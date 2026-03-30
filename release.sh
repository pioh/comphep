#!/bin/bash
#
# release.sh — Build and publish a trusted CompHEP release
#
# Trust levels:
#   1. SHA256 checksums for all artifacts (SHA256SUMS)
#   2. GPG detached signature of checksums (SHA256SUMS.asc)
#   3. GPG-signed git tag + BUILDINFO for reproducibility
#
set -euo pipefail

# ─── Defaults ───────────────────────────────────────────────────
GITVERSE_API="https://api.gitverse.ru"
ACCEPT_HDR="Accept: application/vnd.gitverse.object+json;version=1"
DOCKERFILE="docker/Dockerfile"

: "${GITVERSE_TOKEN:=}"
: "${GIT_REMOTE:=origin}"
: "${SKIP_DOCKER:=0}"
: "${SKIP_PUSH:=0}"

VERSION=""
BRANCH=""
OWNER=""
REPO=""
GPG_KEY=""
PRERELEASE=0
DO_DOCKER_PUSH=0
DO_DEB=0
DO_RPM=0
DEB_DIST="jammy"
DEB_COMP="main"
RPM_GROUP="fedora/43"

# ─── Help ───────────────────────────────────────────────────────
usage() {
    cat <<'HELP'
Usage: ./release.sh [OPTIONS] VERSION

Build, sign, and publish a CompHEP release.

Arguments:
  VERSION                Version: X.Y.Z or X.Y.Z-suffix (e.g. 4.7.0)

Required:
  -t, --token TOKEN      GitVerse personal access token (or GITVERSE_TOKEN env)
  -o, --owner OWNER      GitVerse repo owner (e.g. pioh, comphep)
  -r, --repo REPO        GitVerse repo name (e.g. comphep)

Signing:
  -g, --gpg-key KEY_ID   GPG key for signing (default: git config user.signingkey)
                          The release will include:
                            - SHA256SUMS       checksums of all artifacts
                            - SHA256SUMS.asc   GPG detached signature
                            - BUILDINFO        build reproducibility metadata
                            - GPG-signed git tag (verifiable via git verify-tag)

Git options:
  -b, --branch BRANCH    Branch to release from (default: current)
      --remote REMOTE     Git remote name (default: origin)
      --prerelease        Mark as pre-release on GitVerse

Build options:
      --skip-docker       Skip Docker build (reuse existing dist/)
      --skip-push         Build locally only, don't push or publish
      --dockerfile PATH   Dockerfile path (default: .gena/docker/Dockerfile)

Package options:
      --docker-push       Push Docker image to gitverse.ru registry
      --deb               Create and upload DEB package
      --rpm               Create and upload RPM package
      --packages          Create DEB + RPM (same as --deb --rpm)
      --deb-dist NAME     DEB distribution (default: jammy)
      --deb-comp NAME     DEB component (default: main)
      --rpm-group NAME    RPM group (default: fedora/43)

  -h, --help              Show this help

Examples:
  ./release.sh -t $TOKEN -o comphep -r comphep --packages 4.7.0
  ./release.sh -t $TOKEN -o pioh -r comphep --prerelease --packages 4.7.0-rc1
  ./release.sh -o pioh -r comphep --skip-push 4.7.0

Prerequisites (checked automatically at start):
  docker      — for building binaries in clean environment
  jq          — for JSON processing
  gpg         — for signing artifacts
  curl        — for GitVerse API
  sha256sum   — for checksums
  GPG key     — set via --gpg-key or git config user.signingkey

Artifacts created in dist/:
  comphep-VERSION-linux-ARCH.tar.gz  Binary tarball
  comphep_VERSION-1_ARCH.deb         DEB package (if --deb)
  comphep-VERSION-1.ARCH.rpm         RPM package (if --rpm)
  SHA256SUMS                         Checksums of all artifacts
  SHA256SUMS.asc                     GPG signature of checksums
  BUILDINFO                          Build metadata for reproducibility
HELP
    exit 0
}

# ─── Parse arguments ────────────────────────────────────────────
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)        usage ;;
        -t|--token)       GITVERSE_TOKEN="$2"; shift 2 ;;
        -o|--owner)       OWNER="$2"; shift 2 ;;
        -r|--repo)        REPO="$2"; shift 2 ;;
        -b|--branch)      BRANCH="$2"; shift 2 ;;
        -g|--gpg-key)     GPG_KEY="$2"; shift 2 ;;
        --remote)         GIT_REMOTE="$2"; shift 2 ;;
        --prerelease)     PRERELEASE=1; shift ;;
        --skip-docker)    SKIP_DOCKER=1; shift ;;
        --skip-push)      SKIP_PUSH=1; shift ;;
        --dockerfile)     DOCKERFILE="$2"; shift 2 ;;
        --docker-push)    DO_DOCKER_PUSH=1; shift ;;
        --deb)            DO_DEB=1; shift ;;
        --rpm)            DO_RPM=1; shift ;;
        --packages)       DO_DEB=1; DO_RPM=1; shift ;;
        --deb-dist)       DEB_DIST="$2"; shift 2 ;;
        --deb-comp)       DEB_COMP="$2"; shift 2 ;;
        --rpm-group)      RPM_GROUP="$2"; shift 2 ;;
        -*)               echo "ERROR: Unknown option: $1"; echo "Try: $0 --help"; exit 1 ;;
        *)
            if [ -z "$VERSION" ]; then VERSION="$1"; shift
            else echo "ERROR: Unexpected argument: $1"; exit 1; fi
            ;;
    esac
done

# ═══════════════════════════════════════════════════════════════
# PREFLIGHT CHECKS — verify everything before doing any work
# ═══════════════════════════════════════════════════════════════
echo ""
echo "=== Preflight checks ==="
err=0

# --- Required arguments ---
[ -z "$VERSION" ] && echo "  ✗ VERSION is required" && err=1
[ -n "$VERSION" ] && ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[a-zA-Z0-9.]+)?$ ]] && \
    echo "  ✗ VERSION must be X.Y.Z or X.Y.Z-suffix (got: $VERSION)" && err=1
[ -z "$OWNER" ] && echo "  ✗ --owner is required" && err=1
[ -z "$REPO" ] && echo "  ✗ --repo is required" && err=1

# --- Required tools ---
for tool in docker jq curl sha256sum gpg git; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool — not found"
        err=1
    fi
done

# --- Docker ---
if [ "$SKIP_DOCKER" != "1" ]; then
    if [ -f "$DOCKERFILE" ]; then
        echo "  ✓ Dockerfile: $DOCKERFILE"
    else
        echo "  ✗ Dockerfile not found: $DOCKERFILE"
        err=1
    fi
    if docker info &>/dev/null; then
        echo "  ✓ Docker daemon running"
    else
        echo "  ✗ Docker daemon not running"
        err=1
    fi
fi

# --- GPG key ---
if [ -z "$GPG_KEY" ]; then
    GPG_KEY=$(git config user.signingkey 2>/dev/null || true)
fi
if [ -n "$GPG_KEY" ]; then
    if gpg --list-secret-keys "$GPG_KEY" &>/dev/null; then
        GPG_UID=$(gpg --list-secret-keys --with-colons "$GPG_KEY" 2>/dev/null | grep '^uid' | head -1 | cut -d: -f10)
        echo "  ✓ GPG key: $GPG_KEY ($GPG_UID)"
    else
        echo "  ✗ GPG key $GPG_KEY not found in keyring"
        err=1
    fi
else
    echo "  ✗ No GPG key — set --gpg-key or git config user.signingkey"
    echo "    Generate one: gpg --full-generate-key"
    err=1
fi

# --- GitVerse token ---
if [ "$SKIP_PUSH" != "1" ]; then
    if [ -z "$GITVERSE_TOKEN" ]; then
        echo "  ✗ --token required for publishing (or --skip-push for local)"
        err=1
    else
        code=$(curl -s -w "%{http_code}" -o /dev/null \
            -H "Authorization: Bearer ${GITVERSE_TOKEN}" -H "$ACCEPT_HDR" \
            "${GITVERSE_API}/user" 2>/dev/null || echo "000")
        if [ "$code" = "200" ]; then
            echo "  ✓ GitVerse token valid"
        else
            echo "  ✗ GitVerse token invalid (HTTP $code)"
            echo "    Create at: https://gitverse.ru/settings/applications"
            err=1
        fi
    fi
fi

# --- Git state ---
if [ -n "$VERSION" ]; then
    TAG="v${VERSION}"
    if git tag -l | grep -qx "$TAG"; then
        echo "  ✗ Tag $TAG already exists"
        err=1
    else
        echo "  ✓ Tag $TAG available"
    fi
fi

if [ -n "$(git status --porcelain --ignore-submodules)" ]; then
    echo "  ✗ Working directory not clean"
    git status --short | head -5 | sed 's/^/    /'
    err=1
else
    echo "  ✓ Working directory clean"
fi

# --- Bail if any check failed ---
if [ $err -gt 0 ]; then
    echo ""
    echo "Preflight FAILED — fix the issues above before releasing."
    echo "Try: $0 --help"
    exit 1
fi

echo ""
echo "  All checks passed!"
echo ""

# ─── Computed variables ─────────────────────────────────────────
TAG="v${VERSION}"

HOST_ARCH=$(uname -m)
case "$HOST_ARCH" in
    x86_64)  ARCH="x86_64"; DEB_ARCH="amd64" ;;
    aarch64) ARCH="aarch64"; DEB_ARCH="arm64" ;;
    armv7l)  ARCH="armv7l";  DEB_ARCH="armhf" ;;
    *)       ARCH="$HOST_ARCH"; DEB_ARCH="$HOST_ARCH" ;;
esac

# RPM version: split hyphen suffix (RPM doesn't allow hyphens in Version)
RPM_VERSION="${VERSION%%-*}"           # 4.7.0-rc1 → 4.7.0
RPM_RELEASE="${VERSION#*-}"           # 4.7.0-rc1 → rc1
[ "$RPM_RELEASE" = "$VERSION" ] && RPM_RELEASE="1"  # 4.7.0 → 1
RPM_RELEASE="${RPM_RELEASE//[-.]/_}"  # sanitize: rc.1 → rc_1

BIN_TARBALL="comphep-${VERSION}-linux-${ARCH}.tar.gz"
DEB_FILE="comphep_${VERSION}-1_${DEB_ARCH}.deb"
RPM_FILE="comphep-${VERSION}-1.${ARCH}.rpm"
DIST_DIR="$(pwd)/dist"
BIN_DIST_NAME="comphep-${VERSION}-linux-${ARCH}"
BIN_DIST_DIR="/tmp/${BIN_DIST_NAME}"
DOCKER_IMAGE="comphep-build:${TAG}"
REGISTRY_IMAGE="gitverse.ru/${OWNER}/comphep"

# ─── Switch branch ──────────────────────────────────────────────
if [ -n "$BRANCH" ]; then
    CUR=$(git branch --show-current)
    [ "$CUR" != "$BRANCH" ] && echo "Switching to branch: $BRANCH" && git checkout "$BRANCH"
fi

# ─── Plan ───────────────────────────────────────────────────────
echo "╔═══════════════════════════════════════════════════╗"
echo "║  CompHEP Release: $TAG"
echo "╠═══════════════════════════════════════════════════╣"
echo "║  Owner/Repo  : ${OWNER}/${REPO}"
echo "║  Branch      : $(git branch --show-current)"
echo "║  Architecture: ${ARCH} (deb: ${DEB_ARCH})"
echo "║  GPG key     : ${GPG_KEY}"
printf "║  Packages    : "
pkgs=""
[ "$DO_DOCKER_PUSH" = "1" ] && pkgs="${pkgs}docker "
[ "$DO_DEB" = "1" ] && pkgs="${pkgs}deb "
[ "$DO_RPM" = "1" ] && pkgs="${pkgs}rpm "
echo "${pkgs:-tarball only}"
echo "║  Pre-release : $([ "$PRERELEASE" = "1" ] && echo "yes" || echo "no")"
echo "╚═══════════════════════════════════════════════════╝"
echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 1: Clean + update version
# ═══════════════════════════════════════════════════════════════
echo "=== [1] Preparing ==="
rm -rf "$DIST_DIR"
mkdir -p "$DIST_DIR"

OLD_VERSION=$(cat version)
echo "$VERSION" > version
echo "  Version: $OLD_VERSION → $VERSION"
echo ""

rollback() { echo ""; echo "!!! Failed. Restoring version."; echo "$OLD_VERSION" > version; }
trap rollback ERR

# ═══════════════════════════════════════════════════════════════
# STEP 2: Docker build
# ═══════════════════════════════════════════════════════════════
if [ "$SKIP_DOCKER" != "1" ]; then
    echo "=== [2] Building in Docker ==="
    docker build -f "$DOCKERFILE" -t "$DOCKER_IMAGE" --target build .
    echo ""
    echo "  ✓ Build successful"

    # Record Docker image digest for reproducibility
    DOCKER_DIGEST=$(docker inspect --format='{{.Id}}' "$DOCKER_IMAGE")
    BASE_IMAGE=$(head -20 "$DOCKERFILE" | grep '^FROM' | tail -1 | awk '{print $2}')
    echo ""
else
    echo "=== [2] Docker build — SKIPPED ==="
    DOCKER_DIGEST="skipped"
    BASE_IMAGE="skipped"
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 3: Extract binary distribution + create tarball
# ═══════════════════════════════════════════════════════════════
if [ "$SKIP_DOCKER" != "1" ]; then
    echo "=== [3] Creating binary tarball ==="
    rm -rf "$BIN_DIST_DIR"
    mkdir -p "$BIN_DIST_DIR"

    CID=$(docker create "$DOCKER_IMAGE")

    for dir in bin lib include models help strfun usr doc; do
        docker cp "$CID:/comphep/$dir" "$BIN_DIST_DIR/$dir" 2>/dev/null || true
    done
    for f in Makefile configure INSTALL version README README.md Licence.txt \
             CC CXX CFLAGS CLIBS CLIBS_BASE F77 F77FLAGS F77LIBS RANLIB \
             LHAPDF_DATADIR ROOTFLAGS ROOTLIBS; do
        docker cp "$CID:/comphep/$f" "$BIN_DIST_DIR/$f" 2>/dev/null || true
    done
    # Copy src/*/include/ directories (needed by usr/Makefile for userFun.c)
    for subdir in $(docker run --rm "$DOCKER_IMAGE" find /comphep/src -maxdepth 2 -type d -name include); do
        REL=${subdir#/comphep/}
        mkdir -p "$BIN_DIST_DIR/$(dirname "$REL")"
        docker cp "$CID:/comphep/$REL" "$BIN_DIST_DIR/$REL" 2>/dev/null || true
    done
    # Copy .o files from src/fann/ (linked by usr/Makefile)
    mkdir -p "$BIN_DIST_DIR/src/fann"
    docker cp "$CID:/comphep/src/fann/fannvar.o" "$BIN_DIST_DIR/src/fann/" 2>/dev/null || true
    docker cp "$CID:/comphep/src/fann/unweighter.o" "$BIN_DIST_DIR/src/fann/" 2>/dev/null || true
    docker cp "$CID:/comphep/spack" "$BIN_DIST_DIR/spack" 2>/dev/null || true

    docker rm "$CID" >/dev/null

    (cd /tmp && tar czf "$DIST_DIR/$BIN_TARBALL" "$BIN_DIST_NAME")
    echo "  ✓ $BIN_TARBALL ($(du -h "$DIST_DIR/$BIN_TARBALL" | cut -f1))"
    echo ""
else
    echo "=== [3] Binary tarball — SKIPPED ==="
    if [ -f "$DIST_DIR/$BIN_TARBALL" ] && [ ! -d "$BIN_DIST_DIR" ]; then
        (cd /tmp && tar xzf "$DIST_DIR/$BIN_TARBALL")
    fi
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 4: Create DEB package
# ═══════════════════════════════════════════════════════════════
if [ "$DO_DEB" = "1" ]; then
    echo "=== [4] Creating DEB package ==="
    [ ! -d "$BIN_DIST_DIR" ] && echo "  ERROR: Binary dist not found" && exit 1

    DEB_ROOT="/tmp/comphep-deb-$$"
    rm -rf "$DEB_ROOT"
    mkdir -p "$DEB_ROOT/DEBIAN" "$DEB_ROOT/opt/comphep" "$DEB_ROOT/usr/bin"

    cp -a "$BIN_DIST_DIR"/* "$DEB_ROOT/opt/comphep/"
    ln -sf /opt/comphep/usr/comphep "$DEB_ROOT/usr/bin/comphep"

    cat > "$DEB_ROOT/DEBIAN/control" <<EOF
Package: comphep
Version: ${VERSION}-1
Architecture: ${DEB_ARCH}
Maintainer: CompHEP Collaboration <sherstnv@theory.sinp.msu.ru>
Depends: libx11-6, libxext6, gcc, g++, gfortran, make
Recommends: lhapdf
Section: science
Priority: optional
Homepage: https://gitverse.ru/comphep/comphep
Description: CompHEP - automatic computations in high energy physics
 CompHEP is a package for symbolic and numerical calculations
 of hard scattering processes at colliders.
 .
 After install: make -C /opt/comphep setup WDIR=~/comphep_work
EOF

    docker run --rm -v "/tmp:/tmp" -v "$DIST_DIR:/dist" \
        -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
        ubuntu:22.04 \
        sh -c "dpkg-deb --root-owner-group --build $DEB_ROOT /dist/$DEB_FILE && \
               chown \$HOST_UID:\$HOST_GID /dist/$DEB_FILE"

    rm -rf "$DEB_ROOT" 2>/dev/null || docker run --rm -v "/tmp:/tmp" ubuntu:22.04 rm -rf "$DEB_ROOT"
    echo "  ✓ $DEB_FILE ($(du -h "$DIST_DIR/$DEB_FILE" | cut -f1))"
    echo ""
else
    echo "=== [4] DEB — SKIPPED ==="
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 5: Create RPM package
# ═══════════════════════════════════════════════════════════════
if [ "$DO_RPM" = "1" ]; then
    echo "=== [5] Creating RPM package ==="
    [ ! -d "$BIN_DIST_DIR" ] && echo "  ERROR: Binary dist not found" && exit 1

    RPM_ROOT="/tmp/comphep-rpm-$$"
    rm -rf "$RPM_ROOT"
    mkdir -p "$RPM_ROOT"/{SPECS,SOURCES,BUILD,RPMS,SRPMS,BUILDROOT}
    cp "$DIST_DIR/$BIN_TARBALL" "$RPM_ROOT/SOURCES/"

    cat > "$RPM_ROOT/SPECS/comphep.spec" <<EOF
Name:           comphep
Version:        ${RPM_VERSION}
Release:        ${RPM_RELEASE}
Summary:        CompHEP - automatic computations in high energy physics
License:        CompHEP
URL:            https://gitverse.ru/comphep/comphep
Source0:        ${BIN_TARBALL}
AutoReqProv:    no
%global debug_package %{nil}
Requires:       libX11, libXext, gcc, gcc-c++, gcc-gfortran, make

%description
CompHEP - symbolic and numerical calculations of scattering processes.
After install: make -C /opt/comphep setup WDIR=~/comphep_work

%prep
%setup -q -n ${BIN_DIST_NAME}

%build

%install
mkdir -p %{buildroot}/opt/comphep
cp -a . %{buildroot}/opt/comphep/
mkdir -p %{buildroot}/usr/bin
ln -sf /opt/comphep/usr/comphep %{buildroot}/usr/bin/comphep

%files
/opt/comphep/
/usr/bin/comphep
EOF

    docker run --rm \
        -v "$RPM_ROOT:/rpmbuild" -v "$DIST_DIR:/dist" \
        -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
        fedora:43 \
        sh -c 'dnf install -y rpm-build >/dev/null 2>&1 && \
               rpmbuild -bb --define "_topdir /rpmbuild" /rpmbuild/SPECS/comphep.spec && \
               cp /rpmbuild/RPMS/*/*.rpm /dist/ && \
               chown $HOST_UID:$HOST_GID /dist/*.rpm'

    rm -rf "$RPM_ROOT" 2>/dev/null || docker run --rm -v "/tmp:/tmp" fedora:43 rm -rf "$RPM_ROOT"
    RPM_BUILT=$(find "$DIST_DIR" -name "comphep-*.rpm" -newer "$DIST_DIR/$BIN_TARBALL" 2>/dev/null | head -1)
    if [ -n "$RPM_BUILT" ]; then
        [ "$(basename "$RPM_BUILT")" != "$RPM_FILE" ] && mv "$RPM_BUILT" "$DIST_DIR/$RPM_FILE"
        echo "  ✓ $RPM_FILE ($(du -h "$DIST_DIR/$RPM_FILE" | cut -f1))"
    else
        echo "  ERROR: RPM build produced no output"
    fi
    echo ""
else
    echo "=== [5] RPM — SKIPPED ==="
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 6: SHA256 checksums + GPG signature + BUILDINFO
# ═══════════════════════════════════════════════════════════════
echo "=== [6] Signing artifacts ==="

# --- BUILDINFO: reproducibility metadata ---
BUILDINFO="$DIST_DIR/BUILDINFO.txt"
cat > "$BUILDINFO" <<EOF
CompHEP Build Information
=========================
Version:        ${VERSION}
Tag:            ${TAG}
Commit:         $(git rev-parse HEAD)
Branch:         $(git branch --show-current)
Build date:     $(date -u +%Y-%m-%dT%H:%M:%SZ)
Build host:     $(uname -n)
Build arch:     ${ARCH}
OS:             $(uname -s) $(uname -r)
Docker image:   ${DOCKER_DIGEST}
Base image:     ${BASE_IMAGE}
Dockerfile:     $(sha256sum "$DOCKERFILE" 2>/dev/null | cut -d' ' -f1 || echo "N/A")
GPG key:        ${GPG_KEY}

To reproduce this build:
  git checkout ${TAG}
  docker build -f ${DOCKERFILE} -t comphep-build:${TAG} --target build .
EOF
echo "  ✓ BUILDINFO.txt"

# --- SHA256 checksums ---
CHECKSUMS="$DIST_DIR/SHA256SUMS.txt"
(cd "$DIST_DIR" && sha256sum *.tar.gz *.deb *.rpm BUILDINFO.txt 2>/dev/null > SHA256SUMS.txt)
echo "  ✓ SHA256SUMS.txt"
cat "$CHECKSUMS" | sed 's/^/    /'

# --- GPG detached signature ---
gpg --default-key "$GPG_KEY" --detach-sign --armor --output "$DIST_DIR/SHA256SUMS.txt.sig" "$CHECKSUMS"
echo "  ✓ SHA256SUMS.txt.sig (signed by $GPG_KEY)"

echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 7: GPG-signed git tag + commit
# ═══════════════════════════════════════════════════════════════
echo "=== [7] Commit + signed tag ==="

git add version
git commit -m "Release $TAG"
git tag -s -u "$GPG_KEY" "$TAG" -m "CompHEP $VERSION

Signed-off-by: $(gpg --list-keys --with-colons "$GPG_KEY" 2>/dev/null | grep '^uid' | head -1 | cut -d: -f10)
SHA256SUMS included in release assets."

echo "  ✓ Commit + GPG-signed tag $TAG"
echo "  Verify: git verify-tag $TAG"
echo ""

trap - ERR

# ═══════════════════════════════════════════════════════════════
# STEP 8: Git push
# ═══════════════════════════════════════════════════════════════
if [ "$SKIP_PUSH" != "1" ]; then
    echo "=== [8] Pushing to $GIT_REMOTE ==="
    git push "$GIT_REMOTE" HEAD
    git push "$GIT_REMOTE" "$TAG"
    echo "  ✓ Pushed"
    echo ""
else
    echo "=== [8] Git push — SKIPPED ==="
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 9: GitVerse release + upload assets
# ═══════════════════════════════════════════════════════════════
if [ "$SKIP_PUSH" != "1" ] && [ -n "$GITVERSE_TOKEN" ]; then
    echo "=== [9] Creating GitVerse release ==="

    PRE_JSON="false"; [ "$PRERELEASE" = "1" ] && PRE_JSON="true"

    # Package registry URLs
    PKG_BASE="https://gitverse.ru/${OWNER}/-/packages"
    DEB_PKG_URL="${PKG_BASE}/debian/comphep/${VERSION}-1"
    RPM_PKG_URL="${PKG_BASE}/rpm/comphep/${RPM_VERSION}"

    BODY="## CompHEP ${VERSION}

### Install via package manager

**Debian/Ubuntu (apt):**
\`\`\`bash
sudo curl https://gitverse.ru/api/packages/${OWNER}/debian/repository.key -o /etc/apt/keyrings/gitverse-${OWNER}.asc
echo 'deb [signed-by=/etc/apt/keyrings/gitverse-${OWNER}.asc] https://gitverse.ru/api/packages/${OWNER}/debian ${DEB_DIST} ${DEB_COMP}' | sudo tee /etc/apt/sources.list.d/gitverse.list
sudo apt update && sudo apt install comphep
\`\`\`"

    [ "$DO_DEB" = "1" ] && [ -f "$DIST_DIR/$DEB_FILE" ] && \
        BODY="${BODY}
> [DEB package details](${DEB_PKG_URL})"

    BODY="${BODY}

**Fedora/RHEL (dnf):**
\`\`\`bash
sudo dnf config-manager --add-repo https://gitverse.ru/api/packages/${OWNER}/rpm/${RPM_GROUP}.repo
sudo dnf makecache && sudo dnf install comphep
\`\`\`"

    [ "$DO_RPM" = "1" ] && [ -f "$DIST_DIR/$RPM_FILE" ] && \
        BODY="${BODY}
> [RPM package details](${RPM_PKG_URL})"

    BODY="${BODY}

### Install from binary tarball
\`\`\`bash
tar xzf ${BIN_TARBALL}
cd ${BIN_DIST_NAME}
make setup WDIR=~/comphep_work
cd ~/comphep_work && ./comphep
\`\`\`

### Install from source
\`\`\`bash
tar xzf comphep-*.tar.gz && cd comphep-*/
./configure --with-lhapdf      # optional: requires lhapdf-config in PATH
make && make setup WDIR=~/comphep_work
\`\`\`

### Install via Spack
\`\`\`bash
git clone https://gitverse.ru/${OWNER}/${REPO}.git
spack repo add ${REPO}/spack && spack install comphep@${VERSION} +lhapdf
\`\`\`"

    [ "$DO_DOCKER_PUSH" = "1" ] && \
        BODY="${BODY}

### Docker
\`\`\`bash
docker pull ${REGISTRY_IMAGE}:${VERSION}
\`\`\`"

    BODY="${BODY}

### Verify integrity
\`\`\`bash
sha256sum -c SHA256SUMS.txt
gpg --verify SHA256SUMS.txt.sig SHA256SUMS.txt
git verify-tag ${TAG}
\`\`\`

### Assets
| File | Description |
|------|-------------|
| \`${BIN_TARBALL}\` | Binary (Linux ${ARCH}) |
| \`SHA256SUMS.txt\` | Checksums |
| \`SHA256SUMS.txt.sig\` | GPG signature |
| \`BUILDINFO.txt\` | Build metadata |"

    RESP=$(curl -s -w "\n%{http_code}" -X POST \
        -H "Authorization: Bearer ${GITVERSE_TOKEN}" -H "$ACCEPT_HDR" \
        -H "Content-Type: application/json" \
        "${GITVERSE_API}/repos/${OWNER}/${REPO}/releases" \
        -d "$(printf '%s' "$BODY" | jq -Rs --arg t "$TAG" --arg n "CompHEP $VERSION" --argjson p "$PRE_JSON" \
            '{tag_name:$t, name:$n, body:., draft:false, prerelease:$p}')")

    HTTP=$(echo "$RESP" | tail -1)
    JSON=$(echo "$RESP" | sed '$d')
    RELEASE_ID=$(echo "$JSON" | jq -r '.id // empty')

    if [ -z "$RELEASE_ID" ] || [ "$HTTP" -ge 400 ]; then
        echo "  ERROR: Failed to create release (HTTP $HTTP)"
        echo "  $JSON"
    else
        echo "  ✓ Release created (ID: $RELEASE_ID)"

        # Upload: binary tarball + signing artifacts
        for f in "$DIST_DIR/$BIN_TARBALL" "$DIST_DIR/SHA256SUMS.txt" "$DIST_DIR/SHA256SUMS.txt.sig" "$DIST_DIR/BUILDINFO.txt"; do
            [ ! -f "$f" ] && continue
            NAME=$(basename "$f")
            ENCODED_NAME=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$NAME'))")
            printf "  Uploading %s... " "$NAME"
            URESP=$(curl -s -w "\n%{http_code}" -X POST \
                -H "Authorization: Bearer ${GITVERSE_TOKEN}" -H "$ACCEPT_HDR" \
                -F "attachment=@${f}" \
                "${GITVERSE_API}/repos/${OWNER}/${REPO}/releases/${RELEASE_ID}/assets?name=${ENCODED_NAME}")
            UC=$(echo "$URESP" | tail -1)
            [ "$UC" = "201" ] && echo "OK" || { echo "FAILED (HTTP $UC)"; echo "$URESP" | sed '$d' | head -1; }
        done

        echo ""
        echo "  → https://gitverse.ru/${OWNER}/${REPO}/releases/tag/${TAG}"
    fi
    echo ""
else
    echo "=== [9] GitVerse release — SKIPPED ==="
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 10: Docker push
# ═══════════════════════════════════════════════════════════════
if [ "$DO_DOCKER_PUSH" = "1" ] && [ "$SKIP_PUSH" != "1" ]; then
    echo "=== [10] Pushing Docker image ==="
    echo "${GITVERSE_TOKEN}" | docker login -u "$OWNER" --password-stdin gitverse.ru 2>/dev/null \
        && echo "  ✓ Logged in" \
        || echo "  ERROR: docker login failed"

    DOCKER_TAG="${VERSION}"
    docker tag "$DOCKER_IMAGE" "${REGISTRY_IMAGE}:${DOCKER_TAG}"
    echo "  Pushing ${REGISTRY_IMAGE}:${DOCKER_TAG}..."
    docker push "${REGISTRY_IMAGE}:${DOCKER_TAG}" 2>&1 | tail -3 \
        && echo "  ✓ Pushed :${DOCKER_TAG}" \
        || echo "  ERROR: push failed"

    if [ "$PRERELEASE" != "1" ]; then
        docker tag "$DOCKER_IMAGE" "${REGISTRY_IMAGE}:latest"
        docker push "${REGISTRY_IMAGE}:latest" 2>&1 | tail -3 \
            && echo "  ✓ Pushed :latest" || echo "  ERROR: push failed"
    fi
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 11: Upload DEB to GitVerse packages
# ═══════════════════════════════════════════════════════════════
if [ "$DO_DEB" = "1" ] && [ "$SKIP_PUSH" != "1" ] && [ -f "$DIST_DIR/$DEB_FILE" ]; then
    echo "=== [11] Uploading DEB → packages/${OWNER}/debian ==="
    DRESP=$(curl -s -w "\n%{http_code}" --user "${OWNER}:${GITVERSE_TOKEN}" \
        --upload-file "$DIST_DIR/$DEB_FILE" \
        "https://gitverse.ru/api/packages/${OWNER}/debian/pool/${DEB_DIST}/${DEB_COMP}/upload")
    UC=$(echo "$DRESP" | tail -1)
    [ "$UC" = "201" ] && echo "  ✓ OK" || { echo "  FAILED (HTTP $UC)"; echo "$DRESP" | sed '$d' | head -2; }
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# STEP 12: Upload RPM to GitVerse packages
# ═══════════════════════════════════════════════════════════════
if [ "$DO_RPM" = "1" ] && [ "$SKIP_PUSH" != "1" ] && [ -f "$DIST_DIR/$RPM_FILE" ]; then
    echo "=== [12] Uploading RPM → packages/${OWNER}/rpm ==="
    RRESP=$(curl -s -w "\n%{http_code}" --user "${OWNER}:${GITVERSE_TOKEN}" \
        --upload-file "$DIST_DIR/$RPM_FILE" \
        "https://gitverse.ru/api/packages/${OWNER}/rpm/${RPM_GROUP}/upload")
    UC=$(echo "$RRESP" | tail -1)
    [ "$UC" = "201" ] && echo "  ✓ OK" || { echo "  FAILED (HTTP $UC)"; echo "$RRESP" | sed '$d' | head -2; }
    echo ""
fi

# ═══════════════════════════════════════════════════════════════
# Cleanup
# ═══════════════════════════════════════════════════════════════
rm -rf "$BIN_DIST_DIR"

# ═══════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════
echo "╔═══════════════════════════════════════════════════╗"
echo "║  Release $TAG complete!"
echo "╠═══════════════════════════════════════════════════╣"
for f in "$BIN_TARBALL" "$DEB_FILE" "$RPM_FILE" "SHA256SUMS.txt" "SHA256SUMS.txt.sig" "BUILDINFO.txt"; do
    [ -f "$DIST_DIR/$f" ] && printf "║  %-45s %s\n" "$f" "$(du -h "$DIST_DIR/$f" | cut -f1)"
done
[ "$DO_DOCKER_PUSH" = "1" ] && echo "║  ${REGISTRY_IMAGE}:${VERSION}"
echo "╠═══════════════════════════════════════════════════╣"
echo "║  Verify: sha256sum -c dist/SHA256SUMS.txt"
echo "║          gpg --verify dist/SHA256SUMS.txt.sig"
echo "║          git verify-tag $TAG"
echo "╚═══════════════════════════════════════════════════╝"
