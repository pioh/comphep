# CompHEP build environment with LHAPDF 6
# Layers designed for cache efficiency:
#   1. OS + system deps (rarely changes)
#   2. LHAPDF 6 install (rarely changes)
#   3. PDF data sets (rarely changes)
#   4. CompHEP source copy + build (changes with code edits)

# ============================================================
# Layer 1: Base system with build dependencies
# ============================================================
FROM fedora:43 AS base

RUN dnf install -y \
    gcc gcc-c++ gcc-gfortran make \
    libX11-devel \
    perl \
    lhapdf lhapdf-devel yaml-cpp-devel \
    python3-lhapdf \
    && dnf clean all

# ============================================================
# Layer 2: Download PDF data sets for testing
# ============================================================
FROM base AS pdfdata

RUN lhapdf install cteq6l1 NNPDF40_nlo_as_01180 2>/dev/null || \
    lhapdf get cteq6l1 NNPDF40_nlo_as_01180 2>/dev/null || true

# Verify data installed (non-fatal if download failed)
RUN ls $(lhapdf-config --datadir)/cteq6l1/ && \
    echo "PDF data installed OK" || \
    echo "WARNING: PDF data not available, tests will be limited"

# ============================================================
# Layer 3: CompHEP build
# ============================================================
FROM pdfdata AS build

# Copy source
COPY . /comphep
WORKDIR /comphep

# Clean any prior artifacts
RUN make distclean 2>/dev/null || true

# Configure with LHAPDF 6
RUN ./configure --with-lhapdf --optimise

# Build (--allow-multiple-definition works around pre-existing duplicate
# Fortran symbols in external/FeynHiggs)
RUN sed -i 's/$/ -Wl,--allow-multiple-definition/' CLIBS && \
    sed -i 's/$/ -Wl,--allow-multiple-definition/' CLIBS_BASE && \
    make -j$(nproc)

# Setup working directory
RUN make setup WDIR=/WDIR

# ============================================================
# Layer 4: Test runner (optional)
# ============================================================
FROM build AS test

WORKDIR /WDIR

# Default: interactive shell
CMD ["/bin/bash"]
