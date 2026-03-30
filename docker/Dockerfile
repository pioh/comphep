# CompHEP build environment with LHAPDF 6
#
# Build:  docker build -t comphep-build .
# Usage:  docker run -it comphep-build
#
# Layers:
#   1. base     — OS + build tools + LHAPDF 6
#   2. pdfdata  — download PDF data sets
#   3. build    — compile CompHEP
#   4. test     — working directory for testing
#
# Pin base image by digest for reproducible builds
FROM fedora:43@sha256:781b7642e8bf256e9cf75d2aa58d86f5cc695fd2df113517614e181a5eee9138 AS base

RUN dnf install -y \
    gcc gcc-c++ gcc-gfortran make \
    libX11-devel \
    perl \
    lhapdf lhapdf-devel yaml-cpp-devel \
    python3-lhapdf \
    && dnf clean all

FROM base AS pdfdata

RUN lhapdf install cteq6l1 NNPDF40_nlo_as_01180 2>/dev/null || \
    lhapdf get cteq6l1 NNPDF40_nlo_as_01180 2>/dev/null || true

RUN ls $(lhapdf-config --datadir)/cteq6l1/ && \
    echo "PDF data installed OK" || \
    echo "WARNING: PDF data not available, tests will be limited"

FROM pdfdata AS build

COPY . /comphep
WORKDIR /comphep

RUN make distclean 2>/dev/null || true

RUN ./configure --with-lhapdf --optimise

# -fPIC: make static libs compatible with PIE-default distros (Arch, Ubuntu 17+)
# --allow-multiple-definition: work around duplicate Fortran symbols in external/FeynHiggs
RUN sed -i 's/$/ -fPIC/' CFLAGS && \
    sed -i 's/$/ -fPIC/' F77FLAGS && \
    sed -i 's/$/ -Wl,--allow-multiple-definition/' CLIBS && \
    sed -i 's/$/ -Wl,--allow-multiple-definition/' CLIBS_BASE && \
    make -j$(nproc)

RUN make setup WDIR=/WDIR

FROM build AS test

WORKDIR /WDIR

CMD ["/bin/bash"]
