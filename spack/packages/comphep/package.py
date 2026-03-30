# Copyright CompHEP Collaboration
# SPDX-License-Identifier: LicenseRef-CompHEP

import os

from spack.package import *


class Comphep(MakefilePackage):
    """CompHEP: automatic computation of elementary particle decays
    and collisions at tree level in the Standard Model and beyond.

    CompHEP is a package for symbolic and numerical computations of
    hard scattering processes at colliders. It generates Feynman diagrams,
    calculates squared matrix elements, and performs phase space integration."""

    homepage = "https://gitverse.ru/comphep/comphep"
    git = "https://gitverse.ru/comphep/comphep.git"

    license("LicenseRef-CompHEP")
    maintainers("pioh")

    version("develop", branch="master")
    # Releases — add new versions here (release.sh does this automatically)
    # version("4.7.0", tag="v4.7.0")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")
    depends_on("gmake", type="build")
    depends_on("perl", type="build")
    depends_on("libx11")
    depends_on("libxext")

    variant("lhapdf", default=True, description="Enable LHAPDF 6 PDF sets")
    variant("root", default=False, description="Enable ROOT ntuple output")
    variant("libxml2", default=False, description="Enable HepML via libxml2")
    variant("optimize", default=True, description="Build with -O3 optimization")

    depends_on("lhapdf", when="+lhapdf")
    depends_on("yaml-cpp", when="+lhapdf")
    depends_on("root", when="+root")
    depends_on("libxml2", when="+libxml2")

    def edit(self, spec, prefix):
        """Run CompHEP's custom configure and set up compilers."""

        # Build configure arguments
        args = []
        if "+optimize" in spec:
            args.append("--optimise")
        if "+lhapdf" in spec:
            args.append("--with-lhapdf")
        if "+root" in spec:
            args.append("--with-root")
        if "+libxml2" in spec:
            args.append("--with-libxml")
        if spec.satisfies("target=x86_64:"):
            args.append("--with-m64")

        Executable("./configure")(*args)

        # Override compilers with Spack wrappers (set as env vars by Spack)
        with open("CC", "w") as f:
            f.write(os.environ.get("CC", "gcc"))
        with open("CXX", "w") as f:
            f.write(os.environ.get("CXX", "g++"))
        with open("F77", "w") as f:
            f.write(os.environ.get("FC", os.environ.get("F77", "gfortran")))

        # Work around duplicate Fortran symbols in external/FeynHiggs
        for fname in ["CLIBS", "CLIBS_BASE"]:
            if os.path.isfile(fname):
                with open(fname, "r") as f:
                    content = f.read().rstrip()
                with open(fname, "w") as f:
                    f.write(content + " -Wl,--allow-multiple-definition\n")

    def build(self, spec, prefix):
        make("-j{0}".format(make_jobs))

    def install(self, spec, prefix):
        """Install CompHEP distribution tree.

        CompHEP does not have a standard 'make install'. Instead it uses
        'make setup WDIR=...' to create per-project working directories.

        We install the full distribution tree to prefix so that users can:
          COMPHEP=$(spack location -i comphep)
          make -C $COMPHEP setup WDIR=~/my_project
        """

        # Binaries
        install_tree("bin", prefix.bin)

        # Libraries and headers (needed for user code compilation)
        install_tree("lib", prefix.lib)
        install_tree("include", prefix.include)

        # Data files
        install_tree("models", prefix.join("models"))
        install_tree("help", prefix.join("help"))
        install_tree("strfun", prefix.join("strfun"))
        install_tree("usr", prefix.join("usr"))

        # Source headers needed by user Makefile (usr/Makefile references these)
        src_num_inc = join_path("src", "num", "include")
        if os.path.isdir(src_num_inc):
            mkdirp(prefix.join("src", "num"))
            install_tree(src_num_inc, prefix.join("src", "num", "include"))

        # Documentation
        if os.path.isdir("doc"):
            install_tree("doc", prefix.join("doc"))

        # Build infrastructure (for make setup WDIR=...)
        for f in [
            "Makefile", "configure", "version", "INSTALL", "README",
            "Licence.txt", "CC", "CXX", "CFLAGS", "CLIBS", "CLIBS_BASE",
            "F77", "F77FLAGS", "F77LIBS", "RANLIB",
            "LHAPDF_DATADIR", "ROOTFLAGS", "ROOTLIBS",
        ]:
            if os.path.isfile(f):
                install(f, prefix)
