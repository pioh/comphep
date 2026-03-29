####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Authors: A. Sherstnev, A.Kryukov                  #
####################################################
VERSION=$(shell cat version)
####################################################
# Correct the value of this variable:
#
# The WDIR variable have to point to the 
# user directory where all necessary files 
# for user install will live.
####################################################
WDIR=../comphep_${VERSION}_test

######################################################################################
SHELL=/bin/sh
COMPHEP_DISTR=comphep-${VERSION}
PWD=$(shell pwd)
OSTYPE=$(shell uname)
CHEPTMPDIR="/tmp/comphep$(shell id -u)"
CHEPDIR=.

CC=$(shell cat CC)
CFLAGS=$(shell cat CFLAGS)
CLIBS=$(shell cat CLIBS)
CLIBS_BASE=$(shell cat CLIBS_BASE 2>/dev/null || cat CLIBS)
F77=$(shell cat F77)
F77FLAGS=$(shell cat F77FLAGS)
F77LIBS=$(shell cat F77LIBS)
RANLIB=$(shell cat RANLIB)

LHAPDF_DATADIR=$(shell cat LHAPDF_DATADIR 2>/dev/null)
ROOTFLAGS=
ROOTLIBS=
ifneq (${ROOTSYS},)
ROOTFLAGS=$(shell cat ROOTFLAGS)
ROOTLIBS=$(shell cat ROOTLIBS)
endif

SUBDIR:=service2 chep_crt polynom num plot symb tab mssmlib fann
CS=src
SUB_SUBDIR:=$(addprefix $(CS)/, $(SUBDIR))
COMPHEP_BIN_FILES=s_comphep.exe tab_view.exe archiv.exe slhasuspect.exe \
                  mix.exe diag_viewer.exe diag_viewer2.exe pdf_convertor.exe \
                  translator.exe cascade.exe rtupler.exe lanhep.exe \
                  oneclick.exe
BINDIST_FILES=${COMPHEP_BIN_FILES} models help usr include lib Makefile \
              CC CFLAGS CLIBS CLIBS_BASE F77 F77FLAGS F77LIBS RANLIB ROOTFLAGS ROOTLIBS LHAPDF_DATADIR \
              launch_n_comphep configure \
              INSTALL strfun comphep32.png comphep64.png comphep128.png \
              RELEASE-NOTES README
DOC_FILES=${PWD}/doc/Makefile ${PWD}/doc/*HOWTO ${PWD}/doc/RELEASE-NOTES

######################################################################################
.PHONY: all precompile compile link doc pdf

all: link
	@echo " "
	@echo "******************************************************************"
	@echo "* Binaries for CompHEP-${VERSION} has been successfully prepared *"
	@echo "*                                                                *"
	@echo "*      Create a user working directory using the command         *"
	@echo "*         make setup WDIR=path_to_your_user_work_dir             *"
	@echo "*  Note 1: Do not use '~' to refer to you home directory         *"
	@echo "*          Use the environment variable HOME                     *"
	@echo "******************************************************************"

pdf:
	$(CC)  $(CFLAGS) -o bin/pdf_convertor.exe src/service2/cteq_internal.c $(CLIBS)

link: compile
	if (test ! -f lib/libevents.a  ); then $(RANLIB) lib/libevents.a; fi
	if (test ! -f lib/libexternal.a); then $(RANLIB) lib/libexternal.a; fi
	if (test ! -f lib/libmix.a     ); then $(RANLIB) lib/libmix.a; fi
	if (test ! -f lib/libmssm.a    ); then $(RANLIB) lib/libmssm.a; fi
	if (test ! -f lib/libnum.a     ); then $(RANLIB) lib/libnum.a; fi
	if (test ! -f lib/libpdf.a     ); then $(RANLIB) lib/libpdf.a; fi
	if (test ! -f lib/libserv.a    ); then $(RANLIB) lib/libserv.a; fi
	if (test ! -f lib/libsymb.a    ); then $(RANLIB) lib/libsymb.a; fi
	if (test ! -f lib/libtab.a     ); then $(RANLIB) lib/libtab.a; fi
	if (test ! -f lib/libtranls.a  ); then $(RANLIB) lib/libtranls.a; fi
#	if (test ! -f lib/liblanhep.a  ); then $(RANLIB) lib/liblanhep.a; fi

	$(CC) $(CFLAGS) -o bin/tab_view.exe     src/plot/view_tab.o                    -lnum -lserv                    $(CLIBS)
	$(CC) $(CFLAGS) -o bin/mix.exe          src/num/mix.o                          -lmix -levents -lnum -lserv     $(CLIBS)
	$(CC) $(CFLAGS) -o bin/translator.exe   src/num/translator.o src/num/alphas2.o -ltranls -levents -lpdf -lserv  $(CLIBS)
	$(CC) $(CFLAGS) -o bin/cascade.exe      src/num/cascade.o                      -ltranls -levents -lpdf -lserv  $(CLIBS)
#	$(CC) $(CFLAGS) -o bin/lanhep.exe       src/lanhep/main.o                      -llanhep                        $(CLIBS) $(F77LIBS)
	$(CC) $(CFLAGS) -o bin/s_comphep.exe    src/symb/s_comphep.o                   -lsymb -lserv -lmssm -lexternal $(CLIBS_BASE) $(F77LIBS)
	$(CC) $(CFLAGS) -o bin/archiv.exe       src/symb/archiv.o                      -lsymb -lserv -lmssm -lexternal $(CLIBS_BASE) $(F77LIBS)
	$(CC) $(CFLAGS) -o bin/diag_viewer.exe  src/num/diag_viewer.o                  -lsymb -lserv -lmssm -lexternal $(CLIBS_BASE) $(F77LIBS)
	$(CC) $(CFLAGS) -o bin/diag_viewer2.exe src/num/diag_viewer2.o                 -lsymb -lserv -lmssm -lexternal $(CLIBS_BASE) $(F77LIBS)
	$(CC) $(CFLAGS) -o bin/pdf_reweighter.exe src/num/pdf_reweighter.o             -levents -lnum -lpdf -lserv     $(CLIBS) $(F77LIBS)
ifneq (${ROOTSYS},)
	$(CXX) $(CFLAGS) -o bin/rtupler.exe src/num/rtuple.o -ltranls -lnum -levents -lpdf -lserv $(CLIBS) $(F77LIBS) $(ROOTLIBS)
endif
	$(F77) $(F77FLAGS) -o bin/slhasuspect.exe external/suspect2_call.o -lexternal $(CLIBS) $(F77LIBS)

compile: precompile
	$(MAKE) -C src/service2
	$(MAKE) -C src/chep_crt
	$(MAKE) -C src/num
	$(MAKE) -C src/polynom
	$(MAKE) -C src/symb
	$(MAKE) -C src/mssmlib
	$(MAKE) -C src/plot
	$(MAKE) -C src/tab
#	$(MAKE) -C src/lanhep
	$(MAKE) -C external
	$(MAKE) -C src/fann

precompile:
	@if (test ! -f $(CHEPDIR)/CC); then \
	  echo " ";\
	  echo "*******************************************";\
	  echo " 'make' can not find $(CHEPDIR)/CC file.";\
	  echo " Launch the root 'configure' script first.";\
	  echo "*******************************************";\
	  echo " ";\
	  exit 99; \
	fi
	@if (test ! -f $(CHEPDIR)/CFLAGS); then \
	  echo " ";\
	  echo "*******************************************";\
	  echo " 'make' can not find $(CHEPDIR)/CFLAGS file.";\
	  echo " Launch the root 'configure' script first.";\
	  echo "*******************************************";\
	  echo " ";\
	  exit 99; \
	fi
	@if (test ! -f $(CHEPDIR)/CLIBS); then \
	  echo " ";\
	  echo "*******************************************";\
	  echo " 'make' can not find $(CHEPDIR)/CLIBS file.";\
	  echo " Launch the root 'configure' script first.";\
	  echo "*******************************************";\
	  echo " ";\
	  exit 99; \
	fi
	@if (test ! -f $(CHEPDIR)/F77LIBS); then \
	  echo " ";\
	  echo "*******************************************";\
	  echo " 'make' can not find $(CHEPDIR)/F77LIBS file.";\
	  echo " Launch the root 'configure' script first.";\
	  echo "*******************************************";\
	  echo " ";\
	  exit 99; \
	fi
	@if (test ! -f $(CHEPDIR)/RANLIB); then \
	  echo " ";\
	  echo "*******************************************";\
	  echo " 'make' can not find $(CHEPDIR)/RANLIB file.";\
	  echo " Launch the root 'configure' script first.";\
	  echo "*******************************************";\
	  echo " ";\
	  exit 99; \
	fi

doc:
	$(MAKE) -C doc

######################################################################################
.PHONY: clean myclean execlean inclean distclean totclean

totclean: clean execlean inclean
	@echo "Total cleaning has been done..."

distclean: clean execlean inclean
	${MAKE} -C doc clean 
	@rm -f test/results/* test/tmp/* test/core*
	@echo "Total cleaning has been done..."

inclean:
	@rm -f CC F77 F77FLAGS CFLAGS CLIBS CLIBS_BASE F77LIBS RANLIB ROOTFLAGS ROOTLIBS LHAPDF_DATADIR

execlean:
	@rm -f bin/*.exe lib/*.a
	@rm -f strfun/LHAIndex-comphep.txt

myclean:
	@for dir in $(SUB_SUBDIR); do $(MAKE) -C $$dir clean; done
	@echo "Partial cleaning has been done..."

clean:
	@for dir in $(SUB_SUBDIR); do $(MAKE) -C $$dir clean; done
	${MAKE} -C external clean

######################################################################################
.PHONY: fullsetup setup lhapdf oneclicksetup addsetup

addsetup:
	@if (test ! -d $(WDIR)); then \
	   echo "***** Install WDIR with make setup!"; \
	   exit 97; \
	fi
	@cp usr/nn.c $(WDIR)/usr/nn.c
	@cp usr/unweighter $(WDIR)/unweighter.sh
	@chmod 744 $(WDIR)/unweighter.sh
	$(MAKE) -C $(WDIR)/usr fann CC='$(CC)' CFLAGS='$(CFLAGS)'

fullsetup: setup
	@cp usr/mix            $(WDIR)/mix.sh
	@cp usr/cascade        $(WDIR)/cascade.sh
	@cp usr/translator     $(WDIR)/translator.sh
	@cp usr/addcut         $(WDIR)/addcut.sh
	@cp usr/pdf_reweighter $(WDIR)/pdf_reweighter.sh
	@cp usr/config.reweight $(WDIR)/.
	@chmod 744 $(WDIR)/mix.sh
	@chmod 744 $(WDIR)/cascade.sh
	@chmod 744 $(WDIR)/translator.sh
	@chmod 744 $(WDIR)/addcut.sh
	@chmod 744 $(WDIR)/pdf_reweighter.sh
ifneq (${ROOTSYS},)
	@echo $(ROOTSYS) > $(WDIR)/.rootpath
	@cp usr/rtupler    $(WDIR)/rtupler
	@chmod 744 $(WDIR)/rtupler
endif
	@echo " ";
	@echo "*****************************************************************"
	@echo " Some additional routines have been installed"
ifneq (${ROOTSYS},)
	@echo " Namely: translator, cascade, mix, addcut, pdf_reweighter, rtupler"
else
	@echo " Namely: translator, cascade, mix, addcut, pdf_reweighter"
endif
	@echo "*****************************************************************"

setup: lhapdf
####### Variables test
	@if test "x$(WDIR)" = "x"; then \
	   echo "***** Undefine WDIR in Makefile!"; \
	   exit 99; \
	fi

####### Making directories
	@if (test -d $(WDIR)); then \
	  echo "\n**************************************************************";\
	  echo " The directory WDIR=$(WDIR) exists\n";\
	  echo " Overwrite?"; \
	  echo " Press [Enter] to continue or ^C to stop"; \
	  echo "**************************************************************";\
	  read stdin; \
	  rm -rf $(WDIR);\
	  echo " Old working directory has been deleted"; \
	fi
	@mkdir -p $(WDIR)
	@if (test ! -d $(WDIR)/results); then mkdir -p $(WDIR)/results; fi
	@if (test ! -d $(WDIR)/models); then mkdir -p $(WDIR)/models; fi
	@if (test ! -d $(WDIR)/tmp); then \
	  mkdir -p $(WDIR)/tmp; \
	  cp models/*.mdl $(WDIR)/models/.; \
	fi
	@if (test ! -d $(WDIR)/usr); then mkdir -p $(WDIR)/usr; fi

	@echo $(PWD)           > $(WDIR)/.compheppath
	@cp usr/comphep          $(WDIR)/comphep
	@cp usr/archiv           $(WDIR)/archiv
	@cp usr/comphep.ini      $(WDIR)/comphep.ini
	@cp usr/Makefile         $(WDIR)/usr/Makefile
	@cp usr/mk_tab           $(WDIR)/mk_tab.sh
	@cp usr/tab_view         $(WDIR)/tab_view.sh
	@cp usr/slha.suspect     $(WDIR)/SLHA.sh
	@cp usr/process.dat      $(WDIR)/process.dat
	@cp usr/symb_batch.pl    $(WDIR)/symb_batch.pl
	@cp usr/num_batch.pl     $(WDIR)/num_batch.pl

	@cp usr/Makefile.results $(WDIR)/results/Makefile
	@cp usr/diag_view        $(WDIR)/results/diag_view
	@cp usr/n_comphep        $(WDIR)/results/n_comphep
	@chmod 744 $(WDIR)/archiv
	@chmod 744 $(WDIR)/comphep
	@chmod 744 $(WDIR)/mk_tab.sh
	@chmod 744 $(WDIR)/tab_view.sh
	@chmod 744 $(WDIR)/symb_batch.pl
	@chmod 744 $(WDIR)/num_batch.pl
	@chmod 744 $(WDIR)/SLHA.sh
	@chmod 744 $(WDIR)/results/n_comphep
	@chmod 744 $(WDIR)/results/diag_view

ifneq ($(LHAPDF_DATADIR),)
	@rm -f $(WDIR)/models/beams-lhapdf.mdl $(WDIR)/models/strfuns-lhapdf.mdl
	@echo $(LHAPDF_DATADIR) > $(WDIR)/.lhapdfpath
endif
	@cp usr/userFun* $(WDIR)/usr/.
	@cp usr/userVars* $(WDIR)/usr/.
	$(MAKE) -C $(WDIR)/usr CC='$(CC)' CFLAGS='$(CFLAGS)'

	@echo "\n*****************************************************************"
	@echo " Setup of the user directory for CompHEP-${VERSION} has been"
	@echo " successfully performed"
	@echo " "
	@echo " To start CompHEP go to the user directory ${WDIR}"
	@echo " and enter the command './comphep'"
	@echo "*****************************************************************"

lhapdf:
# LHAPDF 6 is a system library, no symlinks needed

oneclicksetup: lhapdf
	@cp usr/diag_view ../results/diag_view
	@cp usr/n_comphep ../results/n_comphep
	@chmod 744 ../results/n_comphep
	@chmod 744 ../results/diag_view
	@echo $(PWD) > ../.compheppath
ifneq ($(LHAPDF_DATADIR),)
	@echo $(LHAPDF_DATADIR) > ../.lhapdfpath
endif
ifneq (${ROOTSYS},)
	@echo $(ROOTSYS) > ../.rootpath
	@cp usr/rtupler    ../rtupler
	@chmod 744 ../rtupler
endif

######################################################################################
.PHONY: dist bindist

dist: distclean
	@if [ ${COMPHEP_DISTR} != `basename ${PWD}` ]; then \
	  mkdir -p ${CHEPTMPDIR}/${COMPHEP_DISTR}; \
	  cd ${CHEPTMPDIR}; \
	  cp -a ${PWD}/* ${COMPHEP_DISTR}; \
	  mkdir -p ${COMPHEP_DISTR}/doc; \
	  cp -a ${DOC_FILES} ${COMPHEP_DISTR}/doc; \
	  rm -fr `find ${COMPHEP_DISTR} -name CVS`; \
	  rm -fr `find ${COMPHEP_DISTR} -name .svn`; \
	  tar czf ${PWD}/../${COMPHEP_DISTR}.tgz ${COMPHEP_DISTR}; \
	else \
	  cd ..; \
	  tar czf ${COMPHEP_DISTR}.tgz ${COMPHEP_DISTR}; \
	fi
	@cd ${PWD}

bindist:
	@if [ ! -e ${PWD}/bin/s_comphep.exe ]; then \
	  echo " ";\
	  echo "******************************************";\
	  echo " 'make' can not find 's_comphep.exe' file.";\
	  echo "     Please check if you built CompHEP.   ";\
	  echo "******************************************";\
	  echo " ";\
	  exit 99; \
	else \
	  mkdir -p ${CHEPTMPDIR}/${COMPHEP_DISTR}-${OSTYPE}; \
	  cd ${CHEPTMPDIR}; \
	  for x in ${BINDIST_FILES}; do \
		cp -a ${PWD}/$$x ${COMPHEP_DISTR}-${OSTYPE}; \
	  done; \
	  mkdir -p ${COMPHEP_DISTR}-${OSTYPE}/doc; \
	  cp -a ${DOC_FILES} ${COMPHEP_DISTR}-${OSTYPE}/doc; \user
	  cp -a ${DOC_FILES} ${COMPHEP_DISTR}-${OSTYPE}/doc; \
	  mkdir -p ${COMPHEP_DISTR}-${OSTYPE}/src/num; \
	  cp -a ${PWD}/src/num/include ${COMPHEP_DISTR}-${OSTYPE}/src/num; \
	  tar czf ${PWD}/../${COMPHEP_DISTR}-${OSTYPE}.tgz \
	      ${COMPHEP_DISTR}-${OSTYPE}; \
	  rm -fr ${CHEPTMPDIR}/${COMPHEP_DISTR}-${OSTYPE}; \
	  cd ${PWD}; \
	fi
