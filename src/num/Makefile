####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all n_comphep diag_viewer mix translator oneclick

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)

F77=$(shell cat ../../F77)
F77FLAGS=$(shell cat ../../F77FLAGS)

ROOTFLAGS1=$(shell cat ../../ROOTFLAGS)
INS=-I../ -Iinclude -I../../include

OBJS1=alphas2.o alphas_menu.o cs_22.o cut.o decay.o histogram.o \
      phys_val.o q_kin.o roothisto.o regfunal.o regul.o simpson.o \
      spline.o runVegas.o mc_menu.o rw_sess.o vegas.o width_12.o \
      err_code.o kinaux.o kininpt.o num_serv.o param.o strfun.o \
      rd_num.o subproc.o core_data.o userFactor.o

OBJS2=evnt_format1.o evnt_format2.o evnt_format_lhaef.o \
      event_reader.o evnt_tools.o evnt_menu.o compare.o fill.o \
      validity.o tag_routines.o tag_parser.o tag_reader.o \
      tag_writer.o lhef_routines.o lhaef.o addcut_lhef.o addcut.o \
       pdf_reweight_lhef.o
#      weighted_events.o

OBJS3=mix_cpyth1.o mix_cpyth2.o mix_lhef.o 
OBJS4=trans_cpyth1.o trans_cpyth2.o cascade_cpyth1.o \
      cascade_cpyth2.o cascade_cpyth_lhef.o rtuple_cpyth1.o \
      rtuple_cpyth2.o rtuple_lhef.o

OBJS5=lhapdf.o pdf.o sf_epa.o sf_isr.o sf_lsr.o sf_pdf.o sf_lhapdf.o \
      strfun_par.o tools.o

OBJS6=configurator.o rw_batch.o

EXEC=diag_viewer.o diag_viewer2.o n_comphep.o mix.o translator.o \
      cascade.o rtuple.o pdf_reweighter.o

cobjects=$(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(EXEC)
cppobjects=rtuple_routines.o
fobjects=wrap.o 

all: n_comphep mix translator events oneclick pdf $(EXEC)

n_comphep: $(OBJS1) $(OBJS5)
	$(AR) -r ../../lib/libnum.a $(OBJS1)

events: $(OBJS2)
	$(AR) -r ../../lib/libevents.a $(OBJS2)

mix: $(OBJS3)
	$(AR) -r ../../lib/libmix.a $(OBJS3)

translator: $(OBJS4) $(cppobjects)
	$(AR) -r ../../lib/libtranls.a $(OBJS4) $(cppobjects)

pdf: $(OBJS5) $(fobjects)
	$(AR) -r ../../lib/libpdf.a $(OBJS5) $(fobjects)

oneclick: $(OBJS6)
	$(AR) -r ../../lib/liboneclick.a $(OBJS6)

clean:
	@rm -f *.o

$(cobjects): %.o: %.c
	$(CC) -c $(CFLAGS) $(INS) $< -o $@

$(cppobjects): %.o: %.cxx
	$(CXX) $(CFLAGS) $(ROOTFLAGS1) $(INS) -c $< -o $@

$(fobjects): %.o: %.F
	$(F77) $(F77FLAGS) $(INS) -c $<
