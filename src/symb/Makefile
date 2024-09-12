####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-I../ -Iinclude
     
OBJS=chess.o colorf.o cweight.o denominators.o ghosts.o draw_ampl.o \
     model.o m_utils.o physics.o prepdiag.o pre_read.o process.o \
     process_core.o pvars.o r_code.o reader0.o reader_s.o read_mdl.o \
     rfactor.o saveres.o beams.o screen.o sos.o symbolic.o test_wrt.o \
     ini_ses.o combine.o form_code.o l_string.o optimise.o procvar.o \
     out_c.o out_form.o out_math.o out_service.o reader_c.o out_red.o \
     constr.o diagrams.o diaprins.o showgrph.o squar.o batch.o

EXEC=s_comphep.o archiv.o

all: $(OBJS) $(EXEC)
	$(AR) -r ../../lib/libsymb.a $(OBJS)

clean:
	@rm -f *.o

.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
