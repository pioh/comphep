####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-I../ -Iinclude -I$../../include

OBJS=files.o paragraphs.o parser.o read_func.o sets.o tptcmac.o \
     drandXX.o lbl.o syst.o parseref.o f_c.o \
     unix_utils.o getmem.o 4_vector.o kfcodes.o width_func.o \
     lhapdf_update.o

all: $(OBJS) 
	$(AR) -r ../../lib/libserv.a $(OBJS)

clean:
	@rm -f *.o

.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
