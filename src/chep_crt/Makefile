####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS:=-I../ -Iinclude

OBJS=crt.o crt_util.o edittab.o file_scr.o help.o \
     pictures.o tex_util.o X11_crt0.o

all: $(OBJS)
	$(AR) -r ../../lib/libserv.a $(OBJS)

clean:
	@rm -f *.o
	
.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
