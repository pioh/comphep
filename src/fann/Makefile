####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-I../ -Iinclude -I../../include

OBJS=variables.o unweighting_routines.o
EXEC=fannvar.o unweighter.o

all: $(OBJS) $(EXEC)
	$(AR) -r ../../lib/libfann.a $(OBJS)

clean:
	@rm -f *.o

.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
