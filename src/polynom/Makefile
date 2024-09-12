####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: all clean

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-I../ -Iinclude

OBJS=polynom.o spinor.o tensor.o

all: $(OBJS)
	$(AR) -r ../../lib/libsymb.a $(OBJS)

clean:
	@rm -f *.o

.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
