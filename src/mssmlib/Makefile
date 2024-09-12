####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-I../ -Iinclude

OBJS1= fhf3.o sugrac.o
OBJS2= rotate_c.o g2f3_c.o

OBJS= $(OBJS1) $(OBJS2)

all: $(OBJS)
	$(AR) -r ../../lib/libmssm.a $(OBJS)

clean:
	@rm -f *.o
	
.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
