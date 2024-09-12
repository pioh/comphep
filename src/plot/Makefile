####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-I../ -Iinclude

OBJS=plot.o view_tab.o

all: $(OBJS) 
	$(AR) -r ../../lib/libnum.a plot.o 

clean:
	@rm -f *.o
	
.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
