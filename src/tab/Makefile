####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-I../ -Iinclude -I../../include

OBJS=e_tools.o mk_tab.o tab_routines.o
all: $(OBJS)
	$(AR) -r ../../lib/libtab.a $(OBJS)

clean:
	@rm -f *.o

.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
