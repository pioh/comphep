####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean all

F77=$(shell cat ../F77)
F77FLAGS=$(shell cat ../F77FLAGS)
CC=$(shell cat ../CC)
CFLAGS=$(shell cat ../CFLAGS)

TC_OBJS=rotate_f.o g2f3_f.o

OBJS= $(TC_OBJS) bsg.o suspect2.o twoloophiggs.o \
FeynHiggsFast.o FeynHiggsFastCodes.o FeynHiggsFastSub.o delrhoFastSub.o

EXEC=diag_viewer.o n_comphep.o 

all: suspect2_call.o compile

compile: $(OBJS)
	@$(AR) -r ../lib/libexternal.a $(OBJS)

clean:
	@rm -f *.o

.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<

.f.o:
	$(F77) $(F77FLAGS) $(INS) -c $<
