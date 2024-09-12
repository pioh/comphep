####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: all clean fann fann.exe mk_tab.exe addcut.exe

CHEP=$(shell cat ../.compheppath)

CC=$(shell cat $(CHEP)/CC)
CFLAGS=$(shell cat $(CHEP)/CFLAGS)
CLIBS:=$(shell cat $(CHEP)/CLIBS)
INS:=-I$(CHEP)/src

all: mk_tab.exe addcut.exe

fann: fann.exe unweighter.exe

mk_tab.exe: userFun.o
	$(CC) $(CFLAGS) -o mk_tab.exe userFun.o -L$(CHEP)/lib -ltab -lmix -levents -lnum -lserv $(CLIBS)

fann.exe: userVars.o userFun.o
	$(CC) $(CFLAGS) -o fannvars.exe $(CHEP)/src/fann/fannvar.o userVars.o userFun.o -L$(CHEP)/lib -lfann -lmix -levents -lnum -lserv $(CLIBS)

unweighter.exe: nn.o
	$(CC) $(CFLAGS) -o unweighter.exe $(CHEP)/src/fann/unweighter.o $< -L$(CHEP)/lib -lfann -lmix -levents -lnum -lserv $(CLIBS)

addcut.exe: userFun.o
	$(CC) $(CFLAGS) -o addcut.exe userFun.o -L$(CHEP)/lib -ltranls -levents -lserv $(CLIBS)

clean:
	@rm -f *.o *.exe

%.o: %.c
	$(CC) $(CFLAGS) $(INS) -c $<
