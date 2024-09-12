####################################################
# Copyright (C) 2001-2009, CompHEP Collaboration   #
# Author: A. Sherstnev                             #
####################################################
.PHONY: clean

CC=$(shell cat ../../CC)
CFLAGS=$(shell cat ../../CFLAGS)
INS=-Iinclude

FTERM=data.o atoms.o lexic.o write.o lists.o ops.o parse.o itrio.o
FA1=alias.o alg1.o alg1a.o alg1b.o alg1c.o alg1d.o alg1e.o alg1f.o alg1g.o\
     alg1w.o alg1h.o let.o alg1p.o alg1i.o alg1s.o
FA2=alg2.o alg2a.o alg2b.o alg2c.o alg2r.o alg2d.o alg2e.o lagr.o alg2t.o \
    alg2s.o
FPH1=params.o field.o spec.o group.o  fu1.o util.o gauge.o regprt.o color.o\
     redcol.o red56.o imprt.o redspin.o photon.o srules.o angle.o cls.o
FPH2=herm.o mass.o trigon.o tridef.o triheu.o triuti.o trics3.o fu2.o\
     readch.o wrtfa.o

FALL=$(FTERM) $(FA1) $(FA2) $(FPH1) $(FPH2) func.o varver.o
EXEC=main.o term.o

all: $(FALL) $(EXEC)
	$(AR) -r ../../lib/liblanhep.a $(FALL)

term.o: $(FTERM)
	ld -r $(LOPT) -o term.o $(FTERM)

clean:
	@rm -f *.o

.c.o:
	$(CC) $(CFLAGS) $(INS) -c $<
