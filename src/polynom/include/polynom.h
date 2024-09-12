/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __POLYNOM_
#define __POLYNOM_

#include"lnum.h"


#define poly struct monom *
typedef struct monom
  {
    poly next;
    union
      {
	NUM_TYPE num;
	struct
	  {
	    poly re;
	    poly im;
	  }
	complex;
	int type;
      }
    coef;
    union
      {
	unsigned long power[2];
	char tens[2 * sizeof (long)];
	struct
	  {
	    char l;
	    char g5;
	    char g[2 * sizeof (long) - 2];
	  }
	spin;
      }
    tail;
  }
monom;
#undef poly
typedef struct monom *poly;

extern void set_garbage (poly g);
extern poly get_garbage (void);

extern int levi;
extern poly *contracts;
extern int maxLength, monomLength;


extern void newmonom (poly * p);
extern void delpoly (poly * p);
extern void delmonom (poly * p);
extern poly plusone (void);
extern poly copypoly (poly p);
extern void sewpoly (poly * p1, poly * p2);
extern void multpolyint (poly * p, long i);
extern poly multtwopoly (poly q1, poly q2);
extern poly scalarmult (int p1, int p2);
extern void assignsclmult (int p1, int p2, poly p);
extern void deltensor (poly * t);
extern poly copytens (poly t, int ln);
extern void sewtens (poly * t1, poly * t2, int ln);
extern void multtensint (poly * t, long i);
extern void multtensComplexpoly (poly * t, poly re, poly im);
extern void multtenspoly (poly * t, poly p);
extern void (*memoryInfo) (int);
extern void tensRealPart (poly * t);
#endif
