/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __OPTIMISE_
#define __OPTIMISE_

#include "polynom/include/lnum.h"

#define infoptr struct inforec *
#define varptr struct var_rec *
typedef struct var_rec
  {
    varptr next;
    char sgn;
    infoptr coef;
    short vars[STRSIZ];
  }
var_rec;
#undef varptr
#undef infoptr
typedef struct var_rec *varptr;


#define infoptr struct inforec *
typedef struct inforec
  {
    infoptr next;
    unsigned count;
    char name[20];
    enum
      {
	numb, expr, rnumb
      }
    consttype;
    varptr const_;
    NUM_TYPE ival;
    double rval;
  }
inforec;
#undef infoptr
typedef struct inforec *infoptr;


typedef pointer (*smplemit) (varptr ex);
typedef pointer (*vfact) (int ch, int deg, pointer pmult, pointer psum);
typedef pointer (*cfact) (infoptr c, pointer pmult, pointer psum);

extern void initinfo (void);
extern void readpolynom (FILE * fres, varptr * expr_);
extern pointer emitexpr (varptr ex, smplemit smplemitfun,
			 vfact vfactfun, cfact cfactfun);
extern int equalexpr (varptr v1, varptr v2);

#define minvarrec (sizeof(struct var_rec) - sizeof(short)*(STRSIZ-1))

extern infoptr info;

extern int firstVar;

extern int short_strlen (short *);
extern void short_strcpy (short *, short *);
extern int short_strcmp (short *, short *);
extern short *short_strchr (short *, short);
#endif
