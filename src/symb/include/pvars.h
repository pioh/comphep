/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PVARS_
#define __PVARS_

typedef struct varinfo
  {
    unsigned long maxdeg;
    unsigned long zerodeg;
    short num;
    int wordpos;
    char name[7];		/*  used by Editor  */
  }
varinfo;


typedef struct polyvars
  {
    int nvar;			/* <============== Must be initialized to 0  */
    varinfo *vars;
  }
polyvars;

extern polyvars *vardef;

extern void increaseVars (polyvars * v);
extern void clearVars (polyvars * v);
extern void addvar (polyvars * v, char *varname, int deg);
extern void closevars (polyvars * v);
extern void unite_vars (polyvars * result, polyvars * addition);
extern int modelVarPos (char *s);
extern int scalarProductPos (int p1, int p2);

#define ALIG(n) (n>0 ?  8*((n-1)/8+1): 0 )

#endif
