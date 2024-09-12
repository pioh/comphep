/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PREPDIAG_
#define __PREPDIAG_

#include "service2/include/sets.h"
#include "ghosts.h"

typedef setofbyte indset;	/*  1..3*maxVert   */

typedef char momsum[MAXINOUT + 1];

typedef struct
  {
    int g5;
    int lprtcl;
    byte len;
    byte vv[2 * maxvert];
    byte ll[2 * maxvert];
    byte intln[2 * maxvert];
    byte intln2[2 * maxvert];
    char invrt[2 * maxvert];
  }
fermloopstp;


typedef struct
  {
    int r_vert;
    arr4byte subst;
    algvertptr lgrnptr;
  }
vertexhlp;

typedef struct
  {
    byte vrt1;
    byte ln1;
    byte vrt2;
    byte ln2;
  }
linkhlp;

extern vertexhlp vertexes[2 * maxvert];
extern linkhlp massindpos[5 * maxvert];
extern fermloopstp fermloops[maxvert];
extern indset setmassindex;
extern indset setmassindex0;
extern int nloop;
extern int consLow;
extern int fermmap[2 * maxvert];
extern char inoutmasses[MAXINOUT][7];
extern momsum momdep[3 * maxvert];

extern void preperdiagram (void);
extern void coloringvcs (hlpcsptr currentghst);
extern void attachvertexes (void);
extern void findReversVert (void);

#endif
