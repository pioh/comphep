/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __CHESS_
#define __CHESS_

#include "service2/include/sets.h"

typedef setofbyte indvertset;	/*  1..3*maxVert */


typedef struct
  {
    byte weit;
    byte g5;
    byte vlnc;
    indvertset ind;
    byte link[2 * maxvert];
  }
vertinfostr;

extern vertinfostr vertinfo[2 * maxvert];
extern int n_vrt;
extern int prgcode[2 * maxvert][2];


extern void makeprgcode (void);

#define MEMORY_OPTIM  0

#endif
