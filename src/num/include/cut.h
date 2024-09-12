/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __CUT__
#define __CUT__

#include<stdio.h>
#include "chep_crt/include/file_scr.h"
#include "kinaux.h"

extern int fillCutArray (void);
extern int ReadCuts (FILE *);
extern int WriteCuts (FILE *);

extern double calcCutFactor (void);
extern int printDetailsOfCutFactor (void);
extern int get_cutn (void);

extern int rancor (double *vmin, double *vmax, double shift, double fmult, int n);

extern table cutTab;

typedef struct
  {
    char key;
    char lvinvc[PLISTLEN];
    int minon;
    int maxon;
    int exclusive;
    double cvmin;
    double cvmax;
  }
invcut_;

extern invcut_ invcut_1[64];

#endif
