/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef _E_TOOLS_
#define _E_TOOLS_

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

extern int getNames (char * flow);
extern double getTotCS (void);
extern double getCS (void);
extern int getNevents (void);

extern int old_getNames (FILE * flow);
extern int getMasses (FILE * flow);
extern int skipHeadLine (FILE * flow);
extern int skipLine (FILE * flow);
extern int readCS (FILE * flow);
extern int readTotCS (FILE * flow);
extern int readNEvents (FILE * flow);

extern void boost (double *n, double *p);
extern void findBoost (double *p, double *n);
extern int decay2 (double M, double *p1, double *p2);
#endif
