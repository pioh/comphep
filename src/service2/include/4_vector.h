/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#include<stdio.h>

#ifndef __4_VERTOR__
#define __4_VERTOR__
extern double pvect[400];

extern double vsqrt (double a);
extern double vdot4 (int i, int j);
extern double mom4mod (int i);
extern void vsum4 (int i, int j, int k, int isg);
extern void vnull4 (int i);
extern void eps4 (int n1, int n2, int n3, int n4);
extern void pvFill (double mass, double *mom, int pos);
extern void lvtonv (char *lv, int nin, int nvpos);
extern void lorrot (double rapidity, int ntot);
extern void lorenc (double *mom, double *p_frame, double *p_trans);
extern void new_lorenc (double *mom, double *p_frame, double *newmom);

#endif
