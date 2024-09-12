/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __SPLINE__
#define __SPLINE__

extern void progonca (int n, double *f, double *b, double *c, double *d);
extern double spline_for_graph (double iks, double *f, double *b, double *c, double *d);
extern void SPLINE (double koeff, int N, double *Iexp, double *dIexp, double *f, double *b, double *c, double *d);

#endif
