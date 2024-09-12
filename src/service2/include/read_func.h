/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __READER_FUNC__
#define __READER_FUNC__

extern double sort4 (double m1, double m2, double m3, double m4, double dn);
extern int calcExpression (char *s, int (*nameToDouble) (char *, double *), double *p);
extern int rd_num (char *s, double *p);

#endif
