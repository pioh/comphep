/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __ALPHAS2__
#define __ALPHAS2__

extern double alpha_2 (double dscale);

extern void recalc_alphas (void);
extern void setLambda6 (double lam);

extern double QCDLambda (void);
extern int Nflavour (void);
extern int QCDOrder (void);

extern double alpha_em (double dscale);
#endif
