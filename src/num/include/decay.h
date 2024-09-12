/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __DECAY__
#define __DECAY__

extern void decay0 (int nvin, double amm1, double amm2, double *factor);
extern void decay1 (int nvpole, double *hsum, double *hdif);
extern void decay2 (int nvout1, double *parcos);
extern void decay3 (int nvath, double parcos, double parfi, int nvout1, int nvout2);
extern void wrmom (int n);

#endif
