/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __KINAUX__
#define __KINAUX__

extern void coninv_ (char *lv);
extern int eqvect_ (char *lv1, char *lv2);
extern int spole_ (char *lv);
extern int printlv_ (char *lv);
extern void sngpos_ (char *lv, int *ndec, int *nclust, char *lvaux);
extern void lvmirr_ (char *lv);

typedef struct
  {
    char lvin[PLISTLEN];
    char lvout[2][PLISTLEN];
  }
kinmtc_;

extern kinmtc_ kinmtc_1[PLISTLEN];

#endif
