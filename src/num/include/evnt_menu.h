/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __EVNT_MENU__
#define __EVNT_MENU__

#include "vegas.h"
extern void menu_EventGenerator (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, FILE * iprt, int init);
extern void menu_1to2_EventGenerator (double (*func) (double *, double), char * fname, FILE * iprt);

extern int ClearEventMax (void);
extern int WriteEventMax (FILE * f);
extern int ReadEventMax (FILE * f);

int WriteEventSettings (FILE * f);
int ReadEventSettings (FILE * f);
#endif
