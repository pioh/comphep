/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __SF_LSR__
#define __SF_LSR__

#include "service2/include/files.h"
#include "strfun.h"

extern int p_lsr__ (char *p_name__);
extern void n_lsr__ (int i, char * beam, char * pdf);
extern int r_lsr__ (int i, char *name);
extern int b_lsr__ (int i);
extern int m_lsr__ (int i, char * p_name);
extern int i_lsr__ (int i, double *be, double *mass, char * p_name);
extern double c_lsr__ (int i, double x, double q);

extern void info_lsr__ (int i, Str_fun_Info * info);

#endif
