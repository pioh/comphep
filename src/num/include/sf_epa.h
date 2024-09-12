/*
* Copyright (C) 2001-2009, CompHEP Collaboration
*--------------------------------------------------------
*/
#ifndef __SF_EPA__
#define __SF_EPA__

#include "service2/include/files.h"
#include "strfun.h"

extern int p_epa__ (char *p_name__);
extern void n_epa__ (int i, char * beam, char * pdf);
extern int r_epa__ (int i, char *name);
extern int b_epa__ (int i);
extern int m_epa__ (int i, char * p_name);
extern int i_epa__ (int i, double *be, double *mass, char * p_name);
extern double c_epa__ (int i, double x, double q);

extern void info_epa__ (int i, Str_fun_Info * info);
#endif
