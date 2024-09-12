/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* -----------------------------------------------------
*/
#ifndef __IS_ISR__
#define __IS_ISR__

#include<stdio.h>
#include "service2/include/files.h"
#include "strfun.h"

extern int p_isr__ (char *p_name__);
extern void n_isr__ (int i, char * beam, char * pdf);
extern int r_isr__ (int i, char *name);
extern int b_isr__ (int i);
extern int m_isr__ (int i, char * p_name);
extern int i_isr__ (int i, double *be, double *mass, char * p_name);
extern double c_isr__ (int i, double x, double q);

extern void info_isr__ (int i, Str_fun_Info * info);

#endif
