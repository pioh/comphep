/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __SF_PDF__
#define __SF_PDF__

#include "strfun.h"

extern int p_pdf (char *p_name);
extern void n_pdf (int i, char * beam, char * pdf);
extern int r_pdf (int i, char *name);
extern int beam_pdf (int i);
extern int pdf_pdf (int i, char * p_name);
extern int be_pdf (int i, double *be, double *mass, char * p_name);
extern double c_pdf (int i, double x, double q);
extern double alpha_pdf (double q);

extern void info_pdf (int i, Str_fun_Info * info);
extern int pdf_namecmp (void);

#endif
