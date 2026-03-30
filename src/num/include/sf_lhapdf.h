/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __SF_LHAPDF__
#define __SF_LHAPDF__

#include "strfun.h"

extern int p_lhapdf (char * p_name);
extern void n_lhapdf (int i, char * beam, char * pdf);
extern int r_lhapdf (int i, char * name);
extern int beam_lhapdf (int i);
extern int pdf_lhapdf (int i, char * p_name);
extern int be_lhapdf (int i, double * be, double * mass, char * p_name);
extern double c_lhapdf (int i, double x, double q);
extern double alpha_lhapdf (double q);

extern void info_lhapdf (int i, Str_fun_Info * info);
extern int getLHAPDFset (int i);
extern int getLHAPDFmem (int i);
extern int lhapdf_namecmp (void);

#endif
