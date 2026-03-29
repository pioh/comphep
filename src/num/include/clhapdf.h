/*
 * Copyright (C) 2007, Alexander Sherstnev
 * Copyright (C) 2001-2007, CompHEP Collaboration
 *------------------------------------------------------
 * C prototypes for LHAPDF 6 wrapper (lhapdf6_wrapper.cpp)
*/
#ifndef __SF_Clhapdf_
#define __SF_Clhapdf_

void lhapdf6_initpdf(int beam, const char* setname, int member);
void lhapdf6_evolvepdf(int beam, double x, double Q, double* pdf);
double lhapdf6_alphas(int beam, double Q);
double lhapdf6_qcdlambda(int beam);
int lhapdf6_qcdorder(int beam);
int lhapdf6_num_pdfsets(void);
const char* lhapdf6_pdfset_name(int i);
int lhapdf6_num_members(const char* setname);
void lhapdf6_cleanup(void);

#endif
