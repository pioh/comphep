/*
* Copyright (C) 2026, CompHEP Collaboration
* ------------------------------------------------
* C++ shim for LHAPDF 6 API.
* Provides extern "C" functions callable from C code.
*/
/* -DLHAPDF in CFLAGS defines LHAPDF as 1, conflicting with namespace LHAPDF:: */
#ifdef LHAPDF
#undef LHAPDF
#endif
#include "LHAPDF/LHAPDF.h"
#include <cstring>
#include <cstdlib>
#include <vector>
#include <string>

static LHAPDF::PDF* currentPDF[2] = {nullptr, nullptr};

extern "C" {

/* Initialize PDF for beam (0 or 1) */
void lhapdf6_initpdf(int beam, const char* setname, int member) {
    if (currentPDF[beam]) delete currentPDF[beam];
    currentPDF[beam] = LHAPDF::mkPDF(setname, member);
}

/* Compute xf(x,Q) -- fills 13-element array pdf[0..12]
   index: 0=tbar 1=bbar 2=cbar 3=sbar 4=ubar 5=dbar 6=g 7=d 8=u 9=s 10=c 11=b 12=t */
void lhapdf6_evolvepdf(int beam, double x, double Q, double* pdf) {
    if (!currentPDF[beam]) return;
    for (int i = -6; i <= 6; i++) {
        pdf[i+6] = currentPDF[beam]->xfxQ(i, x, Q);
    }
}

/* alpha_s(Q) */
double lhapdf6_alphas(int beam, double Q) {
    if (!currentPDF[beam]) return 0.0;
    return currentPDF[beam]->alphasQ(Q);
}

/* Lambda QCD (Lambda5) from PDF info */
double lhapdf6_qcdlambda(int beam) {
    if (!currentPDF[beam]) return 0.0;
    try {
        return currentPDF[beam]->info().get_entry_as<double>("AlphaS_Lambda5");
    } catch (...) {
        return 0.0;
    }
}

/* QCD order of alpha_s from PDF info */
int lhapdf6_qcdorder(int beam) {
    if (!currentPDF[beam]) return 0;
    try {
        return currentPDF[beam]->info().get_entry_as<int>("AlphaS_OrderQCD");
    } catch (...) {
        return 0;
    }
}

/* Number of available PDF sets */
int lhapdf6_num_pdfsets(void) {
    return (int)LHAPDF::availablePDFSets().size();
}

/* Get name of i-th PDF set */
const char* lhapdf6_pdfset_name(int i) {
    static std::string name;
    const auto& sets = LHAPDF::availablePDFSets();
    if (i < 0 || i >= (int)sets.size()) return "";
    name = sets[i];
    return name.c_str();
}

/* Number of members in a PDF set */
int lhapdf6_num_members(const char* setname) {
    try {
        LHAPDF::PDFSet& pdfset = LHAPDF::getPDFSet(setname);
        return (int)pdfset.size();
    } catch (...) {
        return 0;
    }
}

/* Cleanup */
void lhapdf6_cleanup(void) {
    for (int i = 0; i < 2; i++) {
        if (currentPDF[i]) { delete currentPDF[i]; currentPDF[i] = nullptr; }
    }
}

} /* extern "C" */
