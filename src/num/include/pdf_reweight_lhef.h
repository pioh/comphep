/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __PDF_REWEIGHT_LHEF__
#define __PDF_REWEIGHT_LHEF__

extern int read_config (char configname[]);
extern int pdf_reweight_lhef (char config[], char source[], char target[], int nalphas, long the_seed, int error_used);

#endif
