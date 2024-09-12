/*
* Copyright (C) 2001-2009, CompHEP Collaboration
*---------------------------------------------------
*/
#ifndef __SAVERES_
#define __SAVERES_

#include "pvars.h"
#include "physics.h"
#include "polynom/include/polynom.h"
#include "denominators.h"


extern denom_struct denom[2 * maxvert - 2];

extern byte denrno;

extern void save_analitic_results (FILE * fres, poly rnum, poly factn, poly factd, polyvars * vars1, polyvars * vars2, vcsect vcs,
 int status_key);
extern void save_empty_analitic_results (FILE * fres, vcsect vcs);

#endif
