/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"

#ifdef LHAPDF

#include "clhapdf.h"
#include "lhapdf.h"

static int curparton[2]={0,0};
static double lambda5 = -0.999;
static int decode[10] = {0,5,4,3,-1,-2,0,2,1,0};

void set_QCDLambda (int beam) {
  lambda5 = lhapdf6_qcdlambda(beam);
}

void delLhapdfList (lhapdfList * list) {
  while (list) {
    lhapdfList *next = list->next;
    if (list->pathfile) free (list->pathfile);
    if (list->file) free (list->file);
    if (list->name) free (list->name);
    if (list) free (list);
    list = next;
  }
}

int comphepLhapdfList (lhapdfList ** list) {
  int n, i;
  int pos = 1;

  n = lhapdf6_num_pdfsets();
  for (i = 0; i < n; i++) {
    const char* setname = lhapdf6_pdfset_name(i);
    if (!setname || !setname[0]) continue;

    {
      lhapdfList * new = malloc (sizeof (lhapdfList));
      new->name = malloc (strlen (setname) + 1);
      strcpy (new->name, setname);
      new->set = 0;
      new->mem = 0;
      new->pathfile = NULL;
      new->file = NULL;
      new->position = pos;
      new->next = *list;
      *list = new;
      pos++;
    }
  }

  return 1;
}


void initLHAPDF (int beamnum, const char* setname, int mem, int prt) {

  curparton[beamnum] = prt;
  if (21 == prt) curparton[beamnum] = 0; /* gluon */
  curparton[beamnum] += 6; /* shift in lhapdfVal according to val[-6:6] (FORTRAN) -> val[0:13] (C)*/

  lhapdf6_initpdf(beamnum, setname, mem);
  set_QCDLambda(beamnum);
}


/* return pdf value of i-th parton in a point (x,Q) 
        -6   -5   -4   -3   -2   -1   0 1 2 3 4 5 6
Parton  tbar bbar cbar sbar ubar dbar g d u s c b t
*/
double lhapdfVal (double x, double q, int i) {
  double val;
  double pdf[13];
  int beam = (i < 2) ? i : 0;

  lhapdf6_evolvepdf (beam, x, q, pdf);
  if (i < 2) {
    val = pdf[curparton[i]]/x;
  } else {
/* CompHEP proton: (b,B) (c,C) (s,S) D U G u d
                    2     3      4   5 6 7 8 9 */
    val = pdf[decode[i] + 6]/x;  /* shift val[-6:6] (FORTRAN) -> val[0:13] (C) */
  }

  return val;
}

double lhapdfValCPYTH (double x, double q, int i) {
  double val = 0.;
  double pdf[13];

  lhapdf6_evolvepdf (0, x, q, pdf);
  if (i == 21) {
    val = pdf[6]/x;
  } else {
    val = pdf[6 + i]/x;
  }

  return val;
}

double TESTlhapdfValCPYTH (double x, double q, int i) {
  int j;
  double val = 0.;
  double pdf[13];

  for (j = 0; j < 100; ++j) {
  x = 0.001 + j * 0.001;

  lhapdf6_evolvepdf (0, x, q, pdf);
  if (i == 21) {
    val = pdf[6]/x;
  } else {
    val = pdf[6 + i]/x;
  }

  fprintf (stdout, "x = %f, pdf = %f\n", x, val);
  }
  exit (-1);

  return val;
}

double lhapdf_interAlpha (double q) {
  return lhapdf6_alphas (0, q);
}

double lhapdf_QCDLambda (void) {
  return lambda5;
}

#endif
