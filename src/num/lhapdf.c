/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"

#ifdef LHAPDF

#include "clhapdf.h"
#include "lhapdf.h"

static int curparton[2]={0,0};
static double lambda5 = -0.999;
static int decode[10] = {0,5,4,3,-1,-2,0,2,1,0};

void makelhapdfsilent_(void);

void set_QCDLambda (int set, int mem) {
  lambda5 = lhapdfqcdlam(set, mem);
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

int comphepLhapdfList (char * thepath, char * indexfile, lhapdfList ** list) {
  midstr s;
  midstr file1;
  int pos = 1;
  int set;
  int mem;
  int pdftyp;
  int pdfgup;
  int pdfsup;
  double q2min;
  double q2max;
  double xmin;
  double xmax;

  FILE * f = fopen (indexfile, "r");
  if (!f) {
    return 0;
  }

  while (fgets (s, 1024, f)) {
    int err = sscanf (s," %d %d %d %d %s %d %lf %lf %lf %lf", 
            &set, &pdftyp, &pdfgup, &pdfsup, file1, &mem, &q2min, &q2max, &xmin, &xmax);

    if (10 == err && 0 < strlen (file1)) {
      lhapdfList * new = malloc (sizeof (lhapdfList));

      new->name = malloc (strlen (file1) + 1);
      new->set = set;
      new->mem = mem;
      new->pathfile = malloc (strlen (thepath) + 1);
      new->file = malloc (strlen (file1) + 1);

      strcpy (new->name, file1);
      strcpy (new->pathfile, thepath);
      strcpy (new->file, file1);
      new->position = pos;

      new->next = *list;
      *list = new;
      pos++;
    }
  }
  fclose (f);
  return 1;
}


int setLHAPDFIndexpath (longstr path) {
  FILE * f;
  longstr tmppath;
  int status = 0;

    sprintf (tmppath, "PDFsets.index");
  f = fopen (tmppath, "r");
  if (f) {
    strcpy (path, tmppath);
    fclose(f);
    status = 1;
  } else {
    sprintf (tmppath, "%s/strfun/PDFsets.index", pathtocomphep);
    f = fopen (tmppath, "r");
    if (f) {
      strcpy (path, tmppath);
      fclose(f);
      status = 1;
    } else {
      sprintf (tmppath, "%s/share/lhapdf/PDFsets.index", pathtolhapdf);
      f = fopen (tmppath, "r");
      if (f) {
        strcpy (path, tmppath);
        fclose(f);
        status = 1;
      } else {
        fprintf (stderr,"Warning! can not find the LHAindex-comphep.txt file\n");
      }
    }
  }
  return status;
}


int setLHAPDFIndexpathCPYTH (longstr path) {
  FILE * f;
  longstr tmppath;
  int status = 0;

  sprintf (tmppath, "PDFsets.index");
  f = fopen (tmppath, "r");
  if (f) {
    strcpy (path, tmppath);
    fclose(f);
    status = 1;
  } else {
    sprintf (tmppath, "%s/share/lhapdf/PDFsets.index", pathtolhapdf);
    f = fopen (tmppath, "r");
    if (f) {
      strcpy (path, tmppath);
      fclose(f);
      status = 1;
    } else {
      fprintf (stderr,"Warning! can not find the PDFsets.index file\n");
    }
  }
  return status;
}


void initLHAPDF (int beamnum, char * path, char * file, int set, int mem, int prt) {
  int len;
  longstr pathfile;

  curparton[beamnum] = prt;
  if (21 == prt) curparton[beamnum] = 0; /* gluon */
  curparton[beamnum] += 6; /* shift in lhapdfVal according to val[-6:6] (FORTRAN) -> val[0:13] (C)*/

  sprintf (pathfile, "%s/%s", path, file);
  len = strlen (pathfile);
  makelhapdfsilent_();
  initpdfset(pathfile);
  initpdf(mem);
  set_QCDLambda(set, mem);
}

/* return pdf value of i-th parton in a point (x,Q) 
        -6   -5   -4   -3   -2   -1   0 1 2 3 4 5 6
Parton  tbar bbar cbar sbar ubar dbar g d u s c b t
*/
double lhapdfVal (double x, double q, int i) {
  double val;
  double pdf[13];

  evolvepdf (x,q,pdf);
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

  locevolvepdf (x, q, pdf);
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

  locevolvepdf (x, q, pdf);
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
  double res = lhapdfalphas (q);
  return res;
}

double lhapdf_QCDLambda (void) {
  return lambda5;
}

#endif
