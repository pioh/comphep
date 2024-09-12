/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include<math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "chep_crt/include/chep_crt.h"

#include "model.h"
#include "process_core.h"

static int ninprt;         /* Number of ininial particles */
static int noutprt;        /* Number of ininial particles */
static int ntotprt;        /* Number of ininial particles */
static int n_xprt;         /* Number of ininial particles */

static double sqrts;
static double rapidity;

static shortstr finalstatech;   /* stored (in safe.dat) string with process */
static shortstr processch;      /* stored (in safe.dat) string with process */

static shortstr excl_str;       /* stored (in safe.dat) string with diagrams excluding conditions */
static shortstr keep_str;       /* stored (in safe.dat) string with diagrams keeping conditions */

//hadron hadrons[MAXINOUT];

static int NcInfLimit = 0;

int getntot (void) {return ntotprt;}
int getnin (void)  {return ninprt;}
int getnout (void) {return noutprt;}
int getnx (void) {return n_xprt;}
double getsqrtS (void) { return sqrts;}
double getRapidity (void) { return rapidity;}

void setnin (int num) {
  ninprt = num;
  ntotprt = ninprt + noutprt;
}

void setnout (int num) {
  noutprt = num;
  ntotprt = ninprt + noutprt;
}

void setntot (int num) {
  ntotprt = num;
}

void setnx (int num) {
  n_xprt = num;
}

void setsqrtS (double val) {
  sqrts = val;
}

void setRapidity (double val) {
  rapidity = val;
}

extern void set_sqrts_rap (double e1, double p1, double e2, double p2) {
  sqrts = sqrt ((e1 + e2) * (e1 + e2) - (p1 - p2) * (p1 - p2));
  rapidity = atanh ((p1 - p2) / (e1 + e2));
}

char * getFinalstatech (void) {return finalstatech;}
char * getProcessch (void) {return processch;}
char * getExclprtlist (void) {return excl_str;}
char * getKeepprtlist (void) {return keep_str;}

void setFinalstatech (shortstr s) {
  strcpy(finalstatech, s);
}

void setProcessch (shortstr s) {
  strcpy(processch, s);
}

void setExclprtlist (shortstr s) {
  strcpy(excl_str, s);
}

void setKeepprtlist (shortstr s) {
  strcpy(keep_str, s);
}

int getNcinflimit (void) {
  return NcInfLimit;
}

void setNcinflimit (int num) {
  NcInfLimit = num;
}
