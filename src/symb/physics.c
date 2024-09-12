/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"

#include "process.h"
#include "process_core.h"
#include "physics.h"

vcsect vcs = {0,0,0,0,0,0,{0},{{{0,0,0,0,{0,0},{0,0}}}}};

whohow exclude_list = {{0,0,0}};
whohow keep_list = {{0,0,0}};

multwhohow exclude_composit_list = {{{0},0,0}};
multwhohow keep_composit_list = {{{0},0,0}};

double sqrts = 0.;
double rapidity = 0.;
unsigned subproc_f = 0;
unsigned subproc_sq = 0;

char modelmenu[STRSIZ] = {0};
int maxmodel = 0;

int nsub = 1;
int ndiagr = 0;
static int n_model = 1;


void nilprtcl (whohow p_list) {
  int i;
  for (i = 0; i < whohowMAX; i++) {
    p_list[i].who = 0;
    p_list[i].how = 0;
    p_list[i].type = 1;
  }
}


void getprtcls (char *txt1, prtclsarray pnames) {
  int ps, j;
  char txt[STRSIZ], pnametmp[STRSIZ];
  int numtot = getntot ();

  strcpy (txt, txt1);
  ps = spos ("->", txt);
  txt[ps - 1] = ',';
  txt[ps + 1 - 1] = ' ';
  sbld (txt, " %s", txt);
  
  for (j = numtot + 1; j <= MAXINOUT; j++)
    strcpy (pnames[j - 1], "***");
  j = 1;
  do
    {
      ps = cpos (',', txt);
      if (ps == 0)
	{
	  trim (txt);
	  strncpy (pnames[j - 1], txt, 3); pnames[j - 1][3] = 0;
	  return;
	}
      strcpy (pnametmp, copy (txt, 1, ps - 1));
      strcpy (txt, copy (txt, ps + 1, (int) strlen (txt) - ps));
      trim (pnametmp);
      strncpy (pnames[j - 1], pnametmp, 3); pnames[j - 1][3] = 0;
      ++(j);
    }
  while (TRUE);
}

int getModelNumberSymb (void) {
  return n_model;
}

void setModelNumberSymb (int num) {
  n_model = num;
}
