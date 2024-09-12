/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef ROOTused

#include "service2/include/chep_limits.h"
#include "service2/include/kfcodes.h"
#include "service2/include/4_vector.h"

#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_routines.h"
#include "alphas2.h"
#include "strfun.h"

#include "LesHouches.h"
#include "rtuple_routines.h"

#include "rtuple_cpyth2.h"

static int nin_ = 0;
static int nout_ = 0;
static char pname[MAXNP][20];

static int pbeam[2];
static double ebeam[2];
static int PDFLIBgroup[2];
static int PDFLIBset[2];

static int 
getNames (char * buff)
{
  char * buf;
  char * pch;
  char * temp = malloc (strlen (buff) * sizeof (char));

  nin_ = 0;
  strcpy (temp, buff);
  pch = strtok (buff, " ,->");
  while (pch != NULL) {
    strcpy (pname[nin_], pch);
    pch = strtok (NULL, " ,->");
    ++nin_;
  }

  if (0 == nin_) {
    fprintf (stderr, "mk_tab (error): strange (sub)process name, can't extract parton names and set nin_/nout_\n");
    return -1;
  }

  nout_ = 0;
  buf = strstr(temp, "->") + 3;
  pch = strtok (buf, " ,");
  while (pch != NULL) {
    ++nout_;
    pch = strtok (NULL, " ,");
  }
  nin_ -= nout_;

  return 0;
}

int
rtuple_cpyth2 (char ini_name[], char out_name[], int nevnt)
{
  /* New XML-style fromat */
  int i, j, k;
  int n;
  int nntot;
  int format;
  int itag;
  int mpar;
  int n1, n2, num;
  int Nproc;
  int Nevents = 0;
  int err;

  int * pcolor = NULL;
  int * pscale = NULL;
  int * nin = NULL;
  int * nshft = NULL;
  int * ntot = NULL;
  int * id = NULL;
  int * pnum = NULL;

  midstr buff;
  midstr buf;
  midstr name_proc;

  double cs, cserr;
  double qcdscale;
  char bufff[MAXINOUT + 3][STRSIZ];
  double * pcs = NULL;
  double * pcserr = NULL;
  double * mass = NULL;

  string_comnd com;
  tags * head = init_cap (1);
  FILE * s = fopen (ini_name, "r");

  cup_reader (s, head);
  strcpy (com.name, "Nproc");
  get_tag_with1com (0, head, "total", &com);
  Nproc = atoi (com.value);

  nin      = malloc (Nproc * sizeof (int));
  ntot     = malloc (Nproc * sizeof (int));
  pcs      = malloc (Nproc * sizeof (double));
  pcserr   = malloc (Nproc * sizeof (double));

  /* check file format */
  if (check_cpyth2 (head)) {
    return -1.0;
  }

  /* Get number of events and total CS from the total tag */
  n = get_tag (0, head, "total");
  cs = get_fval (0, "CrosSec", head->tag[n]);
  cserr = get_fval (0, "CrosSecErr", head->tag[n]);

  for (i = 0; i < Nproc; i++) {
    sprintf (com.value, "%i", i + 1);

  /* Get number of partons, parton's names and nin_/nout_ from the process tag */
    strcpy (com.name, "ID");
    n = get_tag_with_exactcom (0, head, "process", com);
    pcs[i]      = get_fval (0, "CrosSec", head->tag[n]);
    pcserr[i]   = get_fval (0, "CrosSecErr", head->tag[n]);
    strcpy (name_proc, get_cval (0, "name", head->tag[n]));
    getNames (name_proc);
    nin[i]  = nin_;
    ntot[i] = nin_ + nout_;
  }

  nshft = malloc (Nproc * sizeof (int));
  nshft[0] = 0;
  for (i = 1; i < Nproc; ++i) {
    nshft[i] = nshft[i - 1] + ntot[i - 1];
  }

  nntot  = nshft[Nproc - 1] + ntot[Nproc - 1];
  mass   = malloc (nntot * sizeof (double));
  id     = malloc (nntot * sizeof (double));
  pnum   = malloc (4 * nntot * sizeof (int));
  pcolor = malloc (Nproc * sizeof (int));
  pscale = malloc (Nproc * sizeof (int));
  for (i = 0; i < Nproc; i++) {
  /* Get partons masses from parton tags */
    strcpy (com.name, "ID");
    sprintf (com.value, "%i", i + 1);
    n = get_tag_with_exactcom (0, head, "process", com);
    strcpy (name_proc, get_cval (0, "name", head->tag[n]));
    getNames (name_proc);
    strcpy (com.name, "IDprocess");
    sprintf (com.value, "%i", i + 1);
    n = -1;
    for (j = 0; j < ntot[i]; j++) {
      n = get_tag_with_exactcom (n + 1, head, "parton", com);
      if (get_ival (0, "in", head->tag[n]) || get_ival (0, "out", head->tag[n]))
        mass[nshft[i] + j] = get_fval (0, "mass", head->tag[n]);
      id[nshft[i] + j] = kfpart (pname[j]);
    }

  /* Construction of event format */
    strcpy (com.name, "IDprocess");
    sprintf (com.value, "%i", i + 1);
    format = get_tag_with_exactcom (0, head, "format", com);
    strcpy (com.name, "p1.3");
    pnum[4 * nshft[i] + 3] = tag_contain_com (&com, head->tag[format]) - 2;
    strcpy (com.name, "p2.3");
    pnum[4 * nshft[i] + 7] = tag_contain_com (&com, head->tag[format]) - 2;
    strcpy (com.name, "color_chain");
    pcolor[i] = tag_contain_com (&com, head->tag[format]) - 2;
    strcpy (com.name, "Qsquared");
    pscale[i] = tag_contain_com (&com, head->tag[format]) - 2;
    for (j = 2; j < ntot[i]; ++j) {
      for (k = 1; k < 4; ++k) {
        sprintf (com.name, "p%i.%i", j + 1, k);
        n = tag_contain_com (&com, head->tag[format]);
        pnum[4 * nshft[i] + 4 * j + k] = n - 2;
      }
    }
  }

  n1 = get_tag ( 0, head, "beam");
  n2 = get_tag (n1 + 1, head, "beam");
  num = get_ival (0, "ID", head->tag[n1]);
  pbeam[num - 1] = get_ival (0, "KF", head->tag[n1]);
  ebeam[num - 1] = get_fval (0, "energy", head->tag[n1]);
  num = get_ival (0, "ID", head->tag[n2]);
  pbeam[num - 1] = get_ival (0, "KF", head->tag[n2]);
  ebeam[num - 1] = get_fval (0, "energy", head->tag[n2]);

  n1 = get_tag ( 0, head, "strfun");
  n2 = get_tag (n1 + 1, head, "strfun");
  num = get_ival (0, "IDbeam", head->tag[n1]);
  PDFLIBset[num - 1]   = get_ival (0, "PDFid", head->tag[n1]);
  PDFLIBgroup[num - 1] = get_ival (0, "PDFgr", head->tag[n1]);
  num = get_ival (0, "IDbeam", head->tag[n2]);
  PDFLIBset[num - 1]   = get_ival (0, "PDFid", head->tag[n2]);
  PDFLIBgroup[num - 1] = get_ival (0, "PDFgr", head->tag[n2]);
  if (0 == PDFLIBgroup[0] || 0 == PDFLIBgroup[1] || 0 == PDFLIBset[0] || 0 == PDFLIBset[1]) {
    num = get_ival (0, "IDbeam", head->tag[n1]);
    PDFLIBset[num - 1]   = get_ival (0, "LHAPDFmem", head->tag[n1]);
    PDFLIBgroup[num - 1] = get_ival (0, "LHAPDFid", head->tag[n1]);
    num = get_ival (0, "IDbeam", head->tag[n2]);
    PDFLIBset[num - 1]   = get_ival (0, "LHAPDFmem", head->tag[n2]);
    PDFLIBgroup[num - 1] = get_ival (0, "LHAPDFid", head->tag[n2]);
  }

  err = book_rtuple (out_name);
  if (err) {
   fprintf (stderr, "rtupler (error): can't create a rtuple-file. exit...\n");
   return -1;
  }

  fgets (buf, 2048, s);
  while (fgets (buf, 2048, s)) {
    int npr;
    int IPSGN = 1;
    int col[MAXINOUT][2];
    eventUP ev;

    sscanf (buf, "%i:%[^\n]", &npr, buff);
    --npr;
    for (i = 0; i < 8; i++) pvect[i] = 0.;
    i = 0;
    char * pch = strtok (buff, ":");
    while (pch != NULL) {
      if (i == pnum[4 * nshft[npr] + 3]) pvect[3] = atof (pch);
      if (i == pnum[4 * nshft[npr] + 7]) pvect[7] = atof (pch);
      if (i == pcolor[npr]) strcpy (bufff[0], pch);
      if (i == pscale[npr]) qcdscale = atof (pch);
      for (j = 2; j < ntot[npr]; ++j) {
        for (k = 1; k < 4; ++k)
          if (i == pnum[4 * nshft[npr] + 4 * j + k]) pvect[4 * j + k] = atof (pch);
      }
      pch = strtok (NULL, ":");
      ++i;
    }
    for (i = 0; i < ntot[npr]; i++)
      pvect[4 * i] = ENERGY (mass[nshft[npr] + i], pvect + 4 * i + 1);

  /* color chains */
    mpar = sscanf (bufff[0], "(%d %d)%[^a]", &n1, &n2, bufff[1]);
    itag = 500;
    for (i = 0; i < ntot[npr]; ++i) col[i][0] = col[i][1] = 0;
    while (3 == mpar || 2 == mpar) {
      --n1;
      --n2;
      if (IPSGN == -1) {
        if (n1 == 0)      n1 = 1;
        else if (n1 == 1) n1 = 0;
        if (n2 == 0)      n2 = 1;
        else if (n2 == 1) n2 = 0;
      }
      if (n1 < 2) col[n1][0] = itag;
      else        col[n1][1] = itag;
      if (n2 < 2) col[n2][1] = itag;
      else        col[n2][0] = itag;
      ++itag;
      mpar = sscanf (bufff[itag - 500], "(%d %d)%[^a]", &n1, &n2, bufff[itag - 499]);
    }

  /* fill the LesHouches structure */
    ev.NpartUP    = ntot[npr];
    ev.IDprocUP   = npr + 1;
    ev.XweightUP  = 1.;
    ev.QscaleUP   = qcdscale;
    ev.QEDalphaUP = alpha_em (qcdscale);
    ev.QCDalphaUP = alpha_2 (qcdscale);

    for (i = 0; i < ntot[npr]; ++i) {
      ev.IDpartUP[i]       = id[nshft[npr] + j];
      ev.statusUP[i]       = 1;
      ev.motherUP[0][i]    = 1;
      ev.motherUP[1][i]    = 2;
      if (2 > i) {
        ev.statusUP[i]     = -1;
        ev.motherUP[0][i]  = 0;
        ev.motherUP[1][i]  = 0;
      }
      ev.colorflowUP[0][i] = col[i][0];
      ev.colorflowUP[1][i] = col[i][1];
      ev.momentumUP[0][i]  = pvect[4 * i + 1];
      ev.momentumUP[1][i]  = pvect[4 * i + 2];
      ev.momentumUP[2][i]  = pvect[4 * i + 3];
      ev.momentumUP[3][i]  = pvect[4 * i];
      ev.momentumUP[4][i]  = mass[nshft[npr] + i];
      ev.timelifeUP[i]     = 0.;
      ev.spinUP[i]         = 9.;
    }
    fill_event (&ev);

    if (Nevents > 0 && 0 == Nevents % 10000) fprintf (stdout, "rtupler (info): %i events read\n", Nevents);
    ++Nevents;
    if (nevnt > 0 && nevnt <= Nevents) {
      break;
    }
  }
  write_rtuple ();

  free (nin);
  free (ntot);
  free (nshft);
  free (id);
  free (mass);
  free (pcs);
  free (pcserr);
  free (pcolor);
  free (pscale);

  fclose (s);

  fprintf (stdout, "rtupler (info): %i events have been written to the rtuple %s\n", Nevents, out_name);
  return 0;
}



#else
int
rtuple_cpyth2 (char ini_name[], char out_name[], int nevnt)
{
  fprintf (stdout, "rtupler (error): ROOT is not linked. Define ROOTSYS and re-compile CompHEP. Exit.\n");
  return 0;
}
#endif
