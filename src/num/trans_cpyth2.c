/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/lbl.h"
#include "service2/include/4_vector.h"
#include "service2/include/kfcodes.h"
#include "service2/include/drandXX.h"

#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_routines.h"
#include "alphas2.h"
#include "strfun.h"
#include "trans_cpyth2.h"

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
  char * temp1 = malloc (strlen (buff) * sizeof (char));

  nin_ = 0;
  strcpy (temp, primetrim (buff));
  strcpy (temp1, temp);
  pch = strtok (temp, " ,->");
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
  buf = strstr(temp1, "->") + 3;
  pch = strtok (buf, " ,");
  while (pch != NULL) {
    ++nout_;
    pch = strtok (NULL, " ,");
  }
  nin_ -= nout_;

  return 0;
}

int 
translate_cpyth2 (const char source[], const char target[]) 
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

  int * pcolor = NULL;
  int * pscale = NULL;
  int * nin = NULL;
  int * nshft = NULL;
  int * ntot = NULL;
  int * the_type = NULL;
  int * id = NULL;
  int * pnum = NULL;

  midstr buff;
  midstr buf;
  midstr name_proc;
  midstr type_proc;

  double cs, cserr;
  double qcdscale;
  char bufff[MAXINOUT + 3][STRSIZ];
  double * pcs = NULL;
  double * pcserr = NULL;
  double * mass = NULL;

  string_comnd com;
  tags * head = init_cap (1);
  FILE * s = fopen (source, "r");
  FILE * t = fopen (target, "w");

  cup_reader (s, head);
  strcpy (com.name, "Nproc");
  get_tag_with1com (0, head, "total", &com);
  Nproc = atoi (com.value);

  nin      = malloc (Nproc * sizeof (int));
  ntot     = malloc (Nproc * sizeof (int));
  the_type = malloc (Nproc * sizeof (int));
  pcs      = malloc (Nproc * sizeof (double));
  pcserr   = malloc (Nproc * sizeof (double));

  /* check file format */
  if (check_cpyth2 (head)) {
    return -1;
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
    if (n < 0) {
      fprintf (stderr, "translator (error): can not fine the neseccary tag \"process\" in the %i-th process. Exit", i+1);
      return -2;
    }
    pcs[i]      = get_fval (0, "CrosSec", head->tag[n]);
    pcserr[i]   = get_fval (0, "CrosSecErr", head->tag[n]);
    strcpy (name_proc, get_cval (0, "name", head->tag[n]));
    strcpy (type_proc, get_cval (0, "type", head->tag[n]));
    the_type[i] = 0;
    if (0 == strcmp (type_proc, "\'decay\'")) the_type[i] = 1;
    if (0 == strcmp (type_proc, "\'scattering\'")) the_type[i] = 2;
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

  for (i = 0; i < Nproc; i++) {
    if (the_type[i] != the_type[0]) {
      fprintf (stderr, "translator (error): diffenrent subprocesses are mixed. Exit");
      return -3;
    }
    if (1 != the_type[i] && 2 != the_type[i]) {
      fprintf (stderr, "translator (error): strange %i-th suprocess: type = %i\n", i + 1, the_type[i]);
      return -4;
    }
  }

  if (2 == the_type[0]) {
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
  } else {
    ebeam[0] = mass[0];
    ebeam[1] = 0.0;
    pbeam[0] = id[0];
    pbeam[1] = 0;
    PDFLIBset[0]   = 0;
    PDFLIBset[1]   = 0;
    PDFLIBgroup[0] = 0;
    PDFLIBgroup[1] = 0;
  }

  fprintf (t, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (t, "<!-- File generated with CompHEP %s -->\n", version ());
  fprintf (t, "<!-- \n"
                    "     The event file is compatible the Les Houches event file format (hep-ph/0609017)\n"
                    "     translated from the CompHEP cpyth2 format, so it has no HepML header\n"
                    "-->\n");
  fprintf (t, "<header>\n");
  fprintf (t, "</header>\n");
  fprintf (t, "<init>\n");
  fprintf (t, "%i %i %17.10E %17.10E %i %i %i %i 3 %i\n",
              pbeam[0], 
              pbeam[1], 
              ebeam[0], 
              ebeam[1], 
              PDFLIBgroup[0], 
              PDFLIBgroup[1], 
              PDFLIBset[0], 
              PDFLIBset[1],
              Nproc);
  for (i = 0; i < Nproc; ++i) {
    fprintf (t, "%17.10E %17.10E %17.10E %i\n", pcs[i], pcserr[i], 1.0, i + 1);
  }
  fprintf (t, "</init>\n");

  fgets (buf, 2048, s);
  while (fgets (buf, 2048, s)) {
    int npr;
    longstr init[2];
    int col[MAXINOUT][2];

    sscanf (buf, "%i:%[^\n]", &npr, buff);
    --npr;
    for (i = 0; i < 8; i++) pvect[i] = 0.;
    i = 0;
    char * pch = strtok (buff, ":");
    while (pch != NULL) {
      if (i == pcolor[npr]) strcpy (bufff[0], pch);
      if (i == pscale[npr]) qcdscale = atof (pch);
      if (2 == nin[npr]) {
        if (i == pnum[4 * nshft[npr] + 3]) pvect[3] = atof (pch);
        if (i == pnum[4 * nshft[npr] + 7]) pvect[7] = atof (pch);
        for (j = 2; j < ntot[npr]; ++j) {
          for (k = 1; k < 4; ++k)
            if (i == pnum[4 * nshft[npr] + 4 * j + k]) pvect[4 * j + k] = atof (pch);
        }
      } else {
        for (j = 1; j < ntot[npr]; ++j) {
          for (k = 1; k < 4; ++k)
            if (i == pnum[4 * nshft[npr] + 4 * j + k]) pvect[4 * j + k] = atof (pch);
        }
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
      if (n1 < nin[npr]) col[n1][0] = itag;
      else               col[n1][1] = itag;
      if (n2 < nin[npr]) col[n2][1] = itag;
      else               col[n2][0] = itag;
      ++itag;
      mpar = sscanf (bufff[itag - 500], "(%d %d)%[^a]", &n1, &n2, bufff[itag - 499]);
    }

  /* print out event */
    fprintf (t, "<event>\n");
    fprintf (t, "%i %i 1.0 %17.10E %17.10E %17.10E\n", ntot[npr], npr + 1, qcdscale, alpha_em (qcdscale), alpha_2 (qcdscale));
    for (j = 0; j < nin[npr]; ++j) {
      int k = 4 * j;
      sprintf (init[j], "%i -1 0 0 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
      id[nshft[npr] + j], col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], mass[nshft[npr] + j]);
    }
    fputs (init[0], t);
    fputs (init[1], t);
    for (j = nin[npr]; j < ntot[npr]; ++j) {
      int k = 4 * j;
      if (2 == nin[npr]) {
        fprintf (t, "%i 1 1 2 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
        id[nshft[npr] + j], col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], mass[nshft[npr] + j]);
      } else {
        fprintf (t, "%i 1 1 0 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
        id[nshft[npr] + j], col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], mass[nshft[npr] + j]);
      }
    }
    fprintf (t, "</event>\n");
    ++Nevents;
    if (0 == Nevents % 20000) fprintf (stdout, "translator (info): %i events processed\n", Nevents);
  }

  free (nin);
  free (ntot);
  free (the_type);
  free (nshft);
//  free (pnum);
  free (id);
  free (mass);
  free (pcs);
  free (pcserr);
  free (pcolor);
  free (pscale);
  fclose (s);
  fclose (t);
  fprintf (stdout, "translator: %i events (cpyth2) have been translated to the LHE format, stored in %s\n", Nevents, target);

  return 0;
}
