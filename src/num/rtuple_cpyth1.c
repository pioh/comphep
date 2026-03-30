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
#include "service2/include/paragraphs.h"
#include "service2/include/lbl.h"
#include "service2/include/4_vector.h"
#include "service2/include/kfcodes.h"
#include "service2/include/drandXX.h"

#include "alphas2.h"
#include "sf_pdf.h"
#include "pdf.h"
#include "sf_lhapdf.h"
#include "strfun_par.h"
#include "LesHouches.h"
#include "rtuple_routines.h"

#include "rtuple_cpyth1.h"

static int nin_ = 0;
static int nout_ = 0;
static double totCS = 0.0;
static char pname[MAXNP][20];
static double pmass[MAXNP];

static struct {
  int (*myParticle) (char *);
  void (*fullName) (int, char *, char *);
  void (*realSTRFUN_info) (int, Str_fun_Info *);
  int (*readName) (int, char *);
  int (*beam_menu) (int);
  int (*menu) (int, char *);
  int (*init) (int, double *, double *, char *);
  double (*val) (int, double, double);
} strFunPdf =
  {
    p_pdf, n_pdf, info_pdf, r_pdf, beam_pdf, pdf_pdf, be_pdf, c_pdf
  };
#ifdef LHAPDF
static struct {
  int (*myParticle) (char *);
  void (*fullName) (int, char *, char *);
  void (*realSTRFUN_info) (int, Str_fun_Info *);
  int (*readName) (int, char *);
  int (*beam_menu) (int);
  int (*menu) (int, char *);
  int (*init) (int, double *, double *, char *);
  double (*val) (int, double, double);
} strFunLha =
  {
    p_lhapdf, n_lhapdf, info_lhapdf, r_lhapdf, beam_lhapdf, pdf_lhapdf, be_lhapdf, c_lhapdf
  };
#endif

int
pinfo (int num, char *name, double *mass)
{
  if (name)
    strcpy (name, pname[num]);
  if (mass)
    *mass = pmass[num];

  return 0;
}

static int
getMasses (FILE * flow)
{
  midstr buff;
  char *pos;
  int i = 0;

  fgets (buff, STRSIZ, flow);
  pos = strtok (buff, " ");
  while (pos)
    {
      sscanf (pos, "%lf", pmass + i);
      i++;
      pos = strtok (NULL, " ");
    }
  return 0;
}

static int
skipLine (FILE * flow)
{
  fscanf (flow, "%*[^\n]\n");
  return 0;
}

static int skipFinalLine (FILE * flow) /* SHOUD RETURN -1 IN ORDER TO END PARAGRAPH READING!!! */
{
  fscanf (flow, "%*[^\n]\n");
  return -1;
}

static int 
getNames (FILE * f)
{
  int ntot_ = 0;
  char * pch;
  char * buf;
  midstr buff;
  midstr temp;

  fgets (buff, STRSIZ, f);
  strcpy (temp, buff);
  pch = strtok (buff, " ,->");
  while (pch != NULL && strcmp(pch, "\n")) {
    strcpy (pname[ntot_], pch);
    pch = strtok (NULL, " ,->");
    ++ntot_;
  }

  if (0 == ntot_) {
    fprintf (stderr, "mk_tab (error): strange (sub)process name, can't extract parton names and set nin_/nout_\n");
    return -1;
  }

  nout_ = 0;
  buf = strstr(temp, "->") + 3;
  pch = strtok (buf, " ,");
  while (pch != NULL && strcmp(pch, "\n")) {
    ++nout_;
    pch = strtok (NULL, " ,");
  }
  nin_ = ntot_ - nout_;
  return 0;
}

static int
readTotCS (FILE * f)
{
  midstr buff;
  fgets (buff, MISSTRLEN, f);
  sscanf (buff, "= %lf", &totCS);
  return 0;
}

static int
readCS (FILE * f)
{
  midstr buff;
  fgets (buff, MISSTRLEN, f);
  sscanf (buff, " %lf", &totCS);
  return 0;
}

static int pbeam[2];
static double ebeam[2];
static int PDFLIBgroup[2];
static int PDFLIBset[2];

static int
readPDF (FILE * f)
{
  int i;
  double SqrtS;
  double Rapidity;
  midstr buff;
  double prt_mass[2] = {0.,0.};
  double be, mass;
  midstr sf_txt[2];

  fgets (buff, MISSTRLEN, f);
  sscanf (buff, "  SQRT(S) %lf", &SqrtS);
  fgets (buff, MISSTRLEN, f);
  sscanf (buff, "  Rapidity(c.m.s) %lf", &Rapidity);

  set_alphaMode (0);
  for (i = 0; i < 2; i++) {
    set_sf_mass (i, 0.);
    set_sf_be (i, 0.);
    set_sf_num (i, 0);
  }

  for (i = 0; i < 2; ++i) {
    double e, p;
    char * pdf;
    char * lha;
    char * bname;
    int len;
    fgets (buff, 1024, f);
    pdf = strstr (buff, "PDF:");
    lha = strstr (buff, "LHA:");
    PDFLIBset[i] = PDFLIBgroup[i] = 0;
    if (NULL != pdf) {
      shortstr pdfver;
      shortstr pdfname;
      strcpy (sf_txt[i], pdf);
      strncpy (pdfname, pdf + 4, 4); pdfname[4] = 0;
      strcpy (pdfver, pdf + 8);
      bname = strstr (pdf, "(");
      len = strlen (bname);
      bname[len - 2] = '\0';
      len = strlen (pdfver) - len;
      pdfver[len] = '\0';
      CERNpdf_number (pdfname, pdfver, &(PDFLIBset[i]), &(PDFLIBgroup[i]));

      if (strFunPdf.myParticle (pname[i]) && strFunPdf.readName (i + 1, sf_txt[i]) && strFunPdf.init (i + 1, &be, &mass, pname[i])) {
        set_sf_num(i, 1);
        set_sf_mass (i, mass);
        set_sf_be (i, be);
        set_alphaMode (1);
      } else {
        fprintf (stderr, "rtuples (warning): the event file has internal CompHEP PDF, but could not open\n");
        fprintf (stderr, "                   the PDF file. So, CompHEP internal alpha_S representation will be used\n");
      }
    }
#ifdef LHAPDF
    if (NULL != lha) {
      strcpy (sf_txt[i], lha);
      bname = strstr (lha, "(");
      len = strlen (bname);
      bname[len - 2] = '\0';
      if (strFunLha.myParticle (pname[i]) && strFunLha.readName (i + 1, sf_txt[i]) && strFunLha.init (i + 1, &be, &mass, pname[i])) {
        set_sf_num(i, 1);
        set_sf_mass (i, mass);
        set_sf_be (i, be);
        set_alphaMode (1);
      } else {
        fprintf (stderr, "rtuples (warning): the event file has LHA PDF, but could not find LHAPDF data\n");
        fprintf (stderr, "                   files. So, CompHEP internal alpha_S representation will be used\n");
      }
    }
#else
    if (NULL != lha) {
      fprintf (stderr, "rtuples (warning): the event file has LHA PDF, but compiled without LHAPDF\n");
      fprintf (stderr, "                   support. So, CompHEP internal alpha_S representation will be used\n");
    }
#endif
    e = (SqrtS * SqrtS + prt_mass[i] * prt_mass[i] - prt_mass[1 - i] * prt_mass[1 - i]) / (2 * SqrtS);
    p = sqrt (e * e - prt_mass[i] * prt_mass[i]);
    p = p * cosh (Rapidity) + e * sinh (Rapidity) * (1 - 2 * i);
    if (p * p < (10.E-10) * SqrtS) p = 0.0;
    ebeam[i] = p;
    pbeam[i] = kfbeam (bname + 1);
  }

  return 0;
}

int
rtuple_cpyth1 (char ini_name[], char out_name[], int nevnt)
{
  int i, j;
  int Nproc;
  int nntot;
  int n1, n2, mpar;
  int itag;
  int stop = 0;
  int Nevents = 0;
  int err;

  int * nin;
  int * nout;
  int * ntot;
  int * nshft;
  int * id;

  double * pcs;
  double * mass;
  double qcdscale;

  midstr buff[MAXINOUT + 3];
  FILE * s = fopen (ini_name, "r");

  rw_paragraph rd_array0[1] = {
    {"PEVLIB_v.1.0", skipFinalLine}
  };
  rw_paragraph rd_array1[7] = {
    {"CompHEP", skipLine},
    {"PROCESS", getNames},
    {"Initial_state", readPDF},
    {"MASSES", getMasses},
    {"Cross_section(Width)", skipLine},
    {"Number_of_events", readCS},
    {"----------------------------------------------------------", skipFinalLine}
  };
  rw_paragraph rd_array2[4] = {
    {"Number_of_subprocesses", skipLine},
    {"Total_cross_section_(pb)", readTotCS},
    {"Events_mixed_and_randomized", skipLine},
    {"Nproc", skipFinalLine}
  };

  /* Get the number of subprocesses in the file */
  while ((NULL == fgets (buff[0], 1024, s)) || !stop) {
    stop = sscanf (buff[0], "#Number_of_subprocesses = %d", &Nproc);
  }
  if (feof (s)) {
    fprintf (stderr, "rtuples (error): strange format, can't find #Number_of_subprocesses string. Exit\n");
    return -1.0;
  }
  if (Nproc < 1) {
    fprintf (stderr, "rtuples (error): strange number of subprocesses. Nproc = %i\n", Nproc);
    return -1.0;
  }

  /* Return to file begin and read first tag: PEVLIB version */
  rewind (s);
  readParagraphs (s, 1, rd_array0);

  /* Read subprocess info: particle names and masses */
  nin   = malloc (Nproc * sizeof (int));
  nout  = malloc (Nproc * sizeof (int));
  ntot  = malloc (Nproc * sizeof (int));
  nshft = malloc (Nproc * sizeof (int));
  pcs   = malloc (Nproc * sizeof (double));
  for (i = 0; i < Nproc; ++i) {
    readParagraphs (s, 7, rd_array1);
    pcs[i]  = totCS;
    nin[i]  = nin_;
    nout[i] = nout_;
    ntot[i] = nin_ + nout_;
  }

  nshft[0] = 0;
  for (i = 1; i < Nproc; ++i) {
    nshft[i] = nshft[i - 1] + ntot[i - 1];
  }

  nntot = nshft[Nproc - 1] + ntot[Nproc - 1];
  mass = malloc (nntot * sizeof (double));
  id   = malloc (nntot * sizeof (double));
  rewind (s);
  readParagraphs (s, 1, rd_array0);
  for (i = 0; i < Nproc; ++i) {
    char name[20];
    double mas;
    readParagraphs (s, 7, rd_array1);
    for (j = 0; j < ntot[i]; ++j) {
      pinfo (j, NULL, &mas);
      mass[nshft[i] + j] = mas;
      pinfo (j, name, NULL);
      id[nshft[i] + j] = kfpart (name);
    }
  }
  readParagraphs (s, 4, rd_array2);

  err = book_rtuple (out_name);
  if (err) {
   fprintf (stderr, "rtupler (error): can't create the rtuple %s. exit...\n", out_name);
   return -1;
  }

  while (fgets (buff[0], 2048, s)) {
    int npr;
    int IPSGN = 1;
    midstr color_string[MAXINOUT + 3];
    int col[MAXINOUT][2];
    eventUP ev;

    sscanf (buff[0], "%i%[^\n]", &npr, buff[1]);
    --npr;
    for (i = 0; i < 8; i++) pvect[i] = 0.;
    sscanf (buff[1], "    %lf %lf%[^a]", pvect + 3, pvect + 7, buff[2]);
    for (i = nin_; i < ntot[npr]; i++) {
      double * p = pvect + 4 * i;
      sscanf (buff[i], " %lf %lf %lf%[^a]", p + 1, p + 2, p + 3, buff[i + 1]);
    }
    for (i = 0; i < ntot[npr]; i++) {
      pvect[4 * i] = ENERGY (mass[nshft[npr] + i], pvect + 4 * i + 1);
    }

  /* color chains */
    sscanf (buff[i], " %lf%[^a]", &qcdscale, color_string[0]);
    mpar = sscanf (color_string[0], "   (%d %d)%[^a]", &n1, &n2, color_string[1]);
    itag = 500;
    for (i = 0; i < ntot[npr]; ++i) col[i][0] = col[i][1] = 0;
    while (2 == mpar || 3 == mpar) {
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
      mpar = sscanf (color_string[itag - 500], "(%d %d)%[^a]", &n1, &n2, color_string[itag - 499]);
    }

  /* fill the LesHouches structure */
    ev.NpartUP    = ntot[npr];
    ev.IDprocUP   = npr + 1;
    ev.XweightUP  = 1.;
    ev.QscaleUP   = qcdscale;
    ev.QEDalphaUP = alpha_em (qcdscale);
    ev.QCDalphaUP = alpha_2 (qcdscale);

    for (i = 0; i < ntot[npr]; ++i) {
      ev.IDpartUP[i]       = id[nshft[npr] + i];
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

  free (mass);
  free (id);
  free (pcs);
  free (nin);
  free (nshft);
  free (nout);
  free (ntot);

  fclose (s);

  fprintf (stdout, "rtupler (info): %i events have been written to the rtuple %s\n", Nevents, out_name);
  return 0;
}



#else
int
rtuple_cpyth1 (char ini_name[], char out_name[], int nevnt)
{
  fprintf (stdout, "rtupler (error): ROOT is not linked. Define ROOTSYS and re-compile CompHEP. Exit.\n");
  return 0;
}
#endif
