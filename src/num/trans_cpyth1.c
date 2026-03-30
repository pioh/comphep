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
#include "trans_cpyth1.h"

#define MAXFUN 4
#define FUNLEN 40

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
  pch = strtok (buff, " ,");
  while (pch != NULL && strcmp(pch, "\n")) {
    if (strcmp ("->", pch)) {
      strcpy (pname[ntot_], pch);
      ++ntot_;
    }
    pch = strtok (NULL, " ,");
  }

  if (0 == ntot_) {
    fprintf (stderr, "traslator (error): strange (sub)process name, can't extract parton names and set nin_/nout_\n");
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
static midstr sf_txt[2];

static int
readPDF (FILE * f)
{
  int i;
  double SqrtS;
  double Rapidity;
  midstr buff;
  double prt_mass[2] = {0.,0.};
  double be, mass;

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

  if (1 == nin_) {
    return 0;
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
        fprintf (stderr, "translator (warning): the event file has internal CompHEP PDF, but traslator could not open\n");
        fprintf (stderr, "                      the PDF file. So, CompHEP internal alpha_S representation will be used\n");
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
        fprintf (stderr, "translator (warning): the event file has LHA PDF, but traslator could not find LHAPDF data\n");
        fprintf (stderr, "                      files. So, CompHEP internal alpha_S representation will be used\n");
      }
    }
#else
    if (NULL != lha) {
      fprintf (stderr, "translator (warning): the event file has LHA PDF, but traslator is compiled without LHAPDF\n");
      fprintf (stderr, "                      support. So, CompHEP internal alpha_S representation will be used\n");
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
translate_cpyth1 (const char source[], const char target[], int mixing_done) 
{
  int i, j;
  int Nproc;
  int nntot;
  int n1, n2, mpar;
  int itag;
  int stop = 0;
  int Nevents = 0;

  int * nin;
  int * nout;
  int * ntot;
  int * nshft;
  int * id;

  double * pcs;
  double * mass;
  double qcdscale;

  midstr buff[MAXINOUT + 3];
  FILE * s = fopen (source, "r");
  FILE * t = fopen (target, "w");

  rw_paragraph rd_array0[1] = {
    {"PEVLIB_v.1.0", skipFinalLine}
  };
  rw_paragraph rd_array1[7] = {
    {"CompHEP", skipLine},
    {"PROCESS", getNames},
    {"Initial_state", readPDF},
    {"MASSES", getMasses},
    {"Cross_section(Width)", readCS},
    {"Number_of_events", skipLine},
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
    fprintf (stderr, "translator (error): strange format, can't find #Number_of_subprocesses string. Exit\n");
    return -1.0;
  }
  if (Nproc < 1) {
    fprintf (stderr, "translator (error): strange number of subprocesses. Nproc = %i\n", Nproc);
    return -1.0;
  }

  /* Return to file begin and read first tag: PEVLIB version */
  rewind (s);
  if (mixing_done) {
    readParagraphs (s, 1, rd_array0);
  }

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

  /* some checks for decays */
    if (0 == i && 1 == nin_) {
      fprintf (stderr, "translator (warning): decay file!\n");
      ebeam[0] = pmass[0];
      ebeam[1] = 0;
      pbeam[0] = kfpart(pname[0]);
      pbeam[1] = 1;
    }
    if (0 < i && 1 == nin_) {
      if ((fabs (ebeam[0] - pmass[0]) > 1e-10) || pbeam[0] != kfpart(pname[0])) {
        fprintf (stderr, "translator (error): this is a decay file, but decay particles are differenc in some subprocesses... Exit\n");
        exit (99);
      }
    }
  }

  nshft[0] = 0;
  for (i = 1; i < Nproc; ++i) {
    nshft[i] = nshft[i - 1] + ntot[i - 1];
  }

  nntot = nshft[Nproc - 1] + ntot[Nproc - 1];
  mass = malloc (nntot * sizeof (double));
  id   = malloc (nntot * sizeof (double));
  rewind (s);
  if (mixing_done) {
    readParagraphs (s, 1, rd_array0);
  }
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

  fprintf (t, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (t, "<!-- File generated with CompHEP %s -->\n", version ());
  fprintf (t, "<!-- \n"
                    "     The event file is compatible the Les Houches event file format (hep-ph/0609017)\n"
                    "     translated from the CompHEP cpyth1 format, so it has no HepML header\n"
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
    fprintf (t, "%17.10E %17.10E %17.10E %i\n", pcs[i], 0.0, 1.0, i + 1);
  }
  fprintf (t, "</init>\n");

  while (fgets (buff[0], 2048, s)) {
    int npr;
    longstr init[2];
    midstr color_string[MAXINOUT + 3];
    int col[MAXINOUT][2];

    sscanf (buff[0], "%i%[^\n]", &npr, buff[1]);
    --npr;
    if (npr < 0) {
      fprintf (stderr, "translator (error): negative process number. Event skipped\n");
      continue;
    }
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
    for (j = 1; j < MAXINOUT + 3; ++j) strcpy (color_string[j], "");
    itag = 500;
    for (j = 0; j < MAXINOUT; ++j) col[j][0] = col[j][1] = 0;
    mpar = sscanf (color_string[0], "   (%d %d)%[^a]", &n1, &n2, color_string[1]);
    while (2 == mpar || 3 == mpar) {
      --n1;
      --n2;
      if (n1 < nin[npr]) col[n1][0] = itag;
      else          col[n1][1] = itag;
      if (n2 < nin[npr]) col[n2][1] = itag;
      else          col[n2][0] = itag;
      ++itag;
      mpar = sscanf (color_string[itag - 500], "(%d %d)%[^a]", &n1, &n2, color_string[itag - 499]);
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
  fprintf (stdout, "translator (info): %i events (cpyth1) have been translated to the LHE format, stored in %s\n", Nevents, target);

  free (mass);
  free (id);
  free (pcs);
  free (nin);
  free (nshft);
  free (nout);
  free (ntot);

  fclose (s);
  fclose (t);

  return 0;
}

static int 
getNames_calchep (FILE * f)
{
  int i = 0;
  char * pch;
  midstr buff;
  midstr temp;

  fgets (buff, STRSIZ, f);
  strcpy (temp, buff);
  pch = strtok (buff, " ");
  while (pch != NULL && strcmp (pch, "\n")) {
    if (strcmp ("->", pch)) {
      int id = 0;
      int err = sscanf (pch, "%d(%s)", &id, pname[i]);
      if (2 != err) {
        fprintf (stderr, "traslator (error): str: %s, id = %i, name: %s\n", pch, id, pname[i]);
      } else {
        pname[i][strlen(pname[i]) - 1] = '\0';
        ++i;
      }
    }
    pch = strtok (NULL, " ,");
  }

  if (0 == i) {
    fprintf (stderr, "traslator (error): strange process name, can't extract parton names\n");
    return -1;
  }

  return 0;
}


static int
readType (FILE * f)
{
  int err;
  int tnout;
  int tnin;
  midstr buff;
  fgets (buff, STRSIZ, f);
  nin_ = nout_ = 0;

  err = sscanf (buff, "  %d -> %d", &tnin, &tnout);
  if (2 == err) {
    nin_ = tnin;
    nout_ = tnout;
  }

  return 0;
}


static int
readPDF_calchep (FILE * f)
{
  int i;
  int err = 0;
  double e1;
  double e2;
  midstr buff;
  char * bname;

  ebeam[0] = ebeam[1] = 0.0;

  fgets (buff, STRSIZ, f);
  err = sscanf (buff, "  P1_3=%lf  P2_3=%lf", &e1, &e2);
  if (2 == err) {
    ebeam[0] = e1;
    ebeam[1] = e2;
  }

  set_alphaMode (0);
  for (i = 0; i < 2; i++) {
    set_sf_mass (i, 0.);
    set_sf_be (i, 0.);
    set_sf_num (i, 0);
  }

  if (1 == nin_) {
    return 0;
  }

  for (i = 0; i < 2; ++i) {
    midstr pdf;
    int idbeam;
    int err = 0;
    char * buf;
    PDFLIBset[i] = PDFLIBgroup[i] = pbeam[i] = 0;
    fgets (buff, 1024, f);
    buf = strstr (buff, "PDT:");
    err = sscanf (buf, "PDT:%s\"%d", pdf, &idbeam);

    if (2 != err) {
      fprintf (stderr, "translator (error): strange PDF string: %s\n", buff);
    }

    if (0 != strlen(pdf)) {
      int len;
      shortstr pdfver;
      shortstr pdfname;
      strcpy (sf_txt[i], pdf);
      strncpy (pdfname, pdf, 4); pdfname[4] = 0;
      strcpy (pdfver, pdf + 4);
      bname = strstr (pdf, "(");
      len = strlen (bname);
      bname[len - 2] = '\0';
      len = strlen (pdfver) - len;
      pdfver[len] = '\0';
      CERNpdf_number (pdfname, pdfver, &(PDFLIBset[i]), &(PDFLIBgroup[i]));
      pbeam[i] = kfbeam (bname + 1);
    }
  }

  return 0;
}


int 
translate_calchep (const char source[], const char target[], int mixing_done) 
{
  int i, j;
  int n1, n2, mpar;
  int itag;
  int Nproc = 1;
  int Nevents = 0;
  int nin;
  int nout;
  int ntot;
  int id[MAXINOUT];
  double pcs;
  double mass[MAXINOUT];
  double qcdscale;
  double be, pmasss;
  midstr buff[MAXINOUT + 3];

  FILE * s = fopen (source, "r");
  FILE * t = fopen (target, "w");

  rw_paragraph rd_array[8] = {
    {"CalcHEP", skipLine},
    {"Type", readType},
    {"Initial_state", readPDF_calchep},
    {"PROCESS", getNames_calchep},
    {"MASSES", getMasses},
    {"Cross_section(Width)", readCS},
    {"Number_of_events", skipLine},
    {"Events", skipFinalLine}
  };

  /* Read subprocess info: particle names and masses */
  readParagraphs (s, 8, rd_array);
  pcs  = totCS;
  nin  = nin_;
  nout = nout_;
  ntot = nin_ + nout_;

  /* if we have decay file */
  if (1 == nin) {
    fprintf (stderr, "translator (info): decay file!\n");
    ebeam[0] = pmass[0];
    ebeam[1] = 0;
    pbeam[0] = kfpart(pname[0]);
    pbeam[1] = 1;
  } else {
    for (i = 0; i < 2; ++i) {
      int done = 0;
      if (strstr (sf_txt[i], "PDF:")) {
        if (strFunPdf.myParticle (pname[i]) && strFunPdf.readName (i + 1, sf_txt[i]) && strFunPdf.init (i + 1, &be, &pmasss, pname[i])) {
          set_sf_num(i, 1);
          set_sf_mass (i, pmasss);
          set_sf_be (i, be);
          set_alphaMode (1);
          done = 1;
        }
      }
#ifdef LHAPDF
      if (!done && strstr (sf_txt[i], "LHA:")) {
        if (strFunLha.myParticle (pname[i]) && strFunLha.readName (i + 1, sf_txt[i]) && strFunLha.init (i + 1, &be, &pmasss, pname[i])) {
          set_sf_num(i, 1);
          set_sf_mass (i, pmasss);
          set_sf_be (i, be);
          set_alphaMode (1);
          done = 1;
        }
      }
#endif
      if (!done) {
        fprintf (stderr, "translator (warning): could not open the PDF file. CompHEP internal alpha_S will be used\n");
      }
    }
  }

  for (j = 0; j < ntot; ++j) {
    double mas;
    char name[16];
    pinfo (j, NULL, &mas);
    mass[j] = mas;
    pinfo (j, name, NULL);
    id[j] = kfpart (name);
  }

  fprintf (t, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (t, "<!-- File generated with CompHEP %s -->\n", version ());
  fprintf (t, "<!-- \n"
                    "     The event file is compatible the Les Houches event file format (hep-ph/0609017)\n"
                    "     translated from the CalcHEP format, so it does not have HepML header\n"
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
  fprintf (t, "%17.10E %17.10E %17.10E %i\n", pcs, 0.0, 1.0, i + 1);
  fprintf (t, "</init>\n");

  while (fgets (buff[0], 2048, s)) {
    int npr;
    longstr init[2];
    midstr color_string[MAXINOUT + 3];
    int col[MAXINOUT][2];

    sscanf (buff[0], "%i%[^\n]", &npr, buff[1]);
    --npr;
    for (i = 0; i < 8; i++) pvect[i] = 0.;
    if (2 == nin) {
      sscanf (buff[1], "    %lf %lf%[^a]", pvect + 3, pvect + 7, buff[2]);
    }
    for (i = nin; i < ntot; i++) {
      double * p = pvect + 4 * i;
      sscanf (buff[i], " %lf %lf %lf%[^a]", p + 1, p + 2, p + 3, buff[i + 1]);
    }
    for (i = 0; i < ntot; i++) {
      pvect[4 * i] = ENERGY (mass[i], pvect + 4 * i + 1);
    }

  /* color chains */
    sscanf (buff[i], "| %lf%[^a]", &qcdscale, color_string[0]);
    for (j = 1; j < MAXINOUT + 3; ++j) strcpy (color_string[j], "");
    mpar = sscanf (color_string[0], "   (%d %d)%[^a]", &n1, &n2, color_string[1]);
    itag = 500;
    for (j = 0; j < MAXINOUT; ++j) col[j][0] = col[j][1] = 0;
    while (2 == mpar || 3 == mpar) {
      --n1;
      --n2;
      if (n1 < nin) col[n1][0] = itag;
      else          col[n1][1] = itag;
      if (n2 < nin) col[n2][1] = itag;
      else          col[n2][0] = itag;
      ++itag;
      mpar = sscanf (color_string[itag - 500], "(%d %d)%[^a]", &n1, &n2, color_string[itag - 499]);
    }

  /* print out event */
    fprintf (t, "<event>\n");
    fprintf (t, "%i %i 1.0 %17.10E %17.10E %17.10E\n", ntot, npr + 1, qcdscale, alpha_em (qcdscale), alpha_2 (qcdscale));
    for (j = 0; j < nin; ++j) {
      int k = 4 * j;
      sprintf (init[j], "%i -1 0 0 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
      id[j], col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], mass[j]);
    }
    fputs (init[0], t);
    fputs (init[1], t);
    for (j = nin; j < ntot; ++j) {
      int k = 4 * j;
      if (2 == nin) {
        fprintf (t, "%i 1 1 2 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
        id[j], col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], mass[j]);
      } else {
        fprintf (t, "%i 1 1 0 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
        id[j], col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], mass[j]);
      }
    }
    fprintf (t, "</event>\n");
    ++Nevents;
    if (0 == Nevents % 20000) fprintf (stdout, "translator (info): %i events processed\n", Nevents);
  }
  fprintf (stdout, "translator (info): %i events (cpyth1) have been translated to the LHE format, stored in %s\n", Nevents, target);

  fclose (s);
  fclose (t);

  return 0;
}
