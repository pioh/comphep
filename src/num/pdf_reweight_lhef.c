/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Author: V.A.Ilyin
* ------------------------------------------------------
*/
#include <unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/files.h"

#include "LesHouches.h"
#include "lhef_routines.h"
#include "lhapdf.h"
#include "event_reader.h"
#include "strfun_par.h"
#include "alphas2.h"

#include "service2/include/syst.h"
#include "pdf_reweight_lhef.h"

static int opdfset, opdfmem;
static int ipdfset, ipdfmem;
static int kf1, kf2;
static int nprup = 0;
static int idwtup;
static double e1, e2;
static double * cs;
static double * cserr;

int read_config (char configname[])
{
  char buff[2048];
  FILE * f = fopen (configname, "r");

  if (NULL == fgets (buff, 2048, f)) return -1;
  if (NULL == fgets (buff, 2048, f)) return -1;
  if (1 != sscanf (buff, "%i", &ipdfset)) return -1;

  if (NULL == fgets (buff, 2048, f)) return -1;
  if (NULL == fgets (buff, 2048, f)) return -1;
  if (1 != sscanf (buff, "%i", &ipdfmem))  return -1;

  if (NULL == fgets (buff, 2048, f)) return -1;
  if (NULL == fgets (buff, 2048, f)) return -1;
  if (1 != sscanf (buff, "%i", &opdfset))  return -1;

  if (NULL == fgets (buff, 2048, f)) return -1;
  if (NULL == fgets (buff, 2048, f)) return -1;
  if (1 != sscanf (buff, "%i", &opdfmem))  return -1;
  fclose (f);

  return 0;
}

#ifdef LHAPDF
static int init_lhapdf (int bnum, int set, int mem, int kfcode) {
  lhapdfList * list = NULL;
  lhapdfList * list_;
  int status = 0;

  /* In LHAPDF 6, set numbers are not used. Init with first available PDF set. */
  comphepLhapdfList (&list);
  list_ = list;

  if (list_) {
    initLHAPDF (bnum, list_->name, mem, kfcode);
    status = 1;
  }
  if (list) {
    delLhapdfList (list);
  }
  return status;
}

static void find_pairs (int num, double * x1, double * x2, int * n1, int * n2) {
  int i, j;
  int k = 0;
  double dx1, dx2, d;

  int * used = malloc ((num + 1) * sizeof (int));
  for (i = 0; i < num; ++i) {
    used[i] = 0;
  }

  for (i = 0; i < num; ++i) {
    if (used[i]) continue;
    int the_n = i;
    double dist = 10.;
    n1[k] = i;
    used[i] = 1;
    for (j = i + 1; j < num; ++j) {
      if (used[j]) continue;
      dx1 = x1[i] - x1[j];
      dx2 = x2[i] - x2[j];
      d = dx1 * dx1 + dx2 * dx2;
      if (d < dist) {dist = d; the_n = j;}
    }
    n2[k] = the_n;
    ++k;
    if (i > 0 && 0 == i % 10000) fprintf (stdout, "pdf_reweighter (info): %i events processed\n", i);
  }
}


int write_event_cap_lhef_cpyth (char fname[], char mode[]) {
  int idwtup = 3;
  long pos;
  int i;

  FILE * outFile = fopen (fname, mode);
  if (!outFile) {
    return 0;
  }

  fprintf (outFile, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (outFile, "<!-- File generated with CompHEP 4.5.0 -->\n");
  fprintf (outFile, "<!-- \n");
  fprintf (outFile, "     Preliminary version, it is compatible with the Les Houches event file\n");
  fprintf (outFile, "     format (hep-ph/0609017), but contains extra tags.\n");
  fprintf (outFile, "-->\n");
  fprintf (outFile, "<header>\n");

/*  fprintf (outFile, "%s", prepare_hepml_header ());*/
  fprintf (outFile, "</header>\n");
  fprintf (outFile, "<init>\n");
  fprintf (outFile, "%i %i %17.10E %17.10E %i %i %i %i %i %i\n",
                    kf1, 
                    kf2, 
                    e1, 
                    e2, 
                    opdfset, 
                    opdfset, 
                    opdfmem, 
                    opdfmem, 
                    idwtup, 
                    nprup);
  for (i = 0; i < nprup; ++i) {
  fprintf (outFile, "%17.10E %17.10E %17.10E %i\n", cs[i], cserr[i], 1.0, i);
  }
  fprintf (outFile, "</init>\n");

  pos = ftell (outFile);
  fclose (outFile);

  return pos;
}

int 
pdf_reweight_lhef (char config[], char source[], char target[], int nalphas, long the_seed, int error_used)
{
  int i;
  int num, nnum;
  double totcs = 0.;
  double totcserr = 0.;
  double * x1;
  double * x2;
  double * pdf_ratio;
  double pdf_ratio_cs;
  double pdf_ratio_cserr;
  int * epos;
  int * inumcount;
  int * onumcount;
  char buff[2048];
  FILE * s = fopen (source, "r");
  FILE * t;
  int set1, set2;
  int mem1, mem2;
  int sposition;
  long pos;
  eventUP evnt;
  double pdfmax = -10000.;
  double pdfmin =  10000.;

  if (read_config (config)) {
    fprintf(stderr,"pdf_reweighter (error): strange config file %s\n", config);
    return -4;
  }

  srand48 (the_seed);

  fgets (buff, 2048, s);
  while (!strstr (buff, "<init>")) {
    fgets (buff, 2048, s);
    if (feof (s)) return -4;
  }
  fgets (buff, 2048, s);
  if (10 != sscanf (buff, " %d %d %le %le %d %d %d %d %d %d", &kf1, &kf2, &e1, &e2, &set1, &set2, &mem1, &mem2, &idwtup, &nprup)) {
    return -4;
  }

  if (4 == set1) {
      switch (mem1) {
        case 34: {set1 = 19150; mem1 = 0; break;} /* cteq4m  */
        case 32: {set1 = 19170; mem1 = 0; break;} /* cteq4l  */
        case 53: {set1 = 19051; mem1 = 0; break;} /* cteq5m1 */
        case 48: {set1 = 19050; mem1 = 0; break;} /* cteq5m  */
        case 46: {set1 = 19070; mem1 = 0; break;} /* cteq5l  */
        case 55: {set1 = 10041; mem1 = 0; break;} /* cteq6l  */
/*        case 56: {set1 = ?????; mem1 = 0; break;}  cteq6d  */
        case 57: {set1 = 10050; mem1 = 0; break;} /* cteq6m  */
        case 58: {set1 = 10042; mem1 = 0; break;} /* cteq6l1 */
        default: {
          fprintf(stderr,"pdf_reweighter (error): can't detect PDF set/mem numbers in the event file (set = %i, mem = %i).\n", set1, mem1);
          return -1;
        }
      }
  }
  if (4 == set2) {
      switch (mem1) {
        case 34: {set2 = 19150; mem2 = 0; break;} /* cteq4m  */
        case 32: {set2 = 19170; mem2 = 0; break;} /* cteq4l  */
        case 53: {set2 = 19051; mem2 = 0; break;} /* cteq5m1 */
        case 48: {set2 = 19050; mem2 = 0; break;} /* cteq5m  */
        case 46: {set2 = 19070; mem2 = 0; break;} /* cteq5l  */
        case 55: {set2 = 10041; mem2 = 0; break;} /* cteq6l  */
/*        case 56: {set2 = ?????; mem2 = 0; break;}  cteq6d  */
        case 57: {set2 = 10050; mem2 = 0; break;} /* cteq6m  */
        case 58: {set2 = 10042; mem2 = 0; break;} /* cteq6l1 */
        default: {
          fprintf(stderr,"pdf_reweighter (error): can't detect PDF set/mem numbers in the event file (set = %i, mem = %i).\n", set1, mem1);
          return -1;
        }
      }
  }

  if (set1 == 0) set1 = ipdfset;
  if (set2 == 0) set2 = ipdfset;
  if (mem1 == 0) mem1 = ipdfmem;
  if (mem2 == 0) mem2 = ipdfmem;

  if (set1 != ipdfset || set2 != ipdfset || mem1 != ipdfmem || mem2 != ipdfmem) {
    fprintf (stdout, "pdf_reweighter (warning): PDF info taken from event file is different then one taken form config file\n");
    fprintf (stdout, "                          event  file: set1 =%i, mem1 - %i, set2 = %i, mem2 = %i\n", set1, mem1, set2, mem2);
    fprintf (stdout, "                          config file: set1 =%i, mem1 - %i, set2 = %i, mem2 = %i\n", ipdfset, ipdfmem, ipdfset, ipdfmem);
  }

  cs    = malloc (nprup * sizeof (double));
  cserr = malloc (nprup * sizeof (double));
  inumcount = malloc (nprup * sizeof (double));
  onumcount = malloc (nprup * sizeof (double));

  for (i = 0; i < nprup; ++i) {
  double wght;
  int nnum;
    fgets (buff, 2048, s);
    if (4 != sscanf (buff, " %le %le %le %i", &cs[i], &cserr[i], &wght, &nnum)) {
      return -4;
    }
    inumcount[i] = 0;
    onumcount[i] = 0;
    totcs += cs[i];
    totcserr += cserr[i] * cserr[i];
  }
  totcserr = sqrt (totcserr);

  fgets (buff, 2048, s);
  while (!strstr (buff, "</init>")) {
    fgets (buff, 2048, s);
    if (feof (s)) return -4;
  }
  sposition = ftell (s);

  num = 0;
  pos = sposition;
  while (!testLHAevent (s, source, pos)) {
    ++num;
    if (num > 0 && 0 == num % 10000) fprintf (stdout, "pdf_reweighter (info): %i events tested\n", num);
    pos = ftell (s);
  }
  fprintf (stdout, "pdf_reweighter (info): %i events read in total\n", num);
  pdf_ratio = malloc ((num + 1) * sizeof (double));
  epos      = malloc ((num + 1) * sizeof (int));

  x1 = malloc ((num + 1) * sizeof (double));
  x2 = malloc ((num + 1) * sizeof (double));

  if (0 == init_lhapdf (0, set1, mem1, kf1) || 0 == init_lhapdf (1, set2, mem2, kf2)) {
    fprintf (stdout, "pdf_reweighter (error): can\'t initialize LHAPDF for PDFset = %i, PDFmem = %i and beams %i and %i\n", 
    set1, mem1, kf1, kf2);
    return -1;
  } else {
    fprintf (stdout, "pdf_reweighter (info): old PDFset = %i, PDFmem = %i\n", set1, mem1);
  }

  if (0 < nalphas) {
    fprintf (stdout, "pdf_reweighter (info): alpha_s power %i is used in re-weighting\n", nalphas);
    set_alphaMode (1);
  }

  num = 0;
  epos[0] = sposition;
  while (0 == getLHAevent (source, s, epos[num], &evnt)) {
    x1[num] = fabs (evnt.momentumUP[2][0] / e1);
    x2[num] = fabs (evnt.momentumUP[2][1] / e2);
    pdf_ratio[num] = lhapdfValCPYTH (x1[num], evnt.QscaleUP, evnt.IDpartUP[0]) * lhapdfValCPYTH (x2[num], evnt.QscaleUP, evnt.IDpartUP[1]);
    if (0 < nalphas) {
      double alpha_s = alpha_2 (evnt.QscaleUP);
      pdf_ratio[num] *=  pow (alpha_s, (double)nalphas);
    }
    ++inumcount[evnt.IDprocUP - 1];
    ++num;
    epos[num] = ftell (s);
  }

  if (0 == init_lhapdf (1, opdfset, opdfmem, kf1) || 0 == init_lhapdf (1, opdfset, opdfmem, kf2)) {
    fprintf (stdout, "pdf_reweighter (err): can\'t initialize LHAPDF for PDFset = %i, PDFmem = %i and beams %i and %i\n", 
    opdfset, opdfmem, kf1, kf2);
    return -1;
  } else {
    fprintf (stdout, "pdf_reweighter (info): new PDFset = %i, PDFmem = %i\n", opdfset, opdfmem);
  }

  num = 0;
  pdf_ratio_cs = 0.;
  while (0 == getLHAevent (source, s, epos[num], &evnt)) {
    if (0. == pdf_ratio[num]) continue;
    pdf_ratio[num] = lhapdfValCPYTH (x1[num], evnt.QscaleUP, evnt.IDpartUP[0]) * lhapdfValCPYTH (x2[num], evnt.QscaleUP, evnt.IDpartUP[1]) / pdf_ratio[num];
    if (0 < nalphas) {
      double alpha_s = alpha_2 (evnt.QscaleUP);
      pdf_ratio[num] *= pow (alpha_s, (double)nalphas);
    }
    pdf_ratio_cs += pdf_ratio[num];
    ++num;
  }
  pdf_ratio_cs /= (double)num;

  nnum = num/2;
  if (error_used) {
    int * n1 = malloc ((nnum + 1) * sizeof (int));
    int * n2 = malloc ((nnum + 1) * sizeof (int));
    double sigma2 = 0.;
    fprintf (stdout, "pdf_reweighter (info): calculating error...\n");
    find_pairs (num, x1, x2, n1, n2);
    for (i = 0; i < nnum; ++i) {
      double dp = pdf_ratio[n1[i]] - pdf_ratio[n2[i]];
      sigma2 += dp * dp;
    }
    pdf_ratio_cserr = sqrt(sigma2 / 2.) / (double)nnum;
    fprintf (stdout, "pdf_reweighter (info): pdf_ratio_cs = %f, pdf_ratio_cserr = %f\n", pdf_ratio_cs, pdf_ratio_cserr);
  } else {
    fprintf (stdout, "pdf_reweighter (info): pdf_ratio_cs = %f\n", pdf_ratio_cs);
  }

  double pdfaver = 0;
  for (i = 0; i < num; ++i) {
    if (pdfmax < pdf_ratio[i]) pdfmax = pdf_ratio[i];
    if (pdfmin > pdf_ratio[i]) pdfmin = pdf_ratio[i];
    pdfaver += pdf_ratio[i];
  }
  fprintf (stdout, "pdf_reweighter (info): weight range: (%f, %f), aver. weight = %f\n", pdfmax, pdfmin, pdfaver / num);

  i = 0;
  write_event_cap_lhef_cpyth (target, "w");
  t = fopen (target, "a+");
  pos = ftell (t);
  while (0 == getLHAevent (source, s, epos[i], &evnt) && i < num) {
    double xrn = drand48 ();
    if (i > 0 && 0 == i % 10000) fprintf (stdout, "pdf_reweighter (info): %i events processed\n", num);
   if (xrn < pdf_ratio[i] / pdfmax) {
      setLHAevent (target, t, pos, &evnt);
      pos = ftell (t);
     ++onumcount[evnt.IDprocUP - 1];
    }
    ++i;
  }
  fclose (s);
  fclose (t);

/* some statistics */
  int nSel = 0;
  for (i = 0; i < nprup; ++i) {
    double k = (double)onumcount[i] / (double)inumcount[i];
    nSel += onumcount[i];
    cs[i] *= k;
    cserr[i] *= k;
  }
  fprintf (stdout, "pdf_reweighter (info): cross section before re-weighting = %f +/- %f\n", totcs, totcserr);
  if (error_used) {
    fprintf (stdout, "pdf_reweighter (info): cross section before re-weighting = %f +/- %f\n", pdf_ratio_cs * totcs, pdf_ratio_cserr * totcs);
  } else {
    fprintf (stdout, "pdf_reweighter (info): cross section before re-weighting = %f +/- unknown\n", pdf_ratio_cs * totcs);
  }
  fprintf (stdout, "pdf_reweighter (info): %i events selected from %i\n", nSel, num);

  write_event_cap_lhef_cpyth (target, "r+");

  free (epos);
  free (pdf_ratio);
  free (x1);
  free (x2);
  free (cs);
  free (cserr);
  free (inumcount);
  free (onumcount);

  return 0;
}
#endif
