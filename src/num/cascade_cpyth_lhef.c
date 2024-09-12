/* 
* Copyright (C) 2008-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ----------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/4_vector.h"

#include "LesHouches.h"
#include "event_reader.h"
#include "lhef_routines.h"
#include "cascade_cpyth_lhef.h"

static int fnl1[1024];
static int fnl2[1024];

static int 
form_new_color (int rpos, int nprd, int * prd1, int * prd2, int ndec, int * dec1, int * dec2) {
  int i;
  int max_prd = 0;
  int min_prd = 10000000;
  int min_dec = 10000000;
  int shift = nprd - 2;
  int ntot = nprd + ndec;
  int nfin = nprd + ndec - 2;
  int excl1 = 0;
  int excl2 = 0;

  if ((prd1[rpos] != 0 && dec1[0] == 0) || (prd1[rpos] == 0 && dec1[0] != 0) ||
      (prd2[rpos] != 0 && dec2[0] == 0) || (prd2[rpos] == 0 && dec2[0] != 0)) {
    return -1;
  }

  for (i = 0; i < ntot; ++i) {
   fnl1[i] = 0;
   fnl2[i] = 0;
  }

  for (i = 0; i < nprd; ++i) {
   if (max_prd < prd1[i]) max_prd = prd1[i];
   if (max_prd < prd2[i]) max_prd = prd2[i];
   if (min_prd > prd1[i] && 0 < prd1[i]) min_prd = prd1[i];
   if (min_prd > prd2[i] && 0 < prd2[i]) min_prd = prd2[i];
  }


  if (dec1[0] != 0) {
    for (i = 1; i < ndec; ++i) {
      if (dec1[0] == dec1[i]) dec1[i] = prd1[rpos] + 10000;
      if (dec1[0] == dec2[i]) dec2[i] = prd1[rpos] + 10000;
    }
    dec1[0] = prd1[rpos] + 10000;
  }

  if (dec2[0] != 0) {
    for (i = 1; i < ndec; ++i) {
      if (dec2[0] == dec1[i]) dec1[i] = prd2[rpos] + 10000;
      if (dec2[0] == dec2[i]) dec2[i] = prd2[rpos] + 10000;
    }
    dec2[0] = prd2[rpos] + 10000;
  }

  for (i = 0; i < ndec; ++i) {
   if (min_dec > dec1[i] && 0 < dec1[i]) min_dec = dec1[i];
   if (min_dec > dec2[i] && 0 < dec2[i]) min_dec = dec2[i];
  }
  for (i = 0; i < ndec; ++i) {
    if (0 < dec1[i] && dec1[i] < 10000) dec1[i] += max_prd - min_dec + 1;
    if (0 < dec2[i] && dec2[i] < 10000) dec2[i] += max_prd - min_dec + 1;
  }

  if (dec1[0] != 0) {
    for (i = 1; i < ndec; ++i) {
      if (dec1[0] == dec1[i]) dec1[i] -= 10000;
      if (dec1[0] == dec2[i]) dec2[i] -= 10000;
    }
    dec1[0] -= 10000;
  }

  if (dec2[0] != 0) {
    for (i = 1; i < ndec; ++i) {
      if (dec2[0] == dec1[i]) dec1[i] -= 10000;
      if (dec2[0] == dec2[i]) dec2[i] -= 10000;
    }
    dec2[0] -= 10000;
  }

/*
  if (prd1[rpos] != 0) {
    
    for (i = 1; i < ndec; ++i) {
      if (dec1[0] == dec1[i]) {
        excl1 = dec1[i];
        dec1[i] = prd1[rpos];
      }
    }
  }
  if (prd2[rpos] != 0) {
    for (i = 1; i < ndec; ++i) {
      if (dec2[0] == dec2[i]) {
        excl2 = dec2[i];
        dec2[i] = prd2[rpos];
      }
    }
  }
*/
  for (i = 0; i < nprd; ++i) {
   if (i < rpos) {
     fnl1[i] = prd1[i];
     fnl2[i] = prd2[i];
   }
   if (i > rpos) {
     fnl1[i - 1] = prd1[i];
     fnl2[i - 1] = prd2[i];
   }
  }

  for (i = 1; i < ndec; ++i) {
   fnl1[shift + i] = dec1[i];
   fnl2[shift + i] = dec2[i];
  }


  if (excl1 != 0) {
    for (i = 0; i < nfin; ++i) {
      if (excl1 < fnl1[i]) --fnl1[i];
      if (excl1 < fnl2[i]) --fnl2[i];
    }
  }
  if (excl2 != 0) {
    for (i = 0; i < nfin; ++i) {
      if (excl2 < fnl1[i]) --fnl1[i];
      if (excl2 < fnl2[i]) --fnl2[i];
    }
  }

  return 0;
}

int cascade_combine (eventUP * prd, eventUP * dec, eventUP * fnl)
{
  int i, j, k;
  int real_decay = 0;
  int shift = prd->NpartUP - 2;
  int res_pos = 0;
  double newrf[4];

  fnl->NpartUP    = prd->NpartUP + dec->NpartUP - 2;
  fnl->IDprocUP   = prd->IDprocUP;
  fnl->XweightUP  = prd->XweightUP;
  fnl->QscaleUP   = prd->QscaleUP;
  fnl->QEDalphaUP = prd->QEDalphaUP;
  fnl->QCDalphaUP = prd->QCDalphaUP;

  j = 0;
  for (i = 0; i < 2; ++i) {
    fnl->IDpartUP[j]    = prd->IDpartUP[i];
    fnl->motherUP[0][j] = prd->motherUP[0][i];
    fnl->motherUP[1][j] = prd->motherUP[1][i];
    fnl->timelifeUP[j]  = prd->timelifeUP[i];
    fnl->spinUP[j]      = prd->spinUP[i];
    fnl->statusUP[j]    = prd->statusUP[i];
    for (k = 0; k < 5; ++k) {
      fnl->momentumUP[k][j] = prd->momentumUP[k][i];
    }
    ++j;
  }

  for (i = 2; i < prd->NpartUP; ++i) {
    if (prd->IDpartUP[i] != dec->IDpartUP[0] || real_decay) {
      fnl->IDpartUP[j]    = prd->IDpartUP[i];
      fnl->motherUP[0][j] = prd->motherUP[0][i];
      fnl->motherUP[1][j] = prd->motherUP[1][i];
      fnl->timelifeUP[j]  = prd->timelifeUP[i];
      fnl->spinUP[j]      = prd->spinUP[i];
      fnl->statusUP[j]    = prd->statusUP[i];
      for (k = 0; k < 5; ++k) {
        fnl->momentumUP[k][j] = prd->momentumUP[k][i];
      }
      ++j;
    } else {
      res_pos = i;
      real_decay = 1;
      newrf[1] = -prd->momentumUP[0][i];
      newrf[2] = -prd->momentumUP[1][i];
      newrf[3] = -prd->momentumUP[2][i];
      newrf[0] = prd->momentumUP[3][i];
    }
  }

  form_new_color (res_pos,  prd->NpartUP, prd->colorflowUP[0], prd->colorflowUP[1], dec->NpartUP, dec->colorflowUP[0], dec->colorflowUP[1]);
  for (i = 0; i < fnl->NpartUP; ++i) {
    fnl->colorflowUP[0][i] = fnl1[i];
    fnl->colorflowUP[1][i] = fnl2[i];
  }

  for (i = 1;  i < dec->NpartUP; ++i) {
    double a[4];
    double b[4];
    fnl->IDpartUP[shift + i]    = dec->IDpartUP[i];
    fnl->motherUP[0][shift + i] = 1;
    fnl->motherUP[1][shift + i] = 2;
    fnl->timelifeUP[shift + i]  = dec->timelifeUP[i];
    fnl->spinUP[shift + i]      = dec->spinUP[i];
    fnl->statusUP[shift + i]    = dec->statusUP[i];

    a[1] = dec->momentumUP[0][i];
    a[2] = dec->momentumUP[1][i];
    a[3] = dec->momentumUP[2][i];
    a[0] = dec->momentumUP[3][i];
    lorenc (a, newrf, b);
    fnl->momentumUP[0][shift + i] = b[1];
    fnl->momentumUP[1][shift + i] = b[2];
    fnl->momentumUP[2][shift + i] = b[3];
    fnl->momentumUP[3][shift + i] = b[0];
    fnl->momentumUP[4][shift + i] = dec->momentumUP[4][i];
  }

  return 0;
}


int
cascade_cpyth_lhef (int regime, char prd_name[], char dec_name[], char out_name[]) 
{
  int i;
  int decnum;
  int dec_id;
  int prd_evnt_used = 0;
  int dec_evnt_used = 0;
  int dec_evnt_num = 1;
  int prd_evnt_num = 0;
  long dec_position;
  long fnl_position;
  long tmp_file_position;
  double xwgt;
  char buff[2048];
  char * stop = NULL;
  eventUP dec_evnt;
  eventUP prd_evnt;
  double dec_rate[MAX_FILE_EVENT];
  double dec_rateerr[MAX_FILE_EVENT];
  analyzeMultiProcLHEfile (prd_name);
  long prd_position = getEventPosition (0);
  int nf = getTotProcNumber ();

  FILE * pfile = fopen (prd_name, "r");
  FILE * dfile = fopen (dec_name, "r");
  FILE * ofile  = fopen (out_name, "w+");

  if (NULL == pfile || NULL == dfile || NULL == ofile) {
    return -1;
  }

/*********************************************************************************************************/
  while (!stop) {
    fgets (buff, 2048, dfile);
    if (feof (dfile)) return -4;
    stop = strstr (buff, "<init>");
  }
  fgets (buff, 2048, dfile);
  fgets (buff, 2048, dfile);
  stop = NULL;
  decnum = 0;
  while (!stop) {
    if (4 != sscanf (buff, "%le %le %le %d", &dec_rate[decnum], &dec_rateerr[decnum], &xwgt, &i)) {
      dec_rate[decnum] = 0.0;
      dec_rateerr[decnum] = 0.0;
    }
    ++decnum;
    fgets (buff, 2048, dfile);
    stop = strstr (buff, "</init>");
  }

  stop = NULL;
  while (!stop) {
    dec_position = ftell (dfile);
    fgets (buff, 2048, dfile);
    if (feof (dfile)) return -4;
    stop = strstr (buff, "<event>");
  }

  if (0 == getLHAevent (dec_name, dfile, dec_position, &dec_evnt)) {
    dec_id = dec_evnt.IDpartUP[0];
  }

  tmp_file_position = dec_position;
  while (0 == getLHAevent (dec_name, dfile, tmp_file_position, &dec_evnt)) {
    tmp_file_position = ftell (dfile);
    int loc_dec_id = dec_evnt.IDpartUP[0];
    if (dec_id != loc_dec_id) {
      fprintf (stderr, "cascade (error): strange decay event file. ");
      fprintf (stderr, "                 It should contain decays of ONE particle! Check the file. Exit\n");
      return -1;
    }
    ++dec_evnt_num;
  }

  tmp_file_position = prd_position;
  while (0 == getLHAevent (prd_name, pfile, tmp_file_position, &prd_evnt)) {
    tmp_file_position = ftell (pfile);
    ++prd_evnt_num;
  }
  fprintf (stdout, "cascade (info): production event file %s contains %i events\n", prd_name, prd_evnt_num);
  fprintf (stdout, "cascade (info): decay event file %s contains %i events\n", dec_name, dec_evnt_num);

/*********************************************************************************************************/
  fprintf (ofile, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (ofile, "<!-- File generated with CompHEP %s -->\n", getCHEPversion ());
#ifdef LIBXML
  fprintf (ofile, "<!-- \n"
                    "     This file is compatible with the Les Houches event file\n"
                    "     format (hep-ph/0609017), but contains extra HepML tags.\n"
                    "-->\n");
  fprintf (ofile, "<header>\n");
  fprintf (ofile, "%s", prepare_hepml_header_libxml2_dynamic ());
#else
  fprintf (ofile, "<!-- \n"
                    "     This file is compatible with the Les Houches event file\n"
                    "     format (hep-ph/0609017).\n"
                    "-->\n");
  fprintf (ofile, "<header>\n");
#endif
  fprintf (ofile, "</header>\n");
  fprintf (ofile, "<init>\n");
  fprintf (ofile, "%i %i %17.10E %17.10E %i %i %i %i 3 %i\n",
                    getPbeam (0), 
                    getPbeam (1), 
                    getEbeam (0), 
                    getEbeam (1), 
                    getPDFLIBgroup (0), 
                    getPDFLIBgroup (1), 
                    getPDFLIBset (0), 
                    getPDFLIBset (1),
                    nf);
  for (i = 0; i < nf; ++i) {
    fprintf (ofile, "%17.10E %17.10E %17.10E %i\n", getCrossSection (i), getCrossSectionErr (i), 1.0, i + 1);
  }
  fprintf (ofile, "</init>\n");
  fnl_position = ftell (ofile);

  prd_evnt_used = 0;
  dec_evnt_used = 0;
  while (0 == getLHAevent (prd_name, pfile, prd_position, &prd_evnt)) {
    int new_event_written = 0;
    eventUP fnl_evnt;
    prd_position = ftell (pfile);
    ++prd_evnt_used;
    for (i = 0; i < prd_evnt.NpartUP; ++i) {
      if (dec_id == prd_evnt.IDpartUP[i]) {
        if (0 == getLHAevent (dec_name, dfile, dec_position, &dec_evnt)) {
          dec_position = ftell (dfile);
          ++dec_evnt_used;
	} else {
          fprintf (stdout, "cascade (info): %i production events have been treated\n", prd_evnt_used);
          fprintf (stdout, "cascade (info): %i decay events have been used\n", dec_evnt_used);
          if (prd_evnt_num > prd_evnt_used) fprintf (stdout, "cascade (warning): "
          "be aware! not all production event have decayed resonances! Tot num = %i, Proc num = %i\n", prd_evnt_num, prd_evnt_used);
          return 0;
        }
        cascade_combine (&prd_evnt, &dec_evnt, &fnl_evnt);
        setLHAevent (out_name, ofile, fnl_position, &fnl_evnt);
        fnl_position = ftell (ofile);
        new_event_written = 1;
        if (0 == regime) break;
      }
    }
    if (0 == new_event_written) {
      setLHAevent (out_name, ofile, fnl_position, &prd_evnt);
      fnl_position = ftell (ofile);
    }
  }
  fprintf (stdout, "cascade (info): %i production events have been treated\n", prd_evnt_used);
  fprintf (stdout, "cascade (info): %i decay events have been used\n", dec_evnt_used);
  if (prd_evnt_num > prd_evnt_used) fprintf (stdout, "cascade (warning): "
  "be aware! not all production event have decayed resonances! Tot num = %i, Proc num = %i\n", prd_evnt_num, prd_evnt_used);
  fclose (pfile);
  fclose (dfile);
  fclose (ofile);

  return 0;
}
