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

#include "service2/include/chep_limits.h"
#include "LesHouches.h"
#include "userFun.h"
#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_routines.h"
#include "event_reader.h"

#include "addcut_lhef.h"

int
addcut_lhef (char ini_name[], char out_name[], int nevnt)
{
  long nposition = 0;
  long tposition = 0;
  int i;
  int Nevents = 0;
  int Nevents_acc = 0;
  eventUP ev;
  double rate = 0.;
  double rate_err = 0.;
  double tmprate = 0.;
  double tmprate_err = 0.;
  double xwgt = 0.;
  char buff[2048];

  FILE * s = fopen (ini_name, "r");
  FILE * t = fopen (out_name, "w");

  fgets (buff, 2048, s);
  while (!strstr (buff, "<init>")) {
    fputs (buff, t);
    fgets (buff, 2048, s);
    if (feof (s)) return -4;
  }

  fgets (buff, 2048, s);
  fputs (buff, t);
  fgets (buff, 2048, s);
  fputs (buff, t);

  while (!strstr (buff, "</init>")) {
    if (4 != sscanf (buff, "%le %le %le %d", &tmprate, &tmprate_err, &xwgt, &i)) {
    tmprate = 0.;
    tmprate_err = 0.;
// message should be issued
    }
    rate += tmprate;
    rate_err += tmprate_err * tmprate_err;
    fputs (buff, t);
    fgets (buff, 2048, s);
  }
  rate_err = sqrt (rate_err);
  tposition = ftell (t);
  nposition = ftell (s);

  while (0 == getLHAevent (ini_name, s, nposition, &ev)) {
    double sel = cutFunction (&ev);
    nposition = ftell (s);
    if (0 == Nevents % 10000) fprintf (stdout, "addcut (info): %i events processed\n", Nevents);
    ++Nevents;
    if (sel) {
      setLHAevent (out_name, t, tposition, &ev);
      tposition = ftell (t);
      ++Nevents_acc;
    }
  }
  fclose (s);

  double ratio = (double)Nevents_acc / (double) Nevents;
  fprintf (stdout, "addcut (info): %i events processed, %i events stored\n", Nevents, Nevents_acc);
  fprintf (stdout, "addcut (info): cs_ini = %f, cs_fin = %f\n", rate, rate * ratio);
  return 0;
}
