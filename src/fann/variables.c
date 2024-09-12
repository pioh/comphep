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
#include "service2/include/tptcmac.h"
#include "num/include/LesHouches.h"
#include "num/include/event_reader.h"

#include "userVars.h"
#include "variables.h"

int
prepare_fann_variables (char ini_name[], char out_name[], int nevnt, int nvarset)
{
  long nposition;
  long evini;
  int Nevents = 0;
  eventUP ev;
  char buff[2048];
  char * stop = NULL;
  double vrs[1024];
  int nvrs = 0;
  char * fnline;

  FILE * s = fopen (ini_name, "r");
  FILE * f = fopen (out_name, "w");

  if (-1 == nevnt) nevnt = 1000000000;
  stop = NULL;
  while (!stop) {
    nposition = ftell (s);
    fgets (buff, 2048, s);
    if (feof (s)) return -4;
    stop = strstr (buff, "<event>");
  }
  evini = nposition;

  fprintf (stdout, "fannvars (info): testing file %s:\n", ini_name);
  while (!testLHAevent (s, ini_name, nposition) && Nevents < nevnt) {
    nposition = ftell (s);
    ++Nevents;
    if (0 == Nevents % 20000) fprintf (stdout, "fannvars (info): %i events tested\n", Nevents);
  }
  fprintf (stdout, "fannvars (info): %i events tested in total\n", Nevents);
  nevnt = Nevents++;

  nposition = evini;
  if (0 == getLHAevent (ini_name, s, nposition, &ev)) {
    int i;
    nposition = ftell (s);
    userVariables (&ev, nvarset, &nvrs, vrs);
    if (NULL == fnline) {
      fnline = malloc (13 * (nvrs + 2) * sizeof (char));
    }
    fprintf (f, "%i %i %i\n", nevnt, nvrs, 1);
    for (i = 0; i < nvrs; ++i) {
      strcat (fnline, scat ("%12.5e ", vrs[i]));
    }
    strcat (fnline, "\n");
    fputs (fnline, f);
    fprintf (f, "%12.5e\n", ev.XweightUP);
    Nevents = 1;

    while (0 == getLHAevent (ini_name, s, nposition, &ev)) {
      nposition = ftell (s);
      userVariables (&ev, nvarset, &nvrs, vrs);
      fnline[0] = '\0';
      for (i = 0; i < nvrs; ++i) {
        strcat (fnline, scat ("%12.5e ", vrs[i]));
      }
      strcat (fnline, "\n");
      fputs (fnline, f);
      fprintf (f, "%12.5e\n", ev.XweightUP);
      ++Nevents;
      if (0 == Nevents % 10000) fprintf (stdout, "fannvars (info): %i events read\n", Nevents);
      if (nevnt > 0 && nevnt <= Nevents) {
        break;
      }
    }
    free (fnline);
    fclose (f);
  } else {
    fprintf (stdout, "rtupler (error): can't read events\n");
    return -1;
  }

  fprintf (stdout, "rtupler (info): %i events have been written to the rtuple %s\n", Nevents, out_name);
  return 0;
}

int
add_fann_variables (char ini_name[], char out_name[], int nevnt, int nvarset)
{
  long nposition = 0L;
  long fposition = 0L;
  long evini;
  int Nevents = 0;
  eventUP ev;
  char buff[2048];
  char * stop = NULL;
  double vrs[1024];
  int nvrs = 0;
  char fnline[2048];

  FILE * s = fopen (ini_name, "r");
  FILE * f = fopen (out_name, "w");

  if (-1 == nevnt) nevnt = 1000000000;
  stop = NULL;
  while (!stop) {
    nposition = ftell (s);
    fposition = ftell (f);
    fgets (buff, 2048, s);
    if (feof (s)) return -4;
    fputs (buff, f);
    stop = strstr (buff, "<event>");
  }
  evini = nposition;

  fprintf (stdout, "fannvars (info): testing file %s:\n", ini_name);
  while (!testLHAevent (s, ini_name, nposition) && Nevents < nevnt) {
    nposition = ftell (s);
    ++Nevents;
    if (0 == Nevents % 20000) fprintf (stdout, "fannvars (info): %i events tested\n", Nevents);
  }
  fprintf (stdout, "fannvars (info): %i events tested in total\n", Nevents);
  nevnt = Nevents++;

  nposition = evini;
  if (0 == getLHAevent (ini_name, s, nposition, &ev)) {
    int i;
    nposition = ftell (s);
    userVariables (&ev, nvarset, &nvrs, vrs);
    fnline[0] = 0;
    strcat (fnline, scat ("#NNnvar %i\n#NNvars ", nvrs));
    for (i = 0; i < nvrs; ++i) {
      strcat (fnline, scat ("%12.5e ", vrs[i]));
    }
    strcat (fnline, scat ("\n#NNwgt %12.5e\n", 1.));
    setLHAevent_with_comments (out_name, f, fposition, &ev, fnline);
    Nevents = 1;

    while (0 == getLHAevent (ini_name, s, nposition, &ev)) {
      nposition = ftell (s);
      fposition = ftell (f);
      setLHAevent (out_name, f, fposition, &ev);
      userVariables (&ev, nvarset, &nvrs, vrs);
      fnline[0] = 0;
      strcat (fnline, scat ("#NNnvar %i\n#NNvars ", nvrs));
      for (i = 0; i < nvrs; ++i) {
        strcat (fnline, scat ("%12.5e ", vrs[i]));
      }
      strcat (fnline, scat ("\n#NNwgt %12.5e\n", 1.));
      setLHAevent_with_comments (out_name, f, fposition, &ev, fnline);
      ++Nevents;
      if (0 == Nevents % 10000) fprintf (stdout, "fannvars (info): %i events written\n", Nevents);
      if (nevnt > 0 && nevnt <= Nevents) {
        break;
      }
    }
    fclose (f);
  } else {
    fprintf (stdout, "rtupler (error): can't read events\n");
    return -1;
  }

  fprintf (stdout, "rtupler (info): %i events have been written to the rtuple %s\n", Nevents, out_name);
  return 0;
}
