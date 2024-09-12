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

extern float rnnfun(double *rin);

int
unweight_events (char ini_name[], char out_name[], int nevnt, long the_seed)
{
  int i;
  long nposition = 0L;
  long fposition = 0L;
  long evini;
  int Nevents = 0;
  int Nsel = 0;
  eventUP ev;
  char buff[2048];
  double max_weight = -1.;
  double avr_weight = 0.;

  FILE * s = fopen (ini_name, "r");
  FILE * f = fopen (out_name, "w");

  srand48 (the_seed);

  if (-1 == nevnt) nevnt = 1000000000;
  fgets (buff, 2048, s);
  while (!strstr (buff, "<event>")) {
    nposition = ftell (s);
    fposition = ftell (f);
    fputs (buff, f);
    if (feof (s)) return -4;
    fgets (buff, 2048, s);
  }
  evini = nposition;

  fprintf (stdout, "unweighter (info): testing file %s:\n", ini_name);
  while (0 == getLHAevent (ini_name, s, nposition, &ev) && Nevents < nevnt) {
    nposition = ftell (s);
    ++Nevents;
    if (max_weight < ev.XweightUP) {
      if (ev.XweightUP < 1.0) max_weight = ev.XweightUP; else max_weight = 1.;
    }
    avr_weight += ev.XweightUP;
    if (0 == Nevents % 20000) fprintf (stdout, "unweighter (info): %i events tested\n", Nevents);
  }
  nevnt = Nevents++;
  fprintf (stdout, "unweighter (info): %i events tested in total, found max = %f, av. weight = %f\n", nevnt, max_weight, avr_weight / (double)Nevents);

  nposition = evini;
  i = 0;
  while (0 == getLHAevent (ini_name, s, nposition, &ev)) {
    double xrn = drand48 ();
    nposition = ftell (s);
    fposition = ftell (f);
    ++i;
    if (ev.XweightUP / max_weight > xrn) {
      char * ev_comments = get_event_comments ();
      ev.XweightUP = 1.;
      setLHAevent (out_name, f, fposition, &ev);
/*      setLHAevent_with_comments (out_name, f, fposition, &ev, ev_comments);*/
      ++Nsel;
    }
  }

  fprintf (stdout, "rtupler (info): %i events have been selected from %i in the file %s\n", Nsel, Nevents, out_name);
  return 0;
}

double 
parse_comments (char cmts[])
{
  double wgt = 0.;
  double thewgt = 0.;
  char * pos = strtok (cmts, "#");
  while (pos)
    {
      if (1 == sscanf (cmts, "NNwgt %lf", &thewgt)) {
        wgt = thewgt;
        break;
      }
      pos = strtok (NULL, " ");
    }

  return wgt;
}

int
unweight_events_with_nn_weights (char ini_name[], char out_name[], int nevnt, long the_seed)
{
  long nposition = 0L;
  long fposition = 0L;
  long evini;
  int Nevents = 0;
  int Nsel = 0;
  eventUP ev;
  char buff[2048];
  double max_weight = -1.;
  double avr_weight = 0.;

  FILE * s = fopen (ini_name, "r");
  FILE * f = fopen (out_name, "w");

  srand48 (the_seed);

  if (-1 == nevnt) nevnt = 1000000000;
  fgets (buff, 2048, s);
  while (!strstr (buff, "<event>")) {
    nposition = ftell (s);
    fposition = ftell (f);
    fputs (buff, f);
    if (feof (s)) return -4;
    fgets (buff, 2048, s);
  }
  evini = nposition;

  fprintf (stdout, "unweighter (info): testing file %s:\n", ini_name);
  max_weight = -1.;
  avr_weight = 0.;
  while (0 == getLHAevent (ini_name, s, nposition, &ev) && Nevents < nevnt) {
    char * ev_comments = get_event_comments ();
    double nn_wight = parse_comments (ev_comments);
    nposition = ftell (s);
    ++Nevents;
    if (max_weight < nn_wight) max_weight = nn_wight;
    avr_weight += nn_wight;
    if (0 == Nevents % 20000) fprintf (stdout, "unweighter (info): %i events tested\n", Nevents);
  }
  fprintf (stdout, "unweighter (info): %i events tested in total, found max = %f, av. weight = %f\n", nevnt, max_weight, avr_weight / (double)Nevents);
  nevnt = Nevents++;

  nposition = evini;
  while (0 == getLHAevent (ini_name, s, nposition, &ev)) {
    double xrn = drand48 ();
    char * ev_comments = get_event_comments ();
    double nn_wight = parse_comments (ev_comments);
    nposition = ftell (s);
    fposition = ftell (f);
    if (nn_wight / max_weight > xrn) {
      setLHAevent_with_comments (out_name, f, fposition, &ev, ev_comments);
      ++Nsel;
    }
  }

  fprintf (stdout, "rtupler (info): %i events have been selected from %i in the file %s\n", Nsel, Nevents, out_name);

  return 0;
}

static int NNnvars = 0;
int 
test_comments (char cmts[], double nnvars[])
{
  int i;
  double v = 0.;
  char postvals[4096];
  char prevals[4096];

  char * pos = strtok (cmts, "#");
  while (pos)
    {
      if (strstr (pos, "NNnvar")) {
	if (1 == sscanf (pos, "NNnvar %d", &i)) {
          if (0 == NNnvars) NNnvars = i;
	  else {
	    if (NNnvars != i) return 0;
	  }
        }
      }
      if (strstr (pos, "NNvars")) {
        sscanf (pos, "NNvars%[^a]", prevals);
        for (i = 0; i < NNnvars && 0 < strlen (prevals); ++i) {
          if (2 == sscanf (prevals, " %lf%[^a]", &v, postvals)) nnvars[i] = v;
          else return 0;
          strcpy (prevals, postvals);
        }
      }
      pos = strtok (NULL, "#");
    }

  return 1;
}

double 
calculate_nnwgt (double vars[])
{
  double nnval = (double) rnnfun(vars);
  return nnval;
}

int
calculate_nn_weights_and_unweight (char ini_name[], char out_name[], int nevnt, long the_seed)
{
  int i;
  int Nevents = 1;
  int Nsel = 0;
  long nposition = 0L;
  long fposition = 0L;
  long evini;
  eventUP ev;
  char buff[2048];
  double max_weight = -1.;
  double avr_weight = 0.;
  double vars[1024];
  double * nnwgts;
  double * mean;
  double * sigma;

  FILE * s = fopen (ini_name, "r");
  FILE * f = fopen (out_name, "w");

  srand48 (the_seed);

  if (-1 == nevnt) nevnt = 1000000000;
  fgets (buff, 2048, s);
  while (!strstr (buff, "<event>")) {
    nposition = ftell (s);
    fposition = ftell (f);
    fputs (buff, f);
    if (feof (s)) return -4;
    fgets (buff, 2048, s);
  }
  evini = nposition;

/* 1st pass: test events and define the number of events */
  if (0 == getLHAevent (ini_name, s, nposition, &ev) && Nevents < nevnt) {
    char * ev_comments = get_event_comments ();
    if (0 == test_comments (ev_comments, vars)) {
      fprintf (stdout, "unweighter (error): problem in event %d. all events should have NN variable for the option! Exit.", Nevents + 1);
      return -1;
    }
  }
  nposition = ftell (s);
  mean = malloc (NNnvars * sizeof (double));
  sigma = malloc (NNnvars * sizeof (double));
  for (i = 0; i < NNnvars ; ++i) {
    mean[i] = vars[i];
    sigma[i] = vars[i] * vars[i];
  }
  fprintf (stdout, "unweighter (info): testing file %s:\n", ini_name);
  while (0 == getLHAevent (ini_name, s, nposition, &ev) && Nevents < nevnt) {
    char * ev_comments = get_event_comments ();
    if (0 == test_comments (ev_comments, vars)) {
      fprintf (stdout, "unweighter (error): problem in event %d. all events should have NN variable for the option! Exit.", Nevents + 1);
      return -1;
    }
    for (i = 0; i < NNnvars ; ++i) {
      mean[i] += vars[i];
      sigma[i] += vars[i] * vars[i];
    }
    nposition = ftell (s);
    ++Nevents;
    if (0 == Nevents % 20000) fprintf (stdout, "unweighter (info): %i events tested\n", Nevents);
  }
  fprintf (stdout, "unweighter (info): %i events found\n", Nevents);
    for (i = 0; i < NNnvars ; ++i) {
      mean[i] = mean[i] / Nevents;
      sigma[i] = sqrt (sigma[i] / Nevents - mean[i] * mean[i]);
    }

/* 2nd pass: calculating NN weights and search for the max weight */
  nnwgts = malloc (Nevents * sizeof (double));
  fprintf (stdout, "unweighter (info): calculation NN weights\n");
  max_weight = -1.;
  avr_weight = 0.;
  nposition = evini;
  Nevents = 0;
  while (0 == getLHAevent (ini_name, s, nposition, &ev) && Nevents < nevnt) {
    char * ev_comments = get_event_comments ();
    test_comments (ev_comments, vars);

    nnwgts[Nevents] = calculate_nnwgt (vars);
    nposition = ftell (s);
    if (max_weight < nnwgts[Nevents]) 
      max_weight = nnwgts[Nevents];
    avr_weight += nnwgts[Nevents];
    ++Nevents;
    if (0 == Nevents % 20000) fprintf (stdout, "unweighter (info): %i NN weights calculated\n", Nevents);
  }
  fprintf (stdout, "unweighter (info): %i events tested in total, found max = %f, av. weight = %f\n", Nevents, max_weight, avr_weight / (double)Nevents);

/* 3rd pass: event unweighting and writing down to the final file */
  nposition = evini;
  Nevents = 0;
  while (0 == getLHAevent (ini_name, s, nposition, &ev)) {
    double xrn = drand48 ();
    char * ev_comments = get_event_comments ();
    nposition = ftell (s);
    if (nnwgts[Nevents] / max_weight > xrn) {
      fposition = ftell (f);
      ev.XweightUP = 1.;
      setLHAevent (out_name, f, fposition, &ev);
/*      setLHAevent_with_comments (out_name, f, fposition, &ev, ev_comments);*/
      ++Nsel;
    }
    ++Nevents;
    if (0 == Nevents % 20000) fprintf (stdout, "unweighter (info): %i NN weights calculated\n", Nevents);
  }
  fprintf (stdout, "rtupler (info): %i events have been selected from %i in the file %s\n", Nsel, Nevents, out_name);
  free (nnwgts);
  free (mean);
  free (sigma);

  return 0;
}
