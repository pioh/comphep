/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <stdio.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "e_tools.h"

#define const
#include "out_ext.h"
#undef const

int nin_  = 0;
int nout_ = 0;
double sum = 0.;
static char p_names[MAXNP][20];
static double p_masses[MAXNP];


int
pinf_ (int nsub, int num, char *name, double *mass)
{
  if (name)
    strcpy (name, p_names[num - 1]);
  if (mass)
    *mass = p_masses[num - 1];
  return 0;
}


void
boost (double *n, double *p)
{
  double M = p[0];		/* particle mass */
  double shY = n[0];		/* h-sine of rapidity */
  double chY = sqrt (1 + shY * shY);
  double f =
    ENERGY (M,
	    p + 1) * shY + (p[1] * n[1] + p[2] * n[2] + p[3] * n[3]) * (chY -
									1);
  int i;

  for (i = 1; i <= 3; i++)
    p[i] += n[i] * f;
}


void
findBoost (double *p, double *n)
{
  double p2 = p[1] * p[1] + p[2] * p[2] + p[3] * p[3];
  double M = p[0];
  int i;

  if (p2 == 0)
    {
      for (i = 0; i <= 3; i++)
	n[i] = 0;
      return;
    }
  p2 = sqrt (p2);
  n[0] = p2 / M;
  for (i = 1; i <= 3; i++)
    n[i] = p[i] / p2;
}


int
getMasses (FILE * flow)
{
  char buff[STRSIZ];
  char *pos;
  int i = 0;

  fgets (buff, STRSIZ, flow);
  pos = strtok (buff, " ");
  while (pos)
    {
      sscanf (pos, "%lf", p_masses + i++);
      pos = strtok (NULL, " ");
    }
  return 0;
}

int
skipHeadLine (FILE * flow)
{
  fscanf (flow, "%*[^\n]\n");
  return -1;
}

int
skipLine (FILE * flow)
{
  fscanf (flow, "%*[^\n]\n");
  return 0;
}

static double CS = 0;
double getCS (void)
{
  return CS;
}

int
readCS (FILE * flow)
{
  fscanf (flow, "%lf", &CS);
  return 0;
}

static double totCS = 0.0;
double getTotCS (void)
{
  return totCS;
}

int
readTotCS (FILE * flow)
{
  char buff[1024];

  fgets (buff, 1024, flow);
  sscanf (buff, "= %lf", &totCS);
  return 0;
}

static long nEvents = 0;
int getNevents (void)
{
  return nEvents;
}

int
readNEvents (FILE * flow)
{
  fscanf (flow, "%ld", &nEvents);
  return 0;
}

int 
getNames (char * buff)
{
  int ntot_ = 0;
  char * buf;
  char * pch;
  char * temp = malloc (strlen (buff) * sizeof (char));

  strcpy (temp, buff);
  pch = strtok (buff, " ,->");
  while (pch != NULL) {
    strcpy (p_names[ntot_], pch);
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
  while (pch != NULL) {
    ++nout_;
    pch = strtok (NULL, " ,");
  }
  nin_ = ntot_ - nout_;

  return 0;
}


int
old_getNames (FILE * flow)
{
  int i;
  int k;
  int arrow = 0;
  char buff[STRSIZ];

  nin_ = 0;
  nout_ = 0;

  fgets (buff, STRSIZ, flow);
  for (i = 0, k = 0; buff[i] && (buff[i] != '\n');)
    {
      if (buff[i] == ' ')
	{
	  i++;
	  continue;
	}
      sscanf (buff + i, "%s", p_names[k]);
      if (strcmp (p_names[k], "->") == 0)
	{
	  arrow = 1;
	  i += 2;
	}
      else
	{
	  if (arrow)
	    nout_++;
	  else
	    nin_++;
	  i += strlen (p_names[k]);
	  k++;
	}
    }
  return 0;
}

int
decay2 (double M, double *p1, double *p2)
{
  double m1 = p1[0];
  double m2 = p2[0];
  double cosX, sinX, phi;
  double P;
  int i;

  if (m1 + m2 >= M)
    return 1;

  P =
    sqrt ((M - m1 - m2) * (M + m1 + m2) * (M - m1 + m2) * (M - m2 +
							   m1)) / (2 * M);

  cosX = 1 - 2 * drand48 ();
  sinX = sqrt (1 - cosX * cosX);

  phi = 2 * M_PI * drand48 ();

  p1[3] = cosX;
  p1[2] = sinX * cos (phi);
  p1[1] = sinX * sin (phi);

  for (i = 1; i <= 3; i++)
    {
      p1[i] = P * p1[i];
      p2[i] = -p1[i];
    }
  return 0;
}
