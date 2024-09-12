/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "service2/include/4_vector.h"
#include "num/include/LesHouches.h"
#include "userFun.h"


double 
userFunction (char *name)
{
/*
    The following commented line is a example of use userFunction.
    It define user variable "UMtAv" that calculate in the process
    of MC integration. This variable can use in the cut-table.
    More details see in Manual :-)
    
*/
/*
  double pt3, pt4, pt5, pt6, qsc;

  pt3 = pvect[8 + 1] * pvect[8 + 1] + pvect[8 + 2] * pvect[8 + 2];
  pt4 = pvect[12 + 1] * pvect[12 + 1] + pvect[12 + 2] * pvect[12 + 2];
  pt5 = pvect[16 + 1] * pvect[16 + 1] + pvect[16 + 2] * pvect[16 + 2];
  pt6 = pvect[20 + 1] * pvect[20 + 1] + pvect[20 + 2] * pvect[20 + 2];
  qsc = sqrt ((175.0 * 175.0) + (pt3 + pt4 + pt5 + pt6) / 4.0);

  if (strcmp (name, "MtAv") == 0)
    {
      if (!finite (qsc))
	fprintf (stdout, " uvar= %18.15g\n", qsc);
      return qsc;
    }
*/

  fprintf (stdout, " ***** userFunction: not defined: %s\n", name);
  exit (10);
  return 0.;
}

/*
typedef struct eventUP
  {
    int NpartUP;                    // number of particle entries in this event
    int IDprocUP;                   // ID of the process for this event
    int IDpartUP[MAXpart];          // particle ID according to Particle Data Group convention
    int statusUP[MAXpart];          // status code:
                                    //   -1 Incoming particle
                                    //   +1 Outgoing final state particle
                                    //   -2 Intermediate space-like propagator defining an 
                                    //   $x$ and $Q^2$ which should be preserved
                                    //   +2 Intermediate resonance, Mass should be preserved
                                    //   +3 Intermediate resonance, for documentation only
                                    //   -9 Incoming beam particles at time $t=-\infty$

    int motherUP[2][MAXpart];       // index of first and last mother
    int colorflowUP[2][MAXpart];    // color flow
                                    //   ICOLUP(1,I) - integer tag for the color flow 
                                    //   line passing through the color of the particle
                                    //   ICOLUP(2,I) - integer tag for the color flow 
                                    //   line passing through the anti-color of the particle

    double XweightUP;               // event weight
    double QscaleUP;                // scale of the event in GeV, as used for calculation of  PDFs
    double QEDalphaUP;              // the QED coupling used for this event
    double QCDalphaUP;              // the QCD coupling used for this event
    double momentumUP[5][MAXpart];  // lab frame momentum $(P_x, P_y, P_z, E, M)$ of particle in GeV
    double timelifeUP[MAXpart];     // invariant lifetime $c\tau$ (distance from 
                                    // production to decay) in mm
    double spinUP[MAXpart];         // cosine of the angle between the spin-vector of 
                                    // particle I and the 3-momentum of the decaying 
                                    // particle, specified in the lab frame
  }
*/

static double get_llprod (double * p, double * q) {
  return p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
}

double 
cutFunction (eventUP * ev)
{
/*
An example:
All leptons are selected. If any invariant mass (l1,l2) 
in the event is less than 10. GeV, the event is rejected

Structure of the eventUP structure is given above
*/
  int i, j;
  int nlep = 0;
  int nmass = 0;
  double selection = 1.;
  double leppx[16];
  double leppy[16];
  double leppz[16];
  double lepp0[16];
  double leppm[16];
  double llmass[16] = {0.};

  for (i = 0; i < ev->NpartUP; ++i) {
    int kf = abs (ev->IDpartUP[i]);
    if (11 == kf || 13 == kf || 15 == kf) {
      leppx[nlep] = ev->momentumUP[0][i];
      leppy[nlep] = ev->momentumUP[1][i];
      leppz[nlep] = ev->momentumUP[2][i];
      lepp0[nlep] = ev->momentumUP[3][i];
      leppm[nlep] = ev->momentumUP[4][i];
      ++nlep;
    }
  }

  for (i = 0; i < nlep; ++i) {
    double p1[4] = {lepp0[i], leppx[i], leppy[i], leppz[i]};
    double ms1 = leppm[i] * leppm[i];
    for (j = i + 1; j < nlep; ++j) {
      double p2[4] = {lepp0[j], leppx[j], leppy[j], leppz[j]};
      double ms2 = leppm[j] * leppm[j];
      llmass [nmass] = sqrt (2. * get_llprod (p1, p2) + ms1 + ms2);
      ++nmass;
    }
  }

  for (i = 0; i < nmass; ++i) {
    if (llmass[i] < 10.)
      selection = 0.;
  }

  return selection;
}
