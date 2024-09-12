/*
* Copyright (C) 2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "num/include/LesHouches.h"
#include "userFun.h"

int
userVariables (eventUP * ev, int nset, int * nvrs, double * vrs)
{
/*
    int NpartUP;
    int IDprocUP;
    int IDpartUP[MAXpart];
    int statusUP[MAXpart];
    int motherUP[2][MAXpart];
    int colorflowUP[2][MAXpart];
    double XweightUP;
    double QscaleUP;
    double QEDalphaUP;
    double QCDalphaUP;
    double momentumUP[5][MAXpart];
    double timelifeUP[MAXpart];
    double spinUP[MAXpart];
*/
  
  switch (nset) {
    case 0:
      {
        int i, j;
        *nvrs = 0;
        for (i = 0; i < ev->NpartUP; ++i) {
            double p_i[4] = {ev->momentumUP[0][i],ev->momentumUP[1][i],ev->momentumUP[2][i],ev->momentumUP[3][i]};
          for (j = i + 1; j < ev->NpartUP; ++j) {
            double p_j[4] = {ev->momentumUP[0][j],ev->momentumUP[1][j],ev->momentumUP[2][j],ev->momentumUP[3][j]};
            double val = p_i[0] * p_j[0] - p_i[1] * p_j[1] - p_i[2] * p_j[2] - p_i[3] * p_j[3];
            if (10e-8 < fabs (val)) {
              vrs[*nvrs] = val;
              ++(*nvrs);
            }
          }
        }
      }
      return 0;
    case 2:
      return 0;
  }

  fprintf (stdout, " ***** userVars: unknown variable set: %i\n", nset);
  return -1;
}
