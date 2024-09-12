/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov, Slava Ilyin 
*------------------------------------------------------
*/
#include <math.h>
#include "4_vector.h"

double pvect[400] = {0.};

void 
lvtonv (char *lv, int nin, int nv)
{
  int i, n;

  vnull4 (nv);
  for (i = 0; (n = lv[i]); i++)
    if (n > nin)
      vsum4 (nv, n, nv, 1);
    else
      vsum4 (nv, n, nv, -1);
}


void
lorenc (double *mom, double *p_frame, double *p_trans)
{
  double v_frame[4];
  double mod_v_frame;
  double G;
  double coef;
  double p_v;

  v_frame[1] = p_frame[1] / p_frame[0];
  v_frame[2] = p_frame[2] / p_frame[0];
  v_frame[3] = p_frame[3] / p_frame[0];
  mod_v_frame = sqrt (v_frame[1] * v_frame[1] + v_frame[2] * v_frame[2] + v_frame[3] * v_frame[3]);
  p_v = mom[1] * v_frame[1] + mom[2] * v_frame[2] + mom[3] * v_frame[3];
  G = 1. / sqrt (1. - mod_v_frame * mod_v_frame);
  coef = G * ((G * p_v) / (G + 1.) - mom[0]);
  p_trans[1] = mom[1] + v_frame[1] * coef;
  p_trans[2] = mom[2] + v_frame[2] * coef;
  p_trans[3] = mom[3] + v_frame[3] * coef;
  p_trans[0] = G * (mom[0] - p_v);
}


void
new_lorenc (double *mom, double *p_frame, double *newmom)
{
  int i;

  if (1.e-10 > fabs (p_frame[0])) {
    for (i = 0; i < 4; ++i)
      newmom[i] = mom[i];
  } else {
    double v_frame[4];
    double mod_v_frame;
    for (i = 1; i < 4; ++i)
     v_frame[i] = p_frame[i] / p_frame[0];

    mod_v_frame = sqrt (v_frame[1] * v_frame[1] + v_frame[2] * v_frame[2] + v_frame[3] * v_frame[3]);
    if (1.e-10 > mod_v_frame) {
      for (i = 0; i < 4; ++i)
        newmom[i] = mom[i];
    } else {
      double p_v = mom[1] * v_frame[1] + mom[2] * v_frame[2] + mom[3] * v_frame[3];
      double G_inversed = sqrt (1. - mod_v_frame * mod_v_frame);
      double coef = (p_v / (G_inversed + 1.) - mom[0]) / G_inversed;

      newmom[0] = (mom[0] - p_v) / G_inversed;
      for (i = 1; i < 4; ++i)
        newmom[i] = mom[i] + v_frame[i] * coef;
    }
  }

  return;
}

/* ******************************** */
/* Improved SQRT function        * */
/* ******************************** */
double 
vsqrt (double a)
{
  if (a < 0)
    return 0;
  else
    return sqrt (a);
}

/* ****************************************** */
/*    Scalar product of two 4-vectors:     * */
/* ****************************************** */
double 
vdot4 (int i, int j)
{
  i = 4 * i - 4;
  j = 4 * j - 4;
  return pvect[i] * pvect[j] - pvect[i + 1] * pvect[j + 1]
    - pvect[i + 2] * pvect[j + 2] - pvect[i + 3] * pvect[j + 3];
}

double 
mom4mod (int i)
{
  i = 4 * i - 4;
  return (pvect[i] - pvect[i + 1]) *(pvect[i] + pvect[i + 1]) - pvect[i + 2] * pvect[i + 2] - pvect[i + 3] * pvect[i + 3];
}

/* ******************************************************* */
/*       SUM or Difference of two 4-vectors:            *  */
/* ISG=1  sum   and  ISG=-1   difference                *  */
/*             P(I) + ISG*P(J)=>P(K)                    *  */
/* ******************************************************* */

void 
vsum4 (int i, int j, int k, int isg)
{
  int l;
  i = 4 * i - 4;
  j = 4 * j - 4;
  k = 4 * k - 4;
  if (isg == 1)
    {
      for (l = 0; l < 4; l++)
	pvect[k + l] = pvect[i + l] + pvect[j + l];
    }
  else
    {
      for (l = 0; l < 4; l++)
	pvect[k + l] = pvect[i + l] - pvect[j + l];
    }
}				/* vsum4_ */

/* ****************************************** */
/* NULLification of a 4-vector:  P(I) => 0 * */
/* ****************************************** */

void 
vnull4 (int i)
{
  int k;
  for (k = 4 * i - 4; k < 4 * i; k++)
    pvect[k] = 0;
}

void 
eps4 (int n1, int n2, int n3, int n4)
{
  double a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33, d1021,
    d1022, d1023, d1122, d1123, d1223;

  n1 *= 4;
  n2 *= 4;
  n3 *= 4;
  n4 *= 4;

  a10 = pvect[n1 - 4];
  a11 = pvect[n1 - 3];
  a12 = pvect[n1 - 2];
  a13 = pvect[n1 - 1];
  a20 = pvect[n2 - 4];
  a21 = pvect[n2 - 3];
  a22 = pvect[n2 - 2];
  a23 = pvect[n2 - 1];
  a30 = pvect[n3 - 4];
  a31 = pvect[n3 - 3];
  a32 = pvect[n3 - 2];
  a33 = pvect[n3 - 1];
/*                               A10  A20  A30  X0 */
/*                               A11  A21  A31  X1 */
/*                               A12  A22  A32  X2 */
/*                               A13  A23  A33  X3 */
  d1021 = a10 * a21 - a20 * a11;
  d1022 = a10 * a22 - a20 * a12;
  d1023 = a10 * a23 - a20 * a13;
  d1122 = a11 * a22 - a21 * a12;
  d1123 = a11 * a23 - a21 * a13;
  d1223 = a12 * a23 - a22 * a13;
  pvect[n4 - 4] = a31 * d1223 - a32 * d1123 + a33 * d1122;
  pvect[n4 - 3] = a30 * d1223 - a32 * d1023 + a33 * d1022;
  pvect[n4 - 2] = -(a30 * d1123 - a31 * d1023 + a33 * d1021);
  pvect[n4 - 1] = a30 * d1122 - a31 * d1022 + a32 * d1021;
}				/* eps4_ */

void 
pvFill (double mass, double *mom4, int pos)
{

  int i, i0 = 4 * (pos - 1);
  pvect[i0] = mass;
  pvect[i0] *= pvect[i0];
  mass = mass * mass;

  for (i = 1; i < 4; i++)
    {
      pvect[i0 + i] = mom4[i];
      pvect[i0] += pvect[i0 + i] * pvect[i0 + i];
    }
  pvect[i0] = sqrt (pvect[i0]);
}

void 
lorrot (double rapidity, int ntot)
{
  static double rapid___ = 0;
  static double sh = 0;
  static double ch = 1;
  double ee, pp;
  int nv;

  if (rapidity != rapid___)
    {
      rapid___ = rapidity;
      sh = sinh (rapidity);
      ch = sqrt (sh * sh + 1);
    }
  if (rapidity)
    for (nv = 4 * ntot - 4; nv >= 0; nv -= 4)
      {
	ee = pvect[nv];
	pp = pvect[nv + 3];
	pvect[nv] = ee * ch + pp * sh;
	pvect[nv + 3] = ee * sh + pp * ch;
      }
}
