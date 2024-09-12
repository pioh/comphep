/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Slava Ilyin 
* ----------------------------------------------------
*/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "service2/include/chep_limits.h"
#include "service2/include/read_func.h"
#include "service2/include/parser.h"
#include "service2/include/syst.h"
#include "chep_crt/include/crt_util.h"
#include "plot/include/plot.h"

#include "strfun_par.h"
#include "strfun.h"
#ifdef LHAPDF
  #include "clhapdf.h"
  #include "lhapdf.h"
  #include "sf_lhapdf.h"
#else
  #include "pdf.h"
  #include "sf_pdf.h"
#endif

#include "alphas2.h"

#define  XMB 4.3
#define  XMC 1.3
#define  XMTOP 175.

static double b0[7];
static double b1[7];
static double b2[7];
static double qcdL6 = 0.1185;

static int re_alpha = 1;


void recalc_alphas (void)
{
  re_alpha = 1;
}

void setLambda6 (double lam)
{
  qcdL6 = lam;
}

/**************************************************
* Transformation of Lambda_QCD from NF to NF-1    *
*    (formula from PDG-94)                        *
**************************************************/
static double 
tonf1 (double qcdl, double xmq, int nf)
{
  double rl, b10n, b10n1, a;

  rl = log (xmq / qcdl);

  b10n = b1[nf] / b0[nf];
  b10n1 = b1[nf - 1] / b0[nf - 1];

  a = qcdl * exp (-1 / b0[nf - 1] *
                  (

                    (b0[nf] - b0[nf - 1]) * rl
                    + (b10n - b10n1) * log (rl * 2.)
                    - b10n1 * log (b0[nf] / b0[nf - 1]) +
                    (
                      b10n * (b10n - b10n1) * log (rl * 2.)

                      + b10n * b10n - b10n1 * b10n1
                      - 0.5 * b2[nf] / b0[nf] + 0.5 * b2[nf - 1] / b0[nf - 1]
                      - 7. / 18.
                    ) / (b0[nf] * rl)
                  )
    );
/* fprintf(stderr,"\na=%f, b10n1=%f, b10n=%f, rl=%f\n",a,b10n1,b10n,rl); */
  return a;
}

double 
alpha_2 (double dscale)
{
  double alphas = 0.;
  static double qcdlf[7] = {0., 0., 0., 0., 0., 0., 0.};

  if (re_alpha) {
    int i;
    for (i = 0; i < 7; ++i) {
      b0[i] = 11. - (2. / 3.) * i;
      b1[i] = 51. - (19. / 3.) * i;
      b2[i] = 2857. - (5033. / 9.) * i + (325. / 27.) * i * i;
    }
    qcdlf[6] = qcdL6;
    qcdlf[5] = tonf1 (qcdlf[6], XMTOP, 6);
    qcdlf[4] = tonf1 (qcdlf[5], XMB, 5);
    qcdlf[3] = tonf1 (qcdlf[4], XMC, 4);

    re_alpha = 0;
  }
/*fprintf(stderr,"\n3=%f  4=%f, 5=%f, 6=%f\n",qcdlf[3],qcdlf[4],qcdlf[5],qcdlf[6]); */

  if (get_alphaMode ()) {
#ifdef LHAPDF
    alphas = alpha_lhapdf (dscale);
#else
    alphas = alpha_pdf (dscale);
#endif
  } else {
    int nf = 3;
    if (dscale > XMC)
      nf = 4;
    if (dscale > XMTOP)
      nf = 6;
    if (dscale > XMB)
      nf = 5;
    {
      double rl = 2 * log (dscale / qcdlf[nf]);
      double d1 = 2 * b1[nf] / (b0[nf] * b0[nf] * rl);
      double d2 = log (rl) - 0.5;
      alphas = 4 * M_PI / (b0[nf] * rl) * (1 - 2 * b1[nf] * log (rl) / (b0[nf] * b0[nf] * rl)
             + d1 * d1 * (d2 * d2 + b2[nf] * b0[nf] / (8 * b1[nf] * b1[nf]) - 1.25));
    }
  }

  return alphas;
}

double QCDLambda (void) {
  double lambda = qcdL6;
  if (get_alphaMode()) {
#ifdef LHAPDF
    lambda = lhapdf_QCDLambda ();
#else
    lambda = pdf_QCDLambda ();
    
#endif
  }
  return lambda;
}

int QCDOrder (void) {
#ifdef LHAPDF
  return lhapdfqcdorder ();
#else
  return 2;
#endif
}

int Nflavour (void) {
  double nfl = 6;
#ifdef LHAPDF
  if (get_alphaMode ()) {
    nfl = 5;
  }
#endif

  return nfl;
}

static double realPi (double r) {
  /* use assymptotic formula */
  if (fabs (r) < 1e-3) {
    return - 1. - 2./3. - log (r);
  }
  /* return zero for large values */
  else if (fabs(r) > 1.e6) {
    return 0.;
  }
  else if(4. * r > 1.) {
    double beta = sqrt (4. * r - 1.);
    double a = acos (1. - 1. / (2. * r));
    return 1./3. -(1. + 2. * r) * (2. - beta * a);
  }
  else {
    double beta = sqrt (1. - 4. * r);
    double a = fabs ((beta - 1.) / (beta + 1.));
    return 1./3. - (1. + 2. * r) * (2. + beta * log (a));
  }
}

double 
alpha_em (double dscale) {
  double alem = 1/137.0;
  double aempi = alem / (3. * M_PI);

  double m1 = 0.000511;
  double m2 = 0.10566;
  double m3 = 1.77699;
  double m4 = 174.3;

  double x = dscale * dscale;

  if (x < 1e-6) {
    return alem;
  } else {
  /* leptonic and t-quark contributions */
    double a1 = 0.0    , b1 = 0.00835, c1 = 1.000;
    double a2 = 0.0    , b2 = 0.00238, c2 = 3.927;
    double a3 = 0.00165, b3 = 0.00299, c3 = 1.000;
    double a4 = 0.00221, b4 = 0.00293, c4 = 1.000;

    double k_e   = realPi (m1 / dscale);
    double k_m   = realPi (m2 / dscale);
    double k_tau = realPi (m3 / dscale);
    double k_top = realPi (m4 / dscale);
    double repigg = aempi * (k_e + k_m + k_tau + k_top);

  /* light quark contribution */
    if (x < 9e-2)      repigg += a1 + b1 * log (1. + c1 * x);
    else if (x < 9.)   repigg += a2 + b2 * log (1. + c2 * x);
    else if (x < 1.e4) repigg += a3 + b3 * log (1. + c3 * x);
    else               repigg += a4 + b4 * log (1. + c4 * x);
    return alem / (1. - repigg);
  }
}
