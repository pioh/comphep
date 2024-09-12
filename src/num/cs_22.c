/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* --------------------------------------------------
*/
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "service2/include/4_vector.h"
#include "chep_crt/include/chep_crt.h"
#include "plot/include/plot.h"
#include "out_ext.h"

#include "num_serv.h"
#include "simpson.h"
#include "const.h"
#include "subproc.h"
#include "mc_menu.h"
#include "err_code.h"
#include "cut.h"
#include "kinaux.h"
#include "cs_22.h"

static double totcoef;
static double cos1, cos2;
static double eps = 0.001;
static int recalc;
static char procname[STRSIZ];

static double pmass[4];
static char pname[4][6];
static double pRestIn, pRestOut;


static void 
findAngleRange (double *cos1ptr, double *cos2ptr)
{
  int ncut;
  char c13[3] =
  {1, 3, 0};
  char c24[3] =
  {2, 4, 0};
  char c23[3] =
  {2, 3, 0};
  char c14[3] =
  {1, 4, 0};

  *cos1ptr = -1.0;
  *cos2ptr = 1.0;

  for (ncut = 0; invcut_1[ncut].key; ncut++)
    if ('C' == invcut_1[ncut].key)
      {

	if (eqvect_ (c13, invcut_1[ncut].lvinvc) ||
	    eqvect_ (c24, invcut_1[ncut].lvinvc))

	  {
	    if (invcut_1[ncut].minon && *cos1ptr < invcut_1[ncut].cvmin)
	      *cos1ptr = invcut_1[ncut].cvmin;
	    if (invcut_1[ncut].maxon && *cos2ptr > invcut_1[ncut].cvmax)
	      *cos2ptr = invcut_1[ncut].cvmax;
	  }

	if (eqvect_ (c23, invcut_1[ncut].lvinvc) ||
	    eqvect_ (c14, invcut_1[ncut].lvinvc))

	  {
	    if (invcut_1[ncut].minon && *cos1ptr < -invcut_1[ncut].cvmax)
	      *cos1ptr = -invcut_1[ncut].cvmax;
	    if (invcut_1[ncut].maxon && *cos2ptr > -invcut_1[ncut].cvmin)
	      *cos2ptr = -invcut_1[ncut].cvmin;
	  }
      }
}

static void 
infotext (void)
{
  goto_xy (1, 4);
  scrcolor (Red, BGmain);
  print (" %-53s\n", "Energy       :");
  print (" Cos(p1,p3): min=                   max=            \n");
  print (" %-53s", "Cross Section: ");
  scrcolor (FGmain, BGmain);
}


static void 
writeinformation (void)
{
  double sqrt_S;
  vinf_ (0, NULL, &sqrt_S);
  scrcolor (FGmain, BGmain);
  goto_xy (18, 4);
  print ("%12f [GeV]    ", sqrt_S);	/*  Energy  */
  goto_xy (18, 5);
  print ("%8.6f", cos1);
  goto_xy (42, 5);
  print ("%8.6f", cos2);
}


static void 
calccoef (void)
{
  int i;
  double lambda12, lambda34, s_, ms, mdiff;
  double sqrt_S;
  vinf_ (0, NULL, &sqrt_S);
  s_ = sqrt_S;


  err_code = 0;

  for (i = 0; i < 4; i++)
    {
      pinf_ (proces_1.nsub, i + 1, NULL, pmass + i);
      if (err_code != 0)
	goto errorexit;
    }

  s_ = sqrt_S * sqrt_S;

  ms = pmass[0] + pmass[1];
  if (ms >= sqrt_S)
    goto errorexit;
  mdiff = pmass[0] - pmass[1];
  lambda12 = sqrt ((s_ - ms * ms) * (s_ - mdiff * mdiff));

  pRestIn = lambda12 / (2 * sqrt_S);

  ms = pmass[2] + pmass[3];
  if (ms >= sqrt_S)
    goto errorexit;
  mdiff = pmass[2] - pmass[3];
  lambda34 = sqrt ((s_ - ms * ms) * (s_ - mdiff * mdiff));

  pRestOut = lambda34 / (2 * sqrt_S);

  totcoef = 3.8937966E8 * lambda34 / (32.0 * M_PI * lambda12 * s_);

  for (i = 0; i < 16; i++)
    pvect[i] = 0;

  pvect[3] = pRestIn;
  pvect[7] = -pRestIn;
  pvect[0] = sqrt (pRestIn * pRestIn + pmass[0] * pmass[0]);
  pvect[4] = sqrt (pRestIn * pRestIn + pmass[1] * pmass[1]);
  pvect[8] = sqrt (pRestOut * pRestOut + pmass[2] * pmass[2]);
  pvect[12] = sqrt (pRestOut * pRestOut + pmass[3] * pmass[3]);

  err_code = 0;
  return;

errorexit:
  if (err_code == 0)
    err_code = 4;
}


static void 
calcscalars (double cos_f)
{
  double sin_f = sqrt (fabs ((1 - cos_f) * (1 + cos_f)));
  pvect[11] = pRestOut * cos_f;
  pvect[15] = -pvect[11];
  pvect[10] = pRestOut * sin_f;
  pvect[14] = -pvect[10];
}


static double 
cross_section (double x)
{
  double r;
  calcscalars (x);
  r = sqme_ (proces_1.nsub, pvect, &err_code) * calcCutFactor ();
  if (err_code != 0)
    return 0;
  return r * totcoef;
}


static int 
fillseq (int n, double *f)
{

  int i;
  double step;

  err_code = 0;
  step = (cos2 - cos1) / (n - 1);
  for (i = 0; i < n; i++)
    {
      f[i] = cross_section (cos1 + i * step);
      if (err_code != 0)
	return FALSE;
    }
  return TRUE;
}


static void 
drawgraph (void)
{
  int n = 101;
  double f[202];


  do
    {
      if (correctInt (56, 8, "Number of points=", &n, 1))
	{
	  if (n < 3)
	    warnanykey (56, 8, "Too few points");
	  if (n > 201)
	    warnanykey (56, 8, "Too many points");
	}
      else
	return;
    }
  while (n < 3 || n > 201);


  if (!fillseq (n, f))
    {
      warnanykey (10, 10, "Error in calculation");
      return;
    }

  plot_histo (cos1, cos2, n, f, NULL, procname, "cos(p1,p3)", "Diff. cross section [pb]");

}


static double 
totcs (void)
{
  double int_val = 0.0;
  calccoef ();
  if (err_code == 0)
    int_val = simpson (cross_section, cos1, cos2, eps);
  return int_val;
}


static double 
acs (void)
{
  double int_val, int_val_a;

  calccoef ();
  if (err_code == 0)
    {
      int_val = simpson (cross_section, 0.0, cos2, eps);
      int_val_a = simpson (cross_section, cos1, 0.0, eps);
      return (int_val - int_val_a) / (int_val + int_val_a);
    }
  else
    return 0.0;
}


static void 
total_cs (void)
{
  double totcs;

  goto_xy (18, 6);
  scrcolor (FGmain, BGmain);
  print ("?                            ");
  goto_xy (18, 6);
  refresh_scr ();
  calccoef ();
  if (err_code)
    print ("incorrect");
  else
    {
      totcs = simpson (cross_section, cos1, cos2, eps);
      if (!err_code)
	print ("%-G [pb]", totcs);
    }
}


int 
cs_numcalc (void)
{
  int k, l;
  void *pscr0 = NULL;
  void *pscr = NULL;
  void *pscr2 = NULL;


  get_text (1, 3, 60, 11, &pscr0);

  for (k = 1; k <= 4; k++)
    {
      pinf_ (proces_1.nsub, k, pname[k - 1], NULL);
      trim (pname[k - 1]);
    }
  sprintf (procname, "%s,%s ->%s,%s", pname[0], pname[1], pname[2], pname[3]);

  findAngleRange (&cos1, &cos2);

  infotext ();
  writeinformation ();
  k = 1;
  l = 1;

  recalc = 1;
  do
    {
      if (recalc)
	{
	  total_cs ();
	  recalc = 0;
	  if (err_code)
	    errormessage ();
	}
      menu1 (54, 4, "", "\030"
	     " Set precision          "
	     " Angular dependence     "
	     " Parameter dependence   ", "n_22_*", &pscr, &k);

      switch (k)
	{
	case 0:
	  break;
	case 1:
	  do
	    {			/* Precision */
	      recalc = correctDouble (1, 23, " Enter precision : ", &eps, 1);
	      if (eps < 1.E-10 || eps > 0.0011)
		warnanykey (10, 12, "Range check error");
	    }
	  while (!(eps >= 1.E-10 && eps <= 0.03));
	  break;
	case 2:
	  if (err_code)
	    errormessage ();
	  else
	    drawgraph ();
	  break;
	case 3:
	  if (err_code)
	    errormessage ();
	  else
	    do
	      {
		l = 1;
		menu1 (54, 7, "", "\030"
		       "   Total Cross Section  "
		       "   Asymmetry            ", "", &pscr2, &l);
		if (l == 1)
		  paramdependence (totcs, procname, "Cross Section [pb]");
		if (l == 2)
		  paramdependence (acs, procname, "Asymmetry");
	      }
	    while (l != 0);
	  break;
	}			/*  switch  */
      if (k > 0)
	writeinformation ();
    }
  while (k != 0);
  put_text (&pscr0);
  return 0;
}
