/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/f_c.h"
#include "chep_crt/include/crt_util.h"

#include "out_ext.h"

#include "subproc.h"
#include "simpson.h"
#include "strfun.h"
#include "tools.h"
#include "sf_isr.h"

#ifdef _WIN32
#include "mathtools.h"
#endif

#define NPOINTS 100
#define   EM     5.1099906e-4
#define   EGAM   0.5772156649	/* Euler constant */
#define  ALPHA   0.0072973530796448189
#define  mmToGeV 5.067E12

#define EPS 1.e-6
#define B2 (1./3.)
#define R3 (1./3.)

static double scale = 100, xy_nm = 560, z_mm = 0.4, qTot = 2.E10;
static int bOn = 0;

static double beta, coeff, b_ncl, b_ips;

static double xi[NPOINTS], yi[NPOINTS];

/* E.A.Kuraev,V.S.Fadin:Sov.J.Nucl.Phys.41(1985)466
   S.Jadach,B.F.L.Ward:Comp.Phys.Commun.56(1990)351 */


static double 
b_h (double etax)
{
  int n;
  double s0 = 0., ds, g[3];

  g[0] = 0.37328217390739632;
  if (etax <= 0.)
    return g[0] * gammai_ (2, b_ncl);
  g[1] = pow (etax, R3) / 1.354117939426404 / 2.;
  g[2] = pow (etax, 2. * R3) / 6.;

  for (n = 1; n < 4; n++)
    {
      s0 += g[n - 1] * gammai_ (n + 1, b_ncl);
    }

  do
    {
      int n3 = (n - 1) % 3;
      g[n3] *= 3 * etax / (n * (n - 1) * (n - 2) * (n - 3));
      ds = g[n3] * gammai_ (++n, b_ncl);
      s0 += ds;
    }
  while (ds > s0 * EPS / 100);

  return s0;
}


static double 
cfbeam (double x)
{
  if (x <= 0.)
    return 0.0;
  else
    {
      double k = 2. / (b_ips * 3.);
      double etax = k * (1 / x - 1);
      if (etax > 50.)
	return 0;
      return pow (k / x, R3) * exp (-etax) / (x) * b_h (etax);
    }
  return 0.0;
}


static double 
cfbeamLog (double y)
{
  return cfbeam (exp (-y)) * pow (divy_ (y), -2. * R3);
}

static double 
cfisr (double x)
{
  return coeff * (x * x + 1 - beta * (log (x) * (3. * x * x + 1.) / 2. + (1. - x) * (1. - x)) / 2.) / 2;
}

static double 
cfisrLog (double y)
{
  return cfisr (exp (-y)) * pow (divy_ (y), beta - 1);
}

int 
p_isr__ (char * p_name)
{
  if (!strcmp (p_name, "e1") ||
      !strcmp (p_name, "E1") ||
      !strcmp (p_name, "e") ||
      !strcmp (p_name, "E"))
    return 1;
  else
    return 0;
}


void 
info_isr__ (int i, Str_fun_Info * info)
{
  double prt_mass;
  char prt_name[20];

  strcpy (info->pdf_name, "ISR");
  pinf_ (proces_1.nsub, i + 1, prt_name, &prt_mass);
  if (!strcmp (prt_name, "e"))
    strcpy (info->prt_name, "electron");
  else
    {
      if (!strcmp (prt_name, "E"))
	strcpy (info->prt_name, "positron");
      else
	{
	  fprintf (stderr, "error: unknown beam particle for electron ISR!");
	  exit (1);
	}
    }
  info->prt_mass = 5.11E-4;
  strcpy (info->version, "");

  if (bOn)
    {
      info->N_extra_commands = 7;
      strcpy (info->extra_commands[1], "BeamStralung=ON");
      snprintf (info->extra_commands[2], 26, "Bunch_x+y_sizes(nm)=%12.5E", xy_nm);
      snprintf (info->extra_commands[3], 30, "Bunch_lenght(mm)=%12.5E", z_mm);
      snprintf (info->extra_commands[4], 33, "Number_of_particles=%12.5E", qTot);
      snprintf (info->extra_commands[5], 18, "N_cl=%12.5E", b_ncl);
      snprintf (info->extra_commands[6], 21, "Upsilon=%12.5E", b_ips);
    }
  else
    {
      info->N_extra_commands = 2;
      strcpy (info->extra_commands[1], "BeamStralung=OFF");
    }
  snprintf (info->extra_commands[0], 22, "ISRscale=%12.5E", scale);
}


void 
n_isr__ (int i, char * beam, char * pdf)
{
  if (bOn) {
    sprintf (pdf, "ISR(%.0f Beamstr.: %.0f,%.2f,%.1E )", scale, xy_nm, z_mm, qTot);
  } else {
    sprintf (pdf, "ISR(%.0f Beamstr.: OFF)", scale);
  }
  sprintf (beam, "electron");
}

static double 
f_test (double x)
{
  return c_isr__ (1, 1 - pow (x, 1 / beta), 0.);
}


static void 
calc_params (void)
{
  double sqrt_S;

  vinf_ (0, NULL, &sqrt_S);
  beta = ALPHA * (2 * log (scale / EM) - 1) / M_PI;
  coeff = exp (beta * (0.75 - EGAM) - lgamma (1 + beta));
  if (bOn)
    {
      b_ncl = 25 * ALPHA * ALPHA * qTot / (12 * EM * (xy_nm * 1.E-6) * mmToGeV);
      b_ips = 5 * ALPHA * qTot * sqrt_S / (12 * EM * EM * EM * z_mm * (xy_nm * 1.E-6) * mmToGeV * mmToGeV);
    }
}


int 
i_isr__ (int ii, double *be, double *mass, char * p_name)
{
  int i;
  static int bOn_old = -1;
  static double beta_old = 0, coeff_old = 0, b_ncl_old = 0, b_ips_old = 0;

  *mass = 5.11E-4;
  calc_params ();

  if (beta == beta_old && coeff == coeff_old && bOn == bOn_old)
    {
      if (!bOn)
	{
	  *be = beta;
	  return 1;
	}
      if (b_ncl == b_ncl_old && b_ips == b_ips_old)
	{
	  *be = beta;
	  return 1;
	}
    }

  for (i = 0; i < NPOINTS; ++i)
    {
      double x = (double) (i) / NPOINTS;
      xi[i] = x;
      x = 1 - x * x * x;
      yi[i] = cfisr (x);
      if (bOn)
	{
	  double lx = -log (x);
	  yi[i] = (yi[i] * (1 - exp (-(b_ncl))) + pow (1 - x, B2) * pow (divy_ (lx), 1 - beta - B2)
		* convol_ (cfisrLog, cfbeamLog, beta, B2, lx, EPS)) / b_ncl;
	}
    }
  printf ("ISR integral %f\n", simpson (f_test, 0., 1., 1.E-8));

  beta_old = beta;
  coeff_old = coeff;
  bOn_old = bOn;
  b_ncl_old = b_ncl;
  b_ips_old = b_ips;

  *be = beta;
  return 1;
}


int 
r_isr__ (int i, char *name)
{
  double z0, z1, z2, z3;

  if (sscanf (name, "ISR(%lf Beamstr.: %lf,%lf,%lf)", &z0, &z1, &z2, &z3) == 4)
    {
      bOn = 1;
      scale = z0;
      xy_nm = z1;
      z_mm = z2;
      qTot = z3;
    }
  else if (sscanf (name, "ISR(%lf Beamstr.: OFF)", &z0) == 1)
    {
      bOn = 0;
      scale = z0;
    }
  else
    return 0;
  return 1;
}

int b_isr__ (int i) {
  return 0;
}


double 
c_isr__ (int i, double x, double q)
{
  x = pow (1 - x, R3);
  return dinter_ (x, NPOINTS, xi, yi);
}


int 
m_isr__ (int i, char * p_name)
{
  void *pscr = NULL;
  int mode = 1;
/*  char xff[128]; */

  for (;;)
    {
      char strmen[] = "\40"
      " ISR scale (GeV)     = XXX      "
      " Beamstralung          ON       "
      " Bunch x+y sizes (nm)= YYY      "
      " Bunch lenght (mm)   = ZZZ      "
      " Number of particles = NNN      "
      "          *     N_cl = NCL      "
      "          *  Upsilon = UPS      ";

      improveStr (strmen, "XXX", "%.1f", scale);
      calc_params ();
      if (bOn)
	{
	  improveStr (strmen, "YYY", "%.0f", xy_nm);
	  improveStr (strmen, "ZZZ", "%.2f", z_mm);
	  improveStr (strmen, "NNN", "%.1e", qTot);
	  improveStr (strmen, "NCL", "%.2f", b_ncl);
	  improveStr (strmen, "UPS", "%.2f", b_ips);
	}
      else
	{
	  improveStr (strmen, "ON", "%3.3s", "OFF");
	  strmen[2 * strmen[0] + 1] = 0;
	}
/*      strFunName(i,xff);
   k=strcmp(xff,"OFF");
   if (!k)
   improveStr(strmen,"XFF","%3.3s","ON ");
   else
   improveStr(strmen,"XFF","%3.3s","OFF");
 */
      menu1 (46, 10, "", strmen, "n_sf_isr", &pscr, &mode);

      switch (mode)
	{
	case 0:
	  return 0;
	case 1:
	  correctDouble (52, 16, "Enter new value ", &scale, 1);
	  break;
	case 2:
	  bOn = !bOn;
	  break;
	case 3:
	  correctDouble (52, 16, "Enter new value ", &xy_nm, 1);
	  break;
	case 4:
	  correctDouble (52, 16, "Enter new value ", &z_mm, 1);
	  break;
	case 5:
	  correctDouble (52, 16, "Enter new value ", &qTot, 1);
	  break;
	case 6:
	case 7:
	  messanykey (10, 10, "This parameter is a function of\n above ones and Sqrt(S)");
	}
    }
  return 99;
}
