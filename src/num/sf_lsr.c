/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Slava Ilyin 
* ------------------------------------------------------
*/
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/f_c.h"

#include "sf_lsr.h"

int 
p_lsr__ (char * p_name)
{
  if (!strcmp (p_name, "A"))
    return 1;
  else
    return 0;
}

void 
info_lsr__ (int i, Str_fun_Info * info)
{
  strcpy (info->pdf_name, "Laser photons");
  strcpy (info->prt_name, "gamma");
  info->prt_mass = 0.0;
  strcpy (info->version, "");
  info->N_extra_commands = 0;
}

void 
n_lsr__ (int i, char * beam, char * pdf)
{
  strcpy (beam, "photon");
  strcpy (pdf, "Laser photon");
}

int 
r_lsr__ (int i, char *name)
{
  if (strcmp (name, "Laser photon(photon)") != 0)
    return 0;
  else
    return 1;
}

int b_lsr__ (int i) {
  return 0;
}

int 
i_lsr__ (int i, double *be, double *mass, char * p_name)
{
  *be = 1.;
  *mass = 0.;
  return 1;
}

int 
m_lsr__ (int i, char * p_name)
{
  return 99;
}

double 
c_lsr__ (int i, double x, double q)
{

  static int first = TRUE;
  static double x0 = 4.82;
  static double xmax;

  /* System generated locals */
  double d__1, d__2;

  /* Local variables */
  static double rnorma;

  xmax = x0 / (1 + x0);
  if (first)
    {
      first = FALSE;
/* Computing 2nd power */
      d__1 = x0;
/* Computing 2nd power */
      d__2 = x0 + 1.;
      rnorma = (1. - 4. / x0 - 8. / (d__1 * d__1)) * log (x0 + 1.) + .5 + 8.
	/ x0 - 1. / (d__2 * d__2 * 2.);
    }
  if (x > xmax)
    return 0.;
  else
    return (1. - x + 1. / (1. - x) * (1. - x * 4. / x0 * (1. - x / (
						 x0 * (1. - x))))) / rnorma;
}
