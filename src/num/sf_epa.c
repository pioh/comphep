/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Slava Ilyin 
* ------------------------------------------------------
*/
#include <math.h>
#include <stdlib.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/f_c.h"
#include "chep_crt/include/crt_util.h"

#include "out_ext.h"

#include "const.h"
#include "strfun.h"
#include "sf_epa.h"

static double xin[2] =
{5.11e-4, 5.11e-4};
static double charge[2] =
{-1., -1.};
static double q_max[2] =
{100., 100.};


/* ********************************************************* */
/*  Equivalent photon approximation structure function.   * */
/*     Improved Weizsaecker-Williams formula              * */
/*      C.F.Weizsaecker, Z.Phys. 88 (1934) 612            * */
/*      E.J.Williams,    Phys.Rev. 45 (1934) 729          * */
/*                                                        * */
/*   V.M.Budnev et al., Phys.Rep. 15C (1975) 181          * */
/* ********************************************************* */


int 
p_epa__ (char * p_name)
{
  if (strcmp (p_name, "A") == 0)
    return 1;
  else
    return 0;
}

void 
info_epa__ (int i, Str_fun_Info * info)
{
  strcpy (info->pdf_name, "WWA");
  strcpy (info->prt_name, "gamma");
  info->prt_mass = xin[i - 1];
  info->N_extra_commands = 2;
  snprintf (info->version, 1, "%s", "");
  snprintf (info->extra_commands[0], 20, "Ch=%12.5E", charge[i]);
  snprintf (info->extra_commands[1], 20, "Q=%12.5E", q_max[i]);
}


void 
n_epa__ (int i, char * beam, char * pdf)
{
  strcpy (beam, "photon");
  sprintf (pdf, "WWA (m=%.10G Ch=%.6G Q=%.8G)", xin[i - 1], charge[i - 1], q_max[i - 1]);
}


int 
r_epa__ (int i, char *name)
{
  double xin__, q_max__, charge__;

  if (3 != sscanf (name, "WWA (m=%lf%*[^=]%*c%lf%*[^=]%*c%lf", &xin__, &charge__, &q_max__))
    goto L10;
  if (xin__ <= 0.)
    goto L10;
  if (q_max__ <= 0.)
    goto L10;
  i--;
  xin[i] = xin__;
  q_max[i] = q_max__;
  charge[i] = charge__;
  return 1;
L10:return 0;
}

int b_epa__ (int i) {
  return 0;
}

int 
i_epa__ (int i, double *be, double *mass, char * p_name)
{
  *be = 1.;
  *mass = xin[i - 1];
  return 1;
}


int 
m_epa__ (int i, char * p_name)
{
  void *pscr = NULL;
  int mode;
  i--;
  for (;;)
    {
      char strmen[] = "\050"
      " Incoming particle mass = XXX           "
      " Incoming particle charge = YYY         "
      " |Q|max = ZZZ                           ";
/*  "123456789 123456789 123456789 123456789 "  */

      improveStr (strmen, "XXX", "%.10G GeV", xin[i]);
      improveStr (strmen, "YYY", "%.6G", charge[i]);
      improveStr (strmen, "ZZZ", "%.8G GeV", q_max[i]);
      menu1 (38, 10, "", strmen, "n_sf_epa", &pscr, &mode);

      switch (mode)
	{
	case 0:
	  return mode;
	case 1:
	  correctDouble (40, 16, "Enter new value ", xin + i, 1);
	  break;
	case 2:
	  correctDouble (40, 16, "Enter new value ", charge + i, 1);
	  break;
	case 3:
	  correctDouble (40, 16, "Enter new value ", q_max + i, 1);
	}
    }
  return 99;
}


double 
c_epa__ (int i, double x, double q)
{
  double alpha = .0072992701;
  double d__1, delt, f;

  i--;
  d__1 = xin[i] / q_max[i];
  delt = d__1 * d__1;
  f = alpha * charge[i] * charge[i] / (M_PI * 2) * (log ((1 - x) / (x * x * delt)) * ((1 - x) * (1 - x) + 1) / x
					- (1 - x - delt * (x * x)) * 2 / x);
  if (f < 0)
    f = 0;
  return f;
}
