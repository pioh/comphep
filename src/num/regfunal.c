/* 
* Copyright (C) 2001-2008, CompHEP Collaboration
* Copyright (C) 1997, Dmitry Kovalenko 
* ------------------------------------------------------
* $Id: regfunal.c,v 1.4 2003/06/25 18:46:41 kryukov Exp $
*
* $Log: regfunal.c,v $
* Revision 1.4  2003/06/25 18:46:41  kryukov
* *** empty log message ***
*
* Revision 1.3  2003/04/22 09:01:19  kryukov
* Apply indent to improve readability of files
*
* Revision 1.2  2003/04/20 17:55:02  kryukov
* Clean up warnings during make
*
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "service2/include/chep_limits.h"
#include "chep_crt/include/chep_crt.h"

#include "kinaux.h"
#include "regul.h"
#include "regfunal.h"

static double 
sing_ (int ityp, double pole, double width, int deg, double x)
{
  double dx = x - pole;

  if (ityp == 1)
    {
    L1:
      if (deg == 2)
	return 1 / (dx * dx);
      else
	return 1 / ABS (dx);
    }
  else
    {
      if (width == 0.)
	goto L1;
      return 1 / (dx * dx + width * width);
    }
}

static double 
singi_ (int ityp, double pole, double width,
	int deg, double x)
{
  double dx = x - pole;

  if (ityp == 1)
    {
      if (deg == 2)
	return -1 / dx;
      else
	{
	  if (dx > 0)
	    return log (dx);
	  else
	    return -log (-dx);
	}
    }
  else if (ityp == 2)
    {
      if (width == 0.)
	{
	  if (deg == 2)
	    return -1 / dx;
	  else
	    {
	      if (dx > 0)
		return log (dx);
	      else
		return -log (-dx);
	    }
	}
      return atan (dx / width) / width;
    }
  printf ("***** singi_: Invalid ityp=%d. ABORT\n", ityp);
  exit (99);
}

static double 
singi1_ (int ityp, double pole, double width,
	 int deg, double x, double xmin)
{
  if (ityp == 1)
    {
    L1:if (deg == 2)
	return pole - 1 / x;
      else
	{
	  if (pole < xmin)
	    return exp (x) + pole;
	  else
	    return pole - exp (-x);
	}
    }
  else
    {
      if (width == 0.)
	goto L1;
      return pole + width * tan (x * width);
    }
}				/* singi1_ */


static void 
regfun_0_ (int factOnly, int itype, int nsing,
	   sing_struct * singar, double xmin, double xmax, double xx,
	   double *xout, double *factor)
{
  double delt0 = 1e-13;
  double ai[101], bi[101];
  double rnorm_0;
  double delt;
  double tintxx;
  int tint0, i, ll, nn, ityp;

  if (nsing == 0)
    {
      if (!factOnly)
	{
	  *xout = xmin * (1 - xx) + xx * (xmax);
	  if (*xout > xmax)
	    *xout = xmax;
	  if (*xout < xmin)
	    *xout = xmin;
	}
      *factor = xmax - xmin;
      return;
    }

  /* Parameter adjustments */

  if (itype > 0)		/* constant is included */
    {
      ityp = itype;
      tint0 = 1;
      rnorm_0 = 1 / (xmax - xmin);
    }
  else
    {
      ityp = -itype;
      tint0 = 0;
      rnorm_0 = 0.;
    }

  for (i = 1; i <= nsing; ++i)
    {
      if (ityp == 1 || (ityp == 2 && singar[i - 1].width == 0.))
	{
	  if (singar[i - 1].pos >= xmin && singar[i - 1].pos <= xmax)
	    {
	      delt = delt0 * (ABS (xmin) + ABS (xmax));
	      if (singar[i - 1].pos - delt <= xmin)
		{
		  singar[i - 1].pos = xmin - delt;
		}
	      else if (singar[i - 1].pos + delt >= xmax)
		{
		  singar[i - 1].pos = xmax + delt;
		}
	      else
		{
		  printf ("%f < x < %f\n", xmin, xmax);
		  printf ("Type= %d  Bad reg= %d\n", itype, i);

		  for (ll = 1; ll <= nsing; ++ll)
		    {
		      printf ("%d %f %f %d\n", ll, singar[ll - 1].pos,
			      singar[ll - 1].width,
			      singar[ll - 1].power);
		    }
		  printf ("ERROR IN REGULARIZATION\n");
		  exit (0);
		}
	    }
	}
    }

  for (i = 1; i <= nsing; ++i)
    {
      ai[i] = singi_ (ityp, singar[i - 1].pos, singar[i - 1].width, singar[i - 1].power, xmin);
      bi[i] = singi_ (ityp, singar[i - 1].pos, singar[i - 1].width, singar[i - 1].power, xmax);
    }

  if (!factOnly)
    {
      tintxx = (tint0 + nsing) * xx;
      if (tint0 > tintxx)
	*xout = xmin + tintxx / rnorm_0;
      else
	{
	  tintxx -= tint0;
	  nn = 1 + tintxx;
	  tintxx -= (nn - 1);
	  *xout = singi1_ (ityp, singar[nn - 1].pos, singar[nn - 1].width,
	   singar[nn - 1].power, tintxx * (bi[nn] - ai[nn]) + ai[nn], xmin);
	}
      if (*xout > xmax)
	*xout = xmax;
      if (*xout < xmin)
	*xout = xmin;
    }
  *factor = rnorm_0;
  for (i = 1; i <= nsing; ++i)
    {
      *factor += sing_ (ityp, singar[i - 1].pos, singar[i - 1].width,
			singar[i - 1].power, *xout) / (bi[i] - ai[i]);
    }
  *factor = (tint0 + nsing) / *factor;
}				/* regfun_ */

void 
regfun_ (int itype, int nsing, sing_struct * singar, double xmin, double xmax,
	 double xx, double *xout, double *factor)
{
  regfun_0_ (0, itype, nsing, singar, xmin, xmax, xx, xout, factor);
}

void 
regfct_ (int itype, int nsing, sing_struct * singar, double xmin, double xmax,
	 double xout, double *factor)
{
  regfun_0_ (1, itype, nsing, singar, xmin, xmax, 0., &xout, factor);
}
