/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "service2/include/chep_limits.h"
#include "service2/include/drandXX.h"
#include "service2/include/4_vector.h"
#include "chep_crt/include/crt_util.h"

#include "vegas.h"

int simplexOn = 1;

static void
drand_arr (int dim, double *x)
{
  int i, j;
  int rest = dim;
  unsigned int umax = UINT_MAX;
  unsigned int randpos = umax * drandXX ();

  for (i = 0; i < dim; i++)
    {
      x[i] = -1.;
    }

  for (i = 0; i < dim; i++)
    {
      int pos = randpos % rest;
      for (j = 0;; ++j)
        {
          if (x[j] < 0)
            {
              if (pos)
                {
                  pos--;
                }
              else
                {
                  x[j] = drandXX ();
                  break;
                }
            }
        }
      randpos /= rest;
      umax /= rest;
      rest--;
      if (umax < rest)
        {
          umax = UINT_MAX;
          randpos = umax * drandXX ();
        }
    }
}


static int Kg[MAX_DIM];
static int Ng[MAX_DIM];
static int Ndim;
static double (*f_) (double *, double);
static int Ndmx;
static double *Xgrid;

#define XG(j,i) Xgrid[(i)+(j)*(Ndmx)]

long
generateVegasCubs (vegasGrid * vegPtr, long nCubs1)
{
  int i;
  double nCubs_ = nCubs1;
  long nCubs = 1;

  Ndim = vegPtr->ndim;
  Ndmx = vegPtr->ndmx + 1;
  Xgrid = vegPtr->x_grid;

  for (i = 0; i < Ndim; i++)
    {
      double nd = Ndim - i;
      Kg[i] = 0;
      Ng[i] = pow (nCubs_, 1. / nd);
      nCubs_ /= Ng[i];
      nCubs *= Ng[i];
    }
  return nCubs;
}


static double
Local2Global (double XLOC[], double XGLOB[], int GRID_LOC[])
{
  int j;
  double JACOB = 1.;

  for (j = 0; j < Ndim; ++j)
    {
      double xlj = (Kg[j] + XLOC[j]) / Ng[j];
      int n = (int) (xlj * (Ndmx - 1));
      double xn2 = XG (j, n + 1);
      double xn1 = 0;
      if (n)
        {
          xn1 = XG (j, n);
        }
      XGLOB[j] = xn1 + (xn2 - xn1) * (xlj * (Ndmx - 1) - n);
      JACOB *= (xn2 - xn1) * (Ndmx - 1);
      if (GRID_LOC) {
        GRID_LOC[j] = n;
      }
    }
  return JACOB;
}

vegasGrid *
vegas_init (int dim, int nd)
{
  vegasGrid *vegPtr;

  if ((dim > MAX_DIM) || (nd > MAX_NDMX))
    {
      return NULL;
    }
  vegPtr = (vegasGrid *) malloc (sizeof (vegasGrid));
  if (vegPtr)
    {
      vegPtr->ndmx = nd;
      vegPtr->ndim = dim;
      vegPtr->x_grid = malloc (dim * (nd + 1) * sizeof (double));
      vegPtr->c_grid = malloc (dim * (nd + 1) * sizeof (double));

      Xgrid = vegPtr->x_grid;
      Ndmx = vegPtr->ndmx + 1;
      if (vegPtr->x_grid && vegPtr->c_grid)
        {
          int i, j;
          double *x_ = vegPtr->x_grid;
          double *c_ = vegPtr->c_grid;
          for (j = 0; j < dim; j++)
            {
              for (i = 0; i <= nd; i++, x_++, c_++)
                {
                  *x_ = i / (double) nd;
                  *c_ = 1. / nd;
                }
            }
        }
      else
        {
          if (vegPtr->x_grid)
            free (vegPtr->x_grid);
          if (vegPtr->c_grid)
            free (vegPtr->c_grid);
          free (vegPtr);
          vegPtr = NULL;
        }
    }
  return vegPtr;
}


void
vegas_finish (vegasGrid * vegPtr)
{
  if (vegPtr)
    {
      free (vegPtr->x_grid);
      free (vegPtr->c_grid);
      free (vegPtr);
    }
}

#define DD(j,i) d[(i)+(j)*Ndmx]

void static
refineGRID (int dim, double *d)
{
  int i, j, k;
  int ndm = Ndmx - 2;
  double alph = 1.5;     /* rate of grid improvement  */

  double r[MAX_NDMX];
  double dt[MAX_DIM];
  double xin[MAX_NDMX];
  double xn, xo, dr;

  for (i = 0; i < dim; ++i)
    {
      xo = DD (i, 0);
      xn = DD (i, 1);
      DD (i, 0) = (xo + xn) / 2.;
      dt[j] = DD (i, 0);
      for (j = 1; j < ndm; ++j)
        {
          DD (i, j) = xo + xn;
          xo = xn;
          xn = DD (i, j + 1);
          DD (i, k) = (DD (i, j) + xn) / 3.;
          dt[i] += DD (i, j);
        }
      DD (i, ndm) = (xn + xo) / 2.;
      dt[i] += DD (i, ndm);
    }

  for (j = 0; j < dim; ++j)
    {
      double rc = 0;
      for (i = 0; i <= ndm; ++i)
        {
          r[i] = 0;
          if (DD (j, i) > 0)
            {
              double xoln = log (dt[j] / DD (j, i));
              if (xoln <= 70.f)
                {
                  r[i] = pow ((1 - exp (-xoln)) / xoln, alph);
                }
              else
                {
                  r[i] = pow (1 / xoln, alph);
                }
            }
          rc += r[i];
        }

      rc /= (Ndmx - 1);
      if (rc)
        {
          for (i = 0, k = 0, xn = 0, dr = 0; i < ndm;)
            {
              do
                {
                  dr += r[k];
                  xo = xn;
                  xn = XG (j, k + 1);
                  k++;
                }
              while (rc > dr);

              do
                {
                  dr -= rc;
                  xin[i] = xn - (xn - xo) * dr / r[k - 1];
                  i++;
                }
              while (rc <= dr);
            }
          for (i = 0; i < ndm; ++i)
            {
              XG (j, i + 1) = xin[i];
            }
          XG (j, ndm + 1) = 1;
          XG (j, 0) = 0;
        }
    }
}

/*
                    *  VEGAS  *
The routine performs NDIM-dimensional Monte-Carlo integration
based on a code by G.P. Lepage, SEPT 1976/(REV)AUG 1979 
Algorithm described in J COMP PHYS 27,192 (1978) 
 */
int
vegas_int (vegasGrid * vegPtr, long ncalls, int nsubpr, int niters, int curitr, double (*fxn) (double *, double),
           double *ti, double *tsi)
{
  int i, j;
  int stop = 0;
  int dim = vegPtr->ndim;
  int fdim = dim * (vegPtr->ndmx + 1);
  long cCub = 0;
  long nCubs = generateVegasCubs (vegPtr, ncalls / 2);
  int npg = ncalls / nCubs;     /* number of calls per Cube! */
  double alph = 1.5;     /* rate of grid improvement  */
  double nCubs2 = (double)nCubs * nCubs;
  double *d = malloc (fdim * sizeof (double));
  FILE * percfile = fopen (".procent", "w");

  *ti = 0.;
  *tsi = 0.;
  for (i = 0; i < fdim; ++i)
    {
      d[i] = 0.;
    }

/* Main interation loop */
  while (cCub < nCubs && !stop)
    {
      int ia[MAX_DIM];
      double fb = 0.;
      double f2b = 0.;
      double xglobal[MAX_DIM];
      double xlocal[MAX_DIM];

      for (i = 0; i < npg; i++)
        {
          double f, f2;
          drand_arr (dim, xlocal);
          f = Local2Global (xlocal, xglobal, ia);
          f *= (*fxn) (xglobal, f);
          fb += f;
          f2 = f * f;
          f2b += f2;
          for (j = 0; j < dim; ++j)
            {
              DD (j, ia[j]) += f2;
            }
        }
      fb /= npg;

      f2b = (f2b / npg - fb * fb) / (npg - 1);
      *ti += fb / nCubs;
      *tsi += f2b / nCubs2;

      for (i = dim - 1; i >= 0; i--)
        {
          if (++Kg[i] < Ng[i])
            break;
          else
            Kg[i] = 0;
        }
      if (NULL != percfile) {
        char dummy[128];
        fseek (percfile, 0L, SEEK_SET);
        for (i = 1; i < nsubpr; ++i) fgets (dummy, 128, percfile);
        fprintf (percfile, "%i: %li %li %i %i\n", nsubpr, cCub, nCubs, niters, curitr);
        fflush (percfile);
      }
      stop = informline (cCub, nCubs);
      ++cCub;
    }
  if (NULL != percfile) {
    fclose (percfile);
  }
  *tsi = sqrt (*tsi);

/* if not stoped -> refine GRID */
  if (!stop)
    {
      double r[MAX_NDMX];
      double dt[MAX_DIM];
      double xin[MAX_NDMX];

      double xn, xo, dr;
      int k, ndm = Ndmx - 2;

      for (j = 0; j < dim; ++j)
        {
          xo = DD (j, 0);
          xn = DD (j, 1);
          DD (j, 0) = (xo + xn) / 2;
          dt[j] = DD (j, 0);
          for (i = 1; i < ndm; ++i)
            {
              DD (j, i) = xo + xn;
              xo = xn;
              xn = DD (j, i + 1);
              DD (j, i) = (DD (j, i) + xn) / 3;
              dt[j] += DD (j, i);
            }
          DD (j, ndm) = (xn + xo) / 2;
          dt[j] += DD (j, ndm);
        }
      for (j = 0; j < dim; ++j)
        {
          double rc = 0;
          for (i = 0; i <= ndm; ++i)
            {
              r[i] = 0;
              if (DD (j, i) > 0)
                {
                  double xoln = log (dt[j] / DD (j, i));
                  if (xoln <= 70.f)
                    r[i] = pow ((1 - exp (-xoln)) / xoln, alph);
                  else
                    r[i] = pow (1 / xoln, alph);
                }
              rc += r[i];
            }

          rc /= (Ndmx - 1);
          if (rc)
            {
              for (i = 0, k = 0, xn = 0, dr = 0; i < ndm;)
                {
                  do
                    {
                      dr += r[k];
                      xo = xn;
                      xn = XG (j, k + 1);
                      k++;
                    }
                  while (rc > dr);
                  do
                    {
                      dr -= rc;
                      xin[i] = xn - (xn - xo) * dr / r[k - 1];
                      i++;
                    }
                  while (rc <= dr);
                }
              for (i = 0; i < ndm; ++i)
                XG (j, i + 1) = xin[i];
              XG (j, ndm + 1) = 1;
              XG (j, 0) = 0;
            }
        }
    }

  free (d);
  return stop;
}
#undef DD

static int sstop = 1;
static double currentPos = 0.;

void
setStopSymb (int s)
{
  sstop = s;
}

int
getCurCub (void)
{
  return currentPos;
}

#define INCUB(x)   ((x)>0   ? (0.5+(x))/((x)+1.) : 0.5/(1.-(x)))
#define OUTCUB(x)  ((x)>0.5 ? (0.5-(x))/((x)-1.) : 1-0.5/x)

static double
amotry (double *p, double *y, int ndim,
        double (*f) (double *), int ilo, double fac)
{
  int i, j;
  double ytry;
  double fac1 = (1. - fac) / ndim;
  double *p_buff = p + (ndim + 1) * ndim;
  double *p_ilo = p + ilo * ndim;

  for (j = 0; j < ndim; j++)
    p_buff[j] = p_ilo[j] * fac;
  for (i = 0; i <= ndim; i++)
    if (i != ilo)
      {
        double *p_i = p + i * ndim;
        for (j = 0; j < ndim; j++)
          p_buff[j] += p_i[j] * fac1;
      }
  ytry = (*f) (p_buff);

  if (ytry > y[ilo])
    {
      for (j = 0; j < ndim; j++)
        p_ilo[j] = p_buff[j];
      y[ilo] = ytry;
    }

  return ytry;
}


static double
amoeba (int type, double *p, double *y, int ndim, double (*f) (double *),
        double eps, int ncalls)
{
  int i, j;
  int ilo, ihi;
  int inlo;
  double ysave, ytry;

  for (;;)
    {
      ihi = 0;
      ilo = y[0] < y[1] ? (inlo = 1, 0) : (inlo = 0, 1);
      for (i = 0; i <= ndim; i++)
        {
          if (y[i] >= y[ihi])
            {
              ihi = i;
            }
          if (y[i] < y[ilo])
            {
              inlo = ilo;
              ilo = i;
            }
          else
            {
              if (y[i] < y[inlo] && i != ilo)
                {
                  inlo = i;
                }
            }
        }

      if (ncalls <= 0
          || 2 * (y[ihi] - y[ilo]) / (fabs (y[ilo]) + fabs (y[ihi])) < eps)
        break;

      ytry = amotry (p, y, ndim, f, ilo, -1.0);
      if (type == 1 && !ytry)
        break;
      --ncalls;
      if (ytry >= y[ihi])
        {
          ytry = amotry (p, y, ndim, f, ilo, 2.);
          if (type == 1 && !ytry)
            break;
          --ncalls;
        }

      else if (ytry <= y[inlo])
        {
          ysave = y[ilo];
          ytry = amotry (p, y, ndim, f, ilo, 0.5);
          if (type == 1 && !ytry)
            break;
          --ncalls;
          if (ytry <= ysave)
            {
              for (i = 0; i <= ndim; i++)
                {
                  double *p_ihi = p + ihi * ndim;
                  if (i != ihi)
                    {
                      double *p_i = p + i * ndim;
                      for (j = 0; j < ndim; j++)
                        p_i[j] = 0.5 * (p_i[j] + p_ihi[j]);
                      y[i] = (*f) (p_i);
                      if (type == 1 && !y[i])
                        {
                          for (i = 0; i <= ndim; i++)
                            if (y[i] >= y[ihi])
                              ihi = i;
                          return y[ihi];
                        }
                    }
                }
              ncalls -= ndim;
            }
        }
    }
  for (i = 0; i <= ndim; i++)
    if (y[i] >= y[ihi])
      ihi = i;
  return y[ihi];
}


static double
f_max (double *x)
{
  int i;
  double f;
  double xglobal[MAX_DIM];
  double xlocal[MAX_DIM];

  for (i = 0; i < Ndim; i++)
    {
      xlocal[i] = INCUB (x[i]);
    }
  f = Local2Global (xlocal, xglobal, NULL);
  f *= f_ (xglobal, f);

  if (f < 0)
    return -f;
  else
    return f;
}


static double
run_amoeba (int type, int ndim, double *xx, double *yy, double step,
            double eps, int ncalls)
{
  int i, j;

  for (i = 1; i <= ndim; ++i)
    {
      double *x_i = xx + i * ndim;
      for (j = 0; j < ndim; ++j)
        {
          x_i[j] = xx[j];
        }
      if (x_i[i - 1] > 0.5)
        {
          x_i[i - 1] -= step;
        }
      else
        {
          x_i[i - 1] += step;
        }
      for (j = 0; j < ndim; j++)
        {
          x_i[j] = OUTCUB (x_i[j]);
        }
      yy[i] = f_max (x_i);
    }
  return amoeba (type, xx, yy, ndim, f_max, eps, ncalls);
}


int
vegas_max (vegasGrid * vegPtr, int nCubsINI, int nPoints,
           double (*fxn) (double *, double), double milk, double *eff,
           float *fmax)
{
  int i, j;
  int dim = vegPtr->ndim;
  int stop = 0;
  int cCub = 0;
  double average = 0;
  double *xx = malloc ((dim + 2) * dim * sizeof (double));
  double *yy = malloc ((dim + 2) * sizeof (double));
  long nCubs = generateVegasCubs (vegPtr, nCubsINI);

  f_ = fxn;

  while (cCub < nCubs && !stop)
    {
      fmax[cCub] = 0.;
      for (i = 0; i < nPoints; ++i)
        {
          double f;
          double xglobal[MAX_DIM];
          double xlocal[MAX_DIM];
          drand_arr (dim, xlocal);
          f = Local2Global (xlocal, xglobal, NULL);
          f *= fabs ((*fxn) (xglobal, f));
          average += f;
          if (f > fmax[cCub])
            {
              fmax[cCub] = f;
              for (j = 0; j < dim; j++)
                xx[j] = xlocal[j];
            }
        }

      if (fmax[cCub] && simplexOn) {
          yy[0] = fmax[cCub];
          fmax[cCub] = run_amoeba (0, dim, xx, yy, 0.1, 0.01, nPoints);
        }

      for (i = dim - 1; i >= 0; i--)
        {
          if (++Kg[i] < Ng[i])
            break;
          else
            Kg[i] = 0;
        }
      stop = informline (cCub, nCubs);
      ++cCub;
    }

  if (!stop)
    {
       double minmax;
       double sum = 0;
       for (i = 0; i < nCubs; ++i)
         sum += fmax[i];
       minmax = sum * milk / nCubs;
       for (i = 0; i < nCubs; ++i)
         if (fmax[i] < minmax)
           fmax[i] = minmax;
      *eff = (average / nPoints) / sum;
    }

  free (xx);
  free (yy);
  return stop;
}

extern int vegas_wgt (
                       vegasGrid * vegPtr, 
		       int nCubsINI, 
		       int nEvents,
                       double gmax,
                       double (*fxn) (double *, double), 
		       void (*out) (long, int, double),
                       float *fmax
) 
{
  int i;
  int stop = 0;
  int dim = vegPtr->ndim;
  long cEvent = 0;
  double sum = 0.;
  float * smax = malloc (nCubsINI * sizeof (float));  /* contribution of each Cube to the total CS */
  double * xx = malloc ((dim + 2) * dim * sizeof (double));
  double * yy = malloc ((dim + 2) * sizeof (double));
  long nCubs = generateVegasCubs (vegPtr, nCubsINI);
  int quant = (nEvents - nEvents % 100) / 100;

  for (i = 0; i < nCubs; ++i)
    {
      sum += fmax[i];
      smax[i] = sum;
    }

  for (i = 0; i < nCubs; ++i)
    {
      smax[i] /= sum;
    }

  f_ = fxn;
  while (cEvent < nEvents && !stop)
    {
      long TheCube;
      long L0;
      double f;
      double xglobal[MAX_DIM];
      double xlocal[MAX_DIM];

  /* find a cube number for stratified sampling */
      double rc = drandXX ();
      long Lmin = 0;
      long Lmax = nCubs - 1;
      while (Lmin + 1 < Lmax)
        {
          TheCube = (Lmin + Lmax) / 2;
          if (smax[TheCube] <= rc)
            Lmin = TheCube;
          else
            Lmax = TheCube;
        }
      if (smax[Lmin] > rc)
        {
          TheCube = Lmin;
        }
      else
        {
          TheCube = Lmax;
        }

  /* mysterious transformation... */
      L0 = TheCube;
      for (i = dim - 1; i >= 0; i--)
        {
          Kg[i] = L0 % Ng[i];
          L0 = L0 / Ng[i];
        }

      drand_arr (dim, xlocal);
      f = Local2Global (xlocal, xglobal, NULL);
      f *= (*fxn) (xglobal, f);

      if (fabs (f) <10e-64) continue;

      if (f < 0) f = -f;

      {
        double f_loc = f / fmax[TheCube] / gmax;
        (*out) (TheCube, 1, f_loc);
      }
/*      (*out) (TheCube, 1, f); */

  /* find new maximum, if the calculated value is bigger thet current maximum */
      if (f > fmax[TheCube])
        {
          fmax[TheCube] = f;
          if (f > fmax[TheCube])
            {
              if (simplexOn)
                {
                  double temp_pvect[400];
                  for (i = 0; i < dim; i++)
                    xx[i] = xlocal[i];
                  yy[0] = f;
                  for (i = 0; i < 4 * MAXINOUT; i++)
                    temp_pvect[i] = pvect[i];
                  fmax[TheCube] = run_amoeba (1, dim, xx, yy, 0.1, 0.01, 50);
                  for (i = 0; i < 4 * MAXINOUT; i++)
                    pvect[i] = temp_pvect[i];
                }
              sum = 0;
              for (i = 0; i < nCubs; ++i)
                {
                  sum += fmax[i];
                  smax[i] = sum;
                }
              for (i = 0; i < nCubs; ++i) 
                {
                  smax[i] /= sum;
                }
            }
        }
      if (0 == cEvent % quant) {
        stop = informline (cEvent, nEvents);
      }
      ++cEvent;
    }
  free (xx);
  free (yy);
  free (smax);

  return stop;
}


static double eff;              /* efficiency */
static double rmax;             /* max reached */
static double mult;             /* partion of multiple events */
static double neg;              /* partion of events with negative weght */

double
get_efficiency (void)
{
  return eff;
}

double
get_rmax (void)
{
  return rmax;
}

double
get_multiplicity (void)
{
  return mult;
}

double
get_negativity (void)
{
  return neg;
}


int
vegas_events (vegasGrid * vegPtr, long nCubsINI, long nEvents, double gmax,
              double (*fxn) (double *, double),
              void (*out) (long, int, double), float * fmax)
{
  int i;
  int stop = 0;
  int dim = vegPtr->ndim;
  long cEvent = 0;
  double sum = 0.;

  float * smax = malloc (nCubsINI * sizeof (float));  /* contribution of each Cube to the total CS */
  double * xx = malloc ((dim + 2) * dim * sizeof (double));
  double * yy = malloc ((dim + 2) * sizeof (double));
  long nCubs = generateVegasCubs (vegPtr, nCubsINI);

  f_ = fxn;

  eff = 0.;
  rmax = 0.;
  mult = 0.;
  neg = 0.;

  for (i = 0; i < nCubs; ++i)
    {
      sum += fmax[i];
      smax[i] = sum;
    }

  for (i = 0; i < nCubs; ++i)
    {
      smax[i] /= sum;
    }

  while (cEvent < nEvents && !stop)
    {
      long TheCube;
      double f;
      int n;
      int sgn = 0;
      double xglobal[MAX_DIM];
      double xlocal[MAX_DIM];

  /* find a cube number for stratified sampling */
      double rc = drandXX ();
      long Lmin = 0;
      long Lmax = nCubs - 1;
      while (Lmin + 1 < Lmax)
        {
          TheCube = (Lmin + Lmax) / 2;
          if (smax[TheCube] <= rc)
            Lmin = TheCube;
          else
            Lmax = TheCube;
        }
      if (smax[Lmin] > rc)
        {
          TheCube = Lmin;
        }
      else
        {
          TheCube = Lmax;
        }

  /* mysterious transformation... */
      {
        long L0 = TheCube;
        for (i = dim - 1; i >= 0; i--)
          {
            Kg[i] = L0 % Ng[i];
            L0 = L0 / Ng[i];
          }
     }

      drand_arr (dim, xlocal);
      f = Local2Global (xlocal, xglobal, NULL);
      f *= (*fxn) (xglobal, f);

      if (fabs (f) <10e-64) continue;

      if (f < 0)
        {
          f = -f;
          sgn = 1;
        }
      ++eff;

  /* MC re-weighting */
      {
        double f_loc = f / fmax[TheCube];
        if (f_loc > rmax)
          rmax = f_loc;
        f_loc /= gmax;
        n = (int) (f_loc);
        f_loc = f_loc - n;
        if (f_loc > drandXX ())
          n++;
        if (n)
          {
            cEvent += n;
            if (cEvent > nEvents)
              {
                n -= (cEvent - nEvents);
              }
            mult += n - 1;
            if (sgn)
              {
                n *= -1;
                neg += n;
              }
            (*out) (TheCube, n, 1.);
          }
      }

  /* find new maximum, if the calculated value is bigger thet current maximum */
      if (f > fmax[TheCube])
        {
          fmax[TheCube] = f;
          if (f > gmax * fmax[TheCube])
            {
              if (simplexOn)
                {
                  double temp_pvect[400];
                  for (i = 0; i < dim; i++)
                    xx[i] = xlocal[i];
                  yy[0] = f;
                  for (i = 0; i < 4 * MAXINOUT; i++)
                    temp_pvect[i] = pvect[i];
                  fmax[TheCube] = run_amoeba (1, dim, xx, yy, 0.1, 0.01, 50);
                  for (i = 0; i < 4 * MAXINOUT; i++)
                    pvect[i] = temp_pvect[i];
                }
              sum = 0;
              for (i = 0; i < nCubs; ++i)
                {
                  sum += fmax[i];
                  smax[i] = sum;
                }
              for (i = 0; i < nCubs; ++i) 
                {
                  smax[i] /= sum;
                }
            }
        }
      stop = informline (cEvent, nEvents);
    }

  neg = neg / cEvent;
  mult = mult / cEvent;
  eff = nEvents / eff;
  free (xx);
  free (yy);
  free (smax);
  return stop;
}


int
vegas_1to2_events (vegasGrid * vegPtr, long nCubsINI, long nEvents, double gmax,
              double (*fxn) (double *, double),
              void (*out) (long, int, double), float * fmax)
{
  int stop = 0;
  long cEvent = 0;

  eff = 1.;
  rmax = 1.;
  mult = 0.;
  neg = 0.;

  while (cEvent < nEvents && !stop) {
    double xglobal[MAX_DIM];
    double test_value = fxn (xglobal, 1.);
    if (test_value < 0.)
      neg = 1.;
    ++cEvent;
    (*out) (1, 1, 1.);
    stop = informline (cEvent, nEvents);
  }

  return stop;
}
