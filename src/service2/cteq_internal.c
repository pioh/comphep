/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static void
skipline (int l)
{
  char c;
  for (; l; l--)
    {
      do
	scanf ("%c", &c);
      while (c != '\n');
    }
}

static double
alpha (int nf, int ord, double lambda, double dscale)
{

  double d__2, d__4;

  double b0 = 11. - (2. / 3.) * nf;
  double b1 = 51. - (19. / 3.) * nf;
  double b2 = 2857. - (5033. / 9.) * nf + (325. / 27.) * nf * nf;

  double rl = 2 * log (dscale / lambda);
  double alpha0 = 4 * M_PI / (b0 * rl);

  d__4 = log (rl) - .5;
  d__2 = 2 * b1 / (b0 * b0 * rl);

  if (ord == 1)
    return alpha0;
  else if (ord == 2)
    return alpha0 * (1 - 2 * b1 * log (rl) / (b0 * b0 * rl));
  else if (ord == 3)
    return alpha0 * (1 - 2 * b1 * log (rl) / (b0 * b0 * rl)
		     + d__2 * d__2 * (d__4 * d__4 + b2 * b0 / (8 * b1 * b1) -
				      1.25));
  else
    {
      fprintf (stderr, "Can not evaluate alpha in so large precision.\n");
      exit (1);
    }
}


int
main (int argc, char **argv)
{
  double buff;
  int i, j, k;
  int nx, nt, NfMx;
  char c;
  double *q;
  char version[30];
  int odr, nf;
  double odrf, nff;
  double lambda;
  char names[4][10] = { "(t,T)", "(b,B)", "(c,C)", "(s,S)" };

  do
    scanf ("%c", &c);
  while (c != ':');
  scanf ("%s", version);
  skipline (2);
  scanf ("%lf %lf %lf", &odrf, &nff, &lambda);
  nf = nff;
  odr = odrf;
  skipline (2);
  scanf ("%d %d %d", &nx, &nt, &NfMx);
  printf ("#distribution \"%s(proton)\"     ", version);
  for (i = 6 - NfMx; i < 4; i++)
    printf (" %s", names[i]);
  printf (" D U G u d \n");

  printf ("#distribution \"%s(anti-proton)\"", version);
  for (i = 6 - NfMx; i < 4; i++)
    printf (" %s", names[i]);
  printf (" d u G U D\n");

  skipline (2);
  scanf ("%lf", &buff);
  printf ("\n#q_min %.5E\n", buff);
  scanf ("%lf", &buff);
  printf ("\n#q_max %.5E\n", buff);
  skipline (1);
  q = (double *) malloc (sizeof (double) * (nt + 1));


  for (i = 0; i <= nt; i++)
    scanf ("%lf", q + i);
  printf ("\n#Q_grid\n");
  for (i = 0, j = 1; i <= nt; i++, j++)
    {
      printf (" %.5E", q[i]);
      if (j == 10)
	{
	  printf ("\n");
	  j = 0;
	}
    }
  printf ("\n");

  printf ("\n#Alpha\n");
  for (i = 0, j = 1; i <= nt; i++, j++)
    {
      printf (" %.5E", alpha (nf, odr, lambda, q[i]));
      if (j == 10)
	{
	  printf ("\n");
	  j = 0;
	}
    }
  printf ("\n");


  skipline (2);
  scanf ("%lf", &buff);
  printf ("\n#x_min %.5E\n", buff);
  skipline (1);
  printf ("\n#X_grid\n");
  for (i = 0, j = 1; i <= nx; i++, j++)
    {
      scanf ("%lf", &buff);
      printf (" %.5E ", buff);
      if (j == 10)
	{
	  printf ("\n");
	  j = 0;
	}
    }
  printf ("\n");

  skipline (2);

  for (k = 0; k < 3 + NfMx; k++)
    {
      printf ("\n#%d-parton\n", k + 1);
      for (j = 0; j < nt + 1; j++)
	{
	  for (i = 0; i < nx + 1; i++)
	    {
	      scanf ("%lf", &buff);
	      printf (" %.5E", buff);
	    }
	  printf ("\n");
	}
    }
  return 0;
}
