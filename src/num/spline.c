/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Copyright (C) 2000, Slava Bunichev
* ------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include "spline.h"

#define N2 (N+1)*(N+1)
#define SET_MEMORY(array,N)  (array)= (double *)malloc((N) *sizeof(double));  for(j=0;j<(N);j++) (array)[j]=0;
#define SET_P(N)  p= (int *)malloc((N) *sizeof(int)); for(j=0;j<(N);j++) p[j]=j*(N);

#define A(i,j)          A[p[(i)]+(j)]
#define C(i,j)          C[p[(i)]+(j)]
#define F3(i,j)         F3[p[(i)]+(j)]
#define matrix_C(i,j)   matrix_C[p[(i)]+(j)]
#define matrix_F3(i,j)  matrix_F3[p[(i)]+(j)]

static int *p;


extern void 
progonca (int N, double *f, double *b, double *c, double *d) 
{
  int j;
  double *alfa, *beta;

  SET_MEMORY (alfa, N + 1) 
  SET_MEMORY (beta, N + 1) 

  alfa[2] = -0.5;
  beta[2] = -0.25 * (f[0] + 4 * f[1] - 5 * f[2]);
  for (j = 2; j <= N - 2; j++)
    {
      alfa[j + 1] = -1. / (4. + alfa[j]);
      beta[j + 1] = (beta[j] + 3 * (f[j - 1] - f[j + 1])) * alfa[j + 1];
    }

  b[N - 1] = (-0.5 * beta[N - 1] + 0.25 * (f[N] + 4 * f[N - 1] - 5 * f[N - 2])) 
              /(1 + 0.5 * alfa[N - 1]);
  for (j = N - 2; j >= 1; j--)
    b[j] = alfa[j + 1] * b[j + 1] + beta[j + 1];
  b[0] = b[2] - 2 * (f[0] - 2 * f[1] + f[2]);
  b[N] = b[N - 2] + 2 * (f[N - 2] - 2 * f[N - 1] + f[N]);

  for (j = 1; j <= N; j++)
    c[j] = 4 * b[j] + 2 * b[j - 1] + 6 * (f[j - 1] - f[j]);

  for (j = 2; j <= N; j++)
    d[j] = c[j] - c[j - 1];
  d[1] = d[2];

  free (alfa);
  free (beta);
}


double 
spline_for_graph (double iks, double *f, double *b, double *c, double *d) 
{
  int i = 0;
  double delta_x;

  while (i < iks)
    i += 1;
  delta_x = iks - i;
  return (f[i] + b[i] * delta_x + c[i] / 2. * delta_x * delta_x + d[i] / 6. * delta_x * delta_x * delta_x);
}


static void 
GAUSS (double psi, int N, double *B, double *C, 
       double *matrix_B, double *matrix_C, double *matrix_F3) 
{
  int i, j, k, l, s, t;
  double ch, *constant;

  SET_MEMORY (constant, N + 1) 
  for (j = 0; j <= N; j++)
    {
      B[j] = matrix_B[j];
      for (i = 0; i <= N; i++)
        C (j, i) = matrix_C (j, i) + psi * matrix_F3 (j, i);
    }
  for (k = 0; k <= N - 1; k++)
    {
      if (k >= (N - 4))
        l = N;
      else
        l = k + 4;

      ch = C (k, k);
      t = k;
      for (j = k + 1; j <= l; j++)
        if (ch < C (j, k))
          {
            ch = C (j, k);
            t = j;
          }
      if (t != k)
        {
          for (i = k; i <= N; i++)
            {
              ch = C (t, i);
              C (t, i) = C (k, i);
              C (k, i) = ch;
            }
          ch = B[t];
          B[t] = B[k];
          B[k] = ch;
        }

      if (k >= (N - 8))
        s = N;
      else
        s = k + 8;
      for (i = s; i >= k; i--)
        constant[i] = C (k, i) / C (k, k);
      for (j = k + 1; j <= l; j++)
        {
          B[j] -= B[k] * C (j, k) / C (k, k);
          for (i = s; i >= k; i--)
            C (j, i) -= C (j, k) * constant[i];
        }
    }

  B[N] /= C (N, N);
  for (j = N - 1; j >= 0; j--)
    {
      if ((N - j) >= 8)
        s = j + 8;
      else
        s = N;
      for (i = s; i >= j + 1; i--)
        B[j] -= C (j, i) * B[i];
      B[j] /= C (j, j);
    }
  free (constant);
}


static void 
NORMIR (int n, double *Inor, double *Iexp) 
{
  int j, N;

  N = 2 * n + 2;
  Inor[0] = Iexp[0];
  for (j = 0; j < n; j++)
    {
      Inor[2 * j + 1] = Iexp[j];
      Inor[2 * j + 2] = Iexp[j] + Iexp[j + 1];
    }
  Inor[N - 1] = Iexp[n];
  Inor[N] = Iexp[n];
  for (j = 0; j <= N; j++)
    if (Inor[j] == 0)
      Inor[j] = 1.e-50;
}


static void 
INTEGRALS (int n, double *Iexp, double *dIexp, double *dIexp_2, 
           double *A, double *matrix_B, double *matrix_C, 
           double *matrix_F3) 
{
  int i, j, N, g, v;
  double *F3, *Inor;

  N = 2 * n + 2;
  SET_MEMORY (F3, N2) 
  SET_MEMORY (Inor, N + 1) 
  NORMIR (n, Inor, Iexp);

  for (j = 0; j <= n; j++)
    {
      if (dIexp[j] == 0)
        dIexp_2[j] = 1.e-50;
      else
        dIexp_2[j] = dIexp[j] * dIexp[j];
    }

  A (0, 0) = 0.375;
  for (j = 0; j < n; j++)
    {
      i = j * 2 + 2;
      A (i - 1, j) = 0.6875;
      A (i - 1, j + 1) = 0.03125;
      A (i, j) = 0.375;
      A (i, j + 1) = 0.375;
      A (i + 1, j) = 0.03125;
      A (i + 1, j + 1) = 0.6875;
    }
  A (N, n) = 0.375;

  for (i = 1; i < N; i++)
    {
       F3 (i - 1, i) = -1. / Inor[i];
       F3 (i, i) = 2. / Inor[i];
       F3 (i + 1, i) = -1. / Inor[i];
       i += 1;
    }
  for (g = 0; g <= N; g++)
    for (v = 0; v <= N; v++)
      for (j = 0; j <= N; j++)
        matrix_F3 (g, v) += F3 (g, j) * F3 (v, j);

  for (i = 0; i < N2; i++)
    F3[i] = 0;
  F3 (0, 0) = 1. / Inor[0];
  F3 (1, 0) = -1. / Inor[0];
  for (i = 2; i <= N - 2; i++)
    {
      F3 (i - 2, i) = -1. / Inor[i];
      F3 (i, i) = 2. / Inor[i];
      F3 (i + 2, i) = -1. / Inor[i];
      i = i + 1;
    }
  F3 (N - 1, N) = -1. / Inor[N];
  F3 (N, N) = 1. / Inor[N];
  for (g = 0; g <= N; g++)
    for (v = 0; v <= N; v++)
      for (j = 0; j <= N; j++)
        matrix_F3 (g, v) += 0.8 * F3 (g, j) * F3 (v, j);

  for (g = 0; g <= N; g++)
    for (v = 0; v <= N; v++)
      for (j = 0; j <= n; j++)
        matrix_C (g, v) += A (g, j) * A (v, j) / dIexp_2[j];
  for (g = 0; g <= N; g++)
    for (j = 0; j <= n; j++)
      matrix_B[g] += A (g, j) * Iexp[j] / dIexp_2[j];

  free (F3);
  free (Inor);
}


static double 
KSI_2 (double psi, int n, double *Iexp, double *dIexp_2, 
       double *A, double *B, double *C, 
       double *matrix_B, double *matrix_C, double *matrix_F3) 
{
  int i, j, s, l, N;
  double *Iteor, ksi2 = 0;

  N = 2 * n + 2;
  SET_MEMORY (Iteor, n + 1) 
  GAUSS (psi, N, B, C, matrix_B, matrix_C, matrix_F3);

  for (j = 0; j <= n; j++)
    {
      if (j == 0)
        s = 0;
      else
        s = 2 * j - 1;
      if (j == n)
        l = N;
      else
        l = 2 * j + 3;
      for (i = s; i <= l; i++)
        Iteor[j] += B[i] * A (i, j);
  }

  for (j = 0; j <= n; j++)
    ksi2 += (Iteor[j] - Iexp[j]) * (Iteor[j] - Iexp[j]) / dIexp_2[j];
  free (Iteor);

  return (ksi2);
}


static void 
dikotomia (double *finale_psi, double *finale_ksi2, int n, 
           double *Iexp, double *dIexp_2, double *A, double *B, 
           double *C, double *matrix_B, double *matrix_C, 
           double *matrix_F3, double koeff) 
{
  double limit_ksi;
  double psi;
  double psi_left, psi_right;
  double ksi_left, ksi_right;
  double delta, delta_ksi, ksi_psi;

  limit_ksi = koeff * (n + 1);

  psi = 1.;

  delta_ksi = limit_ksi * 0.1;

  psi_right = psi;
  psi_left = psi;

  ksi_right = KSI_2 (psi_right, n, Iexp, dIexp_2, A, B, C, matrix_B, matrix_C, matrix_F3);
  ksi_left = ksi_right;

  delta = ksi_right - limit_ksi;

  if (delta < 0)
    delta = -delta;

  if (delta > delta_ksi)
    {
      if (ksi_right < limit_ksi)
        while (ksi_right < limit_ksi)
          {
            psi_left = psi_right;
            psi_right = psi_right * 100.;
            printf ("psi_right= %e\n", psi_right);
            ksi_left = ksi_right;
            ksi_right = KSI_2 (psi_right, n, Iexp, dIexp_2, A, B, C, matrix_B, matrix_C, matrix_F3);
          }
      else
        while (ksi_left > limit_ksi)
          {
            psi_right = psi_left;
            psi_left = psi_left * 0.01;
            printf ("psi_left= %e\n", psi_left);
            ksi_right = ksi_left;
            ksi_left = KSI_2 (psi_left, n, Iexp, dIexp_2, A, B, C, matrix_B, matrix_C, matrix_F3);
          }
    }

    printf ("\npsi_left= %e\n", psi_left);
    printf ("psi_right= %e\n\n", psi_right);

    psi = (psi_left + psi_right) / 2.0;
    printf ("psi= %e\n", psi);
    ksi_psi = KSI_2 (psi, n, Iexp, dIexp_2, A, B, C, matrix_B, matrix_C, matrix_F3);
    delta = ksi_psi - limit_ksi;
    if (delta < 0)
      delta = -delta;

    while (delta > delta_ksi)
      {
        psi = (psi_left + psi_right) / 2.0;
        printf ("psi= %e\n", psi);
        ksi_psi = KSI_2 (psi, n, Iexp, dIexp_2, A, B, C, matrix_B, matrix_C, matrix_F3);
        delta = ksi_psi - limit_ksi;
        if (delta < 0)
          {
            delta = -delta;
            psi_left = psi;
          }
        else
          psi_right = psi;
      }

    *finale_psi = psi;
    *finale_ksi2 = ksi_psi;
}


void 
SPLINE (double koeff, int n, double *Iexp, double *dIexp, double *f, 
        double *b, double *c, double *d) 
{
  int i, j, N;
  double psi, ksi2, *dIexp_2, *A, *B, *C, *matrix_B, *matrix_C, *matrix_F3;

  N = 2 * n + 2;
  SET_P (N + 1) 
  SET_MEMORY (dIexp_2, N + 1) 
  SET_MEMORY (B, N + 1) 
  SET_MEMORY (A, N2) 
  SET_MEMORY (C, N2) 
  SET_MEMORY (matrix_B, N + 1) 
  SET_MEMORY (matrix_C, N2) 
  SET_MEMORY (matrix_F3, N2) 
  INTEGRALS (n, Iexp, dIexp, dIexp_2, A, matrix_B, matrix_C, matrix_F3);
  dikotomia (&psi, &ksi2, n, Iexp, dIexp_2, A, B, C, matrix_B, matrix_C, matrix_F3, koeff);

  /* ksi2= KSI_2(psi,n,Iexp,dIexp_2,A,B,C,matrix_B,matrix_C,matrix_F3); */ 
  printf ("ksi2= %e\n", ksi2);

  for (i = 0; i <= 2 * n; i++)
    {
      f[i] = 0.25 * B[i] + B[i + 1] + 0.25 * B[i + 2];
      b[i] = -0.75 * B[i] + 0.75 * B[i + 2];
      c[i] = 1.5 * B[i] - 3 * B[i + 1] + 1.5 * B[i + 2];
    }

  for (i = 1; i <= 2 * n; i++)
    {
      d[i] = c[i] - c[i - 1];
    }

    free (dIexp_2);
    free (A);
    free (B);
    free (C);
    free (matrix_B);
    free (matrix_C);
    free (matrix_F3);
    free (p);
}
