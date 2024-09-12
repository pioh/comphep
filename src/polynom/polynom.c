/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "service2/include/getmem.h"

#include "polynom.h"

static poly garbage;

int levi = 0;
poly *contracts = NULL;
int monomLength = 0;
int maxLength = 0;

void set_garbage (poly g) {
  garbage = g;
}

poly get_garbage (void) {
  return garbage;
}

void (*memoryInfo) (int) = NULL;

void 
newmonom (poly * p)
{
  poly buffpoly;
  unsigned i /*, imax */ ;

  if (garbage != NULL)
    {
      *p = garbage;
      garbage = garbage->next;
    }
  else
    {
      int brest, memsize;

      memsize = sizeof (monom) + sizeof (long) * (maxLength - 2);
      *p = (poly) getmem_ (memsize);
      brest = blockrest (memsize);

      for (i = 1; i <= brest; i++)
	{
	  buffpoly = garbage;
	  garbage = (poly) getmem_ (memsize);
	  garbage->next = buffpoly;
	}
      if (memoryInfo != NULL)
	(*memoryInfo) (usedmemory);
    }
}

void 
delmonom (poly * p)
{
  (*p)->next = garbage;
  garbage = *p;
  *p = NULL;
}

void 
delpoly (poly * p)
{
  poly m, mm;
  if (*p == NULL)
    return;
  m = *p;
  mm = m->next;
  while (mm != NULL)
    {
      m = mm;
      mm = mm->next;
    }
  m->next = garbage;
  garbage = *p;
  *p = NULL;
}

poly 
plusone (void)
{
  poly p;
  int i;
  newmonom (&p);
  p->next = NULL;
  p->coef.num = NUM_ONE;
  for (i = 0; i < maxLength; i++)
    p->tail.power[i] = 0;
  return p;
}

poly 
copypoly (poly p)
{
  poly pp, qq, copypoly1;
  int i;

  if (p == NULL)
    return NULL;
  else
    {
      newmonom (&pp);
      copypoly1 = pp;
    label_1:
      for (i = 0; i < monomLength; i++)
	pp->tail.power[i] = p->tail.power[i];
      pp->coef.num = p->coef.num;
      p = p->next;
      if (p != NULL)
	{
	  newmonom (&qq);
	  pp->next = qq;
	  pp = qq;
	  goto label_1;
	}
      pp->next = NULL;
    }
  return copypoly1;
}

static void 
comparemonoms (poly m1, poly m2, int *gt, int *eq)
{
  int i = 0;
  int noteq;

  do
    {
      ++(i);
      noteq = (m1->tail.power[i - 1] != m2->tail.power[i - 1]);
    }
  while (!(noteq || i == monomLength));
  *eq = !noteq;
  *gt = (m1->tail.power[i - 1] > m2->tail.power[i - 1]);
}


void 
sewpoly (poly * p1, poly * p2)
{
  poly m, mm, m1, m2;
  int gt, eq;


  /* Nested function: comparemonoms */


  if (*p2 == NULL)
    return;
  if (*p1 == NULL)
    {
      *p1 = *p2;
      *p2 = NULL;
      return;
    }
  m1 = *p1;
  m2 = *p2;
  newmonom (p1);
  (*p1)->next = m1;
  m = *p1;
  *p2 = NULL;

label_1:			/*  (M1,M2,GT,EQ)  */
  comparemonoms (m1, m2, &gt, &eq);

label_2:
  while (gt)
    {
      m = m1;
      m1 = m1->next;
      if (m1 == NULL)
	{
	  m->next = m2;
	  goto label_3;
	}
      /*  (M1,M2,GT,EQ)  */
      comparemonoms (m1, m2, &gt, &eq);
    }
  if (eq)
    {
      m1->coef.num += m2->coef.num;
      mm = m2;
      m2 = m2->next;
      delmonom (&mm);
      if (m1->coef.num == NUM_ZERO)
	{
	  mm = m1;
	  m1 = m1->next;
	  delmonom (&mm);
	  if (m1 == NULL)
	    {
	      m->next = m2;
	      goto label_3;
	    }
	  m->next = m1;
	}
      if (m2 == NULL)
	goto label_3;
      goto label_1;
    }
  mm = m1;
  m1 = m2;
  m2 = mm;
  m->next = m1;
  gt = TRUE;
  goto label_2;

label_3:
  mm = *p1;
  *p1 = (*p1)->next;
  delmonom (&mm);
/*   ;TestPoly(P1)  */
}

void 
multpolyint (poly * p, long i)
{
  poly pp;

  if (i == 0)
    delpoly (p);
  else
    {
      pp = *p;
      while (pp != NULL)
	{
	  pp->coef.num *= i;
	  pp = pp->next;
	}
    }
}

static poly 
multpolymono (poly plnm, poly mono)
{
  poly pp, qq, multpolymono1;
  int i;

  if (plnm == NULL || mono == NULL)
    return NULL;
  newmonom (&pp);
  multpolymono1 = pp;

label_1:
  pp->coef.num = plnm->coef.num * mono->coef.num;
  for (i = 0; i < monomLength; i++)
    pp->tail.power[i] = plnm->tail.power[i] + mono->tail.power[i];
  plnm = plnm->next;
  if (plnm != NULL)
    {
      newmonom (&qq);
      pp->next = qq;
      pp = qq;
      goto label_1;
    }
  pp->next = NULL;
  return multpolymono1;
}

/* -------------------------------------------------- */

poly 
multtwopoly (poly q1, poly q2)
{
  poly mlttwpl, mltplmn, p1, p2;
  mlttwpl = NULL;
  if (q1 != NULL && q2 != NULL)
    {
      p1 = q1;
      p2 = q2;
      for (;;)
	{
	  p1 = p1->next;
	  p2 = p2->next;
	  if (p1 == NULL)
	    {
	      p1 = q2;
	      p2 = q1;
	      break;
	    }
	  if (p2 == NULL)
	    {
	      p1 = q1;
	      p2 = q2;
	      break;
	    }
	}

      mlttwpl = multpolymono (p1, p2);
      for (p2 = p2->next; p2; p2 = p2->next)
	{
	  mltplmn = multpolymono (p1, p2);
	  sewpoly (&mlttwpl, &mltplmn);
	}
    }
  return mlttwpl;
}



/* ---------- Common ----------- */

poly 
scalarmult (int p1, int p2)
{
  unsigned c, cc, n;

  if (p1 < p2)
    {
      c = -p1;
      cc = -p2;
    }
  else
    {
      c = -p2;
      cc = -p1;
    }
  n = cc + c * (c - 1) / 2;
  return contracts[n - 1];
}				/*  ScalarMult  */

void 
assignsclmult (int p1, int p2, poly p)
{
  unsigned c, cc, n;

  if (p1 < p2)
    {
      c = -p1;
      cc = -p2;
    }
  else
    {
      c = -p2;
      cc = -p1;
    }
  n = cc + c * (c - 1) / 2;
  contracts[n - 1] = p;
}				/*  ScalarMult  */

void 
deltensor (poly * t)
{
  poly m, mm;

  if (*t == NULL)
    return;
  mm = *t;
  do
    {
      m = mm;
      delpoly (&m->coef.complex.re);
      delpoly (&m->coef.complex.im);
      mm = m->next;
    }
  while (mm != NULL);
  m->next = garbage;
  garbage = *t;
  *t = NULL;
}

poly 
copytens (poly t, int ln)
{
  poly tt, qq;
  int i;
  poly copytens1;

  if (t == NULL)
    return NULL;
  newmonom (&tt);
  copytens1 = tt;

label_1:
  for (i = 0; i < ln; i++)
    tt->tail.power[i] = t->tail.power[i];
  tt->coef.complex.re = (poly) copypoly (t->coef.complex.re);
  tt->coef.complex.im = (poly) copypoly (t->coef.complex.im);
  t = t->next;
  if (t != NULL)
    {
      newmonom (&qq);
      tt->next = qq;
      tt = qq;
      goto label_1;
    }
  tt->next = NULL;
  return copytens1;
}

static void 
comparetensbody (poly m1, poly m2, int *gt, int *eq,
		 int l)
{
  int i = 0;
  int noteq;

  do
    {
      noteq = (m1->tail.power[i] != m2->tail.power[i]);
      ++(i);
    }
  while (!(noteq || i == l));
  *eq = !noteq;

  *gt = *eq ? FALSE : m1->tail.power[i - 1] > m2->tail.power[i - 1];
}				/*  CompareTensBody  */


void 
sewtens (poly * t1, poly * t2, int ln)
{
  poly m, mm, m1, m2;
  int gt, eq;

  if (*t2 == NULL)
    return;
  if (*t1 == NULL)
    {
      *t1 = *t2;
      *t2 = NULL;
      return;
    }
  m1 = *t1;
  m2 = *t2;
  newmonom (t1);
  (*t1)->next = m1;
  m = *t1;
  *t2 = NULL;

label_1:			/*  (M1,M2,GT,EQ)  */
  comparetensbody (m1, m2, &gt, &eq, ln);

label_2:
  while (gt)
    {
      m = m1;
      m1 = m1->next;
      if (m1 == NULL)
	{
	  m->next = m2;
	  goto label_3;
	}
      /*  (M1,M2,GT,EQ)  */
      comparetensbody (m1, m2, &gt, &eq, ln);
    }
  if (eq)
    {
      sewpoly (&m1->coef.complex.re, &m2->coef.complex.re);
      sewpoly (&m1->coef.complex.im, &m2->coef.complex.im);

      mm = m2;
      m2 = m2->next;
      delmonom (&mm);
      if (m1->coef.complex.re == NULL && m1->coef.complex.im == NULL)
	{
	  mm = m1;
	  m1 = m1->next;
	  delmonom (&mm);
	  if (m1 == NULL)
	    {
	      m->next = m2;
	      goto label_3;
	    }
	  m->next = m1;
	}
      if (m2 == NULL)
	goto label_3;
      goto label_1;
    }
  mm = m1;
  m1 = m2;
  m2 = mm;
  m->next = m1;
  gt = TRUE;
  goto label_2;

label_3:
  mm = *t1;
  *t1 = (*t1)->next;
  delmonom (&mm);
}

void 
multtensint (poly * t, long i)
{
  poly tt;

  if (i == 0)
    deltensor (t);
  else
    {
      tt = *t;
      while (tt != NULL)
	{
	  multpolyint (&tt->coef.complex.re, i);
	  multpolyint (&tt->coef.complex.im, i);
	  tt = tt->next;
	}
    }
}

void 
multtenspoly (poly * t, poly p)
{
  poly tt;
  poly pp;

  if (p == NULL)
    deltensor (t);
  else
    {
      tt = *t;
      while (tt != NULL)
	{
	  pp = (poly) multtwopoly (tt->coef.complex.re, p);
	  delpoly (&tt->coef.complex.re);
	  tt->coef.complex.re = pp;

	  pp = (poly) multtwopoly (tt->coef.complex.im, p);
	  delpoly (&tt->coef.complex.im);
	  tt->coef.complex.im = pp;
	  tt = tt->next;
	}
    }
}

void 
multtensComplexpoly (poly * t, poly re, poly im)
{
  poly tt;
  poly pRe, pIm, qq;

  if (re == NULL && im == NULL)
    deltensor (t);
  else
    {
      tt = *t;
      while (tt != NULL)
	{
	  pRe = (poly) multtwopoly (tt->coef.complex.re, re);
	  qq = (poly) multtwopoly (tt->coef.complex.im, im);
	  multpolyint (&qq, -1);
	  sewpoly (&pRe, &qq);

	  pIm = (poly) multtwopoly (tt->coef.complex.re, im);
	  qq = (poly) multtwopoly (tt->coef.complex.im, re);

	  sewpoly (&pIm, &qq);

	  delpoly (&tt->coef.complex.re);
	  tt->coef.complex.re = pRe;

	  delpoly (&tt->coef.complex.im);
	  tt->coef.complex.im = pIm;

	  tt = tt->next;
	}
    }
}

void 
tensRealPart (poly * t)
{
  poly first, *pred, tt;

  if (!t || !(*t))
    return;
  tt = *t;
  first = tt;
  pred = &first;
  while (tt)
    {
      delpoly (&(tt->coef.complex.im));
      if (tt->coef.complex.re)
	{
	  pred = &(tt->next);
	  tt = tt->next;
	}
      else
	{
	  *pred = tt->next;
	  delmonom (&tt);
	  tt = *pred;
	}
    }
  *t = first;
}
