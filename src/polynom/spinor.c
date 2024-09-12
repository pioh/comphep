/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* ------------------------------------------------------
*/
#include <stdio.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"

#include "polynom.h"
#include "tensor.h"
#include "spinor.h"

int spinLength = 0;

static poly spnr;		/* from calcspur */
static int nused;		/* from calcspur */
static int used[118];		/* from calcspur */
static poly t1_;		/* from multspinmemb */
static poly t2_;		/* from multtwospin */
static int forspur;		/* from multtwospin */
static char indlist[119];	/* from multtwospin */

static poly 
calconespur (void)		/*  (Nused,Used,Spnr), recursion  */
{
  int nu;
  int sign, si, sk;
  int i, k;
  poly ans, r1, r2;

  nu = spnr->tail.spin.l;
  if (nu == nused)		/*  IF last step  */
    {
      newmonom (&ans);
      ans->coef.complex.re = (poly) plusone ();
      ans->coef.complex.im = NULL;
      for (i = 0; i < tensLength; i++)
	ans->tail.power[i] = 0;
      ans->next = NULL;
      return ans;
    }				/*  IF last step  */

  ans = NULL;
  k = 1;
  while (used[k - 1])
    ++(k);
  used[k - 1] = TRUE;
  nused += 2;
  sk = spnr->tail.spin.g[k - 1];
  sign = 1;
  for (i = k + 1; i <= spnr->tail.spin.l; i++)
    if (!used[i - 1])
      {
	si = spnr->tail.spin.g[i - 1];
	used[i - 1] = TRUE;
	r1 = calconespur ();
	used[i - 1] = FALSE;

	if (r1 != NULL)
	  {
	    if (sk < 0 && si < 0)
	      multtenspoly (&r1, scalarmult (sk, si));
	    else
	      {
		r2 = r1;
		while (r2 != NULL)
		  {
		    if (si > 0)
		      r2->tail.tens[si - 1] = sk;
		    if (sk > 0)
		      r2->tail.tens[sk - 1] = si;
		    r2 = r2->next;
		  }
	      }
	  }
	if (sign == -1)
	  multtensint (&r1, -1);
	sign = -sign;
	sewtens (&ans, &r1, tensLength);
      }				/*  For i:=k+1... */
  used[k - 1] = FALSE;
  nused -= 2;
  return ans;
}


poly 
calcspur (poly spnr1)
{
  poly summ, onesp;
  int i, i1, i2, i3, i4, j;
  int c, cc, sgn;
  char ee[4];
  poly tmpsp, tmpsp2;


  /* Nested function: calconespur */
  /*  CalcSpur  */
  spnr = spnr1;
  summ = NULL;
  while (spnr != NULL)
    {
      if (spnr->tail.spin.g5 == 0)
	{
	  nused = 0;
	  for (j = 1; j <= spnr->tail.spin.l; j++)
	    used[j - 1] = FALSE;	/*  for recursion  */
	  /*  (Spnr)  */
	  onesp = calconespur ();
	}
      else
	{
	  if (spnr->tail.spin.l < 4)
	    onesp = NULL;
	  else
	    {
	      onesp = NULL;
	      for (j = 1; j <= spnr->tail.spin.l; j++)
		used[j - 1] = FALSE;
	      for (i1 = 1; i1 <= spnr->tail.spin.l - 3; i1++)
		for (i2 = i1 + 1; i2 <= spnr->tail.spin.l - 2; i2++)
		  for (i3 = i2 + 1; i3 <= spnr->tail.spin.l - 1; i3++)
		    for (i4 = i3 + 1; i4 <= spnr->tail.spin.l; i4++)
		      {
			nused = 4;
			for (j = 1; j <= spnr->tail.spin.l; j++)
			  used[j - 1] = FALSE;
			ee[0] = spnr->tail.spin.g[i1 - 1];
			used[i1 - 1] = TRUE;
			ee[1] = spnr->tail.spin.g[i2 - 1];
			used[i2 - 1] = TRUE;
			ee[2] = spnr->tail.spin.g[i3 - 1];
			used[i3 - 1] = TRUE;
			ee[3] = spnr->tail.spin.g[i4 - 1];
			used[i4 - 1] = TRUE;
			j = 1;
			sgn = ((i1 + i2 + i3 + i4) & 1) == 0 ? 1 : -1;
			do
			  {
			    c = ee[j - 1];
			    cc = ee[j + 1 - 1];
			    if (cc < c)
			      ++(j);
			    else if (cc > c)
			      {
				ee[j - 1] = cc;
				ee[j + 1 - 1] = c;
				sgn = -sgn;
				if (j > 1)
				  --(j);
				else
				  ++(j);
			      }
			    else
			      sgn = 0;
			  }
			while (!(j == 4 || sgn == 0));
			if (sgn != 0)
			  {
			    tmpsp = calconespur ();
			    if (sgn == -1)
			      multtensint (&tmpsp, -1);
			    tmpsp2 = tmpsp;
			    while (tmpsp2 != NULL)
			      {
				for (i = 0; i < 4; i++)
				  tmpsp2->tail.tens[maxIndex + i] = ee[i];
				tmpsp2->coef.complex.im = tmpsp2->coef.complex.re;
				tmpsp2->coef.complex.re = NULL;
				tmpsp2 = tmpsp2->next;
			      }
			    sewtens (&onesp, &tmpsp, tensLength);
			  }
		      }
	    }
	}
      multtensComplexpoly (&onesp, spnr->coef.complex.re, spnr->coef.complex.im);
      sewtens (&summ, &onesp, tensLength);
      onesp = spnr;
      spnr = spnr->next;
      delpoly (&(onesp->coef.complex.re));
      delpoly (&(onesp->coef.complex.im));
      delmonom (&onesp);

    }
/*    MultTensInt(Summ,4); */
  return summ;
}				/*  CalcSpure  */


static void 
kahane (poly * spnrs, int ind)
{
  int m1, m2, nrevol;
  int c, cc, ll, i, j;
  poly spn1, spn2;

  spn1 = *spnrs;
  while (spn1 != NULL)
    {
      m1 = 1;
      while (spn1->tail.spin.g[m1 - 1] != ind)
	++(m1);
      /*  first index position   */
      m2 = m1 + 1;
      while (spn1->tail.spin.g[m2 - 1] != ind)
	++(m2);
      /*  second index position  */

      ll = spn1->tail.spin.l - 2;
      c = spn1->tail.spin.g[m2 - 1 - 1];
      spn1->tail.spin.l = ll;
      for (i = m2 - 1; i <= ll; i++)
	spn1->tail.spin.g[i - 1] = spn1->tail.spin.g[i + 2 - 1];
      spn1->tail.spin.g[ll + 1 - 1] = 0;
      spn1->tail.spin.g[ll + 2 - 1] = 0;

      nrevol = m2 - m1 - 1;
      if (nrevol == 0)
	{
	  multpolyint (&spn1->coef.complex.re, 4);
	  multpolyint (&spn1->coef.complex.im, 4);
	}
      else if ((nrevol & 1) == 0)	/*  even gamma case  */
	{
	  multpolyint (&(spn1->coef.complex.re), 2);
	  multpolyint (&(spn1->coef.complex.im), 2);

	  newmonom (&spn2);
	  spn2->next = spn1->next;
	  spn2->coef.complex.re = (poly) copypoly (spn1->coef.complex.re);
	  spn2->coef.complex.im = (poly) copypoly (spn1->coef.complex.im);
	  for (i = 0; i < spinLength; i++)
	    spn2->tail.power[i] = spn1->tail.power[i];
	  spn2->tail.spin.g[m1 - 1] = c;
	  spn1->next = spn2;

	  i = m1 + 1;
	  j = m2 - 3;
	  while (i < j)
	    {
	      cc = spn1->tail.spin.g[i - 1];
	      spn1->tail.spin.g[i - 1] = spn1->tail.spin.g[j - 1];
	      spn1->tail.spin.g[j - 1] = cc;
	      ++(i);
	      --(j);
	    }
	  spn1->tail.spin.g[m1 - 1] = spn1->tail.spin.g[m2 - 2 - 1];
	  spn1->tail.spin.g[m2 - 2 - 1] = c;
	  spn1 = spn2;		/*  ????  */
	}
      else
	{			/*  odd gamma case  */
	  i = m1 + 1;
	  j = m2 - 2;
	  while (i < j)
	    {
	      cc = spn1->tail.spin.g[i - 1];
	      spn1->tail.spin.g[i - 1] = spn1->tail.spin.g[j - 1];
	      spn1->tail.spin.g[j - 1] = cc;
	      ++(i);
	      --(j);
	    }
	  spn1->tail.spin.g[m1 - 1] = c;
	  multpolyint (&spn1->coef.complex.re, -2);
	  multpolyint (&spn1->coef.complex.im, -2);
/*            Spn2:=Spn1  */
	}
      spn1 = spn1->next;
    }				/*     Kahane  Begin     */
}				/*  Kahane  */

/* (T1_,T2_) */
static poly 
multspinstr (void)
{
  poly ans, ans2, pcoef, q;
  int b1, b2, m1, m2, ls1, ls2, llspin, gg5, sig, nrevol, i;
  int cind, dc, c;

  /*  MultSpinStr  */
  /*    MultSpinStr:=Nil; */

  b1 = 1;
  ls1 = t1_->tail.spin.l;
  b2 = 1;
  ls2 = t2_->tail.spin.l;

  llspin = ls1 + ls2;
  gg5 = (t1_->tail.spin.g5 + t2_->tail.spin.g5) & 1;

  if (forspur)
    {
      if ((llspin & 1) == 1)
	return NULL;
      if (gg5 == 1 && !levi)
	return NULL;
    }

  sig = t2_->tail.spin.g5 * (ls1 & 1);

  if (ls1 > 0 && ls2 > 0)
    {
      c = t1_->tail.spin.g[ls1 - 1];
      if (c < 0 && c == t2_->tail.spin.g[1 - 1])
	{
	  b2 = 2;
	  ls1--;
	  llspin -= 2;
	}
    }

  if (forspur && ls1 > 0 && ls2 >= b2)
    {
      c = t1_->tail.spin.g[0];
      if (c < 0 && c == t2_->tail.spin.g[ls2 - 1])
	{
	  b1 = 2;
	  ls2--;
	  llspin -= 2;
	  if (gg5)
	    sig = !sig;
	}
    }

  if (forspur && gg5 == 1 && llspin - 2 * strlen (indlist) < 4)
    return NULL;

  pcoef = (poly) plusone ();
  if (sig == 1)
    multpolyint (&pcoef, -1);

  if (b1 == 2)
    {
      q = pcoef;
      pcoef = (poly) multtwopoly (q, scalarmult (t1_->tail.spin.g[1 - 1],
						 t1_->tail.spin.g[1 - 1]));
      delpoly (&q);
    }

  if (b2 == 2)
    {
      q = pcoef;
      pcoef = (poly) multtwopoly (q, scalarmult (t2_->tail.spin.g[1 - 1],
						 t2_->tail.spin.g[1 - 1]));
      delpoly (&q);
    }
  if (pcoef == NULL)
    return NULL;

  newmonom (&ans);
  ans->next = NULL;
  ans->coef.complex.re = pcoef;
  ans->coef.complex.im = NULL;
  for (i = 0; i < spinLength; i++)
    ans->tail.power[i] = 0;
  ans->tail.spin.g5 = gg5;

  if (strlen (indlist) == 0)
    {
      m1 = ls1 + 1;
      m2 = b2 - 1;
      dc = 0;
    }
  else
    {
      cind = indlist[0];
      m1 = b1;
      while (t1_->tail.spin.g[m1 - 1] != cind)
	++(m1);
      m2 = b2;
      while (t2_->tail.spin.g[m2 - 1] != cind)
	++(m2);
      llspin -= 2;
      dc = -2;
    }

  c = 1 - b1;
  for (i = b1; i <= m1 - 1; i++)
    ans->tail.spin.g[i + c - 1] = t1_->tail.spin.g[i - 1];
  c += (1 - b2) + ls1 + dc;
  for (i = m2 + 1; i <= ls2; i++)
    ans->tail.spin.g[i + c - 1] = t2_->tail.spin.g[i - 1];
/*      for (i = llspin + 1; i <= 4 * spinlength - 2; i++)
   ans->tail.spin.g[i-1] = 0;
 */
  ans->tail.spin.l = llspin;

  if (strlen (indlist) == 0)
    return ans;

  nrevol = (ls1 - m1) + (m2 - b2);
  if (nrevol == 0)
    {
      multpolyint (&ans->coef.complex.re, 4);
      multpolyint (&ans->coef.complex.im, 4);
    }
  else
    {
      if ((nrevol & 1) == 0)
	{
	  if (m2 == b2)
	    {
	      cind = t1_->tail.spin.g[ls1 - 1];
	      --(ls1);
	    }
	  else
	    {
	      cind = t2_->tail.spin.g[m2 - 1 - 1];
	      --(m2);
	    }
	}
      c = m1 - b1 + m2;		/*  c-(m2-1)= m1+(1-b1)  */
      for (i = b2; i <= m2 - 1; i++)
	ans->tail.spin.g[c - i - 1] = t2_->tail.spin.g[i - 1];
      c = ls1 + m1 + m2 - b2 - b1 + 1;	/*  c-(m1+1)=(m1-b1)+(Ls1-m1)+(m2-b2)  */
      for (i = m1 + 1; i <= ls1; i++)
	ans->tail.spin.g[c - i - 1] = t1_->tail.spin.g[i - 1];

      if ((nrevol & 1) == 1)
	{
	  multpolyint (&ans->coef.complex.re, -2);
	  multpolyint (&ans->coef.complex.im, -2);
	}
      else
	{
	  ans->tail.spin.g[(ls1 - b1) + (m2 - b2) + 1 - 1] = cind;
	  multpolyint (&ans->coef.complex.re, 2);
	  multpolyint (&ans->coef.complex.im, 2);
	  newmonom (&ans2);
	  ans2->next = ans;
	  ans2->coef.complex.re = (poly) copypoly (ans->coef.complex.re);
	  ans2->coef.complex.im = (poly) copypoly (ans->coef.complex.im);
	  for (i = 0; i < spinLength; i++)
	    ans2->tail.power[i] = ans->tail.power[i];

	  ans2->tail.spin.g[m1 + (1 - b1) - 1] = cind;
	  c = ls1 + 2 - b1 - b2;	/*  c=Ls1+(1-b1)+(1-b2)  */
	  for (i = b2; i <= m2 - 1; i++)
	    ans2->tail.spin.g[c + i - 1] = t2_->tail.spin.g[i - 1];
	  c = 1 - b1;
	  for (i = m1 + 1; i <= ls1; i++)
	    ans2->tail.spin.g[c + i - 1] = t1_->tail.spin.g[i - 1];

	  ans = ans2;
	}
    }

  for (i = 2; i <= strlen (indlist); i++)
    kahane (&ans, indlist[i - 1]);

  ans2 = NULL;
  while (ans != NULL)		/*  Sorting  */
    {
      q = ans;
      ans = ans->next;
      q->next = NULL;
      sewtens (&ans2, &q, spinLength);
    }

  return ans2;
}				/*  MultSpinStr  */



/*  (T1,T2_:Poly) */

static poly 
multspinmemb (poly t1)
{
  poly ans, tmpres;

  /* Nested function: multspinstr */

  ans = NULL;
  t1_ = t1;
  do
    {
      /*  (T1_,T2_)  */
      tmpres = multspinstr ();
      multtensComplexpoly (&tmpres, t1_->coef.complex.re, t1_->coef.complex.im);
      sewtens (&ans, &tmpres, spinLength);
      t1_ = t1_->next;
    }
  while (t1_ != NULL);
  multtensComplexpoly (&ans, t2_->coef.complex.re, t2_->coef.complex.im);
  return ans;
}				/*  MultSpinMemb   */


poly 
multtwospin (poly t1, poly t2, int forspur1)
{
  char indlist2[119];
  poly ans, tmpres;
  int i;
  int ch;


  /* Nested function: multspinmemb */

  forspur = forspur1;
  if (t1 == NULL || t2 == NULL)
    return NULL;

  strcpy (indlist2, "");	/*  Seach   coopled  index  */
  for (i = 1; i <= t2->tail.spin.l; i++)
    if (t2->tail.spin.g[i - 1] > 0)
      sprintf (indlist2 + strlen (indlist2), "%c", t2->tail.spin.g[i - 1]);
  strcpy (indlist, "");
  for (i = 1; i <= t1->tail.spin.l; i++)
    if (t1->tail.spin.g[i - 1] > 0)
      {
	ch = t1->tail.spin.g[i - 1];
	if (strchr (indlist2, ch))
	  sprintf (indlist + strlen (indlist), "%c", ch);
      }

  ans = NULL;
  t2_ = t2;
  do
    {
      /*  ( T1,T2_ )  */
      tmpres = (poly) multspinmemb (t1);
      sewtens (&ans, &tmpres, spinLength);
      t2_ = t2_->next;
    }
  while (t2_ != NULL);
  return ans;
}				/*   MultTwoSpin   */


void 
multspintens (poly * spn, poly * tns)
{
  poly spn_, tns_, ans, sum1, sum2, qq;
  char indlist[61];
  int i, n, k;
  int j;

  ans = NULL;
  if (*tns != NULL && *spn != NULL)
    {
      tns_ = *tns;
      do
	{			/*  until Tms_=Nil  */
	  sum1 = NULL; /*   Sum1 will  Spn*Tns_ */ ;
	  strcpy (indlist, "");
	  for (i = 1; i <= maxIndex; i++)
	    if (tns_->tail.tens[i - 1] > i)
	      sprintf (indlist + strlen (indlist), "%c", i);
	  spn_ = *spn;
	  do
	    {
	      newmonom (&sum2);	/*  Sum2 will Snp_*Tns_/Tns_^.coef  */
	      memcpy (sum2->tail.tens, spn_->tail.tens, sizeof (long) * spinLength);
	      sum2->coef.complex.re = (poly) copypoly (spn_->coef.complex.re);
	      sum2->coef.complex.im = (poly) copypoly (spn_->coef.complex.im);
	      sum2->next = NULL;
	      for (i = 1; i <= maxIndex; i++)
		{
		  j = tns_->tail.tens[i - 1];
		  if (j < i && j != 0)
		    {
		      k = 1;
		      while (sum2->tail.spin.g[k - 1] != i)
			++(k);
		      sum2->tail.spin.g[k - 1] = j;
		      if (j < 0)
			{
			  if (k > 1 && j == sum2->tail.spin.g[k - 1 - 1])
			    {
			      for (n = k - 1; n <= sum2->tail.spin.l; n++)
				sum2->tail.spin.g[n - 1] = sum2->tail.spin.g[n + 2 - 1];
			      sum2->tail.spin.l -= 2;

			      qq = sum2->coef.complex.re;
			      sum2->coef.complex.re = (poly) multtwopoly (qq, scalarmult (j, j));
			      delpoly (&qq);

			      qq = sum2->coef.complex.im;
			      sum2->coef.complex.im = (poly) multtwopoly (qq, scalarmult (j, j));
			      delpoly (&qq);
			    }
			  else if (k < sum2->tail.spin.l &&
				   j == sum2->tail.spin.g[k + 1 - 1])
			    {
			      for (n = k; n <= sum2->tail.spin.l; n++)
				sum2->tail.spin.g[n - 1] = sum2->tail.spin.g[n + 2 - 1];
			      sum2->tail.spin.l -= 2;

			      qq = sum2->coef.complex.re;
			      sum2->coef.complex.re = (poly) multtwopoly (qq, scalarmult (j, j));
			      delpoly (&qq);

			      qq = sum2->coef.complex.im;
			      sum2->coef.complex.im = (poly) multtwopoly (qq, scalarmult (j, j));
			      delpoly (&qq);

			    }
			}
		    }
		}
	      if (sum2->coef.complex.re == NULL && sum2->coef.complex.im == NULL)
		delmonom (&sum2);
	      else
		sewtens (&sum1, &sum2, spinLength);
	      spn_ = spn_->next;
	    }
	  while (spn_ != NULL);
	  multtensComplexpoly (&sum1, tns_->coef.complex.re, tns_->coef.complex.im);
	  for (i = 1; i <= strlen (indlist); i++)
	    kahane (&sum1, indlist[i - 1]);
	  sewtens (&ans, &sum1, spinLength);
	  tns_ = tns_->next;
	}
      while (tns_ != NULL);
    }
  deltensor (spn);
  *spn = ans;
}				/* MultSpinTens  */
