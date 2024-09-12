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

int maxIndex = 0;
int tensLength = 0;
int dellevi = 0;

static int sw;			/* from multvectors */
static char ep1[4];		/* from multstructure */
static poly tres;		/* from multstructure */
static poly t1_;		/* from mult_t1_memb_t2 */
static poly t2_;		/* from multtwotens */


static void 
jump (char *s)
{
label_1:
  *s = sw ? t1_->tail.tens[*s - 1] : t2_->tail.tens[*s - 1];
  if (*s < 0)
    return;
  if (tres->tail.tens[*s - 1] == 123)
    {
      tres->tail.tens[*s - 1] = 0;
      sw = !sw;
      goto label_1;
    }
}				/*   Jump   */


static void 
multvectors (void)
{
  poly p, q;
  int i, nloop;
  char s1, s2;

  newmonom (&tres);

  for (i = 0; i < tensLength; i++)
    tres->tail.power[i] = 0;
  for (i = 1; i <= maxIndex; i++)
    tres->tail.tens[i - 1] = t1_->tail.tens[i - 1] == 0 ?
      t2_->tail.tens[i - 1] == 0 ? 0 : 122 :
      t2_->tail.tens[i - 1] == 0 ? 121 : 123;
  nloop = 0;
  p = (poly) plusone ();
  for (i = 1; i <= maxIndex; i++)
    if (tres->tail.tens[i - 1] == 123)
      {
	tres->tail.tens[i - 1] = 0;
	sw = TRUE;
	s1 = i;
	jump (&s1);
	if (s1 == i)
	  ++(nloop);
	else
	  {
	    sw = FALSE;
	    s2 = i;
	    jump (&s2);
	    if (s1 > 0)
	      if (s2 > 0)
		{
		  tres->tail.tens[s1 - 1] = s2;
		  tres->tail.tens[s2 - 1] = s1;
		}
	      else
		tres->tail.tens[s1 - 1] = s2;
	    else if (s2 > 0)
	      tres->tail.tens[s2 - 1] = s1;
	    else
	      {
		q = (poly) multtwopoly (p, scalarmult (s1, s2));
		delpoly (&p);
		p = q;
	      }
	  }
      }
  if (p == NULL)
    delmonom (&tres);
  else
    {
      for (i = 1; i <= maxIndex; i++)
	if (tres->tail.tens[i - 1] == 121)
	  tres->tail.tens[i - 1] = t1_->tail.tens[i - 1];
	else if (tres->tail.tens[i - 1] == 122)
	  tres->tail.tens[i - 1] = t2_->tail.tens[i - 1];
      multpolyint (&p, 1 << (2 * nloop));	/*  loops  */
      tres->coef.complex.re = p;
      tres->coef.complex.im = NULL;
      tres->next = NULL;
    }
}				/*  MultVectors  */

static void 
correcteps (void)
{
  int i;
  int c, cc, sign;

  for (i = maxIndex + 1; i <= maxIndex + 4; i++)
    {
      c = tres->tail.tens[i - 1];
      if (c > 0)
	{
	  if (tres->tail.tens[c - 1] < 0)
	    {
	      tres->tail.tens[i - 1] = tres->tail.tens[c - 1];
	      tres->tail.tens[c - 1] = 0;
	    }
	  else if (tres->tail.tens[c - 1] > 0)
	    {
	      cc = tres->tail.tens[c - 1];
	      tres->tail.tens[i - 1] = cc;
	      tres->tail.tens[c - 1] = 0;
	      tres->tail.tens[cc - 1] = 0;
	    }
	}			/* if c>0 */
    }				/*  For  */
  i = maxIndex + 1;
  sign = 1;
  do
    {
      c = tres->tail.tens[i - 1];
      cc = tres->tail.tens[i + 1 - 1];
      if (c > cc)
	++(i);
      else if (c == cc)
	{
	  deltensor (&tres);
	  return;
	}
      else
	{
	  tres->tail.tens[i - 1] = cc;
	  tres->tail.tens[i + 1 - 1] = c;
	  if (i == maxIndex + 1)
	    ++(i);
	  else
	    --(i);
	  sign = -sign;
	}
    }
  while (i != maxIndex + 4);
  if (sign == -1)
    {
      multpolyint (&tres->coef.complex.re, -1);
      multpolyint (&tres->coef.complex.im, -1);
    }
}				/*   CorrectEps  */


static void 
multeps (void)			/*  Tres:=Tres*Eps1  */
{
  static struct epsdata
  {
    int contr[4];
    int sgn;
  }
  epsdata[24]
  =
  {
    {
      {
	1, 2, 3, 4
      }
      ,1
    }
    ,
    {
      {
	2, 1, 3, 4
      }
      ,-1
    }
    ,
    {
      {
	2, 3, 1, 4
      }
      ,1
    }
    ,
    {
      {
	3, 2, 1, 4
      }
      ,-1
    }
    ,
    {
      {
	3, 1, 2, 4
      }
      ,1
    }
    ,
    {
      {
	1, 3, 2, 4
      }
      ,-1
    }
    ,
    {
      {
	1, 2, 4, 3
      }
      ,-1
    }
    ,
    {
      {
	2, 1, 4, 3
      }
      ,1
    }
    ,
    {
      {
	2, 3, 4, 1
      }
      ,-1
    }
    ,
    {
      {
	3, 2, 4, 1
      }
      ,1
    }
    ,
    {
      {
	3, 1, 4, 2
      }
      ,-1
    }
    ,
    {
      {
	1, 3, 4, 2
      }
      ,1
    }
    ,
    {
      {
	1, 4, 3, 2
      }
      ,-1
    }
    ,
    {
      {
	2, 4, 3, 1
      }
      ,1
    }
    ,
    {
      {
	2, 4, 1, 3
      }
      ,-1
    }
    ,
    {
      {
	3, 4, 1, 2
      }
      ,1
    }
    ,
    {
      {
	3, 4, 2, 1
      }
      ,-1
    }
    ,
    {
      {
	1, 4, 2, 3
      }
      ,1
    }
    ,
    {
      {
	4, 2, 1, 3
      }
      ,1
    }
    ,
    {
      {
	4, 1, 2, 3
      }
      ,-1
    }
    ,
    {
      {
	4, 3, 2, 1
      }
      ,1
    }
    ,
    {
      {
	4, 2, 3, 1
      }
      ,-1
    }
    ,
    {
      {
	4, 1, 3, 2
      }
      ,1
    }
    ,
    {
      {
	4, 3, 1, 2
      }
      ,-1
    }
  };
  static int epsnumb[5] =
  {1, 1, 2, 6, 24};
  static int epsfact[5] =
  {24, 6, 2, 1, 1};

  char ep2[4];
  char ep0[4];
  int i, j, i1, j1, l;
  int c, cc;
  poly tensarr[24];

  for (i = 0; i < 4; i++)
    ep0[i] = 0;
  memcpy (ep2, &(tres->tail.tens[maxIndex]), 4);
  memcpy (&(tres->tail.tens[maxIndex]), ep0, 4);

  l = 4;
  i = 1;
  j = 1;
  c = 1;
  do
    {				/*  delete caupled index  */
      if (ep1[i - 1] > ep2[j - 1])
	++(i);
      else if (ep1[i - 1] < ep2[j - 1])
	++(j);
      else if (ep1[i - 1] > 0)
	{
	  --(l);
	  for (i1 = i; i1 <= l; i1++)
	    ep1[i1 - 1] = ep1[i1];
	  for (j1 = j; j1 <= l; j1++)
	    ep2[j1 - 1] = ep2[j1];
	  if (((i - j) & 1) == 1)
	    c = -c;
	}
    }
  while (!(i > l || j > l || ep1[i - 1] < 0 || ep2[j - 1] < 0));
  if (c < 0)
    multtensint (&tres, -1);
  tensarr[0] = tres;
  for (i = 2; i <= epsnumb[l]; i++)
    tensarr[i - 1] = (poly) copytens (tres, tensLength);
  tres = NULL;
  for (i = 1; i <= epsnumb[l]; i++)
    {
      for (j = 1; j <= l; j++)
	{
	  if (!tensarr[i - 1])
	    break;
	  c = ep1[j - 1];
	  cc = ep2[epsdata[i - 1].contr[j - 1] - 1];
	  if (c > 0)
	    tensarr[i - 1]->tail.tens[c - 1] = cc;
	  if (cc > 0)
	    tensarr[i - 1]->tail.tens[cc - 1] = c;
	  else if (c < 0)
	    multtenspoly (&tensarr[i - 1], scalarmult (c, cc));
	}
      multtensint (&tensarr[i - 1], -epsdata[i - 1].sgn * epsfact[l]);
      sewtens (&tres, &tensarr[i - 1], tensLength);
    }
}				/*  MultEps  */


static poly 
multstructure (void)		/*  T1_,T2_ are extarnal arguments  */
{
  int pseudo1, pseudo2;

  multvectors ();		/*  Tres:=T1_,T2_  Pseudovector parts not done   */
  if (tres == NULL)
    return NULL;
  if (levi)
    {
      pseudo1 = (t1_->tail.tens[maxIndex] != 0);
      pseudo2 = (t2_->tail.tens[maxIndex] != 0);
      if (dellevi && (pseudo1 ^ pseudo2))
	return NULL;
      if (pseudo1)
	{
	  memcpy (&(tres->tail.tens[maxIndex]), &(t1_->tail.tens[maxIndex]), 4);
/*         tres->tail.power[tenslength-1] = t1_->tail.power[tenslength-1]; */
	  correcteps ();	/*  in Tres  */
	  if (tres == NULL)
	    return NULL;
	}
      if (pseudo2)
	{
	  if (pseudo1)
	    memcpy (ep1, &(tres->tail.tens[maxIndex]), 4);
	  memcpy (&(tres->tail.tens[maxIndex]), &(t2_->tail.tens[maxIndex]), 4);
/*            for (i = 1; i <= 4; i++) 
   ep1[i-1] = tres->tail.tens[maxIndex + i-1]; 
   for (k=maxIndex;k<maxIndex+4;k++)  tres->tail.power[k] = t2_->tail.power[k];
 */
	  correcteps ();	/*  in Tres  */
	  if (tres == NULL)
	    return NULL;
	}

      if (pseudo1 && pseudo2)
	multeps ();		/*  Tres:=Tres*Ep1  */
/*      if (!(pseudo1 || pseudo2)) 
   {
   tres->tail.power[tenslength-1] = 0;
   } 
 */
    }				/*  If epsCase  */
  return tres;
}				/*  MultStructure  */


static poly 
mult_t1_memb_t2 (poly t1)	/*  T2 external ; T2 is const  */
{
  poly ans_, tmpres_;


  /* Nested function: multstructure */

  t1_ = t1;
  ans_ = NULL;
  do
    {
      tmpres_ = (poly) multstructure ();	/*  (T1_,T2_)  */
      if (tmpres_ != NULL)
	{
	  multtensComplexpoly (&tmpres_, t1_->coef.complex.re, t1_->coef.complex.im);
	  sewtens (&ans_, &tmpres_, tensLength);
	}
      t1_ = t1_->next;
    }
  while (t1_ != NULL);
  multtensComplexpoly (&ans_, t2_->coef.complex.re, t2_->coef.complex.im);
  return ans_;
}				/*  Mult_T1_Memb_T2  */


poly 
multtwotens (poly t1, poly t2)
{
  poly ans, tmpres;
  int swap;
  poly t_1, t_2;

  t_1 = t1;
  t_2 = t2;

  swap = FALSE;
  while (t_1 != NULL && t_2 != NULL)
    {
      t_1 = t_1->next;
      t_2 = t_2->next;
      if (t_2 == NULL)
	swap = TRUE;
    }

  if (swap)
    {
      t_1 = t1;
      t1 = t2;
      t2 = t_1;
    }

  ans = NULL;
  t2_ = t2;
  if (t1 != NULL && t2 != NULL)
    do
      {
	tmpres = (poly) mult_t1_memb_t2 (t1);	/*  T2 external  */
	sewtens (&ans, &tmpres, tensLength);
	t_2 = t2_;
	t2_ = t2_->next;
	delpoly (&t_2->coef.complex.re);
	delpoly (&t_2->coef.complex.im);
	delmonom (&t_2);
      }
    while (t2_ != NULL);
  else
    deltensor (&t2);

  deltensor (&t1);
  t2 = NULL;
  return ans;
}				/*   MultTwoTens   */
