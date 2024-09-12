/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <limits.h>
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/getmem.h"
#include "service2/include/syst.h"
#include "polynom/include/polynom.h"

#include "physics.h"
#include "sos.h"
#include "pvars.h"

#ifdef STRACE
#include "test_wrt.h"
#endif

polyvars *vardef = NULL;

void 
increaseVars (polyvars * v)
{
  int oldsize = ALIG (v->nvar);
  if (++v->nvar > oldsize)
    v->vars = re_alloc (v->vars, ALIG (v->nvar) * sizeof (*v->vars));
}

void 
clearVars (polyvars * v)
{
  if (v->vars)
    {
      free (v->vars);
      v->vars = NULL;
    }
  v->nvar = 0;
}

void 
unite_vars (polyvars * result, polyvars * addition)
{
  int n, nn;
  char s[STRSIZ];

  for (nn = 0; nn < addition->nvar; nn++)
    {
      strcpy (s, addition->vars[nn].name);

      n = 0;
      while (n < result->nvar &&
	     strcmp (addition->vars[nn].name, result->vars[n].name))
	n++;
      if (n < result->nvar)
	result->vars[n].maxdeg =
	  MAX (result->vars[n].maxdeg, addition->vars[nn].maxdeg);
      else
	{
	  increaseVars (result);
	  strcpy (result->vars[n].name, s);
	  result->vars[n].maxdeg = addition->vars[nn].maxdeg;
	}
    }
}


void 
addvar (polyvars * v, char *varname, int deg)
{
  int n;

  n = 0;
  while (n < v->nvar && strcmp (varname, v->vars[n].name))
    n++;
  if (n < v->nvar)
    v->vars[n].maxdeg += deg;
  else
    {
      increaseVars (v);
      strcpy (v->vars[n].name, varname);
      v->vars[n].maxdeg = deg + 1;
    }
}


int 
scalarProductPos (int p1, int p2)
{
  if (p1 > p2)
    {
      int pp = p1;
      p1 = p2;
      p2 = pp;
    }
  return nmodelvar + 1 + p1 + ((p2 - 1) * (p2 - 2)) / 2;
}



int 
modelVarPos (char *s)
{
  int bt;
  for (bt = 0; bt <= nmodelvar; bt++)
    {
      if (strcmp (modelvars[bt].varname, s) == 0)
	return bt;
    }
  save_sos (14);
  fprintf (stderr, "***** modelVarPos: Reach the end. No return int val.!?\n");
  return 0;
}


void 
closevars (polyvars * v)
{
  int i;
  unsigned long z;
  int wp, p1, p2;


  for (i = 0; i < v->nvar; i++)	/*  numeration  */
    if (cpos ('.', v->vars[i].name) != 0)
      {
	p1 = v->vars[i].name[1] - '0';
	p2 = v->vars[i].name[4] - '0';
	v->vars[i].num = scalarProductPos (p1, p2);
      }
    else
      v->vars[i].num = modelVarPos (v->vars[i].name);

  for (i = 0; i < v->nvar; i++)
    if (!strcmp (v->vars[i].name, "i"))
      {
	v->nvar--;
	break;
      }
  for (; i < v->nvar; i++)
    v->vars[i] = v->vars[i + 1];


  if (v->nvar > 1)		/*  Sorting  */
    {
      i = 1;
      while (i < v->nvar)
	if (v->vars[i - 1].num > v->vars[i].num)
	  i++;
	else
	  {
	    varinfo tmpv = v->vars[i - 1];
	    v->vars[i - 1] = v->vars[i];
	    v->vars[i] = tmpv;
	    if (i == 1)
	      ++i;
	    else
	      --i;
	  }
    }

  monomLength = 1;
  z = 1;
  for (i = 1; i <= v->nvar; i++)
    {
      if (z >= (ULONG_MAX / v->vars[i - 1].maxdeg))
	{
	  monomLength++;
	  z = v->vars[i - 1].maxdeg;
	}
      else
	z *= v->vars[i - 1].maxdeg;

      v->vars[i - 1].wordpos = monomLength;
    }

  z = 1;
  wp = monomLength;
  for (i = v->nvar; i >= 1; i--)
    {
      if (v->vars[i - 1].wordpos == wp)
	{
	  v->vars[i - 1].zerodeg = z;
	  z *= v->vars[i - 1].maxdeg;
	}
      else
	{
	  v->vars[i - 1].zerodeg = 1;
	  z = v->vars[i - 1].maxdeg;
	  --wp;
	}
    }

  set_garbage (NULL);
#ifdef STRACE
  tracePrn ("vars position \n");
  for (i = 0; i < v->nvar; i++)
    {
      tracePrn (" name= %s pos= %d maxdeg=  %lu zerodeg= %lu\n",
      v->vars[i].name, v->vars[i].wordpos, v->vars[i].maxdeg, v->vars[i].zerodeg);
    }
  tracePrn ("monomLength=  %d\n", monomLength);
#endif

}
