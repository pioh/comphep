/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <limits.h>

#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"

#include "physics.h"
#include "out_service.h"
#include "pvars.h"
#include "procvar.h"
#include "optimise.h"

infoptr info = NULL;
int firstVar = 0;

static infoptr infoone;
/* static int ImConjKey = 0; */

void 
initinfo (void)
{
  info = (infoptr) getmem_ (sizeof (struct inforec));
  info->next = NULL;
  strcpy (info->name, "1");
  info->ival = 1;
  info->consttype = numb;
  infoone = info;
}

int 
equalexpr (varptr v1, varptr v2)
{
  while (v1 != NULL || v2 != NULL)
    {
      if (v1 == NULL || v2 == NULL)
	return FALSE;
      if (v1->sgn == v2->sgn && v1->coef == v2->coef &&
	  short_strcmp (v1->vars, v2->vars) == 0)
	{
	  v1 = v1->next;
	  v2 = v2->next;
	}
      else
	return FALSE;
    }
  return TRUE;
}

static void 
readmonom (FILE * file, short *varstr, short *conststr, NUM_TYPE * numc)
{
  int iv = 0, ic = 0;
  int deg, n, k, pos;

  fread (&readBuff->coef.num, readSize, 1, file);
  *numc = readBuff->coef.num;
  if (!*numc)
    return;

  for (n = 0; n < vardef->nvar; n++)
    {
      deg = (readBuff->tail.power[vardef->vars[n].wordpos - 1] /
	     vardef->vars[n].zerodeg) %
	vardef->vars[n].maxdeg;
/*----Slava test --------------------------
      if (deg)
       {
        if(ImConjKey==0) 
        readBuff->coef.num *= -1;   
        printf ("%d*%s",readBuff->coef.num, vararr[vardef->vars[n].num].alias); 
        if (deg > 1) 
          printf ("^%d\n", deg); 
       }   
------------------------------------------*/
      pos = vardef->vars[n].num;
      for (k = 0; k < deg; k++)
	{
	  if (pos >= firstVar)
	    varstr[iv++] = pos;
	  else
	    conststr[ic++] = pos;
	}


    }

  varstr[iv] = 0;
  conststr[ic] = 0;

/*  ImConjKey = 1; */
}				/* ReadMonom */

static void 
addnum (NUM_TYPE n, char *signum, infoptr * ans)
{
  if (n > 0)
    *signum = '+';
  else
    {
      *signum = '-';
      n = -n;
    }
  *ans = info;
  while (*ans != NULL)
    {
      if (((*ans)->consttype == numb) && ((*ans)->ival == n))
	return;
      *ans = (*ans)->next;
    }
  *ans = (infoptr) getmem_ (sizeof (struct inforec));
  (*ans)->next = info;
  (*ans)->consttype = numb;
  (*ans)->ival = n;
  info = *ans;
}

static int 
addtmpconst (varptr tmpconst, char *s, infoptr * coeff)
{
  varptr c;

  revers ((pointer *) & tmpconst);
  if (tmpconst->sgn == '-')
    {
      *s = '-';
      c = tmpconst;
      while (c != NULL)

	{
	  if (c->sgn == '-')
	    c->sgn = '+';
	  else
	    c->sgn = '-';
	  c = c->next;
	}
    }
  else
    *s = '+';
  if ((tmpconst->next == NULL) && (tmpconst->vars[0] == '\0'))
    {
      *coeff = tmpconst->coef;

      return FALSE;
    }

  *coeff = info;
  while (*coeff != NULL)
    if (((*coeff)->consttype == expr) && (equalexpr ((*coeff)->const_, tmpconst)))
      return FALSE;
    else
      *coeff = (*coeff)->next;
  *coeff = (infoptr) getmem_ (sizeof (struct inforec));
  (*coeff)->next = info;
  (*coeff)->const_ = tmpconst;
  (*coeff)->consttype = expr;
  info = *coeff;
  return TRUE;
}				/*  AddTmpConst */


void 
readpolynom (FILE * fres, varptr * expr_)
{
  short varstr[STRSIZ], conststr[STRSIZ];
  NUM_TYPE n;
  pointer pntr;
  varptr tmpconst;
  char s;
  infoptr coeff;
  marktp tmpmark;

  readmonom (fres, varstr, conststr, &n);
  if (!n)
    {
      *expr_ = NULL;
      return;
    }

  *expr_ = (varptr) getmem_ (minvarrec + sizeof (short) * short_strlen (varstr));
  (*expr_)->next = NULL;
  short_strcpy ((*expr_)->vars, varstr);

  addnum (n, &s, &coeff);
  mark_ (&tmpmark);
  tmpconst = (varptr) getmem_ (minvarrec + sizeof (short) * short_strlen (conststr));
  tmpconst->next = NULL;
  short_strcpy (tmpconst->vars, conststr);
  tmpconst->sgn = s;
  tmpconst->coef = coeff;

  while (1)
    {
      readmonom (fres, varstr, conststr, &n);
      if (!n)
	break;
      if (short_strcmp (varstr, (*expr_)->vars) != 0)
	{
	  if (!addtmpconst (tmpconst, &((*expr_)->sgn), &((*expr_)->coef)))
	    release_ (&tmpmark);
	  pntr = (pointer) (*expr_);
	  *expr_ = (varptr) getmem_ (minvarrec + sizeof (short) * short_strlen (varstr));
	  (*expr_)->next = (varptr) pntr;
	  short_strcpy ((*expr_)->vars, varstr);
	  pntr = NULL;
	  addnum (n, &s, &coeff);
	  mark_ (&tmpmark);
	}
      else
	{
	  pntr = (pointer) tmpconst;
	  addnum (n, &s, &coeff);
	}
      tmpconst = (varptr) getmem_ (minvarrec + sizeof (short) * short_strlen (conststr));
      tmpconst->next = (varptr) pntr;
      short_strcpy (tmpconst->vars, conststr);
      tmpconst->sgn = s;
      tmpconst->coef = coeff;
    }
  if (!addtmpconst (tmpconst, &((*expr_)->sgn), &((*expr_)->coef)))
    release_ (&tmpmark);
}				/* ReadPolynom */


static void 
findmaxvar (varptr ex, unsigned *n, short *ch, int *power)
{
  int *nterms;
  int *minpower;
  int k, bt, d, nv;

  nterms = (int *) m_alloc (sizeof (int) * nProcessVar);
  minpower = (int *) m_alloc (sizeof (int) * nProcessVar);

  for (k = firstVar; k < nProcessVar; k++)
    {
      nterms[k] = 0;
      minpower[k] = 0;
    }

  while (ex != NULL)
    {
      if (short_strlen (ex->vars))
	{
	  d = 1;
	  bt = ex->vars[0];
	  for (k = 1; k < short_strlen (ex->vars); k++)
	    {
	      nv = ex->vars[k];
	      if (nv != bt)
		{
		  minpower[bt] = minpower[bt] ? MIN (d, minpower[bt]) : d;
		  nterms[bt]++;
		  d = 1;
		  bt = nv;
		}
	      else
		++(d);
	    }
	  minpower[bt] = minpower[bt] ? MIN (d, minpower[bt]) : d;
	  nterms[bt]++;
	}
      ex = ex->next;
    }

  bt = firstVar;
  *n = nterms[firstVar];
  *power = minpower[firstVar];

  for (k = firstVar + 1; k < nProcessVar; k++)
    if (*n < nterms[k] || (*n == nterms[k] && *power < minpower[k]))
      {
	*n = nterms[k];
	bt = k;
	*power = minpower[k];
      }
  *ch = bt;
  free (nterms);
  free (minpower);
}


static void 
findmaxcoef (varptr ex, infoptr * i, unsigned *n)
{
  varptr jj;
  jj = ex;
  while (jj != NULL)
    {
      (jj->coef)->count = 0;
      jj = jj->next;
    }
  jj = ex;
  while (jj != NULL)
    {
      (jj->coef)->count++;
      jj = jj->next;
    }

  *n = 0;
  infoone = info;
  while (infoone->next != NULL)
    infoone = infoone->next;
  (*i) = infoone;

  (*i)->count = 0;
  jj = ex;
  while (jj != NULL)
    {
      if (*n < (jj->coef)->count)

	{
	  *i = jj->coef;
	  *n = (*i)->count;
	}
      jj = jj->next;
    }
}

static void 
clipvar (varptr ex, short ch, int power, varptr * ex1, varptr * ex2)
{
  varptr exnext;

  var_rec ex1rec;
  var_rec ex2rec;
  varptr ex1_, ex2_;
  short *u;

  ex1_ = &ex1rec;
  ex2_ = &ex2rec;
  while (ex != NULL)
    {
      u = short_strchr (ex->vars, ch);
      exnext = ex->next;
      if (u)
	{
	  int i = 0;
	  do
	    u[i] = u[i + power];
	  while (u[i++]);
	  ex1_->next = ex;
	  ex1_ = ex;
	}
      else
	{
	  ex2_->next = ex;
	  ex2_ = ex;
	}
      ex = exnext;
    }


  ex1_->next = NULL;
  *ex1 = ex1rec.next;
  ex2_->next = NULL;
  *ex2 = ex2rec.next;
}


static void 
clipconst (varptr ex, infoptr i, varptr * ex1, varptr * ex2)
{

  varptr exnext;
  infoptr one;

  var_rec ex1rec, ex2rec;
  varptr ex1_, ex2_;

  one = info;
  while (one->next != NULL)
    one = one->next;

  ex1_ = &ex1rec;
  ex2_ = &ex2rec;

  while (ex != NULL)
    {
      exnext = ex->next;
      if (ex->coef != i)
	{
	  ex2_->next = ex;
	  ex2_ = ex;
	}
      else
	{
	  ex->coef = one;
	  ex1_->next = ex;
	  ex1_ = ex;
	}
      ex = exnext;
    }


  ex1_->next = NULL;
  *ex1 = ex1rec.next;
  ex2_->next = NULL;
  *ex2 = ex2rec.next;

}

pointer 
emitexpr (varptr ex, smplemit smplemitfun, vfact vfactfun,
	  cfact cfactfun)
{
  unsigned nv, nc;
  short ch;
  infoptr i;
  varptr ex1, ex2;
  pointer pmult, psum;
  int deg;

  findmaxcoef (ex, &i, &nc);
  findmaxvar (ex, &nv, &ch, &deg);
  if (nc < 2 && nv < 2)
    return smplemitfun (ex);
  if (nv >= nc)
    {
      clipvar (ex, ch, deg, &ex1, &ex2);
      pmult = emitexpr (ex1, smplemitfun, vfactfun, cfactfun);
      psum = emitexpr (ex2, smplemitfun, vfactfun, cfactfun);
      return vfactfun (ch, deg, pmult, psum);
    }
  else
    {
      clipconst (ex, i, &ex1, &ex2);
      pmult = emitexpr (ex1, smplemitfun, vfactfun, cfactfun);
      psum = emitexpr (ex2, smplemitfun, vfactfun, cfactfun);
      return cfactfun (i, pmult, psum);
    }
}

int 
short_strlen (short *s)
{
  int i = 0;
  while (s[i])
    i++;
  return i;
}

void 
short_strcpy (short *to, short *from)
{
  int i = 0;

  while (from[i])
    {
      to[i] = from[i];
      i++;
    }
  to[i] = 0;
}

int 
short_strcmp (short *s1, short *s2)
{
  int i = 0;
  while (s1[i] && s2[i] && s1[i] == s2[i])
    i++;
  return s1[i] - s2[i];
}


short *
short_strchr (short *str, short s)
{
  int i = 0;
  while (str[i] && str[i] != s)
    i++;
  if (str[i] == s)
    return (str + i);
  else
    return NULL;
}
