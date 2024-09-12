/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/parser.h"
#include "service2/include/syst.h"
#include "chep_crt/include/chep_crt.h"
#include "polynom/include/polynom.h"
#include "polynom/include/tensor.h"

#include "pvars.h"
#include "pre_read.h"
#include "test_wrt.h"

#define xmax 120
#define fle stdout

static int xpos = 1;

static void 
wrtln (void)
{
  fprintf (stdout, "\n");
  xpos = 1;
}


static void 
wrt (char *s)
{
  unsigned l;

  l = strlen (s);
  if (xpos + l <= xmax)
    {
      fprintf (stdout, "%s", s);
      xpos += l;
    }
  else
    {
      l = xmax - xpos;
      while (l > 0 && strchr ("*+-)(^ ", s[l - 1]) == NULL)
	--l;
      if (l == 0)
	{
	  if (xpos == 1)
	    l = xmax;
	  else
	    {
	      wrtln ();
	      wrt (s);
	      return;
	    }
	}
      wrt (copy (s, 1, l));
      wrtln ();
      wrt (copy (s, l + 1, strlen (s) - l));
    }
}


void 
tracePrn (char *format,...)
{
  va_list args;
  char dump[STRSIZ], *beg, *nn;

  va_start (args, format);
  vsprintf (dump, format, args);
  va_end (args);

  beg = dump;
  while (1)
    {
      nn = strchr (beg, '\n');
      if (nn == NULL)
	{
	  wrt (beg);
	  return;
	}
      nn[0] = 0;
      wrt (beg);
      wrtln ();
      beg = nn + 1;
    }
}


static void 
writepoly (polyvars * v, poly p)
{
  char txt[STRSIZ], numtxt[STRSIZ];
  int i, deg;
  int beg, first;
  unsigned long wpower;

  if (!p)
    {
      wrt ("0");
      return;
    }
  beg = 1;

  for (beg = 1; p; p = p->next)
    {
      strcpy (txt, "");
      first = 1;
      wpower = p->tail.power[0];
      for (i = 0; i < v->nvar; i++)
	{
	  deg = (p->tail.power[v->vars[i].wordpos - 1] / v->vars[i].zerodeg) % v->vars[i].maxdeg;

	  if (deg > 0)
	    {
	      if (first)
		first = 0;
	      else
		strcat (txt, "*");
	      strcat (txt, v->vars[i].name);
	      if (deg > 1)
		sprintf (txt + strlen (txt), "^%d", deg);
	    }
	}

      if (strlen (txt) != 0)
	{
	  if (p->coef.num == NUM_ONE || p->coef.num == -NUM_ONE)
	    numtxt[0] = 0;
	  else
	    sprintf (numtxt, "%" NUM_STR "*", p->coef.num > 0 ? p->coef.num : -p->coef.num);
	}
      else
	sprintf (numtxt, "%" NUM_STR, p->coef.num > 0 ? p->coef.num : -p->coef.num);

      if (p->coef.num < NUM_ZERO)
	wrt (" - ");
      if (beg)
	beg = 0;
      else if (p->coef.num > NUM_ZERO)
	wrt (" + ");

      wrt (numtxt);
      wrt (txt);
    }
}


void 
writetens (polyvars * vars, poly p)
{
  char txt[STRSIZ];
  int i, s, l;
  int beg, first;

  if (!p)
    {
      wrt ("0");
      return;
    }

  for (beg = 1; p; p = p->next)
    {
      strcpy (txt, "");
      first = 1;
      for (i = 1; i <= maxIndex; i++)
	{
	  s = p->tail.tens[i - 1];
	  if (s < i && s)
	    {
	      if (first)
		first = 0;
	      else
		strcat (txt, "*");
	      if (s > 0)
		sprintf (txt + strlen (txt), "m%d.m%d", i, s);
	      else
		sprintf (txt + strlen (txt), "p%d.m%d", -s, i);
	    }
	}
      if (levi && p->tail.tens[maxIndex] != 0)
	{
	  if (first)
	    first = 0;
	  else
	    strcat (txt, "*");
	  strcat (txt, "eps(");
	  for (l = 0; l < 4; l++)
	    {
	      s = p->tail.tens[maxIndex + l];
	      if (s < 0)
		sprintf (txt + strlen (txt), "p%d,", -s);
	      else
		sprintf (txt + strlen (txt), "m%d,", s);
	    }
	  txt[strlen (txt) - 1] = ')';
	}
      if (beg)
	beg = 0;
      else
	wrt (" +");
      wrt ("{ ");
      if (p->coef.complex.re)
	writepoly (vars, p->coef.complex.re);
      if (p->coef.complex.im)
	{
	  wrt ("+i*(");
	  writepoly (vars, p->coef.complex.im);
	  wrt (")");
	}
      wrt (" }");

      if (strlen (txt) != 0)
	wrt ("*");
      wrt (txt);
    }
}


void 
writespinor (polyvars * vars, poly p)
{
  char txt[STRSIZ], num1[STRSIZ];
  int i, s, ls;
  int beg;

  if (!p)
    {
      wrt ("0");
      return;
    }
  for (beg = 1; p; p = p->next)
    {
      ls = p->tail.spin.l;
      if (p->tail.spin.g5)
	if (!ls)
	  strcpy (txt, "g(ln,a)");
	else
	  strcpy (txt, "g(ln,a,");
      else if (ls)
	strcpy (txt, "g(ln,");
      else
	strcpy (txt, "");
      for (i = 1; i <= ls; i++)
	{
	  s = p->tail.spin.g[i - 1];
	  if (s < 0)
	    {
	      sprintf (num1, "%d", -s);
	      strcat (txt, "p");
	    }
	  else
	    {
	      sprintf (num1, "%d", s);
	      strcat (txt, "m");
	    }
	  strcat (txt, num1);
	  if (i != ls)
	    strcat (txt, ",");
	}
      if (ls)
	strcat (txt, ")");
      if (beg)
	beg = 0;
      else
	wrt (" + ");
      wrt ("{ ");
      if (p->coef.complex.re)
	writepoly (vars, p->coef.complex.re);
      if (p->coef.complex.im)
	{
	  wrt ("+i*(");
	  writepoly (vars, p->coef.complex.im);
	  wrt (")");
	}
      wrt (" }");

      if (strlen (txt))
	{
	  wrt ("*");
	  wrt (txt);
	}
    }
}

void 
writeexpression (polyvars * vars, poly m)
{
  char s[5];
  int n;
  poly p;

  if (m->coef.type <= tenstp)
    writetens (vars, m->next);
  else if (m->coef.type == spintp)
    writespinor (vars, m->next);
  else if (m->coef.type == vectortp)
    {
      p = m->next;
      if (p == NULL)
	wrt ("0");
      else
	do
	  {
	    n = -p->tail.tens[0];
	    sprintf (s, "p%d*(", n);
	    wrt (s);
	    if (p->coef.complex.re)
	      writepoly (vars, p->coef.complex.re);
	    if (p->coef.complex.im)
	      {
		wrt ("i*(");
		writepoly (vars, p->coef.complex.re);
		wrt (")");
	      }
	    wrt (")");
	    p = p->next;
	    if (p)
	      wrt ("+");
	  }
	while (p != NULL);
    }
  else
    {
      n = m->next->tail.tens[0];
      if (!n)
	n = 1;
      sprintf (s, "l%d", n);
      wrt (s);
    }
}

void 
writevars (polyvars * v)
{
  int i;
  char txt[STRSIZ];
  for (i = 0; i < v->nvar; i++) {
    printf ("var %i: name = %s, num = %i, wordpos = %i, zerodeg = %li, maxdeg = %li\n", i + 1,
    v->vars[i].name, 
    v->vars[i].num, 
    v->vars[i].wordpos, 
    v->vars[i].zerodeg, 
    v->vars[i].maxdeg);
    wrt (txt);
  }
}
