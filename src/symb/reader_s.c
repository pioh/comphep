/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/parser.h"
#include "service2/include/syst.h"
#include "polynom/include/polynom.h"
#include "polynom/include/tensor.h"
#include "polynom/include/spinor.h"

#include "pvars.h"
#include "sos.h"
#include "reader0.h"
#include "reader_s.h"

#ifdef STRACE
#include "test_wrt.h"
#endif

int r_reading2 = FALSE;


void *
rd_symb (char *s)
{
  int i;
  poly p, ans;

  newmonom (&p);

  ans = plusone ();
  p->next = ans;
  ans->coef.complex.re = plusone ();
  ans->coef.complex.im = NULL;


  if (isdigit (s[0]))
    {
      long l;
      sscanf (s, "%ld", &l);
      multtensint (&ans, l);
      p->next = ans;
      p->coef.type = numbertp;
    }
  else
    {
      p->coef.type = polytp;
      if (strlen (s) == 2 && isdigit (s[1]) && s[1] != 0)
	{
	  i = s[1] - '0';
	  switch (s[0])
	    {
	    case 'p':
	    case 'P':
	      p->coef.type = vectortp;
	      ans->tail.tens[0] = -abs (momsubst[i - 1]);
	      if (momsubst[i - 1] < 0)
		(ans->coef.complex.re)->coef.num = -NUM_ONE;
	      break;

	    case 'm':
	    case 'M':
	      p->coef.type = indextp;
	      ans->tail.tens[0] = indsubst[i - 1];
	      if (s[0] == 'M')
		ans->tail.tens[0]--;
	      if (ans->tail.tens[0] == 1)
		ans->tail.tens[0] = 0;
	      else
		ans->tail.tens[ans->tail.tens[0] - 1] = 1;
	    }			/*  Case  */
	  if (strcmp (s, "G5") == 0)
	    {
	      p->coef.type = spintp;
	      ans->tail.spin.g5 = 1;
	    }
	}

      if (p->coef.type == polytp)
	{   
	  if ( strcmp (s, "i") == 0 ) 
	    {
	      ans->coef.complex.im = ans->coef.complex.re;
	      ans->coef.complex.re = NULL;
	    }  
	  else
	    {
	      i = 1;
	      while (strcmp (vardef->vars[i - 1].name, s) != 0)
		if (++i > vardef->nvar)
		  save_sos (14);
	      ans->coef.complex.re->tail.power[vardef->vars[i - 1].wordpos - 1] =
		vardef->vars[i - 1].zerodeg;
	    }
	}
    }
#ifdef STRACE
  tracePrn (" \n rd_ (%s) -> ", s);
  writeexpression (vardef, p);
#endif
  return (void *) p;
}


static void *
uact_symb (char *ch, void *mm)
{
  poly p, m;
  int np;
  int i;

#ifdef STRACE
  tracePrn ("\n uact_symb (%s)\n", ch);
  writeexpression (vardef, (poly) mm);
  tracePrn (" -> ");
#endif

  m = (poly) mm;
  p = m->next;
  if (strcmp (ch, "G") == 0 || strcmp (ch, "g") == 0)
    {
      while (p != NULL)
	{
	  np = p->tail.tens[0];
	  for (i = 0; i < spinLength; i++)
	    p->tail.power[i] = 0;
	  p->tail.spin.l = 1;
	  p->tail.spin.g5 = 0;
	  p->tail.spin.g[0] = np;
	  p = p->next;
	}
      if (m->next->tail.spin.g[0] == 0)
	m->next->tail.spin.g[0] = 1;
      m->coef.type = spintp;
      if (r_reading2)
	m = (poly) uact_symb ("-", m);
    }
  else if (strcmp (ch, "-") == 0)
    multtensint (&p, -1);

#ifdef STRACE
  writeexpression (vardef, (poly) m);
#endif

  return (void *) m;
}


static void *
bact_symb (char ch, void *mm1, void *mm2)
{
  poly p1, p2, p3;
  int t1, t2, t3;
  int i;
  int d;
  void *mm3;

  if (r_reading2 && (ch == '*'))
    {
      mm3 = mm1;
      mm1 = mm2;
      mm2 = mm3;
    }

#ifdef STRACE
  tracePrn ("\n bact_symb (%c)\n", ch);
  writeexpression (vardef, (poly) mm1);
  tracePrn (" |%c| ", ch);
  writeexpression (vardef, (poly) mm2);
  tracePrn (" -> ");
#endif

  p1 = ((poly) mm1)->next;
  p2 = ((poly) mm2)->next;
  t1 = (int) (((poly) mm1)->coef.type);
  t2 = (int) (((poly) mm2)->coef.type);
  delmonom (((poly *) & mm1));
  delmonom (((poly *) & mm2));

  if ((ch == '+' || ch == '*' || ch == '.') && t1 < t2)
    {
      p3 = p1;
      p1 = p2;
      p2 = p3;
      t3 = t1;
      t1 = t2;
      t2 = t3;
    }
  t3 = -1;			/*  Error  */
  switch (ch)
    {
    case '/':
      printf ("unexpected division\n");
      exit (1);

    case '+':
      switch (t1)
	{
	case numbertp:
	case polytp:
	case tenstp:
	case vectortp:
	  sewtens (&p1, &p2, tensLength);
	  p3 = p1;
	  t3 = tenstp;
/* ----- */ break;

	case spintp:
	  if (t2 != spintp)
	    {
	      poly pp = p2;
	      for (; pp; pp = pp->next)
		for (i = tensLength; i < spinLength; i++)
		  pp->tail.power[i] = 0;
	    }
	  sewtens (&p1, &p2, spinLength);
	  p3 = p1;
	  t3 = spintp;

/* ----- */ break;
	}			/*  Case */
      break;

    case '*':
      switch (t1)
	{
	case numbertp:
	case polytp:
	case tenstp:
	case vectortp:
	case indextp:
	  p3 = multtwotens (p1, p2);
	  t3 = t1;
	  break;

	case spintp:
	  if (t2 != spintp)
	    {
	      poly pp = p2;
	      for (; pp; pp = pp->next)
		for (i = tensLength; i < spinLength; i++)
		  pp->tail.power[i] = 0;
	    }
	  p3 = multtwospin (p1, p2, FALSE);
	  deltensor (&p1);
	  deltensor (&p2);
	  t3 = spintp;
	  break;
	}			/* Case */
      break;

    case '^':
      d = (p2->coef.complex.re->coef.num);
      deltensor (&p2);
      p3 = copytens (p1, tensLength);

      for (i = 1; i < d; i++)
	{
	  p2 = copytens (p1, tensLength);
	  p3 = multtwotens (p3, p2);
	}
      deltensor (&p1);
      t3 = polytp;
      break;

    case '.':
      p3 = multtwotens (p1, p2);
      if (t1 == vectortp && t2 == vectortp)
	t3 = polytp;
      else
	t3 = tenstp;
    }				/* Case */

  newmonom (&p1);
  p1->coef.type = t3;
  p1->next = p3;

#ifdef STRACE
  writeexpression (vardef, (poly) p1);
#endif

  return (void *) p1;
}

void *
act_symb (char *name, int n, void **args)
{
  if (n == 1)
    return uact_symb (name, args[0]);
  if (name[0] == '-')
    {
      uact_symb (name, args[0]);
      strcpy (name, "+");
    }
  if (n == 2)
    return bact_symb (name[0], args[0], args[1]);
  if (n == 4 && !strcmp (name, "eps"))
    {
      poly parg[4], p;
      int i;
      poly ans = NULL;

      for (i = 0; i < 4; i++)
	parg[i] = ((poly) args[i])->next;

      for (; parg[0]; parg[0] = parg[0]->next)
	for (; parg[1]; parg[1] = parg[1]->next)
	  for (; parg[2]; parg[2] = parg[2]->next)
	    for (; parg[3]; parg[3] = parg[3]->next)
	      {
		int l[4], sgn;
		poly eps;
		newmonom (&eps);
		eps->next = NULL;
		eps->coef.complex.re = plusone ();
		eps->coef.complex.im = NULL;
		for (i = 0; i < maxIndex; i++)
		  eps->tail.tens[i] = 0;

		for (i = 0; i < 4; i++)
		  {
		    l[i] = parg[i]->tail.tens[0];
		    if (!l[i])
		      l[i] = 1;

		    multtensComplexpoly (&eps, parg[i]->coef.complex.re,
					 parg[i]->coef.complex.im);
		  }

		for (sgn = 1, i = 0; i < 3 && sgn;)
		  {
		    if (l[i + 1] < l[i])
		      i++;
		    else if (l[i + 1] > l[i])
		      {
			int c = l[i];
			l[i] = l[i + 1];
			l[i + 1] = c;
			sgn *= -1;
			if (i)
			  i--;
			else
			  i++;
		      }
		    else
		      sgn = 0;
		  }
		if (!sgn)
		  deltensor (&eps);
		else
		  {
		    for (i = 0; i < 4; i++)
		      eps->tail.tens[maxIndex + i] = l[i];
		    if (sgn == -1)
		      multtensint (&eps, -1);
		    sewtens (&ans, &eps, tensLength);
		  }
	      }
      for (i = 0; i < 4; i++)
	{
	  deltensor (&(((poly) args[i])->next));
	  delmonom ((poly *) (args + i));
	}

      newmonom (&p);
      p->coef.type = tenstp;
      p->next = ans;

      return p;
    }
  return NULL;
}
