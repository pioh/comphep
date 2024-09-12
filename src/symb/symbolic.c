/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/parser.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "polynom/include/polynom.h"
#include "polynom/include/tensor.h"
#include "polynom/include/spinor.h"
#include "chep_crt/include/chep_crt.h"

/*#define STRACE*/

#include "physics.h"
#include "sos.h"
#include "ghosts.h"
#include "cweight.h"
#include "prepdiag.h"
#include "pvars.h"
#include "chess.h"
#include "saveres.h"
#include "pre_read.h"
#include "reader0.h"
#include "reader_s.h"
#include "rfactor.h"
#include "process.h"
#include "process_core.h"
#include "symbolic.h"

#ifdef STRACE
#include "test_wrt.h"
#endif

static int gamma_map[2 * maxvert];
static int maxmomdeg;
static int px, py;
static poly vert[2 * MAXINOUT];
static poly block[2 * MAXINOUT];
static poly fermres[maxvert];
static int calcdiag_sq;

static poly rnum;
static poly factn;
static poly factd;


static poly
polyfactor (polyvars * v, poly p)
{
  int i;
  poly m;
  poly q;
  NUM_TYPE c, num1, num2, numS;

  newmonom (&m);
  m->next = NULL;
  if (p == NULL)
    {
      m->coef.num = 1;
      for (i = 0; i < monomLength; i++)
	m->tail.power[i] = 0;
    }
  else
    {
      num1 = p->coef.num;
      if (num1 < 0)
	{
	  num1 = -num1;
	  numS = -1;
	}
      else
	numS = 1;

      q = p->next;
      while (q)
	{
	  num2 = q->coef.num > 0 ? q->coef.num : -q->coef.num;
	  if (num2 > num1)
	    {
	      c = num1;
	      num1 = num2;
	      num2 = c;
	    }
	  while (num2 != 0)
	    {
	      c = num2;
	      num2 = REST (num1, num2);
	      num1 = c;
	    }
	  q = q->next;
	}
      for (i = 0; i < monomLength; i++)
	m->tail.power[i] = 0;

      for (i = 0; i < v->nvar; i++)
	if (v->vars[i].num <= nmodelvar)
	  {
	    int wrd = v->vars[i].wordpos - 1;
	    unsigned long z_d = v->vars[i].zerodeg;
	    unsigned long m_d = v->vars[i].maxdeg;
	    int deg = (p->tail.power[wrd] / z_d) % m_d;

	    q = p->next;
	    while (q && deg)
	      {
		int deg2 = (q->tail.power[wrd] / z_d) % m_d;
		if (deg2 < deg)
		  deg = deg2;
		q = q->next;
	      }
	    if (deg)
	      m->tail.power[wrd] += deg * z_d;
	  }

      num1 *= numS;
      while (p)
	{
	  for (i = 0; i < monomLength; i++)
	    p->tail.power[i] -= m->tail.power[i];
	  p->coef.num = DIV (p->coef.num, num1);
	  p = p->next;
	}
      m->coef.num = num1;

    }
  return m;
}


static void
memoryInfo_ (int used)
{
  goto_xy (14, 16);
  print ("%d Kb    ", used >> 10);
  if (escpressed ())
    save_sos (-2);
}


static void
wrtoperat (char *s)
{
  scrcolor (Blue, BGmain);
  goto_xy (14, 17);
  clr_eol ();
  print (s);
}


static void
firstvertexreading (polyvars * v)
{
  int i;
  preres m;
  preres m_;

  clearVars (v);

  for (i = 0; i < vcs.sizet; ++i)
    {
      void *agrs[2];
      vardef = v;
      m_ = (preres) readExpression (vertexes[i].lgrnptr->description,
				    rd_pre, act_pre, NULL);
      v = vardef;
      if (rderrcode)
	save_sos (rderrcode);
      if (m_->g5)
	fermloops[fermmap[i] - 1].g5 = 1;
      gamma_map[i] = m_->maxg;
      m_->g5 = 0;
      m_->maxg = 0;
      m_->indlist = setof (_E);
      m_->tp = polytp;
      if (i)
	{
	  agrs[0] = m;
	  agrs[1] = m_;
	  m = (preres) act_pre ("*", 2, agrs);
	}
      else
	{
	  m = m_;
	}
      if (rderrcode)
	save_sos (rderrcode);
    }

  for (i = 0; i < v->nvar; ++i)
    {
      v->vars[i].maxdeg = m->varsdeg[i] + 1;
    }
  maxmomdeg = m->degp;
  levi = 1;
}


static void
calctenslength (void)
{
  int i, j;

  maxIndex = 0;
  for (i = 0; i < vcs.sizet; ++i)
    {
      for (j = 0; j < vcs.valence[i]; ++j)
	{
	  maxIndex = MAX (maxIndex, vcs.vertlist[i][j].lorentz);
	}
    }
  if (maxIndex == 0)
    maxIndex = 1;
  if (levi)
    tensLength = 1 + (maxIndex + 3) / sizeof (long);
  else
    tensLength = (maxIndex + sizeof (long) - 1) / sizeof (long);
  if (tensLength == 0)
    tensLength = 1;
#ifdef STRACE
  tracePrn ("\n tensLength=%d", tensLength);
#endif
}


static void
calcspinlength (void)
{
  int n, currentlen, nl;

  spinLength = 0;
  for (nl = 0; nl < nloop; nl++)
    {
      fermloopstp *with1 = &fermloops[nl];
      currentlen = 0;
      for (n = 0; n < with1->len; n++)
	{
	  currentlen += gamma_map[with1->vv[n] - 1] + 1;
	  if (spinLength < currentlen)
	    spinLength = currentlen;
	  if (with1->intln[n] != 0 &&
	      !insetb (vcs.vertlist[with1->vv[n] - 1]
		       [with1->intln[n] - 1].lorentz, setmassindex))
	    currentlen -= 2;
	}
    }
  spinLength = 1 + (spinLength + 1) / sizeof (long);
}


static void
propagatorsfirstreading (polyvars * v)
{
  int ninout;
  int i, d, n, l, vln;
  char name[7];

  ninout = getntot ();
  for (i = 1; i <= vcs.sizet; ++i)
    {
      vln = vcs.valence[i - 1];
      for (l = 1; l <= vln; l++)
	{
	  edgeinvert *with1 = &vcs.vertlist[i - 1][l - 1];
	  if (with1->moment > 0)
	    {
	      d = prtclbase[with1->partcl - 1].spin;
	      if (d == 2 && !insetb (with1->lorentz, setmassindex))
		d = 0;
	      if (d == 4 && !insetb (with1->lorentz, setmassindex0))
		d = 2;
	      if (d != 0)
		{
		  maxmomdeg += MIN (d, 2);
		  if (with1->moment > ninout)
		    {
		      strcpy (name, prtclbase[with1->partcl - 1].massidnt);
		      if (strcmp (name, "0") != 0)
			{
			  n = 1;
			  while (n <= ninout
				 && strcmp (name, inoutmasses[n - 1]) != 0)
			    ++n;
			  if (n > ninout)
			    addvar (v, name, d);
			  if (d == 4)
			    addvar (v, name, 2);
			}
		    }
		}
	    }
	}
    }
}


static void
addinoutmasses (polyvars * v)
{
  int n, m;
  char massname[7];
  int numtot = getntot ();

  for (n = 1; n <= numtot; n++)
    {
      strcpy (massname, inoutmasses[n - 1]);
      if (strcmp (massname, "0") != 0)
	{
	  m = 1;
	  while (strcmp (massname, inoutmasses[m - 1]) != 0)
	    ++(m);
	  if (n == m)
	    addvar (v, inoutmasses[n - 1], maxmomdeg);
	}
    }
}


static void
addscmult (polyvars * v)
{
  int i, j;
  char pname[MAXINOUT][4];
  int numtot = getntot ();

  px = 0;
  py = 0;
  for (i = 0; i < 3 * maxvert; ++i)
    {
      if (momdep[i][0] == 1)
	++px;			/* For momdep this is a length. V.E. */
      if (momdep[i][0] > 0)
	++py;
    }
  for (i = 0; i < px; ++i)
    sprintf (pname[i], "p%d", i + 1);
  maxmomdeg /= 2;
  for (i = 1; i < px; ++i)
    {
      for (j = 0; j < i; ++j)
	{
	  addvar (v, scat ("%s.%s", pname[j], pname[i]), maxmomdeg);
	}
    }

#ifdef STRACE
  tracePrn ("\n variable before closing\n");
  writevars (vardef);
#endif

  closevars (v);
#ifdef STRACE
  tracePrn ("\n variable after closing\n");
  writevars (vardef);
#endif
  maxLength = MAX (MAX (tensLength, spinLength), monomLength);

#ifdef STRACE
  tracePrn ("\n maxLength= %d", maxLength);
#endif

  dellevi = 0;
  contracts = (poly *) getmem_ ((sizeof (pointer) * py * (py + 1)) / 2);

  vardef = v;
  for (i = 1; i <= numtot; ++i)
    {
      if (strcmp (inoutmasses[i - 1], "0") == 0)
	assignsclmult (-i, -i, NULL);
      else
	{
	  poly ee = (poly) rd_symb (inoutmasses[i - 1]);
	  assignsclmult (-i, -i, multtwopoly (ee->next->coef.complex.re,
					      ee->next->coef.complex.re));
	  delmonom (&(ee->next));
	  delmonom (&ee);
	}
    }

  for (i = 2; i <= px; ++i)
    for (j = 1; j < i; ++j)
      {
	poly ee = (poly) rd_symb (scat ("%s.%s", pname[j - 1], pname[i - 1]));
	assignsclmult (-j, -i, ee->next->coef.complex.re);
	delmonom (&(ee->next));
	delmonom (&ee);
      }
  v = vardef;

  for (i = px + 1; i <= py; ++i)
    for (j = 1; j <= i; ++j)
      if (i != j || i > numtot)
	{
	  poly qq = NULL;
	  int ii, jj;
	  for (ii = 1; ii <= momdep[i - 1][0]; ++ii)
	    {
	      char p1 = (char) momdep[i - 1][ii];
	      for (jj = 1; jj <= momdep[j - 1][0]; ++jj)
		{
		  char p2 = (char) momdep[j - 1][jj];
		  poly pp = copypoly (scalarmult (-abs (p1), -abs (p2)));
		  if (p1 * p2 < 0)
		    multpolyint (&pp, -1);
		  sewpoly (&qq, &pp);
		}
	    }
	  assignsclmult (-j, -i, qq);
	}
}


static void
secondvertexreading (polyvars * vv)
{
  poly m;
  int v, l, ll;
  int sgn;

#ifdef STRACE
  tracePrn ("\n vertex reading \n");
  tracePrn ("tensLength=%d\n", tensLength);
  tracePrn ("spinLength=%d\n", spinLength);
#endif

  for (v = 1; v <= vcs.sizet; v++)
    {
      for (l = 1; l <= MAX (1, vcs.valence[v - 1]); l++)
	{
	  ll = vertexes[v - 1].subst[l - 1];
	  momsubst[l - 1] = vcs.vertlist[v - 1][ll - 1].moment;
	  indsubst[l - 1] = vcs.vertlist[v - 1][ll - 1].lorentz;
	}
      r_reading2 = vertexes[v - 1].r_vert;

      vardef = vv;
      m = (poly) readExpression (vertexes[v - 1].lgrnptr->description,
				 rd_symb, act_symb, NULL);
      vv = vardef;

#ifdef STRACE
      tracePrn ("\n VERTEX %d", v);
      writeexpression (vv, m);
      tracePrn ("\n");
#endif

      vert[v - 1] = m->next;
      delmonom (&m);
      if (vcs.valence[v - 1] == 3 &&
	  prtclbase[vcs.vertlist[v - 1][0].partcl - 1].cdim == 8 &&
	  prtclbase[vcs.vertlist[v - 1][1].partcl - 1].cdim == 8 &&
	  prtclbase[vcs.vertlist[v - 1][2].partcl - 1].cdim == 8)
	{
	  sgn = 1;
	  if (vertexes[v - 1].subst[0] > vertexes[v - 1].subst[1])
	    sgn = -sgn;
	  if (vertexes[v - 1].subst[1] > vertexes[v - 1].subst[2])
	    sgn = -sgn;
	  if (vertexes[v - 1].subst[0] > vertexes[v - 1].subst[2])
	    sgn = -sgn;
	  if (sgn == -1)
	    multtensint (&vert[v - 1], -1);
	}
    }
}


static void
masscalc (int v, int l, int deg, poly * p)
{
  char masstxt[7];
  poly m;

  strcpy (masstxt, prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].massidnt);
  if (strcmp (masstxt, "0") == 0)
    *p = NULL;
  else
    {
      m = (poly) rd_symb (masstxt);
      if (deg == 1)
	*p = m->next->coef.complex.re;
      else
	{
	  *p =
	    multtwopoly (m->next->coef.complex.re, m->next->coef.complex.re);
	  delpoly (&(m->next->coef.complex.re));
	}
      delmonom (&(m->next));
      delmonom (&m);
    }
}


static void
kmorgcalc (int v, int l, poly * p)
{
  int k;
  int i, m;

  newmonom (p);
  (*p)->next = NULL;
  (*p)->coef.complex.re = plusone ();
  (*p)->coef.complex.im = NULL;
  for (i = 0; i < tensLength; i++)
    (*p)->tail.power[i] = 0;
  k = vcs.vertlist[v - 1][l - 1].moment;
  m = vcs.vertlist[v - 1][l - 1].lorentz;
  (*p)->tail.tens[m - 1] = -abs (k);
  if (k < 0)
    multtensint (p, -1);
}


static void
kmsubcalc (int v, int l, poly * p)
{
  int ll;
  int k;
  int i, m;
  poly q;

  m = vcs.vertlist[v - 1][l - 1].lorentz;
  *p = NULL;
  for (ll = 1; ll <= vcs.valence[v - 1]; ll++)
    if (l != ll)
      {
	newmonom (&q);
	q->next = NULL;
	q->coef.complex.re = plusone ();
	q->coef.complex.im = NULL;
	for (i = 0; i < tensLength; i++)
	  q->tail.power[i] = 0;
	k = vcs.vertlist[v - 1][ll - 1].moment;
	q->tail.tens[m - 1] = -abs (k);
	if (k > 0)
	  multtensint (&q, -1);
	sewtens (p, &q, tensLength);
      }
}


static poly
fermpropag (int v, int l)
{
  char proptxt[128];
  char mass[7];
  poly m, p;

  strcpy (mass, prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].massidnt);
  if (1)
    {
      if (prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].hlp == '*')
	strcpy (proptxt, mass);
      else
	{
	  sprintf (proptxt, "G(p%c)", l + '0');
	  if (strcmp (mass, "0"))
	    {
	      strcat (proptxt, "+");
	      strcat (proptxt, mass);
	    }
	  else if (fermionp (vcs.vertlist[v - 1][l - 1].partcl))
	    {
	      if (prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].hlp == 'L')
		strcat (proptxt, "*(1-G5)");
	      else if (prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].hlp ==
		       'R')
		strcat (proptxt, "*(1+G5)");
	    }
	  else
	    {
	      if (prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].hlp == 'R')
		strcat (proptxt, "*(1-G5)");
	      else if (prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].hlp ==
		       'L')
		strcat (proptxt, "*(1+G5)");
	    }
	}
    }
  momsubst[l - 1] = vcs.vertlist[v - 1][l - 1].moment;
  r_reading2 = 0;
  m = (poly) readExpression (proptxt, rd_symb, act_symb, NULL);
  p = (m->next);
  delmonom (&m);
  return p;
}


static int
g5_test (poly sp)
{
  for (; sp; sp = sp->next)
    {
      if (sp->tail.spin.g5)
	return 1;
    }
  return 0;
}


static void
calcfermloops (void)
{
  int v, l, l1, l2, n, lpcount, m, i, nmassind, nv;
  poly mult_1, mult_2, mult_1_, mult_2_, frmprpg, fctmp, km1, km2, m_2,
    m_2_, q, sum;
  poly subv[4], subs[4], mass2[4];
  poly pmem;

/* Main Procedure -- FermPrgEmit */
  for (lpcount = 1; lpcount <= nloop; lpcount++)
    {
      fermloopstp *with1 = &fermloops[lpcount - 1];

      frmprpg = fermpropag (with1->vv[0], with1->ll[0]);

#ifdef STRACE
      tracePrn ("\nfermion loop %d  calculation\n", lpcount);
      tracePrn ("\nvertex=");
      writespinor (vardef, vert[with1->vv[0] - 1]);
      tracePrn ("\n propagator=");
      writespinor (vardef, frmprpg);
#endif
      fctmp = multtwospin (vert[with1->vv[0] - 1], frmprpg, 0);

#ifdef STRACE
      tracePrn ("\n result=");
      writespinor (vardef, fctmp);
#endif
      deltensor (&frmprpg);
      for (n = 2; n <= with1->len; n++)
	{
	  mult_1 = fctmp;
	  frmprpg = fermpropag (with1->vv[n - 1], with1->ll[n - 1]);
#ifdef STRACE
	  tracePrn ("\nvertex=");
	  writespinor (vardef, vert[with1->vv[n - 1] - 1]);
	  tracePrn ("\n propagator=");
	  writespinor (vardef, frmprpg);
#endif

	  mult_2 = multtwospin (vert[with1->vv[n - 1] - 1], frmprpg, 0);
	  deltensor (&frmprpg);
	  fctmp = multtwospin (mult_1, mult_2, n == with1->len);

	  l1 = with1->intln[n - 1];
	  if (l1 &&
	      !insetb (vcs.vertlist[with1->vv[n - 1] - 1][l1 - 1].lorentz,
		       setmassindex))
	    l1 = 0;

	  l2 = with1->intln2[n - 1];
	  if (l2 &&
	      !insetb (vcs.vertlist[with1->vv[n - 1] - 1][l2 - 1].lorentz,
		       setmassindex))
	    l2 = 0;

	  if (!l1 && l2)
	    {
	      l = l1;
	      l1 = l2;
	      l2 = l;
	    }

	  if (l2)
	    {
	      mult_1_ = copytens (mult_1, spinLength);
	      mult_2_ = copytens (mult_2, spinLength);
	    }

	  if (l1)
	    {
	      v = vcs.vertlist[with1->vv[n - 1] - 1][l1 - 1].nextvert.vno;
	      l = vcs.vertlist[with1->vv[n - 1] - 1][l1 - 1].nextvert.edno;
	      masscalc (v, l, 2, &m_2);
	      multtenspoly (&fctmp, m_2);

	      kmsubcalc (v, l, &km1);
	      multspintens (&mult_1, &km1);
	      deltensor (&km1);

	      kmsubcalc (with1->vv[n - 1], l1, &km2);
	      multspintens (&mult_2, &km2);
	      deltensor (&km2);

	      q = multtwospin (mult_1, mult_2, (int) (n == with1->len));
	      sewtens (&fctmp, &q, spinLength);
	      if (l2)
		{
		  v = vcs.vertlist[with1->vv[n - 1] - 1][l2 - 1].nextvert.vno;
		  l =
		    vcs.vertlist[with1->vv[n - 1] - 1][l2 - 1].nextvert.edno;

		  masscalc (v, l, 2, &m_2_);
		  multtenspoly (&fctmp, m_2_);
		  delpoly (&m_2_);

		  kmsubcalc (v, l, &km1);
		  multspintens (&mult_1_, &km1);

		  kmsubcalc (with1->vv[n - 1], l2, &km2);
		  multspintens (&mult_2_, &km2);

		  q = multtwospin (mult_1_, mult_2_, (int) (n == with1->len));
		  multtenspoly (&q, m_2);
		  sewtens (&fctmp, &q, spinLength);

		  multspintens (&mult_1, &km1);
		  multspintens (&mult_2, &km2);
		  q = multtwospin (mult_1, mult_2, (int) (n == with1->len));
		  sewtens (&fctmp, &q, spinLength);

		  deltensor (&km1);
		  deltensor (&km2);
		  deltensor (&mult_1_);
		  deltensor (&mult_2_);
		}
	      delpoly (&m_2);
	    }
	  deltensor (&mult_1);
	  deltensor (&mult_2);
#ifdef STRACE
	  tracePrn ("\n result=");
	  writespinor (vardef, fctmp);
#endif
	}

      for (nv = 1; nv <= strlen (with1->invrt); nv++)
	{
	  v = (with1->invrt[nv - 1]);
	  nmassind = 0;
	  for (l = 1; l <= vcs.valence[v - 1]; l++)
	    {
	      m = vcs.vertlist[v - 1][l - 1].lorentz;
	      if (m && insetb (m, setmassindex))
		{
		  nmassind++;
		  kmorgcalc (v, l, &subv[nmassind - 1]);
		  kmsubcalc (vcs.vertlist[v - 1][l - 1].nextvert.vno,
			     vcs.vertlist[v - 1][l - 1].nextvert.edno,
			     &subs[nmassind - 1]);
		  masscalc (v, l, 2, &mass2[nmassind - 1]);
		}
	    }
	  if (nmassind == 0)
	    multspintens (&fctmp, &vert[v - 1]);
	  else
	    {
	      sum = NULL;
	      for (n = 0; n <= (1 << nmassind) - 1; n++)
		{
		  mult_1 = copytens (vert[v - 1], tensLength);
		  mult_2 = copytens (fctmp, spinLength);
		  for (i = 1; i <= nmassind; i++)
		    if (((1 << (i - 1)) & n) != 0)
		      {
			multspintens (&mult_2, &subs[i - 1]);
			q = mult_1;
			pmem = copytens (subv[i - 1], tensLength);
			mult_1 = multtwotens (q, pmem);
			/*   DelTensor(Q); */
		      }
		    else
		      multtenspoly (&mult_1, mass2[i - 1]);
		  multspintens (&mult_2, &mult_1);
		  deltensor (&mult_1);
		  sewtens (&sum, &mult_2, spinLength);
		}
	      deltensor (&fctmp);
	      fctmp = sum;
	      for (n = 0; n < nmassind; n++)
		{
		  deltensor (&subv[n]);
		  deltensor (&subs[n]);
		  delpoly (&mass2[n]);
		}
	    }
	  fermmap[v - 1] = lpcount;
	}
      with1->g5 = g5_test (fctmp);
      fermres[lpcount - 1] = calcspur (fctmp);
#ifdef STRACE
      tracePrn ("\n Trace calculation = ");
      writetens (vardef, fermres[lpcount - 1]);
#endif
      for (n = 0; n < with1->len; n++)
	deltensor (&vert[with1->vv[n] - 1]);
      for (nv = 1; nv <= strlen (with1->invrt); nv++)
	deltensor (&vert[(with1->invrt[nv - 1]) - 1]);
    }
}


static void
multblocks (int v1, int v2, int saveEps)
{
  poly sub[15], mass2[15];
  poly t1, t2, sum, q, mult_1, mult_2;
  int i, v, l, vmin, nmassind, n;
  indvertset msind, ind1, ind2;
  int lastn;
  poly pmem;
  char messtxt[STRSIZ];

  t1 = block[v1 - 1];
  t2 = block[v2 - 1];


  setofb_cpy (ind1, vertinfo[v1 - 1].ind);
  setofb_cpy (ind2, vertinfo[v2 - 1].ind);
  strcpy (messtxt, "Indices contraction   ");
  for (i = 1; i <= 15; i++)
    if (insetb (i, ind1) && insetb (i, ind2))
      sprintf (messtxt, "%s L%d", messtxt, i);
  strcat (messtxt, "           ");
  wrtoperat (messtxt);

  setofb_cpy (msind, setofb_its (ind1, ind2));
  setofb_cpy (msind, setofb_its (msind, setmassindex));
  if (setofb_eq0 (msind))
    {
      dellevi = !saveEps;
      sum = multtwotens (t1, t2);
      dellevi = 0;
    }
  else
    {
      nmassind = 0;
      for (i = 1; i <= 15; i++)
	if (insetb (i, msind))
	  {
	    ++(nmassind);
	    v = massindpos[i - 1].vrt1;
	    l = massindpos[i - 1].ln1;
	    kmorgcalc (v, l, &sub[nmassind - 1]);
	    masscalc (v, l, 2, &mass2[nmassind - 1]);
	  }
      sum = NULL;
      lastn = (1 << nmassind) - 1;
      for (n = 0; n <= lastn; n++)
	{
	  if (n == lastn)
	    {
	      mult_1 = t1;
	      mult_2 = t2;
	    }
	  else
	    {
	      mult_1 = copytens (t1, tensLength);
	      mult_2 = copytens (t2, tensLength);
	    }
	  for (i = 1; i <= nmassind; i++)
	    if (((1 << (i - 1)) & n) != 0)
	      {
		q = mult_1;
		pmem = copytens (sub[i - 1], tensLength);
		mult_1 = multtwotens (q, pmem);
		q = mult_2;
		pmem = copytens (sub[i - 1], tensLength);
		mult_2 = multtwotens (q, pmem);
		multtensint (&mult_2, -1);
	      }
	    else
	      multtenspoly (&mult_1, mass2[i - 1]);
	  dellevi = !saveEps;
	  q = multtwotens (mult_2, mult_1);
	  dellevi = 0;
	  sewtens (&sum, &q, tensLength);
	}
      for (i = 1; i <= nmassind; i++)
	{
	  deltensor (&sub[i - 1]);
	  delpoly (&mass2[i - 1]);
	}
    }
  vmin = MIN (v1, v2);
  block[vmin - 1] = sum;
  vertinfo[vmin - 1].g5 = vertinfo[v1 - 1].g5 + vertinfo[v2 - 1].g5;
  setofb_cpy (vertinfo[vmin - 1].ind,
	      setofb_aun (setofb_uni (ind1, ind2), setofb_its (ind1, ind2)));
}


static void
del_pp (polyvars * vv, poly * p, poly * fact, long *del)
{
  int d1, n, i, j;

  unsigned long z_d, m_d;

  poly psub;

  poly *pp_powers;

  poly fact1;

  long dmax;

  poly ans;
  poly q;
  int numin = getnin ();

  if (*p == NULL)
    {
      *fact = plusone ();
      *del = 1;
      return;
    }
  n = getntot ();
  psub = scalarmult (-n, -n);
  multpolyint (&psub, -1);
  for (i = 1; i <= n - 1; i++)
    {
      q = scalarmult (-i, -i);
      sewpoly (&psub, &q);
    }

  for (i = 2; i <= n - 1; i++)
    for (j = 1; j <= i - 1; j++)
      if (j != n - 2)
	{
	  q = scalarmult (-i, -j);
	  if (i > numin && j <= numin)
	    multpolyint (&q, -2);
	  else
	    multpolyint (&q, 2);
	  sewpoly (&psub, &q);
	}

  if (getnout () > 2)
    multpolyint (&psub, -1);

  z_d = vv->vars[0].zerodeg;
  m_d = vv->vars[0].maxdeg;

  fact1 = polyfactor (vv, *p);
  dmax = ((*p)->tail.power[0] / z_d) % m_d;

  pp_powers = m_alloc ((dmax + 1) * sizeof (poly));
  pp_powers[0] = plusone ();
  for (i = 1; i <= dmax; i++)
    pp_powers[i] = multtwopoly (pp_powers[i - 1], psub);


  ans = NULL;
  q = *p;
  while (q)
    {
      poly mon = q;
      poly tmp;
      q = q->next;

      mon->next = NULL;
      d1 = (mon->tail.power[0] / z_d) % m_d;
      mon->tail.power[0] -= d1 * z_d;
      tmp = multtwopoly (pp_powers[d1], mon);
      multpolyint (&tmp, 1 << (dmax - d1));
      sewpoly (&ans, &tmp);
      delpoly (&mon);
    }
  q = polyfactor (vv, ans);
  (*fact) = multtwopoly (fact1, q);
  delpoly (&fact1);
  delpoly (&q);
  *del = 1 << dmax;
  (*p) = ans;
  for (i = 0; i <= dmax; i++)
    delpoly (pp_powers + i);
  free (pp_powers);
}


static void
transformfactor (polyvars * var1, polyvars * var2, rmptr * t_fact, poly mon, long del)
{
  int i, j;
  int c;
  long factnum;
  long factdenum;
  char factortxt[STRSIZ];
  rmptr vard;
  rmptr tf_add;
  preres m;
  poly mm;

  sprintf (factortxt, "%" NUM_STR, mon->coef.num);

  for (i = 0; i < var1->nvar; ++i)
    {
      int deg = (mon->tail.power[var1->vars[i].wordpos - 1] /
		 var1->vars[i].zerodeg) % var1->vars[i].maxdeg;
      if (deg != 0)
	{
	  sprintf (factortxt + strlen (factortxt), "*%s", var1->vars[i].name);
	  if (deg > 1)
	    sprintf (factortxt + strlen (factortxt), "^%d", deg);
	}
    }
  vard = (rmptr) read_rmonom (factortxt);
  mult_rptr (t_fact, &vard);

/* ------- Symmetry and Color Factors ------- */
  factnum = vcs.symnum * vcs.clrnum;
  factdenum = vcs.symdenum * vcs.clrdenum * del;

/* ----- average factor --------- */
  c = 1;
  for (i = 0; i < vcs.sizel; ++i)
    {
      for (j = 0; j < vcs.valence[i]; ++j)
	{
	  if (IN_PRTCL & vcs.vertlist[i][j].prop)
	    {
	      int pnum = vcs.vertlist[i][j].partcl - 1;
	      switch (prtclbase[pnum].spin)
		{
		case 1:
		  if (!(prtclbase[pnum].hlp == 'L' ||
			prtclbase[pnum].hlp == 'R'))
		    {
		      c *= 2;
		    }
		  break;
		case 2:
		  if (zeromass (pnum + 1))
		    c *= 2;
		  else
		    c *= 3;
		}
	      c *= abs (prtclbase[pnum].cdim);
	    }
	}
    }
  factdenum *= c;

/* ----- Fermion factor --------- */
  c = 0;
  for (i = 0; i < vcs.sizel; ++i)
    {
      for (j = 0; j < vcs.valence[i]; ++j)
	{
	  int pnum = vcs.vertlist[i][j].partcl - 1;
	  if (1 == prtclbase[pnum].spin &&
	      (IN_PRTCL & vcs.vertlist[i][j].prop))
	    {
	      ++c;
	    }
	}
    }
  c = (c & 1) == 1 ? -1 : 1;
  for (i = 0; i < nloop; ++i)
    c *= -4;
  factnum *= c;

/* ----- Vector/left spinor factor --------- */
  for (i = 0; i < vcs.sizet; ++i)
    {
      for (j = 0; j < vcs.valence[i]; ++j)
	{
	  if (vcs.vertlist[i][j].moment > 0)
	    {
	      int pnum = vcs.vertlist[i][j].partcl - 1;
	      if (prtclbase[pnum].hlp == 'L' || prtclbase[pnum].hlp == 'R')
		factdenum *= 2;
	    }
	}
    }

/* ----------- end of numeric factors --------------- */
  sprintf (factortxt, "%d", (int) factnum);
  tf_add = (rmptr) read_rmonom (factortxt);
  mult_rptr (t_fact, &tf_add);
  sprintf (factortxt, "%d", (int) factdenum);

  for (i = 0; i < vcs.sizet; ++i)
    {
      for (j = 0; j < vcs.valence[i]; ++j)
	{
	  edgeinvert *with1 = &vcs.vertlist[i][j];
	  if (with1->moment > 0 &&
	      (pseudop (with1->partcl) ||
	       insetb (with1->lorentz, setmassindex0)))
	    sprintf (factortxt + strlen (factortxt), "*%s^2",
		     prtclbase[with1->partcl - 1].massidnt);
	}
    }
  sbld (factortxt, "1/(%s)", factortxt);
  tf_add = (rmptr) read_rmonom (factortxt);
  mult_rptr (t_fact, &tf_add);

  clearVars (var2);
  vardef = var2;
  m = (preres) readExpression (rmonomtxt (**t_fact), rd_pre, act_pre, NULL);
  var2 = vardef;
  if (rderrcode != 0)
    save_sos (rderrcode);
  for (i = 0; i < var2->nvar; ++i)
    {
      var2->vars[i].maxdeg = m->varsdeg[i] + 1;
    }
  closevars (var2);
  tensLength = 0;
  maxIndex = 0;
  spinLength = 0;
  levi = 0;
  maxLength = MAX (MAX (tensLength, spinLength), monomLength);
  vardef = var2;
  mm = (poly) readExpression (smonomtxt ((**t_fact).n), rd_symb, act_symb, NULL);
  factn = mm->next->coef.complex.re;
  delmonom (&mm);
  mm = (poly) readExpression (smonomtxt ((**t_fact).d), rd_symb, act_symb, NULL);
  factd = mm->next->coef.complex.re;
  delmonom (&mm);
  var2 = vardef;

  clrvm ((*t_fact)->n.v);
  clrvm ((*t_fact)->d.v);
  free (*t_fact);
}


static void
formblocks (void)
{
  int i, j;
  int count;
  int vertmap[2 * maxvert];

  for (i = 0; i < 2 * maxvert; ++i)
    {
      vertmap[i] = fermmap[i];
    }

  for (i = 0; i < nloop; ++i)
    {
      block[i] = fermres[i];
    }

  count = nloop;
  for (i = 0; i < vcs.sizet; ++i)
    {
      if (!vertmap[i])
	{
	  block[count] = vert[i];
	  vertmap[i] = ++count;
	}
    }

  /*  Begin of Blk filling  */
  n_vrt = count;
  for (i = 0; i < count; ++i)
    {
      vertinfo[i].vlnc = 0;
      vertinfo[i].weit = 1;
      setofb_zero (vertinfo[i].ind);
      vertinfo[i].g5 = 0;
    }

  for (i = 0; i < vcs.sizet; ++i)
    {
      int v1 = vertmap[i];
      if (v1 <= nloop)
	{
	  vertinfo[v1 - 1].weit += 2;
	  if (fermloops[v1 - 1].g5)
	    ++(vertinfo[v1 - 1].weit);
	}

      for (j = 0; j < vcs.valence[i]; ++j)
	{
	  int v2 = vertmap[vcs.vertlist[i][j].nextvert.vno - 1];
	  int lori = vcs.vertlist[i][j].lorentz;
	  if (lori != 0)
	    {
	      if (v1 != v2)
		{
		  int np;
		  ++(vertinfo[v1 - 1].vlnc);
		  vertinfo[v1 - 1].link[vertinfo[v1 - 1].vlnc - 1] = v2;
		  setofb_cpy (vertinfo[v1 - 1].ind,
			      setofb_uni (vertinfo[v1 - 1].ind,
					  setofb (lori, _E)));
		  np = vcs.vertlist[i][j].partcl;
		  if (prtclbase[np - 1].hlp == 't')
		    {
		      ++(vertinfo[v1 - 1].vlnc);
		      vertinfo[v1 - 1].link[vertinfo[v1 - 1].vlnc - 1] = v2;
		      setofb_cpy (vertinfo[v1 - 1].ind,
				  setofb_uni (vertinfo[v1 - 1].ind,
					      setofb (lori - 1, _E)));
		    }
		}
	      vertinfo[v1 - 1].weit += 2;
	      if (insetb (lori, setmassindex))
		++(vertinfo[v1 - 1].weit);
	    }
	}
    }
  for (i = 0; i < nloop; ++i)
    {
      if (fermloops[i].g5)
	vertinfo[i].g5 = 1;
    }
  if (MEMORY_OPTIM)
    {
      for (i = 0; i < n_vrt; ++i)
	vertinfo[i].weit = 0;
      for (i = 0; i < vcs.sizet; ++i)
	vertinfo[vertmap[i] - 1].weit++;
    }
}


static int
delImageryOne (rmptr t_factor)
{
  vmrec rec;
  vmptr m, m1;

  rec.next = (t_factor->n).v;
  m = &rec;
  m1 = m->next;
  while (m1 != NULL && (strcmp (m1->name, "i") != 0)) 
    {
      m = m1;
      m1 = m1->next;
    }

  if (m1)
    {
      m->next = m1->next;
      free (m1);
      m1 = m->next;
      (t_factor->n).v = rec.next;
      return 1;
    }
  else
    {
      return 0;
    }
}


/* returns calculation status: -2/-1/0/1/2 -> outOfMemory/deleted/Rest/calculated/Zero */
static int
symbcalc (polyvars * var1, polyvars * var2, hlpcsptr ghst)
{
  int i, j;
  int sgn;
  int spinl_s;
  int maxmom_s;
  int first;
  long del;
  polyvars tmpvars = { 0, NULL };
  vcsect vcs_copy;
  hlpcsptr gstcopy;
  s_listptr d_facts;
  s_listptr df;
  rmptr t_fact;
  vmptr coefvar;
  poly mfact;
  poly mon;
  poly pp;
  poly mm;

  goto_xy (14, 15);
  clr_eol ();
  wrtoperat ("Factors normalization");
  diagramsrfactors (ghst, &d_facts, &t_fact);
  wrtoperat ("Preparing for calculation");
  vcs_copy = vcs;
  first = 1;

  gstcopy = ghst;
  df = d_facts;

  coloringvcs (ghst);
  attachvertexes ();
  firstvertexreading (var1); /* clear var1 and fill it out */
  propagatorsfirstreading (var1);
  calcspinlength ();
  coefvar = df->monom.v;
  while (coefvar != NULL)
    {
      addvar (var1, coefvar->name, coefvar->deg);
      coefvar = coefvar->next;
    }
  vcs = vcs_copy;
  df = df->next;
  ghst = ghst->next;

  spinl_s = spinLength;
  maxmom_s = maxmomdeg;

  while (ghst != NULL)
    {
      coloringvcs (ghst);
      attachvertexes ();
      firstvertexreading (&tmpvars); /* clear tmpvars and fill it out */
      propagatorsfirstreading (&tmpvars);
      calcspinlength ();
      coefvar = df->monom.v;
      while (coefvar != NULL)
	{
	  addvar (&tmpvars, coefvar->name, coefvar->deg);
	  coefvar = coefvar->next;
	}
      unite_vars (var1, &tmpvars);
      spinl_s = MAX (spinLength, spinl_s);
      maxmom_s = MAX (maxmomdeg, maxmom_s);
      vcs = vcs_copy;
      df = df->next;
      ghst = ghst->next;
    }

  spinLength = spinl_s;
  maxmomdeg = maxmom_s;

  levi = 1;
#ifdef STRACE
  tracePrn ("\n spinLength= %d", spinLength);
#endif

  calctenslength ();
  addinoutmasses (var1);
  addscmult (var1);
  clearpregarbage ();
  ghst = gstcopy;
  df = d_facts;

  rnum = NULL;

  vardef = var1;
  do
    {
      goto_xy (14, 15);
      print ("%d(of %d)  ", ghst->num, ghst->maxnum);
      coloringvcs (ghst);
      attachvertexes ();
      findReversVert ();
      secondvertexreading (var1);
      wrtoperat ("Fermion loops calculation  ");
      calcfermloops ();
      formblocks ();
      makeprgcode ();

#ifdef STRACE
      tracePrn ("n_vrt= %d", n_vrt);
      for (i = 0; i <= n_vrt - 1; i++)
	{
	  tracePrn ("\n vrt=%d:", i);
	  writetens (vardef, block[i]);
	}
#endif

      for (i = n_vrt - 2; i >= 0; i--)
	{
#ifdef STRACE
	  tracePrn ("\n multiplication level %d ", i);
	  tracePrn ("\n T1= ");
	  writetens (vardef, block[prgcode[i][0] - 1]);
	  tracePrn ("\n T2= ");
	  writetens (vardef, block[prgcode[i][1] - 1]);
#endif
	  multblocks (prgcode[i][0], prgcode[i][1], i);
#ifdef STRACE
	  tracePrn ("\n    Final result=");
	  writetens (vardef, block[MIN (prgcode[i][0], prgcode[i][1]) - 1]);
#endif
	}


      for (i = 0; i < vcs.sizet; i++)
	{
	  for (j = 0; j < vcs.valence[i]; j++)
	    {
	      int mom = vcs.vertlist[i][j].moment;
	      int np = vcs.vertlist[i][j].partcl;
	      if ((mom > 0) && (prtclbase[np - 1].hlp == 't'))
		{
		  masscalc (i + 1, j + 1, 2, &mm);
		  if (prtclbase[ghostmother (np) - 1].hlp == '*')
		    {
		      multtenspoly (&block[0], mm);
		      delpoly (&mm);
		    }
		  else
		    {
		      pp = copypoly (scalarmult (-mom, -mom));
		      multpolyint (&pp, -1);
		      sewpoly (&pp, &mm);
		      multtenspoly (&block[0], pp);
		      delpoly (&pp);
		    }
		  if (insetb (vcs.vertlist[i][j].lorentz, setmassindex0))
		    {
		      masscalc (i + 1, j + 1, 2, &mm);
		      multtenspoly (&block[0], mm);
		      delpoly (&mm);
		    }
		}
	    }
	}

      sgn = ghst->sgn;
      for (i = 0; i < vcs.sizet; ++i)
	{
	  sgn *= vertexes[i].lgrnptr->factor;
	}
      if (sgn != 1)
	{
	  multtensint (&block[0], sgn);
	}
      if (block[0] != NULL)
	{
	  char ssss[STRSIZ];
	  strcpy (ssss, smonomtxt (df->monom));
	  mfact = (poly) readExpression (ssss, rd_symb, act_symb, NULL);
	  multtensComplexpoly (&block[0], mfact->next->coef.complex.re,
			       mfact->next->coef.complex.im);
	  delmonom (&(mfact->next));
	  delmonom (&mfact);
	  sewtens (&rnum, &block[0], tensLength);
	}

      vcs = vcs_copy;
      df = df->next;
      ghst = ghst->next;
    }
  while (ghst != NULL);
  eraseslist (d_facts);

  if (!rnum)
    return 2;

  if (delImageryOne (t_fact))
    {
      poly p = plusone ();
      multtensComplexpoly (&rnum, NULL, p);
      delmonom (&p);
    }

  tensRealPart (&rnum);
  if (!rnum)
    return 2;

  if (consLow)
    {
      del_pp (var1, &(rnum->coef.complex.re), &mon, &del);
    }
  else
    {
      mon = polyfactor (var1, rnum->coef.complex.re);
      del = 1;
    }

  wrtoperat ("Total factor calculation");
  transformfactor (var1, var2, &t_fact, mon, del);

  wrtoperat ("Denominator calculation");
  if (rnum->coef.complex.re && factn)
    return 1;
  else
    return 2;
}


static void
calcproc (polyvars * var1, polyvars * var2, csdiagram * csdiagr)
{
  hlpcsptr gstlist;

  transfdiagr (csdiagr, &vcs);
  cwtarg (&vcs);
  if (vcs.clrnum == 0)
    {
      csdiagr->status = 2;
    }
  else
    {
      generateghosts (&vcs, &gstlist);
      if (gstlist == NULL)
	{
	  csdiagr->status = 2;
	}
      else
	{
	  preperdiagram ();
	  csdiagr->status = symbcalc (var1, var2, gstlist);
	}
      eraseghosts (gstlist);
    }
}


static void
writestatistic (unsigned noutmemtot, char *txt)
{
  shortstr numdiag;
  sprintf (numdiag, "%u", calcdiag_sq);
  goto_xy (1, 12);
  print ("%s\n", numdiag);
  goto_xy (1, 13);
  print ("%u", noutmemtot);
  sprintf (numdiag, "%4d", ndiagr);
  goto_xy (17, 14);
  print ("%4d", ndiagr);
  goto_xy (42, 14);
  print ("%s      ", txt);
  goto_xy (14, 17);
  clr_eol ();
}


static void
heap_is_empty (void)
{
  save_sos (-1);
}


void
calcallproc (int firstdiag, int lastdiag)
{
  int ndel, ncalc, nrest;
  int ntotdiag;
  int status_key = 0;
  long nrecord;
  csdiagram csd;
  unsigned noutmemtot;
  shortstr txt;
  marktp heap_beg;
  FILE *archiv = fopen (ARCHIV_NAME, "ab");

  polyvars v1 = { 0, NULL };
  polyvars v2 = { 0, NULL };

  memerror = heap_is_empty;
  memoryInfo = memoryInfo_;

  goto_xy (1, 13);
  scrcolor (FGmain, BGmain);
  print ("0     Out of memory\n");
  print ("current diagram          in (Sub)process \n");
  print ("Subdiagram  :\n");
  print ("Used memory :%3d Kb       \n", (int) (usedmemory / 1000));
  print ("Operation   :\n");
  scrcolor (Yellow, Blue);
  print ("\n");
  print (" Press Esc to halt calculations \n");
  scrcolor (Blue, BGmain);

  catalog = fopen (CATALOG_NAME, "ab");

  calcdiag_sq = 0;
  noutmemtot = 0;
  diagrq = fopen (DIAGRQ_NAME, "rb");
  while (FREAD1 (csd, diagrq))
    {
      switch (csd.status)
	{
	case 1:
	  ++calcdiag_sq;
	  break;
	case -2:
	  ++noutmemtot;
	  break;
	case 2:
	  ++calcdiag_sq;
	}
    }
  fclose (diagrq);

  diagrq = fopen (DIAGRQ_NAME, "r+b");
  menuq = fopen (MENUQ_NAME, "r+b");
  ntotdiag = 1;
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      int naux;
      rd_menu (menuq, 2, nsub, txt, &ndel, &ncalc, &nrest, &nrecord);
      naux = ndel + ncalc + nrest;
      for (ndiagr = 1; ndiagr <= naux; ndiagr++)
	{
	  writestatistic (noutmemtot, txt);
	  fseek (diagrq, sizeof (csd) * (nrecord + ndiagr - 1), SEEK_SET);
	  FREAD1 (csd, diagrq);
/* -------------- koeff status -----------------------------*/
	  if ( (csd.status == 0) || (csd.status == 3) ) 
	    { 
              status_key = csd.status;
              csd.status = 0;
/*---------------------------------------------------------*/
	      mark_ (&heap_beg);
	      if (firstdiag <= ntotdiag && ntotdiag < lastdiag)
		calcproc (&v1, &v2, &csd);

	      fseek (diagrq, sizeof (csd) * (nrecord + ndiagr - 1), SEEK_SET);
	      FWRITE1 (csd, diagrq);
	      if (firstdiag <= ntotdiag && ntotdiag < lastdiag)
		{
		  if (csd.status == 1)
		    {
		      wrtoperat ("Writing result             ");
		      saveent (6);
		  /*    save_analitic_results (archiv, rnum->coef.complex.re,
					     factn, factd, &v1, &v2, vcs);
                  */
                        
                        save_analitic_results (archiv, rnum->coef.complex.re,
					     factn, factd, &v1, &v2, vcs, status_key);
		    }
		}
	      else
		{
		  save_empty_analitic_results (archiv, vcs);
		}
	      release_ (&heap_beg);
	      clearVars (&v1);
	      clearVars (&v2);
	      ncalc++;
	      nrest--;
	      calcdiag_sq++;
	      wrt_menu (menuq, 2, nsub, txt, ndel, ncalc, nrest, nrecord);
	      if (escpressed ())
		goto exi;
	    }
	  ++ntotdiag;
	}
    }
exi:

  fclose (diagrq);
  fclose (menuq);

  fclose (catalog);
  fclose (archiv);

  scrcolor (FGmain, BGmain);
  clrbox (1, 14, 70, 20);
  memerror = NULL;
}
