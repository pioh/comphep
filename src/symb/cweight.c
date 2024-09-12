/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Kryukov 
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"

#include "physics.h"
#include "process.h"
#include "colorf.h"
#include "process_core.h"
#include "cweight.h"


static int 
cerror (int n, char *s)
{
  fprintf (stderr, "comphep (error %i): %s\n", n, s);
  exit (99);
  return 0;
}

/* Return color type of vertex - 06/01/90  */
static vtype 
typev (vert0 v, int valence)
{
  int i;
  int ng = 0;
  int nq = 0;

  for (i = 0; i < valence; ++i)
    {
      if (prtclbase[v[i].partcl - 1].cdim != 1)
	{
	  if (prtclbase[v[i].partcl - 1].cdim == 8)
	    ng++;
	  else
	    nq++;
	}
    }

  switch (ng)
    {
    case 0:
      return nq == 2 ? tv : zv;
    case 1:
      return qg;
    case 2:
      return g2;
    case 3:
      return g3;
    default:
      return cerror (252, "invalid color vertex type");
    }
}


/* Transfer Taranov's representation of graph to Kryukov's representation  */
/* - 07/01/90  */
static void 
t2k2 (vcsect * g, int *nv, cvertex * vl)
{
  int i, j;
  int k;
  int en = 0;
  int maxnv = *nv;
  int maptar[2 * maxvert][MAXVALENCE];

  for (i = 0; i < 2 * maxvert; i++)
    for (k = 0; k < MAXVALENCE; k++)
      maptar[i][k] = 0;

  *nv = 0;
  for (i = 0; i < g->sizet; i++)
    if (typev (g->vertlist[i], g->valence[i]) != zv)
      {
	if (*nv >= maxnv)
	  cerror (251, "too many vertices");
	vl[*nv].e[0] = 0;
	vl[*nv].e[1] = 0;
	vl[*nv].e[2] = 0;
	vl[*nv].vt = typev (g->vertlist[i], g->valence[i]);

	for (k = 0, j = 0; j < g->valence[i]; ++j)
	  {
	    int dim = prtclbase[g->vertlist[i][j].partcl - 1].cdim;
	    if (dim != 1)
	      {
		int l = maptar[i][j];
		if (!l)
		  {
		    l = ++en;
		    maptar[i][j] = l;
		    maptar[g->vertlist[i][j].link.vno]
		      [g->vertlist[i][j].link.edno] = l;
		  }

		if (vl[*nv].vt != g3 && vl[*nv].vt != g2) {
		  switch (dim)
		    {
		    case 8:
		      vl[*nv].e[0] = l;
		      break;
		    case -3:
		      vl[*nv].e[1] = l;
		      break;
		    case 3:
		      vl[*nv].e[2] = l;
		      break;
		    default:
		      cerror (253, "invalid particle color index");
		    }
		} else {
		  vl[*nv].e[k] = l;
		  ++k;
	        }
	      }
	  }
	(*nv)++;
      }
}


static int 
a2k (vampl * g, int nc, int *chains, int *nv, cvertex * vl)
{
  int i, k;
  int ne;
  int en = 0;
  int maxnv = *nv;
  int maptar[maxvert][MAXVALENCE];
  int endl[MAXINOUT];
  int endc[MAXINOUT];

  for (i = 0; i < MAXINOUT; i++)
    {
      endl[i] = 0;
      endc[i] = 0;
    }

  for (i = 0; i < maxvert; i++)
    for (k = 0; k < MAXVALENCE; k++)
      maptar[i][k] = 0;
  *nv = 0;
  for (i = 0; i < g->size; i++)
    if (typev (g->vertlist[i], g->valence[i]) != zv)
      {
	if (*nv >= maxnv)
	  return 1;
	vl[*nv].e[0] = 0;
	vl[*nv].e[1] = 0;
	vl[*nv].e[2] = 0;
	vl[*nv].vt = typev (g->vertlist[i], g->valence[i]);

	for (k = 0, ne = 0; ne < g->valence[i]; ne++)
	  {
	    int dim = prtclbase[g->vertlist[i][ne].partcl - 1].cdim;
	    if (dim != 1)
	      {
		int l = maptar[i][ne];
		if (!l)
		  {
		    l = ++en;
		    if (g->vertlist[i][ne].link.vno != nullvert)
		      {
			maptar[i][ne] = l;
			maptar[g->vertlist[i][ne].link.vno]
			  [g->vertlist[i][ne].link.edno] = l;
		      }
		    else
		      {
			int k = g->vertlist[i][ne].link.edno;
			endl[k] = l;
			endc[k] = prtclbase[g->vertlist[i][ne].partcl - 1].cdim;
		      }
		  }

		if (vl[*nv].vt != g3 && vl[*nv].vt != g2)
		  switch (dim)
		    {
		    case 8:
		      vl[*nv].e[0] = l;
		      break;
		    case -3:
		      vl[*nv].e[1] = l;
		      break;
		    case 3:
		      vl[*nv].e[2] = l;
		      break;
		    default:
		      return 2;
		    }
		else
		  vl[*nv].e[k++] = l;
	      }
	  }
	(*nv)++;
      }

  for (i = 0; i < MAXINOUT; i++)
    if (endc[i] == 8)
      {
	if (*nv >= maxnv)
	  return 1;
	vl[*nv].vt = qg;
	vl[*nv].e[0] = endl[i];
	endl[i] = ++en;
	vl[*nv].e[2] = en;
	vl[*nv].e[1] = en + 1;
	en++;
	(*nv)++;
      }

  for (i = 0; i < nc; i++)
    {
      int from = chains[2 * i + 1];
      int to = chains[2 * i];

      if (*nv >= maxnv)
	return 1;

      vl[*nv].vt = tv;
      vl[*nv].e[0] = 0;

      vl[*nv].e[1] = endl[from];
      if (endc[to] == -3)
	vl[*nv].e[2] = endl[to];
      else
	vl[*nv].e[2] = endl[to] + 1;
      (*nv)++;
    }
  return 0;
}


static int 
maxNcPower (vcsect * g)
{
  int i, j;
  int power = 0;

  for (i = 0; i < g->sizel; i++)
    for (j = 0; j < g->valence[i]; j++)
      if (g->vertlist[i][j].link.vno >= g->sizel)
	{
	  int np = g->vertlist[i][j].partcl;
	  switch (prtclbase[np - 1].cdim)
	    {
	    case 3:
	    case -3:
	      power++;
	      break;
	    case 8:
	      power += 2;
	      break;
	    }
	}
  return power / 2;
}

static void 
getLeadingTerm (factor * f, int maxP, long *num, long *den)
{

  if (maxP < f->len - f->dpow - 1)
    {
      cerror (254, "strange behaviour of the program");
    }

  if (maxP > f->len - f->dpow - 1)
    {
      *num = 0;
      *den = 1;
    }
  else
    {
      *num = f->nc[f->len - 1];
      *den = f->dc;
      while (maxP--)
	(*num) *= 3;
    }
}


void 
c_basis_coef (vampl * g, int pow, int nc, int *chains, long *num, long *den)
{
  int i, nv;
  cvertex vl[3 * MAXINOUT];
  factor *f;

  if (!pow)
    return;

  for (i = 0; i < pow; i++)
    {
      nv = 3 * MAXINOUT;
      a2k (g, nc, chains + 2 * nc * i, &nv, vl);

      f = colorFactor (nv, vl);
      getLeadingTerm (f, nc, num + i, den + i);

      free (f->nc);
      free (f);
    }
}


void 
cwtarg (vcsect * g)
{
  factor *f;
  cvertex vl[2 * MAXINOUT];
  int nv = 2 * MAXINOUT;
  t2k2 (g, &nv, vl);
  f = colorFactor (nv, vl);

  if (getNcinflimit ())
    getLeadingTerm (f, maxNcPower (g), &(g->clrnum), &(g->clrdenum));
  else
    fct_num_calc (f, 3, &(g->clrnum), &(g->clrdenum));

  free (f->nc);
  free (f);
}


static void 
lreduce (long *l1, long *l2)
{
  long c, i1, i2;

  i1 = *l1;
  i2 = *l2;
  if (i1 < 0)
    i1 = -i1;
  if (i2 < 0)
    i2 = -i2;

  if (i2 > i1)
    {
      c = i1;
      i1 = i2;
      i2 = c;
    }
  while (i2 != 0)
    {
      c = i2;
      i2 = i1 % i2;
      i1 = c;
    }
  (*l1) /= i1;
  (*l2) /= i1;
}


int 
generateColorWeights (csdiagram * csdiagr, int cBasisPower, int nC, int *cChains,
		      long *cCoefN, long *cCoefD)
{
  vcsect vcs;
  int i;
  int NcInfLimit_tmp = getNcinflimit ();
  int numtot = getntot ();

  transfdiagr (csdiagr, &vcs);
  cwtarg (&vcs);

  setNcinflimit (1);

  if (vcs.clrnum)
    {
      long n = 1;
      long d = 1;
      long *cCoefNr = malloc (cBasisPower * sizeof (long));
      long *cCoefDr = malloc (cBasisPower * sizeof (long));
      vampl left;
      vampl right;

      decompose (vcs, &left, &right);
      c_basis_coef (&left, cBasisPower, nC, cChains, cCoefN, cCoefD);
      c_basis_coef (&right, cBasisPower, nC, cChains, cCoefNr, cCoefDr);

      for (i = 0; i < nC; i++)
	d *= 3;

      for (i = 0; i < numtot; i++)
	{
	  vertlink Q = left.outer[i];
	  if (8 == prtclbase[left.vertlist[Q.vno][Q.edno].partcl - 1].cdim)
	    n *= 2;
	}

      for (i = 0; i < right.size; i++)
	{
	  int n8 = 0;
	  int j;

	  for (j = 0; j < right.valence[i]; j++)
	    if (8 == prtclbase[right.vertlist[i][j].partcl - 1].cdim)
	      n8++;
	  if (n8 == 3)
	    n *= -1;
	}

      for (i = 0; i < cBasisPower; i++)
	{
	  cCoefN[i] *= cCoefNr[i] * n * vcs.clrdenum;
	  cCoefD[i] *= cCoefDr[i] * d * vcs.clrnum;
	  lreduce (cCoefN + i, cCoefD + i);
	}

      for (i = 0, n = 0, d = 1; i < cBasisPower; i++)
	{
	  n = cCoefN[i] * d + n * cCoefD[i];
	  d *= cCoefD[i];
	  lreduce (&n, &d);
	}

      free (cCoefNr);
      free (cCoefDr);
    }
  else
    for (i = 0; i < cBasisPower; i++)
      {
	cCoefN[i] = 0;
	cCoefD[i] = 0;
      }

  setNcinflimit (NcInfLimit_tmp);
  return vcs.clrnum;
}

static int *used, *perm, *wrt, *pos3, *pos_3;
static int nc_;

static void 
recurGen (int k)
{
  int i;

  if (k == nc_)
    for (i = 0; i < nc_; i++)
      {
	*(wrt++) = pos_3[i];
	*(wrt++) = pos3[perm[i]];
      }
  else
    {
      for (i = 0; i < nc_; i++)
	if (!(used[i] || pos_3[i] == pos3[k]))
	  {
	    used[i] = 1;
	    perm[i] = k;
	    recurGen (k + 1);
	    used[i] = 0;
	  }
    }
}


int 
infCbases (int np, int *cweight, int *nc, int *pow, int **chains)
{
  int n3 = 0, n_3 = 0;
  int i;

  pos3 = (int *) malloc (np * sizeof (int));
  pos_3 = (int *) malloc (np * sizeof (int));

  for (i = 0; i < np; i++)
    switch (cweight[i])
      {
      case -3:
	pos_3[n_3++] = i;
	break;
      case 1:
	break;
      case 3:
	pos3[n3++] = i;
	break;
      case 8:
	pos3[n3++] = i;
	pos_3[n_3++] = i;
	break;
      default:
	return 1;
      }
  if (n3 != n_3)
    return 2;

  nc_ = n3;
  if (nc_)
    {
      int pow_ = 1;
      for (i = 2; i <= nc_; ++i)
        pow_ *= i;

      wrt = (int *) malloc (2 * nc_ * pow_ * sizeof (int));
      *chains = wrt;

      used = (int *) malloc (nc_ * sizeof (int));
      for (i = 0; i < nc_; ++i)
        used[i] = 0;
      perm = (int *) malloc (nc_ * sizeof (int));

      recurGen (0);

      free (perm);
      free (used);
      *pow = (wrt - (*chains)) / (2 * nc_);
      *chains = (int *) realloc (*chains, 2 * nc_ * (*pow) * sizeof (int));
      *nc = nc_;
    }
  else
    {
      *nc = 0;
      *pow = 0;
      *chains = NULL;
    }

  free (pos3);
  free (pos_3);
  return 0;
}




  /* *********************** Modules dependence ************************* */
  /* *                  +------------------+     +---------------+      * */
  /* *               +->|1. CWeight        |<--->|0. "CompHep"   |      * */
  /* *               |  +------------------+     +---------------+      * */
  /* *               |        A    |   A                 A              * */
  /* *               |        |    V   +-------+---------+              * */
  /* * +-----------+ |  +-----A------------+   | +---------------+      * */
  /* * |2. Color   |-+->|3.  Tar2Kr        |<--+-|4. Physics     |      * */
  /* * +-----------+    +------------------+     +---------------+      * */
  /* *                                                                  * */
  /* * 0. Module CWeight imported from "CompHEP" Feynman's graph in     * */
  /* *    Taranov's representation (see module Physics) by use          * */
  /* *    function CWTarG from module CWeight.                          * */
  /* * 1. Module CWeight transfer graph to module Tar2Kr for rebuildung * */
  /* *    in Kryukov's representation (see module Color) by use         * */
  /* *    function T2K, calculated color weight and return result       * */
  /* *    (pair of two int - num. and den. - to "CompHEP".          * */
  /* * 2. Module Color exported types, variables, constants and so on   * */
  /* *    for work with color graph.                                    * */
  /* * 3. Module Tar2Kr tranform graph from Taranov's to Kryukov's      * */
  /* *    representation by use T2K procedure.                          * */
  /* * 4. Module Physics exported types, variables, constants and so    * */
  /* *    to "CompHEP" and module Tar2Kr.                               * */
  /* *                              Good luck!                          * */
  /* ******************************************************************** */
