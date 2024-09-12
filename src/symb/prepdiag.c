/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"

#include "physics.h"
#include "ghosts.h"
#include "process.h"
#include "process_core.h"
#include "prepdiag.h"

vertexhlp vertexes[2 * maxvert] = {{0,{0},NULL}};

linkhlp massindpos[5 * maxvert] = {{0,0,0,0}};

fermloopstp fermloops[maxvert] = {{0,0,0,{0},{0},{0},{0},{0}}};

indset setmassindex = {0};
indset setmassindex0 = {0};
int nloop = 0;
int fermmap[2 * maxvert] = {0};
char inoutmasses[MAXINOUT][7] = {{0}};
momsum momdep[3 * maxvert] = {{0}};

int consLow = 0;

void 
coloringvcs (hlpcsptr currentghst)
{
  int i, j;

  indset setDel;
  setofb_zero (setDel);
  for (i = 0; i < vcs.sizet; i++)
    for (j = 0; j < vcs.valence[i]; j++)
      {
	vcs.vertlist[i][j].partcl += currentghst->hlpcs[i][j];
	switch (currentghst->hlpcs[i][j])
	  {
	  case ghostmark:
	  case antighostmark:
	  case sbosonmark:
	    vcs.vertlist[i][j].lorentz = 0;
	    break;
	  case tbosonmark:
	    setofb_cpy (setDel,
		setofb_uni (setDel, setofb (vcs.vertlist[i][j].lorentz, _E))
	      );
	  }
      }
  setofb_cpy (setmassindex, setofb_aun (setmassindex0, setDel));

}

void 
attachvertexes (void)
{
  int v, l;
  arr4byte vert;

  for (v = 0; v < vcs.sizet; v++)
    {
      for (l = 0; l < vcs.valence[v]; l++)
	vert[l] = vcs.vertlist[v][l].partcl;
      for (l = vcs.valence[v]; l < MAXVALENCE; l++)
	vert[l] = 0;
      vertinlgr (vert, v + 1, vertexes[v].subst, &vertexes[v].lgrnptr);
    }
}

void 
findReversVert (void)
{
  int v, l, lp, k, vin, lin, fin, fout, fl_len;
  for (v = 1; v <= vcs.sizet; v++)
    vertexes[v - 1].r_vert = FALSE;
  for (lp = 1; lp <= nloop; lp++)
    {
      fl_len = fermloops[lp - 1].len;
      for (k = 1; k <= fl_len; k++)
	{
	  l = fermloops[lp - 1].ll[k - 1];
	  v = fermloops[lp - 1].vv[k - 1];
	  if (k > 1)
	    {
	      vin = fermloops[lp - 1].vv[k - 2];
	      lin = fermloops[lp - 1].ll[k - 2];
	    }
	  else
	    {
	      vin = fermloops[lp - 1].vv[fl_len - 1];
	      lin = fermloops[lp - 1].ll[fl_len - 1];
	    }
	  lin = vcs.vertlist[vin - 1][lin - 1].nextvert.edno;
	  fin = 1;
	  while (lin != vertexes[v - 1].subst[fin - 1])
	    fin++;
	  fout = 1;
	  while (l != vertexes[v - 1].subst[fout - 1])
	    fout++;
	  vertexes[v - 1].r_vert = (fin > fout);
	}
    }
}



static void 
findinoutmasses (void)
{
  int v, l, vln;
  int numntot = getntot ();

  for (v = 1; v <= vcs.sizet; v++)
    {
      vln = vcs.valence[v - 1];
      for (l = 1; l <= vln; l++)
	if (vcs.vertlist[v - 1][l - 1].moment > 0 &&
	    vcs.vertlist[v - 1][l - 1].moment <= numntot)
	  strcpy (inoutmasses[vcs.vertlist[v - 1][l - 1].moment - 1],
		  prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].massidnt);
    }
}



static void 
findmassindex (void)
{
  int v, l;

  setofb_zero (setmassindex0);
  for (v = 1; v <= vcs.sizet; v++)
    for (l = 1; l <= vcs.valence[v - 1]; l++)
      if (vcs.vertlist[v - 1][l - 1].lorentz != 0 &&
	  vcs.vertlist[v - 1][l - 1].nextvert.vno < v &&
	  !photonp (vcs.vertlist[v - 1][l - 1].partcl) &&
	  !gaugep (vcs.vertlist[v - 1][l - 1].partcl))
	{
	  setofb_cpy (setmassindex0,
		      setofb_uni (setmassindex0,
			  setofb (vcs.vertlist[v - 1][l - 1].lorentz, _E)));
	  massindpos[vcs.vertlist[v - 1][l - 1].lorentz - 1].vrt1 = v;
	  massindpos[vcs.vertlist[v - 1][l - 1].lorentz - 1].ln1 = l;
	  massindpos[vcs.vertlist[v - 1][l - 1].lorentz - 1].vrt2 =
	    vcs.vertlist[v - 1][l - 1].nextvert.vno;
	  massindpos[vcs.vertlist[v - 1][l - 1].lorentz - 1].ln2 =
	    vcs.vertlist[v - 1][l - 1].nextvert.edno;
	}
}


static int 
vectorslot (int v)
{
  int l;

  l = vcs.valence[v - 1];
  while (l > 0 && vcs.vertlist[v - 1][l - 1].lorentz == 0)
    --(l);
  return l;
}


static int 
vectorslot2 (int v)
{
  int l;

  l = vcs.valence[v - 1];
  while (l > 0 && vcs.vertlist[v - 1][l - 1].lorentz == 0)
    l--;
  if (l != 0)
    {
      l--;
      while (l > 0 && vcs.vertlist[v - 1][l - 1].lorentz == 0)
	l--;
    }
  return l;
}



static void 
nextFerm (int *v_, int *l_)
{
  int l1;
  l1 = vcs.vertlist[(*v_) - 1][(*l_) - 1].nextvert.edno;
  *v_ = vcs.vertlist[(*v_) - 1][(*l_) - 1].nextvert.vno;
  *l_ = 1;
  while ((*l_ == l1) ||
	 (prtclbase[vcs.vertlist[(*v_) - 1][(*l_) - 1].partcl - 1].spin != 1)
    )
    (*l_)++;
}


static void 
findfermcycles (void)
{
  int v, v1, l, l1;
  int count;

  nloop = 0;
  for (v = 1; v <= vcs.sizet; v++)
    fermmap[v - 1] = 0;

  for (v = 1; v <= vcs.sizet; v++)
    {
      if (fermmap[v - 1] == 0)
	{
	  l = vcs.valence[v - 1];
	  while ((l > 0) && (!a_fermionp (vcs.vertlist[v - 1][l - 1].partcl)))
	    l--;
	  if (l != 0)
	    {
	      count = 0;
	      ++(nloop);
	      v1 = v;
	      fermloops[nloop - 1].lprtcl = FALSE;
	      do
		{		/*  Until FermMap[v1]<>0  */
		  ++(count);
		  fermloops[nloop - 1].vv[count - 1] = v1;
		  fermloops[nloop - 1].ll[count - 1] = l;
		  if (strchr ("LR",
		      prtclbase[vcs.vertlist[v1 - 1][l - 1].partcl - 1].hlp)
		      != NULL)
		    fermloops[nloop - 1].lprtcl = TRUE;
		  fermloops[nloop - 1].intln[count - 1] = 0;
		  fermloops[nloop - 1].intln2[count - 1] = 0;
		  fermmap[v1 - 1] = nloop;
		  l1 = vectorslot (v1);
		  if (l1 != 0)
		    {
		      if (fermmap[vcs.vertlist[v1 - 1][l1 - 1].nextvert.vno - 1] == nloop
			)
			fermloops[nloop - 1].intln[count - 1] = l1;
		      l1 = vectorslot2 (v1);
		      if ((l1 != 0) &&
			  (fermmap[vcs.vertlist[v1 - 1][l1 - 1].nextvert.vno - 1] == nloop)
			)
			{
			  if (fermloops[nloop - 1].intln[count - 1] == 0)
			    fermloops[nloop - 1].intln[count - 1] = l1;
			  else
			    fermloops[nloop - 1].intln2[count - 1] = l1;
			}
		    }
		  nextFerm (&v1, &l);

		}
	      while (fermmap[v1 - 1] == 0);
	      fermloops[nloop - 1].g5 = fermloops[nloop - 1].lprtcl;

	      fermloops[nloop - 1].len = count;
	    }
	}
    }
}


static void 
findinnerverts (void)
{
  int v, l, nl;

  for (v = 1; v <= nloop; v++)
    strcpy (fermloops[v - 1].invrt, "");
  for (v = 1; v <= vcs.sizet; v++)
    if (fermmap[v - 1] == 0)
      {
	l = vectorslot (v);
	if (l != 0)
	  {
	    nl = fermmap[vcs.vertlist[v - 1][l - 1].nextvert.vno - 1];
	    if (nl != 0)
	      {
		for (l = l - 1; l >= 1; l--)
		  if (nl != fermmap[vcs.vertlist[v - 1][l - 1].nextvert.vno - 1])
		    goto label_1;
		sbld (fermloops[nl - 1].invrt,
		      "%s%c", fermloops[nl - 1].invrt, v);

	      label_1:;
	      }
	  }
      }
}


static void 
findsubst (int v, int l, char *subst)
{
  momsum frontsubst;
  int i, j, vv, ll;

  subst[0] = 0;			/* strcpy(subst,""); */
  if (				/*(setof(inp,intrp,_E) & vcs.vertlist[v-1][l-1].prop) != setof(_E) */

       vcs.vertlist[v - 1][l - 1].prop & (IN_PRTCL | OUT_PRTCL)

    )
/*      subst[(int)(++subst[0])] = (cahe)(vcs.vertlist[v-1][l-1].moment); */
    subst[(int) (++subst[0])] = (char) (vcs.vertlist[v - 1][l - 1].moment);
  /* sbld(subst,"%s%c",subst,vcs.vertlist[v-1][l-1].moment); */
  else
    for (i = 1; i <= vcs.valence[v - 1]; i++)
      if (i != l)
	{
	  vv = vcs.vertlist[v - 1][i - 1].nextvert.vno;
	  ll = vcs.vertlist[v - 1][i - 1].nextvert.edno;
	  findsubst (vv, ll, frontsubst);
	  for (j = 1; j <= frontsubst[0]; j++)
	    subst[subst[0] + j] = frontsubst[j];
	  subst[0] += frontsubst[0];
	  /* sbld(subst,"%s%s",subst,frontsubst); */
	}
}				/*  FindCond  */


static void 
changesign (char *subst)
{
  int i;

  for (i = 1; i <= /*strlen(subst) */ subst[0]; i++)
    subst[i] = -subst[i];
}


static void 
optimsubst (int v, int l, char *subst)
{
  momsum frontsubst, backsubst;
  int i, vv, ll;

  findsubst (v, l, frontsubst);
  vv = vcs.vertlist[v - 1][l - 1].nextvert.vno;
  ll = vcs.vertlist[v - 1][l - 1].nextvert.edno;
  findsubst (vv, ll, backsubst);
  if (				/*strlen(frontsubst) <= strlen(backsubst) */
       frontsubst[0] <= backsubst[0])
    for (i = 0; i <= frontsubst[0]; i++)
      subst[i] = frontsubst[i];
  /* strcpy(subst,frontsubst); */
  else
    {
      for (i = 0; i <= backsubst[0]; i++)
	subst[i] = backsubst[i];
      /* strcpy(subst,backsubst); */
      changesign (subst);
    }
}				/*  OptimCond  */


static void 
standartsubst (int v, int l, char *subst)
{
  int i, vv, ll;
  int flg1 = FALSE, flg2 = FALSE;
  int ch;
  int numntot = getntot ();
  int numnin = getnin ();

  if (vcs.vertlist[v - 1][l - 1].moment == numntot)
    {
      subst[0] = 0;
      /* strcpy(subst,""); */
      for (i = 1; i <= numnin; i++)
	subst[(int) (++subst[0])] = (char) i;
      /* sbld(subst,"%s%c",subst,i); */
      for (i = numnin + 1; i <= numntot - 1; i++)
	subst[(int) (++subst[0])] = (char) (-i);
      /* sbld(subst,"%s%c",subst,-i); */
      return;
    }

  ch = numntot;
  findsubst (v, l, subst);
  for (i = 1; i <= subst[0]; i++)
    {
      if (subst[i] == ch)
	flg1 = TRUE;
      if (subst[i] == -ch)
	flg2 = TRUE;
    }
  if (flg1 && flg2)
    {
      vv = vcs.vertlist[v - 1][l - 1].nextvert.vno;
      ll = vcs.vertlist[v - 1][l - 1].nextvert.edno;
      findsubst (vv, ll, subst);
      changesign (subst);
    }
}				/*  StandartSubst  */


static void 
findinternalmoments (int indep)
{
  int l, v, vln;
  int m;

  for (m = 1; m <= 3 * maxvert; m++)
    momdep[m - 1][0] = 0;

  for (v = 1; v <= vcs.sizet; v++)
    {
      vln = vcs.valence[v - 1];
      for (l = 1; l <= vln; l++)
	{
	  m = vcs.vertlist[v - 1][l - 1].moment;
	  if (m > 0)
	    {
	      if (indep)
		standartsubst (v, l, momdep[m - 1]);
	      else
		optimsubst (v, l, momdep[m - 1]);
	    }
	}
    }
}


void 
preperdiagram (void)
{
  findinoutmasses ();
  findmassindex ();
  findfermcycles ();
  findinnerverts ();
  findinternalmoments (consLow);
}
