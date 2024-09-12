/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"

#include"saveres.h"
#include"prepdiag.h"
#include"pvars.h"
#include"process.h"
#include"process_core.h"
#include"denominators.h"

static int 
nincount (int v, int l)
{
  int i, vv, ll, summ;

  if (IN_PRTCL & vcs.vertlist[v - 1][l - 1].prop)
    return 1;
  if (OUT_PRTCL & vcs.vertlist[v - 1][l - 1].prop)
    return 0;
  summ = 0;
  vv = vcs.vertlist[v - 1][l - 1].nextvert.vno;
  ll = vcs.vertlist[v - 1][l - 1].nextvert.edno;
  for (i = 1; i <= vcs.valence[vv - 1]; i++)
    if (i != ll)
      summ += nincount (vv, i);
  return summ;
}



int 
ttypepropag (int v, int l)
{
  if (getnin () == 1)
    return FALSE;
  return (nincount (v, l) == 1);
}



/* momdep from prepdiag.h must be calculated before */
void 
calcdenominators (vcsect vcs)
{
  int v, l, k;
  char buff[MAXINOUT + 1];
  int numtot = getntot ();


  denrno = 0;
  for (v = 1; v <= vcs.sizet; v++)
    for (l = 1; l <= vcs.valence[v - 1]; l++)
      {
	edgeinvert *ln = &vcs.vertlist[v - 1][l - 1];
	if (!(ln->moment < 0 || ((IN_PRTCL | OUT_PRTCL) & ln->prop) || pseudop (ln->partcl)))
	  {
	    for (k = 1; k <= momdep[ln->moment - 1][0]; k++)
	      buff[k - 1] = momdep[ln->moment - 1][k];
	    buff[k - 1] = 0;

	    k = -1;
	    while (buff[++k])
	      if (buff[k] < 0)
		buff[k] = -buff[k];
	    if ((2 * strlen (buff) > numtot) ||
		(2 * strlen (buff) == numtot && !strchr (buff, 1))
	      )
	      {
		int ll = 0;
		char buff2[MAXINOUT + 1];
		for (k = 1; k <= numtot; k++)
		  {
		    if (!strchr (buff, k))
		      buff2[ll++] = k;
		  }
		buff2[ll] = 0;
		strcpy (buff, buff2);
	      }
	    k = 0;
	    while (buff[k])
	      {
		if (!k)
		  k++;
		if (buff[k] < buff[k - 1])
		  {
		    int c = buff[k];
		    buff[k] = buff[k - 1];
		    buff[k - 1] = c;
		    k--;
		  }
		else
		  k++;
	      }

	    strcpy (denom[denrno].momStr, buff);
	    denom[denrno].power = 1;

	    denom[denrno].mass = modelVarPos (prtclbase[ln->partcl - 1].massidnt);
	    if (ttypepropag (v, l))
	      denom[denrno].width = 0;
	    else
	      denom[denrno].width = modelVarPos (prtclbase[ln->partcl - 1].imassidnt);

	    for (k = 0; k < denrno; k++)
	      if (!strcmp (denom[denrno].momStr, denom[k].momStr) &&
		  denom[denrno].mass == denom[k].mass &&
		  denom[denrno].width == denom[k].width)
		{
		  ++(denom[k].power);
		  goto label_1;
		}
	    denrno++;
	  label_1:;
	  }
      }
}


void 
denominatorStatistic (FILE * fres, int nsub,
	 int *n_width, int *n_0_width, denlist * allDenominators, FILE * fd)
{
  int i;
  catrec cr;
  denlist den_, den_tmp;
  deninforec dendescript;

  (*n_width) = 0;
  (*n_0_width) = 0;

  den_ = NULL;

  fseek (catalog, 0, SEEK_SET);
  while (FREAD1 (cr, catalog))
    {
      if (cr.nsub_ == nsub)
	{
	  dendescript.cr_pos = ftell (catalog) - sizeof (cr);

	  fseek (fres, cr.denompos, SEEK_SET);
	  readDenominators (fres);
	  dendescript.tot_den = denrno;

	  for (i = 0; i < dendescript.tot_den; i++)
	    {
	      dendescript.denarr[i].power = denom[i].power;
	      dendescript.denarr[i].width = denom[i].width;
	      den_tmp = den_;
	      while (den_tmp != NULL &&
		     (strcmp (denom[i].momStr, den_tmp->momStr)
		      || denom[i].mass != den_tmp->mass
		      || denom[i].width != den_tmp->width))
		den_tmp = den_tmp->next;
	      if (den_tmp == NULL)
		{
		  den_tmp = (denlist) getmem_ ((unsigned) sizeof (denlistrec));
		  den_tmp->next = den_;
		  strcpy (den_tmp->momStr, denom[i].momStr);
		  den_tmp->mass = denom[i].mass;
		  den_tmp->width = denom[i].width;
		  den_ = den_tmp;
		  if (denom[i].width)
		    den_tmp->order_num = ++(*n_width);
		  else
		    den_tmp->order_num = ++(*n_0_width);
		}
	      dendescript.denarr[i].order_num = den_tmp->order_num;
	    }
	  FWRITE1 (dendescript, fd);
	}			/* if CR.nsub_ =nsub */
    }
  *allDenominators = den_;
}


void 
readDenominators (FILE * fres)
{
  int i, m;

  FREAD1 (denrno, fres);	/*  number of demominatirs  */
  for (i = 0; i < denrno; i++)
    {
      FREAD1 (denom[i].power, fres);	/*  power  1 or 2  */
      FREAD1 (denom[i].mass, fres);
      FREAD1 (denom[i].width, fres);
      m = 0;
      do
	FREAD1 (denom[i].momStr[m], fres);
      while (denom[i].momStr[m++]);
    }
}
