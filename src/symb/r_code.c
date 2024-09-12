/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/parser.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/getmem.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"

#include "physics.h"
#include "ghosts.h"
#include "cweight.h"
#include "prepdiag.h"
#include "reader0.h"
#include "diaprins.h"
#include "chess.h"
#include "out_service.h"
#include "rfactor.h"
#include "denominators.h"
#include "process.h"
#include "process_core.h"
#include "r_code.h"

static int vertmap[3 * maxvert];
static int gammaflag;
static char *parsedverts[2 * maxvert];
static int numberg5;


static char *
pexpr (char p)
{
  static char snum[8];
  if (p >= 0)
    sprintf (snum, "P%d", p);
  else
    sprintf (snum, "(-P%d)", -p);
  return snum;
}				/*  Pexpr  */


static char *
iexpr (int i)
{
  static char snum[6];

  sbld (snum, "m%d", i);
  return snum;
}				/*  iexpr  */



static void 
head (void)
{
  int l;
  int numtot = getntot ();

  writeF ("%% ----------- VARIABLES ------------------ \n");
  writeF (" vector  A,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ZERO_;\n");
  writeF (" vector  m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16;\n");
  writeF ("%%\n");
  writeF ("%%--------- Mass shell declarations -----------\n");
  for (l = 1; l <= numtot; l++)
    writeF (" MASS  P%d = %s$  MSHELL P%d$\n", l, inoutmasses[l - 1], l);
  writeF ("%%\n");
  writeF ("operator propDen$\n");
}


static void 
writesubst (void)
{
  int i, j, l;
  int c;

  writeF (
	   "%%-------- Moment substitutions --------\n");
  for (i = 1; i <= 3 * maxvert; i++)
    {
      l = momdep[i - 1][0];	/* For momdep this is a length V.E. */
      if (l > 1)
	{
	  writeF (" Let  p%d = ", i);
	  for (j = 1; j <= l; j++)
	    {
	      c = /* (char) */ momdep[i - 1][j];
	      if (c > 0)
		writeF ("+p%d", c);
	      else
		writeF ("-p%d", -c);
	    }
	  writeF ("$\n");
	}
    }
  writeF (" Let Sqrt2=sqrt(2)$\n");
  writeF ("%%\n");
}

static void 
emitfactors (void)
{
  char num[STRSIZ];
  int c;
  int i, j, s, vln, pnum;

  writeF ("%%---------- Factors ---------------\n");
  /* ------- Symmetry Factor ------- */
  writeF (
	   " SymmFact:=%d/%d$    %% Diagram symmerty factor\n",
	   vcs.symnum, vcs.symdenum);
  /* -----  average factor  --------- */
  c = 1;
  for (i = 1; i <= vcs.sizel; i++)
    {
      vln = vcs.valence[i - 1];
      for (j = 1; j <= vln; j++)
	if (IN_PRTCL & vcs.vertlist[i - 1][j - 1].prop)
	  {
	    pnum = vcs.vertlist[i - 1][j - 1].partcl;
	    s = prtclbase[pnum - 1].spin;
	    switch (s)
	      {
	      case 1:
		if (strchr ("LR", prtclbase[pnum - 1].hlp) == NULL)
		  c *= 2;
		break;

	      case 2:
		if (zeromass (pnum))
		  c *= 2;
		else
		  c *= 3;
	      }			/*  Case  */
	    c *= abs (prtclbase[pnum - 1].cdim);
	  }
    }
  writeF (
	   " AverFact:=1/%d$       %% Normalization factor of polarization average\n", c);
  /* ----- Fermion factor  --------- */
  c = 0;

  for (i = 0; i < vcs.sizel; i++)
    for (j = 0; j < vcs.valence[i]; j++)
      {
	pnum = vcs.vertlist[i][j].partcl;
	if (prtclbase[pnum - 1].spin == 1 && (IN_PRTCL & vcs.vertlist[i][j].prop))
	  ++(c);
      }


  if ((c & 1) == 1)
    strcpy (num, "-1");
  else
    strcpy (num, "1");
  writeF (
    " FermFact:=%s$      %% (-1)**(number of in-fermion particles)\n", num);
  /* ----- Vector factor  --------- */
/*   c = 0;
   for (i = 1; i <= vcs.sizet; i++)
   {
   vln = vcs.valence[i-1];
   for (j = 1; j <= vln; j++)
   if (i < vcs.vertlist[i-1][j-1].nextvert.vno)
   {
   pnum = vcs.vertlist[i-1][j-1].partcl;
   if (prtclbase[pnum-1].spin == 2) ++(c);
   }
   }
   if ((c & 1) == 1)
   strcpy(num,"-1");
   else
   strcpy(num,"1");
   writeF(
   " VectFact:=%s$      %% (-1)**(number of vector propagators)\n",num);
 */
  /* ----- Color factor  --------- */
  writeF (
	   " ColorFact:=%d/%d$    %%  QCD color weight of diagram  \n",
	   vcs.clrnum, vcs.clrdenum);
  writeF ("%%\n");
}


 /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    | Program for emitting denominators of propagators in CrossSection |
    |                        diagramm.                                 |
    |                September 19, 1989                                |
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

static void 
emitdenoms (void)
{
  int j;			/* 1..2 * maxvert */
  int k;			/* 1..4 */
  int denrno;

  /*  EmitDenoms -- main  */

  for (j = 1; j <= vcs.sizet; j++)
    for (k = 1; k <= vcs.valence[j - 1]; k++)
      {
	edgeinvert *ln = &vcs.vertlist[j - 1][k - 1];

	if (ln->moment > 0 && !(ln->prop & (IN_PRTCL | OUT_PRTCL)) && pseudop (ln->partcl))
	  writeF ("totFactor_:=totFactor_/%s**2$\n",
		  prtclbase[ln->partcl - 1].massidnt);

      }

  writeF ("denominator_:=");
  denrno = 0;
  for (j = 1; j <= vcs.sizet; j++)
    for (k = 1; k <= vcs.valence[j - 1]; k++)
      {
	edgeinvert *ln = &vcs.vertlist[j - 1][k - 1];

	if (ln->moment > 0 && !(ln->prop & (IN_PRTCL | OUT_PRTCL)))
	  {
	    if (!pseudop (ln->partcl))
	      {
		char width[10] = "0";
		if (denrno)
		  writeF ("*");
		denrno++;
		if (!ttypepropag (j, k))
		  strcpy (width, prtclbase[ln->partcl - 1].imassidnt);
		writeF ("propDen(p%d,%s,%s)", ln->moment,
			prtclbase[ln->partcl - 1].massidnt, width);
	      }
	  }
      }
  if (!denrno)
    writeF ("1");
  writeF ("$\n");

}				/*   EmitDenoms   */

/*  Programm for constructing vertex-oriented representation for !
   ! amplitude diagram.                                             !
   ! A.Taranov 20.07.89                                            */


static void 
emitindex (indvertset index)
{
  int m = 1, n = 0, l = 1;

  if (setofb_eq0 (index))
    return;
  writeF (" Index ");
  do
    {
      if (insetb (m, index))
	{
	  if (n)
	    writeF (",");
	  n++;
	  writeF ("m%d", l);
	  setofb_cpy (index, setofb_aun (index, setofb (m, _E)));
	}
      l++;
      if (!setofb_eq0 (index))
	m++;
    }
  while (!(n == 1 /*3 */  || setofb_eq0 (index)));
  writeF ("$\n");
}				/*  EmitIndex  */


static void 
reducemult (char *nameres, char *name1, char *name2,
	    indvertset index)
{
  indvertset index1, index2;
  int l, n;

  lvcpy (index1, index);
  lvcpy (index2, index);
  emitindex (index1);
  writeF (" %s:=%s*%s$\n", nameres, name1, name2);
  while (!setofb_eq0 (index1))
    {
      emitindex (index1);
      writeF (" %s:=%s$\n", nameres, nameres);
    }
  if (!setofb_eq0 (index2))
    {
      writeF (" RemInd  ");
      l = 1;
      n = 1;
      while (!setofb_eq0 (index2))
	{
	  if (insetb (l, index2))
	    {
	      if (n != 1)
		writeF (",");
	      ++(n);
	      writeF ("m%d", l);
	    }
	  setofb_cpy (index2, setofb_aun (index2, setofb (l, _E)));
	  ++(l);
	}
      writeF ("$\n");
    }
}


static int 
fermvrt (int v)
{
  int i;

  for (i = 0; i < vcs.valence[v - 1]; i++)
    if (fermionp (vcs.vertlist[v - 1][i].partcl))
      return 1;
  return 0;
}

static char *
subst (int ind, int arg)
{
  int v, l, i;
  static shortstr ans;

  v = massindpos[ind - 1].vrt1;
  if (vertmap[v - 1] == arg)
    l = massindpos[ind - 1].ln1;
  else
    {
      v = massindpos[ind - 1].vrt2;
      l = massindpos[ind - 1].ln2;
    }

  sbld (ans, "%s=>(", iexpr (ind));
  if (gammaflag && fermvrt (v))
    {
      for (i = 1; i <= vcs.valence[v - 1]; i++)
	if (i != l)
	  sbld (ans, "%s+%s", ans, pexpr (vcs.vertlist[v - 1][i - 1].moment));
    }
  else
    sbld (ans, "%s%s", ans, pexpr (-vcs.vertlist[v - 1][l - 1].moment));
  sbld (ans, "%s)/%s", ans, prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].massidnt);
  return ans;
}


static void 
r_mult (int arg1, int arg2, indvertset index1)
{
  char name0[8], name1[8], name2[8];
  indvertset mind;
  unsigned m, mm, maxmm;
  int pstn[15];
  int n, l;
  indvertset index;
  /*  R_mult  */
  setofb_cpy (index, index1);
  sbld (name1, "Vrt_%d", arg1);
  sbld (name2, "Vrt_%d", arg2);
  if (arg1 < arg2)
    strcpy (name0, name1);
  else
    strcpy (name0, name2);
  setofb_cpy (mind, setofb_its (index, setmassindex));
  if (setofb_eq0 (mind))
    reducemult (name0, name1, name2, index);
  else
    {
      reducemult ("Vrt_0", name1, name2, index);
      l = 1;
      n = 0;
      while (!setofb_eq0 (mind))
	{
	  if (insetb (l, mind))
	    pstn[++n - 1] = l;
	  setofb_cpy (mind, setofb_aun (mind, setofb (l, _E)));
	  ++(l);
	}
      maxmm = (1 << n) - 1;
      for (mm = 1; mm <= maxmm; mm++)
	{
	  writeF (
		   " Vrt_L:=%s$   Vrt_R:=%s$\n", name1, name2);
	  setofb_zero (mind);
	  m = mm;
	  n = 1;
	  while (m != 0)
	    {
	      if ((m & 1) != 0 && (l = pstn[n - 1]) != 0)
		{
		  writeF (" Vrt_L:=(Vrt_L where %s)$\n", subst (l, arg1));
		  writeF (" Vrt_R:=(Vrt_R where %s)$\n", subst (l, arg2));
		  setofb_cpy (mind, setofb_uni (mind, setofb (l, _E)));
		}
	      ++(n);
	      m >>= 1;
	    }
	  reducemult ("Vrt_0", "Vrt_0 + Vrt_L", "Vrt_R",
		      setofb_aun (index, mind));
	  writeF (
		   " Clear Vrt_L $    Clear Vrt_R $\n");
	}
      writeF (
	       " %s:=Vrt_0 $    Clear Vrt_0 $\n", name0);
    }
  if (arg1 < arg2)
    {
      writeF (" Clear %s$\n", name2);
      for (n = 1; n <= vcs.sizet; n++)
	if (vertmap[n - 1] == arg2)
	  vertmap[n - 1] = arg1;
    }
  if (arg2 < arg1)
    {
      writeF (" Clear %s$\n", name1);
      for (n = 1; n <= vcs.sizet; n++)
	if (vertmap[n - 1] == arg1)
	  vertmap[n - 1] = arg2;
    }
}				/*  R_mult  */

/*  Emit Reduce prg for ferm loop calculation  */

static char *
fermpropag (int v, int l)
{
  static shortstr ans, mass;
  edgeinvert *ln = &vcs.vertlist[v - 1][l - 1];

  strcpy (mass, prtclbase[ln->partcl - 1].massidnt);
  if (prtclbase[ln->partcl - 1].hlp == '*')
    return mass;

  sprintf (ans, "G(ln,%s)", pexpr (ln->moment));
  if (strcmp (mass, "0") != 0)
    sbld (ans, "(%s+%s)", ans, mass);
  else
    {
      if (fermionp (ln->partcl))
	{
	  if (prtclbase[ln->partcl - 1].hlp == 'L')
	    strcat (ans, "*(1-G(ln,A))/2");
	  else if (prtclbase[ln->partcl - 1].hlp == 'R')
	    strcat (ans, "*(1+G(ln,A))/2");
	}
      else
	{
	  if (prtclbase[ln->partcl - 1].hlp == 'R')
	    strcat (ans, "*(1-G(ln,A))/2");
	  else if (prtclbase[ln->partcl - 1].hlp == 'L')
	    strcat (ans, "*(1+G(ln,A))/2");
	}
    }
  return ans;
}


static void 
fermprgemit (void)
{
  int v, v1, l;
  int i, k, lpcount, nn;
  indvertset indexs;
  char *vertStr, *propStr;

  /* Nested function: fermpropag */
  /* Main Procedure -- FermPrgEmit  */

  numberg5 = 0;
  if (nloop == 0)
    {
      gammaflag = FALSE;
      return;
    }
  gammaflag = TRUE;

  for (v = 1; v <= nloop; v++)
    {
      l = 1;
      fermloops[v - 1].g5 = fermloops[v - 1].lprtcl;
      while (!fermloops[v - 1].g5 && l <= fermloops[v - 1].len)
	{
	  fermloops[v - 1].g5 =
	    spos ("G(ln,A)", parsedverts[fermloops[v - 1].vv[l - 1] - 1]) != 0;
	  ++(l);
	}
      if (fermloops[v - 1].g5)
	++(numberg5);
    }

  for (v = 1; v <= 2 * maxvert; v++)
    vertmap[v - 1] = 0;

  for (lpcount = 1; lpcount <= nloop; lpcount++)
    {
      writeF (
	       "%%------- Fermion loop calculation ------- \n");
      writeF (" NoSpur ln$\n");
      for (nn = 1; nn <= fermloops[lpcount - 1].len; nn++)
	{
	  if (vertexes[fermloops[lpcount - 1].vv[nn - 1] - 1].r_vert)
	    writeF ("%% reversed vertex\n");

	  vertStr = parsedverts[fermloops[lpcount - 1].vv[nn - 1] - 1];
	  propStr = fermpropag (fermloops[lpcount - 1].vv[nn - 1],
				fermloops[lpcount - 1].ll[nn - 1]);
	  writeF (" Vrt_%d:=%s*%s$\n", nn + lpcount - 1, vertStr, propStr);

	  vertmap[fermloops[lpcount - 1].vv[nn - 1] - 1] = nn + lpcount - 1;
	}
/*======================
		if (fermloops[lpcount-1].intln[fermloops[lpcount-1].len-1] == 0)
         setofb_zero(indexs);
		else
			setofb_cpy (indexs,setofb(vcs.vertlist[fermloops[lpcount-1].vv[fermloops[lpcount-1].len-1]-1]
			[fermloops[lpcount-1].
			intln[fermloops[lpcount-1].len-1]-1].lorentz ,_E));
========================*/

      for (v1 = 2; v1 <= fermloops[lpcount - 1].len; v1++)
	{
	  int newInd = fermloops[lpcount - 1].intln[v1 - 1];	/* vector line */
	  if (newInd != 0)
	    newInd =
	      vcs.vertlist[fermloops[lpcount - 1].vv[v1 - 1] - 1][newInd - 1].lorentz;

	  if (newInd == 0)
	    setofb_zero (indexs);
	  else
	    setofb_cpy (indexs, setofb (newInd, _E));
	  newInd = fermloops[lpcount - 1].intln2[v1 - 1];	/* second vector line */
	  if (newInd != 0)
	    newInd =
	      vcs.vertlist[fermloops[lpcount - 1].vv[v1 - 1] - 1][newInd - 1].lorentz;

	  if (newInd != 0)
	    setofb_cpy (indexs, setofb_uni (indexs,
					    setofb (newInd, _E)));

	  r_mult (lpcount, lpcount + v1 - 1, indexs);
	}

      if (numberg5 == 1 && fermloops[lpcount - 1].g5)
	{
	  writeF (
		   " Vrt_%d:=Sub(A=0*A,Vrt_%d)$\n"
	  /* " Vrt_%d:=(Vrt_%d where A=>0*ZERO_)$\n" */
		   ,lpcount, lpcount);
	  numberg5 = 0;
	  fermloops[lpcount - 1].g5 = FALSE;
	}

      for (k = 1; k <= strlen (fermloops[lpcount - 1].invrt); k++)
	{
	  v1 = fermloops[lpcount - 1].invrt[k - 1];
	  setofb_zero (indexs);
	  for (l = 1; l <= vcs.valence[v1 - 1]; l++)
	    {
	      i = vcs.vertlist[v1 - 1][l - 1].lorentz;
	      if (i != 0)
		setofb_cpy (indexs, setofb_uni (indexs, setofb (i, _E)));
	    }
	  vertmap[v1 - 1] = lpcount + 1;
	  writeF ("Vrt_%d:=%s$\n", lpcount + 1, parsedverts[v1 - 1]);
	  r_mult (lpcount, lpcount + 1, indexs);
	}
      writeF (" Spur ln $\n");
      writeF (
	       " Vrt_%d:=-4*Vrt_%d$\n", lpcount, lpcount);
      writeF ("%%\n");
      writeF (
	       " Fl%d:=Vrt_%d$\n", lpcount, lpcount);
      writeF ("%%\n");
    }

  for (v1 = 1; v1 <= 2 * maxvert; v1++)
    fermmap[v1 - 1] = vertmap[v1 - 1];
  writeF ("%%\n");
  gammaflag = FALSE;
}

static void 
formblocks (int *vrtsize)
{
  int v, l, lori, v1, vv, v2, count;
  char vnum[11];

  /*    FormBlocks   */
  for (v = 1; v <= 2 * maxvert; v++)
    vertmap[v - 1] = fermmap[v - 1];	/*  copy fermMap      */
  for (count = 1; count <= nloop; count++)	/*  ferm loop vertex  */
    {
      sbld (vnum, "%d", count);
      writeF (" Vrt_%s:=Fl%s$", vnum, vnum);
      writeF (" Clear Fl%s$\n", vnum);
    }
  count = nloop;
  for (v = 1; v <= vcs.sizet; v++)	/*  simple  vertex  */
    if (vertmap[v - 1] == 0)
      {
	++(count);
	vertmap[v - 1] = count;
	writeF (" Vrt_%d:=%s$\n", count, parsedverts[v - 1]);
      }
  /*  end of vertex filling  */

  /*    Blk.size:=Count;        *//*   Begin of Blk filling  */
  *vrtsize = count;
  for (v = 1; v <= count; v++)
    {
      vertinfo[v - 1].vlnc = 0;
      vertinfo[v - 1].weit = 1;
      setofb_zero (vertinfo[v - 1].ind);
      vertinfo[v - 1].g5 = 0;
    }

  for (v = 1; v <= vcs.sizet; v++)
    {
      v1 = vertmap[v - 1];
      if (v1 <= nloop)
	{
	  vertinfo[v1 - 1].weit += 2;
	  if (fermloops[v1 - 1].g5)
	    {
	      ++(vertinfo[v1 - 1].weit);
	      vertinfo[v1 - 1].g5 = 1;
	    }
	}
      for (l = 1; l <= vcs.valence[v - 1]; l++)
	{
	  vv = vcs.vertlist[v - 1][l - 1].nextvert.vno;
	  v2 = vertmap[vv - 1];
	  lori = vcs.vertlist[v - 1][l - 1].lorentz;
	  if (lori != 0)
	    {
	      if (v1 != v2)
		{
		  int np;
		  ++(vertinfo[v1 - 1].vlnc);
		  vertinfo[v1 - 1].link[vertinfo[v1 - 1].vlnc - 1] = v2;
		  setofb_cpy (vertinfo[v1 - 1].ind,
		      setofb_uni (vertinfo[v1 - 1].ind, setofb (lori, _E)));
		  np = vcs.vertlist[v - 1][l - 1].partcl;
		  if (prtclbase[np - 1].hlp == 't')
		    {
		      ++(vertinfo[v1 - 1].vlnc);
		      vertinfo[v1 - 1].link[vertinfo[v1 - 1].vlnc - 1] = v2;
		      setofb_cpy (vertinfo[v1 - 1].ind,
				  setofb_uni (vertinfo[v1 - 1].ind, setofb (lori - 1, _E)));
		    }
		}
	      vertinfo[v1 - 1].weit += 2;
	      if (insetb (lori, setmassindex))
		++(vertinfo[v1 - 1].weit);
	    }
	}
    }				/*   End of Blk filling  */
}				/*  FormBlocks  */


static void 
emitreducecode (void)
{
  int i, j, v0, v1, v2;
  indvertset ind1, ind2;
  int g5_exist;
  char *pstr;
  int sgn;

  for (i = 1; i <= vcs.sizet; i++)
    {
      for (j = 1; j <= MAX (1, vcs.valence[i - 1]); j++)
	{
	  momsubst[j - 1] =
	    vcs.vertlist[i - 1][vertexes[i - 1].subst[j - 1] - 1].moment;
	  indsubst[j - 1] =
	    vcs.vertlist[i - 1][vertexes[i - 1].subst[j - 1] - 1].lorentz;
	}
      r_reading0 = vertexes[i - 1].r_vert;

      pstr = (char *) readExpression (vertexes[i - 1].lgrnptr->description,
				      rd_rcode, act_rcode, free);
      if (rderrcode != 0)
	{
	  outFileClose ();
	  finish ("");
	  exit (0);
	}

      sgn = 1;
      if (vcs.valence[i - 1] == 3 &&
	  prtclbase[vcs.vertlist[i - 1][0].partcl - 1].cdim == 8 &&
	  prtclbase[vcs.vertlist[i - 1][1].partcl - 1].cdim == 8 &&
	  prtclbase[vcs.vertlist[i - 1][2].partcl - 1].cdim == 8)
	{
	  if (vertexes[i - 1].subst[0] > vertexes[i - 1].subst[1])
	    sgn = -sgn;
	  if (vertexes[i - 1].subst[1] > vertexes[i - 1].subst[2])
	    sgn = -sgn;
	  if (vertexes[i - 1].subst[0] > vertexes[i - 1].subst[2])
	    sgn = -sgn;
	}

      parsedverts[i - 1] = m_alloc (8 + strlen (pstr));

      if (sgn == 1)
	sprintf (parsedverts[i - 1], "(%s)", pstr + 2);
      else
	sprintf (parsedverts[i - 1], "(-1)*(%s)", pstr + 2);
      free (pstr);
    }

  fermprgemit ();
  formblocks (&n_vrt);
  makeprgcode ();
  for (i = 1; i <= vcs.sizet; i++)
    free (parsedverts[i - 1]);


  g5_exist = numberg5 != 0;
  for (i = n_vrt - 1; i >= 1; i--)
    {
      v1 = prgcode[i - 1][0];
      v2 = prgcode[i - 1][1];
      setofb_cpy (ind1, vertinfo[v1 - 1].ind);
      setofb_cpy (ind2, vertinfo[v2 - 1].ind);
      if (g5_exist &&
	  numberg5 == vertinfo[v1 - 1].g5 + vertinfo[v2 - 1].g5)
	{
	  writeF (" Vrt_1r:=(Vrt_%d where i=>0 )$\n", v1);
	  writeF ("Vrt_%d:=Vrt_%d-Vrt_1r $\n", v1, v1);
	  writeF (" Vrt_2r:=(Vrt_%d where i=>0)$\n", v2);
	  writeF ("Vrt_%d:=Vrt_%d-Vrt_2r $\n", v2, v2);
	  r_mult (v1, v2, setofb_its (ind1, ind2));
	  v0 = MIN (v1, v2);
	  writeF ("Vrt_Tmp:=Vrt_%d$\n", v0);
	  writeF (
		   "Vrt_%d:=Vrt_1r $    Clear Vrt_1r $  \n", v1);
	  writeF (
		   "Vrt_%d:=Vrt_2r $    Clear Vrt_2r $  \n", v2);
	  r_mult (v1, v2, setofb_its (ind1, ind2));
	  writeF (
		   "Vrt_%d:=Vrt_%d+Vrt_Tmp $  Clear Vrt_Tmp $\n", v0, v0);
	  g5_exist = FALSE;
	}
      else
	r_mult (v1, v2, setofb_its (ind1, ind2));
      vertinfo[MIN (v1, v2) - 1].g5 = vertinfo[v1 - 1].g5 + vertinfo[v2 - 1].g5;
      setofb_cpy (vertinfo[MIN (v1, v2) - 1].ind,
	     setofb_aun (setofb_uni (ind1, ind2), setofb_its (ind1, ind2)));
    }


  for (i = 0; i < vcs.sizet; i++)
    for (j = 0; j < vcs.valence[i]; j++)
      {
	int mom, np;
	mom = vcs.vertlist[i][j].moment;
	np = vcs.vertlist[i][j].partcl;
	if (mom > 0 && prtclbase[np - 1].hlp == 't')
	  {
	    if (prtclbase[ghostmother (np) - 1].hlp == '*')
	      writeF ("Vrt_1:=Vrt_1*(%s**2)$ \n",
		      prtclbase[np - 1].massidnt);
	    else
	      writeF ("Vrt_1:=Vrt_1*(%s**2-p%d.p%d)$ \n",
		      prtclbase[np - 1].massidnt, mom, mom);
	  }
      }
  writeF ("%%\n");
}


void 
mk_reduceprograms (void)
{
  int ndel, ncalc, nrest;
  long nrecord, naxu;
  csdiagram csd;
  unsigned ncalctot;
  shortstr txt;
  hlpcsptr gstlist, c;
  vcsect vcs_copy;
  int i;
  s_listptr d_facts, df;
  rmptr t_fact;
  int numtot = getntot ();

  goto_xy (1, 21);
  scrcolor (Yellow, Blue);
  print ("  REDUCE code generation \n");
  scrcolor (Red, BGmain);
  print (" Generated........\n");
  print (" current diagram :\n");
  scrcolor (Yellow, Blue);
  print (" Press Esc to halt REDUCE codes generation ");
  scrcolor (FGmain, BGmain);
  diagrq = fopen (DIAGRQ_NAME, "rb");
  ncalctot = 0;
  menuq = fopen (MENUQ_NAME, "rb");

  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      rd_menu (menuq, 2, nsub, txt, &ndel, &ncalc, &nrest, &nrecord);
      fseek (diagrq, nrecord * sizeof (csdiagram), SEEK_SET);
      naxu = ndel + ncalc + nrest;
      for (ndiagr = 1; ndiagr <= naxu; ndiagr++)
	{
	  goto_xy (20, 22);
	  print ("%u", ncalctot);
	  goto_xy (20, 23);
	  print ("%u", ndiagr);
	  clr_eol ();
	  FREAD1 (csd, diagrq);
	  if (csd.status != -1)
	    {

	      outFileOpen (scat ("%sresults%cp%d_%d.red",
				 pathtouser, f_slash, nsub, ndiagr));

	      writeLabel ('%');
	      writeF ("%%\n");

	      transfdiagr (&csd, &vcs);

	      cwtarg (&vcs);
	      if (vcs.clrnum == 0)
		{
		  writeF (
			   "%%-------  Zero color factor --------\n");
		  writeF ("totFactor_:=0$\n");
		  writeF ("numerator_:=0$\n");
		  writeF ("denominator_:=1$\n");
		}
	      else
		{
		  generateghosts (&vcs, &gstlist);
		  if (gstlist == NULL)
		    {
		      writeF ("%%-------  non-existent diagram  --------\n");
		      writeF ("totFactor_:=0$\n");
		      writeF ("numerator_:=0$\n");
		      writeF ("denominator_:=1$\n");
		    }
		  else
		    {
		      goto_xy (40, 23);
		      print ("(%% %4d subdiagrams)", gstlist->maxnum);
		      writeF ("%% The total number of diagrams %d\n", gstlist->maxnum);
		      preperdiagram ();
		      head ();
		      emitfactors ();
		      diagramsrfactors (gstlist, &d_facts, &t_fact);
		      writeF ("totFactor_:=%s$\n", rmonomtxt (*t_fact));

		      writeF ("totFactor_:="
		      "totFactor_*SymmFact*AverFact*FermFact*ColorFact$\n");

		      clrvm (t_fact->n.v);
		      clrvm (t_fact->d.v);
		      free (t_fact);

		      writesubst ();
		      writeF ("numerator_:=0$\n");

		      c = gstlist;
		      df = d_facts;
		      vcs_copy = vcs;
		      while (c != NULL)
			{
			  coloringvcs (c);
			  writeF ("%%  diagram  number =   %d\n", c->num);
			  DiagramToOutFile (&vcs, TRUE, '%');

			  {
			    int k;
			    int sgn = c->sgn;
			    for (k = 0; k < vcs.sizet; k++)
			      sgn *= vertexes[k].lgrnptr->factor;
			    writeF ("  GhostFact:=%d$\n", sgn);
			  }

			  findReversVert ();
			  attachvertexes ();

			  emitreducecode ();

			  writeF (" numerator_:=numerator_ +(%s)*GhostFact*Vrt_1 $\n",
				  smonomtxt (df->monom));
			  writeF (" Clear Vrt_1,GhostFact$\n");
			  writeF ("%%\n");

			  vcs = vcs_copy;
			  c = c->next;
			  df = df->next;
			}

		      eraseslist (d_facts);
		      eraseghosts (gstlist);

		      vcs = vcs_copy;
		      emitdenoms ();
		      writeF (" Clear p%d", numtot + 1);
		      for (i = numtot + 2; i <= 12; i++)
			writeF (",p%d", i);
		      writeF ("$\n");
		      writeF ("%%\n");

		    }
		}
	      writeF ("End$\n");
	      outFileClose ();
	      --(nrest);
	      ++(ncalctot);
	      if (escpressed ())
		goto exi;
	    }
	}
    }

exi:
  fclose (diagrq);
  fclose (menuq);
  clrbox (1, 21, 70, 24);

}

void 
makeghostdiagr (int dnum, char *fname)
{
  csdiagram csd;
  hlpcsptr gstlist, c;
  vcsect vcs_copy;
  FILE *diagrq;
  diagrq = fopen (DIAGRQ_NAME, "rb");
  fseek (diagrq, dnum * sizeof (csdiagram), SEEK_SET);
  FREAD1 (csd, diagrq);
  {
    outFileOpen (fname);
    writeLabel ('%');
    transfdiagr (&csd, &vcs);

    cwtarg (&vcs);
    if (vcs.clrnum == 0)
      writeF ("Color weight equal zero \n");
    else
      {
	generateghosts (&vcs, &gstlist);
	if (gstlist == NULL)
	  {
	    writeF (
		     "Diagrams of this type are absent\n");
	    csd.status = 2;
	    fseek (diagrq, dnum * sizeof (csdiagram), SEEK_SET);
	    FWRITE1 (csd, diagrq);
	  }
	else
	  {
	    emitfactors ();
	    c = gstlist;
	    vcs_copy = vcs;
	    writeF (" total number diagrams of this type is  %u\n", c->maxnum);
	    while (c != NULL)
	      {
		coloringvcs (c);
		{
		  writeF ("  diagrams number =   %u\n", c->num);
		  DiagramToOutFile (&vcs, FALSE, ' ');
		}
		writeF ("  GhostFact:=%d$\n", c->sgn);
		vcs = vcs_copy;
		c = c->next;
	      }

	  }
	eraseghosts (gstlist);
	vcs = vcs_copy;
      }
    writeF ("End$\n");
    outFileClose ();
    if (escpressed ())
      goto exi;
  }

exi:
  fclose (diagrq);
}
