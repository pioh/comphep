/*
* Copyright (C) 2003-2009, CompHEP Collaboration
* Copyright (C) 2003, Slava Bunichev
* ------------------------------------------------------
*/
#include<sys/dir.h>
#include<stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/unix_utils.h"
#include "service2/include/read_func.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "service2/include/parser.h"
#include "service2/include/getmem.h"
#include "chep_crt/include/chep_crt.h"

#include "physics.h"
#include "ghosts.h"
#include "cweight.h"
#include "prepdiag.h"
#include "diaprins.h"
#include "chess.h"
#include "out_service.h"
#include "rfactor.h"
#include "denominators.h"
#include "process.h"
#include "process_core.h"
#include "r_code.h"
#include "procvar.h"
#include "out_c.h"
#include "form_code.h"

char momsubst[9] = { 0 };
char indsubst[9] = { 0 };

static int r_reading0 = FALSE;
static int vertmap[3 * maxvert];
static int gammaflag;
static char *parsedverts[2 * maxvert];
static int numberg5;
static int indicatorFl[20];
static int FermLoopIndicator;


void
run_formprograms (void)
{
  int i;
  DIR *dirp;
  struct direct *directp;
  char form_filename[25];
  char cmd_string[50];

  dirp = opendir ("results/");
  while ((directp = readdir (dirp)) != NULL)
    {
      sscanf (directp->d_name, "%s", form_filename);
      for (i = 0; i < 25; i++)
	if (form_filename[i] == '.')
	  {
	    if ((form_filename[i + 1] == 'f') && (form_filename[i + 2] == 'r')
		&& (form_filename[i + 3] == 'm'))
	      {
		sprintf (cmd_string, "./form ./results/%s\n", form_filename);
		system (cmd_string);
	      }
	  }
    }

  (void) closedir (dirp);
}


static void *
bact0 (char ch, pointer mm1, pointer mm2)
{
  char *m1, *m2, *ans;
  int sgn;

  if (r_reading0 && (ch == '*'))
    {
      m1 = (char *) mm2;
      m2 = (char *) mm1;
    }
  else
    {
      m1 = (char *) mm1;
      m2 = (char *) mm2;
    }
  if (ch == '+' || ch == '-')
    {
      lShift (m1, 2);
      lShift (m2, 2);
    }
  else
    {
      if (m1[0] == 'P' || ch == '^')
	{
	  lShift (m1, 1);
	  m1[0] = '(';
	  strcat (m1, ")");
	}
      else
	lShift (m1, 2);
      if (m2[0] == 'P' || ch == '^')
	{
	  lShift (m2, 1);
	  m2[0] = '(';
	  strcat (m2, ")");
	}
      else
	lShift (m2, 2);
    }

  ans = (char *) m_alloc (strlen (m1) + strlen (m2) + 8);
  switch (ch)
    {
    case '+':
      if (m2[0] == '-')
	sprintf (ans, "P|%s%s", m1, m2);
      else
	sprintf (ans, "P|%s+%s", m1, m2);
      break;

    case '-':
      if (m2[0] == '-')
	sprintf (ans, "P|%s+%s", m1, m2 + 1);
      else
	sprintf (ans, "P|%s-%s", m1, m2);
      break;

    case '*':
      sgn = 1;
      if (m1[0] == '-')
	{
	  lShift (m1, 1);
	  sgn = -sgn;
	}
      if (m2[0] == '-')
	{
	  lShift (m2, 1);
	  sgn = -sgn;
	}
      if (sgn == 1)
	sprintf (ans, "M|%s*%s", m1, m2);
      else
	sprintf (ans, "M|-%s*%s", m1, m2);
      break;

    case '.':
      if (m2[0] != '-')
	sprintf (ans, "M|F(%s,%s)", m1, m2);
      else if (m1[0] != '-')
	sprintf (ans, "M|F(%s,%s)", m2, m1);
      else
	{
	  lShift (m1, 1);
	  lShift (m2, 1);
	  sprintf (ans, "M|F(%s,%s)", m1, m2);
	}
      break;

    case '^':
      sprintf (ans, "M|%s^%s", m1, m2);
    }
  return (pointer) ans;
}


static void *
uact0 (char *ch, pointer mm)
{
  char *m, *ans;

  m = (char *) mm;
  ans = (char *) m_alloc (strlen (m) + 10);

/*   if (strcmp(ch,"-") == 0)*/
  if (m[0] == 'M')
    {
      if (m[2] == '-')
	sprintf (ans, "M|%s", m + 3);
      else
	sprintf (ans, "M|-%s", m + 2);
    }
  else
    sprintf (ans, "M|-(%s)", m + 2);

  if (strcmp (ch, "G") == 0)
    {
      if (r_reading0)
	sprintf (ans, "M|-g_(ln,%s)", m + 2);
      else
	sprintf (ans, "M|g_(ln,%s)", m + 2);
    }
  return (pointer) ans;
}

static void *
act_rcode (char *ch, int n, void **args)
{
  if (n == 1)
    return uact0 (ch, args[0]);
  if (n == 2)
    return bact0 (ch[0], args[0], args[1]);
  if (n == 4 && !strcmp (ch, "eps"))
    {
      int l =
	15 + strlen (args[0]) + strlen (args[1]) + strlen (args[2]) +
	strlen (args[3]);
      char *ans = (char *) m_alloc (l);
      sprintf (ans, "M|eps(%s,%s,%s,%s)", (char *) args[0] + 2, (char *) args[1] + 2,	/* eps(a,b,c,d) */
	       (char *) args[2] + 2, (char *) args[3] + 2);
      return ans;
    }
  return 0;
}


static void *
rd_rcode (char *s)
{
  char *p;
  /*char * GammaString; */
  int num;

  /*GammaString = (char *) m_alloc(12); */
  p = (char *) m_alloc (12);
  p[0] = 0;
  if (strlen (s) == 2 && s[1] > '0' && s[1] <= '9')
    {
      switch (s[0])
	{
	case 'p':
	case 'P':
	  num = s[1] - '0';
	  num = momsubst[num - 1];
	  if (num > 0)
	    sprintf (p, "M|p%d", num);
	  else
	    sprintf (p, "M|-p%d", -num);
	  break;

	case 'm':
	  num = s[1] - '0';
	  num = indsubst[num - 1];
	  sprintf (p, "M|m%d", num);
	  break;
	case 'M':
	  num = s[1] - '0';
	  num = indsubst[num - 1] - 1;
	  sprintf (p, "M|m%d", num);
	}
      if (strcmp (s, "G5") == 0)
	strcpy (p, "M|g_(ln,5_)");
    }
  if (!strlen (p))
    sprintf (p, "M|%s", s);
  return (void *) p;
}


static char *
pexpr (char p)
{
  static char snum[8];
  if (p >= 0)
    sprintf (snum, "p%d", p);
  else
    sprintf (snum, "(-p%d)", -p);
  return snum;
}


static char *
iexpr (int i)
{
  static char snum[6];

  sbld (snum, "m%d", i);
  return snum;
}


static void
var_def (void)
{
  int k, i = 1;
  writeF ("#define dp(n) DP[(n)]\n");
  for (k = 1; k <= nmodelvar; k++)
    {
      varlist modl = modelvars + k;
      if (!modl->func)
	{
	  writeF ("#define %s  va[%d]\n", modl->varname, i);
	  i++;
	}
    }

  for (k = 1; k <= nmodelvar; k++)
    {
      varlist modl = modelvars + k;
      if (modl->func)
	{
	  writeF ("#define %s  va[%d]\n", modl->varname, i);
	  i++;
	}
    }
}


static void
head (void)
{
  writeF ("*----------- VARIABLES ------------------ \n");
  writeF (" S  x,xx,xxx,y,mas,MZ,MW,Mt,Mb,Mc,Mu,Md,Ms,Mtop,Me,MH,");
  writeF (" Mm,SW,EE,GG,Cw,CW,");
  writeF (" Vus,Vcd,Vtb,Vcb,Vud,Vub,Vts,Vtd,Vcs,wtop,wZ,wW,wH,");
  writeF ("  s12,s23,s13,c12,c23,c13,Sqrt2,width,label;\n");
  writeF (" S  d0,...,d27;\n");
  writeF (" Auto S cc;\n");
  writeF (" V  p1,...,p16,ZERO_,z,mom,moment;\n");
  writeF (" I  m1,...,m16,s,t,ind,j,ln;\n");
  writeF (" CFunction dp,ANT,H;\n");
  writeF (" Function F,propDen;\n");
  writeF ("*\n");
  writeF ("Off statistics;\n\n");
  writeF ("#call Optimal()\n");
}



static void
writesubst (void)
{
  int i, j, l;
  int c;
  int ntot = getntot ();

  writeF ("*--- Mass declarations and Moment substitutions ---\n\n");
  writeF ("#procedure Substitution()\n");
  /* writeF("  .sort\n"); */
  writeF ("  repeat;\n");
  for (l = 1; l <= ntot; l++)
    writeF ("    id  p%d.p%d = %s*%s;\n", l, l, inoutmasses[l - 1],
	    inoutmasses[l - 1]);
  writeF ("\n");

  for (i = 1; i <= 3 * maxvert; i++)
    {
      l = momdep[i - 1][0];	/* For momdep this is a length V.E. */
      if (l > 1)
	{
	  writeF ("    id  p%d = ", i);
	  for (j = 1; j <= l; j++)
	    {
	      c = /* (char) */ momdep[i - 1][j];
	      if (c > 0)
		writeF ("+p%d", c);
	      else
		writeF ("-p%d", -c);
	    }
	  writeF (";\n");
	}
    }
  writeF ("  endrepeat;\n");
  writeF ("#endprocedure\n\n\n");
}


static int
MaxDotProd (int i)
{
  if (i == 2)
    return 0;
  if (i == 3)
    return 1;
  if (i == 4)
    return 4;
  if (i == 5)
    return 8;
  if (i == 6)
    return 13;
  if (i == 7)
    return 19;
  if (i == 8)
    return 26;
  return 0;
}


static void
optimal (void)
{
  int i, j, l, max, sign;
  int nin = getnin ();

  writeF ("#procedure Optimal()\n");
  writeF ("  .sort\n");
  for (i = 1; i <= 3 * maxvert; i++)
    {
      l = momdep[i - 1][0];
      if (l > 1)
	{
	  max = i;
	  i = 3 * maxvert + 1;
	}
    }

  writeF ("  L PmaxPmax=");
  sign = 1;
  if ((max - 2) > nin)
    sign *= -1;
  if ((max - 1) > nin)
    sign *= -1;
  if (sign < 0)
    writeF ("-1*(");
  else
    writeF ("(");

  for (i = 1; i < max; i++)
    for (j = i; j < max; j++)
      if ((i != max - 2) || (j != max - 1))
	{
	  sign = 1;
	  if (i > nin)
	    sign *= -1;
	  if (j > nin)
	    sign *= -1;
	  if (sign < 0)
	    writeF ("+");
	  else
	    writeF ("-");
	  writeF ("p%d.p%d", i, j);
	  if (i == j)
	    writeF ("/2");
	}
  writeF ("+p%d.p%d/2);\n", max, max);
  writeF ("  .sort\n");
  writeF ("  skip; nskip PmaxPmax;\n");
  writeF ("  #call Substitution()\n");
  writeF ("  .sort\n");
  /* writeF("  print PmaxPmax;\n"); */
  writeF ("  #$MaxDpMom= p%d.p%d;\n", max - 2, max - 1);
  writeF ("  #$MaxDpSub= PmaxPmax;\n");
  writeF ("  #$MaxDpInd= %d;\n", MaxDotProd (max - 1));
  writeF ("  .sort\n");
  writeF (" Drop PmaxPmax;\n");
  writeF ("  .sort\n");
  writeF ("#endprocedure\n\n\n");
}


static void
emitfactors (void)
{
  char num[STRSIZ];
  int c;
  int i, j, s, vln, pnum;

  writeF ("*---------- Factors ---------------\n");

/* ------- Symmetry Factor ------- */
  writeF (" L SymmFact=%d/%d;    * Diagram symmerty factor\n", vcs.symnum, vcs.symdenum);

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
	      }
	    c *= abs (prtclbase[pnum - 1].cdim);
	  }
    }
  writeF (" L AverFact=1/%d;  * Normalization factor of polarization average\n", c);

/* ----- Fermion factor  --------- */
  c = 0;

  for (i = 0; i < vcs.sizel; i++)
    for (j = 0; j < vcs.valence[i]; j++)
      {
	pnum = vcs.vertlist[i][j].partcl;
	if (prtclbase[pnum - 1].spin == 1
	    && (IN_PRTCL & vcs.vertlist[i][j].prop))
	  ++(c);
      }


  if ((c & 1) == 1)
    strcpy (num, "-1");
  else
    strcpy (num, "1");
  writeF (" L FermFact=%s;      * (-1)**(number of in-fermion particles)\n",
	  num);
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
   writeF(" Local VectFact=%s;      * (-1)**(number of vector propagators)\n",num);
*/
  /* ----- Color factor  --------- */
  writeF (" L ColorFact=%d/%d;    *  QCD color weight of diagram  \n",
	  vcs.clrnum, vcs.clrdenum);
  writeF ("*\n");
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

	if (ln->moment > 0 && !(ln->prop & (IN_PRTCL | OUT_PRTCL))
	    && pseudop (ln->partcl))
	  writeF ("L totFactor=totFactor/%s**2;\n",
		  prtclbase[ln->partcl - 1].massidnt);

      }

  writeF (" L denominator=");
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
  writeF (";\n");

}				/*   EmitDenoms   */

/*  Programm for constructing vertex-oriented representation for !
! amplitude diagram.                                             !
! A.Taranov 20.07.89                                            */

/*
static void  emitindex(indvertset index)
{  int m=1, n=0, l=1;

   if (setofb_eq0(index)) return;
    writeF(" Index ");                                              
   do
   {
      if (insetb(m,index))
      {
         if(n) writeF(",");
         n++;
         writeF("m%d",l);
         setofb_cpy(index,setofb_aun(index,setofb(m,_E)));
      }
      l++;
      if (!setofb_eq0(index)) m++;
   }  while (!(n == 1 || setofb_eq0(index)));
   writeF("$\n");
} 
*/


static void
reducemult (char *nameres, char *name1, char *name2, indvertset index)
{
  indvertset index1, index2;

  lvcpy (index1, index);
  lvcpy (index2, index);
  /* emitindex(index1); */
  writeF (" L %s=%s*%s;\n", nameres, name1, name2);
/*   while (!setofb_eq0(index1))                                           
   {
      emitindex(index1);                                       
      writeF(" %s:=%s$\n",nameres,nameres);      
   }
   if (!setofb_eq0(index2))
   {
      writeF(" RemInd  ");                                        
      l = 1; n = 1;
      while (!setofb_eq0(index2))
      {
         if (insetb(l,index2))
         {
            if (n != 1) writeF(",");                            
            ++(n);
	    writeF("m%d",l);                                   
         }
         setofb_cpy(index2,setofb_aun(index2,setofb(l,_E)));
         ++(l);
      }
      writeF("$\n");
   }
*/

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

  sbld (ans, "%s,(", iexpr (ind));
  if (gammaflag && fermvrt (v))
    {
      for (i = 1; i <= vcs.valence[v - 1]; i++)
	if (i != l)
	  sbld (ans, "%s+%s", ans, pexpr (vcs.vertlist[v - 1][i - 1].moment));
    }
  else
    sbld (ans, "%s%s", ans, pexpr (-vcs.vertlist[v - 1][l - 1].moment));
  sbld (ans, "%s)/%s", ans,
	prtclbase[vcs.vertlist[v - 1][l - 1].partcl - 1].massidnt);
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
  sbld (name1, "Vrt%d", arg1);
  sbld (name2, "Vrt%d", arg2);
  if (arg1 < arg2)
    strcpy (name0, name1);
  else
    strcpy (name0, name2);
  setofb_cpy (mind, setofb_its (index, setmassindex));
  if (setofb_eq0 (mind))
    {
      reducemult (name0, name1, name2, index);
      writeF (" .sort\n");
    }
  else
    {
      reducemult ("Vrt0", name1, name2, index);
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
	  writeF (" L VrtL=%s; L VrtR=%s;\n", name1, name2);
	  setofb_zero (mind);
	  m = mm;
	  n = 1;
	  while (m != 0)
	    {
	      if ((m & 1) != 0 && (l = pstn[n - 1]) != 0)
		{

		  /* writeF(" VrtL=(VrtL where %s);\n",subst(l,arg1));
		     writeF(" VrtR=(VrtR where %s);\n",subst(l,arg2)); */

		  writeF (" #call RepIND(L,%s)\n", subst (l, arg1));
		  writeF (" #call RepIND(R,%s)\n", subst (l, arg2));

		  setofb_cpy (mind, setofb_uni (mind, setofb (l, _E)));
		}
	      ++(n);
	      m >>= 1;
	    }
	  reducemult ("Vrt0", "Vrt0 + VrtL", "VrtR",
		      setofb_aun (index, mind));
	  writeF (" .sort\n");
	  /*  writeF(
	     " Clear Vrt_L $    Clear Vrt_R $\n");  */
	}
      writeF (" L %s=Vrt0;\n", name0);
      writeF (" .sort\n");
    }
  if (arg1 < arg2)
    {
      /* writeF(" Clear %s$\n",name2); */
      for (n = 1; n <= vcs.sizet; n++)
	if (vertmap[n - 1] == arg2)
	  vertmap[n - 1] = arg1;
    }
  if (arg2 < arg1)
    {
      /* writeF(" Clear %s$\n",name1);  */
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
  /*char * GammaString; */

  edgeinvert *ln = &vcs.vertlist[v - 1][l - 1];

  /*GammaString = (char *) m_alloc(18); */

  strcpy (mass, prtclbase[ln->partcl - 1].massidnt);
  if (prtclbase[ln->partcl - 1].hlp == '*')
    return mass;

  sprintf (ans, "g_(ln,%s)", pexpr (ln->moment));
  if (strcmp (mass, "0") != 0)
    sbld (ans, "(%s+%s)", ans, mass);
  else
    {
      if (fermionp (ln->partcl))
	{
	  if (prtclbase[ln->partcl - 1].hlp == 'L')
	    strcat (ans, "*(1-g_(ln,5_))/2");
	  /*{ sprintf(GammaString,"*(1-g_(%d,5_))/2",FermLoopIndicator);
	     strcat(ans,GammaString);
	     } */
	  else if (prtclbase[ln->partcl - 1].hlp == 'R')
	    strcat (ans, "*(1+g_(ln,5_))/2");
	  /*{ sprintf(GammaString,"*(1+g_(%d,5_))/2",FermLoopIndicator);
	     strcat(ans,GammaString);              
	     } */
	}
      else
	{
	  if (prtclbase[ln->partcl - 1].hlp == 'R')
	    strcat (ans, "*(1-g_(ln,5_))/2");
	  /*{ sprintf(GammaString,"*(1-g_(%d,5_))/2",FermLoopIndicator);
	     strcat(ans,GammaString);
	     } */
	  else if (prtclbase[ln->partcl - 1].hlp == 'L')
	    strcat (ans, "*(1+g_(ln,5_))/2");
	  /*{ sprintf(GammaString,"*(1+g_(%d,5_))/2",FermLoopIndicator);
	     strcat(ans,GammaString);
	     } */
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
  char *GammaString;

  GammaString = (char *) m_alloc (12);

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
	  sprintf (GammaString, "g_(ln,5_)");
	  fermloops[v - 1].g5 =
	    spos (GammaString,
		  parsedverts[fermloops[v - 1].vv[l - 1] - 1]) != 0;
	  ++(l);
	}
      if (fermloops[v - 1].g5)
	++(numberg5);
    }

  for (v = 1; v <= 2 * maxvert; v++)
    vertmap[v - 1] = 0;

  for (lpcount = 1; lpcount <= nloop; lpcount++)
    {
      writeF ("*------- Fermion loop calculation ------- \n");
      writeF (" \n");

      FermLoopIndicator++;	/* FermLoopIndicator */

      for (nn = 1; nn <= fermloops[lpcount - 1].len; nn++)
	{
	  if (vertexes[fermloops[lpcount - 1].vv[nn - 1] - 1].r_vert)
	    writeF ("* reversed vertex\n");

	  vertStr = parsedverts[fermloops[lpcount - 1].vv[nn - 1] - 1];
	  propStr = fermpropag (fermloops[lpcount - 1].vv[nn - 1],
				fermloops[lpcount - 1].ll[nn - 1]);
	  writeF (" L Vrt%d=%s*%s;\n", nn + lpcount - 1, vertStr, propStr);

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
      writeF (" .sort\n");

      for (v1 = 2; v1 <= fermloops[lpcount - 1].len; v1++)
	{
	  int newInd = fermloops[lpcount - 1].intln[v1 - 1];	/* vector line */
	  if (newInd != 0)
	    newInd =
	      vcs.vertlist[fermloops[lpcount - 1].vv[v1 - 1] - 1][newInd -
								  1].lorentz;

	  if (newInd == 0)
	    setofb_zero (indexs);
	  else
	    setofb_cpy (indexs, setofb (newInd, _E));
	  newInd = fermloops[lpcount - 1].intln2[v1 - 1];	/* second vector line */
	  if (newInd != 0)
	    newInd =
	      vcs.vertlist[fermloops[lpcount - 1].vv[v1 - 1] - 1][newInd -
								  1].lorentz;

	  if (newInd != 0)
	    setofb_cpy (indexs, setofb_uni (indexs, setofb (newInd, _E)));

	  r_mult (lpcount, lpcount + v1 - 1, indexs);

	}


      if (numberg5 == 1 && fermloops[lpcount - 1].g5)
	{
	  writeF (" #call G5ravno0(%d)\n", lpcount);

	  /* writeF(
	     " Vrt_%d:=Sub(A=0*A,Vrt_%d)$\n"
	     " Vrt_%d:=(Vrt_%d where A=>0*ZERO_)$\n" 
	     ,lpcount,lpcount); */
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
	  writeF (" Local Vrt%d=%s;\n", lpcount + 1, parsedverts[v1 - 1]);

	  writeF (" #call DefineF()\n");

	  r_mult (lpcount, lpcount + 1, indexs);
	}
      /* writeF(" trace4,1;\n");
         writeF(" contract;\n");                                   
         writeF(" .sort\n");
       */


      writeF (" L Fl%d=-Vrt%d;\n", lpcount, lpcount);
      writeF (" .sort\n");
      indicatorFl[lpcount] = 1;

      writeF (" Drop;");
      writeF (" Ndrop");
      for (k = 0; k < 20; k++)
	if (indicatorFl[k] == 1)
	  writeF (" Fl%d,", k);
      writeF (" numerator, GhostFact, totFactor;\n");
      writeF ("  .sort\n");
      writeF (" #call Substitution()\n");
      writeF (" #call RepG(%d,%d)\n", lpcount, FermLoopIndicator);
      /* writeF("  .sort\n"); */
      writeF ("*\n");

    }

  for (v1 = 1; v1 <= 2 * maxvert; v1++)
    fermmap[v1 - 1] = vertmap[v1 - 1];
  writeF ("*\n");
  gammaflag = FALSE;
}


static void
formblocks (int *vrtsize)
{
  int v;
  int l;
  int lori;
  int v1, v2, vv;
  int count;
  char vnum[11];

  /*    FormBlocks   */
  for (v = 1; v <= 2 * maxvert; v++)
    vertmap[v - 1] = fermmap[v - 1];	/*  copy fermMap      */
  for (count = 1; count <= nloop; count++)	/*  ferm loop vertex  */
    {
      sbld (vnum, "%d", count);
      writeF (" L Vrt%s=Fl%s;\n", vnum, vnum);
      indicatorFl[count] = 0;
    }
  count = nloop;
  for (v = 1; v <= vcs.sizet; v++)	/*  simple  vertex  */
    if (vertmap[v - 1] == 0)
      {
	++(count);
	vertmap[v - 1] = count;
	writeF (" L Vrt%d:=%s;\n", count, parsedverts[v - 1]);
      }

  writeF (" #call DefineF()\n");
  /*  end of vertex filling  */

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
			      setofb_uni (vertinfo[v1 - 1].ind,
					  setofb (lori, _E)));
		  np = vcs.vertlist[v - 1][l - 1].partcl;
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
	  finish ("*** Error! CompHEP is restarted");
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
  for (i = 0; i < vcs.sizet; ++i)
    free (parsedverts[i]);

  g5_exist = numberg5 != 0;
  for (i = n_vrt - 1; i >= 1; i--)
    {
      v1 = prgcode[i - 1][0];
      v2 = prgcode[i - 1][1];
      setofb_cpy (ind1, vertinfo[v1 - 1].ind);
      setofb_cpy (ind2, vertinfo[v2 - 1].ind);
      if (g5_exist && numberg5 == vertinfo[v1 - 1].g5 + vertinfo[v2 - 1].g5)
	{
	  writeF (" L Vrt1r=Vrt%d;\n", v1);
	  writeF (" L Vrt2r=Vrt%d;\n", v2);
	  writeF (" #call Iravno0(%d,%d)\n", v1, v2);

	  r_mult (v1, v2, setofb_its (ind1, ind2));
	  v0 = MIN (v1, v2);

	  writeF (" L VrtTmp=Vrt%d; \n", v0);
	  writeF (" .sort\n");
	  writeF (" L Vrt%d=Vrt1r; \n", v1);
	  writeF (" .sort\n");
	  writeF (" L Vrt%d=Vrt2r; \n", v2);
	  writeF (" .sort\n");
	  r_mult (v1, v2, setofb_its (ind1, ind2));
	  writeF (" L Vrt%d=Vrt%d+VrtTmp; \n", v0, v0);
	  writeF (" .sort\n");

	  g5_exist = FALSE;
	}
      else
	r_mult (v1, v2, setofb_its (ind1, ind2));
      vertinfo[MIN (v1, v2) - 1].g5 =
	vertinfo[v1 - 1].g5 + vertinfo[v2 - 1].g5;
      setofb_cpy (vertinfo[MIN (v1, v2) - 1].ind,
		  setofb_aun (setofb_uni (ind1, ind2),
			      setofb_its (ind1, ind2)));
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
	      writeF ("L Vrt1=Vrt1*(%s^2); \n", prtclbase[np - 1].massidnt);
	    else
	      writeF ("L Vrt1=Vrt1*(%s^2-p%d.p%d); \n",
		      prtclbase[np - 1].massidnt, mom, mom);
	  }
      }
  writeF ("*\n");
}


void
mk_formprograms (void)
{
  int k;
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
  int NumberDiagrams = 0;
  int *info_for_sqme;

  int ServiceSqme;

  goto_xy (1, 21);
  scrcolor (Yellow, Blue);
  print ("  FORM code generation \n");
  scrcolor (Red, BGmain);
  print (" Generated........\n");
  print (" current diagram :\n");
  scrcolor (Yellow, Blue);
  print (" Press Esc to halt FORM codes generation ");
  scrcolor (FGmain, BGmain);
  diagrq = fopen (DIAGRQ_NAME, "rb");
  ncalctot = 0;
  menuq = fopen (MENUQ_NAME, "rb");

/*-------------- service files --------------------------------*/
  outFileOpen (scat ("%sresults%cvar_def.h", pathtouser, f_slash));
  var_def ();
  outFileClose ();

/*-------------------------------------------------------------*/
  info_for_sqme = (int *) malloc ((subproc_sq + 1) * sizeof (int));

  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      rd_menu (menuq, 2, nsub, txt, &ndel, &ncalc, &nrest, &nrecord);
      fseek (diagrq, nrecord * sizeof (csdiagram), SEEK_SET);
      naxu = ndel + ncalc + nrest;

      info_for_sqme[nsub] = naxu;

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
	      outFileOpen (scat ("%sresults%cp%d_%d.frm", pathtouser, f_slash, nsub, ndiagr));
	      NumberDiagrams++;
	      for (k = 0; k < 20; k++)
		indicatorFl[k] = 0;
	      writeF ("#-\n");
	      writeF ("#: MaxTermSize 30000\n");
	      writeF ("#: ContinuationLines 30000\n");
	      writeLabel ('*');
	      writeF ("*\n");
	      transfdiagr (&csd, &vcs);
	      writeF ("\n\n#include procedur.prc\n\n");

	      cwtarg (&vcs);
	      if (vcs.clrnum == 0)
		{
		  writeF ("*-------  Zero color factor --------\n");
		  writeF ("L totFactor=0;\n");
		  writeF ("L numerator=0;\n");
		  writeF ("L denominator=1;\n");
		  writeF ("L NumDenominator=1;\n");
		  writeF (".sort\n");
		}
	      else
		{
		  generateghosts (&vcs, &gstlist);
		  if (gstlist == NULL)
		    {
		      writeF ("*-------  non-existent diagram  --------\n");
		      writeF ("L totFactor=0;\n");
		      writeF ("L numerator=0;\n");
		      writeF ("L denominator=1;\n");
		      writeF ("L NumDenominator=1;\n");
		      writeF (".sort\n");
		    }
		  else
		    {
		      goto_xy (40, 23);
		      print ("(%% %4d subdiagrams)", gstlist->maxnum);
		      writeF ("* The total number of diagrams %d\n", gstlist->maxnum);
		      preperdiagram ();
		      writesubst ();
		      optimal ();
		      head ();


		      emitfactors ();
		      diagramsrfactors (gstlist, &d_facts, &t_fact);
		      writeF (" L totFactor=%s;\n", rmonomtxt (*t_fact));
		      writeF (" .sort\n");
		      writeF (" L totFactor="
			      "totFactor*SymmFact*AverFact*FermFact*ColorFact;\n");
		      writeF (" .sort\n\n\n");
		      clrvm (t_fact->n.v);
		      clrvm (t_fact->d.v);
		      free (t_fact);

		      writeF ("L numerator=0;\n");

		      c = gstlist;
		      df = d_facts;
		      vcs_copy = vcs;
		      while (c != NULL)
			{
			  FermLoopIndicator = 0;	/* FermLoopIndicator */
			  coloringvcs (c);
			  writeF ("*  diagram  number =   %d\n", c->num);
			  DiagramToOutFile (&vcs, TRUE, '*');

			  {
			    int k;
			    int sgn = c->sgn;
			    for (k = 0; k < vcs.sizet; k++)
			      sgn *= vertexes[k].lgrnptr->factor;
			    writeF ("L  GhostFact=%d;\n", sgn);
			  }

			  findReversVert ();
			  attachvertexes ();

			  emitreducecode ();
			  writeF
			    (" L numerator=numerator +(%s)*GhostFact*Vrt1;\n",
			     smonomtxt (df->monom));
			  writeF (" .sort\n");
			  writeF (" Drop; Ndrop numerator, totFactor;\n");
			  writeF ("  .sort\n");
			  for (i = 1; i <= FermLoopIndicator; i++)
			    writeF (" trace4,%d;\n", i);
			  writeF (" contract;\n");
			  writeF (" .sort\n");
			  writeF (" #call Substitution()\n");
			  writeF ("  .sort\n");
			  writeF ("*\n");

			  vcs = vcs_copy;
			  c = c->next;
			  df = df->next;
			}

		      eraseslist (d_facts);
		      eraseghosts (gstlist);
		      vcs = vcs_copy;

		      writeF (" #call SubstOptimal()\n");
		      writeF (" #call OtherVar()\n");
		      writeF (" #call OutputTest(%d,%d)\n", nsub, ndiagr);
		      writeF (" #call DotProduct()\n");
		      writeF ("*\n");
		      writeF ("*\n");
		      writeF ("  .sort\n");
		      emitdenoms ();
		      writeF (" #call DefineDenominator()\n");
		    }
		}


	      writeF (" #call OutputC(%d)\n", NumberDiagrams);
	      writeF (".end;\n");
	      outFileClose ();
	      --nrest;
	      ++ncalctot;
	      if (escpressed ())
		goto exi;
	    }
	}
    }

exi:
  fclose (diagrq);
  fclose (menuq);
  clrbox (1, 21, 70, 24);
  ServiceSqme = service_sqme (NumberDiagrams, info_for_sqme);
}
