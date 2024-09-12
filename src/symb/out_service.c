/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <limits.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/parser.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "service2/include/getmem.h"
#include "chep_crt/include/chep_crt.h"
#include "polynom/include/polynom.h"
#include "polynom/include/lnum.h"

#include "physics.h"
#include "diaprins.h"
#include "procvar.h"
#include "pvars.h"
#include "pre_read.h"
#include "saveres.h"
#include "process.h"
#include "process_core.h"
#include "out_service.h"

#define xmax 76

FILE *outFile = NULL;
static int xpos;

static polyvars varsInfo =
{0, NULL};

int readSize = 0;
poly readBuff = NULL;

int outputLanguage = 0;

static void 
wrtln (void)
{
  fprintf (outFile, "\n");
  xpos = 1;
}


static void 
wrt_all (char *s)
{

  int l = strlen (s);

  if (!l)
    return;
  if (xpos > xmax)
    {
      fprintf (outFile, "\n ");
      xpos = 2;
      wrt_all (s);
      return;
    }

  if ((xpos - 1 + l) <= xmax)
    {
      fprintf (outFile, "%s", s);
      xpos += l;
      return;
    }

  l = xmax - xpos + 1;

/*   if(strchr("*=",s[l]) && strchr("*-/+",s[l-1])) l--;  */

  while (l && (!strchr ("*+-)(^=<> ", s[l - 1]) || strchr ("*+-)^=<>;", s[l])))
    l--;


  if (l == 0)
    {
      if (xpos > 2)
	{
	  fprintf (outFile, "\n ");
	  xpos = 2;
	  wrt_all (s);
	  return;
	}
      else
	{
	  l = xmax - xpos + 1;
	  if (strchr ("*=", s[l]) && strchr ("*-/+", s[l - 1]))
	    l++;
	  while (s[l] && strchr ("*+-)(^= ", s[l]) == NULL)
	    l++;
	}
    }

  fprintf (outFile, "%.*s\n ", l, s);
  xpos = 2;
  wrt_all (s + l);
}

static void 
wrt_fort (char *s)
{
  int l;

  l = strlen (s);
  if ((xpos + l) <= 73)
    {
      fprintf (outFile, "%s", s);
      xpos += l;
    }
  else
    {
      l = 73 - xpos;
      fprintf (outFile, "%.*s\n     >", l, s);
      xpos = 7;
      wrt_fort (s + l);
    }
}



void 
writeF (char *format,...)
{
  va_list args;
  char dump[STRSIZ], *beg, *nn;
  void (*wrt_0) (char *);

  va_start (args, format);
  vsprintf (dump, format, args);
  va_end (args);

  if (outputLanguage == 'f')
    wrt_0 = wrt_fort;
  else
    wrt_0 = wrt_all;

  beg = dump;
  while (TRUE)
    {
      nn = strchr (beg, '\n');
      if (nn == NULL)
	{
	  (*wrt_0) (beg);
	  return;
	}
      nn[0] = 0;
      (*wrt_0) (beg);
      wrtln ();
      beg = nn + 1;
    }
}

void 
outFileOpen (char *fname)
{
  outFile = fopen (fname, "w");
  xpos = 1;
}

void 
outFileClose (void)
{
  fclose (outFile);
}



void 
DiagramToOutFile (vcsect * vcs, int label, char comment)
{
  int numtot = getntot ();
/* AK: This is a point to block drawing pg-diagram */
  if (numtot < 7)		/* This is magick V.E. number */
    {
      if (xpos != 1)
	wrtln ();
      writeTextDiagram (vcs, label, comment, outFile);
    }
}

void 
readvardef (FILE * fres)
{
  monom template;
  char *end;

  vardef = &varsInfo;
  fread (&vardef->nvar, sizeof (vardef->nvar), 1, fres);

  if (vardef->nvar)
    {
      vardef->vars = m_alloc (vardef->nvar * sizeof (varinfo));
      fread (vardef->vars, vardef->nvar * sizeof (varinfo), 1, fres);
      end = (char *) &template.tail.power[vardef->vars[vardef->nvar - 1].wordpos];
    }
  else
    {
      end = (char *) &template.tail.power[0];
      vardef->vars = NULL;
    }
  readSize = end - (char *) &template.coef.num;
  if (readBuff)
    free (readBuff);
  readBuff = m_alloc (end - (char *) &template);
}


void 
clearvardef (void)
{
  clearVars (vardef);
  if (readBuff)
    {
      free (readBuff);
      readBuff = NULL;
    }
}


static int 
readmonom (FILE * fres, char *txt)
{
  NUM_TYPE l;

  int deg, n;

  fread (&readBuff->coef.num, readSize, 1, fres);
  if (!txt)
    return readBuff->coef.num;

  l = readBuff->coef.num;

  if (!l)
    return 0;

  if (l == 1 || l == -1)
    txt[0] = 0;
  else
    {
#ifdef NUM_DOUBLE
      if (ABS (l) > 1.E15)
	sprintf (txt, "%+" NUM_STR ".", l);
      else
#endif
	sprintf (txt, "%+" NUM_STR, l);
    }

  for (n = 0; n < vardef->nvar; n++)
    {
      deg = (readBuff->tail.power[vardef->vars[n].wordpos - 1] /
	     vardef->vars[n].zerodeg) %
	vardef->vars[n].maxdeg;
      if (deg)
	{
	  sprintf (txt + strlen (txt), "*%s", vararr[vardef->vars[n].num].alias);
	  if (deg > 1)
	    sprintf (txt + strlen (txt), "^%d", deg);
	}
    }

  if (!txt[0])
    {
#ifdef NUM_DOUBLE
      if (ABS (l) > 1.E15)
	sprintf (txt, "%+" NUM_STR ".", l);
      else
#endif
	sprintf (txt, "%+" NUM_STR, l);
    }
  else
    {
      if (l == 1)
	txt[0] = '+';
      else if (l == -1)
	txt[0] = '-';
    }
  return 1;
}


void 
momentToString (char *momStr, char *outstr)
{
  int m = 0;

  strcpy (outstr, "");
  while (momStr[m])
    {
      if (momStr[m] <= getnin ())
	strcat (outstr, "-");
      else if (m)
	strcat (outstr, "+");
      sprintf (outstr + strlen (outstr), "p%d", momStr[m++]);
    }
}


void 
rewritepolynom (FILE * fsrc)
{
  char monomtxt[STRSIZ];
  if (!readmonom (fsrc, monomtxt))
    {
      wrt_all ("0");
      return;
    }
  if (monomtxt[0] == '+')
    wrt_all (monomtxt + 1);
  else
    wrt_all (monomtxt);
  while (readmonom (fsrc, monomtxt))
    wrt_all (monomtxt);
}


void 
findPrtclNum (char *procName, int *prtclNum)
{
  int k;
  char *pname, txt1[STRSIZ];
  int pnum;

  strcpy (txt1, procName);
  memcpy (strstr (txt1, "->"), ", ", 2);
  pname = strtok (txt1, ",");
  k = 1;
  while (pname != NULL)
    {
      trim (pname);
      pnum = locateinbase (pname);
      prtclNum[k] = pnum;
      k++;
      pname = strtok (NULL, ",");
    }
}

void 
emitconvlow (int *prtclNum)
{
  int i, j, ntot;
  char pmass[MAXINOUT + 1][11];	/* for locateinbase */
  char s[MAXINOUT + 1];
  int numin = getnin ();
  int numout = getnout ();

  ntot = getntot ();
  for (i = 1; i <= ntot; i++)
    strcpy (pmass[i], prtclbase[prtclNum[i] - 1].massidnt);

  for (i = 1; i <= numin; i++)
    s[i] = 1;
  for (i = numin + 1; i <= ntot; i++)
    s[i] = -1;

  switch (outputLanguage)
    {
    case 'R':
      writeF ("\nlet p%d = ", ntot);
      break;
    case 'F':
      writeF ("\nid p%d = ", ntot);
      break;
    case 'M':
      writeF ("\np%d = ", ntot);
      break;
    }
  for (i = 1; i <= numin; i++)
    writeF ("+p%d", i);
  for (i = 1; i <= numout - 1; i++)
    writeF ("-p%d", i + numin);
  writeF (";\n");

  for (i = 1; i <= ntot - 1; i++)
    switch (outputLanguage)
      {
      case 'R':
	writeF ("mass p%d  = %s; Mshell p%d;\n", i, pmass[i], i);
	break;
      case 'F':
	writeF ("id p%d.p%d  = %s^2;\n", i, i, pmass[i]);
	break;
      case 'M':
	writeF ("p%d/: SC[p%d,p%d] =%s^2;\n", i, i, i, pmass[i]);
	break;
      }

  switch (outputLanguage)
    {
    case 'R':
      writeF ("let p%d.p%d = ", ntot - 2, ntot - 1);
      break;
    case 'F':
      writeF ("id p%d.p%d = ", ntot - 2, ntot - 1);
      break;
    case 'M':
      writeF ("p%d/: SC[p%d,p%d] = ", ntot - 2, ntot - 2, ntot - 1);
      break;
    }

  writeF ("%d*(%s^2", s[ntot - 1] * s[ntot - 2], pmass[ntot]);
  for (i = 1; i <= ntot - 1; i++)
    writeF ("-%s^2", pmass[i]);
  for (i = 2; i <= ntot - 1; i++)
    for (j = 1; j <= i - 1; j++)
      if (j < (ntot - 2))
	switch (outputLanguage)
	  {
	  case 'R':
	  case 'F':
	    writeF ("%+d*p%d.p%d", -2 * s[j] * s[i], j, i);
	    break;
	  case 'M':
	    writeF ("%+d*SC[p%d,p%d]", -2 * s[j] * s[i], j, i);
	    break;
	  }
  writeF (")/2;\n");
}

void 
writeLabel (char comment)
{
  if (xpos != 1)
    wrtln ();
  fprintf (outFile, "%c    ==============================\n", comment);
  fprintf (outFile, "%c    *  %s *\n", comment, getname ());
  fprintf (outFile, "%c    ==============================\n", comment);
}



void 
makeOutput (void (*startOutput) (FILE *, int, int *, int),
	    void (*diagramOutput) (FILE *, vcsect *, catrec *),
	    void (*endOutput) (int *)
)
{
  catrec cr;
  int ndel, ncalc, nrest;
  long recpos;
  long count;
  int graphOn;
  shortstr txt;
  int prtclNum[MAXINOUT + 1];
  csdiagram csdiagr;
  vcsect vcs;
  FILE * archiv = fopen (ARCHIV_NAME, "rb");

  informline (0, 1);
  catalog = fopen (CATALOG_NAME, "rb");
  menuq = fopen (MENUQ_NAME, "rb");
  count = 0;
  graphOn = 1;

  if (graphOn)
    diagrq = fopen (DIAGRQ_NAME, "rb");

  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      rd_menu (menuq, 2, nsub, txt, &ndel, &ncalc, &nrest, &recpos);

      findPrtclNum (txt, prtclNum);
      if (ncalc != 0)
	{
	  startOutput (archiv, nsub, prtclNum, ncalc);
	  fseek (catalog, 0, SEEK_SET);
	  while (FREAD1 (cr, catalog))
	    {
	      if (cr.nsub_ == nsub)
		{
		  if (graphOn)
		    {
		      fseek (diagrq, (cr.ndiagr_ + recpos - 1) * sizeof (csdiagram), SEEK_SET);
		      FREAD1 (csdiagr, diagrq);
		      transfdiagr (&csdiagr, &vcs);
		      diagramOutput (archiv, &vcs, &cr);
		    }
		  else
		    diagramOutput (archiv, NULL, &cr);
		  ++(count);
		  if (informline (count, subproc_sq))
		    goto escexit;
		}
	    }
	escexit:
	  endOutput (prtclNum);
	}
    }
  fclose (catalog);
  fclose (archiv);
  fclose (menuq);
  if (graphOn)
    fclose (diagrq);
  informline (subproc_sq, subproc_sq);
}
