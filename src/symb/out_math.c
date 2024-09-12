/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"

#include "physics.h"
#include "procvar.h"
#include "out_service.h"
#include "saveres.h"
#include "process.h"
#include "process_core.h"
#include "out_math.h"


static void 
writeprocessname (int *prtclNum)
{
  int i;
  int numin = getnin ();
  int numtot = getntot ();

  writeF ("  process  ");
  for (i = 1; i <= numin; i++)
    {
      writeF ("%s(p%d)", prtclbase[prtclNum[i] - 1].name, i);
      if (i < numin)
	writeF ("+");
      else
	writeF ("->");
    }

  for (i = numin + 1; i <= numtot; i++)
    {
      writeF ("%s(p%d)", prtclbase[prtclNum[i] - 1].name, i);
      if (i < numtot)
	writeF ("+");
      else
	writeF ("\n");
    }
}


static void 
emitprocessname (int *prtclNum)
{
  int i;
  int numin = getnin ();
  int numtot = getntot ();

  writeF ("inParticles = {");
  for (i = 1; i <= numin; i++)
    {
      writeF ("\"%s\"", prtclbase[prtclNum[i] - 1].name);
      if (i < numin)
	writeF (",");
      else
	writeF ("}\n");
    }

  writeF ("outParticles = {");

  for (i = numin + 1; i <= numtot; i++)
    {
      writeF ("\"%s\"", prtclbase[prtclNum[i] - 1].name);
      if (i < numtot)
	writeF (",");
      else
	writeF ("}\n");
    }

}



static void 
emitexpression (FILE * fsrc, catrec * cr)
{
  int i;

  fseek (fsrc, cr->factpos, SEEK_SET);
  readvardef (fsrc);
  writeF ("totFactor = ((");
  rewritepolynom (fsrc);
  writeF (")/(");
  rewritepolynom (fsrc);
  writeF ("));");
  writeF ("\n");
  clearvardef ();

  fseek (fsrc, cr->rnumpos, SEEK_SET);
  readvardef (fsrc);
  writeF ("numerator =(");
  rewritepolynom (fsrc);
  writeF (");\n");
  clearvardef ();

  fseek (fsrc, cr->denompos, SEEK_SET);
  readDenominators (fsrc);

  writeF ("denominator =");
  for (i = 0; i < denrno; i++)
    {
      char momStr[20];
      if (i)
	writeF ("*");
      else
	writeF ("(");
      momentToString (denom[i].momStr, momStr);
      writeF ("propDen[%s,%s,%s]", momStr, vararr[denom[i].mass].alias,
	      vararr[denom[i].width].alias);
      if (denom[i].power > 1)
	writeF ("^%d", denom[i].power);
    }
  if (i)
    writeF (");\n");
  else
    writeF ("1;\n");
}

static void 
modifyFunc (char *lch)
{
  int c;
  lch[0] = toupper (lch[0]);
  while (lch[0] != '(')
    lch++;
  lch[0] = '[';
  for (c = 1, lch++; c; lch++)
    if (lch[0] == ')')
      c--;
    else if (lch[0] == '(')
      c++;
  lch[-1] = ']';
}


static void 
writeparameters (int nsub)
{
  int k = 0;
  int first = 1;
  char s[STRSIZ];
  char *lch;

  writeF ("\n");
  writeF ("parameters={\n");

  for (k = 1; k <= nmodelvar; k++)
    {
      if (vararr[k].used && !modelvars[k].func)
	{
	  if (first)
	    {
	      first = 0;
	      writeF (" ");
	    }
	  else
	    writeF (",");
	  writeF ("%s -> ", vararr[k].alias);
	  sprintf (s, "%17.11E", vararr[k].tmpvalue);
	  lch = strchr (s, 'E');
	  if (lch)
	    {
	      int d;
	      sscanf (lch + 1, "%d", &d);
	      sprintf (lch, "*10^(%d)", d);
	    }
	  writeF ("%s\n", s);
	}
    }

  writeF ("           };\n");
  writeF ("\n");
  first = 1;
  writeF ("substitutions={\n");

  for (k = nmodelvar; k; k--)
    {
      if (vararr[k].used && modelvars[k].func)
	{
	  sscanf (modelvars[k].func, "%[^|]", s);
	  trim (s);

	  while ((lch = strstr (s, "sqrt(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "sin(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "cos(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "tan(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "asin(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "acos(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "atan(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "exp(")))
	    modifyFunc (lch);
	  while ((lch = strstr (s, "log(")))
	    modifyFunc (lch);

	  if (first)
	    {
	      first = 0;
	      writeF (" ");
	    }
	  else
	    writeF (",");
	  writeF ("%s->%s", vararr[k].alias, s);
	  writeF ("\n");
	}
    }
  writeF ("              };\n");
}


static void 
startMath (FILE * frsc, int nsub, int *prtclNum, int ncalc)
{
  char f_name[STRSIZ];
  outputLanguage = 'M';
  initvararray (frsc, nsub, outputLanguage);

  sprintf (f_name, "%sresults%csymb%d.m", pathtouser, f_slash, nsub);
  outFileOpen (f_name);

  writeF ("(*\n");
  writeLabel (' ');
  writeprocessname (prtclNum);
  writeF ("*)\n");
  writeparameters (nsub);
  writeF ("\n");
  emitprocessname (prtclNum);
  writeF ("\n");
  writeF ("SetAttributes[ SC, Orderless ];\n");
  writeF ("\n");
  writeF ("SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];\n");
  writeF ("\n");
  writeF ("SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]\n");
  writeF ("\n");
  writeF ("\n");

  emitconvlow (prtclNum);

  writeF ("\ninitSum;\n");
}


static void 
diagramMath (FILE * fsrc, vcsect * vcs, catrec * cr)
{
  writeF ("\n(*\n");
  writeF ("  Diagram  %d in subprocess %d\n", cr->ndiagr_, cr->nsub_);
  if (vcs != NULL)
    DiagramToOutFile (vcs, 0, ' ');
  writeF ("*)\n");
  emitexpression (fsrc, cr);

  writeF ("\naddToSum;\n");
}

static void 
endMath (int *prtclNum)
{
  writeF ("\nfinishSum;\n");
  outFileClose ();
}

void 
makeMathOutput (void)
{
  makeOutput (startMath, diagramMath, endMath);
}
