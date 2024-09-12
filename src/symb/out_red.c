/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
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
#include "out_red.h"


static void 
writeprocessname (int *prtclNum)
{
  int i;
  int numin = getnin ();
  int numtot = getntot ();

  writeF ("inParticles:={");
  for (i = 1; i <= numin; i++)
    if (i == numin)
      writeF ("\"%s\"}$\n", prtclbase[prtclNum[i] - 1].name);
    else
      writeF ("\"%s\",", prtclbase[prtclNum[i] - 1].name);

  writeF ("outParticles:={");
  for (i = numin + 1; i <= numtot; i++)
    if (i == numtot)
      writeF ("\"%s\"}$\n", prtclbase[prtclNum[i] - 1].name);
    else
      writeF ("\"%s\",", prtclbase[prtclNum[i] - 1].name);
}

static void 
writeparameters (int nsub)
{

  int k;
  char ch;
  char s[STRSIZ];

  writeF ("%%\n");
  ch = ' ';
  writeF ("parameters:={");

  for (k = 1; k <= nmodelvar; k++)
    {
      if (vararr[k].used && !modelvars[k].func)
	{
	  writeF ("%c%s=>%E", ch, vararr[k].alias, vararr[k].tmpvalue);
	  ch = ',';
	}
    }
  writeF ("}$\n");

  writeF ("%%\n");
  ch = ' ';
  writeF ("substitutions:={");
  for (k = nmodelvar; k; k--)
    {
      if (vararr[k].used && modelvars[k].func)
	{
	  sscanf (modelvars[k].func, "%[^|]", s);
	  trim (s);
	  writeF ("%c%s=>%s", ch, vararr[k].alias, s);
	  ch = ',';
	}
    }
  writeF ("}$\n");
}

static void 
emitexpression (FILE * fsrc, catrec * cr)
{
  int i;

  fseek (fsrc, cr->factpos, SEEK_SET);
  readvardef (fsrc);
  writeF ("totFactor:=(");
  rewritepolynom (fsrc);
  writeF (")/(");
  rewritepolynom (fsrc);
  writeF (")$");
  writeF ("\n");
  clearvardef ();

  fseek (fsrc, cr->rnumpos, SEEK_SET);
  readvardef (fsrc);
  writeF ("numerator:=");
  rewritepolynom (fsrc);
  writeF ("$\n");
  clearvardef ();

  fseek (fsrc, cr->denompos, SEEK_SET);
  readDenominators (fsrc);

  writeF ("denominator:=");
  if (denrno)
    for (i = 0; i < denrno; i++)
      {
	char momStr[20];
	momentToString (denom[i].momStr, momStr);

	if (i)
	  writeF ("*");
	writeF ("propDen(%s,%s,%s)", momStr, vararr[denom[i].mass].alias,
		vararr[denom[i].width].alias);
	if (denom[i].power != 1)
	  writeF ("^%d", denom[i].power);
      }
  else
    writeF ("1");
  writeF ("$\n");
}


static void 
startReduce (FILE * fsrc, int nsub, int *prtclNum, int ncalc)
{
  char f_name[STRSIZ];

  outputLanguage = 'R';
  initvararray (fsrc, nsub, outputLanguage);
  sprintf (f_name, "%sresults%csymb%d.red", pathtouser, f_slash, nsub);
  outFileOpen (f_name);
  writeLabel ('%');
  writeprocessname (prtclNum);
  writeparameters (nsub);
  writeF ("\n\n");
  writeF ("vector p1,p2,p3,p4,p5,p6$\n");
  emitconvlow (prtclNum);

  writeF ("\nvector !=p_,!=q_$\n");
  writeF ("operator propDen$\n");
  writeF ("for all p_,q_,m,w let propDen(0*p_+q_,m,w)=propDen(q_,m,w)$\n");
  writeF ("for all p_,m,w such that ordp(p_,-p_) "
	  "let propDen(p_,m,w)=propDen(-p_,m,w);$\n\n");
  writeF ("initSum();\n");
}

static void 
diagramReduce (FILE * fsrc, vcsect * vcs, catrec * cr)
{
  writeF ("\n");
  writeF ("DiagrNumber:=\"%d_%d\"$\n", cr->nsub_, cr->ndiagr_);
  writeF ("\n");
  if (vcs)
    DiagramToOutFile (vcs, 0, '%');
  emitexpression (fsrc, cr);
  writeF ("\n\n");
  writeF ("addToSum()$\n");
}

static void 
endReduce (int *prtclNum)
{
  writeF ("finishSum();\n");
  writeF ("End$\n");
  outFileClose ();
}

void 
makeReduceOutput (void)
{
  makeOutput (startReduce, diagramReduce, endReduce);
}
