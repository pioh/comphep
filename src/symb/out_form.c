/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/unix_utils.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"

#include "physics.h"
#include "procvar.h"
#include "out_service.h"
#include "saveres.h"
#include "out_form.h"

static void
writeparameters (int nsub)
{
  int k;
  int first;
  char ch;

  writeF ("\n");
  first = TRUE;
  ch = ' ';
  writeF ("Symbols\n");

  for (k = 1; k <= nmodelvar; k++)
    {
      if (vararr[k].used && !modelvars[k].func)
        {
          writeF (" %c%s\n", ch, vararr[k].alias);
          ch = ',';
        }
    }
  writeF (" ;\n");
}


static void
writefunctions (int nsub)
{
  int k;
  char s[STRSIZ];

  writeF ("\n");
  for (k = 1; k <= nmodelvar; k++)
    {
      if (vararr[k].used && modelvars[k].func)
        {
          sscanf (modelvars[k].func, "%[^|]", s);
          trim (s);
          writeF ("id %s=%s;\n", vararr[k].alias, s);
        }
    }
}


static void
startForm (FILE * fsrc, int nsub, int *prtclNum, int ncalc)
{
/*
  char fname[STRSIZ];

  outputLanguage = 'F';
  initvararray (fsrc, nsub, outputLanguage);
  sprintf (fname, "%sresults%cusr%d.frm", pathtouser, f_slash, nsub);
  outFileOpen (fname);

  writeF ("\n#-\nCFunction Sqrt;\n#procedure userWork(nnn)\n");

  writefunctions (nsub);

  emitconvlow (prtclNum);

  writeF ("#if 'nnn' == %d\n.end\n#endif\n#endprocedure\n", ncalc);
  outFileClose ();

  sprintf (fname, "%sresults%csum%d.frm", pathtouser, f_slash, nsub);
  outFileOpen (fname);
  writeLabel ('*');
  writeF ("vector p1,p2,p3,p4,p5,p6;\n");
  writeparameters (nsub);
  writeF ("CFunction den;\n");
  writeF ("Symbol factor;\n");
  writeF ("#include usr%d.frm\n", nsub);
*/
  char fname[STRSIZ];
  char *mstr;

  outputLanguage = 'F';
  initvararray (fsrc, nsub, outputLanguage);
  sprintf (fname, "%sresults%csum%d.frm", pathtouser, f_slash, nsub);
  outFileOpen (fname);
  writeF ("#-\n");
  writeF ("#procedure Output(x,xx)\n");
  writeF (".sort\n");
  writeF ("#$delta=FromC;\n");
  writeF ("#if('$delta'!=0)\n");
  mstr = " %e\",FromC ";
  writeF ("  #write <./results/testresults.txt> \"delta'x'_'xx'= %s\n", mstr);
  writeF ("#endif\n");
  writeF (".sort\n");
  writeF ("#endprocedure\n\n");
  writeF ("*----------- VARIABLES ------------------ \n");
  writeF
    (" Symbols  x,xx,MZ,MW,Mt,Mb,Mc,Mu,Md,Ms,Mtop,Me,MH,Mm,SW,EE,GG,Cw,CW,");
  writeF (" Vus,Vcd,Vtb,Vcb,Vud,Vub,Vts,Vtd,Vcs,wtop,wZ,wW,wH,");
  writeF (" s12,s23,s13,c12,c23,c13,Sqrt2;\n");
  writeF (" Vectors  p1,...,p16;\n");
  writeF (" Function FromForm;\n");
  writeF ("*\n");
  writeF ("Off statistics;\n\n");

  writeF (" #write <./results/testresults.txt> \"file:  delta:\"\n\n");
  /* writeparameters(nsub);       */
}


static void
diagramForm (FILE * fsrc, vcsect * vcs, catrec * cr)
{
  int i;

  writeF ("\n*Diagrama number %d-%d;\n", cr->nsub_, cr->ndiagr_);
  writeF ("#include ./results/t%d_%d.prc\n", cr->nsub_, cr->ndiagr_);

  if (vcs != NULL)
    {
    }

  fseek (fsrc, cr->factpos, SEEK_SET);
  readvardef (fsrc);
  writeF ("Local FACTOR = (");
  rewritepolynom (fsrc);
  writeF (")/(");
  rewritepolynom (fsrc);
  writeF (");\n");
  clearvardef ();
  writeF (".sort\n");
  writeF ("Local FromC = FACTOR*(", cr->ndiagr_);

  fseek (fsrc, cr->rnumpos, SEEK_SET);
  readvardef (fsrc);
  rewritepolynom (fsrc);
  writeF (") -FromForm();");
  clearvardef ();

  fseek (fsrc, cr->denompos, SEEK_SET);
  readDenominators (fsrc);

  for (i = 0; i < denrno; i++)
    {
      char momStr[20];

      momentToString (denom[i].momStr, momStr);
      /* writeF ("\n *den(%s,%s,%d)", momStr,
         vararr[denom[i].mass].alias, -denom[i].power); */
    }
/*  writeF (";\n"); */

  writeF ("\n#call FromForm%dp%d()\n", cr->nsub_, cr->ndiagr_);
  writeF ("#call Output(%d,%d)\n", cr->nsub_, cr->ndiagr_);
  writeF (".sort\n");
  writeF ("Drop;\n");
  writeF (".sort\n");
}

static void
endForm (int *prtclNum)
{
  writeF (".end\n");
  outFileClose ();
}


void
makeFormOutput (void)
{
  makeOutput (startForm, diagramForm, endForm);
}
