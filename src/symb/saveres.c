/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ---------------------------------------------------
*/
#include <limits.h>
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"
#include "chep_crt/include/crt.h"
#include "polynom/include/polynom.h"

#include "physics.h"
#include "pvars.h"
#include "sos.h"
#include "ghosts.h"
#include "rfactor.h"
#include "screen.h"
#include "saveres.h"

denom_struct denom[2 * maxvert - 2] = {{0, 0., 0., {""}}};

byte denrno = 0;

char denStr[2 * maxvert - 2][MAXINOUT] = {{0}};

static void 
wAbort (void)
{
  saveent (menulevel);
  warnanykey (5, 20, "Error in writing on the disk. \n"
	      "Check the existence of the \n"
	      "'tmp' and 'results' directories \n"
	      "or the existence of free disk space");
  finish ("Error...");
  exit (0);			/*  End of work  */
}



static void 
savevars (FILE * file, polyvars * v)
{
  FWRITE1 (v->nvar, file);
  if (v->nvar &&
   fwrite (v->vars, (v->nvar) * sizeof (varinfo), 1, file) != 1)
    wAbort ();
}

static void 
savepoly (FILE * file, poly p, polyvars * v)
{

  int width;
  char *b, *e;
  poly zero = plusone ();

  zero->coef.num = 0;

  b = (char *) &(zero->coef.num);
  if (v->nvar)
    e = (char *) &(zero->tail.power[v->vars[v->nvar - 1].wordpos]);
  else
    e = (char *) &(zero->tail.power[0]);
  width = e - b;

  while (p)
    {
      if (fwrite (&(p->coef.num), width, 1, file) != 1)
	wAbort ();
      p = p->next;
    }

  if (fwrite (&(zero->coef.num), width, 1, file) != 1)
    wAbort ();
}


static void 
printpoly (FILE * file, poly p, polyvars * v)
{
  int deg, n=0, ImConjKey = 0, all_deg = 0;
 
   while (p != NULL)
    {
     for (n = 0; n < v->nvar; n++) 
      { 
        deg = (p->tail.power[v->vars[n].wordpos - 1] /
	     v->vars[n].zerodeg) %
	v->vars[n].maxdeg;
        
        if(deg==2||deg==6||deg==10||deg==14) all_deg = 1; 
        else all_deg = 0;
        if( strncmp(v->vars[n].name,"ImConj",6)==0 && all_deg==1 ) p->coef.num = -(p->coef.num);  
/*        printf ("\n%d*%s^%d\n", p->coef.num, v->vars[n].name, deg); */
        }
/*        if(ImConjKey==1) printf ("\n ImConjKey=%d\n", ImConjKey);*/
/*        if(ImConjKey==1) p->coef.num = -(p->coef.num);  */
        p = p->next;                  
    }

}

void 
save_analitic_results (FILE * fres, poly rnum, poly factn, poly factd, polyvars * vars1, polyvars * vars2, vcsect vcs, 
int status_key)
{
  catrec cr;
  int i;
  int m;

  diskerror = wAbort;

  cr.nsub_ = nsub;
  cr.ndiagr_ = ndiagr;

  cr.factpos = ftell (fres);
  savevars (fres, vars2);
  savepoly (fres, factn, vars2);
  savepoly (fres, factd, vars2);

  cr.rnumpos = ftell (fres);
  printpoly (fres, rnum, vars1);
  savevars (fres, vars1);
  savepoly (fres, rnum, vars1);    

  cr.denompos = ftell (fres);

  if(status_key == 3) cr.status = 4;
  else  cr.status = 1;

  calcdenominators (vcs);
  FWRITE1 (denrno, fres);        /*  number of demominators  */
  for (i = 0; i < denrno; i++)
    {
      FWRITE1 (denom[i].power, fres);      /*  power  1 or 2  */
      FWRITE1 (denom[i].mass, fres);
      FWRITE1 (denom[i].width, fres);
      m = 0;
      do 
        {
          FWRITE1 (denom[i].momStr[m], fres);
        }
      while (denom[i].momStr[m++]);
    }

  FWRITE1 (cr, catalog);
  diskerror = NULL;
}

void 
save_empty_analitic_results (FILE * fres, vcsect vcs)
{
  catrec cr;
  int i;
  int m;

  diskerror = wAbort;

  cr.nsub_ = nsub;
  cr.ndiagr_ = ndiagr;
  cr.factpos = ftell (fres);
  cr.rnumpos = ftell (fres);
  cr.denompos = ftell (fres);
  cr.status = 0;

  calcdenominators (vcs);
  FWRITE1 (denrno, fres);        /*  number of demominators  */
  for (i = 0; i < denrno; i++)
    {
      FWRITE1 (denom[i].power, fres);      /*  power  1 or 2  */
      FWRITE1 (denom[i].mass, fres);
      FWRITE1 (denom[i].width, fres);
      m = 0;
      do 
        {
          FWRITE1 (denom[i].momStr[m], fres);
        }
      while (denom[i].momStr[m++]);
    }

  FWRITE1 (cr, catalog);
  diskerror = NULL;
}
