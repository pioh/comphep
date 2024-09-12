/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "service2/include/parser.h"
#include "service2/include/getmem.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"

#include "physics.h"
#include "pvars.h"
#include "out_service.h"
#include "saveres.h"
#include "process.h"
#include "process_core.h"
#include "procvar.h"

int nProcessVar = 0;

singlevardescription *vararr = NULL;


static void *PP = (void *) "PP";

static void *
rd_hiddenVars (char *s)
{
  int l;
  if (isdigit (*s))
    return PP;
  for (l = 0; l <= nmodelvar; l++)
    {
      if (!strcmp (modelvars[l].varname, s))
	{
	  vararr[l].used = 1;
	  break;
	}
    }
  return PP;
}

static void *
act_hiddenVars (char *ch, int n, void **mm1)
{
  return PP;
}

void 
seekArchiv1 (FILE * f, long n)
{
  fseek (f, n, SEEK_SET);
}


static void 
koeff_check (void)
{
  int l;
 
  for (l = 0; l <= nmodelvar; l++)
    {
      if (strcmp (modelvars[l].varname, "koeff") == 0)
	{
	  vararr[l].used = 1;
	  break;    
	}
    }
}



int 
initvararray (FILE * fres, int nsub, char key)
{
  int i, j, k, kk, l;
  catrec cr;
  int numtot = getntot ();

  polyvars allVars =
  {0, NULL};
  int nvar, nfunc;

  nProcessVar = nmodelvar + 2 + ((MAXINOUT * (MAXINOUT - 1) / 2));
  if (vararr)
    free (vararr);
  vararr = (singlevardescription *) m_alloc (nProcessVar * sizeof (singlevardescription));

  /*sprintf(vararr[0].alias,"0",k); */
  sprintf (vararr[0].alias, "0");
  vararr[0].tmpvalue = 0;
  vararr[0].used = FALSE;

  for (k = 1; k < nProcessVar; k++)
    {
      sprintf (vararr[k].alias, "#%d", k);
      vararr[k].tmpvalue = 0;
      vararr[k].used = FALSE;
    }
  vardef = &allVars;
  fseek (catalog, 0, SEEK_SET);
  while (FREAD1 (cr, catalog))
    {
/*      if ((!nsub || cr.nsub_ == nsub) && 1 == cr.status) */
        if (  ((!nsub || cr.nsub_ == nsub) && 1 == cr.status) ||
              ((!nsub || cr.nsub_ == nsub) && 4 == cr.status)   ) 
	{ 
/*--------check koeff in model------------------*/   
          if(cr.status == 4)  koeff_check ();    
/*-----------------------------------------------*/
	  seekArchiv1 (fres, cr.factpos);
	  readvardef (fres);
	  for (l = 0; l < vardef->nvar; l++)
	    vararr[vardef->vars[l].num].used = TRUE;
	  clearvardef ();
	  seekArchiv1 (fres, cr.rnumpos);
	  readvardef (fres);
	  for (l = 0; l < vardef->nvar; l++)
	    vararr[vardef->vars[l].num].used = TRUE;
	  clearvardef ();
	  seekArchiv1 (fres, cr.denompos);
	  readDenominators (fres);
	  for (l = 0; l < denrno; l++)
	    {
	      if (denom[l].mass)
		vararr[denom[l].mass].used = TRUE;
	      if (denom[l].width)
		vararr[denom[l].width].used = TRUE;
	    }
	  clearvardef ();
	}
    }
  for (k = nmodelvar; k >= 0; k--)
    if (vararr[k].used && modelvars[k].func)
      readExpression (modelvars[k].func, rd_hiddenVars, act_hiddenVars, NULL);

  kk = 0;
  for (i = 2; i <= numtot; i++)
    for (j = 1; j <= i - 1; j++)
      {
	k = scalarProductPos (i, j);
	switch (key)
	  {
	  case 'R':
	  case 'F':
	    sprintf (vararr[k].alias, "p%d.p%d", j, i);
	    break;
	  case 'M':
	    sprintf (vararr[k].alias, "SC[p%d,p%d]", j, i);
	    break;
	  case 'c':
	    sprintf (vararr[k].alias, "DP[%d]", kk);
	    break;
	  case 'f':
	    {
	      char c = kk + 1;
	      if (c < 10)
		sprintf (vararr[k].alias, "P%c", '0' + c);
	      else
		sprintf (vararr[k].alias, "P%c", 'A' + c - 10);
	    }
	    break;
	  }
	vararr[k].used = TRUE;
	kk++;
      }
  nvar = 0;
  for (k = 0; k <= nmodelvar; k++)
    if (vararr[k].used && !modelvars[k].func)
      {
	vararr[k].tmpvalue = modelvars[k].varvalue;

	switch (key)
	  {
	  case 'R':
	  case 'F':
	  case 'M':
	    strcpy (vararr[k].alias, modelvars[k].varname);
	    break;
	  case 'c':
	    sprintf (vararr[k].alias, "va[%d]", ++nvar);
	    break;
	  case 'f':
	    sprintf (vararr[k].alias, "A(%d)", ++nvar);
	    break;
	  }
      }
  nfunc = 0;
  for (k = 0; k <= nmodelvar; k++)
    if (vararr[k].used && modelvars[k].func)
      {
	vararr[k].tmpvalue = modelvars[k].varvalue;
	switch (key)
	  {
	  case 'R':
	  case 'F':
	  case 'M':
	    strcpy (vararr[k].alias, modelvars[k].varname);
	    break;
	  case 'c':
	    sprintf (vararr[k].alias, "va[%d]", ++nfunc + nvar);
	    break;
	  case 'f':
	    sprintf (vararr[k].alias, "A(%d)", ++nfunc + nvar);
	    break;
	  }
      }

  for (k = 1; k <= nmodelvar; k++)
    if (vararr[k].used && !modelvars[k].able)
      return 0;
  return 1;
}
