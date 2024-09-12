/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 2000, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "chep_crt/include/chep_crt.h"
#include "plot/include/plot.h"
#include "out_ext.h"

#include "param.h"

/* ************************************************* */
/* Physics model parameters menu                     */
/* ************************************************* */

void
show_depend (int x, int y)
{
  void *pscr1 = NULL;
  for (;;)
    {
      void *pscr2 = NULL;
      char name1[20];
      int k1 = 0;

      selectParam (&k1, x, y + 1, 0, 0, 1, "Display dependence", &pscr1);

      if (k1 < 0)
	return;
      vinf_ (k1, name1, NULL);
      for (;;)
	{
	  char name2[20];
	  double val2;
	  void *pscr3 = NULL;
	  double xMin, xMax;
	  int nPoints = 100;
	  int k3;

	  int k2 = 0;
	  selectParam (&k2, x, y + 5, 0, 1, 0, "on parameter", &pscr2);
	  if (k2 < 0)
	    break;
	  vinf_ (k2, name2, &val2);
	  xMin = val2 - fabs (val2) / 10;
	  xMax = val2 + fabs (val2) / 10;
	  for (;;)
	    {
	      char strmen[] = "\030 "
		" x-Min = XXX            "
		" x-Max = YYY            "
		" Npoints = NNN          "
		" Display plot           ";

	      improveStr (strmen, "XXX", "%G", xMin);
	      improveStr (strmen, "YYY", "%G", xMax);
	      improveStr (strmen, "NNN", "%d", nPoints);

	      menu1 (x, y + 9, "Plot design", strmen, "", &pscr3, &k3);
	      if (!k3)
		{
		  break;
		}
	      switch (k3)
		{
		case 1:
		  correctDouble (x, y + 12, "xMin = ", &xMin, 1);
		  break;
		case 2:
		  correctDouble (x, y + 12, "xMax = ", &xMax, 1);
		  break;
		case 3:
		  correctInt (x, y + 12, "nPoints = ", &nPoints, 1);
		  break;
		case 4:
		  if (xMax > xMin && nPoints >= 3 && nPoints <= 150)
		    {
		      double dx = (xMax - xMin) / (nPoints - 1);
		      double f[150];
		      int i, NaN = 0;;

		      for (i = 0; i < nPoints; i++)
			{
			  double x = xMin + i * dx;
			  asgn_ (k2, x);
			  NaN = calcFunc ();
			  if (NaN)
			    {
			      char mess[100];
			      sprintf (mess, " NaN is evaluated for %s=%G",
				       name2, x);
			      warnanykey (16, 5, mess);
			      break;
			    }
			  vinf_ (k1, NULL, f + i);
			}
		      asgn_ (k2, val2);
		      if (!NaN)
			{
			  plot_table (xMin, xMax, nPoints, f, NULL,
				  "Plot for Constrain", name2, name1);
			}
		      calcFunc ();
		    }
		  else
		    {
		      messanykey (16, 5,
				  " Correct input is \n"
				  "  xMin<xMax,\n" " 2 < nPoints < 201");
		    }
		  break;
		}
	    }
	}
    }
}


static void 
param_menu (int sqtrS_on, int vars_on, int func_on, char **strmen)
{
  int k;
  double val;
  char name[10];

  int pos;
  int npos = 0;

  if (vars_on)
    npos += nvar_;
  if (func_on)
    npos += nfunc_;

  if (npos == 0)
    {
      *strmen = NULL;
      return;
    }
/*  *strmen = malloc (24 * (npos + 1) + 2); */
  *strmen = malloc (24 * (npos + 2) + 2);
  (*strmen)[0] = '\030';

  pos = 1;
  for (k = !sqtrS_on; k <= nvar_ + nfunc_; ++k)
    {
      if (k == 0 || (vars_on && k <= nvar_) || (func_on && k > nvar_))
        {
          char c = ' ';
          vinf_ (k, name, &val);
          if (k > nvar_ && vars_on)
            c = '*';
          sprintf ((*strmen) + pos, "%c%7s= %-14.5g", c, name, val);
          pos += 24;
        }
    }
 
  char c = ' ';
  val=0;
  sprintf ((*strmen) + pos, "%c  Dummy= %-14.5g", c, val);
  pos += 24;

  (*strmen)[pos] = 0;
}

int WriteConstraints (FILE * f) {
  int i;
  double val;
  char name[10];

 fprintf (f, "\n");
  for (i = 0; i < nfunc_; ++i)
    {
      vinf_ (nvar_ + i + 1, name, &val);
      fprintf (f, "%10s = %.15E\n", name, val);
    }
  return 0;
}

void 
selectParam (int *position, int x, int y, int sqtrS_on, int vars_on, int func_on,
             char *mess, void **pscrPtr)
{
  char *strmen;
  void *pscr;

  if (pscrPtr)
    pscr = *pscrPtr;
  else
    pscr = NULL;

  param_menu (sqtrS_on, vars_on, func_on, &strmen);
  if (strmen)
    {
      menu1 (x, y, mess, strmen, "", &pscr, position);
      free (strmen);

      if (*position == 0)
        {
          if (pscrPtr)
            *pscrPtr = NULL;
          *position = -1;
          return;
        }
      if (pscrPtr)
        *pscrPtr = pscr;
      else
        put_text (&pscr);
      if (sqtrS_on)
        (*position)--;
      if (!vars_on)
        *position += nvar_;
    }
  else
    *position = -1;
}

int change_parameter (int x, int y) {
  double val;
  char name[20];
  int first = 1;
  int returnCode = 0;
  void *pscr = NULL;
  int m;

  do {
    if (!first)
      warnanykey (15, 15, " Wrong parameters\n Can not evaluate constraints");
    for (m = 0; m >= 0;) {
      selectParam (&m, x, y, 0, 1, 0, "Change parameter", &pscr);
      if (m >= 0) {
        vinf_ (m, name, &val);
        strcat (name, " = ");
        if (correctDouble (x, y + 4, name, &val, 1)) {
          asgn_ (m, val);
          returnCode = 1;
        }
      }
    }
    first = 0;
  } while (calcFunc ());
  return returnCode;
}
