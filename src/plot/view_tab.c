/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"
#include "service2/include/syst.h"
#include "service2/include/lbl.h"
#include "chep_crt/include/chep_crt.h"

#include "plot.h"

int
main (int argc, char **argv)
{
  char procName[100];
  char xName[60];
  char yName[60];
  double xMin, xMax;
  double *x = NULL;
  double *f = NULL;
  double *df = NULL;
  int dim;
  int i;
  int n;
  FILE *file;
  FILE *fverion;
  longstr icon_name;
  longstr _pathtocomphep;
  shortstr theversion;
  midstr pathtoversionfile;

  char *p = getenv ("COMPHEP");
  if (!p)
    {
      fprintf (stderr, " Environment variable COMPHEP is not defined.\n");
      exit (0);
    }
  strcpy (_pathtocomphep, p);
  n = strlen (_pathtocomphep) - 1;
  while (n >= 0 && _pathtocomphep[n] != f_slash)
    n--;
  _pathtocomphep[n] = '\0';
  sprintf (pathtocomphep, "%s%c", _pathtocomphep, d_slash);
  sprintf (pathtohelp, "%shelp%c", pathtocomphep, f_slash);

  file = fopen (argv[1], "r");
  if (!file)
    {
      fprintf (stderr, "***Error! File not found: %s\n", argv[1]);
      return 2;
    }

  fscanf (file, "%[^\n]", procName);
  trim (procName);


  fscanf (file, "%*[^\"]%*c%[^\"]%*c%*s%lf%*s%lf%*s%d%*c", xName, &xMin,
          &xMax, &dim);
  x = (double *) malloc (dim * sizeof (double));
  f = (double *) malloc (dim * sizeof (double));
  df = (double *) malloc (dim * sizeof (double));

  fscanf (file, "%[^\n]", yName);
  trim (yName);

  n = fscanf (file, "%lf %lf +/- %lf", x, f, df);
  switch (n)
    {
    case 3:
      {
        for (i = 1; i < dim; i++)
          {
            fscanf (file, "%lf %lf +/- %lf", x + i, f + i, df + i);
          }
        break;
      }
    case 2:
      {
        for (i = 1; i < dim; i++)
          {
            fscanf (file, "%lf %lf", x + i, f + i);
          }
        free (df);
        df = NULL;
        break;
      }
    default:
      {
        fprintf (stderr, "\nErorr! Wrong CompHEP tab format!\n");
        return 0;
      }
    }

  sprintf (icon_name, "%sicon", pathtocomphep);

  sprintf (pathtoversionfile, "%sversion", pathtocomphep);
  fverion = fopen (pathtoversionfile, "r");
  if (fverion != NULL)
    {
      fscanf (fverion, "%s", theversion);
    }
  else
    {
      strcpy (theversion, "unknown");
    }
  setversion (theversion);

  strcpy (theversion, getname ());
  start1 (theversion, icon_name, "comphep.ini;../comphep.ini;[-]comphep.ini");
  clearTypeAhead ();
  plot_histo (xMin, xMax, dim, f, df, procName, xName, yName);
  finish ("");
  return 0;
}
