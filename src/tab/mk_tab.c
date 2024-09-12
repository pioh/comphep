/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <stdio.h>
#include <string.h>

#include "num/include/evnt_tools.h"
#include "num/include/phys_val.h"

#include "e_tools.h"
#include "tab_routines.h"

static void
wrongParam (int N)
{
  if (N) {
    fprintf (stderr, "Wrong parameter %d\n", N);
  } else {
    fprintf (stderr, "Wrong number of parameters \n");
  }

  fprintf (stderr, "Parameters:\n"
           " 1- file name,\n"
           " 2- name of variable,\n"
           " 3- minimum limit,\n"
           " 4- maximum limit,\n" " 5- number of bins(<=300).\n");
}


static void
ShowHelp (int N)
{
  if (N == 0) {
    fprintf (stdout, "Short help:\n"
    "  Using: mk_tab {name_of_file_with_events} Var Param_min Param_max Nbins\n\n"    
    "The program mk_tab constructs a histogram for variable \"Var\" in the limits\n"
    "from \"Param_min\" to \"Param_max\". This histogram will have \"Nbins\" bins\n");
  } else {
    fprintf (stdout, "Long help:\n"
    "Using: mk_tab {name_of_file_with_events} Var Param_min Param_max Nbins\n\n"
    "The program mk_tab constructs a histogram for variable \"Var\" in the limits\n"
    "from \"Param_min\" to \"Param_max\". This histogram will have \"Nbins\" bins"
    "Bla-bla-bla...\n");
  }
}


int
main (int argc, char **argv)
{
  int frmt;
  int nbin;
  double minX, maxX;
  FILE *f;

  if (argc != 6)
    {
      wrongParam (0);
      return 1;
    }
  if (strcmp (argv[1], "") == 0 || strcmp (argv[1], "-h") == 0)
    {
      ShowHelp (0);
      return 1;
    }
  if ((strcmp (argv[1], "-help") == 0) || (strcmp (argv[1], "--help") == 0))
    {
      ShowHelp (1);
      return 1;
    };

  if (sscanf (argv[3], "%lf", &minX) != 1)
    {
      wrongParam (2);
      return 1;
    }
  if (sscanf (argv[4], "%lf", &maxX) != 1 || minX >= maxX)
    {
      wrongParam (3);
      return 1;
    }
  if (sscanf (argv[5], "%d", &nbin) != 1 || nbin <= 0)
    {
      wrongParam (4);
      return 1;
    }

  f = fopen (argv[1], "r");
  if (!f)
    {
      fprintf (stderr, " mk_tab (error): file %s not found\n", argv[1]);
      return -2;
    }

  frmt = CheckFormat (f);
  fclose (f);

  prepare_tab (frmt, argv[1], argv[2], nbin, minX, maxX);

  return 0;
}
