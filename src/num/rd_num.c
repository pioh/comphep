/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <ctype.h>
#include "out_ext.h"

#include "rd_num.h"

int 
rd_num (char *s, double *p)
{
  int i;
  char varname[10];
  double val;

  if (isdigit (*s))
    {
      sscanf (s, "%lf", p);
      return 1;
    }

  for (i = 0; i < nvar_ + nfunc_; ++i)
    {
      vinf_ (i + 1, varname, &val);
      if (strcmp (varname, s) == 0)
	{
	  *p = val;
	  return 1;
	}
    }

  return 0;
}
