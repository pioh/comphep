/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Sergey Shichanin
* ---------------------------------------------------
*/
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "service2/include/chep_limits.h"

#include "file_scr.h"
#include "help.h"

char pathtohelp[256] = "";

int 
show_help (char *fname)
{
  char f_name[STRSIZ];
  FILE *f;
  int z[4];
  sprintf (f_name, "%s%s%s", pathtohelp, fname, ".txt");

  f = fopen (f_name, "r");
  if (f == NULL)
    return 0;
  fgets (f_name, STRSIZ, f);
  sscanf (f_name, "%d %d %d", &z[0], &z[1], &z[2]);
  z[2] = z[2] + z[0] + 1;

  showtext (z[0], z[1], z[2], z[1], fname, f);
  fclose (f);
  return 1;
}
