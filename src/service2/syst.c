/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* ---------------------------------------------------
*/
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>

#include "chep_limits.h"
#include "syst.h"

void (*diskerror) (void) = NULL;

int 
f_printf (FILE * fp, char *format,...)
{
  va_list args;
  char dump[STRSIZ];
  int r;
  va_start (args, format);

  vsprintf (dump, format, args);
  va_end (args);
  r = fputs (dump, fp);
  if (r == EOF)
    {
      if (diskerror)
        (*diskerror) ();
      else
        exit (0);
    }
  return r;
}


size_t 
f_write (void *ptr, size_t size, size_t n, FILE * fp)
{
  size_t nn;
  if ((size == 0) || (n == 0))
    return 0;
  nn = fwrite (ptr, size, n, fp);
  if (nn != n)
    {
      if (diskerror)
        (*diskerror) ();
      else
        exit (0);
    }
  return nn;
}


char *
trim (char *p)
{
  int n1 = 0, n2, k = -1;
  n2 = (int) strlen (p) - 1;
  while (!isgraph (p[n1]) && n1 <= n2)
    n1++;
  while (!isgraph (p[n2]) && n1 < n2)
    n2--;
  while (++k < n2 - n1 + 1)
    p[k] = p[k + n1];
  p[k] = '\0';
  return p;
}

void 
lShift (char *s, int l)
{
  int i;
  if (l > 0)
    {
      int m = strlen (s);
      for (i = 0; i <= m - l; i++)
        s[i] = s[i + l];
    }
  if (l < 0)
    {
      i = strlen (s);
      while (i >= 0)
        {
          s[i - l] = s[i];
          i--;
        }
      for (i = 0; i < -l; i++)
        s[i] = ' ';
    }
}


void 
revers (void **list)
{
  void **q, **p, **r;

  if (*list == NULL)
    return;
  r = (void **) (*list);
  q = (void **) (*r);
  *r = NULL;
  while (q != NULL)
    {
      p = (void **) (*q);
      *q = (void *) r;
      r = q;
      q = p;
    }
  *list = (void *) r;
}

long get_seed (char filename[]) {
  int seed = 123456789;
  FILE * seed_file = fopen (filename, "r");

  if (seed_file) {
    long tmp_seed = 0;
    int rep = fscanf (seed_file, "%li", &tmp_seed);
    if (1 == rep) {
      seed = tmp_seed;
    } else {
      fprintf (stderr, "comphep (warning): can not read seeed from %s\n", filename);
    }
  } else {
    fprintf (stderr, "comphep (warning): seed file %s does not exist. default seed is used\n", filename);
  }

  return seed;
}

int get_sf_info (char sf_info[], char * name, int * set, int * mem)
{
  strcpy (name, sf_info);
  *set = 0;
  *mem = 0;

  if (strstr (sf_info, "LHA") || strstr (sf_info, "PDF")) {
    int num;
    int nstr = 0;
    char * pch = strtok (sf_info, ":");
    while (pch != NULL) {
      if (0 == nstr) sprintf (name, "%s:", pch);
      if (1 == nstr) strcat (name, pch);
      if (2 == nstr) if (1 == sscanf (pch, "%d", &num)) *set = num;
      if (3 == nstr) if (1 == sscanf (pch, "%d", &num)) *mem = num;
      pch = strtok (NULL, ":");
      ++nstr;
    }
  }

  return 1;
}
