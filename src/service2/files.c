/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* -----------------------------------------------------
*/
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>

#include "chep_limits.h" 
#include "chep_crt/include/chep_crt.h"

#include "unix_utils.h"
#include "syst.h"
#include "files.h"

#define BLOCKSIZE 1024

FILE *menup   = 0x0;
FILE *menuq   = 0x0;
FILE *diagrp  = 0x0;			/* file of Adiagram; */
FILE *diagrq  = 0x0;			/* file of CSdiagram; */
FILE *catalog = 0x0;

char mdFls[5][10] =
{"vars", "func", "prtcls", "lgrng", "cpart"};

char *outputDir = "";
longstr pathtocomphep = "";
longstr pathtolhapdf = "";
longstr pathtoresults = "results";
longstr pathtouser = "";

void 
wrt_menu (FILE * men, int menutype, int k, char *txt, int ndel, int ncalc, int nrest,
	  long recpos)
{
  if (menutype == 1)
    {
      fseek (men, (k - 1) * 54 + 2, SEEK_SET);
      fprintf (men, "%4d| %-27.27s|%5d|%5d|%-8d", k, txt, ndel, nrest, (int) recpos);
      fflush(men);
    }
  else
    {
      fseek (men, (k - 1) * 60 + 2, SEEK_SET);
      fprintf (men, "%4d| %-27.27s|%5d|%5d|%5d|%-8d", k, txt, ndel, ncalc, nrest,
              (int) recpos);
      fflush(men);
    }
}


int 
rd_menu (FILE * men, int menutype, int k, char *txt, int *ndel, int *ncalc, int *nrest, long *recpos)
{
  if (menutype == 1)
    {
      fseek (men, (k - 1) * 54 + 2, SEEK_SET);
      if (5 != fscanf (men, "%d| %[^|]%*c%d|%d|%ld", &k, txt, ndel, nrest, recpos))
	return 0;
      *ncalc = 0;
    }
  else
    {
      fseek (men, (k - 1) * (60) + 2, SEEK_SET);
      if (6 != fscanf (men, "%4d| %[^|]%*c%d|%d|%d|%ld", &k, txt, ndel, ncalc, nrest, recpos))
	return 0;
    }
  return 1;
}

void 
copyfile (char *namefrom, char *nameto)
{
  FILE *filefrom;
  FILE *fileto;
  const int _blocksize = BLOCKSIZE;
  char s[BLOCKSIZE];
  int size;
  filefrom = fopen (namefrom, "rb");
  fileto = fopen (nameto, "wb");
  if ((filefrom == NULL) || (fileto == NULL))
    return;

  do
    {
      size = fread (s, 1, _blocksize, filefrom);
      fwrite (s, 1, size, fileto);
    }
  while (size == _blocksize && size != 0);

  fclose (fileto);
  fclose (filefrom);
  return;
}

void 
nextFileName (char *f_name, char *firstname, char *ext)
{
  int tabnum = 0;
  midstr lname;
  FILE *f;

  for (;;)
    {
      tabnum++;
      sprintf (f_name, "%s%s%d", outputDir, firstname, tabnum);
      sprintf (lname, "%s%s", f_name, ext);
      f = fopen (lname, "r");
      if (f)
	fclose (f);
      else
	return;
    }
}
