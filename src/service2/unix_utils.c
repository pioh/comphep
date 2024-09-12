/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Victor Edneral
* ------------------------------------------------------
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>

#include "chep_limits.h"
#include "syst.h"
#include "unix_utils.h"
#include "tptcmac.h"

char f_slash = '/';
char d_slash = '/';
char *defaultPath = ".";

static DIR *dirpff = NULL;
static char nsff[30];
static char esff[10];
static int lnff, leff;

void
init_os ()
{
  char *buf = NULL;
  int buf_size = PATH_MAX;
  for (;;)
    {
      buf = malloc (buf_size);
      assert (buf != NULL);
      if (getcwd (buf, buf_size) == NULL)
	{
	  if (errno == ERANGE)
	    {
	      free (buf);
	      buf_size *= 2;
	    }
	  else
	    {
	      perror ("cannot allocate directory buffer: ");
	      exit (1);
	    }
	}
      else
	{
	  break;
	}
    }
  defaultPath = strdup (buf);
  free (buf);
}

static void
fsplit (char *fname, char *path, char *name, char *ext)
{
  int k = 0, l, m, n;

  trim (fname);

  m = (int) strlen (fname) - 1;
  n = m;
  while (m >= 0 && fname[m] != '.' && fname[m] != f_slash)
    m--;
  if (m >= 0 && fname[m] == '.')
    {
      for (l = m; l <= n; k++, l++)
	ext[k] = fname[l];
      n = --m;
    }
  else
    m = n;
  ext[k] = '\0';

  while (m >= 0 && fname[m] != d_slash && fname[m] != ':')
    m--;
  for (k = 0, l = m + 1; l <= n; k++, l++)
    name[k] = fname[l];
  name[k] = '\0';

  for (k = 0; k <= m; k++)
    path[k] = fname[k];
  path[k] = '\0';
}


int
find_first (char *filename, searchrec * filerec, int attrib)
{
  char ds[STRSIZ];
  int k;

  if (dirpff != NULL)
    closedir (dirpff);

  fsplit (filename, ds, nsff, esff);
  k = (int) strlen (ds);
  ds[k++] = '.';
  ds[k] = '\0';
  lnff = strcmp (nsff, "*") == 0;
  leff = strcmp (esff, ".*") == 0;

  if ((dirpff = opendir (ds)) == NULL)
    return -1;
  return find_next (filerec);
}

int
find_next (searchrec * filerec)
{
  struct dirent *dp;
  char buff[10];
  char ns1[30];
  char es1[10];

  if (dirpff == NULL)
    return -1;
  while ((dp = (struct dirent *) readdir (dirpff)) != NULL)
    if ((dp->d_name)[0] != '.')
      {
	fsplit (dp->d_name, buff, ns1, es1);
	if ((lnff || strcmp (nsff, ns1) == 0)
	    && (leff || strcmp (esff, es1) == 0))
	  {
	    strcpy (filerec->name, dp->d_name);
	    return 0;
	  }
      }
  closedir (dirpff);
  dirpff = NULL;
  return -1;
}


int
chepmkdir (char *p)
{
  searchrec s;

  trim (p);
  if (0 == find_first (p, &s, directory))
    return 1;
  return mkdir (p, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
}
