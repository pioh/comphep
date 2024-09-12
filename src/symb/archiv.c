/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <unistd.h>
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "service2/include/tptcmac.h"

#include "physics.h"
#include "sos.h"
#include "combine.h"

int main (int argc, char **argv) {
  int i;
  int exitlevel;
  int ndiags;
  int err;
  FILE * fverion;
  shortstr theversion;
  longstr com;
  longstr _pathtocomphep;
  longstr pathtocomphep;
  longstr pathtoversionfile;
  FILE** archives;
  FILE** diagpinfo;
  FILE** diagqinfo;
  FILE** catalogs;
  FILE** menups;
  FILE** menuqs;
  FILE** safefiles;
  char** path;
  char * p;
  char * defaultPath = ".";

  int nf = argc - 1;
  safefiles = malloc(nf * sizeof(FILE*));
  archives = malloc(nf * sizeof(FILE*));
  diagpinfo = malloc(nf * sizeof(FILE*));
  diagqinfo = malloc(nf * sizeof(FILE*));
  catalogs = malloc(nf * sizeof(FILE*));
  menups = malloc(nf * sizeof(FILE*));
  menuqs = malloc(nf * sizeof(FILE*));
  path = malloc(nf * sizeof(char*));

  for (i = 1; i < argc; ++i) {
    path[i - 1] = malloc(strlen(argv[i])*sizeof(char));
    strcpy (path[i - 1], argv[i]);
  }

  p = getenv ("COMPHEP");
  if (!p) {
    fprintf (stderr, "Environment variable COMPHEP is not defined.\n");
    exit(-1);
  }
  strcpy (_pathtocomphep, p);
  sprintf (pathtocomphep, "%s/", _pathtocomphep);

  sprintf (pathtouser, "%s/", defaultPath);
  sprintf (pathtoversionfile, "%sversion", pathtocomphep);
  fverion = fopen(pathtoversionfile, "r");
  if (fverion != NULL) {
    fscanf (fverion, "%s", theversion);
  } else {
    strcpy (theversion, "unknown");
  }
  setversion (theversion);

  sprintf (com, "cp %s %stmp/.0_csproces.tp", scat("%s/tmp/csproces.tp", path[0]), pathtouser);
  system (com);
  for (i = 0; i < nf; ++i)
  {
    safefiles[i] = fopen (scat("%s/tmp/safe", path[i]), "rb");
    archives[i]  = fopen (scat("%s/tmp/archive.bt", path[i]), "rb");
    diagqinfo[i] = fopen (scat("%s/tmp/csproces.tp", path[i]), "r+b");
    diagpinfo[i] = fopen (scat("%s/tmp/proces.tp", path[i]), "rb");
    catalogs[i]  = fopen (scat("%s/tmp/catalog.tp", path[i]), "rb");
    menups[i]    = fopen (scat("%s/tmp/menup.ch", path[i]), "rb");
    menuqs[i]    = fopen (scat("%s/tmp/menuq.ch", path[i]), "rb");
  }

  sprintf (com, "cp %s %stmp/.", scat("%s/tmp/safe", path[0]), pathtouser);
  system (com);
  restoreent (&exitlevel);

  ndiags = compare_menuq (nf, nsub, menuqs);

  if (ndiags) {
    err = combine (nf, menuqs[0], archives, diagqinfo, catalogs);
  } else {
    err = -1;
  }

  for (i = 0; i < nf; ++i)
  {
    fclose(safefiles[i]);
    fclose(archives[i]);
    fclose(diagqinfo[i]);
    fclose(diagpinfo[i]);
    fclose(catalogs[i]);
    fclose(menups[i]);
    fclose(menuqs[i]);
  }

  sprintf (com, "cp %s %stmp/.", scat("%s/tmp/csproces.tp", path[0]), pathtouser);
  system (com);
  sprintf (com, "mv %stmp/.0_csproces.tp %s", pathtouser, scat("%s/tmp/csproces.tp", path[0]));
  system (com);
  sprintf (com, "cp %s %stmp/.", scat("%s/tmp/menup.ch", path[0]), pathtouser);
  system (com);
  sprintf (com, "cp %s %stmp/.", scat("%s/tmp/menuq.ch", path[0]), pathtouser);
  system (com);
  sprintf (com, "cp %s %stmp/.", scat("%s/tmp/proces.tp", path[0]), pathtouser);
  system (com);

  for (i = 1; i < argc; ++i) {
    free(path[i - 1]);
  }

  return err;
}
