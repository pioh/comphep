/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "service2/include/syst.h"

#include "evnt_tools.h"
#include "trans_cpyth1.h"
#include "trans_cpyth2.h"

static void
ShowHelp (void)
{
  fprintf (stdout,
  "\n  Using: translator -i {name_of_file_with_cpyth1/2_events} -o {name_of_file_with_lhe_events}\n\n"
  "The program traslates events in the cpyth1/cpyth2 format into events in the LHE format\n");
}

int main (int argc, char** argv)
{
  int i;
  int frmt;
  int hepml = 0;
  midstr source;
  midstr target;
  FILE * inFile;
  FILE * outFile;
  char * p;
  midstr _pathtocomphep;

  for (i = 1; i < argc; ++i) {
    if (!strcmp (argv[i], "-help") || !strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
      ShowHelp ();
      return 0;
    }
  }

  if (3 > argc) {
    fprintf (stderr, "translator (error): Enter file names: input event file and output file (in lhe format)\n");
    return 1;
  }

  p = getenv ("COMPHEP");
  if (!p) {
    fprintf (stderr, " Environment variable COMPHEP is not defined.\n");
    exit (-1);
  }
  strcpy (_pathtocomphep, p);
  sprintf (pathtouser, "%s%c", defaultPath, d_slash);
  sprintf (pathtocomphep, "%s%c", _pathtocomphep, d_slash);

#ifdef LHAPDF
  {
    midstr _pathtolhapdf;
    p = getenv ("LHAPDFPATH");
    if (!p) {
      fprintf (stderr, " Environment variable LHAPDFPATH is not defined.\n");
      exit (-2);
    }
    strcpy (_pathtolhapdf, p);
    sprintf (pathtolhapdf, "%s%c", _pathtolhapdf, d_slash);
  }
#endif

  strcpy (source, "events.txt");
  strcpy (target, "events_lhef.txt");

  for (i = 1; i < argc; ++i) {
    if (strcmp (argv[i], "--hepml") == 0) {
      hepml = 1;
    }
    if (strcmp (argv[i], "-i") == 0) {
      strcpy (source, argv[++i]);
    }
    if (strcmp (argv[i], "-o") == 0) {
      strcpy (target, argv[++i]);
    }
  }

  inFile = fopen (source, "r");
  if (!inFile) {
    fprintf(stderr,"translator (error): file %s not found...\n", source);
    return -2;
  }
  frmt = CheckFormat (inFile);
  fclose (inFile);

  outFile = fopen (target, "r");
  if (outFile) {
    fprintf(stderr,"translator (warning): file %s found and will be rewritten...\n", target);
    fclose (outFile);
  }

  switch (frmt) {
    case 1:
/*      fprintf (stdout, "translator (warning): pass the file via mix, and after that via translator. Exit\n");*/
      fprintf (stdout, "translator (info): format cpyth1 detected in the event file\n");
      translate_cpyth1 (source, target, 0);
      break;
    case 2:
      fprintf (stdout, "translator (info): cpyth2 format detected in the event file\n");
      translate_cpyth2 (source, target);
      break;
    case 3:
      fprintf (stdout, "translator (info): cpyth1 format detected in the event file\n");
      translate_cpyth1 (source, target, 1);
      break;
    case 4:
      fprintf (stdout, "translator (info): lhe format detected and nothing with the file is done\n");
      break;
    case 5:
      fprintf (stdout, "translator (info): CalcHEP format detected in the event file\n");
      translate_calchep (source, target, 0);
      break;
    default:
      fprintf (stderr, "translator (error): unknown format in the event file\n");
      return 4;
    }

  return 0;
}
