/*
* Copyright (C) 2008-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ------------------------------------------------------
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"

#include "evnt_tools.h"
#include "rtuple_cpyth1.h"
#include "rtuple_cpyth2.h"
#include "rtuple_lhef.h"

#ifdef ROOTused
static void print_help (void)
{
  fprintf (stdout, "Help for the rtupler program:\n\n");
  fprintf (stdout, "Usage: ./rtupler -i {ini_event_file} -o {output_rtuple} [-n {nevents}] [-rootcode]\n");
  fprintf (stdout, "     {ini_event_file} - name of the file with events, by default events.pev\n");
  fprintf (stdout, "     {output_rtuple}  - name of the output rtuple, by default events.root\n");
  fprintf (stdout, " Auxiliary options:\n");
  fprintf (stdout, "     -n {nevents}     - number of events in rtuple, by default all events\n");
  fprintf (stdout, "     -rootcode        - by default, rtupler creates rtuples only. A toy analysis\n");
  fprintf (stdout, "                        ROOT code is copied with the option.\n");
}

int
main (int argc, char **argv)
{
  int i;
  int nevt = -1;
  int frmt;
  int err;
  int copy_analysis_code = 0;
  char ini_name[1024];
  char out_name[1024];
  char scriptfilename[1030];
  midstr _pathtocomphep;
  char * p;
  FILE *in, *ot;
  FILE * scriptfile;

  p = getenv ("COMPHEP");
  if (!p) {
    fprintf (stderr, " Environment variable COMPHEP is not defined.\n");
    exit (-1);
  }
  strcpy (_pathtocomphep, p);
  sprintf (pathtocomphep, "%s%c", _pathtocomphep, d_slash);

#ifdef LHAPDF
  {
    FILE *lf = fopen ("../.lhapdfpath", "r");
    if (!lf) lf = fopen ("../../.lhapdfpath", "r");
    if (lf) {
      midstr _pathtolhapdf;
      if (fscanf (lf, "%1023s", _pathtolhapdf) == 1) {
        sprintf (pathtolhapdf, "%s%c", _pathtolhapdf, d_slash);
        setenv ("LHAPDF_DATA_PATH", _pathtolhapdf, 0);
      }
      fclose (lf);
    }
  }
#endif

  strcpy (ini_name, "events.pev");
  strcpy (out_name, "events.root");
  for (i = 1; i < argc; ++i) {
    if (!strcmp (argv[i], "-help") || !strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
      print_help ();
      return 0;
    }
    if (!strcmp (argv[i], "-i")) {
      ++i;
      strcpy (ini_name, argv[i]);
    }
    if (!strcmp (argv[i], "-o")) {
      ++i;
      strcpy (out_name, argv[i]);
    }
    if (!strcmp (argv[i], "-rootcode")) {
      copy_analysis_code = 1;
    }
    if (!strcmp (argv[i], "-n")) {
      ++i;
      err = sscanf (argv[i], "%d", &nevt);
      if (err != 1) {
        fprintf(stderr,"rtupler (warning): can't process the ordered number of events (%s)\n", argv[i]);
        fprintf(stderr,"                   and will use all events in the event file.\n");
        nevt = -1;
      }
    }
  }

  in = fopen (ini_name, "r");
  if (!in) {
    fprintf (stderr, "rtupler (error): initial event file not found: %s\n", ini_name);
    fprintf (stderr, "                 try ./rtupler --help\n");
    return -1;
  }

  ot = fopen (out_name, "a");
  if (!ot) {
    fprintf (stderr, "rtupler (error): can't create output file: %s\n", ini_name);
    fprintf (stderr, "                 try ./rtupler --help\n");
    return -1;
  } else {
    if (0 != ftell (ot)) {
      fprintf(stderr,"rtupler (warning): file %s found and will be rewritten...\n", out_name);
    }
    fclose (ot);
  }

  frmt = CheckFormat (in);
  fclose (in);

  switch (frmt) {
    case 1:
      fprintf (stdout, "rtupler (warning): pass the file via mix, and after that via translator. Exit\n");
      break;
    case 2:
      fprintf (stdout, "rtupler (info): format cpyth2 detected in the event file\n");
      rtuple_cpyth2 (ini_name, out_name, nevt);
      break;
    case 3:
      fprintf (stdout, "rtupler (info): format cpyth1 detected in the event file\n");
      rtuple_cpyth1 (ini_name, out_name, nevt);
      break;
    case 4:
      fprintf (stdout, "rtupler (info): format lhef detected in the event file\n");
      rtuple_cpyth_lhef (ini_name, out_name, nevt);
      break;
    default:
      fprintf(stderr, "rtupler (error): unknown format in event files\n");
      return 4;
    }

  if (copy_analysis_code) {
    p = getenv ("COMPHEP");
    if (!p) {
      fprintf (stderr, "rtupler (warning): Environment variable COMPHEP is not defined.\n");
      fprintf (stderr, "                   I can not find and copy analysis code. Do it by hand!\n");
    } else {
      int err1 = system (scat ("cp %s/usr/Tchep.C .", p));
      int err2 = system (scat ("cp %s/usr/Tchep.h .", p));
      if (-1 == err1 || -1 == err2) {
        fprintf (stderr, "rtupler (warning): can't copy Tchep.C/Tchep.h to the working directory\n");
      }
    }

    nextFileName (scriptfilename, "root_analysis_", ".C");
    strcat (scriptfilename, ".C");
    scriptfile = fopen (scriptfilename, "w");
    if (scriptfile) {
      fputs ("{\ngROOT->ProcessLine (\".L Tchep.C+\");\n\n", scriptfile);
      fprintf (scriptfile, "Tchep * example = new Tchep(\"%s\");\n", out_name);
      fputs ("example->Loop();\n//  delete example;\n}", scriptfile);
      fprintf (stdout, "rtupler (info): analysis code copied, ROOT script %s created\n", scriptfilename);
    } else {
      fprintf (stderr, "rtupler (warning): can not create ROOT script!");
    }
  }

  return 0;
}

#else

int
main (int argc, char **argv)
{
  fprintf (stdout, "rtupler (error): rtupler has been compiled without CompHEP internal ROOT code\n");
  fprintf (stdout, "                 please follow the procedure:\n");
  fprintf (stdout, "                 make distclean\n");
  fprintf (stdout, "                 ./configure --with-root\n");
  fprintf (stdout, "                 make;make setup\n");

  return 0;
}
#endif
