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
#include "addcut_lhef.h"

static void print_help (void)
{
  fprintf (stdout, "Help for the addcut program:\n\n");
  fprintf (stdout, "Usage: ./addcut -i {initial_event_file} -o {output_event_file} [-n {nevents}]\n");
  fprintf (stdout, "     {initial_event_file} - name of the input event file, by default events.lhe\n");
  fprintf (stdout, "     {output_event_file}  - name of the output event, by default events_cut.lhe\n");
  fprintf (stdout, " Auxiliary options:\n");
  fprintf (stdout, "     -n {nevents}     - number of events to be processed, by default all events\n");
}

int
main (int argc, char **argv)
{
  int i;
  int nevt = -1;
  int frmt;
  int err;

  char ini_name[1024];
  char out_name[1024];
  midstr _pathtocomphep;
  char * p;
  FILE *in, *ot;

  p = getenv ("COMPHEP");
  if (!p) {
    fprintf (stderr, "addcut (error): environment variable COMPHEP is not defined.\n");
    fprintf (stderr, "                use script addcut.sh, not addcut.exe!\n");
    exit (-1);
  }
  strcpy (_pathtocomphep, p);
  sprintf (pathtocomphep, "%s%c", _pathtocomphep, d_slash);

  strcpy (ini_name, "events.lhe");
  strcpy (out_name, "events_cut.lhe");
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
    if (!strcmp (argv[i], "-n")) {
      ++i;
      err = sscanf (argv[i], "%d", &nevt);
      if (err != 1) {
        fprintf(stderr,"addcut (warning): can't process the ordered number of events (%s)\n", argv[i]);
        fprintf(stderr,"                  and will use all events in the event file.\n");
        nevt = -1;
      }
    }
  }

  in = fopen (ini_name, "r");
  if (!in) {
    fprintf (stderr, "addcut (error): initial event file not found: %s\n", ini_name);
    fprintf (stderr, "                try ./affcut.sh --help\n");
    return -1;
  }

  ot = fopen (out_name, "a");
  if (!ot) {
    fprintf (stderr, "addcut (error): can't create output file: %s\n", ini_name);
    fprintf (stderr, "                try ./addcut.sh --help\n");
    return -1;
  } else {
    if (0 != ftell (ot)) {
      fprintf(stderr,"addcut (warning): file %s found and will be rewritten...\n", out_name);
    }
    fclose (ot);
  }

  frmt = CheckFormat (in);
  fclose (in);

  switch (frmt) {
    case 1:
      fprintf (stdout, "addcut (error): pass the file via mix, and after that via translator. Exit\n");
      break;
    case 2:
      fprintf (stdout, "addcut (error): format cpyth2 is not supported anymore. Pass the file via translator\n");
      break;
    case 3:
      fprintf (stdout, "addcut (error): format cpyth1 is not supported anymore. Pass the file via translator\n");
      break;
    case 4:
      fprintf (stdout, "addcut (info): format lhef detected in the event file\n");
      addcut_lhef (ini_name, out_name, nevt);
      break;
    default:
      fprintf(stderr, "addcut (error): unknown format in the event file %s\n", ini_name);
      return 4;
    }

  return 0;
}
