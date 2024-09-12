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
#include "service2/include/syst.h"
#include "num/include/evnt_tools.h"

#include "variables.h"

static void print_help (void)
{
  fprintf (stdout, "Help for the fanvarr program:\n");
  fprintf (stdout, "Usage: ./fanvarr -i {ini_event_file} -o {output_file} [-n {nevents}] [-set {setnumber}]\n");
  fprintf (stdout, "     {ini_event_file} - name of the file with events, by default events.txt\n");
  fprintf (stdout, "     {output_file}    - name of the output file, by default fannvars.txt\n");
  fprintf (stdout, " Auxiliary options:\n");
  fprintf (stdout, "     -n {nevents}     - number of events in rtuple, by default all events\n");
  fprintf (stdout, "                        from events.txt are processed\n");
  fprintf (stdout, "     -set {setnumber} - variable set number for the output file, be default 0.\n");
  fprintf (stdout, "     -full            - prepare LHE file with additional line with variable values\n");
}

int
main (int argc, char **argv)
{
  int i;
  int nevt = -1;
  int nvarset = 0;
  int frmt;
  int regime = 0;
  char ini_name[1024];
  char out_name[1024];
  FILE *in;
  FILE *ot;

  strcpy (ini_name, "events.txt");
  strcpy (out_name, "fannvars.txt");
  for (i = 1; i < argc; ++i) {
    if (!strcmp (argv[i], "-help") || !strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
      print_help ();
      return 0;
    }
    if (!strcmp (argv[i], "-i")) {
      strcpy (ini_name, argv[++i]);
    }
    if (!strcmp (argv[i], "-full")) {
      regime = 1;
    }
    if (!strcmp (argv[i], "-o")) {
      strcpy (out_name, argv[++i]);
    }
    if (!strcmp (argv[i], "-n")) {
      int err = sscanf (argv[++i], "%d", &nevt);
      if (err != 1) {
        fprintf(stderr,"fannvar (warning): can't process the ordered number of events (%s)\n", argv[i]);
        fprintf(stderr,"                   and will use all events in the event file.\n");
        nevt = -1;
      }
    }
    if (!strcmp (argv[i], "-set")) {
      int err = sscanf (argv[++i], "%d", &nvarset);
      if (err != 1) {
        fprintf(stderr,"fannvar (warning): can't process the variable set number of events (%s)\n", argv[i]);
        fprintf(stderr,"                   and will use the default value (0)\n");
        nvarset = -1;
      }
    }
  }

  in = fopen (ini_name, "r");
  if (!in) {
    fprintf (stderr, "fannvar (error): event file not found: %s\n", ini_name);
    fprintf (stderr, "                 try ./fannvar --help\n");
    return -1;
  }

  ot = fopen (out_name, "a");
  if (!ot) {
    fprintf (stderr, "fannvar (error): can't create output file: %s\n", ini_name);
    fprintf (stderr, "                 try ./fannvar --help\n");
    return -1;
  } else {
    if (0 != ftell (ot)) {
      fprintf(stderr,"fannvar (warning): file %s found and will be rewritten...\n", out_name);
    }
    fclose (ot);
  }

  frmt = CheckFormat (in);
  fclose (in);

  if (4 != frmt) {
    fprintf (stderr, "fannvar (error): accept event file in the LHE format only");
    fprintf (stderr, "                 use ./traslator to traslate the event file %s to LHE", ini_name);
    return -1;
  }

  fprintf (stdout, "fannvar (info): format LHE detected in the event file\n");
  switch (regime)
    {
      case 0:
        prepare_fann_variables (ini_name, out_name, nevt, nvarset);
        break;
      case 1:
        add_fann_variables (ini_name, out_name, nevt, nvarset);
        break;
      default:
        fprintf (stdout, "fannvar (error): unknown final file regime");
        break;
  }

  return 0;
}
