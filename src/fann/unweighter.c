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

#include "unweighting_routines.h"

static void print_help (void)
{
  fprintf (stdout, "Help for the unweighter program:\n");
  fprintf (stdout, "Usage: ./unweighter -i {ini_event_file} -o {output_file} [-n {nevents}] [-seed] [-nn]\n");
  fprintf (stdout, "     {ini_event_file} - name of the file with events, by default events.txt\n");
  fprintf (stdout, "     {output_file}    - name of the output file, by default unweighted_events.txt\n");
  fprintf (stdout, " Auxiliary options:\n");
  fprintf (stdout, "     -seed            - says the program to use a seed for pseudo-random number generator from\n");
  fprintf (stdout, "                        the file mix_seed.dat. The file should one long int number and should be\n");
  fprintf (stdout, "                        located in the same directory with event files\n");
  fprintf (stdout, "     -n {nevents}     - number of events processed from the initial event file (by default, all events)\n");
  fprintf (stdout, "     -nn              - prepare LHE file with additional lines with information for NN-code\n");
  fprintf (stdout, "     -nncode          - unweight events according NN values from a code generating by MLPfit\n");
}

int
main (int argc, char **argv)
{
  int i;
  int nevt = -1;
  int frmt;
  int seed_file_used = 0;
  int nnweight_used = 0;
  long seed_used = 123456789;
  char ini_name[1024];
  char out_name[1024];
  FILE *in;
  FILE *ot;

  strcpy (ini_name, "events.txt");
  strcpy (out_name, "unweighted_events.txt");
  for (i = 1; i < argc; ++i) {
    if (!strcmp (argv[i], "-help") || !strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
      print_help ();
      return 0;
    }
    if (!strcmp (argv[i], "-i")) {
      strcpy (ini_name, argv[++i]);
    }
    if (!strcmp (argv[i], "-nn")) {
      nnweight_used = 1;
    }
    if (!strcmp (argv[i], "-nncode")) {
      nnweight_used = 2;
    }
    if (!strcmp (argv[i], "-o")) {
      strcpy (out_name, argv[++i]);
    }
    if (!strcmp (argv[i], "-seed")) {
      seed_file_used = 1;
    }
    if (!strcmp (argv[i], "-n")) {
      int err = sscanf (argv[++i], "%d", &nevt);
      if (err != 1) {
        fprintf(stderr,"unweighter (warning): can't process the ordered number of events (%s)\n", argv[i]);
        fprintf(stderr,"                   and will use all events in the event file.\n");
        nevt = -1;
      }
    }
  }

  in = fopen (ini_name, "r");
  if (!in) {
    fprintf (stderr, "unweighter (error): event file not found: %s\n", ini_name);
    fprintf (stderr, "                 try ./unweighter --help\n");
    return -1;
  }

  ot = fopen (out_name, "a");
  if (!ot) {
    fprintf (stderr, "unweighter (error): can't create output file: %s\n", ini_name);
    fprintf (stderr, "                 try ./unweighter --help\n");
    return -1;
  } else {
    if (0 != ftell (ot)) {
      fprintf(stderr,"unweighter (warning): file %s found and will be rewritten...\n", out_name);
    }
    fclose (ot);
  }

  frmt = CheckFormat (in);
  fclose (in);

  if (seed_file_used) {
    FILE * seed_file = fopen ("mix_seed.dat", "r");
    if (seed_file) {
      long tmp_seed = 0;
      int rep = fscanf (seed_file, "%li", &tmp_seed); 
      if (1 == rep) {
        seed_used = tmp_seed;
      } else {
        fprintf (stderr, "mix (warning): can not read seeed from mix_seed.dat\n");
      }
    } else {
      fprintf (stderr, "mix (warning): seed file mix_seed.dat does not exist. default seed is used\n");
    }
  }

  if (4 != frmt) {
    fprintf (stderr, "unweighter (error): accept event file in the LHE format only");
    fprintf (stderr, "                 use ./traslator to traslate the event file %s to LHE", ini_name);
    return -1;
  }

  fprintf (stdout, "unweighter (info): format LHE detected in the event file\n");
  switch (nnweight_used)
    {
      case 0:
        unweight_events (ini_name, out_name, nevt, seed_used);
        break;
      case 1:
        unweight_events_with_nn_weights (ini_name, out_name, nevt, seed_used);
        break;
      case 2:
        calculate_nn_weights_and_unweight (ini_name, out_name, nevt, seed_used);
        break;
      default:
        fprintf (stdout, "unweighter (error): unknown regime");
        break;
  }

  return 0;
}
