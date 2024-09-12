/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/syst.h"
#include "evnt_tools.h"
#include "mix_cpyth1.h"
#include "mix_cpyth2.h"
#include "mix_lhef.h"

static void
ShowHelp (void)
{
  fputs (
  "Help for the mixing program:\n"
  "Usage: ./mix [-seed] {name_of_file_with_events_1} ... {name_of_file_with_events_N} -o {name_of_file_with_mixed_events}\n"
  "     name_of_file_with_events_i     - names of input event files\n"
  "     name_of_file_with_mixed_events - name of the output event file with mixed events,\n"
  "                                      by default Mixed.PEV\n"
  "     -seed - says the program to use a seed for pseudo-random number generator from\n"
  "             the file mix_seed.dat. The file should one long int number and should be\n"
  "             located in the same directory with event files\n"
  "     -zrandom - z-axis randomization (use it if you understand what it means)"
  "\n"
  " The routine prints out some statistic information to stdout. So, user can keep the \n"
  " information redirecting stdout to a file: ./mix.sh ... > file.log\n", stdout);
}

int main (int argc, char** argv)
{
  int i, k;
  int nf, kf;
  int frmt;
  int lenth = 0;
  int first = 1;
  int zrandom_used = 0;
  long seed_used = 123456789;
  char * names;
  char target[2048];
  char user_target[2048] = "";

  nf = kf = argc - 1;
  if (argc < 2) {
    fprintf(stderr, "mix: Enter names of files with events as parameters. \n");
    fprintf(stderr, "     Default name for output event file is Mixed.PEV and\n");
    fprintf(stderr, "     can be overwritten with the option -o {new_file_name}\n");
    return 1;
  }

  for (i = 1; i <= nf; ++i) {
    if (!strcmp (argv[i], "-help") || !strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
      ShowHelp ();
      return 0;
    }

    if (!strcmp (argv[i], "-seed")) {
      seed_used = get_seed ("mix_seed.dat");
      continue;
    }

    if (!strcmp (argv[i], "-zrandom")) {
      zrandom_used = 1;
      continue;
    }

    if (!strcmp (argv[i], "-o")) {
      ++i;
      strcpy (user_target, argv[i]);
      continue;
    }

    {
      int tlenth;
      FILE* inFile = fopen(argv[i], "r");
      if (!inFile) {
        fprintf(stderr,"mix (error): file not found: %s\n", argv[i]);
        return -2;
      } else {
        if (first) {
          frmt = CheckFormat (inFile);
          first = 0;
        }
        if (frmt != CheckFormat (inFile)) {
          fprintf(stderr, "mix (error): different formats in the event files!\n");
          return -3;
        }
        tlenth = strlen (argv[i]);
        if (lenth < tlenth) lenth = tlenth;
      }
      fclose (inFile);
    }
  }

  k = 0;
  names = malloc (nf * lenth * sizeof (char));
  for (i = 1; i <= nf; ++i) {
    if (!strcmp (argv[i], "-zrandom") || !strcmp (argv[i], "-seed")) {
      --kf;
      continue;
    }
    if (!strcmp (argv[i], "-o")) {
      ++i;
      kf -= 2;
      continue;
    }
    strcpy (names + k * lenth, argv[i]);
    k++;
  }

  switch (frmt) {
    case 1:
      strcpy (target, "Mixed.cpyth1");
      break;
    case 2:
      strcpy (target, "Mixed.cpyth2");
      break;
    case 4:
      strcpy (target, "Mixed.lhe");
      break;
    }

  if (0 != strlen (user_target)) strcpy (target, user_target);

  srand48 (seed_used);

  switch (frmt) {
    case 1:
      fprintf (stdout, "mix (info): format cpyth1 detected in all event files\n");
      mix_cpyth1 (kf, target, names, lenth);
      break;
    case 2:
      fprintf (stdout, "mix (info): format cpyth2 detected in all event files\n");
      mix_cpyth2 (kf, target, names, lenth);
      break;
    case 3:
      fprintf (stdout, "mix (info): these files are in the PEVLIB format already. Can't mix them...\n");
      break;
    case 4:
      fprintf (stdout, "mix (info): format lhe detected in all event files\n");
      if (zrandom_used) fprintf (stdout, "mix (info): randomisation along the Z-axis applied\n");
      mix_lhef (kf, target, names, lenth, zrandom_used);
      break;
    default:
      fprintf(stderr, "mix (error): unknown format in event files! Exit\n");
      return 4;
    }
  free (names);
  return 0;
}
