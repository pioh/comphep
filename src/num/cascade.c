/* 
* Copyright (C) 2008-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ----------------------------------------------------
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

#include "evnt_tools.h"
#include "cascade_cpyth1.h"
#include "cascade_cpyth2.h"
#include "cascade_cpyth_lhef.h"

static void print_help (void)
{
  fprintf (stdout, "Help for the cascade program:\n");
  fprintf (stdout, "Using: ./cascade [-all] -p {prod_event_file} -d {decay_event_file} -o {output_event_file}\n");
  fprintf (stdout, "     {prod_event_file}  - name of the file with production events, by default production.pev\n");
  fprintf (stdout, "     {decay_event_file} - name of the file with decay events, by default decay.pev\n");
  fprintf (stdout, "     {output_event_file} - name of the output file, by default cascade.pev\n");
  fprintf (stdout, "     -all - name of the output file, by default cascade.pev\n");
}

int
main (int argc, char **argv)
{
  int i;
  int frmt1, frmt2;
  int regime = 0;
  char prod_name[1024];
  char deca_name[1024];
  char out_name[1024];

  FILE *d, *p, *o;

  strcpy (prod_name, "production.pev");
  strcpy (deca_name, "decay.pev");
  strcpy (out_name, "cascade.pev");
  for (i = 1; i < argc; ++i) {
    if (!strcmp (argv[i], "-help") || !strcmp (argv[i], "--help")) {
      print_help ();
      return 0;
    }
    if (!strcmp (argv[i], "-d")) {
      ++i;
      strcpy (deca_name, argv[i]);
    }
    if (!strcmp (argv[i], "-p")) {
      ++i;
      strcpy (prod_name, argv[i]);
    }
    if (!strcmp (argv[i], "-o")) {
      ++i;
      strcpy (out_name, argv[i]);
    }
    if (!strcmp (argv[i], "-all")) {
      regime = 1;
    }
  }

  p = fopen (prod_name, "r");
  if (!p) {
    fprintf (stderr, "cascade (error): production file not found: %s\n", prod_name);
    fprintf (stderr, "                 try ./cascade --help\n");
    return -1;
  }
  d = fopen (deca_name, "r");
  if (!d) {
    fprintf (stderr, "cascade (error): decay file not found: %s\n", deca_name);
    fprintf (stderr, "                 try ./cascade --help\n");
    return -2;
  }
  o = fopen (out_name, "r");
  if (o) {
    fprintf(stderr,"cascade (warning): file %s found and will be rewritten...\n", out_name);
    fclose (o);
  }

  frmt1 = CheckFormat (p);
  frmt2 = CheckFormat (d);
  fclose (p);
  fclose (d);

  if (frmt1 == frmt2) {
    switch (frmt1) {
      case 1:
        fprintf (stdout, "cascade (warning): pass the file via mix, and after that via translator. Exit\n");
        break;
      case 2:
        fprintf (stdout, "\ncascade (info): format cpyth2 detected in the event file\n");
//        if (regime) fprintf (stdout, "cascade (info): all resonances will decay\n");
        cascade_cpyth2 (regime, prod_name, deca_name, out_name);
        break;
      case 3:
        fprintf (stdout, "\ncascade (info): format cpyth1 detected in the event file\n");
//        if (regime) fprintf (stdout, "cascade (info): all resonances will decay\n");
        cascade_cpyth1 (regime, prod_name, deca_name, out_name);
        break;
      case 4:
        fprintf (stdout, "\ncascade (info): format lhef detected in the event file\n");
        if (regime) fprintf (stdout, "cascade (info): all resonances will decay\n");
        cascade_cpyth_lhef (regime, prod_name, deca_name, out_name);
        break;
      default:
        fprintf(stderr, "\ncascade (error): unknown format in event files\n");
        return 4;
      }
  } else {
    fprintf(stderr, "cascade (error): unknown format in event files or the formats are different\n");
    return 4;
  }

  return 0;
}
