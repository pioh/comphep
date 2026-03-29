/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
* $Id$
*
* $Log$
*/
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/files.h"

#include "evnt_tools.h"
#include "pdf_reweight_lhef.h"

static void
ShowHelp (void)
{
  fputs (
  "Help for the PDF_reweighter program:\n"
  "Usage: ./pdf_reweighter [-seed] -c {name_of_config_file} -i {name_of_initial_event_file} -o {name_of_file_with_reweighted_events} [-nalphas N]\n"
  "     name_of_config_file                 - name of config event file with infom on both PDFs\n"
  "                                           by config.reweigh\n"
  "     name_of_initial_event_file          - name of input event file\n"
  "                                           by default events.lhe\n"
  "     name_of_file_with_reweighted_events - name of the output event file with reweighted events,\n"
  "                                           by default reweighted.lhe\n"
  "     -nalphas N                          - the number of alpha_s in the Matrix Element,\n"
  "                                           by default N = 0\n"
  "     -seed - says the program to use a seed for pseudo-random number generator from\n"
  "             the file mix_seed.dat. The file should one long int number and should be\n"
  "             located in the same directory with event file\n"
  "\n"
  " The routine prints out some statistic information to stdout. So, user can keep the \n"
  " information redirecting stdout to a file: ./pdf_reweighter.sh ... > file.log\n", stdout);
}

static long get_seed (void) {
  int seed = 123456789;
  FILE * seed_file = fopen ("mix_seed.dat", "r");

  if (seed_file) {
    long tmp_seed = 0;
    int rep = fscanf (seed_file, "%li", &tmp_seed); 
    if (1 == rep) {
      seed = tmp_seed;
    } else {
      fprintf (stderr, "mix (warning): can not read seeed from mix_seed.dat\n");
    }
  } else {
    fprintf (stderr, "mix (warning): seed file mix_seed.dat does not exist. default seed is used\n");
  }

  return seed;
}

int main (int argc, char** argv)
{
  int i;
  int frmt;
  long seed_used = 123456789;
  int nalphas = 0;
  int error_used = 0;
  char target[2048];
  char source[2048];
  char config[2048];
  FILE * iniFile;
  FILE * outFile;
  FILE * conFile;

  strcpy (target, "reweighted.lhe");
  strcpy (source, "events.lhe");
  strcpy (config, "config.reweight");

  for (i = 1; i < argc; ++i) {
    if (!strcmp (argv[i], "-help") || !strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
      ShowHelp ();
      return 0;
    }
    if (!strcmp (argv[i], "-o")) {
      strcpy (target, argv[++i]);
    }
    if (!strcmp (argv[i], "-i")) {
      strcpy (source, argv[++i]);
    }
    if (!strcmp (argv[i], "-c")) {
      strcpy (config, argv[++i]);
    }
    if (!strcmp (argv[i], "-error")) {
      error_used = 1;
    }
    if (!strcmp (argv[i], "-nalphas")) {
      int num;
      if (1 == sscanf (argv[++i],"%d", &num)) {
        nalphas = num;
      } else {
        fprintf (stderr,"pdf_reweighter (warning): can't read the number of alpha_s (option -nalphas), so 0 is used\n");
      }
    }
    if (!strcmp (argv[i], "-seed")) {
      seed_used = get_seed ();
    }
  }

  iniFile = fopen (source, "r");
  if (!iniFile) {
    fprintf(stderr,"pdf_reweighter (error): initial event file %s not found\n", source);
    return -1;
  }
  frmt = CheckFormat (iniFile);
  if (4 != frmt) {
    fprintf(stderr, "pdf_reweighter (error): wrong format, the program supports LHE format only\n");
    return -5;
  }
  fclose (iniFile);

  outFile = fopen (target, "r");
  if (outFile) {
    fprintf(stderr,"pdf_reweighter: (warning): final event file %s found and is re-written\n", target);
    fclose (outFile);
  }

  conFile = fopen (config, "r");
  if (!conFile) {
    fprintf(stderr,"pdf_reweighter: (error): config file %s not found\n", config);
    return -3;
  } else {
    fclose (conFile);
    if (read_config (config)) {
      fprintf(stderr,"pdf_reweighter (error): strange config file %s\n", config);
      return -4;
    }
  }

  if (0 == error_used) {
    fprintf(stderr,"pdf_reweighter: (warning): final cross section is not calculated! If you want it, add -error (it will take much time!)\n");
  }

  if (0 == nalphas) {
    fprintf(stderr,"pdf_reweighter: (warning): only PDFs are used in re-weighting! If your process has non-zero power of alpha_s, -nalphas N, where N is the power \n");
  }

#ifdef LHAPDF
  {
    FILE *lf = fopen ("../.lhapdfpath", "r");
    if (!lf) lf = fopen ("../../.lhapdfpath", "r");
    if (lf) {
      midstr _pathtolhapdf;
      if (fscanf (lf, "%1023s", _pathtolhapdf) == 1) {
        sprintf (pathtolhapdf, "%s/", _pathtolhapdf);
        setenv ("LHAPDF_DATA_PATH", _pathtolhapdf, 0);
      }
      fclose (lf);
    }
  }
  pdf_reweight_lhef (config, source, target, nalphas, seed_used, error_used);
#else
  fprintf (stderr, "pdf_reweighter (error) the program can work with LHAPDF only\n");
#endif

  return 0;
}
