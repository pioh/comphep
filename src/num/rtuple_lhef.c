/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef ROOTused

#include "service2/include/chep_limits.h"
#include "LesHouches.h"
#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_routines.h"
#include "event_reader.h"
#include "rtuple_routines.h"

#include "rtuple_lhef.h"

int
rtuple_cpyth_lhef (char ini_name[], char out_name[], int nevnt)
{
  long nposition;
  int Nevents = 0;
  int err;
  eventUP ev;
  char buff[2048];
  char * stop = NULL;

  FILE * s = fopen (ini_name, "r");

  stop = NULL;
  while (!stop) {
    nposition = ftell (s);
    fgets (buff, 2048, s);
    if (feof (s)) return -4;
    stop = strstr (buff, "<event>");
  }

  err = book_rtuple (out_name);
  if (err) {
   fprintf (stderr, "rtupler (error): can't create a rtuple-file. exit...\n");
   return -1;
  }

  while (0 == getLHAevent (ini_name, s, nposition, &ev)) {
    nposition = ftell (s);
    if (Nevents > 0 && 0 == Nevents % 10000) fprintf (stdout, "rtupler (info): %i events read\n", Nevents);
    fill_event (&ev);
    ++Nevents;
    if (nevnt > 0 && nevnt <= Nevents) {
      break;
    }
  }
  write_rtuple ();
  fclose (s);

  fprintf (stdout, "rtupler (info): %i events have been written to the rtuple %s\n", Nevents, out_name);
  return 0;
}

#else
int
rtuple_cpyth_lhef (char ini_name[], char out_name[], int nevnt)
{
  fprintf (stdout, "rtupler (error): ROOT is not linked. Define ROOTSYS and re-compile CompHEP. Exit.\n");
  return 0;
}
#endif
