/*
* Copyright (C) 2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <bzlib.h>

#include "service2/include/chep_limits.h"
#include "service2/include/drandXX.h"
#include "service2/include/lbl.h"
#include "service2/include/4_vector.h"
#include "chep_crt/include/crt_util.h"

#include "out_ext.h"

#include "alphas_menu.h"
#include "core_data.h"
#include "cut.h"
#include "kinaux.h"
#include "runVegas.h"
#include "strfun.h"
#include "subproc.h"
#include "weighted_events.h"

static int nEvents = 1000000;
static int cmprFactor = 1;

extern int the_generator (double (*func) (double *, double), char fname[], int nevnt, int cfactor)
{
  int i;
  int status = 0;
  FILE * events_ = fopen (fname, "w");

  if (events_) {
    int bzerror;
/* Parameter blockSize100k specifies the block size to be used for compression. 
   It should be a value between 1 and 9 inclusive, and the actual block size used 
   is 100000 x this figure. 9 gives the best compression but takes most memory.
*/
    int blockSize100k = cfactor;
/* Parameter verbosity should be set to a number between 0 and 4 inclusive. 
   0 is silent
*/
    int verbosity = 1;

/* Lower values of workFactor reduce the amount of effort the standard algorithm 
    will expend before resorting to the fallback. You should set this parameter 
    carefully; too low, and many inputs will be handled by the fallback algorithm 
    and so compress rather slowly, too high, and your average-to-worst case 
    compression times can become very large. The default value of 30 gives 
    reasonable behaviour over a wide range of circumstances.
*/
    int workFactor = 30;

/*
    BZFILE * cfile = BZ2_bzWriteOpen(&bzerror, events_, blockSize100k, verbosity, workFactor);
    if (cfile) {
      for (i = 0; i < nevnt; ++i) {
        
      }
    } else {
      status = -2;
    }
*/
  } else {
    status = -1;
  }

  return status;
}

void 
menu_generator_weighted_events (double (*func) (double *, double), char fname[]) {
  int mode = 1;
  void * pscr = NULL;

  for (mode = 1;;) {
    char strmen[] = "\030"
      " Number of events=N1    "
      " Compress factor=N2     "
      " Generate events        ";
    improveStr (strmen, "N1", "%d", nEvents);
    improveStr (strmen, "N2", "%d", cmprFactor);
    menu1 (54, 10, "", strmen, "n_gen_*", &pscr, &mode);

    switch (mode) {
      case 0:
        return;
      case 1:
        correctInt (50, 15, "", &nEvents, 1);
        break;
      case 2:
        correctInt (50, 15, "", &cmprFactor, 1);
        break;
      case 3:
        the_generator (func, fname, nEvents, cmprFactor);
        break;
    }
  }
}
