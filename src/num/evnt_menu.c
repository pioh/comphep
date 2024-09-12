/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "chep_crt/include/chep_crt.h"

#include "vegas.h"
#include "evnt_format1.h"
#include "evnt_format2.h"
#include "evnt_format_lhaef.h"
#include "evnt_menu.h"

static float * cubemaxval = NULL;
static int nCubes = 1000;
static long nEvents = 10000;
static double max = 2;
static int nPoints = 500;
static double milk = 0.2;

int ClearEventMax (void) {
  if (cubemaxval) {
    free (cubemaxval);
  }
  cubemaxval = NULL;

  return 1;
}

int WriteEventMax (FILE * f) {
  int l;

  fprintf (f, "nCub=%d \n", nCubes);
  if (cubemaxval) {
    fprintf (f, "Max:\n");
    for (l = 0; l < nCubes; l++) {
      fprintf (f, "%.1E\n", cubemaxval[l]);
    }
  } else {
    fprintf (f, "!Max\n");
  }
  return 0;
}


int ReadEventMax (FILE * f) {
  int l;
  char buff[10];

  fscanf (f, "nCub=%d\n", &nCubes);
  fscanf (f, "%s", buff);
  if (strcmp (buff, "Max:") == 0) {
    cubemaxval = malloc (nCubes * sizeof (float));
    for (l = 0; l < nCubes; l++) {
      fscanf (f, "%f", cubemaxval + l);
    }
  } else {
    if (cubemaxval) {
      free (cubemaxval);
    }
    cubemaxval = NULL;
  }
  return 0;
}


int WriteEventSettings (FILE * f) {
  fprintf (f, "%d %d %f %f %d", nPoints, simplexOn, milk, max, (int) nEvents);
  return 0;
}


int ReadEventSettings (FILE * f) {
  fscanf (f, "%d %d %lf %lf %ld", &nPoints, &simplexOn, &milk, &max, &nEvents);
  return 0;
}


static int grph_event_generator_frmt1 (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, FILE * iprt) {
  int status = prepare_evfile_frmt1 (vegPtr, func, fname, cubemaxval, nEvents, nCubes, max);
  switch (status) {
    case 1: {
        messanykey (25, 15, "---------------\n Events are not saved");
        complete_evfile_frmt1 (fname, 0, 0);
        break;
      }
    case 0: {
        midstr mess;
        double eff  = get_efficiency ();
        double rmax = get_rmax ();
        double mult = get_multiplicity ();
        double neg  = get_negativity ();
        sprintf (mess, "Statistic\n efficiency: %.1E\nReached max: %.1E\n"
             "Mult. events: %.1E \nNeg.events: %.1E \n---------------\n Accept events?", eff, rmax, mult, neg);
        if (mess_y_n (25, 15, mess)) {
          complete_evfile_frmt1 (fname, 1, nEvents);
          fprintf (iprt, "#Event generation\n");
          fprintf (iprt, "  %ld events are stored in '%s'\n", nEvents, fname);
          fprintf (iprt, "  Statistic: eff. = %.1E\n  Max = %.1E\n  Mult.evt. = %.1E\n  Neg.evt = %.1E\n", eff, rmax, mult, neg);
          fflush (iprt);
        } else {
          complete_evfile_frmt1 (fname, 0, 0);
        }
        break;
      }
    case -1:
      messanykey (5, 15, scat ("The file %s is not an event file \nkept in the requested format", fname));
      break;
    case -5:
      messanykey (5, 15, scat ("Can't create event file with name\n%s", fname));
      break;
  }
  return 0;
}

static int grph_event_generator_frmt2 (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, FILE * iprt) {
  int status = prepare_evfile_frmt2 (vegPtr, func, fname, cubemaxval, nEvents, nCubes, max);
  switch (status) {
    case 1: {
        messanykey (25, 15, "---------------\n Events are not saved");
        complete_evfile_frmt2 (fname, 0, 0, 0.0, 0.0);
        break;
      }
    case 0: {
        midstr mess;
        double eff  = get_efficiency ();
        double rmax = get_rmax ();
        double mult = get_multiplicity ();
        double neg  = get_negativity ();
        sprintf (mess, "Statistic\n efficiency: %.1E\nReached max: %.1E\n"
        "Mult. events: %.1E \nNeg.events: %.1E \n---------------\n Accept events?", eff, rmax, mult, neg);
        if (mess_y_n (25, 15, mess)) {
          complete_evfile_frmt2 (fname, 1, nEvents, mult, rmax);
          fprintf (iprt, "#Event generation\n");
          fprintf (iprt, "  %ld events are stored in '%s'\n", nEvents, fname);
          fprintf (iprt, "  Statistic: eff. = %.1E\n  Max = %.1E\n  Mult.evt. = %.1E\n  Neg.evt = %.1E\n", eff, rmax, mult, neg);
          fflush (iprt);
        } else {
          complete_evfile_frmt2 (fname, 0, 0, 0.0, 0.0);
        }
        break;
      }
    case -1:
      messanykey (5, 15, scat ("The file %s is not an event file \nkept in the requested format", fname));
      break;
    case -5:
      messanykey (5, 15, scat ("Can't create event file with name\n%s", fname));
      break;
  }
  return 0;
}

static int grph_event_generator_lhaef (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, FILE * iprt) {
  int status = prepare_evfile_lhaef (vegPtr, func, fname, cubemaxval, nEvents, nCubes, max);
  switch (status) {
    case 1: {
        messanykey (25, 15, "---------------\n Events are not saved");
        complete_evfile_lhaef (fname, 0, 0, 0.0, 0.0);
        break;
      }
    case 0: {
        midstr mess;
        double eff  = get_efficiency ();
        double rmax = get_rmax ();
        double mult = get_multiplicity ();
        double neg  = get_negativity ();
        sprintf (mess, "Statistic\n efficiency: %.1E\nReached max: %.1E\n"
        "Mult. events: %.1E \nNeg.events: %.1E \n---------------\n Accept events?", eff, rmax, mult, neg);
        if (mess_y_n (25, 15, mess)) {
          complete_evfile_lhaef (fname, 1, nEvents, mult, rmax);
          fprintf (iprt, "#Event generation\n");
          fprintf (iprt, "  %ld events are stored in '%s'\n", nEvents, fname);
          fprintf (iprt, "  Statistic: eff. = %.1E\n  Max = %.1E\n  Mult.evt. = %.1E\n  Neg.evt = %.1E\n", eff, rmax, mult, neg);
          fflush (iprt);
        } else {
          complete_evfile_lhaef (fname, 0, 0, 0.0, 0.0);
        }
        break;
      }
    case -1:
      messanykey (5, 15, scat ("The file %s is not an event file \nkept in the requested format", fname));
      break;
    case -2:
      messanykey (5, 15, scat ("The file %s looks like in the LHE format, \nbut does not have the </LesHouchesEvents> tag \nPlease correct it!", fname));
      break;
    case -5:
      messanykey (5, 15, scat ("Can't create event file with name\n%s", fname));
      break;
  }
  return 0;
}

void menu_EventGenerator (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, FILE * iprt, int init) {
  int mode = 1;
  void * pscr = NULL;
  void * pscr_ = NULL;
  double max = 2;
  double eff0;

  if (init) {
    if (cubemaxval) {
      free (cubemaxval);
    }
    cubemaxval = NULL;
  }

  for (;;) {
    if (NULL == cubemaxval) {
      for (mode = 1; mode != 4;) {
        char strmen[] = "\030"
        " sub-cubes = N1         "
        " calls     = N2         "
        " simplex search    |ONN "
        " Start search of maxima ";
        improveStr (strmen, "N1", "%d", nCubes);
        improveStr (strmen, "N2", "%d", nPoints);
        if (simplexOn)
          improveStr (strmen, "|ONN", "|ON ");
        else
          improveStr (strmen, "|ONN", "|OFF");
        menu1 (54, 11, "Preparing of generator", strmen, "n_prep_gen_*", &pscr_, &mode);
        switch (mode) {
          case 0:
            return;
          case 1:
            correctInt (50, 15, "Number of sub-cubes:", &nCubes, 1);
            nCubes = generateVegasCubs (vegPtr, nCubes);
            break;
          case 2:
            correctInt (50, 15, "Calls for each sub-cube:", &nPoints, 1);
            break;
          case 3:
            simplexOn = !simplexOn;
            break;
          case 4: {
              if (nCubes < 1) {
                nCubes = 1;
              }
              cubemaxval  = malloc (nCubes * sizeof (float));
              if (cubemaxval) {
                if (0 == vegas_max (vegPtr, nCubes, nPoints, func, milk, &eff0, cubemaxval)) {
                  midstr mess;
                  sprintf (mess, "Expected efficiency %f", eff0);
                  fprintf (iprt, "#Max\n  Search in %i cubes (%i points per cube) DONE\n", nCubes, nPoints);
                  fprintf (iprt, "  Expected efficiency %f\n", eff0);
                  messanykey (25, 15, mess);
                  put_text (&pscr_);
                } else {
                  mode = 1;
                }
              } else {
                if (cubemaxval) free (cubemaxval);
                cubemaxval = NULL;
                warnanykey (25, 15, "Not enough memory.\nDecrease the number of sub-cubes");
                mode = 1;
              }
            }
            break;
        }
      }
    }

    for (mode = 1;;) {
      char strmen[] = "\031"
        " Number of events=N1     "
        " MAX*N2                  "
        " Generator (new format)  "
        " Generator (old format)  "
        " Generator (LHA format)  "
        " Weighted eventGenerator "
        " New search of maxima    ";
      improveStr (strmen, "N1", "%d", nEvents);
      improveStr (strmen, "N2", "%.2G", max);
      menu1 (54, 10, "", strmen, "n_gen_*", &pscr, &mode);

      switch (mode) {
        case 0:
          return;
        case 1:
          correctLong (50, 15, "", &nEvents, 1);
          break;
        case 2:
          correctDouble (50, 15, "Increase max in ", &max, 1);
          break;
        case 3:
          grph_event_generator_frmt2 (vegPtr, func, fname, iprt);
          break;
        case 4:
          grph_event_generator_frmt1 (vegPtr, func, fname, iprt);
          break;
        case 5:
          set_weighted_flag (0);
          grph_event_generator_lhaef (vegPtr, func, fname, iprt);
          break;
        case 6:
          set_weighted_flag (1);
          grph_event_generator_lhaef (vegPtr, func, fname, iprt);
          break;
      }
      if (7 == mode) {
        free (cubemaxval);
        cubemaxval = NULL;
        put_text (&pscr);
        break;
      }
    }
  }
}

void 
menu_1to2_EventGenerator (double (*func) (double *, double), char * fname, FILE * iprt) {
  int mode = 1;
  void * pscr = NULL;
  vegasGrid * vegPtr = NULL;

  for (mode = 1;;) {
    char strmen[] = "\030"
      "Number of events=N1     "
      " Generator (new format) "
      " Generator (old format) "
      " Generator (LHA format) ";
    improveStr (strmen, "N1", "%d", nEvents);
    menu1 (54, 10, "", strmen, "n_gen_*", &pscr, &mode);
/*
    if (6 == mode) {
      put_text (&pscr);
      break;
    }
*/
    switch (mode) {
      case 0:
        return;
      case 1:
        correctLong (50, 15, "", &nEvents, 1);
        break;
      case 2:
        grph_event_generator_frmt2 (vegPtr, func, fname, iprt);
        break;
      case 3:
        grph_event_generator_frmt1 (vegPtr, func, fname, iprt);
        break;
      case 4:
        grph_event_generator_lhaef (vegPtr, func, fname, iprt);
        break;
    }
  }
}
