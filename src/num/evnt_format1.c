/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "service2/include/chep_limits.h"
#include "service2/include/drandXX.h"
#include "service2/include/lbl.h"
#include "service2/include/4_vector.h"
#include "service2/include/kfcodes.h"

#include "out_ext.h"

#include "alphas_menu.h"
#include "core_data.h"
#include "cut.h"
#include "kinaux.h"
#include "runVegas.h"
#include "strfun.h"
#include "subproc.h"
#include "evnt_tools.h"
#include "evnt_format1.h"

static FILE * events_;
static int * cChains;
static int nC;
static int cBasisPower;


static void write_event_frmt1 (long cCub, int n, double w) {
  int i, k;
  int icc;

  if (cBasisPower) {
    double sum = 0;
    for (i = 0; i < cBasisPower; i++)
      sum += fabs (color_weights[i]);
    sum *= drandXX ();
    for (i = 0; i < cBasisPower; i++){
      sum -= fabs (color_weights[i]);
      if (sum <= 0) {
        break;
      }
    }
    if (i == cBasisPower) {
      i--;
    }
    if (color_weights[i] < 0) {
      n *= -1;
    }
    icc = i;
  }

  for (k = 0; k < n; k++){
    rnd_rotate_momentum (nin_, nout_);
    fprintf (events_, "%2d ", 1);
    if (nin_ == 2) {
      fprintf (events_, "   %17.10E %17.10E", pvect[3], pvect[7]);
    }
    for (i = nin_; i < nout_ + nin_; i++) {
      fprintf (events_, " %17.10E %17.10E %17.10E", pvect[4 * i + 1], pvect[4 * i + 2], pvect[4 * i + 3]);
    }
    fprintf (events_, " %10.3E ", qcd_Scale_chep ());

    if (cBasisPower) {
      int j;
      fprintf (events_, "  ");
      for (j = 0; j < nC; j++)
        fprintf (events_, "(%d %d)", cChains[2 * nC * icc + 2 * j], cChains[2 * nC * icc + 2 * j + 1]);
    }

    fprintf (events_, "\n");
  }
}

static void write_event_cap_frmt1 (void) {
  int i;
  double sqrt_S;
  double cross_section;
  shortstr version1;
  double rap = get_rapidity ();
  vegas_integral in = get_vegas_integral ();

  strcpy (version1, getname());
  fprintf (events_, "#%s\n", version1);

  fprintf (events_, "#PROCESS  ");
  for (i = 1; i <= nin_ + nout_; i++) {
    vshortstr buff;
    pinf_ (proces_1.nsub, i, buff, NULL);
    fprintf (events_, " %-4.4s", buff);
    if (i == nin_) {
      fprintf (events_, "->");
    }
  }
  fprintf (events_, "\n");

  fprintf (events_, "#Initial_state\n");

  vinf_ (0, NULL, &sqrt_S);
  fprintf (events_, "  SQRT(S) %E\n  Rapidity(c.m.s) %E\n", sqrt_S, rap);
  wrt_sf__ (events_);

  fprintf (events_, "#MASSES ");
  for (i = 1; i <= nin_ + nout_; ++i) {
    double q;
    pinf_ (proces_1.nsub, i, NULL, &q);
    fprintf (events_, " %.10E", q);
  }
  fprintf (events_, "\n");

  if (in.n_it) {
    cross_section = in.s1 / in.s0;
    fprintf (events_, "#Cross_section(Width) %E\n", cross_section);
  } else {
    fprintf (events_, "#Cross_section(Width) Unknown\n");
  }
  fprintf (events_, "#Number_of_events %10d\n", 0);
  fprintf (events_,  "#----------------------------------------------------------\n");
  fprintf (events_, "#Number_of_subprocesses = 1\n");
  if (in.n_it) {
    fprintf (events_, "#Total_cross_section_(pb) = %E\n", cross_section);
  } else {
    fprintf (events_, "#Cross_section(Width) Unknown\n");
  }
  fprintf (events_, "%s = %-31d\n", "#Events_mixed_and_randomized", 0);
  fprintf (events_, "#Nproc ================== Events ==========================\n");
}

static long fileEnd;

int prepare_evfile_frmt1 (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, 
   float * cubemaxval, int n_event, int n_cube, double max) {
  int status = 0;

  if (n_cube <= 0) {
   status = -2;
  }

  if (n_event <= 0) {
   status = -3;
  }

  if (max <= 0.0) {
   status = -4;
  }

  if (0 == status) {
    fileEnd = 0;
    events_ = fopen (fname, "a+");
    if (NULL != events_) {
      cStrings (proces_1.nsub, &nC, &cBasisPower, &cChains);
      if (cBasisPower) {
        color_weights = malloc (sizeof (double) * cBasisPower);
      }
      fseek (events_, 0L, SEEK_END);
      if (0 == ftell (events_)) {
        write_event_cap_frmt1 ();
        if (1 == nin_ && 2 == nout_ ) {
          status = vegas_1to2_events (vegPtr, n_cube, n_event, max, func, write_event_frmt1, cubemaxval);
        } else {
          status = vegas_events (vegPtr, n_cube, n_event, max, func, write_event_frmt1, cubemaxval);
        }
        fclose (events_);
        if (1 == status) remove (fname);
      } else {
        int chck = CheckFormat(events_);
        if (1 != chck) {
          status = -1;
        } else {
          fseek (events_, 0L, SEEK_END);
          fileEnd = ftell (events_);
          if (1 == nin_ && 2 == nout_ ) {
            status = vegas_1to2_events (vegPtr, n_cube, n_event, max, func, write_event_frmt1, cubemaxval);
          } else {
            status = vegas_events (vegPtr, n_cube, n_event, max, func, write_event_frmt1, cubemaxval);
          }
        }
        fclose (events_);
      }
      if (cBasisPower) {
        free (color_weights);
        color_weights = NULL;
      }
    } else {
      status = -5;
    }
  }
  return status;
}


int complete_evfile_frmt1 (char * fname, int store, int n_event) {
  int ilen;
  shortstr xfmt;
  shortstr NumEvents;
  long nGenerated = 0;

  if (store) {
    long nEvPos = 0;
    long nEvPos_2 = 0;
    vegas_integral in = get_vegas_integral ();

    in.old = 1;
    set_vegas_integral (in);
    events_ = fopen (fname, "r+");
    while (nEvPos_2 == 0) {
      char ch;
      shortstr word;
      do {
        fscanf (events_, "%c", &ch);
      } while (ch != '#');
      fscanf (events_, "%s", word);
      if (strcmp (word, "Number_of_events") == 0)
        nEvPos = ftell (events_);
      if (strcmp (word, "Events_mixed_and_randomized") == 0)
        nEvPos_2 = ftell (events_);
    }
    fseek (events_, nEvPos, SEEK_SET);
    fscanf (events_, "%ld", &nGenerated);
    nGenerated += n_event;
    fseek (events_, nEvPos, SEEK_SET);
    fprintf (events_, " %10ld", nGenerated);
    sprintf (NumEvents, "%ld", nGenerated);
    ilen = 32 - strlen (NumEvents);
    sprintf (xfmt, "%s%d%s", " = %-", ilen, "ld");
    fseek (events_, nEvPos_2, SEEK_SET);
    fprintf (events_, xfmt, nGenerated);
    fclose (events_);
  } else {
    truncate (fname, fileEnd);
  }
  return 1;
}
