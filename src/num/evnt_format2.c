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
#include "vegas.h"
#include "alphas2.h"
#include "core_data.h"
#include "cut.h"
#include "runVegas.h"
#include "pdf.h"
#include "strfun.h"
#include "subproc.h"
#include "evnt_tools.h"
#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_writer.h"
#include "evnt_format2.h"

static FILE * events_;
static int * cChains;
static int nC;
static int cBasisPower;

static tags * init_cap_info (void) {
  int i, j;
  int prtn = 0;
  int cutn = 0;
  int PDFid = 0;
  int PDFgr = 0;
  int nl;
  int in_lim = 2;
  int tot_tag_num;
  double sqrtS;
  Str_fun_Info strfun[2];
  tags *t;

  for (i = 1; i <= nin_ + nout_; i++) {
    char buf[10];
    pinf_ (proces_1.nsub, i, buf, NULL);
    if (buf != NULL)
      prtn++;
  }

  cutn = get_cutn ();
  if (1 == nin_) {
    t = init_cap (prtn + cutn + 6);
  } else {
    t = init_cap (prtn + cutn + 10);
  }

/*  beams and structure function info */
  for (i = 0; i < 2; i++)
    wrt_sf_NF_ (i, &strfun[i]);
  vinf_ (0, NULL, &sqrtS);
  for (i = 0; nin_ != 1 && i < 2; i++) {
    double p, e;
    double rap = get_rapidity ();

    e = (sqrtS * sqrtS + strfun[i].prt_mass * strfun[i].prt_mass - strfun[1 - i].prt_mass * strfun[1 - i].prt_mass) / (2 * sqrtS);
    p = sqrt (e * e - strfun[i].prt_mass * strfun[i].prt_mass);
    p = p * cosh (rap) + e * sinh (rap) * (1 - 2 * i);
    if (p * p < (10.E-10) * sqrtS) p = 0.0;

    t->tag[i] = init_tag (5);
    sprintf (t->tag[i]->tagname, "beam");
    sprintf (t->tag[i]->commands[0], "ID=%i", i + 1);
    sprintf (t->tag[i]->commands[1], "energy=%.5E", p);
    sprintf (t->tag[i]->commands[2], "KF=%i", kfbeam (strfun[i].prt_name));
    sprintf (t->tag[i]->commands[3], "name=\'%s\'", strfun[i].prt_name);
    sprintf (t->tag[i]->commands[4], "mass=%.5E", strfun[i].prt_mass);

    t->tag[i + 2] = init_tag (3 + strfun[i].N_extra_commands);
    sprintf (t->tag[i + 2]->tagname, "strfun");
    sprintf (t->tag[i + 2]->commands[0], "IDbeam=%i", i + 1);
    sprintf (t->tag[i + 2]->commands[1], "name=\'%s\'", strfun[i].pdf_name);
    sprintf (t->tag[i + 2]->commands[2], "version=\'%s\'", strfun[i].version);
    for (j = 0; j < strfun[i].N_extra_commands; j++) {
      sprintf (t->tag[i + 2]->commands[j + 3], "%s", strfun[i].extra_commands[j]);
    }
  }
  tot_tag_num = 4;
  if (1 == nin_) tot_tag_num = 0;

  t->tag[tot_tag_num] = init_tag (7);
  sprintf (t->tag[tot_tag_num]->tagname, "process");
  sprintf (t->tag[tot_tag_num]->commands[0], "ID=1");
  {
    int k;
    shortstr name;
    strcpy (name, proces_1.proces);
    for (k = 0; '0' != name[k]; ++k) {
      if (',' == name[k]) name[k] = ' ';
    }
    sprintf (t->tag[tot_tag_num]->commands[1], "name=\'%s\'", name);
  }
  sprintf (t->tag[tot_tag_num]->commands[2], "CrosSec=%.5E", 0.0);
  sprintf (t->tag[tot_tag_num]->commands[3], "CrosSecErr=%.5E", 0.0);
  sprintf (t->tag[tot_tag_num]->commands[4], "Nparton=%i", prtn);
  sprintf (t->tag[tot_tag_num]->commands[5], "master=3");
  if (1 == nin_) {
    sprintf (t->tag[tot_tag_num]->commands[6], "type=\'decay\'");
  } else {
    sprintf (t->tag[tot_tag_num]->commands[6], "type=\'scattering\'");
  }

  tot_tag_num += 1;
  t->tag[tot_tag_num] = init_tag (3);
  sprintf (t->tag[tot_tag_num]->tagname, "generator");
  sprintf (t->tag[tot_tag_num]->commands[0], "IDprocess=1");
  sprintf (t->tag[tot_tag_num]->commands[1], "name=\'CompHEP\'");
  sprintf (t->tag[tot_tag_num]->commands[2], "name=\'%s\'", version ());

  tot_tag_num += 1;
  t->tag[tot_tag_num] = init_tag (5);
  sprintf (t->tag[tot_tag_num]->tagname, "n_event");
  sprintf (t->tag[tot_tag_num]->commands[0], "IDprocess=1");
  sprintf (t->tag[tot_tag_num]->commands[1], "N=%*i", 8, 0);
  sprintf (t->tag[tot_tag_num]->commands[2], "mult=%.5E", 0.0);
  sprintf (t->tag[tot_tag_num]->commands[3], "maxW=%.5E", 0.0);
  sprintf (t->tag[tot_tag_num]->commands[4], "CutN=%i", cutn);

/*      partons info    */
  if (1 == nin_) in_lim = 1;
  for (i = 1; i <= prtn; i++) {
    double part_mass;
    char part_name[16];

    t->tag[i + tot_tag_num] = init_tag (6);
    sprintf (t->tag[i + tot_tag_num]->tagname, "parton");
    sprintf (t->tag[i + tot_tag_num]->commands[0], "IDprocess=1");
    if (i <= in_lim) {
      sprintf (t->tag[i + tot_tag_num]->commands[1], "in=%i", i);
      sprintf (t->tag[i + tot_tag_num]->commands[2], "out=0");
    } else {
      sprintf (t->tag[i + tot_tag_num]->commands[1], "in=0");
      sprintf (t->tag[i + tot_tag_num]->commands[2], "out=%i", i);
    }
    pinf_ (proces_1.nsub, i, part_name, &part_mass);
    sprintf (t->tag[i + tot_tag_num]->commands[3], "KF=%i", kfbeam (part_name));
    sprintf (t->tag[i + tot_tag_num]->commands[4], "name=\'%s\'", part_name);
    sprintf (t->tag[i + tot_tag_num]->commands[5], "mass=%.5E", part_mass);
  }

  tot_tag_num += prtn + 1;
  for (i = 0; i < cutn; i++) {
    vshortstr parts;
    int k = 0;

    if (invcut_1[i].key != 'U') {
      if (invcut_1[i].maxon && invcut_1[i].minon) {
        t->tag[i + tot_tag_num] = init_tag (5);
      } else {
        t->tag[i + tot_tag_num] = init_tag (4);
      }
    } else {
      if (invcut_1[i].maxon && invcut_1[i].minon) {
        t->tag[i + tot_tag_num] = init_tag (6);
      } else {
        t->tag[i + tot_tag_num] = init_tag (5);
      }
    }
    sprintf (t->tag[i + tot_tag_num]->tagname, "cut");
    sprintf (t->tag[i + tot_tag_num]->commands[0], "IDprocess=1");
    j = 0;
    while (invcut_1[i].lvinvc[j]) {
      parts[j] = invcut_1[i].lvinvc[j] + '0';
      j++;
    }
    parts[j] = '\0';
    sprintf (t->tag[i + tot_tag_num]->commands[1], "var=\'%c%s\'", invcut_1[i].key, parts);
    if (invcut_1[i].minon) {
      k++;
      sprintf (t->tag[i + tot_tag_num]->commands[1 + k], "min=%.5E", invcut_1[i].cvmin);
    }
    if (invcut_1[i].maxon) {
      k++;
      sprintf (t->tag[i + tot_tag_num]->commands[1 + k], "max=%.5E", invcut_1[i].cvmax);
    }
    if (invcut_1[i].key != 'U') {
      sprintf (t->tag[i + tot_tag_num]->commands[2 + k], "reg=1");
    } else {
      sprintf (t->tag[i + tot_tag_num]->commands[2 + k], "reg=0");
      sprintf (t->tag[i + tot_tag_num]->commands[3 + k], "filename='userFun.c'");
    }
  }
  tot_tag_num += cutn;

  t->tag[tot_tag_num] = init_tag (4);
  sprintf (t->tag[tot_tag_num]->tagname, "QCDinfo");
  sprintf (t->tag[tot_tag_num]->commands[0], "IDprocess=1");
  CERNpdf_number (strfun[0].pdf_name, strfun[0].version, &PDFid, &PDFgr);
  nl=2;
  if(4==PDFgr && (32==PDFid || 46==PDFid || 58==PDFid)) nl=1;
  sprintf (t->tag[tot_tag_num]->commands[1], "NL=%i",nl);
  sprintf (t->tag[tot_tag_num]->commands[2], "Nflavour=%i", Nflavour());
  sprintf (t->tag[tot_tag_num]->commands[3], "QCDLambda=%.5E", QCDLambda());

  tot_tag_num += 1;
  if (2 == nin_) {
    t->tag[tot_tag_num] = init_tag (3 * prtn);
  } else {
    t->tag[tot_tag_num] = init_tag (3 * prtn + 1);
  }
  sprintf (t->tag[tot_tag_num]->tagname, "format");
  sprintf (t->tag[tot_tag_num]->commands[0], "IDprocess=1");
  sprintf (t->tag[tot_tag_num]->commands[1], "ProcNumber='i'");
  if (2 == nin_) {
    sprintf (t->tag[tot_tag_num]->commands[2], "p1.3='17.10E'");
    sprintf (t->tag[tot_tag_num]->commands[3], "p2.3='17.10E'");
    for (i = 2; i < prtn; i++) {
      sprintf (t->tag[tot_tag_num]->commands[3 * i - 2], "p%i.1='17.10E'", i + 1);
      sprintf (t->tag[tot_tag_num]->commands[3 * i - 1], "p%i.2='17.10E'", i + 1);
      sprintf (t->tag[tot_tag_num]->commands[3 * i], "p%i.3='17.10E'", i + 1);
    }
    sprintf (t->tag[tot_tag_num]->commands[3 * prtn - 2], "Qsquared='10.3E'");
    sprintf (t->tag[tot_tag_num]->commands[3 * prtn - 1], "color_chain='string'");
  } else {
    for (i = 1; i < prtn; i++) {
      sprintf (t->tag[tot_tag_num]->commands[3 * i - 1], "p%i.1='17.10E'", i + 1);
      sprintf (t->tag[tot_tag_num]->commands[3 * i    ], "p%i.2='17.10E'", i + 1);
      sprintf (t->tag[tot_tag_num]->commands[3 * i + 1], "p%i.3='17.10E'", i + 1);
    }
    sprintf (t->tag[tot_tag_num]->commands[3 * prtn - 1], "Qsquared='10.3E'");
    sprintf (t->tag[tot_tag_num]->commands[3 * prtn    ], "color_chain='string'");
  }

  tot_tag_num += 1;
  t->tag[tot_tag_num] = init_tag (4);
  sprintf (t->tag[tot_tag_num]->tagname, "total");
  sprintf (t->tag[tot_tag_num]->commands[0], "Nproc=1");
  sprintf (t->tag[tot_tag_num]->commands[1], "Nevent=       0");
  sprintf (t->tag[tot_tag_num]->commands[2], "CrosSec=%.5E", 0.0);
  sprintf (t->tag[tot_tag_num]->commands[3], "CrosSecErr=%.5E", 0.0);
  return t;
}


static void full_free (tags * t) {
  int i, j;
  for (i = 0; i < t->number_of_tags; i++) {
    for (j = 0; j < t->tag[i]->tagsize; j++) {
      if (t->tag[i]->commands[j]) free (t->tag[i]->commands[j]);
    }
    if (t->tag[i]) free (t->tag[i]);
  }
  free (t);
}


static void write_event_frmt2 (long cCub, int n, double w) {
  int i;
  int k;
  int icc;
  int ntot_ = nin_ + nout_;
  double f;

  if (cBasisPower) {
    double sum = 0;
    for (i = 0; i < cBasisPower; i++) {
      sum += fabs (color_weights[i]);
    }
    sum *= drandXX ();
    for (i = 0; i < cBasisPower; i++) {
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

  f = calcCutFactor ();
  if (!f) {
    fprintf (stderr, "Error : generator have written an event \n");
    fprintf (stderr, "        which does not pass cuts.\n");
    printDetailsOfCutFactor ();
  }

  for (k = 0; k < n; k++) {
    rnd_rotate_momentum (nin_, nout_);
    if (2 == nin_) {
      fprintf (events_, " 1:%17.10E:%17.10E", pvect[3], pvect[7]);
    } else {
      fputs (" 1", events_);
    }
    for (i = nin_; i < ntot_; ++i) {
      fprintf (events_, ":%17.10E:%17.10E:%17.10E", pvect[4 * i + 1], pvect[4 * i + 2], pvect[4 * i + 3]);
    }
    fprintf (events_, ":%10.3E:", qcd_Scale_chep ());
    if (cBasisPower) {
      int j;
      for (j = 0; j < nC; j++) {
        fprintf (events_, "(%d %d)", cChains[2 * (nC * icc + j)], cChains[2 * (nC * icc + j) + 1]);
      }
    }
    fprintf (events_, ":\n");
  }
}

static int fileEnd;
static tags * cap_info;

int prepare_evfile_frmt2 (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, 
   float * cubemaxval, int n_event, int n_cube, double max) {
  int status = 0;
  long fileCur;

  cStrings (proces_1.nsub, &nC, &cBasisPower, &cChains);
  if (cBasisPower) {
    color_weights = malloc (sizeof (double) * cBasisPower);
  }

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
    events_ = fopen (fname, "a+");
    fseek (events_, 0L, SEEK_END);
    if (0 == ftell (events_)) {
      rewind (events_);
      cap_info = init_cap_info ();
      cap_writer (events_, cap_info);
      if (1 == nin_ && 2 == nout_) {
        status = vegas_1to2_events (vegPtr, n_cube, n_event, max, func, write_event_frmt2, cubemaxval);
      } else {
        status = vegas_events (vegPtr, n_cube, n_event, max, func, write_event_frmt2, cubemaxval);
      }
    } else {
      int chck = CheckFormat (events_);
      if (2 != chck) {
        status = -1;
      } else {
        rewind (events_);
        cap_info = init_cap (1);
        cup_reader (events_, cap_info);
        fileCur = ftell (events_);
        fseek(events_, 0L, SEEK_END);
        fileEnd = ftell (events_);
        fseek(events_, fileCur, SEEK_SET);
        if (1 == nin_ && 2 == nout_) {
          status = vegas_1to2_events (vegPtr, n_cube, n_event, max, func, write_event_frmt2, cubemaxval);
        } else {
          status = vegas_events (vegPtr, n_cube, n_event, max, func, write_event_frmt2, cubemaxval);
        }
      }
    }
  }
  fclose (events_);
  return status;
}


int complete_evfile_frmt2 (char * fname, int store, int n_event, double mult, double rmax) {
  double cs, cser;

  if (store) {
    vegas_integral in = get_vegas_integral ();
    in.old = 1;
    set_vegas_integral (in);
    if (in.n_it) {
      cs = in.s1 / in.s0;
      cser = sqrt (in.s0) / fabs (in.s1);
    }

    change_cap (cap_info, n_event, mult, rmax, cs, cser);
    events_ = fopen (fname, "r+");
    cap_writer (events_, cap_info);
    fclose (events_);
  } else {
    truncate (fname, fileEnd);
  }
  if (n_event > 10) full_free (cap_info); /* stupid patch!!! */
  return 1;
}
