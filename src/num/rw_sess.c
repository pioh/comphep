/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Slava Ilyin 
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/files.h"
#include "service2/include/drandXX.h"
#include "service2/include/paragraphs.h"
#include "chep_crt/include/crt_util.h"
#include "symb/include/physics.h"
#include "out_ext.h"

#include "cut.h"
#include "kinaux.h"
#include "kininpt.h"
#include "regul.h"
#include "param.h"
#include "alphas_menu.h"
#include "histogram.h"
#include "runVegas.h"
#include "strfun.h"
#include "subproc.h"
#include "evnt_menu.h"
#include "core_data.h"
#include "rw_sess.h"

/*********************************************************/
/* R/W of Statistic string (cross section and its error) */
static int WriteIntegral (FILE * f) {
  vegas_integral in = get_vegas_integral ();
  fprintf (f, " %.17E %.17E %.17E %d %d %d",  
  in.s0, in.s1, in.s2, in.n_it, in.old, (int) in.nCallTot);
  return 0;
}

static int ReadIntegral (FILE * f) {
  vegas_integral in;
  fscanf (f, " %lf %lf %lf %d %d %ld", 
          &(in.s0), &(in.s1), &(in.s2), 
          &(in.n_it), &(in.old), &(in.nCallTot));
  set_vegas_integral (in);
  return 0;
}


/**********************************************/
/* R/W subprocess number and name */
static int WriteSubprocessInfo (FILE * mode) {
  return fprintf (mode, "%d (%s)", proces_1.nsub, proces_1.proces);
}

static int ReadSubprocessInfo (FILE * mode) {
  int status = 0;
  shortstr proc_name;

  fscanf (mode, "%d (%[^)]", &proces_1.nsub, proc_name);
  if (proces_1.nsub > nprc_) {
    fprintf (stderr, "Error! The different number of subprocesses... num_{in file} = %i, nprc_ = %i.\n", proces_1.nsub, nprc_);
    return 2;
  }
  ComposeSubprocessString ();
  if (strcmp (proc_name, proces_1.proces)) {
    fprintf (stderr, "Error! Strange parameter names... name_{in file} = %s, name_{from particle names} = %s.\n", proc_name, proces_1.proces);
    status = 2;
  }
  return status;
}

int ComposeSubprocessString (void) {
  vshortstr buff;
  int pos = 0;
  int i;
  for (i = 1; i <= nin_ + nout_; i++) {
    pinf_ (proces_1.nsub, i, buff, NULL);
    strcpy (proces_1.proces + pos, buff);
    pos += strlen (buff);
    if (i < nin_ + nout_) {
      if (i == nin_)
      {
        strcpy (proces_1.proces + pos, " -> ");
        pos += 4;
      } else {
        strcpy (proces_1.proces + pos, ",");
        pos += 1;
      }
    }
  }
  return 0;
}

static int InitSubprocessString (void) {
  proces_1.nsub = 1;
  ComposeSubprocessString ();
  return 0;
}


/*********************************************************/
/* R/W the Session Number */
static int WriteSessionNumber (FILE * mode) {
  fprintf (mode, "%d", get_nsession ());
  return 0;
}

static int ReadSessionNumber (FILE * mode) {
  int ns;
  int err = fscanf (mode, "%d", &ns);
  if (1 == err) {
    set_nsession (ns);
  } else {
    fprintf (stderr, "comphep (error): can not read the session number\n");
    return -1;
  }
  return 0;
}

static int InitSessionNumber (void) {
  set_nsession (1);
  return 0;
}


/**********************************************/
/* R/W Monte-Carlo parameters */
static int WriteMC (FILE * mode) {
  int ncalls = get_ncalls ();
  int iteras = get_iters ();
  fprintf (mode, "%dx%d", ncalls, iteras);
  return 0;
}

static int ReadMC (FILE * mode) {
  int ncalls;
  int iteras;
  fscanf (mode, "%dx%d", &ncalls, &iteras);
  set_ncalls (ncalls);
  set_iters (iteras);
  return 0;
}


/**********************************************/
/* R/W  */
static int WritePhysParameters (FILE * mode) {
  long i;
  double val;
  vshortstr name;
  fprintf (mode, "\n");
  for (i = 0; i < nvar_; ++i) {
    vinf_ (i + 1, name, &val);
    fprintf (mode, "%10s = %.15E\n", name, val);
  }
  return 0;
}

static int ReadPhysParameters (FILE * mode) {
  static double val;
  long i;
  vshortstr name1;
  vshortstr name2;

  for (i = 1; i <= nvar_; ++i) {
    if (2 != fscanf (mode, "%s = %lf", name1, &val)) {
      fprintf (stderr, "Error! Can not read parameter string...\n");
      return 2;
    }
    vinf_ (i, name2, NULL);
    if (strcmp (name1, name2)) {
      fprintf (stderr, "Error! Strange parameter names... name1 = %s, name2 = %s.\n", name1, name2);
      return 2;
    }
    asgn_ (i, val);
  }
  return 0;
}


/**********************************************/
/* R/W IN-state parameters */
static int WriteInitialState (FILE * mode) {
  double sqrt_S;
  vinf_ (0, NULL, &sqrt_S);
  fprintf (mode, "\n  SQRT(S) %E\n  Rapidity(c.m.s) %E\n", sqrt_S, get_rapidity ());
  wrt_sf__ (mode);
  return 0;
}

static int ReadInitialState (FILE * mode) {
  double sqrt_S, rapidity;
  midstr tmp;
  fgets (tmp, STRSIZ, mode);
  fgets (tmp, STRSIZ, mode);
  sscanf (tmp, "  SQRT(S) %lf", &sqrt_S);
  fgets (tmp, STRSIZ, mode);
  sscanf (tmp, "  Rapidity(c.m.s) %lf", &rapidity);
  set_rapidity (rapidity);
  asgn_ (0, sqrt_S);
  rd_sf__ (mode);
  return 0;
}


/**********************************************/
/* R/W model number */
static int nModel = 0;

static int WriteModelNumber(FILE * mode) {
  fprintf (mode, "%d", nModel);
  return 0;
}

static int ReadModelNumber(FILE * mode) {
  midstr tmp;
  fgets (tmp, STRSIZ, mode);
  sscanf (tmp, "%d", &nModel);
  return 0;
}

int getModelNumber (void) {
  int num = getModelNumberSymb ();
  return num;
}


/**********************************************/
/* R/W initial Random seed */
static int WriteRandom (FILE * f) {
  fprintf (f, "%s\n", seedXX (NULL));
  return 0;
}

static int ReadRandom (FILE * f) {
  static char s[128];
  fscanf (f, "%s", s);
  seedXX (s);
  return 0;
}


/**********************************************/
/* R/W width scheme choice */
static int WriteWidthScheme (FILE * f) {
  fprintf (f, " %i\n", gwidth);
  return 0;
}

static int ReadWidthScheme (FILE * f) {
  fscanf (f, " %i", &gwidth);
  return 0;
}


/**********************************************/
/* Session reader and writer */
int write_session (void) {
  rw_paragraph array[17] =
  {
    {"Subprocess",          WriteSubprocessInfo},
    {"Session_number",      WriteSessionNumber},
    {"Model_number",        WriteModelNumber},
    {"Initial_state",       WriteInitialState},
    {"Physical_Parameters", WritePhysParameters},
    {"Width_scheme",        WriteWidthScheme},
    {"Kinematical_scheme",  WriteKinScheme},
    {"Cuts",                WriteCuts},
    {"Regularization",      WriteRegularisations},
    {"QCD",                 WriteQCDInfo},
    {"Vegas_calls",         WriteMC},
    {"Vegas_integral",      WriteIntegral},
    {"Distributions",       WriteHistograms},
    {"Events",              WriteEventSettings},
    {"Random",              WriteRandom},
    {"VEGAS_Grid",          WriteVegasGrid},
    {"MAX",                 WriteEventMax}
  };
  FILE * f = fopen (scat ("%ssession.dat", outputDir), "w");
  if (f == NULL) {
    return 0;
  }

  if (nin_ == 1) {
    array[3].rw_command = NULL;
  }

  writeParagraphs (f, 17, array);
  fclose (f);
  return 0;
}


int write_prt (FILE * f) {
  rw_paragraph array[12] =
  {
    {"Subprocess",          WriteSubprocessInfo},
    {"Session_number",      WriteSessionNumber},
    {"Model_number",        WriteModelNumber},
    {"Initial_state",       WriteInitialState},
    {"Physical_Parameters", WritePhysParameters},
    {"Constraints",         WriteConstraints},
    {"Width_scheme",        WriteWidthScheme},
    {"Kinematical_scheme",  WriteKinScheme},
    {"Cuts",                WriteCuts},
    {"Regularization",      WriteRegularisations},
    {"QCD",                 WriteQCDInfo},
    {"Vegas_calls",         WriteMC},
  };

  if (nin_ == 1) {
    array[3].rw_command = NULL;
  }
  writeParagraphs (f, 12, array);
  return 0;
}


int read_session (void) {
  rw_paragraph array[17] =
  {
    {"Subprocess",          ReadSubprocessInfo},
    {"Session_number",      ReadSessionNumber},
    {"Model_number",        ReadModelNumber},
    {"Initial_state",       ReadInitialState},
    {"Physical_Parameters", ReadPhysParameters},
    {"Width_scheme",        ReadWidthScheme},
    {"Kinematical_scheme",  ReadKinScheme},
    {"Cuts",                ReadCuts},
    {"Regularization",      ReadRegularisations},
    {"QCD",                 ReadQCDInfo},
    {"Distributions",       ReadHistograms},
    {"Vegas_integral",      ReadIntegral},
    {"Vegas_calls",         ReadMC},
    {"Events",              ReadEventSettings},
    {"Random",              ReadRandom},
    {"VEGAS_Grid",          ReadVegasGrid},
    {"MAX",                 ReadEventMax}
  };
  FILE * f = fopen ("session.dat", "r");
  if (f == NULL) {
    fprintf (stderr, "\nn_comphep (warning): can not open session.dat! Dummy session.dat will be created.\n");
    return 0;
  }
  readParagraphs (f, 17, array);
  fclose (f);
  return 0;
}


int init_session (void) {
  InitSubprocessString ();     /* structure process from the comphep output file */
  InitSessionNumber (); /* session number = 1 */
  InitKinScheme ();
  InitQCDInfo ();
  set_rapidity (0.);
  gwidth = 0;
  return 0;
}


void clearSession (void) {
  cleartab (&cutTab);
  cleartab (&regTab);
  cleartab (&histTab);
  clearHists ();
}
