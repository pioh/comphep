/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ----------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/files.h"
#include "service2/include/paragraphs.h"
#include "service2/include/drandXX.h"
#include "service2/include/unix_utils.h"
#include "chep_crt/include/file_scr.h"
#include "chep_crt/include/crt_util.h"

#include "constr.h"
#include "out_c.h"
#include "beams.h"
#include "process.h"
#include "physics.h"
#include "model.h"
#include "process_core.h"
#include "ini_ses.h"

table cutTab =
{"*** Table ***", " Cuts  ",
 "  Parameter  |> Min bound <|> Max bound <|> Exclusive <|", NULL};

table regTab = { "*** Table ***", " Regularization ",
  " Momentum    |> Mass  <|> Width <| Power|", NULL};

table histTab = { "*** Table ***", "Distributions",
  "  Parameter  |> Min bound <|> Max bound <|> Rest Frame <|", NULL};

static int
i_wPrcName_ (FILE * f)
{
/* initial write subprocess number/name*/
  midstr subprocname;
  get_first_subprocname (subprocname);
  fprintf (f, "1 (%s)", subprocname);
  return 1;
}


static int
i_wSesNumber_ (FILE * f)
{
/* initial write session number*/
  fprintf (f, "%d", 1);
  return 0;
}


static int
i_wInState_ (FILE * f)
{				/*  initial write IN-state parameters */
  int i;
  fprintf (f, "\n  SQRT(S) %E\n  Rapidity(c.m.s) %E\n", getsqrtS (),
	   getRapidity ());

  for (i = 0; i < 2; ++i)
    {
      if (strstr (beam[i].h.sf_name, "LHA") || strstr (beam[i].h.sf_name, "PDF")) {
        construct_full_sf_name (i);
	fprintf (f, "  StrFun%i: %s\n", i + 1, beam[i].h.full_sf_name);
      } else {
        fprintf (f, "  StrFun%i: %s\n", i + 1, beam[i].h.sf_name);
      }
    }
  return 0;
}


static int
i_wMdlInfo_ (FILE * f)
{
  fprintf (f, "\n");
  scanvars (f, 4);
  return 0;
}


static int
i_wGwidth_ (FILE * f)
{
  fprintf (f, " %i\n", 0);
  return 0;
}

static int
i_wKinScheme_ (FILE * f)
{
  int i, j;
  int numtot = getntot ();
  int numin = getnin ();

  if (2 == numin)
    fprintf (f, "\n12");
  else
    fprintf (f, "\n1");

  for (i = numin + 1; i < numtot; i++)
    {
      fprintf (f, " -> %i , ", i);
      for (j = i + 1; j <= numtot; j++)
	fprintf (f, "%d", j);
      fprintf (f, "\n");
      if (i == numtot - 1)
	return 0;
      for (j = i + 1; j <= numtot; j++)
	fprintf (f, "%d", j);
    }
  fprintf (f, "\n");
  return 0;
}


static int
i_wCut_ (FILE * f)
{
  fprintf (f, "\n");
  writetable0 (&cutTab, f);
  return 0;
}


static int
i_wRegul_ (FILE * f)
{
  fprintf (f, "\n");
  writetable0 (&regTab, f);
  return 0;
}


static int
i_wQCDInfo_ (FILE * f)
{
  double Lambda6 = 0.1185;
  shortstr Scale = "91.187";

#ifdef LHAPDF
  fprintf (f, "Lambda5 = %E  Scale = %s", Lambda6, Scale);
#else
  fprintf (f, "Lambda6 = %E  Scale = %s", Lambda6, Scale);
#endif
  return 0;
}


static int
i_wMCInfo_ (FILE * f)
{
/* initial write MC parameters */
/*  Initial values:
    mcintr_.ncall0  = 10000;
    mcintr_.itmx0   = 5;
    mcintr_.wrtEvnt = 0;
*/
  fprintf (f, "%dx%d", 100000, 5);
  return 0;
}


static int
i_wIntegral_ (FILE * f)
{
/*
    integral.s0 = 0.0
    integral.s1 = 0.0
    integral.s2 = 0.0
    integral.nCallTot = 0
    integral.n_it     = 0
    integral.old      = 0
*/
  fprintf (f, " %.17E %.17E %.17E %d %d %d", 0.0, 0.0, 0.0, 0, 0, 0);
  return 0;
}


static int
i_wHistograms_ (FILE * f)
{
  fprintf (f, "\n");
  writetable0 (&histTab, f);
  return 0;
}


static int
i_wEventSettings_ (FILE * f)
{
/*
 nPoints   = 500;
 nEvents   = 10000;
 max       = 2.0;
 milk      = 0.2;
 simplexOn = 1
*/
  fprintf (f, "%d %d %f %f %d", 500, 1, 0.2, 2.0, 100000);
  return 0;
}

static int
i_wModelNumber (FILE * f)
{
  fprintf (f, "%i\n", getModelNumberSymb ());
  return 0;
}


static int
i_wRandom_ (FILE * f)
{
  fprintf (f, "%s\n", seedXX (NULL));
  return 0;
}


static int
i_wVegasGrid_ (FILE * f)
{
  fprintf (f, " Vegas_grid: dim=%d  size=%d\n", 0, 0);
  return 0;
}


static int
i_wEventMax_ (FILE * f)
{
  fprintf (f, "nCub=1000 \n");
  fprintf (f, "!Max\n");
  return 0;
}


int
i_w_sess_ (void)
{
  long pos;
  FILE *f;
  shortstr fname;

  rw_paragraph wrt_array[17] = {
    {"Subprocess", i_wPrcName_},
    {"Session_number", i_wSesNumber_},
    {"Model_number", i_wModelNumber},
    {"Initial_state", i_wInState_},
    {"Physical_Parameters", i_wMdlInfo_},
    {"Width_scheme", i_wGwidth_},
    {"Kinematical_scheme", i_wKinScheme_},
    {"Cuts", i_wCut_},
    {"Regularization", i_wRegul_},
    {"QCD", i_wQCDInfo_},
    {"Vegas_calls", i_wMCInfo_},
    {"Vegas_integral", i_wIntegral_},
    {"Distributions", i_wHistograms_},
    {"Events", i_wEventSettings_},
    {"Random", i_wRandom_},
    {"VEGAS_Grid", i_wVegasGrid_},
    {"MAX", i_wEventMax_}
  };

  if (getnin () == 1)
    wrt_array[3].rw_command = NULL;

  sprintf (fname, "results%csession.dat", f_slash);
  f = fopen (fname, "a");
  pos = ftell (f);
  if (pos)
    return 0;

  writeParagraphs (f, 17, wrt_array);
  fclose (f);
  return 0;
}
