/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ---------------------------------------------------
*/
#include <unistd.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "chep_crt/include/chep_crt.h"

#include "symb/include/sos.h"
#include "symb/include/screen.h"
#include "symb/include/read_mdl.h"
#include "symb/include/physics.h"

static void 
f9_key_prog_diag (int x)
{
  if (mess_y_n (15, 15, " Quit session? "))
    {
      finish ("End of CompHEP diagram viewer");
      exit (0);
    }
}

static void 
view_diagrams (void)
{
  int modelNumber;
  int testmodel;
  modelNumber = getModelNumberSymb();
  testmodel = readModelFiles (modelNumber, "models");
  if (!testmodel) {
    warnanykey (12, 12, "Information about model is not available");
    return;
  }
  testmodel = loadModel (1);
  if (!testmodel) {
    warnanykey (12, 12, "Information about model is not available");
    return;
  } else {
   viewfeyndiag (0);
  }
  return;
}

#ifdef _WIN32
int 
main_n (int argc, char **argv) {
#else
int 
main (int argc, char **argv) {
#endif
  int exitlevel;
  char * p;
  midstr _pathtocomphep;
  FILE * fverion;

  blind = 0;

  p = getenv ("COMPHEP");
  if (!p) {
    fprintf (stderr, " Environment variable COMPHEP is not defined.\n");
    exit (-1);
  }
  strcpy (_pathtocomphep, p);
  sprintf (pathtouser, "%s%c", defaultPath, d_slash);
  sprintf (pathtocomphep, "%s%c", _pathtocomphep, d_slash);
  sprintf (pathtohelp, "%shelp%c", pathtocomphep, f_slash);

#ifdef LHAPDF
  {
    FILE *lf = fopen ("../.lhapdfpath", "r");
    if (!lf) lf = fopen ("../../.lhapdfpath", "r");
    if (lf) {
      midstr _pathtolhapdf;
      if (fscanf (lf, "%1023s", _pathtolhapdf) == 1) {
        sprintf (pathtolhapdf, "%s%c", _pathtolhapdf, d_slash);
        setenv ("LHAPDF_DATA_PATH", _pathtolhapdf, 0);
      }
      fclose (lf);
    }
  }
#endif

  f3_key[6]  = f9_key_prog_diag;
  f3_mess[6] = "Quit";

  restoreent (&exitlevel);
  {
    shortstr inifile;
    shortstr theversion;
    midstr pathtoversionfile;

    sprintf (pathtoversionfile, "%sversion", pathtocomphep);
    fverion = fopen(pathtoversionfile, "r");
    if (fverion != NULL) {
      fscanf (fverion, "%s", theversion);
    } else {
      strcpy (theversion, "unknown");
    }
    setversion (theversion);

    strcpy (inifile, "comphep.ini;../comphep.ini");
    strcpy (theversion, getname ());
    start1 (theversion, scat ("%s%s", pathtocomphep, "icon"), inifile);
  }

  view_diagrams ();
  finish ("End of CompHEP diagram viewer");
  return 0;
}
