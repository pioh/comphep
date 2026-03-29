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
#include "service2/include/syst.h"
#include "chep_crt/include/chep_crt.h"

#include "symb/include/sos.h"
#include "symb/include/screen.h"
#include "symb/include/read_mdl.h"

#include "param.h"
#include "mc_menu.h"
#include "width_12.h"
#include "out_ext.h"
#include "rw_sess.h"
#ifdef LHAPDF
#include "lhapdf.h"
#endif

static int oneclick = 0;

static void 
readtext (char *fname)
{
  FILE *txt;

  trim (fname);
  txt = fopen (fname, "r");
  if (txt == NULL)
    {
      warnanykey (10, 10, " File not found ");
      return;
    }
  showtext (1, 1, 80, maxRow (), "", txt);
  fclose (txt);
}


static void 
viewresults_num (void)
{
  int i, k;
  int doserror;
  char f_name[STRSIZ];
  char menustr[2016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.*");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 2000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
        menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
        menustr[k++] = ' ';
      doserror = find_next (&s);
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
        {
          sprintf (f_name, "%.15s", menustr + k * 15 - 14);
          readtext (f_name);
        }
    }
  while (k != 0);
}


static void 
f6_key_prog_num (int x)
{
  viewresults_num ();
}

static void 
f9_key_prog_num (int x)
{
  if (mess_y_n (15, 15, " Quit session? "))
    {
      write_session ();
      finish ("End of CompHEP numerical session.");
      exit (0);
    }
}

#ifdef LHAPDF
static void f10_refresh_lhapdf (int x) {
  update_lhapdf_mdl ();
  messanykey (10, 10, " LHAPDF PDF list updated ");
}
#endif

static void f4_key_prog_num (int x) {
  system (scat ("%sbin/diag_viewer.exe&", pathtocomphep));
  return;
}

static void f5_key_prog_num (int x) {
  system (scat ("%sbin/diag_viewer2.exe&", pathtocomphep));
  return;
}

static void 
n_comphep (void)
{
  char ExitMessage[] = "Quit session?";
  clr_scr (FGmain, BGmain);

  if (nprc_ > 1)
    {
      goto_xy (4, 2);
      clr_eol ();
      scrcolor (Red, BGmain);
      print ("Process: ");
      scrcolor (Black, BGmain);
      print (processch);
      scrcolor (Red, BGmain);
      print (" (%d subprocesses)", nprc_);
    }

  if (calcFunc ())
    {
      warnanykey (15, 15, "Wrong parameters.\n Can not evaluate constraints");
      change_parameter (54, 7);
    }

  do
    {
      if (oneclick) {
        if (nin_ == 1 && nout_ == 2)
          decay12 ();
        else
          oneclick_monte_carlo_menu ();
      } else {
        if (nin_ == 1 && nout_ == 2)
          decay12 ();
        else
          monte_carlo_menu ();
      }
    }
  while (!mess_y_n (15, 15, ExitMessage));

  write_session ();
}

#ifdef _WIN32
int 
main_n (int argc, char **argv) {
#else
int 
main (int argc, char **argv) {
#endif
  int i, n;
  int exitlevel;
  char * p;
  midstr _pathtocomphep;
  FILE * LOCK;
  FILE * fverion;

  for (n = 1; n < argc; n++) {
    if (strcmp (argv[n], "-blind") == 0) {
      blind = 1;
      inkeyString = argv[n + 1];
    }
    if (strcmp (argv[n], "-nonexpert") == 0) {
      oneclick = 1;
    }
  }

  p = getenv ("TEST_COMPHEP");
  if (!p) {
    LOCK = fopen ("LOCK", "a");
/*    LOCK = fopen ("LOCK", "a+");   FOR TESTS!!!!!!!!!!! */
    if (ftell (LOCK) == 0) {
      fprintf (LOCK, "n_comphep is started here, process number %i\n", getpid ());
      fflush (LOCK);
      fclose (LOCK);
    } else {
      fprintf (stderr, "\ncomphep (warning): program was already launched and not finished normally.\n\n");
      fprintf (stderr, "         Remove the results/LOCK file to start the numerical session anyway");
/*      exit (100);*/
    }
  }

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
  update_lhapdf_mdl ();
#endif

  {
    FILE * percfile = fopen (".procent", "w");
    if (NULL != percfile) {
      for (i = 0; i < nprc_; ++i)  fprintf (percfile, "%i: undefined\n", i + 1);
      fclose (percfile);
    }
  }

  f3_key[1] = f4_key_prog_num;
  if (!oneclick) {
    f3_mess[1] = "Diagrams";
    f3_key[2] = f5_key_prog_num;
    f3_mess[2] = "Squared Diagrams";
    f3_key[3]  = f6_key_prog_num;
  }
  f3_mess[3] = "Results";
  f3_key[6]  = f9_key_prog_num;
  f3_mess[6] = "Quit";
#ifdef LHAPDF
  f3_key[7]  = f10_refresh_lhapdf;
  f3_mess[7] = "LHAPDF Update";
#endif

/*  initialization of the session */
  ComposeSubprocessString ();
  restoreent (&exitlevel);
  init_session ();
  read_session ();

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

  n_comphep ();

  finish ("End of CompHEP numerical session.");
  return 0;
}
