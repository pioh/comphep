/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* ------------------------------------------------------
*/
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/unix_utils.h"
#include "service2/include/read_func.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "chep_crt/include/chep_crt.h"

#include "physics.h"
#include "process_core.h"
#include "screen.h"
#include "constr.h"
#include "read_mdl.h"
#include "sos.h"
#include "batch.h"
#include "process.h"
#include "squar.h"
#include "r_code.h"
#include "out_c.h"
#include "out_red.h"
#include "out_math.h"
#include "out_form.h"
#include "form_code.h"
#include "symbolic.h"
#include "m_utils.h"
#include "ini_ses.h"
#include "amplitudes.h"

static int errorcode = 0;
static int c_make (char *prefix, char *suffix, int n);
static void init_stat (int nfiletot);
static void writestatistic (int nfiletot, int filecount);

#ifdef LHAPDF
static void f10_refresh_lhapdf (int x) {
  update_lhapdf_mdl ();
  messanykey (10, 10, " LHAPDF PDF list updated ");
}
#endif

int
main (int argc, char **argv)
{
  int n;
  int exitlevel;
  long end = -1;
  FILE *LOCK;
  FILE *fverion;
  longstr _pathtocomphep;
/*===================================
if (exitlevel == odd) { don't write chep_label, 
                        the same as a exitlevel--}
0 - Fisrt start
2 - Feynman diagrams menu
4 - Feynman diagrams calculation                        
6 - Squared diagrams menu
8 - Squared diagrams calculation
==========================================================*/
/* 0-Start; 1-Restart; 2-Heap Error,3-Edit Model,4-UserBreak */

  void *pscr1 = NULL;
  void *pscr2 = NULL;
  void *pscr3 = NULL;
  void *pscr4 = NULL;
  void *pscr5 = NULL;
  int k1 = 1, k2 = 1, k3 = 1, k4 = 1, k5 = 1;
  int ProcessChoice = 2;
  int c_prog_ok = 0;
  int diagshow = 0;
  int firstdiag = 1;
  int lastdiag = 100000000;
  char *p;
  shortstr theversion;
  longstr _pathtouser;
  midstr pathtoversionfile;
/* menu variables */

  p = getenv ("COMPHEP");
  if (!p)
    {
      finish ("Environment variable COMPHEP is not defined.\n");
    }
  strcpy (_pathtocomphep, p);
  strcpy (_pathtouser, defaultPath);

  for (n = 1; n < argc; n++)
    {
      if (strcmp (argv[n], "-blind") == 0)
        {
          blind = 1;
          inkeyString = argv[n + 1];
          n++;
        }

      if (strcmp (argv[n], "-diagfirst") == 0)
        {
          int tmp;
          int err = sscanf (argv[n + 1], "%d", &tmp);
          if (1 == err)
            firstdiag = tmp;
          n++;
        }

      if (strcmp (argv[n], "-diaglast") == 0)
        {
          int tmp;
          int err = sscanf (argv[n + 1], "%d", &tmp);
          if (1 == err)
            lastdiag = tmp;
          n++;
        }

      if (strcmp (argv[n], "-diagshow") == 0)
        {
          diagshow = 1;
        }

      if (strcmp (argv[n], "-status") == 0)
        {
          int ncalc = 0;
          int ndel = 0;
          int nrest = 0;
          int ntot = 0;
          long nrecord;
          shortstr txt;
          int nncalc, nndel, nnrest;
          FILE *menuqtmp = fopen (MENUQ_NAME, "rb");
          int subproc = restoreent_dump ();
          for (nsub = 1; nsub <= subproc; nsub++)
            {
              rd_menu (menuqtmp, 2, nsub, txt, &nndel, &nncalc, &nnrest,
                       &nrecord);
              ncalc += nncalc;
              ndel += nndel;
              nrest += nnrest;
              ntot += nndel + nncalc + nnrest;
            }
          fclose (menuqtmp);
          fprintf (stdout,
                   "Status: calculated %i diagrams from %i (deleted %i)\n",
                   ncalc, ntot - ndel, ndel);
          exit (0);
        }
    }

  LOCK = fopen ("LOCK", "a");
  if (LOCK != NULL)
    end = ftell (LOCK);
  if (end == 0)
    {
      fprintf (LOCK, "CompHEP is already started here, process number %i\n",
               getpid ());
      fflush (LOCK);
      fclose (LOCK);
    }
  else
    {
      fprintf (stderr,
               "\nError: CompHEP program was already launched and not finished normaly.\n");
      fprintf (stderr,
               "         Remove the LOCK file to start the CompHEP anyway.\n");
//     exit (100);
    }

  sprintf (pathtouser, "%s%c", _pathtouser, d_slash);
  sprintf (pathtocomphep, "%s%c", _pathtocomphep, d_slash);
  sprintf (pathtohelp, "%shelp%c", pathtocomphep, f_slash);
  outputDir = "results/";

#ifdef LHAPDF
  update_lhapdf_mdl ();
#endif

  sprintf (pathtoversionfile, "%sversion", pathtocomphep);
  fverion = fopen (pathtoversionfile, "r");
  if (fverion != NULL)
    {
      fscanf (fverion, "%s", theversion);
    }
  else
    {
      strcpy (theversion, "unknown");
    }
  setversion (theversion);

  strcpy (theversion, getname ());
  start1 (theversion, scat ("%s%s", pathtocomphep, "icon"),
          "comphep.ini;../comphep.ini;[-]comphep.ini");

  fillModelMenu ();

  f3_key[0] = f3_key_prog;
  f3_mess[0] = "Model";
  f3_key[1] = f4_key_prog;
  f3_mess[1] = "Diagrams";
  f3_key[2] = f5_key_prog;
  f3_mess[2] = "Switches";
  f3_key[3] = f6_key_prog;
  f3_mess[3] = "Results";
  f3_mess[4] = "Del";
  f3_mess[5] = "UnDel";
  f3_key[6] = f9_key_prog;
  f3_mess[6] = "Quit";
#ifdef LHAPDF
  f3_key[7]  = f10_refresh_lhapdf;
  f3_mess[7] = "LHA Upd";
#endif

  restoreent (&exitlevel);
  if (!exitlevel & 1)
    cheplabel ();
  exitlevel = 2 * (exitlevel / 2);
  if (exitlevel >= 2)
    {
      int num = getModelNumberSymb ();
      readModelFiles (num, "models");
      loadModel (FALSE);
      modelinfo ();
      processinfo ();
      diagramsinfo ();
      k1 = num;
    }

  if (exitlevel >= 6)
    {
      sq_diagramsinfo ();
    }

  switch (exitlevel)
    {
    case 2:
      goto label_41;
    case 4:
      goto label_41;
    case 6:
      goto label_50;
    case 8:
      goto restart2;
    }

label_20:                       /* Menu2: ModelMenu */
  f3_key[0] = NULL;             /* models */
  f3_key[1] = NULL;             /* diagrams */
  menulevel = 0;
  menuhelp ();
  for (;;)
    {
      k1 = getModelNumberSymb ();
      menu1 (52, 4, "", modelmenu, "s_1", &pscr1, &k1);
      setModelNumberSymb (k1);
      if (k1 == 0)
        {
          if (mess_y_n (56, 4, "Quit session"))
            {                   /* Exit */
              saveent (menulevel);
              goto exi;
            }
        }
      else if (k1 > maxmodel)
        {
          clrbox (1, 4, 55, 18);
          if (makenewmodel ())
            modelinfo ();
          menuhelp ();
        }
      else if (k1 > 0)
        {
          readModelFiles (k1, "models");
          put_text (&pscr1);
          goto label_30;
        }
    }

label_30:                       /*  Menu3: Enter Process  */
  f3_key[0] = NULL;
  f3_key[1] = NULL;


  menulevel = 0;
  modelinfo ();
  k2 = 2;
/* preload, for proper working of beams table edit and 
           check of ability to work in CompHEP in SUSY models*/
  if (!loadModel (FALSE))
    goto label_20;
  do
    {
      menu1 (51, 6, "", "\032"
             " Enter Decay Process      "
             " Enter Scattering Process "
             " Edit Beams Table         "
             " Edit Str. Functions Table"
             " Edit Model               "
             " Delete Model             ", "s_2_*", &pscr2, &k2);
      switch (k2)
        {
        case 0:         /*  'Esc'  */
          goto_xy (1, 1);
          clr_eol ();
          goto label_20;

        case 1:
          ProcessChoice = 1;
          break;
        case 2:
          ProcessChoice = 2;
          break;
        case 3:
          editBeams ();
          break;
        case 4:
          editStrFuns ();
          break;
        case 5:
          editModel (1);
          break;
        case 6:         /*  Delete  */
          if (deletemodel (getModelNumberSymb ()))
            {
              goto_xy (1, 1);
              clr_eol ();
              setModelNumberSymb (1);
              fillModelMenu ();
              goto label_20;
            }
          else
            {
              int num = getModelNumberSymb ();
              readModelFiles (num, "models");
            }
        }
    }
  while (k2 > 2);

  loadModel (FALSE);


  f3_key[0] = NULL;
  f3_key[1] = NULL;
  menulevel = 0;

label_32:
  errorcode = enter_process (ProcessChoice);    /*  Enter process  */

  if (errorcode != 0)           /*  'Esc' pressed  */
    {
      menuhelp ();
      goto label_30;
    }

  errorcode = construct ();     /*  unSquared diagrams  */
  if (errorcode)
    {
      if (blind)
        {
          batch_error ("Processes of this type are absent", 1);
        }
      else
        {
          clrbox (1, 19, 80, 24);
          goto label_32;        /*  processes are absent  */
        }
    }
  else
    {
      if (!blind)
        {
          if (!viewresults (1))
            goto label_32;
        }
      clr_scr (FGmain, BGmain);
      modelinfo ();
      processinfo ();
      diagramsinfo ();
      goto label_41;
    }

label_40:                       /*  Menu4: Feynman diagram manipulation */
  clr_scr (FGmain, BGmain);
  modelinfo ();
  processinfo ();
  diagramsinfo ();

label_41:
  f3_key[0] = f3_key_prog;
  f3_key[1] = NULL;
  menulevel = 2;
  k3 = 1;
  if (diagshow)
    {
      do
        {
          menu1 (56, 4, "", "\026"      /* Number of columns */
                 " View diagrams        ", "s_squa_*", &pscr3, &k3);
          switch (k3)
            {
            case 0:             /*  Esc  */
              warnanykey (12, 12, "Other menus are unavailable!");
            case 1:             /*  unSquared process menu   */
              viewfeyndiag (1);
              break;
            }
        }
      while (k3 != 2);
    }
  else
    {
      do
        {
          menu1 (56, 4, "", "\026"      /* Number of columns */
                 " View diagrams        "
                 " Square diagrams      ", "s_squa_*", &pscr3, &k3);
          switch (k3)
            {
            case 0:             /*  Esc  */
              clrbox (1, 2, 55, 11);
              menuhelp ();
              goto label_32;
            case 1:             /*  unSquared process menu   */
              viewfeyndiag (1);
            }
        }
      while (k3 != 2);
    }

/*  diagrams  */
  if (blind)
    {
      exclude_feyn_diags (".excluded_digrams");
    }
  if (!squaring ())
    goto label_40;              /*  process are absent  */

  clear_tmp ();

  saveent (menulevel);
  restoreent (&exitlevel);

label_50:                       /*  Menu5: Squared diagram manipulation   */

  f3_key[0] = f3_key_prog;
  f3_key[1] = f4_key_prog;

  menulevel = 6;
  clr_scr (FGmain, BGmain);
  modelinfo ();
  processinfo ();
  diagramsinfo ();
  sq_diagramsinfo ();           /*   ????????   */

  k4 = 1;
  for (;;)
    {
      if (!continuetest ())
        {
          menu1 (56, 4, "", "\026"
                 " View squared diagrams"
                 " Symbolic calculations"
                 " REDUCE program       "
                 " Make n_comphep       "
                 " Prepare process.dat  "
                 " Enter new process    "
 /*                " FORM program         "*/
                 , "s_calc_*", &pscr4, &k4);
          switch (k4)
            {
            case 0:
              if (mess_y_n (50, 3, "Return to previous menu?"))
                {
                  goto label_40;
                }
              break;

            case 1:
              viewsqdiagr ();
              put_text (&pscr4);
              break;

            case 2:             /*  Compute all diagrams   */
            restart2:
              f3_key[0] = f3_key_prog;
              f3_key[1] = f4_key_prog;

              menulevel = 6;
              if (blind)
                {
                  exclude_feyn_csdiags (".excluded_digrams");
                }
              calcallproc (firstdiag, lastdiag);
              sq_diagramsinfo ();
              put_text (&pscr4);
              c_prog_ok = 0;
              k4 = 1;
             break;
            case 3:
              mk_reduceprograms ();
              break;
            case 4:
              saveent (menulevel);
              finish ("\nEnd of CompHEP symbolical session.\n");
              return 24;
            case 5:             /*  unSquared process menu   */
              prepare_process_dat ();
              warnanykey (12, 12,
                          "process_gen.dat has been prepared!\nYou can use it with symb_batch.pl");
              break;
            case 6:
              put_text (&pscr4);
              put_text (&pscr3);
              clrbox (1, 2, 55, 11);
              menuhelp ();
              goto label_30;

            case 7:
/*
              mk_formprograms();
              run_formprograms();
*/
              break;
            }
        }
      else
        {
          saveent (6);
label_60: ;
          char menutxt[] = "\026"
                 " View squared diagrams"
                 " Write results        "
                 " C-compiler           "
                 " Build 1-click tarball"
                 " Enter new process    ";

          menu1 (56, 5, "", menutxt, "s_out_*", &pscr4, &k4);
          switch (k4)
            {
            case 0:
              if (mess_y_n (25, 10, "Return to previous menu?\n All files in results will be removed!"))
                {
                  system ("make -s -C results clean");
                  goto label_40;
                }
              break;

            case 1:
              viewsqdiagr ();
              put_text (&pscr4);
              break;
            case 2:
              {
		FILE * TEST_FILE;
		midstr test_file;
		sprintf (test_file, "%s/results/test_file", pathtouser);
                TEST_FILE = fopen(test_file, "w");
		if (!TEST_FILE) {
		  messanykey (25, 15, "Error: results directory is absent\n Can not write down source files. Create it");
		  break;
		}
		fclose (TEST_FILE);
              }
	      k5 = 1;
              menu1 (56, 6, "", "\026"
                     " C code               "
                     " FORM code            "
                     " REDUCE code          "
                     " MATHEMATICA code     ", NULL, &pscr5, &k5);

              switch (k5)
                {
                case 1: {
                  searchrec s;
                  char * n_chep = scat ("%sresults/n_comphep.exe", pathtouser);
                  if (!find_first (n_chep, &s, 0)) {
                    messanykey (25, 15, "There is n_comphep.exe in Results!\n Clean results at first (Use F6)");
                    break;
                  }
                  c_prog ();
                  {
                    int n_diag = get_diagnumber ();
                    c_make ("f", ".c", n_diag);
                  }
                  c_prog_ok = 1;
                  break;
                }
                case 2:
                  makeFormOutput ();
                  break;
                case 3:
                  makeReduceOutput ();
                  break;
                case 4:
                  makeMathOutput ();
                  break;
                }
              put_text (&pscr5);
              break;

            case 3:
              {
                int err;
                searchrec s;
                char * n_chep = scat ("%sresults/n_comphep.exe", pathtouser);
                char reportlog[128];
		if (blind)
		  strcpy (reportlog, ">>symb_batch.log 2>>symb_batch.log");
		else 
		  strcpy (reportlog, ">>make_gui.log 2>>make_gui.log");
                if (find_first (n_chep, &s, 0))
                  {
                    int i;
                    int j = 0;
                    int n_diag = get_diagnumber ();
                    if (!c_prog_ok)
                      {
                        c_prog ();
                        c_prog_ok = 1;
                      }
                    c_make ("f", ".c", n_diag);
                    if (!blind) {
		      init_stat (n_diag + subproc_sq);
                      err = system ("rm -f make_gui.log");
		    }
                     err = system (scat ("bash -c 'make -C results keep_ffiles %s'", reportlog));
                     for (i = 1; i <= n_diag; ++i) {
                      if (escpressed ()) {
                        goto label_60;
                      }
                      err = system (scat ("bash -c 'make -C results f%i.o %s'", i, reportlog));
                      if (!blind) writestatistic (n_diag + subproc_sq, i);
                      if (0 == i % 100) {
                        ++j;
                        err = system (scat ("bash -c 'make -C results f_compile %s'", reportlog));
                        err = system (scat ("bash -c 'mv results/ffile_%i results/ffile_0 %s'", j, reportlog));
                      }
                    }
                    err = system (scat ("bash -c 'make -C results f_compile %s'", reportlog));
                    for (i = 1; i <= subproc_sq; ++i) {
                      if (escpressed ()) {
                        goto label_60;
                      }
                      err = system (scat ("bash -c 'make -C results d%i.o %s'", i, reportlog));
                      if (!blind) writestatistic (n_diag + subproc_sq, n_diag + i);
                    }
                    err = system (scat ("bash -c 'make -C results d_compile %s'", reportlog));
                    err = system (scat ("bash -c 'make -C results all %s'", reportlog));
                 }
                else
                  {
                    fprintf (stdout, "Old n_comphep.exe has been deleted!");
                    system ("touch results/ffile_0");
                    err = system (scat ("bash -c 'make -C results link %s'", reportlog));
                  }

                if (err != 0)
                  {
                    if (blind)
                      {
                        fprintf (stderr, "\nError during link of n_comphep! See symb_batch.log\n");
                        batch_error ("make error", 1);
                      }
                    else
                      messanykey (25, 15, "Error during link of n_comphep!\n See make_gui.log");
                  }
                else
                  {
                    i_w_sess_ ();
                    if (!blind)
                      {
                        int err = system ("cd results;./n_comphep&");
                        if (err)
                          {
                            messanykey (25, 15, "Can not launch n_comphep automatically!");
                          }
                      }
                  }
              }
              break;
            case 4:
              {
                int prepare_oneclick = 1;
                DIR * d = opendir (scat ("%s/oneclickcode", pathtouser));
                if (NULL != d) {
                  closedir (d);
                  prepare_oneclick = oneclickcode_dir ();
                }
                if (prepare_oneclick) {
                  c_prog ();
                  c_make ("f", ".c", get_diagnumber ());
                  int err = system (scat ("bash -c 'make -C %s oneclick >>make_1click.log 2>>make_1click.log'", pathtoresults));
                  if (err) {
                    messanykey (25, 15, "Can not prepare Oneclick tarball...");
                  } else {
                    messanykey (25, 15, "Oneclick tarball has been prepared!");
                  }
                } else {
                  messanykey (25, 15, "oneclick tarball has not been prepared.");
                }
              }
              break;
            case 5:   /* return to Menu3: Enter Process */
              put_text (&pscr4);
              put_text (&pscr3);
              clrbox (1, 2, 55, 11);
              menuhelp ();
              goto label_30;
            }
        }

    }

exi:
  finish ("End of CompHEP symbolical session.");
  return 0;
}

static int
c_make (char *prefix, char *suffix, int n)
{
  int i, j, max;
  int de, ost;
  int base = 100;
  int num = 0;
  char filename[100];
  char source_files[100 * STRSIZ];
  char ffile[STRSIZ];
  FILE *fi;

  ost = n % base;
  if (!ost)
    de = n / base;
  else
    de = 1 + n / base;

  for (j = 0; j < de; j++)
    {
      sprintf (ffile, "results/ffile_%d", j);
      strcpy (source_files, "");
      max = base;
      if (j == de - 1 && ost)
        {
          max = ost;
        }
      for (i = 0; i < max; i++)
        {
          num++;
          sprintf (filename, "%s%i%s ", prefix, num, suffix);
          strcat (source_files, filename);
        }

      fi = fopen (ffile, "w");
      if (!fi)
        {
          fprintf (stderr, "***** c_make: Cannot open ffile=%s. ABORT\n", ffile);
          exit (99);
        }
      fprintf (fi, "%s", source_files);
      fclose (fi);
    }

  return 0;
}

static void
init_stat (int nfiletot)
{
  goto_xy (1, 17);
  scrcolor (Yellow, Blue);
  print (" C Code Compilation \n");
  scrcolor (Red, BGmain);
  print (" Process.................\n");
  print (" Total number of files...\n");
  print (" Current file............\n");
  print (" Total size Compiled.....\n");
  scrcolor (Yellow, Blue);
  print (" Press Esc to stop    ");
  scrcolor (Black, BGmain);
  goto_xy (25, 18);
  print ("%s", getProcessch ());
  goto_xy (25, 19);
  print ("%4u", nfiletot);
  goto_xy (25, 20);
  print ("   0");
  scrcolor (Yellow, BGmain);
  goto_xy (25, 21);
  print ("   1");
  scrcolor (Yellow, BGmain);
}


static void
writestatistic (int nfiletot, int filecount)
{
  scrcolor (Black, BGmain);
  goto_xy (25, 19);
  print ("%4u", nfiletot);
  goto_xy (25, 20);
  print ("%4u", filecount);
  goto_xy (25, 21);
  print ("%2u (%%)", ((filecount * 100) / nfiletot));
}
