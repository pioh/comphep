/*
* Copyright (C) 2008-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"
#include "service2/include/tptcmac.h"

#include "physics.h"
#include "process_core.h"
#include "diagrams.h"
#include "batch.h"

static char stc[128];
static char ddel1[1];

static int
exclude_feyn_diags_subp (int nsub, char fname[])
{
  int num;
  int stat;
  int nexcd = 0;
  int ndel, nrest, ncalc, ntot;
  long nmin, nmax;
  char line[32];
  FILE * ff = fopen (scat ("%s/tmp/%s", pathtouser, fname), "r");
  FILE * menu = fopen (MENUP_NAME, "r+b");

  ddel1[0] = -1;
  rd_menu (menu, 1, nsub, stc, &ndel, &ncalc, &nrest, &nmin);
  ntot = ndel + ncalc + nrest;
  nmax = nmin + ntot;

  if (ff) {
    fgets (line, 10, ff);
    if (!strcmp (line, "diagrams:")) {
      stat = fscanf (ff, "%d,", &num);
      while (1 == stat) {
        --num;
        if (nmin <= num && num < nmax) ++nexcd;
        stat = fscanf (ff, "%d,", &num);
      }
    } else
      fprintf (stderr, "%s\n", scat ("s_comphep.exe (warning): tmp/%s looks strange [1]. All diagrams are used.", fname));
  } else
    fprintf (stderr, "%s\n", scat ("s_comphep.exe (warning): No tmp/%s. All diagrams are used.", fname));

  if (nexcd > 0) {
    adiagram varc;
    FILE * diagrFile  = fopen (DIAGRP_NAME, "r+b");
    int diagSize   = sizeof (adiagram);
    int delMarkPos = (char *) &(varc.delMark) - (char *) &varc;

    fseek (ff, 0, SEEK_SET);
    fgets (line, 10, ff);
    stat = fscanf (ff, "%d,", &num);
    while (1 == stat) {
      --num;
      if (nmin <= num && num < nmax) {
        fseek (diagrFile, num * diagSize + delMarkPos, SEEK_SET);
        fwrite (ddel1, 1, 1, diagrFile);
        ndel++;
      }
      stat = fscanf (ff, "%d,", &num);
    }
    nrest = ntot - ndel - ncalc;
    wrt_menu (menu, 1, nsub, stc, ndel, ncalc, nrest, nmin);
    fclose (diagrFile);
  }
  if (ff) fclose (ff);
  fclose (menu);

  return 0;
}


static int
exclude_feyn_csdiags_subp (int nsub, char fname[])
{
  int num;
  int stat;
  int nexccsd = 0;
  int ndel, nrest, ncalc, ntot;
  long nmin, nmax;
  long fpos = 0;
  char line[32];
  FILE * ff = fopen (scat ("%s/tmp/%s", pathtouser, fname), "r");
  FILE * menu = fopen (MENUQ_NAME, "r+b");

  ddel1[0] = -1;
  rd_menu (menu, 2, nsub, stc, &ndel, &ncalc, &nrest, &nmin);
  ntot = ndel + ncalc + nrest;
  nmax = nmin + ntot;

  if (ff) {
    fgets (line, 10, ff);
    if (!strcmp (line, "diagrams:")) {
      stat = 1;
      while (1 == stat) {
        stat = fscanf (ff, "%d,", &num);
      }
      fgets (line, 18, ff);
      if (!strcmp (line, "squared diagrams:")) {
        fpos = ftell (ff);
        stat = fscanf (ff, "%d,", &num);
        while (1 == stat) {
          --num;
          if (nmin <= num && num < nmax) ++nexccsd;
          stat = fscanf (ff, "%d,", &num);
        }
      } else
        fprintf (stderr, "%s\n", scat ("s_comphep.exe (warning): tmp/%s looks strange [2]. All diagrams are used.", fname));
    } else
      fprintf (stderr, "%s\n", scat ("s_comphep.exe (warning): tmp/%s looks strange [3]. All diagrams are used.", fname));
  } else
    fprintf (stderr, "%s\n", scat ("s_comphep.exe (warning): No tmp/%s. All diagrams are used.", fname));

  if (nexccsd > 0) {
    csdiagram ars;
    FILE * diagrFile  = fopen (DIAGRQ_NAME, "r+b");
    int diagSize   = sizeof (csdiagram);
    int delMarkPos = (char *) &(ars.status) - (char *) &ars;

    fseek (ff, fpos, SEEK_SET);
    stat = fscanf (ff, "%d,", &num);
    while (1 == stat) {
      --num;
      if (nmin <= num && num < nmax) {
        fseek (diagrFile, num * diagSize + delMarkPos, SEEK_SET);
        fwrite (ddel1, 1, 1, diagrFile);
        ndel++;
      }
      stat = fscanf (ff, "%d,", &num);
    }
    nrest = ntot - ndel - ncalc;
    wrt_menu (menu, 2, nsub, stc, ndel, ncalc, nrest, nmin);
    fclose (diagrFile);
  }
  if (ff) fclose (ff);
  fclose (menu);

  return 0;
}


void
exclude_feyn_diags (char delname[])
{
  int i;

  for (i = 1; i <= subproc_f; ++i)
    exclude_feyn_diags_subp (i, delname);
}

void
exclude_feyn_csdiags (char delname[])
{
  int i;

  for (i = 1; i <= subproc_sq; ++i)
    exclude_feyn_csdiags_subp (i, delname);
}

void form_diag_list (FILE * f) {
  int num = 0;
  int nmax = 0;
  int ncur = 0;
  int new_str = 0;
  adiagram buff;
  FILE * diag = fopen (DIAGRP_NAME, "rb");

  while (FREAD1 (buff, diag))
    {
      if (-1 == buff.delMark)
        ++nmax;
    }
  fseek (diag, 0, SEEK_SET);

  fputs ("excluded diagrams:", f);
  while (FREAD1 (buff, diag))
    {
      ++num;
      if (-1 == buff.delMark) {
        ++ncur;
        ++new_str;
        fprintf (f, "%i,", num);
        if (nmax != ncur && 32 == new_str) {
          fputs ("\nexcluded diagrams:", f);
          new_str = 0;
        }
      }
    }
  fputs ("\n", f);
  fclose (diag);
}


void form_csdiag_list (FILE * f) 
{
  int num = 0;
  int nmax = 0;
  int ncur = 0;
  int new_str = 0;
  csdiagram buff;
  FILE * diag = fopen (DIAGRQ_NAME, "rb");

  while (FREAD1 (buff, diag))
    {
      if (-1 == buff.status)
        ++nmax;
    }
  fseek (diag, 0, SEEK_SET);

  fputs ("excluded squared diagrams:", f);
  while (FREAD1 (buff, diag))
    {
      ++num;
      if (-1 == buff.status) {
        ++ncur;
        ++new_str;
        fprintf (f, "%i,", num);
        if (nmax != ncur && 32 == new_str) {
          fputs ("\nexcluded squared diagrams:", f);
          new_str = 0;
        }
      }
    }
  fputs ("\n", f);
  fclose (diag);
}


void
prepare_process_dat (void)
{
  FILE *ff = fopen (scat ("%sprocess_gen.dat", pathtouser, f_slash), "w");

  fprintf (ff, "############################################################\n");
  fprintf (ff, "#         Data file for symb_script.pl\n");
  fprintf (ff, "#       For the symbolic batch script version 1.0\n");
  fprintf (ff, "############################################################\n\n");
  fprintf (ff, "#  You have to set the model number, which you are going to use.\n");
  fprintf (ff, "#  The model number corresponds to the string number of the model\n");
  fprintf (ff, "#  in the CompHEP model menu in the GUI mode. \n");
  fprintf (ff, "model number: %i\n\n", getModelNumberSymb());
  fprintf (ff, "#  Beam names can be taken from a table of beams \n");
  fprintf (ff, "#  (see CompHEP in the GUI regime). Energy unit is GeV\n");
  fprintf (ff, "beam 1: %s\n", beam[0].h.name);
  fprintf (ff, "beam 2: %s\n", beam[1].h.name);
  fprintf (ff, "beam energy 1: %f\n",   beam[0].energy);
  fprintf (ff, "beam energy 2: %f\n\n", beam[1].energy);
  fprintf (ff, "#  This string defines the final state of your process. Model \n");
  fprintf (ff, "#  particles and composite particles (see the corresponding table)\n");
  fprintf (ff, "#  can be used. \n");
  fprintf (ff, "final state:%s\n\n", getFinalstatech ());
  fprintf (ff, "#  If you'd like to exclude Feynman diagrams with some model\n");
  fprintf (ff, "#  particles (in propagators!), enter the particles here\n");
  fprintf (ff, "exclude diagrams with: %s\n", getExclprtlist ());
  fprintf (ff, "\n");
  fprintf (ff, "#  If you'd like to keep Feynman diagrams with some model\n");
  fprintf (ff, "#  particles (in propagators!), enter the particles here\n");
  fprintf (ff, "keep diagrams with: %s\n\n", getKeepprtlist ());
  fprintf (ff, "#  If you enter no, s_comphep generates diagrams and does not\n");
  fprintf (ff, "#  do symbolic calculations.\n");
  fprintf (ff, "make symbolic calculations(yes/no): yes\n\n");
  fprintf (ff, "#  If you enter no, comphep calculates all squared diagrams, \n");
  fprintf (ff, "#  but n_comphep will not be created.\n");
  fprintf (ff, "make n_comphep generator(yes/no): yes\n");
  fprintf (ff, "\n# excluded Feynman diagrams, simple list of digram numbers\n");
  form_diag_list (ff);
  fprintf (ff, "\n# excluded squared Feynman diagrams, simple list of digram numbers\n");
  form_csdiag_list (ff);
  fclose (ff);
}
