/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "service2/include/syst.h"

#include "tag_reader.h"
#include "tag_writer.h"
#include "tag_parser.h"
#include "fill.h"
#include "validity.h"
#include "structures.h"
#include "compare.h"
#include "mix_cpyth2.h"

int mix_cpyth2 (const int nf, const char rtarget[], const char names[], const int lenth)
{
  int t;
  int nSubtot;
  int nProc;
  int *shift;
  int totwrtevents;
  int distilling_used = 0;
  long i,l;
  long nRectot;
  long maxnRec;
  long num;
  long **pos;
  long fpos;
  long *pos0;
  long *nRec;
  long *real_nRec;
  long *nLeft;
  double *ri, xrn;
  double sigmatot = 0.0;
  char buff[2048];
  char buff1[1024];
  char * etarget;
  char **map;
  char **filename;

  FILE** infile;
  FILE* outFile;
  FILE* evnFile;

  process_ theProcessUP;
  tags **tbase;
  proc_pos *proc_position;

  etarget = malloc ((strlen (rtarget) + 10 ) * sizeof (char));
  strcpy (etarget, rtarget);
  strcat (etarget, ".events");

/* allocate the memory space for work dimensions */
  tbase    = malloc (nf * sizeof (tags *));
  pos0     = malloc (nf * sizeof (long));
  shift    = malloc (nf * sizeof (int));
  infile   = malloc (nf * sizeof(FILE*));
  filename = malloc (nf * sizeof(char*));

  for (i = 0; i < nf; ++i) {
    filename[i] = malloc ((lenth + 1) * sizeof(char));
    strncpy (filename[i], names + i * lenth, lenth);
    filename[i][lenth] = 0;
    trim (filename[i]);
  }

/* check arguments - input event files */
/*  fprintf (stdout, "mix (info): files to mix and randomize:\n");*/
  for (i = 0; i < nf; ++i) {
    infile[i] = fopen (filename[i], "r");
    if (!infile[i]) {
      fprintf (stdout, "mix (error): file %s not found. Exit\n", filename[i]);
      return 2;
    } else {
/*      fprintf (stdout, "%s\n", filename[i]);*/
    }
  }

/* read out headers of event files and check header structure */
  for (i = 0; i < nf; ++i) {
    tbase[i] = init_cap(100);
    cup_reader (infile[i], tbase[i]);
    syntax_validity (tbase[i]);
    structure_validity (tbase[i]);
    physics_validity (tbase[i]);
    pos0[i] = ftell (infile[i]) + 1;  /* position of the 1st event beginning */
  }

/* fill HLE structures and compare subprocesses */
/* BE AWARE!!!! THE NUMBER OF SUBPROCESSES ARE LIMITED IN LHEstructures.h!!!!! */
  fill_LH_structures (nf, tbase, &theProcessUP, names, lenth);
  compare_extra_info (nf, tbase, &theProcessUP);
  nSubtot = theProcessUP.proc_info.NprocRUP;
  proc_position = malloc (nSubtot * sizeof (proc_pos));
  compare_processes (nf, tbase, proc_position);
  fprintf (stdout, "mix (info): the total number of subprocesses in all files is %i\n", nSubtot);

/* allocate memory for dimensions */
  nRec = malloc (nSubtot * sizeof (long));
  real_nRec = malloc (nSubtot * sizeof (long));
  nLeft = malloc (nSubtot * sizeof (long));
  ri = malloc (nSubtot * sizeof (double));
  pos = malloc (nSubtot * sizeof (long *));
  map = malloc (nSubtot * sizeof (char *));

  nRectot = 0;
  maxnRec = 0;
  for (i = 0; i < nSubtot; i++)
    real_nRec[i] = 0;

/* test event files and detect the number of events in the files */
  t = 0;
  for (i = 0; i < nf; i++) {
    int nSubproc = 0;
    int nEvCount = 0;
    int tottagnum = get_tag (0, tbase[i], "total");
    if (tottagnum != -1) {
      nSubproc = get_ival (0, "Nproc", tbase[i]->tag[tottagnum]);
    } else {
      fprintf (stderr, "mix (error): the file %s does not contain the total tag\n", filename[i]);
      exit (2);
    }

    fseek (infile[i], pos0[i], SEEK_SET);
    fprintf (stdout, "mix (info): read file %s:\n", filename[i]);
    while (1) {
      if (!fgets (buff, 1000, infile[i]))
        break;
      sscanf (buff, "%i", &nProc);
      num = t + nProc - 1;
      if (num < nSubtot) {
        ++real_nRec[num];
      } else {
        fprintf (stderr, "mix (error): strange event in %s with subprocess number (%li)\n", filename[i], num);
        fprintf (stderr, "             greater that the total number of subprocesses (%i)\n", nSubtot);
        exit (2);
      }
      if (nEvCount > 0 && 0 == nEvCount % 10000) {
        fprintf (stdout, "mix (info): %i events read\n", nEvCount);
      }
      ++nEvCount;
    }
    fprintf (stdout, "mix (info): %i events read\n", nEvCount);
    t += nSubproc;
  }

/* check the number of events in headers and in files. calculate the total CS */
  sigmatot = 0.0;
  for (i = 0; i < nSubtot; i++) {
    nRec[i] = proc_position[i].Nevents;
    if (real_nRec[i] != nRec[i]) {
      fprintf (stdout, "mix (warning): real number of events for subprocess %li (%li) in file %s does not coincide with header info (%li)\n",
               i + 1, real_nRec[i], filename[proc_position[i].nfile], nRec[i]);
      nRec[i] = real_nRec[i];
    }
    nLeft[i] = nRec[i];
    nRectot = nRectot + nRec[i];
    if (nRec[i] > maxnRec) {
      maxnRec = nRec[i];
    }
    pos[i] = malloc (nRec[i] * sizeof (long));
    map[i] = malloc (nRec[i] * sizeof (char));
    sigmatot += theProcessUP.proc_info.crossecUP[i];
  }

/* build maps of events in files */
  t = 0;
  for (i = 0; i < nSubtot; i++) {
    real_nRec[i] = 0;
  }

  for (i = 0; i < nf; i++) {
    int tottagnum = get_tag (0, tbase[i], "total");
    if (tottagnum != -1) {
      tottagnum = get_ival (0, "Nproc", tbase[i]->tag[tottagnum]);
    } else {
      fprintf (stderr, "mix (error): the file %s does not contain the total tag.\n", filename[i]);
      exit (2);
    }
    shift[i] = t;

    fseek (infile[i], pos0[i], SEEK_SET);
    while (1) {
      fpos = ftell (infile[i]);
      if (!fgets (buff, 1000, infile[i]))
        break;
      sscanf (buff, "%i", &nProc);
      num = t + nProc - 1;
      pos[num][real_nRec[num]] = fpos;
      map[num][real_nRec[num]] = 0;
      real_nRec[num]++;
    }
    t += tottagnum;
  }

/* evaluate map of subprocess weights in the interval [0,1] */
  xrn = 0;
  for (i = 0; i < nSubtot; i++) {
    ri[i] = xrn + theProcessUP.proc_info.crossecUP[i] / sigmatot;
    xrn = ri[i];
  }

/* mix events */
  evnFile = fopen (etarget, "w");
  outFile = fopen (rtarget, "w");
/*
  {
    FILE * tmp = fopen ("test1.log", "w");
    cap_writer (tmp, *tbase);
  }
*/
  write_cap_new (outFile, tbase, nf, theProcessUP);
  fclose (outFile);

  totwrtevents = 0;
  for (i = 0; i < nRectot; i++) {
    int nx;
    int jf = 0;
    xrn = drand48 ();
    while (xrn > ri[jf]) {
      jf++;
    }

    if (nLeft[jf] == 0) break;
    l = 0;
    while (map[jf][l]) {
      l++;
      if (l == nRec[jf]) {
        l = 0;
      }
    }
    fseek (infile[proc_position[jf].nfile], pos[jf][l], SEEK_SET);
    map[jf][l] = 1;
    fgets (buff, 1000, infile[proc_position[jf].nfile]);
    sscanf (buff, "%i:%[^\n]", &nx, buff1);
    fprintf (evnFile, " %i:%s\n", (int) (shift[proc_position[jf].nfile] + nx), buff1);
    --(nLeft[jf]);
    ++totwrtevents;
  }
  fclose (evnFile);

/* print out some statistics after mixing */
  fprintf (stdout, "mix (info): final statistics:\n");
/*
  for (i = 0; i < nf; ++i) {
    fprintf (stdout,  "mix (info): file %s: %li events used from %li events", filename[i], nRec[i] - nLeft[i], nRec[i]);
    if (0 == nLeft[i])  fprintf (stdout, ", file exhausted!\n");
    else fputs ("\n", stdout);
  }
*/
  for (i = 0; i < nSubtot; i++) {
    fprintf (stdout,  "mix (info): file %s: %li events used from %li in subprocess %li", filename[proc_position[i].nfile], nRec[i] - nLeft[i], nRec[i], i);
    if (0 == nLeft[i])  fprintf (stdout, ", file exhausted!\n");
    else fputs ("\n", stdout);
    if (nRec[i] - nLeft[i] == 0) ++distilling_used;
  }
  if (1 < distilling_used) {
    fprintf (stdout,  "mix (warning): distilling needed!\n");
  }

  fprintf (stdout, "mix (info): %i mixed/randomized events written to %s\n", totwrtevents, rtarget);

/* form output event file */
  for (i = 0; i < nSubtot; i++) {
    nLeft[i] = nRec[i] - nLeft[i];
  }
  final_write_cap (rtarget, nLeft, totwrtevents);

  outFile = fopen (rtarget, "a");
  evnFile = fopen (etarget, "r");
  while (!feof (evnFile)) {
    fgets (buff, 2048, evnFile);
    fputs (buff, outFile);
  }

  for (i = 0; i < nSubtot; i++) {
    nLeft[i] += 0;
  }

  free (pos);
  free (map);
  free (proc_position);
  free (tbase);
  free (pos0);
  free (shift);
  free (nRec);
  free (nLeft);
  free (ri);
  for (i = 0; i < nf; ++i) {
    free (filename[i]);
  }
  free (filename);
  remove (etarget);

  if (1 < distilling_used) {
    distilling (rtarget);
  }

  return 1;
}
