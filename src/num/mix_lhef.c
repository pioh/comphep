/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef LIBXML
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"

#include "LesHouches.h"
#include "lhef_routines.h"
#include "event_reader.h"
#include "mix_lhef.h"

int mix_lhef (int nf, const char target[], const char names[], const int lenth, int zrandom_used)
{

  int i, j;
  int mf;
  int final;
  int nRecotot;
  int nLefttot;
  int *nLeft;
  int *nUsed;
  int *file;
  int *xmlerr;
  int *nShift;
  int *idProcShift;

  long endpos;
  long *pos0;
  long *wpos;
  long *pos;

  double *ri;
  double delta;
  double *CSec;
  double *CSerr;
  double xrn;
  double errtot;
  double sigmatot;
  char *map;
  char **filename;

  FILE** infile;
  FILE* outFile;

  filename = malloc (nf * sizeof(char*));
  fprintf (stdout, "mix (info): files to mix and randomize:\n");
  for (i = 0; i < nf; ++i) {
    FILE * f;
    filename[i] = malloc ((lenth + 2) * sizeof(char));
    strncpy (filename[i], names + i * lenth, lenth); filename[i][lenth] = 0;
    trim (filename[i]);
    f = fopen (filename[i], "r");
    if (!f) {
      fprintf (stdout, "mix (warning): file not found: %s\n", filename[i]);
      return 2;
    } else {
      fprintf (stdout, "%s\n", filename[i]);
    }
    fclose (f);
  }

  xmlerr = malloc (nf * sizeof (int));
  pos0   = malloc (nf * sizeof (long));
  CSec   = malloc (nf * sizeof (double));
  CSerr  = malloc (nf * sizeof (double));
  nLeft  = malloc (nf * sizeof (int));
  nUsed  = malloc (nf * sizeof (int));

  idProcShift = malloc (nf * sizeof (int));
  idProcShift[0] = 0;
  for (i = 0; i < nf; ++i) {
#ifdef LIBXML
    xmlerr[i] = formXMLtree (filename[i], i);
#else
    xmlerr[i] = analyzeLHEfile (filename[i], i);
#endif
    pos0[i]   = getEventPosition (i);
    nLeft[i]  = getEventNumber (i);
    CSec[i]   = getCrossSection (i);
    CSerr[i]  = getCrossSectionErr (i);
/*    if (0 < i) idProcShift[i] = idProcShift[i] + getTotProcNumber ();*/
    if (0 < i) idProcShift[i] = idProcShift[i] + 1;
  }

  mf = nf;
  for (i = nf - 1; i >= 0; --i) {
    if (0 > xmlerr[i] || 0 > pos0[i] || 0 > nLeft[i] || 0 > CSec[i]) {
      fprintf (stderr, "mix (error): Remove file %s\n", filename[i]);
      --mf;
      for (j = i; i < mf; ++i) {
        pos0[j]  = pos0 [j + 1];
        nLeft[j] = nLeft [j + 1];
        CSec[j]  = CSec [j + 1];
        strcpy (filename[j], filename[j + 1]);
      }
    }
  }
  nf = mf;

  for (i = 0; i < nf; ++i) {
    nUsed[i] = 0;
  }

  nLefttot = 0;
  sigmatot = 0.0;
  errtot = 0.0;
  for (i = 0; i < nf; ++i) {
    nLefttot += nLeft[i];
    sigmatot += CSec[i];
    delta = 0.0;
    if (0.0 < CSec[i]) delta = CSerr[i]/CSec[i];
    errtot += delta * delta;
  }
  errtot = sigmatot * sqrt (errtot);

  nShift = malloc ((nf + 1) * sizeof (int));
  pos    = malloc ((nLefttot + 1) * sizeof (long));
  map    = malloc ((nLefttot + 1) * sizeof (char));
  for (i = 0; i < nf; ++i) {
    int num = 0;
    FILE * f = fopen (filename[i], "r");

    nShift[0] = 0;
    pos[nShift[i]] = pos0[i];
    map[nShift[i]] = 0;
    fprintf (stdout, "mix (info): reading file %s:\n", filename[i]);
    while (!testLHAevent (f, filename[i], pos[num + nShift[i]]) && num < nLeft[i]) {
      ++num;
      if (0 == num % 10000) fprintf (stdout, "mix (info): %i events read\n", num);
      pos[num + nShift[i]] = ftell (f);
      map[num + nShift[i]] = 0;
    }
    fprintf (stdout, "mix (info): %i events read in total\n", num);
    if (num != nLeft[i]) {
      fprintf (stderr, "mix (warning): the real number of events (%i) in file %s does not correspond reporting number (%i) in header.\n", num, filename[i], nLeft[i]);
      nLeft[i] = num;
    }
    nShift[i + 1] = nLeft[i] + nShift[i];
    fclose (f);
  }

  xrn = 0.0;
  ri = malloc (nf * sizeof (double));
  for (i = 0; i < nf; i++) {
    ri[i] = xrn + CSec[i] / sigmatot;
    xrn = ri[i];
  }

  wpos = malloc (nLefttot * sizeof (long));
  file = malloc (nLefttot * sizeof (int));
  nRecotot = 0;
  final = 0;
  while (nRecotot < nLefttot) {
    xrn = drand48 ();
    i = 0;
    while (xrn > ri[i]) {
      i++;
    }
    file[nRecotot] = i;

    if (0 >= nLeft[i] ) {
      final = i;
      break;
    }

    j = 0;
    while (map[nShift[i] + j]) {
      ++j;
    }
    map[nShift[i] + j] = 1;
    wpos[nRecotot] = pos[nShift[i] + j];
    ++nRecotot;
    --(nLeft[i]);
    ++(nUsed[i]);
    if (0 == nRecotot % 10000) fprintf (stdout, "mix (info): %i events processed\n", nRecotot);
  }
  fprintf (stdout, "mix (info): %i events processed in total\n", nRecotot);

  set_cs (sigmatot, errtot);
  outFile = fopen (target, "w");
  if (!outFile) {
    return 0;
  }

  fprintf (outFile, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (outFile, "<header>\n");


  fprintf (outFile, "<!-- File generated with CompHEP %s -->\n", getCHEPversion ());
#ifdef LIBXML
  fprintf (outFile, "<!-- \n"
                    "     This file is compatible with the Les Houches event file\n"
                    "     format (hep-ph/0609017), but contains extra HepML tags.\n"
                    "-->\n");
  fprintf (outFile, "%s", prepare_hepml_header_libxml2_dynamic ());
#else
  fprintf (outFile, "<!-- \n"
                    "     This file is compatible with the Les Houches event file\n"
                    "-->\n");
#endif

  fprintf (outFile, "</header>\n");
  fprintf (outFile, "<init>\n");
  fprintf (outFile, "%i %i %17.10E %17.10E %i %i %i %i 3 %i\n",
                    getPbeam (0), 
                    getPbeam (1), 
                    getEbeam (0), 
                    getEbeam (1), 
                    getPDFLIBgroup (0), 
                    getPDFLIBgroup (1), 
                    getPDFLIBset (0), 
                    getPDFLIBset (1),
                    nf);
  for (i = 0; i < nf; ++i) {
    fprintf (outFile, "%17.10E %17.10E %17.10E %i\n", CSec[i], CSerr[i], 1.0, i + 1);
  }
  fprintf (outFile, "</init>\n");
  endpos = ftell (outFile);

  infile   = malloc (nf * sizeof (FILE*));
  for (i = 0; i < nf; ++i) {
    infile[i] = fopen (filename[i], "r");
  }
  int num_z_changed = 0;
  for (i = 0; i < nRecotot; ++i) {
    int change_z_axis = 0;
    if (zrandom_used)
      if (0.5 < drand48 ()) {
        change_z_axis = 1;
        ++num_z_changed;
      }
    if (0 == i % 10000) fprintf (stdout, "mix (info): %i events written\n", i + 1);
    moveLHAevent (file[i] + 1, filename[file[i]], infile[file[i]], wpos[i], outFile, endpos, change_z_axis);
    endpos = ftell (outFile);
  }
  fprintf (outFile, "</LesHouchesEvents>\n");
  fclose (outFile);

  fprintf (stdout, "\nmix (info): final statistics:\n");
  for (i = 0; i < nf; ++i) {
    fprintf (stdout,  "mix (info): file %s: %i events used from %i events", filename[i], nUsed[i], nUsed[i] + nLeft[i]);
    if (0 == nLeft[i])  fprintf (stdout, ", file exhausted!\n");
    else fputs ("\n", stdout);
  }
  fprintf (stdout,  "mix (info): %i mixed/randomized events written to %s\n", nRecotot, target);
  if (zrandom_used) fprintf (stdout,  "mix (info): %i events have reversed Z-axis\n", num_z_changed);

  {
    int filesize = 0;
    char check_sum[128];
    struct stat stt;
    char command[1024];

    stat (target, &stt);
    filesize = stt.st_size;
    for (i = 0; i < 128; ++i) check_sum[i] = 0;
    sprintf (command, "md5sum %s > md5sum.log", target);
    system(command);
    FILE * md5 = fopen ("md5sum.log", "r");
    if (md5) {
      int err = fscanf (md5, "%s", check_sum);
      if (1 != err) return -1;
    }
    unlink ("md5sum.log");
    set_final_numbers (nRecotot, filesize, sigmatot, errtot, check_sum);
#ifdef LIBXML
    write_file_header_libxml2 (target, "r+");
#endif
  }

  for (i = 0; i < nf; ++i) {
    free (filename[i]);
  }
  free (filename);

  free (xmlerr);
  free (pos0);
  free (CSec);
  free (CSerr);
  free (nLeft);
  free (nUsed);
  free (nShift);
  free (idProcShift);
  free (pos);
  free (map);
  free (infile);
  free (ri);
  free (wpos);
  free (file);

  return 1;
}

