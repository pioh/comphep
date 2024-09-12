/*
* Copyright (C) 2000-2009, CompHEP Collaboration
* Author: V.A.Ilyin
* ------------------------------------------------------
*/
#include <unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#include "service2/include/syst.h"
#include "mix_cpyth1.h"

/* remove all temp event files, a-la {initial_file_name}___ */
static int removetmpfiles (int nf, char ** filename, int lenth) {
  int i;
  int err = 0;
  char * xtmpfle = malloc((lenth + 2) * sizeof(char));

  for(i = 0; i < nf; ++i) {
    strcpy(xtmpfle, filename[i]);
    strcat(xtmpfle, "___");
    err = unlink(xtmpfle);
    if (err) return err;
  }
  free (xtmpfle);
  return err;
}


/* transform events file of CompHEP format to the PEVLIB format */
static int TotEvPut (const char outFileName[], long posTotEv, long TotEv)
{
  int ilen;
  int nlen;
  char buff[128];
  char xfmt[128];
  char* xmix = "#Events_mixed_and_randomized";

  FILE * outFile = fopen (outFileName, "r+");
  fseek (outFile, posTotEv, SEEK_SET);
  fgets (buff, 1024, outFile);
  ilen = strlen(buff)-32;
  sprintf (buff, "%ld", TotEv);
  nlen = strlen (buff);
  if(ilen - nlen < 0)
  {
    fprintf(stdout, "Error: mixed events more than sum of input events.\n");
    return -1;
  }
  sprintf (xfmt, "%s%d%s", "%s = %-", ilen, "ld");
  sprintf (buff, xfmt, xmix, TotEv);
  fseek (outFile, posTotEv, SEEK_SET);
  fputs (buff, outFile); 
  return 0;
}


/*   transform events file of CompHEP format to the PEVLIB format */
static int toPEVLIB (char * ifile, char * ffile) {
  int version = -1; /* unknown version */
  int nweight;
  long nRec;
  double sigma;
  char buff[1024];
  char buff1[1024];
  char word[100];
  char xstr[72];
  char word_1[8];
  char word_2[7];
  char word_3[5];
  FILE * inFile;
  FILE * outFile;
  char * xdefis =  "#----------------------------------------------------------\n";
  char* xevents = "#Nproc ================== Events ==========================\n";
  char* xnsub   = "#Number_of_subprocesses";
  char* xcross  = "#Total_cross_section_(pb)";
  char* xmix    = "#Events_mixed_and_randomized";
  char* xpevmrk = "#PEVLIB_v.1.0 =============================================\n";

  inFile = fopen (ifile, "r");
  outFile = fopen (ffile, "w");

  strcpy (buff, xpevmrk);
  fputs (buff, outFile);

  fgets(buff1, 30, inFile);
  sscanf(buff1,"%s %s %s", word_1, word_2, word_3);
  fputs(buff1, outFile);

  if (!strcmp(word_1, "#CompHEP")) {
    if (!strcmp(word_3, "41.10")) version = 4110;
    if (!strcmp(word_3, "41.20")) version = 4120;
    if (!strcmp(word_3, "4.2p1")) version = 4120;
  } else {
    fprintf(stderr, "This file is not in CompHEP format!");
    return 1;
  }

// read-write of header
  if (4110 == version) {
    int i;
    for (i = 0; i < 9; ++i) {
      fgets(buff, 1024, inFile);
      sscanf(buff, "%s", word);
      if(!strcmp(word,"#Cross_section(Width)")) {
         sscanf(buff,"#Cross_section(Width) %lf", &sigma);
      }
      if (!strcmp(word,"#Number_of_events")) {
        sscanf(buff,"#Number_of_events  %li", &nRec);
      }
      fputs(buff,outFile);
    }
    strcpy (buff, xdefis);
    fputs (buff, outFile);
    strcpy (xstr, xnsub);
    strcat (xstr," = %d\n");
    sprintf (buff, xstr,1);
    fputs (buff, outFile);
    strcpy (xstr, xcross);
    strcat (xstr," = %E\n");
    sprintf (buff, xstr, sigma);
    fputs (buff, outFile);
    strcpy (xstr, xmix);
    strcat (xstr," = %ld\n");
    sprintf (buff, xstr, nRec);
    fputs (buff, outFile);
    strcpy (buff, xevents);
    fputs (buff, outFile);
    fgets (buff, 1024, inFile);
  } else {
    int while_exit = 0;
    while (while_exit == 0) {
      fgets (buff, 1024, inFile);
      sscanf (buff,"%s", word);
      if(!strcmp (word, "#Nproc")) {
        while_exit = 1;
      }
      fputs (buff, outFile);
    }
  }

// rewriting of events 
  while (1 == fscanf (inFile,"%d", &nweight)) {
    fgets (buff, 1024, inFile);
    while (nweight > 0) {
      if(4110 == version) {
        sprintf(buff1," 1  %s", buff);
      } else { 
        sprintf(buff1," 1%s", buff);
      }
      fputs (buff1, outFile);
      nweight--;
    }
  }

fclose(outFile);
return 0;
}


int mix_cpyth1 (const int nf, const char target[], const char names[], const int lenth)
{
  int j,k;
  int err;
  int nx;
  int nSubtot;
  int nlen;
  long i;
  long l,ipos;
  long nRectot;
  long maxnRec;
  long posTotEv;
  long **pos;
  long *fpos;
  long *pos0;
  long *nsub;
  long *nRec;
  long *nLeft;
  double sigmatot;
  double xrn;
  double *sigma;
  double *ri;
  char buff[1024];
  char buff1[1024];
  char word[128];
  char xstr[128];
  char *fmap;
  char **map;
  char **filename;

  FILE *outFile;
  FILE **infile;

// initialize global variables
  char* xevents  = "#Nproc ================== Events ==========================\n";
  char* xnsub    = "#Number_of_subprocesses";
  char* xcross   = "#Total_cross_section_(pb)";
  char* xmix     = "#Events_mixed_and_randomized";
  char* xpevlib  = "#PEVLIB_v.1.0";
  char* xpevlib2 = "#PEVLIB_v.1.0 =============================================\n";

// allocate the memory space for work dimensions
  pos      = malloc(nf * sizeof (long*));
  map      = malloc(nf * sizeof (char*));
  sigma    = malloc(nf * sizeof (double));
  ri       = malloc(nf * sizeof (double));
  pos0     = malloc(nf * sizeof (long*));
  nsub     = malloc(nf * sizeof (long*));
  nRec     = malloc(nf * sizeof (long*));
  nLeft    = malloc(nf * sizeof (long*));
  infile   = malloc(nf * sizeof (FILE*));
  filename = malloc(nf * sizeof (char*));

  for (i = 0; i < nf; ++i) {
    filename[i] = malloc ((lenth + 1)* sizeof(char));
    strncpy (filename[i], names + i * lenth, lenth);filename[i][lenth] = 0;
    trim(filename[i]);
  }

/****************************************************************/
/* open output files for writing                                */
  outFile = fopen(target, "w");

/****************************************************************/
/* check arguments - input event files */
  fprintf(stdout, "mix (info): files to mix and randomize:\n");
  for(i = 0; i < nf; ++i) {
    infile[i] = fopen (filename[i], "r");
    if (!infile[i]) {
      fprintf (stdout,  "mix (warning): file %s not found!\n", filename[i]);
      exit (2);
    } else {
      fprintf (stdout, "mix (info): %s\n", filename[i]);
    }
  }

/****************************************************************/
/*        check PEVLIB headers of event files                   */
  for (i = 0; i < nf; ++i) {
    fgets(buff, 1024, infile[i]);
    sscanf(buff,"%s", word);
    if(!strcmp(word,"#CompHEP")) {
      char * tmpname = malloc ((lenth + 4)* sizeof(char));
      fclose (infile[i]);
      sprintf (tmpname, "%s___", filename[i]);
      toPEVLIB (filename[i], tmpname);
      infile[i] = fopen (tmpname, "r");
      free(tmpname);
    }
  }

/****************************************************************/
/*        read-write headers of event files                     */
  sigmatot = 0;
  nSubtot = 0;
  strcpy (buff, xpevlib2);
  fputs (buff, outFile);
  for (i = 0; i < nf; ++i) {
    while(1) {
      fgets(buff, 1024, infile[i]);
      sscanf(buff,"%s",word);
      if(0 == strcmp(word,xnsub)) break;
      if(strcmp(word, xpevlib)) fputs(buff, outFile);
    }
    strcpy (xstr, xnsub);
    strcat (xstr," = %ld");
    sscanf (buff, xstr, nsub + i);
    fgets (buff, 1024, infile[i]);
    strcpy (xstr, xcross);
    strcat (xstr," = %lf");
    sscanf (buff, xstr, &sigma[i]);
    fgets (buff, 1024, infile[i]);
    strcpy (xstr, xmix);
    strcat (xstr," = %ld");
    sscanf (buff, xstr, &nRec[i]);
    fgets (buff, 1024, infile[i]);
    pos0[i] = ftell(infile[i]);   /* position of the 1st event beginning */
    nSubtot = nSubtot+nsub[i];
    nsub[i] = nSubtot-nsub[i];
    sigmatot = sigmatot+sigma[i];
  }
   
/****************************************************************/
/*   allocate the memory space for dimensions                   */
  nRectot = 0;
  maxnRec = 0;
  for (i = 0; i < nf;++i) {
    j = 0;
    while (1) {
      if(0 == fgets(buff,1024,infile[i])) break;
      j++;
    }
    nRec[i] = j;
    nLeft[i] = j;
    nRectot = nRectot + nRec[i];
    if(nRec[i]>maxnRec) maxnRec = nRec[i];
  }

  for (i = 0; i < nf; ++i) {
    fpos = malloc(maxnRec * sizeof(long));
    fmap = malloc(maxnRec * sizeof(char));
    fseek(infile[i], pos0[i], SEEK_SET);
    j = 0;
    fprintf (stdout, "mix (info): reading file %s:\n", filename[i]);
    while(1) {
      fpos[j] = ftell(infile[i]);
      fmap[j] = 0;
      if(!fgets(buff, 1024, infile[i])) break;
      ++j;
      if (0 == j % 20000) fprintf (stdout, "mix (info): %i events read\n", j);
    }
    fprintf (stdout, "mix (info): %i events read\n", j);
    pos[i]=fpos;
    map[i]=fmap;
  }

/****************************************************************/
/*     write general info in header                             */
  sprintf(buff,"%s = %d\n", xnsub, nSubtot);
  fputs(buff,outFile);
  sprintf(buff,"%s = %E\n", xcross, sigmatot);
  fputs(buff,outFile);
  sprintf(buff,"%s = %ld\n", xmix, nRectot);
  posTotEv=ftell(outFile);
  fputs(buff,outFile);
  strcpy(buff,xevents);
  fputs(buff,outFile);

/******************************************************************/
/* evaluate the map of subprocesses weights in the interval [0,1] */
   xrn = 0;
   for (i = 0; i < nf; ++i) {
     ri[i] = xrn + sigma[i]/sigmatot;
     xrn = ri[i];
   }

/****************************************************************/
/*         mix events from the files                            */
  for(i = 0; i < nRectot; ++i) {
    xrn = drand48();
    k = 0;
    while (xrn > ri[k]) ++k;
    if (0 == nLeft[k]) break;

    l = drand48() * nRec[k];
    while(map[k][l]) {
      l++;
      if (l == nRec[k]) l = 0;
    }

    fseek(infile[k], pos[k][l], SEEK_SET);
    map[k][l] = 1;
    ipos = ftell(infile[k]);
    fgets(buff,1024,infile[k]);
    sscanf(buff, "%d", &nx);
    sscanf(buff, "%s", word);
    nlen = strlen(word);
    fseek(infile[k], ipos, SEEK_SET);
    fgets(word, nlen + 6, infile[k]);
    fgets(buff1, 1024, infile[k]);
    sprintf(buff, " %-5li%s",nsub[k] + nx, buff1);
    fputs(buff, outFile);
    --(nLeft[k]);
  }
  fclose(outFile);

  fprintf (stdout, "\nmix (info): final statistics:\n");
  for (i = 0; i < nf; ++i) {
    fprintf (stdout,  "mix (info): file %s: %li events used from %li events", filename[i], nRec[i] - nLeft[i], nRec[i]);
    if (0 == nLeft[i])  fprintf (stdout, ", file exhausted!\n");
    else fputs ("\n", stdout);
  }
  fprintf (stdout, "mix (info): %li mixed/randomized events written to %s\n", nRectot, target);

  err = TotEvPut (target, posTotEv, nRectot);
  if (err) return err;
  err = removetmpfiles(nf, filename, lenth);
  if (err) return err;

  free (pos);
  free (map);
  free (sigma);
  free (ri);
  free (pos0);
  free (nsub);
  free (nRec);
  free (nLeft);
  free (infile);
  free (fpos);
  free (fmap);

  return 0;
}
