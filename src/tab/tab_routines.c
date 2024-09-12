/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <stdio.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/paragraphs.h"
#include "service2/include/4_vector.h"

#include "num/include/LesHouches.h"
#include "num/include/tag_reader.h"
#include "num/include/tag_parser.h"
#include "num/include/tag_routines.h"
#include "num/include/event_reader.h"
#include "num/include/lhef_routines.h"
#include "num/include/phys_val.h"
#include "num/include/evnt_tools.h"

#include "e_tools.h"

#define const
#include "out_ext.h"
#undef const

static double  *hist;
static double *dhist;

static int
keep_histo (double *hist, double *dhist, char key, char *plist, double minX, double maxX, int nbin)
{
  int i;
  char xname[200];
  char yname[200];
  char xunits[100];
  double dx;

  for (i = 0; i < nin_ + nout_; i++)
    {
      if (i == nin_)
        fprintf (stdout, " ->");
      else if (i)
        fprintf (stdout, ",");
      pinf_ (1, i + 1, xname, NULL);
      fprintf (stdout, " %s", xname);
    }

  if (nin_ == 2)
    {
      sprintf (yname, "Diff. cross section [pb");
    }
  else
    {
      sprintf (yname, "Diff. width [GeV");
    }
  fprintf (stdout, "\n");

  xName (key, plist, xname, xunits);
  fprintf (stdout, "x-axis: \"%s\"  from %f to %f N_bins= %d\n", xname, minX, maxX, nbin);
  fprintf (stdout, "%s", yname);

  if (xunits[0])
    fprintf (stdout, "/%s", xunits);
  fprintf (stdout, "]");

  /* print histo out */
  dx = (maxX - minX) / nbin;
  for (i = 0; i < nbin; ++i)
    {
      fprintf (stdout, "\n%-12E %-12E", minX + dx * i + dx / 2., hist[i]);
      if (dhist)
        fprintf (stdout, " +/-  %-12E", dhist[i]);
    }
  fprintf (stdout, "\n");

  return 1;
}


static double 
tab_format2 (char fname[], char thevar[], int nbin, double minX, double maxX)
{
  /* New XML-style fromat */
  int i, j, k;
  int n;
  int format;
  int npr;
  int in;
  int out;
  int Nproc;
  int *Npartons;
  int **pnum;
  int nEvents = 0;
  int minnin_, minnout_;

  char key;
  char plist[MAXNP];
  char restfr[10]={0};
  char buff[STRSIZ];
  char buf[STRSIZ];
  midstr name_proc;

  double wght;
  double **mass;

  string_comnd com;
  tags * head = init_cap (1);
  FILE * f = fopen (fname, "r");

  cup_reader (f, head);
  strcpy (com.name, "Nproc");
  get_tag_with1com (0, head, "total", &com);
  Nproc = atoi (com.value);
  Npartons = malloc (Nproc * sizeof (int));
  pnum     = malloc (Nproc * sizeof (int*));
  mass     = malloc (Nproc * sizeof (double*));

  /* check file format and receive nin_ and nout_ */
  if (check_cpyth2 (head)) {
    return -1.0;
  }

  /* Get number of events and total CS from the total tag */
  n = get_tag (0, head, "total");
  wght = get_fval (0, "CrosSec", head->tag[n]);

  minnin_ = 1000;
  minnout_ = 1000;
  for (i = 0; i < Nproc; i++) {
    pnum[i] = malloc (128 * sizeof (int));
    sprintf (com.value, "%i", i + 1);

  /* Get number of partons, parton's names and nin_/nout_ from the process tag */
    strcpy (com.name, "ID");
    n = get_tag_with_exactcom (0, head, "process", com);
    Npartons[i] = get_ival (0, "Nparton", head->tag[n]);
    strcpy (name_proc, get_cval (0, "name", head->tag[n]));
    getNames (name_proc);
    if (minnin_ > nin_) minnin_ = nin_;
    if (minnout_ > nout_) minnout_ = nout_;

  /* Construction of event format */
    strcpy (com.name, "IDprocess");
    sprintf (com.value, "%i", i + 1);
    format = get_tag_with_exactcom (0, head, "format", com);
    strcpy (com.name, "p1.3");
    pnum[i][3] = tag_contain_com (&com, head->tag[format]) - 2;
    strcpy (com.name, "p2.3");
    pnum[i][7] = tag_contain_com (&com, head->tag[format]) - 2;
    for (j = 2; j < Npartons[i]; ++j) {
      for (k = 1; k < 4; ++k) {
        sprintf (com.name, "p%i.%i", j + 1, k);
        n = tag_contain_com (&com, head->tag[format]);
        pnum[i][4 * j + k] = n - 2;
      }
    }

  /* Get partons masses from parton tags */
    strcpy (com.name, "IDprocess");
    sprintf (com.value, "%i", i + 1);
    n = -1;
    mass[i] = malloc (Npartons[i] * sizeof (double));
    for (j = 0; j < Npartons[i]; j++)
      {
        n   = get_tag_with_exactcom (n + 1, head, "parton", com);
        in  = get_ival (0, "in", head->tag[n]);
        out = get_ival (0, "out", head->tag[n]);
        if (in)
          mass[i][in - 1] = get_fval (0, "mass", head->tag[n]);
        if (out)
          mass[i][out - 1] = get_fval (0, "mass", head->tag[n]);
      }
  }
  nin_ = minnin_;
  nout_ = minnout_;

   checkPhysVal (thevar, &key, plist);

/* event manipulation */
  fgets (buf, 1000, f);
  while (fgets (buf, 1000, f)) {
    ++nEvents;
    sscanf (buf, "%i:%[^\n]", &npr, buff);
    npr--;
    for (i = 0; i < 4 * nin_; i++) pvect[i] = 0.0;
    i = 0;
    char * pch = strtok (buff, ":");
    while (pch != NULL) {
      if (i == pnum[npr][3]) pvect[3] = atof (pch);
      if (i == pnum[npr][7]) pvect[7] = atof (pch);
      for (j = 2; j < Npartons[npr]; ++j) {
        for (k = 1; k < 4; ++k)
          if (i == pnum[npr][4 * j + k]) pvect[4 * j + k] = atof (pch);
      }
      pch = strtok (NULL, ":");
      ++i;
    }
    for (i = 0; i < nin_ + nout_; i++)
      pvect[4 * i] = ENERGY (mass[npr][i], pvect + 4 * i + 1);

    i = nbin * (calcPhysVal (key, plist, restfr) - minX) / (maxX - minX);
    if (0 <= i && i < nbin) {
      ++hist[i];
    }
  }

  free (Npartons);
  free (pnum);
  free (mass);

  return wght / nEvents;
}

static double 
tab_format1 (char fname[], char thevar[], int nbin, double minX, double maxX)
{
  /* Old CompHEP_41.10 fromat, but sample contains several subprocesses */
  int i, j;
  int nEvents = 0;
  int stop = 0;
  int Nproc = -1;
  int npr;
  char key;
  int minnin_, minnout_;
  char restfr[10]={0};
  char plist[MAXNP];
  char bufff[MAXINOUT + 3][STRSIZ];
  char buff[STRSIZ];
  double **mass;
  FILE * f = fopen (fname, "r");

  rw_paragraph rd_array0[1] = {
    {"PEVLIB_v.1.0", skipHeadLine}
  };
  rw_paragraph rd_array1[7] = {
    {"CompHEP", NULL},
    {"PROCESS", old_getNames},
    {"Initial_state", skipLine},
    {"MASSES", getMasses},
    {"Cross_section(Width)", skipLine},
    {"Number_of_events", skipLine},
    {"----------------------------------------------------------", skipHeadLine}
  };
  rw_paragraph rd_array2[4] = {
    {"Number_of_subprocesses", skipLine},
    {"Total_cross_section_(pb)", readTotCS},
    {"Events_mixed_and_randomized", skipLine},
    {"Nproc", skipHeadLine}
  };

  /* Get the number of subprocesses in the file */
  while ((NULL == fgets (buff, 1024, f)) || !stop) {
    stop = sscanf (buff, "#Number_of_subprocesses = %d", &Nproc);
  }
  if (Nproc < 1) {
    fprintf (stderr, "mk_tab (error): unknown number of subprocesses! Nproc = %i\n", Nproc);
    return -1.0;
  }

  /* Return to file begin and read first tag: PEVLIB version */
  rewind (f);
  readParagraphs (f, 1, rd_array0);

  /* Read subprocess info: particle names and masses */
  minnin_ = 1000;
  minnout_ = 1000;
  mass  = malloc (Nproc * sizeof (double *));
  for (j = 0; j < Nproc; j++) {
    readParagraphs (f, 7, rd_array1);
    if (minnin_ > nin_) minnin_ = nin_;
    if (minnout_ > nout_) minnout_ = nout_;
    mass[j] = malloc ((nin_ + nout_) * sizeof (double));
    for (i = 0; i < nin_ + nout_; i++) {
      pinf_ (1, i + 1, NULL, mass[j] + i);
    }
  }
  nin_ = minnin_;
  nout_ = minnout_;
  if (2 != nin_) {
    return 0.0;
  }

  /* Read process info: full cross section */
  readParagraphs (f, 4, rd_array2);
  checkPhysVal (thevar, &key, plist);

  /* process events */
  while (fgets (bufff[0], 2048, f)) {
    ++nEvents;
    sscanf (bufff[0], "%i%[^\n]", &npr, bufff[1]);
    if (npr < 1 || npr > Nproc) {
      fprintf (stderr, "mk_tab (error): unknown subprocess (Nproc = %i, maximum Nproc in the sample =%i), don't use the event \n", npr, Nproc);
      continue;
    }
    for (i = 0; i < 4 * nin_; i++) pvect[i] = 0.;

    sscanf (bufff[1], "    %lf %lf%[^a]", pvect + 3, pvect + 7, bufff[2]);

    for (i = nin_; i < nin_ + nout_; i++) {
      double *p = pvect + 4 * i;
      sscanf (bufff[i], " %lf %lf %lf%[^a]", p + 1, p + 2, p + 3, bufff[i + 1]);
    }

    for (i = 0; i < nin_ + nout_; i++) {
      pvect[4 * i] = ENERGY (mass[npr - 1][i], pvect + 4 * i + 1);
    }

    i = nbin * (calcPhysVal (key, plist, restfr) - minX) / (maxX - minX);
    if (0 <= i && i < nbin) {
      ++hist[i];
    }
  }

  for (i = 0; i < Nproc; i++) {
    free (mass[i]);
  }
  free (mass);

  return getTotCS () / nEvents;
}


static double 
tab_format0 (char fname[], char thevar[], int nbin, double minX, double maxX)
{
  /* Old CompHEP_41.10 fromat */
  int i;
  int npr;
  int nEvents = 0;
  char key;
  char restfr[10]={0};
  char plist[MAXNP];
  char buff[STRSIZ];
  FILE * f = fopen (fname, "r");

  double mass[MAXNP];

  rw_paragraph rd_array[11] = {
    {"CompHEP", NULL},
    {"PROCESS", old_getNames},
    {"Initial_state", NULL},
    {"MASSES", getMasses},
    {"Cross_section(Width)", readCS},
    {"Number_of_events", skipLine},
    {"----------------------------------------------------------",
     skipLine},
    {"Number_of_subprocesses", skipLine},
    {"Total_cross_section_(pb)", skipLine},
    {"Events_mixed_and_randomized", skipLine},
    {"Nproc", skipHeadLine}
  };
  readParagraphs (f, 11, rd_array);

  for (i = 0; i < nin_ + nout_; i++) {
    pinf_ (1, i + 1, NULL, mass + i);
  }
  if (2 != nin_) {
    return 0.0;
  }

  checkPhysVal (thevar, &key, plist);
  while (1 == fscanf (f, " %d", &npr)) {
    ++nEvents;
    for (i = 0; i < 4 * nin_; i++) pvect[i] = 0;
    fscanf (f, "      %lf %lf", pvect + 3, pvect + 7);
    for (i = nin_; i < nin_ + nout_; i++) {
      double *p = pvect + 4 * i;
      fscanf (f, " %lf %lf %lf", p + 1, p + 2, p + 3);
    }
    for (i = 0; i < nin_ + nout_; i++) {
      pvect[4 * i] = ENERGY (mass[i], pvect + 4 * i + 1);
    }
    i = nbin * (calcPhysVal (key, plist, restfr) - minX) / (maxX - minX);
    if (0 <= i && i < nbin) {
      ++hist[i];
    }
    fgets (buff, 1000, f);
  }
  fclose (f);

  return getCS () / nEvents;
}

static double 
tab_format_lhef (char fname[], char thevar[], int nbin, double minX, double maxX)
{
  int i, j;
  int err;
  int nEvents = 0;
  char key;
  char restfr[10]={0};
  char plist[MAXNP];
  eventUP evnt;
  long pos0;
  FILE * f = fopen (fname, "r");

#ifdef LIBXML
  formXMLtree (fname, 0);
  pos0 = getEventPosition (0);
  if (0 == pos0) {
    fprintf (stderr, "mk_tab (warning): it seems HepML header lost\n");
    rewind (f);
    setInfoWithoutHEPML (fname, f, 0);
    pos0 = getEventPosition (0);
  }
#else
  fprintf (stderr, "mk_tab (warning): can't process HepML header\n");
  setInfoWithoutHEPML (fname, f, 0);
  pos0 = getEventPosition (0);
#endif
  err = getNames (getSubprocessName ());
  checkPhysVal (thevar, &key, plist);

  while (0 == getLHAevent (fname, f, pos0, &evnt)) {
    pos0 = ftell (f);
    ++nEvents;
    for (i = 0; i < evnt.NpartUP; ++i) {
      for (j = 1; j < 4; ++j)
        pvect[4 * i + j] = evnt.momentumUP[j - 1][i];
      pvect[4 * i] = evnt.momentumUP[3][i];
    }
    i = nbin * (calcPhysVal (key, plist, restfr) - minX) / (maxX - minX);
    if (i >= 0 && i < nbin) {
      hist[i]  += 1.;
    }
  }
  fclose (f);

  return getCrossSection (0) / nEvents;
}


int 
prepare_tab (int frmt, char fname[], char thevar[], int nbin, double minX, double maxX)
{
  int i;
  char key;
  int Npoints = 0;
  double weight;
  char plist[MAXNP];

  /* book histogram */
  hist = malloc (nbin * sizeof (double));
  dhist = malloc (nbin * sizeof (double));
  for (i = 0; i < nbin; i++) {
    hist[i]  = 0.;
    dhist[i] = 0.;
  }

  /* prepare unweighted histogram */
  switch (frmt) {
    case 1:
      weight = tab_format0 (fname, thevar, nbin, minX, maxX);
      break;
    case 2:
      weight = tab_format2 (fname, thevar, nbin, minX, maxX);
      break;
    case 3:
      weight = tab_format1 (fname, thevar, nbin, minX, maxX);
      break;
    case 4:
      weight = tab_format_lhef (fname, thevar, nbin, minX, maxX);
      break;
    default:
      fprintf(stderr, "mk_tab (error): unknown event format in event file, can't prepare histogram\n");
      return 4;
    }

  /* rewighting */
  for (i = 0; i < nbin; i++) {
    Npoints += hist[i];
    hist[i] *= weight;
    dhist[i] = hist[i] * weight;
  }

  if (!Npoints) {
    fprintf(stderr, "mk_tab (error): histogram has no events\n");
    Npoints = 1;
  }
  weight = nbin / (maxX - minX);
  for (i = 0; i < nbin; i++) {
    dhist[i] = weight * sqrt (dhist[i] - hist[i] * hist[i] / Npoints);
    hist[i] *= weight;
  }

  /* build histogram */
  checkPhysVal (thevar, &key, plist);
  keep_histo (hist, dhist, key, plist, minX, maxX, nbin);
  free (hist);
  free (dhist);

  return 0;
}
