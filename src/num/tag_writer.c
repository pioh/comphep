/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_writer.h"

static FILE *evfile = NULL;

void 
tag_writer (FILE * file, elementary_tag * tag)
{
  int i, k;
  int len;

  fprintf (file, "##%s:", tag->tagname);
  len = strlen (tag->tagname);
  for (i = 0; i < tag->tagsize; i++)
    {
      fprintf (file, " %s,", tag->commands[i]);
      if (!((i + 1) % 5) && i != tag->tagsize - 1)
	{
	  fprintf (file, "\n");
	  for (k = 0; k < len + 3; k++)
	    {
	      fputc (' ', file);
	    }
	}
    }
  fprintf (file, ";\n");
}

int
cap_writer (FILE * file, tags * p)
{
  int i;

  for (i = 0; i < p->number_of_tags; i++)
    tag_writer (file, p->tag[i]);
  fprintf (file, ";\n");
  return 0;
}

int
change_cap (tags * t, int nEvents, double mult, double rmax, double cs, double er)
{
  int n;
  int N, Nold;
  double maxW, ml, old_ml;
  string_comnd com;

  n = get_tag (0, t, "total");
  if (n != -1)
    {
      Nold = get_ival (0, "Nevent", t->tag[n]);
      N = nEvents + Nold;
      strcpy (com.name, "Nevent");
      change_ival (N, 8 - intlen (N), com, t->tag[n]);
      strcpy (com.name, "CrosSec");
      change_fval (cs, com, t->tag[n]);
      strcpy (com.name, "CrosSecErr");
      change_fval (er, com, t->tag[n]);
    }
  else
    {
      fprintf (stderr, "***Error! The total tag has not been found...\n");
      return -1;
    }

  n = get_tag (0, t, "process");
  if (n != -1)
    {
      strcpy (com.name, "CrosSec");
      change_fval (cs, com, t->tag[n]);
      strcpy (com.name, "CrosSecErr");
      change_fval (er, com, t->tag[n]);
    }
  else
    {
      fprintf (stderr, "***Error! The process tag has not been found...\n");
      return -1;
    }

  n = get_tag (0, t, "n_event");
  if (n != -1)
    {
      strcpy (com.name, "N");
      change_ival (N, 8 - intlen (N), com, t->tag[n]);
      old_ml = get_fval (0, "mult", t->tag[n]);
      ml = (old_ml * (double) Nold + mult * (double) nEvents) / N;
      strcpy (com.name, "mult");
      change_fval (ml, com, t->tag[n]);
      maxW = get_ival (0, "maxW", t->tag[n]);
      strcpy (com.name, "maxW");
      if (maxW < rmax)
        change_fval (rmax, com, t->tag[n]);
    }
  else
    {
      fprintf (stderr, "***Error! The n_event tag has not been found...\n");
      return -1;
    }
  return 0;
}

static int 
write_mandatory_tag (int m, tags * the_tags, char *tagname)
{
  int n;

  n = get_tag (m, the_tags, tagname);
  if (n == -1)
    {
      fprintf (stderr, "Error: problems with mandatory tags during writing procedure");
      return -1;
    }
  tag_writer (evfile, the_tags->tag[n]);
  return n;
}

int 
write_cap (FILE * file, tags ** the_tags, process_ prUP)
{
  int i, j;
  int npr = 0;
  int nproc = 1;
  int n, nn;
  int Ntot = 0;
  string_comnd com[3];

  evfile = file;
  n = write_mandatory_tag (0, the_tags[0], "beam");
  n = write_mandatory_tag (n + 1, the_tags[0], "beam");
  n = write_mandatory_tag (0, the_tags[0], "strfun");
  n = write_mandatory_tag (n + 1, the_tags[0], "strfun");

  j = -1;
  for (i = 0; i < prUP.proc_info.NprocRUP; i++) {
    sprintf (com[0].value, "%i", nproc);
    sprintf (com[1].value, "%i", i + 1);
    if (i >= npr) {
      nproc = 1;
      j++;
      strcpy (com[0].name, "Nproc");
      n = get_tag_with1com (0, the_tags[j], "total", &(com[0]));
      if (n != -1) {
        npr += atoi (com[0].value);
      } else {
        fprintf (stderr, "1.Erorr during tags writing...\nExit.\n");
        exit (5);
        return -1;
      }
    } else
      nproc++;

    strcpy (com[0].name, "ID");
    sprintf (com[0].value, "%i", nproc);
    strcpy (com[1].name, "ID");

    n = get_tag_with_exactcom (0, the_tags[j], "process", com[0]);
    if (n == -1)
      {
        fprintf (stderr, "2.Erorr during tags writing...\nExit.\n");
        exit (5);
        return -1;
      }

    if (prUP.maps[i].used == 1) {
      replace_com (com[1], the_tags[j]->tag[n]);
      tag_writer (file, the_tags[j]->tag[n]);
      replace_com (com[0], the_tags[j]->tag[n]);

      strcpy (com[0].name, "IDprocess");
      strcpy (com[1].name, "IDprocess");
      n = get_tag_with_exactcom (0, the_tags[j], "generator", com[0]);
      if (n != -1) {
        replace_com (com[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        replace_com (com[0], the_tags[j]->tag[n]);
      }

      n = get_tag_with_exactcom (n + 1, the_tags[j], "n_event", com[0]);
      if (n != -1) {
        int N = get_ival (0, "N", the_tags[j]->tag[n]);
        Ntot += N;
        strcpy (com[2].name, "origN");
        sprintf (com[2].value, "%i", N);
        replace_com (com[1], the_tags[j]->tag[n]);
        nn = replace_com (com[2], the_tags[j]->tag[n]);
        if (nn == -1)
          add_com (com[2], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        replace_com (com[0], the_tags[j]->tag[n]);
      }

      n = get_tag_with_exactcom (n + 1, the_tags[j], "cut", com[0]);
      while (n != -1) {
        replace_com (com[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        replace_com (com[0], the_tags[j]->tag[n]);
        n = get_tag_with_exactcom (n + 1, the_tags[j], "cut", com[0]);
      }

      n = get_tag_with_exactcom (n + 1, the_tags[j], "parton", com[0]);;
      while (n != -1) {
        replace_com (com[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        replace_com (com[0], the_tags[j]->tag[n]);
        n = get_tag_with_exactcom (n + 1, the_tags[j], "parton", com[0]);
      }

      n = get_tag_with_exactcom (n + 1, the_tags[j], "QCDinfo", com[0]);
      if (n != -1) {
        replace_com (com[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        replace_com (com[0], the_tags[j]->tag[n]);
      }

      n = get_tag_with_exactcom (n + 1, the_tags[j], "format", com[0]);
      if (n != -1) {
        replace_com (com[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        replace_com (com[0], the_tags[j]->tag[n]);
      }
    }
  }

  n = get_tag (0, the_tags[0], "total");
  if (n != -1)
    {
      strcpy (com[2].name, "Nproc");
      sprintf (com[2].value, "%i", prUP.proc_info.NprocRUP);
      replace_com (com[2], the_tags[0]->tag[n]);
      strcpy (com[2].name, "Nevent");
      sprintf (com[2].value, "%i", Ntot);
      replace_com (com[2], the_tags[0]->tag[n]);
      tag_writer (file, the_tags[0]->tag[n]);
    }
  fprintf (file, "##:;;\n");

  return 0;
}

int 
write_cap_new (FILE * file, tags ** the_tags, int nf, process_ prUP)
{
  int i, j;
  int nproc = 1;
  int n;
  int Ntot = 0;
  int err = 0;
  int shift = 0;
  string_comnd com[2];
  string_comnd ncom[2];

  evfile = file;

  n = write_mandatory_tag (0, the_tags[0], "beam");
  n = write_mandatory_tag (n + 1, the_tags[0], "beam");
  n = write_mandatory_tag (0, the_tags[0], "strfun");
  n = write_mandatory_tag (n + 1, the_tags[0], "strfun");

  for (j = 0; j < nf && the_tags[j]; ++j) {
    strcpy (com[0].name, "Nproc");
    sprintf (com[0].value, "%i", 0);
    nproc = 0;
    if (get_tag_with1com (0, the_tags[j], "total", &(com[0])) != -1) {
      nproc = atoi (com[0].value);
    } else ++err;

    for (i = 0; i < nproc; i++) {
      strcpy (com[0].name, "ID");
      sprintf (com[0].value, "%i", i + 1);
      strcpy (com[1].name, "IDprocess");
      sprintf (com[1].value, "%i", i + 1);

      strcpy (ncom[0].name, "ID");
      sprintf (ncom[0].value, "%i", i + 1 + shift);
      strcpy (ncom[1].name, "IDprocess");
      sprintf (ncom[1].value, "%i", i + 1 + shift);

      n = get_tag_with_exactcom (0, the_tags[j], "process", com[0]);
      if (n != -1) {
        replace_com (ncom[0], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
      } else ++err;
      n = get_tag_with_exactcom (0, the_tags[j], "generator", com[1]);
      if (n != -1) {
        replace_com (ncom[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
      } else ++err;
      n = get_tag_with_exactcom (0, the_tags[j], "n_event", com[1]);
      if (n != -1) {
        int N = get_ival (0, "N", the_tags[j]->tag[n]);
        string_comnd lcom;
        Ntot += N;
        strcpy (lcom.name, "origN");
        sprintf (lcom.value, "%i", N);
        replace_com (ncom[1], the_tags[j]->tag[n]);
        if (replace_com (lcom, the_tags[j]->tag[n]) == -1)
          add_com (lcom, the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
      } else ++err;
      n = get_tag_with_exactcom (0, the_tags[j], "cut", com[1]);
      while (n != -1) {
        replace_com (ncom[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        n = get_tag_with_exactcom (n + 1, the_tags[j], "cut", com[1]);
      }
      n = get_tag_with_exactcom (0, the_tags[j], "parton", com[1]);;
      while (n != -1) {
        replace_com (ncom[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
        n = get_tag_with_exactcom (n + 1, the_tags[j], "parton", com[1]);
      }
      n = get_tag_with_exactcom (0, the_tags[j], "QCDinfo", com[1]);
      if (n != -1) {
        replace_com (ncom[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
      }
      n = get_tag_with_exactcom (0, the_tags[j], "format", com[1]);
      if (n != -1) {
        replace_com (ncom[1], the_tags[j]->tag[n]);
        tag_writer (file, the_tags[j]->tag[n]);
      } else ++err;
    }
    shift += nproc;
  }

  n = get_tag (0, the_tags[0], "total");
  if (n != -1) {
    strcpy (com[0].name, "Nproc");
    sprintf (com[0].value, "%i", shift);
    strcpy (com[1].name, "Nevent");
    sprintf (com[1].value, "%i", Ntot);

    replace_com (com[0], the_tags[0]->tag[n]);
    replace_com (com[1], the_tags[0]->tag[n]);
    tag_writer (file, the_tags[0]->tag[n]);
  } else ++err;
  fprintf (file, "##:;;\n");

  if (0 != err) {
    fprintf (stderr, "comphep (erorr): can't find necessary cpyth2 tags. Exit\n");
    exit (5);
  }

  return 0;
}


int 
final_write_cap (const char *outFileName, long *nEvent, int tot)
{
  int i, n;
  int Nproc;
  int N;
  int len;
  double *sc;
  double *scerr;
  double sc0 = 0.0;
  double sc1 = 0.0;
  FILE *file;
  tags *tagsbase;
  string_comnd com;

  tagsbase = init_cap (1);
  file = fopen (outFileName, "r+");
  cup_reader (file, tagsbase);
  fclose (file);

  strcpy (com.name, "Nproc");

  n = get_tag_with1com (0, tagsbase, "total", &com);
  Nproc = atoi (com.value);
  sc = malloc (Nproc * sizeof (double));
  scerr = malloc (Nproc * sizeof (double));

  for (i = 0; i < Nproc; i++) {
    strcpy (com.name, "IDprocess");
    sprintf (com.value, "%i", i + 1);
    n = get_tag_with_exactcom (0, tagsbase, "n_event", com);
    if (n != -1) {
      N = get_ival (0, "origN", tagsbase->tag[n]);
      strcpy (com.name, "N");
      sprintf (com.value, "%i", N);
      change_ival (nEvent[i], 0, com, tagsbase->tag[n]);
    }

    strcpy (com.name, "ID");
    sprintf (com.value, "%i", i + 1);
    n = get_tag_with_exactcom (0, tagsbase, "process", com);
    if (n != -1) {
      sc[i] = get_fval (0, "CrosSec", tagsbase->tag[n]);
      scerr[i] = get_fval (0, "CrosSecErr", tagsbase->tag[n]);
    }
  }

  for (i = 0; i < Nproc; i++) {
    sc0 += sc[i];
    sc1 += (scerr[i] * scerr[i]);
  }
  sc1 = sqrt (sc1);
  free (sc);
  free (scerr);

  n = get_tag (0, tagsbase, "total");
  if (n != -1) {
    N = get_ival (0, "Nevent", tagsbase->tag[n]);
    strcpy (com.name, "Nevent");
    sprintf (com.value, "%i", N);
    len = intlen (N) - intlen (tot);
    change_ival (tot, len, com, tagsbase->tag[n]);
    strcpy (com.name, "CrosSec");
    change_fval (sc0, com, tagsbase->tag[n]);
    strcpy (com.name, "CrosSecErr");
    change_fval (sc1, com, tagsbase->tag[n]);
  }

  file = fopen (outFileName, "w");
  for (i = 0; i < tagsbase->number_of_tags; i++) {
    tag_writer (file, tagsbase->tag[i]);
  }
  fprintf (file, "##:;;\n");
  fclose (file);

  return 1;
}

static int 
remove_subproc (tags * base, int num)
{
  int j;
  string_comnd com;

  strcpy (com.name, "IDprocess");
  sprintf (com.value, "%i", num);
  for (j = 0; j < base->number_of_tags; ++j) {
    if (tag_contain_exactcom (com, base->tag[j])) {
      remove_tag_num (j, base);
      --j;
    }
  }

  strcpy (com.name, "ID");
  for (j = 0; j < base->number_of_tags; ++j) {
    if (tag_contain_exactcom (com, base->tag[j]) && strcmp ("beam", base->tag[j]->tagname)) {
      remove_tag_num (j, base);
      --j;
    }
  }

  return 1;
}

int 
distilling (const char * inname)
{
  int i, j, n;
  int nproc = 0;
  int N;
  int num;
  int * newid;
  int * renamed;
  double cs = 0.;
  double cserr = 0.;
  char * outname;
  tags *tagsbase;
  string_comnd com;
  FILE *efile;
  FILE *dfile;

  efile = fopen (inname, "r");
  tagsbase = init_cap (1);
  cup_reader (efile, tagsbase);

  strcpy (com.name, "Nproc");
  if (0 < get_tag_with1com (0, tagsbase, "total", &com));
    nproc = atoi (com.value);

  num = 1;
  newid = malloc (nproc * sizeof (int));
  for (i = 0; i < nproc; i++) {
    strcpy (com.name, "IDprocess");
    sprintf (com.value, "%i", i + 1);
    n = get_tag_with_exactcom (0, tagsbase, "n_event", com);
    if (n != -1) {
      N = get_ival (0, "N", tagsbase->tag[n]);
      if (0 < N) {
        newid[i] = num;
        ++num;
      } else {
/*        fprintf (stdout, "no events in subproc %i\n", i + 1); */
        newid[i] = 0;
      }
    }
  }

  for (i = 0; i < nproc; i++)
    if (0 == newid[i])
      remove_subproc (tagsbase, i + 1);

  remove_null_tags (tagsbase);
  renamed = malloc (tagsbase->number_of_tags * sizeof (int));

  for (i = 0; i < tagsbase->number_of_tags; ++i) renamed[i] = 1;
  strcpy (com.name, "IDprocess");
  for (i = 0; i < tagsbase->number_of_tags; ++i) {
    if (-1 < tag_contain_com (&com, tagsbase->tag[i])) {
      renamed[i] = 0;
    }
  }
  for (i = 0; i < nproc; ++i) {
    if (0 < newid[i]) {
      sprintf (com.value, "%i", i + 1);
      for (j = 0; j < tagsbase->number_of_tags; ++j) {
        if (!renamed[j]) {
          if (tag_contain_exactcom (com, tagsbase->tag[j])) {
            change_ival (newid[i], 0, com, tagsbase->tag[j]);
            renamed[j] = 1;
          }
        }
      }
    }
  }

  for (i = 0; i < tagsbase->number_of_tags; ++i) renamed[i] = 1;
  strcpy (com.name, "ID");
  for (i = 0; i < tagsbase->number_of_tags; ++i) {
    if (-1 < tag_contain_com (&com, tagsbase->tag[i]) && strcmp ("beam", tagsbase->tag[i]->tagname)) {
      renamed[i] = 0;
    }
  }
  for (i = 0; i < nproc; ++i) {
    if (0 < newid[i]) {
      sprintf (com.value, "%i", i + 1);
      for (j = 0; j < tagsbase->number_of_tags; ++j) {
        if (!renamed[j]) {
          if (tag_contain_exactcom (com, tagsbase->tag[j])) {
            change_ival (newid[i], 0, com, tagsbase->tag[j]);
            renamed[j] = 1;
          }
        }
      }
    }
  }

  --num;
  for (i = 0; i < num; ++i) {
    sprintf (com.value, "%i", i + 1);
    n = get_tag_with_exactcom (0, tagsbase, "process", com);
    if (n != -1) {
      double lscerr = get_fval (0, "CrosSecErr", tagsbase->tag[n]);
      cs += get_fval (0, "CrosSec", tagsbase->tag[n]);
      cserr += lscerr * lscerr;
    }
  }
  n = get_tag (0, tagsbase, "total");
  if (n != -1) {
    strcpy (com.name, "Nproc");
    change_ival (num, 0, com, tagsbase->tag[n]);
    strcpy (com.name, "CrosSec");
    change_fval (cs, com, tagsbase->tag[n]);
    strcpy (com.name, "CrosSecErr");
    change_fval (sqrt (cserr), com, tagsbase->tag[n]);
  }

  outname = malloc ((strlen (inname) + 10) * sizeof (char));
  sprintf (outname, "%s.distill", inname);
  dfile = fopen (outname, "w");
  cap_writer (dfile, tagsbase);
  getc (efile);
  {
    int err;
    char oldev[2048];
    char newev[2048];
    while (1) {
      if (!fgets (oldev, 2048, efile))
        break;
      err = sscanf (oldev, "%i%[^a]", &i, newev);
      if (2 == err)
        fprintf (dfile, "%i%s", newid[i - 1], newev);
    }
  }
  fclose (dfile);
  fclose (efile);

  remove (inname);
  rename (outname, inname);
  
  fprintf (stdout, "comphep (info): distilling complete, records for %i subprocesses with no events are deleted\n", nproc - num);

  free (newid);
  free (renamed);
  free (outname);

  return 1;
}
