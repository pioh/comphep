/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/4_vector.h"

#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_routines.h"

elementary_tag **
preparer (FILE * f)
{
  int n, i;
  int Nproc;
  tags *head;
  elementary_tag **format;
  string_comnd com;

  head = init_cap (1);
  cup_reader (f, head);

  strcpy (com.name, "Nproc");
  get_tag_with1com (0, head, "total", &com);
  Nproc = atoi (com.value);

  strcpy (com.name, "IDprocess");
  n = -1;
  for (i = 1; i <= Nproc; ++i) {
    sprintf (com.value, "%i", i);
    n = get_tag_with_exactcom (n + 1, head, "format", com);
    if (n != -1) {
      format[i] = head->tag[n];
    } else {
      fprintf (stderr, "tag_parser (error): the format tag for the %i-th subprocess has been lost.\n", i);
      exit (1);
    }
  }

  return format;
}

char * primetrim (char * p)
{
  int i;
  int len;
  char * val = strchr (p, '\'');  /* remove first prime */
  strcpy (p, val + 1);
  len = strlen (p);
  for (i = len; '\'' != p[i]; --i)
    continue;
  p[i] = 0;  /* remove last prime */
  return p;
}


int check_cpyth2 (tags * head)
{
  int i, j, k;
  int n, m;
  int format;
  int Npartons;
  int Nevents;
  double CSec;
  int Nproc = -1;
  string_comnd com;
  int locnin;
  int locnout;
  int the_type = 0;
  char type_proc[128];

  strcpy (com.name, "Nproc");
  get_tag_with1com (0, head, "total", &com);
  Nproc = atoi (com.value);

  if (0 > Nproc) {
    return -1;
  }

  for (i = 0; i < Nproc; i++) {
    shortstr s;
    sprintf (s, "%i", i + 1);

    strcpy (com.name, "ID");
    strcpy (com.value, s);
    n = get_tag_with_exactcom (0, head, "process", com);
    if (n == -1)
      {
        fprintf (stderr, "tag_parser (error): the process tag is lost (%i-th subprocess)\n", i + 1);
        return -1;
      }
    Npartons = get_ival (0, "Nparton", head->tag[n]);
    CSec = get_fval (0, "CrosSec", head->tag[n]);
    strcpy (type_proc, get_cval (0, "type", head->tag[n]));
    if (0 == strcmp (type_proc, "\'decay\'")) the_type = 1;
    if (0 == strcmp (type_proc, "\'scattering\'")) the_type = 2;

    strcpy (com.name, "IDprocess");
    n = get_tag_with_exactcom (0, head, "n_event", com);
    if (n == -1) {
      fprintf (stderr, "tag_parser (error): the n_event tag is lost (%i-th subprocess)\n", i + 1);
      return -1;
    }
    Nevents = get_ival (0, "N", head->tag[n]);

    if (Npartons < 0 || CSec < 0. || Nevents < 0) {
      fprintf (stderr, "tag_parser (error): strange values in process ot n_event\n");
      return -1;
    }

    sprintf (com.value, "%i", i + 1);
    format = get_tag_with_exactcom (0, head, "format", com);
    if (-1 == format) {
      fprintf (stderr, "tag_parser (error): the format tag is lost (%i-th subprocess)\n", i + 1);
      return -1;
    }

    if (2 == the_type) {
      strcpy (com.name, "p1.3");
      n = tag_contain_com (&com, head->tag[format]);
      strcpy (com.name, "p2.3");
      m = tag_contain_com (&com, head->tag[format]);
      if (n < 2 || m < 2) {
        fprintf (stderr, "tag_parser (error): wrong format tag (%i-th subprocess)\n", i + 1);
        return -1;
      }

      for (k = 2; k < Npartons; ++k) {
        for (j = 0; j < 3; ++j) {
          sprintf (com.name, "p%i.%i", k + 1, j + 1);
          n = tag_contain_com (&com, head->tag[format]);
          if (n < 2) {
            fprintf (stderr, "tag_parser (error): strange format tag (%i-th subprocess)\n", i + 1);
            return -1;
          }
        }
      }
    } else {
      for (k = 1; k < Npartons; ++k) {
        for (j = 0; j < 3; ++j) {
          sprintf (com.name, "p%i.%i", k + 1, j + 1);
          n = tag_contain_com (&com, head->tag[format]);
          if (n < 2) {
            fprintf (stderr, "tag_parser (error): strange format tag (%i-th subprocess)\n", i + 1);
            return -1;
          }
        }
      }
    }
    strcpy (com.name, "IDprocess");
    sprintf (com.value, "%i", i + 1);
    n = -1;
    locnin = locnout = 0;
    for (j = 0; j < Npartons; j++) {
      int in, out;
      n = get_tag_with_exactcom (n + 1, head, "parton", com);
      if (n == -1) {
        fprintf (stderr,
                 "tag_parser (error): The number of parton tags is not equal to the number\n"
                 "                     of parton in process tag (%i-th subprocess).\n", i + 1);
        return -1;
      }
      in  = get_ival (0, "in", head->tag[n]);
      out = get_ival (0, "out", head->tag[n]);
      if ((0 == in && 0 == out) || (0 != in && 0 != out)) {
        fprintf (stderr, "tag_parser (error): strange parton tag (%i-th subprocess).\n", i + 1);
        return -1;
      }
      if (in) ++locnin;
      if (out) ++locnout;
    }
    if (1 > locnin || 2 < locnin || 2 > locnout || 10 < locnout) {
      fprintf (stderr, "tag_parser (error): strange in/out partons: n(in) = %i, n(out) = %i in %i-th subprocess).\n", locnin, locnout, i + 1);
      return -1;
    }
  }

return 0;
}
