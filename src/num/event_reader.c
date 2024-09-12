/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"

#include "LesHouches.h"
#include "tag_reader.h"
#include "tag_parser.h"
#include "event_reader.h"

static char * event_comments;

char * get_event_comments (void)
{
  return event_comments;
}


int getLHAevent (char fname[], FILE * source, long startpos, eventUP * evt) {
  int i;
  int go = 1;
  char buff[4096];
  char * local_event_comments = NULL;

  int nup, idproc, id, stat;
  int moth1, moth2, color1, color2;
  double x, y, z, e, m, time, spin;
  double xwgt, aqed, aqcd, scale;

  fseek (source, startpos, SEEK_SET);

  if (!fgets (buff, 4096, source)) {
    if (feof (source)) return 1;
    fprintf (stderr, "comphep (error): %li - wrong position for LesHouches event in %s (no <event> start-tag)\n", startpos, fname);
    return -1;
  }
  if (!strcmp ("</LesHouchesEvents>\n", buff)) {
    return 1;
  }

  if (!strcmp ("<event>", buff)) {
    fprintf (stderr, "comphep (error): %li - wrong position for LesHouches event in %s (no <event> start-tag)\n", startpos, fname);
    return -2;
  }

  fgets (buff, 4096, source);
  if (feof (source)) {
    goto enderror;
  }
  if (6 != sscanf (buff, "%i %i %lf %le %le %le", &nup, &idproc, &xwgt, &scale, &aqed, &aqcd)) {
    fprintf (stderr, "comphep (error): wrong LesHouches event format in %s (error status 1)\n", fname);
    return -3;
  }
  evt->NpartUP    = nup;
  evt->IDprocUP   = idproc;
  evt->XweightUP  = xwgt;
  evt->QEDalphaUP = aqed;
  evt->QCDalphaUP = aqcd;
  evt->QscaleUP   = scale;

  for (i = 0; i < nup; ++i) {
    fgets (buff, 4096, source);
    if (feof (source)) {
      goto enderror;
    }
    if (13 != sscanf (buff, "%i %i %i %i %i %i %le %le %le %le %le %le %le", &id, &stat, &moth1, &moth2, &color1, &color2, &x, &y, &z, &e, &m, &time, &spin)) {
      fprintf (stderr, "comphep (error): wrong LesHouches event format in %s (error status 2)\n", fname);
      return -4;
    }
    evt->IDpartUP[i]       = id;
    evt->statusUP[i]       = stat;
    evt->motherUP[0][i]    = moth1;
    evt->motherUP[1][i]    = moth2;
    evt->colorflowUP[0][i] = color1;
    evt->colorflowUP[1][i] = color2;

    evt->momentumUP[0][i] = x;
    evt->momentumUP[1][i] = y;
    evt->momentumUP[2][i] = z;
    evt->momentumUP[3][i] = e;
    evt->momentumUP[4][i] = m;
    evt->timelifeUP[i]    = time;
    evt->spinUP[i]        = spin;
  }

  while (go) {
    fgets (buff, 4096, source);
    if (feof (source)) {
      goto enderror;
    }
    if ('#' != buff[0]) go = 0;
    else {
      if (NULL == local_event_comments) {
        local_event_comments = malloc ((strlen (buff) + 1) * sizeof (char));
        local_event_comments[0] = 0;
      } else {
        int len = strlen (local_event_comments) + 1;
        local_event_comments = realloc (local_event_comments, (len + strlen (buff)) * sizeof (char));
      }
      strcat (local_event_comments, buff);
    }
  }

  if (!strcmp ("</event>", buff)) {
    fprintf (stderr, "comphep (error): no </event> end-tag and no #-comment in %s\n", fname);
    return -5;
  }

  if (NULL != event_comments) {
    free (event_comments);
  }
  if (NULL != local_event_comments) {
    int len = strlen (local_event_comments) + 1;
    event_comments = malloc (len * sizeof (char));
    strcpy (event_comments, local_event_comments);
  }

  return 0;
enderror:
    fprintf (stderr, "comphep (error): unexpected end of LHE event file %s \n", fname);
    return -6;
}

int setLHAevent (char fname[], FILE * s, long startpos, eventUP * evt) {
  int i;
  int stat = 1;

  fseek (s, startpos, SEEK_SET);
  fprintf (s, "<event>\n");
  fprintf (s, "%i %i %17.10E %17.10E %17.10E %17.10E\n", 
    evt->NpartUP, evt->IDprocUP, evt->XweightUP, evt->QscaleUP, evt->QEDalphaUP, evt->QCDalphaUP);

  for (i = 0; i < evt->NpartUP; ++i) {
    if (i < 2) stat = -1;
    fprintf (s, "%i %i %i %i %i %i %17.10E %17.10E %17.10E %17.10E %17.10E %17.10E %17.10E\n", 
    evt->IDpartUP[i], 
    evt->statusUP[i], 
    evt->motherUP[0][i], 
    evt->motherUP[1][i], 
    evt->colorflowUP[0][i], 
    evt->colorflowUP[1][i], 
    evt->momentumUP[0][i], 
    evt->momentumUP[1][i], 
    evt->momentumUP[2][i], 
    evt->momentumUP[3][i], 
    evt->momentumUP[4][i], 
    evt->timelifeUP[i], 
    evt->spinUP[i]);
  }
  fprintf (s, "</event>\n");

  return 0;
}

int setLHAevent_with_comments (char fname[], FILE * s, long startpos, eventUP * evt, char comments[]) {
  int i;
  int stat = 1;

  fseek (s, startpos, SEEK_SET);
  fprintf (s, "<event>\n");
  fprintf (s, "%i %i %17.10E %17.10E %17.10E %17.10E\n", 
    evt->NpartUP, evt->IDprocUP, evt->XweightUP, evt->QscaleUP, evt->QEDalphaUP, evt->QCDalphaUP);

  for (i = 0; i < evt->NpartUP; ++i) {
    if (i < 2) stat = -1;
    fprintf (s, "%i %i %i %i %i %i %17.10E %17.10E %17.10E %17.10E %17.10E %17.10E %17.10E\n", 
    evt->IDpartUP[i], 
    evt->statusUP[i], 
    evt->motherUP[0][i], 
    evt->motherUP[1][i], 
    evt->colorflowUP[0][i], 
    evt->colorflowUP[1][i], 
    evt->momentumUP[0][i], 
    evt->momentumUP[1][i], 
    evt->momentumUP[2][i], 
    evt->momentumUP[3][i], 
    evt->momentumUP[4][i], 
    evt->timelifeUP[i], 
    evt->spinUP[i]);
  }
  if (NULL != comments)
    fputs (comments, s);
  fprintf (s, "</event>\n");

  return 0;
}


int moveLHAevent (int num, char fname[], FILE * source, long spos, FILE * target, long tpos, int zrandom) {
  int i;
  int go = 1;
  char buff[4096];

  int nup, idprup;
  double xwgtup, aqedup, aqcdup, scaleup;
  int id, stat;
  int moth1, moth2, color1, color2;
  double x, y, z, e, m, time, spin;

  fseek (source, spos, SEEK_SET);
  fseek (target, tpos, SEEK_SET);

  if (!fgets (buff, 4096, source)) {
    if (feof (source)) return 1;
    fprintf (stderr, "mix: (error) %li - wrong position for LesHouches event in %s (no <event> start-tag)\n", spos, fname);
    return -1;
  }

  if (!strcmp ("</LesHouchesEvents>", buff)) {
    return 1;
  }

  if (!strcmp ("<event>", buff)) {
    fprintf (stderr, "mix: (error) %li - wrong position for LesHouches event in %s (no <event> start-tag)\n", spos, fname);
    return -2;
  }
  fputs (buff, target);

  fgets (buff, 4096, source);
  if (feof (source)) {
    goto enderror;
  }
  if (6 != sscanf (buff, "%i %i %le %le %le %le", &nup, &idprup, &xwgtup, &scaleup, &aqedup, &aqcdup)) {
    fprintf (stderr, "mix: (error) wrong LesHouches event format in %s \n", fname);
    return -4;
  }
  if (num < 0) num = idprup;
  fprintf (target, "%i %i %17.10E %17.10E %17.10E %17.10E\n", nup, num, xwgtup, scaleup, aqedup, aqcdup);
/*
  for (i = 0; i < nup; ++i) {
    fgets (buff, 4096, source);
    if (feof (source)) goto enderror;
    fputs (buff, target);
  }
*/
  char buff2[4096];
  fgets (buff, 4096, source);
  if (feof (source)) goto enderror;
  fgets (buff2, 4096, source);
  if (feof (source)) goto enderror;
  if (zrandom) {
    if (13 != sscanf (buff2, "%i %i %i %i %i %i %le %le %le %le %le %le %le", &id, &stat, &moth1, &moth2, &color1, &color2, &x, &y, &z, &e, &m, &time, &spin)) {
      fprintf (stderr, "comphep (error): wrong LesHouches event format in %s \n", fname);
      return -4;
    }
    fprintf (target, "%i %i %i %i %i %i %17.10E %17.10E %17.10E %17.10E %17.10E %f %f\n", id, stat, moth1, moth2, color1, color2, x, y, -1. * z, e, m, time, spin);
    if (13 != sscanf (buff, "%i %i %i %i %i %i %le %le %le %le %le %le %le", &id, &stat, &moth1, &moth2, &color1, &color2, &x, &y, &z, &e, &m, &time, &spin)) {
      fprintf (stderr, "comphep (error): wrong LesHouches event format in %s \n", fname);
      return -4;
    }
    fprintf (target, "%i %i %i %i %i %i %17.10E %17.10E %17.10E %17.10E %17.10E %f %f\n", id, stat, moth1, moth2, color1, color2, x, y, -1. * z, e, m, time, spin);
  } else {
    if (13 != sscanf (buff, "%i %i %i %i %i %i %le %le %le %le %le %le %le", &id, &stat, &moth1, &moth2, &color1, &color2, &x, &y, &z, &e, &m, &time, &spin)) {
      fprintf (stderr, "comphep (error): wrong LesHouches event format in %s \n", fname);
      return -4;
    }
    fprintf (target, "%i %i %i %i %i %i %17.10E %17.10E %17.10E %17.10E %17.10E %f %f\n", id, stat, moth1, moth2, color1, color2, x, y, z, e, m, time, spin);
    if (13 != sscanf (buff2, "%i %i %i %i %i %i %le %le %le %le %le %le %le", &id, &stat, &moth1, &moth2, &color1, &color2, &x, &y, &z, &e, &m, &time, &spin)) {
      fprintf (stderr, "comphep (error): wrong LesHouches event format in %s \n", fname);
      return -4;
    }
    fprintf (target, "%i %i %i %i %i %i %17.10E %17.10E %17.10E %17.10E %17.10E %f %f\n", id, stat, moth1, moth2, color1, color2, x, y, z, e, m, time, spin);
  }

  for (i = 2; i < nup; ++i) {
    fgets (buff, 4096, source);
    if (feof (source)) {
      goto enderror;
    }
    if (13 != sscanf (buff, "%i %i %i %i %i %i %le %le %le %le %le %le %le", &id, &stat, &moth1, &moth2, &color1, &color2, &x, &y, &z, &e, &m, &time, &spin)) {
      fprintf (stderr, "comphep (error): wrong LesHouches event format in %s \n", fname);
      return -4;
    }
    if (zrandom) z *= -1.;
    fprintf (target, "%i %i %i %i %i %i %17.10E %17.10E %17.10E %17.10E %17.10E %f %f\n", id, stat, moth1, moth2, color1, color2, x, y, z, e, m, time, spin);
  }

  while (go) {
    fgets (buff, 4096, source);
    if (feof (source)) {
      goto enderror;
    }
    if ('#' != buff[0]) {
      go = 0;
    } else {
      fputs (buff, target);
    }
  }

  if (!strcmp ("</event>", buff)) {
    fprintf (stderr, "mix: (error) no </event> end-tag in %s\n", fname);
    return -4;
  }
  fputs (buff, target);

  return 0;
enderror:
    fprintf (stderr, "mix: (error) unexpected end of LHE event file %s \n", fname);
    return -6;
}

int testLHAevent (FILE * source, char fname[], long pos) {
  int i;
  int go = 1;
  char buff[4096];

  int nup, idproc, id, stat;
  int moth1, moth2, color1, color2;
  double x, y, z, e, m, time, spin;
  double xwgt, aqed, aqcd, scale;

  fseek (source, pos, SEEK_SET);
  if (!fgets (buff, 4096, source)) {
    if (feof (source)) return 1;
    fprintf (stderr, "mix: (error) %li - wrong position for LesHouches event in %s (no <event> start-tag)\n", pos, fname);
    return -1;
  }

  if (!strcmp ("</LesHouchesEvents>\n", buff)) {
    return 1;
  }

  if (!strcmp ("<event>", buff)) {
    fprintf (stderr, "mix: (error) %li - wrong position for LesHouches event in %s (no <event> start-tag)\n", pos, fname);
    return -2;
  }

  fgets (buff, 4096, source);
  if (feof (source)) {
    goto enderror;
  }
  if (6 != sscanf (buff, "%i %i %lf %le %le %le", &nup, &idproc, &xwgt, &scale, &aqed, &aqcd)) {
    fprintf (stderr, "mix: (error) wrong LesHouches event format in %s \n", fname);
    return -3;
  }

  for (i = 0; i < nup; ++i) {
    fgets (buff, 4096, source);
    if (feof (source)) {
      goto enderror;
    }
    if (13 != sscanf (buff, "%i %i %i %i %i %i %le %le %le %le %le %le %le", &id, &stat, &moth1, &moth2, &color1, &color2, &x, &y, &z, &e, &m, &time, &spin)) {
      fprintf (stderr, "mix: (error) wrong LesHouches event format in %s \n", fname);
      return -4;
    }
  }

  while (go) {
    fgets (buff, 4096, source);
    if (feof (source)) {
      goto enderror;
    }
    if ('#' != buff[0]) go = 0;
  }

  if (!strcmp ("</event>", buff)) {
    fprintf (stderr, "mix: (error) no </event> end-tag and no #-comment in %s\n", fname);
    return -5;
  }

  return 0;
enderror:
    fprintf (stderr, "mix: (error) unexpected end of LHE event file %s \n", fname);
    return -6;
}
