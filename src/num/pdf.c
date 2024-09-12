/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1999, Alexander Pukhov 
* ----------------------------------------------------
*/
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "service2/include/chep_limits.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"

#include "pdf.h"

static double lambda = -0.999;

int CERNpdf_number (char *pdf, char *ver, int *PDFid, int *PDFgroup) {
  if (0 == *PDFid && 0 == *PDFgroup) {
    if (!strcmp (pdf, "CTEQ") || !strcmp (pdf, "Cteq") || !strcmp (pdf, "cteq")) {
      *PDFgroup = 4;
      if (!strcmp (ver, "4m") || !strcmp (ver, "4M"))  *PDFid = 34;
      if (!strcmp (ver, "4l") || !strcmp (ver, "4L"))  *PDFid = 32;
      if (!strcmp (ver, "5m1")|| !strcmp (ver, "5M1")) *PDFid = 53;
      if (!strcmp (ver, "5m") || !strcmp (ver, "5M"))  *PDFid = 48;
      if (!strcmp (ver, "5l") || !strcmp (ver, "5L"))  *PDFid = 46;
      if (!strcmp (ver, "6d") || !strcmp (ver, "6D"))  *PDFid = 56;
      if (!strcmp (ver, "6l") || !strcmp (ver, "6L"))  *PDFid = 55;
      if (!strcmp (ver, "6m") || !strcmp (ver, "6M"))  *PDFid = 57;
      if (!strcmp (ver, "6l1")|| !strcmp (ver, "6L1")) *PDFid = 58;
      if (0 != *PDFid) return 0;
    }
  } else {
    if (4 == *PDFgroup) {
      strcpy (pdf, "cteq");
      switch (*PDFid) {
        case 34: strcpy (ver, "4m");  break;
        case 32: strcpy (ver, "4l");  break;
        case 53: strcpy (ver, "5m1"); break;
        case 48: strcpy (ver, "5m");  break;
        case 46: strcpy (ver, "5l");  break;
        case 55: strcpy (ver, "6l");  break;
        case 56: strcpy (ver, "6d");  break;
        case 57: strcpy (ver, "6m");  break;
        case 58: strcpy (ver, "6l1"); break;
        default: return 1;
      }
      return 0;
    }
  }
  return 1;
}

void delPdfList (pdfList * list) {
  while (list) {
    pdfList *next = list->next;;
    free (list->pathfile);
    free (list->file);
    free (list->name);
    free (list);
    list = next;
  }
}


void makePdfList (char * pathfile, char * file, char * parton, pdfList ** list) {
  shortstr s;
  shortstr dName;
  int pos;

  FILE * f = fopen (pathfile, "r");
  if (!f) {
    fprintf (stderr,"Warning! can not open file with pdf data: %s\n", pathfile);
    return;
  }

  while (1 == fscanf (f, "%s", s)) {
    if ('#' == s[0]) {
      if (strcmp (s + 1, "distribution")) {   /* NOT distribution*/
        fclose (f);
        return;
      }
      fscanf (f, " \"%[^\"]%*c", dName);
      pos = 1;
    } else {
      char * p;
      if (s[0] == '(')
        while (s[strlen (s) - 1] != ')')
          fscanf (f, "%s", s + strlen (s));

      p = strstr (s, parton);  /* find parton name in pdf name string */
      if (p && (p == s || p[-1] == '(' || p[-1] == ',')) {
        p += strlen (parton);
        if (p[0] == 0 || p[0] == ')' || p[0] == ',') {
          pdfList * new = malloc (sizeof (pdfList));

          new->name = malloc (strlen (dName) + 1);
          new->pathfile = malloc (strlen (pathfile) + 1);
          new->file = malloc (strlen (file) + 1);

          strcpy (new->name, dName);
          strcpy (new->pathfile, pathfile);
          strcpy (new->file, file);
          new->position = pos;

          new->next = *list;
          *list = new;
        }
      }
      pos++;
    }
  }

  fclose (f);
  return;
}

void comphepPdfList (char * p_name, pdfList ** list) {
  int i;
  longstr dNames[3];

  strcpy (dNames[0], "../");
  strcpy (dNames[1], "./");
  sprintf (dNames[2], "%s%s", pathtocomphep, "strfun/");

  for (i = 0; i < 3; i++) {
    int len;
    int err;
    searchrec s;
    longstr pathfile;

    sprintf (pathfile, "%s%s", dNames[i], "*.pdf");
    err = find_first (pathfile, &s, anyfile);

    len = strlen (dNames[i]);
    while (!err) {
      pathfile[len] = 0;
      strcat (pathfile, s.name);
      makePdfList (pathfile, s.name, p_name, list);
      err = find_next (&s);
    }
  }
}


void freePdfData (pdfStr * data) {
  if (data->x_grid);
  data->x_grid = NULL;
  if (data->q_grid);
  data->q_grid = NULL;
  if (data->alpha);
  data->alpha = NULL;
  if (data->strfun);
  data->strfun = NULL;
  if (data->filename);
  data->filename = NULL;
}

int getPdfData (char * pathfile, char * file, int n_parton, pdfStr * data) {
  shortstr pattern;
  shortstr buff;
  char c;
  int nx = 0;
  int nq = 0;
  int xMinOn = 0;
  int xMaxOn = 0;
  int qMinOn = 0;
  int qMaxOn = 0;

  FILE * f = fopen (pathfile, "r");
  if (!f) {
    fprintf (stderr,"Warning! can not open file with pdf data: %s\n", pathfile);
    return -1;
  }

  data->nq = 0;
  data->nx = 0;
  data->number = n_parton;
  data->x_grid = NULL;
  data->q_grid = NULL;
  data->alpha = NULL;
  data->strfun = NULL;
  data->mass = 1;
  data->lin = 0;
  data->filename = NULL;

  data->filename = malloc (strlen (file) * sizeof (char));
  strcpy (data->filename, file);

  sprintf (pattern, "%d-parton", n_parton);

  while (1 == fscanf (f, "%c", &c)){
    if (c == '#') {
      double qq;
      int i;

      fscanf (f, "%s", buff);
      if (!strcmp (buff, "Mass")) {
        if (1 != fscanf (f, "%lf", &data->mass)) {
          goto errexit;
        }
      } else if (!strcmp (buff, "Lambda")) while (fscanf (f, "%lf", &lambda));
       else if (!strcmp (buff, "Q_grid")) {
        long fpos = ftell (f);

        if (data->q_grid || data->strfun) goto errexit;
        while (fscanf (f, "%lf", &qq)) nq++;
        if (nq < 3) goto errexit;

        data->nq = nq;
        data->q_grid = malloc (nq * sizeof (double));
        fseek (f, fpos, SEEK_SET);
        for (i = 0; i < nq; i++) {
          fscanf (f, "%lf", data->q_grid + i);
        }
        if (data->q_grid[0] <= 0) {
          goto errexit;
        }
        data->q_grid[0] = log (data->q_grid[0]);
        for (i = 1; i < nq; i++) {
          if (data->q_grid[i - 1] >= data->q_grid[i]) {
            goto errexit;
          } else {
            data->q_grid[i] = log (data->q_grid[i]);
          }
        }
      }
      else if (!strcmp (buff, "lExtrapolation"))
        data->lin = 1;
      else if (!strcmp (buff, "X_grid"))
        {
          long fpos = ftell (f);

          if (data->x_grid)
            goto errexit;
          while (fscanf (f, "%lf", &qq))
            nx++;
          data->nx = nx;
          if (nx < 3)
            goto errexit;
          data->x_grid = malloc (nx * sizeof (double));
          fseek (f, fpos, SEEK_SET);
          for (i = 0; i < nx; i++)
            fscanf (f, "%lf", data->x_grid + i);
          for (i = 1; i < nx; i++)
            if (data->x_grid[i - 1] >= data->x_grid[i])
              goto errexit;
        }
      else if (!strcmp (buff, "Alpha"))
        {
          if (!data->q_grid)
            goto errexit;
          data->alpha = malloc (nq * sizeof (double));
          for (i = 0; i < nq; i++)
            if (fscanf (f, "%lf", data->alpha + i) != 1)
              goto errexit;
          if (fscanf (f, "%lf", &qq) == 1)
            goto errexit;
        }
      else if (!strcmp (buff, pattern))
        {
          int nn = nq ? nq * nx : nx;

          data->strfun = malloc (nn * sizeof (double));
          for (i = 0; i < nn; i++) {
            if (fscanf (f, "%lf", data->strfun + i) != 1) {
              goto errexit;
            }
          }
          if (fscanf (f, "%lf", &qq) == 1) {
            goto errexit;
          }

          switch (fscanf (f, "powers %lf %lf", &(data->pow0), &(data->pow1)))
            {
            case 0:
              data->pow0 = 0;
            case 1:
              data->pow1 = 0;
            }
          break;
      } else if (!strcmp (buff, "x_min")) {
          if (fscanf (f, "%lf", &(data->x_min)) != 1)
            goto errexit;
          xMinOn = 1;
      } else if (!strcmp (buff, "x_max")) {
          if (fscanf (f, "%lf", &(data->x_max)) != 1)
            goto errexit;
          xMaxOn = 1;
      } else if (!strcmp (buff, "q_min")) {
          if (fscanf (f, "%lf", &(data->q_min)) != 1)
            goto errexit;
          qMinOn = 1;
      } else if (!strcmp (buff, "q_max")) {
          if (fscanf (f, "%lf", &(data->q_max)) != 1)
            goto errexit;
          qMaxOn = 1;
      }
    }
  }
  fclose (f);

  if (!xMinOn) {
    data->x_min = data->x_grid[0];
  }

  if (!xMaxOn) {
    data->x_max = data->x_grid[data->nx - 1];
  }

  if (data->q_grid && !qMinOn) {
    data->q_min = data->q_grid[0];
  }

  if (!qMaxOn) {
    data->q_max = data->q_grid[data->nq - 1];
  }

  return 0;
errexit:
  {
    int errpos = ftell (f);
    fclose (f);
    freePdfData (data);
    return errpos;
  }
}


static double inte2 (double * xa, double * ya, double x) {
  double x0 = x - xa[0];
  double x1 = x - xa[1];
  return (ya[0] * x1 - ya[1] * x0) / (xa[0] - xa[1]);
}


static double inte3 (double * xa, double * ya, double x) {
  double x0 = x - xa[0];
  double x1 = x - xa[1];
  double x2 = x - xa[2];
  double x01 = 1. / (xa[0] - xa[1]);
  double x02 = 1. / (xa[0] - xa[2]);
  double x12 = 1. / (xa[1] - xa[2]);

  return ya[0] * x1 * x2 * x01 * x02 - ya[1] * x0 * x2 * x01 * x12 + ya[2] * x0 * x1 * x02 * x12;
}

static int leftX (int dim, double * xa, double x) {
  int k1, k2, k3;

  if (x < xa[0]) return 0;
  if (x > xa[dim - 1]) return dim - 3;

  k1 = 0;
  k2 = dim - 1;

  while (k2 - k1 > 1) {
    k3 = (k1 + k2) / 2;
    if (xa[k3] > x) {
      k2 = k3;
    } else {
        k1 = k3;
    }
  }

  k3 = k1;
  if (k3 < 0) k3 = 0;
  if (k3 > dim - 2) k3 = dim - 2;

  return k3;
}


double interFunc (double x, double q, pdfStr * W) {
  int px = leftX (W->nx, W->x_grid, x);
  double tmp[3];
  double dat = 0.0;

/* check Q-x range */
  if (W->q_grid && (q > W->q_max || q < W->q_min)) {
    fprintf (stderr, "Warning! Q=%E out of range\n", q);
  }
  if (x > W->x_max || x < W->x_min) {
    fprintf (stderr, "Warning! X=%E out of range\n", x);
  }

  if (W->lin) {
    dat = inte2 (W->x_grid + px, W->strfun + px, x);
    if (W->nq) {
      int i;
      double logQ = log (q);
      int pq = leftX (W->nq, W->q_grid, logQ);

      for (i = 0; i < 2; i++) {
        tmp[i] = inte2 (W->x_grid + px, W->strfun + W->nx * (pq + i) + px, x);
      }
      dat = inte2 (W->q_grid + pq, tmp, logQ);
    }
  } else {
    if (px > W->nx - 3) {
      px = W->nx - 3;
    }
    dat = inte3 (W->x_grid + px, W->strfun + px, x);
    if (W->nq) {
      int i;
      double logQ = log (q);
      int pq = leftX (W->nq, W->q_grid, logQ);

      if (pq > W->nq - 3) {
        pq = W->nq - 3;
      }

      for (i = 0; i < 3; i++) {
        tmp[i] = inte3 (W->x_grid + px, W->strfun + W->nx * (pq + i) + px, x);
      }
      dat = inte3 (W->q_grid + pq, tmp, logQ);
    }
  }

  return dat;
}

double interAlpha (double q, pdfStr * W) {
  double logQ = log (q);
  int pq = leftX (W->nq, W->q_grid, logQ);
  double dat = 0.0;

  if (W->lin) {
      dat = inte2 (W->q_grid + pq, W->alpha, logQ);
  } else {
    if (pq > W->nq - 3) pq = W->nq - 3;
    dat = inte3 (W->q_grid + pq, W->alpha + pq, logQ);
  }

  return dat;
}

double pdf_QCDLambda (void) {
  return lambda;
}
