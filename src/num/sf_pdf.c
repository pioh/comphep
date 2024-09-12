/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* with changes by V.Ilyin for #-models, 6 Sept 2000
* ------------------------------------------------------
*/
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"
#include "service2/include/syst.h"
#include "chep_crt/include/chep_crt.h"

#include "tag_reader.h"
#include "tag_parser.h"
#include "strfun_par.h"
#include "pdf.h"
#include "sf_pdf.h"

static int pdfset[2] = {0, 0};
static int pdfmem[2] = {0, 0};
static shortstr pdfName[2] = {"PDF:OFF", "PDF:OFF"};
static shortstr pdfBeam[2] = {"parton", "parton"};

static pdfStr pdfData[9];       /*  #-mdl */

int p_pdf (char * p_name) {
  int val = 0;
  pdfList * list = NULL;
  comphepPdfList (p_name, &list);
  if (list) {
    delPdfList (list);
    val = 1;
  }
  return val;
}

void info_pdf (int i, Str_fun_Info * info) {
  int PDFid = 0;
  int PDFgr = 0;
  int unknown_pdf;

  strcpy (info->pdf_name, "CTEQ");
  info->prt_mass = 0.938;
  strcpy (info->prt_name, pdfBeam[i]);
  strcpy (info->version, pdfName[i] + 4 + 4); /* remove 'PDF:' and 'cteq' */

  unknown_pdf = CERNpdf_number (info->pdf_name, info->version, &PDFid, &PDFgr);
  if (unknown_pdf) {
    fprintf (stderr, "*** Uknown PDF %s or its version %s", info->pdf_name, info->version);
    exit (1);
  }

  info->N_extra_commands = 2;
  snprintf (info->extra_commands[0], 7 + intlen (PDFid), "PDFid=%i", PDFid);
  snprintf (info->extra_commands[1], 7 + intlen (PDFgr), "PDFgr=%i", PDFgr);
  info->PDFLIBset = PDFid;
  info->PDFLIBgroup = PDFgr;
  info->LHAPDFset = 0;
  info->LHAPDFmember = 0;
}

void n_pdf (int i, char * beam, char * pdf) {
  i--;
  if (strlen(pdfName[i])) {
    strcpy (pdf, pdfName[i]);
      sprintf (pdf, "%s:%i:%i", pdfName[i], pdfset[i], pdfmem[i]);
  } else {
    strcpy (pdfBeam[i], "parton");
    strcpy (pdf, "PDF:OFF:0:0");
  }
  strcpy (beam, pdfBeam[i]);
}

static int 
extract_names (int i, char name[]) {
  int len;

  char * bname = strstr (name, "(");
  strcpy (pdfBeam[i], bname + 1);
  len = strlen (pdfBeam[i]);
  pdfBeam[i][len - 1] = '\0';
  trim (pdfBeam[i]);

  len = strlen (name) - len - 1;
  name[len] = '\0';
  strcpy (pdfName[i], name);
  trim (pdfName[i]);

  return 1;
}

int r_pdf (int i, char * name) {
  --i;
  if (strstr (name, "PDF:")) {
    int set, mem;
    char name1[256];
    extract_names (i, name);
    if (get_sf_info (pdfName[i], name1, &set, &mem)) {
      strcpy (pdfName[i], name1);
      pdfset[i] = set;
      pdfmem[i] = mem;
    }
    return 1;
  }
  return 0;
}


/* initialization function 
  i - beam number
  be - return value (???)
  mass - beam particle mass
  return the string number in PDF menu 
*/

int be_pdf (int i, double * be, double * mass, char * p_name) {
  pdfList *list = NULL;
  pdfList *list_;
  int ret_code = 0;
  static int first[2] = {1, 1};
  int ip, xip;                  /* #-mdl */
  char *xpar;                   /* #-mdl */
  xpar = "anti-proton";         /* #-mdl */
  shortstr tmpname;

  i--;
  if (!first[i])
    freePdfData (pdfData + i);

  comphepPdfList (p_name, &list);
  list_ = list;
  sprintf (tmpname, "%s(%s)", pdfName[i] + 4, pdfBeam[i]);

  while (list_ && ret_code == 0) {
    if (!strcmp (list_->name, tmpname)) {
      if (0 == i || 0 != strlen(pdfName[0]))
        for (ip = 2; ip < 9; ++ip){
          xip = ip;
          if (strstr (pdfBeam[i], xpar)) {
            if (ip == 8) xip = 4;
            if (ip == 4) xip = 8;
            if (ip == 7) xip = 5;
            if (ip == 5) xip = 7;
          }
          getPdfData (list_->pathfile, list_->file, xip, pdfData + ip);   /* #-mdl */
        }
      ret_code = !getPdfData (list_->pathfile, list_->file, list_->position, pdfData + i);
    }
    list_ = list_->next;
  }

  if (list) {
    delPdfList (list);
  }

  if (ret_code) {
    first[i] = 0;
    if (pdfData[i].alpha)
      set_alphaMode (i + 1);
    if (pdfData[i].pow1 < 0) {
      *be = pdfData[i].pow1 + 1;
    } else {
      *be = 1.;
    }
    *mass = pdfData[i].mass;
  }

  return ret_code;
}

/* graphical interface. should be excluded! */
/* return the string number in PDF menu */
int pdf_pdf (int i, char * p_name) {
  int WIDTH = 30;
  int k = 0;
  char * menustring = NULL;
  void *pscr = NULL;
  pdfList *list = NULL;
  pdfList *list_ = NULL;

  comphepPdfList (p_name, &list_);
  if (list_) {
    list = list_;
  } else {
    return 0;
  }

  for (list_ = list; list_; list_ = list_->next) ++k;
  menustring = malloc (2 + WIDTH * k);
  menustring[0] = WIDTH;
  menustring[1] = 0;
  for (list_ = list; list_; list_ = list_->next) {
    sprintf (menustring + strlen (menustring), " PDF:%-*.*s", WIDTH - 5, WIDTH - 5, list_->name);
  }

  k = 0;
  menu1 (80 - 4 - WIDTH, 7, "", menustring, "", &pscr, &k);
  if (k) {
    int j;
    int l = k - 1;
    int sepa, len;
    int id = 0, gr = 0;
    for (list_ = list; l; --l) list_ = list_->next;
    len = strlen(list_->name);
    for (j = 0; j < len; ++j) {
      if ('(' == list_->name[j]) sepa = j;
    }
    --i;
    strcpy (pdfName[i], list_->name);
    sprintf (pdfName[i], "PDF:%s", list_->name);
    pdfName[i][sepa + 4] = 0;
    strcpy (pdfBeam[i], list_->name + sepa + 1);
    pdfBeam[i][len - sepa - 2] = 0;
    trim(pdfName[i]);
    trim(pdfBeam[i]);
    pdfset[i] = 0;
    pdfmem[i] = 0;
    if (!CERNpdf_number ("CTEQ", pdfName[i] + 4 + 4, &id, &gr)) {
      pdfset[i] = gr;
      pdfmem[i] = id;
    }
  }

  free (menustring);
  put_text (&pscr);
  delPdfList (list);
  return k;
}

int beam_pdf (int i) {
  int key = 0;
  void *pscr = NULL;
  shortstr previous;
  char strmen[] = "*" 
                  " parton      "
                  " proton      "
                  " anti-proton ";

  strmen[0] = (strlen (strmen) -1) / 3;
  --i;
  strcpy(previous, pdfBeam[i]);
  menu1 (64, 7, "", strmen, "n_beam_*", &pscr, &key);

  if (key) {
    switch (key) {
      case 1:
        strcpy (pdfName[i], "OFF");
        strcpy (pdfBeam[i], "parton");
        break;
      case 2:
        strcpy (pdfBeam[i], "proton");
        break;
      case 3:
        strcpy (pdfBeam[i], "anti-proton");
        break;
    }
  }
  if (0 == key && !strcmp("parton", pdfBeam[i])) {
    strcpy (pdfName[i], "OFF");
    pdfset[i] = 0;
    pdfmem[i] = 0;
    set_sf_num(i, 0);
  }
  if (!strcmp("proton", pdfBeam[i]) && !strcmp("parton", previous)) {
    pdfset[i] = 0;
    pdfmem[i] = 0;
    strcpy (pdfName[i], "PDF:undefined");
  }

  return 0;
}


/* return pdf value of i-th parton in a point (x,Q) */
/* stability interface to double interFunc (double x, double q, pdfStr * W) from pdf.c */
double c_pdf (int i, double x, double q) {
  double r;

  i--;
  r = interFunc (x, q, pdfData + i);
  if (!r) {
    return 0.;
  }

  if (pdfData[i].pow0) {
    r *= pow (x, pdfData[i].pow0);
  }

  if (pdfData[i].pow1 > 0) {
    r *= pow (1 - x, pdfData[i].pow1);
  }

  if (r < 0) {
    r = 0;
  }

  return r;
}

double alpha_pdf (double q) {
  return interAlpha (q, pdfData + get_alphaMode() - 1);
}

#ifndef LHAPDF
int pdfnamecmp (void) {
  if (!strcmp (pdfName[0], pdfName[1]))
    return 1;
  return 0;
}
#endif
