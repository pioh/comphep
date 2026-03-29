/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/unix_utils.h"
#include "service2/include/files.h"
#include "service2/include/syst.h"
#include "service2/include/kfcodes.h"
#include "chep_crt/include/chep_crt.h"

#include "tag_reader.h"
#include "tag_parser.h"
#include "strfun_par.h"
#include "lhapdf.h"
#include "sf_lhapdf.h"

static int LHAPDFmem[2] = {0, 0};
static shortstr lhapdfName[2] = {"LHA:OFF", "LHA:OFF"};
static shortstr lhapdfBeam[2] = {"parton", "parton"};

int getLHAPDFset (int i) {
  return 0;
}

int getLHAPDFmem (int i) {
  return LHAPDFmem[i];
}

int p_lhapdf (char * p_name) {
  if (fabs(kfpart(p_name)) < 6 || kfpart(p_name) == 21) return 1;
  return 0;
}


void info_lhapdf (int i, Str_fun_Info * info) {
  char * name = strdup (lhapdfName[i] + 4);
  char * dot = strchr (name, '.');
  if (dot) *dot = 0; /* strip .LHpdf/.LHgrid extension if present */

  info->pdf_name[0] = 0;
  if (!strncmp ("H1", name, 2))      {strcpy (info->pdf_name, "H1");      strcpy (info->version, name + 2); }
  if (!strncmp ("GRV", name, 3))     {strcpy (info->pdf_name, "GRV");     strcpy (info->version, name + 3); }
  if (!strncmp ("cteq", name, 4))    {strcpy (info->pdf_name, "CTEQ");    strcpy (info->version, name + 4); }
  if (!strncmp ("MRST", name, 4))    {strcpy (info->pdf_name, "MRST");    strcpy (info->version, name + 4); }
  if (!strncmp ("MSTW", name, 4))    {strcpy (info->pdf_name, "MSTW");    strcpy (info->version, name + 4); }
  if (!strncmp ("a02m", name, 4))    {strcpy (info->pdf_name, "Alekhin"); strcpy (info->version, name + 4); }
  if (!strncmp ("ZEUS", name, 4))    {strcpy (info->pdf_name, "ZEUS");    strcpy (info->version, name + 4); }
  if (!strncmp ("Fermi", name, 5))   {strcpy (info->pdf_name, "Fermi");   strcpy (info->version, name + 5); }
  if (!strncmp ("Botje", name, 5))   {strcpy (info->pdf_name, "Botje");   strcpy (info->version, name + 5); }
  if (!strncmp ("NNPDF", name, 5))   {strcpy (info->pdf_name, "NNPDF");   strcpy (info->version, name + 5); }
  if (!strncmp ("Alekhin", name, 7)) {strcpy (info->pdf_name, "Alekhin"); strcpy (info->version, name + 7); }

  free (name);
  info->prt_mass = 0.938;
  strcpy (info->prt_name, lhapdfBeam[i]);

  info->N_extra_commands = 5;
  sprintf (info->extra_commands[0], "PDFid=%i", 0);
  sprintf (info->extra_commands[1], "PDFgr=%i", 0);
  sprintf (info->extra_commands[2], "LHAPDFid=%i", 0);
  sprintf (info->extra_commands[3], "LHAPDFmem=%i", LHAPDFmem[i]);
  info->PDFLIBset = 0;
  info->PDFLIBgroup = 0;
  info->LHAPDFset = 0;
  info->LHAPDFmember = LHAPDFmem[i];
  snprintf (info->extra_commands[4], strlen (lhapdfName[i]) - 4 + 10, "PDFfile=\'%s\'", lhapdfName[i] + 4);
  strcat (info->extra_commands[4], "\'");
}


void n_lhapdf (int i, char * beam, char * pdf) {
  if (get_sf_num(--i)) {
      sprintf (pdf, "%s:0:%i", lhapdfName[i], LHAPDFmem[i]);
  } else {
    strcpy (lhapdfBeam[i], "parton");
    strcpy (pdf, "LHA:OFF:0:0");
  }
  strcpy (beam, lhapdfBeam[i]);
}

static int 
extract_names (int i, char name[]) {
  int len;

  char * bname = strstr (name, "(");
  strcpy (lhapdfBeam[i], bname + 1);
  len = strlen (lhapdfBeam[i]);
  lhapdfBeam[i][len - 1] = '\0';
  trim (lhapdfBeam[i]);

  len = strlen (name) - len - 1;
  name[len] = '\0';
  strcpy (lhapdfName[i], name);
  trim (lhapdfName[i]);

  return 1;
}

int r_lhapdf (int i, char * name) {
  --i;
  if (strstr (name, "LHA:")) {
    int set, mem;
    char name1[256];
    extract_names (i, name);
    if (get_sf_info (lhapdfName[i], name1, &set, &mem)) {
      strcpy (lhapdfName[i], name1);
      LHAPDFmem[i] = mem;
    }
    return 1;
  }
  return 0;
}

int be_lhapdf (int i, double * be, double * mass, char * p_name) {
  int pdfcode;
  int undefined = 1;
  lhapdfList * list = NULL;
  lhapdfList * list_;

  pdfcode = kfpart(p_name);
  comphepLhapdfList (&list);
  list_ = list;

  --i;
  while (list_) {
    if (0 == strcmp (list_->name, lhapdfName[i]+4) && LHAPDFmem[i] == list_->mem) {
      undefined = 0;
      initLHAPDF (i, list_->name, list_->mem, pdfcode);
      set_alphaMode(i + 1);
    }
    list_ = list_->next;
  }
  if (list) {
    delLhapdfList (list);
  }
  if (undefined) {
    strcpy (lhapdfName[i], "OFF");
  }

  * be = 1.;
  * mass = 1.;
  return 1;
}


int pdf_lhapdf (int i, char * p_name) {
  int k = 0;
  int WIDTH = 42;
  void * pscr = NULL;
  char * menustring = NULL;
  lhapdfList *list = NULL;
  lhapdfList *list_ = NULL;

  comphepLhapdfList (&list_);
  if (list_){
    list = list_;
  } else {
    return 0;
  }
  for (list_ = list; list_; list_ = list_->next) ++k;

  menustring = malloc (2 + WIDTH * k);
  menustring[0] = WIDTH;
  menustring[1] = 0;
  for (list_ = list; list_; list_ = list_->next) {
    char display_name[31];
    if (list_->mem > 99999) {
      fprintf (stderr, "Warning! Too long PDF member number... Should be not longer then 5 symbols! Now it is %i\n", list_->mem);
      list_->mem = 99999;
    }
    strncpy (display_name, list_->name, 30);
    display_name[30] = 0;
    sprintf (menustring + strlen (menustring), " LHA:%-*s  %*i", WIDTH - 12, display_name, 5, list_->mem);
  }

  k = 0;
  menu1 (35, 8, "", menustring, "n_sf_lha", &pscr, &k);

  if (k) {
    int l = k - 1;
    for (list_ = list; l; --l) list_ = list_->next;
    --i;
    sprintf (lhapdfName[i], "LHA:%s", list_->name);
    LHAPDFmem[i] = list_->mem;
    if (!strcmp("parton", lhapdfBeam[i])) strcpy (lhapdfBeam[i], "proton");
  }
  free (menustring);
  put_text (&pscr);
  delLhapdfList (list);
  return k;
}


int beam_lhapdf (int i) {
  int key = 0;
  void *pscr = NULL;
  shortstr previous;
  char strmen[] = "*" 
                  " parton      "
                  " proton      "
                  " anti-proton ";

  strmen[0] = (strlen (strmen) -1) / 3;
  --i;
  strcpy(previous, lhapdfBeam[i]);
  menu1 (64, 7, "", strmen, "n_beam_*", &pscr, &key);
  if (key) {
    switch (key) {
      case 1:
        strcpy (lhapdfName[i], "OFF");
        strcpy (lhapdfBeam[i], "parton");
        break;
      case 2:
        strcpy (lhapdfBeam[i], "proton");
        break;
      case 3:
        strcpy (lhapdfBeam[i], "anti-proton");
        break;
    }
  }
  if (0 == key && !strcmp("parton", lhapdfBeam[i])) {
    strcpy (lhapdfName[i], "OFF");
    LHAPDFmem[i] = 0;
    set_sf_num(i, 0);
  }
  if (!strcmp("proton", lhapdfBeam[i]) && !strcmp("parton", previous)) {
    strcpy (lhapdfName[i], "LHA:undefined");
    LHAPDFmem[i] = 0;
  }
  return 0;
}


double c_lhapdf (int i, double x, double q) {
  double r = lhapdfVal (x, q, i - 1);
  if (r < 0.0) r = 0.0;
  return r;
}


double alpha_lhapdf (double q) {
  return lhapdf_interAlpha (q);
}

int lhapdf_namecmp (void) {
  if (!strcmp (lhapdfName[0], lhapdfName[1]))
    return 1;
  return 0;
}
