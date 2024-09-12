/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* with changes by V.Ilyin for #-mdl, 6 Sept 2000, 8 Oct 2000.
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/f_c.h"
#include "service2/include/syst.h"
#include "chep_crt/include/crt_util.h"
#include "out_ext.h"

#include "subproc.h"
#include "sf_epa.h"
#include "sf_isr.h"
#include "sf_lsr.h"
#include "sf_pdf.h"
#include "sf_lhapdf.h"
#include "strfun_par.h"
#include "strfun.h"

#define MAXFUN 4
#define FUNLEN 40

static int uprt[2];
static int dprt[2];
static int Uprt[2];
static int Dprt[2];
static int beam[2];

static struct {
  int (*myParticle) (char *);
  void (*fullName) (int, char *, char *);
  void (*realSTRFUN_info) (int, Str_fun_Info *);
  int (*readName) (int, char *);
  int (*beam_menu) (int);
  int (*menu) (int, char *);
  int (*init) (int, double *, double *, char *);
  double (*val) (int, double, double);
} strFun[MAXFUN] =
{
  {
    p_epa__, n_epa__, info_epa__, r_epa__, b_epa__, m_epa__, i_epa__, c_epa__
  }
  ,
  {
    p_lsr__, n_lsr__, info_lsr__, r_lsr__, b_lsr__, m_lsr__, i_lsr__, c_lsr__
  }
  ,
  {
    p_isr__, n_isr__, info_isr__, r_isr__, b_isr__, m_isr__, i_isr__, c_isr__
  }
  ,
#ifdef LHAPDF
  {
    p_lhapdf, n_lhapdf, info_lhapdf, r_lhapdf, beam_lhapdf, pdf_lhapdf, be_lhapdf, c_lhapdf
  }
#else
  {
    p_pdf, n_pdf, info_pdf, r_pdf, beam_pdf, pdf_pdf, be_pdf, c_pdf
  }
#endif
};

static double strfun_f (int i, double x, double q) {
  int sf = get_sf_num(i - 1);
  switch (beam[i - 1]) {
    case 0:
      return strFun[sf - 1].val (i, x, q);   /* standard parton */
      break;
    case 1:                     /* f_u(x)+f_c(x) */
      return strFun[sf - 1].val (uprt[i - 1], x, q) + strFun[sf - 1].val (3, x, q);
      break;
    case 2:                     /* f_U(x)+f_C(x) */
      return strFun[sf - 1].val (Uprt[i - 1], x, q) + strFun[sf - 1].val (3, x, q);
      break;
    case 3:                     /* f_d(x)+f_s(x) */
      return strFun[sf - 1].val (dprt[i - 1], x, q) + strFun[sf - 1].val (4, x, q);
      break;
    case 4:                     /* f_U(x)+f_S(x) */
      return strFun[sf - 1].val (Dprt[i - 1], x, q) + strFun[sf - 1].val (4, x, q);
      break;
    case 5:                     /* f_u(x)+f_d(x)+f_s(x)+f_c(x) */
      return
          strFun[sf - 1].val (3, x, q)
        + strFun[sf - 1].val (4, x, q)
        + strFun[sf - 1].val (uprt[i - 1], x, q)
        + strFun[sf - 1].val (dprt[i - 1], x, q);
      break;
    case 6:                     /* f_U(x)+f_D(x)+f_S(x)+f_C(x) */
      return
          strFun[sf - 1].val (3, x, q)
        + strFun[sf - 1].val (4, x, q)
        + strFun[sf - 1].val (Uprt[i - 1], x, q)
        + strFun[sf - 1].val (Dprt[i - 1], x, q);
      break;
  }
  fprintf (stderr, "***Error! function: strfun_f, unknown particle!\n");
  return 0.0;
}

static double strfun_nonf (double x, double y, double q) {
  double sinc2 = 0.221 * 0.221; /* sin of Cabbibo angle squared */
  double cosc2 = 1. - sinc2;
  int sf0 = get_sf_num (0);
  int sf1 = get_sf_num (1);

  switch (beam[1] + 10 * beam[0]) {
    case 0:
      return strFun[sf0 - 1].val (1, x, q) * strFun[sf0 - 1].val (2, y, q);         /* standard parton */
      break;
    case 12:                    /* f_u(x)*f_U(y)+f_c(x)*f_C(y) */
      return
        strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (Uprt[1], y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q);
      break;
    case 21:                    /* f_U(x)*f_u(y)+f_C(x)*f_c(y) */
      return
        strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (uprt[1], y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q);
      break;
    case 34:                    /* f_d(x)*f_D(y)+f_s(x)*f_S(y) */
      return
        strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (Dprt[1], y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q);
      break;
    case 43:                    /* f_D(x)*f_d(y)+f_S(x)*f_s(y) */
      return
        strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (dprt[1], y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q);
      break;
    case 13:                    /* f_u(x)*f_d(y)*cosc2+f_u(x)*f_s(y)*sinc2+f_c(x)*f_d(y)*sinc2+f_c(x)*f_s(y)*cosc2 */
      return
        strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (dprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (4, y, q) * sinc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (dprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (4, y, q) * cosc2;
      break;
    case 31:                    /* f_d(x)*f_u(y)*cosc2+f_s(x)*f_u(y)*sinc2+f_d(x)*f_c(y)*sinc2+f_s(x)*f_c(y)*cosc2 */
      return
        strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (uprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (uprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (3, y, q) * sinc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (3, y, q) * cosc2;
      break;
    case 14:                    /* f_u(x)*f_D(y)*cosc2+f_u(x)*f_S(y)*sinc2+f_c(x)*f_D(y)*sinc2+f_c(x)*f_S(y)*cosc2 */
      return
        strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (Dprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (4, y, q) * sinc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (Dprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (4, y, q) * cosc2;
      break;
    case 41:                    /* f_D(x)*f_u(y)*cosc2+f_S(x)*f_u(y)*sinc2+f_D(x)*f_c(y)*sinc2+f_S(x)*f_c(y)*cosc2 */
      return
        strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (uprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (uprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (3, y, q) * sinc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (3, y, q) * cosc2;
      break;
    case 23:                    /* f_U(x)*f_d(y)*cosc2+f_C(x)*f_d(y)*sinc2+f_U(x)*f_s(y)*sinc2+f_C(x)*f_s(y)*cosc2 */
      return
        strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (dprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (dprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (4, y, q) * sinc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (4, y, q) * cosc2;
      break;
    case 32:                    /* f_d(x)*f_U(y)*cosc2+f_d(x)*f_C(y)*sinc2+f_s(x)*f_U(y)*sinc2+f_s(x)*f_C(y)*cosc2 */
      return
        strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (Uprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (3, y, q) * sinc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (Uprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (3, y, q) * cosc2;
      break;
    case 24:                    /* f_U(x)*f_D(y)*cosc2+f_C(x)*f_D(y)*sinc2+f_U(x)*f_S(y)*sinc2+f_C(x)*f_S(y)*cosc2 */
      return
        strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (Dprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (Dprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (4, y, q) * sinc2
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (4, y, q) * cosc2;
      break;
    case 42:                    /* f_D(x)*f_U(y)*cosc2+f_D(x)*f_C(y)*sinc2+f_S(x)*f_U(y)*sinc2+f_S(x)*f_C(y)*cosc2 */
      return
        strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (Uprt[1], y, q) * cosc2
        + strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (3, y, q) * sinc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (Uprt[1], y, q) * sinc2
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (3, y, q) * cosc2;
      break;
    case 11:                    /* f_u(x)*f_u(y)+f_c(x)*f_c(y) */
      return
        strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (uprt[1], y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q);
      break;
    case 22:                    /* f_U(x)*f_U(y)+f_C(x)*f_c(y) */
      return
        strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (Uprt[1], y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q);
      break;
    case 33:                    /* f_d(x)*f_d(y)+f_s(x)*f_s(y) */
      return
        strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (dprt[1], y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q);
      break;
    case 44:                    /* f_D(x)*f_D(y)+f_S(x)*f_S(y) */
      return
        strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (Dprt[1], y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q);
      break;
    case 56:                    /*f_u(x)*f_U(y)+f_d(x)*f_D(y)+f_c(x)*f_C(y)+f_s(x)*f_S(y) */
      return
        strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (Uprt[1], y, q)
        + strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (Dprt[1], y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q);
      break;
    case 65:                    /*f_U(x)*f_u(y)+f_D(x)*f_d(y)+f_C(x)*f_c(y)+f_S(x)*f_s(y) */
      return
        strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (uprt[1], y, q)
        + strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (dprt[1], y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q);
      break;
    case 55:                    /* f_d(x)*f_d(y)+f_u(x)*f_u(y)+f_s(x)*f_s(y)+f_c(x)*f_c(y) */
      return
        strFun[sf0 - 1].val (uprt[0], x, q) * strFun[sf1 - 1].val (uprt[1], y, q)
        + strFun[sf0 - 1].val (dprt[0], x, q) * strFun[sf1 - 1].val (dprt[1], y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q);
      break;
    case 66:                    /* f_D(x)*f_D(y)+f_U(x)*f_U(y)+f_S(x)*f_S(y)+f_C(x)*f_C(y) */
      return
        strFun[sf0 - 1].val (Uprt[0], x, q) * strFun[sf1 - 1].val (Uprt[1], y, q)
        + strFun[sf0 - 1].val (Dprt[0], x, q) * strFun[sf1 - 1].val (Dprt[1], y, q)
        + strFun[sf0 - 1].val (3, x, q) * strFun[sf1 - 1].val (3, y, q)
        + strFun[sf0 - 1].val (4, x, q) * strFun[sf1 - 1].val (4, y, q);
      break;
  }
  fprintf (stderr, "***Error! function: strfun_nonf, unknown particle!\n");
  return 0.0;
}                               /* strfun_nonf */


double strfun_ (int factr, double x, double y, double q) {
  double xstr, xstr1, xstr2;
  int sf0 = get_sf_num (0);
  int sf1 = get_sf_num (1);

  if (factr) {
    if (sf0) {
      xstr1 = strfun_f (1, x, q);
    } else {
      xstr1 = 1.0;
    }
    if (sf1) {
      xstr2 = strfun_f (2, y, q);
    } else {
      xstr2 = 1.0;
    }
    xstr = xstr1 * xstr2;
  } else {
    if (!sf0 && !sf1) {
      xstr = 1.0;
    } else if (sf0 && sf1) {
        xstr = strfun_nonf (x, y, q);
      } else {
       fprintf (stderr, "***** strfun_: In non-factorize case both parton must have PDF\n");
       xstr = 0.0;
      }
  }
  if (!finite (xstr)) {
    fprintf (stderr, "***Error! x= %18.15g, y=%18.15g, q=%18.15g, xstr=%18.15g\n", x, y, q, xstr);
  }
  return xstr;
}

int initStrFun (char p_name1[], char p_name2[]) {
  int l;
  int i;
  int returnCode = 0;
  vshortstr name[2];

  strcpy(name[0], p_name1);
  strcpy(name[1], p_name2);
  set_alphaMode (0);
  for (i = 0; i < 2; i++) {
    l = get_sf_num (i);
    set_sf_be (i, 0.);
    if (l) {
      double mass;
      double be;
      l--;
      if (!strFun[l].myParticle (name[i]) || !strFun[l].init (i + 1, &be, &mass, name[i])) {
        set_sf_num (i, 0);
        fprintf (stderr, "\n%i-th Stucture function is switched OFF.\n", i + 1);
        returnCode = 2;
      } else {
        set_sf_mass (i, mass);
        set_sf_be (i, be);
        beam[i] = -1;
        if (!strchr (name[i], 35))
          beam[i] = 0;
        else if (!strcmp (name[i], "u#"))
          beam[i] = 1;
        else if (!strcmp (name[i], "U#"))
          beam[i] = 2;
        else if (!strcmp (name[i], "d#"))
          beam[i] = 3;
        else if (!strcmp (name[i], "D#"))
          beam[i] = 4;
        else if (!strcmp (name[i], "q#"))
          beam[i] = 5;
        else if (!strcmp (name[i], "Q#"))
          beam[i] = 6;
        else {
          fprintf (stderr, "***Error! function: initStrFun, unknown initial particle = %s!\n", name[i]);
          exit (2);
        }

        dprt[0] = 9;
        uprt[0] = 8;
        Uprt[0] = 6;
        Dprt[0] = 5;
        if (1 == i) {
          dprt[1] = 9;
          uprt[1] = 8;
          Uprt[1] = 6;
          Dprt[1] = 5;
          if (beam[0] > 0 && beam[1] > 0) {
            if (!pdfnamecmp ()) {
              dprt[1] = 5;
              uprt[1] = 6;
              Uprt[1] = 8;
              Dprt[1] = 9;
            }
          }
        }
      }
    }
  }

  return returnCode;
}

void strFunName (int i, char * pbeam, char * pdf) {
  int num = get_sf_num (i - 1);
  if (num) {
    strFun[num - 1].fullName (i, pbeam, pdf);
  } else {
    strcpy (pdf, "OFF");
    strcpy (pbeam, "parton");
  }

  if (!strcmp (pbeam, "parton") && !strcmp (pdf, "OFF")) {
    strcpy (pdf, "OFF");
    strcpy (pbeam, "parton");
    set_sf_num (i - 1, 0);
  }
}

int beam_menu (int i) {
  int l;
  vshortstr part_name;
  int num = get_sf_num (i - 1);
  if (num) {
    strFun[num - 1].beam_menu (i);
  } else {
    pinf_ (proces_1.nsub, i, part_name, NULL);
    for (l = 0; l < MAXFUN; l++) {
      if (strFun[l].myParticle (part_name)) {
        break;
      }
    }
    set_sf_num (i - 1, l + 1);
    strFun[l].beam_menu (i);
  }
  return 1;
}

int pdf_menu (int ibeam) {
  int k;
  int l;
  int mode;
  int nfun[MAXFUN];
  vshortstr part_name;
  longstr strmen;
  shortstr name;
  void * pscr = NULL;

  strmen[0] = FUNLEN + 1;
  pinf_ (proces_1.nsub, ibeam, part_name, NULL);
  sprintf (strmen + 1, " %-*.*s", FUNLEN, FUNLEN, "OFF");

  k = 0;
  for (l = 0; l < MAXFUN; l++) {
    int ju;
    ju = strFun[l].myParticle (part_name);
    if (ju) {
      shortstr pbeam;
      shortstr pdf;
      nfun[k++] = l;
      strFun[l].fullName (ibeam, pbeam, pdf);
      sprintf (strmen + 1 + (FUNLEN + 1) * k, " %-*.*s", FUNLEN, FUNLEN, pdf);
    }
  }

  if (!k) {
    warnanykey (15, 15, "PDF for this particle\nare not known\n");
    return 0;
  }

  if (strncmp (name, "PDF:", 4)) {       /* Not parton */
    menu1 (77 - FUNLEN, 7, "", strmen, "n_strfun", &pscr, &mode);
    if (mode == 0) {
      return 0;
    }
    if (mode == 1) {
      set_sf_num (ibeam - 1, 0);
    } else {
      strFun[nfun[mode - 2]].menu (ibeam, part_name);
      set_sf_num (ibeam - 1, nfun[mode - 2] + 1);
    }
    put_text (&pscr);
    return 1;
  }

  mode = strFun[nfun[0]].menu (ibeam, part_name);
  switch (mode) {
    case 0:
      return 0;
    case 1:
      set_sf_num (ibeam - 1, 0);
      break;
    default:
      set_sf_num (ibeam - 1, nfun[0] + 1);
  }
  return 1;
}

int rd_sf__ (FILE * mode) {
  midstr sf_txt[2];
  vshortstr part_name[2];
  int l, i;

  for (i = 0; i < 2; i++) {
    set_sf_num (i, 0);
  }

  if (2 != fscanf (mode, "  StrFun1: %[^\n]\n  StrFun2: %[^\n]\n", sf_txt[0], sf_txt[1])) {
    return 1;
  }

  for (i = 0; i < 2; i++)    {
    pinf_ (proces_1.nsub, i + 1, part_name[i], NULL);
    for (l = 0; l < MAXFUN; l++) {
      if (strFun[l].myParticle (part_name[i]) && strFun[l].readName (i + 1, sf_txt[i])) {
        set_sf_num (i, l + 1);
        break;
      }
    }
  }
  initStrFun (part_name[0], part_name[1]);
  return 0;
}

int wrt_sf__ (FILE * mode) {
  char sf_txt[2][STRSIZ];
  int i;

  for (i = 0; i < 2; i++) {
    int num = get_sf_num(i);
    fprintf (mode, "  StrFun%d: ", i + 1);
    if (num) {
      shortstr pbeam;
      shortstr pdf;
      strFun[num - 1].fullName (i + 1, pbeam, pdf);
      sprintf (sf_txt[i], "%s(%s)", pdf, pbeam);
    } else {
      strcpy (sf_txt[i], "OFF(parton)");
    }
    fprintf (mode, "%s\n", sf_txt[i]);
  }
  return 0;
}

void wrt_sf_NF_ (int i, Str_fun_Info * info) {
  int num = get_sf_num (i);
  if (num) {
    strFun[num - 1].realSTRFUN_info (i, info);
  } else {
    double part_mass;
    vshortstr part_name;

    strcpy (info->pdf_name, "OFF");
    pinf_ (proces_1.nsub, i + 1, part_name, &part_mass);
    strcpy (info->prt_name, part_name);
    info->prt_mass = part_mass;
    info->N_extra_commands = 0;
    snprintf (info->version, 1, " ");
  }
}
