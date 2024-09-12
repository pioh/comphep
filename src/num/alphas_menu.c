/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Slava Ilyin 
* ----------------------------------------------------
*/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "service2/include/chep_limits.h"
#include "service2/include/read_func.h"
#include "service2/include/parser.h"
#include "service2/include/syst.h"
#include "chep_crt/include/crt_util.h"
#include "plot/include/plot.h"

#include "out_ext.h"

#ifdef LHAPDF
  #include "clhapdf.h"
  #include "lhapdf.h"
  #include "sf_lhapdf.h"
#else
  #include "pdf.h"
  #include "sf_pdf.h"
#endif

#include "strfun_par.h"
#include "phys_val.h"
#include "alphas2.h"

#define  XMB 4.3
#define  XMC 1.3
#define  XMTOP 175.

static shortstr scale_str = "91.187";

char * get_scale_form (void) 
{
  return scale_str;
}

/********************************************************************
*   This subroutine is called during phase space integration        *
*     (during BASES and SPRING run)                                 *
*     for each phase space point after calculation particles        *
*     momenta and just before calculation of structure functions    *
*     and squared matrix element.                                   *
*                                                                   *
*   It can be used, for example, to provide QCD running strong      *
*     coupling constant. In the program there is used the notation  *
*                          GG                                       *
*     for QCD coupling constant.                                    *
********************************************************************/
void 
alf_ (double q)
{
  int k;
  char xname[7];
/*... search for QCD coupling constant among process parameters */
  for (k = 1; k <= nvar_; k++)
    {
      vinf_ (k, xname, NULL);
      if (strcmp ("GG", xname) == 0)
        {
          asgn_ (k, sqrt (4 * M_PI * alpha_2 (q)));
          return;
        }
    }
}

/***************************************************************************
*  Flavour number matching inALPHAS2:                                      *
*                                                                          *
*        NF              6       5      4        3        2                *
*    Lambda_QCD (GeV)    0.0683  0.135  0.1985   0.2408   0.2293           *
*                        0.1052  0.200  0.28275  0.33078  failed           *
*                                                                          *
*   Comment: for NF=2 and Lambda_QCD(5)=QCDLF(5)=0.2GeV matching is failed *
* because of the term with LOG(LOG(Ms/QCDLF(3))) where Ms<QCDLF(3).        *
***************************************************************************/
static int 
qq0 (char *s, double *v)
{
  char key, plist[20];

  if (isdigit (*s))
    {
      sscanf (s, "%lf", v);
      return 1;
    }
  if (checkPhysVal (s, &key, plist))
    {
      *v = 100;
      return 1;
    }
  else
    return 0;
}

int 
qcdmen_ (void)
{
  void *pscr = NULL;
  int mode;
  int returnCode = 0;
  static double lam = 0.1185;
  vshortstr name[2];

  pinf_ (1, 1, name[0], NULL);
  pinf_ (1, 2, name[1], NULL);
  initStrFun (name[0], name[1]);
L10:
  {
    char strmen[] = "\032"
    " QCD Lambda6= XXX         "
    " Q(GeV) = YYY             "
    " Alpha(Q) plot            ";

    if (get_alphaMode())
      improveStr (strmen, "XXX", "???");
    else
    {
      improveStr (strmen, "XXX", "%-.3lgGeV", lam);
      recalc_alphas ();
      setLambda6 (lam);
    }
    improveStr (strmen, "YYY", "%-.16s", scale_str);
    menu1 (52, 8, "QCD alpha", strmen, "n_alpha", &pscr, &mode);
  }
  switch (mode)
    {
    case 0:
      return returnCode;
    case 1:
      if (get_alphaMode())
        messanykey (20, 20, " alpha(s) is confined "
                    "by PDF structure function");
      else
        {
          goto_xy (50, 12);
          if (correctDouble (3, 15, "Enter new value ", &lam, 1))
            recalc_alphas ();
            setLambda6 (lam);
            returnCode = 1;
        }
      break;
    case 2:
      {
        int npos = 1, rc;
        do
          {
            double dscale;
            goto_xy (2, 12);
            print ("QCD scale: ");
            if (str_redact (scale_str, npos, 60) == KB_ENTER)
              returnCode = 1;
            goto_xy (2, 12);
            clr_eol ();
            rc = calcExpression (scale_str, qq0, &dscale);
            if (rc && rc != cannotevaluate)
              {
                npos = rderrpos;
                if (rc == unknownidentifier)
                  warnanykey (10, 10, " Unknown parameter");
                else
                  warnanykey (10, 10, " Syntax error");
              }
          }
        while (rc && rc != cannotevaluate);
      }
      break;
    case 3:
      {
        void *screen;
        int i;
        double f[150];
        static double qMin = 1, qMax = 100;
        static int nPoints = 100;

        get_text (1, 1, maxCol (), maxRow (), &screen);

        if (correctDouble (40, 15, "Q_min=", &qMin, 0) && qMin >= 0.5
            && correctDouble (40, 16, "Q_max=", &qMax, 0) && qMax > qMin
            && correctInt (33, 17, "number of points=", &nPoints, 0)
            && nPoints > 3 && nPoints <= 150)
          {
            for (i = 0; i < nPoints; i++)
              {
                double theq = qMin + i * (qMax - qMin) / (nPoints - 1);
                f[i] = alpha_2 (theq);
/*              fprintf (stderr, "   f[%i] = %.5E\n", i, f[i]);*/
                fprintf (stdout, " %.5E       %.5E\n", theq, f[i]);
              }
            plot_table (qMin, qMax, nPoints, f, NULL, " ", "Q [GeV]", "Alpha(Q)");
          }
        else
          messanykey (40, 18,
                      " Correct input is \n"
                      " 0.5<= Q_min <Q_max\n"
                      " number of points <=150 and >=4");
        put_text (&screen);
      }

    }
  goto L10;
}

int WriteQCDInfo (FILE * mode) {
  shortstr format;
  double lambda = QCDLambda ();
  strcpy (format, "Lambda6 = \%E  Scale = \%s");
  if (get_alphaMode()) {
#ifdef LHAPDF
    strcpy (format, "Lambda5 = \%E  Scale = \%s");
    lambda = lhapdf_QCDLambda ();
#else
    lambda = pdf_QCDLambda ();
#endif
  }
  fprintf (mode, format, lambda, scale_str);
  return 0;
}

void InitQCDInfo (void) {
  strcpy (scale_str, "91.187");
}

int ReadQCDInfo (FILE * mode) {
  midstr aninf;
  midstr format;
  double lam = 0.;
  fgets (aninf, 1024, mode);
#ifdef LHAPDF
  strcpy (format, " Lambda5 = \%lf  Scale = \%s");
  if (aninf[6] == '6') {
    if (get_alphaMode ())
      fprintf (stderr, "\nn_comphep (warning): the session.dat has been created for \n"
                       "the internal QCD alpha or for PDF in Original CompHEP Format,\n"
                       " QCD Lambda is for Nflavours = 6!!!\n");
    strcpy (format, " Lambda6 = \%lf  Scale = \%s");
  }
#else
  strcpy (format, " Lambda6 = \%lf  Scale = \%s");
  if (aninf[6] == '5') {
    fprintf (stderr, "\nn_comphep (warning): the session.dat has been created for LHAPDF, QCD Lambda is for Nflavours = 5!\n");
    strcpy (format, " Lambda5 = \%lf  Scale = \%s");
  }
#endif
  sscanf (aninf, format, &lam, scale_str);
  recalc_alphas ();
  setLambda6 (lam);
  return 0;
}

/************************************************************
*  This subroutine user can call to evaluate the transfered *
*    momentum scale of the process - could be necessary for *
*    structure function interface.                          *
*                                                           *
*  In the example from STRFUN.F this subroutine is called   *
*    for the evaluation of transfered momentum scale.       *
*  NSCALE = 0 for t-type, and =1 for s-type of              *
*    transfered momentum scale.                             *
************************************************************/
static int 
qq1 (char *s, double *v)
{
  char key;
  char plist[20];

  if (isdigit (*s))
    {
      sscanf (s, "%lf", v);
      return 1;
    }
  checkPhysVal (s, &key, plist);
  *v = calcPhysVal (key, plist, "");
  return 1;
}

double qcd_Scale_chep (void) {
  double dscale;
  if (calcExpression (scale_str, qq1, &dscale))
    {
      fprintf (stderr, " ERROR in evaluation of  QCD scale\n");
      exit (0);
    }

  if (dscale < .3f)
    dscale = .3f;
  return dscale;
}
