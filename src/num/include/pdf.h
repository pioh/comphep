/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/

#ifndef __PDF__
#define __PDF__

typedef struct pdfList {
  struct pdfList * next;
  char * name;                /* title name of distribution             */
  char * pathfile;            /* path to file where it is stored        */
  char * file;                /* name of file where it is stored        */
  int position;               /* ordering number of parton in the list of functions */
} pdfList;

typedef struct pdfStr {
  double mass;                /* mass of composite particle             */
  int number;
  int nq;                     /* number of points in Q-scale grid       */
  int nx;                     /* number of points in X-scale grid       */
  int lin;                    /* flag to detect linear extrapolation    */
  double * x_grid;            /* data for the X-grid                    */
  double * q_grid;            /* log() data for the Q-grid              */
  double x_min;               /* boundaries for x and q arguments       */
  double x_max;               /* boundaries for x and q arguments       */
  double q_min;               /* boundaries for x and q arguments       */
  double q_max;               /* boundaries for x and q arguments       */
  double * alpha;             /* data for QCD-alpha(Q) corresponding to */
                              /* the Q-grid points */
  double * strfun;            /* data for interpolation of distibution corresponding */
                              /* to  (Q-grid)*(X-grid) points;          */
                              /* index of (X-grid) is rinning first     */
  double pow0;                /* factor x^pow *(1-x)^pow1  must be      */
  double pow1;                /* applied after interpolation            */
  int approx;
  char * filename;            /* file name with pdf data                */
} pdfStr;

/*
  build a list od PDF str's
  items for the list are looked for
  in 3 places
  ./ 
  ../
  evn var $COMPHEP
*/
extern void comphepPdfList (char * p_name, pdfList ** list);

/*
  free memory allocated in 'list' 
*/
extern void delPdfList (pdfList * list);

/*
  read 'file' (in 'pathfile') and  fill pdf items 
  in the array data 'data' 
*/
extern int getPdfData (char * pathfile, char *file, int n_parton, pdfStr * data);

/*
  free  memory allocated for data items
  and assign  NULL to them.
*/
extern void freePdfData (pdfStr * data);

/*
  interpolates data for a given (x, q)
  according to information from pdf str W.
  result should be multiplied by the
  power factors later on
*/
extern double interFunc (double x, double q, pdfStr * W);

/*
???
  interpolates data for QCD-alpha(Q)
*/
extern double interAlpha (double q, pdfStr * W);

/*
  return QCD Lambda.
  In current pdf scheme implementation it is just dummy routine!
*/
double pdf_QCDLambda(void);

extern int CERNpdf_number (char *pdf, char *ver, int *PDFid, int *PDFgroup);
#endif
