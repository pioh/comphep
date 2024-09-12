/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Slava Ilyin
* ------------------------------------------------------
*/
#include <unistd.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/files.h"
#include "service2/include/tptcmac.h"
#include "service2/include/unix_utils.h"
#include "service2/include/4_vector.h"
#include "chep_crt/include/crt_util.h"

#include "out_ext.h"

#include "cut.h"
#include "q_kin.h"
#include "regul.h"
#include "rw_sess.h"
#include "subproc.h"
#include "vegas.h"
#include "strfun.h"
#include "alphas_menu.h"
#include "histogram.h"
#include "evnt_menu.h"
#include "core_data.h"
#include "num_serv.h"
#include "runVegas.h"
#include "userFactor.h"

extern double computer_eps;

static int userFactor_key = 0;

static vegasGrid *veg_Ptr = NULL;

int
ClearVegasGrid (void)
{
  vegas_finish (veg_Ptr);
  veg_Ptr = NULL;

  return 1;
}


int
WriteVegasGrid (FILE * f)
{
  if (veg_Ptr)
    {
      int i, j;
      double *x = veg_Ptr->x_grid;
      fprintf (f, " Vegas_grid: dim=%d  size=%d\n", veg_Ptr->ndim,
               veg_Ptr->ndmx);
      for (i = 0; i < veg_Ptr->ndim; i++)
        {
          for (j = 0; j <= veg_Ptr->ndmx; j++)
            fprintf (f, " %.15E", *(x++));
          fprintf (f, "\n");
        }
    }
  else
    {
      fprintf (f, " Vegas_grid: dim=%d  size=%d\n", 0, 0);
    }
  return 0;
}

int
ReadVegasGrid (FILE * f)
{
  int i, j;
  int ndim, ndmx;
  double *x;

  if (veg_Ptr)
    {
      vegas_finish (veg_Ptr);
      veg_Ptr = NULL;
    }
  fscanf (f, " Vegas_grid: dim=%d  size=%d\n", &ndim, &ndmx);
  if (ndim && ndmx)
    {
      veg_Ptr = vegas_init (ndim, ndmx);
      x = veg_Ptr->x_grid;
      for (i = 0; i < ndim; i++) {
        for (j = 0; j <= ndmx; j++) {
          fscanf (f, " %lf", (x++));
        }
      }
    }
  return 0;
}

static int nCall;
static double badPoints;
static double negPoints;
static int hFill = 0;

static double
me2_func (double *x, double wgt)
{
  int ntot = nin_ + nout_;
  double ret_val = 0.;
  double ytilda;
  /* preparation of kinematics for scalar products */
  double factor_0 = mkmom (x, &ytilda);

  nCall++;
  if (factor_0) {
    lorrot (get_rapidity() + ytilda, ntot);
    factor_0 *= calcCutFactor() * userfactorFun(userFactor_key);
    if (factor_0) {
      int err = 0;
  /* call for 'running strong coupling constant' */
      alf_ (qcd_Scale_chep ());
 /*     printf("\n nsub in sqme = %d \n",proces_1.nsub); */
      ret_val = factor_0 * sqme_ (proces_1.nsub, pvect, &err);
      if (err) {
        badPoints += (ret_val > 0 ? ret_val * wgt : -ret_val * wgt);
      }
      if (fabs (ret_val) < computer_eps) {
        negPoints += ret_val * wgt;
      }
    }
  }

  if (hFill) {
    fillHists (ret_val * wgt);
  }

  return ret_val;
}


static int gui_vegas (FILE * prtfile, int nLine) {
  int i;
  int nln = nLine;
  double sd;
  double avgi;
  vegas_integral in = get_vegas_integral ();
  mcintr mc = get_mc_info ();

  if (in.n_it == 0) {
    init_vegas_integral ();
    in = get_vegas_integral ();
  }

  for (i = 0; i < mc.itmx0; ++i) {
    char errtxt[100] = "";

    nCall = 0;
    negPoints = 0;
    badPoints = 0;
    hFill = 1;
    if (vegas_int (veg_Ptr, mc.ncall0, proces_1.nsub, mc.itmx0, i + 1, me2_func, &avgi, &sd)) {
      break;
    }
    in.old = 1;
    negPoints /= nCall;
    badPoints /= nCall;
    in.nCallTot += nCall;
    scrcolor (FGmain, BGmain);
    ++in.n_it;
    printLn (prtfile, &nln, "%3d    %12.4E     %10.2E %8d ", in.n_it, avgi, 100 * sd / fabs (avgi), nCall);
    if (negPoints < 0) {
      sprintf (errtxt, " Negative points %.1G%%;", -100 * negPoints / (avgi - 2 * negPoints));
    }
    if (badPoints) {
      sprintf (errtxt + strlen (errtxt), "Bad Precision %.1G%%;", 100 * badPoints / (avgi - 2 * negPoints));
    }
    if (errtxt[0]) {
      scrcolor (Red, BGmain);
      printLn (prtfile, &nln, "%s", errtxt);
    }
    sd = 1 / (sd * sd);
    in.s0 += sd;
    in.s1 += avgi * sd;
    in.s2 += avgi * avgi * sd;
  }

  if (in.n_it > 1) {
    double ksi_test = (in.s2 - in.s1 * in.s1 / in.s0) / (in.n_it - 1);
    double cs = in.s1 / in.s0;
    double cse = 100. * sqrt (in.s0) / fabs (in.s1);
    scrcolor (FGmain, BGmain);
    printLn (prtfile, &nln, " < >   %12.4E     %10.2E %8d %8.1G",
    cs, cse, in.nCallTot, ksi_test);
  }
  set_vegas_integral (in);
  return nln;
}


static int gui_vegasAll (FILE * prtfile, int nLine) {
  int i;
  int nln = nLine;
  double sd;
  double avgi;

/*---------------------------------*/
  int ndim = imkmom ();
/*  int nses = get_nsession (); */
/*---------------------------------*/

  vegas_integral in = get_vegas_integral ();
  mcintr mc = get_mc_info ();



/*  if (in.n_it == 0) */
 {
    init_vegas_integral ();
    in = get_vegas_integral ();
  }

/*----------------------------------------------------------------------------*/
   vegas_finish (veg_Ptr);
   veg_Ptr = NULL;
   veg_Ptr = vegas_init (ndim, MAX_NDMX);
   correctHistList (0);
   mc.ncall0 = 2 * generateVegasCubs (veg_Ptr, mc.ncall0 / 2.);
   set_mc_info (mc);
/*----------------------------------------------------------------*/


  for (i = 0; i < mc.itmx0; ++i) {
    char errtxt[100] = "";

    nCall = 0;
    negPoints = 0;
    badPoints = 0;
    hFill = 1;
    if (vegas_int (veg_Ptr, mc.ncall0, proces_1.nsub, mc.itmx0, i + 1, me2_func, &avgi, &sd)) {
      break;
    }

    in.old = 1;  

    negPoints /= nCall;
    badPoints /= nCall;
    in.nCallTot += nCall;
    scrcolor (FGmain, BGmain);

    ++in.n_it; 

    printLn (prtfile, &nln, "%3d    %12.4E     %10.2E %8d ", in.n_it, avgi, 100 * sd / fabs (avgi), nCall);
    if (negPoints < 0) {
      sprintf (errtxt, " Negative points %.1G%%;", -100 * negPoints / (avgi - 2 * negPoints));
    }
    if (badPoints) {
      sprintf (errtxt + strlen (errtxt), "Bad Precision %.1G%%;", 100 * badPoints / (avgi - 2 * negPoints));
    }
    if (errtxt[0]) {
      scrcolor (Red, BGmain);
      printLn (prtfile, &nln, "%s", errtxt);
    }
    sd = 1 / (sd * sd);
    in.s0 += sd;
    in.s1 += avgi * sd;
    in.s2 += avgi * avgi * sd;
  }

 if (in.n_it > 1) 
 {
    double ksi_test = (in.s2 - in.s1 * in.s1 / in.s0) / (in.n_it - 1);
    double cs = in.s1 / in.s0;
    double cse = 100. * sqrt (in.s0) / fabs (in.s1);
    scrcolor (FGmain, BGmain);
    printLn (prtfile, &nln, " < >   %12.4E     %10.2E %8d %8.1G",
    cs, cse, in.nCallTot, ksi_test);
  }
  set_vegas_integral (in);
  return nln;
}



static int gui_vegas_new (FILE * prtfile, int nLine) {
  int i;
  int stop;
  int nln = nLine;
  double sd;
  double avgi;
  vegas_integral in = get_vegas_integral ();
  mcintr mc = get_mc_info ();

  if (in.n_it == 0) {
    init_vegas_integral ();
    in = get_vegas_integral ();
  }


  for (i = 0; i < mc.itmx0; ++i) {
    int nCubs;
    char errtxt[128] = "";

    nCall = 0;
    negPoints = 0;
    badPoints = 0;
    hFill = 1;
    nCubs = generateVegasCubs (veg_Ptr, mc.ncall0 / 2);
    setStopSymb (0);

    for (;;) {
      int cCub = (int)(nCubs * getCurCub ());
      stop = informline (cCub, nCubs);
      setStopSymb (stop);
      if (1. == getCurCub () || stop) break;
    }
    if (stop) break;
    negPoints /= nCall;
    badPoints /= nCall;
    in.old = 1;
    in.nCallTot += nCall;
    in.n_it++;
    scrcolor (FGmain, BGmain);
    printLn (prtfile, &nln, "%3d    %12.4E     %10.2E %8d ", in.n_it, avgi, 100 * sd / fabs (avgi), nCall);
    if (negPoints < 0) {
      sprintf (errtxt, " Negative points %.1G%%;", -100. * negPoints / (avgi - 2 * negPoints));
    }
    if (badPoints) {
      sprintf (errtxt + strlen (errtxt), "Bad Precision %.1G%%;", 100. * badPoints / (avgi - 2 * negPoints));
    }

    if (errtxt[0]) {
      scrcolor (Red, BGmain);
      printLn (prtfile, &nln, "%s", errtxt);
    }
    sd = 1 / (sd * sd);
    in.s0 += sd;
    in.s1 += avgi * sd;
    in.s2 += avgi * avgi * sd;
  }

  if (in.n_it > 1) {
    double ksi_test = (in.s2 - in.s1 * in.s1 / in.s0) / (in.n_it - 1);
    double cs = in.s1 / in.s0;
    double cse = 100. * sqrt (in.s0) / fabs (in.s1);
    scrcolor (FGmain, BGmain);
    printLn (prtfile, &nln, " < >   %12.4E     %10.2E %8d %8.1G",
    cs, cse, in.nCallTot, ksi_test);
  }
  set_vegas_integral (in);
  return nln;
}


static double 
crosssec (void)
{   

  int i;
  
  double cs;
  double sd;
  double avgi;
/*---------------------------------*/
  int ndim = imkmom ();
/*  int nses = get_nsession (); */
/*---------------------------------*/


  vegas_integral in = get_vegas_integral ();
  mcintr mc = get_mc_info ();

/*----------------------------------------------------------------------------*/
   vegas_finish (veg_Ptr);
   veg_Ptr = NULL;
   veg_Ptr = vegas_init (ndim, MAX_NDMX);
   correctHistList (0);
   mc.ncall0 = 2 * generateVegasCubs (veg_Ptr, mc.ncall0 / 2.);
   set_mc_info (mc);
/*----------------------------------------------------------------*/


 /* if (in.n_it == 0) */ 
 {
    init_vegas_integral ();
    in = get_vegas_integral ();
  }


  for (i = 0; i < mc.itmx0; ++i) {
    char errtxt[100] = "";

    nCall = 0;
    negPoints = 0;
    badPoints = 0;
    hFill = 1;
/*    printf("\n\n subprocess number = %d\n\n", proces_1.nsub); */
    if (vegas_int (veg_Ptr, mc.ncall0, proces_1.nsub, mc.itmx0, i + 1, me2_func, &avgi, &sd)) {
      break;
    }
 /*   in.old = 1; */

    negPoints /= nCall;
    badPoints /= nCall;
    in.nCallTot += nCall;
    scrcolor (FGmain, BGmain);
   /* ++in.n_it; */

    if (negPoints < 0)  sprintf (errtxt, " Negative points %.1G%%;", -100 * negPoints / (avgi - 2 * negPoints));
    if (badPoints)  sprintf (errtxt + strlen (errtxt), "Bad Precision %.1G%%;", 100 * badPoints / (avgi - 2 * negPoints));
    if (errtxt[0]) scrcolor (Red, BGmain);

    sd = 1 / (sd * sd);
    in.s0 += sd;
    in.s1 += avgi * sd;
    in.s2 += avgi * avgi * sd;
  }

 /* if (in.n_it > 1) */ 
   
    cs = in.s1 / in.s0;
  
  set_vegas_integral (in);
/*  printf("\n cs= %e\n",cs); */
  return cs;
}

int
runVegas (int init)
{
  int subprocnum;
  char mess[64];
  char fname[64];
  FILE *iprt = NULL;
  int mode;
  void *pscr = NULL;
  static int n_Line = 7;
  int ndim = imkmom ();
  int nses = get_nsession ();
  double cross_s;


  vegas_integral in = get_vegas_integral ();
  mcintr mc = get_mc_info ();

  if (2 == nin_) {
    strcpy (mess, " Cross section [pb] ");
  } else {
    strcpy (mess, "  Width     [GeV]   ");
  }

/* save current session parameters */
  write_session ();

/* open protocol file */
  sprintf (fname, "%sprt_%d", outputDir, nses);
  iprt = fopen (fname, "a");
  if (!ftell (iprt)) {
    fprintf (iprt, "    CompHEP kinematics module \n The session parameters:\n");
    write_prt (iprt);
    fprintf (iprt, "===================================\n");
  }

/* init kinematics */
  if (init && veg_Ptr)
    {
      vegas_finish (veg_Ptr);
      veg_Ptr = NULL;
    }
  if (!veg_Ptr)
    veg_Ptr = vegas_init (ndim, MAX_NDMX);

  correctHistList (0);

  mc.ncall0 = 2 * generateVegasCubs (veg_Ptr, mc.ncall0 / 2.);
  set_mc_info (mc);

/* Main loop */
  for (;;)
    {
      char strmen[] = "\030"
        " Itmx   =    N2         "
        " nCall  =    N1         "
        " Use Form Factor   KEY1 "
        " Distributions          "
        " Start integration      "
        " Clear statistic        "
        " Clear  grid            "
        " Generate events        "
        " All subprocess integr. "
        " 2 parameter dependence "
        " 1 parameter dependence ";       

      mc = get_mc_info ();
      improveStr (strmen, "N1", "%d", mc.ncall0);
      improveStr (strmen, "N2", "%d", mc.itmx0);
      
     /* improveStr (strmen, "KEY1", "%d", userFactor_key); */
     if(userFactor_key==0) improveStr (strmen, "KEY1", "NO");
     else improveStr (strmen, "KEY1", "YES");

      menu1 (54, 7, "", strmen, "n_veg_*", &pscr, &mode);
      switch (mode)
        {
        case 0:
          if (iprt) fclose (iprt);
          return 0;

        case 1:
          correctInt (50, 12, "Enter new value ", &mc.itmx0, 1);
          set_mc_info (mc);
          break;

        case 2:
          correctLong (50, 12, "Enter new value ", &mc.ncall0, 1);
          mc.ncall0 = 2 * generateVegasCubs (veg_Ptr, mc.ncall0 / 2.);
          set_mc_info (mc);
          break;

        case 3:
        /*  correctInt (50, 12, "Enter new value ", &userFactor_key, 1); */
        /*  set_mc_info (mc); */
          if(userFactor_key==0) userFactor_key=1;
          else userFactor_key=0;
          break;

        case 4:
          manipulateHists ();
          break;

        case 5:
          if (!in.old || 7 == n_Line) {
            n_Line = 7;
            scrcolor (Blue, BGmain);
            printLn (iprt, &n_Line, "#IT  %20s Error %%    nCall   chi**2", mess);
          }
          n_Line = gui_vegas (iprt, n_Line);
          init = 1;
          messanykey (54, 11, "Integration is over");
          break;

        case 6:
          {
	    init_vegas_integral ();
            clearHists ();
            messanykey (54, 13, "Old results for integral\n"
                        "and distributions\nare deleted.\n"
                        "But grid is saved");
          }
	  break;

        case 7:
          {
	    vegas_finish (veg_Ptr);
            veg_Ptr = vegas_init (ndim, MAX_NDMX);
            init = 1;
            messanykey (57, 11, "OK");
          }
	  break;

        case 8:
          {
            sprintf (fname, "%sevents_%d.txt", outputDir, nses);
            hFill = 0;
            menu_EventGenerator (veg_Ptr, me2_func, fname, iprt, init);
            init = 0;
          }
          break;

        case 9:

         
         for (subprocnum = 1;  subprocnum <= nprc_; subprocnum++)
          { 
            proces_1.nsub = subprocnum;
            ComposeSubprocessString ();
            if (nin_ == 2)
	    {
	     vshortstr name[2];
	     pinf_ (subprocnum, 1, name[0], NULL);
	     pinf_ (subprocnum, 2, name[1], NULL);
	     initStrFun (name[0], name[1]);
	    }
            printf("\n procnum= %d\n",subprocnum);



            if (!in.old || 7 == n_Line) {
               n_Line = 7;
               scrcolor (Blue, BGmain);
               printLn (iprt, &n_Line, "#IT  %20s Error %%    nCall   chi**2", mess);
             }
             n_Line = gui_vegasAll (iprt, n_Line);
             init = 1;
/*             messanykey (54, 11, "Integration is over"); */
           
/*             manipulateHists();-> showHistAll();-> plot_histo(hist->hMin,hist->hMax,nBin,f,df,proces_1.proces,xname,yname);->
                 writeroothisto (xMin, xMax, dim, f, ff, upstr, xstr, ystr);
*/
           showHistAll();

           init_vegas_integral();
           clearHists ();
           vegas_finish (veg_Ptr);
           veg_Ptr = vegas_init (ndim, MAX_NDMX);
           init = 1;  
           }       
         
         break;

        case 10:
          paramtable2 ();
          if (mess_y_n (25, 10, "Run table calculations?")) 
            { 
             paramdependence (crosssec, processch, "Cross section [GeV]");
            }        

/*
          asgn_ (3, 9.1188);
          asgn_ (5, 1.743);
          cross_s = crosssec();
          printf("\n cross_s1 = %e \n", cross_s); 
*/    
 
          init = 1;  
          break;
      

        case 11:
          paramdependence1 (crosssec, processch, "Cross section [GeV]");            
          init = 1;  
          break;

/*        case 8:
          {
            menu_generator_weighted_events (me2_func, fname);
            init = 0;
          }
          break;
*/     

        }
    }
}
