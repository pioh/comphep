/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* Author:  S.Shichanin
* ------------------------------------------------------
*/
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/4_vector.h"
#include "service2/include/files.h"
#include "service2/include/drandXX.h"
#include "chep_crt/include/chep_crt.h"

#include "plot/include/plot.h"

#include "strfun.h"
#include "rw_sess.h"
#include "core_data.h"
#include "subproc.h"
#include "err_code.h"
#include "num_serv.h"
#include "const.h"
#include "param.h"
#include "out_ext.h"
#include "param.h"
#include "evnt_menu.h"
#include "width_12.h"
#include "roothisto.h"

typedef struct canal_12
  {
    char prtclnames[2][6];
    double width;
  }
canal_12;

static canal_12 *allcanal;

static void 
sortcanals (void)
{
  int i;
  canal_12 buff;

  i = 0;
  while (i < nprc_ - 1)
    {
      if (allcanal[i].width < allcanal[i + 1].width)
	{
	  buff = allcanal[i];
	  allcanal[i] = allcanal[i + 1];
	  allcanal[i + 1] = buff;
	  if (i != 0)
	    i--;
	}
      else
	i++;
    }
}


static void 
decay12information (double locw12)
{
  byte xcount, ycount;
  int i;
  char fname[64];
  FILE * iprt;

  sprintf (fname, "%sprt_%d", outputDir, get_nsession ());
  iprt = fopen (fname, "a");
  if (!ftell (iprt)) {
    fprintf (iprt, "    CompHEP kinematics module \n The session parameters:\n");
    write_prt (iprt);
    fprintf (iprt, "===================================\n");
  }

  clrbox (1, 16, 80, 24);
  goto_xy (1, 16);
  scrcolor (Red, BGmain);
  print (" Total width : ");
  scrcolor (FGmain, BGmain);
  fprintf (iprt, "#TotalWidth\n");
  if (err_code == 0) {
    fprintf (iprt, "  Wtot = %.17E GeV\n", locw12);
    print ("%13E GeV", locw12);
  } else {
    fprintf (iprt, "  incorrect\n");
    print (" incorrect       ");
  }
  fprintf (iprt, "#Branchings\n");
  if (locw12 > 1.E-20 && err_code == 0)
    {
      sortcanals ();
      goto_xy (1, 17);
      scrcolor (Red, BGmain);
      print (" Modes and fractions :");
      xcount = 31;
      ycount = 18;
      for (i = 0; i < nprc_; i++)
	{
	  if (ycount < maxRow ())
	    {
	      goto_xy (xcount, ycount - 1);
	      scrcolor (Blue, BGmain);
	      print ("%s %s -  ", allcanal[i].prtclnames[0],
		     allcanal[i].prtclnames[1]);

	      scrcolor (FGmain, BGmain);
	      print ("  %5.2lG%%", 100 * allcanal[i].width / locw12);
	      if (0.0 < allcanal[i].width)
	      fprintf (iprt, "  Br (%s %s) = %f%%\n", 
	        allcanal[i].prtclnames[0],
		allcanal[i].prtclnames[1],
		100. * allcanal[i].width / locw12);
	      xcount += 30;
	      if (xcount > maxCol ())
		{
		  xcount = 1;
		  ++(ycount);
		}
	    }
	}
    }
  scrcolor (FGmain, BGmain);
  if (iprt) fclose (iprt);
}


static double 
calcwidth12_partial (void)
{
  int i, nsub_c;
  double pmass[3];

  err_code = 0;
  double w12 = 0.0;

  nsub_c = proces_1.nsub - 1;
    
      for (i = 0; i < 3; i++) {
	pinf_ (nsub_c + 1, i + 1, NULL, pmass + i);
      }

      if (pmass[0] <= pmass[1] + pmass[2]) {
	allcanal[nsub_c].width = 0.0;
      } else {
	double mm1 = pmass[0] * pmass[0];
	double mm2 = pmass[1] * pmass[1];
	double mm3 = pmass[2] * pmass[2];
	double pRestOut = sqrt ((mm1 - (mm2 + mm3 + 2. * pmass[1] * pmass[2])) *
	        	 (mm1 - (mm2 + mm3 - 2. * pmass[1] * pmass[2]))
	  ) / (2 * pmass[0]);
	double totcof = pRestOut / (8. * M_PI * mm1);

	for (i = 1; i < 12; i++)
	  pvect[i] = 0;
	pvect[0] = pmass[0];
	pvect[7] = pRestOut;
	pvect[4] = sqrt (pvect[7] * pvect[7] + mm2);
	pvect[11] = -pvect[7];
	pvect[8] = sqrt (pvect[7] * pvect[7] + mm3);

	allcanal[nsub_c].width = totcof * sqme_ (nsub_c + 1, pvect, &err_code);
      }
      w12 = allcanal[nsub_c].width;
      pinf_ (nsub_c + 1, 2, allcanal[nsub_c].prtclnames[0], NULL);
      pinf_ (nsub_c + 1, 3, allcanal[nsub_c].prtclnames[1], NULL);
    

  if (err_code != 0) {
    errormessage ();
    return 0.;
  }

  return w12;
}



static double 
calcwidth12 (void)
{
  int i, nsub;
  double pmass[3];

  err_code = 0;
  double w12 = 0.0;

  for (nsub = 0; nsub < nprc_; ++nsub)
    {
      for (i = 0; i < 3; i++) {
	pinf_ (nsub + 1, i + 1, NULL, pmass + i);
      }

      if (pmass[0] <= pmass[1] + pmass[2]) {
	allcanal[nsub].width = 0.0;
      } else {
	double mm1 = pmass[0] * pmass[0];
	double mm2 = pmass[1] * pmass[1];
	double mm3 = pmass[2] * pmass[2];
	double pRestOut = sqrt ((mm1 - (mm2 + mm3 + 2. * pmass[1] * pmass[2])) *
	        	 (mm1 - (mm2 + mm3 - 2. * pmass[1] * pmass[2]))
	  ) / (2 * pmass[0]);
	double totcof = pRestOut / (8. * M_PI * mm1);

	for (i = 1; i < 12; i++)
	  pvect[i] = 0;
	pvect[0] = pmass[0];
	pvect[7] = pRestOut;
	pvect[4] = sqrt (pvect[7] * pvect[7] + mm2);
	pvect[11] = -pvect[7];
	pvect[8] = sqrt (pvect[7] * pvect[7] + mm3);

	allcanal[nsub].width = totcof * sqme_ (nsub + 1, pvect, &err_code);
      }
      w12 += allcanal[nsub].width;
      pinf_ (nsub + 1, 2, allcanal[nsub].prtclnames[0], NULL);
      pinf_ (nsub + 1, 3, allcanal[nsub].prtclnames[1], NULL);
    }

  if (err_code != 0) {
    errormessage ();
    return 0.;
  }

  return w12;
}


static double 
partialwidth (void)
{
  return calcwidth12_partial ();
}


static double 
totwidth (void)
{
  return calcwidth12 ();
}


static double
me2_1to2_func (double *x, double wgt)
{
  int i;
  int err = 0;
  double pmass[3];
  double ret_val = 0.;

  for (i = 0; i < 3; i++) {
    pinf_ (proces_1.nsub, i + 1, NULL, pmass + i);
  }

  if (pmass[0] <= pmass[1] + pmass[2]) {
    ret_val = 0.0;
  } else {
    double mm1 = pmass[0] * pmass[0];
    double mm2 = pmass[1] * pmass[1];
    double mm3 = pmass[2] * pmass[2];
    double pRestOut = sqrt ((mm1 - (mm2 + mm3 + 2. * pmass[1] * pmass[2])) *
   		    (mm1 - (mm2 + mm3 - 2. * pmass[1] * pmass[2]))
      ) / (2 * pmass[0]);
    double totcoef = pRestOut / (8. * M_PI * mm1);

    for (i = 1; i < 12; i++)
      pvect[i] = 0.0;
    pvect[0] = pmass[0];
    pvect[7] = pRestOut;
    pvect[4] = sqrt (pvect[7] * pvect[7] + mm2);
    pvect[11] = -pvect[7];
    pvect[8] = sqrt (pvect[7] * pvect[7] + mm3);

/* 2-d randomization, the last step in 3-d rand. is done in event writing out */
   if (1) {
      double z_rnd = 1 - 2 * drandXX ();
      double r_rnd = sqrt (1 - z_rnd * z_rnd);
      for (i = 1; i < 3; ++i) {
        double ZZ = pvect[4 * i + 3];
        pvect[4 * i + 1] = ZZ * r_rnd;
        pvect[4 * i + 3] = ZZ * z_rnd;
      }
    }
    ret_val = totcoef * sqme_ (proces_1.nsub, pvect, &err);
  }

  return ret_val * wgt;
}

static int
decay_subprocess_menu (void)
{
  int i, j;
  int mode = 0;
  int width;
  int width_;
  char *strmen;
  void *pscr = NULL;
  int nlin;
  char name[16];
  char lin[16];
  char frmt[16];

  sprintf (lin, "%i", nprc_);
  nlin = strlen (lin);

  if (nprc_ > 1)
    {
      width = 0;
      for (i = 0; i < nprc_; ++i)
        {
          for (j = 0; j < nin_ + nout_; ++j)
            {
              pinf_ (i + 1, j + 1, name, NULL);
              width = MAX (width, strlen (name));
            }
        }
      width++;
      width_ = 6 + 2 + nlin + width * (nin_ + nout_);

      strmen = malloc ((2 + nprc_ * width_) * sizeof (char));
      for (i = 1; i <= width_ * nprc_; ++i)
        {
          strmen[i] = ' ';
        }
      strmen[0] = width_;
      strmen[1 + nprc_ * width_] = 0;

      for (i = 0; i < nprc_; ++i)
        {
          sprintf (frmt, "%%-%i", nlin);
          strcat (frmt, "d. ");
          sprintf (lin, frmt, i + 1);
          memcpy (strmen + i * width_ + 1, lin, nlin + 2);
          for (j = 0; j < nin_; ++j)
            {
              pinf_ (i + 1, j + 1, name, NULL);
              memcpy (strmen + i * width_ + j * width + 2 + nlin + 2, name, strlen (name));
            }
          memcpy (strmen + i * width_ + nin_ * width + 2 + nlin + 2, " -> ", 4);
          for (j = nin_; j < nin_ + nout_; ++j)
            {
              pinf_ (i + 1, j + 1, name, NULL);
              memcpy (strmen + i * width_ + j * width + 6 + nlin + 2, name, strlen (name));
            }
        }
      menu1 (54, 4, "Subprocesses", strmen, NULL, &pscr, &mode);
      free (strmen);
    }

  if (mode)
    {
      put_text (&pscr);
      proces_1.nsub = mode;
      ComposeSubprocessString ();
      if (nin_ == 2)
        {
          vshortstr name[2];
          pinf_ (mode, 1, name[0], NULL);
          pinf_ (mode, 2, name[1], NULL);
          initStrFun (name[0], name[1]);
        }
      return 1;
    }
  return 0;
}

void 
decay12 (void)
{
  int k;
  void *pscr = NULL;
  char fname[64];
  allcanal = (canal_12 *) malloc (sizeof (canal_12) * nprc_);

  k = 1;
  double tot_width = calcwidth12 ();
  decay12information (tot_width);

  if (blind) return;
begg: 
  do {
      goto_xy (4, 3);
      clr_eol ();
      scrcolor (Red, BGmain);
      print ("(sub)Process: ");
      scrcolor (FGmain, BGmain);
      print (proces_1.proces);
      goto_xy (4, 4);
      clr_eol ();
      scrcolor (Red, BGmain);
      print ("Monte Carlo session: ");
      scrcolor (Black, BGmain);
      print ("%d", get_nsession ());
      menu1 (56, 6, "", "\026"
	     " Model parameters     "
             " Constraints          "
	     " 2 parameter depend.  "
	     " Sum 2d tables        "
	     " Event generator      ",
	     "n_13_*", &pscr, &k);
      switch (k) {
	case 0:
	  return;
	case 1:
	  if (change_parameter (55, 10))
	    {
              double new_width = calcwidth12 ();
	      set_nsession (get_nsession () + 1);
	      decay12information (new_width);
	    }
	  break;
	case 2:
	  show_depend (54, 7);
	  break;
	case 3:
/*	  paramdependence (totwidth, processch, "Width  [GeV]"); */
/*        paramdependence (partialwidth, processch, "Width  [GeV]");*/
          paramtable2 ();
          if (mess_y_n (25, 10, "Run table calculations?")) 
            { 
             paramdependence (partialwidth, processch, "Width  [GeV]");
            }     
	  break;
	case 4:
          combine_Sum3d ();
	  break;
	case 5:
	  k = 0;
	  break;
	}
  } while (k != 0);

  do
    {
      menu1 (56, 6, "", "\026"
             " Subprocess           "
	     " Events               ",
	     "n_13_*", &pscr, &k);
      switch (k)
	{
	case 0:
	  goto begg;
	case 1: {
            int changed = decay_subprocess_menu ();
	    if (changed) {
	      set_nsession (get_nsession () + 1);
              goto_xy (4, 3);
              clr_eol ();
              scrcolor (Red, BGmain);
              print ("(sub)Process: ");
              scrcolor (FGmain, BGmain);
              print (proces_1.proces);
              goto_xy (4, 4);
              clr_eol ();
              scrcolor (Red, BGmain);
              print ("Monte Carlo session: ");
              scrcolor (Black, BGmain);
              print ("%d", get_nsession ());
            }
	  }
	  break;
	case 2: {
            double * xxx;
	    vegas_integral in = get_vegas_integral ();
            FILE * iprt;
            sprintf (fname, "%sprt_%d", outputDir, get_nsession ());
            iprt = fopen (fname, "a");
            if (!ftell (iprt)) {
              fprintf (iprt, "    CompHEP kinematics module \n The session parameters:\n");
              write_prt (iprt);
              fprintf (iprt, "===================================\n");
            }
            sprintf (fname, "%sevents_%d.txt", outputDir, get_nsession ());
            in.n_it = 1;
            in.s0 = 1./me2_1to2_func (xxx, 1.);
            in.s1 = 1.;
            set_vegas_integral (in);
            menu_1to2_EventGenerator (me2_1to2_func, fname, iprt);
	    fclose (iprt);
	  }
	}
    }
  while (k != 0);
  free (allcanal);
  clrbox (1, 16, 80, 24);
}
