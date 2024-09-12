/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"
#include "plot/include/plot.h"
#include "out_ext.h"

#include "core_data.h"
#include "alphas_menu.h"
#include "kininpt.h"
#include "cut.h"
#include "regul.h"
#include "rw_sess.h"
#include "strfun.h"
#include "strfun_par.h"
#include "cs_22.h"
#include "histogram.h"
#include "q_kin.h"
#include "subproc.h"
#include "param.h"
#include "evnt_menu.h"
#include "vegas.h"
#include "runVegas.h"
#include "mc_menu.h"

static int
subprocess_menu (void)
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
	      memcpy (strmen + i * width_ + j * width + 2 + nlin + 2, name,
		      strlen (name));
	    }
	  memcpy (strmen + i * width_ + nin_ * width + 2 + nlin + 2, " -> ",
		  4);
	  for (j = nin_; j < nin_ + nout_; ++j)
	    {
	      pinf_ (i + 1, j + 1, name, NULL);
	      memcpy (strmen + i * width_ + j * width + 6 + nlin + 2, name,
		      strlen (name));
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


static void
f7_prog (int mode)
{
  int pos = 1;
  double be[2];
  void *pscr = NULL;
  if (mode > 2)
    {
      messanykey (10, 15, " Highlight the corresponding\n"
		  "structure function");
      return;
    }

  if (!get_sf_num (mode - 1))
    {
      return;
    }

  f3_key[4] = NULL;

  for (;;)
    {
      static double xMin = 0.0, xMax = 1.0, scale = 10;
      static int nPoints = 100;

      char strmen[] = "\030 "
	" x-Min = XXX            "
	" x-Max = YYY            "
	" Npoints = NNN          "
	" QCD-scale= QQQ         "
	" Display plot           ";

      improveStr (strmen, "XXX", "%.3f", xMin);
      improveStr (strmen, "YYY", "%.3f", xMax);
      improveStr (strmen, "NNN", "%d", nPoints);
      improveStr (strmen, "QQQ", "%.1fGeV", scale);

      menu1 (54, 14, "", strmen, "n_alpha_view", &pscr, &pos);

      switch (pos)
	{
	case 0:
	  f3_key[4] = f7_prog;
	  return;
	case 1:
	  correctDouble (55, 18, "xMin = ", &xMin, 1);
	  break;
	case 2:
	  correctDouble (55, 18, "xMax = ", &xMax, 1);
	  break;
	case 3:
	  correctInt (50, 18, "nPoints = ", &nPoints, 1);
	  break;
	case 4:
	  correctDouble (50, 18, "QCD-scale = ", &scale, 1);
	  break;
	case 5:
	  if (xMin >= 0. && xMax > xMin && xMax <= 1. && nPoints >= 3
	      && nPoints <= 150 && scale > 0.5)
	    {
	      double f[150];
	      void *screen;
	      double dx = (xMax - xMin) / (2 * nPoints);
	      int i;
	      get_text (1, 1, maxCol (), maxRow (), &screen);
	      for (i = 0; i < nPoints; ++i)
		{
		  double x = xMin + (i + 0.5) * (xMax - xMin) / nPoints;
		  f[i] = strfun_ (0, x, 1, scale) * pow (1 - x, be[mode - 1] - 1.);	/* #-mdl */
		}

	      {
		vshortstr p_name;
		shortstr mess;
		shortstr beam;
		shortstr pdf;
		pinf_ (proces_1.nsub, mode, p_name, NULL);
		trim (p_name);
		strcat (p_name, "(x)");
		strFunName (mode, beam, pdf);
		trim (beam);
		trim (pdf);
		sprintf (mess, "%s(%s)", beam, pdf);
		sprintf (mess + strlen (mess), " [QCD-scale = %.1f GeV]", scale);
		plot_histo (xMin + dx, xMax - dx, nPoints, f, NULL, mess, "x", p_name);
	      }
	      put_text (&screen);
	    }
	  else
	    {
	      messanykey (16, 5, " Correct input is \n"
			  "  0<=xMin<xMax<=1,\n"
			  " QCD-scale > 0.5 GeV\n" " 2 < nPoints < 201");
	    }
	  break;
	}
    }
}


static int
in_setting (void)
{
  int mode = 1;
  int i;
  int returnCode = 0;
  double sqrt_S;
  void *pscr = NULL;
  void (*f7_tmp) (int) = f3_key[4];
  char *f7_mess_tmp = f3_mess[4];
  double e[2], p[2], mass[2];
  double rap = get_rapidity ();

  if (nin_ != 2)
    {
      return returnCode;
    }

  for (i = 0; i < 2; i++)
    {
      if (get_sf_num (i))
	{
	  mass[i] = get_sf_mass (i);
	}
      else
	{
	  pinf_ (proces_1.nsub, i + 1, NULL, mass + i);
	}
    }

  vinf_ (0, NULL, &sqrt_S);
  for (i = 0; i < 2; i++)
    {
      e[i] =
	(sqrt_S * sqrt_S + mass[i] * mass[i] -
	 mass[1 - i] * mass[1 - i]) / (2 * sqrt_S);
      p[i] = sqrt (e[i] * e[i] - mass[i] * mass[i]);
      p[i] = p[i] * cosh (rap) + e[i] * sinh (rap) * (1 - 2 * i);
      if (p[i] * p[i] < (10.E-10) * sqrt_S)
	{
	  p[i] = 0;
	}
    }

  for (;;)
    {
      shortstr beam_txt;
      shortstr pdf_name;
      char strmen[] = "*"
	" Beam particle 1: XXX1                              "
	" Beam particle 2: XXX2                              "
	" Str.Fun.1: YYY1                                    "
	" Str.Fun.2: YYY2                                    "
	" 1 particle momentum[GeV] = PPP1                    "
	" 2 particle momentum[GeV] = PPP2                    ";
      vinf_ (0, NULL, &sqrt_S);

      strmen[0] = (strlen (strmen) - 1) / 6;

      strFunName (1, beam_txt, pdf_name);
      trim (beam_txt);
      trim (pdf_name);
      improveStr (strmen, "XXX1", "%-34.34s", beam_txt);
      improveStr (strmen, "YYY1", "%-34.34s", pdf_name);
      strFunName (2, beam_txt, pdf_name);
      trim (beam_txt);
      trim (pdf_name);
      improveStr (strmen, "XXX2", "%-34.34s", beam_txt);
      improveStr (strmen, "YYY2", "%-34.34s", pdf_name);
      improveStr (strmen, "PPP1", "%-10.4G", p[0]);
      improveStr (strmen, "PPP2", "%-10.4G", p[1]);

      f3_key[4] = f7_prog;
      f3_mess[4] = "Plot";
      menu1 (26, 7, "", strmen, "n_in_*", &pscr, &mode);
      f3_key[4] = f7_tmp;
      f3_mess[4] = f7_mess_tmp;

      switch (mode)
	{
	case 0:
	  for (i = 0; i < 2; i++)
	    {
	      if (get_sf_num (i))
		{
		  mass[i] = get_sf_mass (i);
		}
	      else
		{
		  pinf_ (proces_1.nsub, i + 1, NULL, mass + i);
		}
	      if (mass[i] == 0 && p[i] <= 0)
	        {
	          warnanykey (10, 10,
			      "Momentum must be positive \n"
			      "for a zero mass particle. Momentum is set to 1 GeV");
	          p[i] = 1.;
	        }
	    }

	  for (i = 0; i < 2; i++)
	    e[i] = sqrt (p[i] * p[i] + mass[i] * mass[i]);
	  sqrt_S = sqrt ((e[0] + e[1]) * (e[0] + e[1]) - (p[0] - p[1]) * (p[0] - p[1]));
	  asgn_ (0, sqrt_S);
	  set_rapidity (atanh ((p[0] - p[1]) / (e[0] + e[1])));
	  if (returnCode)
	    {
	      vshortstr name[2];
	      pinf_ (proces_1.nsub, 1, name[0], NULL);
	      pinf_ (proces_1.nsub, 2, name[1], NULL);
	      initStrFun (name[0], name[1]);
	    }
	  return returnCode;
	  break;
	case 1:
	case 2:
	  if (beam_menu (mode))
	    returnCode = 1;
	  break;
	case 3:
	case 4:
	  if (pdf_menu (mode - 2))
	    returnCode = 1;
	  break;
	case 5:
	case 6:
	  if (correctDouble (50, 12, "Enter new value ", p + mode - 5, 1)) 
	    returnCode = 1;
	  break;
	}
    }
}


static int
width_menu ()
{
  int key = 0;
  void *pscr = NULL;
  char strmen[] = "\012" " Fixed    " " Overall  " " Running  ";

  menu1 (64, 7, "", strmen, "n_width_*", &pscr, &key);
  clrbox (1, 5, 160, 50);
  if (key)
    {
      gwidth = key - 1;
      return 1;
    }

  return 0;
}

static int
checkEnergy (void)
{
  int i;
  double mini = 0.0;
  double mfnl = 0.0;
  int to_little_energy = 0;

  for (i = 0; i < nin_ + nout_; ++i)
    {
      double them;
      pinf_ (proces_1.nsub, i + 1, NULL, &them);
      if (i < nin_)
        mini += them;
      else 
        mfnl += them;
    }

  if (nin_ == 1)
    {
      if (mini <= mfnl)
	to_little_energy = 1;
    }
  else
    {
      double sqrt_S;
      vinf_ (0, NULL, &sqrt_S);
      if (sqrt_S <= mfnl || sqrt_S <= mini)
	to_little_energy = 1;
    }

  return to_little_energy;
}

void
monte_carlo_menu (void)
{
  int r = 0;
  int mode = 1;
  void *pscr = NULL;
  void *pscr_mem = NULL;
  char menutxt[] = "\030"
    " Subprocess             "
    " Initial state          "
    " Model parameters       "
    " Constraints            "
    " QCD scale              "
    " Width scheme:  XXXXXXX "
    " Cuts                   "
    " Kinematics             "
    " Regularization         "
    " Numerical Session      " " Simpson                ";
  char ss[3][30] = { "Fixed  ", "Overall", "Running", };
  int gw = gwidth;

  improveStr (menutxt, "XXXXXXX", ss[gwidth]);

  get_text (1, 10, 80, 24, &pscr_mem);
  ComposeSubprocessString ();
  for (;;)
    {
      vegas_integral in = get_vegas_integral ();
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
      if (in.old)
	{
	  print ("(continue)");
	}
      else
	{
	  print ("(begin)");
	  scrcolor (FGmain, BGmain);
	  clrbox (1, 8, 53, maxRow () - 2);
	}

      menutxt[menutxt[0] * 10 + 1] = ' ';  /* menu opton range: Subprocess -  Simpson */
      if (nin_ != 2 || nout_ != 2 || get_sf_num (0) || get_sf_num (1)){
	menutxt[menutxt[0] * 10 + 1] = 0;  /* menu opton range: Subprocess - Vegas */
      }
      menu1 (54, 4, "", menutxt, "n_mc_*", &pscr, &mode);

      switch (mode)
	{
	case 0:
	  put_text (&pscr_mem);
	  return;
	case 1:
	  r = r | (3 * subprocess_menu ());
	  break;
	case 2:
	  if (nin_ == 2)
	    r = r | in_setting ();
	  break;
	case 3:
	  r = r | change_parameter (54, 7);
	  break;
	case 4:
	  show_depend (54, 7);
	  break;
	case 5:
	  r = r | qcdmen_ ();
	  break;
	case 6:
	  {
	    r = r | width_menu ();
	    improveStr (menutxt, ss[gw], ss[gwidth]);
	    gw = gwidth;
	    break;
	  }
	case 7:
	  do
	    r = r | (3 * edittable (1, 4, &cutTab, 1, "n_cut", 0));
	  while (fillCutArray ());
	  break;
	case 8:
	  r = r | 3 * EnterKinScheme ();
	  break;
	case 9:
	  do
	    r = r | (3 * edittable (1, 4, &regTab, 1, "n_reg", 0));
	  while (fillRegArray ());
	  break;
	case 10:
	case 11:
	  {
	    if (calcFunc ())
	      {
		warnanykey (15, 15, "Can not evaluate constraints");
		break;
	      }
	    if (fillCutArray ())
	      {
		warnanykey (15, 15, "Can not evaluate cuts limits");
		break;
	      }
	    if (fillRegArray ())
	      {
		warnanykey (15, 15,
			    "Can not evaluate regularization parameters");
		break;
	      }
	    if (checkEnergy ())
	      {
		warnanykey (15, 15, "Energy is too small");
		break;
	      }

	    if (mode == 11)
	      {			/* executed if Simpson is in the menu only! */
		cs_numcalc ();
	      }
	    else
	      {
		if (r & 2)
		  {
		    runVegas (1);
		  }
		else
		  {
		    runVegas (0);
		  }
	      }
	    r = 0;
	    break;
	  }
	}

      in = get_vegas_integral ();
      if ((r & 1) && in.old)
	{
	  messanykey (15, 15,
		      "Some parameters where changed.\nSo integral and statictics for\n"
		      "distribushions is forgotten!\nSession number is increased.");
	  clearHists ();
	  init_vegas_integral ();
	  ClearEventMax ();
	  ClearVegasGrid ();
	  set_nsession (get_nsession () + 1);
	  write_session ();
	}
    }
}

void
oneclick_monte_carlo_menu (void)
{
  int mode = 1;
  void *pscr = NULL;
  void *pscr_mem = NULL;
  char ss[3][32] = { "Fixed  ", "Overall", "Running", };
  int gw = gwidth;

  get_text (1, 10, 80, 24, &pscr_mem);
  ComposeSubprocessString ();
  for (;;) {
    mcintr mc = get_mc_info ();
    char menutxt[] = "\030"
      " Subprocess             "
      " Initial state          "
      " Model parameters       "
      " Constraints            "
      " QCD scale              "
      " Width scheme:  XXXXXXX "
      " Cuts                   "
      " Itmx   =    N2         "
      " nCall  =    N1         "
      " New batch              "
      ;

      vegas_integral in = get_vegas_integral ();
      goto_xy (4, 3);
      clr_eol ();
      scrcolor (Red, BGmain);
      print ("(sub)Process: ");
      scrcolor (FGmain, BGmain);
      print (proces_1.proces);
      goto_xy (4, 4);
      clr_eol ();

      if (in.old)
	{
	  print ("(continue)");
	}
      else
	{
	  print ("(begin)");
	  scrcolor (FGmain, BGmain);
	  clrbox (1, 8, 53, maxRow () - 2);
	}

      improveStr (menutxt, "XXXXXXX", ss[gwidth]);
      improveStr (menutxt, "N1", "%d", mc.ncall0);
      improveStr (menutxt, "N2", "%d", mc.itmx0);

      menu1 (54, 4, "", menutxt, "n_mc_*", &pscr, &mode);
      switch (mode)
	{
	case 0:
	  put_text (&pscr_mem);
	  return;
	case 1:
	  subprocess_menu ();
	  break;
	case 2:
	  if (nin_ == 2)
	    in_setting ();
	  break;
	case 3:
	  change_parameter (54, 7);
	  break;
	case 4:
	  show_depend (54, 7);
	  break;
	case 5:
	  qcdmen_ ();
	  break;
	case 6:
	  {
	    width_menu ();
	    improveStr (menutxt, ss[gw], ss[gwidth]);
	    gw = gwidth;
	    break;
	  }
	case 7:
	  do
	    edittable (1, 4, &cutTab, 1, "n_cut", 0);
	  while (fillCutArray ());
	  break;
        case 8:
          correctInt (50, 12, "Enter new value ", &mc.itmx0, 1);
          set_mc_info (mc);
          break;

        case 9:
          correctLong (50, 12, "Enter new value ", &mc.ncall0, 1);
          set_mc_info (mc);
          break;

	case 10:
          {
	    int err = 0;
	    FILE * batchfile = fopen (scat ("%s/batch.dat", pathtoresults), "a");
	    write_session ();
	    system (scat ("cp session.dat %s/.", pathtoresults));
	    if (0 == ftell(batchfile)) {
               err = system (scat ("%s/num_batch.pl -d %s", pathtouser, pathtoresults));
	       messanykey (15, 15, scat ("There is no %s/batch.dat The file has been created.\n"
	                                 "Configure CompHEP with %s/n_comphep and num_batch.pl\n", 
	                                  pathtoresults, pathtoresults));
	    } else {
              err = system (scat ("%s/num_batch.pl -d %s --geninfo", pathtouser, pathtoresults));
	    }
	    if (err)
	      messanykey (15, 15, "Alas! There is a problem with num_batch.pl. \n"
	                          "          Consult CompHEP authors!\n");
	  }
	  break;
	}
    }
}
