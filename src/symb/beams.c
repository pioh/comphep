/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ----------------------------------------------
*/
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/files.h"
#include "service2/include/syst.h"
#include "chep_crt/include/chep_crt.h"

#include "sos.h"
#include "screen.h"
#include "physics.h"
#include "process.h"
#include "read_mdl.h"
#include "beams.h"

static int
In_PDFbase_num (int num)
{
  int i;

  for (i = 0; i < n_strfun; i++)
    {
      if (num == strfunbase[i].num)
	{
	  return i;
	}
    }
  return -1;
}


static int
In_PDFbase_name (shortstr name)
{
  int i;

  for (i = 0; i < n_strfun; i++)
    {
      if (!strcmp (name, strfunbase[i].name))
	{
	  return i;
	}
    }
  return -1;
}


static int
In_hadronbase (shortstr name)
{
  int i;

  for (i = 0; i <= n_hadron; i++)
    {
      if (!strcmp (name, hadronbase[i].name))
	{
	  return i;
	}
    }
  return -1;
}

void construct_full_sf_name (int i)
{
  int len;
  shortstr beamname;
  shortstr sfname;
  char * bname = strstr (beam[i].h.sf_name, "(");

  if (NULL != bname) {
    strcpy (beamname, bname + 1);
    strcpy (sfname, strdup (beam[i].h.sf_name));
    len = strlen (beamname);
    beamname[len - 1] = '\0';
    sfname[strlen (beam[i].h.sf_name) - len - 1] = '\0';

    trim (beamname);
    trim (sfname);

    sprintf (beam[i].h.full_sf_name, "%s:%i:%i(%s)", sfname, beam[0].h.sf_set, beam[0].h.sf_mem, beamname);
  } else {
    strcpy (beam[i].h.full_sf_name, beam[i].h.sf_name);
  }
}

int
In_composite_particle_base (shortstr name)
{
  int i;

  for (i = 0; i <= n_cpart; i++)
    {
      if (!strcmp (name, cpartbase[i].name))
	{
	  return i;
	}
    }
  return -1;
}


int
copy_hadrons (hadron * dist, hadron src)
{
  int i, err = 0;

  strcpy (dist->name, src.name);
  strcpy (dist->sf_name, src.sf_name);
  dist->sf_set = src.sf_set;
  dist->sf_mem = src.sf_mem;
  dist->mass = src.mass;
  dist->how = src.how;
  for (i = 0; i < src.how; i++)
    dist->parton[i] = src.parton[i];

  return err;
}

int
get_hadroncontent (hadron h, shortstr res)
{
  int i, err = 0;
  vshortstr nm;
  shortstr tmp;

  if (!h.how)
    {
      strcpy (res, "");
      return 1;
    }

  res[0] = '\0';
  for (i = 0; i < h.how - 1; i++)
    {
      prtclname (h.parton[i], nm);
      sprintf (tmp, "%s,", nm);
      strcat (res, tmp);
    }
  prtclname (h.parton[i], nm);
  sprintf (tmp, "%s", nm);
  strcat (res, tmp);
  return err;
}


static int nTotB;
static int nFirstB;

static void
beamlist (int key)
{
  int i, k;
  int tabMax;
  int tabSz;
  int n_len, n_max_len;
  shortstr content;
  longstr str;

  tabMax = ycons - 7;
  if (-1 == key)
    {
      key = 0;
      for (k = 2; k <= 24; k++)
	{
	  goto_xy (1, k);
	  clr_eol ();
	}
    }

  if (key == 0)
    {
      scrcolor (FGmain, BGmain);
      for (k = 2; k < tabMax + 5; k++)
	{
	  goto_xy (1, k);
	  clr_eol ();
	}
      goto_xy (10, 3);
      scrcolor (Blue, BGmain);
      print ("List of beams (switch to the particle list by F3)");
      nTotB = 0;
      nFirstB = 0;
    }
  else
    {
      if (nTotB <= tabMax)
	return;
      switch (key)
	{
	case KB_DOWN:
	  nFirstB++;
	  break;
	case KB_UP:
	  nFirstB--;
	  break;
	case KB_PAGED:
	  nFirstB += tabMax;
	  break;
	case KB_PAGEU:
	  nFirstB -= tabMax;
	  break;
	}
      if (nFirstB < 0)
	nFirstB = 0;
      if (nTotB - nFirstB < tabMax)
	nFirstB = nTotB - tabMax;
      clrbox (1, 4, 79, 5 + tabMax);
    }
  goto_xy (3, 5);
  scrcolor (FGmain, BGmain);
  n_max_len = 0;
  for (k = 0; k < n_hadron; k++)
    {
      if (k >= nFirstB && (k - nFirstB) < tabMax)
	{
	  get_hadroncontent (hadronbase[k], content);
	  n_len = strlen (hadronbase[k].name);
	  n_len += strlen (content);
	  if (n_len > n_max_len)
	    n_max_len = n_len;
	}
    }

  for (k = 0; k < n_hadron; k++)
    {
      if (k >= nFirstB && (k - nFirstB) < tabMax)
	{
	  midstr bksp;
	  get_hadroncontent (hadronbase[k], content);
	  n_len = strlen (hadronbase[k].name);
	  n_len += strlen (content);
	  for (i = 0; i < n_max_len - n_len; i++)
	    bksp[i] = ' ';
	  bksp[i++] = '\0';
          if (strstr (hadronbase[k].sf_name, "LHA") || strstr (hadronbase[k].sf_name, "PDF")) {
	    sprintf (str, "Name: %s (%s) %s%s:%i:%i", hadronbase[k].name,
		   content, bksp, hadronbase[k].sf_name, hadronbase[k].sf_set, hadronbase[k].sf_mem);
	  } else {
	    sprintf (str, "Name: %s (%s) %sPDF: %s", hadronbase[k].name,
		   content, bksp, hadronbase[k].sf_name);
	  }
	  print ("%s", str);
	  str[0] = '\0';
	  goto_xy (3, where_y () + 1);
	}
    }
  nTotB = n_hadron;
  tabSz = MIN ((nTotB + 2), tabMax);
  chepbox (1, 4, 79, 5 + tabSz);

  if (nFirstB > 0)
    {
      goto_xy (72, 4);
      print ("PgUp");
    }
  if (nFirstB + tabSz < nTotB)
    {
      goto_xy (72, 5 + tabMax);
      print ("PgDn");
    }
  scrcolor (FGmain, BGmain);
}


static int nTotS;
static int nFirstS;

static void
strfunlist (int key)
{
  int k;
  int tabMax;
  int tabSz;

  tabMax = ycons - 7;
  if (key == 0)
    {
      scrcolor (FGmain, BGmain);
      for (k = 2; k <= 18; k++)
	{
	  goto_xy (1, k);
	  clr_eol ();
	}
      goto_xy (24, 3);
      scrcolor (Blue, BGmain);
      print ("List of structure functions");
      nTotS = 0;
      nFirstS = 0;
    }
  else
    {
      if (nTotS <= tabMax)
	return;
      switch (key)
	{
	case KB_DOWN:
	  nFirstS++;
	  break;
	case KB_UP:
	  nFirstS--;
	  break;
	case KB_PAGED:
	  nFirstS += tabMax;
	  break;
	case KB_PAGEU:
	  nFirstS -= tabMax;
	  break;
	}
      if (nFirstS < 0)
	nFirstS = 0;
      if (nTotS - nFirstS < tabMax)
	nFirstS = nTotS - tabMax;
      clrbox (1, 4, 79, 5 + tabMax);
    }
  goto_xy (3, 5);
  scrcolor (FGmain, BGmain);
  for (k = 0; k < n_strfun; k++)
    {
      if (k >= nFirstS && (k - nFirstS) < tabMax)
	{
          if (strstr (strfunbase[k].name, "LHA") || strstr (strfunbase[k].name, "PDF")) {
	    print ("%i: %s:%i:%i ", strfunbase[k].num, strfunbase[k].name, strfunbase[k].set, strfunbase[k].mem);
	  } else {
            print ("%i: %s ", strfunbase[k].num, strfunbase[k].name);
	  }
	  goto_xy (3, where_y () + 1);
	}
    }
  nTotS = n_strfun;
  tabSz = MIN (nTotS + 2, tabMax);
  chepbox (1, 4, 79, 5 + tabSz);

  if (nFirstS > 0)
    {
      goto_xy (72, 4);
      print ("PgUp");
    }
  if (nFirstS + tabSz < nTotS)
    {
      goto_xy (72, 5 + tabMax);
      print ("PgDn");
    }
  scrcolor (FGmain, BGmain);
}


static int
enter_beam (void)
{
  int i, j, jj;
  int y0;
  int redres;
  int screen_choise;
  int bm;
  int same_beams = 0, model_beam[2];
  int MaxR = maxRow ();
  int strfun_tmp;
  double sqrt_tmp;
  shortstr scrch1, scrch2, scrch3, scrch4, scrch5, scrch6;
  shortstr buff;
  midstr errTxt;

  scrcolor (Red, BGmain);
  y0 = ycons;
  redres = -1;
  screen_choise = 1;

/*** First Beam Name ***/
  strcpy (scrch1, beam[0].h.name);
again1:
  if (1 == screen_choise) beamlist (redres);
  else                    prtcllist (redres,1);
  if (!blind) {
    redres = input (y0, "Enter 1st Beam: ", scrch1, 1, 70);
  } else {
    redres = KB_ENTER;
    redres = inkey();
  }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_DOWN:
    case KB_F2:
      goto again1;
    case KB_F1:
      show_help ("s_ent_1");
      redres = 0;
      goto again1;
    case KB_F3:
      screen_choise*=-1;
      redres = 0;
      goto again1;
    case KB_ESC:
      clrbox (1, 2, maxCol (), 24);
      return 1;
    }
  trim (scrch1);

  if (strlen (scrch1) > 20)
    {
      warnanykey (12, 12, "Too long the 1st beam name");
      if (blind)
	{
	  sprintf (errTxt, "1st beam: too long the beam name %s", scrch1);
	  batch_error (errTxt, 1);
	}
      goto again1;
    }

  j = locateinbase (scrch1);
  jj = pseudop (j);
  model_beam[0] = 0;
  if (j && !jj)
    {
      model_beam[0] = 1;
      strcpy (beam[0].h.name, scrch1);
      strcpy (beam[0].h.sf_name, "OFF");
      beam[0].h.sf_set = 0;
      beam[0].h.sf_mem = 0;
      beam[0].h.how = 1;
      beam[0].h.parton[0] = j;
      beam[0].h.mass = prtclbase[j - 1].mass;
    }
  else
    {
      bm = In_hadronbase (scrch1);
      if (bm < 0)
	{
	  sprintf (errTxt, "%s is not beam or particle name in the model.", scrch1);
	  warnanykey (12, 12, errTxt);
          if (blind)
	    {
	      sprintf (errTxt, "1st beam: %s is not beam or particle name in the model.", scrch1);
	      batch_error (errTxt, 1);
	    }
	  goto again1;
	}
      else
	{
	  copy_hadrons (&(beam[0].h), hadronbase[bm]);
	}
    }
  if (y0 < MaxR)
    y0++;
  redres = 0;

/*   Enter first beam energy */
  sprintf (scrch2, "%f", beam[0].energy);
again2:
  prtcllist (redres,0);
  if (!blind) {
    redres = input (y0, "Enter 1st Beam Energy (GeV) : ", scrch2, 1, 70);
  } else {
    redres = KB_ENTER;
    redres = inkey();
  }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_DOWN:
    case KB_F2:
    case KB_F3:
      goto again2;
    case KB_F1:
      show_help ("s_ent_3");
      redres = 0;
      goto again2;
    case KB_ESC:
      goto_xy (1, y0);
      clr_eol ();
      redres = 0;
      y0--;
      goto again1;
    }

  if (1 != sscanf (scrch2, "%lf%s", &sqrt_tmp, buff))
    {
      warnanykey(12, 12, "Wrong string for the energy");
      if (blind)
	{
          sprintf(errTxt, "1st beam: wrong string for the beam energy");
	  batch_error (errTxt, 1);
	}
      goto again2;
    }
  if (0 >= sqrt_tmp)
    {
      warnanykey(12, 12, "Energy has to be greater than 0");
      if (blind)
	{
	  sprintf(errTxt, "1st beam: energy has to be greater than 0");
	  batch_error (errTxt, 1);
	}
      goto again2;
    }
  if (beam[0].h.mass > sqrt_tmp)
    {
      warnanykey(12, 12, "The beam particle mass is larger than the beam energy");
      if (blind)
	{
          sprintf(errTxt, "1st beam: beam particle mass is larger than the beam energy");
	  batch_error (errTxt, 1);
	}
      goto again2;
    }
  beam[0].energy = sqrt_tmp;
  if (y0 < MaxR)
    y0++;
  redres = 0;
  screen_choise = 1;

/*** Second Beam Name ***/
  strcpy (scrch3, beam[1].h.name);
again3:
  if (1 == screen_choise) beamlist (redres);
  else                    prtcllist (redres,1);

  if (!blind) {
    redres = input (y0, "Enter 2nd Beam: ", scrch3, 1, 70);
  } else {
    redres = KB_ENTER;
    redres = inkey();
  }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_DOWN:
    case KB_F2:
      goto again3;
    case KB_F1:
      show_help ("s_ent_1");
      redres = 0;
      goto again3;
    case KB_F3:
      screen_choise*=-1;
      redres = 0;
      goto again3;
    case KB_ESC:
      goto_xy (1, y0);
      clr_eol ();
      redres = 0;
      y0--;
      goto again2;
    }
  trim (scrch3);

  if (strlen (scrch3) > 20)
    {
      warnanykey (12, 12, "Too long beam or paricle name");
      if (blind)
	{
	  sprintf (errTxt, "2nd beam: too long beam or particle name %s", scrch3);
	  batch_error (errTxt, 1);
	}
      goto again3;
    }

  j = locateinbase (scrch3);
  jj = pseudop (j);
  model_beam[1] = 0;
  if (j && !jj)
    {
      model_beam[1] = 1;
      strcpy (beam[1].h.name, scrch3);
      strcpy (beam[1].h.sf_name, "OFF");
      beam[1].h.sf_set = 0;
      beam[1].h.sf_mem = 0;
      beam[1].h.how = 1;
      beam[1].h.parton[0] = j;
      beam[1].h.mass = prtclbase[j - 1].mass;
    }
  else
    {
      bm = In_hadronbase (scrch3);
      if (bm < 0)
	{
	  sprintf (errTxt, "%s is not beam or particle name in the model.", scrch3);
	  warnanykey (12, 12, errTxt);
          if (blind)
	    {
	      sprintf (errTxt, "2nd beam: %s is not beam or particle name in the model.", scrch3);
	      batch_error (errTxt, 1);
	    }
	  goto again3;
	}
      else
	{
	  copy_hadrons (&(beam[1].h), hadronbase[bm]);
	}
    }
  if (y0 < MaxR)
    {
      y0++;
    }
  redres = 0;


/*** Enter energy ***/
  sprintf (scrch4, "%f", beam[1].energy);
again4:
  prtcllist (redres,0);
  if (!blind) {
    redres = input (y0, "Enter 2nd Beam Energy (GeV) : ", scrch4, 1, 70);
  } else {
    redres = KB_ENTER;
    redres = inkey();
  }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_DOWN:
    case KB_F2:
    case KB_F3:
      goto again4;
    case KB_F1:
      show_help ("s_ent_3");
      redres = 0;
      goto again4;
    case KB_ESC:
      goto_xy (1, y0);
      clr_eol ();
      redres = 0;
      y0--;
      goto again3;
    }

  if (1 != sscanf (scrch4, "%lf%s", &sqrt_tmp, buff))
    {
      warnanykey(12, 12, "Wrong string for the energy");
      if (blind)
	{
          sprintf(errTxt, "2nd beam: wrong string for the beam energy: %s", scrch4);
	  batch_error (errTxt, 1);
	}
      goto again4;
    }
  if (0 >= sqrt_tmp)
    {
      warnanykey(12, 12, "Energy has to be greater than 0");
      if (blind)
	{
	  sprintf(errTxt, "snd beam: energy has to be greater than 0: %f", sqrt_tmp);
	  batch_error (errTxt, 1);
	}
      goto again4;
    }
  if (beam[1].h.mass > sqrt_tmp)
    {
      warnanykey(12, 12, "The beam particle mass is larger than the beam energy");
      if (blind)
	{
          sprintf(errTxt, "2nd beam: beam particle mass is larger than the beam energy: %f", sqrt_tmp);
	  batch_error (errTxt, 1);
	}
      goto again4;
    }
  beam[1].energy = sqrt_tmp;

  if (y0 < MaxR)
    y0++;
  redres = 0;

/*   Enter structure functions */
  if (0 == strcmp (scrch1, scrch3))
    same_beams = 1;

  if (model_beam[0]) {
      if (blind) { // one brach should be thrown away
        redres = inkey();
      }
  } else {
      i = In_PDFbase_name (beam[0].h.sf_name);
      if (-1 == i)
        i = 0;
      sprintf (scrch5, "%i", strfunbase[i].num);

again5:
      strfunlist (redres);
      if (!blind) {
        if (same_beams)
    	  redres = input (y0, "Enter PDF number : ", scrch5, 1, 70);
    	else
    	  redres = input (y0, "Enter 1st beam PDF number: ", scrch5, 1, 70);
      } else {
        redres = KB_ENTER;
        redres = inkey();
      }
      switch (redres)
    	{
    	case KB_PAGED:
    	case KB_PAGEU:
    	case KB_UP:
    	case KB_DOWN:
    	case KB_F2:
    	case KB_F3:
    	  goto again5;
    	case KB_F1:
    	  show_help ("s_ent_5");
    	  redres = 0;
    	  goto again5;
    	case KB_ESC:
    	  goto_xy (1, y0);
    	  clr_eol ();
    	  redres = 0;
    	  y0--;
    	  goto again4;
    	}

      trim (scrch5);
      if (1 != sscanf (scrch5, "%i%s", &strfun_tmp, buff))
    	{
    	  warnanykey (12, 12, "Wrong string for the PDF number");
    	  if (blind)
    	    {
    	      sprintf (errTxt, "1st beam: wrong string for the PDF number: %s", scrch5);
    	      batch_error (errTxt, 1);
    	    }
    	  goto again5;
    	}
    
      i = In_PDFbase_num (strfun_tmp);
      if (-1 == i)
    	{
    	  warnanykey (12, 11, "Unknown PDF number");
    	  if (blind)
    	    {
	      sprintf (errTxt, "1st beam: unknown PDF number: %s", scrch5);
    	      batch_error (errTxt, 1);
    	    }
    	  goto again5;
    	}
      else
       {
    	  strcpy (beam[0].h.sf_name, strfunbase[i].name);
    	  beam[0].h.sf_set = strfunbase[i].set;
    	  beam[0].h.sf_mem = strfunbase[i].mem;
        }
      if (y0 < MaxR)
      y0++;
    }

  if (same_beams)
    {
      strcpy (beam[1].h.sf_name, beam[0].h.sf_name);
      beam[1].h.sf_set = beam[0].h.sf_set;
      beam[1].h.sf_mem = beam[0].h.sf_mem;
      if (blind) { // one brach should be thrown away
        redres = inkey();
      }
    }
  else
    {
      if (model_beam[1]) {
        if (blind) { // one brach should be thrown away
          redres = inkey();
	}
      } else {
          i = In_PDFbase_name (beam[1].h.sf_name);
          if (-1 == i) {i = 0;}
          sprintf (scrch6, "%i", strfunbase[i].num);
          redres = 0;
again6:
          strfunlist (redres);
          if (!blind) {
	     redres = input (y0, "Enter 2nd beam PDF number: ", scrch6, 1, 70);
          } else {
	    redres = KB_ENTER;
            redres = inkey();
          }
	  switch (redres)
	    {
	    case KB_PAGED:
	    case KB_PAGEU:
	    case KB_UP:
	    case KB_DOWN:
            case KB_F2:
            case KB_F3:
	      goto again6;
	    case KB_F1:
	      show_help ("s_ent_5");
	      redres = 0;
	      goto again6;
	    case KB_ESC:
	      goto_xy (1, y0);
	      clr_eol ();
	      redres = 0;
	      y0--;
	      if (model_beam[0])
	        goto again4;
	      else
	        goto again5;
	     }

          trim (scrch6);
          if (1 != sscanf (scrch6, "%i%s", &strfun_tmp, buff))
	    {
	      warnanykey (12, 11, "Wrong string for the PDF number");
	      if (blind)
	        {
	          sprintf (errTxt, "2nd beam: wrong string for the PDF number: %s", scrch6);
		  batch_error (errTxt, 1);
		}
	      goto again6;
	    }

          i = In_PDFbase_num (strfun_tmp);
          if (-1 == i)
	    {
    	      warnanykey (12, 11, "Unknown PDF number");
    	      if (blind)
    	        {
	          sprintf (errTxt, "2nd beam: unknown PDF number: %s", scrch6);
    	          batch_error (errTxt, 1);
    	        }
	      goto again6;
	    }
          else
	    {
	      strcpy (beam[1].h.sf_name, strfunbase[i].name);
	      beam[1].h.sf_set = strfunbase[i].set;
	      beam[1].h.sf_mem = strfunbase[i].mem;
	    }
      }
    }
  return 0;
}


int
enter_beams (void)
{
  int errorcode;

  errorcode = enter_beam ();	/*   enter beams */
  if (errorcode != 0)
    {
      menuhelp ();
      return 1;			/*  'Esc' pressed  */
    }

  return 0;
}
