/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/parser.h"
#include "service2/include/read_func.h"
#include "service2/include/syst.h"
#include "chep_crt/include/crt_util.h"
#include "out_ext.h"

#include "rd_num.h"
#include "const.h"
#include "phys_val.h"
#include "cut.h"

invcut_ invcut_1[64] = {{0,{0},0,0,0,0.,0.}};

table cutTab =
{"*** Table ***", " Cuts  ",
 "  Parameter  |> Min bound <|> Max bound <|> Exclusive <|", NULL};

int cutnumber = 0;

static void 
cutInit (void)
{
  invcut_1[0].key = 0;
  cutnumber = 0;
}


static int 
addcut (char key, char *lv, int minonw, int maxonw, double cvminw, double cvmaxw, int excl)
{
  int ncut;

  for (ncut = 0; invcut_1[ncut].key; ncut++)
    {
      if (key == invcut_1[ncut].key && !strcmp (lv, invcut_1[ncut].lvinvc))
        {
          if (minonw && (!invcut_1[ncut].minon || cvminw > invcut_1[ncut].cvmin))
            {
              invcut_1[ncut].minon = 1;
              invcut_1[ncut].cvmin = cvminw;
            }
          if (maxonw && (!invcut_1[ncut].maxon || cvmaxw < invcut_1[ncut].cvmax))
            {
              invcut_1[ncut].maxon = 1;
              invcut_1[ncut].cvmax = cvmaxw;
            }
          return 0;
        }
    }
  if (ncut >= 59)
    return 2;

  invcut_1[ncut].key = key;
  strcpy (invcut_1[ncut].lvinvc, lv);
  invcut_1[ncut].minon = minonw;
  invcut_1[ncut].maxon = maxonw;
  invcut_1[ncut].cvmin = cvminw;
  invcut_1[ncut].cvmax = cvmaxw;
  invcut_1[ncut].exclusive = excl;
  invcut_1[ncut + 1].key = 0;
  return 0;
}


int WriteCuts (FILE * nchan) {
  fprintf (nchan, "\n");
  writetable0 (&cutTab, nchan);
  return 0;
}


int ReadCuts (FILE * nchan) {
  fscanf (nchan, "\n");
  return readtable0 (&cutTab, nchan);
}

int get_cutn (void) {
  return cutnumber;
}


int fillCutArray (void) {
  linelist ln = cutTab.strings;
  int lineNum = 0;
  int minOn;
  int maxOn;
  int exclOn;
  char keyChar;
  midstr cutStr;
  midstr minStr;
  midstr maxStr;
  midstr excStr;
  shortstr fieldName;

  cutInit ();
  while (ln != NULL)
    {
      double min_ = 0.0;
      double max_ = 0.0;
      vshortstr lv;
      int i, k;

      cutnumber++;
      cutStr[0] = 0;
      minStr[0] = 0;
      maxStr[0] = 0;
      lineNum++;
      sscanf (ln->line, "%[^|]%*c%[^|]%*c%[^|]%*c%[^|]", cutStr, minStr, maxStr, excStr);

/*============ Parameter ===========*/
      strcpy (fieldName, "Wrong field 'Parameter'");
      if (!checkPhysVal (cutStr, &keyChar, lv))
        goto errorExit;

/*================ MIN bound ============*/
      strcpy (fieldName, "Wrong field 'Min. bound'");
      i = 0;
      while (minStr[i] == ' ')
        i++;
      minOn = (minStr[i] != 0);
      if (minOn && calcExpression (minStr, rd_num, &min_))
        goto errorExit;

/*================== MAX bound ==========*/
      strcpy (fieldName, "Wrong field 'Max bound'");
      i = 0;
      while (maxStr[i] == ' ')
        i++;
      maxOn = (maxStr[i] != 0);
      if (maxOn && calcExpression (maxStr, rd_num, &max_))
        goto errorExit;

/*================== Exclusive condition ==========*/
      strcpy (fieldName, "Wrong field 'Exclusive cond.'");
      i = 0;
      while (excStr[i] == ' ')
        i++;
      exclOn = (excStr[i] != 0);
      if ((!maxOn || !minOn) && exclOn)
        goto errorExit;

/* =========== fill array ========== */

      if (keyChar != 'U')
        SORTARR (lv, strlen (lv));
      if (keyChar == 'A')
        {
          int tmpOn = minOn;
          double min = max_;
          double max = min_;
          keyChar = 'C';
          minOn = maxOn;
          maxOn = tmpOn;
          if (minOn)
            min_ = cos (min * (M_PI) / 180);
          if (maxOn)
            max_ = cos (max * (M_PI) / 180);
        }
      if (keyChar == 'M')
        {
          keyChar = 'S';
          if (minOn)
            min_ = min_ * fabs (min_);
          if (maxOn)
            max_ = max_ * fabs (max_);
        }
      if ((keyChar == 'S'))
        coninv_ (lv);

      if (minOn || maxOn)
        {
          k = addcut (keyChar, lv, minOn, maxOn, min_, max_, exclOn);
          if (k == 2)
            {
              strcpy (fieldName, "Too many cuts ");
              goto errorExit;
            }
        }
      ln = ln->next;
    }
  return 0;

errorExit:
  sprintf (errorText, " Error in Cut table line %d .\n%s.", lineNum, fieldName);
  warnanykey (2, 10, errorText);
  return 1;
}


int rancor (double *vmin, double *vmax, double shift, double fmult, int n) {
  static double vnew;

  if (n < 0)
    return 0;

  if (!invcut_1[n].exclusive) {
    if (invcut_1[n].minon) {
      vnew = (invcut_1[n].cvmin - shift) * fmult;
      if (fmult > 0.)
        *vmin = MAX (*vmin, vnew);
      else
        *vmax = MIN (*vmax, vnew);
    }

    if (invcut_1[n].maxon) {
      vnew = (invcut_1[n].cvmax - shift) * fmult;
      if (fmult > 0.)
        *vmax = MIN (*vmax, vnew);
      else
        *vmin = MAX (*vmin, vnew);
    }
  }

  return 0;
}


double calcCutFactor (void) {
  int i;
  double val;

  for (i = 0; invcut_1[i].key; i++) {
    val = calcPhysVal (invcut_1[i].key, invcut_1[i].lvinvc, "");
    if (invcut_1[i].exclusive) {
      if (invcut_1[i].minon && invcut_1[i].maxon && (val > invcut_1[i].cvmin && val < invcut_1[i].cvmax))
        return 0.;
    } else {
      if (invcut_1[i].minon && (val < invcut_1[i].cvmin))
        return 0.;
      if (invcut_1[i].maxon && (val > invcut_1[i].cvmax))
        return 0.;
    }
  }

  return 1.;
}


int printDetailsOfCutFactor (void){
  int i;
  double val;

  for (i = 0; invcut_1[i].key; i++)
    {
      val = calcPhysVal (invcut_1[i].key, invcut_1[i].lvinvc, "");
      if (invcut_1[i].minon)
        {
          fprintf (stderr, "%c: min = %f, val = %f\n", invcut_1[i].key, invcut_1[i].cvmin, val);
        }
      if (invcut_1[i].maxon)
        {
          fprintf (stderr, "%c: max = %f, val = %f\n", invcut_1[i].key, invcut_1[i].cvmax, val);
        }
    }
  return 0;
}
