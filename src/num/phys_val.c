/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/

#include <math.h>
#include <ctype.h>

#include "service2/include/chep_limits.h"
#include "service2/include/4_vector.h"

#include "out_ext.h"

#include "LesHouches.h"
#include "subproc.h"
#include "kinaux.h"
#include "userFun.h"
#include "const.h"
#include "phys_val.h"

double 
calcPhysVal (char key, char *lv, char *restframe)
{
  int i, j;
  int ntot_ = nin_ + nout_;
  int np1 = 4 * (lv[0] - 1);
  int np2 = 4 * (lv[1] - 1);

  double p1, p2, p3, q1, q2, q3, mp, mq, cs, dl, ep, eq;
  double nu1,nu2;
  double s = 0;
  double pp[4] = {0., 0., 0., 0.};
  double newrf[4] = {0., 0., 0., 0.};
  double lab[4] = {0., 0., 0., 0.};
  double boost[4] = {0, 0, 0, 0};
  double pvect_boosted[400]={0};

/* ----- New Rest Frame -----------*/
  if ( 0 != restframe[0]) {
    i = 0;
    while (restframe[i]) {
      for (j = 0; j < 4; ++j)
        newrf[j] += pvect[4 * (restframe[i] - 1) + j];
      ++i;
    }

    for(i = 0; i < ntot_; ++i) {
      for (j = 0; j < 4; ++j)
        lab[j] = pvect[4 * i + j];
      lorenc (lab, newrf, boost);
      for (j = 0; j < 4; ++j)
       pvect_boosted[4 * i + j] = boost[j];
    }
  } else {
    for (i = 0; i < ntot_; ++i)
      for (j = 0; j < 4; ++j)
       pvect_boosted[4 * i + j] = pvect[4 * i + j];
  }
  i = 0;
/*---------------------------------*/

  switch (key)
    {
    case 'A':
    case 'C':
      p1 = pvect_boosted[np1 + 1];
      p2 = pvect_boosted[np1 + 2];
      p3 = pvect_boosted[np1 + 3];

      q1 = pvect_boosted[np2 + 1];
      q2 = pvect_boosted[np2 + 2];
      q3 = pvect_boosted[np2 + 3];

      cs = (p1 * q1 + p2 * q2 + p3 * q3) /
	sqrt ((p1 * p1 + p2 * p2 + p3 * p3) * (q1 * q1 + q2 * q2 + q3 * q3));
      if (key == 'A')
	return acos (cs) * 180 / M_PI;
      return cs;
    case 'B':
      ep = pvect_boosted[np1 + 0];
      p1 = pvect_boosted[np1 + 1];
      p2 = pvect_boosted[np1 + 2];
      p3 = pvect_boosted[np1 + 3];

      eq = pvect_boosted[np2 + 0];
      q1 = pvect_boosted[np2 + 1];
      q2 = pvect_boosted[np2 + 2];
      q3 = pvect_boosted[np2 + 3];

      return (p1 * q1 + p2 * q2 + p3 * q3) / ep / sqrt (q1 * q1 + q2 * q2 + q3 * q3);
    case 'D':
      ep = pvect_boosted[np1 + 0];
      p1 = pvect_boosted[np1 + 1];
      p2 = pvect_boosted[np1 + 2];
      p3 = pvect_boosted[np1 + 3];

      eq = pvect_boosted[np2 + 0];
      q1 = pvect_boosted[np2 + 1];
      q2 = pvect_boosted[np2 + 2];
      q3 = pvect_boosted[np2 + 3];

      return (p1 * q1 + p2 * q2 + p3 * q3) / eq / ep;
    case 'J':
      p1 = pvect_boosted[np1 + 1];
      p2 = pvect_boosted[np1 + 2];
      p3 = pvect_boosted[np1 + 3];
      mp = sqrt (p1 * p1 + p2 * p2 + p3 * p3);

      q1 = pvect_boosted[np2 + 1];
      q2 = pvect_boosted[np2 + 2];
      q3 = pvect_boosted[np2 + 3];
      mq = sqrt (q1 * q1 + q2 * q2 + q3 * q3);

      cs = (p1 * q1 + p2 * q2) / sqrt ((p1 * p1 + p2 * p2) * (q1 * q1 + q2 * q2));
      cs = acos (cs);

      dl = (mp + p3) * (mq - q3) / (mp - p3) / (mq + q3);
      dl = 0.5 * log (dl);

      return sqrt (dl * dl + cs * cs);

    case 'P':
      {
	double mtot, mtot2, ms, md, pcm, p;

	pinf_ (proces_1.nsub, lv[0], NULL, &mp);
	pinf_ (proces_1.nsub, lv[1], NULL, &mq);

	for (j = 0; j < 4; j++)
	  pp[j] += pvect_boosted[np1 + j] + pvect_boosted[np2 + j];

	p = pp[1] * pp[1] + pp[2] * pp[2] + pp[3] * pp[3];
	mtot2 = pp[0] * pp[0] - p;
	mtot = sqrt (mtot2);

	ms = mp + mq;
	md = mp - mq;

	pcm = sqrt ((mtot2 - ms * ms) * (mtot2 - md * md)) / (2 * mtot);
	s = pp[1] * pvect_boosted[np1 + 1] + pp[2] * pvect_boosted[np1 + 2] + pp[3] * pvect_boosted[np1 + 3];

	return (s * pp[0] - pvect_boosted[np1] * p) / (sqrt (p) * mtot * pcm);

      }
    case 'E':
      {
       while (lv[i] != 0)
       	s += pvect_boosted[(lv[i++] << 2) - 4];
      return s;
      } 
    case 'T':
      {
	i = 0;
	do
	  for (j = 1; j < 3; j++)
	    pp[j] += pvect_boosted[4 * (lv[i] - 1) + j];
	while (lv[++i]);
	return sqrt (pp[1] * pp[1] + pp[2] * pp[2]);
      }
    case 'S':
    case 'M':
      {
	do
	  {
	    if (lv[i] > nin_)
	      for (j = 0; j < 4; j++)
		pp[j] += pvect_boosted[4 * (lv[i] - 1) + j];
	    else
	      for (j = 0; j < 4; j++)
		pp[j] -= pvect_boosted[4 * (lv[i] - 1) + j];
	  }
	while (lv[++i]);
	s = pp[0] * pp[0];
	for (j = 1; j < 4; j++)
	  s -= pp[j] * pp[j];
	if (key == 'M')
	  return sqrt (s);
	return s;
      }
    case 'Y':
      do
	for (j = 0; j < 4; j += 3)
	  pp[j] += pvect_boosted[4 * (lv[i] - 1) + j];
      while (lv[++i]);
      return log ((pp[0] + pp[3]) / (pp[0] - pp[3])) / 2;
    case 'N':
      do
	for (j = 0; j < 4; j++)
	  pp[j] += pvect_boosted[4 * (lv[i] - 1) + j];
      while (lv[++i]);
      mp = sqrt (pp[1] * pp[1] + pp[2] * pp[2] + pp[3] * pp[3]);
      return log ((mp + pp[3]) / (mp - pp[3])) / 2;


    case 'X':
   	for (j=0;j<4;j++) pp[j] += pvect_boosted[np1 + j];
        mp = sqrt (pp[1] * pp[1] + pp[2] * pp[2] + pp[3] * pp[3]);
        nu1= log ((mp + pp[3]) / (mp - pp[3])) /2;
      
        for (j=0;j<4;j++) pp[j] += pvect_boosted[np2 + j];
        mp = sqrt (pp[1] * pp[1] + pp[2] * pp[2] + pp[3] * pp[3]);
        nu2= log ((mp + pp[3]) / (mp - pp[3])) /2;

        return (nu2-nu1);


      case 'U':
       return userFunction (lv); 
 

    }
  return 0;
}

int 
checkPhysVal (char *name, char *key, char *plist)
{
  int i = 0;
  int j = 0;
  int n, k;

  while (name[i] == ' ' && name[i] != 0)
    i++;
  *key = name[i++];

  if (*key == 0)
    return 0;
  *key = toupper (*key);
  if (strchr ("ABDCEJMPSTUYNX", *key) == NULL)
    return 0;


  if (*key == 'U')
    {
      for (; name[i] && name[i] != ' ' && i < 6; i++)
	plist[j++] = name[i];
      plist[j] = 0;
      for (; name[i]; i++)
	if (name[i] != ' ')
	  return 0;
      return 1;
    }


  for (; name[i] && name[i] != ' '; i++)
    {
      n = name[i] - '0';
      if (n <= 0 || n > nin_ + nout_)
	return 0;
      for (k = 0; k < j; k++)
	if (plist[k] == n)
	  return 0;
      plist[j++] = n;
    }
  plist[j] = 0;
  for (; name[i]; i++)
    if (name[i] != ' ')
      return 0;


  if (strchr ("BCDAJP", *key) != NULL && strlen (plist) != 2)
    return 0;

  if (strchr ("MS", *key) != NULL && strlen (plist) < 2)
    return 0;

  if (strchr ("JPT", *key) != NULL)
    for (i = 0; i < strlen (plist); i++)
      {
	if (plist[i] <= nin_)
	  return 0;
      }

  if (strchr ("MNYX", *key) && !spole_ (plist))
    return 0;
/*
  if (nin_ == 1)
    {
      if (strchr ("TYN", *key))
	return 0;
      if (strchr ("ACP", *key) && (plist[0] == 1 || plist[1] == 1))
	return 0;
    }
*/
  return 1;
}

void 
xName (char key, char *plist, char *xname, char *units)
{
  int i;

  switch (key)
    {
    case 'A':
      sprintf (units, "Deg");
      sprintf (xname, "Angle(p%d,p%d)", (int) labs (plist[0]), (int) labs (plist[1]));
      break;
    case 'C':
      strcpy (units, "");
      sprintf (xname, "Cosine(p%d,p%d)", (int) labs (plist[0]), (int) labs (plist[1]));
      break;
    case 'E':
      sprintf (units, "GeV");
      sprintf (xname, "Energy E%d", (int) labs (plist[0]));
      i = 1;
      while (plist[i])
	sprintf (xname + strlen (xname), "+E%d", (int) labs (plist[i++]));
      break;
    case 'J':
      strcpy (units, "");
      sprintf (xname, "J(p%d,p%d)", (int) labs (plist[0]), (int) labs (plist[1]));
      break;
    case 'M':
      sprintf (units, "GeV");
      sprintf (xname, "Mass{p%d", (int) labs (plist[0]));
      i = 1;
      while (plist[i])
	sprintf (xname + strlen (xname), "+p%d", (int) labs (plist[i++]));
      strcat (xname, "}");
      break;
    case 'P':
      strcpy (units, "");
      sprintf (xname, "S.M.Cosine(p%d,p%d)", (int) labs (plist[0]), (int) labs (plist[1]));
      break;
    case 'S':
      sprintf (units, "GeV^2");
      sprintf (xname, "(p%d", plist[0]);
      i = 1;
      while (plist[i])
	sprintf (xname + strlen (xname), "+p%d", (int) labs (plist[i++]));
      strcat (xname, ")^2");
      break;
    case 'T':
      sprintf (units, "GeV");
      sprintf (xname, "Transverse momentum Pt%d", (int) labs (plist[0]));
      i = 1;
      while (plist[i])
	sprintf (xname + strlen (xname), "+Pt%d", (int) labs (plist[i++]));
      break;
    case 'Y':
      strcpy (units, "");
      sprintf (xname, "Rapidity_");
      for (i = 0; plist[i]; i++)
	sprintf (xname + strlen (xname), "%d", (int) labs (plist[i]));
      break;
    case 'N':
      strcpy (units, "");
      sprintf (xname, "pseudo-rapidity_");
      for (i = 0; plist[i]; i++)
	sprintf (xname + strlen (xname), "%d", (int) labs (plist[i]));
      break;
   case 'X':
      strcpy (units, "");
      sprintf (xname, "delta pseudo-rapidity_");
      for (i = 0; plist[i]; i++)
	sprintf (xname + strlen (xname), "%d", (int) labs (plist[i]));
      break;
    case 'U':
      sprintf (units, "?");
      strcpy (xname, plist);
      break;
    }
}
