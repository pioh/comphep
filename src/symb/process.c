/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* ------------------------------------------------------
*/
#include <math.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/files.h"
#include "service2/include/syst.h"
#include "chep_crt/include/chep_crt.h"

#include "sos.h"
#include "screen.h"
#include "read_mdl.h"
#include "beams.h"
#include "physics.h"
#include "model.h"
#include "m_utils.h"
#include "process_core.h"
#include "process.h"

static int hadrons_are_initialized = 0;
static int nTotP;
static int nFirstP;
static int nTotCP;
static int nFirstCP;


char **
stritems (char *format, char *s)
{
  int i, space = 1, item = 0, l = strlen (s);
  char **res = NULL;

  for (i = 0; i < l; i++)
    if (!strchr (format, s[i]))
      {
        if (space)
          {
            space = 0;
            item++;
          }
      }
    else
      space = 1;

  res = malloc ((item + 1) * sizeof (char *));
  res[item] = NULL;

  item = 0;
  space = 1;
  for (i = 0; i < l; i++)
    if (!strchr (format, s[i]))
      {
        if (space)
          {
            space = 0;
            res[item++] = s + i;
          }
      }
    else
      space = 1;

  return res;
}


int
input (int y0, char *directive, char *text, int cur, int lim)
{
  int err;

  goto_xy (1, y0);
  scrcolor (FGmain, BGmain);
  print ("%s", directive);
  err = str_redact (text, cur, lim);
  return err;
}


static void
addlim (whohow p_list, int j, int k, int type)
{
  int i;

  if (prtclbase[j - 1].anti < j)
    j = prtclbase[j - 1].anti;
  for (i = 0; i < whohowMAX - 1; i++)
    {
      if (p_list[i].who == j && type == 1)
        return;
      if (p_list[i].who == 0)
        {
          p_list[i].who = j;
          p_list[i].how = k;
          p_list[i].type = type;
          p_list[i + 1].who = 0;
          return;
        }
    }
}


int
enter_hadron (int *y, char *nme, int num)
{
  int i;
  int j = 0;
  char **items;
  shortstr errTxt;

/* If the particle is a model one, the hadron is formed by this particle only */
  j = locateinbase (nme);
  if (j && !pseudop (j))
    {
      strcpy (hadrons[num].name, nme);
      strcpy (hadrons[num].sf_name, "OFF");
      hadrons[num].sf_set = 0;
      hadrons[num].sf_mem = 0;
      hadrons[num].parton[0] = j;
      hadrons[num].how = 1;
      return 0;
    }

/* if the particle counsides with previons one, the new hadron is formed by the previous hadron */
  for (i = 0; i < num; i++)
    {
      if (!strcmp (hadrons[i].name, nme))
        {
          strcpy (hadrons[num].name, nme);
          strcpy (hadrons[num].sf_name, hadrons[i].sf_name);
          hadrons[num].sf_set = hadrons[i].sf_set;
          hadrons[num].sf_mem = hadrons[i].sf_mem;
          hadrons[num].how = hadrons[i].how;
          for (j = 0; j < hadrons[num].how; j++)
            {
              hadrons[num].parton[j] = hadrons[i].parton[j];
            }
          return 0;
        }
    }

  if (0 == strlen (nme) || strlen (nme) > 10)
    {
      fprintf (stderr, "CompHEP: internal error %i. Please, inform authors.",
               2);
      return -1;
    }

/* if composite particles DB contains the particle with name the hadron is formed by this cparticle */
  if (0 > In_composite_particle_base (nme))
    {
      int stop = 1;
      do
        {
          sprintf (errTxt,
                   "The model does not contain the composite particle %s\n"
                   "Do you want to introduce new particle?", nme);
          if (mess_y_n (12, 11, errTxt))
            {
              int n_model;
              do
                {
                  if (blind)
                    {
                      sprintf (errTxt, "Unknown composite particle %s", nme);
                      batch_error (errTxt, 1);
                    }
                  edittable (1, 1, &modelTab[4], 1, "s_mdl_5", 0);
                  n_model = getModelNumberSymb ();
                  writeModelFiles (n_model, "models");
                  if (loadModel (TRUE))
                    stop = 0;
                }
              while (stop);
            }
          else
            {
              stop = 0;
            }
        }
      while (stop);
    }

  j = In_composite_particle_base (nme);
  items = stritems (" ,", cpartbase[j].cpart);
  for (i = 0, hadrons[num].how = 0; items[i]; i++)
    {
      shortstr pname;
      sscanf (items[i], "%[^ ,]", pname);
      j = locateinbase (pname);
      if (!j || pseudop (j))
        {
          if (blind)
            {
              sprintf (errTxt,
                       "Unknown constituent %s in composite particle %s",
                       pname, nme);
              batch_error (errTxt, 1);
            }
          warnanykey (12, 12, "Unknown model particle");
          return 1;
        }
      hadrons[num].parton[hadrons[num].how++] = j;
    }
  strcpy (hadrons[num].name, nme);
  strcpy (hadrons[num].sf_name, "OFF");
  hadrons[num].sf_set = 0;
  hadrons[num].sf_mem = 0;

  return 0;
}

void
prtcllist (int key, int switch_list)
{
  int i, j, pnum;
  int tabMax, tabSz;

  shortstr buf1, buf2, buf3, buf4, hlp, p1, p2, fullname;
  longstr str;
  linelist ln;

  tabMax = ycons - 7;
  if (-1 == key)
    {
      key = 0;
      for (i = 2; i < 24; i++)
        {
          goto_xy (1, i);
          clr_eol ();
        }
    }

  if (key == 0)
    {
      scrcolor (FGmain, BGmain);
      for (i = 2; i <= 18; i++)
        {
          goto_xy (1, i);
          clr_eol ();
        }
/* switch_list = 1(in beam menu), 2(in decay menu), 0(other cases) */
      switch (switch_list)
        {
        case 1:
          {
            goto_xy (10, 3);
            scrcolor (Blue, BGmain);
            print ("List of (anti)particles (switch to the beam list by F3)");
            break;
          }
        case 2:
          {
            goto_xy (5, 3);
            scrcolor (Blue, BGmain);
            print
              ("List of (anti)particles (switch to the composite particle list by F3)");
            break;
          }
        case 0:
          {
            goto_xy (24, 3);
            scrcolor (Blue, BGmain);
            print ("List of (anti)particles");
            break;
          }
        }
      nTotP = 0;
      nFirstP = 1;
    }
  else
    {
      if (nTotP <= 3 * tabMax)
        return;
      switch (key)
        {
        case KB_DOWN:
          nFirstP += 3;
          break;
        case KB_UP:
          nFirstP -= 3;
          break;
        case KB_PAGED:
          nFirstP += 3 * tabMax;
          break;
        case KB_PAGEU:
          nFirstP -= 3 * tabMax;
          break;
        }
      if (nFirstP < 1)
        nFirstP = 1;
      if (nTotP - nFirstP + 3 < 3 * tabMax)
        nFirstP = 1 + 3 * ((nTotP + 2) / 3) - 3 * tabMax;
      clrbox (1, 4, 79, 5 + tabMax);
    }
  goto_xy (3, 5);
  i = 0;
  ln = prtcls_tab.strings;
  scrcolor (FGmain, BGmain);
  while (ln != NULL)
    {
      sscanf (ln->line,
              "%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]",
              fullname, p1, p2, buf1, buf2, buf3, buf4, hlp);
      trim (p1);
      trim (p2);
      trim (fullname);
      pnum = locateinbase (p1);
      trim (hlp);
      if (prtclbase[pnum - 1].top != NULL && strcmp (hlp, "*") != 0)
        {
          i++;
          if (i >= nFirstP && (i - nFirstP) / 3 < tabMax)
            {
              switch (strlen (p1))
                {
                case 1:
                  sprintf (str, "%s(%s)     %s", p1, p2, fullname);
                  break;
                case 2:
                  sprintf (str, "%s(%s)   %s", p1, p2, fullname);
                  break;
                case 3:
                  sprintf (str, "%s(%s) %s", p1, p2, fullname);
                  break;
                default:
                  fprintf (stderr,
                           "***** prtcllist: Too long name p1=%s. ABORT\n",
                           p1);
                  exit (99);
                }
              print ("%s", str);
              str[0] = '\0';
              j = i % 3;
              if (j == 0)
                goto_xy (3, where_y () + 1);
              else
                goto_xy (3 + 26 * j, where_y ());
            }
        }
      ln = ln->next;
    }
  nTotP = i;
  tabSz = MIN ((nTotP + 2) / 3, tabMax);
  chepbox (1, 4, 79, 5 + tabSz);

  if (nFirstP > 1)
    {
      goto_xy (72, 4);
      print ("PgUp");
    }
  if (nFirstP + 3 * tabSz <= nTotP)
    {
      goto_xy (72, 5 + tabMax);
      print ("PgDn");
    }
  scrcolor (FGmain, BGmain);
}


void
composite_prtcllist (int key)
{
  int i;
  int tabMax, tabSz;

  tabMax = ycons - 7;
  if (-1 == key)
    {
      key = 0;
      for (i = 2; i < 24; i++)
        {
          goto_xy (1, i);
          clr_eol ();
        }
    }

  if (key == 0)
    {
      scrcolor (FGmain, BGmain);
      for (i = 2; i <= 18; i++)
        {
          goto_xy (1, i);
          clr_eol ();
        }
      goto_xy (10, 3);
      scrcolor (Blue, BGmain);
      print
        ("List of composite particles (switch to the particle list by F3)");
      nTotCP = 0;
      nFirstCP = 0;
    }
  else
    {
      if (nTotCP <= tabMax)
        return;
      switch (key)
        {
        case KB_DOWN:
          nFirstCP++;
          break;
        case KB_UP:
          nFirstCP--;
          break;
        case KB_PAGED:
          nFirstCP += tabMax;
          break;
        case KB_PAGEU:
          nFirstCP -= tabMax;
          break;
        }
      if (nFirstCP < 0)
        nFirstCP = 0;
      if (nTotCP - nFirstCP < tabMax)
        nFirstCP = nTotCP - tabMax;
      clrbox (1, 4, 79, 5 + tabMax);
    }
  goto_xy (3, 5);
  scrcolor (FGmain, BGmain);
  for (i = 0; i < n_cpart; i++)
    {
      if (i >= nFirstCP && (i - nFirstCP) < tabMax)
        {
          print ("Name: %s (%s)", cpartbase[i].name, cpartbase[i].cpart);
          goto_xy (3, where_y () + 1);
        }
    }
  nTotCP = n_cpart;
  tabSz = MIN (nTotCP + 2, tabMax);
  chepbox (1, 4, 79, 5 + tabSz);

  if (nFirstCP > 0)
    {
      goto_xy (72, 4);
      print ("PgUp");
    }
  if (nFirstCP + tabSz <= nTotCP)
    {
      goto_xy (72, 5 + tabMax);
      print ("PgDn");
    }
  scrcolor (FGmain, BGmain);
}


int
exclude_diagrams (int *y0)
{
  int i;
  int jj, j;
  int k, ntot;
  int redres;
  int type = 1;
  char **items;
  char *n;
  shortstr frgm;
  shortstr errTxt;

  if (*y0 < maxRow ())
    (*y0)++;
  redres = 0;

again:
  prtcllist (redres, 0);
  if (!blind)
    {
      shortstr tmpexcl_str;
      strcpy (tmpexcl_str, getExclprtlist ());
      redres = input (*y0, "Exclude diagrams with ", tmpexcl_str, 1, 70);
      setExclprtlist (tmpexcl_str);
    }
  else
    {
      redres = KB_ENTER;
      redres = inkey ();
    }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_DOWN:
    case KB_F2:
    case KB_F3:
      goto again;
    case KB_F1:
      show_help ("s_ent_4");
      redres = 0;
      goto again;
    case KB_ESC:
      clrbox (1, 2, maxCol (), 24);
      (*y0) = (*y0) - 2;
      return 0;
    }

  nilprtcl (exclude_list);
  ntot = 0;
  items = stritems (",", getExclprtlist ());
  for (i = 0; items[i]; i++)
    {
      sscanf (items[i], "%[^,]", frgm); /*  read of condition for excluding */

/*
type defines logical connective for excuding/keeping:
  type =  2 -> exclude/keep;
  type =  1 -> operation ">";
  type = -1 -> operation "!=";
  type = -2 -> operation "=";
  type = -3 -> operation "<";
*/
      n = strchr (frgm, '!');
      if (n)
        type = -1;
      else
        {
          n = strchr (frgm, '=');
          if (n)
            type = -2;
          else
            {
              n = strchr (frgm, '<');
              if (n)
                type = -3;
              else
                {
                  n = strchr (frgm, '>');
                  if (n)
                    type = 1;
                  else
                    type = 2;
                }
            }
        }

      if (n)
        {
          n[0] = 0;
          k = 0;
          sscanf (n + 1, "%d", &k);
          if (k < 0)
            {
              warnanykey (12, 12, "Wrong number");
              if (blind)
                {
                  sprintf (errTxt, "Diagram excluding: wrong number %i", k);
                  batch_error (errTxt, 1);
                }
              goto again;
            }
        }
      else
        k = 1;
      trim (frgm);

      j = locateinbase (frgm);
      jj = pseudop (j);
      if (!j || jj)
        {
          warnanykey (12, 12, "Unknown model particle");
          if (blind)
            {
              sprintf (errTxt,
                       "Unknown model particle %s in the excluding list",
                       frgm);
              batch_error (errTxt, 1);
            }
          goto again;
        }
      ntot++;
      if (ntot >= whohowMAX)
        {
          warnanykey (12, 12, "Too many items");
          if (blind)
            {
              batch_error ("Too many particles in the excluding list", 1);
            }
          goto again;
        }
      addlim (exclude_list, j, k, type);
    }
  free (items);
  return 1;
}


int
keep_diagrams (int *y0)
{
  int i;
  int jj, j;
  int k, ntot;
  int redres;
  int type = 1;
  char *n;
  char **items;
  shortstr errTxt;
  shortstr frgm;

  if (*y0 < maxRow ())
    (*y0)++;
  redres = 0;

again:
  prtcllist (redres, 0);
  if (!blind)
    {
      shortstr tmpkeep_str;
      strcpy (tmpkeep_str, getKeepprtlist ());
      redres = input (*y0, "Keep diagrams with ", tmpkeep_str, 1, 70);
      setKeepprtlist (tmpkeep_str);
    }
  else
    {
      redres = KB_ENTER;
      redres = inkey ();
    }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_DOWN:
    case KB_F2:
    case KB_F3:
      goto again;
    case KB_F1:
      show_help ("s_ent_5");
      redres = 0;
      goto again;
    case KB_ESC:
      clrbox (1, 2, maxCol (), 24);
      (*y0) = (*y0) - 2;
      return 0;
    }

  nilprtcl (keep_list);
  ntot = 0;
  items = stritems (",", getKeepprtlist ());
  for (i = 0; items[i]; i++)
    {
      sscanf (items[i], "%[^,]", frgm); /*  read condition for keeping */
      n = strchr (frgm, '!');
      if (n)
        type = -1;
      else
        {
          n = strchr (frgm, '=');
          if (n)
            type = -2;
          else
            {
              n = strchr (frgm, '<');
              if (n)
                type = -3;
              else
                {
                  n = strchr (frgm, '>');
                  if (n)
                    type = 1;
                  else
                    type = 2;
                }
            }
        }

      if (n)
        {
          n[0] = 0;
          k = 0;
          sscanf (n + 1, "%d", &k);
          if (k < 0)
            {
              warnanykey (12, 12, "wrong number");
              if (blind)
                {
                  sprintf (errTxt, "Diagram keeping: wrong number: %i", k);
                  batch_error (errTxt, 1);
                }
              goto again;
            }
        }
      else
        k = 1;
      trim (frgm);

      j = locateinbase (frgm);
      jj = pseudop (j);
      if (!j || jj)
        {
          warnanykey (12, 12, "Unknown model particle");
          if (blind)
            {
              sprintf (errTxt, "Unknown model particle %s in keep list",
                       frgm);
              batch_error (errTxt, 1);
            }
          goto again;
        }
      ntot++;
      if (ntot >= whohowMAX)
        {
          warnanykey (12, 12, "Too many items");
          if (blind)
            {
              batch_error ("Too many particles in the keep list", 1);
            }
          goto again;
        }
      addlim (keep_list, j, k, type);
    }
  free (items);
  return 1;
}


static int
enter_decay_particle (void)
{
  int i, y0;
  int redres;
  int screen_choice;
  double m, m0;
  shortstr scrch;
  shortstr name;
  shortstr errTxt;

  y0 = ycons;
  if (!hadrons_are_initialized)
    {
      hadrons_are_initialized = 1;
      for (i = 0; i < MAXINOUT; i++)
        hadrons[i].name[0] = 0;
    }
  y0 = ycons;
  redres = -1;
  screen_choice = 1;

  strcpy (scrch, beam[0].h.name);
/* artificial change of the default name e -> Z */
  if (0 == strcmp (scrch, "e"))
    strcpy (scrch, "Z");
again:
  if (1 == screen_choice)
    prtcllist (redres, 2);
  else
    composite_prtcllist (redres);

  if (!blind)
    {
      redres = input (y0, "Enter decayed particle: ", scrch, 1, 70);
    }
  else
    {
      redres = KB_ENTER;
      redres = inkey ();
    }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_DOWN:
    case KB_F2:
      goto again;
    case KB_F3:
      screen_choice *= -1;
      redres = 0;
      goto again;
    case KB_F1:
      show_help ("s_ent_1");
      redres = 0;
      goto again;
    case KB_ESC:
      clrbox (1, 2, maxCol (), 24);
      return 0;
    }

  trim (scrch);

/* Check of the final particles (from the final state string)*/
  sscanf (scrch, "%[^, ]", name);
  if (strlen (name) > 6)
    {
      warnanykey (12, 12, "Too long name");
      if (blind)
        {
          sprintf (errTxt, "Initial state: too long particle name %s", name);
          batch_error (errTxt, 1);
        }
      goto again;
    }
  if (strlen (name) == 3
      && name[1] == '*'
      && (name[2] == 'X' || name[2] == 'x') && isdigit (name[0]))
    {
/* X-particles... */

      int locn_x = getnx () + name[0] - '0';
      int locnout = getnout () + name[0] - '0';
      setnout (locnout);
      setnx (locn_x);
    }
  else
    {
/* Ordinary particles */
      if (enter_hadron (&y0, name, 0))
        {
          warnanykey (12, 12, "Error in final state");
          goto again;
        }
    }

  copy_hadrons (&(beam[0].h), hadrons[0]);
  copy_hadrons (&(beam[1].h), hadrons[0]);
  m0 = prtclbase[hadrons[0].parton[0] - 1].mass;
  for (i = 1; i < hadrons[0].how; i++)
    {
      m = prtclbase[hadrons[0].parton[i] - 1].mass;
      if (m < m0)
        m0 = m;
    }
  beam[0].energy = m0;

  return y0++;
}


static int
enter_final_state (int *y0)
{
  int i, y00;
  int redres;
  int curh;
  int screen_choice;
  double inMass, outMass, m, m0;
  double e[2], p[2];
  shortstr name;
  shortstr instatech;
  shortstr errTxt;
  int locn_x;
  int locnout;

  char **items = NULL;

  scrcolor (FGmain, BGmain);
  for (i = *y0 + 1; i < 25; i++)
    {
      goto_xy (1, i);
      clr_eol ();
    }
  scrcolor (Red, BGmain);
  y00 = *y0;
  if (*y0 < maxRow ())
    {
      (*y0)++;
      y00++;
    }

  if (!hadrons_are_initialized)
    {
      hadrons_are_initialized = 1;
      for (i = 0; i < MAXINOUT; i++)
        {
          hadrons[i].name[0] = 0;
          hadrons[i].sf_name[0] = 0;
          hadrons[i].sf_set = 0;
          hadrons[i].sf_mem = 0;
        }
    }
  screen_choice = 1;
  redres = 0;

again:
  if (1 == screen_choice)
    prtcllist (redres, 2);
  else
    composite_prtcllist (redres);

  if (2 == getnin ())
    sprintf (instatech, "Enter Final State: %s,%s -> ", beam[0].h.name,
             beam[1].h.name);
  else
    sprintf (instatech, "Enter  Final State: %s -> ", beam[0].h.name);

  if (!blind)
    {
      shortstr tmpfinalstatech;
      strcpy (tmpfinalstatech, getFinalstatech ());
      redres = input (y00, instatech, tmpfinalstatech, 1, 70);
      setFinalstatech (tmpfinalstatech);
    }
  else
    {
      redres = KB_ENTER;
      redres = inkey ();
    }
  switch (redres)
    {
    case KB_PAGED:
    case KB_PAGEU:
    case KB_UP:
    case KB_F2:
    case KB_DOWN:
      goto again;
    case KB_F1:
      show_help ("s_ent_1");
      goto again;
    case KB_F3:
      screen_choice *= -1;
      redres = 0;
      goto again;
    case KB_ESC:
      clrbox (1, 2, maxCol (), 24);
      (*y0) = (*y0) - 2;
      return 0;
    }

  if (!items)
    free (items);
  items = stritems (" ,", getFinalstatech ());

  curh = 2;
  if (1 == getnin ())
    {
      curh = 1;
      setsqrtS (0.0);
      setRapidity (0.0);
    }

/* Check of the final particles (from the final state string)*/
  locn_x = 0;
  locnout = 0;
  setnout (locnout);
  setnx (locn_x);
  for (i = 0; items[i]; i++)
    {
      sscanf (items[i], "%[^, ]", name);
      if (strlen (name) > 7)
        {
          warnanykey (12, 12, "Too long name");
          if (blind)
            {
              sprintf (errTxt, "Final state: too long particle name %s",
                       name);
              batch_error (errTxt, 1);
            }
          goto again;
        }

      if (strlen (name) == 3
          && name[1] == '*'
          && (name[2] == 'X' || name[2] == 'x') && isdigit (name[0]))
        {
/* X-particles... */
          locn_x += name[0] - '0';
          locnout += name[0] - '0';
          setnout (locnout);
          setnx (locn_x);
        }
      else
        {
/* Ordinary particles */
          if (curh == MAXINOUT - 1)
            {
              warnanykey (12, 12, "Too many particles in the process");
              if (blind)
                {
                  sprintf (errTxt,
                           "Final state: too many particles in the process");
                  batch_error (errTxt, 1);
                }
              goto again;
            }
          if (enter_hadron (y0, name, curh))
            {
              warnanykey (12, 12, "Error in final state");
              (*y0) = y00;
              goto again;
            }
          locnout++;
          setnout (locnout);
          curh++;
        }
    }
  free (items);

  for (i = 0; i < getnin (); i++)
    {
      copy_hadrons (&(hadrons[i]), beam[i].h);
      e[i] = beam[i].energy;
      p[i] = sqrt (e[i] * e[i] - hadrons[i].mass * hadrons[i].mass);
    }

  if (getnout () < 2)
    {
      warnanykey (12, 12,
                  "The number of outgoing particles has to be greater than 1");
      if (blind)
        {
          sprintf (errTxt,
                   "Final state: the number of outgoing particles has to be greater than 1");
          batch_error (errTxt, 1);
        }
      (*y0) = y00;
      goto again;
    }

  for (curh = 0, inMass = 0, outMass = 0; curh < getntot () - getnx ();
       curh++)
    {
      m0 = prtclbase[hadrons[curh].parton[0] - 1].mass;
      for (i = 1; i < hadrons[curh].how; i++)
        {
          m = prtclbase[hadrons[curh].parton[i] - 1].mass;
          if (m < m0)
            m0 = m;
        }
      if (curh < getnin ())
        inMass += m0;
      else
        outMass += m0;
    }

  if (2 == getnin ())
    {
      double locrapidity;
      double locsqrts =
        sqrt ((e[0] + e[1]) * (e[0] + e[1]) - (p[0] - p[1]) * (p[0] - p[1]));
      setsqrtS (locsqrts);
      if (locsqrts < 0.0)
        {
          locsqrts = hadrons[0].mass + hadrons[1].mass;
          setsqrtS (locsqrts);
        }
      locrapidity = atanh ((p[0] - p[1]) / (e[0] + e[1]));
      setRapidity (locrapidity);
      if (locsqrts <= inMass || locsqrts <= outMass)
        {
          warnanykey (12, 12, "Energy is not enough");
          if (blind)
            {
              sprintf (errTxt, "Final state: Energy is not enough");
              batch_error (errTxt, 1);
            }
          (*y0) = y00;
          goto again;
        }
    }
  else
    {
      if (inMass <= outMass)
        {
          warnanykey (12, 12, "This decay mode is forbidden");
          if (blind)
            {
              sprintf (errTxt, "Final state: this decay mode is forbidden");
              batch_error (errTxt, 1);
            }
          (*y0) = y00;
          goto again;
        }
    }
  return 1;
}


static int
enter_scattering_process (void)
{
  shortstr tmps;
  int y0 = ycons - 1;

  setnin (2);

again:
  if (!enter_final_state (&y0))
    return 1;

  if (!exclude_diagrams (&y0))
    {
      y0 = ycons - 1;
      goto again;
    }

  if (!keep_diagrams (&y0))
    {
      y0 = ycons - 1;
      goto again;
    }

  sprintf (tmps, "%s,%s -> %s", hadrons[0].name, hadrons[1].name,
           getFinalstatech ());
  setProcessch (tmps);
  return 0;
}


static int
enter_decay_process (void)
{
  shortstr tmps;
  int y0;

  setnin (1);

again0:
  y0 = enter_decay_particle ();
  if (!y0)
    return 1;

again1:
  if (!enter_final_state (&y0))
    goto again0;                /*  'Esc' pressed  */

  if (!exclude_diagrams (&y0))
    {
      y0 = ycons;
      goto again1;              /*  'Esc' pressed  */
    }

  if (!keep_diagrams (&y0))
    {
      y0 = ycons;
      goto again1;              /*  'Esc' pressed  */
    }

  sprintf (tmps, "%s -> %s", hadrons[0].name, getFinalstatech ());
  setProcessch (tmps);
  return 0;
}


int
enter_process (int ProcessChoice)
{
  int res = 1;

  if (1 == ProcessChoice)
    {
      res = enter_decay_process ();
    }
  else if (2 == ProcessChoice)
    {
    again:
      res = enter_beams ();
      if (0 != res)
        {
          menuhelp ();
          return res;           /*  'Esc' pressed  */
        }
      res = enter_scattering_process ();
      if (0 != res)
        {
          menuhelp ();
          goto again;           /*  'Esc' pressed  */
        }
    }
  return res;
}
