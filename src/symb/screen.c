/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* ------------------------------------------------------
*/
#include <unistd.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"

#include "m_utils.h"
#include "sos.h"
#include "showgrph.h"
#include "r_code.h"
#include "read_mdl.h"
#include "physics.h"
#include "screen.h"
#include "prepdiag.h"
#include "process.h"
#include "process_core.h"
#include "cweight.h"
int menulevel = 0;

#define tcol Green
#define mpos 7
#define graphpos 8


void
diag_stat (int type, int *n_sub, int *n_del, int *n_calc, int *n_rest)
{
  int ndel, ncalc, nrest;
  long nn;
  shortstr buff;
  FILE *tmp, **menu;

  *n_sub = 0, *n_del = 0, *n_calc = 0, *n_rest = 0;

  if (type == 1)
    {
      menu = &menup;
      strcpy (buff, MENUP_NAME);
    }
  else
    {
      menu = &menuq;
      strcpy (buff, MENUQ_NAME);
    }

  tmp = *menu;
  *menu = fopen (buff, "rb");
  if (*menu)
    {
      while (rd_menu
             (*menu, type, *n_sub + 1, buff, &ndel, &ncalc, &nrest, &nn))
        {
          (*n_sub)++;
          *n_del += ndel;
          *n_calc += ncalc;
          *n_rest += nrest;
        }
      fclose (*menu);
    }
  *menu = tmp;
  return;
}


static void
readtext (char *fname)
{
  FILE *txt;
  trim (fname);
  txt = fopen (fname, "r");
  if (txt == NULL)
    {
      warnanykey (10, 10, " File not found");
      return;
    }
  showtext (1, 1, 80, 1, "", txt);
  fclose (txt);
}


void
editModel (int edit)
{
  int n = 1, i, j;
  void *pscr = NULL;
  int edited = FALSE;
  char tabMenu[STRSIZ], tabName[80];
  char menuName[30];
  char tabhelp[4][10] = { "s_mdl_1", "s_mdl_2", "s_mdl_3", "s_mdl_4" };

cont:
  do
    {
      if (edit)
        strcpy (menuName, "Edit Model");
      else
        strcpy (menuName, "View Model");
      strcpy (tabMenu, "\017");
      for (i = 0; i < 5; i++)
        {
          strcpy (tabName, modelTab[i].headln);
          trim (tabName);
          sbld (tabName, " %s", tabName);
          for (j = strlen (tabName); j < 15; j++)
            tabName[j] = ' ';
          tabName[15] = 0;
          strcat (tabMenu, tabName);
        }
      if (edit)
        menu1 (17, 15, "", tabMenu, "s_mdl_e", &pscr, &n);
      else
        menu1 (17, 15, "", tabMenu, "s_mdl_v", &pscr, &n);
      if (n > 0 && n <= 5)
        edited = edittable (1, 1, &modelTab[n - 1], 1, tabhelp[n - 1], !edit)
          || edited;
    }
  while (n != 0);

  if (edited)
    {
      int num = getModelNumberSymb ();
      if (mess_y_n (15, 19, " Save corrections ?"))
        {
          if (loadModel (TRUE))
            {
              writeModelFiles (num, "models");
            }
          else
            {
              goto cont;
            }
        }
      else
        {
          readModelFiles (num, "models");
        }
    }
}


void
menuhelp (void)
{
  scrcolor (Red, BGmain);
  goto_xy (23, 4);
  print ("Abstract");
  scrcolor (FGmain, BGmain);
  goto_xy (1, 6);
  print ("    CompHEP package is created for calculation   \n");
  print (" of decay and high energy collision processes of \n");
  print (" elementary particles in the tree approximation. \n");
  print ("    The main idea put into the CompHEP was to    \n");
  print (" make available passing from the Lagrangian to   \n");
  print (" the final distributions effectively, with the   \n");
  print (" high level of automatization.                   \n");
  print ("    Use the F2 key to get the information about  \n");
  print (" interface facilities and the F1 key to get      \n");
  print (" online help.                                    \n");

  scrcolor (Black, BGmain);
  chepbox (1, 5, 50, 16);
  scrcolor (FGmain, BGmain);
}


void
modelinfo (void)
{
  int num = getModelNumberSymb ();
  goto_xy (5, 1);
  scrcolor (Red, BGmain);
  print ("   Model:  ");
  scrcolor (FGmain, BGmain);
  print ("%s", copy (modelmenu, num * 22 - 20, 22));
}


void
processinfo (void)
{
  goto_xy (5, 3);
  scrcolor (Red, BGmain);
  print (" Process:  ");
  scrcolor (FGmain, BGmain);
  print ("%s", getProcessch ());
}


void
diagramsinfo (void)
{
  int n_sub, n_del, n_calc, n_rest;

  diag_stat (1, &n_sub, &n_del, &n_calc, &n_rest);
  if (!n_sub)
    return;
  goto_xy (15, 5);
  scrcolor (Red, BGmain);
  print (" Feynman diagrams \n");
  scrcolor (FGmain, BGmain);
  print ("      diagrams in      subprocesses are constructed.\n");
  print ("      diagrams are deleted.");
  scrcolor (Blue, BGmain);
  goto_xy (1, 6);
  print ("%d", n_del + n_calc + n_rest);
  goto_xy (20, 6);
  print ("%d", n_sub);
  goto_xy (1, 7);
  print ("%u", n_del);
  while (where_y () < 7)
    print (" ");

}


void
sq_diagramsinfo (void)
{
  int n_sub, n_del, n_calc, n_rest;


  diag_stat (2, &n_sub, &n_del, &n_calc, &n_rest);
  if (!n_sub)
    return;

  goto_xy (15, 9);
  scrcolor (Red, BGmain);
  print (" Squared diagrams \n");
  scrcolor (FGmain, BGmain);
  print ("      diagrams in      subprocesses are constructed.\n");
  print ("      diagrams are deleted.\n");
  print ("      diagrams are calculated.");
  scrcolor (Blue, BGmain);
  goto_xy (1, 10);
  print ("%d", n_del + n_calc + n_rest);
  goto_xy (20, 10);
  print ("%d", n_sub);
  goto_xy (1, 11);
  print ("%d", n_del);
  while (where_x () < 7)
    print (" ");
  goto_xy (1, 12);
  print ("%d", n_calc);
  goto_xy (1, 13);

}

static void
f7_key_prog (int x)
{
  static char delstr[5] = "}D{";
  inkeyString = delstr;
}

static void
f8_key_prog (int x)
{
  static char delstr[5] = "}R{";
  inkeyString = delstr;
}


static void
menu_f (int col, int row, char *label, char *f_name, char *help,
        void *hscr, int *kk)
{
  FILE *f;
  int nline, i;

  char *menustr;
  char ch1, ch2;

  f = fopen (f_name, "r");
  if (f == NULL)
    return;

  fread (&ch1, 1, 1, f);
  fread (&ch2, 1, 1, f);

  fseek (f, 0, SEEK_END);
  nline = ftell (f) / ch2;

  menustr = (char *) malloc (2 + nline * ch1);
  fseek (f, 2, SEEK_SET);

  menustr[0] = ch1;
  for (i = 1; i <= nline; i++)
    {
      fread (menustr + 1 + (i - 1) * ch1, ch1, 1, f);
      fseek (f, ch2 - ch1, SEEK_CUR);
    }
  fclose (f);
  menustr[1 + nline * ch1] = 0;
  menu1 (col, row, label, menustr, help, hscr, kk);

}


void
sqdiagrmenu (void)              /*  form menu of squared diagrams  */
{
  void *pscr = NULL;
  shortstr buff;

  if (subproc_sq == 1)
    {
      nsub = 1;
      return;
    }

  strcpy (buff, MENUQ_NAME);
  menu_f (9, 16, "NN      Subprocess                Del   Calc  Rest ", buff,
          "s_sq_proc", &pscr, &nsub);
  if (nsub)
    {
      put_text (&pscr);
    }
}


void
viewsqdiagr (void)              /*  View squared diagrams  */
{
  nsub = 1;
  do
    {
      sqdiagrmenu ();
      if (nsub != 0)
        showgraphs (2);
      sq_diagramsinfo ();
    }
  while (!(nsub == 0 || subproc_sq == 1));      /*  Esc  */
  sq_diagramsinfo ();
}


void
viewfeyndiag (int del_mode)
{
  void *pscr = NULL;
  int upr = del_mode ? 1 : -1;

  if (del_mode)
    {
      f3_key[4] = f7_key_prog;
      f3_key[5] = f8_key_prog;
    }
  nsub = 1;
  do
    {
      if (subproc_f == 1)
        nsub = 1;
      else
        {
          menu_f (9, 11, "NN        Subprocess              Del   Rest ",
                  MENUP_NAME, "s_proc", &pscr, &nsub);
        }
      if (nsub > 0)
        {
          showgraphs (upr);
          if (del_mode)
            diagramsinfo ();
        }
    }
  while (!(nsub == 0 || subproc_f == 1));

  if (del_mode)
    {
      f3_key[4] = NULL;
      f3_key[5] = NULL;
    }
}

int
viewresults (int toDelete)
{
  int i, k = 1;
  int doserror;
  void *pscr = NULL;
  void *pscr2 = NULL;
  void *pscr3 = NULL;

  shortstr newname;
  shortstr oldname;
  shortstr f_name;
  midstr menustr;
  searchrec s;

  menustr[0] = 15;
  doserror = find_first (scat ("%s/%s/*.*", pathtouser, pathtoresults), &s, anyfile);
  while (0 == doserror && k <= 2048)
    {
      if (strcmp (s.name, "Makefile") && strcmp (s.name, "n_comphep")
          && strcmp (s.name, "diag_view"))
        {
          for (i = 0; (i < strlen (s.name)) && (i < 15); ++i)
            {
              menustr[k++] = s.name[i];
            }
          i = strlen (s.name);
          for (; i < 15; ++i)
            {
              menustr[k++] = ' ';
            }
        }
      doserror = find_next (&s);
    }
  menustr[k] = 0;

  if (menustr[1] == 0)
    {
      if (!toDelete)
        messanykey (10, 15, "Results is empty");
      return 1;
    }

  if (toDelete)
    messanykey (10, 10,
                "There are files in Results.\n Clear or rename this directory");

  while (1)
    {
      int kmenu = 1;

      menu1 (10, 10, "", "\010"
             " View   "
	     " Delete "
	     " Rename ",
	     "s_res", &pscr, &kmenu);

      switch (kmenu)
        {
        case 0:
          return 0;
        case 1:
          k = 1;
          do
            {
              menu1 (10, 10, "", menustr, "", &pscr2, &k);
              if (k > 0)
                {
                  sprintf (f_name, "%s/%s/%.15s", pathtouser, pathtoresults, menustr + k * 15 - 14);
                  readtext (f_name);
                }
            }
          while (k != 0);
          break;

        case 2:
          if (mess_y_n (10, 13, " Delete files "))
            {
              system (scat ("bash -c 'make -C %s clean'", pathtoresults));
              return 1;
            }

        case 3:
          strcpy (newname, " ");
          while (1)
            {
              get_text (1, maxRow (), maxCol (), maxRow (), &pscr3);
              goto_xy (1, maxRow ());
              print ("Enter new name: ");
              k = str_redact (newname, 1, 30);
              if (k == KB_ESC)
                {
                  goto_xy (1, 24);
                  clr_eol ();
                  break;
                }
              if (k == KB_ENTER)
                {
                  trim (newname);
                  sbld (oldname, "%s/%s", pathtouser, pathtoresults);
                  sbld (newname, "%s/%s", pathtouser, newname);
                  if (!rename (oldname, newname))
                    {
                      chepmkdir (oldname);
                      system (scat ("bash -c 'cp %s/Makefile %s/Makefile'",   newname, oldname));
                      system (scat ("bash -c 'cp %s/n_comphep %s/n_comphep'", newname, oldname));
                      system (scat ("bash -c 'cp %s/diag_view %s/diag_view'", newname, oldname));
                      put_text (&pscr);
                      put_text (&pscr3);
                      return 1;
                    }
                  else {
                    warnanykey (10, 15, " Can't rename the directory.\nTry another name");
		  }
                }
              put_text (&pscr3);
            }
        }
    }
}


void
f3_key_prog (int x)
{
  int i;
  for (i = 0; i < 8 && f3_key_prog != f3_key[i]; i++);
  if (i < 8)
    f3_key[i] = NULL;           /* LOCK */

  editModel (FALSE);

  if (i < 8)
    f3_key[i] = f3_key_prog;    /* UNLOCK */

}



void
f4_key_prog (int x)
{
  int nsubtmp;
  int i;
  for (i = 0; i < 8 && f4_key_prog != f3_key[i]; i++);
  if (i < 8)
    f3_key[i] = NULL;           /* LOCK */
  nsubtmp = nsub;
  viewfeyndiag (FALSE);
  nsub = nsubtmp;
  if (i < 8)
    f3_key[i] = f4_key_prog;    /* UNLOCK */
}

void
f5_key_prog (int x)
{
  int kmenu = 1;
  void *pscr = NULL;

  int nfun;
  for (nfun = 0; nfun < 8 && f5_key_prog != f3_key[nfun]; nfun++);
  if (nfun < 8)
    f3_key[nfun] = NULL;



  while (kmenu)
    {
      char strmen[] = "\040"
        " Symb. momentum conservation OF1" " Number of QCD colors =   Nc    ";

      if (consLow)
        improveStr (strmen, "OF1", "ON ");
      else
        improveStr (strmen, "OF1", "OFF");

      if (getNcinflimit ())
        improveStr (strmen, "Nc", "Inf");
      else
        improveStr (strmen, "Nc", "3");


      menu1 (20, 18, "Switches", strmen, "s_res", &pscr, &kmenu);
      switch (kmenu)
        {
        case 1:
          consLow = !consLow;
          break;
        case 2:
          {
            int num = getNcinflimit ();
            setNcinflimit (!num);
            break;
          }
        }

    }
  if (nfun < 8)
    f3_key[nfun] = f5_key_prog;

}


void
f6_key_prog (int x)
{
  int i;
  for (i = 0; i < 8 && f6_key_prog != f3_key[i]; i++);
  if (i < 8)
    f3_key[i] = NULL;
  viewresults (0);
  if (i < 8)
    f3_key[i] = f6_key_prog;
}

void
f9_key_prog (int x)
{
  if (mess_y_n (56, maxRow () - 5, " Quit session?  "))
    {
      saveent (menulevel);
      finish ("End of CompHEP symbolical session.");
      exit (0);
    }
}

int oneclickcode_dir (void) 
{
  int i, k = 1;
  int doserror;
  void *pscr = NULL;
  void *pscr2 = NULL;
  void *pscr3 = NULL;

  shortstr newname;
  shortstr oldname;
  shortstr f_name;
  midstr menustr;
  searchrec s;

  menustr[0] = 15;
  doserror = find_first (scat ("%s/oneclickcode/*.*", pathtouser), &s, anyfile);
  while (0 == doserror && k <= 2048)
    {
      if (strcmp (s.name, "Makefile") && strcmp (s.name, "n_comphep")
          && strcmp (s.name, "diag_view"))
        {
          for (i = 0; (i < strlen (s.name)) && (i < 15); ++i)
            {
              menustr[k++] = s.name[i];
            }
          i = strlen (s.name);
          for (; i < 15; ++i)
            {
              menustr[k++] = ' ';
            }
        }
      doserror = find_next (&s);
    }
  menustr[k] = 0;

  messanykey (10, 10, "The oneclickcode exists.\n Clear or rename this directory");
  while (1) {
    int kmenu = 1;
    menu1 (10, 10, "", "\010"
                       " View   "
                       " Delete "
                       " Rename ",
                       "s_res", &pscr, &kmenu);

    switch (kmenu) {
      case 0:
        return 0;
      case 1:
        k = 1;
        do {
          menu1 (10, 10, "", menustr, "", &pscr2, &k);
          if (k > 0) {
   	    sprintf (f_name, "%s/oneclickcode/%.15s", pathtouser, menustr + k * 15 - 14);
   	    readtext (f_name);
          }
        } while (k != 0);
        break;

      case 2:
        if (mess_y_n (10, 13, " Delete files ")) {
          system ("rm -rf oneclickcode");
          return 1;
        }
        break;

      case 3:
        strcpy (newname, " ");
        while (1) {
          get_text (1, maxRow (), maxCol (), maxRow (), &pscr3);
          goto_xy (1, maxRow ());
          print ("Enter new name: ");
          k = str_redact (newname, 1, 30);
          if (k == KB_ESC) {
            goto_xy (1, 24);
            clr_eol ();
            return 0;
          }
          if (k == KB_ENTER) {
            trim (newname);
            sbld (oldname, "%s/oneclickcode", pathtouser);
            sbld (newname, "%s/%s", pathtouser, newname);
            if (!rename (oldname, newname)) {
              put_text (&pscr);
              put_text (&pscr3);
              return 1;
            } else {
              warnanykey (10, 15, " Can't rename the directory.\nTry another name");
            }
          }
          put_text (&pscr3);
        }
        break;
    }
  }

  return 0;
}
