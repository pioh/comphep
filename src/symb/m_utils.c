/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"

#include "physics.h"
#include "screen.h"
#include "read_mdl.h"
#include "m_utils.h"

#define nhardmdl 5
#define newmodeltxt " CREATE NEW MODEL     "

#define MAX_MODEL_NUM  ((STRSIZ-2)/22 -1)

/*
1. Otkrivaet vse files var*.mdl v dir models i chitaet nazvanie modeli iz pervoy stroki.
   Nazvaniya sbrasivautsya v string modelmenu. K kontsu modelmenu prisobachivaetsya
   string newmodeltxt.
   
2. Pitaetsya otkrit' files cpart (kompozitnie chastitsi). Esli ne poluchaetsya,
   to sozdaet etot file v modeli pustim.

3. Pitaetsya otkrit' file beams.mdl (info o vozmozhnih beams). Esli ne poluchaetsya,
   to sozdaet pustoy file beams.mdl.

4. Pitaetsya otkrit' file strfuns.mdl (info o vozmozhnih str. function). Esli ne poluchaetsya,
   to sozdaet pustoy file strfuns.mdl.
*/
void 
fillModelMenu (void)
{
  int i, j;
  FILE *txt;
  char name[80];
  char filepath[STRSIZ];
  strcpy (modelmenu, "\026");
  maxmodel = 0;
  for (i = 1; i <= MAX_MODEL_NUM; i++)
    {
      txt = fopen (scat ("%smodels%c%s%d.mdl", pathtouser, f_slash, mdFls[0], i), "r");
      if (txt == NULL)
	goto exi;
      fgets (name, 60, txt);
      trim (name);
      for (j = strlen (name); j < 21; j++)
	name[j] = ' ';
      name[21] = 0;
      strcat (modelmenu, " ");
      strcat (modelmenu, name);
      fclose (txt);
      maxmodel++;
      strcpy (filepath, scat ("%smodels%c%s%d.mdl", pathtouser, f_slash, mdFls[4], i));
      txt = fopen (filepath, "r");
      if (txt == NULL)
	{
	  txt = fopen (filepath, "w");
	  fprintf (txt, "%s", name);
	  fprintf (txt, "\n Composite \n");
	  fprintf (txt, " Abr  |> elementary particles               <|> Comment                        <|\n");
	  fclose (txt);
	}
      else
	fclose (txt);
    }
exi:
  strcpy (filepath, scat ("%smodels%cbeams.mdl", pathtouser, f_slash));
  txt = fopen (filepath, "r");
  if (txt == NULL)
    {
      txt = fopen (filepath, "w");
      fprintf (txt, "\n Beams \n");
      fprintf (txt, " Name     |> Mass    <|> Content                <|> Structure Functions      <|\n");
      fclose (txt);
    }
  else
    {
      fclose (txt);
    }
  strcpy (filepath, scat ("%smodels%cstrfuns.mdl", pathtouser, f_slash));
  txt = fopen (filepath, "r");
  if (txt == NULL)
    {
      txt = fopen (filepath, "w");
      fprintf (txt, "\n Strfuns \n");
      fprintf (txt, " Name                   <|> Number    <|\n");
      fclose (txt);
    }
  else
    {
      fclose (txt);
    }
  if (maxmodel < MAX_MODEL_NUM)
    strcat (modelmenu, newmodeltxt);
}


int 
deletemodel (int n)
{
  int i;
  char from[STRSIZ], to[STRSIZ];
  searchrec s;


  sprintf (from, "%smodels%c%s%d.mdl", pathtocomphep, f_slash, mdFls[0], n);

  if (find_first (from, &s, archive) == 0)
    {
      if (mess_y_n (44, 10, "This model can not be deleted.\nRestore original version?"))
	for (i = 0; i < 4; i++)
	  {
	    sprintf (from, "%smodels%c%s%d.mdl", pathtocomphep, f_slash, mdFls[i], n);
	    sprintf (to, "%smodels%c%s%d.mdl", pathtouser, f_slash, mdFls[i], n);
	    copyfile (from, to);
	  }
      else
	return FALSE;
    }
  else if (mess_y_n (56, 10, "Delete model?"))
    {
      for (i = 0; i < 4; i++)
	{
	  sprintf (to, "%smodels%c%s%d.mdl", pathtouser, f_slash, mdFls[i], n);
	  unlink (to);
	}

      if (n < maxmodel)
	{
	  for (i = 0; i < 4; i++)
	    {
	      sprintf (from, "%smodels%c%s%d.mdl", pathtouser, f_slash, mdFls[i], maxmodel);
	      sprintf (to, "%smodels%c%s%d.mdl", pathtouser, f_slash, mdFls[i], n);
	      rename (from, to);
	    }
	}
      return TRUE;
    }
  return FALSE;
}



int 
makenewmodel (void)
{
  int i;
  int num;
  int key;
  int nmdl;
  int npos = 1;
  shortstr newName;
  char oldModels[STRSIZ];
  void  *pscr = NULL;

  strcpy (newName, " ");
  do {
    goto_xy (1, 7);
    print ("Model name :");
    key = str_redact (newName, npos, 25);
  } while (key != KB_ENTER && key != KB_ESC);
  trim (newName);
  if (key == KB_ESC || strlen (newName) == 0) {
    return FALSE;
  }
  sbld (newName, "_%s", newName);
  nmdl = 1;
  strcpy (oldModels, copy (modelmenu, 1, (int) strlen (modelmenu) - 22));
  menu1 (5, 10, "Choose a template", oldModels, NULL, &pscr, &nmdl);
  if (nmdl == 0)
    return FALSE;
  put_text (&pscr);
  readModelFiles (nmdl, "models");
  for (i = 0; i < 4; i++) {
    strcpy (modelTab[i].mdlName, newName);
  }
  num = getModelNumberSymb();
  writeModelFiles (num, "models");
  fillModelMenu ();
  clrbox (5, 10, 30, 12);
  return TRUE;
}


int 
continuetest (void)
{
  shortstr txt;
  int k, ndel, ncalc, nrest;
  long recpos;

  menuq = fopen (MENUQ_NAME, "rb");

  for (k = 1; k <= subproc_sq; k++)
    {
      rd_menu (menuq, 2, k, txt, &ndel, &ncalc, &nrest, &recpos);
      if (nrest != 0)
	{
	  fclose (menuq);
	  return FALSE;
	}
    }
  fclose (menuq);
  return TRUE;
}


void 
clear_tmp (void)
{
  unlink (CATALOG_NAME);
  unlink (ARCHIV_NAME);
}


int
writeHadrons(char * path)
{
  longstr fname;

  sprintf (fname, "%s%s%cbeams.mdl", pathtouser, path, f_slash);
  writetable (&hadron_tab, fname);
  return 1;
}


int
writeStrFuns(char * path)
{
  longstr fname;
  
  sprintf (fname, "%s%s%cstrfuns.mdl", pathtouser, path, f_slash);
  writetable (&strfun_tab, fname);
  return 1;
}


int 
writeModelFiles (int l, char * path)
{
  int i;
  longstr fName;

  sprintf (fName, "%s%s%c%%s%d.mdl", pathtouser, path, f_slash, l);

  for (i = 0; i < 5; i++)
    {
      longstr sname;
      FILE * ftest;
      sprintf (sname, fName, mdFls[i]);
      ftest = fopen(sname, "a");
      if (ftest == NULL) {
        longstr mess;
        sprintf (mess, "I can't store the file %s", sname);
        messanykey (10, 15, mess);
        return 0;
      }
      else fclose(ftest);

      writetable (&modelTab[i], scat (fName, mdFls[i]));
    }
  return 1;
}

int 
editBeams (void)
{
  int edit = 1;
  int edited = 0;
  int exit = 1;
  char tabhelp[10] =  {"s_beam_1"};

  do
    {
      edited = edittable (1, 1, &hadron_tab, 1, tabhelp, !edit) || edited;
    
      if (readhadrons (1))
      {
        writeHadrons("models");
        exit = 0;
      }
    }
  while (exit);
  return 1;
}


int 
editStrFuns (void)
{
  int edit = 1;
  int edited = 0;
  int exit = 1;
  char tabhelp[10] =  {"s_strfn_1"};

  do
    {
      edited = edittable (1, 1, &strfun_tab, 1, tabhelp, !edit) || edited;
    
      if (readstrfuns (1))
      {
        writeStrFuns("models");
        exit = 0;
      }
    }
  while (exit);
  return 1;
}
