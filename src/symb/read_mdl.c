/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov
* --------------------------------------------------
*/
#include <limits.h>
#include <math.h>
#include <ctype.h>

#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/read_func.h"
#include "service2/include/parser.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "chep_crt/include/chep_crt.h"
#include "chep_crt/include/file_scr.h"

#include "physics.h"
#include "process.h"
#include "model.h"
#include "sos.h"
#include "pre_read.h"
#include "pvars.h"
#include "read_mdl.h"

#define minvarmem  ((unsigned)sizeof(struct varrec)  + 1 - STRSIZ)
#define minlagrmem ((unsigned)sizeof(struct algvert) + 1 - STRSIZ)

table modelTab[5] =
{
  {"", "", "", NULL},
  {"", "", "", NULL},
  {"", "", "", NULL},
  {"", "", "", NULL},
  {"", "", "", NULL}
};

table hadron_tab =  {"", "", "", NULL};
table strfun_tab =  {"", "", "", NULL};

static int lastModel = 0;
static vshortstr Sqrt2FuncTxt = "sqrt(2)";
static midstr tabName;
static int nLine;


static void 
errorMessage (char *fieldName, char *errComment)
{
  char *mess;
  switch (rderrcode)
    {
    case braketexpected:
      mess = "Bracket expected";
      break;
    case unexpectedcharacter:
      mess = "Unexpected character";
      break;
    case operationexpected:
      mess = "Operation expected";
      break;
    case toolongidentifier:
      mess = "Too long identifier/number";
      break;
    case toomanyagruments:
      mess = "Too many arguments";
      break;
    case cannotread:
      mess = "Undefined identifier";
      break;
    case cannotevaluate:
      mess = "Can not evaluate function";
      break;
    case unknownidentifier:
      mess = "Unknown identifier";
      break;
    case unexpectedoperation:
      mess = "Unexpected operation";
      break;
    case unknownfunction:
      mess = "Unknown function";
      break;
    case wrongnumberofagr:
      mess = "Wrong number of arguments";
      break;
    case typemismatch:
      mess = "Type mismatch";
      break;
    case naninoperation:
      mess = "NAN in operation";
      break;
    case toolargenumber:
      mess = "Too large number";
      break;
    case indexuncompatibility:
      mess = "Index uncompatibility";
      break;
    case rangecheckerror:
      mess = "Range check error";
      break;
    default:
      mess = "Unknown error code";
      break;
    }

  if (strcmp (errComment, "*") == 0)
    sprintf (errorText, "Error in table '%s' line %d field '%s'\nposition %u: %s",
	     tabName, nLine, fieldName, rderrpos, mess);
  else
    sprintf (errorText, "Error in table '%s' line %d field '%s' \n %s",
	     tabName, nLine, fieldName, errComment);
  warnanykey (2, 10, errorText);
#undef MAXERRMSG
}


static int 
tabCharPos (char *str, int n)
{
  int k = 0;
  int nn = 0;
  if (n == 0)
    return 0;
  while (str[k] != 0)
    {
      if (str[k] == '|')
	{
	  nn++;
	  if (nn == n)
	    return (k + 1);
	}
      k++;
    }
  return k;
}


static int 
isVarName (char *s)
{
  int i;
  if (!isalpha (s[0]))
    return 0;
  for (i = 1; s[i] != 0; i++)
    if (!isalnum (s[i]))
      return 0;

  if (s[0] == 'p' || s[0] == 'm' || s[0] == 'M')
    {
      i = 1;
      while (isdigit (s[i]))
	i++;
      if (s[i] == 0)
	return 0;
    }

  if (!strcmp (s, "G5"))
    return 0;
  if (!strcmp (s, "g5"))
    return 0;
  if (!strcmp (s, "I"))
    return 0;
  if (!strcmp (s, "i"))
    return 0; 
  return 1;
}

static int 
isOriginName (char *s)
{
  char newName[10], oldName[10];
  int i = -1, k;

  strcpy (newName, s);
  for (i = 0; newName[i] != 0; i++)
    newName[i] = toupper (newName[i]);

  for (k = 1; k == nmodelvar; k++)
    {
      strcpy (oldName, modelvars[k].varname);
      for (i = 0; oldName[i] != 0; i++)
	oldName[i] = toupper (oldName[i]);
      if (!strcmp (newName, oldName))
	return 0;
    }
  return 1;
}


static int 
rd_3 (char *s, double *p)
{
  int i;

  if (isdigit (*s))
    {
      if (strchr (s, '.'))
	return 0;
      sscanf (s, "%lf", p);
      return 1;
    }

  if (strcmp (s, strongconst) == 0 || strcmp (s, "i") == 0)
    return 0;

  for (i = 1; i <= nmodelvar; i++)
    {
      if (strcmp (modelvars[i].varname, s) == 0)
	{
	  *p = modelvars[i].varvalue;
	  if (modelvars[i].able)
	    return 1;
	  else
	    return -1;
	}
    }
  return 0;
}


static int 
readvars (int check)
{
  char numtxt[60];
  char name[60];
  char *ss;
  char *endstr;
  int funcShift;
  varlist mvars;
  double varvalue_tmp;
  linelist ln;
  int ggOn = 0;
  int nv = 3;			/* Plus 0,i,Sqrt2 */

  ln = vars_tab.strings;
  while (ln != NULL)
    {
      ln = ln->next;
      nv++;
    }
  ln = func_tab.strings;
  while (ln != NULL)
    {
      ln = ln->next;
      nv++;
    }

  strcpy(tabName, vars_tab.headln);
  if (modelvars)
    free (modelvars);
  modelvars = m_alloc (nv * sizeof (*modelvars));
  nmodelvar = 0;

  strcpy (modelvars[0].varname, "0");
  modelvars[0].varvalue = 0.;
  modelvars[0].func = NULL;
  modelvars[0].able = 1;

  nLine = 1;
  ln = vars_tab.strings;
  while (ln)
    {
      ss = ln->line;

      sscanf (ss, "%[^|]%*c%[^|]", name, numtxt);
      trim (name);
      trim (numtxt);

      if (check && (!isVarName (name)))
	{
	  errorMessage ("Name", scat ("incorrect name '%s'", name));
	  goto errExi1;
	}

      if (check && (!isOriginName (name)))
	{
	  errorMessage ("Name", scat ("duplicate name '%s'", name));
	  goto errExi1;
	}

      varvalue_tmp = strtod (trim (numtxt), &endstr);

      if (check && (endstr != numtxt + strlen (numtxt)))
	{
	  errorMessage ("Value", scat ("wrong number '%s'", numtxt));
	  goto errExi1;
	}

      if (strcmp (name, strongconst) == 0)
	{
	  if (ggOn)
	    {
	      errorMessage ("Name", scat ("duplicate name '%s'", name));
	      goto errExi1;
	    }
	  ggOn = 1;
	  mvars = modelvars + nv - 1;
	}
      else
	{
	  nmodelvar++;
	  mvars = modelvars + nmodelvar;
	  if (nmodelvar > USHRT_MAX - 1)
	    {
	      errorMessage ("Name", "too many parameters");
	      goto errExi1;
	    }
	}

      strcpy (mvars->varname, name);
      mvars->varvalue = varvalue_tmp;
      mvars->func = NULL;
      mvars->able = 1;

      ln = ln->next;
      nLine++;
    }

  nmodelvar++;
  mvars = modelvars + nmodelvar;
  strcpy (mvars->varname, "i");
  mvars->func = NULL;
  mvars->varvalue = 0.;
  mvars->able = 1;
  mvars++;
  nmodelvar++;
  strcpy (mvars->varname, "Sqrt2");
  mvars->func = Sqrt2FuncTxt;
  mvars->varvalue = sqrt (2.0);
  mvars->able = 1;
  nLine = 1;
  ln = func_tab.strings;
  strcpy(tabName, func_tab.headln);
  funcShift = tabCharPos (func_tab.format, 1);

  while (ln)
    {
      ss = ln->line;
      sscanf (ss, "%[^|]", name);
      trim (name);

      if (!isVarName (name))
	{
	  errorMessage ("Name", scat ("incorrect name '%s'", name));
	  goto errExi1;
	}

      if (check && (!isOriginName (name)))
	{
	  errorMessage ("Name", scat ("duplicate name '%s'", name));
	  goto errExi1;
	}

      calcExpression (ln->line + funcShift, rd_3, &varvalue_tmp);
      if (check && rderrcode && (rderrcode != unknownfunction))
	{
	  errorMessage ("Expression", "*");
	  goto errExi1;
	}

      nmodelvar++;
      if (nmodelvar > USHRT_MAX - 1)
	{
	  errorMessage ("Name", "too many parameters");
	  goto errExi1;
	}

      mvars = modelvars + nmodelvar;
      if (rderrcode == unknownfunction)
	{
	  mvars->able = 0;
	  varvalue_tmp = 0;
	}
      else
	mvars->able = 1;
      mvars->varvalue = varvalue_tmp;
      strcpy (mvars->varname, name);
      mvars->func = ln->line + funcShift;
      ln = ln->next;
      nLine++;
    }
  nmodelvar = nv - 1;
  return 1;

errExi1:
  free (modelvars);
  modelvars = NULL;
  return 0;
}


static void 
findvar (char *txt, double *num, int *err)
{
  int i;
  trim (txt);
  for (i = 1; i <= nmodelvar; i++)
    {
      if (strcmp (txt, modelvars[i].varname) == 0)
	{
	  *num = modelvars[i].varvalue;
	  *err = 0;
	  return;
	}
    }
  *err = -1;
}


static int 
ghostaddition (void)
{
  int i, nPrim;
  nPrim = nparticles;
  if ((prtclbase[nPrim - 1].spin == 2) && (prtclbase[nPrim - 1].cdim != 1))
    {
      nPrim++;
      prtclbase[nPrim - 1] = prtclbase[nparticles - 1];
      prtclbase[nparticles - 1].hlp = 't';
      prtclbase[nparticles - 1].spin = 4;
      nparticles++;
    }

  if (gaugep (nPrim))
    {
      nparticles++;
      prtclbase[nparticles - 1] = prtclbase[nPrim - 1];
      prtclbase[nparticles - 1].hlp = 'c';
      nparticles++;
      prtclbase[nparticles - 1] = prtclbase[nPrim - 1];
      prtclbase[nparticles - 1].hlp = 'C';

      if (strcmp (prtclbase[nPrim - 1].massidnt, "0") != 0)
	{
	  nparticles++;
	  prtclbase[nparticles - 1] = prtclbase[nPrim - 1];
	  prtclbase[nparticles - 1].hlp = 'f';
	}
      for (i = nPrim + 1; i <= nparticles; i++)
	prtclbase[i - 1].spin = 0;
    }
  return nPrim;
}


static void 
clearlgrgn (void)
{
  algvertptr l1;
  while (lgrgn != NULL)
    {
      l1 = lgrgn;
      lgrgn = lgrgn->next;
      free (l1);
    }
}


static void 
cleardecaylist (void)
{
  decaylink v1, v2;
  int j;

  if (!prtclbase)
    return;
  for (j = 0; j < nparticles; j++)
    {
      v1 = prtclbase[j].top;
      while (v1 != NULL)
	{
	  v2 = v1;
	  v1 = v1->next;
	  free (v2);
	}
      prtclbase[j].top = NULL;
    }
}

static void 
clearLatexNames (void)
{
  int j;
  if (!prtclbase)
    return;
  for (j = 0; j < nparticles; j++)
    {
      if (!strchr ("fcCt", prtclbase[j].hlp))
	free (prtclbase[j].latex);
    }
}


static int 
isPrtclName (char *p)
{
  if (strlen (p) > 3)
    return 0;
  if (strlen (p) > 2)
    if (p[1] != '_' && p[0] != '~')
      return 0;
  return 1;
}


static int 
readparticles (int check)
{
  char *ss;
  char fullname[60];
  char massname[60];
  char imassname[60];
  char p1[60];
  char p2[60];
  char latex[STRSIZ];
  char latex_[STRSIZ];
  int i, j;
  char s[60], c[60];
  char *endstr;
  char chlp[40];
  int itmp;
  int errcode;
  int np1, np2;
  int nparticleLimit = 128;
  linelist ln;

  ln = prtcls_tab.strings;
  strcpy(tabName, prtcls_tab.headln);

  if (prtclbase)
    {
      cleardecaylist ();
      clearLatexNames ();
      free (prtclbase);
    }

  prtclbase = (prtcl_base *) malloc (nparticleLimit * sizeof (prtcl_base));
  nparticles = 0;

  for (i = nparticles; i < nparticleLimit; i++)
    {
      prtclbase[i].top = NULL;
      prtclbase[i].latex = NULL;
    }

  nLine = 1;
  while (ln != NULL)
    {
      ss = ln->line;
      if (nparticles >= nparticleLimit - 8)
	{
	  nparticleLimit += 128;
	  prtclbase = re_alloc (prtclbase, nparticleLimit * sizeof (prtcl_base));
	  if (!prtclbase)
	    {
	      errorMessage (" P ", "too many particles");
	      return 0;
	    }
	  for (i = nparticles; i < nparticleLimit; i++)
	    {
	      prtclbase[i].top = NULL;
	      prtclbase[i].latex = NULL;
	    }
	}

      sscanf (ss, "%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]",
	  fullname, p1, p2, s, massname, imassname, c, chlp, latex, latex_);
      trim (p1);
      trim (p2);
      trim (latex);
      trim (latex_);
      {
	static char fldName[2][5] =
	{" P ", " aP"};
	char *pName[2];
	pName[0] = p1;
	pName[1] = p2;

	for (i = 0; i <= 1; i++)
	  {
	    if (check && (!isPrtclName (pName[i])))
	      {
		errorMessage (fldName[i], scat ("incorrect particle name '%s'", pName[i]));
		return 0;
	      }
	    if (check)
	      {
		j = locateinbase (pName[i]);
		if (j != 0)
		  {
		    errorMessage (fldName[i], scat ("duplicate particle name '%s'", pName[i]));
		    return 0;
		  }
	      }
	  }
      }
      nparticles++;
      strcpy (prtclbase[nparticles - 1].name, p1);

      itmp = strtol (trim (s), &endstr, 10);
      if (check)
	{
	  if (s + strlen (s) != endstr)
	    {
	      errorMessage ("2*spin", "number expected");
	      return 0;
	    }
	  if ((itmp != 0) && (itmp != 1) && (itmp != 2))
	    {
	      errorMessage ("2*spin", "value out of range");
	      return 0;
	    }
	}
      prtclbase[nparticles - 1].spin = itmp;
      trim (massname);
      if (strcmp (massname, "0") == 0)
	prtclbase[nparticles - 1].mass = 0.0;
      else
	{
	  findvar (massname, &(prtclbase[nparticles - 1].mass), &errcode);
	  if (check && (errcode != 0))
	    {
	      errorMessage ("mass", scat ("unknown variable %s", massname));
	      return 0;
	    }
	  if (prtclbase[nparticles - 1].mass < 0)
	    prtclbase[nparticles - 1].mass *= -1;
	}
      strcpy (prtclbase[nparticles - 1].massidnt, massname);

      trim (imassname);
      strcpy (prtclbase[nparticles - 1].imassidnt, imassname);
      if (check && (strcmp (imassname, "0") != 0))
	{
	  double r;
	  findvar (imassname, &r, &errcode);
	  if (errcode != 0)
	    {
	      errorMessage ("width", scat ("unknown variable %s", imassname));
	      return 0;
	    }
	}
      itmp = strtol (trim (c), &endstr, 10);
      if (check)
	{
	  if (c + strlen (c) != endstr)
	    {
	      errorMessage ("color", "number expected");
	      return 0;
	    }
	  if (((itmp != 1) && (itmp != 3) && (itmp != 8)) || ((itmp == 3) && (strcmp (p1, p2) == 0)))
	    {
	      errorMessage ("color", "value out of range");
	      return 0;
	    }
	}
      prtclbase[nparticles - 1].cdim = itmp;
      trim (chlp);
      if (strcmp (chlp, "") == 0)
	strcpy (chlp, " ");
      prtclbase[nparticles - 1].hlp = toupper (chlp[0]);
      if (check)
	{
	  int ner;
	  ner = 1;
	  switch (prtclbase[nparticles - 1].hlp)
	    {
	    case ' ':
	      if (prtclbase[nparticles - 1].spin == 2 && !strcmp (prtclbase[nparticles - 1].massidnt, "0"))
		{
		  errorMessage ("aux", "Massless vector boson must\n be a gauge particle");
		  return 0;
		}
	      break;
	    case 'L':
	    case 'R':
	      if ((prtclbase[nparticles - 1].spin != 1)
		  || ((prtclbase[nparticles - 1].massidnt[0]) != '0')
		  || (strcmp (p1, p2) == 0))
		ner = 0;
	      break;
	    case '*':
	      if (prtclbase[nparticles - 1].massidnt[0] == '0')
		ner = 0;
	      break;
	    case 'G':
	      if (prtclbase[nparticles - 1].spin != 2)
		ner = 0;
	      break;
	    default:
	      ner = 0;
	    }
	  if (!ner)
	    {
	      if (prtclbase[nparticles - 1].hlp == ' ')
		errorMessage ("aux", "unexpeted character");
	      return 0;
	    }
	}
      prtclbase[nparticles - 1].latex = malloc (1 + strlen (latex));
      strcpy (prtclbase[nparticles - 1].latex, latex);

      np1 = ghostaddition ();
      if (strcmp (p1, p2) == 0)
	prtclbase[np1 - 1].anti = np1;
      else
	{
	  ++(nparticles);
	  prtclbase[nparticles - 1] = prtclbase[np1 - 1];
	  strcpy (prtclbase[nparticles - 1].name, p2);
	  prtclbase[nparticles - 1].latex = malloc (1 + strlen (latex_));
	  strcpy (prtclbase[nparticles - 1].latex, latex_);
	  if (prtclbase[np1 - 1].cdim == 3)
	    prtclbase[nparticles - 1].cdim = -3;
	  np2 = ghostaddition ();
	  prtclbase[np1 - 1].anti = np2;
	  prtclbase[np2 - 1].anti = np1;
	}
      ln = ln->next;
      nLine++;
    }

  for (i = 1; i <= nparticles; i++)
    {
      prtcl_base *with1 = &prtclbase[i - 1];

      with1->top = NULL;
      if (strchr ("fcCt", with1->hlp) != NULL)
	{
	  sbld (with1->name, "%s.%c", with1->name, with1->hlp);
	  j = prtclbase[ghostmother (i) - 1].anti;
	  switch (with1->hlp)
	    {
	    case 'c':
	      with1->anti = prtclbase[i - 1 - 1].anti + 2;
	      break;
	    case 'C':
	      with1->anti = prtclbase[i - 1 - 2].anti + 1;
	      break;
	    case 'f':
	      with1->anti = prtclbase[i - 1 - 3].anti + 3;
	      break;
	    case 't':
	      with1->anti = prtclbase[i - 1 + 1].anti - 1;
	    }
	}
    }
  return 1;
}


static int 
testLgrgn (algvertptr lgrgn)
{
  preres m;
  int n;

  m = (preres) readExpression (lgrgn->comcoef, rd_pre, act_preF, NULL);
  if (rderrcode)
    {
      errorMessage ("Factor", "*");
      return 0;
    }
  m->free = 1;

  if (m->tp > rationtp)
    {
      errorMessage ("Factor", "scalar expected");
      return 0;
    }

  if (m->maxp > 0)
    {
      errorMessage ("Factor", scat ("moments p%d are not permitable here", m->maxp));
      return 0;
    }

  for (n = 0; n < vardef->nvar; n++)
    {
      int err;
      double val;
      findvar (vardef->vars[n].name, &val, &err);
      if (err)
	{
	  errorMessage ("Factor", scat ("unknown variable '%s'", vardef->vars[n].name));
	  return 0;
	}
    }

  clearVars (vardef);
  m = (preres) readExpression (lgrgn->description, rd_pre, act_pre, NULL);
  if (rderrcode != 0)
    {
      errorMessage ("Lorentz part", "*");
      return 0;
    }
  m->free = 1;

  if (m->tp == rationtp)
    {
      errorMessage ("Lorentz part", "division is not permited here");
      return 0;
    }

  if ((m->tp == spintp) && ((prtclbase[lgrgn->fields[0] - 1].spin != 1) && (prtclbase[lgrgn->fields[1] - 1].spin != 1)
			    && (prtclbase[lgrgn->fields[2] - 1].spin != 1)))
    {
      errorMessage ("Lorentz part", "Dirac's gamma matrix not expected");
      return 0;
    }

  if ((m->tp == tenstp) && ((prtclbase[lgrgn->fields[0] - 1].spin == 1) || (prtclbase[lgrgn->fields[1] - 1].spin == 1)
			    || (prtclbase[lgrgn->fields[2] - 1].spin == 1)))
    {
      errorMessage ("Lorentz part",
		    scat ("structure as m2.m3 is not permited here.\n", " use identity G(m2)*G(m3)+G(m3)*G(m2) = 2*m2.m3 "));
      return 0;
    }

  if ((m->maxp == 4) && (lgrgn->fields[3] == 0))
    {
      errorMessage ("Lorentz part", "p4 are not permited here");
      return 0;
    }

  for (n = 0; n < vardef->nvar; n++)
    {
      int err;
      double val;
      findvar (vardef->vars[n].name, &val, &err);
      if (err)
	{
	  errorMessage ("Lorentz part", scat ("unknown variable '%s'", vardef->vars[n].name));
	  return 0;
	}
    }

  clearVars (vardef);

  for (n = 0; n <= 3; n++)
    {
      int ind1, ind2, np;

      ind1 = 0;
      np = lgrgn->fields[n];
      if (np != 0)
	switch (prtclbase[np - 1].spin)
	  {
	  case 2:
	    ind1 = 1;
	    break;
	  case 4:
	    ind1 = 3;
	  }

      ind2 = 0;
      if (inset (n + 1, m->indlist))
	ind2 += 1;
      if (inset (n + 1 + 4, m->indlist))
	ind2 += 2;
      if (ind1 != ind2)
	{
	  errorMessage ("Lorentz part", scat ("index 'm%d'  unbalanced", n + 1));
	  return 0;
	}
    }
  return 1;
}


static int 
readlagrangian (int check)
{
  algvertptr lgrgn1, lgrgn2;
  int i, j, mm;
  char p1[60], p2[60], p3[60], p4[60];
  char *ss;
  char *pPtr[4];
  int factorShift, lorentzShift;
  arr4byte f_copy;
  int mLine, totcolor, color, spinorNumb;
  linelist ln;
  static char fName[4][5] =
  {"P1", "P2", "P3", "P4"};
  polyvars var_testing =
  {0, NULL};

  vardef = &(var_testing);

  clearlgrgn ();
  nLine = 1;
  ln = lgrng_tab.strings;
  factorShift = tabCharPos (lgrng_tab.format, 4);
  lorentzShift = tabCharPos (lgrng_tab.format, 5);
  strcpy(tabName, lgrng_tab.headln);
  while (ln != NULL)
    {
      ss = ln->line;
      sscanf (ss, "%[^|]%*c%[^|]%*c%[^|]%*c%[^|]", p1, p2, p3, p4);
      pPtr[0] = p1;
      pPtr[1] = p2;
      pPtr[2] = p3;
      pPtr[3] = p4;

      lgrgn1 = (algvertptr) m_alloc (sizeof (*lgrgn1));
      lgrgn1->next = lgrgn;
      lgrgn = lgrgn1;
      lgrgn->comcoef = ln->line + factorShift;
      lgrgn->description = ln->line + lorentzShift;
      for (i = 0; i < 4; i++)
	{
	  trim (pPtr[i]);
	  if (pPtr[i][0] != 0)
	    {
	      j = locateinbase (pPtr[i]);
	      if (check && j == 0)
		{
		  errorMessage (fName[i], scat (" unknown particle %s", pPtr[i]));
		  return 0;
		}
	      lgrgn->fields[i] = j;
	    }
	  else
	    {
	      if (i == 3)
		lgrgn->fields[i] = 0;
	      else
		{
		  errorMessage (fName[i], "particle name is expected");
		  return 0;
		}
	    }
	}
      if (check)
	{
	  totcolor = 1;
	  for (mm = 0; ((mm < 4) && (lgrgn->fields[mm] != 0)); mm++)
	    {
	      color = prtclbase[lgrgn->fields[mm] - 1].cdim;
	      if (color == -3)
		color = 5;
	      totcolor = totcolor * color;
	    }
	  if ((totcolor != 1) && (totcolor != 15) && (totcolor != 64) && (totcolor != 120) && (totcolor != 512))
	    {
	      errorMessage ("Lorentz part", "wrong color structure");
	      return 0;
	    }
	  spinorNumb = 0;
	  for (mm = 0; ((mm < 4) && (lgrgn->fields[mm] != 0)); mm++)
	    {
	      if (prtclbase[lgrgn->fields[mm] - 1].spin == 1)
		spinorNumb++;
	    }
	  if ((spinorNumb != 0) && (spinorNumb != 2))
	    {
	      errorMessage ("Lorentz part", "wrong spinor structure");
	      return 0;
	    }
	}
      if (!testLgrgn (lgrgn))
	{
	  clearVars (vardef);
	  return 0;
	}
      ln = ln->next;
      nLine++;
    }

  clearVars (vardef);
  clearpregarbage ();
  lgrgn1 = lgrgn;		/* Sorting */
  do
    {
      lgrgn1->factor = 1;
      for (i = 0; i < 4 && lgrgn1->fields[i]; i++)
	{
	  int hlp = prtclbase[lgrgn1->fields[i] - 1].hlp;
	  if (hlp == 'C')
	    break;
	  else if (hlp == 'c')
	    {
	      lgrgn1->factor = -1;
	      break;
	    }
	}
      for (i = 1; i <= 4; i++)
	lgrgn1->perm[i - 1] = i;
      i = 1;
      while (i < 4)
	if (lgrgn1->fields[i - 1] >= lgrgn1->fields[i + 1 - 1])
	  ++(i);
	else
	  {
	    mm = lgrgn1->fields[i - 1];
	    lgrgn1->fields[i - 1] = lgrgn1->fields[i + 1 - 1];
	    lgrgn1->fields[i + 1 - 1] = mm;
	    mm = lgrgn1->perm[i - 1];
	    lgrgn1->perm[i - 1] = lgrgn1->perm[i + 1 - 1];
	    lgrgn1->perm[i + 1 - 1] = mm;
	    if (i == 1)
	      ++(i);
	    else
	      --(i);
	  }
      lgrgn1 = lgrgn1->next;
    }
  while (lgrgn1 != NULL);

  if (check)
    {
      mLine = nLine;
      lgrgn1 = lgrgn;		/* check1 */
      do
	{
	  nLine--;
	  lgrgn2 = lgrgn1->next;
	  while (lgrgn2 != NULL)
	    {
	      if ((lgrgn1->fields[0] == lgrgn2->fields[0]) &&
		  (lgrgn1->fields[1] == lgrgn2->fields[1]) &&
		  (lgrgn1->fields[2] == lgrgn2->fields[2]) &&
		  (lgrgn1->fields[3] == lgrgn2->fields[3]))
		{
		  errorMessage ("P1,P2,P3,P4", "duplicate vertex");
		  return 0;
		}
	      lgrgn2 = lgrgn2->next;
	    }
	  lgrgn1 = lgrgn1->next;
	}
      while (lgrgn1 != NULL);

      nLine = mLine;
      lgrgn1 = lgrgn;		/* check2 */
      do
	{
	  nLine--;
	  for (i = 0; i < 4; i++)
	    {
	      f_copy[i] = lgrgn1->fields[i];
	      if (f_copy[i] != 0)
		{
		  mm = ghostmother (f_copy[i]);
		  f_copy[i] = prtclbase[mm - 1].anti + f_copy[i] - mm;
		}
	    }
	  i = 1;
	  while (i < 4)
	    if (f_copy[i - 1] >= f_copy[i])
	      ++(i);
	    else
	      {
		mm = f_copy[i - 1];
		f_copy[i - 1] = f_copy[i];
		f_copy[i] = mm;
		if (i == 1)
		  ++(i);
		else
		  --(i);
	      }
	  lgrgn2 = lgrgn;
	  while ((lgrgn2 != NULL) && ((f_copy[0] != lgrgn2->fields[0]) ||
				      (f_copy[1] != lgrgn2->fields[1]) ||
				      (f_copy[2] != lgrgn2->fields[2]) ||
				      (f_copy[3] != lgrgn2->fields[3])))
	    {
	      lgrgn2 = lgrgn2->next;
	    }
	  if (lgrgn2 == NULL)
	    {
	      char sss[10];
	      strcpy (sss, "");
	      for (i = 0; i < 3; i++)
		{
		  strcat (sss, prtclbase[lgrgn1->fields[i] - 1].name);
		  strcat (sss, " ");
		}
	      if (lgrgn1->fields[3] != 0)
		strcat (sss, prtclbase[lgrgn1->fields[3] - 1].name);
	      errorMessage ("P1,P2,P3,P4", scat ("conjugated vertex for %s not found", sss));
	      return 0;
	    }
	  lgrgn1 = lgrgn1->next;
	}
      while (lgrgn1 != NULL);
    }
  return 1;
}


static void 
filldecaylist (void)
{
  algvertptr lgrgn1;
  int i, j, k, n;
  particleNumType pn[5], cc[3];
  decaylink kk, qq;

  lgrgn1 = lgrgn;
  do
    {
      for (i = 1; i <= 4; i++)
	{
	  pn[i - 1] = ghostmother (lgrgn1->fields[i - 1]);
	  if (pn[i - 1] != 0)
	    pn[i - 1] = prtclbase[pn[i - 1] - 1].anti;
	}
      pn[4] = 0;

      for (i = 1; i <= 4; i++)
	if (pn[i - 1] != pn[i + 1 - 1] && pn[i - 1] != 0)
	  {
	    j = 1;
	    for (k = 1; k <= 4; k++)
	      if (k != i)
		{
		  cc[j - 1] = pn[k - 1];
		  j++;
		}
	    n = prtclbase[pn[i - 1] - 1].anti;

	    if (prtclbase[n - 1].top == NULL)
	      {
		prtclbase[n - 1].top = (decaylink) m_alloc (sizeof (modeofdecay));
		prtclbase[n - 1].top->next = NULL;
		memcpy (prtclbase[n - 1].top->part, cc, 3 * sizeof (particleNumType));
	      }
	    else
	      {
		qq = prtclbase[n - 1].top;
		while (1)
		  {
		    k = 1;
		    while (k < 4 && qq->part[k - 1] == cc[k - 1])
		      k++;
		    if (k == 4)
		      goto exi;
		    if (qq->part[k - 1] > cc[k - 1])
		      {
			kk = (decaylink) m_alloc (sizeof (modeofdecay));
			kk->next = qq->next;
			qq->next = kk;
			memcpy (kk->part, qq->part, 3 * sizeof (particleNumType));
			memcpy (qq->part, cc, 3 * sizeof (particleNumType));
			goto exi;
		      }
		    if (qq->next == NULL)
		      {
			kk = (decaylink) m_alloc (sizeof (modeofdecay));
			kk->next = qq->next;
			qq->next = kk;
			memcpy (kk->part, cc, 3 * sizeof (particleNumType));
			goto exi;
		      }
		    qq = qq->next;
		  }
	      exi:;
	      }
	  }
      lgrgn1 = lgrgn1->next;
    }
  while (lgrgn1 != NULL);
}


static int 
isDublicateName (shortstr s)
{
  shortstr name;
  int k;

  for (k = 0; k < n_cpart; k++)
    {
      strcpy (name, cpartbase[k].name);
      if (!strcmp (name, s))
	return 0;
    }
  return 1;
}


static int 
isOrigName (shortstr s)
{
  shortstr name;
  int k;

  for (k = 0; k < nparticles; k++)
    {
      strcpy (name, prtclbase[k].name);
      if (!strcmp (name, s))
	return 0;
    }
  return 1;
}


static int 
cpart_check (shortstr str, char *wpart, int *cont)
{
  int i = 0;
  int n = 0;
  int num = 0;
  int len;

  len = strlen (str);
  if (!isalpha (str[0]) && '~'!=str[0])
    return 0;

  for (i = 0; i < len; i++)
    {
      if (str[i] == ',' || str[i] == ' ')
	{
	  wpart[n] = 0;
	  cont[num] = locateinbase (wpart);
	  if (n && !(cont[num])) return 0;
	  num++;
	  n = 0;
	}
      else
	{
	  wpart[n] = str[i], n++;
	}
    }

  if (0 != n)
    {
      wpart[n]=0;
      cont[num] = locateinbase (wpart);
      if (len && !(cont[num])) return 0;
      num++;
    }
  return num;
}


int 
readcpart (int check)
{
  shortstr cpart_tmp;
  shortstr name, wrong;
  midstr s;
  linelist ln;
  int nv = 0;
  int ncp;
  int content[1024];

  ln = cpart_tab.strings;
  while (ln)
    {
      ln = ln->next;
      nv++;
    }
    
  strcpy(tabName, cpart_tab.headln);
  if (cpartbase)
    {
      free (cpartbase);
    }
  cpartbase = m_alloc (nv * sizeof (*cpartbase));

  n_cpart = 0;

  ln = cpart_tab.strings;
  nLine = 1;
  while (ln)
    {
      strcpy (s, ln->line);
      sscanf (s, "%[^|]%*c%[^|]", name, cpart_tmp);
      trim (name);
      trim (cpart_tmp);

      if (check && (!isDublicateName (name)))
	{
	  errorMessage ("Abr", scat ("identical composite particles '%s'", name));
	  goto errExi;
	}

      if (check && (!isOrigName (name)))
	{
	  errorMessage ("Abr", scat ("The name %s coincides with the name of model particle", name));
	  goto errExi;
	}

      ncp = cpart_check (cpart_tmp, wrong, content);
      if (ncp > 40)
	{
	  errorMessage ("elementary particles", "too many model particles");
	  goto errExi;
	}
      if (check && !ncp)
	{
	  errorMessage ("elementary particles", scat ("wrong model particle '%s'", wrong));
	  goto errExi;
	}

      n_cpart++;
      if (n_cpart > nv)
        {
          cpartbase = realloc (cpartbase, n_cpart * sizeof (comppart));
	}
             cpartbase[n_cpart - 1].how = ncp;
      strcpy (cpartbase[n_cpart - 1].name,  name);
      strcpy (cpartbase[n_cpart - 1].cpart, cpart_tmp);

      if (n_cpart > USHRT_MAX - 1)
	{
	  errorMessage ("Name", "too many composite particles");
	  goto errExi;
	}

      ln = ln->next;
      nLine++;
    }
  return 1;

errExi:
  free (cpartbase);
  cpartbase = NULL;
  return 0;
}


int
readhadrons(int check)
{
  int i;
  int nbp;
  int nv = 0;
  int content[1024];
  double ms;
  shortstr name;
  shortstr mass;
  shortstr bm;
  shortstr wrong;
  midstr s;
  linelist ln;

  ln = hadron_tab.strings;
  while (ln)
    {
      ln = ln->next;
      nv++;
    }
  strcpy(tabName, hadron_tab.headln);
  if (hadronbase)
    {
      free (hadronbase);
    }

  hadronbase = m_alloc (nv * sizeof (*hadronbase));

  n_hadron = 0;

  ln = hadron_tab.strings;
  nLine = 1;
  while (ln)
    {
      shortstr namefrommdl;
      strcpy (s, ln->line);
      sscanf (s, "%[^|]%*c%[^|]%*c%[^|]%*c%[^|]", name, mass, bm, namefrommdl);
      trim (name);
      trim (mass);
      trim (bm);
      trim (namefrommdl);

      if (check && (!isOrigName (name)))
	{
	  errorMessage ("Abr", scat ("The name %s coincides with the name of model particle", name));
	  goto errExi;
	}

      if (1 != sscanf (mass, "%lf", &ms))
	{
	  errorMessage ("Particle mass", scat ("Unknown string for mass of particle %s", name));
	  goto errExi;
	}
      
      nbp = cpart_check (bm, wrong, content);
      if (nbp > 100)
	{
	  errorMessage ("elementary particles", "too many model particles");
	  goto errExi;
	}
      if (nbp)
	{
          int set, mem;
	  shortstr name1;
          n_hadron++;
          if (n_hadron > nv)
            {
              hadronbase = re_alloc (hadronbase, n_hadron * sizeof (hadron));
	    }
          strcpy (hadronbase[n_hadron - 1].name, name);
	  hadronbase[n_hadron - 1].mass = ms;
          if (get_sf_info (namefrommdl, name1, &set, &mem)) {
	    strcpy (hadronbase[n_hadron - 1].sf_name, name1);
	    hadronbase[n_hadron - 1].sf_set = set;
	    hadronbase[n_hadron - 1].sf_mem = mem;
	  }
          hadronbase[n_hadron - 1].how = nbp;
	  for(i=0;i<nbp;i++) hadronbase[n_hadron - 1].parton[i] = content[i];

          if (n_hadron > USHRT_MAX - 1)    /*Const SHRT_MAX = 32767 (from limits.h)*/
	    {
	      errorMessage ("Name", "too many beam particles");
	      goto errExi;
	    }
        }
      else
	{
/*	  errorMessage ("elementary particles", scat ("wrong model particle '%s'", wrong));
	  goto errExi;
*/
	}
      ln = ln->next;
      nLine++;
    }
  return 1;

errExi:
  free (hadronbase);
  hadronbase = NULL;
  return 0;
}

int
readstrfuns(int check)
{
  int nv = 0;
  int the_set;
  int the_mem;
  double nm;
  shortstr name;
  shortstr num;
  shortstr set;
  shortstr mem;
  midstr s;
  linelist ln;

  ln = strfun_tab.strings;
  while (ln)
    {
      ln = ln->next;
      nv++;
    }
  strcpy(tabName, strfun_tab.headln);

  if (strfunbase)
    {
      free (strfunbase);
    }

  strfunbase = m_alloc (nv * sizeof (*strfunbase));

  n_strfun = 0;

  ln = strfun_tab.strings;
  nLine = 1;
  while (ln)
    {
      strcpy (s, ln->line);
      strcpy (mem, "0");
      sscanf (s, "%[^|]%*c%[^|]%*c%[^|]%*c%[^|]", name, num, set, mem);
      trim (name);
      trim (num);
      trim (set);
      trim (mem);

      if (1 != sscanf (num, "%lf", &nm) || 1 != sscanf (set, "%d", &the_set) || 1 != sscanf (mem, "%d", &the_mem))
	{
	  errorMessage ("Abr", scat ("Unparsed string in the structure function table:\n %s", s));
	  goto errExi;
	}

      n_strfun++;
      if (n_strfun > nv)
        {
          strfunbase = re_alloc (strfunbase, n_strfun * sizeof (strfun));
        }
      strcpy (strfunbase[n_strfun - 1].name, name);
      strfunbase[n_strfun - 1].num = nm;
      strfunbase[n_strfun - 1].set = the_set;
      strfunbase[n_strfun - 1].mem = the_mem;

      if (n_strfun > USHRT_MAX - 1)
        {
          errorMessage ("Name", "too many structure functions");
          goto errExi;
        }
      ln = ln->next;
      nLine++;
    }
  return 1;

errExi:
  free (strfunbase);
  strfunbase = NULL;
  return 0;
}


int 
loadModel (int check)
{
  int num = getModelNumberSymb();
  errorText[0] = 0;
  if ((!check) && (lastModel == num))
    return 1;

  if (!readvars (check))
    return 0;
  if (!readparticles (check))
    return 0;
  if (!readlagrangian (check))
    return 0;
  filldecaylist ();

  if (!readcpart (check))
    return 0;
  if (!readhadrons (check))
    return 0;
  if (!readstrfuns (check))
    return 0;

  lastModel = num;
  return 1;
}


int 
readModelFiles (int l, char * path)
{
  FILE * test = NULL;
  char tmpPath[STRSIZ];
  
  cleartab (&vars_tab);
  cleartab (&func_tab);
  cleartab (&prtcls_tab);
  cleartab (&lgrng_tab);
  cleartab (&cpart_tab);
  cleartab (&hadron_tab);
  cleartab (&strfun_tab);

  lastModel = 0;

  strcpy (tmpPath, scat ("%s%s%cvars%d.mdl", pathtouser, path, f_slash, l));
  test = fopen(tmpPath, "r");
  if (test == NULL) return 0;
  readtable (&vars_tab, tmpPath);

  strcpy (tmpPath, scat ("%s%s%cfunc%d.mdl", pathtouser, path, f_slash, l));
  test = fopen(tmpPath, "r");
  if (test == NULL) return 0;
  readtable (&func_tab, tmpPath);

  strcpy (tmpPath, scat ("%s%s%cprtcls%d.mdl", pathtouser, path, f_slash, l));
  test = fopen(tmpPath, "r");
  if (test == NULL) return 0;
  readtable (&prtcls_tab, tmpPath);

  strcpy (tmpPath, scat ("%s%s%clgrng%d.mdl", pathtouser, path, f_slash, l));
  test = fopen(tmpPath, "r");
  if (test == NULL) return 0;
  readtable (&lgrng_tab, tmpPath);

  strcpy (tmpPath, scat ("%s%s%ccpart%d.mdl", pathtouser, path, f_slash, l));
  test = fopen(tmpPath, "r");
  if (test == NULL) return 0;
  readtable (&cpart_tab, tmpPath);

  strcpy (tmpPath, scat ("%s%s%cbeams.mdl", pathtouser, path, f_slash));
  test = fopen(tmpPath, "r");
  if (test == NULL) return 0;
  readtable (&hadron_tab, tmpPath);

  strcpy (tmpPath, scat ("%s%s%cstrfuns.mdl", pathtouser, path, f_slash));
  test = fopen(tmpPath, "r");
  if (test == NULL) return 0;
  readtable (&strfun_tab, tmpPath);

return 1;
}

int 
copyModelFiles (int l, char * path)
{
  int result = 1;
  char tmpComm[STRSIZ];
  
  strcpy (tmpComm, scat ("cp %smodels%cvars%d.mdl   %sresults%c%s%c.", pathtouser, f_slash, l, pathtouser, f_slash, path, f_slash));
  result = result * system (tmpComm);
  strcpy (tmpComm, scat ("cp %smodels%cfunc%d.mdl   %sresults%c%s%c.", pathtouser, f_slash, l, pathtouser, f_slash, path, f_slash));
  result = result * system (tmpComm);
  strcpy (tmpComm, scat ("cp %smodels%cprtcls%d.mdl %sresults%c%s%c.", pathtouser, f_slash, l, pathtouser, f_slash, path, f_slash));
  result = result * system (tmpComm);
  strcpy (tmpComm, scat ("cp %smodels%clgrng%d.mdl  %sresults%c%s%c.", pathtouser, f_slash, l, pathtouser, f_slash, path, f_slash));
  result = result * system (tmpComm);
  strcpy (tmpComm, scat ("cp %smodels%ccpart%d.mdl  %sresults%c%s%c.", pathtouser, f_slash, l, pathtouser, f_slash, path, f_slash));
  result = result * system (tmpComm);
  strcpy (tmpComm, scat ("cp %smodels%cbeams.mdl    %sresults%c%s%c.", pathtouser, f_slash, pathtouser, f_slash, path, f_slash));
  result = result * system (tmpComm);
  strcpy (tmpComm, scat ("cp %smodels%cstrfuns.mdl  %sresults%c%s%c.", pathtouser, f_slash, pathtouser, f_slash, path, f_slash));
  result = result * system (tmpComm);
  
  return result;
}
