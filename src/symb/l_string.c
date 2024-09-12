/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <limits.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/getmem.h"
#include "service2/include/syst.h"

#include "physics.h"
#include "procvar.h"
#include "optimise.h"
#include "out_service.h"
#include "l_string.h"

#define buffsize MIN(1290,STRSIZ-3)

typedef struct llongstr
  {
    int len;
    char txt[buffsize + 1];
  }
llongstr;

typedef llongstr *longstrptr;



typedef struct degreerec
  {
    struct degreerec *next;
    char name[10];
    int deg;
  }
degreerec;

static degreerec **degnamesptr = NULL;

static int degnamecount;
static int tmpNameNum;
static int c_style;
static int tmpNameMax = 0;

static void (*writelongstr) (char *name, longstrptr longs);


void 
initdegnames (void)
{
  int k;
  degnamesptr = (degreerec **) m_alloc (sizeof (degreerec **) * (nProcessVar));
  for (k = 0; k < nProcessVar; k++)
    degnamesptr[k] = NULL;
  degnamecount = 0;
}

int 
cleardegnames (void)
{
  int k;
  degreerec *p, *q;
  if (!degnamesptr)
    return 0;
  for (k = 0; k < nProcessVar; k++)
    {
      q = degnamesptr[k];
      while (q != NULL)
	{
	  p = q;
	  q = p->next;
	  free (p);
	}
    }
  free (degnamesptr);
  degnamesptr = NULL;
  return degnamecount;
}


static void 
writelongstr_c (char *name, longstrptr longs)
{
  writeF ("%s=", name);
  if (longs == NULL || longs->len == 0)
    {
      writeF ("0;\n");
      return;
    }

  writeF ("%.*s;\n", longs->len, longs->txt);
}


static void 
writelongstr_f (char *name, longstrptr longs)
{
  writeF ("      %s=", name);
  if (longs == NULL || longs->len == 0)
    {
      writeF ("0\n");
      return;
    }
  writeF ("%.*s\n", longs->len, longs->txt);
}



static void 
addstring (longstrptr longs, char *s)
{
  int i, l, ll;
  char name[STRSIZ];
  ll = longs->len;
  l = strlen (s);
  if (ll + l > buffsize)
    {
      if (c_style)
	sprintf (name, "tmp[%d]", tmpNameNum++);
      else
	sprintf (name, "tmp%d", ++tmpNameNum);
      writelongstr (name, longs);
      longs->len = 0;
      addstring (longs, name);
      ll = longs->len;
    }

  for (i = 0; i < l; i++)
    longs->txt[ll + i] = s[i];
  longs->len += l;
}


static pointer 
gorner (char *s, longstrptr pmult, longstrptr psum)
{
  char name[STRSIZ], name2[STRSIZ];
  longstrptr ans;
  pointer pchange;

  if (pmult == NULL)
    return (pointer) psum;
  ans = (longstrptr) m_alloc (sizeof (longstr));
  ans->len = 0;
  addstring (ans, s);

  if (3 + ans->len + pmult->len > buffsize)
    {
      if (c_style)
	sprintf (name, "tmp[%d]", tmpNameNum++);
      else
	sprintf (name, "tmp%d", ++tmpNameNum);
      writelongstr (name, pmult);
      addstring (ans, scat ("*%s", name));
    }
  else
    {
      if (pmult->txt[0] == '+')
	{
	  pmult->txt[0] = '(';
	  addstring (ans, "*");
	}
      else
	addstring (ans, "*(");

      memcpy (&(ans->txt[ans->len]), pmult->txt, pmult->len);
      ans->len += pmult->len;
      addstring (ans, ")");
    }
  free (pmult);

  if (psum == NULL)
    return (pointer) ans;
  if (ans->len + psum->len > buffsize)
    {
      if (c_style)
	sprintf (name, "tmp[%d]", tmpNameNum++);
      else
	sprintf (name, "tmp%d", ++tmpNameNum);
      if (ans->len > psum->len)
	{
	  pchange = (pointer) ans;
	  ans = psum;
	  psum = (longstrptr) pchange;
	}
      writelongstr (name, psum);
      if (ans->len + strlen (name) >= buffsize)
	{
	  if (c_style)
	    sprintf (name2, "tmp[%d]", tmpNameNum++);
	  else
	    sprintf (name2, "tmp%d", ++tmpNameNum);
	  writelongstr (name2, ans);
	  ans->len = 0;
	  addstring (ans, scat ("+%s+%s", name, name2));
	}
      else
	addstring (ans, scat ("+%s", name));

    }
  else
    {
      memcpy (&(ans->txt[ans->len]), psum->txt, psum->len);
      ans->len += psum->len;
    }
  free (psum);
  return (pointer) ans;
}


static char *
writevardeg (int nv, int deg)
{
  degreerec *p;
  static char namest[21];
  if (deg == 1)
    return vararr[nv].alias;
  p = degnamesptr[nv];
  while (p != NULL)
    if (p->deg == deg)
      return p->name;
    else
      p = p->next;
  p = (degreerec *) m_alloc (sizeof (degreerec));
  p->deg = deg;
  if (c_style)
    sprintf (namest, "S[%d]", degnamecount++);
  else
    sprintf (namest, "S%d", ++degnamecount);
  strcpy (p->name, namest);
  p->next = degnamesptr[nv];
  degnamesptr[nv] = p;
  if (c_style)
    {
      int k;
      writeF ("%s=%s", namest, vararr[nv].alias);
      for (k = 1; k < deg; k++)
	writeF ("*%s", vararr[nv].alias);
      writeF (";\n");
    }
  else
    f_printf (outFile, "      %s=%s**%d\n", namest, vararr[nv].alias, deg);
  return namest;
}


static pointer 
smpl_emit (varptr ex)
{
  longstrptr ans;
  char s[STRSIZ];
  int k, bt, nv, deg;
  int star;
  varptr ex_, exbeg;

  if (ex == NULL)
    return NULL;

  if (ex->sgn == '-')
    {
      ex_ = ex;

      while (ex_->next != NULL && ex_->next->sgn == '-')
	ex_ = ex_->next;

      if (ex_->next != NULL)
	{
	  exbeg = ex_->next;
	  ex_->next = exbeg->next;
	  exbeg->next = ex;
	  ex = exbeg;
	}

    }

  ans = (longstrptr) m_alloc (sizeof (llongstr));
  ans->len = 0;
  while (ex != NULL)
    {
      sprintf (s, "%c", ex->sgn);
      star = (strcmp ((ex->coef)->name, "1") != 0);
      if (star || short_strlen (ex->vars) == 0)
	strcat (s, (ex->coef)->name);
      if (short_strlen (ex->vars) != 0)
	{
	  bt = (ex->vars[0]);
	  deg = 1;
	  for (k = 2; k <= short_strlen (ex->vars); k++)
	    {
	      nv = (ex->vars[k - 1]);
	      if (bt != nv)
		{
		  if (star)
		    strcat (s, "*");
		  else
		    star = TRUE;
		  strcat (s, writevardeg (bt, deg));
		  deg = 1;
		  bt = nv;
		}
	      else
		++(deg);
	    }
	  if (star)
	    strcat (s, "*");
	  else
	    star = TRUE;
	  strcat (s, writevardeg (bt, deg));
	}
      addstring (ans, s);
      ex = ex->next;
    }
  return (pointer) ans;
}

static pointer 
v_gorner (int ch, int deg, pointer pmult, pointer psum)
{
  char b[STRSIZ];
  sprintf (b, "+%s", writevardeg (ch, deg));
  return gorner (b, pmult, psum);
}

static pointer 
c_gorner (infoptr i, pointer pmult, pointer psum)
{
  char b[STRSIZ];
  sprintf (b, "+%s", i->name);
  return gorner (b, pmult, psum);
}

void 
fortwriter (char *name, varptr fortformula)
{
  longstrptr tmp;
  tmpNameNum = 0;
  tmp = (longstrptr) emitexpr (fortformula, smpl_emit, v_gorner, c_gorner);
  writelongstr (name, tmp);
  if (tmp != NULL)
    free (tmp);
  if (tmpNameMax < tmpNameNum)
    tmpNameMax = tmpNameNum;
}


int 
write_const (void)
{
  infoptr i;
  int firstVarTmp;
  int constcount = 0;


  firstVarTmp = firstVar;
  firstVar = 1;

  for (i = info; i; i = i->next)
    {
      if (i->consttype == numb)
	{
	  if (ABS (i->ival) >= LONG_MAX)
	    sprintf (i->name, "%" NUM_STR ".", i->ival);
	  else
	    sprintf (i->name, "%" NUM_STR, i->ival);
	}
      else
	sprintf (i->name, "C[%d]", constcount++);
    }

  for (i = info; i; i = i->next)
    {
      if (i->consttype != numb)
	fortwriter (i->name, (pointer) i->const_);
    }

  firstVar = firstVarTmp;

  return constcount;
}				/*  WriteConst  */


void 
initfortwriting (char style)
{
  c_style = (style == 'c');
  if (c_style)
    writelongstr = &writelongstr_c;
  else
    writelongstr = &writelongstr_f;
  initinfo ();
}

int gettmpNameMax (void) {
  return tmpNameMax;
}

void settmpNameMax (int num) {
  tmpNameMax = num;
}
