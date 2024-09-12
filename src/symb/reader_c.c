/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/getmem.h"
#include "service2/include/syst.h"

#include "symb/include/physics.h"
#include "symb/include/procvar.h"
#include "symb/include/reader_c.h"

FILE *ext_h = NULL;

static void *
bact5 (char ch, void *mm1, void *mm2)
{
  char *m1;
  char *m2;
  char *ans;
  char r_n, p_m;

  m1 = (char *) mm1;
  m2 = (char *) mm2;

  if (ch == '+' || ch == '-')
    p_m = 'P';
  else
    p_m = 'M';

  r_n = 'R';

  if (m1[0] == 'M' || ch == '+' || ch == '-')
    lShift (m1, 3);
  else
    {
      lShift (m1, 2);
      m1[0] = '(';
      strcat (m1, ")");
    }
  if ((m2[0] == 'M' || ch == '+') && ch != '/')
    lShift (m2, 3);
  else
    {
      lShift (m2, 2);
      m2[0] = '(';
      strcat (m2, ")");
    }

  ans = m_alloc (strlen (m1) + strlen (m2) + 30);

  switch (ch)
    {
    case '+':
      if (m1[0] == '-')
	sprintf (ans, "%c%c|%s%s", p_m, r_n, m2, m1);
      else if (m2[0] == '-')
	sprintf (ans, "%c%c|%s%s", p_m, r_n, m1, m2);
      else
	sprintf (ans, "%c%c|%s+%s", p_m, r_n, m1, m2);
      break;

    case '-':
      if (m2[0] == '-')
	sprintf (ans, "%c%c|%s+%s", p_m, r_n, m1, m2 + 1);
      else
	sprintf (ans, "%c%c|%s-%s", p_m, r_n, m1, m2);
      break;


    case '*':
      if (m2[0] != '-')
	sprintf (ans, "%c%c|%s*%s", p_m, r_n, m1, m2);
      else if (m1[0] != '-')
	sprintf (ans, "%c%c|%s*%s", p_m, r_n, m2, m1);
      else
	sprintf (ans, "%c%c|%s*%s", p_m, r_n, m1 + 1, m2 + 1);
      break;

    case '/':
      if (m2[0] != '-')
	sprintf (ans, "%c%c|%s/%s", p_m, r_n, m1, m2);
      else
	{
	  if (m1[0] == '-')
	    sprintf (ans, "%c%c|%s/%s", p_m, r_n, m1 + 1, m2 + 1);
	  else
	    sprintf (ans, "%c%c|-%s/%s", p_m, r_n, m1, m2 + 1);
	}
      break;

    case '^':
      sprintf (ans, "%c%c|pow(%s,%s)", p_m, r_n, m1, m2);

    }				/* Case */
  return (void *) ans;
}


static void *
uact5 (char *ch, void *mm)
{
  char *m, *ans;
  m = (char *) mm;
  ans = m_alloc (strlen (m) + 30);

  if (strcmp (ch, "-") == 0)
    {
      if (m[0] == 'M')
	{
	  if (m[3] == '-')
	    sprintf (ans, "M%c|%s", m[1], m + 4);
	  else
	    sprintf (ans, "M%c|-%s", m[1], m + 3);
	}
      else
	sprintf (ans, "M%c|-(%s)", m[1], m + 3);
    }

  if (strcmp (ch, "sqrt") == 0)
    {
      if (m[1] == 'N')
	sprintf (ans, "MR|sqrt_e((double)(%s),&err)", m + 3);
      else
	sprintf (ans, "MR|sqrt_e(%s,&err)", m + 3);
    }
  return (void *) ans;
}

void *
act_c (char *name, int n, void **args)
{
  int l, i;
  char *ans;

  if (!isalpha (name[0]))
    {
      if (n == 1)
	return uact5 (name, args[0]);
      if (n == 2)
	return bact5 (name[0], args[0], args[1]);
    }
  else
    {
      if (strcmp (name, "sqrt") && strcmp (name, "sin") && strcmp (name, "cos")
	  && strcmp (name, "tan") && strcmp (name, "asin") && strcmp (name, "acos")
	  && strcmp (name, "atan") && strcmp (name, "exp") && strcmp (name, "log")
	  && strcmp (name, "fabs") && strcmp (name, "atan2") && strcmp (name, "if")
	)
	{
	  fprintf (ext_h, " extern double %s(", name);
	  for (i = 1; i < n; i++)
	    fprintf (ext_h, "double,");
	  fprintf (ext_h, "double);\n");
	}
    }
  l = n + 10 + strlen (name);
  for (i = 0; i < n; i++)
    l += strlen ((char *) args[i]);
  ans = m_alloc (l);
  if (!strcmp (name, "if") && n == 3)
    sprintf (ans, "MR|(%s>0 ? %s : %s)", (char *) args[0] + 3, (char *) args[1] + 3,
	     (char *) args[2] + 3);
  else
    {
      sprintf (ans, "MR|%s(", name);
      for (i = 0; i < n; i++)
	{
	  strcat (ans, (char *) args[i] + 3);
	  strcat (ans, ",");
	}
      ans[strlen (ans) - 1] = ')';
    }
  return ans;
}

void *
rd_c (char *s)
{
  char *p;
  int l;
  p = m_alloc (40);
  if ('0' <= s[0] && s[0] <= '9')
    sprintf (p, "MN|(double)%s", s);
  else
    {
      for (l = 1; l <= nmodelvar; l++)
	{
	  if (!strcmp (s, modelvars[l].varname))
	    {
	      sprintf (p, "MR|%s", vararr[l].alias);
	      return (void *) p;
	    }
	}
    }
  return (void *) p;
}
