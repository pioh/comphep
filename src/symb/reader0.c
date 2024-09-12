/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"

#include "physics.h"
#include "reader0.h"

char momsubst[9] = {0};
char indsubst[9] = {0};
int r_reading0 = FALSE;

static void *
bact0 (char ch, pointer mm1, pointer mm2)
{
  char *m1, *m2, *ans;
  int sgn;

  if (r_reading0 && (ch == '*'))
    {
      m1 = (char *) mm2;
      m2 = (char *) mm1;
    }
  else
    {
      m1 = (char *) mm1;
      m2 = (char *) mm2;
    }
  if (ch == '+' || ch == '-')
    {
      lShift (m1, 2);
      lShift (m2, 2);
    }
  else
    {
      if (m1[0] == 'P' || ch == '^')
	{
	  lShift (m1, 1);
	  m1[0] = '(';
	  strcat (m1, ")");
	}
      else
	lShift (m1, 2);
      if (m2[0] == 'P' || ch == '^')
	{
	  lShift (m2, 1);
	  m2[0] = '(';
	  strcat (m2, ")");
	}
      else
	lShift (m2, 2);
    }

  ans = (char *) m_alloc (strlen (m1) + strlen (m2) + 8);
  switch (ch)
    {
    case '+':
      if (m2[0] == '-')
	sprintf (ans, "P|%s%s", m1, m2);
      else
	sprintf (ans, "P|%s+%s", m1, m2);
      break;

    case '-':
      if (m2[0] == '-')
	sprintf (ans, "P|%s+%s", m1, m2 + 1);
      else
	sprintf (ans, "P|%s-%s", m1, m2);
      break;

    case '*':
      sgn = 1;
      if (m1[0] == '-')
	{
	  lShift (m1, 1);
	  sgn = -sgn;
	}
      if (m2[0] == '-')
	{
	  lShift (m2, 1);
	  sgn = -sgn;
	}
      if (sgn == 1)
	sprintf (ans, "M|%s*%s", m1, m2);
      else
	sprintf (ans, "M|-%s*%s", m1, m2);
      break;

    case '.':
      if (m2[0] != '-')
	sprintf (ans, "M|%s.%s", m1, m2);
      else if (m1[0] != '-')
	sprintf (ans, "M|%s.%s", m2, m1);
      else
	{
	  lShift (m1, 1);
	  lShift (m2, 1);
	  sprintf (ans, "M|%s.%s", m1, m2);
	}
      break;

    case '^':
      sprintf (ans, "M|%s^%s", m1, m2);
    }
  return (pointer) ans;
}


static void *
uact0 (char *ch, pointer mm)
{
  char *m;
  char *ans;

  m = (char *) mm;
  ans = (char *) m_alloc (strlen (m) + 10);

  if (strcmp (ch, "-") == 0)
    {
      if (m[0] == 'M')
	{
	  if (m[2] == '-')
	    sprintf (ans, "M|%s", m + 3);
	  else
	    sprintf (ans, "M|-%s", m + 2);
	}
      else
	sprintf (ans, "M|-(%s)", m + 2);
    }
  else if (strcmp (ch, "G") == 0)
	{
	  if (r_reading0)
	    sprintf (ans, "M|-G(ln,%s)", m + 2);
	  else
	    sprintf (ans, "M|G(ln,%s)", m + 2);
	}
  else
    {
       printf("***** uact0: invalid symbol '%s'. Fatal error", ch);
       exit(99);
    }
  return (pointer) ans;
}

void *
act_rcode (char *ch, int n, void **args)
{
  if (n == 1)
    return uact0 (ch, args[0]);
  if (n == 2)
    return bact0 (ch[0], args[0], args[1]);
  if (n == 4 && !strcmp (ch, "eps"))
    {
      int l = 15 + strlen (args[0]) + strlen (args[1]) + strlen (args[2]) + strlen (args[3]);
      char *ans = (char *) m_alloc (l);
      sprintf (ans, "M|eps(%s,%s,%s,%s)", (char *) args[0] + 2, (char *) args[1] + 2,
	       (char *) args[2] + 2, (char *) args[3] + 2);
      return ans;
    }
  fprintf (stderr, "***** act_rcode: Reach the end. No return void* val.!?\n");
  return 0;
}


void *
rd_rcode (char *s)
{
  char *p;
  int num;

  p = (char *) m_alloc (12);
  p[0] = 0;
  if (strlen (s) == 2 && s[1] > '0' && s[1] <= '9')
    {
      switch (s[0])
	{
	case 'p':
	case 'P':
	  num = s[1] - '0';
	  num = momsubst[num - 1];
	  if (num > 0)
	    sprintf (p, "M|p%d", num);
	  else
	    sprintf (p, "M|-p%d", -num);
	  break;

	case 'm':
	  num = s[1] - '0';
	  num = indsubst[num - 1];
	  sprintf (p, "M|m%d", num);
	  break;
	case 'M':
	  num = s[1] - '0';
	  num = indsubst[num - 1] - 1;
	  sprintf (p, "M|m%d", num);
	}
      if (strcmp (s, "G5") == 0)
	strcpy (p, "M|G(ln,A)");
    }
  if (!strlen (p))
    sprintf (p, "M|%s", s);
  return (void *) p;
}
