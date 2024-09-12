/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/

#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "chep_crt/include/crt_util.h"

#include "mssmlib/include/sugrac.h"
#include "mssmlib/include/fhf3.h"
#include "mssmlib/include/fhf3.h"
#include "mssmlib/include/g2f3_c.h"
#include "mssmlib/include/rotate_c.h"

#include "getmem.h"
#include "f_c.h"
#include "parser.h"
#include "read_func.h"

#include "width_func.h"

#define MAX_ARGS (30)

static int (*nameToVal) (char *, double *);
static int isAble;



double 
sort4 (double m1, double m2, double m3, double m4, double dn)
{
  int n;
  int i, f = 0;
  double m[4];
  n = (int) floor (dn - 0.5);
  m[0] = m1;
  m[1] = m2;
  m[2] = m3;
  m[3] = m4;

  do
    {
      f = 0;
      for (i = 0; i < 3; i++)
	if (fabs (m[i]) > fabs (m[i + 1]))
	  {
	    double tmp;
	    tmp = m[i];
	    m[i] = m[i + 1];
	    m[i + 1] = tmp;
	    f = 1;
	  }

    }
  while (f);

  return m[n];
}

static void *
rd_num_local (char *s)
{
  int rc;
  double * p = (double *) getmem_ (sizeof (double));

  if (!nameToVal)
    {
      fprintf (stderr, "Error in programming, nameToVal==NULL\n");
      exit (1);
    }

  rc = (*nameToVal) (s, p);
  if (!rc)
    {
      if (isdigit (*s))
	rderrcode = typemismatch;
      else
	rderrcode = unknownidentifier;
      return NULL;
    }
  else if (rc == -1)
    isAble = 0;
  return (void *) p;
}

static void *
act_num (char *ch, int n, void **args)
{
  int i;
  double pp[MAX_ARGS];

  if (!isAble)
    {
      *(double *) args[0] = 0.;
      return args[0];
    }

  if (n >= MAX_ARGS)
    {
      fprintf (stderr, "***** act_num: n=%d more MAX_ARGS=%d\n", n, MAX_ARGS);
      exit (99);
    }

  for (i = 0; i < n; i++)
    {
      pp[i] = *(double *) args[i];
    }

  switch (ch[0])
    {
    case '-':
      if (n == 1)
	pp[0] *= -1;
      else
	pp[0] -= pp[1];
      break;
    case '+':
      pp[0] += pp[1];
      break;
    case '*':
      pp[0] *= pp[1];
      break;
    case '/':
      if (pp[1] == 0.)
	rderrcode = naninoperation;
      else
	pp[0] /= pp[1];
      break;
    case '^':
      if (pp[1] == floor (pp[1]))
	pp[0] = pow_dl (pp[0], floor (pp[1]));
      else
	pp[0] = pow (pp[0], pp[1]);
      break;
    case '.':
      rderrcode = typemismatch;
      break;
    default:
      switch (n)
	{
	case 1:
	  if (!strcmp (ch, "sqrt"))
	    pp[0] = sqrt (pp[0]);
	  else if (!strcmp (ch, "sin"))
	    pp[0] = sin (pp[0]);
	  else if (!strcmp (ch, "cos"))
	    pp[0] = cos (pp[0]);
	  else if (!strcmp (ch, "tan"))
	    pp[0] = tan (pp[0]);
	  else if (!strcmp (ch, "asin"))
	    pp[0] = asin (pp[0]);
	  else if (!strcmp (ch, "acos"))
	    pp[0] = acos (pp[0]);
	  else if (!strcmp (ch, "atan"))
	    pp[0] = atan (pp[0]);
	  else if (!strcmp (ch, "exp"))
	    pp[0] = exp (pp[0]);
	  else if (!strcmp (ch, "log"))
	    pp[0] = log (pp[0]);
	  else if (!strcmp (ch, "fabs"))
	    pp[0] = fabs (pp[0]);
	  else
	    isAble = 0;
	  break;
	case 2:
	  if (!strcmp (ch, "atan2"))
	    pp[0] = atan2 (pp[0], pp[1]);
	  else if (!strcmp (ch, "feynhiggs"))
	    pp[0] = feynhiggs (pp[0], pp[1]);
	  else
	    isAble = 0;
	  break;
	case 3:
	  if (!strcmp (ch, "if"))
	    {
	      if (pp[0] > 0)
		pp[0] = pp[1];
	      else
		pp[0] = pp[2];
	    }
	  else
	    isAble = 0;
	  break;
	case 4:
	  if (!strcmp (ch, "pmass2"))
	    {
	      pp[0] = pmass2 (pp[0], pp[1], pp[2], pp[3]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	case 5:
	  if (!strcmp (ch, "rotate2")) 
	    {
	      pp[0] = rotate2 (pp[0], pp[1], pp[2], pp[3], pp[4]);
	    }
	  else if (!strcmp (ch, "sort4")) 
	    {
	      pp[0] = sort4 (pp[0], pp[1], pp[2], pp[3], pp[4]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 6:
	case 7:
	  if (!strcmp (ch, "g2f3"))
	    {
	      pp[0] = g2f3 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], 
	              pp[6]);
	    }
	  else if (!strcmp (ch, "pmass3"))
	    {
	      pp[0] = pmass3 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], 
	              pp[6]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 8:
	  if (!strcmp (ch, "rotate3"))
	    {
	      pp[0] = rotate3 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], 
	              pp[6], pp[7]);
	    }
	  else if (!strcmp (ch, "sugra"))
	    {
	      pp[0] = sugra (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], 
	              pp[6], pp[7]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 9:
	  if (!strcmp (ch, "gmsb"))
	    {
              pp[0] = gmsb (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], 
	              pp[6], pp[7], pp[8]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 10:
	case 11:
	  if (!strcmp (ch, "pmass4"))
	    {
	      pp[0] = pmass4 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], 
	              pp[6], pp[7], pp[8], pp[9], pp[10]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 12:
	  if(!strcmp (ch, "rotate4"))
	    {
	      pp[0] = rotate4 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5],
	              pp[6], pp[7], pp[8], pp[9], pp[10], pp[11]);
	    }
	  else if(!strcmp (ch, "feynhiggs1"))
	    {
	      pp[0] = feynhiggs1 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5],
	              pp[6], pp[7], pp[8], pp[9], pp[10], pp[11]);
	    }
	  else if(!strcmp (ch, "feynhiggs2"))
	    {
	      pp[0] = feynhiggs2 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5],
	              pp[6], pp[7], pp[8], pp[9], pp[10], pp[11]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 13:
	  if (!strcmp (ch, "fhf2"))
	    {
	      pp[0] = fhf2 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5],
	              pp[6], pp[7], pp[8], pp[9], pp[10], pp[11], pp[12]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 16:
	  if (!strcmp (ch, "pmass5"))
	    {
	      pp[0] = pmass5 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5],
	              pp[6], pp[7], pp[8], pp[9], pp[10], pp[11], pp[12], 
	              pp[13], pp[14], pp[15]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	case 17:
	  if (!strcmp (ch, "rotate5"))
	    {
	      pp[0] = rotate5 (pp[0], pp[1], pp[2], pp[3], pp[4], pp[5],
	              pp[6], pp[7], pp[8], pp[9], pp[10], pp[11], pp[12], 
	              pp[13], pp[14], pp[15], pp[16]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
        case 18:
	  if (!strcmp (ch, "width1"))
	    {
	      pp[0] = width1 (pp[0], pp[1]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
        case 19:
	  if (!strcmp (ch, "width2"))
	    {
	      pp[0] = width2 (pp[0], pp[1]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
        case 20:
	  if (!strcmp (ch, "width3"))
	    {
	      pp[0] = width3 (pp[0], pp[1]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
        case 21:
	  if (!strcmp (ch, "width4"))
	    {
	      pp[0] = width4 (pp[0], pp[1]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
        case 22:
	  if (!strcmp (ch, "width5"))
	    {
	      pp[0] = width5 (pp[0], pp[1]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
        case 23:
	  if (!strcmp (ch, "koeff1"))
	    {
	      pp[0] = width5 (pp[0], pp[1]);
	    }
	  else
	    {
	      fprintf (stderr, "***** act_num: Undefine function ch=%s\n", ch);
	      exit (99);
	    }
	  break;
	default:
	  isAble = 0;
	}
    }

  if (rderrcode)
    return NULL;
  if (!isAble)
    pp[0] = 0;
  *(double *) args[0] = pp[0];

  return args[0];
}

int 
calcExpression (char *s, int (*nameToDouble) (char *, double *), double *p)
{
  marktp heapbeg;
  double *r;
  isAble = 1;

  nameToVal = nameToDouble;

  mark_ (&heapbeg);
  r = (double *) readExpression (s, rd_num_local, act_num, NULL);
  if (!isAble)
    {
      rderrcode = unknownfunction;
    }
  else
    {
      if (!rderrcode)
	*p = *r;
      if (!rderrcode && !finite (*r))
	rderrcode = cannotevaluate;
    }
  release_ (&heapbeg);
  return rderrcode;
}
