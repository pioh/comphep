/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Victor Edneral
* ------------------------------------------------------------
*/
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "chep_limits.h"
#include "sets.h"

/* construct and return a set of the specified character values */
int 
setof (int i,...)
{
  va_list v;
  int k = 0;
  int r = 0;

  va_start (v, i);

  while (i != _E)
    {
      if (i > 15)
	{
	  fprintf (stderr, "Too large setof element!! Hit enter:");
	  getchar ();
	  exit (-1);
	}
      if (i == UpTo)
	i = va_arg (v, int);
      else
	k = i;
      for (; k <= i; k++)
	r |= 1 << i;
      i = va_arg (v, int);
    }
  va_end (v);
  return r;
}

/*
 *   int inset(int ex,int setrec)
 *      predicate returns TRUE if expression ex is a member of
 *      the set parameter
 *   (It is suppoused enum type setrec takes two bytes here; a,b,ex < 16) */
int 
inset (int a, int sp)
{
  if (a > 15)
    {
      fprintf (stderr, "Too large setof element!! Hit enter:");
      getchar ();
      exit (-1);
    }
  return ((1 << a) & sp) != 0;
}


/*
 *   setofbyte setofb(int a,int b,...,_E)
 *   int insetb(unsigned ex, setofbyte setrec)
 *   Here a,b,ex<256; setofbyte=unsigned[16].
 */
unsigned *
setofb (int i,...)
{
  va_list v;
  int k = 0;
  static setofbyte r;

  va_start (v, i);
  memset ((char *) r, 0, sizeof (setofbyte));
  while (i != _E)
    {
      if (i > 255)
	{
	  fprintf (stderr, "Too large setofb element!! Hit enter:");
	  getchar ();
	  return 0;
	}
      if (i == UpTo)
	i = va_arg (v, int);
      else
	k = i;
      for (; k <= i; k++)
	r[k >> 4] |= 1 << (k & 0xF);
      i = va_arg (v, int);
    }
  va_end (v);
  return r;
}

unsigned *
setofb_add1 (setofbyte a, int k)
{
  a[k >> 4] |= 1 << (k & 0xF);
  return a;
}


void 
setofb_zero (setofbyte sp)
{
  memset ((char *) sp, 0, sizeof (setofbyte));
}


int 
insetb (unsigned a, setofbyte sp)
{
  if (a > 255)
    {
      fprintf (stderr, "Too large setofb element!! Hit enter:");
      getchar ();
      exit (-1);
    }
  return ((1 << (a & 0xF)) & sp[a >> 4]) != 0;
}

void 
setofb_cpy (setofbyte dest, setofbyte source)
{
  memcpy ((char *) dest, (char *) source, sizeof (setofbyte));
}

unsigned *
setofb_uni (setofbyte a, setofbyte b)
{
  static setofbyte res;
  int k;
  for (k = 0; k < 16; k++)
    res[k] = a[k] | b[k];
  return res;
}

unsigned *
setofb_aun (setofbyte a, setofbyte b)
{
  static setofbyte res;
  int k;
  for (k = 0; k < 16; k++)
    res[k] = a[k] & (~b[k]);
  return res;
}

unsigned *
setofb_its (setofbyte a, setofbyte b)
{
  static setofbyte res;
  int k;
  for (k = 0; k < 16; k++)
    res[k] = a[k] & b[k];
  return res;
}


int 
setofb_eql (setofbyte a, setofbyte b)
{
  int k;
  for (k = 0; k < 16; k++)
    if (a[k] != b[k])
      return FALSE;
  return TRUE;
}

int 
setofb_eq0 (setofbyte a)
{
  int k;
  for (k = 0; k < 16; k++)
    if (a[k] != 0)
      return FALSE;
  return TRUE;
}

void 
setofb_dpl (setofbyte a)
{
  int k;
  for (k = 0; k < 64; k++)
    {
      if (insetb (k, a))
	fprintf (stderr, "%d ", k);
    }
  fprintf (stderr, " ");
  for (; k < 128; k++)
    {
      if (insetb (k, a))
	fprintf (stderr, "%d ", k);
    }
  fprintf (stderr, " ");
  for (; k < 192; k++)
    {
      if (insetb (k, a))
	fprintf (stderr, "%d ", k);
    }
  fprintf (stderr, " ");
  for (; k < 256; k++)
    {
      if (insetb (k, a))
	fprintf (stderr, "%d ", k);
    }
  fprintf (stderr, " ");
}
