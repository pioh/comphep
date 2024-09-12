/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
*------------------------------------------------------
*/
#include <stdlib.h>
#include <stdio.h>
#include "drandXX.h"

static double float48 = 1. / (((double) 0x10000) * ((double) 0x10000) * ((double) 0x10000));

static unsigned long Xlong = 0x1234ABCD;
static unsigned long Xshort = 0x330E;
#define Along   0x5DEEC
#define Ashort  0xE66D
#define Along16 0xDEEC0000
#define  Clong  0
#define  Cshort 0xB
#define FIRST16 0xFFFF


double 
drandXX (void)
{
  unsigned long bot = Xshort * Ashort;
  unsigned long top = Clong + (bot >> 16);

  bot = (bot & FIRST16) + Cshort;

  Xlong = top + (bot >> 16) + Along * Xshort + Ashort * Xlong + ((Along16 * Xlong));
  Xshort = bot & FIRST16;

  return ((double) Xshort + (Xlong & 0xFFFFFFFF) * (double) 0x10000) * float48;
}

char *
seedXX (char *init)
{
  unsigned long Xlong_, Xshort_;
  static char cbuff[128];
  static char frmt1[128];
  static char frmt2[128];
  int lenth = 2*sizeof(unsigned long);

  sprintf (frmt1, "%s%i%s", "%0", lenth, "lX%04lX");
  sprintf (frmt2, "%s%i%s", "%", lenth, "lX%4lX");

  sprintf (cbuff, frmt1, Xlong, Xshort);

  if (init)
    {
      if (sscanf (init, frmt2, &Xlong_, &Xshort_) == 2)
	{
	  Xlong = Xlong_;
	  Xshort = Xshort_;
	}
      else
	return NULL;
    }
  return cbuff;
}
void get_long_seed (unsigned long * seed1, unsigned long * seed2)
{
  *seed1 = Xlong;
  *seed2 = Xshort;
}

void set_long_seed (unsigned long seed1, unsigned long seed2)
{
  Xlong = seed1;
  Xshort = seed2;
}
