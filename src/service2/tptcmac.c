/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Victor Edneral
* ---------------------------------------------------
*/
#include "chep_limits.h"
#include "tptcmac.h"

char *
copy (char *str, int from, int len)
{
  static char buf[STRSIZ];
  buf[0] = '\0';
  if (from > strlen (str))	/* copy past end gives null string */
    return buf;
  strcpy (buf, str + from - 1);	/* skip over first part of string */
  if (len < STRSIZ)
    buf[len] = '\0';		/* truncate after len characters */
  return buf;
}


/* String/character concatenation function
 *
 * This function takes a sprintf-like control string, a variable number of
 * parameters, and returns a pointer a static location where the processed
 * string is to be stored. */
static int scatBufNum = 0;

char *
scat (char *control,...)
{
  va_list args;
  static char buf[4][STRSIZ];

  va_start (args, control);	/* get variable arg pointer */
  scatBufNum++;
  if (scatBufNum > 3)
    scatBufNum = 0;
  vsprintf (buf[scatBufNum], control, args);	/* format into buf with variable args */
  va_end (args);		/* finish the arglist */

  return buf[scatBufNum];	/* return a pointer to the string */
}

/* string build - like scat, sprintf, but will not over-write any input parameters */
void 
sbld (char *dest, char *control,...)
{
  va_list args;
  char buf[STRSIZ];
  va_start (args, control);	/* get variable arg pointer */
  vsprintf (buf, control, args);	/* format into buf with variable args */
  va_end (args);		/* finish the arglist */
  strcpy (dest, buf);		/* copy result */
}


/* spos(str1,str2) - returns index of first occurence of str1 within str2;
 *    1=first char of str2
 *    0=nomatch */
int 
spos (char *str1, char *str2)
{
  int i, k, n;
  n = (int) strlen (str2) - (int) strlen (str1);
  for (k = 0; k <= n; k++)
    {
      i = 0;
      while (str1[i] != '\0' && str1[i] == str2[i + k])
	i++;
      if (str1[i] == '\0')
	return k + 1;
    }
  return 0;
}


/* cpos(str1,str2) - returns index of first occurence of c within str2;
 *    1=first char of str2
 *    0=nomatch */
int 
cpos (char c, char *str2)
{
  char *res;
  res = strchr (str2, c);
  if (res == NULL)
    return 0;
  else
    return (int) (res - str2 + 1);
}
