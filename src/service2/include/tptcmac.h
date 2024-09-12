/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __TPTCMAC_
#define __TPTCMAC_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <stdarg.h>
typedef void *pointer;

typedef unsigned char byte;



#define lvcpy(d,s)      memcpy((char *)d,(char *)s,(size_t)sizeof(d))
/*It copies SIZEOF(dest) bytes of src to the dest address. */

/* Definition of constans */


#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif


/* String's  department */

/* char *copy(char *str,int from,int len):
 * copy len bytes from the dynamic string dstr
 * starting at position from.
 *
 * String/character concatenation function
 * char *scat(char *control, ...):
 * This function takes a sprintf-like control string, a variable number of
 * parameters, and returns a pointer a static location where the processed
 * string is to be stored.
 *
 * void sbld(char *dest,char *control, ...):
 * string build - like scat, sprintf
 , but will not over-write any
 *                input parameters.
 *
 * int spos(char *str1,char *str2):
 * returns index of first occurence of str1 within str2;
 *    1=first char of str2
 *    0=nomatch
 *
 * int cpos(char c,char *str2):
 * returns index of first occurence of c within str2;
 *    1=first char of str2
 *    0=nomatch
 */

extern char *copy (char *str, int from, int len);
extern char *scat (char *control,...);
extern void sbld (char *dest, char *control,...);
extern int spos (char *str1, char *str2);
extern int cpos (char c, char *str2);

#endif
