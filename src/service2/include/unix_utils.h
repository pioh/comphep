/*
* Copyright (C) 1997-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __OS_
#define __OS_
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define anyfile       0x00
#define readonly      0x01
#define hidden        0x02
#define sysfile       0x04
#define volumeid      0x08
#define directory     0x10
#define archive       0x20

extern char f_slash;
extern char d_slash;
extern char *defaultPath;
typedef struct
  {
    char fill[21];		/* Reserved by TOS */
    unsigned char attr;		/* Attribute found */
    unsigned time;		/* File time */
    unsigned date;		/* File date */
    long size;			/* File size */
    char name[50];		/* File name found */
#ifdef _WIN32			/* Windows95 */
    HANDLE h;
    WIN32_FIND_DATA w32FD;
#endif
  }
searchrec;

extern int find_first (char *filename, searchrec * filerec, int attrib);
/* It sets the first file at mask filename with attrib. Result in fileRec */
extern int find_next (searchrec * filerec);
extern int find_close (searchrec * filerec);
/* It works with findfirst for the loop in appropriate files */
extern int unlink (const char *path);
/* It erases the filename file, it is used by "erase" */
extern int chepmkdir (char *path);
extern void init_os(void);

#endif
