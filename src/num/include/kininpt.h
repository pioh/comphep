/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __KININPT__
#define __KININPT__
#include<stdio.h>

extern int EnterKinScheme (void);

extern void InitKinScheme (void);
extern int WriteKinScheme (FILE * nchan);
extern int ReadKinScheme (FILE * nchan);
#endif
