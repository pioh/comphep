/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __HIST__
#define __HIST__

#include"chep_crt/include/file_scr.h"

extern table histTab;

extern void manipulateHists (void);
extern int correctHistList (int mode);

extern int WriteHistograms (FILE * nchan);
extern int ReadHistograms (FILE * nchan);

extern void showHistAll (void);
extern void clearHists (void);
extern void fillHists (double w);
#endif
