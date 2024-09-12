/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __COMBINE__
#define __COMBINE__

extern int combine (int nfiletot, FILE * menuq, FILE ** archives, FILE ** diaginfo, FILE ** catalogs);
extern int compare_menuq (int nfiletot, int nsub, FILE ** menuqs);
#endif
