/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __PHYS_VAL__
#define __PHYS_VAL__

/*extern double calcPhysVal (char key, char *lv);*/
extern double calcPhysVal (char key, char *lv, char *restframe);
extern int checkPhysVal (char *name, char *key, char *plist);
extern void xName (char key, char *plist, char *xname, char *units);

#endif
