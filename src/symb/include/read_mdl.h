/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __READ_MDL_
#define __READ_MDL_

#include "chep_crt/include/file_scr.h"
extern table modelTab[5];

#define vars_tab modelTab[0]
#define func_tab modelTab[1]
#define prtcls_tab modelTab[2]
#define lgrng_tab modelTab[3]
#define cpart_tab modelTab[4]

extern table hadron_tab;
extern table strfun_tab;

extern void read_beams(void);
extern int readModelFiles (int l, char * path);
extern int copyModelFiles (int l, char * path);

extern int loadModel (int check);
extern int readhadrons(int check);
extern int readstrfuns(int check);

#endif
