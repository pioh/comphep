/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __ALPHAS_MENU__
#define __ALPHAS_MENU__

extern int WriteQCDInfo (FILE * mode);
extern int ReadQCDInfo (FILE * mode);
extern void InitQCDInfo (void);

extern int qcdmen_ (void);
extern double qcd_Scale_chep (void);
extern char * get_scale_form (void);

extern void alf_ (double q);

#endif
