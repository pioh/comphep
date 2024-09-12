/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __RW_SESS__
#define __RW_SESS__
#include<stdio.h>

extern int write_prt (FILE * mode);

extern int write_session (void);
extern int read_session (void);
extern int init_session (void);

extern void clearSession (void);	/* for interpreter */
extern int ComposeSubprocessString (void);
extern int getModelNumber (void);
#endif
