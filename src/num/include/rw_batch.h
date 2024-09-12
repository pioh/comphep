/*
* Copyright (C) 2008-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __RW_BATCH__
#define __RW_BATCH__
#include<stdio.h>

extern int w_batch_session (FILE * mode);
extern int r_batch_session (FILE * mode);
extern int init_batch_session (void);
extern int get_process_name_safe (char * procname);
#endif
