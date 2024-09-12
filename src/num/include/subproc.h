/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __SUBPROC__
#define __SUBPROC__

typedef struct
  {
    int nsub;
    shortstr proces;
  }
proces_;
extern proces_ proces_1;

char * get_subproc_name (void);
char * get_final_state_name (void);

#endif
