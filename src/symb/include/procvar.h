/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PROCVAR_
#define __PROCVAR_

#include "physics.h"

typedef struct
  {
    char alias[10];
    double tmpvalue;
    int used;
  }
singlevardescription;

extern singlevardescription *vararr;
extern int nProcessVar;

extern int initvararray (FILE * fres, int nsub, char key);


#endif
