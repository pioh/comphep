/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __AMPL__
#define __AMPL__

#include"diagrams.h"


typedef struct vamplExt
  {
    struct vamplExt *next;
    int status;
    int nsub;
    vampl Ampl;
  }
vamplExt;


extern void Amplitudes (void);

#endif
