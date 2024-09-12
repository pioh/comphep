/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PRE_READ_
#define __PRE_READ_
#include "service2/include/sets.h"

typedef setrec
{
  abracadabra
}
typeindlist;			/* [1..9], used by Editor  */

typedef struct preresrecord
  {
    struct preresrecord *next;
    int free;
    int tp;
    long num;
    unsigned maxp;		/*  How many pulses: used by Editor   */
    unsigned degp;
    int g5;
    unsigned maxg;
    typeindlist indlist;
    int nvar;
    unsigned *varsdeg;
  }
preresrecord;
typedef struct preresrecord *preres;

extern preres pregarbage;

extern void clearpregarbage (void);
extern void *act_pre (char *ch, int n, void **agrs);
extern void *act_preF (char *ch, int n, void **agrs);
extern void *rd_pre (char *s);

#endif
