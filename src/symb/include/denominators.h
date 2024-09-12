/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef  __DENOMINATORS__
#define  __DENOMINATORS__
#include"physics.h"

typedef struct
  {
    int power;
    short mass, width;
    char momStr[MAXINOUT];
  }
denom_struct;


typedef struct denlistrec
  {
    struct denlistrec *next;
    int order_num;
    int mass, width;
    char momStr[MAXINOUT];
    double val0, val1, val2;	/*for  Interpreter denominators in power 0,1 
				   and 2 correspondingly      */
  }
denlistrec;
typedef struct denlistrec *denlist;


typedef struct deninforec
  {
    long cr_pos;
    int tot_den;
    struct denarr
      {
	int width;
	int power;
	unsigned order_num;
      }
    denarr[2 * (MAXINOUT - 3)];
  }
deninforec;

extern void calcdenominators (vcsect vcs);
extern void denominatorStatistic (FILE * fres, int nsub,
	int *n_width, int *n_0_width, denlist * allDenominators, FILE * fd);
extern int ttypepropag (int v, int l);

extern void readDenominators (FILE * fres);

#endif
