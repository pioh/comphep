/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __CORE_DATA__
#define __CORE_DATA__

extern long get_ncalls (void);
extern int get_iters (void);
extern void set_ncalls (long ncalls);
extern void set_iters (int iteras);

extern int get_nsession (void);
extern void set_nsession (int ns);

extern double get_rapidity (void);
extern void set_rapidity (double ra);


typedef struct vegas_integral
  {
    double s0, s1, s2;
    long nCallTot;
    int n_it;
    int old;
  }
vegas_integral;

typedef struct mcintr
  {
    long ncall0;
    int itmx0;
    int wrtEvnt;
  }
mcintr;

extern vegas_integral get_vegas_integral (void);
extern void set_vegas_integral (vegas_integral i);
extern void init_vegas_integral (void);

extern mcintr get_mc_info (void);
extern void set_mc_info (mcintr mc);
#endif
