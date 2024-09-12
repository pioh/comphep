/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* author: Alexander Sherstnev
* ------------------------------------------------------
*/
#include "core_data.h"

static int n_session = 1;
static double rapidity = 0.;
static vegas_integral integral = { 0, 0, 0, 0, 0, 0 };
static mcintr mcintr_1 = { 60000, 5, 0 };

mcintr get_mc_info (void) {
  return mcintr_1;
}

void set_mc_info (mcintr mc) {
  mcintr_1.itmx0   = mc.itmx0  ;
  mcintr_1.ncall0  = mc.ncall0 ;
  mcintr_1.wrtEvnt = mc.wrtEvnt;
}

long get_ncalls (void) {
  return mcintr_1.ncall0;
}

int get_iters (void) {
  return mcintr_1.itmx0;
}

int get_nsession (void) {
  return n_session;
}

void set_ncalls (long ncalls) {
  mcintr_1.ncall0 = ncalls;
}

void set_iters (int iteras) {
  mcintr_1.itmx0 = iteras;
}

void set_nsession (int ns) {
  n_session = ns;
}

double get_rapidity(void) {
  return rapidity;
}

void set_rapidity(double ra) {
  rapidity = ra;
}

vegas_integral get_vegas_integral (void) {
  return integral;
}

void set_vegas_integral (vegas_integral in) {
    integral.s0 = in.s0;
    integral.s1 = in.s1;
    integral.s2 = in.s2;
    integral.nCallTot = in.nCallTot;
    integral.n_it     = in.n_it    ;
    integral.old      = in.old     ;
}

void init_vegas_integral (void) {
    integral.s0 = 0.0;
    integral.s1 = 0.0;
    integral.s2 = 0.0;
    integral.nCallTot = 0;
    integral.n_it     = 0;
    integral.old      = 0;
}

