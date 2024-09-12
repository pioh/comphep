/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __EVENT_FORMAT2__
#define __EVENT_FORMAT2__

#include<stdio.h>
extern int prepare_evfile_frmt2 (
 vegasGrid * vegPtr,
 double (*func) (double *, double), 
 char * fname,
 float * cubemaxval, 
 int n_event, 
 int n_cube, 
 double max);

extern int complete_evfile_frmt2 (char * fname, int store, int n_event, double mult, double rmax);
#endif
