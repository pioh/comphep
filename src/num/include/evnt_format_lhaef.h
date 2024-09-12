/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __EVENT_FORMAT_LHAEF__
#define __EVENT_FORMAT_LHAEF__

extern int prepare_evfile_lhaef (
 vegasGrid * vegPtr,
 double (*func) (double *, double), 
 char * fname,
 float * cubemaxval, 
 int n_event, 
 int n_cube, 
 double max);

extern int complete_evfile_lhaef (char * fname, int store, int n_event, double mult, double rmax);

extern void set_weighted_flag (int i);

#endif
