/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __UNWEIGHTING_ROUTINES__
#define __UNWEIGHTING_ROUTINES__

extern int
unweight_events (char ini_name[], char out_name[], int nevnt, long the_seed);
extern int
unweight_events_with_nn_weights (char ini_name[], char out_name[], int nevnt, long the_seed);
extern int
calculate_nn_weights_and_unweight (char ini_name[], char out_name[], int nevnt, long the_seed);

#endif
