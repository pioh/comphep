/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __BEAMS_
#define __BEAMS_

#include"model.h"

extern int In_composite_particle_base (shortstr name);

extern int enter_hadron (int *y, char *name, int num);

extern int copy_hadrons(hadron* dist,hadron src);

extern int get_hadroncontent(hadron h, shortstr res);

extern int enter_beams(void);

extern void construct_full_sf_name (int i);

#endif
