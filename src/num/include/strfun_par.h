/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __STRFUN_PAR__
#define __STRFUN_PAR__

extern void set_alphaMode (int mode);
extern int get_alphaMode (void);
extern void set_sf_mass (int i, double mass);
extern double get_sf_mass (int i);
extern void set_sf_be (int i, double be);
extern double get_sf_be (int i);
extern void set_sf_num (int i, int num);
extern int get_sf_num (int i);

#endif
