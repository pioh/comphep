/*
* Copyright (C) 2001-2008, CompHEP Collaboration
* ------------------------------------------------------
* $Id$
*
* $Log$
*/
#ifndef __REGFUNALL_
#define __REGFUNALL_

extern void regfun_ (int itype, int nsing, sing_struct * singar, double xmin, double xmax,
                     double xx, double *xout, double *factor);
extern void regfct_ (int itype, int nsing, sing_struct * singar, double xmin, double xmax,
                     double xout, double *factor);

#endif
