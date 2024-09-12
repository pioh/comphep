/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __SPINOR_
#define __SPINOR_


extern poly multtwospin (poly t1,
			 poly t2,
			 int forspur);

extern poly calcspur (poly spnr);

extern void multspintens (poly * spn,
			  poly * tns);

extern int spinLength;
#endif
