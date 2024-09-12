/*
 * Copyright (C) 2007, Alexander Sherstnev 
 * Copyright (C) 2001-2007, CompHEP Collaboration
 *------------------------------------------------------
 * $Id$
 *
 * $Log$
*/
#ifndef __SF_Clhapdf_
#define __SF_Clhapdf_

#define f2cFortran

#include "cfortran.h"

PROTOCCALLSFSUB1(initpdfset,initpdfset,STRINGV)
#define initpdfset(name) \
CCALLSFSUB1(initpdfset,initpdfset,STRINGV,name) 

PROTOCCALLSFSUB1(initpdf,initpdf,INT)
#define initpdf(mem) \
CCALLSFSUB1(initpdf,initpdf,INT,mem) 

PROTOCCALLSFSUB2(getlam5,getlam5,INT,DOUBLE)
#define getlam5(mem,qcdl5) \
CCALLSFSUB2(getlam5,getlam5,INT,DOUBLE,mem,qcdl5)

PROTOCCALLSFSUB3(evolvepdf,evolvepdf,DOUBLE,DOUBLE,DOUBLEV)
#define evolvepdf(x,Q,f) \
CCALLSFSUB3(evolvepdf,evolvepdf,DOUBLE,DOUBLE,DOUBLEV,x,Q,f) 

PROTOCCALLSFSUB3(locevolvepdf,locevolvepdf,DOUBLE,DOUBLE,DOUBLEV)
#define locevolvepdf(x,Q,pdf) \
CCALLSFSUB3(locevolvepdf,locevolvepdf,DOUBLE,DOUBLE,DOUBLEV,x,Q,pdf) 
/*
PROTOCCALLSFSUB1(getorderas,getorderas,INT)
#define getorderas(mem) \
CCALLSFSUB1(getorderas,getorderas,INT,order) 
*/
/* PROTOCCALLSFSUBn is optional for C, but mandatory for C++.

PROTOCCALLSFSUBn(ROUTINE_NAME,routine_name,argtype_1,...,argtype_n)
#define     Routine_name(argname_1,..,argname_n)               \
CCALLSFSUBn(ROUTINE_NAME,routine_name,argtype_1,...,argtype_n, argname_1,..,argname_n) 

PROTOCCALLSFFUNn(routine_type,ROUTINE_NAME,routine_name,argtype_1,...,argtype_n)
#define     Routine_name(argname_1,..,argname_n)               \
CCALLSFFUNn(ROUTINE_NAME,routine_name,argtype_1,...,argtype_n, argname_1,..,argname_n) 

*/

PROTOCCALLSFFUN1(DOUBLE,LHAPDFALPHAS,lhapdfalphas,DOUBLE)
#define lhapdfalphas(Q) \
CCALLSFFUN1(LHAPDFALPHAS,lhapdfalphas,DOUBLE,Q)

PROTOCCALLSFFUN2(DOUBLE,LHAPDFQCDLAM,lhapdfqcdlam,INT,INT)
#define lhapdfqcdlam(set,mem) \
CCALLSFFUN2(LHAPDFQCDLAM,lhapdfqcdlam,INT,INT,set,mem)

PROTOCCALLSFFUN0(INT,LHAPDFQCDORDER,lhapdfqcdorder)
#define lhapdfqcdorder() \
CCALLSFFUN0(LHAPDFQCDORDER,lhapdfqcdorder)

#endif
