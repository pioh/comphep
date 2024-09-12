/*
* Copyright (C) 2009-2010, CompHEP Collaboration
*------------------------------------------------------
*/
#ifndef __ROTATE__
#define __ROTATE__

/* 17 parameters */
extern 
double rotate5 (double v11, double v12, double v13, double v14, double v15, 
		double v22, double v23,double v24,  double v25, 
		double v33, double v34, double v35, 
		double v44, double v45, 
		double v55, 
		double ri, double rj);

/* 12 parameters */
extern 
double rotate4 (double v11, double v12, double v13, double v14, 
		double v22, double v23, double v24, 
		double v33, double v34, 
		double v44, 
		double ri, double rj);

/* 8 parameters */
extern 
double rotate3 (double v11, double v12, double v13, 
		double v22, double v23, 
		double v33, 
		double ri, double rj);

/* 5 parameters */
extern 
double rotate2 (double v11, double v12, 
		double v22, 
		double ri, double rj);

/* 16 parameters */
extern 
double pmass5 ( double v11, double v12, double v13, double v14, double v15, 
		double v22, double v23, double v24, double v25, 
		double v33, double v34, double v35, 
		double v44, double v45, 
		double v55, 
		double ri);

/* 11 parameters */
extern 
double pmass4 ( double v11, double v12, double v13, double v14, 
		double v22, double v23, double v24, 
		double v33, double v34, 
		double v44, 
		double ri);

/* 7 parameters */
extern 
double pmass3 ( double v11, double v12, double v13, 
		double v22, double v23, 
		double v33, 
		double ri);

/* 4 parameters */
extern 
double pmass2 ( double v11, double v12, 
		double v22, 
		double ri);

#endif
