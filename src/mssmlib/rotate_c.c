extern double rotate5_ (double *v11, double *v12, double *v13, double *v14, double *v15, 
			double *v22, double *v23,double *v24, double *v25, 
			double *v33, double *v34, double *v35, 
			double *v44, double *v45, 
			double *v55, 
			double *ri, double *rj);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v13 = 1.;
			*v14 = 1.;
			*v15 = 1.;
			*v22 = 1.;
			*v23 = 1.;
			*v24 = 1.;
			*v25 = 1.;
			*v33 = 1.;
			*v34 = 1.;
			*v35 = 1.;
			*v44 = 1.;
			*v45 = 1.;
			*v55 = 1.;

			*ri  = 1.;
			*rj  = 1.;
			  return 1.;
			};
*/

extern double rotate4_ (double *v11, double *v12, double *v13, double *v14, 
			double *v22, double *v23, double *v24, 
			double *v33, double *v34, 
			double *v44, 
			double *ri, double *rj);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v13 = 1.;
			*v14 = 1.;
			*v22 = 1.;
			*v23 = 1.;
			*v24 = 1.;
			*v33 = 1.;
			*v34 = 1.;
			*v44 = 1.;

			*ri  = 1.;
			*rj  = 1.;
			  return 1.;
			};
*/

extern double rotate3_ (double *v11, double *v12, double *v13, 
			double *v22, double *v23, 
			double *v33, 
			double *ri, double *rj);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v13 = 1.;
			*v22 = 1.;
			*v23 = 1.;
			*v33 = 1.;

			*ri  = 1.;
			*rj  = 1.;
			  return 1.;
			};
*/

extern double  rotate2_(double *v11, double *v12, 
			double *v22, 
			double *ri, double *rj);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v22 = 1.;

			*ri  = 1.;
			*rj  = 1.;
			  return 1.;
			};
*/

extern double pmass5_ ( double *v11, double *v12, double *v13, double *v14, double *v15, 
			double *v22, double *v23,double *v24, double *v25,
			double *v33, double *v34, double *v35, 
			double *v44, double *v45, 
			double *v55, 
			double *ri);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v13 = 1.;
			*v14 = 1.;
			*v15 = 1.;
			*v22 = 1.;
			*v23 = 1.;
			*v24 = 1.;
			*v25 = 1.;
			*v33 = 1.;
			*v34 = 1.;
			*v35 = 1.;
			*v44 = 1.;
			*v45 = 1.;
			*v55 = 1.;

			*ri  = 1.;
			  return 1.;
			};
*/

extern double pmass4_ ( double *v11, double *v12, double *v13, double *v14, 
			double *v22, double *v23, double *v24, 
			double *v33, double *v34, 
			double *v44, 
			double *ri);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v13 = 1.;
			*v14 = 1.;
			*v22 = 1.;
			*v23 = 1.;
			*v24 = 1.;
			*v33 = 1.;
			*v34 = 1.;
			*v44 = 1.;

			*ri  = 1.;
			  return 1.;
			};
*/

extern double pmass3_ ( double *v11, double *v12, double *v13, 
			double *v22, double *v23, 
			double *v33, 
			double *ri);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v13 = 1.;
			*v22 = 1.;
			*v23 = 1.;
			*v33 = 1.;

			*ri  = 1.;
			  return 1.;
			};
*/

extern double pmass2_ ( double *v11, double *v12, 
			double *v22, 
			double *ri);
/*			{
                        *v11 = 1.;
			*v12 = 1.;
			*v22 = 1.;

			*ri  = 1.;
			  return 1.;
			};
*/

/*********************************** INTERFACES ****************************************************/
/* 17 parameters */
double rotate5 (double v11, double v12, double v13, double v14, double v15, 
		double v22, double v23,double v24,  double v25, 
		double v33, double v34, double v35, 
		double v44, double v45, 
		double v55, 
		double ri, double rj) {
		  return rotate5_(&v11, &v12, &v13, &v14, &v15, &v22, &v23,&v24, &v25, &v33, &v34, &v35, &v44, &v45, &v55, &ri, &rj);
		}

/* 12 parameters */
double rotate4 (double v11, double v12, double v13, double v14, 
		double v22, double v23, double v24, 
		double v33, double v34, 
		double v44, 
		double ri, double rj) {
		  return rotate4_(&v11, &v12, &v13, &v14, &v22, &v23, &v24, &v33, &v34, &v44, &ri, &rj);
		}

/* 8 parameters */
double rotate3 (double v11, double v12, double v13, 
		double v22, double v23, 
		double v33, 
		double ri, double rj) {
		  return rotate3_(&v11, &v12, &v13, &v22, &v23, &v33, &ri, &rj);
		}

/* 5 parameters */
double rotate2 (double v11, double v12, 
		double v22, 
		double ri, double rj) {
		  return rotate2_(&v11, &v12, &v22, &ri, &rj);
		}

/* 16 parameters */
double pmass5 ( double v11, double v12, double v13, double v14, double v15, 
		double v22, double v23, double v24, double v25, 
		double v33, double v34, double v35, 
		double v44, double v45, 
		double v55, 
		double ri) {
		  return pmass5_ (&v11, &v12, &v13, &v14, &v15, &v22, &v23,&v24, &v25, &v33, &v34, &v35, &v44, &v45, &v55, &ri);
		}

/* 11 parameters */
double pmass4 ( double v11, double v12, double v13, double v14, 
		double v22, double v23, double v24, 
		double v33, double v34, 
		double v44, 
		double ri) {
		 return pmass4_ (&v11, &v12, &v13, &v14, &v22, &v23, &v24, &v33, &v34, &v44, &ri);
		 }

/* 7 parameters */
double pmass3 ( double v11, double v12, double v13, 
		double v22, double v23, 
		double v33, 
		double ri) {
	        return pmass3_ (&v11, &v12, &v13, &v22, &v23, &v33, &ri);
	      }

/* 4 parameters */
double pmass2 ( double v11, double v12, 
		double v22, 
		double ri)
		{
		  return pmass2_(&v11, &v12, &v22, &ri);
		}
