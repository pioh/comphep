/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Author: A.Kryukov 
* ------------------------------------------------------
*/
#ifndef __CWEIGHT_
#define __CWEIGHT_

extern void cwtarg (vcsect * g);
extern void c_basis_coef (vampl * g, int pow, int nc, int *chains, long *num, long *den);

extern int generateColorWeights (csdiagram * csdiagr,
				 int cBasisPower, int nC, int *cChains,
				 long *cCoefN, long *cCoefD);

extern int infCbases (int np,	/* number of particles */
		      int *cweight,	/* array of particle color weights */
		      int *nc,	/* number of color chains */
		      int *pow,	/* power of basis */
		      int **chains	/* returns array   which descibes   
					   (*pow) basis elements, 
					   each of them contains  (*nc)  chains,
					   each of them is a couple of 
					   particle numbers 
					 */
);

#endif
