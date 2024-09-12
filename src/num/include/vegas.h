/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __VEGAS__
#define __VEGAS__

extern int simplexOn;

typedef struct vegasGrid
  {

    int ndim,                   /* number of dimensions */
      ndmx;
    double *x_grid;
    double *c_grid;
  }
vegasGrid;

extern vegasGrid *vegas_init
  (int dim,                     /* number of dimensions */
   int ndmx                     /* size of grid */
);

extern void vegas_finish (vegasGrid * vegPtr);

extern int vegas_int (vegasGrid * vegPtr,
                      long ncalls,      /* number of integrand calls */
                      int  nsubpr,      /* subprocess number */
                      int  niters,      /* number of iterations */
                      int  curitr,      /* current iteration number */
                      double (*fxn) (double *, double),         /* integrand */
                      double *ti,       /* integral estimation */
                      double *tsi       /* standard deviation */
);

extern int vegas_int_noGUI (vegasGrid * vegPtr,
                      long ncall0,      /* number of integrand calls */
                      double alph,      /* rate of grid improvement  */
                      double (*fxn) (double *, double),         /* integrand */
                      double *ti,       /* integral estimation */
                      double *tsi       /* standard deviation */
);


extern int vegas_max (
                       vegasGrid * vegPtr,
                       int nCubs,
                       int nPoints,
                       double (*fxn) (double *, double),
                       double milk,
                       double *eff,
                       float *fmax
);


extern int vegas_events (
                          vegasGrid * vegPtr,
                          long nCubs,
                          long nEvents,
                          double gmax,
                          double (*fxn) (double *, double),
                          void (*out) (long, int, double),
                          float *fmax
);

extern int vegas_1to2_events (
                          vegasGrid * vegPtr,
                          long nCubs,
                          long nEvents,
                          double gmax,
                          double (*fxn) (double *, double),
                          void (*out) (long, int, double),
                          float *fmax
);

extern int vegas_wgt (
                       vegasGrid * vegPtr, 
		       int nCubsINI, 
		       int nEvents,
                       double gmax,
                       double (*fxn) (double *, double), 
		       void (*out) (long, int, double),
                       float *fmax
);

extern long generateVegasCubs (vegasGrid * vegPtr, long nCubs);

extern double get_efficiency (void);
extern double get_rmax (void);
extern double get_multiplicity (void);
extern double get_negativity (void);

extern void setStopSymb (int s);
extern int getCurCub (void);
extern int 
vegas_int_nongui (vegasGrid * vegPtr, long Ncalls, double alph,
           double (*fxn) (double *, double), double *ti, double *tsi);

#endif
