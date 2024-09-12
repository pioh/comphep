/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PLOT__
#define __PLOT__

extern void plot_histo (double xMin, double xMax,
		    int dim, double *f, double *ff,
		    char *upstr, char *xstr, char *ystr);

extern void plot_table (double xMin, double xMax,
		    int dim, double *f, double *ff,
		    char *upstr, char *xstr, char *ystr);
#endif
