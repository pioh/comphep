/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __SYST_
#define __SYST_

extern void (*diskerror) (void);
extern int f_printf (FILE * fp, char *format,...);
extern size_t f_write (void *ptr, size_t size, size_t n, FILE * fp);
extern char *trim (char *);

extern void revers (void **list);
extern void lShift (char *s, int l);

extern long get_seed (char filename[]);
extern int get_sf_info (char sf_info[], char * name, int * set, int * mem);

#ifdef LHAPDF
extern int update_lhapdf_mdl (void);
#endif

#endif
