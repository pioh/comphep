/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __RANDOM__
#define __RANDOM__
extern double drandXX (void);
extern char *seedXX (char *init);

extern void get_long_seed (unsigned long * sd1, unsigned long * sd2);
extern void set_long_seed (unsigned long sd1, unsigned long sd2);


#endif
