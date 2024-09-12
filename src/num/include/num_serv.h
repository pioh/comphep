/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __NUM_SERV_
#define __NUM_SERV_
typedef double (*r_func) (void);
extern void paramtable2 (void);
extern void paramdependence (r_func ff, char *procname, char *resultname);
extern void paramdependence1 (r_func ff, char *procname, char *resultname);
#endif
