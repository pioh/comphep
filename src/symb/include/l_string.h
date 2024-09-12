/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __L_STRING_
#define __L_STRING_

extern void fortwriter (char *name, varptr fortformula);
extern void initfortwriting (char type);
extern int write_const (void);
extern int cleardegnames (void);
extern void initdegnames (void);
int gettmpNameMax (void);
void settmpNameMax (int num);

#endif
