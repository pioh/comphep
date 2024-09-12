/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __PARAM__
#define __PARAM__

extern void selectParam (int *position, int x, int y, int sqtrS_on, int vars_on,
			 int func_on, char *mess, void **pscrPrt);
extern int change_parameter (int x, int y);
extern void show_depend (int x, int y);
extern int WriteConstraints (FILE * f);

#endif
