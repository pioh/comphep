/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PROCESS_
#define __PROCESS_

#define ycons 19

extern int enter_process (int ProcessChoice);
extern int enter_beams (void);

extern char ** stritems (char *format, char *s);
extern int input (int y0, char *directive, char *text, int cur, int lim);
extern char ** stritems (char *format, char *s);
void prtcllist (int key, int beam);

#endif
