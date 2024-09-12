/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __SOS_
#define __SOS_

extern void save_sos (int ercode);
extern void restoreent (int *exitlevel);
extern void saveent (int exitlevel);

extern int restoreent_dump (void);
extern int batch_composite (char *name, shortstr h);

#endif
