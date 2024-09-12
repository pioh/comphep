/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __RUNVEGAS__
#define __RUNVEGAS__

extern int runVegas (int init);

extern int ClearVegasGrid (void);
extern int WriteVegasGrid (FILE * f);
extern int ReadVegasGrid (FILE * f);

#endif
