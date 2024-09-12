/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __OUT_SERVICE
#define __OUT_SERVICE

#include "service2/include/files.h"
#include "physics.h"
#include "polynom/include/polynom.h"

extern int outputLanguage;
extern void writeF (char *format,...);
extern void outFileOpen (char *fname);
extern void outFileClose (void);

extern void momentToString (char *moment, char *outstr);

extern void rewritepolynom (FILE * fsrc);
extern void findPrtclNum (char *procName, int *prtclNum);
extern void emitconvlow (int *prtclNum);
extern void writeLabel (char comment);
extern void DiagramToOutFile (vcsect * vcs, int label, char comment);

extern void makeOutput (void (*startOutput) (FILE * , int, int *, int),
			void (*diagramOutput) (FILE * , vcsect *, catrec *),
			void (*endOutput) (int *)
);

extern void readvardef (FILE * fres);
extern void clearvardef (void);

extern FILE *outFile;
extern poly readBuff;
extern int readSize;

#define wrtBuffSize 2000
#endif
