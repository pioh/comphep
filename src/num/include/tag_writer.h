/*
* Copyright (C) 2002-2008, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef _CWRITER_H_
#define _CWRITER_H_
#include "tag_reader.h"
#include "LesHouches.h"

#define MAXP 10

extern int intlen (int src);

extern void tag_writer (FILE * file, elementary_tag * tag);

extern int cap_writer (FILE * file, tags * p);
extern int change_cap (tags * t, int nEvents, double mult, double rmax, double cs, double er);

extern int final_write_cap (const char *outFileName, long *nEvent, int tot);
extern int write_cap (FILE * file, tags ** the_tags, process_ prUP);
extern int write_cap_new (FILE * file, tags ** the_tags, int nf, process_ prUP);
extern int distilling (const char * inname);

#endif /* reader.h  */
