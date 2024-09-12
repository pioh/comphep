/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __EVENTREADER__
#define __EVENTREADER__

extern int getLHAevent (char fname[], FILE * source, long startpos, eventUP * evt);
extern int setLHAevent (char fname[], FILE * s, long startpos, eventUP * evt);
extern int setLHAevent_with_comments (char fname[], FILE * s, long startpos, eventUP * evt, char comments[]);
extern int moveLHAevent (int num, char fname[], FILE * source, long spos, FILE * target, long tpos, int zrandom);
extern int testLHAevent (FILE * f, char fname[], long pos);
extern char * get_event_comments (void);

#endif
