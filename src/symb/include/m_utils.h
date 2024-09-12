/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __M_UTILS_
#define __M_UTILS_

extern int clearresults (void);
extern void fillModelMenu (void);
extern int deletemodel (int n);
extern int makenewmodel (void);
extern int continuetest (void);
extern void changeexitmenu (void);
extern void new_user (void);
extern void clear_tmp (void);

extern int writeModelFiles (int l, char * path);
extern int writeHadrons(char * path);
extern int writeStrFuns(char * path);

extern int editBeams (void);
extern int editStrFuns (void);

#endif
