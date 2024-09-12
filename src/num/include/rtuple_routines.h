/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __RTUPLE_ROUTINES__
#define __RTUPLE_ROUTINES__

#ifdef __cplusplus
extern "C" {
#endif
int
book_rtuple (char rtuple_name[]);

int
write_rtuple (void);

int
fill_event (eventUP * ev);
#ifdef __cplusplus
}
#endif

#define RTPLMAXINOUT 25

#endif
