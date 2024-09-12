/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __FILL__
#define __FILL__
#include "LesHouches.h"
#include "structures.h"
#include "tag_reader.h"

extern int validity (tags ** the_tags);
extern int fill_LH_structures (int nf, tags ** tbase, process_ * prUP, const char names[], int lenth);

int fill_heprup (processUP * prUP, tags * head);

int prefill_hepeup (eventUP * ev, tags * head);
int fill_hepeup (eventUP * ev, char *evstr, int pnum);

#endif
