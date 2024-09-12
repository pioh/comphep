/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <string.h>

#include "service2/include/chep_limits.h"
#include "subproc.h"

proces_ proces_1 =
{1, ""};

char * get_subproc_name (void) {
  return proces_1.proces;
}

char * get_final_state_name (void) {
  return strstr (proces_1.proces, "->") + 3;
}
