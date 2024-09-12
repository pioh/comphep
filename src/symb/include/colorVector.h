/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __COLOR_VECTOR__
#define __COLOR_VECTOR__

#include "amplitudes.h"

typedef struct
  {
    double n, d;
  }
rat;

extern void generateColorVectors (vamplExt * ans, int *nvect, rat ** c_vectors);

#endif
