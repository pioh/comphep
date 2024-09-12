/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "service2/include/drandXX.h"
#include "service2/include/4_vector.h"

#include "const.h"
#include "evnt_tools.h"

char new_symbol (FILE * file) {
  static int char_switch = 0;
  char sym = fgetc (file);

  if (39 == sym) {
    char_switch++;
    char_switch = char_switch % 2;
  }

  if (char_switch) {
    return sym;
  }
  while (32 == sym || 10 == sym) {
    sym = fgetc (file);
  }
  return sym;
}


int CheckFormat (FILE * f) {
  int NF = -1;
  char ch;

  rewind (f);
  ch = new_symbol (f);
  if (ch == '#') {
    ch = new_symbol (f);
    if (ch == '#') {
      NF = 2;     /* CompHEP cpyth2 event file */
    } else if (ch == 'C') {
      ch = new_symbol (f);
      if (ch == 'a')
        NF = 5;     /* CalcHEP event file */
      else
        NF = 1;     /* CompHEP cpyth1 event file */
    } else if (ch == 'P') {
      NF = 3;
    }
  } else {
    if (ch == '<') {
      char word[128];
      fscanf (f, "%s", word);
      if (0 == strcmp (word, "LesHouchesEvents")) {
        NF = 4;     /* CompHEP lhe event file */
      }
    }
  }
  rewind (f);

  return NF;
}


void rnd_rotate_momentum (int nnin, int nnout) {
  int i;
  int nntot = nnin + nnout;
  double phi_rnd = 2 * M_PI * drandXX ();
  double cs = cos (phi_rnd);
  double sn = sin (phi_rnd);
  for (i = nnin; i < nntot; ++i) {
    double XX = pvect[4 * i + 1];
    double YY = pvect[4 * i + 2];
    pvect[4 * i + 1] = XX * cs - YY * sn;
    pvect[4 * i + 2] = XX * sn + YY * cs;
  }
}
