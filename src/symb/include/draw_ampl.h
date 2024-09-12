/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __DRAW_AMPL__
#define __DRAW_AMPL__

#include"diagrams.h"

extern int amplitudeFrameX (int nout);
extern int amplitudeFrameY (int nout);
extern void drawAmpitudeDiagram (vampl * diagr, int x, int y);
void drawSquaredDiagram (csdiagram * diagr, int x_init, int x_final, int y);
extern void setPictureScale (int upr_, int *xn, int *ynu);

#endif
