/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/

#include "strfun_par.h"

static int sf_num[2]     = {0, 0};
static double sf_mass[2] = {0, 0};
static double sf_be[2]   = {0, 0};
static int alphaMode = 0;


void set_alphaMode (int mode)
{
  alphaMode = mode;
}

int get_alphaMode (void)
{
  return alphaMode;
}

void set_sf_mass (int i, double mass)
{
  sf_mass[i] = mass;
}

double get_sf_mass (int i)
{
  return sf_mass[i];
}

void set_sf_be (int i, double be)
{
  sf_be[i] = be;
}

double get_sf_be (int i)
{
  return sf_be[i];
}

void set_sf_num (int i, int num)
{
  sf_num[i] = num;
}

int get_sf_num (int i)
{
  return sf_num[i];
}
