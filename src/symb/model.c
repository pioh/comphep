/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include<stdlib.h>
#include<string.h>

#include "service2/include/chep_limits.h"

#include"model.h"

int nmodelvar = 0;
varlist modelvars = NULL;
prtcl_base *prtclbase = NULL;

int nparticles = 0;                     /*  Number particles in model */
algvertptr lgrgn = NULL;

int n_cpart = 0;
cparticles cpartbase = NULL;

int n_hadron = 0;
_hadron_  hadronbase = NULL;
hadron hadrons[MAXINOUT] = {{"",0.,0,{0},""}};
__beam__ beam[2]= {{{"",0.,0,{0},""}, 0.}};

int n_strfun = 0;
_strfun_ strfunbase = NULL;


int 
locateinbase (char *name)
{
  int i;
  for (i = 0; i < nparticles; i++)
    {
      if (strcmp (name, prtclbase[i].name) == 0)
	return i + 1;
    }
  return 0;
}


int
prtclname (int n, vshortstr nm)
{
  if (n > nparticles-1 || n < 1) strcpy(nm,"");
  else strcpy(nm,prtclbase[n-1].name);
  return 0;
}



int 
pseudop (int p)
{
  return (p == 0 || p > nparticles || prtclbase[p - 1].hlp == '*');
}
int 


fermionp (int p)
{
  return (prtclbase[p - 1].spin % 2 == 1 && p <= prtclbase[p - 1].anti);
}


int 
a_fermionp (int p)
{
  return (prtclbase[p - 1].spin % 2 == 1 && p >= prtclbase[p - 1].anti);
}


int 
bosonp (int p)
{
  return (prtclbase[p - 1].spin % 2 == 0);
}


int 
vectorp (int p)
{
  return (prtclbase[p - 1].spin == 2);
}


int 
zeromass (int p)
{
  return (strcmp (prtclbase[p - 1].massidnt, "0") == 0);
}
int 


photonp (int p)
{
  return ((vectorp (p) && zeromass (p)) || pseudop (p));
}


int 
ghostp (int p)
{
  return (prtclbase[p - 1].hlp == 'c' || prtclbase[p - 1].hlp == 'C' || prtclbase[p - 1].hlp == 'f');
}


int 
gaugep (int p)
{
  return (p != 0 && prtclbase[p - 1].hlp == 'G');
}


int 
ghostmother (int j)
{
  if (j == 0)
    return 0;
  switch (prtclbase[j - 1].hlp)
    {
    case 'c':
      return j - 1;
    case 'C':
      return j - 2;
    case 'f':
      return j - 3;
    case 't':
      return j + 1;
    default:
      return j;
    }
}
