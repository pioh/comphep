/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"

#include "physics.h"
#include "chess.h"

vertinfostr vertinfo[2 * maxvert] = {{0,0,0,{0},{0}}};

int n_vrt = 0;
int prgcode[2 * maxvert][2] = {{0}};


static int ncode;


static int 
setpower (indvertset ss)
{
  int sp = 0, i = 1;
  indvertset s;

  setofb_cpy (s, ss);
  while (!setofb_eq0 (s))
    {
      if (insetb (i, s))
	{
	  ++(sp);
	  setofb_cpy (s, setofb_aun (s, setofb (i, _E)));
	}
      ++(i);
    }
  return sp;
}


static int 
firstelement (indvertset s)
{
  int fe;

  if (setofb_eq0 (s))
    return 0;
  fe = 1;
  while (!insetb (fe, s))
    ++(fe);
  return fe;
}


static void 
split (indvertset s, indvertset s1, indvertset s2)
{
  int nsub = 0;
  indvertset ss, s1_, s2_;
  int w = 0, w_;
  int i, j, l = 0, cross, cross_, dim;
  int *nextelem;

  dim = 1;
  for (i = 1; i <= 2 * MAXINOUT - 5; i++)
    dim = 2 * dim;
  nextelem = m_alloc (dim * sizeof (int));
  setofb_cpy (ss, s);
  while (!setofb_eq0 (ss))
    {
      if (insetb (l, ss))
	{
	  w += vertinfo[l - 1].weit;
	  setofb_cpy (ss, setofb_aun (ss, setofb (l, _E)));
	  if (!setofb_eq0 (ss))
	    {
	      ++(nsub);
	      nextelem[nsub - 1] = l;
	      for (i = 1; i <= nsub - 1; i++)
		nextelem[nsub + i - 1] = -nextelem[nsub - i - 1];
	      nsub = 2 * nsub - 1;
	    }
	}
      ++(l);
    }
  l = nextelem[0];
  w -= 2 * vertinfo[l - 1].weit;
  setofb_cpy (s2, setofb (l, _E));
  setofb_cpy (s1, setofb_aun (s, s2));
  cross = 0;
  for (j = 1; j <= vertinfo[l - 1].vlnc; j++)
    if (insetb (vertinfo[l - 1].link[j - 1], s1))
      ++(cross);
  setofb_cpy (s1_, s1);
  setofb_cpy (s2_, s2);
  w_ = w;
  cross_ = cross;
  for (i = 2; i <= nsub; i++)
    {
      l = abs (nextelem[i - 1]);
      if (nextelem[i - 1] > 0)
	{
	  setofb_cpy (s1_, setofb_aun (s1_, setofb (l, _E)));
	  for (j = 1; j <= vertinfo[l - 1].vlnc; j++)
	    if (insetb (vertinfo[l - 1].link[j - 1], s1_))
	      ++(cross_);
	    else if (insetb (vertinfo[l - 1].link[j - 1], s2_))
	      --(cross_);
	  setofb_cpy (s2_, setofb_uni (s2_, setofb (l, _E)));
	  w_ -= 2 * vertinfo[l - 1].weit;
	}
      else
	{
	  setofb_cpy (s2_, setofb_aun (s2_, setofb (l, _E)));
	  for (j = 1; j <= vertinfo[l - 1].vlnc; j++)
	    if (insetb (vertinfo[l - 1].link[j - 1], s2_))
	      ++(cross_);
	    else if (insetb (vertinfo[l - 1].link[j - 1], s1_))
	      --(cross_);
	  setofb_cpy (s1_, setofb_uni (s1_, setofb (l, _E)));
	  w_ += 2 * vertinfo[l - 1].weit;
	}

      if (MEMORY_OPTIM)
	{
	  if (abs (w_) < abs (w) ||
	      (abs (w_) == abs (w) && cross_ < cross))
	    {
	      setofb_cpy (s1, s1_);
	      setofb_cpy (s2, s2_);
	      w = w_;
	      cross = cross_;
	    }
	}
      else
	{
	  if (cross_ < cross ||
	      (cross_ == cross && abs (w_) < abs (w)))
	    {
	      setofb_cpy (s1, s1_);
	      setofb_cpy (s2, s2_);
	      w = w_;
	      cross = cross_;
	    }
	}
    }
  free (nextelem);
}


static void 
programer (indvertset s)
{
  indvertset s1, s2;

  if (setpower (s) > 1)
    {
      split (s, s1, s2);
      ++(ncode);
      prgcode[ncode - 1][0] = firstelement (s1);
      prgcode[ncode - 1][1] = firstelement (s2);
      programer (s1);
      programer (s2);
    }
}

void 
makeprgcode (void)
{
  indvertset ss;
  setofb_cpy (ss, setofb (1, UpTo, n_vrt, _E));
  ncode = 0;
  programer (ss);
}
