/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structures.h"
#include "fill.h"
#include "tag_writer.h"
#include "tag_parser.h"
#include "compare.h"

static void quicksort (int * array, int llimit, int rlimit)
{
  int temp;
  int left = llimit;
  int right = rlimit;
  int pivot = (left + right) / 2;	// find the median
  int median = array[pivot];

  do
    {
      while ((array[left] < median) && (left < rlimit))
	{
	  left++;
	}
      while ((median < array[right]) && (right > llimit))
	{
	  right--;
	}
      if (left <= right)
	{
	  temp = array[left];
	  array[left] = array[right];
	  array[right] = temp;
	  left++;
	  right--;
	}
    }
  while (left <= right);
  if (llimit < right)
    {
      quicksort (array, llimit, right);
    }
  if (left < rlimit)
    {
      quicksort (array, left, rlimit);
    }
}


/****************************************************/
/* Compare two process using their particle content */
/* 1 - the processes coinside                       */
/* 0 - the processes are different                  */
/****************************************************/
static int
proc (process_info f, process_info s)
{
  int i, j, k;
  int num;
  int fkf[2];
  int skf[2];
  int *pr1;
  int *pr2;

  if (f.Nparton != s.Nparton) return 0;

  for (i = 0; i < f.Nparton; i++)
    {
      if (1 == f.parton[i].in) fkf[0] = f.parton[i].KF;
      if (2 == s.parton[i].in) fkf[1] = f.parton[i].KF;
      if (1 == s.parton[i].in) skf[0] = s.parton[i].KF;
      if (2 == s.parton[i].in) skf[1] = s.parton[i].KF;
    }
  if ((skf[0] != fkf[0]) || (skf[1] != fkf[1])) return 0;

  pr1 = malloc (f.Nparton * sizeof(int));
  pr2 = malloc (f.Nparton * sizeof(int));
  j = k = 0;
  for (i = 0; i < f.Nparton; i++)
    {
      if (0 != f.parton[i].out) {pr1[j] = f.parton[i].KF;++j;}
      if (0 != s.parton[i].out) {pr2[k] = s.parton[i].KF;++k;}
    }

  if (j != k) {
    free (pr1);
    free (pr2);
    return 0;
  }
  num = j;

  for (i = 0; i < num; i++)
    {
      if (0 == pr1[i])
  	{
  	  fprintf (stderr, "mix: there is unkonwn particle '%s' in process %s\n", f.parton[i].name, f.name);
  	  pr1[i] = 1000000000;
  	}
      if (0 == pr2[i])
  	{
  	  fprintf (stderr, "mix: there is unkonwn particle '%s' in process %s\n", s.parton[i].name, s.name);
  	  pr2[i] = 1000000000;
  	}
    }

  quicksort (pr1, 0, num - 1);
  quicksort (pr2, 0, num - 1);

  for (i = 0; i < num; ++i)
    {
      if (pr1[i] != pr2[i]) {
        free (pr1);
        free (pr2);
        return 0;
      }
    }

  free (pr1);
  free (pr2);
  return 1;
}

int
compare_processes (int nf, tags ** tagbase, proc_pos * pr_pos)
{
  int i, j, k, n;
  int subtot;
  int result;
  int Nproc;
  string_comnd com;
  process_info *pr;

  if (com.name == NULL || com.value == NULL || pr == NULL)
    {
      fprintf (stderr, "Error: problems with memory allocation");
      return 0;
    }

  subtot = 0;
  for (i = 0; i < nf; ++i) {
    strcpy (com.name, "Nproc");
    get_tag_with1com (0, tagbase[i], "total", &com);
    subtot += atoi (com.value);
  }
  pr = malloc (subtot * sizeof (process_info));

  subtot = 0;
  for (i = 0; i < nf; ++i)
    {
      strcpy (com.name, "Nproc");
      get_tag_with1com (0, tagbase[i], "total", &com);
      Nproc = atoi (com.value);
      for (j = 0; j < Nproc; j++)
        {
          strcpy (com.name, "ID");
          sprintf (com.value, "%i", j + 1);
          n = get_tag_with_exactcom (0, tagbase[i], "process", com);
          if (n != -1)
            {
              pr[j + subtot].Nparton = get_ival (0, "Nparton", tagbase[i]->tag[n]);
              strcpy (pr[j + subtot].name, get_cval (0, "name", tagbase[i]->tag[n]));
            }
          n = -1;
          strcpy (com.name, "IDprocess");
          for (k = 0; k < pr[j + subtot].Nparton; k++)
            {
              n = get_tag_with_exactcom (n + 1, tagbase[i], "parton", com);
              if (n != 1)
                {
                  pr[j + subtot].parton[k].in = get_ival (0, "in", tagbase[i]->tag[n]);
                  pr[j + subtot].parton[k].out = get_ival (0, "out", tagbase[i]->tag[n]);
                  pr[j + subtot].parton[k].in = get_ival (0, "in", tagbase[i]->tag[n]);
                  pr[j + subtot].parton[k].out = get_ival (0, "out", tagbase[i]->tag[n]);
                  pr[j + subtot].parton[k].KF = get_ival (0, "KF", tagbase[i]->tag[n]);
                  strcpy (pr[j + subtot].parton[k].name, get_cval (0, "name", tagbase[i]->tag[n]));
                }
            }

          pr_pos[j + subtot].nfile = i;
          n = get_tag_with_exactcom (0, tagbase[i], "n_event", com);
          if (n != -1)
            pr_pos[j + subtot].Nevents = get_ival (0, "N", tagbase[i]->tag[n]);
        }
      subtot += Nproc;
    }

  for (i = 0; i < subtot - 1; i++)
    {
      for (j = i + 1; j < subtot; j++)
        {
          result = proc (pr[i], pr[j]);
          if (result)
            {
              fprintf (stderr, "mix: You are mixing identical processes!\n");
              fprintf (stderr, "     processes: %i and %i\n", i + 1, j + 1);
            }
        }
    }
  free (pr);
  return 0;
}

static int
proto_crucial (elementary_tag f, elementary_tag s)
{
  return 1;
}

static int
ele_crucial (elementary_tag f, elementary_tag s)
{
  return 1;
}

static int
photo_crucial (elementary_tag f, elementary_tag s)
{
  return 1;
}

int
compare_extra_info (int nf, tags ** t, process_ * pr)
{
  int l, i, j, k;
  int res, n;
  string_comnd com;
  elementary_tag *beam[nf][2];
  elementary_tag *strfun[nf][2];

  if (com.name == NULL || com.value == NULL)
    {
      fprintf (stderr, "Error: problems with memory allocation");
      return 0;
    }
  for (j = 0; j < nf; j++)
    {
      for (i = 0; i < 2; i++)
        {
          strcpy (com.name, "ID");
          sprintf (com.value, "%i", i + 1);
          n = get_tag_with_exactcom (0, t[j], "beam", com);
          if (n != -1)
            {
              beam[j][i] = t[j]->tag[n];
            }
          strcpy (com.name, "IDbeam");
          n = get_tag_with_exactcom (0, t[j], "strfun", com);
          if (n != -1)
            {
              strfun[j][i] = t[j]->tag[n];
            }
        }
    }

  for (l = 0; l < 2; l++)
    {
      if (abs (pr->proc_info.IDBeamUP[l]) == 2212)
        {
          for (i = 0; i < 2; i++)
            {
              for (j = 0; j < nf - 1; j++)
                {
                  for (k = j + 1; k < nf; k++)
                    {
                      res = proto_crucial (*strfun[j][i], *strfun[k][i]);
                      if (!res)
                        {
                          fprintf (stderr, "\nError: proton extra commands\n");
                          exit (2);
                        }
                    }
                }
            }
          return 1;
        }
      else if (abs (pr->proc_info.IDBeamUP[l]) == 11)
        {
          for (i = 0; i < 2; i++)
            {
              for (j = 0; j < nf - 1; j++)
                {
                  for (k = j + 1; k < nf; k++)
                    {
                      res = ele_crucial (*strfun[j][i], *strfun[k][i]);
                      if (!res)
                        {
                          fprintf (stderr, "\nError: electron extra commands\n");
                          exit (2);
                        }
                    }
                }
            }
          return 1;
        }
      else if (abs (pr->proc_info.IDBeamUP[l]) == 22)
        {
          for (i = 0; i < 2; i++)
            {
              for (j = 0; j < nf - 1; j++)
                {
                  for (k = j + 1; k < nf; k++)
                    {
                      res = photo_crucial (*strfun[j][i], *strfun[k][i]);
                      if (!res)
                        {
                          fprintf (stderr, "\nError: photon extra commands\n");
                          exit (2);
                        }
                    }
                }
            }
          return 1;
        }
      else
        {
          return 0;
        }
    }
  return 0;
}
