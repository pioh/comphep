/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structures.h"
#include "tag_reader.h"
#include "tag_parser.h"
#include "fill.h"


static int 
number_proc (tags * t)
{
  int res;
  string_comnd com;

  strcpy (com.name, "Nproc");
  get_tag_with1com (0, t, "total", &com);
  res = atoi (com.value);
  return res;
}

static void 
get_gen (tags * t, int *IDBeam, double *energy, int *PDFgroup, int *PDFid)
{
  int i, n;
  string_comnd com;

  for (i = 0; i < 2; i++)
    {
      strcpy (com.name, "ID");
      sprintf (com.value, "%i", i + 1);
      n = get_tag_with_exactcom (0, t, "beam", com);
      if (n == -1)
	{
	  fprintf (stderr, "\nError: unknown beam tag\n");
	  exit (3);
	}
      IDBeam[i] = get_ival (0, "KF", t->tag[n]);
      energy[i] = get_fval (0, "energy", t->tag[n]);

      strcpy (com.name, "IDbeam");
      n = get_tag_with_exactcom (0, t, "strfun", com);
      if (n == -1)
	{
	  fprintf (stderr, "\nError: unknown strfun tag\n");
	  exit (3);
	}
      PDFid[i] = get_ival (0, "PDFid", t->tag[n]);
      PDFgroup[i] = get_ival (0, "PDFgr", t->tag[n]);
      if (abs (IDBeam[i]) != 2212)
	{
	  PDFgroup[i] = -1;
	  PDFid[i] = -1;
	}
    }
}

int 
fill_LH_structures (int nf, tags ** tbase, process_ * prUP, const char names[], int lenth)
{
  int i, j;
  int n, len;
  int Nproc;
  int subtot;
  int IDBeam[2], PDFgroup[2], PDFid[2];
  int t_IDBeam[2], t_PDFgroup[2], t_PDFid[2];
  double energy[2], t_energy[2];
  string_comnd com;

  get_gen (tbase[0], IDBeam, energy, PDFgroup, PDFid);
  for (i = 0; i < 2; i++)
    {
      prUP->proc_info.IDBeamUP[i] = IDBeam[i];
      prUP->proc_info.PDFgroupUP[i] = PDFgroup[i];
      prUP->proc_info.PDFidUP[i] = PDFid[i];
      prUP->proc_info.energyUP[i] = energy[i];
    }

  for (i = 1; i < nf; i++)
    {
      get_gen (tbase[i], t_IDBeam, t_energy, t_PDFgroup, t_PDFid);
      for (j = 0; j < 2; j++)
	{
	  if ((t_IDBeam[j] != IDBeam[j]) ||
	      (t_energy[j] != energy[j]) ||
	      (t_PDFgroup[j] != PDFgroup[j]) ||
	      (t_PDFid[j] != PDFid[j]))
	    {
	      fprintf (stderr, "Diffirent general tasks in %i and %i files:\n", 1, i + 1);
	      fprintf (stderr, "   file %i: IDBeam[0]=%i, IDBeam[1]=%i\n", i + 1, t_IDBeam[0], t_IDBeam[1]);
	      fprintf (stderr, "	     energy[0]=%f, energy[1]=%f\n", t_energy[0], t_energy[1]);
	      fprintf (stderr, "	     PDFgroupUP[0]=%i PDFgroupUP[1]=%i\n",
		       t_PDFgroup[0], t_PDFgroup[1]);
	      fprintf (stderr, "	     PDFid[0]=%i, PDFid[1]=%i\n", t_PDFid[0], t_PDFid[1]);
	      fprintf (stderr, " -----------------------------------\n");
	      fprintf (stderr, "   file %i: IDBeam[0]=%i, IDBeam[1]=%i\n", 1, IDBeam[0], IDBeam[1]);
	      fprintf (stderr, "	     energy[0]=%f, energy[1]=%f\n", energy[0], energy[1]);
	      fprintf (stderr, "	     PDFgroupUP[0]=%i PDFgroupUP[1]=%i\n",
		       PDFgroup[0], PDFgroup[1]);
	      fprintf (stderr, "	     PDFid[0]=%i, PDFid[1]=%i\n", PDFid[0], PDFid[1]);
//         return -1;  
	    }
	}
    }

  subtot = 0;
  for (i = 0; i < nf; i++)
    {
      subtot += number_proc (tbase[i]);
    }
  prUP->maps = malloc (subtot * sizeof (process_map));

  subtot = 0;
  strcpy (com.name, "ID");
  for (i = 0; i < nf; i++)
    {
      len = strlen (names + i * lenth) + 5;
      Nproc = number_proc (tbase[i]);
      for (j = 0; j < Nproc; j++)
	{
	  if (MAXproc < j + subtot) {
	    fprintf (stderr, "tag_reader (error): Number of subprocesses is too large! \n");
	    fprintf (stderr, "                    Enlarge the paramerter MAXproc in comphep-XXX/src/num/include/LesHouches.h\n");
	    fprintf (stderr, "                    and recompile CompHEP\n");
	    exit (-10);
	  }
	  sprintf (com.value, "%i", j + 1);
	  n = get_tag_with_exactcom (0, tbase[i], "process", com);

	  prUP->maps[j + subtot].file = i;
	  prUP->maps[j + subtot].number = j;
	  prUP->maps[j + subtot].filename = malloc (len * sizeof (char));
	  strcpy (prUP->maps[j + subtot].filename, names + i * lenth);
	  prUP->maps[j + subtot].used = 1;

	  prUP->proc_info.switchUP = get_ival (0, "master", tbase[i]->tag[n]);
	  prUP->proc_info.listprocRUP[j + subtot] = j + subtot;
	  prUP->proc_info.crossecUP[j + subtot] = get_fval (0, "CrosSec", tbase[i]->tag[n]);
	  prUP->proc_info.errorcsUP[j + subtot] = get_fval (0, "CrosSecErr", tbase[i]->tag[n]);
	  prUP->proc_info.maxweightUP[j + subtot] = 1.0;
	  strcpy (prUP->proc_info.PRnameUP[j + subtot], get_cval (0, "name", tbase[i]->tag[n]));
	}
      subtot += Nproc;
    }
  prUP->proc_info.NprocRUP = subtot;

  return 0;
}



int 
fill_heprup (processUP * prUP, tags * head)
{
  int i, n;
  int switch1;
  int IDBeam[2], PDFgroup[2], PDFid[2];
  double energy[2];
  string_comnd com;

  get_gen (head, IDBeam, energy, PDFgroup, PDFid);
  for (i = 0; i < 2; i++)
    {
      prUP->IDBeamUP[i] = IDBeam[i];
      prUP->PDFgroupUP[i] = PDFgroup[i];
      prUP->PDFidUP[i] = PDFid[i];
      prUP->energyUP[i] = energy[i];
    }

  strcpy (com.name, "ID");

  prUP->NprocRUP = number_proc (head);
  for (i = 0; i < prUP->NprocRUP; i++)
    {
      sprintf (com.value, "%i", i + 1);
      n = get_tag_with_exactcom (0, head, "process", com);
      switch1 = get_ival (0, "master", head->tag[n]);
      if (0 == i)
	prUP->switchUP = switch1;
      else
	{
	  if (prUP->switchUP != switch1)
	    {
	      fprintf (stderr, "Warning: switch parameter for subprocess %i differ\n"
		       "         from one for subprocess 1.", i + 1);
	    }
	}
      prUP->listprocRUP[i] = i;
      prUP->crossecUP[i] = get_fval (0, "CrosSec", head->tag[n]);
      prUP->errorcsUP[i] = get_fval (0, "CrosSecErr", head->tag[n]);
      prUP->maxweightUP[i] = 1.0;
      strcpy (prUP->PRnameUP[i], get_cval (0, "name", head->tag[n]));
    }

  return 0;
}

int 
prefill_hepeup (eventUP * ev, tags * head)
{

  return 0;
}


int 
fill_hepeup (eventUP * ev, char *evstr, int pnum)
{

  return 0;
}
