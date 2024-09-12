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
#include "tag_reader.h"
#include "tag_writer.h"
#include "tag_parser.h"
#include "validity.h"

int 
syntax_validity (tags * tagbase)
{
  return 1;
}

int 
structure_validity (tags * tagbase)
{
  int i, n;
  int total_tag = 0;
  int beam_tag = 0;
  int strfun_tag = 0;
  int nproc;
  int Nparton;
  int parton = 0;
  int n_event = 0;
  int qcd = 0;
  char *tagname;
  string_comnd com;

/****************************************************************/
/*              Check presence of mandatory tags                */
  for (i = 0; i < tagbase->number_of_tags; i++)
    {
      tagname = malloc (strlen (tagbase->tag[i]->tagname) * sizeof (char));
      strcpy (tagname, tagbase->tag[i]->tagname);
      if (!strcmp (tagname, "total"))
	total_tag++;
      if (!strcmp (tagname, "beam"))
	beam_tag++;
      if (!strcmp (tagname, "strfun"))
	strfun_tag++;
      free (tagname);
    }

  if (total_tag != 1)
    {
      fprintf (stderr, "Validity error: number of total tags %i\n", total_tag);
      return -1;
    }
  if (beam_tag != 2)
    {
      fprintf (stderr, "Validity error: number of beam tags %i\n", beam_tag);
      return -1;
    }
  if (strfun_tag != 2)
    {
      fprintf (stderr, "Validity error: number of strfun tags %i\n", strfun_tag);
      return -1;
    }

/****************************************************************/
/*              Check consistency of a tagbase                  */
  strcpy (com.name, "Nproc");
  n = get_tag_with1com (0, tagbase, "total", &com);
  if (n != -1)
    {
      nproc = atoi (com.value);
    }
  else
    {
      fprintf (stderr, "Validity error: total tag does not contain Nproc command.\n");
      return -1;
    }

  for (i = 1; i <= nproc; i++)
    {
      strcpy (com.name, "ID");
      sprintf (com.value, "%i", i);
      n = get_tag_with_exactcom (0, tagbase, "process", com);
      if (n == -1)
	{
	  fprintf (stderr, "Validity error: problems with process ID numbers\n");
	  return -1;
	}
      else
	{
	  Nparton = get_ival (0, "Nparton", tagbase->tag[n]);
	}

      strcpy (com.name, "IDprocess");

      parton = 0;
      n_event = 0;
      qcd = 0;
      for (i = 0; i < tagbase->number_of_tags; i++)
	{
	  tagname = malloc (strlen (tagbase->tag[i]->tagname) * sizeof (char));
	  strcpy (tagname, tagbase->tag[i]->tagname);
	  if (!strcmp (tagname, "parton"))
	    {
	      if (tag_contain_exactcom (com, tagbase->tag[i]))
		parton++;
	    }
	  if (!strcmp (tagname, "n_event"))
	    {
	      if (tag_contain_exactcom (com, tagbase->tag[i]))
		n_event++;
	    }
	  if (!strcmp (tagname, "QCDinfo"))
	    {
	      if (tag_contain_exactcom (com, tagbase->tag[i]))
		qcd++;
	    }
	  free (tagname);
	}

      if (Nparton != parton)
	{
	  fprintf (stderr, "Validity error: number of parton tags and Nparton do not");
	  fprintf (stderr, "              : coincide! parton tags=%i, Nparton=%i\n", parton, Nparton);
	  return -1;
	}
      if (n_event != 1)
	{
	  fprintf (stderr, "Validity error: number of n_event tags in process %i is %i\n", i, n_event);
	  return -1;
	}
      if (qcd != 1)
	{
	  fprintf (stderr, "Validity error: number of QCDinfo tags in process %i is %i\n", i, qcd);
	  return -1;
	}
    }
  return 1;
}

int 
physics_validity (tags * tagbase)
{
  return 1;
}
