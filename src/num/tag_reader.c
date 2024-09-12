/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"

#include "evnt_tools.h"
#include "tag_reader.h"

static FILE *evfile = NULL;

int char_switch = 0;

static int
command_reader (char * command_buffer)
{
  char symb = '0';
  char buff[MAX_COMMAND_LENTH];
  int i;

  for (i = 0;; i++)
    {
      symb = new_symbol (evfile);
      if (i == MAX_COMMAND_LENTH - 1)
        {
          fprintf (stderr, "\nWarning: size of the command exceeds limit (%i tokens)! Command: %s\n",
                   MAX_COMMAND_LENTH, command_buffer);
        }

      if (symb != ',' && i < MAX_COMMAND_LENTH)
        {
          buff[i] = symb;
        }
      else
        {
          strncpy (command_buffer, buff, i); command_buffer[i] = 0;
          return 0;
        }
    }
}


static int
tagname_reader (midstr tagname)
{
  int i;
  char symb;

  int max_tg_l = sizeof (midstr) / sizeof (char);
  for (i = 0;; i++)
    {
      symb = new_symbol (evfile);
      if (i == max_tg_l)
        {
          fprintf (stderr, "\nWarning: size of the tag name exceeds limit (%i tokens)! Tag name: %s\n",
                   max_tg_l, tagname);
        }
      if (symb != ':' && i <= max_tg_l)
        {
          tagname[i] = symb;
        }
      else
        {
          tagname[i] = '\0';
          return 0;
        }
    }
}


static void
tag_reader (elementary_tag * tag)
{
  int i;

  if (tagname_reader (tag->tagname))
    {
      fprintf (stderr, "\n***Error! function: command_reader, error with tag syntax in tag %s\n", tag->tagname);
      exit (1);
    }

  for (i = 0; i < MAX_COMMAND_NUMBER; i++)
    {
      char symb = new_symbol (evfile);
      if (symb == ';')
        {
          tag->tagsize = i;
          return;
        }
      else
        {
          fseek (evfile, -1, 1);
        }
      {
        char * com = malloc (MAX_COMMAND_LENTH * sizeof (char));
        if (command_reader (com)) {
          fprintf (stderr, "\n***Error! function: command_reader, error with command syntax in tag %s\n", tag->tagname);
          exit (1);
        }
        strcpy (tag->commands[i], com);
        free (com);
      }
      tag->commands = realloc (tag->commands, (2 + i) * sizeof (char *));
      tag->commands[i + 1] = malloc (MAX_COMMAND_LENTH * sizeof (char));
    }

  fprintf (stderr, "\nWarning: number of commands in the tag %s exceeds limit (%i commands)!\n",
           tag->tagname, MAX_COMMAND_NUMBER);
  for (;;)
    {
      char symb = new_symbol (evfile);
      if (symb == ';')
        return;
    }
}


extern int
cup_reader (FILE * inFile, tags * all_tags)
{
  int j = 0, y;

  evfile = inFile;
  for (y = 0;; y++)
    {
      char symbol = new_symbol (evfile);
      if (symbol == '#') {
        symbol = new_symbol (evfile);
        if (symbol != '#') {
          fprintf (stderr, "CPYTH format error ...");
          return -1;
        }
      } else {
        fprintf (stderr, "This file is not in CPYTH format!");
        return -2;
      }
      {
        elementary_tag * tag = init_tag (1);
        tag_reader (tag);
        if (strlen (tag->tagname) > 0) {
          if (NULL != all_tags->tag[j]) all_tags->tag[j] = NULL;
          all_tags->tag[j] = tag;
          j++;
          all_tags->tag = realloc (all_tags->tag, (j + 1) * sizeof (elementary_tag));
        } else {
          free_tag (tag);
        }
      }
      symbol = new_symbol (evfile);
      if (symbol == ';' || y > 10000)
        {
          all_tags->number_of_tags = j;
          return 0;
        }
      else
        fseek (evfile, -1, SEEK_CUR);
    }
}


elementary_tag *
init_tag (int n)
{
  int i;
  elementary_tag *tag;

  tag = malloc (sizeof (elementary_tag));
  tag->used = 0;
  tag->tagsize = n;
  tag->commands = malloc (n * sizeof (char *));
  for (i = 0; i < n; i++)
    {
      tag->commands[i] = malloc (MAX_COMMAND_LENTH * sizeof (char));
    }
  return tag;
}

void 
free_tag (elementary_tag * t)
{
  int i;

  for (i = 0; i < t->tagsize; i++) {
    if (NULL != t->commands[i]) free (t->commands[i]);
  }
  free (t->commands);
  return;
}

extern elementary_tag *
init_tag_unallocated (int n)
{
  int i;
  elementary_tag * tag = malloc (sizeof (elementary_tag));
  tag->used = 0;
  tag->tagsize = n;
  tag->commands = malloc (n * sizeof (char *));
  for (i = 0; i < n; i++)
    {
      tag->commands[i] = NULL;
    }
  return tag;
}

tags *
init_cap (int in)
{
  tags *the_tags;
  int i;

  the_tags = malloc (sizeof (tags));
  the_tags->number_of_tags = in;
  the_tags->tag = malloc (in * sizeof (elementary_tag));
  for (i = 0; i < in; i++)
    the_tags->tag[i] = init_tag (1);

  return the_tags;
}
