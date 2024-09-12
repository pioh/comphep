/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tag_reader.h"
#include "tag_parser.h"

/**************************************************/
int
intlen (int src)
{
  int i = 0;
  int digit = 1;

  while (abs (src / digit) >= 1.0)
    {
      digit *= 10;
      i++;
    }
  return i;
}


static char * ttrim (char * p)
{
  char * pos;
  char * f = malloc (strlen (p) * sizeof (char));
  f[0] = 0;
  pos = strtok (p, " ");
  while (pos) {
    strcpy (f, pos);
    pos = strtok (NULL, " ");
  }
  strcpy (p, f);
  free (f);
  return p;
}


/**************************************************/
static int 
command_parser (const char * s, string_comnd * com)
{
  char * val = strchr (s, '=');
  if (NULL != val) {
    int len = strlen (s) - strlen (val);
    char * shtr = strchr (val + 1, '\'');
    if (NULL != shtr) {
       strcpy (com->value, shtr);
/*
       int i;
       int len1;
       strcpy (com->value, shtr + 1);              // remove first prime
       len1 = strlen (com->value);
       for (i = len1; '\'' != com->value[i]; --i)
         continue;
       com->value[i] = 0;                   // remove last prime
*/
    } else {
      strcpy (com->value, val + 1);
      if (NULL != strchr (com->value, ' ')) ttrim (com->value);
    }
    strncpy (com->name, s, len); com->name[len] = 0;
    if (NULL != strchr (com->name, ' ')) ttrim (com->name);
    return 0;
  }
  return -1;
}

/**************************************************/
int
get_tag (int first, tags * base, char *name)
{
  int i, n;

  n = base->number_of_tags;
  if (first <= n)
    {
      for (i = first; i < n; i++)
        {
          if (!strcmp (base->tag[i]->tagname, name))
            {
              return i;
            }
        }
    }
  return -1;
}


/**************************************************/
int
get_tag_with1com (int first, tags * base, char tagname[2048], string_comnd * com)
{
  int i, j;
  int n = base->number_of_tags;
  if (first <= n)
    {
      for (i = first; i < n; i++)
        {
          if (!strcmp (base->tag[i]->tagname, tagname))
            {
              for (j = 0; j < base->tag[i]->tagsize; j++)
                {
                  string_comnd temp;
                  command_parser (base->tag[i]->commands[j], &temp);
                  if (!strcmp (com->name, temp.name))
                    {
                      strcpy (com->value, temp.value);
                      return i;
                    }
                }
            }
        }
    }
  return -1;
}


/**************************************************/
/*
static int
get_tag_with_exactcom_old (int first, tags * base, char tagname[2048], string_comnd com)
{
  int i, j;
  int n = base->number_of_tags;
  int len1 = strlen (com.name);
  int len2 = strlen (com.value);
  if (first <= n)
    {
      for (i = first; i < n; i++)
        {
          if (!strcmp (base->tag[i]->tagname, tagname))
            {
              for (j = 0; j < base->tag[i]->tagsize; j++)
                {
                  string_comnd temp;
                  command_parser (base->tag[i]->commands[j], &temp);
                  if (!strncmp (com.name, temp.name, len1) && !strncmp (com.value, temp.value, len2))
                    {
                      return i;
                    }
                }
            }
        }
    }
  return -1;
}
*/
/**************************************************/
int
get_tag_with_exactcom (int first, tags * base, char tagname[2048], string_comnd com)
{
  int i, j;
  int n = base->number_of_tags;
  if (first <= n)
    {
      for (i = first; i < n; i++)
        {
          if (!strcmp (base->tag[i]->tagname, tagname))
            {
              for (j = 0; j < base->tag[i]->tagsize; j++)
                {
                  string_comnd temp;
                  command_parser (base->tag[i]->commands[j], &temp);
                  if (!strcmp (com.name, temp.name) && !strcmp (com.value, temp.value))
                    {
                      return i;
                    }
                }
            }
        }
    }
  return -1;
}

/**************************************************/
int
get_tag_with_exact2com (int first, tags * base, char tagname[2048], string_comnd * com1, string_comnd * com2)
{
  int i, j;
  int in[2];
  int stop = 1;
  int n = base->number_of_tags;
  if (first > n)
    {
      return -1;
    }

  i = first - 1;
  while (stop)
    {
      i++;
      if (!strcmp (base->tag[i]->tagname, tagname))
        {
          in[0] = 0;
          in[1] = 0;
          for (j = 0; j < base->tag[i]->tagsize; j++)
            {
              string_comnd temp;
              command_parser (base->tag[i]->commands[j], &temp);
              if (!strcmp (com1->name, temp.name) && !strcmp (com1->value, temp.value))
                in[0] = i;
              if (!strcmp (com2->name, temp.name) && !strcmp (com2->value, temp.value))
                in[1] = i;
            }
          if ((in[0] && in[1]) || i == n)
            {
              stop = 0;
            }
        }
    }

  if (in[0] && in[1])
    {
      return in[0];
    }
  else
    {
      return -1;
    }
}


/**************************************************/
int
tag_contain_com (string_comnd * com, elementary_tag * t)
{
  int i;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (com->name, temp.name))
        {
          strcpy (com->value, temp.value);
          return i;
        }
    }
  return -1;
}


/**************************************************/
int
tag_contain_exactcom (string_comnd com, elementary_tag * t)
{
  int i;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (com.name, temp.name) && !strcmp (com.value, temp.value))
        {
          return 1;
        }
    }
  return 0;
}


/**************************************************/
int
add_com (string_comnd com, elementary_tag * t)
{
  int i;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (com.name, temp.name))
        {
          return -1;
        }
    }
  i = strlen (com.name) + strlen (com.value) + 3;
  t->tagsize += 1;
  t->commands = realloc (t->commands, t->tagsize * sizeof (char *));
  t->commands[t->tagsize - 1] = malloc (i * sizeof (char));
  sprintf (t->commands[t->tagsize - 1], "%s=%s", com.name, com.value);
  return 1;
}


/**************************************************/
int
remove_com (char *cname, elementary_tag * t)
{
  int i, j;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (cname, temp.name))
        {
          t->tagsize -= 1;
          for (j = i; j < t->tagsize; j++)
            {
              t->commands[j] = realloc (t->commands[j], strlen (t->commands[j + 1]) * sizeof (char));
              strcpy (t->commands[i], t->commands[j + 1]);
            }
          t->commands = realloc (t->commands, t->tagsize * sizeof (char *));
          return 1;
        }
    }
  return -1;
}


/**************************************************/
int
add_tag (elementary_tag * tag, tags * base)
{
  int i;

  for (i = 0; i < base->number_of_tags; i++)
    {
      if (!strcmp (base->tag[i]->tagname, tag->tagname))
        {
          return -1;
        }
    }
  base->number_of_tags += 1;
  base->tag = realloc (base->tag, base->number_of_tags * sizeof (elementary_tag));
  base->tag[base->number_of_tags - 1] = tag;
  return 1;
}


/**************************************************/
int
remove_tag (char *tname, tags * base)
{
  int i, j;

  for (i = 0; i < base->number_of_tags; i++)
    {
      if (!strcmp (base->tag[i]->tagname, tname))
        {
          base->number_of_tags -= 1;
          for (j = i; j < base->number_of_tags; j++)
            {
              base->tag[j] = base->tag[j + 1];
            }
          base->tag = realloc (base->tag, base->number_of_tags * sizeof (elementary_tag));
          return 1;
        }
    }
  return -1;
}

/**************************************************/
int
remove_tag_num (int num, tags * base)
{
  int i;

  if (0 > num || num > base->number_of_tags - 1) return -1;

  for (i = num; i < base->number_of_tags; ++i)
    {
      base->tag[i] = base->tag[i + 1];
    }
  --(base->number_of_tags);
  base->tag = realloc (base->tag, base->number_of_tags * sizeof (elementary_tag));

  return 1;
}

/**************************************************/
int
remove_null_tags (tags * base)
{
  int i, j;

  for (i = 0; i < base->number_of_tags; i++)
    {
      if (0 == strlen (base->tag[i]->tagname))
        {
          base->number_of_tags -= 1;
          for (j = i; j < base->number_of_tags; j++)
            {
              base->tag[j] = base->tag[j + 1];
            }
          base->tag = realloc (base->tag, base->number_of_tags * sizeof (elementary_tag));
        }
    }

  return 1;
}

/**************************************************/
int
get_ival (int n, char *cname, elementary_tag * t)
{
  int i;
  int res;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (cname, temp.name))
        {
          res = atoi (temp.value);
          return res;
        }
    }
  return -1;
}


/**************************************************/
double
get_fval (int n, char *cname, elementary_tag * t)
{
  int i;
  double res;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (cname, temp.name))
        {
          res = atof (temp.value);
          return res;
        }
    }
  return -1.0;
}


/**************************************************/
char *
get_cval (int n, char * cname, elementary_tag * t)
{
  int i;
  char * dummy = "";

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (cname, temp.name))
        {
          return temp.value;
        }
    }
  return dummy;

}


/**************************************************/
int
replace_com (string_comnd com, elementary_tag * t)
{
  int i;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (com.name, temp.name))
        {
          if (!(t->commands[i]))
            free (t->commands[i]);
          t->commands[i] = malloc ((strlen (com.name) + strlen (com.value) + 2) * sizeof (char));
          sprintf (t->commands[i], "%s=%s", com.name, com.value);
          return i;
        }
    }
  return -1;
}


/**************************************************/
int
change_ival (int newval, int shift, string_comnd com, elementary_tag * t)
{
  int i, j;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (com.name, temp.name))
        {
          char *format;
          int len = strlen (com.name) + shift + intlen (newval) + 1;
          if (!(t->commands[i]))
            free (t->commands[i]);
          t->commands[i] = malloc (len * sizeof (char));
          format = malloc ((shift + 5) * sizeof (char));
          strcpy (format, "%s=");
          for (j = 0; j < shift; j++)
            format[j + 3] = ' ';
          format[shift + 3] = '\0';
          strcat (format, "%i");
          sprintf (t->commands[i], format, temp.name, newval);
          return 1;
        }
    }
  return -1;
}


/**************************************************/
int
change_fval (double newval, string_comnd com, elementary_tag * t)
{
  int i;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (com.name, temp.name))
        {
          if (!(t->commands[i]))
            free (t->commands[i]);
          t->commands[i] = malloc ((strlen (com.name) + 13) * sizeof (char));
          sprintf (t->commands[i], "%s=%.5E", com.name, newval);
          return 1;
        }
    }
  return -1;
}


/**************************************************/
int
change_cval (char *newval, string_comnd com, elementary_tag * t)
{
  int i;

  for (i = 0; i < t->tagsize; i++)
    {
      string_comnd temp;
      command_parser (t->commands[i], &temp);
      if (!strcmp (com.name, temp.name))
        {
          if (!(t->commands[i]))
            free (t->commands[i]);
          t->commands[i] = malloc ((strlen (com.name) + strlen (newval) + 2) * sizeof (char));
          sprintf (t->commands[i], "%s=%s", com.name, newval);
          return 1;
        }
    }
  return -1;
}

/**************************************************/
int
event_parser (char *in, char **eventword)
{
  int i, n = 0;
  int len, fulllen;
  int last = 0;

  fulllen = strlen (in);
  if (13 == in[fulllen-1]) in[fulllen-1] = '\0';
  fulllen = strlen (in);
  while (last < fulllen)
    {
      char *st = strchr (in, ':');
      if (!st)
        {
          fprintf (stderr, "Error! an event does not contain colon!\n");
          return 1;
        }
      len = fulllen - strlen (st);
      eventword[n] = malloc ((len - last + 1) * sizeof (char));

      for (i = last; i < len; i++)
        eventword[n][i - last] = in[i];
      eventword[n][len - last + 1] = '\0';

      last = len + 1;
      in[len] = '_';

      n++;
    }
  return 0;
}
