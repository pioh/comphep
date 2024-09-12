/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef _CREADER_H_
#define _CREADER_H_

#define MAX_COMMAND_LENTH     1024
#define MAX_COMMAND_NUMBER    256
#define MAX_NUMBER_KNOWN_TAGS 256


typedef char middlestr[MAX_COMMAND_LENTH];

typedef struct string_comnd
  {
    middlestr name;
    middlestr value;
  }
string_comnd;

typedef struct elementary_tag
  {
    int tagsize;
    int used;
    middlestr tagname;
    char **commands;
  }
elementary_tag;

typedef struct tags
  {
    int number_of_tags;
    elementary_tag **tag;
  }
tags;

extern elementary_tag * init_tag (int n);
extern void free_tag (elementary_tag * t);

extern elementary_tag * init_tag_unallocated (int n);
extern tags *init_cap (int in);

extern int cup_reader (FILE * inFile, tags * all_tags);

#endif /* c-reader.h  */
