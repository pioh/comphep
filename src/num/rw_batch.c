/*
* Copyright (C) 2008-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ------------------------------------------------------
*/
#include <unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#include "service2/include/syst.h"

#include "rw_sess.h"
#include "rw_batch.h"

int 
w_batch_session (FILE * mode)
{
  fprintf (stdout,  "configurator: still not realised\n");
  return 0;
}


int 
r_batch_session (FILE * mode)
{
  fprintf (stdout,  "configurator: still not realised\n");
  return 0;
}

int 
init_batch_session (void)
{
  fprintf (stdout,  "configurator: still not realised\n");
  return 0;
}

int 
get_process_name_safe (char * procname)
{
  return 0;
}
