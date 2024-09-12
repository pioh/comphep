/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#include "LesHouches.h"

#ifndef _STRUCTURES_H_
#define _STRUCTURES_H_

#define MAX_COMMAND_LENTH     1024
#define MAX_COMMAND_NUMBER    256
#define MAX_NUMBER_KNOWN_TAGS 256
#define STRLEN 2048

/********************************************/
typedef struct proc_pos
  {
    int nfile, Nevents;
  }
proc_pos;

typedef struct _beam_
  {
    double energy, mass;
    int KF, IDbeam;
    char name[30];
    int N_Extra_Com;
    char ExtraInfo[50][50];
  }
_beam_;

typedef struct _strfun_
  {
    char name[30], version[30];
    int N_Extra_Com, PDFid, PDFgr;
    char ExtraInfo[50][50];
  }
_strfun_;

typedef struct _nevent_
  {
    int IDprocess, N, multN, origN;
    double maxW;
    char orig_file[200];
  }
_nevent_;

typedef struct _parton_
  {
    int IDprocess, in, out, KF;
    char name[10];
    double mass;
  }
_parton_;

typedef struct _QCDinfo_
  {
    int IDprocess, NL, Nflavour;
    double QCDLambda;
  }
_QCDinfo_;


/********************************************/
/*typedef struct general_info
   {
   _beam_ beam[2];
   _strfun_ strfun[2];
   } general_info; */

typedef struct process_info
  {
    int ID;
    char name[200];
    char generator[20];
    char version[20];
    double CrosSec, CrosSecErr;
    int Nparton, master;

    _nevent_ n_event;
    _parton_ parton[10];
    _QCDinfo_ QCDinfo;
  }
process_info;

#endif /* structures.h  */
