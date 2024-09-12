/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PHYSICS_
#define __PHYSICS_

#include "service2/include/tptcmac.h"
#include "model.h"
#include "diagrams.h"

#define whohowMAX 32
#define mult_whohowMAX 16
typedef struct __mult_whohow__
  {
    int who[64];		/* Type chastitsi, obichno nomer chastitsi v localbase */
    int how;			/* Colichestvo chastits, obichno nomer chastitsi v localbase */
    int type;			/* Special parameter. Vveden dlya peredachi usloviy excluding/keeping */
  }
__mult_whohow__;

/*
type defines logical connective for excuding/keeping:
  type =  2 -> exclude/keep;
  type =  1 -> operation ">";
  type =  0 -> the hadron does not exist;
  type = -1 -> operation "!=";
  type = -2 -> operation "=";
  type = -3 -> operation "<";
*/

typedef struct __whohow__
  {
    int who;			/* Type chastitsi, obichno nomer chastitsi v localbase */
    int how;			/* Colichestvo chastits, obichno nomer chastitsi v localbase */
    int type;			/* Special parameter. Vveden dlya peredachi usloviy excluding/keeping */
  }
__whohow__;

typedef __mult_whohow__ multwhohow[mult_whohowMAX];
typedef __whohow__ whohow[whohowMAX];

extern multwhohow exclude_composit_list;
extern multwhohow keep_composit_list;

extern whohow exclude_list;
extern whohow keep_list;

extern __beam__ beam[2];

extern void nilprtcl (whohow p_list);

/* ======================================================== */
typedef char shortname[6];
typedef shortname prtclsarray[MAXINOUT];
extern void getprtcls (char *txt, prtclsarray pnames);
/*===========================================================*/

extern unsigned subproc_f, subproc_sq;

extern char modelmenu[STRSIZ];
extern int maxmodel;

extern int nsub, ndiagr;

extern int getModelNumberSymb (void);
extern void setModelNumberSymb (int num);

#endif
