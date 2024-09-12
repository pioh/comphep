/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __DIAGRAMS__
#define __DIAGRAMS__

#include"model.h"

/*========== Pukhov representation ================== */

#define ldiagram (2*MAXINOUT-3)
/*  maximum number of particles in diagram  */

typedef particleNumType decayDiagram[ldiagram];

typedef struct adiagram
  {
    decayDiagram dgrm0;
    char delMark;
    int nsub;
  }
adiagram;

typedef struct csdiagram
  {
    decayDiagram dgrm1, dgrm2;
    int lnk[MAXINOUT];
    int mult;
    unsigned del;
    char status;		/* -2-outOfMemory,-1-deleted,
				   0-Rest, 1-calculated,2-Zero  */
    int nsub;
  }
csdiagram;

#define maxvert (MAXINOUT-2)	/* maximal # of verteces in amplitude */
#define nullvert 253		/*# of next vert for  unused edge     */

/*================== Taranov representation =========== */

#define IN_PRTCL (int)1
#define OUT_PRTCL (int)2

typedef struct vertlink
  {
    int vno, edno;		/* # of vert, # of edge in vert (slot) */
  }
vertlink;

typedef struct edgeinvert
  {
    int lorentz;
    int moment;
    int prop;  /* 1 - incoming prtcl, 2 - outgoing prtcl */
    int partcl;
    vertlink nextvert;
    vertlink link;
  }
edgeinvert;

typedef edgeinvert vert0[MAXVALENCE];


typedef struct vampl
  {
    int size, outno;		/*  how many  verts and outgoing edges  */
    int valence[maxvert];
    vertlink outer[MAXINOUT];	/*  adresses of external edges */
    vert0 vertlist[maxvert];	/*  array of verts  */
  }
vampl;


typedef struct vcsect
  {
    int sizel;
    int sizet;
    long symnum;
    long symdenum;
    long clrnum;
    long clrdenum;
    int valence[2 * maxvert];	/* 1..4 */
    vert0 vertlist[2 * maxvert];
  }
vcsect;

extern vcsect vcs;

extern void transfdiagr (csdiagram * diag, vcsect * vcs);
extern void mkverts (decayDiagram diag1, vampl * vlist1);
extern void InOutPrtclsNumb (decayDiagram a, int *numb, int sort);
extern void proccessName (decayDiagram a, char *txt);
extern void decompose (vcsect vcs, vampl * left, vampl * right);
extern void mkcsections (csdiagram * diagr, vcsect * vcs);
extern void printDiagram (vampl * vlist);
extern void printCsDiagram (vcsect * vlist);
#endif
