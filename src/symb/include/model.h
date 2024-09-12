/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __MODEL__
#define __MODEL__


/* ================== variables ==================== */
#define VAR_NAME_SIZE 7
#define strongconst "GG"

typedef struct varrec
  {
    char varname[VAR_NAME_SIZE];
    int able;
    double varvalue;
    char *func;
  }
varrec;
typedef struct varrec *varlist;

extern int nmodelvar;
extern varlist modelvars;

/*=================== particles ==================== */

#define P_NAME_SIZE 6
#define MAXVALENCE 4

typedef short particleNumType;
typedef struct modeofdecay
  {
    struct modeofdecay *next;
    particleNumType part[MAXVALENCE - 1];
  }
modeofdecay;
typedef struct modeofdecay *decaylink;

typedef struct prtcl_base
  {
    char name[P_NAME_SIZE];
    int anti, spin;
    double mass;
    char massidnt[VAR_NAME_SIZE], imassidnt[VAR_NAME_SIZE];
    int cdim;
    int hlp;
    char *latex;
    decaylink top;
  }
prtcl_base;

extern int nparticles;		/* Number particles in model */
extern prtcl_base *prtclbase;

extern int pseudop (int np);
extern int fermionp (int p);
extern int a_fermionp (int p);
extern int bosonp (int p);
extern int vectorp (int p);
extern int zeromass (int p);
extern int photonp (int p);
extern int ghostp (int p);
extern int ghostmother (int j);
extern int gaugep (int j);
extern int prtclname (int n, vshortstr nm);

extern int locateinbase (char *name);

/* ================== lagrangian ================ */
typedef int arr4byte[4];
typedef struct algvert
  {
    struct algvert *next;
    int fields[MAXVALENCE];
    int perm[MAXVALENCE];
    int factor;
    char *comcoef;
    char *description;
  }
algvert;
typedef struct algvert *algvertptr;

extern algvertptr lgrgn;

/* ================== composite particles ================ */
typedef struct comppart
  {
    shortstr name;     /* composite particle name */
    int how;           /* how many partons is in the composite particle */
    shortstr cpart;    /* string with model particles */
  }
comppart;
typedef struct comppart *cparticles;

extern int n_cpart;
extern cparticles cpartbase;

/* ================== beams  ================ */
typedef struct hadron
  {
    shortstr name;          /* hadron name */
    double mass;            /* hadron mass */
    int how;                /* how many partons is in the hadron */
    int parton[64];         /* parton is model paticle number in particle base of the model */
    shortstr sf_name;       /* stucture function name */
    shortstr full_sf_name;  /* stucture function name */
    int sf_set;             /* stucture function set number */
    int sf_mem;             /* stucture function mem number */
  }
hadron;

typedef struct hadron * _hadron_;

extern int n_hadron;
extern hadron hadrons[MAXINOUT];
extern _hadron_  hadronbase;

typedef struct __beam__
{
  hadron h;
  double energy;
} __beam__;

extern __beam__ beam[2];

/* ================== beams  ================ */
typedef struct strfun
  {
    shortstr name;     /* stucture function name */
    int set;           /* stucture function LHAPDF set number */
    int mem;           /* stucture function LHAPDF mem number */
    int num;           /* stucture function internal CompHEP number */
    shortstr file;     /* stucture function file */
  }
strfun;

typedef struct strfun * _strfun_;

extern int n_strfun;
extern _strfun_  strfunbase;

#endif
