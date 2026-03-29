/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __LHAPDF__
#define __LHAPDF__

typedef struct lhapdfList {
    struct lhapdfList  *next;
    char * name;                 /* title name of distribution             */
    int set;                     /* set number                             */
    int mem;                     /* member number                          */
    char * pathfile;             /* path to file where it is stored        */
    char * file;                 /* name of file where it is stored        */
    int position;                /* ordering number of parton in the list of functions */
} lhapdfList;

extern void delLhapdfList (lhapdfList * list);

extern int comphepLhapdfList (lhapdfList ** list);

/* free  memory allocated for 'list' */
extern void delLhapdfList (lhapdfList * list);

/* interpolates data for given x and q    */
/* according to information stored in W.  */
/* result should be multiplied by the     */
/* power factors later on                 */
extern double lhapdfVal (double x, double q, int i);
extern double lhapdfValCPYTH (double x, double q, int i);
extern double TESTlhapdfValCPYTH (double x, double q, int i);

/* interpolates data for QCD-alpha(Q)     */
extern double lhapdf_interAlpha (double q);

extern void set_QCDLambda (int beam);
extern double lhapdf_QCDLambda(void);
extern void initLHAPDF (int beamnum, const char* setname, int mem, int prt);
#endif
