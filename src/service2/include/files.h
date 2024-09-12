/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __FILES_
#define __FILES_

typedef struct catrec
  {
    int nsub_;
    int status;
    int ndiagr_;
    long factpos;
    long rnumpos;
    long denompos;
  }
catrec;

extern FILE *menup;
extern FILE *menuq;
extern FILE *diagrp;		/* file of Adiagram; */
extern FILE *diagrq;		/* file of CSdiagram; */
extern FILE *catalog;

extern char *outputDir;
extern longstr pathtocomphep;
extern longstr pathtolhapdf;
extern longstr pathtouser;
extern longstr pathtoresults;
extern char mdFls[5][10];

#define MENUQ_NAME   scat("%stmp%cmenuq.ch",pathtouser,f_slash)
#define MENUP_NAME   scat("%stmp%cmenup.ch",pathtouser,f_slash)
#define DIAGRP_NAME  scat("%stmp%cproces.tp",pathtouser,f_slash)
#define DIAGRQ_NAME  scat("%stmp%ccsproces.tp",pathtouser,f_slash)
#define ARCHIV_NAME  scat("%stmp%carchive.bt",pathtouser,f_slash)
#define CATALOG_NAME scat("%stmp%ccatalog.tp",pathtouser,f_slash)


extern void wrt_menu (FILE * men, int menutype, int k,
		    char *txt, int ndel, int ncalc, int nrest, long recpos);
extern int rd_menu (FILE * men, int menutype, int k,
		char *txt, int *ndel, int *ncalc, int *nrest, long *recpos);

extern void copyfile (char *namefrom, char *nameto);
extern void nextFileName (char *f_name, char *firstname, char *ext);

#define  FREAD1(d,f)   fread(&(d),sizeof(d),1,f)
#define  FWRITE1(d,f)  f_write(&(d),sizeof(d),1,f)
#endif
