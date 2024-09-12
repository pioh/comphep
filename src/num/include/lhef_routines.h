/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __MIX_LHEFF_ROUTINES__
#define __MIX_LHEFF_ROUTINES__

#define MAX_FILE_EVENT 2048

#ifdef LIBXML
  extern int formXMLtree (const char fname[], int i);
  extern int write_file_header_libxml2 (const char fname[], const char mode[]);
  extern char * prepare_hepml_header_libxml2_dynamic (void);
#else
  extern int analyzeLHEfile (const char fname[], int i);
#endif
extern int analyzeMultiProcLHEfile (const char fname[]);

int setInfoWithoutHEPML (char fname[], FILE * f, int procnum);

extern int set_final_numbers (int n, int size, double cs, double cserr, char check[128]);
extern void set_cs (double cs, double cserr);

extern long getEventPosition (int samplenum);
extern int getEventNumber (int samplenum);
extern double getCrossSection (int samplenum);
extern double getCrossSectionErr (int samplenum);
extern char * getFinalStateName (void);
extern char * getSubprocessName (void);

extern int getTotProcNumber (void);
extern int getPbeam (int i);
extern double getEbeam (int i);
extern int getPDFLIBgroup (int i);
extern int getPDFLIBset (int i);
extern char * getCHEPversion (void);
#endif
