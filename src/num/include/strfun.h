/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __STRFUN__
#define __STRFUN__

extern int pdf_menu (int ibeam);
extern int beam_menu (int ibeam);
extern double strfun_ (int factr, double x, double y, double q);        /* #-mdl */

extern int rd_sf__ (FILE * mode);
extern int wrt_sf__ (FILE * mode);

extern int initStrFun (char p_name1[], char p_name2[]);

extern void strFunName (int i, char * beam, char * pdf);

typedef struct Str_fun_Info
  {
    char pdf_name[20];          /* Name of structure function */
    char version[20];           /* Version of structure function */
    char prt_name[20];          /* Name of beam particle in structure function */
    double prt_mass;            /* Mass of beam parlicle in structure function */
    int PDFLIBset;
    int PDFLIBgroup;
    int LHAPDFset;
    int LHAPDFmember;
    int N_extra_commands;
   char extra_commands[50][50];        /* Extra info about structure function */
  }
Str_fun_Info;


extern void wrt_sf_NF_ (int i, Str_fun_Info * info);
#endif
