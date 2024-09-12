/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __LHAEF__
#define __LHAEF__

#define NAMELENGTH 128
#define TITLLENGTH 1024
#define DESCLENGTH 4096

typedef struct general_type
  {
    char * name;
  }
general_type;

/* defines all necessary notations for visualization */
typedef struct lhaef_math_notation
  {
    char plain[TITLLENGTH];         // in plain texts
    char html[TITLLENGTH];          // in html documents
    char latex[NAMELENGTH];         // in LaTeX documents
    char mathml[TITLLENGTH];       // in MathML documents
  }
lhaef_math_notation;

/* defines an elementary cut */
/* cut_logic means whether the interval (min_value, max_value) should be excluded or kept */
typedef struct lhaef_cut
{
    char name[NAMELENGTH];
    lhaef_math_notation * notation;
    double min_value;
    double max_value;
    char cut_logic[NAMELENGTH];
  }
lhaef_cut;

/* defines formal information about an event files some general parameters */
typedef struct lhaef_file
  {
    char location_path[TITLLENGTH];
    char location_type[TITLLENGTH];
    char machine_name[TITLLENGTH];
    int nevents;                            // the numner of event in the file
    int size;                               // the file size (in bytes)
    double cross_section;                   // the total cross section
    double cross_section_upper_error;
    double cross_section_lower_error;
    char cross_section_unit[NAMELENGTH];    // this parameter must be equal to a value of cross_section_type
    char check_sum[128];                    // check sum parameter
    char check_sum_type[TITLLENGTH];
  }
lhaef_file;

/* a short description of the experiment the event has been prepared to */
typedef struct lhapdf_experiment
  {
    char experiment[TITLLENGTH];
    char group[TITLLENGTH];
    char responsiblePerson[TITLLENGTH];
    char description[DESCLENGTH];
  }
lhapdf_experiment;

/* the most general information on the sample, required by MCDB and can be omitted */
typedef struct lhaef_general_description
  {
    char title[TITLLENGTH];
    char abstruct[DESCLENGTH];
    char comments[DESCLENGTH];
    lhapdf_experiment * experiment;
  }
lhaef_general_description;

/* define parton density functions, pdf_set/pdf_group and lhapdf_set/lhapdf_member can't be define at once */
typedef struct lhaef_pdf
  {
    char name[TITLLENGTH];
    char version[TITLLENGTH];
    int pdflib_set;
    int pdflib_group;
    int lhapdf_set;
    int lhapdf_member;
    char lhapdf_filename[TITLLENGTH];
  }
lhaef_pdf;

/* name and short description of the physical model applied */
typedef struct lhaef_model
  {
    char name[TITLLENGTH];            /* Name of the model */
    char description[DESCLENGTH];      /* Short description of the model */
  }
lhaef_model;

/* description of model parameters */
typedef struct lhaef_parameter
  {
    char name[TITLLENGTH];            /* Name of the parameter */
    char value[TITLLENGTH];           /* string with the parameter value */
    char description[DESCLENGTH];      /* Short description of the parameter */
    lhaef_math_notation * notation;
  }
lhaef_parameter;

/* QCD Lambda based definition of the QCD coupling, can be included into PDF or/and process */
typedef struct lhaef_alphas
  {
    int nloops;
    int nflavours;
    double lambda;
    char lambda_unit[NAMELENGTH];
    char description[DESCLENGTH];
  }
lhaef_alphas;

/* definion of beam, a beam can be alone if a fixed target is considered */
typedef struct lhaef_beam
  {
    char particle_name[NAMELENGTH];
    int particle_pdgcode;
    double energy;
    char energy_unit[NAMELENGTH];
  }
lhaef_beam;

/* general information on process. At least one subprocess should associated with the process */
typedef struct lhaef_process
  {
    char final_state[NAMELENGTH];
    lhaef_math_notation * final_state_notation;
    double cross_section;
    double cross_section_upper_error;
    double cross_section_lower_error;
    char cross_section_unit[NAMELENGTH];
  }
lhaef_process;

/* some information on author */
typedef struct lhaef_author
  {
    char firstName[NAMELENGTH];
    char lastName[NAMELENGTH];
    char email[TITLLENGTH];
    char experiment[TITLLENGTH];
    char group[TITLLENGTH];
    char organization[TITLLENGTH];
  }
lhaef_author;

/* description of the Monte-Carlo generator */
typedef struct lhaef_generator
  {
    char name[NAMELENGTH];
    char version[NAMELENGTH];
    char homepage[TITLLENGTH];
    char description[DESCLENGTH];
  }
lhaef_generator;

/* information on a subprocess: name, the total cross section, factorisation and renormalisation scales */
typedef struct lhaef_subprocess
  {
    char name[TITLLENGTH];
    double cross_section;
    double cross_section_upper_error;
    double cross_section_lower_error;
    char cross_section_unit[NAMELENGTH];
    lhaef_math_notation * FactorisationScale;
    lhaef_math_notation * RenormalisationScale;
  }
lhaef_subprocess;

/***********************************************************************************************/

void * init_lhaef_document (void);

int set_general_description (void * doc, lhaef_general_description * sample);
int set_file_description (void * doc, lhaef_file * file);
char * get_file_description (void);

int set_generator (void * doc, lhaef_generator * gen);

int set_model (void * doc, lhaef_model * model);
int add_model_parameter (void * doc, lhaef_parameter * par);

int add_cutset (void * doc, int cutsetnumber);
int add_cut (void * doc, lhaef_cut * cut, int minon, int maxon);

int set_author_record (void * doc);
int add_author (void * doc, lhaef_author * author);

int create_process (void * doc, lhaef_process * proc);
int add_process_beam (void * doc, int i, lhaef_beam * beam);
int add_process_beam_pdf (void * doc, int i, lhaef_pdf * pdf);
int add_process_beam_pdf_alphas (void * doc, lhaef_alphas * alphas);
int add_process_alphas (void * doc, lhaef_alphas * alphas);

int add_subprocess (void * doc, lhaef_subprocess * subproc, int cutsetnumber);

char * form_hepml_document (void * doc);

#endif
