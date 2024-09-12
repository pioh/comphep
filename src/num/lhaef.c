/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lhaef.h"

static char * gen_desc_str  = NULL;
static char * file_str      = NULL;
static char * generator_str = NULL;
static char * model_str     = NULL;
static char * model_par_str = NULL;
static char * cuts_str      = NULL;
static char * authors_str   = NULL;
static char * process_str   = NULL;
static char * beam_str[2]   = {NULL, NULL};
static char * pdf_str[2]    = {NULL, NULL};
static char * pdf_qcd_str   = NULL;
static char * proc_qcd_str  = NULL;
static char * subproc_str   = NULL;

static int gen_desc_used  = 0;
static int file_used      = 0;
static int generator_used = 0;
static int model_used     = 0;
static int model_par_used = 0;
static int cuts_used      = 0;
static int authors_used   = 0;
static int process_used   = 0;
static int beam_used[2]   = {0, 0};
static int pdf_used[2]    = {0, 0};
static int pdf_qcd_used   = 0;
static int proc_qcd_used  = 0;
static int subproc_used   = 0;

void * init_lhaef_document (void)
{
  return 0;
}

int set_general_description (void * doc, lhaef_general_description * smpl)
{
  if (!gen_desc_used) {
    gen_desc_str = malloc ((
    1024 + 
    strlen (smpl->title) + 
    strlen (smpl->abstruct) + 
    strlen (smpl->comments) + 
    strlen (smpl->experiment->experiment) + 
    strlen (smpl->experiment->group) + 
    strlen (smpl->experiment->responsiblePerson) + 
    strlen (smpl->experiment->description)
    )* sizeof (char));

    sprintf (gen_desc_str, 
    "        <title>\%s</title>\n"                            
    "        <abstract>\%s</abstract>\n"                      
    "        <authorComments>\%s</authorComments>\n"          
    "        <experimentGroup>\n"                             
    "            <experiment>\%s</experiment>\n"              
    "            <group>\%s</group>\n"                        
    "            <responsiblePerson>\%s</responsiblePerson>\n"
    "            <description>\%s</description>\n"            
    "        </experimentGroup>\n"                            
    , 
    smpl->title, 
    smpl->abstruct, 
    smpl->comments, 
    smpl->experiment->experiment, 
    smpl->experiment->group, 
    smpl->experiment->responsiblePerson, 
    smpl->experiment->description);
    gen_desc_used = 1;
  }

  return 0;
}


int set_file_description (void * doc, lhaef_file * file)
{
  if (!file_used) {
    file_str = malloc ((
    1024 + 
    strlen (file->cross_section_unit) + 
    strlen (file->check_sum) + 
    strlen (file->check_sum_type)
    ) * sizeof (char));

    sprintf (file_str, 
    "            <eventsNumber>\%10i</eventsNumber>\n" 
    "            <crossSection unit=\"\%s\" "
    "errorMinus=\"\%11.4E\" errorPlus=\"\%11.4E\">"
    "\%11.4E</crossSection>\n"
    "            <fileSize>\%13i</fileSize>\n"
    "            <checksum type=\"\%s\">\%s</checksum>\n"
    "            <comments></comments>\n"
    "            <location>\n"
    "              <path/>\n"
    "            </location>\n"
    , 
    file->nevents, 
    file->cross_section_unit, 
    file->cross_section_lower_error, 
    file->cross_section_upper_error, 
    file->cross_section, 
    file->size, 
    file->check_sum_type, 
    file->check_sum);
    file_used = 1;
  }

  return 0;
}


char * get_file_description (void) 
{
  if (file_used) {
    return file_str;
  }
  return 0;
}

int set_generator (void * doc, lhaef_generator * gen)
{
  if (!generator_used) {
    generator_str = malloc ((
    512 + 
    strlen (gen->name) +
    strlen (gen->version) +
    strlen (gen->homepage) +
    strlen (gen->description)
    ) * sizeof (char));

    sprintf (generator_str,
    "            <name>\%s</name>\n"
    "            <version>\%s</version>\n"
    "            <homepage>\%s</homepage>\n"
    "            <description>\n"
    "\%s\n"
    "            </description>\n"
    , 
    gen->name, 
    gen->version, 
    gen->homepage, 
    gen->description );
    generator_used = 1;
  }

  return 0;
}


int set_model (void * doc, lhaef_model * model)
{
  if (!model_used) {
    model_str = malloc ((
    512 + 
    strlen (model->name) + 
    strlen (model->description)
    ) * sizeof (char));

    sprintf (model_str, 
    "             <name>%s</name>\n"
    "             <description>%s</description>\n"
    ,
    model->name, 
    model->description);
    model_used = 1;
  }

  return 0;
}


int add_model_parameter (void * doc, lhaef_parameter * par)
{
  int len;
  char * tmp_par = malloc ((
  512 + 
  strlen (par->name            ) +
  strlen (par->value           ) +
  strlen (par->notation->plain ) +
  strlen (par->notation->latex ) +
  strlen (par->description     )
  ) * sizeof (char));

  sprintf (tmp_par, 
  "                <parameter>\n"                       
  "                    <name>\%s</name>\n"              
  "                    <value>\%s</value>\n"            
  "                    <notation>\n"                    
  "                        <plain>\%s</plain>\n"        
  "                        <Latex>\%s</Latex>\n"        
  "                    </notation>\n"                   
  "                    <description>\%s</description>\n"
  "                </parameter>\n"                      
  ,
  par->name, 
  par->value, 
  par->notation->plain, 
  par->notation->latex, 
  par->description);

  len = strlen (tmp_par) + 1;
  if (NULL != model_par_str) {
    len += strlen (model_par_str);
    model_par_str = realloc (model_par_str, len * sizeof(char));
  } else {
    model_par_str = malloc (len * sizeof(char));
    model_par_str[0] = '\0';
    model_par_used = 1;
  }

  strcat (model_par_str, tmp_par);
  free(tmp_par);

  return 0;
}


int add_cutset (void * doc, int cutsetnumber)
{
  if (!cuts_used) {
    cuts_str = malloc (256 * sizeof (char));
    sprintf (cuts_str, "            <cutSet cutset_id=\"%i\">\n", cutsetnumber);
    cuts_used = 1;
  } else {
    char tmp[256];
    int len = 256 + strlen (cuts_str);
    cuts_str = realloc (cuts_str, len * sizeof (char));
    strcat (cuts_str , "            </cutSet>\n");
    sprintf (tmp, "            <cutSet cutset_id=\"%i\">\n", cutsetnumber);
    strcat (cuts_str , tmp);
  }

  return 0;
}

int add_cut (void * doc, lhaef_cut * cut, int minon, int maxon)
{
  char * tmp_cut;
  char frmt[1024];
  int len;

  frmt[0] = '\0';
  strcat (frmt, "                <cut>\n"                           );
  strcat (frmt, "                    <object>\n"                    );
  strcat (frmt, "                        <name>\%s</name>\n"        );
  strcat (frmt, "                        <notation>\n"              );
  strcat (frmt, "                            <plain>\%s</plain>\n"  );
  strcat (frmt, "                            <Latex>\%s</Latex>\n"  );
  strcat (frmt, "                        </notation>\n"             );
  strcat (frmt, "                    </object>\n"                   );
  if (minon) strcat (frmt, "                    <minValue>%11.4E</minValue>\n");
  if (maxon) strcat (frmt, "                    <maxValue>%11.4E</maxValue>\n");
  strcat (frmt, "                    <logic>\%s</logic>\n"          );
  strcat (frmt, "                </cut>\n"                          );

  tmp_cut = malloc ((
  1024 + 
  strlen (cut->name            ) +
  strlen (cut->cut_logic       ) +
  strlen (cut->notation->plain ) +
  strlen (cut->notation->latex )
  ) * sizeof (char));

  if (minon && maxon) {
    sprintf (tmp_cut, frmt, 
    cut->name, 
    cut->notation->plain, 
    cut->notation->latex, 
    cut->min_value, 
    cut->max_value, 
    cut->cut_logic 
    );
  }

  if (minon && !maxon) {
    sprintf (tmp_cut, frmt, 
    cut->name, 
    cut->notation->plain, 
    cut->notation->latex, 
    cut->min_value, 
    cut->cut_logic 
    );
  }

  if (!minon && maxon) {
    sprintf (tmp_cut, frmt, 
    cut->name, 
    cut->notation->plain, 
    cut->notation->latex, 
    cut->max_value, 
    cut->cut_logic 
    );
  }

  len = strlen (tmp_cut) + 1;
  if (NULL != cuts_str) {
    len += strlen (cuts_str);
    cuts_str = realloc (cuts_str, len * sizeof(char));
  } else {
    cuts_str = malloc (len * sizeof(char));
    cuts_str[0] = '\0';
    cuts_used = 1;
  }

  cuts_str = strcat (cuts_str, tmp_cut);
  free(tmp_cut);
  return 0;
}

int set_author_record (void * doc)
{
  return 0;
}


int add_author (void * doc, lhaef_author * author)
{
  int len = 512 + 
  strlen (author->firstName   ) +
  strlen (author->lastName    ) +
  strlen (author->email       ) +
  strlen (author->experiment  ) +
  strlen (author->group       ) +
  strlen (author->organization);
  char * tmp_author = malloc (len * sizeof (char));

  sprintf (tmp_author, 
  "            <author>\n"
  "                <firstName>\%s</firstName>\n"
  "                <lastName>\%s</lastName>\n"
  "                <email>\%s</email>\n"
  "                <experiment>\%s</experiment>\n"
  "                <group>\%s</group>\n"
  "                <organization>\%s</organization>\n"
  "            </author>\n"
  ,
  author->firstName, 
  author->lastName, 
  author->email, 
  author->experiment, 
  author->group, 
  author->organization 
  );

  len = strlen (tmp_author) + 2;
  if (NULL != authors_str) {
    len += strlen (authors_str);
    authors_str = realloc (authors_str, len * sizeof(char));
  } else {
    authors_str = malloc (len * sizeof(char));
    authors_str[0] = '\0';
    authors_used = 1;
  }

  authors_str = strcat (authors_str, tmp_author);
  free (tmp_author);

  return 0;
}


int create_process (void * doc, lhaef_process * proc)
{
  if (!process_used) {
    process_str = malloc ((
    1024 + 
    strlen (proc->final_state                 ) +
    strlen (proc->final_state_notation->plain ) +
    strlen (proc->final_state_notation->latex ) +
    strlen (proc->cross_section_unit          ) 
    ) * sizeof (char));

    sprintf (process_str, 
    "            <finalState>\n"
    "                <state>\%s</state>\n"
    "                <notation>\n"
    "                    <plain>\%s</plain>\n"
    "                    <Latex>\%s</Latex>\n"
    "                </notation>\n"
    "            </finalState>\n"
    "            <crossSection unit=\"\%s\">\%11.4E</crossSection>\n"
    ,
    proc->final_state, 
    proc->final_state_notation->plain, 
    proc->final_state_notation->latex, 
    proc->cross_section_unit, 
    proc->cross_section );
    process_used = 1;
  }

  return 0;
}


int add_process_beam (void * doc, int i, lhaef_beam * beam)
{
  if (!beam_used[i]) {
    beam_str[i] = malloc ((1024 + strlen (beam->particle_name)) * sizeof (char));

    sprintf (beam_str[i], 
    "                <particle pdgcode=\"%i\">%s</particle>\n"
    "                <energy unit=\"GeV\">%f</energy>\n"
    , 
    beam->particle_pdgcode, 
    beam->particle_name, 
    beam->energy);
    beam_used[i] = 1;
  }

  return 0;
}


int add_process_beam_pdf (void * doc, int i, lhaef_pdf * pdf)
{
  if (!pdf_used[i]) {
    pdf_str[i] = malloc ((
    1024 + 
    strlen (pdf->name) +
    strlen (pdf->version) +
    strlen (pdf->lhapdf_filename)
    ) * sizeof (char));

    sprintf (pdf_str[i], 
    "                <pdf name=\"%s\" version=\"%s\" "
    "PDFLIBset=\"%i\" PDFLIBgroup=\"%i\" LHAPDFset=\"%i\" "
    "LHAPDFmember=\"%i\" LHAPDFfile=\"%s\" />\n"
    , 
    pdf->name, 
    pdf->version, 
    pdf->pdflib_set, 
    pdf->pdflib_group, 
    pdf->lhapdf_set, 
    pdf->lhapdf_member, 
    pdf->lhapdf_filename );
    pdf_used[i] = 1;
  }

  return 0;
}


int add_process_beam_pdf_alphas (void * doc, lhaef_alphas * alphas)
{
  if (!pdf_qcd_used) {
    pdf_qcd_str = malloc ((1024 + strlen (alphas->description)) * sizeof (char));

    sprintf (pdf_qcd_str, 
    "                <QCDCoupling>\n"
    "                    <Lambda unit=\"GeV\">\%f</Lambda>\n"
    "                    <NFlavours>\%i</NFlavours>\n"
    "                    <NLoopsAlphaS>\%i</NLoopsAlphaS>\n"
    "                </QCDCoupling>\n"
    , 
    alphas->lambda, 
    alphas->nflavours, 
    alphas->nloops);
    pdf_qcd_used = 1;
  }

  return 0;
}


int add_process_alphas (void * doc, lhaef_alphas * alphas)
{
  if (!proc_qcd_used) {
    proc_qcd_str = malloc ((1024 + strlen (alphas->description)) * sizeof (char));

    sprintf (proc_qcd_str, 
    "            <QCDCoupling>\n"
    "                <Lambda unit=\"GeV\">\%f</Lambda>\n"
    "                <NFlavours>\%i</NFlavours>\n"
    "                <NLoopsAlphaS>\%i</NLoopsAlphaS>\n"
    "            </QCDCoupling>\n"
    , 
    alphas->lambda, 
    alphas->nflavours, 
    alphas->nloops);
    proc_qcd_used = 1;
  }

return 0;
}


int add_subprocess (void * doc, lhaef_subprocess * sub, int cutsetnumber)
{
  int len;
  char * tmp = malloc ((
  1024 + 
  strlen (sub->name                        ) +
  strlen (sub->cross_section_unit          ) +
  strlen (sub->FactorisationScale->plain   ) +
  strlen (sub->FactorisationScale->latex   ) +
  strlen (sub->RenormalisationScale->plain ) +
  strlen (sub->RenormalisationScale->latex ) 
  ) * sizeof (char));

  sprintf (tmp, 
  "                <subprocess cutset_id=\"\%i\">\n"
  "                    <notation>\%s</notation>\n"
  "                    <crossSection unit=\"\%s\" "
  "errorMinus=\"\%11.4E\" errorPlus=\"\%11.4E\">"
  "\%11.4E</crossSection>\n"
  "                    <FactorisationScale>\n"
  "                        <plain>\%s</plain>\n"
  "                        <Latex>\%s</Latex>\n"
  "                    </FactorisationScale>\n"
  "                    <RenormalisationScale>\n"
  "                        <plain>\%s</plain>\n"
  "                        <Latex>\%s</Latex>\n"
  "                    </RenormalisationScale>\n"
  "                </subprocess>\n"
  ,
  cutsetnumber,
  sub->name, 
  sub->cross_section_unit, 
  sub->cross_section_lower_error, 
  sub->cross_section_upper_error, 
  sub->cross_section, 
  sub->FactorisationScale->plain, 
  sub->FactorisationScale->latex, 
  sub->RenormalisationScale->plain, 
  sub->RenormalisationScale->latex 
  );

  len = strlen (tmp) + 1;
  if (NULL != subproc_str) {
    len += strlen (subproc_str);
    subproc_str = realloc (subproc_str, len * sizeof(char));
  } else {
    subproc_str = malloc (len * sizeof(char));
    subproc_str[0] = '\0';
    subproc_used = 1;
  }

  subproc_str = strcat (subproc_str, tmp);
  free(tmp);

  return 0;
}

char * form_hepml_document (void * doc)
{
  int i;
  char * the_hepml_doc = NULL;
  int len[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int llen = 1024;

  if (NULL != file_str)      len[ 0] = strlen (file_str);
  if (NULL != gen_desc_str)  len[ 1] = strlen (gen_desc_str);
  if (NULL != generator_str) len[ 2] = strlen (generator_str);
  if (NULL != model_str)     len[ 3] = strlen (model_str);
  if (NULL != model_par_str) len[ 4] = strlen (model_par_str);
  if (NULL != beam_str[0])   len[ 5] = strlen (beam_str[0]);
  if (NULL != pdf_str[0])    len[ 6] = strlen (pdf_str[0]);
  if (NULL != pdf_qcd_str)   len[ 7] = strlen (pdf_qcd_str);
  if (NULL != beam_str[1])   len[ 8] = strlen (beam_str[1]);
  if (NULL != pdf_str[1])    len[ 9] = strlen (pdf_str[1]);
  if (NULL != pdf_qcd_str)   len[10] = strlen (pdf_qcd_str);
  if (NULL != proc_qcd_str)  len[11] = strlen (proc_qcd_str);
  if (NULL != process_str)   len[12] = strlen (process_str);
  if (NULL != subproc_str)   len[13] = strlen (subproc_str);
  if (NULL != cuts_str)      len[14] = strlen (cuts_str);
  if (NULL != authors_str)   len[15] = strlen (authors_str);
  for (i = 0; i < 16; ++i) llen += len[i];

  the_hepml_doc = malloc (llen * sizeof (char));
  the_hepml_doc[0] = '\0';

  strcat (the_hepml_doc, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  strcat (the_hepml_doc, "<samples xmlns=\"http://mcdb.cern.ch/hepml/0.2/\"\n");
  strcat (the_hepml_doc, "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
  strcat (the_hepml_doc, "    xsi:schemaLocation=\"http://mcdb.cern.ch/hepml/0.2/ file:/home/dudko/MCDB/hepml/schemas/0.2/hepml.xsd\">\n");
  strcat (the_hepml_doc, "    <files>\n");
  strcat (the_hepml_doc, "        <file>\n");
  if (1 == file_used) strcat (the_hepml_doc, file_str);
  strcat (the_hepml_doc, "        </file>\n");
  strcat (the_hepml_doc, "    </files>\n");
  strcat (the_hepml_doc, "    <description>\n");
  if (1 == gen_desc_used) strcat (the_hepml_doc, gen_desc_str);
  strcat (the_hepml_doc, "        <generator>\n");
  if (1 == generator_used) strcat (the_hepml_doc, generator_str);
  strcat (the_hepml_doc, "        </generator>\n");
  strcat (the_hepml_doc, "        <process>\n");
  strcat (the_hepml_doc, "            <beam1>\n");
  if (1 == beam_used[0]) strcat (the_hepml_doc, beam_str[0]);
  if (1 == pdf_used[0])  strcat (the_hepml_doc, pdf_str[0]);
  if (1 == pdf_qcd_used) strcat (the_hepml_doc, pdf_qcd_str);
  strcat (the_hepml_doc, "            </beam1>\n");
  strcat (the_hepml_doc, "            <beam2>\n");
  if (1 == beam_used[1]) strcat (the_hepml_doc, beam_str[1]);
  if (1 == pdf_used[1])  strcat (the_hepml_doc, pdf_str[1]);
  if (1 == pdf_qcd_used) strcat (the_hepml_doc, pdf_qcd_str);
  strcat (the_hepml_doc, "            </beam2>\n");
  if (1 == proc_qcd_used) strcat (the_hepml_doc, proc_qcd_str);
  if (1 == process_used) strcat (the_hepml_doc, process_str);
  strcat (the_hepml_doc, "            <subprocesses>\n");
  if (1 == subproc_used) strcat (the_hepml_doc, subproc_str);
  strcat (the_hepml_doc, "            </subprocesses>\n");
  strcat (the_hepml_doc, "        </process>\n");
  strcat (the_hepml_doc, "        <model>\n");
  if (1 == model_used) strcat (the_hepml_doc, model_str);
  strcat (the_hepml_doc, "            <parameters>\n");
  if (1 == model_par_used) strcat (the_hepml_doc, model_par_str);
  strcat (the_hepml_doc, "            </parameters>\n");
  strcat (the_hepml_doc, "        </model>\n");
  if (1 == cuts_used) {
  strcat (the_hepml_doc, "        <cuts>\n");
    strcat (the_hepml_doc, cuts_str);
    strcat (the_hepml_doc, "            </cutSet>\n");
  strcat (the_hepml_doc, "        </cuts>\n");
  } else {
    strcat (the_hepml_doc, "        <cuts/>\n");
  }
  strcat (the_hepml_doc, "        <authors>\n");
  if (1 == authors_used) strcat (the_hepml_doc, authors_str);
  strcat (the_hepml_doc, "        </authors>\n");
  strcat (the_hepml_doc, "        <relatedPapers></relatedPapers>\n");
  strcat (the_hepml_doc, "    </description>\n");
  strcat (the_hepml_doc, "</samples>\n");

  if (NULL != file_str     ) free (file_str     );
  if (NULL != gen_desc_str ) free (gen_desc_str );
  if (NULL != generator_str) free (generator_str);
  if (NULL != model_str    ) free (model_str    );
  if (NULL != model_par_str) free (model_par_str);
  if (NULL != beam_str[0]  ) free (beam_str[0]  );
  if (NULL != pdf_str[0]   ) free (pdf_str[0]   );
  if (NULL != beam_str[1]  ) free (beam_str[1]  );
  if (NULL != pdf_str[1]   ) free (pdf_str[1]   );
  if (NULL != pdf_qcd_str  ) free (pdf_qcd_str  );
  if (NULL != proc_qcd_str ) free (proc_qcd_str );
  if (NULL != process_str  ) free (process_str  );
  if (NULL != subproc_str  ) free (subproc_str  );
  if (NULL != cuts_str     ) free (cuts_str     );
  if (NULL != authors_str  ) free (authors_str  );

  file_str      = NULL;
  gen_desc_str  = NULL;
  generator_str = NULL;
  model_str     = NULL;
  model_par_str = NULL;
  beam_str[0]   = NULL;
  pdf_str[0]    = NULL;
  beam_str[1]   = NULL;
  pdf_str[1]    = NULL;
  pdf_qcd_str   = NULL;
  proc_qcd_str  = NULL;
  process_str   = NULL;
  subproc_str   = NULL;
  cuts_str      = NULL;
  authors_str   = NULL;

  file_used      = 0;
  gen_desc_used  = 0;
  generator_used = 0;
  model_used     = 0;
  model_par_used = 0;
  beam_used[0]   = 0;
  pdf_used[0]    = 0;
  beam_used[1]   = 0;
  pdf_used[1]    = 0;
  pdf_qcd_used   = 0;
  proc_qcd_used  = 0;
  process_used   = 0;
  subproc_used   = 0;
  cuts_used      = 0;
  authors_used   = 0;

  return the_hepml_doc;
}
