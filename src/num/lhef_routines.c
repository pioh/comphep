/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
*------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#ifdef LIBXML
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#include "service2/include/kfcodes.h"

#include "LesHouches.h"
#include "lhaef.h"
#include "tag_reader.h"
#include "tag_parser.h"
#include "event_reader.h"
#include "lhef_routines.h"

static int nevents[MAX_FILE_EVENT];
static double cs[MAX_FILE_EVENT];
static double cserr[MAX_FILE_EVENT];
static long nposition[MAX_FILE_EVENT];
static char final_state_name[1024];
static char subproc_name[1024];

static int tot_proc_num;
static int the_pbeam[2];
static double the_ebeam[2];
static int the_PDFLIBgroup[2];
static int the_PDFLIBset[2];
static char the_version[256];

static double tot_cs = 0.0;
static double tot_cserr = 0.0;
void set_cs (double c1, double c2) {
  tot_cs = c1;
  tot_cserr = c2;
  return;
}

int set_final_numbers (int n, int size, double cs1, double cserr1, char check[128]) {
  lhaef_file lhafile;
  void * hepml_header = init_lhaef_document ();

  tot_cs = cs1;
  tot_cserr = cserr1;
  lhafile.nevents = n;
  lhafile.size = size;
  lhafile.cross_section = cs1;
  lhafile.cross_section_upper_error = cserr1;
  lhafile.cross_section_lower_error = cserr1;
  strcpy (lhafile.cross_section_unit, "pb");
  strcpy (lhafile.check_sum, check);
  strcpy (lhafile.check_sum_type, "md5");

  return set_file_description (hepml_header, &lhafile);;
}

#ifdef LIBXML
static int fill_base_hepml = 1;
static const xmlChar * nameStr        = (const xmlChar *)"name";
static const xmlChar * valueStr       = (const xmlChar *)"value";
static const xmlChar * plainStr       = (const xmlChar *)"plain";
static const xmlChar * latexStr       = (const xmlChar *)"latex";
static const xmlChar * symbStr        = (const xmlChar *)"notation";
static const xmlChar * subproclistStr = (const xmlChar *)"subprocesses";
static const xmlChar * subprocStr     = (const xmlChar *)"subprocess";
static const xmlChar * cutSetStr      = (const xmlChar *)"cutSet";
static const xmlChar * descriptionStr = (const xmlChar *)"description";
static const xmlChar * processStr     = (const xmlChar *)"process";
static const xmlChar * filelistStr    = (const xmlChar *)"files";
static const xmlChar * fileStr        = (const xmlChar *)"file";
static const xmlChar * cutlistStr     = (const xmlChar *)"cuts";
static const xmlChar * cutStr         = (const xmlChar *)"cut";
static const xmlChar * modelStr       = (const xmlChar *)"model";
static const xmlChar * parlistStr     = (const xmlChar *)"parameters";
static const xmlChar * parStr         = (const xmlChar *)"parameter";
static const xmlChar * stateStr       = (const xmlChar *)"state";
static const xmlChar * fstateStr      = (const xmlChar *)"finalState";
static const xmlChar * errorPlusStr   = (const xmlChar *)"errorPlus";
static const xmlChar * errorMinusStr  = (const xmlChar *)"errorMinus";
static const xmlChar * cutsetidStr = (const xmlChar *)"cutset_id";
static const xmlChar * beamStr[2] = {(const xmlChar *)"beam1", (const xmlChar *)"beam2"};

static int prepare_hepml_header_libxml2_static (const xmlNodePtr desc);

static xmlNodePtr mdlparlist = NULL;
static xmlNodePtr sublist = NULL;
static xmlNodePtr cutlist = NULL;

static xmlDocPtr parse_doc (char fname[]) {
  xmlDocPtr doc = NULL;
  xmlParserCtxtPtr parser = xmlNewParserCtxt ();
  if (NULL == parser) {
    fprintf (stderr, "mix: (error) failed to allocate libxml2 parser for %s\n", fname);
  } else {
    doc = xmlCtxtReadFile (parser, fname, NULL, XML_PARSE_PEDANTIC);
    xmlFreeParserCtxt (parser);
  }
  return doc;
}


static int get_son_number (xmlNodePtr parent, const xmlChar * name) {
  int num = 0;
  xmlNodePtr son = parent->children;
  while (son) {
    if (!xmlStrcmp(son->name, name)) {
      ++num;
    }
    son = son->next;
  }
  return num;
}


static xmlNodePtr get_first_son (xmlNodePtr parent, const xmlChar * name) {
  if (!parent) return NULL;

  xmlNodePtr son = parent->children;
  while (son) {
    if (!xmlStrcmp (son->name, name)) {
      return xmlCopyNode(son, 1);
    }
    son = son->next;
  }
  return NULL;
}


static xmlNodePtr get_first_grandson (xmlNodePtr parent, const xmlChar * sonname, const xmlChar * grandsonname) {
  if (!parent) return NULL;

  xmlNodePtr son = parent->children;
  while (son) {
    if (!xmlStrcmp (son->name, sonname)) {
      xmlNodePtr grandson = son->children;
      while (grandson) {
        if (!xmlStrcmp (grandson->name, grandsonname)) {
          return xmlCopyNode(grandson, 1);
        }
        grandson = grandson->next;
      }
    }
    son = son->next;
  }
  return NULL;
}


static int modelMerger (const xmlNodePtr pars) {
  xmlNodePtr son = NULL;
  xmlNodePtr mdldiff = xmlNewNode (NULL, (const xmlChar *)"model_difference");

  if (!pars) return -1;
  if (!mdlparlist) {
    mdlparlist = xmlCopyNode (pars, 1);
    return 0;
  }

  son = pars->children;
  while (son) {
    if (!xmlStrcmp(son->name, parStr)) {
      xmlNodePtr t1 = xmlCopyNode (get_first_son (son, nameStr), 1);
      xmlNodePtr t2 = xmlCopyNode (get_first_son (son, valueStr), 1);
      if (t1 && t2) {
        int unused_par = 1;
        const xmlChar * name = xmlNodeGetContent (t1);
        const xmlChar * valu = xmlNodeGetContent (t2);
        xmlNodePtr tmpson = mdlparlist->children;
        while (tmpson) {
          if (!xmlStrcmp(tmpson->name, parStr)) {
            xmlNodePtr t3 = xmlCopyNode (get_first_son (tmpson, nameStr), 1);
            xmlNodePtr t4 = xmlCopyNode (get_first_son (tmpson, valueStr), 1);
            if (t1 && t2) {
              const xmlChar * tmpname = xmlNodeGetContent (t3);
              const xmlChar * tmpvalu = xmlNodeGetContent (t4);
              if (!strcmp((char *)name, (char *)tmpname)) {
                if (!strcmp((char *)valu, (char *)tmpvalu)) {
                  unused_par = 0;
                  break;
                }
              }
            }
          }
          tmpson = tmpson->next;
        }
        if (unused_par) {
          if (!xmlAddChild (mdldiff, son)) return -2;
        }
      }
    }
    son = son->next;
  }

  son = mdldiff->children;
  while (son) {
    if (!xmlAddChild (mdlparlist, son)) return -3;
    son = son->next;
  }

  return 0;
}

static int curnum = 0;
static int xml_add_subprocess (const xmlNodePtr subprocesses, const xmlNodePtr cuts) {
  if (!subprocesses || !cuts) return -1;

  if (!sublist) {
    sublist = xmlNewNode (NULL, (const xmlChar *)"internal_subprocesses");
    curnum = 0;
  }

  if (!cutlist) {
    cutlist = xmlNewNode (NULL, (const xmlChar *)"internal_cuts");
  }
  if (!sublist || !cutlist) return -2;

//////////////////////////////////////////////////////////////////////////////////////////////////
  {
    char numStr[8];
    xmlNodePtr son = NULL;
    int subnum = get_son_number (subprocesses, subprocStr);
    int cutnum = get_son_number (cuts, cutSetStr);
    if (1 != subnum || 1 != cutnum) return -3;

    son = subprocesses->children;
    while (son) {
      if (!xmlStrcmp (son->name, subprocStr)) {
        xmlChar * numstr = xmlGetProp (son, cutsetidStr);
        if (numstr) {
          int n;
          if (1 == sscanf ((char *)numstr, "%i", &n)) {
            sprintf (numStr, "%i", n + curnum);
            xmlSetProp (son, cutsetidStr, (const xmlChar *)numStr);
            if (!xmlAddChild (sublist, son)) return -4;
            ++curnum;
          }
        }
      }
      son = son->next;
    }

    son = cuts->children;
    while (son) {
      if (!xmlStrcmp (son->name, cutSetStr)) {
        if (xmlGetProp (son, cutsetidStr)) {
          xmlSetProp (son, cutsetidStr, (const xmlChar *)numStr);
          if (!xmlAddChild (cutlist, son)) return -5;
        } else {
          return -6;
        }
      }
      son = son->next;
    }
  }
  return 0;
}


static long formXMLdocument (const char fname[], const char xml[]) {
  long pos;
  char buff[2048];
  char * event = NULL;
  FILE * source = fopen (fname, "r");
  FILE * target = fopen (xml, "w");

  if (!source || !target) return -1;
  fgets (buff, 2048, source);
  if (!strstr (buff, "<LesHouchesEvents version=\"1.0\">")) return -2;

  while (1) {
    pos = ftell (source);
    fgets (buff, 2048, source);
    char * tmp = strstr (buff, "<header>");
    if (tmp) {
      break;
    }
  }
  fseek (source, pos, SEEK_SET);

  while (1) {
    char c1 = fgetc(source);
    if (c1 == '<') {
      char c2 = fgetc(source);
      if (c2 == '?') {
        fgets (buff, 4, source);
        if (!strcmp (buff, "xml")) {
          pos = ftell (source);
          fseek (source, pos - 5, SEEK_SET);
          break;
        } else {
          pos = ftell (source);
          fseek (source, pos - 4, SEEK_SET);
        }
      } else {
        pos = ftell (source);
        fseek (source, pos - 1, SEEK_SET);
      }
    }
    if (c1 == EOF) return -3;
  }

  while (1) {
    pos = ftell (source);
    fgets (buff, 2048, source);
    char * tmp = strstr (buff, "</samples>");
    if (tmp) {
      fputs (tmp, target);
      break;
    }
    fputs (buff, target);
  }

  while (!event) {
    pos = ftell (source);
    fgets (buff, 2048, source);
    if (feof (source)) return -4;
    event = strstr (buff, "<event>");
  }
  fclose (target);
  fclose (source);

  return pos;
}


int formXMLtree (const char fname[], int procnum) {
  char * xml = malloc ((strlen(fname) + 5) * sizeof (char));
  int n;

  if (procnum > MAX_FILE_EVENT) {
    fprintf (stderr, "mix: (error) too many event files! Increase MAX_FILE_EVENT in src/num/include/lhef_routines.h!\n");
    exit (-1);
  }

  sprintf (xml, "%s.xml", fname);
  n = formXMLdocument (fname, xml);
  switch (n) {
    case -1:
      fprintf (stderr, "lhe module (warning): can not open file %s or %s\n", fname, xml);
      goto error;
    case -2:
      fprintf (stderr, "lhe module (warning): file %s is not LesHouchesEvents file\n", fname);
      goto error;
    case -3:
      fprintf (stderr, "lhe module (warning): no XML block in file %s\n", fname);
      goto error;
    case -4:
     fprintf (stderr, "lhe module (warning): no events in file %s\n", fname);
      goto error;
  }
  nposition[procnum] = n;

  xmlDocPtr doc = parse_doc (xml);
  if (!doc) {
    fprintf (stderr, "lhe module (error): failed to parse XML block in %s\n", xml);
      goto error;
  }

  {
    xmlNodePtr rnode = xmlDocGetRootElement (doc);
    xmlNodePtr xml_proc   = xmlCopyNode (get_first_grandson (rnode, descriptionStr, processStr), 1);
    xmlNodePtr xml_list   = xmlCopyNode (get_first_son (xml_proc, subproclistStr), 1);
    xmlNodePtr xml_fstate = xmlCopyNode (get_first_grandson (xml_proc, fstateStr, stateStr), 1);
    xmlNodePtr xml_subprc = xmlCopyNode (get_first_grandson (xml_proc, subproclistStr, subprocStr), 1);
    xmlNodePtr xml_sbname = xmlCopyNode (get_first_son (xml_subprc, symbStr), 1);
    xmlNodePtr xml_cuts   = xmlCopyNode (get_first_grandson (rnode, descriptionStr, cutlistStr), 1);
    xmlNodePtr xml_files  = xmlCopyNode (get_first_son (rnode, filelistStr), 1);
    xmlNodePtr xml_model  = xmlCopyNode (get_first_grandson (rnode, descriptionStr, modelStr), 1);
    xmlNodePtr xml_pars   = xmlCopyNode (get_first_son (xml_model, parlistStr), 1);
    xmlNodePtr xml_event  = xmlCopyNode (get_first_grandson (xml_files, fileStr, (const xmlChar *)"eventsNumber"), 1);
    xmlNodePtr xml_crosc  = xmlCopyNode (get_first_grandson (xml_files, fileStr, (const xmlChar *)"crossSection"), 1);

    if (fill_base_hepml) {
      xmlNodePtr xml_tmp = xmlCopyNode (get_first_son (rnode, descriptionStr), 1);
      prepare_hepml_header_libxml2_static (xml_tmp);
      fill_base_hepml = 0;
    }

    if (!xml_list || !xml_cuts || !xml_event || !xml_crosc || !xml_pars) {
      fprintf (stderr, "mix: (error) failed to find tags \"subrocesses\", \"cuts\", \"files\" in %s. Remove the file\n", xml);
      goto error;
    }
    if (1 != sscanf ((char *)xmlNodeGetContent (xml_event), "%i", &(nevents[procnum]))) {
        nevents[procnum] = -1;
      }
    if (1 != sscanf ((char *)xmlNodeGetContent (xml_crosc), "%le", &(cs[procnum]))) {
      cs[procnum] = -1.0;
    }
    if (1 != sscanf ((char *)xmlGetProp (xml_crosc, errorPlusStr), "%le", &(cserr[procnum]))) {
      cserr[procnum] = -1.0;
    }
    if (xml_fstate) {
      strcpy (final_state_name, (char *)xmlNodeGetContent (xml_fstate));
    }
    if (xml_sbname) {
      strcpy (subproc_name, (char *)xmlNodeGetContent (xml_sbname));
    }
    xml_add_subprocess (xml_list, xml_cuts);
    modelMerger (xml_pars);
  }

  xmlFreeDoc (doc);
  unlink (xml);
  free (xml);

  {
    char buff[2048];
    FILE * source = fopen (fname, "r");
    while (1) {
      if (feof (source)) {
        fclose (source);
        return -4;
      }
      fgets (buff, 2048, source);
      if (strstr (buff, "<init>")) {
        break;
      }
    }

    tot_proc_num = 0;
    fgets (buff, 1024, source);
    if (1) {
      int i1, i2, i3, i4, i5, i6, i7;
      double d1, d2;
      if (9 == sscanf (buff, "%d %d %le %le %d %d %d %d 3 %d", &i1, &i2, &d1, &d2, &i3, &i4, &i5, &i6, &i7)) {
        tot_proc_num = i7;
      }
    }
    fclose (source);
  }

  return 0;

error:
  free (xml);
  return -1;
}


static int prepare_hepml_header_libxml2_static (const xmlNodePtr desc) {
  int err = 0;
  void * hepml_header = init_lhaef_document ();
  xmlNodePtr xml_notation_plain;
  xmlNodePtr xml_notation_latex;

////////////////////////////////////////////////////////////////////
  lhaef_file lhafile;
  lhafile.nevents = 0;
  lhafile.size = 0;
  lhafile.cross_section = 0.0;
  lhafile.cross_section_upper_error = 0.0;
  lhafile.cross_section_lower_error = 0.0;
  strcpy (lhafile.cross_section_unit, "pb");
  strcpy (lhafile.check_sum, "00000000000000000000000000000000");  // 32-bit md5 code!!!!!
  strcpy (lhafile.check_sum_type, "md5");
  err = set_file_description (hepml_header, &lhafile);

////////////////////////////////////////////////////////////////////
  {
    lhaef_general_description gendesc;
    lhapdf_experiment exp;
    xmlNodePtr xml_title    = xmlCopyNode (get_first_son (desc, (const xmlChar *)"title"), 1);
    xmlNodePtr xml_abstruct = xmlCopyNode (get_first_son (desc, (const xmlChar *)"abstract"), 1);
    xmlNodePtr xml_comments = xmlCopyNode (get_first_son (desc, (const xmlChar *)"authorComments"), 1);
    xmlNodePtr xml_ExpGroup = xmlCopyNode (get_first_son (desc, (const xmlChar *)"experimentGroup"), 1);
    xmlNodePtr xml_ExpGroup_exp    = xmlCopyNode (get_first_son (xml_ExpGroup, (const xmlChar *)"experiment"), 1);
    xmlNodePtr xml_ExpGroup_group  = xmlCopyNode (get_first_son (xml_ExpGroup, (const xmlChar *)"group"), 1);
    xmlNodePtr xml_ExpGroup_person = xmlCopyNode (get_first_son (xml_ExpGroup, (const xmlChar *)"responsiblePerson"), 1);
    xmlNodePtr xml_ExpGroup_desc   = xmlCopyNode (get_first_son (xml_ExpGroup, descriptionStr), 1);

    gendesc.title[0] = 0;
    gendesc.abstruct[0] = 0;
    gendesc.comments[0] = 0;
    exp.experiment[0] = 0;
    exp.group[0] = 0;
    exp.responsiblePerson[0] = 0;
    exp.description[0] = 0;

    if (xml_title)           strcpy (gendesc.title,         xmlNodeGetContent (xml_title));
    if (xml_abstruct)        strcpy (gendesc.abstruct,      xmlNodeGetContent (xml_abstruct));
    if (xml_comments)        strcpy (gendesc.comments,      xmlNodeGetContent (xml_comments));
    if (xml_ExpGroup_exp)    strcpy (exp.experiment,        xmlNodeGetContent (xml_ExpGroup_exp));
    if (xml_ExpGroup_group)  strcpy (exp.group,             xmlNodeGetContent (xml_ExpGroup_group));
    if (xml_ExpGroup_person) strcpy (exp.responsiblePerson, xmlNodeGetContent (xml_ExpGroup_person));
    if (xml_ExpGroup_desc)   strcpy (exp.description,       xmlNodeGetContent (xml_ExpGroup_desc));

    gendesc.experiment = &exp;
    err *= set_general_description (hepml_header, &gendesc);
  }

////////////////////////////////////////////////////////////////////
  {
    lhaef_generator gen;
    xmlNodePtr xml_gen = xmlCopyNode (get_first_son (desc, (const xmlChar *)"generator"), 1);
    xmlNodePtr xml_gen_name     = xmlCopyNode (get_first_son (xml_gen, nameStr), 1);
    xmlNodePtr xml_gen_version  = xmlCopyNode (get_first_son (xml_gen, (const xmlChar *)"version"), 1);
    xmlNodePtr xml_gen_homepage = xmlCopyNode (get_first_son (xml_gen, (const xmlChar *)"homepage"), 1);
    xmlNodePtr xml_gen_desc     = xmlCopyNode (get_first_son (xml_gen, descriptionStr), 1);

    gen.name[0] = 0;
    the_version[0] = 0;
    gen.version[0] = 0;
    gen.homepage[0] = 0;
    gen.description[0] = 0;

    if (xml_gen_name)     strcpy (gen.name,        xmlNodeGetContent (xml_gen_name));
    if (xml_gen_version)  strcpy (the_version,     xmlNodeGetContent (xml_gen_version));
    if (xml_gen_version)  strcpy (gen.version,     the_version);
    if (xml_gen_homepage) strcpy (gen.homepage,    xmlNodeGetContent (xml_gen_homepage));
    if (xml_gen_desc)     {
      xmlChar * text = xmlNodeGetContent (xml_gen_desc);
      strcpy (gen.description, text);
    }
    err *= set_generator (hepml_header, &gen);
  }

////////////////////////////////////////////////////////////////////
  {
    xmlNodePtr xml_model = xmlCopyNode (get_first_son (desc, modelStr), 1);
    xmlNodePtr xml_model_name = xmlCopyNode (get_first_son (xml_model, nameStr), 1);
    xmlNodePtr xml_model_desc = xmlCopyNode (get_first_son (xml_model, descriptionStr), 1);
    lhaef_model model;
  
    model.name[0] = 0;
    model.description[0] = 0;
    if (xml_model_name) strcpy (model.name, xmlNodeGetContent (xml_model_name));
    if (xml_model_desc) strcpy (model.description, xmlNodeGetContent (xml_model_desc));
    err *= set_model (hepml_header, &model);
  }

////////////////////////////////////////////////////////////////////
  {
    xmlNodePtr xml_authors = xmlCopyNode (get_first_son (desc, (const xmlChar *)"authors"), 1);
    xmlNodePtr son = xml_authors->children;
    err = set_author_record (hepml_header);
    while (son) {
      if (!xmlStrcmp (son->name, (const xmlChar *)"author")) {
        lhaef_author author;
        xmlNodePtr xml_author_fname = xmlCopyNode (get_first_son (son, (const xmlChar *)"firstName"), 1);
        xmlNodePtr xml_author_lname = xmlCopyNode (get_first_son (son, (const xmlChar *)"lastName"), 1);
        xmlNodePtr xml_author_email = xmlCopyNode (get_first_son (son, (const xmlChar *)"email"), 1);
        xmlNodePtr xml_author_exprm = xmlCopyNode (get_first_son (son, (const xmlChar *)"experiment"), 1);
        xmlNodePtr xml_author_group = xmlCopyNode (get_first_son (son, (const xmlChar *)"group"), 1);
        xmlNodePtr xml_author_organ = xmlCopyNode (get_first_son (son, (const xmlChar *)"organization"), 1);

        author.firstName[0]    = 0;
        author.lastName[0]     = 0;
        author.email[0]        = 0;
        author.experiment[0]   = 0;
        author.group[0]        = 0;
        author.organization[0] = 0;

        if (xml_author_fname) strcpy (author.firstName,    xmlNodeGetContent (xml_author_fname));
        if (xml_author_lname) strcpy (author.lastName,     xmlNodeGetContent (xml_author_lname));
        if (xml_author_email) strcpy (author.email,        xmlNodeGetContent (xml_author_email));
        if (xml_author_exprm) strcpy (author.experiment,   xmlNodeGetContent (xml_author_exprm));
        if (xml_author_group) strcpy (author.group,        xmlNodeGetContent (xml_author_group));
        if (xml_author_organ) strcpy (author.organization, xmlNodeGetContent (xml_author_organ));
        err = add_author(hepml_header, &author);
      }
      son = son->next;
    }
  }

////////////////////////////////////////////////////////////////////
  {
    lhaef_process proc;
    lhaef_math_notation fs_notation;
    xmlNodePtr xml_process     = xmlCopyNode (get_first_son (desc, processStr), 1);
    xmlNodePtr xml_cs          = xmlCopyNode (get_first_son (xml_process, (const xmlChar *)"crossSection"), 1);
    xmlNodePtr xml_final       = xmlCopyNode (get_first_son (xml_process, (const xmlChar *)"finalState"), 1);
    xmlNodePtr xml_final_state = xmlCopyNode (get_first_son (xml_final, (const xmlChar *)"state"), 1);
    xml_notation_plain         = xmlCopyNode (get_first_grandson (xml_final, symbStr, plainStr), 1);
    xml_notation_latex         = xmlCopyNode (get_first_grandson (xml_final, symbStr, latexStr), 1);

    proc.final_state[0] = 0;
    fs_notation.plain[0] = 0;
    fs_notation.latex[0] = 0;

    if (xml_final_state) strcpy (proc.final_state, xmlNodeGetContent (xml_final_state));
    if (xml_notation_plain) strcpy (fs_notation.plain, xmlNodeGetContent (xml_notation_plain));
    if (xml_notation_latex) strcpy (fs_notation.latex, xmlNodeGetContent (xml_notation_latex));
    proc.final_state_notation = &fs_notation;

    strcpy (proc.cross_section_unit, "pb");
    if (xml_cs) {
      xmlChar * prop1 = xmlGetProp (xml_cs, errorPlusStr);
      xmlChar * prop2 = xmlGetProp (xml_cs, errorMinusStr);
      if (prop1) {
        if (1 != sscanf (prop1, "%le", &(proc.cross_section_upper_error))) {
          proc.cross_section_upper_error = -1.0;
        }
      }
      if (prop2) {
        if (1 != sscanf (prop2, "%le", &(proc.cross_section_lower_error))) {
          proc.cross_section_lower_error = -1.0;
        }
      }
      if (1 != sscanf (xmlNodeGetContent (xml_cs), "%le", &(proc.cross_section))) {
        proc.cross_section = 0;
      }
    }
    err *= create_process(hepml_header, &proc);
  }

////////////////////////////////////////////////////////////////////
  {
    int i;
    for (i = 0; i < 2; ++i) {
      lhaef_beam beam;
      lhaef_pdf pdf;
      lhaef_alphas alphas;

      xmlNodePtr xml_process = xmlCopyNode (get_first_son (desc, processStr), 1);
      xmlNodePtr xml_beam1   = xmlCopyNode (get_first_son (xml_process, beamStr[i]), 1);
      xmlNodePtr xml_part    = xmlCopyNode (get_first_son (xml_beam1, (const xmlChar *)"particle"), 1);
      xmlNodePtr xml_energy  = xmlCopyNode (get_first_son (xml_beam1, (const xmlChar *)"energy"), 1);
      xmlNodePtr xml_pdf     = xmlCopyNode (get_first_son (xml_beam1, (const xmlChar *)"pdf"), 1);
      xmlNodePtr xml_qcd     = xmlCopyNode (get_first_son (xml_beam1, (const xmlChar *)"QCDCoupling"), 1);

      xmlNodePtr xml_qcd_lambda = xmlCopyNode (get_first_son (xml_qcd, (const xmlChar *)"Lambda"), 1);
      xmlNodePtr xml_qcd_nflv   = xmlCopyNode (get_first_son (xml_qcd, (const xmlChar *)"NFlavours"), 1);
      xmlNodePtr xml_qcd_nloops = xmlCopyNode (get_first_son (xml_qcd, (const xmlChar *)"NLoopsAlphaS"), 1);

      if (xml_part) strcpy (beam.particle_name, xmlNodeGetContent (xml_part));
      if (xml_part) {
        if (1 != sscanf (xmlGetProp (xml_part, (const xmlChar *)"pdgcode"), "%i", &(beam.particle_pdgcode))) beam.particle_pdgcode = 0;
      }
      if (xml_energy) {
        if (1 != sscanf (xmlNodeGetContent (xml_energy), "%lf", &(beam.energy))) beam.energy = 0;
      }
      strcpy (beam.energy_unit, "GeV");

      the_pbeam[i] = beam.particle_pdgcode;
      the_ebeam[i] = beam.energy;

      strcpy (alphas.lambda_unit, "GeV");
      alphas.description[0] = 0;
      if (xml_qcd_nloops) {
        if (1 != sscanf (xmlNodeGetContent (xml_qcd_nloops), "%i", &(alphas.nloops))) alphas.nloops = 0;
      }
      if (xml_qcd_nflv) {
        if (1 != sscanf (xmlNodeGetContent (xml_qcd_nflv), "%i", &(alphas.nflavours))) alphas.nflavours = 0;
      }
      if (xml_qcd_lambda) {
        if (1 != sscanf (xmlNodeGetContent (xml_qcd_lambda), "%lf", &(alphas.lambda))) alphas.lambda = 0;
      }

      pdf.name[0] = 0;
      pdf.version[0] = 0;
      pdf.lhapdf_filename[0] = 0;
      pdf.pdflib_set = -1;
      pdf.pdflib_group = -1;
      pdf.lhapdf_set = -1;
      pdf.lhapdf_member = -1;

      if (2212 == abs(beam.particle_pdgcode)) {
        xmlChar * prop1 = xmlGetProp (xml_pdf, nameStr);
        xmlChar * prop2 = xmlGetProp (xml_pdf, (const xmlChar *)"version");
        xmlChar * prop3 = xmlGetProp (xml_pdf, (const xmlChar *)"PDFLIBset");
        xmlChar * prop4 = xmlGetProp (xml_pdf, (const xmlChar *)"PDFLIBgroup");
        xmlChar * prop5 = xmlGetProp (xml_pdf, (const xmlChar *)"LHAPDFset");
        xmlChar * prop6 = xmlGetProp (xml_pdf, (const xmlChar *)"LHAPDFmember");

        pdf.name[0] = 0;
        pdf.version[0] = 0;
        if (prop1) strcpy (pdf.name,    prop1);
        if (prop2) strcpy (pdf.version, prop2);
        if (prop3) if (1 != sscanf (prop3, "%i", &(pdf.pdflib_set)))    pdf.pdflib_set = -1;
        if (prop4) if (1 != sscanf (prop4, "%i", &(pdf.pdflib_group)))  pdf.pdflib_group = -1;
        if (prop5) if (1 != sscanf (prop5, "%i", &(pdf.lhapdf_set)))    pdf.lhapdf_set = -1;
        if (prop6) if (1 != sscanf (prop6, "%i", &(pdf.lhapdf_member))) pdf.lhapdf_member = -1;
        err *= add_process_beam (hepml_header, i, &beam);
        err *= add_process_beam_pdf (hepml_header, i, &pdf);
        err *= add_process_beam_pdf_alphas (hepml_header, &alphas);
      } else {
        err *= add_process_beam (hepml_header, i, &beam);
        err *= add_process_alphas (hepml_header, &alphas);
      }

      the_PDFLIBgroup[i] = pdf.pdflib_group;
      the_PDFLIBset[i]   = pdf.pdflib_set;
    }
  }

  return err;
}

char * prepare_hepml_header_libxml2_dynamic (void) {
  int err = 0;
  void * hepml_header = init_lhaef_document ();

////////////////////////////////////////////////////////////////////
  xmlNodePtr son = mdlparlist->children;
  while (son) {
    if (!xmlStrcmp (son->name, parStr)) {
      lhaef_parameter par;
      lhaef_math_notation not;
      xmlNodePtr xml_model_par_name = xmlCopyNode (get_first_son (son, nameStr), 1);
      xmlNodePtr xml_model_par_valu = xmlCopyNode (get_first_son (son, valueStr), 1);
      xmlNodePtr xml_model_par_desc = xmlCopyNode (get_first_son (son, descriptionStr), 1);
      xmlNodePtr xml_notation_plain = xmlCopyNode (get_first_grandson (son, symbStr, plainStr), 1);
      xmlNodePtr xml_notation_latex = xmlCopyNode (get_first_grandson (son, symbStr, latexStr), 1);

      par.name[0]        = 0;
      par.value[0]       = 0;
      par.description[0] = 0;
      not.plain[0]       = 0;
      not.latex[0]       = 0;

      if (xml_model_par_name) strcpy (par.name,        xmlNodeGetContent (xml_model_par_name));
      if (xml_model_par_valu) strcpy (par.value,       xmlNodeGetContent (xml_model_par_valu));
      if (xml_model_par_desc) strcpy (par.description, xmlNodeGetContent (xml_model_par_desc));
      if (xml_notation_plain) strcpy (not.plain,       xmlNodeGetContent (xml_notation_plain));
      if (xml_notation_latex) strcpy (not.latex,       xmlNodeGetContent (xml_notation_latex));
      par.notation = &not;
      err = add_model_parameter (hepml_header, &par);
    }
    son = son->next;
  }

////////////////////////////////////////////////////////////////////
  son = cutlist->children;
  while (son) {
    if (!xmlStrcmp (son->name, (const xmlChar *)"cutSet")) {
      int cutsetid = 0;
      xmlNodePtr grandson = son->children;
      xmlChar * prop = xmlGetProp (son, cutsetidStr);
      if (NULL != prop) {
        if (1 != sscanf (prop, "%i", &cutsetid)) {
          cutsetid = 0;
        }
      }
      err = add_cutset (hepml_header, cutsetid);
      while (grandson) {
        if (!xmlStrcmp (grandson->name, cutStr)) {
          int minon = 0;
          int maxon = 0;
          lhaef_cut cut;
          lhaef_math_notation not;
          xmlNodePtr xml_cut_object      = xmlCopyNode (get_first_son (grandson, (const xmlChar *)"object"), 1);
          xmlNodePtr xml_cut_min         = xmlCopyNode (get_first_son (grandson, (const xmlChar *)"minValue"), 1);
          xmlNodePtr xml_cut_max         = xmlCopyNode (get_first_son (grandson, (const xmlChar *)"maxValue"), 1);
          xmlNodePtr xml_cut_logic       = xmlCopyNode (get_first_son (grandson, (const xmlChar *)"logic"), 1);
          xmlNodePtr xml_cut_object_name = xmlCopyNode (get_first_son (xml_cut_object, nameStr), 1);
          xmlNodePtr xml_notation_plain  = xmlCopyNode (get_first_grandson (xml_cut_object, symbStr, plainStr), 1);
          xmlNodePtr xml_notation_latex  = xmlCopyNode (get_first_grandson (xml_cut_object, symbStr, latexStr), 1);

          cut.name[0] = 0;
          cut.cut_logic[0] = 0;
          not.plain[0] = 0;
          not.latex[0] = 0;
          if (xml_cut_object_name) {
            strcpy (cut.name, xmlNodeGetContent (xml_cut_object_name));
          }
          if (xml_cut_min) {
            minon = 1;
            if (1 != sscanf (xmlNodeGetContent (xml_cut_min), "%lf", &(cut.min_value))) {
              cut.min_value = 0;
            }
          }
          if (xml_cut_max) {
            maxon = 1;
            if (1 != sscanf (xmlNodeGetContent (xml_cut_max), "%lf", &(cut.max_value))) {
              cut.max_value = 0;
            }
          }
          if (xml_cut_logic)      strcpy (cut.cut_logic, xmlNodeGetContent (xml_cut_logic));
          if (xml_notation_plain) strcpy (not.plain, xmlNodeGetContent (xml_notation_plain));
          if (xml_notation_latex) strcpy (not.latex, xmlNodeGetContent (xml_notation_latex));
          cut.notation = &not;

          err = add_cut (hepml_header, &cut, minon, maxon);
        }
        grandson = grandson->next;
      }
    }
    son = son->next;
  }

////////////////////////////////////////////////////////////////////
  son = sublist->children;
  while (son) {
    int the_num;
    if (!xmlStrcmp (son->name, subprocStr)) {
      lhaef_subprocess subproc;
      lhaef_math_notation scale_notation;
      if (1 != sscanf (xmlGetProp (son, cutsetidStr), "%i", &(the_num))) {
        the_num = 0;
      }

      xmlNodePtr xml_sub_name   = xmlCopyNode (get_first_son (son, symbStr), 1);
      xmlNodePtr xml_sub_cs     = xmlCopyNode (get_first_son (son, (const xmlChar *)"crossSection"), 1);
      xmlNodePtr xml_sub_fscale = xmlCopyNode (get_first_son (son, (const xmlChar *)"FactorisationScale"), 1);
      xmlNodePtr xml_sub_rscale = xmlCopyNode (get_first_son (son, (const xmlChar *)"RenormalisationScale"), 1);

      subproc.name[0] = 0;
      if (xml_sub_name) strcpy (subproc.name, xmlNodeGetContent (xml_sub_name));
      if (1 != sscanf (xmlNodeGetContent (xml_sub_cs), "%le", &(subproc.cross_section))) {
        subproc.cross_section = 0;
      }
      strcpy (subproc.cross_section_unit, "pb");

      if (1 != sscanf (xmlGetProp (xml_sub_cs, errorPlusStr), "%le", &(subproc.cross_section_upper_error))) {
        subproc.cross_section_upper_error = -1.0;
      }
      if (1 != sscanf (xmlGetProp (xml_sub_cs, errorMinusStr), "%le", &(subproc.cross_section_lower_error))) {
        subproc.cross_section_lower_error = -1.0;
      }

      scale_notation.plain[0] = 0;
      scale_notation.latex[0] = 0;
      xmlNodePtr xml_notation_plain = xmlCopyNode (get_first_son (xml_sub_fscale, plainStr), 1);
      xmlNodePtr xml_notation_latex = xmlCopyNode (get_first_son (xml_sub_fscale, latexStr), 1);
      if (xml_notation_plain) strcpy (scale_notation.plain, xmlNodeGetContent (xml_notation_plain));
      if (xml_notation_latex) strcpy (scale_notation.latex, xmlNodeGetContent (xml_notation_latex));
      subproc.FactorisationScale = &scale_notation;

      scale_notation.plain[0] = 0;
      scale_notation.latex[0] = 0;
      xmlNodePtr xml_notation_plain2 = xmlCopyNode (get_first_son (xml_sub_rscale, plainStr), 1);
      xmlNodePtr xml_notation_latex2 = xmlCopyNode (get_first_son (xml_sub_rscale, latexStr), 1);
      if (xml_notation_plain2) strcpy (scale_notation.plain, xmlNodeGetContent (xml_notation_plain2));
      if (xml_notation_latex2) strcpy (scale_notation.latex, xmlNodeGetContent (xml_notation_latex2));
      subproc.RenormalisationScale = &scale_notation;

      err = add_subprocess (hepml_header, &subproc, the_num);
    }
    son = son->next;
  }

  return form_hepml_document (hepml_header);
}


static char * prepare_file_header_libxml2_dynamic (void) {
  char * filestr = malloc ((512 + strlen (get_file_description ())) * sizeof (char));
  strcpy (filestr, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  strcat (filestr, "<samples xmlns=\"http://mcdb.cern.ch/hepml/0.2/\"\n");
  strcat (filestr, "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
  strcat (filestr, "    xsi:schemaLocation=\"http://mcdb.cern.ch/hepml/0.2/ file:/home/dudko/MCDB/hepml/schemas/0.2/hepml.xsd\">\n");
  strcat (filestr, "    <files>\n");
  strcat (filestr, "        <file>\n");
  strcat (filestr, get_file_description ());
  strcat (filestr, "        </file>\n");
  strcat (filestr, "    </files>\n");

  return filestr;
}


int write_file_header_libxml2 (const char fname[], const char mode[]) {

  FILE * outFile = fopen (fname, mode);
  if (!outFile) {
    return 0;
  }

  fprintf (outFile, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (outFile, "<!-- File generated with CompHEP %s -->\n", the_version);
  fprintf (outFile, "<!-- \n"
                    "     This file is compatible with the Les Houches event file\n"
                    "     format (hep-ph/0609017), but contains extra HepML tags.\n"
                    "-->\n");
  fprintf (outFile, "<header>\n");
  fprintf (outFile, "%s", prepare_file_header_libxml2_dynamic ());
  fclose (outFile);

  return 0;
}
#endif


int analyzeLHEfile (const char fname[], int i) {
  int j;
  long pos;
  char buff[2048];
  char * event = NULL;
  FILE * source = fopen (fname, "r");
  char * tmp;

  tot_proc_num = 0;
  if (!source) return -1;
  fgets (buff, 2048, source);
  if (!strstr (buff, "<LesHouchesEvents version=\"1.0\">")) return -2;

  cs[i]    = -1.;
  cserr[i] = -1.;
  nevents[i] = -1;
  while (1) {
    pos = ftell (source);
    fgets (buff, 2048, source);
    if (strstr (buff, "<!-- File generated with CompHEP")) {
       int len;
       char * ttt;
       char ver[1024];
       strcpy (ver, buff + 33);
       ttt = strstr(ver, "-->");
       len = strlen(ver) - strlen (ttt);
       strncpy (the_version, ver, len);
       the_version[len] = 0;
    }
    tmp = strstr (buff, "<crossSection unit=\"pb\"");
    if (tmp || strstr (buff, "<init>")) {
      break;
    }
  }
  if (tmp) {
    double d1, d2, d3;
    if (3 == sscanf (buff, "            <crossSection unit=\"pb\" errorMinus=\"%le\" errorPlus=\"%le\">%le</crossSection>", &d1, &d2, &d3)) {
      cs[i]    = d3;
      cserr[i] = d1;
    }
  }

  rewind (source);
  while (1) {
    pos = ftell (source);
    fgets (buff, 2048, source);
    tmp = strstr (buff, "<eventsNumber>");
    if (tmp || strstr (buff, "<init>")) {
      break;
    }
  }
  if (tmp) {
    int i1;
    if (1 == sscanf (buff, "            <eventsNumber>%d</eventsNumber>", &i1)) {
      nevents[i] = i1;
    }
  }

  fseek (source, pos, SEEK_SET);
  while (!event) {
    fgets (buff, 2048, source);
    if (feof (source)) {
      fclose (source);
      return -4;
    }
    event = strstr (buff, "<init>");
  }
  fgets (buff, 1024, source);
  if (1) {
    int i1, i2, i3, i4, i5, i6, i7;
    double d1, d2;
    if (9 == sscanf (buff, "%d %d %le %le %d %d %d %d 3 %i", &i1, &i2, &d1, &d2, &i3, &i4, &i5, &i6, &i7)) {
      the_pbeam[0]       = i1;
      the_pbeam[1]       = i2;
      the_ebeam[0]       = d1;
      the_ebeam[1]       = d2;
      the_PDFLIBgroup[0] = i3;
      the_PDFLIBgroup[1] = i4;
      the_PDFLIBset[0]   = i5;
      the_PDFLIBset[1]   = i6;
      tot_proc_num       = i7;
    }
  }

  double cserr2 = 0.;
  cs[i]    = 0.;
  cserr[i] = 0.;
  for (j = 0; j < tot_proc_num; ++j) {
    int i1;
    double d1, d2, d3;
    fgets (buff, 1024, source);
    if (4 == sscanf (buff, "%le %le %le %d\n", &d1, &d2, &d3, &i1)) {
      cs[i]  += d1;
      cserr2 += d2 * d2;
    }
    cserr[i] = sqrt (cserr2);
  }

  event = NULL;
  while (!event) {
    pos = ftell (source);
    fgets (buff, 2048, source);
    if (feof (source)) {
      fclose (source);
      return -4;
    }
    event = strstr (buff, "<event>");
  }
  nposition[i] = pos;

  if (nevents[i] < 0) {
    nevents[i] = 0;
    event = NULL;
    while (!feof (source)) {
      fgets (buff, 2048, source);
      event = strstr (buff, "</event>");
      if (event) {
        ++nevents[i];
      }
      buff[0] = 0;
    }
  }

  fclose (source);
  return 0;
}


int analyzeMultiProcLHEfile (const char fname[]) {
  long pos;
  char buff[2048];
  char * event = NULL;
  FILE * source = fopen (fname, "r");
  int numproc = 0;

  tot_proc_num = 0;
  if (!source) return -1;
  fgets (buff, 2048, source);
  if (!strstr (buff, "<LesHouchesEvents version=\"1.0\">")) return -2;

  while (1) {
    if (feof (source)) {
      fclose (source);
      return -4;
    }
    fgets (buff, 2048, source);
    if (strstr (buff, "<init>")) {
      break;
    }
  }

  fgets (buff, 1024, source);
  if (1) {
    int i1, i2, i3, i4, i5, i6, i7;
    double d1, d2;
    if (9 == sscanf (buff, "%d %d %le %le %d %d %d %d 3 %d", &i1, &i2, &d1, &d2, &i3, &i4, &i5, &i6, &i7)) {
      the_pbeam[0]       = i1;
      the_pbeam[1]       = i2;
      the_ebeam[0]       = d1;
      the_ebeam[1]       = d2;
      the_PDFLIBgroup[0] = i3;
      the_PDFLIBgroup[1] = i4;
      the_PDFLIBset[0]   = i5;
      the_PDFLIBset[1]   = i6;
      tot_proc_num       = i7;
    }
  }

  fgets (buff, 2048, source);
  while (numproc < tot_proc_num) {
    if (feof (source)) {
      fclose (source);
      return -5;
    }
    if (numproc > MAX_FILE_EVENT) {
      fclose (source);
      return -6;
    }
    int i1;
    double d1, d2, d3;
    if (4 == sscanf (buff, "%le %le %le %d\n", &d1, &d2, &d3, &i1)) {
      cs[numproc]    = d1;
      cserr[numproc] = d2;
      ++numproc;
    }
    fgets (buff, 2048, source);
    if (strstr (buff, "</init>")) {
      break;
    }
  }

  event = NULL;
  while (!event) {
    pos = ftell (source);
    fgets (buff, 2048, source);
    if (feof (source)) {
      fclose (source);
      return -4;
    }
    event = strstr (buff, "<event>");
  }
  nposition[0] = pos;

  fclose (source);
  return 0;
}


int setInfoWithoutHEPML (char fname[], FILE * f, int procnum) {
  int n, i;
  double xwgt;
  eventUP evnt;
  char buff[2048];
  char * stop = NULL;

  while (!stop) {
    fgets (buff, 2048, f);
    if (feof (f)) return -4;
    stop = strstr (buff, "<init>");
  }
  fgets (buff, 2048, f);
  fgets (buff, 2048, f);
  if (4 != sscanf (buff, "%le %le %le %d", &cs[procnum], &cserr[procnum], &xwgt, &n)) {
    cs[procnum] = 0.0;
    cserr[procnum] = 0.0;
  }

  stop = NULL;
  while (!stop) {
    nposition[procnum] = ftell (f);
    fgets (buff, 2048, f);
    if (feof (f)) return -4;
    stop = strstr (buff, "<event>");
  }

  if (0 == getLHAevent (fname, f, nposition[procnum], &evnt)) {
      sprintf (subproc_name, "%s, %s -> ", kfname(evnt.IDpartUP[0]), kfname(evnt.IDpartUP[1]));
    for (i = 2; i < evnt.NpartUP; ++i) {
      strcat (subproc_name, kfname(evnt.IDpartUP[i]));
      if (i < evnt.NpartUP - 1) strcat (subproc_name, ", ");
    }
  }

  return 0;
}

long getEventPosition (int procnum) {
  return nposition[procnum];
}

char * getFinalStateName (void) {
  return final_state_name;
}

char * getSubprocessName (void) {
  return subproc_name;
}

int getEventNumber (int procnum) {
  return nevents[procnum];
}

double getCrossSection (int procnum) {
  return cs[procnum];
}

double getCrossSectionErr (int procnum) {
  return cserr[procnum];
}

int getPbeam (int i) {
  return the_pbeam[i];
}

double getEbeam (int i) {
  return the_ebeam[i];
}

int getPDFLIBgroup (int i) {
  return the_PDFLIBgroup[i];
}

int getPDFLIBset (int i) {
  return the_PDFLIBset[i];
}

char * getCHEPversion (void) {
  return the_version;
}

int getTotProcNumber (void) {
  return tot_proc_num;
}

