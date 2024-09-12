/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pwd.h>

#include "service2/include/chep_limits.h"
#include "service2/include/drandXX.h"
#include "service2/include/syst.h"
#include "service2/include/lbl.h"
#include "service2/include/parseref.h"
#include "service2/include/4_vector.h"
#include "service2/include/kfcodes.h"

#include "out_ext.h"

#include "alphas_menu.h"
#include "alphas2.h"
#include "cut.h"
#include "core_data.h"
#include "lhaef.h"
#include "kinaux.h"
#include "runVegas.h"
#include "strfun.h"
#include "subproc.h"
#include "vegas.h"
#include "evnt_tools.h"
#include "evnt_format_lhaef.h"

static FILE * events_;
static int * cChains;
static int nC;
static int cBasisPower;

static int pbeam[2];
static double ebeam[2];
static Str_fun_Info strfun[2];

static int nevents = 0;
static double cs = 0.0;
static double cserr = 0.0;
static int filesize = 0;
static char check_sum[128];

static char * prepare_hepml_header (void)
{
  int i;
  int err = 0;
  void * hepml_header = init_lhaef_document ();
  char * machine_name = "";

  lhaef_file lhafile;
  strcpy (lhafile.location_path, "");
  strcpy (lhafile.location_type, "local");
  strcpy (lhafile.machine_name, machine_name);
  lhafile.nevents = nevents;
  lhafile.size = filesize;
  lhafile.cross_section = cs;
  lhafile.cross_section_upper_error = cserr * cs;
  lhafile.cross_section_lower_error = cserr * cs;
  strcpy (lhafile.cross_section_unit, "pb");
  strcpy (lhafile.check_sum, check_sum);
  strcpy (lhafile.check_sum_type, "md5");
  err = set_file_description(hepml_header, &lhafile);

  lhaef_general_description gendesc;
  lhapdf_experiment exp;
  strcpy (gendesc.title, "");
  strcpy (gendesc.abstruct, "");
  strcpy (gendesc.comments, "");
  strcpy (exp.experiment, get_author_experiment(0));
  strcpy (exp.group, get_author_group(0));
  strcpy (exp.responsiblePerson, "");
  strcpy (exp.description, "");
  gendesc.experiment = &exp;
  err = set_general_description(hepml_header, &gendesc);

  lhaef_generator gen;
  strcpy (gen.name, "CompHEP");
  strcpy (gen.version, version());
  strcpy (gen.homepage, "http://comphep.sinp.msu.ru");
  strcpy (gen.description, 
    "                CompHEP: a package for evaluation of Feynman diagrams, integration over multi-particle \n"
    "                phase space and event generation (supported in part by RFBR grants 96-02-19773-a, \n"
    "                99-02-04011-a, 01-02-16710-a, 04-02-17448-a). CompHEP Collaboration: \n"
    "                E.Boos, V.Bunichev, M.Dubinin, L.Dudko, V.Edneral, V.ILyin, A.Kryukov, V.Savrin \n"
    "                (SINP MSU, Moscow, Russia), A.Semenov (JINR, Dubna, Russia), A.Sherstnev \n"
    "                (SINP MSU, Moscow, Russia and University of Cambridge, UK) ");
  err = set_generator(hepml_header, &gen);

  lhaef_model model;
  strcpy (model.name, "");
  strcpy (model.description, "");
  err = set_model(hepml_header, &model);
  for (i = 0; i < nvar_ + nfunc_ + 1; ++i) {
    char name[128];
    char value[128];
    double val;
    lhaef_parameter par;
    lhaef_math_notation notation;
    vinf_ (i, name, &val);
    sprintf (value, "%f", val);
    strcpy (par.name, name);
    strcpy (par.value, value);
    strcpy (par.description, "");
    strcpy (notation.plain, name);
    strcpy (notation.html, "");
    strcpy (notation.latex, "");
    strcpy (notation.mathml, "");
    par.notation = &notation;
    err = add_model_parameter(hepml_header, &par);
  }

  err = add_cutset (hepml_header, 1);
  for (i = 0; invcut_1[i].key; ++i) {
    int j;
    vshortstr parts;
    lhaef_cut cut;
    lhaef_math_notation notation;
    for(j = 0; invcut_1[i].lvinvc[j]; ++j) {
      parts[j] = invcut_1[i].lvinvc[j] + '0';
    }
    parts[j] = '\0';
    sprintf (cut.name, "%c%s", invcut_1[i].key, parts);
    if (invcut_1[i].minon) cut.min_value = invcut_1[i].cvmin;
    if (invcut_1[i].maxon) cut.max_value = invcut_1[i].cvmax;
    strcpy (cut.cut_logic, "include");
    strcpy (notation.plain, cut.name);
    strcpy (notation.html, "");
    strcpy (notation.latex, "");
    strcpy (notation.mathml, "");
    cut.notation = &notation;
    err = add_cut (hepml_header, &cut, invcut_1[i].minon, invcut_1[i].maxon);
    cut.notation = NULL;
  }
/*
  err = set_author_record (hepml_header);
  for (i = 0; i < getnauthor (); ++i) {
    lhaef_author author;
    strcpy (author.firstName, get_author_firstname(i));
    strcpy (author.lastName, get_author_lastname(i));
    strcpy (author.email, get_author_email(i));
    strcpy (author.experiment, get_author_experiment(i));
    strcpy (author.group, get_author_group(i));
    strcpy (author.organization, get_author_organization(i));
    err = add_author (hepml_header, &author);
  }
*/
/* Adopt simpler solution */
  err = set_author_record (hepml_header);
  {
    lhaef_author author;
    struct passwd * userrecord = getpwuid (getuid());
    strcpy (author.firstName, userrecord->pw_gecos);
    strcpy (author.lastName, userrecord->pw_gecos);
    strcpy (author.email, "");
    strcpy (author.experiment, "");
    strcpy (author.group, "");
    strcpy (author.organization, "");
    err = add_author (hepml_header, &author);
  }

  lhaef_process proc;
  lhaef_math_notation fs_notation;
  strcpy (proc.final_state, get_final_state_name ());
  strcpy (fs_notation.plain, get_final_state_name ());
  strcpy (fs_notation.html, "");
  strcpy (fs_notation.latex, "");
  strcpy (fs_notation.mathml, "");
  proc.final_state_notation = &fs_notation;
  proc.cross_section = cs;
  proc.cross_section_upper_error = cserr * cs;
  proc.cross_section_lower_error = cserr * cs;
  strcpy (proc.cross_section_unit, "pb");
  err = create_process(hepml_header, &proc);

  for (i = 0; i < 2; ++i) {
    wrt_sf_NF_ (i, &strfun[i]);
  }

  for (i = 0; i < 2; ++i) {
    lhaef_beam beam;
    lhaef_pdf pdf;
    lhaef_alphas alphas;

    if (2 == nin_) {
      double e, p, sqrtS;
      double rap = get_rapidity ();
      vinf_ (0, NULL, &sqrtS);
      e = (sqrtS * sqrtS + strfun[i].prt_mass * strfun[i].prt_mass - strfun[1 - i].prt_mass * strfun[1 - i].prt_mass) / (2 * sqrtS);
      p = sqrt (e * e - strfun[i].prt_mass * strfun[i].prt_mass);
      p = p * cosh (rap) + e * sinh (rap) * (1 - 2 * i);
      if (p * p < (10.E-10) * sqrtS) p = 0.0;
      ebeam[i] = e;
      pbeam[i] = kfbeam (strfun[i].prt_name);
    } else {
      vshortstr n;
      double q;
      pinf_ (proces_1.nsub, i, NULL, &q);
      ebeam[0] = q;
      ebeam[1] = 0.;
      pinf_ (proces_1.nsub, i, n, NULL);
      pbeam[0] = kfpart (n);
      pbeam[1] = 0;
    }

    strcpy (beam.particle_name, strfun[i].prt_name);
    beam.particle_pdgcode = pbeam[i];
    beam.energy = ebeam[i];
    strcpy (beam.energy_unit, "GeV");
    err = add_process_beam (hepml_header, i, &beam);

    strcpy (alphas.lambda_unit, "GeV");
    strcpy (alphas.description, "");
    alphas.nloops = 2;
    alphas.nflavours = Nflavour ();
    alphas.lambda = QCDLambda ();
    if (2212 == abs(pbeam[i])) {
      strcpy (pdf.name, strfun[i].pdf_name);
      strcpy (pdf.version, strfun[i].version);
      pdf.pdflib_set = strfun[i].PDFLIBset;
      pdf.pdflib_group = strfun[i].PDFLIBgroup;
      pdf.lhapdf_set = strfun[i].LHAPDFset;
      pdf.lhapdf_member = strfun[i].LHAPDFmember;
      strcpy (pdf.lhapdf_filename, "");
      alphas.nloops = QCDOrder ();
      err = add_process_beam_pdf(hepml_header, i, &pdf);
      err = add_process_beam_pdf_alphas(hepml_header, &alphas);
    } else {
      err = add_process_alphas(hepml_header, &alphas);
    }
  }

  lhaef_subprocess subproc;
  lhaef_math_notation scale_notation;
  strcpy (subproc.name, get_subproc_name());
  subproc.cross_section = cs;
  subproc.cross_section_upper_error = cserr * cs;
  subproc.cross_section_lower_error = cserr * cs;
  strcpy (subproc.cross_section_unit, "pb");
  strcpy (scale_notation.plain, get_scale_form());
  strcpy (scale_notation.html, "");
  strcpy (scale_notation.latex, "");
  strcpy (scale_notation.mathml, "");
  subproc.FactorisationScale = &scale_notation;
  subproc.RenormalisationScale = &scale_notation;
  err = add_subprocess(hepml_header, &subproc, 1);

  return form_hepml_document (hepml_header);
}


static void write_event_lhaef (long cCub, int n, double w) {
  int i;
  int j;
  int icc;
  double f;
  int ntot_ = nout_ + nin_;
  double alphaem = 1/128.;
  int col[MAXINOUT][2];
  int IPSGN = 1;

  for (j = 0; j < MAXINOUT; ++j) col[j][0] = col[j][1] = 0;

  if (cBasisPower) {
    int itag = 500;
    double sum = 0;
    for (i = 0; i < cBasisPower; i++) {
      sum += fabs (color_weights[i]);
    }
    sum *= drandXX ();
    for (i = 0; i < cBasisPower; i++) {
      sum -= fabs (color_weights[i]);
      if (sum <= 0) {
        break;
      }
    }
    if (i == cBasisPower) {
      i--;
    }
    if (color_weights[i] < 0) {
      n *= -1;
    }
    icc = i;

    for (j = 0; j < nC; ++j) {
      int n1 = cChains[2 * (nC * icc + j)] - 1;
      int n2 = cChains[2 * (nC * icc + j) + 1] - 1;

      if (2 == nin_) {
        if (IPSGN == -1) {
          if (n1 == 0)      n1 = 1;
          else if (n1 == 1) n1 = 0;
          if (n2 == 0)      n2 = 1;
          else if (n2 == 1) n2 = 0;
        }

        if (n1 < 2) col[n1][0] = itag;
        else        col[n1][1] = itag;
        if (n2 < 2) col[n2][1] = itag;
        else        col[n2][0] = itag;
      } else {
        if (n1 < 1) col[n1][0] = itag;
        else        col[n1][1] = itag;
        if (n2 < 1) col[n2][1] = itag;
        else        col[n2][0] = itag;
      }
      ++itag;
    }
  }

  f = calcCutFactor ();
  if (!f) {
    fprintf (stderr, "Error : generator is trying to keep an event \n");
    fprintf (stderr, "        which does not pass cuts.\n");
    printDetailsOfCutFactor ();
  }

  for (i = 0; i < n; ++i) {
    int id;
    double m;
    vshortstr pname;
    longstr init[2];
    double qcdscale = qcd_Scale_chep ();
    rnd_rotate_momentum (nin_, nout_);

    fprintf (events_, "<event>\n");
    fprintf (events_, "%i 1 %17.10E %17.10E %17.10E %17.10E\n", ntot_, w, qcdscale, alpha_2 (qcdscale), alphaem);
    for (j = 0; j < nin_; ++j) {
      int k = 4 * j;
      pinf_ (proces_1.nsub, j + 1, pname, &m);
      id = kfpart (pname);
      pvect[k + 3] *= IPSGN;
      sprintf (init[j], "%i -1 0 0 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
      id, col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], m);
    }
    if (1 < nin_) {
      if (-1 == IPSGN) {
        fputs (init[1], events_);
        fputs (init[0], events_);
      } else {
        fputs (init[0], events_);
        fputs (init[1], events_);
      }
    } else {
      fputs (init[0], events_);
    }
    for (j = nin_; j < ntot_; ++j) {
      int k = 4 * j;
      pinf_ (proces_1.nsub, j + 1, pname, &m);
      id = kfpart (pname);
      pvect[k + 3] *= IPSGN;
      if (2 == nin_) {
        fprintf (events_, "%i 1 1 2 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
        id, col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], m);
      } else {
        fprintf (events_, "%i 1 1 0 %i %i %17.10E %17.10E %17.10E %17.10E %17.10E 0.0 9.0\n", 
        id, col[j][0], col[j][1], pvect[k + 1], pvect[k + 2], pvect[k + 3], pvect[k], m);
      }
    }
    fprintf (events_, "</event>\n");
  }
}


int write_event_cap_lhef (char fname[], char mode[]) {
  int nprup = 1;
  int idwtup = 3;
  long pos;

  FILE * outFile = fopen (fname, mode);
  if (!outFile) {
    return 0;
  }

  fprintf (outFile, "<LesHouchesEvents version=\"1.0\">\n");
  fprintf (outFile, "<!-- File generated with CompHEP %s -->\n", version ());
  fprintf (outFile, "<!-- \n"
                    "     Preliminary version, it is compatible with the Les Houches event file\n"
                    "     format (hep-ph/0609017), but contains extra tags.\n"
                    "-->\n");
  fprintf (outFile, "<header>\n");

  fprintf (outFile, "%s", prepare_hepml_header ());
  fprintf (outFile, "</header>\n");
  fprintf (outFile, "<init>\n");
  fprintf (outFile, "%i %i %17.10E %17.10E %i %i %i %i %i %i\n",
                    pbeam[0], 
                    pbeam[1], 
                    ebeam[0], 
                    ebeam[1], 
#ifdef LHAPDF
                    strfun[0].LHAPDFset, 
                    strfun[1].LHAPDFset, 
                    strfun[0].LHAPDFmember, 
                    strfun[1].LHAPDFmember, 
#else
                    strfun[0].PDFLIBgroup, 
                    strfun[1].PDFLIBgroup, 
                    strfun[0].PDFLIBset, 
                    strfun[1].PDFLIBset, 
#endif
                    idwtup, 
                    nprup);
  fprintf (outFile, "%17.10E %17.10E %17.10E %i\n", cs, cserr * cs, 1.0, 1);
  fprintf (outFile, "</init>\n");

  pos = ftell (outFile);
  fclose (outFile);

  return pos;
}

static long fileEnd;

static int read_event_cap_lhef ()
{
  int n = 0;
  char word[1024];
  fscanf (events_, "%s", word);
  while (!strstr (word, "</header>")) {
    if (NULL != strstr (word, "<eventsNumber>")) {
      fscanf (events_, "%s", word);
      sscanf (word, "%i</eventsNumber>", &n);
    }
    fscanf (events_, "%s", word);
  }
  return n;
}

static long find_last_event_position ()
{
  midstr buff;
  long finpos = ftell (events_);

  fgets (buff, 1024, events_);
  while (!strstr (buff, "</LesHouchesEvents>")) {
    if (feof (events_)) return -1;
    finpos = ftell (events_);
    fgets (buff, 1024, events_);
  }
  return finpos;
}

static int weighted_event = 0;

void set_weighted_flag (int i)
{
  weighted_event = i;
}


int prepare_evfile_lhaef (vegasGrid * vegPtr, double (*func) (double *, double), char * fname, 
   float * cubemaxval, int n_event, int n_cube, double max) {
  int status = 0;

  cStrings (proces_1.nsub, &nC, &cBasisPower, &cChains);
  if (cBasisPower) {
    color_weights = malloc (sizeof (double) * cBasisPower);
  }

  if (n_cube <= 0) {
   status = -2;
  }

  if (n_event <= 0) {
   status = -3;
  }

  if (max <= 0.0) {
   status = -4;
  }

  if (0 == status) {
    events_ = fopen (fname, "a+");
    if (NULL != events_) {
      fseek (events_, 0L, SEEK_END);
      if (0 == ftell (events_)) {
        fclose (events_);
        write_event_cap_lhef (fname, "w");
        events_ = fopen (fname, "a+");
        if (1 == nin_ && 2 == nout_ ) {
          status = vegas_1to2_events (vegPtr, n_cube, n_event, max, func, write_event_lhaef, cubemaxval);
        } else {
          if (weighted_event) {
            status = vegas_wgt (vegPtr, n_cube, n_event, max, func, write_event_lhaef, cubemaxval);
          } else {
            status = vegas_events (vegPtr, n_cube, n_event, max, func, write_event_lhaef, cubemaxval);
          }
        }
        fclose (events_);
        if (0 == status) nevents = n_event;
      } else {
        int chck = CheckFormat (events_);
        if (4 != chck) {
          status = -1;
        } else {
          long flepos = 0L;
          fseek (events_, 0L, SEEK_SET);
          nevents = read_event_cap_lhef ();
          flepos = find_last_event_position ();
          if (0L < flepos) {
            truncate (fname, flepos);
            fileEnd = flepos;
            if (1 == nin_ && 2 == nout_ ) {
              status = vegas_1to2_events (vegPtr, n_cube, n_event, max, func, write_event_lhaef, cubemaxval);
            } else {
              if (weighted_event) {
                status = vegas_wgt (vegPtr, n_cube, n_event, max, func, write_event_lhaef, cubemaxval);
              } else {
                status = vegas_events (vegPtr, n_cube, n_event, max, func, write_event_lhaef, cubemaxval);
              }
            }
            if (0 == status) nevents += n_event;
            fclose (events_);
          } else {
            status = -2;
          }
        }
      }
    } else {
      status = -5;
    }
  }
  return status;
}

int complete_evfile_lhaef (char * fname, int store, int n_event, double mult, double rmax)
{
  struct stat stt;
  if (store) {
    vegas_integral in = get_vegas_integral ();
    in.old = 1;
    set_vegas_integral (in);
    if (in.n_it) {
      cs = in.s1 / in.s0;
      cserr = sqrt (in.s0) / fabs (in.s1);
    }

    stat(fname, &stt);
    filesize = stt.st_size;
    write_event_cap_lhef (fname, "r+");
  } else {
    truncate (fname, fileEnd);
  }
  events_ = fopen (fname, "a+");
  fseek (events_, 0L, SEEK_END);
  if (0 != ftell (events_)) {
    fprintf (events_, "</LesHouchesEvents>\n");
    fclose (events_);
  } else {
    fclose (events_);
    remove (fname);
  }
  return 1;
}
