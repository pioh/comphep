/*
* Copyright (C) 2008-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include "LesHouches.h"

#ifdef ROOTused
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TObject.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"

#include "rtuple_routines.h"


static TTree * tree;
static TFile * file;

static Int_t nprt;
static Int_t procid;
static Int_t kf[RTPLMAXINOUT];
static Int_t status[RTPLMAXINOUT];
static Int_t mother1[RTPLMAXINOUT];
static Int_t mother2[RTPLMAXINOUT];
static Int_t color1[RTPLMAXINOUT];
static Int_t color2[RTPLMAXINOUT];

static Double_t weight;
static Double_t scale;
static Double_t QEDalpha;
static Double_t QCDalpha;
static Double_t px[RTPLMAXINOUT];
static Double_t py[RTPLMAXINOUT];
static Double_t pz[RTPLMAXINOUT];
static Double_t e[RTPLMAXINOUT];
static Double_t m[RTPLMAXINOUT];
static Double_t lifetime[RTPLMAXINOUT];
static Double_t spin[RTPLMAXINOUT];

using namespace std;
using namespace ROOT;

#ifdef __cplusplus
extern "C" {
#endif
int
book_rtuple (char rtuple_name[])
{
  int nnet = 0;
  int comp = 1;
  // Authorize Trees up to 2 Terabytes (if the system can do it)
  TTree::SetMaxTreeSize (1000 * Long64_t (2000000000));

  if (nnet) {
    file = new TNetFile (rtuple_name, "RECREATE", "CompHEP ROOT file");
//    file->UseCache (10);
  } else {
    file = new TFile (rtuple_name, "RECREATE", "CompHEP ROOT file");
  }
  file->SetCompressionLevel (comp);

  tree = new TTree ("Tchep", "CompHEP LHA-I COMMON BLOCK based ROOT tree");
//  tree->SetAutoSave (1048576);   // auto-save as 1Mb written
  tree->SetDirectory (file);

  tree->Branch ("nprt",     &nprt,     "nprt/I");
  tree->Branch ("procid",   &procid,   "procid/I");
  tree->Branch ("kf",        kf,       "kf[nprt]/I");
  tree->Branch ("status",    status,   "status[nprt]/I");
  tree->Branch ("mother1",   mother1,  "mother1[nprt]/I");
  tree->Branch ("mother2",   mother2,  "mother2[nprt]/I");
  tree->Branch ("color1",    color1,   "color1[nprt]/I");
  tree->Branch ("color2",    color2,   "color2[nprt]/I");
  tree->Branch ("weight",   &weight,   "weight/D");
  tree->Branch ("scale",    &scale,    "scale/D");
  tree->Branch ("QEDalpha", &QEDalpha, "QEDalpha/D");
  tree->Branch ("QCDalpha", &QCDalpha, "QCDalpha/D");
  tree->Branch ("px",        px,       "px[nprt]/D");
  tree->Branch ("py",        py,       "py[nprt]/D");
  tree->Branch ("pz",        pz,       "pz[nprt]/D");
  tree->Branch (" e",        e,        "e[nprt]/D");
  tree->Branch (" m",        m,        "m[nprt]/D");
  tree->Branch ("lifetime",  lifetime, "lifetime[nprt]/D");
  tree->Branch ("spin",      spin,     "spin[nprt]/D");

  return 0;
}

int 
write_rtuple (void)
{
  file = tree->GetCurrentFile ();
  file->Write ();
//  tree->Print ();
  return 0;
}


int
fill_event (eventUP * ev)
{
  int i;

  nprt     = ev->NpartUP;
  procid   = ev->IDprocUP;
  weight   = ev->XweightUP;
  scale    = ev->QscaleUP;
  QEDalpha = ev->QEDalphaUP;
  QCDalpha = ev->QCDalphaUP;

  for (i = 0; i < nprt; ++i) {
    kf[i]       = ev->IDpartUP[i];
    status[i]   = ev->statusUP[i];
    mother1[i]  = ev->motherUP[0][i];
    mother2[i]  = ev->motherUP[1][i];
    color1[i]   = ev->colorflowUP[0][i];
    color2[i]   = ev->colorflowUP[1][i];
    px[i]       = ev->momentumUP[0][i];
    py[i]       = ev->momentumUP[1][i];
    pz[i]       = ev->momentumUP[2][i];
    e[i]        = ev->momentumUP[3][i];
    m[i]        = ev->momentumUP[4][i];
    lifetime[i] = ev->timelifeUP[i];
    spin[i]     = ev->spinUP[i];
  }
  tree->Fill ();

  return 0;
}
}


#else
int
book_rtuple (char rtuple_name[])
{
  fprintf (stdout, "rtupler (error): ROOT is not linked. Define ROOTSYS and re-compile CompHEP. Exit.\n");
  return 0;
}

int 
write_rtuple (void)
{
  fprintf (stdout, "rtupler (error): ROOT is not linked. Define ROOTSYS and re-compile CompHEP. Exit.\n");
  return 0;
}


int
fill_event (eventUP * ev)
{
  fprintf (stdout, "rtupler (error): ROOT is not linked. Define ROOTSYS and re-compile CompHEP. Exit.\n");
  return 0;
}
#endif
