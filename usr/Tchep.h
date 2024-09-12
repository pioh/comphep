/*
* Copyright (C) 2008-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef Tchep_h
#define Tchep_h

#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

class Tchep {
public :
   TTree *fChain;
   Int_t fCurrent;
   vector<TLorentzVector> prt;

   // Declaration of leaf types
   Int_t     nprt;
   Int_t     procid;
   Int_t     kf[16];
   Int_t     status[16];
   Int_t     mother1[16];
   Int_t     mother2[16];
   Int_t     color1[16];
   Int_t     color2[16];
   Double_t  weight;
   Double_t  scale;
   Double_t  QEDalpha;
   Double_t  QCDalpha;
   Double_t  px[16];
   Double_t  py[16];
   Double_t  pz[16];
   Double_t  e[16];
   Double_t  m[16];
   Double_t  lifetime[16];
   Double_t  spin[16];

   // List of branches
   TBranch * b_nprt;
   TBranch * b_procid;
   TBranch * b_kf;
   TBranch * b_status;
   TBranch * b_mother1;
   TBranch * b_mother2;
   TBranch * b_color1;
   TBranch * b_color2;
   TBranch * b_weight;
   TBranch * b_scale;
   TBranch * b_QEDalpha;
   TBranch * b_QCDalpha;
   TBranch * b_px;
   TBranch * b_py;
   TBranch * b_pz;
   TBranch * b_e;
   TBranch * b_m;
   TBranch * b_lifetime;
   TBranch * b_spin;

   Tchep(Char_t * file_name);
   virtual ~Tchep();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Loop();
   virtual Bool_t   Notify();
};

#endif

#ifdef Tchep_cxx
Tchep::Tchep (Char_t * file_name)
{
 TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject (file_name);
  if (!f) {
    f = new TFile (file_name);
  }
  fChain = (TTree*)gDirectory->Get ("Tchep");

  fCurrent = -1;
  fChain->SetMakeClass (1);

  fChain->SetBranchAddress ("nprt",     &nprt,     &b_nprt);
  fChain->SetBranchAddress ("procid",   &procid,   &b_procid);
  fChain->SetBranchAddress ("kf",       kf,        &b_kf);
  fChain->SetBranchAddress ("status",   status,    &b_status);
  fChain->SetBranchAddress ("mother1",  mother1,   &b_mother1);
  fChain->SetBranchAddress ("mother2",  mother2,   &b_mother2);
  fChain->SetBranchAddress ("color1",   color1,    &b_color1);
  fChain->SetBranchAddress ("color2",   color2,    &b_color2);
  fChain->SetBranchAddress ("weight",   &weight,   &b_weight);
  fChain->SetBranchAddress ("scale",    &scale,    &b_scale);
  fChain->SetBranchAddress ("QEDalpha", &QEDalpha, &b_QEDalpha);
  fChain->SetBranchAddress ("QCDalpha", &QCDalpha, &b_QCDalpha);
  fChain->SetBranchAddress ("px",       px,        &b_px);
  fChain->SetBranchAddress ("py",       py,        &b_py);
  fChain->SetBranchAddress ("pz",       pz,        &b_pz);
  fChain->SetBranchAddress (" e",       e,         &b_e);
  fChain->SetBranchAddress (" m",       m,         &b_m);
  fChain->SetBranchAddress ("lifetime", lifetime,  &b_lifetime);
  fChain->SetBranchAddress ("spin",     spin,      &b_spin);
  Notify();
}

Tchep::~Tchep()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile ();
}

Int_t Tchep::GetEntry(Long64_t entry)
{
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t Tchep::LoadTree(Long64_t entry)
{
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

Bool_t Tchep::Notify()
{
  return kTRUE;
}

#endif // #ifdef Tchep_cxx
