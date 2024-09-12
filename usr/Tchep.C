/*
* Copyright (C) 2008-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ------------------------------------------------------
*/
#define Tchep_cxx

#include <iostream>
#include <fstream>
#include <iomanip>

#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "Tchep.h"

void Tchep::Loop ()
{
  if (fChain == 0) return;

  TH1F * hist = new TH1F ("1", "T_P(t)", 100, 0.0, 100.0);

  Long64_t nentries = fChain->GetEntriesFast();

  cout << "Test analysis code: Pt of the 3d particle in range (0 - 100) GeV is built " << endl;
  for (Long64_t j = 0; j < nentries; ++j) {
    prt.clear ();
    if (0 > LoadTree (j)) break;
    if (0 == j % 10000) cout << "Main Loop: processed event N = " << j << endl;

    fChain->GetEntry (j);
    for (Int_t i = 0; i < nprt; ++i) {
      TLorentzVector p (px[i], py[i], pz[i], e[i]);
      prt.push_back (p);
    }
    hist->Fill (prt[3].Pt ());
  }

  hist->Draw();
}
