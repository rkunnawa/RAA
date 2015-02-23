// Raghav Kunnawalkam Elayavalli
// Feb 21 2014
// Rugers 

//
// Macro to study the effect of jet id cuts. this makes ratio plots from data and MC so we can study the effect. 
//
//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <fstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"


#include "Headers/plot.h"
 
static const int nbins_pt = 39;
static const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 272, 300, 
  330, 362, 395, 430,
  468, 507, 548, 592,
  638, 686, 1000 
};


using namespace std;

void RAA_plot_jetid_checks_ratio_plots(int radius = 3, char *algo = "Pu", char *jet_type = "PF"){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  TFile *fData = TFile::Open("/Users/keraghav/WORK/RAA/Output/PbPb_jetid_checks_akPuPF_20150220.root");
  TFile *fMC = TFile::Open("/Users/keraghav/WORK/RAA/Output/PbPb_mc_chMaxjtpt_norawptcut_spectra_akPuPF_20150220.root");

  //get the histograms: 
  TH1F* hData_Jet80 = (TH1F*)fData->Get("hpbpb_Jet80");
  TH1F* hData_Jet80_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p01");
  TH1F* hData_Jet80_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p02");
  TH1F* hData_Jet80_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p03");
  TH1F* hData_Jet80_chMaxJtpt0p04 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p04");
  TH1F* hData_Jet80_chMaxJtpt0p05 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p05");
  TH1F* hData_Jet80_eMaxJtpt0p1 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p1");
  TH1F* hData_Jet80_eMaxJtpt0p2 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p2");
  TH1F* hData_Jet80_eMaxJtpt0p3 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p3");
  TH1F* hData_Jet80_eMaxJtpt0p4 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p4");
  TH1F* hData_Jet80_eMaxJtpt0p5 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5");
  TH1F* hData_Jet80_eMaxJtpt0p6 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6");
  TH1F* hData_Jet80_eMaxJtpt0p7 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7");
  TH1F* hData_Jet80_eMaxJtpt0p8 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p8");
  TH1F* hData_Jet80_eMaxJtpt0p9 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p9");
  TH1F* hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hData_Jet65 = (TH1F*)fData->Get("hpbpb_Jet65");
  TH1F* hData_Jet65_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p01");
  TH1F* hData_Jet65_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p02");
  TH1F* hData_Jet65_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p03");
  TH1F* hData_Jet65_chMaxJtpt0p04 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p04");
  TH1F* hData_Jet65_chMaxJtpt0p05 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p05");
  TH1F* hData_Jet65_eMaxJtpt0p1 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p1");
  TH1F* hData_Jet65_eMaxJtpt0p2 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p2");
  TH1F* hData_Jet65_eMaxJtpt0p3 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p3");
  TH1F* hData_Jet65_eMaxJtpt0p4 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p4");
  TH1F* hData_Jet65_eMaxJtpt0p5 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5");
  TH1F* hData_Jet65_eMaxJtpt0p6 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6");
  TH1F* hData_Jet65_eMaxJtpt0p7 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7");
  TH1F* hData_Jet65_eMaxJtpt0p8 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p8");
  TH1F* hData_Jet65_eMaxJtpt0p9 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p9");
  TH1F* hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hData_Jet55 = (TH1F*)fData->Get("hpbpb_Jet55");
  TH1F* hData_Jet55_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p01");
  TH1F* hData_Jet55_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p02");
  TH1F* hData_Jet55_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p03");
  TH1F* hData_Jet55_chMaxJtpt0p04 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p04");
  TH1F* hData_Jet55_chMaxJtpt0p05 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p05");
  TH1F* hData_Jet55_eMaxJtpt0p1 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p1");
  TH1F* hData_Jet55_eMaxJtpt0p2 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p2");
  TH1F* hData_Jet55_eMaxJtpt0p3 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p3");
  TH1F* hData_Jet55_eMaxJtpt0p4 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p4");
  TH1F* hData_Jet55_eMaxJtpt0p5 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5");
  TH1F* hData_Jet55_eMaxJtpt0p6 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6");
  TH1F* hData_Jet55_eMaxJtpt0p7 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7");
  TH1F* hData_Jet55_eMaxJtpt0p8 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p8");
  TH1F* hData_Jet55_eMaxJtpt0p9 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p9");
  TH1F* hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03");

  TH2F *hData_chMaxJtpt_jtpt = (TH2F*)fData->Get("hpbpb_chMaxJtpt_jtpt");
  TH2F *hData_eMaxJtpt_jtpt = (TH2F*)fData->Get("hpbpb_eMaxJtpt_jtpt");
  TH2F *hData_eMaxJtpt_chMaxJtpt = (TH2F*)fData->Get("hpbpb_eMaxJtpt_chMaxJtpt");


  //get the histograms: 
  TH1F* hMC_Jet80 = (TH1F*)fMC->Get("hpbpb_Jet80");
  TH1F* hMC_Jet80_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p01");
  TH1F* hMC_Jet80_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p02");
  TH1F* hMC_Jet80_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p03");
  TH1F* hMC_Jet80_chMaxJtpt0p04 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p04");
  TH1F* hMC_Jet80_chMaxJtpt0p05 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p05");
  TH1F* hMC_Jet80_eMaxJtpt0p1 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p1");
  TH1F* hMC_Jet80_eMaxJtpt0p2 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p2");
  TH1F* hMC_Jet80_eMaxJtpt0p3 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p3");
  TH1F* hMC_Jet80_eMaxJtpt0p4 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p4");
  TH1F* hMC_Jet80_eMaxJtpt0p5 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5");
  TH1F* hMC_Jet80_eMaxJtpt0p6 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6");
  TH1F* hMC_Jet80_eMaxJtpt0p7 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7");
  TH1F* hMC_Jet80_eMaxJtpt0p8 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p8");
  TH1F* hMC_Jet80_eMaxJtpt0p9 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p9");
  TH1F* hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hMC_Jet65 = (TH1F*)fMC->Get("hpbpb_Jet65");
  TH1F* hMC_Jet65_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p01");
  TH1F* hMC_Jet65_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p02");
  TH1F* hMC_Jet65_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p03");
  TH1F* hMC_Jet65_chMaxJtpt0p04 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p04");
  TH1F* hMC_Jet65_chMaxJtpt0p05 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p05");
  TH1F* hMC_Jet65_eMaxJtpt0p1 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p1");
  TH1F* hMC_Jet65_eMaxJtpt0p2 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p2");
  TH1F* hMC_Jet65_eMaxJtpt0p3 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p3");
  TH1F* hMC_Jet65_eMaxJtpt0p4 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p4");
  TH1F* hMC_Jet65_eMaxJtpt0p5 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5");
  TH1F* hMC_Jet65_eMaxJtpt0p6 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6");
  TH1F* hMC_Jet65_eMaxJtpt0p7 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7");
  TH1F* hMC_Jet65_eMaxJtpt0p8 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p8");
  TH1F* hMC_Jet65_eMaxJtpt0p9 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p9");
  TH1F* hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hMC_Jet55 = (TH1F*)fMC->Get("hpbpb_Jet55");
  TH1F* hMC_Jet55_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p01");
  TH1F* hMC_Jet55_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p02");
  TH1F* hMC_Jet55_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p03");
  TH1F* hMC_Jet55_chMaxJtpt0p04 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p04");
  TH1F* hMC_Jet55_chMaxJtpt0p05 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p05");
  TH1F* hMC_Jet55_eMaxJtpt0p1 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p1");
  TH1F* hMC_Jet55_eMaxJtpt0p2 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p2");
  TH1F* hMC_Jet55_eMaxJtpt0p3 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p3");
  TH1F* hMC_Jet55_eMaxJtpt0p4 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p4");
  TH1F* hMC_Jet55_eMaxJtpt0p5 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5");
  TH1F* hMC_Jet55_eMaxJtpt0p6 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6");
  TH1F* hMC_Jet55_eMaxJtpt0p7 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7");
  TH1F* hMC_Jet55_eMaxJtpt0p8 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p8");
  TH1F* hMC_Jet55_eMaxJtpt0p9 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p9");
  TH1F* hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03");

  TH2F *hMC_chMaxJtpt_jtpt = (TH2F*)fMC->Get("hpbpb_chMaxJtpt_jtpt");
  TH2F *hMC_eMaxJtpt_jtpt = (TH2F*)fMC->Get("hpbpb_eMaxJtpt_jtpt");
  TH2F *hMC_eMaxJtpt_chMaxJtpt = (TH2F*)fMC->Get("hpbpb_eMaxJtpt_chMaxJtpt");


  hData_Jet80_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p03->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p04->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p05->Divide(hData_Jet80);

  hData_Jet80_eMaxJtpt0p1->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p2->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p3->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p4->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p7->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p8->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p9->Divide(hData_Jet80);

  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hData_Jet80);

  hData_Jet65_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p03->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p04->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p05->Divide(hData_Jet65);

  hData_Jet65_eMaxJtpt0p1->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p2->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p3->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p4->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p7->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p8->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p9->Divide(hData_Jet65);

  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hData_Jet65);

  hData_Jet55_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p03->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p04->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p05->Divide(hData_Jet55);

  hData_Jet55_eMaxJtpt0p1->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p2->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p3->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p4->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p7->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p8->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p9->Divide(hData_Jet55);

  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hData_Jet55);

  // lets start making the plots 3 plots, data on top MC at the bottom: 
  // 3 triggers x 3 different Jet iD, plots with 3x1 showing the   
  TCanvas * cJet55 = new TCanvas("cJet55","",1200,1000);
  makeMultiPanelCanvas(cJet55,3,3,0.0,0.0,0.2,0.15,0.07);
  TLine *line = new TLine(30,1,300,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  
  cJet55->cd(1);
  hData_Jet55_chMaxJtpt0p01->SetTitle(" ");
  hData_Jet55_chMaxJtpt0p01->SetXTitle(" Jet p_{T} (GeV/c)");
  hData_Jet55_chMaxJtpt0p01->SetYTitle(" Cut effiiency (with Cut/without Cut) ");
  hData_Jet55_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet55_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet55_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet55_chMaxJtpt0p01->Draw();
  drawText("Data",0.5,0.8,14);

  hData_Jet55_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_chMaxJtpt0p02->Draw("same");

  hData_Jet55_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_chMaxJtpt0p03->Draw("same");

  hData_Jet55_chMaxJtpt0p04->Rebin(5);
  hData_Jet55_chMaxJtpt0p04->Scale(1./5);
  hData_Jet55_chMaxJtpt0p04->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p04->SetMarkerColor(4);
  hData_Jet55_chMaxJtpt0p04->Draw("same");

  hData_Jet55_chMaxJtpt0p05->Rebin(5);
  hData_Jet55_chMaxJtpt0p05->Scale(1./5);
  hData_Jet55_chMaxJtpt0p05->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p05->SetMarkerColor(5);
  hData_Jet55_chMaxJtpt0p05->Draw("same");
  line->Draw();

  cJet55->cd(2);
  hMC_Jet55_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet55_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_chMaxJtpt0p01->Draw();
  drawText("MC",0.2,0.8,14); 

  hMC_Jet55_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_chMaxJtpt0p03->Draw("same");

  hMC_Jet55_chMaxJtpt0p04->Rebin(5);
  hMC_Jet55_chMaxJtpt0p04->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p04->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p04->SetMarkerColor(4);
  hMC_Jet55_chMaxJtpt0p04->Draw("same");

  hMC_Jet55_chMaxJtpt0p05->Rebin(5);
  hMC_Jet55_chMaxJtpt0p05->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p05->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p05->SetMarkerColor(5);
  hMC_Jet55_chMaxJtpt0p05->Draw("same");
  line->Draw();

  cJet55->cd(3);
  drawText("#frac{charged Max}{Jet p_{T}} > X ",0.2,0.7,16);
  drawText("Jet 55 trigger (does not include higher triggers)", 0.2,0.8,16);
  TLegend *Jet55_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p01,"0.01","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p02,"0.02","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p03,"0.03","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p04,"0.04","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p05,"0.05","pl");
  Jet55_chMax->Draw();

  cJet55->cd(4);
  hData_Jet55_eMaxJtpt0p1->Rebin(5);
  hData_Jet55_eMaxJtpt0p1->Scale(1./5);
  hData_Jet55_eMaxJtpt0p1->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p1->SetMarkerColor(1);
  hData_Jet55_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hData_Jet55_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hData_Jet55_eMaxJtpt0p1->Draw();

  hData_Jet55_eMaxJtpt0p2->Rebin(5);
  hData_Jet55_eMaxJtpt0p2->Scale(1./5);
  hData_Jet55_eMaxJtpt0p2->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p2->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p2->Draw("same");

  hData_Jet55_eMaxJtpt0p3->Rebin(5);
  hData_Jet55_eMaxJtpt0p3->Scale(1./5);
  hData_Jet55_eMaxJtpt0p3->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p3->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p3->Draw("same");

  hData_Jet55_eMaxJtpt0p4->Rebin(5);
  hData_Jet55_eMaxJtpt0p4->Scale(1./5);
  hData_Jet55_eMaxJtpt0p4->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p4->SetMarkerColor(4);
  hData_Jet55_eMaxJtpt0p4->Draw("same");

  hData_Jet55_eMaxJtpt0p5->Rebin(5);
  hData_Jet55_eMaxJtpt0p5->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p5->SetMarkerColor(5);
  hData_Jet55_eMaxJtpt0p5->Draw("same");

  hData_Jet55_eMaxJtpt0p6->Rebin(5);
  hData_Jet55_eMaxJtpt0p6->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p6->SetMarkerColor(6);
  hData_Jet55_eMaxJtpt0p6->Draw("same");

  hData_Jet55_eMaxJtpt0p7->Rebin(5);
  hData_Jet55_eMaxJtpt0p7->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7->SetMarkerColor(7);
  hData_Jet55_eMaxJtpt0p7->Draw("same");

  hData_Jet55_eMaxJtpt0p8->Rebin(5);
  hData_Jet55_eMaxJtpt0p8->Scale(1./5);
  hData_Jet55_eMaxJtpt0p8->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p8->SetMarkerColor(8);
  hData_Jet55_eMaxJtpt0p8->Draw("same");

  hData_Jet55_eMaxJtpt0p9->Rebin(5);
  hData_Jet55_eMaxJtpt0p9->Scale(1./5);
  hData_Jet55_eMaxJtpt0p9->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p9->SetMarkerColor(9);
  hData_Jet55_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet55->cd(5);

  hMC_Jet55_eMaxJtpt0p1->Rebin(5);
  hMC_Jet55_eMaxJtpt0p1->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p1->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p1->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hMC_Jet55_eMaxJtpt0p1->Draw();

  hMC_Jet55_eMaxJtpt0p2->Rebin(5);
  hMC_Jet55_eMaxJtpt0p2->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p2->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p2->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p2->Draw("same");

  hMC_Jet55_eMaxJtpt0p3->Rebin(5);
  hMC_Jet55_eMaxJtpt0p3->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p3->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p3->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p3->Draw("same");

  hMC_Jet55_eMaxJtpt0p4->Rebin(5);
  hMC_Jet55_eMaxJtpt0p4->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p4->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p4->SetMarkerColor(4);
  hMC_Jet55_eMaxJtpt0p4->Draw("same");

  hMC_Jet55_eMaxJtpt0p5->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p5->SetMarkerColor(5);
  hMC_Jet55_eMaxJtpt0p5->Draw("same");

  hMC_Jet55_eMaxJtpt0p6->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p6->SetMarkerColor(6);
  hMC_Jet55_eMaxJtpt0p6->Draw("same");

  hMC_Jet55_eMaxJtpt0p7->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7->SetMarkerColor(7);
  hMC_Jet55_eMaxJtpt0p7->Draw("same");

  hMC_Jet55_eMaxJtpt0p8->Rebin(5);
  hMC_Jet55_eMaxJtpt0p8->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p8->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p8->SetMarkerColor(8);
  hMC_Jet55_eMaxJtpt0p8->Draw("same");

  hMC_Jet55_eMaxJtpt0p9->Rebin(5);
  hMC_Jet55_eMaxJtpt0p9->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p9->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p9->SetMarkerColor(9);
  hMC_Jet55_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet55->cd(6);
  drawText("#frac{electron Max}{Jet p_{T}} < X ",0.2,0.7,16);
  TLegend *Jet55_eMax = myLegend(0.2,0.2,0.6,0.6);
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p1," 0.1 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p2," 0.2 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p3," 0.3 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p4," 0.4 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p5," 0.5 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p6," 0.6 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p7," 0.7 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p8," 0.8 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p9," 0.9 ","pl");
  Jet55_eMax->Draw();

  cJet55->cd(7);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetXTitle("Jet p_{T} (GeV/c)");
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet55->cd(8);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetXTitle("Jet p_{T} (GeV/c)");
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet55->cd(9);
  drawText("#frac{electron Max}{Jet p_{T}} < X && #frac{charged Max}{Jet p_{T}} < Y",0.2,0.7,16);
  TLegend *Jet55_eMax_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01,"X = 7, Y = 1","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02,"X = 7, Y = 2","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03,"X = 7, Y = 3","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01,"X = 6, Y = 1","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02,"X = 6, Y = 2","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03,"X = 6, Y = 3","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01,"X = 5, Y = 1","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02,"X = 5, Y = 2","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03,"X = 5, Y = 3","pl");
  Jet55_eMax_chMax->Draw();

  cJet55->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jet55_spectra_ratio_jetid.pdf","RECREATE");



  TCanvas * cJet65 = new TCanvas("cJet65","",1200,1000);
  makeMultiPanelCanvas(cJet65,3,3,0.0,0.0,0.2,0.15,0.07);

  cJet65->cd(1);
  hData_Jet65_chMaxJtpt0p01->SetTitle(" ");
  hData_Jet65_chMaxJtpt0p01->SetXTitle(" Jet p_{T} (GeV/c)");
  hData_Jet65_chMaxJtpt0p01->SetYTitle(" Cut effiiency (with Cut/without Cut) ");
  hData_Jet65_chMaxJtpt0p01->SetAxisRange(0,2,"Y");
  hData_Jet65_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet65_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_chMaxJtpt0p01->Draw();

  hData_Jet65_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_chMaxJtpt0p02->Draw("same");

  hData_Jet65_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_chMaxJtpt0p03->Draw("same");

  hData_Jet65_chMaxJtpt0p04->Rebin(5);
  hData_Jet65_chMaxJtpt0p04->Scale(1./5);
  hData_Jet65_chMaxJtpt0p04->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p04->SetMarkerColor(4);
  hData_Jet65_chMaxJtpt0p04->Draw("same");

  hData_Jet65_chMaxJtpt0p05->Rebin(5);
  hData_Jet65_chMaxJtpt0p05->Scale(1./5);
  hData_Jet65_chMaxJtpt0p05->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p05->SetMarkerColor(5);
  hData_Jet65_chMaxJtpt0p05->Draw("same");
  line->Draw();
  drawText("Data",0.2,0.8,16);

  cJet65->cd(2);
  hMC_Jet65_chMaxJtpt0p01->SetAxisRange(0,2,"Y");
  hMC_Jet65_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet65_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_chMaxJtpt0p01->Draw();

  hMC_Jet65_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_chMaxJtpt0p03->Draw("same");

  hMC_Jet65_chMaxJtpt0p04->Rebin(5);
  hMC_Jet65_chMaxJtpt0p04->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p04->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p04->SetMarkerColor(4);
  hMC_Jet65_chMaxJtpt0p04->Draw("same");

  hMC_Jet65_chMaxJtpt0p05->Rebin(5);
  hMC_Jet65_chMaxJtpt0p05->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p05->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p05->SetMarkerColor(5);
  hMC_Jet65_chMaxJtpt0p05->Draw("same");
  line->Draw();

  cJet65->cd(3);
  drawText("Jet 65 Trigger",0.3,0.8,16);
  drawText("#frac{charged Max}{Jet p_{T}} > X ",0.2,0.7,16);
  TLegend *Jet65_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p01,"0.01","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p02,"0.02","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p03,"0.03","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p04,"0.04","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p05,"0.05","pl");
  Jet65_chMax->Draw();

  cJet65->cd(4);
  hData_Jet65_eMaxJtpt0p1->Rebin(5);
  hData_Jet65_eMaxJtpt0p1->Scale(1./5);
  hData_Jet65_eMaxJtpt0p1->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p1->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p1->Draw();

  hData_Jet65_eMaxJtpt0p2->Rebin(5);
  hData_Jet65_eMaxJtpt0p2->Scale(1./5);
  hData_Jet65_eMaxJtpt0p2->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p2->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p2->Draw("same");

  hData_Jet65_eMaxJtpt0p3->Rebin(5);
  hData_Jet65_eMaxJtpt0p3->Scale(1./5);
  hData_Jet65_eMaxJtpt0p3->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p3->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p3->Draw("same");

  hData_Jet65_eMaxJtpt0p4->Rebin(5);
  hData_Jet65_eMaxJtpt0p4->Scale(1./5);
  hData_Jet65_eMaxJtpt0p4->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p4->SetMarkerColor(4);
  hData_Jet65_eMaxJtpt0p4->Draw("same");

  hData_Jet65_eMaxJtpt0p5->Rebin(5);
  hData_Jet65_eMaxJtpt0p5->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p5->SetMarkerColor(5);
  hData_Jet65_eMaxJtpt0p5->Draw("same");

  hData_Jet65_eMaxJtpt0p6->Rebin(5);
  hData_Jet65_eMaxJtpt0p6->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p6->SetMarkerColor(6);
  hData_Jet65_eMaxJtpt0p6->Draw("same");

  hData_Jet65_eMaxJtpt0p7->Rebin(5);
  hData_Jet65_eMaxJtpt0p7->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7->SetMarkerColor(7);
  hData_Jet65_eMaxJtpt0p7->Draw("same");

  hData_Jet65_eMaxJtpt0p8->Rebin(5);
  hData_Jet65_eMaxJtpt0p8->Scale(1./5);
  hData_Jet65_eMaxJtpt0p8->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p8->SetMarkerColor(8);
  hData_Jet65_eMaxJtpt0p8->Draw("same");

  hData_Jet65_eMaxJtpt0p9->Rebin(5);
  hData_Jet65_eMaxJtpt0p9->Scale(1./5);
  hData_Jet65_eMaxJtpt0p9->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p9->SetMarkerColor(9);
  hData_Jet65_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet65->cd(5);

  hMC_Jet65_eMaxJtpt0p1->Rebin(5);
  hMC_Jet65_eMaxJtpt0p1->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p1->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p1->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p1->Draw();

  hMC_Jet65_eMaxJtpt0p2->Rebin(5);
  hMC_Jet65_eMaxJtpt0p2->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p2->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p2->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p2->Draw("same");

  hMC_Jet65_eMaxJtpt0p3->Rebin(5);
  hMC_Jet65_eMaxJtpt0p3->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p3->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p3->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p3->Draw("same");

  hMC_Jet65_eMaxJtpt0p4->Rebin(5);
  hMC_Jet65_eMaxJtpt0p4->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p4->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p4->SetMarkerColor(4);
  hMC_Jet65_eMaxJtpt0p4->Draw("same");

  hMC_Jet65_eMaxJtpt0p5->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p5->SetMarkerColor(5);
  hMC_Jet65_eMaxJtpt0p5->Draw("same");

  hMC_Jet65_eMaxJtpt0p6->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p6->SetMarkerColor(6);
  hMC_Jet65_eMaxJtpt0p6->Draw("same");

  hMC_Jet65_eMaxJtpt0p7->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7->SetMarkerColor(7);
  hMC_Jet65_eMaxJtpt0p7->Draw("same");

  hMC_Jet65_eMaxJtpt0p8->Rebin(5);
  hMC_Jet65_eMaxJtpt0p8->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p8->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p8->SetMarkerColor(8);
  hMC_Jet65_eMaxJtpt0p8->Draw("same");

  hMC_Jet65_eMaxJtpt0p9->Rebin(5);
  hMC_Jet65_eMaxJtpt0p9->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p9->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p9->SetMarkerColor(9);
  hMC_Jet65_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet65->cd(6);
  drawText("#frac{electron Max}{Jet p_{T}} < X ",0.2,0.7,16);
  TLegend *Jet65_eMax = myLegend(0.2,0.2,0.6,0.6);
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p1," 0.1 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p2," 0.2 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p3," 0.3 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p4," 0.4 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p5," 0.5 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p6," 0.6 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p7," 0.7 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p8," 0.8 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p9," 0.9 ","pl");
  Jet65_eMax->Draw();

  cJet65->cd(7);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet65->cd(8);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet65->cd(9);
  drawText("#frac{electron Max}{Jet p_{T}} < X && #frac{charged Max}{Jet p_{T}} < Y",0.2,0.7,16);
  TLegend *Jet65_eMax_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01,"X = 7, Y = 1","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02,"X = 7, Y = 2","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03,"X = 7, Y = 3","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01,"X = 6, Y = 1","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02,"X = 6, Y = 2","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03,"X = 6, Y = 3","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01,"X = 5, Y = 1","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02,"X = 5, Y = 2","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03,"X = 5, Y = 3","pl");
  Jet65_eMax_chMax->Draw();

  cJet65->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jet65_spectra_ratio_jetid.pdf","RECREATE");

  
  TCanvas * cJet80 = new TCanvas("cJet80","",1200,1000);
  makeMultiPanelCanvas(cJet80,3,3,0.0,0.0,0.2,0.15,0.07);

  cJet80->cd(1);
  hData_Jet80_chMaxJtpt0p01->SetTitle(" ");
  hData_Jet80_chMaxJtpt0p01->SetXTitle(" Jet p_{T} (GeV/c)");
  hData_Jet80_chMaxJtpt0p01->SetYTitle(" Cut effiiency (with Cut/without Cut) ");
  hData_Jet80_chMaxJtpt0p01->SetAxisRange(0,2,"Y");
  hData_Jet80_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet80_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_chMaxJtpt0p01->Draw();

  hData_Jet80_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_chMaxJtpt0p02->Draw("same");

  hData_Jet80_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_chMaxJtpt0p03->Draw("same");

  hData_Jet80_chMaxJtpt0p04->Rebin(5);
  hData_Jet80_chMaxJtpt0p04->Scale(1./5);
  hData_Jet80_chMaxJtpt0p04->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p04->SetMarkerColor(4);
  hData_Jet80_chMaxJtpt0p04->Draw("same");

  hData_Jet80_chMaxJtpt0p05->Rebin(5);
  hData_Jet80_chMaxJtpt0p05->Scale(1./5);
  hData_Jet80_chMaxJtpt0p05->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p05->SetMarkerColor(5);
  hData_Jet80_chMaxJtpt0p05->Draw("same");
  line->Draw();

  cJet80->cd(2);
  hMC_Jet80_chMaxJtpt0p01->SetAxisRange(0,2,"Y");
  hMC_Jet80_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet80_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_chMaxJtpt0p01->Draw();

  hMC_Jet80_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_chMaxJtpt0p03->Draw("same");

  hMC_Jet80_chMaxJtpt0p04->Rebin(5);
  hMC_Jet80_chMaxJtpt0p04->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p04->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p04->SetMarkerColor(4);
  hMC_Jet80_chMaxJtpt0p04->Draw("same");

  hMC_Jet80_chMaxJtpt0p05->Rebin(5);
  hMC_Jet80_chMaxJtpt0p05->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p05->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p05->SetMarkerColor(5);
  hMC_Jet80_chMaxJtpt0p05->Draw("same");
  line->Draw();

  cJet80->cd(3);
  drawText("#frac{charged Max}{Jet p_{T}} > X ",0.2,0.7,16);
  TLegend *Jet80_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p01,"0.01","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p02,"0.02","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p03,"0.03","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p04,"0.04","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p05,"0.05","pl");
  Jet80_chMax->Draw();

  cJet80->cd(4);
  hData_Jet80_eMaxJtpt0p1->Rebin(5);
  hData_Jet80_eMaxJtpt0p1->Scale(1./5);
  hData_Jet80_eMaxJtpt0p1->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p1->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p1->Draw();

  hData_Jet80_eMaxJtpt0p2->Rebin(5);
  hData_Jet80_eMaxJtpt0p2->Scale(1./5);
  hData_Jet80_eMaxJtpt0p2->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p2->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p2->Draw("same");

  hData_Jet80_eMaxJtpt0p3->Rebin(5);
  hData_Jet80_eMaxJtpt0p3->Scale(1./5);
  hData_Jet80_eMaxJtpt0p3->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p3->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p3->Draw("same");

  hData_Jet80_eMaxJtpt0p4->Rebin(5);
  hData_Jet80_eMaxJtpt0p4->Scale(1./5);
  hData_Jet80_eMaxJtpt0p4->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p4->SetMarkerColor(4);
  hData_Jet80_eMaxJtpt0p4->Draw("same");

  hData_Jet80_eMaxJtpt0p5->Rebin(5);
  hData_Jet80_eMaxJtpt0p5->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p5->SetMarkerColor(5);
  hData_Jet80_eMaxJtpt0p5->Draw("same");

  hData_Jet80_eMaxJtpt0p6->Rebin(5);
  hData_Jet80_eMaxJtpt0p6->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p6->SetMarkerColor(6);
  hData_Jet80_eMaxJtpt0p6->Draw("same");

  hData_Jet80_eMaxJtpt0p7->Rebin(5);
  hData_Jet80_eMaxJtpt0p7->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7->SetMarkerColor(7);
  hData_Jet80_eMaxJtpt0p7->Draw("same");

  hData_Jet80_eMaxJtpt0p8->Rebin(5);
  hData_Jet80_eMaxJtpt0p8->Scale(1./5);
  hData_Jet80_eMaxJtpt0p8->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p8->SetMarkerColor(8);
  hData_Jet80_eMaxJtpt0p8->Draw("same");

  hData_Jet80_eMaxJtpt0p9->Rebin(5);
  hData_Jet80_eMaxJtpt0p9->Scale(1./5);
  hData_Jet80_eMaxJtpt0p9->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p9->SetMarkerColor(9);
  hData_Jet80_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet80->cd(5);

  hMC_Jet80_eMaxJtpt0p1->Rebin(5);
  hMC_Jet80_eMaxJtpt0p1->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p1->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p1->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p1->Draw();

  hMC_Jet80_eMaxJtpt0p2->Rebin(5);
  hMC_Jet80_eMaxJtpt0p2->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p2->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p2->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p2->Draw("same");

  hMC_Jet80_eMaxJtpt0p3->Rebin(5);
  hMC_Jet80_eMaxJtpt0p3->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p3->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p3->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p3->Draw("same");

  hMC_Jet80_eMaxJtpt0p4->Rebin(5);
  hMC_Jet80_eMaxJtpt0p4->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p4->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p4->SetMarkerColor(4);
  hMC_Jet80_eMaxJtpt0p4->Draw("same");

  hMC_Jet80_eMaxJtpt0p5->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p5->SetMarkerColor(5);
  hMC_Jet80_eMaxJtpt0p5->Draw("same");

  hMC_Jet80_eMaxJtpt0p6->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p6->SetMarkerColor(6);
  hMC_Jet80_eMaxJtpt0p6->Draw("same");

  hMC_Jet80_eMaxJtpt0p7->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7->SetMarkerColor(7);
  hMC_Jet80_eMaxJtpt0p7->Draw("same");

  hMC_Jet80_eMaxJtpt0p8->Rebin(5);
  hMC_Jet80_eMaxJtpt0p8->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p8->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p8->SetMarkerColor(8);
  hMC_Jet80_eMaxJtpt0p8->Draw("same");

  hMC_Jet80_eMaxJtpt0p9->Rebin(5);
  hMC_Jet80_eMaxJtpt0p9->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p9->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p9->SetMarkerColor(9);
  hMC_Jet80_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet80->cd(6);
  drawText("#frac{electron Max}{Jet p_{T}} < X ",0.2,0.7,16);
  TLegend *Jet80_eMax = myLegend(0.2,0.2,0.6,0.6);
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p1," 0.1 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p2," 0.2 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p3," 0.3 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p4," 0.4 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p5," 0.5 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p6," 0.6 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p7," 0.7 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p8," 0.8 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p9," 0.9 ","pl");
  Jet80_eMax->Draw();

  cJet80->cd(7);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet80->cd(8);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet80->cd(9);
  drawText("#frac{electron Max}{Jet p_{T}} < X && #frac{charged Max}{Jet p_{T}} > Y",0.2,0.7,16);
  TLegend *Jet80_eMax_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01,"X = 7, Y = 1","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02,"X = 7, Y = 2","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03,"X = 7, Y = 3","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01,"X = 6, Y = 1","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02,"X = 6, Y = 2","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03,"X = 6, Y = 3","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01,"X = 5, Y = 1","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02,"X = 5, Y = 2","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03,"X = 5, Y = 3","pl");
  Jet80_eMax_chMax->Draw();

  cJet80->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jet80_spectra_ratio_jetid.pdf","RECREATE");



  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
