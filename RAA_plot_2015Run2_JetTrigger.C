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


using namespace std;

void RAA_plot_2015Run2_JetTrigger(int radius = 3, char * algo= "Pu", char * jet_type = "Calo"){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // get the scaling factor from the files at 2.76 TeV and 5.02 TeV pythia+hydjet pthat=30
  TFile *f5p02 = TFile::Open("/mnt/hadoop/cms/store/user/luck/L1Emulator/PyquenUnquenched_DiJet_pt30_PbPb_5020GeV_actuallyEmbedded_HiForest.root");
  TFile *f2p76 = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat30_Track9_Jet30_matchEqR_merged_forest_0.root");

  // get the respective trees.

  TTree * jet2p76 = (TTree*)f2p76->Get("akPu3CaloJetAnalyzer/t");
  TTree * evt2p76 = (TTree*)f2p76->Get("hiEvtAnalyzer/HiTree");
  TTree * hlt2p76 = (TTree*)f2p76->Get("hltanalysis/HltTree");
  TTree * skm2p76 = (TTree*)f2p76->Get("skimanalysis/HltTree");

  jet2p76->AddFriend(evt2p76);
  jet2p76->AddFriend(hlt2p76);
  jet2p76->AddFriend(skm2p76); 
  
  TTree * jet5p02 = (TTree*)f5p02->Get("akPu3CaloJetAnalyzer/t");
  TTree * evt5p02 = (TTree*)f5p02->Get("hiEvtAnalyzer/HiTree");
  TTree * hlt5p02 = (TTree*)f5p02->Get("hltanalysis/HltTree");
  TTree * skm5p02 = (TTree*)f5p02->Get("skimanalysis/HltTree");

  jet5p02->AddFriend(evt5p02);
  jet5p02->AddFriend(hlt5p02);
  jet5p02->AddFriend(skm5p02);

  // declare the histograms to get the scale factor 
  TH1F * h2p76 = new TH1F("h2p76","",200,0,400);
  TH1F * h5p02 = new TH1F("h5p02","",200,0,400);
  
  jet2p76->Draw("jtpt>>h2p76","pcollisionEventSelection && vz<15 && vz>-15 && jteta<2 && jteta>-2","goff");
  Float_t nEvents_2p76 = jet2p76->GetEntries("pcollisionEventSelection && vz<15 && vz>-15 && jteta<2 && jteta>-2");
  Float_t nJets_2p76 = h2p76->GetEntries();
  cout<<"2p76 Number of events = "<<nEvents_2p76<<endl;
  cout<<"2p76 Number of Jets   = "<<nJets_2p76<<endl;
  h2p76->Scale(1./nEvents_2p76);
  
  jet5p02->Draw("jtpt>>h5p02","pcollisionEventSelection && vz<15 && vz>-15 && jteta<2 && jteta>-2","goff");
  Float_t nEvents_5p02 = jet5p02->GetEntries("pcollisionEventSelection && vz<15 && vz>-15 && jteta<2 && jteta>-2");
  Float_t nJets_5p02 = h5p02->GetEntries();
  cout<<"5p02 Number of events = "<<nEvents_5p02<<endl;
  cout<<"5p02 Number of Jets   = "<<nJets_5p02<<endl;
  h5p02->Scale(1./nEvents_5p02);

  Float_t ScaleFactor = (Float_t) nEvents_5p02/nEvents_2p76;
  cout<<" Scale factor in terms of number of events = "<<ScaleFactor<<endl;
  Float_t ScaleFactor_Jets = (Float_t) nJets_5p02/nJets_2p76;
  cout<<" Scale factor in terms of number of Jets   = "<<ScaleFactor_Jets<<endl;
  
  TH1F * hRatio = (TH1F*)h5p02->Clone("hRatio");
  hRatio->Divide(h2p76);

  TCanvas * cRatio = new TCanvas("cRatio","",800,600);
  makeHistTitle(hRatio," "," akPu3Calo Jet p_{T} (GeV/c)"," #frac{5.02 TeV #hat{p_{T}}=30}{2.76 TeV #hat{p_{T}}=30}");
  hRatio->Draw();
  cRatio->SaveAs("../../Plots/scalefactor_2p75_5p02_ratio_jetpT_pthat30.pdf","RECREATE");
    


}


