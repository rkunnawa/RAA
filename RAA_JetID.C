// Raghav Kunnawalkam Elayavalli
// Nov 18th 2014
// Rutgers 
// raghav.k.e at CERN dot CH 

// Jet ID macro - To study the issue of fake jets in the Jet Reconstruction for the RAA analysis. First idea is that the electron reconstructed pT is off. Going to use analysis cuts : neMax/(neMax+chMax+phMax)<0.975 and muMax/(neMax+chMax+phMax)<0.975 and testing out this cut in Data and MC: eMax/(neMax+chMax+phMax)<0.975. 

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

static const int nbins_eta = 15;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0}, {-2.0,+2.0}, {-3.0,+3.0},
  {-3.0,-2.5}, {-2.5,-2.0}, {-2.0,-1.5}, 
  {-1.5,-1.0}, {-1.0,-0.5}, {-0.5,0}, {0,+0.5}, 
  {+0.5,+1.0}, {+1.0,+1.5}, {+1.5,+2.0}, 
  {+2.0,+2.5}, {+2.5,+3.0}
};

static const double delta_eta[nbins_eta] = {
  2.0, 4.0, 6.0, 
  0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 
  0.5, 0.5
};

static const char etaWidth [nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20","n30_eta_p30",
  "n30_eta_n25","n25_eta_n20","n20_eta_n15",
  "n15_eta_n10","n10_eta_n05","n05_eta_0","0_eta_p05",
  "p05_eta_p10","p10_eta_p15","p15_eta_p20",
  "p20_eta_p25","p25_eta_p30"
};

static const int no_radius = 7; 
static const int list_radius[no_radius] = {1,2,3,4,5,6,7};

// divide by bin width
void divideBinWidth(TH1 *h)
{
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    if(val!=0){
      val/=h->GetBinWidth(i);
      valErr/=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }  
  }

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

using namespace std;


void RAA_JetID(int radius = 3, char *algo = "Pu", char *jet_type = "PF"){
  
  TH1::SetDefaultSumw2();

  TStopwatch timer;
  timer.Start();

  TDatime date;

  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = true;

  // Get the input files for Data and MC including the trees that are necessary 
  TFile *fDatain = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_data_ak%s%s_testComb4_cut1_20141114.root",algo,jet_type));
  TFile *fMCin = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_pp_mc_nocut_ak%s%s_20141118.root",algo,jet_type));

  TTree *jetData = (TTree*)fDatain->Get("jets_ID");
  TTree *jetMC = (TTree*)fMCin->Get("jets_ID");

  //check if the trees are filled. 
  if(printDebug)cout<<"jetData entries = "<<jetData->GetEntries()<<endl;
  if(printDebug)cout<<"jetMC entries   = "<<jetMC->GetEntries()<<endl;

  // lets think about what we want to plot here: 
  // i want to plot the ratio of jet spectra (from the specific triggers Jet55 or Jet65 or Jet80) with cut over without cut. 
  // both the trees already have the basic event selection here: pHBHENoiseFilter (Data), pcollisionEventSelection (MC and Data) and |vz|<15 and |eta|<2 for both data and MC. 
  // 

  TH1F *hMC_Jet55_withCuts = new TH1F("hMC_Jet55_withCuts","",1000,0,1000);

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}
