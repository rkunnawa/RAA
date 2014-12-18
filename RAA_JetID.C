// Raghav Kunnawalkam Elayavalli
// Nov 18th 2014
// Rutgers 
// raghav.k.e at CERN dot CH 

// Jet ID macro - To study the issue of fake jets in the Jet Reconstruction for the RAA analysis. First idea is that the electron reconstructed pT is off. Going to use analysis cuts : neMax/(neMax+chMax+phMax)<0.975 and muMax/(neMax+chMax+phMax)<0.975 and testing out this cut in Data and MC: eMax/(neMax+chMax+phMax)<0.975. 

// Dec 9th making the plots for PU since RAA is finally going to PU Jets 
// Dec 10th change the value of the cuts to reflect + adding new histograms 

// Dec 11th - first include the plots for individual pfcandidate variables like chMax, chSum etc... and also make the 2D scatter plot of their ratio with the pT vs pT. for the individual candidates.  

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

static const int TrigValue = 3;
static const int CutValue = 2;
static const char TrigName [TrigValue][256] = {"HLT55","HLT65","HLT80"};
static const char isJetID [CutValue][256] = {"without","with"};

static const int nbins_cent = 6;
static const double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
static const char centWidth[nbins_cent+1][256] = {"0 - 5 %","5 - 10 %","10 - 30 %","30 - 50 %","50 - 70 %","70 - 90 %","0 - 100 %"};

static const double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};

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
  TFile *fDatain = TFile::Open(Form("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/PbPb_jetntuple_withEvtCuts_SuperNovaRejected_ak%s%d%s_20141209.root",algo,radius,jet_type));
  TFile *fMCin = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_mc_nocut_ak%s%s_20141210.root",algo,jet_type));

  TTree *jetData = (TTree*)fDatain->Get("jets_ID");
  TTree *jetMC = (TTree*)fMCin->Get("jets_ID");

  //check if the trees are filled. 
  if(printDebug)cout<<"jetData no of jets = "<<jetData->GetEntries()<<endl;
  if(printDebug)cout<<"jetMC no of jets  = "<<jetMC->GetEntries()<<endl;

  // lets think about what we want to plot here: 
  // i want to plot the ratio of jet spectra (from the specific triggers Jet55 or Jet65 or Jet80) with cut over without cut. 
  // both the trees already have the basic event selection here: pHBHENoiseFilter (Data), pcollisionEventSelection (MC and Data) and |vz|<15 and |eta|<2 for both data and MC. 
  // 

  //TH1F *hMC_Jet55_withCuts = new TH1F("hMC_Jet55_withCuts","",1000,0,1000);

  // setup of data histograms: 2D arrays, [JetTrigger][isCut]
  // JetTrigger is either 0 - Jet55; 1 - Jet65; 2 - Jet80. 
  // isJetID is either 0: No Jet ID cuts, 1: Yes Jet ID Cuts
  // necessary Jet ID cuts: 

  TCut Jet55 = "jet55 && l1sj36 && !jet65 && !jet80";
  TCut Jet65 = "jet65 && l1sj36 && !jet80";
  TCut Jet80 = "jet80 && l1sj52";
  TCut elRjc = "eMax/jtpt<0.3";
  TCut neRjc = "neMax/(chMax+neMax+phMax)<0.9";
  TCut phRjc = "phMax/(chMax+neMax+phMax)<0.9";
  TCut chRjc = "chMax/(chMax+neMax+phMax)<0.9";
  TCut chID  = "chMax/jtpt>0.05";
  TCut muRjc = "muMax/(chMax+neMax+phMax)<0.9";
  //TCut cent[nbins_cent] = {"0","1","2","3","4","5"};

  //Get the Spectra in the centrality bins as well, which would mean that we have to add in the centrality loop

  TH1F *hData[3][2][nbins_cent+1];
  TH1F *hData_Ratio[TrigValue][nbins_cent+1];
  TH1F *hMC[3][2][nbins_cent+1];
  TH1F *hMC_Ratio[TrigValue][nbins_cent+1];

  // need the ratio histograms from each time we add the Jet ID cuts for the full spectra. 
  // need to make 2d histograms like what Yetkin showed plotting the ratio of species fraction/jtpt vs jtpt and for genpt for MC - done
  
  TH1F *hData_chMax[nbins_cent+1], *hData_chSum[nbins_cent+1],*hData_neMax[nbins_cent+1], *hData_neSum[nbins_cent+1],*hData_phMax[nbins_cent+1], *hData_phSum[nbins_cent+1],*hData_eMax[nbins_cent+1], *hData_eSum[nbins_cent+1],*hData_muMax[nbins_cent+1], *hData_muSum[nbins_cent+1];

  TH1F *hMC_chMax[nbins_cent+1], *hMC_chSum[nbins_cent+1],*hMC_neMax[nbins_cent+1], *hMC_neSum[nbins_cent+1],*hMC_phMax[nbins_cent+1], *hMC_phSum[nbins_cent+1],*hMC_eMax[nbins_cent+1], *hMC_eSum[nbins_cent+1],*hMC_muMax[nbins_cent+1], *hMC_muSum[nbins_cent+1];
  
  TH2F *hData_chMaxJtPt_Pt[nbins_cent+1], *hData_chSumJtPt_Pt[nbins_cent+1], *hData_neMaxJtPt_Pt[nbins_cent+1], *hData_neSumJtPt_Pt[nbins_cent+1], *hData_phMaxJtPt_Pt[nbins_cent+1], *hData_phSumJtPt_Pt[nbins_cent+1], *hData_eMaxJtPt_Pt[nbins_cent+1], *hData_eSumJtPt_Pt[nbins_cent+1], *hData_muMaxJtPt_Pt[nbins_cent+1], *hData_muSumJtPt_Pt[nbins_cent+1];
  
  TH2F *hMC_chMaxJtPt_Pt[nbins_cent+1], *hMC_chSumJtPt_Pt[nbins_cent+1], *hMC_neMaxJtPt_Pt[nbins_cent+1], *hMC_neSumJtPt_Pt[nbins_cent+1], *hMC_phMaxJtPt_Pt[nbins_cent+1], *hMC_phSumJtPt_Pt[nbins_cent+1], *hMC_eMaxJtPt_Pt[nbins_cent+1], *hMC_eSumJtPt_Pt[nbins_cent+1], *hMC_muMaxJtPt_Pt[nbins_cent+1], *hMC_muSumJtPt_Pt[nbins_cent+1];

  TH2F *hMC_chMaxGenPt_Pt[nbins_cent+1], *hMC_chSumGenPt_Pt[nbins_cent+1], *hMC_neMaxGenPt_Pt[nbins_cent+1], *hMC_neSumGenPt_Pt[nbins_cent+1], *hMC_phMaxGenPt_Pt[nbins_cent+1], *hMC_phSumGenPt_Pt[nbins_cent+1], *hMC_eMaxGenPt_Pt[nbins_cent+1], *hMC_eSumGenPt_Pt[nbins_cent+1], *hMC_muMaxGenPt_Pt[nbins_cent+1], *hMC_muSumGenPt_Pt[nbins_cent+1];

  TH2F *hData_chMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_chSumSumMaxChNePh_Pt[nbins_cent+1], *hData_neMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_neSumSumMaxChNePh_Pt[nbins_cent+1], *hData_phMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_phSumSumMaxChNePh_Pt[nbins_cent+1], *hData_eMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_eSumSumMaxChNePh_Pt[nbins_cent+1], *hData_muMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_muSumSumMaxChNePh_Pt[nbins_cent+1];
  
  TH2F *hMC_chMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_chSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_neMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_neSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_phMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_phSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_eMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_eSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_muMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_muSumSumMaxChNePh_Pt[nbins_cent+1];

  TH2F *hMC_chMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_chSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_neMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_neSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_phMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_phSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_eMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_eSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_muMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_muSumSumMaxChNePh_GenPt[nbins_cent+1];

#if 0
  TH2F *hData_chMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_chSumSumMaxChNePh_Pt[nbins_cent+1], *hData_neMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_neSumSumMaxChNePh_Pt[nbins_cent+1], *hData_phMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_phSumSumMaxChNePh_Pt[nbins_cent+1], *hData_eMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_eSumSumMaxChNePh_Pt[nbins_cent+1], *hData_muMaxSumMaxChNePh_Pt[nbins_cent+1], *hData_muSumSumMaxChNePh_Pt[nbins_cent+1];
  
  TH2F *hMC_chMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_chSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_neMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_neSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_phMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_phSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_eMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_eSumSumMaxChNePh_Pt[nbins_cent+1], *hMC_muMaxSumMaxChNePh_Pt[nbins_cent+1], *hMC_muSumSumMaxChNePh_Pt[nbins_cent+1];

  TH2F *hMC_chMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_chSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_neMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_neSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_phMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_phSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_eMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_eSumSumMaxChNePh_GenPt[nbins_cent+1], *hMC_muMaxSumMaxChNePh_GenPt[nbins_cent+1], *hMC_muSumSumMaxChNePh_GenPt[nbins_cent+1];
#endif

  //We are also going to check for the individual application of the cuts rather than the sequencial application. for both data and MC. and we are also going to need the ratio plots 


  TH1F *hData_chMax_SumMaxChNePh[CutValue][nbins_cent+1], *hData_neMax_SumMaxChNePh[CutValue][nbins_cent+1], *hData_phMax_SumMaxChNePh[CutValue][nbins_cent+1], *hData_eMax_SumMaxChNePh[CutValue][nbins_cent+1], *hData_muMax_SumMaxChNePh[CutValue][nbins_cent+1]; 

  TH1F *hMC_chMax_SumMaxChNePh[CutValue][nbins_cent+1], *hMC_neMax_SumMaxChNePh[CutValue][nbins_cent+1], *hMC_phMax_SumMaxChNePh[CutValue][nbins_cent+1], *hMC_eMax_SumMaxChNePh[CutValue][nbins_cent+1], *hMC_muMax_SumMaxChNePh[CutValue][nbins_cent+1]; 

  TH1F *hData_Ratio_chMax_SumMaxChNePh[nbins_cent+1], *hData_Ratio_neMax_SumMaxChNePh[nbins_cent+1], *hData_Ratio_phMax_SumMaxChNePh[nbins_cent+1], *hData_Ratio_muMax_SumMaxChNePh[nbins_cent+1], *hData_Ratio_eMax_SumMaxChNePh[nbins_cent+1];

  TH1F *hMC_Ratio_chMax_SumMaxChNePh[nbins_cent+1], *hMC_Ratio_neMax_SumMaxChNePh[nbins_cent+1], *hMC_Ratio_phMax_SumMaxChNePh[nbins_cent+1], *hMC_Ratio_muMax_SumMaxChNePh[nbins_cent+1], *hMC_Ratio_eMax_SumMaxChNePh[nbins_cent+1];

  TH1F *hData_chMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hData_neMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hData_phMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hData_eMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hData_muMax_SumSumChPhNeEMu[CutValue][nbins_cent+1]; 

  TH1F *hMC_chMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hMC_neMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hMC_phMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hMC_eMax_SumSumChPhNeEMu[CutValue][nbins_cent+1], *hMC_muMax_SumSumChPhNeEMu[CutValue][nbins_cent+1]; 

  TH1F *hData_Ratio_chMax_SumSumChPhNeEMu[nbins_cent+1], *hData_Ratio_neMax_SumSumChPhNeEMu[nbins_cent+1], *hData_Ratio_phMax_SumSumChPhNeEMu[nbins_cent+1], *hData_Ratio_muMax_SumSumChPhNeEMu[nbins_cent+1], *hData_Ratio_eMax_SumSumChPhNeEMu[nbins_cent+1];

  TH1F *hMC_Ratio_chMax_SumSumChPhNeEMu[nbins_cent+1], *hMC_Ratio_neMax_SumSumChPhNeEMu[nbins_cent+1], *hMC_Ratio_phMax_SumSumChPhNeEMu[nbins_cent+1], *hMC_Ratio_muMax_SumSumChPhNeEMu[nbins_cent+1], *hMC_Ratio_eMax_SumSumChPhNeEMu[nbins_cent+1];

  TH1F * hData_eMaxJtPt_0p9[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p8[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p7[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p6[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p5[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p4[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p3[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p2[nbins_cent+1];
  TH1F * hData_eMaxJtPt_0p1[nbins_cent+1];

  TH1F * hData_noCut[nbins_cent+1];

  TH1F * hData_Ratio_eMaxJtPt_0p9[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p8[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p7[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p6[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p5[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p4[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p3[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p2[nbins_cent+1];
  TH1F * hData_Ratio_eMaxJtPt_0p1[nbins_cent+1];

  TH1F * hMC_eMaxJtPt_0p9[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p8[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p7[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p6[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p5[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p4[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p3[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p2[nbins_cent+1];
  TH1F * hMC_eMaxJtPt_0p1[nbins_cent+1];

  TH1F * hMC_noCut[nbins_cent+1];

  TH1F * hMC_Ratio_eMaxJtPt_0p9[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p8[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p7[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p6[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p5[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p4[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p3[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p2[nbins_cent+1];
  TH1F * hMC_Ratio_eMaxJtPt_0p1[nbins_cent+1];



  for(int i = 0;i<=nbins_cent;i++){

    hData_chMax[i] = new TH1F(Form("hData_chMax_cent%d",i),"",1000,0,1000);
    hData_chSum[i] = new TH1F(Form("hData_chSum_cent%d",i),"",1000,0,1000);
    hData_neMax[i] = new TH1F(Form("hData_neMax_cent%d",i),"",1000,0,1000);
    hData_neSum[i] = new TH1F(Form("hData_neSum_cent%d",i),"",1000,0,1000);
    hData_phMax[i] = new TH1F(Form("hData_phMax_cent%d",i),"",1000,0,1000);
    hData_phSum[i] = new TH1F(Form("hData_phSum_cent%d",i),"",1000,0,1000);
    hData_eMax[i]  = new TH1F(Form("hData_eMax_cent%d",i),"",1000,0,1000);
    hData_eSum[i]  = new TH1F(Form("hData_eSum_cent%d",i),"",1000,0,1000);
    hData_muMax[i] = new TH1F(Form("hData_muMax_cent%d",i),"",1000,0,1000);
    hData_muSum[i] = new TH1F(Form("hData_muSum_cent%d",i),"",1000,0,1000);

    hMC_chMax[i] = new TH1F(Form("hMC_chMax_cent%d",i),"",1000,0,1000);
    hMC_chSum[i] = new TH1F(Form("hMC_chSum_cent%d",i),"",1000,0,1000);
    hMC_neMax[i] = new TH1F(Form("hMC_neMax_cent%d",i),"",1000,0,1000);
    hMC_neSum[i] = new TH1F(Form("hMC_neSum_cent%d",i),"",1000,0,1000);
    hMC_phMax[i] = new TH1F(Form("hMC_phMax_cent%d",i),"",1000,0,1000);
    hMC_phSum[i] = new TH1F(Form("hMC_phSum_cent%d",i),"",1000,0,1000);
    hMC_eMax[i]  = new TH1F(Form("hMC_eMax_cent%d",i),"",1000,0,1000);
    hMC_eSum[i]  = new TH1F(Form("hMC_eSum_cent%d",i),"",1000,0,1000);
    hMC_muMax[i] = new TH1F(Form("hMC_muMax_cent%d",i),"",1000,0,1000);
    hMC_muSum[i] = new TH1F(Form("hMC_muSum_cent%d",i),"",1000,0,1000);    

    hData_chMaxJtPt_Pt[i] = new TH2F(Form("hData_chMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_chSumJtPt_Pt[i] = new TH2F(Form("hData_chSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_phMaxJtPt_Pt[i] = new TH2F(Form("hData_phMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_phSumJtPt_Pt[i] = new TH2F(Form("hData_phSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_neMaxJtPt_Pt[i] = new TH2F(Form("hData_neMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_neSumJtPt_Pt[i] = new TH2F(Form("hData_neSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_eMaxJtPt_Pt[i]  = new TH2F(Form("hData_eMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_eSumJtPt_Pt[i]  = new TH2F(Form("hData_eSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_muMaxJtPt_Pt[i] = new TH2F(Form("hData_muMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_muSumJtPt_Pt[i] = new TH2F(Form("hData_muSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);

    hMC_chMaxJtPt_Pt[i] = new TH2F(Form("hMC_chMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_chSumJtPt_Pt[i] = new TH2F(Form("hMC_chSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phMaxJtPt_Pt[i] = new TH2F(Form("hMC_phMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phSumJtPt_Pt[i] = new TH2F(Form("hMC_phSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neMaxJtPt_Pt[i] = new TH2F(Form("hMC_neMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neSumJtPt_Pt[i] = new TH2F(Form("hMC_neSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eMaxJtPt_Pt[i]  = new TH2F(Form("hMC_eMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eSumJtPt_Pt[i]  = new TH2F(Form("hMC_eSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muMaxJtPt_Pt[i] = new TH2F(Form("hMC_muMaxJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muSumJtPt_Pt[i] = new TH2F(Form("hMC_muSumJtPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);

    hMC_chMaxGenPt_Pt[i] = new TH2F(Form("hMC_chMaxGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_chSumGenPt_Pt[i] = new TH2F(Form("hMC_chSumGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phMaxGenPt_Pt[i] = new TH2F(Form("hMC_phMaxGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phSumGenPt_Pt[i] = new TH2F(Form("hMC_phSumGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neMaxGenPt_Pt[i] = new TH2F(Form("hMC_neMaxGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neSumGenPt_Pt[i] = new TH2F(Form("hMC_neSumGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eMaxGenPt_Pt[i]  = new TH2F(Form("hMC_eMaxGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eSumGenPt_Pt[i]  = new TH2F(Form("hMC_eSumGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muMaxGenPt_Pt[i] = new TH2F(Form("hMC_muMaxGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muSumGenPt_Pt[i] = new TH2F(Form("hMC_muSumGenPt_Pt_cent%d",i),"",100,0,10,1000,0,1000);

    hData_chMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hData_chMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_chSumSumMaxChNePh_Pt[i] = new TH2F(Form("hData_chSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_phMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hData_phMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_phSumSumMaxChNePh_Pt[i] = new TH2F(Form("hData_phSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_neMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hData_neMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_neSumSumMaxChNePh_Pt[i] = new TH2F(Form("hData_neSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_eMaxSumMaxChNePh_Pt[i]  = new TH2F(Form("hData_eMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_eSumSumMaxChNePh_Pt[i]  = new TH2F(Form("hData_eSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_muMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hData_muMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hData_muSumSumMaxChNePh_Pt[i] = new TH2F(Form("hData_muSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);

    hMC_chMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_chMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_chSumSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_chSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_phMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phSumSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_phSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_neMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neSumSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_neSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eMaxSumMaxChNePh_Pt[i]  = new TH2F(Form("hMC_eMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eSumSumMaxChNePh_Pt[i]  = new TH2F(Form("hMC_eSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muMaxSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_muMaxSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muSumSumMaxChNePh_Pt[i] = new TH2F(Form("hMC_muSumSumMaxChNePh_Pt_cent%d",i),"",100,0,10,1000,0,1000);

    hMC_chMaxSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_chMaxSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_chSumSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_chSumSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phMaxSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_phMaxSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_phSumSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_phSumSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neMaxSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_neMaxSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_neSumSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_neSumSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eMaxSumMaxChNePh_GenPt[i]  = new TH2F(Form("hMC_eMaxSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_eSumSumMaxChNePh_GenPt[i]  = new TH2F(Form("hMC_eSumSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muMaxSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_muMaxSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);
    hMC_muSumSumMaxChNePh_GenPt[i] = new TH2F(Form("hMC_muSumSumMaxChNePh_GenPt_cent%d",i),"",100,0,10,1000,0,1000);

    for(int a = 0;a<TrigValue;a++){

      for(int b = 0;b<CutValue;b++){

        hData[a][b][i] = new TH1F(Form("hData_%s_%s_JetID_cent%d",TrigName[a],isJetID[b],i),Form("Data Jet Spectra from %s trigger %s JetID cuts in the centrality bin %s",TrigName[a],isJetID[b],centWidth[i]),1000,0,1000);
	hMC[a][b][i] = new TH1F(Form("hMC_%s_%s_JetID_cent%d",TrigName[a],isJetID[b],i),Form("MC Jet Spectra from %s trigger %s JetID cuts in the centrality bin %s",TrigName[a],isJetID[b],centWidth[i]),1000,0,1000);

      }
      
    }
    
  }

  for(int i = 0;i<nbins_cent+1;i++){

    for(int b = 0;b<CutValue;b++){

      hData_chMax_SumMaxChNePh[b][i] = new TH1F(Form("hData_chMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hMC_chMax_SumMaxChNePh[b][i] = new TH1F(Form("hMC_chMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hData_phMax_SumMaxChNePh[b][i] = new TH1F(Form("hData_phMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hMC_phMax_SumMaxChNePh[b][i] = new TH1F(Form("hMC_phMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hData_neMax_SumMaxChNePh[b][i] = new TH1F(Form("hData_neMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hMC_neMax_SumMaxChNePh[b][i] = new TH1F(Form("hMC_neMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hData_muMax_SumMaxChNePh[b][i] = new TH1F(Form("hData_muMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hMC_muMax_SumMaxChNePh[b][i] = new TH1F(Form("hMC_muMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hData_eMax_SumMaxChNePh[b][i] = new TH1F(Form("hData_eMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);
      hMC_eMax_SumMaxChNePh[b][i] = new TH1F(Form("hMC_eMax_SumMaxChNePh_%s_JetID_cent%d",isJetID[b],i),"",1000,0,1000);

    }

    hData_noCut[i] = new TH1F(Form("hData_noCut_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p9[i] = new TH1F(Form("hData_eMaxJtPt_0p9_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p8[i] = new TH1F(Form("hData_eMaxJtPt_0p8_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p7[i] = new TH1F(Form("hData_eMaxJtPt_0p7_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p6[i] = new TH1F(Form("hData_eMaxJtPt_0p6_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p5[i] = new TH1F(Form("hData_eMaxJtPt_0p5_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p4[i] = new TH1F(Form("hData_eMaxJtPt_0p4_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p3[i] = new TH1F(Form("hData_eMaxJtPt_0p3_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p2[i] = new TH1F(Form("hData_eMaxJtPt_0p2_cent%d",i),"",1000,0,1000);
    hData_eMaxJtPt_0p1[i] = new TH1F(Form("hData_eMaxJtPt_0p1_cent%d",i),"",1000,0,1000);

    hMC_noCut[i] = new TH1F(Form("hMC_noCut_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p9[i] = new TH1F(Form("hMC_eMaxJtPt_0p9_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p8[i] = new TH1F(Form("hMC_eMaxJtPt_0p8_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p7[i] = new TH1F(Form("hMC_eMaxJtPt_0p7_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p6[i] = new TH1F(Form("hMC_eMaxJtPt_0p6_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p5[i] = new TH1F(Form("hMC_eMaxJtPt_0p5_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p4[i] = new TH1F(Form("hMC_eMaxJtPt_0p4_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p3[i] = new TH1F(Form("hMC_eMaxJtPt_0p3_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p2[i] = new TH1F(Form("hMC_eMaxJtPt_0p2_cent%d",i),"",1000,0,1000);
    hMC_eMaxJtPt_0p1[i] = new TH1F(Form("hMC_eMaxJtPt_0p1_cent%d",i),"",1000,0,1000);

  }

  // Lets start to get the histograms from the Tree - Project is a complete waste of time since it has to run over the whole tree every time. 
  // best would be to set branch addresses: 
  // int L1_MB_1;
  // int L1_MB_p_1;
  Float_t L1_sj36_1;
  Float_t L1_sj52_1;
  // Float_t jetMB_1;
  Float_t jet55_1;
  Float_t jet65_1;
  Float_t jet80_1;
  // Float_t jetMB_p_1;
  Float_t L1_sj36_p_1;
  Float_t L1_sj52_p_1;
  Float_t jet55_p_1;
  Float_t jet65_p_1;
  Float_t jet80_p_1;
  float trgObjpt_1;
  float cent_1;
  float jtpt_1;
  float raw_1;
  // float eta_1;
  // float eta_1_CM;
  // float phi_1;
  float chMax_1;
  // float trkMax_1;
  float chSum_1;
  float phSum_1;
  float neSum_1;
  // float trkSum_1;
  float phMax_1;
  float neMax_1;
  float eMax_1;
  float muMax_1;
  float eSum_1;
  float muSum_1;
  float jtpu_1;

  jetData->SetBranchAddress("rawpt",&raw_1);
  jetData->SetBranchAddress("jtpt",&jtpt_1);
  jetData->SetBranchAddress("jtpu",&jtpu_1);
  jetData->SetBranchAddress("l1sj36",&L1_sj36_1);
  jetData->SetBranchAddress("l1sj36_prescl",&L1_sj36_p_1);
  jetData->SetBranchAddress("l1sj52",&L1_sj52_1);
  jetData->SetBranchAddress("l1sj52_prescl",&L1_sj52_p_1);
  jetData->SetBranchAddress("jet55",&jet55_1);
  jetData->SetBranchAddress("jet55_prescl",&jet55_p_1);
  jetData->SetBranchAddress("jet65",&jet65_1);
  jetData->SetBranchAddress("jet65_prescl",&jet65_p_1);
  jetData->SetBranchAddress("jet80",&jet80_1);
  jetData->SetBranchAddress("jet80_prescl",&jet80_p_1);
  jetData->SetBranchAddress("trgObjpt",&trgObjpt_1);
  jetData->SetBranchAddress("cent",&cent_1);
  jetData->SetBranchAddress("chMax",&chMax_1);
  jetData->SetBranchAddress("chSum",&chSum_1);
  jetData->SetBranchAddress("phMax",&phMax_1);
  jetData->SetBranchAddress("phSum",&phSum_1);
  jetData->SetBranchAddress("neMax",&neMax_1);
  jetData->SetBranchAddress("neSum",&neSum_1);
  jetData->SetBranchAddress("muMax",&muMax_1);
  jetData->SetBranchAddress("muSum",&muSum_1);
  jetData->SetBranchAddress("eMax",&eMax_1);
  jetData->SetBranchAddress("eSum",&eSum_1);
  Float_t effecPrescl = 2.047507;
  int centBin_1 = 0;

  // Set Variables for MC. (_2) 
  Float_t jet55_2;
  Float_t jet65_2;
  Float_t jet80_2;
  Float_t jet55_p_2;
  Float_t jet65_p_2;
  Float_t jet80_p_2;
  float cent_2;
  float jtpt_2;
  float raw_2;
  float refpt_2;
  float scale_2;
  float weight_vz_2;
  float weight_cent_2;
  float subid_2;
  float chMax_2;
  float chSum_2;
  float phSum_2;
  float neSum_2;
  float phMax_2;
  float neMax_2;
  float eMax_2;
  float muMax_2;
  float eSum_2;
  float muSum_2;
  float jtpu_2;

  jetMC->SetBranchAddress("rawpt",&raw_2);
  jetMC->SetBranchAddress("jtpt",&jtpt_2);
  jetMC->SetBranchAddress("jtpu",&jtpu_2);
  jetMC->SetBranchAddress("refpt",&refpt_2);
  jetMC->SetBranchAddress("scale",&scale_2);
  jetMC->SetBranchAddress("weight_cent",&weight_cent_2);
  jetMC->SetBranchAddress("weight_vz",&weight_vz_2);
  jetMC->SetBranchAddress("subid",&subid_2);
  jetMC->SetBranchAddress("jet55",&jet55_2);
  jetMC->SetBranchAddress("jet55_prescl",&jet55_p_2);
  jetMC->SetBranchAddress("jet65",&jet65_2);
  jetMC->SetBranchAddress("jet65_prescl",&jet65_p_2);
  jetMC->SetBranchAddress("jet80",&jet80_2);
  jetMC->SetBranchAddress("jet80_prescl",&jet80_p_2);
  jetMC->SetBranchAddress("cent",&cent_2);
  jetMC->SetBranchAddress("chMax",&chMax_2);
  jetMC->SetBranchAddress("chSum",&chSum_2);
  jetMC->SetBranchAddress("phMax",&phMax_2);
  jetMC->SetBranchAddress("phSum",&phSum_2);
  jetMC->SetBranchAddress("neMax",&neMax_2);
  jetMC->SetBranchAddress("neSum",&neSum_2);
  jetMC->SetBranchAddress("muMax",&muMax_2);
  jetMC->SetBranchAddress("muSum",&muSum_2);
  jetMC->SetBranchAddress("eMax",&eMax_2);
  jetMC->SetBranchAddress("eSum",&eSum_2);
  int centBin_2 = 0;

#if 0 
  // Data loop
  for(int jentry = 0;jentry<jetData->GetEntries();jentry++){

    jetData->GetEntry(jentry);
    centBin_1 = (int)cent_1;
    if(jentry%1000000==0)cout<<"Data "<<jentry<<" of "<<jetData->GetEntries()<<endl;

    hData_noCut[centBin_1]->Fill(jtpt_1);
    hData_noCut[nbins_cent]->Fill(jtpt_1);

    if((Float_t)eMax_1/jtpt_1<0.9) {
      hData_eMaxJtPt_0p9[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p9[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.8) {
      hData_eMaxJtPt_0p8[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p8[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.7) {
      hData_eMaxJtPt_0p7[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p7[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.6) {
      hData_eMaxJtPt_0p6[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p6[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.5) {
      hData_eMaxJtPt_0p5[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p5[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.4) {
      hData_eMaxJtPt_0p4[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p4[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.3) {
      hData_eMaxJtPt_0p3[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p3[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.2) {
      hData_eMaxJtPt_0p2[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p2[nbins_cent]->Fill(jtpt_1);
    }

    if((Float_t)eMax_1/jtpt_1<0.1) {
      hData_eMaxJtPt_0p1[centBin_1]->Fill(jtpt_1);
      hData_eMaxJtPt_0p1[nbins_cent]->Fill(jtpt_1);
    }

    hData_chMaxSumMaxChNePh_Pt[centBin_1]->Fill((Float_t)chMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_phMaxSumMaxChNePh_Pt[centBin_1]->Fill((Float_t)phMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_neMaxSumMaxChNePh_Pt[centBin_1]->Fill((Float_t)neMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_muMaxSumMaxChNePh_Pt[centBin_1]->Fill((Float_t)muMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_eMaxSumMaxChNePh_Pt[centBin_1]->Fill((Float_t)eMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);

    hData_chMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)chMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_phMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)phMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_neMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)neMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_muMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)muMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);
    hData_eMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)eMax_1/(chMax_1+neMax_1+phMax_1),jtpt_1);


    hData_chMax[centBin_1]->Fill(chMax_1);
    hData_chMax[nbins_cent]->Fill(chMax_1);
    hData_chMaxJtPt_Pt[centBin_1]->Fill((Float_t)chMax_1/jtpt_1,jtpt_1);
    hData_chMaxJtPt_Pt[nbins_cent]->Fill((Float_t)chMax_1/jtpt_1,jtpt_1);

    hData_chSum[centBin_1]->Fill(chSum_1);
    hData_chSum[nbins_cent]->Fill(chSum_1);
    hData_chSumJtPt_Pt[centBin_1]->Fill((Float_t)chSum_1/jtpt_1,jtpt_1);
    hData_chSumJtPt_Pt[nbins_cent]->Fill((Float_t)chSum_1/jtpt_1,jtpt_1);

    hData_phMax[centBin_1]->Fill(phMax_1);
    hData_phMax[nbins_cent]->Fill(phMax_1);
    hData_phMaxJtPt_Pt[centBin_1]->Fill((Float_t)phMax_1/jtpt_1,jtpt_1);
    hData_phMaxJtPt_Pt[nbins_cent]->Fill((Float_t)phMax_1/jtpt_1,jtpt_1);

    hData_phSum[centBin_1]->Fill(phSum_1);
    hData_phSum[nbins_cent]->Fill(phSum_1);
    hData_phSumJtPt_Pt[centBin_1]->Fill((Float_t)phSum_1/jtpt_1,jtpt_1);
    hData_phSumJtPt_Pt[nbins_cent]->Fill((Float_t)phSum_1/jtpt_1,jtpt_1);

    hData_neMax[centBin_1]->Fill(neMax_1);
    hData_neMax[nbins_cent]->Fill(neMax_1);
    hData_neMaxJtPt_Pt[centBin_1]->Fill((Float_t)neMax_1/jtpt_1,jtpt_1);
    hData_neMaxJtPt_Pt[nbins_cent]->Fill((Float_t)neMax_1/jtpt_1,jtpt_1);

    hData_neSum[centBin_1]->Fill(neSum_1);
    hData_neSum[nbins_cent]->Fill(neSum_1);
    hData_neSumJtPt_Pt[centBin_1]->Fill((Float_t)neSum_1/jtpt_1,jtpt_1);
    hData_neSumJtPt_Pt[nbins_cent]->Fill((Float_t)neSum_1/jtpt_1,jtpt_1);

    hData_muMax[centBin_1]->Fill(muMax_1);
    hData_muMax[nbins_cent]->Fill(muMax_1);
    hData_muMaxJtPt_Pt[centBin_1]->Fill((Float_t)muMax_1/jtpt_1,jtpt_1);
    hData_muMaxJtPt_Pt[nbins_cent]->Fill((Float_t)muMax_1/jtpt_1,jtpt_1);

    hData_muSum[centBin_1]->Fill(muSum_1);
    hData_muSum[nbins_cent]->Fill(muSum_1);
    hData_muSumJtPt_Pt[centBin_1]->Fill((Float_t)muSum_1/jtpt_1,jtpt_1);
    hData_muSumJtPt_Pt[nbins_cent]->Fill((Float_t)muSum_1/jtpt_1,jtpt_1);

    hData_eMax[centBin_1]->Fill(eMax_1);
    hData_eMax[nbins_cent]->Fill(eMax_1);
    hData_eMaxJtPt_Pt[centBin_1]->Fill((Float_t)eMax_1/jtpt_1,jtpt_1);
    hData_eMaxJtPt_Pt[nbins_cent]->Fill((Float_t)eMax_1/jtpt_1,jtpt_1);

    hData_eSum[centBin_1]->Fill(eSum_1);
    hData_eSum[nbins_cent]->Fill(eSum_1);
    hData_eSumJtPt_Pt[centBin_1]->Fill((Float_t)eSum_1/jtpt_1,jtpt_1);
    hData_eSumJtPt_Pt[nbins_cent]->Fill((Float_t)eSum_1/jtpt_1,jtpt_1);
 
    if(jet80_1 && L1_sj52_1){

      if((neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[2][0][centBin_1]->Fill(jtpt_1);
	hData[2][0][nbins_cent]->Fill(jtpt_1);
      }
      if((eMax_1/jtpt_1<0.7) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[2][1][centBin_1]->Fill(jtpt_1);
	hData[2][1][nbins_cent]->Fill(jtpt_1);
      }

    }else if(jet65_1 && L1_sj36_1 && !jet80_1){
  
      if((neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[1][0][centBin_1]->Fill(jtpt_1);
	hData[1][0][nbins_cent]->Fill(jtpt_1);
      }
      if((eMax_1/jtpt_1<0.7) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) &&  (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[1][1][centBin_1]->Fill(jtpt_1);
	hData[1][1][nbins_cent]->Fill(jtpt_1);
      }

    }else if(jet55_1 && L1_sj36_1 && !jet65_1 && !jet80_1){

      if((neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[0][0][centBin_1]->Fill(jtpt_1,effecPrescl);
	hData[0][0][nbins_cent]->Fill(jtpt_1,effecPrescl);
      }
      if((eMax_1/jtpt_1<0.7) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) &&  (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[0][1][centBin_1]->Fill(jtpt_1,effecPrescl);
	hData[0][1][nbins_cent]->Fill(jtpt_1,effecPrescl);
      }

    }

  
  // for(int i = 0;i<nbins_cent;i++){
    
  //   jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[0],isJetID[0],i),"jtpt",Jet55&&cent[i]);
  //   hData[0][0][i]->Print("base");
  //   jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[0],isJetID[1],i),"jtpt",Jet55&&elRjc&&neRjc&&phRjc&&chID&&muRjc &&cent[i]);
  //   hData[0][1][i]->Print("base");
    
  //   jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[1],isJetID[0],i),"jtpt",Jet65 && cent[i]);
  //   hData[1][0][i]->Print("base");
  //   jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[1],isJetID[1],i),"jtpt",Jet65 && elRjc && neRjc && phRjc && chID && muRjc && cent[i]);
  //   hData[1][1][i]->Print("base");
    
  //   jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[2],isJetID[0],i),"jtpt",Jet80 && cent[i]);  
  //   hData[2][0][i]->Print("base");
  //   jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[2],isJetID[1],i),"jtpt",Jet80 && elRjc && neRjc && phRjc && chID && muRjc && cent[i]);
  //   hData[2][1][i]->Print("base");
    
  // }

  // jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[0],isJetID[0],nbins_cent),"jtpt",Jet55);
  // hData[0][0][nbins_cent]->Print("base");
  // jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[0],isJetID[1],nbins_cent),"jtpt",Jet55 && elRjc && neRjc && phRjc && chID && muRjc);
  // hData[0][1][nbins_cent]->Print("base");

  // jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[1],isJetID[0],nbins_cent),"jtpt",Jet65);
  // hData[1][0][nbins_cent]->Print("base");
  // jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[1],isJetID[1],nbins_cent),"jtpt",Jet65 && elRjc && neRjc && phRjc && chID && muRjc);
  // hData[1][1][nbins_cent]->Print("base");

  // jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[2],isJetID[0],nbins_cent),"jtpt",Jet80);  
  // hData[2][0][nbins_cent]->Print("base");
  // jetData->Project(Form("hData_%s_%s_JetID_cent%d",TrigName[2],isJetID[1],nbins_cent),"jtpt",Jet80 && elRjc && neRjc && phRjc && chID && muRjc);
  // hData[2][1][nbins_cent]->Print("base");

  }// entry loop

#endif 
 
  Float_t weight = 0;
  
  // MC loop
  for(int jentry = 0;jentry<jetMC->GetEntries();jentry++){

    jetMC->GetEntry(jentry);
    centBin_2 = (int)cent_2;
    if(jentry%1000000==0)cout<<"MC "<<jentry<<" of "<<jetMC->GetEntries()<<endl;

    if(subid_2!=0)continue;

    weight = scale_2 * weight_vz_2 * weight_cent_2;

    hMC_noCut[centBin_2]->Fill(jtpt_2,weight);
    hMC_noCut[nbins_cent]->Fill(jtpt_2,weight);

    hMC_chMaxSumMaxChNePh_Pt[centBin_2]->Fill((Float_t)chMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_phMaxSumMaxChNePh_Pt[centBin_2]->Fill((Float_t)phMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_neMaxSumMaxChNePh_Pt[centBin_2]->Fill((Float_t)neMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_muMaxSumMaxChNePh_Pt[centBin_2]->Fill((Float_t)muMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_eMaxSumMaxChNePh_Pt[centBin_2]->Fill((Float_t)eMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);

    hMC_chMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)chMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_phMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)phMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_neMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)neMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_muMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)muMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);
    hMC_eMaxSumMaxChNePh_Pt[nbins_cent]->Fill((Float_t)eMax_2/(chMax_2+neMax_2+phMax_2),jtpt_2, weight);

    hMC_chMaxSumMaxChNePh_GenPt[centBin_2]->Fill((Float_t)chMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_phMaxSumMaxChNePh_GenPt[centBin_2]->Fill((Float_t)phMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_neMaxSumMaxChNePh_GenPt[centBin_2]->Fill((Float_t)neMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_muMaxSumMaxChNePh_GenPt[centBin_2]->Fill((Float_t)muMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_eMaxSumMaxChNePh_GenPt[centBin_2]->Fill((Float_t)eMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);

    hMC_chMaxSumMaxChNePh_GenPt[nbins_cent]->Fill((Float_t)chMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_phMaxSumMaxChNePh_GenPt[nbins_cent]->Fill((Float_t)phMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_neMaxSumMaxChNePh_GenPt[nbins_cent]->Fill((Float_t)neMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_muMaxSumMaxChNePh_GenPt[nbins_cent]->Fill((Float_t)muMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);
    hMC_eMaxSumMaxChNePh_GenPt[nbins_cent]->Fill((Float_t)eMax_2/(chMax_2+neMax_2+phMax_2),refpt_2, weight);

    if((Float_t)eMax_2/jtpt_2<0.9) {
      hMC_eMaxJtPt_0p9[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p9[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.8) {
      hMC_eMaxJtPt_0p8[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p8[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.7) {
      hMC_eMaxJtPt_0p7[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p7[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.6) {
      hMC_eMaxJtPt_0p6[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p6[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.5) {
      hMC_eMaxJtPt_0p5[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p5[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.4) {
      hMC_eMaxJtPt_0p4[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p4[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.3) {
      hMC_eMaxJtPt_0p3[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p3[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.2) {
      hMC_eMaxJtPt_0p2[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p2[nbins_cent]->Fill(jtpt_2,weight);
    }

    if((Float_t)eMax_2/jtpt_2<0.1) {
      hMC_eMaxJtPt_0p1[centBin_2]->Fill(jtpt_2,weight);
      hMC_eMaxJtPt_0p1[nbins_cent]->Fill(jtpt_2,weight);
    }



    hMC_chMax[centBin_2]->Fill(chMax_2, weight);
    hMC_chMax[nbins_cent]->Fill(chMax_2, weight);
    hMC_chMaxJtPt_Pt[centBin_2]->Fill((Float_t)chMax_2/jtpt_2,jtpt_2, weight);
    hMC_chMaxJtPt_Pt[nbins_cent]->Fill((Float_t)chMax_2/jtpt_2,jtpt_2, weight);

    hMC_chSum[centBin_2]->Fill(chSum_2, weight);
    hMC_chSum[nbins_cent]->Fill(chSum_2, weight);
    hMC_chSumJtPt_Pt[centBin_2]->Fill((Float_t)chSum_2/jtpt_2,jtpt_2, weight);
    hMC_chSumJtPt_Pt[nbins_cent]->Fill((Float_t)chSum_2/jtpt_2,jtpt_2, weight);

    hMC_phMax[centBin_2]->Fill(phMax_2, weight);
    hMC_phMax[nbins_cent]->Fill(phMax_2, weight);
    hMC_phMaxJtPt_Pt[centBin_2]->Fill((Float_t)phMax_2/jtpt_2,jtpt_2, weight);
    hMC_phMaxJtPt_Pt[nbins_cent]->Fill((Float_t)phMax_2/jtpt_2,jtpt_2, weight);

    hMC_phSum[centBin_2]->Fill(phSum_2, weight);
    hMC_phSum[nbins_cent]->Fill(phSum_2, weight);
    hMC_phSumJtPt_Pt[centBin_2]->Fill((Float_t)phSum_2/jtpt_2,jtpt_2, weight);
    hMC_phSumJtPt_Pt[nbins_cent]->Fill((Float_t)phSum_2/jtpt_2,jtpt_2, weight);

    hMC_neMax[centBin_2]->Fill(neMax_2, weight);
    hMC_neMax[nbins_cent]->Fill(neMax_2, weight);
    hMC_neMaxJtPt_Pt[centBin_2]->Fill((Float_t)neMax_2/jtpt_2,jtpt_2, weight);
    hMC_neMaxJtPt_Pt[nbins_cent]->Fill((Float_t)neMax_2/jtpt_2,jtpt_2, weight);

    hMC_neSum[centBin_2]->Fill(neSum_2, weight);
    hMC_neSum[nbins_cent]->Fill(neSum_2, weight);
    hMC_neSumJtPt_Pt[centBin_2]->Fill((Float_t)neSum_2/jtpt_2,jtpt_2, weight);
    hMC_neSumJtPt_Pt[nbins_cent]->Fill((Float_t)neSum_2/jtpt_2,jtpt_2, weight);

    hMC_muMax[centBin_2]->Fill(muMax_2, weight);
    hMC_muMax[nbins_cent]->Fill(muMax_2, weight);
    hMC_muMaxJtPt_Pt[centBin_2]->Fill((Float_t)muMax_2/jtpt_2,jtpt_2, weight);
    hMC_muMaxJtPt_Pt[nbins_cent]->Fill((Float_t)muMax_2/jtpt_2,jtpt_2, weight);

    hMC_muSum[centBin_2]->Fill(muSum_2, weight);
    hMC_muSum[nbins_cent]->Fill(muSum_2, weight);
    hMC_muSumJtPt_Pt[centBin_2]->Fill((Float_t)muSum_2/jtpt_2,jtpt_2, weight);
    hMC_muSumJtPt_Pt[nbins_cent]->Fill((Float_t)muSum_2/jtpt_2,jtpt_2, weight);

    hMC_eMax[centBin_2]->Fill(eMax_2, weight);
    hMC_eMax[nbins_cent]->Fill(eMax_2, weight);
    hMC_eMaxJtPt_Pt[centBin_2]->Fill((Float_t)eMax_2/jtpt_2,jtpt_2, weight);
    hMC_eMaxJtPt_Pt[nbins_cent]->Fill((Float_t)eMax_2/jtpt_2,jtpt_2, weight);

    hMC_eSum[centBin_2]->Fill(eSum_2, weight);
    hMC_eSum[nbins_cent]->Fill(eSum_2, weight);
    hMC_eSumJtPt_Pt[centBin_2]->Fill((Float_t)eSum_2/jtpt_2,jtpt_2, weight);
    hMC_eSumJtPt_Pt[nbins_cent]->Fill((Float_t)eSum_2/jtpt_2,jtpt_2, weight);

    // histograms with refpt instead of jtpt. 

    hMC_chMax[centBin_2]->Fill(chMax_2, weight);
    hMC_chMax[nbins_cent]->Fill(chMax_2, weight);
    hMC_chMaxGenPt_Pt[centBin_2]->Fill((Float_t)chMax_2/refpt_2,refpt_2, weight);
    hMC_chMaxGenPt_Pt[nbins_cent]->Fill((Float_t)chMax_2/refpt_2,refpt_2, weight);

    hMC_chSum[centBin_2]->Fill(chSum_2, weight);
    hMC_chSum[nbins_cent]->Fill(chSum_2, weight);
    hMC_chSumGenPt_Pt[centBin_2]->Fill((Float_t)chSum_2/refpt_2,refpt_2, weight);
    hMC_chSumGenPt_Pt[nbins_cent]->Fill((Float_t)chSum_2/refpt_2,refpt_2, weight);

    hMC_phMax[centBin_2]->Fill(phMax_2, weight);
    hMC_phMax[nbins_cent]->Fill(phMax_2, weight);
    hMC_phMaxGenPt_Pt[centBin_2]->Fill((Float_t)phMax_2/refpt_2,refpt_2, weight);
    hMC_phMaxGenPt_Pt[nbins_cent]->Fill((Float_t)phMax_2/refpt_2,refpt_2, weight);

    hMC_phSum[centBin_2]->Fill(phSum_2, weight);
    hMC_phSum[nbins_cent]->Fill(phSum_2, weight);
    hMC_phSumGenPt_Pt[centBin_2]->Fill((Float_t)phSum_2/refpt_2,refpt_2, weight);
    hMC_phSumGenPt_Pt[nbins_cent]->Fill((Float_t)phSum_2/refpt_2,refpt_2, weight);

    hMC_neMax[centBin_2]->Fill(neMax_2, weight);
    hMC_neMax[nbins_cent]->Fill(neMax_2, weight);
    hMC_neMaxGenPt_Pt[centBin_2]->Fill((Float_t)neMax_2/refpt_2,refpt_2, weight);
    hMC_neMaxGenPt_Pt[nbins_cent]->Fill((Float_t)neMax_2/refpt_2,refpt_2, weight);

    hMC_neSum[centBin_2]->Fill(neSum_2, weight);
    hMC_neSum[nbins_cent]->Fill(neSum_2, weight);
    hMC_neSumGenPt_Pt[centBin_2]->Fill((Float_t)neSum_2/refpt_2,refpt_2, weight);
    hMC_neSumGenPt_Pt[nbins_cent]->Fill((Float_t)neSum_2/refpt_2,refpt_2, weight);

    hMC_muMax[centBin_2]->Fill(muMax_2, weight);
    hMC_muMax[nbins_cent]->Fill(muMax_2, weight);
    hMC_muMaxGenPt_Pt[centBin_2]->Fill((Float_t)muMax_2/refpt_2,refpt_2, weight);
    hMC_muMaxGenPt_Pt[nbins_cent]->Fill((Float_t)muMax_2/refpt_2,refpt_2, weight);

    hMC_muSum[centBin_2]->Fill(muSum_2, weight);
    hMC_muSum[nbins_cent]->Fill(muSum_2, weight);
    hMC_muSumGenPt_Pt[centBin_2]->Fill((Float_t)muSum_2/refpt_2,refpt_2, weight);
    hMC_muSumGenPt_Pt[nbins_cent]->Fill((Float_t)muSum_2/refpt_2,refpt_2, weight);

    hMC_eMax[centBin_2]->Fill(eMax_2, weight);
    hMC_eMax[nbins_cent]->Fill(eMax_2, weight);
    hMC_eMaxGenPt_Pt[centBin_2]->Fill((Float_t)eMax_2/refpt_2,refpt_2, weight);
    hMC_eMaxGenPt_Pt[nbins_cent]->Fill((Float_t)eMax_2/refpt_2,refpt_2, weight);

    hMC_eSum[centBin_2]->Fill(eSum_2, weight);
    hMC_eSum[nbins_cent]->Fill(eSum_2, weight);
    hMC_eSumGenPt_Pt[centBin_2]->Fill((Float_t)eSum_2/refpt_2,refpt_2, weight);
    hMC_eSumGenPt_Pt[nbins_cent]->Fill((Float_t)eSum_2/refpt_2,refpt_2, weight);
    

    if(jet80_2){

      if((neMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (phMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/jtpt_2>0.05) && (muMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/(chMax_2+neMax_2+phMax_2)<0.9)){
	hMC[2][0][centBin_2]->Fill(jtpt_2,weight);
	hMC[2][0][nbins_cent]->Fill(jtpt_2,weight);
      }
      if((eMax_2/jtpt_2<0.7) && (neMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (phMax_2/(chMax_2+neMax_2+phMax_2)<0.9) &&(chMax_2/jtpt_2>0.05) && (muMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/(chMax_2+neMax_2+phMax_2)<0.9)){
	hMC[2][1][centBin_2]->Fill(jtpt_2,weight);
	hMC[2][1][nbins_cent]->Fill(jtpt_2,weight);
      }
      
    }else if(jet65_2 && !jet80_2){
      
      if((neMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (phMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/jtpt_2>0.05) && (muMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/(chMax_2+neMax_2+phMax_2)<0.9)){
        hMC[1][0][centBin_2]->Fill(jtpt_2,weight);
	hMC[1][0][nbins_cent]->Fill(jtpt_2,weight);
      }
      if((eMax_2/jtpt_2<0.7) && (neMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (phMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/jtpt_2>0.05) && (muMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/(chMax_2+neMax_2+phMax_2)<0.9)){
	hMC[1][1][centBin_2]->Fill(jtpt_2,weight);
	hMC[1][1][nbins_cent]->Fill(jtpt_2,weight);
      }

    }else if(jet55_2 && !jet65_2 && !jet80_2){

      if((neMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (phMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/jtpt_2>0.05) && (muMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/(chMax_2+neMax_2+phMax_2)<0.9)){
	hMC[0][0][centBin_2]->Fill(jtpt_2,weight);
	hMC[0][0][nbins_cent]->Fill(jtpt_2,weight);
      }
      if((eMax_2/jtpt_2<0.7) && (neMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (phMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/jtpt_2>0.05) && (muMax_2/(chMax_2+neMax_2+phMax_2)<0.9) && (chMax_2/(chMax_2+neMax_2+phMax_2)<0.9)){
	hMC[0][1][centBin_2]->Fill(jtpt_2,weight);
	hMC[0][1][nbins_cent]->Fill(jtpt_2,weight);
      }

    }
  }
  
  // Divide the histograms to get an idea of the effect to the Jet Spectra
  
  for(int i = 0;i<nbins_cent+1;i++){
    for(int a = 0;a<TrigValue;a++){

      hData_Ratio[a][i] = (TH1F*)hData[a][1][i]->Clone(Form("hData_Ratio_ID_over_noIDCut_%s_cent%d",TrigName[a],i));
      hData_Ratio[a][i]->Divide(hData[a][0][i]);
      hMC_Ratio[a][i] = (TH1F*)hMC[a][1][i]->Clone(Form("hMC_Ratio_ID_over_noIDCut_%s_cent%d",TrigName[a],i));
      hMC_Ratio[a][i]->Divide(hMC[a][0][i]);
    }
    hData_Ratio_eMaxJtPt_0p9[i] = (TH1F*)hData_eMaxJtPt_0p9[i]->Clone(Form("hData_Ratio_emaxJtPt_0p9_cent%d",i));
    hData_Ratio_eMaxJtPt_0p9[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p8[i] = (TH1F*)hData_eMaxJtPt_0p8[i]->Clone(Form("hData_Ratio_emaxJtPt_0p8_cent%d",i));
    hData_Ratio_eMaxJtPt_0p8[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p7[i] = (TH1F*)hData_eMaxJtPt_0p7[i]->Clone(Form("hData_Ratio_emaxJtPt_0p7_cent%d",i));
    hData_Ratio_eMaxJtPt_0p7[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p6[i] = (TH1F*)hData_eMaxJtPt_0p6[i]->Clone(Form("hData_Ratio_emaxJtPt_0p6_cent%d",i));
    hData_Ratio_eMaxJtPt_0p6[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p5[i] = (TH1F*)hData_eMaxJtPt_0p5[i]->Clone(Form("hData_Ratio_emaxJtPt_0p5_cent%d",i));
    hData_Ratio_eMaxJtPt_0p5[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p4[i] = (TH1F*)hData_eMaxJtPt_0p4[i]->Clone(Form("hData_Ratio_emaxJtPt_0p4_cent%d",i));
    hData_Ratio_eMaxJtPt_0p4[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p3[i] = (TH1F*)hData_eMaxJtPt_0p3[i]->Clone(Form("hData_Ratio_emaxJtPt_0p3_cent%d",i));
    hData_Ratio_eMaxJtPt_0p3[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p2[i] = (TH1F*)hData_eMaxJtPt_0p2[i]->Clone(Form("hData_Ratio_emaxJtPt_0p2_cent%d",i));
    hData_Ratio_eMaxJtPt_0p2[i]->Divide(hData_noCut[i]);
    hData_Ratio_eMaxJtPt_0p1[i] = (TH1F*)hData_eMaxJtPt_0p1[i]->Clone(Form("hData_Ratio_emaxJtPt_0p1_cent%d",i));
    hData_Ratio_eMaxJtPt_0p1[i]->Divide(hData_noCut[i]);

    hMC_Ratio_eMaxJtPt_0p9[i] = (TH1F*)hMC_eMaxJtPt_0p9[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p9_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p9[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p8[i] = (TH1F*)hMC_eMaxJtPt_0p8[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p8_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p8[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p7[i] = (TH1F*)hMC_eMaxJtPt_0p7[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p7_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p7[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p6[i] = (TH1F*)hMC_eMaxJtPt_0p6[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p6_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p6[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p5[i] = (TH1F*)hMC_eMaxJtPt_0p5[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p5_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p5[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p4[i] = (TH1F*)hMC_eMaxJtPt_0p4[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p4_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p4[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p3[i] = (TH1F*)hMC_eMaxJtPt_0p3[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p3_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p3[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p2[i] = (TH1F*)hMC_eMaxJtPt_0p2[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p2_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p2[i]->Divide(hMC_noCut[i]);
    hMC_Ratio_eMaxJtPt_0p1[i] = (TH1F*)hMC_eMaxJtPt_0p1[i]->Clone(Form("hMC_Ratio_emaxJtPt_0p1_cent%d",i));
    hMC_Ratio_eMaxJtPt_0p1[i]->Divide(hMC_noCut[i]);

  }

  //TFile f("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/RAA_JetID_MC_withCutHasAllExceptElecRejection.root","RECREATE");
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/RAA_JetID_mc_subid0_ElecCutRejectionStudy_with2DplotsSumMaxChNePhCut_%s%d%s_%d.root",algo,radius,jet_type,date.GetDate()),"RECREATE");
  f.cd();

  for(int i = 0;i<=nbins_cent;i++){
    for(int a = 0;a<TrigValue;a++){
      for(int b = 0;b<CutValue;b++){
	hData[a][b][i]->Write();
	hData[a][b][i]->Print("base");
	hMC[a][b][i]->Write();
	hMC[a][b][i]->Print("base");
      }
      hData_Ratio[a][i]->Write();
      hData_Ratio[a][i]->Print("base");
      hMC_Ratio[a][i]->Write();
      hMC_Ratio[a][i]->Print("base");
    } 

    hData_noCut[i]->Write();
    hMC_noCut[i]->Write();
    
    hData_chMax[i]->Write();
    hData_chSum[i]->Write();
    hData_neMax[i]->Write();
    hData_neSum[i]->Write();
    hData_phMax[i]->Write();
    hData_phSum[i]->Write();
    hData_muMax[i]->Write();
    hData_muSum[i]->Write();
    hData_eMax[i]->Write();
    hData_eSum[i]->Write();

    hMC_chMax[i]->Write();
    hMC_chSum[i]->Write();
    hMC_neMax[i]->Write();
    hMC_neSum[i]->Write();
    hMC_phMax[i]->Write();
    hMC_phSum[i]->Write();
    hMC_muMax[i]->Write();
    hMC_muSum[i]->Write();
    hMC_eMax[i]->Write();
    hMC_eSum[i]->Write();

    hData_chMaxJtPt_Pt[i]->Write();
    hData_chSumJtPt_Pt[i]->Write();
    hData_neMaxJtPt_Pt[i]->Write();
    hData_neSumJtPt_Pt[i]->Write();
    hData_phMaxJtPt_Pt[i]->Write();
    hData_phSumJtPt_Pt[i]->Write();
    hData_muMaxJtPt_Pt[i]->Write();
    hData_muSumJtPt_Pt[i]->Write();
    hData_eMaxJtPt_Pt[i]->Write();
    hData_eSumJtPt_Pt[i]->Write();

    hMC_chMaxJtPt_Pt[i]->Write();
    hMC_chSumJtPt_Pt[i]->Write();
    hMC_neMaxJtPt_Pt[i]->Write();
    hMC_neSumJtPt_Pt[i]->Write();
    hMC_phMaxJtPt_Pt[i]->Write();
    hMC_phSumJtPt_Pt[i]->Write();
    hMC_muMaxJtPt_Pt[i]->Write();
    hMC_muSumJtPt_Pt[i]->Write();
    hMC_eMaxJtPt_Pt[i]->Write();
    hMC_eSumJtPt_Pt[i]->Write();

    hMC_chMaxGenPt_Pt[i]->Write();
    hMC_chSumGenPt_Pt[i]->Write();
    hMC_neMaxGenPt_Pt[i]->Write();
    hMC_neSumGenPt_Pt[i]->Write();
    hMC_phMaxGenPt_Pt[i]->Write();
    hMC_phSumGenPt_Pt[i]->Write();
    hMC_muMaxGenPt_Pt[i]->Write();
    hMC_muSumGenPt_Pt[i]->Write();
    hMC_eMaxGenPt_Pt[i]->Write();
    hMC_eSumGenPt_Pt[i]->Write();    

    hData_chMaxSumMaxChNePh_Pt[i]->Write();
    hData_phMaxSumMaxChNePh_Pt[i]->Write();
    hData_neMaxSumMaxChNePh_Pt[i]->Write();
    hData_muMaxSumMaxChNePh_Pt[i]->Write();
    hData_eMaxSumMaxChNePh_Pt[i]->Write();

    hMC_chMaxSumMaxChNePh_Pt[i]->Write();
    hMC_phMaxSumMaxChNePh_Pt[i]->Write();
    hMC_neMaxSumMaxChNePh_Pt[i]->Write();
    hMC_muMaxSumMaxChNePh_Pt[i]->Write();
    hMC_eMaxSumMaxChNePh_Pt[i]->Write();

    hMC_chMaxSumMaxChNePh_GenPt[i]->Write();
    hMC_phMaxSumMaxChNePh_GenPt[i]->Write();
    hMC_neMaxSumMaxChNePh_GenPt[i]->Write();
    hMC_muMaxSumMaxChNePh_GenPt[i]->Write();
    hMC_eMaxSumMaxChNePh_GenPt[i]->Write();

    hData_eMaxJtPt_0p9[i]->Write();
    hData_eMaxJtPt_0p8[i]->Write();
    hData_eMaxJtPt_0p7[i]->Write();
    hData_eMaxJtPt_0p6[i]->Write();
    hData_eMaxJtPt_0p5[i]->Write();
    hData_eMaxJtPt_0p4[i]->Write();
    hData_eMaxJtPt_0p3[i]->Write();
    hData_eMaxJtPt_0p2[i]->Write();
    hData_eMaxJtPt_0p1[i]->Write();

    hData_Ratio_eMaxJtPt_0p9[i]->Write();
    hData_Ratio_eMaxJtPt_0p8[i]->Write();
    hData_Ratio_eMaxJtPt_0p7[i]->Write();
    hData_Ratio_eMaxJtPt_0p6[i]->Write();
    hData_Ratio_eMaxJtPt_0p5[i]->Write();
    hData_Ratio_eMaxJtPt_0p4[i]->Write();
    hData_Ratio_eMaxJtPt_0p3[i]->Write();
    hData_Ratio_eMaxJtPt_0p2[i]->Write();
    hData_Ratio_eMaxJtPt_0p1[i]->Write();

    hMC_eMaxJtPt_0p9[i]->Write();
    hMC_eMaxJtPt_0p8[i]->Write();
    hMC_eMaxJtPt_0p7[i]->Write();
    hMC_eMaxJtPt_0p6[i]->Write();
    hMC_eMaxJtPt_0p5[i]->Write();
    hMC_eMaxJtPt_0p4[i]->Write();
    hMC_eMaxJtPt_0p3[i]->Write();
    hMC_eMaxJtPt_0p2[i]->Write();
    hMC_eMaxJtPt_0p1[i]->Write();

    hMC_Ratio_eMaxJtPt_0p9[i]->Write();
    hMC_Ratio_eMaxJtPt_0p8[i]->Write();
    hMC_Ratio_eMaxJtPt_0p7[i]->Write();
    hMC_Ratio_eMaxJtPt_0p6[i]->Write();
    hMC_Ratio_eMaxJtPt_0p5[i]->Write();
    hMC_Ratio_eMaxJtPt_0p4[i]->Write();
    hMC_Ratio_eMaxJtPt_0p3[i]->Write();
    hMC_Ratio_eMaxJtPt_0p2[i]->Write();
    hMC_Ratio_eMaxJtPt_0p1[i]->Write();

  }

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}
