// Raghav Kunnawalkam Elayavalli
// Nov 18th 2014
// Rutgers 
// raghav.k.e at CERN dot CH 

// Jet ID macro - To study the issue of fake jets in the Jet Reconstruction for the RAA analysis. First idea is that the electron reconstructed pT is off. Going to use analysis cuts : neMax/(neMax+chMax+phMax)<0.975 and muMax/(neMax+chMax+phMax)<0.975 and testing out this cut in Data and MC: eMax/(neMax+chMax+phMax)<0.975. 

// Dec 9th making the plots for PU since RAA is finally going to PU Jets 

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
  TFile *fDatain = TFile::Open(Form("/Users/keraghav/WORK/RAA/Output/PbPb_jetntuple_withEvtCuts_SuperNovaRejected_ak%s%d%s_20141209.root",algo,radius,jet_type));
  //TFile *fMCin = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_pp_mc_nocut_ak%s%s_20141118.root",algo,jet_type));

  TTree *jetData = (TTree*)fDatain->Get("jets_ID");
  //TTree *jetMC = (TTree*)fMCin->Get("jets_ID");

  //check if the trees are filled. 
  if(printDebug)cout<<"jetData entries = "<<jetData->GetEntries()<<endl;
  //if(printDebug)cout<<"jetMC entries   = "<<jetMC->GetEntries()<<endl;

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

  for(int i = 0;i<=nbins_cent;i++){

    for(int a = 0;a<TrigValue;a++){

      for(int b = 0;b<CutValue;b++){

        hData[a][b][i] = new TH1F(Form("hData_%s_%s_JetID_cent%d",TrigName[a],isJetID[b],i),Form("Jet Spectra from %s trigger %s JetID cuts in the centrality bin %s",TrigName[a],isJetID[b],centWidth[i]),1000,0,1000);

      }
      
    }
    
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
  float cent;
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
  jetData->SetBranchAddress("cent",&cent);
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

  int centBin = cent;

  for(int jentry = 0;jentry<jetData->GetEntries();jentry++){

    jetData->GetEntry(jentry);

    if(jentry%10000==0)cout<<jentry<<" of "<<jetData->GetEntries()<<endl;

    if(jet80_1 && L1_sj52_1){

      hData[2][0][centBin]->Fill(jtpt_1);
      hData[2][0][nbins_cent]->Fill(jtpt_1);
      
      if((eMax_1/jtpt_1<0.3) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[2][1][centBin]->Fill(jtpt_1);
	hData[2][1][nbins_cent]->Fill(jtpt_1);
      }

    }else if(jet65_1 && L1_sj36_1 && !jet80_1){

      hData[1][0][centBin]->Fill(jtpt_1);
      hData[1][0][nbins_cent]->Fill(jtpt_1);
      
      if((eMax_1/jtpt_1<0.3) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[1][1][centBin]->Fill(jtpt_1);
	hData[1][1][nbins_cent]->Fill(jtpt_1);
      }

    }else if(jet55_1 && L1_sj36_1 && !jet65_1 && !jet80_1){

      hData[0][0][centBin]->Fill(jtpt_1,effecPrescl);
      hData[0][0][nbins_cent]->Fill(jtpt_1,effecPrescl);
      
      if((eMax_1/jtpt_1<0.3) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[0][1][centBin]->Fill(jtpt_1,effecPrescl);
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

  // Divide the histograms to get an idea of the effect to the Jet Spectra
  
  for(int i = 0;i<nbins_cent+1;i++){
  
    for(int a = 0;a<TrigValue;a++){

      hData_Ratio[a][i] = (TH1F*)hData[a][1][i]->Clone(Form("hRatio_ID_over_noIDCut_%s_cent%d",TrigName[a],i));
      hData_Ratio[a][i]->Divide(hData[a][0][i]);
    }
    
  }

  TFile f("/Users/keraghav/WORK/RAA/Output/RAA_JetID_output.root","RECREATE");
  f.cd();

  for(int i = 0;i<=nbins_cent;i++){
    for(int a = 0;a<TrigValue;a++){
      for(int b = 0;b<CutValue;b++){
	hData[a][b][i]->Write();
	hData[a][b][i]->Print("base");
      }
      hData_Ratio[a][i]->Write();
      hData_Ratio[a][i]->Print("base");
    }
  }

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}
