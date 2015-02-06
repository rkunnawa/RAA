// Raghav Kunnawalkam Elayavalli
// Rutgers 
// Jan 24th 2015
// for questions or comments: raghav.k.e at CERN dot CH 

//
// macro to load in the pp data ntuple and process the spectra required to do the analysis. 
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
static const char TrigName [TrigValue][256] = {"HLT40","HLT60","HLT80"};
static const char isJetID [CutValue][256] = {"without","with"};

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


void RAA_JetID_pp(int radius = 3, char *jet_type = "PF"){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  timer.Start();

  TDatime date;

  cout<<"Running for "<<jet_type<<" pp jets "<<endl;
  
  bool printDebug = true;

  TFile *fDatain = TFile::Open(Form("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/pp_data_ak%s_20150123.root",jet_type));
  TTree* jetData = (TTree*)fDatain->Get("jets_ID");

  if(printDebug)cout<<"jetData no of jets = "<<jetData->GetEntries()<<endl;

  TH1F *hData[3][2];
  for(int a = 0;a<TrigValue;a++){
    for(int b = 0;b<CutValue;b++){
      hData[a][b] = new TH1F(Form("hData_%s_%s_JetID",TrigName[a],isJetID[b]),Form("Data Jet Spectra from %s trigger %s JetID cuts",TrigName[a],isJetID[b]),1000,0,1000);
    }   
  }
  
  // set the branch addresses: 
  Float_t jet40_1;
  Float_t jet60_1;
  Float_t jet80_1;
  Float_t jet40_p_1;
  Float_t jet60_p_1;
  Float_t jet80_p_1;
  float trgObjpt_1;
  float raw_1;
  float chMax_1;
  float chSum_1;
  float phSum_1;
  float neSum_1;
  float phMax_1;
  float neMax_1;
  float eMax_1;
  float muMax_1;
  float eSum_1;
  float muSum_1;
  float jtpt_1;

  jetData->SetBranchAddress("rawpt",&raw_1);
  jetData->SetBranchAddress("jtpt",&jtpt_1);
  jetData->SetBranchAddress("jet40",&jet40_1);
  jetData->SetBranchAddress("jet40_prescl",&jet40_p_1);
  jetData->SetBranchAddress("jet60",&jet60_1);
  jetData->SetBranchAddress("jet60_prescl",&jet60_p_1);
  jetData->SetBranchAddress("jet80",&jet80_1);
  jetData->SetBranchAddress("jet80_prescl",&jet80_p_1);
  jetData->SetBranchAddress("trgObjpt",&trgObjpt_1);
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

  for(int jentry = 0; jentry<jetData->GetEntries(); jentry++){

    jetData->GetEntry(jentry);
    if(jentry%1000000==0)cout<<"Data "<<jentry<<" of "<<jetData->GetEntries()<<endl;

    if(raw_1 < 30) continue;
    
    if(jet80_1 && trgObjpt_1 >= 80){

      if(1>0){
	hData[2][0]->Fill(jtpt_1);
	hData[2][0]->Fill(jtpt_1);
      }
      if((eSum_1/(chSum_1+neSum_1+phSum_1+muSum_1)<0.7) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[2][1]->Fill(jtpt_1);
	hData[2][1]->Fill(jtpt_1);
      }

    }else if(jet60_1 && !jet80_1 && trgObjpt_1 >= 65 && trgObjpt_1 < 80){
  
      if(1>0/*(neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)*/){
	hData[1][0]->Fill(jtpt_1,jet60_p_1);
	hData[1][0]->Fill(jtpt_1,jet60_p_1);
      }
      if((eSum_1/(chSum_1+neSum_1+phSum_1+muSum_1)<0.7) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) &&  (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[1][1]->Fill(jtpt_1,jet60_p_1);
	hData[1][1]->Fill(jtpt_1,jet60_p_1);
      }

    }else if(jet40_1 && !jet60_1 && !jet80_1 && trgObjpt_1 >= 55 && trgObjpt_1<65){

      if(1>0/*(neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) && (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)*/){
	hData[0][0]->Fill(jtpt_1,jet40_p_1);
	hData[0][0]->Fill(jtpt_1,jet40_p_1);
      }
      if((eSum_1/(chSum_1+neSum_1+phSum_1+muSum_1)<0.7) && (neMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (phMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/jtpt_1>0.05) &&  (muMax_1/(chMax_1+neMax_1+phMax_1)<0.9) && (chMax_1/(chMax_1+neMax_1+phMax_1)<0.9)){
	hData[0][1]->Fill(jtpt_1,jet40_p_1);
	hData[0][1]->Fill(jtpt_1,jet40_p_1);
      }

    }

  }// entry loop

  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/RAA_JetID_pp_data_final_jetIDcut_ptGreater30_%d%s_%d.root",radius,jet_type,date.GetDate()),"RECREATE");
  f.cd();

  for(int a = 0;a<TrigValue;a++){
    for(int b = 0;b<CutValue;b++){
      hData[a][b]->Write();
      hData[a][b]->Print("base");
     }
  }
  

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;

}// macro end
