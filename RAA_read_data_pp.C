// Raghav Kunnawalam Elayavalli
// June 23rd 2014
// CERN

// 
// Macro separated from RAA_read_data to just read in the pp datasets - due to condor submission details. look at the read_data_pbpb macro to get more information. 
// 

// old prescale info: 
  // get the prescl factor information. 
  // Float_t presclpbpb3 = (Float_t)jetpbpb1_v2->GetEntries("jet80")/jetpbpb1_v2->GetEntries("jet55&&jet80");
  // cout<<"pbpb prescl3 = "<<presclpbpb3<<endl;//1.99871
  // Float_t presclpp3 = (Float_t)jetpp1_v2->GetEntries("jet80")/jetpp1_v2->GetEntries("jet40&&jet80");
  // cout<<"pp prescl3 = "<<presclpp3<<endl; //9.24968

// Aug 28th 2014 - read in data from the forest at CERN. will be changed to read in Aditya's ntuples once they are ready. 

#include <iostream>
#include <stdio.h>
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

// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}


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

static const int nbins_eta = 14;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0}, {-2.0,+2.0}, {-3.0,+3.0},
  {-3.0,-2.5}, {-2.5,-2.0}, {-2.0,-1.5}, 
  {-1.5,-1.0}, {-1.0,-0.5}, {-0.5,+0.5}, 
  {+0.5,+1.0}, {+1.0,+1.5}, {+1.5,+2.0}, 
  {+2.0,+2.5}, {+2.5,+3.0}
};

static const double delta_eta[nbins_eta] = {
  2.0, 4.0, 6.0, 
  0.5, 0.5, 0.5, 
  0.5, 0.5, 1.0, 
  0.5, 0.5, 0.5, 
  0.5, 0.5
};

static const char etaWidth [nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20","n30_eta_p30",
  "n30_eta_n25","n25_eta_n20","n20_eta_n15",
  "n15_eta_n10","n10_eta_n05","n05_eta_p05",
  "p05_eta_p10","p10_eta_p15","p15_eta_p20",
  "p20_eta_p25","p25_eta_p30"
};

static const int no_radius = 3;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {3,4,5};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 7; 
//static const int list_radius[no_radius] = {1,2,3,4,5,6,7};

using namespace std;

void RAA_read_data_pp(char *jet_type = "Calo"){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TStopwatch timer;
  timer.Start();

  bool printDebug = false;

  if(printDebug) cout<<"Reading PP 2013 data"<<endl;

   // data files - pp 
  TFile *fpp40or60 = TFile::Open("root://eoscms.cern.ch://store/group/phys_heavyions/yjlee/pp2013/promptReco/PP2013_HiForest_PromptReco_JSon_Jet40Jet60_ppTrack_forestv84.root");
  TFile *fpp80or100 = TFile::Open("root://eoscms.cern.ch://store/caf/user/yjlee/pp2013/promptReco/PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root");

  TTree *jetpp80or100[no_radius];
  TTree *jetpp40or60[no_radius];

  TTree *evtpp80or100 = (TTree*)fpp80or100->Get("hiEvtAnalyzer/HiTree");
  TTree *evtpp40or60 = (TTree*)fpp40or60->Get("hiEvtAnalyzer/HiTree");
  
  TTree *hltpp80or100 = (TTree*)fpp80or100->Get("hltanalysis/HltTree");
  TTree *hltpp40or60 = (TTree*)fpp40or60->Get("hltanalysis/HltTree");
 
  TTree *skmpp80or100 = (TTree*)fpp80or100->Get("skimanalysis/HltTree");
  TTree *skmpp40or60 = (TTree*)fpp40or60->Get("skimanalysis/HltTree");

  for(int k = 0;k<no_radius;k++){

    jetpp80or100[k] = (TTree*)fpp80or100->Get(Form("ak%d%sJetAnalyzer/t",list_radius[k],jet_type));
    jetpp40or60[k] = (TTree*)fpp40or60->Get(Form("ak%d%sJetAnalyzer/t",list_radius[k],jet_type));

    jetpp80or100[k]->AddFriend(evtpp80or100);
    jetpp40or60[k]->AddFriend(evtpp40or60);

    jetpp80or100[k]->AddFriend(hltpp80or100);
    jetpp40or60[k]->AddFriend(hltpp40or60);

    jetpp80or100[k]->AddFriend(skmpp80or100);
    jetpp40or60[k]->AddFriend(skmpp40or60);

  }//radius loop

  if(printDebug) cout<<"no of entries in the Jet80or100 tree = "<<jetpp80or100[0]->GetEntries()<<endl;
  if(printDebug) cout<<"no of entries in the Jet40or60 tree = "<<jetpp40or60[0]->GetEntries()<<endl;

  //declare the histograms here. currently they are made uisng the jet trigger and 12003 combination method. Once we have the new forests with the trigger objects we will change it to use the 14007 method. 

  TH1F *hpp_Trg80[no_radius][nbins_eta];
  TH1F *hpp_Trg60[no_radius][nbins_eta];
  TH1F *hpp_Trg40[no_radius][nbins_eta];
  TH1F *hpp_TrgComb[no_radius][nbins_eta];

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      hpp_Trg80[k][j] = new TH1F(Form("hpp_Trg80_R%d_%s",list_radius[k],etaWidth[j]),Form("Spectra from Jet Trigger 80 for R%d and %s",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      hpp_Trg60[k][j] = new TH1F(Form("hpp_Trg60_R%d_%s",list_radius[k],etaWidth[j]),Form("Spectra from Jet Trigger 60 for R%d and %s",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      hpp_Trg40[k][j] = new TH1F(Form("hpp_Trg40_R%d_%s",list_radius[k],etaWidth[j]),Form("Spectra from Jet Trigger 40 for R%d and %s",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      hpp_TrgComb[k][j] = new TH1F(Form("hpp_TrgComb_R%d_%s",list_radius[k],etaWidth[j]),Form("Combined Spectra from Jet Triggers for R%d and %s",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      
    }// eta bin loop
    
  }// radius loop


  //setting branch addresses: 

  //file 1: Jet80or100
  // jet tree
  int nrefe_1;
  float pt_1[1000];
  // float old_pt3[1000];
  float raw_1[1000];
  float eta_1[1000];
  // float eta_1_CM[1000];
  float phi_1[1000];
  float chMax_1[1000];
  float trkMax_1[1000];
  float chSum_1[1000];
  float phSum_1[1000];
  float neSum_1[1000];
  float trkSum_1[1000];
  float phMax_1[1000];
  float neMax_1[1000];

  // event tree
  int evt_1;
  int run_1;
  int lumi_1;
  float vx_1;
  float vy_1;
  float vz_1;
  // int hiNtracks_1;
  // float hiHFminus_1;
  // float hiHFplus_1;
  // float hiHFplusEta4_1;
  // float hiHFminusEta4_1;
  int pPAcollisionEventSelectionPA_1;
  int pHBHENoiseFilter_1;
  // int pprimaryvertexFilter_1;
  // int pVertexFilterCutGplus_1;

  // trigger tree
  // int L1_MB_1;
  // int L1_MB_p_1;
  // int jetMB_1;
  int jet40_1;
  int jet60_1;
  int jet80_1;
  // int jetMB_p_1;
  int jet40_p_1;
  int jet60_p_1;
  int jet80_p_1;


  // trigger object tree - this contains the maximum value of the particular trigger object. 
  // float trgObj_id_1;
  // float trgObj_pt_1;
  // float trgObj_eta_1;
  // float trgObj_phi_1;
  // float trgObj_mass_1;
  
  //file 2: Jet40or60
  // jet tree
  int nrefe_2;
  float pt_2[1000];
  //float old_pt3[1000];
  float raw_2[1000];
  float eta_2[1000];
  // float eta_2_CM[1000];
  float phi_2[1000];
  float chMax_2[1000];
  float trkMax_2[1000];
  float chSum_2[1000];
  float phSum_2[1000];
  float neSum_2[1000];
  float trkSum_2[1000];
  float phMax_2[1000];
  float neMax_2[1000];

  // event tree
  int evt_2;
  int run_2;
  int lumi_2;
  float vx_2;
  float vy_2;
  float vz_2;
  // int hiNtracks_2;
  // float hiHFminus_2;
  // float hiHFplus_2;
  // float hiHFplusEta4_2;
  // float hiHFminusEta4_2;
  int pPAcollisionEventSelectionPA_2;
  int pHBHENoiseFilter_2;
  // int pprimaryvertexFilter_2;
  // int pVertexFilterCutGplus_2;

  // trigger tree
  // int L1_MB_2;
  // int L1_MB_p_2;
  // int jetMB_2;
  int jet40_2;
  int jet60_2;
  int jet80_2;
  // int jetMB_p_2;
  int jet40_p_2;
  int jet60_p_2;
  int jet80_p_2;


  // trigger object tree - this contains the maximum value of the particular trigger object. 
  // float trgObj_id_2;
  // float trgObj_pt_2;
  // float trgObj_eta_2;
  // float trgObj_phi_2;
  // float trgObj_mass_2;

  
  
  for(int k = 0;k<no_radius;k++){
    //set the branch addresses:  - one of the most boring parts of the code: 
    jetpp80or100[k]->SetBranchAddress("evt",&evt_1);
    jetpp80or100[k]->SetBranchAddress("run",&run_1);
    jetpp80or100[k]->SetBranchAddress("lumi",&lumi_1);
    jetpp80or100[k]->SetBranchAddress("vz",&vz_1);
    jetpp80or100[k]->SetBranchAddress("vx",&vx_1);
    jetpp80or100[k]->SetBranchAddress("vy",&vy_1);
    // jetpp80or100[k]->SetBranchAddress("hiNtracks",&hiNtracks_1);
    // jetpp80or100[k]->SetBranchAddress("hiHFminus",&hiHFminus_1);
    // jetpp80or100[k]->SetBranchAddress("hiHFplus",&hiHFplus_1);
    // jetpp80or100[k]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_1);
    // jetpp80or100[k]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_1);
    jetpp80or100[k]->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA_1);
    jetpp80or100[k]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_1);
    // jetpp80or100[k]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_1);
    // jetpp80or100[k]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_1);
  
    jetpp80or100[k]->SetBranchAddress("nref",&nrefe_1);
    jetpp80or100[k]->SetBranchAddress("jtpt",&pt_1);
    jetpp80or100[k]->SetBranchAddress("jteta",&eta_1);
    jetpp80or100[k]->SetBranchAddress("jtphi",&phi_1);
    jetpp80or100[k]->SetBranchAddress("rawpt",&raw_1);
    jetpp80or100[k]->SetBranchAddress("chargedMax",&chMax_1);
    jetpp80or100[k]->SetBranchAddress("chargedSum",&chSum_1);
    jetpp80or100[k]->SetBranchAddress("trackMax",&trkMax_1);
    jetpp80or100[k]->SetBranchAddress("trackSum",&trkSum_1);
    jetpp80or100[k]->SetBranchAddress("photonMax",&phMax_1);
    jetpp80or100[k]->SetBranchAddress("photonSum",&phSum_1);
    jetpp80or100[k]->SetBranchAddress("neutralMax",&neMax_1);
    jetpp80or100[k]->SetBranchAddress("neutralSum",&neSum_1);

    //jetpp80or100[k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
    //jetpp80or100[k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
    //jetpp80or100[k]->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
    //jetpp80or100[k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
    jetpp80or100[k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40_1);
    jetpp80or100[k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p_1);
    jetpp80or100[k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60_1);
    jetpp80or100[k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p_1);
    jetpp80or100[k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80_1);
    jetpp80or100[k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p_1);

    // jetpp80or100[k]->SetBranchAddress("id",&trgObj_id_1);
    // jetpp80or100[k]->SetBranchAddress("pt",&trgObj_pt_1);
    // jetpp80or100[k]->SetBranchAddress("eta",&trgObj_eta_1);
    // jetpp80or100[k]->SetBranchAddress("phi",&trgObj_phi_1);
    // jetpp80or100[k]->SetBranchAddress("mass",&trgObj_mass_1);

    //set the branch addresses:  - one of the most boring parts of the code: 
    jetpp40or60[k]->SetBranchAddress("evt",&evt_2);
    jetpp40or60[k]->SetBranchAddress("run",&run_2);
    jetpp40or60[k]->SetBranchAddress("lumi",&lumi_2);
    jetpp40or60[k]->SetBranchAddress("vz",&vz_2);
    jetpp40or60[k]->SetBranchAddress("vx",&vx_2);
    jetpp40or60[k]->SetBranchAddress("vy",&vy_2);
    // jetpp40or60[k]->SetBranchAddress("hiNtracks",&hiNtracks_2);
    // jetpp40or60[k]->SetBranchAddress("hiHFminus",&hiHFminus_2);
    // jetpp40or60[k]->SetBranchAddress("hiHFplus",&hiHFplus_2);
    // jetpp40or60[k]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_2);
    // jetpp40or60[k]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_2);
    jetpp40or60[k]->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA_2);
    jetpp40or60[k]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_2);
    // jetpp40or60[k]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_2);
    // jetpp40or60[k]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_2);
  
    jetpp40or60[k]->SetBranchAddress("nref",&nrefe_2);
    jetpp40or60[k]->SetBranchAddress("jtpt",&pt_2);
    jetpp40or60[k]->SetBranchAddress("jteta",&eta_2);
    jetpp40or60[k]->SetBranchAddress("jtphi",&phi_2);
    jetpp40or60[k]->SetBranchAddress("rawpt",&raw_2);
    jetpp40or60[k]->SetBranchAddress("chargedMax",&chMax_2);
    jetpp40or60[k]->SetBranchAddress("chargedSum",&chSum_2);
    jetpp40or60[k]->SetBranchAddress("trackMax",&trkMax_2);
    jetpp40or60[k]->SetBranchAddress("trackSum",&trkSum_2);
    jetpp40or60[k]->SetBranchAddress("photonMax",&phMax_2);
    jetpp40or60[k]->SetBranchAddress("photonSum",&phSum_2);
    jetpp40or60[k]->SetBranchAddress("neutralMax",&neMax_2);
    jetpp40or60[k]->SetBranchAddress("neutralSum",&neSum_2);

    // jetpp40or60[k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_2);
    // jetpp40or60[k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_2);
    // jetpp40or60[k]->SetBranchAddress("L1_ZeroBias",&L1_MB_2);
    // jetpp40or60[k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_2);
    jetpp40or60[k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40_2);
    jetpp40or60[k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p_2);
    jetpp40or60[k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60_2);
    jetpp40or60[k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p_2);
    jetpp40or60[k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80_2);
    jetpp40or60[k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p_2);

    // jetpp40or60[k]->SetBranchAddress("id",&trgObj_id_2);
    // jetpp40or60[k]->SetBranchAddress("pt",&trgObj_pt_2);
    // jetpp40or60[k]->SetBranchAddress("eta",&trgObj_eta_2);
    // jetpp40or60[k]->SetBranchAddress("phi",&trgObj_phi_2);
    // jetpp40or60[k]->SetBranchAddress("mass",&trgObj_mass_2);
  
  }//radius loop

  //get all the pp spectra here: 
  TCut pp3 = "jet40&&!jet60&&!jet80&&(chMax/pt)>0.01&&(TMath::Max(neMax,chMax)/TMath::Max(chSum,neSum))>=0.975";
  
  for(int k = 0;k<no_radius;k++){

    if(printDebug) cout<<"Reading data for R = "<<list_radius[k]<<endl;
    Long64_t nentries_jet80or100 = jetpp80or100[k]->GetEntries();
    if(printDebug) nentries_jet80or100 = 2;
    
    for(int jentry = 0;jentry<nentries_jet80or100;jentry++){

      jetpp80or100[k]->GetEntry(jentry);
      if(printDebug && jentry%100000==0)cout<<"Jet 80or100 file"<<endl;
      if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_1<<"; run = "<<run_1<<endl;

      for(int j = 0;j<nbins_eta;j++){
	
	if(pHBHENoiseFilter_1 && pPAcollisionEventSelectionPA_1 &&(vz_1<15)){

	  for(int g = 0;g<nrefe_1;g++){
	    
	    if(/* put your favourite QA cut here */1>0){
	      
	      if(eta_1[g]>=boundaries_eta[j][0] && eta_1[g]<boundaries_eta[j][1]){
		
		if(jet80_1) hpp_Trg80[k][j]->Fill(pt_1[g],jet80_p_1);
		
	      }// eta selection
	      
	    }// qa cut selection
	    
	  }// jet loop
	  
	}// event selection cuts

      }// eta bin loop

    }// event loop for jet80or100

    Long64_t nentries_jet40or60 = jetpp40or60[k]->GetEntries();
    if(printDebug) nentries_jet40or60 = 2;
    
    for(int jentry = 0;jentry<nentries_jet40or60;jentry++){
      
      jetpp40or60[k]->GetEntry(jentry);
      if(printDebug && jentry%100000==0)cout<<"Jet 40or60 file"<<endl;
      if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_2<<"; run = "<<run_2<<endl;

      for(int j = 0;j<nbins_eta;j++){
	
	if(pHBHENoiseFilter_2 && pPAcollisionEventSelectionPA_2 &&(vz_2<15)){

	  for(int g = 0;g<nrefe_2;g++){
	    
	    if(/* put your favourite QA cut here */1>0){
	      
	      if(eta_2[g]>=boundaries_eta[j][0] && eta_2[g]<boundaries_eta[j][1]){
		
		if(jet60_2) hpp_Trg60[k][j]->Fill(pt_2[g],jet60_p_2);

		if(jet40_2) hpp_Trg40[k][j]->Fill(pt_2[g],jet40_p_2);
		
	      }// eta selection
	      
	    }// qa cut selection
	    
	  }// jet loop
	  
	}// event selection cuts

      }// eta bin loop      

    }// event loop for jet40or60

  }// radius loop


  // output file declaration

  TDatime date;

  TFile f(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_data_ak%s_%d.root",jet_type,date.GetDate()),"RECREATE");
  f.cd();

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      hpp_Trg80[k][j]->Scale(1./5300e6/delta_eta[j]);//5300e6 is the lumi of the dataset. 
      hpp_Trg60[k][j]->Scale(1./5300e6/delta_eta[j]);
      hpp_Trg40[k][j]->Scale(1./5300e6/delta_eta[j]);
      
      divideBinWidth(hpp_Trg80[k][j]);// divide by delta pt. 
      divideBinWidth(hpp_Trg60[k][j]);
      divideBinWidth(hpp_Trg40[k][j]);

      hpp_TrgComb[k][j]->Add(hpp_Trg80[k][j]);
      hpp_TrgComb[k][j]->Add(hpp_Trg60[k][j]);
      hpp_TrgComb[k][j]->Add(hpp_Trg40[k][j]);

      hpp_Trg80[k][j]->Write();
      if(printDebug)hpp_Trg80[k][j]->Print();
      hpp_Trg60[k][j]->Write();
      if(printDebug)hpp_Trg60[k][j]->Print();
      hpp_Trg40[k][j]->Write();
      if(printDebug)hpp_Trg40[k][j]->Print();
      hpp_TrgComb[k][j]->Write();
      if(printDebug)hpp_TrgComb[k][j]->Print();

    }// eta bin loop

  }// radius loop
  
  f.Write();

  f.Close();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;
  

}
