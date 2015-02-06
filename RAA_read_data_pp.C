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

// Jan 23th 2015 - running it on the new files in MIT v29 JEC produced by pawan. 

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
#include <TNtuple.h>
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

// static const int nbins_eta = 14;
// static const double boundaries_eta[nbins_eta][2] = {
//   {-1.0,+1.0}, {-2.0,+2.0}, {-3.0,+3.0},
//   {-3.0,-2.5}, {-2.5,-2.0}, {-2.0,-1.5}, 
//   {-1.5,-1.0}, {-1.0,-0.5}, {-0.5,+0.5}, 
//   {+0.5,+1.0}, {+1.0,+1.5}, {+1.5,+2.0}, 
//   {+2.0,+2.5}, {+2.5,+3.0}
// };

// static const double delta_eta[nbins_eta] = {
//   2.0, 4.0, 6.0, 
//   0.5, 0.5, 0.5, 
//   0.5, 0.5, 1.0, 
//   0.5, 0.5, 0.5, 
//   0.5, 0.5
// };

// static const char etaWidth [nbins_eta][256] = {
//   "n10_eta_p10","n20_eta_p20","n30_eta_p30",
//   "n30_eta_n25","n25_eta_n20","n20_eta_n15",
//   "n15_eta_n10","n10_eta_n05","n05_eta_p05",
//   "p05_eta_p10","p10_eta_p15","p15_eta_p20",
//   "p20_eta_p25","p25_eta_p30"
// };

static const int nbins_eta = 1;
static const double boundaries_eta[nbins_eta][2] = {
  {-2.0,+2.0}
};

static const double delta_eta[nbins_eta] = {
  4.0
};

static const char etaWidth [nbins_eta][256] = {
  "n20_eta_p20"
};

static const int no_radius = 3;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {2,3,4};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 7; 
//static const int list_radius[no_radius] = {1,2,3,4,5,6,7};

using namespace std;

void RAA_read_data_pp(int startfile = 0,int endfile = 1,char *jet_type = "PF"){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TStopwatch timer;
  timer.Start();

  bool printDebug = false;
  
  if(printDebug) cout<<"Reading PP 2013 data"<<endl;

   // data files - pp 
  std::string infile1;
  infile1 = "jetRAA_pp_data_forest.txt";
  
  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;
  
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }

  const int N = 5;
  TChain * jetPP[N][no_radius];
  string dir[N][no_radius];

  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%d%sJetAnalyzer",list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    dir[4][k] = "hltobject";
  }
  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree",
    "jetObjTree",
  };

  //this loop is to assign the tree values before we go into the file loop. 
  for(int k = 0;k<no_radius;k++){
    for(int t = 0;t<N;t++){
      jetPP[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data());
    }//tree loop ends
  }// radius loop ends
  
  for(int ifile = startfile;ifile<endfile;ifile++){
    instr1>>filename1;
    if(printDebug)cout<<"File: "<<filename1<<endl;
    for(int k = 0;k<no_radius;k++){
      for(int t = 0;t<N;t++){	
	//jetpbpb1[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data()) ;
	jetPP[t][k]->Add(filename1.c_str());
	if(printDebug)cout << "Tree loaded  " << string(dir[t][k]+"/"+trees[t]).data() << endl;
	if(printDebug)cout << "Entries : " << jetPP[t][k]->GetEntries() << endl;
      }// tree loop ends
    }// radius loop ends    
  }// file loop ends

    for(int k = 0;k<no_radius;k++){
    jetPP[2][k]->AddFriend(jetPP[0][k]);
    jetPP[2][k]->AddFriend(jetPP[1][k]);
    jetPP[2][k]->AddFriend(jetPP[3][k]);
    jetPP[2][k]->AddFriend(jetPP[4][k]);
  }// radius loop ends
  
  if(printDebug) cout<<"no of entries in the forest = "<<jetPP[0][0]->GetEntries()<<endl;

  //declare the histograms here. currently they are made uisng the jet trigger and 12003 combination method. Once we have the new forests with the trigger objects we will change it to use the 14007 method. 

  TH1F *hpp_Trg80[no_radius][nbins_eta];
  TH1F *hpp_Trg60[no_radius][nbins_eta];
  TH1F *hpp_Trg40[no_radius][nbins_eta];
  TH1F *hpp_TrgComb[no_radius][nbins_eta];

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      hpp_Trg80[k][j] = new TH1F(Form("hpp_Trg80_R%d_%s",list_radius[k],etaWidth[j]),Form("Spectra from Jet Trigger 80 for R%d and %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Trg60[k][j] = new TH1F(Form("hpp_Trg60_R%d_%s",list_radius[k],etaWidth[j]),Form("Spectra from Jet Trigger 60 for R%d and %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Trg40[k][j] = new TH1F(Form("hpp_Trg40_R%d_%s",list_radius[k],etaWidth[j]),Form("Spectra from Jet Trigger 40 for R%d and %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_TrgComb[k][j] = new TH1F(Form("hpp_TrgComb_R%d_%s",list_radius[k],etaWidth[j]),Form("Combined Spectra from Jet Triggers for R%d and %s",list_radius[k],etaWidth[j]),1000,0,1000);
      
    }// eta bin loop
    
  }// radius loop


  //setting branch addresses: 

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
  float eMax_1[1000];
  float muMax_1[1000];
  float eSum_1[1000];
  float muSum_1[1000];

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
  float trgObj_id_1;
  float trgObj_pt_1;
  float trgObj_eta_1;
  float trgObj_phi_1;
  float trgObj_mass_1;

  
  
  for(int k = 0;k<no_radius;k++){
    //set the branch addresses:  - one of the most boring parts of the code: 
    jetPP[2][k]->SetBranchAddress("evt",&evt_1);
    jetPP[2][k]->SetBranchAddress("run",&run_1);
    jetPP[2][k]->SetBranchAddress("lumi",&lumi_1);
    jetPP[2][k]->SetBranchAddress("vz",&vz_1);
    jetPP[2][k]->SetBranchAddress("vx",&vx_1);
    jetPP[2][k]->SetBranchAddress("vy",&vy_1);
    // jetPP[2][k]->SetBranchAddress("hiNtracks",&hiNtracks_1);
    // jetPP[2][k]->SetBranchAddress("hiHFminus",&hiHFminus_1);
    // jetPP[2][k]->SetBranchAddress("hiHFplus",&hiHFplus_1);
    // jetPP[2][k]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_1);
    // jetPP[2][k]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_1);
    jetPP[2][k]->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA_1);
    jetPP[2][k]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_1);
    // jetPP[2][k]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_1);
    // jetPP[2][k]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_1);
  
    jetPP[2][k]->SetBranchAddress("nref",&nrefe_1);
    jetPP[2][k]->SetBranchAddress("jtpt",&pt_1);
    jetPP[2][k]->SetBranchAddress("jteta",&eta_1);
    jetPP[2][k]->SetBranchAddress("jtphi",&phi_1);
    jetPP[2][k]->SetBranchAddress("rawpt",&raw_1);
    jetPP[2][k]->SetBranchAddress("chargedMax",&chMax_1);
    jetPP[2][k]->SetBranchAddress("chargedSum",&chSum_1);
    jetPP[2][k]->SetBranchAddress("trackMax",&trkMax_1);
    jetPP[2][k]->SetBranchAddress("trackSum",&trkSum_1);
    jetPP[2][k]->SetBranchAddress("photonMax",&phMax_1);
    jetPP[2][k]->SetBranchAddress("photonSum",&phSum_1);
    jetPP[2][k]->SetBranchAddress("eMax",&eMax_1);
    jetPP[2][k]->SetBranchAddress("eSum",&eSum_1);
    jetPP[2][k]->SetBranchAddress("neutralMax",&neMax_1);
    jetPP[2][k]->SetBranchAddress("neutralSum",&neSum_1);
    jetPP[2][k]->SetBranchAddress("muMax",&muMax_1);
    jetPP[2][k]->SetBranchAddress("muSum",&muSum_1);

    //jetPP[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
    //jetPP[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
    //jetPP[2][k]->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
    //jetPP[2][k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
    jetPP[2][k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40_1);
    jetPP[2][k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p_1);
    jetPP[2][k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60_1);
    jetPP[2][k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p_1);
    jetPP[2][k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80_1);
    jetPP[2][k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p_1);

    jetPP[2][k]->SetBranchAddress("id",&trgObj_id_1);
    jetPP[2][k]->SetBranchAddress("pt",&trgObj_pt_1);
    jetPP[2][k]->SetBranchAddress("eta",&trgObj_eta_1);
    jetPP[2][k]->SetBranchAddress("phi",&trgObj_phi_1);
    jetPP[2][k]->SetBranchAddress("mass",&trgObj_mass_1);

  }//radius loop

  TH1F *hEvents_HLT80 = new TH1F("hEvents_HLT80","",4,0,2);
  TH1F *hEvents_HLT60 = new TH1F("hEvents_HLT60","",4,0,2);
  TH1F *hEvents_HLT40 = new TH1F("hEvents_HLT40","",4,0,2);

  //get all the pp spectra here: 
  //TCut pp3 = "jet40&&!jet60&&!jet80&&(chMax/pt)>0.01&&(TMath::Max(neMax,chMax)/TMath::Max(chSum,neSum))>=0.975";
  //TNtuple *jets_ID[no_radius];
  //for(int k = 0;k<no_radius;k++){
  //  jets_ID[k] = new TNtuple(Form("jets_R%d_ID",list_radius[k]),"","rawpt:jtpt:jet40:jet40_prescl:jet60:jet60_prescl:jet80:jet80_prescl:trgObjpt:chMax:chSum:phMax:phSum:neMax:neSum:muMax:muSum:eMax:eSum:trkMax:trkSum");
  //}  
  //Float_t arrayValues[21];
  
  for(int k = 0;k<no_radius;k++){

    if(printDebug) cout<<"Reading data for R = "<<list_radius[k]<<endl;
    Long64_t nentries_jetPP = jetPP[2][k]->GetEntries();
    if(printDebug) nentries_jetPP = 2;
    
    for(int jentry = 0;jentry<nentries_jetPP;jentry++){

      jetPP[0][k]->GetEntry(jentry);
      jetPP[1][k]->GetEntry(jentry);
      jetPP[2][k]->GetEntry(jentry);
      jetPP[3][k]->GetEntry(jentry);
      jetPP[4][k]->GetEntry(jentry);
      if(printDebug && jentry%100000==0)cout<<"pp data file"<<endl;
      if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_1<<"; run = "<<run_1<<endl;

      if(pHBHENoiseFilter_1 == 0 || pPAcollisionEventSelectionPA_1==0) continue;
      if(fabs(vz_1) > 15) continue;
      if(jet40_1 == 0 && jet60_1 == 0 && jet80_1 == 0) continue;

      if(jet80_1) hEvents_HLT80->Fill(1);
      if(jet60_1 && !jet80_1) hEvents_HLT60->Fill(1,jet60_p_1);
      if(jet60_1 && !jet80_1) hEvents_HLT60->Fill(0);
      if(jet40_1 && !jet60_1 && !jet80_1) hEvents_HLT40->Fill(1,jet40_p_1);
      if(jet40_1 && !jet60_1 && !jet80_1) hEvents_HLT40->Fill(0);

      for(int j = 0;j<nbins_eta;j++){

	for(int g = 0;g<nrefe_1;g++){

	  if(eta_1[g]<boundaries_eta[j][0] || eta_1[g]>=boundaries_eta[j][1]) continue;

#if 0
	  arrayValues[0] = raw_1[g];
	  arrayValues[1] = pt_1[g];
	  arrayValues[2] = jet40_1;
	  arrayValues[3] = jet40_p_1;
	  arrayValues[4] = jet60_1;
	  arrayValues[5] = jet60_p_1;
	  arrayValues[6] = jet80_1;
	  arrayValues[7] = jet80_p_1;
	  arrayValues[8] = trgObj_pt_1;
	  arrayValues[9] = chMax_1[g];
	  arrayValues[10] = chSum_1[g];
	  arrayValues[11] = phMax_1[g];
	  arrayValues[12] = phSum_1[g];
	  arrayValues[13] = neMax_1[g];
	  arrayValues[14] = neSum_1[g];
	  arrayValues[15] = muMax_1[g];
	  arrayValues[16] = muSum_1[g];
	  arrayValues[17] = eMax_1[g];
	  arrayValues[18] = eSum_1[g];
	  arrayValues[19] = trkMax_1[g];
	  arrayValues[20] = trkSum_1[g];

	  jets_ID[k]->Fill(arrayValues);
#endif

	  //if(raw_1[g] < 30) continue;
	  //if((neMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9) && (phMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9) && (chMax_1[g]/pt_1[g]>0.05) && (muMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9) && (chMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9)){
	  //if(eSum_1[g]/(chSum_1[g]+neSum_1[g]+phSum_1[g]+muSum_1[g])<0.7){
	  if(chMax_1[g]/pt_1[g]>0.05){
	  //if(1>0){
	  
	    if(jet80_1 && trgObj_pt_1>=80){
	      hpp_Trg80[k][j]->Fill(pt_1[g],jet80_p_1);
	    }
	    if(jet60_1==1 && trgObj_pt_1>=65 && trgObj_pt_1<80){
	      hpp_Trg60[k][j]->Fill(pt_1[g],jet60_p_1);
	    }
	    if(jet40_1==1 && trgObj_pt_1>=55 && trgObj_pt_1<65){
	      hpp_Trg40[k][j]->Fill(pt_1[g],jet40_p_1);
	    }
	    
	  }// qa condition

	}// jet loop
	
      }// eta bins
    }// event loop for jet80or100
  }// radius loop


  // output file declaration

  TDatime date;

  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/pp_data_spectra_trgObj_chMaxjtpt_norawptcut_ak%s_%d_%d.root",jet_type,date.GetDate(),endfile),"RECREATE");
  f.cd();

  hEvents_HLT80->Write();
  hEvents_HLT60->Write();
  hEvents_HLT40->Write();

  //for(int k = 0;k<no_radius;k++) jets_ID[k]->Write();

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      // normalize all the spectra to barns / delta pt delta eta

      // hpp_Trg80[k][j]->Scale(1./(5.3*1e12*delta_eta[j]));//5.3 pico barns^-1 is the lumi of the dataset. 
      // hpp_Trg60[k][j]->Scale(1./(5.3*1e12*delta_eta[j]));
      // hpp_Trg40[k][j]->Scale(1./(5.3*1e12*delta_eta[j]));
      
      // divideBinWidth(hpp_Trg80[k][j]);// divide by delta pt. 
      // divideBinWidth(hpp_Trg60[k][j]);
      // divideBinWidth(hpp_Trg40[k][j]);

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
