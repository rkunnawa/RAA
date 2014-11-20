// Raghav Kunnawalkam Elayavalli
// Nov 19th 2014
// Rutgers 
// raghav.k.e at CERN dot CH 

// Macro to make histograms for the HFVs validation study. 

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

/*
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
*/


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


static const int nbins_cent = 6;
static Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 5 to get your actual centrality
static Double_t ncoll[nbins_cent+1] = { 1660, 1310, 745, 251, 62.8, 10.8 ,362.24}; //last one is for 0-200 bin. 

//static const int no_radius = 3;//necessary for the RAA analysis  
//static const int list_radius[no_radius] = {2,3,4};

static const int no_radius = 1;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {3};

using namespace std;

void RAA_HFVsValidation(char *algo= "Vs", char *jet_type = "PF", int startfile = 0, int endfile = 1){

  TH1::SetDefaultSumw2();

  
  TStopwatch timer;
  timer.Start();

  TDatime date;
  
  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = true;
  
  // number convension:
  // 0 - MB
  // 1 - 55 or 65
  // 2 - 80 or 95 
  // 80 is the unprescaled trigger - yes
  // 
  
  // Now im going to change the file reading here for PbPb to look at the unmerged files through condor. 
  std::string infile1;
  //infile1 = "jet55or65_filelist.txt";
  infile1 = "jetRAA_PbPb_data_forest.txt";
  
  //std::string infile2;
  //infile2 = "jet80_filelist.txt";

  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;

  //std::ifstream instr2(infile2.c_str(),std::ifstream::in);
  //std::string filename2;

  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }

  //for(int ifile = 0;ifile<boundaries_fileno_job[startfile];ifile++){
  //  instr2>>filename2;
  //}
 
  const int N = 6;
    
  TChain *jetpbpb1[N][no_radius];
  //TChain *jetpbpb2[N][no_radius];
  
  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%d%sJetAnalyzer",algo,list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    dir[4][k] = "hltobject";
    dir[5][k] = "pfcandAnalyzer";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree",
    "jetObjTree",
    "pfTree"
  };
  
  //this loop is to assign the tree values before we go into the file loop. 
  for(int k = 0;k<no_radius;k++){
    for(int t = 0;t<N;t++){
      jetpbpb1[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data());
    }//tree loop ends
  }// radius loop ends
  
  for(int ifile = startfile;ifile<endfile;ifile++){
    
    instr1>>filename1;
    if(printDebug)cout<<"File: "<<filename1<<endl;
    for(int k = 0;k<no_radius;k++){

      for(int t = 0;t<N;t++){
	
	//jetpbpb1[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data()) ;
	jetpbpb1[t][k]->Add(filename1.c_str());
	if(printDebug)cout << "Tree loaded  " << string(dir[t][k]+"/"+trees[t]).data() << endl;
	if(printDebug)cout << "Entries : " << jetpbpb1[t][k]->GetEntries() << endl;
	
      }// tree loop ends
      
      //jetpbpb1->Add(filename.c_str());
      //hltpbpb1->Add(filename.c_str());
      //skmpbpb1->Add(filename.c_str());
      //hltobjpbpb1->Add(filename.c_str());
      //evtpbpb1->Add(filename.c_str());
      //cout<<"Entries = "<<jetpbpb1->GetEntries()<<endl;
    }// radius loop ends
    
  }// file loop ends
  

  for(int k = 0;k<no_radius;k++){
    jetpbpb1[2][k]->AddFriend(jetpbpb1[0][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[1][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[3][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[4][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[5][k]);
  }// radius loop ends
 
  
  //file 1: 
  // jet tree
  int nrefe_1;
  float pt_1[1000];
  //float old_pt3[1000];
  float raw_1[1000];
  float eta_1[1000];
  float eta_1_CM[1000];
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
  float jtpu_1[1000];

  // event tree
  int evt_1;
  int run_1;
  int lumi_1;
  int hiBin_1;
  float vx_1;
  float vy_1;
  float vz_1;
  int hiNpix_1;
  int hiNtracks_1;
  float hiHFminus_1;
  float hiHF_1;
  float hiHFplus_1;
  float hiHFplusEta4_1;
  float hiHFminusEta4_1;
  int pcollisionEventSelection_1;
  int pHBHENoiseFilter_1;
  int pprimaryvertexFilter_1;
  int pVertexFilterCutGplus_1;

  // trigger tree
  int L1_MB_1;
  int L1_MB_p_1;
  int L1_sj36_1;
  int L1_sj52_1;
  int jetMB_1;
  int jet55_1;
  int jet65_1;
  int jet80_1;
  int jetMB_p_1;
  int L1_sj36_p_1;
  int L1_sj52_p_1;
  int jet55_p_1;
  int jet65_p_1;
  int jet80_p_1;

  // trigger object tree - this contains the maximum value of the particular trigger object. 
  float trgObj_id_1;
  float trgObj_pt_1;
  float trgObj_eta_1;
  float trgObj_phi_1;
  float trgObj_mass_1;
  
  //from the pfcandanalyzer tree:
  Int_t nPFpart;
  Int_t pfID[10000];
  Float_t pfPt[10000];
  Float_t pfVsPt[10000];
  Float_t pfVsPtInitial[10000];
  Float_t pfArea[10000];
  Float_t pfEta[10000];
  Float_t pfPhi[10000];
  Float_t v_n[5][15];
  Float_t psi_n[5][15];
  Float_t sumpT[15];


  for(int k = 0;k<no_radius;k++){
    //set the branch addresses:  - one of the most boring parts of the code: 

    jetpbpb1[2][k]->SetBranchAddress("evt",&evt_1);
    jetpbpb1[2][k]->SetBranchAddress("run",&run_1);
    jetpbpb1[2][k]->SetBranchAddress("lumi",&lumi_1);
    jetpbpb1[2][k]->SetBranchAddress("hiBin",&hiBin_1);
    jetpbpb1[2][k]->SetBranchAddress("vz",&vz_1);
    jetpbpb1[2][k]->SetBranchAddress("vx",&vx_1);
    jetpbpb1[2][k]->SetBranchAddress("vy",&vy_1);
    jetpbpb1[2][k]->SetBranchAddress("hiNpix",&hiNpix_1);
    jetpbpb1[2][k]->SetBranchAddress("hiNtracks",&hiNtracks_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFminus",&hiHFminus_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHF",&hiHF_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFplus",&hiHFplus_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_1);
    jetpbpb1[2][k]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_1);
    jetpbpb1[2][k]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_1);
    //jetpbpb1[2][k]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_1);
    //jetpbpb1[2][k]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_1);
  
    jetpbpb1[2][k]->SetBranchAddress("nref",&nrefe_1);
    jetpbpb1[2][k]->SetBranchAddress("jtpt",&pt_1);
    jetpbpb1[2][k]->SetBranchAddress("jteta",&eta_1);
    jetpbpb1[2][k]->SetBranchAddress("jtphi",&phi_1);
    jetpbpb1[2][k]->SetBranchAddress("rawpt",&raw_1);
    jetpbpb1[2][k]->SetBranchAddress("jtpu",&jtpu_1);
    jetpbpb1[2][k]->SetBranchAddress("chargedMax",&chMax_1);
    jetpbpb1[2][k]->SetBranchAddress("chargedSum",&chSum_1);
    jetpbpb1[2][k]->SetBranchAddress("trackMax",&trkMax_1);
    jetpbpb1[2][k]->SetBranchAddress("trackSum",&trkSum_1);
    jetpbpb1[2][k]->SetBranchAddress("photonMax",&phMax_1);
    jetpbpb1[2][k]->SetBranchAddress("photonSum",&phSum_1);
    jetpbpb1[2][k]->SetBranchAddress("neutralMax",&neMax_1);
    jetpbpb1[2][k]->SetBranchAddress("neutralSum",&neSum_1);
    jetpbpb1[2][k]->SetBranchAddress("eSum",&eSum_1);
    jetpbpb1[2][k]->SetBranchAddress("eMax",&eMax_1);
    jetpbpb1[2][k]->SetBranchAddress("muSum",&muSum_1);
    jetpbpb1[2][k]->SetBranchAddress("muMax",&muMax_1);

    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v1",&jet55_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v1",&jet65_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v1",&jet80_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_1);

    jetpbpb1[2][k]->SetBranchAddress("id",&trgObj_id_1);
    jetpbpb1[2][k]->SetBranchAddress("pt",&trgObj_pt_1);
    jetpbpb1[2][k]->SetBranchAddress("eta",&trgObj_eta_1);
    jetpbpb1[2][k]->SetBranchAddress("phi",&trgObj_phi_1);
    jetpbpb1[2][k]->SetBranchAddress("mass",&trgObj_mass_1);

    jetpbpb1[2][k]->SetBranchAddress("nPFpart",&nPFpart);
    jetpbpb1[2][k]->SetBranchAddress("pfId",&pfID);
    jetpbpb1[2][k]->SetBranchAddress("pfPt",&pfPt);
    jetpbpb1[2][k]->SetBranchAddress("pfVsPt",&pfVsPt);
    jetpbpb1[2][k]->SetBranchAddress("pfVsPtInitial",&pfVsPtInitial);
    jetpbpb1[2][k]->SetBranchAddress("pfArea",&pfArea);
    jetpbpb1[2][k]->SetBranchAddress("pfEta",&pfEta);
    jetpbpb1[2][k]->SetBranchAddress("pfPhi",&pfPhi);
    jetpbpb1[2][k]->SetBranchAddress("vn",&v_n);
    jetpbpb1[2][k]->SetBranchAddress("psin",&psi_n);
    jetpbpb1[2][k]->SetBranchAddress("sumpt",&sumpT);


  }//radius loop

  // lets start making the 2D histograms necessary: 
  // 1) sumpT vs centrality
  // 2) Event Plane, official vs PF:

  // information about the objects in the pfcandanalyzer: 
  // The 15 eta regions are delimited by a detector-technology-motivated
  // split into the following edges, where the first ([0]) and last ([14])
  // are the HF- and HF+:
  // -5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522, 0.522,
  // 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191
  // vn[5][15] contains the flow magnitude from n = 0 to 4. 


  TH2F *hSumpTvsHF[15];
  TH2F *hSumpTvshiBin[15];
  
  TH2F *hPsivsHF[5][15];
  
  // declare the histograms in loops:
  for(int a = 0;a<15;a++){

    hSumpTvsHF[a] = new TH2F(Form("hSumpT_eta%d_vsHF",a),"",5000,0,10000,5000,0,10000);
    hSumpTvshiBin[a] = new TH2F(Form("hSumpT_eta%d_vshiBin",a),"",5000,0,10000,200,0,200);

    for(int b = 0;b<5;b++){

      hPsivsHF[b][a] = new TH2F(Form("hPsi%d_etabin%d_HF",b,a),"",630,-3.15,+3.15,5000,0,10000);

    }// no of flow components

  }// eta bin

  for(int k = 0;k<no_radius;k++){
    for(int ievt = 0;ievt<jetpbpb1[2][k]->GetEntries();ievt++){

      jetpbpb1[2][k]->GetEntry(ievt);
      if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 
      if(fabs(vz_1)>15) continue;

      for(int a = 0;a<15;a++){

	hSumpTvsHF[a]->Fill(sumpT[a],hiHF_1);
	hSumpTvshiBin[a]->Fill(sumpT[a],hiBin_1);

	for(int b = 0;b<5;b++){
	
	  hPsivsHF[b][a]->Fill(psi_n[b][a],hiHF_1);
	
	}// no of flow components 

      }// eta bin

    }// event loop

  }// radius loop
  
  // declare the output root file
  TFile fout(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_HFVsvalidation_plots_ak%s%s_R3_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  fout.cd();
  for(int a = 0;a<15;a++){
    hSumpTvsHF[a]->Write();
    if(printDebug)hSumpTvsHF[a]->Print("base");
    hSumpTvshiBin[a]->Write();
    if(printDebug)hSumpTvshiBin[a]->Print("base");
    for(int b = 0;b<5;b++){
      hPsivsHF[b][a]->Write();
      if(printDebug)hPsivsHF[b][a]->Print("base");
    }
  }
  fout.Close();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;


}// main function
