// Raghav Kunnawalkam Elayavalli
// April 19th 2015
// Rutgers 

//
// Macro to read in 2011 data and study how the trigger will be affected for the amount of statistics we need.
// we also need to understand how each trigger would need to be prescaled inorder to have the same error bars as we cross each trigger threshold 
//


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

static const int no_radius = 1;
static const int list_radius[no_radius] = {3};

using namespace std;

void RAA_2015Run2_JetTrigger(int startfile = 0, int endfile = 1, char *algo = "Pu", char *jet_type = "Calo"){

  TH1::SetDefaultSumw2();
  
  TStopwatch timer;
  timer.Start();
  
  TDatime date;
  
  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = false;
  
  std::string infile1;
  //infile1 = "jetRAA_MinBiasUPC_forest.txt";
  infile1 = "PbPb_HYDJETMB_5p02_forest.txt";
  //infile1 = "PbPb_HydjetMinBias_forest.txt";
  
  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;
  
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }  

  const int N = 5;
  
  TChain *jetpbpb1[N][no_radius];

  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%d%sJetAnalyzer",algo,list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    dir[4][k] = "hltobject";
    //dir[6][k] = "pfcandAnalyzer";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree",
    "jetObjTree",
    //"pfTree"
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
      
    }// radius loop ends
    
  }// file loop ends

  for(int k = 0;k<no_radius;k++){
    jetpbpb1[2][k]->AddFriend(jetpbpb1[0][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[1][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[3][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[4][k]);

  }// radius loop ends
 
  // Ok this should now work. 

  cout<<"total no of entries in the combined forest files = "<<jetpbpb1[2][0]->GetEntries()<<endl;
  
  //file 1: 
  // jet tree 1
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
  float hiHF_1;
  float hiHFminus_1;
  float hiHFplus_1;
  float hiHFplusEta4_1;
  float hiHFminusEta4_1;
  int pcollisionEventSelection_1;
  int pHBHENoiseFilter_1;
  int pprimaryvertexFilter_1;
  int pVertexFilterCutGplus_1;

  // trigger tree
  int L1_sj36_1;
  int L1_sj52_1;
  int jet55_1;
  int jet65_1;
  int jet80_1;
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
  
  /*
  Int_t nPFpart;
  Int_t pfId[NOBJECT_MAX];
  Float_t pfPt[NOBJECT_MAX];
  Float_t pfVsPtInitial[NOBJECT_MAX];
  Float_t pfVsPt[NOBJECT_MAX];
  Float_t pfEta[NOBJECT_MAX];
  Float_t pfPhi[NOBJECT_MAX];
  Float_t pfArea[NOBJECT_MAX];
  Float_t v_n[5][15];
  Float_t psi_n[5][15];
  Float_t sumpT[15];
  */
  
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

    /*
    jetpbpb1[2][k]->SetBranchAddress("nPFpart", &nPFpart);
    jetpbpb1[2][k]->SetBranchAddress("pfId", pfId);
    jetpbpb1[2][k]->SetBranchAddress("pfPt", pfPt);
    jetpbpb1[2][k]->SetBranchAddress("pfVsPtInitial", pfVsPtInitial);
    jetpbpb1[2][k]->SetBranchAddress("pfVsPt", pfVsPt);
    jetpbpb1[2][k]->SetBranchAddress("pfEta", pfEta);
    jetpbpb1[2][k]->SetBranchAddress("pfPhi", pfPhi);
    jetpbpb1[2][k]->SetBranchAddress("pfArea", pfArea);
    jetpbpb1[2][k]->SetBranchAddress("vn",&v_n);
    jetpbpb1[2][k]->SetBranchAddress("psin",&psi_n);
    jetpbpb1[2][k]->SetBranchAddress("sumpt",&sumpT);
    */
  }//radius loop

  // Declare the output file 
  TFile fout(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Inclusive_minbias_JetTrigger_2015Run2_5p02HYDJETMinBias_akPu3Calo_%d_%d.root",date.GetDate(), endfile),"RECREATE");
  fout.cd();
  
  // declare the histograms 
  TH1F * hJet80_HLT        = new TH1F("hJet80_HLT","Jets from HLT_HIJet80 trigger",250,0,500);
  TH1F * hJet65_HLT        = new TH1F("hJet65_HLT","Jets from HLT_HIJet65 trigger",250,0,500);
  TH1F * hJet55_HLT        = new TH1F("hJet55_HLT","Jets from HLT_HIJet55 trigger",250,0,500);
  TH1F * hJet80_HLT_prescl = new TH1F("hJet80_HLT_prescl","Jets from HLT_HIJet80 trigger with prescl in forest",250,0,500);
  TH1F * hJet65_HLT_prescl = new TH1F("hJet65_HLT_prescl","Jets from HLT_HIJet65 trigger with prescl in forest",250,0,500);
  TH1F * hJet55_HLT_prescl = new TH1F("hJet55_HLT_prescl","Jets from HLT_HIJet55 trigger with prescl in forest",250,0,500);
  TH1F * hJet80            = new TH1F("hJet80","Jets from events with leading jet > 80 GeV, prescl = 1",250,0,500);
  TH1F * hJet70            = new TH1F("hJet70","Jets from events with leading jet > 70 GeV, prescl = 1",250,0,500);
  TH1F * hJet60            = new TH1F("hJet60","Jets from events with leading jet > 60 GeV, prescl = 1",250,0,500);
  TH1F * hJet50            = new TH1F("hJet50","Jets from events with leading jet > 50 GeV, prescl = 1",250,0,500);
  TH1F * hJet40            = new TH1F("hJet40","Jets from events with leading jet > 40 GeV, prescl = 1",250,0,500);

  TH1F * hEvents           = new TH1F("hEvent","Total no of events",10,0,10);
  TH1F * hEvents_LeadJet   = new TH1F("hEvents_LeadJet","Leading Jet pT",250,0,500);
  TH1F * hEvents_LeadJet80 = new TH1F("hEvents_LeadJet80","leading Jet pt > 80 GeV",250,0,500);
  TH1F * hEvents_LeadJet70 = new TH1F("hEvents_LeadJet70","leading Jet pt > 70 GeV",250,0,500);
  TH1F * hEvents_LeadJet60 = new TH1F("hEvents_LeadJet60","leading Jet pt > 60 GeV",250,0,500);
  TH1F * hEvents_LeadJet50 = new TH1F("hEvents_LeadJet50","leading Jet pt > 50 GeV",250,0,500);
  TH1F * hEvents_LeadJet40 = new TH1F("hEvents_LeadJet40","leading Jet pt > 40 GeV",250,0,500);

  for(int k = 0; k<no_radius; k++){

    for(int nentry = 0; nentry < jetpbpb1[2][k]->GetEntries(); ++nentry){

      //if(nentry % 10000 == 0) cout<<nentry<<"/"<<jetpbpb1[2][k]->GetEntries()<<endl;

      jetpbpb1[0][k]->GetEntry(nentry);
      jetpbpb1[1][k]->GetEntry(nentry);
      jetpbpb1[2][k]->GetEntry(nentry);    
      jetpbpb1[3][k]->GetEntry(nentry);
      jetpbpb1[4][k]->GetEntry(nentry);

      hEvents->Fill(0);

      if(pHBHENoiseFilter_1 == 0 || pcollisionEventSelection_1 == 0 || TMath::Abs(vz_1) > 15) continue;

      hEvents->Fill(1);

      // get the from the trigger 
      for(int jentry = 0; jentry < nrefe_1; ++jentry){
	if(jet80_1){
	  hEvents->Fill(8);
	  hJet80_HLT->Fill(pt_1[jentry]);
	  hJet80_HLT_prescl->Fill(pt_1[jentry], jet80_p_1);
	}
	if(jet65_1){
	  hEvents->Fill(6);
	  hJet65_HLT->Fill(pt_1[jentry]);
	  hJet65_HLT_prescl->Fill(pt_1[jentry], jet65_p_1);
	}	
	if(jet55_1){
	  hEvents->Fill(5);
	  hJet55_HLT->Fill(pt_1[jentry]);
	  hJet55_HLT_prescl->Fill(pt_1[jentry], jet55_p_1);
	}
      }

      hEvents_LeadJet->Fill(pt_1[0]);
      
      if(pt_1[0] >= 80) {
	hEvents_LeadJet80->Fill(pt_1[0]);
	for(int jentry = 0; jentry < nrefe_1; ++jentry)
	  hJet80->Fill(pt_1[jentry]);
      }

      if(pt_1[0] >= 70) {
	hEvents_LeadJet70->Fill(pt_1[0]);
	for(int jentry = 0; jentry < nrefe_1; ++jentry)
	  hJet70->Fill(pt_1[jentry]);
      }

      if(pt_1[0] >= 60) {
	hEvents_LeadJet60->Fill(pt_1[0]);
	for(int jentry = 0; jentry < nrefe_1; ++jentry)
	  hJet60->Fill(pt_1[jentry]);
      }

      if(pt_1[0] >= 50) {
	hEvents_LeadJet50->Fill(pt_1[0]);
	for(int jentry = 0; jentry < nrefe_1; ++jentry)
	  hJet50->Fill(pt_1[jentry]);
      }
      
      if(pt_1[0] >= 40) {
	hEvents_LeadJet40->Fill(pt_1[0]);
	for(int jentry = 0; jentry < nrefe_1; ++jentry)
	  hJet40->Fill(pt_1[jentry]);
      }

    } // event loop

  } // radius loop

  fout.Write();
  fout.Close();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;



}







