// Raghav Kunnawalkam Elayavalli
// Dec 19th 2014
// Rutgers 
// comments and questions: raghav.k.e at CERN dot CH

//
// macro to read in PbPb minbias data and make ntuples to study the pf electron problems 
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

//static const int nbins_pt = 29;
//static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

//only for the jet80 merged file - not necessary now 
//static const int job_no = 11;
//static const double boundaries_loopno_job[job_no+1] = {0,100000,200000,300000,400000,50000,600000,700000,800000,900000,1000000,1152308};

//static const double boundaries_fileno_job[job_no+1] = {0, 413, 826, 1239, 1652, 2065, 2478, 2891, 3304, 3717, 4130, 4542};

#define NOBJECT_MAX 16384

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

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 7; 
//static const int list_radius[no_radius] = {1,2,3,4,5,6,7};


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

int findBin(int hiBin){
  int binNo = 0;

  for(int i = 0;i<nbins_cent;i++){
    if(hiBin>=5*boundaries_cent[i] && hiBin<5*boundaries_cent[i+1]) {
      binNo = i;
      break;
    }
  }

  return binNo;
}


using namespace std;

void RAA_read_MinBias(int startfile = 0, int endfile = 1, char *algo = "Pu", char *jet_type = "PF"){

  TH1::SetDefaultSumw2();
  
  TStopwatch timer;
  timer.Start();
  
  TDatime date;
  
  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = false;

  // Now im going to change the file reading here for PbPb to look at the unmerged files through condor. 
  std::string infile1;
  //infile1 = "jet55or65_filelist.txt";
  infile1 = "PbPb_HydjetMinBias_forest.txt";
  
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
    //dir[3][k] = Form("akPu%d%sJetAnalyzer",list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    dir[4][k] = "hltobject";
    //dir[6][k] = "pfcandAnalyzer";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    //"t",
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


  // we need to add the histograms to find the jet spectra from normal and failure mode- infact just add them to the ntuples per event the value of the HFSumpT*vn*cos/sin(n*psi_n) so we can plot the spectra at the final stage. this would make things easier. 
  TNtuple *jets_ID = new TNtuple("jets_ID","","rawpt:jtpt:jtpu:l1sj36:l1sj36_prescl:l1sj52:l1sj52_prescl:jet55:jet55_prescl:jet65:jet65_prescl:jet80:jet80_prescl:trgObjpt:cent:chMax:chSum:phMax:phSum:neMax:neSum:muMax:muSum:eMax:eSum");
  Float_t arrayValues[25];
 
  for(int k = 0;k<no_radius;k++){

    if(printDebug)cout<<"Running data reading for R = "<<list_radius[k]<<endl;
    // loop for the jetpbpb1[2] tree 
    Long64_t nentries_jet55or65 = jetpbpb1[2][k]->GetEntries();
    if(printDebug)cout<<"nentries_jet55or65or80 = "<<nentries_jet55or65<<endl;
    //if(printDebug)nentries_jet55or65 = 2;
   
    for(int jentry = 0;jentry<nentries_jet55or65;jentry++){

      jetpbpb1[0][k]->GetEntry(jentry);
      jetpbpb1[1][k]->GetEntry(jentry);
      jetpbpb1[2][k]->GetEntry(jentry);    
      jetpbpb1[3][k]->GetEntry(jentry);
      jetpbpb1[4][k]->GetEntry(jentry);

      //if(printDebug && jentry%100000==0)cout<<"Jet 55or65 file"<<endl;
      if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_1<<"; run = "<<run_1<<endl;
      
      int centBin = findBin(hiBin_1);//tells us the centrality of the event. 
      
      if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 

      if(fabs(vz_1)>15) continue;
      
      int jetCounter = 0;//counts jets which are going to be used in the supernova cut rejection. 
      
      //if(algo=="Vs"){
      
      for(int j = 0;j<nbins_eta;j++){
	  
	for(int g = 0;g<nrefe_1;g++){
	    
	  if(eta_1[g]>=boundaries_eta[j][0] && eta_1[g]<boundaries_eta[j][1]){
	      
	    if(pt_1[g]>=50) jetCounter++;
	    
	  }//eta selection cut
	  
	}// jet loop
	
      }//eta bins loop
      
      // if(printDebug)cout<<"pixel hit = "<<hiNpix_1<<", jet counter = "<<jetCounter<<endl;
      
      // hpbpb_Npix_before_cut[k][centBin]->Fill(jetCounter,hiNpix_1);
      // hpbpb_Npix_before_cut[k][nbins_cent]->Fill(jetCounter,hiNpix_1);	
      
      //if(hiBin_1 >= 0 && hiBin_1 < 1) hpbpb_Npix_before_cut[k][nbins_cent+1]->Fill(jetCounter,hiNpix_1);	

      // apply the correct supernova selection cut rejection here: 
      if(hiNpix_1 > 38000 - 500*jetCounter){
       	if(printDebug) cout<<"removed this supernova event"<<endl;
      	continue;
      }
      
#if 0
      Float_t Vs_0_x_minus = sumpT[0]*v_n[0][0]*TMath::Cos(0*psi_n[0][0]);
      Float_t Vs_0_x_plus = sumpT[14]*v_n[0][14]*TMath::Cos(0*psi_n[0][14]);
      Float_t Vs_0_y_minus = sumpT[0]*v_n[0][0]*TMath::Sin(0*psi_n[0][0]);
      Float_t Vs_0_y_plus = sumpT[14]*v_n[0][14]*TMath::Sin(0*psi_n[0][14]);
      Float_t Vs_0_x = Vs_0_x_minus + Vs_0_x_plus;
      Float_t Vs_0_y = Vs_0_y_minus + Vs_0_y_plus;
      //if(printDebug)std::cout<<"Vs_0_x = "<<Vs_0_x<<"; Vs_0_y =  "<<Vs_0_y<<std::endl;
      //if(TMath::Abs(Vs_0_x)> v0_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v0 > 5200"<<std::endl;
      //if(TMath::Abs(Vs_0_y)> v0_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v0 > 5200"<<std::endl;
      if(TMath::Abs(Vs_0_x)> v0_tight || TMath::Abs(Vs_0_y)> v0_tight) hEvents_Diverge_v0_tight->Fill(1);
      if(TMath::Abs(Vs_0_x)> v0_loose || TMath::Abs(Vs_0_y)> v0_loose) hEvents_Diverge_v0_loose->Fill(1);
      
      Float_t Vs_1_x_minus = sumpT[0]*v_n[1][0]*TMath::Cos(1*psi_n[1][0]);
      Float_t Vs_1_x_plus = sumpT[14]*v_n[1][14]*TMath::Cos(1*psi_n[1][14]);
      Float_t Vs_1_y_minus = sumpT[0]*v_n[1][0]*TMath::Sin(1*psi_n[1][0]);
      Float_t Vs_1_y_plus = sumpT[14]*v_n[1][14]*TMath::Sin(1*psi_n[1][14]);
      Float_t Vs_1_x = Vs_1_x_minus + Vs_1_x_plus;
      Float_t Vs_1_y = Vs_1_y_minus + Vs_1_y_plus;
      //if(printDebug)std::cout<<"Vs_1_x = "<<Vs_1_x<<"; Vs_1_y =  "<<Vs_1_y<<std::endl;
      //if(TMath::Abs(Vs_1_x)> v1_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v1 > 100"<<std::endl;
      //if(TMath::Abs(Vs_1_y)> v1_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v1 > 100"<<std::endl;
      if(TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight) hEvents_Diverge_v1_tight->Fill(1);
      if(TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose) hEvents_Diverge_v1_loose->Fill(1);

      Float_t Vs_2_x_minus = sumpT[0]*v_n[2][0]*TMath::Cos(2*psi_n[2][0]);
      Float_t Vs_2_x_plus = sumpT[14]*v_n[2][14]*TMath::Cos(2*psi_n[2][14]);
      Float_t Vs_2_y_minus = sumpT[0]*v_n[2][0]*TMath::Sin(2*psi_n[2][0]);
      Float_t Vs_2_y_plus = sumpT[14]*v_n[2][14]*TMath::Sin(2*psi_n[2][14]);
      Float_t Vs_2_x = Vs_2_x_minus + Vs_2_x_plus;
      Float_t Vs_2_y = Vs_2_y_minus + Vs_2_y_plus;
      //if(printDebug)std::cout<<"Vs_2_x = "<<Vs_2_x<<"; Vs_2_y =  "<<Vs_2_y<<std::endl;
      //if(TMath::Abs(Vs_2_x)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v2 > 140"<<std::endl;
      //if(TMath::Abs(Vs_2_y)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v2 > 140"<<std::endl;
      if(TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight) hEvents_Diverge_v2_tight->Fill(1);
      if(TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose) hEvents_Diverge_v2_loose->Fill(1);

      Float_t Vs_3_x_minus = sumpT[0]*v_n[3][0]*TMath::Cos(3*psi_n[3][0]);
      Float_t Vs_3_x_plus = sumpT[14]*v_n[3][14]*TMath::Cos(3*psi_n[3][14]);
      Float_t Vs_3_y_minus = sumpT[0]*v_n[3][0]*TMath::Sin(3*psi_n[3][0]);
      Float_t Vs_3_y_plus = sumpT[14]*v_n[3][14]*TMath::Sin(3*psi_n[3][14]);
      Float_t Vs_3_x = Vs_3_x_minus + Vs_3_x_plus;
      Float_t Vs_3_y = Vs_3_y_minus + Vs_3_y_plus;
      //if(printDebug)std::cout<<"Vs_3_x = "<<Vs_3_x<<"; Vs_3_y =  "<<Vs_3_y<<std::endl;
      //if(TMath::Abs(Vs_3_x)>120 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v3 > 120"<<std::endl;
      //if(TMath::Abs(Vs_3_y)>120 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v3 > 120"<<std::endl;
      if(TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight) hEvents_Diverge_v3_tight->Fill(1);
      if(TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose) hEvents_Diverge_v3_loose->Fill(1);

      Float_t Vs_4_x_minus = sumpT[0]*v_n[4][0]*TMath::Cos(4*psi_n[4][0]);
      Float_t Vs_4_x_plus = sumpT[14]*v_n[4][14]*TMath::Cos(4*psi_n[4][14]);
      Float_t Vs_4_y_minus = sumpT[0]*v_n[4][0]*TMath::Sin(4*psi_n[4][0]);
      Float_t Vs_4_y_plus = sumpT[14]*v_n[4][14]*TMath::Sin(4*psi_n[4][14]);
      Float_t Vs_4_x = Vs_4_x_minus + Vs_4_x_plus;
      Float_t Vs_4_y = Vs_4_y_minus + Vs_4_y_plus;

      //if(printDebug)std::cout<<"Vs_4_x = "<<Vs_4_x<<"; Vs_4_y =  "<<Vs_4_y<<std::endl;
      //if(TMath::Abs(Vs_4_x)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v4 > 140"<<std::endl;
      //if(TMath::Abs(Vs_4_y)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v4 > 140"<<std::endl;
      if(TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) hEvents_Diverge_v4_tight->Fill(1);
      if(TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) hEvents_Diverge_v4_loose->Fill(1);

      if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
	hEvents_Diverge_tight->Fill(1);
	isDivergeTight = true;
      }
      if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
	hEvents_Diverge_loose->Fill(1);
	isDivergeLoose = true;
      }
      

#endif

      for(int j = 0;j<nbins_eta;j++){

	for(int g = 0;g<nrefe_1;g++){ // this is the loop for the  Jets we are interested in.  
	  
	  //if((chMax_1[g]/pt_1[g]>0.05) && (chMax_1[g]/pt_1[g]<1.00)){
	  
	  if(eta_1[g]<boundaries_eta[j][0] || eta_1[g]>=boundaries_eta[j][1]) continue;
	  
	  //if(raw_1[g] < 20) continue;
	  
	  // if(printDebug)cout<<"jet variables:  "<<endl;
	  // if(printDebug)cout<<"chSum = "<<chSum_1[g]<<endl;
	  // if(printDebug)cout<<"phSum = "<<phSum_1[g]<<endl;
	  // if(printDebug)cout<<"neSum = "<<neSum_1[g]<<endl;
	  // if(printDebug)cout<<"muSum = "<<muSum_1[g]<<endl;
	  // if(printDebug)cout<<"eSum  = "<<eSum_1[g]<<endl;
	  // if(printDebug)cout<<"jtpt  = "<<pt_1[g]<<endl;
	  // if(printDebug)cout<<"rawpT = "<<raw_1[g]<<endl;
	  // if(printDebug)cout<<"jtpu  = "<<jtpu_1[g]<<endl;
	  // cut3 = (float)(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g] + eSum_1[g] - jtpu_1[g])/(pt_1[g]);
	  // cut1 = (float)chMax_1[g]/(pt_1[g]);
	  // cut2 = (float)neMax_1[g]/TMath::Max(chSum_1[g],neSum_1[g]);
	  // cut2b = (float)phMax_1[g]/TMath::Max(chSum_1[g],neSum_1[g]);
	  // cut4 = (float)(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g] + eSum_1[g])/(0.5*raw_1[g]);
	  // cut5 = (float)neMax_1[g]/(neMax_1[g] + chMax_1[g] + phMax_1[g]);
	  // cut6 = (float)muMax_1[g]/(neMax_1[g] + chMax_1[g] + phMax_1[g]);
	  // cut7 = (float)eMax_1[g]/(neMax_1[g] + chMax_1[g] + phMax_1[g]);
	  // if(printDebug)cout<<"Cut1 value = "<<cut1<<endl;
	  // if(printDebug)cout<<"Cut2 value = "<<cut2<<endl;
	  // if(printDebug)cout<<"Cut3 value = "<<cut3<<endl;
	  // if(printDebug)cout<<"Cut4 value = "<<cut4<<endl;
	  // if(printDebug)cout<<"Cut4 value = "<<cut5<<endl;
	  
	  //hCut3->Fill(cut3);
	  //hCut5->Fill(cut5);
	  //hCut4->Fill(cut4);
	  //hCut2->Fill(cut2);
	  //hCut1->Fill(cut1);
	  
	  
	  arrayValues[0] = raw_1[g];
	  arrayValues[1] = pt_1[g];
	  arrayValues[2] = jtpu_1[g];
	  arrayValues[3] = L1_sj36_1;
	  arrayValues[4] = L1_sj36_p_1;
	  arrayValues[5] = L1_sj52_1;
	  arrayValues[6] = L1_sj52_p_1;
	  arrayValues[7] = jet55_1;
	  arrayValues[8] = jet55_p_1;
	  arrayValues[9] = jet65_1;
	  arrayValues[10] = jet65_p_1;
	  arrayValues[11] = jet80_1;
	  arrayValues[12] = jet80_p_1;
	  arrayValues[13] = trgObj_pt_1;
	  arrayValues[14] = centBin;
	  arrayValues[15] = chMax_1[g];
	  arrayValues[16] = chSum_1[g];
	  arrayValues[17] = phMax_1[g];
	  arrayValues[18] = phSum_1[g];
	  arrayValues[19] = neMax_1[g];
	  arrayValues[20] = neSum_1[g];
	  arrayValues[21] = muMax_1[g];
	  arrayValues[22] = muSum_1[g];
	  arrayValues[23] = eMax_1[g];
	  arrayValues[24] = eSum_1[g];

	  jets_ID->Fill(arrayValues);
	
	  // going to use the effective prescl for the Jet55 trigger: 2.0475075465
	  Float_t effecPrescl = 2.047507;

#if 0

	  //if(cut5<0.95 && cut6<0.95 && cut7<0.95 && cut1>0.05 && cut2<0.95 && cut2b<0.95){
	  if(1>0){

	    //	    hJets->Fill(1);
	    //	    if(raw_1[g] < 20) continue;

	    if(jet80_1==1 && L1_sj52_1==1) {
	      if(trgObj_pt_1>=80){
		hpbpb_TrgObj80[k][j][centBin]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		hpbpb_TrgObj80[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		
		if(jetCounter>=7){
		  hpbpb_TrgObj80_nJet_g7[k][j][centBin]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		  hpbpb_TrgObj80_nJet_g7[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		}
		if(jetCounter<7){
		  hpbpb_TrgObj80_nJet_l7[k][j][centBin]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		  hpbpb_TrgObj80_nJet_l7[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		}
	      }
	      
	      if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
		hpbpb_Jet80[1][1][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[1][1][nbins_cent]->Fill(pt_1[g]);
	      }
	      if(TMath::Abs(Vs_0_x)<v0_tight && TMath::Abs(Vs_0_y)<v0_tight && TMath::Abs(Vs_1_x)<v1_tight && TMath::Abs(Vs_1_y)<v1_tight && TMath::Abs(Vs_2_x)<v2_tight && TMath::Abs(Vs_2_y)<v2_tight && TMath::Abs(Vs_3_x)<v3_tight && TMath::Abs(Vs_3_y)<v3_tight && TMath::Abs(Vs_4_x)<v4_tight && TMath::Abs(Vs_4_y)<v4_tight){
		hpbpb_Jet80[0][1][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[0][1][nbins_cent]->Fill(pt_1[g]);
	      }

	      if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
		hpbpb_Jet80[1][0][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[1][0][nbins_cent]->Fill(pt_1[g]);
	      }
	      if(TMath::Abs(Vs_0_x)<v0_loose && TMath::Abs(Vs_0_y)<v0_loose && TMath::Abs(Vs_1_x)<v1_loose && TMath::Abs(Vs_1_y)<v1_loose && TMath::Abs(Vs_2_x)<v2_loose && TMath::Abs(Vs_2_y)<v2_loose && TMath::Abs(Vs_3_x)<v3_loose && TMath::Abs(Vs_3_y)<v3_loose && TMath::Abs(Vs_4_x)<v4_loose && TMath::Abs(Vs_4_y)<v4_loose){
		hpbpb_Jet80[0][0][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[0][0][nbins_cent]->Fill(pt_1[g]);
	      }

	      hpbpb_FullJet80[centBin]->Fill(pt_1[g]);
	      hpbpb_FullJet80[nbins_cent]->Fill(pt_1[g]);

	    }else if(jet65_1==1 && L1_sj36_1==1) {


	      // if(trgObj_pt_1>=65 && trgObj_pt_1<80){
	      // 	hpbpb_TrgObj65[k][j][centBin]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	hpbpb_TrgObj65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
		
	      // 	if(jetCounter>=7){
	      // 	  hpbpb_TrgObj65_nJet_g7[k][j][centBin]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	  hpbpb_TrgObj65_nJet_g7[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	}
	      // 	if(jetCounter<7){
	      // 	  hpbpb_TrgObj65_nJet_l7[k][j][centBin]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	  hpbpb_TrgObj65_nJet_l7[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	}
	      // }
	      // hpbpb_Jet65[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
	      // hpbpb_Jet65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);

	      if(jet80_1==0){
		
		if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
		  hpbpb_Jet65[1][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[1][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_tight && TMath::Abs(Vs_0_y)<v0_tight && TMath::Abs(Vs_1_x)<v1_tight && TMath::Abs(Vs_1_y)<v1_tight && TMath::Abs(Vs_2_x)<v2_tight && TMath::Abs(Vs_2_y)<v2_tight && TMath::Abs(Vs_3_x)<v3_tight && TMath::Abs(Vs_3_y)<v3_tight && TMath::Abs(Vs_4_x)<v4_tight && TMath::Abs(Vs_4_y)<v4_tight){
		  hpbpb_Jet65[0][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[0][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
		  hpbpb_Jet65[1][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[1][0][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_loose && TMath::Abs(Vs_0_y)<v0_loose && TMath::Abs(Vs_1_x)<v1_loose && TMath::Abs(Vs_1_y)<v1_loose && TMath::Abs(Vs_2_x)<v2_loose && TMath::Abs(Vs_2_y)<v2_loose && TMath::Abs(Vs_3_x)<v3_loose && TMath::Abs(Vs_3_y)<v3_loose && TMath::Abs(Vs_4_x)<v4_loose && TMath::Abs(Vs_4_y)<v4_loose){
		  hpbpb_Jet65[0][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[0][0][nbins_cent]->Fill(pt_1[g]);
		}
	      
		hpbpb_FullJet65[centBin]->Fill(pt_1[g]);
		hpbpb_FullJet65[nbins_cent]->Fill(pt_1[g]);
	      }

	    }else if(jet55_1==1 && L1_sj36_1==1) { // passes the jet55 trigger
#if 0
	      if(trgObj_pt_1>=55 && trgObj_pt_1<65){ // check for the trigger object pt to lie inbetween the two trigger values 
		hpbpb_TrgObj55[k][j][centBin]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		hpbpb_TrgObj55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		
		if(jetCounter>=7){
		  hpbpb_TrgObj55_nJet_g7[k][j][centBin]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		  hpbpb_TrgObj55_nJet_g7[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		}
		if(jetCounter<7){
		  hpbpb_TrgObj55_nJet_l7[k][j][centBin]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		  hpbpb_TrgObj55_nJet_l7[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		}
		//if(pt_1[g]>3*trgObj_pt_1){ // get an idea of the event information from these large pt jets 
		//  fHLT_high[k]<<evt_1<<" "<<run_1<<" "<<lumi_1<<" "<<vz_1<<" "<<centBin<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<pt_1[g]<<" "<<raw_1[g]<<" "<<eta_1[g]<<" "<<phi_1[g]<<" "<<endl;
		//}
	      }
	      // hpbpb_Jet55[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
	      // hpbpb_Jet55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
#endif
	      if((jet65_1==0) && (jet80_1==0)){ // this is to just check

		if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
		  hpbpb_Jet55[1][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[1][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_tight && TMath::Abs(Vs_0_y)<v0_tight && TMath::Abs(Vs_1_x)<v1_tight && TMath::Abs(Vs_1_y)<v1_tight && TMath::Abs(Vs_2_x)<v2_tight && TMath::Abs(Vs_2_y)<v2_tight && TMath::Abs(Vs_3_x)<v3_tight && TMath::Abs(Vs_3_y)<v3_tight && TMath::Abs(Vs_4_x)<v4_tight && TMath::Abs(Vs_4_y)<v4_tight){
		  hpbpb_Jet55[0][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[0][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
		  hpbpb_Jet55[1][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[1][0][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_loose && TMath::Abs(Vs_0_y)<v0_loose && TMath::Abs(Vs_1_x)<v1_loose && TMath::Abs(Vs_1_y)<v1_loose && TMath::Abs(Vs_2_x)<v2_loose && TMath::Abs(Vs_2_y)<v2_loose && TMath::Abs(Vs_3_x)<v3_loose && TMath::Abs(Vs_3_y)<v3_loose && TMath::Abs(Vs_4_x)<v4_loose && TMath::Abs(Vs_4_y)<v4_loose){
		  hpbpb_Jet55[0][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[0][0][nbins_cent]->Fill(pt_1[g]);
		}
	      
		hpbpb_FullJet55[centBin]->Fill(pt_1[g]);
		hpbpb_FullJet55[nbins_cent]->Fill(pt_1[g]);
	      
	      }
	    }

	    // if(jet80_1) {
	    //   hpbpb_Jet80[k][j][centBin]->Fill(pt_1[g],jet80_p_1);
	    //   hpbpb_Jet80[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1);
	    // }//Jet80 trigger selection

	    // if(jet65_1) {
	    //   hpbpb_Jet65[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
	    //   hpbpb_Jet65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
	    // }//Jet65 trigger selection

	    // if(jet55_1) {
	    //   hpbpb_Jet55[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
	    //   hpbpb_Jet55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
	    // }//Jet55 trigger selection

	    // if(jet65_1 && !jet80_1){
	    //   hpbpb_Jet65_noJet80[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
	    //   hpbpb_Jet65_noJet80[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
	    // }//jet65 and not Jet80 selection 

	    // if(jet55_1 && !jet65_1 && !jet80_1){
	    //   hpbpb_Jet55_noJet65_80[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
	    //   hpbpb_Jet55_noJet65_80[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
	    // }//jet55 and not Jet65 and not Jet80 selection 

	      
	  }//qa cut selection
#endif
	    
	}//jet loop
	  
      }//eta bin loop
    
    }//nentries_jet55or65 loop
    
    
    // //loop for jetpbpb2[2] tree
    //   Long64_t nentries_jet80or95 = jetpbpb2[2][k]->GetEntries();
    //   if(printDebug)cout<<"nentries_jet80or95 = "<<nentries_jet80or95<<endl;
    //   if(printDebug)nentries_jet80or95 = 2;

    //   for(int jentry = 0;jentry<nentries_jet80or95;jentry++){
    
    //     jetpbpb2[2][k]->GetEntry(jentry);
    //     if(printDebug && jentry%100000==0)cout<<"Jet 80or95 file"<<endl;
    //     if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_2<<", run = "<<run_2<<endl;

    //     //
    //     int centBin = findBin(hiBin_2);//tells us the centrality of the event
    //     if(printDebug)cout<<"centrality bin = "<<5*boundaries_cent[centBin]<< " to "<<5*boundaries_cent[centBin+1]<<endl;

    //     for(int j = 0;j<nbins_eta;j++){

    //       if(pHBHENoiseFilter_2 && pcollisionEventSelection_2 && (vz_2<15)) {
    
    //         for(int g = 0;g<nrefe_2;g++){

    //           if(/*put your favourite QA cut here*/chMax_2[g]/pt_2[g]>0.01){

    //             if(eta_2[g]>=boundaries_eta[j][0] && eta_2[g]<boundaries_eta[j][1]){

    //               if(jet80_2 && trgObj_pt_2>=80){

    //                 hpbpb_TrgObj80[k][j][centBin]->Fill(pt_2[g],jet80_p_2);
    //                 hpbpb_TrgObj80[k][j][nbins_cent]->Fill(pt_2[g],jet80_p_2);

    //               }//80 trg obj selection

    //             }//eta selection

    //           }//qa cut selection

    //         }//jet loop

    //       }//event selection cuts

    //     }//eta bin loop  
    
    //   }//nentries_jet80or95 loop

    //fVs_failure[k].close();
    
  }//radius loop. 
 
 
  /*
  //declare the output file
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_data_ak%s%s_12003_effPres_severalcut_%d_endfile_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  //TFile f(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_data_ak%s%s_%d_endfile_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  
  f.cd();
  
  for(int k = 0;k<no_radius;k++){
    
  for(int j = 0;j<nbins_eta;j++){
      
  for(int i = 0;i<=nbins_cent;i++){
	
  //divide by delta eta, delta pt, luminosity seen by the triggers and ncoll ... and that will give us cross section/ncoll. 
  // lumi seen by the respective triggers. - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HIJetRAA#Luminosity_Cross_check_for_HIJet
  // once we include the prescl then we have to normalize w.r.t the hightest lumi trigger which would have a prescl of 1 - in this case its the jet80 trigger. 

  //not dividing by ncoll here. will do that in the plotting macro under RAA. 
	
  //hpbpb_TrgObj80[k][j][i]->Scale(1./(delta_eta[j]*145.156*1e6));
  //hpbpb_TrgObj65[k][j][i]->Scale(1./(delta_eta[j]*145.156*1e6));
  //hpbpb_TrgObj55[k][j][i]->Scale(1./(delta_eta[j]*145.156*1e6));
	
  // divide by bin width is doing something very weird. Its making histograms which are empty get 40 entries from somewhere - fixed that issue due to empty bins being divided by bin width. 

  divideBinWidth(hpbpb_TrgObj80[k][j][i]);
  divideBinWidth(hpbpb_TrgObj65[k][j][i]);
  divideBinWidth(hpbpb_TrgObj55[k][j][i]);
	
  hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj80[k][j][i]);
  hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj65[k][j][i]);
  hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj55[k][j][i]);

  hpbpb_TrgObjComb_nJet_g7[k][j][i]->Add(hpbpb_TrgObj80_nJet_g7[k][j][i]);
  hpbpb_TrgObjComb_nJet_g7[k][j][i]->Add(hpbpb_TrgObj65_nJet_g7[k][j][i]);
  hpbpb_TrgObjComb_nJet_g7[k][j][i]->Add(hpbpb_TrgObj55_nJet_g7[k][j][i]);

  hpbpb_TrgObjComb_nJet_l7[k][j][i]->Add(hpbpb_TrgObj80_nJet_l7[k][j][i]);
  hpbpb_TrgObjComb_nJet_l7[k][j][i]->Add(hpbpb_TrgObj65_nJet_l7[k][j][i]);
  hpbpb_TrgObjComb_nJet_l7[k][j][i]->Add(hpbpb_TrgObj55_nJet_l7[k][j][i]);
	
  // old method histograms, just as a test: 
	
  // hpbpb_Jet80[k][j][i]->Scale(1./delta_eta[j]);
  // hpbpb_Jet65[k][j][i]->Scale(1./delta_eta[j]);
  // hpbpb_Jet55[k][j][i]->Scale(1./delta_eta[j]);

  // hpbpb_Jet65_noJet80[k][j][i]->Scale(1./delta_eta[j]);
  // hpbpb_Jet55_noJet65_80[k][j][i]->Scale(1./delta_eta[j]);
	
  divideBinWidth(hpbpb_Jet80[k][j][i]);
  // divideBinWidth(hpbpb_Jet65[k][j][i]);
  // divideBinWidth(hpbpb_Jet55[k][j][i]);

  divideBinWidth(hpbpb_Jet65_noJet80[k][j][i]);
  divideBinWidth(hpbpb_Jet55_noJet65_80[k][j][i]);
	
  hpbpb_JetComb_old[k][j][i]->Add(hpbpb_Jet80[k][j][i]);
  hpbpb_JetComb_old[k][j][i]->Add(hpbpb_Jet65_noJet80[k][j][i]);
  hpbpb_JetComb_old[k][j][i]->Add(hpbpb_Jet55_noJet65_80[k][j][i]);
	
  hpbpb_TrgObjComb[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObjComb[k][j][i]->Print();
  hpbpb_TrgObj80[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj80[k][j][i]->Print();
  hpbpb_TrgObj65[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj65[k][j][i]->Print();
  hpbpb_TrgObj55[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj55[k][j][i]->Print();

  hpbpb_TrgObjComb_nJet_g7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObjComb_nJet_g7[k][j][i]->Print();
  hpbpb_TrgObj80_nJet_g7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj80_nJet_g7[k][j][i]->Print();
  hpbpb_TrgObj65_nJet_g7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj65_nJet_g7[k][j][i]->Print();
  hpbpb_TrgObj55_nJet_g7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj55_nJet_g7[k][j][i]->Print();

  hpbpb_TrgObjComb_nJet_l7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObjComb_nJet_l7[k][j][i]->Print();
  hpbpb_TrgObj80_nJet_l7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj80_nJet_l7[k][j][i]->Print();
  hpbpb_TrgObj65_nJet_l7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj65_nJet_l7[k][j][i]->Print();
  hpbpb_TrgObj55_nJet_l7[k][j][i]->Write();
  if(printDebug)hpbpb_TrgObj55_nJet_l7[k][j][i]->Print();	

  hpbpb_Jet80[k][j][i]->Write();
  if(printDebug)hpbpb_Jet80[k][j][i]->Print();
  // hpbpb_Jet65[k][j][i]->Write();
  // if(printDebug)hpbpb_Jet65[k][j][i]->Print();
  // hpbpb_Jet55[k][j][i]->Write();
  // if(printDebug)hpbpb_Jet55[k][j][i]->Print();

  hpbpb_Jet65_noJet80[k][j][i]->Write();
  if(printDebug)hpbpb_Jet65_noJet80[k][j][i]->Print();
  hpbpb_Jet55_noJet65_80[k][j][i]->Write();
  if(printDebug)hpbpb_Jet55_noJet65_80[k][j][i]->Print();
  hpbpb_JetComb_old[k][j][i]->Write();
  if(printDebug)hpbpb_JetComb_old[k][j][i]->Print();

  }//cent bins loop

  }//eta bins loop

  for(int i = 0;i<=nbins_cent;i++){

  if(printDebug)hpbpb_Npix_before_cut[k][i]->Print("base");
  hpbpb_Npix_before_cut[k][i]->Write();
  if(printDebug)hpbpb_Npix_after_cut[k][i]->Print("base");
  hpbpb_Npix_after_cut[k][i]->Write();
      
  }
  if(printDebug)hpbpb_Npix_before_cut[k][nbins_cent+1]->Print("base");
  hpbpb_Npix_before_cut[k][nbins_cent+1]->Write();

  hpbpb_cent[k]->Write();
  if(printDebug)hpbpb_cent[k]->Print("base");
  hpbpb_vz[k]->Write();
  if(printDebug)hpbpb_vz[k]->Print("base");
  hpbpb_vx[k]->Write();
  if(printDebug)hpbpb_vx[k]->Print("base");
  hpbpb_vy[k]->Write();
  if(printDebug)hpbpb_vy[k]->Print("base");

  }//radius loop

  //hCut1->Write();
  //hCut2->Write();
  //hCut3->Write();
  //hCut4->Write();
  //hCut5->Write();

  hEvents->Write();
  hJets->Write();
  //jets_ID->Write();

  //if(printDebug)jets_ID->Print();
  if(printDebug)hEvents->Print();
  if(printDebug)hJets->Print();

  f.Write();

  f.Close();
  */
  

  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_jetntuple_withEvtCuts_SuperNovaRejected_ak%s3%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  f.cd();

  jets_ID->Write();


  // hEvents->Write();
  // hEvents_Diverge_tight->Write(); 
  // hEvents_Diverge_loose->Write(); 
  // hEvents_Diverge_v0_tight->Write(); 
  // hEvents_Diverge_v0_loose->Write(); 
  // hEvents_Diverge_v1_tight->Write(); 
  // hEvents_Diverge_v1_loose->Write();
  // hEvents_Diverge_v2_tight->Write(); 
  // hEvents_Diverge_v2_loose->Write();
  // hEvents_Diverge_v3_tight->Write(); 
  // hEvents_Diverge_v3_loose->Write();
  // hEvents_Diverge_v4_tight->Write(); 
  // hEvents_Diverge_v4_loose->Write();

  // for(int i = 0;i<nbins_cent+1;i++){
  
  //   for(int a = 0;a<2;a++){

  //     for(int b = 0;b<2;b++){

  // 	hpbpb_Jet80[b][a][i]->Write();
  // 	if(printDebug)hpbpb_Jet80[b][a][i]->Print();
  // 	hpbpb_Jet65[b][a][i]->Write();
  // 	if(printDebug)hpbpb_Jet65[b][a][i]->Print();
  // 	hpbpb_Jet55[b][a][i]->Write();
  // 	if(printDebug)hpbpb_Jet55[b][a][i]->Print();		

  // 	hpbpb_Pu_Jet80[b][a][i]->Write();
  // 	if(printDebug)hpbpb_Pu_Jet80[b][a][i]->Print();
  // 	hpbpb_Pu_Jet65[b][a][i]->Write();
  // 	if(printDebug)hpbpb_Pu_Jet65[b][a][i]->Print();
  // 	hpbpb_Pu_Jet55[b][a][i]->Write();
  // 	if(printDebug)hpbpb_Pu_Jet55[b][a][i]->Print();

  //     }
      
  //   }

  //   hpbpb_FullJet80[i]->Write();
  //   if(printDebug)hpbpb_FullJet80[i]->Print();
  //   hpbpb_FullJet65[i]->Write();
  //   if(printDebug)hpbpb_FullJet65[i]->Print();
  //   hpbpb_FullJet55[i]->Write();
  //   if(printDebug)hpbpb_FullJet55[i]->Print();
  
  //   hpbpb_Pu_FullJet80[i]->Write();
  //   if(printDebug)hpbpb_Pu_FullJet80[i]->Print();
  //   hpbpb_Pu_FullJet65[i]->Write();
  //   if(printDebug)hpbpb_Pu_FullJet65[i]->Print();
  //   hpbpb_Pu_FullJet55[i]->Write();
  //   if(printDebug)hpbpb_Pu_FullJet55[i]->Print();
    
  // }

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;

}
