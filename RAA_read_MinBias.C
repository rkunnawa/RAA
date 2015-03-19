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
#define pi 3.14159265

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
static Double_t ncoll[nbins_cent+1] = { 1660, 1310, 745, 251, 62.8, 10.8, 362.24}; //last one is for 0-200 bin. 

static const int no_radius = 3;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {2,3,4};

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
  int binNo = -1;

  for(int i = 0;i<nbins_cent;i++){
    if(hiBin>=5*boundaries_cent[i] && hiBin<5*boundaries_cent[i+1]) {
      binNo = i;
      break;
    }
  }

  return binNo;
}

double Calc_deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  double dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}

using namespace std;

void RAA_read_MinBias(int startfile = 0, int endfile = 1, char *algo = "Pu", char *jet_type = "PF", char *Type = "Data"){

  TH1::SetDefaultSumw2();
  
  TStopwatch timer;
  timer.Start();
  
  TDatime date;
  
  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = true;

  // Now im going to change the file reading here for PbPb to look at the unmerged files through condor. 
  std::string infile1;
  if(Type=="MC")infile1 = "PbPb_HydjetMinBias_forest.txt";
  //if(Type=="Data")infile1 = "PbPb_MinBiasUPC_forest.txt";
  if(Type=="Data")infile1 = "jetRAA_MinBiasUPC_forest.txt";
  //if(Type=="Data")infile1 = "14010_MinBiasUPC_forest.txt";
  
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
    dir[4][k] = "hltobject"; //- not there in MinBias files 
    //dir[6][k] = "pfcandAnalyzer";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    //"t",
    "HiTree",
    "jetObjTree"
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
  float refpt_1[1000];
  float subid_1[1000];
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
  int l1MB_1;
  int L1_MB_p_1;
  int l1MB_p_1;
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
    if(Type=="MC")jetpbpb1[2][k]->SetBranchAddress("refpt",&refpt_1);
    if(Type=="MC")jetpbpb1[2][k]->SetBranchAddress("subid",&subid_1);
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
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1",&jetMB_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1_Prescl",&jetMB_p_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
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

    jetpbpb1[4][k]->SetBranchAddress("pt",&trgObj_pt_1);
    jetpbpb1[4][k]->SetBranchAddress("eta",&trgObj_eta_1);
    jetpbpb1[4][k]->SetBranchAddress("phi",&trgObj_phi_1);
    
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
  /*
  TNtuple *jets_ID;
  if(Type=="Data") jets_ID= new TNtuple("jets_ID","","rawpt:jtpt:jtpu:L1_MB:L1_MB_prescl:HLTZeroBias_MB:HLTZeroBias_MB_prescl:cent:chMax:chSum:phMax:phSum:neMax:neSum:muMax:muSum:eMax:eSum");
  if(Type=="MC")  jets_ID= new TNtuple("jets_ID","","rawpt:jtpt:jtpu:L1_MB:L1_MB_prescl:HLTZeroBias_MB:HLTZeroBias_MB_prescl:cent:chMax:chSum:phMax:phSum:neMax:neSum:muMax:muSum:eMax:eSum:refpt:subid");
  Float_t arrayValues_Data[18];
  Float_t arrayValues_MC[20];
  */
  
  TH1F *hEvents  = new TH1F("hEvents","",2,0,1);
  TH1F *hJets = new TH1F("hJets","",2,0,1);

  TH1F *hJetMB_80_1 = new TH1F("hJetMB_80_1","",200,0,200);
  TH1F *hJetMB_65_1 = new TH1F("hJetMB_65_1_","",200,0,200);
  TH1F *hJetMB_55_1 = new TH1F("hJetMB_55_1","",200,0,200);
  TH1F *hJetMBSpectra[no_radius][nbins_cent];
  for(int k = 0;k<no_radius;++k){
    for(int i = 0;i<nbins_cent;++i)
      hJetMBSpectra[k][i] = new TH1F(Form("hJetMBSpectra_R%d_cent%d",list_radius[k],i),"Data from MB trigger alone to add to the Jet triggered data",1000,0,1000);
  }
  
  TH1F *hL1MB = new TH1F("hL1MB","",200,0,200);
  TH1F *hJet80 = new TH1F("hJet80","",200,0,200);
  TH1F *hJet80_Trig = new TH1F("hJet80_Trig","",200,0,200);
  TH1F *hJet65 = new TH1F("hJet65","",200,0,200);
  TH1F *hJet55 = new TH1F("hJet55","",200,0,200);
  TH1F *hL1SJ36 = new TH1F("hL1SJ36","",200,0,200);
  TH1F *hL1SJ52 = new TH1F("hL1SJ52","",200,0,200);
  TH1F *hL1MB_JetMB = new TH1F("hL1MB_JetMB","",200,0,200);
  TH1F *hJet80_JetMB = new TH1F("hJet80_JetMB","",200,0,200);
  TH1F *hJet65_JetMB = new TH1F("hJet65_JetMB","",200,0,200);
  TH1F *hJet55_JetMB = new TH1F("hJet55_JetMB","",200,0,200);
  TH1F *hJet65_Jet80 = new TH1F("hJet65_Jet80","",200,0,200);
  TH1F *hJet55_Jet80 = new TH1F("hJet55_Jet80","",200,0,200);
  TH1F *hL1SJ36_JetMB = new TH1F("hL1SJ36_JetMB","",200,0,200);
  TH1F *hL1SJ52_JetMB = new TH1F("hL1SJ52_JetMB","",200,0,200);

  TH1F * hJet80_Prescl = new TH1F("hJet80_Prescl","",200,0,10);
  TH1F * hJet65_Prescl = new TH1F("hJet65_Prescl","",200,0,10);
  TH1F * hJet55_Prescl = new TH1F("hJet55_Prescl","",200,0,10);
  
    
  // jetpbpb1[2][0]->Draw("jtpt[0]>>hJet80_JetMB","(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&HLT_HIJet80_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&HLT_HIJet80_v1_Prescl==1)","goff");
  // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet80_Trig","(pcollisionEventSelection&&HLT_HIJet80_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&jetObjTree.pt>80)","goff");
  // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet65_JetMB","(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&HLT_HIJet65_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&HLT_HIJet65_v1_Prescl==1)","goff");
  // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet55_JetMB","(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&HLT_HIJet55_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&HLT_HIJet55_v1_Prescl==1)","goff");
  // // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet65_Jet80","(pcollisionEventSelection&&HLT_HIJet65_v1&&HLT_HIJet80_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&pt>80)","goff");
  // // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet55_Jet80","(pcollisionEventSelection&&HLT_HIJet80_v1&&HLT_HIJet55_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&pt>80)","goff");
  // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJetMB_80_1","(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&HLT_HIJet80_v1_Prescl==1)","goff");
  // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJetMB_65_1","(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&HLT_HIJet65_v1_Prescl==1)","goff");
  // jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJetMB_55_1","(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&HLT_HIJet55_v1_Prescl==1)","goff");

  // for(int k = 0;k<no_radius;++k){

  // for(int i = 0;i<nbins_cent;i++){
  //   //cout<<Form("%f %f",5*boundaries_cent[i],5*boundaries_cent[i+1])<<endl;
  //   //cout<<Form("38.695*(pcollisionEventSelection && HLT_HIMinBiasHfOrBSC_v1 && abs(jteta)<2 && abs(vz)<15 && pHBHENoiseFilter && !HLT_HIJet80_v1 && !HLT_HIJet65_v1 && !HLT_HIJet55_v1 && hiBin>=%f && hiBin<%f)",5*boundaries_cent[i],5*boundaries_cent[i+1])<<endl;
  //   jetpbpb1[2][k]->Draw(Form("jtpt>>hJetMBSpectra_R%d_cent%d",list_radius[k],i),Form("38.695*(pcollisionEventSelection && HLT_HIMinBiasHfOrBSC_v1 && abs(jteta)<2 && abs(vz)<15 && pHBHENoiseFilter && !HLT_HIJet80_v1 && !HLT_HIJet65_v1 && !HLT_HIJet55_v1 && hiBin>=%f && hiBin<%f && chargedMax/jtpt>0.02 && eMax/jtpt<0.6)",5*boundaries_cent[i],5*boundaries_cent[i+1]),"goff");
  //   //hJetMBSpectra[k][i]->Print("base");
  // }
    
  // TTree* trigger_info;
  // Int_t evt_value;
  // Int_t run_value;
  // Int_t lumi_value;
  // Int_t nref;
  // Float_t jetpt[1000];
  // Float_t jeteta[1000];
  // Int_t Jet80;
  // Int_t Jet80_prescl;
  // Int_t Jet65;
  // Int_t Jet65_prescl;
  // Int_t Jet55;
  // Int_t Jet55_prescl;
  // Int_t HIMinBias;
  // Int_t HIMinBias_prescl;
  
  // trigger_info = new TTree("trigger_info","");
  // trigger_info->Branch("run_value",&run_value,"run_value/I");
  // trigger_info->Branch("evt_value",&evt_value,"evt_value/I");
  // trigger_info->Branch("lumi_value",&lumi_value,"lumi_value/I");
  // trigger_info->Branch("nref",&nref,"nref/I");
  // trigger_info->Branch("jetpt",&jetpt,"jetpt[nref]/F");
  // trigger_info->Branch("jeteta",&jeteta,"jeteta[nref]/F");
  // trigger_info->Branch("Jet80",&Jet80,"Jet80/I");
  // trigger_info->Branch("Jet80_prescl",&Jet80_prescl,"Jet80_prescl/I");
  // trigger_info->Branch("Jet65",&Jet65,"Jet65/I");
  // trigger_info->Branch("Jet65_prescl",&Jet65_prescl,"Jet65_prescl/I");
  // trigger_info->Branch("Jet55",&Jet55,"Jet55/I");
  // trigger_info->Branch("Jet55_prescl",&Jet55_prescl,"Jet55_prescl/I");
  // trigger_info->Branch("HIMinBias",&HIMinBias,"HIMinBias/I");
  // trigger_info->Branch("HIMinBias_prescl",&HIMinBias_prescl,"HIMinBias_prescl/I");  
  
  // prescl is calculated as the ratio of the 
  //minbias prescl = 38.695
  // jet80  prescl = 1 (same in minbias sample as well)
  // jet65  prescl = 1.11398 (1.11287)
  // jet55  prescl = 2.0475 (2.0292)

  
  TH1F* hDenominator[no_radius]; 
  TH1F* hNumerator_80[no_radius];
  TH1F* hNumerator_65[no_radius];
  TH1F* hNumerator_55[no_radius];

  TH1F * hNumerator_55_n65_n80[no_radius];
  TH1F * hNumerator_65_n80[no_radius];

  for(int k = 0;k<no_radius;k++){

    hDenominator[k] =new TH1F(Form("hDenominator_R%d",list_radius[k]),"",nbins_pt,boundaries_pt);
    hNumerator_80[k] =new TH1F(Form("hNumerator_80_R%d",list_radius[k]),"",nbins_pt,boundaries_pt);
    hNumerator_65[k] =new TH1F(Form("hNumerator_65_R%d",list_radius[k]),"",nbins_pt,boundaries_pt);
    hNumerator_55[k] =new TH1F(Form("hNumerator_55_R%d",list_radius[k]),"",nbins_pt,boundaries_pt);
    hNumerator_55_n65_n80[k] =new TH1F(Form("hNumerator_55_n65_n80_R%d",list_radius[k]),"",nbins_pt,boundaries_pt);
    hNumerator_65_n80[k] =new TH1F(Form("hNumerator_65_n80_R%d",list_radius[k]),"",nbins_pt,boundaries_pt);

  }
  
  // but these prescl values are derived from the triggered data sample, need to check if they are the same in the minbias sample 
  // cout<<"total number of events with minbias trigger = "<<jetpbpb1[2][0]->GetEntries("HLT_HIMinBiasHfOrBSC_v1 && HLT_HIJet80_v1_Prescl==1")<<endl;
  // cout<<" numerator   "<<jetpbpb1[2][0]->GetEntries("HLT_HIMinBiasHfOrBSC_v1 && HLT_HIJet80_v1_Prescl==1 && HLT_HIJet80_v1 && jtpt[0]>110 && TMath::Abs(jteta[0])<2")<<endl;
  // cout<<" denominator "<<jetpbpb1[2][0]->GetEntries("HLT_HIMinBiasHfOrBSC_v1 && HLT_HIJet80_v1_Prescl==1 && jtpt[0]>110 && TMath::Abs(jteta[0])<2")<<endl;

  for(int k = 0;k<no_radius;k++){

    if(printDebug)cout<<"Running data reading for R = "<<list_radius[k]<<endl;
    // loop for the jetpbpb1[2] tree 
    Long64_t nentries_MB = jetpbpb1[2][k]->GetEntries();
    if(printDebug)cout<<"nentries_MB = "<<nentries_MB<<endl;
    
    for(int jentry = 0;jentry<nentries_MB;jentry++){

      jetpbpb1[0][k]->GetEntry(jentry);
      jetpbpb1[1][k]->GetEntry(jentry);
      jetpbpb1[2][k]->GetEntry(jentry);    
      jetpbpb1[3][k]->GetEntry(jentry);
      //jetpbpb1[4][k]->GetEntry(jentry);

      //if(printDebug && jentry%100000==0)cout<<"Jet 55or65 file"<<endl;
      if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_1<<"; run = "<<run_1<<endl;
      
      int centBin = findBin(hiBin_1);//tells us the centrality of the event. 

      if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 

      if(fabs(vz_1)>15) continue;
      
      if(jetMB_1==0 || jet80_p_1==0 || jet65_p_1 || jet55_p_1)  continue;

      hEvents->Fill(1);

      if(nrefe_1 > 40) continue;
      
      // Float_t largepT = pt_1[0];
      // //cout<<"large pT = "<<largepT<<endl;
      // Float_t deltaR = 0;
      // Float_t small_deltaR = 100.;
      // Int_t small_position = 0;
      // // find the jet thats matched with the trigger object in delta R. 
      // for(int j = 0;j<nrefe_1;++j){

      // 	deltaR = Calc_deltaR(eta_1[j],phi_1[j], trgObj_eta_1, trgObj_phi_1);
	
      // 	if(deltaR < small_deltaR){
	  
      // 	  small_deltaR = deltaR;
      // 	  small_position = j;
	  
      // 	}

      // }

      // Float_t deltaR_leading = Calc_deltaR(eta_1[0],phi_1[0], trgObj_eta_1, trgObj_phi_1);

      //if(trgObj_pt_1 != -99)cout<<" JetMB = "<<jetMB_1<<" Jet80 = "<<jet80_1<<" Jet65 = "<<jet65_1<<" Jet55 = "<<jet55_1<<endl;
      //if(trgObj_pt_1 != -99)cout<<" Trigger object pt =  "<<trgObj_pt_1<<", matched jet pt = "<<pt_1[small_position]<<", deltaR = "<<small_deltaR<<", leading jet pt = "<<pt_1[0]<<", delta R with leading jet = "<<deltaR_leading<<endl;

      // hDenominator[k]->Fill(pt_1[small_position]);
      
      // if(jet80_1==1) hNumerator_80[k]->Fill(pt_1[small_position], jet80_p_1); 
      // if(jet65_1==1) hNumerator_65[k]->Fill(pt_1[small_position], jet65_p_1); 
      // if(jet55_1==1) hNumerator_55[k]->Fill(pt_1[small_position], jet55_p_1);

      // now the jet which is matched with the trigger object is the jet at position 
      for(int j = 0;j<nrefe_1;++j){

	//if(chMax_1[j]/pt_1[j] < 0.02 || eMax_1[j]/pt_1[j]>0.6) continue; 
	
	hDenominator[k]->Fill(pt_1[j]);
	
	if(jet80_1==1) hNumerator_80[k]->Fill(pt_1[j]); 
	if(jet65_1==1) hNumerator_65[k]->Fill(pt_1[j]); 
	if(jet55_1==1) hNumerator_55[k]->Fill(pt_1[j]);
	
	//if(jet55_1==1 && jet65_1==0 && jet80_1==0) hNumerator_55_n65_n80[k]->Fill(pt_1[j]);
	//if(jet65_1==1 && jet80_1==0) hNumerator_65_n80[k]->Fill(pt_1[j]);
	

	// if(jetMB_1==1 && jet80_p_1==1) {
	
	//   hDenominator[k]->Fill(pt_1[j]);
	  
	//   if(jet80_1==1) hNumerator_80[k]->Fill(pt_1[j]); 
	//   // if(jet65_1==1) hNumerator_65[k]->Fill(pt_1[0]); 
	//   // if(jet55_1==1) hNumerator_55[k]->Fill(pt_1[0]);
	  
	//   // if(jet55_1==1 && jet65_1==0 && jet80_1==0) hNumerator_55_n65_n80[k]->Fill(pt_1[0]);
	//   // if(jet65_1==1 && jet80_1==0) hNumerator_65_n80[k]->Fill(pt_1[0]);
	  
	// }
	
	
      }
	
      //trigger_info->Fill();


      
      //cout<<"after eta selection large pT = "<<largepT<<endl;
      
      //if(largepT!=test)cout<<"large pT is different in |eta|<2"<<endl;

      // hJet80_Prescl->Fill(jet80_p_1);
      // hJet65_Prescl->Fill(jet65_p_1);
      // hJet55_Prescl->Fill(jet55_p_1);

      // //if(l1MB_1 && jetMB_1) hJetMB->Fill(pt_1[0]);
      // if(jetMB_1==1)hJetMB_80_1->Fill(largepT);
      // if(jetMB_1==1 && jet65_p_1 == 1)hJetMB_65_1->Fill(largepT);
      // if(jetMB_1==1 && jet55_p_1 == 1)hJetMB_55_1->Fill(largepT);

      // if(jetMB_1==1 && jet80_1==1) hJet80_JetMB->Fill(largepT);
      // if(jetMB_1==1 && jet65_1==1 && jet65_p_1 ==1 ) hJet65_JetMB->Fill(largepT);
      // if(jetMB_1==1 && jet55_1==1 && jet55_p_1 ==1 ) hJet55_JetMB->Fill(largepT);
      
      //if(printDebug && jentry%100==0) cout<<"MB prescl = "<<jetMB_p_1<<endl;
      //if(printDebug && jentry%100==0) cout<<"Jet 55 prescl = "<<jet55_p_1<<endl;
      //if(printDebug && jentry%100==0) cout<<"Jet 65 prescl = "<<jet65_p_1<<endl;
      //if(printDebug && jentry%100==0) cout<<"Jet 80 prescl = "<<jet80_p_1<<endl;

    }//nentries_jet55or65 loop
    
  }//radius loop.
  
  // TH1F * hJet55_TrigTurnon = (TH1F*)hJet55_JetMB->Clone("hJet55_TrigTurnon");
  // hJet55_TrigTurnon->Divide(hJetMB);

  // TH1F * hJet65_TrigTurnon = (TH1F*)hJet65_JetMB->Clone("hJet65_TrigTurnon");
  // hJet65_TrigTurnon->Divide(hJetMB);

  // TH1F * hJet80_TrigTurnon = (TH1F*)hJet80_JetMB->Clone("hJet80_TrigTurnon");
  // hJet80_TrigTurnon->Divide(hJetMB);

  // TH1F * hL1SJ36_TrigTurnon = (TH1F*)hL1SJ36_JetMB->Clone("hL1SJ36_TrigTurnon");
  // hL1SJ36_TrigTurnon->Divide(hJetMB);

  // TH1F * hL1SJ52_TrigTurnon = (TH1F*)hL1SJ52_JetMB->Clone("hL1SJ52_TrigTurnon");
  // hL1SJ52_TrigTurnon->Divide(hJetMB); 
  
  if(Type=="Data"){
    TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_MinBiasUPC_ak%s%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
    f.cd();

    //jets_ID->Write();
    //hJets->Write();
    //hEvents->Write();

    //trigger_info->Write();
    for(int k = 0;k<no_radius;k++){
      hDenominator[k]->Write();
      hNumerator_80[k]->Write();
      hNumerator_65[k]->Write();
      hNumerator_55[k]->Write();
      hNumerator_55_n65_n80[k]->Write();
      hNumerator_65_n80[k]->Write();
    }
    
    
#if 0
    hJetMB_55_1->Write();
    hJetMB_55_1->Print("base");
    hJetMB_65_1->Write();
    hJetMB_65_1->Print("base");
    hJetMB_80_1->Write();
    hJetMB_80_1->Print("base");

    // for(int k = 0;k<no_radius;++k){
    //   for(int i = 0;i<nbins_cent;++i){
    // 	hJetMBSpectra[k][i]->Write();
    // 	hJetMBSpectra[k][i]->Print("base");
    //   }
    // }
    
    //hJet80->Write();
    //hJet65->Write();
    //hJet55->Write();
    //hL1MB->Write();
    //hL1SJ36->Write();
    //hL1SJ52->Write();
    hJet80_JetMB->Write();
    hJet80_JetMB->Print("base");
    hJet65_JetMB->Write();
    hJet65_JetMB->Print("base");
    hJet55_JetMB->Write();
    hJet55_JetMB->Print("base");
    hJet80_Prescl->Write();
    hJet80_Prescl->Print("base");
    hJet65_Prescl->Write();
    hJet65_Prescl->Print("base");
    hJet55_Prescl->Write();
    hJet55_Prescl->Print("base");
    // hJet65_Jet80->Write();
    // hJet65_Jet80->Print("base");
    // hJet55_Jet80->Write();
    // hJet55_Jet80->Print("base");
    // hJet80_Trig->Write();
    // hJet80_Trig->Print("base");
    //hL1MB_JetMB->Write();
    //hL1SJ36_JetMB->Write();
    //hL1SJ52_JetMB->Write();
    // hJet55_TrigTurnon->Write();
    // hJet65_TrigTurnon->Write();
    // hJet80_TrigTurnon->Write();
    // hL1SJ36_TrigTurnon->Write();
    // hL1SJ52_TrigTurnon->Write();
#endif
  }

  if(Type=="MC"){
    TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_HydjetMinBias_jetntuple_withEvtCuts_ak%s3%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
    f.cd();

    //jets_ID->Write();
    //hJets->Write();
    //hEvents->Write();
  }
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;

}
