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

void RAA_read_MinBias(int startfile = 0, int endfile = 1, char *algo = "Pu", char *jet_type = "Calo", char *Type = "Data"){

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
  
  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;
  
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }
  
  const int N = 4;
  
  TChain *jetpbpb1[N][no_radius];

  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%d%sJetAnalyzer",algo,list_radius[k],jet_type);
    //dir[3][k] = Form("akPu%d%sJetAnalyzer",list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    //dir[4][k] = "hltobject"; //- not there in MinBias files 
    //dir[6][k] = "pfcandAnalyzer";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    //"t",
    "HiTree",
    // "jetObjTree",
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
    //jetpbpb1[2][k]->AddFriend(jetpbpb1[4][k]);
  }// radius loop ends

  cout<<"total no of entries in the combined forest files = "<<jetpbpb1[2][0]->GetEntries()<<endl;

  #if 0
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
    jetpbpb1[2][k]->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1",&l1MB_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl",&l1MB_p_1);
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

#endif

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

  TH1F *hJetMB = new TH1F("hJetMB","",200,0,200);
  TH1F *hJetMBSpectra[no_radius][nbins_cent];
  for(int k = 0;k<no_radius;++k){
    for(int i = 0;i<nbins_cent;++i)
      hJetMBSpectra[k][i] = new TH1F(Form("hJetMBSpectra_R%d_cent%d",list_radius[k],i),"Data from MB trigger alone to add to the Jet triggered data",1000,0,1000);
  }
  
  TH1F *hL1MB = new TH1F("hL1MB","",200,0,200);
  TH1F *hJet80 = new TH1F("hJet80","",200,0,200);
  TH1F *hJet65 = new TH1F("hJet65","",200,0,200);
  TH1F *hJet55 = new TH1F("hJet55","",200,0,200);
  TH1F *hL1SJ36 = new TH1F("hL1SJ36","",200,0,200);
  TH1F *hL1SJ52 = new TH1F("hL1SJ52","",200,0,200);
  TH1F *hL1MB_JetMB = new TH1F("hL1MB_JetMB","",200,0,200);
  TH1F *hJet80_JetMB = new TH1F("hJet80_JetMB","",200,0,200);
  TH1F *hJet65_JetMB = new TH1F("hJet65_JetMB","",200,0,200);
  TH1F *hJet55_JetMB = new TH1F("hJet55_JetMB","",200,0,200);
  TH1F *hL1SJ36_JetMB = new TH1F("hL1SJ36_JetMB","",200,0,200);
  TH1F *hL1SJ52_JetMB = new TH1F("hL1SJ52_JetMB","",200,0,200);
  
  for(int k = 0;k<no_radius;++k){
    
    jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet80_JetMB","1*(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&HLT_HIJet80_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&L1_SingleJet52_BptxAND)","goff");
    jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet65_JetMB","1.11287*(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&HLT_HIJet65_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&L1_SingleJet36_BptxAND)","goff");
    jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJet55_JetMB","2.0292*(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&HLT_HIJet55_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter&&L1_SingleJet36_BptxAND)","goff");
    jetpbpb1[2][0]->Draw("Max$(jtpt)>>hJetMB","(pcollisionEventSelection&&HLT_HIMinBiasHfOrBSC_v1&&abs(jteta)<2&&abs(vz)<15&&pHBHENoiseFilter)","goff");
    
    // for(int i = 0;i<nbins_cent;i++){
    //   //cout<<Form("%f %f",5*boundaries_cent[i],5*boundaries_cent[i+1])<<endl;
    //   //cout<<Form("38.695*(pcollisionEventSelection && HLT_HIMinBiasHfOrBSC_v1 && abs(jteta)<2 && abs(vz)<15 && pHBHENoiseFilter && !HLT_HIJet80_v1 && !HLT_HIJet65_v1 && !HLT_HIJet55_v1 && hiBin>=%f && hiBin<%f)",5*boundaries_cent[i],5*boundaries_cent[i+1])<<endl;
    //   jetpbpb1[2][k]->Draw(Form("jtpt>>hJetMBSpectra_R%d_cent%d",list_radius[k],i),Form("38.695*(pcollisionEventSelection && HLT_HIMinBiasHfOrBSC_v1 && abs(jteta)<2 && abs(vz)<15 && pHBHENoiseFilter && !HLT_HIJet80_v1 && !HLT_HIJet65_v1 && !HLT_HIJet55_v1 && hiBin>=%f && hiBin<%f && (chargedMax/jtpt > 0.02 || eMax/jtpt<0.6))",5*boundaries_cent[i],5*boundaries_cent[i+1]),"goff");
    //   //hJetMBSpectra[i]->Print("base");
    // }
    
  }
  
  // prescl is calculated as the ratio of the 
  //minbias prescl = 38.695
  // jet80  prescl = 1 (same in minbias sample as well)
  // jet65  prescl = 1.11398 (1.11287)
  // jet55  prescl = 2.0475 (2.0292)

  // but these prescl values are derived from the triggered data sample, need to check if they are the same in the minbias sample 

#if 0
  for(int k = 0;k<no_radius;k++){

    if(printDebug)cout<<"Running data reading for R = "<<list_radius[k]<<endl;
    // loop for the jetpbpb1[2] tree 
    Long64_t nentries_MB = jetpbpb1[2][k]->GetEntries();
    if(printDebug)cout<<"nentries_MB = "<<nentries_MB<<endl;

    cout<<"No of events with HLT_HIMinBiasHfOrBSC_v1 = "<<jetpbpb1[2][k]->GetEntries("HLT_HIMinBiasHfOrBSC_v1")<<endl;
    cout<<"No of events with HLT_HIMinBiasHfOrBSC_v1 = "<<jetpbpb1[2][k]->GetEntries("HLT_HIMinBiasHfOrBSC_v1")<<endl;
    cout<<"No of events with L1_ZeroBias = "<<jetpbpb1[2][k]->GetEntries("L1_ZeroBias")<<endl;
    cout<<"No of events with HLT_HIJet55_v1 = "<<jetpbpb1[2][k]->GetEntries("HLT_HIJet55_v1")<<endl;
    cout<<"No of events with HLT_HIJet65_v1 = "<<jetpbpb1[2][k]->GetEntries("HLT_HIJet65_v1")<<endl;
    cout<<"No of events with HLT_HIJet80_v1 = "<<jetpbpb1[2][k]->GetEntries("HLT_HIJet80_v1")<<endl;
    cout<<"No of events with L1_SingleJet36_BptxAND = "<<jetpbpb1[2][k]->GetEntries("L1_SingleJet36_BptxAND")<<endl;
    cout<<"No of events with L1_SingleJet52_BptxAND = "<<jetpbpb1[2][k]->GetEntries("L1_SingleJet52_BptxAND")<<endl;
    
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

      hEvents->Fill(1);
      
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

      
      if(eta_1[0]<-2 || eta_1[0]>=2) continue;
      if(pt_1[0] < 30) continue;
      /*
	if(jetMB_1) hJetMB->Fill(pt_1[g],jetMB_p_1);
	if(L1_MB_1) hL1MB->Fill(pt_1[g],L1_MB_p_1);
	if(L1_sj36_1) hL1SJ36->Fill(pt_1[g],L1_sj36_1);
	if(L1_sj52_1) hL1SJ52->Fill(pt_1[g],L1_sj52_1);
	if(jet80_1) hJet80->Fill(pt_1[g],jet80_p_1);
	if(jet65_1) hJet65->Fill(pt_1[g],jet65_p_1);
	if(jet55_1) hJet55->Fill(pt_1[g],jet55_p_1);
	if(jetMB_1 && L1_MB_1) hL1MB_JetMB->Fill(pt_1[g],L1_MB_p_1*jetMB_p_1);
	if(jetMB_1 && L1_sj36_1) hL1SJ36_JetMB->Fill(pt_1[g],L1_sj36_1*jetMB_p_1);
	if(jetMB_1 && L1_sj52_1) hL1SJ52_JetMB->Fill(pt_1[g],L1_sj52_1*jetMB_p_1);
	if(jetMB_1 && jet80_1) hJet80_JetMB->Fill(pt_1[g],jet80_p_1*jetMB_p_1);
	if(jetMB_1 && jet65_1) hJet65_JetMB->Fill(pt_1[g],jet65_p_1*jetMB_p_1);
	if(jetMB_1 && jet55_1) hJet55_JetMB->Fill(pt_1[g],jet55_p_1*jetMB_p_1);	  
      */
	  
      //if(l1MB_1 && jetMB_1) hJetMB->Fill(pt_1[0]);
      hJetMB->Fill(pt_1[0]);

      if(jetMB_1 && !jet80_1 && !jet65_1 && !jet55_1) hJetMBSpectra->Fill(pt_1[0],jetMB_p_1);

      if(L1_MB_1) hL1MB->Fill(pt_1[0]);

      if(L1_sj36_1) hL1SJ36->Fill(pt_1[0]);

      if(L1_sj52_1) hL1SJ52->Fill(pt_1[0]);

      if(jet80_1) hJet80->Fill(pt_1[0]);

      if(jet65_1) hJet65->Fill(pt_1[0]);

      if(jet55_1) hJet55->Fill(pt_1[0]);

      if(jetMB_1 && L1_MB_1) hL1MB_JetMB->Fill(pt_1[0]);

      if(jetMB_1 && L1_sj36_1) hL1SJ36_JetMB->Fill(pt_1[0]);

      if(jetMB_1 && L1_sj52_1) hL1SJ52_JetMB->Fill(pt_1[0]);

      //if(l1MB_1 && jetMB_1 && jet80_1) hJet80_JetMB->Fill(pt_1[0]);
      if(jet80_1) hJet80_JetMB->Fill(pt_1[0]);

      //if(l1MB_1 && jetMB_1 && jet65_1) hJet65_JetMB->Fill(pt_1[0]);
      if(jet65_1) hJet65_JetMB->Fill(pt_1[0]);

      //if(l1MB_1 && jetMB_1 && jet55_1) hJet55_JetMB->Fill(pt_1[0]);
      if(jet55_1) hJet55_JetMB->Fill(pt_1[0]);

      
      //if(printDebug && jentry%100==0) cout<<"MB prescl = "<<jetMB_p_1<<endl;
      //if(printDebug && jentry%100==0) cout<<"Jet 55 prescl = "<<jet55_p_1<<endl;
      //if(printDebug && jentry%100==0) cout<<"Jet 65 prescl = "<<jet65_p_1<<endl;
      //if(printDebug && jentry%100==0) cout<<"Jet 80 prescl = "<<jet80_p_1<<endl;

#if 0
      for(int j = 0;j<nbins_eta;j++){
	
	for(int g = 0;g<nrefe_1;g++){ // this is the loop for the  Jets we are interested in.  
	  
	  if(eta_1[g]<boundaries_eta[j][0] || eta_1[g]>=boundaries_eta[j][1]) continue;
	  
	  /*
	  if(jetMB_1) hJetMB->Fill(pt_1[g],jetMB_p_1);
	  if(L1_MB_1) hL1MB->Fill(pt_1[g],L1_MB_p_1);
	  if(L1_sj36_1) hL1SJ36->Fill(pt_1[g],L1_sj36_1);
	  if(L1_sj52_1) hL1SJ52->Fill(pt_1[g],L1_sj52_1);
	  if(jet80_1) hJet80->Fill(pt_1[g],jet80_p_1);
	  if(jet65_1) hJet65->Fill(pt_1[g],jet65_p_1);
	  if(jet55_1) hJet55->Fill(pt_1[g],jet55_p_1);
	  if(jetMB_1 && L1_MB_1) hL1MB_JetMB->Fill(pt_1[g],L1_MB_p_1*jetMB_p_1);
	  if(jetMB_1 && L1_sj36_1) hL1SJ36_JetMB->Fill(pt_1[g],L1_sj36_1*jetMB_p_1);
	  if(jetMB_1 && L1_sj52_1) hL1SJ52_JetMB->Fill(pt_1[g],L1_sj52_1*jetMB_p_1);
	  if(jetMB_1 && jet80_1) hJet80_JetMB->Fill(pt_1[g],jet80_p_1*jetMB_p_1);
	  if(jetMB_1 && jet65_1) hJet65_JetMB->Fill(pt_1[g],jet65_p_1*jetMB_p_1);
	  if(jetMB_1 && jet55_1) hJet55_JetMB->Fill(pt_1[g],jet55_p_1*jetMB_p_1);	  
	  */
	  
	  
	  if(jetMB_1) hJetMB->Fill(pt_1[g]);

	  if(jetMB_1 && !jet80_1 && !jet65_1 && !jet55_1) hJetMBSpectra->Fill(pt_1[g],jetMB_p_1);

	  if(L1_MB_1) hL1MB->Fill(pt_1[g]);

	  if(L1_sj36_1) hL1SJ36->Fill(pt_1[g]);

	  if(L1_sj52_1) hL1SJ52->Fill(pt_1[g]);

	  if(jet80_1) hJet80->Fill(pt_1[g]);

	  if(jet65_1) hJet65->Fill(pt_1[g]);

	  if(jet55_1) hJet55->Fill(pt_1[g]);

	  if(jetMB_1 && L1_MB_1) hL1MB_JetMB->Fill(pt_1[g]);

	  if(jetMB_1 && L1_sj36_1) hL1SJ36_JetMB->Fill(pt_1[g]);

	  if(jetMB_1 && L1_sj52_1) hL1SJ52_JetMB->Fill(pt_1[g]);

	  if(jetMB_1 && jet80_1) hJet80_JetMB->Fill(pt_1[g]);

	  if(jetMB_1 && jet65_1) hJet65_JetMB->Fill(pt_1[g]);

	  if(jetMB_1 && jet55_1) hJet55_JetMB->Fill(pt_1[g]);

	  
	  /*
	  if(Type=="Data"){
	    arrayValues_Data[0] = raw_1[g];
	    arrayValues_Data[1] = pt_1[g];
	    arrayValues_Data[2] = jtpu_1[g];
	    arrayValues_Data[3] = L1_MB_1; 
	    arrayValues_Data[4] = L1_MB_p_1; 
	    arrayValues_Data[5] = jetMB_1; 
	    arrayValues_Data[6] = jetMB_p_1; 
	    arrayValues_Data[7] = centBin;
	    arrayValues_Data[8] = chMax_1[g];
	    arrayValues_Data[9] = chSum_1[g];
	    arrayValues_Data[10] = phMax_1[g];
	    arrayValues_Data[11] = phSum_1[g];
	    arrayValues_Data[12] = neMax_1[g];
	    arrayValues_Data[13] = neSum_1[g];
	    arrayValues_Data[14] = muMax_1[g];
	    arrayValues_Data[15] = muSum_1[g];
	    arrayValues_Data[16] = eMax_1[g];
	    arrayValues_Data[17] = eSum_1[g];

	    jets_ID->Fill(arrayValues_Data);
	  }

	  if(Type=="MC"){
	    arrayValues_MC[0] = raw_1[g];
	    arrayValues_MC[1] = pt_1[g];
	    arrayValues_MC[2] = jtpu_1[g];
	    arrayValues_MC[3] = L1_MB_1; 
	    arrayValues_MC[4] = L1_MB_p_1; 
	    arrayValues_MC[5] = jetMB_1; 
	    arrayValues_MC[6] = jetMB_p_1; 
	    arrayValues_MC[7] = centBin;
	    arrayValues_MC[8] = chMax_1[g];
	    arrayValues_MC[9] = chSum_1[g];
	    arrayValues_MC[10] = phMax_1[g];
	    arrayValues_MC[11] = phSum_1[g];
	    arrayValues_MC[12] = neMax_1[g];
	    arrayValues_MC[13] = neSum_1[g];
	    arrayValues_MC[14] = muMax_1[g];
	    arrayValues_MC[15] = muSum_1[g];
	    arrayValues_MC[16] = eMax_1[g];
	    arrayValues_MC[17] = eSum_1[g];
	    arrayValues_MC[18] = refpt_1[g];
	    arrayValues_MC[19] = subid_1[g];
	   
	    jets_ID->Fill(arrayValues_MC);

	  }
	  */

	  hJets->Fill(1);
	  
	    
	}//jet loop
	  
      }//eta bin loop

#endif
      
    }//nentries_jet55or65 loop
    
  }//radius loop.

#endif
  
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
    TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_MinBiasUPC_trigger_turnoncurves_SuperNovaRejected_ak%s%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
    f.cd();

    //jets_ID->Write();
    //hJets->Write();
    //hEvents->Write();
    
    hJetMB->Write();
    hJetMB->Print("base");
    for(int k = 0;k<no_radius;++k){
      for(int i = 0;i<nbins_cent;++i){
	hJetMBSpectra[k][i]->Write();
	hJetMBSpectra[k][i]->Print("base");
      }
    }
    
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
    //hL1MB_JetMB->Write();
    //hL1SJ36_JetMB->Write();
    //hL1SJ52_JetMB->Write();
    // hJet55_TrigTurnon->Write();
    // hJet65_TrigTurnon->Write();
    // hJet80_TrigTurnon->Write();
    // hL1SJ36_TrigTurnon->Write();
    // hL1SJ52_TrigTurnon->Write();
    
  }

  if(Type=="MC"){
    TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_HydjetMinBias_jetntuple_withEvtCuts_SuperNovaRejected_ak%s3%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
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
