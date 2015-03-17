// Raghav Kunnawalkam Elayavalli
// June 5th 2014
// CERN
// for questions or comments: raghav.k.e at CERN dot CH

// 
// read all the MC files for PbPb and pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 
// 

// July 19 - all pp histograms will have 2D arrays with [radius][eta_bin]. the PbPb histograms will be defined by 3D arrays with [radius][eta_bin][centrality]. 
// July 20 - the loop structure(s) are defined as follows for the several histograms 
//            Radius    : iteration variable: k;       number of iterations: no_radius;                                 values for the radii: list_radius
//            Eta bins  : iteration variable: j;       number of iterations: nbins_eta;                                 values of the bins  : boundaries_eta
//            Centrality: iteration variable: i;       number of iterations: nbins_cent +1;                             values of the bins  : boundaries_cent (till nbins_cent) + 0-200 (the whole range) for the final iteration. 
//            p_T Hats  : iteration variable: h;       number of iterations: nbins_pthat (PbPb) and nbinsPP_pthat (pp); values of the bins  : boundaries_pthat (PbPb) and boundariesPP_pthat (pp)  
//            jets      : iteration variable: g;       number of iterations: no of jets in Data[k][h];
//            p_T       : defined just below as nbins_pt with 39 bins. to match our NLO and jet RpA analysis bins. 

// Oct 23 - removed the cuts from the MC -> like the noisefilter etc... 

// Nov 4th - added the supernova event cut rejection based on the no of hits in the pixel. 

// Dec 9th - going to PU for the Jet RAA. 

// Dec 17th - changing the file list to smaller 50k files on which JEC were derived to check for PF electron problems, requested by Marguerite.

// Jan 13th 2015 - adding in the official pp mc (from Dragos) 
//               - this is going to be a bit tricky since each file is split up into 4 smaller files. so each pthat will have a TChain!


// Feb 12th - cleaned up the macro to make it usable (hopefuly) by others. 

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
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
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"


//static const int nbins_pt = 29; //old bins with slight difference in the low and high pt ranges. 
//static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

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

static const char etaWidth[nbins_eta][256] = {
  "n20_eta_p20"
};


static const int no_radius = 1;//testing purposes 
static const int list_radius[no_radius] = {3};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 3; 
//static const int list_radius[no_radius] = {2,3,4};

static const int nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 
static const Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
static const int trigValue = 4;
static const char trigName [trigValue][256] = {"HLT55","HLT65","HLT80","Combined"};
static const Float_t effecPrescl = 2.047507;

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


// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

class JetData
{
public:
  JetData(char *fileName, char *jetTree, char *genJetTree, bool loadGenJet = 0,bool isPbPb = 0) {
    cout <<"Open "<<fileName<<endl;
    tFile = new TFile(fileName,"read");
    tEvt = (TTree*)tFile->Get("hiEvtAnalyzer/HiTree");
    tSkim = (TTree*)tFile->Get("skimanalysis/HltTree");
    tHlt = (TTree*)tFile->Get("hltanalysis/HltTree");
    tpfCand = (TTree*)tFile->Get("pfcandAnalyzer/pfTree");
    tJet = (TTree*)tFile->Get(jetTree);
    tJet->SetBranchAddress("jtpt" , jtpt );
    if(isPbPb)tJet->SetBranchAddress("jtpu", jtpu );
    if(isPbPb)tEvt->SetBranchAddress("hiNpix",&hiNpix);
    if(isPbPb)tEvt->SetBranchAddress("run",&run);
    if(isPbPb)tEvt->SetBranchAddress("evt",&evt);
    if(isPbPb)tEvt->SetBranchAddress("lumi",&lumi);
    if(isPbPb)tEvt->SetBranchAddress("hiNtracks",&hiNtracks);
    if(isPbPb)tEvt->SetBranchAddress("hiHF",&hiHF);

    tpfCand->SetBranchAddress("nPFpart",&nPFpart);
    tpfCand->SetBranchAddress("pfId",&pfID);
    tpfCand->SetBranchAddress("pfPt",&pfPt);
    tpfCand->SetBranchAddress("pfVsPt",&pfVsPt);
    tpfCand->SetBranchAddress("pfVsPtInitial",&pfVsPtInitial);
    tpfCand->SetBranchAddress("pfArea",&pfArea);
    tpfCand->SetBranchAddress("pfEta",&pfEta);
    tpfCand->SetBranchAddress("pfPhi",&pfPhi);
    tpfCand->SetBranchAddress("vn",&v_n);
    tpfCand->SetBranchAddress("psin",&psi_n);
    tpfCand->SetBranchAddress("sumpt",&sumpT);
    
    tJet->SetBranchAddress("rawpt", rawpt);
    tJet->SetBranchAddress("trackMax" , trackMax );
    tJet->SetBranchAddress("chargedMax",chargedMax);
    tJet->SetBranchAddress("chargedSum",chargedSum);
    tJet->SetBranchAddress("neutralMax",neutralMax);
    tJet->SetBranchAddress("neutralSum",neutralSum);
    tJet->SetBranchAddress("photonSum",photonSum);
    tJet->SetBranchAddress("photonMax",photonMax);
    tJet->SetBranchAddress("eSum",eSum);
    tJet->SetBranchAddress("eMax",eMax);
    tJet->SetBranchAddress("muSum",muSum);
    tJet->SetBranchAddress("muMax",muMax);
    tJet->SetBranchAddress("refpt", refpt);
    tJet->SetBranchAddress("nref" ,&njets);
    tJet->SetBranchAddress("jteta", jteta);
    tJet->SetBranchAddress("jtphi", jtphi);
    tJet->SetBranchAddress("jtm", jtmass);
    tJet->SetBranchAddress("pthat", &pthat);
    if(isPbPb) tJet->SetBranchAddress("subid",&subid);
    if (loadGenJet) tGenJet = (TTree*)tFile->Get(genJetTree);
    if (loadGenJet) tGenJet->SetBranchAddress("ngen" ,&ngen);
    if (loadGenJet) tGenJet->SetBranchAddress("genpt", genpt);
    if (loadGenJet) tGenJet->SetBranchAddress("gensubid", gensubid);
    tEvt->SetBranchAddress("hiBin",&bin);
    tEvt->SetBranchAddress("vz",&vz);
    tEvt->SetBranchAddress("vx",&vx);
    tEvt->SetBranchAddress("vy",&vy);
    tSkim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    if(isPbPb) tSkim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    else tSkim->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
    if(isPbPb) tHlt->SetBranchAddress("HLT_HIJet80_v7",&jet80_1);
    if(isPbPb) tHlt->SetBranchAddress("HLT_HIJet65_v7",&jet65_1);
    if(isPbPb) tHlt->SetBranchAddress("HLT_HIJet55_v7",&jet55_1);
    if(isPbPb) tHlt->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_p_1);
    if(isPbPb) tHlt->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_p_1);
    if(isPbPb) tHlt->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_p_1);
    tJet->AddFriend(tEvt);
    tJet->AddFriend(tSkim);
    tJet->AddFriend(tHlt);
    tJet->AddFriend(tpfCand);

  };
  
  TFile *tFile;
  TTree *tJet;
  TTree *tGenJet;
  TTree *tEvt;
  TTree *tHlt;
  TTree *tSkim;
  TTree *tpfCand;
 
  //event varianbles
  int run;
  int evt;
  int lumi;
  int hiNTracks;
  Float_t hiHF;
  float vz;
  float vx;
  float vy;
  float pthat;
  int hiNpix;
  int hiNtracks;
  int hiNevtPlane;
  Float_t hiEvtPlanes[1000];


  //jet variables 
  float jtpt[1000];
  float jtpu[1000];
  float rawpt[1000];
  float refpt[1000];
  float jteta[1000];
  float jtphi[1000];
  float jtmass[1000];
  float trackMax[1000];
  float chargedMax[1000];
  float neutralMax[1000];
  float chargedSum[1000];
  float neutralSum[1000];
  float photonSum[1000];
  float eSum[1000];
  float muSum[1000];
  float photonMax[1000];
  float eMax[1000];
  float muMax[1000];
  float genpt[1000];
  int gensubid[1000];
  float subid[1000];


  int njets;
  int ngen;
  int bin;     
  int pHBHENoiseFilter;
  int pPAcollisionEventSelectionPA;
  int pcollisionEventSelection;
  int jet55_1;
  int jet65_1;
  int jet80_1;
  int jet55_p_1;
  int jet65_p_1;
  int jet80_p_1;

  // //from the pfcandanalyzer tree:
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


};


using namespace std;

void RAA_read_mc(int startfile = 0, int endfile = 9, char *algo = "Pu", char *jet_type = "PF", int sub_id = 0){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  cout<<"Running for Algorithm "<<algo<<" "<<jet_type<<endl;
 
  bool printDebug = true;
  TDatime date;

  const int nbins_pthat = 9;
  Double_t boundaries_pthat[nbins_pthat+1];
  char *fileName_pthat[nbins_pthat+1];
  Double_t xsection[nbins_pthat+1];
  //Double_t entries[nbins_pthat]; 
  //there are two ways in which we can select the no of events we use to scale - it has to be between the pthat range. 
  //first file name - partial 50K statistics, second one is full statistics sample. 
  //similarly the entries number is for the small statistics. 
  
  // refer this twiki for the data and MC files: http://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForestPA2014#PYTHIA_HYDJET_embedded_sample

  boundaries_pthat[0]=15;
  fileName_pthat[0] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat15_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[0] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT15_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[0] = "/export/d00/scratch/dav2105/badjets/bad15.root";
  xsection[0]= 2.034e-01;
  //entries[0] = ;//total - 48588
  
  boundaries_pthat[1]=30;
  fileName_pthat[1] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat30_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[1] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT30_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[1] = "/export/d00/scratch/dav2105/badjets/bad30.root";
  xsection[1]= 1.075e-02;
  // entries[1] = ;//total - 48428
  
  boundaries_pthat[2]=50;
  fileName_pthat[2] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat50_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[2] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT50_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[2] = "/export/d00/scratch/dav2105/badjets/bad50.root";
  xsection[2]= 1.025e-03;
  // entries[2] = ;//total - 50000
  
  boundaries_pthat[3]=80;
  fileName_pthat[3] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[3] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT80_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[3] = "/export/d00/scratch/dav2105/badjets/bad80.root";
  xsection[3]= 9.865e-05;
  // entries[3] = ;//total - 49500
  
  boundaries_pthat[4]=120;
  fileName_pthat[4] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat120_Track9_Jet30_matchEqR_merged_forest_0.root";  
  //fileName_pthat[4] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT120_Hydjet1p8_STARTHI53_LV1_v15_full.root";  
  //fileName_pthat[4] = "/export/d00/scratch/dav2105/badjets/bad120.root";
  xsection[4]= 1.129e-05;
  // entries[4] = ;//total - 49500

  boundaries_pthat[5]=170;
  fileName_pthat[5] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat170_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[5] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT170_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[5] = "/export/d00/scratch/dav2105/badjets/bad120.root";
  xsection[5]= 1.465e-06;
  // entries[5] = ;//total - 49444

  boundaries_pthat[6]=220;
  fileName_pthat[6] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[6] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT220_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[6] = "/export/d00/scratch/dav2105/badjets/bad220.root";
  xsection[6]= 2.837e-07;
  // entries[6] = ;//total - 49460

  boundaries_pthat[7]=280;
  fileName_pthat[7] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat280_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[7] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT280_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[7] = "/export/d00/scratch/dav2105/badjets/bad280.root";
  xsection[7]= 5.323e-08;
  // entries[7] = ;//total - 49541

  boundaries_pthat[8]=370;
  fileName_pthat[8] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat370_Track9_Jet30_matchEqR_merged_forest_0.root";
  //fileName_pthat[8] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT370_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  //fileName_pthat[8] = "/export/d00/scratch/dav2105/badjets/bad370.root";
  xsection[8]= 5.934e-09;
  // entries[8] = ;//total - 19031

  boundaries_pthat[9] = 2000;
  xsection[9] = 0.0;

  // Vertex & centrality reweighting for PbPb
  TF1 *fVz;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);

  //get the centrality weight from the root file created in the plotting macro. 
  TFile *fcentin = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_DataMC_cent_ratio_20141117.root");
  TH1F *hCentWeight = (TH1F*)fcentin->Get("hCentRatio");

  // lets declare all the histograms here. 

  TH1F *hpbpb_gen[no_radius][nbins_eta][nbins_cent+1],*hpbpb_reco[no_radius][nbins_eta][nbins_cent+1];
  TH2F *hpbpb_matrix[no_radius][nbins_eta][nbins_cent+1];
  TH2F *hpbpb_matrix_HLT[no_radius][nbins_eta][nbins_cent+1];
  TH2F *hpbpb_mcclosure_matrix[no_radius][nbins_eta][nbins_cent+1];
  TH2F *hpbpb_mcclosure_matrix_HLT[no_radius][nbins_eta][nbins_cent+1];
  //TH2F *hpbpb_response[nbins_cent+1];
  TH1F *hpbpb_mcclosure_JetComb_data[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_data[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_Jet80_data[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_Jet65_data[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_Jet55_data[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_gen[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_JetComb_gen[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_Jet80_gen[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_Jet65_gen[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_mcclosure_Jet55_gen[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_jtpu[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_jtpu_noScale[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_vx[no_radius];
  TH1F *hpbpb_vy[no_radius];
  TH1F *hpbpb_vz[no_radius];
  TH1F *hpbpb_cent[no_radius];
  TH1F *hpbpb_RecoOverRaw[no_radius][nbins_eta][nbins_cent+1];
  TH2F *hpbpb_RecoOverRaw_jtpt[no_radius][nbins_eta][nbins_cent+1];
  //TH2F *hpbpb_RecoOverref_refpt[no_radius][nbins_eta][nbins_cent+1];


  //get the spectra with the specific trigger object from the different files. 
  TH1F *hpbpb_Jet80_gen[no_radius][nbins_eta][nbins_cent+1],*hpbpb_Jet80_reco[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet65_gen[no_radius][nbins_eta][nbins_cent+1],*hpbpb_Jet65_reco[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet55_gen[no_radius][nbins_eta][nbins_cent+1],*hpbpb_Jet55_reco[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_JetComb_gen[no_radius][nbins_eta][nbins_cent+1],*hpbpb_JetComb_reco[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet80_raw[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet65_raw[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet55_raw[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_JetComb_raw[no_radius][nbins_eta][nbins_cent+1];

  //TH1F *hpbpb_eta[no_radius][nbins_eta][nbins_cent+1], *hpbpb_phi[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eta_full[no_radius], *hpbpb_phi_full[no_radius];
  TH1F *hpbpb_eta_full_noScale[no_radius], *hpbpb_phi_full_noScale[no_radius];
  
  TH1F *hpbpb_chMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_chSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eSum[no_radius][nbins_eta][nbins_cent+1];

  TH1F *hpbpb_chMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_chSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  
  TH1F *hCentMC[no_radius];
  TH1F *hVzMC[no_radius];
  TH1F *hPbPb_pthat_fine[no_radius];
  TH1F *hPbPb_pthat_fine_noScale[no_radius];
  TH1F *hPtHat[no_radius];
  TH1F *hPtHatRaw[no_radius];
  
  // histograms for the supernova cut rejection 
  TH2F *hpbpb_Npix_before_cut[no_radius][nbins_cent+2];// the last cent is for ultra central events.  
  TH2F *hpbpb_Npix_after_cut[no_radius][nbins_cent+1]; 
  
  // histograms to check the contribution of the nJets>7 (with pT>50 and |eta|<2) to the nJets < 7 
  TH1F *hpbpb_pt_Njet_g7[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_pt_Njet_l7[no_radius][nbins_eta][nbins_cent+1];
  
  
  
  //declare the output file 
  //TNtuple *jets_ID = new TNtuple("jets_ID","","rawpt:refpt:jtpt:jtpu:jet55:jet55_prescl:jet65:jet65_prescl:jet80:jet80_prescl:scale:weight_vz:weight_cent:cent:subid:chMax:chSum:phMax:phSum:neMax:neSum:muMax:muSum:eMax:eSum");
  //Float_t arrayValues[25];

  for(int k = 0;k<no_radius;k++){
    //cout<<"radius = "<<list_radius[k]<<endl;
    for(int j = 0;j<nbins_eta;j++){
      //cout<<"eta bin = "<<j<<endl;
      for(int i = 0;i<nbins_cent;i++){
	//cout<<"cent bin = "<<i<<endl;

        hpbpb_gen[k][j][i] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Gen refpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	//cout<<"A"<<endl;
	hpbpb_reco[k][j][i] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Reco jtpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	//cout<<"B"<<endl;
	hpbpb_matrix[k][j][i] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Matrix refpt jtpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
	hpbpb_matrix_HLT[k][j][i] = new TH2F(Form("hpbpb_matrix_HLT_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
	hpbpb_mcclosure_matrix[k][j][i] = new TH2F(Form("hpbpb_mcclosure_matrix_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Matrix for mcclosure refpt jtpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
	hpbpb_mcclosure_matrix_HLT[k][j][i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
	//cout<<"C"<<endl;
	hpbpb_mcclosure_data[k][j][i] = new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("data for unfolding mc closure test R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_JetComb_data[k][j][i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_Jet80_data[k][j][i] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("data for unfolding mc closure test trigger 80  R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_Jet65_data[k][j][i] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("data for unfolding mc closure test trigger 65  R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_Jet55_data[k][j][i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

	hpbpb_mcclosure_gen[k][j][i] = new TH1F(Form("hpbpb_mcclosure_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("gen spectra for unfolding mc closure test R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_JetComb_gen[k][j][i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_Jet80_gen[k][j][i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet80_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_Jet65_gen[k][j][i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet65_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("gen spectra for unfolding mc closure test trigger 65 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_mcclosure_Jet55_gen[k][j][i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

	hpbpb_jtpu[k][j][i] = new TH1F(Form("hpbpb_jtpu_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("jtpu Vs algorithm R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,500);
	hpbpb_jtpu_noScale[k][j][i] = new TH1F(Form("hpbpb_jtpu_noScale_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("jtpu Vs algorithm not Scaled R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,500);

	hpbpb_pt_Njet_g7[k][j][i] = new TH1F(Form("hpbpb_pt_Njet_g7_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("jet spectra from the events with Njet(pT>50)>7 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_pt_Njet_l7[k][j][i] = new TH1F(Form("hpbpb_pt_Njet_l7_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("jet spectra from the events with Njet(pT>50)<7 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

	//hpbpb_eta[k][j][i] = new TH1F(Form("hpbpb_eta_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eta distribution R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],8*boundaries_cent[i],5*boundaries_cent[i+1]),60,-4,+4);
	//hpbpb_phi[k][j][i] = new TH1F(Form("hpbpb_phi_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phi distribution R%d in eta widths %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],8*boundaries_cent[i],5*boundaries_cent[i+1]),60,-3.2,+3.2);
	
	//hpbpb_response[h] = new TH2F(Form("hpbpb_response_cent%d",i),Form("response jtpt refpt %2.0f - %2.0f cent",5*boundaries_cent[h],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
        hpbpb_JetComb_gen[k][j][i] = new TH1F(Form("hpbpb_JetComb_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
        hpbpb_JetComb_reco[k][j][i] = new TH1F(Form("hpbpb_JetComb_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
        hpbpb_Jet80_gen[k][j][i] = new TH1F(Form("hpbpb_Jet80_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
        hpbpb_Jet80_reco[k][j][i] = new TH1F(Form("hpbpb_Jet80_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_Jet65_gen[k][j][i] = new TH1F(Form("hpbpb_Jet65_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
        hpbpb_Jet65_reco[k][j][i] = new TH1F(Form("hpbpb_Jet65_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_Jet55_gen[k][j][i] = new TH1F(Form("hpbpb_Jet55_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
        hpbpb_Jet55_reco[k][j][i] = new TH1F(Form("hpbpb_Jet55_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	
	hpbpb_RecoOverRaw[k][j][i] = new TH1F(Form("hpbpb_RecoOverRaw_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco over raw ratio R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,10);
	hpbpb_RecoOverRaw_jtpt[k][j][i] = new TH2F(Form("hpbpb_RecoOverRaw_jtpt_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco over raw ratio versus jtpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,100,0,10);
	//hpbpb_RecoOverRaw[k][j][i] = new TH2F(Form("hpbpb_RecoOverRaw_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco over raw ratio R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,10);
	hpbpb_chMax[k][j][i] = new TH1F(Form("hpbpb_chMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phMax[k][j][i] = new TH1F(Form("hpbpb_phMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neMax[k][j][i] = new TH1F(Form("hpbpb_neMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muMax[k][j][i] = new TH1F(Form("hpbpb_muMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eMax[k][j][i] = new TH1F(Form("hpbpb_eMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);

	hpbpb_chSum[k][j][i] = new TH1F(Form("hpbpb_chSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phSum[k][j][i] = new TH1F(Form("hpbpb_phSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neSum[k][j][i] = new TH1F(Form("hpbpb_neSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muSum[k][j][i] = new TH1F(Form("hpbpb_muSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eSum[k][j][i] = new TH1F(Form("hpbpb_eSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);

	hpbpb_chMax_withCut[k][j][i] = new TH1F(Form("hpbpb_chMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phMax_withCut[k][j][i] = new TH1F(Form("hpbpb_phMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neMax_withCut[k][j][i] = new TH1F(Form("hpbpb_neMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muMax_withCut[k][j][i] = new TH1F(Form("hpbpb_muMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eMax_withCut[k][j][i] = new TH1F(Form("hpbpb_eMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);

	hpbpb_chSum_withCut[k][j][i] = new TH1F(Form("hpbpb_chSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phSum_withCut[k][j][i] = new TH1F(Form("hpbpb_phSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neSum_withCut[k][j][i] = new TH1F(Form("hpbpb_neSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muSum_withCut[k][j][i] = new TH1F(Form("hpbpb_muSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eSum_withCut[k][j][i] = new TH1F(Form("hpbpb_eSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	
      }// centrality bin loop

      hpbpb_pt_Njet_g7[k][j][nbins_cent] = new TH1F(Form("hpbpb_pt_Njet_g7_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("jet spectra from the events with Njet(pT>50)>7 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_pt_Njet_l7[k][j][nbins_cent] = new TH1F(Form("hpbpb_pt_Njet_l7_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("jet spectra from the events with Njet(pT>50)<7 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);

      
      hpbpb_chMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_chMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("chMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_phMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_phMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("phMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_neMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_neMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("neMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_muMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_muMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("muMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_eMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_eMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("eMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      
      hpbpb_chSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_chSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("chSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_phSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_phSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("phSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_neSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_neSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("neSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_muSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_muSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("muSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_eSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_eSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("eSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      
      hpbpb_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Gen refpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_reco[k][j][nbins_cent] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Reco jtpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_matrix[k][j][nbins_cent] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Matrix refpt jtpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      hpbpb_matrix_HLT[k][j][nbins_cent] = new TH2F(Form("hpbpb_matrix_HLT_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Matrix refpt jtpt from HLT combination R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      hpbpb_mcclosure_matrix[k][j][nbins_cent] = new TH2F(Form("hpbpb_mcclosure_matrix_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Matrix for mcclosure refpt jtpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      hpbpb_mcclosure_matrix_HLT[k][j][nbins_cent] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Matrix for mcclosure refpt jtpt trigger spectra R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);

      hpbpb_JetComb_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_JetComb_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Gen refpt from HLT trigger combined R%d %s 0-200",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_JetComb_reco[k][j][nbins_cent] = new TH1F(Form("hpbpb_JetComb_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco jtpt from HLT trigger combined R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_Jet80_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet80_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Gen refpt from Jet80 trigger R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_Jet80_reco[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet80_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco jtpt from Jet80 trigger R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_Jet65_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet65_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_Jet65_reco[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet65_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_Jet55_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet55_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_Jet55_reco[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet55_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      
      hpbpb_mcclosure_data[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("data for unfolding mc closure test R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_JetComb_data[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("data for unfolding mc closure test trigger combined R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_Jet80_data[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("data for unfolding mc closure test trigger 80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_Jet65_data[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("data for unfolding mc closure test trigger 65 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_Jet55_data[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("data for unfolding mc closure test trigger 55 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);

      hpbpb_mcclosure_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("gen for unfolding mc closure test R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_JetComb_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_JetComb_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("gen for unfolding mc closure test trigger combined R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_Jet80_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_Jet80_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("gen for unfolding mc closure test trigger 80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_Jet65_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_Jet65_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("gen for unfolding mc closure test trigger 65 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_mcclosure_Jet55_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_Jet55_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("gen for unfolding mc closure test trigger 55 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);


      hpbpb_jtpu[k][j][nbins_cent] = new TH1F(Form("hpbpb_jtpu_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("jtpu Vs algorithm R%d %s",list_radius[k],etaWidth[j]),1000,0,500);
      hpbpb_jtpu_noScale[k][j][nbins_cent] = new TH1F(Form("hpbpb_jtpu_noScale_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("jtpu Vs algorithm not Scaled R%d %s",list_radius[k],etaWidth[j]),1000,0,500);
      //hpbpb_response[nbins_cent] = new TH2F(Form("hpbpb_response_cent%d",nbins_cent),"response jtpt refpt 0-200 cent",1000,0,1000,1000,0,1000);

      hpbpb_RecoOverRaw[k][j][nbins_cent] = new TH1F(Form("hpbpb_RecoOverRaw_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco over raw ratio R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,10);
      hpbpb_RecoOverRaw_jtpt[k][j][nbins_cent] = new TH2F(Form("hpbpb_RecoOverRaw_jtpt_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco over raw ratio versus jtpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000,100,0,10);
#if 0
      hpp_gen[k][j] = new TH1F(Form("hpp_gen_R%d_%s",list_radius[k],etaWidth[j]),Form("gen refpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_reco[k][j] = new TH1F(Form("hpp_reco_R%d_%s",list_radius[k],etaWidth[j]),Form("reco jtpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_matrix[k][j] = new TH2F(Form("hpp_matrix_R%d_%s",list_radius[k],etaWidth[j]),Form("matrix refpt jtpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      hpp_mcclosure_matrix[k][j] = new TH2F(Form("hpp_mcclosure_matrix_R%d_%s",list_radius[k],etaWidth[j]),Form("matrix mcclosure refpt jtpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      //TH2F* hpp_response = new TH2F("hpp_response","response jtpt refpt",1000,0,1000,1000,0,1000);
      hpp_mcclosure_data[k][j] = new TH1F(Form("hpp_mcclosure_data_R%d_%s",list_radius[k],etaWidth[j]),Form("data for unfolding mc closure test pp R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
#endif
    }// eta bin loop
    
    hVzMC[k] = new TH1F(Form("hVzMC_R%d",list_radius[k]),Form("PbPb MC Vz R%d",list_radius[k]),60,-15,+15);    
    hCentMC[k] = new TH1F(Form("hCentMC_R%d",list_radius[k]),"",100,0,200);
    hPtHat[k] = new TH1F(Form("hPtHat_R%d",list_radius[k]),"",nbins_pthat,boundaries_pthat);
    hPtHatRaw[k] = new TH1F(Form("hPtHatRaw_R%d",list_radius[k]),"",nbins_pthat,boundaries_pthat);
    hPbPb_pthat_fine[k] = new TH1F(Form("hPbPb_pthat_fine_R%d",list_radius[k]),Form("PbPb pthat distribution for R=0.%d",list_radius[k]),700,0,700);
    hPbPb_pthat_fine_noScale[k] = new TH1F(Form("hPbPb_pthat_fine_noScale_R%d",list_radius[k]),Form("PbPb pthat distribution (unscaled) for R=0.%d",list_radius[k]),700,0,700);
    hpbpb_eta_full[k] = new TH1F(Form("hpbpb_eta_full_R%d",list_radius[k]),Form("PbPb eta distribution for R=%d",list_radius[k]),400,-4,+4);
    hpbpb_phi_full[k] = new TH1F(Form("hpbpb_phi_full_R%d",list_radius[k]),Form("PbPb phi distribution for R=%d",list_radius[k]),400,-4,+4);
    hpbpb_eta_full_noScale[k] = new TH1F(Form("hpbpb_eta_full_noScale_R%d",list_radius[k]),Form("PbPb eta distribution noScale for R=%d",list_radius[k]),400,-4,+4);
    hpbpb_phi_full_noScale[k] = new TH1F(Form("hpbpb_phi_full_noScale_R%d",list_radius[k]),Form("PbPb phi distribution noScale for R=%d",list_radius[k]),400,-4,+4);
    hpbpb_cent[k] = new TH1F(Form("hpbpb_cent_R%d",list_radius[k]),Form("centrality distributions R%d",list_radius[k]),200,0,200);
    hpbpb_vz[k] = new TH1F(Form("hpbpb_vz_R%d",list_radius[k]),Form("vz distribution R%d",list_radius[k]),60,-15,15);
    hpbpb_vx[k] = new TH1F(Form("hpbpb_vx_R%d",list_radius[k]),Form("vx distribution R%d",list_radius[k]),60,-15,15);
    hpbpb_vy[k] = new TH1F(Form("hpbpb_vy_R%d",list_radius[k]),Form("vy distribution R%d",list_radius[k]),60,-15,15);

    
    for(int i = 0;i<nbins_cent;i++){
      hpbpb_Npix_before_cut[k][i] = new TH2F(Form("hpbpb_Npix_before_cut_R%d_n20_eta_p20_cent%d",list_radius[k],i),Form("Number of pixels hit per no of jets pT>50 before cut R%d n20_eta_p20 %2.0f - %2.0f cent",list_radius[k],5*boundaries_cent[i],5*boundaries_cent[i+1]),50,0,50,100,0,60000);
      hpbpb_Npix_after_cut[k][i] = new TH2F(Form("hpbpb_Npix_after_cut_R%d_n20_eta_p20_cent%d",list_radius[k],i),Form("Number of pixels hit per no of jets pT>50 after cut R%d n20_eta_p20 %2.0f - %2.0f cent",list_radius[k],5*boundaries_cent[i],5*boundaries_cent[i+1]),50,0,50,100,0,60000);
    }

    hpbpb_Npix_before_cut[k][nbins_cent] = new TH2F(Form("hpbpb_Npix_before_cut_R%d_n20_eta_p20_cent%d",list_radius[k],nbins_cent),Form("Number of pixels hit per no of jets pT>50 before cut R%d n20_eta_p20 0-200cent",list_radius[k]),50,0,50,100,0,60000);
    hpbpb_Npix_before_cut[k][nbins_cent+1] = new TH2F(Form("hpbpb_Npix_before_cut_R%d_n20_eta_p20_cent%d",list_radius[k],nbins_cent+1),Form("Number of pixels hit per no of jets pT>50 before cut R%d n20_eta_p20 ultraCentral 0 to 1 centrality",list_radius[k]),50,0,50,100,0,60000);
    hpbpb_Npix_after_cut[k][nbins_cent] = new TH2F(Form("hpbpb_Npix_after_cut_R%d_n20_eta_p20_cent%d",list_radius[k],nbins_cent),Form("Number of pixels hit per no of jets pT>50 after cut R%d n20_eta_p20 0-200cent",list_radius[k]),50,0,50,100,0,60000);

  }// radii loop
  
  
  TH1F *hJets_noTrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_3TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_5TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_7TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_10TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];

  for(int k = 0;k<no_radius;k++){
    for(int l = 0;l<trigValue;l++){
      for(int j = 0;j<nbins_eta;j++){
	for(int i = 0;i<nbins_cent;i++){
	  hJets_noTrkpTCut[k][j][i][l] = new TH1F(Form("hJets_noTrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with no fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_3TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_3TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 3GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_5TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_5TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 5GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_7TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_7TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 7GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_10TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_10TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 10GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  
	}// cent loop
	
      }// eta loop

    }// trigger loop

  }// radius loop

  // TTree *evt_electron_failure[no_radius];
  // TTree *evt_electron_good[no_radius];

  // Int_t event_value;
  // Int_t run_value;
  // Int_t lumi_value;
  // Int_t pthat_value;

  // for(int k = 0;k<no_radius;k++){

  //   evt_electron_failure[k] = new TTree (Form("evt_electron_failure_R%d",list_radius[k]),"");
  //   evt_electron_failure[k]->Branch("run_value",&run_value,"run_value/I");
  //   evt_electron_failure[k]->Branch("lumi_value",&lumi_value,"lumi_value/I");
  //   evt_electron_failure[k]->Branch("event_value",&event_value,"event_value/I");
  //   evt_electron_failure[k]->Branch("pthat_value",&pthat_value,"pthat_value/I");
  //   evt_electron_good[k] = new TTree (Form("evt_electron_good_R%d",list_radius[k]),"");
  //   evt_electron_good[k]->Branch("run_value",&run_value,"run_value/I");
  //   evt_electron_good[k]->Branch("lumi_value",&lumi_value,"lumi_value/I");
  //   evt_electron_good[k]->Branch("event_value",&event_value,"event_value/I");
  //   evt_electron_good[k]->Branch("pthat_value",&pthat_value,"pthat_value/I");
  // }


  // Setup jet data branches - this will be 2D with [radius][pthat-file], but the histogram here is just 1D with [radius]
  JetData *data[no_radius][nbins_pthat]; 
  for(int k = 0;k<no_radius;k++){
    if(printDebug)cout<<"Radius = "<<list_radius[k]<<endl;
    if(printDebug)cout<<"reading all the pbpb mc files"<<endl;
    for (int h=0;h<endfile;h++) {
      //cout<<Form("ak%s%dJetAnalyzer/t",algo,list_radius[k])<<endl;
      data[k][h] = new JetData(fileName_pthat[h],Form("ak%s%d%sJetAnalyzer/t",algo,list_radius[k],jet_type),Form("ak%s%d%sJetAnalyzer/t",algo,list_radius[k],jet_type),0,1);
      //cout<<"A"<<endl;
      TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbins_pthat,boundaries_pthat);
      //cout<<"B"<<endl;
      data[k][h]->tJet->Project("hPtHatTmp","pthat");
      //cout<<"C"<<endl;
      hPtHatRaw[k]->Add(hPtHatTmp);
      //cout<<"D"<<endl;
      delete hPtHatTmp;
    }// pthat loop
    
  }//radius loop

  // checking the histograms to see if something is filled. 
  if(printDebug)hPtHatRaw[0]->Print("base");
  // if(printDebug)hPtHatRawPP[0]->Print("base");

#if 0

  //plots for the jet ID (hopefully final check) // will be the same in Data and MC
  TH1F *hpbpb_Jet80 = new TH1F("hpbpb_Jet80","",100,0,300);
  TH1F* hpbpb_Jet80_chMaxJtpt0p01 = new TH1F("hpbpb_Jet80_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet80_chMaxJtpt0p02 = new TH1F("hpbpb_Jet80_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet80_chMaxJtpt0p03 = new TH1F("hpbpb_Jet80_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet80_chMaxJtpt0p04 = new TH1F("hpbpb_Jet80_chMaxJtpt0p04","",100,0,300);
  TH1F* hpbpb_Jet80_chMaxJtpt0p05 = new TH1F("hpbpb_Jet80_chMaxJtpt0p05","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p1 = new TH1F("hpbpb_Jet80_eMaxJtpt0p1","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p2 = new TH1F("hpbpb_Jet80_eMaxJtpt0p2","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p3 = new TH1F("hpbpb_Jet80_eMaxJtpt0p3","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p4 = new TH1F("hpbpb_Jet80_eMaxJtpt0p4","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p5 = new TH1F("hpbpb_Jet80_eMaxJtpt0p5","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p6 = new TH1F("hpbpb_Jet80_eMaxJtpt0p6","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p7 = new TH1F("hpbpb_Jet80_eMaxJtpt0p7","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p8 = new TH1F("hpbpb_Jet80_eMaxJtpt0p8","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p9 = new TH1F("hpbpb_Jet80_eMaxJtpt0p9","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p1 = new TH1F("hpbpb_Jet80_eMaxSumcand0p1","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p2 = new TH1F("hpbpb_Jet80_eMaxSumcand0p2","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p3 = new TH1F("hpbpb_Jet80_eMaxSumcand0p3","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p4 = new TH1F("hpbpb_Jet80_eMaxSumcand0p4","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p5 = new TH1F("hpbpb_Jet80_eMaxSumcand0p5","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p6 = new TH1F("hpbpb_Jet80_eMaxSumcand0p6","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p7 = new TH1F("hpbpb_Jet80_eMaxSumcand0p7","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p8 = new TH1F("hpbpb_Jet80_eMaxSumcand0p8","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxSumcand0p9 = new TH1F("hpbpb_Jet80_eMaxSumcand0p9","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01 = new TH1F("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02 = new TH1F("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03 = new TH1F("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01 = new TH1F("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02 = new TH1F("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03 = new TH1F("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01 = new TH1F("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02 = new TH1F("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03 = new TH1F("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03","",100,0,300);

  TH1F *hpbpb_Jet65 = new TH1F("hpbpb_Jet65","",100,0,300);
  TH1F* hpbpb_Jet65_chMaxJtpt0p01 = new TH1F("hpbpb_Jet65_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet65_chMaxJtpt0p02 = new TH1F("hpbpb_Jet65_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet65_chMaxJtpt0p03 = new TH1F("hpbpb_Jet65_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet65_chMaxJtpt0p04 = new TH1F("hpbpb_Jet65_chMaxJtpt0p04","",100,0,300);
  TH1F* hpbpb_Jet65_chMaxJtpt0p05 = new TH1F("hpbpb_Jet65_chMaxJtpt0p05","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p1 = new TH1F("hpbpb_Jet65_eMaxJtpt0p1","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p2 = new TH1F("hpbpb_Jet65_eMaxJtpt0p2","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p3 = new TH1F("hpbpb_Jet65_eMaxJtpt0p3","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p4 = new TH1F("hpbpb_Jet65_eMaxJtpt0p4","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p5 = new TH1F("hpbpb_Jet65_eMaxJtpt0p5","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p6 = new TH1F("hpbpb_Jet65_eMaxJtpt0p6","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p7 = new TH1F("hpbpb_Jet65_eMaxJtpt0p7","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p8 = new TH1F("hpbpb_Jet65_eMaxJtpt0p8","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p9 = new TH1F("hpbpb_Jet65_eMaxJtpt0p9","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p1 = new TH1F("hpbpb_Jet65_eMaxSumcand0p1","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p2 = new TH1F("hpbpb_Jet65_eMaxSumcand0p2","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p3 = new TH1F("hpbpb_Jet65_eMaxSumcand0p3","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p4 = new TH1F("hpbpb_Jet65_eMaxSumcand0p4","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p5 = new TH1F("hpbpb_Jet65_eMaxSumcand0p5","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p6 = new TH1F("hpbpb_Jet65_eMaxSumcand0p6","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p7 = new TH1F("hpbpb_Jet65_eMaxSumcand0p7","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p8 = new TH1F("hpbpb_Jet65_eMaxSumcand0p8","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxSumcand0p9 = new TH1F("hpbpb_Jet65_eMaxSumcand0p9","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01 = new TH1F("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02 = new TH1F("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03 = new TH1F("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01 = new TH1F("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02 = new TH1F("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03 = new TH1F("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01 = new TH1F("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02 = new TH1F("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03 = new TH1F("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03","",100,0,300);

  TH1F *hpbpb_Jet55 = new TH1F("hpbpb_Jet55","",100,0,300);
  TH1F* hpbpb_Jet55_chMaxJtpt0p01 = new TH1F("hpbpb_Jet55_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet55_chMaxJtpt0p02 = new TH1F("hpbpb_Jet55_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet55_chMaxJtpt0p03 = new TH1F("hpbpb_Jet55_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet55_chMaxJtpt0p04 = new TH1F("hpbpb_Jet55_chMaxJtpt0p04","",100,0,300);
  TH1F* hpbpb_Jet55_chMaxJtpt0p05 = new TH1F("hpbpb_Jet55_chMaxJtpt0p05","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p1 = new TH1F("hpbpb_Jet55_eMaxJtpt0p1","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p2 = new TH1F("hpbpb_Jet55_eMaxJtpt0p2","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p3 = new TH1F("hpbpb_Jet55_eMaxJtpt0p3","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p4 = new TH1F("hpbpb_Jet55_eMaxJtpt0p4","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p5 = new TH1F("hpbpb_Jet55_eMaxJtpt0p5","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p6 = new TH1F("hpbpb_Jet55_eMaxJtpt0p6","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p7 = new TH1F("hpbpb_Jet55_eMaxJtpt0p7","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p8 = new TH1F("hpbpb_Jet55_eMaxJtpt0p8","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p9 = new TH1F("hpbpb_Jet55_eMaxJtpt0p9","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p1 = new TH1F("hpbpb_Jet55_eMaxSumcand0p1","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p2 = new TH1F("hpbpb_Jet55_eMaxSumcand0p2","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p3 = new TH1F("hpbpb_Jet55_eMaxSumcand0p3","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p4 = new TH1F("hpbpb_Jet55_eMaxSumcand0p4","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p5 = new TH1F("hpbpb_Jet55_eMaxSumcand0p5","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p6 = new TH1F("hpbpb_Jet55_eMaxSumcand0p6","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p7 = new TH1F("hpbpb_Jet55_eMaxSumcand0p7","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p8 = new TH1F("hpbpb_Jet55_eMaxSumcand0p8","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxSumcand0p9 = new TH1F("hpbpb_Jet55_eMaxSumcand0p9","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01 = new TH1F("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02 = new TH1F("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03 = new TH1F("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01 = new TH1F("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02 = new TH1F("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03 = new TH1F("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01 = new TH1F("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02 = new TH1F("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02","",100,0,300);
  TH1F* hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03 = new TH1F("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03","",100,0,300);

#endif
  TH2F* hpbpb_chMaxJtpt_jtpt = new TH2F("hpbpb_chMaxJtpt_jtpt","",100,0,300,500,0,5);
  TH2F* hpbpb_eMaxJtpt_jtpt = new TH2F("hpbpb_eMaxJtpt_jtpt","",100,0,300,500,0,5);
  TH2F* hpbpb_eMaxSumcand_jtpt = new TH2F("hpbpb_eMaxSumcand_jtpt","",100,0,300,1000,0,10);
  TH2F* hpbpb_eMaxSumcand_refpt = new TH2F("hpbpb_eMaxSumcand_refpt","",100,0,300,1000,0,10);
  TH2F* hpbpb_eMaxJtpt_chMaxJtpt = new TH2F("hpbpb_eMaxJtpt_chMaxJtpt","",500,0,5,500,0,5);
  TH2F* hpbpb_eMaxSumcand_chMaxJtpt = new TH2F("hpbpb_eMaxSumcand_chMaxJtpt","",500,0,5,1000,0,10);

  static const int ptSelection = 19;
  static const int ptBoundary[ptSelection+1] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
  
  TH2F* hpbpb_chMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_phMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_neMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_muMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_jtpt_refptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_jtpt_refANDjtptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_refpt_refptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_chMaxJtpt_refANDjtptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_eMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_eMaxJtpt_refANDjtptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_neMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_phMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_muMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];

  for(int a = 0;a<ptSelection;++a){
    for(int t = 0;t<trigValue;++t){
      hpbpb_chMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_chMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_eMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_phMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_phMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_neMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_neMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_muMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_muMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_eMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_phMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_phMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_neMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_neMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_muMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_muMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_eMaxSumcand_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,1000,0,10);
      hpbpb_eMaxSumcand_chMaxJtpt_refANDjtptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_chMaxJtpt_%d_refANDjtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,1000,0,10);
      hpbpb_eMaxSumcand_eMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_eMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,1000,0,10);
      hpbpb_eMaxSumcand_eMaxJtpt_refANDjtptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_eMaxJtpt_%d_refANDjtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,1000,0,10);
      hpbpb_eMaxSumcand_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,1000,0,10);
      hpbpb_eMaxSumcand_jtpt_refptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_jtpt_%d_refpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,1000,0,10);
      hpbpb_eMaxSumcand_jtpt_refANDjtptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_jtpt_%d_refANDjtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,1000,0,10);
      hpbpb_eMaxSumcand_refpt_refptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_refpt_%d_refpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,1000,0,10);

    }

  }
  
  for(int k = 0;k<no_radius;k++){
    if(printDebug)cout<<"Filling MC for radius = "<<list_radius[k]<<endl;
    // fill PbPb MC 
    if(printDebug)cout<<"Filling PbPb MC"<<endl;

    for (int h=0;h<endfile;h++) {
      int goodCounter = 0;

      if (xsection[h]==0) continue;
      if(printDebug)cout <<"Loading pthat"<<boundaries_pthat[h]<<" sample, cross section = "<<xsection[h]<< Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[h],boundaries_pthat[h+1])<<endl;
      
      //from Pawan's code: /net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/combinePtHatBins/pbpbJEC2014/condor/CondorPbPbCalJec.C
      TEventList *el = new TEventList("el","el");
      //double pthat_event = data[k][h]->pthat;
      //double pthat_lower = boundaries_pthat[h];
      double pthat_upper = boundaries_pthat[h+1];
      stringstream selection; selection<<"pthat<"<<pthat_upper;
      
      data[k][h]->tJet->Draw(">>el",selection.str().c_str());
      double fentries = el->GetN();
      if(printDebug)cout<<"tree entries: "<<data[k][h]->tJet->GetEntries()<<" elist: "<<fentries<<endl;
      delete el;
      int test_counter = 0; 

      Int_t nEntries = data[k][h]->tJet->GetEntries();
      //if(printDebug) nEntries = 10;

      for (Long64_t jentry=0; jentry<nEntries;jentry++) {
	
        data[k][h]->tEvt->GetEntry(jentry);
        data[k][h]->tJet->GetEntry(jentry);
        data[k][h]->tSkim->GetEntry(jentry);
        
	if(data[k][h]->pthat > boundaries_pthat[h] && data[k][h]->pthat < boundaries_pthat[h+1]) test_counter++;

        //remember this cut is there because there was some rediculous values of pthats of -1 in the private production forests
        //if(jentry%100==0)cout<<"pthat of that event = "<<data[k][h]->pthat<<endl;
      
        int pthatBin = hPtHat[k]->FindBin(data[k][h]->pthat);

	double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/fentries;

	if(!data[k][h]->pcollisionEventSelection) continue;
        int cBin = findBin(data[k][h]->bin);
        //int cBin = nbins_cent-1;
        double weight_cent=1;
        double wegiht_pt=1;
        double weight_vz=1;
	Float_t cent = cBin;

        weight_cent = hCentWeight->GetBinContent(hCentWeight->FindBin(data[k][h]->bin));
	if(fabs(data[k][h]->vz)>15) continue;

        weight_vz = fVz->Eval(data[k][h]->vz);

	hpbpb_vz[k]->Fill(data[k][h]->vz,weight_vz);
	hpbpb_vx[k]->Fill(data[k][h]->vx);
	hpbpb_vy[k]->Fill(data[k][h]->vy);

	hpbpb_cent[k]->Fill(data[k][h]->bin,weight_cent);

	if(scale*weight_cent*weight_vz <=0 ) {
	  cout<<"RED FLAG RED FLAG RED FLAG"<<endl;
	  cout<<"pthat file = "<<boundaries_pthat[h]<<endl;
	  continue;
	}

	int jetCounter = 0;//counts jets which are going to be used in the supernova cut rejection. 
	
	for(int j = 0;j<nbins_eta;j++){

	  for(int g = 0;g<data[k][h]->njets;g++){

	    if(data[k][h]->jteta[g] >= boundaries_eta[j][0] && data[k][h]->jteta[g] < boundaries_eta[j][1]){
	      if(data[k][h]->jtpt[g]>=50) jetCounter++;
	      
	    }// eta selection loop

	  }//jet loop

	}//eta bins loop

	// apply the supernova events cut rejection here: 
	if(data[k][h]->hiNpix > 38000 - 500*jetCounter){
	  if(printDebug) cout<<"removed this supernova event"<<endl;
	  continue;
	}

        if (cBin>=nbins_cent) continue;
        if (cBin==-1) continue;
	// hPtHat[k]->Fill(data[k][h]->pthat,scale*weight_cent*weight_vz);
	
        //cout<<"scale = "<<scale<<endl;
	
#if 0    
	int hasLeadingJet = 0;
	for (int k= 0; k < data[h]->njets; k++) { 
	if ( data[h]->jteta[k]  > 2. || data[h]->jteta[k] < -2. ) continue;
	if ( data[h]->jtpt[k]>100) {
	hasLeadingJet = 1;
      }
	break;
	
      }
	if (hasLeadingJet == 0) continue;
#endif


	// if(data[k][h]->jet55_1 && (data[k][h]->eMax[0]/data[k][h]->jtpt[0]) > 0.8 && data[k][h]->jtpt[0]>80){
	//   event_value = data[k][h]->evt;
	//   run_value = data[k][h]->run;
	//   lumi_value = data[k][h]->lumi;
	//   pthat_value = boundaries_pthat[h];
	//   evt_electron_failure[k]->Fill();
	// }
	// if(data[k][h]->jet55_1 && (data[k][h]->eMax[0]/data[k][h]->jtpt[0]) < 0.4 && data[k][h]->jtpt[0]>80){
	//   event_value = data[k][h]->evt;
	//   run_value = data[k][h]->run;
	//   lumi_value = data[k][h]->lumi;
	//   pthat_value = boundaries_pthat[h];
	//   evt_electron_good[k]->Fill();
	// }
	
	for (int g = 0; g < data[k][h]->njets; g++) {

	  //removed the following for no cut histogram: 
	  //if ( data[k][h]->rawpt[g] < 30. ) continue;
	  if ( data[k][h]->subid[g] != sub_id ) continue;
	  if ( data[k][h]->jtpt[g] > 2.*data[k][h]->pthat) continue;
	  // jet quality cuts here:
	  //if ( data[k][h]->neutralMax[g]/TMath::Max(data[h]->chargedSum[k],data[h]->neutralSum[k]) < 0.975)continue;
	  //if ( data[k][h]->neutralMax[g]/(data[k][h]->chargedMax[g] + data[k][h]->photonMax[g] + data[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( data[k][h]->photonMax[g]/(data[k][h]->chargedMax[g] + data[k][h]->photonMax[g] + data[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( data[k][h]->chargedMax[g]/(data[k][h]->chargedMax[g] + data[k][h]->photonMax[g] + data[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( data[k][h]->muMax[g]/(data[k][h]->chargedMax[g] + data[k][h]->photonMax[g] + data[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( data[k][h]->eSum[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]) > 0.7 )continue;
	  //if(((data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->eSum[g] + data[k][h]->muSum[g])/data[k][h]->jtpt[g])>1.01) continue;

	  hpbpb_eta_full[k]->Fill(data[k][h]->jteta[g],scale*weight_vz*weight_cent);
	  hpbpb_phi_full[k]->Fill(data[k][h]->jtphi[g],scale*weight_vz*weight_cent);
	  hpbpb_eta_full_noScale[k]->Fill(data[k][h]->jteta[g]);
	  hpbpb_phi_full_noScale[k]->Fill(data[k][h]->jtphi[g]);
  
          for(int j = 0;j<nbins_eta;j++){

#if 0
	    vector <Float_t> inJetPFcand_pT;

            //int subEvt=-1;
	    
            if ( data[k][h]->jteta[g]  > boundaries_eta[j][1] || data[k][h]->jteta[g] < boundaries_eta[j][0] ) continue;
	    
	    // have to run through all the particle flow candidates:
	    for(int ipf = 0; ipf<data[k][h]->nPFpart; ipf++){

	      if(TMath::Sqrt((data[k][h]->jteta[g] - data[k][h]->pfEta[ipf])*(data[k][h]->jteta[g] - data[k][h]->pfEta[ipf]) + (data[k][h]->jtphi[g] - data[k][h]->pfPhi[ipf])*(data[k][h]->jtphi[g] - data[k][h]->pfPhi[ipf])) <= (float) list_radius[k]/10){

		inJetPFcand_pT.push_back(data[k][h]->pfPt[ipf]);
	      
	      }// checking for candidates to be inside the jet cone 

	    } // pf flow candidate loop

	    // Now find the pT of the smallest candidate inside the jet.
	    float large = inJetPFcand_pT[0];
	    for(int a = 1;a<inJetPFcand_pT.size();a++){
	      if(large < inJetPFcand_pT[a])
		large = inJetPFcand_pT[a];
	    }
	  
	    //for(int a = 0;a<inJetPFcand_pT.size();a++) cout<<inJetPFcand_pT[a]<<", ";
	    //cout<<endl;

	    // test checking the sum of the pf candidate to equal the jet pt
	    if(printDebug) cout<<"jtpt = "<<data[k][h]->pfPt[g]<<", large candidate = "<<large<<endl;
	  
	    // now that we have the smallest candidate we can check for the different fragmentation biases:
	 	  
	    if(data[k][h]->jet80_1==1) {

	      hJets_noTrkpTCut[k][j][cBin][2]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 3) hJets_3TrkpTCut[k][j][cBin][2]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 5) hJets_5TrkpTCut[k][j][cBin][2]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 7) hJets_7TrkpTCut[k][j][cBin][2]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 10) hJets_10TrkpTCut[k][j][cBin][2]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      
	    }

	    if(data[k][h]->jet65_1==1){

	      hJets_noTrkpTCut[k][j][cBin][1]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 3) hJets_3TrkpTCut[k][j][cBin][1]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 5) hJets_5TrkpTCut[k][j][cBin][1]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 7) hJets_7TrkpTCut[k][j][cBin][1]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 10) hJets_10TrkpTCut[k][j][cBin][1]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);

	    }

	    if(data[k][h]->jet55_1==1 && data[k][h]->jet65_1==0 && data[k][h]->jet80_1==0){
	    
	      hJets_noTrkpTCut[k][j][cBin][0]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 3) hJets_3TrkpTCut[k][j][cBin][0]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 5) hJets_5TrkpTCut[k][j][cBin][0]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 7) hJets_7TrkpTCut[k][j][cBin][0]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);
	      if(large > 10) hJets_10TrkpTCut[k][j][cBin][0]->Fill(data[k][h]->pfPt[g],scale*weight_vz*weight_cent);

	    }
	
	    for(int t = 0;t<trigValue-1;t++) {

	      hJets_noTrkpTCut[k][j][cBin][trigValue-1]->Add(hJets_noTrkpTCut[k][j][cBin][t]);
	      hJets_3TrkpTCut[k][j][cBin][trigValue-1]->Add(hJets_3TrkpTCut[k][j][cBin][t]);
	      hJets_5TrkpTCut[k][j][cBin][trigValue-1]->Add(hJets_5TrkpTCut[k][j][cBin][t]);
	      hJets_7TrkpTCut[k][j][cBin][trigValue-1]->Add(hJets_7TrkpTCut[k][j][cBin][t]);
	      hJets_10TrkpTCut[k][j][cBin][trigValue-1]->Add(hJets_10TrkpTCut[k][j][cBin][t]);
	    
	    }

#endif
	    //get the spectra histograms to make the ratios for the different Jet ID cuts. 0-30%
	    //if(cBin==3 || cBin==4 || cBin==5) continue;

	    if(k==0){
	      hpbpb_chMaxJtpt_jtpt->Fill(data[k][h]->jtpt[g],data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]);
	      hpbpb_eMaxJtpt_jtpt->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/data[k][h]->jtpt[g]);
	      hpbpb_eMaxSumcand_chMaxJtpt->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
	      hpbpb_eMaxSumcand_jtpt->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
	      hpbpb_eMaxJtpt_chMaxJtpt->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/data[k][h]->jtpt[g]);

	      for(int a = 0;a<ptSelection;a++) {
		
		if(data[k][h]->jet80_1){
		  if(data[k][h]->jtpt[g] >= ptBoundary[a] && data[k][h]->jtpt[g] < ptBoundary[a+1]){
		    hpbpb_eMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_neMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->neutralMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_phMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->photonMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_muMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->muMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_eMaxSumcand_chMaxJtpt_ptselection[2][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxSumcand_jtpt_ptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxJtpt_jtpt_ptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_chMaxJtpt_jtpt_ptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->chargedMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_phMaxJtpt_jtpt_ptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->photonMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_neMaxJtpt_jtpt_ptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->neutralMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_muMaxJtpt_jtpt_ptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->muMax[g]/(data[k][h]->jtpt[g]));
		    
		    hpbpb_eMaxSumcand_eMaxJtpt_ptselection[2][a]->Fill(data[k][h]->eMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    if(data[k][h]->refpt[g] >= ptBoundary[a] && data[k][h]->refpt[g] < ptBoundary[a+1]){
		      hpbpb_eMaxSumcand_chMaxJtpt_refANDjtptselection[2][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		      hpbpb_eMaxSumcand_eMaxJtpt_refANDjtptselection[2][a]->Fill(data[k][h]->eMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		      hpbpb_eMaxSumcand_jtpt_refANDjtptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    }
		   
		  }
		  if(data[k][h]->refpt[g] >= ptBoundary[a] && data[k][h]->refpt[g] < ptBoundary[a+1]){
		    hpbpb_eMaxSumcand_jtpt_refptselection[2][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		  }
		  
		}

		if(data[k][h]->jet65_1 && !data[k][h]->jet80_1){

		  if(data[k][h]->jtpt[g] >= ptBoundary[a] && data[k][h]->jtpt[g] < ptBoundary[a+1]){
		    hpbpb_eMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_neMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->neutralMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_phMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->photonMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_muMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->muMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_eMaxSumcand_chMaxJtpt_ptselection[1][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxSumcand_jtpt_ptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxSumcand_jtpt_ptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxJtpt_jtpt_ptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_chMaxJtpt_jtpt_ptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->chargedMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_phMaxJtpt_jtpt_ptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->photonMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_neMaxJtpt_jtpt_ptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->neutralMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_muMaxJtpt_jtpt_ptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->muMax[g]/(data[k][h]->jtpt[g]));

		    hpbpb_eMaxSumcand_eMaxJtpt_ptselection[1][a]->Fill(data[k][h]->eMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    if(data[k][h]->refpt[g] >= ptBoundary[a] && data[k][h]->refpt[g] < ptBoundary[a+1]){
		      hpbpb_eMaxSumcand_chMaxJtpt_refANDjtptselection[1][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		      hpbpb_eMaxSumcand_eMaxJtpt_refANDjtptselection[1][a]->Fill(data[k][h]->eMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		      hpbpb_eMaxSumcand_jtpt_refANDjtptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    }
		    
		  }
		  if(data[k][h]->refpt[g] >= ptBoundary[a] && data[k][h]->refpt[g] < ptBoundary[a+1]){
		    hpbpb_eMaxSumcand_jtpt_refptselection[1][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		  }
		}


		if(data[k][h]->jet55_1 && !data[k][h]->jet65_1 && !data[k][h]->jet80_1){

		  if(data[k][h]->jtpt[g] >= ptBoundary[a] && data[k][h]->jtpt[g] < ptBoundary[a+1]){
		    hpbpb_eMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_neMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->neutralMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_phMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->photonMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_muMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->muMax[g]/data[k][h]->jtpt[g]);
		    hpbpb_eMaxSumcand_chMaxJtpt_ptselection[0][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxSumcand_jtpt_ptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxSumcand_jtpt_ptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    hpbpb_eMaxJtpt_jtpt_ptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_chMaxJtpt_jtpt_ptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->chargedMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_phMaxJtpt_jtpt_ptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->photonMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_neMaxJtpt_jtpt_ptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->neutralMax[g]/(data[k][h]->jtpt[g]));
		    hpbpb_muMaxJtpt_jtpt_ptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->muMax[g]/(data[k][h]->jtpt[g]));

		    hpbpb_eMaxSumcand_eMaxJtpt_ptselection[0][a]->Fill(data[k][h]->eMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    if(data[k][h]->refpt[g] >= ptBoundary[a] && data[k][h]->refpt[g] < ptBoundary[a+1]){
		      hpbpb_eMaxSumcand_chMaxJtpt_refANDjtptselection[0][a]->Fill(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		      hpbpb_eMaxSumcand_eMaxJtpt_refANDjtptselection[0][a]->Fill(data[k][h]->eMax[g]/data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		      hpbpb_eMaxSumcand_jtpt_refANDjtptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		    }
		  }
		  if(data[k][h]->refpt[g] >= ptBoundary[a] && data[k][h]->refpt[g] < ptBoundary[a+1]){
		    hpbpb_eMaxSumcand_jtpt_refptselection[0][a]->Fill(data[k][h]->jtpt[g],data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g]));
		  }
		}

	      }

#if 0
	      if(data[k][h]->jet80_1){
		hpbpb_Jet80->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01)hpbpb_Jet80_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02)hpbpb_Jet80_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03)hpbpb_Jet80_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.04)hpbpb_Jet80_chMaxJtpt0p04->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.05)hpbpb_Jet80_chMaxJtpt0p05->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);


		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.1)hpbpb_Jet80_eMaxSumcand0p1->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.2)hpbpb_Jet80_eMaxSumcand0p2->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.3)hpbpb_Jet80_eMaxSumcand0p3->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.4)hpbpb_Jet80_eMaxSumcand0p4->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.5)hpbpb_Jet80_eMaxSumcand0p5->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.6)hpbpb_Jet80_eMaxSumcand0p6->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.7)hpbpb_Jet80_eMaxSumcand0p7->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.8)hpbpb_Jet80_eMaxSumcand0p8->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.9)hpbpb_Jet80_eMaxSumcand0p9->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);


		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.1)hpbpb_Jet80_eMaxJtpt0p1->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.2)hpbpb_Jet80_eMaxJtpt0p2->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.3)hpbpb_Jet80_eMaxJtpt0p3->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.4)hpbpb_Jet80_eMaxJtpt0p4->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet80_eMaxJtpt0p5->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet80_eMaxJtpt0p6->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet80_eMaxJtpt0p7->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.8)hpbpb_Jet80_eMaxJtpt0p8->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.9)hpbpb_Jet80_eMaxJtpt0p9->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      }
	      if(data[k][h]->jet65_1 && !data[k][h]->jet80_1){
		hpbpb_Jet65->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01)hpbpb_Jet65_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02)hpbpb_Jet65_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03)hpbpb_Jet65_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.04)hpbpb_Jet65_chMaxJtpt0p04->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.05)hpbpb_Jet65_chMaxJtpt0p05->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.1)hpbpb_Jet65_eMaxJtpt0p1->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.2)hpbpb_Jet65_eMaxJtpt0p2->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.3)hpbpb_Jet65_eMaxJtpt0p3->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.4)hpbpb_Jet65_eMaxJtpt0p4->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet65_eMaxJtpt0p5->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet65_eMaxJtpt0p6->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet65_eMaxJtpt0p7->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.8)hpbpb_Jet65_eMaxJtpt0p8->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.9)hpbpb_Jet65_eMaxJtpt0p9->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

		
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.1)hpbpb_Jet65_eMaxSumcand0p1->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.2)hpbpb_Jet65_eMaxSumcand0p2->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.3)hpbpb_Jet65_eMaxSumcand0p3->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.4)hpbpb_Jet65_eMaxSumcand0p4->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.5)hpbpb_Jet65_eMaxSumcand0p5->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.6)hpbpb_Jet65_eMaxSumcand0p6->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.7)hpbpb_Jet65_eMaxSumcand0p7->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.8)hpbpb_Jet65_eMaxSumcand0p8->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.9)hpbpb_Jet65_eMaxSumcand0p9->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);


		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      }
	      if(data[k][h]->jet55_1 && !data[k][h]->jet65_1 && !data[k][h]->jet80_1){
		hpbpb_Jet55->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01)hpbpb_Jet55_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02)hpbpb_Jet55_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03)hpbpb_Jet55_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.04)hpbpb_Jet55_chMaxJtpt0p04->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.05)hpbpb_Jet55_chMaxJtpt0p05->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.1)hpbpb_Jet55_eMaxSumcand0p1->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.2)hpbpb_Jet55_eMaxSumcand0p2->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.3)hpbpb_Jet55_eMaxSumcand0p3->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.4)hpbpb_Jet55_eMaxSumcand0p4->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.5)hpbpb_Jet55_eMaxSumcand0p5->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.6)hpbpb_Jet55_eMaxSumcand0p6->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.7)hpbpb_Jet55_eMaxSumcand0p7->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.8)hpbpb_Jet55_eMaxSumcand0p8->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/(data[k][h]->chargedSum[g] + data[k][h]->photonSum[g] + data[k][h]->neutralSum[g] + data[k][h]->muSum[g])<0.9)hpbpb_Jet55_eMaxSumcand0p9->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.1)hpbpb_Jet55_eMaxJtpt0p1->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.2)hpbpb_Jet55_eMaxJtpt0p2->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.3)hpbpb_Jet55_eMaxJtpt0p3->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.4)hpbpb_Jet55_eMaxJtpt0p4->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet55_eMaxJtpt0p5->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet55_eMaxJtpt0p6->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet55_eMaxJtpt0p7->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.8)hpbpb_Jet55_eMaxJtpt0p8->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.9)hpbpb_Jet55_eMaxJtpt0p9->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.7)hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.6)hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.01 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.02 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		if(data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]>0.03 && data[k][h]->eMax[g]/data[k][h]->jtpt[g]<0.5)hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      }
#endif
	    }
	  
#if 0
	    
	    hpbpb_RecoOverRaw[k][j][cBin]->Fill((Float_t)data[k][h]->jtpt[g]/data[k][h]->rawpt[g]);
	    hpbpb_RecoOverRaw[k][j][nbins_cent]->Fill((Float_t)data[k][h]->jtpt[g]/data[k][h]->rawpt[g]);
	  
	    hpbpb_RecoOverRaw_jtpt[k][j][cBin]->Fill(data[k][h]->jtpt[g],(Float_t)data[k][h]->jtpt[g]/data[k][h]->rawpt[g]);
	    hpbpb_RecoOverRaw_jtpt[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],(Float_t)data[k][h]->jtpt[g]/data[k][h]->rawpt[g]);
#endif

	    
	    hpbpb_chMax[k][j][cBin]->Fill(data[k][h]->chargedMax[g],weight_cent*weight_vz);
	    hpbpb_phMax[k][j][cBin]->Fill(data[k][h]->photonMax[g],weight_cent*weight_vz);
	    hpbpb_neMax[k][j][cBin]->Fill(data[k][h]->neutralMax[g],weight_cent*weight_vz);
	    hpbpb_muMax[k][j][cBin]->Fill(data[k][h]->muMax[g],weight_cent*weight_vz);
	    hpbpb_eMax[k][j][cBin]->Fill(data[k][h]->eMax[g],weight_cent*weight_vz);
	    
	    hpbpb_chSum[k][j][cBin]->Fill(data[k][h]->chargedSum[g],weight_cent*weight_vz);
	    hpbpb_phSum[k][j][cBin]->Fill(data[k][h]->photonSum[g],weight_cent*weight_vz);
	    hpbpb_neSum[k][j][cBin]->Fill(data[k][h]->neutralSum[g],weight_cent*weight_vz);
	    hpbpb_muSum[k][j][cBin]->Fill(data[k][h]->muSum[g],weight_cent*weight_vz);
	    hpbpb_eSum[k][j][cBin]->Fill(data[k][h]->eSum[g],weight_cent*weight_vz);

	    // hpbpb_chMax[k][j][nbins_cent]->Fill(data[k][h]->chargedMax[g],weight_cent*weight_vz);
	    // hpbpb_phMax[k][j][nbins_cent]->Fill(data[k][h]->photonMax[g],weight_cent*weight_vz);
	    // hpbpb_neMax[k][j][nbins_cent]->Fill(data[k][h]->neutralMax[g],weight_cent*weight_vz);
	    // hpbpb_muMax[k][j][nbins_cent]->Fill(data[k][h]->muMax[g],weight_cent*weight_vz);
	    // hpbpb_eMax[k][j][nbins_cent]->Fill(data[k][h]->eMax[g],weight_cent*weight_vz);
	    
	    // hpbpb_chSum[k][j][nbins_cent]->Fill(data[k][h]->chargedSum[g],weight_cent*weight_vz);
	    // hpbpb_phSum[k][j][nbins_cent]->Fill(data[k][h]->photonSum[g],weight_cent*weight_vz);
	    // hpbpb_neSum[k][j][nbins_cent]->Fill(data[k][h]->neutralSum[g],weight_cent*weight_vz);
	    // hpbpb_muSum[k][j][nbins_cent]->Fill(data[k][h]->muSum[g],weight_cent*weight_vz);
	    // hpbpb_eSum[k][j][nbins_cent]->Fill(data[k][h]->eSum[g],weight_cent*weight_vz);

	    if ( data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]<0.02 || data[k][h]->eMax[g]/data[k][h]->jtpt[g]>0.6 ) continue;

	    hpbpb_chMax_withCut[k][j][cBin]->Fill(data[k][h]->chargedMax[g],weight_cent*weight_vz);
	    hpbpb_phMax_withCut[k][j][cBin]->Fill(data[k][h]->photonMax[g],weight_cent*weight_vz);
	    hpbpb_neMax_withCut[k][j][cBin]->Fill(data[k][h]->neutralMax[g],weight_cent*weight_vz);
	    hpbpb_muMax_withCut[k][j][cBin]->Fill(data[k][h]->muMax[g],weight_cent*weight_vz);
	    hpbpb_eMax_withCut[k][j][cBin]->Fill(data[k][h]->eMax[g],weight_cent*weight_vz);

	    hpbpb_chSum_withCut[k][j][cBin]->Fill(data[k][h]->chargedSum[g],weight_cent*weight_vz);
	    hpbpb_phSum_withCut[k][j][cBin]->Fill(data[k][h]->photonSum[g],weight_cent*weight_vz);
	    hpbpb_neSum_withCut[k][j][cBin]->Fill(data[k][h]->neutralSum[g],weight_cent*weight_vz);
	    hpbpb_muSum_withCut[k][j][cBin]->Fill(data[k][h]->muSum[g],weight_cent*weight_vz);
	    hpbpb_eSum_withCut[k][j][cBin]->Fill(data[k][h]->eSum[g],weight_cent*weight_vz);

	    // hpbpb_chMax_withCut[k][j][nbins_cent]->Fill(data[k][h]->chargedMax[g],weight_cent*weight_vz);
	    // hpbpb_phMax_withCut[k][j][nbins_cent]->Fill(data[k][h]->photonMax[g],weight_cent*weight_vz);
	    // hpbpb_neMax_withCut[k][j][nbins_cent]->Fill(data[k][h]->neutralMax[g],weight_cent*weight_vz);
	    // hpbpb_muMax_withCut[k][j][nbins_cent]->Fill(data[k][h]->muMax[g],weight_cent*weight_vz);
	    // hpbpb_eMax_withCut[k][j][nbins_cent]->Fill(data[k][h]->eMax[g],weight_cent*weight_vz);

	    // hpbpb_chSum_withCut[k][j][nbins_cent]->Fill(data[k][h]->chargedSum[g],weight_cent*weight_vz);
	    // hpbpb_phSum_withCut[k][j][nbins_cent]->Fill(data[k][h]->photonSum[g],weight_cent*weight_vz);
	    // hpbpb_neSum_withCut[k][j][nbins_cent]->Fill(data[k][h]->neutralSum[g],weight_cent*weight_vz);
	    // hpbpb_muSum_withCut[k][j][nbins_cent]->Fill(data[k][h]->muSum[g],weight_cent*weight_vz);
	    // hpbpb_eSum_withCut[k][j][nbins_cent]->Fill(data[k][h]->eSum[g],weight_cent*weight_vz);

	    if(jetCounter>=7) hpbpb_pt_Njet_g7[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	    if(jetCounter<7) hpbpb_pt_Njet_l7[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	    hpbpb_matrix[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	    hpbpb_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	    hpbpb_reco[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	    hpbpb_jtpu[k][j][cBin]->Fill(data[k][h]->jtpu[g],scale*weight_vz*weight_cent);	    
	    hpbpb_jtpu_noScale[k][j][cBin]->Fill(data[k][h]->jtpu[g],1);

	    //hpbpb_response[nbins_cent]->Fill(data[h]->jtpt[k],data[h]->refpt[k],scale*weight_vz);
	    hpbpb_matrix[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	    hpbpb_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	    hpbpb_reco[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	    
	    hpbpb_jtpu[k][j][nbins_cent]->Fill(data[k][h]->jtpu[g],scale*weight_vz*weight_cent);
	    hpbpb_jtpu_noScale[k][j][nbins_cent]->Fill(data[k][h]->jtpu[g],1);		    

	    if(jentry%2==0) {
	      hpbpb_mcclosure_data[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_mcclosure_data[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	      if(data[k][h]->jet80_1){
		hpbpb_mcclosure_Jet80_data[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_Jet80_data[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	      }else if(data[k][h]->jet65_1 && !data[k][h]->jet80_1){

		hpbpb_mcclosure_Jet65_data[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_Jet65_data[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	      }else if(data[k][h]->jet55_1 && !data[k][h]->jet65_1 && !data[k][h]->jet80_1){

		hpbpb_mcclosure_Jet55_data[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_Jet55_data[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	      }

	    }

	    if(jentry%2==1) {
	      hpbpb_mcclosure_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_mcclosure_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_mcclosure_matrix[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_mcclosure_matrix[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	      if(data[k][h]->jet80_1){
		hpbpb_mcclosure_Jet80_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_Jet80_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_matrix_HLT[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_matrix_HLT[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      }else if(data[k][h]->jet65_1 && !data[k][h]->jet80_1){

		hpbpb_mcclosure_Jet65_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_Jet65_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_matrix_HLT[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_matrix_HLT[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      }else if(data[k][h]->jet55_1 && !data[k][h]->jet65_1 && !data[k][h]->jet80_1){
		
		hpbpb_mcclosure_Jet55_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_Jet55_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_matrix_HLT[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
		hpbpb_mcclosure_matrix_HLT[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      }
	      
	    }

	    if(data[k][h]->jet80_1){

	      hpbpb_Jet80_reco[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet80_reco[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet80_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet80_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_matrix_HLT[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_matrix_HLT[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	    
	    }else if(data[k][h]->jet65_1 && !data[k][h]->jet80_1){

	      hpbpb_Jet65_reco[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet65_reco[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet65_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet65_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_matrix_HLT[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_matrix_HLT[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz*weight_cent);

	    }else if(data[k][h]->jet55_1 && !data[k][h]->jet65_1 && !data[k][h]->jet80_1){

	      hpbpb_Jet55_reco[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet55_reco[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet55_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_Jet55_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz*weight_cent);
	      hpbpb_matrix_HLT[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],data[k][h]->jet55_p_1*scale*weight_vz*weight_cent);
	      hpbpb_matrix_HLT[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],data[k][h]->jet55_p_1*scale*weight_vz*weight_cent);

	    }
	    
          }// eta bins loop
	      
        }//njets loop
	
      }//nentry loop

      if(printDebug)cout<<"no of events inbetween pthat "<<boundaries_pthat[h]<<" and "<<boundaries_pthat[h+1]<<" = "<<test_counter<<endl;

    }//ptbins loop
    

  }// radius loop

  for(int k = 0;k<no_radius;k++){
    for(int j = 0;j<nbins_eta;j++){
      for(int i = 0;i<nbins_cent;i++){
	hpbpb_JetComb_gen[k][j][i]->Add(hpbpb_Jet80_gen[k][j][i]);
	hpbpb_JetComb_gen[k][j][i]->Add(hpbpb_Jet65_gen[k][j][i]);
	hpbpb_JetComb_gen[k][j][i]->Add(hpbpb_Jet55_gen[k][j][i]);

	hpbpb_JetComb_reco[k][j][i]->Add(hpbpb_Jet80_reco[k][j][i]);
	hpbpb_JetComb_reco[k][j][i]->Add(hpbpb_Jet65_reco[k][j][i]);
	hpbpb_JetComb_reco[k][j][i]->Add(hpbpb_Jet55_reco[k][j][i]);

	hpbpb_mcclosure_JetComb_gen[k][j][i]->Add(hpbpb_mcclosure_Jet80_gen[k][j][i]);
	hpbpb_mcclosure_JetComb_gen[k][j][i]->Add(hpbpb_mcclosure_Jet65_gen[k][j][i]);
	hpbpb_mcclosure_JetComb_gen[k][j][i]->Add(hpbpb_mcclosure_Jet55_gen[k][j][i]);

	hpbpb_mcclosure_JetComb_data[k][j][i]->Add(hpbpb_mcclosure_Jet80_data[k][j][i]);
	hpbpb_mcclosure_JetComb_data[k][j][i]->Add(hpbpb_mcclosure_Jet65_data[k][j][i]);
	hpbpb_mcclosure_JetComb_data[k][j][i]->Add(hpbpb_mcclosure_Jet55_data[k][j][i]);
		
      }
    }
  }
  

  TFile f(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_mc_ak%s%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  f.cd();
#if 0
  hpbpb_Jet80_chMaxJtpt0p01->Divide(hpbpb_Jet80);
  hpbpb_Jet80_chMaxJtpt0p01->Write();
  hpbpb_Jet80_chMaxJtpt0p02->Divide(hpbpb_Jet80);
  hpbpb_Jet80_chMaxJtpt0p02->Write();
  hpbpb_Jet80_chMaxJtpt0p03->Divide(hpbpb_Jet80);
  hpbpb_Jet80_chMaxJtpt0p03->Write();
  hpbpb_Jet80_chMaxJtpt0p04->Divide(hpbpb_Jet80);
  hpbpb_Jet80_chMaxJtpt0p04->Write();
  hpbpb_Jet80_chMaxJtpt0p05->Divide(hpbpb_Jet80);
  hpbpb_Jet80_chMaxJtpt0p05->Write();

  hpbpb_Jet80_eMaxJtpt0p1->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p1->Write();
  hpbpb_Jet80_eMaxJtpt0p2->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p2->Write();
  hpbpb_Jet80_eMaxJtpt0p3->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p3->Write();
  hpbpb_Jet80_eMaxJtpt0p4->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p4->Write();
  hpbpb_Jet80_eMaxJtpt0p5->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p5->Write();
  hpbpb_Jet80_eMaxJtpt0p6->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p6->Write();
  hpbpb_Jet80_eMaxJtpt0p7->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p7->Write();
  hpbpb_Jet80_eMaxJtpt0p8->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p8->Write();
  hpbpb_Jet80_eMaxJtpt0p9->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p9->Write();

  hpbpb_Jet80_eMaxSumcand0p1->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p1->Write();
  hpbpb_Jet80_eMaxSumcand0p2->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p2->Write();
  hpbpb_Jet80_eMaxSumcand0p3->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p3->Write();
  hpbpb_Jet80_eMaxSumcand0p4->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p4->Write();
  hpbpb_Jet80_eMaxSumcand0p5->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p5->Write();
  hpbpb_Jet80_eMaxSumcand0p6->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p6->Write();
  hpbpb_Jet80_eMaxSumcand0p7->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p7->Write();
  hpbpb_Jet80_eMaxSumcand0p8->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p8->Write();
  hpbpb_Jet80_eMaxSumcand0p9->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxSumcand0p9->Write();

  hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Write();
  hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Write();
  hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Write();
  hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Write();
  hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Write();
  hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Write();
  hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Write();
  hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Write();
  hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hpbpb_Jet80);
  hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Write();

  
  hpbpb_Jet65_chMaxJtpt0p01->Divide(hpbpb_Jet65);
  hpbpb_Jet65_chMaxJtpt0p01->Write();
  hpbpb_Jet65_chMaxJtpt0p02->Divide(hpbpb_Jet65);
  hpbpb_Jet65_chMaxJtpt0p02->Write();
  hpbpb_Jet65_chMaxJtpt0p03->Divide(hpbpb_Jet65);
  hpbpb_Jet65_chMaxJtpt0p03->Write();
  hpbpb_Jet65_chMaxJtpt0p04->Divide(hpbpb_Jet65);
  hpbpb_Jet65_chMaxJtpt0p04->Write();
  hpbpb_Jet65_chMaxJtpt0p05->Divide(hpbpb_Jet65);
  hpbpb_Jet65_chMaxJtpt0p05->Write();

  hpbpb_Jet65_eMaxJtpt0p1->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p1->Write();
  hpbpb_Jet65_eMaxJtpt0p2->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p2->Write();
  hpbpb_Jet65_eMaxJtpt0p3->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p3->Write();
  hpbpb_Jet65_eMaxJtpt0p4->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p4->Write();
  hpbpb_Jet65_eMaxJtpt0p5->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p5->Write();
  hpbpb_Jet65_eMaxJtpt0p6->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p6->Write();
  hpbpb_Jet65_eMaxJtpt0p7->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p7->Write();
  hpbpb_Jet65_eMaxJtpt0p8->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p8->Write();
  hpbpb_Jet65_eMaxJtpt0p9->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p9->Write();

  hpbpb_Jet65_eMaxSumcand0p1->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p1->Write();
  hpbpb_Jet65_eMaxSumcand0p2->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p2->Write();
  hpbpb_Jet65_eMaxSumcand0p3->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p3->Write();
  hpbpb_Jet65_eMaxSumcand0p4->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p4->Write();
  hpbpb_Jet65_eMaxSumcand0p5->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p5->Write();
  hpbpb_Jet65_eMaxSumcand0p6->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p6->Write();
  hpbpb_Jet65_eMaxSumcand0p7->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p7->Write();
  hpbpb_Jet65_eMaxSumcand0p8->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p8->Write();
  hpbpb_Jet65_eMaxSumcand0p9->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxSumcand0p9->Write();

  hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Write();
  hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Write();
  hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Write();
  hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Write();
  hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Write();
  hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Write();
  hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Write();
  hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Write();
  hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hpbpb_Jet65);
  hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Write();

  
  hpbpb_Jet55_chMaxJtpt0p01->Divide(hpbpb_Jet55);
  hpbpb_Jet55_chMaxJtpt0p01->Write();
  hpbpb_Jet55_chMaxJtpt0p02->Divide(hpbpb_Jet55);
  hpbpb_Jet55_chMaxJtpt0p02->Write();
  hpbpb_Jet55_chMaxJtpt0p03->Divide(hpbpb_Jet55);
  hpbpb_Jet55_chMaxJtpt0p03->Write();
  hpbpb_Jet55_chMaxJtpt0p04->Divide(hpbpb_Jet55);
  hpbpb_Jet55_chMaxJtpt0p04->Write();
  hpbpb_Jet55_chMaxJtpt0p05->Divide(hpbpb_Jet55);
  hpbpb_Jet55_chMaxJtpt0p05->Write();

  hpbpb_Jet55_eMaxJtpt0p1->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p1->Write();
  hpbpb_Jet55_eMaxJtpt0p2->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p2->Write();
  hpbpb_Jet55_eMaxJtpt0p3->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p3->Write();
  hpbpb_Jet55_eMaxJtpt0p4->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p4->Write();
  hpbpb_Jet55_eMaxJtpt0p5->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p5->Write();
  hpbpb_Jet55_eMaxJtpt0p6->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p6->Write();
  hpbpb_Jet55_eMaxJtpt0p7->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p7->Write();
  hpbpb_Jet55_eMaxJtpt0p8->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p8->Write();
  hpbpb_Jet55_eMaxJtpt0p9->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p9->Write();

  hpbpb_Jet55_eMaxSumcand0p1->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p1->Write();
  hpbpb_Jet55_eMaxSumcand0p2->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p2->Write();
  hpbpb_Jet55_eMaxSumcand0p3->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p3->Write();
  hpbpb_Jet55_eMaxSumcand0p4->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p4->Write();
  hpbpb_Jet55_eMaxSumcand0p5->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p5->Write();
  hpbpb_Jet55_eMaxSumcand0p6->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p6->Write();
  hpbpb_Jet55_eMaxSumcand0p7->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p7->Write();
  hpbpb_Jet55_eMaxSumcand0p8->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p8->Write();
  hpbpb_Jet55_eMaxSumcand0p9->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxSumcand0p9->Write();

  hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Write();
  hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Write();
  hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Write();
  hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Write();
  hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Write();
  hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Write();
  hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Write();
  hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Write();
  hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hpbpb_Jet55);
  hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Write();

  #endif
  hpbpb_chMaxJtpt_jtpt->Write();
  hpbpb_eMaxJtpt_jtpt->Write();
  hpbpb_eMaxJtpt_chMaxJtpt->Write();
  hpbpb_eMaxSumcand_jtpt->Write();

  for(int a = 0;a<ptSelection;++a){
    for(int t = 0;t<trigValue-1;++t){
      hpbpb_eMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_phMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_neMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_muMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_eMaxSumcand_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_eMaxSumcand_jtpt_ptselection[t][a]->Write();
      hpbpb_chMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_phMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_neMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_muMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_eMaxJtpt_jtpt_ptselection[t][a]->Write();
      //hpbpb_eMaxSumcand_chMaxJtpt_refptselection[t][a]->Write();
      hpbpb_eMaxSumcand_jtpt_refptselection[t][a]->Write();
      hpbpb_eMaxSumcand_refpt_refptselection[t][a]->Write();
      hpbpb_eMaxSumcand_eMaxJtpt_refANDjtptselection[t][a]->Write();
      hpbpb_eMaxSumcand_eMaxJtpt_ptselection[t][a]->Write();
      hpbpb_eMaxSumcand_chMaxJtpt_refANDjtptselection[t][a]->Write();
      hpbpb_eMaxSumcand_jtpt_refANDjtptselection[t][a]->Write();
    }
  }

  for(int k = 0;k<no_radius;k++){
    for(int j = 0;j<nbins_eta;j++){
      for(int i = 0;i<nbins_cent;i++){
	hpbpb_chMax[k][j][i]->Write();
	hpbpb_phMax[k][j][i]->Write();
	hpbpb_neMax[k][j][i]->Write();
	hpbpb_muMax[k][j][i]->Write();
	hpbpb_eMax[k][j][i]->Write();
	hpbpb_chSum[k][j][i]->Write();
	hpbpb_phSum[k][j][i]->Write();
	hpbpb_neSum[k][j][i]->Write();
	hpbpb_muSum[k][j][i]->Write();
	hpbpb_eSum[k][j][i]->Write();

	hpbpb_chMax_withCut[k][j][i]->Write();
	hpbpb_phMax_withCut[k][j][i]->Write();
	hpbpb_neMax_withCut[k][j][i]->Write();
	hpbpb_muMax_withCut[k][j][i]->Write();
	hpbpb_eMax_withCut[k][j][i]->Write();
	hpbpb_chSum_withCut[k][j][i]->Write();
	hpbpb_phSum_withCut[k][j][i]->Write();
	hpbpb_neSum_withCut[k][j][i]->Write();
	hpbpb_muSum_withCut[k][j][i]->Write();
	hpbpb_eSum_withCut[k][j][i]->Write();

      }
    }
  }
  
#if 0

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      for(int i = 0;i<nbins_cent;i++){

	for(int t = 0;t<trigValue;t++){

	  hJets_noTrkpTCut[k][j][i][t]->Write();
	  hJets_3TrkpTCut[k][j][i][t]->Write();
	  hJets_5TrkpTCut[k][j][i][t]->Write();
	  hJets_7TrkpTCut[k][j][i][t]->Write();
	  hJets_10TrkpTCut[k][j][i][t]->Write();

	}

      }

    }
  }
#endif

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      for(int i = 0;i<nbins_cent;i++){

	hpbpb_JetComb_gen[k][j][i]->Write();
	hpbpb_Jet80_gen[k][j][i]->Write();
	hpbpb_Jet65_gen[k][j][i]->Write();
	hpbpb_Jet55_gen[k][j][i]->Write();
	hpbpb_JetComb_reco[k][j][i]->Write();
	hpbpb_Jet80_reco[k][j][i]->Write();
	hpbpb_Jet65_reco[k][j][i]->Write();
	hpbpb_Jet55_reco[k][j][i]->Write();
	
	hpbpb_mcclosure_JetComb_gen[k][j][i]->Write();
	hpbpb_mcclosure_Jet80_gen[k][j][i]->Write();
	hpbpb_mcclosure_Jet65_gen[k][j][i]->Write();
	hpbpb_mcclosure_Jet55_gen[k][j][i]->Write();
	hpbpb_mcclosure_JetComb_data[k][j][i]->Write();
	hpbpb_mcclosure_Jet80_data[k][j][i]->Write();
	hpbpb_mcclosure_Jet65_data[k][j][i]->Write();
	hpbpb_mcclosure_Jet55_data[k][j][i]->Write();
	
	hpbpb_matrix_HLT[k][j][i]->Write();
	hpbpb_matrix[k][j][i]->Write();
	hpbpb_gen[k][j][i]->Write();
	hpbpb_reco[k][j][i]->Write();
	hpbpb_mcclosure_data[k][j][i]->Write();
	hpbpb_mcclosure_gen[k][j][i]->Write();
	hpbpb_mcclosure_matrix[k][j][i]->Write();
	
      }
      
    }
  }

  f.Write();
  f.Close();

  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
