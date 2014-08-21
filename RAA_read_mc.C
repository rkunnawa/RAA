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
static const int nbins_eta = 10;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0},
  {-2.0,+2.0},
  {-2.5,-2.0},
  {-2.0,-1.5},
  {-1.5,-1.0},
  {-1.0,-0.5},
  {-0.5,+0.5},
  {+0.5,+1.0},
  {+1.0,+1.5},
  {+1.5,+2.0}
};

static const char etaWidth [nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
  "n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
  "p10_eta_p15","p15_eta_p20"
};

*/

static const int nbins_eta = 2;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0},
  {-2.0,+2.0}
};

static const char etaWidth[nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20"
};

static const int no_radius = 2;//testing purposes 
static const int list_radius[no_radius] = {3,4};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 7; 
//static const int list_radius[no_radius] = {1,2,3,4,5,6,7};

static const int nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 
static const Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };

//static const int nAlgos = 9;
//static const int BinLabelN = 11;
//remember to change this to run akPu3PF for pPb and akVs3PF for Pbpb datasets. or just create a separate header file which will be way easier. 
//static const char *algoName[nAlgos] = { "", "icPu5", "akPu2PF", "akPu3PF", "akPu4PF", "akPu5PF" , "akPu2Calo", "akPu3Calo", "akPu4Calo" };
//static const char *algoNamePP[nAlgos] = { "", "icPu5", "ak2PF", "ak3PF", "ak4PF", "ak5PF" , "ak2Calo", "ak3Calo", "ak4Calo" };
//static const char *algoNameGen[nAlgos] = { "", "icPu5", "akPu2PF", "akVs3PF", "akPu4PF", "akPu2PF", "akPu3PF", "akPu4PF" };
//static const char *BinLabel[BinLabelN] = {"100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180", "180-200", "200-240","240-300" };


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
    tJet = (TTree*)tFile->Get(jetTree);
    tJet->SetBranchAddress("jtpt" , jtpt );
    tJet->SetBranchAddress("rawpt", rawpt);
    tJet->SetBranchAddress("trackMax" , trackMax );
    tJet->SetBranchAddress("chargedMax",chargedMax);
    tJet->SetBranchAddress("chargedSum",chargedSum);
    tJet->SetBranchAddress("neutralMax",neutralMax);
    tJet->SetBranchAddress("neutralSum",neutralSum);
    tJet->SetBranchAddress("refpt", refpt);
    tJet->SetBranchAddress("nref" ,&njets);
    tJet->SetBranchAddress("jteta", jteta);
    tJet->SetBranchAddress("jtm",jtmass);
    tJet->SetBranchAddress("pthat",&pthat);
    if(isPbPb) tJet->SetBranchAddress("subid",&subid);
    if (loadGenJet) tGenJet = (TTree*)tFile->Get(genJetTree);
    if (loadGenJet) tGenJet->SetBranchAddress("ngen" ,&ngen);
    if (loadGenJet) tGenJet->SetBranchAddress("genpt", genpt);
    if (loadGenJet) tGenJet->SetBranchAddress("gensubid", gensubid);
    tEvt->SetBranchAddress("hiBin",&bin);
    tEvt->SetBranchAddress("vz",&vz);
    tSkim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    if(isPbPb) tSkim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    else tSkim->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
    tJet->AddFriend(tEvt);
    tJet->AddFriend(tSkim);
  };
  TFile *tFile;
  TTree *tJet;
  TTree *tGenJet;
  TTree *tEvt;
  TTree* tSkim;
  float jtpt[1000];
  float rawpt[1000];
  float refpt[1000];
  float jteta[1000];
  float jtmass[1000];
  float trackMax[1000];
  float chargedMax[1000];
  float neutralMax[1000];
  float chargedSum[1000];
  float neutralSum[1000];
  float genpt[1000];
  int gensubid[1000];
  float subid[1000];
  float vz;
  float pthat;
  int njets;
  int ngen;
  int bin;     
  int pHBHENoiseFilter;
  int pPAcollisionEventSelectionPA;
  int pcollisionEventSelection;
};


using namespace std;

void RAA_read_mc(char *algo = "Vs", char *jet_type = "Calo"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  cout<<"Running for Algorithm "<<algo<<" "<<jet_type<<endl;
 
  bool printDebug = true;

  const int nbins_pthat = 9;
  Double_t boundaries_pthat[nbins_pthat+1];
  char *fileName_pthat[nbins_pthat+1];
  Double_t xsection[nbins_pthat+1];
  Double_t entries[nbins_pthat]; 
  //there are two ways in which we can select the no of events we use to scale - it has to be between the pthat range. 
  //first file name - partial 50K statistics, second one is full statistics sample. 
  //similarly the entries number is for the small statistics. 

  boundaries_pthat[0]=15;
  //fileName_pthat[0] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15/hiForest_DijetpT15_Hydjet1p8_STARTHI53_LV1_Track9_Jet30_v15.root";
  fileName_pthat[0] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT15_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[0]= 2.034e-01;
  //entries[0] = ;//total - 48588
  
  boundaries_pthat[1]=30;
  fileName_pthat[1] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT30_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[1]= 1.075e-02;
  //entries[1] = ;//total - 48428
  
  boundaries_pthat[2]=50;
  fileName_pthat[2] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT50_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[2]= 1.025e-03;
  //entries[2] = ;//total - 50000
  
  boundaries_pthat[3]=80;
  fileName_pthat[3] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT80_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[3]= 9.865e-05;
  //entries[3] = ;//total - 49500
  
  boundaries_pthat[4]=120;
  fileName_pthat[4] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT120_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[4]= 1.129e-05;
  //entries[4] = ;//total - 49500

  boundaries_pthat[5]=170;
  fileName_pthat[5] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT170_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[5]= 1.465e-06;
  //entries[5] = ;//total - 49444

  boundaries_pthat[6]=220;
  fileName_pthat[6] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT220_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[6]= 2.837e-07;
  //entries[6] = ;//total - 49460

  boundaries_pthat[7]=280;
  fileName_pthat[7] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT280_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[7]= 5.323e-08;
  //entries[7] = ;//total - 49541

  boundaries_pthat[8]=370;
  fileName_pthat[8] = "/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/hiForest_DijetpT370_Hydjet1p8_STARTHI53_LV1_v15_full.root";
  xsection[8]= 5.934e-09;
  //entries[8] = ;//total - 19031

  boundaries_pthat[9] = 2000;
  xsection[9] = 0.0;

  // Vertex & centrality reweighting for PbPb
  TF1 *fVz;
  TF1* fCentralityWeight;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);
  fCentralityWeight = new TF1("fCentralityWeight","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
  fCentralityWeight->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04);

  /* 
     const int nbinsPP_pthat = 11;
     Double_t boundariesPP_pthat[nbinsPP_pthat+1];
     char *fileNamePP_pthat[nbinsPP_pthat+1];
     Double_t xsectionPP[nbinsPP_pthat+1];
  
     boundariesPP_pthat[0]=15;
     fileNamePP_pthat[0]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt15/HiForest_v81_merged01/pt15_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[0]= 1.079e-02;
     xsectionPP[0]= 2.034e-01;
  
     boundariesPP_pthat[1]=30;
     fileNamePP_pthat[1]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt30/HiForest_v81_merged01/pt30_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[1]= 1.021e-03;
     xsectionPP[1]= 1.075e-02;
  
     boundariesPP_pthat[2]=50;
     fileNamePP_pthat[2]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt50/HiForest_v81_merged01/pt50_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[2]= 9.913e-05;
     xsectionPP[2]= 1.025e-03;
  
     boundariesPP_pthat[3]=80;
     fileNamePP_pthat[3]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt80/HiForest_v81_merged01/pt80_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[3]= 1.128e-05;
     xsectionPP[3]= 9.865e-05;
  
     boundariesPP_pthat[4]=120;
     fileNamePP_pthat[4]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt120/HiForest_v81_merged01/pt120_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[4]= 1.470e-06;
     xsectionPP[4]= 1.129e-05;
  
     boundariesPP_pthat[5]=170;
     fileNamePP_pthat[5]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt170/HiForest_v81_merged01/pt170_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[5]= 5.310e-07;
     xsectionPP[5]= 1.465e-06;
  
     boundariesPP_pthat[6]=220;
     fileNamePP_pthat[6]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt220/HiForest_v81_merged01/pt220_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[6]= 1.192e-07;	
     xsectionPP[6]= 2.837e-07;
  
     boundariesPP_pthat[7]=280;
     fileNamePP_pthat[7]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt280/HiForest_v81_merged01/pt280_pp2013_P01_prod22_v81_merged_forest_0.root";
     //xsectionPP[7]= 3.176e-08;
     xsectionPP[7]= 5.323e-08;
 
     xsectionPP[8] = 0;
     boundariesPP_pthat[8]=1000;
  */

  const int nbinsPP_pthat = 11;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  char *fileNamePP_pthat[nbinsPP_pthat+1];
  Double_t xsectionPP[nbinsPP_pthat+1];
  
  boundariesPP_pthat[0]=15;
  fileNamePP_pthat[0] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat15_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[0]= 0.2034;
  //  entries[0] = 71680;  
  
  boundariesPP_pthat[1]=30;
  fileNamePP_pthat[1] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat30_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[1]= 0.01075;
  //entries[1] = 52160;
  
  boundariesPP_pthat[2]=50;
  fileNamePP_pthat[2] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat50_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[2]= 0.001025;
  // entries[2] = 50240;
  
  boundariesPP_pthat[3]=80;
  fileNamePP_pthat[3] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat80_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[3]= 9.8650e-05;
  // entries[3] = 52160;
  
  boundariesPP_pthat[4]=120;
  fileNamePP_pthat[4] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat120_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[4]= 1.1290e-05;
  // entries[4] = 53760;

  boundariesPP_pthat[5] = 170;
  fileNamePP_pthat[5] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat170_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[5]= 1.4650e-06;
  //entries[5] = 53120;
  
  boundariesPP_pthat[6]=220;
  fileNamePP_pthat[6] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat220_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[6]= 2.8370e-07;
  // entries[6] = 54080;
  
  boundariesPP_pthat[7]=280;
  fileNamePP_pthat[7] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat280_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[7]= 5.3230e-08;
  // entries[7] = 53120;
  
  boundariesPP_pthat[8]=370;
  fileNamePP_pthat[8] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat370_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[8]= 5.9340e-09;
  //entries[8] = 52800;
  
  boundariesPP_pthat[9]=460;
  fileNamePP_pthat[9] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat460_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[9]= 8.1250e-10;
  //entries[9] = 54080;
  
  boundariesPP_pthat[10]=540;
  fileNamePP_pthat[10] = "/mnt/hadoop/cms/store/user/rkunnawa/53X_Production/pp_official_MC_merged_files/HiForest_pp_official_MC_pthat540_53X_STARTHI53_V28_5_3_16_trk8_Jet28_merged.root";
  xsectionPP[10]= 1.4670e-10;
  //entries[10] = 53440;  
  
  xsectionPP[11] = 0;
  boundariesPP_pthat[11]=2000; 
  
  // lets declare all the histograms here. 

  TH1F *hpbpb_gen[no_radius][nbins_eta][nbins_cent+1],*hpbpb_reco[no_radius][nbins_eta][nbins_cent+1];
  TH2F *hpbpb_matrix[no_radius][nbins_eta][nbins_cent+1];
  //TH2F *hpbpb_response[nbins_cent+1];
  TH1F *hpbpb_mcclosure_data[no_radius][nbins_eta][nbins_cent+1];

  TH1F *hpp_gen[no_radius][nbins_eta];
  TH1F *hpp_reco[no_radius][nbins_eta];
  TH2F *hpp_matrix[no_radius][nbins_eta];
  TH1F *hpp_mcclosure_data[no_radius][nbins_eta];

  TH1F *hCentMC[no_radius];
  
  TH1F *hVzMC[no_radius];
  TH1F *hVzPPMC[no_radius];

  TH1F *hPtHat[no_radius];
  TH1F *hPtHatRaw[no_radius];
  TH1F *hPtHatPP[no_radius];
  TH1F *hPtHatRawPP[no_radius];

  TH1F *hPbPb_pthat_fine[no_radius];
  TH1F *hPP_pthat_fine[no_radius];

  TH1F *hPbPb_pthat_fine_noScale[no_radius];
  TH1F *hPP_pthat_fine_noScale[no_radius];
  

  for(int k = 0;k<no_radius;k++){
    //cout<<"radius = "<<list_radius[k]<<endl;
    for(int j = 0;j<nbins_eta;j++){
      //cout<<"eta bin = "<<j<<endl;
      for(int i = 0;i<nbins_cent;i++){
	//cout<<"cent bin = "<<i<<endl;
        hpbpb_gen[k][j][i] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Gen refpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
	//cout<<"A"<<endl;
	hpbpb_reco[k][j][i] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Reco jtpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
	//cout<<"B"<<endl;
	hpbpb_matrix[k][j][i] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Matrix refpt jtpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
	//cout<<"C"<<endl;
	hpbpb_mcclosure_data[k][j][i] = new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("data for unfolding mc closure test R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
	//cout<<"D"<<endl;
	//hpbpb_response[h] = new TH2F(Form("hpbpb_response_cent%d",i),Form("response jtpt refpt %2.0f - %2.0f cent",5*boundaries_cent[h],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
      }// centrality bin loop
      
      hpbpb_gen[k][j][nbins_cent] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Gen refpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpbpb_reco[k][j][nbins_cent] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Reco jtpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpbpb_matrix[k][j][nbins_cent] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Matrix refpt jtpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
      hpbpb_mcclosure_data[k][j][nbins_cent] = new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("data for unfolding mc closure test R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      //hpbpb_response[nbins_cent] = new TH2F(Form("hpbpb_response_cent%d",nbins_cent),"response jtpt refpt 0-200 cent",1000,0,1000,1000,0,1000);

      hpp_gen[k][j] = new TH1F(Form("hpp_gen_R%d_%s",list_radius[k],etaWidth[j]),Form("gen refpt R%d %s",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpp_reco[k][j] = new TH1F(Form("hpp_reco_R%d_%s",list_radius[k],etaWidth[j]),Form("reco jtpt R%d %s",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpp_matrix[k][j] = new TH2F(Form("hpp_matrix_R%d_%s",list_radius[k],etaWidth[j]),Form("matrix refpt jtpt R%d %s",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
      //TH2F* hpp_response = new TH2F("hpp_response","response jtpt refpt",1000,0,1000,1000,0,1000);
      hpp_mcclosure_data[k][j] = new TH1F(Form("hpp_mcclosure_data_R%d_%s",list_radius[k],etaWidth[j]),Form("data for unfolding mc closure test pp R%d %s",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      
    }// eta bin loop
    
    hVzMC[k] = new TH1F(Form("hVzMC_R%d",list_radius[k]),Form("PbPb MC Vz R%d",list_radius[k]),60,-15,+15);
    hVzPPMC[k] = new TH1F(Form("hVzPPMC_R%d",list_radius[k]),Form("PP MC Vz R%d",list_radius[k]),60,-15,+15);
    
    hCentMC[k] = new TH1F(Form("hCentMC_R%d",list_radius[k]),"",nbins_cent,boundaries_cent);
    
    hPtHat[k] = new TH1F(Form("hPtHat_R%d",list_radius[k]),"",nbins_pthat,boundaries_pthat);
    hPtHatRaw[k] = new TH1F(Form("hPtHatRaw_R%d",list_radius[k]),"",nbins_pthat,boundaries_pthat);
    hPtHatPP[k] = new TH1F(Form("hPtHatPP_R%d",list_radius[k]),"",nbinsPP_pthat,boundariesPP_pthat);
    hPtHatRawPP[k] = new TH1F(Form("hPtHatRawPP_R%d",list_radius[k]),"",nbinsPP_pthat,boundariesPP_pthat);

    hPbPb_pthat_fine[k] = new TH1F(Form("hPbPb_pthat_fine_R%d",list_radius[k]),Form("PbPb pthat distribution for R=0.%d",list_radius[k]),700,0,700);
    hPP_pthat_fine[k] = new TH1F(Form("hPP_pthat_fine_R%d",list_radius[k]),Form("pp pthat distribution for R=0.%d",list_radius[k]),1000,0,1000);

    hPbPb_pthat_fine_noScale[k] = new TH1F(Form("hPbPb_pthat_fine_noScale_R%d",list_radius[k]),Form("PbPb pthat distribution (unscaled) for R=0.%d",list_radius[k]),700,0,700);
    hPP_pthat_fine_noScale[k] = new TH1F(Form("hPP_pthat_fine_noScale_R%d",list_radius[k]),Form("PP pthat distribution (unscaled) for R=0.%d",list_radius[k]),1000,0,1000);

  }// radii loop

  // Setup jet data branches - this will be 2D with [radius][pthat-file], but the histogram here is just 1D with [radius]
  JetData *data[no_radius][nbins_pthat]; 
  JetData *dataPP[no_radius][nbinsPP_pthat];
  for(int k = 0;k<no_radius;k++){
    //cout<<"Radius = "<<list_radius[k]<<endl;
    //cout<<"reading all the pbpb mc files"<<endl;
    for (int h=0;h<nbins_pthat;h++) {
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
    cout<<"reading all the pp mc files"<<endl;
    for (int h=0;h<nbinsPP_pthat;h++){ 
      dataPP[k][h] = new JetData(fileNamePP_pthat[h],Form("ak%d%sJetAnalyzer/t",list_radius[k],jet_type),Form("ak%d%sJetAnalyzer/t",list_radius[k],jet_type),0,0);
      TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbinsPP_pthat,boundariesPP_pthat);
      dataPP[k][h]->tJet->Project("hPtHatTmp","pthat");
      hPtHatRawPP[k]->Add(hPtHatTmp);
      delete hPtHatTmp;
    }//pthatpp loop
  }//radius loop

  // checking the histograms to see if something is filled. 
  hPtHatRaw[1]->Print("base");
  hPtHatRawPP[1]->Print("base");
  
  for(int k = 0;k<no_radius;k++){
    cout<<"Filling MC for radius = "<<list_radius[k]<<endl;
    // fill PbPb MC 
    cout<<"Filling PbPb MC"<<endl;
    for (int h=0;h<nbins_pthat;h++) {
      if (xsection[h]==0) continue;
      cout <<"Loading pthat"<<boundaries_pthat[h]<<" sample, cross section = "<<xsection[h]<< Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[h],boundaries_pthat[h+1])<<endl;
      cout<<data[k][h]->tJet->GetEntries()<<endl;
      for (Long64_t jentry=0; jentry<data[k][h]->tJet->GetEntries();jentry++) {
	//for (Long64_t jentry=0; jentry<100;jentry++) {

        //cout<<"hi"<<endl;
        data[k][h]->tEvt->GetEntry(jentry);
        data[k][h]->tJet->GetEntry(jentry);
        //data[k][h]->tGenJet->GetEntry(jentry);
        //if(data[k][h]->pthat<boundaries_pthat[h] || data[k][h]->pthat>boundaries_pthat[h+1]) continue;
        //remember this cut is there because there was some rediculous values of pthats of -1 in the private production forests
        //if(jentry%100==0)cout<<"pthat of that event = "<<data[k][h]->pthat<<endl;
      
        int pthatBin = hPtHat[k]->FindBin(data[k][h]->pthat);
	//cout<<"pthat = "<<data[k][h]->pthat<<", pthatBin = "<<pthatBin<<", boundaries_pthat[pthatBin] = "<<boundaries_pthat[pthatBin]<<endl;
        //if(jentry2%100==0)cout<<"pthatBin = "<<pthatBin<<endl;
      
        //cout<<xsection[pthatBin-1]-xsection[pthatBin]<<endl;
        //cout<<"nentries = "<<hPtHatRaw->GetBinContent(pthatBin)<<endl;
        double scale_old = (double)(xsection[pthatBin-1]-xsection[pthatBin])/hPtHatRaw[k]->GetBinContent(pthatBin);
	//cout<<"scale = "<<scale<<endl;

	//from Pawan's code: /net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/combinePtHatBins/pbpbJEC2014/condor/CondorPbPbCalJec.C
	TEventList *el = new TEventList("el","el");
	double pthat_event = data[k][h]->pthat;
	double pthat_lower = boundaries_pthat[h];
	double pthat_upper = boundaries_pthat[h+1];
	stringstream selection; selection<<"pthat_lower<"<<pthat_upper;

	data[k][h]->tJet->Draw(">>el",selection.str().c_str());
	double fentries = el->GetN();
	if(jentry==0)cout<<"tree entries: "<<data[k][h]->tJet->GetEntries()<<" elist: "<<fentries<<endl;
	delete el;

	//double fentries = data[k][h]->tJet->GetEntries(data[k][h]->pthat>=boundaries_pthat[h] && data[k][h]->pthat<boundaries_pthat[h+1]);
	//if(jentry==0)cout<<fentries<<endl;
	double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/fentries;

	//cout<<"xsection[pthatBin-1] = "<<xsection[pthatBin-1]<<", xsection[pthatBin] = "<<xsection[pthatBin]<<", bin content = "<<hPtHatRaw[k]->GetBinContent(pthatBin)<<endl;
        //double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/entries[h];
	
        if(fabs(data[k][h]->vz)>15) continue;
        int cBin = hCentMC[k]->FindBin(data[k][h]->bin)-1;
        //int cBin = nbins_cent-1;
        double weight_cent=1;
        double weight_pt=1;
        double weight_vz=1;
	
        //weight_cent = fCentralityWeight->Eval(data[h]->bin);
        //weight_vz = fVz->Eval(data[k][h]->vz);
	
        if(scale*weight_cent*weight_vz <=0 ) {
	  cout<<"RED FLAG RED FLAG RED FLAG"<<endl;
	  continue;
	}
	
	hPbPb_pthat_fine[k]->Fill(data[k][h]->pthat,weight_vz*scale);
	hPbPb_pthat_fine_noScale[k]->Fill(data[k][h]->pthat);
        hCentMC[k]->Fill(data[k][h]->bin,scale*weight_cent*weight_vz);
        hVzMC[k]->Fill(data[k][h]->vz,scale*weight_cent*weight_vz);
        if (cBin>=nbins_cent) continue;
        if (cBin==-1) continue;
        hPtHat[k]->Fill(data[k][h]->pthat,scale*weight_cent*weight_vz);
	
        //cout<<"scale = "<<scale<<endl;
	
        /*
	  int hasLeadingJet = 0;
	  for (int k= 0; k < data[h]->njets; k++) { 
	  if ( data[h]->jteta[k]  > 2. || data[h]->jteta[k] < -2. ) continue;
	  if ( data[h]->jtpt[k]>100) {
	  hasLeadingJet = 1;
	  }
	  break;
				 
	  }
	  if (hasLeadingJet == 0) continue;
        */

        for (int g = 0; g < data[k][h]->njets; g++) {
  
          for(int j = 0;j<nbins_eta;j++){

            //int subEvt=-1;
	    if ( data[k][h]->subid[g] != 0 ) continue;
            if ( data[k][h]->rawpt[g]  <= 10. ) continue;
	    if ( data[k][h]->jtpt[g] > 2.*data[k][h]->pthat) continue;
            if ( data[k][h]->jteta[g]  > boundaries_eta[j][1] || data[k][h]->jteta[g] < boundaries_eta[j][0] ) continue;
	    
            // jet quality cuts here: 
            if ( data[k][h]->chargedMax[g]/data[k][h]->jtpt[g]<0.01) continue;
	    //if ( data[k][h]->neutralMax[g]/TMath::Max(data[h]->chargedSum[k],data[h]->neutralSum[k]) < 0.975)continue;

	    //for (int l= 0; l< data[h]->ngen;l++) {
	    //  if (data[h]->refpt[k]==data[h]->genpt[l]) {
	    //    subEvt = data[h]->gensubid[l];
	    //    break;
	    //  } 
	    //}
	    //if (subEvt!=0) continue;
	    //if (uhist[cBin]->hMeasMatch!=0) {
	    //   int ptBinNumber = uhist[cBin]->hMeasMatch->FindBin(data[h]->jtpt[k]);
	    //   int ratio = uhist[cBin]->hMeasMatch->GetBinContent(ptBinNumber);
	    //if (ratio!=0) weight_pt = 1./ratio;
	    //}
	    //if (!isMC||jentry2<data[h]->tJet->GetEntries()/2.) {
	    //cout<<"going to fill the histograms now"<<endl;
	    //cout<<"fvz = "<<weight_vz<<endl;
	    
	    //hpbpb_response[cBin]->Fill(data[h]->jtpt[k],data[h]->refpt[k],scale*weight_vz);
	    hpbpb_matrix[k][j][cBin]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz);
	    hpbpb_gen[k][j][cBin]->Fill(data[k][h]->refpt[g],scale*weight_vz);
	    hpbpb_reco[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz);
	    
	    //hpbpb_response[nbins_cent]->Fill(data[h]->jtpt[k],data[h]->refpt[k],scale*weight_vz);
	    hpbpb_matrix[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],data[k][h]->jtpt[g],scale*weight_vz);
	    hpbpb_gen[k][j][nbins_cent]->Fill(data[k][h]->refpt[g],scale*weight_vz);
	    hpbpb_reco[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz);
	    
	    if (jentry>data[k][h]->tJet->GetEntries()/2.) {
	      hpbpb_mcclosure_data[k][j][cBin]->Fill(data[k][h]->jtpt[g],scale*weight_vz);
	      hpbpb_mcclosure_data[k][j][nbins_cent]->Fill(data[k][h]->jtpt[g],scale*weight_vz);
	    }
	 
	    //uhist[cBin]-> hMeasJECSys->Fill(data[h]->jtpt[k]*(1.+0.02/nbins_cent*(nbins_cent-i)),scale*weight_cent*weight_pt*weight_vz); 
	
	  }// eta bins loop
	      
        }//njets loop
      
      }//nentry loop
    
    }//ptbins loop
  
  
    // Vertex reweighting for pp
    TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
    fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
  
    cout<<"Filling PP MC"<<endl;
    // fill pp MC
    for (int h=0;h<nbinsPP_pthat;h++) {
      if (xsectionPP[h]==0) continue;
      //float scale=(xsectionPP[h]-xsectionPP[i+1])/dataPP[k][h]->tJet->GetEntries(Form("pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[h],boundariesPP_pthat[i+1])); 
      cout <<"Loading PP pthat"<<boundariesPP_pthat[h]<<" sample, cross section = "<<xsectionPP[h]<< Form(" pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[h],boundariesPP_pthat[h+1])<<endl;
      //cout<<""<<endl;
      for (Long64_t jentry=0; jentry<dataPP[k][h]->tJet->GetEntries();jentry++) {
	//for (Long64_t jentry=0; jentry<10;jentry++) {
        dataPP[k][h]->tEvt->GetEntry(jentry);
	dataPP[k][h]->tJet->GetEntry(jentry);

	//dataPP[k][h]->tGenJet->GetEntry(jentry);
	//if(dataPP[k][h]->pthat<boundariesPP_pthat[h] || dataPP[k][h]->pthat>boundariesPP_pthat[i+1]) continue;
        //if(dataPP[k][h]->bin<=28) continue;
        int pthatBin = hPtHatPP[k]->FindBin(dataPP[k][h]->pthat);
        float scalepp = (xsectionPP[pthatBin-1]-xsectionPP[pthatBin])/hPtHatRawPP[k]->GetBinContent(pthatBin);
        if(fabs(dataPP[k][h]->vz)>15) continue;
        double weight_cent=1;
        double weight_pt=1;
        double weight_vz=1;
        //if(!dataPP[k][h]->pPAcollisionEventSelectionPA || !dataPP[k][h]->pHBHENoiseFilter) continue;
	//for now the MC doesnt have pPAcollisionEventSelectionPA so dont search for it. 
	if(!dataPP[k][h]->pHBHENoiseFilter) continue;

        weight_vz = fVzPP->Eval(dataPP[k][h]->vz);
        //if (weight_vz>5||weight_vz<0.5) cout <<dataPP[k][h]->vz<<" "<<weight_vz<<endl;
        //weight_vz = 1;
	hPP_pthat_fine[k]->Fill(dataPP[k][h]->pthat,scalepp*weight_vz);
	hPP_pthat_fine_noScale[k]->Fill(dataPP[k][h]->pthat);
        hPtHatPP[k]->Fill(dataPP[k][h]->pthat,scalepp*weight_vz);
        int hasLeadingJet = 0;
        hVzPPMC[k]->Fill(dataPP[k][h]->vz,scalepp*weight_vz);
        /*
	  for (int k= 0; k < dataPP[k][h]->njets; k++) { 
	  if ( dataPP[k][h]->jteta[k]  > 2. || dataPP[k][h]->jteta[k] < -2. ) continue;
	  if ( dataPP[k][h]->jtpt[k]>100) {
	  hasLeadingJet = 1;
	  }
	  break;
	
	  }
	  if (hasLeadingJet == 0) continue;
        */

        for (int g= 0; g< dataPP[k][h]->njets; g++) { 

          for(int j = 0;j<nbins_eta;j++){
            
            int subEvt=-1;
            if ( dataPP[k][h]->rawpt[g]  <= 10. ) continue;
            if ( dataPP[k][h]->jteta[g]  > boundaries_eta[j][1] || dataPP[k][h]->jteta[g] < boundaries_eta[j][0] ) continue;

            // jet QA cuts: 
            if ( dataPP[k][h]->chargedMax[g]/dataPP[k][h]->jtpt[g]<0.01) continue;
            //if ( dataPP[k][h]->neutralMax[g]/TMath::Max(dataPP[k][h]->chargedSum[g],dataPP[k][h]->neutralSum[g]) < 0.975)continue;
            //if ( dataPP[k][h]->neu)

            //if (uhist[nbins_cent]->hMeasMatch!=0) {
            //   int ptBinNumber = uhist[nbins_cent]->hMeasMatch->FindBin(dataPP[k][h]->jtpt[k]);
            //   int ratio = uhist[nbins_cent]->hMeasMatch->GetBinContent(ptBinNumber);
            //if (ratio!=0) weight_pt = 1./ratio;
            //}
          
            //if (!isMC||jentry<dataPP[k][h]->tJet->GetEntries()/2.) {
          
            //hpp_response->Fill(dataPP[k][h]->jtpt[k],dataPP[k][h]->refpt[k],scalepp*weight_vz);
            hpp_matrix[k][j]->Fill(dataPP[k][h]->refpt[g],dataPP[k][h]->jtpt[g],scalepp*weight_vz);
            hpp_gen[k][j]->Fill(dataPP[k][h]->refpt[g],scalepp*weight_vz);   
            hpp_reco[k][j]->Fill(dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	    
	    
            if (jentry>dataPP[k][h]->tJet->GetEntries()/2.)
              hpp_mcclosure_data[k][j]->Fill(dataPP[k][h]->jtpt[g],scalepp*weight_vz);
            
          }//eta loop
	              
        }//njet loop     
      
      }//nentry loop
    
    }//ptbins loop

  }// radius loop

  TDatime date;

  //declare the output file 
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_mc_ak%s%s_%d.root",algo,jet_type,date.GetDate()),"RECREATE");
  f.cd();

  for(int k = 0;k<no_radius;k++){
    
    for(int j=0;j<nbins_eta;j++){
      
      for(int i = 0;i<=nbins_cent;i++){
    
        divideBinWidth(hpbpb_gen[k][j][i]);
        divideBinWidth(hpbpb_reco[k][j][i]);
        divideBinWidth(hpbpb_mcclosure_data[k][j][i]);
        hpbpb_gen[k][j][i]->Write();
        hpbpb_gen[k][j][i]->Print("base");
        hpbpb_reco[k][j][i]->Write();
        hpbpb_reco[k][j][i]->Print("base");
        hpbpb_matrix[k][j][i]->Write();
        hpbpb_matrix[k][j][i]->Print("base");
        hpbpb_mcclosure_data[k][j][i]->Write();
        hpbpb_mcclosure_data[k][j][i]->Print("base");
    
      }// cent loop 
      
      divideBinWidth(hpp_gen[k][j]);
      divideBinWidth(hpp_reco[k][j]);
      divideBinWidth(hpp_mcclosure_data[k][j]);

      hpp_gen[k][j]->Write();
      hpp_gen[k][j]->Print("base");
      hpp_reco[k][j]->Write();
      hpp_reco[k][j]->Print("base");
      hpp_matrix[k][j]->Write();
      hpp_matrix[k][j]->Print("base");
      hpp_mcclosure_data[k][j]->Write();
      hpp_mcclosure_data[k][j]->Print("base");

    }//eta loop
    //just the check the Pthat distributions for PbPb and pp. should be fine. 

    hCentMC[k]->Print("base");
    hCentMC[k]->Write();
    hPtHat[k]->Print("base");
    hPtHat[k]->Write();
    hPtHatPP[k]->Print("base");
    hPtHatPP[k]->Write();
    hPbPb_pthat_fine[k]->Print("base");
    hPbPb_pthat_fine[k]->Write();
    hPP_pthat_fine[k]->Print("base");
    hPP_pthat_fine[k]->Write();

    hPbPb_pthat_fine_noScale[k]->Print("base");
    hPbPb_pthat_fine_noScale[k]->Write();
    hPP_pthat_fine_noScale[k]->Print("base");
    hPP_pthat_fine_noScale[k]->Write();

  }// radius loop

  
  f.Write();
  f.Close();
  
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
