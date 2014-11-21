// Raghav Kunnawalkam Elayavalli
// Nov 20th 2014
// Rutgers
// for questions or comments: raghav.k.e at CERN dot CH


//
// Read in the MC files to figure out whats going wrong with the Voronoi subtraction:
//


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

//static const int nAlgos = 9;
//static const int BinLabelN = 11;
//remember to change this to run akPu3PF for pPb and akVs3PF for Pbpb datasets. or just create a separate header file which will be way easier. 
//static const char *algoName[nAlgos] = { "", "icPu5", "akPu2PF", "akPu3PF", "akPu4PF", "akPu5PF" , "akPu2Calo", "akPu3Calo", "akPu4Calo" };
//static const char *algoNamePP[nAlgos] = { "", "icPu5", "ak2PF", "ak3PF", "ak4PF", "ak5PF" , "ak2Calo", "ak3Calo", "ak4Calo" };
//static const char *algoNameGen[nAlgos] = { "", "icPu5", "akPu2PF", "akVs3PF", "akPu4PF", "akPu2PF", "akPu3PF", "akPu4PF" };
//static const char *BinLabel[BinLabelN] = {"100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180", "180-200", "200-240","240-300" };

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

    if(isPbPb){
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
      
      tEvt->SetBranchAddress("hiNevtPlane",&hiNevtPlane);
      tEvt->SetBranchAddress("hiEvtPlanes",&hiEvtPlanes);
    }

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


};


using namespace std;

void RAA_HFVsValidation_MC(char *algo = "Vs", char *jet_type = "PF"){


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
  Double_t entries[nbins_pthat]; 
  //there are two ways in which we can select the no of events we use to scale - it has to be between the pthat range. 
  //first file name - partial 50K statistics, second one is full statistics sample. 
  //similarly the entries number is for the small statistics. 
  
  // refer this twiki for the data and MC files: http://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForestPA2014#PYTHIA_HYDJET_embedded_sample

  boundaries_pthat[0]=15;
  //fileName_pthat[0] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat15_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[0] = "/export/d00/scratch/dav2105/badjets/bad15.root";
  xsection[0]= 2.034e-01;
  //entries[0] = ;//total - 48588
  
  boundaries_pthat[1]=30;
  //fileName_pthat[1] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat30_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[1] = "/export/d00/scratch/dav2105/badjets/bad30.root";
  xsection[1]= 1.075e-02;
  // entries[1] = ;//total - 48428
  
  boundaries_pthat[2]=50;
  //fileName_pthat[2] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat50_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[2] = "/export/d00/scratch/dav2105/badjets/bad50.root";
  xsection[2]= 1.025e-03;
  // entries[2] = ;//total - 50000
  
  boundaries_pthat[3]=80;
  //fileName_pthat[3] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[3] = "/export/d00/scratch/dav2105/badjets/bad80.root";
  xsection[3]= 9.865e-05;
  // entries[3] = ;//total - 49500
  
  boundaries_pthat[4]=120;
  //fileName_pthat[4] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat120_Track9_Jet30_matchEqR_merged_forest_0.root";  
  fileName_pthat[4] = "/export/d00/scratch/dav2105/badjets/bad120.root";
  xsection[4]= 1.129e-05;
  // entries[4] = ;//total - 49500

  boundaries_pthat[5]=170;
  //fileName_pthat[5] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat170_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[5] = "/export/d00/scratch/dav2105/badjets/bad120.root";
  xsection[5]= 1.465e-06;
  // entries[5] = ;//total - 49444

  boundaries_pthat[6]=220;
  //fileName_pthat[6] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[6] = "/export/d00/scratch/dav2105/badjets/bad220.root";
  xsection[6]= 2.837e-07;
  // entries[6] = ;//total - 49460

  boundaries_pthat[7]=280;
  //fileName_pthat[7] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat280_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[7] = "/export/d00/scratch/dav2105/badjets/bad280.root";
  xsection[7]= 5.323e-08;
  // entries[7] = ;//total - 49541

  boundaries_pthat[8]=370;
  //fileName_pthat[8] = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat370_Track9_Jet30_matchEqR_merged_forest_0.root";
  fileName_pthat[8] = "/export/d00/scratch/dav2105/badjets/bad370.root";
  xsection[8]= 5.934e-09;
  // entries[8] = ;//total - 19031

  boundaries_pthat[9] = 2000;
  xsection[9] = 0.0;

  // Vertex & centrality reweighting for PbPb
  TF1 *fVz;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);

  //get the centrality weight from the root file created in the plotting macro. 
  TFile *fcentin = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_DataMC_cent_ratio_20141117.root");
  TH1F *hCentWeight = (TH1F*)fcentin->Get("hCentRatio");


  TH2F *hSumpTvsHF[15];
  TH2F *hSumpTvsHF_cent[nbins_cent];
  TH2F *hSumpTvshiBin[15];
  TH2F *hSumpTvshiBin_cent[nbins_cent];
  TH2F *hNJetsvsSumpT[nbins_cent];
  TH1F *hSumpT[nbins_cent];
  TH1F *hPsi[5][15][nbins_cent];
  //TH2F *hPsi_2d[5][15][nbins_cent];

  TH2F *hPsivsHF[5][15];
  TH2F *hPsivsHF_cent[5][nbins_cent];
  TH2F *hvnvsHF[5][15];
  TH2F *hvnvscent[5][nbins_cent];


  // UseFull histograms: Event Plane from HF for the official versus HF/Vs algorithm calculation of the Event Plane.
  // first [3] array elements - Psi_2, Psi_3, Psi_4  only in the HF for tonight. 
  // [3] - total, p - positive and n - negative in eta; based on what we have from the https://github.com/CmsHI/cmssw/blob/forest_CMSSW_5_3_20/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h#L31 
  TH1F *hEP_HF_Official[3][3][nbins_cent+1];
  TH1F *hEP_HF_Vs[3][3][nbins_cent+1];
  TH2F *hEP_HF[3][3][nbins_cent+1];
  TH1F *hAngle_2_HF_Vs[nbins_cent+1];
  TH1F *hAngle_3_HF_Vs[nbins_cent+1];
  TH1F *hAngle_4_HF_Vs[nbins_cent+1];
  TH1F *hAngle_2_HF_Official[nbins_cent+1];
  TH1F *hAngle_3_HF_Official[nbins_cent+1];
  TH1F *hAngle_4_HF_Official[nbins_cent+1];

  for(int i = 0;i<nbins_cent+1;i++){
    
    for(int z = 0;z<3;z++){
      for(int x = 0;x<3;x++){
	hEP_HF_Official[z][x][i] = new TH1F(Form("hEP_HF_Official_Psi%d_%d_cent%d",z,x,i),"",630,-3.15,3.15);
	hEP_HF_Vs[z][x][i] = new TH1F(Form("hEP_HF_Vs_Psi%d_%d_cent%d",z,x,i),"",630,-3.15,3.15);
	hEP_HF[z][x][i] = new TH2F(Form("hEP_HF_Psi%d_%d_cent%d",z,x,i),"",630,-3.15,3.15,630,-3.15,3.15);
      }
    }

    hSumpTvsHF_cent[i] = new TH2F(Form("hSumpT_vsHF_cent%d",i),"",5000,0,100000,5000,0,100000);
    hSumpTvshiBin_cent[i] = new TH2F(Form("hSumpT_vshiBin_cent%d",i),"",5000,0,100000,200,0,200);
    hNJetsvsSumpT[i] = new TH2F(Form("hNJetsvsSumpT_cent%d",i),"",50,0,50,5000,0,100000);
    hAngle_2_HF_Vs[i] = new TH1F(Form("hAngle_2_HF_Vs_cent%d",i),"",630,-6.30,+6.30);
    hAngle_3_HF_Vs[i] = new TH1F(Form("hAngle_3_HF_Vs_cent%d",i),"",630,-6.30,+6.30);
    hAngle_4_HF_Vs[i] = new TH1F(Form("hAngle_4_HF_Vs_cent%d",i),"",630,-6.30,+6.30);
    hAngle_2_HF_Official[i] = new TH1F(Form("hAngle_2_HF_Official_cent%d",i),"",630,-6.30,+6.30);
    hAngle_3_HF_Official[i] = new TH1F(Form("hAngle_3_HF_Official_cent%d",i),"",630,-6.30,+6.30);
    hAngle_4_HF_Official[i] = new TH1F(Form("hAngle_4_HF_Official_cent%d",i),"",630,-6.30,+6.30);
    hSumpT[i] = new TH1F(Form("hSumpT_cent%d",i),"",5000,0,100000);
   
    /*
    for(int b = 0;b<5;b++){

      hPsivsHF_cent[b][i] = new TH2F(Form("hPsi%d_cent%d",b,i),"",630,-3.15,+3.15,5000,0,100000);
      hvnvsHF[b][i] = new TH2F(Form("hvn%d_HF_cent%d",b,i),"",100,0,1,5000,0,100000);
      //hvnvscent[b][i] = new TH2F
      for(int a = 0;a<15;a++){

	//hPsi[b][a][i] = new TH1F(Form("hPsi_n%d_eta%d_cent%d",b,a,i),"",630,-3.15,3.15);
	
      }
      
    }
    */
  }
  
  /*  
  // declare the histograms in loops:
  for(int a = 0;a<15;a++){
    
    hSumpTvsHF[a] = new TH2F(Form("hSumpT_eta%d_vsHF",a),"",5000,0,100000,5000,0,100000);
    hSumpTvshiBin[a] = new TH2F(Form("hSumpT_eta%d_vshiBin",a),"",5000,0,100000,200,0,200);
    
    for(int b = 0;b<5;b++){
      
      hPsivsHF[b][a] = new TH2F(Form("hPsi%d_etabin%d_HF",b,a),"",630,-3.15,+3.15,5000,0,100000);
      hvnvsHF[b][a] = new TH2F(Form("hvn%d_etabin%d_HF",b,a),"",100,0,1,5000,0,100000);
      
    }// no of flow components
    
  }// eta bin
  */
  TH1F *hPtHat[no_radius];
  TH1F *hPtHatRaw[no_radius];
  // Setup jet data branches - this will be 2D with [radius][pthat-file], but the histogram here is just 1D with [radius]
  JetData *data[no_radius][nbins_pthat]; 
  for(int k = 0;k<no_radius;k++){
    hPtHatRaw[k] = new TH1F(Form("hPtHatRaw_R%d",list_radius[k]),"",nbins_pthat,boundaries_pthat);
    hPtHat[k] = new TH1F(Form("hPtHat_R%d",list_radius[k]),"",nbins_pthat,boundaries_pthat);
    if(printDebug)cout<<"Radius = "<<list_radius[k]<<endl;
    if(printDebug)cout<<"reading all the pbpb mc files"<<endl;
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
  }
  
  
  for(int k = 0;k<no_radius;k++){
    if(printDebug)cout<<"Filling MC for radius = "<<list_radius[k]<<endl;
    // fill PbPb MC 
    if(printDebug)cout<<"Filling PbPb MC"<<endl;
    
    for (int h=0;h<nbins_pthat;h++) {
      if (xsection[h]==0) continue;
      if(printDebug)cout <<"Loading pthat"<<boundaries_pthat[h]<<" sample, cross section = "<<xsection[h]<< Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[h],boundaries_pthat[h+1])<<endl;
      
      //TCut pthatcut = Form("pthat>%d && pthat<%d",boundaries_pthat[h],boundaries_pthat[h+1]);
      //double fentries_test = data[k][h]->tJet->GetEntries(pthatcut);
      //cout<<"fentries_test = "<<fentries_test<<endl;
      
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
      
      for (Long64_t jentry=0; jentry<data[k][h]->tJet->GetEntries();jentry++) {
	//for (Long64_t jentry=0; jentry<10;jentry++) {
	
        //cout<<"hi"<<endl;
        data[k][h]->tEvt->GetEntry(jentry);
        data[k][h]->tJet->GetEntry(jentry);
        data[k][h]->tSkim->GetEntry(jentry);
        	
        int pthatBin = hPtHat[k]->FindBin(data[k][h]->pthat);
	
	double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/fentries;
	
	if(!data[k][h]->pcollisionEventSelection) continue;
        int cBin = findBin(data[k][h]->bin);
        //int cBin = nbins_cent-1;
        double weight_cent=1;
        double weight_pt=1;
        double weight_vz=1;
	Float_t cent = cBin;
	
        weight_cent = hCentWeight->GetBinContent(hCentWeight->FindBin(data[k][h]->bin));
	if(fabs(data[k][h]->vz)>15) continue;
	
        weight_vz = fVz->Eval(data[k][h]->vz);
	
        if(scale*weight_cent*weight_vz <=0 ) {
	  cout<<"RED FLAG RED FLAG RED FLAG"<<endl;
	  continue;
	}
	
	int jetCounter = 0;//counts jets which are going to be used in the supernova cut rejection. 
	
	for(int j = 0;j<nbins_eta;j++){
	  
	  for(int g = 0;g<data[k][h]->njets;g++){
	    
	    if(data[k][h]->jteta[g] >= boundaries_eta[j][0] && data[k][h]->jteta[g] < boundaries_eta[j][1]){
	      //cout<<"jtpt = "<<data[k][h]->jtpt[g]<<endl;
	      if(data[k][h]->jtpt[g]>=50) jetCounter++;
	      
	    }// eta selection loop
	    
	  }//jet loop
	  
	}//eta bins loop

	//constructing the x-,y- and x+,y+ based on the following equation: 
	// x- = sumpt[0]*vn[2][0]*cos(2*vpsi[2][0]) y- = sumpt[0]*vn[2][0]*sin(2*vpsi[2][0])
	Float_t Vs_2_x_minus = data[k][h]->sumpT[0]*data[k][h]->v_n[2][0]*TMath::Cos(2*data[k][h]->psi_n[2][0]);
	Float_t Vs_2_x_plus = data[k][h]->sumpT[14]*data[k][h]->v_n[2][14]*TMath::Cos(2*data[k][h]->psi_n[2][14]);
	Float_t Vs_2_y_minus = data[k][h]->sumpT[0]*data[k][h]->v_n[2][0]*TMath::Cos(2*data[k][h]->psi_n[2][0]);
	Float_t Vs_2_y_plus = data[k][h]->sumpT[14]*data[k][h]->v_n[2][14]*TMath::Cos(2*data[k][h]->psi_n[2][14]);
	Float_t Vs_2_x = Vs_2_x_minus + Vs_2_x_plus;
	Float_t Vs_2_y = Vs_2_y_minus + Vs_2_y_plus;
	Float_t Vs_2_angle = TMath::ATan2(Vs_2_y,Vs_2_x);
	hAngle_2_HF_Vs[cBin]->Fill(Vs_2_angle);
	hAngle_2_HF_Vs[nbins_cent]->Fill(Vs_2_angle);
	cout<<"Vs_2_angle = "<<Vs_2_angle<<endl;

	Float_t Vs_3_x_minus = data[k][h]->sumpT[0]*data[k][h]->v_n[3][0]*TMath::Cos(2*data[k][h]->psi_n[3][0]);
	Float_t Vs_3_x_plus = data[k][h]->sumpT[14]*data[k][h]->v_n[3][14]*TMath::Cos(2*data[k][h]->psi_n[3][14]);
	Float_t Vs_3_y_minus = data[k][h]->sumpT[0]*data[k][h]->v_n[3][0]*TMath::Cos(2*data[k][h]->psi_n[3][0]);
	Float_t Vs_3_y_plus = data[k][h]->sumpT[14]*data[k][h]->v_n[3][14]*TMath::Cos(2*data[k][h]->psi_n[3][14]);
	Float_t Vs_3_x = Vs_3_x_minus + Vs_3_x_plus;
	Float_t Vs_3_y = Vs_3_y_minus + Vs_3_y_plus;
	Float_t Vs_3_angle = TMath::ATan2(Vs_3_y,Vs_3_x);
	hAngle_3_HF_Vs[cBin]->Fill(Vs_3_angle);
	hAngle_3_HF_Vs[nbins_cent]->Fill(Vs_3_angle);
	cout<<"Vs_3_angle = "<<Vs_3_angle<<endl;

	Float_t Vs_4_x_minus = data[k][h]->sumpT[0]*data[k][h]->v_n[4][0]*TMath::Cos(2*data[k][h]->psi_n[4][0]);
	Float_t Vs_4_x_plus = data[k][h]->sumpT[14]*data[k][h]->v_n[4][14]*TMath::Cos(2*data[k][h]->psi_n[4][14]);
	Float_t Vs_4_y_minus = data[k][h]->sumpT[0]*data[k][h]->v_n[4][0]*TMath::Cos(2*data[k][h]->psi_n[4][0]);
	Float_t Vs_4_y_plus = data[k][h]->sumpT[14]*data[k][h]->v_n[4][14]*TMath::Cos(2*data[k][h]->psi_n[4][14]);
	Float_t Vs_4_x = Vs_4_x_minus + Vs_4_x_plus;
	Float_t Vs_4_y = Vs_4_y_minus + Vs_4_y_plus;
	Float_t Vs_4_angle = TMath::ATan2(Vs_4_y,Vs_4_x);
	hAngle_4_HF_Vs[cBin]->Fill(Vs_4_angle);
	hAngle_4_HF_Vs[nbins_cent]->Fill(Vs_4_angle);
      	cout<<"Vs_3_angle = "<<Vs_3_angle<<endl;

	//getting the event plane information for the Official and Vs calculations. 
	/*
	hEP_HF_Official[0][0][cBin]->Fill(data[k][h]->hiEvtPlanes[21],scale*weight_vz*weight_cent);
	hEP_HF_Official[0][1][cBin]->Fill(data[k][h]->hiEvtPlanes[22],scale*weight_vz*weight_cent);
	hEP_HF_Official[0][2][cBin]->Fill(data[k][h]->hiEvtPlanes[23],scale*weight_vz*weight_cent);

	hEP_HF_Official[1][0][cBin]->Fill(data[k][h]->hiEvtPlanes[24],scale*weight_vz*weight_cent);
	hEP_HF_Official[1][1][cBin]->Fill(data[k][h]->hiEvtPlanes[25],scale*weight_vz*weight_cent);
	hEP_HF_Official[1][2][cBin]->Fill(data[k][h]->hiEvtPlanes[26],scale*weight_vz*weight_cent);

	hEP_HF_Official[2][0][cBin]->Fill(data[k][h]->hiEvtPlanes[27],scale*weight_vz*weight_cent);
	hEP_HF_Official[2][1][cBin]->Fill(data[k][h]->hiEvtPlanes[28],scale*weight_vz*weight_cent);
	hEP_HF_Official[2][2][cBin]->Fill(data[k][h]->hiEvtPlanes[29],scale*weight_vz*weight_cent);

	hEP_HF_Vs[0][0][cBin]->Fill(data[k][h]->psi_n[2][0],scale*weight_vz*weight_cent);
	hEP_HF_Vs[0][0][cBin]->Fill(data[k][h]->psi_n[2][14],scale*weight_vz*weight_cent);
	hEP_HF_Vs[0][1][cBin]->Fill(data[k][h]->psi_n[2][0],scale*weight_vz*weight_cent);
	hEP_HF_Vs[0][2][cBin]->Fill(data[k][h]->psi_n[2][14],scale*weight_vz*weight_cent);

	hEP_HF_Vs[1][0][cBin]->Fill(data[k][h]->psi_n[3][0],scale*weight_vz*weight_cent);
	hEP_HF_Vs[1][0][cBin]->Fill(data[k][h]->psi_n[3][14],scale*weight_vz*weight_cent);
	hEP_HF_Vs[1][1][cBin]->Fill(data[k][h]->psi_n[3][0],scale*weight_vz*weight_cent);
	hEP_HF_Vs[1][2][cBin]->Fill(data[k][h]->psi_n[3][14],scale*weight_vz*weight_cent);

	hEP_HF_Vs[2][0][cBin]->Fill(data[k][h]->psi_n[4][0],scale*weight_vz*weight_cent);
	hEP_HF_Vs[2][0][cBin]->Fill(data[k][h]->psi_n[4][14],scale*weight_vz*weight_cent);
	hEP_HF_Vs[2][1][cBin]->Fill(data[k][h]->psi_n[4][0],scale*weight_vz*weight_cent);
	hEP_HF_Vs[2][2][cBin]->Fill(data[k][h]->psi_n[4][14],scale*weight_vz*weight_cent);
	
	hEP_HF[0][0][cBin]->Fill(data[k][h]->hiEvtPlanes[21],data[k][h]->psi_n[2][0],scale*weight_vz*weight_cent);
	hEP_HF[0][0][cBin]->Fill(data[k][h]->hiEvtPlanes[21],data[k][h]->psi_n[2][14],scale*weight_vz*weight_cent);
	hEP_HF[0][1][cBin]->Fill(data[k][h]->hiEvtPlanes[22],data[k][h]->psi_n[2][0],scale*weight_vz*weight_cent);
	hEP_HF[0][2][cBin]->Fill(data[k][h]->hiEvtPlanes[23],data[k][h]->psi_n[2][14],scale*weight_vz*weight_cent);

	hEP_HF[1][0][cBin]->Fill(data[k][h]->hiEvtPlanes[24],data[k][h]->psi_n[3][0],scale*weight_vz*weight_cent);
	hEP_HF[1][0][cBin]->Fill(data[k][h]->hiEvtPlanes[24],data[k][h]->psi_n[3][14],scale*weight_vz*weight_cent);
	hEP_HF[1][1][cBin]->Fill(data[k][h]->hiEvtPlanes[25],data[k][h]->psi_n[3][0],scale*weight_vz*weight_cent);
	hEP_HF[1][2][cBin]->Fill(data[k][h]->hiEvtPlanes[26],data[k][h]->psi_n[3][14],scale*weight_vz*weight_cent);

	hEP_HF[2][0][cBin]->Fill(data[k][h]->hiEvtPlanes[27],data[k][h]->psi_n[4][0],scale*weight_vz*weight_cent);
	hEP_HF[2][0][cBin]->Fill(data[k][h]->hiEvtPlanes[27],data[k][h]->psi_n[4][14],scale*weight_vz*weight_cent);
	hEP_HF[2][1][cBin]->Fill(data[k][h]->hiEvtPlanes[28],data[k][h]->psi_n[4][0],scale*weight_vz*weight_cent);
	hEP_HF[2][2][cBin]->Fill(data[k][h]->hiEvtPlanes[29],data[k][h]->psi_n[4][14],scale*weight_vz*weight_cent);
	*/

	float sumpTtotal = 0;
	for(int a = 0;a<15;a++){

	  //hSumpTvsHF[a]->Fill(data[k][h]->sumpT[a],data[k][h]->hiHF,scale*weight_vz*weight_cent);
	  //hSumpTvshiBin[a]->Fill(data[k][h]->sumpT[a],data[k][h]->bin,scale*weight_vz*weight_cent);
	  sumpTtotal+= data[k][h]->sumpT[a];
	  //for(int b = 0;b<5;b++){
	    
	  //  hPsivsHF[b][a]->Fill(data[k][h]->psi_n[b][a],data[k][h]->hiHF,scale*weight_vz*weight_cent);
	  //  hvnvsHF[b][a]->Fill(data[k][h]->v_n[b][a],data[k][h]->hiHF,scale*weight_vz*weight_cent);
	    //hPsi[cBin]->Fill(data[k][h]->psi_n[b][a]);

	  //}// no of flow components 
	}// eta bin
	
	hNJetsvsSumpT[cBin]->Fill(jetCounter,sumpTtotal,scale*weight_vz*weight_cent);
	hNJetsvsSumpT[nbins_cent]->Fill(jetCounter,sumpTtotal,scale*weight_vz*weight_cent);
	hSumpTvsHF_cent[cBin]->Fill(sumpTtotal,data[k][h]->hiHF,scale*weight_vz*weight_cent);
	hSumpTvsHF_cent[nbins_cent]->Fill(sumpTtotal,data[k][h]->hiHF,scale*weight_vz*weight_cent);
	//hSumpTvshiBin_cent[cBin]->Fill(data[k][h]->sumpT[a],data[k][h]->hiHF);
	hSumpT[cBin]->Fill(sumpTtotal,scale*weight_vz*weight_cent);
	hSumpT[nbins_cent]->Fill(sumpTtotal,scale*weight_vz*weight_cent);
       
	
      }// event entry loop

    }// pthat loop

  }// radius loop

  TFile fout(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_MC_HFVsValidation_failure_histograms_ak%s%s_%d.root",algo,jet_type,date.GetDate()),"RECREATE");
  fout.cd();

  for(int i = 0;i<=nbins_cent;i++){
    hNJetsvsSumpT[i]->Write();
    hSumpTvsHF_cent[i]->Write();
    hSumpT[i]->Write();
    hAngle_2_HF_Vs[i]->Write();
    hAngle_3_HF_Vs[i]->Write();
    hAngle_4_HF_Vs[i]->Write();
  }

}
