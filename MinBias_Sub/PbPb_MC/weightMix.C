#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TChain.h"

#include "TCut.h"
#include "TNtuple.h"

#include "THStack.h"

using namespace std;

static bool subtract = 0;
static int QID = 3;

static TString weightString;

static bool normLead = 0;

static int mixColor = 2;
static int dataColor = 1;
static int ppColor = 4;

static int centralBin = 8;
static int leadCut = 120;
static int subleadCut = 30;

static double sideMin = 0.1;
static double sideMax = TMath::Pi()/3 + 0.1;

static double sideCorrect = 1;

static const char* LUM = "#int L dt=150";

static bool plotSubtraction = 0;

static bool reweightCentrality = 1;

static const double pi = TMath::Pi();


void weightMix(){


  TH1::SetDefaultSumw2();

  int Npt = 10;

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForestPA2013#Dijet_Cross_Sections_for_reweigh
  double pthatBins2011[] = {15,30,50,80,120,170,220,280,370,9999};
  double xs2011[] = {2.034e-01, 1.075E-02, 1.025E-03, 9.865E-05, 1.129E-05, 1.465E-06, 2.837E-07, 5.323e-08, 5.934e-09, 0};

  double *xs;
  double *pthats;
  int nev[20];
  int n[20];

  Npt = 9;
  xs = xs2011;
  pthats = pthatBins2011;

  for(int i = 0; i < Npt-1; ++i){
    xs[i] -= xs[i+1];
  }

  TChain* nt;
  TChain* evt;
  TFile* outf;

  nt = new TChain("akPu3PFJetAnalyzer/t");
  evt = new TChain("hiEvtAnalyzer/HiTree");
  // Order of files matters a lot!
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat15_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat30_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat50_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat120_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat170_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat280_Track9_Jet30_matchEqR_merged_forest_0.root");
  nt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat370_Track9_Jet30_matchEqR_merged_forest_0.root");

  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat15_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat30_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat50_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat120_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat170_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat280_Track9_Jet30_matchEqR_merged_forest_0.root");
  evt->AddFile("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat370_Track9_Jet30_matchEqR_merged_forest_0.root");

  nt->AddFriend(evt);

  outf = new TFile(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/weights.root"),"recreate");
  TNtuple* ntw = new TNtuple("weights","","weight:ptweight:cweight:varpthat:vzweight");
  float varpthat, varbin, varvz;

  TF1 *fVz=0;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);

  nt->SetBranchAddress("pthat",&varpthat);
  evt->SetBranchAddress("hiBin",&varbin);
  evt->SetBranchAddress("vz",&varvz);

  for(int i = 0; i < Npt; ++i){
    TCut pthatCut(Form("pthat >= %f && pthat < %f",pthats[i],pthats[i+1]));
    n[i] = nt->GetEntries(pthatCut);
  }

  for(int ie = 0; ie < nt->GetEntries(); ++ie){

    if(ie%50000 == 0) cout<<ie<<"/"<<nt->GetEntries();
    
    nt->GetEntry(ie);
    double pthatweight = 0;
    
    for(int i = 0; i < Npt; ++i){
      if(n[i] > 0 && varpthat >= pthats[i]) pthatweight = xs[i]/n[i];
    }

    double cweight = exp(-  pow(varbin+1.11957e+01,2)  /    pow(1.34120e+01,2) / 2);
    double vzweight = fVz->Eval(varvz);    
    ntw->Fill(pthatweight*cweight*vzweight,pthatweight,cweight,varpthat,vzweight);
  }

  outf->cd();
  ntw->Write();
  outf->Write();

}

