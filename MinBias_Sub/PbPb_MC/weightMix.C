#include <iostream>
#include <fstream>
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
static const int Nfiles = 9; 

void weightMix(){

  TH1::SetDefaultSumw2();

  int Npt = 10;

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForestPA2013#Dijet_Cross_Sections_for_reweigh
  double pthatBins2011[] = {15,30,50,80,120,170,220,280,370,9999};
  double xs2011[] = {2.034e-01, 1.075E-02, 1.025E-03, 9.865E-05, 1.129E-05, 1.465E-06, 2.837E-07, 5.323e-08, 5.934e-09, 0};

  double *xs;
  double *pthats;
  int nev[20];
  int n[9] = {315269, 243155, 379848, 362476, 359581, 359272, 217873, 88991, 25252};
// no of events in pthat = 15 = 315269
// no of events in pthat = 30 = 243155
// no of events in pthat = 50 = 379848
// no of events in pthat = 80 = 362476
// no of events in pthat = 120 = 359581
// no of events in pthat = 170 = 359272
// no of events in pthat = 220 = 217873
// no of events in pthat = 280 = 88991
// no of events in pthat = 370 = 25252
  Npt = 9;
  xs = xs2011;
  pthats = pthatBins2011;

  for(int i = 0; i < Npt-1; ++i){
    xs[i] -= xs[i+1];
  }

  //TChain* nt;
  //TChain* evt;
  TFile* outf[Nfiles];

  //nt->AddFriend(evt);

  std::string infile_Forest;
  infile_Forest = "jetRAA_PbPb_mc_forest.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  std::string filename_array[Nfiles];

  for(int ifile = 0; ifile<Nfiles; ++ifile){
    instr_Forest>>filename_Forest;
    //nt->AddFile(filename_Forest.c_str());
    filename_array[ifile] = filename_Forest.c_str();
  }

  for(int i = 0; i<Nfiles; ++i){

    cout<<"Now looping over individual file : "<<filename_array[i].c_str()<<endl;

    TFile * fin = new TFile(filename_array[i].c_str());
    TTree * newTree = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
    TTree * evtTree = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");

    double pthatweight = 0;
    float varpthat;
    int hiBin, evnt, lumi;
    float vz;
    
    newTree->SetBranchAddress("pthat",&varpthat);
    evtTree->SetBranchAddress("hiBin",&hiBin);
    evtTree->SetBranchAddress("evt", &evnt);
    evtTree->SetBranchAddress("lumi", &lumi);
    evtTree->SetBranchAddress("vz", &vz);

    outf[i] = new TFile(Form("weights_pbpb_%d.root",i+1),"recreate");
    
    TTree *ntw = new TTree("weights","");

    ntw->Branch("pthatweight",&pthatweight,"pthatweight/D");
    ntw->Branch("hiBin",&hiBin,"hiBin/I");
    ntw->Branch("evt",&evnt,"evt/I");
    ntw->Branch("lumi",&lumi,"lumi/I");
    ntw->Branch("vz",&vz,"vz/F");
  
    cout<<"loaded all the files"<<endl;
   
    for(int ie = 0; ie < newTree->GetEntries(); ++ie){

      if(ie%50000 == 0) cout<<ie<<"/"<<newTree->GetEntries()<<endl;
    
      newTree->GetEntry(ie);    
    
      for(int i = 0; i < Npt; ++i){
	if(n[i] > 0 && varpthat >= pthats[i]) pthatweight = xs[i]/n[i];
      }
      ntw->Fill();

    }

    outf[i]->cd();
    ntw->Write();
    outf[i]->Write();
    outf[i]->Close();
    
  }

}

