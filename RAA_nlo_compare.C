// Raghav Kunnawalkam Elayavalli
// July 1st 2014
// CERN

//
// Macro to readin the NLO files and compare it with the pp data we have. 
// 

// I need to decide if this macro will also plot the reuqired histograms or would i have to create a new macro for that
// going by what ive done before i think ill create a new macro. 

// July 11 - we also have the NP corrections here /net/hisrv0001/home/rkunnawa/WORK/NP_factors/test/txtfiles

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

  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};


using namespace std;

static const int dir = 50;

void RAA_nlo_compare(int radius = 3, int energy = 2760){

  TStopwatch timer;
  timer.Start();
  
  TDatime date;
  
  gStyle->SetOptStat(0);
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  cout<<"Starting the comparison of the NLO with data"<<endl;
  cout<<"Running for Energy = "<<energy<<" and Radius = "<<radius<<endl;
  
  //nlo files: right now they are not for all the eta ranges. they are only present 
  TFile *fNLO_nnpdf = TFile::Open("nlo_files/input_rootfiles/fnl4350a_nnpdf21-nlo_aspdf_new.root");
  TFile *fNLO_cteq = TFile::Open("nlo_files/input_rootfiles/fnl4350a_cteq66-nlo_aspdf_all_new.root");
  TFile *fNLO_ct10n = TFile::Open("nlo_files/input_rootfiles/fnl4350a_ct10n-nlo_aspdf_new.root");
  TFile *fNLO_hera = TFile::Open("nlo_files/input_rootfiles/fnl4350a_hera15all-nlo_aspdf_new.root");

  //For the moment, the NLO histograms will have the radius as the array counter. (2,3,4). later i will include the eta bins as the other array counter. which will make these 
  TH1F *hPP_nnpdf[3],*hPP_cteq[3],*hPP_ct10n[3],*hPP_hera[3];

  TH1F *hPP_err[3];

  for(int i = 0;i<3;i++){

    hPP_nnpdf[i] = (TH1F*)fNLO_nnpdf->Get(Form("h100%d00",i+1));
    hPP_cteq[i] = (TH1F*)fNLO_cteq->Get(Form("h100%d00",i+1));
    hPP_ct10n[i] = (TH1F*)fNLO_ct10n->Get(Form("h100%d00",i+1));
    hPP_hera[i] = (TH1F*)fNLO_hera->Get(Form("h100%d00",i+1));
    
    hPP_err[i] = (TH1F*)fNLO_cteq->Get(Form("h100%d03",i+1));

    for(int j = 0;j<hPP_nnpdf[0]->GetNbinsX();j++){
      //Float_t valErr = hPP_err[i]->GetBinError(j);
      Float_t valErr = hPP_err[i]->GetBinContent(j);
      hPP_nnpdf[i]->SetBinError(j,valErr);
      hPP_cteq[i]->SetBinError(j,valErr);
      hPP_ct10n[i]->SetBinError(j,valErr);
      hPP_hera[i]->SetBinError(j,valErr);
      
    }

  }
  
  //npc files: 
  // for now only take the central eta file for the different radii (# 1 + i*10 in the array). then we can take all the rest. 
  /*
  char dirName[dir][256] = {
    "ak2GenJetSpectrum_n10_p10","ak2GenJetSpectrum_n20_p20","ak2GenJetSpectrum_n25_n20","ak2GenJetSpectrum_n20_n15",
    "ak2GenJetSpectrum_n15_n10","ak2GenJetSpectrum_n10_n05","ak2GenJetSpectrum_n05_p05","ak2GenJetSpectrum_p05_p10",
    "ak2GenJetSpectrum_p10_p15","ak2GenJetSpectrum_p15_p20",
    "ak3GenJetSpectrum_n10_p10","ak3GenJetSpectrum_n20_p20","ak3GenJetSpectrum_n25_n20","ak3GenJetSpectrum_n20_n15",
    "ak3GenJetSpectrum_n15_n10","ak3GenJetSpectrum_n10_n05","ak3GenJetSpectrum_n05_p05","ak3GenJetSpectrum_p05_p10",
    "ak3GenJetSpectrum_p10_p15","ak3GenJetSpectrum_p15_p20",
    "ak4GenJetSpectrum_n10_p10","ak4GenJetSpectrum_n20_p20","ak4GenJetSpectrum_n25_n20","ak4GenJetSpectrum_n20_n15",
    "ak4GenJetSpectrum_n15_n10","ak4GenJetSpectrum_n10_n05","ak4GenJetSpectrum_n05_p05","ak4GenJetSpectrum_p05_p10",
    "ak4GenJetSpectrum_p10_p15","ak4GenJetSpectrum_p15_p20",
    "ak5GenJetSpectrum_n10_p10","ak5GenJetSpectrum_n20_p20","ak5GenJetSpectrum_n25_n20","ak5GenJetSpectrum_n20_n15",
    "ak5GenJetSpectrum_n15_n10","ak5GenJetSpectrum_n10_n05","ak5GenJetSpectrum_n05_p05","ak5GenJetSpectrum_p05_p10",
    "ak5GenJetSpectrum_p10_p15","ak5GenJetSpectrum_p15_p20",
    "ak7GenJetSpectrum_n10_p10","ak7GenJetSpectrum_n20_p20","ak7GenJetSpectrum_n25_n20","ak7GenJetSpectrum_n20_n15",
    "ak7GenJetSpectrum_n15_n10","ak7GenJetSpectrum_n10_n05","ak7GenJetSpectrum_n05_p05","ak7GenJetSpectrum_p05_p10",
    "ak7GenJetSpectrum_p10_p15","ak7GenJetSpectrum_p15_p20"};
  */
  char etaWidth[dir][256] = {
    "n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
    "n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
    "p10_eta_p15","p15_eta_p20",
    "n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
    "n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
    "p10_eta_p15","p15_eta_p20",
    "n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
    "n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
    "p10_eta_p15","p15_eta_p20",
    "n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
    "n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
    "p10_eta_p15","p15_eta_p20",
    "n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
    "n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
    "p10_eta_p15","p15_eta_p20"
  };
  
  char radius_lable[dir][256] = {
    "R2","R2","R2","R2","R2","R2","R2","R2","R2","R2",
    "R3","R3","R3","R3","R3","R3","R3","R3","R3","R3",
    "R4","R4","R4","R4","R4","R4","R4","R4","R4","R4",
    "R5","R5","R5","R5","R5","R5","R5","R5","R5","R5",
    "R7","R7","R7","R7","R7","R7","R7","R7","R7","R7"
  };
  
  ifstream fin_txt[dir];
  
  for(int i = 0;i<dir;i++){
    
    //ostringstream filename;
    //filename<<"/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Macros/RAA/nlo_files/input_np_txtfiles/NPC_ak_"<<radius_lable[i]<<etaWidth[i]<<"_energy"<<energy<<".txt";
    fin_txt[i].open(Form("nlo_files/input_np_txtfiles/NPC_ak_%s_%s_energy%d.txt",radius_lable[i],etaWidth[i],energy));
    
  }

  //apply the correction factors to the NLO histograms 

  Float_t npc = 0;
  Float_t bin = 0;

  //int counter = 1;
  /*
  while(1){

    for(int i = 0;i<hPP_nnpdf[0]->GetNbinsX();i++){
      fin_txt[i]>>bin>>npc;
      if(!fin_txt[i].good())break;
      for(int j = 0;j<3;j++){
	hPP_nnpdf[j]->SetBinContent(i,npc*hPP_nnpdf[j]->GetBinContent(i));
	hPP_ct10n[j]->SetBinContent(i,npc*hPP_ct10n[j]->GetBinContent(i));
	hPP_cteq[j]->SetBinContent(i,npc*hPP_cteq[j]->GetBinContent(i));
	hPP_hera[j]->SetBinContent(i,npc*hPP_hera[j]->GetBinContent(i));
      }
    }
  }
  */

  //have to think about the best way to get the np factors from the files to the histograms without falling into a loop mess. 
  // so we have histograms compared with txt files for those histograms. total of 50 txt files and (ideally should have 50 histograms as well) with only 3 histograms, central eta range for R=0.2,0.3,0.4 for the different nlo. 
  hPP_nnpdf[0]->Print("base");

  int counter = 0;

  for(int i = 0;i<3;i++){

    while(1){

      fin_txt[i*10+1]>>bin>>npc;
      if(!fin_txt[i*10+1].good()) break;

      hPP_nnpdf[i]->SetBinContent(counter,npc*hPP_nnpdf[i]->GetBinContent(counter));
      hPP_cteq[i]->SetBinContent(counter,npc*hPP_cteq[i]->GetBinContent(counter));
      hPP_ct10n[i]->SetBinContent(counter,npc*hPP_ct10n[i]->GetBinContent(counter));
      hPP_hera[i]->SetBinContent(counter,npc*hPP_hera[i]->GetBinContent(counter));

      counter++;
    }


  }

  hPP_nnpdf[0]->Print("base");

  //the nlo histograms are now ready

  // lets load the data here: 
  

  timer.Stop();
  float rtime  = timer.RealTime();
  float ctime  = timer.CpuTime();
  
  cout<<"\t"<<endl;
  cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<endl;
  cout<<"\t"<<endl;
  cout<<"Good bye : " <<"\t"<<endl;
  



}









































