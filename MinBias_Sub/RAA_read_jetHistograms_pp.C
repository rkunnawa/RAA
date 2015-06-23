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
#include "../Headers/plot.h"


// static const int nbins_pt = 32;
// static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 272, 300, 330, 362, 395, 501};

// static const int nbins_pt = 30;
// static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};

static const int nbins_pt = 17;
static const double boundaries_pt[nbins_pt+1] = { 32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330};

// the following bins is the cms pp nlo pt bins
// static const int nbins_pt = 29;
// static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638,790,967};

// the following bins is the atlas spectra pt binning
// static const int nbins_pt = 12;
// static const double boundaries_pt[nbins_pt+1] = {31., 39., 50., 63., 79., 100., 125., 158., 199., 251., 316., 398., 501};

// the following//  bins is the atlas Rcp pt binning
// static const int nbins_pt = 12;
// static const double boundaries_pt[nbins_pt+1] = {38.36, 44.21, 50.94, 58.7, 67.64 , 77.94 , 89.81, 103.5, 119.3, 137.4 , 158.3, 182.5,  210.3};

const double kdelrcut=0.3;

using namespace std;

void RAA_read_jetHistograms_pp(char * etaWidth = (char*)"20_eta_20", 
			       Int_t radius = 3, 
			       Int_t etaLow = 20, 
			       Int_t etaHigh = 20,
			       char * ptbins = (char*)"finebinscut")
{

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // char * etaWidth = (char*)Form("%d_eta_%d",etaLow, etaHigh);
  cout<<"etaWidth = "<<etaWidth<<endl;
  cout<<"Radius = "<<radius<<endl;

  bool isSymm = false;
  if(etaLow == etaHigh) isSymm = true;
  char * ntuple = (char*)"Pawan"; //  or "Pawan"
  
  // the cut is a 3 step cut based on the different value of the calopt/pfpt - copy the following lines into your loop (with the corresponding branch address set)
  // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) hGood->Fill();
  // if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) ) hGood->Fill();
  // if(calopt/pfpt > 0.85 & eMax/Sumcand > 0.9) hGood->Fill();

  TFile * fData, * fMC;
  TTree * Data_matched, * Data_unmatched, * MC_matched, * MC_unmatched; 

  if(ntuple == "Pawan"){
    fData = TFile::Open("/mnt/hadoop/cms/store/user/pawan/trees/JetRaaTree_ak234_pp_DATA.root");
    fMC = TFile::Open("/mnt/hadoop/cms/store/user/pawan/trees/JetRaaTree_ak234_pp_MC.root");

    Data_matched = (TTree*)fData->Get(Form("ak%dJetAnalyzer/jetTree",radius));
    MC_matched = (TTree*)fMC->Get(Form("ak%dJetAnalyzer/jetTree",radius));
  }
  if(ntuple == "Raghav"){
    // Raghav's ntuples - running to find the difference in spectra between Pawan's and mine.
    fData = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/pp_Data_calo_pf_jet_correlation_deltaR_0p2_ak%d_20150331.root",radius));
    fMC = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/pp_MC_calo_pf_jet_correlation_deltaR_0p2_ak%d_20150331.root",radius));

    Data_matched= (TTree*)fData->Get("matchedJets");
    Data_unmatched = (TTree*)fData->Get("unmatchedPFJets");

    MC_matched = (TTree*)fMC->Get("matchedJets");
    MC_unmatched = (TTree*)fMC->Get("unmatchedPFJets");
  }
  
  // setup the residual correction factors 
  // TF1 = 1 - [0]/pow(x,[1]), in the pt range till 180. after 180, all the factors are 1. 

  TF1 * fResidual = new TF1("fResidual","1 - [0]/pow(x,[1])");
  if(radius == 2){
    fResidual->SetParameter(0, -0.05756);
    fResidual->SetParameter(1,  0.42750);
  }
  if(radius == 3){
    fResidual->SetParameter(0, -0.66229);
    fResidual->SetParameter(1,  1.02119);
  }
  if(radius == 4){
    fResidual->SetParameter(0, -1.28178);
    fResidual->SetParameter(1,  1.17348);
  }

  // Vertex reweighting for pp
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
  

  TH1F * hMC_Jet40_noCut = new TH1F("hMC_Jet40_noCut","data from matched jets without any jet ID cut",501,0,501);
  TH1F * hMC_Jet40_CutA = new TH1F("hMC_Jet40_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",501,0,501);
  TH1F * hMC_Jet40_CutA_rej = new TH1F("hMC_Jet40_CutA_rej","data from matched jets rejected by Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",501,0,501);

  TH1F * hMC_Jet60_noCut = new TH1F("hMC_Jet60_noCut","data from matched jets without any jet ID cut",501,0,501);
  TH1F * hMC_Jet60_CutA = new TH1F("hMC_Jet60_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",501,0,501);
  TH1F * hMC_Jet60_CutA_rej = new TH1F("hMC_Jet60_CutA_rej","data from matched jets rejected with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",501,0,501);

  TH1F * hMC_Jet80_noCut = new TH1F("hMC_Jet80_noCut","data from matched jets without any jet ID cut",501,0,501);
  TH1F * hMC_Jet80_CutA = new TH1F("hMC_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",501,0,501);
  TH1F * hMC_Jet80_CutA_rej = new TH1F("hMC_Jet80_CutA_rej","data from matched jets rejected with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",501,0,501);

  TH1F * hData_Jet40_noCut = new TH1F("hData_Jet40_noCut","data from matched jets without any jet ID cut",501,0,501);
  TH1F * hData_Jet40_CutA = new TH1F("hData_Jet40_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",501,0,501);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet40_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet40_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet40_CutB->Fill()} 
  TH1F * hData_Jet40_CutA_rej = new TH1F("hData_Jet40_CutA_rej","",501,0,501);

  TH1F * hData_Jet60_noCut = new TH1F("hData_Jet60_noCut","data from matched jets without any jet ID cut",501,0,501);
  TH1F * hData_Jet60_CutA = new TH1F("hData_Jet60_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",501,0,501);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet60_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet60_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet60_CutB->Fill()} 
  TH1F * hData_Jet60_CutA_rej = new TH1F("hData_Jet60_CutA_rej","",501,0,501);
  
  TH1F * hData_Jet80_noCut = new TH1F("hData_Jet80_noCut","data from matched jets without any jet ID cut",501,0,501);
  TH1F * hData_Jet80_CutA = new TH1F("hData_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",501,0,501);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet80_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet80_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet80_CutB->Fill()} 
  TH1F * hData_Jet80_CutA_rej = new TH1F("hData_Jet80_CutA_rej","",501,0,501);

  TH1F * hData_unmatched_Jet80_noCut = new TH1F("hData_unmatched_Jet80_noCut","",501,0,501);
  TH1F * hData_unmatched_Jet80_CutA = new TH1F("hData_unmatched_Jet80_CutA","",501,0,501);
  TH1F * hData_unmatched_Jet80_CutA_rej = new TH1F("hData_unmatched_Jet80_CutA_rej","",501,0,501);
  TH1F * hData_unmatched_Jet60_noCut = new TH1F("hData_unmatched_Jet60_noCut","",501,0,501);
  TH1F * hData_unmatched_Jet60_CutA = new TH1F("hData_unmatched_Jet60_CutA","",501,0,501);
  TH1F * hData_unmatched_Jet60_CutA_rej = new TH1F("hData_unmatched_Jet60_CutA_rej","",501,0,501);
  TH1F * hData_unmatched_Jet40_noCut = new TH1F("hData_unmatched_Jet40_noCut","",501,0,501);
  TH1F * hData_unmatched_Jet40_CutA = new TH1F("hData_unmatched_Jet40_CutA","",501,0,501);
  TH1F * hData_unmatched_Jet40_CutA_rej = new TH1F("hData_unmatched_Jet40_CutA_rej","",501,0,501);

  TH1F * hMC_unmatched_Jet80_noCut = new TH1F("hMC_unmatched_Jet80_noCut","",501,0,501);
  TH1F * hMC_unmatched_Jet80_CutA = new TH1F("hMC_unmatched_Jet80_CutA","",501,0,501);
  TH1F * hMC_unmatched_Jet80_CutA_rej = new TH1F("hMC_unmatched_Jet80_CutA_rej","",501,0,501);
  TH1F * hMC_unmatched_Jet60_noCut = new TH1F("hMC_unmatched_Jet60_noCut","",501,0,501);
  TH1F * hMC_unmatched_Jet60_CutA = new TH1F("hMC_unmatched_Jet60_CutA","",501,0,501);
  TH1F * hMC_unmatched_Jet60_CutA_rej = new TH1F("hMC_unmatched_Jet60_CutA_rej","",501,0,501);
  TH1F * hMC_unmatched_Jet40_noCut = new TH1F("hMC_unmatched_Jet40_noCut","",501,0,501);
  TH1F * hMC_unmatched_Jet40_CutA = new TH1F("hMC_unmatched_Jet40_CutA","",501,0,501);
  TH1F * hMC_unmatched_Jet40_CutA_rej = new TH1F("hMC_unmatched_Jet40_CutA_rej","",501,0,501);

  
  //get the spectra with the specific trigger object from the different files. 
  TH1F *hpp_Jet80_gen,*hpp_Jet80_reco;
  TH1F *hpp_Jet60_gen,*hpp_Jet60_reco;
  TH1F *hpp_Jet40_gen,*hpp_Jet40_reco;
  TH1F *hpp_JetComb_gen,*hpp_JetComb_reco;

  TH1F *hpp_anaBin_Jet80_gen,*hpp_anaBin_Jet80_reco;
  TH1F *hpp_anaBin_Jet60_gen,*hpp_anaBin_Jet60_reco;
  TH1F *hpp_anaBin_Jet40_gen,*hpp_anaBin_Jet40_reco;
  TH1F *hpp_anaBin_JetComb_gen,*hpp_anaBin_JetComb_reco;
  
  TH1F *hpp_gen,*hpp_reco;
  TH2F *hpp_matrix;
  TH2F *hpp_matrix_HLT;
  TH2F *hpp_anaBin_matrix_HLT;
  TH2F *hpp_mcclosure_matrix;
  TH2F *hpp_mcclosure_matrix_HLT;
  //TH2F *hpp_response;
  TH1F *hpp_mcclosure_JetComb_data;
  TH1F *hpp_mcclosure_data;
  TH1F *hpp_mcclosure_data_train;
  TH1F *hpp_mcclosure_JetComb_data_train;
  TH1F *hpp_mcclosure_Jet80_data_train;
  TH1F *hpp_mcclosure_Jet60_data_train;
  TH1F *hpp_mcclosure_Jet40_data_train;
  TH1F *hpp_mcclosure_Jet80_data;
  TH1F *hpp_mcclosure_Jet60_data;
  TH1F *hpp_mcclosure_Jet40_data;
  TH1F *hpp_mcclosure_gen;
  TH1F *hpp_mcclosure_JetComb_gen;
  TH1F *hpp_mcclosure_Jet80_gen;
  TH1F *hpp_mcclosure_Jet60_gen;
  TH1F *hpp_mcclosure_Jet40_gen;

  TH1F *hpp_TrgObj80;
  TH1F *hpp_TrgObj60;
  TH1F *hpp_TrgObj40;
  TH1F *hpp_TrgObjComb;

  TH1F *hpp_JEC_TrgObj80;
  TH1F *hpp_JEC_TrgObj60;
  TH1F *hpp_JEC_TrgObj40;
  TH1F *hpp_JEC_TrgObjComb;
  
  TH1F *hpp_Smear_TrgObj80;
  TH1F *hpp_Smear_TrgObj60;
  TH1F *hpp_Smear_TrgObj40;
  TH1F *hpp_Smear_TrgObjComb;
  
  TH1F *hpp_anaBin_TrgObj80;
  TH1F *hpp_anaBin_TrgObj60;
  TH1F *hpp_anaBin_TrgObj40;
  TH1F *hpp_anaBin_TrgObjComb;
  
  TH1F * hpp_Data_Jet80_noCut = new TH1F("hpp_Data_Jet80_noCut","",501,0,501);
  TH1F * hpp_Data_Jet60_noCut = new TH1F("hpp_Data_Jet60_noCut","",501,0,501);
  TH1F * hpp_Data_Jet40_noCut = new TH1F("hpp_Data_Jet40_noCut","",501,0,501);
  TH1F * hpp_Data_Comb_noCut = new TH1F("hpp_Data_Comb_noCut","",501,0,501);

  TH1F * hpp_MC_Jet80_noCut = new TH1F("hpp_MC_Jet80_noCut","",501,0,501);
  TH1F * hpp_MC_Jet60_noCut = new TH1F("hpp_MC_Jet60_noCut","",501,0,501);
  TH1F * hpp_MC_Jet40_noCut = new TH1F("hpp_MC_Jet40_noCut","",501,0,501);
  TH1F * hpp_MC_Comb_noCut = new TH1F("hpp_MC_Comb_noCut","",501,0,501);


  hpp_TrgObj80 = new TH1F(Form("hpp_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_TrgObj60 = new TH1F(Form("hpp_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_TrgObj40 = new TH1F(Form("hpp_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_TrgObjComb = new TH1F(Form("hpp_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),501,0,501);

  hpp_JEC_TrgObj80 = new TH1F(Form("hpp_JEC_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_JEC_TrgObj60 = new TH1F(Form("hpp_JEC_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_JEC_TrgObj40 = new TH1F(Form("hpp_JEC_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_JEC_TrgObjComb = new TH1F(Form("hpp_JEC_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),501,0,501);

  hpp_Smear_TrgObj80 = new TH1F(Form("hpp_Smear_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_Smear_TrgObj60 = new TH1F(Form("hpp_Smear_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_Smear_TrgObj40 = new TH1F(Form("hpp_Smear_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_Smear_TrgObjComb = new TH1F(Form("hpp_Smear_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),501,0,501);

  hpp_anaBin_TrgObj80 = new TH1F(Form("hpp_anaBin_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_TrgObj60 = new TH1F(Form("hpp_anaBin_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_TrgObj40 = new TH1F(Form("hpp_anaBin_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_TrgObjComb = new TH1F(Form("hpp_anaBin_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  
  hpp_gen = new TH1F(Form("hpp_gen_R%d_%s",radius,etaWidth),Form("Gen refpt R%d %s ",radius,etaWidth),501,0,501);
  //cout<<"A"<<endl;
  hpp_reco = new TH1F(Form("hpp_reco_R%d_%s",radius,etaWidth),Form("Reco jtpt R%d %s ",radius,etaWidth),501,0,501);
  //cout<<"B"<<endl;
  hpp_matrix = new TH2F(Form("hpp_matrix_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt R%d %s ",radius,etaWidth),501,0,501,501,0,501);
  hpp_matrix_HLT = new TH2F(Form("hpp_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),501,0,501,501,0,501);
  hpp_anaBin_matrix_HLT = new TH2F(Form("hpp_anaBin_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix refpt jtpt from trigger addition R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
  hpp_mcclosure_matrix = new TH2F(Form("hpp_mcclosure_matrix_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt R%d %s ",radius,etaWidth),501,0,501,501,0,501);
  hpp_mcclosure_matrix_HLT = new TH2F(Form("hpp_mcclosure_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),501,0,501,501,0,501);
  //cout<<"C"<<endl;
  hpp_mcclosure_data = new TH1F(Form("hpp_mcclosure_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_JetComb_data = new TH1F(Form("hpp_mcclosure_JetComb_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger combined  R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet80_data = new TH1F(Form("hpp_mcclosure_Jet80_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 80  R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet60_data = new TH1F(Form("hpp_mcclosure_Jet60_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 60  R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet40_data = new TH1F(Form("hpp_mcclosure_Jet40_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 40  R%d %s ",radius,etaWidth),501,0,501);

  hpp_mcclosure_data_train = new TH1F(Form("hpp_mcclosure_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_JetComb_data_train = new TH1F(Form("hpp_mcclosure_JetComb_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger combined  R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet80_data_train = new TH1F(Form("hpp_mcclosure_Jet80_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger 80  R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet60_data_train = new TH1F(Form("hpp_mcclosure_Jet60_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger 60  R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet40_data_train = new TH1F(Form("hpp_mcclosure_Jet40_data_train_R%d_%s",radius,etaWidth),Form("data_train for unfolding mc closure test trigger 40  R%d %s ",radius,etaWidth),501,0,501);

  
  hpp_mcclosure_gen = new TH1F(Form("hpp_mcclosure_gen_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_JetComb_gen = new TH1F(Form("hpp_mcclosure_gen_JetComb_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger combined R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet80_gen = new TH1F(Form("hpp_mcclosure_gen_Jet80_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet60_gen = new TH1F(Form("hpp_mcclosure_gen_Jet60_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 60 R%d %s ",radius,etaWidth),501,0,501);
  hpp_mcclosure_Jet40_gen = new TH1F(Form("hpp_mcclosure_gen_Jet40_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 40 R%d %s ",radius,etaWidth),501,0,501);

  hpp_JetComb_gen = new TH1F(Form("hpp_JetComb_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from HLT trigger combined R%d %s ",radius,etaWidth),501,0,501);
  hpp_JetComb_reco = new TH1F(Form("hpp_JetComb_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from HLT trigger combined R%d %s ",radius,etaWidth),501,0,501);
  hpp_Jet80_gen = new TH1F(Form("hpp_Jet80_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet80 trigger R%d %s ",radius,etaWidth),501,0,501);
  hpp_Jet80_reco = new TH1F(Form("hpp_Jet80_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet80 trigger R%d %s ",radius,etaWidth),501,0,501);
  hpp_Jet60_gen = new TH1F(Form("hpp_Jet60_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),501,0,501);
  hpp_Jet60_reco = new TH1F(Form("hpp_Jet60_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),501,0,501);
  hpp_Jet40_gen = new TH1F(Form("hpp_Jet40_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),501,0,501);
  hpp_Jet40_reco = new TH1F(Form("hpp_Jet40_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),501,0,501);

  hpp_anaBin_JetComb_gen = new TH1F(Form("hpp_anaBin_JetComb_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from HLT trigger combined R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_JetComb_reco = new TH1F(Form("hpp_anaBin_JetComb_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from HLT trigger combined R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_Jet80_gen = new TH1F(Form("hpp_anaBin_Jet80_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_Jet80_reco = new TH1F(Form("hpp_anaBin_Jet80_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_Jet60_gen = new TH1F(Form("hpp_anaBin_Jet60_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_Jet60_reco = new TH1F(Form("hpp_anaBin_Jet60_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_Jet40_gen = new TH1F(Form("hpp_anaBin_Jet40_gen_R%d_%s",radius,etaWidth),Form("Gen refpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_Jet40_reco = new TH1F(Form("hpp_anaBin_Jet40_reco_R%d_%s",radius,etaWidth),Form("reco jtpt from Jet40 && !Jet60 && !Jet80 trigger R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);  


  // Define all the histograms necessary for the analysis: 
  
  // 1 - Data, 2 - MC
  Float_t pfpt_1[1000], pfpt_2[1000];
  Float_t pfrefpt_2[1000];
  Float_t pthat_2;
  Float_t calopt_1[1000], calopt_2[1000];
  Int_t npf_1, npf_2;
  Float_t eMax_1[1000], eMax_2[1000];
  Float_t chMax_1[1000], chMax_2[1000];
  Float_t chSum_1[1000], chSum_2[1000];
  Float_t phSum_1[1000], phSum_2[1000];
  Float_t neSum_1[1000], neSum_2[1000];
  Float_t muSum_1[1000], muSum_2[1000];
  Float_t phMax_1[1000], phMax_2[1000];
  Float_t neMax_1[1000], neMax_2[1000];
  Float_t muMax_1[1000], muMax_2[1000];
  Int_t jet40_1, jet60_1, jet80_1;
  Int_t jet40_p_1, jet60_p_1, jet80_p_1;
  Int_t jet40_2, jet60_2, jet80_2;
  Int_t jet40_p_2;
  Double_t weight;
  Float_t vz; 
  Int_t subid_2[1000];
  Float_t eta_1[1000], eta_2[1000];
  Int_t nref_2;
  Int_t pfrefidx_2[1000]; //! REF -> PF match index
  Float_t refdrjt_2[1000];

  Data_matched->SetBranchAddress("calopt",&calopt_1);
  Data_matched->SetBranchAddress("npf", &npf_1);
  Data_matched->SetBranchAddress("pfpt",&pfpt_1);
  Data_matched->SetBranchAddress("eMax",&eMax_1);
  Data_matched->SetBranchAddress("chMax",&chMax_1);
  Data_matched->SetBranchAddress("chSum",&chSum_1);
  Data_matched->SetBranchAddress("phSum",&phSum_1);
  Data_matched->SetBranchAddress("neSum",&neSum_1);
  Data_matched->SetBranchAddress("muSum",&muSum_1);
  Data_matched->SetBranchAddress("phMax",&phMax_1);
  Data_matched->SetBranchAddress("neMax",&neMax_1);
  Data_matched->SetBranchAddress("muMax",&muMax_1);
  Data_matched->SetBranchAddress("jet40",&jet40_1);
  Data_matched->SetBranchAddress("jet60",&jet60_1);
  Data_matched->SetBranchAddress("jet80",&jet80_1);
  Data_matched->SetBranchAddress("jet40_prescl",&jet40_p_1);
  Data_matched->SetBranchAddress("jet60_prescl",&jet60_p_1);
  Data_matched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  Data_matched->SetBranchAddress("pfeta",&eta_1);

#if 0
  Data_unmatched->SetBranchAddress("pfpt",&pfpt_1);
  Data_unmatched->SetBranchAddress("eMax",&eMax_1);
  Data_unmatched->SetBranchAddress("chMax",&chMax_1);
  Data_unmatched->SetBranchAddress("chSum",&chSum_1);
  Data_unmatched->SetBranchAddress("phSum",&phSum_1);
  Data_unmatched->SetBranchAddress("neSum",&neSum_1);
  Data_unmatched->SetBranchAddress("muSum",&muSum_1);
  Data_unmatched->SetBranchAddress("jet40",&jet40_1);
  Data_unmatched->SetBranchAddress("jet60",&jet60_1);
  Data_unmatched->SetBranchAddress("jet80",&jet80_1);
  Data_unmatched->SetBranchAddress("jet40_prescl",&jet40_p_1);
  Data_unmatched->SetBranchAddress("jet60_prescl",&jet60_p_1);
  Data_unmatched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  Data_unmatched->SetBranchAddress("pfeta",&eta_1);
#endif
  MC_matched->SetBranchAddress("npf", &npf_2);
  MC_matched->SetBranchAddress("nref", &nref_2);
  MC_matched->SetBranchAddress("pfrefidx",&pfrefidx_2);
  MC_matched->SetBranchAddress("refdrjt",&refdrjt_2);
  MC_matched->SetBranchAddress("pthat", &pthat_2);
  MC_matched->SetBranchAddress("calopt",&calopt_2);
  MC_matched->SetBranchAddress("pfpt",&pfpt_2);
  MC_matched->SetBranchAddress("eMax",&eMax_2);
  MC_matched->SetBranchAddress("vz",&vz);
  MC_matched->SetBranchAddress("chMax",&chMax_2);
  MC_matched->SetBranchAddress("phMax",&phMax_2);
  MC_matched->SetBranchAddress("neMax",&neMax_2);
  MC_matched->SetBranchAddress("muMax",&muMax_2);
  MC_matched->SetBranchAddress("chSum",&chSum_2);
  MC_matched->SetBranchAddress("phSum",&phSum_2);
  MC_matched->SetBranchAddress("neSum",&neSum_2);
  MC_matched->SetBranchAddress("muSum",&muSum_2);
  if(ntuple == "Pawan")  MC_matched->SetBranchAddress("refpt",&pfrefpt_2);
  //if(ntuple == "Raghav") MC_matched->SetBranchAddress("pfrefpt",&pfrefpt_2);  MC_matched->SetBranchAddress("jet40",&jet40_2);
  MC_matched->SetBranchAddress("jet60",&jet60_2);
  MC_matched->SetBranchAddress("jet80",&jet80_2);
  MC_matched->SetBranchAddress("weight", &weight);
  MC_matched->SetBranchAddress("subid", &subid_2);
  MC_matched->SetBranchAddress("jet40_prescl",&jet40_p_2);
  MC_matched->SetBranchAddress("pfeta",&eta_2);

#if 0
  MC_unmatched->SetBranchAddress("pfpt",&pfpt_2);
  MC_unmatched->SetBranchAddress("eMax",&eMax_2);
  MC_unmatched->SetBranchAddress("chMax",&chMax_2);
  MC_unmatched->SetBranchAddress("chSum",&chSum_2);
  MC_unmatched->SetBranchAddress("phSum",&phSum_2);
  MC_unmatched->SetBranchAddress("neSum",&neSum_2);
  MC_unmatched->SetBranchAddress("muSum",&muSum_2);
  if(ntuple == "Pawan")  MC_matched->SetBranchAddress("refpt",&pfrefpt_2);
  if(ntuple == "Raghav") MC_matched->SetBranchAddress("pfrefpt",&pfrefpt_2);
  MC_unmatched->SetBranchAddress("jet40",&jet40_2);
  MC_unmatched->SetBranchAddress("jet60",&jet60_2);
  MC_unmatched->SetBranchAddress("jet80",&jet80_2);
  MC_unmatched->SetBranchAddress("weight", & weight);
  MC_unmatched->SetBranchAddress("subid", &subid_2);
  MC_unmatched->SetBranchAddress("jet40_prescl",&jet40_p_2);
  MC_unmatched->SetBranchAddress("pfeta",&eta_2);
  MC_unmatched->SetBranchAddress("vz",&vz);
#endif
  // random value for smear systematics. value per jet = rnd.Gaus(0,1);
  TRandom rnd; 

  // data loop
  long entries = Data_matched->GetEntries();
  //entries = 1000;
  cout<<"matched Data ntuple "<<endl;

  Float_t Jet40_prescl = 9.275;
  
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    Data_matched->GetEntry(nentry);

    for(int g = 0; g<npf_1; ++g){

      //if(muMax_1[g]/(chMax_1[g]+phMax_1[g]+neMax_1[g]+muMax_1[g]+eMax_1[g]) > 0.975) continue;

      Float_t Sumcand = chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g];

      if(isSymm && TMath::Abs(eta_1[g]) > (Float_t)etaHigh/10) continue;       
      if(!isSymm && (TMath::Abs(eta_1[g]) < (Float_t)etaLow/10 || TMath::Abs(eta_1[g]) > (Float_t)etaHigh/10)) continue;

      if(pfpt_1[g] > 40 && pfpt_1[g] <= 180) 
	pfpt_1[g] = fResidual->Eval(pfpt_1[g]) * pfpt_1[g];
    
      if(jet40_1 == 1 && jet60_1==0 && jet80_1==0) {
      
	hData_Jet40_noCut->Fill(pfpt_1[g], Jet40_prescl);
	hpp_Data_Jet40_noCut->Fill(pfpt_1[g], Jet40_prescl);

	//if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) {
	hData_Jet40_CutA->Fill(pfpt_1[g], Jet40_prescl);
	hpp_TrgObj40->Fill(pfpt_1[g], Jet40_prescl);
	hpp_JEC_TrgObj40->Fill(pfpt_1[g]*1.005, Jet40_prescl);
	hpp_Smear_TrgObj40->Fill(pfpt_1[g]+ rnd.Gaus(0,1), Jet40_prescl);
	hpp_anaBin_TrgObj40->Fill(pfpt_1[g], Jet40_prescl);

	//}
	//if(calopt_1[g]/pfpt_1[g] > 0.85){
	// hData_Jet40_CutA->Fill(pfpt_1[g], Jet40_prescl);
	// hpp_TrgObj40->Fill(pfpt_1[g], Jet40_prescl);
	// //}
	// //if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1/Sumcand < 0.05) {
	// hData_Jet40_CutA->Fill(pfpt_1[g], Jet40_prescl);
	// hpp_TrgObj40->Fill(pfpt_1[g], Jet40_prescl);
	// //}

	// if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet40_CutA_rej->Fill(pfpt_1[g], Jet40_prescl);
	// if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) hData_Jet40_CutA_rej->Fill(pfpt_1[g], Jet40_prescl);
      
      }
    
      if(jet60_1 == 1 && jet80_1==0) {
      
	hData_Jet60_noCut->Fill(pfpt_1[g]);
	hpp_Data_Jet60_noCut->Fill(pfpt_1[g]);

	//if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)){
	hData_Jet60_CutA->Fill(pfpt_1[g]);
	hpp_TrgObj60->Fill(pfpt_1[g]);
	hpp_JEC_TrgObj60->Fill(pfpt_1[g]*1.005);
	hpp_Smear_TrgObj60->Fill(pfpt_1[g]+rnd.Gaus(0,1));
	hpp_anaBin_TrgObj60->Fill(pfpt_1[g]);
	//}
	// //if(calopt_1[g]/pfpt_1[g] > 0.85) {
	// hData_Jet60_CutA->Fill(pfpt_1[g]);
	// hpp_TrgObj60->Fill(pfpt_1[g]);
	// //}
	// //if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1/Sumcand < 0.05) {
	// hData_Jet60_CutA->Fill(pfpt_1[g]);
	// hpp_TrgObj60->Fill(pfpt_1[g]);
	// //}
	// if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet60_CutA_rej->Fill(pfpt_1[g]);
	// if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) hData_Jet60_CutA_rej->Fill(pfpt_1[g]);
      
      }

      if(jet80_1 == 1) {
    
	hData_Jet80_noCut->Fill(pfpt_1[g]);
	hpp_Data_Jet80_noCut->Fill(pfpt_1[g]);

	//if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) {
	hData_Jet80_CutA->Fill(pfpt_1[g]);
	hpp_TrgObj80->Fill(pfpt_1[g]);
	hpp_JEC_TrgObj80->Fill(pfpt_1[g]*1.005);
	hpp_Smear_TrgObj80->Fill(pfpt_1[g]+rnd.Gaus(0,1));
	hpp_anaBin_TrgObj80->Fill(pfpt_1[g]);
	//}
	// //if(calopt_1[g]/pfpt_1[g] > 0.85){
	// hData_Jet80_CutA->Fill(pfpt_1[g]);
	// hpp_TrgObj80->Fill(pfpt_1[g]);
	// //}
	// //if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1/Sumcand < 0.05){
	// hData_Jet80_CutA->Fill(pfpt_1[g]);
	// hpp_TrgObj80->Fill(pfpt_1[g]);
	// //}
	// if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet80_CutA_rej->Fill(pfpt_1[g]);
	// if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) hData_Jet80_CutA_rej->Fill(pfpt_1[g]);
      
      }

    }// jet ntuple 
    
  }// data ntuple loop

#if 0
  // data unmatched loop:
  entries = Data_unmatched->GetEntries();
  //entries = 1000;
  cout<<"Unmatched Data ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    Data_unmatched->GetEntry(nentry);

    if(pfpt_1 > 40 && pfpt_1 <= 180) 
      pfpt_1 = fResidual->Eval(pfpt_1) * pfpt_1;

    Float_t Sumcand = chSum_1 + phSum_1 + neSum_1 + muSum_1;

    if(isSymm && TMath::Abs(eta_1) > (Float_t)etaHigh/10) continue;       
    if(!isSymm && (TMath::Abs(eta_1) < (Float_t)etaLow/10 || TMath::Abs(eta_1) > (Float_t)etaHigh/10)) continue;
    
    if(jet40_1 == 1 && jet60_1 == 0 && jet80_1 == 0) {
    
      hData_unmatched_Jet40_noCut->Fill(pfpt_1, Jet40_prescl);
      hpp_Data_Jet40_noCut->Fill(pfpt_1, Jet40_prescl);
      //if(eMax_1/Sumcand < 0.05 ){
      hpp_TrgObj40->Fill(pfpt_1, Jet40_prescl);
      hpp_JEC_TrgObj40->Fill(pfpt_1*1.005, Jet40_prescl);
      hpp_Smear_TrgObj40->Fill(pfpt_1+rnd.Gaus(0,1), Jet40_prescl);
      hpp_anaBin_TrgObj40->Fill(pfpt_1, Jet40_prescl);
      hData_unmatched_Jet40_CutA->Fill(pfpt_1, Jet40_prescl);
      //}else hData_unmatched_Jet40_CutA_rej->Fill(pfpt_1, Jet40_prescl);
      
    }

    if(jet60_1 == 1 && jet80_1==0) {

      hData_unmatched_Jet60_noCut->Fill(pfpt_1);
      hpp_Data_Jet60_noCut->Fill(pfpt_1);
      
      //if(eMax_1/Sumcand < 0.05  ){
	hpp_TrgObj60->Fill(pfpt_1);
	hpp_JEC_TrgObj60->Fill(pfpt_1*1.005);
	hpp_Smear_TrgObj60->Fill(pfpt_1+rnd.Gaus(0,1));
	hpp_anaBin_TrgObj60->Fill(pfpt_1);
	hData_unmatched_Jet60_CutA->Fill(pfpt_1);
	//}else hData_unmatched_Jet60_CutA_rej->Fill(pfpt_1);
      
    }

    if(jet80_1 == 1) {
    
      hData_unmatched_Jet80_noCut->Fill(pfpt_1);
      hpp_Data_Jet80_noCut->Fill(pfpt_1);
      //if(eMax_1/Sumcand < 0.05  ){
      hpp_TrgObj80->Fill(pfpt_1);
      hpp_JEC_TrgObj80->Fill(pfpt_1*1.005);
      hpp_Smear_TrgObj80->Fill(pfpt_1+rnd.Gaus(0,1));
      hpp_anaBin_TrgObj80->Fill(pfpt_1);
      hData_unmatched_Jet80_CutA->Fill(pfpt_1);
      //}else hData_unmatched_Jet80_CutA_rej->Fill(pfpt_1);
      
    }
    
  }// data unmatched ntuple loop
#endif

  entries = MC_matched->GetEntries();
  //entries = 1000;
  // MC loop
  cout<<" looking at matched MC ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    MC_matched->GetEntry(nentry);

    for(int g = 0; g<npf_2; ++g){
    
      Float_t Sumcand = chSum_2[g] + phSum_2[g] + neSum_2[g] + muSum_2[g];
      int refid = -1; 
      
      if(pfpt_2[g] > 2 * pthat_2) continue;
      refid = pfrefidx_2[g];
      if(subid_2[refid] != 0 || fabs(refdrjt_2[refid]) > kdelrcut) continue;
      if(refid < 0) continue;
      //if(muMax_1[g]/(chMax_2[g]+phMax_2[g]+neMax_2[g]+muMax_2[g]+eMax_2[g]) > 0.975) continue;

      if(isSymm && TMath::Abs(eta_2[g]) > (Float_t)etaHigh/10) continue;       
      if(!isSymm && (TMath::Abs(eta_2[g]) < (Float_t)etaLow/10 || TMath::Abs(eta_2[g]) > (Float_t)etaHigh/10)) continue;
    
      //weight = (Float_t)weight * fVzPP->Eval(vz) * 1e-5;

      hpp_gen->Fill(pfrefpt_2[refid], weight);
      hpp_reco->Fill(pfpt_2[g], weight);
      hpp_matrix->Fill(pfrefpt_2[refid], pfpt_2[g], weight);

      if(nentry%2 == 0){
	hpp_mcclosure_data->Fill(pfpt_2[g], weight);
      }else {
	hpp_mcclosure_matrix->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	hpp_mcclosure_gen->Fill(pfrefpt_2[refid], weight);
	hpp_mcclosure_data_train->Fill(pfpt_2[g], weight);
      }
      
      if(jet40_2 == 1 && jet60_2==0 && jet80_2 == 0){
      
	hMC_Jet40_noCut->Fill(pfrefpt_2[refid], weight);
	hpp_MC_Jet40_noCut->Fill(pfrefpt_2[refid], weight);
	//if(calopt_2[g]/pfpt_2[g] > 0.5 && calopt_2[g]/pfpt_2[g] <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2[g]/pfpt_2[g] - (Float_t)9/7)){
	hMC_Jet40_CutA->Fill(pfrefpt_2[refid], weight);
	hpp_Jet40_gen->Fill(pfrefpt_2[refid], weight);
	hpp_Jet40_reco->Fill(pfpt_2[g], weight);
	hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	
	hpp_anaBin_Jet40_gen->Fill(pfrefpt_2[refid], weight);
	hpp_anaBin_Jet40_reco->Fill(pfpt_2[g], weight);
	hpp_anaBin_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	
	if(nentry%2==0) {
	  hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	  hpp_mcclosure_Jet40_gen->Fill(pfrefpt_2[refid], weight);
	  hpp_mcclosure_Jet40_data_train->Fill(pfpt_2[g], weight);
	}
	if(nentry%2==1) {
	  hpp_mcclosure_Jet40_data->Fill(pfpt_2[g], weight);
	}
	
	// //}
	// //if(calopt_2[g]/pfpt_2[g] > 0.85) {
	// hMC_Jet40_CutA->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet40_gen->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet40_reco->Fill(pfpt_2[g], weight);
	// hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);

	// if(nentry%2==0) {
	//   hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	//   hpp_mcclosure_Jet40_gen->Fill(pfrefpt_2[refid], weight);
	// }
	// if(nentry%2==1) {
	//   hpp_mcclosure_Jet40_data->Fill(pfpt_2[g], weight);
	// }
	
	// //}
	// //if(calopt_2[g]/pfpt_2[g] <= 0.5 && eMax_2/Sumcand < 0.05) {
	// hMC_Jet40_CutA->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet40_gen->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet40_reco->Fill(pfpt_2[g], weight);
	// hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);

	// if(nentry%2==0) {
	//   hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	//   hpp_mcclosure_Jet40_gen->Fill(pfrefpt_2[refid], weight);
	// }
	// if(nentry%2==1) {
	//   hpp_mcclosure_Jet40_data->Fill(pfpt_2[g], weight);
	// }
	
	//}
	// if(calopt_2[g]/pfpt_2[g] <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet40_CutA_rej->Fill(pfrefpt_2[refid], weight);
	// if(calopt_2[g]/pfpt_2[g] > 0.5 && calopt_2[g]/pfpt_2[g] <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2[g]/pfpt_2[g] - (Float_t)9/7)) hMC_Jet40_CutA_rej->Fill(pfrefpt_2[refid], weight);
      
      }
    
      if(jet60_2 == 1 && jet80_2 == 0){
      

	hMC_Jet60_noCut->Fill(pfrefpt_2[refid], weight);
	hpp_MC_Jet60_noCut->Fill(pfrefpt_2[refid], weight);

	//if(calopt_2[g]/pfpt_2[g] > 0.5 && calopt_2[g]/pfpt_2[g] <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2[g]/pfpt_2[g] - (Float_t)9/7)) {
	hMC_Jet60_CutA->Fill(pfrefpt_2[refid], weight);
	hpp_Jet60_gen->Fill(pfrefpt_2[refid], weight);
	hpp_Jet60_reco->Fill(pfpt_2[g], weight);
	hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);

	hpp_anaBin_Jet60_gen->Fill(pfrefpt_2[refid], weight);
	hpp_anaBin_Jet60_reco->Fill(pfpt_2[g], weight);
	hpp_anaBin_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	
	if(nentry%2==0) {
	  hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	  hpp_mcclosure_Jet60_gen->Fill(pfrefpt_2[refid], weight);
	  hpp_mcclosure_Jet60_data_train->Fill(pfpt_2[g], weight);
	}
	if(nentry%2==1) {
	  hpp_mcclosure_Jet60_data->Fill(pfpt_2[g], weight);
	}
	
	// //}
	// //if(calopt_2[g]/pfpt_2[g] > 0.85) {
	// hMC_Jet60_CutA->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet60_gen->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet60_reco->Fill(pfpt_2[g], weight);
	// hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);

	// if(nentry%2==0) {
	//   hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	//   hpp_mcclosure_Jet60_gen->Fill(pfrefpt_2[refid], weight);
	// }
	// if(nentry%2==1) {
	//   hpp_mcclosure_Jet60_data->Fill(pfpt_2[g], weight);
	// }
	// //}
	// //if(calopt_2[g]/pfpt_2[g] <= 0.5 && eMax_2/Sumcand < 0.05) {
	// hMC_Jet60_CutA->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet60_gen->Fill(pfrefpt_2[refid], weight);
	// hpp_Jet60_reco->Fill(pfpt_2[g], weight);
	// hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);

	// if(nentry%2==0) {
	//   hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	//   hpp_mcclosure_Jet60_gen->Fill(pfrefpt_2[refid], weight);
	// }
	// if(nentry%2==1) {
	//   hpp_mcclosure_Jet60_data->Fill(pfpt_2[g], weight);
	// }
	// //}

	// if(calopt_2[g]/pfpt_2[g] <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet60_CutA_rej->Fill(pfrefpt_2[refid], weight);
	// if(calopt_2[g]/pfpt_2[g] > 0.5 && calopt_2[g]/pfpt_2[g] <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2[g]/pfpt_2[g] - (Float_t)9/7)) hMC_Jet60_CutA_rej->Fill(pfrefpt_2[refid], weight);

      }

    
      if(jet80_2 == 1){

	hMC_Jet80_noCut->Fill(pfrefpt_2[refid], weight);
	hpp_MC_Jet80_noCut->Fill(pfrefpt_2[refid], weight);

	//if(calopt_2[g]/pfpt_2[g] > 0.5 && calopt_2[g]/pfpt_2[g] <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2[g]/pfpt_2[g] - (Float_t)9/7)) {
	hMC_Jet80_CutA->Fill(pfrefpt_2[refid], weight);
	hpp_Jet80_gen->Fill(pfrefpt_2[refid], weight);
	hpp_Jet80_reco->Fill(pfpt_2[g], weight);
	hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);

	hpp_anaBin_Jet80_gen->Fill(pfrefpt_2[refid], weight);
	hpp_anaBin_Jet80_reco->Fill(pfpt_2[g], weight);
	hpp_anaBin_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	
	if(nentry%2==0) {
	  hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight);
	  hpp_mcclosure_Jet80_gen->Fill(pfrefpt_2[refid], weight);
	  hpp_mcclosure_Jet80_data_train->Fill(pfpt_2[g], weight);
	}
	if(nentry%2==1) {
	  hpp_mcclosure_Jet80_data->Fill(pfpt_2[g], weight);
	}
	//}
	// //if(calopt_2[g]/pfpt_2[g] > 0.85) {
	// hMC_Jet80_CutA->Fill(pfrefpt_2[refid], weight[g]);
	// hpp_Jet80_gen->Fill(pfrefpt_2[refid], weight[g]);
	// hpp_Jet80_reco->Fill(pfpt_2[g], weight[g]);
	// hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight[g]);

	// if(nentry%2==0) {
	//   hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight[g]);
	//   hpp_mcclosure_Jet80_gen->Fill(pfrefpt_2[refid], weight[g]);
	// }
	// if(nentry%2==1) {
	//   hpp_mcclosure_Jet80_data->Fill(pfpt_2[g], weight[g]);
	// }
	// //}
	// //if(calopt_2[g]/pfpt_2[g] <= 0.5 && eMax_2/Sumcand < 0.05) {
	// hMC_Jet80_CutA->Fill(pfrefpt_2[refid], weight[g]);
	// hpp_Jet80_gen->Fill(pfrefpt_2[refid], weight[g]);
	// hpp_Jet80_reco->Fill(pfpt_2[g], weight[g]);
	// hpp_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight[g]);

	// if(nentry%2==0) {
	//   hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2[refid], pfpt_2[g], weight[g]);
	//   hpp_mcclosure_Jet80_gen->Fill(pfrefpt_2[refid], weight[g]);
	// }
	// if(nentry%2==1) {
	//   hpp_mcclosure_Jet80_data->Fill(pfpt_2[g], weight[g]);
	// }
	//}
      
	// if(calopt_2[g]/pfpt_2[g] <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet80_CutA_rej->Fill(pfrefpt_2[refid], weight[g]);
	// if(calopt_2[g]/pfpt_2[g] > 0.5 && calopt_2[g]/pfpt_2[g] <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2[g]/pfpt_2[g] - (Float_t)9/7)) hMC_Jet80_CutA_rej->Fill(pfrefpt_2[refid], weight[g]);

      }//jet80 selection

    }//jet loop
    
    
  }// mc ntuple loop

#if 0
  entries = MC_unmatched->GetEntries();
  //entries = 1000;
  // MC loop
  cout<<" looking at unmatched MC ntuple"<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    MC_unmatched->GetEntry(nentry);

    if(subid_2 != 0) continue;

    Float_t Sumcand = chSum_2 + phSum_2 + neSum_2 + muSum_2;

    if(isSymm && TMath::Abs(eta_2) > (Float_t)etaHigh/10) continue;       
    if(!isSymm && (TMath::Abs(eta_2) < (Float_t)etaLow/10 || TMath::Abs(eta_2) > (Float_t)etaHigh/10)) continue;
    
    weight = (Float_t)weight * fVzPP->Eval(vz) * 1e-5;

    hpp_gen->Fill(pfrefpt_2, weight);
    hpp_reco->Fill(pfpt_2, weight);
    hpp_matrix->Fill(pfrefpt_2, pfpt_2, weight);
    
    if(jet40_2 == 1 && jet60_2==0 && jet80_2==0){
      
      hMC_unmatched_Jet40_noCut->Fill(pfrefpt_2, weight);
      hpp_MC_Jet40_noCut->Fill(pfrefpt_2, weight);

      //if(eMax_2/Sumcand < 0.05  ){
	hpp_Jet40_gen->Fill(pfrefpt_2, weight);
	hpp_Jet40_reco->Fill(pfpt_2, weight);
	hpp_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);

	hpp_anaBin_Jet40_gen->Fill(pfrefpt_2, weight);
	hpp_anaBin_Jet40_reco->Fill(pfpt_2, weight);
	hpp_anaBin_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);
	
	if(nentry%2==0) {
	  hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);
	  hpp_mcclosure_Jet40_gen->Fill(pfrefpt_2, weight);
	}
	if(nentry%2==1) {
	  hpp_mcclosure_Jet40_data->Fill(pfpt_2, weight);
	}
      
	hMC_unmatched_Jet40_CutA->Fill(pfrefpt_2, weight);
	//}else hMC_unmatched_Jet40_CutA_rej->Fill(pfrefpt_2, weight);
      
    }

    
    if(jet60_2 == 1 && jet80_2==0){

      hpp_MC_Jet60_noCut->Fill(pfrefpt_2, weight);
      hMC_unmatched_Jet60_noCut->Fill(pfrefpt_2, weight);
      //if(eMax_2/Sumcand < 0.05  ){
	hpp_Jet60_gen->Fill(pfrefpt_2, weight);
	hpp_Jet60_reco->Fill(pfpt_2, weight);
	hpp_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);

	hpp_anaBin_Jet60_gen->Fill(pfrefpt_2, weight);
	hpp_anaBin_Jet60_reco->Fill(pfpt_2, weight);
	hpp_anaBin_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);
	
	if(nentry%2==0) {
	  hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);
	  hpp_mcclosure_Jet60_gen->Fill(pfrefpt_2, weight);
	}
	if(nentry%2==1) {
	  hpp_mcclosure_Jet60_data->Fill(pfpt_2, weight);
	}
      
	hMC_unmatched_Jet60_CutA->Fill(pfrefpt_2, weight);
	//}else hMC_unmatched_Jet60_CutA_rej->Fill(pfrefpt_2);
      
    }

    
    if(jet80_2 == 1){

      hpp_MC_Jet80_noCut->Fill(pfrefpt_2, weight);
      hMC_unmatched_Jet80_noCut->Fill(pfrefpt_2, weight);
      //if(eMax_2/Sumcand < 0.05  ){
	hpp_Jet80_gen->Fill(pfrefpt_2, weight);
	hpp_Jet80_reco->Fill(pfpt_2, weight);
	hpp_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);

	hpp_anaBin_Jet80_gen->Fill(pfrefpt_2, weight);
	hpp_anaBin_Jet80_reco->Fill(pfpt_2, weight);
	hpp_anaBin_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);
	
	if(nentry%2==0) {
	  hpp_mcclosure_matrix_HLT->Fill(pfrefpt_2, pfpt_2, weight);
	  hpp_mcclosure_Jet80_gen->Fill(pfrefpt_2, weight);
	}
	if(nentry%2==1) {
	  hpp_mcclosure_Jet80_data->Fill(pfpt_2, weight);
	}
      
	hMC_unmatched_Jet80_CutA->Fill(pfrefpt_2, weight);
	//}else hMC_unmatched_Jet80_CutA_rej->Fill(pfrefpt_2); 
    }
    
  }// mc unmatched  ntuple loop
#endif

  TFile fout(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/%s_TTree_PP_Data_MC_spectra_residualFactor_%s_%s_R0p%d.root",ntuple, ptbins, etaWidth, radius),"RECREATE");
  fout.cd();
  
  hpp_TrgObjComb->Add(hpp_TrgObj80);
  hpp_TrgObjComb->Add(hpp_TrgObj60);
  hpp_TrgObjComb->Add(hpp_TrgObj40);

  divideBinWidth(hpp_TrgObjComb);
  divideBinWidth(hpp_TrgObj80);
  divideBinWidth(hpp_TrgObj60);
  divideBinWidth(hpp_TrgObj40);
  
  hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj80);
  hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj60);
  hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj40);

  divideBinWidth(hpp_JEC_TrgObjComb);
  divideBinWidth(hpp_JEC_TrgObj80);
  divideBinWidth(hpp_JEC_TrgObj60);
  divideBinWidth(hpp_JEC_TrgObj40);

  hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj80);
  hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj60);
  hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj40);

  divideBinWidth(hpp_Smear_TrgObjComb);
  divideBinWidth(hpp_Smear_TrgObj80);
  divideBinWidth(hpp_Smear_TrgObj60);
  divideBinWidth(hpp_Smear_TrgObj40);


  hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj80);
  hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj60);
  hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj40);

  divideBinWidth(hpp_anaBin_TrgObjComb);
  divideBinWidth(hpp_anaBin_TrgObj80);
  divideBinWidth(hpp_anaBin_TrgObj60);
  divideBinWidth(hpp_anaBin_TrgObj40);

  hpp_MC_Comb_noCut->Add(hpp_MC_Jet80_noCut);
  hpp_MC_Comb_noCut->Add(hpp_MC_Jet60_noCut);
  hpp_MC_Comb_noCut->Add(hpp_MC_Jet40_noCut);
  
  hpp_Data_Comb_noCut->Add(hpp_Data_Jet80_noCut);
  hpp_Data_Comb_noCut->Add(hpp_Data_Jet60_noCut);
  hpp_Data_Comb_noCut->Add(hpp_Data_Jet40_noCut);
  
  hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet80_data);
  hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet60_data);
  hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet40_data);

  divideBinWidth(hpp_mcclosure_JetComb_data);

  hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet80_gen);
  hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet60_gen);
  hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet40_gen);

  divideBinWidth(hpp_mcclosure_JetComb_gen);
  
  hpp_JetComb_reco->Add(hpp_Jet80_reco);
  hpp_JetComb_reco->Add(hpp_Jet60_reco);
  hpp_JetComb_reco->Add(hpp_Jet40_reco);

  divideBinWidth(hpp_JetComb_reco);
  
  hpp_JetComb_gen->Add(hpp_Jet80_gen);
  hpp_JetComb_gen->Add(hpp_Jet60_gen);
  hpp_JetComb_gen->Add(hpp_Jet40_gen);

  divideBinWidth(hpp_JetComb_gen);

  hpp_anaBin_JetComb_reco->Add(hpp_anaBin_Jet80_reco);
  hpp_anaBin_JetComb_reco->Add(hpp_anaBin_Jet60_reco);
  hpp_anaBin_JetComb_reco->Add(hpp_anaBin_Jet40_reco);

  divideBinWidth(hpp_anaBin_JetComb_reco);
  
  hpp_anaBin_JetComb_gen->Add(hpp_anaBin_Jet80_gen);
  hpp_anaBin_JetComb_gen->Add(hpp_anaBin_Jet60_gen);
  hpp_anaBin_JetComb_gen->Add(hpp_anaBin_Jet40_gen);

  divideBinWidth(hpp_anaBin_JetComb_gen);

  hpp_MC_Comb_noCut->Write();
  hpp_Data_Comb_noCut->Write();

  hpp_reco->Write();
  hpp_gen->Write();
  hpp_matrix->Write();
  
  hpp_TrgObjComb->Write();
  hpp_TrgObj80->Write();
  hpp_TrgObj60->Write();
  hpp_TrgObj40->Write();

  hpp_JEC_TrgObjComb->Write();
  hpp_JEC_TrgObj80->Write();
  hpp_JEC_TrgObj60->Write();
  hpp_JEC_TrgObj40->Write();

  hpp_Smear_TrgObjComb->Write();
  hpp_Smear_TrgObj80->Write();
  hpp_Smear_TrgObj60->Write();
  hpp_Smear_TrgObj40->Write();

  hpp_matrix_HLT->Write();
  hpp_mcclosure_matrix_HLT->Write();
  hpp_mcclosure_JetComb_data->Write();
  hpp_mcclosure_Jet80_data->Write();
  hpp_mcclosure_Jet60_data->Write();
  hpp_mcclosure_Jet40_data->Write();
  hpp_mcclosure_JetComb_gen->Write();    
  hpp_mcclosure_Jet80_gen->Write();
  hpp_mcclosure_Jet60_gen->Write();
  hpp_mcclosure_Jet40_gen->Write();

  hpp_JetComb_reco->Write();
  hpp_Jet80_reco->Write();
  hpp_Jet60_reco->Write();
  hpp_Jet40_reco->Write();
  hpp_JetComb_gen->Write();
  hpp_Jet80_gen->Write();
  hpp_Jet60_gen->Write();
  hpp_Jet40_gen->Write();

  hpp_anaBin_TrgObjComb->Write();
  hpp_anaBin_TrgObj80->Write();
  hpp_anaBin_TrgObj60->Write();
  hpp_anaBin_TrgObj40->Write();
  hpp_anaBin_matrix_HLT->Write();
  hpp_anaBin_JetComb_reco->Write();
  hpp_anaBin_Jet80_reco->Write();
  hpp_anaBin_Jet60_reco->Write();
  hpp_anaBin_Jet40_reco->Write();
  hpp_anaBin_JetComb_gen->Write();
  hpp_anaBin_Jet80_gen->Write();
  hpp_anaBin_Jet60_gen->Write();
  hpp_anaBin_Jet40_gen->Write();
 

  // add the unmatched histograms to the matched ones to get the final cut efficiency
  hData_Jet40_noCut->Add(hData_unmatched_Jet40_noCut);
  hData_Jet60_noCut->Add(hData_unmatched_Jet60_noCut);
  hData_Jet80_noCut->Add(hData_unmatched_Jet80_noCut);
  
  hData_Jet40_CutA->Add(hData_unmatched_Jet40_CutA);
  hData_Jet60_CutA->Add(hData_unmatched_Jet60_CutA);
  hData_Jet80_CutA->Add(hData_unmatched_Jet80_CutA);

  hData_Jet40_CutA_rej->Add(hData_unmatched_Jet40_CutA_rej);
  hData_Jet60_CutA_rej->Add(hData_unmatched_Jet60_CutA_rej);
  hData_Jet80_CutA_rej->Add(hData_unmatched_Jet80_CutA_rej);

  hMC_Jet40_noCut->Add(hMC_unmatched_Jet40_noCut);
  hMC_Jet60_noCut->Add(hMC_unmatched_Jet60_noCut);
  hMC_Jet80_noCut->Add(hMC_unmatched_Jet80_noCut);
  
  hMC_Jet40_CutA->Add(hMC_unmatched_Jet40_CutA);
  hMC_Jet60_CutA->Add(hMC_unmatched_Jet60_CutA);
  hMC_Jet80_CutA->Add(hMC_unmatched_Jet80_CutA);

  hMC_Jet40_CutA_rej->Add(hMC_unmatched_Jet40_CutA_rej);
  hMC_Jet60_CutA_rej->Add(hMC_unmatched_Jet60_CutA_rej);
  hMC_Jet80_CutA_rej->Add(hMC_unmatched_Jet80_CutA_rej);

  

  hData_Jet60_noCut->Write();
  hData_Jet60_CutA->Write();
  // hData_Jet60_CutB->Write();
  TH1F * hData_Jet60_CutA_eff = (TH1F*)hData_Jet60_CutA->Clone("hData_Jet60_CutA_eff");
  hData_Jet60_CutA_eff->Divide(hData_Jet60_noCut);
  hData_Jet60_CutA_eff->Write();
  // TH1F * hData_Jet60_CutB_eff = (TH1F*)hData_Jet60_CutB->Clone("hData_Jet60_CutB_eff");
  // hData_Jet60_CutB_eff->Divide(hData_Jet60_noCut);
  // hData_Jet60_CutB_eff->Write();
  hData_Jet60_CutA_rej->Write();
  // hData_Jet60_CutB_rej->Write();
  TH1F * hData_Jet60_CutA_rej_eff = (TH1F*)hData_Jet60_CutA_rej->Clone("hData_Jet60_CutA_rej_eff");
  hData_Jet60_CutA_rej_eff->Divide(hData_Jet60_noCut);
  hData_Jet60_CutA_rej_eff->Write();
  // hData_Jet60_CutB_rej_eff = (TH1F*)hData_Jet60_CutB_rej->Clone("hData_Jet60_CutB_rej_eff");
  // hData_Jet60_CutB_rej_eff->Divide(hData_Jet60_noCut);
  // hData_Jet60_CutB_rej_eff->Write();

  hData_Jet40_noCut->Write();
  hData_Jet40_CutA->Write();
  // hData_Jet40_CutB->Write();
  TH1F * hData_Jet40_CutA_eff = (TH1F*)hData_Jet40_CutA->Clone("hData_Jet40_CutA_eff");
  hData_Jet40_CutA_eff->Divide(hData_Jet40_noCut);
  hData_Jet40_CutA_eff->Write();
  // TH1F * hData_Jet40_CutB_eff = (TH1F*)hData_Jet40_CutB->Clone("hData_Jet40_CutB_eff");
  // hData_Jet40_CutB_eff->Divide(hData_Jet40_noCut);
  // hData_Jet40_CutB_eff->Write();
  hData_Jet40_CutA_rej->Write();
  // hData_Jet40_CutB_rej->Write();
  TH1F * hData_Jet40_CutA_rej_eff = (TH1F*)hData_Jet40_CutA_rej->Clone("hData_Jet40_CutA_rej_eff");
  hData_Jet40_CutA_rej_eff->Divide(hData_Jet40_noCut);
  hData_Jet40_CutA_rej_eff->Write();
  // hData_Jet40_CutB_rej_eff = (TH1F*)hData_Jet40_CutB_rej->Clone("hData_Jet40_CutB_rej_eff");
  // hData_Jet40_CutB_rej_eff->Divide(hData_Jet40_noCut);
  // hData_Jet40_CutB_rej_eff->Write();

  hData_Jet80_noCut->Write();
  hData_Jet80_CutA->Write();
  //hData_Jet80_CutB->Write();
  TH1F * hData_Jet80_CutA_eff = (TH1F*)hData_Jet80_CutA->Clone("hData_Jet80_CutA_eff");
  hData_Jet80_CutA_eff->Divide(hData_Jet80_noCut);
  hData_Jet80_CutA_eff->Write();
  // TH1F * hData_Jet80_CutB_eff = (TH1F*)hData_Jet80_CutB->Clone("hData_Jet80_CutB_eff");
  // hData_Jet80_CutB_eff->Divide(hData_Jet80_noCut);
  // hData_Jet80_CutB_eff->Write();
  hData_Jet80_CutA_rej->Write();
  // hData_Jet80_CutB_rej->Write();
  TH1F * hData_Jet80_CutA_rej_eff = (TH1F*)hData_Jet80_CutA_rej->Clone("hData_Jet80_CutA_rej_eff");
  hData_Jet80_CutA_rej_eff->Divide(hData_Jet80_noCut);
  hData_Jet80_CutA_rej_eff->Write();
  // hData_Jet80_CutB_rej_eff = (TH1F*)hData_Jet80_CutB_rej->Clone("hData_Jet80_CutB_rej_eff");
  // hData_Jet80_CutB_rej_eff->Divide(hData_Jet80_noCut);
  // hData_Jet80_CutB_rej_eff->Write();
  
  
  hMC_Jet80_noCut->Write();
  hMC_Jet80_CutA->Write();
  // hMC_Jet80_CutB->Write();
  TH1F * hMC_Jet80_CutA_eff = (TH1F*)hMC_Jet80_CutA->Clone("hMC_Jet80_CutA_eff");
  hMC_Jet80_CutA_eff->Divide(hMC_Jet80_noCut);
  hMC_Jet80_CutA_eff->Write();
  // TH1F * hMC_Jet80_CutB_eff = (TH1F*)hMC_Jet80_CutB->Clone("hMC_Jet80_CutB_eff");
  // hMC_Jet80_CutB_eff->Divide(hMC_Jet80_noCut);
  // hMC_Jet80_CutB_eff->Write();
  
  hMC_Jet60_noCut->Write();
  hMC_Jet60_CutA->Write();
  // hMC_Jet60_CutB->Write();
  TH1F * hMC_Jet60_CutA_eff = (TH1F*)hMC_Jet60_CutA->Clone("hMC_Jet60_CutA_eff");
  hMC_Jet60_CutA_eff->Divide(hMC_Jet60_noCut);
  hMC_Jet60_CutA_eff->Write();
  // TH1F * hMC_Jet60_CutB_eff = (TH1F*)hMC_Jet60_CutB->Clone("hMC_Jet60_CutB_eff");
  // hMC_Jet60_CutB_eff->Divide(hMC_Jet60_noCut);
  // hMC_Jet60_CutB_eff->Write();
  
  hMC_Jet40_noCut->Write();
  hMC_Jet40_CutA->Write();
  // hMC_Jet40_CutB->Write();
  TH1F * hMC_Jet40_CutA_eff = (TH1F*)hMC_Jet40_CutA->Clone("hMC_Jet40_CutA_eff");
  hMC_Jet40_CutA_eff->Divide(hMC_Jet40_noCut);
  hMC_Jet40_CutA_eff->Write();
  // TH1F * hMC_Jet40_CutB_eff = (TH1F*)hMC_Jet40_CutB->Clone("hMC_Jet40_CutB_eff");
  // hMC_Jet40_CutB_eff->Divide(hMC_Jet40_noCut);
  // hMC_Jet40_CutB_eff->Write();

  // save the unmatched histograms as well:
  hData_unmatched_Jet80_noCut->Write();
  hData_unmatched_Jet80_CutA->Write();
  hData_unmatched_Jet80_CutA_rej->Write();
  hData_unmatched_Jet60_noCut->Write();
  hData_unmatched_Jet60_CutA->Write();
  hData_unmatched_Jet60_CutA_rej->Write();
  hData_unmatched_Jet40_noCut->Write();
  hData_unmatched_Jet40_CutA->Write();
  hData_unmatched_Jet40_CutA_rej->Write();

  hMC_unmatched_Jet80_noCut->Write();
  hMC_unmatched_Jet80_CutA->Write();
  hMC_unmatched_Jet80_CutA_rej->Write();
  hMC_unmatched_Jet60_noCut->Write();
  hMC_unmatched_Jet60_CutA->Write();
  hMC_unmatched_Jet60_CutA_rej->Write();
  hMC_unmatched_Jet40_noCut->Write();
  hMC_unmatched_Jet40_CutA->Write();
  hMC_unmatched_Jet40_CutA_rej->Write();
  
  
  TCanvas * cJet80_CutEfficiency_Jet80 = new TCanvas("cJet80_CutEfficiency_Jet80","",1000,800);
  cJet80_CutEfficiency_Jet80->Divide(2,1);
  cJet80_CutEfficiency_Jet80->cd(1);

  hData_Jet80_CutA_eff->Rebin(5);
  hData_Jet80_CutA_eff->Scale(1./5);
  hData_Jet80_CutA_eff->SetMarkerColor(kRed);
  hData_Jet80_CutA_eff->SetMarkerStyle(24);
  hData_Jet80_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet80_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet80_CutA_eff->SetXTitle(Form("ak%dPF p_{T}",radius));
  hData_Jet80_CutA_eff->SetTitle("Data");
  hData_Jet80_CutA_eff->SetYTitle("Jet80_Cut efficiency");
  hData_Jet80_CutA_eff->Draw();
  // hData_Jet80_CutB_eff->Rebin(20);
  // hData_Jet80_CutB_eff->Scale(1./20);
  // hData_Jet80_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet80_CutB_eff->SetMarkerStyle(33);
  // hData_Jet80_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet80_CutB_eff->Draw("same");

  cJet80_CutEfficiency_Jet80->cd(2);
  hMC_Jet80_CutA_eff->Rebin(5);
  hMC_Jet80_CutA_eff->Scale(1./5);
  hMC_Jet80_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet80_CutA_eff->SetMarkerStyle(24);
  hMC_Jet80_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet80_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet80_CutA_eff->SetXTitle(Form("ak%dPF ref p_{T}",radius));
  hMC_Jet80_CutA_eff->SetTitle("MC");
  hMC_Jet80_CutA_eff->SetYTitle("Jet80_Cut efficiency");
  hMC_Jet80_CutA_eff->Draw();
  // hMC_Jet80_CutB_eff->Rebin(20);
  // hMC_Jet80_CutB_eff->Scale(1./20);
  // hMC_Jet80_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet80_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet80_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet80_CutB_eff->Draw("same");

  cJet80_CutEfficiency_Jet80->SaveAs(Form("Pp_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_R0p%d_zoomed.pdf",radius),"RECREATE");

  // TCanvas * cCutRejection_Jet80 = new TCanvas("cCutRejection_Jet80","",1000,800);

  // hData_Jet80_CutA_rej->Rebin(5);
  // hData_Jet80_CutA_rej->Scale(1./5);
  // hData_Jet80_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet80_CutA_rej->SetMarkerStyle(24);
  // hData_Jet80_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet80_CutA_rej->SetAxisRange(0,1.2,"Y");
  // hData_Jet80_CutA_rej->SetXTitle("matched ak3PF p_{T}");
  // hData_Jet80_CutA_rej->SetTitle("Data");
  // hData_Jet80_CutA_rej->SetYTitle("Jet80_Cut Rejection");
  // hData_Jet80_CutA_rej->Draw();
  // // hData_Jet80_CutB_rej->Rebin(20);
  // // hData_Jet80_CutB_rej->Scale(1./20);
  // // hData_Jet80_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet80_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet80_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet80_CutB_rej->Draw("same");

  // cCutRejection_Jet80->SaveAs("Pp_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection.pdf","RECREATE");

  TCanvas * cJet40_CutEfficiency_Jet40 = new TCanvas("cJet40_CutEfficiency_Jet40","",1000,800);
  cJet40_CutEfficiency_Jet40->Divide(2,1);
  cJet40_CutEfficiency_Jet40->cd(1);

  hData_Jet40_CutA_eff->Rebin(5);
  hData_Jet40_CutA_eff->Scale(1./5);
  hData_Jet40_CutA_eff->SetMarkerColor(kRed);
  hData_Jet40_CutA_eff->SetMarkerStyle(24);
  hData_Jet40_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet40_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet40_CutA_eff->SetXTitle(Form("ak%dPF p_{T}",radius));
  hData_Jet40_CutA_eff->SetTitle("Data");
  hData_Jet40_CutA_eff->SetYTitle("Jet40_Cut efficiency");
  hData_Jet40_CutA_eff->Draw();
  // hData_Jet40_CutB_eff->Rebin(20);
  // hData_Jet40_CutB_eff->Scale(1./20);
  // hData_Jet40_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet40_CutB_eff->SetMarkerStyle(33);
  // hData_Jet40_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet40_CutB_eff->SetAxisRange(0,1.2,"Y");
  // hData_Jet40_CutB_eff->Draw("same");

  cJet40_CutEfficiency_Jet40->cd(2);
  hMC_Jet40_CutA_eff->Rebin(5);
  hMC_Jet40_CutA_eff->Scale(1./5);
  hMC_Jet40_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet40_CutA_eff->SetMarkerStyle(24);
  hMC_Jet40_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet40_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet40_CutA_eff->SetXTitle(Form("ak%dPF ref p_{T}",radius));
  hMC_Jet40_CutA_eff->SetTitle("MC");
  hMC_Jet40_CutA_eff->SetYTitle("Jet40_Cut efficiency");
  hMC_Jet40_CutA_eff->Draw();
  // hMC_Jet40_CutB_eff->Rebin(20);
  // hMC_Jet40_CutB_eff->Scale(1./20);
  // hMC_Jet40_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet40_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet40_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet40_CutB_eff->Draw("same");

  cJet40_CutEfficiency_Jet40->SaveAs(Form("Pp_YetkinCuts_Jet40_noJet60noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched_R0p%d_zoomed.pdf",radius),"RECREATE");

  // TCanvas * cCutRejection_Jet40 = new TCanvas("cCutRejection_Jet40","",1000,800);

  // hData_Jet40_CutA_rej->Rebin(5);
  // hData_Jet40_CutA_rej->Scale(1./5);
  // hData_Jet40_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet40_CutA_rej->SetMarkerStyle(24);
  // hData_Jet40_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet40_CutA_rej->SetAxisRange(0,1.2,"Y");
  // hData_Jet40_CutA_rej->SetXTitle("matched ak3PF p_{T}");
  // hData_Jet40_CutA_rej->SetTitle("Data");
  // hData_Jet40_CutA_rej->SetYTitle("Jet40_Cut Rejection");
  // hData_Jet40_CutA_rej->Draw();
  // // hData_Jet40_CutB_rej->Rebin(20);
  // // hData_Jet40_CutB_rej->Scale(1./20);
  // // hData_Jet40_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet40_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet40_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet40_CutB_rej->Draw("same");

  // cCutRejection_Jet40->SaveAs("Pp_YetkinCuts_Jet40_noJet60noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");

  TCanvas * cJet60_CutEfficiency_Jet60 = new TCanvas("cJet60_CutEfficiency_Jet60","",1000,800);
  cJet60_CutEfficiency_Jet60->Divide(2,1);
  cJet60_CutEfficiency_Jet60->cd(1);

  hData_Jet60_CutA_eff->Rebin(5);
  hData_Jet60_CutA_eff->Scale(1./5);
  hData_Jet60_CutA_eff->SetMarkerColor(kRed);
  hData_Jet60_CutA_eff->SetMarkerStyle(24);
  hData_Jet60_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet60_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet60_CutA_eff->SetXTitle(Form("matched ak%dPF p_{T}",radius));
  hData_Jet60_CutA_eff->SetTitle("Data");
  hData_Jet60_CutA_eff->SetYTitle("Jet60_Cut efficiency");
  hData_Jet60_CutA_eff->Draw();
  // hData_Jet60_CutB_eff->Rebin(20);
  // hData_Jet60_CutB_eff->Scale(1./20);
  // hData_Jet60_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet60_CutB_eff->SetMarkerStyle(33);
  // hData_Jet60_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet60_CutB_eff->Draw("same");

  cJet60_CutEfficiency_Jet60->cd(2);
  hMC_Jet60_CutA_eff->Rebin(5);
  hMC_Jet60_CutA_eff->Scale(1./5);
  hMC_Jet60_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet60_CutA_eff->SetMarkerStyle(24);
  hMC_Jet60_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet60_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet60_CutA_eff->SetXTitle(Form("matched ak%dPF ref p_{T}",radius));
  hMC_Jet60_CutA_eff->SetTitle("MC");
  hMC_Jet60_CutA_eff->SetYTitle("Jet60_Cut efficiency");
  hMC_Jet60_CutA_eff->Draw();
  // hMC_Jet60_CutB_eff->Rebin(20);
  // hMC_Jet60_CutB_eff->Scale(1./20);
  // hMC_Jet60_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet60_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet60_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet60_CutB_eff->Draw("same");

  cJet60_CutEfficiency_Jet60->SaveAs(Form("Pp_YetkinCuts_Jet60_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched_R0p%d_zoomed.pdf",radius),"RECREATE");

  // TCanvas * cCutRejection_Jet60 = new TCanvas("cCutRejection_Jet60","",1000,800);

  // hData_Jet60_CutA_rej->Rebin(5);
  // hData_Jet60_CutA_rej->Scale(1./5);
  // hData_Jet60_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet60_CutA_rej->SetMarkerStyle(24);
  // hData_Jet60_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet60_CutA_rej->SetAxisRange(0, 1.2,"Y");
  // hData_Jet60_CutA_rej->SetXTitle("matched ak3PF p_{T}");
  // hData_Jet60_CutA_rej->SetTitle("Data");
  // hData_Jet60_CutA_rej->SetYTitle("Jet60_Cut Rejection");
  // hData_Jet60_CutA_rej->Draw();
  // // hData_Jet60_CutB_rej->Rebin(20);
  // // hData_Jet60_CutB_rej->Scale(1./20);
  // // hData_Jet60_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet60_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet60_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet60_CutB_rej->Draw("same");

  // cCutRejection_Jet60->SaveAs("Pp_YetkinCuts_Jet60_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");


  // plot the trigger combination from this, and the total cut efficiency:
  
  TH1F * hData_Jet80 = (TH1F*)hData_Jet80_CutA->Clone("hData_Jet80");
  TH1F * hData_Jet60 = (TH1F*)hData_Jet60_CutA->Clone("hData_Jet60");
  TH1F * hData_Jet40 = (TH1F*)hData_Jet40_CutA->Clone("hData_Jet40");

  TH1F * hData_Combined = (TH1F*)hData_Jet80->Clone("hData_Combined");
  hData_Combined->Add(hData_Jet60);
  hData_Combined->Add(hData_Jet40);

  hData_Combined->Print("base");
  
  TH1F * hData_noCut_Jet80 = (TH1F*)hData_Jet80_noCut->Clone("hData_noCut_Jet80");
  TH1F * hData_noCut_Jet60 = (TH1F*)hData_Jet60_noCut->Clone("hData_noCut_Jet60");
  TH1F * hData_noCut_Jet40 = (TH1F*)hData_Jet40_noCut->Clone("hData_noCut_Jet40");

  TH1F * hData_noCut_Combined = (TH1F*)hData_noCut_Jet80->Clone("hData_noCut_Combined");
  hData_noCut_Combined->Add(hData_noCut_Jet60);
  hData_noCut_Combined->Add(hData_noCut_Jet40);
  
  hData_noCut_Combined->Print("base");
  
  TH1F * hMC_Jet80 = (TH1F*)hMC_Jet80_CutA->Clone("hMC_Jet80");
  TH1F * hMC_Jet60 = (TH1F*)hMC_Jet60_CutA->Clone("hMC_Jet60");
  TH1F * hMC_Jet40 = (TH1F*)hMC_Jet40_CutA->Clone("hMC_Jet40");

  TH1F * hMC_Combined = (TH1F*)hMC_Jet80->Clone("hMC_Combined");
  hMC_Combined->Add(hMC_Jet60);
  hMC_Combined->Add(hMC_Jet40);

  hMC_Combined->Print("base");
  
  TH1F * hMC_noCut_Jet80 = (TH1F*)hMC_Jet80_noCut->Clone("hMC_noCut_Jet80");
  TH1F * hMC_noCut_Jet60 = (TH1F*)hMC_Jet60_noCut->Clone("hMC_noCut_Jet60");
  TH1F * hMC_noCut_Jet40 = (TH1F*)hMC_Jet40_noCut->Clone("hMC_noCut_Jet40");

  TH1F * hMC_noCut_Combined = (TH1F*)hMC_noCut_Jet80->Clone("hMC_noCut_Combined");
  hMC_noCut_Combined->Add(hMC_noCut_Jet60);
  hMC_noCut_Combined->Add(hMC_noCut_Jet40);

  hMC_noCut_Combined->Print("base");
  
  TH1F * hData_Combined_Efficiency = (TH1F*)hData_Combined->Clone("hData_Combined_Efficiency");
  hData_Combined_Efficiency->Divide(hData_noCut_Combined);
  hData_Combined_Efficiency->Print("base");
  
  TH1F * hMC_Combined_Efficiency = (TH1F*)hMC_Combined->Clone("hMC_Combined_Efficiency");
  hMC_Combined_Efficiency->Divide(hMC_noCut_Combined);
  hMC_Combined_Efficiency->Print("base");
  
  TCanvas * cCombinedEff = new TCanvas("cCombinedEff","",800,600);
  hData_Combined_Efficiency->SetXTitle("Jet p_{T}");
  hData_Combined_Efficiency->SetYTitle("Combined Jet ID cut efficiency");
  hData_Combined_Efficiency->SetMarkerStyle(20);
  hData_Combined_Efficiency->SetMarkerColor(kBlack);
  hData_Combined_Efficiency->Rebin(10);
  hData_Combined_Efficiency->Scale(1/10);
  hData_Combined_Efficiency->SetAxisRange(30,350,"X");
  hData_Combined_Efficiency->Draw();
  
  hMC_Combined_Efficiency->SetMarkerStyle(24);
  hMC_Combined_Efficiency->SetMarkerColor(kRed);
  hMC_Combined_Efficiency->Rebin(10);
  hMC_Combined_Efficiency->Scale(1/10);
  hMC_Combined_Efficiency->Draw("same");

  cCombinedEff->SaveAs(Form("Combined_trigger_pp_efficiency_YetkinCut_R0p%d.pdf",radius),"RECREATE");
  
  TCanvas * cTriggerCombination = new TCanvas("cTriggerCombination","",800,600);
  cTriggerCombination->SetLogy();

  hData_Combined->SetMarkerColor(kBlack);
  hData_Combined->SetMarkerStyle(25);
  hData_Combined->SetAxisRange(30,350,"X");
  //hData_Combined->SetXTitle("");
  //hData_Combined->Draw();

  hData_Jet80->SetMarkerColor(kRed);
  hData_Jet80->SetMarkerStyle(20);
  hData_Jet80->Draw("same");
  
  hData_Jet60->SetMarkerColor(kBlue);
  hData_Jet60->SetMarkerStyle(20);
  hData_Jet60->Draw("same");

  hData_Jet40->SetMarkerColor(kGreen);
  hData_Jet40->SetMarkerStyle(20);
  hData_Jet40->Draw("same");

  //drawText()

  cTriggerCombination->SaveAs(Form("TriggerCombination_pp_YetkinCuts_R0p%d.pdf",radius),"RECREATE");


}
