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


#include "Headers/plot.h"



static const int nbins_cent = 7;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36,40};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 

int findBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins               
  if(bin<=10)ibin=0; //! 0-5%
  else if(bin>10  && bin<=20 )ibin=1; //! 5-10%
  else if(bin>20  && bin<=60 )ibin=2;  //! 10-30%
  else if(bin>60  && bin<=100)ibin=3;  //! 30-50%
  else if(bin>100 && bin<=140)ibin=4;  //! 50-70%
  else if(bin>140 && bin<=180)ibin=5;  //! 70-90%
  else if(bin>180 && bin<=200)ibin=6;  //! 90-100%
  return ibin;
}

static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};


using namespace std;


void RAA_plot_yetkinCutEfficiency(Int_t radius = 3){

  TH1::SetDefaultSumw2();
  
  // the cut is a 3 step cut based on the different value of the calopt/pfpt - copy the following lines into your loop (with the corresponding branch address set)
  // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) hGood->Fill();
  // if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) ) hGood->Fill();
  // if(calopt/pfpt > 0.85 & eMax/Sumcand > 0.9) hGood->Fill();
  
  
  char * etaWidth = (char*)"n20_eta_p20";
  TFile * fData, * fMC; 

  if(radius == 2) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu2_20150331.root");
  if(radius == 3) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150331.root");
  //if(radius == 3) fData = TFile::Open("/export/d00/scratch/pawan/condorfiles/pbpb/ntuples/Merged_JetRAA_akPu3_PbPb_Data.root");
  if(radius == 4) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu4_20150331.root");

  if(radius == 2) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu2_20150331.root");
  if(radius == 3) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150331.root");
  if(radius == 4) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu4_20150331.root");

  TTree * Data_matched = (TTree*)fData->Get("matchedJets");
  TTree * Data_unmatched = (TTree*)fData->Get("unmatchedPFJets");

  TTree * MC_matched = (TTree*)fMC->Get("matchedJets");
  TTree * MC_unmatched = (TTree*)fMC->Get("unmatchedPFJets");


  TH1F * hMC_Jet55_noCut = new TH1F("hMC_Jet55_noCut","data from matched jets without any jet ID cut",nbins_pt,boundaries_pt);
  TH1F * hMC_Jet55_CutA = new TH1F("hMC_Jet55_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);
  TH1F * hMC_Jet55_CutA_rej = new TH1F("hMC_Jet55_CutA_rej","data from matched jets rejected by Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);

  TH1F * hMC_Jet65_noCut = new TH1F("hMC_Jet65_noCut","data from matched jets without any jet ID cut",nbins_pt,boundaries_pt);
  TH1F * hMC_Jet65_CutA = new TH1F("hMC_Jet65_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);
  TH1F * hMC_Jet65_CutA_rej = new TH1F("hMC_Jet65_CutA_rej","data from matched jets rejected with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);

  TH1F * hMC_Jet80_noCut = new TH1F("hMC_Jet80_noCut","data from matched jets without any jet ID cut",nbins_pt,boundaries_pt);
  TH1F * hMC_Jet80_CutA = new TH1F("hMC_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);
  TH1F * hMC_Jet80_CutA_rej = new TH1F("hMC_Jet80_CutA_rej","data from matched jets rejected with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);

  TH1F * hData_Jet55_noCut = new TH1F("hData_Jet55_noCut","data from matched jets without any jet ID cut",nbins_pt,boundaries_pt);
  TH1F * hData_Jet55_CutA = new TH1F("hData_Jet55_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet55_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet55_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet55_CutB->Fill()} 
  TH1F * hData_Jet55_CutA_rej = new TH1F("hData_Jet55_CutA_rej","",nbins_pt,boundaries_pt);

  TH1F * hData_Jet65_noCut = new TH1F("hData_Jet65_noCut","data from matched jets without any jet ID cut",nbins_pt,boundaries_pt);
  TH1F * hData_Jet65_CutA = new TH1F("hData_Jet65_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet65_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet65_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet65_CutB->Fill()} 
  TH1F * hData_Jet65_CutA_rej = new TH1F("hData_Jet65_CutA_rej","",nbins_pt,boundaries_pt);
  
  TH1F * hData_Jet80_noCut = new TH1F("hData_Jet80_noCut","data from matched jets without any jet ID cut",nbins_pt,boundaries_pt);
  TH1F * hData_Jet80_CutA = new TH1F("hData_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",nbins_pt,boundaries_pt);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet80_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet80_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet80_CutB->Fill()} 
  TH1F * hData_Jet80_CutA_rej = new TH1F("hData_Jet80_CutA_rej","",nbins_pt,boundaries_pt);

  TH1F * hData_unmatched_Jet80_noCut = new TH1F("hData_unmatched_Jet80_noCut","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet80_CutA = new TH1F("hData_unmatched_Jet80_CutA","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet80_CutA_rej = new TH1F("hData_unmatched_Jet80_CutA_rej","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet65_noCut = new TH1F("hData_unmatched_Jet65_noCut","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet65_CutA = new TH1F("hData_unmatched_Jet65_CutA","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet65_CutA_rej = new TH1F("hData_unmatched_Jet65_CutA_rej","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet55_noCut = new TH1F("hData_unmatched_Jet55_noCut","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet55_CutA = new TH1F("hData_unmatched_Jet55_CutA","",nbins_pt,boundaries_pt);
  TH1F * hData_unmatched_Jet55_CutA_rej = new TH1F("hData_unmatched_Jet55_CutA_rej","",nbins_pt,boundaries_pt);

  TH1F * hMC_unmatched_Jet80_noCut = new TH1F("hMC_unmatched_Jet80_noCut","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet80_CutA = new TH1F("hMC_unmatched_Jet80_CutA","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet80_CutA_rej = new TH1F("hMC_unmatched_Jet80_CutA_rej","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet65_noCut = new TH1F("hMC_unmatched_Jet65_noCut","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet65_CutA = new TH1F("hMC_unmatched_Jet65_CutA","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet65_CutA_rej = new TH1F("hMC_unmatched_Jet65_CutA_rej","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet55_noCut = new TH1F("hMC_unmatched_Jet55_noCut","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet55_CutA = new TH1F("hMC_unmatched_Jet55_CutA","",nbins_pt,boundaries_pt);
  TH1F * hMC_unmatched_Jet55_CutA_rej = new TH1F("hMC_unmatched_Jet55_CutA_rej","",nbins_pt,boundaries_pt);

  //get the spectra with the specific trigger object from the different files. 
  TH1F *hpbpb_Jet80_gen[nbins_cent],*hpbpb_Jet80_reco[nbins_cent];
  TH1F *hpbpb_Jet65_gen[nbins_cent],*hpbpb_Jet65_reco[nbins_cent];
  TH1F *hpbpb_Jet55_gen[nbins_cent],*hpbpb_Jet55_reco[nbins_cent];
  TH1F *hpbpb_JetComb_gen[nbins_cent],*hpbpb_JetComb_reco[nbins_cent];

  TH1F * hpbpb_MC_noCut[nbins_cent];
  TH1F * hpbpb_MC_Comb_noCut[nbins_cent];
  TH1F * hpbpb_MC_Jet80_noCut[nbins_cent];
  TH1F * hpbpb_MC_Jet65_noCut[nbins_cent];
  TH1F * hpbpb_MC_Jet55_noCut[nbins_cent];
  TH1F * hpbpb_Data_Comb_noCut[nbins_cent];
  TH1F * hpbpb_Data_Jet55_noCut[nbins_cent];
  TH1F * hpbpb_Data_Jet65_noCut[nbins_cent];
  TH1F * hpbpb_Data_Jet80_noCut[nbins_cent];

  TH1F *hpbpb_gen[nbins_cent],*hpbpb_reco[nbins_cent];
  TH2F *hpbpb_matrix[nbins_cent];
  TH2F *hpbpb_matrix_HLT[nbins_cent];
  TH2F *hpbpb_mcclosure_matrix[nbins_cent];
  TH2F *hpbpb_mcclosure_matrix_HLT[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_data[nbins_cent];
  TH1F *hpbpb_mcclosure_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet80_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet65_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_data[nbins_cent];
  TH1F *hpbpb_mcclosure_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet80_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet65_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_gen[nbins_cent];

  TH1F *hpbpb_TrgObj80[nbins_cent];
  TH1F *hpbpb_TrgObj65[nbins_cent];
  TH1F *hpbpb_TrgObj55[nbins_cent];
  TH1F *hpbpb_TrgObjComb[nbins_cent];

  TH1F *hpbpb_JEC_TrgObj80[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj65[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj55[nbins_cent];
  TH1F *hpbpb_JEC_TrgObjComb[nbins_cent];

  TH1F *hpbpb_Smear_TrgObj80[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj65[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj55[nbins_cent];
  TH1F *hpbpb_Smear_TrgObjComb[nbins_cent];
  
  for(int i = 0;i<nbins_cent;++i){
    //cout<<"cent bin = "<<i<<endl;

    hpbpb_MC_noCut[i] = new TH1F(Form("hpbpb_MC_noCut_cent%d",i),"",nbins_pt,boundaries_pt);
    
    hpbpb_MC_Jet55_noCut[i] = new TH1F(Form("hpbpb_MC_Jet55_noCut_cent%d",i),"",nbins_pt,boundaries_pt);
    hpbpb_MC_Jet65_noCut[i] = new TH1F(Form("hpbpb_MC_Jet65_noCut_cent%d",i),"",nbins_pt,boundaries_pt);
    hpbpb_MC_Jet80_noCut[i] = new TH1F(Form("hpbpb_MC_Jet80_noCut_cent%d",i),"",nbins_pt,boundaries_pt);
    hpbpb_MC_Comb_noCut[i] = new TH1F(Form("hpbpb_MC_Comb_noCut_cent%d",i),"",nbins_pt,boundaries_pt);

    hpbpb_Data_Jet55_noCut[i] = new TH1F(Form("hpbpb_Data_Jet55_noCut_cent%d",i),"",nbins_pt,boundaries_pt);
    hpbpb_Data_Jet65_noCut[i] = new TH1F(Form("hpbpb_Data_Jet65_noCut_cent%d",i),"",nbins_pt,boundaries_pt);
    hpbpb_Data_Jet80_noCut[i] = new TH1F(Form("hpbpb_Data_Jet80_noCut_cent%d",i),"",nbins_pt,boundaries_pt);
    hpbpb_Data_Comb_noCut[i] = new TH1F(Form("hpbpb_Data_Comb_noCut_cent%d",i),"",nbins_pt,boundaries_pt);

    hpbpb_TrgObj80[i] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_TrgObj65[i] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_TrgObj55[i] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_TrgObjComb[i] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);

    hpbpb_JEC_TrgObj80[i] = new TH1F(Form("hpbpb_JEC_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_JEC_TrgObj65[i] = new TH1F(Form("hpbpb_JEC_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_JEC_TrgObj55[i] = new TH1F(Form("hpbpb_JEC_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_JEC_TrgObjComb[i] = new TH1F(Form("hpbpb_JEC_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);

    hpbpb_Smear_TrgObj80[i] = new TH1F(Form("hpbpb_Smear_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Smear_TrgObj65[i] = new TH1F(Form("hpbpb_Smear_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Smear_TrgObj55[i] = new TH1F(Form("hpbpb_Smear_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Smear_TrgObjComb[i] = new TH1F(Form("hpbpb_Smear_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);

    hpbpb_gen[i] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    //cout<<"A"<<endl;
    hpbpb_reco[i] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("Reco jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    //cout<<"B"<<endl;
    hpbpb_matrix[i] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    hpbpb_matrix_HLT[i] = new TH2F(Form("hpbpb_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    hpbpb_mcclosure_matrix[i] = new TH2F(Form("hpbpb_mcclosure_matrix_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    hpbpb_mcclosure_matrix_HLT[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    //cout<<"C"<<endl;
    hpbpb_mcclosure_data[i] = new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_JetComb_data[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_Jet80_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 80  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_Jet65_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 65  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_Jet55_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);

    hpbpb_mcclosure_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_JetComb_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_Jet80_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet80_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_Jet65_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet65_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 65 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_mcclosure_Jet55_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);

    hpbpb_JetComb_gen[i] = new TH1F(Form("hpbpb_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_JetComb_reco[i] = new TH1F(Form("hpbpb_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Jet80_gen[i] = new TH1F(Form("hpbpb_Jet80_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Jet80_reco[i] = new TH1F(Form("hpbpb_Jet80_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Jet65_gen[i] = new TH1F(Form("hpbpb_Jet65_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Jet65_reco[i] = new TH1F(Form("hpbpb_Jet65_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Jet55_gen[i] = new TH1F(Form("hpbpb_Jet55_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_Jet55_reco[i] = new TH1F(Form("hpbpb_Jet55_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);

	
  }// centrality bin loop
  


  // Define all the histograms necessary for the analysis: 
  
  // 1 - Data, 2 - MC
  Float_t pfpt_1, pfpt_2;
  Float_t pfrefpt_2;
  Float_t calopt_1, calopt_2;
  Float_t eMax_1, eMax_2;
  Float_t chMax_1, chMax_2;
  Float_t chSum_1, chSum_2;
  Float_t phSum_1, phSum_2;
  Float_t neSum_1, neSum_2;
  Float_t muSum_1, muSum_2;
  Int_t jet55_1, jet65_1, jet80_1;
  Int_t jet55_p_1, jet65_p_1, jet80_p_1;
  Int_t jet55_2, jet65_2, jet80_2;
  Int_t jet55_p_2;
  Float_t weight;
  Int_t subid_2;
  Int_t hiBin_1, hiBin_2;

  Data_matched->SetBranchAddress("calopt",&calopt_1);
  Data_matched->SetBranchAddress("pfpt",&pfpt_1);
  Data_matched->SetBranchAddress("eMax",&eMax_1);
  Data_matched->SetBranchAddress("chMax",&chMax_1);
  Data_matched->SetBranchAddress("chSum",&chSum_1);
  Data_matched->SetBranchAddress("phSum",&phSum_1);
  Data_matched->SetBranchAddress("neSum",&neSum_1);
  Data_matched->SetBranchAddress("muSum",&muSum_1);
  Data_matched->SetBranchAddress("hiBin",&hiBin_1);
  Data_matched->SetBranchAddress("jet55",&jet55_1);
  Data_matched->SetBranchAddress("jet65",&jet65_1);
  Data_matched->SetBranchAddress("jet80",&jet80_1);
  Data_matched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_matched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_matched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  
  Data_unmatched->SetBranchAddress("pfpt",&pfpt_1);
  Data_unmatched->SetBranchAddress("eMax",&eMax_1);
  Data_unmatched->SetBranchAddress("chMax",&chMax_1);
  Data_unmatched->SetBranchAddress("chSum",&chSum_1);
  Data_unmatched->SetBranchAddress("phSum",&phSum_1);
  Data_unmatched->SetBranchAddress("neSum",&neSum_1);
  Data_unmatched->SetBranchAddress("muSum",&muSum_1);
  Data_unmatched->SetBranchAddress("hiBin",&hiBin_1);
  Data_unmatched->SetBranchAddress("jet55",&jet55_1);
  Data_unmatched->SetBranchAddress("jet65",&jet65_1);
  Data_unmatched->SetBranchAddress("jet80",&jet80_1);
  Data_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_unmatched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_unmatched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  
  MC_matched->SetBranchAddress("calopt",&calopt_2);
  MC_matched->SetBranchAddress("pfpt",&pfpt_2);
  MC_matched->SetBranchAddress("eMax",&eMax_2);
  MC_matched->SetBranchAddress("chMax",&chMax_2);
  MC_matched->SetBranchAddress("chSum",&chSum_2);
  MC_matched->SetBranchAddress("phSum",&phSum_2);
  MC_matched->SetBranchAddress("neSum",&neSum_2);
  MC_matched->SetBranchAddress("muSum",&muSum_2);
  MC_matched->SetBranchAddress("hiBin",&hiBin_2);
  MC_matched->SetBranchAddress("pfrefpt",&pfrefpt_2);
  MC_matched->SetBranchAddress("jet55",&jet55_2);
  MC_matched->SetBranchAddress("jet65",&jet65_2);
  MC_matched->SetBranchAddress("jet80",&jet80_2);
  MC_matched->SetBranchAddress("weight", &weight);
  MC_matched->SetBranchAddress("subid", &subid_2);
  MC_matched->SetBranchAddress("jet55_prescl",&jet55_p_2);
  
  MC_unmatched->SetBranchAddress("pfpt",&pfpt_2);
  MC_unmatched->SetBranchAddress("eMax",&eMax_2);
  MC_unmatched->SetBranchAddress("chMax",&chMax_2);
  MC_unmatched->SetBranchAddress("chSum",&chSum_2);
  MC_unmatched->SetBranchAddress("phSum",&phSum_2);
  MC_unmatched->SetBranchAddress("neSum",&neSum_2);
  MC_unmatched->SetBranchAddress("muSum",&muSum_2);
  MC_unmatched->SetBranchAddress("hiBin",&hiBin_2);
  MC_unmatched->SetBranchAddress("pfrefpt",&pfrefpt_2);
  MC_unmatched->SetBranchAddress("jet55",&jet55_2);
  MC_unmatched->SetBranchAddress("jet65",&jet65_2);
  MC_unmatched->SetBranchAddress("jet80",&jet80_2);
  MC_unmatched->SetBranchAddress("weight", & weight);
  MC_unmatched->SetBranchAddress("subid", &subid_2);
  MC_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_2);

  // data loop
  long entries = Data_matched->GetEntries();
  //entries = 1000;
  Float_t Jet55_prescl = 2.0475;

  // get the random value for smear systematics, TRandom rnd, value per jet = rnd.Gaus(0,1);
  TRandom rnd; 
  
  cout<<"matched Data ntuple "<<endl;
  
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    Data_matched->GetEntry(nentry);

    Int_t cBin = findBin(hiBin_1);
    if(cBin == -1 || cBin >= nbins_cent) continue;
    
    Float_t Sumcand = chSum_1 + phSum_1 + neSum_1 + muSum_1;

    if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {
      
      hData_Jet55_noCut->Fill(pfpt_1, Jet55_prescl);
      hpbpb_Data_Jet55_noCut[cBin]->Fill(pfpt_1, Jet55_prescl);

      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) {
	hData_Jet55_CutA->Fill(pfpt_1, Jet55_prescl);
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl);
      }
      if(calopt_1/pfpt_1 > 0.85){
	hData_Jet55_CutA->Fill(pfpt_1, Jet55_prescl);
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand < 0.05) {
	hData_Jet55_CutA->Fill(pfpt_1, Jet55_prescl);
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet55_CutA_rej->Fill(pfpt_1, Jet55_prescl);
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) hData_Jet55_CutA_rej->Fill(pfpt_1, Jet55_prescl);
      
    }
    
    if(jet65_1 == 1 && jet80_1 == 0 ) {
      
      hData_Jet65_noCut->Fill(pfpt_1);
      hpbpb_Data_Jet65_noCut[cBin]->Fill(pfpt_1);
	    
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)){
	hData_Jet65_CutA->Fill(pfpt_1);
	hpbpb_TrgObj65[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
      }
      if(calopt_1/pfpt_1 > 0.85) {
	hData_Jet65_CutA->Fill(pfpt_1);
	hpbpb_TrgObj65[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand < 0.05) {
	hData_Jet65_CutA->Fill(pfpt_1);
	hpbpb_TrgObj65[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet65_CutA_rej->Fill(pfpt_1);
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) hData_Jet65_CutA_rej->Fill(pfpt_1);
      
    }

    if(jet80_1 == 1) {
    
      hData_Jet80_noCut->Fill(pfpt_1);
      hpbpb_Data_Jet80_noCut[cBin]->Fill(pfpt_1);   
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) {
	hData_Jet80_CutA->Fill(pfpt_1);
	hpbpb_TrgObj80[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
      }
      if(calopt_1/pfpt_1 > 0.85){
	hData_Jet80_CutA->Fill(pfpt_1);
	hpbpb_TrgObj80[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand < 0.05){
	hData_Jet80_CutA->Fill(pfpt_1);
	hpbpb_TrgObj80[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet80_CutA_rej->Fill(pfpt_1);
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) hData_Jet80_CutA_rej->Fill(pfpt_1);
      
    }
    
  }// data ntuple loop

  // data unmatched loop:
  entries = Data_unmatched->GetEntries();
  //entries = 1000;
  cout<<"Unmatched Data ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    Data_unmatched->GetEntry(nentry);
    Int_t cBin = findBin(hiBin_1);
    if(cBin == -1 || cBin >= nbins_cent) continue;
    
    Float_t Sumcand = chSum_1 + phSum_1 + neSum_1 + muSum_1;

    if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {
    
      hData_unmatched_Jet55_noCut->Fill(pfpt_1, Jet55_prescl);
      hpbpb_Data_Jet55_noCut[cBin]->Fill(pfpt_1, Jet55_prescl);

      if(eMax_1/Sumcand < 0.05 ){
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1);
	hData_unmatched_Jet55_CutA->Fill(pfpt_1, Jet55_prescl);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl);

      }
      else hData_unmatched_Jet55_CutA_rej->Fill(pfpt_1, Jet55_prescl);
      
    }

    if(jet65_1 == 1 && jet80_1 == 0 ) {

      hData_unmatched_Jet65_noCut->Fill(pfpt_1);
      hpbpb_Data_Jet65_noCut[cBin]->Fill(pfpt_1);

      if(eMax_1/Sumcand < 0.05  ){hpbpb_TrgObj65[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
	hData_unmatched_Jet65_CutA->Fill(pfpt_1);}
      else hData_unmatched_Jet65_CutA_rej->Fill(pfpt_1);
      
    }

    if(jet80_1 == 1) {
    
      hData_unmatched_Jet80_noCut->Fill(pfpt_1);
      hpbpb_Data_Jet80_noCut[cBin]->Fill(pfpt_1);

      if(eMax_1/Sumcand < 0.05  ){hpbpb_TrgObj80[cBin]->Fill(pfpt_1);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)));
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1));
	hData_unmatched_Jet80_CutA->Fill(pfpt_1);}
      else hData_unmatched_Jet80_CutA_rej->Fill(pfpt_1);
      
    }
    
  }// data ntuple loop

  entries = MC_matched->GetEntries();
  //entries = 1000;
  // MC loop
  cout<<" looking at matched MC ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    MC_matched->GetEntry(nentry);
    
    Float_t Sumcand = chSum_2 + phSum_2 + neSum_2 + muSum_2;
    if(subid_2 != 0) continue;

    Int_t cBin = findBin(hiBin_2);
    if(cBin == -1 || cBin >= nbins_cent) continue;    
    
    if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)){
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
    }
    if(calopt_2/pfpt_2 > 0.85) {
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
    }
    if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
    }

    hpbpb_MC_noCut[cBin]->Fill(pfrefpt_2, weight);
    
    if(jet55_2 == 1 && jet65_2==0 && jet80_2 == 0){
      
      hMC_Jet55_noCut->Fill(pfrefpt_2, weight);
      hpbpb_MC_Jet55_noCut[cBin]->Fill(pfrefpt_2,weight);

      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)){
	hMC_Jet55_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

      }
      if(calopt_2/pfpt_2 > 0.85) {
	hMC_Jet55_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

	
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
	hMC_Jet55_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

	
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet55_CutA_rej->Fill(pfrefpt_2, weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) hMC_Jet55_CutA_rej->Fill(pfrefpt_2, weight);
      
    }
    
    if(jet65_2 == 1 && jet80_2 == 0){

      hMC_Jet65_noCut->Fill(pfrefpt_2, weight);
      hpbpb_MC_Jet65_noCut[cBin]->Fill(pfrefpt_2,weight);

      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) {
	hMC_Jet65_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
      }
      if(calopt_2/pfpt_2 > 0.85) {
	hMC_Jet65_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
	hMC_Jet65_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

      }

      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet65_CutA_rej->Fill(pfrefpt_2, weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) hMC_Jet65_CutA_rej->Fill(pfrefpt_2, weight);

    }

    
    if(jet80_2 == 1){

      hMC_Jet80_noCut->Fill(pfrefpt_2, weight);
      hpbpb_MC_Jet80_noCut[cBin]->Fill(pfrefpt_2,weight);
 
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) {
	hMC_Jet80_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

      }
      if(calopt_2/pfpt_2 > 0.85) {
	hMC_Jet80_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
	hMC_Jet80_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

      }
      
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet80_CutA_rej->Fill(pfrefpt_2, weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) hMC_Jet80_CutA_rej->Fill(pfrefpt_2, weight);

    }
    
    
  }// mc ntuple loop


  entries = MC_unmatched->GetEntries();
  //entries = 1000;
  // MC loop
  cout<<" looking at unmatched MC ntuple"<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    MC_unmatched->GetEntry(nentry);

    if(subid_2 != 0) continue;

    Int_t cBin = findBin(hiBin_2);
    if(cBin == -1 || cBin >= nbins_cent) continue;
    
    Float_t Sumcand = chSum_2 + phSum_2 + neSum_2 + muSum_2;

    if(eMax_2/Sumcand < 0.05  ){
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

    }

    hpbpb_MC_noCut[cBin]->Fill(pfrefpt_2, weight);

    if(jet55_2 == 1 && jet65_2==0 && jet80_2 == 0){
      hpbpb_MC_Jet55_noCut[cBin]->Fill(pfrefpt_2,weight);
      hMC_unmatched_Jet55_noCut->Fill(pfrefpt_2, weight);
      if(eMax_2/Sumcand < 0.05  ){
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

      
	hMC_unmatched_Jet55_CutA->Fill(pfrefpt_2, weight);}
      else hMC_unmatched_Jet55_CutA_rej->Fill(pfrefpt_2, weight);
      
    }

    
    if(jet65_2 == 1 && jet80_2 == 0){
      hpbpb_MC_Jet65_noCut[cBin]->Fill(pfrefpt_2,weight);
      hMC_unmatched_Jet65_noCut->Fill(pfrefpt_2, weight);
      if(eMax_2/Sumcand < 0.05  ){
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

	hMC_unmatched_Jet65_CutA->Fill(pfrefpt_2, weight);}
      else hMC_unmatched_Jet65_CutA_rej->Fill(pfrefpt_2);
      
    }

    
    if(jet80_2 == 1){
      hpbpb_MC_Jet80_noCut[cBin]->Fill(pfrefpt_2,weight);
      hMC_unmatched_Jet80_noCut->Fill(pfrefpt_2, weight);
      if(eMax_2/Sumcand < 0.05  ){
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);

	hMC_unmatched_Jet80_CutA->Fill(pfrefpt_2, weight);}
      
      else hMC_unmatched_Jet80_CutA_rej->Fill(pfrefpt_2, weight);
      
    }
    
    
  }// mc unmatched  ntuple loop

  TFile fout(Form("../../Output/PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_rebinned_eMaxSumcand_A_R0p%d.root",radius),"RECREATE");
  fout.cd();
  
  for(int i = 0;i<nbins_cent;++i){

    hpbpb_MC_Comb_noCut[i]->Add(hpbpb_MC_Jet80_noCut[i]);
    hpbpb_MC_Comb_noCut[i]->Add(hpbpb_MC_Jet65_noCut[i]);
    hpbpb_MC_Comb_noCut[i]->Add(hpbpb_MC_Jet55_noCut[i]);

    divideBinWidth(hpbpb_MC_Comb_noCut[i]);
    divideBinWidth(hpbpb_MC_Jet80_noCut[i]);
    divideBinWidth(hpbpb_MC_Jet65_noCut[i]);
    divideBinWidth(hpbpb_MC_Jet55_noCut[i]);

    hpbpb_Data_Comb_noCut[i]->Add(hpbpb_Data_Jet80_noCut[i]);
    hpbpb_Data_Comb_noCut[i]->Add(hpbpb_Data_Jet65_noCut[i]);
    hpbpb_Data_Comb_noCut[i]->Add(hpbpb_Data_Jet55_noCut[i]);

    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj80[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj65[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj55[i]);

    divideBinWidth(hpbpb_TrgObjComb[i]);
    divideBinWidth(hpbpb_TrgObj80[i]);
    divideBinWidth(hpbpb_TrgObj65[i]);
    divideBinWidth(hpbpb_TrgObj55[i]);

    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj80[i]);
    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj65[i]);
    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj55[i]);

    divideBinWidth(hpbpb_Smear_TrgObjComb[i]);

    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj80[i]);
    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj65[i]);
    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj55[i]);

    divideBinWidth(hpbpb_JEC_TrgObjComb[i]);
    
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet55_data[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_data[i]);
    
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet55_gen[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_gen[i]);	
	
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet80_reco[i]);
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet65_reco[i]);
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet55_reco[i]);
    
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet80_gen[i]);
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet65_gen[i]);
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet55_gen[i]);

    divideBinWidth(hpbpb_JetComb_gen[i]);
    divideBinWidth(hpbpb_JetComb_reco[i]);
    
  }

  for(int i = 0; i<nbins_cent;++i){

    hpbpb_MC_noCut[i]->Write();
    
    hpbpb_MC_Comb_noCut[i]->Write();
    hpbpb_Data_Comb_noCut[i]->Write();
    

    hpbpb_TrgObjComb[i]->Write();
    hpbpb_TrgObj80[i]->Write();
    hpbpb_TrgObj65[i]->Write();
    hpbpb_TrgObj55[i]->Write();

    hpbpb_JEC_TrgObjComb[i]->Write();
    hpbpb_JEC_TrgObj80[i]->Write();
    hpbpb_JEC_TrgObj65[i]->Write();
    hpbpb_JEC_TrgObj55[i]->Write();

    hpbpb_Smear_TrgObjComb[i]->Write();
    hpbpb_Smear_TrgObj80[i]->Write();
    hpbpb_Smear_TrgObj65[i]->Write();
    hpbpb_Smear_TrgObj55[i]->Write();
    
    hpbpb_matrix_HLT[i]->Write();
    hpbpb_mcclosure_matrix_HLT[i]->Write();
    hpbpb_mcclosure_JetComb_data[i]->Write();
    hpbpb_mcclosure_Jet80_data[i]->Write();
    hpbpb_mcclosure_Jet65_data[i]->Write();
    hpbpb_mcclosure_Jet55_data[i]->Write();
    hpbpb_mcclosure_JetComb_gen[i]->Write();    
    hpbpb_mcclosure_Jet80_gen[i]->Write();
    hpbpb_mcclosure_Jet65_gen[i]->Write();
    hpbpb_mcclosure_Jet55_gen[i]->Write();

    hpbpb_JetComb_reco[i]->Write();
    hpbpb_Jet80_reco[i]->Write();
    hpbpb_Jet65_reco[i]->Write();
    hpbpb_Jet55_reco[i]->Write();
    hpbpb_JetComb_gen[i]->Write();
    hpbpb_Jet80_gen[i]->Write();
    hpbpb_Jet65_gen[i]->Write();
    hpbpb_Jet55_gen[i]->Write();
 
    hpbpb_reco[i]->Write();
    hpbpb_gen[i]->Write();
    hpbpb_matrix[i]->Write();

  }

 
  // add the unmatched histograms to the matched ones to get the final cut efficiency
  hData_Jet55_noCut->Add(hData_unmatched_Jet55_noCut);
  hData_Jet65_noCut->Add(hData_unmatched_Jet65_noCut);
  hData_Jet80_noCut->Add(hData_unmatched_Jet80_noCut);
  
  hData_Jet55_CutA->Add(hData_unmatched_Jet55_CutA);
  hData_Jet65_CutA->Add(hData_unmatched_Jet65_CutA);
  hData_Jet80_CutA->Add(hData_unmatched_Jet80_CutA);

  hData_Jet55_CutA_rej->Add(hData_unmatched_Jet55_CutA_rej);
  hData_Jet65_CutA_rej->Add(hData_unmatched_Jet65_CutA_rej);
  hData_Jet80_CutA_rej->Add(hData_unmatched_Jet80_CutA_rej);

  hMC_Jet55_noCut->Add(hMC_unmatched_Jet55_noCut);
  hMC_Jet65_noCut->Add(hMC_unmatched_Jet65_noCut);
  hMC_Jet80_noCut->Add(hMC_unmatched_Jet80_noCut);
  
  hMC_Jet55_CutA->Add(hMC_unmatched_Jet55_CutA);
  hMC_Jet65_CutA->Add(hMC_unmatched_Jet65_CutA);
  hMC_Jet80_CutA->Add(hMC_unmatched_Jet80_CutA);

  hMC_Jet55_CutA_rej->Add(hMC_unmatched_Jet55_CutA_rej);
  hMC_Jet65_CutA_rej->Add(hMC_unmatched_Jet65_CutA_rej);
  hMC_Jet80_CutA_rej->Add(hMC_unmatched_Jet80_CutA_rej);

  

  hData_Jet65_noCut->Write();
  hData_Jet65_CutA->Write();
  // hData_Jet65_CutB->Write();
  TH1F * hData_Jet65_CutA_eff = (TH1F*)hData_Jet65_CutA->Clone("hData_Jet65_CutA_eff");
  hData_Jet65_CutA_eff->Divide(hData_Jet65_noCut);
  hData_Jet65_CutA_eff->Write();
  // TH1F * hData_Jet65_CutB_eff = (TH1F*)hData_Jet65_CutB->Clone("hData_Jet65_CutB_eff");
  // hData_Jet65_CutB_eff->Divide(hData_Jet65_noCut);
  // hData_Jet65_CutB_eff->Write();
  hData_Jet65_CutA_rej->Write();
  // hData_Jet65_CutB_rej->Write();
  TH1F * hData_Jet65_CutA_rej_eff = (TH1F*)hData_Jet65_CutA_rej->Clone("hData_Jet65_CutA_rej_eff");
  hData_Jet65_CutA_rej_eff->Divide(hData_Jet65_noCut);
  hData_Jet65_CutA_rej_eff->Write();
  // hData_Jet65_CutB_rej_eff = (TH1F*)hData_Jet65_CutB_rej->Clone("hData_Jet65_CutB_rej_eff");
  // hData_Jet65_CutB_rej_eff->Divide(hData_Jet65_noCut);
  // hData_Jet65_CutB_rej_eff->Write();

  hData_Jet55_noCut->Write();
  hData_Jet55_CutA->Write();
  // hData_Jet55_CutB->Write();
  TH1F * hData_Jet55_CutA_eff = (TH1F*)hData_Jet55_CutA->Clone("hData_Jet55_CutA_eff");
  hData_Jet55_CutA_eff->Divide(hData_Jet55_noCut);
  hData_Jet55_CutA_eff->Write();
  // TH1F * hData_Jet55_CutB_eff = (TH1F*)hData_Jet55_CutB->Clone("hData_Jet55_CutB_eff");
  // hData_Jet55_CutB_eff->Divide(hData_Jet55_noCut);
  // hData_Jet55_CutB_eff->Write();
  hData_Jet55_CutA_rej->Write();
  // hData_Jet55_CutB_rej->Write();
  TH1F * hData_Jet55_CutA_rej_eff = (TH1F*)hData_Jet55_CutA_rej->Clone("hData_Jet55_CutA_rej_eff");
  hData_Jet55_CutA_rej_eff->Divide(hData_Jet55_noCut);
  hData_Jet55_CutA_rej_eff->Write();
  // hData_Jet55_CutB_rej_eff = (TH1F*)hData_Jet55_CutB_rej->Clone("hData_Jet55_CutB_rej_eff");
  // hData_Jet55_CutB_rej_eff->Divide(hData_Jet55_noCut);
  // hData_Jet55_CutB_rej_eff->Write();

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
  
  hMC_Jet65_noCut->Write();
  hMC_Jet65_CutA->Write();
  // hMC_Jet65_CutB->Write();
  TH1F * hMC_Jet65_CutA_eff = (TH1F*)hMC_Jet65_CutA->Clone("hMC_Jet65_CutA_eff");
  hMC_Jet65_CutA_eff->Divide(hMC_Jet65_noCut);
  hMC_Jet65_CutA_eff->Write();
  // TH1F * hMC_Jet65_CutB_eff = (TH1F*)hMC_Jet65_CutB->Clone("hMC_Jet65_CutB_eff");
  // hMC_Jet65_CutB_eff->Divide(hMC_Jet65_noCut);
  // hMC_Jet65_CutB_eff->Write();
  
  hMC_Jet55_noCut->Write();
  hMC_Jet55_CutA->Write();
  // hMC_Jet55_CutB->Write();
  TH1F * hMC_Jet55_CutA_eff = (TH1F*)hMC_Jet55_CutA->Clone("hMC_Jet55_CutA_eff");
  hMC_Jet55_CutA_eff->Divide(hMC_Jet55_noCut);
  hMC_Jet55_CutA_eff->Write();
  // TH1F * hMC_Jet55_CutB_eff = (TH1F*)hMC_Jet55_CutB->Clone("hMC_Jet55_CutB_eff");
  // hMC_Jet55_CutB_eff->Divide(hMC_Jet55_noCut);
  // hMC_Jet55_CutB_eff->Write();

  // save the unmatched histograms as well:
  hData_unmatched_Jet80_noCut->Write();
  hData_unmatched_Jet80_CutA->Write();
  hData_unmatched_Jet80_CutA_rej->Write();
  hData_unmatched_Jet65_noCut->Write();
  hData_unmatched_Jet65_CutA->Write();
  hData_unmatched_Jet65_CutA_rej->Write();
  hData_unmatched_Jet55_noCut->Write();
  hData_unmatched_Jet55_CutA->Write();
  hData_unmatched_Jet55_CutA_rej->Write();

  hMC_unmatched_Jet80_noCut->Write();
  hMC_unmatched_Jet80_CutA->Write();
  hMC_unmatched_Jet80_CutA_rej->Write();
  hMC_unmatched_Jet65_noCut->Write();
  hMC_unmatched_Jet65_CutA->Write();
  hMC_unmatched_Jet65_CutA_rej->Write();
  hMC_unmatched_Jet55_noCut->Write();
  hMC_unmatched_Jet55_CutA->Write();
  hMC_unmatched_Jet55_CutA_rej->Write();
  
  
  TCanvas * cJet80_CutEfficiency_Jet80 = new TCanvas("cJet80_CutEfficiency_Jet80","",1000,800);
  cJet80_CutEfficiency_Jet80->Divide(2,1);
  cJet80_CutEfficiency_Jet80->cd(1);

  hData_Jet80_CutA_eff->Rebin(5);
  hData_Jet80_CutA_eff->Scale(1./5);
  hData_Jet80_CutA_eff->SetMarkerColor(kRed);
  hData_Jet80_CutA_eff->SetMarkerStyle(24);
  hData_Jet80_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet80_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet80_CutA_eff->SetXTitle(Form("akPu%dPF p_{T}",radius));
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
  hMC_Jet80_CutA_eff->SetXTitle(Form("akPu%dPF ref p_{T}",radius));
  hMC_Jet80_CutA_eff->SetTitle("MC");
  hMC_Jet80_CutA_eff->SetYTitle("Jet80_Cut efficiency");
  hMC_Jet80_CutA_eff->Draw();
  // hMC_Jet80_CutB_eff->Rebin(20);
  // hMC_Jet80_CutB_eff->Scale(1./20);
  // hMC_Jet80_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet80_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet80_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet80_CutB_eff->Draw("same");

  cJet80_CutEfficiency_Jet80->SaveAs(Form("PbPb_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_R0p%d_zoomed.pdf",radius),"RECREATE");

  // TCanvas * cCutRejection_Jet80 = new TCanvas("cCutRejection_Jet80","",1000,800);

  // hData_Jet80_CutA_rej->Rebin(5);
  // hData_Jet80_CutA_rej->Scale(1./5);
  // hData_Jet80_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet80_CutA_rej->SetMarkerStyle(24);
  // hData_Jet80_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet80_CutA_rej->SetAxisRange(0,1.2,"Y");
  // hData_Jet80_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  // hData_Jet80_CutA_rej->SetTitle("Data");
  // hData_Jet80_CutA_rej->SetYTitle("Jet80_Cut Rejection");
  // hData_Jet80_CutA_rej->Draw();
  // // hData_Jet80_CutB_rej->Rebin(20);
  // // hData_Jet80_CutB_rej->Scale(1./20);
  // // hData_Jet80_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet80_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet80_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet80_CutB_rej->Draw("same");

  // cCutRejection_Jet80->SaveAs("PbPb_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection.pdf","RECREATE");

  TCanvas * cJet55_CutEfficiency_Jet55 = new TCanvas("cJet55_CutEfficiency_Jet55","",1000,800);
  cJet55_CutEfficiency_Jet55->Divide(2,1);
  cJet55_CutEfficiency_Jet55->cd(1);

  hData_Jet55_CutA_eff->Rebin(5);
  hData_Jet55_CutA_eff->Scale(1./5);
  hData_Jet55_CutA_eff->SetMarkerColor(kRed);
  hData_Jet55_CutA_eff->SetMarkerStyle(24);
  hData_Jet55_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet55_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet55_CutA_eff->SetXTitle(Form("akPu%dPF p_{T}",radius));
  hData_Jet55_CutA_eff->SetTitle("Data");
  hData_Jet55_CutA_eff->SetYTitle("Jet55_Cut efficiency");
  hData_Jet55_CutA_eff->Draw();
  // hData_Jet55_CutB_eff->Rebin(20);
  // hData_Jet55_CutB_eff->Scale(1./20);
  // hData_Jet55_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet55_CutB_eff->SetMarkerStyle(33);
  // hData_Jet55_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet55_CutB_eff->SetAxisRange(0,1.2,"Y");
  // hData_Jet55_CutB_eff->Draw("same");

  cJet55_CutEfficiency_Jet55->cd(2);
  hMC_Jet55_CutA_eff->Rebin(5);
  hMC_Jet55_CutA_eff->Scale(1./5);
  hMC_Jet55_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet55_CutA_eff->SetMarkerStyle(24);
  hMC_Jet55_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet55_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet55_CutA_eff->SetXTitle(Form("akPu%dPF ref p_{T}",radius));
  hMC_Jet55_CutA_eff->SetTitle("MC");
  hMC_Jet55_CutA_eff->SetYTitle("Jet55_Cut efficiency");
  hMC_Jet55_CutA_eff->Draw();
  // hMC_Jet55_CutB_eff->Rebin(20);
  // hMC_Jet55_CutB_eff->Scale(1./20);
  // hMC_Jet55_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet55_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet55_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet55_CutB_eff->Draw("same");

  cJet55_CutEfficiency_Jet55->SaveAs(Form("PbPb_YetkinCuts_Jet55_noJet65noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched_R0p%d_zoomed.pdf",radius),"RECREATE");

  // TCanvas * cCutRejection_Jet55 = new TCanvas("cCutRejection_Jet55","",1000,800);

  // hData_Jet55_CutA_rej->Rebin(5);
  // hData_Jet55_CutA_rej->Scale(1./5);
  // hData_Jet55_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet55_CutA_rej->SetMarkerStyle(24);
  // hData_Jet55_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet55_CutA_rej->SetAxisRange(0,1.2,"Y");
  // hData_Jet55_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  // hData_Jet55_CutA_rej->SetTitle("Data");
  // hData_Jet55_CutA_rej->SetYTitle("Jet55_Cut Rejection");
  // hData_Jet55_CutA_rej->Draw();
  // // hData_Jet55_CutB_rej->Rebin(20);
  // // hData_Jet55_CutB_rej->Scale(1./20);
  // // hData_Jet55_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet55_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet55_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet55_CutB_rej->Draw("same");

  // cCutRejection_Jet55->SaveAs("PbPb_YetkinCuts_Jet55_noJet65noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");

  TCanvas * cJet65_CutEfficiency_Jet65 = new TCanvas("cJet65_CutEfficiency_Jet65","",1000,800);
  cJet65_CutEfficiency_Jet65->Divide(2,1);
  cJet65_CutEfficiency_Jet65->cd(1);

  hData_Jet65_CutA_eff->Rebin(5);
  hData_Jet65_CutA_eff->Scale(1./5);
  hData_Jet65_CutA_eff->SetMarkerColor(kRed);
  hData_Jet65_CutA_eff->SetMarkerStyle(24);
  hData_Jet65_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet65_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet65_CutA_eff->SetXTitle(Form("matched akPu%dPF p_{T}",radius));
  hData_Jet65_CutA_eff->SetTitle("Data");
  hData_Jet65_CutA_eff->SetYTitle("Jet65_Cut efficiency");
  hData_Jet65_CutA_eff->Draw();
  // hData_Jet65_CutB_eff->Rebin(20);
  // hData_Jet65_CutB_eff->Scale(1./20);
  // hData_Jet65_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet65_CutB_eff->SetMarkerStyle(33);
  // hData_Jet65_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet65_CutB_eff->Draw("same");

  cJet65_CutEfficiency_Jet65->cd(2);
  hMC_Jet65_CutA_eff->Rebin(5);
  hMC_Jet65_CutA_eff->Scale(1./5);
  hMC_Jet65_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet65_CutA_eff->SetMarkerStyle(24);
  hMC_Jet65_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet65_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet65_CutA_eff->SetXTitle(Form("matched akPu%dPF ref p_{T}",radius));
  hMC_Jet65_CutA_eff->SetTitle("MC");
  hMC_Jet65_CutA_eff->SetYTitle("Jet65_Cut efficiency");
  hMC_Jet65_CutA_eff->Draw();
  // hMC_Jet65_CutB_eff->Rebin(20);
  // hMC_Jet65_CutB_eff->Scale(1./20);
  // hMC_Jet65_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet65_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet65_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet65_CutB_eff->Draw("same");

  cJet65_CutEfficiency_Jet65->SaveAs(Form("PbPb_YetkinCuts_Jet65_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched_R0p%d_zoomed.pdf",radius),"RECREATE");

  // TCanvas * cCutRejection_Jet65 = new TCanvas("cCutRejection_Jet65","",1000,800);

  // hData_Jet65_CutA_rej->Rebin(5);
  // hData_Jet65_CutA_rej->Scale(1./5);
  // hData_Jet65_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet65_CutA_rej->SetMarkerStyle(24);
  // hData_Jet65_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet65_CutA_rej->SetAxisRange(0, 1.2,"Y");
  // hData_Jet65_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  // hData_Jet65_CutA_rej->SetTitle("Data");
  // hData_Jet65_CutA_rej->SetYTitle("Jet65_Cut Rejection");
  // hData_Jet65_CutA_rej->Draw();
  // // hData_Jet65_CutB_rej->Rebin(20);
  // // hData_Jet65_CutB_rej->Scale(1./20);
  // // hData_Jet65_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet65_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet65_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet65_CutB_rej->Draw("same");

  // cCutRejection_Jet65->SaveAs("PbPb_YetkinCuts_Jet65_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");


  // plot the trigger combination from this, and the total cut efficiency:
  
  TH1F * hData_Jet80 = (TH1F*)hData_Jet80_CutA->Clone("hData_Jet80");
  TH1F * hData_Jet65 = (TH1F*)hData_Jet65_CutA->Clone("hData_Jet65");
  TH1F * hData_Jet55 = (TH1F*)hData_Jet55_CutA->Clone("hData_Jet55");

  TH1F * hData_Combined = (TH1F*)hData_Jet80->Clone("hData_Combined");
  hData_Combined->Add(hData_Jet65);
  hData_Combined->Add(hData_Jet55);

  hData_Combined->Print("base");
  
  TH1F * hData_noCut_Jet80 = (TH1F*)hData_Jet80_noCut->Clone("hData_noCut_Jet80");
  TH1F * hData_noCut_Jet65 = (TH1F*)hData_Jet65_noCut->Clone("hData_noCut_Jet65");
  TH1F * hData_noCut_Jet55 = (TH1F*)hData_Jet55_noCut->Clone("hData_noCut_Jet55");

  TH1F * hData_noCut_Combined = (TH1F*)hData_noCut_Jet80->Clone("hData_noCut_Combined");
  hData_noCut_Combined->Add(hData_noCut_Jet65);
  hData_noCut_Combined->Add(hData_noCut_Jet55);
  
  hData_noCut_Combined->Print("base");
  
  TH1F * hMC_Jet80 = (TH1F*)hMC_Jet80_CutA->Clone("hMC_Jet80");
  TH1F * hMC_Jet65 = (TH1F*)hMC_Jet65_CutA->Clone("hMC_Jet65");
  TH1F * hMC_Jet55 = (TH1F*)hMC_Jet55_CutA->Clone("hMC_Jet55");

  TH1F * hMC_Combined = (TH1F*)hMC_Jet80->Clone("hMC_Combined");
  hMC_Combined->Add(hMC_Jet65);
  hMC_Combined->Add(hMC_Jet55);

  hMC_Combined->Print("base");
  
  TH1F * hMC_noCut_Jet80 = (TH1F*)hMC_Jet80_noCut->Clone("hMC_noCut_Jet80");
  TH1F * hMC_noCut_Jet65 = (TH1F*)hMC_Jet65_noCut->Clone("hMC_noCut_Jet65");
  TH1F * hMC_noCut_Jet55 = (TH1F*)hMC_Jet55_noCut->Clone("hMC_noCut_Jet55");

  TH1F * hMC_noCut_Combined = (TH1F*)hMC_noCut_Jet80->Clone("hMC_noCut_Combined");
  hMC_noCut_Combined->Add(hMC_noCut_Jet65);
  hMC_noCut_Combined->Add(hMC_noCut_Jet55);

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

  cCombinedEff->SaveAs(Form("Combined_trigger_efficiency_YetkinCut_R0p%d.pdf",radius),"RECREATE");
  
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
  
  hData_Jet65->SetMarkerColor(kBlue);
  hData_Jet65->SetMarkerStyle(20);
  hData_Jet65->Draw("same");

  hData_Jet55->SetMarkerColor(kGreen);
  hData_Jet55->SetMarkerStyle(20);
  hData_Jet55->Draw("same");

  //drawText()

  cTriggerCombination->SaveAs(Form("TriggerCombination_YetkinCuts_R0p%d.pdf",radius),"RECREATE");


}
