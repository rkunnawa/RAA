// Raghav Kunnawalkam Elayavalli
// April 13 2015
// Rutgers
// raghav.k.e at CERN dot CH

//
// Macro to plot the systematics from the RAA analysis.
//


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <fstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TColor.h>
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
#include "../Headers/utilities.h"

// May29th after HIN-W meeting:  im going to plot the Rcp here and increase the npart plots to 2 more pT ranges. 

TLegend *getLegend(double x1, double y1, double x2, double y2)
{
  TLegend *leg = new TLegend(x1,y1,x2,y2,NULL,"BRNDC");
  leg->SetHeader("");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(1001);
  return leg;
}


// all the pt bins are declared in the utilities.h header file to make things easier. 

using namespace std;

void RAA_plot_Systematics(int radius = 2, char *algo = "Pu", char *jet_type = "PF", int unfoldingCut = 40){
  
  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  char * etaWidth = (char*) "20_eta_20";
  Float_t etaBoundary = 2.0; 

  TFile *fin_R2, *fin_R3, *fin_R4; 
  fin_R2 = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_SevilFakeMBnoJet80Cut_unfold_mcclosure_oppside_trgMC_fullfinebins_noSmear_%s_%dGeVCut_ak%s_20150603.root",2,etaWidth,30,jet_type));
  //fin_R3 = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_SevilNewMBSubtraction_unfold_mcclosure_oppside_trgMC_noSmear_%s_%dGeVCut_ak%s_20150522.root",3,etaWidth,40,jet_type));
  fin_R3 = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_SevilFakeMBnoJet80Cut_unfold_mcclosure_oppside_trgMC_fullfinebins_noSmear_%s_%dGeVCut_ak%s_20150603.root",3,etaWidth,40,jet_type));
  //fin_R4 = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_SevilFakeMBnoJet80Cut_unfold_mcclosure_oppside_trgMC_noSmear_%s_%dGeVCut_ak%s_20150528.root",4,etaWidth,50,jet_type));
  fin_R4 = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_SevilFakeMBnoJet80Cut_unfold_mcclosure_oppside_trgMC_fullfinebins_noSmear_%s_%dGeVCut_ak%s_20150603.root",4,etaWidth,50,jet_type));
  //fin_R4 = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_without80FakeRemoval_unfold_mcclosure_oppside_trgMC_noSmear_%s_%dGeVCut_ak%s_20150521.root",4,etaWidth,50,jet_type));
  
  
  // // get the unfolded error correction files and histograms
  TFile * fError_R2, * fError_R3, * fError_R4;
  // fError_R2 = TFile::Open(Form("Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_%s_unfoldingCut_%d_MinBiasFakeCut_NoJet80_data_driven_correction_akPu%s.root",2,2, etaWidth, 40, jet_type));
  // fError_R3 = TFile::Open(Form("Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_%s_unfoldingCut_%d_MinBiasFakeCut_NoJet80_data_driven_correction_akPu%s.root",3,3, etaWidth, 40, jet_type));
  // fError_R4 = TFile::Open(Form("Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_%s_unfoldingCut_%d_MinBiasFakeCut_NoJet80_data_driven_correction_akPu%s.root",4,4, etaWidth, 40, jet_type));

  fError_R2 = TFile::Open(Form("Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_fullfinebin_%s_unfoldingCut_%d_SevilFakeMBnoJet80Cut_data_driven_correction_akPu%s.root",2,2, etaWidth, 40, jet_type));
  //fError_R3 = TFile::Open(Form("Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_bin_%s_unfoldingCut_%d_SevilfakeCut_data_driven_correction_akPu%s.root",3,3, etaWidth, 40, jet_type));
  fError_R3 = TFile::Open(Form("Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_fullfinebin_%s_unfoldingCut_%d_SevilFakeMBnoJet80Cut_data_driven_correction_akPu%s.root",3,3, etaWidth, 40, jet_type));
  fError_R4 = TFile::Open(Form("Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_fullfinebin_%s_unfoldingCut_%d_SevilFakeMBnoJet80Cut_data_driven_correction_akPu%s.root",4,4, etaWidth, 40, jet_type));
  
  TH1F * hPbPb_R2_ErrorFix[nbins_cent], * hPbPb_R3_ErrorFix[nbins_cent], * hPbPb_R4_ErrorFix[nbins_cent];
  TH1F * hPbPb_R2_measured[nbins_cent], * hPbPb_R3_measured[nbins_cent], * hPbPb_R4_measured[nbins_cent];
  
  // get the histograms.
  TH1F * uPbPb_R2_Bayes[nbins_cent], * uPP_R2_Bayes, * uPbPb_R3_Bayes[nbins_cent], * uPP_R3_Bayes, * uPbPb_R4_Bayes[nbins_cent], * uPP_R4_Bayes;
  TH1F * mPbPb_R2[nbins_cent], * mPP_R2, * mPbPb_R3[nbins_cent], * mPP_R3, * mPbPb_R4[nbins_cent], * mPP_R4;
  
  TH1F * RAA_R2_Bayes[nbins_cent], * RAA_R3_Bayes[nbins_cent], * RAA_R4_Bayes[nbins_cent];
  TH1F * RAA_R2_BinByBin[nbins_cent], * RAA_R3_BinByBin[nbins_cent], * RAA_R4_BinByBin[nbins_cent];
  TH1F * RAA_R2_Meas[nbins_cent], * RAA_R3_Meas[nbins_cent], * RAA_R4_Meas[nbins_cent];

  uPP_R2_Bayes = (TH1F*)fin_R2->Get("PP_bayesian_unfolded_spectra");
  uPP_R2_Bayes->Print("base");
  uPP_R3_Bayes = (TH1F*)fin_R3->Get("PP_bayesian_unfolded_spectra");
  uPP_R3_Bayes->Print("base");
  uPP_R4_Bayes = (TH1F*)fin_R4->Get("PP_bayesian_unfolded_spectra");
  uPP_R4_Bayes->Print("base");

  mPP_R2 = (TH1F*)fin_R2->Get("PP_Gen_spectra_refpt");
  mPP_R2->Print("base");
  mPP_R3 = (TH1F*)fin_R3->Get("PP_Gen_spectra_refpt");
  mPP_R3->Print("base");
  mPP_R4 = (TH1F*)fin_R4->Get("PP_Gen_spectra_refpt");
  mPP_R4->Print("base");

  TH1F * Rcp_R2_Bayes[2], * Rcp_R3_Bayes[2], * Rcp_R4_Bayes[2];
  TH1F * Rcp3_vs_Rcp2[2], * Rcp4_vs_Rcp2[2]; 

  TH1F * RAA_R2_Bayes_fineBin[nbins_cent],* RAA_R3_Bayes_fineBin[nbins_cent],* RAA_R4_Bayes_fineBin[nbins_cent];
  TH1F * RAA_R2_Bayes_atlasBin[nbins_cent],* RAA_R3_Bayes_atlasBin[nbins_cent],* RAA_R4_Bayes_atlasBin[nbins_cent];
  TH1F * RAA_R2_Bayes_atlasRcpBin[nbins_cent],* RAA_R3_Bayes_atlasRcpBin[nbins_cent],* RAA_R4_Bayes_atlasRcpBin[nbins_cent];

  
  for(int i = 0; i<nbins_cent; ++i){
    cout<<i<<endl;
    hPbPb_R2_ErrorFix[i] = (TH1F*)fError_R2->Get(Form("PbPb_BayesianUnfolded_cent%d",i));
    hPbPb_R2_ErrorFix[i]->Print("base");
    hPbPb_R3_ErrorFix[i] = (TH1F*)fError_R3->Get(Form("PbPb_BayesianUnfolded_cent%d",i));
    hPbPb_R3_ErrorFix[i]->Print("base");
    hPbPb_R4_ErrorFix[i] = (TH1F*)fError_R4->Get(Form("PbPb_BayesianUnfolded_cent%d",i));
    hPbPb_R4_ErrorFix[i]->Print("base");

    hPbPb_R2_measured[i] = (TH1F*)fError_R2->Get(Form("PbPb_data_minbiasSub_cent%d",i));
    hPbPb_R2_measured[i]->Print("base");
    hPbPb_R3_measured[i] = (TH1F*)fError_R3->Get(Form("PbPb_data_minbiasSub_cent%d",i));
    hPbPb_R3_measured[i]->Print("base");
    hPbPb_R4_measured[i] = (TH1F*)fError_R4->Get(Form("PbPb_data_minbiasSub_cent%d",i));
    hPbPb_R4_measured[i]->Print("base");

    uPbPb_R2_Bayes[i] = (TH1F*)fin_R2->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R2_Bayes[i]->Print("base");
    uPbPb_R3_Bayes[i] = (TH1F*)fin_R3->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R3_Bayes[i]->Print("base");
    uPbPb_R4_Bayes[i] = (TH1F*)fin_R4->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R4_Bayes[i]->Print("base");
 
    mPbPb_R2[i] = (TH1F*)fin_R2->Get(Form("PbPb_Gen_spectra_refpt_cent%d",i));
    mPbPb_R2[i]->Print("base");
    mPbPb_R3[i] = (TH1F*)fin_R3->Get(Form("PbPb_Gen_spectra_refpt_cent%d",i));
    mPbPb_R3[i]->Print("base");
    mPbPb_R4[i] = (TH1F*)fin_R4->Get(Form("PbPb_Gen_spectra_refpt_cent%d",i));
    mPbPb_R4[i]->Print("base");
    
    RAA_R2_Bayes[i]   = (TH1F*)fin_R2->Get(Form("RAA_bayesian_cent%d",i));  
    RAA_R2_Bayes[i]->Print("base");
    RAA_R3_Bayes[i]   = (TH1F*)fin_R3->Get(Form("RAA_bayesian_cent%d",i));  
    RAA_R3_Bayes[i]->Print("base");
    RAA_R4_Bayes[i]   = (TH1F*)fin_R4->Get(Form("RAA_bayesian_cent%d",i));  
    RAA_R4_Bayes[i]->Print("base");
    
    // RAA_R2_Bayes[i]   = (TH1F*)RAA_R2_Bayes_fineBin[i]->Rebin(nbins_pt,Form("RAA_R2_analysisbins_bayesian_cent%d",i),boundaries_pt);  
    // divideBinWidth(RAA_R2_Bayes[i]);
    // RAA_R2_Bayes[i]->Print("base");
    // RAA_R3_Bayes[i]   = (TH1F*)RAA_R3_Bayes_fineBin[i]->Rebin(nbins_pt,Form("RAA_R3_analysisbins_bayesian_cent%d",i),boundaries_pt);  
    // divideBinWidth(RAA_R3_Bayes[i]);
    // RAA_R3_Bayes[i]->Print("base");
    // RAA_R4_Bayes[i]   = (TH1F*)RAA_R4_Bayes_fineBin[i]->Rebin(nbins_pt,Form("RAA_R4_analysisbins_bayesian_cent%d",i),boundaries_pt);  
    // divideBinWidth(RAA_R4_Bayes[i]);
    // RAA_R4_Bayes[i]->Print("base");

    // RAA_R2_Bayes_atlasBin[i]   = (TH1F*)RAA_R2_Bayes_fineBin[i]->Rebin(nbins_atlas,Form("RAA_R2_analysisbins_bayesian_cent%d",i),boundaries_atlas);  
    // divideBinWidth(RAA_R2_Bayes_atlasBin[i]);
    // RAA_R2_Bayes_atlasBin[i]->Print("base");
    // RAA_R3_Bayes_atlasBin[i]   = (TH1F*)RAA_R3_Bayes_fineBin[i]->Rebin(nbins_atlas,Form("RAA_R3_analysisbins_bayesian_cent%d",i),boundaries_atlas);  
    // divideBinWidth(RAA_R3_Bayes_atlasBin[i]);
    // RAA_R3_Bayes_atlasBin[i]->Print("base");
    // RAA_R4_Bayes_atlasBin[i]   = (TH1F*)RAA_R4_Bayes_fineBin[i]->Rebin(nbins_atlas,Form("RAA_R4_analysisbins_bayesian_cent%d",i),boundaries_atlas);  
    // divideBinWidth(RAA_R4_Bayes_atlasBin[i]);
    // RAA_R4_Bayes_atlasBin[i]->Print("base");

    // RAA_R2_Bayes_atlasRcpBin[i]   = (TH1F*)RAA_R2_Bayes_fineBin[i]->Rebin(nbins_atlasRcp,Form("RAA_R2_analysisbins_bayesian_cent%d",i),boundaries_atlasRcp);  
    // divideBinWidth(RAA_R2_Bayes_atlasRcpBin[i]);
    // RAA_R2_Bayes_atlasRcpBin[i]->Print("base");
    // RAA_R3_Bayes_atlasRcpBin[i]   = (TH1F*)RAA_R3_Bayes_fineBin[i]->Rebin(nbins_atlasRcp,Form("RAA_R3_analysisbins_bayesian_cent%d",i),boundaries_atlasRcp);  
    // divideBinWidth(RAA_R3_Bayes_atlasRcpBin[i]);
    // RAA_R3_Bayes_atlasRcpBin[i]->Print("base");
    // RAA_R4_Bayes_atlasRcpBin[i]   = (TH1F*)RAA_R4_Bayes_fineBin[i]->Rebin(nbins_atlasRcp,Form("RAA_R4_analysisbins_bayesian_cent%d",i),boundaries_atlasRcp);  
    // divideBinWidth(RAA_R4_Bayes_atlasRcpBin[i]);
    // RAA_R4_Bayes_atlasRcpBin[i]->Print("base");

    
    RAA_R2_BinByBin[i]   = (TH1F*)fin_R2->Get(Form("RAA_binbybin_cent%d",i));  
    RAA_R3_BinByBin[i]   = (TH1F*)fin_R3->Get(Form("RAA_binbybin_cent%d",i));  
    RAA_R4_BinByBin[i]   = (TH1F*)fin_R4->Get(Form("RAA_binbybin_cent%d",i));  
    
    RAA_R2_Meas[i]   = (TH1F*)fin_R2->Get(Form("RAA_measured_cent%d",i));  
    RAA_R3_Meas[i]   = (TH1F*)fin_R3->Get(Form("RAA_measured_cent%d",i));  
    RAA_R4_Meas[i]   = (TH1F*)fin_R4->Get(Form("RAA_measured_cent%d",i));  
    
  }
  
  // declare the systematics for the RAA
  SysData systematics_R2;
  SysData systematics_R3;
  SysData systematics_R4;

  // declare the systematics for the PbPb Spectra, and pp spectra
  // it turns out that they are already located in the systematics_R2 but im going to keep them seperate anyway 
  SysData systematics_PbPb_R2;
  SysData systematics_PbPb_R3;
  SysData systematics_PbPb_R4;

  SysData systematics_PP_R2;
  SysData systematics_PP_R3;
  SysData systematics_PP_R4;
  
  // ncoll uncertainty
  prepareNcollUnc(nbins_pt, 300.);
  
  // get the necessary histograms
  const int Iterations = 6; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_BayesianIter_R2[nbins_cent+1][Iterations];
  TH1F *uPP_BayesianIter_R2[Iterations];
  TH1F *uPbPb_BayesianIter_R3[nbins_cent+1][Iterations];
  TH1F *uPP_BayesianIter_R3[Iterations];
  TH1F *uPbPb_BayesianIter_R4[nbins_cent+1][Iterations];
  TH1F *uPP_BayesianIter_R4[Iterations];
  
  for(int j = 0;j<Iterations;j++){
    for(int i = 0;i<nbins_cent;++i){
      uPbPb_BayesianIter_R2[i][j] = (TH1F*)fin_R2->Get(Form("uPbPb_BayesianIter%d_cent%d",j+1,i));
      uPbPb_BayesianIter_R3[i][j] = (TH1F*)fin_R3->Get(Form("uPbPb_BayesianIter%d_cent%d",j+1,i));
      uPbPb_BayesianIter_R4[i][j] = (TH1F*)fin_R4->Get(Form("uPbPb_BayesianIter%d_cent%d",j+1,i));
    }
    uPP_BayesianIter_R2[j] = (TH1F*)fin_R2->Get(Form("uPP_BayesianIter%d",j+1));
    uPP_BayesianIter_R3[j] = (TH1F*)fin_R3->Get(Form("uPP_BayesianIter%d",j+1));
    uPP_BayesianIter_R4[j] = (TH1F*)fin_R4->Get(Form("uPP_BayesianIter%d",j+1));
  }
  
  //TH1F * uPbPb_Bayes[nbins_cent], * uPbPb_JEC_Bayes[nbins_cent], * uPbPb_Smear_Bayes[nbins_cent];
  //TH1F * uPP_Bayes, dPP_Bayes, dPP_BinByBin;
  TH1F *RAA_measured_R2[nbins_cent+1];
  TH1F *RAA_binbybin_R2[nbins_cent+1];
  TH1F *RAA_bayesian_R2[nbins_cent+1];
  TH1F *RAA_JEC_bayesian_R2[nbins_cent+1];
  TH1F *RAA_Smear_bayesian_R2[nbins_cent+1];
  
  TH1F *RAA_measured_R3[nbins_cent+1];
  TH1F *RAA_binbybin_R3[nbins_cent+1];
  TH1F *RAA_bayesian_R3[nbins_cent+1];
  TH1F *RAA_JEC_bayesian_R3[nbins_cent+1];
  TH1F *RAA_Smear_bayesian_R3[nbins_cent+1];
  
  TH1F *RAA_measured_R4[nbins_cent+1];
  TH1F *RAA_binbybin_R4[nbins_cent+1];
  TH1F *RAA_bayesian_R4[nbins_cent+1];
  TH1F *RAA_JEC_bayesian_R4[nbins_cent+1];
  TH1F *RAA_Smear_bayesian_R4[nbins_cent+1];
  

  for(int i = 0; i<nbins_cent; ++i){

    RAA_bayesian_R2[i] = (TH1F*)fin_R2->Get(Form("RAA_bayesian_cent%d",i));
    RAA_JEC_bayesian_R2[i] = (TH1F*)fin_R2->Get(Form("RAA_JEC_bayesian_cent%d",i));
    RAA_Smear_bayesian_R2[i] = (TH1F*)fin_R2->Get(Form("RAA_Smear_bayesian_cent%d",i));
    RAA_binbybin_R2[i] = (TH1F*)fin_R2->Get(Form("RAA_binbybin_cent%d",i));
    RAA_measured_R2[i] = (TH1F*)fin_R2->Get(Form("RAA_measured_cent%d",i));

    RAA_bayesian_R3[i] = (TH1F*)fin_R3->Get(Form("RAA_bayesian_cent%d",i));
    RAA_JEC_bayesian_R3[i] = (TH1F*)fin_R3->Get(Form("RAA_JEC_bayesian_cent%d",i));
    RAA_Smear_bayesian_R3[i] = (TH1F*)fin_R3->Get(Form("RAA_Smear_bayesian_cent%d",i));
    RAA_binbybin_R3[i] = (TH1F*)fin_R3->Get(Form("RAA_binbybin_cent%d",i));
    RAA_measured_R3[i] = (TH1F*)fin_R3->Get(Form("RAA_measured_cent%d",i));

    RAA_bayesian_R4[i] = (TH1F*)fin_R4->Get(Form("RAA_bayesian_cent%d",i));
    RAA_JEC_bayesian_R4[i] = (TH1F*)fin_R4->Get(Form("RAA_JEC_bayesian_cent%d",i));
    RAA_Smear_bayesian_R4[i] = (TH1F*)fin_R4->Get(Form("RAA_Smear_bayesian_cent%d",i));
    RAA_binbybin_R4[i] = (TH1F*)fin_R4->Get(Form("RAA_binbybin_cent%d",i));
    RAA_measured_R4[i] = (TH1F*)fin_R4->Get(Form("RAA_measured_cent%d",i));
    // uPbPb_Bayes[i] = (TH1F*)fin->Get(Form("PbPb_bayesian_unfolded_spectra_combined__cent%d",i));
    // uPbPb_Bayes[i]->Print("base");
  }
  
  // make the ratios for each iteration with respect to iteration 4
  TH1F * hDenominator_R2[nbins_cent];
  TH1F * hNumerator_R2[Iterations][nbins_cent];
  TH1F * hDenominatorPP_R2;
  TH1F * hNumeratorPP_R2[Iterations];
  
  TH1F * hDenominator_R3[nbins_cent];
  TH1F * hNumerator_R3[Iterations][nbins_cent];
  TH1F * hDenominatorPP_R3;
  TH1F * hNumeratorPP_R3[Iterations];

  TH1F * hDenominator_R4[nbins_cent];
  TH1F * hNumerator_R4[Iterations][nbins_cent];
  TH1F * hDenominatorPP_R4;
  TH1F * hNumeratorPP_R4[Iterations];
  
  for(int i = 0; i<nbins_cent; ++i){
    hDenominator_R2[i] = (TH1F*)uPbPb_BayesianIter_R2[i][BayesIter-1]->Clone(Form("Denominator_PbPb_R2_iter%d_cent%d",BayesIter,i));
    for(int j = 0; j<Iterations; ++j){
      hNumerator_R2[j][i] = (TH1F*)uPbPb_BayesianIter_R2[i][j]->Clone(Form("Numerator_PbPb_R2_iter%d_cent%d",j+1,i));
      hNumerator_R2[j][i]->Divide(hDenominator_R2[i]);
    }
  }

  hDenominatorPP_R2 = (TH1F*)uPP_BayesianIter_R2[BayesIter-1]->Clone(Form("Denominator_PP_R2_iter%d",BayesIter));
  for(int j = 0; j<Iterations; ++j){
    hNumeratorPP_R2[j] = (TH1F*)uPP_BayesianIter_R2[j]->Clone(Form("Numerator_PP_R2_iter%d",j+1));
    hNumeratorPP_R2[j]->Divide(hDenominatorPP_R2);
  }  

  for(int i = 0; i<nbins_cent; ++i){
    hDenominator_R3[i] = (TH1F*)uPbPb_BayesianIter_R3[i][BayesIter-1]->Clone(Form("Denominator_PbPb_R3_iter%d_cent%d",BayesIter,i));
    for(int j = 0; j<Iterations; ++j){
      hNumerator_R3[j][i] = (TH1F*)uPbPb_BayesianIter_R3[i][j]->Clone(Form("Numerator_PbPb_R3_iter%d_cent%d",j+1,i));
      hNumerator_R3[j][i]->Divide(hDenominator_R3[i]);
    }
  }

  hDenominatorPP_R3 = (TH1F*)uPP_BayesianIter_R3[BayesIter-1]->Clone(Form("Denominator_PP_R3_iter%d",BayesIter));
  for(int j = 0; j<Iterations; ++j){
    hNumeratorPP_R3[j] = (TH1F*)uPP_BayesianIter_R3[j]->Clone(Form("Numerator_PP_R3_iter%d",j+1));
    hNumeratorPP_R3[j]->Divide(hDenominatorPP_R3);
  }  

  for(int i = 0; i<nbins_cent; ++i){
    hDenominator_R4[i] = (TH1F*)uPbPb_BayesianIter_R4[i][BayesIter-1]->Clone(Form("Denominator_PbPb_R4_iter%d_cent%d",BayesIter,i));
    for(int j = 0; j<Iterations; ++j){
      hNumerator_R4[j][i] = (TH1F*)uPbPb_BayesianIter_R4[i][j]->Clone(Form("Numerator_PbPb_R4_iter%d_cent%d",j+1,i));
      hNumerator_R4[j][i]->Divide(hDenominator_R4[i]);
    }
  }

  hDenominatorPP_R4 = (TH1F*)uPP_BayesianIter_R4[BayesIter-1]->Clone(Form("Denominator_PP_R4_iter%d",BayesIter));
  for(int j = 0; j<Iterations; ++j){
    hNumeratorPP_R4[j] = (TH1F*)uPP_BayesianIter_R4[j]->Clone(Form("Numerator_PP_R4_iter%d",j+1));
    hNumeratorPP_R4[j]->Divide(hDenominatorPP_R4);
  }  

  // define the TLines for the JEC upper and lower bounds. 
  TLine * bounds_up[nbins_cent], *bounds_low[nbins_cent];
  bounds_up[0] = new TLine(65,1.12,299,1.12); 
  bounds_up[1] = new TLine(65,1.10,299,1.10); 
  bounds_up[2] = new TLine(65,1.08,299,1.08); 
  bounds_up[3] = new TLine(65,1.05,299,1.05); 
  bounds_up[4] = new TLine(65,1.03,299,1.03); 
  bounds_up[5] = new TLine(65,1.02,299,1.02); 
  
  bounds_low[0] = new TLine(65,0.88,299,0.88); 
  bounds_low[1] = new TLine(65,0.90,299,0.90); 
  bounds_low[2] = new TLine(65,0.92,299,0.92); 
  bounds_low[3] = new TLine(65,0.95,299,0.95); 
  bounds_low[4] = new TLine(65,0.97,299,0.97); 
  bounds_low[5] = new TLine(65,0.98,299,0.98); 
  
  for(int i = 0; i<nbins_cent; ++i){

    bounds_up[i]->SetLineWidth(2);
    bounds_up[i]->SetLineStyle(2);
    bounds_up[i]->SetLineColor(kRed);
    bounds_low[i]->SetLineWidth(2);
    bounds_low[i]->SetLineStyle(2);
    bounds_low[i]->SetLineColor(kRed);    

  }

  
  //make the canvas
  // line at 1
  TLine *linePbPb_iter = new TLine(65,1,299,1);
  linePbPb_iter->SetLineStyle(2);
  linePbPb_iter->SetLineWidth(2);
  

  TCanvas * cPbPb_Its_R2 = new TCanvas("cPbPb_Its_R2","",1200,1000);
  makeMultiPanelCanvas(cPbPb_Its_R2,3,2,0.0,0.0,0.2,0.15,0.07);
  
  TLegend *PbPb_itersys_R2 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R2->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R2[j][i]->SetMarkerStyle(33);
      hNumerator_R2[j][i]->SetMarkerColor(j+1);
      hNumerator_R2[j][i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
      hNumerator_R2[j][i]->SetAxisRange(0,2,"Y");

      if(j==1){
	makeHistTitle(hNumerator_R2[j][i]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
	hNumerator_R2[j][i]->Draw();
      }
      hNumerator_R2[j][i]->Draw("same");
	
      if(i==0) PbPb_itersys_R2->AddEntry(hNumerator_R2[j][i],Form("Iteration %d",j+1),"pl");
      checkMaximumSys(systematics_R2.hSysIter[i], hNumerator_R2[j][i], 0, 1.05);
      checkMaximumSys(systematics_PbPb_R2.hSysIter[i], hNumerator_R2[j][i], 0, 1.05);
      
    }

    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    linePbPb_iter->Draw();
    drawEnvelope(systematics_R2.hSysIter[i], "hist same");

  }
  PbPb_itersys_R2->Draw();

  cPbPb_Its_R2->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, 2),0.2,0.23,20);
  drawText(Form("|#eta|< %2.0f, |vz|<15",etaBoundary),0.6,0.31,20);

  cPbPb_Its_R2->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR2%s_%d_%s_pawan_ntuple.pdf",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPbPb_Its_R2->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR2%s_%d_%s_pawan_ntuple.C",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPbPb_Its_R2->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR2%s_%d_%s_pawan_ntuple.root",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  

  TCanvas * cPP_Its_R2 = new TCanvas("cPP_Its_R2","",800,600);
  // line at 1
  TLine *linePP_iter_R2 = new TLine(65,1,299,1);
  linePP_iter_R2->SetLineStyle(2);
  linePP_iter_R2->SetLineWidth(2);
  
  TLegend *PP_itersys_R2 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int j = 1;j<Iterations;j++){
    hNumeratorPP_R2[j]->SetMarkerStyle(33);
    hNumeratorPP_R2[j]->SetMarkerColor(j+1);
    hNumeratorPP_R2[j]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    hNumeratorPP_R2[j]->SetAxisRange(0,2,"Y");

    if(j==1){
      makeHistTitle(hNumeratorPP_R2[j]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
      hNumeratorPP_R2[j]->Draw();
    }
    hNumeratorPP_R2[j]->Draw("same");

    checkMaximumSys(systematics_R2.hSysIter[nbins_cent], hNumeratorPP_R2[j],0, 1.05);
    checkMaximumSys(systematics_PP_R2.hSysIter[nbins_cent], hNumeratorPP_R2[j],0, 1.05);
    PP_itersys_R2->AddEntry(hNumeratorPP_R2[j],Form("Iteration %d",j+1),"pl");

  }

  linePP_iter_R2->Draw();  
  PP_itersys_R2->Draw();
  drawEnvelope(systematics_R2.hSysIter[nbins_cent],"hist same");

  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, 2),0.2,0.23,20);
  drawText(Form("|#eta|< %2.0f, |vz|<15",etaBoundary),0.6,0.31,20);

  cPP_Its_R2->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR2%s_%d_%s_pawan_ntuple.pdf",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPP_Its_R2->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR2%s_%d_%s_pawan_ntuple.C",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPP_Its_R2->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR2%s_%d_%s_pawan_ntuple.root",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");

  // define the Jet ID efficiency as a function of jet pT: taken from the above plots: 
  TF1 * fPol_x2 = new TF1("fPol_x2","1-[0]/pow(x,[1])");
  TF1 * fPol_x1 = new TF1("fPol_x1","[0]+[1]*x");
  
  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R2 = new TCanvas("cRAA_JEC_sys_R2","",1200,1000);
  makeMultiPanelCanvas(cRAA_JEC_sys_R2,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R2->cd(nbins_cent-i);
    RAA_JEC_bayesian_R2[i]->Divide(RAA_JEC_bayesian_R2[i], RAA_bayesian_R2[i], 1,1, "B");
    RAA_JEC_bayesian_R2[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R2[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R2[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1], "X");
    //RAA_JEC_bayesian_R2[i]->Fit("fPol_x2","","",60,299);


    RAA_JEC_bayesian_R2[i]->Draw("p");
    //checkMaximumSys(systematics_R2.hSysJEC[i], functionHist(fPol_x2, systematics_R2.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();

    bounds_up[i]->Draw();
    //bounds_low[i]->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }

  cRAA_JEC_sys_R2->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R0.2",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|< %2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_JEC_sys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRAA_JEC_sys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRAA_JEC_sys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

  // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R2 = new TCanvas("cRAA_Smear_sys_R2","",1200,1000);
  makeMultiPanelCanvas(cRAA_Smear_sys_R2,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R2->cd(nbins_cent-i);
    RAA_Smear_bayesian_R2[i]->Divide(RAA_Smear_bayesian_R2[i], RAA_bayesian_R2[i], 1,1, "B");
    RAA_Smear_bayesian_R2[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R2[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R2[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1], "X");
    RAA_Smear_bayesian_R2[i]->Fit("fPol_x1","","",60,299);
    RAA_Smear_bayesian_R2[i]->Draw("p");
    checkMaximumSys(systematics_R2.hSysSmear[i], functionHist(fPol_x1, systematics_R2.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));
    checkMaximumSys(systematics_PbPb_R2.hSysSmear[i], functionHist(fPol_x1, systematics_R2.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }

  cRAA_Smear_sys_R2->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.2",algo, jet_type),0.2,0.23,20);
  if(etaWidth=="n16_eta_p16")drawText("|#eta|<1.6, |vz|<15",0.6,0.31,20);
  if(etaWidth=="n20_eta_p20")drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_Smear_sys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));
  cRAA_Smear_sys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R2%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));
  cRAA_Smear_sys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R2%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));


  // draw the total systematics
  TCanvas * cSys_R2 = new TCanvas("cSys_R2","Total Systematics",1200,1200);
  makeMultiPanelCanvas(cSys_R2,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0; i<nbins_cent; ++i){

    cSys_R2->cd(nbins_cent-i);
    cSys_R2->cd(nbins_cent-i)->SetGridy();
    cSys_R2->cd(nbins_cent-i)->SetGridx();
    systematics_R2.DrawComponent(i);

    //drawPanelLabel(i);
    //TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    //title_->AddEntry(RAA_bayesian_R2[i],Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    //title_->SetTextSize(0.06);
    //title_->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  cSys_R2->cd(1);
  putCMSPrel();
  linePbPb_iter->Draw();
  cSys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cSys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cSys_R2->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");


  TCanvas * cPbPb_Its_R3 = new TCanvas("cPbPb_Its_R3","",1200,1000);
  makeMultiPanelCanvas(cPbPb_Its_R3,3,2,0.0,0.0,0.2,0.15,0.07);
  
  TLegend *PbPb_itersys_R3 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R3->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R3[j][i]->SetMarkerStyle(33);
      hNumerator_R3[j][i]->SetMarkerColor(j+1);
      hNumerator_R3[j][i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
      hNumerator_R3[j][i]->SetAxisRange(0,2,"Y");

      if(j==1){
	makeHistTitle(hNumerator_R3[j][i]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
	hNumerator_R3[j][i]->Draw();
      }
      hNumerator_R3[j][i]->Draw("same");
	
      if(i==0) PbPb_itersys_R3->AddEntry(hNumerator_R3[j][i],Form("Iteration %d",j+1),"pl");
      checkMaximumSys(systematics_R3.hSysIter[i], hNumerator_R3[j][i], 0, 1.05);
      
    }

    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    linePbPb_iter->Draw();
    drawEnvelope(systematics_R3.hSysIter[i], "hist same");

  }
  PbPb_itersys_R3->Draw();

  cPbPb_Its_R3->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, 3),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  
  cPbPb_Its_R3->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR3%s_%d_%s_pawan_ntuple.pdf",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPbPb_Its_R3->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR3%s_%d_%s_pawan_ntuple.C",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPbPb_Its_R3->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR3%s_%d_%s_pawan_ntuple.root",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  

  TCanvas * cPP_Its_R3 = new TCanvas("cPP_Its_R3","",800,600);
  // line at 1
  TLine *linePP_iter_R3 = new TLine(65,1,299,1);
  linePP_iter_R3->SetLineStyle(2);
  linePP_iter_R3->SetLineWidth(2);
  
  TLegend *PP_itersys_R3 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int j = 1;j<Iterations;j++){
    hNumeratorPP_R3[j]->SetMarkerStyle(33);
    hNumeratorPP_R3[j]->SetMarkerColor(j+1);
    hNumeratorPP_R3[j]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    hNumeratorPP_R3[j]->SetAxisRange(0,2,"Y");

    if(j==1){
      makeHistTitle(hNumeratorPP_R3[j]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
      hNumeratorPP_R3[j]->Draw();
    }
    hNumeratorPP_R3[j]->Draw("same");

    checkMaximumSys(systematics_R3.hSysIter[nbins_cent], hNumeratorPP_R3[j],0, 1.05);
    PP_itersys_R3->AddEntry(hNumeratorPP_R3[j],Form("Iteration %d",j+1),"pl");

  }

  linePP_iter_R3->Draw();  
  PP_itersys_R3->Draw();
  drawEnvelope(systematics_R3.hSysIter[nbins_cent],"hist same");

  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, 3),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  
  cPP_Its_R3->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR3%s_%d_%s_pawan_ntuple.pdf",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPP_Its_R3->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR3%s_%d_%s_pawan_ntuple.C",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPP_Its_R3->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR3%s_%d_%s_pawan_ntuple.root",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");

  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R3 = new TCanvas("cRAA_JEC_sys_R3","",1200,1000);
  makeMultiPanelCanvas(cRAA_JEC_sys_R3,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R3->cd(nbins_cent-i);
    RAA_JEC_bayesian_R3[i]->Divide(RAA_JEC_bayesian_R3[i], RAA_bayesian_R3[i], 1,1, "B");
    RAA_JEC_bayesian_R3[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R3[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R3[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1], "X");
    //RAA_JEC_bayesian_R3[i]->Fit("fPol_x2","","",60,299);

    RAA_JEC_bayesian_R3[i]->Draw("p");
    //checkMaximumSys(systematics_R3.hSysJEC[i], functionHist(fPol_x2, systematics_R3.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

    bounds_up[i]->Draw();
    //bounds_low[i]->Draw();

  }

  cRAA_JEC_sys_R3->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.3",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_JEC_sys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRAA_JEC_sys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRAA_JEC_sys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

  // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R3 = new TCanvas("cRAA_Smear_sys_R3","",1200,1000);
  makeMultiPanelCanvas(cRAA_Smear_sys_R3,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R3->cd(nbins_cent-i);
    RAA_Smear_bayesian_R3[i]->Divide(RAA_Smear_bayesian_R3[i], RAA_bayesian_R3[i], 1,1, "B");
    RAA_Smear_bayesian_R3[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R3[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R3[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1], "X");
    RAA_Smear_bayesian_R3[i]->Fit("fPol_x1","","",60,299);
    RAA_Smear_bayesian_R3[i]->Draw("p");
    checkMaximumSys(systematics_R3.hSysSmear[i], functionHist(fPol_x1, systematics_R3.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }

  cRAA_Smear_sys_R3->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.3",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_Smear_sys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));
  cRAA_Smear_sys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R3%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));
  cRAA_Smear_sys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R3%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));


  // draw the total systematics
  TCanvas * cSys_R3 = new TCanvas("cSys_R3","Total Systematics",1200,1200);
  makeMultiPanelCanvas(cSys_R3,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0; i<nbins_cent; ++i){

    cSys_R3->cd(nbins_cent-i);
    cSys_R3->cd(nbins_cent-i)->SetGridy();
    cSys_R3->cd(nbins_cent-i)->SetGridx();
    systematics_R3.DrawComponent(i);

    //drawPanelLabel(i);
    //TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    //title_->AddEntry(RAA_bayesian_R3[i],Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    //title_->SetTextSize(0.06);
    //title_->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  cSys_R3->cd(1);
  putCMSPrel();
  linePbPb_iter->Draw();
  cSys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cSys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cSys_R3->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");



  TCanvas * cPbPb_Its_R4 = new TCanvas("cPbPb_Its_R4","",1200,1000);
  makeMultiPanelCanvas(cPbPb_Its_R4,3,2,0.0,0.0,0.2,0.15,0.07);
  
  TLegend *PbPb_itersys_R4 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R4->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R4[j][i]->SetMarkerStyle(33);
      hNumerator_R4[j][i]->SetMarkerColor(j+1);
      hNumerator_R4[j][i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
      hNumerator_R4[j][i]->SetAxisRange(0,2,"Y");

      if(j==1){
	makeHistTitle(hNumerator_R4[j][i]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
	hNumerator_R4[j][i]->Draw();
      }
      hNumerator_R4[j][i]->Draw("same");
	
      if(i==0) PbPb_itersys_R4->AddEntry(hNumerator_R4[j][i],Form("Iteration %d",j+1),"pl");
      checkMaximumSys(systematics_R4.hSysIter[i], hNumerator_R4[j][i], 0, 1.05);
      
    }

    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    linePbPb_iter->Draw();
    drawEnvelope(systematics_R4.hSysIter[i], "hist same");

  }
  PbPb_itersys_R4->Draw();

  cPbPb_Its_R4->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, 4),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cPbPb_Its_R4->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR4%s_%d_%s_pawan_ntuple.pdf",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPbPb_Its_R4->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR4%s_%d_%s_pawan_ntuple.C",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPbPb_Its_R4->SaveAs(Form("June4/Analysis_full_bins/PbPb_unfoldng_iterations_%d_systematics_%dGeVCut_ak%sR4%s_%d_%s_pawan_ntuple.root",Iterations, unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  

  TCanvas * cPP_Its_R4 = new TCanvas("cPP_Its_R4","",800,600);
  // line at 1
  TLine *linePP_iter_R4 = new TLine(65,1,299,1);
  linePP_iter_R4->SetLineStyle(2);
  linePP_iter_R4->SetLineWidth(2);
  
  TLegend *PP_itersys_R4 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int j = 1;j<Iterations;j++){
    hNumeratorPP_R4[j]->SetMarkerStyle(33);
    hNumeratorPP_R4[j]->SetMarkerColor(j+1);
    hNumeratorPP_R4[j]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    hNumeratorPP_R4[j]->SetAxisRange(0,2,"Y");

    if(j==1){
      makeHistTitle(hNumeratorPP_R4[j]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
      hNumeratorPP_R4[j]->Draw();
    }
    hNumeratorPP_R4[j]->Draw("same");

    checkMaximumSys(systematics_R4.hSysIter[nbins_cent], hNumeratorPP_R4[j],0, 1.05);
    PP_itersys_R4->AddEntry(hNumeratorPP_R4[j],Form("Iteration %d",j+1),"pl");

  }

  linePP_iter_R4->Draw();  
  PP_itersys_R4->Draw();
  drawEnvelope(systematics_R4.hSysIter[nbins_cent],"hist same");

  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, 4),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cPP_Its_R4->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR4%s_%d_%s_pawan_ntuple.pdf",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPP_Its_R4->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR4%s_%d_%s_pawan_ntuple.C",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");
  cPP_Its_R4->SaveAs(Form("June4/Analysis_full_bins/PP_unfoldng_iterations_%d_systematics_%dGeVCut_akR4%s_%d_%s_pawan_ntuple.root",Iterations, unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");

  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R4 = new TCanvas("cRAA_JEC_sys_R4","",1200,1000);
  makeMultiPanelCanvas(cRAA_JEC_sys_R4,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R4->cd(nbins_cent-i);
    RAA_JEC_bayesian_R4[i]->Divide(RAA_JEC_bayesian_R4[i], RAA_bayesian_R4[i], 1,1, "B");
    RAA_JEC_bayesian_R4[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R4[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R4[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1], "X");
    //RAA_JEC_bayesian_R4[i]->Fit("fPol_x2","","",60,299);

    RAA_JEC_bayesian_R4[i]->Draw("p");
    //checkMaximumSys(systematics_R4.hSysJEC[i], functionHist(fPol_x2, systematics_R4.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();

    bounds_up[i]->Draw();
    //bounds_low[i]->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }

  cRAA_JEC_sys_R4->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.4",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_JEC_sys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRAA_JEC_sys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRAA_JEC_sys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

  // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R4 = new TCanvas("cRAA_Smear_sys_R4","",1200,1000);
  makeMultiPanelCanvas(cRAA_Smear_sys_R4,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R4->cd(nbins_cent-i);
    RAA_Smear_bayesian_R4[i]->Divide(RAA_Smear_bayesian_R4[i], RAA_bayesian_R4[i], 1,1, "B");
    RAA_Smear_bayesian_R4[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R4[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R4[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1], "X");
    RAA_Smear_bayesian_R4[i]->Fit("fPol_x1","","",60,299);
    RAA_Smear_bayesian_R4[i]->Draw("p");
    checkMaximumSys(systematics_R4.hSysSmear[i], functionHist(fPol_x1, systematics_R4.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }

  cRAA_Smear_sys_R4->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.4",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_Smear_sys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));
  cRAA_Smear_sys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R4%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));
  cRAA_Smear_sys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R4%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));


  // draw the total systematics
  TCanvas * cSys_R4 = new TCanvas("cSys_R4","Total Systematics",1200,1200);
  makeMultiPanelCanvas(cSys_R4,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0; i<nbins_cent; ++i){

    cSys_R4->cd(nbins_cent-i);
    cSys_R4->cd(nbins_cent-i)->SetGridy();
    cSys_R4->cd(nbins_cent-i)->SetGridx();
    systematics_R4.DrawComponent(i);

    //drawPanelLabel(i);
    //TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    //title_->AddEntry(RAA_bayesian_R4[i],Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    //title_->SetTextSize(0.06);
    //title_->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  cSys_R4->cd(1);
  putCMSPrel();
  linePbPb_iter->Draw();
  cSys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cSys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cSys_R4->SaveAs(Form("June4/Analysis_full_bins/RAA_total_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");


  // ATLAS Rcp points: 

  // Plot: p8170_d2x1y4 Rcp R=0.2, 0-10% 
  double p8170_d2x1y1_xval[] = { 41.285, 47.575, 54.82, 63.17, 72.79, 83.875, 96.655, 111.4, 128.35, 
    147.85, 170.4, 196.4 };
  double p8170_d2x1y1_xerrminus[] = { 2.924999999999997, 3.365000000000002, 3.8800000000000026, 4.469999999999999, 5.150000000000006, 5.935000000000002, 6.844999999999999, 7.900000000000006, 9.049999999999997, 
    10.449999999999989, 12.099999999999994, 13.900000000000006 };
  double p8170_d2x1y1_xerrplus[] = { 2.9250000000000043, 3.364999999999995, 3.8800000000000026, 4.469999999999999, 5.1499999999999915, 5.935000000000002, 6.844999999999999, 7.8999999999999915, 9.050000000000011, 
    10.450000000000017, 12.099999999999994, 13.900000000000006 };
  double p8170_d2x1y1_yval[] = { 0.4704, 0.4343, 0.4209, 0.4396, 0.4619, 0.481, 0.4918, 0.4962, 0.4937, 
    0.4939, 0.4908, 0.4872 };
  double p8170_d2x1y1_yerrminus[] = { 0.10432473340488342, 0.05235350991098878, 0.049562889342733035, 0.05282508873631923, 0.05807934228277728, 0.0634736165662553, 0.06800889647685808, 0.07199187454150642, 0.07530166000826277, 
    0.07942128178265571, 0.08351850094440154, 0.08881756583018925 };
  double p8170_d2x1y1_yerrplus[] = { 0.08469226647103029, 0.056619607911040856, 0.057403222906035514, 0.05767616492104863, 0.06102220907177976, 0.0657116428039963, 0.06940244952449444, 0.06852342665103665, 0.07217818229908536, 
    0.07699045395372078, 0.08104936767180852, 0.0844627728647361 };
  double p8170_d2x1y1_ystatminus[] = { 0.0114, 0.0112, 0.0112, 0.0129, 0.0161, 0.0208, 0.027, 0.0345, 0.042, 
    0.0497, 0.0562, 0.0614 };
  double p8170_d2x1y1_ystatplus[] = { 0.0114, 0.0112, 0.0112, 0.0129, 0.0161, 0.0208, 0.027, 0.0345, 0.042, 
    0.0497, 0.0562, 0.0614 };
  int p8170_d2x1y1_numpoints = 12;
  TGraphAsymmErrors p8170_d2x1y1(p8170_d2x1y1_numpoints, p8170_d2x1y1_xval, p8170_d2x1y1_yval, p8170_d2x1y1_xerrminus, p8170_d2x1y1_xerrplus, p8170_d2x1y1_yerrminus, p8170_d2x1y1_yerrplus);
  p8170_d2x1y1.SetName("/HepData/8170/d2x1y1");
  p8170_d2x1y1.SetTitle("/HepData/8170/d2x1y1");

  // Rcp R=0.3, 0-10% 
  double p8170_d2x1y2_xval[] = { 41.285, 47.575, 54.82, 63.17, 72.79, 83.875, 96.655, 111.4, 128.35, 
    147.85, 170.4, 196.4 };
  double p8170_d2x1y2_xerrminus[] = { 2.924999999999997, 3.365000000000002, 3.8800000000000026, 4.469999999999999, 5.150000000000006, 5.935000000000002, 6.844999999999999, 7.900000000000006, 9.049999999999997, 
    10.449999999999989, 12.099999999999994, 13.900000000000006 };
  double p8170_d2x1y2_xerrplus[] = { 2.9250000000000043, 3.364999999999995, 3.8800000000000026, 4.469999999999999, 5.1499999999999915, 5.935000000000002, 6.844999999999999, 7.8999999999999915, 9.050000000000011, 
    10.450000000000017, 12.099999999999994, 13.900000000000006 };
  double p8170_d2x1y2_yval[] = { 0.5157, 0.4779, 0.458, 0.4595, 0.4749, 0.4919, 0.5027, 0.5098, 0.5115, 
    0.5129, 0.5132, 0.5122 };
  double p8170_d2x1y2_yerrminus[] = { 0.10010319675215172, 0.08186946927884656, 0.06598014852969035, 0.0605403171448581, 0.061047768181973695, 0.06320711985211792, 0.06533153909100871, 0.06801602752293022, 0.07133764784459885, 
    0.07709597914288396, 0.08411991440794504, 0.09128143294230213 };
  double p8170_d2x1y2_yerrplus[] = { 0.1029421682305167, 0.08473051398404237, 0.06911881075365808, 0.06372142496837307, 0.0641477980915947, 0.066283180370287, 0.06869308553267935, 0.06712868239433871, 0.07115173926194637, 
    0.0761998687662912, 0.0817672917492074, 0.08721960788721765 };
  double p8170_d2x1y2_ystatminus[] = { 0.0128, 0.0119, 0.0115, 0.0121, 0.0143, 0.018, 0.0232, 0.0303, 0.0384, 
    0.0471, 0.0556, 0.063 };
  double p8170_d2x1y2_ystatplus[] = { 0.0128, 0.0119, 0.0115, 0.0121, 0.0143, 0.018, 0.0232, 0.0303, 0.0384, 
    0.0471, 0.0556, 0.063 };
  int p8170_d2x1y2_numpoints = 12;
  TGraphAsymmErrors p8170_d2x1y2(p8170_d2x1y2_numpoints, p8170_d2x1y2_xval, p8170_d2x1y2_yval, p8170_d2x1y2_xerrminus, p8170_d2x1y2_xerrplus, p8170_d2x1y2_yerrminus, p8170_d2x1y2_yerrplus);
  p8170_d2x1y2.SetName("/HepData/8170/d2x1y2");
  p8170_d2x1y2.SetTitle("/HepData/8170/d2x1y2");

  // Rcp R=0.4, 0-10%
  double p8170_d2x1y3_xval[] = { 41.285, 47.575, 54.82, 63.17, 72.79, 83.875, 96.655, 111.4, 128.35, 
    147.85, 170.4, 196.4 };
  double p8170_d2x1y3_xerrminus[] = { 2.924999999999997, 3.365000000000002, 3.8800000000000026, 4.469999999999999, 5.150000000000006, 5.935000000000002, 6.844999999999999, 7.900000000000006, 9.049999999999997, 
    10.449999999999989, 12.099999999999994, 13.900000000000006 };
  double p8170_d2x1y3_xerrplus[] = { 2.9250000000000043, 3.364999999999995, 3.8800000000000026, 4.469999999999999, 5.1499999999999915, 5.935000000000002, 6.844999999999999, 7.8999999999999915, 9.050000000000011, 
    10.450000000000017, 12.099999999999994, 13.900000000000006 };
  double p8170_d2x1y3_yval[] = { 0.528, 0.5556, 0.5594, 0.5527, 0.5502, 0.549, 0.5398, 0.5324, 0.5199, 
    0.5096, 0.4987, 0.4888 };
  double p8170_d2x1y3_yerrminus[] = { 0.10688096182201955, 0.10246155376530262, 0.08521772116173959, 0.07754282687650742, 0.06963404339832636, 0.06771078200700388, 0.06720892797835716, 0.06921719439561243, 0.0723074684939253, 
    0.07731836780481077, 0.08354244430228265, 0.09071460742350154 };
  double p8170_d2x1y3_yerrplus[] = { 0.11300769885277728, 0.1054047437262669, 0.09213500963260383, 0.08474184326529605, 0.07993209618169661, 0.07867127811342588, 0.07771312630437667, 0.07457921962584484, 0.07786334182399314, 
    0.08311732671350781, 0.08935636519017545, 0.09581429955909504 };
  double p8170_d2x1y3_ystatminus[] = { 0.0107, 0.013, 0.0137, 0.0142, 0.016, 0.0193, 0.0238, 0.0305, 0.0382, 
    0.0464, 0.0542, 0.0609 };
  double p8170_d2x1y3_ystatplus[] = { 0.0107, 0.013, 0.0137, 0.0142, 0.016, 0.0193, 0.0238, 0.0305, 0.0382, 
    0.0464, 0.0542, 0.0609 };
  int p8170_d2x1y3_numpoints = 12;
  TGraphAsymmErrors p8170_d2x1y3(p8170_d2x1y3_numpoints, p8170_d2x1y3_xval, p8170_d2x1y3_yval, p8170_d2x1y3_xerrminus, p8170_d2x1y3_xerrplus, p8170_d2x1y3_yerrminus, p8170_d2x1y3_yerrplus);
  p8170_d2x1y3.SetName("/HepData/8170/d2x1y3");
  p8170_d2x1y3.SetTitle("/HepData/8170/d2x1y3");


  // these ones are for Rcp (0.3) / Rcp(0.2) 0-10% 
  // Plot: p8170_d44x1y3
  double p8170_d44x1y1_xval[] = { 41.285, 47.575, 54.82, 63.17, 72.79, 83.875, 96.655, 111.4, 128.35, 
    147.85, 170.4, 196.4 };
  double p8170_d44x1y1_xerrminus[] = { 2.924999999999997, 3.365000000000002, 3.8800000000000026, 4.469999999999999, 5.150000000000006, 5.935000000000002, 6.844999999999999, 7.900000000000006, 9.049999999999997, 
    10.449999999999989, 12.099999999999994, 13.900000000000006 };
  double p8170_d44x1y1_xerrplus[] = { 2.9250000000000043, 3.364999999999995, 3.8800000000000026, 4.469999999999999, 5.1499999999999915, 5.935000000000002, 6.844999999999999, 7.8999999999999915, 9.050000000000011, 
    10.450000000000017, 12.099999999999994, 13.900000000000006 };
  double p8170_d44x1y1_yval[] = { 1.096, 1.1, 1.088, 1.045, 1.028, 1.023, 1.022, 1.028, 1.036, 
    1.039, 1.045, 1.051 };
  double p8170_d44x1y1_yerrminus[] = { 0.1582234495894967, 0.10879963235232, 0.08648098056798385, 0.05816201509576504, 0.050735687637007545, 0.0560280286999284, 0.06196499011538693, 0.06764066824034191, 0.0761121540885554, 
    0.08579026751327914, 0.10094622330726395, 0.12271470979471044 };
  double p8170_d44x1y1_yerrplus[] = { 0.04773237894762841, 0.10873237788257921, 0.10165008607964875, 0.06138517736392068, 0.040830135929237364, 0.04144128376389901, 0.04421244168783262, 0.05035136542339244, 0.06193367420071249, 
    0.0754761551749955, 0.0897681457979388, 0.10377485244508904 };
  double p8170_d44x1y1_ystatminus[] = { 0.0221, 0.0194, 0.0188, 0.0193, 0.0223, 0.0277, 0.0359, 0.0474, 0.0607, 
    0.0744, 0.0882, 0.1005 };
  double p8170_d44x1y1_ystatplus[] = { 0.0221, 0.0194, 0.0188, 0.0193, 0.0223, 0.0277, 0.0359, 0.0474, 0.0607, 
    0.0744, 0.0882, 0.1005 };
  int p8170_d44x1y1_numpoints = 12;
  TGraphAsymmErrors p8170_d44x1y1(p8170_d44x1y1_numpoints, p8170_d44x1y1_xval, p8170_d44x1y1_yval, p8170_d44x1y1_xerrminus, p8170_d44x1y1_xerrplus, p8170_d44x1y1_yerrminus, p8170_d44x1y1_yerrplus);
  p8170_d44x1y1.SetName("/HepData/8170/d44x1y1");
  p8170_d44x1y1.SetTitle("/HepData/8170/d44x1y1");

  // these ones are for Rcp (0.4) / Rcp(0.2) 0-10% 
  double p8170_d44x1y2_xval[] = { 41.285, 47.575, 54.82, 63.17, 72.79, 83.875, 96.655, 111.4, 128.35, 
    147.85, 170.4, 196.4 };
  double p8170_d44x1y2_xerrminus[] = { 2.924999999999997, 3.365000000000002, 3.8800000000000026, 4.469999999999999, 5.150000000000006, 5.935000000000002, 6.844999999999999, 7.900000000000006, 9.049999999999997, 
    10.449999999999989, 12.099999999999994, 13.900000000000006 };
  double p8170_d44x1y2_xerrplus[] = { 2.9250000000000043, 3.364999999999995, 3.8800000000000026, 4.469999999999999, 5.1499999999999915, 5.935000000000002, 6.844999999999999, 7.8999999999999915, 9.050000000000011, 
    10.450000000000017, 12.099999999999994, 13.900000000000006 };
  double p8170_d44x1y2_yval[] = { 1.122, 1.279, 1.329, 1.257, 1.191, 1.141, 1.098, 1.073, 1.053, 
    1.032, 1.016, 1.003 };
  double p8170_d44x1y2_yerrminus[] = { 0.1643743288959684, 0.14402090126089337, 0.12131957797486768, 0.0759939471273864, 0.05351382998814419, 0.057808476887044866, 0.06572853261712147, 0.07394112522811645, 0.08437825549275121, 
    0.09866215079755762, 0.11775852410759911, 0.141549037439327 };
  double p8170_d44x1y2_yerrplus[] = { 0.06375186271788456, 0.14826290837562847, 0.15079502644318213, 0.1100721581509148, 0.08116384663136661, 0.07695316757613035, 0.07434762941748714, 0.0772503721673883, 0.08790955579457788, 
    0.10392203808625002, 0.12225898739969998, 0.1402864569372254 };
  double p8170_d44x1y2_ystatminus[] = { 0.0254, 0.0269, 0.026, 0.0248, 0.0263, 0.0307, 0.038, 0.0493, 0.0632, 
    0.0785, 0.0937, 0.1072 };
  double p8170_d44x1y2_ystatplus[] = { 0.0254, 0.0269, 0.026, 0.0248, 0.0263, 0.0307, 0.038, 0.0493, 0.0632, 
    0.0785, 0.0937, 0.1072 };
  int p8170_d44x1y2_numpoints = 12;
  TGraphAsymmErrors p8170_d44x1y2(p8170_d44x1y2_numpoints, p8170_d44x1y2_xval, p8170_d44x1y2_yval, p8170_d44x1y2_xerrminus, p8170_d44x1y2_xerrplus, p8170_d44x1y2_yerrminus, p8170_d44x1y2_yerrplus);
  p8170_d44x1y2.SetName("/HepData/8170/d44x1y2");
  p8170_d44x1y2.SetTitle("/HepData/8170/d44x1y2");

  // these ones are for Rcp (0.5) / Rcp(0.2) 0-10% 
  // double p8170_d44x1y3_xval[] = { 41.285, 47.575, 54.82, 63.17, 72.79, 83.875, 96.655, 111.4, 128.35, 
  //   147.85, 170.4, 196.4 };
  // double p8170_d44x1y3_xerrminus[] = { 2.924999999999997, 3.365000000000002, 3.8800000000000026, 4.469999999999999, 5.150000000000006, 5.935000000000002, 6.844999999999999, 7.900000000000006, 9.049999999999997, 
  //   10.449999999999989, 12.099999999999994, 13.900000000000006 };
  // double p8170_d44x1y3_xerrplus[] = { 2.9250000000000043, 3.364999999999995, 3.8800000000000026, 4.469999999999999, 5.1499999999999915, 5.935000000000002, 6.844999999999999, 7.8999999999999915, 9.050000000000011, 
  //   10.450000000000017, 12.099999999999994, 13.900000000000006 };
  // double p8170_d44x1y3_yval[] = { 1.186, 1.331, 1.471, 1.452, 1.415, 1.358, 1.307, 1.269, 1.243, 
  //   1.219, 1.201, 1.188 };
  // double p8170_d44x1y3_yerrminus[] = { 0.20627045353128015, 0.19308073440921028, 0.18799029230255482, 0.15606207098459254, 0.13660179354605853, 0.10531695020270954, 0.08157070552594234, 0.08430996382397515, 0.09882580634631827, 
  //   0.1297209697774419, 0.16574061059378295, 0.20571847267564478 };
  // double p8170_d44x1y3_yerrplus[] = { 0.1749122637209867, 0.23427831739194302, 0.2676643794007712, 0.25556104945785457, 0.22740784946874637, 0.19165432424028422, 0.1711231720136113, 0.1616201720083233, 0.15815574602271015, 
  //   0.16805493149562734, 0.19107197596717315, 0.22119079998951133 };
  // double p8170_d44x1y3_ystatminus[] = { 0.0314, 0.0345, 0.0345, 0.0318, 0.0326, 0.0371, 0.0451, 0.0566, 0.0702, 
  //   0.0846, 0.0989, 0.112 };
  // double p8170_d44x1y3_ystatplus[] = { 0.0314, 0.0345, 0.0345, 0.0318, 0.0326, 0.0371, 0.0451, 0.0566, 0.0702, 
  //   0.0846, 0.0989, 0.112 };
  // int p8170_d44x1y3_numpoints = 12;
  // p8170_d44x1y3 = TGraphAsymmErrors(p8170_d44x1y3_numpoints, p8170_d44x1y3_xval, p8170_d44x1y3_yval, p8170_d44x1y3_xerrminus, p8170_d44x1y3_xerrplus, p8170_d44x1y3_yerrminus, p8170_d44x1y3_yerrplus);
  // p8170_d44x1y3.SetName("/HepData/8170/d44x1y3");
  // p8170_d44x1y3.SetTitle("/HepData/8170/d44x1y3");
  // lets plot the Rcp here

  // set up the error bars correctly: 
  Rcp_R2_Bayes[0] = (TH1F*)RAA_R2_Bayes[0]->Clone("Rcp_R2_Bayes_cent0_cent5");
  Rcp_R2_Bayes[0]->Divide(RAA_R2_Bayes[5]);
  Rcp_R2_Bayes[1] = (TH1F*)RAA_R2_Bayes[1]->Clone("Rcp_R2_Bayes_cent1_cent5");
  Rcp_R2_Bayes[1]->Divide(RAA_R2_Bayes[5]);

  Rcp_R3_Bayes[0] = (TH1F*)RAA_R3_Bayes[0]->Clone("Rcp_R3_Bayes_cent0_cent5");
  Rcp_R3_Bayes[0]->Divide(RAA_R3_Bayes[5]);
  Rcp_R3_Bayes[1] = (TH1F*)RAA_R3_Bayes[1]->Clone("Rcp_R3_Bayes_cent1_cent5");
  Rcp_R3_Bayes[1]->Divide(RAA_R3_Bayes[5]);

  Rcp_R4_Bayes[0] = (TH1F*)RAA_R4_Bayes[0]->Clone("Rcp_R4_Bayes_cent0_cent5");
  Rcp_R4_Bayes[0]->Divide(RAA_R4_Bayes[5]);
  Rcp_R4_Bayes[1] = (TH1F*)RAA_R4_Bayes[1]->Clone("Rcp_R4_Bayes_cent1_cent5");
  Rcp_R4_Bayes[1]->Divide(RAA_R4_Bayes[5]);

  Rcp3_vs_Rcp2[0] = (TH1F*)Rcp_R3_Bayes[0]->Clone("Rcp3_vs_Rcp2_cent0_cent5");
  Rcp3_vs_Rcp2[0]->Divide(Rcp_R2_Bayes[0]);
  Rcp3_vs_Rcp2[1] = (TH1F*)Rcp_R3_Bayes[1]->Clone("Rcp3_vs_Rcp2_cent1_cent5");
  Rcp3_vs_Rcp2[1]->Divide(Rcp_R2_Bayes[1]);  

  Rcp4_vs_Rcp2[0] = (TH1F*)Rcp_R4_Bayes[0]->Clone("Rcp4_vs_Rcp2_cent0_cent5");
  Rcp4_vs_Rcp2[0]->Divide(Rcp_R2_Bayes[0]);
  Rcp4_vs_Rcp2[1] = (TH1F*)Rcp_R4_Bayes[1]->Clone("Rcp4_vs_Rcp2_cent1_cent5");
  Rcp4_vs_Rcp2[1]->Divide(Rcp_R2_Bayes[1]);

  TCanvas * cRcp = new TCanvas("cRcp","",1000,800);
  makeMultiPanelCanvas(cRcp,3,1,0.0,0.0,0.2,0.15,0.07);
  TLegend * leg_rcp = myLegend(0.6,0.7,0.8,0.9);
  
  cRcp->cd(1);
  for(int i = 0; i<2; i++){

    makeHistTitle(Rcp_R2_Bayes[i]," ","  Jet p_{T} (GeV/c) ","R_{cp}");
    Rcp_R2_Bayes[i]->SetMarkerStyle(33+i);
    Rcp_R2_Bayes[i]->SetMarkerColor(2+i);
    Rcp_R2_Bayes[i]->SetAxisRange(50,299,"X");
    Rcp_R2_Bayes[i]->SetAxisRange(0,2,"Y");
    p8170_d2x1y1.SetMaximum(2);
    p8170_d2x1y1.SetMaximum(0);
    p8170_d2x1y1.GetXaxis()->SetLimits( 50, 299); 
    if(i == 0)p8170_d2x1y1.Draw("AP");
    if(i == 0)Rcp_R2_Bayes[i]->Draw("same p");
    if(i == 1)Rcp_R2_Bayes[i]->Draw("same p");
    if(i == 0)leg_rcp->AddEntry(Rcp_R2_Bayes[i],"#frac{0-5%}{70-90%}","pl");
    if(i == 1)leg_rcp->AddEntry(Rcp_R2_Bayes[i],"#frac{5-10%}{70-90%}","pl");
    
  }
  leg_rcp->SetTextSize(0.04);
  leg_rcp->Draw();
  drawText("R=0.2",0.25,0.25,14);
  cRcp->cd(2);
  for(int i = 0; i<2; i++){

    makeHistTitle(Rcp_R3_Bayes[i]," ","  Jet p_{T} (GeV/c)","R_{cp}");
    Rcp_R3_Bayes[i]->SetMarkerStyle(33+i);
    Rcp_R3_Bayes[i]->SetMarkerColor(2+i);
    Rcp_R3_Bayes[i]->SetAxisRange(50,299,"X");
    Rcp_R3_Bayes[i]->SetAxisRange(0,2,"Y");
    // p8170_d2x1y2.SetMaximum(2);
    // p8170_d2x1y2.SetMaximum(0);
    // p8170_d2x1y2.GetXaxis()->SetLimits( 50, 299); 
    if(i == 0)p8170_d2x1y2.Draw("AP");
    if(i == 0)Rcp_R3_Bayes[i]->Draw("same p");
    if(i == 1)Rcp_R3_Bayes[i]->Draw("same p");
    
  }
  drawText("R=0.3",0.2,0.25,14);
  cRcp->cd(3);
  for(int i = 0; i<2; i++){

    makeHistTitle(Rcp_R4_Bayes[i]," ","  Jet p_{T} (GeV/c)","R_{cp}");
    Rcp_R4_Bayes[i]->SetMarkerStyle(33+i);
    Rcp_R4_Bayes[i]->SetMarkerColor(2+i);
    Rcp_R4_Bayes[i]->SetAxisRange(50,299,"X");
    Rcp_R4_Bayes[i]->SetAxisRange(0,2,"Y");
    p8170_d2x1y3.SetMaximum(2);
    p8170_d2x1y3.SetMaximum(0);
    p8170_d2x1y3.GetXaxis()->SetLimits( 50, 299); 
    if(i == 0)p8170_d2x1y3.Draw("AP");
    if(i == 0)Rcp_R4_Bayes[i]->Draw("same p");
    if(i == 1)Rcp_R4_Bayes[i]->Draw("same p");
    
  }
  drawText("R=0.4",0.2,0.25,14);
  cRcp->cd(1);
  putCMSPrel();
  cRcp->SaveAs(Form("June4/Analysis_full_bins/Rcp_05_510_ratioto_5060_R234_%dGevCut_ak%s%s_%d_%s.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRcp->SaveAs(Form("June4/Analysis_full_bins/Rcp_05_510_ratioto_5060_R234_%dGevCut_ak%s%s_%d_%s.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRcp->SaveAs(Form("June4/Analysis_full_bins/Rcp_05_510_ratioto_5060_R234_%dGevCut_ak%s%s_%d_%s.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

  // now to make the Rcp double ratio: Rcp3 vs Rcp2
  TCanvas * cRcpvs2 = new TCanvas("cRcpvs2","",1000,800);

  TLegend * leg_Rcpvs2 = myLegend(0.6,0.7,0.8,0.9);
  
  for(int i = 0; i<2; i++){

    cRcpvs2->SetLogx();
    makeHistTitle(Rcp3_vs_Rcp2[i]," "," Jet p_{T} (GeV/c) ","#frac{R_{cp} R=0.X}{R_{cp} R=0.2}");
    Rcp3_vs_Rcp2[i]->SetMarkerStyle(33+i);
    Rcp3_vs_Rcp2[i]->SetMarkerColor(2+i);
    Rcp3_vs_Rcp2[i]->SetAxisRange(60,200,"X");
    Rcp3_vs_Rcp2[i]->SetAxisRange(0.6,2,"Y");
    p8170_d44x1y1.SetMaximum(2);
    p8170_d44x1y1.SetMaximum(0);
    p8170_d44x1y1.GetXaxis()->SetLimits( 50, 200); 
    if(i==0)p8170_d44x1y1.Draw("AP");
    if(i==0)Rcp3_vs_Rcp2[i]->Draw("same p");
    if(i==1)Rcp3_vs_Rcp2[i]->Draw("same p");

    Rcp4_vs_Rcp2[i]->SetMarkerStyle(24+i);
    Rcp4_vs_Rcp2[i]->SetMarkerColor(6+i);
    Rcp4_vs_Rcp2[i]->SetAxisRange(50,299,"X");
    if(i==0)p8170_d44x1y2.Draw("P same");
    Rcp4_vs_Rcp2[i]->Draw("same p");
    Rcp4_vs_Rcp2[i]->Draw("same p");
    
  }

  leg_Rcpvs2->AddEntry(Rcp3_vs_Rcp2[0],"X = 0.3, c = 0-5%","pl");
  leg_Rcpvs2->AddEntry(Rcp3_vs_Rcp2[1],"X = 0.3, c = 5-10%","pl");
  leg_Rcpvs2->AddEntry(Rcp4_vs_Rcp2[0],"X = 0.4, c = 0-5%","pl");
  leg_Rcpvs2->AddEntry(Rcp4_vs_Rcp2[1],"X = 0.4, c = 5-10%","pl");
  leg_Rcpvs2->SetTextSize(0.03);
  leg_Rcpvs2->Draw();

  putCMSPrel();
  cRcpvs2->SaveAs(Form("June4/Analysis_full_bins/Rcp34_vs_Rcp2_05_510_ratioto_5060_%dGevCut_ak%s%s_%d_%s.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRcpvs2->SaveAs(Form("June4/Analysis_full_bins/Rcp34_vs_Rcp2_05_510_ratioto_5060_%dGevCut_ak%s%s_%d_%s.C",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  cRcpvs2->SaveAs(Form("June4/Analysis_full_bins/Rcp34_vs_Rcp2_05_510_ratioto_5060_%dGevCut_ak%s%s_%d_%s.root",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");
  
  // im  plotting the raa vs npart plot here, since it needs the systematics
  // plot 3 - RAA as a function of npart - taken from http://dde.web.cern.ch/dde/glauber_lhc.htm for 84 < pT < 97 in PbPb,PP
  //       - need to decide if we have to unfold this? or if we can just take that respective pt ranges from the already existing RAA histograms.  this is bin number 16 from the centrality classes weve measured.

  // // get the statistical Error bars from the unfolding data driven check
  // TFile * fR2_err, * fR3_err, * fR4_err;
  // fR2_err = TFile::Open(Form("PbPb_R2_pp_R2_n20_eta_p20_unfoldingCut_%d_data_driven_correction_akPuPF_20150417.root",unfoldingCut));
  // fR3_err = TFile::Open(Form("PbPb_R3_pp_R3_n20_eta_p20_unfoldingCut_%d_data_driven_correction_akPuPF_20150417.root",unfoldingCut));
  // fR4_err = TFile::Open(Form("PbPb_R4_pp_R4_n20_eta_p20_unfoldingCut_%d_data_driven_correction_akPuPF_20150417.root",unfoldingCut));

  // TH1F * R2_error[nbins_cent], * R3_error[nbins_cent], * R4_error[nbins_cent];

  // get the responsible histograms for this.
  TH1F * hRAA_R2_npart = new TH1F("hRAA_R2_npart","",450, 0, 450);
  //hRAA_R2_npart->LabelsOption(">","X");
  TH1F * hRAA_R3_npart = new TH1F("hRAA_R3_npart","",450, 0, 450);
  //hRAA_R3_npart->LabelsOption(">","X");
  TH1F * hRAA_R4_npart = new TH1F("hRAA_R4_npart","",450, 0, 450);
  //hRAA_R4_npart->LabelsOption(">","X");

  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart->SetBinContent(hRAA_R2_npart->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinContent(RAA_R2_Bayes[i]->FindBin(100)));
    hRAA_R2_npart->SetBinError(hRAA_R2_npart->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinError(RAA_R2_Bayes[i]->FindBin(100)));
    
    hRAA_R3_npart->SetBinContent(hRAA_R3_npart->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(100)));
    hRAA_R3_npart->SetBinError(hRAA_R3_npart->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinError(RAA_R3_Bayes[i]->FindBin(100)));
    
    hRAA_R4_npart->SetBinContent(hRAA_R4_npart->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinContent(RAA_R4_Bayes[i]->FindBin(100)));    
    hRAA_R4_npart->SetBinError(hRAA_R4_npart->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinError(RAA_R4_Bayes[i]->FindBin(100)));    
  }

  TLine *lineRAA_R2_npart = new TLine(0,1,450,1);
  lineRAA_R2_npart->SetLineStyle(2);
  lineRAA_R2_npart->SetLineWidth(2);

  TCanvas * cRAA_npart = new TCanvas("cRAA_npart","",800,600);
  //cRAA_npart->SetGridy();
  //cRAA_npart->SetGridx();

  hRAA_R2_npart->SetTitle(" ");
  hRAA_R2_npart->SetXTitle(" N_{part} ");
  hRAA_R2_npart->SetYTitle(" R_{AA} ");
  hRAA_R2_npart->SetMarkerColor(kRed);
  hRAA_R2_npart->SetLineColor(kRed);
  hRAA_R2_npart->SetMarkerStyle(20);
  hRAA_R2_npart->SetAxisRange(0,1.8, "Y");
  hRAA_R2_npart->Draw("E1");
  hRAA_R3_npart->SetMarkerColor(kBlack);
  hRAA_R3_npart->SetLineColor(kBlack);
  hRAA_R3_npart->SetMarkerStyle(20);
  hRAA_R3_npart->Draw("same E1");
  hRAA_R4_npart->SetMarkerColor(kBlue);
  hRAA_R4_npart->SetLineColor(kBlue);
  hRAA_R4_npart->SetMarkerStyle(20);
  hRAA_R4_npart->Draw("same E1");

  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.5,0.25,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.25,0.2,16);
  drawText("97 < Jet p_{T} < 114",0.25,0.25,16);
  
  TLegend * npart1 = myLegend(0.6,0.7,0.8,0.9);
  npart1->AddEntry(hRAA_R2_npart,"R=0.2", "pl");
  npart1->AddEntry(hRAA_R3_npart,"R=0.3", "pl");
  npart1->AddEntry(hRAA_R4_npart,"R=0.4", "pl");
  npart1->SetTextSize(0.04);
  npart1->Draw();
 
  for(int i = 0; i<nbins_cent; ++i){
    //systematics_R3.calcTotalSys(i);
    systematics_R3.DrawNpartSys(RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(100)),i,npart[i], 100);
  }

  hRAA_R2_npart->Draw("same E1");
  hRAA_R3_npart->Draw("same E1");
  hRAA_R4_npart->Draw("same E1");
  lineRAA_R2_npart->Draw();
  putCMSPrel(0.2,0.85,0.025);
  putPbPbLumi(0.2,0.8,0.025);
  putPPLumi(0.2,0.75,0.025);
  
  cRAA_npart->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_97_pT_114_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_97_pT_114_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_97_pT_114_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");

  // plot npart raa for dufferent pt range around 130
  
  // get the responsible histograms for this.
  TH1F * hRAA_R2_npart2 = new TH1F("hRAA_R2_npart2","",450, 0, 450);
  //hRAA_R2_npart2->LabelsOption(">","X");
  TH1F * hRAA_R3_npart2 = new TH1F("hRAA_R3_npart2","",450, 0, 450);
  //hRAA_R3_npart2->LabelsOption(">","X");
  TH1F * hRAA_R4_npart2 = new TH1F("hRAA_R4_npart2","",450, 0, 450);
  //hRAA_R4_npart2->LabelsOption(">","X");

  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart2->SetBinContent(hRAA_R2_npart2->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinContent(RAA_R2_Bayes[i]->FindBin(135)));
    hRAA_R2_npart2->SetBinError(hRAA_R2_npart2->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinError(RAA_R2_Bayes[i]->FindBin(135)));
    
    hRAA_R3_npart2->SetBinContent(hRAA_R3_npart2->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(135)));
    hRAA_R3_npart2->SetBinError(hRAA_R3_npart2->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinError(RAA_R3_Bayes[i]->FindBin(135)));
    
    hRAA_R4_npart2->SetBinContent(hRAA_R4_npart2->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinContent(RAA_R4_Bayes[i]->FindBin(135)));    
    hRAA_R4_npart2->SetBinError(hRAA_R4_npart2->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinError(RAA_R4_Bayes[i]->FindBin(135)));    
  }


  TCanvas * cRAA_npart2 = new TCanvas("cRAA_npart2","",800,600);
  //cRAA_npart2->SetGridy();
  //cRAA_npart2->SetGridx();

  hRAA_R2_npart2->SetTitle(" ");
  hRAA_R2_npart2->SetXTitle(" N_{part} ");
  hRAA_R2_npart2->SetYTitle(" R_{AA} ");
  hRAA_R2_npart2->SetMarkerColor(kRed);
  hRAA_R2_npart2->SetLineColor(kRed);
  hRAA_R2_npart2->SetMarkerStyle(20);
  hRAA_R2_npart2->SetAxisRange(0,1.8, "Y");
  hRAA_R2_npart2->Draw("E1");
  hRAA_R3_npart2->SetMarkerColor(kBlack);
  hRAA_R3_npart2->SetLineColor(kBlack);
  hRAA_R3_npart2->SetMarkerStyle(20);
  hRAA_R3_npart2->Draw("same E1");
  hRAA_R4_npart2->SetMarkerColor(kBlue);
  hRAA_R4_npart2->SetLineColor(kBlue);
  hRAA_R4_npart2->SetMarkerStyle(20);
  hRAA_R4_npart2->Draw("same E1");
  lineRAA_R2_npart->Draw();

  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f",etaBoundary),0.5,0.25,16);
  //drawText("|vz|<15, HBHEfilter, pCES",0.25,0.2,16);
  drawText("133 < Jet p_{T} < 153", 0.25,0.25,16);
  
  TLegend * npart2 = myLegend(0.6,0.7,0.8,0.9);
  npart2->AddEntry(hRAA_R2_npart2,"R=0.2", "pl");
  npart2->AddEntry(hRAA_R3_npart2,"R=0.3", "pl");
  npart2->AddEntry(hRAA_R4_npart2,"R=0.4", "pl");
  npart2->SetTextSize(0.04);
  npart2->Draw();
 
  for(int i = 0; i<nbins_cent; ++i){
    //systematics_R3.calcTotalSys(i);
    systematics_R3.DrawNpartSys(RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(135)),i,npart[i], 135);
  }

  hRAA_R2_npart2->Draw("same E1");
  hRAA_R3_npart2->Draw("same E1");
  hRAA_R4_npart2->Draw("same E1");
  putCMSPrel(0.2,0.85,0.025);
  putPbPbLumi(0.2,0.8,0.025);
  putPPLumi(0.2,0.75,0.025);

  cRAA_npart2->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_133_pT_153_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart2->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_133_pT_153_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart2->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_133_pT_153_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");
  
  // plot RAA vs Npart for a smaller pT : 74 - 84
  
  // get the responsible histograms for this.
  TH1F * hRAA_R2_npart3 = new TH1F("hRAA_R2_npart3","",450, 0, 450);
  //hRAA_R2_npart3->LabelsOption(">","X");
  TH1F * hRAA_R3_npart3 = new TH1F("hRAA_R3_npart3","",450, 0, 450);
  //hRAA_R3_npart3->LabelsOption(">","X");
  TH1F * hRAA_R4_npart3 = new TH1F("hRAA_R4_npart3","",450, 0, 450);
  //hRAA_R4_npart3->LabelsOption(">","X");
  
  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart3->SetBinContent(hRAA_R2_npart3->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinContent(RAA_R2_Bayes[i]->FindBin(80)));
    hRAA_R2_npart3->SetBinError(hRAA_R2_npart3->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinError(RAA_R2_Bayes[i]->FindBin(80)));
    
    hRAA_R3_npart3->SetBinContent(hRAA_R3_npart3->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(80)));
    hRAA_R3_npart3->SetBinError(hRAA_R3_npart3->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinError(RAA_R3_Bayes[i]->FindBin(80)));
    
    hRAA_R4_npart3->SetBinContent(hRAA_R4_npart3->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinContent(RAA_R4_Bayes[i]->FindBin(80)));    
    hRAA_R4_npart3->SetBinError(hRAA_R4_npart3->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinError(RAA_R4_Bayes[i]->FindBin(80)));    
  }


  TCanvas * cRAA_npart3 = new TCanvas("cRAA_npart3","",800,600);
  //cRAA_npart3->SetGridy();
  //cRAA_npart3->SetGridx();

  hRAA_R2_npart3->SetTitle(" ");
  hRAA_R2_npart3->SetXTitle(" N_{part} ");
  hRAA_R2_npart3->SetYTitle(" R_{AA} ");
  hRAA_R2_npart3->SetMarkerColor(kRed);
  hRAA_R2_npart3->SetLineColor(kRed);
  hRAA_R2_npart3->SetMarkerStyle(20);
  hRAA_R2_npart3->SetAxisRange(0,1.8, "Y");
  hRAA_R2_npart3->Draw("E1");
  hRAA_R3_npart3->SetMarkerColor(kBlack);
  hRAA_R3_npart3->SetLineColor(kBlack);
  hRAA_R3_npart3->SetMarkerStyle(20);
  hRAA_R3_npart3->Draw("same E1");
  hRAA_R4_npart3->SetMarkerColor(kBlue);
  hRAA_R4_npart3->SetLineColor(kBlue);
  hRAA_R4_npart3->SetMarkerStyle(20);
  hRAA_R4_npart3->Draw("same E1");
  lineRAA_R2_npart->Draw();

  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.5,0.25,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.25,0.2,16);
  drawText("74 < Jet p_{T} < 84", 0.25,0.25,16);
  
  TLegend * npart3 = myLegend(0.6,0.7,0.8,0.9);
  npart3->AddEntry(hRAA_R2_npart3,"R=0.2", "pl");
  npart3->AddEntry(hRAA_R3_npart3,"R=0.3", "pl");
  npart3->AddEntry(hRAA_R4_npart3,"R=0.4", "pl");
  npart3->SetTextSize(0.04);
  npart3->Draw();
 
  for(int i = 0; i<nbins_cent; ++i){
    //systematics_R3.calcTotalSys(i);
    systematics_R3.DrawNpartSys(RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(80)),i,npart[i], 80);
  }

  hRAA_R2_npart3->Draw("same E1");
  hRAA_R3_npart3->Draw("same E1");
  hRAA_R4_npart3->Draw("same E1");

  putCMSPrel(0.2,0.85,0.025);
  putPbPbLumi(0.2,0.8,0.025);
  putPPLumi(0.2,0.75,0.025);

  cRAA_npart3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_74_pT_84_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_74_pT_84_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_74_pT_84_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");
  
  // draw the RAA vs npart 
  // plot RAA vs Npart for a smaller pT : 64 - 74
  
  // get the responsible histograms for this.
  TH1F * hRAA_R2_npart4 = new TH1F("hRAA_R2_npart4","",450, 0, 450);
  //hRAA_R2_npart4->LabelsOption(">","X");
  TH1F * hRAA_R3_npart4 = new TH1F("hRAA_R3_npart4","",450, 0, 450);
  //hRAA_R3_npart4->LabelsOption(">","X");
  TH1F * hRAA_R4_npart4 = new TH1F("hRAA_R4_npart4","",450, 0, 450);
  //hRAA_R4_npart4->LabelsOption(">","X");
  
  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart4->SetBinContent(hRAA_R2_npart4->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinContent(RAA_R2_Bayes[i]->FindBin(70)));
    hRAA_R2_npart4->SetBinError(hRAA_R2_npart4->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinError(RAA_R2_Bayes[i]->FindBin(70)));
    
    hRAA_R3_npart4->SetBinContent(hRAA_R3_npart4->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(70)));
    hRAA_R3_npart4->SetBinError(hRAA_R3_npart4->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinError(RAA_R3_Bayes[i]->FindBin(70)));
    
    hRAA_R4_npart4->SetBinContent(hRAA_R4_npart4->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinContent(RAA_R4_Bayes[i]->FindBin(70)));    
    hRAA_R4_npart4->SetBinError(hRAA_R4_npart4->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinError(RAA_R4_Bayes[i]->FindBin(70)));    
  }


  TCanvas * cRAA_npart4 = new TCanvas("cRAA_npart4","",800,600);
  //cRAA_npart4->SetGridy();
  //cRAA_npart4->SetGridx();

  hRAA_R2_npart4->SetTitle(" ");
  hRAA_R2_npart4->SetXTitle(" N_{part} ");
  hRAA_R2_npart4->SetYTitle(" R_{AA} ");
  hRAA_R2_npart4->SetMarkerColor(kRed);
  hRAA_R2_npart4->SetLineColor(kRed);
  hRAA_R2_npart4->SetMarkerStyle(20);
  hRAA_R2_npart4->SetAxisRange(0,1.8, "Y");
  hRAA_R2_npart4->Draw("E1");
  hRAA_R3_npart4->SetMarkerColor(kBlack);
  hRAA_R3_npart4->SetLineColor(kBlack);
  hRAA_R3_npart4->SetMarkerStyle(20);
  hRAA_R3_npart4->Draw("same E1");
  hRAA_R4_npart4->SetMarkerColor(kBlue);
  hRAA_R4_npart4->SetLineColor(kBlue);
  hRAA_R4_npart4->SetMarkerStyle(20);
  hRAA_R4_npart4->Draw("same E1");
  lineRAA_R2_npart->Draw();

  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.5,0.25,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.25,0.2,16);
  drawText("64 < Jet p_{T} < 74", 0.25,0.25,16);
  
  TLegend * npart4 = myLegend(0.6,0.7,0.8,0.9);
  npart4->AddEntry(hRAA_R2_npart4,"R=0.2", "pl");
  npart4->AddEntry(hRAA_R3_npart4,"R=0.3", "pl");
  npart4->AddEntry(hRAA_R4_npart4,"R=0.4", "pl");
  npart4->SetTextSize(0.04);
  npart4->Draw();
 
  for(int i = 0; i<nbins_cent; ++i){
    //systematics_R3.calcTotalSys(i);
    systematics_R3.DrawNpartSys(RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(70)),i,npart[i], 70);
  }

  hRAA_R2_npart4->Draw("same E1");
  hRAA_R3_npart4->Draw("same E1");
  hRAA_R4_npart4->Draw("same E1");

  putCMSPrel(0.2,0.85,0.025);
  putPbPbLumi(0.2,0.8,0.025);
  putPPLumi(0.2,0.75,0.025);

  cRAA_npart4->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_64_pT_74_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart4->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_64_pT_74_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart4->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_64_pT_74_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");

  // plot RAA vs Npart for a smaller pT : 84 - 97 
  
  // get the responsible histograms for this.
  TH1F * hRAA_R2_npart5 = new TH1F("hRAA_R2_npart5","",450, 0, 450);
  //hRAA_R2_npart5->LabelsOption(">","X");
  TH1F * hRAA_R3_npart5 = new TH1F("hRAA_R3_npart5","",450, 0, 450);
  //hRAA_R3_npart5->LabelsOption(">","X");
  TH1F * hRAA_R4_npart5 = new TH1F("hRAA_R4_npart5","",450, 0, 450);
  //hRAA_R4_npart5->LabelsOption(">","X");
  
  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart5->SetBinContent(hRAA_R2_npart5->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinContent(RAA_R2_Bayes[i]->FindBin(90)));
    hRAA_R2_npart5->SetBinError(hRAA_R2_npart5->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinError(RAA_R2_Bayes[i]->FindBin(90)));
    
    hRAA_R3_npart5->SetBinContent(hRAA_R3_npart5->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(90)));
    hRAA_R3_npart5->SetBinError(hRAA_R3_npart5->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinError(RAA_R3_Bayes[i]->FindBin(90)));
    
    hRAA_R4_npart5->SetBinContent(hRAA_R4_npart5->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinContent(RAA_R4_Bayes[i]->FindBin(90)));    
    hRAA_R4_npart5->SetBinError(hRAA_R4_npart5->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinError(RAA_R4_Bayes[i]->FindBin(90)));    
  }


  TCanvas * cRAA_npart5 = new TCanvas("cRAA_npart5","",800,600);
  //cRAA_npart5->SetGridy();
  //cRAA_npart5->SetGridx();

  hRAA_R2_npart5->SetTitle(" ");
  hRAA_R2_npart5->SetXTitle(" N_{part} ");
  hRAA_R2_npart5->SetYTitle(" R_{AA} ");
  hRAA_R2_npart5->SetMarkerColor(kRed);
  hRAA_R2_npart5->SetLineColor(kRed);
  hRAA_R2_npart5->SetMarkerStyle(20);
  hRAA_R2_npart5->SetAxisRange(0,1.8, "Y");
  hRAA_R2_npart5->Draw("E1");
  hRAA_R3_npart5->SetMarkerColor(kBlack);
  hRAA_R3_npart5->SetLineColor(kBlack);
  hRAA_R3_npart5->SetMarkerStyle(20);
  hRAA_R3_npart5->Draw("same E1");
  hRAA_R4_npart5->SetMarkerColor(kBlue);
  hRAA_R4_npart5->SetLineColor(kBlue);
  hRAA_R4_npart5->SetMarkerStyle(20);
  hRAA_R4_npart5->Draw("same E1");
  lineRAA_R2_npart->Draw();

  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.5,0.25,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.25,0.2,16);
  drawText("84 < Jet p_{T} < 97", 0.25,0.25,16);
  
  TLegend * npart5 = myLegend(0.6,0.7,0.8,0.9);
  npart5->AddEntry(hRAA_R2_npart5,"R=0.2", "pl");
  npart5->AddEntry(hRAA_R3_npart5,"R=0.3", "pl");
  npart5->AddEntry(hRAA_R4_npart5,"R=0.4", "pl");
  npart5->SetTextSize(0.04);
  npart5->Draw();
 
  for(int i = 0; i<nbins_cent; ++i){
    //systematics_R3.calcTotalSys(i);
    systematics_R3.DrawNpartSys(RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(90)),i,npart[i], 90);
  }

  hRAA_R2_npart5->Draw("same E1");
  hRAA_R3_npart5->Draw("same E1");
  hRAA_R4_npart5->Draw("same E1");

  putCMSPrel(0.2,0.85,0.025);
  putPbPbLumi(0.2,0.8,0.025);
  putPPLumi(0.2,0.75,0.025);

  cRAA_npart5->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_84_pT_97_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart5->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_84_pT_97_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart5->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_84_pT_97_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");

  // plot RAA vs Npart for a smaller pT : 84 - 97 
  
  // get the responsible histograms for this.
  TH1F * hRAA_R2_npart6 = new TH1F("hRAA_R2_npart6","",450, 0, 450);
  //hRAA_R2_npart6->LabelsOption(">","X");
  TH1F * hRAA_R3_npart6 = new TH1F("hRAA_R3_npart6","",450, 0, 450);
  //hRAA_R3_npart6->LabelsOption(">","X");
  TH1F * hRAA_R4_npart6 = new TH1F("hRAA_R4_npart6","",450, 0, 450);
  //hRAA_R4_npart6->LabelsOption(">","X");
  
  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart6->SetBinContent(hRAA_R2_npart6->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinContent(RAA_R2_Bayes[i]->FindBin(120)));
    hRAA_R2_npart6->SetBinError(hRAA_R2_npart6->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinError(RAA_R2_Bayes[i]->FindBin(120)));
    
    hRAA_R3_npart6->SetBinContent(hRAA_R3_npart6->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(120)));
    hRAA_R3_npart6->SetBinError(hRAA_R3_npart6->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinError(RAA_R3_Bayes[i]->FindBin(120)));
    
    hRAA_R4_npart6->SetBinContent(hRAA_R4_npart6->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinContent(RAA_R4_Bayes[i]->FindBin(120)));    
    hRAA_R4_npart6->SetBinError(hRAA_R4_npart6->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinError(RAA_R4_Bayes[i]->FindBin(120)));    
  }


  TCanvas * cRAA_npart6 = new TCanvas("cRAA_npart6","",800,600);
  //cRAA_npart6->SetGridy();
  //cRAA_npart6->SetGridx();

  hRAA_R2_npart6->SetTitle(" ");
  hRAA_R2_npart6->SetXTitle(" N_{part} ");
  hRAA_R2_npart6->SetYTitle(" R_{AA} ");
  hRAA_R2_npart6->SetMarkerColor(kRed);
  hRAA_R2_npart6->SetLineColor(kRed);
  hRAA_R2_npart6->SetMarkerStyle(20);
  hRAA_R2_npart6->SetAxisRange(0,1.8, "Y");
  hRAA_R2_npart6->Draw("E1");
  hRAA_R3_npart6->SetMarkerColor(kBlack);
  hRAA_R3_npart6->SetLineColor(kBlack);
  hRAA_R3_npart6->SetMarkerStyle(20);
  hRAA_R3_npart6->Draw("same E1");
  hRAA_R4_npart6->SetMarkerColor(kBlue);
  hRAA_R4_npart6->SetLineColor(kBlue);
  hRAA_R4_npart6->SetMarkerStyle(20);
  hRAA_R4_npart6->Draw("same E1");
  lineRAA_R2_npart->Draw();

  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.5,0.25,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.25,0.2,16);
  drawText("114 < Jet p_{T} < 133", 0.25,0.25,16);
  
  TLegend * npart6 = myLegend(0.6,0.7,0.8,0.9);
  npart6->AddEntry(hRAA_R2_npart6,"R=0.2", "pl");
  npart6->AddEntry(hRAA_R3_npart6,"R=0.3", "pl");
  npart6->AddEntry(hRAA_R4_npart6,"R=0.4", "pl");
  npart6->SetTextSize(0.04);
  npart6->Draw();
 
  for(int i = 0; i<nbins_cent; ++i){
    //systematics_R3.calcTotalSys(i);
    systematics_R3.DrawNpartSys(RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(120)),i,npart[i], 120);
  }

  hRAA_R2_npart6->Draw("same E1");
  hRAA_R3_npart6->Draw("same E1");
  hRAA_R4_npart6->Draw("same E1");

  putCMSPrel(0.2,0.85,0.025);
  putPbPbLumi(0.2,0.8,0.025);
  putPPLumi(0.2,0.75,0.025);

  cRAA_npart6->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_114_pT_133_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart6->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_114_pT_133_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_npart6->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_npart_114_pT_133_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");




  
  
  // plot the RAA for each radii here along with the systematics shown in yellow, for measured, bayesian and bin by bin unfoldeding.
  TCanvas *cRAA_R2 = new TCanvas("cRAA_R2","RAA",1400,1200);
  makeMultiPanelCanvas(cRAA_R2,3,2,0.0,0.0,0.2,0.15,0.07);
  
  TLegend *tRAA_R2 = myLegend(0.35,0.55,0.65,0.7);
  TLine *lineRAA_R2 = new TLine(unfoldingCut+30,1,299,1);
  lineRAA_R2->SetLineStyle(2);
  lineRAA_R2->SetLineWidth(2);
  
  for(int i = 0;i<nbins_cent;++i){

    cRAA_R2->cd(nbins_cent-i);

    // if(i==0){
    //   for(int j = RAA_R2_Bayes[i]->FindBin(unfoldingCut); j<RAA_R2_Bayes[i]->FindBin(100); ++j){
    // 	RAA_R2_Meas[i]->SetBinContent(j,0);
    // 	RAA_R2_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R2_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R2_Meas[i]->SetBinError(j,0);
    // 	RAA_R2_Bayes[i]->SetBinError(j,0);
    // 	RAA_R2_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==5 || i==4 || i==3){
    //   for(int j = RAA_R2_Bayes[i]->FindBin(240); j<RAA_R2_Bayes[i]->FindBin(300); ++j){
    // 	RAA_R2_Meas[i]->SetBinContent(j,0);
    // 	RAA_R2_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R2_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R2_Meas[i]->SetBinError(j,0);
    // 	RAA_R2_Bayes[i]->SetBinError(j,0);
    // 	RAA_R2_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==1 || i==2){
    //   for(int j = RAA_R2_Bayes[i]->FindBin(unfoldingCut); j<RAA_R2_Bayes[i]->FindBin(unfoldingCut+30); ++j){
    // 	RAA_R2_Meas[i]->SetBinContent(j,0);
    // 	RAA_R2_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R2_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R2_Meas[i]->SetBinError(j,0);
    // 	RAA_R2_Bayes[i]->SetBinError(j,0);
    // 	RAA_R2_BinByBin[i]->SetBinError(j,0);
    //   }
    // }

    RAA_R2_Meas[i]->SetMarkerColor(kBlack);
    RAA_R2_Meas[i]->SetMarkerStyle(24);
    makeHistTitle(RAA_R2_Meas[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    RAA_R2_Meas[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    RAA_R2_Meas[i]->SetAxisRange(0,2,"Y");
    RAA_R2_Meas[i]->Draw("E1");

    RAA_R2_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R2_Bayes[i]->SetMarkerStyle(20);
    RAA_R2_Bayes[i]->Draw("same E1");

    RAA_R2_BinByBin[i]->SetMarkerStyle(33);
    RAA_R2_BinByBin[i]->SetMarkerColor(kRed);
    RAA_R2_BinByBin[i]->Draw("same E1");

    lineRAA_R2->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.85,20);

    systematics_R2.calcTotalSys(i);
    systematics_R2.Draw(RAA_R2_Bayes[i],i,2);
    
    RAA_R2_Meas[i]->Draw("same E1");
    RAA_R2_BinByBin[i]->Draw("same E1");
    RAA_R2_Bayes[i]->Draw("same E1");

  }
    
  tRAA_R2->AddEntry(RAA_R2_Meas[0],"Measured","pl");
  tRAA_R2->AddEntry(RAA_R2_Bayes[0],"Bayesian","pl");
  tRAA_R2->AddEntry(RAA_R2_BinByBin[0],"bin-by-bin","pl");
  tRAA_R2->SetTextSize(0.04);

  cRAA_R2->cd(1);
  tRAA_R2->Draw();
  cRAA_R2->cd(1);
  putCMSPrel();
  cRAA_R2->cd(2);
  putPbPbLumi();
  cRAA_R2->cd(3);
  putPPLumi();
  cRAA_R2->cd(1);
  drawText(Form("Anti-k_{T} %s R=0.2 %s Jets",algo,jet_type),0.25,0.20,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_R2->cd(2);
  drawText(Form("Jet ID cut, |#eta|<%2.0f",etaBoundary),0.1,0.2,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.1,16);
  cRAA_R2->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.2,16);
  drawText("Pile up rejection cut applied",0.1,0.1,16);
  cRAA_R2->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R2_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_R2->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R2_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_R2->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R2_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");


  // draw it for R=0.3
  TCanvas *cRAA_onlyB_R3 = new TCanvas("cRAA_onlyB_R3","RAA",1400,1200);
  makeMultiPanelCanvas(cRAA_onlyB_R3,3,2,0.0,0.0,0.2,0.15,0.07);

  TLegend *tRAA_onlyB_R3 = myLegend(0.35,0.55,0.65,0.7);
  TLine *lineRAA_onlyB_R3 = new TLine(unfoldingCut+30,1,299,1);
  lineRAA_onlyB_R3->SetLineStyle(2);
  lineRAA_onlyB_R3->SetLineWidth(2);
  
  for(int i = 0;i<nbins_cent;++i){

    cRAA_onlyB_R3->cd(nbins_cent-i);

    // if(i==0){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(100); ++j){
    // 	RAA_R3_Meas[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R3_Meas[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==5 || i==4 || i==3){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(240); j<RAA_R3_Bayes[i]->FindBin(300); ++j){
    // 	RAA_R3_Meas[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R3_Meas[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==1 || i==2){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(unfoldingCut+30); ++j){
    // 	RAA_R3_Meas[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R3_Meas[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinError(j,0);
    //   }
    // }

    RAA_R3_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R3_Bayes[i]->SetMarkerStyle(24);
    makeHistTitle(RAA_R3_Bayes[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    RAA_R3_Bayes[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    RAA_R3_Bayes[i]->SetAxisRange(0,2,"Y");
    RAA_R3_Bayes[i]->Draw("E1");

    // RAA_R3_Bayes[i]->SetMarkerColor(kBlack);
    // RAA_R3_Bayes[i]->SetMarkerStyle(20);
    // RAA_R3_Bayes[i]->Draw("same E1");

    // RAA_R3_BinByBin[i]->SetMarkerStyle(33);
    // RAA_R3_BinByBin[i]->SetMarkerColor(kRed);
    // RAA_R3_BinByBin[i]->Draw("same E1");

    lineRAA_onlyB_R3->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.85,20);

    systematics_R3.calcTotalSys(i);
    systematics_R3.Draw(RAA_R3_Bayes[i],i,2);
    
    //RAA_R3_Meas[i]->Draw("same E1");
    //RAA_R3_BinByBin[i]->Draw("same E1");
    RAA_R3_Bayes[i]->Draw("same E1");

  }
    
  //tRAA_onlyB_R3->AddEntry(RAA_R3_Meas[0],"Measured","pl");
  tRAA_onlyB_R3->AddEntry(RAA_R3_Bayes[0],"Bayesian","pl");
  //tRAA_R3->AddEntry(RAA_R3_BinByBin[0],"bin-by-bin","pl");
  tRAA_onlyB_R3->SetTextSize(0.04);

  cRAA_onlyB_R3->cd(1);
  tRAA_onlyB_R3->Draw();
  cRAA_onlyB_R3->cd(1);
  putCMSPrel();
  cRAA_onlyB_R3->cd(2);
  putPbPbLumi();
  cRAA_onlyB_R3->cd(3);
  putPPLumi();
  cRAA_onlyB_R3->cd(1);
  drawText(Form("Anti-k_{T} %s R=0.3 %s Jets",algo,jet_type),0.25,0.20,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_onlyB_R3->cd(2);
  drawText(Form("|#eta|<%2.0f",etaBoundary),0.1,0.2,16);
  //drawText("|vz|<15, HBHEfilter, pCES",0.1,0.1,16);
  //cRAA_onlyB_R3->cd(3);
  //drawText("Jet RAA dataset, trigger combined",0.1,0.2,16);
  //drawText("Pile up rejection cut applied",0.1,0.1,16);

  cRAA_onlyB_R3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R3_Bayes_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_onlyB_R3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R3_Bayes_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_onlyB_R3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R3_Bayes_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");


  ///
  // Plot the comparison with the PAS here with the systematic boxes 
  
  TH1F *PASRAA_bayesian[3][nbins_cent];

  Double_t xAxis2[12] = {100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300}; 

  PASRAA_bayesian[1][5] = new TH1F("PASRAA_bayesian_R3_cent5","",11,xAxis2);
  PASRAA_bayesian[1][5]->SetBinContent(0,0.7225377);
  PASRAA_bayesian[1][5]->SetBinContent(1,0.8600517);
  PASRAA_bayesian[1][5]->SetBinContent(2,0.8100158);
  PASRAA_bayesian[1][5]->SetBinContent(3,0.749716);
  PASRAA_bayesian[1][5]->SetBinContent(4,0.728152);
  PASRAA_bayesian[1][5]->SetBinContent(5,0.8079244);
  PASRAA_bayesian[1][5]->SetBinContent(6,1.045846);
  PASRAA_bayesian[1][5]->SetBinContent(7,0.99492);
  PASRAA_bayesian[1][5]->SetBinContent(8,0.8501027);
  PASRAA_bayesian[1][5]->SetBinContent(9,0.7283186);
  PASRAA_bayesian[1][5]->SetBinContent(10,0.8501291);
  PASRAA_bayesian[1][5]->SetBinContent(11,0.8562201);
  PASRAA_bayesian[1][5]->SetBinContent(12,2.485746);

  PASRAA_bayesian[1][5]->SetBinError(0,0.006807087);
  PASRAA_bayesian[1][5]->SetBinError(1,0.03660645);
  PASRAA_bayesian[1][5]->SetBinError(2,0.04424438);
  PASRAA_bayesian[1][5]->SetBinError(3,0.05483107);
  PASRAA_bayesian[1][5]->SetBinError(4,0.06882145);
  PASRAA_bayesian[1][5]->SetBinError(5,0.100899);
  PASRAA_bayesian[1][5]->SetBinError(6,0.1601293);
  PASRAA_bayesian[1][5]->SetBinError(7,0.181309);
  PASRAA_bayesian[1][5]->SetBinError(8,0.1830885);
  PASRAA_bayesian[1][5]->SetBinError(9,0.1512963);
  PASRAA_bayesian[1][5]->SetBinError(10,0.1969083);
  PASRAA_bayesian[1][5]->SetBinError(11,0.3481827);
  PASRAA_bayesian[1][5]->SetBinError(12,0.7989296);
   
  // PASRAA_bayesian[1][5]->SetBinError(0,0.006807087);
  // PASRAA_bayesian[1][5]->SetBinError(1,0.01931266);
  // PASRAA_bayesian[1][5]->SetBinError(2,0.02226945);
  // PASRAA_bayesian[1][5]->SetBinError(3,0.02646403);
  // PASRAA_bayesian[1][5]->SetBinError(4,0.03199702);
  // PASRAA_bayesian[1][5]->SetBinError(5,0.04537526);
  // PASRAA_bayesian[1][5]->SetBinError(6,0.06991723);
  // PASRAA_bayesian[1][5]->SetBinError(7,0.07712956);
  // PASRAA_bayesian[1][5]->SetBinError(8,0.07612922);
  // PASRAA_bayesian[1][5]->SetBinError(9,0.06113506);
  // PASRAA_bayesian[1][5]->SetBinError(10,0.07653623);
  // PASRAA_bayesian[1][5]->SetBinError(11,0.1332746);
  // PASRAA_bayesian[1][5]->SetBinError(12,0.7989296);
  PASRAA_bayesian[1][5]->SetMinimum(0);
  PASRAA_bayesian[1][5]->SetMaximum(2);
  PASRAA_bayesian[1][5]->SetEntries(1770.971);
  
  PASRAA_bayesian[1][4] = new TH1F("PASRAA_bayesian_R3_cent4","",11,xAxis2);
  PASRAA_bayesian[1][4]->SetBinContent(0,0.7299647);
  PASRAA_bayesian[1][4]->SetBinContent(1,0.7481229);
  PASRAA_bayesian[1][4]->SetBinContent(2,0.7609972);
  PASRAA_bayesian[1][4]->SetBinContent(3,0.7759888);
  PASRAA_bayesian[1][4]->SetBinContent(4,0.7765765);
  PASRAA_bayesian[1][4]->SetBinContent(5,0.7793122);
  PASRAA_bayesian[1][4]->SetBinContent(6,0.7750514);
  PASRAA_bayesian[1][4]->SetBinContent(7,0.8105651);
  PASRAA_bayesian[1][4]->SetBinContent(8,0.778167);
  PASRAA_bayesian[1][4]->SetBinContent(9,0.6663016);
  PASRAA_bayesian[1][4]->SetBinContent(10,0.736169);
  PASRAA_bayesian[1][4]->SetBinContent(11,0.8040183);
  PASRAA_bayesian[1][4]->SetBinContent(12,0.9631981);

  PASRAA_bayesian[1][4]->SetBinError(0,0.003463908);
  PASRAA_bayesian[1][4]->SetBinError(1,0.01638808);
  PASRAA_bayesian[1][4]->SetBinError(2,0.02150757);
  PASRAA_bayesian[1][4]->SetBinError(3,0.02926648);
  PASRAA_bayesian[1][4]->SetBinError(4,0.0388617);
  PASRAA_bayesian[1][4]->SetBinError(5,0.05095474);
  PASRAA_bayesian[1][4]->SetBinError(6,0.0616764);
  PASRAA_bayesian[1][4]->SetBinError(7,0.07817987);
  PASRAA_bayesian[1][4]->SetBinError(8,0.09560783);
  PASRAA_bayesian[1][4]->SetBinError(9,0.07394511);
  PASRAA_bayesian[1][4]->SetBinError(10,0.09310688);
  PASRAA_bayesian[1][4]->SetBinError(11,0.1812772);
  PASRAA_bayesian[1][4]->SetBinError(12,0.2031483);
   
  // PASRAA_bayesian[1][4]->SetBinError(0,0.003463908);
  // PASRAA_bayesian[1][4]->SetBinError(1,0.008645942);
  // PASRAA_bayesian[1][4]->SetBinError(2,0.01082537);
  // PASRAA_bayesian[1][4]->SetBinError(3,0.01412537);
  // PASRAA_bayesian[1][4]->SetBinError(4,0.01806789);
  // PASRAA_bayesian[1][4]->SetBinError(5,0.02291484);
  // PASRAA_bayesian[1][4]->SetBinError(6,0.02692976);
  // PASRAA_bayesian[1][4]->SetBinError(7,0.03325802);
  // PASRAA_bayesian[1][4]->SetBinError(8,0.03975427);
  // PASRAA_bayesian[1][4]->SetBinError(9,0.02987937);
  // PASRAA_bayesian[1][4]->SetBinError(10,0.03618969);
  // PASRAA_bayesian[1][4]->SetBinError(11,0.06938782);
  // PASRAA_bayesian[1][4]->SetBinError(12,0.2031483);
  PASRAA_bayesian[1][4]->SetMinimum(0);
  PASRAA_bayesian[1][4]->SetMaximum(2);
  PASRAA_bayesian[1][4]->SetEntries(6061.544);

  PASRAA_bayesian[1][3] = new TH1F("PASRAA_bayesian_R3_cent3","",11,xAxis2);
  PASRAA_bayesian[1][3]->SetBinContent(0,0.7279794);
  PASRAA_bayesian[1][3]->SetBinContent(1,0.6681795);
  PASRAA_bayesian[1][3]->SetBinContent(2,0.6799733);
  PASRAA_bayesian[1][3]->SetBinContent(3,0.6649472);
  PASRAA_bayesian[1][3]->SetBinContent(4,0.6252716);
  PASRAA_bayesian[1][3]->SetBinContent(5,0.6117619);
  PASRAA_bayesian[1][3]->SetBinContent(6,0.6540294);
  PASRAA_bayesian[1][3]->SetBinContent(7,0.7157548);
  PASRAA_bayesian[1][3]->SetBinContent(8,0.6989863);
  PASRAA_bayesian[1][3]->SetBinContent(9,0.5692854);
  PASRAA_bayesian[1][3]->SetBinContent(10,0.6007913);
  PASRAA_bayesian[1][3]->SetBinContent(11,0.8094043);
  PASRAA_bayesian[1][3]->SetBinContent(12,2.021915);

  PASRAA_bayesian[1][3]->SetBinError(0,0.002593759);
  PASRAA_bayesian[1][3]->SetBinError(1,0.01113931);
  PASRAA_bayesian[1][3]->SetBinError(2,0.01478439);
  PASRAA_bayesian[1][3]->SetBinError(3,0.0194451);
  PASRAA_bayesian[1][3]->SetBinError(4,0.02369379);
  PASRAA_bayesian[1][3]->SetBinError(5,0.03078356);
  PASRAA_bayesian[1][3]->SetBinError(6,0.03982553);
  PASRAA_bayesian[1][3]->SetBinError(7,0.05365841);
  PASRAA_bayesian[1][3]->SetBinError(8,0.0640314);
  PASRAA_bayesian[1][3]->SetBinError(9,0.04796186);
  PASRAA_bayesian[1][3]->SetBinError(10,0.05829843);
  PASRAA_bayesian[1][3]->SetBinError(11,0.1392753);
  PASRAA_bayesian[1][3]->SetBinError(12,0.2970788);

  // PASRAA_bayesian[1][3]->SetBinError(0,0.002593759);
  // PASRAA_bayesian[1][3]->SetBinError(1,0.005876821);
  // PASRAA_bayesian[1][3]->SetBinError(2,0.007441405);
  // PASRAA_bayesian[1][3]->SetBinError(3,0.009385111);
  // PASRAA_bayesian[1][3]->SetBinError(4,0.01101591);
  // PASRAA_bayesian[1][3]->SetBinError(5,0.01384366);
  // PASRAA_bayesian[1][3]->SetBinError(6,0.01738902);
  // PASRAA_bayesian[1][3]->SetBinError(7,0.0228265);
  // PASRAA_bayesian[1][3]->SetBinError(8,0.02662461);
  // PASRAA_bayesian[1][3]->SetBinError(9,0.01938019);
  // PASRAA_bayesian[1][3]->SetBinError(10,0.02266);
  // PASRAA_bayesian[1][3]->SetBinError(11,0.05331068);
  // PASRAA_bayesian[1][3]->SetBinError(12,0.2970788);
  PASRAA_bayesian[1][3]->SetMinimum(0);
  PASRAA_bayesian[1][3]->SetMaximum(2);
  PASRAA_bayesian[1][3]->SetEntries(9256.671);

  PASRAA_bayesian[1][2] = new TH1F("PASRAA_bayesian_R3_cent2","",11,xAxis2);
  PASRAA_bayesian[1][2]->SetBinContent(0,0.8238835);
  PASRAA_bayesian[1][2]->SetBinContent(1,0.5604139);
  PASRAA_bayesian[1][2]->SetBinContent(2,0.5365105);
  PASRAA_bayesian[1][2]->SetBinContent(3,0.5556849);
  PASRAA_bayesian[1][2]->SetBinContent(4,0.5274078);
  PASRAA_bayesian[1][2]->SetBinContent(5,0.5573406);
  PASRAA_bayesian[1][2]->SetBinContent(6,0.6040812);
  PASRAA_bayesian[1][2]->SetBinContent(7,0.6451674);
  PASRAA_bayesian[1][2]->SetBinContent(8,0.5992499);
  PASRAA_bayesian[1][2]->SetBinContent(9,0.4717564);
  PASRAA_bayesian[1][2]->SetBinContent(10,0.4812902);
  PASRAA_bayesian[1][2]->SetBinContent(11,0.6858578);
  PASRAA_bayesian[1][2]->SetBinContent(12,1.051401);

  PASRAA_bayesian[1][2]->SetBinError(0,0.002676996);
  PASRAA_bayesian[1][2]->SetBinError(1,0.008546944);
  PASRAA_bayesian[1][2]->SetBinError(2,0.01060751);
  PASRAA_bayesian[1][2]->SetBinError(3,0.01473257);
  PASRAA_bayesian[1][2]->SetBinError(4,0.01825233);
  PASRAA_bayesian[1][2]->SetBinError(5,0.02542376);
  PASRAA_bayesian[1][2]->SetBinError(6,0.03367376);
  PASRAA_bayesian[1][2]->SetBinError(7,0.04421273);
  PASRAA_bayesian[1][2]->SetBinError(8,0.04964768);
  PASRAA_bayesian[1][2]->SetBinError(9,0.03594492);
  PASRAA_bayesian[1][2]->SetBinError(10,0.04249006);
  PASRAA_bayesian[1][2]->SetBinError(11,0.1093235);
  PASRAA_bayesian[1][2]->SetBinError(12,0.1496494);
  
  // PASRAA_bayesian[1][2]->SetBinError(0,0.002676996);
  // PASRAA_bayesian[1][2]->SetBinError(1,0.004509155);
  // PASRAA_bayesian[1][2]->SetBinError(2,0.005339063);
  // PASRAA_bayesian[1][2]->SetBinError(3,0.007110624);
  // PASRAA_bayesian[1][2]->SetBinError(4,0.008486022);
  // PASRAA_bayesian[1][2]->SetBinError(5,0.01143331);
  // PASRAA_bayesian[1][2]->SetBinError(6,0.01470297);
  // PASRAA_bayesian[1][2]->SetBinError(7,0.01880827);
  // PASRAA_bayesian[1][2]->SetBinError(8,0.02064378);
  // PASRAA_bayesian[1][2]->SetBinError(9,0.01452445);
  // PASRAA_bayesian[1][2]->SetBinError(10,0.01651545);
  // PASRAA_bayesian[1][2]->SetBinError(11,0.04184596);
  // PASRAA_bayesian[1][2]->SetBinError(12,0.1496494);
  PASRAA_bayesian[1][2]->SetMinimum(0);
  PASRAA_bayesian[1][2]->SetMaximum(2);
  PASRAA_bayesian[1][2]->SetEntries(10967.25);

  PASRAA_bayesian[1][1] = new TH1F("PASRAA_bayesian_R3_cent1","",11,xAxis2);
  PASRAA_bayesian[1][1]->SetBinContent(0,1.031868);
  PASRAA_bayesian[1][1]->SetBinContent(1,0.4557194);
  PASRAA_bayesian[1][1]->SetBinContent(2,0.480005);
  PASRAA_bayesian[1][1]->SetBinContent(3,0.4997296);
  PASRAA_bayesian[1][1]->SetBinContent(4,0.444041);
  PASRAA_bayesian[1][1]->SetBinContent(5,0.4247413);
  PASRAA_bayesian[1][1]->SetBinContent(6,0.4454679);
  PASRAA_bayesian[1][1]->SetBinContent(7,0.4785669);
  PASRAA_bayesian[1][1]->SetBinContent(8,0.4834603);
  PASRAA_bayesian[1][1]->SetBinContent(9,0.4274004);
  PASRAA_bayesian[1][1]->SetBinContent(10,0.4192492);
  PASRAA_bayesian[1][1]->SetBinContent(11,0.5767608);
  PASRAA_bayesian[1][1]->SetBinContent(12,1.223037);

  PASRAA_bayesian[1][1]->SetBinError(0,0.003550127);
  PASRAA_bayesian[1][1]->SetBinError(1,0.007539356);
  PASRAA_bayesian[1][1]->SetBinError(2,0.01035939);
  PASRAA_bayesian[1][1]->SetBinError(3,0.01458809);
  PASRAA_bayesian[1][1]->SetBinError(4,0.01679069);
  PASRAA_bayesian[1][1]->SetBinError(5,0.02080943);
  PASRAA_bayesian[1][1]->SetBinError(6,0.02718877);
  PASRAA_bayesian[1][1]->SetBinError(7,0.03569299);
  PASRAA_bayesian[1][1]->SetBinError(8,0.04427906);
  PASRAA_bayesian[1][1]->SetBinError(9,0.03604142);
  PASRAA_bayesian[1][1]->SetBinError(10,0.0406791);
  PASRAA_bayesian[1][1]->SetBinError(11,0.09958546);
  PASRAA_bayesian[1][1]->SetBinError(12,0.1827313);
  
  // PASRAA_bayesian[1][1]->SetBinError(0,0.003550127);
  // PASRAA_bayesian[1][1]->SetBinError(1,0.003977577);
  // PASRAA_bayesian[1][1]->SetBinError(2,0.005214178);
  // PASRAA_bayesian[1][1]->SetBinError(3,0.00704089);
  // PASRAA_bayesian[1][1]->SetBinError(4,0.007806462);
  // PASRAA_bayesian[1][1]->SetBinError(5,0.0093582);
  // PASRAA_bayesian[1][1]->SetBinError(6,0.01187143);
  // PASRAA_bayesian[1][1]->SetBinError(7,0.01518394);
  // PASRAA_bayesian[1][1]->SetBinError(8,0.01841148);
  // PASRAA_bayesian[1][1]->SetBinError(9,0.01456344);
  // PASRAA_bayesian[1][1]->SetBinError(10,0.01581155);
  // PASRAA_bayesian[1][1]->SetBinError(11,0.03811852);
  // PASRAA_bayesian[1][1]->SetBinError(12,0.1827313);
  PASRAA_bayesian[1][1]->SetMinimum(0);
  PASRAA_bayesian[1][1]->SetMaximum(2);
  PASRAA_bayesian[1][1]->SetEntries(9198.666);

  PASRAA_bayesian[1][0] = new TH1F("PASRAA_bayesian_R3_cent0","",11,xAxis2);
  PASRAA_bayesian[1][0]->SetBinContent(0,1.418427);
  PASRAA_bayesian[1][0]->SetBinContent(1,0.4702078);
  PASRAA_bayesian[1][0]->SetBinContent(2,0.4736779);
  PASRAA_bayesian[1][0]->SetBinContent(3,0.4640111);
  PASRAA_bayesian[1][0]->SetBinContent(4,0.4361698);
  PASRAA_bayesian[1][0]->SetBinContent(5,0.4342382);
  PASRAA_bayesian[1][0]->SetBinContent(6,0.4725478);
  PASRAA_bayesian[1][0]->SetBinContent(7,0.5191602);
  PASRAA_bayesian[1][0]->SetBinContent(8,0.5139099);
  PASRAA_bayesian[1][0]->SetBinContent(9,0.4391448);
  PASRAA_bayesian[1][0]->SetBinContent(10,0.4229113);
  PASRAA_bayesian[1][0]->SetBinContent(11,0.5843745);
  PASRAA_bayesian[1][0]->SetBinContent(12,0.7831888);

  PASRAA_bayesian[1][0]->SetBinError(0,0.004730243);
  PASRAA_bayesian[1][0]->SetBinError(1,0.007550343);
  PASRAA_bayesian[1][0]->SetBinError(2,0.009840934);
  PASRAA_bayesian[1][0]->SetBinError(3,0.01303877);
  PASRAA_bayesian[1][0]->SetBinError(4,0.01589791);
  PASRAA_bayesian[1][0]->SetBinError(5,0.02094939);
  PASRAA_bayesian[1][0]->SetBinError(6,0.02765753);
  PASRAA_bayesian[1][0]->SetBinError(7,0.03752524);
  PASRAA_bayesian[1][0]->SetBinError(8,0.04478431);
  PASRAA_bayesian[1][0]->SetBinError(9,0.03525832);
  PASRAA_bayesian[1][0]->SetBinError(10,0.03945337);
  PASRAA_bayesian[1][0]->SetBinError(11,0.09843785);
  PASRAA_bayesian[1][0]->SetBinError(12,0.1174314);
  
  // PASRAA_bayesian[1][0]->SetBinError(0,0.004730243);
  // PASRAA_bayesian[1][0]->SetBinError(1,0.003983373);
  // PASRAA_bayesian[1][0]->SetBinError(2,0.004953222);
  // PASRAA_bayesian[1][0]->SetBinError(3,0.006293119);
  // PASRAA_bayesian[1][0]->SetBinError(4,0.007391385);
  // PASRAA_bayesian[1][0]->SetBinError(5,0.009421142);
  // PASRAA_bayesian[1][0]->SetBinError(6,0.0120761);
  // PASRAA_bayesian[1][0]->SetBinError(7,0.01596338);
  // PASRAA_bayesian[1][0]->SetBinError(8,0.01862156);
  // PASRAA_bayesian[1][0]->SetBinError(9,0.01424701);
  // PASRAA_bayesian[1][0]->SetBinError(10,0.01533512);
  // PASRAA_bayesian[1][0]->SetBinError(11,0.03767925);
  // PASRAA_bayesian[1][0]->SetBinError(12,0.1174314);
  PASRAA_bayesian[1][0]->SetMinimum(0);
  PASRAA_bayesian[1][0]->SetMaximum(2);
  PASRAA_bayesian[1][0]->SetEntries(9671.122);

  //now that we have taken the necessary histogram, lets start to make plots:

  //plot1 - 6 panel plot of RAA, at R=0.3, in the old PAS and now.
  TCanvas *cRAA_PASComp = new TCanvas("cRAA_PASComp","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA_PASComp,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA_PASComp = myLegend(0.15,0.75,0.85,0.9);
  TLine *lineRAA = new TLine(65,1,299,1);
  lineRAA->SetLineStyle(2);
  lineRAA->SetLineWidth(2);

  int ci;
  TBox *box;
  for(int i = 0;i<nbins_cent;i++){

    cRAA_PASComp->cd(nbins_cent-i);

    RAA_R3_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R3_Bayes[i]->SetMarkerStyle(20);
    makeHistTitle( RAA_R3_Bayes[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    RAA_R3_Bayes[i]->SetAxisRange(60,299,"X");
    RAA_R3_Bayes[i]->SetAxisRange(0,2,"Y");
    RAA_R3_Bayes[i]->Draw("E0");

    PASRAA_bayesian[1][i]->SetMarkerColor(kRed);
    PASRAA_bayesian[1][i]->SetMarkerStyle(33);
    PASRAA_bayesian[1][i]->Draw("same E0");

    //RAA_binbybin[1][i]->SetMarkerStyle(29);
    //RAA_binbybin[1][i]->SetMarkerColor(kBlue);
    //RAA_binbybin[1][i]->Draw("same");

    lineRAA->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

    switch ( i ) {
	
    case 0:
      box = new TBox(100,0.3992887,110,0.5411268);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(110,0.4004558,120,0.5469);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(120,0.3797207,130,0.5483015);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(130,0.3460133,140,0.5263263);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(140,0.3479742,150,0.5205022);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(150,0.3878594,160,0.5572361);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(160,0.4375477,170,0.6007726);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(170,0.432544,180,0.5952758);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(180,0.3688692,200,0.5094203);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(200,0.3537775,240,0.492045);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(240,0.485418,300,0.6833309);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
	
      break;

    case 1:
      box = new TBox(100,0.3922236,110,0.5192152);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(110,0.4114913,120,0.5485187);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(120,0.4142097,130,0.5852495);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(130,0.3565539,140,0.5315281);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(140,0.3448709,150,0.5046117);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(150,0.3711734,160,0.5197624);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(160,0.4105411,170,0.5465928);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(170,0.4144865,180,0.5524341);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(180,0.3660878,200,0.4887129);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(200,0.3584384,240,0.4800599);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(240,0.4915343,300,0.6619874);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      break;

    case 2:
      box = new TBox(100,0.4885437,110,0.6322841);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(110,0.4659843,120,0.6070368);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(120,0.4659718,130,0.6453979);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(130,0.4280813,140,0.6267344);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(140,0.4578278,150,0.6568533);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(150,0.5101032,160,0.6980592);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(160,0.5623943,170,0.7279404);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(170,0.5223327,180,0.6761671);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(180,0.4111619,200,0.5323508);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(200,0.4193837,240,0.5431968);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(240,0.5974073,300,0.7743082);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();

      break;

    case 3:

      box = new TBox(100,0.5888569,110,0.7475021);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(110,0.5972088,120,0.7627378);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(120,0.5630347,130,0.7668597);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(130,0.5120569,140,0.7384863);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(140,0.507423,150,0.7161008);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(150,0.5585642,160,0.7494945);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(160,0.6327194,170,0.7987902);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(170,0.6181627,180,0.7798098);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(180,0.5037794,200,0.6347914);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(200,0.532314,240,0.6692686);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(240,0.7185097,300,0.9002989);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
	
      break;

    case 4:
      box = new TBox(100,0.6653466,110,0.8308992);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(110,0.6746653,120,0.8473291);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(120,0.662311,130,0.8896665);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(130,0.6405509,140,0.9126022);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(140,0.651443,150,0.9071814);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(150,0.6679778,160,0.882125);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(160,0.7248008,170,0.8963293);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(170,0.6962772,180,0.8600568);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(180,0.596691,200,0.7359122);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(200,0.6601083,240,0.8122298);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(240,0.7215662,300,0.8864704);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();

      break;

    case 5:
      box = new TBox(100,0.7705583,110,0.9495451);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(110,0.7227901,120,0.8972415);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(120,0.6428426,130,0.8565895);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(130,0.602727,140,0.853577);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(140,0.6775657,150,0.9382831);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(150,0.9043127,160,1.187379);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(160,0.8928045,170,1.097036);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(170,0.762919,180,0.9372863);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(180,0.6536744,200,0.8029628);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(200,0.762955,240,0.9373031);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
      box = new TBox(240,0.7678421,300,0.9445981);

      ci = TColor::GetColor("#cccccc");
      box->SetFillColor(ci);

      ci = TColor::GetColor("#cccccc");
      box->SetLineColor(ci);
      box->Draw();
		
      break;
	
    }// switch


    RAA_R3_Bayes[i]->Draw("same E1");
    systematics_R3.Draw(RAA_R3_Bayes[i],i,1);

    PASRAA_bayesian[1][i]->Draw("same E1");

  }// centrality loop

  tRAA_PASComp->AddEntry( RAA_R3_Bayes[0],"13-005, New JetID Cut,","pl");
  tRAA_PASComp->AddEntry(PASRAA_bayesian[1][0],"12-004, trkMax/jtpt > 0.01","pl");
  tRAA_PASComp->SetTextSize(0.04);

  cRAA_PASComp->cd(1);
  tRAA_PASComp->Draw();
  drawText("Bayesian Unfolding, 4 iterations",0.2,0.7,16);
  cRAA_PASComp->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.3",algo,jet_type),0.2,0.23,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_PASComp->cd(2);
  drawText("|#eta|<2",0.1,0.3,16);
  //drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  //cRAA_PASComp->cd(3);
  //drawText("Trig Combined with MinBias&&!Jet80 subtracted",0.06,0.2,16);

  cRAA_PASComp->SaveAs(Form("June4/Analysis_full_bins/RAA_newComparingwithPAS_ak%s3%s_%d.pdf",algo,jet_type,date.GetDate()),"RECREATE");
  cRAA_PASComp->SaveAs(Form("June4/Analysis_full_bins/RAA_newComparingwithPAS_ak%s3%s_%d.C",algo,jet_type,date.GetDate()),"RECREATE");
  cRAA_PASComp->SaveAs(Form("June4/Analysis_full_bins/RAA_newComparingwithPAS_ak%s3%s_%d.root",algo,jet_type,date.GetDate()),"RECREATE");
  
  


  ///


    
  // draw it for R=0.3
  TCanvas *cRAA_R3 = new TCanvas("cRAA_R3","RAA",1400,1200);
  makeMultiPanelCanvas(cRAA_R3,3,2,0.0,0.0,0.2,0.15,0.07);

  TLegend *tRAA_R3 = myLegend(0.35,0.55,0.65,0.7);
  TLine *lineRAA_R3 = new TLine(unfoldingCut+30,1,299,1);
  lineRAA_R3->SetLineStyle(2);
  lineRAA_R3->SetLineWidth(2);
    
  for(int i = 0;i<nbins_cent;++i){

    cRAA_R3->cd(nbins_cent-i);

    // if(i==0){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(100); ++j){
    // 	RAA_R3_Meas[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R3_Meas[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==5 || i==4 || i==3){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(240); j<RAA_R3_Bayes[i]->FindBin(300); ++j){
    // 	RAA_R3_Meas[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R3_Meas[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==1 || i==2){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(unfoldingCut+30); ++j){
    // 	RAA_R3_Meas[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R3_Meas[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_BinByBin[i]->SetBinError(j,0);
    //   }
    // }

    RAA_R3_Meas[i]->SetMarkerColor(kBlack);
    RAA_R3_Meas[i]->SetMarkerStyle(24);
    makeHistTitle(RAA_R3_Meas[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    RAA_R3_Meas[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    RAA_R3_Meas[i]->SetAxisRange(0,2,"Y");
    RAA_R3_Meas[i]->Draw("E1");

    RAA_R3_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R3_Bayes[i]->SetMarkerStyle(20);
    RAA_R3_Bayes[i]->Draw("same E1");

    RAA_R3_BinByBin[i]->SetMarkerStyle(33);
    RAA_R3_BinByBin[i]->SetMarkerColor(kRed);
    RAA_R3_BinByBin[i]->Draw("same E1");

    lineRAA_R3->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.85,20);

    systematics_R3.calcTotalSys(i);
    systematics_R3.Draw(RAA_R3_Bayes[i],i,2);
    
    RAA_R3_Meas[i]->Draw("same E1");
    RAA_R3_BinByBin[i]->Draw("same E1");
    RAA_R3_Bayes[i]->Draw("same E1");

  }
    
  tRAA_R3->AddEntry(RAA_R3_Meas[0],"Measured","pl");
  tRAA_R3->AddEntry(RAA_R3_Bayes[0],"Bayesian","pl");
  tRAA_R3->AddEntry(RAA_R3_BinByBin[0],"bin-by-bin","pl");
  tRAA_R3->SetTextSize(0.04);

  cRAA_R3->cd(1);
  tRAA_R3->Draw();
  cRAA_R3->cd(1);
  putCMSPrel();
  cRAA_R3->cd(2);
  putPbPbLumi();
  cRAA_R3->cd(3);
  putPPLumi();
  cRAA_R3->cd(1);
  drawText(Form("Anti-k_{T} %s R=0.3 %s Jets",algo,jet_type),0.25,0.20,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_R3->cd(2);
  drawText(Form("Jet ID cut, |#eta|<%2.0f",etaBoundary),0.1,0.2,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.1,16);
  cRAA_R3->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.2,16);
  drawText("Pile up rejection cut applied",0.1,0.1,16);

  cRAA_R3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R3_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_R3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R3_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_R3->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R3_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");

  // plot it for R=0.4
  TCanvas *cRAA_R4 = new TCanvas("cRAA_R4","RAA",1400,1200);
  makeMultiPanelCanvas(cRAA_R4,3,2,0.0,0.0,0.2,0.15,0.07);

  TLegend *tRAA_R4 = myLegend(0.35,0.55,0.65,0.7);
  TLine *lineRAA_R4 = new TLine(65,1,299,1);
  lineRAA_R4->SetLineStyle(2);
  lineRAA_R4->SetLineWidth(2);
    
  for(int i = 0;i<nbins_cent;++i){

    cRAA_R4->cd(nbins_cent-i);

    // if(i==0){
    //   for(int j = RAA_R4_Bayes[i]->FindBin(unfoldingCut); j<RAA_R4_Bayes[i]->FindBin(100); ++j){
    // 	RAA_R4_Meas[i]->SetBinContent(j,0);
    // 	RAA_R4_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R4_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R4_Meas[i]->SetBinError(j,0);
    // 	RAA_R4_Bayes[i]->SetBinError(j,0);
    // 	RAA_R4_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==5 || i==4 || i==3){
    //   for(int j = RAA_R4_Bayes[i]->FindBin(240); j<RAA_R4_Bayes[i]->FindBin(300); ++j){
    // 	RAA_R4_Meas[i]->SetBinContent(j,0);
    // 	RAA_R4_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R4_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R4_Meas[i]->SetBinError(j,0);
    // 	RAA_R4_Bayes[i]->SetBinError(j,0);
    // 	RAA_R4_BinByBin[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==1 || i==2){
    //   for(int j = RAA_R4_Bayes[i]->FindBin(unfoldingCut); j<RAA_R4_Bayes[i]->FindBin(unfoldingCut+30); ++j){
    // 	RAA_R4_Meas[i]->SetBinContent(j,0);
    // 	RAA_R4_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R4_BinByBin[i]->SetBinContent(j,0);
    // 	RAA_R4_Meas[i]->SetBinError(j,0);
    // 	RAA_R4_Bayes[i]->SetBinError(j,0);
    // 	RAA_R4_BinByBin[i]->SetBinError(j,0);
    //   }
    // }

    RAA_R4_Meas[i]->SetMarkerColor(kBlack);
    RAA_R4_Meas[i]->SetMarkerStyle(24);
    makeHistTitle(RAA_R4_Meas[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    RAA_R4_Meas[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    RAA_R4_Meas[i]->SetAxisRange(0,2,"Y");
    RAA_R4_Meas[i]->Draw("E1");

    RAA_R4_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R4_Bayes[i]->SetMarkerStyle(20);
    RAA_R4_Bayes[i]->Draw("same E1");

    RAA_R4_BinByBin[i]->SetMarkerStyle(33);
    RAA_R4_BinByBin[i]->SetMarkerColor(kRed);
    RAA_R4_BinByBin[i]->Draw("same E1");

    lineRAA_R4->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.85,20);

    systematics_R4.calcTotalSys(i);
    systematics_R4.Draw(RAA_R4_Bayes[i],i,2);
    
    RAA_R4_Meas[i]->Draw("same E1");
    RAA_R4_BinByBin[i]->Draw("same E1");
    RAA_R4_Bayes[i]->Draw("same E1");

  }
    
  tRAA_R4->AddEntry(RAA_R4_Meas[0],"Measured","pl");
  tRAA_R4->AddEntry(RAA_R4_Bayes[0],"Bayesian","pl");
  tRAA_R4->AddEntry(RAA_R4_BinByBin[0],"bin-by-bin","pl");
  tRAA_R4->SetTextSize(0.04);

  cRAA_R4->cd(1);
  tRAA_R4->Draw();
  cRAA_R4->cd(1);
  putCMSPrel();
  cRAA_R4->cd(2);
  putPbPbLumi();
  cRAA_R4->cd(3);
  putPPLumi();
  cRAA_R4->cd(1);
  drawText(Form("Anti-k_{T} %s R=0.4 %s Jets",algo,jet_type),0.25,0.20,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_R4->cd(2);
  drawText(Form("Jet ID cut, |#eta|<%2.0f",etaBoundary),0.1,0.2,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.1,16);
  cRAA_R4->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.2,16);
  drawText("Pile up rejection cut applied",0.1,0.1,16);

  cRAA_R4->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R4_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA_R4->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R4_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA_R4->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_R4_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");
  
  // plot 2 - Bayesian unfolded RAA as a function of pT for the different radii
  //        - regular 6 panel plot 
  
  // again this will be a 6 panel plot. showing measured, unfolded Bayesian, and unfolded Bin By Bin methods. 
  TCanvas *cRAA = new TCanvas("cRAA","RAA",1400,1200);
  makeMultiPanelCanvas(cRAA,3,2,0.0,0.0,0.2,0.15,0.07);

  TLegend *tRAA = myLegend(0.35,0.55,0.65,0.7);
  //TLine *lineRAA = new TLine(65,1,299,1);
  lineRAA->SetLineStyle(2);
  lineRAA->SetLineWidth(2);

  TLine *lUnfoldingCut = new TLine(unfoldingCut+30,0,unfoldingCut+30,2);
  lUnfoldingCut->SetLineStyle(4);
  lUnfoldingCut->SetLineWidth(2);
    
  for(int i = 0;i<nbins_cent;++i){

    cRAA->cd(nbins_cent-i);

    // if(i==0){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(100); ++j){
    // 	RAA_R2_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R4_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R2_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R4_Bayes[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==5 || i==4 || i==3){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(240); j<RAA_R3_Bayes[i]->FindBin(300); ++j){
    // 	RAA_R2_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R4_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R2_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R4_Bayes[i]->SetBinError(j,0);
    //   }
    // }
    // if(i==1 || i==2){
    //   for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(unfoldingCut+30); ++j){
    // 	RAA_R2_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R3_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R4_Bayes[i]->SetBinContent(j,0);
    // 	RAA_R2_Bayes[i]->SetBinError(j,0);
    // 	RAA_R3_Bayes[i]->SetBinError(j,0);
    // 	RAA_R4_Bayes[i]->SetBinError(j,0);
    //   }
    // }

    RAA_R2_Bayes[i]->SetMarkerColor(kRed);
    RAA_R2_Bayes[i]->SetMarkerStyle(20);
    makeHistTitle(RAA_R2_Bayes[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    RAA_R2_Bayes[i]->SetAxisRange(65, boundaries_pt[nbins_pt+1],"X");
    RAA_R2_Bayes[i]->SetAxisRange(0,2,"Y");
    RAA_R2_Bayes[i]->Draw("E1");

    RAA_R3_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R3_Bayes[i]->SetMarkerStyle(20);
    RAA_R3_Bayes[i]->Draw("same E1");

    RAA_R4_Bayes[i]->SetMarkerStyle(20);
    RAA_R4_Bayes[i]->SetMarkerColor(kBlue);
    RAA_R4_Bayes[i]->Draw("same E1");

    lineRAA->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.85,20);

    systematics_R2.calcTotalSys(i);
    systematics_R2.Draw(RAA_R2_Bayes[i],i,2);
    systematics_R3.calcTotalSys(i);
    systematics_R3.Draw(RAA_R3_Bayes[i],i,1);
    systematics_R4.calcTotalSys(i);
    systematics_R4.Draw(RAA_R4_Bayes[i],i,4);
    
    RAA_R2_Bayes[i]->Draw("same E1");
    RAA_R4_Bayes[i]->Draw("same E1");
    RAA_R3_Bayes[i]->Draw("same E1");

  }
    
  tRAA->AddEntry(RAA_R2_Bayes[0],"R=0.2","pl");
  tRAA->AddEntry(RAA_R3_Bayes[0],"R=0.3","pl");
  tRAA->AddEntry(RAA_R4_Bayes[0],"R=0.4","pl");
  tRAA->SetTextSize(0.04);

  cRAA->cd(1);
  tRAA->Draw();
  cRAA->cd(1);
  putCMSPrel();
  cRAA->cd(2);
  putPbPbLumi();
  cRAA->cd(3);
  putPPLumi();
  cRAA->cd(1);
  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.25,0.20,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA->cd(2);
  drawText(Form("Jet ID cut, |#eta|<%2.0f",etaBoundary),0.1,0.2,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.1,16);
  cRAA->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.2,16);
  drawText("Pile up rejection cut applied",0.1,0.1,16);

  cRAA->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cRAA->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cRAA->SaveAs(Form("June4/Analysis_full_bins/Final_paper_plots_RAA_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");
    

  
  // plot - comparison with ATLAS
  // get the ATLAS most central bin 0-10% RAA from the hepdata http://hepdata.cedar.ac.uk/view/ins1326911/next

  
  // Plot: p8719_d27x1y1 - 0-10%
  double p8719_d27x1y1_xval[] = { 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5, 357.0 };
  double p8719_d27x1y1_xerrminus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d27x1y1_xerrplus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d27x1y1_yval[] = { 0.472, 0.491, 0.453, 0.478, 0.511, 0.565, 0.598, 0.595, 0.563 };
  double p8719_d27x1y1_yerrminus[] = { 0.05435587740070064, 0.08343822336315654, 0.06473622013988768, 0.052501723552660626, 0.05534391863610671, 0.06402723990615243, 0.07539070218534909, 0.09135526996293097, 0.13564678562354507 };
  double p8719_d27x1y1_yerrplus[] = { 0.05435587740070064, 0.08343822336315654, 0.06473622013988768, 0.052501723552660626, 0.05534391863610671, 0.06402723990615243, 0.07539070218534909, 0.09135526996293097, 0.13564678562354507 };
  double p8719_d27x1y1_ystatminus[] = { 0.011799999999999998, 0.010802, 0.008154, 0.007647999999999999, 0.008176000000000001, 0.0113, 0.017342, 0.031534999999999994, 0.053485 };
  double p8719_d27x1y1_ystatplus[] = { 0.011799999999999998, 0.010802, 0.008154, 0.007647999999999999, 0.008176000000000001, 0.0113, 0.017342, 0.031534999999999994, 0.053485 };
  int p8719_d27x1y1_numpoints = 9;
  TGraphAsymmErrors * p8719_d27x1y1 = new TGraphAsymmErrors(p8719_d27x1y1_numpoints, p8719_d27x1y1_xval, p8719_d27x1y1_yval, p8719_d27x1y1_xerrminus, p8719_d27x1y1_xerrplus, p8719_d27x1y1_yerrminus, p8719_d27x1y1_yerrplus);
  p8719_d27x1y1->SetName("/HepData/8719/d27x1y1");
  p8719_d27x1y1->SetTitle(" ");
  p8719_d27x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");
  p8719_d27x1y1->GetYaxis()->SetTitle("R_{AA}");
  

  // Plot: p8719_d28x1y1 - 10-20% 
  double p8719_d28x1y1_xval[] = { 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d28x1y1_xerrminus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d28x1y1_xerrplus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d28x1y1_yval[] = { 0.525, 0.529, 0.496, 0.543, 0.58, 0.651, 0.653, 0.657 };
  double p8719_d28x1y1_yerrminus[] = { 0.07913406346195044, 0.07710708064762925, 0.04466754705600028, 0.0552073737194589, 0.05947475094525406, 0.0699484277521661, 0.0784823419948717, 0.10211014442747597 };
  double p8719_d28x1y1_yerrplus[] = { 0.07913406346195044, 0.07710708064762925, 0.04466754705600028, 0.0552073737194589, 0.05947475094525406, 0.0699484277521661, 0.0784823419948717, 0.10211014442747597 };
  double p8719_d28x1y1_ystatminus[] = { 0.0126, 0.013754000000000002, 0.008928, 0.008688000000000001, 0.009859999999999999, 0.014322000000000001, 0.022202, 0.040077 };
  double p8719_d28x1y1_ystatplus[] = { 0.0126, 0.013754000000000002, 0.008928, 0.008688000000000001, 0.009859999999999999, 0.014322000000000001, 0.022202, 0.040077 };
  int p8719_d28x1y1_numpoints = 8;
  TGraphAsymmErrors * p8719_d28x1y1 = new TGraphAsymmErrors(p8719_d28x1y1_numpoints, p8719_d28x1y1_xval, p8719_d28x1y1_yval, p8719_d28x1y1_xerrminus, p8719_d28x1y1_xerrplus, p8719_d28x1y1_yerrminus, p8719_d28x1y1_yerrplus);
  p8719_d28x1y1->SetName("/HepData/8719/d28x1y1");
  p8719_d28x1y1->SetTitle(" ");
  p8719_d28x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");
  p8719_d28x1y1->GetYaxis()->SetTitle("R_{AA}");

  // Plot: p8719_d29x1y1 20-30% 
  double p8719_d29x1y1_xval[] = { 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d29x1y1_xerrminus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d29x1y1_xerrplus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d29x1y1_yval[] = { 0.525, 0.612, 0.569, 0.603, 0.646, 0.668, 0.678, 0.693 };
  double p8719_d29x1y1_yerrminus[] = { 0.04394033596821945, 0.04984325783894949, 0.052726992205131516, 0.05972135866672828, 0.06435405180717063, 0.07073864527964895, 0.07955745719918404, 0.09429384914192442 };
  double p8719_d29x1y1_yerrplus[] = { 0.04394033596821945, 0.04984325783894949, 0.052726992205131516, 0.05972135866672828, 0.06435405180717063, 0.07073864527964895, 0.07955745719918404, 0.09429384914192442 };
  double p8719_d29x1y1_ystatminus[] = { 0.0105, 0.011016, 0.008535, 0.007839, 0.010336000000000001, 0.014696, 0.025764000000000002, 0.043658999999999996 };
  double p8719_d29x1y1_ystatplus[] = { 0.0105, 0.011016, 0.008535, 0.007839, 0.010336000000000001, 0.014696, 0.025764000000000002, 0.043658999999999996 };
  int p8719_d29x1y1_numpoints = 8;
  TGraphAsymmErrors * p8719_d29x1y1 = new TGraphAsymmErrors(p8719_d29x1y1_numpoints, p8719_d29x1y1_xval, p8719_d29x1y1_yval, p8719_d29x1y1_xerrminus, p8719_d29x1y1_xerrplus, p8719_d29x1y1_yerrminus, p8719_d29x1y1_yerrplus);
  p8719_d29x1y1->SetName("/HepData/8719/d29x1y1");
  p8719_d29x1y1->SetTitle(" ");
  p8719_d29x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");
  p8719_d29x1y1->GetYaxis()->SetTitle("R_{AA}");

  // Plot: p8719_d30x1y1 30-40% 
  double p8719_d30x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d30x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d30x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d30x1y1_yval[] = { 0.542, 0.636, 0.676, 0.646, 0.676, 0.726, 0.768, 0.781, 0.786 };
  double p8719_d30x1y1_yerrminus[] = { 0.050102029699404395, 0.052299158234143696, 0.06069712111789158, 0.05365690155049955, 0.05949952161152222, 0.06826330577989904, 0.08064, 0.09072047052898259, 0.10676299077864015 };
  double p8719_d30x1y1_yerrplus[] = { 0.050102029699404395, 0.052299158234143696, 0.06069712111789158, 0.05365690155049955, 0.05949952161152222, 0.06826330577989904, 0.08064, 0.09072047052898259, 0.10676299077864015 };
  double p8719_d30x1y1_ystatminus[] = { 0.01084, 0.010812, 0.012168000000000002, 0.008398000000000001, 0.008788, 0.011616000000000001, 0.019200000000000002, 0.032021, 0.049518 };
  double p8719_d30x1y1_ystatplus[] = { 0.01084, 0.010812, 0.012168000000000002, 0.008398000000000001, 0.008788, 0.011616000000000001, 0.019200000000000002, 0.032021, 0.049518 };
  int p8719_d30x1y1_numpoints = 9;
  TGraphAsymmErrors * p8719_d30x1y1 = new TGraphAsymmErrors(p8719_d30x1y1_numpoints, p8719_d30x1y1_xval, p8719_d30x1y1_yval, p8719_d30x1y1_xerrminus, p8719_d30x1y1_xerrplus, p8719_d30x1y1_yerrminus, p8719_d30x1y1_yerrplus);
  p8719_d30x1y1->SetName("/HepData/8719/d30x1y1");
  p8719_d30x1y1->SetTitle(" ");
  p8719_d30x1y1->GetYaxis()->SetTitle("R_{AA}");
  p8719_d30x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");
  
  // Plot: p8719_d31x1y1 40-50% 
  double p8719_d31x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d31x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d31x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d31x1y1_yval[] = { 0.633, 0.707, 0.747, 0.716, 0.745, 0.776, 0.766, 0.725, 0.679 };
  double p8719_d31x1y1_yerrminus[] = { 0.058046479660699486, 0.05817196687408807, 0.06816144756385387, 0.060733515623583, 0.06801160452452214, 0.07559132267661416, 0.08290512519742071, 0.10722649509333035, 0.11789403304238939 };
  double p8719_d31x1y1_yerrplus[] = { 0.058046479660699486, 0.05817196687408807, 0.06816144756385387, 0.060733515623583, 0.06801160452452214, 0.07559132267661416, 0.08290512519742071, 0.10722649509333035, 0.11789403304238939 };
  double p8719_d31x1y1_ystatminus[] = { 0.01266, 0.013432999999999999, 0.013446, 0.01074, 0.011175000000000001, 0.01552, 0.024512000000000003, 0.041325, 0.061789000000000004 };
  double p8719_d31x1y1_ystatplus[] = { 0.01266, 0.013432999999999999, 0.013446, 0.01074, 0.011175000000000001, 0.01552, 0.024512000000000003, 0.041325, 0.061789000000000004 };
  int p8719_d31x1y1_numpoints = 9;
  TGraphAsymmErrors * p8719_d31x1y1 = new TGraphAsymmErrors(p8719_d31x1y1_numpoints, p8719_d31x1y1_xval, p8719_d31x1y1_yval, p8719_d31x1y1_xerrminus, p8719_d31x1y1_xerrplus, p8719_d31x1y1_yerrminus, p8719_d31x1y1_yerrplus);
  p8719_d31x1y1->SetName("/HepData/8719/d31x1y1");
  p8719_d31x1y1->SetTitle(" ");
  p8719_d31x1y1->GetYaxis()->SetTitle("R_{AA}");
  p8719_d31x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");

  // Plot: p8719_d32x1y1 50-60% 
  double p8719_d32x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d32x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d32x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d32x1y1_yval[] = { 0.696, 0.762, 0.797, 0.756, 0.818, 0.809, 0.759, 0.677, 0.623 };
  double p8719_d32x1y1_yerrminus[] = { 0.07291077549992181, 0.06435630234872106, 0.0687132337836024, 0.06296622294532203, 0.07997178997621598, 0.08381255617746067, 0.09553154035186497, 0.12123026731802583, 0.1683506930131266 };
  double p8719_d32x1y1_yerrplus[] = { 0.07291077549992181, 0.06435630234872106, 0.0687132337836024, 0.06296622294532203, 0.07997178997621598, 0.08381255617746067, 0.09553154035186497, 0.12123026731802583, 0.1683506930131266 };
  double p8719_d32x1y1_ystatminus[] = { 0.015312, 0.016764, 0.017534, 0.013608, 0.015541999999999999, 0.021843, 0.033396, 0.05348300000000001, 0.081613 };
  double p8719_d32x1y1_ystatplus[] = { 0.015312, 0.016764, 0.017534, 0.013608, 0.015541999999999999, 0.021843, 0.033396, 0.05348300000000001, 0.081613 };
  int p8719_d32x1y1_numpoints = 9;
  TGraphAsymmErrors * p8719_d32x1y1 = new TGraphAsymmErrors(p8719_d32x1y1_numpoints, p8719_d32x1y1_xval, p8719_d32x1y1_yval, p8719_d32x1y1_xerrminus, p8719_d32x1y1_xerrplus, p8719_d32x1y1_yerrminus, p8719_d32x1y1_yerrplus);
  p8719_d32x1y1->SetName("/HepData/8719/d32x1y1");
  p8719_d32x1y1->SetTitle(" ");
  p8719_d32x1y1->GetYaxis()->SetTitle("R_{AA}");
  p8719_d32x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");

  // Plot: p8719_d33x1y1 60-70% 
  double p8719_d33x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0 };
  double p8719_d33x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0 };
  double p8719_d33x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0 };
  double p8719_d33x1y1_yval[] = { 0.762, 0.773, 0.806, 0.791, 0.809, 0.825, 0.812, 0.768 };
  double p8719_d33x1y1_yerrminus[] = { 0.0842587336007372, 0.06169489225211436, 0.06025624950160772, 0.0649680231036777, 0.08105356425228935, 0.08651886355009525, 0.0994559132882505, 0.13203794072917072 };
  double p8719_d33x1y1_yerrplus[] = { 0.0842587336007372, 0.06169489225211436, 0.06025624950160772, 0.0649680231036777, 0.08105356425228935, 0.08651886355009525, 0.0994559132882505, 0.13203794072917072 };
  double p8719_d33x1y1_ystatminus[] = { 0.017526, 0.018552, 0.017732000000000005, 0.016611, 0.020225, 0.030525000000000004, 0.049532, 0.077568 };
  double p8719_d33x1y1_ystatplus[] = { 0.017526, 0.018552, 0.017732000000000005, 0.016611, 0.020225, 0.030525000000000004, 0.049532, 0.077568 };
  int p8719_d33x1y1_numpoints = 8;
  TGraphAsymmErrors * p8719_d33x1y1 = new TGraphAsymmErrors(p8719_d33x1y1_numpoints, p8719_d33x1y1_xval, p8719_d33x1y1_yval, p8719_d33x1y1_xerrminus, p8719_d33x1y1_xerrplus, p8719_d33x1y1_yerrminus, p8719_d33x1y1_yerrplus);
  p8719_d33x1y1->SetName("/HepData/8719/d33x1y1");
  p8719_d33x1y1->SetTitle(" ");
  p8719_d33x1y1->GetYaxis()->SetTitle("R_{AA}");
  p8719_d33x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");

  // Plot: p8719_d34x1y1 70-80% 
  double p8719_d34x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5 };
  double p8719_d34x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5 };
  double p8719_d34x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5 };
  double p8719_d34x1y1_yval[] = { 0.823, 0.795, 0.812, 0.809, 0.814, 0.819, 0.801 };
  double p8719_d34x1y1_yerrminus[] = { 0.10341675352668928, 0.07952384642357285, 0.07091415012534523, 0.06899304539589479, 0.08790069371739907, 0.09622465819113102, 0.11634094161558088 };
  double p8719_d34x1y1_yerrplus[] = { 0.10341675352668928, 0.07952384642357285, 0.07091415012534523, 0.06899304539589479, 0.08790069371739907, 0.09622465819113102, 0.11634094161558088 };
  double p8719_d34x1y1_ystatminus[] = { 0.024689999999999997, 0.026235, 0.023548, 0.02427, 0.030931999999999998, 0.047501999999999996, 0.078498 };
  double p8719_d34x1y1_ystatplus[] = { 0.024689999999999997, 0.026235, 0.023548, 0.02427, 0.030931999999999998, 0.047501999999999996, 0.078498 };
  int p8719_d34x1y1_numpoints = 7;
  TGraphAsymmErrors * p8719_d34x1y1 = new TGraphAsymmErrors(p8719_d34x1y1_numpoints, p8719_d34x1y1_xval, p8719_d34x1y1_yval, p8719_d34x1y1_xerrminus, p8719_d34x1y1_xerrplus, p8719_d34x1y1_yerrminus, p8719_d34x1y1_yerrplus);
  p8719_d34x1y1->SetName("/HepData/8719/d34x1y1");
  p8719_d34x1y1->SetTitle(" ");
  p8719_d34x1y1->GetYaxis()->SetTitle("R_{AA}");
  p8719_d34x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");

  TCanvas * cATLAS = new TCanvas("cATLAS","",1200,1000);
  makeMultiPanelCanvas(cATLAS,3,2,0.0,0.0,0.2,0.15,0.07);
  cATLAS->cd(6);

  // centrality bin - 0-10
  p8719_d27x1y1->SetMaximum(2);
  p8719_d27x1y1->SetMinimum(0);
  p8719_d27x1y1->SetMarkerStyle(33);
  p8719_d27x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d27x1y1->Draw("ap");
  systematics_R4.calcTotalSys(0);
  systematics_R4.Draw(RAA_R4_Bayes[0],0,4);
  systematics_R4.calcTotalSys(1);
  systematics_R4.Draw(RAA_R4_Bayes[1],0,4);
  p8719_d27x1y1->Draw("psame");
  RAA_R4_Bayes[0]->SetMarkerStyle(24);
  RAA_R4_Bayes[0]->Draw("same E1");
  RAA_R4_Bayes[1]->SetMarkerStyle(24);
  RAA_R4_Bayes[1]->SetMarkerColor(kRed);
  RAA_R4_Bayes[1]->Draw("same E1");
  //drawText("0-10%", 0.2,0.7,16);
  
  TLegend *leg1 = getLegend(0.20,0.65,0.4,0.85);
  leg1->AddEntry(p8719_d27x1y1,"ATLAS 0-10%","p");
  leg1->AddEntry(RAA_R4_Bayes[0],"CMS 0-5%","p");
  leg1->AddEntry(RAA_R4_Bayes[1],"CMS 5-10%","p");
  leg1->SetTextSize(0.04);
  leg1->Draw();
  
  // TLegend * comp0 = myLegend(0.2,0.2,0.4,0.4);
  // comp0->AddEntry(p8719_d27x1y1,"ATLAS 0-10%","pl");
  // comp0->AddEntry(RAA_R4_Bayes[0],"CMS 0-5%","pl");
  // comp0->Draw();

  // centrality bin 10-30
  cATLAS->cd(5);
  p8719_d28x1y1->SetMaximum(2);
  p8719_d28x1y1->SetMinimum(0);
  p8719_d28x1y1->SetMarkerStyle(33);
  p8719_d28x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d28x1y1->Draw("AP");
  RAA_R4_Bayes[2]->Draw("same E1");
  RAA_R4_Bayes[2]->SetMarkerStyle(24);
  systematics_R4.calcTotalSys(2);
  systematics_R4.Draw(RAA_R4_Bayes[2],2,4);
  p8719_d28x1y1->Draw("psame");
  p8719_d29x1y1->SetMarkerStyle(34);
  p8719_d29x1y1->Draw("psame");
  RAA_R4_Bayes[2]->Draw("same E1");
  //drawText("10-30%", 0.2,0.7,16);
  TLegend *leg2 = getLegend(0.20,0.65,0.4,0.85);
  leg2->AddEntry(p8719_d28x1y1,"ATLAS 10-20%","p");
  leg2->AddEntry(p8719_d29x1y1,"ATLAS 20-30%","p");
  leg2->AddEntry(RAA_R4_Bayes[2],"CMS 10-30%","p");
  leg2->SetTextSize(0.04);
  leg2->Draw();
  // TLegend * comp2 = myLegend(0.2,0.2,0.4,0.4);
  // comp2->AddEntry(p8719_d27x1y1,"ATLAS 10-20%","pl");
  // comp2->AddEntry(RAA_R4_Bayes[2],"CMS 10-30%","pl");
  // comp2->Draw();
  
  // centrality bin 30-50 %
  cATLAS->cd(4);
  p8719_d30x1y1->SetMaximum(2);
  p8719_d30x1y1->SetMinimum(0);
  p8719_d30x1y1->SetMarkerStyle(33);
  p8719_d30x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d30x1y1->Draw("AP");
  RAA_R4_Bayes[3]->Draw("same E1");
  systematics_R4.calcTotalSys(3);
  systematics_R4.Draw(RAA_R4_Bayes[3],3,4);
  RAA_R4_Bayes[3]->SetMarkerStyle(24);
  p8719_d30x1y1->Draw("psame");
  p8719_d31x1y1->SetMarkerStyle(34);
  p8719_d31x1y1->Draw("psame");
  RAA_R4_Bayes[3]->Draw("same E1");
  //drawText("30-50%", 0.2,0.7,16);
  TLegend *leg3 = getLegend(0.20,0.65,0.4,0.85);
  leg3->AddEntry(p8719_d30x1y1,"ATLAS 30-40%","p");
  leg3->AddEntry(p8719_d31x1y1,"ATLAS 40-50%","p");
  leg3->AddEntry(RAA_R4_Bayes[3],"CMS 30-50%","p");
  leg3->SetTextSize(0.04);
  leg3->Draw();

  
  // centrality bin 50-70 %
  cATLAS->cd(3);
  p8719_d32x1y1->SetMaximum(2);
  p8719_d32x1y1->SetMinimum(0);
  p8719_d32x1y1->SetMarkerStyle(33);
  p8719_d32x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d32x1y1->Draw("AP");
  RAA_R4_Bayes[4]->Draw("same E1");
  systematics_R4.calcTotalSys(4);
  systematics_R4.Draw(RAA_R4_Bayes[4],4,4);
  RAA_R4_Bayes[4]->SetMarkerStyle(24);
  p8719_d32x1y1->Draw("psame");
  p8719_d33x1y1->SetMarkerStyle(34);
  p8719_d33x1y1->Draw("psame");
  RAA_R4_Bayes[4]->Draw("same E1");
  TLegend *leg4 = getLegend(0.20,0.65,0.4,0.85);
  leg4->AddEntry(p8719_d32x1y1,"ATLAS 50-60%","p");
  leg4->AddEntry(p8719_d33x1y1,"ATLAS 60-70%","p");
  leg4->AddEntry(RAA_R4_Bayes[4],"CMS 50-70%","p");
  leg4->SetTextSize(0.04);
  leg4->Draw();

  
  // centrality bin 50-70 %
  cATLAS->cd(2);
  p8719_d34x1y1->SetMaximum(2);
  p8719_d34x1y1->SetMinimum(0);
  p8719_d34x1y1->SetMarkerStyle(33);
  p8719_d34x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d34x1y1->Draw("AP");
  RAA_R4_Bayes[5]->Draw("same E1");
  systematics_R4.calcTotalSys(5);
  systematics_R4.Draw(RAA_R4_Bayes[5],5,4);
  RAA_R4_Bayes[5]->SetMarkerStyle(24);
  p8719_d34x1y1->Draw("psame");
  RAA_R4_Bayes[5]->Draw("same E1");
  TLegend *leg5 = getLegend(0.20,0.65,0.4,0.85);
  leg5->AddEntry(p8719_d32x1y1,"ATLAS 70-80%","p");
  leg5->AddEntry(RAA_R4_Bayes[5],"CMS 70-90%","p");
  leg5->SetTextSize(0.04);
  leg5->Draw();
  
  
  cATLAS->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_RAA_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cATLAS->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_RAA_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cATLAS->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_RAA_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");


  // Lets compare the PbPb spectra here:
  // ATLAS pp spectra: R=0.4
  // Plot: p8719_d2x1y1
  double p8719_d2x1y1_xval[] = { 35.0, 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 
				 283.5, 357.0, 449.5 };
  double p8719_d2x1y1_xerrminus[] = { 4.0, 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 
				      32.5, 41.0, 51.5 };
  double p8719_d2x1y1_xerrplus[] = { 4.0, 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 
				     32.5, 41.0, 51.5 };
  double p8719_d2x1y1_yval[] = { 180.0, 55.7, 16.9, 4.85, 1.42, 0.364, 0.0882, 0.0197, 0.00406, 
				 7.35E-4, 1.14E-4, 1.41E-5 };
  double p8719_d2x1y1_yerrminus[] = { 29.95905205442923, 7.828376905208385, 2.6254355333925075, 0.627695686618922, 0.18780540141327137, 0.047724273069372145, 0.011038045881404914, 0.002292152355320213, 4.82252125345239E-4, 
				      8.981748117710716E-5, 1.442494214893079E-5, 1.9885499490835025E-6 };
  double p8719_d2x1y1_yerrplus[] = { 29.95905205442923, 7.828376905208385, 2.6254355333925075, 0.627695686618922, 0.18780540141327137, 0.047724273069372145, 0.011038045881404914, 0.002292152355320213, 4.82252125345239E-4, 
				     8.981748117710716E-5, 1.442494214893079E-5, 1.9885499490835025E-6 };
  double p8719_d2x1y1_ystatminus[] = { 1.9800000000000002, 0.6684, 0.16899999999999998, 0.048499999999999995, 0.01136, 0.00182, 5.292E-4, 1.773E-4, 5.2780000000000006E-5, 
				       1.6904999999999998E-5, 4.674E-6, 8.883E-7 };
  double p8719_d2x1y1_ystatplus[] = { 1.9800000000000002, 0.6684, 0.16899999999999998, 0.048499999999999995, 0.01136, 0.00182, 5.292E-4, 1.773E-4, 5.2780000000000006E-5, 
				      1.6904999999999998E-5, 4.674E-6, 8.883E-7 };
  int p8719_d2x1y1_numpoints = 12;
  TGraphAsymmErrors * p8719_d2x1y1 = new TGraphAsymmErrors(p8719_d2x1y1_numpoints, p8719_d2x1y1_xval, p8719_d2x1y1_yval, p8719_d2x1y1_xerrminus, p8719_d2x1y1_xerrplus, p8719_d2x1y1_yerrminus, p8719_d2x1y1_yerrplus);
  p8719_d2x1y1->SetName("/HepData/8719/d2x1y1");
  p8719_d2x1y1->SetTitle(" ");


  // i took care of all the scaling in the RAA_analyze marco 
  //uPP_R4_Bayes->Scale(1./4/5.3/1e3);// to get cross sections in nano barns
  //uPP_R4_Bayes->Scale(1./4/5.429/1e3/0.82698/0.25249);
  //uPP_R4_Bayes->Scale(1./(4 * 5.429 * 1e3 * 0.82698));
  // the 5.429 is from the ntuple and the other factor of 0.82698 is due to event losses from the ntuple making with jet matching
    
  TCanvas * cATLAS_pp = new TCanvas("cATLAS_pp","",1200,1000);
  cATLAS_pp->SetLogy();
  p8719_d2x1y1->GetXaxis()->SetTitle("ak R=0.4 Jet p_{T} (GeV/c)");
  p8719_d2x1y1->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T} d#eta} nb");  
  p8719_d2x1y1->SetMarkerStyle(33);
  p8719_d2x1y1->SetMarkerColor(kBlack);
  p8719_d2x1y1->GetXaxis()->SetLimits(50,450);
  p8719_d2x1y1->Draw("ap");
  uPP_R4_Bayes->SetMarkerStyle(24);
  uPP_R4_Bayes->SetMarkerColor(kBlue);
  uPP_R4_Bayes->Draw("same E1");
  systematics_R4.calcTotalSys(6);
  systematics_R4.Draw(uPP_R4_Bayes,6,2);
  
  TLegend *legpp_1 = getLegend(0.40,0.65,0.6,0.85);
  legpp_1->AddEntry(p8719_d2x1y1,"ATLAS pp","p");
  legpp_1->AddEntry(uPP_R4_Bayes,"CMS pp","p");
  legpp_1->SetTextSize(0.04);
  legpp_1->Draw();
  
  cATLAS_pp->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_pp_spectra_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cATLAS_pp->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_pp_spectra_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cATLAS_pp->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_pp_spectra_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");
  
  // TAA for atlas plots:// which is in inverse milli barns 
  float_t TAA[8] = {23.45, 14.43, 8.73, 5.04, 2.7, 1.33, 0.59, 0.24};
  
  // PbPb spectra for R=0.4, 0-10%
  // Plot: p8719_d7x1y1
  double p8719_d7x1y1_xval[] = { 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5, 357.0 };
  double p8719_d7x1y1_xerrminus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d7x1y1_xerrplus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d7x1y1_yval[] = { 1.87E-4, 5.58E-5, 1.51E-5, 4.08E-6, 1.06E-6, 2.61E-7, 5.7E-8, 1.03E-8, 1.5E-9 };
  double p8719_d7x1y1_yerrminus[] = { 3.043220327219178E-5, 1.1792957400075692E-5, 3.363912646012081E-6, 7.439709671754671E-7, 1.809684204495359E-7, 4.3891478443998676E-8, 1.0058844764683469E-8, 2.0244937564734545E-9, 3.6489484786716294E-10 };
  double p8719_d7x1y1_yerrplus[] = { 3.043220327219178E-5, 1.1792957400075692E-5, 3.363912646012081E-6, 7.439709671754671E-7, 1.809684204495359E-7, 4.3891478443998676E-8, 1.0058844764683469E-8, 2.0244937564734545E-9, 3.6489484786716294E-10 };
  double p8719_d7x1y1_ystatminus[] = { 4.114E-6, 1.116E-6, 2.416E-7, 6.12E-8, 1.59E-8, 4.698000000000001E-9, 1.482E-9, 4.841E-10, 1.29E-10 };
  double p8719_d7x1y1_ystatplus[] = { 4.114E-6, 1.116E-6, 2.416E-7, 6.12E-8, 1.59E-8, 4.698000000000001E-9, 1.482E-9, 4.841E-10, 1.29E-10 };
  int p8719_d7x1y1_numpoints = 9;
  for(int i = 0; i<p8719_d7x1y1_numpoints; ++i) p8719_d7x1y1_yval[i] = p8719_d7x1y1_yval[i] * (1./TAA[0]) * 1e6;
  TGraphAsymmErrors * p8719_d7x1y1 = new TGraphAsymmErrors(p8719_d7x1y1_numpoints, p8719_d7x1y1_xval, p8719_d7x1y1_yval, p8719_d7x1y1_xerrminus, p8719_d7x1y1_xerrplus, p8719_d7x1y1_yerrminus, p8719_d7x1y1_yerrplus);
  p8719_d7x1y1->SetName("/HepData/8719/d7x1y1");
  p8719_d7x1y1->SetTitle("  ");

  // Plot: p8719_d8x1y1, 10-20% 

  double p8719_d8x1y1_xval[] = { 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5, 357.0 };
  double p8719_d8x1y1_xerrminus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d8x1y1_xerrplus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d8x1y1_yval[] = { 1.28E-4, 3.7E-5, 1.02E-5, 2.85E-6, 7.38E-7, 1.86E-7, 3.82E-8, 6.97E-9, 1.07E-9 };
  double p8719_d8x1y1_yerrminus[] = { 2.365369349594266E-5, 7.71625219909251E-6, 1.7929093340155267E-6, 5.048523818899936E-7, 1.2186165641414859E-7, 3.001466041786913E-8, 6.654369827414162E-9, 1.3162416403533207E-9, 2.3045773169932924E-10 };
  double p8719_d8x1y1_yerrplus[] = { 2.365369349594266E-5, 7.71625219909251E-6, 1.7929093340155267E-6, 5.048523818899936E-7, 1.2186165641414859E-7, 3.001466041786913E-8, 6.654369827414162E-9, 1.3162416403533207E-9, 2.3045773169932924E-10 };
  double p8719_d8x1y1_ystatminus[] = { 2.8160000000000002E-6, 8.88E-7, 1.734E-7, 4.2749999999999996E-8, 1.1808E-8, 3.72E-9, 1.1842E-9, 3.9729000000000003E-10, 1.0058000000000002E-10 };
  double p8719_d8x1y1_ystatplus[] = { 2.8160000000000002E-6, 8.88E-7, 1.734E-7, 4.2749999999999996E-8, 1.1808E-8, 3.72E-9, 1.1842E-9, 3.9729000000000003E-10, 1.0058000000000002E-10 };
  int p8719_d8x1y1_numpoints = 9;
  for(int i = 0; i<p8719_d8x1y1_numpoints; ++i) p8719_d8x1y1_yval[i] = p8719_d8x1y1_yval[i] * (1./TAA[1]) * 1e6;
  TGraphAsymmErrors * p8719_d8x1y1 = new TGraphAsymmErrors(p8719_d8x1y1_numpoints, p8719_d8x1y1_xval, p8719_d8x1y1_yval, p8719_d8x1y1_xerrminus, p8719_d8x1y1_xerrplus, p8719_d8x1y1_yerrminus, p8719_d8x1y1_yerrplus);
  p8719_d8x1y1->SetName("/HepData/8719/d8x1y1");
  p8719_d8x1y1->SetTitle("  ");

  // Plot: p8719_d9x1y1, 20-30% 
  double p8719_d9x1y1_xval[] = { 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5, 357.0 };
  double p8719_d9x1y1_xerrminus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d9x1y1_xerrplus[] = { 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 41.0 };
  double p8719_d9x1y1_yval[] = { 7.77E-5, 2.59E-5, 7.05E-6, 1.92E-6, 4.98E-7, 1.15E-7, 2.41E-8, 4.45E-9, 7.23E-10 };
  double p8719_d9x1y1_yerrminus[] = { 1.3550577037897685E-5, 4.4971380554748375E-6, 1.258834702214711E-6, 3.4389445822810234E-7, 8.247878345368583E-8, 1.9121149285542434E-8, 4.183832188317309E-9, 8.232018904861674E-10, 1.4669596951518472E-10 };
  double p8719_d9x1y1_yerrplus[] = { 1.3550577037897685E-5, 4.4971380554748375E-6, 1.258834702214711E-6, 3.4389445822810234E-7, 8.247878345368583E-8, 1.9121149285542434E-8, 4.183832188317309E-9, 8.232018904861674E-10, 1.4669596951518472E-10 };
  double p8719_d9x1y1_ystatminus[] = { 1.3209E-6, 3.8849999999999996E-7, 9.165E-8, 2.304E-8, 7.47E-9, 2.415E-9, 8.435E-10, 2.6255000000000004E-10, 6.0732E-11 };
  double p8719_d9x1y1_ystatplus[] = { 1.3209E-6, 3.8849999999999996E-7, 9.165E-8, 2.304E-8, 7.47E-9, 2.415E-9, 8.435E-10, 2.6255000000000004E-10, 6.0732E-11 };
  int p8719_d9x1y1_numpoints = 9;
  for(int i = 0; i<p8719_d9x1y1_numpoints; ++i) p8719_d9x1y1_yval[i] = p8719_d9x1y1_yval[i] * (1./TAA[2]) * 1e6;
  TGraphAsymmErrors * p8719_d9x1y1 = new TGraphAsymmErrors(p8719_d9x1y1_numpoints, p8719_d9x1y1_xval, p8719_d9x1y1_yval, p8719_d9x1y1_xerrminus, p8719_d9x1y1_xerrplus, p8719_d9x1y1_yerrminus, p8719_d9x1y1_yerrplus);
  p8719_d9x1y1->SetName("/HepData/8719/d9x1y1");
  p8719_d9x1y1->SetTitle("  ");

  // Plot: p8719_d10x1y1, 30-40% 
  double p8719_d10x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5, 
				  357.0 };
  double p8719_d10x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 
				       41.0 };
  double p8719_d10x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5, 
				      41.0 };
  double p8719_d10x1y1_yval[] = { 1.52E-4, 5.42E-5, 1.65E-5, 4.62E-6, 1.24E-6, 3.22E-7, 7.64E-8, 1.6E-8, 2.91E-9, 
				  4.65E-10 };
  double p8719_d10x1y1_yerrminus[] = { 2.361534145423267E-5, 9.526565496547013E-6, 3.137951242451036E-6, 8.184965991377117E-7, 2.0816124903545332E-7, 5.382410839020002E-8, 1.3172964125055529E-8, 2.8537021568481882E-9, 5.391437661700263E-10, 
				       8.978956147570829E-11 };
  double p8719_d10x1y1_yerrplus[] = { 2.361534145423267E-5, 9.526565496547013E-6, 3.137951242451036E-6, 8.184965991377117E-7, 2.0816124903545332E-7, 5.382410839020002E-8, 1.3172964125055529E-8, 2.8537021568481882E-9, 5.391437661700263E-10, 
				      8.978956147570829E-11 };
  double p8719_d10x1y1_ystatminus[] = { 2.28E-6, 7.588E-7, 2.31E-7, 5.082E-8, 1.488E-8, 4.83E-9, 1.7571999999999997E-9, 6.24E-10, 1.6878E-10, 
					3.5805E-11 };
  double p8719_d10x1y1_ystatplus[] = { 2.28E-6, 7.588E-7, 2.31E-7, 5.082E-8, 1.488E-8, 4.83E-9, 1.7571999999999997E-9, 6.24E-10, 1.6878E-10, 
				       3.5805E-11 };
  int p8719_d10x1y1_numpoints = 10;
  for(int i = 0; i<p8719_d10x1y1_numpoints; ++i) p8719_d10x1y1_yval[i] = p8719_d10x1y1_yval[i] * (1./TAA[3]) * 1e6;
  TGraphAsymmErrors * p8719_d10x1y1 = new TGraphAsymmErrors(p8719_d10x1y1_numpoints, p8719_d10x1y1_xval, p8719_d10x1y1_yval, p8719_d10x1y1_xerrminus, p8719_d10x1y1_xerrplus, p8719_d10x1y1_yerrminus, p8719_d10x1y1_yerrplus);
  p8719_d10x1y1->SetName("/HepData/8719/d10x1y1");
  p8719_d10x1y1->SetTitle("  ");

  // Plot: p8719_d11x1y1, 40-50% 
  double p8719_d11x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d11x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d11x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d11x1y1_yval[] = { 9.53E-5, 3.24E-5, 9.8E-6, 2.75E-6, 7.33E-7, 1.85E-7, 4.09E-8, 7.96E-9, 1.35E-9 };
  double p8719_d11x1y1_yerrminus[] = { 1.6077480010561357E-5, 6.086803732666266E-6, 1.9425271169278434E-6, 5.289802099322809E-7, 1.3585240943759519E-7, 3.3748178469363346E-8, 7.763785335131312E-9, 1.5941278779320057E-9, 3.0047130478633063E-10 };
  double p8719_d11x1y1_yerrplus[] = { 1.6077480010561357E-5, 6.086803732666266E-6, 1.9425271169278434E-6, 5.289802099322809E-7, 1.3585240943759519E-7, 3.3748178469363346E-8, 7.763785335131312E-9, 1.5941278779320057E-9, 3.0047130478633063E-10 };
  double p8719_d11x1y1_ystatminus[] = { 1.5248000000000002E-6, 5.184E-7, 1.4699999999999998E-7, 3.3E-8, 1.0262E-8, 3.515E-9, 1.2679000000000001E-9, 4.3780000000000003E-10, 1.1880000000000002E-10 };
  double p8719_d11x1y1_ystatplus[] = { 1.5248000000000002E-6, 5.184E-7, 1.4699999999999998E-7, 3.3E-8, 1.0262E-8, 3.515E-9, 1.2679000000000001E-9, 4.3780000000000003E-10, 1.1880000000000002E-10 };
  int p8719_d11x1y1_numpoints = 9;
  for(int i = 0; i<p8719_d11x1y1_numpoints; ++i) p8719_d11x1y1_yval[i] = p8719_d11x1y1_yval[i] * (1./TAA[4]) * 1e6;
  TGraphAsymmErrors * p8719_d11x1y1 = new TGraphAsymmErrors(p8719_d11x1y1_numpoints, p8719_d11x1y1_xval, p8719_d11x1y1_yval, p8719_d11x1y1_xerrminus, p8719_d11x1y1_xerrplus, p8719_d11x1y1_yerrminus, p8719_d11x1y1_yerrplus);
  p8719_d11x1y1->SetName("/HepData/8719/d11x1y1");
  p8719_d11x1y1->SetTitle("  ");

  // Plot: p8719_d12x1y1, 50-60% 
  double p8719_d12x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d12x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d12x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d12x1y1_yval[] = { 5.16E-5, 1.72E-5, 5.13E-6, 1.43E-6, 3.96E-7, 9.49E-8, 1.99E-8, 3.66E-9, 6.08E-10 };
  double p8719_d12x1y1_yerrminus[] = { 1.0819524618022736E-5, 3.5765248384430385E-6, 1.081213676800289E-6, 2.96365559402573E-7, 8.163557048248024E-8, 1.968254369130169E-8, 4.313941809992341E-9, 8.801672222935821E-10, 1.7106963559907412E-10 };
  double p8719_d12x1y1_yerrplus[] = { 1.0819524618022736E-5, 3.5765248384430385E-6, 1.081213676800289E-6, 2.96365559402573E-7, 8.163557048248024E-8, 1.968254369130169E-8, 4.313941809992341E-9, 8.801672222935821E-10, 1.7106963559907412E-10 };
  double p8719_d12x1y1_ystatminus[] = { 9.288E-7, 3.268E-7, 1.026E-7, 2.2880000000000004E-8, 7.524E-9, 2.4674E-9, 8.557E-10, 2.8548E-10, 7.8432E-11 };
  double p8719_d12x1y1_ystatplus[] = { 9.288E-7, 3.268E-7, 1.026E-7, 2.2880000000000004E-8, 7.524E-9, 2.4674E-9, 8.557E-10, 2.8548E-10, 7.8432E-11 };
  int p8719_d12x1y1_numpoints = 9;
  for(int i = 0; i<p8719_d12x1y1_numpoints; ++i) p8719_d12x1y1_yval[i] = p8719_d12x1y1_yval[i] * (1./TAA[5]) * 1e6;  
  TGraphAsymmErrors * p8719_d12x1y1 = new TGraphAsymmErrors(p8719_d12x1y1_numpoints, p8719_d12x1y1_xval, p8719_d12x1y1_yval, p8719_d12x1y1_xerrminus, p8719_d12x1y1_xerrplus, p8719_d12x1y1_yerrminus, p8719_d12x1y1_yerrplus);
  p8719_d12x1y1->SetName("/HepData/8719/d12x1y1");
  p8719_d12x1y1->SetTitle("  ");

  
  // Plot: p8719_d13x1y1, 60-70% 
  double p8719_d13x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0, 283.5 };
  double p8719_d13x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d13x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0, 32.5 };
  double p8719_d13x1y1_yval[] = { 2.52E-5, 7.75E-6, 2.31E-6, 6.66E-7, 1.74E-7, 4.31E-8, 9.49E-9, 1.85E-9, 3.25E-10 };
  double p8719_d13x1y1_yerrminus[] = { 5.791889847018847E-6, 1.9076652942798953E-6, 5.322733258392721E-7, 1.5522801339964382E-7, 4.1300357916124654E-8, 1.0076462920588751E-8, 2.284016482536849E-9, 4.806084945150263E-10, 9.45839342859029E-11 };
  double p8719_d13x1y1_yerrplus[] = { 5.791889847018847E-6, 1.9076652942798953E-6, 5.322733258392721E-7, 1.5522801339964382E-7, 4.1300357916124654E-8, 1.0076462920588751E-8, 2.284016482536849E-9, 4.806084945150263E-10, 9.45839342859029E-11 };
  double p8719_d13x1y1_ystatminus[] = { 5.04E-7, 1.6275E-7, 4.389E-8, 1.332E-8, 4.35E-9, 1.5947E-9, 5.693999999999999E-10, 1.85E-10, 4.81E-11 };
  double p8719_d13x1y1_ystatplus[] = { 5.04E-7, 1.6275E-7, 4.389E-8, 1.332E-8, 4.35E-9, 1.5947E-9, 5.693999999999999E-10, 1.85E-10, 4.81E-11 };
  int p8719_d13x1y1_numpoints = 9;
  for(int i = 0; i<p8719_d13x1y1_numpoints; ++i) p8719_d13x1y1_yval[i] = p8719_d13x1y1_yval[i] * (1./TAA[6]) * 1e6;
  TGraphAsymmErrors * p8719_d13x1y1 = new TGraphAsymmErrors(p8719_d13x1y1_numpoints, p8719_d13x1y1_xval, p8719_d13x1y1_yval, p8719_d13x1y1_xerrminus, p8719_d13x1y1_xerrplus, p8719_d13x1y1_yerrminus, p8719_d13x1y1_yerrplus);
  p8719_d13x1y1->SetName("/HepData/8719/d13x1y1");
  p8719_d13x1y1->SetTitle("  ");

  // Plot: p8719_d14x1y1, 70-80% 
  double p8719_d14x1y1_xval[] = { 44.5, 56.5, 71.0, 89.5, 112.5, 141.5, 178.5, 225.0 };
  double p8719_d14x1y1_xerrminus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0 };
  double p8719_d14x1y1_xerrplus[] = { 5.5, 6.5, 8.0, 10.5, 12.5, 16.5, 20.5, 26.0 };
  double p8719_d14x1y1_yval[] = { 1.08E-5, 3.16E-6, 9.25E-7, 2.7E-7, 6.97E-8, 1.7E-8, 3.72E-9, 7.31E-10 };
  double p8719_d14x1y1_yerrminus[] = { 2.994659085772536E-6, 9.421146647834327E-7, 2.559744262128543E-7, 7.343751834042325E-8, 1.943530477687448E-8, 4.7535295307802604E-9, 1.0803830454056562E-9, 2.3182152646162954E-10 };
  double p8719_d14x1y1_yerrplus[] = { 2.994659085772536E-6, 9.421146647834327E-7, 2.559744262128543E-7, 7.343751834042325E-8, 1.943530477687448E-8, 4.7535295307802604E-9, 1.0803830454056562E-9, 2.3182152646162954E-10 };
  double p8719_d14x1y1_ystatminus[] = { 2.9160000000000004E-7, 9.796E-8, 2.4975E-8, 7.83E-9, 2.6485999999999997E-9, 9.69E-10, 3.6083999999999997E-10, 1.1769100000000002E-10 };
  double p8719_d14x1y1_ystatplus[] = { 2.9160000000000004E-7, 9.796E-8, 2.4975E-8, 7.83E-9, 2.6485999999999997E-9, 9.69E-10, 3.6083999999999997E-10, 1.1769100000000002E-10 };
  int p8719_d14x1y1_numpoints = 8;
  for(int i = 0; i<p8719_d14x1y1_numpoints; ++i) p8719_d14x1y1_yval[i] = p8719_d14x1y1_yval[i] * (1./TAA[7]) * 1e6;

  TGraphAsymmErrors * p8719_d14x1y1 = new TGraphAsymmErrors(p8719_d14x1y1_numpoints, p8719_d14x1y1_xval, p8719_d14x1y1_yval, p8719_d14x1y1_xerrminus, p8719_d14x1y1_xerrplus, p8719_d14x1y1_yerrminus, p8719_d14x1y1_yerrplus);
  p8719_d14x1y1->SetName("/HepData/8719/d14x1y1");
  p8719_d14x1y1->SetTitle(" ");


  /*
    uPbPb_Bayes[i]->Scale(1./deltaEta);// delta eta
    //uPbPb_Bayes[i]->Scale(1./145.156/1e6);// Jet 80 luminosity
    //uPbPb_Bayes[i]->Scale(1./1.1153/1e6);// equivalent no of minbias events 
    uPbPb_Bayes[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    uPbPb_Bayes[i]->Scale(1./145.156);
    //uPbPb_Bayes[i]->Scale(1./161.939);
    uPbPb_Bayes[i]->Scale(1./(7.65*1e6));
    uPbPb_Bayes[i]->Scale(64.*1e9/(ncoll[i]*1e3));
    uPbPb_Bayes[i] = (TH1F*)uPbPb_Bayes[i]->Rebin(nbins_pt,Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i),boundaries_pt);
    divideBinWidth(uPbPb_Bayes[i]);
    uPbPb_Bayes[i]->Write();
  */
  // scale our PbPb spectra to 1/(Nevt = 2.1966046e+07) to get yield.

  // April 22 - update from Nevt calculation:
  // High pT data:
  // total number of events = 21965900
  // Jet80 number of events =  1115080
  // MinBias data:
  // total number of events = 58867850
  // Jet80 number of events =    26423

  // ratio of jet80/total# in minbias = 0.00045
  // There total number of events = 21965900 * 1/0.00045 = 48813111111 = 4.88e+10

  // for(int i = 0; i<nbins_cent;++i){
  //   uPbPb_R4_Bayes[i]->Scale(64/(0.787 * 7.65 * ncoll[i] * 4 * 0.967 * 160.521 * (0.025*(boundaries_cent[i+1] - boundaries_cent[i]))));
  //   // for 0.787 is the difference between the event calculation versus lumi approach. 
  //   //uPbPb_R4_Bayes[i]->Scale(1./(9.44 4 * (0.025*(boundaries_cent[i+1] - boundaries_cent[i])))); 
  // }

  TCanvas * cATLAS_pbpb = new TCanvas("cATLAS_pbpb","",1200,1000);
  makeMultiPanelCanvas(cATLAS_pbpb,3,2,0.0,0.0,0.2,0.15,0.07);

  // centrality bin 0-10% 
  cATLAS_pbpb->cd(6);
  cATLAS_pbpb->cd(6)->SetLogy();
  p8719_d7x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d7x1y1->SetMaximum(1e2);
  p8719_d7x1y1->SetMinimum(1e-4);
  p8719_d7x1y1->SetMarkerStyle(33);
  //  p8719_d7x1y1->Scale(1./TAA[0]);
  p8719_d7x1y1->Draw("ap");
  uPbPb_R4_Bayes[0]->SetMarkerStyle(24);
  uPbPb_R4_Bayes[0]->SetMarkerColor(kBlue);
  uPbPb_R4_Bayes[0]->Draw("same");
  uPbPb_R4_Bayes[1]->SetMarkerStyle(24);
  uPbPb_R4_Bayes[1]->SetMarkerColor(kRed);
  uPbPb_R4_Bayes[1]->Draw("same");
  systematics_R4.calcTotalSys(0);
  systematics_R4.Draw(uPbPb_R4_Bayes[0],0,2);
  systematics_R4.calcTotalSys(1);
  systematics_R4.Draw(uPbPb_R4_Bayes[1],1,3);
  
  TLegend *legpbpb_1 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_1->AddEntry(p8719_d7x1y1,"ATLAS 0-10%","p");
  legpbpb_1->AddEntry(uPbPb_R4_Bayes[0],"CMS 0-5%","p");
  legpbpb_1->AddEntry(uPbPb_R4_Bayes[1],"CMS 5-10%","p");
  legpbpb_1->SetTextSize(0.04);
  legpbpb_1->Draw();

  // centrality bin 10-30% 
  cATLAS_pbpb->cd(5);
  cATLAS_pbpb->cd(5)->SetLogy();
  p8719_d8x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d8x1y1->SetMarkerStyle(33);
  p8719_d8x1y1->SetMaximum(1e2);
  p8719_d8x1y1->SetMinimum(1e-4);
  // p8719_d8x1y1->Scale(1./TAA[1]);
  p8719_d8x1y1->Draw("ap");
  p8719_d9x1y1->SetMarkerStyle(34);
  //p8719_d9x1y1->Scale(1./TAA[2]);
  p8719_d9x1y1->Draw("psame");
  uPbPb_R4_Bayes[2]->SetMarkerStyle(24);
  uPbPb_R4_Bayes[2]->SetMarkerColor(kBlue);
  uPbPb_R4_Bayes[2]->Draw("same");
  systematics_R4.calcTotalSys(2);
  systematics_R4.Draw(uPbPb_R4_Bayes[2],2,2);
  
  TLegend *legpbpb_2 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_2->AddEntry(p8719_d8x1y1,"ATLAS 10-20%","p");
  legpbpb_2->AddEntry(p8719_d9x1y1,"ATLAS 20-30%","p");
  legpbpb_2->AddEntry(uPbPb_R4_Bayes[2],"CMS 10-30%","p");
  legpbpb_2->SetTextSize(0.04);
  legpbpb_2->Draw();
  
  // centrality bin 30-50% 
  cATLAS_pbpb->cd(4);
  cATLAS_pbpb->cd(4)->SetLogy();
  p8719_d10x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d10x1y1->SetMarkerStyle(33);
  p8719_d10x1y1->SetMaximum(1e2);
  p8719_d10x1y1->SetMinimum(1e-4);
  //p8719_d10x1y1->Scale(1./TAA[3]);
  p8719_d10x1y1->Draw("ap");
  p8719_d11x1y1->SetMarkerStyle(34);
  //p8719_d11x1y1->Scale(1./TAA[4]);
  p8719_d11x1y1->Draw("psame");
  uPbPb_R4_Bayes[3]->SetMarkerStyle(24);
  uPbPb_R4_Bayes[3]->SetMarkerColor(kBlue);
  uPbPb_R4_Bayes[3]->Draw("same");
  systematics_R4.calcTotalSys(3);
  systematics_R4.Draw(uPbPb_R4_Bayes[3],3,2);
  
  TLegend *legpbpb_3 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_3->AddEntry(p8719_d10x1y1,"ATLAS 30-40%","p");
  legpbpb_3->AddEntry(p8719_d11x1y1,"ATLAS 40-50%","p");
  legpbpb_3->AddEntry(uPbPb_R4_Bayes[3],"CMS 30-50%","p");
  legpbpb_3->SetTextSize(0.04);
  legpbpb_3->Draw();

  // centrality bin 50-70% 
  cATLAS_pbpb->cd(3);
  cATLAS_pbpb->cd(3)->SetLogy();
  p8719_d12x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d12x1y1->SetMarkerStyle(33);
  p8719_d12x1y1->SetMaximum(1e2);
  p8719_d12x1y1->SetMinimum(1e-4);
  //p8719_d12x1y1->Scale(1./TAA[5]);
  p8719_d12x1y1->Draw("ap");
  p8719_d13x1y1->SetMarkerStyle(34);
  //p8719_d13x1y1->Scale(1./TAA[6]);
  p8719_d13x1y1->Draw("psame");
  uPbPb_R4_Bayes[4]->SetMarkerStyle(24);
  uPbPb_R4_Bayes[4]->SetMarkerColor(kBlue);
  uPbPb_R4_Bayes[4]->Draw("same");
  systematics_R4.calcTotalSys(4);
  systematics_R4.Draw(uPbPb_R4_Bayes[4],4,2);
  
  TLegend *legpbpb_4 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_4->AddEntry(p8719_d10x1y1,"ATLAS 50-60%","p");
  legpbpb_4->AddEntry(p8719_d11x1y1,"ATLAS 60-70%","p");
  legpbpb_4->AddEntry(uPbPb_R4_Bayes[4],"CMS 50-70%","p");
  legpbpb_4->SetTextSize(0.04);
  legpbpb_4->Draw();

  // centrality bin 70-90% 
  cATLAS_pbpb->cd(2);
  cATLAS_pbpb->cd(2)->SetLogy();
  p8719_d14x1y1->GetXaxis()->SetLimits(65, boundaries_pt[nbins_pt+1]);
  p8719_d14x1y1->SetMarkerStyle(33);
  p8719_d14x1y1->SetMaximum(1e2);
  p8719_d14x1y1->SetMinimum(1e-4);
  //p8719_d14x1y1->Scale(1./TAA[7]);
  p8719_d14x1y1->Draw("ap");
  uPbPb_R4_Bayes[5]->SetMarkerStyle(24);
  uPbPb_R4_Bayes[5]->SetMarkerColor(kBlue);
  uPbPb_R4_Bayes[5]->Draw("same");
  systematics_R4.calcTotalSys(5);
  systematics_R4.Draw(uPbPb_R4_Bayes[5],5,2);
  
  TLegend *legpbpb_5 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_5->AddEntry(p8719_d14x1y1,"ATLAS 70-80%","p");
  legpbpb_5->AddEntry(uPbPb_R4_Bayes[5],"CMS 70-90%","p");
  legpbpb_5->SetTextSize(0.04);
  legpbpb_5->Draw();

  cATLAS_pbpb->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_pbpb_spectra_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  cATLAS_pbpb->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_pbpb_spectra_%d_%s_pawan_ntuple.C",date.GetDate(),etaWidth),"RECREATE");
  cATLAS_pbpb->SaveAs(Form("June4/Analysis_full_bins/comparison_with_ATLAS_pbpb_spectra_%d_%s_pawan_ntuple.root",date.GetDate(),etaWidth),"RECREATE");


  // do the statistical error check
  // we have three histograms here:
  // 1) input data spectra.
  // 2) bayesian unfolding after 4 iterations - take it from the iteration systematics. 
  // 3) error fixed unfolding.

  // make those Error histograms:
  TH1F * hError_R2_Meas[nbins_cent], * hError_R2_4Iter[nbins_cent], * hError_R2_Fixed[nbins_cent];
  TH1F * hError_R3_Meas[nbins_cent], * hError_R3_4Iter[nbins_cent], * hError_R3_Fixed[nbins_cent];
  TH1F * hError_R4_Meas[nbins_cent], * hError_R4_4Iter[nbins_cent], * hError_R4_Fixed[nbins_cent];

  for(int i = 0; i<nbins_cent; ++i){
    cout<<i<<endl;
    hError_R2_Meas[i] = new TH1F(Form("hError_R2_Meas_cent%d",i),"",nbins_pt, boundaries_pt); 
    hError_R2_4Iter[i] = new TH1F(Form("hError_R2_4Iter_cent%d",i),"",nbins_pt, boundaries_pt); 
    hError_R2_Fixed[i] = new TH1F(Form("hError_R2_Fixed_cent%d",i),"",nbins_pt, boundaries_pt); 

    hError_R3_Meas[i] = new TH1F(Form("hError_R3_Meas_cent%d",i),"",nbins_pt, boundaries_pt); 
    hError_R3_4Iter[i] = new TH1F(Form("hError_R3_4Iter_cent%d",i),"",nbins_pt, boundaries_pt); 
    hError_R3_Fixed[i] = new TH1F(Form("hError_R3_Fixed_cent%d",i),"",nbins_pt, boundaries_pt); 

    hError_R4_Meas[i] = new TH1F(Form("hError_R4_Meas_cent%d",i),"",nbins_pt, boundaries_pt); 
    hError_R4_4Iter[i] = new TH1F(Form("hError_R4_4Iter_cent%d",i),"",nbins_pt, boundaries_pt); 
    hError_R4_Fixed[i] = new TH1F(Form("hError_R4_Fixed_cent%d",i),"",nbins_pt, boundaries_pt); 
    
    for(int j = 1; j<=nbins_pt; ++j){
      
      hError_R2_Meas[i]->SetBinContent(j, hPbPb_R2_measured[i]->GetBinError(j));
      hError_R2_Meas[i]->SetTitle(" ");
      hError_R2_Meas[i]->SetXTitle("Jet p_{T} (GeV/c)");
      hError_R2_Meas[i]->SetYTitle("Statistical Error Bar height");
      hError_R2_4Iter[i]->SetBinContent(j, uPbPb_BayesianIter_R2[i][4]->GetBinError(j));
      hError_R2_Fixed[i]->SetBinContent(j, hPbPb_R2_ErrorFix[i]->GetBinError(j));
      
      hError_R3_Meas[i]->SetBinContent(j, hPbPb_R3_measured[i]->GetBinError(j));
      hError_R3_Meas[i]->SetTitle(" ");
      hError_R3_Meas[i]->SetXTitle("Jet p_{T} (GeV/c)");
      hError_R3_Meas[i]->SetYTitle("Statistical Error Bar height");
      hError_R3_4Iter[i]->SetBinContent(j, uPbPb_BayesianIter_R3[i][4]->GetBinError(j));
      hError_R3_Fixed[i]->SetBinContent(j, hPbPb_R3_ErrorFix[i]->GetBinError(j));
      
      hError_R4_Meas[i]->SetBinContent(j, hPbPb_R4_measured[i]->GetBinError(j));
      hError_R4_Meas[i]->SetTitle(" ");
      hError_R4_Meas[i]->SetXTitle("Jet p_{T} (GeV/c)");
      hError_R4_Meas[i]->SetYTitle("Statistical Error Bar height");
      hError_R4_4Iter[i]->SetBinContent(j, uPbPb_BayesianIter_R4[i][4]->GetBinError(j));
      hError_R4_Fixed[i]->SetBinContent(j, hPbPb_R4_ErrorFix[i]->GetBinError(j));
      
    }
    
  }
 
  TCanvas * cErrorFix_R2 = new TCanvas("cErrorFix_R2","",1200,1000);
  makeMultiPanelCanvas(cErrorFix_R2,3,2,0.0,0.0,0.2,0.15,0.07);  
  for(int i = 0; i<nbins_cent; ++i){
    
    cErrorFix_R2->cd(nbins_cent-i);
    cErrorFix_R2->cd(nbins_cent-i)->SetLogy();

    hError_R2_Meas[i]->SetMarkerStyle(24);
    hError_R2_Meas[i]->SetMarkerColor(kBlack);
    hError_R2_Meas[i]->SetAxisRange(50, 299, "X");
    hError_R2_Meas[i]->Draw("p");

    hError_R2_4Iter[i]->SetMarkerStyle(25);
    hError_R2_4Iter[i]->SetMarkerColor(kRed);
    hError_R2_4Iter[i]->Draw("psame");
    hError_R2_Fixed[i]->SetMarkerStyle(33);
    hError_R2_Fixed[i]->SetMarkerColor(kRed);
    hError_R2_Fixed[i]->Draw("psame");
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  TLegend * err_R2 = myLegend(0.55,0.55,0.75,0.75);
  cErrorFix_R2->cd(1);
  putCMSPrel();
  err_R2->AddEntry(hError_R2_Meas[0],"Measured","pl");
  err_R2->AddEntry(hError_R2_4Iter[0],"4 Bayes Iterations","pl");
  err_R2->AddEntry(hError_R2_Fixed[0],"Data Driven Correction","pl");
  err_R2->SetTextSize(0.04);
  err_R2->Draw();
  
  cErrorFix_R2->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.pdf",etaWidth,2,date.GetDate()),"RECREATE");
  cErrorFix_R2->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.C",etaWidth,2,date.GetDate()),"RECREATE");
  cErrorFix_R2->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.root",etaWidth,2,date.GetDate()),"RECREATE");

  TCanvas * cErrorFix_R3 = new TCanvas("cErrorFix_R3","",1200,1000);
  makeMultiPanelCanvas(cErrorFix_R3,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0; i<nbins_cent; ++i){
    
    cErrorFix_R3->cd(nbins_cent-i);
    cErrorFix_R3->cd(nbins_cent-i)->SetLogy();

    hError_R3_Meas[i]->SetMarkerStyle(24);
    hError_R3_Meas[i]->SetMarkerColor(kBlack);
    hError_R3_Meas[i]->SetAxisRange(50, 299, "X");
    hError_R3_Meas[i]->Draw("p");

    hError_R3_4Iter[i]->SetMarkerStyle(25);
    hError_R3_4Iter[i]->SetMarkerColor(kRed);
    hError_R3_4Iter[i]->Draw("psame");
    hError_R3_Fixed[i]->SetMarkerStyle(33);
    hError_R3_Fixed[i]->SetMarkerColor(kRed);
    hError_R3_Fixed[i]->Draw("psame");
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  TLegend * err_R3 = myLegend(0.55,0.55,0.75,0.75);
  cErrorFix_R3->cd(1);
  putCMSPrel();
  err_R3->AddEntry(hError_R3_Meas[0],"Measured","pl");
  err_R3->AddEntry(hError_R3_4Iter[0],"4 Bayes Iterations","pl");
  err_R3->AddEntry(hError_R3_Fixed[0],"Data Driven Correction","pl");
  err_R3->SetTextSize(0.04);
  err_R3->Draw();
  
  cErrorFix_R3->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.pdf",etaWidth,3,date.GetDate()),"RECREATE");
  cErrorFix_R3->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.C",etaWidth,3,date.GetDate()),"RECREATE");
  cErrorFix_R3->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.root",etaWidth,3,date.GetDate()),"RECREATE");

  TCanvas * cErrorFix_R4 = new TCanvas("cErrorFix_R4","",1200,1000);
  makeMultiPanelCanvas(cErrorFix_R4,3,2,0.0,0.0,0.2,0.15,0.07);
  
  for(int i = 0; i<nbins_cent; ++i){
    
    cErrorFix_R4->cd(nbins_cent-i);
    cErrorFix_R4->cd(nbins_cent-i)->SetLogy();

    hError_R4_Meas[i]->SetMarkerStyle(24);
    hError_R4_Meas[i]->SetMarkerColor(kBlack);
    hError_R4_Meas[i]->SetAxisRange(50, 299, "X");
    hError_R4_Meas[i]->Draw("p");

    hError_R4_4Iter[i]->SetMarkerStyle(25);
    hError_R4_4Iter[i]->SetMarkerColor(kRed);
    hError_R4_4Iter[i]->Draw("psame");
    hError_R4_Fixed[i]->SetMarkerStyle(33);
    hError_R4_Fixed[i]->SetMarkerColor(kRed);
    hError_R4_Fixed[i]->Draw("psame");
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

    
  }
  TLegend * err_R4 = myLegend(0.55,0.55,0.75,0.75);
  cErrorFix_R4->cd(1);
  putCMSPrel();
  err_R4->AddEntry(hError_R4_Meas[0],"Measured","pl");
  err_R4->AddEntry(hError_R4_4Iter[0],"4 Bayes Iterations","pl");
  err_R4->AddEntry(hError_R4_Fixed[0],"Data Driven Correction","pl");
  err_R4->SetTextSize(0.04);
  err_R4->Draw();
  
  cErrorFix_R4->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.pdf",etaWidth,4,date.GetDate()),"RECREATE");
  cErrorFix_R4->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.C",etaWidth,4,date.GetDate()),"RECREATE");
  cErrorFix_R4->SaveAs(Form("June4/Analysis_full_bins/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.root",etaWidth,4,date.GetDate()),"RECREATE");
  // make plot to compare the bayesian unfolded spectra with 4 iterations and with the data driven error corrections. 
  
}
