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
#include "Headers/utilities.h"

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

using namespace std;

void RAA_plot_Systematics(int radius = 2, char *algo = "Pu", char *jet_type = "PF", int unfoldingCut = 40){
  
  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  char * etaWidth = (char*) "10_eta_10";
  Float_t etaBoundary = 2.0; 

  TFile *fin_R2, *fin_R3, *fin_R4; 
  fin_R2 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_noFakeWeight_unfold_mcclosure_oppside_trgMC_%s_%dGeVCut_ak%s_20150506.root",2,etaWidth,unfoldingCut,jet_type));
  fin_R3 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_noFakeWeight_unfold_mcclosure_oppside_trgMC_%s_%dGeVCut_ak%s_20150506.root",3,etaWidth,unfoldingCut,jet_type));
  fin_R4 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_noFakeWeight_unfold_mcclosure_oppside_trgMC_%s_%dGeVCut_ak%s_20150506.root",4,etaWidth,unfoldingCut,jet_type));

  // // get the unfolded error correction files and histograms
  TFile * fError_R2, * fError_R3, * fError_R4;
  fError_R2 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_%s_unfoldingCut_%d_noFakeWeight_data_driven_correction_akPu%s_20150506.root",2,2, etaWidth, unfoldingCut, jet_type));
  fError_R3 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_%s_unfoldingCut_%d_noFakeWeight_data_driven_correction_akPu%s_20150506.root",3,3, etaWidth, unfoldingCut, jet_type));
  fError_R4 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pawan_ntuple_PbPb_R%d_pp_R%d_noJetID_%s_unfoldingCut_%d_noFakeWeight_data_driven_correction_akPu%s_20150506.root",4,4, etaWidth, unfoldingCut, jet_type));


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
  
  for(int i = 0; i<nbins_cent; ++i){
    cout<<i<<endl;
    hPbPb_R2_ErrorFix[i] = (TH1F*)fError_R2->Get(Form("PbPb_BayesianUnfolded_cent%d",i));
    hPbPb_R3_ErrorFix[i] = (TH1F*)fError_R3->Get(Form("PbPb_BayesianUnfolded_cent%d",i));
    hPbPb_R4_ErrorFix[i] = (TH1F*)fError_R4->Get(Form("PbPb_BayesianUnfolded_cent%d",i));

    hPbPb_R2_measured[i] = (TH1F*)fError_R2->Get(Form("hpbpb_HLTComb_R2_%s_cent%d",etaWidth,i));
    hPbPb_R3_measured[i] = (TH1F*)fError_R3->Get(Form("hpbpb_HLTComb_R3_%s_cent%d",etaWidth,i));
    hPbPb_R4_measured[i] = (TH1F*)fError_R4->Get(Form("hpbpb_HLTComb_R4_%s_cent%d",etaWidth,i));

    uPbPb_R2_Bayes[i] = (TH1F*)fin_R2->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R2_Bayes[i]->Print("base");
    uPbPb_R3_Bayes[i] = (TH1F*)fin_R3->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R3_Bayes[i]->Print("base");
    uPbPb_R4_Bayes[i] = (TH1F*)fin_R4->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R3_Bayes[i]->Print("base");
 
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
    
    RAA_R2_BinByBin[i]   = (TH1F*)fin_R2->Get(Form("RAA_binbybin_cent%d",i));  
    RAA_R3_BinByBin[i]   = (TH1F*)fin_R3->Get(Form("RAA_binbybin_cent%d",i));  
    RAA_R4_BinByBin[i]   = (TH1F*)fin_R4->Get(Form("RAA_binbybin_cent%d",i));  
    
    RAA_R2_Meas[i]   = (TH1F*)fin_R2->Get(Form("RAA_measured_cent%d",i));  
    RAA_R3_Meas[i]   = (TH1F*)fin_R3->Get(Form("RAA_measured_cent%d",i));  
    RAA_R4_Meas[i]   = (TH1F*)fin_R4->Get(Form("RAA_measured_cent%d",i));  
    
  }
  // declare the systematics
  SysData systematics_R2;
  SysData systematics_R3;
  SysData systematics_R4;

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

  //make the canvas
  // line at 1
  TLine *linePbPb_iter = new TLine(unfoldingCut,1,299,1);
  linePbPb_iter->SetLineStyle(2);
  linePbPb_iter->SetLineWidth(2);
  

  TCanvas * cPbPb_Its_R2 = new TCanvas("cPbPb_Its_R2","",1000,800);
  makeMultiPanelCanvasWithGap(cPbPb_Its_R2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  TLegend *PbPb_itersys_R2 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R2->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R2[j][i]->SetMarkerStyle(33);
      hNumerator_R2[j][i]->SetMarkerColor(j+1);
      hNumerator_R2[j][i]->SetAxisRange(unfoldingCut,299,"X");
      hNumerator_R2[j][i]->SetAxisRange(0,2,"Y");

      if(j==1){
	makeHistTitle(hNumerator_R2[j][i]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
	hNumerator_R2[j][i]->Draw();
      }
      hNumerator_R2[j][i]->Draw("same");
	
      if(i==0) PbPb_itersys_R2->AddEntry(hNumerator_R2[j][i],Form("Iteration %d",j+1),"pl");
      checkMaximumSys(systematics_R2.hSysIter[i], hNumerator_R2[j][i], 0, 1.05);
      
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

  cPbPb_Its_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%sR2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  

  TCanvas * cPP_Its_R2 = new TCanvas("cPP_Its_R2","",600,400);
  // line at 1
  TLine *linePP_iter_R2 = new TLine(unfoldingCut,1,299,1);
  linePP_iter_R2->SetLineStyle(2);
  linePP_iter_R2->SetLineWidth(2);
  
  TLegend *PP_itersys_R2 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int j = 1;j<Iterations;j++){
    hNumeratorPP_R2[j]->SetMarkerStyle(33);
    hNumeratorPP_R2[j]->SetMarkerColor(j+1);
    hNumeratorPP_R2[j]->SetAxisRange(unfoldingCut,299,"X");
    hNumeratorPP_R2[j]->SetAxisRange(0,2,"Y");

    if(j==1){
      makeHistTitle(hNumeratorPP_R2[j]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
      hNumeratorPP_R2[j]->Draw();
    }
    hNumeratorPP_R2[j]->Draw("same");

    checkMaximumSys(systematics_R2.hSysIter[nbins_cent], hNumeratorPP_R2[j],0, 1.05);
    PP_itersys_R2->AddEntry(hNumeratorPP_R2[j],Form("Iteration %d",j+1),"pl");

  }

  linePP_iter_R2->Draw();  
  PP_itersys_R2->Draw();
  drawEnvelope(systematics_R2.hSysIter[nbins_cent],"hist same");

  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, 2),0.2,0.23,20);
  drawText(Form("|#eta|< %2.0f, |vz|<15",etaBoundary),0.6,0.31,20);

  cPP_Its_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_akR2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");

  // define the Jet ID efficiency as a function of jet pT: taken from the above plots: 
  TF1 * fPol = new TF1("fPol","[0]+[1]*x");
  
  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R2 = new TCanvas("cRAA_JEC_sys_R2","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_JEC_sys_R2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R2->cd(nbins_cent-i);
    RAA_JEC_bayesian_R2[i]->Divide(RAA_JEC_bayesian_R2[i], RAA_bayesian_R2[i], 1,1, "B");
    RAA_JEC_bayesian_R2[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R2[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R2[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_JEC_bayesian_R2[i]->Fit("fPol","","",60,299);

    RAA_JEC_bayesian_R2[i]->Draw("p");
    checkMaximumSys(systematics_R2.hSysJEC[i], functionHist(fPol, systematics_R2.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_JEC_sys_R2->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R0.2",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|< %2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_JEC_sys_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

    // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R2 = new TCanvas("cRAA_Smear_sys_R2","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Smear_sys_R2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R2->cd(nbins_cent-i);
    RAA_Smear_bayesian_R2[i]->Divide(RAA_Smear_bayesian_R2[i], RAA_bayesian_R2[i], 1,1, "B");
    RAA_Smear_bayesian_R2[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R2[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R2[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_Smear_bayesian_R2[i]->Fit("fPol","","",60,299);
    RAA_Smear_bayesian_R2[i]->Draw("p");
    checkMaximumSys(systematics_R2.hSysSmear[i], functionHist(fPol, systematics_R2.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_Smear_sys_R2->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R2",algo, jet_type, radius),0.2,0.23,20);
  if(etaWidth=="n16_eta_p16")drawText("|#eta|<1.6, |vz|<15",0.6,0.31,20);
  if(etaWidth=="n20_eta_p20")drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_Smear_sys_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));


  // draw the total systematics
  TCanvas * cSys_R2 = new TCanvas("cSys_R2","Total Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cSys_R2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0; i<nbins_cent; ++i){

    cSys_R2->cd(nbins_cent-i);
    cSys_R2->cd(nbins_cent-i)->SetGridy();
    cSys_R2->cd(nbins_cent-i)->SetGridx();
    systematics_R2.DrawComponent(i);

    //drawPanelLabel(i);
    TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    title_->AddEntry(RAA_bayesian_R2[i],Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    title_->SetTextSize(0.06);
    title_->Draw();
  
  }
  cSys_R2->cd(1);
  putCMSPrel();
  linePbPb_iter->Draw();
  cSys_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_total_systematics_unfoldingCut%dGeV_ak%sR2%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");


  TCanvas * cPbPb_Its_R3 = new TCanvas("cPbPb_Its_R3","",1000,800);
  makeMultiPanelCanvasWithGap(cPbPb_Its_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  TLegend *PbPb_itersys_R3 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R3->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R3[j][i]->SetMarkerStyle(33);
      hNumerator_R3[j][i]->SetMarkerColor(j+1);
      hNumerator_R3[j][i]->SetAxisRange(unfoldingCut,299,"X");
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
  
  cPbPb_Its_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%sR3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  

  TCanvas * cPP_Its_R3 = new TCanvas("cPP_Its_R3","",600,400);
  // line at 1
  TLine *linePP_iter_R3 = new TLine(unfoldingCut,1,299,1);
  linePP_iter_R3->SetLineStyle(2);
  linePP_iter_R3->SetLineWidth(2);
  
  TLegend *PP_itersys_R3 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int j = 1;j<Iterations;j++){
    hNumeratorPP_R3[j]->SetMarkerStyle(33);
    hNumeratorPP_R3[j]->SetMarkerColor(j+1);
    hNumeratorPP_R3[j]->SetAxisRange(unfoldingCut,299,"X");
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
  
  cPP_Its_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_akR3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");

  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R3 = new TCanvas("cRAA_JEC_sys_R3","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_JEC_sys_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R3->cd(nbins_cent-i);
    RAA_JEC_bayesian_R3[i]->Divide(RAA_JEC_bayesian_R3[i], RAA_bayesian_R3[i], 1,1, "B");
    RAA_JEC_bayesian_R3[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R3[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R3[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_JEC_bayesian_R3[i]->Fit("fPol","","",60,299);

    RAA_JEC_bayesian_R3[i]->Draw("p");
    checkMaximumSys(systematics_R3.hSysJEC[i], functionHist(fPol, systematics_R3.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_JEC_sys_R3->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.3",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_JEC_sys_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

    // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R3 = new TCanvas("cRAA_Smear_sys_R3","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Smear_sys_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R3->cd(nbins_cent-i);
    RAA_Smear_bayesian_R3[i]->Divide(RAA_Smear_bayesian_R3[i], RAA_bayesian_R3[i], 1,1, "B");
    RAA_Smear_bayesian_R3[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R3[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R3[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_Smear_bayesian_R3[i]->Fit("fPol","","",60,299);
    RAA_Smear_bayesian_R3[i]->Draw("p");
    checkMaximumSys(systematics_R3.hSysSmear[i], functionHist(fPol, systematics_R3.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_Smear_sys_R3->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.3",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_Smear_sys_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));


  // draw the total systematics
  TCanvas * cSys_R3 = new TCanvas("cSys_R3","Total Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cSys_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0; i<nbins_cent; ++i){

    cSys_R3->cd(nbins_cent-i);
    cSys_R3->cd(nbins_cent-i)->SetGridy();
    cSys_R3->cd(nbins_cent-i)->SetGridx();
    systematics_R3.DrawComponent(i);

    //drawPanelLabel(i);
    TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    title_->AddEntry(RAA_bayesian_R3[i],Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    title_->SetTextSize(0.06);
    title_->Draw();
  
  }
  cSys_R3->cd(1);
  putCMSPrel();
  linePbPb_iter->Draw();
  cSys_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_total_systematics_unfoldingCut%dGeV_ak%sR3%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");



  TCanvas * cPbPb_Its_R4 = new TCanvas("cPbPb_Its_R4","",1000,800);
  makeMultiPanelCanvasWithGap(cPbPb_Its_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  TLegend *PbPb_itersys_R4 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R4->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R4[j][i]->SetMarkerStyle(33);
      hNumerator_R4[j][i]->SetMarkerColor(j+1);
      hNumerator_R4[j][i]->SetAxisRange(unfoldingCut,299,"X");
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
  cPbPb_Its_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%sR4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut,algo,jet_type,date.GetDate(),etaWidth),"RECREATE");
  

  TCanvas * cPP_Its_R4 = new TCanvas("cPP_Its_R4","",600,400);
  // line at 1
  TLine *linePP_iter_R4 = new TLine(unfoldingCut,1,299,1);
  linePP_iter_R4->SetLineStyle(2);
  linePP_iter_R4->SetLineWidth(2);
  
  TLegend *PP_itersys_R4 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int j = 1;j<Iterations;j++){
    hNumeratorPP_R4[j]->SetMarkerStyle(33);
    hNumeratorPP_R4[j]->SetMarkerColor(j+1);
    hNumeratorPP_R4[j]->SetAxisRange(unfoldingCut,299,"X");
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
  cPP_Its_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_akR4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut,jet_type,date.GetDate(),etaWidth),"RECREATE");

  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R4 = new TCanvas("cRAA_JEC_sys_R4","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_JEC_sys_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R4->cd(nbins_cent-i);
    RAA_JEC_bayesian_R4[i]->Divide(RAA_JEC_bayesian_R4[i], RAA_bayesian_R4[i], 1,1, "B");
    RAA_JEC_bayesian_R4[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R4[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R4[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_JEC_bayesian_R4[i]->Fit("fPol","","",60,299);

    RAA_JEC_bayesian_R4[i]->Draw("p");
    checkMaximumSys(systematics_R4.hSysJEC[i], functionHist(fPol, systematics_R4.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_JEC_sys_R4->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.4",algo, jet_type),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_JEC_sys_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

    // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R4 = new TCanvas("cRAA_Smear_sys_R4","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Smear_sys_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R4->cd(nbins_cent-i);
    RAA_Smear_bayesian_R4[i]->Divide(RAA_Smear_bayesian_R4[i], RAA_bayesian_R4[i], 1,1, "B");
    RAA_Smear_bayesian_R4[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R4[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R4[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_Smear_bayesian_R4[i]->Fit("fPol","","",60,299);
    RAA_Smear_bayesian_R4[i]->Draw("p");
    checkMaximumSys(systematics_R4.hSysSmear[i], functionHist(fPol, systematics_R4.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_Smear_sys_R4->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R4",algo, jet_type, radius),0.2,0.23,20);
  drawText(Form("|#eta|<%2.0f, |vz|<15",etaBoundary),0.6,0.31,20);
  cRAA_Smear_sys_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth));


  // draw the total systematics
  TCanvas * cSys_R4 = new TCanvas("cSys_R4","Total Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cSys_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0; i<nbins_cent; ++i){

    cSys_R4->cd(nbins_cent-i);
    cSys_R4->cd(nbins_cent-i)->SetGridy();
    cSys_R4->cd(nbins_cent-i)->SetGridx();
    systematics_R4.DrawComponent(i);

    //drawPanelLabel(i);
    TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    title_->AddEntry(RAA_bayesian_R4[i],Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    title_->SetTextSize(0.06);
    title_->Draw();
  
  }
  cSys_R4->cd(1);
  putCMSPrel();
  linePbPb_iter->Draw();
  cSys_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_total_systematics_unfoldingCut%dGeV_ak%sR4%s_%d_%s_pawan_ntuple.pdf",unfoldingCut, algo, jet_type, date.GetDate(),etaWidth),"RECREATE");

  
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
  TLine *lineRAA_R2 = new TLine(unfoldingCut+30,1,299,1);
  lineRAA_R2->SetLineStyle(2);
  lineRAA_R2->SetLineWidth(2);

  TCanvas * cRAA_npart = new TCanvas("cRAA_npart","",600,400);
  cRAA_npart->SetGridy();
  cRAA_npart->SetGridx();

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
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.6,0.31,20);
  drawText("|vz|<15, HBHEfilter, pCES",0.15,0.2,16);
  drawText("97 < Jet p_{T} < 114", 0.15,0.8,16);
  
  TLegend * npart1 = myLegend(0.7,0.7,0.9,0.9);
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
  lineRAA_R2->Draw();
  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  
  cRAA_npart->SaveAs(Form("../../Plots/Final_paper_plots_RAA_npart_97_pT_114_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");

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


  TCanvas * cRAA_npart2 = new TCanvas("cRAA_npart2","",600,400);
  cRAA_npart2->SetGridy();
  cRAA_npart2->SetGridx();

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
  lineRAA_R2->Draw();
  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.6,0.31,20);
  drawText("|vz|<15, HBHEfilter, pCES",0.15,0.2,16);
  drawText("133 < Jet p_{T} < 153", 0.15,0.8,16);
  
  TLegend * npart2 = myLegend(0.7,0.7,0.9,0.9);
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

  
  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  
  cRAA_npart2->SaveAs(Form("../../Plots/Final_paper_plots_RAA_npart_133_pT_153_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");

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


  TCanvas * cRAA_npart3 = new TCanvas("cRAA_npart3","",600,400);
  cRAA_npart3->SetGridy();
  cRAA_npart3->SetGridx();

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
  lineRAA_R2->Draw();
  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.6,0.2,16);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.6,0.31,20);
  drawText("|vz|<15, HBHEfilter, pCES",0.15,0.2,16);
  drawText("74 < Jet p_{T} < 84", 0.15,0.8,16);
  
  TLegend * npart3 = myLegend(0.7,0.7,0.9,0.9);
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

  
  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  
  cRAA_npart3->SaveAs(Form("../../Plots/Final_paper_plots_RAA_npart_74_pT_84_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");

  
  
  // plot the RAA for each radii here along with the systematics shown in yellow, for measured, bayesian and bin by bin unfoldeding.
 TCanvas *cRAA_R2 = new TCanvas("cRAA_R2","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA_R2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA_R2 = myLegend(0.45,0.75,0.85,0.9);

    
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
    RAA_R2_Meas[i]->SetAxisRange(unfoldingCut+30,299,"X");
    RAA_R2_Meas[i]->SetAxisRange(0,2,"Y");
    RAA_R2_Meas[i]->Draw("E1");

    RAA_R2_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R2_Bayes[i]->SetMarkerStyle(20);
    RAA_R2_Bayes[i]->Draw("same E1");

    RAA_R2_BinByBin[i]->SetMarkerStyle(25);
    RAA_R2_BinByBin[i]->SetMarkerColor(kRed);
    RAA_R2_BinByBin[i]->Draw("same E1");

    lineRAA_R2->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

    systematics_R2.calcTotalSys(i);
    systematics_R2.Draw(RAA_R2_Bayes[i],i,2);
    
    RAA_R2_Meas[i]->Draw("same E1");
    RAA_R2_BinByBin[i]->Draw("same E1");
    RAA_R2_Bayes[i]->Draw("same E1");

  }
    
  tRAA_R2->AddEntry(RAA_R2_Meas[0],"Measured","pl");
  tRAA_R2->AddEntry(RAA_R2_Bayes[0],"Bayesian","pl");
  tRAA_R2->AddEntry(RAA_R2_BinByBin[0],"BinByBin","pl");
  tRAA_R2->SetTextSize(0.04);

  cRAA_R2->cd(1);
  tRAA_R2->Draw();
  cRAA_R2->cd(1);
  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  drawText(Form("Anti-k_{T} %s R=0.2 %s Jets",algo,jet_type),0.2,0.23,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_R2->cd(2);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.6,0.31,20);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  cRAA_R2->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.3,16);
  drawText("Pile up rejection cut applied",0.1,0.2,16);

  cRAA_R2->SaveAs(Form("../../Plots/Final_paper_plots_RAA_R2_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
    
  // draw it for R=0.3
  TCanvas *cRAA_R3 = new TCanvas("cRAA_R3","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA_R3 = myLegend(0.45,0.75,0.85,0.9);
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
    RAA_R3_Meas[i]->SetAxisRange(unfoldingCut+30,299,"X");
    RAA_R3_Meas[i]->SetAxisRange(0,2,"Y");
    RAA_R3_Meas[i]->Draw("E1");

    RAA_R3_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R3_Bayes[i]->SetMarkerStyle(20);
    RAA_R3_Bayes[i]->Draw("same E1");

    RAA_R3_BinByBin[i]->SetMarkerStyle(25);
    RAA_R3_BinByBin[i]->SetMarkerColor(kRed);
    RAA_R3_BinByBin[i]->Draw("same E1");

    lineRAA_R3->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

    systematics_R3.calcTotalSys(i);
    systematics_R3.Draw(RAA_R3_Bayes[i],i,2);
    
    RAA_R3_Meas[i]->Draw("same E1");
    RAA_R3_BinByBin[i]->Draw("same E1");
    RAA_R3_Bayes[i]->Draw("same E1");

  }
    
  tRAA_R3->AddEntry(RAA_R3_Meas[0],"Measured","pl");
  tRAA_R3->AddEntry(RAA_R3_Bayes[0],"Bayesian","pl");
  tRAA_R3->AddEntry(RAA_R3_BinByBin[0],"BinByBin","pl");
  tRAA_R3->SetTextSize(0.04);

  cRAA_R3->cd(1);
  tRAA_R3->Draw();
  cRAA_R3->cd(1);
  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  drawText(Form("Anti-k_{T} %s R=0.3 %s Jets",algo,jet_type),0.2,0.23,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_R3->cd(2);
  drawText(Form("|#eta|<%2.0f, |vz|<15, Jet ID Cut",etaBoundary),0.6,0.31,20);

  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  cRAA_R3->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.3,16);
  drawText("Pile up rejection cut applied",0.1,0.2,16);

  cRAA_R3->SaveAs(Form("../../Plots/Final_paper_plots_RAA_R3_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");

  // plot it for R=0.4
   TCanvas *cRAA_R4 = new TCanvas("cRAA_R4","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA_R4 = myLegend(0.45,0.75,0.85,0.9);
  TLine *lineRAA_R4 = new TLine(unfoldingCut+30,1,299,1);
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
    RAA_R4_Meas[i]->SetAxisRange(unfoldingCut+30,299,"X");
    RAA_R4_Meas[i]->SetAxisRange(0,2,"Y");
    RAA_R4_Meas[i]->Draw("E1");

    RAA_R4_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R4_Bayes[i]->SetMarkerStyle(20);
    RAA_R4_Bayes[i]->Draw("same E1");

    RAA_R4_BinByBin[i]->SetMarkerStyle(25);
    RAA_R4_BinByBin[i]->SetMarkerColor(kRed);
    RAA_R4_BinByBin[i]->Draw("same E1");

    lineRAA_R4->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

    systematics_R4.calcTotalSys(i);
    systematics_R4.Draw(RAA_R4_Bayes[i],i,2);
    
    RAA_R4_Meas[i]->Draw("same E1");
    RAA_R4_BinByBin[i]->Draw("same E1");
    RAA_R4_Bayes[i]->Draw("same E1");

  }
    
  tRAA_R4->AddEntry(RAA_R4_Meas[0],"Measured","pl");
  tRAA_R4->AddEntry(RAA_R4_Bayes[0],"Bayesian","pl");
  tRAA_R4->AddEntry(RAA_R4_BinByBin[0],"BinByBin","pl");
  tRAA_R4->SetTextSize(0.04);

  cRAA_R4->cd(1);
  tRAA_R4->Draw();
  cRAA_R4->cd(1);
  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  drawText(Form("Anti-k_{T} %s R=0.4 %s Jets",algo,jet_type),0.2,0.23,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA_R4->cd(2);
  drawText(Form("|#eta|<%2.0f , Jet ID Cut",etaBoundary),0.6,0.31,20);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  cRAA_R4->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.3,16);
  drawText("Pile up rejection cut applied",0.1,0.2,16);

  cRAA_R4->SaveAs(Form("../../Plots/Final_paper_plots_RAA_R4_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");


  
  // plot 2 - Bayesian unfolded RAA as a function of pT for the different radii
  //        - regular 6 panel plot 
  
  // again this will be a 6 panel plot. showing measured, unfolded Bayesian, and unfolded Bin By Bin methods. 
  TCanvas *cRAA = new TCanvas("cRAA","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA = myLegend(0.45,0.75,0.85,0.9);
  TLine *lineRAA = new TLine(unfoldingCut,1,299,1);
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
    RAA_R2_Bayes[i]->SetAxisRange(unfoldingCut+30,299,"X");
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
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

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
  putPbPbLumi();
  putPPLumi();
  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.2,0.23,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA->cd(2);
  drawText(Form("Jet ID cut, |#eta|<2=%2.0f",etaBoundary),0.1,0.3,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  cRAA->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.3,16);
  drawText("Pile up rejection cut applied",0.1,0.2,16);

  cRAA->SaveAs(Form("../../Plots/Final_paper_plots_RAA_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
    

  
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

  TCanvas * cATLAS = new TCanvas("cATLAS","",800,600);
  makeMultiPanelCanvasWithGap(cATLAS,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  cATLAS->cd(6);

  // centrality bin - 0-10
  p8719_d27x1y1->SetMaximum(2);
  p8719_d27x1y1->SetMinimum(0);
  p8719_d27x1y1->SetMarkerStyle(33);
  p8719_d27x1y1->GetXaxis()->SetLimits(50,299);
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
  p8719_d28x1y1->GetXaxis()->SetLimits(50,299);
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
  p8719_d30x1y1->GetXaxis()->SetLimits(50,299);
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
  p8719_d32x1y1->GetXaxis()->SetLimits(50,299);
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
  p8719_d34x1y1->GetXaxis()->SetLimits(50,299);
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
 
  
  cATLAS->SaveAs(Form("../../Plots/comparison_with_ATLAS_RAA_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");



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

  TCanvas * cATLAS_pp = new TCanvas("cATLAS_pp","",800,600);
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

  TLegend *legpp_1 = getLegend(0.20,0.65,0.4,0.85);
  legpp_1->AddEntry(p8719_d2x1y1,"ATLAS pp","p");
  legpp_1->AddEntry(uPP_R4_Bayes,"CMS pp","p");
  legpp_1->SetTextSize(0.04);
  legpp_1->Draw();
  
  cATLAS_pp->SaveAs(Form("../../Plots/comparison_with_ATLAS_pp_spectra_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");
  
  
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
  double p8719_d8x1y1_yval[] = { 1.28E-4, 3.7E-5, 1.02E-5, 2.85E-5, 7.38E-6, 1.86E-7, 3.82E-7, 6.97E-8, 1.07E-9 };
  double p8719_d8x1y1_yerrminus[] = { 2.365369349594266E-5, 7.71625219909251E-6, 1.7929093340155267E-6, 5.048523818899937E-6, 1.2186165641414859E-6, 3.001466041786913E-8, 6.654369827414165E-8, 1.3162416403533205E-8, 2.3045773169932924E-10 };
  double p8719_d8x1y1_yerrplus[] = { 2.365369349594266E-5, 7.71625219909251E-6, 1.7929093340155267E-6, 5.048523818899937E-6, 1.2186165641414859E-6, 3.001466041786913E-8, 6.654369827414165E-8, 1.3162416403533205E-8, 2.3045773169932924E-10 };
  double p8719_d8x1y1_ystatminus[] = { 2.8160000000000002E-6, 8.88E-7, 1.734E-7, 4.2750000000000004E-7, 1.1808000000000001E-7, 3.72E-9, 1.1842000000000001E-8, 3.9729E-9, 1.0058000000000002E-10 };
  double p8719_d8x1y1_ystatplus[] = { 2.8160000000000002E-6, 8.88E-7, 1.734E-7, 4.2750000000000004E-7, 1.1808000000000001E-7, 3.72E-9, 1.1842000000000001E-8, 3.9729E-9, 1.0058000000000002E-10 };
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

  //for(int i = 0; i<nbins_cent;++i){
    // uPbPb_R4_Bayes[i]->Scale(145.156); //remove lumi scaling
    // uPbPb_R4_Bayes[i]->Scale(7.46*1e6); // remove pp sigma
    // uPbPb_R4_Bayes[i]->Scale((ncoll[i]*1e3)/(64./1e9)); // remove  ncoll and other pp related normalization
    // uPbPb_R4_Bayes[i]->Scale(1./4.88e10); // apply scaling by no of events to get yield. 
    // uPbPb_R4_Bayes[i]->Scale(1./1e5);
  //}

  TCanvas * cATLAS_pbpb = new TCanvas("cATLAS_pbpb","",800,600);
  makeMultiPanelCanvasWithGap(cATLAS_pbpb,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  // centrality bin 0-10% 
  cATLAS_pbpb->cd(6);
  cATLAS_pbpb->cd(6)->SetLogy();
  p8719_d7x1y1->GetXaxis()->SetLimits(50,299);
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

  TLegend *legpbpb_1 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_1->AddEntry(p8719_d7x1y1,"ATLAS 0-10%","p");
  legpbpb_1->AddEntry(uPbPb_R4_Bayes[0],"CMS 0-5%","p");
  legpbpb_1->AddEntry(uPbPb_R4_Bayes[1],"CMS 5-10%","p");
  legpbpb_1->SetTextSize(0.04);
  legpbpb_1->Draw();

  // centrality bin 10-30% 
  cATLAS_pbpb->cd(5);
  cATLAS_pbpb->cd(5)->SetLogy();
  p8719_d8x1y1->GetXaxis()->SetLimits(50,299);
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
  
  TLegend *legpbpb_2 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_2->AddEntry(p8719_d8x1y1,"ATLAS 10-20%","p");
  legpbpb_2->AddEntry(p8719_d9x1y1,"ATLAS 20-30%","p");
  legpbpb_2->AddEntry(uPbPb_R4_Bayes[2],"CMS 10-30%","p");
  legpbpb_2->SetTextSize(0.04);
  legpbpb_2->Draw();
  
  // centrality bin 30-50% 
  cATLAS_pbpb->cd(4);
  cATLAS_pbpb->cd(4)->SetLogy();
  p8719_d10x1y1->GetXaxis()->SetLimits(50,299);
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
  
  TLegend *legpbpb_3 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_3->AddEntry(p8719_d10x1y1,"ATLAS 30-40%","p");
  legpbpb_3->AddEntry(p8719_d11x1y1,"ATLAS 40-50%","p");
  legpbpb_3->AddEntry(uPbPb_R4_Bayes[3],"CMS 30-50%","p");
  legpbpb_3->SetTextSize(0.04);
  legpbpb_3->Draw();

  // centrality bin 50-70% 
  cATLAS_pbpb->cd(3);
  cATLAS_pbpb->cd(3)->SetLogy();
  p8719_d12x1y1->GetXaxis()->SetLimits(50,299);
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
  
  TLegend *legpbpb_4 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_4->AddEntry(p8719_d10x1y1,"ATLAS 50-60%","p");
  legpbpb_4->AddEntry(p8719_d11x1y1,"ATLAS 60-70%","p");
  legpbpb_4->AddEntry(uPbPb_R4_Bayes[4],"CMS 50-70%","p");
  legpbpb_4->SetTextSize(0.04);
  legpbpb_4->Draw();

  // centrality bin 70-90% 
  cATLAS_pbpb->cd(2);
  cATLAS_pbpb->cd(2)->SetLogy();
  p8719_d14x1y1->GetXaxis()->SetLimits(50,299);
  p8719_d14x1y1->SetMarkerStyle(33);
  p8719_d14x1y1->SetMaximum(1e2);
  p8719_d14x1y1->SetMinimum(1e-4);
  //p8719_d14x1y1->Scale(1./TAA[7]);
  p8719_d14x1y1->Draw("ap");
  uPbPb_R4_Bayes[5]->SetMarkerStyle(24);
  uPbPb_R4_Bayes[5]->SetMarkerColor(kBlue);
  uPbPb_R4_Bayes[5]->Draw("same");
  
  TLegend *legpbpb_5 = getLegend(0.60,0.65,0.8,0.85);
  legpbpb_5->AddEntry(p8719_d14x1y1,"ATLAS 70-80%","p");
  legpbpb_5->AddEntry(uPbPb_R4_Bayes[5],"CMS 70-90%","p");
  legpbpb_5->SetTextSize(0.04);
  legpbpb_5->Draw();

  cATLAS_pbpb->SaveAs(Form("../../Plots/comparison_with_ATLAS_pbpb_spectra_%d_%s_pawan_ntuple.pdf",date.GetDate(),etaWidth),"RECREATE");


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
 
  TCanvas * cErrorFix_R2 = new TCanvas("cErrorFix_R2","",800,600);
  makeMultiPanelCanvasWithGap(cErrorFix_R2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);  
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
  
  cErrorFix_R2->SaveAs(Form("../../Plots/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.pdf",etaWidth,2,date.GetDate()),"RECREATE");

  TCanvas * cErrorFix_R3 = new TCanvas("cErrorFix_R3","",800,600);
  makeMultiPanelCanvasWithGap(cErrorFix_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
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
  
  cErrorFix_R3->SaveAs(Form("../../Plots/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.pdf",etaWidth,3,date.GetDate()),"RECREATE");

  TCanvas * cErrorFix_R4 = new TCanvas("cErrorFix_R4","",800,600);
  makeMultiPanelCanvasWithGap(cErrorFix_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
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
  
  cErrorFix_R4->SaveAs(Form("../../Plots/UnfoldingErrorFix_PbPb_%s_R%d_%d_pawan_ntuple.pdf",etaWidth,4,date.GetDate()),"RECREATE");
  // make plot to compare the bayesian unfolded spectra with 4 iterations and with the data driven error corrections. 

  
}
