// Raghav Kunnawalkam Elayavalli
// April 13 2014
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


using namespace std;

void RAA_plot_Systematics(int radius = 4, char *algo = "Pu", char *jet_type = "PF", int unfoldingCut = 30){

  
  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();



  TFile *fin_R2, *fin_R3, *fin_R4; 
  fin_R2 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_unfold_mcclosure_oppside_trgMC_n20_eta_p20_%dGeVCut_ak%s_20150417.root",2,unfoldingCut,jet_type));
  fin_R3 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_unfold_mcclosure_oppside_trgMC_n20_eta_p20_%dGeVCut_ak%s_20150417.root",3,unfoldingCut,jet_type));
  fin_R4 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_unfold_mcclosure_oppside_trgMC_n20_eta_p20_%dGeVCut_ak%s_20150417.root",4,unfoldingCut,jet_type));

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
      hNumerator_R2[j][i]->SetAxisRange(0.4,1.6,"Y");

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
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPbPb_Its_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%sR2%s_%d.pdf",unfoldingCut,algo,jet_type,date.GetDate()),"RECREATE");
  

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
    hNumeratorPP_R2[j]->SetAxisRange(0.4,1.6,"Y");

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
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPP_Its_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_akR2%s_%d.pdf",unfoldingCut,jet_type,date.GetDate()),"RECREATE");

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
    RAA_JEC_bayesian_R2[i]->Fit("fPol","");

    RAA_JEC_bayesian_R2[i]->Draw("p");
    checkMaximumSys(systematics_R2.hSysJEC[i], functionHist(fPol, systematics_R2.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_JEC_sys_R2->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R2",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_JEC_sys_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR2%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()),"RECREATE");

    // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R2 = new TCanvas("cRAA_Smear_sys_R2","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Smear_sys_R2,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R2->cd(nbins_cent-i);
    RAA_Smear_bayesian_R2[i]->Divide(RAA_Smear_bayesian_R2[i], RAA_bayesian_R2[i], 1,1, "B");
    RAA_Smear_bayesian_R2[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R2[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R2[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_Smear_bayesian_R2[i]->Fit("fPol","");
    RAA_Smear_bayesian_R2[i]->Draw("p");
    checkMaximumSys(systematics_R2.hSysSmear[i], functionHist(fPol, systematics_R2.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_Smear_sys_R2->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R2",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_Smear_sys_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R2%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()));


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
  cSys_R2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_total_systematics_unfoldingCut%dGeV_ak%sR2%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()),"RECREATE");


  TCanvas * cPbPb_Its_R3 = new TCanvas("cPbPb_Its_R3","",1000,800);
  makeMultiPanelCanvasWithGap(cPbPb_Its_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  TLegend *PbPb_itersys_R3 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R3->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R3[j][i]->SetMarkerStyle(33);
      hNumerator_R3[j][i]->SetMarkerColor(j+1);
      hNumerator_R3[j][i]->SetAxisRange(unfoldingCut,299,"X");
      hNumerator_R3[j][i]->SetAxisRange(0.4,1.6,"Y");

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
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPbPb_Its_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%sR3%s_%d.pdf",unfoldingCut,algo,jet_type,date.GetDate()),"RECREATE");
  

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
    hNumeratorPP_R3[j]->SetAxisRange(0.4,1.6,"Y");

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
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPP_Its_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_akR3%s_%d.pdf",unfoldingCut,jet_type,date.GetDate()),"RECREATE");

  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R3 = new TCanvas("cRAA_JEC_sys_R3","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_JEC_sys_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R3->cd(nbins_cent-i);
    RAA_JEC_bayesian_R3[i]->Divide(RAA_JEC_bayesian_R3[i], RAA_bayesian_R3[i], 1,1, "B");
    RAA_JEC_bayesian_R3[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R3[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R3[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_JEC_bayesian_R3[i]->Fit("fPol","");

    RAA_JEC_bayesian_R3[i]->Draw("p");
    checkMaximumSys(systematics_R3.hSysJEC[i], functionHist(fPol, systematics_R3.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_JEC_sys_R3->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R3",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_JEC_sys_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR3%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()),"RECREATE");

    // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R3 = new TCanvas("cRAA_Smear_sys_R3","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Smear_sys_R3,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R3->cd(nbins_cent-i);
    RAA_Smear_bayesian_R3[i]->Divide(RAA_Smear_bayesian_R3[i], RAA_bayesian_R3[i], 1,1, "B");
    RAA_Smear_bayesian_R3[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R3[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R3[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_Smear_bayesian_R3[i]->Fit("fPol","");
    RAA_Smear_bayesian_R3[i]->Draw("p");
    checkMaximumSys(systematics_R3.hSysSmear[i], functionHist(fPol, systematics_R3.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_Smear_sys_R3->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R3",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_Smear_sys_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R3%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()));


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
  cSys_R3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_total_systematics_unfoldingCut%dGeV_ak%sR3%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()),"RECREATE");



  TCanvas * cPbPb_Its_R4 = new TCanvas("cPbPb_Its_R4","",1000,800);
  makeMultiPanelCanvasWithGap(cPbPb_Its_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  TLegend *PbPb_itersys_R4 = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its_R4->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator_R4[j][i]->SetMarkerStyle(33);
      hNumerator_R4[j][i]->SetMarkerColor(j+1);
      hNumerator_R4[j][i]->SetAxisRange(unfoldingCut,299,"X");
      hNumerator_R4[j][i]->SetAxisRange(0.4,1.6,"Y");

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
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPbPb_Its_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%sR4%s_%d.pdf",unfoldingCut,algo,jet_type,date.GetDate()),"RECREATE");
  

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
    hNumeratorPP_R4[j]->SetAxisRange(0.4,1.6,"Y");

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
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPP_Its_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_akR4%s_%d.pdf",unfoldingCut,jet_type,date.GetDate()),"RECREATE");

  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys_R4 = new TCanvas("cRAA_JEC_sys_R4","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_JEC_sys_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys_R4->cd(nbins_cent-i);
    RAA_JEC_bayesian_R4[i]->Divide(RAA_JEC_bayesian_R4[i], RAA_bayesian_R4[i], 1,1, "B");
    RAA_JEC_bayesian_R4[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian_R4[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian_R4[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_JEC_bayesian_R4[i]->Fit("fPol","");

    RAA_JEC_bayesian_R4[i]->Draw("p");
    checkMaximumSys(systematics_R4.hSysJEC[i], functionHist(fPol, systematics_R4.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_JEC_sys_R4->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R4",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_JEC_sys_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_JEC_systematics_unfoldingCut%dGeV_ak%sR4%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()),"RECREATE");

    // plot for the Smear sys

  TCanvas * cRAA_Smear_sys_R4 = new TCanvas("cRAA_Smear_sys_R4","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Smear_sys_R4,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys_R4->cd(nbins_cent-i);
    RAA_Smear_bayesian_R4[i]->Divide(RAA_Smear_bayesian_R4[i], RAA_bayesian_R4[i], 1,1, "B");
    RAA_Smear_bayesian_R4[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian_R4[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian_R4[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_Smear_bayesian_R4[i]->Fit("fPol","");
    RAA_Smear_bayesian_R4[i]->Draw("p");
    checkMaximumSys(systematics_R4.hSysSmear[i], functionHist(fPol, systematics_R4.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_Smear_sys_R4->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R4",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_Smear_sys_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s_R4%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()));


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
  cSys_R4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_total_systematics_unfoldingCut%dGeV_ak%sR4%s_%d.pdf",unfoldingCut, algo, jet_type, date.GetDate()),"RECREATE");

  
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

    for(int j = 1; j<RAA_measured_R2[i]->GetNbinsX(); ++j){
      
      RAA_R2_Bayes[i]->SetBinError(j, RAA_measured_R2[i]->GetBinError(j));
      RAA_R3_Bayes[i]->SetBinError(j, RAA_measured_R3[i]->GetBinError(j));
      RAA_R4_Bayes[i]->SetBinError(j, RAA_measured_R4[i]->GetBinError(j));
      
    }

  }

  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart->SetBinContent(hRAA_R2_npart->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinContent(RAA_R2_Bayes[i]->FindBin(100)));
    hRAA_R2_npart->SetBinError(hRAA_R2_npart->FindBin(npart[i])-3, RAA_R2_Bayes[i]->GetBinError(RAA_R2_Bayes[i]->FindBin(100)));
    
    hRAA_R3_npart->SetBinContent(hRAA_R3_npart->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(100)));
    hRAA_R3_npart->SetBinError(hRAA_R3_npart->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinError(RAA_R3_Bayes[i]->FindBin(100)));
    
    hRAA_R4_npart->SetBinContent(hRAA_R4_npart->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinContent(RAA_R4_Bayes[i]->FindBin(100)));    
    hRAA_R4_npart->SetBinError(hRAA_R4_npart->FindBin(npart[i])+3, RAA_R4_Bayes[i]->GetBinError(RAA_R4_Bayes[i]->FindBin(100)));    
  }


  TCanvas * cRAA_npart = new TCanvas("cRAA_npart","",600,400);
  cRAA_npart->SetGridy();
  cRAA_npart->SetGridx();

  hRAA_R2_npart->SetTitle(" ");
  hRAA_R2_npart->SetXTitle(" N_{part} ");
  hRAA_R2_npart->SetYTitle(" R_{AA} ");
  hRAA_R2_npart->SetMarkerColor(kRed);
  hRAA_R2_npart->SetLineColor(kRed);
  hRAA_R2_npart->SetMarkerStyle(20);
  hRAA_R2_npart->SetAxisRange(0.0,1.1, "Y");
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
  drawText("Jet ID cut, |#eta|<2",0.15,0.3,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.15,0.2,16);
  drawText("97 < Jet p_{T} < 114", 0.15,0.8,16);
  
  TLegend * npart1 = myLegend(0.7,0.7,0.9,0.9);
  npart1->AddEntry(hRAA_R2_npart,"R=0.2", "pl");
  npart1->AddEntry(hRAA_R3_npart,"R=0.3", "pl");
  npart1->AddEntry(hRAA_R4_npart,"R=0.4", "pl");
  npart1->SetTextSize(0.04);
  npart1->Draw();
 
  for(int i = 0; i<nbins_cent; ++i){
    systematics_R3.calcTotalSys(i);
    systematics_R3.DrawNpartSys(RAA_R3_Bayes[i]->GetBinContent(RAA_R3_Bayes[i]->FindBin(100)),i,npart[i]);
  }

  hRAA_R2_npart->Draw("same E1");
  hRAA_R3_npart->Draw("same E1");
  hRAA_R4_npart->Draw("same E1");

  
  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  
  cRAA_npart->SaveAs(Form("../../Plots/Final_paper_plots_RAA_npart_%d.pdf",date.GetDate()),"RECREATE");
  
  
  // plotting the RAA here: 


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

    if(i==0){
      for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(100); ++j){
	RAA_R2_Bayes[i]->SetBinContent(j,0);
	RAA_R3_Bayes[i]->SetBinContent(j,0);
	RAA_R4_Bayes[i]->SetBinContent(j,0);
	RAA_R2_Bayes[i]->SetBinError(j,0);
	RAA_R3_Bayes[i]->SetBinError(j,0);
	RAA_R4_Bayes[i]->SetBinError(j,0);
      }
    }
    if(i==5 || i==4 || i==3){
      for(int j = RAA_R3_Bayes[i]->FindBin(240); j<RAA_R3_Bayes[i]->FindBin(300); ++j){
	RAA_R2_Bayes[i]->SetBinContent(j,0);
	RAA_R3_Bayes[i]->SetBinContent(j,0);
	RAA_R4_Bayes[i]->SetBinContent(j,0);
	RAA_R2_Bayes[i]->SetBinError(j,0);
	RAA_R3_Bayes[i]->SetBinError(j,0);
	RAA_R4_Bayes[i]->SetBinError(j,0);
      }
    }
    if(i==1 || i==2){
      for(int j = RAA_R3_Bayes[i]->FindBin(unfoldingCut); j<RAA_R3_Bayes[i]->FindBin(unfoldingCut+30); ++j){
	RAA_R2_Bayes[i]->SetBinContent(j,0);
	RAA_R3_Bayes[i]->SetBinContent(j,0);
	RAA_R4_Bayes[i]->SetBinContent(j,0);
	RAA_R2_Bayes[i]->SetBinError(j,0);
	RAA_R3_Bayes[i]->SetBinError(j,0);
	RAA_R4_Bayes[i]->SetBinError(j,0);
      }
    }

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
    RAA_R3_Bayes[i]->Draw("same E1");
    RAA_R4_Bayes[i]->Draw("same E1");

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
  drawText("Jet ID cut, |#eta|<2",0.1,0.3,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  cRAA->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.3,16);
  drawText("Pile up rejection cut applied",0.1,0.2,16);

  cRAA->SaveAs(Form("../../Plots/Final_paper_plots_RAA_%d.pdf",date.GetDate()),"RECREATE");
    
  
  
}
