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

void RAA_plot_Systematics(int radius = 3, char *algo = "Pu", char *jet_type = "PF", int unfoldingCut = 60){

  
  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  TFile *fin; 
  fin= TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_fullMC_n20_eta_p20_%dGeVCut_ak%s_20150414.root",radius,unfoldingCut,jet_type));

  // declare the systematics
  SysData systematics;
  // ncoll uncertainty
  prepareNcollUnc(nbins_pt, 300.);
  
  
  // get the necessary histograms
  const int Iterations = 6; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_BayesianIter[nbins_cent+1][Iterations];
  TH1F *uPP_BayesianIter[Iterations];

  for(int j = 0;j<Iterations;j++){
    for(int i = 0;i<nbins_cent;++i)
      uPbPb_BayesianIter[i][j] = (TH1F*)fin->Get(Form("uPbPb_BayesianIter%d_cent%d",j+1,i));
    uPP_BayesianIter[j] = (TH1F*)fin->Get(Form("uPP_BayesianIter%d",j+1));
  }

  //TH1F * uPbPb_Bayes[nbins_cent], * uPbPb_JEC_Bayes[nbins_cent], * uPbPb_Smear_Bayes[nbins_cent];
  //TH1F * uPP_Bayes, dPP_Bayes, dPP_BinByBin;
  TH1F *RAA_measured[nbins_cent+1];
  TH1F *RAA_binbybin[nbins_cent+1];
  TH1F *RAA_bayesian[nbins_cent+1];
  TH1F *RAA_JEC_bayesian[nbins_cent+1];
  TH1F *RAA_Smear_bayesian[nbins_cent+1];
  
  for(int i = 0; i<nbins_cent; ++i){

    RAA_bayesian[i] = (TH1F*)fin->Get(Form("RAA_bayesian_cent%d",i));
    RAA_JEC_bayesian[i] = (TH1F*)fin->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Smear_bayesian[i] = (TH1F*)fin->Get(Form("RAA_bayesian_cent%d",i));
    RAA_binbybin[i] = (TH1F*)fin->Get(Form("RAA_binbybin_cent%d",i));
    RAA_measured[i] = (TH1F*)fin->Get(Form("RAA_measured_cent%d",i));
    // uPbPb_Bayes[i] = (TH1F*)fin->Get(Form("PbPb_bayesian_unfolded_spectra_combined__cent%d",i));
    // uPbPb_Bayes[i]->Print("base");
  }
  
  // make the ratios for each iteration with respect to iteration 4
  TH1F * hDenominator[nbins_cent];
  TH1F * hNumerator[Iterations][nbins_cent];
  TH1F * hDenominatorPP;
  TH1F * hNumeratorPP[Iterations];
  
  for(int i = 0; i<nbins_cent; ++i){
    hDenominator[i] = (TH1F*)uPbPb_BayesianIter[i][BayesIter-1]->Clone(Form("Denominator_PbPb_iter%d_cent%d",BayesIter,i));
    for(int j = 0; j<Iterations; ++j){
      hNumerator[j][i] = (TH1F*)uPbPb_BayesianIter[i][j]->Clone(Form("Numerator_PbPb_iter%d_cent%d",j+1,i));
      hNumerator[j][i]->Divide(hDenominator[i]);
    }
  }

  hDenominatorPP = (TH1F*)uPP_BayesianIter[BayesIter-1]->Clone(Form("Denominator_PP_iter%d",BayesIter));
  for(int j = 0; j<Iterations; ++j){
    hNumeratorPP[j] = (TH1F*)uPP_BayesianIter[j]->Clone(Form("Numerator_PP_iter%d",j+1));
    hNumeratorPP[j]->Divide(hDenominatorPP);
  }  

  //make the canvas
  
  TCanvas * cPbPb_Its = new TCanvas("cPbPb_Its","",1000,800);
  makeMultiPanelCanvasWithGap(cPbPb_Its,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  // line at 1
  TLine *linePbPb_iter = new TLine(unfoldingCut,1,299,1);
  linePbPb_iter->SetLineStyle(2);
  linePbPb_iter->SetLineWidth(2);
  
  TLegend *PbPb_itersys = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cPbPb_Its->cd(nbins_cent-i);

    for(int j = 1;j<Iterations;j++){
      hNumerator[j][i]->SetMarkerStyle(33);
      hNumerator[j][i]->SetMarkerColor(j+1);
      hNumerator[j][i]->SetAxisRange(unfoldingCut,299,"X");
      hNumerator[j][i]->SetAxisRange(0.4,1.6,"Y");

      if(j==1){
	makeHistTitle(hNumerator[j][i]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
	hNumerator[j][i]->Draw();
      }
      hNumerator[j][i]->Draw("same");
	
      if(i==0) PbPb_itersys->AddEntry(hNumerator[j][i],Form("Iteration %d",j+1),"pl");
      checkMaximumSys(systematics.hSysIter[i], hNumerator[j][i], 0, 1.05);
      
    }

    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    linePbPb_iter->Draw();
    drawEnvelope(systematics.hSysIter[i], "hist same");

  }
  PbPb_itersys->Draw();

  cPbPb_Its->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPbPb_Its->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%s%d%s_%d.pdf",unfoldingCut,algo,radius,jet_type,date.GetDate()),"RECREATE");
  

  TCanvas * cPP_Its = new TCanvas("cPP_Its","",600,400);
  // line at 1
  TLine *linePP_iter = new TLine(unfoldingCut,1,299,1);
  linePP_iter->SetLineStyle(2);
  linePP_iter->SetLineWidth(2);
  
  TLegend *PP_itersys = myLegend(0.53,0.65,0.85,0.9);
  
  for(int j = 1;j<Iterations;j++){
    hNumeratorPP[j]->SetMarkerStyle(33);
    hNumeratorPP[j]->SetMarkerColor(j+1);
    hNumeratorPP[j]->SetAxisRange(unfoldingCut,299,"X");
    hNumeratorPP[j]->SetAxisRange(0.4,1.6,"Y");

    if(j==1){
      makeHistTitle(hNumeratorPP[j]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
      hNumeratorPP[j]->Draw();
    }
    hNumeratorPP[j]->Draw("same");

    checkMaximumSys(systematics.hSysIter[nbins_cent], hNumeratorPP[j],0, 1.05);
    PP_itersys->AddEntry(hNumeratorPP[j],Form("Iteration %d",j+1),"pl");

  }

  linePP_iter->Draw();  
  PP_itersys->Draw();
  drawEnvelope(systematics.hSysIter[nbins_cent],"hist same");

  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPP_Its->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_ak%d%s_%d.pdf",unfoldingCut,radius,jet_type,date.GetDate()),"RECREATE");

  // define the Jet ID efficiency as a function of jet pT: taken from the above plots: 
  
  // plot for the JEC sys

  TCanvas * cRAA_JEC_sys = new TCanvas("cRAA_JEC_sys","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_JEC_sys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  TF1 * fPol = new TF1("fPol","[0]+[1]*x");
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_JEC_sys->cd(nbins_cent-i);
    RAA_JEC_bayesian[i]->Divide(RAA_JEC_bayesian[i], RAA_bayesian[i], 1,1, "B");
    RAA_JEC_bayesian[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_JEC_bayesian[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_JEC_bayesian[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_JEC_bayesian[i]->Fit("fPol","");

    RAA_JEC_bayesian[i]->Draw("p");
    checkMaximumSys(systematics.hSysJEC[i], functionHist(fPol, systematics.hSysJEC[i], Form("hist_sysJEC_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_JEC_sys->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_JEC_sys->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_JEC_systematics_unfoldingCut%dGeV_ak%s%d%s_%d.pdf",unfoldingCut, algo, radius, jet_type, date.GetDate()),"RECREATE");

    // plot for the Smear sys

  TCanvas * cRAA_Smear_sys = new TCanvas("cRAA_Smear_sys","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Smear_sys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRAA_Smear_sys->cd(nbins_cent-i);
    RAA_Smear_bayesian[i]->Divide(RAA_Smear_bayesian[i], RAA_bayesian[i], 1,1, "B");
    RAA_Smear_bayesian[i]->SetAxisRange(0.61, 1.39, "Y");
    makeHistTitle(RAA_Smear_bayesian[i], "", "Jet p_{T} (GeV/c)","Ratio",2);
    RAA_Smear_bayesian[i]->SetAxisRange(unfoldingCut, 299, "X");
    RAA_Smear_bayesian[i]->Fit("fPol","");
    RAA_Smear_bayesian[i]->Draw("p");
    checkMaximumSys(systematics.hSysSmear[i], functionHist(fPol, systematics.hSysSmear[i], Form("hist_sysSmear_cent%d",i)));

    linePbPb_iter->Draw();
  }

  cRAA_Smear_sys->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cRAA_Smear_sys->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_Smear_systematics_unfoldingCut%dGeV_ak%s%d%s_%d.pdf",unfoldingCut, algo, radius, jet_type, date.GetDate()));


  // draw the total systematics
  TCanvas * cSys = new TCanvas("cSys","Total Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0; i<nbins_cent; ++i){

    cSys->cd(nbins_cent-i);
    cSys->cd(nbins_cent-i)->SetGridy();
    cSys->cd(nbins_cent-i)->SetGridx();
    systematics.DrawComponent(i);

    //drawPanelLabel(i);
    TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    title_->AddEntry(RAA_bayesian[i],Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    title_->SetTextSize(0.06);
    title_->Draw();
  
  }
  cSys->cd(1);
  putCMSPrel();
  linePbPb_iter->Draw();
  cSys->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_total_systematics_unfoldingCut%dGeV_ak%s%d%s_%d.pdf",unfoldingCut, algo, radius, jet_type, date.GetDate()),"RECREATE");
  
  
}
