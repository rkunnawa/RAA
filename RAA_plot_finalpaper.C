// Raghav Kunnawalkam Elayavalli
// April 14 2014
// Rutgers
// raghav.k.e at CERN dot CH

//
// Macro to plot the final paper plots.  
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

void RAA_plot_finalpaper(Int_t unfoldingCut = 30 , char *algo = "Pu", char *jet_type = "PF"){
    
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
  
  // plot 1 - spectra plot showing pp and 6 different centrality classes PbPb spectra
  //        - have a 3 panel plot for the different radii, with each of them scaled by two orders of magnitude 

  // first we need to scale the MC to the level of Data:
  // PbPb Data scaling:
  //   uPbPb_Bayes[i]->Scale(1./deltaEta);// delta eta
  //   //uPbPb_Bayes[i]->Scale(1./145.156/1e6);// Jet 80 luminosity
  //   //uPbPb_Bayes[i]->Scale(1./1.1153/1e6);// equivalent no of minbias events 
  //   uPbPb_Bayes[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
  //   //uPbPb_Bayes[i]->Scale(1./145.156);
  //   //uPbPb_Bayes[i]->Scale(1./161.939);
  //   uPbPb_Bayes[i]->Scale(1./(7.65*1e6));
  //   uPbPb_Bayes[i]->Scale(64.*1e9/(ncoll[i]*1e3));
  //   uPbPb_Bayes[i] = (TH1F*)uPbPb_Bayes[i]->Rebin(nbins_pt,Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i),boundaries_pt);
  //   divideBinWidth(uPbPb_Bayes[i]);
  //   uPbPb_Bayes[i]->Write();
  //   So finally PbPb is in 

  // PbPb MC scaling: is already in sigma (mb) / (dEta dpT)
  //   mPbPb_Reco[i]->Scale(1./deltaEta);// delta eta
  //   mPbPb_Reco[i] = (TH1F*)mPbPb_Reco[i]->Rebin(nbins_pt,Form("PbPb_Reco_spectra_refpt_cent%d",i),boundaries_pt);
  //   divideBinWidth(mPbPb_Reco[i]);
  //   mPbPb_Reco[i]->Write();

  // take MC to nano barns from milli barns 
  mPP_R2->Scale(1e6);
  mPP_R3->Scale(1e6);
  mPP_R4->Scale(1e6);

  // divide by delta eta = 4. 
  mPP_R2->Scale(1./4);
  mPP_R3->Scale(1./4);
  mPP_R4->Scale(1./4);

  for(int i = 0; i<nbins_cent; ++i){

    mPbPb_R2[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    mPbPb_R2[i]->Scale(64.*1e9/(ncoll[i]*1e3));
    mPbPb_R2[i]->Scale(1./(7.65));

    mPbPb_R3[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    mPbPb_R3[i]->Scale(64.*1e9/(ncoll[i]*1e3));
    mPbPb_R3[i]->Scale(1./(7.65));

    mPbPb_R4[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    mPbPb_R4[i]->Scale(64.*1e9/(ncoll[i]*1e3));
    mPbPb_R4[i]->Scale(1./(7.65));
    
  }
  
  Double_t ScaleFactor[nbins_cent+2] = {1, 1e2, 1e4, 1e6, 1e8, 1e10, 1e12, 1e14};  
  
  TCanvas * cSpectra = new TCanvas("cSpectra","",1200,1000);
  //makeMultiPanelCanvas(cSpectra,3,1,0.0,0.0,0.2,0.15,0.07);
  
#if 0
  cSpectra->cd(1);
  cSpectra->cd(1)->SetLogy();
  cSpectra->cd(1)->SetLogx();

  uPP_R2_Bayes->Scale(ScaleFactor[0]);
  uPP_R2_Bayes->SetMarkerStyle(20);
  uPP_R2_Bayes->SetMarkerColor(kBlack);
  makeHistTitle(uPP_R2_Bayes," "," Jet p_{T} (GeV/c)","cross section");
  uPP_R2_Bayes->SetAxisRange(unfoldingCut, 299, "X");
  uPP_R2_Bayes->SetAxisRange(1e-5, 1e8, "Y");
  uPP_R2_Bayes->Draw();
  
  // draw the MC
  mPP_R2->Scale(ScaleFactor[0]);
  mPP_R2->SetLineStyle(2);
  mPP_R2->SetLineColor(kBlack);
  mPP_R2->Draw("same C");
  
  for(int i = 0; i<nbins_cent; ++i){
    uPbPb_R2_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R2_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R2_Bayes[i]->SetMarkerColor(kBlack);
    uPbPb_R2_Bayes[i]->Draw("same");

    mPbPb_R2[i]->Scale(ScaleFactor[i+1]);
    mPbPb_R2[i]->SetLineStyle(2);
    mPbPb_R2[i]->SetLineColor(kBlack);
    mPbPb_R2[i]->Draw("same C");
    
  }
#endif

  //cSpectra->cd(2);
  cSpectra->SetLogy();
  cSpectra->SetGridy();
  cSpectra->SetLogx();

  uPP_R3_Bayes->Scale(ScaleFactor[0]);
  uPP_R3_Bayes->SetMarkerStyle(20);
  uPP_R3_Bayes->SetMarkerColor(kBlack);
  makeHistTitle(uPP_R3_Bayes," "," Jet p_{T} (GeV/c)","#frac{d #sigma}{T_{AA} dp_{T} d#eta} nb");
  uPP_R3_Bayes->SetAxisRange(unfoldingCut, 299, "X");
  uPP_R3_Bayes->SetAxisRange(1e-4, 1e14, "Y");
  uPP_R3_Bayes->Draw();

  // draw the MC
  mPP_R3->Scale(ScaleFactor[0]);
  //mPP_R3->SetLineStyle(2);
  mPP_R3->SetLineColor(kBlack);
  mPP_R3->SetAxisRange(unfoldingCut, 299, "X");
  mPP_R3->Draw("same Lhist");
  
  for(int i = 0; i<nbins_cent; ++i){
    uPbPb_R3_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R3_Bayes[i]->SetMarkerStyle(33);
    uPbPb_R3_Bayes[i]->SetMarkerColor(kRed);
    uPbPb_R3_Bayes[i]->SetAxisRange(unfoldingCut, 299, "X");
    uPbPb_R3_Bayes[i]->Draw("same");

    mPbPb_R3[i]->Scale(ScaleFactor[i+2]);
    //mPbPb_R3[i]->SetLineStyle(2);
    mPbPb_R3[i]->SetLineColor(kRed);
    mPbPb_R3[i]->SetAxisRange(unfoldingCut, 299, "X");
    mPbPb_R3[i]->Draw("same Lhist");
  }

  TLegend * leg1 = myLegend(0.15,0.1,0.25,0.2);
  leg1->AddEntry(uPP_R3_Bayes,"PP Data","pl");
  leg1->AddEntry(mPP_R3,"PYTHIA","pl");
  leg1->SetTextSize(0.02);
  leg1->Draw();


  TLegend * leg2 = myLegend(0.75,0.8,0.85,0.9);
  leg2->AddEntry(uPbPb_R3_Bayes[0],"PbPb Data","pl");
  leg2->AddEntry(mPbPb_R3[0],"PYTHIA+HYDJET","pl");
  leg2->SetTextSize(0.02);
  leg2->Draw();
  
  
  drawText("R=0.3, anti k_{T} PF Jets", 0.15,0.2,16);
  drawText("R=0.3, anti k_{T} Pu PF Jets", 0.75,0.78,16);

#if 0
  cSpectra->cd(3);
  cSpectra->cd(3)->SetLogy();
  cSpectra->cd(3)->SetLogx();

  uPP_R4_Bayes->Scale(ScaleFactor[0]);
  uPP_R4_Bayes->SetMarkerStyle(20);
  uPP_R4_Bayes->SetMarkerColor(kBlack);
  makeHistTitle(uPP_R4_Bayes," "," Jet p_{T} (GeV/c)","cross section");
  uPP_R4_Bayes->SetAxisRange(unfoldingCut, 299, "X");
  uPP_R4_Bayes->SetAxisRange(1e-5, 1e8, "Y");
  uPP_R4_Bayes->Draw();

  // draw the MC
  mPP_R4->Scale(ScaleFactor[0]);
  mPP_R4->SetLineStyle(2);
  mPP_R4->SetLineColor(kBlack);
  mPP_R4->Draw("same C");
  
  uPP_R4_Bayes->Draw();
  for(int i = 0; i<nbins_cent; ++i){
    uPbPb_R4_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R4_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R4_Bayes[i]->SetMarkerColor(kBlack);
    uPbPb_R4_Bayes[i]->Draw("same");
    
    mPbPb_R4[i]->Scale(ScaleFactor[i]);
    mPbPb_R4[i]->SetLineStyle(2);
    mPbPb_R4[i]->SetLineColor(kBlack);
    mPbPb_R4[i]->Draw("same C");
  }
#endif

  putCMSPrel();
  putPbPbLumi();
  putPPLumi();
  
  //drawText("pp", 0.7,0.15,16);
  drawText("0-5% x 10^{2}", 0.8,0.20,16);
  drawText("5-10% x 10^{4}", 0.8,0.28,16);
  drawText("10-30% x 10^{6}", 0.8,0.38,16);
  drawText("30-50% x 10^{8}", 0.8,0.47,16);
  drawText("50-70% x 10^{10}", 0.8,0.54,16);
  drawText("70-90% x 10^{12}", 0.8,0.63,16);

  cSpectra->SaveAs(Form("../../Plots/Final_paper_plots_spectra_%d.pdf",date.GetDate()),"RECREATE");

}
