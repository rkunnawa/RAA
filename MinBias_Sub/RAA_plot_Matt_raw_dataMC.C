// Raghav Kunnawalkam Elayavalli
// May 6th 2015
// Rutgers

// contact: raghav dot k dot e at CERN dot ch 


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


static const int nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};


  
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


void RAA_plot_Matt_raw_dataMC(char* etaWidth = (char*)"20_eta_20", Int_t radius = 2){

  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  TDatime date;

  // char * etaWidth = (char*)Form("%d_eta_%d",etaLow, etaHigh);
  cout<<"etaWidth = "<<etaWidth<<endl;

  Float_t etaLow = 0.0;
  Float_t etaHigh = 1.0; 
  
  if(etaWidth == "20_eta_20") {
    etaLow = 0.0;
    etaHigh = 2.0;
  }
  if(etaWidth == "10_eta_18") {
    etaLow = 1.0;
    etaHigh = 1.8;
  }
  if(etaWidth == "10_eta_10") {
    etaLow = 0.0;
    etaHigh = 1.0;
  }
  
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  int unfoldingCut = 40; 

  if(radius == 2) unfoldingCut = 30;
  if(radius == 3) unfoldingCut = 40;
  if(radius == 4) unfoldingCut = 50;
  
  // Make Matt's requested plots about the unfolding prior distributions.
  // need to make ratio histograms for RAA, and PbPb spectra

  TFile * fNoSmear = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_without80FakeRemoval_unfold_mcclosure_oppside_trgMC_noSmear_%s_%dGeVCut_akPF_20150521.root",radius,etaWidth,unfoldingCut));
  TFile * fGenSmear = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_without80FakeRemoval_unfold_mcclosure_oppside_trgMC_GenSmear_%s_%dGeVCut_akPF_20150519.root",radius,etaWidth,unfoldingCut));
  TFile * fgen2pSmear = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_without80FakeRemoval_unfold_mcclosure_oppside_trgMC_gen2pSmear_%s_%dGeVCut_akPF_20150519.root",radius,etaWidth,unfoldingCut));
  TFile * fRecoSmear = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_without80FakeRemoval_unfold_mcclosure_oppside_trgMC_RecoSmear_%s_%dGeVCut_akPF_20150519.root",radius,etaWidth,unfoldingCut));
  TFile * fBothSmear = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_without80FakeRemoval_unfold_mcclosure_oppside_trgMC_BothSmear_%s_%dGeVCut_akPF_20150519.root",radius,etaWidth,unfoldingCut));
  
  TH1F * RAA_Bayesian[nbins_cent], * RAA_Bayesian_gen2pSmear[nbins_cent], * RAA_Bayesian_GenSmear[nbins_cent], * RAA_Bayesian_RecoSmear[nbins_cent], * RAA_Bayesian_BothSmear[nbins_cent];
  TH1F * hPbPb_Bayesian[nbins_cent], * hPbPb_Bayesian_gen2pSmear[nbins_cent], * hPbPb_Bayesian_GenSmear[nbins_cent], * hPbPb_Bayesian_RecoSmear[nbins_cent], * hPbPb_Bayesian_BothSmear[nbins_cent];

  TH1F * RAA_Sys_GenSmear[nbins_cent], * RAA_Sys_gen2pSmear[nbins_cent], * RAA_Sys_RecoSmear[nbins_cent],* RAA_Sys_BothSmear[nbins_cent];
  TH1F * hPbPb_Sys_GenSmear[nbins_cent], * hPbPb_Sys_gen2pSmear[nbins_cent], * hPbPb_Sys_RecoSmear[nbins_cent],* hPbPb_Sys_BothSmear[nbins_cent];
  
  // get the input files:
  TFile * fin = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_ntuple_PbPb_data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius));

  // get the histograms
  TH1F * hData[nbins_cent], * hData_unf[nbins_cent], * hMC[nbins_cent], * hMC_gen[nbins_cent];
  TH1F * hRatio[nbins_cent], * hRatio_unf[nbins_cent], * hRatio_gen[nbins_cent], * hRatio_unf_gen[nbins_cent];

  for(int i = 0; i<nbins_cent; ++i){

    RAA_Bayesian[i] = (TH1F*)fNoSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian[i]->Print("base");
    RAA_Bayesian_gen2pSmear[i] = (TH1F*)fgen2pSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian_gen2pSmear[i]->Print("base");
    RAA_Bayesian_GenSmear[i] = (TH1F*)fGenSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian_GenSmear[i]->Print("base");
    RAA_Bayesian_RecoSmear[i] = (TH1F*)fRecoSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian_RecoSmear[i]->Print("base");
    RAA_Bayesian_BothSmear[i] = (TH1F*)fBothSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian_BothSmear[i]->Print("base");

    RAA_Sys_GenSmear[i] = (TH1F*)RAA_Bayesian_GenSmear[i]->Clone(Form("RAA_Sys_GenSmear_cent%d",i));
    RAA_Sys_GenSmear[i]->Divide(RAA_Bayesian[i]);
    RAA_Sys_gen2pSmear[i] = (TH1F*)RAA_Bayesian_gen2pSmear[i]->Clone(Form("RAA_Sys_gen2pSmear_cent%d",i));
    RAA_Sys_gen2pSmear[i]->Divide(RAA_Bayesian[i]);
    RAA_Sys_RecoSmear[i] = (TH1F*)RAA_Bayesian_RecoSmear[i]->Clone(Form("RAA_Sys_RecoSmear_cent%d",i));
    RAA_Sys_RecoSmear[i]->Divide(RAA_Bayesian[i]);
    RAA_Sys_BothSmear[i] = (TH1F*)RAA_Bayesian_BothSmear[i]->Clone(Form("RAA_Sys_BothSmear_cent%d",i));
    RAA_Sys_BothSmear[i]->Divide(RAA_Bayesian[i]);
    
    hPbPb_Bayesian[i] = (TH1F*)fNoSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian[i]->Print("base");
    hPbPb_Bayesian_GenSmear[i] = (TH1F*)fGenSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian_GenSmear[i]->Print("base");
    hPbPb_Bayesian_gen2pSmear[i] = (TH1F*)fgen2pSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian_gen2pSmear[i]->Print("base");
    hPbPb_Bayesian_RecoSmear[i] = (TH1F*)fRecoSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian_RecoSmear[i]->Print("base");
    hPbPb_Bayesian_BothSmear[i] = (TH1F*)fBothSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian_BothSmear[i]->Print("base");

    hPbPb_Sys_GenSmear[i] = (TH1F*)hPbPb_Bayesian_GenSmear[i]->Clone(Form("hPbPb_Sys_GenSmear_cent%d",i));
    hPbPb_Sys_GenSmear[i]->Divide(hPbPb_Bayesian[i]);
    hPbPb_Sys_gen2pSmear[i] = (TH1F*)hPbPb_Bayesian_gen2pSmear[i]->Clone(Form("hPbPb_Sys_gen2pSmear_cent%d",i));
    hPbPb_Sys_gen2pSmear[i]->Divide(hPbPb_Bayesian[i]);
    hPbPb_Sys_RecoSmear[i] = (TH1F*)hPbPb_Bayesian_RecoSmear[i]->Clone(Form("hPbPb_Sys_RecoSmear_cent%d",i));
    hPbPb_Sys_RecoSmear[i]->Divide(hPbPb_Bayesian[i]);
    hPbPb_Sys_BothSmear[i] = (TH1F*)hPbPb_Bayesian_BothSmear[i]->Clone(Form("hPbPb_Sys_BothSmear_cent%d",i));
    hPbPb_Sys_BothSmear[i]->Divide(hPbPb_Bayesian[i]);
    
    hData[i] = (TH1F*)fNoSmear->Get(Form("PbPb_measured_spectra_combined_cent%d",i));
    hData[i]->Print("base");
    hData[i]->Scale((0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    hData[i]->Scale(ncoll[i]/(64.*1e3));

    hData_unf[i] = (TH1F*)fNoSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hData_unf[i]->Print("base");
    hData_unf[i]->Scale((0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    hData_unf[i]->Scale(ncoll[i]/(64.*1e3));
    
    hMC[i] = (TH1F*)fNoSmear->Get(Form("PbPb_Reco_spectra_jtpt_cent%d",i));
    hMC[i]->Print("base");

    hMC_gen[i] = (TH1F*)fNoSmear->Get(Form("PbPb_Gen_spectra_refpt_cent%d",i));
    hMC_gen[i]->Print("base");
    
    hRatio[i] = (TH1F*)hData[i]->Clone(Form("Data_MC_jtpT_ratio_cent%d",i));
    hRatio[i]->Divide(hMC[i]);
    hRatio[i]->Print("base");

    hRatio_unf[i] = (TH1F*)hData_unf[i]->Clone(Form("Data_unf_MC_jtpT_ratio_cent%d",i));
    hRatio_unf[i]->Divide(hMC[i]);
    hRatio_unf[i]->Print("base");

    hRatio_gen[i] = (TH1F*)hData[i]->Clone(Form("Data_MC_refpT_ratio_cent%d",i));
    hRatio_gen[i]->Divide(hMC_gen[i]);
    hRatio_gen[i]->Print("base");

    hRatio_unf_gen[i] = (TH1F*)hData_unf[i]->Clone(Form("Data_unf_MC_refpT_ratio_cent%d",i));
    hRatio_unf_gen[i]->Divide(hMC_gen[i]);
    hRatio_unf_gen[i]->Print("base");
    
    for(int j = 1; j<hPbPb_Sys_GenSmear[i]->GetNbinsX(); ++j){

      RAA_Sys_BothSmear[i]->SetBinError(j,0);
      RAA_Sys_GenSmear[i]->SetBinError(j,0);
      RAA_Sys_gen2pSmear[i]->SetBinError(j,0);
      RAA_Sys_RecoSmear[i]->SetBinError(j,0);
      hPbPb_Sys_GenSmear[i]->SetBinError(j,0);
      hPbPb_Sys_gen2pSmear[i]->SetBinError(j,0);
      hPbPb_Sys_RecoSmear[i]->SetBinError(j,0);
      hPbPb_Sys_BothSmear[i]->SetBinError(j,0);

    }
    
  }

  TLine * line = new TLine(65,1,299,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
#if 0
  // make the bayesian unfolding prior systematics plot
  TCanvas * cRAA_Sys = new TCanvas("cRAA_Sys","",1000,1000);
  makeMultiPanelCanvas(cRAA_Sys,3,2,0.0,0.0,0.2,0.15,0.07);

  TLegend * lRAA_Sys = myLegend(0.25,0.7,0.4,0.9);
  for(int i = 0; i<nbins_cent; ++i){

    cRAA_Sys->cd(nbins_cent-i);
  
    RAA_Sys_GenSmear[i]->SetMarkerStyle(24);
    RAA_Sys_GenSmear[i]->SetMarkerColor(kRed);
    RAA_Sys_GenSmear[i]->SetAxisRange(65,299,"X");
    RAA_Sys_GenSmear[i]->SetAxisRange(0.5,1.5,"Y");
    makeHistTitle(RAA_Sys_GenSmear[i]," ","Jet p_{T} (GeV/c)","RAA ratio smear/ no smear");
    RAA_Sys_GenSmear[i]->Draw("p");

    RAA_Sys_RecoSmear[i]->SetMarkerStyle(24);
    RAA_Sys_RecoSmear[i]->SetMarkerColor(kBlue);
    RAA_Sys_RecoSmear[i]->Draw("same p");

    RAA_Sys_gen2pSmear[i]->SetMarkerStyle(24);
    RAA_Sys_gen2pSmear[i]->SetMarkerColor(kOrange);
    RAA_Sys_gen2pSmear[i]->Draw("same p");

    RAA_Sys_BothSmear[i]->SetMarkerStyle(24);
    RAA_Sys_BothSmear[i]->SetMarkerColor(kGreen);
    RAA_Sys_BothSmear[i]->Draw("same p");

    drawText(Form("%2.0f-%.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    line->Draw();
    RAA_Sys_GenSmear[i]->Draw("same p");
    RAA_Sys_RecoSmear[i]->Draw("same p");
    RAA_Sys_gen2pSmear[i]->Draw("same p");
    RAA_Sys_BothSmear[i]->Draw("same p");

  }

  cRAA_Sys->cd(1);
  putCMSPrel();
  
  lRAA_Sys->AddEntry(RAA_Sys_GenSmear[0],"Gen pT smeared","pl");
  lRAA_Sys->AddEntry(RAA_Sys_gen2pSmear[0],"Gen pT 2% resolution smear","pl");
  lRAA_Sys->AddEntry(RAA_Sys_RecoSmear[0],"Reco pT smeared","pl");
  lRAA_Sys->AddEntry(RAA_Sys_BothSmear[0],"Gen and Reco pT smeared","pl");
  lRAA_Sys->SetTextSize(0.04);
  lRAA_Sys->Draw();
  drawText(Form("ak Pu R = 0.%d PF jets, %2.0f < |eta| < %2.0f",radius, etaLow, etaHigh), 0.25,0.2,15);

  cRAA_Sys->SaveAs(Form("May19/RAA_bayesian_unfolded_prior_smear_shift_systematics_%s_R%d_%d.pdf",etaWidth,radius,date.GetDate()),"RECREATE");
#endif
#if 0
  // plot for spectra
  TCanvas * chPbPb_Sys = new TCanvas("chPbPb_Sys","",1000,800);
  makeMultiPanelCanvasWithGap(chPbPb_Sys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend * lhPbPb_Sys = myLegend(0.2,0.6,0.4,0.8);
  for(int i = 0; i<nbins_cent; ++i){

    chPbPb_Sys->cd(nbins_cent-i);

    hPbPb_Sys_GenSmear[i]->SetMarkerStyle(33);
    hPbPb_Sys_GenSmear[i]->SetMarkerColor(kRed);
    hPbPb_Sys_GenSmear[i]->SetAxisRange(40,299,"X");
    hPbPb_Sys_GenSmear[i]->SetAxisRange(0,2,"Y");
    makeHistTitle(hPbPb_Sys_GenSmear[i]," ","Jet p_{T} (GeV/c)","PbPb spectra ratio smear/ no smear");
    hPbPb_Sys_GenSmear[i]->Draw();

    hPbPb_Sys_RecoSmear[i]->SetMarkerStyle(33);
    hPbPb_Sys_RecoSmear[i]->SetMarkerColor(kBlue);
    hPbPb_Sys_RecoSmear[i]->Draw("same");

    hPbPb_Sys_gen2pSmear[i]->SetMarkerStyle(33);
    hPbPb_Sys_gen2pSmear[i]->SetMarkerColor(kPink);
    hPbPb_Sys_gen2pSmear[i]->Draw("same");

    hPbPb_Sys_BothSmear[i]->SetMarkerStyle(33);
    hPbPb_Sys_BothSmear[i]->SetMarkerColor(kGreen);
    hPbPb_Sys_BothSmear[i]->Draw("same");

    drawText(Form("%2.0f-%.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    line->Draw();
    
  }

  chPbPb_Sys->cd(1);
  putCMSPrel();
  
  lhPbPb_Sys->AddEntry(hPbPb_Sys_GenSmear[0],"Gen pT smeared","pl");
  lhPbPb_Sys->AddEntry(hPbPb_Sys_gen2pSmear[0],"Gen pT 2% resolution smear","pl");
  lhPbPb_Sys->AddEntry(hPbPb_Sys_RecoSmear[0],"Reco pT smeared","pl");
  lhPbPb_Sys->AddEntry(hPbPb_Sys_BothSmear[0],"Gen and Reco pT smeared","pl");
  lhPbPb_Sys->SetTextSize(0.04);
  lhPbPb_Sys->Draw();
  drawText(Form("ak Pu R = 0.%d PF jets, %2.0f < |eta| < %2.0f",radius, etaLow, etaHigh), 0.2,0.2,15);

  chPbPb_Sys->SaveAs(Form("May19/PbPb_spectra_bayesian_unfolded_prior_smear_shift_systematics_%s_R%d_%d.pdf",etaWidth,radius,date.GetDate()),"RECREATE");

#endif
  // make the plots - 3x2 plots.
  TCanvas * cDataMC_raw = new TCanvas("cDataMC_raw","",1000,800);
  makeMultiPanelCanvas(cDataMC_raw,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0; i<nbins_cent; ++i){

    cDataMC_raw->cd(nbins_cent-i);
    cDataMC_raw->cd(nbins_cent-i)->SetLogy();

    hRatio[i]->SetMarkerStyle(33);
    hRatio[i]->SetMarkerColor(kBlack);
    makeHistTitle(hRatio[i]," ","Jet p_{T} (GeV/c)","Data/MC arbitraty scale");
    hRatio[i] = (TH1F*)hRatio[i]->Rebin(nbins_pt, Form("Data_MC_jtpT_ratio_cent%d",i), boundaries_pt);
    divideBinWidth(hRatio[i]);
    hRatio_unf[i] = (TH1F*)hRatio_unf[i]->Rebin(nbins_pt, Form("Data_unf_MC_jtpT_ratio_cent%d",i), boundaries_pt);
    divideBinWidth(hRatio_unf[i]);
    hRatio_gen[i] = (TH1F*)hRatio_gen[i]->Rebin(nbins_pt, Form("Data_MC_refpT_ratio_cent%d",i), boundaries_pt);
    divideBinWidth(hRatio_gen[i]);
    hRatio_unf_gen[i] = (TH1F*)hRatio_unf_gen[i]->Rebin(nbins_pt, Form("Data_unf_MC_refpT_ratio_cent%d",i), boundaries_pt);
    divideBinWidth(hRatio_unf_gen[i]);
    hRatio[i]->SetAxisRange(30,299,"X");
    //hRatio[i]->GetYaxis()->SetNdivisions(2,kFALSE);
    hRatio[i]->Draw(" p");

    hRatio_gen[i]->SetMarkerStyle(33);
    hRatio_gen[i]->SetMarkerColor(kRed);
    hRatio_gen[i]->Draw("same p");
      
    hRatio_unf[i]->SetMarkerStyle(33);
    hRatio_unf[i]->SetMarkerColor(kGreen);
    hRatio_unf[i]->Draw("same p");
    
    hRatio_unf_gen[i]->SetMarkerStyle(33);
    hRatio_unf_gen[i]->SetMarkerColor(kBlue);
    hRatio_unf_gen[i]->Draw("same p");
    
    drawText(Form("%2.0f-%.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  cDataMC_raw->cd(1);
  putCMSPrel();
  drawText(Form("ak Pu R = 0.%d PF jets, %2.0f < |eta| < %2.0f",radius, etaLow, etaHigh), 0.2,0.13,15);

  TLegend * leg = myLegend(0.2,0.2,0.4,0.4);
  leg->AddEntry(hRatio[0],"Data Raw / MC Reco","pl");
  leg->AddEntry(hRatio_unf[0],"Data Unfo / MC Reco","pl");
  leg->AddEntry(hRatio_gen[0],"Data Raw / MC Gen","pl");
  leg->AddEntry(hRatio_unf_gen[0],"Data Unfo / MC Gen","pl");
  leg->SetTextSize(0.04);
  leg->Draw();
  
  cDataMC_raw->SaveAs(Form("May22/PbPb_DataMC_ratio_%s_R%d_%d.pdf",etaWidth,radius,date.GetDate()),"RECREATE");


  
}
