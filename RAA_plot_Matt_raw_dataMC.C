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
#include "Headers/plot.h"


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


void RAA_plot_Matt_raw_dataMC(char* etaWidth = (char*)"10_eta_18", Int_t radius = 4){

  
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
  
  // Make Matt's requested plots about the unfolding prior distributions.
  // need to make ratio histograms for RAA, and PbPb spectra

  TFile * fNoSmear = TFile::Open(Form("../../Output/Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_noFakeWeight_unfold_mcclosure_oppside_trgMC_%s_40GeVCut_akPF_20150506.root",radius,etaWidth));
  TFile * fGenSmear = TFile::Open(Form("../../Output/Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_noFakeWeight_unfold_mcclosure_oppside_trgMC_GenSmear_%s_40GeVCut_akPF_20150507.root",radius,etaWidth));
  TFile * fRecoSmear = TFile::Open(Form("../../Output/Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_noFakeWeight_unfold_mcclosure_oppside_trgMC_RecoSmear_%s_40GeVCut_akPF_20150507.root",radius,etaWidth));
  TFile * fBothSmear = TFile::Open(Form("../../Output/Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_noFakeWeight_unfold_mcclosure_oppside_trgMC_BothSmear_%s_40GeVCut_akPF_20150507.root",radius,etaWidth));
  
  TH1F * RAA_Bayesian[nbins_cent], * RAA_Bayesian_GenSmear[nbins_cent], * RAA_Bayesian_RecoSmear[nbins_cent], * RAA_Bayesian_BothSmear[nbins_cent];
  TH1F * hPbPb_Bayesian[nbins_cent], * hPbPb_Bayesian_GenSmear[nbins_cent], * hPbPb_Bayesian_RecoSmear[nbins_cent], * hPbPb_Bayesian_BothSmear[nbins_cent];

  TH1F * RAA_Sys_GenSmear[nbins_cent], * RAA_Sys_RecoSmear[nbins_cent],* RAA_Sys_BothSmear[nbins_cent];
  TH1F * hPbPb_Sys_GenSmear[nbins_cent], * hPbPb_Sys_RecoSmear[nbins_cent],* hPbPb_Sys_BothSmear[nbins_cent];
  
  // get the input files:
  TFile * fin = TFile::Open(Form("../../Output/Pawan_ntuple_PbPb_data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius));

  // get the histograms
  TH1F * hData[nbins_cent], * hMC[nbins_cent];
  TH1F * hRatio[nbins_cent];

  for(int i = 0; i<nbins_cent; ++i){

    RAA_Bayesian[i] = (TH1F*)fNoSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian[i]->Print("base");
    RAA_Bayesian_GenSmear[i] = (TH1F*)fGenSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian_GenSmear[i]->Print("base");
    RAA_Bayesian_RecoSmear[i] = (TH1F*)fRecoSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian_RecoSmear[i]->Print("base");
    RAA_Bayesian_BothSmear[i] = (TH1F*)fBothSmear->Get(Form("RAA_bayesian_cent%d",i));
    RAA_Bayesian_BothSmear[i]->Print("base");

    RAA_Sys_GenSmear[i] = (TH1F*)RAA_Bayesian_GenSmear[i]->Clone(Form("RAA_Sys_GenSmear_cent%d",i));
    RAA_Sys_GenSmear[i]->Divide(RAA_Bayesian[i]);
    RAA_Sys_RecoSmear[i] = (TH1F*)RAA_Bayesian_RecoSmear[i]->Clone(Form("RAA_Sys_RecoSmear_cent%d",i));
    RAA_Sys_RecoSmear[i]->Divide(RAA_Bayesian[i]);
    RAA_Sys_BothSmear[i] = (TH1F*)RAA_Bayesian_BothSmear[i]->Clone(Form("RAA_Sys_BothSmear_cent%d",i));
    RAA_Sys_BothSmear[i]->Divide(RAA_Bayesian[i]);
    
    hPbPb_Bayesian[i] = (TH1F*)fNoSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian[i]->Print("base");
    hPbPb_Bayesian_GenSmear[i] = (TH1F*)fGenSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian_GenSmear[i]->Print("base");
    hPbPb_Bayesian_RecoSmear[i] = (TH1F*)fRecoSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian_RecoSmear[i]->Print("base");
    hPbPb_Bayesian_BothSmear[i] = (TH1F*)fBothSmear->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    hPbPb_Bayesian_BothSmear[i]->Print("base");

    hPbPb_Sys_GenSmear[i] = (TH1F*)hPbPb_Bayesian_GenSmear[i]->Clone(Form("hPbPb_Sys_GenSmear_cent%d",i));
    hPbPb_Sys_GenSmear[i]->Divide(hPbPb_Bayesian[i]);
    hPbPb_Sys_RecoSmear[i] = (TH1F*)hPbPb_Bayesian_RecoSmear[i]->Clone(Form("hPbPb_Sys_RecoSmear_cent%d",i));
    hPbPb_Sys_RecoSmear[i]->Divide(hPbPb_Bayesian[i]);
    hPbPb_Sys_BothSmear[i] = (TH1F*)hPbPb_Bayesian_BothSmear[i]->Clone(Form("hPbPb_Sys_BothSmear_cent%d",i));
    hPbPb_Sys_BothSmear[i]->Divide(hPbPb_Bayesian[i]);
    
    hData[i] = (TH1F*)fin->Get(Form("hpbpb_raw_HLTComb_R%d_%s_cent%d",radius, etaWidth, i));
    hData[i]->Print("base");
    hMC[i] = (TH1F*)fin->Get(Form("hpbpb_JetComb_raw_R%d_%s_cent%d",radius, etaWidth, i));
    hMC[i]->Print("base");
    hRatio[i] = (TH1F*)hData[i]->Clone(Form("Data_MC_rawpT_ratio_cent%d",i));
    hRatio[i]->Divide(hMC[i]);
    hRatio[i]->Print("base");
    
  }

  TLine * line = new TLine(40,1,299,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  
  // make the bayesian unfolding prior systematics plot
  TCanvas * cRAA_Sys = new TCanvas("cRAA_Sys","",1000,800);
  makeMultiPanelCanvasWithGap(cRAA_Sys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend * lRAA_Sys = myLegend(0.2,0.6,0.4,0.8);
  for(int i = 0; i<nbins_cent; ++i){

    cRAA_Sys->cd(nbins_cent-i);
  
    RAA_Sys_GenSmear[i]->SetMarkerStyle(33);
    RAA_Sys_GenSmear[i]->SetMarkerColor(kRed);
    RAA_Sys_GenSmear[i]->SetAxisRange(40,299,"X");
    RAA_Sys_GenSmear[i]->SetAxisRange(0,2,"Y");
    makeHistTitle(RAA_Sys_GenSmear[i]," ","Jet p_{T} (GeV/c)","RAA ratio smear/ no smear");
    RAA_Sys_GenSmear[i]->Draw();

    RAA_Sys_RecoSmear[i]->SetMarkerStyle(33);
    RAA_Sys_RecoSmear[i]->SetMarkerColor(kBlue);
    RAA_Sys_RecoSmear[i]->Draw("same");

    RAA_Sys_BothSmear[i]->SetMarkerStyle(33);
    RAA_Sys_BothSmear[i]->SetMarkerColor(kGreen);
    RAA_Sys_BothSmear[i]->Draw("same");

    drawText(Form("%2.0f-%.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    line->Draw();
    
  }

  cRAA_Sys->cd(1);
  putCMSPrel();
  
  lRAA_Sys->AddEntry(RAA_Sys_GenSmear[0],"Gen pT smeared","pl");
  lRAA_Sys->AddEntry(RAA_Sys_RecoSmear[0],"Reco pT smeared","pl");
  lRAA_Sys->AddEntry(RAA_Sys_BothSmear[0],"Gen and Reco pT smeared","pl");
  lRAA_Sys->SetTextSize(0.04);
  lRAA_Sys->Draw();
  drawText(Form("ak Pu R = 0.%d PF jets, %2.0f < |eta| < %2.0f",radius, etaLow, etaHigh), 0.2,0.2,15);

  cRAA_Sys->SaveAs(Form("RAA_bayesian_unfolded_prior_smear_shift_systematics_%s_R%d_%d.pdf",etaWidth,radius,date.GetDate()),"RECREATE");

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

    hPbPb_Sys_BothSmear[i]->SetMarkerStyle(33);
    hPbPb_Sys_BothSmear[i]->SetMarkerColor(kGreen);
    hPbPb_Sys_BothSmear[i]->Draw("same");

    drawText(Form("%2.0f-%.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
    line->Draw();
    
  }

  chPbPb_Sys->cd(1);
  putCMSPrel();
  
  lhPbPb_Sys->AddEntry(hPbPb_Sys_GenSmear[0],"Gen pT smeared","pl");
  lhPbPb_Sys->AddEntry(hPbPb_Sys_RecoSmear[0],"Reco pT smeared","pl");
  lhPbPb_Sys->AddEntry(hPbPb_Sys_BothSmear[0],"Gen and Reco pT smeared","pl");
  lhPbPb_Sys->SetTextSize(0.04);
  lhPbPb_Sys->Draw();
  drawText(Form("ak Pu R = 0.%d PF jets, %2.0f < |eta| < %2.0f",radius, etaLow, etaHigh), 0.2,0.2,15);

  chPbPb_Sys->SaveAs(Form("PbPb_spectra_bayesian_unfolded_prior_smear_shift_systematics_%s_R%d_%d.pdf",etaWidth,radius,date.GetDate()),"RECREATE");


  // make the plots - 3x2 plots.
  TCanvas * cDataMC_raw = new TCanvas("cDataMC_raw","",1000,800);
  makeMultiPanelCanvasWithGap(cDataMC_raw,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0; i<nbins_cent; ++i){

    cDataMC_raw->cd(nbins_cent-i);
    cDataMC_raw->cd(nbins_cent-i)->SetLogy();

    hRatio[i]->SetMarkerStyle(20);
    hRatio[i]->SetMarkerColor(kBlack);
    makeHistTitle(hRatio[i]," ","Raw Jet p_{T} (GeV/c)","Data/MC arbitraty scale");
    hRatio[i] = (TH1F*)hRatio[i]->Rebin(nbins_pt, Form("Data_MC_rawpT_ratio_cent%d",i), boundaries_pt);
    divideBinWidth(hRatio[i]);
    hRatio[i]->SetAxisRange(30,299,"X");
    hRatio[i]->GetYaxis()->SetNdivisions(2,kFALSE);
    hRatio[i]->Draw();
    drawText(Form("%2.0f-%.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  cDataMC_raw->cd(1);
  putCMSPrel();
  drawText(Form("ak Pu R = 0.%d PF jets, %2.0f < |eta| < %2.0f",radius, etaLow, etaHigh), 0.2,0.2,15);

  cDataMC_raw->SaveAs(Form("../../Plots/PbPb_DataMC_rawpT_%s_R%d_%d.pdf",etaWidth,radius,date.GetDate()),"RECREATE");


  



  
}
