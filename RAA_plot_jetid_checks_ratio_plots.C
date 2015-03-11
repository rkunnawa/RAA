// Raghav Kunnawalkam Elayavalli
// Feb 21 2014
// Rugers 

//
// Macro to study the effect of jet id cuts. this makes ratio plots from data and MC so we can study the effect. 
//
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
 
static const int nbins_pt = 39;
static const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 272, 300, 
  330, 362, 395, 430,
  468, 507, 548, 592,
  638, 686, 1000 
};


using namespace std;

void RAA_plot_jetid_checks_ratio_plots(int radius = 3, char *algo = "Pu", char *jet_type = "PF"){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  TFile *fData = TFile::Open("/Users/keraghav/WORK/RAA/Output/PbPb_data_akPuPF_20150308.root");
  TFile *fMC = TFile::Open("/Users/keraghav/WORK/RAA/Output/PbPb_mc_akPuPF_20150308.root");

  //get the 2d histograms 
  static const int ptSelection = 19;
  static const int ptBoundary[ptSelection+1] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};


#if 0
  TH2F* hData_chMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hData_eMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hData_phMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hData_neMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hData_muMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hData_eMaxSumcand_jtpt_ptselection[ptSelection];
  //TH2F* hData_eMaxSumcand_jtpt_refptselection[ptSelection];
  //TH2F* hData_eMaxSumcand_refpt_refptselection[ptSelection];
  TH2F* hData_eMaxSumcand_chMaxJtpt_ptselection[ptSelection];
  //TH2F* hData_eMaxSumcand_chMaxJtpt_refptselection[ptSelection];
  TH2F* hData_eMaxJtpt_chMaxJtpt_ptselection[ptSelection];
  TH2F* hData_neMaxJtpt_chMaxJtpt_ptselection[ptSelection];
  TH2F* hData_phMaxJtpt_chMaxJtpt_ptselection[ptSelection];
  TH2F* hData_muMaxJtpt_chMaxJtpt_ptselection[ptSelection];

  
  TH2F* hMC_chMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hMC_eMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hMC_phMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hMC_neMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hMC_muMaxJtpt_jtpt_ptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_jtpt_ptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_jtpt_refptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_refpt_refptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_chMaxJtpt_ptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_chMaxJtpt_refptselection[ptSelection];
  TH2F* hMC_eMaxJtpt_chMaxJtpt_ptselection[ptSelection];
  TH2F* hMC_neMaxJtpt_chMaxJtpt_ptselection[ptSelection];
  TH2F* hMC_phMaxJtpt_chMaxJtpt_ptselection[ptSelection];
  TH2F* hMC_muMaxJtpt_chMaxJtpt_ptselection[ptSelection];
  
  //TCanvas *cchMaxJtpt_jtpt[ptSelection], *cphMaxJtpt_jtpt[ptSelection], *cneMaxJtpt_jtpt[ptSelection], *cmuMaxJtpt_jtpt[ptSelection], *ceMaxJtpt_jtpt[ptSelection], *ceMaxSumcand_jtpt[ptSelection];
  TCanvas *ceMaxSumcand_chMaxJtpt[ptSelection], *ceMaxJtpt_chMaxJtpt[ptSelection], *cphMaxJtpt_chMaxJtpt[ptSelection], *cneMaxJtpt_chMaxJtpt[ptSelection], *cmuMaxJtpt_chMaxJtpt[ptSelection];
  TH2F *hData_chMaxJtpt_jtpt, *hData_phMaxJtpt_jtpt, *hData_neMaxJtpt_jtpt, *hData_muMaxJtpt_jtpt, *hData_eMaxJtpt_jtpt, *hData_eMaxSumcand_jtpt;
  TH2F *hMC_chMaxJtpt_jtpt, *hMC_phMaxJtpt_jtpt, *hMC_neMaxJtpt_jtpt, *hMC_muMaxJtpt_jtpt, *hMC_eMaxJtpt_jtpt, *hMC_eMaxSumcand_jtpt_jtpt, *hMC_eMaxSumcand_jtpt_refpt, *hMC_eMaxSumcand_refpt_refpt;
  
  for(int a = 0;a<ptSelection;a++){

    hData_chMaxJtpt_jtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_chMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_phMaxJtpt_jtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_phMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_neMaxJtpt_jtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_neMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_muMaxJtpt_jtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_muMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_eMaxJtpt_jtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_eMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_eMaxSumcand_jtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_eMaxSumcand_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_eMaxSumcand_chMaxJtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_eMaxSumcand_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_eMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_eMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_neMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_neMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_phMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_phMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hData_muMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fData->Get(Form("hpbpb_muMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));

    hMC_chMaxJtpt_jtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_chMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_phMaxJtpt_jtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_phMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_neMaxJtpt_jtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_neMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_muMaxJtpt_jtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_muMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_eMaxJtpt_jtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_eMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_eMaxSumcand_jtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_eMaxSumcand_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_eMaxSumcand_jtpt_refptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_eMaxSumcand_jtpt_%d_refpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_eMaxSumcand_refpt_refptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_eMaxSumcand_refpt_%d_refpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_eMaxSumcand_chMaxJtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_eMaxSumcand_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_eMaxSumcand_chMaxJtpt_refptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_eMaxSumcand_chMaxJtpt_%d_refpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_eMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_eMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_neMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_neMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_phMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_phMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
    hMC_muMaxJtpt_chMaxJtpt_ptselection[a] = (TH2F*)fMC->Get(Form("hpbpb_muMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));

    if(a == 0) {
      hData_chMaxJtpt_jtpt = (TH2F*)hData_chMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_chMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hData_phMaxJtpt_jtpt = (TH2F*)hData_phMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_phMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hData_neMaxJtpt_jtpt = (TH2F*)hData_neMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_neMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hData_muMaxJtpt_jtpt = (TH2F*)hData_muMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_muMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hData_eMaxJtpt_jtpt = (TH2F*)hData_eMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_eMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hData_eMaxSumcand_jtpt = (TH2F*)hData_eMaxSumcand_jtpt_ptselection[a]->Clone(Form("hpbpb_eMaxSumcand_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_chMaxJtpt_jtpt = (TH2F*)hMC_chMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_chMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_phMaxJtpt_jtpt = (TH2F*)hMC_phMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_phMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_neMaxJtpt_jtpt = (TH2F*)hMC_neMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_neMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_muMaxJtpt_jtpt = (TH2F*)hMC_muMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_muMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_eMaxJtpt_jtpt = (TH2F*)hMC_eMaxJtpt_jtpt_ptselection[a]->Clone(Form("hpbpb_eMaxJtpt_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_eMaxSumcand_jtpt_jtpt = (TH2F*)hMC_eMaxSumcand_jtpt_ptselection[a]->Clone(Form("hpbpb_eMaxSumcand_jtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_eMaxSumcand_jtpt_refpt = (TH2F*)hMC_eMaxSumcand_jtpt_refptselection[a]->Clone(Form("hpbpb_eMaxSumcand_jtpt_%d_refpt_%d",ptBoundary[a],ptBoundary[a+1]));
      hMC_eMaxSumcand_refpt_refpt = (TH2F*)hMC_eMaxSumcand_refpt_refptselection[a]->Clone(Form("hpbpb_eMaxSumcand_refpt_%d_refpt_%d",ptBoundary[a],ptBoundary[a+1]));

    } else {
      hData_chMaxJtpt_jtpt->Add(hData_chMaxJtpt_jtpt_ptselection[a]);
      hData_phMaxJtpt_jtpt->Add(hData_phMaxJtpt_jtpt_ptselection[a]);
      hData_neMaxJtpt_jtpt->Add(hData_neMaxJtpt_jtpt_ptselection[a]);
      hData_muMaxJtpt_jtpt->Add(hData_muMaxJtpt_jtpt_ptselection[a]);
      hData_eMaxJtpt_jtpt->Add(hData_eMaxJtpt_jtpt_ptselection[a]);
      hData_eMaxSumcand_jtpt->Add(hData_eMaxSumcand_jtpt_ptselection[a]);
      hMC_chMaxJtpt_jtpt->Add(hMC_chMaxJtpt_jtpt_ptselection[a]);
      hMC_phMaxJtpt_jtpt->Add(hMC_phMaxJtpt_jtpt_ptselection[a]);
      hMC_neMaxJtpt_jtpt->Add(hMC_neMaxJtpt_jtpt_ptselection[a]);
      hMC_muMaxJtpt_jtpt->Add(hMC_muMaxJtpt_jtpt_ptselection[a]);
      hMC_eMaxJtpt_jtpt->Add(hMC_eMaxJtpt_jtpt_ptselection[a]);
      hMC_eMaxSumcand_jtpt_jtpt->Add(hMC_eMaxSumcand_jtpt_ptselection[a]);
      hMC_eMaxSumcand_jtpt_refpt->Add(hMC_eMaxSumcand_jtpt_refptselection[a]);
      hMC_eMaxSumcand_refpt_refpt->Add(hMC_eMaxSumcand_refpt_refptselection[a]);
    }

  }


  TCanvas * ceMaxSumcand_jtpt = new TCanvas("ceMaxSumcand_jtpt","",1200,800);
  makeMultiPanelCanvas(ceMaxSumcand_jtpt,4,1,0.0,0.0,0.2,0.15,0.07);
  ceMaxSumcand_jtpt->cd(1)->SetLogz();
  hData_eMaxSumcand_jtpt->SetTitle(" ");
  hData_eMaxSumcand_jtpt->SetXTitle("Reco jet p_{T}");
  hData_eMaxSumcand_jtpt->SetYTitle("eMax/Sumcand");
  hData_eMaxSumcand_jtpt->Draw("colz");
  drawText("Data 0-30% 10 < jtpt < 200",0.3,0.8,14);
  ceMaxSumcand_jtpt->cd(2)->SetLogz();
  hMC_eMaxSumcand_jtpt_jtpt->SetTitle(" ");
  hMC_eMaxSumcand_jtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hMC_eMaxSumcand_jtpt_jtpt->SetYTitle("eMax/Sumcand");
  hMC_eMaxSumcand_jtpt_jtpt->Draw("colz");
  drawText("MC 0-30%, 10 < jtpt < 200",0.3,0.8,14);
  ceMaxSumcand_jtpt->cd(3)->SetLogz();
  hMC_eMaxSumcand_jtpt_refpt->SetTitle(" ");
  hMC_eMaxSumcand_jtpt_refpt->SetXTitle("Reco jet p_{T}");
  hMC_eMaxSumcand_jtpt_refpt->SetYTitle("eMax/Sumcand");
  hMC_eMaxSumcand_jtpt_refpt->Draw("colz");
  drawText("MC 0-30%, 10 < refpt < 200",0.3,0.8,14);
  ceMaxSumcand_jtpt->cd(4)->SetLogz();
  hMC_eMaxSumcand_refpt_refpt->SetTitle(" ");
  hMC_eMaxSumcand_refpt_refpt->SetXTitle("Gen jet p_{T}");
  hMC_eMaxSumcand_refpt_refpt->SetYTitle("eMax/Sumcand");
  hMC_eMaxSumcand_refpt_refpt->Draw("colz");
  drawText("MC 0-30%, 10 < refpt < 200",0.3,0.8,14);
  ceMaxSumcand_jtpt->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_eMaxSumcand_jtpt_10_pt_200.pdf","RECREATE");

  TCanvas * ceMaxJtpt_jtpt = new TCanvas("ceMaxJtpt_jtpt","",1200,800);
  makeMultiPanelCanvas(ceMaxJtpt_jtpt,2,1,0.0,0.0,0.2,0.15,0.07);
  ceMaxJtpt_jtpt->cd(1)->SetLogz();
  hData_eMaxJtpt_jtpt->SetTitle(" ");
  hData_eMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hData_eMaxJtpt_jtpt->SetYTitle("eMax/Jtpt");
  hData_eMaxJtpt_jtpt->Draw("colz");
  drawText("Data 0-30% 10 < jtpt < 200",0.3,0.8,14);
  ceMaxJtpt_jtpt->cd(2)->SetLogz();
  hMC_eMaxJtpt_jtpt->SetTitle(" ");
  hMC_eMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hMC_eMaxJtpt_jtpt->SetYTitle("eMax/Sumcand");
  hMC_eMaxJtpt_jtpt->Draw("colz");
  drawText("MC 0-30%, 10 < jtpt < 200",0.3,0.8,14);
  ceMaxJtpt_jtpt->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_eMaxJtpt_jtpt_10_pt_200.pdf","RECREATE");

  TCanvas * cchMaxJtpt_jtpt = new TCanvas("cchMaxJtpt_jtpt","",1200,800);
  makeMultiPanelCanvas(cchMaxJtpt_jtpt,2,1,0.0,0.0,0.2,0.15,0.07);
  cchMaxJtpt_jtpt->cd(1)->SetLogz();
  hData_chMaxJtpt_jtpt->SetTitle(" ");
  hData_chMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hData_chMaxJtpt_jtpt->SetYTitle("chMax/Jtpt");
  hData_chMaxJtpt_jtpt->Draw("colz");
  drawText("Data 0-30% 10 < jtpt < 200",0.3,0.8,14);
  cchMaxJtpt_jtpt->cd(2)->SetLogz();
  hMC_chMaxJtpt_jtpt->SetTitle(" ");
  hMC_chMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hMC_chMaxJtpt_jtpt->SetYTitle("eMax/Sumcand");
  hMC_chMaxJtpt_jtpt->Draw("colz");
  drawText("MC 0-30%, 10 < jtpt < 200",0.3,0.8,14);
  cchMaxJtpt_jtpt->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_chMaxJtpt_jtpt_10_pt_200.pdf","RECREATE");

  TCanvas * cneMaxJtpt_jtpt = new TCanvas("cneMaxJtpt_jtpt","",1200,800);
  makeMultiPanelCanvas(cneMaxJtpt_jtpt,2,1,0.0,0.0,0.2,0.15,0.07);
  cneMaxJtpt_jtpt->cd(1)->SetLogz();
  hData_neMaxJtpt_jtpt->SetTitle(" ");
  hData_neMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hData_neMaxJtpt_jtpt->SetYTitle("neMax/Jtpt");
  hData_neMaxJtpt_jtpt->Draw("colz");
  drawText("Data 0-30% 10 < jtpt < 200",0.3,0.8,14);
  cneMaxJtpt_jtpt->cd(2)->SetLogz();
  hMC_neMaxJtpt_jtpt->SetTitle(" ");
  hMC_neMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hMC_neMaxJtpt_jtpt->SetYTitle("eMax/Sumcand");
  hMC_neMaxJtpt_jtpt->Draw("colz");
  drawText("MC 0-30%, 10 < jtpt < 200",0.3,0.8,14);
  cneMaxJtpt_jtpt->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_neMaxJtpt_jtpt_10_pt_200.pdf","RECREATE");

  TCanvas * cphMaxJtpt_jtpt = new TCanvas("cphMaxJtpt_jtpt","",1200,800);
  makeMultiPanelCanvas(cphMaxJtpt_jtpt,2,1,0.0,0.0,0.2,0.15,0.07);
  cphMaxJtpt_jtpt->cd(1)->SetLogz();
  hData_phMaxJtpt_jtpt->SetTitle(" ");
  hData_phMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hData_phMaxJtpt_jtpt->SetYTitle("phMax/Jtpt");
  hData_phMaxJtpt_jtpt->Draw("colz");
  drawText("Data 0-30% 10 < jtpt < 200",0.3,0.8,14);
  cphMaxJtpt_jtpt->cd(2)->SetLogz();
  hMC_phMaxJtpt_jtpt->SetTitle(" ");
  hMC_phMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hMC_phMaxJtpt_jtpt->SetYTitle("eMax/Sumcand");
  hMC_phMaxJtpt_jtpt->Draw("colz");
  drawText("MC 0-30%, 10 < jtpt < 200",0.3,0.8,14);
  cphMaxJtpt_jtpt->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_phMaxJtpt_jtpt_10_pt_200.pdf","RECREATE");
  
  TCanvas * cmuMaxJtpt_jtpt = new TCanvas("cmuMaxJtpt_jtpt","",1200,800);
  makeMultiPanelCanvas(cmuMaxJtpt_jtpt,2,1,0.0,0.0,0.2,0.15,0.07);
  cmuMaxJtpt_jtpt->cd(1)->SetLogz();
  hData_muMaxJtpt_jtpt->SetTitle(" ");
  hData_muMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hData_muMaxJtpt_jtpt->SetYTitle("muMax/Jtpt");
  hData_muMaxJtpt_jtpt->Draw("colz");
  drawText("Data 0-30% 10 < jtpt < 200",0.3,0.8,14);
  cmuMaxJtpt_jtpt->cd(2)->SetLogz();
  hMC_muMaxJtpt_jtpt->SetTitle(" ");
  hMC_muMaxJtpt_jtpt->SetXTitle("Reco jet p_{T}");
  hMC_muMaxJtpt_jtpt->SetYTitle("eMax/Sumcand");
  hMC_muMaxJtpt_jtpt->Draw("colz");
  drawText("MC 0-30%, 10 < jtpt < 200",0.3,0.8,14);
  cmuMaxJtpt_jtpt->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_muMaxJtpt_jtpt_10_pt_200.pdf","RECREATE");

  
  for(int a = 3;a<ptSelection;a++){
    ceMaxSumcand_chMaxJtpt[a] = new TCanvas(Form("ceMaxSumcand_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]),"",1200,800);
    makeMultiPanelCanvas(ceMaxSumcand_chMaxJtpt[a],3,1,0.0,0.0,0.2,0.15,0.07);
    ceMaxSumcand_chMaxJtpt[a]->cd(1)->SetLogz();
    hData_eMaxSumcand_chMaxJtpt_ptselection[a]->SetTitle(" ");
    hData_eMaxSumcand_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hData_eMaxSumcand_chMaxJtpt_ptselection[a]->SetYTitle("eMax/Sumcand");
    hData_eMaxSumcand_chMaxJtpt_ptselection[a]->RebinX(5);
    hData_eMaxSumcand_chMaxJtpt_ptselection[a]->RebinY(10);
    hData_eMaxSumcand_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hData_eMaxSumcand_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText("Data",0.45,0.7,15);
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.3,0.8,15);
    ceMaxSumcand_chMaxJtpt[a]->cd(2)->SetLogz();
    hMC_eMaxSumcand_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hMC_eMaxSumcand_chMaxJtpt_ptselection[a]->SetYTitle("eMax/Sumcand");
    hMC_eMaxSumcand_chMaxJtpt_ptselection[a]->RebinX(5);
    hMC_eMaxSumcand_chMaxJtpt_ptselection[a]->RebinY(10);
    hMC_eMaxSumcand_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hMC_eMaxSumcand_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.2,0.8,15);
    drawText("MC",0.45,0.7,15);
    ceMaxSumcand_chMaxJtpt[a]->cd(3)->SetLogz();
    hMC_eMaxSumcand_chMaxJtpt_refptselection[a]->SetXTitle("chMax/jtpt");
    hMC_eMaxSumcand_chMaxJtpt_refptselection[a]->SetYTitle("eMax/Sumcand");
    hMC_eMaxSumcand_chMaxJtpt_refptselection[a]->RebinX(5);
    hMC_eMaxSumcand_chMaxJtpt_refptselection[a]->RebinY(10);
    hMC_eMaxSumcand_chMaxJtpt_refptselection[a]->SetAxisRange(0,1.5,"X");
    hMC_eMaxSumcand_chMaxJtpt_refptselection[a]->Draw("colz");
    drawText(Form("%d < ref pt < %d",ptBoundary[a],ptBoundary[a+1]),0.2,0.8,15);
    drawText("MC",0.45,0.7,15);
    ceMaxSumcand_chMaxJtpt[a]->SaveAs(Form("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_eMaxSumcand_chMaxJtpt_%d_jtpt_%d.pdf",ptBoundary[a],ptBoundary[a+1]),"RECREATE");

    ceMaxJtpt_chMaxJtpt[a] = new TCanvas(Form("ceMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]),"",1200,800);
    makeMultiPanelCanvas(ceMaxJtpt_chMaxJtpt[a],2,1,0.0,0.0,0.2,0.15,0.07);
    ceMaxJtpt_chMaxJtpt[a]->cd(1)->SetLogz();
    hData_eMaxJtpt_chMaxJtpt_ptselection[a]->SetTitle(" ");
    hData_eMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hData_eMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("eMax/jtpt");
    hData_eMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(1);
    hData_eMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hData_eMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,3,"Y");
    hData_eMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText("Data",0.45,0.7,15);
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.3,0.8,15);
    ceMaxJtpt_chMaxJtpt[a]->cd(2)->SetLogz();
    hMC_eMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hMC_eMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("eMax/jtpt");
    hMC_eMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(1);
    hMC_eMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hMC_eMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,3,"Y");
    hMC_eMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.2,0.8,15);
    drawText("MC",0.45,0.7,15);
    // ceMaxJtpt_chMaxJtpt[a]->cd(3)->SetLogz();
    // hMC_eMaxJtpt_chMaxJtpt_refptselection[a]->SetXTitle("chMax/jtpt");
    // hMC_eMaxJtpt_chMaxJtpt_refptselection[a]->SetYTitle("eMax/Sumcand");
    // hMC_eMaxJtpt_chMaxJtpt_refptselection[a]->Draw("colz");
    // drawText(Form("%d < ref pt < %d",ptSelection[a],ptSelection[a+1]),0.2,0.8,15);
    // drawText("MC",0.45,0.7,15);
    ceMaxJtpt_chMaxJtpt[a]->SaveAs(Form("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_eMaxJtpt_chMaxJtpt_%d_jtpt_%d.pdf",ptBoundary[a],ptBoundary[a+1]),"RECREATE");

    cneMaxJtpt_chMaxJtpt[a] = new TCanvas(Form("cneMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]),"",1200,800);
    makeMultiPanelCanvas(cneMaxJtpt_chMaxJtpt[a],2,1,0.0,0.0,0.2,0.15,0.07);
    cneMaxJtpt_chMaxJtpt[a]->cd(1)->SetLogz();
    hData_neMaxJtpt_chMaxJtpt_ptselection[a]->SetTitle(" ");
    hData_neMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hData_neMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("neMax/jtpt");
    hData_neMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(5);
    hData_neMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hData_neMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"Y");
    hData_neMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText("Data",0.45,0.7,15);
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.3,0.8,15);
    cneMaxJtpt_chMaxJtpt[a]->cd(2)->SetLogz();
    hMC_neMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hMC_neMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("neMax/jtpt");
    hMC_neMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(5);
    hMC_neMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hMC_neMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"Y");
    hMC_neMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.2,0.8,15);
    drawText("MC",0.45,0.7,15);
    cneMaxJtpt_chMaxJtpt[a]->SaveAs(Form("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_neMaxJtpt_chMaxJtpt_%d_jtpt_%d.pdf",ptBoundary[a],ptBoundary[a+1]),"RECREATE");


    cphMaxJtpt_chMaxJtpt[a] = new TCanvas(Form("cphMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]),"",1200,800);
    makeMultiPanelCanvas(cphMaxJtpt_chMaxJtpt[a],2,1,0.0,0.0,0.2,0.15,0.07);
    cphMaxJtpt_chMaxJtpt[a]->cd(1)->SetLogz();
    hData_phMaxJtpt_chMaxJtpt_ptselection[a]->SetTitle(" ");
    hData_phMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hData_phMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("phMax/jtpt");
    hData_phMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(5);
    hData_phMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hData_phMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,2,"Y");
    hData_phMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText("Data",0.45,0.7,15);
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.3,0.8,15);
    cphMaxJtpt_chMaxJtpt[a]->cd(2)->SetLogz();
    hMC_phMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hMC_phMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("phMax/jtpt");
    hMC_phMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(5);
    hMC_phMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hMC_phMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,2,"Y");
    hMC_phMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.2,0.8,15);
    drawText("MC",0.45,0.7,15);
    cphMaxJtpt_chMaxJtpt[a]->SaveAs(Form("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_phMaxJtpt_chMaxJtpt_%d_jtpt_%d.pdf",ptBoundary[a],ptBoundary[a+1]),"RECREATE");

    //cmuMaxJtpt_chMaxJtpt[a] = new TCanvas(Form("cmuMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]),"",1200,800);
    cmuMaxJtpt_chMaxJtpt[a] = new TCanvas(Form("cmuMaxJtpt_chMaxJtpt_%d_jtpt_%d",ptBoundary[a],ptBoundary[a+1]),"",1200,800);
    makeMultiPanelCanvas(cmuMaxJtpt_chMaxJtpt[a],2,1,0.0,0.0,0.2,0.15,0.07);
    cmuMaxJtpt_chMaxJtpt[a]->cd(1)->SetLogz();
    hData_muMaxJtpt_chMaxJtpt_ptselection[a]->SetTitle(" ");
    hData_muMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hData_muMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("muMax/jtpt");
    hData_muMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(5);
    hData_muMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hData_muMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"Y");
    hData_muMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText("Data",0.45,0.7,15);
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.3,0.8,15);
    cmuMaxJtpt_chMaxJtpt[a]->cd(2)->SetLogz();
    hMC_muMaxJtpt_chMaxJtpt_ptselection[a]->SetXTitle("chMax/jtpt");
    hMC_muMaxJtpt_chMaxJtpt_ptselection[a]->SetYTitle("muMax/Sumcand");
    hMC_muMaxJtpt_chMaxJtpt_ptselection[a]->Rebin2D(5);
    hMC_muMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"X");
    hMC_muMaxJtpt_chMaxJtpt_ptselection[a]->SetAxisRange(0,1.5,"Y");
    hMC_muMaxJtpt_chMaxJtpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < reco pt < %d",ptBoundary[a],ptBoundary[a+1]),0.2,0.8,15);
    drawText("MC",0.45,0.7,15);
    cmuMaxJtpt_chMaxJtpt[a]->SaveAs(Form("/Users/keraghav/WORK/RAA/Plots/PbPb_jetid_muMaxJtpt_chMaxJtpt_%d_jtpt_%d.pdf",ptBoundary[a],ptBoundary[a+1]),"RECREATE");
    
  }

  //TCanvas *cchMaxJtpt_jtpt


  
  //get the histograms:
  TH1F* hData_Jet80 = (TH1F*)fData->Get("hpbpb_Jet80");
  TH1F* hData_Jet80_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p01");
  TH1F* hData_Jet80_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p02");
  TH1F* hData_Jet80_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p03");
  TH1F* hData_Jet80_chMaxJtpt0p04 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p04");
  TH1F* hData_Jet80_chMaxJtpt0p05 = (TH1F*)fData->Get("hpbpb_Jet80_chMaxJtpt0p05");
  TH1F* hData_Jet80_eMaxJtpt0p1 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p1");
  TH1F* hData_Jet80_eMaxJtpt0p2 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p2");
  TH1F* hData_Jet80_eMaxJtpt0p3 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p3");
  TH1F* hData_Jet80_eMaxJtpt0p4 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p4");
  TH1F* hData_Jet80_eMaxJtpt0p5 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5");
  TH1F* hData_Jet80_eMaxJtpt0p6 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6");
  TH1F* hData_Jet80_eMaxJtpt0p7 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7");
  TH1F* hData_Jet80_eMaxJtpt0p8 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p8");
  TH1F* hData_Jet80_eMaxJtpt0p9 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p9");
  TH1F* hData_Jet80_eMaxSumcand0p1 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p1");
  TH1F* hData_Jet80_eMaxSumcand0p2 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p2");
  TH1F* hData_Jet80_eMaxSumcand0p3 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p3");
  TH1F* hData_Jet80_eMaxSumcand0p4 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p4");
  TH1F* hData_Jet80_eMaxSumcand0p5 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p5");
  TH1F* hData_Jet80_eMaxSumcand0p6 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p6");
  TH1F* hData_Jet80_eMaxSumcand0p7 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p7");
  TH1F* hData_Jet80_eMaxSumcand0p8 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p8");
  TH1F* hData_Jet80_eMaxSumcand0p9 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxSumcand0p9");
  TH1F* hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hData_Jet65 = (TH1F*)fData->Get("hpbpb_Jet65");
  TH1F* hData_Jet65_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p01");
  TH1F* hData_Jet65_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p02");
  TH1F* hData_Jet65_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p03");
  TH1F* hData_Jet65_chMaxJtpt0p04 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p04");
  TH1F* hData_Jet65_chMaxJtpt0p05 = (TH1F*)fData->Get("hpbpb_Jet65_chMaxJtpt0p05");
  TH1F* hData_Jet65_eMaxJtpt0p1 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p1");
  TH1F* hData_Jet65_eMaxJtpt0p2 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p2");
  TH1F* hData_Jet65_eMaxJtpt0p3 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p3");
  TH1F* hData_Jet65_eMaxJtpt0p4 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p4");
  TH1F* hData_Jet65_eMaxJtpt0p5 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5");
  TH1F* hData_Jet65_eMaxJtpt0p6 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6");
  TH1F* hData_Jet65_eMaxJtpt0p7 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7");
  TH1F* hData_Jet65_eMaxJtpt0p8 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p8");
  TH1F* hData_Jet65_eMaxJtpt0p9 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p9");

  TH1F* hData_Jet65_eMaxSumcand0p1 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p1");
  TH1F* hData_Jet65_eMaxSumcand0p2 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p2");
  TH1F* hData_Jet65_eMaxSumcand0p3 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p3");
  TH1F* hData_Jet65_eMaxSumcand0p4 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p4");
  TH1F* hData_Jet65_eMaxSumcand0p5 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p5");
  TH1F* hData_Jet65_eMaxSumcand0p6 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p6");
  TH1F* hData_Jet65_eMaxSumcand0p7 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p7");
  TH1F* hData_Jet65_eMaxSumcand0p8 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p8");
  TH1F* hData_Jet65_eMaxSumcand0p9 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxSumcand0p9");
  TH1F* hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hData_Jet55 = (TH1F*)fData->Get("hpbpb_Jet55");
  TH1F* hData_Jet55_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p01");
  TH1F* hData_Jet55_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p02");
  TH1F* hData_Jet55_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p03");
  TH1F* hData_Jet55_chMaxJtpt0p04 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p04");
  TH1F* hData_Jet55_chMaxJtpt0p05 = (TH1F*)fData->Get("hpbpb_Jet55_chMaxJtpt0p05");
  TH1F* hData_Jet55_eMaxJtpt0p1 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p1");
  TH1F* hData_Jet55_eMaxJtpt0p2 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p2");
  TH1F* hData_Jet55_eMaxJtpt0p3 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p3");
  TH1F* hData_Jet55_eMaxJtpt0p4 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p4");
  TH1F* hData_Jet55_eMaxJtpt0p5 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5");
  TH1F* hData_Jet55_eMaxJtpt0p6 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6");
  TH1F* hData_Jet55_eMaxJtpt0p7 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7");
  TH1F* hData_Jet55_eMaxJtpt0p8 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p8");
  TH1F* hData_Jet55_eMaxJtpt0p9 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p9");

  TH1F* hData_Jet55_eMaxSumcand0p1 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p1");
  TH1F* hData_Jet55_eMaxSumcand0p2 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p2");
  TH1F* hData_Jet55_eMaxSumcand0p3 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p3");
  TH1F* hData_Jet55_eMaxSumcand0p4 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p4");
  TH1F* hData_Jet55_eMaxSumcand0p5 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p5");
  TH1F* hData_Jet55_eMaxSumcand0p6 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p6");
  TH1F* hData_Jet55_eMaxSumcand0p7 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p7");
  TH1F* hData_Jet55_eMaxSumcand0p8 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p8");
  TH1F* hData_Jet55_eMaxSumcand0p9 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxSumcand0p9");
  TH1F* hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fData->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03");

  //TH2F *hData_chMaxJtpt_jtpt = (TH2F*)fData->Get("hpbpb_chMaxJtpt_jtpt");
  //TH2F *hData_eMaxJtpt_jtpt = (TH2F*)fData->Get("hpbpb_eMaxJtpt_jtpt");
  //TH2F *hData_eMaxJtpt_chMaxJtpt = (TH2F*)fData->Get("hpbpb_eMaxJtpt_chMaxJtpt");


  //get the histograms: 
  TH1F* hMC_Jet80 = (TH1F*)fMC->Get("hpbpb_Jet80");
  TH1F* hMC_Jet80_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p01");
  TH1F* hMC_Jet80_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p02");
  TH1F* hMC_Jet80_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p03");
  TH1F* hMC_Jet80_chMaxJtpt0p04 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p04");
  TH1F* hMC_Jet80_chMaxJtpt0p05 = (TH1F*)fMC->Get("hpbpb_Jet80_chMaxJtpt0p05");

  TH1F* hMC_Jet80_eMaxSumcand0p1 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p1");
  TH1F* hMC_Jet80_eMaxSumcand0p2 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p2");
  TH1F* hMC_Jet80_eMaxSumcand0p3 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p3");
  TH1F* hMC_Jet80_eMaxSumcand0p4 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p4");
  TH1F* hMC_Jet80_eMaxSumcand0p5 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p5");
  TH1F* hMC_Jet80_eMaxSumcand0p6 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p6");
  TH1F* hMC_Jet80_eMaxSumcand0p7 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p7");
  TH1F* hMC_Jet80_eMaxSumcand0p8 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p8");
  TH1F* hMC_Jet80_eMaxSumcand0p9 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxSumcand0p9");
  TH1F* hMC_Jet80_eMaxJtpt0p1 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p1");
  TH1F* hMC_Jet80_eMaxJtpt0p2 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p2");
  TH1F* hMC_Jet80_eMaxJtpt0p3 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p3");
  TH1F* hMC_Jet80_eMaxJtpt0p4 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p4");
  TH1F* hMC_Jet80_eMaxJtpt0p5 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5");
  TH1F* hMC_Jet80_eMaxJtpt0p6 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6");
  TH1F* hMC_Jet80_eMaxJtpt0p7 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7");
  TH1F* hMC_Jet80_eMaxJtpt0p8 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p8");
  TH1F* hMC_Jet80_eMaxJtpt0p9 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p9");
  TH1F* hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet80_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hMC_Jet65 = (TH1F*)fMC->Get("hpbpb_Jet65");
  TH1F* hMC_Jet65_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p01");
  TH1F* hMC_Jet65_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p02");
  TH1F* hMC_Jet65_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p03");
  TH1F* hMC_Jet65_chMaxJtpt0p04 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p04");
  TH1F* hMC_Jet65_chMaxJtpt0p05 = (TH1F*)fMC->Get("hpbpb_Jet65_chMaxJtpt0p05");
  TH1F* hMC_Jet65_eMaxJtpt0p1 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p1");
  TH1F* hMC_Jet65_eMaxJtpt0p2 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p2");
  TH1F* hMC_Jet65_eMaxJtpt0p3 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p3");
  TH1F* hMC_Jet65_eMaxJtpt0p4 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p4");
  TH1F* hMC_Jet65_eMaxJtpt0p5 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5");
  TH1F* hMC_Jet65_eMaxJtpt0p6 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6");
  TH1F* hMC_Jet65_eMaxJtpt0p7 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7");
  TH1F* hMC_Jet65_eMaxJtpt0p8 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p8");
  TH1F* hMC_Jet65_eMaxJtpt0p9 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p9");

  TH1F* hMC_Jet65_eMaxSumcand0p1 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p1");
  TH1F* hMC_Jet65_eMaxSumcand0p2 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p2");
  TH1F* hMC_Jet65_eMaxSumcand0p3 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p3");
  TH1F* hMC_Jet65_eMaxSumcand0p4 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p4");
  TH1F* hMC_Jet65_eMaxSumcand0p5 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p5");
  TH1F* hMC_Jet65_eMaxSumcand0p6 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p6");
  TH1F* hMC_Jet65_eMaxSumcand0p7 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p7");
  TH1F* hMC_Jet65_eMaxSumcand0p8 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p8");
  TH1F* hMC_Jet65_eMaxSumcand0p9 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxSumcand0p9");
  TH1F* hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet65_eMaxJtpt0p5_chMaxJtpt0p03");

  TH1F* hMC_Jet55 = (TH1F*)fMC->Get("hpbpb_Jet55");
  TH1F* hMC_Jet55_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p01");
  TH1F* hMC_Jet55_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p02");
  TH1F* hMC_Jet55_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p03");
  TH1F* hMC_Jet55_chMaxJtpt0p04 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p04");
  TH1F* hMC_Jet55_chMaxJtpt0p05 = (TH1F*)fMC->Get("hpbpb_Jet55_chMaxJtpt0p05");
  TH1F* hMC_Jet55_eMaxJtpt0p1 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p1");
  TH1F* hMC_Jet55_eMaxJtpt0p2 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p2");
  TH1F* hMC_Jet55_eMaxJtpt0p3 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p3");
  TH1F* hMC_Jet55_eMaxJtpt0p4 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p4");
  TH1F* hMC_Jet55_eMaxJtpt0p5 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5");
  TH1F* hMC_Jet55_eMaxJtpt0p6 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6");
  TH1F* hMC_Jet55_eMaxJtpt0p7 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7");
  TH1F* hMC_Jet55_eMaxJtpt0p8 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p8");
  TH1F* hMC_Jet55_eMaxJtpt0p9 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p9");

  TH1F* hMC_Jet55_eMaxSumcand0p1 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p1");
  TH1F* hMC_Jet55_eMaxSumcand0p2 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p2");
  TH1F* hMC_Jet55_eMaxSumcand0p3 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p3");
  TH1F* hMC_Jet55_eMaxSumcand0p4 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p4");
  TH1F* hMC_Jet55_eMaxSumcand0p5 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p5");
  TH1F* hMC_Jet55_eMaxSumcand0p6 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p6");
  TH1F* hMC_Jet55_eMaxSumcand0p7 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p7");
  TH1F* hMC_Jet55_eMaxSumcand0p8 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p8");
  TH1F* hMC_Jet55_eMaxSumcand0p9 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxSumcand0p9");
  TH1F* hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p01");
  TH1F* hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p02");
  TH1F* hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p7_chMaxJtpt0p03");
  TH1F* hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p01");
  TH1F* hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p02");
  TH1F* hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p6_chMaxJtpt0p03");
  TH1F* hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p01");
  TH1F* hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p02");
  TH1F* hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03 = (TH1F*)fMC->Get("hpbpb_Jet55_eMaxJtpt0p5_chMaxJtpt0p03");

  //TH2F *hMC_chMaxJtpt_jtpt = (TH2F*)fMC->Get("hpbpb_chMaxJtpt_jtpt");
  //TH2F *hMC_eMaxJtpt_jtpt = (TH2F*)fMC->Get("hpbpb_eMaxJtpt_jtpt");
  //TH2F *hMC_eMaxJtpt_chMaxJtpt = (TH2F*)fMC->Get("hpbpb_eMaxJtpt_chMaxJtpt");


  hData_Jet80_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p03->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p04->Divide(hData_Jet80);
  hData_Jet80_chMaxJtpt0p05->Divide(hData_Jet80);

  hData_Jet80_eMaxJtpt0p1->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p2->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p3->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p4->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p7->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p8->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p9->Divide(hData_Jet80);


  hData_Jet80_eMaxSumcand0p1->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p2->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p3->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p4->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p5->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p6->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p7->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p8->Divide(hData_Jet80);
  hData_Jet80_eMaxSumcand0p9->Divide(hData_Jet80);


  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hData_Jet80);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hData_Jet80);

  hData_Jet65_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p03->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p04->Divide(hData_Jet65);
  hData_Jet65_chMaxJtpt0p05->Divide(hData_Jet65);

  hData_Jet65_eMaxSumcand0p1->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p2->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p3->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p4->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p5->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p6->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p7->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p8->Divide(hData_Jet65);
  hData_Jet65_eMaxSumcand0p9->Divide(hData_Jet65);

  hData_Jet65_eMaxJtpt0p1->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p2->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p3->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p4->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p7->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p8->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p9->Divide(hData_Jet65);

  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hData_Jet65);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hData_Jet65);

  hData_Jet55_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p03->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p04->Divide(hData_Jet55);
  hData_Jet55_chMaxJtpt0p05->Divide(hData_Jet55);

  hData_Jet55_eMaxSumcand0p1->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p2->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p3->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p4->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p5->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p6->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p7->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p8->Divide(hData_Jet55);
  hData_Jet55_eMaxSumcand0p9->Divide(hData_Jet55);

  hData_Jet55_eMaxJtpt0p1->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p2->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p3->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p4->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p7->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p8->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p9->Divide(hData_Jet55);

  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Divide(hData_Jet55);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Divide(hData_Jet55);

  TLine *line = new TLine(30,1,300,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TLine *line3per = new TLine(30,0.97,300,0.97);
  line3per->SetLineStyle(2);
  line3per->SetLineWidth(2);

  //make the plot from the sumcand - 2x2 plot
  TCanvas *cSumcandData = new TCanvas("cSumcandData","",1200,1000);
  makeMultiPanelCanvas(cSumcandData,2,2,0.0,0.0,0.2,0.15,0.07);
  cSumcandData->cd(1);
  
  hData_Jet55_eMaxSumcand0p1->SetTitle(" ");
  hData_Jet55_eMaxSumcand0p1->SetXTitle("Jet p_{T} (GeV/c)");
  hData_Jet55_eMaxSumcand0p1->SetYTitle(" Cut efficiency with Cut/without Cut ");
  hData_Jet55_eMaxSumcand0p1->Rebin(5);
  hData_Jet55_eMaxSumcand0p1->Scale(1./5);
  hData_Jet55_eMaxSumcand0p1->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p1->SetMarkerColor(1);
  hData_Jet55_eMaxSumcand0p1->SetAxisRange(30,300,"X");
  hData_Jet55_eMaxSumcand0p1->SetAxisRange(0,1.2,"Y");
  hData_Jet55_eMaxSumcand0p1->Draw();
  drawText("Data 0-30%",0.5,0.8,14);
  drawText("Jet 55",0.2,0.8,14);

  hData_Jet55_eMaxSumcand0p2->Rebin(5);
  hData_Jet55_eMaxSumcand0p2->Scale(1./5);
  hData_Jet55_eMaxSumcand0p2->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p2->SetMarkerColor(2);
  hData_Jet55_eMaxSumcand0p2->Draw("same");

  hData_Jet55_eMaxSumcand0p3->Rebin(5);
  hData_Jet55_eMaxSumcand0p3->Scale(1./5);
  hData_Jet55_eMaxSumcand0p3->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p3->SetMarkerColor(3);
  hData_Jet55_eMaxSumcand0p3->Draw("same");

  hData_Jet55_eMaxSumcand0p4->Rebin(5);
  hData_Jet55_eMaxSumcand0p4->Scale(1./5);
  hData_Jet55_eMaxSumcand0p4->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p4->SetMarkerColor(4);
  hData_Jet55_eMaxSumcand0p4->Draw("same");

  hData_Jet55_eMaxSumcand0p5->Rebin(5);
  hData_Jet55_eMaxSumcand0p5->Scale(1./5);
  hData_Jet55_eMaxSumcand0p5->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p5->SetMarkerColor(5);
  hData_Jet55_eMaxSumcand0p5->Draw("same");

  hData_Jet55_eMaxSumcand0p6->Rebin(5);
  hData_Jet55_eMaxSumcand0p6->Scale(1./5);
  hData_Jet55_eMaxSumcand0p6->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p6->SetMarkerColor(6);
  hData_Jet55_eMaxSumcand0p6->Draw("same");

  hData_Jet55_eMaxSumcand0p7->Rebin(5);
  hData_Jet55_eMaxSumcand0p7->Scale(1./5);
  hData_Jet55_eMaxSumcand0p7->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p7->SetMarkerColor(7);
  hData_Jet55_eMaxSumcand0p7->Draw("same");

  hData_Jet55_eMaxSumcand0p8->Rebin(5);
  hData_Jet55_eMaxSumcand0p8->Scale(1./5);
  hData_Jet55_eMaxSumcand0p8->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p8->SetMarkerColor(8);
  hData_Jet55_eMaxSumcand0p8->Draw("same");

  hData_Jet55_eMaxSumcand0p9->Rebin(5);
  hData_Jet55_eMaxSumcand0p9->Scale(1./5);
  hData_Jet55_eMaxSumcand0p9->SetMarkerStyle(33);
  hData_Jet55_eMaxSumcand0p9->SetMarkerColor(9);
  hData_Jet55_eMaxSumcand0p9->Draw("same");
  line->Draw();

  cSumcandData->cd(2);

  hData_Jet65_eMaxSumcand0p1->SetTitle(" ");
  hData_Jet65_eMaxSumcand0p1->SetXTitle("Jet p_{T} (GeV/c)");
  hData_Jet65_eMaxSumcand0p1->SetYTitle(" Cut efficiency with Cut/without Cut ");
  hData_Jet65_eMaxSumcand0p1->Rebin(5);
  hData_Jet65_eMaxSumcand0p1->Scale(1./5);
  hData_Jet65_eMaxSumcand0p1->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p1->SetMarkerColor(1);
  hData_Jet65_eMaxSumcand0p1->SetAxisRange(30,300,"X");
  hData_Jet65_eMaxSumcand0p1->SetAxisRange(0,1.2,"Y");
  hData_Jet65_eMaxSumcand0p1->Draw();
  drawText("Jet65",0.2,0.8,14);

  hData_Jet65_eMaxSumcand0p2->Rebin(5);
  hData_Jet65_eMaxSumcand0p2->Scale(1./5);
  hData_Jet65_eMaxSumcand0p2->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p2->SetMarkerColor(2);
  hData_Jet65_eMaxSumcand0p2->Draw("same");

  hData_Jet65_eMaxSumcand0p3->Rebin(5);
  hData_Jet65_eMaxSumcand0p3->Scale(1./5);
  hData_Jet65_eMaxSumcand0p3->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p3->SetMarkerColor(3);
  hData_Jet65_eMaxSumcand0p3->Draw("same");

  hData_Jet65_eMaxSumcand0p4->Rebin(5);
  hData_Jet65_eMaxSumcand0p4->Scale(1./5);
  hData_Jet65_eMaxSumcand0p4->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p4->SetMarkerColor(4);
  hData_Jet65_eMaxSumcand0p4->Draw("same");

  hData_Jet65_eMaxSumcand0p5->Rebin(5);
  hData_Jet65_eMaxSumcand0p5->Scale(1./5);
  hData_Jet65_eMaxSumcand0p5->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p5->SetMarkerColor(5);
  hData_Jet65_eMaxSumcand0p5->Draw("same");

  hData_Jet65_eMaxSumcand0p6->Rebin(5);
  hData_Jet65_eMaxSumcand0p6->Scale(1./5);
  hData_Jet65_eMaxSumcand0p6->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p6->SetMarkerColor(6);
  hData_Jet65_eMaxSumcand0p6->Draw("same");

  hData_Jet65_eMaxSumcand0p7->Rebin(5);
  hData_Jet65_eMaxSumcand0p7->Scale(1./5);
  hData_Jet65_eMaxSumcand0p7->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p7->SetMarkerColor(7);
  hData_Jet65_eMaxSumcand0p7->Draw("same");

  hData_Jet65_eMaxSumcand0p8->Rebin(5);
  hData_Jet65_eMaxSumcand0p8->Scale(1./5);
  hData_Jet65_eMaxSumcand0p8->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p8->SetMarkerColor(8);
  hData_Jet65_eMaxSumcand0p8->Draw("same");

  hData_Jet65_eMaxSumcand0p9->Rebin(5);
  hData_Jet65_eMaxSumcand0p9->Scale(1./5);
  hData_Jet65_eMaxSumcand0p9->SetMarkerStyle(33);
  hData_Jet65_eMaxSumcand0p9->SetMarkerColor(9);
  hData_Jet65_eMaxSumcand0p9->Draw("same");
  line->Draw();
  
  cSumcandData->cd(3);

  hData_Jet80_eMaxSumcand0p1->SetTitle(" ");
  hData_Jet80_eMaxSumcand0p1->SetXTitle("Jet p_{T} (GeV/c)");
  hData_Jet80_eMaxSumcand0p1->SetYTitle(" Cut efficiency with Cut/without Cut ");
  hData_Jet80_eMaxSumcand0p1->Rebin(5);
  hData_Jet80_eMaxSumcand0p1->Scale(1./5);
  hData_Jet80_eMaxSumcand0p1->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p1->SetMarkerColor(1);
  hData_Jet80_eMaxSumcand0p1->SetAxisRange(30,300,"X");
  hData_Jet80_eMaxSumcand0p1->SetAxisRange(0,1.2,"Y");
  hData_Jet80_eMaxSumcand0p1->Draw();
  drawText("Jet80",0.2,0.8,14);

  hData_Jet80_eMaxSumcand0p2->Rebin(5);
  hData_Jet80_eMaxSumcand0p2->Scale(1./5);
  hData_Jet80_eMaxSumcand0p2->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p2->SetMarkerColor(2);
  hData_Jet80_eMaxSumcand0p2->Draw("same");

  hData_Jet80_eMaxSumcand0p3->Rebin(5);
  hData_Jet80_eMaxSumcand0p3->Scale(1./5);
  hData_Jet80_eMaxSumcand0p3->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p3->SetMarkerColor(3);
  hData_Jet80_eMaxSumcand0p3->Draw("same");

  hData_Jet80_eMaxSumcand0p4->Rebin(5);
  hData_Jet80_eMaxSumcand0p4->Scale(1./5);
  hData_Jet80_eMaxSumcand0p4->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p4->SetMarkerColor(4);
  hData_Jet80_eMaxSumcand0p4->Draw("same");

  hData_Jet80_eMaxSumcand0p5->Rebin(5);
  hData_Jet80_eMaxSumcand0p5->Scale(1./5);
  hData_Jet80_eMaxSumcand0p5->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p5->SetMarkerColor(5);
  hData_Jet80_eMaxSumcand0p5->Draw("same");

  hData_Jet80_eMaxSumcand0p6->Rebin(5);
  hData_Jet80_eMaxSumcand0p6->Scale(1./5);
  hData_Jet80_eMaxSumcand0p6->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p6->SetMarkerColor(6);
  hData_Jet80_eMaxSumcand0p6->Draw("same");

  hData_Jet80_eMaxSumcand0p7->Rebin(5);
  hData_Jet80_eMaxSumcand0p7->Scale(1./5);
  hData_Jet80_eMaxSumcand0p7->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p7->SetMarkerColor(7);
  hData_Jet80_eMaxSumcand0p7->Draw("same");

  hData_Jet80_eMaxSumcand0p8->Rebin(5);
  hData_Jet80_eMaxSumcand0p8->Scale(1./5);
  hData_Jet80_eMaxSumcand0p8->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p8->SetMarkerColor(8);
  hData_Jet80_eMaxSumcand0p8->Draw("same");

  hData_Jet80_eMaxSumcand0p9->Rebin(5);
  hData_Jet80_eMaxSumcand0p9->Scale(1./5);
  hData_Jet80_eMaxSumcand0p9->SetMarkerStyle(33);
  hData_Jet80_eMaxSumcand0p9->SetMarkerColor(9);
  hData_Jet80_eMaxSumcand0p9->Draw("same");
  line->Draw();

  cSumcandData->cd(4);
  drawText("#frac{electron Max}{#Sigma h^{#pm} + #Sigma #gamma + #Sigma h^{0} + #Sigma #mu} < X ",0.2,0.7,16);
  TLegend *Jet55_SumcandData = myLegend(0.2,0.2,0.6,0.6);
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p1," 0.1 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p2," 0.2 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p3," 0.3 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p4," 0.4 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p5," 0.5 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p6," 0.6 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p7," 0.7 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p8," 0.8 ","pl");
  Jet55_SumcandData->AddEntry(hData_Jet55_eMaxSumcand0p9," 0.9 ","pl");
  Jet55_SumcandData->Draw();

  cSumcandData->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_data_eMaxSumcand_ratio_plot.pdf","RECREATE");

  //make the plot from the sumcand - 2x2 plot
  TCanvas *cSumcandMC = new TCanvas("cSumcandMC","",1200,1000);
  makeMultiPanelCanvas(cSumcandMC,2,2,0.0,0.0,0.2,0.15,0.07);
  cSumcandMC->cd(1);
  
  hMC_Jet55_eMaxSumcand0p1->SetTitle(" ");
  hMC_Jet55_eMaxSumcand0p1->SetXTitle("Jet p_{T} (GeV/c)");
  hMC_Jet55_eMaxSumcand0p1->SetYTitle(" Cut efficiency with Cut/without Cut ");
  hMC_Jet55_eMaxSumcand0p1->Rebin(5);
  hMC_Jet55_eMaxSumcand0p1->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p1->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p1->SetMarkerColor(1);
  hMC_Jet55_eMaxSumcand0p1->SetAxisRange(30,300,"X");
  hMC_Jet55_eMaxSumcand0p1->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_eMaxSumcand0p1->Draw();
  drawText("MC 0-30%",0.5,0.8,14);
  drawText("Jet 55",0.2,0.8,14);

  hMC_Jet55_eMaxSumcand0p2->Rebin(5);
  hMC_Jet55_eMaxSumcand0p2->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p2->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p2->SetMarkerColor(2);
  hMC_Jet55_eMaxSumcand0p2->Draw("same");

  hMC_Jet55_eMaxSumcand0p3->Rebin(5);
  hMC_Jet55_eMaxSumcand0p3->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p3->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p3->SetMarkerColor(3);
  hMC_Jet55_eMaxSumcand0p3->Draw("same");

  hMC_Jet55_eMaxSumcand0p4->Rebin(5);
  hMC_Jet55_eMaxSumcand0p4->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p4->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p4->SetMarkerColor(4);
  hMC_Jet55_eMaxSumcand0p4->Draw("same");

  hMC_Jet55_eMaxSumcand0p5->Rebin(5);
  hMC_Jet55_eMaxSumcand0p5->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p5->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p5->SetMarkerColor(5);
  hMC_Jet55_eMaxSumcand0p5->Draw("same");

  hMC_Jet55_eMaxSumcand0p6->Rebin(5);
  hMC_Jet55_eMaxSumcand0p6->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p6->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p6->SetMarkerColor(6);
  hMC_Jet55_eMaxSumcand0p6->Draw("same");

  hMC_Jet55_eMaxSumcand0p7->Rebin(5);
  hMC_Jet55_eMaxSumcand0p7->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p7->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p7->SetMarkerColor(7);
  hMC_Jet55_eMaxSumcand0p7->Draw("same");

  hMC_Jet55_eMaxSumcand0p8->Rebin(5);
  hMC_Jet55_eMaxSumcand0p8->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p8->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p8->SetMarkerColor(8);
  hMC_Jet55_eMaxSumcand0p8->Draw("same");

  hMC_Jet55_eMaxSumcand0p9->Rebin(5);
  hMC_Jet55_eMaxSumcand0p9->Scale(1./5);
  hMC_Jet55_eMaxSumcand0p9->SetMarkerStyle(33);
  hMC_Jet55_eMaxSumcand0p9->SetMarkerColor(9);
  hMC_Jet55_eMaxSumcand0p9->Draw("same");
  line->Draw();

  cSumcandMC->cd(2);

  hMC_Jet65_eMaxSumcand0p1->SetTitle(" ");
  hMC_Jet65_eMaxSumcand0p1->SetXTitle("Jet p_{T} (GeV/c)");
  hMC_Jet65_eMaxSumcand0p1->SetYTitle(" Cut efficiency with Cut/without Cut ");
  hMC_Jet65_eMaxSumcand0p1->Rebin(5);
  hMC_Jet65_eMaxSumcand0p1->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p1->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p1->SetMarkerColor(1);
  hMC_Jet65_eMaxSumcand0p1->SetAxisRange(30,300,"X");
  hMC_Jet65_eMaxSumcand0p1->SetAxisRange(0,1.2,"Y");
  hMC_Jet65_eMaxSumcand0p1->Draw();
  drawText("Jet65",0.2,0.8,14);

  hMC_Jet65_eMaxSumcand0p2->Rebin(5);
  hMC_Jet65_eMaxSumcand0p2->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p2->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p2->SetMarkerColor(2);
  hMC_Jet65_eMaxSumcand0p2->Draw("same");

  hMC_Jet65_eMaxSumcand0p3->Rebin(5);
  hMC_Jet65_eMaxSumcand0p3->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p3->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p3->SetMarkerColor(3);
  hMC_Jet65_eMaxSumcand0p3->Draw("same");

  hMC_Jet65_eMaxSumcand0p4->Rebin(5);
  hMC_Jet65_eMaxSumcand0p4->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p4->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p4->SetMarkerColor(4);
  hMC_Jet65_eMaxSumcand0p4->Draw("same");

  hMC_Jet65_eMaxSumcand0p5->Rebin(5);
  hMC_Jet65_eMaxSumcand0p5->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p5->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p5->SetMarkerColor(5);
  hMC_Jet65_eMaxSumcand0p5->Draw("same");

  hMC_Jet65_eMaxSumcand0p6->Rebin(5);
  hMC_Jet65_eMaxSumcand0p6->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p6->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p6->SetMarkerColor(6);
  hMC_Jet65_eMaxSumcand0p6->Draw("same");

  hMC_Jet65_eMaxSumcand0p7->Rebin(5);
  hMC_Jet65_eMaxSumcand0p7->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p7->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p7->SetMarkerColor(7);
  hMC_Jet65_eMaxSumcand0p7->Draw("same");

  hMC_Jet65_eMaxSumcand0p8->Rebin(5);
  hMC_Jet65_eMaxSumcand0p8->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p8->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p8->SetMarkerColor(8);
  hMC_Jet65_eMaxSumcand0p8->Draw("same");

  hMC_Jet65_eMaxSumcand0p9->Rebin(5);
  hMC_Jet65_eMaxSumcand0p9->Scale(1./5);
  hMC_Jet65_eMaxSumcand0p9->SetMarkerStyle(33);
  hMC_Jet65_eMaxSumcand0p9->SetMarkerColor(9);
  hMC_Jet65_eMaxSumcand0p9->Draw("same");
  line->Draw();
  
  cSumcandMC->cd(3);

  hMC_Jet80_eMaxSumcand0p1->SetTitle(" ");
  hMC_Jet80_eMaxSumcand0p1->SetXTitle("Jet p_{T} (GeV/c)");
  hMC_Jet80_eMaxSumcand0p1->SetYTitle(" Cut efficiency with Cut/without Cut ");
  hMC_Jet80_eMaxSumcand0p1->Rebin(5);
  hMC_Jet80_eMaxSumcand0p1->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p1->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p1->SetMarkerColor(1);
  hMC_Jet80_eMaxSumcand0p1->SetAxisRange(30,300,"X");
  hMC_Jet80_eMaxSumcand0p1->SetAxisRange(0,1.2,"Y");
  hMC_Jet80_eMaxSumcand0p1->Draw();
  drawText("Jet80",0.2,0.8,14);

  hMC_Jet80_eMaxSumcand0p2->Rebin(5);
  hMC_Jet80_eMaxSumcand0p2->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p2->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p2->SetMarkerColor(2);
  hMC_Jet80_eMaxSumcand0p2->Draw("same");

  hMC_Jet80_eMaxSumcand0p3->Rebin(5);
  hMC_Jet80_eMaxSumcand0p3->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p3->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p3->SetMarkerColor(3);
  hMC_Jet80_eMaxSumcand0p3->Draw("same");

  hMC_Jet80_eMaxSumcand0p4->Rebin(5);
  hMC_Jet80_eMaxSumcand0p4->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p4->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p4->SetMarkerColor(4);
  hMC_Jet80_eMaxSumcand0p4->Draw("same");

  hMC_Jet80_eMaxSumcand0p5->Rebin(5);
  hMC_Jet80_eMaxSumcand0p5->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p5->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p5->SetMarkerColor(5);
  hMC_Jet80_eMaxSumcand0p5->Draw("same");

  hMC_Jet80_eMaxSumcand0p6->Rebin(5);
  hMC_Jet80_eMaxSumcand0p6->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p6->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p6->SetMarkerColor(6);
  hMC_Jet80_eMaxSumcand0p6->Draw("same");

  hMC_Jet80_eMaxSumcand0p7->Rebin(5);
  hMC_Jet80_eMaxSumcand0p7->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p7->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p7->SetMarkerColor(7);
  hMC_Jet80_eMaxSumcand0p7->Draw("same");

  hMC_Jet80_eMaxSumcand0p8->Rebin(5);
  hMC_Jet80_eMaxSumcand0p8->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p8->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p8->SetMarkerColor(8);
  hMC_Jet80_eMaxSumcand0p8->Draw("same");

  hMC_Jet80_eMaxSumcand0p9->Rebin(5);
  hMC_Jet80_eMaxSumcand0p9->Scale(1./5);
  hMC_Jet80_eMaxSumcand0p9->SetMarkerStyle(33);
  hMC_Jet80_eMaxSumcand0p9->SetMarkerColor(9);
  hMC_Jet80_eMaxSumcand0p9->Draw("same");
  line->Draw();

  cSumcandMC->cd(4);
  drawText("#frac{electron Max}{#Sigma h^{#pm} + #Sigma #gamma + #Sigma neutral + #Sigma #mu} < X ",0.2,0.7,16);
  TLegend *Jet55_SumcandMC = myLegend(0.2,0.2,0.6,0.6);
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p1," 0.1 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p2," 0.2 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p3," 0.3 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p4," 0.4 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p5," 0.5 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p6," 0.6 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p7," 0.7 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p8," 0.8 ","pl");
  Jet55_SumcandMC->AddEntry(hMC_Jet55_eMaxSumcand0p9," 0.9 ","pl");
  Jet55_SumcandMC->Draw();

  cSumcandMC->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_MC_eMaxSumcand_ratio_plot.pdf","RECREATE");


  // lets start making the plots 3 plots, data on top MC at the bottom: 
  // 3 triggers x 3 different Jet iD, plots with 3x1 showing the   
  TCanvas * cJet55 = new TCanvas("cJet55","",1200,1000);
  makeMultiPanelCanvas(cJet55,3,3,0.0,0.0,0.2,0.15,0.07);
  
  cJet55->cd(1);
  hData_Jet55_chMaxJtpt0p01->SetTitle(" ");
  hData_Jet55_chMaxJtpt0p01->SetXTitle(" Jet p_{T} (GeV/c)");
  hData_Jet55_chMaxJtpt0p01->SetYTitle(" Cut effiiency (with Cut/without Cut) ");
  hData_Jet55_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet55_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet55_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet55_chMaxJtpt0p01->Draw();
  drawText("Data 0-30%",0.5,0.8,14);

  hData_Jet55_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_chMaxJtpt0p02->Draw("same");

  hData_Jet55_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_chMaxJtpt0p03->Draw("same");

  hData_Jet55_chMaxJtpt0p04->Rebin(5);
  hData_Jet55_chMaxJtpt0p04->Scale(1./5);
  hData_Jet55_chMaxJtpt0p04->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p04->SetMarkerColor(4);
  hData_Jet55_chMaxJtpt0p04->Draw("same");

  hData_Jet55_chMaxJtpt0p05->Rebin(5);
  hData_Jet55_chMaxJtpt0p05->Scale(1./5);
  hData_Jet55_chMaxJtpt0p05->SetMarkerStyle(33);
  hData_Jet55_chMaxJtpt0p05->SetMarkerColor(5);
  hData_Jet55_chMaxJtpt0p05->Draw("same");
  line->Draw();

  cJet55->cd(2);
  hMC_Jet55_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet55_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_chMaxJtpt0p01->Draw();
  drawText("MC 0-30%",0.2,0.8,14); 

  hMC_Jet55_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_chMaxJtpt0p03->Draw("same");

  hMC_Jet55_chMaxJtpt0p04->Rebin(5);
  hMC_Jet55_chMaxJtpt0p04->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p04->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p04->SetMarkerColor(4);
  hMC_Jet55_chMaxJtpt0p04->Draw("same");

  hMC_Jet55_chMaxJtpt0p05->Rebin(5);
  hMC_Jet55_chMaxJtpt0p05->Scale(1./5);
  hMC_Jet55_chMaxJtpt0p05->SetMarkerStyle(33);
  hMC_Jet55_chMaxJtpt0p05->SetMarkerColor(5);
  hMC_Jet55_chMaxJtpt0p05->Draw("same");
  line->Draw();

  cJet55->cd(3);
  drawText("#frac{charged Max}{Jet p_{T}} > X ",0.2,0.7,16);
  drawText("Jet 55 trigger (does not include higher triggers)", 0.1,0.8,15);
  TLegend *Jet55_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p01,"0.01","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p02,"0.02","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p03,"0.03","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p04,"0.04","pl");
  Jet55_chMax->AddEntry(hData_Jet55_chMaxJtpt0p05,"0.05","pl");
  Jet55_chMax->Draw();

  cJet55->cd(4);
  hData_Jet55_eMaxJtpt0p1->Rebin(5);
  hData_Jet55_eMaxJtpt0p1->Scale(1./5);
  hData_Jet55_eMaxJtpt0p1->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p1->SetMarkerColor(1);
  hData_Jet55_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hData_Jet55_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hData_Jet55_eMaxJtpt0p1->Draw();

  hData_Jet55_eMaxJtpt0p2->Rebin(5);
  hData_Jet55_eMaxJtpt0p2->Scale(1./5);
  hData_Jet55_eMaxJtpt0p2->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p2->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p2->Draw("same");

  hData_Jet55_eMaxJtpt0p3->Rebin(5);
  hData_Jet55_eMaxJtpt0p3->Scale(1./5);
  hData_Jet55_eMaxJtpt0p3->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p3->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p3->Draw("same");

  hData_Jet55_eMaxJtpt0p4->Rebin(5);
  hData_Jet55_eMaxJtpt0p4->Scale(1./5);
  hData_Jet55_eMaxJtpt0p4->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p4->SetMarkerColor(4);
  hData_Jet55_eMaxJtpt0p4->Draw("same");

  hData_Jet55_eMaxJtpt0p5->Rebin(5);
  hData_Jet55_eMaxJtpt0p5->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p5->SetMarkerColor(5);
  hData_Jet55_eMaxJtpt0p5->Draw("same");

  hData_Jet55_eMaxJtpt0p6->Rebin(5);
  hData_Jet55_eMaxJtpt0p6->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p6->SetMarkerColor(6);
  hData_Jet55_eMaxJtpt0p6->Draw("same");

  hData_Jet55_eMaxJtpt0p7->Rebin(5);
  hData_Jet55_eMaxJtpt0p7->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7->SetMarkerColor(7);
  hData_Jet55_eMaxJtpt0p7->Draw("same");

  hData_Jet55_eMaxJtpt0p8->Rebin(5);
  hData_Jet55_eMaxJtpt0p8->Scale(1./5);
  hData_Jet55_eMaxJtpt0p8->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p8->SetMarkerColor(8);
  hData_Jet55_eMaxJtpt0p8->Draw("same");

  hData_Jet55_eMaxJtpt0p9->Rebin(5);
  hData_Jet55_eMaxJtpt0p9->Scale(1./5);
  hData_Jet55_eMaxJtpt0p9->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p9->SetMarkerColor(9);
  hData_Jet55_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet55->cd(5);

  hMC_Jet55_eMaxJtpt0p1->Rebin(5);
  hMC_Jet55_eMaxJtpt0p1->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p1->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p1->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hMC_Jet55_eMaxJtpt0p1->Draw();

  hMC_Jet55_eMaxJtpt0p2->Rebin(5);
  hMC_Jet55_eMaxJtpt0p2->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p2->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p2->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p2->Draw("same");

  hMC_Jet55_eMaxJtpt0p3->Rebin(5);
  hMC_Jet55_eMaxJtpt0p3->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p3->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p3->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p3->Draw("same");

  hMC_Jet55_eMaxJtpt0p4->Rebin(5);
  hMC_Jet55_eMaxJtpt0p4->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p4->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p4->SetMarkerColor(4);
  hMC_Jet55_eMaxJtpt0p4->Draw("same");

  hMC_Jet55_eMaxJtpt0p5->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p5->SetMarkerColor(5);
  hMC_Jet55_eMaxJtpt0p5->Draw("same");

  hMC_Jet55_eMaxJtpt0p6->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p6->SetMarkerColor(6);
  hMC_Jet55_eMaxJtpt0p6->Draw("same");

  hMC_Jet55_eMaxJtpt0p7->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7->SetMarkerColor(7);
  hMC_Jet55_eMaxJtpt0p7->Draw("same");

  hMC_Jet55_eMaxJtpt0p8->Rebin(5);
  hMC_Jet55_eMaxJtpt0p8->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p8->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p8->SetMarkerColor(8);
  hMC_Jet55_eMaxJtpt0p8->Draw("same");

  hMC_Jet55_eMaxJtpt0p9->Rebin(5);
  hMC_Jet55_eMaxJtpt0p9->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p9->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p9->SetMarkerColor(9);
  hMC_Jet55_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet55->cd(6);
  drawText("#frac{electron Max}{Jet p_{T}} < X ",0.2,0.7,16);
  TLegend *Jet55_eMax = myLegend(0.2,0.2,0.6,0.6);
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p1," 0.1 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p2," 0.2 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p3," 0.3 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p4," 0.4 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p5," 0.5 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p6," 0.6 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p7," 0.7 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p8," 0.8 ","pl");
  Jet55_eMax->AddEntry(hMC_Jet55_eMaxJtpt0p9," 0.9 ","pl");
  Jet55_eMax->Draw();

  cJet55->cd(7);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetXTitle("Jet p_{T} (GeV/c)");
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);  
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet55->cd(8);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetXTitle("Jet p_{T} (GeV/c)");
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet55->cd(9);
  drawText("#frac{electron Max}{Jet p_{T}} < 0.X && #frac{charged Max}{Jet p_{T}} > 0.0Y",0.2,0.7,16);
  TLegend *Jet55_eMax_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p01,"X = 7, Y = 1","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p02,"X = 7, Y = 2","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p7_chMaxJtpt0p03,"X = 7, Y = 3","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p01,"X = 6, Y = 1","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p02,"X = 6, Y = 2","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p6_chMaxJtpt0p03,"X = 6, Y = 3","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p01,"X = 5, Y = 1","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p02,"X = 5, Y = 2","pl");
  Jet55_eMax_chMax->AddEntry(hMC_Jet55_eMaxJtpt0p5_chMaxJtpt0p03,"X = 5, Y = 3","pl");
  Jet55_eMax_chMax->Draw();

  cJet55->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jet55_spectra_ratio_jetid.pdf","RECREATE");



  TCanvas * cJet65 = new TCanvas("cJet65","",1200,1000);
  makeMultiPanelCanvas(cJet65,3,3,0.0,0.0,0.2,0.15,0.07);

  cJet65->cd(1);
  hData_Jet65_chMaxJtpt0p01->SetTitle(" ");
  hData_Jet65_chMaxJtpt0p01->SetXTitle(" Jet p_{T} (GeV/c)");
  hData_Jet65_chMaxJtpt0p01->SetYTitle(" Cut effiiency (with Cut/without Cut) ");
  hData_Jet65_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet65_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet65_chMaxJtpt0p01->Draw();

  hData_Jet65_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_chMaxJtpt0p02->Draw("same");

  hData_Jet65_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_chMaxJtpt0p03->Draw("same");

  hData_Jet65_chMaxJtpt0p04->Rebin(5);
  hData_Jet65_chMaxJtpt0p04->Scale(1./5);
  hData_Jet65_chMaxJtpt0p04->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p04->SetMarkerColor(4);
  hData_Jet65_chMaxJtpt0p04->Draw("same");

  hData_Jet65_chMaxJtpt0p05->Rebin(5);
  hData_Jet65_chMaxJtpt0p05->Scale(1./5);
  hData_Jet65_chMaxJtpt0p05->SetMarkerStyle(33);
  hData_Jet65_chMaxJtpt0p05->SetMarkerColor(5);
  hData_Jet65_chMaxJtpt0p05->Draw("same");
  line->Draw();
  drawText("Data 0-30%",0.5,0.8,16);

  cJet65->cd(2);
  hMC_Jet65_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet65_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet65_chMaxJtpt0p01->Draw();

  hMC_Jet65_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_chMaxJtpt0p03->Draw("same");

  hMC_Jet65_chMaxJtpt0p04->Rebin(5);
  hMC_Jet65_chMaxJtpt0p04->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p04->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p04->SetMarkerColor(4);
  hMC_Jet65_chMaxJtpt0p04->Draw("same");

  hMC_Jet65_chMaxJtpt0p05->Rebin(5);
  hMC_Jet65_chMaxJtpt0p05->Scale(1./5);
  hMC_Jet65_chMaxJtpt0p05->SetMarkerStyle(33);
  hMC_Jet65_chMaxJtpt0p05->SetMarkerColor(5);
  hMC_Jet65_chMaxJtpt0p05->Draw("same");
  line->Draw();
  drawText("MC 0-30%",0.5,0.8,16);

  cJet65->cd(3);
  drawText("Jet 65 Trigger without Jet 80",0.1,0.8,15);
  drawText("#frac{charged Max}{Jet p_{T}} > X ",0.2,0.7,16);
  TLegend *Jet65_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p01,"0.01","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p02,"0.02","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p03,"0.03","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p04,"0.04","pl");
  Jet65_chMax->AddEntry(hData_Jet65_chMaxJtpt0p05,"0.05","pl");
  Jet65_chMax->Draw();

  cJet65->cd(4);
  hData_Jet65_eMaxJtpt0p1->Rebin(5);
  hData_Jet65_eMaxJtpt0p1->Scale(1./5);
  hData_Jet65_eMaxJtpt0p1->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p1->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hData_Jet65_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hData_Jet65_eMaxJtpt0p1->Draw();

  hData_Jet65_eMaxJtpt0p2->Rebin(5);
  hData_Jet65_eMaxJtpt0p2->Scale(1./5);
  hData_Jet65_eMaxJtpt0p2->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p2->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p2->Draw("same");

  hData_Jet65_eMaxJtpt0p3->Rebin(5);
  hData_Jet65_eMaxJtpt0p3->Scale(1./5);
  hData_Jet65_eMaxJtpt0p3->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p3->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p3->Draw("same");

  hData_Jet65_eMaxJtpt0p4->Rebin(5);
  hData_Jet65_eMaxJtpt0p4->Scale(1./5);
  hData_Jet65_eMaxJtpt0p4->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p4->SetMarkerColor(4);
  hData_Jet65_eMaxJtpt0p4->Draw("same");

  hData_Jet65_eMaxJtpt0p5->Rebin(5);
  hData_Jet65_eMaxJtpt0p5->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p5->SetMarkerColor(5);
  hData_Jet65_eMaxJtpt0p5->Draw("same");

  hData_Jet65_eMaxJtpt0p6->Rebin(5);
  hData_Jet65_eMaxJtpt0p6->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p6->SetMarkerColor(6);
  hData_Jet65_eMaxJtpt0p6->Draw("same");

  hData_Jet65_eMaxJtpt0p7->Rebin(5);
  hData_Jet65_eMaxJtpt0p7->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7->SetMarkerColor(7);
  hData_Jet65_eMaxJtpt0p7->Draw("same");

  hData_Jet65_eMaxJtpt0p8->Rebin(5);
  hData_Jet65_eMaxJtpt0p8->Scale(1./5);
  hData_Jet65_eMaxJtpt0p8->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p8->SetMarkerColor(8);
  hData_Jet65_eMaxJtpt0p8->Draw("same");

  hData_Jet65_eMaxJtpt0p9->Rebin(5);
  hData_Jet65_eMaxJtpt0p9->Scale(1./5);
  hData_Jet65_eMaxJtpt0p9->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p9->SetMarkerColor(9);
  hData_Jet65_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet65->cd(5);

  hMC_Jet65_eMaxJtpt0p1->Rebin(5);
  hMC_Jet65_eMaxJtpt0p1->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p1->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p1->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hMC_Jet65_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hMC_Jet65_eMaxJtpt0p1->Draw();

  hMC_Jet65_eMaxJtpt0p2->Rebin(5);
  hMC_Jet65_eMaxJtpt0p2->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p2->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p2->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p2->Draw("same");

  hMC_Jet65_eMaxJtpt0p3->Rebin(5);
  hMC_Jet65_eMaxJtpt0p3->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p3->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p3->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p3->Draw("same");

  hMC_Jet65_eMaxJtpt0p4->Rebin(5);
  hMC_Jet65_eMaxJtpt0p4->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p4->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p4->SetMarkerColor(4);
  hMC_Jet65_eMaxJtpt0p4->Draw("same");

  hMC_Jet65_eMaxJtpt0p5->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p5->SetMarkerColor(5);
  hMC_Jet65_eMaxJtpt0p5->Draw("same");

  hMC_Jet65_eMaxJtpt0p6->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p6->SetMarkerColor(6);
  hMC_Jet65_eMaxJtpt0p6->Draw("same");

  hMC_Jet65_eMaxJtpt0p7->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7->SetMarkerColor(7);
  hMC_Jet65_eMaxJtpt0p7->Draw("same");

  hMC_Jet65_eMaxJtpt0p8->Rebin(5);
  hMC_Jet65_eMaxJtpt0p8->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p8->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p8->SetMarkerColor(8);
  hMC_Jet65_eMaxJtpt0p8->Draw("same");

  hMC_Jet65_eMaxJtpt0p9->Rebin(5);
  hMC_Jet65_eMaxJtpt0p9->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p9->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p9->SetMarkerColor(9);
  hMC_Jet65_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet65->cd(6);
  drawText("#frac{electron Max}{Jet p_{T}} < X ",0.2,0.7,16);
  TLegend *Jet65_eMax = myLegend(0.2,0.2,0.6,0.6);
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p1," 0.1 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p2," 0.2 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p3," 0.3 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p4," 0.4 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p5," 0.5 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p6," 0.6 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p7," 0.7 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p8," 0.8 ","pl");
  Jet65_eMax->AddEntry(hMC_Jet65_eMaxJtpt0p9," 0.9 ","pl");
  Jet65_eMax->Draw();

  cJet65->cd(7);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet65->cd(8);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet65->cd(9);
  drawText("#frac{electron Max}{Jet p_{T}} < 0.X && #frac{charged Max}{Jet p_{T}} > 0.0Y",0.2,0.7,16);
  TLegend *Jet65_eMax_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p01,"X = 7, Y = 1","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p02,"X = 7, Y = 2","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p7_chMaxJtpt0p03,"X = 7, Y = 3","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p01,"X = 6, Y = 1","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p02,"X = 6, Y = 2","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p6_chMaxJtpt0p03,"X = 6, Y = 3","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p01,"X = 5, Y = 1","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p02,"X = 5, Y = 2","pl");
  Jet65_eMax_chMax->AddEntry(hMC_Jet65_eMaxJtpt0p5_chMaxJtpt0p03,"X = 5, Y = 3","pl");
  Jet65_eMax_chMax->Draw();

  cJet65->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jet65_spectra_ratio_jetid.pdf","RECREATE");

  
  TCanvas * cJet80 = new TCanvas("cJet80","",1200,1000);
  makeMultiPanelCanvas(cJet80,3,3,0.0,0.0,0.2,0.15,0.07);

  cJet80->cd(1);
  hData_Jet80_chMaxJtpt0p01->SetTitle(" ");
  hData_Jet80_chMaxJtpt0p01->SetXTitle(" Jet p_{T} (GeV/c)");
  hData_Jet80_chMaxJtpt0p01->SetYTitle(" Cut effiiency (with Cut/without Cut) ");
  hData_Jet80_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet80_chMaxJtpt0p01->SetAxisRange(30,300,"X");  
  hData_Jet80_chMaxJtpt0p01->Draw();

  hData_Jet80_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_chMaxJtpt0p02->Draw("same");

  hData_Jet80_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_chMaxJtpt0p03->Draw("same");

  hData_Jet80_chMaxJtpt0p04->Rebin(5);
  hData_Jet80_chMaxJtpt0p04->Scale(1./5);
  hData_Jet80_chMaxJtpt0p04->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p04->SetMarkerColor(4);
  hData_Jet80_chMaxJtpt0p04->Draw("same");

  hData_Jet80_chMaxJtpt0p05->Rebin(5);
  hData_Jet80_chMaxJtpt0p05->Scale(1./5);
  hData_Jet80_chMaxJtpt0p05->SetMarkerStyle(33);
  hData_Jet80_chMaxJtpt0p05->SetMarkerColor(5);
  hData_Jet80_chMaxJtpt0p05->Draw("same");
  line->Draw();
  drawText("Data 0-30%",0.5,0.8,15);

  cJet80->cd(2);
  hMC_Jet80_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet80_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet80_chMaxJtpt0p01->Draw();

  hMC_Jet80_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_chMaxJtpt0p03->Draw("same");

  hMC_Jet80_chMaxJtpt0p04->Rebin(5);
  hMC_Jet80_chMaxJtpt0p04->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p04->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p04->SetMarkerColor(4);
  hMC_Jet80_chMaxJtpt0p04->Draw("same");

  hMC_Jet80_chMaxJtpt0p05->Rebin(5);
  hMC_Jet80_chMaxJtpt0p05->Scale(1./5);
  hMC_Jet80_chMaxJtpt0p05->SetMarkerStyle(33);
  hMC_Jet80_chMaxJtpt0p05->SetMarkerColor(5);
  hMC_Jet80_chMaxJtpt0p05->Draw("same");
  line->Draw();
  drawText("MC 0-30%",0.5,0.8,15);

  cJet80->cd(3);
  drawText("#frac{charged Max}{Jet p_{T}} > X ",0.2,0.7,16);
  drawText("Jet 80 trigger",0.2,0.8,16);
  TLegend *Jet80_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p01,"0.01","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p02,"0.02","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p03,"0.03","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p04,"0.04","pl");
  Jet80_chMax->AddEntry(hData_Jet80_chMaxJtpt0p05,"0.05","pl");
  Jet80_chMax->Draw();

  cJet80->cd(4);
  hData_Jet80_eMaxJtpt0p1->Rebin(5);
  hData_Jet80_eMaxJtpt0p1->Scale(1./5);
  hData_Jet80_eMaxJtpt0p1->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p1->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hData_Jet80_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hData_Jet80_eMaxJtpt0p1->Draw();

  hData_Jet80_eMaxJtpt0p2->Rebin(5);
  hData_Jet80_eMaxJtpt0p2->Scale(1./5);
  hData_Jet80_eMaxJtpt0p2->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p2->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p2->Draw("same");

  hData_Jet80_eMaxJtpt0p3->Rebin(5);
  hData_Jet80_eMaxJtpt0p3->Scale(1./5);
  hData_Jet80_eMaxJtpt0p3->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p3->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p3->Draw("same");

  hData_Jet80_eMaxJtpt0p4->Rebin(5);
  hData_Jet80_eMaxJtpt0p4->Scale(1./5);
  hData_Jet80_eMaxJtpt0p4->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p4->SetMarkerColor(4);
  hData_Jet80_eMaxJtpt0p4->Draw("same");

  hData_Jet80_eMaxJtpt0p5->Rebin(5);
  hData_Jet80_eMaxJtpt0p5->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p5->SetMarkerColor(5);
  hData_Jet80_eMaxJtpt0p5->Draw("same");

  hData_Jet80_eMaxJtpt0p6->Rebin(5);
  hData_Jet80_eMaxJtpt0p6->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p6->SetMarkerColor(6);
  hData_Jet80_eMaxJtpt0p6->Draw("same");

  hData_Jet80_eMaxJtpt0p7->Rebin(5);
  hData_Jet80_eMaxJtpt0p7->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7->SetMarkerColor(7);
  hData_Jet80_eMaxJtpt0p7->Draw("same");

  hData_Jet80_eMaxJtpt0p8->Rebin(5);
  hData_Jet80_eMaxJtpt0p8->Scale(1./5);
  hData_Jet80_eMaxJtpt0p8->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p8->SetMarkerColor(8);
  hData_Jet80_eMaxJtpt0p8->Draw("same");

  hData_Jet80_eMaxJtpt0p9->Rebin(5);
  hData_Jet80_eMaxJtpt0p9->Scale(1./5);
  hData_Jet80_eMaxJtpt0p9->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p9->SetMarkerColor(9);
  hData_Jet80_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet80->cd(5);

  hMC_Jet80_eMaxJtpt0p1->Rebin(5);
  hMC_Jet80_eMaxJtpt0p1->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p1->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p1->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p1->SetAxisRange(0,1.2,"Y");
  hMC_Jet80_eMaxJtpt0p1->SetAxisRange(30,300,"X");
  hMC_Jet80_eMaxJtpt0p1->Draw();

  hMC_Jet80_eMaxJtpt0p2->Rebin(5);
  hMC_Jet80_eMaxJtpt0p2->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p2->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p2->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p2->Draw("same");

  hMC_Jet80_eMaxJtpt0p3->Rebin(5);
  hMC_Jet80_eMaxJtpt0p3->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p3->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p3->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p3->Draw("same");

  hMC_Jet80_eMaxJtpt0p4->Rebin(5);
  hMC_Jet80_eMaxJtpt0p4->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p4->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p4->SetMarkerColor(4);
  hMC_Jet80_eMaxJtpt0p4->Draw("same");

  hMC_Jet80_eMaxJtpt0p5->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p5->SetMarkerColor(5);
  hMC_Jet80_eMaxJtpt0p5->Draw("same");

  hMC_Jet80_eMaxJtpt0p6->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p6->SetMarkerColor(6);
  hMC_Jet80_eMaxJtpt0p6->Draw("same");

  hMC_Jet80_eMaxJtpt0p7->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7->SetMarkerColor(7);
  hMC_Jet80_eMaxJtpt0p7->Draw("same");

  hMC_Jet80_eMaxJtpt0p8->Rebin(5);
  hMC_Jet80_eMaxJtpt0p8->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p8->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p8->SetMarkerColor(8);
  hMC_Jet80_eMaxJtpt0p8->Draw("same");

  hMC_Jet80_eMaxJtpt0p9->Rebin(5);
  hMC_Jet80_eMaxJtpt0p9->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p9->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p9->SetMarkerColor(9);
  hMC_Jet80_eMaxJtpt0p9->Draw("same");
  line->Draw();

  cJet80->cd(6);
  drawText("#frac{electron Max}{Jet p_{T}} < X ",0.2,0.7,16);
  TLegend *Jet80_eMax = myLegend(0.2,0.2,0.6,0.6);
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p1," 0.1 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p2," 0.2 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p3," 0.3 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p4," 0.4 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p5," 0.5 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p6," 0.6 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p7," 0.7 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p8," 0.8 ","pl");
  Jet80_eMax->AddEntry(hMC_Jet80_eMaxJtpt0p9," 0.9 ","pl");
  Jet80_eMax->Draw();

  cJet80->cd(7);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hData_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet80->cd(8);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(0,1.2,"Y");
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->SetAxisRange(30,300,"X");
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01->Draw();

  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerStyle(33);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03->Draw("same");

  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerStyle(25);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01->Draw("same");

  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerStyle(25);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerStyle(25);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03->Draw("same");

  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerStyle(26);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->SetMarkerColor(1);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01->Draw("same");

  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerStyle(26);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->SetMarkerColor(2);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02->Draw("same");

  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Rebin(5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Scale(1./5);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerStyle(26);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->SetMarkerColor(3);
  hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03->Draw("same");
  line->Draw();

  cJet80->cd(9);
  drawText("#frac{electron Max}{Jet p_{T}} < 0.X && #frac{charged Max}{Jet p_{T}} > 0.0Y",0.2,0.7,16);
  TLegend *Jet80_eMax_chMax = myLegend(0.2,0.2,0.6,0.6);
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p01,"X = 7, Y = 1","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p02,"X = 7, Y = 2","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p7_chMaxJtpt0p03,"X = 7, Y = 3","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p01,"X = 6, Y = 1","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p02,"X = 6, Y = 2","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p6_chMaxJtpt0p03,"X = 6, Y = 3","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p01,"X = 5, Y = 1","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p02,"X = 5, Y = 2","pl");
  Jet80_eMax_chMax->AddEntry(hMC_Jet80_eMaxJtpt0p5_chMaxJtpt0p03,"X = 5, Y = 3","pl");
  Jet80_eMax_chMax->Draw();

  cJet80->SaveAs("/Users/keraghav/WORK/RAA/Plots/PbPb_jet80_spectra_ratio_jetid.pdf","RECREATE");

#endif 
  
  // get the pfcand variables before and after the cut and take ratio:

  TH1F* hMC_chMax[nbins_cent], * hMC_phMax[nbins_cent], * hMC_neMax[nbins_cent], * hMC_muMax[nbins_cent], * hMC_eMax[nbins_cent];
  TH1F* hMC_chSum[nbins_cent], * hMC_phSum[nbins_cent], * hMC_neSum[nbins_cent], * hMC_muSum[nbins_cent], * hMC_eSum[nbins_cent];

  TH1F* hMC_chMax_withCut[nbins_cent], * hMC_phMax_withCut[nbins_cent], * hMC_neMax_withCut[nbins_cent], * hMC_muMax_withCut[nbins_cent], * hMC_eMax_withCut[nbins_cent];
  TH1F* hMC_chSum_withCut[nbins_cent], * hMC_phSum_withCut[nbins_cent], * hMC_neSum_withCut[nbins_cent], * hMC_muSum_withCut[nbins_cent], * hMC_eSum_withCut[nbins_cent];

  TH1F* hData_chMax[nbins_cent], * hData_phMax[nbins_cent], * hData_neMax[nbins_cent], * hData_muMax[nbins_cent], * hData_eMax[nbins_cent];
  TH1F* hData_chSum[nbins_cent], * hData_phSum[nbins_cent], * hData_neSum[nbins_cent], * hData_muSum[nbins_cent], * hData_eSum[nbins_cent];

  TH1F* hData_chMax_withCut[nbins_cent], * hData_phMax_withCut[nbins_cent], * hData_neMax_withCut[nbins_cent], * hData_muMax_withCut[nbins_cent], * hData_eMax_withCut[nbins_cent];
  TH1F* hData_chSum_withCut[nbins_cent], * hData_phSum_withCut[nbins_cent], * hData_neSum_withCut[nbins_cent], * hData_muSum_withCut[nbins_cent], * hData_eSum_withCut[nbins_cent];

  TH1F * hRatio_MC_chMax[nbins_cent], * hRatio_MC_phMax[nbins_cent], * hRatio_MC_neMax[nbins_cent], * hRatio_MC_muMax[nbins_cent], * hRatio_MC_eMax[nbins_cent];
  TH1F * hRatio_MC_chSum[nbins_cent], * hRatio_MC_phSum[nbins_cent], * hRatio_MC_neSum[nbins_cent], * hRatio_MC_muSum[nbins_cent], * hRatio_MC_eSum[nbins_cent];
  TH1F * hRatio_Data_chMax[nbins_cent], * hRatio_Data_phMax[nbins_cent], * hRatio_Data_neMax[nbins_cent], * hRatio_Data_muMax[nbins_cent], * hRatio_Data_eMax[nbins_cent];
  TH1F * hRatio_Data_chSum[nbins_cent], * hRatio_Data_phSum[nbins_cent], * hRatio_Data_neSum[nbins_cent], * hRatio_Data_muSum[nbins_cent], * hRatio_Data_eSum[nbins_cent];
  
  for(int i = 0;i<nbins_cent;++i){

    hMC_chMax[i] = (TH1F*)fMC->Get(Form("hpbpb_chMax_R3_n20_eta_p20_cent%d",i));
    hMC_phMax[i] = (TH1F*)fMC->Get(Form("hpbpb_phMax_R3_n20_eta_p20_cent%d",i));
    hMC_neMax[i] = (TH1F*)fMC->Get(Form("hpbpb_neMax_R3_n20_eta_p20_cent%d",i));
    hMC_muMax[i] = (TH1F*)fMC->Get(Form("hpbpb_muMax_R3_n20_eta_p20_cent%d",i));
    hMC_eMax[i] = (TH1F*)fMC->Get(Form("hpbpb_eMax_R3_n20_eta_p20_cent%d",i));

    hMC_chSum[i] = (TH1F*)fMC->Get(Form("hpbpb_chSum_R3_n20_eta_p20_cent%d",i));
    hMC_phSum[i] = (TH1F*)fMC->Get(Form("hpbpb_phSum_R3_n20_eta_p20_cent%d",i));
    hMC_neSum[i] = (TH1F*)fMC->Get(Form("hpbpb_neSum_R3_n20_eta_p20_cent%d",i));
    hMC_muSum[i] = (TH1F*)fMC->Get(Form("hpbpb_muSum_R3_n20_eta_p20_cent%d",i));
    hMC_eSum[i] = (TH1F*)fMC->Get(Form("hpbpb_eSum_R3_n20_eta_p20_cent%d",i));

    hMC_chMax_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_chMax_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_phMax_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_phMax_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_neMax_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_neMax_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_muMax_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_muMax_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_eMax_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_eMax_withCut_R3_n20_eta_p20_cent%d",i));

    hMC_chSum_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_chSum_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_phSum_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_phSum_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_neSum_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_neSum_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_muSum_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_muSum_withCut_R3_n20_eta_p20_cent%d",i));
    hMC_eSum_withCut[i] = (TH1F*)fMC->Get(Form("hpbpb_eSum_withCut_R3_n20_eta_p20_cent%d",i));

    
    hData_chMax[i] = (TH1F*)fData->Get(Form("hpbpb_chMax_R3_n20_eta_p20_cent%d",i));
    hData_phMax[i] = (TH1F*)fData->Get(Form("hpbpb_phMax_R3_n20_eta_p20_cent%d",i));
    hData_neMax[i] = (TH1F*)fData->Get(Form("hpbpb_neMax_R3_n20_eta_p20_cent%d",i));
    hData_muMax[i] = (TH1F*)fData->Get(Form("hpbpb_muMax_R3_n20_eta_p20_cent%d",i));
    hData_eMax[i] = (TH1F*)fData->Get(Form("hpbpb_eMax_R3_n20_eta_p20_cent%d",i));

    hData_chSum[i] = (TH1F*)fData->Get(Form("hpbpb_chSum_R3_n20_eta_p20_cent%d",i));
    hData_phSum[i] = (TH1F*)fData->Get(Form("hpbpb_phSum_R3_n20_eta_p20_cent%d",i));
    hData_neSum[i] = (TH1F*)fData->Get(Form("hpbpb_neSum_R3_n20_eta_p20_cent%d",i));
    hData_muSum[i] = (TH1F*)fData->Get(Form("hpbpb_muSum_R3_n20_eta_p20_cent%d",i));
    hData_eSum[i] = (TH1F*)fData->Get(Form("hpbpb_eSum_R3_n20_eta_p20_cent%d",i));

    hData_chMax_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_chMax_withCut_R3_n20_eta_p20_cent%d",i));
    hData_phMax_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_phMax_withCut_R3_n20_eta_p20_cent%d",i));
    hData_neMax_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_neMax_withCut_R3_n20_eta_p20_cent%d",i));
    hData_muMax_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_muMax_withCut_R3_n20_eta_p20_cent%d",i));
    hData_eMax_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_eMax_withCut_R3_n20_eta_p20_cent%d",i));

    hData_chSum_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_chSum_withCut_R3_n20_eta_p20_cent%d",i));
    hData_phSum_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_phSum_withCut_R3_n20_eta_p20_cent%d",i));
    hData_neSum_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_neSum_withCut_R3_n20_eta_p20_cent%d",i));
    hData_muSum_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_muSum_withCut_R3_n20_eta_p20_cent%d",i));
    hData_eSum_withCut[i] = (TH1F*)fData->Get(Form("hpbpb_eSum_withCut_R3_n20_eta_p20_cent%d",i));

    hRatio_MC_chMax[i] = (TH1F*)hMC_chMax_withCut[i]->Clone(Form("hRatio_MC_chMax_cent%d",i));
    hRatio_MC_chMax[i]->Divide(hMC_chMax[i]);
    hRatio_MC_phMax[i] = (TH1F*)hMC_phMax_withCut[i]->Clone(Form("hRatio_MC_phMax_cent%d",i));
    hRatio_MC_phMax[i]->Divide(hMC_phMax[i]);
    hRatio_MC_neMax[i] = (TH1F*)hMC_neMax_withCut[i]->Clone(Form("hRatio_MC_neMax_cent%d",i));
    hRatio_MC_neMax[i]->Divide(hMC_neMax[i]);
    hRatio_MC_muMax[i] = (TH1F*)hMC_muMax_withCut[i]->Clone(Form("hRatio_MC_muMax_cent%d",i));
    hRatio_MC_muMax[i]->Divide(hMC_muMax[i]);
    hRatio_MC_eMax[i] = (TH1F*)hMC_eMax_withCut[i]->Clone(Form("hRatio_MC_eMax_cent%d",i));
    hRatio_MC_eMax[i]->Divide(hMC_eMax[i]);

    hRatio_MC_chSum[i] = (TH1F*)hMC_chSum_withCut[i]->Clone(Form("hRatio_MC_chSum_cent%d",i));
    hRatio_MC_chSum[i]->Divide(hMC_chSum[i]);
    hRatio_MC_phSum[i] = (TH1F*)hMC_phSum_withCut[i]->Clone(Form("hRatio_MC_phSum_cent%d",i));
    hRatio_MC_phSum[i]->Divide(hMC_phSum[i]);
    hRatio_MC_neSum[i] = (TH1F*)hMC_neSum_withCut[i]->Clone(Form("hRatio_MC_neSum_cent%d",i));
    hRatio_MC_neSum[i]->Divide(hMC_neSum[i]);
    hRatio_MC_muSum[i] = (TH1F*)hMC_muSum_withCut[i]->Clone(Form("hRatio_MC_muSum_cent%d",i));
    hRatio_MC_muSum[i]->Divide(hMC_muSum[i]);
    hRatio_MC_eSum[i] = (TH1F*)hMC_eSum_withCut[i]->Clone(Form("hRatio_MC_eSum_cent%d",i));
    hRatio_MC_eSum[i]->Divide(hMC_eSum[i]);
    
    hRatio_Data_chMax[i] = (TH1F*)hData_chMax_withCut[i]->Clone(Form("hRatio_Data_chMax_cent%d",i));
    hRatio_Data_chMax[i]->Divide(hData_chMax[i]);
    hRatio_Data_phMax[i] = (TH1F*)hData_phMax_withCut[i]->Clone(Form("hRatio_Data_phMax_cent%d",i));
    hRatio_Data_phMax[i]->Divide(hData_phMax[i]);
    hRatio_Data_neMax[i] = (TH1F*)hData_neMax_withCut[i]->Clone(Form("hRatio_Data_neMax_cent%d",i));
    hRatio_Data_neMax[i]->Divide(hData_neMax[i]);
    hRatio_Data_muMax[i] = (TH1F*)hData_muMax_withCut[i]->Clone(Form("hRatio_Data_muMax_cent%d",i));
    hRatio_Data_muMax[i]->Divide(hData_muMax[i]);
    hRatio_Data_eMax[i] = (TH1F*)hData_eMax_withCut[i]->Clone(Form("hRatio_Data_eMax_cent%d",i));
    hRatio_Data_eMax[i]->Divide(hData_eMax[i]);

    hRatio_Data_chSum[i] = (TH1F*)hData_chSum_withCut[i]->Clone(Form("hRatio_Data_chSum_cent%d",i));
    hRatio_Data_chSum[i]->Divide(hData_chSum[i]);
    hRatio_Data_phSum[i] = (TH1F*)hData_phSum_withCut[i]->Clone(Form("hRatio_Data_phSum_cent%d",i));
    hRatio_Data_phSum[i]->Divide(hData_phSum[i]);
    hRatio_Data_neSum[i] = (TH1F*)hData_neSum_withCut[i]->Clone(Form("hRatio_Data_neSum_cent%d",i));
    hRatio_Data_neSum[i]->Divide(hData_neSum[i]);
    hRatio_Data_muSum[i] = (TH1F*)hData_muSum_withCut[i]->Clone(Form("hRatio_Data_muSum_cent%d",i));
    hRatio_Data_muSum[i]->Divide(hData_muSum[i]);
    hRatio_Data_eSum[i] = (TH1F*)hData_eSum_withCut[i]->Clone(Form("hRatio_Data_eSum_cent%d",i));
    hRatio_Data_eSum[i]->Divide(hData_eSum[i]);
      
  }

  TCanvas * cMC_chMax = new TCanvas("cMC_chMax","",800,600);
  makeMultiPanelCanvas(cMC_chMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_phMax = new TCanvas("cMC_phMax","",800,600);
  makeMultiPanelCanvas(cMC_phMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_neMax = new TCanvas("cMC_neMax","",800,600);
  makeMultiPanelCanvas(cMC_neMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_muMax = new TCanvas("cMC_muMax","",800,600);
  makeMultiPanelCanvas(cMC_muMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_eMax = new TCanvas("cMC_eMax","",800,600);
  makeMultiPanelCanvas(cMC_eMax,3,3,0.0,0.0,0.2,0.15,0.07);

  TCanvas * cMC_chMax_ratio = new TCanvas("cMC_chMax_ratio","",800,600);
  makeMultiPanelCanvas(cMC_chMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_phMax_ratio = new TCanvas("cMC_phMax_ratio","",800,600);
  makeMultiPanelCanvas(cMC_phMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_neMax_ratio = new TCanvas("cMC_neMax_ratio","",800,600);
  makeMultiPanelCanvas(cMC_neMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_muMax_ratio = new TCanvas("cMC_muMax_ratio","",800,600);
  makeMultiPanelCanvas(cMC_muMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_eMax_ratio = new TCanvas("cMC_eMax_ratio","",800,600);
  makeMultiPanelCanvas(cMC_eMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);


  TCanvas * cMC_chSum = new TCanvas("cMC_chSum","",800,600);
  makeMultiPanelCanvas(cMC_chSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_phSum = new TCanvas("cMC_phSum","",800,600);
  makeMultiPanelCanvas(cMC_phSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_neSum = new TCanvas("cMC_neSum","",800,600);
  makeMultiPanelCanvas(cMC_neSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_muSum = new TCanvas("cMC_muSum","",800,600);
  makeMultiPanelCanvas(cMC_muSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_eSum = new TCanvas("cMC_eSum","",800,600);
  makeMultiPanelCanvas(cMC_eSum,3,3,0.0,0.0,0.2,0.15,0.07);

  TCanvas * cMC_chSum_ratio = new TCanvas("cMC_chSum_ratio","",800,600);
  makeMultiPanelCanvas(cMC_chSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_phSum_ratio = new TCanvas("cMC_phSum_ratio","",800,600);
  makeMultiPanelCanvas(cMC_phSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_neSum_ratio = new TCanvas("cMC_neSum_ratio","",800,600);
  makeMultiPanelCanvas(cMC_neSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_muSum_ratio = new TCanvas("cMC_muSum_ratio","",800,600);
  makeMultiPanelCanvas(cMC_muSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cMC_eSum_ratio = new TCanvas("cMC_eSum_ratio","",800,600);
  makeMultiPanelCanvas(cMC_eSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);

  
  TCanvas * cData_chMax = new TCanvas("cData_chMax","",800,600);
  makeMultiPanelCanvas(cData_chMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_phMax = new TCanvas("cData_phMax","",800,600);
  makeMultiPanelCanvas(cData_phMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_neMax = new TCanvas("cData_neMax","",800,600);
  makeMultiPanelCanvas(cData_neMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_muMax = new TCanvas("cData_muMax","",800,600);
  makeMultiPanelCanvas(cData_muMax,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_eMax = new TCanvas("cData_eMax","",800,600);
  makeMultiPanelCanvas(cData_eMax,3,3,0.0,0.0,0.2,0.15,0.07);

  TCanvas * cData_chMax_ratio = new TCanvas("cData_chMax_ratio","",800,600);
  makeMultiPanelCanvas(cData_chMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_phMax_ratio = new TCanvas("cData_phMax_ratio","",800,600);
  makeMultiPanelCanvas(cData_phMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_neMax_ratio = new TCanvas("cData_neMax_ratio","",800,600);
  makeMultiPanelCanvas(cData_neMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_muMax_ratio = new TCanvas("cData_muMax_ratio","",800,600);
  makeMultiPanelCanvas(cData_muMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_eMax_ratio = new TCanvas("cData_eMax_ratio","",800,600);
  makeMultiPanelCanvas(cData_eMax_ratio,3,3,0.0,0.0,0.2,0.15,0.07);


  TCanvas * cData_chSum = new TCanvas("cData_chSum","",800,600);
  makeMultiPanelCanvas(cData_chSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_phSum = new TCanvas("cData_phSum","",800,600);
  makeMultiPanelCanvas(cData_phSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_neSum = new TCanvas("cData_neSum","",800,600);
  makeMultiPanelCanvas(cData_neSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_muSum = new TCanvas("cData_muSum","",800,600);
  makeMultiPanelCanvas(cData_muSum,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_eSum = new TCanvas("cData_eSum","",800,600);
  makeMultiPanelCanvas(cData_eSum,3,3,0.0,0.0,0.2,0.15,0.07);

  TCanvas * cData_chSum_ratio = new TCanvas("cData_chSum_ratio","",800,600);
  makeMultiPanelCanvas(cData_chSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_phSum_ratio = new TCanvas("cData_phSum_ratio","",800,600);
  makeMultiPanelCanvas(cData_phSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_neSum_ratio = new TCanvas("cData_neSum_ratio","",800,600);
  makeMultiPanelCanvas(cData_neSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_muSum_ratio = new TCanvas("cData_muSum_ratio","",800,600);
  makeMultiPanelCanvas(cData_muSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);
  TCanvas * cData_eSum_ratio = new TCanvas("cData_eSum_ratio","",800,600);
  makeMultiPanelCanvas(cData_eSum_ratio,3,3,0.0,0.0,0.2,0.15,0.07);

  
  for(int i = 0;i<nbins_cent;++i){

    cMC_chMax->cd(nbins_cent-i);
    cMC_chMax->cd(nbins_cent-i)->SetLogy();

    hMC_chMax[i]->SetYTitle("counts");
    hMC_chMax[i]->SetXTitle("chMax");
    hMC_chMax[i]->SetMarkerColor(kBlack);
    hMC_chMax[i]->SetMarkerStyle(24);
    hMC_chMax[i]->Draw();

    hMC_chMax_withCut[i]->SetMarkerColor(kRed);
    hMC_chMax_withCut[i]->SetMarkerStyle(24);
    hMC_chMax_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_chMax = myLegend(0.4,0.4,0.6,0.6);
      l_chMax->AddEntry(hMC_chMax[i],"noCut","pl");
      l_chMax->AddEntry(hMC_chMax_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_chMax->SetTextSize(0.04);
      l_chMax->Draw();
    }

    cMC_phMax->cd(nbins_cent-i);
    cMC_phMax->cd(nbins_cent-i)->SetLogy();

    hMC_phMax[i]->SetYTitle("counts");
    hMC_phMax[i]->SetXTitle("chMax");
    hMC_phMax[i]->SetMarkerColor(kBlack);
    hMC_phMax[i]->SetMarkerStyle(24);
    hMC_phMax[i]->Draw();

    hMC_phMax_withCut[i]->SetMarkerColor(kRed);
    hMC_phMax_withCut[i]->SetMarkerStyle(24);
    hMC_phMax_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_phMax = myLegend(0.4,0.4,0.6,0.6);
      l_phMax->AddEntry(hMC_phMax[i],"noCut","pl");
      l_phMax->AddEntry(hMC_phMax_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_phMax->SetTextSize(0.04);
      l_phMax->Draw();
    }

    cMC_neMax->cd(nbins_cent-i);
    cMC_neMax->cd(nbins_cent-i)->SetLogy();

    hMC_neMax[i]->SetYTitle("counts");
    hMC_neMax[i]->SetXTitle("chMax");
    hMC_neMax[i]->SetMarkerColor(kBlack);
    hMC_neMax[i]->SetMarkerStyle(24);
    hMC_neMax[i]->Draw();

    hMC_neMax_withCut[i]->SetMarkerColor(kRed);
    hMC_neMax_withCut[i]->SetMarkerStyle(24);
    hMC_neMax_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_neMax = myLegend(0.4,0.4,0.6,0.6);
      l_neMax->AddEntry(hMC_neMax[i],"noCut","pl");
      l_neMax->AddEntry(hMC_neMax_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_neMax->SetTextSize(0.04);
      l_neMax->Draw();
    }

    cMC_muMax->cd(nbins_cent-i);
    cMC_muMax->cd(nbins_cent-i)->SetLogy();

    hMC_muMax[i]->SetYTitle("counts");
    hMC_muMax[i]->SetXTitle("chMax");
    hMC_muMax[i]->SetMarkerColor(kBlack);
    hMC_muMax[i]->SetMarkerStyle(24);
    hMC_muMax[i]->Draw();

    hMC_muMax_withCut[i]->SetMarkerColor(kRed);
    hMC_muMax_withCut[i]->SetMarkerStyle(24);
    hMC_muMax_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_muMax = myLegend(0.4,0.4,0.6,0.6);
      l_muMax->AddEntry(hMC_muMax[i],"noCut","pl");
      l_muMax->AddEntry(hMC_muMax_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_muMax->SetTextSize(0.04);
      l_muMax->Draw();
    }

    cMC_eMax->cd(nbins_cent-i);
    cMC_eMax->cd(nbins_cent-i)->SetLogy();

    hMC_eMax[i]->SetYTitle("counts");
    hMC_eMax[i]->SetXTitle("chMax");
    hMC_eMax[i]->SetMarkerColor(kBlack);
    hMC_eMax[i]->SetMarkerStyle(24);
    hMC_eMax[i]->Draw();

    hMC_eMax_withCut[i]->SetMarkerColor(kRed);
    hMC_eMax_withCut[i]->SetMarkerStyle(24);
    hMC_eMax_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_eMax = myLegend(0.4,0.4,0.6,0.6);
      l_eMax->AddEntry(hMC_eMax[i],"noCut","pl");
      l_eMax->AddEntry(hMC_eMax_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_eMax->SetTextSize(0.04);
      l_eMax->Draw();
    }

    cMC_chSum->cd(nbins_cent-i);
    cMC_chSum->cd(nbins_cent-i)->SetLogy();

    hMC_chSum[i]->SetYTitle("counts");
    hMC_chSum[i]->SetXTitle("chSum");
    hMC_chSum[i]->SetMarkerColor(kBlack);
    hMC_chSum[i]->SetMarkerStyle(24);
    hMC_chSum[i]->Draw();

    hMC_chSum_withCut[i]->SetMarkerColor(kRed);
    hMC_chSum_withCut[i]->SetMarkerStyle(24);
    hMC_chSum_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_chSum = myLegend(0.4,0.4,0.6,0.6);
      l_chSum->AddEntry(hMC_chSum[i],"noCut","pl");
      l_chSum->AddEntry(hMC_chSum_withCut[i],"chSum/jtpt >0.02 && eSum/jtpt<0.6","pl");
      l_chSum->SetTextSize(0.04);
      l_chSum->Draw();
    }

    cMC_phSum->cd(nbins_cent-i);
    cMC_phSum->cd(nbins_cent-i)->SetLogy();

    hMC_phSum[i]->SetYTitle("counts");
    hMC_phSum[i]->SetXTitle("chSum");
    hMC_phSum[i]->SetMarkerColor(kBlack);
    hMC_phSum[i]->SetMarkerStyle(24);
    hMC_phSum[i]->Draw();

    hMC_phSum_withCut[i]->SetMarkerColor(kRed);
    hMC_phSum_withCut[i]->SetMarkerStyle(24);
    hMC_phSum_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_phSum = myLegend(0.4,0.4,0.6,0.6);
      l_phSum->AddEntry(hMC_phSum[i],"noCut","pl");
      l_phSum->AddEntry(hMC_phSum_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_phSum->SetTextSize(0.04);
      l_phSum->Draw();
    }

    cMC_neSum->cd(nbins_cent-i);
    cMC_neSum->cd(nbins_cent-i)->SetLogy();

    hMC_neSum[i]->SetYTitle("counts");
    hMC_neSum[i]->SetXTitle("chSum");
    hMC_neSum[i]->SetMarkerColor(kBlack);
    hMC_neSum[i]->SetMarkerStyle(24);
    hMC_neSum[i]->Draw();

    hMC_neSum_withCut[i]->SetMarkerColor(kRed);
    hMC_neSum_withCut[i]->SetMarkerStyle(24);
    hMC_neSum_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_neSum = myLegend(0.4,0.4,0.6,0.6);
      l_neSum->AddEntry(hMC_neSum[i],"noCut","pl");
      l_neSum->AddEntry(hMC_neSum_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_neSum->SetTextSize(0.04);
      l_neSum->Draw();
    }

    cMC_muSum->cd(nbins_cent-i);
    cMC_muSum->cd(nbins_cent-i)->SetLogy();

    hMC_muSum[i]->SetYTitle("counts");
    hMC_muSum[i]->SetXTitle("chSum");
    hMC_muSum[i]->SetMarkerColor(kBlack);
    hMC_muSum[i]->SetMarkerStyle(24);
    hMC_muSum[i]->Draw();

    hMC_muSum_withCut[i]->SetMarkerColor(kRed);
    hMC_muSum_withCut[i]->SetMarkerStyle(24);
    hMC_muSum_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_muSum = myLegend(0.4,0.4,0.6,0.6);
      l_muSum->AddEntry(hMC_muSum[i],"noCut","pl");
      l_muSum->AddEntry(hMC_muSum_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_muSum->SetTextSize(0.04);
      l_muSum->Draw();
    }

    cMC_eSum->cd(nbins_cent-i);
    cMC_eSum->cd(nbins_cent-i)->SetLogy();

    hMC_eSum[i]->SetYTitle("counts");
    hMC_eSum[i]->SetXTitle("chSum");
    hMC_eSum[i]->SetMarkerColor(kBlack);
    hMC_eSum[i]->SetMarkerStyle(24);
    hMC_eSum[i]->Draw();

    hMC_eSum_withCut[i]->SetMarkerColor(kRed);
    hMC_eSum_withCut[i]->SetMarkerStyle(24);
    hMC_eSum_withCut[i]->Draw("same");
    
    if(i==nbins_cent-1){
      TLegend * l_eSum = myLegend(0.4,0.4,0.6,0.6);
      l_eSum->AddEntry(hMC_eSum[i],"noCut","pl");
      l_eSum->AddEntry(hMC_eSum_withCut[i],"chMax/jtpt >0.02 && eMax/jtpt<0.6","pl");
      l_eSum->SetTextSize(0.04);
      l_eSum->Draw();
    }

  }
  
  cMC_chMax->SaveAs("","RECREATE");
  cMC_phMax->SaveAs("","RECREATE");
  cMC_neMax->SaveAs("","RECREATE");
  cMC_muMax->SaveAs("","RECREATE");
  cMC_eMax->SaveAs("","RECREATE");

  cMC_chMax_ratio->SaveAs("","RECREATE");
  cMC_phMax_ratio->SaveAs("","RECREATE");
  cMC_neMax_ratio->SaveAs("","RECREATE");
  cMC_muMax_ratio->SaveAs("","RECREATE");
  cMC_eMax_ratio->SaveAs("","RECREATE");

  cMC_chSum->SaveAs("","RECREATE");
  cMC_phSum->SaveAs("","RECREATE");
  cMC_neSum->SaveAs("","RECREATE");
  cMC_muSum->SaveAs("","RECREATE");
  cMC_eSum->SaveAs("","RECREATE");

  cMC_chSum_ratio->SaveAs("","RECREATE");
  cMC_phSum_ratio->SaveAs("","RECREATE");
  cMC_neSum_ratio->SaveAs("","RECREATE");
  cMC_muSum_ratio->SaveAs("","RECREATE");
  cMC_eSum_ratio->SaveAs("","RECREATE");

  cData_chMax->SaveAs("","RECREATE");
  cData_phMax->SaveAs("","RECREATE");
  cData_neMax->SaveAs("","RECREATE");
  cData_muMax->SaveAs("","RECREATE");
  cData_eMax->SaveAs("","RECREATE");

  cData_chMax_ratio->SaveAs("","RECREATE");
  cData_phMax_ratio->SaveAs("","RECREATE");
  cData_neMax_ratio->SaveAs("","RECREATE");
  cData_muMax_ratio->SaveAs("","RECREATE");
  cData_eMax_ratio->SaveAs("","RECREATE");

  cData_chSum->SaveAs("","RECREATE");
  cData_phSum->SaveAs("","RECREATE");
  cData_neSum->SaveAs("","RECREATE");
  cData_muSum->SaveAs("","RECREATE");
  cData_eSum->SaveAs("","RECREATE");

  cData_chSum_ratio->SaveAs("","RECREATE");
  cData_phSum_ratio->SaveAs("","RECREATE");
  cData_neSum_ratio->SaveAs("","RECREATE");
  cData_muSum_ratio->SaveAs("","RECREATE");
  cData_eSum_ratio->SaveAs("","RECREATE");
  
  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
