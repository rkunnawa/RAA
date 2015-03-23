// Raghav Kunnawalkam Elayavalli
// Rutgers University
// Feb 25th 2015
// please contact raghav.k.e at CERN dot CH 

//
// Macro to plot the ALICE plot showing the fragmentation bias 
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

static const int trigValue = 4;
static const char trigNamePbPb [trigValue][256] = {"HLT55","HLT65","HLT80","Combined"};
static const char trigNamePP [trigValue][256] = {"HLT40","HLT60","HLT80","Combined"};
static const int nbins_eta = 1;

static const int no_radius = 1;
static const int list_radius[no_radius] = {2};


//put the radius loop here: 

using namespace std;

void RAA_plot_fragbiascheck(){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};

  // get the input files 
  TFile * fPbPbData = TFile::Open("/Users/raghavke/WORK/RAA/Output/PbPb_fragBiasCheck_akPuPF_20150222.root");
  TFile * fPbPbMC   = TFile::Open("/Users/raghavke/WORK/RAA/Output/PbPb_mc_fragbiascheck_akPuPF_20150223.root");
  TFile * fPPData   = TFile::Open("/Users/raghavke/WORK/RAA/Output/pp_data_fragbiascheck_akPF_20150223.root");
  TFile * fPPMC;

  if(list_radius[0]==2) fPPMC = TFile::Open("/Users/raghavke/WORK/RAA/Output/pp_mc_fragbiascheck_ak2PF_20150226.root");
  if(list_radius[0]==3) fPPMC = TFile::Open("/Users/raghavke/WORK/RAA/Output/pp_mc_fragBiasCheck_ak3PF_20150224.root");

  TH1F *hJets_PbPb_Data_noTrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_Data_3TrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_Data_5TrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_Data_7TrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_Data_10TrkpTCut[no_radius][nbins_cent][trigValue];

  TH1F *hJets_PbPb_MC_noTrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_MC_3TrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_MC_5TrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_MC_7TrkpTCut[no_radius][nbins_cent][trigValue];
  TH1F *hJets_PbPb_MC_10TrkpTCut[no_radius][nbins_cent][trigValue];

  TH1F *hJets_pp_Data_noTrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_Data_noTrkpTCut_ratiowith5[no_radius][trigValue];
  TH1F *hJets_pp_Data_3TrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_Data_5TrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_Data_7TrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_Data_10TrkpTCut[no_radius][trigValue];

  TH1F *hJets_pp_MC_noTrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_MC_noTrkpTCut_ratiowith5[no_radius][trigValue];
  TH1F *hJets_pp_MC_3TrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_MC_5TrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_MC_7TrkpTCut[no_radius][trigValue];
  TH1F *hJets_pp_MC_10TrkpTCut[no_radius][trigValue];

  for(int k = 0;k<no_radius;k++){
    cout<<"radius = "<<k<<endl;
    for(int t = 0;t<trigValue;t++){
      for(int i = 0;i<nbins_cent;i++){

	hJets_PbPb_Data_noTrkpTCut[k][t][i] = (TH1F*)fPbPbData->Get(Form("hJets_noTrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_Data_3TrkpTCut[k][t][i] = (TH1F*)fPbPbData->Get(Form("hJets_3TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_Data_5TrkpTCut[k][t][i] = (TH1F*)fPbPbData->Get(Form("hJets_5TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_Data_7TrkpTCut[k][t][i] = (TH1F*)fPbPbData->Get(Form("hJets_7TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_Data_10TrkpTCut[k][t][i] = (TH1F*)fPbPbData->Get(Form("hJets_10TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));

	hJets_PbPb_Data_noTrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_Data_3TrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_Data_5TrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_Data_7TrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_Data_10TrkpTCut[k][t][i]->Print("base");
		
	hJets_PbPb_Data_noTrkpTCut[k][t][i]->Divide(hJets_PbPb_Data_5TrkpTCut[k][t][i]);
	hJets_PbPb_Data_3TrkpTCut[k][t][i]->Divide(hJets_PbPb_Data_5TrkpTCut[k][t][i]);
	hJets_PbPb_Data_7TrkpTCut[k][t][i]->Divide(hJets_PbPb_Data_5TrkpTCut[k][t][i]);
	hJets_PbPb_Data_10TrkpTCut[k][t][i]->Divide(hJets_PbPb_Data_5TrkpTCut[k][t][i]);

	hJets_PbPb_Data_noTrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_Data_noTrkpTCut[k][t][i]->Scale(1./10);
	hJets_PbPb_Data_3TrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_Data_3TrkpTCut[k][t][i]->Scale(1./10);
	hJets_PbPb_Data_7TrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_Data_7TrkpTCut[k][t][i]->Scale(1./10);
	hJets_PbPb_Data_10TrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_Data_10TrkpTCut[k][t][i]->Scale(1./10);      

	hJets_PbPb_MC_noTrkpTCut[k][t][i] = (TH1F*)fPbPbMC->Get(Form("hJets_noTrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_MC_3TrkpTCut[k][t][i] = (TH1F*)fPbPbMC->Get(Form("hJets_3TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_MC_5TrkpTCut[k][t][i] = (TH1F*)fPbPbMC->Get(Form("hJets_5TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_MC_7TrkpTCut[k][t][i] = (TH1F*)fPbPbMC->Get(Form("hJets_7TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	hJets_PbPb_MC_10TrkpTCut[k][t][i] = (TH1F*)fPbPbMC->Get(Form("hJets_10TrkpTCut_%s_R%d_n20_eta_p20_cent%d",trigNamePbPb[t],list_radius[k],i));
	
	hJets_PbPb_MC_noTrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_MC_3TrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_MC_5TrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_MC_7TrkpTCut[k][t][i]->Print("base");
	hJets_PbPb_MC_10TrkpTCut[k][t][i]->Print("base");
		
	hJets_PbPb_MC_noTrkpTCut[k][t][i]->Divide(hJets_PbPb_MC_5TrkpTCut[k][t][i]);
	hJets_PbPb_MC_3TrkpTCut[k][t][i]->Divide(hJets_PbPb_MC_5TrkpTCut[k][t][i]);
	hJets_PbPb_MC_7TrkpTCut[k][t][i]->Divide(hJets_PbPb_MC_5TrkpTCut[k][t][i]);
	hJets_PbPb_MC_10TrkpTCut[k][t][i]->Divide(hJets_PbPb_MC_5TrkpTCut[k][t][i]);
	
	hJets_PbPb_MC_noTrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_MC_noTrkpTCut[k][t][i]->Scale(1./10);
	hJets_PbPb_MC_3TrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_MC_3TrkpTCut[k][t][i]->Scale(1./10);
	hJets_PbPb_MC_7TrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_MC_7TrkpTCut[k][t][i]->Scale(1./10);
	hJets_PbPb_MC_10TrkpTCut[k][t][i]->Rebin(10);
	hJets_PbPb_MC_10TrkpTCut[k][t][i]->Scale(1./10);   
	
      }
      
      cout<<"trigger value = "<<trigNamePP[t]<<endl;

      hJets_pp_Data_noTrkpTCut[k][t] = (TH1F*)fPPData->Get(Form("hJets_noTrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_Data_noTrkpTCut_ratiowith5[k][t] = (TH1F*)fPPData->Get(Form("hJets_5TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_Data_3TrkpTCut[k][t] = (TH1F*)fPPData->Get(Form("hJets_3TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_Data_5TrkpTCut[k][t] = (TH1F*)fPPData->Get(Form("hJets_5TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_Data_7TrkpTCut[k][t] = (TH1F*)fPPData->Get(Form("hJets_7TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_Data_10TrkpTCut[k][t] = (TH1F*)fPPData->Get(Form("hJets_10TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_Data_noTrkpTCut[k][t]->Print("base");
      hJets_pp_Data_3TrkpTCut[k][t]->Print("base");
      hJets_pp_Data_5TrkpTCut[k][t]->Print("base");
      hJets_pp_Data_7TrkpTCut[k][t]->Print("base");
      hJets_pp_Data_10TrkpTCut[k][t]->Print("base");

      hJets_pp_Data_noTrkpTCut_ratiowith5[k][t]->Divide(hJets_pp_Data_noTrkpTCut[k][t]);
      hJets_pp_Data_noTrkpTCut[k][t]->Divide(hJets_pp_Data_5TrkpTCut[k][t]);
      hJets_pp_Data_3TrkpTCut[k][t]->Divide(hJets_pp_Data_5TrkpTCut[k][t]);      
      hJets_pp_Data_7TrkpTCut[k][t]->Divide(hJets_pp_Data_5TrkpTCut[k][t]);
      hJets_pp_Data_10TrkpTCut[k][t]->Divide(hJets_pp_Data_5TrkpTCut[k][t]);

      hJets_pp_Data_noTrkpTCut_ratiowith5[k][t]->Rebin(10);
      hJets_pp_Data_noTrkpTCut_ratiowith5[k][t]->Scale(1./10);
      hJets_pp_Data_noTrkpTCut[k][t]->Rebin(10);
      hJets_pp_Data_noTrkpTCut[k][t]->Scale(1./10);
      hJets_pp_Data_3TrkpTCut[k][t]->Rebin(10);
      hJets_pp_Data_3TrkpTCut[k][t]->Scale(1./10);
      hJets_pp_Data_7TrkpTCut[k][t]->Rebin(10);
      hJets_pp_Data_7TrkpTCut[k][t]->Scale(1./10);
      hJets_pp_Data_10TrkpTCut[k][t]->Rebin(10);
      hJets_pp_Data_10TrkpTCut[k][t]->Scale(1./10);

      hJets_pp_MC_noTrkpTCut[k][t] = (TH1F*)fPPMC->Get(Form("hJets_noTrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_MC_noTrkpTCut_ratiowith5[k][t] = (TH1F*)fPPMC->Get(Form("hJets_5TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_MC_3TrkpTCut[k][t] = (TH1F*)fPPMC->Get(Form("hJets_3TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_MC_5TrkpTCut[k][t] = (TH1F*)fPPMC->Get(Form("hJets_5TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_MC_7TrkpTCut[k][t] = (TH1F*)fPPMC->Get(Form("hJets_7TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));
      hJets_pp_MC_10TrkpTCut[k][t] = (TH1F*)fPPMC->Get(Form("hJets_10TrkpTCut_%s_R%d_0_absEta_20",trigNamePP[t],list_radius[k]));

      hJets_pp_MC_noTrkpTCut[k][t]->Print("base");
      hJets_pp_MC_3TrkpTCut[k][t]->Print("base");
      hJets_pp_MC_5TrkpTCut[k][t]->Print("base");
      hJets_pp_MC_7TrkpTCut[k][t]->Print("base");
      hJets_pp_MC_10TrkpTCut[k][t]->Print("base");

      hJets_pp_MC_noTrkpTCut_ratiowith5[k][t]->Divide(hJets_pp_MC_noTrkpTCut[k][t]);
      hJets_pp_MC_noTrkpTCut[k][t]->Divide(hJets_pp_MC_5TrkpTCut[k][t]);
      hJets_pp_MC_3TrkpTCut[k][t]->Divide(hJets_pp_MC_5TrkpTCut[k][t]);
      hJets_pp_MC_7TrkpTCut[k][t]->Divide(hJets_pp_MC_5TrkpTCut[k][t]);
      hJets_pp_MC_10TrkpTCut[k][t]->Divide(hJets_pp_MC_5TrkpTCut[k][t]);

      hJets_pp_MC_noTrkpTCut_ratiowith5[k][t]->Rebin(10);
      hJets_pp_MC_noTrkpTCut_ratiowith5[k][t]->Scale(1./10);
      hJets_pp_MC_noTrkpTCut[k][t]->Rebin(10);
      hJets_pp_MC_noTrkpTCut[k][t]->Scale(1./10);
      hJets_pp_MC_3TrkpTCut[k][t]->Rebin(10);
      hJets_pp_MC_3TrkpTCut[k][t]->Scale(1./10);
      hJets_pp_MC_7TrkpTCut[k][t]->Rebin(10);
      hJets_pp_MC_7TrkpTCut[k][t]->Scale(1./10);
      hJets_pp_MC_10TrkpTCut[k][t]->Rebin(10);
      hJets_pp_MC_10TrkpTCut[k][t]->Scale(1./10);
      
    }// trigger value loop
    
  }// radius loop

  cout<<"Done gettng all the histograms "<<endl;
  
  // make the plots, one plot for each centrality bin and each plot for PbPb will have 4 panels for each of the plots showing data and mc. and one separate plot for pp. 

  
  TCanvas *cPbPb[nbins_cent];
  
  for(int i = 0;i<nbins_cent;i++){

    cPbPb[i] = new TCanvas(Form("cPbPb_cent%d",i),"",800,600);
    makeMultiPanelCanvas(cPbPb[i],2,2,0.0,0.0,0.2,0.15,0.07);

    cPbPb[i]->cd(1);
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->Print("base");
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->SetMarkerStyle(20);
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->SetMarkerColor(kBlack);
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->SetYTitle("#frac{jets candidiate pT > X}{jets candidate pT > 5}");
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->SetXTitle("Jet p_{T} (GeV/c)");
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->SetTitle(" ");
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->SetAxisRange(0,140,"X");
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->SetAxisRange(0,1.5,"Y");
    hJets_PbPb_Data_noTrkpTCut[0][3][i]->Draw();

    hJets_PbPb_MC_noTrkpTCut[0][3][i]->SetMarkerStyle(24);
    hJets_PbPb_MC_noTrkpTCut[0][3][i]->SetMarkerColor(kBlack);
    hJets_PbPb_MC_noTrkpTCut[0][3][i]->Draw("same");
    putCMSPrel();
    drawText(Form("akPu%dPF, |#eta|<2 ",list_radius[0]),0.7,0.2,14);
    drawText("X = 0",0.6,0.8,14);

    cPbPb[i]->cd(2);
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->SetMarkerStyle(20);
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->SetMarkerColor(kBlack);
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->SetYTitle("#frac{#jets candidiate pT > X}{# jets candidate pT > 5}");
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->SetXTitle("Jet p_{T} (GeV/c)");
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->SetTitle(" ");
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->SetAxisRange(0,140,"X");
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->SetAxisRange(0,1.5,"Y");
    hJets_PbPb_Data_3TrkpTCut[0][3][i]->Draw();

    hJets_PbPb_MC_3TrkpTCut[0][3][i]->SetMarkerStyle(24);
    hJets_PbPb_MC_3TrkpTCut[0][3][i]->SetMarkerColor(kBlack);
    hJets_PbPb_MC_3TrkpTCut[0][3][i]->Draw("same");
    drawText("X = 3",0.6,0.8,14);
    drawText(Form("%2.0f-%2.0f%% ",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.4,0.8,16);

    TLegend * PbPb = myLegend(0.3,0.2,0.6,0.5);
    PbPb->AddEntry(hJets_PbPb_Data_3TrkpTCut[0][3][i],"Data","pl");
    PbPb->AddEntry(hJets_PbPb_MC_3TrkpTCut[0][3][i],"MC","pl");
    PbPb->SetTextSize(0.04);
    PbPb->Draw();

    cPbPb[i]->cd(3);
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->SetMarkerStyle(20);
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->SetMarkerColor(kBlack);
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->SetYTitle("#frac{#jets candidiate pT > X}{# jets candidate pT > 5}");
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->SetXTitle("Jet p_{T} (GeV/c)");
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->SetTitle(" ");
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->SetAxisRange(0,140,"X");
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->SetAxisRange(0,1.5,"Y");
    hJets_PbPb_Data_7TrkpTCut[0][3][i]->Draw();

    hJets_PbPb_MC_7TrkpTCut[0][3][i]->SetMarkerStyle(24);
    hJets_PbPb_MC_7TrkpTCut[0][3][i]->SetMarkerColor(kBlack);
    hJets_PbPb_MC_7TrkpTCut[0][3][i]->Draw("same");
    drawText("X = 7",0.6,0.8,14);

    cPbPb[i]->cd(4);
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->SetMarkerStyle(20);
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->SetMarkerColor(kBlack);
  
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->SetYTitle("#frac{#jets candidiate pT > X}{# jets candidate pT > 5}");
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->SetXTitle("Jet p_{T} (GeV/c)");
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->SetTitle(" ");
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->SetAxisRange(0,140,"X");
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->SetAxisRange(0,1.5,"Y");
    hJets_PbPb_Data_10TrkpTCut[0][3][i]->Draw();

    hJets_PbPb_MC_10TrkpTCut[0][3][i]->SetMarkerStyle(24);
    hJets_PbPb_MC_10TrkpTCut[0][3][i]->SetMarkerColor(kBlack);
    hJets_PbPb_MC_10TrkpTCut[0][3][i]->Draw("same");
    drawText("X = 10",0.6,0.8,14);


    cPbPb[i]->SaveAs(Form("PbPb_fragmentation_bias_check_R%d_akPuPF_cent%d.pdf",list_radius[0],i),"RECREATE");
  }

  
  // make pp plots: 
  TCanvas * cPP = new TCanvas("cPP","",800,600);
  makeMultiPanelCanvas(cPP,2,2,0.0,0.0,0.2,0.15,0.07);
  cPP->cd(1);
  hJets_pp_Data_noTrkpTCut[0][3]->Print("base");
  hJets_pp_Data_noTrkpTCut[0][3]->SetMarkerStyle(20);
  hJets_pp_Data_noTrkpTCut[0][3]->SetMarkerColor(kBlack);
  hJets_pp_Data_noTrkpTCut[0][3]->SetYTitle("#frac{jets candidiate pT > X}{jets candidate pT > 5}");
  hJets_pp_Data_noTrkpTCut[0][3]->SetXTitle("Jet p_{T} (GeV/c)");
  hJets_pp_Data_noTrkpTCut[0][3]->SetTitle(" ");
  hJets_pp_Data_noTrkpTCut[0][3]->SetAxisRange(0,140,"X");
  hJets_pp_Data_noTrkpTCut[0][3]->SetAxisRange(0,1.5,"Y");
  hJets_pp_Data_noTrkpTCut[0][3]->Draw();

  hJets_pp_MC_noTrkpTCut[0][3]->SetMarkerStyle(24);
  hJets_pp_MC_noTrkpTCut[0][3]->SetMarkerColor(kBlack);
  hJets_pp_MC_noTrkpTCut[0][3]->Draw("same");
  putCMSPrel();
  drawText(Form("ak%dPF",list_radius[0]),0.7,0.2,14);
  drawText("X = 0",0.6,0.8,14);

  cPP->cd(2);
  hJets_pp_Data_3TrkpTCut[0][3]->SetMarkerStyle(20);
  hJets_pp_Data_3TrkpTCut[0][3]->SetMarkerColor(kBlack);
  hJets_pp_Data_3TrkpTCut[0][3]->SetYTitle("#frac{#jets candidiate pT > X}{# jets candidate pT > 5}");
  hJets_pp_Data_3TrkpTCut[0][3]->SetXTitle("Jet p_{T} (GeV/c)");
  hJets_pp_Data_3TrkpTCut[0][3]->SetTitle(" ");
  hJets_pp_Data_3TrkpTCut[0][3]->SetAxisRange(0,140,"X");
  hJets_pp_Data_3TrkpTCut[0][3]->SetAxisRange(0,1.5,"Y");
  hJets_pp_Data_3TrkpTCut[0][3]->Draw();

  hJets_pp_MC_3TrkpTCut[0][3]->SetMarkerStyle(24);
  hJets_pp_MC_3TrkpTCut[0][3]->SetMarkerColor(kBlack);
  hJets_pp_MC_3TrkpTCut[0][3]->Draw("same");
  drawText("X = 3",0.6,0.8,14);

  TLegend * pp = myLegend(0.3,0.2,0.6,0.5);
  pp->AddEntry(hJets_pp_Data_3TrkpTCut[0][3],"Data","pl");
  pp->AddEntry(hJets_pp_MC_3TrkpTCut[0][3],"MC","pl");
  pp->SetTextSize(0.04);
  pp->Draw();

  cPP->cd(3);
  hJets_pp_Data_7TrkpTCut[0][3]->SetMarkerStyle(20);
  hJets_pp_Data_7TrkpTCut[0][3]->SetMarkerColor(kBlack);
  hJets_pp_Data_7TrkpTCut[0][3]->SetYTitle("#frac{#jets candidiate pT > X}{# jets candidate pT > 5}");
  hJets_pp_Data_7TrkpTCut[0][3]->SetXTitle("Jet p_{T} (GeV/c)");
  hJets_pp_Data_7TrkpTCut[0][3]->SetTitle(" ");
  hJets_pp_Data_7TrkpTCut[0][3]->SetAxisRange(0,140,"X");
  hJets_pp_Data_7TrkpTCut[0][3]->SetAxisRange(0,1.5,"Y");
  hJets_pp_Data_7TrkpTCut[0][3]->Draw();

  hJets_pp_MC_7TrkpTCut[0][3]->SetMarkerStyle(24);
  hJets_pp_MC_7TrkpTCut[0][3]->SetMarkerColor(kBlack);
  hJets_pp_MC_7TrkpTCut[0][3]->Draw("same");
  drawText("X = 7",0.6,0.8,14);

  cPP->cd(4);
  hJets_pp_Data_10TrkpTCut[0][3]->SetMarkerStyle(20);
  hJets_pp_Data_10TrkpTCut[0][3]->SetMarkerColor(kBlack);
  
  hJets_pp_Data_10TrkpTCut[0][3]->SetYTitle("#frac{#jets candidiate pT > X}{# jets candidate pT > 5}");
  hJets_pp_Data_10TrkpTCut[0][3]->SetXTitle("Jet p_{T} (GeV/c)");
  hJets_pp_Data_10TrkpTCut[0][3]->SetTitle(" ");
  hJets_pp_Data_10TrkpTCut[0][3]->SetAxisRange(0,140,"X");
  hJets_pp_Data_10TrkpTCut[0][3]->SetAxisRange(0,1.5,"Y");
  hJets_pp_Data_10TrkpTCut[0][3]->Draw();

  hJets_pp_MC_10TrkpTCut[0][3]->SetMarkerStyle(24);
  hJets_pp_MC_10TrkpTCut[0][3]->SetMarkerColor(kBlack);
  hJets_pp_MC_10TrkpTCut[0][3]->Draw("same");
  drawText("X = 10",0.6,0.8,14);

  cPP->SaveAs(Form("pp_fragmentation_bias_check_R%d_akPF.pdf",list_radius[0]),"RECREATE");
  

  TCanvas *cPP_1panel = new TCanvas("cPP_1panel","",800,600);
 
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->SetMarkerStyle(20);
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->SetMarkerColor(kBlack);
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->SetYTitle("#frac{jets candidate pT > 5}{jets candidiate pT > 0}");
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->SetXTitle("Jet p_{T} (GeV/c)");
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->SetTitle(" ");
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->SetAxisRange(0,140,"X");
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->SetAxisRange(0,1.5,"Y");
  hJets_pp_Data_noTrkpTCut_ratiowith5[0][3]->Draw();

  hJets_pp_MC_noTrkpTCut_ratiowith5[0][3]->SetMarkerStyle(24);
  hJets_pp_MC_noTrkpTCut_ratiowith5[0][3]->SetMarkerColor(kBlack);
  hJets_pp_MC_noTrkpTCut_ratiowith5[0][3]->Draw("same");

  putCMSPrel();
  drawText(Form("ak%dPF",list_radius[0]),0.2,0.3,16);
  cPP_1panel->SaveAs(Form("pp_fragmentation_bias_check_inverseratio_R%d_akPF.pdf",list_radius[0]),"RECREATE");

  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;


}
