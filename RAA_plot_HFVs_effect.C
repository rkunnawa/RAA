// Raghav Kunnawalkam Elayavalli
// Dec 7th 2014
// Rutgers
// raghav.k.eat CERN dot CH

//
// plot the effect of the HF/Vs polynomial divergece events on the Jet Spectra. 
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


static const int nbins_eta = 15;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0}, {-2.0,+2.0}, {-3.0,+3.0},
  {-3.0,-2.5}, {-2.5,-2.0}, {-2.0,-1.5}, 
  {-1.5,-1.0}, {-1.0,-0.5}, {-0.5,0}, {0,+0.5}, 
  {+0.5,+1.0}, {+1.0,+1.5}, {+1.5,+2.0}, 
  {+2.0,+2.5}, {+2.5,+3.0}
};

static const double delta_eta[nbins_eta] = {
  2.0, 4.0, 6.0, 
  0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 
  0.5, 0.5
};

static const char etaWidth [nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20","n30_eta_p30",
  "n30_eta_n25","n25_eta_n20","n20_eta_n15",
  "n15_eta_n10","n10_eta_n05","n05_eta_0","0_eta_p05",
  "p05_eta_p10","p10_eta_p15","p15_eta_p20",
  "p20_eta_p25","p25_eta_p30"
};


//these are the only radii we are interested for the RAA analysis: 2,3,4,5
static const int no_radius = 7; 
static const int list_radius[no_radius] = {1,2,3,4,5,6,7};


using namespace std;


void RAA_plot_HFVs_effect(){
  
  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  
  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  TFile *fin = TFile::Open("/Users/raghavke/WORK/RAA/Output/PbPb_HF_divergence_events_spectra_VsPF_20141207.root");
  
  TH1F *hpbpb_Jet80[2][2][nbins_cent+1];
  TH1F *hpbpb_Jet65[2][2][nbins_cent+1];
  TH1F *hpbpb_Jet55[2][2][nbins_cent+1];
  TH1F *hpbpb_JetComb[2][2][nbins_cent+1];  

  TH1F *hpbpb_FullJet80[nbins_cent+1];
  TH1F *hpbpb_FullJet65[nbins_cent+1];
  TH1F *hpbpb_FullJet55[nbins_cent+1];
  TH1F *hpbpb_FullJetComb[nbins_cent+1];

  //Get the ratio histograms w.r.t the full spectra. 
  TH1F *hpbpb_Jet80_Ratio[2][2][nbins_cent+1];
  TH1F *hpbpb_Jet65_Ratio[2][2][nbins_cent+1];
  TH1F *hpbpb_Jet55_Ratio[2][2][nbins_cent+1];

  for(int i = 0;i<nbins_cent+1;i++){

    for(int a = 0;a<2;a++){

      for(int b = 0;b<2;b++){

	hpbpb_Jet80[b][a][i] = (TH1F*)fin->Get(Form("hpbpb_Jet80_isDiverge_%d_isTightCut_%d_cent%d",b,a,i));
	hpbpb_Jet80[b][a][i] = (TH1F*)hpbpb_Jet80[b][a][i]->Rebin(nbins_pt,Form("hpbpb_Jet80_isDiverge_%d_isTightCut_%d_cent%d",b,a,i),boundaries_pt);
	divideBinWidth(hpbpb_Jet80[b][a][i]);
	hpbpb_Jet65[b][a][i] = (TH1F*)fin->Get(Form("hpbpb_Jet65_isDiverge_%d_isTightCut_%d_cent%d",b,a,i));
	hpbpb_Jet65[b][a][i] = (TH1F*)hpbpb_Jet65[b][a][i]->Rebin(nbins_pt,Form("hpbpb_Jet65_isDiverge_%d_isTightCut_%d_cent%d",b,a,i),boundaries_pt);
	divideBinWidth(hpbpb_Jet65[b][a][i]);
	hpbpb_Jet55[b][a][i] = (TH1F*)fin->Get(Form("hpbpb_Jet55_isDiverge_%d_isTightCut_%d_cent%d",b,a,i));
	hpbpb_Jet55[b][a][i] = (TH1F*)hpbpb_Jet55[b][a][i]->Rebin(nbins_pt,Form("hpbpb_Jet55_isDiverge_%d_isTightCut_%d_cent%d",b,a,i),boundaries_pt);
	divideBinWidth(hpbpb_Jet55[b][a][i]);

      }

    }

    hpbpb_FullJet80[i] = (TH1F*)fin->Get(Form("hpbpb_Jet80_cent%d",i));
    hpbpb_FullJet80[i] = (TH1F*)hpbpb_FullJet80[i]->Rebin(nbins_pt,Form("hpbpb_Jet80_cent%d",i),boundaries_pt);
    divideBinWidth(hpbpb_FullJet80[i]);
    hpbpb_FullJet65[i] = (TH1F*)fin->Get(Form("hpbpb_Jet65_cent%d",i));
    hpbpb_FullJet65[i] = (TH1F*)hpbpb_FullJet65[i]->Rebin(nbins_pt,Form("hpbpb_Jet65_cent%d",i),boundaries_pt);
    divideBinWidth(hpbpb_FullJet65[i]);
    hpbpb_FullJet55[i] = (TH1F*)fin->Get(Form("hpbpb_Jet55_cent%d",i));
    hpbpb_FullJet55[i] = (TH1F*)hpbpb_FullJet55[i]->Rebin(nbins_pt,Form("hpbpb_Jet55_cent%d",i),boundaries_pt);
    divideBinWidth(hpbpb_FullJet55[i]);

    for(int a = 0;a<2;a++){

      for(int b = 0;b<2;b++){

	hpbpb_Jet80_Ratio[b][a][i] = (TH1F*)hpbpb_Jet80[b][a][i]->Clone(Form("hpbpb_Jet80_isDiverge_%d_isTightCut_%d_cent%d_Ratio_with_fullSpectra_Jet80",b,a,i));
	hpbpb_Jet80_Ratio[b][a][i]->Divide(hpbpb_FullJet80[i]);
	hpbpb_Jet65_Ratio[b][a][i] = (TH1F*)hpbpb_Jet65[b][a][i]->Clone(Form("hpbpb_Jet65_isDiverge_%d_isTightCut_%d_cent%d_Ratio_with_fullSpectra_Jet65",b,a,i));
	hpbpb_Jet65_Ratio[b][a][i]->Divide(hpbpb_FullJet65[i]);
	hpbpb_Jet55_Ratio[b][a][i] = (TH1F*)hpbpb_Jet55[b][a][i]->Clone(Form("hpbpb_Jet80_isDiverge_%d_isTightCut_%d_cent%d_Ratio_with_fullSpectra_Jet55",b,a,i));
	hpbpb_Jet55_Ratio[b][a][i]->Divide(hpbpb_FullJet55[i]);
      }

    }
    
  }// cent bin
  
  // Lets start making the plots. These plots will all be a 3x3 panel plot for the different centrality bins. 
  // plots 1,2 and 3 - show the spectra for isdiverge 1 and with and without tight cut for each Jet80 and Jet65 and Jet55 in each centrality  
  // plots 4,5 and 6 - similarly show the ratio histogram for each Jet80 and Jet65 and Jet55 in each centrailty. 

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // plot 1 - is diverge with and without tight cuts plus the whole events list. 

  TCanvas *cJet80_Spectra = new TCanvas("cJet80_Spectra","",800,600);
  makeMultiPanelCanvas(cJet80_Spectra,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0;i<nbins_cent;i++){

    cJet80_Spectra->cd(nbins_cent-i);
    cJet80_Spectra->cd(nbins_cent-i)->SetLogy();
    
    hpbpb_FullJet80[i]->SetMarkerStyle(24);
    hpbpb_FullJet80[i]->SetMarkerColor(kBlack);
    hpbpb_FullJet80[i]->SetTitle(" ");
    hpbpb_FullJet80[i]->SetXTitle("Jet p_{T} (GeV/c)");
    hpbpb_FullJet80[i]->SetYTitle("#frac{dN}{dp_{T}}");
    hpbpb_FullJet80[i]->SetAxisRange(0,500,"X");
    hpbpb_FullJet80[i]->Draw();

    hpbpb_Jet80[1][0][i]->SetMarkerStyle(25);
    hpbpb_Jet80[1][0][i]->SetMarkerColor(kRed);
    hpbpb_Jet80[1][0][i]->Draw("same");

    hpbpb_Jet80[1][1][i]->SetMarkerStyle(33);
    hpbpb_Jet80[1][1][i]->SetMarkerColor(kBlue);
    hpbpb_Jet80[1][1][i]->Draw("same");
    
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);

  }

  cJet80_Spectra->cd(1);
  putCMSPrel();

  TLegend *Jet80_spectra = myLegend(0.5,0.55,0.9,0.8);
  Jet80_spectra->AddEntry(hpbpb_FullJet80[0],"Full Spectra","pl");
  Jet80_spectra->AddEntry(hpbpb_Jet80[1][0][0],"Diverging with Loose Cut","pl");
  Jet80_spectra->AddEntry(hpbpb_Jet80[1][1][0],"Diverging with Tight Cut","pl");
  Jet80_spectra->SetTextSize(0.04);
  Jet80_spectra->Draw();
  
  drawText("Jet80 Trigger",0.5,0.4,14);
  
  cJet80_Spectra->cd(2);
  drawText("pCES,HBHE,|vz|<15",0.2,0.8,14);
  
  cJet80_Spectra->cd(3);
  drawText("supernova diagonal cut",0.2,0.8,14);
  
  cJet80_Spectra->SaveAs(Form("/Users/raghavke/WORK/RAA/Plots/HFVS_RAA_Effect_Jet80_spectra_%d.pdf",date.GetDate()),"RECREATE");
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // plot 2 - is diverge with and without tight cuts plus the whole events list - Jet65

  TCanvas *cJet65_Spectra = new TCanvas("cJet65_Spectra","",800,600);
  makeMultiPanelCanvas(cJet65_Spectra,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0;i<nbins_cent;i++){

    cJet65_Spectra->cd(nbins_cent-i);
    cJet65_Spectra->cd(nbins_cent-i)->SetLogy();
    
    hpbpb_FullJet65[i]->SetMarkerStyle(24);
    hpbpb_FullJet65[i]->SetMarkerColor(kBlack);
    hpbpb_FullJet65[i]->SetTitle(" ");
    hpbpb_FullJet65[i]->SetXTitle("Jet p_{T} (GeV/c)");
    hpbpb_FullJet65[i]->SetYTitle("#frac{dN}{dp_{T}}");
    hpbpb_FullJet65[i]->SetAxisRange(0,500,"X");
    hpbpb_FullJet65[i]->Draw();
    
    hpbpb_Jet65[1][0][i]->SetMarkerStyle(25);
    hpbpb_Jet65[1][0][i]->SetMarkerColor(kRed);
    hpbpb_Jet65[1][0][i]->Draw("same");
    
    hpbpb_Jet65[1][1][i]->SetMarkerStyle(33);
    hpbpb_Jet65[1][1][i]->SetMarkerColor(kBlue);
    hpbpb_Jet65[1][1][i]->Draw("same");
    
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);
    
  }

  cJet65_Spectra->cd(1);
  putCMSPrel();

  TLegend *Jet65_spectra = myLegend(0.5,0.55,0.9,0.8);
  Jet65_spectra->AddEntry(hpbpb_FullJet65[0],"Full Spectra","pl");
  Jet65_spectra->AddEntry(hpbpb_Jet65[1][0][0],"Diverging with Loose Cut","pl");
  Jet65_spectra->AddEntry(hpbpb_Jet65[1][1][0],"Diverging with Tight Cut","pl");
  Jet65_spectra->SetTextSize(0.04);
  Jet65_spectra->Draw();
  
  drawText("Jet65 Trigger",0.5,0.4,14);
  
  cJet65_Spectra->cd(2);
  drawText("pCES,HBHE,|vz|<15",0.2,0.8,14);
  
  cJet65_Spectra->cd(3);
  drawText("supernova diagonal cut",0.2,0.8,14);
  
  cJet65_Spectra->SaveAs(Form("/Users/raghavke/WORK/RAA/Plots/HFVS_RAA_Effect_Jet65_spectra_%d.pdf",date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // plot 3 - is diverge with and without tight cuts plus the whole events list - Jet65

  TCanvas *cJet55_Spectra = new TCanvas("cJet55_Spectra","",800,600);
  makeMultiPanelCanvas(cJet55_Spectra,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0;i<nbins_cent;i++){

    cJet55_Spectra->cd(nbins_cent-i);
    cJet55_Spectra->cd(nbins_cent-i)->SetLogy();
    
    hpbpb_FullJet55[i]->SetMarkerStyle(24);
    hpbpb_FullJet55[i]->SetMarkerColor(kBlack);
    hpbpb_FullJet55[i]->SetTitle(" ");
    hpbpb_FullJet55[i]->SetXTitle("Jet p_{T} (GeV/c)");
    hpbpb_FullJet55[i]->SetYTitle("#frac{dN}{dp_{T}}");
    hpbpb_FullJet55[i]->SetAxisRange(0,500,"X");
    hpbpb_FullJet55[i]->Draw();
    
    hpbpb_Jet55[1][0][i]->SetMarkerStyle(25);
    hpbpb_Jet55[1][0][i]->SetMarkerColor(kRed);
    hpbpb_Jet55[1][0][i]->Draw("same");
    
    hpbpb_Jet55[1][1][i]->SetMarkerStyle(33);
    hpbpb_Jet55[1][1][i]->SetMarkerColor(kBlue);
    hpbpb_Jet55[1][1][i]->Draw("same");
    
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);
    
  }

  cJet55_Spectra->cd(1);
  putCMSPrel();

  TLegend *Jet55_spectra = myLegend(0.5,0.55,0.9,0.8);
  Jet55_spectra->AddEntry(hpbpb_FullJet55[0],"Full Spectra","pl");
  Jet55_spectra->AddEntry(hpbpb_Jet55[1][0][0],"Diverging with Loose Cut","pl");
  Jet55_spectra->AddEntry(hpbpb_Jet55[1][1][0],"Diverging with Tight Cut","pl");
  Jet55_spectra->SetTextSize(0.04);
  Jet55_spectra->Draw();
  
  drawText("Jet55 Trigger",0.5,0.4,14);
  
  cJet55_Spectra->cd(2);
  drawText("pCES,HBHE,|vz|<15",0.2,0.8,14);
  
  cJet55_Spectra->cd(3);
  drawText("supernova diagonal cut",0.2,0.8,14);
  
  cJet55_Spectra->SaveAs(Form("/Users/raghavke/WORK/RAA/Plots/HFVS_RAA_Effect_Jet55_spectra_%d.pdf",date.GetDate()),"RECREATE");

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // plot 4 - Ratio plot: is diverge with and without tight cuts plus the whole events list - Jet65

  TCanvas *cJet80_Ratio = new TCanvas("cJet80_Ratio","",800,600);
  makeMultiPanelCanvas(cJet80_Ratio,3,2,0.0,0.0,0.2,0.15,0.07);

  TLine *line80 = new TLine(0,1,500,1);
  line80->SetLineStyle(2);
  line80->SetLineWidth(2);

  for(int i = 0;i<nbins_cent;i++){

    cJet80_Ratio->cd(nbins_cent-i);
    //cJet80_Ratio->cd(nbins_cent-i)->SetLogy();
    
    hpbpb_Jet80_Ratio[1][0][i]->SetMarkerStyle(25);
    hpbpb_Jet80_Ratio[1][0][i]->SetMarkerColor(kRed);
    hpbpb_Jet80_Ratio[1][0][i]->SetXTitle("Jet p_{T} (GeV/c)");
    hpbpb_Jet80_Ratio[1][0][i]->SetYTitle("Ratio with full spectra");
    hpbpb_Jet80_Ratio[1][0][i]->SetTitle(" ");
    hpbpb_Jet80_Ratio[1][0][i]->SetAxisRange(0,500,"X");
    hpbpb_Jet80_Ratio[1][0][i]->SetAxisRange(0,1,"Y");
    hpbpb_Jet80_Ratio[1][0][i]->Draw("same");
    
    hpbpb_Jet80_Ratio[1][1][i]->SetMarkerStyle(33);
    hpbpb_Jet80_Ratio[1][1][i]->SetMarkerColor(kBlue);
    hpbpb_Jet80_Ratio[1][1][i]->Draw("same");
 
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);
    
    line80->Draw();
  }

  cJet80_Ratio->cd(1);
  putCMSPrel();

  TLegend *Jet80_ratio = myLegend(0.5,0.55,0.9,0.8);
  Jet80_ratio->AddEntry(hpbpb_Jet80[1][0][0],"Diverging with Loose Cut","pl");
  Jet80_ratio->AddEntry(hpbpb_Jet80[1][1][0],"Diverging with Tight Cut","pl");
  Jet80_ratio->SetTextSize(0.04);
  Jet80_ratio->Draw();
  
  drawText("Jet80 Trigger",0.5,0.4,14);
  
  cJet80_Ratio->cd(2);
  drawText("pCES,HBHE,|vz|<15",0.2,0.8,14);
  drawText("Ratio plots with full spectra",0.15,0.7,14);
  
  cJet80_Ratio->cd(3);
  drawText("supernova diagonal cut",0.2,0.8,14);
  
  cJet80_Ratio->SaveAs(Form("/Users/raghavke/WORK/RAA/Plots/HFVS_RAA_Effect_Jet80_ratio_%d.pdf",date.GetDate()),"RECREATE");
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // plot 5 - Ratio plot: is diverge with and without tight cuts plus the whole events list - Jet65

  TCanvas *cJet65_Ratio = new TCanvas("cJet65_Ratio","",800,600);
  makeMultiPanelCanvas(cJet65_Ratio,3,2,0.0,0.0,0.2,0.15,0.07);

  TLine *line65 = new TLine(0,1,500,1);
  line65->SetLineStyle(2);
  line65->SetLineWidth(2);

  for(int i = 0;i<nbins_cent;i++){

    cJet65_Ratio->cd(nbins_cent-i);
    //cJet65_Ratio->cd(nbins_cent-i)->SetLogy();
    
    hpbpb_Jet65_Ratio[1][0][i]->SetMarkerStyle(25);
    hpbpb_Jet65_Ratio[1][0][i]->SetMarkerColor(kRed);
    hpbpb_Jet65_Ratio[1][0][i]->SetXTitle("Jet p_{T} (GeV/c)");
    hpbpb_Jet65_Ratio[1][0][i]->SetYTitle("Ratio with full spectra");
    hpbpb_Jet65_Ratio[1][0][i]->SetTitle(" ");
    hpbpb_Jet65_Ratio[1][0][i]->SetAxisRange(0,500,"X");
    hpbpb_Jet65_Ratio[1][0][i]->SetAxisRange(0,1,"Y");
    hpbpb_Jet65_Ratio[1][0][i]->Draw("same");
    
    hpbpb_Jet65_Ratio[1][1][i]->SetMarkerStyle(33);
    hpbpb_Jet65_Ratio[1][1][i]->SetMarkerColor(kBlue);
    hpbpb_Jet65_Ratio[1][1][i]->Draw("same");
 
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);
    
    line65->Draw();
  }

  cJet65_Ratio->cd(1);
  putCMSPrel();

  TLegend *Jet65_ratio = myLegend(0.5,0.55,0.9,0.8);
  Jet65_ratio->AddEntry(hpbpb_Jet65[1][0][0],"Diverging with Loose Cut","pl");
  Jet65_ratio->AddEntry(hpbpb_Jet65[1][1][0],"Diverging with Tight Cut","pl");
  Jet65_ratio->SetTextSize(0.04);
  Jet65_ratio->Draw();
  
  drawText("Jet65 Trigger",0.5,0.4,14);
  
  cJet65_Ratio->cd(2);
  drawText("pCES,HBHE,|vz|<15",0.2,0.8,14);
  drawText("Ratio plots with full spectra",0.15,0.7,14);
  
  cJet65_Ratio->cd(3);
  drawText("supernova diagonal cut",0.2,0.8,14);
  
  cJet65_Ratio->SaveAs(Form("/Users/raghavke/WORK/RAA/Plots/HFVS_RAA_Effect_Jet65_ratio_%d.pdf",date.GetDate()),"RECREATE");


  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // plot 6 - Ratio plot: is diverge with and without tight cuts plus the whole events list - Jet55

  TCanvas *cJet55_Ratio = new TCanvas("cJet55_Ratio","",800,600);
  makeMultiPanelCanvas(cJet55_Ratio,3,2,0.0,0.0,0.2,0.15,0.07);

  TLine *line55 = new TLine(0,1,500,1);
  line55->SetLineStyle(2);
  line55->SetLineWidth(2);

  for(int i = 0;i<nbins_cent;i++){

    cJet55_Ratio->cd(nbins_cent-i);
    //cJet55_Ratio->cd(nbins_cent-i)->SetLogy();
    
    hpbpb_Jet55_Ratio[1][0][i]->SetMarkerStyle(25);
    hpbpb_Jet55_Ratio[1][0][i]->SetMarkerColor(kRed);
    hpbpb_Jet55_Ratio[1][0][i]->SetXTitle("Jet p_{T} (GeV/c)");
    hpbpb_Jet55_Ratio[1][0][i]->SetYTitle("Ratio with full spectra");
    hpbpb_Jet55_Ratio[1][0][i]->SetTitle(" ");
    hpbpb_Jet55_Ratio[1][0][i]->SetAxisRange(0,500,"X");
    hpbpb_Jet55_Ratio[1][0][i]->SetAxisRange(0,1.1,"Y");
    hpbpb_Jet55_Ratio[1][0][i]->Draw("same");
    
    hpbpb_Jet55_Ratio[1][1][i]->SetMarkerStyle(33);
    hpbpb_Jet55_Ratio[1][1][i]->SetMarkerColor(kBlue);
    hpbpb_Jet55_Ratio[1][1][i]->Draw("same");
 
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.7,0.85,20);
    
    line55->Draw();
  }

  cJet55_Ratio->cd(1);
  putCMSPrel();

  TLegend *Jet55_ratio = myLegend(0.5,0.55,0.9,0.8);
  Jet55_ratio->AddEntry(hpbpb_Jet55[1][0][0],"Diverging with Loose Cut","pl");
  Jet55_ratio->AddEntry(hpbpb_Jet55[1][1][0],"Diverging with Tight Cut","pl");
  Jet55_ratio->SetTextSize(0.04);
  Jet55_ratio->Draw();
  
  drawText("Jet55 Trigger",0.5,0.4,14);
  
  cJet55_Ratio->cd(2);
  drawText("pCES,HBHE,|vz|<15",0.2,0.8,14);
  drawText("Ratio plots with full spectra",0.15,0.7,14);
  
  cJet55_Ratio->cd(3);
  drawText("supernova diagonal cut",0.2,0.8,14);
  
  cJet55_Ratio->SaveAs(Form("/Users/raghavke/WORK/RAA/Plots/HFVS_RAA_Effect_Jet55_ratio_%d.pdf",date.GetDate()),"RECREATE");

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  

  //
  timer.Stop();
  cout<<" Total time taken CPU(mins) = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<" Total time taken Real(mins) = "<<(Float_t)timer.RealTime()/60<<endl;


}// macro main

