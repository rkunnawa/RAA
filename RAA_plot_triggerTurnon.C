// Raghav Kunnawalkam Elayavalli
// Feb 26th 2015
// Rutgers
// questions or comments: raghav.k.e at CERN dot CH

//
// Macro to plot the trigger turn on curve  
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

void RAA_plot_triggerTurnon(){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TFile * fin = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_MinBiasUPC_trigger_turnoncurves_SuperNovaRejected_akPu3Calo_20150302.root");

  TH1F * hJet80_MB = (TH1F*)fin->Get("hJet80_JetMB");
  TH1F * hJet65_MB = (TH1F*)fin->Get("hJet65_JetMB");
  TH1F * hJet55_MB = (TH1F*)fin->Get("hJet55_JetMB");

  TH1F * hJetMB = (TH1F*)fin->Get("hJetMB");
  //TH1F * hJet80 = (TH1F*)fin->Get("hJet80");
  //TH1F * hJet65 = (TH1F*)fin->Get("hJet65");
  //TH1F * hJet55 = (TH1F*)fin->Get("hJet55");

  TGraphAsymmErrors *hJet80_Turnon = new TGraphAsymmErrors;
  hJet80_Turnon->BayesDivide(hJet80_MB,hJetMB);
  TGraphAsymmErrors *hJet65_Turnon = new TGraphAsymmErrors;
  hJet65_Turnon->BayesDivide(hJet65_MB,hJetMB);
  TGraphAsymmErrors *hJet55_Turnon = new TGraphAsymmErrors;
  hJet55_Turnon->BayesDivide(hJet55_MB,hJetMB);

  //TH1F * hTurnon_Jet55 = (TH1F*)hJet55_MB->Clone("hTurnon_Jet55");
  //hTurnon_Jet55->Divide(hJetMB);
 
  
  //plot the turn on curves:

  TCanvas * cturnon = new TCanvas ("cTurnon","",800,600);
  cturnon->SetLogy();
  //hJet80_Turnon->SetTitle(" ");
  //hJet80_Turnon->SetXTitle("Jet p_{T} (GeV/c)");
  //hJet80_Turnon->SetYTitle("Trigger Efficiency");
  //hJet80_Turnon->SetMarkerStyle(20);
  //hJet80_Turnon->SetMarkerColor(kBlack);
  //hJet80_Turnon->Rebin(5);
  //hJet80_Turnon->Scale(1./5);
  //hJet80_Turnon->SetAxisRange(20,170,"X");
  //hJet80_Turnon->SetAxisRange(0,1.2,"Y");
  hJetMB->SetMarkerStyle(24);
  hJetMB->Draw();
  hJet80_MB->SetMarkerStyle(25);
  hJet80_MB->SetMarkerColor(kBlue);
  hJet80_MB->Draw("same");
  hJet65_MB->SetMarkerStyle(26);
  hJet65_MB->SetMarkerColor(kRed);
  hJet65_MB->Draw("same");
  hJet55_MB->SetMarkerStyle(27);
  hJet55_MB->SetMarkerColor(kGreen);
  hJet55_MB->Draw("same");
  

  //hJet65_Turnon->SetMarkerStyle(20);
  //hJet65_Turnon->SetMarkerColor(kGreen);
  // //hJet65_Turnon->Rebin(5);
  // //hJet65_Turnon->Scale(1./5);
  //hJet65_Turnon->Draw("Ap");

  // Jet55_Turnon->SetMarkerStyle(20);
  // Jet55_Turnon->SetMarkerColor(kBlue);
  // // Jet55_Turnon->Rebin(5);
  // // Jet55_Turnon->Scale(1./5);
  // Jet55_Turnon->Draw("Ap");


  TLine * line = new TLine(20,1,200,1);
  line->SetLineWidth(2);
  //line->Draw();

  TLegend * turnon = myLegend(0.5,0.2,0.7,0.4);
  turnon->AddEntry(hJetMB,"MB","pl");
  turnon->AddEntry(hJet80_MB,"80 MB","pl");
  turnon->AddEntry(hJet65_MB,"65 MB","pl");
  turnon->AddEntry(hJet55_MB,"55 MB","pl");
  //turnon->AddEntry(hJet55_Turnon,"HLT_HIJet55_v1","pl");
  //turnon->AddEntry(hJet65_Turnon,"HLT_HIJet65_v1","pl");
  //turnon->AddEntry(hJet80_Turnon,"HLT_HIJet80_v1","pl");
  turnon->SetTextSize(0.04);
  turnon->Draw();
  
  putCMSPrel();
  drawText("HLT turnon curves",0.2,0.8,16);

  cturnon->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_triggerturnon_%d.pdf",date.GetDate()),"RECREATE");
  

}
