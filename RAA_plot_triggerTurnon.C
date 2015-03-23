
#include <iostream>
#include <stdio.h>
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
#include <TMultiGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

using namespace std;


void RAA_plot_triggerTurnon(){

  TH1::SetDefaultSumw2();

  TFile * fin = TFile::Open("/Users/keraghav/WORK/RAA/Output/PbPb_MinBiasUPC_ntuple_SuperNovaRejected_akPuCalo_20150311.root");

  TTree * trig = (TTree*)fin->Get("trigger_info");

  int run, evt, lumi, nref, jet80, jet80_prescl, jet65, jet65_prescl, jet55, jet55_prescl, hiMB, hiMB_prescl;
  float jetpt[1000], jeteta[1000];

  trig->SetBranchAddress("evt_value",&evt);
  trig->SetBranchAddress("run_value",&run);
  trig->SetBranchAddress("lumi_value",&lumi);
  trig->SetBranchAddress("nref",&nref);
  trig->SetBranchAddress("jetpt",&jetpt);
  trig->SetBranchAddress("jeteta",&jeteta);
  trig->SetBranchAddress("Jet80",&jet80);
  trig->SetBranchAddress("Jet80_prescl",&jet80_prescl);
  trig->SetBranchAddress("Jet65",&jet65);
  trig->SetBranchAddress("Jet65_prescl",&jet65_prescl);
  trig->SetBranchAddress("Jet55",&jet55);
  trig->SetBranchAddress("Jet55_prescl",&jet55_prescl);
  trig->SetBranchAddress("HIMinBias",&hiMB);
  trig->SetBranchAddress("HIMinBias_prescl",&hiMB_prescl);

  TH1F* hDenominator_80 = new TH1F("hDenominator_80","",14,0,140);
  TH1F* hDenominator_65 = new TH1F("hDenominator_65","",14,0,140);
  TH1F* hDenominator_55 = new TH1F("hDenominator_55","",14,0,140);
  TH1F* hNumerator_80 = new TH1F("hNumerator_80","",14,0,140);
  TH1F* hNumerator_65 = new TH1F("hNumerator_65","",14,0,140);
  TH1F* hNumerator_55 = new TH1F("hNumerator_55","",14,0,140);

  
  Long_t  nentries = trig->GetEntries();
  cout<<nentries<<endl;
  // nentries = 10000000;
  for(Long_t i = 0;i<nentries;i++){
    trig->GetEntry(i);
    if(i%100000 == 0) cout<<i<<"/"<<nentries<<endl;
   
    if(jet80_prescl==1 && jet65_prescl==1 && jet55_prescl==1){

      hDenominator_80->Fill(jetpt[0]);
      hDenominator_65->Fill(jetpt[0]);
      hDenominator_55->Fill(jetpt[0]);
      
      if(jet80)
	hNumerator_80->Fill(jetpt[0]);
      if(jet65)
	hNumerator_65->Fill(jetpt[0]);
      if(jet55)
	hNumerator_55->Fill(jetpt[0]);
      
    }

  }
  
  hDenominator_80->Print("base");
  hDenominator_65->Print("base");
  hDenominator_55->Print("base");
  hNumerator_80->Print("base");
  hNumerator_65->Print("base");
  hNumerator_55->Print("base");
  
  TGraphAsymmErrors * Jet80 = new TGraphAsymmErrors();
  Jet80->BayesDivide(hNumerator_80, hDenominator_80);
  Jet80->SetMarkerColor(kRed);
  Jet80->SetMarkerSize(2.0);
  TGraphAsymmErrors * Jet65 = new TGraphAsymmErrors();
  Jet65->BayesDivide(hNumerator_65, hDenominator_65);
  Jet65->SetMarkerColor(kGreen);
  Jet65->SetMarkerSize(2.0);
  TGraphAsymmErrors * Jet55 = new TGraphAsymmErrors();
  Jet55->BayesDivide(hNumerator_55, hDenominator_55);
  Jet55->SetMarkerColor(kBlue);
  Jet55->SetMarkerSize(2.0);
  
  TCanvas * cSpectra = new TCanvas ("cSpectra","",1000,600);
  cSpectra->Divide(2,1);
  cSpectra->cd(1)->SetLogy();

  hDenominator_80->SetMarkerColor(kRed);
  hDenominator_80->SetMarkerStyle(25);
  hDenominator_80->Draw();
  
  hDenominator_65->SetMarkerColor(kGreen);
  hDenominator_65->SetMarkerStyle(24);
  hDenominator_65->Draw("same");
  
  hDenominator_55->SetMarkerColor(kBlue);
  hDenominator_55->SetMarkerStyle(27);
  hDenominator_55->Draw("same");

  hNumerator_80->SetMarkerColor(kRed);
  hNumerator_80->SetMarkerStyle(21);
  hNumerator_80->Draw("same");
  hNumerator_65->SetMarkerColor(kGreen);
  hNumerator_65->SetMarkerStyle(20);
  hNumerator_65->Draw("same");
  hNumerator_55->SetMarkerColor(kBlue);
  hNumerator_55->SetMarkerStyle(33);
  hNumerator_55->Draw("same");

  cSpectra->cd(2);
    
  TMultiGraph * mg = new TMultiGraph();
  mg->SetTitle(";lead jet p_{T}^{reco};Efficiency");
  mg->Add(Jet80,"");
  mg->Add(Jet65,"");
  mg->Add(Jet55,"");

  mg->Draw("AP");

  cSpectra->SaveAs("Trigger_spectra_MinBias.pdf","RECREATE");
  

}
