// Raghav Kunnawalkam Elayavalliu
// Dec 4th 2014
// CERN

// macro to check the rerecoed hf/Vs algorithm by plotting the two jet spectras on top of each other. 
//


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

#include "Headers/plot.h"

using namespace std;

void RAA_VsValidationCheck(){

  TH1F::SetDefaultSumw2();

  // get the input files: 
  TFile *fOldData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_data_akPuPF_testComb4_cut1_20141112.root");
  TFile *fReReco = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/event_list_data_akvs3pf_v0linv2_unequalized.root");

  TH1F *hOldJet55 = new TH1F("hOldJet55","",500,0,1000);
  TH1F *hReReco55 = new TH1F("hReReco55","",500,0,1000);

  TTree *jetOldData = (TTree*)fOldData->Get("jets_ID");
  TTree *jetReReco = (TTree*)fReReco->Get("t");

  jetOldData->Project("hOldJet55","pt","jet55 && cent==0");
  jetReReco->Project("hReReco55","jtpt");

  TCanvas *c = new TCanvas("c","",600,400);
  makeMultiPanelCanvas(c,1,2,0.0,0.0,0.2,0.15,0.07);

  c->cd(1);
  c->cd(1)->SetLogy();
  hOldJet55->Print("base");
  hReReco55->Print("base");
  
  //hOldJet55->Scale(1./hOldJet55->Integral());
  //hReReco55->Scale(1./hReReco55->Integral());

  hOldJet55->SetMarkerStyle(20);
  hOldJet55->SetMarkerColor(kBlack);
  hOldJet55->SetXTitle("Jet55, central bin 0-5%, Jet pT");
  hOldJet55->SetYTitle("normalized to hist integral");
  hOldJet55->Draw();
  
  hReReco55->SetMarkerStyle(24);
  hReReco55->SetMarkerColor(kBlack);
  hReReco55->Draw("same");
  
  TLegend *PbPbSpectra = myLegend(0.35,0.65,0.85,0.9);
  PbPbSpectra->AddEntry(hOldJet55,"Old Jet 55","pl");
  PbPbSpectra->AddEntry(hReReco55,"New ReReco sample, Poly v0, Linv2","pl");
  PbPbSpectra->SetTextSize(0.04);
  PbPbSpectra->Draw();

  c->cd(2);
  TH1F *hRatio = (TH1F*)hOldJet55->Clone("hRatio");
  hRatio->Divide(hReReco55);
  hRatio->SetMarkerStyle(33);
  hRatio->SetYTitle("Existing/new ReReco");
  hRatio->SetXTitle("Jet pT");
  hRatio->Draw();
  
  c->SaveAs("oldvsrereco.pdf","RECREATE");
  c->SaveAs("oldvsrereco.C","RECREATE");

}
