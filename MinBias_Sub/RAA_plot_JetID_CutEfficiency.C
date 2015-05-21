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

using namespace std;

void RAA_plot_JetID_CutEfficiency(  Int_t radius = 4,
				    char * etaWidth = (char*)"20_eta_p20"){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  char * Coll = "PbPb";
  //char * Coll = "Pp";

  TDatime date;

  Int_t trigger[3]; 

  if(Coll == "PbPb") {trigger[0] = 55; trigger[1] = 65; trigger[2] = 80;}
  if(Coll == "Pp") {trigger[0] = 40; trigger[1] = 60; trigger[2] =  80;}

  char * jetType; 

  if(Coll == "PbPb") jetType = Form("akPu%dPF",radius);
  if(Coll == "Pp") jetType = Form("ak%dPF",radius);

  const int nbins_pt = 38;
  const double boundaries_pt[nbins_pt+1] = {
    3, 4, 5, 7, 9, 12, 
    15, 18, 21, 24, 28,
    32, 37, 43, 49, 56,
    64, 74, 84, 97, 114,
    133, 153, 174, 196,
    220, 245, 300, 
    330, 362, 395, 430,
    468, 507, 548, 592,
    638, 686, 1000 
  };
  
  // Pawan's files:
  TFile * fIn = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_ntuple_PbPb_data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius));
//= TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pawan_ntuple_PbPb_data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth, radius)); 
  TFile * fPP_in = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_ntuple_PP_data_MC_spectra_residualFactor_finebins_%s_R0p%d.root",etaWidth, radius));
  
  // get the histograms, the Cut Efficiency is plotted from the MC histograms 
  TH1F * hMC_noCut[3], * hMC_Cut[3], * hMC_unm_noCut[3], * hMC_unm_Cut[3]; 

  TH1F * hMC_Denominator, *hMC_Numerator;

  for(int i = 0; i<3; i++){

    hMC_noCut[i] = (TH1F*)fIn->Get(Form("hMC_Jet%d_noCut", trigger[i]));
    hMC_Cut[i] = (TH1F*)fIn->Get(Form("hMC_Jet%d_CutA", trigger[i]));
    hMC_unm_noCut[i] = (TH1F*)fIn->Get(Form("hMC_unmatched_Jet%d_noCut", trigger[i]));
    hMC_unm_Cut[i] = (TH1F*)fIn->Get(Form("hMC_unmatched_Jet%d_CutA", trigger[i]));

  }

  hMC_Denominator = (TH1F*)hMC_noCut[0]->Clone("hMC_Denominator");
  hMC_Denominator->Add(hMC_noCut[1]);
  hMC_Denominator->Add(hMC_noCut[2]);
  hMC_Denominator->Add(hMC_unm_noCut[0]);
  hMC_Denominator->Add(hMC_unm_noCut[1]);
  hMC_Denominator->Add(hMC_unm_noCut[2]);
  
  hMC_Numerator = (TH1F*)hMC_Cut[0]->Clone("hMC_Numerator");
  hMC_Numerator->Add(hMC_Cut[1]);
  hMC_Numerator->Add(hMC_Cut[2]);
  hMC_Numerator->Add(hMC_unm_Cut[0]);
  hMC_Numerator->Add(hMC_unm_Cut[1]);
  hMC_Numerator->Add(hMC_unm_Cut[2]);

  TH1F * hCutEff = (TH1F*)hMC_Numerator->Clone("hCutEff");
  hCutEff->Divide(hMC_Denominator);
  hCutEff = (TH1F*)hCutEff->Rebin(nbins_pt, "hCutEff", boundaries_pt);
  divideBinWidth(hCutEff);

  // line at 1

  
  TCanvas * cCutEff = new TCanvas("cCutEff","",800,600);
  hCutEff->SetXTitle(Form("%s Gen p_{T} (GeV/c)",jetType));
  hCutEff->SetYTitle("Jet ID Cut Efficiency");
  hCutEff->SetTitle(" ");
  hCutEff->SetAxisRange(40, 300, "X");
  hCutEff->SetAxisRange(0.9, 1.05, "Y");
  hCutEff->SetMarkerStyle(20);
  hCutEff->SetMarkerColor(kBlack);
  hCutEff->Draw();

  putCMSPrel();
  drawText("|#eta|<2, |vz|<15", 0.6, 0.31, 20);
  drawText(Form("%s",Coll), 0.6, 0.21, 20);

  cCutEff->SaveAs(Form("May20/%s_Combined_CutEfficiency_R0p%d_%s_%d.pdf",Coll,radius,etaWidth,date.GetDate()),"RECREATE");
  


  // plot the curves for the different centrality classes. get the combined spectra and do the ratio 
  
  const int nbins_cent = 6; 

  TH1F * hMC_Denominator[nbins_cent];
  TH1F * hMC_Numerator[nbins_cent]; 

  TH1F * hData_Denominator[nbins_cent];
  TH1F * hData_Numerator[nbins_cent]; 

  TH1F * hCutEff[nbins_cent];
  TH1F * hCutEff_Data[nbins_cent];

  TH1F * hPP_MC_Denominator = (TH1F*)fIn_PP->Get("hpp_MC_Comb_noCut");
  TH1F * hPP_MC_Numerator = (TH1F*)fIn_PP->Get(Form("hpp_JetComb_gen_R%d_%s",radius,etaWidth)); 

  TH1F * hPP_Data_Denominator = (TH1F*)fIn_PP->Get("hpp_Data_Comb_noCut");
  TH1F * hPP_Data_Numerator = (TH1F*)fIn_PP->Get(Form("hpp_HLTComb_R%d_%s",radius,etaWidth)); 

  TH1F * hPP_CutEff = (TH1F*)hPP_MC_Numerator->Clone("hPP_CutEff");
  hPP_CutEff->Divide(hPP_MC_Denominator);
  hPP_CutEff->Rebin(10);
  hPP_CutEff->Scale(1./10);
  TH1F * hPP_CutEff_Data = (TH1F*)hPP_Data_Numerator->Clone("hPP_CutEff_Data");
  hPP_CutEff_Data->Divide(hPP_Data_Denominator);
  hPP_CutEff_Data->Rebin(10);
  hPP_CutEff_Data->Scale(1./10);
  
  TLine *line = new TLine(60,1,299,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  for(int i = 0; i<nbins_cent; ++i){
    
    hMC_Denominator[i] = (TH1F*)fIn->Get(Form("hpbpb_MC_noCut_cent%d",i));
    hMC_Numerator[i] = (TH1F*)fIn->Get(Form("hpbpb_gen_R%d_%s_cent%d",radius,etaWidth, i));
    hData_Denominator[i] = (TH1F*)fIn->Get(Form("hpbpb_Data_Comb_noCut_cent%d",i));
    hData_Numerator[i] = (TH1F*)fIn->Get(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i));

    hCutEff[i] = (TH1F*)hMC_Numerator[i]->Clone(Form("CutEff_cent%d",i));
    hCutEff[i]->Divide(hMC_Denominator[i]);
    hCutEff[i]->Rebin(20);
    hCutEff[i]->Scale(1./20);

    hCutEff_Data[i] = (TH1F*)hData_Numerator[i]->Clone(Form("CutEff_Data_cent%d",i));
    hCutEff_Data[i]->Divide(hData_Denominator[i]);
    hCutEff_Data[i]->Rebin(20);
    hCutEff_Data[i]->Scale(1./20);

  }

  TCanvas * cCutEff = new TCanvas("cCutEff","",800,600);
  makeMultiPanelCanvas(cCutEff,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0; i<nbins_cent; ++i){
    cCutEff->cd(nbins_cent-i);
    hCutEff[i]->SetXTitle(Form("akPu%dPF Gen p_{T} (GeV/c)", radius));
    hCutEff[i]->SetYTitle("Jet ID Cut Efficiency");
    hCutEff[i]->SetTitle(" ");
    hCutEff[i]->SetAxisRange(50, 299, "X");
    hCutEff[i]->SetAxisRange(0.95, 1.05, "Y");
    hCutEff[i]->SetMarkerStyle(20);
    hCutEff[i]->SetMarkerColor(kBlack);
    hCutEff[i]->Draw();
    line->Draw();
  }

  putCMSPrel();
  drawText("|#eta|<2, |vz|<15", 0.6, 0.31, 20);
  drawText(Form("%s",Coll), 0.6, 0.21, 20);

  cCutEff->SaveAs(Form("May20/%s_Combined_CutEfficiency_centralityClass_R0p%d_%s_%d.pdf",Coll,radius,etaWidth,date.GetDate()),"RECREATE");

  TCanvas * cCutEff_Data = new TCanvas("cCutEff_Data","",800,600);
  makeMultiPanelCanvas(cCutEff_Data,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i = 0; i<nbins_cent; ++i){
    cCutEff_Data->cd(nbins_cent-i);
    hCutEff_Data[i]->SetXTitle(Form("akPu%dPF Data reco p_{T} (GeV/c)", radius));
    hCutEff_Data[i]->SetYTitle("Jet ID Cut Efficiency");
    hCutEff_Data[i]->SetTitle(" ");
    hCutEff_Data[i]->SetAxisRange(60, 299, "X");
    hCutEff_Data[i]->SetAxisRange(0.2, 1.1, "Y");
    hCutEff_Data[i]->SetMarkerStyle(20);
    hCutEff_Data[i]->SetMarkerColor(kBlack);
    hCutEff_Data[i]->Draw();
    line->Draw();
  }

  putCMSPrel();
  drawText("|#eta|<2, |vz|<15", 0.6, 0.31, 20);
  drawText(Form("%s",Coll), 0.6, 0.21, 20);

  cCutEff_Data->SaveAs(Form("May20/%s_Combined_CutEfficiency_Data_centralityClass_R0p%d_%s_%d.pdf",Coll,radius,etaWidth,date.GetDate()),"RECREATE");

  
  TCanvas * cPP_CutEff = new TCanvas("cPP_CutEff","",800,600);
  hPP_CutEff->SetXTitle(Form("ak%dPF Gen p_{T} (GeV/c)", radius));
  hPP_CutEff->SetYTitle("Jet ID Cut Efficiency");
  hPP_CutEff->SetTitle(" ");
  hPP_CutEff->SetAxisRange(60, 299, "X");
  hPP_CutEff->SetAxisRange(0.8,1.1,"Y");
  hPP_CutEff->SetMarkerStyle(20);
  hPP_CutEff->Draw();

  cPP_CutEff->SaveAs(Form("May20/PP_Combined_CutEfficiency_R0p%d_%d.pdf",radius, date.GetDate()));

  TCanvas * cPP_CutEff_Data = new TCanvas("cPP_CutEff_Data","",800,600);
  hPP_CutEff_Data->SetXTitle(Form(" ak%dPF Data reco p_{T} (GeV/c)", radius));
  hPP_CutEff_Data->SetYTitle("Jet ID Cut Efficiency");
  hPP_CutEff_Data->SetTitle(" ");
  hPP_CutEff_Data->SetAxisRange(60, 299, "X");
  hPP_CutEff_Data->SetAxisRange(0.8,1.1,"Y");
  hPP_CutEff_Data->SetMarkerStyle(20);
  hPP_CutEff_Data->Draw();

  cPP_CutEff_Data->SaveAs(Form("May20/PP_Combined_CutEfficiency_Data_R0p%d_%s_%d.pdf",radius,etaWidth, date.GetDate()));
  
}

