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
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 

static const int nRadii = 3;
static const Int_t list_Radii[nRadii] = {2,3,4};

using namespace std;


void RAA_plot_PawanvsRaghav_ntuple(char* etaWidth = (char*)"20_eta_20"){

  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  TDatime date;

  
  // get the input files.
  TFile * fPawan_PbPb[nRadii], * fRaghav_PbPb[nRadii];
  TFile * fPawan_PP[nRadii], * fRaghav_PP[nRadii];

  for(int k = 0; k<nRadii; ++k){
    fPawan_PbPb[k]  = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_ntuple_PbPb_Data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root", etaWidth, list_Radii[k]));
    fRaghav_PbPb[k] = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Raghav_ntuple_PbPb_Data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root", etaWidth, list_Radii[k]));
    fPawan_PP[k]  = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_ntuple_PP_data_MC_spectra_residualFactor_finebins_%s_R0p%d.root", etaWidth, list_Radii[k]));
    fRaghav_PP[k] = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Raghav_ntuple_PP_data_MC_spectra_residualFactor_finebins_%s_R0p%d.root", etaWidth, list_Radii[k]));
  }
  
  // declare the histograms and get them from the root files

  TH1F * hPbPb_Pawan[nRadii][nbins_cent], * hPbPb_Raghav[nRadii][nbins_cent], * hPbPb_Ratio[nRadii][nbins_cent]; // ratio = pawan/raghav
  TH1F * hPP_Pawan[nRadii], * hPP_Raghav[nRadii], * hPP_Ratio[nRadii]; // ratio = pawan/raghav

  char * histType_PbPb = "Data_Comb_noCut";
  // can be JetComb_gen
  //        JetComb_reco
  //        HLTComb

  char * histType_pp = "HLTComb";
  // can be JetComb_gen
  //        JetComb_reco
  //        HLTComb
  
  
  for(int k = 0; k<nRadii; ++k){
    cout<<"Radius = "<<list_Radii[k]<<endl;
    hPP_Pawan[k] = (TH1F*)fPawan_PP[k]->Get(Form("hpp_%s_R%d_%s",histType_pp, list_Radii[k], etaWidth));
    hPP_Pawan[k]->Print("base");
    hPP_Raghav[k] = (TH1F*)fRaghav_PP[k]->Get(Form("hpp_%s_R%d_%s",histType_pp,list_Radii[k], etaWidth));
    hPP_Raghav[k]->Print("base");
    hPP_Ratio[k] = (TH1F*)hPP_Pawan[k]->Clone(Form("hPP_PawanOverRaghav_R%d",list_Radii[k]));
    hPP_Ratio[k]->Divide(hPP_Raghav[k]);
    hPP_Ratio[k] = (TH1F*)hPP_Ratio[k]->Rebin(10);
    hPP_Ratio[k]->Scale(1./10);
    
    for(int i = 0; i<nbins_cent; ++i){
      cout<<"cent = "<<i<<endl;
      hPbPb_Pawan[k][i] = (TH1F*)fPawan_PbPb[k]->Get(Form("hpbpb_%s_cent%d",histType_PbPb,i));
      hPbPb_Pawan[k][i]->Print("base");
      hPbPb_Raghav[k][i] = (TH1F*)fRaghav_PbPb[k]->Get(Form("hpbpb_%s_cent%d",histType_PbPb,i));
      hPbPb_Raghav[k][i]->Print("base");
      hPbPb_Ratio[k][i] = (TH1F*)hPbPb_Pawan[k][i]->Clone(Form("hPbPb_PawanOverRaghav_R%d_cent%d",list_Radii[k],i));
      hPbPb_Ratio[k][i]->Divide(hPbPb_Raghav[k][i]);
      hPbPb_Ratio[k][i] = (TH1F*)hPbPb_Ratio[k][i]->Rebin(10);
      hPbPb_Ratio[k][i]->Scale(1./10);

    }
  }

  // make the plots

  TCanvas * cRatio = new TCanvas("cRatio","",1000,800);
  makeMultiPanelCanvas(cRatio,3,2,0.0,0.0,0.2,0.15,0.07);
  TLine * line = new TLine(15,1,399,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TLegend *leg = myLegend(0.6,0.7,0.8,0.9);
  
  for(int i = 0; i<nbins_cent; ++i){
    cRatio->cd(nbins_cent-i);

    for(int k = 0; k<nRadii; ++k){
      makeHistTitle(hPbPb_Ratio[k][i], " ", "Jet p_{T} (GeV/c)", "#frac{Pawan}{Raghav}");
      hPbPb_Ratio[k][i]->SetMarkerStyle(24+k);
      hPbPb_Ratio[k][i]->SetMarkerColor(2+k);
      hPbPb_Ratio[k][i]->SetAxisRange(15,299,"X");
      hPbPb_Ratio[k][i]->SetAxisRange(0,2,"Y");
      if(k==0) hPbPb_Ratio[k][i]->Draw("p");
      else hPbPb_Ratio[k][i]->Draw("same p");
      if(i==0)
	leg->AddEntry(hPbPb_Ratio[k][i],Form("R=0.%d",list_Radii[k]),"pl");
      
    }
    if(i==nbins_cent-1) leg->Draw();
    line->Draw();
    drawText(Form("%2.0f-%.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }
  
  cRatio->SaveAs(Form("Pawan_vs_Raghav_ntuple_PbPb_spectra_%s_ratio_%s_%d.pdf", histType_PbPb, etaWidth, date.GetDate()),"RECREATE");

  TCanvas * cRatio_pp = new TCanvas("cRatio_pp","",600,600);
  TLegend *leg_pp = myLegend(0.6,0.7,0.8,0.9);
  for(int k = 0; k<nRadii; ++k){
    makeHistTitle(hPP_Ratio[k], " ", "Jet p_{T} (GeV/c)", "#frac{Pawan}{Raghav}");
    hPP_Ratio[k]->SetMarkerStyle(24+k);
    hPP_Ratio[k]->SetMarkerColor(2+k);
    hPP_Ratio[k]->SetAxisRange(15,299,"X");
    hPP_Ratio[k]->SetAxisRange(0,2,"Y");
    if(k==0) hPP_Ratio[k]->Draw("p");
      else hPP_Ratio[k]->Draw("same p");
    leg_pp->AddEntry(hPP_Ratio[k],Form("R=0.%d",list_Radii[k]),"pl");
    
  }
  leg->Draw();
  line->Draw();
  cRatio_pp->SaveAs(Form("Pawan_vs_Raghav_ntuple_PP_spectra_%s_ratio_%s_%d.pdf",histType_pp,etaWidth, date.GetDate()),"RECREATE");

}
