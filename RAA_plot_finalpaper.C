// Raghav Kunnawalkam Elayavalli
// April 14 2014
// Rutgers
// raghav.k.e at CERN dot CH

//
// Macro to plot the final paper plots.  
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
#include "Headers/utilities.h"


using namespace std;

void RAA_plot_finalpaper(Int_t unfoldingCut = 60 , char *algo = "Pu", char *jet_type = "PF"){

    
  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  double npart[nbins_cent+1] = {389.84, 307.65, 223.95, 107.5, 41.65, 11.55, 112.9};
  
  const int nbins_pt = 30;
  const double boundaries_pt[nbins_pt+1] = {
    3, 4, 5, 7, 9, 12, 
    15, 18, 21, 24, 28,
    32, 37, 43, 49, 56,
    64, 74, 84, 97, 114,
    133, 153, 174, 196,
    220, 245, 300, 
    330, 362, 395
  };

  
  TFile *fin_R2, *fin_R3, *fin_R4; 
  fin= TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_fullMC_n20_eta_p20_%dGeVCut_ak%s_20150414.root",2,unfoldingCut,jet_type));
  fin= TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_fullMC_n20_eta_p20_%dGeVCut_ak%s_20150414.root",3,unfoldingCut,jet_type));
  fin= TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_fullMC_n20_eta_p20_%dGeVCut_ak%s_20150414.root",4,unfoldingCut,jet_type));

  // get the histograms.
  TH1F * uPbPb_R2_Bayes[nbins_cent], * uPP_R2_Bayes;
  TH1F * RAA_R2_Bayes[nbins_cent], * RAA_R3_Bayes[nbins_cent], * RAA_R4_Bayes[nbins_cent];
  TH1F * RAA_R2_BinByBin[nbins_cent], * RAA_R3_BinByBin[nbins_cent], * RAA_R4_BinByBin[nbins_cent];
  TH1F * RAA_R2_Meas[nbins_cent], * RAA_R3_Meas[nbins_cent], * RAA_R4_Meas[nbins_cent];

  
  // plot 1 - spectra plot showing pp and 6 different centrality classes PbPb spectra
  //        - have a 3 panel plot for the different radii, with each of them scaled by two orders of magnitude 

  Double_t ScaleFactor[nbins_cent+1] = {1, 1e2, 1e4, 1e6, 1e8, 1e10, 1e12};  

  TCanvas * cSpectra = new TCanvas("cSpectra","",1200,1000);
  makeMultiPanelCanvasWithGap(cSpectra,3,1,0.01,0.01,0.16,0.2,0.04,0.04);

  cSpectra->cd(1);
  cSpectra->cd(1)->SetLogy();
  cSpectra->cd(1)->SetLogx();

  uPP_R2_Bayes->Scale(ScaleFactor[0]);
  uPP_R2_Bayes->SetMarkerStyle(20);
  uPP_R2_Bayes->SetMarkerColor(2);
  makeHistTitle(uPP_R2_Bayes," "," Jet p_{T} (GeV/c)","cross section");
  uPP_R2_Bayes->SetAxisRange(unfoldingCut, 299, "X");
  
  uPP_R2_Bayes->Draw();
  for(int i = 0. i<nbins_cent; ++i){
    uPbPb_R2_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R2_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R2_Bayes[i]->SetMarkerColor(i+3);
    uPbPb_R2_Bayes[i]->Draw("same");
  }
  
  cSpectra->cd(2);
  cSpectra->cd(2)->SetLogy();
  cSpectra->cd(2)->SetLogx();

  uPP_R3_Bayes->Scale(ScaleFactor[0]);
  uPP_R3_Bayes->SetMarkerStyle(20);
  uPP_R3_Bayes->SetMarkerColor(2);
  makeHistTitle(uPP_R3_Bayes," "," Jet p_{T} (GeV/c)","cross section");
  uPP_R3_Bayes->SetAxisRange(unfoldingCut, 299, "X");
  
  uPP_R3_Bayes->Draw();
  for(int i = 0. i<nbins_cent; ++i){
    uPbPb_R3_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R3_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R3_Bayes[i]->SetMarkerColor(i+3);
    uPbPb_R3_Bayes[i]->Draw("same");
  }

  cSpectra->cd(3);
  cSpectra->cd(3)->SetLogy();
  cSpectra->cd(3)->SetLogx();

  uPP_R4_Bayes->Scale(ScaleFactor[0]);
  uPP_R4_Bayes->SetMarkerStyle(20);
  uPP_R4_Bayes->SetMarkerColor(2);
  makeHistTitle(uPP_R4_Bayes," "," Jet p_{T} (GeV/c)","cross section");
  uPP_R4_Bayes->SetAxisRange(unfoldingCut, 299, "X");
  
  uPP_R4_Bayes->Draw();
  for(int i = 0. i<nbins_cent; ++i){
    uPbPb_R4_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R4_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R4_Bayes[i]->SetMarkerColor(i+3);
    uPbPb_R4_Bayes[i]->Draw("same");
  }

  // also need to draw the MC as lines here: in all the 3 canvas. 
  
  cSpectra->SaveAs("../../Plots/Final_paper_plots_spectra.pdf","RECREATE");
  

  // plot 2 - Bayesian unfolded RAA as a function of pT for the different radii
  //        - regular 6 panel plot 
  
  // again this will be a 6 panel plot. showing measured, unfolded Bayesian, and unfolded Bin By Bin methods. 
  TCanvas *cRAA = new TCanvas("cRAA","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA = myLegend(0.45,0.75,0.85,0.9);
  TLine *lineRAA = new TLine(unfoldingCut,1,299,1);
  lineRAA->SetLineStyle(2);
  lineRAA->SetLineWidth(2);

  TLine *lUnfoldingCut = new TLine(unfoldingCut+30,0,unfoldingCut+30,2);
  lUnfoldingCut->SetLineStyle(4);
  lUnfoldingCut->SetLineWidth(2);
    
  for(int i = 0;i<nbins_cent;++i){

    cRAA->cd(nbins_cent-i);

    RAA_R2_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R2_Bayes[i]->SetMarkerStyle(24);
    makeHistTitle(RAA_R2_Bayes[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    RAA_R2_Bayes[i]->SetAxisRange(unfoldingCut,299,"X");
    RAA_R2_Bayes[i]->SetAxisRange(0,2,"Y");
    RAA_R2_Bayes[i]->Draw("E1");

    RAA_R3_Bayes[i]->SetMarkerColor(kBlack);
    RAA_R3_Bayes[i]->SetMarkerStyle(20);
    RAA_R3_Bayes[i]->Draw("same E1");

    RAA_R4_Bayes[i]->SetMarkerStyle(33);
    RAA_R4_Bayes[i]->SetMarkerColor(kRed);
    RAA_R4_Bayes[i]->Draw("same E1");

    lineRAA->Draw();
    //lUnfoldingCut->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

  }
    
  tRAA->AddEntry(RAA_R2_Bayes[0],"No Unfolding","pl");
  tRAA->AddEntry(RAA_R3_Bayes[0],"Bayesian","pl");
  tRAA->AddEntry(RAA_R4_Bayes[0],"BinbyBin","pl");
  tRAA->SetTextSize(0.04);

  cRAA->cd(1);
  tRAA->Draw();
  cRAA->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.2,0.23,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA->cd(2);
  drawText("Jet ID cut, |#eta|<2",0.1,0.3,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  cRAA->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.3,16);
  drawText("Pile up rejection cut applied",0.1,0.2,16);

  cRAA->SaveAs("../../Plots/Final_paper_plots_RAA.pdf","RECREATE");
    
  // plot 3 - RAA as a function of npart - taken from http://dde.web.cern.ch/dde/glauber_lhc.htm for 84 < pT < 97 in PbPb,PP
  //        - need to decide if we have to unfold this? or if we can just take that respective pt ranges from the already existing RAA histograms.  this is bin number 16 from the centrality classes weve measured.
  
  TCanvas * cRAA_npart = new TCanvas("cRAA_npart","",800,600);

  // get the responsible histograms for this.
  TH1F * hRAA_R2_npart = new TH1F("hRAA_R2_npart","",450, 0, 450);
  //hRAA_R2_npart->LabelsOption(">","X");
  TH1F * hRAA_R3_npart = new TH1F("hRAA_R3_npart","",450, 0, 450);
  //hRAA_R3_npart->LabelsOption(">","X");
  TH1F * hRAA_R4_npart = new TH1F("hRAA_R4_npart","",450, 0, 450);
  //hRAA_R4_npart->LabelsOption(">","X");

  for(int i = 0; i<nbins_cent; ++i){
    hRAA_R2_npart->SetBinContent(hRAA_R2_npart->FindBin(npart[i]), RAA_R2_Bayes[i]->GetBinContent[16]);
    hRAA_R3_npart->SetBinContent(hRAA_R3_npart->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent[16]);
    hRAA_R4_npart->SetBinContent(hRAA_R4_npart->FindBin(npart[i]), RAA_R4_Bayes[i]->GetBinContent[16]);    
  }

  makeHistTitle(hRAA_R2_npart, " ", "N_{part} - number of participating nucleons ", "R_{AA{");
  hRAA_R2_npart->SetMarkerColor(kRed);
  hRAA_R2_npart->SetMarkerStyle(20);
  hRAA_R2_npart->Draw();
  hRAA_R3_npart->SetMarkerColor(kBlack);
  hRAA_R3_npart->SetMarkerStyle(20);
  hRAA_R3_npart->Draw("same");
  hRAA_R4_npart->SetMarkerColor(kBlue);
  hRAA_R4_npart->SetMarkerStyle(20);
  hRAA_R4_npart->Draw("same");

  putCMSPrel();
  cRAA_npart->SaveAs("","RECREATE");
  
}
