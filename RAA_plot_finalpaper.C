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

  TFile *fin_R2, *fin_R3, *fin_R4; 
  fin_R2 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_trgMC_n20_eta_p20_%dGeVCut_ak%s_20150415.root",2,unfoldingCut,jet_type));
  fin_R3 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_trgMC_n20_eta_p20_%dGeVCut_ak%s_20150415.root",3,unfoldingCut,jet_type));
  fin_R4 = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_trgMC_n20_eta_p20_%dGeVCut_ak%s_20150415.root",4,unfoldingCut,jet_type));

  // get the histograms.
  TH1F * uPbPb_R2_Bayes[nbins_cent], * uPP_R2_Bayes, * uPbPb_R3_Bayes[nbins_cent], * uPP_R3_Bayes, * uPbPb_R4_Bayes[nbins_cent], * uPP_R4_Bayes;
  TH1F * mPbPb_R2[nbins_cent], * mPP_R2, * mPbPb_R3[nbins_cent], * mPP_R3, * mPbPb_R4[nbins_cent], * mPP_R4;
  
  TH1F * RAA_R2_Bayes[nbins_cent], * RAA_R3_Bayes[nbins_cent], * RAA_R4_Bayes[nbins_cent];
  TH1F * RAA_R2_BinByBin[nbins_cent], * RAA_R3_BinByBin[nbins_cent], * RAA_R4_BinByBin[nbins_cent];
  TH1F * RAA_R2_Meas[nbins_cent], * RAA_R3_Meas[nbins_cent], * RAA_R4_Meas[nbins_cent];

  uPP_R2_Bayes = (TH1F*)fin_R2->Get("PP_byesian_unfoded_spectra");
  uPP_R3_Bayes = (TH1F*)fin_R3->Get("PP_byesian_unfoded_spectra");
  uPP_R4_Bayes = (TH1F*)fin_R4->Get("PP_byesian_unfoded_spectra");

  mPP_R2 = (TH1F*)fin_R2->Get("PP_Reco_spectra_refpt");
  mPP_R3 = (TH1F*)fin_R3->Get("PP_Reco_spectra_refpt");
  mPP_R4 = (TH1F*)fin_R4->Get("PP_Reco_spectra_refpt");
  
  for(int i = 0; i<nbins_cent; ++i){

    uPbPb_R2_Bayes[i] = (TH1F*)fin_R2->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R3_Bayes[i] = (TH1F*)fin_R3->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_R4_Bayes[i] = (TH1F*)fin_R4->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));

    mPbPb_R2[i] = (TH1F*)fin_R2->Get(Form("PbPb_Reco_spectra_refpt_cent%d",i));
    mPbPb_R3[i] = (TH1F*)fin_R3->Get(Form("PbPb_Reco_spectra_refpt_cent%d",i));
    mPbPb_R4[i] = (TH1F*)fin_R4->Get(Form("PbPb_Reco_spectra_refpt_cent%d",i));

    RAA_R2_Bayes[i]   = (TH1F*)fin_R2->Get(Form("RAA_bayesian_cent%d",i));  
    RAA_R3_Bayes[i]   = (TH1F*)fin_R3->Get(Form("RAA_bayesian_cent%d",i));  
    RAA_R4_Bayes[i]   = (TH1F*)fin_R4->Get(Form("RAA_bayesian_cent%d",i));  
    
    RAA_R2_BinByBin[i]   = (TH1F*)fin_R2->Get(Form("RAA_binbybin_cent%d",i));  
    RAA_R3_BinByBin[i]   = (TH1F*)fin_R3->Get(Form("RAA_binbybin_cent%d",i));  
    RAA_R4_BinByBin[i]   = (TH1F*)fin_R4->Get(Form("RAA_binbybin_cent%d",i));  
    
    RAA_R2_Meas[i]   = (TH1F*)fin_R2->Get(Form("RAA_measured_cent%d",i));  
    RAA_R3_Meas[i]   = (TH1F*)fin_R3->Get(Form("RAA_measured_cent%d",i));  
    RAA_R4_Meas[i]   = (TH1F*)fin_R4->Get(Form("RAA_measured_cent%d",i));  
    
  }
  
  // plot 1 - spectra plot showing pp and 6 different centrality classes PbPb spectra
  //        - have a 3 panel plot for the different radii, with each of them scaled by two orders of magnitude 

  // first we need to scale the MC to the level of Data:
  // PbPb Data scaling:
  //   uPbPb_Bayes[i]->Scale(1./deltaEta);// delta eta
  //   //uPbPb_Bayes[i]->Scale(1./145.156/1e6);// Jet 80 luminosity
  //   //uPbPb_Bayes[i]->Scale(1./1.1153/1e6);// equivalent no of minbias events 
  //   uPbPb_Bayes[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
  //   //uPbPb_Bayes[i]->Scale(1./145.156);
  //   //uPbPb_Bayes[i]->Scale(1./161.939);
  //   uPbPb_Bayes[i]->Scale(1./(7.65*1e6));
  //   uPbPb_Bayes[i]->Scale(64.*1e9/(ncoll[i]*1e3));
  //   uPbPb_Bayes[i] = (TH1F*)uPbPb_Bayes[i]->Rebin(nbins_pt,Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i),boundaries_pt);
  //   divideBinWidth(uPbPb_Bayes[i]);
  //   uPbPb_Bayes[i]->Write();
  //   So finally PbPb is in 

  // PbPb MC scaling: is already in sigma (mb) / (dEta dpT)
  //   mPbPb_Reco[i]->Scale(1./deltaEta);// delta eta
  //   mPbPb_Reco[i] = (TH1F*)mPbPb_Reco[i]->Rebin(nbins_pt,Form("PbPb_Reco_spectra_refpt_cent%d",i),boundaries_pt);
  //   divideBinWidth(mPbPb_Reco[i]);
  //   mPbPb_Reco[i]->Write();

  mPP_R2->Scale(1e6);
  mPP_R3->Scale(1e6);
  mPP_R4->Scale(1e6);

  for(int i = 0; i<nbins_cent; ++i){

    mPbPb_R2[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    mPbPb_R2[i]->Scale(64.*1e9/(ncoll[i]*1e3));
    mPbPb_R2[i]->Scale(1./(7.65));

    mPbPb_R3[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    mPbPb_R3[i]->Scale(64.*1e9/(ncoll[i]*1e3));
    mPbPb_R3[i]->Scale(1./(7.65));

    mPbPb_R4[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
    mPbPb_R4[i]->Scale(64.*1e9/(ncoll[i]*1e3));
    mPbPb_R4[i]->Scale(1./(7.65));
    
  }
  
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
  
  // draw the MC
  mPP_R2->Scale(ScaleFactor[0]);
  mPP_R2->SetLineStyle(2);
  mPP_R2->SetLineColor(2);
  mPP_R2->Draw("same L");
  
  uPP_R2_Bayes->Draw();
  for(int i = 0; i<nbins_cent; ++i){
    uPbPb_R2_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R2_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R2_Bayes[i]->SetMarkerColor(i+3);
    uPbPb_R2_Bayes[i]->Draw("same");

    mPbPb_R2[i]->Scale(ScaleFactor[i]);
    mPbPb_R2[i]->SetLineStyle(2);
    mPbPb_R2[i]->SetLineColor(i+3);
    mPbPb_R2[i]->Draw("same L");
    
  }
  
  cSpectra->cd(2);
  cSpectra->cd(2)->SetLogy();
  cSpectra->cd(2)->SetLogx();

  uPP_R3_Bayes->Scale(ScaleFactor[0]);
  uPP_R3_Bayes->SetMarkerStyle(20);
  uPP_R3_Bayes->SetMarkerColor(2);
  makeHistTitle(uPP_R3_Bayes," "," Jet p_{T} (GeV/c)","cross section");
  uPP_R3_Bayes->SetAxisRange(unfoldingCut, 299, "X");
  
  // draw the MC
  mPP_R3->Scale(ScaleFactor[0]);
  mPP_R3->SetLineStyle(2);
  mPP_R3->SetLineColor(2);
  mPP_R3->Draw("same L");
  
  uPP_R3_Bayes->Draw();
  for(int i = 0; i<nbins_cent; ++i){
    uPbPb_R3_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R3_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R3_Bayes[i]->SetMarkerColor(i+3);
    uPbPb_R3_Bayes[i]->Draw("same");

    mPbPb_R3[i]->Scale(ScaleFactor[i]);
    mPbPb_R3[i]->SetLineStyle(2);
    mPbPb_R3[i]->SetLineColor(i+3);
    mPbPb_R3[i]->Draw("same L");
  }

  cSpectra->cd(3);
  cSpectra->cd(3)->SetLogy();
  cSpectra->cd(3)->SetLogx();

  uPP_R4_Bayes->Scale(ScaleFactor[0]);
  uPP_R4_Bayes->SetMarkerStyle(20);
  uPP_R4_Bayes->SetMarkerColor(2);
  makeHistTitle(uPP_R4_Bayes," "," Jet p_{T} (GeV/c)","cross section");
  uPP_R4_Bayes->SetAxisRange(unfoldingCut, 299, "X");

  // draw the MC
  mPP_R4->Scale(ScaleFactor[0]);
  mPP_R4->SetLineStyle(2);
  mPP_R4->SetLineColor(2);
  mPP_R4->Draw("same L");
  
  uPP_R4_Bayes->Draw();
  for(int i = 0; i<nbins_cent; ++i){
    uPbPb_R4_Bayes[i]->Scale(ScaleFactor[i+1]);
    uPbPb_R4_Bayes[i]->SetMarkerStyle(20);
    uPbPb_R4_Bayes[i]->SetMarkerColor(i+3);
    uPbPb_R4_Bayes[i]->Draw("same");
    
    mPbPb_R4[i]->Scale(ScaleFactor[i]);
    mPbPb_R4[i]->SetLineStyle(2);
    mPbPb_R4[i]->SetLineColor(i+3);
    mPbPb_R4[i]->Draw("same L");
  }

  // also need to draw the MC as lines here: in all the 3 canvas. 
  
  cSpectra->SaveAs(Form("../../Plots/Final_paper_plots_spectra_%d.pdf",date.GetDate()),"RECREATE");
  

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
    
  tRAA->AddEntry(RAA_R2_Bayes[0],"R=0.2","pl");
  tRAA->AddEntry(RAA_R3_Bayes[0],"R=0.3","pl");
  tRAA->AddEntry(RAA_R4_Bayes[0],"R=0.4","pl");
  tRAA->SetTextSize(0.04);

  cRAA->cd(1);
  tRAA->Draw();
  cRAA->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets",algo,jet_type),0.2,0.23,16);
  //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
  cRAA->cd(2);
  drawText("Jet ID cut, |#eta|<2",0.1,0.3,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
  cRAA->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.3,16);
  drawText("Pile up rejection cut applied",0.1,0.2,16);

  cRAA->SaveAs(Form("../../Plots/Final_paper_plots_RAA_%d.pdf",date.GetDate()),"RECREATE");
    
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
    hRAA_R2_npart->SetBinContent(hRAA_R2_npart->FindBin(npart[i]), RAA_R2_Bayes[i]->GetBinContent(16));
    hRAA_R3_npart->SetBinContent(hRAA_R3_npart->FindBin(npart[i]), RAA_R3_Bayes[i]->GetBinContent(16));
    hRAA_R4_npart->SetBinContent(hRAA_R4_npart->FindBin(npart[i]), RAA_R4_Bayes[i]->GetBinContent(16));    
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
  cRAA_npart->SaveAs(Form("../../Plots/Final_paper_plots_RAA_npart_%d.pdf_",date.GetDate()),"RECREATE");
  
}
