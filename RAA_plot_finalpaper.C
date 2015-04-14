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
  
  // plot 1 - spectra plot showing pp and 6 different centrality classes PbPb spectra
  //        - have a 3 panel plot for the different radii, with each of them scaled by two orders of magnitude 

  Double_t ScaleFactor[nbins_cent+1] = {1, 1e2, 1e4, 1e6, 1e8, 1e10, 1e12};  

  TCanvas * cSpectra = new TCanvas("cSpectra","",1200,1000);
  makeMultiPanelCanvasWithGap(cSpectra,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  c

  // plot 2 - Bayesian unfolded RAA as a function of pT for the different radii
  //        - regular 6 panel plot 


    
  // plot 3 - RAA as a function of npart - taken from http://dde.web.cern.ch/dde/glauber_lhc.htm for 84 < pT < 114 in PbPb,PP
  //        - need to decide if we have to unfold this? or if we can just take that respective pt ranges from the already existing RAA histograms.  this is bin number 16 and 17 from the centrality classes weve measured.
  


}
