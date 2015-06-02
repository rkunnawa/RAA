// Raghav Kunnawalkam Elayavalli
// May 25 2015
// Rutgers
// raghav.k.e at CERN dot CH

//
// Macro to plot the ratio between Data and MC in the measured spectra before unfolding - to get them in the same units. 
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
#include <TColor.h>
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
#include "../Headers/utilities.h"


using namespace std;

void RAA_plot_datavsMC_check(int radius = 2,
			     chat * etaWidth = (chat*) "20_eta_20",
			     char *jet_type = "PF",
			     int unfoldingCut = 40)
{

  
  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  TFile *fin = TFile::Open(Form("Pawan_ntuple_PbPb_pp_calopfpt_ppNoJetidcut_R0p%d_SevilNewMBSubtraction_unfold_mcclosure_oppside_trgMC_noSmear_%s_%dGeVCut_ak%s_20150522.root",3,etaWidth,40,jet_type));

  TH1F * hData[nbins_cent], * hMC[nbins_cent], * hRatio[nbins_cent];
  
  
  

}
