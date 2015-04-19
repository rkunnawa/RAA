// Raghav Kunnawalkam Elayavalli
// Feb 5th 2014
// Rutgers

//
// Macro to read in the RAA/spectra from the PAS and make nice comparison plots to show in the HIN-W meeting on Feb 6th
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

/*
static const int nbins_eta = 2;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0},
  {-2.0,+2.0}
};

static const double delta_eta[nbins_eta] = {
  2.0,4.0
};

static const char etaWidth[nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20"
};
*/

//static const int no_radius = 2;//testing purposes 
//static const int list_radius[no_radius] = {3,4};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
static const int no_radius = 7; 
static const int list_radius[no_radius] = {1,2,3,4,5,6,7};
/*
// divide by bin width
void divideBinWidth(TH1 *h)
{
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    if(val!=0){
      val/=h->GetBinWidth(i);
      valErr/=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }  
  }
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}
*/

using namespace std;

void RAA_NEWcomparisonwithPAS(char *algo = "Pu", char *jet_type = "PF"){
  
  TStopwatch timer;
  timer.Start();

  cout<<"started the macro"<<endl;
  
  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};

  // load in the necessary files:
  TFile *f_unfold_R2 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_ppNoJetidcut_R0p2_unfold_mcclosure_oppside_trgMC_n20_eta_p20_30GeVCut_akPF_20150417.root");
  TFile *f_unfold_R3 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_ppNoJetidcut_R0p3_unfold_mcclosure_oppside_trgMC_n20_eta_p20_30GeVCut_akPF_20150417.root");
  TFile *f_unfold_R4 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_ppNoJetidcut_R0p4_unfold_mcclosure_oppside_trgMC_n20_eta_p20_30GeVCut_akPF_20150417.root");

  TFile *fPbPb_data_R2 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p2.root");
  TFile *fPP_data_R2 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pp_CutEfficiency_noJetID_exclusionhighertriggers_A_R0p2.root");

  TFile *fPbPb_data_R3 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p3.root");
  TFile *fPP_data_R3 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pp_CutEfficiency_noJetID_exclusionhighertriggers_A_R0p3.root");

  TFile *fPbPb_data_R4 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p4.root");
  TFile *fPP_data_R4 = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/Pp_CutEfficiency_noJetID_exclusionhighertriggers_A_R0p4.root");

  cout<<"loaded the files "<<endl;

  // get the latest histograms 
  // first array index corresponds to radius 
  TH1F *PbPb_bayesian[3][nbins_cent], *PbPb_binbybin[3][nbins_cent], *PbPb_measured[3][nbins_cent];
  TH1F *PbPb_measured_fine[3][nbins_cent];
  TH1F *PP_bayesian[3], *PP_binbybin[3], *PP_measured[3];
  TH1F *PP_measured_fine[3];
  TH1F *RAA_measured[3][nbins_cent];
  TH1F *RAA_binbybin[3][nbins_cent];
  TH1F *RAA_bayesian[3][nbins_cent];

  Double_t xAxis1[38] = {30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400};

  for(int i = 0;i<nbins_cent;i++){

    cout<<i<<endl;

    PbPb_bayesian[0][i] = (TH1F*)f_unfold_R2->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    PbPb_measured[0][i] = (TH1F*)f_unfold_R2->Get(Form("PbPb_measured_spectra_combined_cent%d",i));
    PbPb_binbybin[0][i] = (TH1F*)f_unfold_R2->Get(Form("PbPb_BinByBin_unfolded_spectra_combined_cent%d",i));
    PbPb_bayesian[1][i] = (TH1F*)f_unfold_R3->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    PbPb_measured[1][i] = (TH1F*)f_unfold_R3->Get(Form("PbPb_measured_spectra_combined_cent%d",i));
    PbPb_binbybin[1][i] = (TH1F*)f_unfold_R3->Get(Form("PbPb_BinByBin_unfolded_spectra_combined_cent%d",i));

    PbPb_bayesian[2][i] = (TH1F*)f_unfold_R4->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    PbPb_measured[2][i] = (TH1F*)f_unfold_R4->Get(Form("PbPb_measured_spectra_combined_cent%d",i));
    PbPb_binbybin[2][i] = (TH1F*)f_unfold_R4->Get(Form("PbPb_BinByBin_unfolded_spectra_combined_cent%d",i));

    // this can be HLT80 or HLTComb, check with both: 
    PbPb_measured_fine[0][i] = (TH1F*)fPbPb_data_R2->Get(Form("hpbpb_HLTComb_R2_n20_eta_p20_cent%d",i));
    PbPb_measured_fine[0][i] = (TH1F*)PbPb_measured_fine[0][i]->Rebin(37,Form("PbPb_measured_rebin_R2",i),xAxis1);
    //divideBinWidth(PbPb_measured_fine[0][i]);
    PbPb_measured_fine[1][i] = (TH1F*)fPbPb_data_R3->Get(Form("hpbpb_HLTComb_R3_n20_eta_p20_cent%d",i));
    PbPb_measured_fine[1][i] = (TH1F*)PbPb_measured_fine[1][i]->Rebin(37,Form("PbPb_measured_rebin_R3",i),xAxis1);
    //divideBinWidth(PbPb_measured_fine[1][i]);
    PbPb_measured_fine[2][i] = (TH1F*)fPbPb_data_R4->Get(Form("hpbpb_HLTComb_R4_n20_eta_p20_cent%d",i));
    PbPb_measured_fine[2][i] = (TH1F*)PbPb_measured_fine[2][i]->Rebin(37,Form("PbPb_measured_rebin_R4",i),xAxis1);
    //divideBinWidth(PbPb_measured_fine[2][i]);
    
    RAA_bayesian[0][i] = (TH1F*)f_unfold_R2->Get(Form("RAA_bayesian_cent%d",i));
    RAA_measured[0][i] = (TH1F*)f_unfold_R2->Get(Form("RAA_measured_cent%d",i));
    RAA_binbybin[0][i] = (TH1F*)f_unfold_R2->Get(Form("RAA_binbybin_cent%d",i));
    
    RAA_bayesian[1][i] = (TH1F*)f_unfold_R3->Get(Form("RAA_bayesian_cent%d",i));
    RAA_measured[1][i] = (TH1F*)f_unfold_R3->Get(Form("RAA_measured_cent%d",i));
    RAA_binbybin[1][i] = (TH1F*)f_unfold_R3->Get(Form("RAA_binbybin_cent%d",i));
    
    RAA_bayesian[2][i] = (TH1F*)f_unfold_R4->Get(Form("RAA_bayesian_cent%d",i));
    RAA_measured[2][i] = (TH1F*)f_unfold_R4->Get(Form("RAA_measured_cent%d",i));
    RAA_binbybin[2][i] = (TH1F*)f_unfold_R4->Get(Form("RAA_binbybin_cent%d",i));
    
    // for(int j = 1; j<=RAA_bayesian[0][i]->GetNbinsX(); ++j){
      
    //   cout<<j<<endl;
    //   RAA_bayesian[0][i]->SetBinError(j, RAA_measured[0][i]->GetBinError(j));
    //   RAA_bayesian[1][i]->SetBinError(j, RAA_measured[1][i]->GetBinError(j));
    //   RAA_bayesian[2][i]->SetBinError(j, RAA_measured[2][i]->GetBinError(j));
      
    // }
    
  }
  
  cout<<"loaded PbPb histograms "<<endl;
  
  PP_bayesian[0] = (TH1F*)f_unfold_R2->Get("PP_bayesian_unfolded_spectra");
  PP_measured[0] = (TH1F*)f_unfold_R2->Get("PP_measured_unfolded_spectra");
  PP_binbybin[0] = (TH1F*)f_unfold_R2->Get("PP_binbybin_unfolded_spectra");
  
  PP_bayesian[1] = (TH1F*)f_unfold_R3->Get("PP_bayesian_unfolded_spectra");
  PP_measured[1] = (TH1F*)f_unfold_R3->Get("PP_measured_unfolded_spectra");
  PP_binbybin[1] = (TH1F*)f_unfold_R3->Get("PP_binbybin_unfolded_spectra");
  
  PP_bayesian[2] = (TH1F*)f_unfold_R4->Get("PP_bayesian_unfolded_spectra");
  PP_measured[2] = (TH1F*)f_unfold_R4->Get("PP_measured_unfolded_spectra");
  PP_binbybin[2] = (TH1F*)f_unfold_R4->Get("PP_binbybin_unfolded_spectra");
  
  PP_measured_fine[0] = (TH1F*)fPP_data_R2->Get("hpp_HLTComb_R2_n20_eta_p20");
  PP_measured_fine[0] = (TH1F*)PP_measured_fine[0]->Rebin(37,"PP_measured_rebin_R2",xAxis1);
  //divideBinWidth(PP_measured_fine[0]);
  PP_measured_fine[1] = (TH1F*)fPP_data_R3->Get("hpp_HLTComb_R3_n20_eta_p20");
  PP_measured_fine[1] = (TH1F*)PP_measured_fine[1]->Rebin(37,"PP_measured_rebin_R3",xAxis1);
  //divideBinWidth(PP_measured_fine[1]);
  PP_measured_fine[2] = (TH1F*)fPP_data_R4->Get("hpp_HLTComb_R4_n20_eta_p20");
  PP_measured_fine[2] = (TH1F*)PP_measured_fine[2]->Rebin(37,"PP_measured_rebin_R4",xAxis1);
  //divideBinWidth(PP_measured_fine[2]);
  
  cout<<"loaded pp histograms "<<endl;
  
  
  //lets get the pas histograms

  TH1F *PASPbPb_measured[nbins_cent], *PASPP_measured;

  PASPbPb_measured[0] = new TH1F("PASPbPb_measured_cent0","",47, xAxis1);
  PASPbPb_measured[0]->SetBinContent(0,1.152652e+07);
  PASPbPb_measured[0]->SetBinContent(1,157758);
  PASPbPb_measured[0]->SetBinContent(2,36422);
  PASPbPb_measured[0]->SetBinContent(3,21272);
  PASPbPb_measured[0]->SetBinContent(4,23997);
  PASPbPb_measured[0]->SetBinContent(5,27940);
  PASPbPb_measured[0]->SetBinContent(6,24859);
  PASPbPb_measured[0]->SetBinContent(7,17870);
  PASPbPb_measured[0]->SetBinContent(8,11265);
  PASPbPb_measured[0]->SetBinContent(9,6837);
  PASPbPb_measured[0]->SetBinContent(10,4196);
  PASPbPb_measured[0]->SetBinContent(11,2621);
  PASPbPb_measured[0]->SetBinContent(12,1709);
  PASPbPb_measured[0]->SetBinContent(13,1160);
  PASPbPb_measured[0]->SetBinContent(14,762);
  PASPbPb_measured[0]->SetBinContent(15,511);
  PASPbPb_measured[0]->SetBinContent(16,337);
  PASPbPb_measured[0]->SetBinContent(17,266);
  PASPbPb_measured[0]->SetBinContent(18,180);
  PASPbPb_measured[0]->SetBinContent(19,119);
  PASPbPb_measured[0]->SetBinContent(20,101);
  PASPbPb_measured[0]->SetBinContent(21,69);
  PASPbPb_measured[0]->SetBinContent(22,46);
  PASPbPb_measured[0]->SetBinContent(23,34);
  PASPbPb_measured[0]->SetBinContent(24,33);
  PASPbPb_measured[0]->SetBinContent(25,26);
  PASPbPb_measured[0]->SetBinContent(26,17);
  PASPbPb_measured[0]->SetBinContent(27,14);
  PASPbPb_measured[0]->SetBinContent(28,7);
  PASPbPb_measured[0]->SetBinContent(29,8);
  PASPbPb_measured[0]->SetBinContent(30,7);
  PASPbPb_measured[0]->SetBinContent(31,1);
  PASPbPb_measured[0]->SetBinContent(32,3);
  PASPbPb_measured[0]->SetBinContent(33,3);
  PASPbPb_measured[0]->SetBinContent(34,3);
  PASPbPb_measured[0]->SetBinContent(35,1);
  PASPbPb_measured[0]->SetBinContent(36,1);
  PASPbPb_measured[0]->SetBinContent(38,1);
  PASPbPb_measured[0]->SetBinError(0,3395.073);
  PASPbPb_measured[0]->SetBinError(1,397.1876);
  PASPbPb_measured[0]->SetBinError(2,190.8455);
  PASPbPb_measured[0]->SetBinError(3,145.8492);
  PASPbPb_measured[0]->SetBinError(4,154.9097);
  PASPbPb_measured[0]->SetBinError(5,167.1526);
  PASPbPb_measured[0]->SetBinError(6,157.6674);
  PASPbPb_measured[0]->SetBinError(7,133.6787);
  PASPbPb_measured[0]->SetBinError(8,106.1367);
  PASPbPb_measured[0]->SetBinError(9,82.68615);
  PASPbPb_measured[0]->SetBinError(10,64.77654);
  PASPbPb_measured[0]->SetBinError(11,51.1957);
  PASPbPb_measured[0]->SetBinError(12,41.34005);
  PASPbPb_measured[0]->SetBinError(13,34.05877);
  PASPbPb_measured[0]->SetBinError(14,27.60435);
  PASPbPb_measured[0]->SetBinError(15,22.60531);
  PASPbPb_measured[0]->SetBinError(16,18.35756);
  PASPbPb_measured[0]->SetBinError(17,16.30951);
  PASPbPb_measured[0]->SetBinError(18,13.41641);
  PASPbPb_measured[0]->SetBinError(19,10.90871);
  PASPbPb_measured[0]->SetBinError(20,10.04988);
  PASPbPb_measured[0]->SetBinError(21,8.306624);
  PASPbPb_measured[0]->SetBinError(22,6.78233);
  PASPbPb_measured[0]->SetBinError(23,5.830952);
  PASPbPb_measured[0]->SetBinError(24,5.744563);
  PASPbPb_measured[0]->SetBinError(25,5.09902);
  PASPbPb_measured[0]->SetBinError(26,4.123106);
  PASPbPb_measured[0]->SetBinError(27,3.741657);
  PASPbPb_measured[0]->SetBinError(28,2.645751);
  PASPbPb_measured[0]->SetBinError(29,2.828427);
  PASPbPb_measured[0]->SetBinError(30,2.645751);
  PASPbPb_measured[0]->SetBinError(31,1);
  PASPbPb_measured[0]->SetBinError(32,1.732051);
  PASPbPb_measured[0]->SetBinError(33,1.732051);
  PASPbPb_measured[0]->SetBinError(34,1.732051);
  PASPbPb_measured[0]->SetBinError(35,1);
  PASPbPb_measured[0]->SetBinError(36,1);
  PASPbPb_measured[0]->SetBinError(38,1);
  PASPbPb_measured[0]->SetEntries(1.186697e+07);

  
  PASPbPb_measured[1] = new TH1F("PASPbPb_measured_cent1","",47, xAxis1);
  PASPbPb_measured[1]->SetBinContent(0,8384057);
  PASPbPb_measured[1]->SetBinContent(1,68266);
  PASPbPb_measured[1]->SetBinContent(2,19519);
  PASPbPb_measured[1]->SetBinContent(3,14704);
  PASPbPb_measured[1]->SetBinContent(4,18007);
  PASPbPb_measured[1]->SetBinContent(5,22146);
  PASPbPb_measured[1]->SetBinContent(6,20369);
  PASPbPb_measured[1]->SetBinContent(7,14461);
  PASPbPb_measured[1]->SetBinContent(8,8903);
  PASPbPb_measured[1]->SetBinContent(9,5530);
  PASPbPb_measured[1]->SetBinContent(10,3333);
  PASPbPb_measured[1]->SetBinContent(11,2159);
  PASPbPb_measured[1]->SetBinContent(12,1292);
  PASPbPb_measured[1]->SetBinContent(13,897);
  PASPbPb_measured[1]->SetBinContent(14,614);
  PASPbPb_measured[1]->SetBinContent(15,465);
  PASPbPb_measured[1]->SetBinContent(16,283);
  PASPbPb_measured[1]->SetBinContent(17,196);
  PASPbPb_measured[1]->SetBinContent(18,136);
  PASPbPb_measured[1]->SetBinContent(19,86);
  PASPbPb_measured[1]->SetBinContent(20,65);
  PASPbPb_measured[1]->SetBinContent(21,49);
  PASPbPb_measured[1]->SetBinContent(22,37);
  PASPbPb_measured[1]->SetBinContent(23,25);
  PASPbPb_measured[1]->SetBinContent(24,23);
  PASPbPb_measured[1]->SetBinContent(25,17);
  PASPbPb_measured[1]->SetBinContent(26,18);
  PASPbPb_measured[1]->SetBinContent(27,11);
  PASPbPb_measured[1]->SetBinContent(28,4);
  PASPbPb_measured[1]->SetBinContent(29,6);
  PASPbPb_measured[1]->SetBinContent(30,5);
  PASPbPb_measured[1]->SetBinContent(31,10);
  PASPbPb_measured[1]->SetBinContent(32,2);
  PASPbPb_measured[1]->SetBinContent(33,3);
  PASPbPb_measured[1]->SetBinContent(34,2);
  PASPbPb_measured[1]->SetBinContent(35,2);
  PASPbPb_measured[1]->SetBinContent(36,1);
  PASPbPb_measured[1]->SetBinContent(38,1);
  PASPbPb_measured[1]->SetBinError(0,2895.524);
  PASPbPb_measured[1]->SetBinError(1,261.2776);
  PASPbPb_measured[1]->SetBinError(2,139.7104);
  PASPbPb_measured[1]->SetBinError(3,121.2601);
  PASPbPb_measured[1]->SetBinError(4,134.1902);
  PASPbPb_measured[1]->SetBinError(5,148.8153);
  PASPbPb_measured[1]->SetBinError(6,142.72);
  PASPbPb_measured[1]->SetBinError(7,120.2539);
  PASPbPb_measured[1]->SetBinError(8,94.35571);
  PASPbPb_measured[1]->SetBinError(9,74.36397);
  PASPbPb_measured[1]->SetBinError(10,57.73214);
  PASPbPb_measured[1]->SetBinError(11,46.46504);
  PASPbPb_measured[1]->SetBinError(12,35.9444);
  PASPbPb_measured[1]->SetBinError(13,29.94996);
  PASPbPb_measured[1]->SetBinError(14,24.77902);
  PASPbPb_measured[1]->SetBinError(15,21.56386);
  PASPbPb_measured[1]->SetBinError(16,16.8226);
  PASPbPb_measured[1]->SetBinError(17,14);
  PASPbPb_measured[1]->SetBinError(18,11.6619);
  PASPbPb_measured[1]->SetBinError(19,9.273618);
  PASPbPb_measured[1]->SetBinError(20,8.062258);
  PASPbPb_measured[1]->SetBinError(21,7);
  PASPbPb_measured[1]->SetBinError(22,6.082763);
  PASPbPb_measured[1]->SetBinError(23,5);
  PASPbPb_measured[1]->SetBinError(24,4.795832);
  PASPbPb_measured[1]->SetBinError(25,4.123106);
  PASPbPb_measured[1]->SetBinError(26,4.242641);
  PASPbPb_measured[1]->SetBinError(27,3.316625);
  PASPbPb_measured[1]->SetBinError(28,2);
  PASPbPb_measured[1]->SetBinError(29,2.44949);
  PASPbPb_measured[1]->SetBinError(30,2.236068);
  PASPbPb_measured[1]->SetBinError(31,3.162278);
  PASPbPb_measured[1]->SetBinError(32,1.414214);
  PASPbPb_measured[1]->SetBinError(33,1.732051);
  PASPbPb_measured[1]->SetBinError(34,1.414214);
  PASPbPb_measured[1]->SetBinError(35,1.414214);
  PASPbPb_measured[1]->SetBinError(36,1);
  PASPbPb_measured[1]->SetBinError(38,1);
  PASPbPb_measured[1]->SetEntries(8585704);

  PASPbPb_measured[2] = new TH1F("PASPbPb_measured_cent2","",47, xAxis1);
  PASPbPb_measured[2]->SetBinContent(0,1.677722e+07);
  PASPbPb_measured[2]->SetBinContent(1,79071);
  PASPbPb_measured[2]->SetBinContent(2,35711);
  PASPbPb_measured[2]->SetBinContent(3,30933);
  PASPbPb_measured[2]->SetBinContent(4,37773);
  PASPbPb_measured[2]->SetBinContent(5,51709);
  PASPbPb_measured[2]->SetBinContent(6,51987);
  PASPbPb_measured[2]->SetBinContent(7,37244);
  PASPbPb_measured[2]->SetBinContent(8,23141);
  PASPbPb_measured[2]->SetBinContent(9,13787);
  PASPbPb_measured[2]->SetBinContent(10,8444);
  PASPbPb_measured[2]->SetBinContent(11,5412);
  PASPbPb_measured[2]->SetBinContent(12,3335);
  PASPbPb_measured[2]->SetBinContent(13,2344);
  PASPbPb_measured[2]->SetBinContent(14,1551);
  PASPbPb_measured[2]->SetBinContent(15,1019);
  PASPbPb_measured[2]->SetBinContent(16,716);
  PASPbPb_measured[2]->SetBinContent(17,535);
  PASPbPb_measured[2]->SetBinContent(18,352);
  PASPbPb_measured[2]->SetBinContent(19,229);
  PASPbPb_measured[2]->SetBinContent(20,188);
  PASPbPb_measured[2]->SetBinContent(21,144);
  PASPbPb_measured[2]->SetBinContent(22,126);
  PASPbPb_measured[2]->SetBinContent(23,79);
  PASPbPb_measured[2]->SetBinContent(24,48);
  PASPbPb_measured[2]->SetBinContent(25,26);
  PASPbPb_measured[2]->SetBinContent(26,36);
  PASPbPb_measured[2]->SetBinContent(27,32);
  PASPbPb_measured[2]->SetBinContent(28,19);
  PASPbPb_measured[2]->SetBinContent(29,20);
  PASPbPb_measured[2]->SetBinContent(30,15);
  PASPbPb_measured[2]->SetBinContent(31,8);
  PASPbPb_measured[2]->SetBinContent(32,4);
  PASPbPb_measured[2]->SetBinContent(33,3);
  PASPbPb_measured[2]->SetBinContent(34,2);
  PASPbPb_measured[2]->SetBinContent(35,6);
  PASPbPb_measured[2]->SetBinContent(36,2);
  PASPbPb_measured[2]->SetBinContent(37,1);
  PASPbPb_measured[2]->SetBinContent(38,2);
  PASPbPb_measured[2]->SetBinContent(39,1);
  PASPbPb_measured[2]->SetBinContent(41,1);
  PASPbPb_measured[2]->SetBinContent(44,2);
  PASPbPb_measured[2]->SetBinContent(46,1);
  PASPbPb_measured[2]->SetBinContent(48,5);
  PASPbPb_measured[2]->SetBinError(0,4144.642);
  PASPbPb_measured[2]->SetBinError(1,281.1957);
  PASPbPb_measured[2]->SetBinError(2,188.9735);
  PASPbPb_measured[2]->SetBinError(3,175.8778);
  PASPbPb_measured[2]->SetBinError(4,194.3528);
  PASPbPb_measured[2]->SetBinError(5,227.3961);
  PASPbPb_measured[2]->SetBinError(6,228.0066);
  PASPbPb_measured[2]->SetBinError(7,192.987);
  PASPbPb_measured[2]->SetBinError(8,152.1217);
  PASPbPb_measured[2]->SetBinError(9,117.4181);
  PASPbPb_measured[2]->SetBinError(10,91.89124);
  PASPbPb_measured[2]->SetBinError(11,73.5663);
  PASPbPb_measured[2]->SetBinError(12,57.74946);
  PASPbPb_measured[2]->SetBinError(13,48.41487);
  PASPbPb_measured[2]->SetBinError(14,39.38274);
  PASPbPb_measured[2]->SetBinError(15,31.92178);
  PASPbPb_measured[2]->SetBinError(16,26.75818);
  PASPbPb_measured[2]->SetBinError(17,23.13007);
  PASPbPb_measured[2]->SetBinError(18,18.76166);
  PASPbPb_measured[2]->SetBinError(19,15.13275);
  PASPbPb_measured[2]->SetBinError(20,13.71131);
  PASPbPb_measured[2]->SetBinError(21,12);
  PASPbPb_measured[2]->SetBinError(22,11.22497);
  PASPbPb_measured[2]->SetBinError(23,8.888194);
  PASPbPb_measured[2]->SetBinError(24,6.928203);
  PASPbPb_measured[2]->SetBinError(25,5.09902);
  PASPbPb_measured[2]->SetBinError(26,6);
  PASPbPb_measured[2]->SetBinError(27,5.656854);
  PASPbPb_measured[2]->SetBinError(28,4.358899);
  PASPbPb_measured[2]->SetBinError(29,4.472136);
  PASPbPb_measured[2]->SetBinError(30,3.872983);
  PASPbPb_measured[2]->SetBinError(31,2.828427);
  PASPbPb_measured[2]->SetBinError(32,2);
  PASPbPb_measured[2]->SetBinError(33,1.732051);
  PASPbPb_measured[2]->SetBinError(34,1.414214);
  PASPbPb_measured[2]->SetBinError(35,2.44949);
  PASPbPb_measured[2]->SetBinError(36,1.414214);
  PASPbPb_measured[2]->SetBinError(37,1);
  PASPbPb_measured[2]->SetBinError(38,1.414214);
  PASPbPb_measured[2]->SetBinError(39,1);
  PASPbPb_measured[2]->SetBinError(41,1);
  PASPbPb_measured[2]->SetBinError(44,1.414214);
  PASPbPb_measured[2]->SetBinError(46,1);
  PASPbPb_measured[2]->SetBinError(48,2.236068);
  PASPbPb_measured[2]->SetEntries(1.756412e+07);

  
  PASPbPb_measured[3] = new TH1F("PASPbPb_measured_cent3","",47, xAxis1);
  PASPbPb_measured[3]->SetBinContent(0,4552914);
  PASPbPb_measured[3]->SetBinContent(1,14329);
  PASPbPb_measured[3]->SetBinContent(2,10658);
  PASPbPb_measured[3]->SetBinContent(3,9834);
  PASPbPb_measured[3]->SetBinContent(4,10871);
  PASPbPb_measured[3]->SetBinContent(5,16906);
  PASPbPb_measured[3]->SetBinContent(6,19614);
  PASPbPb_measured[3]->SetBinContent(7,15027);
  PASPbPb_measured[3]->SetBinContent(8,9252);
  PASPbPb_measured[3]->SetBinContent(9,5603);
  PASPbPb_measured[3]->SetBinContent(10,3456);
  PASPbPb_measured[3]->SetBinContent(11,2105);
  PASPbPb_measured[3]->SetBinContent(12,1330);
  PASPbPb_measured[3]->SetBinContent(13,929);
  PASPbPb_measured[3]->SetBinContent(14,600);
  PASPbPb_measured[3]->SetBinContent(15,437);
  PASPbPb_measured[3]->SetBinContent(16,301);
  PASPbPb_measured[3]->SetBinContent(17,186);
  PASPbPb_measured[3]->SetBinContent(18,144);
  PASPbPb_measured[3]->SetBinContent(19,103);
  PASPbPb_measured[3]->SetBinContent(20,66);
  PASPbPb_measured[3]->SetBinContent(21,63);
  PASPbPb_measured[3]->SetBinContent(22,39);
  PASPbPb_measured[3]->SetBinContent(23,34);
  PASPbPb_measured[3]->SetBinContent(24,19);
  PASPbPb_measured[3]->SetBinContent(25,10);
  PASPbPb_measured[3]->SetBinContent(26,16);
  PASPbPb_measured[3]->SetBinContent(27,14);
  PASPbPb_measured[3]->SetBinContent(28,11);
  PASPbPb_measured[3]->SetBinContent(29,8);
  PASPbPb_measured[3]->SetBinContent(30,6);
  PASPbPb_measured[3]->SetBinContent(31,4);
  PASPbPb_measured[3]->SetBinContent(32,5);
  PASPbPb_measured[3]->SetBinContent(33,5);
  PASPbPb_measured[3]->SetBinContent(34,2);
  PASPbPb_measured[3]->SetBinContent(35,2);
  PASPbPb_measured[3]->SetBinContent(40,1);
  PASPbPb_measured[3]->SetBinContent(43,3);
  PASPbPb_measured[3]->SetBinContent(48,1);
  PASPbPb_measured[3]->SetBinError(0,2133.756);
  PASPbPb_measured[3]->SetBinError(1,119.7038);
  PASPbPb_measured[3]->SetBinError(2,103.2376);
  PASPbPb_measured[3]->SetBinError(3,99.16653);
  PASPbPb_measured[3]->SetBinError(4,104.2641);
  PASPbPb_measured[3]->SetBinError(5,130.0231);
  PASPbPb_measured[3]->SetBinError(6,140.05);
  PASPbPb_measured[3]->SetBinError(7,122.5847);
  PASPbPb_measured[3]->SetBinError(8,96.18732);
  PASPbPb_measured[3]->SetBinError(9,74.85319);
  PASPbPb_measured[3]->SetBinError(10,58.78775);
  PASPbPb_measured[3]->SetBinError(11,45.88028);
  PASPbPb_measured[3]->SetBinError(12,36.46917);
  PASPbPb_measured[3]->SetBinError(13,30.4795);
  PASPbPb_measured[3]->SetBinError(14,24.4949);
  PASPbPb_measured[3]->SetBinError(15,20.90454);
  PASPbPb_measured[3]->SetBinError(16,17.34935);
  PASPbPb_measured[3]->SetBinError(17,13.63818);
  PASPbPb_measured[3]->SetBinError(18,12);
  PASPbPb_measured[3]->SetBinError(19,10.14889);
  PASPbPb_measured[3]->SetBinError(20,8.124038);
  PASPbPb_measured[3]->SetBinError(21,7.937254);
  PASPbPb_measured[3]->SetBinError(22,6.244998);
  PASPbPb_measured[3]->SetBinError(23,5.830952);
  PASPbPb_measured[3]->SetBinError(24,4.358899);
  PASPbPb_measured[3]->SetBinError(25,3.162278);
  PASPbPb_measured[3]->SetBinError(26,4);
  PASPbPb_measured[3]->SetBinError(27,3.741657);
  PASPbPb_measured[3]->SetBinError(28,3.316625);
  PASPbPb_measured[3]->SetBinError(29,2.828427);
  PASPbPb_measured[3]->SetBinError(30,2.44949);
  PASPbPb_measured[3]->SetBinError(31,2);
  PASPbPb_measured[3]->SetBinError(32,2.236068);
  PASPbPb_measured[3]->SetBinError(33,2.236068);
  PASPbPb_measured[3]->SetBinError(34,1.414214);
  PASPbPb_measured[3]->SetBinError(35,1.414214);
  PASPbPb_measured[3]->SetBinError(40,1);
  PASPbPb_measured[3]->SetBinError(43,1.732051);
  PASPbPb_measured[3]->SetBinError(48,1);
  PASPbPb_measured[3]->SetEntries(4674908);

  PASPbPb_measured[4] = new TH1F("PASPbPb_measured_cent4","",47, xAxis1);
  PASPbPb_measured[4]->SetBinContent(0,733773);
  PASPbPb_measured[4]->SetBinContent(1,3350);
  PASPbPb_measured[4]->SetBinContent(2,2719);
  PASPbPb_measured[4]->SetBinContent(3,2564);
  PASPbPb_measured[4]->SetBinContent(4,2616);
  PASPbPb_measured[4]->SetBinContent(5,4024);
  PASPbPb_measured[4]->SetBinContent(6,5391);
  PASPbPb_measured[4]->SetBinContent(7,4199);
  PASPbPb_measured[4]->SetBinContent(8,2593);
  PASPbPb_measured[4]->SetBinContent(9,1548);
  PASPbPb_measured[4]->SetBinContent(10,951);
  PASPbPb_measured[4]->SetBinContent(11,627);
  PASPbPb_measured[4]->SetBinContent(12,423);
  PASPbPb_measured[4]->SetBinContent(13,230);
  PASPbPb_measured[4]->SetBinContent(14,164);
  PASPbPb_measured[4]->SetBinContent(15,108);
  PASPbPb_measured[4]->SetBinContent(16,77);
  PASPbPb_measured[4]->SetBinContent(17,71);
  PASPbPb_measured[4]->SetBinContent(18,37);
  PASPbPb_measured[4]->SetBinContent(19,17);
  PASPbPb_measured[4]->SetBinContent(20,31);
  PASPbPb_measured[4]->SetBinContent(21,18);
  PASPbPb_measured[4]->SetBinContent(22,7);
  PASPbPb_measured[4]->SetBinContent(23,11);
  PASPbPb_measured[4]->SetBinContent(24,5);
  PASPbPb_measured[4]->SetBinContent(25,5);
  PASPbPb_measured[4]->SetBinContent(27,4);
  PASPbPb_measured[4]->SetBinContent(28,2);
  PASPbPb_measured[4]->SetBinContent(29,1);
  PASPbPb_measured[4]->SetBinContent(30,1);
  PASPbPb_measured[4]->SetBinContent(32,3);
  PASPbPb_measured[4]->SetBinError(0,856.6055);
  PASPbPb_measured[4]->SetBinError(1,57.87918);
  PASPbPb_measured[4]->SetBinError(2,52.14403);
  PASPbPb_measured[4]->SetBinError(3,50.63596);
  PASPbPb_measured[4]->SetBinError(4,51.14685);
  PASPbPb_measured[4]->SetBinError(5,63.43501);
  PASPbPb_measured[4]->SetBinError(6,73.42343);
  PASPbPb_measured[4]->SetBinError(7,64.79969);
  PASPbPb_measured[4]->SetBinError(8,50.92151);
  PASPbPb_measured[4]->SetBinError(9,39.34463);
  PASPbPb_measured[4]->SetBinError(10,30.83829);
  PASPbPb_measured[4]->SetBinError(11,25.03997);
  PASPbPb_measured[4]->SetBinError(12,20.56696);
  PASPbPb_measured[4]->SetBinError(13,15.16575);
  PASPbPb_measured[4]->SetBinError(14,12.80625);
  PASPbPb_measured[4]->SetBinError(15,10.3923);
  PASPbPb_measured[4]->SetBinError(16,8.774964);
  PASPbPb_measured[4]->SetBinError(17,8.42615);
  PASPbPb_measured[4]->SetBinError(18,6.082763);
  PASPbPb_measured[4]->SetBinError(19,4.123106);
  PASPbPb_measured[4]->SetBinError(20,5.567764);
  PASPbPb_measured[4]->SetBinError(21,4.242641);
  PASPbPb_measured[4]->SetBinError(22,2.645751);
  PASPbPb_measured[4]->SetBinError(23,3.316625);
  PASPbPb_measured[4]->SetBinError(24,2.236068);
  PASPbPb_measured[4]->SetBinError(25,2.236068);
  PASPbPb_measured[4]->SetBinError(27,2);
  PASPbPb_measured[4]->SetBinError(28,1.414214);
  PASPbPb_measured[4]->SetBinError(29,1);
  PASPbPb_measured[4]->SetBinError(30,1);
  PASPbPb_measured[4]->SetBinError(32,1.732051);
  PASPbPb_measured[4]->SetEntries(765570);

  PASPbPb_measured[5] = new TH1F("PASPbPb_measured_cent5","",47, xAxis1);
  PASPbPb_measured[5]->SetBinContent(0,54371);
  PASPbPb_measured[5]->SetBinContent(1,560);
  PASPbPb_measured[5]->SetBinContent(2,422);
  PASPbPb_measured[5]->SetBinContent(3,458);
  PASPbPb_measured[5]->SetBinContent(4,440);
  PASPbPb_measured[5]->SetBinContent(5,649);
  PASPbPb_measured[5]->SetBinContent(6,924);
  PASPbPb_measured[5]->SetBinContent(7,744);
  PASPbPb_measured[5]->SetBinContent(8,510);
  PASPbPb_measured[5]->SetBinContent(9,314);
  PASPbPb_measured[5]->SetBinContent(10,168);
  PASPbPb_measured[5]->SetBinContent(11,101);
  PASPbPb_measured[5]->SetBinContent(12,67);
  PASPbPb_measured[5]->SetBinContent(13,63);
  PASPbPb_measured[5]->SetBinContent(14,33);
  PASPbPb_measured[5]->SetBinContent(15,26);
  PASPbPb_measured[5]->SetBinContent(16,10);
  PASPbPb_measured[5]->SetBinContent(17,12);
  PASPbPb_measured[5]->SetBinContent(18,9);
  PASPbPb_measured[5]->SetBinContent(19,5);
  PASPbPb_measured[5]->SetBinContent(20,6);
  PASPbPb_measured[5]->SetBinContent(21,1);
  PASPbPb_measured[5]->SetBinContent(22,4);
  PASPbPb_measured[5]->SetBinContent(23,1);
  PASPbPb_measured[5]->SetBinContent(24,2);
  PASPbPb_measured[5]->SetBinContent(27,1);
  PASPbPb_measured[5]->SetBinContent(30,2);
  PASPbPb_measured[5]->SetBinError(0,233.1759);
  PASPbPb_measured[5]->SetBinError(1,23.66432);
  PASPbPb_measured[5]->SetBinError(2,20.54264);
  PASPbPb_measured[5]->SetBinError(3,21.40093);
  PASPbPb_measured[5]->SetBinError(4,20.97618);
  PASPbPb_measured[5]->SetBinError(5,25.47548);
  PASPbPb_measured[5]->SetBinError(6,30.39737);
  PASPbPb_measured[5]->SetBinError(7,27.27636);
  PASPbPb_measured[5]->SetBinError(8,22.58318);
  PASPbPb_measured[5]->SetBinError(9,17.72005);
  PASPbPb_measured[5]->SetBinError(10,12.96148);
  PASPbPb_measured[5]->SetBinError(11,10.04988);
  PASPbPb_measured[5]->SetBinError(12,8.185353);
  PASPbPb_measured[5]->SetBinError(13,7.937254);
  PASPbPb_measured[5]->SetBinError(14,5.744563);
  PASPbPb_measured[5]->SetBinError(15,5.09902);
  PASPbPb_measured[5]->SetBinError(16,3.162278);
  PASPbPb_measured[5]->SetBinError(17,3.464102);
  PASPbPb_measured[5]->SetBinError(18,3);
  PASPbPb_measured[5]->SetBinError(19,2.236068);
  PASPbPb_measured[5]->SetBinError(20,2.44949);
  PASPbPb_measured[5]->SetBinError(21,1);
  PASPbPb_measured[5]->SetBinError(22,2);
  PASPbPb_measured[5]->SetBinError(23,1);
  PASPbPb_measured[5]->SetBinError(24,1.414214);
  PASPbPb_measured[5]->SetBinError(27,1);
  PASPbPb_measured[5]->SetBinError(30,1.414214);
  PASPbPb_measured[5]->SetEntries(59903);

  PASPP_measured = new TH1F("PASPP_measured","",47,xAxis1);
  PASPP_measured->SetBinContent(0,14639);
  PASPP_measured->SetBinContent(1,4909);
  PASPP_measured->SetBinContent(2,3996);
  PASPP_measured->SetBinContent(3,3742);
  PASPP_measured->SetBinContent(4,3820);
  PASPP_measured->SetBinContent(5,5043);
  PASPP_measured->SetBinContent(6,7599);
  PASPP_measured->SetBinContent(7,6703);
  PASPP_measured->SetBinContent(8,4182);
  PASPP_measured->SetBinContent(9,2566);
  PASPP_measured->SetBinContent(10,1518);
  PASPP_measured->SetBinContent(11,937);
  PASPP_measured->SetBinContent(12,609);
  PASPP_measured->SetBinContent(13,390);
  PASPP_measured->SetBinContent(14,270);
  PASPP_measured->SetBinContent(15,143);
  PASPP_measured->SetBinContent(16,139);
  PASPP_measured->SetBinContent(17,97);
  PASPP_measured->SetBinContent(18,72);
  PASPP_measured->SetBinContent(19,45);
  PASPP_measured->SetBinContent(20,39);
  PASPP_measured->SetBinContent(21,17);
  PASPP_measured->SetBinContent(22,11);
  PASPP_measured->SetBinContent(23,16);
  PASPP_measured->SetBinContent(24,7);
  PASPP_measured->SetBinContent(25,5);
  PASPP_measured->SetBinContent(26,6);
  PASPP_measured->SetBinContent(27,4);
  PASPP_measured->SetBinContent(28,2);
  PASPP_measured->SetBinContent(29,3);
  PASPP_measured->SetBinContent(32,1);
  PASPP_measured->SetBinContent(33,1);
  PASPP_measured->SetBinContent(35,2);
  PASPP_measured->SetBinContent(48,2);
  PASPP_measured->SetBinError(0,120.9917);
  PASPP_measured->SetBinError(1,70.06426);
  PASPP_measured->SetBinError(2,63.21392);
  PASPP_measured->SetBinError(3,61.17189);
  PASPP_measured->SetBinError(4,61.80615);
  PASPP_measured->SetBinError(5,71.01408);
  PASPP_measured->SetBinError(6,87.17224);
  PASPP_measured->SetBinError(7,81.87185);
  PASPP_measured->SetBinError(8,64.66838);
  PASPP_measured->SetBinError(9,50.6557);
  PASPP_measured->SetBinError(10,38.96152);
  PASPP_measured->SetBinError(11,30.61046);
  PASPP_measured->SetBinError(12,24.67793);
  PASPP_measured->SetBinError(13,19.74842);
  PASPP_measured->SetBinError(14,16.43168);
  PASPP_measured->SetBinError(15,11.95826);
  PASPP_measured->SetBinError(16,11.78983);
  PASPP_measured->SetBinError(17,9.848858);
  PASPP_measured->SetBinError(18,8.485281);
  PASPP_measured->SetBinError(19,6.708204);
  PASPP_measured->SetBinError(20,6.244998);
  PASPP_measured->SetBinError(21,4.123106);
  PASPP_measured->SetBinError(22,3.316625);
  PASPP_measured->SetBinError(23,4);
  PASPP_measured->SetBinError(24,2.645751);
  PASPP_measured->SetBinError(25,2.236068);
  PASPP_measured->SetBinError(26,2.44949);
  PASPP_measured->SetBinError(27,2);
  PASPP_measured->SetBinError(28,1.414214);
  PASPP_measured->SetBinError(29,1.732051);
  PASPP_measured->SetBinError(32,1);
  PASPP_measured->SetBinError(33,1);
  PASPP_measured->SetBinError(35,1.414214);
  PASPP_measured->SetBinError(48,1.414214);
  PASPP_measured->SetEntries(61535);

  
  TH1F *PASRAA_measured[3][nbins_cent];
  TH1F *PASRAA_bayesian[3][nbins_cent];
  TH1F *PASRAA_binbybin[3][nbins_cent];

  Double_t xAxis2[12] = {100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300}; 

  PASRAA_bayesian[1][5] = new TH1F("PASRAA_bayesian_R3_cent5","",11,xAxis2);
  PASRAA_bayesian[1][5]->SetBinContent(0,0.7225377);
  PASRAA_bayesian[1][5]->SetBinContent(1,0.8600517);
  PASRAA_bayesian[1][5]->SetBinContent(2,0.8100158);
  PASRAA_bayesian[1][5]->SetBinContent(3,0.749716);
  PASRAA_bayesian[1][5]->SetBinContent(4,0.728152);
  PASRAA_bayesian[1][5]->SetBinContent(5,0.8079244);
  PASRAA_bayesian[1][5]->SetBinContent(6,1.045846);
  PASRAA_bayesian[1][5]->SetBinContent(7,0.99492);
  PASRAA_bayesian[1][5]->SetBinContent(8,0.8501027);
  PASRAA_bayesian[1][5]->SetBinContent(9,0.7283186);
  PASRAA_bayesian[1][5]->SetBinContent(10,0.8501291);
  PASRAA_bayesian[1][5]->SetBinContent(11,0.8562201);
  PASRAA_bayesian[1][5]->SetBinContent(12,2.485746);
  PASRAA_bayesian[1][5]->SetBinError(0,0.006807087);
  PASRAA_bayesian[1][5]->SetBinError(1,0.01931266);
  PASRAA_bayesian[1][5]->SetBinError(2,0.02226945);
  PASRAA_bayesian[1][5]->SetBinError(3,0.02646403);
  PASRAA_bayesian[1][5]->SetBinError(4,0.03199702);
  PASRAA_bayesian[1][5]->SetBinError(5,0.04537526);
  PASRAA_bayesian[1][5]->SetBinError(6,0.06991723);
  PASRAA_bayesian[1][5]->SetBinError(7,0.07712956);
  PASRAA_bayesian[1][5]->SetBinError(8,0.07612922);
  PASRAA_bayesian[1][5]->SetBinError(9,0.06113506);
  PASRAA_bayesian[1][5]->SetBinError(10,0.07653623);
  PASRAA_bayesian[1][5]->SetBinError(11,0.1332746);
  PASRAA_bayesian[1][5]->SetBinError(12,0.7989296);
  PASRAA_bayesian[1][5]->SetMinimum(0);
  PASRAA_bayesian[1][5]->SetMaximum(2);
  PASRAA_bayesian[1][5]->SetEntries(1770.971);
  
  PASRAA_bayesian[1][4] = new TH1F("PASRAA_bayesian_R3_cent4","",11,xAxis2);
  PASRAA_bayesian[1][4]->SetBinContent(0,0.7299647);
  PASRAA_bayesian[1][4]->SetBinContent(1,0.7481229);
  PASRAA_bayesian[1][4]->SetBinContent(2,0.7609972);
  PASRAA_bayesian[1][4]->SetBinContent(3,0.7759888);
  PASRAA_bayesian[1][4]->SetBinContent(4,0.7765765);
  PASRAA_bayesian[1][4]->SetBinContent(5,0.7793122);
  PASRAA_bayesian[1][4]->SetBinContent(6,0.7750514);
  PASRAA_bayesian[1][4]->SetBinContent(7,0.8105651);
  PASRAA_bayesian[1][4]->SetBinContent(8,0.778167);
  PASRAA_bayesian[1][4]->SetBinContent(9,0.6663016);
  PASRAA_bayesian[1][4]->SetBinContent(10,0.736169);
  PASRAA_bayesian[1][4]->SetBinContent(11,0.8040183);
  PASRAA_bayesian[1][4]->SetBinContent(12,0.9631981);
  PASRAA_bayesian[1][4]->SetBinError(0,0.003463908);
  PASRAA_bayesian[1][4]->SetBinError(1,0.008645942);
  PASRAA_bayesian[1][4]->SetBinError(2,0.01082537);
  PASRAA_bayesian[1][4]->SetBinError(3,0.01412537);
  PASRAA_bayesian[1][4]->SetBinError(4,0.01806789);
  PASRAA_bayesian[1][4]->SetBinError(5,0.02291484);
  PASRAA_bayesian[1][4]->SetBinError(6,0.02692976);
  PASRAA_bayesian[1][4]->SetBinError(7,0.03325802);
  PASRAA_bayesian[1][4]->SetBinError(8,0.03975427);
  PASRAA_bayesian[1][4]->SetBinError(9,0.02987937);
  PASRAA_bayesian[1][4]->SetBinError(10,0.03618969);
  PASRAA_bayesian[1][4]->SetBinError(11,0.06938782);
  PASRAA_bayesian[1][4]->SetBinError(12,0.2031483);
  PASRAA_bayesian[1][4]->SetMinimum(0);
  PASRAA_bayesian[1][4]->SetMaximum(2);
  PASRAA_bayesian[1][4]->SetEntries(6061.544);

  PASRAA_bayesian[1][3] = new TH1F("PASRAA_bayesian_R3_cent3","",11,xAxis2);
  PASRAA_bayesian[1][3]->SetBinContent(0,0.7279794);
  PASRAA_bayesian[1][3]->SetBinContent(1,0.6681795);
  PASRAA_bayesian[1][3]->SetBinContent(2,0.6799733);
  PASRAA_bayesian[1][3]->SetBinContent(3,0.6649472);
  PASRAA_bayesian[1][3]->SetBinContent(4,0.6252716);
  PASRAA_bayesian[1][3]->SetBinContent(5,0.6117619);
  PASRAA_bayesian[1][3]->SetBinContent(6,0.6540294);
  PASRAA_bayesian[1][3]->SetBinContent(7,0.7157548);
  PASRAA_bayesian[1][3]->SetBinContent(8,0.6989863);
  PASRAA_bayesian[1][3]->SetBinContent(9,0.5692854);
  PASRAA_bayesian[1][3]->SetBinContent(10,0.6007913);
  PASRAA_bayesian[1][3]->SetBinContent(11,0.8094043);
  PASRAA_bayesian[1][3]->SetBinContent(12,2.021915);
  PASRAA_bayesian[1][3]->SetBinError(0,0.002593759);
  PASRAA_bayesian[1][3]->SetBinError(1,0.005876821);
  PASRAA_bayesian[1][3]->SetBinError(2,0.007441405);
  PASRAA_bayesian[1][3]->SetBinError(3,0.009385111);
  PASRAA_bayesian[1][3]->SetBinError(4,0.01101591);
  PASRAA_bayesian[1][3]->SetBinError(5,0.01384366);
  PASRAA_bayesian[1][3]->SetBinError(6,0.01738902);
  PASRAA_bayesian[1][3]->SetBinError(7,0.0228265);
  PASRAA_bayesian[1][3]->SetBinError(8,0.02662461);
  PASRAA_bayesian[1][3]->SetBinError(9,0.01938019);
  PASRAA_bayesian[1][3]->SetBinError(10,0.02266);
  PASRAA_bayesian[1][3]->SetBinError(11,0.05331068);
  PASRAA_bayesian[1][3]->SetBinError(12,0.2970788);
  PASRAA_bayesian[1][3]->SetMinimum(0);
  PASRAA_bayesian[1][3]->SetMaximum(2);
  PASRAA_bayesian[1][3]->SetEntries(9256.671);

  PASRAA_bayesian[1][2] = new TH1F("PASRAA_bayesian_R3_cent2","",11,xAxis2);
  PASRAA_bayesian[1][2]->SetBinContent(0,0.8238835);
  PASRAA_bayesian[1][2]->SetBinContent(1,0.5604139);
  PASRAA_bayesian[1][2]->SetBinContent(2,0.5365105);
  PASRAA_bayesian[1][2]->SetBinContent(3,0.5556849);
  PASRAA_bayesian[1][2]->SetBinContent(4,0.5274078);
  PASRAA_bayesian[1][2]->SetBinContent(5,0.5573406);
  PASRAA_bayesian[1][2]->SetBinContent(6,0.6040812);
  PASRAA_bayesian[1][2]->SetBinContent(7,0.6451674);
  PASRAA_bayesian[1][2]->SetBinContent(8,0.5992499);
  PASRAA_bayesian[1][2]->SetBinContent(9,0.4717564);
  PASRAA_bayesian[1][2]->SetBinContent(10,0.4812902);
  PASRAA_bayesian[1][2]->SetBinContent(11,0.6858578);
  PASRAA_bayesian[1][2]->SetBinContent(12,1.051401);
  PASRAA_bayesian[1][2]->SetBinError(0,0.002676996);
  PASRAA_bayesian[1][2]->SetBinError(1,0.004509155);
  PASRAA_bayesian[1][2]->SetBinError(2,0.005339063);
  PASRAA_bayesian[1][2]->SetBinError(3,0.007110624);
  PASRAA_bayesian[1][2]->SetBinError(4,0.008486022);
  PASRAA_bayesian[1][2]->SetBinError(5,0.01143331);
  PASRAA_bayesian[1][2]->SetBinError(6,0.01470297);
  PASRAA_bayesian[1][2]->SetBinError(7,0.01880827);
  PASRAA_bayesian[1][2]->SetBinError(8,0.02064378);
  PASRAA_bayesian[1][2]->SetBinError(9,0.01452445);
  PASRAA_bayesian[1][2]->SetBinError(10,0.01651545);
  PASRAA_bayesian[1][2]->SetBinError(11,0.04184596);
  PASRAA_bayesian[1][2]->SetBinError(12,0.1496494);
  PASRAA_bayesian[1][2]->SetMinimum(0);
  PASRAA_bayesian[1][2]->SetMaximum(2);
  PASRAA_bayesian[1][2]->SetEntries(10967.25);

  PASRAA_bayesian[1][1] = new TH1F("PASRAA_bayesian_R3_cent1","",11,xAxis2);
  PASRAA_bayesian[1][1]->SetBinContent(0,1.031868);
  PASRAA_bayesian[1][1]->SetBinContent(1,0.4557194);
  PASRAA_bayesian[1][1]->SetBinContent(2,0.480005);
  PASRAA_bayesian[1][1]->SetBinContent(3,0.4997296);
  PASRAA_bayesian[1][1]->SetBinContent(4,0.444041);
  PASRAA_bayesian[1][1]->SetBinContent(5,0.4247413);
  PASRAA_bayesian[1][1]->SetBinContent(6,0.4454679);
  PASRAA_bayesian[1][1]->SetBinContent(7,0.4785669);
  PASRAA_bayesian[1][1]->SetBinContent(8,0.4834603);
  PASRAA_bayesian[1][1]->SetBinContent(9,0.4274004);
  PASRAA_bayesian[1][1]->SetBinContent(10,0.4192492);
  PASRAA_bayesian[1][1]->SetBinContent(11,0.5767608);
  PASRAA_bayesian[1][1]->SetBinContent(12,1.223037);
  PASRAA_bayesian[1][1]->SetBinError(0,0.003550127);
  PASRAA_bayesian[1][1]->SetBinError(1,0.003977577);
  PASRAA_bayesian[1][1]->SetBinError(2,0.005214178);
  PASRAA_bayesian[1][1]->SetBinError(3,0.00704089);
  PASRAA_bayesian[1][1]->SetBinError(4,0.007806462);
  PASRAA_bayesian[1][1]->SetBinError(5,0.0093582);
  PASRAA_bayesian[1][1]->SetBinError(6,0.01187143);
  PASRAA_bayesian[1][1]->SetBinError(7,0.01518394);
  PASRAA_bayesian[1][1]->SetBinError(8,0.01841148);
  PASRAA_bayesian[1][1]->SetBinError(9,0.01456344);
  PASRAA_bayesian[1][1]->SetBinError(10,0.01581155);
  PASRAA_bayesian[1][1]->SetBinError(11,0.03811852);
  PASRAA_bayesian[1][1]->SetBinError(12,0.1827313);
  PASRAA_bayesian[1][1]->SetMinimum(0);
  PASRAA_bayesian[1][1]->SetMaximum(2);
  PASRAA_bayesian[1][1]->SetEntries(9198.666);

  PASRAA_bayesian[1][0] = new TH1F("PASRAA_bayesian_R3_cent0","",11,xAxis2);
  PASRAA_bayesian[1][0]->SetBinContent(0,1.418427);
  PASRAA_bayesian[1][0]->SetBinContent(1,0.4702078);
  PASRAA_bayesian[1][0]->SetBinContent(2,0.4736779);
  PASRAA_bayesian[1][0]->SetBinContent(3,0.4640111);
  PASRAA_bayesian[1][0]->SetBinContent(4,0.4361698);
  PASRAA_bayesian[1][0]->SetBinContent(5,0.4342382);
  PASRAA_bayesian[1][0]->SetBinContent(6,0.4725478);
  PASRAA_bayesian[1][0]->SetBinContent(7,0.5191602);
  PASRAA_bayesian[1][0]->SetBinContent(8,0.5139099);
  PASRAA_bayesian[1][0]->SetBinContent(9,0.4391448);
  PASRAA_bayesian[1][0]->SetBinContent(10,0.4229113);
  PASRAA_bayesian[1][0]->SetBinContent(11,0.5843745);
  PASRAA_bayesian[1][0]->SetBinContent(12,0.7831888);
  PASRAA_bayesian[1][0]->SetBinError(0,0.004730243);
  PASRAA_bayesian[1][0]->SetBinError(1,0.003983373);
  PASRAA_bayesian[1][0]->SetBinError(2,0.004953222);
  PASRAA_bayesian[1][0]->SetBinError(3,0.006293119);
  PASRAA_bayesian[1][0]->SetBinError(4,0.007391385);
  PASRAA_bayesian[1][0]->SetBinError(5,0.009421142);
  PASRAA_bayesian[1][0]->SetBinError(6,0.0120761);
  PASRAA_bayesian[1][0]->SetBinError(7,0.01596338);
  PASRAA_bayesian[1][0]->SetBinError(8,0.01862156);
  PASRAA_bayesian[1][0]->SetBinError(9,0.01424701);
  PASRAA_bayesian[1][0]->SetBinError(10,0.01533512);
  PASRAA_bayesian[1][0]->SetBinError(11,0.03767925);
  PASRAA_bayesian[1][0]->SetBinError(12,0.1174314);
  PASRAA_bayesian[1][0]->SetMinimum(0);
  PASRAA_bayesian[1][0]->SetMaximum(2);
  PASRAA_bayesian[1][0]->SetEntries(9671.122);

  //now that we have taken the necessary histogram, lets start to make plots:

  bool doRAA = true;
  bool doPbPbSpectra = false;
  bool doPPSpectra = false;

  if(doRAA){
    //plot1 - 6 panel plot of RAA, at R=0.3, in the old PAS and now.
    TCanvas *cRAA = new TCanvas("cRAA","RAA",1200,800);
    makeMultiPanelCanvasWithGap(cRAA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

    TLegend *tRAA = myLegend(0.15,0.75,0.85,0.9);
    TLine *lineRAA = new TLine(100,1,299,1);
    lineRAA->SetLineStyle(2);
    lineRAA->SetLineWidth(2);

    int ci;
    TBox *box;
    for(int i = 0;i<nbins_cent;i++){

      cRAA->cd(nbins_cent-i);

      RAA_bayesian[1][i]->SetMarkerColor(kBlack);
      RAA_bayesian[1][i]->SetMarkerStyle(20);
      makeHistTitle(RAA_bayesian[1][i],"","Jet p_{T} (GeV/c)","R_{AA}");
      RAA_bayesian[1][i]->SetAxisRange(60,299,"X");
      RAA_bayesian[1][i]->SetAxisRange(0,2,"Y");
      RAA_bayesian[1][i]->Draw("E0");

      PASRAA_bayesian[1][i]->SetMarkerColor(kBlack);
      PASRAA_bayesian[1][i]->SetMarkerStyle(24);
      PASRAA_bayesian[1][i]->Draw("same E0");

      //RAA_binbybin[1][i]->SetMarkerStyle(29);
      //RAA_binbybin[1][i]->SetMarkerColor(kBlue);
      //RAA_binbybin[1][i]->Draw("same");

      lineRAA->Draw();
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

      switch ( i ) {
	
      case 0:
	box = new TBox(100,0.3992887,110,0.5411268);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(110,0.4004558,120,0.5469);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(120,0.3797207,130,0.5483015);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(130,0.3460133,140,0.5263263);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(140,0.3479742,150,0.5205022);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(150,0.3878594,160,0.5572361);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(160,0.4375477,170,0.6007726);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(170,0.432544,180,0.5952758);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(180,0.3688692,200,0.5094203);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(200,0.3537775,240,0.492045);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(240,0.485418,300,0.6833309);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	
	break;

      case 1:
	box = new TBox(100,0.3922236,110,0.5192152);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(110,0.4114913,120,0.5485187);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(120,0.4142097,130,0.5852495);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(130,0.3565539,140,0.5315281);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(140,0.3448709,150,0.5046117);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(150,0.3711734,160,0.5197624);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(160,0.4105411,170,0.5465928);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(170,0.4144865,180,0.5524341);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(180,0.3660878,200,0.4887129);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(200,0.3584384,240,0.4800599);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(240,0.4915343,300,0.6619874);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	break;

      case 2:
	box = new TBox(100,0.4885437,110,0.6322841);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(110,0.4659843,120,0.6070368);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(120,0.4659718,130,0.6453979);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(130,0.4280813,140,0.6267344);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(140,0.4578278,150,0.6568533);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(150,0.5101032,160,0.6980592);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(160,0.5623943,170,0.7279404);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(170,0.5223327,180,0.6761671);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(180,0.4111619,200,0.5323508);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(200,0.4193837,240,0.5431968);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(240,0.5974073,300,0.7743082);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();

	break;

      case 3:

	box = new TBox(100,0.5888569,110,0.7475021);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(110,0.5972088,120,0.7627378);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(120,0.5630347,130,0.7668597);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(130,0.5120569,140,0.7384863);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(140,0.507423,150,0.7161008);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(150,0.5585642,160,0.7494945);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(160,0.6327194,170,0.7987902);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(170,0.6181627,180,0.7798098);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(180,0.5037794,200,0.6347914);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(200,0.532314,240,0.6692686);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(240,0.7185097,300,0.9002989);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	
	break;

      case 4:
	box = new TBox(100,0.6653466,110,0.8308992);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(110,0.6746653,120,0.8473291);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(120,0.662311,130,0.8896665);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(130,0.6405509,140,0.9126022);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(140,0.651443,150,0.9071814);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(150,0.6679778,160,0.882125);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(160,0.7248008,170,0.8963293);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(170,0.6962772,180,0.8600568);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(180,0.596691,200,0.7359122);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(200,0.6601083,240,0.8122298);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(240,0.7215662,300,0.8864704);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();

	break;

      case 5:
	box = new TBox(100,0.7705583,110,0.9495451);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(110,0.7227901,120,0.8972415);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(120,0.6428426,130,0.8565895);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(130,0.602727,140,0.853577);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(140,0.6775657,150,0.9382831);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(150,0.9043127,160,1.187379);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(160,0.8928045,170,1.097036);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(170,0.762919,180,0.9372863);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(180,0.6536744,200,0.8029628);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(200,0.762955,240,0.9373031);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
	box = new TBox(240,0.7678421,300,0.9445981);

	ci = TColor::GetColor("#cccccc");
	box->SetFillColor(ci);

	ci = TColor::GetColor("#cccccc");
	box->SetLineColor(ci);
	box->Draw();
		
	break;
	
      }// switch

      RAA_bayesian[1][i]->Draw("same E1");
      PASRAA_bayesian[1][i]->Draw("same E1");

    }// centrality loop

    tRAA->AddEntry(RAA_bayesian[1][0],"13-005, chMax/jtpt > 0.05,","pl");
    tRAA->AddEntry(PASRAA_bayesian[1][0],"12-004, trkMax/jtpt > 0.01","pl");
    tRAA->SetTextSize(0.04);

    cRAA->cd(1);
    tRAA->Draw();
    drawText("Bayesian Unfolding, 4 iterations",0.2,0.7,16);
    cRAA->cd(1);
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.3",algo,jet_type),0.2,0.23,16);
    //drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
    cRAA->cd(2);
    drawText("|#eta|<2",0.1,0.3,16);
    drawText("|vz|<15, HBHEfilter, pCES",0.1,0.2,16);
    //cRAA->cd(3);
    //drawText("Marguerite file, trigger HLT80",0.1,0.2,16);

    cRAA->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_newComparingwithPAS_ak%s3%s_%d.pdf",algo,jet_type,date.GetDate()),"RECREATE");
  }
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doPbPbSpectra){
    // plot 2 - spectra comparison: measured for PbPb - 6 panel plot.
    TCanvas *cPbPb_sigma = new TCanvas("cPbPb_sigma","PbPb inclusive jet invariant cross section",1200,800);
    makeMultiPanelCanvas(cPbPb_sigma,3,2,0.0,0.0,0.2,0.15,0.07); 

    for(int i = 0;i<nbins_cent;i++){

      cPbPb_sigma->cd(nbins_cent-i);
      cPbPb_sigma->cd(nbins_cent-i)->SetLogy();
      cPbPb_sigma->cd(nbins_cent-i)->SetLogx();

      makeHistTitle(PbPb_measured_fine[1][i]," ","Jet p_{T} (GeV/c)","arbitrary units");
      PbPb_measured_fine[1][i]->SetMarkerStyle(24);
      PbPb_measured_fine[1][i]->SetMarkerColor(kBlack);
      PbPb_measured_fine[1][i]->SetAxisRange(1e-1,1e6,"Y");
      PbPb_measured_fine[1][i]->SetAxisRange(30,500,"X");
      PbPb_measured_fine[1][i]->Draw();
    
      //PASPbPb_measured[i]->Scale(1./4);
      PASPbPb_measured[i]->SetMarkerStyle(20);
      PASPbPb_measured[i]->SetMarkerColor(kBlack);
      PASPbPb_measured[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.8,20);
    }

    cPbPb_sigma->cd(1);
    TLegend *PbPb_sigma = myLegend(0.25,0.6,0.5,0.9);
    PbPb_sigma->AddEntry(PbPb_measured[1][0],"13-005, chMax/jtpt  > 0.05","pl");
    PbPb_sigma->AddEntry(PASPbPb_measured[0],"12-004, trkMax/jtpt > 0.01","pl");
    PbPb_sigma->SetTextSize(0.04);
    PbPb_sigma->Draw();

    putCMSPrel();
    putPbPbLumi();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.3",algo,jet_type),0.3,0.2,16);
    //drawText("Spectra with NJet(p_{T}>50) < 7",0.55,0.55,16);

    cPbPb_sigma->cd(2);
    drawText("|#eta|<2",0.15,0.35,16);
    drawText("|vz|<15, pCES, HBHE",0.15,0.25,16);
    //drawText("hiNpix_1 > 38000 - 500*NJet",0.15,0.15,16);
    //cPbPb_sigma->cd(3);
    drawText("measured spectra, HLT80",0.15,0.15,16);

    cPbPb_sigma->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_spectra_newComparingwithPAS_ak%s3%s_%d.pdf",algo,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doPPSpectra){
    TCanvas *cPP_sigma = new TCanvas("cPP_sigma","PP inclusive jet invariant cross section",600,400);
    cPP_sigma->SetLogy();
    cPP_sigma->SetLogx();

    PP_measured_fine[1]->SetMarkerStyle(24);
    PP_measured_fine[1]->SetMarkerColor(kBlack);
    PP_measured_fine[1]->SetTitle(" ");
    PP_measured_fine[1]->SetXTitle("Jet p_{T} (GeV/c)");
    PP_measured_fine[1]->SetYTitle("arbitrary units");
    PP_measured_fine[1]->SetAxisRange(30,500,"X");
    PP_measured_fine[1]->Draw();

    PASPP_measured->SetMarkerColor(kBlack);
    PASPP_measured->SetMarkerStyle(20);
    PASPP_measured->Draw("same");

    TLegend *PP_sigma = myLegend(0.43,0.65,0.75,0.9);
    PP_sigma->AddEntry(PP_measured_fine[1]," latest Jet ID cut","pl");
    PP_sigma->AddEntry(PASPP_measured,"2012 PAS, trxMax/jtpt > 0.1","pl");
    PP_sigma->SetTextSize(0.04);
    PP_sigma->Draw();

    putCMSPrel();
    putPPLumi();
    drawText(Form("Anti-k_{T} %s Jets R=0.3",jet_type),0.15,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.15,0.33,16);
    drawText("measured spectra",0.15,0.16,16);

    cPP_sigma->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_spectra_newComparingwithPAS_ak3%s_%d.pdf",jet_type,date.GetDate()),"RECREATE");

  }
  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
