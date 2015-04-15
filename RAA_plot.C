// Raghav Kunnawalkam Elayavalli
// June 11 2014
// CERN

// RAA - plotting macro. pretty much all the neccessary plots will be made here in this macro. 
// will take input from the read macro and the analysis macro. Maybe not from the read macro - i will have to decide on that. 

// June 25th - just started working on the macro. will plot the iteration systematics, MC closure test, create setup for Unfolded vs measured, RAA and Normalized Response matrix plots as well. 

// July 1st - finished the macro. it would have been very easy and efficient to put them all in 2 segments, 
//            one for PbPb (with the centrality loop) and one for pp. 
//          - But ive decided to keep it this way since we can easily delete any segments which are complete and remake any plots individually. 

// Nov 14th - Fully active plotting macro with boolean variables to make certain plots alone. 

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

void RAA_plot(int radius = 2, char *algo = "Pu", char *jet_type = "PF", int unfoldingCut = 60){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  char *location = "MIT";

  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  
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
  
  //Boolean variable to make the small data histogram only output file. 
  bool makeRootFile = false;

  //Boolean Variables for the several plots: 
  bool doSpectra = false;

  bool doPbPbIterSys = false;
  bool doPPIterSys = false;
  bool doRAA = true;
  bool doPbPbMCClosure = true;
  bool doPPMCClosure = true;
  bool doPbPbDatavsMC = false;
  bool doPPDatavsMC = false;
  bool doPbPbNormRes = true;
  bool doPPNormRes = true;
  bool doPbPbsigma = true;
  bool doPPsigma = true;
  bool doGenSpectra = false;
  bool doPbPbTrgComb = true;
  bool doPbPb12003TrgComb = false;
  bool doPPTrgComb = true;
  bool doPPTrgContribution = false;

  bool doRandomCone = false;
  bool doSupernovaData = false;
  bool doSupernovaMC = false;
  bool doSupernovaDataMCUltraCentral = false;
  bool doSupernovaContributionData = false;
  bool doSupernovaContributionMC = false;
  bool doCentDatavsMC = false;
  bool doVertexDatavsMC = false;
  bool doPFElectronCheck = false;
  bool doJetVariablesCheck = false;
  bool doJetID = false;
  bool do2DCutpTplot = false;
  bool do2DCut1 = false;
  bool do2DCut2 = false;
  bool do2DCut3 = false;
  bool do2DCut4 = false;
  bool do2DCut5 = false;
  
  TFile *fin; 
  
  //if(location=="MIT") 
  fin= TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_calopfpt_jetidcut_R0p%d_unfold_mcclosure_oppside_trgMC_n20_eta_p20_%dGeVCut_ak%s_20150415.root",radius,unfoldingCut,jet_type));
  //fin= TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_data_ak%s%s_testComb4_cut1_20141111.root",algo,jet_type));
  //if(location=="CERN")fin= TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_unfo_ak%s%d%s_20140911.root",algo,radius,jet_type));
  //if(location=="MPB") fin= TFile::Open(Form(""))

  
  TH1F *dPbPb_TrgComb[nbins_cent+1], *dPbPb_Comb[nbins_cent+1], *dPbPb_Trg80[nbins_cent+1], *dPbPb_Trg65[nbins_cent+1], *dPbPb_Trg55[nbins_cent+1], *dPbPb_1[nbins_cent+1], *dPbPb_2[nbins_cent+1], *dPbPb_3[nbins_cent+1], *dPbPb_80[nbins_cent+1], *dPbPb_65[nbins_cent+1], *dPbPb_55[nbins_cent+1], *dPbPb_MinBias[nbins_cent+1];
  
  TH1F *mPbPb_Gen[nbins_cent+1], *mPbPb_Reco[nbins_cent+1];
  TH2F *mPbPb_Matrix[nbins_cent+1], *mPbPb_Response[nbins_cent+1], *mPbPb_ResponseNorm[nbins_cent+1];
  TH1F *mPbPb_mcclosure_data[nbins_cent+1];
  TH1F *mPbPb_mcclosure_gen[nbins_cent+1];
  TH1F *dPP_1, *dPP_2, *dPP_3, *dPP_Comb;
  
  TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix, *mPP_Response;
  TH2F *mPP_ResponseNorm;
  TH1F *mPP_mcclosure_data;
  TH1F *mPP_mcclosure_gen;
  
  const int Iterations = 10; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1], *uPbPb_BinByBin[nbins_cent+1]; 
  TH1F *uPbPb_BayesianIter[nbins_cent+1][Iterations];
  TH1F *uPbPb_MC_Bayes[nbins_cent];
  TH1F *uPbPb_MC_BayesianIter[nbins_cent+1][Iterations];
  TH1F *uPbPb_MC_BinByBin[nbins_cent+1];
  
  TH1F *uPP_Bayes, *uPP_BinByBin;
  TH1F *uPP_BayesianIter[Iterations];
  TH1F *uPP_MC_Bayes, *uPP_MC_BinByBin;
  TH1F *uPP_MC_BayesianIter[Iterations];
  
  TH1F *RAA_measured[nbins_cent+1];
  TH1F *RAA_binbybin[nbins_cent+1];
  TH1F *RAA_bayesian[nbins_cent+1];
  
  for(int i = 0;i<nbins_cent;++i){
    
    cout<<"cent = "<<i<<endl;
    dPbPb_TrgComb[i] = (TH1F*)fin->Get(Form("PbPb_measured_spectra_combined_cent%d",i));
    dPbPb_TrgComb[i]->Print("base");
    dPbPb_Trg80[i] = (TH1F*)fin->Get(Form("PbPb_measured_spectra_jet80_cent%d",i));
    dPbPb_Trg65[i] = (TH1F*)fin->Get(Form("PbPb_measured_spectra_jet65_cent%d",i));
    dPbPb_Trg55[i] = (TH1F*)fin->Get(Form("PbPb_measured_spectra_jet55_cent%d",i));    
    dPbPb_MinBias[i] = (TH1F*)fin->Get(Form("PbPb_measured_spectra_MinBias_cent%d",i));    
    
    mPbPb_Reco[i] = (TH1F*)fin->Get(Form("hpbpb_reco_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Gen[i] = (TH1F*)fin->Get(Form("hpbpb_gen_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_ResponseNorm[i] = (TH2F*)fin->Get(Form("PbPb_normalized_response_matrix_cent%d",i));
    mPbPb_Response[i] = (TH2F*)fin->Get(Form("mPbPb_Response_cent%d",i));
    uPbPb_Bayes[i] = (TH1F*)fin->Get(Form("PbPb_bayesian_unfolded_spectra_combined_cent%d",i));
    uPbPb_MC_Bayes[i] = (TH1F*)fin->Get(Form("uPbPb_MC_Bayes_cent%d",i));
    uPbPb_BinByBin[i] = (TH1F*)fin->Get(Form("PbPb_BinByBin_unfolded_spectra_combined_cent%d",i));
    uPbPb_MC_BinByBin[i] = (TH1F*)fin->Get(Form("uPbPb_MC_BinByBin_cent%d",i));
    RAA_bayesian[i] = (TH1F*)fin->Get(Form("RAA_bayesian_cent%d",i));
    RAA_binbybin[i] = (TH1F*)fin->Get(Form("RAA_binbybin_cent%d",i));
    RAA_measured[i] = (TH1F*)fin->Get(Form("RAA_measured_cent%d",i));
    mPbPb_mcclosure_data[i] = (TH1F*)fin->Get(Form("mPbPb_mclosure_data_cent%d",radius,i));
    mPbPb_mcclosure_gen[i] = (TH1F*)fin->Get(Form("hpbpb_mcclosure_gen_JetComb_R%d_n20_eta_p20_cent%d",radius,i));
    
    for(int j = 0;j<Iterations;j++){
      uPbPb_BayesianIter[i][j] = (TH1F*)fin->Get(Form("uPbPb_BayesianIter%d_cent%d",j+1,i));
      //uPbPb_BayesianIter[j][i]->Rebin(10);
      //uPbPb_BayesianIter[j][i]->Scale(1./10);
      uPbPb_MC_BayesianIter[i][j] = (TH1F*)fin->Get(Form("uPbPb_MC_BayesianIter%d_cent%d",j+1,i));
    }
  }
  
  dPP_Comb = (TH1F*)fin->Get("pp_measured_spectra_combined");
  dPP_Comb->Print("base");
  dPP_1 = (TH1F*)fin->Get("pp_measured_spectra_jet80");
  dPP_1->Print("base");
  dPP_2 = (TH1F*)fin->Get("pp_measured_spectra_jet60");
  dPP_2->Print("base");
  dPP_3 = (TH1F*)fin->Get("pp_measured_spectra_jet40");
  dPP_3->Print("base");
  
  mPP_ResponseNorm = (TH2F*)fin->Get("PP_normalized_response_matrix");	
  mPP_mcclosure_data = (TH1F*)fin->Get("mPP_mcclosure_data");
  mPP_mcclosure_gen = (TH1F*)fin->Get("mPP_mcclosure_gen");
  
  uPP_Bayes = (TH1F*)fin->Get("PP_bayesian_unfolded_spectra");
  uPP_Bayes->Print("base");
  uPP_BinByBin = (TH1F*)fin->Get("PP_BinByBin_unfolded_spectra");
  uPP_MC_Bayes = (TH1F*)fin->Get("PP_MC_Bayes_unfolded_spectra");
  uPP_MC_BinByBin = (TH1F*)fin->Get("PP_MC_BinByBin_unfolded_spectra");
  //mPP_Gen = (TH1F*)fin->Get(Form("hpp_gen_R%d_n20_eta_p20",radius));
  //mPP_Reco = (TH1F*)fin->Get(Form("hpp_reco_R%d_n20_eta_p20",radius));
  
  for(int i = 0;i<Iterations;i++){

    uPP_BayesianIter[i] = (TH1F*)fin->Get(Form("uPP_BayesianIter%d",i+1));
    //uPP_BayesianIter[i]->Rebin(10);
    //uPP_BayesianIter[i]->Scale(1./10);
    uPP_MC_BayesianIter[i] = (TH1F*)fin->Get(Form("uPP_MC_BayesianIter%d",i+1));
    
  }
  
  
  TFile * fPP_in = TFile::Open(Form("../../Output/Pp_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));
  TFile * fPP_in_trigger = TFile::Open(Form("../../Output/Pp_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d_for_trigger.root",radius));
  fPP_in->Print();
  
  TH1F * dPPComb = (TH1F*)fPP_in->Get(Form("hpp_HLTComb_R%d_n20_eta_p20",radius)); 
  
  TH1F * hPP_80Cont = (TH1F*)fPP_in_trigger->Get(Form("hpp_HLT80_R%d_n20_eta_p20",radius));
  hPP_80Cont->Divide(dPPComb);
  TH1F * hPP_60Cont = (TH1F*)fPP_in_trigger->Get(Form("hpp_HLT60_R%d_n20_eta_p20",radius));
  hPP_60Cont->Divide(dPPComb);
  TH1F * hPP_40Cont = (TH1F*)fPP_in_trigger->Get(Form("hpp_HLT40_R%d_n20_eta_p20",radius));
  hPP_40Cont->Divide(dPPComb);
  
  hPP_80Cont = (TH1F*)hPP_80Cont->Rebin(nbins_pt, "hPP_80", boundaries_pt);
  divideBinWidth(hPP_80Cont);
  hPP_60Cont = (TH1F*)hPP_60Cont->Rebin(nbins_pt, "hPP_60", boundaries_pt);
  divideBinWidth(hPP_60Cont);
  hPP_40Cont = (TH1F*)hPP_40Cont->Rebin(nbins_pt, "hPP_40", boundaries_pt);
  divideBinWidth(hPP_40Cont);  

  /*
  // get histograms from the MC file. 
  TFile *fMCin = TFile::Open(Form("/Users/keraghav/WORK/RAA/Output/PbPb_mc_nocut_ak%s%s_20141210.root",algo,jet_type));
  TFile *fDatain = TFile::Open(Form("/Users/keraghav/WORK/RAA/Output/PbPb_jetntuple_withEvtCuts_SuperNovaRejected_ak%s%d%s_20141209.root",algo,radius,jet_type));
  
  TH1F *hPbPb_MC_jtpu[no_radius][nbins_eta][nbins_cent+1];
  
  TH1F *hPbPb12003Jet80[nbins_cent+1], *hPbPb12003Jet65[nbins_cent+1], *hPbPb12003Jet55[nbins_cent+1], *hPbPb12003Comb[nbins_cent+1];
  
  TH1F *hpbpb_TrgObj80_nJet_g7[nbins_cent+1];
  TH1F *hpbpb_TrgObj65_nJet_g7[nbins_cent+1];
  TH1F *hpbpb_TrgObj55_nJet_g7[nbins_cent+1];
  TH1F *hpbpb_TrgObjMB_nJet_g7[nbins_cent+1];
  TH1F *hpbpb_TrgObjComb_nJet_g7[nbins_cent+1];
  
  TH1F *hpbpb_TrgObj80_nJet_l7[nbins_cent+1];
  TH1F *hpbpb_TrgObj65_nJet_l7[nbins_cent+1];
  TH1F *hpbpb_TrgObj55_nJet_l7[nbins_cent+1];
  TH1F *hpbpb_TrgObjMB_nJet_l7[nbins_cent+1];
  TH1F *hpbpb_TrgObjComb_nJet_l7[nbins_cent+1];
  
  // histograms for the supernova cut rejection 
  TH2F *hpbpb_Npix_before_cut_data[nbins_cent+2]; 
  TH2F *hpbpb_Npix_after_cut_data[nbins_cent+1]; 
  TH2F *hpbpb_Npix_before_cut_mc[nbins_cent+2]; 
  TH2F *hpbpb_Npix_after_cut_mc[nbins_cent+1]; 
  
  // histograms to check the contribution of the nJets>7 (with pT>50 and |eta|<2) to the nJets < 7 
  TH1F *hpbpb_pt_Njet_g7[nbins_cent+1];
  TH1F *hpbpb_pt_Njet_l7[nbins_cent+1];
  
  // hpbpb_Npix_before_cut_data[nbins_cent+1] = (TH2F*)fDatain->Get(Form("hpbpb_Npix_before_cut_R%d_n20_eta_p20_cent%d",radius,nbins_cent+1));
  // hpbpb_Npix_before_cut_data[nbins_cent+1]->Print("base");
  
  // hpbpb_Npix_before_cut_mc[nbins_cent+1] = (TH2F*)fMCin->Get(Form("hpbpb_Npix_before_cut_R%d_n20_eta_p20_cent%d",radius,nbins_cent+1));
  // hpbpb_Npix_before_cut_mc[nbins_cent+1]->Print("base");
  
  // I changed it from <= to < due to a problem with the data macro which is now fixed for the hpbpb_JetComb_oldR%d_n20_eta_p20_cent%d should be hpbpb_JetComb_old_R%d_n20_eta_p20_cent%d
  /*
  for(int i = 0;i<nbins_cent;i++){

    cout<<"cent = "<<i<<endl;
    dPbPb_TrgComb[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObjComb_R%d_n20_eta_p20_cent%d",radius,i));
    dPbPb_TrgComb[i]->Print("base");
    dPbPb_Trg80[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj80_R%d_n20_eta_p20_cent%d",radius,i));
    dPbPb_Trg65[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj65_R%d_n20_eta_p20_cent%d",radius,i));
    dPbPb_Trg55[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj55_R%d_n20_eta_p20_cent%d",radius,i));   

    hPbPb12003Jet80[i] = (TH1F*)fDatain->Get(Form("hpbpb_Jet80_R%d_n20_eta_p20_cent%d",radius,i));
    hPbPb12003Jet80[i]->Print("base");
    hPbPb12003Jet65[i] = (TH1F*)fDatain->Get(Form("hpbpb_Jet65_noJet80_R%d_n20_eta_p20_cent%d",radius,i));
    hPbPb12003Jet65[i]->Print("base");
    hPbPb12003Jet55[i] = (TH1F*)fDatain->Get(Form("hpbpb_Jet55_noJet65_80_R%d_n20_eta_p20_cent%d",radius,i));
    hPbPb12003Jet55[i]->Print("base");
    hPbPb12003Comb[i] = (TH1F*)fDatain->Get(Form("hpbpb_JetComb_oldR%d_n20_eta_p20_cent%d",radius,i));
    hPbPb12003Comb[i]->Print("base");

    hpbpb_Npix_before_cut_data[i] = (TH2F*)fDatain->Get(Form("hpbpb_Npix_before_cut_R%d_n20_eta_p20_cent%d",radius,i));
    hpbpb_Npix_after_cut_data[i] = (TH2F*)fDatain->Get(Form("hpbpb_Npix_after_cut_R%d_n20_eta_p20_cent%d",radius,i));
    hpbpb_Npix_before_cut_data[i]->Print("base");
    hpbpb_Npix_after_cut_data[i]->Print("base");
    
    hpbpb_Npix_before_cut_mc[i] = (TH2F*)fMCin->Get(Form("hpbpb_Npix_before_cut_R%d_n20_eta_p20_cent%d",radius,i));
    hpbpb_Npix_before_cut_mc[i]->Print("base");
    hpbpb_Npix_after_cut_mc[i] = (TH2F*)fMCin->Get(Form("hpbpb_Npix_after_cut_R%d_n20_eta_p20_cent%d",radius,i));
    hpbpb_Npix_after_cut_mc[i]->Print("base");
    
    hpbpb_pt_Njet_g7[i] = (TH1F*)fMCin->Get(Form("hpbpb_pt_Njet_g7_R%d_n20_eta_p20_cent%d",radius,i));
    hpbpb_pt_Njet_g7[i]->Print("base");
    hpbpb_pt_Njet_l7[i] = (TH1F*)fMCin->Get(Form("hpbpb_pt_Njet_l7_R%d_n20_eta_p20_cent%d",radius,i));
    hpbpb_pt_Njet_l7[i]->Print("base");
    
    
    // hpbpb_TrgObj80_nJet_g7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj80_nJet_g7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObj80_nJet_g7[i]->Print("base");
    // hpbpb_TrgObj65_nJet_g7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj65_nJet_g7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObj65_nJet_g7[i]->Print("base");            
    // hpbpb_TrgObj55_nJet_g7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj55_nJet_g7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObj55_nJet_g7[i]->Print("base");
    // hpbpb_TrgObjComb_nJet_g7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObjComb_nJet_g7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObjComb_nJet_g7[i]->Print("base");

    // hpbpb_TrgObj80_nJet_l7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj80_nJet_l7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObj80_nJet_l7[i]->Print("base");
    // hpbpb_TrgObj65_nJet_l7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj65_nJet_l7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObj65_nJet_l7[i]->Print("base");
    // hpbpb_TrgObj55_nJet_l7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObj55_nJet_l7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObj55_nJet_l7[i]->Print("base");
    // hpbpb_TrgObjComb_nJet_l7[i] = (TH1F*)fDatain->Get(Form("hpbpb_TrgObjComb_nJet_l7_R%d_n20_eta_p20_cent%d",radius,i));
    // hpbpb_TrgObjComb_nJet_l7[i]->Print("base");
    
    
    for(int k = 0;k<no_radius;k++){
      
      for(int j = 0;j<nbins_eta;j++){
	
	//hPbPb_MC_jtpu[k][j][i] = (TH1F*)fMCin->Get(Form("hpbpb_jtpu_noScale_R%d_%s_cent%d",list_radius[k],etaWidth[j],i));
	
      }//eta bins loop  
    }//radius loop
    
  }//centrality loop

  */
  
  if(makeRootFile){
    //intermediate (small) output file for running the unfolding: 
    TFile fout(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_data_ak%s%s_testComb4_cut1_%d.root",algo,jet_type,date.GetDate()),"RECREATE");
    fout.cd();
    for(int i = 0;i<=nbins_cent;i++){
      dPbPb_TrgComb[i]->Write();
      dPbPb_Trg80[i]->Write();
      dPbPb_Trg65[i]->Write();
      dPbPb_Trg55[i]->Write();
    }
    fout.Write();
    fout.Close();
  }

  
  //Ok now that we have loaded all the histograms we need - lets start making the plots 
  Double_t scaleFactor[nbins_cent+1] = {1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6};


  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  if(doSpectra){  
    // plot 0 - PbPb and pp Unfolded Spectra compared with Data(measured) and MC Spectra. in a 2 by 1 canvas. with the PbPb doing the multiply by powers of 10. unfolded in black circles, MC in black (dotted) line, measured in red open boxes. 
  
    TCanvas *cSpectra = new TCanvas("cSpectra","PbPb and PP spectra",1000,1000);
    cSpectra->Divide(2,1);
  
    cSpectra->cd(1);
    cSpectra->cd(1)->SetLogy();
    cSpectra->cd(1)->SetLogx();
    //cSpectra->Set
  
    TLegend *PbPbSpectra = myLegend(0.4,0.65,0.85,0.9);
  
    for(int i = 0;i<=nbins_cent;i++){
    
      uPbPb_Bayes[i]->Scale(1./(scaleFactor[i]));
      dPbPb_TrgComb[i]->Scale(1./(scaleFactor[i]));
      mPbPb_Reco[i]->Scale(1./(scaleFactor[i]));

      uPbPb_Bayes[i]->SetMarkerStyle(20);
      uPbPb_Bayes[i]->SetMarkerColor(i+1);
      if(i==0){
	uPbPb_Bayes[i]->SetAxisRange(1e-15,1e5,"Y");
	uPbPb_Bayes[i]->SetAxisRange(30,600,"X");
	makeHistTitle(uPbPb_Bayes[i],"","Jet p_{T} (GeV/c)","arbitrary for now");
	uPbPb_Bayes[i]->Draw();
      }else 
	uPbPb_Bayes[i]->Draw("same");

      dPbPb_TrgComb[i]->SetMarkerStyle(25);
      dPbPb_TrgComb[i]->SetMarkerColor(i+1);
      dPbPb_TrgComb[i]->Draw("same");

      mPbPb_Reco[i]->SetMarkerStyle(27);
      mPbPb_Reco[i]->SetMarkerColor(i+1);
      mPbPb_Reco[i]->Draw("same");

    }

    for(int i = nbins_cent;i>=0;i--){
      if(i<nbins_cent)
	PbPbSpectra->AddEntry(uPbPb_Bayes[i],Form("%2.0f - %2.0f cent * 10^{%d}",5*boundaries_cent[i],5*boundaries_cent[i+1],i),"pl");
      else 
	PbPbSpectra->AddEntry(uPbPb_Bayes[i],Form("0 - 200 cent x 10^{%d} x 10^{6}",i),"pl");
    }

    PbPbSpectra->SetTextSize(0.04);
    PbPbSpectra->Draw();
    putCMSPrel();
    drawText("PbPb",0.2,0.8,20);
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.5,0.9,20);
    drawText("#frac{chMax}{jtpt}>0.01",0.15,0.2,20);

    cSpectra->cd(2);						
    cSpectra->cd(2)->SetLogy();
    cSpectra->cd(2)->SetLogx();
    makeHistTitle(uPP_Bayes,"","Jet p_{T} (GeV/c)","arbitrary for now");
    uPP_Bayes->SetMarkerStyle(20);
    uPP_Bayes->SetMarkerColor(kBlack);
    uPP_Bayes->SetAxisRange(30,600,"X");
    uPP_Bayes->Draw();

    drawText("#frac{chMax}{jtpt}>0.01",0.15,0.2,20);
    drawText(Form("Anti-k_{T} %s Jets R=0.%d",jet_type,radius),0.45,0.9,20);
    drawText("2.76 TeV, |#eta|<2, |vz|<15",0,0.9,20);
    drawText("pp",0.3,0.8,20);
    
    dPP_Comb->SetMarkerStyle(25);
    dPP_Comb->SetMarkerColor(kBlack);
    dPP_Comb->Draw("same");
    
    mPP_Reco->SetMarkerStyle(27);
    mPP_Reco->SetMarkerColor(kBlack);
    mPP_Reco->Draw("same");
  
    TLegend *Spectra = myLegend(0.53,0.65,0.85,0.9);
    Spectra->AddEntry(uPP_Bayes,"Bayesian","pl");
    Spectra->AddEntry(dPP_Comb,"Measured","pl");
    Spectra->AddEntry(mPP_Reco,"MC Reco","pl");
    Spectra->SetTextSize(0.04);
    Spectra->Draw();
  
    cSpectra->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/Unfolded_spectra_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

  
    //rescale the histograms which were scaled: 
  

    for(int i = 0;i<=nbins_cent;i++){
    
      uPbPb_Bayes[i]->Scale(scaleFactor[i]);
      dPbPb_TrgComb[i]->Scale(scaleFactor[i]);
      mPbPb_Reco[i]->Scale(scaleFactor[i]);

    }

  }
 
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(doPbPbIterSys){
    //plot 1 - PbPb iteration systematics. 
    // this will be a 3 by 2 panel plot showing bayesian for PbPb. per centrality bin. 
    // divide different unfolding iterations with iteration 4 - the nominal one. 
    TCanvas *cIterSysPbPb = new TCanvas("cIterSysPbPb","PbPb Iteration systematics",1200,800);
    makeMultiPanelCanvasWithGap(cIterSysPbPb,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

    // line at 1
    TLine *linePbPb_iter = new TLine(30,1,299,1);
    linePbPb_iter->SetLineStyle(2);
    linePbPb_iter->SetLineWidth(2);
  
    TLegend *PbPb_itersys = myLegend(0.53,0.65,0.85,0.9);
  
    for(int i = 0;i<nbins_cent;i++){
      cIterSysPbPb->cd(nbins_cent-i);

      for(int j = 1;j<Iterations-3;j++){
	//uPbPb_BayesianIter[i][j]->Divide(uPbPb_BayesianIter[i][3]);
	uPbPb_BayesianIter[i][j]->SetMarkerStyle(33);
	uPbPb_BayesianIter[i][j]->SetMarkerColor(j+1);
	uPbPb_BayesianIter[i][j]->SetAxisRange(unfoldingCut,299,"X");
	uPbPb_BayesianIter[i][j]->SetAxisRange(0,2,"Y");

	if(j==0){
	  makeHistTitle(uPbPb_BayesianIter[i][j]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal(4 iterations))");
	  uPbPb_BayesianIter[i][j]->Draw();
	}
        uPbPb_BayesianIter[i][j]->Draw("same");
	
	if(i==0) PbPb_itersys->AddEntry(uPbPb_BayesianIter[i][j],Form("Iteration %d",j+1),"pl");

      }

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    
      linePbPb_iter->Draw();

    }
    PbPb_itersys->Draw();

    cIterSysPbPb->cd(1);
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
    drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

    cIterSysPbPb->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_%dGeVCut_ak%s%d%s_%d.pdf",unfoldingCut,algo,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPPIterSys){
    //plot 2 - pp iteration systematics 
    // this is just one panel plot showing bayesian for pp. 
    TCanvas *cIterSysPP = new TCanvas("cIterSysPP","PP Iteration systematics",600,400);

    TLegend *PP_itersys = myLegend(0.53,0.65,0.85,0.9);
    TLine *linePP_iter = new TLine(30,1,299,1);
    linePP_iter->SetLineStyle(2);
    linePP_iter->SetLineWidth(2);

    for(int i = 0;i<Iterations; i++){
      uPP_BayesianIter[i]->Divide(uPP_BayesianIter[3]);
      uPP_BayesianIter[i]->SetMarkerStyle(33);
      uPP_BayesianIter[i]->SetMarkerColor(i);
      uPP_BayesianIter[i]->SetAxisRange(30,299,"X");
      uPP_BayesianIter[i]->SetAxisRange(0,2,"Y");
      uPP_BayesianIter[i]->SetXTitle("Jet p_{T} (GeV/c)");
      uPP_BayesianIter[i]->SetYTitle("Ratio (unfolded/Nominal)");
      uPP_BayesianIter[i]->SetTitle(" ");
      uPP_BayesianIter[i]->GetXaxis()->CenterTitle();
      uPP_BayesianIter[i]->GetYaxis()->CenterTitle();

      if(i==0){
	//makeHistTitle(uPP_BayesianIter[i]," ","Jet p_{T} (GeV/c)","Ratio (unfolded/Nominal)");
	uPP_BayesianIter[i]->Draw();
      }
      uPP_BayesianIter[i]->Draw("same");
      PP_itersys->AddEntry(uPP_BayesianIter[i],Form("Iteration %d",i),"pl");
    }
    linePP_iter->Draw();
    PP_itersys->Draw();
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type,radius),0.2,0.23,20);
    drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

    cIterSysPP->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_%dGeVCut_ak%d%s_%d.pdf",unfoldingCut,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doRAA){
    //plot 3 - RAA 
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

      RAA_measured[i]->SetMarkerColor(kBlack);
      RAA_measured[i]->SetMarkerStyle(24);
      makeHistTitle(RAA_measured[i],"","Jet p_{T} (GeV/c)","R_{AA}");
      RAA_measured[i]->SetAxisRange(unfoldingCut+5,299,"X");
      RAA_measured[i]->SetAxisRange(0,2,"Y");
      RAA_measured[i]->Draw("E1");

      RAA_bayesian[i]->SetMarkerColor(kBlack);
      RAA_bayesian[i]->SetMarkerStyle(20);
      RAA_bayesian[i]->Draw("same E1");

      RAA_binbybin[i]->SetMarkerStyle(33);
      RAA_binbybin[i]->SetMarkerColor(kRed);
      RAA_binbybin[i]->Draw("same E1");

      lineRAA->Draw();
      lUnfoldingCut->Draw();
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.9,20);

    }
    
    tRAA->AddEntry(RAA_measured[0],"No Unfolding","pl");
    tRAA->AddEntry(RAA_bayesian[0],"Bayesian","pl");
    tRAA->AddEntry(RAA_binbybin[0],"BinbyBin","pl");
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

    cRAA->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_%dGeVCut_ak%s%d%s_%d.pdf",unfoldingCut,algo,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPbPbMCClosure){
    //plot 4 - PbPb MC closure test 
    // this will also be a 6 panel plot showing bayesian iteration 4 and binbybin divided by measured which is actually the MC
    TCanvas *cPbPbMCclosure = new TCanvas("cPbPbMCclosure","PbPb MC closure test",1200,800);
    makeMultiPanelCanvas(cPbPbMCclosure,3,2,0.0,0.0,0.2,0.15,0.07);

    TLine *linePbPb_mcclosure = new TLine(40,1,300,1);
    linePbPb_mcclosure->SetLineStyle(2);
    linePbPb_mcclosure->SetLineWidth(2);

    TH1F *hMCClosurePbPb_Meas[nbins_cent+1], *hMCClosurePbPb_Bayesian[nbins_cent+1], *hMCClosurePbPb_BinByBin[nbins_cent+1];
    for(int i = 0;i<nbins_cent;i++){

      cPbPbMCclosure->cd(nbins_cent-i);
      //hMCClosurePbPb_Meas[i] = (TH1F*)mPbPb_Reco[i]->Clone(Form("hMCClosurePbPb_Meas_cent%d",i));
      hMCClosurePbPb_Meas[i] = (TH1F*)mPbPb_mcclosure_data[i]->Rebin(nbins_pt,Form("hMCClosurePbPb_Meas_cent%d",i),boundaries_pt);
      hMCClosurePbPb_Bayesian[i] = (TH1F*)uPbPb_MC_Bayes[i]->Rebin(nbins_pt,Form("hMCClosurePbPb_Bayesian_cent%d",i),boundaries_pt);
      hMCClosurePbPb_BinByBin[i] = (TH1F*)uPbPb_MC_BinByBin[i]->Rebin(nbins_pt,Form("hMCClosurePbPb_BinByBin_cent%d",i),boundaries_pt);

      mPbPb_mcclosure_gen[i] = (TH1F*)mPbPb_mcclosure_gen[i]->Rebin(nbins_pt,Form("mc_pbpb_mcclosure_gen_cent%d",i),boundaries_pt);

      hMCClosurePbPb_Meas[i]->Divide(mPbPb_mcclosure_gen[i]);
      hMCClosurePbPb_Bayesian[i]->Divide(mPbPb_mcclosure_gen[i]);
      hMCClosurePbPb_BinByBin[i]->Divide(mPbPb_mcclosure_gen[i]);

      // hMCClosurePbPb_Meas[i]->Scale(2);
      // hMCClosurePbPb_Bayesian[i]->Scale(2);
      // hMCClosurePbPb_BinByBin[i]->Scale(2);

      // makeHistTitle(hMCClosurePbPb_Meas[i]," ","Jet p_{T} (GeV/c)","Reco/Truth");
      // hMCClosurePbPb_Meas[i]->SetMarkerStyle(24);
      // hMCClosurePbPb_Meas[i]->SetMarkerColor(kBlack);
      // hMCClosurePbPb_Meas[i]->SetAxisRange(30,500,"X");
      // hMCClosurePbPb_Meas[i]->SetAxisRange(0,2,"Y");
      //hMCClosurePbPb_Meas[i]->Draw();

      makeHistTitle(hMCClosurePbPb_Bayesian[i]," ","Jet p_{T} (GeV/c)","Reco/Truth");
      hMCClosurePbPb_Bayesian[i]->SetMarkerStyle(20);
      hMCClosurePbPb_Bayesian[i]->SetMarkerColor(kBlack);
      hMCClosurePbPb_Bayesian[i]->SetAxisRange(40,300,"X");
      hMCClosurePbPb_Bayesian[i]->SetAxisRange(0,2,"Y");
      hMCClosurePbPb_Bayesian[i]->Draw();

      hMCClosurePbPb_BinByBin[i]->SetMarkerColor(kRed);
      hMCClosurePbPb_BinByBin[i]->SetMarkerStyle(33);
      hMCClosurePbPb_BinByBin[i]->Draw("same");

      linePbPb_mcclosure->Draw();
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.6,0.33,16);

    }

    cPbPbMCclosure->cd(1);
    TLegend *pbpbmcclosure = myLegend(0.53,0.65,0.85,0.9);
    //pbpbmcclosure->AddEntry(hMCClosurePbPb_Meas[0],"PbPb no unfolding","pl");
    pbpbmcclosure->AddEntry(hMCClosurePbPb_Bayesian[0],"PbPb Bayesian 4 Iter","pl");
    pbpbmcclosure->AddEntry(hMCClosurePbPb_BinByBin[0],"PbPb BinbyBin","pl");
    pbpbmcclosure->Draw();
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("Data and gen, opp half statistics",0.3,0.13,16);
  
    cPbPbMCclosure->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_mc_closure_opphalf_%dGeVCut_ak%s%d%s_%d.pdf",unfoldingCut,algo,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPPMCClosure){
    //plot 5 - PP MC closure test
    TCanvas * cPPMCclosure = new TCanvas("cPPMCclosure","PP MC closure test",600,600);
    //TH1F *hMCClosurePP_Meas = (TH1F*)mPP_Reco->Clone("hMCClosurePP_Meas");
    mPP_mcclosure_gen = (TH1F*)mPP_mcclosure_gen->Rebin(nbins_pt,"mPP_mcclosure_gen",boundaries_pt);
    divideBinWidth(mPP_mcclosure_gen);
    TH1F *hMCClosurePP_Meas = (TH1F*)mPP_mcclosure_data->Clone("hMCClosurePP_Meas");
    TH1F *hMCClosurePP_Bayesian = (TH1F*)uPP_MC_Bayes->Clone("hMCClosurePP_Bayesian");
    TH1F *hMCClosurePP_BinbyBin = (TH1F*)uPP_MC_BinByBin->Clone("hMCClosurePP_BinbyBin");

    TLine *linePP_mcclosure = new TLine(40,1,300,1);
    linePP_mcclosure->SetLineStyle(2);
    linePP_mcclosure->SetLineWidth(2);

    hMCClosurePP_Bayesian->Divide(mPP_mcclosure_gen);
    hMCClosurePP_Meas->Divide(mPP_mcclosure_gen);
    hMCClosurePP_BinbyBin->Divide(mPP_mcclosure_gen);

    // hMCClosurePP_Meas->SetAxisRange(30,500,"X");
    // hMCClosurePP_Meas->SetAxisRange(0,2,"Y");

    // hMCClosurePP_Meas->SetMarkerStyle(24);
    // hMCClosurePP_Meas->SetMarkerColor(kBlack);

    // hMCClosurePP_Meas->SetTitle(" ");
    // hMCClosurePP_Meas->SetXTitle("Jet p_{T} (GeV/c)");
    // hMCClosurePP_Meas->SetYTitle("Reco/Truth");
  
    // hMCClosurePP_Meas->GetXaxis()->CenterTitle();
    // hMCClosurePP_Meas->GetYaxis()->CenterTitle();

    //hMCClosurePP_Meas->Draw();

    hMCClosurePP_Bayesian->SetTitle(" ");
    hMCClosurePP_Bayesian->SetXTitle("Jet p_{T} (GeV/c)");
    hMCClosurePP_Bayesian->SetYTitle("Reco/Truth");
  
    hMCClosurePP_Bayesian->GetXaxis()->CenterTitle();
    hMCClosurePP_Bayesian->GetYaxis()->CenterTitle();

    hMCClosurePP_Bayesian->SetAxisRange(unfoldingCut,299,"X");
    hMCClosurePP_Bayesian->SetAxisRange(0,2,"Y");

    hMCClosurePP_Bayesian->SetMarkerStyle(33);
    hMCClosurePP_Bayesian->SetMarkerColor(kRed);
    hMCClosurePP_Bayesian->Draw();

    hMCClosurePP_BinbyBin->SetMarkerStyle(29);
    hMCClosurePP_BinbyBin->SetMarkerColor(kBlue);
    hMCClosurePP_BinbyBin->Draw("same");

    TLegend *ppmcclosure = myLegend(0.53,0.65,0.85,0.9);
    //ppmcclosure->AddEntry(hMCClosurePP_Meas,"pp no unfolding","pl");
    ppmcclosure->AddEntry(hMCClosurePP_Bayesian,"pp Bayesian","pl");
    ppmcclosure->AddEntry(hMCClosurePP_BinbyBin,"pp BinbyBin","pl");
    ppmcclosure->Draw();

    linePP_mcclosure->Draw();

    putCMSPrel();
    drawText(Form("Anti-k_{T} %s Jets R=0.%d",jet_type,radius),0.2,0.23,20);
    drawText("|#eta|<2, |vz|<15",0.6,0.31,22);
  
    cPPMCclosure->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_mc_closure_test_%dGeVCut_ak%d%s_%d.pdf",unfoldingCut,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPbPbDatavsMC){
    //plot 6 - data vs MC for PbPb 
    // i want to make it like ratio plots for each centrality class. and plot the unfolded ratio as well. 
    // so i want a line at 1, and 2 overlayed plots showing measured/Gen ratio and unfo/gen ratio. 
    TCanvas *cPbPb_data_vs_mc = new TCanvas("cPbPb_data_vs_mc","PbPb data vs mc",1200,800);
    makeMultiPanelCanvasWithGap(cPbPb_data_vs_mc,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

    TH1F *PbPb_meas_vs_mc[nbins_cent+1];
    TH1F *PbPb_unfo_vs_mc[nbins_cent+1];

    TLine *linePbPb_data_vs_mc = new TLine(30,1,500,1);
    linePbPb_data_vs_mc->SetLineStyle(2);
    linePbPb_data_vs_mc->SetLineWidth(2);

    for(int i = 0;i<nbins_cent;i++){

      cPbPb_data_vs_mc->cd(nbins_cent-i);
      cPbPb_data_vs_mc->cd(nbins_cent-i)->SetLogy();

      PbPb_meas_vs_mc[i] = (TH1F*)dPbPb_TrgComb[i]->Clone("PbPb_meas_vs_mc");
      PbPb_unfo_vs_mc[i] = (TH1F*)uPbPb_Bayes[i]->Clone("PbPb_unfo_vs_mc");

      PbPb_meas_vs_mc[i]->Divide(mPbPb_Reco[i]);
      PbPb_unfo_vs_mc[i]->Divide(mPbPb_Reco[i]);

      makeHistTitle(PbPb_meas_vs_mc[i]," ","Jet p_{T} (GeV/c)", "#frac{data}{MC RECO}");
      PbPb_meas_vs_mc[i]->SetMarkerStyle(25);
      PbPb_meas_vs_mc[i]->SetMarkerColor(kBlack);
      PbPb_meas_vs_mc[i]->SetAxisRange(30,800,"X");
      PbPb_meas_vs_mc[i]->SetAxisRange(0,1e8,"Y");
      PbPb_meas_vs_mc[i]->Draw();

      PbPb_unfo_vs_mc[i]->SetMarkerStyle(27);
      PbPb_unfo_vs_mc[i]->SetMarkerColor(kBlue);
      PbPb_unfo_vs_mc[i]->Draw("same");

      linePbPb_data_vs_mc->Draw();
    
      drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

    }

    cPbPb_data_vs_mc->cd(1);
    TLegend *PbPb_datavsmc = myLegend(0.53,0.65,0.85,0.9);
    PbPb_datavsmc->AddEntry(PbPb_meas_vs_mc[0],"no unfolding","pl");
    PbPb_datavsmc->AddEntry(PbPb_unfo_vs_mc[0],"Bayesian 4 Iter","pl");
    PbPb_datavsmc->Draw();

    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.2,0.23,20);
    drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

    cPbPb_data_vs_mc->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_data_vs_mc_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPPDatavsMC){
    //plot7 - data vs MC for pp
    // made similar to above. 
    TCanvas *cPP_data_vs_mc = new TCanvas("cPP_data_vs_mc","PP data vs MC",600,600);

    TH1F *PP_meas_vs_mc = (TH1F*)dPP_Comb->Clone("PP_meas_vs_mc");
    TH1F *PP_unfo_vs_mc = (TH1F*)uPP_Bayes->Clone("PP_unfo_vs_mc");
  
    TLine *linePP_data_vs_mc = new TLine(30,1,500,1);
    linePP_data_vs_mc->SetLineStyle(2);
    linePP_data_vs_mc->SetLineWidth(2);

    PP_meas_vs_mc->Divide(mPP_Reco);
    PP_unfo_vs_mc->Divide(mPP_Reco);

    //makeHistTitle(PP_meas_vs_mc," ","Jet p_{T} (GeV/c)","#frac{data}{MC RECO}");

    PP_meas_vs_mc->SetTitle(" ");
    PP_meas_vs_mc->SetXTitle("Jet p_{T} (GeV/c)");
    PP_meas_vs_mc->SetYTitle("#frac{data}{MC RECO}");
    PP_meas_vs_mc->GetXaxis()->CenterTitle();
    PP_meas_vs_mc->GetYaxis()->CenterTitle();  

    PP_meas_vs_mc->SetMarkerStyle(25);
    PP_meas_vs_mc->SetMarkerColor(kBlack);
    PP_meas_vs_mc->Draw();
    PP_unfo_vs_mc->SetMarkerStyle(27);
    PP_unfo_vs_mc->SetMarkerColor(kBlue);
    PP_unfo_vs_mc->Draw("same");

    linePP_data_vs_mc->Draw();

    TLegend *PP_datavsmc = myLegend(0.53,0.65,0.85,0.9);
    PP_datavsmc->AddEntry(PP_meas_vs_mc,"no unfolding","pl");
    PP_datavsmc->AddEntry(PP_unfo_vs_mc,"Bayesian 4 Iter","pl");
    PP_datavsmc->SetTextSize(0.02);
    PP_datavsmc->Draw();

    putCMSPrel();
    drawText(Form("Anti-k_{T} %s Jets R=0.%d",jet_type,radius),0.2,0.23,20);
    drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

    cPP_data_vs_mc->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_data_vs_mc_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  
  if(doPbPbNormRes){
    //plot8 - normalized Response Matrix for PbPb. 
    TCanvas *cPbPb_NormResMat = new TCanvas("cPbPb_NormResMat","Normalized Response Matrix for PbPb",1200,800);
    makeMultiPanelCanvas(cPbPb_NormResMat,3,2,0.0,0.0,0.2,0.15,0.07);
  
    TLine *linePbPb_Res = new TLine(40,40,300,300);
    linePbPb_Res->SetLineStyle(2);
    linePbPb_Res->SetLineWidth(2);

    for(int i = 0;i<nbins_cent;i++){
      cPbPb_NormResMat->cd(nbins_cent-i);
      //cPbPb_NormResMat->cd(nbins_cent-i)->SetLogy();
      //cPbPb_NormResMat->cd(nbins_cent-i)->SetLogx();
      //cPbPb_NormResMat->cd(nbins_cent-i)->SetLogz();

      makeHistTitle(mPbPb_ResponseNorm[i]," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");

      mPbPb_ResponseNorm[i] = (TH2F*) mPbPb_ResponseNorm[i]->Rebin2D(10);

      mPbPb_ResponseNorm[i]->SetAxisRange(1e-10,1,"Z");
      mPbPb_ResponseNorm[i]->SetAxisRange(40,300,"X");
      mPbPb_ResponseNorm[i]->SetAxisRange(40,300,"Y");
   
      mPbPb_ResponseNorm[i]->Draw("textcolz");
    
      linePbPb_Res->Draw();
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.6,0.8,20);

    }
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.2,0.7,20);
    drawText("|#eta|<2, |vz|<15",0.65,0.31,16);
    //drawText("chMax/jtpt>0.05",0.65,0.4,16);
    //drawText("opp-side statistics",0.2,0.6,16);
    cPbPb_NormResMat->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_normalized_response_matrix_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPPNormRes){
    //plot9 - normalized response matrix for pp
    TCanvas *cPP_NormResMat = new TCanvas("cPP_NormResMat","Normalized Response Matrix fr PP",600,600);
    //cPP_NormResMat->SetLogy();
    //cPP_NormResMat->SetLogx();
    //cPP_NormResMat->SetLogz();

    TLine *linePP_Res = new TLine(40,40,300,300);
    linePP_Res->SetLineStyle(2);
    linePP_Res->SetLineWidth(2);

    //makeHistTitle(mPP_ResponseNorm," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");

    mPP_ResponseNorm = (TH2F*) mPP_ResponseNorm->Rebin2D(10);

    mPP_ResponseNorm->SetAxisRange(1e-10,1,"Z");
    mPP_ResponseNorm->SetAxisRange(40,300,"X");
    mPP_ResponseNorm->SetAxisRange(40,300,"Y");

    mPP_ResponseNorm->SetTitle(" ");
    mPP_ResponseNorm->SetXTitle("Gen p_{T} (GeV/c)");
    mPP_ResponseNorm->SetYTitle("Reco p_{T} (GeV/c)");
    mPP_ResponseNorm->GetXaxis()->CenterTitle();
    mPP_ResponseNorm->GetYaxis()->CenterTitle();

    mPP_ResponseNorm->Draw("textcolz");

    linePP_Res->Draw();

    putCMSPrel();
    drawText(Form("Anti-k_{T} %s Jets R=0.%d",jet_type,radius),0.2,0.23,20);
    drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
    cPP_NormResMat->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_normalized_response_matrix_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPbPbsigma){
    //plot10 - Cross section for PbPb - in proper units per centrality bin.  
    TCanvas *cPbPb_sigma = new TCanvas("cPbPb_sigma","PbPb inclusive jet invariant cross section",1200,800);
    makeMultiPanelCanvas(cPbPb_sigma,3,2,0.0,0.0,0.2,0.15,0.07); 

    for(int i = 0;i<nbins_cent;i++){

      cPbPb_sigma->cd(nbins_cent-i);
      cPbPb_sigma->cd(nbins_cent-i)->SetLogy();
      cPbPb_sigma->cd(nbins_cent-i)->SetLogx();

      makeHistTitle(dPbPb_TrgComb[i]," ","Jet p_{T} (GeV/c)","TAA cross section in nano barns");
      //dPbPb_TrgComb[i]->Scale(1e-6)
      dPbPb_TrgComb[i]->SetMarkerStyle(24);
      dPbPb_TrgComb[i]->SetMarkerColor(kBlack);
      dPbPb_TrgComb[i]->SetAxisRange(1e-6,100,"Y");
      //dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins_pt,Form("rebin_measured_cent%d",i),boundaries_pt);
      //divideBinWidth(dPbPb_TrgComb[i]);
      dPbPb_TrgComb[i]->SetAxisRange(50,500,"X");
      dPbPb_TrgComb[i]->Draw();
    
      uPbPb_Bayes[i]->SetMarkerStyle(20);
      uPbPb_Bayes[i]->SetMarkerColor(kBlack);
      //uPbPb_Bayes[i] = (TH1F*)uPbPb_Bayes[i]->Rebin(nbins_pt,Form("rebin_bayes_cent%d",i),boundaries_pt);
      //divideBinWidth(uPbPb_Bayes[i]);
      uPbPb_Bayes[i]->Draw("same");

      //uPbPb_BinByBin[i]->Scale(1e-6);
      uPbPb_BinByBin[i]->SetMarkerStyle(33);
      uPbPb_BinByBin[i]->SetMarkerColor(kRed);
      //uPbPb_BinByBin[i] = (TH1F*)uPbPb_BinByBin[i]->Rebin(nbins_pt,Form("rebin_binbybin_cent%d",i),boundaries_pt);
      //divideBinWidth(uPbPb_BinByBin[i]);
      uPbPb_BinByBin[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.8,20);
    }

    cPbPb_sigma->cd(1);
    TLegend *PbPb_sigma = myLegend(0.33,0.3,0.5,0.5);
    PbPb_sigma->AddEntry(dPbPb_TrgComb[0],"No Unfolding","pl");
    PbPb_sigma->AddEntry(uPbPb_Bayes[0],"Bayesian 4 Iter","pl");
    PbPb_sigma->AddEntry(uPbPb_BinByBin[0],"BinByBin","pl");
    PbPb_sigma->SetTextSize(0.04);
    PbPb_sigma->Draw();

    putCMSPrel();
    putPbPbLumi();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.3,0.2,16);
    //drawText("Spectra with NJet(p_{T}>50) < 7",0.55,0.55,16);

    cPbPb_sigma->cd(2);
    drawText("|#eta|<2, Jet ID cut",0.15,0.35,16);
    drawText("|vz|<15, pCES, HBHE",0.15,0.25,16);
    drawText("Pile up rejection cut",0.15,0.15,16);
    cPbPb_sigma->cd(3);
    //drawText("Marguerite file, ",0.1,0.2,16);

    cPbPb_sigma->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_invariant_cross_section_combined_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  if(doPPsigma){
    //plot 11 - cross section for PP 
    TCanvas *cPP_sigma = new TCanvas("cPP_sigma","PP inclusive jet invariant cross section",800,600);
    cPP_sigma->SetLogy();
    cPP_sigma->SetLogx();

    dPP_Comb->SetMarkerStyle(24);
    dPP_Comb->SetMarkerColor(kBlack);
    dPP_Comb->SetTitle(" ");
    dPP_Comb->SetXTitle("Jet p_{T} (GeV/c)");
    dPP_Comb->SetYTitle("#frac{d^2 #sigma}{d#eta dp_{T}} nano barns");
    dPP_Comb->SetAxisRange(40,300,"X");
    dPP_Comb->Draw();

    uPP_Bayes->SetMarkerColor(kBlack);
    uPP_Bayes->SetMarkerStyle(20);
    //uPP_Bayes->Scale(1./1e3);
    uPP_Bayes->Draw("same");

    uPP_BinByBin->SetMarkerStyle(33);
    uPP_BinByBin->SetMarkerColor(kRed);
    //uPP_BinByBin->Scale(1./1e3);
    uPP_BinByBin->Draw("same");

    TLegend *PP_sigma = myLegend(0.53,0.65,0.85,0.9);
    PP_sigma->AddEntry(dPP_Comb,"No Unfolding","pl");
    PP_sigma->AddEntry(uPP_Bayes,"Bayesian 4 Iter","pl");
    PP_sigma->AddEntry(uPP_BinByBin,"BinByBin","pl");
    PP_sigma->SetTextSize(0.04);
    PP_sigma->Draw();

    putCMSPrel();
    putPPLumi();
    drawText(Form("Anti-k_{T} %s Jets R=0.%d",jet_type,radius),0.15,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.15,0.33,16);
    cPP_sigma->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_invariant_cross_section_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  if(doGenSpectra){
    // plot 12 - Gen Spectra for PbPb and pp. 
  
    TCanvas *cGenSpectra = new TCanvas("cGenSpectra","Generator Level Spectra for PbPb and pp",1400,1000);
  
    cGenSpectra->Divide(1,1);
  
    cGenSpectra->cd(1);
    cGenSpectra->cd(1)->SetLogy();
    cGenSpectra->cd(1)->SetLogx();
    //cSpectra->Set
    //Double_t scaleFactor[nbins_cent+1] = {1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6};

    TLegend *PbPbGenSpectra = myLegend(0.4,0.65,0.85,0.9);

    for(int i = 0;i<=nbins_cent;i++){
    
      mPbPb_Gen[i]->Scale(1./scaleFactor[i]);
      mPbPb_Reco[i]->Scale(1./scaleFactor[i]); 

      mPbPb_Gen[i]->SetMarkerStyle(28);
      mPbPb_Gen[i]->SetMarkerColor(i+1);

      if(i==0){
	mPbPb_Gen[i]->SetAxisRange(1e-15,1e2,"Y");
	mPbPb_Gen[i]->SetAxisRange(15,600,"X");
	makeHistTitle(mPbPb_Gen[i],"","MC: Jet p_{T} (GeV/c)","arbitrary for now");
	mPbPb_Gen[i]->Draw();
      }else 
	mPbPb_Gen[i]->Draw("same");

      mPbPb_Reco[i]->SetMarkerStyle(29);
      mPbPb_Reco[i]->SetMarkerColor(i+1);
      mPbPb_Reco[i]->Draw("same");

    }

    for(int i = nbins_cent;i>=0;i--){
      if(i<nbins_cent)
	PbPbGenSpectra->AddEntry(mPbPb_Gen[i],Form("%2.0f - %2.0f cent * 10^{%d}",5*boundaries_cent[i],5*boundaries_cent[i+1],i),"pl");
      else 
	PbPbGenSpectra->AddEntry(mPbPb_Gen[i],Form("0-200 cent * 10^{%d}",i),"pl");
    }

    PbPbGenSpectra->SetTextSize(0.04);
    PbPbGenSpectra->Draw();
    putCMSSim();
    drawText("PbPb",0.2,0.8,20);
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.5,0.9,20);
    drawText("#frac{chMax}{jtpt}>0.01",0.15,0.2,20);
  
    // cGenSpectra->cd(2);						
    // cGenSpectra->cd(2)->SetLogy();
    // cGenSpectra->cd(2)->SetLogx();
    // makeHistTitle(mPP_Gen,"","MC: Jet p_{T} (GeV/c)","arbitrary for now");
    // mPP_Gen->SetMarkerStyle(28);
    // mPP_Gen->SetMarkerColor(kBlack);
    // mPP_Gen->SetAxisRange(1e-14,1,"Y");
    // mPP_Gen->SetAxisRange(15,600,"X");
    // mPP_Gen->Draw();
  
    // drawText("#frac{chMax}{jtpt}>0.01",0.15,0.2,20);
    // drawText(Form("Anti-k_{T} %s Jets R=0.%d",jet_type,radius),0.45,0.9,20);
  
    drawText("2.76 TeV, |#eta|<2, |vz|<15",0,0.9,20);
    //drawText("pp",0.3,0.8,20);

    // 
    // mPP_Reco->SetMarkerStyle(29);
    // mPP_Reco->SetMarkerColor(kBlack);
    // mPP_Reco->Draw("same");
  
    // TLegend *GenSpectra = myLegend(0.15,0.3,0.35,0.5);
    // GenSpectra->AddEntry(mPP_Gen,"MC GenJet - refpt","pl");
    // GenSpectra->AddEntry(mPP_Reco,"MC RecoJet - jtpt","pl");
    // GenSpectra->SetTextSize(0.04);
    // GenSpectra->Draw();
    // 

    cGenSpectra->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/MC_spectra_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
 
    //unscale the gen and reco histograms: 
    for(int i = 0;i<=nbins_cent;i++){ 
      mPbPb_Gen[i]->Scale(1./scaleFactor[i]);
      mPbPb_Reco[i]->Scale(1./scaleFactor[i]); 
    }
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPbPbTrgComb){
    // Plotting the trigObj spectra from individual triggers vs the combined triggers. 
    TCanvas *cDataMerge = new TCanvas("cDataMerge","",1000,800);
    makeMultiPanelCanvas(cDataMerge,3,2,0.0,0.0,0.2,0.15,0.07);

    for(int i = 0;i<nbins_cent;i++){

      cDataMerge->cd(nbins_cent-i);
      cDataMerge->cd(nbins_cent-i)->SetLogy();
      cDataMerge->cd(nbins_cent-i)->SetLogx();

      makeHistTitle(dPbPb_TrgComb[i]," ","Jet p_{T} (GeV/c)","TAA scaled #sigma in nano barns");
      //dPbPb_TrgComb[i]->Scale(1e-6)
      //dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins_pt,Form("rebin_pbpb_meas_cent%d",i),boundaries_pt);
      //divideBinWidth(dPbPb_TrgComb[i]);
      dPbPb_TrgComb[i]->SetMarkerStyle(25);
      dPbPb_TrgComb[i]->SetMarkerColor(kBlack);
      dPbPb_TrgComb[i]->SetAxisRange(1e-6,100,"Y");
      dPbPb_TrgComb[i]->SetAxisRange(50,500,"X");
      dPbPb_TrgComb[i]->Draw();

      //dPbPb_MinBias[i]->SetMarkerStyle(20);
      //dPbPb_MinBias[i]->SetMarkerColor(kBlack);
      //dPbPb_MinBias[i] = (TH1F*)dPbPb_MinBias[i]->Rebin(nbins_pt,Form("rebin_measured_80_cent%d",i),boundaries_pt);
      //divideBinWidth(dPbPb_MinBias[i]);
      //dPbPb_MinBias[i]->Draw("same");
      
      dPbPb_Trg80[i]->SetMarkerStyle(20);
      dPbPb_Trg80[i]->SetMarkerColor(kRed);
      //dPbPb_Trg80[i] = (TH1F*)dPbPb_Trg80[i]->Rebin(nbins_pt,Form("rebin_measured_80_cent%d",i),boundaries_pt);
      //divideBinWidth(dPbPb_Trg80[i]);
      dPbPb_Trg80[i]->Draw("same");

      dPbPb_Trg65[i]->SetMarkerStyle(20);
      dPbPb_Trg65[i]->SetMarkerColor(kGreen);
      //dPbPb_Trg65[i] = (TH1F*)dPbPb_Trg65[i]->Rebin(nbins_pt,Form("rebin_measured_65_cent%d",i),boundaries_pt);
      //divideBinWidth(dPbPb_Trg65[i]);
      dPbPb_Trg65[i]->Draw("same");

      dPbPb_Trg55[i]->SetMarkerStyle(20);
      dPbPb_Trg55[i]->SetMarkerColor(kBlue);
      //dPbPb_Trg55[i] = (TH1F*)dPbPb_Trg55[i]->Rebin(nbins_pt,Form("rebin_measured_55_cent%d",i),boundaries_pt);
      //divideBinWidth(dPbPb_Trg55[i]);
      dPbPb_Trg55[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.87,16);
    
    }

    cDataMerge->cd(1);
    TLegend *PbPb_dataMerge = myLegend(0.75,0.65,0.90,0.85);
    PbPb_dataMerge->AddEntry(dPbPb_TrgComb[0],"Combined","pl");
    //PbPb_dataMerge->AddEntry(dPbPb_MinBias[0],"JetMB","pl");
    PbPb_dataMerge->AddEntry(dPbPb_Trg80[0],"Jet80","pl");
    PbPb_dataMerge->AddEntry(dPbPb_Trg65[0],"Jet65","pl");
    PbPb_dataMerge->AddEntry(dPbPb_Trg55[0],"Jet55","pl");
    PbPb_dataMerge->SetTextSize(0.04);
    PbPb_dataMerge->Draw();

    putCMSPrel();
    putPbPbLumi();
    drawText("|#eta|<2, |vz|<15",0.25,0.21,16);
    drawText("pCES, HBHE",0.25,0.11,16);
    cDataMerge->cd(2);
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.1,0.95,16);
    //drawText("L1_SingleJet36 added to jet55 and 65 triggers",0.22,0.83,15);
    //drawText("cut5: #frac{neMax}{neMax + chMax + phMax}",0.22,0.73,15);
    //drawText("cut4: #frac{#sum #left(ch+nu+ph+mu+e#right)}{0.5*raw p_{T}}",0.22,0.83,15);
    //drawText("chMax/jtpt>0.02 && eMax/jtpt<0.6",0.15,0.11,16);
    cDataMerge->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_data_trigger_merging_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doPbPb12003TrgComb){

        // Plotting the trigObj spectra from individual triggers vs the combined triggers. 
    TCanvas *cData12003Merge = new TCanvas("cData12003Merge","",1000,800);
    makeMultiPanelCanvas(cData12003Merge,3,2,0.0,0.0,0.2,0.15,0.07);

    for(int i = 0;i<nbins_cent;i++){

      cData12003Merge->cd(nbins_cent-i);
      cData12003Merge->cd(nbins_cent-i)->SetLogy();
      cData12003Merge->cd(nbins_cent-i)->SetLogx();

      makeHistTitle(hPbPb12003Comb[i]," ","Jet p_{T} (GeV/c)","#frac{d^{2} N}{dp_{T}}");
      //hPbPb12003Comb[i]->Scale(1e-6)
      hPbPb12003Comb[i] = (TH1F*)hPbPb12003Comb[i]->Rebin(nbins_pt,Form("rebin_pbpb_meas_cent%d",i),boundaries_pt);
      divideBinWidth(hPbPb12003Comb[i]);
      hPbPb12003Comb[i]->SetMarkerStyle(25);
      hPbPb12003Comb[i]->SetMarkerColor(kBlack);
      hPbPb12003Comb[i]->SetAxisRange(1e-2,1e8,"Y");
      hPbPb12003Comb[i]->SetAxisRange(30,500,"X");
      hPbPb12003Comb[i]->Draw();

      hPbPb12003Jet80[i]->SetMarkerStyle(20);
      hPbPb12003Jet80[i]->SetMarkerColor(kRed);
      hPbPb12003Jet80[i] = (TH1F*)hPbPb12003Jet80[i]->Rebin(nbins_pt,Form("rebin_measured_80_cent%d",i),boundaries_pt);
      divideBinWidth(hPbPb12003Jet80[i]);
      hPbPb12003Jet80[i]->Draw("same");

      hPbPb12003Jet65[i]->SetMarkerStyle(20);
      hPbPb12003Jet65[i]->SetMarkerColor(kBlue);
      hPbPb12003Jet65[i] = (TH1F*)hPbPb12003Jet65[i]->Rebin(nbins_pt,Form("rebin_measured_65_cent%d",i),boundaries_pt);
      divideBinWidth(hPbPb12003Jet65[i]);
      hPbPb12003Jet65[i]->Draw("same");

      hPbPb12003Jet55[i]->SetMarkerStyle(20);
      hPbPb12003Jet55[i]->SetMarkerColor(kGreen);
      hPbPb12003Jet55[i] = (TH1F*)hPbPb12003Jet55[i]->Rebin(nbins_pt,Form("rebin_measured_55_cent%d",i),boundaries_pt);
      divideBinWidth(hPbPb12003Jet55[i]);
      hPbPb12003Jet55[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.7,20);
    
    }

    cData12003Merge->cd(1);
    TLegend *PbPb_data12003Merge = myLegend(0.53,0.75,0.8,0.9);
    PbPb_data12003Merge->AddEntry(hPbPb12003Comb[0],"Combined","pl");
    PbPb_data12003Merge->AddEntry(hPbPb12003Jet80[0],"Jet80","pl");
    PbPb_data12003Merge->AddEntry(hPbPb12003Jet65[0],"Jet65","pl");
    PbPb_data12003Merge->AddEntry(hPbPb12003Jet55[0],"Jet55","pl");
    PbPb_data12003Merge->SetTextSize(0.04);
    PbPb_data12003Merge->Draw();

    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.5,0.95,15);
    drawText("|#eta|<2, |vz|<15, cut5<0.975",0.25,0.61,16);
    drawText("pCES, HBHE",0.55,0.51,16);
    cData12003Merge->cd(2);
    drawText("L1_SingleJet36, 52 added to triggers",0.22,0.83,15);
    //drawText("cut4: #frac{#sum #left(ch+nu+ph+mu+e#right)}{0.5*raw p_{T}}",0.22,0.83,15);
    //drawText("cut1: #frac{chMax}{jet p_{T}}",0.22,0.73,15);
    drawText("cut5: #frac{neMax}{neMax + chMax + phMax}",0.22,0.73,15);
    drawText("12003 Method",0.22,0.63,15);
    cData12003Merge->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_data_trigger_merging_12003_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doPPTrgComb){
    // Plotting the trigger combination for the pp dataset 
    TCanvas *cPPMerge = new TCanvas("cPPMerge","",800,600);
    
    cPPMerge->SetLogy();
    cPPMerge->SetLogx();
    
    makeHistTitle(dPP_Comb," ","Jet p_{T} (GeV/c)","#sigma in nano barns");
    //dPP_Comb->Scale(1e-6)
    //dPP_Comb = (TH1F*)dPP_Comb->Rebin(nbins_pt,Form("rebin_pbpb_meas_cent%d",i),boundaries_pt);
    //divideBinWidth(dPP_Comb);
    dPP_Comb->SetMarkerStyle(25);
    dPP_Comb->SetMarkerColor(kBlack);
    dPP_Comb->SetAxisRange(1e-6,100,"Y");
    dPP_Comb->SetAxisRange(50,500,"X");
    dPP_Comb->Draw();

    dPP_1->SetMarkerStyle(20);
    dPP_1->SetMarkerColor(kRed);
    //dPP_1 = (TH1F*)dPP_1->Rebin(nbins_pt,Form("rebin_measured_80_cent%d",i),boundaries_pt);
    //divideBinWidth(dPP_1);
    dPP_1->Draw("same");

    dPP_2->SetMarkerStyle(20);
    dPP_2->SetMarkerColor(kBlue);
    //dPP_2 = (TH1F*)dPP_2->Rebin(nbins_pt,Form("rebin_measured_65_cent%d",i),boundaries_pt);
    //divideBinWidth(dPP_2);
    dPP_2->Draw("same");

    dPP_3->SetMarkerStyle(20);
    dPP_3->SetMarkerColor(kGreen);
    //dPP_3 = (TH1F*)dPP_3->Rebin(nbins_pt,Form("rebin_measured_55_cent%d",i),boundaries_pt);
    //divideBinWidth(dPP_3);
    dPP_3->Draw("same");

    TLegend *PP_dataMerge = myLegend(0.65,0.65,0.80,0.85);
    PP_dataMerge->AddEntry(dPP_Comb,"Combined","pl");
    PP_dataMerge->AddEntry(dPP_1,"Jet80","pl");
    PP_dataMerge->AddEntry(dPP_2,"Jet60","pl");
    PP_dataMerge->AddEntry(dPP_3,"Jet40","pl");
    PP_dataMerge->SetTextSize(0.04);
    PP_dataMerge->Draw();

    putCMSPrel();
    putPPLumi();
    drawText("|#eta|<2, |vz|<15",0.15,0.26,16);
    drawText("pCES, HBHE",0.15,0.21,16);
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.15,0.15,16);
    //drawText("L1_SingleJet36 added to jet55 and 65 triggers",0.22,0.83,15);
    //drawText("cut5: #frac{neMax}{neMax + chMax + phMax}",0.22,0.73,15);
    //drawText("cut4: #frac{#sum #left(ch+nu+ph+mu+e#right)}{0.5*raw p_{T}}",0.22,0.83,15);
    //drawText("chMax/jtpt>0.05",0.65,0.60,16);
    cPPMerge->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_data_trigger_merging_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");

  }
  
    
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doPPTrgContribution){
    // plot the pp trigger contribution 
    // get the data sepctra from the root files before the unfolding: 
    TLine *line = new TLine(40,1,140,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);

    TCanvas * cPPTrgCont = new TCanvas("cPPTrgCont","",800,600);

    hPP_80Cont->SetMarkerStyle(20);
    hPP_80Cont->SetMarkerColor(kBlack);
    hPP_80Cont->SetXTitle(Form("reco R=0.%d Jet p_{T} (GeV/c)",radius));
    hPP_80Cont->SetYTitle("Trigger Contribution");
    hPP_80Cont->SetTitle(" ");

    // hPP_80Cont = (TH1F*)hPP_80Cont->Rebin(nbins_pt, "hPP_80Cont", boundaries_pt);
    // divideBinWidth(hPP_80Cont);
    // hPP_60Cont = (TH1F*)hPP_60Cont->Rebin(nbins_pt, "hPP_60Cont", boundaries_pt);
    // divideBinWidth(hPP_60Cont);
    // hPP_40Cont = (TH1F*)hPP_40Cont->Rebin(nbins_pt, "hPP_40Cont", boundaries_pt);
    // divideBinWidth(hPP_40Cont);

    hPP_80Cont->SetAxisRange(40, 140, "X");
    hPP_80Cont->Draw();

    hPP_60Cont->SetMarkerStyle(20);
    hPP_60Cont->SetMarkerColor(kGreen);

    hPP_60Cont->Draw("same");

    hPP_40Cont->SetMarkerStyle(20);
    hPP_40Cont->SetMarkerColor(kBlue);

    hPP_40Cont->Draw("same");
    hPP_80Cont->Draw("same");

    line->Draw();

    TLegend *PP_Contribution = myLegend(0.65,0.25,0.80,0.45);
    PP_Contribution->AddEntry(hPP_80Cont,"Jet80","pl");
    PP_Contribution->AddEntry(hPP_60Cont,"Jet60","pl");
    PP_Contribution->AddEntry(hPP_40Cont,"Jet40","pl");
    PP_Contribution->SetTextSize(0.04);
    PP_Contribution->Draw();

    putCMSPrel();
    putPPLumi();
    drawText("|#eta|<2, |vz|<15",0.15,0.26,16);
    drawText("pCES, HBHE",0.15,0.21,16);
    
    cPPTrgCont->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_data_trigger_contribution_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");

  }

 
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  if(doRandomCone){

    // Plotting the average energy subtracted in the Vs cone 
    if(algo=="Vs"){


      TCanvas *cPbPb_MC_jtpu[nbins_eta];

      TF1 *fMC[no_radius][nbins_eta][nbins_cent];
      double c1FixMinX = 0.2;
      double c1xMax = 5;

      double fit_mean[no_radius][nbins_eta][nbins_cent];
      double fit_mean_err[no_radius][nbins_eta][nbins_cent];

      for(int j = 0;j<nbins_eta;j++){

	cPbPb_MC_jtpu[j] = new TCanvas(Form("cPbPb_MC_jtpu_%s",etaWidth[j]),Form("energy subtracted from Jets in the Vs algorithmin the range %s",etaWidth[j]),1400,1200);
	makeMultiPanelCanvas(cPbPb_MC_jtpu[j],3,2,0.0,0.0,0.2,0.15,0.07);
	TLegend *LPbPb_MC_jtpu[nbins_cent];

	for(int i = 0;i<nbins_cent;i++){

	  LPbPb_MC_jtpu[i] = myLegend(0.55,0.70,0.85,0.90);

	  cPbPb_MC_jtpu[j]->cd(nbins_cent-i)->SetLogy();
	  //cPbPb_MC_jtpu[j]->cd(i+1)->SetLogx();
	  //cPbPb_MC_jtpu[j]->cd(i+1)->SetGridy();
	  cPbPb_MC_jtpu[j]->cd(nbins_cent-i)->SetGridx();
	  cPbPb_MC_jtpu[j]->cd(nbins_cent-i);

	  for(int k = 1;k<5;k++){

	    hPbPb_MC_jtpu[k][j][i]->SetMarkerStyle(20+k);
	    hPbPb_MC_jtpu[k][j][i]->SetMarkerColor(k);
	    hPbPb_MC_jtpu[k][j][i]->SetAxisRange(5,250,"X");
	    hPbPb_MC_jtpu[k][j][i]->SetAxisRange(1,3.5e3,"Y");
	    makeHistTitle(hPbPb_MC_jtpu[k][j][i]," ","Subtracted p_{T} (GeV/c) from Jet","counts");
	 
	    fMC[k][j][i] = new TF1(Form("fMC_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),"[0]*exp(-0.5*((x-[1])/[2])^2)",5,250);
	    //fMC[k][j][i]->SetParameters(0.3,hPbPb_MC_jtpu[k][j][i]->GetMean(),10);
	    //fMC[k][j][i]->FixParameter(1,hPbPb_MC_jtpu[k][j][i]->GetMean());

	    hPbPb_MC_jtpu[k][j][i]->Fit(Form("fMC_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),"0");
	  
	    fit_mean[k][j][i] = fMC[k][j][i]->GetParameter(2);
	    fit_mean_err[k][j][i] = fMC[k][j][i]->GetParError(2);

	    // this fit values doesnt look reasonable??? 

	    if(k==1)
	      hPbPb_MC_jtpu[k][j][i]->Draw();
	    else 
	      hPbPb_MC_jtpu[k][j][i]->Draw("same");
	  
	    //LPbPb_MC_jtpu[i]->AddEntry(,"","");
	    LPbPb_MC_jtpu[i]->AddEntry(hPbPb_MC_jtpu[k][j][i],Form("R=0.%d #bar{p_{T}}:%5.2f #pm %5.2f",list_radius[k],hPbPb_MC_jtpu[k][j][i]->GetMean(),hPbPb_MC_jtpu[k][j][i]->GetMeanError()),"pl");
	    //LPbPb_MC_jtpu[i]->AddEntry("",Form("Fit mean:%5.2f #pm %5.2f",fit_mean[k][j][i],fit_mean_err[k][j][i]),"");

	  }//radius loop
	
	  LPbPb_MC_jtpu[i]->SetTextSize(0.03);
	  LPbPb_MC_jtpu[i]->Draw();
	  drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.87,0.95,20);

	}//centrality loop
	cPbPb_MC_jtpu[j]->cd(1);
	putCMSSim();
	cPbPb_MC_jtpu[j]->cd(nbins_cent);
	drawText(Form("Anti-k_{T} %s %s",algo,jet_type),0.5,0.55,16);
	drawText(Form("%s |vz|<15",etaWidth[j]),0.5,0.60,16);
      
	cPbPb_MC_jtpu[j]->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/bkg_energy_subtracted_%s_ak%s_%s_%d.pdf",etaWidth[j],algo,jet_type,date.GetDate()),"RECREATE");

      }//eta bins loop for the canvas

    
      // TCanvas *cMeanVsEta[4];
      // TH1F *hPbPb_MC_jtpu_mean_vs_eta[4];
      // hPbPb_MC_jtpu_mean_vs_eta[k] = new TH1F(Form("hPbPb_MC_jtpu_mean_vs_eta_R%d",list_radius[k]),"",12,-3,+3);

      // for(int k = 1;k<5;k++){
      
      //   cMeanVsEta[k-1] = new TCanvas(Form("cMeanVsEta_R%d",list_radius[k]),"",800,600);
      //   for(int n = 0;n<12;n++){
      // 	hPbPb_MC_jtpu_mean_vs_eta[k]->SetBin;
      //   }
  

    // plot from the Hydjet and MB data the jtpu variable as a function of centrality 
    
    }// random cone plotting if statement only for Vs so far

  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(doSupernovaData){
    // plotting the hiNpix vs NJets(pT>50 and |eta|<2) for each centrality bin. 
    TCanvas *cSupernova_data = new TCanvas("cSupernova_data","",1000,800);
    makeMultiPanelCanvas(cSupernova_data,3,2,0.0,0.0,0.2,0.15,0.07);
    //makeMultiPanelCanvasWithGap(cSupernova_data,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

    //TLine *lineSupernova_data = new TLine(0,38000,50,35500);
    TLine *lineSupernova_data = new TLine(0,38000,50,13000);
    lineSupernova_data->SetLineStyle(2);
    lineSupernova_data->SetLineWidth(2);

    TLine *lineSupernova_data_2 = new TLine(0,38000,50,35500);
    lineSupernova_data_2->SetLineStyle(2);
    lineSupernova_data_2->SetLineWidth(2);
    lineSupernova_data_2->SetLineColor(kRed);

    for(int i = 0;i<nbins_cent;i++){

      cSupernova_data->cd(nbins_cent-i);
    
      hpbpb_Npix_before_cut_data[i]->SetYTitle("No of Pixels");
      hpbpb_Npix_before_cut_data[i]->SetXTitle("No of Jets (p_{T}>50 GeV/c, |#eta|<2)");
      hpbpb_Npix_before_cut_data[i]->SetTitle(" ");
      hpbpb_Npix_before_cut_data[i]->Draw("col");

      lineSupernova_data->Draw();
      //lineSupernova_data_2->Draw();
      drawText("Data",0.6,0.2,16);
      drawText("|vz|<15, pCES, HBHE",0.45,0.75,16);
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.2,16);

    }
    cSupernova_data->cd(1);
    putCMSPrel();
    drawText("Full dataset Jet55 or Jet65 or Jet80",0.3,0.85,16);
    cSupernova_data->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_data_supernova_events_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doSupernovaMC){
    // plotting the hiNpix vs NJets(pT>50 and |eta|<2) for each centrality bin. for MC
    TCanvas *cSupernova_mc = new TCanvas("cSupernova_mc","",1000,800);
    makeMultiPanelCanvas(cSupernova_mc,3,2,0.0,0.0,0.2,0.15,0.07);

    //TLine *lineSupernova_mc = new TLine(0,38000,50,35500);
    TLine *lineSupernova_mc = new TLine(0,38000,50,13000);
    lineSupernova_mc->SetLineStyle(2);
    lineSupernova_mc->SetLineWidth(2);

    TLine *lineSupernova_mc_2 = new TLine(0,38000,50,35500);
    lineSupernova_mc_2->SetLineStyle(2);
    lineSupernova_mc_2->SetLineWidth(2);
    lineSupernova_mc_2->SetLineColor(kRed);

    for(int i = 0;i<nbins_cent;i++){

      cSupernova_mc->cd(nbins_cent-i);
    
      hpbpb_Npix_before_cut_mc[i]->SetYTitle("No of Pixels");
      hpbpb_Npix_before_cut_mc[i]->SetTitle(" ");
      hpbpb_Npix_before_cut_mc[i]->SetXTitle("No of Jets (p_{T}>50 GeV/c, |#eta|<2)");
      hpbpb_Npix_before_cut_mc[i]->Draw("col");

      lineSupernova_mc->Draw();
      //lineSupernova_mc_2->Draw();
      drawText("MC",0.6,0.2,16);
      drawText("|vz|<15, pCES",0.45,0.75,16);
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.2,16);
    }
    cSupernova_mc->cd(1);
    putCMSSim();
    drawText("Pythia + HYDJET",0.3,0.85,16);
    cSupernova_mc->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_mc_supernova_events_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  
  }
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doSupernovaDataMCUltraCentral){

    TCanvas *cSupernovaDataMC = new TCanvas("cSupernovaDataMC","",1000,800);
    makeMultiPanelCanvas(cSupernovaDataMC,2,1,0.0,0.0,0.2,0.15,0.07);
    
    //TLine *lineSupernova_mc = new TLine(0,38000,50,35500);
    TLine *lineSupernova_mc = new TLine(0,38000,50,13000);
    lineSupernova_mc->SetLineStyle(2);
    lineSupernova_mc->SetLineWidth(2);
    
    cSupernovaDataMC->cd(1);
    
    hpbpb_Npix_before_cut_data[7]->SetYTitle("No of Pixels");
    hpbpb_Npix_before_cut_data[7]->SetTitle(" ");
    hpbpb_Npix_before_cut_data[7]->SetXTitle("No of Jets (p_{T}>50 GeV/c, |#eta|<2)");
    hpbpb_Npix_before_cut_data[7]->Draw("col");
    
    lineSupernova_mc->Draw();
    drawText("Data",0.6,0.2,16);
    //drawText("|vz|<15, pCES",0.55,0.75,16);
    putCMSPrel();
    
    cSupernovaDataMC->cd(2);
    
    hpbpb_Npix_before_cut_mc[7]->SetYTitle("No of Pixels");
    hpbpb_Npix_before_cut_mc[7]->SetTitle(" ");
    hpbpb_Npix_before_cut_mc[7]->SetXTitle("No of Jets (p_{T}>50 GeV/c, |#eta|<2)");
    hpbpb_Npix_before_cut_mc[7]->Draw("col");
    
    lineSupernova_mc->Draw();
    drawText("Pythia + HYDJET",0.3,0.2,16);
    drawText(" 0 - 0.5% central events",0.3,0.85,16);
    drawText("|vz|<15, pCES",0.55,0.65,16);

    cSupernovaDataMC->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_DataMC_supernova_ultraCentral_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");    
    
  }
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doSupernovaContributionData){
    // make the plot which shows the contribution of the cut. 
    TCanvas *cpT_njet_gl7[nbins_cent+1];

    for(int i = 0;i<=nbins_cent;i++){
      cpT_njet_gl7[i] = new TCanvas(Form("cpT_njet_gl7_cent%d",i),"",1000,1000);
    
      makeMultiPanelCanvas(cpT_njet_gl7[i],2,2,0.0,0.0,0.2,0.15,0.07);    
      cpT_njet_gl7[i]->cd(1);
      cpT_njet_gl7[i]->cd(1)->SetLogy();
      cpT_njet_gl7[i]->cd(1)->SetGridx();

      makeHistTitle(dPbPb_Trg55[i]," ","Jet p_{T} (GeV/c)","counts/dp_{T}");
      dPbPb_Trg55[i] = (TH1F*)dPbPb_Trg55[i]->Rebin(nbins_pt,Form("rebin_trgobj55_cent%d",i),boundaries_pt);
      divideBinWidth(dPbPb_Trg55[i]);
      hpbpb_TrgObj55_nJet_l7[i] = (TH1F*)hpbpb_TrgObj55_nJet_l7[i]->Rebin(nbins_pt,Form("rebin_trgobj55_l7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObj55_nJet_l7[i]);
      hpbpb_TrgObj55_nJet_g7[i] = (TH1F*)hpbpb_TrgObj55_nJet_g7[i]->Rebin(nbins_pt,Form("rebin_trgobj55_g7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObj55_nJet_g7[i]);

      dPbPb_Trg55[i]->SetMarkerColor(kBlack);
      dPbPb_Trg55[i]->SetMarkerStyle(25);
      dPbPb_Trg55[i]->SetAxisRange(30,400,"X");
      dPbPb_Trg55[i]->SetAxisRange(1e-2,1e8,"Y");
      dPbPb_Trg55[i]->Draw();

      hpbpb_TrgObj55_nJet_l7[i]->SetMarkerColor(kBlack);
      hpbpb_TrgObj55_nJet_l7[i]->SetMarkerStyle(20);
      hpbpb_TrgObj55_nJet_l7[i]->Draw("same");
    
      hpbpb_TrgObj55_nJet_g7[i]->SetMarkerColor(kRed);
      hpbpb_TrgObj55_nJet_g7[i]->SetMarkerStyle(20);
      hpbpb_TrgObj55_nJet_g7[i]->Draw("same");
    
      drawText("Jet55, 55<TrgObjpT<65",0.5,0.65,16);

      TLegend *PbPb_lg7 = myLegend(0.5,0.75,0.75,0.9);
      PbPb_lg7->AddEntry(dPbPb_Trg55[i],"Total","pl");
      PbPb_lg7->AddEntry(hpbpb_TrgObj55_nJet_l7[i],"NJet < 7 (Jet pT>50)","pl");
      PbPb_lg7->AddEntry(hpbpb_TrgObj55_nJet_g7[i],"NJet > 7 (Jet pT>50)","pl");
      PbPb_lg7->SetTextSize(0.04);
      PbPb_lg7->Draw();

      putCMSPrel();
      if(i!=nbins_cent)drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);
      else drawText("0-100 %",0.75,0.9,20);
      drawText("|#eta|<2, chMax/jtpt>0.01",0.5,0.55,16);
      drawText("|vz|<15, pCES, HBHE",0.5,0.45,16);
      drawText("hiNpix_1 > 38000 - 500*NJet",0.5,0.35,16);

      cpT_njet_gl7[i]->cd(2);
      cpT_njet_gl7[i]->cd(2)->SetLogy();
      cpT_njet_gl7[i]->cd(2)->SetGridx();

      makeHistTitle(dPbPb_Trg65[i]," ","Jet p_{T} (GeV/c)","counts/dp_{T}");
      dPbPb_Trg65[i] = (TH1F*)dPbPb_Trg65[i]->Rebin(nbins_pt,Form("rebin_trgobj65_cent%d",i),boundaries_pt);
      divideBinWidth(dPbPb_Trg65[i]);
      hpbpb_TrgObj65_nJet_l7[i] = (TH1F*)hpbpb_TrgObj65_nJet_l7[i]->Rebin(nbins_pt,Form("rebin_trgobj65_l7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObj65_nJet_l7[i]);
      hpbpb_TrgObj65_nJet_g7[i] = (TH1F*)hpbpb_TrgObj65_nJet_g7[i]->Rebin(nbins_pt,Form("rebin_trgobj65_g7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObj65_nJet_g7[i]);

      dPbPb_Trg65[i]->SetMarkerColor(kBlack);
      dPbPb_Trg65[i]->SetMarkerStyle(25);
      dPbPb_Trg65[i]->SetAxisRange(30,400,"X");
      dPbPb_Trg65[i]->SetAxisRange(1e-2,1e8,"Y");
      dPbPb_Trg65[i]->Draw();

      hpbpb_TrgObj65_nJet_l7[i]->SetMarkerColor(kBlack);
      hpbpb_TrgObj65_nJet_l7[i]->SetMarkerStyle(20);
      hpbpb_TrgObj65_nJet_l7[i]->Draw("same");
    
      hpbpb_TrgObj65_nJet_g7[i]->SetMarkerColor(kRed);
      hpbpb_TrgObj65_nJet_g7[i]->SetMarkerStyle(20);
      hpbpb_TrgObj65_nJet_g7[i]->Draw("same");

      drawText("Jet65, 65<TrgObjpT<80",0.5,0.65,16);
      drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type,radius),0.5,0.55,16);

      cpT_njet_gl7[i]->cd(3);
      cpT_njet_gl7[i]->cd(3)->SetLogy();
      cpT_njet_gl7[i]->cd(3)->SetGridx();

      makeHistTitle(dPbPb_Trg80[i]," ","Jet p_{T} (GeV/c)","counts/dp_{T}");
      dPbPb_Trg80[i] = (TH1F*)dPbPb_Trg80[i]->Rebin(nbins_pt,Form("rebin_trgobj80_cent%d",i),boundaries_pt);
      divideBinWidth(dPbPb_Trg80[i]);
      hpbpb_TrgObj80_nJet_l7[i] = (TH1F*)hpbpb_TrgObj80_nJet_l7[i]->Rebin(nbins_pt,Form("rebin_trgobj80_l7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObj80_nJet_l7[i]);
      hpbpb_TrgObj80_nJet_g7[i] = (TH1F*)hpbpb_TrgObj80_nJet_g7[i]->Rebin(nbins_pt,Form("rebin_trgobj80_g7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObj80_nJet_g7[i]);

      dPbPb_Trg80[i]->SetMarkerColor(kBlack);
      dPbPb_Trg80[i]->SetMarkerStyle(25);
      dPbPb_Trg80[i]->SetAxisRange(30,400,"X");
      dPbPb_Trg80[i]->SetAxisRange(1e-2,1e8,"Y");
      dPbPb_Trg80[i]->Draw();

      hpbpb_TrgObj80_nJet_l7[i]->SetMarkerColor(kBlack);
      hpbpb_TrgObj80_nJet_l7[i]->SetMarkerStyle(20);
      hpbpb_TrgObj80_nJet_l7[i]->Draw("same");
    
      hpbpb_TrgObj80_nJet_g7[i]->SetMarkerColor(kRed);
      hpbpb_TrgObj80_nJet_g7[i]->SetMarkerStyle(20);
      hpbpb_TrgObj80_nJet_g7[i]->Draw("same");

      drawText("Jet80, 80<TrgObjpT",0.5,0.65,16);

      cpT_njet_gl7[i]->cd(4);
      cpT_njet_gl7[i]->cd(4)->SetLogy();
      cpT_njet_gl7[i]->cd(4)->SetGridx();

      makeHistTitle(dPbPb_TrgComb[i]," ","Jet p_{T} (GeV/c)","counts/dp_{T}");
      dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins_pt,Form("rebin_trgobjComb_cent%d",i),boundaries_pt);
      divideBinWidth(dPbPb_TrgComb[i]);
      hpbpb_TrgObjComb_nJet_l7[i] = (TH1F*)hpbpb_TrgObjComb_nJet_l7[i]->Rebin(nbins_pt,Form("rebin_trgobjComb_l7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObjComb_nJet_l7[i]);
      hpbpb_TrgObjComb_nJet_g7[i] = (TH1F*)hpbpb_TrgObjComb_nJet_g7[i]->Rebin(nbins_pt,Form("rebin_trgobjComb_g7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_TrgObjComb_nJet_g7[i]);

      dPbPb_TrgComb[i]->SetMarkerColor(kBlack);
      dPbPb_TrgComb[i]->SetMarkerStyle(25);
      dPbPb_TrgComb[i]->SetAxisRange(30,400,"X");
      dPbPb_TrgComb[i]->SetAxisRange(1e-2,1e8,"Y");
      dPbPb_TrgComb[i]->Draw();

      hpbpb_TrgObjComb_nJet_l7[i]->SetMarkerColor(kBlack);
      hpbpb_TrgObjComb_nJet_l7[i]->SetMarkerStyle(20);
      hpbpb_TrgObjComb_nJet_l7[i]->Draw("same");
    
      hpbpb_TrgObjComb_nJet_g7[i]->SetMarkerColor(kRed);
      hpbpb_TrgObjComb_nJet_g7[i]->SetMarkerStyle(20);
      hpbpb_TrgObjComb_nJet_g7[i]->Draw("same");

      drawText("Jet Trigger Combined",0.5,0.65,16);

      cpT_njet_gl7[i]->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_data_spectra_nJetcut_ak%s%d%s_cent%d_%d.pdf",algo,radius,jet_type,i,date.GetDate()),"RECREATE");

    }//cent loop
  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doSupernovaContributionMC){
    // make the same plot as above, showing the contributions from the NJet>7 vs NJet<7 for MC.
  
    TCanvas *cMC_njet_lg7 = new TCanvas("nMC_njet_lg7","",1000,800);
    makeMultiPanelCanvas(cMC_njet_lg7,3,2,0.0,0.0,0.2,0.15,0.07);

    for(int i = 0;i<nbins_cent;i++){

      cMC_njet_lg7->cd(nbins_cent-i);
      cMC_njet_lg7->cd(nbins_cent-i)->SetLogy();
    
      mPbPb_Reco[i] = (TH1F*)mPbPb_Reco[i]->Rebin(nbins_pt,Form("rebin_mc_jtpt_cent%d",i),boundaries_pt);
      divideBinWidth(mPbPb_Reco[i]);
      mPbPb_Reco[i]->Scale(4);
      hpbpb_pt_Njet_l7[i] = (TH1F*)hpbpb_pt_Njet_l7[i]->Rebin(nbins_pt,Form("rebin_mc_jtpt_l7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_pt_Njet_l7[i]);
      hpbpb_pt_Njet_g7[i] = (TH1F*)hpbpb_pt_Njet_g7[i]->Rebin(nbins_pt,Form("rebin_mc_jtpt_g7_cent%d",i),boundaries_pt);
      divideBinWidth(hpbpb_pt_Njet_g7[i]);
    
      makeHistTitle(mPbPb_Reco[i]," ","Jet RECO p_{T} (GeV/c)","d#sigma/dp_{T}");
      mPbPb_Reco[i]->SetMarkerColor(kBlack);
      mPbPb_Reco[i]->SetMarkerStyle(25);
      mPbPb_Reco[i]->SetAxisRange(30,500,"X");
      mPbPb_Reco[i]->SetAxisRange(1e-12,1e-2,"Y");
      mPbPb_Reco[i]->Draw();

      hpbpb_pt_Njet_l7[i]->SetMarkerColor(kBlack);
      hpbpb_pt_Njet_l7[i]->SetMarkerStyle(20);
      hpbpb_pt_Njet_l7[i]->Draw("same");

      hpbpb_pt_Njet_g7[i]->SetMarkerColor(kRed);
      hpbpb_pt_Njet_g7[i]->SetMarkerStyle(20);
      hpbpb_pt_Njet_g7[i]->Draw("same");
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.85,20);

    }

    cMC_njet_lg7->cd(1);
    
    putCMSSim();
    TLegend *PbPb_mc_lg7 = myLegend(0.5,0.75,0.75,0.9);
    PbPb_mc_lg7->AddEntry(mPbPb_Reco[0],"Total","pl");
    PbPb_mc_lg7->AddEntry(hpbpb_pt_Njet_l7[0],"NJet < 7 (Jet pT>50)","pl");
    PbPb_mc_lg7->AddEntry(hpbpb_pt_Njet_g7[0],"NJet > 7 (Jet pT>50)","pl");
    PbPb_mc_lg7->SetTextSize(0.04);
    PbPb_mc_lg7->Draw();

    drawText("|#eta|<2, chMax/jtpt>0.01",0.5,0.55,16);
    drawText("|vz|<15, pCES",0.5,0.45,16);
    //drawText("hiNpix_1 > 38000 - 500*NJet",0.5,0.35,16);
    
    cMC_njet_lg7->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_mc_spectra_nJetcut_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  
  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doCentDatavsMC){

    TH1F *hDataCent = (TH1F*)fDatain->Get(Form("hpbpb_cent_R%d",radius));
    TH1F *hMCCent = (TH1F*)fMCin->Get(Form("hpbpb_cent_R%d",radius));

    TCanvas *cDataMCCent = new TCanvas("cDataMCCent","",800,600);
    makeMultiPanelCanvas(cDataMCCent,1,2,0.0,0.0,0.2,0.15,0.07);
	
    cDataMCCent->cd(1)->SetLogy();
    hDataCent->Scale(1./hDataCent->Integral());
    hMCCent->Scale(1./hMCCent->Integral());

    makeHistTitle(hDataCent," ","hiBin","Event Fraction");
    hDataCent->SetMarkerStyle(24);
    hDataCent->SetMarkerColor(kBlack);
    hDataCent->Draw();
    
    hMCCent->SetLineColor(kRed);
    hMCCent->SetLineStyle(1);
    hMCCent->Draw("same L");

    TLegend *pbpbCent = myLegend(0.77,0.65,0.9,0.9);
    pbpbCent->AddEntry(hDataCent,"Data","p");
    pbpbCent->AddEntry(hMCCent,"MC","l");
    pbpbCent->SetTextSize(0.04);
    pbpbCent->Draw();
    putCMSPrel();    
    //drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    //drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);
    
    cDataMCCent->cd(2);
    cDataMCCent->cd(2)->SetLogy();
    TH1F* hCentRatio = (TH1F*)hDataCent->Clone("hCentRatio");
    hCentRatio->Divide(hMCCent);

    TF1 *fCentralityWeight = new TF1("fCentralityWeight","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,200);

    hCentRatio->SetMarkerStyle(20);
    hCentRatio->SetMarkerColor(kBlack);
    makeHistTitle(hCentRatio," ","hiBin","Data/MC");
    hCentRatio->Fit("fCentralityWeight","LL","",0,200);
    hCentRatio->Draw();
    
    cDataMCCent->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_datavsmc_centrality_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    // TFile fout(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_DataMC_cent_ratio_%d.root",date.GetDate()),"RECREATE");
    // fout.cd();
    // hCentRatio->Write();
    // fout.Close();

  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doVertexDatavsMC){

    TH1F *hDataVx = (TH1F*)fDatain->Get(Form("hpbpb_vx_R%d",radius));
    TH1F *hMCVx = (TH1F*)fMCin->Get(Form("hpbpb_vx_R%d",radius));  

    TH1F *hDataVy = (TH1F*)fDatain->Get(Form("hpbpb_vy_R%d",radius));
    TH1F *hMCVy = (TH1F*)fMCin->Get(Form("hpbpb_vy_R%d",radius));

    TH1F *hDataVz = (TH1F*)fDatain->Get(Form("hpbpb_vz_R%d",radius));
    TH1F *hMCVz = (TH1F*)fMCin->Get(Form("hpbpb_vz_R%d",radius));

    TCanvas *cDataMCVertex = new TCanvas("cDataMCVertex","",600,600);
    makeMultiPanelCanvas(cDataMCVertex,1,2,0.0,0.0,0.2,0.15,0.07);

    hDataVx->Scale(1./hDataVx->Integral());
    hDataVy->Scale(1./hDataVy->Integral());
    hDataVz->Scale(1./hDataVz->Integral());

    hMCVx->Scale(1./hMCVx->Integral());
    hMCVy->Scale(1./hMCVy->Integral());
    hMCVz->Scale(1./hMCVz->Integral());
    /*
    cDataMCVertex->cd(1);
    cDataMCVertex->cd(1)->SetLogy();
    makeHistTitle(hDataVx," ","Vx","Event Fraction");
    hDataVx->SetMarkerStyle(24);
    hDataVx->SetMarkerColor(kBlack);
    hDataVx->Draw();

    hMCVx->SetLineColor(kRed);
    hMCVx->SetLineStyle(1);
    hMCVx->Draw("same L");
    
    TLegend *pbpbVx = myLegend(0.77,0.65,0.9,0.9);
    pbpbVx->AddEntry(hDataVx,"Data","p");
    pbpbVx->AddEntry(hMCVx,"MC","l");
    pbpbVx->SetTextSize(0.04);
    pbpbVx->Draw();
    putCMSPrel();    
    //drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    //drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);

    cDataMCVertex->cd(4);
    TH1F *hVxRatio = (TH1F*)hDataVx->Clone("hVxRatio");
    hVxRatio->Divide(hMCVx);
    hVxRatio->SetMarkerColor(kBlack);
    hVxRatio->SetMarkerStyle(20);
    hVxRatio->SetAxisRange(0,2,"X");
    makeHistTitle(hVxRatio," ","Vx","Data/MC");
    hVxRatio->Draw();

    cDataMCVertex->cd(2);
    cDataMCVertex->cd(2)->SetLogy();
    makeHistTitle(hDataVy," ","Vy","Event Fraction");
    hDataVy->SetMarkerStyle(24);
    hDataVy->SetMarkerColor(kBlack);
    hDataVy->Draw();

    hMCVy->SetLineColor(kRed);
    hMCVy->SetLineStyle(1);
    hMCVy->Draw("same L");

    cDataMCVertex->cd(5);
    TH1F *hVyRatio = (TH1F*)hDataVy->Clone("hVyRatio");
    hVyRatio->Divide(hMCVy);
    hVyRatio->SetMarkerColor(kBlack);
    hVyRatio->SetMarkerStyle(20);
    hVyRatio->SetAxisRange(0,2,"X");
    makeHistTitle(hVyRatio," ","Vy","Data/MC");
    hVyRatio->Draw();
    */
    cDataMCVertex->cd(1);
    cDataMCVertex->cd(1)->SetLogy();
    makeHistTitle(hDataVz," ","Vz","Event Fraction");
    hDataVz->SetMarkerStyle(24);
    hDataVz->SetMarkerColor(kBlack);
    hDataVz->Draw();

    hMCVz->SetLineColor(kRed);
    hMCVz->SetLineStyle(1);
    hMCVz->Draw("same L");

    cDataMCVertex->cd(2);
    TH1F *hVzRatio = (TH1F*)hDataVz->Clone("hVzRatio");
    hVzRatio->Divide(hMCVz);

    TF1 *fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-15,15);
    hVzRatio->Fit("fVz","","",-15,15);
    
    hVzRatio->SetMarkerColor(kBlack);
    hVzRatio->SetMarkerStyle(20);
    //hVzRatio->SetAxisRange(0,2,"X");
    makeHistTitle(hVzRatio," ","Vz","Data/MC");
    hVzRatio->Draw();
    
    cDataMCVertex->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_datavsmc_vertex_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
 
  
  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if(doPFElectronCheck){

    TFile * fhistin = TFile::Open("/Users/keraghav/WORK/RAA/Output/RAA_JetID_Data_mc_ElecCutRejectionStudy_with2DplotsSumMaxChNePhCut_Pu3PF_20141216.root");
    TFile * fMCsubidin = TFile::Open("/Users/keraghav/WORK/RAA/Output/RAA_JetID_mc_subid0_ElecCutRejectionStudy_with2DplotsSumMaxChNePhCut_Pu3PF_20141217.root");
    cout<<"checkpoint 1"<<endl;
    // [3] is 0: Data, 1: MC no subid cut, 2: MC with subid==0 cut; [9] is eMax/jtpt<0.1 to 0.9  
    TH1F *hElecFracRatio[3][9][nbins_cent+1];
    //TH1F *hData_nocut[3][nbins_cent+1];

    for(int i = 0;i<nbins_cent+1;i++){
      for(int q = 0;q<9;q++){
	hElecFracRatio[0][q][i] = (TH1F*)fhistin->Get(Form("hData_Ratio_emaxJtPt_0p%d_cent%d",q+1,i));
	hElecFracRatio[1][q][i] = (TH1F*)fhistin->Get(Form("hMC_Ratio_emaxJtPt_0p%d_cent%d",q+1,i));
	hElecFracRatio[2][q][i] = (TH1F*)fMCsubidin->Get(Form("hMC_Ratio_emaxJtPt_0p%d_cent%d",q+1,i));

	hElecFracRatio[0][q][i] = (TH1F*)hElecFracRatio[0][q][i]->Rebin(nbins_pt,Form("hData_Ratio_emaxJtPt_0p%d_cent%d",q+1,i),boundaries_pt);
	hElecFracRatio[1][q][i] = (TH1F*)hElecFracRatio[1][q][i]->Rebin(nbins_pt,Form("hMC_Ratio_emaxJtPt_0p%d_cent%d",q+1,i),boundaries_pt);
	hElecFracRatio[2][q][i] = (TH1F*)hElecFracRatio[2][q][i]->Rebin(nbins_pt,Form("hMC_subid_Ratio_emaxJtPt_0p%d_cent%d",q+1,i),boundaries_pt);

	divideBinWidth(hElecFracRatio[0][q][i]);
	divideBinWidth(hElecFracRatio[1][q][i]);
	divideBinWidth(hElecFracRatio[2][q][i]);

      }

    }

    cout<<"hi check"<<endl;

    //make the plots: 

    TCanvas *cData = new TCanvas("cData","",800,600);
    makeMultiPanelCanvas(cData,2,2,0.0,0.0,0.2,0.15,0.07);

    cData->cd(1);
    for(int q = 0;q<9;q++){
      makeHistTitle(hElecFracRatio[0][q][nbins_cent],"","Jet p_{T} (GeV/c)","With Cut / without cut");
      hElecFracRatio[0][q][nbins_cent]->SetMarkerStyle(29);
      hElecFracRatio[0][q][nbins_cent]->SetMarkerColor(q+1);
      hElecFracRatio[0][q][nbins_cent]->SetAxisRange(0,550,"X");
      hElecFracRatio[0][q][nbins_cent]->SetAxisRange(0.4,1.1,"Y");
      if(q==1)hElecFracRatio[0][q][nbins_cent]->Draw();
      else hElecFracRatio[0][q][nbins_cent]->Draw("same");
    }
    putCMSPrel();
    drawText("Data",0.3,0.2,20);

    cData->cd(2);
    for(int q = 0;q<9;q++){
      makeHistTitle(hElecFracRatio[1][q][nbins_cent],"","Jet p_{T} (GeV/c)","With Cut / without cut");
      hElecFracRatio[1][q][nbins_cent]->SetMarkerStyle(29);
      hElecFracRatio[1][q][nbins_cent]->SetMarkerColor(q+1);
      hElecFracRatio[1][q][nbins_cent]->SetAxisRange(0,550,"X");
      hElecFracRatio[1][q][nbins_cent]->SetAxisRange(0.4,1.1,"Y");
      if(q==1)hElecFracRatio[1][q][nbins_cent]->Draw();
      else hElecFracRatio[1][q][nbins_cent]->Draw("same");
    }
    //purCMSPrel();
    drawText("MC - no SubID cut",0.3,0.2,20);
    
    cData->cd(3);
    for(int q = 0;q<9;q++){
      makeHistTitle(hElecFracRatio[2][q][nbins_cent],"","Jet p_{T} (GeV/c)","With Cut / without cut");
      hElecFracRatio[2][q][nbins_cent]->SetMarkerStyle(29);
      hElecFracRatio[2][q][nbins_cent]->SetMarkerColor(q+1);
      hElecFracRatio[2][q][nbins_cent]->SetAxisRange(0,550,"X");
      hElecFracRatio[2][q][nbins_cent]->SetAxisRange(0.8,1.1,"Y");
      if(q==1)hElecFracRatio[2][q][nbins_cent]->Draw();
      else hElecFracRatio[2][q][nbins_cent]->Draw("same");
    }
    //purCMSPrel();
    drawText("MC - SubID==0 cut",0.3,0.2,20);

    cData->cd(4);
    TLegend *Lelec = myLegend(0.15,0.1,0.45,0.9);
    TH1F *htest = new TH1F("htest","",10,0,10);
    //Lelec->AddEntry(htest,"eMax/jtpt < ","");
    for(int q = 8;q>=0;q--) Lelec->AddEntry(hElecFracRatio[0][q][nbins_cent],Form("eMax/jtpt < 0.%d",q+1),"pl");
    Lelec->SetTextSize(0.04);
    Lelec->Draw();

    drawText("0 - 100 %",0.75,0.9,16);
    drawText(Form("ak%s%d%s Jets",algo,radius,jet_type),0.5,0.8,20);
    drawText("|#eta|<2, |vz|<15",0.5,0.7,16);
    drawText("pCES, HBHE",0.5,0.6,16);
    drawText("Npix Diag cut",0.5,0.5,16);
    
    cData->SaveAs(Form("/Users/keraghav/WORK/RAA/Plots/PbPb_eMaxJtPt_ratio_centFull_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if(doJetVariablesCheck){

    cout<<endl<<"Plotting JetVarialbles Check"<<endl;

    //get the jetID trees from the MC and data files 
    // TTree *jetID_data = (TTree*)fDatain->Get("jets_ID");
    // TTree *jetID_mc = (TTree*)fMCin->Get("jets_ID");

    TFile * fData = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src//Output/PbPb_marguerite_chMaxjtpt_norawptcut_electronFailureNtuple_file_spectra_histograms_akPuPF_20150205.root");
    TFile * fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_mc_chMaxjtpt_norawptcut_spectra_akPuPF_20150210.root");
    
    // declare the histograms
    TH1F *hchMax[2][nbins_cent+1], *hchSum[2][nbins_cent+1], *hphMax[2][nbins_cent+1], *hphSum[2][nbins_cent+1], *hneMax[2][nbins_cent+1], *hneSum[2][nbins_cent+1], *hmuMax[2][nbins_cent+1], *hmuSum[2][nbins_cent+1], *heMax[2][nbins_cent+1], *heSum[2][nbins_cent+1], *hjtpu[2][nbins_cent+1]; 
    
    for(int i = 0;i<=nbins_cent;i++){
      cout<<"cent = "<<i<<endl;
      hchMax[0][i] = (TH1F*)fData->Get(Form("hpbpb_chMax_R3_n20_eta_p20_cent%d",i));
      hchSum[0][i] = (TH1F*)fData->Get(Form("hpbpb_chSum_R3_n20_eta_p20_cent%d",i));
      hphMax[0][i] = (TH1F*)fData->Get(Form("hpbpb_phMax_R3_n20_eta_p20_cent%d",i));
      hphSum[0][i] = (TH1F*)fData->Get(Form("hpbpb_phSum_R3_n20_eta_p20_cent%d",i));
      hneMax[0][i] = (TH1F*)fData->Get(Form("hpbpb_neMax_R3_n20_eta_p20_cent%d",i));
      hneSum[0][i] = (TH1F*)fData->Get(Form("hpbpb_neSum_R3_n20_eta_p20_cent%d",i));
      hmuMax[0][i] = (TH1F*)fData->Get(Form("hpbpb_muMax_R3_n20_eta_p20_cent%d",i));
      hmuSum[0][i] = (TH1F*)fData->Get(Form("hpbpb_muSum_R3_n20_eta_p20_cent%d",i));
      heMax[0][i]  = (TH1F*)fData->Get(Form("hpbpb_eMax_R3_n20_eta_p20_cent%d",i));
      heSum[0][i]  = (TH1F*)fData->Get(Form("hpbpb_eSum_R3_n20_eta_p20_cent%d",i));
      // hjtpu[0][i] = new TH1F(Form("hpbpb_jtpu_R3_n20_eta_p20_cent%d",i));
      
      hchMax[1][i] = (TH1F*)fMC->Get(Form("hpbpb_chMax_R3_n20_eta_p20_cent%d",i));
      hchSum[1][i] = (TH1F*)fMC->Get(Form("hpbpb_chSum_R3_n20_eta_p20_cent%d",i));
      hphMax[1][i] = (TH1F*)fMC->Get(Form("hpbpb_phMax_R3_n20_eta_p20_cent%d",i));
      hphSum[1][i] = (TH1F*)fMC->Get(Form("hpbpb_phSum_R3_n20_eta_p20_cent%d",i));
      hneMax[1][i] = (TH1F*)fMC->Get(Form("hpbpb_neMax_R3_n20_eta_p20_cent%d",i));
      hneSum[1][i] = (TH1F*)fMC->Get(Form("hpbpb_neSum_R3_n20_eta_p20_cent%d",i));
      hmuMax[1][i] = (TH1F*)fMC->Get(Form("hpbpb_muMax_R3_n20_eta_p20_cent%d",i));
      hmuSum[1][i] = (TH1F*)fMC->Get(Form("hpbpb_muSum_R3_n20_eta_p20_cent%d",i));
      heMax[1][i]  = (TH1F*)fMC->Get(Form("hpbpb_eMax_R3_n20_eta_p20_cent%d",i));
      heSum[1][i]  = (TH1F*)fMC->Get(Form("hpbpb_eSum_R3_n20_eta_p20_cent%d",i));
      //hjtpu[1][i] = new TH1F(Form("hpbpb_jtpu_R3_n20_eta_p20_cent%d",i));

      hchMax[0][i] = (TH1F*)hchMax[0][i]->Rebin(nbins_pt,Form("hpbpb_chMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hphMax[0][i] = (TH1F*)hphMax[0][i]->Rebin(nbins_pt,Form("hpbpb_phMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hneMax[0][i] = (TH1F*)hneMax[0][i]->Rebin(nbins_pt,Form("hpbpb_neMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hmuMax[0][i] = (TH1F*)hmuMax[0][i]->Rebin(nbins_pt,Form("hpbpb_muMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      heMax[0][i] = (TH1F*)heMax[0][i]->Rebin(nbins_pt,Form("hpbpb_eMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hchSum[0][i] = (TH1F*)hchSum[0][i]->Rebin(nbins_pt,Form("hpbpb_chSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hphSum[0][i] = (TH1F*)hphSum[0][i]->Rebin(nbins_pt,Form("hpbpb_phSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hneSum[0][i] = (TH1F*)hneSum[0][i]->Rebin(nbins_pt,Form("hpbpb_neSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hmuSum[0][i] = (TH1F*)hmuSum[0][i]->Rebin(nbins_pt,Form("hpbpb_muSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      heSum[0][i] = (TH1F*)heSum[0][i]->Rebin(nbins_pt,Form("hpbpb_eSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);

      hchMax[1][i] = (TH1F*)hchMax[1][i]->Rebin(nbins_pt,Form("hpbpb_chMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hphMax[1][i] = (TH1F*)hphMax[1][i]->Rebin(nbins_pt,Form("hpbpb_phMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hneMax[1][i] = (TH1F*)hneMax[1][i]->Rebin(nbins_pt,Form("hpbpb_neMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hmuMax[1][i] = (TH1F*)hmuMax[1][i]->Rebin(nbins_pt,Form("hpbpb_muMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      heMax[1][i] = (TH1F*)heMax[1][i]->Rebin(nbins_pt,Form("hpbpb_eMax_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hchSum[1][i] = (TH1F*)hchSum[1][i]->Rebin(nbins_pt,Form("hpbpb_chSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hphSum[1][i] = (TH1F*)hphSum[1][i]->Rebin(nbins_pt,Form("hpbpb_phSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hneSum[1][i] = (TH1F*)hneSum[1][i]->Rebin(nbins_pt,Form("hpbpb_neSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      hmuSum[1][i] = (TH1F*)hmuSum[1][i]->Rebin(nbins_pt,Form("hpbpb_muSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);
      heSum[1][i] = (TH1F*)heSum[1][i]->Rebin(nbins_pt,Form("hpbpb_eSum_R3_n20_eta_p20_cent%d",i),boundaries_pt);

      divideBinWidth(hchMax[0][i]);
      divideBinWidth(hphMax[0][i]);
      divideBinWidth(hneMax[0][i]);
      divideBinWidth(hmuMax[0][i]);
      divideBinWidth(heMax[0][i]);
      divideBinWidth(hchSum[0][i]);
      divideBinWidth(hphSum[0][i]);
      divideBinWidth(hneSum[0][i]);
      divideBinWidth(hmuSum[0][i]);
      divideBinWidth(heSum[0][i]);

      divideBinWidth(hchMax[1][i]);
      divideBinWidth(hphMax[1][i]);
      divideBinWidth(hneMax[1][i]);
      divideBinWidth(hmuMax[1][i]);
      divideBinWidth(heMax[1][i]);
      divideBinWidth(hchSum[1][i]);
      divideBinWidth(hphSum[1][i]);
      divideBinWidth(hneSum[1][i]);
      divideBinWidth(hmuSum[1][i]);
      divideBinWidth(heSum[1][i]);
            
      /*
      //draw the histograms from the trees: 
      jetID_data->Draw(Form("chMax>>hchMax_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_data->Draw(Form("chSum>>hchSum_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("chMax>>hchMax_mc_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("chSum>>hchSum_mc_cent%d",i),Form("cent==%d",i),"goff");   
      cout<<"got charged histograms for data and MC"<<endl;
      jetID_data->Draw(Form("phMax>>hphMax_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_data->Draw(Form("phSum>>hphSum_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("phMax>>hphMax_mc_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("phSum>>hphSum_mc_cent%d",i),Form("cent==%d",i),"goff");   
      cout<<"got photon histograms for data and MC"<<endl;
      jetID_data->Draw(Form("neMax>>hneMax_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_data->Draw(Form("neSum>>hneSum_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("neMax>>hneMax_mc_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("neSum>>hneSum_mc_cent%d",i),Form("cent==%d",i),"goff");       
      cout<<"got neutral histograms for data and MC"<<endl;
      jetID_data->Draw(Form("muMax>>hmuMax_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_data->Draw(Form("muSum>>hmuSum_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("muMax>>hmuMax_mc_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("muSum>>hmuSum_mc_cent%d",i),Form("cent==%d",i),"goff"); 
      cout<<"got muon histograms for data and MC"<<endl;
      jetID_data->Draw(Form("eMax>>heMax_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_data->Draw(Form("eSum>>heSum_data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("eMax>>heMax_mc_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("eSum>>heSum_mc_cent%d",i),Form("cent==%d",i),"goff"); 
      cout<<"got electron histograms for data and MC"<<endl;
      //jetID_data->Draw(Form("jtpu>>hjtpu_data_cent%d",i),Form("cent==%d",i),"goff");
      //jetID_mc->Draw(Form("jtpu>>hjtpu_mc_cent%d",i),Form("cent==%d",i),"goff"); 
      //cout<<"got charged histograms for data and MC"<<endl;
      */
      
      hchMax[0][i]->Scale(1./hchMax[0][i]->GetEntries());   
      hchSum[0][i]->Scale(1./hchSum[0][i]->GetEntries());
      //hchSum[1][i]->Scale(1./hchSum[1][i]->GetEntries());
      //hchMax[1][i]->Scale(1./hchMax[1][i]->GetEntries());
      hphSum[0][i]->Scale(1./hphSum[0][i]->GetEntries());
      hphMax[0][i]->Scale(1./hphMax[0][i]->GetEntries());
      //hphSum[1][i]->Scale(1./hphSum[1][i]->GetEntries());
      //hphMax[1][i]->Scale(1./hphMax[1][i]->GetEntries());
      hneSum[0][i]->Scale(1./hneSum[0][i]->GetEntries());
      hneMax[0][i]->Scale(1./hneMax[0][i]->GetEntries());
      //hneSum[1][i]->Scale(1./hneSum[1][i]->GetEntries());
      //hneMax[1][i]->Scale(1./hneMax[1][i]->GetEntries());
      hmuSum[0][i]->Scale(1./hmuSum[0][i]->GetEntries());
      hmuMax[0][i]->Scale(1./hmuMax[0][i]->GetEntries());
      //hmuSum[1][i]->Scale(1./hmuSum[1][i]->GetEntries());
      //hmuMax[1][i]->Scale(1./hmuMax[1][i]->GetEntries());
      heSum[0][i]->Scale(1./heSum[0][i]->GetEntries());
      heMax[0][i]->Scale(1./heMax[0][i]->GetEntries());
      //heSum[1][i]->Scale(1./heSum[1][i]->GetEntries());
      //heMax[1][i]->Scale(1./heMax[1][i]->GetEntries());
      
      //hjtpu[0][i]->Scale(1./hjtpu[0][i]->GetEntries());
      //hjtpu[1][i]->Scale(1./hjtpu[1][i]->GetEntries());

    }


    //make the plots:
    // before doing it in differnet centralities lets do it for the whole range. so two 3x3 plot (with gap) with data and mc in each panel where panel 1 is chMax, 2 is phMax, 3 is neMax, 4 is muMax and 5 is eMax and 6 is where we draw the legend and the useful information about the cuts etc... the second plot will do the same thing for Sum variables. 

    TCanvas *cMaxVariables = new TCanvas("cMaxVariables","",1000,800);
    makeMultiPanelCanvasWithGap(cMaxVariables,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

    cMaxVariables->cd(1);
    cMaxVariables->cd(1)->SetLogy();
    makeHistTitle(hchMax[0][nbins_cent],"","chMax p_{T} (GeV/c)","Event Fraction");
    hchMax[0][nbins_cent]->SetMarkerStyle(24);
    hchMax[0][nbins_cent]->SetMarkerColor(kBlack);
    hchMax[0][nbins_cent]->SetAxisRange(0,500,"X");
    hchMax[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hchMax[0][nbins_cent]->Draw();
   
    hchMax[1][nbins_cent]->SetLineStyle(1);
    hchMax[1][nbins_cent]->SetLineColor(kRed);
    hchMax[1][nbins_cent]->Draw("same L");
    putCMSPrel();
    
    cMaxVariables->cd(2);
    cMaxVariables->cd(2)->SetLogy();
    makeHistTitle(hphMax[0][nbins_cent],"","phMax p_{T} (GeV/c)","Event Fraction");
    hphMax[0][nbins_cent]->SetMarkerStyle(24);
    hphMax[0][nbins_cent]->SetMarkerColor(kBlack);
    hphMax[0][nbins_cent]->SetAxisRange(0,500,"X");
    hphMax[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hphMax[0][nbins_cent]->Draw();
    
    hphMax[1][nbins_cent]->SetLineStyle(1);
    hphMax[1][nbins_cent]->SetLineColor(kRed);
    hphMax[1][nbins_cent]->Draw("same L");

    cMaxVariables->cd(3);
    cMaxVariables->cd(3)->SetLogy();
    makeHistTitle(hneMax[0][nbins_cent],"","neMax p_{T} (GeV/c)","Event Fraction");
    hneMax[0][nbins_cent]->SetMarkerStyle(24);
    hneMax[0][nbins_cent]->SetMarkerColor(kBlack);
    hneMax[0][nbins_cent]->SetAxisRange(0,500,"X");
    hneMax[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hneMax[0][nbins_cent]->Draw();
    
    hneMax[1][nbins_cent]->SetLineStyle(1);
    hneMax[1][nbins_cent]->SetLineColor(kRed);
    hneMax[1][nbins_cent]->Draw("same L");

    cMaxVariables->cd(4);
    cMaxVariables->cd(4)->SetLogy();
    makeHistTitle(hmuMax[0][nbins_cent],"","muMax p_{T} (GeV/c)","Event Fraction");
    hmuMax[0][nbins_cent]->SetMarkerStyle(24);
    hmuMax[0][nbins_cent]->SetMarkerColor(kBlack);
    hmuMax[0][nbins_cent]->SetAxisRange(0,500,"X");
    hmuMax[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hmuMax[0][nbins_cent]->Draw();
    
    hmuMax[1][nbins_cent]->SetLineStyle(33);
    hmuMax[1][nbins_cent]->SetLineColor(kRed);
    hmuMax[1][nbins_cent]->Draw("same L");

    cMaxVariables->cd(5);
    cMaxVariables->cd(5)->SetLogy();
    makeHistTitle(heMax[0][nbins_cent],"","eMax p_{T} (GeV/c)","Event Fraction");
    heMax[0][nbins_cent]->SetMarkerStyle(24);
    heMax[0][nbins_cent]->SetMarkerColor(kBlack);
    heMax[0][nbins_cent]->SetAxisRange(0,500,"X");
    heMax[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    heMax[0][nbins_cent]->Draw();
    
    heMax[1][nbins_cent]->SetLineStyle(33);
    heMax[1][nbins_cent]->SetLineColor(kRed);
    heMax[1][nbins_cent]->Draw("same L");

    cMaxVariables->cd(6);
    drawText("0-100 %",0.55,0.9,18);
    TLegend *pbpbMax = myLegend(0.5,0.65,0.9,0.9);
    pbpbMax->AddEntry(hchMax[0][nbins_cent],"Data","p");
    pbpbMax->AddEntry(hchMax[1][nbins_cent],"MC","l");
    pbpbMax->SetTextSize(0.06);
    pbpbMax->Draw();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);

    cMaxVariables->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_Max_centFull_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
      
    //

    TCanvas *cSumVariables = new TCanvas("cSumVariables","",1000,800);
    makeMultiPanelCanvasWithGap(cSumVariables,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

    cSumVariables->cd(1);
    cSumVariables->cd(1)->SetLogy();
    makeHistTitle(hchSum[0][nbins_cent],"","chSum p_{T} (GeV/c)","Event Fraction");
    hchSum[0][nbins_cent]->SetMarkerStyle(24);
    hchSum[0][nbins_cent]->SetMarkerColor(kBlack);
    hchSum[0][nbins_cent]->SetAxisRange(0,500,"X");
    hchSum[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hchSum[0][nbins_cent]->Draw();
    
    hchSum[1][nbins_cent]->SetLineStyle(33);
    hchSum[1][nbins_cent]->SetLineColor(kRed);
    hchSum[1][nbins_cent]->Draw("same L");
    putCMSPrel();    

    cSumVariables->cd(2);
    cSumVariables->cd(2)->SetLogy();
    makeHistTitle(hphSum[0][nbins_cent],"","phSum p_{T} (GeV/c)","Event Fraction");
    hphSum[0][nbins_cent]->SetMarkerStyle(24);
    hphSum[0][nbins_cent]->SetMarkerColor(kBlack);
    hphSum[0][nbins_cent]->SetAxisRange(0,500,"X");
    hphSum[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hphSum[0][nbins_cent]->Draw();
    
    hphSum[1][nbins_cent]->SetLineStyle(33);
    hphSum[1][nbins_cent]->SetLineColor(kRed);
    hphSum[1][nbins_cent]->Draw("same L");

    cSumVariables->cd(3);
    cSumVariables->cd(3)->SetLogy();
    makeHistTitle(hneSum[0][nbins_cent],"","neSum p_{T} (GeV/c)","Event Fraction");
    hneSum[0][nbins_cent]->SetMarkerStyle(24);
    hneSum[0][nbins_cent]->SetMarkerColor(kBlack);
    hneSum[0][nbins_cent]->SetAxisRange(0,500,"X");
    hneSum[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hneSum[0][nbins_cent]->Draw();
    
    hneSum[1][nbins_cent]->SetLineStyle(33);
    hneSum[1][nbins_cent]->SetLineColor(kRed);
    hneSum[1][nbins_cent]->Draw("same L");

    cSumVariables->cd(4);
    cSumVariables->cd(4)->SetLogy();
    makeHistTitle(hmuSum[0][nbins_cent],"","muSum p_{T} (GeV/c)","Event Fraction");
    hmuSum[0][nbins_cent]->SetMarkerStyle(24);
    hmuSum[0][nbins_cent]->SetMarkerColor(kBlack);
    hmuSum[0][nbins_cent]->SetAxisRange(0,500,"X");
    hmuSum[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    hmuSum[0][nbins_cent]->Draw();
    
    hmuSum[1][nbins_cent]->SetLineStyle(33);
    hmuSum[1][nbins_cent]->SetLineColor(kRed);
    hmuSum[1][nbins_cent]->Draw("same L");

    cSumVariables->cd(5);
    cSumVariables->cd(5)->SetLogy();
    makeHistTitle(heSum[0][nbins_cent],"","eSum p_{T} (GeV/c)","Event Fraction");
    heSum[0][nbins_cent]->SetMarkerStyle(24);
    heSum[0][nbins_cent]->SetMarkerColor(kBlack);
    heSum[0][nbins_cent]->SetAxisRange(0,500,"X");
    heSum[0][nbins_cent]->SetAxisRange(1e-11,1,"Y");
    heSum[0][nbins_cent]->Draw();
    
    heSum[1][nbins_cent]->SetLineStyle(33);
    heSum[1][nbins_cent]->SetLineColor(kRed);
    heSum[1][nbins_cent]->Draw("same L");

    cSumVariables->cd(6);
    drawText("0-100 %",0.55,0.9,18);
    TLegend *pbpbSum = myLegend(0.5,0.65,0.9,0.9);
    pbpbSum->AddEntry(hchSum[0][nbins_cent],"Data","p");
    pbpbSum->AddEntry(hchSum[1][nbins_cent],"MC","l");
    pbpbSum->SetTextSize(0.06);
    pbpbSum->Draw();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);

    cSumVariables->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_Sum_centFull_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    //
    //
    //

#if 0
    TCanvas *cchMax = new TCanvas("cchMax","",1000,800);
    // make plot for chMax
    makeMultiPanelCanvas(cchMax,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){
      cchMax->cd(nbins_cent-i);
      cchMax->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hchMax[0][i]," ","chMax p_{T} (GeV/c)","Event Fraction");
      hchMax[0][i]->SetMarkerStyle(24);
      hchMax[0][i]->SetMarkerColor(kBlack);
      hchMax[0][i]->SetAxisRange(0,500,"X");
      hchMax[0][i]->Draw();

      hchMax[1][i]->SetLineStyle(1);
      hchMax[1][i]->SetLineColor(kRed);
      hchMax[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cchMax->cd(1);
    TLegend *pbpbchMax = myLegend(0.77,0.65,0.9,0.9);
    pbpbchMax->AddEntry(hchMax[0][0],"Data","p");
    pbpbchMax->AddEntry(hchMax[1][0],"MC","l");
    pbpbchMax->SetTextSize(0.04);
    pbpbchMax->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cchMax->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_chMax_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    
    // make plot for chSum
    TCanvas *cchSum = new TCanvas("cchSum","",1000,800);
    makeMultiPanelCanvas(cchSum,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cchSum->cd(nbins_cent-i);
      cchSum->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hchSum[0][i]," ","chSum p_{T} (GeV/c)","Event Fraction");
      hchSum[0][i]->SetMarkerStyle(24);
      hchSum[0][i]->SetMarkerColor(kBlack);
      hchSum[0][i]->SetAxisRange(0,500,"X");
      hchSum[0][i]->Draw();

      hchSum[1][i]->SetLineStyle(1);
      hchSum[1][i]->SetLineColor(kRed);
      hchSum[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cchSum->cd(1);
    TLegend *pbpbchSum = myLegend(0.77,0.65,0.9,0.9);
    pbpbchSum->AddEntry(hchSum[0][0],"Data","p");
    pbpbchSum->AddEntry(hchSum[1][0],"MC","l");
    pbpbchSum->SetTextSize(0.04);
    pbpbchSum->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cchSum->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_chSum_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");


        
    // make plot for phMax
    TCanvas *cphMax = new TCanvas("cphMax","",1000,800);
    makeMultiPanelCanvas(cphMax,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cphMax->cd(nbins_cent-i);
      cphMax->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hphMax[0][i]," ","phMax p_{T} (GeV/c)","Event Fraction");
      hphMax[0][i]->SetMarkerStyle(24);
      hphMax[0][i]->SetMarkerColor(kBlack);
      hphMax[0][i]->SetAxisRange(0,500,"X");
      hphMax[0][i]->Draw();

      hphMax[1][i]->SetLineStyle(1);
      hphMax[1][i]->SetLineColor(kRed);
      hphMax[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cphMax->cd(1);
    TLegend *pbpbphMax = myLegend(0.77,0.65,0.9,0.9);
    pbpbphMax->AddEntry(hphMax[0][0],"Data","p");
    pbpbphMax->AddEntry(hphMax[1][0],"MC","l");
    pbpbphMax->SetTextSize(0.04);
    pbpbphMax->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cphMax->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_phMax_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    // make plot for phSum
    TCanvas *cphSum = new TCanvas("cphSum","",1000,800);
    makeMultiPanelCanvas(cphSum,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cphSum->cd(nbins_cent-i);
      cphSum->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hphSum[0][i]," ","phSum p_{T} (GeV/c)","Event Fraction");
      hphSum[0][i]->SetMarkerStyle(24);
      hphSum[0][i]->SetMarkerColor(kBlack);
      hphSum[0][i]->SetAxisRange(0,500,"X");
      hphSum[0][i]->Draw();

      hphSum[1][i]->SetLineStyle(1);
      hphSum[1][i]->SetLineColor(kRed);
      hphSum[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cphSum->cd(1);
    TLegend *pbpbphSum = myLegend(0.77,0.65,0.9,0.9);
    pbpbphSum->AddEntry(hphSum[0][0],"Data","p");
    pbpbphSum->AddEntry(hphSum[1][0],"MC","l");
    pbpbphSum->SetTextSize(0.04);
    pbpbphSum->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cphSum->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_phSum_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

        
    // make plot for neMax
    TCanvas *cneMax = new TCanvas("cneMax","",1000,800);
    makeMultiPanelCanvas(cneMax,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cneMax->cd(nbins_cent-i);
      cneMax->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hneMax[0][i]," ","neMax p_{T} (GeV/c)","Event Fraction");
      hneMax[0][i]->SetMarkerStyle(24);
      hneMax[0][i]->SetMarkerColor(kBlack);
      hneMax[0][i]->SetAxisRange(0,500,"X");
      hneMax[0][i]->Draw();

      hneMax[1][i]->SetLineStyle(1);
      hneMax[1][i]->SetLineColor(kRed);
      hneMax[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cneMax->cd(1);
    TLegend *pbpbneMax = myLegend(0.77,0.65,0.9,0.9);
    pbpbneMax->AddEntry(hneMax[0][0],"Data","p");
    pbpbneMax->AddEntry(hneMax[1][0],"MC","l");
    pbpbneMax->SetTextSize(0.04);
    pbpbneMax->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cneMax->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_neMax_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    // make plot for neSum
    TCanvas *cneSum = new TCanvas("cneSum","",1000,800);
    makeMultiPanelCanvas(cneSum,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cneSum->cd(nbins_cent-i);
      cneSum->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hneSum[0][i]," ","neSum p_{T} (GeV/c)","Event Fraction");
      hneSum[0][i]->SetMarkerStyle(24);
      hneSum[0][i]->SetMarkerColor(kBlack);
      hneSum[0][i]->SetAxisRange(0,500,"X");
      hneSum[0][i]->Draw();

      hneSum[1][i]->SetLineStyle(1);
      hneSum[1][i]->SetLineColor(kRed);
      hneSum[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cneSum->cd(1);
    TLegend *pbpbneSum = myLegend(0.77,0.65,0.9,0.9);
    pbpbneSum->AddEntry(hneSum[0][0],"Data","p");
    pbpbneSum->AddEntry(hneSum[1][0],"MC","l");
    pbpbneSum->SetTextSize(0.04);
    pbpbneSum->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cneSum->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_neSum_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");


        
    // make plot for muMax
    TCanvas *cmuMax = new TCanvas("cmuMax","",1000,800);
    makeMultiPanelCanvas(cmuMax,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cmuMax->cd(nbins_cent-i);
      cmuMax->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hmuMax[0][i]," ","muMax p_{T} (GeV/c)","Event Fraction");
      hmuMax[0][i]->SetMarkerStyle(24);
      hmuMax[0][i]->SetMarkerColor(kBlack);
      hmuMax[0][i]->SetAxisRange(0,500,"X");
      hmuMax[0][i]->Draw();

      hmuMax[1][i]->SetLineStyle(1);
      hmuMax[1][i]->SetLineColor(kRed);
      hmuMax[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cmuMax->cd(1);
    TLegend *pbpbmuMax = myLegend(0.77,0.65,0.9,0.9);
    pbpbmuMax->AddEntry(hmuMax[0][0],"Data","p");
    pbpbmuMax->AddEntry(hmuMax[1][0],"MC","l");
    pbpbmuMax->SetTextSize(0.04);
    pbpbmuMax->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cmuMax->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_muMax_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    // make plot for muSum
    TCanvas *cmuSum = new TCanvas("cmuSum","",1000,800);
    makeMultiPanelCanvas(cmuSum,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cmuSum->cd(nbins_cent-i);
      cmuSum->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hmuSum[0][i]," ","muSum p_{T} (GeV/c)","Event Fraction");
      hmuSum[0][i]->SetMarkerStyle(24);
      hmuSum[0][i]->SetMarkerColor(kBlack);
      hmuSum[0][i]->SetAxisRange(0,500,"X");
      hmuSum[0][i]->Draw();

      hmuSum[1][i]->SetLineStyle(1);
      hmuSum[1][i]->SetLineColor(kRed);
      hmuSum[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cmuSum->cd(1);
    TLegend *pbpbmuSum = myLegend(0.77,0.65,0.9,0.9);
    pbpbmuSum->AddEntry(hmuSum[0][0],"Data","p");
    pbpbmuSum->AddEntry(hmuSum[1][0],"MC","l");
    pbpbmuSum->SetTextSize(0.04);
    pbpbmuSum->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    cmuSum->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_muSum_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");


   // make plot for eMax
    TCanvas *ceMax = new TCanvas("ceMax","",1000,800);
    makeMultiPanelCanvas(ceMax,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      ceMax->cd(nbins_cent-i);
      ceMax->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(heMax[0][i]," ","eMax p_{T} (GeV/c)","Event Fraction");
      heMax[0][i]->SetMarkerStyle(24);
      heMax[0][i]->SetMarkerColor(kBlack);
      heMax[0][i]->SetAxisRange(0,500,"X");
      heMax[0][i]->Draw();

      heMax[1][i]->SetLineStyle(1);
      heMax[1][i]->SetLineColor(kRed);
      heMax[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    ceMax->cd(1);
    TLegend *pbpbeMax = myLegend(0.77,0.65,0.9,0.9);
    pbpbeMax->AddEntry(heMax[0][0],"Data","p");
    pbpbeMax->AddEntry(heMax[1][0],"MC","l");
    pbpbeMax->SetTextSize(0.04);
    pbpbeMax->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    ceMax->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_eMax_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    // make plot for eSum
    TCanvas *ceSum = new TCanvas("ceSum","",1000,800);
    makeMultiPanelCanvas(ceSum,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      ceSum->cd(nbins_cent-i);
      ceSum->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(heSum[0][i]," ","eSum p_{T} (GeV/c)","Event Fraction");
      heSum[0][i]->SetMarkerStyle(24);
      heSum[0][i]->SetMarkerColor(kBlack);
      heSum[0][i]->SetAxisRange(0,500,"X");
      heSum[0][i]->Draw();

      heSum[1][i]->SetLineStyle(1);
      heSum[1][i]->SetLineColor(kRed);
      heSum[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    ceSum->cd(1);
    TLegend *pbpbeSum = myLegend(0.77,0.65,0.9,0.9);
    pbpbeSum->AddEntry(heSum[0][0],"Data","p");
    pbpbeSum->AddEntry(heSum[1][0],"MC","l");
    pbpbeSum->SetTextSize(0.04);
    pbpbeSum->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.4,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.5,0.33,16);
    drawText("pCES, HBHE(data)",0.5,0.43,16);
    ceSum->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_eSum_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    
    /*
    // make plot for jtpu
    makeMultiPanelCanvas(cjtpu,3,2,0.0,0.0,0.2,0.15,0.07);
    for(int i = 0;i<nbins_cent;i++){

      cjtpu->cd(nbins_cent-i);
      cjtpu->cd(nbins_cent-i)->SetLogy();

      makeHistTitle(hmuSum[0][i]," ","jtpu p_{T} (GeV/c)","Event Fraction");
      hjtpu[0][i]->SetMarkerStyle(20);
      hjtpu[0][i]->SetMarkerColor(kBlack);
      hjtpu[0][i]->Draw();

      hjtpu[1][i]->SetLineStyle(1);
      hjtpu[1][i]->SetLineColor(kRed);
      hjtpu[1][i]->Draw("same L");
      
      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cmuSum->cd(1);
    TLegend *pbpbjtpu = myLegend(0.77,0.65,0.9,0.9);
    pbpbjtpu->AddEntry(hjtpu[0][0],"Data","p");
    pbpbjtpu->AddEntry(hjtpu[1][0],"MC","l");
    pbpbjtpu->SetTextSize(0.04);
    pbpbjtpu->Draw();
    putCMSPrel();    
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);
    cjtpu->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetVariables_jtpu_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    */

#endif 
  }

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  

  if(doJetID){

    cout<<endl<<"Plotting JetID histograms"<<endl;

    //get the jetID trees from the MC and data files 
    TTree *jetID_data = (TTree*)fDatain->Get("jets_ID");
    TTree *jetID_mc = (TTree*)fMCin->Get("jets_ID");

    //declare histogram arrays
    TH1F* hCut1_Data[nbins_cent];
    TH1F* hCut1_MC[nbins_cent]; 
    TH1F* hCut2_Data[nbins_cent];
    TH1F* hCut2_MC[nbins_cent]; 
    TH1F* hCut3_Data[nbins_cent];
    TH1F* hCut3_MC[nbins_cent]; 
    TH1F* hCut4_Data[nbins_cent];
    TH1F* hCut4_MC[nbins_cent]; 
    TH1F* hCut5_Data[nbins_cent];
    TH1F* hCut5_MC[nbins_cent]; 

    for(int i = 0;i<nbins_cent;i++){
      cout<<"cent = "<<i<<endl;
      //declare the histograms 
      hCut1_Data[i] = new TH1F(Form("hCut1_Data_cent%d",i),"",21,0,1.05);
      hCut1_MC[i] = new TH1F(Form("hCut1_MC_cent%d",i),"",21,0,1.05);
      hCut2_Data[i] = new TH1F(Form("hCut2_Data_cent%d",i),"",21,0,1.05);
      hCut2_MC[i] = new TH1F(Form("hCut2_MC_cent%d",i),"",21,0,1.05);
      hCut3_Data[i] = new TH1F(Form("hCut3_Data_cent%d",i),"",100,0,5);
      hCut3_MC[i] = new TH1F(Form("hCut3_MC_cent%d",i),"",100,0,5);
      hCut4_Data[i] = new TH1F(Form("hCut4_Data_cent%d",i),"",100,0,15);
      hCut4_MC[i] = new TH1F(Form("hCut4_MC_cent%d",i),"",100,0,15);
      hCut5_Data[i] = new TH1F(Form("hCut5_Data_cent%d",i),"",21,0,1.05);
      hCut5_MC[i] = new TH1F(Form("hCut5_MC_cent%d",i),"",21,0,1.05);

      //get the histograms from the trees 
      jetID_data->Draw(Form("cut1>>hCut1_Data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("cut1>>hCut1_MC_cent%d",i),Form("cent==%d",i),"goff");
      hCut1_Data[i]->Print("base");
      hCut1_MC[i]->Print("base");
      jetID_data->Draw(Form("cut2>>hCut2_Data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("cut2>>hCut2_MC_cent%d",i),Form("cent==%d",i),"goff");
      hCut2_Data[i]->Print("base");
      hCut2_MC[i]->Print("base");
      jetID_data->Draw(Form("cut3>>hCut3_Data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("cut3>>hCut3_MC_cent%d",i),Form("cent==%d",i),"goff");
      hCut3_Data[i]->Print("base");
      hCut3_MC[i]->Print("base");
      jetID_data->Draw(Form("cut4>>hCut4_Data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("cut4>>hCut4_MC_cent%d",i),Form("cent==%d",i),"goff");
      hCut4_Data[i]->Print("base");
      hCut4_MC[i]->Print("base");
      jetID_data->Draw(Form("cut5>>hCut5_Data_cent%d",i),Form("cent==%d",i),"goff");
      jetID_mc->Draw(Form("cut5>>hCut5_MC_cent%d",i),Form("cent==%d",i),"goff");      
      hCut5_Data[i]->Print("base");
      hCut5_MC[i]->Print("base");

      //normalize w.r.t no of jets for each. (should ideally normalize to the number of events in the dataset)
      hCut1_Data[i]->Scale(1./hCut1_Data[i]->Integral());
      hCut1_MC[i]->Scale(1./hCut1_MC[i]->Integral());
      hCut2_Data[i]->Scale(1./hCut2_Data[i]->Integral());
      hCut2_MC[i]->Scale(1./hCut2_MC[i]->Integral());
      hCut3_Data[i]->Scale(1./hCut3_Data[i]->Integral());
      hCut3_MC[i]->Scale(1./hCut3_MC[i]->Integral());
      hCut4_Data[i]->Scale(1./hCut4_Data[i]->Integral());
      hCut4_MC[i]->Scale(1./hCut4_MC[i]->Integral());
      hCut5_Data[i]->Scale(1./hCut5_Data[i]->Integral());
      hCut5_MC[i]->Scale(1./hCut5_MC[i]->Integral());

      //hCut1_Data[i]->Print("base");
      //hCut1_MC[i]->Print("base");
    
    }
    TCanvas *cCut1 = new TCanvas("cCut1","",1200,800);
    makeMultiPanelCanvas(cCut1,3,2,0.0,0.0,0.2,0.15,0.07);
    
    for(int i = 0;i<nbins_cent;i++){
      cCut1->cd(nbins_cent-i);
      cCut1->cd(nbins_cent-i)->SetLogy();
      
      makeHistTitle(hCut1_MC[i]," ","#frac{chMax}{Jet p_{T}}","Event Fraction");
      hCut1_MC[i]->SetFillStyle(3005);
      hCut1_MC[i]->SetFillColor(kCyan);
      hCut1_MC[i]->Draw("L");      

      hCut1_Data[i]->SetMarkerStyle(20);
      hCut1_Data[i]->SetMarkerColor(kBlack);
      hCut1_Data[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cCut1->cd(1);
    TLegend *pbpbCut1 = myLegend(0.77,0.65,0.9,0.9);
    pbpbCut1->AddEntry(hCut1_Data[0],"Data","p");
    pbpbCut1->AddEntry(hCut1_MC[0],"MC","f");
    pbpbCut1->SetTextSize(0.04);
    pbpbCut1->Draw();
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);

    cCut1->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut1_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    TCanvas *cCut2 = new TCanvas("cCut2","",1200,800);
    makeMultiPanelCanvas(cCut2,3,2,0.0,0.0,0.2,0.15,0.07);
    
    for(int i = 0;i<nbins_cent;i++){
      cCut2->cd(nbins_cent-i);
      cCut2->cd(nbins_cent-i)->SetLogy();
      
      makeHistTitle(hCut2_MC[i]," ","#frac{neMax}{Maximum(chSum,neSum)}","Event Fraction");
      hCut2_MC[i]->SetFillStyle(3005);
      hCut2_MC[i]->SetFillColor(kCyan);
      hCut2_MC[i]->Draw("L");      

      hCut2_Data[i]->SetMarkerStyle(20);
      hCut2_Data[i]->SetMarkerColor(kBlack);
      hCut2_Data[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cCut2->cd(1);
    TLegend *pbpbCut2 = myLegend(0.77,0.65,0.9,0.9);
    pbpbCut2->AddEntry(hCut2_Data[0],"Data","p");
    pbpbCut2->AddEntry(hCut2_MC[0],"MC","f");
    pbpbCut2->SetTextSize(0.04);
    pbpbCut2->Draw();
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);

    cCut2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut2_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    TCanvas *cCut3 = new TCanvas("cCut3","",1200,800);
    makeMultiPanelCanvas(cCut3,3,2,0.0,0.0,0.2,0.15,0.07);
    
    for(int i = 0;i<nbins_cent;i++){
      cCut3->cd(nbins_cent-i);
      cCut3->cd(nbins_cent-i)->SetLogy();
      
      makeHistTitle(hCut3_MC[i]," ","#frac{chSum+phsum+neSum+muSum+eSum–jtpu}{Jet p_{T}}","Event Fraction");
      hCut3_MC[i]->SetFillStyle(3005);
      hCut3_MC[i]->SetFillColor(kCyan);
      hCut3_MC[i]->Draw("L");      

      hCut3_Data[i]->SetMarkerStyle(20);
      hCut3_Data[i]->SetMarkerColor(kBlack);
      hCut3_Data[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cCut3->cd(1);
    TLegend *pbpbCut3 = myLegend(0.77,0.65,0.9,0.9);
    pbpbCut3->AddEntry(hCut3_Data[0],"Data","p");
    pbpbCut3->AddEntry(hCut3_MC[0],"MC","f");
    pbpbCut3->SetTextSize(0.04);
    pbpbCut3->Draw();
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);

    cCut3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut3_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    TCanvas *cCut4 = new TCanvas("cCut4","",1200,800);
    makeMultiPanelCanvas(cCut4,3,2,0.0,0.0,0.2,0.15,0.07);
    
    for(int i = 0;i<nbins_cent;i++){
      cCut4->cd(nbins_cent-i);
      cCut4->cd(nbins_cent-i)->SetLogy();
      
      makeHistTitle(hCut4_MC[i]," ","#frac{chSum+phsum+neSum+muSum+eSum}{0.5*raw p_{T}}","Event Fraction");
      hCut4_MC[i]->SetFillStyle(3005);
      hCut4_MC[i]->SetFillColor(kCyan);
      hCut4_MC[i]->Draw("L");      

      hCut4_Data[i]->SetMarkerStyle(20);
      hCut4_Data[i]->SetMarkerColor(kBlack);
      hCut4_Data[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cCut4->cd(1);
    TLegend *pbpbCut4 = myLegend(0.77,0.65,0.9,0.9);
    pbpbCut4->AddEntry(hCut4_Data[0],"Data","p");
    pbpbCut4->AddEntry(hCut4_MC[0],"MC","f");
    pbpbCut4->SetTextSize(0.04);
    pbpbCut4->Draw();
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);

    cCut4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut4_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

    TCanvas *cCut5 = new TCanvas("cCut5","",1200,800);
    makeMultiPanelCanvas(cCut5,3,2,0.0,0.0,0.2,0.15,0.07);
    
    for(int i = 0;i<nbins_cent;i++){
      cCut5->cd(nbins_cent-i);
      cCut5->cd(nbins_cent-i)->SetLogy();
      
      makeHistTitle(hCut5_MC[i]," ","#frac{neMax}{neMax+chMax+phMax}","Event Fraction");
      hCut5_MC[i]->SetFillStyle(3005);
      hCut5_MC[i]->SetFillColor(kCyan);
      hCut5_MC[i]->Draw("L");      

      hCut5_Data[i]->SetMarkerStyle(20);
      hCut5_Data[i]->SetMarkerColor(kBlack);
      hCut5_Data[i]->Draw("same");

      drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.55,0.83,16);

    }
    cCut5->cd(1);
    TLegend *pbpbCut5 = myLegend(0.77,0.65,0.9,0.9);
    pbpbCut5->AddEntry(hCut5_Data[0],"Data","p");
    pbpbCut5->AddEntry(hCut5_MC[0],"MC","f");
    pbpbCut5->SetTextSize(0.04);
    pbpbCut5->Draw();
    putCMSPrel();
    drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo,jet_type, radius),0.3,0.23,16);
    drawText("|#eta|<2, |vz|<15",0.3,0.33,16);
    drawText("pCES, HBHE",0.3,0.43,16);

    cCut5->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut5_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");

  }//do jetID


  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  

  if(do2DCutpTplot){

    //get the jetID trees from the MC and data files 
    TTree *jetID_data = (TTree*)fDatain->Get("jets_ID");
    TTree *jetID_mc = (TTree*)fMCin->Get("jets_ID");

    //each cut gets a canvas. the way these plots are going to be plotted as 4 (Jet80, Jet65, Jet55 and MC) x 6 centrality 
    TH2F *hCut1_pt[4][nbins_cent], *hCut2_pt[4][nbins_cent], *hCut3_pt[4][nbins_cent] ,*hCut4_pt[4][nbins_cent], *hCut5_pt[4][nbins_cent];
    
    for(int i = 0;i<nbins_cent;i++){
      for(int j = 0;j<4;j++){
	hCut1_pt[i][j] = new TH2F(Form("hCut1_pt_%d_cent%d",j,i),"",40,0,10,400,0,400);
	hCut2_pt[i][j] = new TH2F(Form("hCut2_pt_%d_cent%d",j,i),"",40,0,10,400,0,400);
	hCut3_pt[i][j] = new TH2F(Form("hCut3_pt_%d_cent%d",j,i),"",40,0,10,400,0,400);
	hCut4_pt[i][j] = new TH2F(Form("hCut4_pt_%d_cent%d",j,i),"",40,0,10,400,0,400);
	hCut5_pt[i][j] = new TH2F(Form("hCut5_pt_%d_cent%d",j,i),"",40,0,10,400,0,400);
      }
    }
    
    if(do2DCut1){

      jetID_data->Draw("(cut1:pt)>>hCut1_pt_0_cent0","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_0_cent1","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_0_cent2","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_0_cent3","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_0_cent4","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_0_cent5","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==5 && pt<500","goff");    

      jetID_data->Draw("(cut1:pt)>>hCut1_pt_1_cent0","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_1_cent1","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_1_cent2","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_1_cent3","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_1_cent4","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_1_cent5","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut1:pt)>>hCut1_pt_2_cent0","jet80 && l1sl52 && trgObjpt>80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_2_cent1","jet80 && l1sl52 && trgObjpt>80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_2_cent2","jet80 && l1sl52 && trgObjpt>80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_2_cent3","jet80 && l1sl52 && trgObjpt>80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_2_cent4","jet80 && l1sl52 && trgObjpt>80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut1:pt)>>hCut1_pt_2_cent5","jet80 && l1sl52 && trgObjpt>80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut1:jtpt)>>hCut1_pt_2_cent0","refpt<500 && cent==0 && subid==0","goff");
      jetID_data->Draw("(cut1:jtpt)>>hCut1_pt_2_cent1","refpt<500 && cent==1 && subid==0","goff");
      jetID_data->Draw("(cut1:jtpt)>>hCut1_pt_2_cent2","refpt<500 && cent==2 && subid==0","goff");
      jetID_data->Draw("(cut1:jtpt)>>hCut1_pt_2_cent3","refpt<500 && cent==3 && subid==0","goff");
      jetID_data->Draw("(cut1:jtpt)>>hCut1_pt_2_cent4","refpt<500 && cent==4 && subid==0","goff");
      jetID_data->Draw("(cut1:jtpt)>>hCut1_pt_2_cent5","refpt<500 && cent==5 && subid==0","goff");

      TCanvas *cCut1_pt = new TCanvas("cCut1_pt","",1000,800);
      makeMultiPanelCanvas(cCut1_pt,4,6,0.0,0.0,0.2,0.15,0.07);
      cCut1_pt->cd(1);cCut1_pt->cd(1)->SetLogy();hCut1_pt[0][0]->Draw("colz");
      drawText("0-5%",0.5,0.6,15);
      drawText("Jet55",0.2,0.8,15);
      putCMSPrel();
      cCut1_pt->cd(2);cCut1_pt->cd(2)->SetLogy();hCut1_pt[0][1]->Draw("colz");
      drawText("5-10%",0.5,0.6,15);
      cCut1_pt->cd(3);cCut1_pt->cd(3)->SetLogy();hCut1_pt[0][2]->Draw("colz");
      drawText("10-30%",0.5,0.6,15);
      cCut1_pt->cd(4);cCut1_pt->cd(4)->SetLogy();hCut1_pt[0][3]->Draw("colz");
      drawText("30-50%",0.5,0.6,15);
      cCut1_pt->cd(5);cCut1_pt->cd(5)->SetLogy();hCut1_pt[0][4]->Draw("colz");
      drawText("50-70%",0.5,0.6,15);
      cCut1_pt->cd(6);cCut1_pt->cd(6)->SetLogy();hCut1_pt[0][5]->Draw("colz");
      drawText("70-90%",0.5,0.6,15);
    
      cCut1_pt->cd(7);cCut1_pt->cd(7)->SetLogy();hCut1_pt[1][0]->Draw("colz");
      drawText("Jet65",0.2,0.8,15);
      cCut1_pt->cd(8);cCut1_pt->cd(8)->SetLogy();hCut1_pt[1][1]->Draw("colz");
      cCut1_pt->cd(9);cCut1_pt->cd(9)->SetLogy();hCut1_pt[1][2]->Draw("colz");
      cCut1_pt->cd(10);cCut1_pt->cd(10)->SetLogy();hCut1_pt[1][3]->Draw("colz");
      cCut1_pt->cd(11);cCut1_pt->cd(11)->SetLogy();hCut1_pt[1][4]->Draw("colz");
      cCut1_pt->cd(12);cCut1_pt->cd(12)->SetLogy();hCut1_pt[1][5]->Draw("colz");
    
      cCut1_pt->cd(13);cCut1_pt->cd(13)->SetLogy();hCut1_pt[2][0]->Draw("colz");
      drawText("Jet80",0.2,0.8,15);
      cCut1_pt->cd(14);cCut1_pt->cd(14)->SetLogy();hCut1_pt[2][1]->Draw("colz");
      cCut1_pt->cd(15);cCut1_pt->cd(15)->SetLogy();hCut1_pt[2][2]->Draw("colz");
      cCut1_pt->cd(16);cCut1_pt->cd(16)->SetLogy();hCut1_pt[2][3]->Draw("colz");
      cCut1_pt->cd(17);cCut1_pt->cd(17)->SetLogy();hCut1_pt[2][4]->Draw("colz");
      cCut1_pt->cd(18);cCut1_pt->cd(18)->SetLogy();hCut1_pt[2][5]->Draw("colz");
    
      cCut1_pt->cd(19);cCut1_pt->cd(19)->SetLogy();hCut1_pt[3][0]->Draw("colz");
      drawText("MC",0.2,0.8,15);
      cCut1_pt->cd(20);cCut1_pt->cd(20)->SetLogy();hCut1_pt[3][1]->Draw("colz");
      cCut1_pt->cd(21);cCut1_pt->cd(21)->SetLogy();hCut1_pt[3][2]->Draw("colz");
      cCut1_pt->cd(22);cCut1_pt->cd(22)->SetLogy();hCut1_pt[3][3]->Draw("colz");
      cCut1_pt->cd(23);cCut1_pt->cd(23)->SetLogy();hCut1_pt[3][4]->Draw("colz");
      cCut1_pt->cd(24);cCut1_pt->cd(24)->SetLogy();hCut1_pt[3][5]->Draw("colz");
    
      cCut1_pt->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut1_pt_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    }


    if(do2DCut2){    

      jetID_data->Draw("(cut2:pt)>>hCut2_pt_0_cent0","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_0_cent1","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_0_cent2","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_0_cent3","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_0_cent4","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_0_cent5","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==5 && pt<500","goff");    

      jetID_data->Draw("(cut2:pt)>>hCut2_pt_1_cent0","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_1_cent1","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_1_cent2","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_1_cent3","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_1_cent4","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_1_cent5","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut2:pt)>>hCut2_pt_2_cent0","jet80 && l1sl52 && trgObjpt>80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_2_cent1","jet80 && l1sl52 && trgObjpt>80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_2_cent2","jet80 && l1sl52 && trgObjpt>80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_2_cent3","jet80 && l1sl52 && trgObjpt>80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_2_cent4","jet80 && l1sl52 && trgObjpt>80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut2:pt)>>hCut2_pt_2_cent5","jet80 && l1sl52 && trgObjpt>80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut2:jtpt)>>hCut2_pt_2_cent0","refpt<500 && cent==0 && subid==0","goff");
      jetID_data->Draw("(cut2:jtpt)>>hCut2_pt_2_cent1","refpt<500 && cent==1 && subid==0","goff");
      jetID_data->Draw("(cut2:jtpt)>>hCut2_pt_2_cent2","refpt<500 && cent==2 && subid==0","goff");
      jetID_data->Draw("(cut2:jtpt)>>hCut2_pt_2_cent3","refpt<500 && cent==3 && subid==0","goff");
      jetID_data->Draw("(cut2:jtpt)>>hCut2_pt_2_cent4","refpt<500 && cent==4 && subid==0","goff");
      jetID_data->Draw("(cut2:jtpt)>>hCut2_pt_2_cent5","refpt<500 && cent==5 && subid==0","goff");

      TCanvas *cCut2_pt = new TCanvas("cCut2_pt","",1000,800);
      makeMultiPanelCanvas(cCut2_pt,4,6,0.0,0.0,0.2,0.15,0.07);
      cCut2_pt->cd(1);cCut2_pt->cd(1)->SetLogy();hCut2_pt[0][0]->Draw("colz");
      drawText("0-5%",0.5,0.6,15);
      drawText("Jet55",0.2,0.8,15);
      putCMSPrel();
      cCut2_pt->cd(2);cCut2_pt->cd(2)->SetLogy();hCut2_pt[0][1]->Draw("colz");
      drawText("5-10%",0.5,0.6,15);
      cCut2_pt->cd(3);cCut2_pt->cd(3)->SetLogy();hCut2_pt[0][2]->Draw("colz");
      drawText("10-30%",0.5,0.6,15);
      cCut2_pt->cd(4);cCut2_pt->cd(4)->SetLogy();hCut2_pt[0][3]->Draw("colz");
      drawText("30-50%",0.5,0.6,15);
      cCut2_pt->cd(5);cCut2_pt->cd(5)->SetLogy();hCut2_pt[0][4]->Draw("colz");
      drawText("50-70%",0.5,0.6,15);
      cCut2_pt->cd(6);cCut2_pt->cd(6)->SetLogy();hCut2_pt[0][5]->Draw("colz");
      drawText("70-90%",0.5,0.6,15);
    
      cCut2_pt->cd(7);cCut2_pt->cd(7)->SetLogy();hCut2_pt[1][0]->Draw("colz");
      drawText("Jet65",0.2,0.8,15);
      cCut2_pt->cd(8);cCut2_pt->cd(8)->SetLogy();hCut2_pt[1][1]->Draw("colz");
      cCut2_pt->cd(9);cCut2_pt->cd(9)->SetLogy();hCut2_pt[1][2]->Draw("colz");
      cCut2_pt->cd(10);cCut2_pt->cd(10)->SetLogy();hCut2_pt[1][3]->Draw("colz");
      cCut2_pt->cd(11);cCut2_pt->cd(11)->SetLogy();hCut2_pt[1][4]->Draw("colz");
      cCut2_pt->cd(12);cCut2_pt->cd(12)->SetLogy();hCut2_pt[1][5]->Draw("colz");
    
      cCut2_pt->cd(13);cCut2_pt->cd(13)->SetLogy();hCut2_pt[2][0]->Draw("colz");
      drawText("Jet80",0.2,0.8,15);
      cCut2_pt->cd(14);cCut2_pt->cd(14)->SetLogy();hCut2_pt[2][1]->Draw("colz");
      cCut2_pt->cd(15);cCut2_pt->cd(15)->SetLogy();hCut2_pt[2][2]->Draw("colz");
      cCut2_pt->cd(16);cCut2_pt->cd(16)->SetLogy();hCut2_pt[2][3]->Draw("colz");
      cCut2_pt->cd(17);cCut2_pt->cd(17)->SetLogy();hCut2_pt[2][4]->Draw("colz");
      cCut2_pt->cd(18);cCut2_pt->cd(18)->SetLogy();hCut2_pt[2][5]->Draw("colz");
    
      cCut2_pt->cd(19);cCut2_pt->cd(19)->SetLogy();hCut2_pt[3][0]->Draw("colz");
      drawText("MC",0.2,0.8,15);
      cCut2_pt->cd(20);cCut2_pt->cd(20)->SetLogy();hCut2_pt[3][1]->Draw("colz");
      cCut2_pt->cd(21);cCut2_pt->cd(21)->SetLogy();hCut2_pt[3][2]->Draw("colz");
      cCut2_pt->cd(22);cCut2_pt->cd(22)->SetLogy();hCut2_pt[3][3]->Draw("colz");
      cCut2_pt->cd(23);cCut2_pt->cd(23)->SetLogy();hCut2_pt[3][4]->Draw("colz");
      cCut2_pt->cd(24);cCut2_pt->cd(24)->SetLogy();hCut2_pt[3][5]->Draw("colz");
    
      cCut2_pt->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut2_pt_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    }
    
    if(do2DCut3){

      jetID_data->Draw("(cut3:pt)>>hCut3_pt_0_cent0","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_0_cent1","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_0_cent2","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_0_cent3","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_0_cent4","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_0_cent5","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==5 && pt<500","goff");    

      jetID_data->Draw("(cut3:pt)>>hCut3_pt_1_cent0","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_1_cent1","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_1_cent2","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_1_cent3","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_1_cent4","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_1_cent5","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut3:pt)>>hCut3_pt_2_cent0","jet80 && l1sl52 && trgObjpt>80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_2_cent1","jet80 && l1sl52 && trgObjpt>80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_2_cent2","jet80 && l1sl52 && trgObjpt>80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_2_cent3","jet80 && l1sl52 && trgObjpt>80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_2_cent4","jet80 && l1sl52 && trgObjpt>80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut3:pt)>>hCut3_pt_2_cent5","jet80 && l1sl52 && trgObjpt>80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut3:jtpt)>>hCut3_pt_2_cent0","refpt<500 && cent==0 && subid==0","goff");
      jetID_data->Draw("(cut3:jtpt)>>hCut3_pt_2_cent1","refpt<500 && cent==1 && subid==0","goff");
      jetID_data->Draw("(cut3:jtpt)>>hCut3_pt_2_cent2","refpt<500 && cent==2 && subid==0","goff");
      jetID_data->Draw("(cut3:jtpt)>>hCut3_pt_2_cent3","refpt<500 && cent==3 && subid==0","goff");
      jetID_data->Draw("(cut3:jtpt)>>hCut3_pt_2_cent4","refpt<500 && cent==4 && subid==0","goff");
      jetID_data->Draw("(cut3:jtpt)>>hCut3_pt_2_cent5","refpt<500 && cent==5 && subid==0","goff");

      TCanvas *cCut3_pt = new TCanvas("cCut3_pt","",1000,800);
      makeMultiPanelCanvas(cCut3_pt,4,6,0.0,0.0,0.2,0.15,0.07);
      cCut3_pt->cd(1);cCut3_pt->cd(1)->SetLogy();hCut3_pt[0][0]->Draw("colz");
      drawText("0-5%",0.5,0.6,15);
      drawText("Jet55",0.2,0.8,15);
      putCMSPrel();
      cCut3_pt->cd(2);cCut3_pt->cd(2)->SetLogy();hCut3_pt[0][1]->Draw("colz");
      drawText("5-10%",0.5,0.6,15);
      cCut3_pt->cd(3);cCut3_pt->cd(3)->SetLogy();hCut3_pt[0][2]->Draw("colz");
      drawText("10-30%",0.5,0.6,15);
      cCut3_pt->cd(4);cCut3_pt->cd(4)->SetLogy();hCut3_pt[0][3]->Draw("colz");
      drawText("30-50%",0.5,0.6,15);
      cCut3_pt->cd(5);cCut3_pt->cd(5)->SetLogy();hCut3_pt[0][4]->Draw("colz");
      drawText("50-70%",0.5,0.6,15);
      cCut3_pt->cd(6);cCut3_pt->cd(6)->SetLogy();hCut3_pt[0][5]->Draw("colz");
      drawText("70-90%",0.5,0.6,15);
    
      cCut3_pt->cd(7);cCut3_pt->cd(7)->SetLogy();hCut3_pt[1][0]->Draw("colz");
      drawText("Jet65",0.2,0.8,15);
      cCut3_pt->cd(8);cCut3_pt->cd(8)->SetLogy();hCut3_pt[1][1]->Draw("colz");
      cCut3_pt->cd(9);cCut3_pt->cd(9)->SetLogy();hCut3_pt[1][2]->Draw("colz");
      cCut3_pt->cd(10);cCut3_pt->cd(10)->SetLogy();hCut3_pt[1][3]->Draw("colz");
      cCut3_pt->cd(11);cCut3_pt->cd(11)->SetLogy();hCut3_pt[1][4]->Draw("colz");
      cCut3_pt->cd(12);cCut3_pt->cd(12)->SetLogy();hCut3_pt[1][5]->Draw("colz");
    
      cCut3_pt->cd(13);cCut3_pt->cd(13)->SetLogy();hCut3_pt[2][0]->Draw("colz");
      drawText("Jet80",0.2,0.8,15);
      cCut3_pt->cd(14);cCut3_pt->cd(14)->SetLogy();hCut3_pt[2][1]->Draw("colz");
      cCut3_pt->cd(15);cCut3_pt->cd(15)->SetLogy();hCut3_pt[2][2]->Draw("colz");
      cCut3_pt->cd(16);cCut3_pt->cd(16)->SetLogy();hCut3_pt[2][3]->Draw("colz");
      cCut3_pt->cd(17);cCut3_pt->cd(17)->SetLogy();hCut3_pt[2][4]->Draw("colz");
      cCut3_pt->cd(18);cCut3_pt->cd(18)->SetLogy();hCut3_pt[2][5]->Draw("colz");
    
      cCut3_pt->cd(19);cCut3_pt->cd(19)->SetLogy();hCut3_pt[3][0]->Draw("colz");
      drawText("MC",0.2,0.8,15);
      cCut3_pt->cd(20);cCut3_pt->cd(20)->SetLogy();hCut3_pt[3][1]->Draw("colz");
      cCut3_pt->cd(21);cCut3_pt->cd(21)->SetLogy();hCut3_pt[3][2]->Draw("colz");
      cCut3_pt->cd(22);cCut3_pt->cd(22)->SetLogy();hCut3_pt[3][3]->Draw("colz");
      cCut3_pt->cd(23);cCut3_pt->cd(23)->SetLogy();hCut3_pt[3][4]->Draw("colz");
      cCut3_pt->cd(24);cCut3_pt->cd(24)->SetLogy();hCut3_pt[3][5]->Draw("colz");
    
      cCut3_pt->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut3_pt_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    }

    if(do2DCut4){
 
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_0_cent0","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_0_cent1","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_0_cent2","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_0_cent3","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_0_cent4","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_0_cent5","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==5 && pt<500","goff");    

      jetID_data->Draw("(cut4:pt)>>hCut4_pt_1_cent0","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_1_cent1","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_1_cent2","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_1_cent3","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_1_cent4","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_1_cent5","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut4:pt)>>hCut4_pt_2_cent0","jet80 && l1sl52 && trgObjpt>80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_2_cent1","jet80 && l1sl52 && trgObjpt>80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_2_cent2","jet80 && l1sl52 && trgObjpt>80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_2_cent3","jet80 && l1sl52 && trgObjpt>80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_2_cent4","jet80 && l1sl52 && trgObjpt>80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut4:pt)>>hCut4_pt_2_cent5","jet80 && l1sl52 && trgObjpt>80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut4:jtpt)>>hCut4_pt_2_cent0","refpt<500 && cent==0 && subid==0","goff");
      jetID_data->Draw("(cut4:jtpt)>>hCut4_pt_2_cent1","refpt<500 && cent==1 && subid==0","goff");
      jetID_data->Draw("(cut4:jtpt)>>hCut4_pt_2_cent2","refpt<500 && cent==2 && subid==0","goff");
      jetID_data->Draw("(cut4:jtpt)>>hCut4_pt_2_cent3","refpt<500 && cent==3 && subid==0","goff");
      jetID_data->Draw("(cut4:jtpt)>>hCut4_pt_2_cent4","refpt<500 && cent==4 && subid==0","goff");
      jetID_data->Draw("(cut4:jtpt)>>hCut4_pt_2_cent5","refpt<500 && cent==5 && subid==0","goff");

      TCanvas *cCut4_pt = new TCanvas("cCut4_pt","",1000,800);
      makeMultiPanelCanvas(cCut4_pt,4,6,0.0,0.0,0.2,0.15,0.07);
      cCut4_pt->cd(1);cCut4_pt->cd(1)->SetLogy();hCut4_pt[0][0]->Draw("colz");
      drawText("0-5%",0.5,0.6,15);
      drawText("Jet55",0.2,0.8,15);
      putCMSPrel();
      cCut4_pt->cd(2);cCut4_pt->cd(2)->SetLogy();hCut4_pt[0][1]->Draw("colz");
      drawText("5-10%",0.5,0.6,15);
      cCut4_pt->cd(3);cCut4_pt->cd(3)->SetLogy();hCut4_pt[0][2]->Draw("colz");
      drawText("10-30%",0.5,0.6,15);
      cCut4_pt->cd(4);cCut4_pt->cd(4)->SetLogy();hCut4_pt[0][3]->Draw("colz");
      drawText("30-50%",0.5,0.6,15);
      cCut4_pt->cd(5);cCut4_pt->cd(5)->SetLogy();hCut4_pt[0][4]->Draw("colz");
      drawText("50-70%",0.5,0.6,15);
      cCut4_pt->cd(6);cCut4_pt->cd(6)->SetLogy();hCut4_pt[0][5]->Draw("colz");
      drawText("70-90%",0.5,0.6,15);
    
      cCut4_pt->cd(7);cCut4_pt->cd(7)->SetLogy();hCut4_pt[1][0]->Draw("colz");
      drawText("Jet65",0.2,0.8,15);
      cCut4_pt->cd(8);cCut4_pt->cd(8)->SetLogy();hCut4_pt[1][1]->Draw("colz");
      cCut4_pt->cd(9);cCut4_pt->cd(9)->SetLogy();hCut4_pt[1][2]->Draw("colz");
      cCut4_pt->cd(10);cCut4_pt->cd(10)->SetLogy();hCut4_pt[1][3]->Draw("colz");
      cCut4_pt->cd(11);cCut4_pt->cd(11)->SetLogy();hCut4_pt[1][4]->Draw("colz");
      cCut4_pt->cd(12);cCut4_pt->cd(12)->SetLogy();hCut4_pt[1][5]->Draw("colz");
    
      cCut4_pt->cd(13);cCut4_pt->cd(13)->SetLogy();hCut4_pt[2][0]->Draw("colz");
      drawText("Jet80",0.2,0.8,15);
      cCut4_pt->cd(14);cCut4_pt->cd(14)->SetLogy();hCut4_pt[2][1]->Draw("colz");
      cCut4_pt->cd(15);cCut4_pt->cd(15)->SetLogy();hCut4_pt[2][2]->Draw("colz");
      cCut4_pt->cd(16);cCut4_pt->cd(16)->SetLogy();hCut4_pt[2][3]->Draw("colz");
      cCut4_pt->cd(17);cCut4_pt->cd(17)->SetLogy();hCut4_pt[2][4]->Draw("colz");
      cCut4_pt->cd(18);cCut4_pt->cd(18)->SetLogy();hCut4_pt[2][5]->Draw("colz");
    
      cCut4_pt->cd(19);cCut4_pt->cd(19)->SetLogy();hCut4_pt[3][0]->Draw("colz");
      drawText("MC",0.2,0.8,15);
      cCut4_pt->cd(20);cCut4_pt->cd(20)->SetLogy();hCut4_pt[3][1]->Draw("colz");
      cCut4_pt->cd(21);cCut4_pt->cd(21)->SetLogy();hCut4_pt[3][2]->Draw("colz");
      cCut4_pt->cd(22);cCut4_pt->cd(22)->SetLogy();hCut4_pt[3][3]->Draw("colz");
      cCut4_pt->cd(23);cCut4_pt->cd(23)->SetLogy();hCut4_pt[3][4]->Draw("colz");
      cCut4_pt->cd(24);cCut4_pt->cd(24)->SetLogy();hCut4_pt[3][5]->Draw("colz");
    
      cCut4_pt->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut4_pt_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    }

    if(do2DCut5){

      jetID_data->Draw("(cut5:pt)>>hCut5_pt_0_cent0","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_0_cent1","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_0_cent2","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_0_cent3","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_0_cent4","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_0_cent5","jet55 && !jet65 && !jet80 && l1sl36 && trgObjpt>55 && trgObjpt<65 && cent==5 && pt<500","goff");    

      jetID_data->Draw("(cut5:pt)>>hCut5_pt_1_cent0","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_1_cent1","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_1_cent2","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_1_cent3","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_1_cent4","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_1_cent5","jet65 && !jet80 && l1sl36 && trgObjpt>65 && trgObjpt<80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut5:pt)>>hCut5_pt_2_cent0","jet80 && l1sl52 && trgObjpt>80 && cent==0 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_2_cent1","jet80 && l1sl52 && trgObjpt>80 && cent==1 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_2_cent2","jet80 && l1sl52 && trgObjpt>80 && cent==2 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_2_cent3","jet80 && l1sl52 && trgObjpt>80 && cent==3 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_2_cent4","jet80 && l1sl52 && trgObjpt>80 && cent==4 && pt<500","goff");    
      jetID_data->Draw("(cut5:pt)>>hCut5_pt_2_cent5","jet80 && l1sl52 && trgObjpt>80 && cent==5 && pt<500","goff"); 

      jetID_data->Draw("(cut5:jtpt)>>hCut5_pt_2_cent0","refpt<500 && cent==0 && subid==0","goff");
      jetID_data->Draw("(cut5:jtpt)>>hCut5_pt_2_cent1","refpt<500 && cent==1 && subid==0","goff");
      jetID_data->Draw("(cut5:jtpt)>>hCut5_pt_2_cent2","refpt<500 && cent==2 && subid==0","goff");
      jetID_data->Draw("(cut5:jtpt)>>hCut5_pt_2_cent3","refpt<500 && cent==3 && subid==0","goff");
      jetID_data->Draw("(cut5:jtpt)>>hCut5_pt_2_cent4","refpt<500 && cent==4 && subid==0","goff");
      jetID_data->Draw("(cut5:jtpt)>>hCut5_pt_2_cent5","refpt<500 && cent==5 && subid==0","goff");

      TCanvas *cCut5_pt = new TCanvas("cCut5_pt","",1000,800);
      makeMultiPanelCanvas(cCut5_pt,4,6,0.0,0.0,0.2,0.15,0.07);
      cCut5_pt->cd(1);cCut5_pt->cd(1)->SetLogy();hCut5_pt[0][0]->Draw("colz");
      drawText("0-5%",0.5,0.6,15);
      drawText("Jet55",0.2,0.8,15);
      putCMSPrel();
      cCut5_pt->cd(2);cCut5_pt->cd(2)->SetLogy();hCut5_pt[0][1]->Draw("colz");
      drawText("5-10%",0.5,0.6,15);
      cCut5_pt->cd(3);cCut5_pt->cd(3)->SetLogy();hCut5_pt[0][2]->Draw("colz");
      drawText("10-30%",0.5,0.6,15);
      cCut5_pt->cd(4);cCut5_pt->cd(4)->SetLogy();hCut5_pt[0][3]->Draw("colz");
      drawText("30-50%",0.5,0.6,15);
      cCut5_pt->cd(5);cCut5_pt->cd(5)->SetLogy();hCut5_pt[0][4]->Draw("colz");
      drawText("50-70%",0.5,0.6,15);
      cCut5_pt->cd(6);cCut5_pt->cd(6)->SetLogy();hCut5_pt[0][5]->Draw("colz");
      drawText("70-90%",0.5,0.6,15);
    
      cCut5_pt->cd(7);cCut5_pt->cd(7)->SetLogy();hCut5_pt[1][0]->Draw("colz");
      drawText("Jet65",0.2,0.8,15);
      cCut5_pt->cd(8);cCut5_pt->cd(8)->SetLogy();hCut5_pt[1][1]->Draw("colz");
      cCut5_pt->cd(9);cCut5_pt->cd(9)->SetLogy();hCut5_pt[1][2]->Draw("colz");
      cCut5_pt->cd(10);cCut5_pt->cd(10)->SetLogy();hCut5_pt[1][3]->Draw("colz");
      cCut5_pt->cd(11);cCut5_pt->cd(11)->SetLogy();hCut5_pt[1][4]->Draw("colz");
      cCut5_pt->cd(12);cCut5_pt->cd(12)->SetLogy();hCut5_pt[1][5]->Draw("colz");
    
      cCut5_pt->cd(13);cCut5_pt->cd(13)->SetLogy();hCut5_pt[2][0]->Draw("colz");
      drawText("Jet80",0.2,0.8,15);
      cCut5_pt->cd(14);cCut5_pt->cd(14)->SetLogy();hCut5_pt[2][1]->Draw("colz");
      cCut5_pt->cd(15);cCut5_pt->cd(15)->SetLogy();hCut5_pt[2][2]->Draw("colz");
      cCut5_pt->cd(16);cCut5_pt->cd(16)->SetLogy();hCut5_pt[2][3]->Draw("colz");
      cCut5_pt->cd(17);cCut5_pt->cd(17)->SetLogy();hCut5_pt[2][4]->Draw("colz");
      cCut5_pt->cd(18);cCut5_pt->cd(18)->SetLogy();hCut5_pt[2][5]->Draw("colz");
    
      cCut5_pt->cd(19);cCut5_pt->cd(19)->SetLogy();hCut5_pt[3][0]->Draw("colz");
      drawText("MC",0.2,0.8,15);
      cCut5_pt->cd(20);cCut5_pt->cd(20)->SetLogy();hCut5_pt[3][1]->Draw("colz");
      cCut5_pt->cd(21);cCut5_pt->cd(21)->SetLogy();hCut5_pt[3][2]->Draw("colz");
      cCut5_pt->cd(22);cCut5_pt->cd(22)->SetLogy();hCut5_pt[3][3]->Draw("colz");
      cCut5_pt->cd(23);cCut5_pt->cd(23)->SetLogy();hCut5_pt[3][4]->Draw("colz");
      cCut5_pt->cd(24);cCut5_pt->cd(24)->SetLogy();hCut5_pt[3][5]->Draw("colz");
    
      cCut5_pt->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_JetID_cut5_pt_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
    }
  }

  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;


}























