// Raghav Kunnawalkam Elayavalli
// April 22th 2015
// CERN
// raghav.k.e AT CERN DOT CH

// Macro to perform svd unfolding.


#include <iostream>
#include <stdio.h>

#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/prior.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/bayesianUnfold.h"
#include "TStopwatch.h"


static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};



void RAA_unfold_svd(int radius = 3, char* algo = (char*) "Pu", char *jet_type = (char*) "PF", int unfoldingCut = 30, char* etaWidth = (char*) "n20_eta_p20", double deltaEta = 4.0){

  TStopwatch timer; 
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  bool printDebug = true;

  // get the data and mc histograms from the output of the read macro. 
  
  TDatime date;//this is just here to get them to run optimized. 

  // Raghav's files: 
  TFile * fPbPb_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pawan_ntuple_PbPb_data_MC_spectra_JetID_CutA_analysisbins_%s_R0p%d.root",etaWidth,radius));
  //TFile * fPP_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pp_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_rebinned_eMaxSumcand_A_R0p%d.root",radius));
  // no jet ID cut. 
  TFile * fPP_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pawan_ntuple_PP_data_MC_spectra_residualFactor_analysisbins_%s_R0p%d.root",etaWidth,radius));
  
  
  // Pawan's files:
  // TFile * fPbPb_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pawan_ntuplehistograms/PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));
  // //TFile * fPP_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pp_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));
  // TFile * fPP_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pawan_ntuplehistograms/Pp_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));

  TFile * fPbPb_MB_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_MinBiasUPC_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));

  // we also need to get the files for the MinBias closure histograms.
  TFile * fPbPb_MCclosure_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_MC_closure_histogram_deltaR_0p2_akPu%d_20150423.root",radius));
  TFile * fPP_MCclosure_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_MC_closure_histogram_deltaR_0p2_ak%d_20150423.root",radius));
  //TFile * fPP_MCclosure_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_MC_calo_pf_jet_correlation_mcclosure_test_sameside_deltaR_0p2_ak%d_20150413.root",radius));
  
  cout<<"after input file declaration"<<endl;
  // need to make sure that the file names are in prefect order so that i can run them one after another. 
  // for the above condition, i might have to play with the date stamp. 
  
  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  // histogram declarations with the following initial appendage: d - Data, m - MC, u- Unfolded
  // for the MC closure test, ive kept separate 
  
  // setup the radius and the eta bin loop here later. not for the time being. Aug 20th. only run the -2 < eta < 2 with the differenent centrality bins 
  
  TH1F *dPbPb_TrgComb[nbins_cent];
  TH2F *mPbPb_Matrix[nbins_cent];
  TH1F *uPbPb_SVD[nbins_cent];
  TH1F *mPbPb_mcclosure_data[nbins_cent];
  TH2F *mPbPb_mcclosure_Matrix[nbins_cent];
  TH1F *uPbPb_MC_SVD[nbins_cent];
  
  TH1F *dPP_Comb;
  TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix;
  TH1F *mPP_mcclosure_data;
  TH2F *mPP_mcclosure_Matrix;
  TH1F *uPP_SVD;
  TH1F *uPP_MC_SVD;
  
  // would be better to read in the histograms and rebin them. come to think of it, it would be better to have them already rebinned (and properly scaled - to the level of differential cross section in what ever barns (inverse micro barns) but keep it consistent) from the read macro. 
      
  // get PbPb data
  for(int i = 0;i<nbins_cent;++i){
    if(printDebug) cout<<"cent_"<<i<<endl;
    dPbPb_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i));
    dPbPb_TrgComb[i]->Print("base");
    //dPbPb_TrgComb[i]->Scale(1./1e3/145.156);
    mPbPb_Matrix[i] = (TH2F*)fPbPb_in->Get(Form("hpbpb_matrix_R%d_%s_cent%d",radius,etaWidth,i));
    mPbPb_Matrix[i]->Print("base");
    // mPbPb_mcclosure_data[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_reco_R%d_%s_cent%d",radius,etaWidth,i));
    // mPbPb_mcclosure_data[i]->Print("base");
    // mPbPb_mcclosure_Matrix[i] = (TH2F*)fPbPb_in->Get(Form("hpbpb_matrix_R%d_%s_cent%d",radius,etaWidth,i));
    // mPbPb_mcclosure_Matrix[i]->Print("base");
  }

  dPP_Comb = (TH1F*)fPP_in->Get(Form("hpp_HLTComb_R%d_%s",radius,etaWidth));
  dPP_Comb->Print("base");
  //dPP_Comb->Scale(1./5.3/1e9);
  mPP_Matrix = (TH2F*)fPP_in->Get(Form("hpp_matrix_HLT_R%d_%s",radius,etaWidth));
  mPP_Matrix->Print("base");
  // mPP_mcclosure_data = (TH1F*)fPP_in->Get(Form("hpp_JetComb_reco_R%d_%s",radius,etaWidth));
  // mPP_mcclosure_data->Print("base");
  // mPP_mcclosure_Matrix = (TH2F*)fPP_in->Get(Form("hpp_matrix_HLT_R%d_%s",radius,etaWidth));
  // mPP_mcclosure_Matrix->Print("base");
  
  Int_t nSVDIter = 1;
  
  for(int j = 4; j<=20; ++j){
    nSVDIter = j;

    for(int i = 0;i<nbins_cent;++i){
      // rebin histograms before we send it to unfolding 
      // histograms are already binned
      // need to scale the Data to be in the same units as the prior. MC is in mb and data was in counts. 
      
      //since SVD is very straight forward, lets do it rignt here:
      //get the SVD response matrix:
      RooUnfoldResponse ruResponse(mPbPb_Matrix[i]->ProjectionY(),mPbPb_Matrix[i]->ProjectionX(), mPbPb_Matrix[i],"","");
      //regularization parameter definition: 
      RooUnfoldSvd unfoldSvd(&ruResponse, dPbPb_TrgComb[i], nSVDIter);
      uPbPb_SVD[i] = (TH1F*)unfoldSvd.Hreco();
      uPbPb_SVD[i]->SetName(Form("PbPb_SVD_unfolding_cent%d",i));

      // RooUnfoldResponse ruResponseMC(mPbPb_mcclosure_Matrix[i]->ProjectionY(),mPbPb_mcclosure_Matrix[i]->ProjectionX(), mPbPb_mcclosure_Matrix[i],"","");
      // //regularization parameter definition: 
      // RooUnfoldSvd unfoldSvdMC(&ruResponseMC, mPbPb_mcclosure_data[i], nSVDIter);
      // uPbPb_MC_SVD[i] = (TH1F*)unfoldSvdMC.Hreco();
      // uPbPb_MC_SVD[i]->SetName(Form("PbPb_MC_SVD_unfolding_cent%d",i));
      
    }// centrality bin loop

    // RooUnfoldResponse ruResponseMCPP(mPP_mcclosure_Matrix->ProjectionY(),mPP_mcclosure_Matrix->ProjectionX(), mPP_mcclosure_Matrix,"","");
    // //regularization parameter definition: 
    // RooUnfoldSvd unfoldSvdMCPP(&ruResponseMCPP, mPP_mcclosure_data, nSVDIter);
    // uPP_MC_SVD = (TH1F*)unfoldSvdMCPP.Hreco();
    // uPP_MC_SVD->SetName("PP_MC_SVD_unfolding");
    
    RooUnfoldResponse ruResponsePP(mPP_Matrix->ProjectionY(),mPP_Matrix->ProjectionX(), mPP_Matrix,"","");
    //regularization parameter definition: 
    RooUnfoldSvd unfoldSvdPP(&ruResponsePP, dPP_Comb, nSVDIter);
    uPP_SVD = (TH1F*)unfoldSvdPP.Hreco();
    uPP_SVD->SetName("PP_SVD_unfolding");
  
    TFile fout(Form("../../Output/svd_unfolding_matrix_param%d_%s_R%d_%d.root",nSVDIter,etaWidth,radius,date.GetDate()),"RECREATE");
    for(int i = 0; i<nbins_cent;i++){

      dPbPb_TrgComb[i]->Write();
      uPbPb_SVD[i]->Write();
      mPbPb_Matrix[i]->Write();
      // mPbPb_mcclosure_data[i]->Write();
      // uPbPb_MC_SVD[i]->Write();
      // mPbPb_mcclosure_Matrix[i]->Write();
  
    }

    uPP_SVD->Write();
    dPP_Comb->Write();
    mPP_Matrix->Write();
    // mPP_mcclosure_data->Write();
    // uPP_MC_SVD->Write();
    // mPP_mcclosure_Matrix->Write();
    
    fout.Close();
    
  }// svd unfolding iterations loop
  
  timer.Stop();
  if(printDebug) cout<<"CPU time (mins) = "<<(Float_t)timer.CpuTime()/60<<endl;
  if(printDebug) cout<<"Real tile (mins) = "<<(Float_t)timer.RealTime()/60<<endl;
  

}
