// Raghav Kunnawalkam Elayavalli
// Monday March 25th
// Rutgers

//
// Macro for data driven unfolding error correction since just plain unfolding doesnt fix the error bars correctly.
// here is the method:
//1) For each pT bin in your measured spectra, generate a gaussian distribution with the BinContent as mean and BinError as sigma. 
//2) generate a large number of spectra, say a thousand, by getting a random value from the gaussian for each pt bin. 
//3) now unfold these 1000 spectra, with the same response matrix (derived from recopt-genpt distributions)
//4) Now we have 1000 values of bin contents and bin errors for each pt bin. Fill the bin contents into a histogram and find the mean (which is our BinContent for our final unfolded histogram) and RMS (which the BinError) per pt bin. 
//
//


#include <iostream>
#include <stdio.h>

#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/prior.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/bayesianUnfold.h"
#include "TStopwatch.h"
#include "TRandom3.h"

//static const int nbins_pt = 29;
//static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

static const int nbins_pt = 400;

static const int nbins_pt_coarse = 30;
static const double boundaries_pt_coarse[nbins_pt_coarse+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 300, 
  330, 362, 395
};


// Remove bins with error > central value
void cleanup(TH1F *h){
  for (int i=1;i<=h->GetNbinsX();i++){
    double val1 = h->GetBinContent(i);
    double valErr1 = h->GetBinError(i);
    if (valErr1>=val1) {
      h->SetBinContent(i,0);
      h->SetBinError(i,0);
    }
  }   
}

// Remove error 
void removeError(TH1F *h){
  for (int i=1;i<=h->GetNbinsX();i++){
    h->SetBinError(i,0);
  } 
}  

// Remove Zero
void removeZero(TH1 *h){
  double min = 0;
  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i)>min&&h->GetBinContent(i)>0)
      min = h->GetBinContent(i);
  }
  
  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i) == 0){
      h->SetBinContent(i,min/10.);
      h->SetBinError(i,min/10.);
    }
  }
}

// make a histogram from TF1 function
TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname){
  TH1F *hF = (TH1F*)h->Clone(fHistname);
  for (int i=1;i<=h->GetNbinsX();i++){
    double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
    hF->SetBinContent(i,var);
    hF->SetBinError(i,0);
  }
  return hF;
}

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

void RAA_dataDrivenUnfoldingErrorCheck(int radius = 4, int radiusPP = 4, char* algo = (char*) "Pu", char *jet_type = (char*) "PF", int unfoldingCut = 60, char* etaWidth = (char*) "n20_eta_p20", double deltaEta = 4.0){

  
  TStopwatch timer; 
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  bool printDebug = true;

  // get the data and mc histograms from the output of the read macro. 
  
  TDatime date;//this is just here to get them to run optimized. 

  // //TFile* fMC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_mc_final_chMaxjtpt_ak%s%s_20150202.root",algo,jet_type));
  // //TFile* fMC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_mc_allcutswithouteSum_ak%s%s_20150203.root",algo,jet_type));
  // //TFile* fMC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_mc_nocut_spectra_ak%s%s_20150203.root",algo,jet_type));
  // TFile* fMC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_mc_ak%s%s_20150304.root",algo,jet_type));
  // TFile *fMC_pp_in;
  // //fMC_pp_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_mc_akPF_20150320.root");
  // fMC_pp_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_mc_spectra_chMaxjtpt0p02_eMaxjtpt0p6_ak3PF_20150225.root");
  
  // //TFile *fData_PbPb_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_spectra_histograms_akPuPF_20150128.root");
  // //TFile *fData_PbPb_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_doga_file_spectra_histograms_akPuPF_20150202.root");
  // //TFile *fData_PbPb_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_marguerite_file_spectra_histograms_akPuPF_20150203.root");
  // //TFile *fData_PbPb_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_marguerite_allcutswithouteSum_file_spectra_histograms_akPuPF_20150203.root");
  // TFile *fData_PbPb_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_data_akPuPF_20150304.root");
  // TFile * fData_MinBias_PbPb = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_MinBiasUPC_trigger_turnoncurves_SuperNovaRejected_akPuPF_20150304.root");
  
  // //TFile *fData_PbPb_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_marguerite_nocut_file_spectra_histograms_akPuPF_20150203.root");
  // //TFile *fData_pp_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_data_spectra_trgObj_chMaxjtpt_akPF_20150203.root");
  // //TFile *fData_pp_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_data_spectra_trgObj_allcutwithouteSum_akPF_20150203.root");
  // //TFile *fData_pp_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_data_spectra_trgObj_chMaxjtpt0p02_eMaxjtpt0p6_akPF_20150226.root");
  // TFile *fData_pp_in = TFile::Open("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_data_spectra_trgObj_chMaxjtpt0p02_eMaxjtpt0p6_akPF_20150319.root");
  
  TFile * fPbPb_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));
  TFile * fPP_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/Pp_withresiduals_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));
  TFile * fPbPb_MB_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_MinBiasUPC_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_R0p%d.root",radius));
  
  //TH1F * htest = new TH1F("htest","",nbins_pt, boundaries_pt);
  //Int_t unfoldingCutBin = htest->FindBin(unfoldingCut);
  
  cout<<"after input file declaration"<<endl;
  // need to make sure that the file names are in prefect order so that i can run them one after another. 
  // for the above condition, i might have to play with the date stamp. 
  
  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  // histogram declarations with the following initial appendage: d - Data, m - MC, u- Unfolded
  // for the MC closure test, ive kept separate 

  // setup the radius and the eta bin loop here later. not for the time being. Aug 20th. only run the -2 < eta < 2 with the differenent centrality bins 

  TH1F *dPbPb_TrgComb[nbins_cent+1], *dPbPb_Comb[nbins_cent+1], *dPbPb_Trg80[nbins_cent+1], *dPbPb_Trg65[nbins_cent+1], *dPbPb_Trg55[nbins_cent+1], *dPbPb_1[nbins_cent+1], *dPbPb_2[nbins_cent+1], *dPbPb_3[nbins_cent+1], *dPbPb_80[nbins_cent+1], *dPbPb_65[nbins_cent+1], *dPbPb_55[nbins_cent+1];
  
  TH1F *mPbPb_Gen[nbins_cent+1], *mPbPb_Reco[nbins_cent+1];
  TH2F *mPbPb_Matrix[nbins_cent+1], *mPbPb_Response[nbins_cent+1], *mPbPb_ResponseNorm[nbins_cent+1];
  TH1F *mPbPb_mcclosure_data[nbins_cent+1];
  TH2F *mPbPb_mcclosure_Matrix[nbins_cent+1],*mPbPb_mcclosure_Response[nbins_cent+1], *mPbPb_mcclosure_ResponseNorm[nbins_cent+1];
  TH1F *mPbPb_mcclosure_gen[nbins_cent+1];
  const int Iterations = 20; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1], *uPbPb_BinByBin[nbins_cent+1], *uPbPb_SVD[nbins_cent+1]; 
  TH1F *uPbPb_BayesianIter[nbins_cent+1][Iterations];
  TH1F *dPbPb_MinBias[nbins_cent];
  
  TH1F *dPP_1, *dPP_2, *dPP_3, *dPP_Comb;
  TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix, *mPP_Response,*mPP_ResponseNorm;
  TH1F *mPP_mcclosure_data;
  TH2F *mPP_mcclosure_Matrix, *mPP_mcclosure_Response,*mPP_mcclosure_ResponseNorm;
  TH1F *mPP_mcclosure_Gen;
  TH1F *uPP_Bayes, *uPP_BinByBin, *uPP_SVD;
  TH1F *uPP_BayesianIter[Iterations];

  // would be better to read in the histograms and rebin them. come to think of it, it would be better to have them already rebinned (and properly scaled - to the level of differential cross section in what ever barns (inverse micro barns) but keep it consistent) from the read macro. 

  // get PbPb data
  for(int i = 0;i<nbins_cent;i++){
    if(printDebug) cout<<"cent_"<<i<<endl;
    dPbPb_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLTComb_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_TrgComb[i]->Scale(4*145.156*1e6);
    dPbPb_TrgComb[i]->Print("base");
    dPbPb_Trg80[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT80_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_Trg80[i]->Scale(4*145.156*1e6);
    dPbPb_Trg80[i]->Print("base");
    dPbPb_Trg65[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT65_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_Trg65[i]->Scale(4*145.156*1e6);
    dPbPb_Trg65[i]->Print("base");
    dPbPb_Trg55[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT55_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_Trg55[i]->Scale(4*145.156*1e6);
    dPbPb_Trg55[i]->Print("base");
    //dPbPb_TrgComb[i] = (TH1F*)dPbPb_Trg80[i]->Clone(Form("Jet_80_triggered_spectra_data_PbPb_cent%d",i));
    
    //dPbPb_MinBias[i] = (TH1F*)fPbPb_MB_in->Get(Form("hpbpb_HLTComb_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_MinBias[i]->Print("base");
    dPbPb_TrgComb[i]->Scale(1./(145.156 * 1e9));
    //dPbPb_MinBias[i]->Scale(1./(161.939 * 1e9));
    
    //dPbPb_TrgComb[i]->Add(dPbPb_MinBias[i]);
    
    for(int k = 1;k<=unfoldingCut;k++) {
      dPbPb_TrgComb[i]->SetBinContent(k,0);
      dPbPb_Trg80[i]->SetBinContent(k,0);
      dPbPb_Trg65[i]->SetBinContent(k,0);
      dPbPb_Trg55[i]->SetBinContent(k,0);
    }
    
  }
  
  //Int_t nSVDIter = 4;
  
  if(printDebug)cout<<"loaded the data histograms PbPb"<<endl;
  // get PbPb MC
  for(int i = 0;i<nbins_cent;i++){
    
    mPbPb_Gen[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_JetComb_gen_R%d_n20_eta_p20_cent%d",radius,i));
    //mPbPb_Gen[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_gen_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Gen[i]->Print("base");
    mPbPb_Reco[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_JetComb_reco_R%d_n20_eta_p20_cent%d",radius,i));
    //mPbPb_Reco[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_reco_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Reco[i]->Print("base");
    mPbPb_Matrix[i] = (TH2F*)fPbPb_in->Get(Form("hpbpb_matrix_HLT_R%d_n20_eta_p20_cent%d",radius,i));
    //mPbPb_Matrix[i] = (TH2F*)fPbPb_in->Get(Form("hpbpb_matrix_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Matrix[i]->Print("base");
    mPbPb_mcclosure_data[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_mcclosure_JetComb_data_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_mcclosure_data[i]->Print("base");
    mPbPb_mcclosure_gen[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_mcclosure_gen_JetComb_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_mcclosure_gen[i]->Print("base");
    mPbPb_mcclosure_Matrix[i] = (TH2F*)fPbPb_in->Get(Form("hpbpb_mcclosure_matrix_HLT_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_mcclosure_Matrix[i]->Print("base");
    
    //since SVD is very straight forward, lets do it rignt here:
    //get the SVD response matrix:
    //RooUnfoldResponse ruResponse(mPbPb_Matrix[i]->ProjectionY(),mPbPb_Matrix[i]->ProjectionX(), mPbPb_Matrix[i],"","");
    //regularization parameter definition: 
    //RooUnfoldSvd unfoldSvd(&ruResponse, dPbPb_TrgComb[i], nSVDIter);
    //uPbPb_SVD[i] = (TH1F*)unfoldSvd.Hreco();
  
    
    // for(int k = 1;k<=unfoldingCut;k++){

    //   mPbPb_Gen[i]->SetBinContent(k,0);
    //   mPbPb_Reco[i]->SetBinContent(k,0);
    //   mPbPb_mcclosure_data[i]->SetBinContent(k,0);
    //   mPbPb_mcclosure_gen[i]->SetBinContent(k,0);
    //   for(int l = 1;l<=1000;l++){
    // 	mPbPb_Matrix[i]->SetBinContent(k,l,0);
    // 	mPbPb_mcclosure_Matrix[i]->SetBinContent(k,l,0);
    // 	mPbPb_Matrix[i]->SetBinContent(l,k,0);
    // 	mPbPb_mcclosure_Matrix[i]->SetBinContent(l,k,0);	
    //   }
    // }
    
    //mPbPb_Response[i] = new TH2F(Form("mPbPb_Response_cent%d",i),"Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    //mPbPb_ResponseNorm[i] = new TH2F(Form("mPbPb_ResponseNorm_cent%d",i),"Normalized Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  }
  
  if(printDebug) cout<<"loaded the data and mc PbPb histograms from the files"<<endl;

  // get PP data
  if(printDebug) cout<<"Getting PP data and MC"<<endl;
  dPP_1 = (TH1F*)fPP_in->Get(Form("hpp_HLT80_R%d_%s",radiusPP,etaWidth)); 
  dPP_1->Print("base");
  dPP_2 = (TH1F*)fPP_in->Get(Form("hpp_HLT60_R%d_%s",radiusPP,etaWidth));
  dPP_2->Print("base");
  dPP_3 = (TH1F*)fPP_in->Get(Form("hpp_HLT40_R%d_%s",radiusPP,etaWidth));
  dPP_3->Print("base");
  dPP_Comb = (TH1F*)fPP_in->Get(Form("hpp_HLTComb_R%d_%s",radiusPP,etaWidth));   
  //dPP_Comb = (TH1F*)dPP_1->Clone(Form("hpp_TrgComb_R%d_n20_eta_p20",radiusPP,etaWidth));   
  dPP_Comb->Print("base");

  dPP_Comb->Scale(1./(5.3 * 1e9));
  
  for(int k = 1;k<=unfoldingCut;k++) {
    dPP_Comb->SetBinContent(k,0);
    dPP_1->SetBinContent(k,0);
    dPP_2->SetBinContent(k,0);
    dPP_3->SetBinContent(k,0);
  }
  
  // get PP MC
  mPP_Gen = (TH1F*)fPP_in->Get(Form("hpp_JetComb_gen_R%d_%s",radiusPP,etaWidth));
  mPP_Gen->Print("base");
  mPP_Reco = (TH1F*)fPP_in->Get(Form("hpp_JetComb_reco_R%d_%s",radiusPP,etaWidth));
  mPP_Reco->Print("base");
  mPP_Matrix = (TH2F*)fPP_in->Get(Form("hpp_matrix_HLT_R%d_%s",radiusPP,etaWidth));
  mPP_Matrix->Print("base");
  mPP_mcclosure_data = (TH1F*)fPP_in->Get(Form("hpp_mcclosure_JetComb_data_R%d_%s",radiusPP,etaWidth));
  mPP_mcclosure_data->Print("base");
  mPP_mcclosure_Matrix = (TH2F*)fPP_in->Get(Form("hpp_mcclosure_matrix_HLT_R%d_%s",radiusPP,etaWidth));
  mPP_mcclosure_Matrix->Print("base");

  //RooUnfoldResponse ruResponsePP(mPP_Matrix->ProjectionY(),mPP_Matrix->ProjectionX(), mPP_Matrix,"","");
  //regularization parameter definition: 
  //RooUnfoldSvd unfoldSvdPP(&ruResponsePP, dPP_Comb, nSVDIter);
  //uPP_SVD = (TH1F*)unfoldSvdPP.Hreco();

  
  // for(int k = 1;k<=unfoldingCut;k++){
  //   mPP_Gen->SetBinContent(k,0);
  //   mPP_Reco->SetBinContent(k,0);
  //   mPP_mcclosure_data->SetBinContent(k,0);
  //   for(int l = 1;l<=1000;l++){
  //     mPP_Matrix->SetBinContent(k,l,0);
  //     mPP_mcclosure_Matrix->SetBinContent(k,l,0);
  //     mPP_Matrix->SetBinContent(l,k,0);
  //     mPP_mcclosure_Matrix->SetBinContent(l,k,0);
  //   }
  // }

  
  if(printDebug) cout<<"Filling the PbPb response Matrix"<<endl;

  // response matrix and unfolding for PbPb 
  // going to try it the way kurt has it. 

  for(int i = 0;i<nbins_cent;i++){
    if(printDebug) cout<<"centrality bin iteration = "<<i<<endl;
    TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
    f->SetParameters(1e10,-8.8,40);
    // TH1F *hGenSpectraCorr = (TH1F*)mPbPb_Matrix[i]->ProjectionX()->Clone(Form("hGenSpectraCorr_cent%d",i));
    // hGenSpectraCorr->Fit("f"," ");
    // hGenSpectraCorr->Fit("f","","");
    // hGenSpectraCorr->Fit("f","LL");
    // TH1F *fHist = functionHist(f,hGenSpectraCorr,Form("fHist_cent%d",i));// function that you get from the fitting 
    // hGenSpectraCorr->Divide(fHist);
    for (int y=1;y<=mPbPb_Matrix[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_Matrix[i]->GetNbinsX();x++) {
	if (mPbPb_Matrix[i]->GetBinContent(x,y)<=1*mPbPb_Matrix[i]->GetBinError(x,y)) {
	  //in the above line mine had 0*getbinerror while Kurt's had 1*. 
	  mPbPb_Matrix[i]->SetBinContent(x,y,0);
	  mPbPb_Matrix[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_Matrix[i]->GetBinContent(x,y);
      }
      
      for (int x=1;x<=mPbPb_Matrix[i]->GetNbinsX();x++) {	   
	double ratio = 1;
	// if (hGenSpectraCorr->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorr->GetBinContent(x);
	mPbPb_Matrix[i]->SetBinContent(x,y,mPbPb_Matrix[i]->GetBinContent(x,y)*ratio);
	mPbPb_Matrix[i]->SetBinError(x,y,mPbPb_Matrix[i]->GetBinError(x,y)*ratio);
      }
    }
    //mPbPb_Matrix[i]->Smooth(0);
    // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  
    
    mPbPb_Response[i] = (TH2F*)mPbPb_Matrix[i]->Clone(Form("mPbPb_Response_cent%d",i));
    TH1F *hProj = (TH1F*)mPbPb_Response[i]->ProjectionY()->Clone(Form("hProj_cent%d",i));

    for (int y=1;y<=mPbPb_Response[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_Response[i]->GetNbinsX();x++) {
	if (mPbPb_Response[i]->GetBinContent(x,y)<=1*mPbPb_Response[i]->GetBinError(x,y)) {
	  // in the above if loop, kurt has 1*error and my old had 0*error
	  mPbPb_Response[i]->SetBinContent(x,y,0);
	  mPbPb_Response[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_Response[i]->GetBinContent(x,y);
      }
      
      for (int x=1;x<=mPbPb_Response[i]->GetNbinsX();x++) {  	
	if (sum==0) continue;
	double ratio = 1;
	//if(dPbPb_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
	// else ratio = dPbPb_TrgComb[i]->GetBinContent(y)/sum
	ratio = 1./sum;
	if (hProj->GetBinContent(y)==0) ratio = 1e-100/sum;
	else ratio = hProj->GetBinContent(y)/sum;
	mPbPb_Response[i]->SetBinContent(x,y,mPbPb_Response[i]->GetBinContent(x,y)*ratio);
	mPbPb_Response[i]->SetBinError(x,y,mPbPb_Response[i]->GetBinError(x,y)*ratio);
      }
    }
    
    mPbPb_ResponseNorm[i] = (TH2F*)mPbPb_Matrix[i]->Clone(Form("mPbPb_ResponseNorm_cent%d",i));
    for (int x=1;x<=mPbPb_ResponseNorm[i]->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=mPbPb_ResponseNorm[i]->GetNbinsY();y++) {
	if (mPbPb_ResponseNorm[i]->GetBinContent(x,y)<=1*mPbPb_ResponseNorm[i]->GetBinError(x,y)) {
	  mPbPb_ResponseNorm[i]->SetBinContent(x,y,0);
	  mPbPb_ResponseNorm[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_ResponseNorm[i]->GetBinContent(x,y);
      }
      
      for (int y=1;y<=mPbPb_ResponseNorm[i]->GetNbinsY();y++) {  	
	if (sum==0) continue;
	double ratio = 1./sum;
	mPbPb_ResponseNorm[i]->SetBinContent(x,y,mPbPb_ResponseNorm[i]->GetBinContent(x,y)*ratio);
	mPbPb_ResponseNorm[i]->SetBinError(x,y,mPbPb_ResponseNorm[i]->GetBinError(x,y)*ratio);
      }
      
    }
    
    
  }

  
  if(printDebug) cout<<"Filling PP response Matrix"<<endl;

  // response matrix for pp.  
  // Kurt doesnt have this whole hGenSpectraCorr thing in his macro. need to check why the difference exists between out codes
  
  TF1 *fpp = new TF1("fpp","[0]*pow(x+[2],[1])");
  fpp->SetParameters(1e10,-8.8,40);
  // if(printDebug) cout<<"before getting the gen spectra corr matrix"<<endl;
  // TH1F *hGenSpectraCorrPP = (TH1F*)mPP_Matrix->ProjectionX()->Clone("hGenSpectraCorrPP");
  // if(printDebug) cout<<"after gettign the gen spectra corr matrix"<<endl;
  // hGenSpectraCorrPP->Fit("f"," ");
  // hGenSpectraCorrPP->Fit("f","","");
  // hGenSpectraCorrPP->Fit("f","LL");
  // TH1F *fHistPP = functionHist(fpp,hGenSpectraCorrPP,"fHistPP");// that the function that you get from the fitting 
  // hGenSpectraCorrPP->Divide(fHistPP);
  
  for (int y=1;y<=mPP_Matrix->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_Matrix->GetNbinsX();x++) {
      if (mPP_Matrix->GetBinContent(x,y)<=1*mPP_Matrix->GetBinError(x,y)) {
	mPP_Matrix->SetBinContent(x,y,0);
	mPP_Matrix->SetBinError(x,y,0);
      }
      sum+=mPP_Matrix->GetBinContent(x,y);
    }
    
    for (int x=1;x<=mPP_Matrix->GetNbinsX();x++) {	   
      double ratio = 1;
      // if (hGenSpectraCorrPP->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorrPP->GetBinContent(x);
      mPP_Matrix->SetBinContent(x,y,mPP_Matrix->GetBinContent(x,y)*ratio);
      mPP_Matrix->SetBinError(x,y,mPP_Matrix->GetBinError(x,y)*ratio);
    }
  }
  // mPbPb_Matrix[i]->Smooth(0);
  
  // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  

  if(printDebug) cout<<"getting the response matrix"<<endl;

  mPP_Response = (TH2F*)mPP_Matrix->Clone("mPP_Response");
  TH1F *hProjPP = (TH1F*)mPP_Response->ProjectionY()->Clone("hProjPP");
  
  
  for (int y=1;y<=mPP_Response->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_Response->GetNbinsX();x++) {
      if (mPP_Response->GetBinContent(x,y)<=1*mPP_Response->GetBinError(x,y)) {
	// in the above if statement, kurt has 1*error and my old has 0*error
	mPP_Response->SetBinContent(x,y,0);
	mPP_Response->SetBinError(x,y,0);
      }
      sum+=mPP_Response->GetBinContent(x,y);
    }
    
    for (int x=1;x<=mPP_Response->GetNbinsX();x++) {  	
      if (sum==0) continue;
      double ratio = 1;
      //if(dPbPb_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
      // else ratio = dPbPb_TrgComb[i]->GetBinContent(y)/sum
      ratio = 1./sum;
      if (hProjPP->GetBinContent(y)==0) ratio = 1e-100/sum;
      else ratio = hProjPP->GetBinContent(y)/sum;
      mPP_Response->SetBinContent(x,y,mPP_Response->GetBinContent(x,y)*ratio);
      mPP_Response->SetBinError(x,y,mPP_Response->GetBinError(x,y)*ratio);
    }
  }
  if(printDebug) cout<<"getting the normalized response matrix"<<endl;
  mPP_ResponseNorm = (TH2F*)mPP_Matrix->Clone("mPP_ResponseNorm");
  for (int x=1;x<=mPP_ResponseNorm->GetNbinsX();x++) {
    double sum=0;
    for (int y=1;y<=mPP_ResponseNorm->GetNbinsY();y++) {
      if (mPP_ResponseNorm->GetBinContent(x,y)<=1*mPP_ResponseNorm->GetBinError(x,y)) {
	mPP_ResponseNorm->SetBinContent(x,y,0);
	mPP_ResponseNorm->SetBinError(x,y,0);
      }
      sum+=mPP_ResponseNorm->GetBinContent(x,y);
    }
    
    for (int y=1;y<=mPP_ResponseNorm->GetNbinsY();y++) {  	
      if (sum==0) continue;
      double ratio = 1./sum;
      mPP_ResponseNorm->SetBinContent(x,y,mPP_ResponseNorm->GetBinContent(x,y)*ratio);
      mPP_ResponseNorm->SetBinError(x,y,mPP_ResponseNorm->GetBinError(x,y)*ratio);
    }
    
    
  }
  
  // scale the spectra to the respective units

  // for(int i = 0;i<nbins_cent;++i){
  //   dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins_pt,Form("PbPb_measured_spectra_combined_cent%d",i),boundaries_pt);
  //   divideBinWidth(dPbPb_TrgComb[i]);
  // }

  // dPP_Comb = (TH1F*)dPP_Comb->Rebin(nbins_pt,"pp_measured_spectra_combined",boundaries_pt);
  // divideBinWidth(dPP_Comb);
  // dPP_Comb->Scale(1./ dPP_Comb->GetBinContent(nbins_pt));
  
  // Now that we have all the response matrix for the 6 centralities in PbPb and one pp spectra lets start doing the steps:
  // we have 39 pt bins, so we need 1000 gaussian functions for each pt bin.
  
  Int_t unfoldingTrials = 100;
  Double_t meanMeasPbPb[nbins_pt][nbins_cent], sigmaMeasPbPb[nbins_pt][nbins_cent];
  Double_t meanMeasPP[nbins_pt], sigmaMeasPP[nbins_pt];
  Double_t meanUnfoldPbPb[nbins_pt][nbins_cent][unfoldingTrials], sigmaUnfoldPbPb[nbins_pt][nbins_cent][unfoldingTrials];
  Double_t meanUnfoldPP[nbins_pt][unfoldingTrials], sigmaUnfoldPP[nbins_pt][unfoldingTrials]; 
  
  TRandom3 *random = new TRandom3(0);

  for(int u = 0;u<unfoldingTrials;++u){
    cout<<"unfolding trial no = "<<u+1<<endl;
  
    for(int j = 0;j<nbins_pt;++j){
      for(int i = 0;i<nbins_cent;++i){
      
	meanMeasPbPb[j][i] = dPbPb_TrgComb[i]->GetBinContent(j+1);
	sigmaMeasPbPb[j][i] = dPbPb_TrgComb[i]->GetBinError(j+1);

      }// centrality loop

      meanMeasPP[j] = dPP_Comb->GetBinContent(j+1);
      sigmaMeasPP[j] = dPP_Comb->GetBinContent(j+1);
      
    }// nbins_pt loop

    // now proceed to unfolding for each trial.

    for(int i = 0;i<nbins_cent;++i){
      //cout<<"centrality = "<<i<<endl;

      TH1F * hPreUnfoldingSpectra = new TH1F("hPreUnfoldingSpectra","",nbins_pt,0,nbins_pt);
      TH1F * hAfterUnfoldingSpectra;

      for(int j = 0;j<nbins_pt;++j){
	
	hPreUnfoldingSpectra->SetBinContent(j+1, random->Gaus(meanMeasPbPb[j][i], sigmaMeasPbPb[j][i]));
	hPreUnfoldingSpectra->SetBinError(j+1, sigmaMeasPbPb[j][i]/sqrt(unfoldingTrials));
        //if(j==100)cout << " before unfolding bin " << j << " value = " << hPreUnfoldingSpectra->GetBinContent(j+1)<<endl;
        //if(j==100)cout << " before unfolding bin " << j << " error = " << hPreUnfoldingSpectra->GetBinError(j+1)<<endl;
	
      }// nbins_pt loop

      TH1F* hMCGen          = (TH1F*)mPbPb_Response[i]->ProjectionX();
      removeZero(hMCGen);
      //cout << " MC bin " << 100 << " value = " << hMCGen->GetBinContent(100)<<endl;
      bayesianUnfold myUnfoldingMulti(mPbPb_Matrix[i], hMCGen, 0);
      myUnfoldingMulti.unfold(hPreUnfoldingSpectra, BayesIter);

      hAfterUnfoldingSpectra = (TH1F*) myUnfoldingMulti.hPrior->Clone("hAfterUnfoldingSpectra");

      for(int j = 0;j<nbins_pt;++j){
	
	//if(j==100)cout << " before unfolding bin " << j << " value = " << hPreUnfoldingSpectra->GetBinContent(j+1)<<endl;
	//if(j==100)cout << " after  unfolding bin " << j << " value = " << hAfterUnfoldingSpectra->GetBinContent(j+1)<<endl;
	
	meanUnfoldPbPb[j][i][u] = hAfterUnfoldingSpectra->GetBinContent(j+1);
	sigmaUnfoldPbPb[j][i][u] = hAfterUnfoldingSpectra->GetBinError(j+1);

	// cout << "after unfolding meanUnfoldPbPb[" << j << "][" << i << "][" << u<< "] = " <<meanUnfoldPbPb[j][i][u]<<"    ";
	// cout << "after unfolding meanUnfoldPbPb[" << j << "][" << i << "][" << u<< "] = " <<sigmaUnfoldPbPb[j][i][u]<<endl;
	
      }// nbins_pt loop
      
      //hPreUnfoldingSpectra->Print("base");
      //hAfterUnfoldingSpectra->Print("base");
      
      delete hPreUnfoldingSpectra;
      delete hAfterUnfoldingSpectra;
      delete hMCGen; 
      
    }// centrality loop

    cout<<"pp "<<endl;

    // now do it for the pp:
    TH1F * hPreUnfoldingSpectraPP = new TH1F("hPreUnfoldingSpectraPP","",nbins_pt,0,nbins_pt);
    TH1F * hAfterUnfoldingSpectraPP;
    
    for(int j = 0;j<nbins_pt;++j){
	
      hPreUnfoldingSpectraPP->SetBinContent(j+1, random->Gaus(meanMeasPP[j], sigmaMeasPP[j]));
      hPreUnfoldingSpectraPP->SetBinError(j+1, sigmaMeasPP[j]/sqrt(unfoldingTrials));
        
    }// nbins_pt loop
    TH1F* hMCGenPP          = (TH1F*)mPP_Response->ProjectionX();
    removeZero(hMCGenPP);
    bayesianUnfold myUnfoldingMultiPP(mPP_Matrix, hMCGenPP, 0);
    myUnfoldingMultiPP.unfold(hPreUnfoldingSpectraPP, BayesIter);

    hAfterUnfoldingSpectraPP = (TH1F*) myUnfoldingMultiPP.hPrior->Clone("hAfterUnfoldingSpectraPP");

    for(int j = 0;j<nbins_pt;++j){

      meanUnfoldPP[j][u] = hAfterUnfoldingSpectraPP->GetBinContent(j+1);
      sigmaUnfoldPP[j][u] = hAfterUnfoldingSpectraPP->GetBinError(j+1);

    }// nbins_pt loop

    delete hPreUnfoldingSpectraPP;
    delete hAfterUnfoldingSpectraPP;
    delete hMCGenPP; 
    
  }// unfolding trials loop


  // Now that we have all the necesary values we need, lets proceed to fill a histogram with the mean values for each ptbin and get the corrected values.
  TH1F * hAfterUnfoldingptBinDistribution[nbins_pt];
  TH1F * hCorrUnfoldingPbPb[nbins_cent];
  
  for(int i = 0;i<nbins_cent;++i){

    hCorrUnfoldingPbPb[i] = new TH1F(Form("PbPb_BayesianUnfolded_cent%d",i),"Spectra after correction", nbins_pt, 0, nbins_pt);

    for(int j = 0;j<nbins_pt;++j){
      
      //hAfterUnfoldingptBinDistribution[j] = new TH1F(Form("hAfterUnfoldingptBinDistribution_ptBin%d",j),"",100,	(meanMeasPbPb[j][i]-10) * sigmaMeasPbPb[j][i], (meanMeasPbPb[j][i]+10) * sigmaMeasPbPb[j][i]);
      hAfterUnfoldingptBinDistribution[j] = new TH1F(Form("hAfterUnfoldingptBinDistribution_ptBin%d",j),"",100,	0, 1);
      for(int u = 0;u<unfoldingTrials;++u){

	hAfterUnfoldingptBinDistribution[j]->Fill(meanUnfoldPbPb[j][i][u]);

	//if(j==100) cout<< "unfolding_trial = " << u+1 << " mean unfold value = "<< meanUnfoldPbPb[j][i][u] <<endl;

      }// unfolding trials loop

      //if(j==100) cout<<"Mean of that value for pt=100 = "<< (Float_t)hAfterUnfoldingptBinDistribution[j]->GetMean() <<endl;      
      hCorrUnfoldingPbPb[i]->SetBinContent(j+1, hAfterUnfoldingptBinDistribution[j]->GetMean());
      hCorrUnfoldingPbPb[i]->SetBinError(j+1, hAfterUnfoldingptBinDistribution[j]->GetRMS());

      delete hAfterUnfoldingptBinDistribution[j];
      
    }// nbins_pt loop

  }// centrality loop

  // similar for the pp:
  TH1F * hAfterUnfoldingptBinDistributionPP[nbins_pt];
  TH1F * hCorrUnfoldingPP;
  
  hCorrUnfoldingPP = new TH1F("PP_BayesianUnfolded","Spectra after unfolding error correction",nbins_pt, 0, nbins_pt);
  
  for(int j = 0;j<nbins_pt;++j){
    
    //hAfterUnfoldingptBinDistributionPP[j] = new TH1F(Form("hAfterUnfoldingptBinDistributionPP_ptBin%d",j),"",1000,(meanMeasPP[j]-10) * sigmaMeasPP[j], (meanMeasPP[j]+10) * sigmaMeasPP[j]);
    hAfterUnfoldingptBinDistributionPP[j] = new TH1F(Form("hAfterUnfoldingptBinDistributionPP_ptBin%d",j),"",100, 0, 1);
    for(int u = 0;u<unfoldingTrials;++u){
      
      hAfterUnfoldingptBinDistributionPP[j]->Fill(meanUnfoldPP[j][u]);
      
    }// unfolding trials loop
    
    hCorrUnfoldingPP->SetBinContent(j+1, hAfterUnfoldingptBinDistributionPP[j]->GetMean());
    hCorrUnfoldingPP->SetBinError(j+1, hAfterUnfoldingptBinDistributionPP[j]->GetRMS());
    
    delete hAfterUnfoldingptBinDistributionPP[j];
    
  }// nbins_pt loop
    
  TFile f(Form("../../Output/PbPb_R%d_pp_withresiduals_R%d_%s_unfoldingCut_%d_data_driven_correction_ak%s%s_%d.root",radius, radiusPP, etaWidth ,unfoldingCut,algo,jet_type,date.GetDate()),"RECREATE");
  f.cd();

  for(int i = 0;i<nbins_cent;i++) {
    hCorrUnfoldingPbPb[i]->Write();
    hCorrUnfoldingPbPb[i]->Print("base");

    dPbPb_TrgComb[i]->Write();
    dPbPb_TrgComb[i]->Print("base");
    
  }

  hCorrUnfoldingPP->Write();
  hCorrUnfoldingPP->Print("base");
  dPP_Comb->Write();
  dPP_Comb->Print("base");
  
  f.Write();
  f.Close();

  timer.Stop();
  if(printDebug) cout<<"CPU time (mins) = "<<(Float_t)timer.CpuTime()/60<<endl;
  if(printDebug) cout<<"Real tile (mins) = "<<(Float_t)timer.RealTime()/60<<endl;
  

}
