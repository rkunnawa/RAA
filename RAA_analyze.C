// Raghav Kunnawalkam Elayavalli
// June 10th 2014
// CERN
// raghav.k.e AT CERN DOT CH

// RAA - analyze macro. following the steps of having a separate macro to do the work. 

// this will be taken mostly from Unfold_RAA_V0.C from MIT. Im just keeping it here so that it will be easier to run. 
// important point in the code is that in the loops for pbpb, the array index corrsponding to nbins_cent represents the 0-100 centrality bin. 
// PbPb and pp are separated here - as they should be. 

// June 25th - add separate files for reading in data from pp and PbPb. 
//           - Run the full macro to check if it works. 
//           - Wooo Hoooo! Macro runs fully without any error and produces the histograms 
//           - Now to go ahead and check if they contain meaningful stuff. 
//           - obviously till i get the full dataset, i cant get the lumi normalization required to add the triggers. 

// June 30th - changed the names of hitograms to match what they mean. especially the unfolded ones. 
//           - Ok there is a slight issue in the naming of histograms for the bayesian iteration unfolding. 
//             - the histograms are called [nbins_cent][iterations] while they are saved ad [iterations][nbins_cent]
//             - Just be careful when loading the histograms, especially in the plotting macro 

// July 19th - started adding in the radius loop and the eta bins loop.

// Aug 20th - continued editing the radius and eta bin loop. 

// Nov 4th - have the unfolding closure test working. 

// Nov 5th - working on the data bayesin closure using the fine bins technique. 
 

#include <iostream>
#include <stdio.h>

#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/prior.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/bayesianUnfold.h"
#include "TStopwatch.h"

//static const int nbins_pt = 29;
//static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

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


void RAA_analyze(int radius = 3, char* algo = "Pu", char *jet_type = "PF"){

  TStopwatch timer; 
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  bool printDebug = true;

  // get the data and mc histograms from the output of the read macro. 
  
  TDatime date;//this is just here to get them to run optimized. 

  TFile* fData_PbPb_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_data_ak%s%s_testComb2_cut3_test_20141110.root",algo,jet_type));
  TFile *fData_pp_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/pp_data_ak%s_20140829.root",jet_type));
  TFile* fMC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_mc_closure_oppside_ak%s%s_20141111.root",algo,jet_type));

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
  
  TH1F *dPP_1, *dPP_2, *dPP_3, *dPP_Comb;
  
  TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix, *mPP_Response,*mPP_ResponseNorm;
  TH1F *mPP_mcclosure_data;
  TH2F *mPP_mcclosure_Matrix, *mPP_mcclosure_Response,*mPP_mcclosure_ResponseNorm;
  TH1F *mPP_mcclosure_Gen;

  const int Iterations = 20; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1], *uPbPb_BinByBin[nbins_cent+1]; 
  TH1F *uPbPb_BayesianIter[nbins_cent+1][Iterations];

  TH1F *uPP_Bayes, *uPP_BinByBin;
  TH1F *uPP_BayesianIter[Iterations];
  
  
  // would be better to read in the histograms and rebin them. come to think of it, it would be better to have them already rebinned (and properly scaled - to the level of differential cross section in what ever barns (inverse micro barns) but keep it consistent) from the read macro. 
  
  // get PbPb data
  for(int i = 0;i<=nbins_cent;i++){
    if(printDebug) cout<<"cent_"<<i<<endl;
    dPbPb_TrgComb[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb_TrgObjComb_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_TrgComb[i]->Scale(4*145.156*1e6);
    dPbPb_TrgComb[i]->Print("base");
    dPbPb_Trg80[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb_TrgObj80_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_Trg80[i]->Scale(4*145.156*1e6);
    dPbPb_Trg80[i]->Print("base");
    dPbPb_Trg65[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb_TrgObj65_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_Trg65[i]->Scale(4*145.156*1e6);
    dPbPb_Trg65[i]->Print("base");
    dPbPb_Trg55[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb_TrgObj55_R%d_n20_eta_p20_cent%d",radius,i));
    //dPbPb_Trg55[i]->Scale(4*145.156*1e6);
    dPbPb_Trg55[i]->Print("base");

    /*
    dPbPb_Comb[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpbComb_cent%d",i));
    dPbPb_Comb[i]->Print("base");
    dPbPb_1[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb1_cent%d",i));
    dPbPb_1[i]->Print("base");
    dPbPb_2[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb2_cent%d",i));
    dPbPb_2[i]->Print("base");
    dPbPb_3[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb3_cent%d",i));
    dPbPb_3[i]->Print("base");

    dPbPb_80[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb_80_cent%d",i));
    dPbPb_80[i]->Print("base");
    dPbPb_65[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb_65_cent%d",i));
    dPbPb_65[i]->Print("base");
    dPbPb_55[i] = (TH1F*)fData_PbPb_in->Get(Form("hpbpb_55_cent%d",i));
    dPbPb_55[i]->Print("base");
    */

  }

  // get PbPb MC
  for(int i = 0;i<=nbins_cent;i++){
    
    mPbPb_Gen[i] = (TH1F*)fMC_in->Get(Form("hpbpb_gen_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Gen[i]->Print("base");
    mPbPb_Reco[i] = (TH1F*)fMC_in->Get(Form("hpbpb_reco_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Reco[i]->Print("base");
    mPbPb_Matrix[i] = (TH2F*)fMC_in->Get(Form("hpbpb_matrix_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Matrix[i]->Print("base");
    mPbPb_mcclosure_data[i] = (TH1F*)fMC_in->Get(Form("hpbpb_mcclosure_data_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_mcclosure_data[i]->Print("base");
    mPbPb_mcclosure_gen[i] = (TH1F*)fMC_in->Get(Form("hpbpb_mcclosure_gen_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_mcclosure_gen[i]->Print("base");
    mPbPb_mcclosure_Matrix[i] = (TH2F*)fMC_in->Get(Form("hpbpb_mcclosure_matrix_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_mcclosure_Matrix[i]->Print("base");
    
    //mPbPb_Response[i] = new TH2F(Form("mPbPb_Response_cent%d",i),"Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    //mPbPb_ResponseNorm[i] = new TH2F(Form("mPbPb_ResponseNorm_cent%d",i),"Normalized Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  }
  
  // get PP data
  if(printDebug) cout<<"Getting PP data and MC"<<endl;
  dPP_1 = (TH1F*)fData_pp_in->Get(Form("hpp_Trg80_R%d_n20_eta_p20",radius));
  dPP_1->Print("base");
  dPP_2 = (TH1F*)fData_pp_in->Get(Form("hpp_Trg60_R%d_n20_eta_p20",radius));
  dPP_2->Print("base");
  dPP_3 = (TH1F*)fData_pp_in->Get(Form("hpp_Trg40_R%d_n20_eta_p20",radius));
  dPP_3->Print("base");
  dPP_Comb = (TH1F*)fData_pp_in->Get(Form("hpp_TrgComb_R%d_n20_eta_p20",radius));
  dPP_Comb->Print("base");
  
  // get PP MC
  mPP_Gen = (TH1F*)fMC_in->Get(Form("hpp_gen_R%d_n20_eta_p20",radius));
  mPP_Gen->Print("base");
  mPP_Reco = (TH1F*)fMC_in->Get(Form("hpp_reco_R%d_n20_eta_p20",radius));
  mPP_Reco->Print("base");
  mPP_Matrix = (TH2F*)fMC_in->Get(Form("hpp_matrix_R%d_n20_eta_p20",radius));
  mPP_Matrix->Print("base");
  mPP_mcclosure_data = (TH1F*)fMC_in->Get(Form("hpp_mcclosure_data_R%d_n20_eta_p20",radius));
  mPP_mcclosure_data->Print("base");
  mPP_mcclosure_Matrix = (TH2F*)fMC_in->Get(Form("hpp_mcclosure_matrix_R%d_n20_eta_p20",radius));
  mPP_mcclosure_Matrix->Print("base");
  //mPP_Matrix->Print("base");
  //mPP_Response = (TH2F*)fMc_in->Get("hpp_gen");
  
  // make the response matrix.
  // here since we dont have the simple nature of the uhist histograms which makes debugging hard, we have to run the response matrix and unfolding separately for PbPb and pp which makes code ugly and not efficient but easy to debug at the same time. 
  
  if(printDebug) cout<<"Filling the PbPb response Matrix"<<endl;

  // response matrix and unfolding for PbPb 
  // going to try it the way kurt has it. 

  for(int i = 0;i<=nbins_cent;i++){
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

  

  if(printDebug) cout<<"Filling the PbPb mcclosure response Matrix"<<endl;

  // response matrix and unfolding for PbPb for the closure test 
  // going to try it the way kurt has it. 

  for(int i = 0;i<=nbins_cent;i++){
    if(printDebug) cout<<"centrality bin iteration = "<<i<<endl;
    TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
    f->SetParameters(1e10,-8.8,40);
    // TH1F *hGenSpectraCorr = (TH1F*)mPbPb_mcclosure_Matrix[i]->ProjectionX()->Clone(Form("hGenSpectraCorr_cent%d",i));
    // hGenSpectraCorr->Fit("f"," ");
    // hGenSpectraCorr->Fit("f","","");
    // hGenSpectraCorr->Fit("f","LL");
    // TH1F *fHist = functionHist(f,hGenSpectraCorr,Form("fHist_cent%d",i));// function that you get from the fitting 
    // hGenSpectraCorr->Divide(fHist);
    for (int y=1;y<=mPbPb_mcclosure_Matrix[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_mcclosure_Matrix[i]->GetNbinsX();x++) {
	if (mPbPb_mcclosure_Matrix[i]->GetBinContent(x,y)<=1*mPbPb_mcclosure_Matrix[i]->GetBinError(x,y)) {
	  //in the above line mine had 0*getbinerror while Kurt's had 1*. 
	  mPbPb_mcclosure_Matrix[i]->SetBinContent(x,y,0);
	  mPbPb_mcclosure_Matrix[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_mcclosure_Matrix[i]->GetBinContent(x,y);
      }
      
      for (int x=1;x<=mPbPb_mcclosure_Matrix[i]->GetNbinsX();x++) {	   
	double ratio = 1;
	// if (hGenSpectraCorr->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorr->GetBinContent(x);
	mPbPb_mcclosure_Matrix[i]->SetBinContent(x,y,mPbPb_mcclosure_Matrix[i]->GetBinContent(x,y)*ratio);
	mPbPb_mcclosure_Matrix[i]->SetBinError(x,y,mPbPb_mcclosure_Matrix[i]->GetBinError(x,y)*ratio);
      }
    }
    //mPbPb_mcclosure_Matrix[i]->Smooth(0);
    // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  
    
    mPbPb_mcclosure_Response[i] = (TH2F*)mPbPb_mcclosure_Matrix[i]->Clone(Form("mPbPb_mcclosure_Response_cent%d",i));
    TH1F *hProj = (TH1F*)mPbPb_mcclosure_Response[i]->ProjectionY()->Clone(Form("hProj_cent%d",i));

    for (int y=1;y<=mPbPb_mcclosure_Response[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_mcclosure_Response[i]->GetNbinsX();x++) {
	if (mPbPb_mcclosure_Response[i]->GetBinContent(x,y)<=1*mPbPb_mcclosure_Response[i]->GetBinError(x,y)) {
	  // in the above if loop, kurt has 1*error and my old had 0*error
	  mPbPb_mcclosure_Response[i]->SetBinContent(x,y,0);
	  mPbPb_mcclosure_Response[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_mcclosure_Response[i]->GetBinContent(x,y);
      }
      
      for (int x=1;x<=mPbPb_mcclosure_Response[i]->GetNbinsX();x++) {  	
	if (sum==0) continue;
	double ratio = 1;
	//if(dPbPb_mcclosure_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
	// else ratio = dPbPb_mcclosure_TrgComb[i]->GetBinContent(y)/sum
	ratio = 1./sum;
	if (hProj->GetBinContent(y)==0) ratio = 1e-100/sum;
	else ratio = hProj->GetBinContent(y)/sum;
	mPbPb_mcclosure_Response[i]->SetBinContent(x,y,mPbPb_mcclosure_Response[i]->GetBinContent(x,y)*ratio);
	mPbPb_mcclosure_Response[i]->SetBinError(x,y,mPbPb_mcclosure_Response[i]->GetBinError(x,y)*ratio);
      }
    }
    
    mPbPb_mcclosure_ResponseNorm[i] = (TH2F*)mPbPb_mcclosure_Matrix[i]->Clone(Form("mPbPb_mcclosure_ResponseNorm_cent%d",i));
    for (int x=1;x<=mPbPb_mcclosure_ResponseNorm[i]->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=mPbPb_mcclosure_ResponseNorm[i]->GetNbinsY();y++) {
	if (mPbPb_mcclosure_ResponseNorm[i]->GetBinContent(x,y)<=1*mPbPb_mcclosure_ResponseNorm[i]->GetBinError(x,y)) {
	  mPbPb_mcclosure_ResponseNorm[i]->SetBinContent(x,y,0);
	  mPbPb_mcclosure_ResponseNorm[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_mcclosure_ResponseNorm[i]->GetBinContent(x,y);
      }
      
      for (int y=1;y<=mPbPb_mcclosure_ResponseNorm[i]->GetNbinsY();y++) {  	
	if (sum==0) continue;
	double ratio = 1./sum;
	mPbPb_mcclosure_ResponseNorm[i]->SetBinContent(x,y,mPbPb_mcclosure_ResponseNorm[i]->GetBinContent(x,y)*ratio);
	mPbPb_mcclosure_ResponseNorm[i]->SetBinError(x,y,mPbPb_mcclosure_ResponseNorm[i]->GetBinError(x,y)*ratio);
      }
      
    }
  }


  
  if(printDebug) cout<<"Filling PP response Matrix"<<endl;

  // response matrix for pp.  
  // Kurt doesnt have this whole hGenSpectraCorr thing in his macro. need to check why the difference exits between out codes
  
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
  


  if(printDebug) cout<<"Filling PP mcclosure response Matrix"<<endl;

  // response matrix for pp.  
  // Kurt doesnt have this whole hGenSpectraCorr thing in his macro. need to check why the difference exits between out codes
  
  TF1 *fpp_mcclosure = new TF1("fpp_mcclosure","[0]*pow(x+[2],[1])");
  fpp_mcclosure->SetParameters(1e10,-8.8,40);
  // if(printDebug) cout<<"before getting the gen spectra corr matrix"<<endl;
  // TH1F *hGenSpectraCorrPP = (TH1F*)mPP_mcclosure_Matrix->ProjectionX()->Clone("hGenSpectraCorrPP");
  // if(printDebug) cout<<"after gettign the gen spectra corr matrix"<<endl;
  // hGenSpectraCorrPP->Fit("f"," ");
  // hGenSpectraCorrPP->Fit("f","","");
  // hGenSpectraCorrPP->Fit("f","LL");
  // TH1F *fHistPP = functionHist(fpp,hGenSpectraCorrPP,"fHistPP");// that the function that you get from the fitting 
  // hGenSpectraCorrPP->Divide(fHistPP);
  
  for (int y=1;y<=mPP_mcclosure_Matrix->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_mcclosure_Matrix->GetNbinsX();x++) {
      if (mPP_mcclosure_Matrix->GetBinContent(x,y)<=1*mPP_mcclosure_Matrix->GetBinError(x,y)) {
	mPP_mcclosure_Matrix->SetBinContent(x,y,0);
	mPP_mcclosure_Matrix->SetBinError(x,y,0);
      }
      sum+=mPP_mcclosure_Matrix->GetBinContent(x,y);
    }
    
    for (int x=1;x<=mPP_mcclosure_Matrix->GetNbinsX();x++) {	   
      double ratio = 1;
      // if (hGenSpectraCorrPP->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorrPP->GetBinContent(x);
      mPP_mcclosure_Matrix->SetBinContent(x,y,mPP_mcclosure_Matrix->GetBinContent(x,y)*ratio);
      mPP_mcclosure_Matrix->SetBinError(x,y,mPP_mcclosure_Matrix->GetBinError(x,y)*ratio);
    }
  }
  // mPbPb_Matrix[i]->Smooth(0);
  
  // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  

  if(printDebug) cout<<"getting the response matrix"<<endl;

  mPP_mcclosure_Response = (TH2F*)mPP_mcclosure_Matrix->Clone("mPP_mcclosure_Response");
  TH1F *hProjPP_mcclosure = (TH1F*)mPP_mcclosure_Response->ProjectionY()->Clone("hProjPP_mcclosure");
  
  
  for (int y=1;y<=mPP_mcclosure_Response->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_mcclosure_Response->GetNbinsX();x++) {
      if (mPP_mcclosure_Response->GetBinContent(x,y)<=1*mPP_mcclosure_Response->GetBinError(x,y)) {
	// in the above if statement, kurt has 1*error and my old has 0*error
	mPP_mcclosure_Response->SetBinContent(x,y,0);
	mPP_mcclosure_Response->SetBinError(x,y,0);
      }
      sum+=mPP_mcclosure_Response->GetBinContent(x,y);
    }
    
    for (int x=1;x<=mPP_mcclosure_Response->GetNbinsX();x++) {  	
      if (sum==0) continue;
      double ratio = 1;
      //if(dPbPb_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
      // else ratio = dPbPb_TrgComb[i]->GetBinContent(y)/sum
      ratio = 1./sum;
      if (hProjPP_mcclosure->GetBinContent(y)==0) ratio = 1e-100/sum;
      else ratio = hProjPP_mcclosure->GetBinContent(y)/sum;
      mPP_mcclosure_Response->SetBinContent(x,y,mPP_mcclosure_Response->GetBinContent(x,y)*ratio);
      mPP_mcclosure_Response->SetBinError(x,y,mPP_mcclosure_Response->GetBinError(x,y)*ratio);
    }
  }
  if(printDebug) cout<<"getting the normalized response matrix"<<endl;
  mPP_mcclosure_ResponseNorm = (TH2F*)mPP_mcclosure_Matrix->Clone("mPP_mcclosure_ResponseNorm");
  for (int x=1;x<=mPP_mcclosure_ResponseNorm->GetNbinsX();x++) {
    double sum=0;
    for (int y=1;y<=mPP_mcclosure_ResponseNorm->GetNbinsY();y++) {
      if (mPP_mcclosure_ResponseNorm->GetBinContent(x,y)<=1*mPP_mcclosure_ResponseNorm->GetBinError(x,y)) {
	mPP_mcclosure_ResponseNorm->SetBinContent(x,y,0);
	mPP_mcclosure_ResponseNorm->SetBinError(x,y,0);
      }
      sum+=mPP_mcclosure_ResponseNorm->GetBinContent(x,y);
    }
    
    for (int y=1;y<=mPP_mcclosure_ResponseNorm->GetNbinsY();y++) {  	
      if (sum==0) continue;
      double ratio = 1./sum;
      mPP_mcclosure_ResponseNorm->SetBinContent(x,y,mPP_mcclosure_ResponseNorm->GetBinContent(x,y)*ratio);
      mPP_mcclosure_ResponseNorm->SetBinError(x,y,mPP_mcclosure_ResponseNorm->GetBinError(x,y)*ratio);
    }
    
    
  }

  

  if(printDebug) cout<<"finished with all the response matrix. now going for unfolding"<<endl;
  
  // do the unfolding - including the iteration systematics (ofcourse). Similar to the above version, we have to create 2 separate unfolding for PbPb and pp. 

  
 
  // first for PbPb
  for (int i=0;i<=nbins_cent;i++) {

    // Do Bin-by-bin
    if(printDebug) cout<<"doing bin by bin unfolding for PbPb for centrality = "<<i<<endl;

    TH1F* hBinByBinCorRaw = (TH1F*)mPbPb_Response[i]->ProjectionY();
    TH1F* hMCGen          = (TH1F*)mPbPb_Response[i]->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",50,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    //TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));

    if(printDebug) cout<<"passed here after binbybincorrraw"<<endl;
    dPbPb_TrgComb[i]->Print("base");
    uPbPb_BinByBin[i] = (TH1F*)dPbPb_TrgComb[i]->Clone(Form("uPbPb_BinByBin_cent%d",i));
    if(printDebug) cout<<"maybe here?"<<endl;
    uPbPb_BinByBin[i]->Divide(hBinByBinCor);
    //      uPbPb_BinByBin[i] = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
    if(printDebug) cout<<"passed A"<<endl;
    // Do unfolding
    //prior myPrior(mPbPb_Matrix[i],dPbPb_TrgComb[i],0);
    //myPrior.unfold(dPbPb_TrgComb[i],1);
    TH1F* hPrior = (TH1F*)hMCGen->Clone("hPrior");
    removeZero(hPrior);
    //hPrior->Scale(dPbPb_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
    if(printDebug) cout<<"passed B"<<endl;
    //TH1F *hReweightedPbPb = (TH1F*)mPbPb_Response[i]->ProjectionY()->Clone(Form("hReweightedPbPb_cent%d",i));

    //bayesianUnfold myUnfoldingJECSys(mPbPb_Matrix[i],hPrior,0);
    //myUnfoldingJECSys.unfold(dPbPb_TrgComb[i]JECSys,nBayesianIter);
    //bayesianUnfold myUnfoldingSmearSys(mPbPb_Matrix[i],hPrior,0);
    //myUnfoldingSmearSys.unfold(dPbPb_TrgComb[i]SmearSys,nBayesianIter);

    bayesianUnfold myUnfolding(mPbPb_Matrix[i],hPrior,0);
    myUnfolding.unfold(dPbPb_TrgComb[i],BayesIter);

    if(printDebug) cout <<"Unfolding bin "<<i<<endl;

    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<Iterations;j++){
      bayesianUnfold myUnfoldingSys(mPbPb_Matrix[i],hPrior,0);
      myUnfoldingSys.unfold(dPbPb_TrgComb[i],j);
      uPbPb_BayesianIter[i][j]= (TH1F*) myUnfoldingSys.hPrior->Clone(Form("uPbPb_BayesianIter%d_cent%d",j,i));
      uPbPb_BayesianIter[i][j]->Print("base");
      if(i!=6) uPbPb_BayesianIter[i][j]->SetTitle(Form("Unfolded PbPb Bayesian iteration %d with %2.0f - %2.0f cent",j,5*boundaries_cent[i],5*boundaries_cent[i+1]));
      else uPbPb_BayesianIter[i][j]->SetTitle(Form("Unfolded PbPb Bayesian iteration %d with 0-200 cent",j));
    }
    if(printDebug) cout<<"passed iteration sys"<<endl;
    uPbPb_Bayes[i]        = (TH1F*) uPbPb_BayesianIter[i][BayesIter]->Clone(Form("uPbPb_Bayes_cent%i",i));
    //uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJECSys_cent%i",i));
    //uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    //uPbPb_BinByBin[i] = (TH1F*) unfold2.Hreco();
    uPbPb_BinByBin[i]->SetName(Form("uPbPb_BinByBin_cent%i",i));
    if(i!=6){
      uPbPb_BinByBin[i]->SetTitle(Form("Unfolded PbPb bin by bin %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]));
      uPbPb_Bayes[i]->SetTitle(Form("Unfolded PbPb Bayes %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]));
    }else {
      uPbPb_BinByBin[i]->SetTitle("Unfolded PbPb bin by bin 0-200 cent");
      uPbPb_Bayes[i]->SetTitle("Unfolded PbPb Bayes 0-200 cent");
    }

    delete hPrior;
    

    // // Do unfolding
    // prior myPrior(mPbPb_Matrix[i],dPbPb_TrgComb[i],0.0);
    // myPrior.unfold(dPbPb_TrgComb[i],1);
    // TH1F *hPrior = (TH1F*)mPbPb_Matrix[i]->ProjectionX()->Clone(Form("hPrior_cent%d",i));
    // hPrior->Scale(dPbPb_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
    // bayesianUnfold myUnfoldingJECSys(mPbPb_Matrix[i],hPrior,0);
    // myUnfoldingJECSys.unfold(dPbPb_TrgComb[i]JECSys,nBayesianIter);
    // bayesianUnfold myUnfoldingSmearSys(mPbPb_Matrix[i],hPrior,0);
    // myUnfoldingSmearSys.unfold(dPbPb_TrgComb[i]SmearSys,nBayesianIter);
    // if(printDebug) cout <<"Unfolding bin "<<i<<endl;
    // // Iteration Systematics
    // for (int j=2;j<7;j++)
    // {
    // bayesianUnfold myUnfoldingSys(mPbPb_Matrix[i],hPrior,0);
    // myUnfoldingSys.unfold(dPbPb_TrgComb[i],j);
    // uhist[i]->hRecoIterSys[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
    // }
		
    // // Do Bin-by-bin
    // TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY(); 
    // TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
    // hBinByBinCorRaw->Divide(hMCGen);
    // TF1 *f = new TF1("f","[0]+[1]*x");
    // hBinByBinCorRaw->Fit("f","LL ","",100,300);
    // TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    // //      TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    // delete hBinByBinCorRaw,hMCGen;
    // uPbPb_BinByBin[i] = (TH1F*) dPbPb_TrgComb[i]->Clone(Form("hRecoBinByBin_cent%d",i));
    // uPbPb_BinByBin[i]->Divide(hBinByBinCor);
    // //      uPbPb_BinByBin[i] = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
		
    // uhist[i]->hReco         = (TH1F*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    // uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJeCSys_cent%i",i));
    // uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    // //uPbPb_BinByBin[i] = (TH1F*) unfold2.Hreco();
    // uPbPb_BinByBin[i]->SetName(Form("UnfoldedBinByBin_cent%i",i));
    
    

    // if (doToy) {
    //   TCanvas *cToy = new TCanvas("cToy","toy",600,600);
    //   int nExp=1000;
    //   TH1F *hTmp[nbins_truth+1];
    //   TH1F *hTmp2[nbins_truth+1];
    //   for (int j=1;j<=nbins_truth;j++) {
    // 	hTmp[j] = new TH1F(Form("hTmp%d",j),"",200,0,10.+uhist[i]->hReco->GetBinContent(j)*2);
    // 	hTmp2[j] = new TH1F(Form("hTmp2%d",j),"",200,0,10.+uPbPb_BinByBin[i]->GetBinContent(j)*2);
    //   }
    //   for (int exp =0; exp<nExp; exp++) {
    // 	TH1F *hToy = (TH1F*)dPbPb_TrgComb[i]->Clone();   
    // 	TH2F *hMatrixToy = (TH2F*)mPbPb_Matrix[i]->Clone();
    // 	hToy->SetName("hToy");
    // 	if (exp%100==0) if(printDebug) cout <<"Pseudo-experiment "<<exp<<endl;
    // 	for (int j=1;j<=hToy->GetNbinsX();j++) {
    // 	  double value = gRandom->Poisson(dPbPb_TrgComb[i]->GetBinContent(j));
    // 	  hToy->SetBinContent(j,value);
    // 	}

    // 	for (int j=1;j<=hMatrixToy->GetNbinsX();j++) {
    // 	  for (int k=1;k<=hMatrixToy->GetNbinsY();k++) {
    // 	    double value = gRandom->Gaus(mPbPb_Matrix[i]->GetBinContent(j,k),mPbPb_Matrix[i]->GetBinError(j,k));
    // 	    hMatrixToy->SetBinContent(j,k,value);
    // 	  }
    // 	}
    // 	//RooUnfoldBayes unfoldToy(response[i],hToy,2);
    // 	prior myPriorToy(hMatrixToy,hToy,0.0);
    // 	myPriorToy.unfold(hToy,1);
    // 	bayesianUnfold myUnfoldingToy(hMatrixToy,hPrior,0.0);
    // 	myUnfoldingToy.unfold(hToy,nBayesianIter);
    // 	RooUnfoldBinByBin unfoldToy2(response[i],hToy);
    // 	TH1F *hRecoTmp = (TH1F*) myUnfoldingToy.hPrior->Clone();
    // 	TH1F *hRecoTmp2 = (TH1F*) unfoldToy2.Hreco();

    // 	for (int j=1;j<=hRecoTmp->GetNbinsX();j++) {
    // 	  hTmp[j]->Fill(hRecoTmp->GetBinContent(j));
    // 	  hTmp2[j]->Fill(hRecoTmp2->GetBinContent(j));
    // 	}
    // 	delete hToy;
    // 	delete hRecoTmp;
    // 	delete hRecoTmp2;
    // 	delete hMatrixToy;
    //   }
    //   TF1 *f = new TF1("f","[0]*TMath::Gaus(x,[1],[2])");
    //   for (int j=1;j<=nbins_truth;j++)
    // 	{
    // 	  f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());

    // 	  if (hTmp[j]->GetMean()>0) {
    // 	    hTmp[j]->Fit("f","LL Q ");
    // 	    hTmp[j]->Fit("f","LL Q ");
    // 	    //	       cToy->SaveAs(Form("toy/cent-%d-pt-%.0f.gif",i,uhist[i]->hReco->GetBinCenter(j)));
    // 	    //     	       if(printDebug) cout <<j<<" "<<f->GetParameter(2)<<endl;
    // 	    uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
    // 	  }	       
    // 	  f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
    // 	  if (hTmp2[j]->GetMean()>0) {
    // 	    hTmp2[j]->Fit("f","LL Q ");
    // 	    hTmp2[j]->Fit("f","LL Q ");
    // 	    //cToy->SaveAs(Form("toy/cent2-%d-pt-%.0f.gif",i,uhist[i]->hReco->GetBinCenter(j)));
    // 	    //if(printDebug) cout <<j<<" "<<f->GetParameter(2)<<endl;
    // 	    uPbPb_BinByBin[i]->SetBinError(j,f->GetParameter(2));
    // 	  }	       
    // 	  delete hTmp[j];
    // 	  delete hTmp2[j];
    // 	}
    // }
 
  }
  
  // do the pp unfolding. 
  // Do Bin-by-bin
  if(printDebug) cout<<"doing bin by bin unfolding - PP"<<endl;

  TH1F* hBinByBinCorRawPP = (TH1F*)mPP_Response->ProjectionY();
  TH1F* hMCGenPP          = (TH1F*)mPP_Response->ProjectionX(); // gen
  hBinByBinCorRawPP->Divide(hMCGenPP);
  TF1 *fPP = new TF1("fPP","[0]+[1]*x");
  hBinByBinCorRawPP->Fit("fPP","LL ","",50,300);
  TH1F* hBinByBinCorPP = (TH1F*)hBinByBinCorRawPP->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCorPP");
  //TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCorPP");

  uPP_BinByBin = (TH1F*)dPP_Comb->Clone("uPP_BinByBin");
  uPP_BinByBin->Divide(hBinByBinCorPP);
  //      uPP_BinByBin[i] = (TH1F*) hMCReco->Clone("hRecoBinByBinPP");
  
  // Do unfolding
  //prior myPrior(mPP_Matrix[i],dPP_TrgComb[i],0);
  //myPrior.unfold(dPP_TrgComb[i],1);
  TH1F* hPriorPP = (TH1F*)hMCGenPP->Clone("hPriorPP");
  removeZero(hPriorPP);
  //hPrior->Scale(dPP_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
  
  //TH1F *hReweightedPP = (TH1F*)mPP_Response->ProjectionY("hReweightedPP");
  
  //bayesianUnfold myUnfoldingJECSys(mPP_Matrix[i],hPrior,0);
  //myUnfoldingJECSys.unfold(dPP_TrgComb[i]JECSys,nBayesianIter);
  //bayesianUnfold myUnfoldingSmearSys(mPP_Matrix[i],hPrior,0);
  //myUnfoldingSmearSys.unfold(dPP_TrgComb[i]SmearSys,nBayesianIter);
  
  bayesianUnfold myUnfolding(mPP_Matrix,hPriorPP,0);
  myUnfolding.unfold(dPP_Comb,BayesIter);
  
  delete hBinByBinCorRawPP;
  delete hMCGenPP;
  
  // Iteration Systematics
  for (int j=2;j<Iterations;j++){
    bayesianUnfold myUnfoldingSys(mPP_Matrix,hPriorPP,0);
    myUnfoldingSys.unfold(dPP_Comb,j);
    uPP_BayesianIter[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("uPP_BayesianIter%d",j));
    uPP_BayesianIter[j]->Print("base");
    uPP_BayesianIter[j]->SetTitle(Form("Unfolded pp Bayesian iteration %d",j));
  }
  
  uPP_Bayes        = (TH1F*) uPP_BayesianIter[BayesIter]->Clone("uPP_Bayes");
  //uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone("UnfoldedJECSys_cent%i",i));
  //uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone("UnfoldedSmearSys_cent%i",i));
  //uPP_BinByBin[i] = (TH1F*) unfold2.Hreco();
  uPP_BinByBin->SetName("uPP_BinByBin");  
  
  uPP_Bayes->SetTitle("Unfolded pp Bayes");
  uPP_BinByBin->SetTitle("Unfolded pp BinByBin");

  delete hPriorPP;
  
  // calculate the RAA 
  // Scale pp by 1./64 - sigma pp 
  // scale PbPb by 1./ncoll[i]
  // scale PbPb by 1./7.65 sigma inelastic
  // and just divide each other. 
  // RAA = (dsigma ^PbPb / dp_T deta) / (ncoll * dsigma ^PP / dp_T deta)
  
  TH1F *RAA_bayesian[nbins_cent+1];
  TH1F *RAA_binbybin[nbins_cent+1];
  TH1F *RAA_measured[nbins_cent+1];
  //uPP_Bayes->Scale(1./64);
  for(int i = 0;i<=nbins_cent;i++){

    //uPbPb_Bayes[i]->Scale(1./ncoll[i]);
    //uPbPb_Bayes[i]->Scale(1./7.65);
    RAA_bayesian[i] = (TH1F*)uPbPb_Bayes[i]->Clone(Form("RAA_bayesian_cent%d",i));
    RAA_binbybin[i] = (TH1F*)uPbPb_BinByBin[i]->Clone(Form("RAA_binbybin_cent%d",i));
    RAA_measured[i] = (TH1F*)dPbPb_TrgComb[i]->Clone(Form("RAA_measured_cent%d",i));
    RAA_bayesian[i]->Divide(uPP_Bayes);
    RAA_measured[i]->Divide(dPP_Comb);
    RAA_binbybin[i]->Divide(uPP_BinByBin);

    RAA_bayesian[i]->Scale(64./(ncoll[i]*7.65));
    RAA_measured[i]->Scale(64./(ncoll[i]*7.65));
    RAA_binbybin[i]->Scale(64./(ncoll[i]*7.65));
    
    if(i!=6){
      RAA_bayesian[i]->SetTitle(Form("RAA ak%s%d%s bayesian unfolded %2.0f - %2.0f cent",algo,radius,jet_type,5*boundaries_cent[i], 5*boundaries_cent[i+1]));
      RAA_measured[i]->SetTitle(Form("RAA ak%s%d%s measured unfolded %2.0f - %2.0f cent",algo,radius,jet_type,5*boundaries_cent[i], 5*boundaries_cent[i+1]));
      RAA_binbybin[i]->SetTitle(Form("RAA ak%s%d%s binbybin unfolded %2.0f - %2.0f cent",algo,radius,jet_type,5*boundaries_cent[i], 5*boundaries_cent[i+1]));
    }else 
      RAA_bayesian[i]->SetTitle(Form("RAA ak%s%d%s bayesian unfolded 0 - 200 cent",algo,radius,jet_type));
      RAA_measured[i]->SetTitle(Form("RAA ak%s%d%s measured unfolded 0 - 200 cent",algo,radius,jet_type));
      RAA_binbybin[i]->SetTitle(Form("RAA ak%s%d%s binbybin unfolded 0 - 200 cent",algo,radius,jet_type));
  }
  
  // think if there are any other systematic checks to be performed. 
  // so i have iteration systematics, what about Unfolding closure using the MC. For that i need to take in half of the MC as data and then unfold it using the same response matrix. I think it would be better if i create a new macro for that. Or coming to think of it, i have all the unfolding here so might as well just make a copy of it. 

  // unfolding for the MC closure test. 

  // create the necessary histograms.  
  TH1F* uPbPb_MC_Bayes[nbins_cent+1];
  TH1F* uPbPb_MC_BinByBin[nbins_cent+1];
  TH1F* uPbPb_MC_BayesianIter[nbins_cent+1][Iterations];

  // first for PbPb
  for (int i=0;i<=nbins_cent;i++) {

    // Do Bin-by-bin
    if(printDebug) cout<<"doing bin by bin unfolding for PbPb MC closure test for centrality = "<<i<<endl;

    TH1F* hBinByBinCorRaw = (TH1F*)mPbPb_mcclosure_Response[i]->ProjectionY();
    TH1F* hMCGen          = (TH1F*)mPbPb_mcclosure_Response[i]->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",50,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    //TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));

    //if(printDebug) cout<<"passed here after binbybincorrraw"<<endl;
    //dPbPb_TrgComb[i]->Print("base");
    uPbPb_MC_BinByBin[i] = (TH1F*)mPbPb_mcclosure_data[i]->Clone(Form("uPbPb_MC_BinByBin_cent%d",i));
    //if(printDebug) cout<<"maybe here?"<<endl;
    uPbPb_MC_BinByBin[i]->Divide(hBinByBinCor);
    //uPbPb_BinByBin[i] = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
    //if(printDebug) cout<<"passed A"<<endl;
    
    // Do unfolding
    //prior myPrior(mPbPb_Matrix[i],dPbPb_TrgComb[i],0);
    //myPrior.unfold(dPbPb_TrgComb[i],1);
    TH1F* hPriorMC = (TH1F*)hMCGen->Clone("hPriorMC");
    //mPbPb_mcclosure_Gen[i] = (TH1F*)hMCGen->Clone(Form("hPbPb_mcclosure_gen_R%d_n20_eta_p20_cent%d",radius,i));

    removeZero(hPriorMC);
    //hPrior->Scale(dPbPb_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
    //if(printDebug) cout<<"passed B"<<endl;
    //TH1F *hReweightedPbPb = (TH1F*)mPbPb_Response[i]->ProjectionY()->Clone(Form("hReweightedPbPb_cent%d",i));

    //bayesianUnfold myUnfoldingJECSys(mPbPb_Matrix[i],hPrior,0);
    //myUnfoldingJECSys.unfold(dPbPb_TrgComb[i]JECSys,nBayesianIter);
    //bayesianUnfold myUnfoldingSmearSys(mPbPb_Matrix[i],hPrior,0);
    //myUnfoldingSmearSys.unfold(dPbPb_TrgComb[i]SmearSys,nBayesianIter);

    bayesianUnfold myUnfoldingMC(mPbPb_mcclosure_Matrix[i],hPriorMC,0);
    myUnfoldingMC.unfold(mPbPb_mcclosure_data[i],BayesIter);
    
    mPbPb_mcclosure_data[i]->SetTitle(Form("PbPb MC closure test data cent%d",i));
    mPbPb_mcclosure_data[i]->SetName(Form("mPbPb_mclosure_data_cent%d",i));

    if(printDebug) cout <<"Unfolding bin "<<i<<endl;

    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<Iterations;j++){
      bayesianUnfold myUnfoldingSys(mPbPb_mcclosure_Matrix[i],hPriorMC,0);
      myUnfoldingSys.unfold(mPbPb_mcclosure_data[i],j);
      uPbPb_MC_BayesianIter[i][j]  = (TH1F*) myUnfoldingSys.hUnfolded->Clone(Form("uPbPb_MC_BayesianIter%d_cent%d",j,i));
      uPbPb_MC_BayesianIter[i][j] ->Print("base");
      if(i<6) uPbPb_MC_BayesianIter[i][j]->SetTitle(Form("Unfolded PbPb MC closure test Bayesian iteration %d with %2.0f - %2.0f cent",j,5*boundaries_cent[i],5*boundaries_cent[i+1]));
      else uPbPb_MC_BayesianIter[i][j]->SetTitle(Form("Unfolded PbPb MC closure test Bayesian iteration %d with 0-200 cent",j));
    }
    if(printDebug) cout<<"passed iteration sys"<<endl;
    uPbPb_MC_Bayes[i]        = (TH1F*) uPbPb_MC_BayesianIter[i][BayesIter]->Clone(Form("uPbPb_MC_Bayes_cent%i",i));
    //uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJECSys_cent%i",i));
    //uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    uPbPb_MC_BinByBin[i]->SetName(Form("uPbPb_MC_BinByBin_cent%i",i));

    if(i==6){
      uPbPb_MC_Bayes[i]->SetTitle("Unfolded PbPb MC Closure test Bayes 0-200 cent");
      uPbPb_MC_BinByBin[i]->SetTitle("Unfolded PbPb MC closure test BinByBin 0-200 cent");
    }else {
      uPbPb_MC_Bayes[i]->SetTitle(Form("Unfolded PbPb MC closure test Bayes %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]));
      uPbPb_MC_BinByBin[i]->SetTitle(Form("Unfolded PbPb MC closure test BinByBin %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]));
    }

    delete hPriorMC;
    
  }

  
  // do pp mc unfolding closure test 
  // Do Bin-by-bin
  if(printDebug) cout<<"doing bin by bin unfolding for MC closure test - PP"<<endl;

  TH1F *uPP_MC_BinByBin, *uPP_MC_Bayes;
  TH1F *uPP_MC_BayesianIter[Iterations];

  TH1F* hBinByBinCorRawPPMC = (TH1F*)mPP_mcclosure_Response->ProjectionY();
  TH1F* hMCGenPPMC          = (TH1F*)mPP_mcclosure_Response->ProjectionX(); // gen
  hBinByBinCorRawPPMC->Divide(hMCGenPPMC);
  TF1 *fPPMC = new TF1("fPPMC","[0]+[1]*x");
  hBinByBinCorRawPPMC->Fit("fPPMC","LL ","",50,300);
  TH1F* hBinByBinCorPPMC = (TH1F*)hBinByBinCorRawPPMC->Clone("hBinByBinCorPPMC");//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCorPP");
  //TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCorPP");

  uPP_MC_BinByBin = (TH1F*)mPP_mcclosure_data->Clone("uPP_MC_BinByBin");
  uPP_MC_BinByBin->Divide(hBinByBinCorPPMC);
  //      uPP_BinByBin[i] = (TH1F*) hMCReco->Clone("hRecoBinByBinPP");
  
  // Do unfolding
  //prior myPrior(mPP_Matrix[i],dPP_TrgComb[i],0);
  //myPrior.unfold(dPP_TrgComb[i],1);
  TH1F* hPriorPPMC = (TH1F*)hMCGenPPMC->Clone("hPriorPPMC");
  removeZero(hPriorPPMC);
  //hPrior->Scale(dPP_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
  
  //TH1F *hReweightedPP = (TH1F*)mPP_Response->ProjectionY("hReweightedPP");
  
  //bayesianUnfold myUnfoldingJECSys(mPP_Matrix[i],hPrior,0);
  //myUnfoldingJECSys.unfold(dPP_TrgComb[i]JECSys,nBayesianIter);
  //bayesianUnfold myUnfoldingSmearSys(mPP_Matrix[i],hPrior,0);
  //myUnfoldingSmearSys.unfold(dPP_TrgComb[i]SmearSys,nBayesianIter);
  
  bayesianUnfold myUnfoldingPPMC(mPP_mcclosure_Matrix,hPriorPPMC,0);
  myUnfoldingPPMC.unfold(mPP_mcclosure_data,BayesIter);
  
  mPP_mcclosure_data->SetName("mPP_mcclosure_data");
  mPP_mcclosure_data->SetTitle("PP mcclosure data");
  
  delete hBinByBinCorRawPPMC;
  delete hMCGenPPMC;
  
  // Iteration Systematics
  for (int j=2;j<Iterations;j++){
    bayesianUnfold myUnfoldingSys(mPP_mcclosure_Matrix,hPriorPPMC,0);
    myUnfoldingSys.unfold(mPP_mcclosure_data,j);
    uPP_MC_BayesianIter[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("uPP_MC_BayesianIter%d",j));
    uPP_MC_BayesianIter[j]->Print("base");
    uPP_MC_BayesianIter[j]->SetTitle(Form("Unfolded PP MC closure test iteration %d",j));
  }
  
  uPP_MC_Bayes        = (TH1F*) uPP_MC_BayesianIter[BayesIter]->Clone("uPP_MC_Bayes");
  //uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone("UnfoldedJECSys_cent%i",i));
  //uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone("UnfoldedSmearSys_cent%i",i));
  //uPP_BinByBin[i] = (TH1F*) unfold2.Hreco();
  uPP_MC_BinByBin->SetName("uPP_MC_BinByBin");  

  uPP_MC_Bayes->SetTitle("Unfolded pp MC closure test Bayes");
  uPP_MC_BinByBin->SetTitle("Unfolded pp MC closure test BinByBin");

  delete hPriorPPMC;
  
  

  // write it to the output file
  
  cout<<"writing to output file"<<endl;
    
  TFile fout(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_unfo_mcclosure_oppside_ak%s%d%s_%d_test.root",algo,radius,jet_type,date.GetDate()),"RECREATE");
  fout.cd();

  for(int i = 0;i<=nbins_cent;i++){
    
    uPbPb_Bayes[i]->Write();
    uPbPb_BinByBin[i]->Write();
    dPbPb_TrgComb[i]->Write();
    dPbPb_Trg80[i]->Write();
    dPbPb_Trg65[i]->Write();
    dPbPb_Trg55[i]->Write();

    mPbPb_ResponseNorm[i]->Write();
    mPbPb_Response[i]->Write();
    mPbPb_Matrix[i]->Write();

    mPbPb_mcclosure_ResponseNorm[i]->Write();
    mPbPb_mcclosure_Response[i]->Write();
    mPbPb_mcclosure_Matrix[i]->Write();

    mPbPb_Gen[i]->Write();
    mPbPb_Reco[i]->Write();

    uPbPb_MC_Bayes[i]->Write();
    uPbPb_MC_BinByBin[i]->Write();

    mPbPb_mcclosure_data[i]->Write();
    mPbPb_mcclosure_gen[i]->Write();
    //RAA_measured[i]->Write();
    //RAA_binbybin[i]->Write();
    //RAA_bayesian[i]->Write();
    
    for(int j = 2;j<Iterations;j++){
      uPbPb_BayesianIter[i][j]->Write();
      uPbPb_MC_BayesianIter[i][j]->Write();
    }
  }//cent bin loop
  
  dPP_Comb->Write();
  mPP_ResponseNorm->Write();
  mPP_Response->Write();
  mPP_mcclosure_ResponseNorm->Write();
  mPP_mcclosure_Response->Write();
  mPP_Reco->Write();
  mPP_Gen->Write();
  mPP_Matrix->Write();
  mPP_mcclosure_Matrix->Write();
  mPP_mcclosure_data->Write();
  
  for(int i= 2;i<Iterations;i++){
    uPP_BayesianIter[i]->Write();
    //uPP_MC_BayesianIter[i]->Write();
  }

  uPP_Bayes->Write();
  //uPP_MC_Bayes->Write();
  uPP_BinByBin->Write();
  //uPP_MC_BinByBin->Write();
  
  //fout.Write();
  fout.Close();
  
  timer.Stop();
  if(printDebug) cout<<"CPU time = "<<timer.CpuTime()<<endl;
  if(printDebug) cout<<"Real tile = "<<timer.RealTime()<<endl;
  
  
}
