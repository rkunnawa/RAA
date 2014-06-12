// Raghav Kunnawalkam Elayavalli
// June 10th 2014
// CERN

// RAA - analyze macro. following the steps of having a separate macro to do the work. 

// this will be taken mostly from Unfold_RAA_V0.C from MIT. Im just keeping it here so that it will be easier to run. 
// importnat point in the code is that in the loops for pbpb, the array index corrsponding to nbins_cent represents the 0-100 centrality bin. 
// PbPb and pp are separated here - as they should be. 

#include <iostream>
#include <stdio.h>

#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/prior.h"
#include "/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Headers/bayesianUnfold.h"

  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};


void RAA_analysis(int radius = 3, char* algo = "Pu"){

  TStopwatch timer; 
  timer.Start();
  

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  // get the data and mc histograms from the output of the read macro. 
  
  TDatime date;//this is just here to get them to run optimized. 
  
  TFile* fData_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_data_ak%d_%s_%d_chMax_12003cut.root",radius,algo,date.GetDate()));
  TFile* fMC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_mc_ak%d_%s_%d_chMax_12003cut.root",radius,algo,date.GetDate()));
  // need to make sure that the file names are in prefect order so that i can run them one after another. 
  // for the above condition, i might have to play with the date stamp. 
  
  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  
  // histogram declarations with the following initial appendage: d - Data, m - MC, u- Unfolded
  
  TH1F *dPbPb_TrgComb[nbins_cent+1], *dPbPb_Comb[nbins_cent+1], *dPbPb_Trg80[nbins_cent+1], *dPbPb_Trg65[nbins_cent+1], *dPbPb_Trg55[nbins_cent+1], *dPbPb_1[nbins_cent+1], *dPbPb_2[nbins_cent+1], *dPbPb_3[nbins_cent+1], *dPbPb_80[nbins_cent+1], *dPbPb_65[nbins_cent+1], *dPbPb_55[nbins_cent+1];
  
  TH1F *mPbPb_Gen[nbins_cent+1], *mPbPb_Reco[nbins_cent+1];
  TH2F *mPbPb_Matrix[nbins_cent+1], *mPbPb_Response[nbins_cent+1], *mPbPb_ResponseNorm[nbins_cent+1];
  
  TH1F *dPP_1, *dPP_2, *dPP_3, *dPP_Comb;
  
  TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix, *mPP_Response, *mPP_ResponseNorm;
  
  static const int Iterations = 20;
  static const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1][Iterations], *uPbPb_BinByBin[nbins_cent+1];
  
  // would be better to read in the histograms and rebin them. come to think of it, it would be better to have them already rebinned (and properly scaled - to the level of differential cross section in what ever barns (inverse micro barns) but keep it consistent) from the read macro. 
  
  
  for(int i = 0;i<=nbins_cent;i++){
    
    dPbPb_TrgComb[i] = (TH1F*)fData_in->Get(Form("hpbpb_TrgComb_%d",i));
    dPbPb_Trg80[i] = (TH1F*)fData_in->Get(Form("hpbpb_Trg80_%d",i));
    dPbPb_Trg65[i] = (TH1F*)fData_in->Get(Form("hpbpb_Trg65_%d",i));
    dPbPb_Trg55[i] = (TH1F*)fData_in->Get(Form("hpbpb_Trg55_%d",i));

    dPbPb_Comb[i] = (TH1F*)fData_in->Get(Form("hpbpbComb_%d",i));
    dPbPb_1[i] = (TH1F*)fData_in->Get(Form("hpbpb1_%d",i));
    dPbPb_2[i] = (TH1F*)fData_in->Get(Form("hpbpb2_%d",i));
    dPbPb_3[i] = (TH1F*)fData_in->Get(Form("hpbpb3_%d",i));

    dPbPb_80[i] = (TH1F*)fData_in->Get(Form("hpbpb_80_%d",i));
    dPbPb_65[i] = (TH1F*)fData_in->Get(Form("hpbpb_65_%d",i));
    dPbPb_55[i] = (TH1F*)fData_in->Get(Form("hpbpb_55_%d",i));

  }

  for(int i = 0;i<=nbins_cent;i++){

    mPbPb_Gen[i] = (TH1F*)fMc_in->Get(Form("hpbpb_gen_%d",i));
    mPbPb_Reco[i] = (TH1F*)fMc_in->Get(Form("hpbpb_reco_%d",i));
    mPbPb_Matrix[i] = (TH2F*)fMc_in->Get(Form("hpbpb_matrix_%d",i));
    //mPbPb_Response[i] = new TH2F(Form("mPbPb_Response_%d",i),"Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    //mPbPb_ResponseNorm[i] = new TH2F(Form("mPbPb_ResponseNorm_%d",i),"Normalized Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    
  }

  dPP_1 = (TH1F*)fData_in->Get("hpp1");
  dPP_2 = (TH1F*)fData_in->Get("hpp2");
  dPP_3 = (TH1F*)fData_in->Get("hpp3");
  dPPComb = (TH1F*)fData_in->Get("hppComb");

  mPP_Gen = (TH1F*)fMc_in->Get("hpp_gen");
  mPP_Reco = (TH1F*)fMc_in->Get("hpp_gen");
  mPP_Matrix = (TH2F*)fMc_in->Get("hpp_gen");
  mPP_Response = (TH2F*)fMc_in->Get("hpp_gen");


  // make the response matrix.
  // here since we dont have the simple nature of the uhist histograms which makes debugging a pain in the ass, we have to run the response matrix and unfolding separately for PbPb and pp which makes code ugly but easy to debugg at the same time. 

  // response matrix and unfolding for PbPb 
  for(int i = 0;i<=nbins_cent;i++){
    TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
    f->SetParameters(1e10,-8.8,40);
    TH1F *hGenSpectraCorr = (TH1F*)mPbPb_Matrix[i]->ProjectionX()->Clone(Form("hGenSpectraCorr_cent%d",i));
    hGenSpectraCorr->Fit("f"," ");
    hGenSpectraCorr->Fit("f","","");
    hGenSpectraCorr->Fit("f","LL");
    TH1F *fHist = functionHist(f,hGenSpectraCorr,Form("fHist_cent%d",i));// that the function that you get from the fitting 
    hGenSpectraCorr->Divide(fHist);
    for (int y=1;y<=mPbPb_Matrix[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_Matrix[i]->GetNbinsX();x++) {
	if (mPbPb_Matrix[i]->GetBinContent(x,y)<=0*mPbPb_Matrix[i]->GetBinError(x,y)) {
	  mPbPb_Matrix[i]->SetBinContent(x,y,0);
	  mPbPb_Matrix[i]->SetBinError(x,y,0);
	}
	sum+=mPbPb_Matrix[i]->GetBinContent(x,y);
      }
      
      for (int x=1;x<=mPbPb_Matrix[i]->GetNbinsX();x++) {	   
	double ratio = 1;
	if (hGenSpectraCorr->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorr->GetBinContent(x);
	mPbPb_Matrix[i]->SetBinContent(x,y,mPbPb_Matrix[i]->GetBinContent(x,y)*ratio);
	mPbPb_Matrix[i]->SetBinError(x,y,mPbPb_Matrix[i]->GetBinError(x,y)*ratio);
      }
    }
    //uhist[i]->hMatrix->Smooth(0);

    // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  
    
    mPbPb_Response[i] = (TH2F*)mPbPb_Matrix[i]->Clone(Form("mPbPb_Response_cent%d",i));
    TH1F *hProj = (TH1F*)mPbPb_Response[i]->ProjectionY(Form("hProj_cent%d",i));


    for (int y=1;y<=mPbPb_Response[i]->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=mPbPb_Response[i]->GetNbinsX();x++) {
	if (mPbPb_Response[i]->GetBinContent(x,y)<=1*mPbPb_Response[i]->GetBinError(x,y)) {
	  // in the above if loop, kurt has 1*error and my old has 0*error
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
	if (mPbPb_ResponseNorm[i]->GetBinContent(x,y)<=0*mPbPb_ResponseNorm[i]->GetBinError(x,y)) {
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
  
  
  
  // response matrix for pp.  
  
  TF1 *fpp = new TF1("fpp","[0]*pow(x+[2],[1])");
  fpp->SetParameters(1e10,-8.8,40);
  TH1F *hGenSpectraCorrPP = (TH1F*)mPP_Matrix->ProjectionX()->Clone("hGenSpectraCorrPP");
  hGenSpectraCorrPP->Fit("f"," ");
  hGenSpectraCorrPP->Fit("f","","");
  hGenSpectraCorrPP->Fit("f","LL");
  TH1F *fHistPP = functionHist(fpp,hGenSpectraCorrPP,"fHistPP");// that the function that you get from the fitting 
  hGenSpectraCorrPP->Divide(fHistPP);
  for (int y=1;y<=mPP_Matrix->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_Matrix->GetNbinsX();x++) {
      if (mPP_Matrix->GetBinContent(x,y)<=0*mPP_Matrix->GetBinError(x,y)) {
	mPP_Matrix->SetBinContent(x,y,0);
	mPP_Matrix->SetBinError(x,y,0);
      }
      sum+=mPP_Matrix->GetBinContent(x,y);
    }
    
    for (int x=1;x<=mPP_Matrix->GetNbinsX();x++) {	   
      double ratio = 1;
      if (hGenSpectraCorrPP->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorrPP->GetBinContent(x);
      mPP_Matrix->SetBinContent(x,y,mPP_Matrix->GetBinContent(x,y)*ratio);
      mPP_Matrix->SetBinError(x,y,mPP_Matrix->GetBinError(x,y)*ratio);
    }
  }
  //uhist[i]->hMatrix->Smooth(0);
  
  // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  
  
  mPP_Response = (TH2F*)mPP_Matrix->Clone("mPbPb_Response");
  TH1F *hProjPP = (TH1F*)mPP_Response->ProjectionY("hProjPP");
  
  
  for (int y=1;y<=mPP_Response->GetNbinsY();y++) {
    double sum=0;
    for (int x=1;x<=mPP_Response->GetNbinsX();x++) {
      if (mPP_Response->GetBinContent(x,y)<=1*mPP_Response->GetBinError(x,y)) {
	// in the above if loop, kurt has 1*error and my old has 0*error
	mPP_Response->SetBinContent(x,y,0);
	mPP_Response->SetBinError(x,y,0);
      }
      sum+=mPP_Response->GetBinContent(x,y);
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
      if (mPbPb_ResponseNorm[i]->GetBinContent(x,y)<=0*mPbPb_ResponseNorm[i]->GetBinError(x,y)) {
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
  
  
  
  
  // do the unfolding - including the iteration systematics (ofcourse)
  
  
  
  // think if there are any other systematic checks to be performed. 
  

  
  // write it to the output file
  
  

  TDatime date;

  TFile fout(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_unfo_ak%d_%s_%d_chMax_12003cut.root",radius,algo,data.GetDate()),"RECREATE");
  fout.cd();

  fout.Write();
  fout.Close();

  timer.Stop();
  cout<<"CPU time = "<<timer.CpuTime()<<endl;
  cout<<"Real tile = "<<timer.RealTime()<<endl;
  

}
