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
  
  TFile* fData_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_data_ak%d_%s_cent%d_chMax_12003cut.root",radius,algo,date.GetDate()));
  TFile* fMC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_mc_ak%d_%s_cent%d_chMax_12003cut.root",radius,algo,date.GetDate()));
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
  
  static const int Iterations = 20; //for unfolding systematics. 
  static const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1], *uPbPb_BinByBin[nbins_cent+1]; 
  TH1F *uPbPb_BayesianIter[nbins_cent+1][Iterations];
  
  // would be better to read in the histograms and rebin them. come to think of it, it would be better to have them already rebinned (and properly scaled - to the level of differential cross section in what ever barns (inverse micro barns) but keep it consistent) from the read macro. 
  
  for(int i = 0;i<=nbins_cent;i++){
    
    dPbPb_TrgComb[i] = (TH1F*)fData_in->Get(Form("hpbpb_TrgComb_cent%d",i));
    dPbPb_Trg80[i] = (TH1F*)fData_in->Get(Form("hpbpb_Trg80_cent%d",i));
    dPbPb_Trg65[i] = (TH1F*)fData_in->Get(Form("hpbpb_Trg65_cent%d",i));
    dPbPb_Trg55[i] = (TH1F*)fData_in->Get(Form("hpbpb_Trg55_cent%d",i));

    dPbPb_Comb[i] = (TH1F*)fData_in->Get(Form("hpbpbComb_cent%d",i));
    dPbPb_1[i] = (TH1F*)fData_in->Get(Form("hpbpb1_cent%d",i));
    dPbPb_2[i] = (TH1F*)fData_in->Get(Form("hpbpb2_cent%d",i));
    dPbPb_3[i] = (TH1F*)fData_in->Get(Form("hpbpb3_cent%d",i));

    dPbPb_80[i] = (TH1F*)fData_in->Get(Form("hpbpb_80_cent%d",i));
    dPbPb_65[i] = (TH1F*)fData_in->Get(Form("hpbpb_65_cent%d",i));
    dPbPb_55[i] = (TH1F*)fData_in->Get(Form("hpbpb_55_cent%d",i));

  }

  for(int i = 0;i<=nbins_cent;i++){

    mPbPb_Gen[i] = (TH1F*)fMc_in->Get(Form("hpbpb_gen_cent%d",i));
    mPbPb_Reco[i] = (TH1F*)fMc_in->Get(Form("hpbpb_reco_cent%d",i));
    mPbPb_Matrix[i] = (TH2F*)fMc_in->Get(Form("hpbpb_matrix_cent%d",i));
    //mPbPb_Response[i] = new TH2F(Form("mPbPb_Response_cent%d",i),"Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    //mPbPb_ResponseNorm[i] = new TH2F(Form("mPbPb_ResponseNorm_cent%d",i),"Normalized Response Matrix",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    
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
    //mPbPb_Matrix[i]->Smooth(0);

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
  //mPbPb_Matrix[i]->Smooth(0);
  
  // Ok major differences here between my code and Kurt in b-jet Tools under Unfold - lines 469 and above.  
  
  mPP_Response = (TH2F*)mPP_Matrix->Clone("mPP_Response");
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
    
    for (int x=1;x<=mPP_Response->GetNbinsX();x++) {  	
      if (sum==0) continue;
      double ratio = 1;
      //if(dPbPb_TrgComb[i]->GetBinContent(y)==0) ratio = 1e-100/sum;
      // else ratio = dPbPb_TrgComb[i]->GetBinContent(y)/sum
      ratio = 1./sum;
      if (hProj->GetBinContent(y)==0) ratio = 1e-100/sum;
      else ratio = hProj->GetBinContent(y)/sum;
      mPP_Response->SetBinContent(x,y,mPP_Response->GetBinContent(x,y)*ratio);
      mPP_Response->SetBinError(x,y,mPP_Response->GetBinError(x,y)*ratio);
    }
  }
  
  mPP_ResponseNorm = (TH2F*)mPP_Matrix[i]->Clone(Form("mPP_ResponseNorm_cent%d",i));
  for (int x=1;x<=mPP_ResponseNorm->GetNbinsX();x++) {
    double sum=0;
    for (int y=1;y<=mPP_ResponseNorm->GetNbinsY();y++) {
      if (mPP_ResponseNorm->GetBinContent(x,y)<=0*mPP_ResponseNorm->GetBinError(x,y)) {
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
  
  // do the unfolding - including the iteration systematics (ofcourse). Similar to the above version, we have to create 2 separate unfolding for PbPb and pp. 

  // first for PbPb
  for (int i=0;i<=nbins_cent;i++) {

    // Do Bin-by-bin
    cout<<"doing bin by bin unfolding for PbPb for centrality = "<<i<<endl;

    TH1F* hBinByBinCorRaw = (TH1F*)mPbPb_Response[i]->ProjectionY();
    TH1F* hMCGen          = (TH1F*)mPbPb_Response[i]->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",50,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    //TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));

    uPbPb_BinByBin[i] = (TH1F*) dPbPb_TrgComb[i]->Clone(Form("uPbPb_BinByBin_cent%d",i));
    uPbPb_BinByBin[i]->Divide(hBinByBinCor);
    //      uPbPb_BinByBin[i] = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
    
    // Do unfolding
    //prior myPrior(mPbPb_Matrix[i],dPbPb_TrgComb[i],0);
    //myPrior.unfold(dPbPb_TrgComb[i],1);
    TH1F *hPrior = (TH1F*)hMCGen->Clone("hPrior");
    removeZero(hPrior);
    //hPrior->Scale(dPbPb_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
    
    TH1F *hReweighted = (TH1F*)mPbPb_Response[i]->ProjectionY(Form("hReweighted_cent%d",i));

    //bayesianUnfold myUnfoldingJECSys(mPbPb_Matrix[i],hPrior,0);
    //myUnfoldingJECSys.unfold(dPbPb_TrgComb[i]JECSys,nBayesianIter);
    //bayesianUnfold myUnfoldingSmearSys(mPbPb_Matrix[i],hPrior,0);
    //myUnfoldingSmearSys.unfold(dPbPb_TrgComb[i]SmearSys,nBayesianIter);

    bayesianUnfold myUnfolding(mPbPb_Matrix[i],hPrior,0);
    myUnfolding.unfold(dPbPb_TrgComb[i],BayesIter);

    cout <<"Unfolding bin "<<i<<endl;

    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<Iterations;j++){
      bayesianUnfold myUnfoldingSys(mPbPb_Matrix[i],hPrior,0);
      myUnfoldingSys.unfold(dPbPb_TrgComb[i],j);
      uPbPb_BayesianIter[i][j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("uPbPb_IterSys%d_cent%d",j,i));
      uPbPb_BayesianIter[i][j] ->Print("base");
    }

    uPbPb_Bayes[i]        = (TH1F*) uPbPb_BayesianIterSys[i][BayesIter]->Clone(Form("uPbPb_Iter_cent%i",i));
    //uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJECSys_cent%i",i));
    //uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    //uPbPb_BinByBin[i] = (TH1F*) unfold2.Hreco();
    uPbPb_BinByBin[i]->SetName(Form("uPbPb_BinByBin_cent%i",i));

    /*
    // Do unfolding
    prior myPrior(mPbPb_Matrix[i],dPbPb_TrgComb[i],0.0);
    myPrior.unfold(dPbPb_TrgComb[i],1);
    TH1F *hPrior = (TH1F*)mPbPb_Matrix[i]->ProjectionX()->Clone(Form("hPrior_cent%d",i));
    hPrior->Scale(dPbPb_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
    bayesianUnfold myUnfoldingJECSys(mPbPb_Matrix[i],hPrior,0);
    myUnfoldingJECSys.unfold(dPbPb_TrgComb[i]JECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(mPbPb_Matrix[i],hPrior,0);
    myUnfoldingSmearSys.unfold(dPbPb_TrgComb[i]SmearSys,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;
    // Iteration Systematics
    for (int j=2;j<7;j++)
    {
    bayesianUnfold myUnfoldingSys(mPbPb_Matrix[i],hPrior,0);
    myUnfoldingSys.unfold(dPbPb_TrgComb[i],j);
    uhist[i]->hRecoIterSys[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
    }
		
    // Do Bin-by-bin
    TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY(); 
    TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",100,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    //      TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    delete hBinByBinCorRaw,hMCGen;
    uPbPb_BinByBin[i] = (TH1F*) dPbPb_TrgComb[i]->Clone(Form("hRecoBinByBin_cent%d",i));
    uPbPb_BinByBin[i]->Divide(hBinByBinCor);
    //      uPbPb_BinByBin[i] = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
		
    uhist[i]->hReco         = (TH1F*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJeCSys_cent%i",i));
    uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    //uPbPb_BinByBin[i] = (TH1F*) unfold2.Hreco();
    uPbPb_BinByBin[i]->SetName(Form("UnfoldedBinByBin_cent%i",i));
    */
    

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
    // 	if (exp%100==0) cout <<"Pseudo-experiment "<<exp<<endl;
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
    // 	    //     	       cout <<j<<" "<<f->GetParameter(2)<<endl;
    // 	    uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
    // 	  }	       
    // 	  f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
    // 	  if (hTmp2[j]->GetMean()>0) {
    // 	    hTmp2[j]->Fit("f","LL Q ");
    // 	    hTmp2[j]->Fit("f","LL Q ");
    // 	    //cToy->SaveAs(Form("toy/cent2-%d-pt-%.0f.gif",i,uhist[i]->hReco->GetBinCenter(j)));
    // 	    //cout <<j<<" "<<f->GetParameter(2)<<endl;
    // 	    uPbPb_BinByBin[i]->SetBinError(j,f->GetParameter(2));
    // 	  }	       
    // 	  delete hTmp[j];
    // 	  delete hTmp2[j];
    // 	}
    // }
 
  }

  // do the pp unfolding. 
      // Do Bin-by-bin
    cout<<"doing bin by bin unfolding - PP"<<endl;

    TH1F* hBinByBinCorRaw = (TH1F*)mPP_Response->ProjectionY();
    TH1F* hMCGen          = (TH1F*)mPP_Response->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",50,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCorPP");
    //TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCorPP");

    uPP_BinByBin = (TH1F*) dPPComb->Clone("uPP_BinByBinPP");
    uPP_BinByBin->Divide(hBinByBinCor);
    //      uPP_BinByBin[i] = (TH1F*) hMCReco->Clone("hRecoBinByBinPP");
    
    // Do unfolding
    //prior myPrior(mPP_Matrix[i],dPP_TrgComb[i],0);
    //myPrior.unfold(dPP_TrgComb[i],1);
    TH1F *hPrior = (TH1F*)hMCGen->Clone("hPrior");
    removeZero(hPrior);
    //hPrior->Scale(dPP_TrgComb[i]->Integral(0,1000)/hPrior->Integral(0,1000));
    
    TH1F *hReweighted = (TH1F*)mPP_Response->ProjectionY("hReweightedPP");

    //bayesianUnfold myUnfoldingJECSys(mPP_Matrix[i],hPrior,0);
    //myUnfoldingJECSys.unfold(dPP_TrgComb[i]JECSys,nBayesianIter);
    //bayesianUnfold myUnfoldingSmearSys(mPP_Matrix[i],hPrior,0);
    //myUnfoldingSmearSys.unfold(dPP_TrgComb[i]SmearSys,nBayesianIter);

    bayesianUnfold myUnfolding(mPP_Matrix,hPrior,0);
    myUnfolding.unfold(dPPComb,BayesIter);

    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<Iterations;j++){
      bayesianUnfold myUnfoldingSys(mPP_Matrix,hPrior,0);
      myUnfoldingSys.unfold(dPPComb,j);
      uPP_BayesianIter[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone("uPP_IterSys%d",j));
      uPP_BayesianIter[j] ->Print("base");
    }

    uPP_Bayes        = (TH1F*) uPP_BayesianIterSys[BayesIter]->Clone("uPP_IterPP");
    //uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone("UnfoldedJECSys_cent%i",i));
    //uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone("UnfoldedSmearSys_cent%i",i));
    //uPP_BinByBin[i] = (TH1F*) unfold2.Hreco();
    uPP_BinByBin->SetName("uPP_BinByBin");
  
  // think if there are any other systematic checks to be performed. 
  

  
  // write it to the output file
  
  

  TDatime date;

  TFile fout(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_unfo_ak%d_%s_cent%d_chMax_12003cut.root",radius,algo,data.GetDate()),"RECREATE");
  fout.cd();

  fout.Write();
  fout.Close();

  timer.Stop();
  cout<<"CPU time = "<<timer.CpuTime()<<endl;
  cout<<"Real tile = "<<timer.RealTime()<<endl;
  

}
