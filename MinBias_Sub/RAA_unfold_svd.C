// Raghav Kunnawalkam Elayavalli
// April 22th 2015
// CERN
// raghav.k.e AT CERN DOT CH

// Macro to perform svd unfolding.


#include <iostream>
#include <stdio.h>

#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldResponse.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldBayes.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldSvd.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldBinByBin.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/prior.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/bayesianUnfold.h"

#include "TStopwatch.h"


static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};

// divide by bin width
void divideBinWidth(TH1 *h)
{
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();++i){
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

using namespace std;

void RAA_unfold_svd(int radius = 3,
		    char* algo = (char*) "Pu",
		    char *jet_type = (char*) "PF",
		    int unfoldingCut = 60,
		    char* etaWidth = (char*) "20_eta_20",
		    double deltaEta = 4.0)
{

  TStopwatch timer; 
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  bool printDebug = true;

  // get the data and mc histograms from the output of the read macro. 
  
  TDatime date;//this is just here to get them to run optimized. 

  
  // Pawan's files:
  TFile * fPbPb_in = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_TTree_PbPb_Data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius));
  TFile * fPbPb_MC_in = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_TTree_PbPb_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius));
  TFile * fPP_in = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_TTree_PP_data_MC_spectra_residualFactor_finebins_%s_R0p%d.root",etaWidth, radius));


  // we also need to get the files for the MC closure histograms.
  TFile * fMCClosure = TFile::Open("Histogram_pp_PbPb_unfoldMatrix.root");
  
  TFile * fMinBias = TFile::Open(Form("Pawan_ntuple_PbPb_MinBiasData_spectra_JetID_CutA_finebins_CentralityWeightedAllMB_%s_R0p%d.root",etaWidth,radius)); //MinBias File 

  cout<<"after input file declaration"<<endl;
  // need to make sure that the file names are in prefect order so that i can run them one after another. 
  // for the above condition, i might have to play with the date stamp. 
  
  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
  
  // histogram declarations with the following initial appendage: d - Data, m - MC, u- Unfolded
  // for the MC closure test, ive kept separate 
  
  // setup the radius and the eta bin loop here later. not for the time being. Aug 20th. only run the -2 < eta < 2 with the differenent centrality bins 
  
  TH1F *uPbPb_SVD[nbins_cent];
  TH1F *uPbPb_MC_SVD[nbins_cent];
  
  TH1F *dPbPb_TrgComb[nbins_cent+1];
  
  TH1F *mPbPb_Gen[nbins_cent+1], *mPbPb_Reco[nbins_cent+1];
  TH2F *mPbPb_Matrix[nbins_cent+1];
  TH1F *mPbPb_mcclosure_data[nbins_cent+1];
  TH2F *mPbPb_mcclosure_Matrix[nbins_cent+1];
  TH1F *mPbPb_mcclosure_gen[nbins_cent+1];
  
  TH1F *dPP_Comb;
  TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix;
  TH1F *mPP_mcclosure_data;
  TH2F *mPP_mcclosure_Matrix;
  TH1F *mPP_mcclosure_gen;

  TH1F *hMinBias[nbins_cent+1];
  TH1F *uPP_SVD;
  TH1F *uPP_MC_SVD;

  
  TH1F * htest = new TH1F("htest","",nbins_pt, boundaries_pt);
  Int_t unfoldingCutBin = htest->FindBin(unfoldingCut);

  Float_t scaleFactor_dataMC[nbins_cent] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  //Float_t scaleFactor_dataMC[nbins_cent] = {1./145.156/1e5, 0.0, 0.0, 0.0, 0.0, 0.0};
  Float_t scaleFactor_dataMC_pp = 0.0;
  
  // get PbPb data
  for(int i = 0;i<nbins_cent;++i){
    if(printDebug) cout<<"cent_"<<i<<endl;

    //  hMinBias[i]     = (TH1F*)fMinBias->Get(Form("hpbpb_noTrg_R%d_%s_cent%d",radius,etaWidth,i)); //MinBias Histo
    hMinBias[i]     = (TH1F*)fMinBias->Get(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i)); //MinBias Histo
    dPbPb_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i));
    // //dPbPb_TrgComb[i]->Scale(4*145.156*1e6);
    dPbPb_TrgComb[i]->Print("base");

    TH1F * hMinBias_test = new TH1F("hMinBias_test","",501,0,501);
    for(int j = 0; j<hMinBias[i]->GetNbinsX();++j){
      hMinBias_test->SetBinContent(j+1, hMinBias[i]->GetBinContent(j+1));
      hMinBias_test->SetBinError(j+1, hMinBias[i]->GetBinError(j+1));
    }
    hMinBias_test->Print("base");
    hMinBias[i]->Print("base");

    Float_t   bin_no = dPbPb_TrgComb[i]->FindBin(15);
    Float_t bin_end=dPbPb_TrgComb[i]->FindBin(25);
      
    Float_t   bin_nomb = hMinBias_test->FindBin(15);
    Float_t bin_endmb=hMinBias_test->FindBin(25);
      
    float scalerangeweight=dPbPb_TrgComb[i]->Integral(bin_no,bin_end)/hMinBias_test->Integral(bin_nomb,bin_endmb);

    for(int j = 0; j<hMinBias_test->GetNbinsX(); ++j)
      hMinBias_test->SetBinError(j+1,hMinBias_test->GetBinError(j+1)/scalerangeweight);

    hMinBias_test->Scale(scalerangeweight);

    hMinBias[i]->Scale(scalerangeweight);
    dPbPb_TrgComb[i]->Add(hMinBias_test, -1);

    delete hMinBias_test;
    
    dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins_pt, Form("PbPb_data_minbiasSub_cent%d",i), boundaries_pt);
    divideBinWidth(dPbPb_TrgComb[i]);

    //dPbPb_TrgComb[i]->Scale(1./(145.156 * 1e3)); // scale it to milli barns which is the luminosity of the MC 

    mPbPb_Gen[i] = (TH1F*)fPbPb_MC_in->Get(Form("hpbpb_anaBin_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i));
    mPbPb_Gen[i]->Print("base");
    mPbPb_Reco[i] = (TH1F*)fPbPb_MC_in->Get(Form("hpbpb_anaBin_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i));
    mPbPb_Reco[i]->Print("base");
    mPbPb_Matrix[i] = (TH2F*)fPbPb_MC_in->Get(Form("hpbpb_anaBin_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i));
    mPbPb_Matrix[i]->Print("base");

    cout<<"matrix bin content 100 before scaling = "<<mPbPb_Matrix[i]->GetBinContent(22,22)<<endl;
    cout<<"matrix bin Error   100 before scaling = "<<mPbPb_Matrix[i]->GetBinError(22,22)<<endl;
    
    scaleFactor_dataMC[i] = (Float_t) mPbPb_Reco[i]->GetBinContent(mPbPb_Reco[i]->FindBin(120)) / dPbPb_TrgComb[i]->GetBinContent(dPbPb_TrgComb[i]->FindBin(120));
    cout<<"scale factor for cent "<<i<<" = "<<scaleFactor_dataMC[i]<<endl;
    // dPbPb_TrgComb[i]->Scale(scaleFactor_dataMC[i]);
    mPbPb_Gen[i]->Scale(1./scaleFactor_dataMC[i]);
    mPbPb_Reco[i]->Scale(1./scaleFactor_dataMC[i]);
    mPbPb_Matrix[i]->Scale(1./scaleFactor_dataMC[i]);

    cout<<"matrix bin content 100 after scaling = "<<mPbPb_Matrix[i]->GetBinContent(22,22)<<endl;
    cout<<"matrix bin Error   100 after scaling = "<<mPbPb_Matrix[i]->GetBinError(22,22)<<endl;
    
    mPbPb_mcclosure_data[i] = (TH1F*)fMCClosure->Get(Form("akPu%dJetAnalyzer/hrec_c_pbpb_akPu%d_1_%d",radius, radius, i));
    mPbPb_mcclosure_data[i]->Print("base");
    mPbPb_mcclosure_gen[i] = (TH1F*)fMCClosure->Get(Form("akPu%dJetAnalyzer/hgen_pbpb_akPu%d_1_%d",radius, radius, i));
    mPbPb_mcclosure_gen[i]->Print("base");
    mPbPb_mcclosure_Matrix[i] = (TH2F*)fMCClosure->Get(Form("akPu%dJetAnalyzer/hmatrix_pbpb_akPu%d_1_%d",radius, radius, i));
    mPbPb_mcclosure_Matrix[i]->Print("base");
    
    for(int k = 1;k<=unfoldingCutBin;k++){
      dPbPb_TrgComb[i]->SetBinContent(k,0);
      mPbPb_Gen[i]->SetBinContent(k,0);
      mPbPb_Reco[i]->SetBinContent(k,0);
      dPbPb_TrgComb[i]->SetBinError(k,0);
      mPbPb_Gen[i]->SetBinError(k,0);
      mPbPb_Reco[i]->SetBinError(k,0);
      for(int l = 1;l<=nbins_pt;l++){
    	//mPbPb_Matrix[i]->SetBinContent(k,l,0);
    	mPbPb_Matrix[i]->SetBinContent(l,k,0);
    	//mPbPb_Matrix[i]->SetBinError(k,l,0);
    	mPbPb_Matrix[i]->SetBinError(l,k,0);	
      }
    }

    // for(int k = 28;k<=nbins_pt;k++){
    //   dPbPb_TrgComb[i]->SetBinContent(k,0);
    //   mPbPb_Gen[i]->SetBinContent(k,0);
    //   mPbPb_Reco[i]->SetBinContent(k,0);
    //   dPbPb_TrgComb[i]->SetBinError(k,0);
    //   mPbPb_Gen[i]->SetBinError(k,0);
    //   mPbPb_Reco[i]->SetBinError(k,0);
    //   for(int l = 1;l<=nbins_pt;l++){
    // 	mPbPb_Matrix[i]->SetBinContent(k,l,0);
    // 	mPbPb_Matrix[i]->SetBinContent(l,k,0);
    // 	mPbPb_Matrix[i]->SetBinError(k,l,0);
    // 	mPbPb_Matrix[i]->SetBinError(l,k,0);	
    //   }
    // }
  }

  dPP_Comb = (TH1F*)fPP_in->Get(Form("hpp_anaBin_HLTComb_R%d_%s",radius,etaWidth));
  // dPP_Comb->Scale(1./(5.3 * 1e3)); // scale it to luminosity of the MC but there is something going wrong here.
  // get PP MC
  mPP_Gen = (TH1F*)fPP_in->Get(Form("hpp_anaBin_JetComb_gen_R%d_%s",radius,etaWidth));
  mPP_Gen->Print("base");
  mPP_Reco = (TH1F*)fPP_in->Get(Form("hpp_anaBin_JetComb_reco_R%d_%s",radius,etaWidth));
  mPP_Reco->Print("base");
  mPP_Matrix = (TH2F*)fPP_in->Get(Form("hpp_anaBin_matrix_HLT_R%d_%s",radius,etaWidth));
  mPP_Matrix->Print("base");

  scaleFactor_dataMC_pp = (Float_t)mPP_Reco->GetBinContent(mPP_Reco->FindBin(120))/dPP_Comb->GetBinContent(dPP_Comb->FindBin(120));
  // dPP_Comb->Scale(scaleFactor_dataMC_pp);
  mPP_Gen->Scale(1./scaleFactor_dataMC_pp);
  mPP_Reco->Scale(1./scaleFactor_dataMC_pp);
  mPP_Matrix->Scale(1./scaleFactor_dataMC_pp);
  
  // get the mc closure test histograms
  mPP_mcclosure_data = (TH1F*)fMCClosure->Get(Form("ak%dJetAnalyzer/hrec_c_pp_ak%d_0_7",radius,radius));
  mPP_mcclosure_data->Print("base");
  mPP_mcclosure_gen = (TH1F*)fMCClosure->Get(Form("ak%dJetAnalyzer/hgen_pp_ak%d_0_7",radius,radius));
  mPP_mcclosure_gen->Print("base");
  mPP_mcclosure_Matrix = (TH2F*)fMCClosure->Get(Form("ak%dJetAnalyzer/hmatrix_pp_ak%d_0_7",radius,radius));
  mPP_mcclosure_Matrix->Print("base");

  
  // would be better to read in the histograms and rebin them. come to think of it, it would be better to have them already rebinned (and properly scaled - to the level of differential cross section in what ever barns (inverse micro barns) but keep it consistent) from the read macro. 

  
  Int_t nSVDIter = 1.0;
  
  for(int j = 3; j<=7; ++j){
    nSVDIter = j;
    
    for(int i = 0;i<nbins_cent;++i){
      // rebin histograms before we send it to unfolding 
      // histograms are already binned
      // need to scale the Data to be in the same units as the prior. MC is in mb and data was in counts. 
      cout<<"cent = "<<i<<endl;
      //since SVD is very straight forward, lets do it rignt here:
      //get the SVD response matrix:
      RooUnfoldResponse ruResponse(mPbPb_Reco[i],mPbPb_Gen[i], mPbPb_Matrix[i],"","");
      //regularization parameter definition: 
      RooUnfoldSvd unfoldSvd(&ruResponse, dPbPb_TrgComb[i], nSVDIter);
      uPbPb_SVD[i] = (TH1F*)unfoldSvd.Hreco();
      uPbPb_SVD[i]->SetName(Form("PbPb_SVD_unfolding_cent%d",i));

      // RooUnfoldResponse ruResponseMC(mPbPb_mcclosure_Matrix[i]->ProjectionY(),mPbPb_mcclosure_Matrix[i]->ProjectionX(), mPbPb_mcclosure_Matrix[i],"","");
      // // //regularization parameter definition: 
      // RooUnfoldSvd unfoldSvdMC(&ruResponseMC, mPbPb_mcclosure_data[i], nSVDIter);
      // uPbPb_MC_SVD[i] = (TH1F*)unfoldSvdMC.Hreco();
      // uPbPb_MC_SVD[i]->SetName(Form("PbPb_MC_SVD_unfolding_cent%d",i));
      
    }// centrality bin loop

    // RooUnfoldResponse ruResponseMCPP(mPP_Reco, mPP_Gen, mPP_mcclosure_Matrix,"","");
    // //regularization parameter definition: 
    // RooUnfoldSvd unfoldSvdMCPP(&ruResponseMCPP, mPP_mcclosure_data, nSVDIter);
    // uPP_MC_SVD = (TH1F*)unfoldSvdMCPP.Hreco();
    // uPP_MC_SVD->SetName("PP_MC_SVD_unfolding");
    
    RooUnfoldResponse ruResponsePP(mPP_Reco, mPP_Gen, mPP_Matrix,"","");
    // //regularization parameter definition: 
    RooUnfoldSvd unfoldSvdPP(&ruResponsePP, dPP_Comb, nSVDIter);
    uPP_SVD = (TH1F*)unfoldSvdPP.Hreco();
    uPP_SVD->SetName("PP_SVD_unfolding");
  
    TFile fout(Form("svd_unfolding_matrix_param%d_%s_%dGeVCut_R%d_%d.root",j,etaWidth,unfoldingCut,radius,date.GetDate()),"RECREATE");

    for(int i = 0; i<nbins_cent;i++){

      //dPbPb_TrgComb[i]->Scale(1./scaleFactor_dataMC[i]);
      //uPbPb_SVD[i]->Scale(1./scaleFactor_dataMC[i]);
      dPbPb_TrgComb[i]->Write();
      uPbPb_SVD[i]->Write();
      mPbPb_Matrix[i]->Write();
      mPbPb_Gen[i]->Write();
      // mPbPb_mcclosure_data[i]->Write();
      // uPbPb_MC_SVD[i]->Write();
      // mPbPb_mcclosure_Matrix[i]->Write();
      // mPbPb_mcclosure_gen[i]->Write();
    }

    //uPP_SVD->Scale(1./scaleFactor_dataMC_pp);
    //dPP_Comb->Scale(1./scaleFactor_dataMC_pp);
    uPP_SVD->Write();
    dPP_Comb->Write();
    mPP_Matrix->Write();
    mPP_Gen->Write();
    // mPP_mcclosure_data->Write();
    // uPP_MC_SVD->Write();
    // mPP_mcclosure_Matrix->Write();
    // mPP_mcclosure_gen->Write();
    
    fout.Close();
    
  }// svd unfolding iterations loop
  
  timer.Stop();
  if(printDebug) cout<<"CPU time (mins) = "<<(Float_t)timer.CpuTime()/60<<endl;
  if(printDebug) cout<<"Real tile (mins) = "<<(Float_t)timer.RealTime()/60<<endl;
  

}
