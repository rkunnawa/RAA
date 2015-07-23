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

#include "../Headers/plot.h"
#include "../Headers/utilities.h"

void RAA_dataDrivenUnfoldingErrorCheck_new(int radius = 3, bool isATLASCut = true)
{

  TStopwatch timer; 
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  bool printDebug = true;
  bool dofakeremove = false;
  bool do10GeVBins = true;

  char * scale = (char*)"NeqScale";
  // 944Scale
  // NeqScale
  // NeqScalePerCent
  
  int unfoldingCut = 20;
  char * outLocation = (char*) "July20/";
  if(isATLASCut) outLocation = (char*)"July20/ATLASCut/";
  
  // get the data and mc histograms from the output of the read macro. 
  
  TDatime date;//this is just here to get them to run optimized. 

  // Pawan's files:
  //TFile * fPbPb_in = TFile::Open(Form("Pawan_TTree_PbPb_Data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root", etaWidth, radius));
  //TFile * fPbPb_in = TFile::Open(Form("Pawan_TTree_PbPb_Data_MC_subid0_spectra_JetID_CutA_muMaxOverSumcandMaxLT0p975_finebins_%s_R0p%d.root", etaWidth, radius));
  //TFile * fPbPb_in = TFile::Open(Form("Pawan_TTree_PbPb_Data_MC_subid0_spectra_JetID_CutA_trkMaxOverpfptGT0p02_finebins_%s_R0p%d.root", etaWidth, radius));
  //TFile * fPbPb_in = TFile::Open(Form("Pawan_TTree_PbPb_Data_noPrescl_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root", etaWidth, radius));
  //TFile * fPbPb_in = TFile::Open(Form("Pawan_TTree_PbPb_Data_noPrescl_MC_subid0_spectra_JetID_CutA_7GeVTrackCut_finebinscut_%s_R0p%d.root", etaWidth, radius));
  // get the files to perform MC closure from the fixed ntuples
  TFile * fPbPb_in, *fMinBias, *fPbPb_MC_in, * fTrig;
  
  if(isATLASCut) fPbPb_in    = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_Data_histograms_FromForest_trkMax7OrNeMax8GeVCut_akPu%d_20_eta_20.root",radius),"r");
  if(!isATLASCut) fPbPb_in    = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_Data_histograms_FromForest_akPu%d_20_eta_20.root",radius),"r");
  if(isATLASCut) fMinBias    = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_MB_Data_histograms_FromForest_trkMax7OrNeMax8GeVCut_fix_pt15GeVCut_akPu%d_20_eta_20.root",radius));
  if(!isATLASCut) fMinBias    = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_MB_Data_histograms_FromForest_fix_pt15GeVCut_akPu%d_20_eta_20.root",radius));
  //TFile * fPbPb_MC_in = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_MC_histograms_FromForest_pthat50andabove_akPu%d_20_eta_20.root",radius));
  if(isATLASCut) fPbPb_MC_in = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_MC_histograms_FromForest_trkMax7OrNeMax8GeVCut_akPu%d_20_eta_20.root",radius),"r");
  if(!isATLASCut) fPbPb_MC_in = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_MC_histograms_FromForest_akPu%d_20_eta_20.root",radius),"r");
  TFile * fPP_in      = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/pp_Data_histograms_FromForest_ak%d_20_eta_20.root",radius),"r");
  TFile * fPP_MC_in   = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/pp_MC_histograms_FromForest_ak%d_20_eta_20.root",radius),"r");
  if(isATLASCut) fTrig       = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_Data_nofkSub_trkMax7OrNeMax8GeVCut_new_MC_turnonCurves_R%d_20_eta_20_20150715.root",radius),"r");
  if(!isATLASCut) fTrig       = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/PbPb_Data_nofkSub_new_MC_turnonCurves_R%d_20_eta_20_20150715.root",radius),"r");
  TFile * fTrig_pp = new TFile(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/Jun29/pp_MC_turnonCurves_R%d_20_eta_20_20150702.root", radius),"r");
  // TFile * fPPMCTrig = new TFile(Form("PP_Data_MC_turnonCurves_R%d_20_eta_20_20150618.root",radius),"r");

  TH1F * hMC_turnon[nbins_cent+1], * hData_turnon[nbins_cent+1];
  cout<<"after input file declaration"<<endl;
  
  
  // histogram declarations with the following initial appendage: d - Data, m - MC, u- Unfolded
  // for the MC closure test, ive kept separate 

  TH1F *dPbPb_TrgComb[nbins_cent+1], *dPbPb_TrgCombInput[nbins_cent+1], *dPbPb_Comb[nbins_cent+1], *dPbPb_Trg80[nbins_cent+1], *dPbPb_Trg65[nbins_cent+1], *dPbPb_Trg55[nbins_cent+1], *dPbPb_1[nbins_cent+1], *dPbPb_2[nbins_cent+1], *dPbPb_3[nbins_cent+1], *dPbPb_80[nbins_cent+1], *dPbPb_65[nbins_cent+1], *dPbPb_55[nbins_cent+1];
  
  TH1F *mPbPb_GenInput[nbins_cent+1], *mPbPb_RecoInput[nbins_cent+1];
  TH1F *mPbPb_Gen[nbins_cent+1], *mPbPb_Reco[nbins_cent+1];
  TH2F *mPbPb_Matrix[nbins_cent+1], * mPbPb_MatrixInput[nbins_cent+1], *mPbPb_Response[nbins_cent+1], *mPbPb_ResponseNorm[nbins_cent+1];
  TH1F *mPbPb_mcclosure_data[nbins_cent+1];
  TH2F *mPbPb_mcclosure_Matrix[nbins_cent+1],*mPbPb_mcclosure_Response[nbins_cent+1], *mPbPb_mcclosure_ResponseNorm[nbins_cent+1];
  TH1F *mPbPb_mcclosure_gen[nbins_cent+1];
  const int Iterations = 20; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1], *uPbPb_BinByBin[nbins_cent+1], *uPbPb_SVD[nbins_cent+1]; 
  TH1F *uPbPb_BayesianIter[nbins_cent+1][Iterations];
  TH1F *dPbPb_MinBias[nbins_cent];
  TH1F *hMinBias[nbins_cent];
  
  TH1F *dPP_1, *dPP_2, *dPP_3, *dPP_Comb, * dPP_CombInput;
  TH1F *mPP_Gen, *mPP_Reco, *mPP_GenInput, *mPP_RecoInput;
  TH2F *mPP_Matrix, *mPP_MatrixInput, *mPP_Response,*mPP_ResponseNorm;
  TH1F *mPP_mcclosure_data;
  TH2F *mPP_mcclosure_Matrix, *mPP_mcclosure_Response,*mPP_mcclosure_ResponseNorm;
  TH1F *mPP_mcclosure_Gen;
  TH1F *uPP_Bayes, *uPP_BinByBin, *uPP_SVD;
  TH1F *uPP_BayesianIter[Iterations];

  TH1F * hData_FaketoSub_fullbin[nbins_cent+1];
  
  // would be better to read in the histograms and rebin them. come to think of it, it would be better to have them already rebinned (and properly scaled - to the level of differential cross section in what ever barns (inverse micro barns) but keep it consistent) from the read macro. 

  if(radius == 2) unfoldingCut = unfoldingCut_R2;
  if(radius == 3) unfoldingCut = unfoldingCut_R3;
  if(radius == 4) unfoldingCut = unfoldingCut_R4;
  
  // TH1F * htest = new TH1F("htest","",nbins, ptbins_long);
  // Int_t unfoldingCutBin = htest->FindBin(unfoldingCut);
  
  //Float_t cutarray[6]={50,50,40,35,35,35};
  float cutarray[nbins_cent] = {0.0,0.0,0.0,0.0,0.0,0.0};
  
  if(radius == 2){
    cutarray[0] = 50;
    cutarray[1] = 40;
    cutarray[2] = 40;
    cutarray[3] = 40;
    cutarray[4] = 30;
    cutarray[5] = 30;
  }
  if(radius == 3){
    cutarray[0] = 55;
    cutarray[1] = 50;
    cutarray[2] = 50;
    cutarray[3] = 40;
    cutarray[4] = 35;
    cutarray[5] = 40;
  }
    
  if(radius == 4){
    cutarray[0] = 70;
    cutarray[1] = 60;
    cutarray[2] = 60;
    cutarray[3] = 45;
    cutarray[4] = 40;
    cutarray[5] = 30;
  }
    
  TH1F * hDataBeforeSub[nbins_cent], * hDataAfterSub[nbins_cent];
  
  // get PbPb data
  for(int i = 0;i<nbins_cent;++i){
    if(printDebug) cout<<"cent_"<<i<<endl;

    hData_turnon[i] = (TH1F*)fTrig->Get(Form("hHist_Data_Turnon_cent%d",i));
    
    hData_turnon[i] = (TH1F*)hData_turnon[i]->Rebin(nbins_short, Form("hData_turnon_cent%d",i), boundaries_short);
    //hData_turnon[i] = (TH1F*)hData_turnon[i]->Rebin(10);

    divideBinWidth(hData_turnon[i]);

     //  hMinBias[i]     = (TH1F*)fMinBias->Get(Form("hpbpb_noTrg_R%d_%s_cent%d",radius,etaWidth,i)); //MinBias Histo
    hMinBias[i]     = (TH1F*)fMinBias->Get(Form("hpbpb_HLTMBwoLJSbJ_R%d_%s_cent%d",radius,etaWidth,i)); //MinBias Histo
    hMinBias[i]->Print("base");

    dPbPb_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLTComb_%s_R%d_%s_cent%d",scale, radius,etaWidth,i));
    //dPbPb_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i));
    //dPbPb_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i));
    //dPbPb_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i));
    // dPbPb_TrgComb[i]->Scale(1./(1+0.898+0.494)/1e16);
    // //dPbPb_TrgComb[i]->Scale(4*145.156*1e6);
    dPbPb_TrgComb[i]->Print("base");

    hDataBeforeSub[i] = (TH1F*)dPbPb_TrgComb[i]->Clone(Form("hData_Before_Sub_cent%d",i));
    
    //dPbPb_TrgComb[i] = (TH1F*)fPbPb_DataCorr_in->Get(Form("Data_TrigEffCorrected_FakeSub_cent%d",i));
    //dPbPb_TrgComb[i]->Print("base");
    
    // dPbPb_JEC_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_JEC_HLTComb_R%d_%s_cent%d",radius,etaWidth,i));
    // // //dPbPb_TrgComb[i]->Scale(4*145.156*1e6);
    // dPbPb_JEC_TrgComb[i]->Print("base");
    // dPbPb_Smear_TrgComb[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_Smear_HLTComb_R%d_%s_cent%d",radius,etaWidth,i));
    // // //dPbPb_TrgComb[i]->Scale(4*145.156*1e6);
    // dPbPb_Smear_TrgComb[i]->Print("base");
    // dPbPb_Trg80[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i));
    // //dPbPb_Trg80[i]->Scale(4*145.156*1e6);
    // dPbPb_Trg80[i]->Print("base");
    // dPbPb_Trg65[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i));
    // //dPbPb_Trg65[i]->Scale(4*145.156*1e6);
    // dPbPb_Trg65[i]->Print("base");
    // dPbPb_Trg55[i] = (TH1F*)fPbPb_in->Get(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i));
    // //dPbPb_Trg55[i]->Scale(4*145.156*1e6);
    // dPbPb_Trg55[i]->Print("base");

    // if(dotrigcor){
      
    //   doTrigCorr(dPbPb_TrgComb[i], hData_turnon[i]);
    //   doTrigCorr(dPbPb_JEC_TrgComb[i], hData_turnon[i]);
    //   doTrigCorr(dPbPb_Smear_TrgComb[i], hData_turnon[i]);

    // }
    
    //Lets do the subtraction here _Sevil 

      // Float_t bincon=cutarray[i]; 
      // Int_t bincut= hMinBias[i]->FindBin(bincon); 
      
      // for(int k = bincut;k<=hMinBias[i]->GetNbinsX();k++) { 
      // 	hMinBias[i]->SetBinContent(k,0); 
      // 	hMinBias[i]->SetBinError(k,0);
      // } 

      // for(int k = 1;k<=15;k++) { 
      // 	hMinBias[i]->SetBinContent(k,0); 
      // 	hMinBias[i]->SetBinError(k,0);
      // }
      
      Float_t   bin_no = dPbPb_TrgComb[i]->FindBin(15);
      Float_t bin_end=dPbPb_TrgComb[i]->FindBin(25);
      
      Float_t   bin_nomb = hMinBias[i]->FindBin(15);
      Float_t bin_endmb=hMinBias[i]->FindBin(25);
      
      float scalerangeweight=dPbPb_TrgComb[i]->Integral(bin_no,bin_end)/hMinBias[i]->Integral(bin_nomb,bin_endmb);

      // for(int j = 0; j<hMinBias[i]->GetNbinsX(); ++j)
      // 	hMinBias[i]->SetBinError(j+1, (Float_t)hMinBias[i]->GetBinError(j+1)/scalerangeweight);
      
      hMinBias[i]->Scale(scalerangeweight);
      if(dofakeremove) dPbPb_TrgComb[i]->Add(hMinBias[i], -1);

      hDataAfterSub[i] = (TH1F*)dPbPb_TrgComb[i]->Clone(Form("hData_After_Sub_cent%d",i));
      
      // dPbPb_JEC_TrgComb[i]->Add(hMinBias[i], -1);
      // dPbPb_Smear_TrgComb[i]->Add(hMinBias[i], -1);

      // for(int j = 1; j<dPbPb_TrgComb[i]->GetNbinsX(); ++j){

      // 	if(dPbPb_TrgComb[i]->GetBinContent(j) <= 0 ||dPbPb_JEC_TrgComb[i]->GetBinContent(j) <= 0||dPbPb_Smear_TrgComb[i]->GetBinContent(j) <= 0){
      // 	  dPbPb_TrgComb[i]->SetBinContent(j, 0);
      // 	  dPbPb_JEC_TrgComb[i]->SetBinContent(j, 0);
      // 	  dPbPb_Smear_TrgComb[i]->SetBinContent(j, 0);
      // 	  dPbPb_TrgComb[i]->SetBinError(j, 0);
      // 	  dPbPb_JEC_TrgComb[i]->SetBinError(j, 0);
      // 	  dPbPb_Smear_TrgComb[i]->SetBinError(j, 0);
      // 	}
      // }
      
    

    // // lets truncate the histograms here:
    // cout<<" going to truncate Data histogram here cent "<<i<<endl;
    // dPbPb_TrgCombInput[i]->Print("base");

    // dPbPb_TrgComb[i] = new TH1F(Form("PbPb_data_minbiasSub_cent%d",i),"",365, 30, 395);
    // Truncate1D(dPbPb_TrgCombInput[i], dPbPb_TrgComb[i]);    
    
    // // dPbPb_TrgComb[i] = (TH1F*)Truncate1D(dPbPb_TrgComb[i], 340, unfoldingCutBin, 395);
    // // dPbPb_TrgComb[i]->Print("base");
    
    //dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins, Form("PbPb_data_minbiasSub_cent%d",i), ptbins_long);
    dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(10);
    dPbPb_TrgComb[i]->SetName(Form("PbPb_data_minbiasSub_cent%d",i));
    divideBinWidth(dPbPb_TrgComb[i]);
    
    hMinBias[i] = (TH1F*)hMinBias[i]->Rebin(10);
    hDataAfterSub[i] = (TH1F*)hDataAfterSub[i]->Rebin(10);
    hDataBeforeSub[i] = (TH1F*)hDataBeforeSub[i]->Rebin(10);

    divideBinWidth(hMinBias[i]);
    divideBinWidth(hDataAfterSub[i]);
    divideBinWidth(hDataBeforeSub[i]);
    
    dPbPb_TrgComb[i]->Scale(1./(166 * 1e9));
    
    // dPbPb_TrgComb[i]->Print("base");
    
  }

  if(printDebug)cout<<"loaded the data histograms PbPb"<<endl;
  // get PbPb MC
  for(int i = 0;i<nbins_cent;i++){
    
    // mPbPb_GenInput[i] = (TH1F*)fPbPb_MC_in->Get(Form("hpbpb_anaBin_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i));
    // mPbPb_GenInput[i]->Print("base");
    // mPbPb_RecoInput[i] = (TH1F*)fPbPb_MC_in->Get(Form("hpbpb_anaBin_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i));
    // mPbPb_RecoInput[i]->Print("base");
    // mPbPb_MatrixInput[i] = (TH2F*)fPbPb_MC_in->Get(Form("hpbpb_anaBin_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i));
    // mPbPb_MatrixInput[i]->Print("base");

    mPbPb_Gen[i] = (TH1F*)fPbPb_MC_in->Get(Form("hpbpb_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i));
    //mPbPb_Gen[i]->Rebin(nbins, Form("mPbPb_Gen_cent%d",i), ptbins_long);
    mPbPb_Gen[i]->Rebin(10);
    divideBinWidth(mPbPb_Gen[i]);
    mPbPb_Gen[i]->Print("base");
    mPbPb_Reco[i] = (TH1F*)fPbPb_MC_in->Get(Form("hpbpb_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i));
    //mPbPb_Reco[i]->Rebin(nbins, Form("mPbPb_Reco_cent%d",i), ptbins_long);
    mPbPb_Reco[i]->Rebin(10);
    divideBinWidth(mPbPb_Reco[i]);
    mPbPb_Reco[i]->Print("base");
    mPbPb_Matrix[i] = (TH2F*)fPbPb_MC_in->Get(Form("hpbpb_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i));
    //mPbPb_Matrix[i] = (TH2F*)fPbPb_MC_in->Get(Form("hpbpb_anaBin_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i));
    mPbPb_Matrix[i]->Rebin2D(10, 10);
    mPbPb_Matrix[i]->Print("base");
    
    // if(etaWidth == "10_eta_10"){
    //   if(i == 0 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 1 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 2 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 3 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 4 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==2) unfoldingCutBin = htest->FindBin(30);

    //   if(i == 0 && radius==3) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 1 && radius==3) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 2 && radius==3) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 3 && radius==3) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 4 && radius==3) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==3) unfoldingCutBin = htest->FindBin(30);

    //   if(i == 0 && radius==4) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 1 && radius==4) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 2 && radius==4) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 3 && radius==4) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 4 && radius==4) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==4) unfoldingCutBin = htest->FindBin(30);
    // }

    // if(etaWidth == "10_eta_18"){
    //   if(i == 0 && radius==2) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 1 && radius==2) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 2 && radius==2) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 3 && radius==2) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 4 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==2) unfoldingCutBin = htest->FindBin(30);

    //   if(i == 0 && radius==3) unfoldingCutBin = htest->FindBin(60);
    //   if(i == 1 && radius==3) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 2 && radius==3) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 3 && radius==3) unfoldingCutBin = htest->FindBin(40);
    //   if(i == 4 && radius==3) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==3) unfoldingCutBin = htest->FindBin(30);

    //   if(i == 0 && radius==4) unfoldingCutBin = htest->FindBin(70);
    //   if(i == 1 && radius==4) unfoldingCutBin = htest->FindBin(60);
    //   if(i == 2 && radius==4) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 3 && radius==4) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 4 && radius==4) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==4) unfoldingCutBin = htest->FindBin(30);
    // }

    // if(etaWidth == "20_eta_20"){
    //   if(i == 0 && radius==2) unfoldingCutBin = htest->FindBin(70);
    //   if(i == 1 && radius==2) unfoldingCutBin = htest->FindBin(60);
    //   if(i == 2 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 3 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 4 && radius==2) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==2) unfoldingCutBin = htest->FindBin(30);

    //   if(i == 0 && radius==3) unfoldingCutBin = htest->FindBin(70);
    //   if(i == 1 && radius==3) unfoldingCutBin = htest->FindBin(60);
    //   if(i == 2 && radius==3) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 3 && radius==3) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 4 && radius==3) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==3) unfoldingCutBin = htest->FindBin(30);

    //   if(i == 0 && radius==4) unfoldingCutBin = htest->FindBin(80);
    //   if(i == 1 && radius==4) unfoldingCutBin = htest->FindBin(60);
    //   if(i == 2 && radius==4) unfoldingCutBin = htest->FindBin(50);
    //   if(i == 3 && radius==4) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 4 && radius==4) unfoldingCutBin = htest->FindBin(30);
    //   if(i == 5 && radius==4) unfoldingCutBin = htest->FindBin(30);
    // }

    int bincut = mPbPb_Gen[i]->FindBin(50);
    for(int k = 1;k<=bincut;k++){

    //   mPbPb_Gen[i]->SetBinContent(k,0);
    //   mPbPb_Reco[i]->SetBinContent(k,0);
    //   mPbPb_Gen[i]->SetBinError(k,0);
    //   mPbPb_Reco[i]->SetBinError(k,0);
    //   // set bin content matrix l,k works 
      for(int l = 1;l<=mPbPb_Gen[i]->GetNbinsX();l++){
	mPbPb_Matrix[i]->SetBinContent(l,k,0);
	mPbPb_Matrix[i]->SetBinError(l,k,0);
      }
      
    }

    SetUnfoldBins1D(dPbPb_TrgComb[i], 50, 350);
    
    // cout<<"going to truncate the MC histograms here."<<endl;
    
    // // mPbPb_Reco[i]->Print("base");
    // // mPbPb_Reco[i] = (TH1F*)Truncate1D_anaBin(mPbPb_Reco[i], nbins, ptbins_long);
    // // mPbPb_Reco[i]->Print("base");

    // mPbPb_GenInput[i]->Print("base");
    // mPbPb_Gen[i] = new TH1F(Form("mPbPb_Gen_spectra_cent%d",i),"",nbins, ptbins_long);
    // Truncate1D(mPbPb_GenInput[i], mPbPb_Gen[i]);
    // mPbPb_Gen[i]->Print("base");
    
    // mPbPb_RecoInput[i]->Print("base");
    // mPbPb_Reco[i] = new TH1F(Form("mPbPb_REco_spectra_cent%d",i),"",nbins, ptbins_long);
    // Truncate1D(mPbPb_RecoInput[i], mPbPb_Reco[i]);
    // mPbPb_Reco[i]->Print("base");
    
    // mPbPb_MatrixInput[i]->Print("base");
    // mPbPb_Matrix[i] = new TH2F(Form("mPbPb_Response_Matrix_cent%d",i),"",nbins_truncated, ptbins_long_truncated, nbins_truncated, ptbins_long_truncated);
    // Truncate2D(mPbPb_MatrixInput[i], mPbPb_Matrix[i]);
    // mPbPb_Matrix[i]->Print("base");
 
  }
  
  if(printDebug) cout<<"loaded the data and mc PbPb histograms from the files"<<endl;

  // get PP data
  if(printDebug) cout<<"Getting PP data and MC"<<endl;

  //fPP_in->ls();

  // dPP_1 = (TH1F*)fPP_in->Get(Form("hpp_HLT80_R%d_%s",radius,etaWidth)); 
  // dPP_1->Print("base");
  // dPP_2 = (TH1F*)fPP_in->Get(Form("hpp_HLT60_R%d_%s",radius,etaWidth));
  // dPP_2->Print("base");
  // dPP_3 = (TH1F*)fPP_in->Get(Form("hpp_HLT40_R%d_%s",radius,etaWidth));
  // dPP_3->Print("base");
  dPP_Comb = (TH1F*)fPP_in->Get(Form("hpp_HLTComb_R%d_%s",radius,etaWidth));
  //dPP_Comb = (TH1F*)dPP_1->Clone(Form("hpp_TrgComb_R%d_n20_eta_p20",radius,etaWidth));   
  //dPP_CombInput->Print("base");
  dPP_Comb->Scale(1./(5.3 * 1e9));
  
  // dPP_Comb = new TH1F("PP_MeasuredSpectra","",365, 30, 395);
  // Truncate1D(dPP_CombInput, dPP_Comb);

  //dPP_Comb = (TH1F*)dPP_Comb->Rebin(nbins, "PP_MeasuredSpectra", ptbins_long);
  dPP_Comb = (TH1F*)dPP_Comb->Rebin(10);
  dPP_Comb->SetName("PP_MeasuredSpectra");
  divideBinWidth(dPP_Comb);
  dPP_Comb->Print("base");
  
  // get PP MC
  // mPP_GenInput = (TH1F*)fPP_MC_in->Get(Form("hpp_anaBin_JetComb_gen_R%d_%s",radius,etaWidth));
  // mPP_GenInput->Print("base");
  // mPP_Gen = new TH1F("mPP_Gen_spectra","",nbins, ptbins_long);
  // Truncate1D(mPP_GenInput, mPP_Gen);
  // mPP_Gen->Print("base");
  
  // mPP_RecoInput = (TH1F*)fPP_MC_in->Get(Form("hpp_anaBin_JetComb_reco_R%d_%s",radius,etaWidth));
  // mPP_RecoInput->Print("base");
  // mPP_Reco = new TH1F("mPP_Reco_spectra","",nbins, ptbins_long);
  // Truncate1D(mPP_RecoInput, mPP_Reco);
  // mPP_Reco->Print("base");
  
  // mPP_MatrixInput = (TH2F*)fPP_MC_in->Get(Form("hpp_anaBin_matrix_HLT_R%d_%s",radius,etaWidth));
  // mPP_MatrixInput->Print("base");
  // mPP_Matrix = new TH2F("mPP_response_Matrix","",nbins_truncated, ptbins_long_truncated, nbins_truncated, ptbins_long_truncated);
  // Truncate2D(mPP_MatrixInput, mPP_Matrix);
  // mPP_Matrix->Print("base");

  // get PP MC
  // change from fPP_MC_in to fPP_in to run finebinscut
  //mPP_Gen = (TH1F*)fPP_MC_in->Get(Form("hpp_anaBin_JetComb_gen_R%d_20_eta_20",radius));
  mPP_Gen = (TH1F*)fPP_MC_in->Get(Form("hpp_JetComb_gen_R%d_20_eta_20",radius));
  //mPP_Gen->Rebin(nbins, "mPP_Gen", ptbins_long);
  mPP_Gen->Rebin(10);
  divideBinWidth(mPP_Gen);
  mPP_Gen->Print("base");

  //mPP_Reco = (TH1F*)fPP_MC_in->Get(Form("hpp_anaBin_JetComb_reco_R%d_20_eta_20",radius));
  mPP_Reco = (TH1F*)fPP_MC_in->Get(Form("hpp_JetComb_reco_R%d_20_eta_20",radius));
  //mPP_Gen->Rebin(nbins, "mPP_Gen", ptbins_long);
  mPP_Reco->Rebin(10);
  divideBinWidth(mPP_Reco);
  mPP_Reco->Print("base");

  //mPP_Matrix = (TH2F*)fPP_MC_in->Get(Form("hpp_anaBin_matrix_HLT_R%d_20_eta_20",radius));
  mPP_Matrix = (TH2F*)fPP_MC_in->Get(Form("hpp_matrix_HLT_R%d_20_eta_20",radius));
  mPP_Matrix->Rebin2D(10, 10);
  mPP_Matrix->Print("base");

  if(printDebug) cout<<"Filling the PbPb response Matrix"<<endl;

  // response matrix and unfolding for PbPb 
  // going to try it the way kurt has its. 

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
  //   dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins,Form("PbPb_measured_spectra_combined_cent%d",i),ptbins_long);
  //   divideBinWidth(dPbPb_TrgComb[i]);
  // }

  // dPP_Comb = (TH1F*)dPP_Comb->Rebin(nbins,"pp_measured_spectra_combined",ptbins_long);
  // divideBinWidth(dPP_Comb);
  // dPP_Comb->Scale(1./ dPP_Comb->GetBinContent(nbins));
  
  // Now that we have all the response matrix for the 6 centralities in PbPb and one pp spectra lets start doing the steps:
  // we have 39 pt bins, so we need 1000 gaussian functions for each pt bin.
  
  Int_t unfoldingTrials = 1000;
  Double_t meanMeasPbPb[nbins][nbins_cent], sigmaMeasPbPb[nbins][nbins_cent];
  Double_t meanMeasPP[nbins], sigmaMeasPP[nbins];
  Double_t meanUnfoldPbPb[nbins][nbins_cent][unfoldingTrials], sigmaUnfoldPbPb[nbins][nbins_cent][unfoldingTrials];
  Double_t meanUnfoldPP[nbins][unfoldingTrials], sigmaUnfoldPP[nbins][unfoldingTrials]; 
  
  TRandom3 *random = new TRandom3(0);

  TH1F * hPbPb_beforeUnfold_Gaussian_pt150[nbins_cent];
  TH1F * hPP_beforeUnfold_Gaussian_pt150; 
  hPP_beforeUnfold_Gaussian_pt150 = new TH1F("hPP_beforeUnfold_Gaussian_pt150","",1000, 0.1 * dPP_Comb->GetBinContent(dPP_Comb->FindBin(150)) , 1.9 * dPP_Comb->GetBinContent(dPP_Comb->FindBin(150)));
  
  for(int i = 0; i<nbins_cent; ++i)
    hPbPb_beforeUnfold_Gaussian_pt150[i] = new TH1F(Form("hPbPb_beforeUnfold_Gaussian_pt150_cent%d",i),"Before Unfolding pt bin at 150 value spectra",1000, 0.1 * dPbPb_TrgComb[i]->GetBinContent(dPbPb_TrgComb[i]->FindBin(150)), 1.9 * dPbPb_TrgComb[i]->GetBinContent(dPbPb_TrgComb[i]->FindBin(150)));
  

  for(int u = 0;u<unfoldingTrials;++u){
    cout<<"unfolding trial no = "<<u+1<<endl;
  
    for(int j = 0;j<nbins;++j){
      for(int i = 0;i<nbins_cent;++i){
      
	meanMeasPbPb[j][i] = dPbPb_TrgComb[i]->GetBinContent(j+1);
	sigmaMeasPbPb[j][i] = dPbPb_TrgComb[i]->GetBinError(j+1);

      }// centrality loop

      meanMeasPP[j] = dPP_Comb->GetBinContent(j+1);
      sigmaMeasPP[j] = dPP_Comb->GetBinError(j+1);
      
    }// nbins loop

    // now proceed to unfolding for each trial.

    for(int i = 0;i<nbins_cent;++i){

      TH1F * hPreUnfoldingSpectra = new TH1F("hPreUnfoldingSpectra","",nbins,0, 1000);
      TH1F * hAfterUnfoldingSpectra;

      for(int j = 0;j<nbins;++j){
	
	hPreUnfoldingSpectra->SetBinContent(j+1, random->Gaus(meanMeasPbPb[j][i], sigmaMeasPbPb[j][i]));
	hPreUnfoldingSpectra->SetBinError(j+1, sigmaMeasPbPb[j][i]/sqrt(unfoldingTrials));
	if(j+1 == dPbPb_TrgComb[i]->FindBin(150)) hPbPb_beforeUnfold_Gaussian_pt150[i]->Fill(random->Gaus(meanMeasPbPb[j][i], sigmaMeasPbPb[j][i]));
	
      }// nbins loop

      TH1F* hMCGen          = (TH1F*)mPbPb_Response[i]->ProjectionX();
      removeZero(hMCGen);
      bayesianUnfold myUnfoldingMulti(mPbPb_Matrix[i], hMCGen, 0);
      myUnfoldingMulti.unfold(hPreUnfoldingSpectra, BayesIter);

      hAfterUnfoldingSpectra = (TH1F*) myUnfoldingMulti.hPrior->Clone("hAfterUnfoldingSpectra");

      for(int j = 0;j<nbins;++j){

	meanUnfoldPbPb[j][i][u] = hAfterUnfoldingSpectra->GetBinContent(j+1);
	sigmaUnfoldPbPb[j][i][u] = hAfterUnfoldingSpectra->GetBinError(j+1);

      }// nbins loop
      
      delete hPreUnfoldingSpectra;
      delete hAfterUnfoldingSpectra;
      delete hMCGen; 
      
    }// centrality loop

    cout<<"pp "<<endl;

    // now do it for the pp:
    TH1F * hPreUnfoldingSpectraPP = new TH1F("hPreUnfoldingSpectraPP","",nbins,0, 1000);
    TH1F * hAfterUnfoldingSpectraPP;
    
    for(int j = 0;j<nbins;++j){
	
      hPreUnfoldingSpectraPP->SetBinContent(j+1, random->Gaus(meanMeasPP[j], sigmaMeasPP[j]));
      hPreUnfoldingSpectraPP->SetBinError(j+1, sigmaMeasPP[j]/sqrt(unfoldingTrials));
      if(j+1 == dPP_Comb->FindBin(150)) hPP_beforeUnfold_Gaussian_pt150->Fill(random->Gaus(meanMeasPP[j], sigmaMeasPP[j]));
      
    }// nbins loop
    TH1F* hMCGenPP          = (TH1F*)mPP_Response->ProjectionX();
    removeZero(hMCGenPP);
    bayesianUnfold myUnfoldingMultiPP(mPP_Matrix, hMCGenPP, 0);
    myUnfoldingMultiPP.unfold(hPreUnfoldingSpectraPP, BayesIter);

    hAfterUnfoldingSpectraPP = (TH1F*) myUnfoldingMultiPP.hPrior->Clone("hAfterUnfoldingSpectraPP");

    for(int j = 0;j<nbins;++j){

      meanUnfoldPP[j][u] = hAfterUnfoldingSpectraPP->GetBinContent(j+1);
      sigmaUnfoldPP[j][u] = hAfterUnfoldingSpectraPP->GetBinError(j+1);

    }// nbins loop

    delete hPreUnfoldingSpectraPP;
    delete hAfterUnfoldingSpectraPP;
    delete hMCGenPP; 
    
  }// unfolding trials loop


  // Now that we have all the necesary values we need, lets proceed to fill a histogram with the mean values for each ptbin and get the corrected values.
  TH1F * hAfterUnfoldingptBinDistribution[nbins];
  TH1F * hCorrUnfoldingPbPb[nbins_cent];

  // we need to store one gaussian histogram in the root file which we can plot 
  TH1F * hPbPb_Gaussian_pt150[nbins_cent];
  TH1F * hPP_Gaussian_pt150;


  for(int i = 0;i<nbins_cent;++i){

    hCorrUnfoldingPbPb[i] = new TH1F(Form("PbPb_BayesianUnfolded_cent%d",i),"Spectra after correction", nbins,0, 1000);
    hPbPb_Gaussian_pt150[i] = new TH1F(Form("PbPb_Gaussian_pt150_cent%d",i),"gaussian distribution of values at pt bin at 150",1000, 0.1 * dPbPb_TrgComb[i]->GetBinContent(dPbPb_TrgComb[i]->FindBin(150)), 1.9 * dPbPb_TrgComb[i]->GetBinContent(dPbPb_TrgComb[i]->FindBin(150)));

    for(int j = 0;j<nbins;++j){
      
      hAfterUnfoldingptBinDistribution[j] = new TH1F(Form("hAfterUnfoldingptBinDistribution_ptBin%d",j),"",100,	0, 1);
      for(int u = 0;u<unfoldingTrials;++u){

	hAfterUnfoldingptBinDistribution[j]->Fill(meanUnfoldPbPb[j][i][u]);
	if(j+1 == dPbPb_TrgComb[i]->FindBin(150)) hPbPb_Gaussian_pt150[i]->Fill(meanUnfoldPbPb[j][i][u]);

      }// unfolding trials loop

      hCorrUnfoldingPbPb[i]->SetBinContent(j+1, hAfterUnfoldingptBinDistribution[j]->GetMean());
      hCorrUnfoldingPbPb[i]->SetBinError(j+1, hAfterUnfoldingptBinDistribution[j]->GetRMS());

      delete hAfterUnfoldingptBinDistribution[j];
      
    }// nbins loop

  }// centrality loop

  // similar for the pp:
  TH1F * hAfterUnfoldingptBinDistributionPP[nbins];
  TH1F * hCorrUnfoldingPP;
  
  hCorrUnfoldingPP = new TH1F("PP_BayesianUnfolded","Spectra after unfolding error correction",nbins,0, 1000);
  hPP_Gaussian_pt150 = new TH1F("PP_Gaussian_pt100","gaussian distribution of values at pt bin at 150",1000, 0.1 * dPP_Comb->GetBinContent(dPP_Comb->FindBin(150)) , 1.9 * dPP_Comb->GetBinContent(dPP_Comb->FindBin(150)));
  
  for(int j = 0;j<nbins;++j){
    
    hAfterUnfoldingptBinDistributionPP[j] = new TH1F(Form("hAfterUnfoldingptBinDistributionPP_ptBin%d",j),"",100, 0, 1);
    for(int u = 0;u<unfoldingTrials;++u){
      
      hAfterUnfoldingptBinDistributionPP[j]->Fill(meanUnfoldPP[j][u]);
      if(j+1 == dPP_Comb->FindBin(150)) hPP_Gaussian_pt150->Fill(meanUnfoldPP[j][u]);
      
    }// unfolding trials loop
    
    hCorrUnfoldingPP->SetBinContent(j+1, hAfterUnfoldingptBinDistributionPP[j]->GetMean());
    hCorrUnfoldingPP->SetBinError(j+1, hAfterUnfoldingptBinDistributionPP[j]->GetRMS());
    
    delete hAfterUnfoldingptBinDistributionPP[j];
    
  }// nbins loop
    
  TFile f(Form("July20/HiForest_%disATLASCut_%ddo10GeVBins_data_driven_correction_ak%d.root" , isATLASCut, do10GeVBins, radius),"RECREATE");
  f.cd();

  for(int i = 0;i<nbins_cent;i++) {

    //hCorrUnfoldingPbPb[i] = (TH1F*)hCorrUnfoldingPbPb[i]->Rebin(nbins_coarse, Form("PbPb_BayesianUnfolded_cent%d",i), ptbins_long_coarse);
    //divideBinWidth(hCorrUnfoldingPbPb[i]);
    //dPbPb_TrgComb[i] = (TH1F*)dPbPb_TrgComb[i]->Rebin(nbins_coarse, Form("PbPb_measured_cent%d",i), ptbins_long_coarse);
    //divideBinWidth(dPbPb_TrgComb[i]);

    hMinBias[i]->Write();
    hDataBeforeSub[i]->Write();
    hDataAfterSub[i]->Write();
    
    hCorrUnfoldingPbPb[i]->Scale(166 * 1e9);
    hCorrUnfoldingPbPb[i]->Write();
    hCorrUnfoldingPbPb[i]->Print("base");

    dPbPb_TrgComb[i]->Scale(166 * 1e9);
    dPbPb_TrgComb[i]->SetName(Form("PbPb_data_minbiasSub_cent%d",i));
    //dPbPb_TrgComb[i]->Scale(145.156 * 1e9);
    dPbPb_TrgComb[i]->Write();
    dPbPb_TrgComb[i]->Print("base");

    hPbPb_beforeUnfold_Gaussian_pt150[i]->Write();
    hPbPb_beforeUnfold_Gaussian_pt150[i]->Print("base");
    
    hPbPb_Gaussian_pt150[i]->Write();
    hPbPb_Gaussian_pt150[i]->Print("base");

    mPbPb_Matrix[i]->Write();
    
  }

  //hCorrUnfoldingPP = (TH1F*)hCorrUnfoldingPP->Rebin(nbins_coarse, "PP_BayesianUnfolded", ptbins_long_coarse);
  //divideBinWidth(hCorrUnfoldingPP);
  //dPP_Comb = (TH1F*)dPP_Comb->Rebin(nbins_coarse, "PP_measured", ptbins_long_coarse);  
  //divideBinWidth(dPP_Comb);
  
  hCorrUnfoldingPP->Scale(5.3 * 1e9);
  hCorrUnfoldingPP->Write();
  hCorrUnfoldingPP->Print("base");

  dPP_Comb->Scale(5.3 * 1e9);
  dPP_Comb->Write();
  dPP_Comb->Print("base");

  hPP_beforeUnfold_Gaussian_pt150->Write();
  hPP_beforeUnfold_Gaussian_pt150->Print("base");
  
  hPP_Gaussian_pt150->Write();
  hPP_Gaussian_pt150->Print("base");
  mPP_Matrix->Write();

  f.Write();
  f.Close();

  // make the data driven Error correction histograms and plots here:
  
  TH1F * hError_Meas[nbins_cent+1], * hError_Fixed[nbins_cent+1];
  for(int i = 0; i<nbins_cent+1; ++i){
    cout<<"centrality "<<i<<endl;
    if(i < nbins_cent){
      hError_Meas[i]  = new TH1F(Form("hError_Meas_cent%d",i),"",nbins, 0, 1000);
      hError_Fixed[i] = new TH1F(Form("hError_Fixed_cent%d",i),"",nbins, 0, 1000);
    }
    if(i == nbins_cent){
      hError_Meas[i]  = new TH1F(Form("hError_PP_Meas_cent%d",i),"",nbins, 0, 1000);
      hError_Fixed[i] = new TH1F(Form("hError_PP_Fixed_cent%d",i),"",nbins, 0, 1000);
    }

    for(int j = 1; j<=nbins; ++j){
      //cout<<"ptbins "<<j<<endl;
      if(i < nbins_cent){
	if(dPbPb_TrgComb[i]->GetBinContent(j)!=0) hError_Meas[i]->SetBinContent(j, (float)dPbPb_TrgComb[i]->GetBinError(j)/dPbPb_TrgComb[i]->GetBinContent(j));
        //hError_Meas[i]->SetBinContent(j, (float)dPbPb_TrgComb[i]->GetBinError(j));
	if(hCorrUnfoldingPbPb[i]->GetBinContent(j)!=0)hError_Fixed[i]->SetBinContent(j, (float)hCorrUnfoldingPbPb[i]->GetBinError(j)/hCorrUnfoldingPbPb[i]->GetBinContent(j));
	cout<<j<<" "<<hError_Fixed[i]->GetBinContent(j)<<endl;
	//hError_Fixed[i]->SetBinContent(j, (float)hCorrUnfoldingPbPb[i]->GetBinError(j));
      }
      if(i == nbins_cent){
	hError_Meas[i]->SetBinContent(j, (float)dPP_Comb->GetBinError(j)/dPP_Comb->GetBinContent(j));
	//hError_Meas[i]->SetBinContent(j, (float)dPP_Comb->GetBinError(j));
	hError_Fixed[i]->SetBinContent(j, (float)hCorrUnfoldingPP->GetBinError(j)/hCorrUnfoldingPP->GetBinContent(j));
	//hError_Fixed[i]->SetBinContent(j, (float)hCorrUnfoldingPP->GetBinError(j));
      }
    }
    
    hError_Meas[i]->SetAxisRange(50, 299, "X");
    //hError_Meas[i]->Print("base");
    //hError_Fixed[i]->Print("base");
    //hError_Meas[i]->SetAxisRange(1e-12, 1, "Y");
  }
  //cout<<" passed the loop"<<endl;

  TCanvas * cSpectra = new TCanvas("cSpectra","",1200,1000);
  makeMultiPanelCanvas(cSpectra,3,2,0.0,0.0,0.2,0.15,0.07);  
  for(int i = 0; i<nbins_cent; ++i){
    //cout<<i<<endl;
    cSpectra->cd(nbins_cent-i);
    cSpectra->cd(nbins_cent-i)->SetLogy();
    dPbPb_TrgComb[i]->SetMarkerStyle(24);
    dPbPb_TrgComb[i]->SetMarkerColor(kBlack);
    makeHistTitle(dPbPb_TrgComb[i]," ","jet pT","dN/dpT");
    dPbPb_TrgComb[i]->SetAxisRange(50, 299, "X");
    dPbPb_TrgComb[i]->Draw("p");

    hCorrUnfoldingPbPb[i]->SetMarkerStyle(33);
    hCorrUnfoldingPbPb[i]->SetMarkerColor(kRed);
    hCorrUnfoldingPbPb[i]->Draw("psame");

  }
  TLegend * Spec = myLegend(0.55,0.55,0.75,0.75);
  cSpectra->cd(1);
  putCMSPrel();
  Spec->AddEntry(dPbPb_TrgComb[0],"Measured","pl");
  Spec->AddEntry(hCorrUnfoldingPbPb[0],"Data Driven Correction","pl");
  Spec->SetTextSize(0.04);
  Spec->Draw();
  
  cSpectra->SaveAs(Form("%sUnfoldingSpectra_fromDataDrivenMacro_PbPb_%s_R%d_%d_hiForest_%dGeVCut.pdf",outLocation,etaWidth,radius,date.GetDate(),unfoldingCut),"RECREATE");

  
  TCanvas * cErrorFix = new TCanvas("cErrorFix","",1200,1000);
  makeMultiPanelCanvas(cErrorFix,3,2,0.0,0.0,0.2,0.15,0.07);  
  for(int i = 0; i<nbins_cent; ++i){
    //cout<<i<<endl;
    cErrorFix->cd(nbins_cent-i);
    cErrorFix->cd(nbins_cent-i)->SetLogy();
    makeHistTitle(hError_Meas[i]," ","jet pT","Error/Content");
    hError_Meas[i]->SetMarkerStyle(24);
    hError_Meas[i]->SetMarkerColor(kBlack);
    hError_Meas[i]->Draw("p");

    hError_Fixed[i]->SetMarkerStyle(33);
    hError_Fixed[i]->SetMarkerColor(kRed);
    hError_Fixed[i]->Draw("psame");

  }
  TLegend * err = myLegend(0.55,0.55,0.75,0.75);
  cErrorFix->cd(1);
  putCMSPrel();
  err->AddEntry(hError_Meas[0],"Measured","pl");
  err->AddEntry(hError_Fixed[0],"Data Driven Correction","pl");
  err->SetTextSize(0.04);
  err->Draw();
  
  cErrorFix->SaveAs(Form("%sUnfoldingErrorFix_fromDataDrivenMacro_PbPb_%s_R%d_%d_hiForest_%dGeVCut.pdf",outLocation,etaWidth,radius,date.GetDate(),unfoldingCut),"RECREATE");

  TCanvas * cErrorFixPP = new TCanvas("cErrorFixPP","",800,600);
  cErrorFixPP->SetLogy();
  hError_Meas[nbins_cent]->SetMarkerStyle(24);
  hError_Meas[nbins_cent]->SetMarkerColor(kBlack);
  hError_Meas[nbins_cent]->Draw("p");
  
  hError_Fixed[nbins_cent]->SetMarkerStyle(33);
  hError_Fixed[nbins_cent]->SetMarkerColor(kRed);
  hError_Fixed[nbins_cent]->Draw("psame");
  TLegend * errPP = myLegend(0.55,0.55,0.75,0.75);
  putCMSPrel();
  errPP->AddEntry(hError_Meas[nbins_cent],"Measured","pl");
  errPP->AddEntry(hError_Fixed[nbins_cent],"Data Driven Correction","pl");
  errPP->SetTextSize(0.04);
  errPP->Draw();
  
  cErrorFixPP->SaveAs(Form("%sUnfoldingErrorFix_fromDataDrivenMacro_PP_%s_R%d_%d_hiForest_%dGeVCut.pdf",outLocation,etaWidth,radius,date.GetDate(),unfoldingCut),"RECREATE");
  
  timer.Stop();
  if(printDebug) cout<<"CPU time (mins) = "<<(Float_t)timer.CpuTime()/60<<endl;
  if(printDebug) cout<<"Real tile (mins) = "<<(Float_t)timer.RealTime()/60<<endl;


}
