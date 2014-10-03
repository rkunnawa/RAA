// Raghav Kunnawalkam Elayavalli
// June 11 2014
// CERN

// RAA - plotting macro. pretty much all the neccessary plots will be made here in this macro. 
// will take input from the read macro and the analysis macro. Maybe not from the read macro - i will have to decide on that. 

// June 25th - just started working on the macro. will plot the iteration systematics, MC closure test, create setup for Unfolded vs measured, RAA and Normalized Response matrix plots as well. 

// July 1st - finished the macro. it would have been very easy and efficient to put them all in 2 segments, 
//            one for PbPb (with the centrality loop) and one for pp. 
//          - But ive decided to keep it this way since we can easily delete any segments which are complete and remake any plots individually. 

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


using namespace std;

void RAA_plot(int radius = 3, char *algo = "Vs", char *jet_type = "Calo"){

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

  TFile *fin; 
  /*
  //if(location=="MIT") 
  fin= TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_unfo_ak%s%d%s_20140912.root",algo,radius,jet_type));
  //if(location=="CERN")fin= TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_unfo_ak%s%d%s_20140911.root",algo,radius,jet_type));
  //if(location=="MPB") fin= TFile::Open(Form(""))
  

  TH1F *dPbPb_TrgComb[nbins_cent+1], *dPbPb_Comb[nbins_cent+1], *dPbPb_Trg80[nbins_cent+1], *dPbPb_Trg65[nbins_cent+1], *dPbPb_Trg55[nbins_cent+1], *dPbPb_1[nbins_cent+1], *dPbPb_2[nbins_cent+1], *dPbPb_3[nbins_cent+1], *dPbPb_80[nbins_cent+1], *dPbPb_65[nbins_cent+1], *dPbPb_55[nbins_cent+1];
  
  TH1F *mPbPb_Gen[nbins_cent+1], *mPbPb_Reco[nbins_cent+1];
  TH2F *mPbPb_Matrix[nbins_cent+1], *mPbPb_Response[nbins_cent+1], *mPbPb_ResponseNorm[nbins_cent+1];
  TH1F *mPbPb_mcclosure_data[nbins_cent+1];
  
  TH1F *dPP_1, *dPP_2, *dPP_3, *dPP_Comb;
  
  TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix, *mPP_Response;
  TH2F *mPP_ResponseNorm;
  TH1F *mPP_mcclosure_data;
  
  const int Iterations = 20; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1], *uPbPb_BinByBin[nbins_cent+1]; 
  TH1F *uPbPb_BayesianIter[Iterations][nbins_cent+1];
  TH1F *uPbPb_MC_Bayes[nbins_cent];
  TH1F *uPbPb_MC_BayesianIter[Iterations][nbins_cent+1];
  TH1F *uPbPb_MC_BinByBin[nbins_cent+1];

  TH1F *uPP_Bayes, *uPP_BinByBin;
  TH1F *uPP_BayesianIter[Iterations];
  TH1F *uPP_MC_Bayes, *uPP_MC_BinByBin;
  TH1F *uPP_MC_BayesianIter[Iterations];

  TH1F *RAA_measured[nbins_cent+1];
  TH1F *RAA_binbybin[nbins_cent+1];
  TH1F *RAA_bayesian[nbins_cent+1];

  for(int i = 0;i<=nbins_cent;i++){
    
    cout<<"cent = "<<i<<endl;
    dPbPb_TrgComb[i] = (TH1F*)fin->Get(Form("hpbpb_TrgObjComb_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Reco[i] = (TH1F*)fin->Get(Form("hpbpb_reco_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_Gen[i] = (TH1F*)fin->Get(Form("hpbpb_gen_R%d_n20_eta_p20_cent%d",radius,i));
    mPbPb_ResponseNorm[i] = (TH2F*)fin->Get(Form("mPbPb_ResponseNorm_cent%d",i));
    mPbPb_Response[i] = (TH2F*)fin->Get(Form("mPbPb_Response_cent%d",i));
    uPbPb_Bayes[i] = (TH1F*)fin->Get(Form("uPbPb_Bayes_cent%d",i));
    uPbPb_MC_Bayes[i] = (TH1F*)fin->Get(Form("uPbPb_MC_Bayes_cent%d",i));
    uPbPb_BinByBin[i] = (TH1F*)fin->Get(Form("uPbPb_BinByBin_cent%d",i));
    uPbPb_MC_BinByBin[i] = (TH1F*)fin->Get(Form("uPbPb_MC_BinByBin_cent%d",i));
    RAA_bayesian[i] = (TH1F*)fin->Get(Form("RAA_bayesian_cent%d",i));
    RAA_binbybin[i] = (TH1F*)fin->Get(Form("RAA_binbybin_cent%d",i));
    RAA_measured[i] = (TH1F*)fin->Get(Form("RAA_measured_cent%d",i));
    mPbPb_mcclosure_data[i] = (TH1F*)fin->Get(Form("mPbPb_mclosure_data_cent%d",i));
    
    for(int j = 0;j<Iterations;j++){
      uPbPb_BayesianIter[j][i] = (TH1F*)fin->Get(Form("uPbPb_BayesianIter%d_cent%d",j,i));
      uPbPb_MC_BayesianIter[j][i] = (TH1F*)fin->Get(Form("uPbPb_MC_BayesianIter%d_cent%d",j,i));
    }

  }

  dPP_Comb = (TH1F*)fin->Get(Form("hpp_TrgComb_R%d_n20_eta_p20",radius));
  mPP_ResponseNorm = (TH2F*)fin->Get(Form("mPP_ResponseNorm",radius));	
  mPP_mcclosure_data = (TH1F*)fin->Get(Form("mPP_mcclosure_data",radius));

  uPP_Bayes = (TH1F*)fin->Get("uPP_Bayes");
  uPP_BinByBin = (TH1F*)fin->Get("uPP_BinByBin");
  uPP_MC_Bayes = (TH1F*)fin->Get("uPP_MC_Bayes");
  uPP_MC_BinByBin = (TH1F*)fin->Get("uPP_MC_BinByBin");
  mPP_Gen = (TH1F*)fin->Get(Form("hpp_gen_R%d_n20_eta_p20",radius));
  mPP_Reco = (TH1F*)fin->Get(Form("hpp_reco_R%d_n20_eta_p20",radius));

  for(int i = 0;i<Iterations;i++){

  	uPP_BayesianIter[i] = (TH1F*)fin->Get(Form("uPP_BayesianIter%d",i));
  	uPP_MC_BayesianIter[i] = (TH1F*)fin->Get(Form("uPP_MC_BayesianIter%d",i));

  }


  //Ok now that we have loaded all the histograms we need - lets start making the plots 
  /*
  // line at 1
  TLine *line = new TLine(50,1,300,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // plot 0 - PbPb and pp Unfolded Spectra compared with Data(measured) and MC Spectra. in a 2 by 1 canvas. with the PbPb doing the multiply by powers of 10. unfolded in black circles, MC in black (dotted) line, measured in red open boxes. 
  TCanvas *cSpectra = new TCanvas("cSpectra","PbPb and PP spectra",1000,1000);
  cSpectra->Divide(2,1);
  
  cSpectra->cd(1);
  cSpectra->cd(1)->SetLogy();
  cSpectra->cd(1)->SetLogx();
  //cSpectra->Set
  Double_t scaleFactor[nbins_cent+1] = {1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6};
  
  TLegend *PbPbSpectra = myLegend(0.4,0.65,0.85,0.9);
  
  for(int i = 0;i<=nbins_cent;i++){
    
    uPbPb_Bayes[i]->Scale(1./(scaleFactor[i]*1e6));
    dPbPb_TrgComb[i]->Scale(1./(scaleFactor[i]*1e6));
    mPbPb_Reco[i]->Scale(1./scaleFactor[i]);

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
      PbPbSpectra->AddEntry(uPbPb_Bayes[i],Form("0 - 200 cent * 10^{%d}",i),"pl");
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
  
  
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  //plot 1 - PbPb iteration systematics. 
  // this will be a 3 by 2 panel plot showing bayesian for PbPb. per centrality bin. 
  // divide different unfolding iterations with iteration 4 - the nominal one. 
  TCanvas *cIterSysPbPb = new TCanvas("cIterSysPbPb","PbPb Iteration systematics",1200,800);
  makeMultiPanelCanvasWithGap(cIterSysPbPb,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
 
  TLegend *PbPb_itersys = myLegend(0.53,0.65,0.85,0.9);
  
  for(int i = 0;i<nbins_cent;i++){
    cIterSysPbPb->cd(nbins_cent-i);
    line->Draw();

    drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

    for(int j = 2;j<7;j++){
  		
      uPbPb_BayesianIter[j][i]->Divide(uPbPb_Bayes[i]);
      uPbPb_BayesianIter[j][i]->SetMarkerStyle(2);
      uPbPb_BayesianIter[j][i]->SetMarkerColor(j);
      if(j==2){
	makeHistTitle(uPbPb_BayesianIter[j][i]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal)");
	uPbPb_BayesianIter[j][i]->Draw();
      }else if(j==4);
      else {
	uPbPb_BayesianIter[j][i]->Draw("same");
      }
      if(i==0) PbPb_itersys->AddEntry(uPbPb_BayesianIter[j][i],Form("Iteration %d",j),"pl");

    }


  }
  PbPb_itersys->Draw();

  cIterSysPbPb->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s %s Jets R=0.%d",algo, jet_type, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cIterSysPbPb->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_iteration_systematics_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  
  //plot 2 - pp iteration systematics 
  // this is just one panel plot showing bayesian for pp. 
  TCanvas *cIterSysPP = new TCanvas("cIterSysPP","PP Iteration systematics",600,400);

  TLegend *PP_itersys = myLegend(0.53,0.65,0.85,0.9);
  line->Draw();

  for(int i = 2;i<7;i++){
    uPP_BayesianIter[i]->Divide(uPP_Bayes);
    uPP_BayesianIter[i]->SetMarkerStyle(2);
    uPP_BayesianIter[i]->SetMarkerColor(i);
    if(i==2){
      makeHistTitle(uPP_BayesianIter[i]," ","Jet p_{T} (GeV/c)","Ratio (unfolded/Nominal)");
      uPP_BayesianIter[i]->Draw();
    }else if(i==4);
    else uPP_BayesianIter[i]->Draw("same");
    PP_itersys->AddEntry(uPP_BayesianIter[i],Form("Iteration %d",i),"pl");
  }

  PP_itersys->Draw();
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Jets R=0.%d", jet_type,radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cIterSysPP->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_iteration_systematics_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  
  //plot 3 - RAA 
  // again this will be a 6 panel plot. showing measured, unfolded Bayesian, and unfolded Bin By Bin methods. 
  TCanvas *cRAA = new TCanvas("cRAA","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA = myLegend(0.53,0.65,0.85,0.9);

  for(int i = 0;i<nbins_cent;i++){

    cRAA->cd(nbins_cent-i);
    drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

    RAA_measured[i]->SetMarkerColor(kBlack);
    RAA_measured[i]->SetMarkerStyle(24);
    RAA_measured[i]->Draw();

    RAA_bayesian[i]->SetMarkerColor(kRed);
    RAA_bayesian[i]->SetMarkerStyle(33);
    RAA_bayesian[i]->Draw("same");

    RAA_binbybin[i]->SetMarkerStyle(29);
    RAA_binbybin[i]->SetMarkerColor(kBlue);
    RAA_binbybin[i]->Draw("same");

    line->Draw();

  }

  tRAA->AddEntry(RAA_measured[0],"No Unfolding","pl");
  tRAA->AddEntry(RAA_bayesian[0],"Bayesian","pl");
  tRAA->AddEntry(RAA_binbybin[0],"BinbyBin","pl");

  cRAA->cd(0);
  tRAA->Draw();
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Particle Flow Jets R=0.%d",algo, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cRAA->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/RAA_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 4 - PbPb MC closure test 
  // this will also be a 6 panel plot showing bayesian iteration 4 and binbybin divided by measured which is actually the MC
  TCanvas *cPbPbMCclosure = new TCanvas("cPbPbMCclosure","PbPb MC closure test",1200,800);
  makeMultiPanelCanvasWithGap(cPbPbMCclosure,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TH1F *hMCClosurePbPb_Meas[nbins_cent+1], *hMCClosurePbPb_Bayesian[nbins_cent+1], *hMCClosurePbPb_BinByBin[nbins_cent+1];
  for(int i = 0;i<nbins_cent;i++){

    cPbPbMCclosure->cd(nbins_cent-i);
    drawText(Form("%2.0f-%2.0f",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);
    //hMCClosurePbPb_Meas[i] = (TH1F*)mPbPb_Reco[i]->Clone(Form("hMCClosurePbPb_Meas_cent%d",i));
    hMCClosurePbPb_Meas[i] = (TH1F*)mPbPb_mcclosure_data[i]->Clone(Form("hMCClosurePbPb_Meas_cent%d",i));
    hMCClosurePbPb_Bayesian[i] = (TH1F*)uPbPb_MC_Bayes[i]->Clone(Form("hMCClosurePbPb_Bayesian_cent%d",i));
    hMCClosurePbPb_BinByBin[i] = (TH1F*)uPbPb_MC_BinByBin[i]->Clone(Form("hMCClosurePbPb_BinByBin_cent%d",i));

    hMCClosurePbPb_Meas[i]->Divide(mPbPb_mcclosure_data[i]);
    hMCClosurePbPb_Bayesian[i]->Divide(mPbPb_mcclosure_data[i]);
    hMCClosurePbPb_BinByBin[i]->Divide(mPbPb_mcclosure_data[i]);

    makeHistTitle(hMCClosurePbPb_Meas[i]," ","Jet p_{T} (GeV/c)","Reco/Truth");
    hMCClosurePbPb_Meas[i]->SetMarkerStyle(24);
    hMCClosurePbPb_Meas[i]->SetMarkerColor(kBlack);
    hMCClosurePbPb_Meas[i]->Draw();

    hMCClosurePbPb_Bayesian[i]->SetMarkerStyle(33);
    hMCClosurePbPb_Bayesian[i]->SetMarkerColor(kRed);
    hMCClosurePbPb_Bayesian[i]->Draw("same");

    hMCClosurePbPb_BinByBin[i]->SetMarkerColor(kBlue);
    hMCClosurePbPb_BinByBin[i]->SetMarkerStyle(29);
    hMCClosurePbPb_BinByBin[i]->Draw("same");

    line->Draw();

  }

  cPbPbMCclosure->cd(1);
  TLegend *pbpbmcclosure = myLegend(0.53,0.65,0.85,0.9);
  pbpbmcclosure->AddEntry(hMCClosurePbPb_Meas[0],"PbPb no unfolding","pl");
  pbpbmcclosure->AddEntry(hMCClosurePbPb_Bayesian[0],"PbPb Bayesian 4 Iter","pl");
  pbpbmcclosure->AddEntry(hMCClosurePbPb_BinByBin[0],"PbPb BinbyBin","pl");
  pbpbmcclosure->Draw();
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Particle Flow Jets R=0.%d",algo, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  
  cPbPbMCclosure->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_unfoldng_mc_closure_test_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 5 - PP MC closure test
  TCanvas * cPPMCclosure = new TCanvas("cPPMCclosure","PP MC closure test",600,400);
  //TH1F *hMCClosurePP_Meas = (TH1F*)mPP_Reco->Clone("hMCClosurePP_Meas");
  TH1F *hMCClosurePP_Meas = (TH1F*)mPP_mcclosure_data->Clone("hMCClosurePP_Meas");
  TH1F *hMCClosurePP_Bayesian = (TH1F*)uPP_MC_Bayes->Clone("hMCClosurePP_Bayesian");
  TH1F *hMCClosurePP_BinbyBin = (TH1F*)uPP_MC_BinByBin->Clone("hMCClosurePP_BinbyBin");

  hMCClosurePP_Bayesian->Divide(mPP_mcclosure_data);
  hMCClosurePP_Meas->Divide(mPP_mcclosure_data);
  hMCClosurePP_BinbyBin->Divide(mPP_mcclosure_data);

  makeHistTitle(hMCClosurePP_Meas," ","Jet p_{T} (GeV/c)","Reco/Truth");
  hMCClosurePP_Meas->SetAxisRange(50,300,"X");
  hMCClosurePP_Meas->SetAxisRange(0,2,"Y");

  hMCClosurePP_Meas->SetMarkerStyle(24);
  hMCClosurePP_Meas->SetMarkerColor(kBlack);
  hMCClosurePP_Meas->Draw();

  hMCClosurePP_Bayesian->SetMarkerStyle(33);
  hMCClosurePP_Bayesian->SetMarkerColor(kRed);
  hMCClosurePP_Bayesian->Draw("same");

  hMCClosurePP_BinbyBin->SetMarkerStyle(29);
  hMCClosurePP_BinbyBin->SetMarkerColor(kBlue);
  hMCClosurePP_BinbyBin->Draw("same");

  line->Draw();

  TLegend *ppmcclosure = myLegend(0.53,0.65,0.85,0.9);
  ppmcclosure->AddEntry(hMCClosurePP_Meas,"pp no unfolding","pl");
  ppmcclosure->AddEntry(hMCClosurePP_Bayesian,"pp Bayesian","pl");
  ppmcclosure->AddEntry(hMCClosurePP_BinbyBin,"pp BinbyBin","pl");
  ppmcclosure->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,22);
  
  cPPMCclosure->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_unfoldng_mc_closure_test_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  
  //plot 6 - data vs MC for PbPb 
  // i want to make it like ratio plots for each centrality class. and plot the unfolded ratio as well. 
  // so i want a line at 1, and 2 overlayed plots showing measured/Gen ratio and unfo/gen ratio. 
  TCanvas *cPbPb_data_vs_mc = new TCanvas("cPbPb_data_vs_mc","PbPb data vs mc",1200,800);
  makeMultiPanelCanvasWithGap(cPbPb_data_vs_mc,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TH1F *PbPb_meas_vs_mc[nbins_cent+1];
  TH1F *PbPb_unfo_vs_mc[nbins_cent+1];

  for(int i = 0;i<nbins_cent;i++){

    cPbPb_data_vs_mc->cd(nbins_cent-i);

    PbPb_meas_vs_mc[i] = (TH1F*)dPbPb_TrgComb[i]->Clone("PbPb_meas_vs_mc");
    PbPb_unfo_vs_mc[i] = (TH1F*)uPbPb_Bayes[i]->Clone("PbPb_unfo_vs_mc");

    PbPb_meas_vs_mc[i]->Divide(mPbPb_Reco[i]);
    PbPb_unfo_vs_mc[i]->Divide(mPbPb_Reco[i]);

    makeHistTitle(PbPb_meas_vs_mc[i]," ","Jet p_{T} (GeV/c)", "data/mc");
    PbPb_meas_vs_mc[i]->SetMarkerStyle(25);
    PbPb_meas_vs_mc[i]->SetMarkerColor(kBlack);
    PbPb_meas_vs_mc[i]->Draw();

    PbPb_unfo_vs_mc[i]->SetMarkerStyle(27);
    PbPb_unfo_vs_mc[i]->SetMarkerColor(kBlue);
    PbPb_unfo_vs_mc[i]->Draw("same");

    line->Draw();

    drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

  }

  cPbPb_data_vs_mc->cd(1);
  TLegend *PbPb_datavsmc = myLegend(0.53,0.65,0.85,0.9);
  PbPb_datavsmc->AddEntry(PbPb_meas_vs_mc[0],"no unfolding","pl");
  PbPb_datavsmc->AddEntry(PbPb_unfo_vs_mc[0],"Bayesian 4 Iter","pl");
  PbPb_datavsmc->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Particle Flow Jets R=0.%d",algo, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPbPb_data_vs_mc->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_data_vs_mc_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot7 - data vs MC for pp
  // made similar to above. 
  TCanvas *cPP_data_vs_mc = new TCanvas("cPP_data_vs_mc","PP data vs MC",600,400);

  TH1F *PP_meas_vs_mc = (TH1F*)dPP_Comb->Clone("PP_meas_vs_mc");
  TH1F *PP_unfo_vs_mc = (TH1F*)uPP_Bayes->Clone("PP_unfo_vs_mc");

  PP_meas_vs_mc->Divide(mPP_Reco);
  PP_unfo_vs_mc->Divide(mPP_Reco);

  makeHistTitle(PP_meas_vs_mc," ","Jet p_{T} (GeV/c)","data/mc");
  PP_meas_vs_mc->SetMarkerStyle(25);
  PP_meas_vs_mc->SetMarkerColor(kBlack);
  PP_meas_vs_mc->Draw();
  PP_unfo_vs_mc->SetMarkerStyle(27);
  PP_unfo_vs_mc->SetMarkerColor(kBlue);
  PP_unfo_vs_mc->Draw("same");

  line->Draw();

  TLegend *PP_datavsmc = myLegend(0.53,0.65,0.85,0.9);
  PP_datavsmc->AddEntry(PP_meas_vs_mc,"no unfolding","pl");
  PP_datavsmc->AddEntry(PP_unfo_vs_mc,"Bayesian 4 Iter","pl");
  PP_datavsmc->SetTextSize(0.02);
  PP_datavsmc->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPP_data_vs_mc->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_data_vs_mc_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot8 - normalized Response Matrix for PbPb. 
  TCanvas *cPbPb_NormResMat = new TCanvas("cPbPb_NormResMat","Normalized Response Matrix for PbPb",1200,800);
  makeMultiPanelCanvasWithGap(cPbPb_NormResMat,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0;i<nbins_cent;i++){
    cPbPb_NormResMat->cd(nbins_cent-i);

    makeHistTitle(mPbPb_ResponseNorm[i]," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");
    drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

    mPbPb_ResponseNorm[i]->SetAxisRange(1e-10,1,"Z");
    mPbPb_ResponseNorm[i]->Draw("colz");

  }
  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow %s Jets R=0.%d",algo,radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPbPb_NormResMat->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_normalized_response_matrix_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot9 - normalized response matrix for pp
  TCanvas *cPP_NormResMat = new TCanvas("cPP_NormResMat","Normalized Response Matrix fr PP",600,400);

  makeHistTitle(mPP_ResponseNorm," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");

  mPP_ResponseNorm->SetAxisRange(1e-10,1,"Z");
  mPP_ResponseNorm->Draw("colz");

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPP_NormResMat->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_normalized_response_matrix_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  

  //plot10 - Cross section for PbPb - in proper units per centrality bin.  
  TCanvas *cPbPb_sigma = new TCanvas("cPbPb_sigma","PbPb inclusive jet invariant cross section",1200,800);
  makeMultiPanelCanvasWithGap(cPbPb_sigma,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0;i<nbins_cent;i++){

    cPbPb_sigma->cd(nbins_cent-i);
    cPbPb_sigma->cd(nbins_cent-i)->SetLogy();
    cPbPb_sigma->cd(nbins_cent-i)->SetLogx();
    
    drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.85,0.8,20);

    makeHistTitle(dPbPb_TrgComb[i]," ","Jet p_{T} (GeV/c)","#frac{d^2 #sigma}{d#eta dp_{T}} micro barns");
    dPbPb_TrgComb[i]->SetMarkerStyle(24);
    dPbPb_TrgComb[i]->SetMarkerColor(kBlack);
    dPbPb_TrgComb[i]->SetAxisRange(20,600,"X");
    dPbPb_TrgComb[i]->Draw();

    uPbPb_Bayes[i]->SetMarkerStyle(33);
    uPbPb_Bayes[i]->SetMarkerColor(kRed);
    uPbPb_Bayes[i]->Draw("same");

    uPbPb_BinByBin[i]->SetMarkerStyle(29);
    uPbPb_BinByBin[i]->SetMarkerColor(kBlue);
    uPbPb_BinByBin[i]->Draw("same");

  }

  cPbPb_sigma->cd(1);
  TLegend *PbPb_sigma = myLegend(0.53,0.65,0.85,0.9);
  PbPb_sigma->AddEntry(dPbPb_TrgComb[0],"No Unfolding","pl");
  PbPb_sigma->AddEntry(uPbPb_Bayes[0],"Bayesian 4 Iter","pl");
  PbPb_sigma->AddEntry(uPbPb_BinByBin[0],"BinByBin","pl");
  PbPb_sigma->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow %s Jets R=0.%d",algo,radius),0.3,0.95,15);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPbPb_sigma->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_invariant_cross_section_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  

  //plot 11 - cross section for PP 
  TCanvas *cPP_sigma = new TCanvas("cPP_sigma","PP inclusive jet invariant cross section",600,400);
  cPP_sigma->SetLogy();

  makeHistTitle(dPP_Comb,"","Jet p_{T} (GeV/c)","#frac{d^2 #sigma}{d#eta dp_{T}} nano barns");
  dPP_Comb->SetMarkerStyle(24);
  dPP_Comb->SetMarkerColor(kBlack);
  dPP_Comb->Draw();

  uPP_Bayes->SetMarkerColor(kRed);
  uPP_Bayes->SetMarkerStyle(33);
  uPP_Bayes->Draw("same");

  uPP_BinByBin->SetMarkerStyle(29);
  uPP_BinByBin->SetMarkerColor(kBlue);
  uPP_BinByBin->Draw("same");

  TLegend *PP_sigma = myLegend(0.53,0.65,0.85,0.9);
  PP_sigma->AddEntry(dPP_Comb,"No Unfolding","pl");
  PP_sigma->AddEntry(uPP_Bayes,"Bayesian 4 Iter","pl");
  PP_sigma->AddEntry(uPP_BinByBin,"BinByBin","pl");
  PP_sigma->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPbPb_sigma->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PP_invariant_cross_section_ak%d%s_%d.pdf",radius,jet_type,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  

  // plot 12 - Gen Spectra for PbPb and pp. 
  TCanvas *cGenSpectra = new TCanvas("cGenSpectra","Generator Level Spectra for PbPb and pp",1400,1000);
  
  cGenSpectra->Divide(2,1);
  
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
      mPbPb_Gen[i]->SetAxisRange(1e-15,1e5,"Y");
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

  cGenSpectra->cd(2);						
  cGenSpectra->cd(2)->SetLogy();
  cGenSpectra->cd(2)->SetLogx();
  makeHistTitle(mPP_Gen,"","MC: Jet p_{T} (GeV/c)","arbitrary for now");
  mPP_Gen->SetMarkerStyle(28);
  mPP_Gen->SetMarkerColor(kBlack);
  mPP_Gen->SetAxisRange(1e-14,1,"Y");
  mPP_Gen->SetAxisRange(15,600,"X");
  mPP_Gen->Draw();

  drawText("#frac{chMax}{jtpt}>0.01",0.15,0.2,20);
  drawText(Form("Anti-k_{T} %s Jets R=0.%d",jet_type,radius),0.45,0.9,20);
  drawText("2.76 TeV, |#eta|<2, |vz|<15",0,0.9,20);
  drawText("pp",0.3,0.8,20);

  mPP_Reco->SetMarkerStyle(29);
  mPP_Reco->SetMarkerColor(kBlack);
  mPP_Reco->Draw("same");
  
  TLegend *GenSpectra = myLegend(0.15,0.3,0.35,0.5);
  GenSpectra->AddEntry(mPP_Gen,"MC GenJet - refpt","pl");
  GenSpectra->AddEntry(mPP_Reco,"MC RecoJet - jtpt","pl");
  GenSpectra->SetTextSize(0.04);
  GenSpectra->Draw();
  
  cGenSpectra->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/MC_spectra_ak%s%d%s_%d.pdf",algo,radius,jet_type,date.GetDate()),"RECREATE");
 
  //unscale the gen and reco histograms: 
  for(int i = 0;i<=nbins_cent;i++){ 
    mPbPb_Gen[i]->Scale(1./scaleFactor[i]);
    mPbPb_Reco[i]->Scale(1./scaleFactor[i]); 
  }
 
  */
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Plotting the average energy subtracted in the Vs cone 
  if(algo=="Vs"){

    TFile *fMCin = TFile::Open("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_pp_mc_akVsCalo_20141003.root");
    TH1F *hPbPb_jtpu[no_radius][nbins_eta][nbins_cent+1];
    
    for(int i = 0;i<=nbins_cent;i++){

      for(int j = 0;j<nbins_eta;j++){

	for(int k = 0;k<no_radius;k++){

	  hPbPb_jtpu[i][j][k] = (TH1F*)fMCin->Get(Form("hpbpb_jtpu_R%d_%s_cent%d",list_radius[k],etaWidth[j],i));

	}// radius loop
	
      }//eta bins loop

    }//centrality loop

    TCanvas *cPbPb_MC_jtpu[nbins_eta];
    
    for(int j = 0;j<nbins_eta;j++){

      cPbPb_MC_jtpu[j] = new TCanvas(Form("cPbPb_MC_jtpu_%s",etaWidth[j]),FOrm("energy subtracted from Jets in the Vs algorithmin the range %s",etaWidth[j]),1000,800);
      makeMultiPanelCanvasWithGap(cPbPb_MC_jtpu[j],3,2,0.01,0.01,0.16,0.2,0.04,0.04);
      
      TLegend *LPbPb_MC_jtpu = myLegend(0.53,0.65,0.85,0.9);
      
      for(int i = 0;i<nbins_cent;i++){

	cPbPb_MC_jtpu[j]->cd(i);
	cPbPb_MC_jtpu[j]->cd(i)->SetLogy();

	for(int k = 1;k<4;k++){

	  hPbPb_jtpu[k][j][i]->Set

	}//radius loop

      }//centrality loop

      cPbPb_MC_jtpu[j]->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/bkg_energy_subtracted_%s_ak%s_%s_%d.pdf",algo,jet_type,date.GetDate()),"RECREATE");

    }//eta bins loop for the canvas

  }// random cone plotting if statement only for Vs so far

 
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;


}





































