// Raghav Kunnawalkam Elayavalli
// April 22 2015
// CERN




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

static const int nbins_cent = 6;
static const Int_t nSVD =3;
double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};

using namespace std;

void RAA_plot_SVD_unfolding_check(int radius = 3){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // get the input SVD files and histograms
  //Int_t nSVDIter[nSVD] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  Int_t nSVDIter[nSVD] = {4,5,6};
  TFile * fin[nSVD];
  TFile * fPbPb_MCin;
  TFile * fPP_MCin;
  fPbPb_MCin = TFile::Open(Form("../../Output/PbPb_MC_closure_histogram_deltaR_0p2_akPu%d_20150423.root",radius));
  fPP_MCin = TFile::Open(Form("../../Output/pp_MC_closure_histogram_deltaR_0p2_ak%d_20150423.root",radius));
  TH1F * hSVD[nSVD][nbins_cent], * hData[nbins_cent];
  TH1F * hSVD_pp[nSVD], * hData_pp;
  TH1F * hSVD_MC[nSVD][nbins_cent], * hGen[nbins_cent];
  TH1F * hSVD_MC_pp[nSVD], * hGen_pp;
  TH2F * hmatrix[nbins_cent];
  TH1F * hMC_Closure[nSVD][nbins_cent];
  TH1F * hMC_Closure_pp[nSVD];

  TFile * fPbPb_in = TFile::Open(Form("../../Output/PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_rebinned_eMaxSumcand_A_R0p%d.root",radius));
  TFile * fPP_in = TFile::Open(Form("../../Output/Pp_CutEfficiency_YetkinCuts_exclusionhighertriggers_rebinned_eMaxSumcand_A_R0p%d.root",radius));
  
  for(int j = 0; j<nSVD; ++j){
    cout<<"iteration j = "<<nSVDIter[j]<<endl;
    for(int i = 0; i<nbins_cent; ++i){
      cout<<"centrality bin = "<<i<<endl;
      fin[j] = TFile::Open(Form("../../Output/svd_unfolding_sameside_matrix_param%d_R%d_20150424.root",nSVDIter[j],radius));
      if(j==0) hmatrix[i] = (TH2F*)fPbPb_in->Get(Form("hpbpb_matrix_R%d_n20_eta_p20_cent%d",radius,i));
      if(j==0) hData[i] = (TH1F*)fin[j]->Get(Form("hpbpb_HLTComb_R%d_n20_eta_p20_cent%d",radius,i));
      if(j==0) hGen[i] = (TH1F*)hmatrix[i]->ProjectionX();
      hSVD[j][i] = (TH1F*)fin[j]->Get(Form("PbPb_SVD_unfolding_cent%d",i));
      hSVD_MC[j][i] = (TH1F*)fin[j]->Get(Form("PbPb_MC_SVD_unfolding_cent%d",i));
      hMC_Closure[j][i] = (TH1F*)hSVD_MC[j][i]->Clone(Form("hSVD_Closure_iter%d_cent%d",j,i));
      hMC_Closure[j][i]->Divide(hGen[i]);

    }

    hSVD_pp[j] = (TH1F*)fin[j]->Get("PP_SVD_unfolding_cent%d");
    hSVD_MC_pp[j] = (TH1F*)fin[j]->Get("PP_MC_SVD_unfolding_cent%d");
    if(j==0)hData_pp = (TH1F*)fin[j]->Get(Form("hpp_HLTComb_R%d_n20_eta_p20",radius));
    if(j==0)hGen_pp = (TH1F*)fPP_MCin->Get(Form("hpp_mcclosure_gen_JetComb_R%d_n20_eta_p20",radius));
    hMC_Closure_pp[j] = (TH1F*)hSVD_MC_pp[j]->Clone(Form("hSVD_Closure_PP_iter%d",j));
    hMC_Closure_pp[j]->Divide(hGen_pp);
    
  }

  // Start plotting the curves in 6 different centrality bins.
  TCanvas * cSVD_check = new TCanvas("cSVD_check","",1200,1000);
  makeMultiPanelCanvasWithGap(cSVD_check,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend * leg = myLegend(0.5,0.55,0.7,0.95);
  for(int i = 0; i<nbins_cent; ++i){

    cSVD_check->cd(nbins_cent-i);
    cSVD_check->cd(nbins_cent-i)->SetLogy();
    cSVD_check->cd(nbins_cent-i)->SetLogx();

    hData[i]->SetMarkerStyle(20);
    hData[i]->SetMarkerColor(kBlack);
    hData[i]->SetAxisRange(40,300,"X");
    hData[i]->SetAxisRange(1e-2,1e7,"Y");
    makeHistTitle(hData[i]," ", " ak R=0.3 Jet p_{T} (GeV/c) " ,"counts");
    hData[i]->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);
    for(int j = 0; j<nSVD; ++j){

      hSVD[j][i]->SetMarkerStyle(24+j);
      hSVD[j][i]->SetMarkerColor(j+2);
      hSVD[j][i]->Draw("same");

      if(i==0)leg->AddEntry(hSVD[j][i],Form("k = %d",nSVDIter[j]),"pl");

    }

  }

  cSVD_check->cd(1);
  putCMSPrel();
  leg->Draw();

  cSVD_check->cd(2);
  drawText("SVD unfolding",0.2,0.2,14);

  cSVD_check->SaveAs(Form("../../Plots/SVD_unfolding_%d_iteration_check_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");
  
  // plot the MC closure for PbPb

  TCanvas * cClosure_PbPb = new TCanvas("cClosure_PbPb","",1200,1000);
  makeMultiPanelCanvasWithGap(cClosure_PbPb,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  TLegend * mcclo = myLegend(0.4,0.1,0.7,0.9);
  for(int i = 0; i<nbins_cent; ++i){

    cClosure_PbPb->cd(nbins_cent-i);
    for(int j = 0; j<nSVD; ++j){

      hMC_Closure[j][i]->SetMarkerStyle(24);
      hMC_Closure[j][i]->SetMarkerColor(j+1);
      hMC_Closure[j][i]->SetAxisRange(0,2,"Y");
      hMC_Closure[j][i]->SetAxisRange(30,300,"X");
      makeHistTitle(hMC_Closure[j][i]," ", " ak R=0.3 Jet p_{T} (GeV/c) " ,"Reco/Truth");

      if(j==0) hMC_Closure[j][i]->Draw();
      if(i==0) mcclo->AddEntry(hMC_Closure[j][i],Form("k = %d",nSVDIter[j]),"pl");
      hMC_Closure[j][i]->Draw("same");
    }
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);
    
  }
  cClosure_PbPb->cd(1);
  putCMSPrel();
  mcclo->Draw();
  
  cClosure_PbPb->cd(2);
  drawText("SVD Unfolding",0.2,0.2,14);

  cClosure_PbPb->SaveAs(Form("../../Plots/SVD_unfolding_%d_mcclosure_same_check_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");


  // Do the spectra check for pp and PP mc closure.
  // Start plotting the curves in 6 different centrality bins.
  TCanvas * cSVD_check_pp = new TCanvas("cSVD_check_pp","",1200,1000);
  TLegend * leg_pp = myLegend(0.5,0.55,0.7,0.95);
  cSVD_check_pp->SetLogy();
  cSVD_check_pp->SetLogx();
  hData_pp->SetMarkerStyle(20);
  hData_pp->SetMarkerColor(kBlack);
  hData_pp->SetAxisRange(40,400,"X");
  hData_pp->SetAxisRange(1e-2,1e7,"Y");
  makeHistTitle(hData_pp," ", " ak R=0.3 Jet p_{T} (GeV/c) " ,"counts");
  hData_pp->Draw();
  for(int j = 0; j<nSVD; ++j){

    hSVD_pp[j]->SetMarkerStyle(24+j);
    hSVD_pp[j]->SetMarkerColor(j+2);
    hSVD_pp[j]->Draw("same");

    leg->AddEntry(hSVD_pp[j],Form("k = %d",nSVDIter[j]),"pl");

  }
  putCMSPrel();
  leg_pp->Draw();

  drawText("SVD unfolding pp",0.2,0.2,14);

  cSVD_check_pp->SaveAs(Form("../../Plots/SVD_unfolding_%d_iteration_check_pp_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");
  
  // plot the MC closure for PbPb

  TCanvas * cClosure_pp = new TCanvas("cClosure_pp","",1200,1000);
  TLegend * mcclo_pp = myLegend(0.4,0.1,0.7,0.9);

  for(int j = 0; j<nSVD; ++j){
    
    hMC_Closure_pp[j]->SetMarkerStyle(24);
    hMC_Closure_pp[j]->SetMarkerColor(j+1);
    hMC_Closure_pp[j]->SetAxisRange(0,2,"Y");
    hMC_Closure_pp[j]->SetAxisRange(30,300,"X");
    makeHistTitle(hMC_Closure_pp[j]," ", " ak R=0.3 Jet p_{T} (GeV/c) " ,"Reco/Truth");
    
    if(j==0) hMC_Closure_pp[j]->Draw();
    mcclo_pp->AddEntry(hMC_Closure_pp[j],Form("k = %d",nSVDIter[j]),"pl");
    hMC_Closure_pp[j]->Draw("same");
  }
  
  putCMSPrel();
  mcclo_pp->Draw();
  
  drawText("SVD Unfolding pp",0.2,0.2,14);

  cClosure_pp->SaveAs(Form("../../Plots/SVD_unfolding_%d_mcclosure_same_check_pp_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");

  
  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
  
