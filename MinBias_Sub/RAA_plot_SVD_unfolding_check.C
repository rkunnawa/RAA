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


#include "../Headers/plot.h"

static const int nbins_cent = 6;
static const Int_t nSVD = 19;
double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};
static const double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};


using namespace std;

void RAA_plot_SVD_unfolding_check(int radius = 3,
				  char * etaWidth = (char*)"20_eta_20",
				  Float_t etaBoundary = 2.0)
{

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // get the input SVD files and histograms
  Int_t nSVDIter[nSVD] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  TFile * fin[nSVD];
  TH1F * hSVD[nSVD][nbins_cent], * hData[nbins_cent];
  TH1F * hSVD_pp[nSVD], * hData_pp;
  TH1F * hSVD_MC[nSVD][nbins_cent], * hGen[nbins_cent];
  TH1F * hSVD_MC_pp[nSVD], * hGen_pp;
  TH1F * hMC_Closure[nSVD][nbins_cent];
  TH1F * hMC_Closure_pp[nSVD];
  TH1F * hclosure_gen[nSVD][nbins_cent];
  TH1F * hclosure_gen_pp[nSVD];

  TH1F * hRAA_SVD[nSVD][nbins_cent];
  TH1F * hRAA_meas[nbins_cent];
  
  for(int j = 0; j<nSVD; ++j){

    for(int i = 0; i<nbins_cent; ++i)
      fin[j] = TFile::Open(Form("svd_unfolding_matrix_param%d_%s_60GeVCut_R%d_20150527.root",nSVDIter[j], etaWidth,radius));
    
    
    cout<<"resolution parameter j = "<<nSVDIter[j]<<endl;
    
    hSVD_pp[j] = (TH1F*)fin[j]->Get("PP_SVD_unfolding");
    hSVD_pp[j]->Print("base");
    // hSVD_MC_pp[j] = (TH1F*)fin[j]->Get("PP_MC_SVD_unfolding");
    // hSVD_MC_pp[j]->Print("base");
    if(j==0){
      hData_pp = (TH1F*)fin[j]->Get(Form("hpp_anaBin_HLTComb_R%d_%s",radius,etaWidth));
      hData_pp->Print("base");
      hGen_pp = (TH1F*)fin[j]->Get(Form("hpp_anaBin_JetComb_gen_R%d_%s",radius,etaWidth));
      hGen_pp->Print("base");
      // hData_pp->Scale(5.3 * 1e3);
      
    }
    // hclosure_gen_pp[j] = (TH1F*)fin[j]->Get(Form("hgen_pp_ak%d_0_7",radius));
    // hclosure_gen_pp[j]->Print("base");
    // hMC_Closure_pp[j] = (TH1F*)hSVD_MC_pp[j]->Clone(Form("hSVD_Closure_PP_iter%d",j));
    // hMC_Closure_pp[j]->Divide(hclosure_gen_pp[j]);
    
    // hSVD_pp[j]->Scale(5.3 * 1e3);
    for(int l = 1; l<=hSVD_pp[j]->GetNbinsX(); ++l){
      hSVD_pp[j]->SetBinError(l,hData_pp->GetBinError(l));
      // hSVD_MC_pp[j]->SetBinError(l,hData_pp->GetBinError(l));
    }
    
    for(int i = 0; i<nbins_cent; ++i){
      cout<<"centrality bin = "<<i<<endl;

      if(j==0) {

	hData[i] = (TH1F*)fin[j]->Get(Form("PbPb_data_minbiasSub_cent%d",i));
	hData[i]->Print("base");
	hGen[i] = (TH1F*)fin[j]->Get(Form("hpbpb_anaBin_JetComb_gen_R%d_%s_cent%d", radius, etaWidth, i));
	hGen[i]->Print("base");
	// hData[i]->Scale(145.156 * 1e3);
	
	hRAA_meas[i] = (TH1F*)hData[i]->Clone(Form("Measured_RAA_cent%d",i));
	hRAA_meas[i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
	hRAA_meas[i]->Scale(1./145.156); // triggered value
	hRAA_meas[i]->Scale(1./(7.65*1e6));
	hRAA_meas[i]->Scale(64.*1e9/(ncoll[i]*1e3));
	hRAA_meas[i]->Scale(5.3*1e3);
	hRAA_meas[i]->Divide(hData_pp);

      }

      // hclosure_gen[j][i] = (TH1F*)fin[j]->Get(Form("hgen_pbpb_akPu%d_1_%d",radius,i));
      // hclosure_gen[j][i]->Print("base");
      hSVD[j][i] = (TH1F*)fin[j]->Get(Form("PbPb_SVD_unfolding_cent%d",i));
      hSVD[j][i]->Print("base");
      // hSVD_MC[j][i] = (TH1F*)fin[j]->Get(Form("PbPb_MC_SVD_unfolding_cent%d",i));
      // hSVD_MC[j][i]->Print("base");
      // hMC_Closure[j][i] = (TH1F*)hSVD_MC[j][i]->Clone(Form("hSVD_Closure_iter%d_cent%d",j,i));
      // hMC_Closure[j][i]->Divide(hclosure_gen[j][i]);

      // hSVD[j][i]->Scale(145.156 * 1e3);
      for(int l = 1; l<=hSVD[j][i]->GetNbinsX(); ++l){
	hSVD[j][i]->SetBinError(l,hData[i]->GetBinError(l));
	// hSVD_MC[j][i]->SetBinError(l,hData[i]->GetBinError(l));
	
      }

      hRAA_SVD[j][i] = (TH1F*)hSVD[j][i]->Clone(Form("SVD_RAA_param%d_cent%d",j,i));
      hRAA_SVD[j][i]->Scale(1./(0.025*(boundaries_cent[i+1] - boundaries_cent[i])));
      hRAA_SVD[j][i]->Scale(1./145.156); // triggered value
      hRAA_SVD[j][i]->Scale(1./(7.65*1e6));
      hRAA_SVD[j][i]->Scale(64.*1e9/(ncoll[i]*1e3));
      hRAA_SVD[j][i]->Scale(5.3*1e3);
      hRAA_SVD[j][i]->Divide(hSVD_pp[j]);
      
    }

    
  }

  // make the RAA plot including SVD and measured.
  TCanvas * cSVD_RAA = new TCanvas("cSVD_RAA","",1200,1000);
  makeMultiPanelCanvas(cSVD_RAA,3,2,0.0,0.0,0.2,0.15,0.07);

  TLegend *tRAA = myLegend(0.45,0.45,0.65,0.9);
  TLine *lineRAA = new TLine(65,1,299,1);
  lineRAA->SetLineStyle(2);
  lineRAA->SetLineWidth(2);

  for(int i = 0;i<nbins_cent;++i){

    cSVD_RAA->cd(nbins_cent-i);

    hRAA_meas[i]->SetMarkerColor(kBlack);
    hRAA_meas[i]->SetMarkerStyle(24);
    makeHistTitle(hRAA_meas[i],"","Jet p_{T} (GeV/c)","R_{AA}");
    hRAA_meas[i]->SetAxisRange(65,299,"X");
    hRAA_meas[i]->SetAxisRange(0,2,"Y");
    hRAA_meas[i]->Draw("p");
    
    if(i == 0) tRAA->AddEntry(hRAA_meas[i],"Measured","pl");
    
    for(int j = 0; j<nSVD; ++j){

      // plot the different regularization parameter RAA
      hRAA_SVD[j][i]->SetMarkerStyle(20+j);
      hRAA_SVD[j][i]->SetMarkerColor(j+2);
      hRAA_SVD[j][i]->Draw("same p");
      
      if(i == 0){
	tRAA->AddEntry(hRAA_SVD[j][i],Form("SVD, k = %d",nSVDIter[j]),"pl");
      }
      
    }

    lineRAA->Draw();
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.85,20);

  }
  
  tRAA->SetTextSize(0.04);
  cSVD_RAA->cd(1);
  tRAA->Draw();
  cSVD_RAA->cd(1);
  putCMSPrel();
  cSVD_RAA->cd(2);
  putPbPbLumi();
  cSVD_RAA->cd(3);
  putPPLumi();
  cSVD_RAA->cd(1);
  drawText(Form("Anti-k_{T} Pu R=0.%d PF Jets", radius),0.25,0.20,16);
  cSVD_RAA->cd(2);
  drawText(Form("Jet ID cut, |#eta|<%2.0f",etaBoundary),0.1,0.2,16);
  drawText("|vz|<15, HBHEfilter, pCES",0.1,0.1,16);
  cSVD_RAA->cd(3);
  drawText("Jet RAA dataset, trigger combined",0.1,0.2,16);
  drawText("Pile up rejection cut applied",0.1,0.1,16);

  cSVD_RAA->SaveAs(Form("Final_paper_plots_RAA_SVD__20__%d_R%d_%s_pawan_ntuple.pdf",date.GetDate(),radius,etaWidth),"RECREATE");
  
  
  // Start plotting the curves in 6 different centrality bins.
  TCanvas * cSVD_check = new TCanvas("cSVD_check","",1200,1000);
  makeMultiPanelCanvasWithGap(cSVD_check,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend * leg = myLegend(0.6,0.4,0.8,0.9);
  for(int i = 0; i<nbins_cent; ++i){

    cSVD_check->cd(nbins_cent-i);
    cSVD_check->cd(nbins_cent-i)->SetLogy();
    cSVD_check->cd(nbins_cent-i)->SetLogx();

    hData[i]->SetMarkerStyle(20);
    hData[i]->SetMarkerColor(kBlack);
    hData[i]->SetAxisRange(40,300,"X");
    hData[i]->SetAxisRange(1e-2,1e7,"Y");
    makeHistTitle(hData[i]," ", " ak R=0.3 Jet p_{T} (GeV/c) " ,"counts");
    hData[i]->Draw("p");
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);
    
    for(int j = 0; j<nSVD; ++j){

      hSVD[j][i]->SetMarkerStyle(24+j);
      hSVD[j][i]->SetMarkerColor(j+2);
      hSVD[j][i]->Draw("same p");

      if(i==0)leg->AddEntry(hSVD[j][i],Form("SVD k = %d",nSVDIter[j]),"pl");

    }

  }

  cSVD_check->cd(1);
  putCMSPrel();
  leg->AddEntry(hData[0],"Measured","pl");
  leg->SetTextSize(0.04);
  
  leg->Draw();

  cSVD_check->cd(2);
  drawText("SVD unfolding",0.2,0.2,14);

  cSVD_check->SaveAs(Form("SVD__20__unfolding_%d_iteration_check_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");

#if 0
  
  // plot the MC closure for PbPb

  TCanvas * cClosure_PbPb = new TCanvas("cClosure_PbPb","",1200,1000);
  makeMultiPanelCanvasWithGap(cClosure_PbPb,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
  TLegend * mcclo = myLegend(0.6,0.7,0.8,0.9);
  for(int i = 0; i<nbins_cent; ++i){

    cClosure_PbPb->cd(nbins_cent-i);
    for(int j = 0; j<nSVD; ++j){

      hMC_Closure[j][i]->SetMarkerStyle(24);
      hMC_Closure[j][i]->SetMarkerColor(j+1);
      hMC_Closure[j][i]->SetAxisRange(0,2,"Y");
      hMC_Closure[j][i]->SetAxisRange(30,300,"X");
      makeHistTitle(hMC_Closure[j][i]," ", " ak R=0.3 Jet p_{T} (GeV/c) " ,"Reco/Truth");

      if(j==0) hMC_Closure[j][i]->Draw("p");
      if(i==0) mcclo->AddEntry(hMC_Closure[j][i],Form("k = %d",nSVDIter[j]),"pl");
      hMC_Closure[j][i]->Draw("same p");
    }
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);
    
  }
  
  cClosure_PbPb->cd(1);
  putCMSPrel();
  mcclo->Draw();
  
  cClosure_PbPb->cd(2);
  drawText("SVD Unfolding",0.2,0.2,14);
  cClosure_PbPb->SaveAs(Form("SVD__5__unfolding_%d_mcclosure_same_check_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");
  
#endif
  
  // Do the spectra check for pp and PP mc closure.
  TCanvas * cSVD_check_pp = new TCanvas("cSVD_check_pp","",800,700);
  TLegend * leg_pp = myLegend(0.6,0.4,0.8,0.9);
  cSVD_check_pp->SetLogy();
  cSVD_check_pp->SetLogx();
  hData_pp->SetMarkerStyle(20);
  hData_pp->SetMarkerColor(kBlack);
  hData_pp->SetAxisRange(60,299,"X");
  hData_pp->SetAxisRange(1e-2,1e7,"Y");
  makeHistTitle(hData_pp," ", Form(" ak R=0.%d Jet p_{T} (GeV/c) ",radius) ,"counts");
  hData_pp->Draw("p");
  for(int j = 0; j<nSVD; ++j){

    hSVD_pp[j]->SetMarkerStyle(24+j);
    hSVD_pp[j]->SetMarkerColor(j+2);
    hSVD_pp[j]->Draw("same p");

    leg_pp->AddEntry(hSVD_pp[j],Form("SVD k = %d",nSVDIter[j]),"pl");

  }
  putCMSPrel();
  leg_pp->AddEntry(hData_pp,"Measured","pl");
  leg_pp->SetTextSize(0.04);
  leg_pp->Draw();

  drawText("SVD unfolding pp",0.2,0.2,14);

  cSVD_check_pp->SaveAs(Form("SVD__20__unfolding_%d_iteration_check_pp_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");

#if 0
  // plot the MC closure for pp
  TCanvas * cClosure_pp = new TCanvas("cClosure_pp","",1200,1000);
  TLegend * mcclo_pp = myLegend(0.6,0.4,0.8,0.9);

  for(int j = 0; j<nSVD; ++j){
    
    hMC_Closure_pp[j]->SetMarkerStyle(24);
    hMC_Closure_pp[j]->SetMarkerColor(j+1);
    hMC_Closure_pp[j]->SetAxisRange(0.6,1.4,"Y");
    hMC_Closure_pp[j]->SetAxisRange(60,299,"X");
    makeHistTitle(hMC_Closure_pp[j]," ", Form(" ak R=0.%d Jet p_{T} (GeV/c) ",radius) ,"Reco/Truth");
    
    if(j==0) hMC_Closure_pp[j]->Draw("p");
    mcclo_pp->AddEntry(hMC_Closure_pp[j],Form("SVD k = %d",nSVDIter[j]),"pl");
    hMC_Closure_pp[j]->Draw("same p");
  }
  
  putCMSPrel();
  mcclo_pp->Draw();
  
  drawText("SVD Unfolding pp",0.2,0.2,14);

  cClosure_pp->SaveAs(Form("SVD__5__unfolding_%d_mcclosure_same_check_pp_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");
#endif  
  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
  
