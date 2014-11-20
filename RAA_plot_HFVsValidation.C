// Raghav Kunnawalkam Elayavalli
// Nov 19 2014
// Rutgers
// raghav.k.e at CERN dot CH

//
// Plotting macro for HF Vs algorithm validation. 
// 


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


static const int nbins_eta = 15;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0}, {-2.0,+2.0}, {-3.0,+3.0},
  {-3.0,-2.5}, {-2.5,-2.0}, {-2.0,-1.5}, 
  {-1.5,-1.0}, {-1.0,-0.5}, {-0.5,0}, {0,+0.5}, 
  {+0.5,+1.0}, {+1.0,+1.5}, {+1.5,+2.0}, 
  {+2.0,+2.5}, {+2.5,+3.0}
};

static const double delta_eta[nbins_eta] = {
  2.0, 4.0, 6.0, 
  0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 
  0.5, 0.5
};

static const char etaWidth [nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20","n30_eta_p30",
  "n30_eta_n25","n25_eta_n20","n20_eta_n15",
  "n15_eta_n10","n10_eta_n05","n05_eta_0","0_eta_p05",
  "p05_eta_p10","p10_eta_p15","p15_eta_p20",
  "p20_eta_p25","p25_eta_p30"
};


/*
static const int nbins_eta = 2;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0},
  {-2.0,+2.0}
};

static const double delta_eta[nbins_eta] = {
  2.0,4.0
};

static const char etaWidth[nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20"
};
*/

//static const int no_radius = 2;//testing purposes 
//static const int list_radius[no_radius] = {3,4};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
static const int no_radius = 7; 
static const int list_radius[no_radius] = {1,2,3,4,5,6,7};
/*
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
*/

using namespace std;


void RAA_plot_HFVsValidation(int radius = 3, char *algo = "Vs", char *jet_type="PF"){


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
  
  //get the histograms from the MC file. 
  TFile *fMCin = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/PbPb_pp_mc_nocut_scaletest_ak%s%s_20141120.root",algo,jet_type));
  
  // UseFull histograms: Event Plane from HF for the official versus HF/Vs algorithm calculation of the Event Plane.
  // first [3] array elements - Psi_2, Psi_3, Psi_4  only in the HF for tonight. 
  // [3] - total, p - positive and n - negative in eta; based on what we have from the https://github.com/CmsHI/cmssw/blob/forest_CMSSW_5_3_20/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h#L31 
  TH1F *hEP_HF_Official[3][3][nbins_cent];
  TH1F *hEP_HF_Vs[3][3][nbins_cent];
  TH2F *hEP_HF[3][3][nbins_cent];
  TH1F *hNJetsvsSumpT[nbins_cent];
  TH1F *hSumpT[nbins_cent];
  TH2F *hSumpTvsHF[15];

  for(int i = 0;i<nbins_cent;i++){
    hNJetsvsSumpT[i] = (TH1F*)fMCin->Get(Form("hNJetsvsSumpT_cent%d",i));
    hSumpT[i] = (TH1F*)fMCin->Get(Form("hSumpT_cent%d",i));
    for(int z = 0;z<3;z++){
      for(int x = 0;x<3;x++){
	hEP_HF_Official[z][x][i] = (TH1F*)fMCin->Get(Form("hEP_HF_Official_Psi%d_%d_cent%d",z,x,i));
	hEP_HF_Vs[z][x][i] = (TH1F*)fMCin->Get(Form("hEP_HF_Vs_Psi%d_%d_cent%d",z,x,i));
	hEP_HF[z][x][i] = (TH2F*)fMCin->Get(Form("hEP_HF_Psi%d_%d_cent%d",z,x,i));
      }
    }
  }

  for(int a = 0;a<15;a++){
    hSumpTvsHF[a] = (TH2F*)fMCin->Get(Form("hsumpT_eta%d_vsHF",a));
  }
  
  // time to make plots! 
  
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Plot showing the comparison bewteen official and Vs calculation of the event plane. 
  // this is going to be a 3 by 2 set of 3 plots, each for v2, v3 and v4. 
  /*
  TCanvas *cEP_Psi[3][3]; 
  
  for(int z = 0;z<3;z++){
    for(int x = 0;x<3;x++){

      cEP_Psi[z][x] = new TCanvas(Form("cEP_Psi_%d_%d",z,x),"",1000,800);
      makeMultiPanelCanvas(cEP_Psi[z][x],3,2,0.0,0.0,0.2,0.15,0.07); 
      
      for(int i = 0;i<nbins_cent;i++){

	cEP_Psi[z][x]->cd(nbins_cent-i);
	cEP_Psi[z][x]->cd(nbins_cent-i)->SetLogy();
	makeHistTitle(hEP_HF_Vs[z][x][i]," ",Form("#Psi_{%d}",z+2),"Event Fraction");

	hEP_HF_Vs[z][x][i]->Rebin(30);
	divideBinWidth(hEP_HF_Vs[z][x][i]);
	hEP_HF_Vs[z][x][i]->Scale(1./hEP_HF_Vs[z][x][i]->Integral());
	hEP_HF_Vs[z][x][i]->SetMarkerStyle(24);
	hEP_HF_Vs[z][x][i]->SetMarkerColor(kRed);
	hEP_HF_Vs[z][x][i]->Draw();

	hEP_HF_Official[z][x][i]->Rebin(30);
	divideBinWidth(hEP_HF_Official[z][x][i]);
	hEP_HF_Official[z][x][i]->Scale(1./hEP_HF_Official[z][x][i]->Integral());
	hEP_HF_Official[z][x][i]->SetMarkerStyle(24);
	hEP_HF_Official[z][x][i]->SetMarkerColor(kBlack);
	hEP_HF_Official[z][x][i]->Draw("same");

	drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.7,20);

      }
   
      cEP_Psi[z][x]->cd(1);
      TLegend *PbPb_Psi = myLegend(0.3,0.15,0.5,0.35);
      PbPb_Psi->AddEntry(hEP_HF_Official[0][0][0],"Official EP","pl");
      PbPb_Psi->AddEntry(hEP_HF_Vs[0][0][0],"Vs EP","pl");
      PbPb_Psi->SetTextSize(0.04);
      PbPb_Psi->Draw();
      
      putCMSSim();
      cEP_Psi[z][x]->cd(2);
      if(x==0)drawText("Event Plane Both HF",0.2,0.2,16);
      if(x==1)drawText("Event Plane forward HF",0.2,0.2,16);
      if(x==2)drawText("Event Plane backward HF",0.2,0.2,16);
      cEP_Psi[z][x]->cd(3);
      drawText("|vz|<15, pCES, HBHE",0.23,0.85,16);
      cEP_Psi[z][x]->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_HFVs_validation_scaletest_ak%s%s_EventPlane_Psi%d_%d_%d.pdf",algo,jet_type,z,x,date.GetDate()),"RECREATE");
    }
  }
  */

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Plot showing the comparison bewteen official and Vs calculation of the event plane. 
  // this is going to be a 6 by 3 set of 3 plots, each for v2, v3 and v4. 
  /*
  TCanvas *cEP_Psi2D[3][3]; 

  for(int z = 0;z<3;z++){
    for(int x = 0;x<3;x++){

      cEP_Psi2D[z][x] = new TCanvas(Form("cEP_Psi2D_%d_%d",z,x),"",1000,800);
      makeMultiPanelCanvas(cEP_Psi2D[z][x],3,2,0.0,0.0,0.2,0.15,0.07); 
      
      for(int i = 0;i<nbins_cent;i++){

	cEP_Psi2D[z][x]->cd(nbins_cent-i);
	cEP_Psi2D[z][x]->cd(nbins_cent-i)->SetLogz();
	makeHistTitle(hEP_HF[z][x][i]," ",Form("Official #Psi_{%d}",z+2),Form("Vs #Psi_{%d}",z+2));

	//hEP_HF[z][x][i]->Rebin2D(30);
	//divideBinWidth(hEP_HF_Vs[z][x][i]);
	hEP_HF[z][x][i]->Scale(1./hEP_HF[z][x][i]->Integral());
	hEP_HF[z][x][i]->Draw("colz");

	drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.7,20);

      }
   
      cEP_Psi2D[z][x]->cd(1);

      putCMSSim();
      cEP_Psi2D[z][x]->cd(2);
      if(x==0)drawText("Event Plane Both HF",0.2,0.2,16);
      if(x==1)drawText("Event Plane forward HF",0.2,0.2,16);
      if(x==2)drawText("Event Plane backward HF",0.2,0.2,16);
      cEP_Psi2D[z][x]->cd(3);
      drawText("|vz|<15, pCES, HBHE",0.23,0.85,16);
      cEP_Psi2D[z][x]->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_HFVs_validation_scaletest_ak%s%s_EventPlane2D_Psi%d_%d_%d.pdf",algo,jet_type,z,x,date.GetDate()),"RECREATE");
    }
  }
  */
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Plot showing the NJets vs SumpT for each centralty class. 
  TCanvas *cSumpT = new TCanvas("cSumpT","",1000,800);
  makeMultiPanelCanvas(cSumpT,3,2,0.0,0.0,0.2,0.15,0.07); 

  for(int i = 0;i<nbins_cent;i++){

    cSumpT->cd(nbins_cent-i);
    cSumpT->cd(nbins_cent-i)->SetLogy();
    cSumpT->cd(nbins_cent-i)->SetLogz();

    makeHistTitle(hNJetsvsSumpT[i],"","No Jets (pT>50GeV)","SumpT from all eta bins");
    hNJetsvsSumpT[i]->SetAxisRange(1,1e6,"Y");
    hNJetsvsSumpT[i]->Draw("colz");
    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.75,0.7,20);

  }
  cSumpT->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Plots/PbPb_HFVs_validation_scaletest_ak%s%s_sumpT_vs_NJets_%d.pdf",algo,jet_type,date.GetDate()),"RECREATE");

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // 2d plot showing the correlation of the two event planes on an event by event level.


  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // plot showing the SumpT vs HF energy deposited, in all the eta bins. 
  //TH2F *hSumpTTotvsHF = (TH2F*)hSumpTvsHF[0]->Clone("hSumpTTotvsHF");
  

  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
