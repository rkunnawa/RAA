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
static const Int_t nSVD =1;
double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};

using namespace std;

void RAA_plot_SVD_unfolding_check(int radius = 3){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // get the input SVD files and histograms
  Int_t nSVDIter[nSVD] = {5};
  TFile * fin[nSVD];
  TH1F * hSVD[nSVD][nbins_cent], * hData[nbins_cent];

  for(int j = 0; j<nSVD; ++j){

    for(int i = 0; i<nbins_cent; ++i){

      fin[j] = TFile::Open(Form("../../Output/svd_unfolding_noHLT_param%d_test_R%d_20150422.root",nSVDIter[j],radius));

      if(j==0) hData[i] = (TH1F*)fin[j]->Get(Form("hpbpb_HLTComb_R%d_n20_eta_p20_cent%d",radius,i));
      hSVD[j][i] = (TH1F*)fin[j]->Get(Form("PbPb_SVD_unfolding_cent%d",i));

    }
    
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
  drawText("",0.2,0.2,14);

  cSVD_check->SaveAs(Form("../../Plots/SVD_unfolding_%d_iteration_check_R%d_%d.pdf",nSVD,radius,date.GetDate()),"RECREATE");
  
  //
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
  
