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



static const int nbins_cent = 7;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36,40};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 

int findBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins               
  if(bin<=10)ibin=0; //! 0-5%
  else if(bin>10  && bin<=20 )ibin=1; //! 5-10%
  else if(bin>20  && bin<=60 )ibin=2;  //! 10-30%
  else if(bin>60  && bin<=100)ibin=3;  //! 30-50%
  else if(bin>100 && bin<=140)ibin=4;  //! 50-70%
  else if(bin>140 && bin<=180)ibin=5;  //! 70-90%
  else if(bin>180 && bin<=200)ibin=6;  //! 90-100%
  return ibin;
}

static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};



using namespace std;

void RAA_plot_YenJie_CutfromMinBias(char* etaWidth = (char*)"10_eta_18", Int_t radius = 3)
{

  TStopwatch timer;
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  cout<<"etaWidth = "<<etaWidth<<endl;
  bool isAnalysisBin = false; 

  TDatime date;
  
  TFile * fMinBias, * fJetTrig;

  fMinBias = TFile::Open(Form("../../Output/Pawan_ntuple_PbPb_MinBiasData_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius));
  fJetTrig = TFile::Open(Form("../../Output/Pawan_ntuple_PbPb_data_MC_subid0_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius));

  // get the histograms
  TH1F * hJetTrigComb[nbins_cent], * hMinBias[nbins_cent];
  TH1F * hSubtracted[nbins_cent];
  Float_t bin_no = 1;

  for(int i = 0; i<nbins_cent; ++i){

    hJetTrigComb[i] = (TH1F*)fJetTrig->Get(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i));
    hMinBias[i]     = (TH1F*)fMinBias->Get(Form("hpbpb_noTrg_R%d_%s_cent%d",radius,etaWidth,i));

    // rebin the histograms (maybe, commented out for now, since we are doing the analysis in fine bins) and normalize them here to the bin width of the first pT bin for the analysis: around 15 GeV

    // the input to the datadrivenunfolding takes in analysis bins input, so for that you need to rebin these histograms. but then for RAA_analyze.C the input is fine bins, so you need to rerun this without the rebinning. 
    if(isAnalysisBin){
      hJetTrigComb[i] = (TH1F*)hJetTrigComb[i]->Rebin(nbins_cent, Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i), boundaries_cent);
      divideBinWidth(hJetTrigComb[i]);
      hMinBias[i] = (TH1F*)hMinBias[i]->Rebin(nbins_cent, Form("hpbpb_noTrg_R%d_%s_cent%d",radius,etaWidth,i), boundaries_cent);
      divideBinWidth(hMinBias[i]);
    }
    // find the bin number of 15 GeV
    bin_no = hJetTrigComb[i]->FindBin(15);
    // set all the entries before bin no 15 to be zero.
    for(int j = 1; j<bin_no; ++j){
      hJetTrigComb[i]->SetBinContent(j,0);
      hJetTrigComb[i]->SetBinError(j,0);
      hMinBias[i]->SetBinContent(j,0);
      hMinBias[i]->SetBinError(j,0);
    }

    // normalize the histogram by the jet entries in bin of pT = 15.
    hJetTrigComb[i]->Scale(1./hJetTrigComb[i]->GetBinContent(bin_no));
    hMinBias[i]->Scale(1./hMinBias[i]->GetBinContent(bin_no));
    
    hSubtracted[i]  = (TH1F*)hJetTrigComb[i]->Clone(Form("hFakeMinBias_SubtractedFrom_Data_R%d_%s_cent%d",radius, etaWidth, i));
    hSubtracted[i]->Add(hMinBias[i], -1);
    
  }

  // make the plots: just a 3x2 panel plot showing all the above three curves and write to output.
  
  TCanvas * cFakeSub = new TCanvas("cFakeSub","",1000,800);
  makeMultiPanelCanvasWithGap(cFakeSub,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend * lSub = myLegend(0.5,0.5,0.7,0.7);
  for(int i = 0; i<nbins_cent; ++i){

    cFakeSub->cd(nbins_cent-i);
    cFakeSub->cd(nbins_cent-i)->SetLogy();

    hJetTrigComb[i]->SetMarkerStyle(24);    
    hJetTrigComb[i]->SetMarkerColor(kBlack);
    hJetTrigComb[i]->SetAxisRange(13,299,"X");
    makeHistTitle(hJetTrigComb[i]," ","Jet p_{T} (GeV/c)","Counts normalized to ptbin = 15 GeV");
    hJetTrigComb[i]->Draw();

    hMinBias[i]->SetMarkerStyle(24);    
    hMinBias[i]->SetMarkerColor(kRed);
    hMinBias[i]->Draw("same");

    hSubtracted[i]->SetMarkerStyle(20);    
    hSubtracted[i]->SetMarkerColor(kBlack);
    hSubtracted[i]->Draw("same");

    drawText(Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.8,0.85,20);    

  }

  cFakeSub->cd(1);
  putCMSPrel();
  
  lSub->AddEntry(hJetTrigComb[0],"Jet Triggered Data","pl");
  lSub->AddEntry(hMinBias[0],"MinBias data","pl");
  lSub->AddEntry(hSubtracted[0],"Jet Triggered - MinBias","pl");
  lSub->SetTextSize(0.04);
  lSub->Draw();
  
  cFakeSub->SaveAs(Form("../../Plots/PbPb_FakeMinBias_SubtractedFrom_Data_R%d_%s_%d.pdf",radius, etaWidth, date.GetDate()),"RECREATE");
  
  // create output root file to write the output histogram. 
  TFile fout(Form("../../Output/PbPb_FakeMinBias_SubtractedFrom_Data_R%d_%s_%d.root",radius, etaWidth, date.GetDate()),"RECREATE");
  fout.cd();

  for(int i = 0; i<nbins_cent; ++i){

    hSubtracted[i]->Write();
    hJetTrigComb[i]->Write();
    hMinBias[i]->Write();

  }

  fout.Close();

  // macro end. 
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}
