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

using namespace std;

void RAA_plot_matched_unmatched_eMaxSumcand(){

  Int_t radius = 3;
  TFile * fData, * fMC; 
  /*
  if(radius == 2) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu2_20150331.root");
  if(radius == 3) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150331.root");
  //if(radius == 3) fData = TFile::Open("/export/d00/scratch/pawan/condorfiles/pbpb/ntuples/Merged_JetRAA_akPu3_PbPb_Data.root");
  if(radius == 4) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu4_20150331.root");

  if(radius == 2) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu2_20150331.root");
  if(radius == 3) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150331.root");
  if(radius == 4) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu4_20150331.root");
  */

  if(radius == 2) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/pp_Data_calo_pf_jet_correlation_deltaR_0p2_ak2_20150331.root");
  if(radius == 3) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/pp_Data_calo_pf_jet_correlation_deltaR_0p2_ak3_20150331.root");
  if(radius == 4) fData = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/pp_Data_calo_pf_jet_correlation_deltaR_0p2_ak4_20150331.root");

  if(radius == 2) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/pp_MC_calo_pf_jet_correlation_deltaR_0p2_ak2_20150331.root");
  if(radius == 3) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/pp_MC_calo_pf_jet_correlation_deltaR_0p2_ak3_20150331.root");
  if(radius == 4) fMC = TFile::Open("/export/d00/scratch/rkunnawa/rootfiles/pp_MC_calo_pf_jet_correlation_deltaR_0p2_ak4_20150331.root");

  TTree * Data_matched = (TTree*)fData->Get("matchedJets");
  TTree * Data_unmatched = (TTree*)fData->Get("unmatchedPFJets");

  TTree * MC_matched = (TTree*)fMC->Get("matchedJets");
  TTree * MC_unmatched = (TTree*)fMC->Get("unmatchedPFJets");

  const int ptSelection = 17;
  const int ptBoundary[ptSelection+1] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
  
  TH2F* hData_eMaxSumcand_calopfpt_ptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_calopfpt_ptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_calopfpt_refptselection[ptSelection];
  TCanvas * ceMaxSumcand [ptSelection];

  TLine * CutA_1 = new TLine(0.5, 0.05, 0.85, 0.9);
  CutA_1->SetLineStyle(2);
  CutA_1->SetLineWidth(2);
  CutA_1->SetLineColor(kRed);
  TLine * CutA_2 = new TLine(0.85, 0.9, 0.85, 5);
  CutA_2->SetLineStyle(2);
  CutA_2->SetLineWidth(2);
  CutA_2->SetLineColor(kRed);
  TLine * CutA_3 = new TLine(0, 0.05, 0.5, 0.05);
  CutA_3->SetLineStyle(2);
  CutA_3->SetLineWidth(2);
  CutA_3->SetLineColor(kRed);

  TH2F* hData_unmatch_eMaxSumcand_calopfpt_ptselection[ptSelection];
  TH2F* hMC_unmatch_eMaxSumcand_calopfpt_ptselection[ptSelection];
  TH2F* hMC_unmatch_eMaxSumcand_calopfpt_refptselection[ptSelection];
 
  for(int a = 0;a<ptSelection; ++a) {

    hData_eMaxSumcand_calopfpt_ptselection[a] = new TH2F(Form("hData_eMaxSumcand_calopfpt_ptselection_%d",a),"",100, 0, 2.5, 100, 0, 5);
    hMC_eMaxSumcand_calopfpt_ptselection[a] = new TH2F(Form("hMC_eMaxSumcand_calopfpt_ptselection_%d",a),"",100, 0, 2.5, 100, 0, 5);
    hMC_eMaxSumcand_calopfpt_refptselection[a] = new TH2F(Form("hMC_eMaxSumcand_calopfpt_refptselection_%d",a),"",100, 0, 2.5, 100, 0, 5);
    
    hData_unmatch_eMaxSumcand_calopfpt_ptselection[a] = new TH2F(Form("hData_unmatch_eMaxSumcand_calopfpt_ptselection_%d",a),"",100, 0, 2.5, 100, 0, 5);
    hMC_unmatch_eMaxSumcand_calopfpt_ptselection[a] = new TH2F(Form("hMC_unmatch_eMaxSumcand_calopfpt_ptselection_%d",a),"",100, 0, 2.5, 100, 0, 5);
    hMC_unmatch_eMaxSumcand_calopfpt_refptselection[a] = new TH2F(Form("hMC_unmatch_eMaxSumcand_calopfpt_refptselection_%d",a),"",100, 0, 2.5, 100, 0, 5);
    
    Data_matched->Draw(Form("eMax/(chSum+phSum+neSum+muSum):calopt/pfpt>>hData_eMaxSumcand_calopfpt_ptselection_%d",a),Form("pfpt > %d && pfpt < %d && jet40", ptBoundary[a], ptBoundary[a+1]),"goff");
    MC_matched->Draw(Form("eMax/(chSum+phSum+neSum+muSum):calopt/pfpt>>hMC_eMaxSumcand_calopfpt_ptselection_%d",a),Form("pfpt > %d && pfpt < %d && jet40", ptBoundary[a], ptBoundary[a+1]),"goff");
    MC_matched->Draw(Form("eMax/(chSum+phSum+neSum+muSum):calopt/pfpt>>hMC_eMaxSumcand_calopfpt_refptselection_%d",a),Form("pfpt > %d && pfpt < %d && pfrefpt > %d && pfrefpt < %d && jet40",  ptBoundary[a], ptBoundary[a+1],  ptBoundary[a], ptBoundary[a+1]),"goff");

    Data_unmatched->Draw(Form("eMax/(chSum+phSum+neSum+muSum):(pfpt-pfpt)>>hData_unmatch_eMaxSumcand_calopfpt_ptselection_%d",a),Form("pfpt > %d && pfpt < %d && jet40", ptBoundary[a], ptBoundary[a+1]),"goff");
    MC_unmatched->Draw(Form("eMax/(chSum+phSum+neSum+muSum):(pfpt-pfpt)>>hMC_unmatch_eMaxSumcand_calopfpt_ptselection_%d",a),Form("pfpt > %d && pfpt < %d && jet40", ptBoundary[a], ptBoundary[a+1]),"goff");
    MC_unmatched->Draw(Form("eMax/(chSum+phSum+neSum+muSum):(pfpt-pfpt)>>hMC_unmatch_eMaxSumcand_calopfpt_refptselection_%d",a),Form("pfpt > %d && pfpt < %d && pfrefpt > %d && pfrefpt < %d && jet40",  ptBoundary[a], ptBoundary[a+1],  ptBoundary[a], ptBoundary[a+1]),"goff");

    hData_eMaxSumcand_calopfpt_ptselection[a]->Add(hData_unmatch_eMaxSumcand_calopfpt_ptselection[a]);
    hMC_eMaxSumcand_calopfpt_ptselection[a]->Add(hMC_unmatch_eMaxSumcand_calopfpt_ptselection[a]);
    hMC_eMaxSumcand_calopfpt_refptselection[a]->Add(hMC_unmatch_eMaxSumcand_calopfpt_refptselection[a]);
    
    ceMaxSumcand[a] = new TCanvas(Form("ceMaxSumcand_%d",a),"",1000,800);
    ceMaxSumcand[a]->Divide(3,1);
    ceMaxSumcand[a]->cd(1);
    ceMaxSumcand[a]->cd(1)->SetLogz();
    hData_eMaxSumcand_calopfpt_ptselection[a]->SetXTitle("calopt/pfpt");
    hData_eMaxSumcand_calopfpt_ptselection[a]->SetYTitle("eMax/(chSum+phSum+neSum+muSum)");
    hData_eMaxSumcand_calopfpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < pfpt < %d", ptBoundary[a], ptBoundary[a+1]),0.2,0.7,14);
    drawText("Data - Jet40",0.25,0.8,14);
    CutA_1->Draw();
    CutA_2->Draw();
    CutA_3->Draw();
    
    ceMaxSumcand[a]->cd(2);
    ceMaxSumcand[a]->cd(2)->SetLogz();
    hMC_eMaxSumcand_calopfpt_ptselection[a]->SetXTitle("calopt/pfpt");
    hMC_eMaxSumcand_calopfpt_ptselection[a]->SetYTitle("eMax/(chSum+phSum+neSum+muSum)");
    hMC_eMaxSumcand_calopfpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < pfpt < %d", ptBoundary[a], ptBoundary[a+1]),0.2,0.7,14);
    drawText("MC - Jet40",0.25,0.8,14);
    CutA_1->Draw();
    CutA_2->Draw();
    CutA_3->Draw();

    ceMaxSumcand[a]->cd(3);
    ceMaxSumcand[a]->cd(3)->SetLogz();
    hMC_eMaxSumcand_calopfpt_refptselection[a]->SetXTitle("calopt/pfpt");
    hMC_eMaxSumcand_calopfpt_refptselection[a]->SetYTitle("eMax/(chSum+phSum+neSum+muSum)");
    hMC_eMaxSumcand_calopfpt_refptselection[a]->Draw("colz");
    drawText(Form("%d < pfpt, pfrefpt < %d", ptBoundary[a], ptBoundary[a+1]),0.2,0.7,14);
    drawText("MC - Jet40",0.25,0.8,14);
    CutA_1->Draw();
    CutA_2->Draw();
    CutA_3->Draw();

    ceMaxSumcand[a]->SaveAs(Form("../../Output/eMaxOverSumCand_calpfpt_PPJet40_showingCut_fullstat_R0p%d_%d_ptrange_%d.pdf",radius, ptBoundary[a], ptBoundary[a+1]),"RECREATE");
    
  }








}
