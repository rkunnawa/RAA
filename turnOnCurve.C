//#include "HiForestAnalysis/hiForest.h"
#include <iostream>
#include "TProfile.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "TAttMarker.h"
#include "TStyle.h"
#include "TAttLine.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TColor.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
//#include "JEC7tev/get7tevPt.h"

void turnOnCurve(int mode = 0)
{
  TH1::SetDefaultSumw2();

  const int nbins = 6; 
  float bins[nbins+2] = {30,50,80,120,170,220,10000,10000};

  //pp cross sections for 2.76 tev
  float crossSection[nbins+1] = {1.075E-02,1.025E-03,9.865E-05,1.129E-05,1.465E-06,2.837E-07,0};
  //pPb cross sections for 5.02 tev
  
  /*crossSection[0] = 5.335E-01;
    crossSection[1] = 3.378E-02;
    crossSection[2] = 3.778E-03;
    crossSection[3] = 4.412E-04;
    crossSection[4] = 6.147E-05;
    crossSection[5] = 1.018E-05;
    crossSection[6] = 2.477E-06;
    crossSection[7] = 6.160E-07;
    crossSection[8] = 0;
    }*/
  
  /*HiForest *  h[nbins];
    for(int file = 0; file<nbins; file++)
    {
    h[file] = new HiForest(Form("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat%d_Track9_Jet30_matchEqR_merged_forest_0.root",(int)bins[file]),"forest",cPbPb,1);
    //h[file] = new HiForest("/mnt/hadoop/cms/store/user/dgulhan/HIMinBiasUPC-HIRun2011-14Mar2014-v2_tag_HI_MatchEqR_DatabaseJEC_merged2/HIMinBiasUPC-HIRun2011-14Mar2014-v2_tag_HI_MatchEqR_DatabaseJEC.root","forest",cPbPb,1);
    h[file]->PrintStatus(); 
    h[file]->LoadNoTrees();
    h[file]->hasAkVs2CaloJetTree = true;
    h[file]->hasSkimTree = true;
    h[file]->hasHltTree = true;
    h[file]->hasEvtTree = true;
    //h[file]->PrintStatus();
    }*/

  //TFile * inf = new TFile("/mnt/hadoop/cms/store/user/dgulhan/HIMinBiasUPC-HIRun2011-14Mar2014-v2_tag_HI_MatchEqR_DatabaseJEC_merged2/HIMinBiasUPC-HIRun2011-14Mar2014-v2_tag_HI_MatchEqR_DatabaseJEC.root","read");

  int trigger = -1;
  int noise = -1;
  int pcoll = -1;
  float vz = -20;
  float jetpt2[1000]={0};
  float jeteta2[1000]={0};
  float jetpt3[1000]={0};
  float jeteta3[1000]={0};
  float jetpt4[1000]={0};
  float jeteta4[1000]={0};
  float jetpt5[1000]={0};
  float jeteta5[1000]={0};
  int nref2=0;
  int nref3=0;
  int nref4=0;
  int nref5=0;
 
  TH1D * turnon2 = new TH1D("turnon2","",50,0,200);
  TH1D * turnon3 = new TH1D("turnon3","",50,0,200);
  TH1D * turnon4 = new TH1D("turnon4","",50,0,200);
  TH1D * turnon5 = new TH1D("turnon5","",50,0,200);
  TH1D * denom2 = new TH1D("denom2","",50,0,200);
  TH1D * denom3 = new TH1D("denom3","",50,0,200);
  TH1D * denom4 = new TH1D("denom4","",50,0,200);
  TH1D * denom5 = new TH1D("denom5","",50,0,200);

  for(int ifile = 0; ifile<6; ifile++)
    {
      TFile * inf = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat%d_Track9_Jet30_matchEqR_merged_forest_0.root",(int)bins[ifile]),"read");
      TTree * skim = (TTree*)inf->Get("skimanalysis/HltTree");
      TTree * evt  = (TTree*)inf->Get("hiEvtAnalyzer/HiTree"); 
      TTree * jet2 = (TTree*)inf->Get("akPu2PFJetAnalyzer/t");
      TTree * jet3 = (TTree*)inf->Get("akPu3PFJetAnalyzer/t");
      TTree * jet4 = (TTree*)inf->Get("akPu4PFJetAnalyzer/t");
      TTree * jet5 = (TTree*)inf->Get("akPu5PFJetAnalyzer/t");
      TTree * hlt =  (TTree*)inf->Get("hltanalysis/HltTree");

      skim->AddFriend(evt);
      skim->AddFriend(jet2);
      skim->AddFriend(jet3);
      skim->AddFriend(jet4);
      skim->AddFriend(jet5);
      skim->AddFriend(hlt);

      evt->SetBranchAddress("vz",&vz);
      skim->SetBranchAddress("pcollisionEventSelection",&pcoll);
      skim->SetBranchAddress("pHBHENoiseFilter",&noise);
      hlt->SetBranchAddress("HLT_HIJet80_v7",&trigger);
      jet2->SetBranchAddress("jtpt",&jetpt2);
      jet2->SetBranchAddress("jteta",&jeteta2);
      jet3->SetBranchAddress("jtpt",&jetpt3);
      jet3->SetBranchAddress("jteta",&jeteta3);
      jet4->SetBranchAddress("jtpt",&jetpt4);
      jet4->SetBranchAddress("jteta",&jeteta4);
      jet5->SetBranchAddress("jtpt",&jetpt5);
      jet5->SetBranchAddress("jteta",&jeteta5); 
      jet2->SetBranchAddress("nref",&nref2);
      jet3->SetBranchAddress("nref",&nref3);
      jet4->SetBranchAddress("nref",&nref4);
      jet5->SetBranchAddress("nref",&nref5);


      int nEntries = skim->GetEntries();
      //int nEntries = 100;
      for(int i = 0; i<nEntries; i++)
	{
	  if(i%1000==0) std::cout << i << "/" << skim->GetEntries() << std::endl;
	  skim->GetEntry(i);
	  if(pcoll==0) continue;
	  if(vz>15 || vz<-15) continue;
	  //if(noise==0) continue;
	  //std::cout << noise << " " << pcoll << " " << trigger << " " << std::endl;
	  //std::cout << jetpt3[0] << " " << jetpt4[0] << " " << jetpt5[0] << std::endl;
	  //std::cout << nref3 << nref4 << nref5 << std::endl;
	  if(nref2!=0)
	    {
	      if(TMath::Abs(jeteta2[0])<2)
		{
		  denom2->Fill(jetpt2[0]);
		  if(trigger) turnon2->Fill(jetpt2[0]);
		}
	    }
	  if(nref3!=0)//TMath::Abs(jeteta3[0])<1.6)
	    { 
	      if(TMath::Abs(jeteta3[0])<2)
		{
		  denom3->Fill(jetpt3[0]);
		  if(trigger) turnon3->Fill(jetpt3[0]);
		}
	    }
	  if(nref4!=0)//TMath::Abs(jeteta4[0])<1.6)
	    { 
	      if(TMath::Abs(jeteta4[0])<2)
		{
		  denom4->Fill(jetpt4[0]);
		  if(trigger) turnon4->Fill(jetpt4[0]);
		}
	    }
	  if(nref5!=0)//TMath::Abs(jeteta5[0])<1.6)
	    {
	      if(TMath::Abs(jeteta5[0])<2)
		{
		  denom5->Fill(jetpt5[0]);
		  if(trigger) turnon5->Fill(jetpt5[0]);
		}
	    }
	}
      inf->Close();
    }
  /*int maxFile = nbins;
    for(int ifile = 0; ifile<maxFile; ifile++)
    {
    std::cout << h[ifile] << std::endl;
    int nEntries = h[ifile]->GetEntries();
    //int  nEntries = 5000;
    std::cout << nEntries << std::endl;
    for(int i = 0; i<50; i++)
    {
    h[ifile]->GetEntry(i);
    std::cout << h[ifile]->skim.pcollisionEventSelection << std::endl;
    if(i%1000 == 0){ std::cout << i << "/" << nEntries << std::endl;
    std::cout << h[ifile]->akVs2Calo.jteta[0] << std::endl;}     

    if(!(h[ifile]->skim.pcollisionEventSelection == 1)||(TMath::Abs(h[ifile]->evt.vz)>15)) continue; 
  
    if(TMath::Abs(h[ifile]->akVs2Calo.jteta[0])<1.6)
    {
    denom2->Fill(h[ifile]->akVs2Calo.jtpt[0]);
    if(h[ifile]->hlt.HLT_HIJet80_v7) turnon2->Fill(h[ifile]->akVs2Calo.jtpt[0]);
    }
    h[ifile]->hasAkVs2CaloJetTree = false;
    h[ifile]->hasAkVs3CaloJetTree = true;
    h[ifile]->GetEntry(i);

    if(i%10000 == 0){std::cout << h[ifile]->akVs2Calo.jteta[0] << std::endl;}
    if(TMath::Abs(h[ifile]->akVs3Calo.jteta[0])<1.6)
    {
    denom3->Fill(h[ifile]->akVs3Calo.jtpt[0]);
    if(h[ifile]->hlt.HLT_HIJet80_v7) turnon3->Fill(h[ifile]->akVs3Calo.jtpt[0]);
    }
    h[ifile]->hasAkVs3CaloJetTree = false;
    h[ifile]->hasAkVs4CaloJetTree = true;
    h[ifile]->GetEntry(i);

    if(TMath::Abs(h[ifile]->akVs4Calo.jteta[0])<1.6)
    {
    denom4->Fill(h[ifile]->akVs4Calo.jtpt[0]);
    if(h[ifile]->hlt.HLT_HIJet80_v7) turnon4->Fill(h[ifile]->akVs4Calo.jtpt[0]);
    }
    h[ifile]->hasAkVs4CaloJetTree = false;
    h[ifile]->hasAkVs5CaloJetTree = true;
    h[ifile]->GetEntry(i);

    if(TMath::Abs(h[ifile]->akVs5Calo.jteta[0])<1.6)
    {
    denom5->Fill(h[ifile]->akVs5Calo.jtpt[0]);
    if(h[ifile]->hlt.HLT_HIJet80_v7) turnon5->Fill(h[ifile]->akVs5Calo.jtpt[0]);
    }
    h[ifile]->hasAkVs5CaloJetTree = false;
    h[ifile]->hasAkVs2CaloJetTree = true;
    }
    }
  */
  TGraphAsymmErrors * turnon2Asym;// = new TGraphAsymmErrors(turnon40);
  TGraphAsymmErrors * turnon3Asym;// = new TGraphAsymmErrors(turnon80);
  TGraphAsymmErrors * turnon4Asym;// = new TGraphAsymmErrors(turnon30);
  TGraphAsymmErrors * turnon5Asym;// = new TGraphAsymmErrors(turnon60);
 

  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  TMultiGraph * mg = new TMultiGraph();
  mg->SetTitle(";lead jet p_{T}^{reco};Efficiency");

  turnon2Asym = new TGraphAsymmErrors();
  turnon2Asym->SetName("turnon2Asym");
  turnon2Asym->BayesDivide(turnon2,denom2);
  for(int i =1; i<41; i++)
    {
      turnon2Asym->SetPointEXlow(i,0);
      turnon2Asym->SetPointEXhigh(i,0);
    }
  turnon2Asym->SetLineColor(1);
  turnon2Asym->SetMarkerSize(0.8);

  turnon3Asym = new TGraphAsymmErrors();
  turnon3Asym->SetName("turnon3Asym");
  turnon3Asym->BayesDivide(turnon3,denom3);
  for(int i =1; i<41; i++)
    {
      turnon3Asym->SetPointEXlow(i,0);
      turnon3Asym->SetPointEXhigh(i,0);
    }
  turnon3Asym->SetLineColor(kRed+1);
  turnon3Asym->SetMarkerColor(kRed+1);
  turnon3Asym->SetMarkerSize(0.8);

  turnon4Asym = new TGraphAsymmErrors();
  turnon4Asym->SetName("turnon4Asym");
  turnon4Asym->BayesDivide(turnon4,denom4);
  for(int i =1; i<41; i++)
    {
      turnon4Asym->SetPointEXlow(i,0);
      turnon4Asym->SetPointEXhigh(i,0);
    }
  turnon4Asym->SetLineColor(kBlue+1);
  turnon4Asym->SetMarkerColor(kBlue+1);
  turnon4Asym->SetMarkerSize(0.8);

  turnon5Asym = new TGraphAsymmErrors();
  turnon5Asym->SetName("turnon5Asym");
  turnon5Asym->BayesDivide(turnon5,denom5);
  for(int i =1; i<41; i++)
    {
      turnon5Asym->SetPointEXlow(i,0);
      turnon5Asym->SetPointEXhigh(i,0);
    }
  turnon5Asym->SetLineColor(kGreen+1);
  turnon5Asym->SetMarkerColor(kGreen+1);
  turnon5Asym->SetMarkerSize(0.8);

  mg->Add(turnon2Asym,"");
  mg->Add(turnon3Asym,"");
  mg->Add(turnon4Asym,"");
  mg->Add(turnon5Asym,"");
  mg->Draw("AP");

  TLegend * leg = new TLegend(0.80,0.2,0.9,0.4);
  leg->AddEntry((TObject*)0,"2.76 TeV PYTHIA+HYJET","");
  leg->AddEntry(turnon2Asym,"R=0.2 Jet80 Trigger","p");
  leg->AddEntry(turnon3Asym,"R=0.3 Jet80 Trigger","p");
  leg->AddEntry(turnon4Asym,"R=0.4 Jet80 Trigger","p");
  leg->AddEntry(turnon5Asym,"R=0.5 Jet80 Trigger","p");
  leg->Draw("same");
  c2->SaveAs("TriggerTurnOn_PbPbMC_etaLT2_Jet80.png");
  c2->SaveAs("TriggerTurnOn_PbPbMC_etaLT2_jet80.pdf");
  
}
