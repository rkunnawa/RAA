// Raghav Kunnawalam Elayavalli
// June 23rd 2014
// CERN

// 
// Macro separated from RAA_read_data to just read in the pp datasets - due to condor submission details. look at the read_data_pbpb macro to get more information. 
// 

#include <iostream>
#include <stdio.h>
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

static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}

using namespace std;


void RAA_read_data_pp(int radius = 3){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TStopwatch timer;
  timer.Start();

  cout<<"Reading PP 2013 data"<<endl;
  cout<<"Running for Radius = "<<radius<<endl;

   // data files - pp 
  TFile *fpp1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet80_v2.root");
  TFile *fpp2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet40_v2.root");


    //do it for the pp - need to check up on this. 
  TTree *jetpp1_v2 = (TTree*)fpp1_v2->Get(Form("jetR%d",radius));
  TTree *jetpp2_v2 = (TTree*)fpp2_v2->Get(Form("jetR%d",radius));

  TTree *evtpp1_v2 = (TTree*)fpp1_v2->Get("evt");
  TTree *evtpp2_v2 = (TTree*)fpp2_v2->Get("evt");

  jetpp1_v2->AddFriend(evtpp1_v2);
  jetpp2_v2->AddFriend(evtpp2_v2);


  //get all the pp spectra here: 
  TCut pp3 = "abs(eta)<2&&jet40&&!jet60&&!jet80&&(chMax/pt)>0.01&&(TMath::Max(neMax,chMax)/TMath::Max(chSum,neSum))>=0.975";
  
  TH1F *hpp1 = new TH1F("hpp1","Spectra from Jet 80",1000,0,1000);
  TH1F *hpp2 = new TH1F("hpp2","Spectra from Jet 60 & !Jet80",1000,0,1000);
  TH1F *hpp3 = new TH1F("hpp3","Spectra from Jet 40 & !Jet60 & !Jet80",1000,0,1000);
  TH1F *hppComb = new TH1F("hppComb","Combined spectra according to 12003 method",1000,0,1000);
  
  //get the prescl factor information. 
  //Float_t presclpbpb3 = (Float_t)jetpbpb1_v2->GetEntries("jet80")/jetpbpb1_v2->GetEntries("jet55&&jet80");
  //cout<<"pbpb prescl3 = "<<presclpbpb3<<endl;//1.99871
  //Float_t presclpp3 = (Float_t)jetpp1_v2->GetEntries("jet80")/jetpp1_v2->GetEntries("jet40&&jet80");
  //cout<<"pp prescl3 = "<<presclpp3<<endl; //9.24968

  //root [9] (Float_t)jet->GetEntries("HLT_HIJet80_v1")/jet->GetEntries("HLT_HIJet80_v1&&HLT_HIJet55_v1")
  //(double)2.34995051108819153e+00
  //ive commented this below - to just check for the pbpb histograms to load. 
  
  // include the neutralMax/ max(chargedSum, neutralSum)>0.975 cut here

  // because whatever cut we use in PbPb, we have to use in pp. 

  jetpp1_v2->Project("hpp1","pt","abs(eta)<2&&jet80&&(chMax/pt)>0.01&&(TMath::Max(chMax,neMax)/TMath::Max(chSum,neSum))>=0.975");
  hpp1->Print("base");
 
  jetpp2_v2->Project("hpp2","pt","abs(eta)<2&&jet60&&!jet80&&(chMax/pt)>0.01&&(TMath::Max(chMax,neMax)/TMath::Max(chSum,neSum))>=0.975");
  hpp2->Print("base");

  jetpp2_v2->Project("hpp3","pt","9.25038"*pp3);
  // 9.25038 was the value. 
  //jetpp2_v2->Project("hpp3","pt","jet40_p"*pp3);
  hpp3->Print("base");
 
  
  hpp1->Scale(1./5300e6);//pp lumi
  hpp2->Scale(1./5300e6);
  hpp3->Scale(1./5300e6);

  hpp1->Scale(1./4);//delta eta
  hpp2->Scale(1./4);
  hpp3->Scale(1./4);

  hppComb->Add(hpp1,1);
  hppComb->Add(hpp2,1);
  hppComb->Add(hpp3,1);
  hppComb->Print("base");

  hppComb = (TH1F*)hppComb->Rebin(nbins_pt,"hppComb",boundaries_pt);
  hpp3 = (TH1F*)hpp3->Rebin(nbins_pt,"hpp3",boundaries_pt);
  hpp2 = (TH1F*)hpp2->Rebin(nbins_pt,"hpp2",boundaries_pt);
  hpp1 = (TH1F*)hpp1->Rebin(nbins_pt,"hpp1",boundaries_pt);

  divideBinWidth(hppComb);
  divideBinWidth(hpp1);
  divideBinWidth(hpp2);
  divideBinWidth(hpp3);

  // output file declaration

  TDatime date;

  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Output/pp_data_ak%d_%d_chMax_12003cut.root",radius,date.GetDate()),"RECREATE");
  f.cd();

  hpp1->Write();
  hpp2->Write();
  hpp3->Write();
  hppComb->Write();
  
  f.Write();

  f.Close();


  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<timer.CpuTime()<<endl;
  cout<<"Real time (min) = "<<timer.RealTime()<<endl;
  

}
