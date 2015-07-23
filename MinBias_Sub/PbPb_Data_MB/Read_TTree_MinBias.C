//Based on Raghav's  RAA_read_jetHistograms.C code

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
#include "TLorentzVector.h"
#include "TVector3.h"
#include "../Headers/plot.h"


static const int nbins_cent = 7;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
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

static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 355, 400};


float deltaphi(float phi1, float phi2); 
float deltaR(float eta1, float phi1, float eta2, float phi2);
Double_t GetdPhi(Double_t mphi,Double_t vphi);

using namespace std;


void Read_TTree_MinBias(char* etaWidth = (char*)"20_eta_20",
			Int_t radius = 3,
			Int_t etaLow = 20,
			Int_t etaHigh = 20)
{
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // char * etaWidth = (char*)Form("%d_eta_%d",etaLow, etaHigh);
  cout<<"etaWidth = "<<etaWidth<<endl;
  cout<<"Radius = "<<radius<<endl;

  bool isSymm = false;
  if(etaLow == etaHigh) isSymm = true;
   
  // the cut is a 3 step cut based on the different value of the calopt/pfpt - copy the following lines into your loop (with the corresponding branch address set)
  // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) hGood->Fill();
  // if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) ) hGood->Fill();
  // if(calopt/pfpt > 0.85 & eMax/Sumcand > 0.9) hGood->Fill();
  
  TFile * fData, * fMC;
  TTree * Data_matched, * Data_unmatched, *MC_matched, * MC_unmatched;
  fData = TFile::Open("/mnt/hadoop/cms/store/user/pawan/trees/JetRaaTree_ak234_pp_DATA.root");
  //  fData = TFile::Open("/mnt/hadoop/cms/store/user/pawan/trees/JetRaaTree_akPu234_PbPb_DATA.root");
  //   fData = TFile::Open("/mnt/hadoop/cms/store/user/pawan/trees/JetRaaTree_akPu234_PbPb_MinBias_DATA.root");
  //fMC = TFile::Open("/mnt/hadoop/cms/store/user/pawan/trees/JetRaaTree_akPu234_PbPb_MC.root");
  // Data_matched = (TTree*)fData->Get(Form("akPu%dJetAnalyzer/jetTree",radius));
  Data_matched = (TTree*)fData->Get(Form("ak%dJetAnalyzer/jetTree",radius));
  //  MC_matched = (TTree*)fMC->Get(Form("akPu%dJetAnalyzer/jetTree",radius));
  TFile *fcentin = TFile::Open("MinBiasHLTRatioWeightsTTree.root");
  TH1F *hCentWeight = (TH1F*)fcentin->Get("hRatio");
  TFile* fDataFlag = TFile::Open(Form("Pawan_TTree_PbPb_Data_MC_DupEvents_flag_R0p%d.root", radius));
  // TFile* fDataFlag = TFile::Open(Form("/export/d00/scratch/rkunnawa/rootfiles/RAA/Pawan_TTree_PbPb_MBData_DupEvents_flag_R0p%d.root", radius));
  TTree* flagt=(TTree*)fDataFlag->Get("jetDup_Data");

  Data_matched->AddFriend(flagt);
  
  // TFile fout1(Form("Pawan_TTree_PbPb_Data_MC_DupEvents_flag_R0p%d.root", radius),"RECREATE");
  // fout1.cd();
  
  //get the spectra with the specific trigger object from the different files. 
  TH1F *hpbpb_TrgObj80[nbins_cent];
  TH1F *hpbpb_TrgObj65[nbins_cent];
  TH1F *hpbpb_TrgObj55[nbins_cent];

  TH1F *hpbpb_TrgObjComb[nbins_cent];
  TH1F *hpbpb_TrgObjMB[nbins_cent];
  TH1F *hpbpb_TrgObjMBArray[nbins_cent];
  TH1F *hpbpb_TrgObjMBwoLJ[nbins_cent];
  TH1F *hpbpb_TrgObjMBLJ[nbins_cent];
  TH1F *hpbpb_TrgObjMBwoHT[nbins_cent];
  
  TH2F *hdphiptcent[nbins_cent];
  //TH2F *hdphiptMBwoHLTcent[nbins_cent];
  //  TH1F *hBin_[nbins_cent];
  TH1F *hEvent_Vz_[nbins_cent];
  
  TH1F * hEvent_Vz = new TH1F("hEvent_Vz","Primary Vertex Z",400,-20,20);
  TH1F * hBin = new TH1F("hBin","Centrality Bins",200,0,200);

  TH1F *heta = new TH1F("heta","eta distribution",100,-2.5,2.5);
  TH1F *hphi = new TH1F("hphi","phi distribution",100,-3.5,3.5);
  TH1F *hdphi = new TH1F("hdphi","delta phi distribution",70,0,3.5);
  TH2F *hdphipt=new TH2F("hdphipt","pt vs delta phi distribution",100,0,200,70,0,3.5);
  TH1F *hdptratio = new TH1F("hdptratio","pt ratio distribution",100,0,5);

  
  for(int i = 0;i<nbins_cent;++i){
    //cout<<"cent bin = "<<i<<endl;
    //  hBin_[i] = new TH1F(Form("hBin_cent_%d",i),Form("HBin Values %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),200,0,200);
    hEvent_Vz_[i] = new TH1F(Form("hVz_cent_%d",i),Form("HVz Values %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),400,-20,20);
      
    hpbpb_TrgObj80[i] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObj65[i] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObj55[i] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjComb[i] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjMB[i] = new TH1F(Form("hpbpb_HLTMB_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjMBArray[i] = new TH1F(Form("hpbpb_HLTMBArray_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB ArrayR%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjMBwoLJ[i] = new TH1F(Form("hpbpb_HLTMBwoLJ_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB without LJ R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjMBLJ[i] = new TH1F(Form("hpbpb_HLTMBLJ_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB Leading Jet R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    
    hpbpb_TrgObjMBwoHT[i] = new TH1F(Form("hpbpb_HLTMBwoHT_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB without HT R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);


    
    hdphiptcent[i]=new TH2F(Form("hdphipt_cent%d",i),Form("hdphipt%2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200,160,0,3.2);
    // hdphiptwoHLTcent[i]=new TH2F(Form("hdphiptwoHLT_cent%d",i),Form("hdphipt%2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200,70,0,3.5);
  }// centrality bin loop
  

  // Define all the histograms necessary for the analysis: 
  
  // 1 - Data,
  Float_t pfpt_1[1000];
  Float_t calopt_1[1000]; 
  Int_t npf_1;
  Int_t isCaloMatch_1[1000];
  Float_t eMax_1[1000];
  Float_t trMax_1[1000];
  Float_t chMax_1[1000];
  Float_t chSum_1[1000];
  Float_t phSum_1[1000];
  Float_t neSum_1[1000];
  Float_t muSum_1[1000];
  Int_t jet55_1, jet65_1, jet80_1, jetMB_1;
  Int_t jet55_p_1, jet65_p_1, jet80_p_1;
  Double_t weight;
  Int_t hiBin_1;
  Float_t vz_1;
  Float_t eta_1[1000];
  Float_t phi_1[1000];
  Float_t pfrawpt_1[1000];
  Int_t run, evt, lumi;
  Int_t dup[1000];

  flagt->SetBranchAddress("jetDup",&dup);
  
  Data_matched->SetBranchAddress("run_value", &run);
  Data_matched->SetBranchAddress("evt_value", &evt);
  Data_matched->SetBranchAddress("lumi_value", &lumi);
 
  Data_matched->SetBranchAddress("vz", &vz_1);
  Data_matched->SetBranchAddress("npf", &npf_1);
  Data_matched->SetBranchAddress("isCaloMatch", &isCaloMatch_1);
  Data_matched->SetBranchAddress("calopt",&calopt_1);
  Data_matched->SetBranchAddress("pfpt",&pfpt_1);
  Data_matched->SetBranchAddress("eMax",&eMax_1);
  Data_matched->SetBranchAddress("chMax",&chMax_1);
  Data_matched->SetBranchAddress("trMax",&trMax_1);
  Data_matched->SetBranchAddress("chSum",&chSum_1);
  Data_matched->SetBranchAddress("phSum",&phSum_1);
  Data_matched->SetBranchAddress("neSum",&neSum_1);
  Data_matched->SetBranchAddress("muSum",&muSum_1);
  Data_matched->SetBranchAddress("hiBin",&hiBin_1);
  Data_matched->SetBranchAddress("jet55",&jet55_1);
  Data_matched->SetBranchAddress("jet65",&jet65_1);
  Data_matched->SetBranchAddress("jet80",&jet80_1);
  Data_matched->SetBranchAddress("jetMB",&jetMB_1);

  // Data_matched->SetBranchAddress("jet40",&jet40_1);
  //Data_matched->SetBranchAddress("jet60",&jet60_1);
  
  // Data_matched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  //Data_matched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  //Data_matched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  Data_matched->SetBranchAddress("pfeta",&eta_1);
  Data_matched->SetBranchAddress("pfphi",&phi_1);
  //  Data_matched->SetBranchAddress("pfrawpt",&pfrawpt_1);
  
  
  //TTree * jetDup_Data = new TTree("jetDup_Data","");
  // Int_t jetDup[1000];
  //Int_t nPFjets;
  // jetDup_Data->Branch("nPFjets",&nPFjets,"nPFjets/I");
  // jetDup_Data->Branch("jetDup",&jetDup, "jetDup[nPFjets]/I");
  
  // data loop
  // long entries = Data_matched->GetEntries();
  long entries = 200000;
  Float_t Jet55_prescl = 2.0475;

  // get the random value for smear systematics, TRandom rnd, value per jet = rnd.Gaus(0,1);
  TRandom rnd; 
  TH1F * htest = new TH1F("htest","",nbins_pt, boundaries_pt);

  // cout<<"matched Data ntuple "<<endl;

  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;

    //  cout<<nentry<<"/"<<entries<<endl;
    Data_matched->GetEntry(nentry);
    // Int_t cBin =1;

  
    // nPFjets = npf_1;
    
    Int_t cBin = findBin(hiBin_1);
    if(cBin == -1 || cBin >= nbins_cent) {
      cout<<" problematic centrality bin!!!!  "<< hiBin_1<<endl;
      continue;
    }
    if(jetMB_1 == 1)
      //  if(jet55_1 == 1 || jet65_1 == 1 || jet80_1 == 1)
      {
	hBin->Fill(hiBin_1);
	hEvent_Vz->Fill(vz_1);
	hEvent_Vz_[cBin]->Fill(vz_1);
      }
    

    Float_t weight_cent = hCentWeight->GetBinContent(hCentWeight->FindBin(hiBin_1));
    Float_t wght = weight_cent;
       
       
    //Float_t wght=1;
    // if(etaBoundary == 1.0){ 
    //   if(radius == 2) wght = R2nCorr[cBin][htest->FindBin(pfpt_1)];
    //   if(radius == 3) wght = R3nCorr[cBin][htest->FindBin(pfpt_1)];
    //   if(radius == 4) wght = R4nCorr[cBin][htest->FindBin(pfpt_1)];  
    // }

    // continue;
    TLorentzVector jetarray[1000];
    for(int i=0;i<1000;i++){
      jetarray[i].SetPtEtaPhiM(0.001,0,0,0);
    }


    //  TVector3 jetarray[1000];
    //for(int i=0;i<1000;i++){
    //  jetarray[i].SetPtEtaPhi(0.001,0,0);
    // }
    



    /* vector <TLorentzVector> jetVecA;
       TLorentzVector TV;
       TV.SetPtEtaPhiM(0.,0.,0.,0.);
       jetVecA.push_back(TV);
    */


    
    //   vector  <int> indj;
    //vector  <float> indpt;
    //indj.push_back(0);
    //indpt.push_back(0.);

    // int pfdup_loc = 0; 
      
    if(npf_1>1000)cout<<"npf_1 bigger than 1000"<<endl;
    
    for(int g = 0; g<npf_1; ++g){ // jet loop
      //if(dup[g]>0)    cout<<"dup "<<dup[g]<<"    "<<g<<endl;     
      if(dup[g]==1) continue;
      
      if(isSymm && TMath::Abs(eta_1[g]) > (Float_t)etaHigh/10) continue;       
      if(!isSymm && (TMath::Abs(eta_1[g]) < (Float_t)etaLow/10 || TMath::Abs(eta_1[g]) > (Float_t)etaHigh/10)) continue;

      // jetDup[g] = 0;
      
      if(pfpt_1[g] <15) continue;
      
      Float_t Sumcand = chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g];
      //  cout<<"testing the pt values   "<<pfpt_1[g]<<endl;
      if(isCaloMatch_1[g] == 1){
	jetarray[g].SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {
	  
	  if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1[g]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) {
	    hpbpb_TrgObj55[cBin]->Fill(pfpt_1[g], Jet55_prescl* wght);
	  }
	  if(calopt_1[g]/pfpt_1[g] > 0.85){
	    hpbpb_TrgObj55[cBin]->Fill(pfpt_1[g], Jet55_prescl* wght);
	  }
	  if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1[g]/Sumcand < 0.05) {
		  
	    hpbpb_TrgObj55[cBin]->Fill(pfpt_1[g], Jet55_prescl* wght);
	  }
	}
	
	if(jet65_1 == 1 && jet80_1 == 0 ) {
	  
	  if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1[g]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)){
	    
	    hpbpb_TrgObj65[cBin]->Fill(pfpt_1[g], wght);
	  }
	  if(calopt_1[g]/pfpt_1[g] > 0.85) {
	    hpbpb_TrgObj65[cBin]->Fill(pfpt_1[g], wght);
	  }
	  if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1[g]/Sumcand < 0.05) {
	    hpbpb_TrgObj65[cBin]->Fill(pfpt_1[g], wght);
	  }
	}
	
	if(jet80_1 == 1) {
	  if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1[g]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) {
	    hpbpb_TrgObj80[cBin]->Fill(pfpt_1[g], wght);
	  }
	  if(calopt_1[g]/pfpt_1[g] > 0.85){
	    hpbpb_TrgObj80[cBin]->Fill(pfpt_1[g], wght);
	  }
	  if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1[g]/Sumcand < 0.05){
	    hpbpb_TrgObj80[cBin]->Fill(pfpt_1[g], wght);
	  }
	}
	
	
	if(jetMB_1==1){// || jet65_1 == 1 || jet55_1 == 1   ) {
	  
	  //  if (nentry==5703 ||nentry==18159 || nentry==21380) cout<<"no cuts  "<<pfpt_1[g]<<"   "<<eta_1[g]<<"   "<<phi_1[g]<<"   "<<g<<" "<<run<<"   "<<evt<<"   "<<lumi<<endl; 
	  
	  
	  
	  if(calopt_1[g]/pfpt_1[g] > 0.5 && calopt_1[g]/pfpt_1[g] <= 0.85 && eMax_1[g]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1[g]/pfpt_1[g] - (Float_t)9/7)) {
	    //if(g!=0) indj.push_back(g);
	    //indpt.push_back(pfpt_1[g]);
	    // jetarray[g].SetPerp(pfpt_1[g]);
	    //   jetarray[g].SetPtEtaPhi(pfpt_1[g],eta_1[g],phi_1[g]);
	    jetarray[g].SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	    
	    // if (nentry==5703 ||nentry==18159 || nentry==21380) cout<<pfpt_1[g]<<"   "<<eta_1[g]<<"   "<<g<<""<<endl; 
	    
	    //  jetVecA.at(g).SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	    if(trMax_1[g]>8) hpbpb_TrgObjMBwoHT[cBin]->Fill(pfpt_1[g], wght);
	    hpbpb_TrgObjMB[cBin]->Fill(pfpt_1[g], wght);
	  }
	  if(calopt_1[g]/pfpt_1[g] > 0.85){
	    // if (nentry==5703 ||nentry==18159 || nentry==21380) cout<<"matched loop 2th cut"<<endl;
	    //if(g!=0) indj.push_back(g);
	    //indpt.push_back(pfpt_1[g]);
	    // jetarray[g].SetPtEtaPhi(pfpt_1[g],eta_1[g],phi_1[g]);
	    jetarray[g].SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	    // jetVecA.at(g).SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	    if(trMax_1[g]>8) hpbpb_TrgObjMBwoHT[cBin]->Fill(pfpt_1[g], wght);
	    hpbpb_TrgObjMB[cBin]->Fill(pfpt_1[g], wght);
	    //  if (nentry==5703 ||nentry==18159 || nentry==21380) cout<<pfpt_1[g]<<"   "<<eta_1[g]<<"   "<<g<<endl; 
	    
	  }
	  if(calopt_1[g]/pfpt_1[g] <= 0.5 && eMax_1[g]/Sumcand < 0.05){
	    
	    //  if (nentry==5703 ||nentry==18159 || nentry==21380) cout<<"matched loop 3rd cut"<<endl;
	    //  if(g!=0)  indj.push_back(g);
	    //indpt.push_back(pfpt_1[g]);
	    // jetarray[g].SetPtEtaPhi(pfpt_1[g],eta_1[g],phi_1[g]);
	    jetarray[g].SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0); 
	    // jetVecA.at(g).SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	    if(trMax_1[g]>8) hpbpb_TrgObjMBwoHT[cBin]->Fill(pfpt_1[g], wght);
	    hpbpb_TrgObjMB[cBin]->Fill(pfpt_1[g], wght);
	    
	    
	    //  if (nentry==5703 ||nentry==18159 || nentry==21380) cout<<pfpt_1[g]<<"   "<<eta_1[g]<<"   "<<g<<endl;  
	  }
	}
	
      }
      
      // fill the histograms with unmatched Jets 
      if(isCaloMatch_1[g] ==0) {
	//	jetarray[g].SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {    

	  if(eMax_1[g]/Sumcand < 0.05 ){
	    hpbpb_TrgObj55[cBin]->Fill(pfpt_1[g], Jet55_prescl*wght);
	  }
	}

	if(jet65_1 == 1 && jet80_1 == 0 ) {
	  if(eMax_1[g]/Sumcand < 0.05  ){
	    
	    hpbpb_TrgObj65[cBin]->Fill(pfpt_1[g], wght);
	  }
	}

	if(jet80_1 == 1) {
	  if(eMax_1[g]/Sumcand < 0.05  ){
	    hpbpb_TrgObj80[cBin]->Fill(pfpt_1[g], wght);
	  }
	}
	
	if(jetMB_1==1){// || jet65_1 == 1 || jet55_1 == 1 ) {
	  if(eMax_1[g]/Sumcand < 0.05  ){
	    if(trMax_1[g]>8) hpbpb_TrgObjMBwoHT[cBin]->Fill(pfpt_1[g], wght);
	    hpbpb_TrgObjMB[cBin]->Fill(pfpt_1[g], wght);
	    
	    // if (nentry==5703 ||nentry==18159 || nentry==21380) cout<<pfpt_1[g]<<"   "<<eta_1[g]<<"   "<<g<<endl; 
	    //if(g!=0)  indj.push_back(g);
	    //indpt.push_back(pfpt_1[g]);
	    // jetarray[g].SetPtEtaPhi(pfpt_1[g],eta_1[g],phi_1[g]);
	    jetarray[g].SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	    // jetVecA.at(g).SetPtEtaPhiM(pfpt_1[g],eta_1[g],phi_1[g],0);
	  }
	}
		
		
      }
      
      
      // if(jetarray[g].Phi()==2.98746 ||jetarray[g].Phi()==2.61922 || jetarray[g].Phi()==0.994431)
      //  if(nentry==5703 || nentry==18159 || nentry==21380)	cout<<jetarray[g].Perp()<<"   "<<jetarray[g].Eta()<<"   "<<jetarray[g].Phi()<<"    "<<g<<endl;

      //if(nentry==19000 ) cout<<"next event"<<endl;
      // if(nentry==20000 ) cout<<"next event"<<endl;
    }// jet loop
    //if(nentry==19000 ) cout<<"next event"<<endl;
    //  if(nentry==20000 ) cout<<"next event"<<endl;
    int flagt=1;
    float temppt=0;
    float tempeta=0;
    float tempphi=0;
    
    for(int g = 1; (g<1000) &&flagt; ++g){ // TLorentzVector Loop
      flagt=0;

      // Float_t ptvalue=jetarray[g].Perp(); 
      //      if(ptvalue!=0)cout<<"  ptvalue from lorentz vector  "<<ptvalue<<endl;
      
      for(int h = 0; h<1000; ++h){ // TLorentzVector Loop

	Float_t ptvalue=jetarray[h+1].Perp(); 
	if(ptvalue>jetarray[h].Perp()) {
	  temppt=jetarray[h].Perp();
	  tempeta=jetarray[h].Eta();
	  tempphi=jetarray[h].Phi();
	  jetarray[h]=jetarray[h+1];
	  
	  //  jetarray[h+1].SetPtEtaPhi(temppt,tempeta,tempphi);
	  jetarray[h+1].SetPtEtaPhiM(temppt,tempeta,tempphi,0);
	  flagt=1;
	  // cout<<"found one   "<<ptvalue<<"  h "<<h<<"   "<<g<<endl; 
	}
	//	else continue;
		
      }
    }

    
    /*   for(int g = 0; g<npf_1; ++g){

	 if(jetarray[0].Perp() == pfpt_1[g]){
	 pfdup_loc = g;
	 break;
	 }
      
	 }*/

  
    for(int g = 1; (g<=npf_1); ++g){

      //if(jetarray[g].Perp()<15) continue;
      
      if(jetarray[0].Perp()<jetarray[g].Perp()) cout<<"problem with ordering pt bigger than the leading jet"<<endl;

      // cout<<jetarray[g-1].Phi()<<"   "<<g-1<<endl;
      //	Float_t deltaval=0;
      Float_t deltaval=deltaphi(jetarray[0].Phi(),jetarray[g].Phi());
      // Float_t deltaval=GetdPhi(jetarray[0].Phi(),jetarray[g].Phi());
      hpbpb_TrgObjMBwoLJ[cBin]->Fill(jetarray[g].Perp(),wght);
      
      //if(jetarray[0].Perp()>35 && jetarray[g].Perp()>30){
      if(deltaval<0.01 && jetarray[0].Perp()==jetarray[g].Perp() && jetarray[g].Phi()!=0 ) {
	//	jetDup[g] = 1;

       	cout<<jetarray[g].Perp() <<"   "<<jetarray[0].Phi()<<"   "<<jetarray[g].Phi()<<"   "<<g<<"   "<<deltaval<<"    "<<jetarray[0].Eta()<<"    "<<jetarray[g].Eta()<<"   "<<"   "<<g<<" "<<run<<"   "<<evt<<"   "<<lumi<<endl;
      }
      hdphi->Fill(deltaval,wght);
      hdphipt->Fill(jetarray[g].Perp(),deltaval,wght);
      hdphiptcent[cBin]->Fill(jetarray[g].Perp(),deltaval,wght);
      
      //}
      
    }

    // jetDup_Data->Fill();
    
    // if(jetarray[0].Perp()>15) hpbpb_TrgObjMBLJ[cBin]->Fill(jetarray[0].Perp(),wght);
     
    
    //  std::sort(jetVecA.begin(), jetVecA.end());


    /*  cout << "DEBUG before: a = ("; 
	for (size_t i = 0; i != indj.size(); ++i) {
	cout << indj.at(i);
	if (i != (indj.size()-1))
	cout << ", ";
	} // for (i)
	cout << ")" << endl;
    */

    /*  cout << "DEBUG before: a = ("; 
	for (size_t i = 0; i != indpt.size(); ++i) {
	cout << indpt.at(i);
	if (i != (indpt.size()-1))
	cout << ", ";
	} // for (i)
	cout << ")" << endl;
    */

    
    
    //  std::sort(indpt.rbegin(), indpt.rend());

    //if(indpt.at(0)!=0)hpbpb_TrgObjMBLJ[cBin]->Fill(indpt.at(0));
  
    /* for (size_t i = 0; i != indj.size(); ++i) {
       if(i!=0 && indj.at(i)!=0) hpbpb_TrgObjMBwoLJ[cBin]->Fill(indj.at(i));
       if(indj.at(i)!=0)hpbpb_TrgObjMBArray[cBin]->Fill(indj.at(0));
       }*/

    // cout<<"event "<<endl;
    
  }//// data ntuple loop

  //fout1.cd();
  //jetDup_Data->Write();
  //fout1.Close();


  /*
    TFile *fout=new TFile("test.root","RECREATE");
    //  TFile *fout=new TFile(Form("TTree_PbPb_MBwoHLT_spectra_finebins_wohtp_%s_R0p%d.root",etaWidth,radius),"RECREATE");
    // TFile *fout=new TFile(Form("TTree_PbPb_HLT_spectra_finebins_%s_R0p%d.root",etaWidth,radius),"RECREATE");
    fout->cd();
 
    hEvent_Vz->Write();
    hBin->Write();
    hdphi->Write();
    hdphipt->Write();
    for(int i = 0;i<nbins_cent;++i){
    
    hpbpb_TrgObj80[i]->Write();
    hpbpb_TrgObj65[i]->Write();
    hpbpb_TrgObj55[i]->Write();
    
    hpbpb_TrgObjComb[i]->Write();
    hpbpb_TrgObjMB[i]->Write();
    hpbpb_TrgObjMBArray[i]->Write();
    hpbpb_TrgObjMBwoLJ[i]->Write();
    hpbpb_TrgObjMBLJ[i]->Write();
    hEvent_Vz_[i]->Write();
    hdphiptcent[i]->Write();
    hpbpb_TrgObjMBLJ[i]->Write();
    hpbpb_TrgObjMBwoLJ[i]->Write();
    hpbpb_TrgObjMBwoHT[i]->Write();
      
    }

    fout->Close();*/
  
}

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float pi=TMath::Pi();
  
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}

float deltaphi(float phi1, float phi2)
{
  float pi=TMath::Pi();
 
  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}


Double_t GetdPhi(Double_t mphi,Double_t vphi)
{
  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());

  return dphi;
}
