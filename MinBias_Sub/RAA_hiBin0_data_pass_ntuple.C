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

static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};


using namespace std;

const char *ksp="pbpb";

int RAA_hiBin0_data_pass_ntuple(std::string kAlgName="akPu3",
				char* etaWidth = (char*)"20_eta_20",
				Int_t radius = 4)
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // char * etaWidth = (char*)Form("%d_eta_%d",etaLow, etaHigh);
  cout<<"etaWidth = "<<etaWidth<<endl;

  //  bool isSymm = false;
  //  if(etaLow == etaHigh) isSymm = true;
  
  // the cut is a 3 step cut based on the different value of the calopt/pfpt - copy the following lines into your loop (with the corresponding branch address set)
  // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) hGood->Fill();
  // if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) ) hGood->Fill();
  // if(calopt/pfpt > 0.85 & eMax/Sumcand > 0.9) hGood->Fill();
  
  TFile * fData, * fMC; 

  fData = TFile::Open("/mnt/hadoop/cms/store/user/pawan/ntuples/JetRaa_akPu234_PbPb_Data.root");
  fMC = TFile::Open("/mnt/hadoop/cms/store/user/pawan/ntuples/JetRaa_akPu234_PbPb_MC.root");

  TTree * Data_matched= (TTree*)fData->Get(Form("akPu%dJetAnalyzer/matchedJets",radius));
  TTree * Data_unmatched = (TTree*)fData->Get(Form("akPu%dJetAnalyzer/unmatchedPFJets",radius));

  //  TTree * MC_matched = (TTree*)fMC->Get(Form("akPu%dJetAnalyzer/matchedJets",radius));
  //  TTree * MC_unmatched = (TTree*)fMC->Get(Form("akPu%dJetAnalyzer/unmatchedPFJets",radius));
  
  TFile *fout = new TFile("/export/d00/scratch/rkunnawa/rootfiles/RAA/hiBin0_PbPb_data.root","RECREATE");  
  fout->cd();

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("PbPbData : %s %s ",ksp, kAlgName.c_str())<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  fout->mkdir(Form("%sJetAnalyzer",kAlgName.c_str()));
  fout->cd(Form("%sJetAnalyzer"   ,kAlgName.c_str()));
  
  TTree *Data_matched_hiBin0 = new TTree("Data_matched_hiBin0","data matched hiBin0 tree");
  TTree *Data_unmatched_hiBin0 = new TTree("Data_unmatched_hiBin0","data unmatched hiBin0 tree");

  // 1 - Data, 3 - after hiBin selection
  Float_t pfpt_1, pfpt_3;
  Float_t calopt_1, calopt_3;
  Float_t eMax_1, eMax_3;
  Float_t chMax_1, chMax_3;
  Float_t chSum_1, chSum_3;
  Float_t phSum_1, phSum_3;
  Float_t neSum_1, neSum_3;
  Float_t muSum_1, muSum_3;
  Int_t jet55_1, jet65_1, jet80_1;
  Int_t jet55_p_1, jet65_p_1, jet80_p_1;
  Int_t jet55_p_3, jet65_p_3, jet80_p_3;
  Int_t jet55_3, jet65_3, jet80_3;
  Int_t subid_3;
  Int_t hiBin_1, hiBin_3;
  Float_t eta_1, eta_3;
  Float_t phi_1, phi_3;
  Float_t pfrawpt_1, pfrawpt_3; 

  Data_matched->SetBranchAddress("calopt",&calopt_1);
  Data_matched->SetBranchAddress("pfpt",&pfpt_1);
  Data_matched->SetBranchAddress("eMax",&eMax_1);
  Data_matched->SetBranchAddress("chMax",&chMax_1);
  Data_matched->SetBranchAddress("chSum",&chSum_1);
  Data_matched->SetBranchAddress("phSum",&phSum_1);
  Data_matched->SetBranchAddress("neSum",&neSum_1);
  Data_matched->SetBranchAddress("muSum",&muSum_1);
  Data_matched->SetBranchAddress("hiBin",&hiBin_1);
  Data_matched->SetBranchAddress("jet55",&jet55_1);
  Data_matched->SetBranchAddress("jet65",&jet65_1);
  Data_matched->SetBranchAddress("jet80",&jet80_1);
  Data_matched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_matched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_matched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  Data_matched->SetBranchAddress("pfeta",&eta_1);
  Data_matched->SetBranchAddress("pfphi",&phi_1);
  Data_matched->SetBranchAddress("pfrawpt",&pfrawpt_1);
  
  Data_unmatched->SetBranchAddress("pfpt",&pfpt_1);
  Data_unmatched->SetBranchAddress("eMax",&eMax_1);
  Data_unmatched->SetBranchAddress("chMax",&chMax_1);
  Data_unmatched->SetBranchAddress("chSum",&chSum_1);
  Data_unmatched->SetBranchAddress("phSum",&phSum_1);
  Data_unmatched->SetBranchAddress("neSum",&neSum_1);
  Data_unmatched->SetBranchAddress("muSum",&muSum_1);
  Data_unmatched->SetBranchAddress("hiBin",&hiBin_1);
  Data_unmatched->SetBranchAddress("jet55",&jet55_1);
  Data_unmatched->SetBranchAddress("jet65",&jet65_1);
  Data_unmatched->SetBranchAddress("jet80",&jet80_1);
  Data_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_unmatched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_unmatched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  Data_unmatched->SetBranchAddress("pfeta",&eta_1);
  Data_unmatched->SetBranchAddress("pfphi",&phi_1);
  Data_unmatched->SetBranchAddress("pfrawpt",&pfrawpt_1);

  // hiBin0 tree
  Data_matched_hiBin0->Branch("calopt",&calopt_3,"calopt/F");
  Data_matched_hiBin0->Branch("pfpt",&pfpt_3,"pfpt/F");
  Data_matched_hiBin0->Branch("eMax",&eMax_3,"eMax/F");
  Data_matched_hiBin0->Branch("chMax",&chMax_3,"chMax/F");
  Data_matched_hiBin0->Branch("chSum",&chSum_3,"chSum/F");
  Data_matched_hiBin0->Branch("phSum",&phSum_3,"phSum/F");
  Data_matched_hiBin0->Branch("neSum",&neSum_3,"neSum/F");
  Data_matched_hiBin0->Branch("muSum",&muSum_3,"muSum/F");
  Data_matched_hiBin0->Branch("jet55",&jet55_3,"jet55/I");
  Data_matched_hiBin0->Branch("jet65",&jet65_3,"jet65/I");
  Data_matched_hiBin0->Branch("jet80",&jet80_3,"jet80/I");
  Data_matched_hiBin0->Branch("jet55_prescl",&jet55_p_3,"jet55_prescl/I");
  Data_matched_hiBin0->Branch("jet65_prescl",&jet65_p_3,"jet65_prescl/I");
  Data_matched_hiBin0->Branch("jet80_prescl",&jet80_p_3,"jet80_prescl/I");
  Data_matched_hiBin0->Branch("pfeta",&eta_3,"pfeta/F");
  Data_matched_hiBin0->Branch("pfphi",&phi_3,"pfphi/F");
  Data_matched_hiBin0->Branch("pfrawpt",&pfrawpt_3,"pfrawpt/F");
  
  Data_unmatched_hiBin0->Branch("pfpt",&pfpt_3,"pfpt/F");
  Data_unmatched_hiBin0->Branch("eMax",&eMax_3,"eMax/F");
  Data_unmatched_hiBin0->Branch("chMax",&chMax_3,"chMax/F");
  Data_unmatched_hiBin0->Branch("chSum",&chSum_3,"chSum/F");
  Data_unmatched_hiBin0->Branch("phSum",&phSum_3,"phSum/F");
  Data_unmatched_hiBin0->Branch("neSum",&neSum_3,"neSum/F");
  Data_unmatched_hiBin0->Branch("muSum",&muSum_3,"muSum/F");
  Data_unmatched_hiBin0->Branch("jet55",&jet55_3,"jet55/I");
  Data_unmatched_hiBin0->Branch("jet65",&jet65_3,"jet65/I");
  Data_unmatched_hiBin0->Branch("jet80",&jet80_3,"jet80/I");
  Data_unmatched_hiBin0->Branch("jet55_prescl",&jet55_p_3,"jet55_prescl/I");
  Data_unmatched_hiBin0->Branch("jet65_prescl",&jet65_p_3,"jet65_prescl/I");
  Data_unmatched_hiBin0->Branch("jet80_prescl",&jet80_p_3,"jet80_prescl/I");
  Data_unmatched_hiBin0->Branch("pfeta",&eta_3,"pfeta/F");
  Data_unmatched_hiBin0->Branch("pfphi",&phi_3,"pfphi/F");
  Data_unmatched_hiBin0->Branch("pfrawpt",&pfrawpt_3,"pfrawpt/F");
  
  //   MC_matched->SetBranchAddress("calopt",&calopt_3);
  //   MC_matched->SetBranchAddress("pfpt",&pfpt_3);
  //   MC_matched->SetBranchAddress("eMax",&eMax_3);
  //   MC_matched->SetBranchAddress("chMax",&chMax_3);
  //   MC_matched->SetBranchAddress("chSum",&chSum_3);
  //   MC_matched->SetBranchAddress("phSum",&phSum_3);
  //   MC_matched->SetBranchAddress("neSum",&neSum_3);
  //   MC_matched->SetBranchAddress("muSum",&muSum_3);
  //   MC_matched->SetBranchAddress("hiBin",&hiBin_3);
  //   MC_matched->SetBranchAddress("refpt",&pfrefpt_3);
  //   MC_matched->SetBranchAddress("jet55",&jet55_3);
  //   MC_matched->SetBranchAddress("jet65",&jet65_3);
  //   MC_matched->SetBranchAddress("jet80",&jet80_3);
  //   MC_matched->SetBranchAddress("weight", &weight);
  //   MC_matched->SetBranchAddress("subid", &subid_3);
  //   MC_matched->SetBranchAddress("jet55_prescl",&jet55_p_3);
  //   MC_matched->SetBranchAddress("pfeta",&eta_3);
  //   MC_matched->SetBranchAddress("pfrawpt",&pfrawpt_3);
  // 
  //   MC_unmatched->SetBranchAddress("pfpt",&pfpt_3);
  //   MC_unmatched->SetBranchAddress("eMax",&eMax_3);
  //   MC_unmatched->SetBranchAddress("chMax",&chMax_3);
  //   MC_unmatched->SetBranchAddress("chSum",&chSum_3);
  //   MC_unmatched->SetBranchAddress("phSum",&phSum_3);
  //   MC_unmatched->SetBranchAddress("neSum",&neSum_3);
  //   MC_unmatched->SetBranchAddress("muSum",&muSum_3);
  //   MC_unmatched->SetBranchAddress("hiBin",&hiBin_3);
  //   MC_unmatched->SetBranchAddress("refpt",&pfrefpt_3);
  //   MC_unmatched->SetBranchAddress("jet55",&jet55_3);
  //   MC_unmatched->SetBranchAddress("jet65",&jet65_3);
  //   MC_unmatched->SetBranchAddress("jet80",&jet80_3);
  //   MC_unmatched->SetBranchAddress("weight", & weight);
  //   MC_unmatched->SetBranchAddress("subid", &subid_3);
  //   MC_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_3);
  //   MC_unmatched->SetBranchAddress("pfeta",&eta_3);
  //   MC_unmatched->SetBranchAddress("pfrawpt",&pfrawpt_3);
  //   
  cout<<"matched Data ntuple "<<endl;
  long entries = Data_matched->GetEntries();
  //entries = 1000;
  Float_t Jet55_prescl = 2.0475;
  
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    Data_matched->GetEntry(nentry);

    Int_t cBin = findBin(hiBin_1);
    if(cBin == -1 || cBin >= nbins_cent) continue;
    if(hiBin_1==0){
      calopt_3=calopt_1;
      pfpt_3=pfpt_1;
      eMax_3=eMax_1;
      chMax_3=chMax_1;
      chSum_3=chSum_1;
      phSum_3=phSum_1;
      neSum_3=neSum_1;
      muSum_3=muSum_1;
      jet55_3=jet55_1;
      jet65_3=jet65_1;
      jet80_3=jet80_1;
      jet55_p_3=jet55_p_1;
      jet65_p_3=jet65_p_1;
      jet80_p_3=jet80_p_1;
      eta_3=eta_1;
      phi_3=phi_1;
      pfrawpt_3=pfrawpt_1;

      Data_matched_hiBin0->Fill();
      
    }
  }
  // data ntuple loop

  // data unmatched loop:
  entries = Data_unmatched->GetEntries();
  //entries = 1000;
  cout<<"Unmatched Data ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    Data_unmatched->GetEntry(nentry);
    Int_t cBin = findBin(hiBin_1);
    if(cBin == -1 || cBin >= nbins_cent) continue;
    if(hiBin_1==0){
      pfpt_3=pfpt_1;
      eMax_3=eMax_1;
      chMax_3=chMax_1;
      chSum_3=chSum_1;
      phSum_3=phSum_1;
      neSum_3=neSum_1;
      muSum_3=muSum_1;
      jet55_3=jet55_1;
      jet65_3=jet65_1;
      jet80_3=jet80_1;
      jet55_p_3=jet55_p_1;
      jet65_p_3=jet65_p_1;
      jet80_p_3=jet80_p_1;
      eta_3=eta_1;
      phi_3=phi_1;
      pfrawpt_3=pfrawpt_1;
      Data_unmatched_hiBin0->Fill();
    }   
  }  

  Data_matched_hiBin0->Write();
  Data_unmatched_hiBin0->Write();
  //fout->Write();
  fout->Close();   
  return 1;

}
