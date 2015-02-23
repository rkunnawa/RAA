// Raghav Kunnawalkam Elayavalli
// Jan 23rd 2015
// Rutgers
// for questions or comments: raghav.k.e at CERN dot CH

//
// Read the pp MC files. each pthat is split into 4 files (except pthat 540 which has 5 files). Have a TChain for each pthat that gets the branch address set properly. 
//
// similar structure to the read mc macro except it loads a file list for each pthat. similar to whats necessary in the Jet ID workshop! 
//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
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
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

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

static const int nbins_eta = 4;
static const double boundaries_eta[nbins_eta][2] = {
  {0.0,0.5}, {0.5,1.0}, {1.0,1.5}, {1.5,2.0}
};

static const double delta_eta[nbins_eta] = {
  1.0, 1.0, 1.0, 1.0
};

static const char etaWidth [nbins_eta][256] = {
  "0_absEta_05","05_absEta_10","10_absEta_15","15_absEta_20"
};


static const int no_radius = 1;//testing purposes 
static const int list_radius[no_radius] = {5};

// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

static const int N = 5;

class JetData
{
public:
  JetData(string infile, string jetTree, int nofiles) {

    std::ifstream instr(infile.c_str(),std::ifstream::in);

    dir[0] = "hltanalysis";
    dir[1] = "skimanalysis";
    dir[2] = jetTree;
    dir[3] = "hiEvtAnalyzer";
    dir[4] = "hltobject";

    string trees[N] = {
      "HltTree",
      "HltTree",
      "t",
      "HiTree",
      "jetObjTree",
    };
    
    for(int t = 0;t<N;t++){
      jetPP[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
    }//tree loop ends

    for(int i = 0;i<nofiles;i++){
      instr>>filename;
      cout<<"File: "<<filename<<endl;
      for(int t = 0;t<N;t++){
	jetPP[t]->Add(filename.c_str());
      }
    }

    jetPP[2]->AddFriend(jetPP[0]);
    jetPP[2]->AddFriend(jetPP[1]);
    jetPP[2]->AddFriend(jetPP[3]);
    jetPP[2]->AddFriend(jetPP[3]);
    
    jetPP[2]->SetBranchAddress("jtpt" , jtpt );
    jetPP[2]->SetBranchAddress("rawpt", rawpt);
    jetPP[2]->SetBranchAddress("trackMax" , trackMax );
    jetPP[2]->SetBranchAddress("chargedMax",chargedMax);
    jetPP[2]->SetBranchAddress("chargedSum",chargedSum);
    jetPP[2]->SetBranchAddress("neutralMax",neutralMax);
    jetPP[2]->SetBranchAddress("neutralSum",neutralSum);
    jetPP[2]->SetBranchAddress("photonSum",photonSum);
    jetPP[2]->SetBranchAddress("photonMax",photonMax);
    jetPP[2]->SetBranchAddress("eSum",eSum);
    jetPP[2]->SetBranchAddress("eMax",eMax);
    jetPP[2]->SetBranchAddress("muSum",muSum);
    jetPP[2]->SetBranchAddress("muMax",muMax);
    jetPP[2]->SetBranchAddress("refpt", refpt);
    jetPP[2]->SetBranchAddress("nref" ,&njets);
    jetPP[2]->SetBranchAddress("jteta", jteta);
    jetPP[2]->SetBranchAddress("jtphi", jtphi);
    jetPP[2]->SetBranchAddress("jtm", jtmass);
    jetPP[2]->SetBranchAddress("pthat", &pthat);
    jetPP[2]->SetBranchAddress("vz",&vz);
    jetPP[2]->SetBranchAddress("vx",&vx);
    jetPP[2]->SetBranchAddress("vy",&vy);
    jetPP[2]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    jetPP[2]->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
<<<<<<< HEAD
    jetPP[2]->SetBranchAddress("HLT_PAJet40_noJetID_v1",&jet40_1);
=======
    jetPP[2]->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40_1);
    jetPP[2]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p_1);
    jetPP[2]->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60_1);
    jetPP[2]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p_1);
    jetPP[2]->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80_1);
    jetPP[2]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p_1);
>>>>>>> ebe8d1a62d11fe500bdf74e154536b8b51c0ee17

  };

  string filename;
  TChain *jetPP[N];
  string dir[N];
  
  //event varianbles
  int run;
  int evt;
  int lumi;
  float vz;
  float vx;
  float vy;
  float pthat;

  //jet variables 
  float jtpt[1000];
  float rawpt[1000];
  float refpt[1000];
  float jteta[1000];
  float jtphi[1000];
  float jtmass[1000];
  float trackMax[1000];
  float chargedMax[1000];
  float neutralMax[1000];
  float chargedSum[1000];
  float neutralSum[1000];
  float photonSum[1000];
  float eSum[1000];
  float muSum[1000];
  float photonMax[1000];
  float eMax[1000];
  float muMax[1000];
  float genpt[1000];
  
  int njets;
  int pHBHENoiseFilter;
  int pPAcollisionEventSelectionPA;
  int jet40_1;
  int jet60_1;
  int jet80_1;
  int jet40_p_1;
  int jet60_p_1;
  int jet80_p_1;

};

using namespace std;

void RAA_read_mc_pp(char *jet_type="PF", int radius = 5){

  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  cout<<"Running for "<<jet_type<<" jets, R = "<<radius<<endl;
 
  bool printDebug = true;
  TDatime date;
  
  const int nbinsPP_pthat = 11;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  int nofiles_pthat[nbinsPP_pthat];
  string filelistPP_pthat[nbinsPP_pthat];
  Double_t xsectionPP[nbinsPP_pthat+1];

  boundariesPP_pthat[0]=15;
  filelistPP_pthat[0] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat15_forest.txt";
  nofiles_pthat[0] = 4;
  xsectionPP[0]= 0.2034;
  
  boundariesPP_pthat[1]=30;
  filelistPP_pthat[1] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat30_forest.txt";
  nofiles_pthat[1] = 4;
  xsectionPP[1]= 0.01075;
  
  boundariesPP_pthat[2]=50;
  filelistPP_pthat[2] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat50_forest.txt";
  nofiles_pthat[2] = 4;
  xsectionPP[2]= 0.001025;
  
  boundariesPP_pthat[3]=80;
  filelistPP_pthat[3] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat80_forest.txt";
  nofiles_pthat[3] = 4;
  xsectionPP[3]= 9.8650e-05;
  
  boundariesPP_pthat[4]=120;
  filelistPP_pthat[4] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat120_forest.txt";
  nofiles_pthat[4] = 4;
  xsectionPP[4]= 1.1290e-05;

  boundariesPP_pthat[5] = 170;
  filelistPP_pthat[5] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat170_forest.txt";
  nofiles_pthat[5] = 4;
  xsectionPP[5]= 1.4650e-06;
  
  boundariesPP_pthat[6]=220;
  filelistPP_pthat[6] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat220_forest.txt";
  nofiles_pthat[6] = 4;
  xsectionPP[6]= 2.8370e-07;
  
  boundariesPP_pthat[7]=280;
  filelistPP_pthat[7] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat280_forest.txt";
  nofiles_pthat[7] = 4;
  xsectionPP[7]= 5.3230e-08;
  
  boundariesPP_pthat[8]=370;
  filelistPP_pthat[8] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat370_forest.txt";
  nofiles_pthat[8] = 4;
  xsectionPP[8]= 5.9340e-09;
  
  boundariesPP_pthat[9]=460;
  filelistPP_pthat[9] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat460_forest.txt";
  nofiles_pthat[9] = 4;
  xsectionPP[9]= 8.1250e-10;
  
  boundariesPP_pthat[10]=540;
  filelistPP_pthat[10] = "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/jetRAA_pp_mc_pthat540_forest.txt";
  nofiles_pthat[10] = 5;
  xsectionPP[10]= 1.4670e-10;
  
  xsectionPP[11] = 0;
  boundariesPP_pthat[11]=2000;

  // declare the histograms: 
  TH1F *hpp_gen[no_radius][nbins_eta];
  TH1F *hpp_reco[no_radius][nbins_eta];
  TH2F *hpp_matrix[no_radius][nbins_eta];
  TH1F *hpp_JetComb_gen[no_radius][nbins_eta];
  TH1F *hpp_JetComb_reco[no_radius][nbins_eta];
  TH1F *hpp_Jet80_gen[no_radius][nbins_eta];
  TH1F *hpp_Jet80_reco[no_radius][nbins_eta];
  TH1F *hpp_Jet60_gen[no_radius][nbins_eta];
  TH1F *hpp_Jet60_reco[no_radius][nbins_eta];
  TH1F *hpp_Jet40_gen[no_radius][nbins_eta];
  TH1F *hpp_Jet40_reco[no_radius][nbins_eta];
  TH2F *hpp_matrix_HLT[no_radius][nbins_eta];
  TH2F *hpp_mcclosure_matrix[no_radius][nbins_eta];
  TH1F *hpp_mcclosure_data[no_radius][nbins_eta];
  //TH1F *hpp_eta[no_radius][nbins_eta], *hpp_phi[no_radius][nbins_eta];
  TH1F *hpp_eta_full[no_radius], *hpp_phi_full[no_radius];
  TH1F *hpp_eta_full_noScale[no_radius], *hpp_phi_full_noScale[no_radius];
  TH1F *hVzPPMC[no_radius];
  TH1F *hPP_pthat_fine_noScale[no_radius];
  TH1F *hPP_pthat_fine[no_radius];
  TH1F *hPtHatPP[no_radius];
  TH1F *hPtHatRawPP[no_radius];

  for(int k = 0;k<no_radius;k++){
    for(int j = 0;j<nbins_eta;j++){
      hpp_gen[k][j] = new TH1F(Form("hpp_gen_R%d_%s",list_radius[k],etaWidth[j]),Form("gen refpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_reco[k][j] = new TH1F(Form("hpp_reco_R%d_%s",list_radius[k],etaWidth[j]),Form("reco jtpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      
      hpp_JetComb_gen[k][j] = new TH1F(Form("hpp_JetComb_gen_R%d_%s",list_radius[k],etaWidth[j]),Form("gen refpt from JetComb trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_JetComb_reco[k][j] = new TH1F(Form("hpp_JetComb_reco_R%d_%s",list_radius[k],etaWidth[j]),Form("reco jtpt from JetComb trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Jet80_gen[k][j] = new TH1F(Form("hpp_Jet80_gen_R%d_%s",list_radius[k],etaWidth[j]),Form("gen refpt from Jet80 trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Jet80_reco[k][j] = new TH1F(Form("hpp_Jet80_reco_R%d_%s",list_radius[k],etaWidth[j]),Form("reco jtpt from Jet80 trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Jet60_gen[k][j] = new TH1F(Form("hpp_Jet60_gen_R%d_%s",list_radius[k],etaWidth[j]),Form("gen refpt from Jet60 && !Jet80 trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Jet60_reco[k][j] = new TH1F(Form("hpp_Jet60_reco_R%d_%s",list_radius[k],etaWidth[j]),Form("reco jtpt from Jet60 && !jet80  trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Jet40_gen[k][j] = new TH1F(Form("hpp_Jet40_gen_R%d_%s",list_radius[k],etaWidth[j]),Form("gen refpt from Jet40 && !jet65 && !jet55 trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_Jet40_reco[k][j] = new TH1F(Form("hpp_Jet40_reco_R%d_%s",list_radius[k],etaWidth[j]),Form("reco jtpt from Jet40 && !jet65 && !jet55 trigger R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
      hpp_matrix[k][j] = new TH2F(Form("hpp_matrix_R%d_%s",list_radius[k],etaWidth[j]),Form("matrix refpt jtpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      hpp_matrix_HLT[k][j] = new TH2F(Form("hpp_matrix_HLT_R%d_%s",list_radius[k],etaWidth[j]),Form("matrix refpt jtpt from the HLT triggers combined R%d %s",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      hpp_mcclosure_matrix[k][j] = new TH2F(Form("hpp_mcclosure_matrix_R%d_%s",list_radius[k],etaWidth[j]),Form("matrix mcclosure refpt jtpt R%d %s",list_radius[k],etaWidth[j]),1000,0,1000,1000,0,1000);
      //TH2F* hpp_response = new TH2F("hpp_response","response jtpt refpt",1000,0,1000,1000,0,1000);
      hpp_mcclosure_data[k][j] = new TH1F(Form("hpp_mcclosure_data_R%d_%s",list_radius[k],etaWidth[j]),Form("data for unfolding mc closure test pp R%d %s",list_radius[k],etaWidth[j]),1000,0,1000);
    }
    
    hVzPPMC[k] = new TH1F(Form("hVzPPMC_R%d",list_radius[k]),Form("PP MC Vz R%d",list_radius[k]),60,-15,+15);
    hPtHatPP[k] = new TH1F(Form("hPtHatPP_R%d",list_radius[k]),"",nbinsPP_pthat,boundariesPP_pthat);
    hPtHatRawPP[k] = new TH1F(Form("hPtHatRawPP_R%d",list_radius[k]),"",nbinsPP_pthat,boundariesPP_pthat);
    hPP_pthat_fine[k] = new TH1F(Form("hPP_pthat_fine_R%d",list_radius[k]),Form("pp pthat distribution for R=0.%d",list_radius[k]),1000,0,1000);
    hPP_pthat_fine_noScale[k] = new TH1F(Form("hPP_pthat_fine_noScale_R%d",list_radius[k]),Form("PP pthat distribution (unscaled) for R=0.%d",list_radius[k]),1000,0,1000);    
    hpp_eta_full[k] = new TH1F(Form("hpp_eta_full_R%d",list_radius[k]),Form("PP eta distribution for R=%d",list_radius[k]),400,-4,+4);
    hpp_phi_full[k] = new TH1F(Form("hpp_phi_full_R%d",list_radius[k]),Form("PP phi distribution for R=%d",list_radius[k]),400,-4,+4); 
    hpp_eta_full_noScale[k] = new TH1F(Form("hpp_eta_full_noScale_R%d",list_radius[k]),Form("PP eta distribution for noScale R=%d",list_radius[k]),400,-4,+4);
    hpp_phi_full_noScale[k] = new TH1F(Form("hpp_phi_full_noScale_R%d",list_radius[k]),Form("PP phi distribution for noScale R=%d",list_radius[k]),400,-4,+4);

  }

  // setup the jet Data branches:
  JetData *dataPP[no_radius][nbinsPP_pthat];
  for(int k = 0;k<no_radius;k++){
    for(int h = 0;h<nbinsPP_pthat;h++){
      dataPP[k][h] = new JetData(filelistPP_pthat[h],Form("ak%d%sJetAnalyzer",list_radius[k],jet_type),nofiles_pthat[h]);
      TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbinsPP_pthat,boundariesPP_pthat);
      dataPP[k][h]->jetPP[2]->Project("hPtHatTmp","pthat");
      hPtHatRawPP[k]->Add(hPtHatTmp);
      delete hPtHatTmp;
    }
  }
  // Vertex reweighting for pp
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
  
  if(printDebug)hPtHatRawPP[0]->Print("base");
  if(printDebug)cout<<"Filling PP MC"<<endl;

  for(int k = 0;k<no_radius;k++){
    
    for (int h=0;h<nbinsPP_pthat;h++) {
      if (xsectionPP[h]==0) continue;
      if(printDebug)cout <<"Loading PP pthat"<<boundariesPP_pthat[h]<<" sample, cross section = "<<xsectionPP[h]<< Form(" pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[h],boundariesPP_pthat[h+1])<<endl;
     
      //from Pawan's code: /net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/combinePtHatBins/pbpbJEC2014/condor/CondorPbPbCalJec.C
      TEventList *el = new TEventList("el","el");
      double pthat_upper = boundariesPP_pthat[h+1];
      stringstream selection; selection<<"pthat<"<<pthat_upper;
      
      dataPP[k][h]->jetPP[2]->Draw(">>el",selection.str().c_str());
      double fentries = el->GetN();
      if(printDebug)cout<<"tree entries: "<<dataPP[k][h]->jetPP[2]->GetEntries()<<" elist: "<<fentries<<endl;
      delete el;
      
      for (Long64_t jentry=0; jentry<dataPP[k][h]->jetPP[2]->GetEntries();jentry++) {
        dataPP[k][h]->jetPP[2]->GetEntry(jentry);
	dataPP[k][h]->jetPP[2]->GetEntry(jentry);

        int pthatBin = hPtHatPP[k]->FindBin(dataPP[k][h]->pthat);
        float scalepp = (xsectionPP[pthatBin-1]-xsectionPP[pthatBin])/fentries;

        if(fabs(dataPP[k][h]->vz)>15) continue;
        double weight_pt=1;
        double weight_vz=1;
        
	if(!dataPP[k][h]->pPAcollisionEventSelectionPA) continue;
	//if(!dataPP[k][h]->pHBHENoiseFilter) continue;
       
        weight_vz = fVzPP->Eval(dataPP[k][h]->vz);
        //if (weight_vz>5||weight_vz<0.5) cout <<dataPP[k][h]->vz<<" "<<weight_vz<<endl;
        //weight_vz = 1;
	hPP_pthat_fine[k]->Fill(dataPP[k][h]->pthat,scalepp*weight_vz);
	hPP_pthat_fine_noScale[k]->Fill(dataPP[k][h]->pthat);
        hPtHatPP[k]->Fill(dataPP[k][h]->pthat,scalepp*weight_vz);
        int hasLeadingJet = 0;
        hVzPPMC[k]->Fill(dataPP[k][h]->vz,scalepp);
	
        for (int g= 0; g< dataPP[k][h]->njets; g++) { 

	  hpp_eta_full_noScale[k]->Fill(dataPP[k][h]->jteta[g]);
	  hpp_phi_full_noScale[k]->Fill(dataPP[k][h]->jtphi[g]);

	  //if ( dataPP[k][h]->rawpt[g] < 30. ) continue;
	  if ( dataPP[k][h]->jtpt[g] > 2.*dataPP[k][h]->pthat) continue;
	  
	  // jet QA cuts: 
	  if ( dataPP[k][h]->chargedMax[g]/dataPP[k][h]->jtpt[g]<0.05) continue;
	  //if ( dataPP[k][h]->neutralMax[g]/(dataPP[k][h]->chargedMax[g] + dataPP[k][h]->photonMax[g] + dataPP[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( dataPP[k][h]->photonMax[g]/(dataPP[k][h]->chargedMax[g] + dataPP[k][h]->photonMax[g] + dataPP[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( dataPP[k][h]->chargedMax[g]/(dataPP[k][h]->chargedMax[g] + dataPP[k][h]->photonMax[g] + dataPP[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( dataPP[k][h]->muMax[g]/(dataPP[k][h]->chargedMax[g] + dataPP[k][h]->photonMax[g] + dataPP[k][h]->neutralMax[g]) > 0.9 )continue;
	  //if ( dataPP[k][h]->eSum[g]/(dataPP[k][h]->chargedSum[g] + dataPP[k][h]->photonSum[g] + dataPP[k][h]->neutralSum[g] + dataPP[k][h]->muSum[g]) > 0.7 )continue;
	  
	  hpp_eta_full[k]->Fill(dataPP[k][h]->jteta[g],scalepp*weight_vz);
	  hpp_phi_full[k]->Fill(dataPP[k][h]->jtphi[g],scalepp*weight_vz);
	  
	  
          for(int j = 0;j<nbins_eta;j++){
	    
            int subEvt=-1;
	    
            if ( TMath::Abs(dataPP[k][h]->jteta[g])  > boundaries_eta[j][1] || TMath::Abs(dataPP[k][h]->jteta[g]) < boundaries_eta[j][0] ) continue;
	    
            //hpp_response->Fill(dataPP[k][h]->jtpt[k],dataPP[k][h]->refpt[k],scalepp*weight_vz);
            hpp_matrix[k][j]->Fill(dataPP[k][h]->refpt[g],dataPP[k][h]->jtpt[g],scalepp*weight_vz);
            hpp_gen[k][j]->Fill(dataPP[k][h]->refpt[g],scalepp*weight_vz);   
            hpp_reco[k][j]->Fill(dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	    
            if(jentry%2==0) {
	      hpp_mcclosure_data[k][j]->Fill(dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	    }
	    if(jentry%2==1) {
	      hpp_mcclosure_matrix[k][j]->Fill(dataPP[k][h]->refpt[g],dataPP[k][h]->jtpt[g],scalepp*weight_vz);	      
	    }

	    if(dataPP[k][h]->jet80_1){
	      hpp_Jet80_gen[k][j]->Fill(dataPP[k][h]->refpt[g],scalepp*weight_vz);
	      hpp_Jet80_reco[k][j]->Fill(dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	      hpp_matrix_HLT[k][j]->Fill(dataPP[k][h]->refpt[g],dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	    }else if(dataPP[k][h]->jet60_1 && !dataPP[k][h]->jet80_1){
	      hpp_Jet60_gen[k][j]->Fill(dataPP[k][h]->refpt[g],scalepp*weight_vz);
	      hpp_Jet60_reco[k][j]->Fill(dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	      hpp_matrix_HLT[k][j]->Fill(dataPP[k][h]->refpt[g],dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	    }else if(dataPP[k][h]->jet40_1 && !dataPP[k][h]->jet60_1 && !dataPP[k][h]->jet80_1){
	      hpp_Jet40_gen[k][j]->Fill(dataPP[k][h]->refpt[g],scalepp*weight_vz);
	      hpp_Jet40_reco[k][j]->Fill(dataPP[k][h]->jtpt[g],scalepp*weight_vz);
	      hpp_matrix_HLT[k][j]->Fill(dataPP[k][h]->refpt[g],dataPP[k][h]->jtpt[g],dataPP[k][h]->jet40_p_1*scalepp*weight_vz);
	    }
	    
          }//eta loop
	  
        }//njet loop     
	
      }//nentry loop
      
    }//ptbins loop
    
  }//radius loop
  
  TFile f(Form("/export/d00/scratch/rkunnawa/rootfiles/pp_mc_chMaxjtpt_norawpt_absEta_ak%d%s_%d.root",radius,jet_type,date.GetDate()),"RECREATE");
  f.cd();
  
  for(int k = 0;k<no_radius;k++){
    
    for(int j=0;j<nbins_eta;j++){

      divideBinWidth(hpp_gen[k][j]);
      divideBinWidth(hpp_reco[k][j]);
      divideBinWidth(hpp_mcclosure_data[k][j]);

      //hpp_gen[k][j]->Scale(1./(delta_eta[j]));
      //hpp_reco[k][j]->Scale(1./(delta_eta[j]));
      //hpp_mcclosure_data[k][j]->Scale(1./(delta_eta[j]));

      hpp_Jet80_gen[k][j]->Write();
      hpp_Jet80_reco[k][j]->Write();

      hpp_Jet60_gen[k][j]->Write();
      hpp_Jet60_reco[k][j]->Write();

      hpp_Jet40_gen[k][j]->Write();
      hpp_Jet40_reco[k][j]->Write();

      hpp_JetComb_gen[k][j]->Add(hpp_Jet80_gen[k][j]);
      hpp_JetComb_gen[k][j]->Add(hpp_Jet60_gen[k][j]);
      hpp_JetComb_gen[k][j]->Add(hpp_Jet40_gen[k][j]);
      hpp_JetComb_gen[k][j]->Write();

      hpp_JetComb_reco[k][j]->Add(hpp_Jet80_reco[k][j]);
      hpp_JetComb_reco[k][j]->Add(hpp_Jet60_reco[k][j]);
      hpp_JetComb_reco[k][j]->Add(hpp_Jet40_reco[k][j]);
      hpp_JetComb_reco[k][j]->Write();

      hpp_matrix_HLT[k][j]->Write();

      hpp_gen[k][j]->Write();
      if(printDebug)hpp_gen[k][j]->Print("base");
      hpp_reco[k][j]->Write();
      if(printDebug)hpp_reco[k][j]->Print("base");
      hpp_matrix[k][j]->Write();
      if(printDebug)hpp_matrix[k][j]->Print("base");
      hpp_mcclosure_matrix[k][j]->Write();
      if(printDebug)hpp_mcclosure_matrix[k][j]->Print("base");
      hpp_mcclosure_data[k][j]->Write();
      if(printDebug)hpp_mcclosure_data[k][j]->Print("base");

    }
  
    if(printDebug)hPtHatPP[k]->Print("base");
    hPtHatPP[k]->Write();
    if(printDebug)hPP_pthat_fine[k]->Print("base");
    hPP_pthat_fine[k]->Write();
    hpp_eta_full[k]->Write();
    if(printDebug)hpp_eta_full[k]->Print("base");
    hpp_phi_full[k]->Write();
    if(printDebug)hpp_phi_full[k]->Print("base");
    hpp_eta_full_noScale[k]->Write();
    if(printDebug)hpp_eta_full_noScale[k]->Print("base");
    hpp_phi_full_noScale[k]->Write();
    if(printDebug)hpp_phi_full_noScale[k]->Print("base");
    if(printDebug)hPP_pthat_fine_noScale[k]->Print("base");
    hPP_pthat_fine_noScale[k]->Write();
    hVzPPMC[k]->Write();

  }

  f.Write();
  f.Close();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;

}



