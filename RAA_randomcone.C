// Raghav Kunnawalkam Elayavalli
// Oct 6th 2014
// Rutgers
// For questions or comments: raghav.k.e at CERN dot CH

// 
// Macro from Jorge to do the random cone study in the pf candidates 
// one of the main comments is to look at the distribution in centrality and eta with pt on the z axis. 
//



#include <TROOT.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <TVector3.h> 
#include <stdio.h>
#include <string.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TRandom3.h>

//#include "hiForest.h"
//#define PI 3.14159;

//#define nRan 100;

using namespace std;

void RAA_randomcone(int rad=3, const char* jet_type="PF", const char *algo="Vs",const char *type="MC"){
  
  TDatime date;

  TFile *FileA; 
  if(type=="data")
    FileA = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_lumi2_FOREST_TRY2merged/0.root");
  else if(type=="MC")
    FileA = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/HIMC/MB/Track8_Jet26_STARTHI53_LV1/merged2/HiForest_HYDJET_Track8_Jet26_STARTHI53_LV1_merged_forest_0.root");
  
  TFile* outf = new TFile(Form("test_randomcone_%s_ak%s%d%s_%d.root",type,algo,rad,jet_type,date.GetDate()),"recreate"); 
			  
  //TFile *FileA = TFile::Open(Form("/net/hisrv0001/home/icali/hadoop/HIMinBiasUPC_skimmed/MinBias-reTracking-merged/MinBias_Merged_tracking_all.root"));
  //TString outname = "dataAKSkimNtupleRandomConeRings_v4_TkpTCut0_ak3dataMB.root"; 
  //TFile* outf = new TFile(outname,"recreate");

  cout<<"The radius is : "<<rad<<endl;
  char *COLLTYPE = "PbPb";
  cout<<"   COLLTYPE : "<<COLLTYPE<<endl;
  
  //****** Analysis knobs *****************
  Float_t JETRADIUS = (Float_t)rad/10;
  //Float_t PTCUT1 = 120;      //leading jet
  //Float_t PTCUT2 = 50;       //away jet
  Float_t PTCUT3 = 20;       //perp jets
  Float_t ETACUT = 2.0;      //eta acceptance
  //Float_t PHICUT = 0.05;     //dijets phi angle (in PI radians)
  Float_t TRACKPTCUT = 0.0; //
  Float_t TRACKETACUT = 2.0; //eta acceptance
  Float_t PF_TRACKPTCUT = 0.0; //
  Float_t PF_TRACKETACUT = 1.0; //eta acceptance
  Float_t Tower_TRACKPTCUT = 0.0; //
  Float_t Tower_TRACKETACUT = 2.0; //eta acceptance
  //bool debug = true;
  double LOW_dR[] = {   0, 0.05, 0.10, 0.15, 0.20, 0.25};
  double HI_dR[]  = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30};
  const int nRings = 6;      // rings for JetShapes 
  //const Int_t nRan = 100;    // number of random cones
  int nRandom = 100;
  // int nRan= 100; - not used here 
  int Nevents = 1000;
  int nMarks = 50;
  //****************************************
  
  TTree* ak      = (TTree*)FileA->Get(Form("ak%s%d%sJetAnalyzer/t",algo,rad,jet_type));
  // - TTree* ic      = (TTree*)FileA->Get("icPu5JetAnalyzer/t"); // - not there in the new forest 
  TTree* hlt     = (TTree*)FileA->Get("hltanalysis/HltTree");
  TTree* skim    = (TTree*)FileA->Get("skimanalysis/HltTree");
  TTree* hiEv    = (TTree*)FileA->Get("hiEvtAnalyzer/HiTree");
  //TTree* rhEB    = (TTree*)FileA->Get("rechitanalyzer/eb");
  //TTree* rhEE    = (TTree*)FileA->Get("rechitanalyzer/ee");
  //TTree* rhHBHE  = (TTree*)FileA->Get("rechitanalyzer/hbhe");
  TTree* rhtower = (TTree*)FileA->Get("rechitanalyzer/tower");
  //TTree* hcal    = (TTree*)FileA->Get("hcalNoise/hcalNoise"); //- not there in the new forests 
  //if(COLLTYPE=="pA")   TTree* ttrack  = (TTree*)FileA->Get("ppTrack/trackTree");
  //if(COLLTYPE=="PbPb") TTree* ttrack  = (TTree*)FileA->Get("mergedTrack/trackTree"); // - also not needed here 
  TTree* pflow   = (TTree*)FileA->Get("pfcandAnalyzer/pfTree");
  
  // variables for the input file: 

  Float_t vz;
  Int_t nAKJets;
  Float_t jteta[1000];
  Float_t jtphi[1000];
  Float_t jtpt[1000];
  Float_t jtpu[1000];
  Float_t rawpt[1000];
  /*
  Int_t nTrk;
  Float_t trkEta[10000];
  Float_t trkPhi[10000];
  Float_t trkPt[10000];
  Int_t HBHEn;
  Float_t HBHEeta[10000];
  Float_t HBHEphi[10000];
  Float_t HBHEet[10000];
  Int_t EEn;
  Float_t EEeta[10000000];
  Float_t EEphi[10000000];
  Float_t EEet[10000000];
  Int_t EBn;
  Float_t EBeta[10000000];
  Float_t EBphi[10000000];
  Float_t EBet[10000000];
  */

  Int_t PF_n;
  Float_t PF_eta[10000000];
  Float_t PF_phi[10000000];
  Float_t PF_pt[10000000];
  Float_t PF_vspt[10000000];
  Float_t PF_sumpt[10000000];
  Int_t PF_id[10000000];

  Int_t Tower_n;
  Float_t Tower_eta[10000000];
  Float_t Tower_phi[10000000];
  Float_t Tower_pt[10000000];
  Float_t Tower_emEt[10000000];
  Float_t Tower_hadEt[10000000];
  Float_t Tower_sumpt[10000000];
  Float_t Tower_vspt[10000000];

  Int_t eventN;
  Int_t runN;
  Int_t lumiS;
  Int_t cBin;
  Int_t pCES;
  //Int_t phfPosFilter1;
  //Int_t phfNegFilter1;
  //Int_t phltPixelClusterShapeFilter;
  Int_t pprimaryvertexFilter;
  Int_t pHBHENoiseFilter;
  //Int_t HCALmaxhpdhits;
  //Int_t HCALmaxrbxhits;
  //Int_t HCALntrianglenoise;
  //Int_t HCALnspikenoise;
  //Bool_t HCALhasBadRBXTS4TS5;

  //Take the variables from the incoming TTrees
  skim->SetBranchAddress("pcollisionEventSelection",&pCES);
  skim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  
  //if(COLLTYPE=="pA")skim->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
  //if(COLLTYPE=="pA")skim->SetBranchAddress("phltPixelClusterShapeFilter",&phltPixelClusterShapeFilter);
  //if(COLLTYPE=="pA")skim->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
  //if(COLLTYPE=="pA")skim->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
  
  hlt->SetBranchAddress("Run",&runN);
  hlt->SetBranchAddress("Event",&eventN);
  hlt->SetBranchAddress("LumiBlock",&lumiS);
  hiEv->SetBranchAddress("hiBin",&cBin);
  hiEv->SetBranchAddress("vz",&vz);
  
  ak->SetBranchAddress("jteta",&jteta);
  ak->SetBranchAddress("jtpt",&jtpt);
  ak->SetBranchAddress("jtphi",&jtphi);
  ak->SetBranchAddress("jtpu",&jtpu);
  ak->SetBranchAddress("rawpt",&rawpt);
  ak->SetBranchAddress("nref",&nAKJets);

  // ttrack->SetBranchAddress("trkEta",&trkEta);
  // ttrack->SetBranchAddress("trkPt",&trkPt);
  // ttrack->SetBranchAddress("trkPhi",&trkPhi);
  // ttrack->SetBranchAddress("nTrk",&nTrk);
  // rhHBHE->SetBranchAddress("eta",&HBHEeta);
  // rhHBHE->SetBranchAddress("phi",&HBHEphi);
  // rhHBHE->SetBranchAddress("et",&HBHEet);
  // rhHBHE->SetBranchAddress("n",&HBHEn);
  // rhEB->SetBranchAddress("eta",&EBeta);
  // rhEB->SetBranchAddress("phi",&EBphi);
  // rhEB->SetBranchAddress("et",&EBet);
  // rhEB->SetBranchAddress("n",&EBn);
  // rhEE->SetBranchAddress("eta",&EEeta);
  // rhEE->SetBranchAddress("phi",&EEphi);
  // rhEE->SetBranchAddress("et",&EEet);
  // rhEE->SetBranchAddress("n",&EEn);
  // hcal->SetBranchAddress("nspikenoise",&HCALnspikenoise);
  // hcal->SetBranchAddress("hasBadRBXTS4TS5",&HCALhasBadRBXTS4TS5);
  // hcal->SetBranchAddress("maxhpdhits",&HCALmaxhpdhits);
  // hcal->SetBranchAddress("maxrbxhits",&HCALmaxrbxhits);
  // hcal->SetBranchAddress("ntrianglenoise",&HCALntrianglenoise);
  
  pflow->SetBranchAddress("nPFpart",&PF_n);
  pflow->SetBranchAddress("pfEta",&PF_eta);
  pflow->SetBranchAddress("pfPhi",&PF_phi);
  pflow->SetBranchAddress("pfPt",&PF_pt);
  pflow->SetBranchAddress("pfVsPt",&PF_vspt);
  pflow->SetBranchAddress("sumpt",&PF_sumpt);
  pflow->SetBranchAddress("pfId",&PF_id);
  
  rhtower->SetBranchAddress("n",&Tower_n);
  rhtower->SetBranchAddress("eta",&Tower_eta);
  rhtower->SetBranchAddress("phi",&Tower_phi);
  rhtower->SetBranchAddress("et",&Tower_pt);
  rhtower->SetBranchAddress("emEt",&Tower_emEt);
  rhtower->SetBranchAddress("hadEt",&Tower_hadEt);
  rhtower->SetBranchAddress("sumpt",&Tower_sumpt);
  rhtower->SetBranchAddress("vsPt",&Tower_vspt);

  Int_t nJets;
  Float_t jeta[1000];
  Float_t jphi[1000];
  Float_t jpt[1000];
  Float_t jpu[1000];
  Float_t jrawpt[1000];
  Int_t nTracks;
//   Float_t treta[1000];
//   Float_t trphi[1000];
//   Float_t trpt[1000];
  
  // Int_t nTrInJet[1000];
  // Float_t trkSumPtInJet[1000];
  // Float_t dr_TrJet[10000][1000];
  // Int_t nHBHErhits[1000];
  // Float_t HBHEsumEt[1000];
  // Int_t nEErhits[1000];
  // Float_t EEsumEt[1000];
  // Int_t nEBrhits[1000];
  // Float_t EBsumEt[1000];
  
  Int_t nTower_objs[1000];
  Float_t Tower_sumEt[1000];
  Int_t nTower_objs_ring0[1000];
  Int_t nTower_objs_ring1[1000];
  Int_t nTower_objs_ring2[1000];
  Int_t nTower_objs_ring3[1000];
  Int_t nTower_objs_ring4[1000];
  Int_t nTower_objs_ring5[1000];
  Float_t Tower_sumEt_ring0[1000];
  Float_t Tower_sumEt_ring1[1000];
  Float_t Tower_sumEt_ring2[1000];
  Float_t Tower_sumEt_ring3[1000];
  Float_t Tower_sumEt_ring4[1000];
  Float_t Tower_sumEt_ring5[1000];
  // we dont have charge identification in the towers 

  Int_t nPFobjs[1000];
  Float_t PFsumEt[1000];
  Int_t nPFobjs_ring0[1000];
  Int_t nPFobjs_ring1[1000];
  Int_t nPFobjs_ring2[1000];
  Int_t nPFobjs_ring3[1000];
  Int_t nPFobjs_ring4[1000];
  Int_t nPFobjs_ring5[1000];
  Float_t PFsumEt_ring0[1000];
  Float_t PFsumEt_ring1[1000];
  Float_t PFsumEt_ring2[1000];
  Float_t PFsumEt_ring3[1000];
  Float_t PFsumEt_ring4[1000];
  Float_t PFsumEt_ring5[1000];
  Int_t nPFchObjs_ring0[1000];
  Int_t nPFchObjs_ring1[1000];
  Int_t nPFchObjs_ring2[1000];
  Int_t nPFchObjs_ring3[1000];
  Int_t nPFchObjs_ring4[1000];
  Int_t nPFchObjs_ring5[1000];
  Float_t PFchSumEt_ring0[1000];
  Float_t PFchSumEt_ring1[1000];
  Float_t PFchSumEt_ring2[1000];
  Float_t PFchSumEt_ring3[1000];
  Float_t PFchSumEt_ring4[1000];
  Float_t PFchSumEt_ring5[1000];

  // not using tracks in the analysis 
  // Int_t nTks_ring0[10000];
  // Int_t nTks_ring1[10000];
  // Int_t nTks_ring2[10000];
  // Int_t nTks_ring3[10000];
  // Int_t nTks_ring4[10000];
  // Int_t nTks_ring5[10000];
  // Float_t TkSumEt_ring0[10000];
  // Float_t TkSumEt_ring1[10000];
  // Float_t TkSumEt_ring2[10000];
  // Float_t TkSumEt_ring3[10000];
  // Float_t TkSumEt_ring4[10000];
  // Float_t TkSumEt_ring5[10000];
  // Int_t nRanTks_ring0[1000];
  // Int_t nRanTks_ring1[1000];
  // Int_t nRanTks_ring2[1000];
  // Int_t nRanTks_ring3[1000];
  // Int_t nRanTks_ring4[1000];
  // Int_t nRanTks_ring5[1000];
  // Float_t ranTkSumEt_ring0[1000];
  // Float_t ranTkSumEt_ring1[1000];
  // Float_t ranTkSumEt_ring2[1000];
  // Float_t ranTkSumEt_ring3[1000];
  // Float_t ranTkSumEt_ring4[1000];
  // Float_t ranTkSumEt_ring5[1000];

  Int_t run ;
  Int_t event ;
  Int_t lumi ;
  Float_t bin;
  Int_t CES;
//   Float_t dPhi_TT, dPhi_JT, dPhi_JJ;
//   Float_t dEta_TT, dEta_JT, dEta_JJ;
  // Int_t maxHPDhits;
  // Int_t maxRBXhits;
  // Int_t ntrianglenoise;
  // Int_t nspikenoise;
  // Bool_t hasBadRBXTS4TS5;
  // //const Int_t InTracks_set = 3; 
  // Float_t inTreta;
  // Float_t inTrphi;
  // Float_t inTrpt;
  // Float_t outTreta;
  // Float_t outTrphi;
  // Float_t outTrpt;
  // Float_t dRinCone[1000][1000];

 
  Float_t ranConeEta[1000];
  Float_t ranConePhi[1000];
  // Float_t ranTrkSumPt[1000];
  // Int_t ranNtracks[1000];
  // Float_t ranHBHEsumEt[1000];
  // Int_t ranNHBHErhits[1000];
  // Float_t ranEBsumEt[1000];
  // Int_t ranNEBrhits[1000];
  // Float_t ranEEsumEt[1000];
  // Int_t ranNEErhits[1000];

  
  Float_t ranTower_sumEt[1000];
  Int_t ranNTower_objs[1000];
  Int_t ranNTower_objsRing[1000];
  Float_t ranTower_ringSumEt[1000];

  Int_t ranNTower_objs_ring0[1000];
  Float_t ranTower_sumEt_ring0[1000];
  Int_t ranNTower_objs_ring1[1000];
  Float_t ranTower_sumEt_ring1[1000];
  Int_t ranNTower_objs_ring2[1000];
  Float_t ranTower_sumEt_ring2[1000];
  Int_t ranNTower_objs_ring3[1000];
  Float_t ranTower_sumEt_ring3[1000];
  Int_t ranNTower_objs_ring4[1000];
  Float_t ranTower_sumEt_ring4[1000];
  Int_t ranNTower_objs_ring5[1000];
  Float_t ranTower_sumEt_ring5[1000];

  
  Float_t ranPFsumEt[1000];
  Int_t ranNPFobjs[1000];
  Int_t ranNPFobjsRing[1000];
  Float_t ranPFringSumEt[1000];

  Int_t ranNPFobjs_ring0[1000];
  Float_t ranPFsumEt_ring0[1000];
  Int_t ranNPFobjs_ring1[1000];
  Float_t ranPFsumEt_ring1[1000];
  Int_t ranNPFobjs_ring2[1000];
  Float_t ranPFsumEt_ring2[1000];
  Int_t ranNPFobjs_ring3[1000];
  Float_t ranPFsumEt_ring3[1000];
  Int_t ranNPFobjs_ring4[1000];
  Float_t ranPFsumEt_ring4[1000];
  Int_t ranNPFobjs_ring5[1000];
  Float_t ranPFsumEt_ring5[1000];

  Int_t ranNPFchObjs_ring0[1000];
  Float_t ranPFchSumEt_ring0[1000];
  Int_t ranNPFchObjs_ring1[1000];
  Float_t ranPFchSumEt_ring1[1000];
  Int_t ranNPFchObjs_ring2[1000];
  Float_t ranPFchSumEt_ring2[1000];
  Int_t ranNPFchObjs_ring3[1000];
  Float_t ranPFchSumEt_ring3[1000];
  Int_t ranNPFchObjs_ring4[1000];
  Float_t ranPFchSumEt_ring4[1000];
  Int_t ranNPFchObjs_ring5[1000];
  Float_t ranPFchSumEt_ring5[1000];
  
 //Define the branches of the output TTree
  TTree* nt = new TTree("nt","nt");
  nt->Branch("run",&run,"run/I");
  nt->Branch("lumi",&lumi,"lumi/I");
  nt->Branch("event",&event,"event/I");
  nt->Branch("bin", &bin, "bin/F");
  nt->Branch("CES",&CES,"CES/I");
  //nt->Branch("maxHPDhits",&maxHPDhits,"maxHPDhits/I");
  //nt->Branch("maxRBXhits",&maxRBXhits,"maxRBXhits/I");
  //nt->Branch("ntrianglenoise",&ntrianglenoise,"ntrianglenoise/I");
  //nt->Branch("nspikenoise",&nspikenoise,"nspikenoise/I");
  //nt->Branch("hasBadRBXTS4TS5",&hasBadRBXTS4TS5,"hasBadRBXTS4TS5/O");
  nt->Branch("nJets",&nJets,"nJets/I");
  nt->Branch("jeta",jeta,"jeta[nJets]/F");
  nt->Branch("jphi",jphi,"jphi[nJets]/F");
  nt->Branch("jpt",jpt,"jpt[nJets]/F");
  nt->Branch("jpu",jpu,"jpu[nJets]/F");
  nt->Branch("jrawpt",jrawpt,"jrawpt[nJets]/F");
  //nt->Branch("nTracks",&nTracks,"nTracks/I");
  //nt->Branch("treta",treta,"treta[nTracks]/F");
  //nt->Branch("trphi",trphi,"trphi[nTracks]/F");
  //nt->Branch("trpt",trpt,"trpt[nTracks]/F");
  //nt->Branch("nTrInJet",nTrInJet,"nTrInJet[nJets]/I");
  //nt->Branch("trkSumPtInJet",trkSumPtInJet,"trkSumPtInJet[nJets]/F");

  //nt->Branch("dr_TrJet",dr_TrJet,"dr_TrJet[nJets][nTracks]/F");
  //nt->Branch("nHBHErhits",nHBHErhits,"nHBHErhits[nJets]/I");    
  //nt->Branch("HBHEsumEt",HBHEsumEt,"HBHEsumEt[nJets]/F");
  //nt->Branch("nEBrhits",nEBrhits,"nEBrhits[nJets]/I");    
  //nt->Branch("EBsumEt",EBsumEt,"EBsumEt[nJets]/F");
  //nt->Branch("nEErhits",nEErhits,"nEErhits[nJets]/I");    
  //nt->Branch("EEsumEt",EEsumEt,"EEsumEt[nJets]/F");

  if(jet_type=="Calo"){
    nt->Branch("nTower_objs",nTower_objs,"nTower_objs[nJets]/I");    
    nt->Branch("Tower_sumEt",Tower_sumEt,"Tower_sumEt[nJets]/F");

    nt->Branch("nTower_objs_ring0",nTower_objs_ring0,"nTower_objs_ring0[nJets]/I");    
    nt->Branch("Tower_sumEt_ring0",Tower_sumEt_ring0,"Tower_sumEt_ring0[nJets]/F");
    nt->Branch("nTower_objs_ring1",nTower_objs_ring1,"nTower_objs_ring1[nJets]/I");    
    nt->Branch("Tower_sumEt_ring1",Tower_sumEt_ring1,"Tower_sumEt_ring1[nJets]/F");
    nt->Branch("nTower_objs_ring2",nTower_objs_ring2,"nTower_objs_ring2[nJets]/I");    
    nt->Branch("Tower_sumEt_ring2",Tower_sumEt_ring2,"Tower_sumEt_ring2[nJets]/F");
    nt->Branch("nTower_objs_ring3",nTower_objs_ring3,"nTower_objs_ring3[nJets]/I");    
    nt->Branch("Tower_sumEt_ring3",Tower_sumEt_ring3,"Tower_sumEt_ring3[nJets]/F");
    nt->Branch("nTower_objs_ring4",nTower_objs_ring4,"nTower_objs_ring4[nJets]/I");    
    nt->Branch("Tower_sumEt_ring4",Tower_sumEt_ring4,"Tower_sumEt_ring4[nJets]/F");
    nt->Branch("nTower_objs_ring5",nTower_objs_ring5,"nTower_objs_ring5[nJets]/I");    
    nt->Branch("Tower_sumEt_ring5",Tower_sumEt_ring5,"Tower_sumEt_ring5[nJets]/F");
  }else if(jet_type=="PF"){
  
    nt->Branch("nPFobjs",nPFobjs,"nPFobjs[nJets]/I");    
    nt->Branch("PFsumEt",PFsumEt,"PFsumEt[nJets]/F");

    nt->Branch("nPFobjs_ring0",nPFobjs_ring0,"nPFobjs_ring0[nJets]/I");    
    nt->Branch("PFsumEt_ring0",PFsumEt_ring0,"PFsumEt_ring0[nJets]/F");
    nt->Branch("nPFobjs_ring1",nPFobjs_ring1,"nPFobjs_ring1[nJets]/I");    
    nt->Branch("PFsumEt_ring1",PFsumEt_ring1,"PFsumEt_ring1[nJets]/F");
    nt->Branch("nPFobjs_ring2",nPFobjs_ring2,"nPFobjs_ring2[nJets]/I");    
    nt->Branch("PFsumEt_ring2",PFsumEt_ring2,"PFsumEt_ring2[nJets]/F");
    nt->Branch("nPFobjs_ring3",nPFobjs_ring3,"nPFobjs_ring3[nJets]/I");    
    nt->Branch("PFsumEt_ring3",PFsumEt_ring3,"PFsumEt_ring3[nJets]/F");
    nt->Branch("nPFobjs_ring4",nPFobjs_ring4,"nPFobjs_ring4[nJets]/I");    
    nt->Branch("PFsumEt_ring4",PFsumEt_ring4,"PFsumEt_ring4[nJets]/F");
    nt->Branch("nPFobjs_ring5",nPFobjs_ring5,"nPFobjs_ring5[nJets]/I");    
    nt->Branch("PFsumEt_ring5",PFsumEt_ring5,"PFsumEt_ring5[nJets]/F");

    nt->Branch("nPFchObjs_ring0",nPFchObjs_ring0,"nPFchObjs_ring0[nJets]/I");    
    nt->Branch("PFchSumEt_ring0",PFchSumEt_ring0,"PFchSumEt_ring0[nJets]/F");
    nt->Branch("nPFchObjs_ring1",nPFchObjs_ring1,"nPFchObjs_ring1[nJets]/I");    
    nt->Branch("PFchSumEt_ring1",PFchSumEt_ring1,"PFchSumEt_ring1[nJets]/F");
    nt->Branch("nPFchObjs_ring2",nPFchObjs_ring2,"nPFchObjs_ring2[nJets]/I");    
    nt->Branch("PFchSumEt_ring2",PFchSumEt_ring2,"PFchSumEt_ring2[nJets]/F");
    nt->Branch("nPFchObjs_ring3",nPFchObjs_ring3,"nPFchObjs_ring3[nJets]/I");    
    nt->Branch("PFchSumEt_ring3",PFchSumEt_ring3,"PFchSumEt_ring3[nJets]/F");
    nt->Branch("nPFchObjs_ring4",nPFchObjs_ring4,"nPFchObjs_ring4[nJets]/I");    
    nt->Branch("PFchSumEt_ring4",PFchSumEt_ring4,"PFchSumEt_ring4[nJets]/F");
    nt->Branch("nPFchObjs_ring5",nPFchObjs_ring5,"nPFchObjs_ring5[nJets]/I");    
    nt->Branch("PFchSumEt_ring5",PFchSumEt_ring5,"PFchSumEt_ring5[nJets]/F");
  }
  /*
  nt->Branch("nTks_ring0",nTks_ring0,"nTks_ring0[nJets]/I");    
  nt->Branch("TkSumEt_ring0",TkSumEt_ring0,"TkSumEt_ring0[nJets]/F");
  nt->Branch("nTks_ring1",nTks_ring1,"nTks_ring1[nJets]/I");    
  nt->Branch("TkSumEt_ring1",TkSumEt_ring1,"TkSumEt_ring1[nJets]/F");
  nt->Branch("nTks_ring2",nTks_ring2,"nTks_ring2[nJets]/I");    
  nt->Branch("TkSumEt_ring2",TkSumEt_ring2,"TkSumEt_ring2[nJets]/F");
  nt->Branch("nTks_ring3",nTks_ring3,"nTks_ring3[nJets]/I");    
  nt->Branch("TkSumEt_ring3",TkSumEt_ring3,"TkSumEt_ring3[nJets]/F");
  nt->Branch("nTks_ring4",nTks_ring4,"nTks_ring4[nJets]/I");    
  nt->Branch("TkSumEt_ring4",TkSumEt_ring4,"TkSumEt_ring4[nJets]/F");
  nt->Branch("nTks_ring5",nTks_ring5,"nTks_ring5[nJets]/I");    
  nt->Branch("TkSumEt_ring5",TkSumEt_ring5,"TkSumEt_ring5[nJets]/F");


  nt->Branch("nRanTks_ring0",nRanTks_ring0,"nRanTks_ring0[200]/I");    
  nt->Branch("ranTkSumEt_ring0",ranTkSumEt_ring0,"ranTkSumEt_ring0[200]/F");
  nt->Branch("nRanTks_ring1",nRanTks_ring1,"nRanTks_ring1[200]/I");    
  nt->Branch("ranTkSumEt_ring1",ranTkSumEt_ring1,"ranTkSumEt_ring1[200]/F");
  nt->Branch("nRanTks_ring2",nRanTks_ring2,"nRanTks_ring2[200]/I");    
  nt->Branch("ranTkSumEt_ring2",ranTkSumEt_ring2,"ranTkSumEt_ring2[200]/F");
  nt->Branch("nRanTks_ring3",nRanTks_ring3,"nRanTks_ring3[200]/I");    
  nt->Branch("ranTkSumEt_ring3",ranTkSumEt_ring3,"ranTkSumEt_ring3[200]/F");
  nt->Branch("nRanTks_ring4",nRanTks_ring4,"nRanTks_ring4[200]/I");    
  nt->Branch("ranTkSumEt_ring4",ranTkSumEt_ring4,"ranTkSumEt_ring4[200]/F");
  nt->Branch("nRanTks_ring5",nRanTks_ring5,"nRanTks_ring5[200]/I");    
  nt->Branch("ranTkSumEt_ring5",ranTkSumEt_ring5,"ranTkSumEt_ring5[200]/F");
//   nt->Branch("dPhi_JJ",&dPhi_JJ,"dPhi_JJ/F");
//   nt->Branch("dPhi_JT",&dPhi_JT,"dPhi_JT/F");
//   nt->Branch("dPhi_TT",&dPhi_TT,"dPhi_TT/F");
//   nt->Branch("dEta_JJ",&dEta_JJ,"dEta_JJ/F");
//   nt->Branch("dEta_JT",&dEta_JT,"dEta_JT/F");
//   nt->Branch("dEta_TT",&dEta_TT,"dEta_TT/F");
  nt->Branch("inTreta",&inTreta,"inTreta/F");
  nt->Branch("inTrphi",&inTrphi,"inTrphi/F");
  nt->Branch("inTrpt",&inTrpt,"inTrpt/F");
  nt->Branch("outTreta",&outTreta,"outTreta/F");
  nt->Branch("outTrphi",&outTrphi,"outTrphi/F");
  nt->Branch("outTrpt",&outTrpt,"outTrpt/F");
  nt->Branch("dRinCone",dRinCone,"dRinCone[nJets][1000]");
  */

  nt->Branch("ranConeEta",ranConeEta,"ranConeEta[200]/F");
  nt->Branch("ranConePhi",ranConePhi,"ranConePhi[200]/F");
  //nt->Branch("ranTrkSumPt",ranTrkSumPt,"ranTrkSumPt[200]/F");
  //nt->Branch("ranNtracks",ranNtracks,"ranNtracks[200]/I");
  /*
  nt->Branch("ranHBHEsumEt",ranHBHEsumEt,"ranHBHEsumEt[200]/F");
  nt->Branch("ranNHBHErhits",ranNHBHErhits,"ranNHBHErhits[200]/I");
  nt->Branch("ranEBsumEt",ranEBsumEt,"ranEBsumEt[200]/F");
  nt->Branch("ranNEBrhits",ranNEBrhits,"ranNEBrhits[200]/I");
  nt->Branch("ranEEsumEt",ranEEsumEt,"ranEEsumEt[200]/F");
  nt->Branch("ranNEErhits",ranNEErhits,"ranNEErhits[200]/I");
  */

  if(jet_type=="Calo"){
    nt->Branch("ranTower_sumEt",ranTower_sumEt,"ranTower_sumEt[200]/F");
    nt->Branch("ranNTower_objs",ranNTower_objs,"ranNTower_objs[200]/I");
  
    nt->Branch("ranNTower_objs_ring0",ranNTower_objs_ring0,"ranNTower_objs_ring0[200]/I");    
    nt->Branch("ranTower_sumEt_ring0",ranTower_sumEt_ring0,"ranTower_sumEt_ring0[200]/F");
    nt->Branch("ranNTower_objs_ring1",ranNTower_objs_ring1,"ranNTower_objs_ring1[200]/I");    
    nt->Branch("ranTower_sumEt_ring1",ranTower_sumEt_ring1,"ranTower_sumEt_ring1[200]/F");
    nt->Branch("ranNTower_objs_ring2",ranNTower_objs_ring2,"ranNTower_objs_ring2[200]/I");    
    nt->Branch("ranTower_sumEt_ring2",ranTower_sumEt_ring2,"ranTower_sumEt_ring2[200]/F");
    nt->Branch("ranNTower_objs_ring3",ranNTower_objs_ring3,"ranNTower_objs_ring3[200]/I");    
    nt->Branch("ranTower_sumEt_ring3",ranTower_sumEt_ring3,"ranTower_sumEt_ring3[200]/F");
    nt->Branch("ranNTower_objs_ring4",ranNTower_objs_ring4,"ranNTower_objs_ring4[200]/I");    
    nt->Branch("ranTower_sumEt_ring4",ranTower_sumEt_ring4,"ranTower_sumEt_ring4[200]/F");
    nt->Branch("ranNTower_objs_ring5",ranNTower_objs_ring5,"ranNTower_objs_ring5[200]/I");    
    nt->Branch("ranTower_sumEt_ring5",ranTower_sumEt_ring5,"ranTower_sumEt_ring5[200]/F");
  }else if(jet_type=="PF"){
    nt->Branch("ranPFsumEt",ranPFsumEt,"ranPFsumEt[200]/F");
    nt->Branch("ranNPFobjs",ranNPFobjs,"ranNPFobjs[200]/I");
  
    nt->Branch("ranNPFobjs_ring0",ranNPFobjs_ring0,"ranNPFobjs_ring0[200]/I");    
    nt->Branch("ranPFsumEt_ring0",ranPFsumEt_ring0,"ranPFsumEt_ring0[200]/F");
    nt->Branch("ranNPFobjs_ring1",ranNPFobjs_ring1,"ranNPFobjs_ring1[200]/I");    
    nt->Branch("ranPFsumEt_ring1",ranPFsumEt_ring1,"ranPFsumEt_ring1[200]/F");
    nt->Branch("ranNPFobjs_ring2",ranNPFobjs_ring2,"ranNPFobjs_ring2[200]/I");    
    nt->Branch("ranPFsumEt_ring2",ranPFsumEt_ring2,"ranPFsumEt_ring2[200]/F");
    nt->Branch("ranNPFobjs_ring3",ranNPFobjs_ring3,"ranNPFobjs_ring3[200]/I");    
    nt->Branch("ranPFsumEt_ring3",ranPFsumEt_ring3,"ranPFsumEt_ring3[200]/F");
    nt->Branch("ranNPFobjs_ring4",ranNPFobjs_ring4,"ranNPFobjs_ring4[200]/I");    
    nt->Branch("ranPFsumEt_ring4",ranPFsumEt_ring4,"ranPFsumEt_ring4[200]/F");
    nt->Branch("ranNPFobjs_ring5",ranNPFobjs_ring5,"ranNPFobjs_ring5[200]/I");    
    nt->Branch("ranPFsumEt_ring5",ranPFsumEt_ring5,"ranPFsumEt_ring5[200]/F");

    nt->Branch("ranNPFchObjs_ring0",ranNPFchObjs_ring0,"ranNPFchObjs_ring0[200]/I");    
    nt->Branch("ranPFchSumEt_ring0",ranPFchSumEt_ring0,"ranPFchSumEt_ring0[200]/F");
    nt->Branch("ranNPFchObjs_ring1",ranNPFchObjs_ring1,"ranNPFchObjs_ring1[200]/I");    
    nt->Branch("ranPFchSumEt_ring1",ranPFchSumEt_ring1,"ranPFchSumEt_ring1[200]/F");
    nt->Branch("ranNPFchObjs_ring2",ranNPFchObjs_ring2,"ranNPFchObjs_ring2[200]/I");    
    nt->Branch("ranPFchSumEt_ring2",ranPFchSumEt_ring2,"ranPFchSumEt_ring2[200]/F");
    nt->Branch("ranNPFchObjs_ring3",ranNPFchObjs_ring3,"ranNPFchObjs_ring3[200]/I");    
    nt->Branch("ranPFchSumEt_ring3",ranPFchSumEt_ring3,"ranPFchSumEt_ring3[200]/F");
    nt->Branch("ranNPFchObjs_ring4",ranNPFchObjs_ring4,"ranNPFchObjs_ring4[200]/I");    
    nt->Branch("ranPFchSumEt_ring4",ranPFchSumEt_ring4,"ranPFchSumEt_ring4[200]/F");
    nt->Branch("ranNPFchObjs_ring5",ranNPFchObjs_ring5,"ranNPFchObjs_ring5[200]/I");    
    nt->Branch("ranPFchSumEt_ring5",ranPFchSumEt_ring5,"ranPFchSumEt_ring5[200]/F");
  }
  // Int_t ranNPFchObjs_ring0[1000];
  // Float_t ranPFchSumEt_ring0[1000];

  TRandom3 myDice(314);

  cout<<"Right before the loop"<<endl;
  //int Nevents = ak->GetEntries();

  for(int iev = 0; iev < Nevents; ++iev)//loop over all the entries (events)
    {
 
      if (iev%nMarks==0) 
	{
	cout<<"event: "<<iev<<" / "<<Nevents<<endl;
	} 
      
      //----------------------------------------------------------------------
      //Initialize some stuff and get the actual event
      //----------------------------------------------------------------------
      nJets = 0;
      hlt->GetEntry(iev);
      ak->GetEntry(iev);
      // ic->GetEntry(iev);
      skim->GetEntry(iev);
      hiEv->GetEntry(iev);
      // rhHBHE->GetEntry(iev);
      // rhEE->GetEntry(iev);
      // rhEB->GetEntry(iev);
      // hcal->GetEntry(iev);
      // ttrack->GetEntry(iev);
      rhtower->GetEntry(iev);
      pflow->GetEntry(iev);
      //----------------------------------------------------------------------
      // Event variables, and HCAL noise cleaning
      //----------------------------------------------------------------------
      
      CES = pCES;
      if (CES==0 || pHBHENoiseFilter==0 || fabs(vz)>15) continue; // since this is only for PbPb
      
      //if ( (COLLTYPE=="PbPb") && (CES==0)  ) continue;
      
      //if ( (COLLTYPE=="pA") && (pHBHENoiseFilter==0 ||  phltPixelClusterShapeFilter==0 || pprimaryvertexFilter == 0|| pHBHENoiseFilter==0 || phfPosFilter1==0 || phfNegFilter1==0)) continue;
      /*
      if (COLLTYPE=="pA"){
      	if (fabs(vz) > 15)continue; 
      	if (TYPE=="DATA"){
      	  if (runN!=202792) continue;
      	  if (lumiS<317 || lumiS>1014) continue; 
      	}
      }

      */
      
      run = runN;
      event = eventN;
      lumi = lumiS;
      bin = cBin;
      
      //cout<<"run = "<<run<<endl;
      //cout<<"event = "<<event<<endl;
      
      // maxHPDhits = HCALmaxhpdhits;
      // maxRBXhits = HCALmaxrbxhits;
      // ntrianglenoise = HCALntrianglenoise;
      // nspikenoise = HCALnspikenoise;
      // hasBadRBXTS4TS5 = HCALhasBadRBXTS4TS5;
      
      nJets = nAKJets;
      // nTracks = nTrk;
      //nJets = 1; //for dijet studies
      for (int t=0; t<=nJets; t++){//Initialize and clear arrays
	jpt[t] = -99;
	jeta[t] = -99;
	jphi[t] = -99;
	jpu[t] = -99;
	jrawpt[t] = -99;
	//nTrInJet[t] = 0;
	//trkSumPtInJet[t] = -99;
	// HBHEsumEt[t] = -99;
 	// nHBHErhits[t] = 0;
 	// EBsumEt[t] = -99;
 	// nEBrhits[t] = 0;
 	// EEsumEt[t] = -99;
 	// nEErhits[t] = 0;
	
	nTower_objs[t] = 0;
	Tower_sumEt[t] =-99;
	Tower_sumEt_ring0[t] = -99;
	nTower_objs_ring0[t] = 0;
	Tower_sumEt_ring1[t] = -99;
	nTower_objs_ring1[t] = 0;
	Tower_sumEt_ring2[t] = -99;
	nTower_objs_ring2[t] = 0;
	Tower_sumEt_ring3[t] = -99;
	nTower_objs_ring3[t] = 0;
	Tower_sumEt_ring4[t] = -99;
	nTower_objs_ring4[t] = 0;
	Tower_sumEt_ring5[t] = -99;
	nTower_objs_ring5[t] = 0;
	
	
	nPFobjs[t] = 0;
	PFsumEt[t] =-99;
	PFsumEt_ring0[t] = -99;
	nPFobjs_ring0[t] = 0;
	PFsumEt_ring1[t] = -99;
	nPFobjs_ring1[t] = 0;
	PFsumEt_ring2[t] = -99;
	nPFobjs_ring2[t] = 0;
	PFsumEt_ring3[t] = -99;
	nPFobjs_ring3[t] = 0;
	PFsumEt_ring4[t] = -99;
	nPFobjs_ring4[t] = 0;
	PFsumEt_ring5[t] = -99;
	nPFobjs_ring5[t] = 0;
	PFchSumEt_ring0[t] = -99;
	nPFchObjs_ring0[t] = 0;
	PFchSumEt_ring1[t] = -99;
	nPFchObjs_ring1[t] = 0;
	PFchSumEt_ring2[t] = -99;
	nPFchObjs_ring2[t] = 0;
	PFchSumEt_ring3[t] = -99;
	nPFchObjs_ring3[t] = 0;
	PFchSumEt_ring4[t] = -99;
	nPFchObjs_ring4[t] = 0;
	PFchSumEt_ring5[t] = -99;
	nPFchObjs_ring5[t] = 0;
	// TkSumEt_ring0[t] = -99;
	// nTks_ring0[t] = 0;
	// TkSumEt_ring1[t] = -99;
	// nTks_ring1[t] = 0;
	// TkSumEt_ring2[t] = -99;
	// nTks_ring2[t] = 0;
	// TkSumEt_ring3[t] = -99;
	// nTks_ring3[t] = 0;
	// TkSumEt_ring4[t] = -99;
	// nTks_ring4[t] = 0;
	// TkSumEt_ring5[t] = -99;
	// nTks_ring5[t] = 0;
      }
   
      for (int qq =0; qq<nRandom; qq++){
	ranConeEta[qq] =-99;
	ranConePhi[qq] =-99;
	// ranTrkSumPt[qq] =-99;
	// ranNtracks[qq] =0;
 	// ranHBHEsumEt[qq] =-99;
 	// ranNHBHErhits[qq] =0;
 	// ranEBsumEt[qq] =-99;
 	// ranNEBrhits[qq] =0;
 	// ranEEsumEt[qq] =-99;
 	// ranNEErhits[qq] =0;

	ranNTower_objs[qq] = 0;
	ranTower_sumEt[qq] =-99;
	ranTower_sumEt_ring0[qq] = -99;
	ranNTower_objs_ring0[qq] = 0;
	ranTower_sumEt_ring1[qq] = -99;
	ranNTower_objs_ring1[qq] = 0;
	ranTower_sumEt_ring2[qq] = -99;
	ranNTower_objs_ring2[qq] = 0;
	ranTower_sumEt_ring3[qq] = -99;
	ranNTower_objs_ring3[qq] = 0;
	ranTower_sumEt_ring4[qq] = -99;
	ranNTower_objs_ring4[qq] = 0;
	ranTower_sumEt_ring5[qq] = -99;
	ranNTower_objs_ring5[qq] = 0;
	
	ranNPFobjs[qq] = 0;
	ranPFsumEt[qq] =-99;
	ranPFsumEt_ring0[qq] = -99;
	ranNPFobjs_ring0[qq] = 0;
	ranPFsumEt_ring1[qq] = -99;
	ranNPFobjs_ring1[qq] = 0;
	ranPFsumEt_ring2[qq] = -99;
	ranNPFobjs_ring2[qq] = 0;
	ranPFsumEt_ring3[qq] = -99;
	ranNPFobjs_ring3[qq] = 0;
	ranPFsumEt_ring4[qq] = -99;
	ranNPFobjs_ring4[qq] = 0;
	ranPFsumEt_ring5[qq] = -99;
	ranNPFobjs_ring5[qq] = 0;
	ranPFchSumEt_ring0[qq] = -99;
	ranNPFchObjs_ring0[qq] = 0;
	ranPFchSumEt_ring1[qq] = -99;
	ranNPFchObjs_ring1[qq] = 0;
	ranPFchSumEt_ring2[qq] = -99;
	ranNPFchObjs_ring2[qq] = 0;
	ranPFchSumEt_ring3[qq] = -99;
	ranNPFchObjs_ring3[qq] = 0;
	ranPFchSumEt_ring4[qq] = -99;
	ranNPFchObjs_ring4[qq] = 0;
	ranPFchSumEt_ring5[qq] = -99;
	ranNPFchObjs_ring5[qq] = 0;
	// ranTkSumEt_ring0[qq] = -99;
	// nRanTks_ring0[qq] = 0;
	// ranTkSumEt_ring1[qq] = -99;
	// nRanTks_ring1[qq] = 0;
	// ranTkSumEt_ring2[qq] = -99;
	// nRanTks_ring2[qq] = 0;
	// ranTkSumEt_ring3[qq] = -99;
	// nRanTks_ring3[qq] = 0;
	// ranTkSumEt_ring4[qq] = -99;
	// nRanTks_ring4[qq] = 0;
	// ranTkSumEt_ring5[qq] = -99;
	// nRanTks_ring5[qq] = 0;
      }
      
      
      //----------------------------------------------------------------------
      //  Now loop over the Jets
      //----------------------------------------------------------------------
      for (int iJet=0; iJet<nJets; iJet++)
	{ 
	  jeta[iJet]   = jteta[iJet];
	  jpt[iJet]    = jtpt[iJet]; 
	  jphi[iJet]   = jtphi[iJet];
	  jpu[iJet]    = jtpu[iJet];
	  jrawpt[iJet] = rawpt[iJet];
	  //cout<<"HBEH eta: "<<HBHEeta[iJet]<<"EEeta: "<<EEeta[iJet]<<endl;
	  //----------------------------------------------------------------------
	  // Use kinematic cuts
	  //----------------------------------------------------------------------
	  if ( jteta[iJet] > fabs(ETACUT) ) continue;
	  if ( jtpt[iJet] < PTCUT3 ) continue;

	  
	  //cout<<"jet pt = "<<jpt[iJet]<<endl;
	  //cout<<"jet pu = "<<jpu[iJet]<<endl;
	  //cout<<"raw pt = "<<jrawpt[iJet]<<endl;
	  //cout<<"eta: "<<jeta[iJet]<<", phi: "<<jphi[iJet]<<endl;

	  /*

	  //----------------------------------------------------------------------
	  //  Now loop over the tracks, and get the ones inside the jet
	  //----------------------------------------------------------------------
	  Int_t nTrkInJetCounter=0;
	  Float_t trkSumPt = 0;
	  double maxPtIn = 0; //to identify the track with the largest pT, inside the cone
	  double maxPtOut = 0;//to identify the track with the largest pT, outside the cone
	  int iOut = 0;
	  int iIn = 0;

	  vector<int> iInTracks;
	  vector<int> iOutTracks;
	  double ringSumEtTk[nRings] ={0,0,0,0,0,0};
	  Int_t ringCounterTk[nRings] = {0,0,0,0,0,0};
	  for (int iTrk = 0; iTrk<nTracks ; iTrk++)
	    {
	      if (trkPt[iTrk] < TRACKPTCUT) continue;
	      //float iieta = trkEta[iTrk];
	      if (fabs(trkEta[iTrk]) > TRACKETACUT ) continue;
	      double dr_tj = sqrt(((trkEta[iTrk]-jteta[iJet])*(trkEta[iTrk]-jteta[iJet])) + (((trkPhi[iTrk]-jtphi[iJet])*(trkPhi[iTrk]-jtphi[iJet]))));
	      //cout<<"the radius of theis track to the center of teh jet is :"<<dr_tj1<<endl; 
	      //dr_TrJet[iJet][iTrk] = dr_tj;
	      if (dr_tj <= JETRADIUS )
		{
		  nTrkInJetCounter +=1;
		  trkSumPt += trkPt[iTrk];
		  //if(iJet==0)cout<<"tracks around(in) this jeta have pt: "<<trkPt[iTrk]<<endl;
		}
	      dRinCone[iJet][iTrk] = dr_tj;
	      //For the JetShapes analysis get the energy in each ring in dR(\eta,\phi)
	      for (int iRing = 0; iRing<=nRings; iRing++)
		{
		  if (dr_tj> LOW_dR[iRing] && dr_tj<= HI_dR[iRing])
		    {
		      ringSumEtTk[iRing] += trkPt[iTrk];
		      ringCounterTk[iRing] +=1;

		    } 
		}	  
	      
	      //Now gather up the sums for each jet
	      
	      TkSumEt_ring0[iTrk] = ringSumEtTk[0];
	      TkSumEt_ring1[iTrk] = ringSumEtTk[1];
	      TkSumEt_ring2[iTrk] = ringSumEtTk[2];
	      TkSumEt_ring3[iTrk] = ringSumEtTk[3];
	      TkSumEt_ring4[iTrk] = ringSumEtTk[4];
	      TkSumEt_ring5[iTrk] = ringSumEtTk[5];
	      
	      nTks_ring0[iTrk] = ringCounterTk[0];
	      nTks_ring1[iTrk] = ringCounterTk[1];
	      nTks_ring2[iTrk] = ringCounterTk[2];
	      nTks_ring3[iTrk] = ringCounterTk[3];
	      nTks_ring4[iTrk] = ringCounterTk[4];
	      nTks_ring5[iTrk] = ringCounterTk[5];
	      //----------------------------------------------------------------------
	      //  Now, some correlations:: leadingJet-subLeadingJet, 
	      //  leadingJet-leadingAwayTrack, leadingInsideJetTrack-leadingAwayTrack	
	      //----------------------------------------------------------------------
	      if (trkPt[iTrk] > maxPtOut && dr_tj > JETRADIUS ) {
		//if current track has a gretaer pT that previous and is inside the 
		//jet-cone radius re-value the maxpT and keep the iterator value
		maxPtOut = trkPt[iTrk];
		iOut = iTrk;
		
	      } 
	      if (trkPt[iTrk] > maxPtIn && dr_tj <= JETRADIUS ) {
		maxPtIn = trkPt[iTrk];	  
		iIn =iTrk;
	      }
	    }

	  if (iJet==0){// catch only the correlations to the leading jet
	    //cout<<"**** and i kept the pt : "<<trkPt[iIn]<<"  element: "<<iIn<<endl;
	    inTreta = trkEta[iIn];
	    inTrphi = trkPhi[iIn];
	    inTrpt = trkPt[iIn];
	    outTreta = trkEta[iOut];
	    outTrphi = trkPhi[iOut];
	    outTrpt = trkPt[iOut];
	  }

	  nTrInJet[iJet] = nTrkInJetCounter;
	  trkSumPtInJet[iJet] = trkSumPt;
	  //cout<<"event: "<<iev<<" Jet: "<<iJet<<" with eta: "<<jteta[iJet]<<" and "<<nTrkInJetCounter<<" tracks"<<endl;	
	  if (nTrkInJetCounter>2500) cout<<"Alert: event: "<<iev<<" Jet: "<<iJet<<" with eta: "<<jteta[iJet]<<" and "<<nTrkInJetCounter<<" tracks"<<endl;
	
	  //
	  */

	  if(jet_type=="PF"){
	    //----------------------------------------------------------------------
	    //  Now loop over the Pflow candidates, and get the ones inside the jet	
	    //----------------------------------------------------------------------
	    double sumEtPF=0;
	    Int_t PFCounter = 0;	  
	    double ringSumEt[nRings]={0,0,0,0,0,0};
	    Int_t ringCounter[nRings] ={0,0,0,0,0,0};
	    double ringSumEtCh[nRings]={0,0,0,0,0,0};
	    Int_t ringCounterCh[nRings] ={0,0,0,0,0,0};	  
	    //cout<<"no of PF candidates = "<<PF_n<<endl;
	    for(int iPF = 0; iPF<PF_n; iPF++)
	      {
		if (PF_pt[iPF] < PF_TRACKPTCUT) continue;
		if ( fabs(PF_eta[iPF] ) >= PF_TRACKETACUT ) continue;
		double dRhitJetPF = sqrt(((PF_eta[iPF]-jteta[iJet])*(PF_eta[iPF]-jteta[iJet])) + ((PF_phi[iPF]-jtphi[iJet])*(PF_phi[iPF]-jtphi[iJet])));  
		if (dRhitJetPF<=JETRADIUS)
		  {
		    sumEtPF += PF_pt[iPF];
		    PFCounter +=1;
		    //cout<<"PF pt ["<<iPF<<"] = "<<PF_pt[iPF]<<" with current sum: "<<sumEtPF<<endl;
		    //cout<<"PF eta["<<iPF<<"] = "<<PF_eta[iPF]<<", PF phi["<<iPF<<"] = "<<PF_phi[iPF]<<endl;
		    //cout<<"delta R = "<<dRhitJetPF<<endl;
		  }	      
		//For the JetShapes analysis get the energy in each ring in dR(\eta,\phi)
		for (int iRing = 0; iRing<=nRings; iRing++)
		  {
		    if (dRhitJetPF> LOW_dR[iRing] && dRhitJetPF<= HI_dR[iRing])
		      {
			ringSumEt[iRing] += PF_pt[iPF];
			ringCounter[iRing] +=1;
			//these are the charged tracks
			if ( PF_id[iPF] ==1 || PF_id[iPF]==2 || PF_id[iPF]==3 )
			  {
			    ringSumEtCh[iRing]+= PF_pt[iPF];
			    ringCounterCh[iRing]+=1;
			  }
		      } 
		  }	      
	      }
	  
	    //Now gather up the sums for each jet	
	    PFsumEt_ring0[iJet] = ringSumEt[0];
	    PFsumEt_ring1[iJet] = ringSumEt[1];
	    PFsumEt_ring2[iJet] = ringSumEt[2];
	    PFsumEt_ring3[iJet] = ringSumEt[3];
	    PFsumEt_ring4[iJet] = ringSumEt[4];
	    PFsumEt_ring5[iJet] = ringSumEt[5];
	    nPFobjs_ring0[iJet] = ringCounter[0];
	    nPFobjs_ring1[iJet] = ringCounter[1];
	    nPFobjs_ring2[iJet] = ringCounter[2];
	    nPFobjs_ring3[iJet] = ringCounter[3];
	    nPFobjs_ring4[iJet] = ringCounter[4];
	    nPFobjs_ring5[iJet] = ringCounter[5];
	    PFchSumEt_ring0[iJet] = ringSumEtCh[0];
	    PFchSumEt_ring1[iJet] = ringSumEtCh[1];
	    PFchSumEt_ring2[iJet] = ringSumEtCh[2];
	    PFchSumEt_ring3[iJet] = ringSumEtCh[3];
	    PFchSumEt_ring4[iJet] = ringSumEtCh[4];
	    PFchSumEt_ring5[iJet] = ringSumEtCh[5];
	    nPFchObjs_ring0[iJet] = ringCounterCh[0];
	    nPFchObjs_ring1[iJet] = ringCounterCh[1];
	    nPFchObjs_ring2[iJet] = ringCounterCh[2];
	    nPFchObjs_ring3[iJet] = ringCounterCh[3];
	    nPFchObjs_ring4[iJet] = ringCounterCh[4];
	    nPFchObjs_ring5[iJet] = ringCounterCh[5];
	    PFsumEt[iJet] = sumEtPF;
	    //cout<<"Event no = "<<iJet<<", sumEtPF = "<<sumEtPF<<endl;
	    nPFobjs[iJet] = PFCounter;

	  }else if(jet_type=="Calo"){

	    //	  //----------------------------------------------------------------------
	    // 	  //  Now loop over the recHits(Towers ), and get the ones inside the jet	
	    // 	  //----------------------------------------------------------------------
	    double sumEtTower=0;
	    Int_t TowerCounter = 0;	  
	    double ringTowerSumEt[nRings]={0,0,0,0,0,0};
	    Int_t ringTowerCounter[nRings] ={0,0,0,0,0,0};
	    for(int iTower = 0; iTower<Tower_n; iTower++)
	      {
		if (Tower_pt[iTower] < Tower_TRACKPTCUT) continue;
		if ( fabs(Tower_eta[iTower] ) >= Tower_TRACKETACUT ) continue;
		double dRhitJetTower = sqrt(((Tower_eta[iTower]-jteta[iJet])*(Tower_eta[iTower]-jteta[iJet])) + ((Tower_phi[iTower]-jtphi[iJet])*(Tower_phi[iTower]-jtphi[iJet])));	  
		if (dRhitJetTower<=JETRADIUS)
		  {
		    sumEtTower += Tower_pt[iTower];
		    TowerCounter +=1;
		    //cout<<"Toweret ["<<iTower<<"] = "<<Toweret[iTower]<<" with current sum: "<<sumEtTower<<endl;
		  }	      
		//For the JetShapes analysis get the energy in each ring in dR(\eta,\phi)
		for (int iRing = 0; iRing<=nRings; iRing++)
		  {
		    if (dRhitJetTower> LOW_dR[iRing] && dRhitJetTower<= HI_dR[iRing])
		      {
			ringTowerSumEt[iRing] += Tower_pt[iTower];
			ringTowerCounter[iRing] +=1;
		      } 
		  }	      
	      }
	  
	    //Now gather up the sums for each jet	
	    Tower_sumEt_ring0[iJet] = ringTowerSumEt[0];
	    Tower_sumEt_ring1[iJet] = ringTowerSumEt[1];
	    Tower_sumEt_ring2[iJet] = ringTowerSumEt[2];
	    Tower_sumEt_ring3[iJet] = ringTowerSumEt[3];
	    Tower_sumEt_ring4[iJet] = ringTowerSumEt[4];
	    Tower_sumEt_ring5[iJet] = ringTowerSumEt[5];
	    nTower_objs_ring0[iJet] = ringTowerCounter[0];
	    nTower_objs_ring1[iJet] = ringTowerCounter[1];
	    nTower_objs_ring2[iJet] = ringTowerCounter[2];
	    nTower_objs_ring3[iJet] = ringTowerCounter[3];
	    nTower_objs_ring4[iJet] = ringTowerCounter[4];
	    nTower_objs_ring5[iJet] = ringTowerCounter[5];

	    Tower_sumEt[iJet] = sumEtTower;
	    nTower_objs[iJet] = TowerCounter;

	  }
// 	  //----------------------------------------------------------------------
// 	  //  Now loop over the recHits(HBHE), and get the ones inside the jet	
// 	  //----------------------------------------------------------------------
// 	  double sumEtHBHE=0;
// 	  Int_t recHitHBHECounter = 0;
// 	  for(int iHBHE = 0; iHBHE<HBHEn; iHBHE++)
// 	    {
// 	      double dRhitJetHBHE = sqrt(((HBHEeta[iHBHE]-jteta[iJet])*(HBHEeta[iHBHE]-jteta[iJet])) + (((HBHEphi[iHBHE]-jtphi[iJet])*(HBHEphi[iHBHE]-jtphi[iJet]))));	  
// 	      if (dRhitJetHBHE<=JETRADIUS)
// 		{
// 		  sumEtHBHE += HBHEet[iHBHE];
// 		  recHitHBHECounter +=1;
// 		  //cout<<"HBHEet ["<<iHBHE<<"] = "<<HBHEet[iHBHE]<<" with current sum: "<<sumEtHBHE<<endl;
// 		}
// 	    }
// 	  HBHEsumEt[iJet] = sumEtHBHE;
// 	  nHBHErhits[iJet] = recHitHBHECounter;
// 	  //----------------------------------------------------------------------
// 	  //  Now loop over the recHits(EE), and get the ones inside the jet	
// 	  //----------------------------------------------------------------------
// 	  double sumEtEE=0;
// 	  Int_t recHitEECounter = 0;
// 	  for(int iEE = 0; iEE<EEn; iEE++)
// 	    {
// 	      double dRhitJetEE = sqrt(((EEeta[iEE]-jteta[iJet])*(EEeta[iEE]-jteta[iJet])) + (((EEphi[iEE]-jtphi[iJet])*(EEphi[iEE]-jtphi[iJet]))));	  
// 	      if (dRhitJetEE<=JETRADIUS)
// 		{
// 		  sumEtEE += EEet[iEE];
// 		  //cout<<"EEet ["<<iEE<<"] = "<<EEet[iEE]<<" with current sum: "<<sumEtEE<<endl;
// 		  recHitEECounter +=1;
// 		}
// 	    }
// 	  EEsumEt[iJet] = sumEtEE;
// 	  nEErhits[iJet] = recHitEECounter;
// 	  //----------------------------------------------------------------------
// 	  //  Now loop over the recHits(EB), and get the ones inside the jet	
// 	  //----------------------------------------------------------------------
// 	  double sumEtEB=0;
// 	  Int_t recHitEBCounter = 0;
// 	  for(int iEB = 0; iEB<EBn; iEB++)
// 	    {
// 	      double dRhitJetEB = sqrt(((EBeta[iEB]-jteta[iJet])*(EBeta[iEB]-jteta[iJet])) + (((EBphi[iEB]-jtphi[iJet])*(EBphi[iEB]-jtphi[iJet]))));	  
// 	      if (dRhitJetEB<=JETRADIUS)
// 		{
// 		  sumEtEB += EBet[iEB];
// 		  recHitEBCounter +=1;
// 		}
// 	    }
// 	  EBsumEt[iJet] = sumEtEB;
// 	  nEBrhits[iJet] = recHitEBCounter;
	  
	}

      
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //  Get random distribution in eta and phi
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------  

    
      for (int iRan =0 ; iRan<nRandom; iRan++ )
	{
	  float ranEta = myDice.Uniform(-2,2);
	  float ranPhi = myDice.Uniform(-3.14159,3.13159);

	  //cout<<"ranEta = "<<ranEta<<", ranPhi = "<<ranPhi<<endl;

	  /*

	  //----------------------------------------------------------------------
	  //  Now loop over the tracks, and get the ones inside random cones
	  //----------------------------------------------------------------------
	  int ranConeTrkCounter = 0;
	  double  ranConeTrkSumPt = 0;
	  double ranRingSumEtTk[nRings]={0,0,0,0,0,0};
	  Int_t ranRingCounterTk[nRings] ={0,0,0,0,0,0};
	  for (int jTrk = 0; jTrk<nTracks ; jTrk++)
	    {
	      if (trkPt[jTrk] < TRACKPTCUT) continue;
	      //float iiieta = trkEta[iTrk];
	      if (fabs(trkEta[jTrk]) > TRACKETACUT ) continue;
	      double dr1 = sqrt(((trkEta[jTrk]-ranEta)*(trkEta[jTrk]-ranEta)) + (((trkPhi[jTrk]-ranPhi)*(trkPhi[jTrk]-ranPhi))));
	      if (dr1 <= JETRADIUS )
		{
		  ranConeTrkCounter +=1;
		  ranConeTrkSumPt += trkPt[jTrk];
		}
	      ranNtracks[iRan] = ranConeTrkCounter;
	      ranTrkSumPt[iRan] = ranConeTrkSumPt;
	      ranConeEta[iRan] = ranEta;
	      ranConePhi[iRan] = ranPhi;

	      //For the JetShapes analysis get the energy in each ring in dR(\eta,\phi)
	      for (int iRing = 0; iRing<=nRings; iRing++)
		{
		  if (dr1> LOW_dR[iRing] && dr1<= HI_dR[iRing])
		    {
		      ranRingSumEtTk[iRing] += trkPt[jTrk];
		      ranRingCounterTk[iRing] +=1;

		    } 
		}	  
	    }
	  //Now gather up the sums for each jet
	  
	  ranTkSumEt_ring0[iRan] = ranRingSumEtTk[0];
	  ranTkSumEt_ring1[iRan] = ranRingSumEtTk[1];
	  ranTkSumEt_ring2[iRan] = ranRingSumEtTk[2];
	  ranTkSumEt_ring3[iRan] = ranRingSumEtTk[3];
	  ranTkSumEt_ring4[iRan] = ranRingSumEtTk[4];
	  ranTkSumEt_ring5[iRan] = ranRingSumEtTk[5];
	  
	  nRanTks_ring0[iRan] = ranRingCounterTk[0];
	  nRanTks_ring1[iRan] = ranRingCounterTk[1];
	  nRanTks_ring2[iRan] = ranRingCounterTk[2];
	  nRanTks_ring3[iRan] = ranRingCounterTk[3];
	  nRanTks_ring4[iRan] = ranRingCounterTk[4];
	  nRanTks_ring5[iRan] = ranRingCounterTk[5];
	  
	  */
	
	  //----------------------------------------------------------------------
	  //  Now loop over the PFlow object, and get the ones inside the random cone	
	  //----------------------------------------------------------------------

	  if(jet_type=="PF"){

	    double ranSumEtPF=0;
	    Int_t ranPFCounter = 0;
	    double ranRingSumEt[nRings] ={0,0,0,0,0,0};
	    Int_t ranRingCounter[nRings] ={0,0,0,0,0,0};
	    double ranRingSumEtCh[nRings] ={0,0,0,0,0,0};
	    Int_t ranRingCounterCh[nRings] ={0,0,0,0,0,0};
	    for(int iPF = 0; iPF<PF_n; iPF++)
	      {
		if (PF_pt[iPF] < PF_TRACKPTCUT) continue;
		if (fabs(PF_eta[iPF]) > PF_TRACKETACUT ) continue;
		double dR12 = sqrt(((PF_eta[iPF]-ranEta)*(PF_eta[iPF]-ranEta)) + (((PF_phi[iPF]-ranPhi)*(PF_phi[iPF]-ranPhi))));	  
		if (dR12<=JETRADIUS)
		  {
		    ranSumEtPF += PF_pt[iPF];
		    ranPFCounter +=1;
		    //cout<<"random cone PF_pt ["<<iPF<<"] = "<<PF_pt[iPF]<<" with current sum: "<<ranSumEtPF<<endl;
		  }
		//For the JetShapes analysis get the energy in each ring in dR(\eta,\phi)
		for (int iRing = 0; iRing<=nRings; iRing++)
		  {
		    if (dR12> LOW_dR[iRing] && dR12<= HI_dR[iRing])
		      {
			ranRingSumEt[iRing] += PF_pt[iPF];
			ranRingCounter[iRing] +=1;
			if (PF_id[iPF] ==1 || PF_id[iPF]==2 || PF_id[iPF]==3)//these are the charged tracks
			  {
			    ranRingSumEtCh[iRing]+= PF_pt[iPF];
			    ranRingCounterCh[iRing]+=1;
			  }
		      }
		  }
	      }
	    //Now gather up the sums for each jet

	    ranPFsumEt_ring0[iRan] = ranRingSumEt[0];
	    ranPFsumEt_ring1[iRan] = ranRingSumEt[1];
	    ranPFsumEt_ring2[iRan] = ranRingSumEt[2];
	    ranPFsumEt_ring3[iRan] = ranRingSumEt[3];
	    ranPFsumEt_ring4[iRan] = ranRingSumEt[4];
	    ranPFsumEt_ring5[iRan] = ranRingSumEt[5];

	    ranNPFobjs_ring0[iRan] = ranRingCounter[0];
	    ranNPFobjs_ring1[iRan] = ranRingCounter[1];
	    ranNPFobjs_ring2[iRan] = ranRingCounter[2];
	    ranNPFobjs_ring3[iRan] = ranRingCounter[3];
	    ranNPFobjs_ring4[iRan] = ranRingCounter[4];
	    ranNPFobjs_ring5[iRan] = ranRingCounter[5];

	    ranPFsumEt[iRan] = ranSumEtPF;
	    ranNPFobjs[iRan] = ranPFCounter;

	    //cout<<"Event no = "<<iev<<", sum of random PF candidates in the R=0."<<rad<<" = "<<ranSumEtPF<<endl;

	    ranPFchSumEt_ring0[iRan] = ranRingSumEtCh[0];
	    ranPFchSumEt_ring1[iRan] = ranRingSumEtCh[1];
	    ranPFchSumEt_ring2[iRan] = ranRingSumEtCh[2];
	    ranPFchSumEt_ring3[iRan] = ranRingSumEtCh[3];
	    ranPFchSumEt_ring4[iRan] = ranRingSumEtCh[4];
	    ranPFchSumEt_ring5[iRan] = ranRingSumEtCh[5];

	    ranNPFchObjs_ring0[iRan] = ranRingCounterCh[0];
	    ranNPFchObjs_ring1[iRan] = ranRingCounterCh[1];
	    ranNPFchObjs_ring2[iRan] = ranRingCounterCh[2];
	    ranNPFchObjs_ring3[iRan] = ranRingCounterCh[3];
	    ranNPFchObjs_ring4[iRan] = ranRingCounterCh[4];
	    ranNPFchObjs_ring5[iRan] = ranRingCounterCh[5];


	  }else if(jet_type=="CALO"){


	    double ranSumEtTower=0;
	    Int_t ranTowerCounter = 0;
	    double ranRingSumEt[nRings] ={0,0,0,0,0,0};
	    Int_t ranRingCounter[nRings] ={0,0,0,0,0,0};
	    double ranRingSumEtCh[nRings] ={0,0,0,0,0,0};
	    Int_t ranRingCounterCh[nRings] ={0,0,0,0,0,0};
	    for(int iTower = 0; iTower<Tower_n; iTower++)
	      {
		if (Tower_pt[iTower] < Tower_TRACKPTCUT) continue;
		if (fabs(Tower_eta[iTower]) > Tower_TRACKETACUT ) continue;
		double dR12 = sqrt(((Tower_eta[iTower]-ranEta)*(Tower_eta[iTower]-ranEta)) + (((Tower_phi[iTower]-ranPhi)*(Tower_phi[iTower]-ranPhi))));	  
		if (dR12<=JETRADIUS)
		  {
		    ranSumEtTower += Tower_pt[iTower];
		    ranTowerCounter +=1;
		    //cout<<"Toweret ["<<iTower<<"] = "<<Toweret[iTower]<<" with current sum: "<<sumEtTower<<endl;
		  }
		//For the JetShapes analysis get the energy in each ring in dR(\eta,\phi)
		for (int iRing = 0; iRing<=nRings; iRing++)
		  {
		    if (dR12> LOW_dR[iRing] && dR12<= HI_dR[iRing])
		      {
			ranRingSumEt[iRing] += Tower_pt[iTower];
			ranRingCounter[iRing] +=1;
		      }
		  }
	      }
	    //Now gather up the sums for each jet

	    ranTower_sumEt_ring0[iRan] = ranRingSumEt[0];
	    ranTower_sumEt_ring1[iRan] = ranRingSumEt[1];
	    ranTower_sumEt_ring2[iRan] = ranRingSumEt[2];
	    ranTower_sumEt_ring3[iRan] = ranRingSumEt[3];
	    ranTower_sumEt_ring4[iRan] = ranRingSumEt[4];
	    ranTower_sumEt_ring5[iRan] = ranRingSumEt[5];

	    ranNTower_objs_ring0[iRan] = ranRingCounter[0];
	    ranNTower_objs_ring1[iRan] = ranRingCounter[1];
	    ranNTower_objs_ring2[iRan] = ranRingCounter[2];
	    ranNTower_objs_ring3[iRan] = ranRingCounter[3];
	    ranNTower_objs_ring4[iRan] = ranRingCounter[4];
	    ranNTower_objs_ring5[iRan] = ranRingCounter[5];

	    ranTower_sumEt[iRan] = ranSumEtTower;
	    ranNTower_objs[iRan] = ranTowerCounter;



	  }
// 	  //----------------------------------------------------------------------
// 	  //  Now loop over the recHits(HBHE), and get the ones inside the random cone	
// 	  //----------------------------------------------------------------------
// 	  double ranSumEtHBHE=0;
// 	  Int_t ranRecHitHBHECounter = 0;
// 	  for(int kHBHE = 0; kHBHE<HBHEn; kHBHE++)
// 	    {
// 	      double dr2 = sqrt(((HBHEeta[kHBHE]-ranEta)*(HBHEeta[kHBHE]-ranEta)) + (((HBHEphi[kHBHE]-ranPhi)*(HBHEphi[kHBHE]-ranPhi))));	  
// 	      if (dr2<=JETRADIUS)
// 		{
// 		  ranSumEtHBHE += HBHEet[kHBHE];
// 		  ranRecHitHBHECounter +=1;
// 		}
// 	    }
// 	  ranHBHEsumEt[iRan] = ranSumEtHBHE;
// 	  ranNHBHErhits[iRan] = ranRecHitHBHECounter;
// 	  //----------------------------------------------------------------------
// 	  //  Now loop over the recHits(EB), and get the ones inside the random cone	
// 	  //----------------------------------------------------------------------
// 	  double ranSumEtEB=0;
// 	  Int_t ranRecHitEBCounter = 0;
// 	  for(int kEB = 0; kEB<EBn; kEB++)
// 	    {
// 	      double dr2 = sqrt(((EBeta[kEB]-ranEta)*(EBeta[kEB]-ranEta)) + (((EBphi[kEB]-ranPhi)*(EBphi[kEB]-ranPhi))));	  
// 	      if (dr2<=JETRADIUS)
// 		{
// 		  ranSumEtEB += EBet[kEB];
// 		  ranRecHitEBCounter +=1;
// 		}
// 	    }
// 	  ranEBsumEt[iRan] = ranSumEtEB;
// 	  ranNEBrhits[iRan] = ranRecHitEBCounter;
// 	  //----------------------------------------------------------------------
// 	  //  Now loop over the recHits(EE), and get the ones inside the random cone	
// 	  //----------------------------------------------------------------------
// 	  double ranSumEtEE=0;
// 	  Int_t ranRecHitEECounter = 0;
// 	  for(int kEE = 0; kEE<EEn; kEE++)
// 	    {
// 	      double dr2 = sqrt(((EEeta[kEE]-ranEta)*(EEeta[kEE]-ranEta)) + (((EEphi[kEE]-ranPhi)*(EEphi[kEE]-ranPhi))));	  
// 	      if (dr2<=JETRADIUS)
// 		{
// 		  ranSumEtEE += EEet[kEE];
// 		  ranRecHitEECounter +=1;
// 		}
// 	    }
// 	  ranEEsumEt[iRan] = ranSumEtEE;
// 	  ranNEErhits[iRan] = ranRecHitEECounter;
	  
	}//Random number loop
          
      nt->Fill(); 
      //cout<<"done with this event "<<iev<<endl;
    }  
  cout<<"out of the loop"<<endl;  
  outf->Write();
  outf->Close();

}
