// Raghav Kunnawalkam Elayavalli
// June 18th 2014
// CERN

// 
// Macro to look specifically at the fake or supernova jets in the Jet55or65 PbPb data forest.
// 
// This macro will look at specific jet quality cuts (or if possible track cuts) aimed at understanding where these crazy jets come from and help to remove them. 
// list of cuts: 
// 1) chargedMax / jtpt > 0.01
// 2) TMath::Max(charged Max, neutral Max) / TMath::Max(charged Sum, neutral Sum) < 0.975 - from 12003 and 14007
// 3) jtpt/trgObj_pt < 2. this should also remove all weird jets with much greater pt than the triggered object. 

// Modify the script to run on condor! for each file. otherwise it takes forever to run on all those files in TChain.
// - Done 

// better question:  should i include the duplicate event calculation here? since it requires the loop structure. 
// naa i'll add it in a separate include. 



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

using namespace std;
  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};



//static const double boundaries_run_job[job_no+1] = {};

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

void RAA_fakecheck(int startfile = 0, int endfile = 1, int radius = 3, char *algo = "Vs"){

  TH1::SetDefaultSumw2();
  TStopwatch timer;
  timer.Start();

  cout<<"Radius = "<<radius<<" and Algo = "<<algo<<endl;

  // Change to macro to run on condor since its taking a freaking long time. 
  // 

  std::string infile;
  infile = "jet55or65_filelist.txt";
  
  std::ifstream instr(infile.c_str(),std::ifstream::in);
  std::string filename;
  //int nFiles = 11;

  //just to read the files till the start number
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr>>filename;
  }

  const int N = 5;
  //Create chain
  TChain* ch[N];

  string dir[N] = {
    "hltanalysis",
    "skimanalysis",
    //"hcalNoise",
    Form("ak%s%dPFJetAnalyzer",algo,radius),
    //"akPu5PFJetAnalyzer",
    //"multiPhotonAnalyzer",
    //"ppTrack",
    //"pfcandAnalyzer",
    //"anaMET",
    //"muonTree",
    "hiEvtAnalyzer",
    "hltobject"
  };
    
  string trees[N] = {
    "HltTree",
    "HltTree",
    //"hbhenoise",
    "t",
    //"t",
    //"photon",
    //"trackTree",
    //"pfTree",
    //"metTree",
    //"HLTMuTree",
    "HiTree",
    "jetObjTree"
  };

  // data files for PbPb Jet55or66 are a list. So we need to TChain them first and then we should be able to use those chaing to get the events. they will have the same names as we used here. 

  //TChain *jetpbpb1 = new TChain(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  //TChain *evtpbpb1 = new TChain("hiEvtAnalyzer/HiTree");
  //TChain *hltpbpb1 = new TChain("hltanalysis/HltTree");
  //TChain *skmpbpb1 = new TChain("skimanalysis/HltTree");
  //TChain *hltobjpbpb1 = new TChain("hltobject/jetObjTree");
  
  //for(int i = 0;i<N;i++)       ch[i] = new TChain(string(dir[i]+"/"+trees[i]).data());


  for(int ifile = startfile;ifile<endfile;ifile++){
    
    instr>>filename;
    cout<<"File: "<<ifile<<" = "<<filename<<endl;

    for(int i = 0;i<N;i++){
      ch[i] = new TChain(string(dir[i]+"/"+trees[i]).data());
      ch[i]->Add(filename.c_str());
      cout << "Tree loaded  " << string(dir[i]+"/"+trees[i]).data() << endl;
      cout << "Entries : " << ch[i]->GetEntries() << endl;

    }

    //jetpbpb1->Add(filename.c_str());
    //hltpbpb1->Add(filename.c_str());
    //skmpbpb1->Add(filename.c_str());
    //hltobjpbpb1->Add(filename.c_str());
    //evtpbpb1->Add(filename.c_str());
    //cout<<"Entries = "<<jetpbpb1->GetEntries()<<endl;

  }


  //these ones were used when i ran the macro standalone in root. 
//   // file list - are all named like this: /mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet55or65/HIRun2011-14Mar2014-v2-6lumi-jet55or65-forest-v9/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_##.root
//   // where ## goes from 0-11

//   for(int i = 0;i<=11;i++){
//     if(i==9)continue;
//     jetpbpb1->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet55or65/HIRun2011-14Mar2014-v2-6lumi-jet55or65-forest-v9/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_%d.root",i));
//     hltpbpb1->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet55or65/HIRun2011-14Mar2014-v2-6lumi-jet55or65-forest-v9/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_%d.root",i));
//     skmpbpb1->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet55or65/HIRun2011-14Mar2014-v2-6lumi-jet55or65-forest-v9/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_%d.root",i));
//     evtpbpb1->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet55or65/HIRun2011-14Mar2014-v2-6lumi-jet55or65-forest-v9/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_%d.root",i));
//     hltobjpbpb1->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet55or65/HIRun2011-14Mar2014-v2-6lumi-jet55or65-forest-v9/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_%d.root",i));
  
//   }

  //jetpbpb1->AddFriend(evtpbpb1);
  //jetpbpb1->AddFriend(hltpbpb1);
  //jetpbpb1->AddFriend(skmpbpb1);
  //jetpbpb1->AddFriend(hltobjpbpb1);

  ch[2]->AddFriend(ch[0]);
  ch[2]->AddFriend(ch[1]);
  ch[2]->AddFriend(ch[3]);
  ch[2]->AddFriend(ch[4]);

  ch[3]->AddFriend(ch[0]);
  ch[3]->AddFriend(ch[1]);
  ch[3]->AddFriend(ch[2]);
  ch[3]->AddFriend(ch[4]);

  Float_t nEntries = (Float_t)ch[2]->GetEntries();
  cout<<"Total number of entries = "<<nEntries<<endl;

  //jetpbpb1->Print();
  //jetpbpb1->GetListOfBranches();
  //jetpbpb1->GetListOfFriends();

  //jetpbpb1->AddFriend(hltpbpb1);
  
  //cout<<"# of events which satisfy the Jet55 criteria = "<<ch[2]->GetEntries("HLT_HIJet55_v1")<<endl;
  //cout<<"# of events which satisfy the Jet65 criteria = "<<ch[2]->GetEntries("HLT_HIJet65_v1")<<endl;
  //cout<<"testing if HLT tree knows that branch, no of entries there = "<<hltpbpb1->GetEntries("HLT_HIJet55_v1")<<endl;
  //cout<<"# of events which satisfy Jet55 but fail Jet65 and Jet80 = "<<ch[2]->GetEntries("HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1")<<endl;
  

  

  // declare the output file here: 
  TFile fout(Form("pbpb_ak%d_%s_fakejet_histos_%d.root",radius,algo,endfile),"RECREATE");
  fout.cd();
  
  
  //histogram from the HLT_HIJet55 trigger alone
  TH1F *hJet55 = new TH1F("hJet55","Jets pt which are selected by HLT_HIJet55 trigger (along with eventsel)",1000,0,1000);
  TH1F *hJet55_QA1 = new TH1F("hJet55_QA1","HLT_HIJet55 jets with chMax/jtpt>0.01 (along with eventsel)",1000,0,1000);
  TH1F *hJet55_QA2 = new TH1F("hJet55_QA2","HLT_HIJet55 jets with Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975 (along with eventsel)",1000,0,1000);
  TH1F *hJet55_QA1_2 = new TH1F("hJet55_QA1_2","HLT_HIJet55 jets with chMax/jtpt>0.01 and Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975 (along with eventsel)",1000,0,1000);
  TH1F *hJet55_QA3 = new TH1F("hJet55_QA3","HLT_HIJet55 jets with jtpt/trgobjpt<3 (along with eventsel)",1000,0,1000);

  TH1F *hJet65 = new TH1F("hJet65","Jets pt which are selected by HLT_HIJet65 trigger (along with eventsel)",1000,0,1000);
  TH1F *hJet65_QA1 = new TH1F("hJet65_QA1","HLT_HIJet65 jets with chMax/jtpt>0.01 (along with eventsel)",1000,0,1000);
  TH1F *hJet65_QA2 = new TH1F("hJet65_QA2","HLT_HIJet65 jets with Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975 (along with eventsel)",1000,0,1000);
  TH1F *hJet65_QA1_2 = new TH1F("hJet65_QA1_2","HLT_HIJet65 jets with chMax/jtpt>0.01 and Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975 (along with eventsel)",1000,0,1000);
  TH1F *hJet65_QA3 = new TH1F("hJet65_QA3","HLT_HIJet65 jets with jtpt/trgobjpt<3 (along with eventsel)",1000,0,1000);

  //TH1F *hJet55_QA4 = new TH1F("hJet55_QA4","HLT_HIJet55 jets with chMax/jtpt>0.01 (along with eventsel)",1000,0,1000);
  TH3F *hJet55_3D = new TH3F("hJet55_3D","3D lego histogram of jets with Jet55 trigger, pt vs eta vs phi",60,-2.5,2.5,60,-3,3,1000,0,1000);
  TH1F *hJet55Fake = new TH1F("hJet55Fake","jets with pt>80 which pass Jet55 and fail higher triggers",1000,0,1000);
  // the above one is just to get a feel for all the jets in the Jet55 trigger. 

  TH1F *hJet55_trg = new TH1F("hJet55_trg","HLT_HIJet55 and trgpt>=55 and <65 along with eventsel",1000,0,1000);
  TH1F *hJet55_trg_QA1 = new TH1F("hJet55_trg_QA1","HLT_HI_Jet55 and trgpt>=55 and <65 eventsel and chMax/jtpt>0.01",1000,0,1000);
  TH1F *hJet55_trg_QA2 = new TH1F("hJet55_trg_QA2","HLT_HI_Jet55 and trgpt>=55 and <65 eventsel and  Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975",1000,0,1000);
  TH1F *hJet55_trg_QA1_2 = new TH1F("hJet55_trg_QA1_2","HLT_HI_Jet55 and trgpt>=55 and <65 eventsel and  Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975 and chMax/jtpt>0.01",1000,0,1000);
  TH1F *hJet55_trg_QA3 = new TH1F("hJet55_trg_QA3","HLT_HI_Jet55 and trgpt>=55 and <65 eventsel and jtpt/trgobjpt<3",1000,0,1000);
  TH1F *hJet55_only = new TH1F("hJet55_only","HLT_HIJet55 jets with no HLT_HIJet65 and no 80 along with eventsel",1000,0,1000);
 
  TH1F *hJet65_only = new TH1F("hJet65_only","HLT_HIJet65 jets with no HLT_HIJet80 along with eventsel",1000,0,1000);
  TH1F *hJet65_trg = new TH1F("hJet65_trg","HLT_HIJet65 and trgpt>=65 and <80 along with eventsel",1000,0,1000);
  TH1F *hJet65_trg_QA1 = new TH1F("hJet65_trg_QA1","HLT_HIJet65 and trgpt>=65 and <80 with chMax/jtpt>0.01 and with eventsel",1000,0,1000);
  TH1F *hJet65_trg_QA2 = new TH1F("hJet65_trg_QA2","HLT_HIJet65 and trgpt>=65 and <80 with evetsel and  Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975",1000,0,1000);
  TH1F *hJet65_trg_QA1_2 =  new TH1F("hJet65_trg_QA1_2","HLT_HI_Jet65 and trgpt>=65 and <80 eventsel and  Max(neutralMax,chargedMax)/Max(chargedSum,neutralSum)<0.975 and chMax/jtpt>0.01",1000,0,1000);
  TH1F *hJet65_trg_QA3 = new TH1F("hJet65_trg_QA3","HLT_HIJet 65 and trgpt>=65 and <80 with jtpt/trgobjpt<3 and eventsel",1000,0,1000);
  

  // should i make one with QA1 and QA3? since QA2 wont work in R=0.2 
  
  TH1F *hJet55_trg_QA1_3 = new TH1F("hJet55_trg_QA1_3","HLT_HIJet55 jets with trgpt>=55 and <65 with eventsel and QA1 and QA3",1000,0,1000);
  TH1F *hJet65_trg_QA1_3 = new TH1F("hJet65_trg_QA1_3","HLT_HIJet65 jets with trgpt>=65 and <80 with eventsel and QA1 and QA3",1000,0,1000);

  TH1F *hJet55_QA1_3 = new TH1F("hJet55_QA1_3","HLT_HIJet55 jets with eventsel and QA1 and QA3",1000,0,1000);
  TH1F *hJet65_QA1_3 = new TH1F("hJet65_QA1_3","HLT_HIJet65 jets with eventsel and QA1 and QA3",1000,0,1000);
  
  /*
  
  //include the jet80 tree just to find the numbers the fake cut removes. 
  TFile *fin80 = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HIRun2011-14Mar2014-v2-6lumi-jet80-forest-v4ANDv9-merged/0.root");
  TTree *Jet80 = (TTree*)fin80->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  //Jet80->Print();
  TTree *evt80 = (TTree*)fin80->Get("hiEvtAnalyzer/HiTree");
  //evt80->Print();
  TTree *hlt80 = (TTree*)fin80->Get("hltanalysis/HltTree");
  //hlt80->Print();
  TTree *skm80 = (TTree*)fin80->Get("skimanalysis/HltTree");
  //skm80->Print();
  TTree *Trg80 = (TTree*)fin80->Get("hltobject/jetObjTree");
  //Trg80->Print();
  
  Jet80->AddFriend(evt80);
  Jet80->AddFriend(hlt80);
  Jet80->AddFriend(skm80);
  Jet80->AddFriend(Trg80);
  */

  TCut jet55 = "HLT_HIJet55_v1";
  TCut jet65 = "HLT_HIJet65_v1";
  TCut jet80 = "HLT_HIJet80_v1";
  TCut evtSel = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2";
  TCut qalCut1 = "chargedMax/jtpt>0.01";
  TCut qalCut2 = "TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975";
  TCut qalCut3 = "jtpt/pt<3";
  //TCut qalCut4 = "";
  TCut jet55only = "HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1";
  TCut jet65only = "HLT_HIJet65_v1&&!HLT_HIJet80_v1";
  TCut largejetpt = "jtpt>80";
  TCut trg55 = "pt>=55&&pt<65";
  TCut trg65 = "pt>=65&&pt<80";
  TCut trg80 = "pt>=80";

  //TH1F total,evtSel,jet55,jet55_evtsel,jet55only_ectsel,jet55_evtSel_QA1,jet55_evtSel_QA2,jet55_evtSel_QA3,jet55_evtSel_QA1_2,jet55_evtSel_QA1_3,jet55_evtSel_QA3_2;
  //TH1F jet55_trg55,jet55_trg_evtSel,jet55only_trg55_evtSel,jet55_trg55_evtSel_QA1,jet55_trg55_evtSel_QA2,jet55_trg55_evtSel_QA3,jet55_trg55_evtSel_qalCut1_qalCut2,jet55_trg55_evtSel_qalCut1_qalCut3,jet55_trg55_evtSel_qalCut3_qalCut2;
  //TH1F jet65,jet65_evtSel,jet65only_evtSel,jet65_evtSel_qalCut1,jet65_evtSel_qalCut3,jet65_evtSel_qalCut3,jet65_evtSel_qalCut1_qalCut3,jet65_evtSel_qalCut1_qalCut2,jet65_evtSel_qalCut3_qalCut2;
  //TH1F jet65_

  
  //cout<<"total # of entries = "<<ch[2]->GetEntries()<<endl;
  /*
  cout<<"total # of entries = "<<ch[2]->GetEntries()<<endl;
  ch[3]->Draw("hiBin>>total");
  cout<<"with evtSel added  = "<<ch[2]->GetEntries(evtSel)<<endl;
  //ch[2]->Draw("evtSel","hiBin",);
  cout<<"with evtSel added  = "<<ch[2]->GetEntries(evtSel)<<endl;
  ch[3]->Draw("hiBin>>evtSel",evtSel);
  cout<<"with jet55"<<" = "<<ch[2]->GetEntries(jet55)<<endl;
  ch[3]->Draw("hiBin>>jet55",jet55);
  cout<<"with evtSel and Jet55 added = "<<ch[2]->GetEntries(jet55&&evtSel)<<endl;
  ch[3]->Draw("hiBin>>jet55_evtsel",jet55&&evtSel);
  cout<<"with Jet55_only and evtSel added = "<<ch[2]->GetEntries(jet55only&&evtSel)<<endl;
  ch[3]->Draw("hiBin>>jet55only_evtsel",jet55only&&evtSel);
  */
  //ch[2]->Project("total","hiBin");
  
  cout<<"with evtSel added  = "<<ch[2]->GetEntries(evtSel)<<endl;
  //ch[2]->Project("evtSel","hiBin",);
  cout<<"with evtSel added  = "<<ch[2]->GetEntries(evtSel)<<endl;
  //ch[2]->Project("evtSel","hiBin",evtSel);
  cout<<"with jet55"<<" = "<<ch[2]->GetEntries(jet55)<<endl;
  //ch[2]->Project("jet55","hiBin",jet55);
  cout<<"with evtSel and Jet55 added = "<<ch[2]->GetEntries(jet55&&evtSel)<<endl;
  //ch[2]->Project("jet55_evtsel","hiBin",jet55&&evtSel);
  cout<<"with Jet55_only and evtSel added = "<<ch[2]->GetEntries(jet55only&&evtSel)<<endl;
  //ch[2]->Project("jet55only_evtsel","hiBin",jet55only&&evtSel);

  cout<<"jet55 && evtSel && qalCut1"<<" = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut1)<<endl;
  //ch[2]->Project("jet55_evtSel_QA1","hiBin",jet55&&evtSel&&qalCut1);
  cout<<"jet55 && evtSel && qalCut2"<<" = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut2)<<endl;
  //ch[2]->Project("jet55_evtSel_QA2","hiBin",jet55&&evtSel&&qalCut2);
  cout<<"jet55 && evtSel && qalCut3"<<" = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut3)<<endl;
  //ch[2]->Project("jet55_evtSel_QA3","hiBin",jet55&&evtSel&&qalCut3);

  cout<<"jet55 and evtSel and QA1_2 = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut1&&qalCut2)<<endl;
  //ch[2]->Project("jet55_evtSel_QA1_QA2","hiBin",jet55&&evtSel&&qalCut1&&qalCut2);
  cout<<"jet55 and evtSel and QA1_3 = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut1&&qalCut3)<<endl;
  //ch[2]->Project("jet55_evtSel_QA1_QA3","hiBin",jet55&&evtSel&&qalCut1&&qalCut3);
  cout<<"jet55 and evtSel and QA2_3 = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut2&&qalCut3)<<endl;
  //ch[2]->Project("jet55_evtSel_QA3_QA2","hiBin",jet55&&evtSel&&qalCut3&&qalCut2);


  cout<<"with jet55 and trgObjpt selection"<<" = "<<ch[2]->GetEntries(jet55&&trg55)<<endl;
  //ch[2]->Project("jet55_trg55","hiBin",jet55&&trg55);
  cout<<"with evtSel and Jet55 and trgobjpt added = "<<ch[2]->GetEntries(jet55&&trg55&&evtSel)<<endl;
  //ch[2]->Project("jet55_trg55_evtSel","hiBin",jet55&&trg55&&evtSel);
  cout<<"with Jet55_only and trgObjpt and evtSel added = "<<ch[2]->GetEntries(jet55only&&trg55&&evtSel)<<endl;
  //ch[2]->Project("jet55only_trg55_evtSel","hiBin",jet55only&&trg55&&evtSel);

  cout<<"jet55 and trgObjpt && evtSel && qalCut1"<<" = "<<ch[2]->GetEntries(jet55&&trg55&&evtSel&&qalCut1)<<endl;
  //ch[2]->Project("jet55_trg55_evtSel_QA1","hiBin",jet55&&trg55&&evtSel&&qalCut1);
  cout<<"jet55 and trgObjpt && evtSel && qalCut2"<<" = "<<ch[2]->GetEntries(jet55&&trg55&&evtSel&&qalCut2)<<endl;
  //ch[2]->Project("jet55_trg55_evtSel_QA2","hiBin",jet55&&trg55&&evtSel&&qalCut2);
  cout<<"jet55 and trgObjpt && evtSel && qalCut3"<<" = "<<ch[2]->GetEntries(jet55&&trg55&&evtSel&&qalCut3)<<endl;
  //ch[2]->Project("jet55_trg55_evtSel_QA3","hiBin",jet55&&trg55&&evtSel&&qalCut3);

  cout<<"jet55 and trgObjpt and evtSel and QA1_2 = "<<ch[2]->GetEntries(jet55&&trg55&&evtSel&&qalCut1&&qalCut2)<<endl;
  //ch[2]->Project("jet55_trg55_evtSel_qalCut1_qalCut2","hiBin",jet55&&trg55&&evtSel&&qalCut1&&qalCut2);
  cout<<"jet55 and trgObjpt and evtSel and QA1_3 = "<<ch[2]->GetEntries(jet55&&trg55&&evtSel&&qalCut1&&qalCut3)<<endl;
  //ch[2]->Project("jet55_trg55_evtSel_qalCut1_qalCut3","hiBin",jet55&&trg55&&evtSel&&qalCut1&&qalCut3);
  cout<<"jet55 and trgObjpt and evtSel and QA2_3 = "<<ch[2]->GetEntries(jet55&&trg55&&evtSel&&qalCut2&&qalCut3)<<endl;
  //ch[2]->Project("jet55_trg55_evtSel_qalCut1_qa3Cut2","hiBin",jet55&&trg55&&evtSel&&qalCut3&&qalCut2);


  cout<<"with jet65"<<" = "<<ch[2]->GetEntries(jet65)<<endl;
  //ch[2]->Project("jet65","hiBin",jet65);
  cout<<"with evtSel and Jet65 added = "<<ch[2]->GetEntries(jet65&&evtSel)<<endl;
  //ch[2]->Project("jet65_evtSel","hiBin",jet65&&evtSel);
  cout<<"with Jet65_only and evtSel added = "<<ch[2]->GetEntries(jet65only&&evtSel)<<endl;
  //ch[2]->Project("jet65only_evtSel","hiBin",jet65only&&evtSel);

  cout<<"jet65 && evtSel && qalCut1"<<" = "<<ch[2]->GetEntries(jet65&&evtSel&&qalCut1)<<endl;
  //ch[2]->Project("jet65_evtSel_qalCut1","hiBin",jet65&&evtSel&&qalCut1);
  cout<<"jet65 && evtSel && qalCut2"<<" = "<<ch[2]->GetEntries(jet65&&evtSel&&qalCut2)<<endl;
  //ch[2]->Project("jet65_evtSel_qalCut2","hiBin",jet65&&evtSel&&qalCut2);
  cout<<"jet65 && evtSel && qalCut3"<<" = "<<ch[2]->GetEntries(jet65&&evtSel&&qalCut3)<<endl;
  //ch[2]->Project("jet65_evtSel_qalCut3","hiBin",jet65&&evtSel&&qalCut3);

  cout<<"jet65 and evtSel and QA1_2 = "<<ch[2]->GetEntries(jet65&&evtSel&&qalCut1&&qalCut2)<<endl;
  //ch[2]->Project("jet65_evtSel_qalCut1_qalCut2","hiBin",jet65&&evtSel&&qalCut1&&qalCut2);
  cout<<"jet65 and evtSel and QA1_3 = "<<ch[2]->GetEntries(jet65&&evtSel&&qalCut1&&qalCut3)<<endl;
  //ch[2]->Project("jet65_evtSel_qalCut1_qalCut3","hiBin",jet65&&evtSel&&qalCut1&&qalCut3);
  cout<<"jet65 and evtSel and QA2_3 = "<<ch[2]->GetEntries(jet65&&evtSel&&qalCut2&&qalCut3)<<endl;
  //ch[2]->Project("jet65_evtSel_qalCut3_qalCut2","hiBin",jet65&&evtSel&&qalCut3&&qalCut2);

  cout<<"with jet65 and trgObjpt selection"<<" = "<<ch[2]->GetEntries(jet65&&trg65)<<endl;
  //ch[2]->Project("jet65_trg65","hiBin",jet65&&trg65);
  cout<<"with evtSel and Jet65 and trgobjpt added = "<<ch[2]->GetEntries(jet65&&trg65&&evtSel)<<endl;
  //ch[2]->Project("jet65_trg65_evtSel","hiBin",jet65&&trg65&&evtSel);
  cout<<"with Jet65_only and trgObjpt and evtSel added = "<<ch[2]->GetEntries(jet65only&&trg65&&evtSel)<<endl;
  //ch[2]->Project("jet65only_trg65_evtSel","hiBin",jet65only&&trg65&&evtSel);

  cout<<"jet65 and trgObjpt && evtSel && qalCut1"<<" = "<<ch[2]->GetEntries(jet65&&trg65&&evtSel&&qalCut1)<<endl;
  //ch[2]->Project("jet65_trg65_evtSel_qalCut1","hiBin",jet65&&trg65&&evtSel&&qalCut1);
  cout<<"jet65 and trgObjpt && evtSel && qalCut2"<<" = "<<ch[2]->GetEntries(jet65&&trg65&&evtSel&&qalCut2)<<endl;
  //ch[2]->Project("jet65_trg65_evtSel_qalCut2","hiBin",jet65&&trg65&&evtSel&&qalCut2);
  cout<<"jet65 and trgObjpt && evtSel && qalCut3"<<" = "<<ch[2]->GetEntries(jet65&&trg65&&evtSel&&qalCut3)<<endl;
  //ch[2]->Project("jet65_trg65_evtSel_qalCut3","hiBin",jet65&&trg65&&evtSel&&qalCut3);

  cout<<"jet65 and trgObjpt and evtSel and QA1_2 = "<<ch[2]->GetEntries(jet65&&trg65&&evtSel&&qalCut1&&qalCut2)<<endl;
  //ch[2]->Project("jet65_trg65_evtSel_qalCut1_qalCut2","hiBin",jet65&&trg65&&evtSel&&qalCut1&&qalCut2);
  cout<<"jet65 and trgObjpt and evtSel and QA1_3 = "<<ch[2]->GetEntries(jet65&&trg65&&evtSel&&qalCut1&&qalCut3)<<endl;
  //ch[2]->Project("jet65_trg65_evtSel_qalCut1_qalCut3","hiBin",jet65&&trg65&&evtSel&&qalCut1&&qalCut3);
  cout<<"jet65 and trgObjpt and evtSel and QA2_3 = "<<ch[2]->GetEntries(jet65&&trg65&&evtSel&&qalCut2&&qalCut3)<<endl;
  //ch[2]->Project("jet65_trg65_evtSel_qalCut3_qalCut2","hiBin",jet65&&trg65&&evtSel&&qalCut3&&qalCut2);
  

  /*
  cout<<"with jet80"<<" = "<<Jet80->GetEntries(jet80)<<endl;
  cout<<"with evtSel and Jet80 added = "<<Jet80->GetEntries(jet80&&evtSel)<<endl;

  cout<<"jet80 && evtSel && qalCut1"<<" = "<<Jet80->GetEntries(jet80&&evtSel&&qalCut1)<<endl;
  cout<<"jet80 && evtSel && qalCut2"<<" = "<<Jet80->GetEntries(jet80&&evtSel&&qalCut2)<<endl;
  cout<<"jet80 && evtSel && qalCut3"<<" = "<<Jet80->GetEntries(jet80&&evtSel&&qalCut3)<<endl;

  cout<<"jet80 and evtSel and QA1_2 = "<<Jet80->GetEntries(jet80&&evtSel&&qalCut1&&qalCut2)<<endl;
  cout<<"jet80 and evtSel and QA1_3 = "<<Jet80->GetEntries(jet80&&evtSel&&qalCut1&&qalCut3)<<endl;
  cout<<"jet80 and evtSel and QA2_3 = "<<Jet80->GetEntries(jet80&&evtSel&&qalCut2&&qalCut3)<<endl;
  cout<<"with jet80 and trgObjpt selection"<<" = "<<Jet80->GetEntries(jet80&&trg80)<<endl;
  cout<<"with evtSel and Jet680 and trgobjpt added = "<<Jet80->GetEntries(jet80&&trg80&&evtSel)<<endl;
  //cout<<"with Jet65_only and trgObjpt and evtSel added = "<<ch[2]->GetEntries(jet65only&&trg65&&evtSel)<<endl;

  cout<<"jet80 and trgObjpt && evtSel && qalCut1"<<" = "<<Jet80->GetEntries(jet80&&trg80&&evtSel&&qalCut1)<<endl;
  cout<<"jet80 and trgObjpt && evtSel && qalCut2"<<" = "<<Jet80->GetEntries(jet80&&trg80&&evtSel&&qalCut2)<<endl;
  cout<<"jet80 and trgObjpt && evtSel && qalCut3"<<" = "<<Jet80->GetEntries(jet80&&trg80&&evtSel&&qalCut3)<<endl;

  cout<<"jet80 and trgObjpt and evtSel and QA1_2 = "<<Jet80->GetEntries(jet80&&trg80&&evtSel&&qalCut1&&qalCut2)<<endl;
  cout<<"jet80 and trgObjpt and evtSel and QA1_3 = "<<Jet80->GetEntries(jet80&&trg80&&evtSel&&qalCut1&&qalCut3)<<endl;
  cout<<"jet80 and trgObjpt and evtSel and QA2_3 = "<<Jet80->GetEntries(jet80&&trg80&&evtSel&&qalCut2&&qalCut3)<<endl;
  */



  cout<<"jet55 only and trig 55 and eventSel and largejetpt"<<" = "<<ch[2]->GetEntries(evtSel&&jet55only&&trg55&&largejetpt)<<endl;

  

  //ch[2]->Print();
  
  ch[2]->Project("hJet55_trg_QA1_3","jtpt",jet55&&trg55&&evtSel&&qalCut1&&qalCut3);
  hJet55_trg_QA1_3->Print("base");
  ch[2]->Project("hJet65_trg_QA1_3","jtpt",jet65&&trg65&&evtSel&&qalCut1&&qalCut3);
  hJet65_trg_QA1_3->Print("base");
  ch[2]->Project("hJet55_QA1_3","jtpt",jet55&&evtSel&&qalCut1&&qalCut3);
  hJet55_QA1_3->Print("base");
  ch[2]->Project("hJet65_QA1_3","jtpt",jet65&&evtSel&&qalCut1&&qalCut3);
  hJet65_QA1_3->Print("base");

  ch[2]->Project("hJet55","jtpt",jet55&&evtSel);
  hJet55->Print("base");
  ch[2]->Project("hJet55_QA1","jtpt",evtSel&&jet55&&qalCut1);
  hJet55_QA1->Print("base");
  ch[2]->Project("hJet55_QA2","jtpt",evtSel&&jet55&&qalCut2);
  hJet55_QA2->Print("base");
  ch[2]->Project("hJet55_QA1_2","jtpt",evtSel&&jet55&&qalCut1&&qalCut2);
  hJet55_QA1_2->Print("base");
  ch[2]->Project("hJet55_QA3","jtpt",evtSel&&jet55&&qalCut3);
  hJet55_QA3->Print("base");
  ch[2]->Project("hJet55Fake","jtpt",evtSel&&jet55only&&largejetpt);
  hJet55Fake->Print("base");
  ch[2]->Project("hJet55_only","jtpt",jet55only&&evtSel);
  hJet55_only->Print("base");

  ch[2]->Project("hJet55_3D","jteta:jtphi:jtpt",evtSel&&jet55only&&largejetpt);
  hJet55_3D->Print("base");

  ch[2]->Project("hJet55_trg","jtpt",evtSel&&jet55&&trg55);
  hJet55_trg->Print("base");
  ch[2]->Project("hJet55_trg_QA1","jtpt",evtSel&&jet55&&trg55&&qalCut1);
  hJet55_trg_QA1->Print("base");
  ch[2]->Project("hJet55_trg_QA2","jtpt",evtSel&&jet55&&trg55&&qalCut2);
  hJet55_trg_QA2->Print("base");
  ch[2]->Project("hJet55_trg_QA1_2","jtpt",evtSel&&jet55&&trg55&&qalCut1&&qalCut2);
  hJet55_trg_QA1_2->Print("base");
  ch[2]->Project("hJet55_trg_QA3","jtpt",evtSel&&jet55&&trg55&&qalCut3);
  hJet55_trg_QA3->Print("base");

  ch[2]->Project("hJet65","jtpt",jet65&&evtSel);
  hJet65->Print("base");
  ch[2]->Project("hJet65_QA1","jtpt",evtSel&&jet65&&qalCut1);
  hJet65_QA1->Print("base");
  ch[2]->Project("hJet65_QA2","jtpt",evtSel&&jet65&&qalCut2);
  hJet65_QA2->Print("base");
  ch[2]->Project("hJet65_QA1_2","jtpt",evtSel&&jet65&&qalCut1&&qalCut2);
  hJet65_QA1_2->Print("base");
  ch[2]->Project("hJet65_QA3","jtpt",evtSel&&jet65&&qalCut3);
  hJet65_QA3->Print("base");

  ch[2]->Project("hJet65_only","jtpt",evtSel&&jet65only);
  hJet65_only->Print("base");
  ch[2]->Project("hJet65_trg","jtpt",evtSel&&jet65&&trg65);
  hJet65_trg->Print("base");
  ch[2]->Project("hJet65_trg_QA1","jtpt",evtSel&&jet65&&trg65&&qalCut1);
  hJet65_trg_QA1->Print("base");
  ch[2]->Project("hJet65_trg_QA2","jtpt",evtSel&&jet65&&trg65&&qalCut2);
  hJet65_trg_QA2->Print("base");
  ch[2]->Project("hJet65_trg_QA3","jtpt",evtSel&&jet65&&trg65&&qalCut3);
  hJet65_trg_QA3->Print("base");
  ch[2]->Project("hJet65_trg_QA1_2","jtpt",evtSel&&jet65&&trg65&&qalCut1&&qalCut2);
  hJet65_trg_QA1_2->Print("base");
  

  hJet55 = (TH1F*)hJet55->Rebin(nbins_pt,"hJet55",boundaries_pt);
  divideBinWidth(hJet55);
  hJet55_QA1 = (TH1F*)hJet55_QA1->Rebin(nbins_pt,"hJet55_QA1",boundaries_pt);
  divideBinWidth(hJet55_QA1);
  hJet55_QA2 = (TH1F*)hJet55_QA2->Rebin(nbins_pt,"hJet55_QA2",boundaries_pt);
  divideBinWidth(hJet55_QA2);
  hJet55_QA1_2 = (TH1F*)hJet55_QA1_2->Rebin(nbins_pt,"hJet55_QA1_2",boundaries_pt);
  divideBinWidth(hJet55_QA1_2);
  hJet55_QA3 = (TH1F*)hJet55_QA3->Rebin(nbins_pt,"hJet55_QA3",boundaries_pt);
  divideBinWidth(hJet55_QA3);
  hJet55_only = (TH1F*)hJet55_only->Rebin(nbins_pt,"hJet55_only",boundaries_pt);
  divideBinWidth(hJet55_only);
  hJet55Fake = (TH1F*)hJet55Fake->Rebin(nbins_pt,"hJet55Fake",boundaries_pt);
  divideBinWidth(hJet55Fake);
  hJet55_trg = (TH1F*)hJet55_trg->Rebin(nbins_pt,"hJet55_trg",boundaries_pt);
  divideBinWidth(hJet55_trg);
  hJet55_trg_QA1 = (TH1F*)hJet55_trg_QA1->Rebin(nbins_pt,"hJet55_trg_QA1",boundaries_pt);
  divideBinWidth(hJet55_trg_QA1);
  hJet55_trg_QA2 = (TH1F*)hJet55_trg_QA2->Rebin(nbins_pt,"hJet55_trg_QA2",boundaries_pt);
  divideBinWidth(hJet55_trg_QA2);
  hJet55_trg_QA1_2 = (TH1F*)hJet55_trg_QA1_2->Rebin(nbins_pt,"hJet55_trg_QA1_2",boundaries_pt);
  divideBinWidth(hJet55_trg_QA1_2);
  hJet55_trg_QA3 = (TH1F*)hJet55_trg_QA3->Rebin(nbins_pt,"hJet55_trg_QA3",boundaries_pt);
  divideBinWidth(hJet55_trg_QA3);

  hJet65 = (TH1F*)hJet65->Rebin(nbins_pt,"hJet65",boundaries_pt);
  divideBinWidth(hJet65);
  hJet65_QA1 = (TH1F*)hJet65_QA1->Rebin(nbins_pt,"hJet65_QA1",boundaries_pt);
  divideBinWidth(hJet65_QA1);
  hJet65_QA2 = (TH1F*)hJet65_QA2->Rebin(nbins_pt,"hJet65_QA2",boundaries_pt);
  divideBinWidth(hJet65_QA2);
  hJet65_QA1_2 = (TH1F*)hJet65_QA1_2->Rebin(nbins_pt,"hJet65_QA1_2",boundaries_pt);
  divideBinWidth(hJet65_QA1_2);
  hJet65_QA3 = (TH1F*)hJet65_QA3->Rebin(nbins_pt,"hJet65_QA3",boundaries_pt);
  divideBinWidth(hJet65_QA3);
  hJet65_only = (TH1F*)hJet65_only->Rebin(nbins_pt,"hJet65_only",boundaries_pt);
  divideBinWidth(hJet65_only);
  //hJet65_fake = (TH1F*)hJet65_fake->Rebin(nbins_pt,"hJet65_fake",boundaries_pt);
  //divideBinWidth(hJet65_fake);
  hJet65_trg = (TH1F*)hJet65_trg->Rebin(nbins_pt,"hJet65_trg",boundaries_pt);
  divideBinWidth(hJet65_trg);
  hJet65_trg_QA1 = (TH1F*)hJet65_trg_QA1->Rebin(nbins_pt,"hJet65_trg_QA1",boundaries_pt);
  divideBinWidth(hJet65_trg_QA1);
  hJet65_trg_QA2 = (TH1F*)hJet65_trg_QA2->Rebin(nbins_pt,"hJet65_trg_QA2",boundaries_pt);
  divideBinWidth(hJet65_trg_QA2);
  hJet65_trg_QA1_2 = (TH1F*)hJet65_trg_QA1_2->Rebin(nbins_pt,"hJet65_trg_QA1_2",boundaries_pt);
  divideBinWidth(hJet65_trg_QA1_2);
  hJet65_trg_QA3 = (TH1F*)hJet65_trg_QA3->Rebin(nbins_pt,"hJet65_trg_QA3",boundaries_pt);
  divideBinWidth(hJet65_trg_QA3);

  hJet55_trg_QA1_3 = (TH1F*)hJet55_trg_QA1_3->Rebin(nbins_pt,"hJet55_trg_QA1_3",boundaries_pt);
  divideBinWidth(hJet55_trg_QA1_3);
  hJet65_trg_QA1_3 = (TH1F*)hJet65_trg_QA1_3->Rebin(nbins_pt,"hJet65_trg_QA1_3",boundaries_pt);
  divideBinWidth(hJet65_trg_QA1_3);
  hJet55_QA1_3 = (TH1F*)hJet55_QA1_3->Rebin(nbins_pt,"hJet55_QA1_3",boundaries_pt);
  divideBinWidth(hJet55_QA1_3);
  hJet65_QA1_3 = (TH1F*)hJet65_QA1_3->Rebin(nbins_pt,"hJet65_QA1_3",boundaries_pt);
  divideBinWidth(hJet65_QA1_3);
  

  fout.Write();
  fout.Close();
  
  timer.Stop();
  cout<<"Real time(sec) = "<<timer.RealTime()<<endl;
  cout<<"CPU time(sec)  = "<<timer.CpuTime()<<endl;
  cout<<"All done"<<endl;

}
