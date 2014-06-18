// Raghav Kunnawalkam Elayavalli
// June 18th 2014
// CERN

// 
// Macro to look specifically at the fake or supernova jets in the Jet55or65 PbPb data forest.
// 
// This macro will look at specific jet quality cuts (or if possible track cuts) aimed at understanding where these crazy jets come from and help to remove them. 
// list of cuts: 
// 1) chargedMax / jtpt > 0.01
// 2) neutral Max / Max(charged Sum, neutral Sum) >= 0.975 - from 12003 
// 3) 

// Modify the script to run on condor! for each file. otherwise it takes forever to run on all those files in TChain. 



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
  
  for(int ifile = startfile;ifile<endfile;ifile++){
    
    instr>>filename;
    cout<<"File: "<<filename<<endl;

    for(int i = 0;i<N;i++){

      ch[i] = new TChain(string(dir[i]+"/"+trees[i]).data()) ;
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


  Float_t nEntries = (Float_t)ch[2]->GetEntries();
  cout<<"Total number of entries = "<<nEntries<<endl;

  //jetpbpb1->Print();
  //jetpbpb1->GetListOfBranches();
  //jetpbpb1->GetListOfFriends();

  //jetpbpb1->AddFriend(hltpbpb1);
  
  cout<<"# of events which satisfy the Jet55 criteria = "<<ch[2]->GetEntries("HLT_HIJet55_v1")<<endl;
  cout<<"# of events which satisfy the Jet65 criteria = "<<ch[2]->GetEntries("HLT_HIJet65_v1")<<endl;
  //cout<<"testing if HLT tree knows that branch, no of entries there = "<<hltpbpb1->GetEntries("HLT_HIJet55_v1")<<endl;
  cout<<"# of events which satisfy Jet55 but fail Jet65 and Jet80 = "<<ch[2]->GetEntries("HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1")<<endl;


  // declare the output file here: 
  TFile fout(Form("pbpb_ak%d_%s_fakejet_histos_%d.root",radius,algo,endfile),"RECREATE");
  fout.cd();


  //histogram from the HLT_HIJet55 trigger alone
  TH1F *hJet55 = new TH1F("hJet55","Jets pt which are selected by HLT_HIJet55 trigger",1000,0,1000);
  TH1F *hJet55_QA1 = new TH1F("hJet55_QA1","HLT_HIJet55 jets with chMax/jtpt>0.01 (along with eventsel)",1000,0,1000);
  TH1F *hJet55_QA2 = new TH1F("hJet55_QA2","HLT_HIJet55 jets with neutralMax/Max(chargedSum,neutralSum) (along with eventsel)",1000,0,1000);
  //TH1F *hJet55_QA3 = new TH1F("hJet55_QA3","HLT_HIJet55 jets with  (along with eventsel)",1000,0,1000);
  //TH1F *hJet55_QA4 = new TH1F("hJet55_QA4","HLT_HIJet55 jets with chMax/jtpt>0.01 (along with eventsel)",1000,0,1000);
  TH3F *hJet55_3D = new TH3F("hJet55_3D","3D lego histogram of jets with Jet55 trigger, pt vs eta vs phi",60,-2.5,2.5,60,-3,3,1000,0,1000);
  TH1F *hJet55Fake = new TH1F("hJet55Fake","jets with pt>80 which pass Jet55 and fail higher triggers",1000,0,1000);
  // the above one is just to get a feel for all the jets in the Jet55 trigger. 

  TCut jet55 = "HLT_HIJet55_v1";
  TCut jet65 = "HLT_HIJet65_v1";
  TCut jet80 = "HLT_HIJet80_v1";
  TCut evtSel = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2";
  TCut qalCut1 = "chargedMax/jtpt>0.01";
  TCut qalCut2 = "neutralMax/TMath::Max(chargedSum,neutralSum)>=0.975";
  //TCut qalCut3 = "";
  //TCut qalCut4 = "";
  TCut fakeJet = "HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1&&jtpt>80";

  cout<<"jet55"<<" = "<<ch[2]->GetEntries(jet55)<<endl;
  cout<<"jet65"<<" = "<<ch[2]->GetEntries(jet65)<<endl;
  cout<<"jet80"<<" = "<<ch[2]->GetEntries(jet80)<<endl;
  cout<<"jet55 && evtSel"<<" = "<<ch[2]->GetEntries(jet55&&evtSel)<<endl;
  cout<<"jet55 && evtSel && qalCut1"<<" = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut1)<<endl;
  cout<<"jet55 && evtSel && qalCut2"<<" = "<<ch[2]->GetEntries(jet55&&evtSel&&qalCut2)<<endl;
  cout<<"fakeJet && evtSel"<<" = "<<ch[2]->GetEntries(evtSel&&fakeJet)<<endl;


  ch[2]->Project("hJet55","jtpt",jet55);
  hJet55->Print("base");
  ch[2]->Project("hJet55_QA1","jtpt",evtSel&&jet55&&qalCut1);
  hJet55_QA1->Print("base");
  ch[2]->Project("hJet55_QA2","jtpt",evtSel&&jet55&&qalCut2);
  hJet55_QA2->Print("base");
  ch[2]->Project("hJet55Fake","jtpt",evtSel&&fakeJet);
  hJet55Fake->Print("base");

  ch[2]->Project("hJet55_3D","jteta:jtphi:jtpt",evtSel&&fakeJet);
  hJet55_3D->Print("base");
  
  fout.Write();
  fout.Close();
  
  timer.Stop();
  cout<<"Real time(sec) = "<<timer.RealTime()<<endl;
  cout<<"CPU time(sec)  = "<<timer.CpuTime()<<endl;
  cout<<"All done"<<endl;

}
