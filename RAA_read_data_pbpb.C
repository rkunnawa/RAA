// Raghav Kunnawalkam Elayavalli
// June 23rd 2014
// CERN

// also now that we have the new Jet55or65 dataset, need the lumimask to find out what lumi that particular trigger saw which is necessary when adding the triggers. 

// Now we have to make the macro able to run on condor jobs so we have to split the required files into 11 jobs to match the jet55or65 file list. 
// for this I have to create an array for the jet80 file showing the no of events and each job will go from one array element to the other. which means that i have to stop doing the ttree project method. 
// I would have to think of another way to do that for the project. 

// the best way i think would be to just use the unmerged files for the condor script. that way once we have a TChain of all the files everything below would work. 

// I have the lumi mask. will get the required numbers from that - which are as follows:  
// CORRECT Lumi seen by the respective triggers 
// Jet65 - 139.571 ub-1
// Jet55 - 75.086 ub-1 
// Jet80 - 149.427 ub-1 

// And it looks like we have to separate the pp side from this macro due to condor job submission - Done
// July 25 - add in the capability to read different radii and eta widths in the same macro. 

// Oct 21 - Now that we have the new data forest files which are with made with the triggers: jet55 or jet65 or jet80 we only need one file loop.

// OCt 27th - still have to run the lumi calculator for the triggers using the lumi mask files. Since we lost some files during the processing (~1%).

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

//static const int nbins_pt = 29;
//static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

//only for the jet80 merged file - not necessary now 
//static const int job_no = 11;
//static const double boundaries_loopno_job[job_no+1] = {0,100000,200000,300000,400000,50000,600000,700000,800000,900000,1000000,1152308};

//static const double boundaries_fileno_job[job_no+1] = {0, 413, 826, 1239, 1652, 2065, 2478, 2891, 3304, 3717, 4130, 4542};

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


static const int nbins_eta = 15;
static const double boundaries_eta[nbins_eta][2] = {
  {-1.0,+1.0}, {-2.0,+2.0}, {-3.0,+3.0},
  {-3.0,-2.5}, {-2.5,-2.0}, {-2.0,-1.5}, 
  {-1.5,-1.0}, {-1.0,-0.5}, {-0.5,0}, {0,+0.5}, 
  {+0.5,+1.0}, {+1.0,+1.5}, {+1.5,+2.0}, 
  {+2.0,+2.5}, {+2.5,+3.0}
};

static const double delta_eta[nbins_eta] = {
  2.0, 4.0, 6.0, 
  0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 0.5, 
  0.5, 0.5, 0.5, 
  0.5, 0.5
};

static const char etaWidth [nbins_eta][256] = {
  "n10_eta_p10","n20_eta_p20","n30_eta_p30",
  "n30_eta_n25","n25_eta_n20","n20_eta_n15",
  "n15_eta_n10","n10_eta_n05","n05_eta_0","0_eta_p05",
  "p05_eta_p10","p10_eta_p15","p15_eta_p20",
  "p20_eta_p25","p25_eta_p30"
};
/*

static const int nbins_eta = 1;
static const double boundaries_eta[nbins_eta][2] = {
  {-2.0,+2.0}
};

static const double delta_eta[nbins_eta] = {
  4.0
};

static const char etaWidth [nbins_eta][256] = {
  "n20_eta_p20"
};

*/
static const int nbins_cent = 6;
static Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 5 to get your actual centrality
static Double_t ncoll[nbins_cent+1] = { 1660, 1310, 745, 251, 62.8, 10.8 ,362.24}; //last one is for 0-200 bin. 

static const int no_radius = 3;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {2,3,4};

//static const int no_radius = 1;//necessary for the RAA analysis  
//static const int list_radius[no_radius] = {3};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 7; 
//static const int list_radius[no_radius] = {1,2,3,4,5,6,7};


// divide by bin width
void divideBinWidth(TH1 *h)
{
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    if(val!=0){
      val/=h->GetBinWidth(i);
      valErr/=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }  
  }
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

class JetData
{
public:
  JetData(char *fileName, char *jetTree, char *genJetTree, bool loadGenJet = 0,bool isPbPb = 0) {
    cout <<"Open "<<fileName<<endl;
    tFile = new TFile(fileName,"read");
    tEvt = (TTree*)tFile->Get("hiEvtAnalyzer/HiTree");
    tSkim = (TTree*)tFile->Get("skimanalysis/HltTree");
    tJet = (TTree*)tFile->Get(jetTree);
    tJet->SetBranchAddress("jtpt" , jtpt );
    tJet->SetBranchAddress("trackMax" , trackMax );
    tJet->SetBranchAddress("chargedMax",chargedMax);
    tJet->SetBranchAddress("refpt", refpt);
    tJet->SetBranchAddress("nref" ,&njets);
    tJet->SetBranchAddress("jteta", jteta);
    tJet->SetBranchAddress("jtm",jtmass);
    tJet->SetBranchAddress("pthat",&pthat);
    if (loadGenJet) tGenJet = (TTree*)tFile->Get(genJetTree);
    if (loadGenJet) tGenJet->SetBranchAddress("ngen" ,&ngen);
    if (loadGenJet) tGenJet->SetBranchAddress("genpt", genpt);
    if (loadGenJet) tGenJet->SetBranchAddress("gensubid", gensubid);
    tEvt->SetBranchAddress("hiBin",&bin);
    tEvt->SetBranchAddress("vz",&vz);
    tSkim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    if(isPbPb) tSkim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    else tSkim->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
    tJet->AddFriend(tEvt);
    tJet->AddFriend(tSkim);
  };
  TFile *tFile;
  TTree *tJet;
  TTree *tGenJet;
  TTree *tEvt;
  TTree* tSkim;
  float jtpt[1000];
  float refpt[1000];
  float jteta[1000];
  float jtmass[1000];
  float trackMax[1000];
  float chargedMax[1000];
  float genpt[1000];
  int gensubid[1000];
  float vz;
  float pthat;
  int njets;
  int ngen;
  int bin;     
  int pHBHENoiseFilter;
  int pPAcollisionEventSelectionPA;
  int pcollisionEventSelection;
};

int findBin(int hiBin){
  int binNo = 0;

  for(int i = 0;i<nbins_cent;i++){
    if(hiBin>=5*boundaries_cent[i] && hiBin<5*boundaries_cent[i+1]) {
      binNo = i;
      break;
    }
  }

  return binNo;
}


using namespace std;

void RAA_read_data_pbpb(int startfile = 0, int endfile = 1, char *algo = "Vs", char *jet_type = "PF"){

  TH1::SetDefaultSumw2();
  //gStyle->SetOptStat(0);
  
  TStopwatch timer;
  timer.Start();

  TDatime date;
  
  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = false;
  
  // number convension:
  // 0 - MB
  // 1 - 55 or 65
  // 2 - 80 or 95 
  // 80 is the unprescaled trigger - yes
  // 
  
  // Now im going to change the file reading here for PbPb to look at the unmerged files through condor. 
  std::string infile1;
  //infile1 = "jet55or65_filelist.txt";
  infile1 = "jetRAA_PbPb_data_forest.txt";
  
  //std::string infile2;
  //infile2 = "jet80_filelist.txt";

  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;

  //std::ifstream instr2(infile2.c_str(),std::ifstream::in);
  //std::string filename2;

  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }

  //for(int ifile = 0;ifile<boundaries_fileno_job[startfile];ifile++){
  //  instr2>>filename2;
  //}
 
  const int N = 5;
    
  TChain *jetpbpb1[N][no_radius];
  //TChain *jetpbpb2[N][no_radius];
  
  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%d%sJetAnalyzer",algo,list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    dir[4][k] = "hltobject";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree",
    "jetObjTree"
  };
  
  //this loop is to assign the tree values before we go into the file loop. 
  for(int k = 0;k<no_radius;k++){
    for(int t = 0;t<N;t++){
      jetpbpb1[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data());
    }//tree loop ends
  }// radius loop ends
  
  for(int ifile = startfile;ifile<endfile;ifile++){
    
    instr1>>filename1;
    if(printDebug)cout<<"File: "<<filename1<<endl;
    for(int k = 0;k<no_radius;k++){

      for(int t = 0;t<N;t++){
	
	//jetpbpb1[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data()) ;
	jetpbpb1[t][k]->Add(filename1.c_str());
	if(printDebug)cout << "Tree loaded  " << string(dir[t][k]+"/"+trees[t]).data() << endl;
	if(printDebug)cout << "Entries : " << jetpbpb1[t][k]->GetEntries() << endl;
	
      }// tree loop ends
      
      //jetpbpb1->Add(filename.c_str());
      //hltpbpb1->Add(filename.c_str());
      //skmpbpb1->Add(filename.c_str());
      //hltobjpbpb1->Add(filename.c_str());
      //evtpbpb1->Add(filename.c_str());
      //cout<<"Entries = "<<jetpbpb1->GetEntries()<<endl;
    }// radius loop ends
    
  }// file loop ends
  
  

  // for(int k = 0;k<no_radius;k++){
  //   for(int t = 0;t<N;t++){
  //     jetpbpb2[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data());
  //   }//tree loop ends
  // }// radius loop ends

  // for(int ifile = boundaries_fileno_job[startfile];ifile<boundaries_fileno_job[endfile];ifile++){
  // //for(int ifile = 0;ifile<10;ifile++){
    
  //   instr2>>filename2;
  //   if(printDebug)cout<<"File: "<<filename2<<endl;

  //   for(int k = 0;k<no_radius;k++){

  //     for(int t = 0;t<N;t++){

  //     //jetpbpb2[i] = new TChain(string(dir[i]+"/"+trees[i]).data()) ;
  //     jetpbpb2[t][k]->Add(filename2.c_str());
  //     //cout << "Tree loaded  " << string(dir[i]+"/"+trees[i]).data() << endl;
  //     //cout << "Entries : " << jetpbpb2[i]->GetEntries() << endl;

  //     }//tree loop ends
  //   }// radius loop ends
  //   //jetpbpb1->Add(filename.c_str());
  //   //hltpbpb1->Add(filename.c_str());
  //   //skmpbpb1->Add(filename.c_str());
  //   //hltobjpbpb1->Add(filename.c_str());
  //   //evtpbpb1->Add(filename.c_str());
  //   //cout<<"Entries = "<<jetpbpb1->GetEntries()<<endl;
  // }// file loop ends

  for(int k = 0;k<no_radius;k++){
    jetpbpb1[2][k]->AddFriend(jetpbpb1[0][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[1][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[3][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[4][k]);

    // jetpbpb2[2][k]->AddFriend(jetpbpb2[0][k]);
    // jetpbpb2[2][k]->AddFriend(jetpbpb2[1][k]);
    // jetpbpb2[2][k]->AddFriend(jetpbpb2[3][k]);
    // jetpbpb2[2][k]->AddFriend(jetpbpb2[4][k]);

  }// radius loop ends
 
  // Ok this should now work. 

  cout<<"total no of entries in the combined forest files = "<<jetpbpb1[2][0]->GetEntries()<<endl;
  //  if(printDebug)cout<<"total no of entries in the Jet80 Tree     = "<<jetpbpb2[2][0]->GetEntries()<<endl;

  //these were for doing it from the forests directly without the proper JEC's 
  //add the centrality cuts: 

  //const int nbins_cent = 1;
  //Double_t boundaries_cent[nbins_cent+1] = {0,40};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
  //Double_t ncoll[nbins_cent] = { 362.24};


  //Double_t jet55or65CentWeight[nbins_cent] = {0.3734,0.2509,0.3222,0.0352,0.0066,0.0010}; //total = 0.9983
  //Double_t jet80or90CentWeight[nbins_cent] = {0.2334,0.1764,0.3601,0.1117,0.0283,0.0054}; //total = 0.9153


  // ok so this is pretty important here: 
  // the structure of the macro realies heavily on the centrality loop. so histograms arrays from 0 to nbins_cent-1 will have the spectra for the pbpb at different centrality classes. the one at nbins_cent contains the spectra for the 0-200 bin. so its the full spectra. this was done as a cross check rather than just adding the other histograms. 

  /*
  
  TCut pbpb0 = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIZeroBiasPizel_SingleTrack_v1&&chargedMax/jtpt>0.01";//this is just for the MB file. not really used here so far. 

  TCut pbpb1[nbins_cent+1];
  TCut pbpb2[nbins_cent+1];
  TCut pbpb3[nbins_cent+1];

  TCut pbpb80[nbins_cent+1];
  TCut pbpb65[nbins_cent+1];
  TCut pbpb55[nbins_cent+1];

  TH1F *hpbpb1[nbins_cent+1],*hpbpb2[nbins_cent+1],*hpbpb3[nbins_cent+1];
  TH1F *hpbpbComb[nbins_cent+1];
  //TH1F* htest = new TH1F("htest","",1000,0,1000);
  TH1F *hpbpb_80[nbins_cent+1],*hpbpb_65[nbins_cent+1],*hpbpb_55[nbins_cent+1]; //histos to check the separate spectra, weighted by event by event prescl 
  // I should also add the trigger objects merging method. 

  //old way of finding trigger turn on which didnt work since we dont have a good MB sample. 
  
  TH1F* hTurnon80_old = new TH1F("hTurnon80_old","",150,0,150);
  TH1F* hTurnon65_old = new TH1F("hTurnon65_old","",150,0,150);
  TH1F* hTurnon55_old = new TH1F("hTurnon55_old","",150,0,150);

  TH1F* hTriggerMerged_old = new TH1F("hTriggerMerged_old","",150,0,150);
  TH1F* htest80_old = new TH1F("htest80_old","",150,0,150);
  TH1F* htest65_old = new TH1F("htest65_old","",150,0,150);
  TH1F* htest55_old = new TH1F("htest55_old","",150,0,150);

  // check the trigger turn on curve from the MB file. 
  TH1F* hMB_old = new TH1F("hMB_old","",150,0,150);
  
  //jetpbpb2->Project("htest80","jtpt","HLT_Hjentryet80_v1");
  //htest80->Print("base");
  //jetpbpb1->Project("htest65","jtpt","HLT_Hjentryet65_v1&&!HLT_Hjentryet80_v1");
  //htest65->Print("base");
  //jetpbpb1->Project("htest55","jtpt","HLT_Hjentryet55_v1_Prescl*(HLT_Hjentryet55_v1&&!HLT_Hjentryet65_v1&&!HLT_Hjentryet80_v1)");
  //htest55->Print("base");
  
  jetpbpb0_old->Project("hTurnon80_old","jtpt","HLT_Hjentryet80_v1_Prescl*HLT_Hjentryet80_v1");
  jetpbpb0_old->Project("hTurnon65_old","jtpt","HLT_Hjentryet65_v1_Prescl*HLT_Hjentryet65_v1");
  jetpbpb0_old->Project("hTurnon55_old","jtpt","HLT_Hjentryet55_v1_Prescl*HLT_Hjentryet55_v1");

  TCut MB_prescl = "HLT_HIMinBiasHfOrBSC_v1_Prescl*HLT_HIMinBiasHfOrBSC_v1";
  jetpbpb0_old->Project("hMB_old","jtpt","30"*MB_prescl);

  //hTurnon80->Print("base");
  //hTurnon65->Print("base");
  //hTurnon55->Print("base");

  //hTriggerMerged->Add(htest80);
  //hTriggerMerged->Add(htest65);
  //hTriggerMerged->Add(htest55);

  hTurnon80_old->Divide(hMB_old);
  hTurnon65_old->Divide(hMB_old);
  hTurnon55_old->Divide(hMB_old);
  

  
  //centrality loop for the pbpb files/histograms 
  for(int i = 0;i<nbins_cent;i++){

    cout<<"centrality boundary = "<<boundaries_cent[i]*5<<" - "<<boundaries_cent[i+1]*5<<endl;

   
    //  Ok so i screwed up the convention here. the tcuts and the histogram has the number appendage going from 
    //  1 - 80
    //  2 - 65
    //  3 - 55

    //  but thats the reverse for the jet trees. 
    //  2 - HLT_80 or HLT_95
    //  1 - HLT_55 or HLT_65

      so be careful when moving through the code. 
    

    //list of cuts we are using to remove all the fake jets. 
    // chMax/jtpt>0.01 - studied
    // neMax/jtpt>0.01 - just what im trying out now - needs to be studied
    // TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975 used by kurt in 14007 and in 12003 as well.  
    // 

    // lets talk about the actual events/jets which are irritating here. the so called supernova events. 
    // they come from smaller jet triggers but still have larger jet pt. 

    pbpb1[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb2[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet65_v1&&!HLT_Hjentryet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb3[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet55_v1&&!HLT_Hjentryet65_v1&&!HLT_Hjentryet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);

    pbpb80[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb65[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet65_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb55[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet55_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);

    hpbpb1[i] = new TH1F(Form("hpbpb1_cent%d",i),Form("spectra from Jet 80 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    //hpbpb1[i]->Print("base");
    hpbpb2[i] = new TH1F(Form("hpbpb2_cent%d",i),Form("spectra from jet 65 & !jet 80 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb3[i] = new TH1F(Form("hpbpb3_cent%d",i),Form("spectra from jet 55 & !jet 65 & !jet 80 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpbComb[i] = new TH1F(Form("hpbpbComb_cent%d",i),Form("Spectra Combined using 12003 method %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_80[i] = new TH1F(Form("hpbpb_80_cent%d",i),Form("Spectra from Jet80 trigger alone %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_65[i] = new TH1F(Form("hpbpb_65_cent%d",i),Form("Spectra from Jet65 trigger alone %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_55[i] = new TH1F(Form("hpbpb_55_cent%d",i),Form("Spectra from Jet55 trigger alone %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    jetpbpb2[2]->Project(Form("hpbpb1_cent%d",i),"jtpt",pbpb1[i]);
    hpbpb1[i]->Print("base");
    //divideBinWidth(hpbpb1);
    
    jetpbpb1[2]->Project(Form("hpbpb2_cent%d",i),"jtpt",pbpb2[i]);
    hpbpb2[i]->Print("base");
    //divideBinWidth(hpbpb2);
    
    jetpbpb1[2]->Project(Form("hpbpb3_cent%d",i),"jtpt","HLT_Hjentryet55_v1_Prescl"*pbpb3[i]);
    //jetpbpb1[2]->Project("hpbpb3","jtpt","2.34995"*pbpb3);
    hpbpb3[i]->Print("base");
    //divideBinWidth(hpbpb3);

    jetpbpb2[2]->Project(Form("hpbpb_80_cent%d",i),"jtpt","HLT_Hjentryet80_v1_Prescl"*pbpb80[i]);
    hpbpb_80[i]->Print("base");
    jetpbpb1[2]->Project(Form("hpbpb_65_cent%d",i),"jtpt","HLT_Hjentryet65_v1_Prescl"*pbpb65[i]);
    hpbpb_65[i]->Print("base");
    jetpbpb1[2]->Project(Form("hpbpb_55_cent%d",i),"jtpt","HLT_Hjentryet55_v1_Prescl"*pbpb55[i]);
    hpbpb_55[i]->Print("base");

    //scale the PbPb histograms before adding them
    //we have to scale them according to the lumi of the Jet80 file. 
    // HLT file  |   Lumi inverse micro barns 
    // HLT_80    |   149.382 
    // HLT_65    |   3.195
    // HLT_55    |   2.734
    // 

    // the files which im using now is only a fraction of events of that:
    // for the PbPb 55 or 65 file its 0.977
    // for the PbPb 80 file its 0.304
    // now using the full sample file 
    
    hpbpb1[i]->Scale(1./149.382e6);//respective lumi seen by the trigger all in inverse micro barns 
    hpbpb2[i]->Scale(1./3.195e6);
    hpbpb3[i]->Scale(1./2.734e6);
    
    hpbpb1[i]->Scale(1./4);//delta eta
    hpbpb2[i]->Scale(1./4);
    hpbpb3[i]->Scale(1./4);

    //hpbpb1[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i])/jet80or95CentWeight[i]);//centrality bin width and scaling by the centrality events fraction. 
    //hpbpb2[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i])/jet55or65CentWeight[i]);
    //hpbpb3[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i])/jet55or65CentWeight[i]);

    hpbpb1[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));//centrality bin width 
    hpbpb2[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));
    hpbpb3[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));

    //might have to end up adding a centrality weight - the ratio of events per centrality class. 
    
    //add the histograms  
    hpbpbComb[i]->Add(hpbpb1[i]);
    hpbpbComb[i]->Add(hpbpb2[i]);
    hpbpbComb[i]->Add(hpbpb3[i]);
    hpbpbComb[i]->Print("base");

    hpbpbComb[i] = (TH1F*)hpbpbComb[i]->Rebin(nbins_pt,Form("hpbpbComb_cent%d",i),boundaries_pt);
    hpbpb3[i] = (TH1F*)hpbpb3[i]->Rebin(nbins_pt,Form("hpbpb3_cent%d",i),boundaries_pt);
    hpbpb2[i] = (TH1F*)hpbpb2[i]->Rebin(nbins_pt,Form("hpbpb2_cent%d",i),boundaries_pt);
    hpbpb1[i] = (TH1F*)hpbpb1[i]->Rebin(nbins_pt,Form("hpbpb1_cent%d",i),boundaries_pt);

    divideBinWidth(hpbpbComb[i]);
    divideBinWidth(hpbpb3[i]);
    divideBinWidth(hpbpb2[i]);
    divideBinWidth(hpbpb1[i]);
    
  }

  // doing it the 0-200 centrality bin 

  pbpb1[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet80_v1&&(chargedMax/jtpt)>0.01&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975";
  pbpb2[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet65_v1&&!HLT_Hjentryet80_v1&&(chargedMax/jtpt)>0.01&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975";
  pbpb3[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet55_v1&&!HLT_Hjentryet65_v1&&!HLT_Hjentryet80_v1&&(chargedMax/jtpt)>0.01&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975";
  
  pbpb80[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet80_v1&&(chargedMax/jtpt>0.01)&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975";
  pbpb65[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet65_v1&&(chargedMax/jtpt>0.01&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975";
  pbpb55[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_Hjentryet55_v1&&(chargedMax/jtpt)>0.01&&TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975";

  hpbpb1[nbins_cent] = new TH1F(Form("hpbpb1_cent%d",nbins_cent),"Spectra from Jet80 0-200 cent",1000,0,1000);
  //hpbpb1[i]->Print("base");
  hpbpb2[nbins_cent] = new TH1F(Form("hpbpb2_cent%d",nbins_cent),"Spectra from Jet 65 & !Jet80 0-200 cent",1000,0,1000);
  hpbpb3[nbins_cent] = new TH1F(Form("hpbpb3_cent%d",nbins_cent),"Spectra from Jet 55 & !Jet65 & !Jet80 0-200 cent",1000,0,1000);
  hpbpbComb[nbins_cent] = new TH1F(Form("hpbpbComb_cent%d",nbins_cent),"Combined Jet spectra 12003 method 0-200 cent",1000,0,1000);
  
  hpbpb_80[nbins_cent] = new TH1F(Form("hpbpb_80_cent%d",nbins_cent),"Spectra from Jet 80 alone 0-200 cent",1000,0,1000);
  hpbpb_65[nbins_cent] = new TH1F(Form("hpbpb_65_cent%d",nbins_cent),"Spectra from Jet 65 alone 0-200 cent",1000,0,1000);
  hpbpb_55[nbins_cent] = new TH1F(Form("hpbpb_55_cent%d",nbins_cent),"Spectra from Jet 55 alone 0-200 cent",1000,0,1000);
  
  jetpbpb2[2]->Project(Form("hpbpb1_cent%d",nbins_cent),"jtpt",pbpb1[nbins_cent]);
  hpbpb1[nbins_cent]->Print("base");
  //divideBinWidth(hpbpb1);
    
  jetpbpb1[2]->Project(Form("hpbpb2_cent%d",nbins_cent),"jtpt",pbpb2[nbins_cent]);
  hpbpb2[nbins_cent]->Print("base");
  //divideBinWidth(hpbpb2);
  
  jetpbpb1[2]->Project(Form("hpbpb3_cent%d",nbins_cent),"jtpt","HLT_Hjentryet55_v1_Prescl"*pbpb3[nbins_cent]);
  //jetpbpb1[2]->Project("hpbpb3","jtpt","2.34995"*pbpb3);
  hpbpb3[nbins_cent]->Print("base");
  //divideBinWidth(hpbpb3);
  
  //following histograms are for the trigger turnon curve. no cuts apart from the trigger selection. 
  jetpbpb2[2]->Project(Form("hpbpb_80_cent%d",nbins_cent),"jtpt","HLT_Hjentryet80_v1_Prescl*HLT_Hjentryet80_v1");
  hpbpb_80[nbins_cent]->Print("base");
  jetpbpb1[2]->Project(Form("hpbpb_65_cent%d",nbins_cent),"jtpt","HLT_Hjentryet65_v1_Prescl*HLT_Hjentryet65_v1");
  hpbpb_65[nbins_cent]->Print("base");
  jetpbpb1[2]->Project(Form("hpbpb_55_cent%d",nbins_cent),"jtpt","HLT_Hjentryet55_v1_Prescl*HLT_Hjentryet55_v1");
  hpbpb_55[nbins_cent]->Print("base");
  
  //scale the PbPb histograms before adding them
  //we have to scale them according to the lumi of the Jet80 file. 
  // HLT file  |   Lumi inverse micro barns 
  // HLT_80    |   149.382 
  // HLT_65    |   3.195
  // HLT_55    |   2.734
  // 
  
  // the files which im using now is only a fraction of events of that:
  // for the PbPb 55 or 65 file its 0.977
  // for the PbPb 80 file its 0.304
  // now using the full sample file 
  
  hpbpb1[nbins_cent]->Scale(1./149.382e6);//respective lumi seen by the trigger all in inverse micro barns 
  hpbpb2[nbins_cent]->Scale(1./3.195e6);
  hpbpb3[nbins_cent]->Scale(1./2.734e6);
  
  hpbpb1[nbins_cent]->Scale(1./4);//delta eta
  hpbpb2[nbins_cent]->Scale(1./4);
  hpbpb3[nbins_cent]->Scale(1./4);
  
  //might have to end up adding a centrality weight - the ratio of events per centrality class. 
  
  //add the histograms  
  hpbpbComb[nbins_cent]->Add(hpbpb1[nbins_cent]);
  hpbpbComb[nbins_cent]->Add(hpbpb2[nbins_cent]);
  hpbpbComb[nbins_cent]->Add(hpbpb3[nbins_cent]);
  hpbpbComb[nbins_cent]->Print("base");
  
  hpbpbComb[nbins_cent] = (TH1F*)hpbpbComb[nbins_cent]->Rebin(nbins_pt,Form("hpbpbComb_cent%d",nbins_cent),boundaries_pt);
  hpbpb3[nbins_cent] = (TH1F*)hpbpb3[nbins_cent]->Rebin(nbins_pt,Form("hpbpb3_cent%d",nbins_cent),boundaries_pt);
  hpbpb2[nbins_cent] = (TH1F*)hpbpb2[nbins_cent]->Rebin(nbins_pt,Form("hpbpb2_cent%d",nbins_cent),boundaries_pt);
  hpbpb1[nbins_cent] = (TH1F*)hpbpb1[nbins_cent]->Rebin(nbins_pt,Form("hpbpb1_cent%d",nbins_cent),boundaries_pt);
  
  divideBinWidth(hpbpbComb[nbins_cent]);
  divideBinWidth(hpbpb3[nbins_cent]);
  divideBinWidth(hpbpb2[nbins_cent]);
  divideBinWidth(hpbpb1[nbins_cent]);
  
  //ok now we have the spectra for the 0-200% centrality. 
  
  */

  
  // do the trigger object merging here: 
  // this has to be done in the event loop which means that we have to get the 
  // create the trees and set the branch address
  // jet tree

  // similarly here 0 - MB file, 1 - 55or65, 2 - 80or95
  
//   //file 0:
//   // jet tree
//   int nrefe_0;
//   float pt_0[1000];
//   //float old_pt3[1000];
//   float raw_0[1000];
//   float eta_0[1000];
//   float eta_0_CM[1000];
//   float phi_0[1000];
//   float chMax_0[1000];
//   float trkMax_0[1000];
//   float chSum_0[1000];
//   float phSum_0[1000];
//   float neSum_0[1000];
//   float trkSum_0[1000];
//   float phMax_0[1000];
//   float neMax_0[1000];

//   // event tree
//   int evt_0;
//   int run_0;
//   int lumi_0;
//   int hiBin_0;
//   float vx_0;
//   float vy_0;
//   float vz_0;
//   int hiNtracks_0;
//   float hiHFminus_0;
//   float hiHFplus_0;
//   float hiHFplusEta4_0;
//   float hiHFminusEta4_0;
//   int pcollisionEventSelection_0;
//   int pHBHENoiseFilter_0;
//   int pprimaryvertexFilter_0;
//   int pVertexFilterCutGplus_0;

//   // trigger tree
//   int L1_MB_0;
//   int L1_MB_p_0;
//   int jetMB_0;
//   int jet55_0;
//   int jet65_0;
//   int jet80_0;
//   int jetMB_p_0;
//   int jet55_p_0;
//   int jet65_p_0;
//   int jet80_p_0;

//   // trigger object tree - this contains the maximum value of the particular trigger object. 
//   float trgObj_id_0;
//   float trgObj_pt_0;
//   float trgObj_eta_0;
//   float trgObj_phi_0;
//   float trgObj_mass_0;


  //file 1: 
  // jet tree
  int nrefe_1;
  float pt_1[1000];
  //float old_pt3[1000];
  float raw_1[1000];
  float eta_1[1000];
  float eta_1_CM[1000];
  float phi_1[1000];
  float chMax_1[1000];
  float trkMax_1[1000];
  float chSum_1[1000];
  float phSum_1[1000];
  float neSum_1[1000];
  float trkSum_1[1000];
  float phMax_1[1000];
  float neMax_1[1000];

  // event tree
  int evt_1;
  int run_1;
  int lumi_1;
  int hiBin_1;
  float vx_1;
  float vy_1;
  float vz_1;
  int hiNtracks_1;
  float hiHFminus_1;
  float hiHFplus_1;
  float hiHFplusEta4_1;
  float hiHFminusEta4_1;
  int pcollisionEventSelection_1;
  int pHBHENoiseFilter_1;
  int pprimaryvertexFilter_1;
  int pVertexFilterCutGplus_1;

  // trigger tree
  int L1_MB_1;
  int L1_MB_p_1;
  int jetMB_1;
  int jet55_1;
  int jet65_1;
  int jet80_1;
  int jetMB_p_1;
  int jet55_p_1;
  int jet65_p_1;
  int jet80_p_1;


  // trigger object tree - this contains the maximum value of the particular trigger object. 
  float trgObj_id_1;
  float trgObj_pt_1;
  float trgObj_eta_1;
  float trgObj_phi_1;
  float trgObj_mass_1;
  
  // //file 2: 
  // // jet tree
  // int nrefe_2;
  // float pt_2[1000];
  // //float old_pt3[1000];
  // float raw_2[1000];
  // float eta_2[1000];
  // float eta_2_CM[1000];
  // float phi_2[1000];
  // float chMax_2[1000];
  // float trkMax_2[1000];
  // float chSum_2[1000];
  // float phSum_2[1000];
  // float neSum_2[1000];
  // float trkSum_2[1000];
  // float phMax_2[1000];
  // float neMax_2[1000];

  // // event tree
  // int evt_2;
  // int run_2;
  // int lumi_2;
  // int hiBin_2;
  // float vx_2;
  // float vy_2;
  // float vz_2;
  // int hiNtracks_2;
  // float hiHFminus_2;
  // float hiHFplus_2;
  // float hiHFplusEta4_2;
  // float hiHFminusEta4_2;
  // int pcollisionEventSelection_2;
  // int pHBHENoiseFilter_2;
  // int pprimaryvertexFilter_2;
  // int pVertexFilterCutGplus_2;

  // // trigger tree
  // int L1_MB_2;
  // int L1_MB_p_2;
  // int jetMB_2;
  // int jet55_2;
  // int jet65_2;
  // int jet80_2;
  // int jetMB_p_2;
  // int jet55_p_2;
  // int jet65_p_2;
  // int jet80_p_2;


  // // trigger object tree - this contains the maximum value of the particular trigger object. 
  // float trgObj_id_2;
  // float trgObj_pt_2;
  // float trgObj_eta_2;
  // float trgObj_phi_2;
  // float trgObj_mass_2;
  
  
//   //set the branch addresses:  - one of the most boring parts of the code: 
//   jetpbpb0_old->SetBranchAddress("evt",&evt_0);
//   jetpbpb0_old->SetBranchAddress("run",&run_0);
//   jetpbpb0_old->SetBranchAddress("lumi",&lumi_0);
//   jetpbpb0_old->SetBranchAddress("hiBin",&hiBin_0);
//   jetpbpb0_old->SetBranchAddress("vz",&vz_0);
//   jetpbpb0_old->SetBranchAddress("vx",&vx_0);
//   jetpbpb0_old->SetBranchAddress("vy",&vy_0);
//   jetpbpb0_old->SetBranchAddress("hiNtracks",&hiNtracks_0);
//   jetpbpb0_old->SetBranchAddress("hiHFminus",&hiHFminus_0);
//   jetpbpb0_old->SetBranchAddress("hiHFplus",&hiHFplus_0);
//   jetpbpb0_old->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_0);
//   jetpbpb0_old->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_0);
//   jetpbpb0_old->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_0);
//   jetpbpb0_old->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_0);
//   //jetpbpb0_old->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_0);
//   //jetpbpb0_old->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_0);
  
//   jetpbpb0_old->SetBranchAddress("nref",&nrefe_0);
//   jetpbpb0_old->SetBranchAddress("jtpt",&pt_0);
//   jetpbpb0_old->SetBranchAddress("jteta",&eta_0);
//   jetpbpb0_old->SetBranchAddress("jtphi",&phi_0);
//   jetpbpb0_old->SetBranchAddress("rawpt",&raw_0);
//   jetpbpb0_old->SetBranchAddress("chargedMax",&chMax_0);
//   jetpbpb0_old->SetBranchAddress("chargedSum",&chSum_0);
//   jetpbpb0_old->SetBranchAddress("trackMax",&trkMax_0);
//   jetpbpb0_old->SetBranchAddress("trackSum",&trkSum_0);
//   jetpbpb0_old->SetBranchAddress("photonMax",&phMax_0);
//   jetpbpb0_old->SetBranchAddress("photonSum",&phSum_0);
//   jetpbpb0_old->SetBranchAddress("neutralMax",&neMax_0);
//   jetpbpb0_old->SetBranchAddress("neutralSum",&neSum_0);
//   /*
//     jetTree->SetBranchAddress("nTrk",&nTrack);
//     jetTree->SetBranchAddress("trkPt",&trkPt);
//     jetTree->SetBranchAddress("trkEta",&trkEta);
//     jetTree->SetBranchAddress("trkPhi",&trkPhi);
//     jetTree->SetBranchAddress("highPurity",&highPurity);
//     jetTree->SetBranchAddress("trkDz1",&trkDz1);
//     jetTree->SetBranchAddress("trkDzError1",&trkDzError1);
//     jetTree->SetBranchAddress("trkDxy1",&trkDxy1);
//     jetTree->SetBranchAddress("trkDxyError1",&trkDxyError1);
//   */
//   //jetpbpb0_old->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_0);
//   //jetpbpb0_old->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_0);
//   //jetpbpb0_old->SetBranchAddress("L1_ZeroBias",&L1_MB_0);
//   //jetpbpb0_old->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_0);
//   jetpbpb0_old->SetBranchAddress("HLT_Hjentryet55_v1",&jet55_0);
//   jetpbpb0_old->SetBranchAddress("HLT_Hjentryet55_v1_Prescl",&jet55_p_0);
//   jetpbpb0_old->SetBranchAddress("HLT_Hjentryet65_v1",&jet65_0);
//   jetpbpb0_old->SetBranchAddress("HLT_Hjentryet65_v1_Prescl",&jet65_p_0);
//   jetpbpb0_old->SetBranchAddress("HLT_Hjentryet80_v1",&jet80_0);
//   jetpbpb0_old->SetBranchAddress("HLT_Hjentryet80_v1_Prescl",&jet80_p_0);
//   /*
//   jetpbpb0_old->SetBranchAddress("id",&trgObj_id_0);
//   jetpbpb0_old->SetBranchAddress("pt",&trgObj_pt_0);
//   jetpbpb0_old->SetBranchAddress("eta",&trgObj_eta_0);
//   jetpbpb0_old->SetBranchAddress("phi",&trgObj_phi_0);
//   jetpbpb0_old->SetBranchAddress("mass",&trgObj_mass_0);
//   */

  for(int k = 0;k<no_radius;k++){
    //set the branch addresses:  - one of the most boring parts of the code: 
    jetpbpb1[2][k]->SetBranchAddress("evt",&evt_1);
    jetpbpb1[2][k]->SetBranchAddress("run",&run_1);
    jetpbpb1[2][k]->SetBranchAddress("lumi",&lumi_1);
    jetpbpb1[2][k]->SetBranchAddress("hiBin",&hiBin_1);
    jetpbpb1[2][k]->SetBranchAddress("vz",&vz_1);
    jetpbpb1[2][k]->SetBranchAddress("vx",&vx_1);
    jetpbpb1[2][k]->SetBranchAddress("vy",&vy_1);
    jetpbpb1[2][k]->SetBranchAddress("hiNtracks",&hiNtracks_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFminus",&hiHFminus_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFplus",&hiHFplus_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_1);
    jetpbpb1[2][k]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_1);
    jetpbpb1[2][k]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_1);
    //jetpbpb1[2][k]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_1);
    //jetpbpb1[2][k]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_1);
  
    jetpbpb1[2][k]->SetBranchAddress("nref",&nrefe_1);
    jetpbpb1[2][k]->SetBranchAddress("jtpt",&pt_1);
    jetpbpb1[2][k]->SetBranchAddress("jteta",&eta_1);
    jetpbpb1[2][k]->SetBranchAddress("jtphi",&phi_1);
    jetpbpb1[2][k]->SetBranchAddress("rawpt",&raw_1);
    jetpbpb1[2][k]->SetBranchAddress("chargedMax",&chMax_1);
    jetpbpb1[2][k]->SetBranchAddress("chargedSum",&chSum_1);
    jetpbpb1[2][k]->SetBranchAddress("trackMax",&trkMax_1);
    jetpbpb1[2][k]->SetBranchAddress("trackSum",&trkSum_1);
    jetpbpb1[2][k]->SetBranchAddress("photonMax",&phMax_1);
    jetpbpb1[2][k]->SetBranchAddress("photonSum",&phSum_1);
    jetpbpb1[2][k]->SetBranchAddress("neutralMax",&neMax_1);
    jetpbpb1[2][k]->SetBranchAddress("neutralSum",&neSum_1);

    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v1",&jet55_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v1",&jet65_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v1",&jet80_1);
    jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_1);

    jetpbpb1[2][k]->SetBranchAddress("id",&trgObj_id_1);
    jetpbpb1[2][k]->SetBranchAddress("pt",&trgObj_pt_1);
    jetpbpb1[2][k]->SetBranchAddress("eta",&trgObj_eta_1);
    jetpbpb1[2][k]->SetBranchAddress("phi",&trgObj_phi_1);
    jetpbpb1[2][k]->SetBranchAddress("mass",&trgObj_mass_1);

    // //set the branch addresses:  - one of the most boring parts of the code: 
    // jetpbpb2[2][k]->SetBranchAddress("evt",&evt_2);
    // jetpbpb2[2][k]->SetBranchAddress("run",&run_2);
    // jetpbpb2[2][k]->SetBranchAddress("lumi",&lumi_2);
    // jetpbpb2[2][k]->SetBranchAddress("hiBin",&hiBin_2);
    // jetpbpb2[2][k]->SetBranchAddress("vz",&vz_2);
    // jetpbpb2[2][k]->SetBranchAddress("vx",&vx_2);
    // jetpbpb2[2][k]->SetBranchAddress("vy",&vy_2);
    // jetpbpb2[2][k]->SetBranchAddress("hiNtracks",&hiNtracks_2);
    // jetpbpb2[2][k]->SetBranchAddress("hiHFminus",&hiHFminus_2);
    // jetpbpb2[2][k]->SetBranchAddress("hiHFplus",&hiHFplus_2);
    // jetpbpb2[2][k]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_2);
    // jetpbpb2[2][k]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_2);
    // jetpbpb2[2][k]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_2);
    // jetpbpb2[2][k]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_2);
    // //jetpbpb2[2][k]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_2);
    // //jetpbpb2[2][k]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_2);
  
    // jetpbpb2[2][k]->SetBranchAddress("nref",&nrefe_2);
    // jetpbpb2[2][k]->SetBranchAddress("jtpt",&pt_2);
    // jetpbpb2[2][k]->SetBranchAddress("jteta",&eta_2);
    // jetpbpb2[2][k]->SetBranchAddress("jtphi",&phi_2);
    // jetpbpb2[2][k]->SetBranchAddress("rawpt",&raw_2);
    // jetpbpb2[2][k]->SetBranchAddress("chargedMax",&chMax_2);
    // jetpbpb2[2][k]->SetBranchAddress("chargedSum",&chSum_2);
    // jetpbpb2[2][k]->SetBranchAddress("trackMax",&trkMax_2);
    // jetpbpb2[2][k]->SetBranchAddress("trackSum",&trkSum_2);
    // jetpbpb2[2][k]->SetBranchAddress("photonMax",&phMax_2);
    // jetpbpb2[2][k]->SetBranchAddress("photonSum",&phSum_2);
    // jetpbpb2[2][k]->SetBranchAddress("neutralMax",&neMax_2);
    // jetpbpb2[2][k]->SetBranchAddress("neutralSum",&neSum_2);

    // //jetpbpb2[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_2);
    // //jetpbpb2[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_2);
    // //jetpbpb2[2][k]->SetBranchAddress("L1_ZeroBias",&L1_MB_2);
    // //jetpbpb2[2][k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_2);
    // jetpbpb2[2][k]->SetBranchAddress("HLT_HIJet55_v1",&jet55_2);
    // jetpbpb2[2][k]->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_2);
    // jetpbpb2[2][k]->SetBranchAddress("HLT_HIJet65_v1",&jet65_2);
    // jetpbpb2[2][k]->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_2);
    // jetpbpb2[2][k]->SetBranchAddress("HLT_HIJet80_v1",&jet80_2);
    // jetpbpb2[2][k]->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_2);

    // jetpbpb2[2][k]->SetBranchAddress("id",&trgObj_id_2);
    // jetpbpb2[2][k]->SetBranchAddress("pt",&trgObj_pt_2);
    // jetpbpb2[2][k]->SetBranchAddress("eta",&trgObj_eta_2);
    // jetpbpb2[2][k]->SetBranchAddress("phi",&trgObj_phi_2);
    // jetpbpb2[2][k]->SetBranchAddress("mass",&trgObj_mass_2);
  
  }//radius loop

  //now that we have all the branch addresses set up we can start going through the loop to look at the trigger objects 
  
  //before we go and do the trigger object merging lets do the text file information here: 
  //we need one text file per spectra, which means we need 3 -> one for each trigger object 
  //and we also need one for the high pt jets with low trigger objects 
  // and these files will have the following structure: 
  // run lumi evt HLTobjpt HLTobjeta HLTobjphi hibin jtpt jteta jtphi

  //ofstream fHLT_80,fHLT_65,fHLT_55;
  ofstream fHLT_high[no_radius];
  //fHLT_80.open(Form("pbpb_%s_R%d_max_trgObj_80_jets.txt",algo,radius));
  //fHLT_65.open(Form("pbpb_%s_R%d_80_max_trgObj_65_jets.txt",algo,radius));
  //fHLT_55.open(Form("pbpb_%s_R%d_65_max_trgObj_55_jets.txt",algo,radius));
  for(int k = 0;k<no_radius;k++){
    fHLT_high[k].open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/pbpb_%s_R%d_abnormal_jets_%d_%d.txt",algo,list_radius[k],date.GetDate(),endfile));
  }
  // the actual cut we are going to be using is 
  // TMath::Max(chargedMax,neutralMax)/(TMath::Max(chargedSum,neutralSum))<0.975
  // which is pretty much the same as doing the one mentioned above.  
  // ok the cut mentioned above is not the one we are using for the analysis now. We have to check if that cut actually removes a lot of the fake events for us to have a meaningful combined jet spectra. 
  

  // TMath::Max(chargedMax,neutralMax)/TMath::Max(chargedSum,neutralSum)<0.975 im going to add this from Kurt used in 12003 here to see if it sort of helps. he also used this in 14007 analysis. 


  //declare the histograms needed for the hpbpb_TrigObj80, hpbpb_TrigObj65, hpbpb_TrigObj55, hpbpb_TrigObjMB, hpbpb_TrigComb;

  TH1F *hpbpb_TrgObj80[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_TrgObj65[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_TrgObj55[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_TrgObjMB[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_TrgObjComb[no_radius][nbins_eta][nbins_cent+1];
  //test histograms for the spectra alone 
  TH1F *hpbpb_Jet80[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet65[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet55[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet65_noJet80[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_Jet55_noJet65_80[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_JetComb_old[no_radius][nbins_eta][nbins_cent+1];

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      for(int i = 0;i<nbins_cent;i++){

        hpbpb_TrgObj80[k][j][i] = new TH1F(Form("hpbpb_TrgObj80_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from Trig Object > 80 and Jet 80 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
        hpbpb_TrgObj65[k][j][i] = new TH1F(Form("hpbpb_TrgObj65_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from Trig Object > 65 and < 80 and Jet 65 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
        hpbpb_TrgObj55[k][j][i] = new TH1F(Form("hpbpb_TrgObj55_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from Trig Object > 55 and < 65 and Jet 55 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
        hpbpb_TrgObjMB[k][j][i] = new TH1F(Form("hpbpb_TrgObjMB_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from MB file R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
        hpbpb_TrgObjComb[k][j][i] = new TH1F(Form("hpbpb_TrgObjComb_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Trig Combined Spectra using 14007 method R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
	
	hpbpb_Jet80[k][j][i] = new TH1F(Form("hpbpb_Jet80_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from Jet 80 trigger alone R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
	hpbpb_Jet65[k][j][i] = new TH1F(Form("hpbpb_Jet65_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from Jet 65 trigger alone R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
	hpbpb_Jet55[k][j][i] = new TH1F(Form("hpbpb_Jet55_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from Jet 55 trigger alone R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
	hpbpb_Jet65_noJet80[k][j][i] = new TH1F(Form("hpbpb_Jet65_noJet80_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra for Jet65 and not Jet80 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
	hpbpb_Jet55_noJet65_80[k][j][i] = new TH1F(Form("hpbpb_Jet55_noJet65_80_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra for Jet55 and not Jet65 and not Jet80 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
	hpbpb_JetComb_old[k][j][i] = new TH1F(Form("hpbpb_JetComb_oldR%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Combined spectra using old method R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);

      }//cent bin loop

      hpbpb_TrgObj80[k][j][nbins_cent] = new TH1F(Form("hpbpb_TrgObj80_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra from Trig Object > 80 and Jet 80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpbpb_TrgObj65[k][j][nbins_cent] = new TH1F(Form("hpbpb_TrgObj65_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra from Trig Object > 65 and < 80 and Jet 65 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpbpb_TrgObj55[k][j][nbins_cent] = new TH1F(Form("hpbpb_TrgObj55_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra from Trig Object > 55 and < 65 and Jet 55 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpbpb_TrgObjMB[k][j][nbins_cent] = new TH1F(Form("hpbpb_TrgObjMB_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra from MB file R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
      hpbpb_TrgObjComb[k][j][nbins_cent] = new TH1F(Form("hpbpb_TrgObjComb_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Trig combined spectra using 14007 method R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt,boundaries_pt);
    
      hpbpb_Jet80[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet80_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra form Jet80 only R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      hpbpb_Jet65[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet65_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra form Jet65 only R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      hpbpb_Jet55[k][j][nbins_cent] = new TH1F(Form("hpbpb_Jet55_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra form Jet55 only R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);

      hpbpb_Jet65_noJet80[k][j][nbins_cent] = new TH1F(Form("hpbbp_Jet65_noJet80_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra for Jet65 and not Jet80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      hpbpb_Jet55_noJet65_80[k][j][nbins_cent] = new TH1F(Form("hpbbp_Jet55_noJet65_80_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra for Jet55 and not Jet65 and not Jet80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);
      hpbpb_JetComb_old[k][j][nbins_cent] = new TH1F(Form("hpbpb_JetComb_old_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Combined spectra using the old method R%d %s 0-200 cent",list_radius[k],etaWidth[j]),nbins_pt, boundaries_pt);

    }//eta bin loop  

  }//radius loop

  //loop for the MB tree. 
  //empty for now - need to fix that once we have the MB data.  
  
  if(printDebug)hpbpb_TrgObj80[0][0][0]->Print("base");
  
  for(int k = 0;k<no_radius;k++){

    if(printDebug)cout<<"Running data reading for R = "<<list_radius[k]<<endl;
    // loop for the jetpbpb1[2] tree 
    Long64_t nentries_jet55or65 = jetpbpb1[2][k]->GetEntries();
    if(printDebug)cout<<"nentries_jet55or65or80 = "<<nentries_jet55or65<<endl;
    if(printDebug)nentries_jet55or65 = 2;

    for(int jentry = 0;jentry<nentries_jet55or65;jentry++){
    
      jetpbpb1[2][k]->GetEntry(jentry);
      //if(printDebug && jentry%100000==0)cout<<"Jet 55or65 file"<<endl;
      if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_1<<"; run = "<<run_1<<endl;
      
      // get the stuff required for the trigger turn on curve later. in a separate loop till i understand how to put this in here. 

      int centBin = findBin(hiBin_1);//tells us the centrality of the event. 
      if(printDebug)cout<<"cent bin = "<<centBin<<endl;
      if(printDebug)cout<<"centrality bin = "<<5*boundaries_cent[centBin]<< " to "<<5*boundaries_cent[centBin+1]<<endl;

      for(int j = 0;j<nbins_eta;j++){

        if(pHBHENoiseFilter_1 && pcollisionEventSelection_1 && (vz_1<15)) {

	  if(printDebug)cout<<" trigger object pt =  "<<trgObj_pt_1<<endl;
	  if(printDebug)cout<<" jet80 = "<<jet80_1<<endl;
	  if(printDebug)cout<<" jet65 = "<<jet65_1<<endl;
	  if(printDebug)cout<<" jet55 = "<<jet55_1<<endl;
	  if(printDebug)cout<<" nrefe = "<<nrefe_1<<endl;

          for(int g = 0;g<nrefe_1;g++){

            if(/*put your favourite QA cut here*/(chMax_1[g]/pt_1[g]>0.01) /*&& (pt_1[g]<3*trgObj_pt_1)*/){

              if(eta_1[g]>=boundaries_eta[j][0] && eta_1[g]<boundaries_eta[j][1]){
		
                if(jet55_1==1) { // passes the jet55 trigger
		  if(trgObj_pt_1>=55 && trgObj_pt_1<65){ // check for the trigger object pt to lie inbetween the two trigger values 
		    hpbpb_TrgObj55[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
		    hpbpb_TrgObj55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
		    
		    if(pt_1[g]>3*trgObj_pt_1){ // get an idea of the event information from these large pt jets 
		      fHLT_high[k]<<evt_1<<" "<<run_1<<" "<<lumi_1<<" "<<vz_1<<" "<<trgObj_pt_1<<" "<<pt_1[g]<<" "<<eta_1[g]<<" "<<endl;
		    }

		  }
		  hpbpb_Jet55[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
		  hpbpb_Jet55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
		  if((jet65_1==0) && (jet80_1==0)){ // this is to just check
		    hpbpb_Jet55_noJet65_80[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
		    hpbpb_Jet55_noJet65_80[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
		  }
                }else if(jet65_1==1) {
		  if(trgObj_pt_1>=65 && trgObj_pt_1<80){
		    hpbpb_TrgObj65[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
		    hpbpb_TrgObj65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
		  }
		  hpbpb_Jet65[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
		  hpbpb_Jet65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
		  if(jet80_1==0){
		    hpbpb_Jet65_noJet80[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
		    hpbpb_Jet65_noJet80[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
		  }
                }else if(jet80_1==1) {
		  if(trgObj_pt_1>=80){
		    hpbpb_TrgObj80[k][j][centBin]->Fill(pt_1[g],jet80_p_1);
		    hpbpb_TrgObj80[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1);
		  }
		  hpbpb_Jet80[k][j][centBin]->Fill(pt_1[g],jet80_p_1);
		  hpbpb_Jet80[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1);
		}

		// if(jet80_1) {
		//   hpbpb_Jet80[k][j][centBin]->Fill(pt_1[g],jet80_p_1);
		//   hpbpb_Jet80[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1);
		// }//Jet80 trigger selection

		// if(jet65_1) {
		//   hpbpb_Jet65[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
		//   hpbpb_Jet65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
		// }//Jet65 trigger selection

		// if(jet55_1) {
		//   hpbpb_Jet55[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
		//   hpbpb_Jet55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
		// }//Jet55 trigger selection

		// if(jet65_1 && !jet80_1){
		//   hpbpb_Jet65_noJet80[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
		//   hpbpb_Jet65_noJet80[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
		// }//jet65 and not Jet80 selection 

		// if(jet55_1 && !jet65_1 && !jet80_1){
		//   hpbpb_Jet55_noJet65_80[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
		//   hpbpb_Jet55_noJet65_80[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
		// }//jet55 and not Jet65 and not Jet80 selection 

              }//eta selection

            }//qa cut selection

          }//jet loop

        }//event selection cuts

      }//eta bin loop
    
    }//nentries_jet55or65 loop
    
    
    // //loop for jetpbpb2[2] tree
  //   Long64_t nentries_jet80or95 = jetpbpb2[2][k]->GetEntries();
  //   if(printDebug)cout<<"nentries_jet80or95 = "<<nentries_jet80or95<<endl;
  //   if(printDebug)nentries_jet80or95 = 2;

  //   for(int jentry = 0;jentry<nentries_jet80or95;jentry++){
    
  //     jetpbpb2[2][k]->GetEntry(jentry);
  //     if(printDebug && jentry%100000==0)cout<<"Jet 80or95 file"<<endl;
  //     if(printDebug && jentry%100000==0)cout<<jentry<<": event = "<<evt_2<<", run = "<<run_2<<endl;

  //     //
  //     int centBin = findBin(hiBin_2);//tells us the centrality of the event
  //     if(printDebug)cout<<"centrality bin = "<<5*boundaries_cent[centBin]<< " to "<<5*boundaries_cent[centBin+1]<<endl;

  //     for(int j = 0;j<nbins_eta;j++){

  //       if(pHBHENoiseFilter_2 && pcollisionEventSelection_2 && (vz_2<15)) {
    
  //         for(int g = 0;g<nrefe_2;g++){

  //           if(/*put your favourite QA cut here*/chMax_2[g]/pt_2[g]>0.01){

  //             if(eta_2[g]>=boundaries_eta[j][0] && eta_2[g]<boundaries_eta[j][1]){

  //               if(jet80_2 && trgObj_pt_2>=80){

  //                 hpbpb_TrgObj80[k][j][centBin]->Fill(pt_2[g],jet80_p_2);
  //                 hpbpb_TrgObj80[k][j][nbins_cent]->Fill(pt_2[g],jet80_p_2);

  //               }//80 trg obj selection

  //             }//eta selection

  //           }//qa cut selection

  //         }//jet loop

  //       }//event selection cuts

  //     }//eta bin loop  
    
  //   }//nentries_jet80or95 loop

    fHLT_high[k].close();
    
  }//radius loop. 
  

  //declare the output file
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_data_ak%s%s_%d_endfile_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  
  f.cd();

  for(int k = 0;k<no_radius;k++){

    for(int j = 0;j<nbins_eta;j++){

      for(int i = 0;i<=nbins_cent;i++){

	//divide by delta eta, delta pt, luminosity seen by the triggers and ncoll ... and that will give us cross section/ncoll. 
	// lumi seen by the respective triggers. - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HIJetRAA#Luminosity_Cross_check_for_HIJet
	// once we include the prescl then we have to normalize w.r.t the hightest lumi trigger which would have a prescl of 1 - in this case its the jet80 trigger. 

	//not dividing by ncoll here. will do that in the plotting macro under RAA. 
	
	hpbpb_TrgObj80[k][j][i]->Scale(1./(delta_eta[j]*145.156*1e6));
	hpbpb_TrgObj65[k][j][i]->Scale(1./(delta_eta[j]*145.156*1e6));
	hpbpb_TrgObj55[k][j][i]->Scale(1./(delta_eta[j]*145.156*1e6));
	
	// divide by bin width is doing something very weird. Its making histograms which are empty get 40 entries from somewhere. (TH1.Print Name  = hpbpb_TrgObj80_R3_n20_eta_p20_cent2, Entries= 40, Total sum= 0 example)

	divideBinWidth(hpbpb_TrgObj80[k][j][i]);
	divideBinWidth(hpbpb_TrgObj65[k][j][i]);
	divideBinWidth(hpbpb_TrgObj55[k][j][i]);
	
        hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj80[k][j][i]);
        hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj65[k][j][i]);
        hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj55[k][j][i]);

	// old method histograms, just as a test: 
	
	hpbpb_Jet80[k][j][i]->Scale(1./delta_eta[j]);
	hpbpb_Jet65[k][j][i]->Scale(1./delta_eta[j]);
	hpbpb_Jet55[k][j][i]->Scale(1./delta_eta[j]);

	hpbpb_Jet65_noJet80[k][j][i]->Scale(1./delta_eta[j]);
	hpbpb_Jet55_noJet65_80[k][j][i]->Scale(1./delta_eta[j]);
	
	divideBinWidth(hpbpb_Jet80[k][j][i]);
	divideBinWidth(hpbpb_Jet65[k][j][i]);
	divideBinWidth(hpbpb_Jet55[k][j][i]);

	divideBinWidth(hpbpb_Jet65_noJet80[k][j][i]);
	divideBinWidth(hpbpb_Jet55_noJet65_80[k][j][i]);
	
	hpbpb_JetComb_old[k][j][i]->Add(hpbpb_Jet80[k][j][i]);
	hpbpb_JetComb_old[k][j][i]->Add(hpbpb_Jet65_noJet80[k][j][i]);
	hpbpb_JetComb_old[k][j][i]->Add(hpbpb_Jet55_noJet65_80[k][j][i]);
	
	
        hpbpb_TrgObjComb[k][j][i]->Write();
        if(printDebug)hpbpb_TrgObjComb[k][j][i]->Print();
        hpbpb_TrgObj80[k][j][i]->Write();
        if(printDebug)hpbpb_TrgObj80[k][j][i]->Print();
        hpbpb_TrgObj65[k][j][i]->Write();
        if(printDebug)hpbpb_TrgObj65[k][j][i]->Print();
        hpbpb_TrgObj55[k][j][i]->Write();
        if(printDebug)hpbpb_TrgObj55[k][j][i]->Print();

	hpbpb_Jet80[k][j][i]->Write();
	if(printDebug)hpbpb_Jet80[k][j][i]->Print();
	hpbpb_Jet65[k][j][i]->Write();
	if(printDebug)hpbpb_Jet65[k][j][i]->Print();
	hpbpb_Jet55[k][j][i]->Write();
	if(printDebug)hpbpb_Jet55[k][j][i]->Print();

	hpbpb_Jet65_noJet80[k][j][i]->Write();
	if(printDebug)hpbpb_Jet65_noJet80[k][j][i]->Print();
	hpbpb_Jet55_noJet65_80[k][j][i]->Write();
	if(printDebug)hpbpb_Jet55_noJet65_80[k][j][i]->Print();
	hpbpb_JetComb_old[k][j][i]->Write();
	if(printDebug)hpbpb_JetComb_old[k][j][i]->Print();

      }//cent bins loop

    }//eta bins loop

  }//radius loop
  
  f.Write();

  f.Close();

  //fHLT_high.write();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;

}
