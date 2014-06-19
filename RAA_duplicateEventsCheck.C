// Raghav Kunnawalkam Elayavalli
// June 19th 2014
// CERN

// 
// Macro to check for duplicate events in the latest Jet55or65 forest production and send the output to a text file with information about those particular events. 
// 

// make this macro also run via condor which will be faster. 

#include <iostream>
#include <stdio.h>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>


class DuplicateEvents {
public:
  DuplicateEvents(std::string infname){
    inf = TFile::Open(infname.c_str());
    t = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  };
  ~DuplicateEvents(){
    delete inf;
  }
  
  void MakeList(){
    cout<<"starting make list to check for duplicate events"<<endl;
    evts.clear();
    occurence.clear();
    int run,evt;
    t->SetBranchAddress("run",&run);
    t->SetBranchAddress("evt",&evt);
    for(int i = 0;i<t->GetEntries();i++){
      t->GetEntry(i);
      if(i%100000==0) cout<<i<<" / "<<t->GetEntries()<<" run: "<<run<<" evt: "<<evt<<endl;
      int occur = (int)FindOccurences(run,evt);
      if(occur==0) occurence.push_back(1);
      else occurence.push_back(2);
      evts.push_back(std::make_pair(run,evt));
    }
  }
  int FindOccurences(int run, int evt){
    int noccur = count(evts.begin(),evts.end(),std::make_pair(run,evt));
    return noccur;
  }
  TFile* inf;
  TTree* t;
  vector <pair<int,int> > evts;
  vector <int> occurence;
};


using namespace std;

void RAA_duplicateEventsCheck( int startfile = 0, int endfile = 1, int radius = 3, char *algo = "Vs"){
  
  //TH1::SetDefaultSumw2();
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
  
  for(int ifile = startfile;ifile<endfile;ifile++){
    instr>>filename;
    cout<<"File: "<<filename<<endl;  
  }  
  
  DuplicateEvents dupEvt(filename.c_str());
  dupEvt.MakeList();

 

  //i dont need any of the following stuff since the DuplicateEvents class will do all that for me.  
  // on second thought looks like i need that. :) 
  
  TFile *fin = TFile::Open(filename.c_str());
  
  TTree *jetpbpb1 = (TTree*)fin->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *evtpbpb1 = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  TTree *hltpbpb1 = (TTree*)fin->Get("hltanalysis/HltTree");
  TTree *skmpbpb1 = (TTree*)fin->Get("skimanalysis/HltTree");
  TTree *trgpbpb1 = (TTree*)fin->Get("hltobject/jetObjTree");

  jetpbpb1->AddFriend(evtpbpb1);
  jetpbpb1->AddFriend(hltpbpb1);
  jetpbpb1->AddFriend(skmpbpb1);
  jetpbpb1->AddFriend(trgpbpb1);

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


  //set the branch addresses:  - one of the most boring parts of the code: 
  jetpbpb1->SetBranchAddress("evt",&evt_1);
  jetpbpb1->SetBranchAddress("run",&run_1);
  jetpbpb1->SetBranchAddress("lumi",&lumi_1);
  jetpbpb1->SetBranchAddress("hiBin",&hiBin_1);
  jetpbpb1->SetBranchAddress("vz",&vz_1);
  jetpbpb1->SetBranchAddress("vx",&vx_1);
  jetpbpb1->SetBranchAddress("vy",&vy_1);
  jetpbpb1->SetBranchAddress("hiNtracks",&hiNtracks_1);
  jetpbpb1->SetBranchAddress("hiHFminus",&hiHFminus_1);
  jetpbpb1->SetBranchAddress("hiHFplus",&hiHFplus_1);
  jetpbpb1->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_1);
  jetpbpb1->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_1);
  jetpbpb1->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_1);
  jetpbpb1->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_1);
  //jetpbpb1->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_1);
  //jetpbpb1->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_1);
  
  jetpbpb1->SetBranchAddress("nref",&nrefe_1);
  jetpbpb1->SetBranchAddress("jtpt",&pt_1);
  jetpbpb1->SetBranchAddress("jteta",&eta_1);
  jetpbpb1->SetBranchAddress("jtphi",&phi_1);
  jetpbpb1->SetBranchAddress("rawpt",&raw_1);
  jetpbpb1->SetBranchAddress("chargedMax",&chMax_1);
  jetpbpb1->SetBranchAddress("chargedSum",&chSum_1);
  jetpbpb1->SetBranchAddress("trackMax",&trkMax_1);
  jetpbpb1->SetBranchAddress("trackSum",&trkSum_1);
  jetpbpb1->SetBranchAddress("photonMax",&phMax_1);
  jetpbpb1->SetBranchAddress("photonSum",&phSum_1);
  jetpbpb1->SetBranchAddress("neutralMax",&neMax_1);
  jetpbpb1->SetBranchAddress("neutralSum",&neSum_1);

  //jetpbpb1->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
  //jetpbpb1->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
  //jetpbpb1->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
  //jetpbpb1->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet55_v1",&jet55_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet65_v1",&jet65_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet80_v1",&jet80_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_1);

  jetpbpb1->SetBranchAddress("id",&trgObj_id_1);
  jetpbpb1->SetBranchAddress("pt",&trgObj_pt_1);
  jetpbpb1->SetBranchAddress("eta",&trgObj_eta_1);
  jetpbpb1->SetBranchAddress("phi",&trgObj_phi_1);
  jetpbpb1->SetBranchAddress("mass",&trgObj_mass_1);


  // write the duplicate events to file. 
  
  ofstream outfile;
  outfile.open(Form("pbpb_jet55or65_duplicate_events_loop_run_lumi_event_%d.txt",endfile));
  
  Long64_t nentries = jetpbpb1->GetEntries();

  for(int ievt = 0;ievt<nentries;ievt++){

    jetpbpb1->GetEntry(ievt);

    if(dupEvt.occurence[ievt] >= 2) {

      outfile<<ievt<<" "<<run_1<<" "<<lumi_1<<" "<<evt_1<<" "<<endl;

    }
    
  }

  outfile.close();


}
