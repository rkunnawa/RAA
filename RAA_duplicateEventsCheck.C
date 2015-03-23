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
#include <TNtuple.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include <vector>

class DuplicateEvents {
  
public:
  DuplicateEvents(TTree * t){
    event = (TTree*)t;
  };
  ~DuplicateEvents(){
    delete event;
  }
  
  void MakeList(){
    cout<<"starting make list to check for duplicate events"<<endl;
    evts.clear();
    occurence.clear();
    int run,evt;
    //float vz;
    event->SetBranchAddress("run",&run);
    event->SetBranchAddress("evt",&evt);
    //t->SetBranchAddress("vz",&vz);
    for(int i = 0;i<event->GetEntries();i++){
      event->GetEntry(i);
      if(i%100000==0) cout<<i<<" / "<<event->GetEntries()<<" run: "<<run<<" evt: "<<evt<<endl;
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
  TTree * event;
  vector <pair<int,int> > evts;
  vector <int> occurence;
};

static const int no_radius = 1;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {3};

using namespace std;

void RAA_duplicateEventsCheck( int startfile = 0, int endfile = 902, int radius = 3, char *algo = "Vs"){
  
  //TH1::SetDefaultSumw2();
  TStopwatch timer;
  timer.Start();
  
  // Change to macro to run on condor since its taking a freaking long time. 
  // done! 

  bool printDebug = false;
  
  std::string infile;
  //infile = "jetRAA_PbPb_filelist.txt";
  infile = "jetRAA_MinBiasUPC_forest.txt";
  
  std::ifstream instr(infile.c_str(),std::ifstream::in);
  std::string filename;
  
  //just to read the files till the start number
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr>>filename;
  }

  //filename = "/mnt/hadoop/cms/store/user/dgulhan/HIMC/MB/Track8_Jet26_STARTHI53_LV1/merged2/HiForest_HYDJET_Track8_Jet26_STARTHI53_LV1_merged_forest_0.root";
  
 

  TChain *jetpbpb1;
  string dir = "hiEvtAnalyzer";
  string trees = "HiTree";
  
  //this loop is to assign the tree values before we go into the file loop. 
  jetpbpb1 = new TChain(string(dir+"/"+trees).data());
  
  for(int ifile = startfile;ifile<endfile;ifile++){
    
    instr>>filename;
    if(printDebug)cout<<"File: "<<filename<<endl;
    
    jetpbpb1->Add(filename.c_str());
    if(printDebug)cout << "Tree loaded  " << string(dir+"/"+trees).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb1->GetEntries() << endl;
    
  }// file loop ends

  DuplicateEvents dupEvt(jetpbpb1);
  dupEvt.MakeList();
  
  // event tree
  int evt_1;
  int run_1;
  int lumi_1;

  //set the branch addresses:  - one of the most boring parts of the code: 
  jetpbpb1->SetBranchAddress("evt",&evt_1);
  jetpbpb1->SetBranchAddress("run",&run_1);
  jetpbpb1->SetBranchAddress("lumi",&lumi_1);

  // write the duplicate events to file. 
  
  ofstream outfile;
  outfile.open(Form("/export/d00/scratch/rkunnawa/pbpb_new_minbiasUPF_duplicate_events_run_lumi_event_%d.txt",endfile));
  
  Long64_t nentries = jetpbpb1->GetEntries();

  for(int ievt = 0;ievt<nentries;ievt++){

    jetpbpb1->GetEntry(ievt);

    if(dupEvt.occurence[ievt] >= 2) {
      
      outfile<<ievt<<" "<<run_1<<" "<<lumi_1<<" "<<evt_1<<" "<<endl;
      
    }
    
  }
  
  outfile.close();
  
  
}
