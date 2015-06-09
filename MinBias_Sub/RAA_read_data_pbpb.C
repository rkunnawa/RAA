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

// Oct 27th - still have to run the lumi calculator for the triggers using the lumi mask files. Since we lost some files during the processing (~1%) - done

// Nov 4th - implementing supernova event rejection https://twiki.cern.ch/twiki/pub/CMS/HighPt2014/141102-yenjie-highPt-jetSelection.pdf 

// Dec 7th - added in histograms to show the contamination of the Vs Algorithm to the Jet spectra. 

// Dec 9th - Jet RAA is going to Pu Jets! make the macro to load in the jets and make the ntuple to study the jet ID cuts. 
//           Only going to make the ntuple first. and then we can look at the histograms. 

// Jan 28th forget the ntuples, taking forever. am going directly to the spectra. 
// Feb 2nd - make the spectra from Doga's jet80 spectra. 

// Feb 7th making the change to a TTree from TNtuple at the request of Alex ($!@#$%^&*^%$%^)

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
#include <TNtuple.h>
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

#define NOBJECT_MAX 16384

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

/*
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
*/


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

static const int trigValue = 4;
static const char trigName [trigValue][256] = {"HLT55","HLT65","HLT80","Combined"};
static const Float_t effecPrescl = 2.047507;
static const int nbins_cent = 6;
static Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 5 to get your actual centrality
static Double_t ncoll[nbins_cent+1] = { 1660, 1310, 745, 251, 62.8, 10.8 ,362.24}; //last one is for 0-200 bin. 

static const int no_radius = 1;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {3};

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

int findBin(int bin){

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


using namespace std;

void RAA_read_data_pbpb(int startfile = 0, int endfile = 22, char *algo = "Pu", char *jet_type = "PF"){

  TH1::SetDefaultSumw2();
  //gStyle->SetOptStat(0);
  
  TStopwatch timer;
  timer.Start();
  
  TDatime date;
  
  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = false;
  int istightCut = 0;

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
  //infile1 = "doga_jet80_pbpb_data.txt";
  
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

  //jetpbpb2 and dir2 are all for akPu3PF  - removed since we are only looking for akPuPF jets 

  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%d%sJetAnalyzer",algo,list_radius[k],jet_type);
    //dir[3][k] = Form("akPu%d%sJetAnalyzer",list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    dir[4][k] = "hltobject";
    //dir[6][k] = "pfcandAnalyzer";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    //"t",
    "HiTree",
    "jetObjTree",
    //"pfTree"
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
      
    }// radius loop ends
    
  }// file loop ends
  
  


  for(int k = 0;k<no_radius;k++){
    jetpbpb1[2][k]->AddFriend(jetpbpb1[0][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[1][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[3][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[4][k]);

  }// radius loop ends
 
  // Ok this should now work. 

  cout<<"total no of entries in the combined forest files = "<<jetpbpb1[2][0]->GetEntries()<<endl;
  //  if(printDebug)cout<<"total no of entries in the Jet80 Tree     = "<<jetpbpb2[2][0]->GetEntries()<<endl;


  
  //file 1: 
  // jet tree 1
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
  float eMax_1[1000];
  float muMax_1[1000];
  float eSum_1[1000];
  float muSum_1[1000];
  float jtpu_1[1000];
  /*
  // jet tree 2
  int nrefe_2;
  float pt_2[1000];
  //float old_2pt3[1000];
  float raw_2[1000];
  float eta_2[1000];
  float eta_2_CM[1000];
  float phi_2[1000];
  float chMax_2[1000];
  float trkMax_2[1000];
  float chSum_2[1000];
  float phSum_2[1000];
  float neSum_2[1000];
  float trkSum_2[1000];
  float phMax_2[1000];
  float neMax_2[1000];
  float eMax_2[1000];
  float muMax_2[1000];
  float eSum_2[1000];
  float muSum_2[1000];
  float jtpu_2[1000];
  */
  // event tree
  int evt_1;
  int run_1;
  int lumi_1;
  int hiBin_1;
  float vx_1;
  float vy_1;
  float vz_1;
  int hiNpix_1;
  int hiNtracks_1;
  float hiHF_1;
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
  int L1_sj36_1;
  int L1_sj52_1;
  int jetMB_1;
  int jet55_1;
  int jet65_1;
  int jet80_1;
  int jetMB_p_1;
  int L1_sj36_p_1;
  int L1_sj52_p_1;
  int jet55_p_1;
  int jet65_p_1;
  int jet80_p_1;


  // trigger object tree - this contains the maximum value of the particular trigger object. 
  float trgObj_id_1;
  float trgObj_pt_1;
  float trgObj_eta_1;
  float trgObj_phi_1;
  float trgObj_mass_1;
  
  /*
  Int_t nPFpart;
  Int_t pfId[NOBJECT_MAX];
  Float_t pfPt[NOBJECT_MAX];
  Float_t pfVsPtInitial[NOBJECT_MAX];
  Float_t pfVsPt[NOBJECT_MAX];
  Float_t pfEta[NOBJECT_MAX];
  Float_t pfPhi[NOBJECT_MAX];
  Float_t pfArea[NOBJECT_MAX];
  Float_t v_n[5][15];
  Float_t psi_n[5][15];
  Float_t sumpT[15];
  */
  
  for(int k = 0;k<no_radius;k++){
    //set the branch addresses:  - one of the most boring parts of the code: 
    jetpbpb1[2][k]->SetBranchAddress("evt",&evt_1);
    jetpbpb1[2][k]->SetBranchAddress("run",&run_1);
    jetpbpb1[2][k]->SetBranchAddress("lumi",&lumi_1);
    jetpbpb1[2][k]->SetBranchAddress("hiBin",&hiBin_1);
    jetpbpb1[2][k]->SetBranchAddress("vz",&vz_1);
    jetpbpb1[2][k]->SetBranchAddress("vx",&vx_1);
    jetpbpb1[2][k]->SetBranchAddress("vy",&vy_1);
    jetpbpb1[2][k]->SetBranchAddress("hiNpix",&hiNpix_1);
    jetpbpb1[2][k]->SetBranchAddress("hiNtracks",&hiNtracks_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHFminus",&hiHFminus_1);
    jetpbpb1[2][k]->SetBranchAddress("hiHF",&hiHF_1);
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
    jetpbpb1[2][k]->SetBranchAddress("jtpu",&jtpu_1);
    jetpbpb1[2][k]->SetBranchAddress("chargedMax",&chMax_1);
    jetpbpb1[2][k]->SetBranchAddress("chargedSum",&chSum_1);
    jetpbpb1[2][k]->SetBranchAddress("trackMax",&trkMax_1);
    jetpbpb1[2][k]->SetBranchAddress("trackSum",&trkSum_1);
    jetpbpb1[2][k]->SetBranchAddress("photonMax",&phMax_1);
    jetpbpb1[2][k]->SetBranchAddress("photonSum",&phSum_1);
    jetpbpb1[2][k]->SetBranchAddress("neutralMax",&neMax_1);
    jetpbpb1[2][k]->SetBranchAddress("neutralSum",&neSum_1);
    jetpbpb1[2][k]->SetBranchAddress("eSum",&eSum_1);
    jetpbpb1[2][k]->SetBranchAddress("eMax",&eMax_1);
    jetpbpb1[2][k]->SetBranchAddress("muSum",&muSum_1);
    jetpbpb1[2][k]->SetBranchAddress("muMax",&muMax_1);

    /*
    jetpbpb1[3][k]->SetBranchAddress("nref",&nrefe_2);
    jetpbpb1[3][k]->SetBranchAddress("jtpt",&pt_2);
    jetpbpb1[3][k]->SetBranchAddress("jteta",&eta_2);
    jetpbpb1[3][k]->SetBranchAddress("jtphi",&phi_2);
    jetpbpb1[3][k]->SetBranchAddress("rawpt",&raw_2);
    jetpbpb1[3][k]->SetBranchAddress("jtpu",&jtpu_2);
    jetpbpb1[3][k]->SetBranchAddress("chargedMax",&chMax_2);
    jetpbpb1[3][k]->SetBranchAddress("chargedSum",&chSum_2);
    jetpbpb1[3][k]->SetBranchAddress("trackMax",&trkMax_2);
    jetpbpb1[3][k]->SetBranchAddress("trackSum",&trkSum_2);
    jetpbpb1[3][k]->SetBranchAddress("photonMax",&phMax_2);
    jetpbpb1[3][k]->SetBranchAddress("photonSum",&phSum_2);
    jetpbpb1[3][k]->SetBranchAddress("neutralMax",&neMax_2);
    jetpbpb1[3][k]->SetBranchAddress("neutralSum",&neSum_2);
    jetpbpb1[3][k]->SetBranchAddress("eSum",&eSum_2);
    jetpbpb1[3][k]->SetBranchAddress("eMax",&eMax_2);
    jetpbpb1[3][k]->SetBranchAddress("muSum",&muSum_2);
    jetpbpb1[3][k]->SetBranchAddress("muMax",&muMax_2);
    */
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
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_1);

    jetpbpb1[2][k]->SetBranchAddress("id",&trgObj_id_1);
    jetpbpb1[2][k]->SetBranchAddress("pt",&trgObj_pt_1);
    jetpbpb1[2][k]->SetBranchAddress("eta",&trgObj_eta_1);
    jetpbpb1[2][k]->SetBranchAddress("phi",&trgObj_phi_1);
    jetpbpb1[2][k]->SetBranchAddress("mass",&trgObj_mass_1);
    /*
    jetpbpb1[2][k]->SetBranchAddress("nPFpart", &nPFpart);
    jetpbpb1[2][k]->SetBranchAddress("pfId", pfId);
    jetpbpb1[2][k]->SetBranchAddress("pfPt", pfPt);
    jetpbpb1[2][k]->SetBranchAddress("pfVsPtInitial", pfVsPtInitial);
    jetpbpb1[2][k]->SetBranchAddress("pfVsPt", pfVsPt);
    jetpbpb1[2][k]->SetBranchAddress("pfEta", pfEta);
    jetpbpb1[2][k]->SetBranchAddress("pfPhi", pfPhi);
    jetpbpb1[2][k]->SetBranchAddress("pfArea", pfArea);
    jetpbpb1[2][k]->SetBranchAddress("vn",&v_n);
    jetpbpb1[2][k]->SetBranchAddress("psin",&psi_n);
    jetpbpb1[2][k]->SetBranchAddress("sumpt",&sumpT);
    */
  }//radius loop

  //now that we have all the branch addresses set up we can start going through the loop to look at the trigger objects 
  
  //before we go and do the trigger object merging lets do the text file information here: 
  //we need one text file per spectra, which means we need 3 -> one for each trigger object 
  //and we also need one for the high pt jets with low trigger objects 
  // and these files will have the following structure: 
  // run lumi evt HLTobjpt HLTobjeta HLTobjphi hibin jtpt jteta jtphi

  //ofstream fHLT_80,fHLT_65,fHLT_55;

  
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_data_vz_cent_ak%s%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  f.cd();
  // need to keep all the input histograms in scope of the output file. 

  TH1F *hEvents_HLT80 = new TH1F("hEvents_HLT80","",4,0,2);
  TH1F *hEvents_HLT65 = new TH1F("hEvents_HLT65","",4,0,2);
  TH1F *hEvents_HLT55 = new TH1F("hEvents_HLT55","",4,0,2);
  TH1F *hEvents = new TH1F("hEvents","",4,0,2);
  TH1F *hCentEvents = new TH1F("hCentEvents","",10,0,10);

#if 0
  TH1F *hEvents_pCES = new TH1F("hEvents_pCES","",4,0,2);
  TH1F *hEvents_pHBHE = new TH1F("hEvents_pHBHE","",4,0,2);
  TH1F *hEvents_vz15 = new TH1F("hEvents_vz15","",4,0,2);
  TH1F *hEvents_eta2 = new TH1F("hEvents_eta2","",4,0,2);
  TH1F *hEvents_supernova = new TH1F("hEvents_supernova","",4,0,2);
  TH1F *hEvents_chMaxjtpt = new TH1F("hEvents_chMaxjtpt","",4,0,2);
  TH1F *hpbpb_TrgObj80[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_TrgObj65[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_TrgObj55[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_TrgObjComb[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_RecoOverRaw[no_radius][nbins_eta][nbins_cent+1];
  TH2F *hpbpb_RecoOverRaw_jtpt[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_chMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eMax[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_chSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eSum[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_chMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eMax_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_chSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_phSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_neSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_muSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  TH1F *hpbpb_eSum_withCut[no_radius][nbins_eta][nbins_cent+1];
  for(int k = 0;k<no_radius;k++){
    for(int j = 0;j<nbins_eta;j++){
      for(int i = 0;i<nbins_cent;i++){
        hpbpb_TrgObj80[k][j][i] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
        hpbpb_TrgObj65[k][j][i] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
        hpbpb_TrgObj55[k][j][i] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_TrgObjComb[k][j][i] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	hpbpb_RecoOverRaw[k][j][i] = new TH1F(Form("hpbpb_RecoOverRaw_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco over raw ratio R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,10);
	hpbpb_RecoOverRaw_jtpt[k][j][i] = new TH2F(Form("hpbpb_RecoOverRaw_jtpt_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("reco over raw ratio versus jtpt R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,100,0,10);
	hpbpb_chMax[k][j][i] = new TH1F(Form("hpbpb_chMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phMax[k][j][i] = new TH1F(Form("hpbpb_phMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neMax[k][j][i] = new TH1F(Form("hpbpb_neMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muMax[k][j][i] = new TH1F(Form("hpbpb_muMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eMax[k][j][i] = new TH1F(Form("hpbpb_eMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eMax variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_chSum[k][j][i] = new TH1F(Form("hpbpb_chSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phSum[k][j][i] = new TH1F(Form("hpbpb_phSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neSum[k][j][i] = new TH1F(Form("hpbpb_neSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muSum[k][j][i] = new TH1F(Form("hpbpb_muSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eSum[k][j][i] = new TH1F(Form("hpbpb_eSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eSum variable for R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_chMax_withCut[k][j][i] = new TH1F(Form("hpbpb_chMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phMax_withCut[k][j][i] = new TH1F(Form("hpbpb_phMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neMax_withCut[k][j][i] = new TH1F(Form("hpbpb_neMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muMax_withCut[k][j][i] = new TH1F(Form("hpbpb_muMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eMax_withCut[k][j][i] = new TH1F(Form("hpbpb_eMax_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eMax variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_chSum_withCut[k][j][i] = new TH1F(Form("hpbpb_chSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("chSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_phSum_withCut[k][j][i] = new TH1F(Form("hpbpb_phSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("phSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_neSum_withCut[k][j][i] = new TH1F(Form("hpbpb_neSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("neSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_muSum_withCut[k][j][i] = new TH1F(Form("hpbpb_muSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("muSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);
	hpbpb_eSum_withCut[k][j][i] = new TH1F(Form("hpbpb_eSum_withCut_R%d_%s_cent%d",list_radius[k],etaWidth[j],i),Form("eSum variable for withCut_R%d %s %2.0f - %2.0f cent",list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200);	
      }//cent bin loop
      hpbpb_TrgObj80[k][j][nbins_cent] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra from Jet 80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_TrgObj65[k][j][nbins_cent] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra from Jet 65 && !jet80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_TrgObj55[k][j][nbins_cent] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_TrgObjComb[k][j][nbins_cent] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("Trig combined spectra R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000);
      hpbpb_RecoOverRaw[k][j][nbins_cent] = new TH1F(Form("hpbpb_RecoOverRaw_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco over raw ratio R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,10);
      hpbpb_RecoOverRaw_jtpt[k][j][nbins_cent] = new TH2F(Form("hpbpb_RecoOverRaw_jtpt_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("reco over raw ratio versus jtpt R%d %s 0-200 cent",list_radius[k],etaWidth[j]),1000,0,1000,100,0,10);
      
      hpbpb_chMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_chMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("chMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_phMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_phMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("phMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_neMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_neMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("neMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_muMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_muMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("muMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_eMax[k][j][nbins_cent] = new TH1F(Form("hpbpb_eMax_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("eMax variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      
      hpbpb_chSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_chSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("chSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_phSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_phSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("phSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_neSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_neSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("neSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_muSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_muSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("muSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
      hpbpb_eSum[k][j][nbins_cent] = new TH1F(Form("hpbpb_eSum_R%d_%s_cent%d",list_radius[k],etaWidth[j],nbins_cent),Form("eSum variable for R%d %s 0-200 cent",list_radius[k],etaWidth[j]),100,0,200);
    }
   }
  TH1F * hpbpb_Jet[trigValue-1][nbins_cent];
  TH1F * hpbpb_Jet_chMaxJtpt[trigValue-1][5][nbins_cent];
  TH1F * hpbpb_Jet_eMaxJtpt[trigValue-1][9][nbins_cent];
  //TH1F * hpbpb_Jet_eMaxSumcand[trigValue-1][9][nbins_cent];
  TH1F * hpbpb_Jet_eMaxJtpt_chMaxJtpt[trigValue-1][3][3][nbins_cent];
  for(int i = 0;i<nbins_cent;++i){
    
    for(int t = 0;t<trigValue-1;t++){
      hpbpb_Jet[t][i] = new TH1F(Form("hpbpb_%s_cent%d",trigName[t], i),"",100,0,300);
      for(int a = 0;a<5;++a)
	hpbpb_Jet_chMaxJtpt[t][a][i] = new TH1F(Form("hpbpb_%s_chMaxJtpt0p0%d_cent%d",trigName[t], a+1, i),"",100,0,300);
      for(int a = 0;a<9;++a)
	hpbpb_Jet_eMaxJtpt[t][a][i] = new TH1F(Form("hpbpb_%s_eMaxJtpt0p%d_cent%d", trigName[t], a+1, i),"",100,0,300);
      
      for(int a = 0;a<3;++a)
	for(int b = 0;b<3;++b)  
	  hpbpb_Jet_eMaxJtpt_chMaxJtpt[t][a][b][i] = new TH1F(Form("hpbpb_%s_eMaxJtpt0p0%d_chMaxJtpt0p%d_cent%d",trigName[t], a+5, b+1, i),"",100,0,300);
    }
  }
  
  static const int ptSelection = 19;
  static const int ptBoundary[ptSelection+1] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
  TH2F* hpbpb_chMaxJtpt_jtpt = new TH2F("hpbpb_chMaxJtpt_jtpt","",100,0,300,500,0,5);
  TH2F* hpbpb_eMaxJtpt_jtpt = new TH2F("hpbpb_eMaxJtpt_jtpt","",100,0,300,500,0,5);
  TH2F* hpbpb_eMaxSumcand_jtpt = new TH2F("hpbpb_eMaxSumcand_jtpt","",100,0,300,1000,0,10);
  TH2F* hpbpb_eMaxSumcand_chMaxJtpt = new TH2F("hpbpb_eMaxSumcand_chMaxJtpt","",500,0,5,1000,0,10);
  TH2F* hpbpb_eMaxJtpt_chMaxJtpt = new TH2F("hpbpb_eMaxJtpt_chMaxJtpt","",500,0,5,500,0,5);
  TH2F* hpbpb_chMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_phMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_neMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_muMaxJtpt_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_jtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxSumcand_eMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_eMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_neMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_phMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];
  TH2F* hpbpb_muMaxJtpt_chMaxJtpt_ptselection[trigValue][ptSelection];
  for(int a = 0;a<ptSelection;a++){
    for(int t = 0;t<trigValue;++t){
      hpbpb_chMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_chMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_eMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_phMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_phMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_neMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_neMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_muMaxJtpt_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_muMaxJtpt_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,500,0,5);
      hpbpb_eMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_phMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_phMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_neMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_neMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_muMaxJtpt_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_muMaxJtpt_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,500,0,5);
      hpbpb_eMaxSumcand_chMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_chMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,1000,0,10);
      hpbpb_eMaxSumcand_eMaxJtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_eMaxJtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",500,0,5,1000,0,10);
      hpbpb_eMaxSumcand_jtpt_ptselection[t][a] = new TH2F(Form("hpbpb_eMaxSumcand_jtpt_%d_jtpt_%d_%s",ptBoundary[a],ptBoundary[a+1],trigName[t]),"",100,0,300,1000,0,10);
    }
  }
  TH1F *hpbpb_vx[no_radius];
  TH1F *hpbpb_vy[no_radius];
  TH1F *hpbpb_vz[no_radius];
  TH1F *hpbpb_cent[no_radius];
#endif
  
  for(int k = 0;k<no_radius;k++){

    // hpbpb_cent[k] = new TH1F(Form("hpbpb_cent_R%d",list_radius[k]),Form("centrality distributions R%d",list_radius[k]),200,0,200);
    // hpbpb_vz[k] = new TH1F(Form("hpbpb_vz_R%d",list_radius[k]),Form("vz distribution R%d",list_radius[k]),60,-15,15);
    // hpbpb_vx[k] = new TH1F(Form("hpbpb_vx_R%d",list_radius[k]),Form("vx distribution R%d",list_radius[k]),60,-15,15);
    // hpbpb_vy[k] = new TH1F(Form("hpbpb_vy_R%d",list_radius[k]),Form("vy distribution R%d",list_radius[k]),60,-15,15);

    if(printDebug)cout<<"Running data reading for R = "<<list_radius[k]<<endl;
    // loop for the jetpbpb1[2] tree 
    Long64_t nentries_jet55or65 = jetpbpb1[2][k]->GetEntries();
    if(printDebug)cout<<"nentries_jet55or65or80 = "<<nentries_jet55or65<<endl;
    if(printDebug)nentries_jet55or65 = 10;

    // fill the no of pixel vs no of jets histogram here 
    //jetpbpb1[2][k]->Draw(Form("hiNpix:Sum$(jtpt>50&&abs(jteta)<2)>>hpbpb_Npix_cut_R%d_n20_eta_p20_cent%d",list_radius[k],centBin),"pcollisionEventSelection&&pHBHENoiseFilter","goff");
    //jetpbpb1[2][k]->Project(Form("hpbpb_Npix_cut_R%d_n20_eta_p20_cent%d",list_radius[k],nbins_cent),"hiNpix:Sum$(jtpt>50&&abs(jteta)<2)","pcollisionEventSelection&&pHBHENoiseFilter");

    for(int jentry = 0;jentry<nentries_jet55or65;jentry++){
      jetpbpb1[0][k]->GetEntry(jentry);
      jetpbpb1[1][k]->GetEntry(jentry);
      jetpbpb1[2][k]->GetEntry(jentry);    
      jetpbpb1[3][k]->GetEntry(jentry);
      jetpbpb1[4][k]->GetEntry(jentry);

      //if(printDebug && jentry%100000==0)cout<<"Jet 55or65 file"<<endl;
      if(printDebug)cout<<jentry<<": event = "<<evt_1<<"; run = "<<run_1<<"; lumi = "<<lumi_1<<endl;
      
      // get the stuff required for the trigger turn on curve later. in a separate loop till i understand how to put this in here. 
      
      int centBin = findBin(hiBin_1);//tells us the centrality of the event. 
      if(centBin==-1 || centBin==nbins_cent) continue;

      //if(printDebug)cout<<"cent bin = "<<centBin<<endl;
      //if(printDebug)cout<<"centrality bin = "<<5*boundaries_cent[centBin]<< " to "<<5*boundaries_cent[centBin+1]<<endl;
      // if(k==1)hEvents->Fill(1);
      // if(k==1 && pcollisionEventSelection_1 == 1) hEvents_pCES->Fill(1);
      // if(k==1 && pcollisionEventSelection_1 == 1 && pHBHENoiseFilter_1 ==1) hEvents_pHBHE->Fill(1);
      // if(k==1 && pcollisionEventSelection_1 == 1 && pHBHENoiseFilter_1 ==1 && fabs(vz_1)<15) hEvents_vz15->Fill(1);
      
      if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 

      if(TMath::Abs(vz_1)>15) continue;
      
      // if(printDebug)cout<<" trigger object pt =  "<<trgObj_pt_1<<endl;
      // if(printDebug)cout<<" jet80 = "<<jet80_1<<endl;
      // if(printDebug)cout<<" jet65 = "<<jet65_1<<endl;
      // if(printDebug)cout<<" jet55 = "<<jet55_1<<endl;
      // if(printDebug)cout<<" nrefe = "<<nrefe_1<<endl;
      
      int jetCounter = 0;//counts jets which are going to be used in the supernova cut rejection. 
      
      //if(algo=="Vs"){
     
      for(int j = 0;j<nbins_eta;j++){
	  
	for(int g = 0;g<nrefe_1;g++){
	    
	  if(eta_1[g]>=boundaries_eta[j][0] && eta_1[g]<boundaries_eta[j][1]){
	      
	    if(pt_1[g]>=50) jetCounter++;
	    
	  }//eta selection cut
	  
	}// jet loop
	
      }//eta bins loop
      

      // if(printDebug)cout<<"pixel hit = "<<hiNpix_1<<", jet counter = "<<jetCounter<<endl;
      
      // hpbpb_Npix_before_cut[k][centBin]->Fill(jetCounter,hiNpix_1);
      // hpbpb_Npix_before_cut[k][nbins_cent]->Fill(jetCounter,hiNpix_1);	
      
      //if(hiBin_1 >= 0 && hiBin_1 < 1) hpbpb_Npix_before_cut[k][nbins_cent+1]->Fill(jetCounter,hiNpix_1);	

      // if(k==1 && pcollisionEventSelection_1 == 1 && pHBHENoiseFilter_1 ==1 && fabs(vz_1)<15 && (hiNpix_1 < 38000 - 500*jetCounter)) hEvents_supernova->Fill(1);
      //if(k==1 && pcollisionEventSelection_1 == 1 && pHBHENoiseFilter_1 ==1 && fabs(vz_1)<15 && (hiNpix_1 < 38000 - 500*jetCounter) && eta_1[0] > -2 && eta_1[0] < 2) hEvents_eta2->Fill(1);
      
      // apply the correct supernova selection cut rejection here: 
      if(hiNpix_1 > 38000 - 500*jetCounter){
       	if(printDebug) cout<<"removed this supernova event"<<endl;
      	continue;
      }

      hEvents->Fill(1);
      hCentEvents->Fill(centBin);
      if(jet80_1)hEvents_HLT80->Fill(1);
      
      // if(chMax_1[0]/pt_1[0] < 0.02 || eMax_1[0]/pt_1[0] > 0.6) continue;

      // hpbpb_cent[k]->Fill(hiBin_1);
      // hpbpb_vz[k]->Fill(vz_1);
      // hpbpb_vx[k]->Fill(vx_1);
      // hpbpb_vy[k]->Fill(vy_1);
#if 0
      //get the failed events after the beam scrapping and pile up rejection cuts. 
      //if(jetCounter>10){
      //put the failure mode event variables here: 
      //fVs_failure[k]<<run_1<<" "<<lumi_1<<" "<<evt_1<<" "<<vz_1<<" "<<hiHF_1<<" "<<hiNpix_1<<" "<<hiNtracks_1<<endl;
      //      }
      
      //hpbpb_Npix_after_cut[k][centBin]->Fill(jetCounter,hiNpix_1);
      //hpbpb_Npix_after_cut[k][nbins_cent]->Fill(jetCounter,hiNpix_1);
      //}//if Vs search for supernova events. 
      
      // hEvents->Fill(1);
      // if(jet80_1) hEvents_HLT80->Fill(1);
      // if(jet65_1 && !jet80_1) hEvents_HLT65->Fill(1,jet65_p_1);
      // if(jet65_1 && !jet80_1) hEvents_HLT65->Fill(0);
      // if(jet55_1 && !jet65_1 && !jet80_1) hEvents_HLT55->Fill(1,jet55_p_1);
      // if(jet55_1 && !jet65_1 && !jet80_1) hEvents_HLT55->Fill(0);
      
      
#if 0
      Float_t Vs_0_x_minus = sumpT[0]*v_n[0][0]*TMath::Cos(0*psi_n[0][0]);
      Float_t Vs_0_x_plus = sumpT[14]*v_n[0][14]*TMath::Cos(0*psi_n[0][14]);
      Float_t Vs_0_y_minus = sumpT[0]*v_n[0][0]*TMath::Sin(0*psi_n[0][0]);
      Float_t Vs_0_y_plus = sumpT[14]*v_n[0][14]*TMath::Sin(0*psi_n[0][14]);
      Float_t Vs_0_x = Vs_0_x_minus + Vs_0_x_plus;
      Float_t Vs_0_y = Vs_0_y_minus + Vs_0_y_plus;
      //if(printDebug)std::cout<<"Vs_0_x = "<<Vs_0_x<<"; Vs_0_y =  "<<Vs_0_y<<std::endl;
      //if(TMath::Abs(Vs_0_x)> v0_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v0 > 5200"<<std::endl;
      //if(TMath::Abs(Vs_0_y)> v0_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v0 > 5200"<<std::endl;
      if(TMath::Abs(Vs_0_x)> v0_tight || TMath::Abs(Vs_0_y)> v0_tight) hEvents_Diverge_v0_tight->Fill(1);
      if(TMath::Abs(Vs_0_x)> v0_loose || TMath::Abs(Vs_0_y)> v0_loose) hEvents_Diverge_v0_loose->Fill(1);
      
      Float_t Vs_1_x_minus = sumpT[0]*v_n[1][0]*TMath::Cos(1*psi_n[1][0]);
      Float_t Vs_1_x_plus = sumpT[14]*v_n[1][14]*TMath::Cos(1*psi_n[1][14]);
      Float_t Vs_1_y_minus = sumpT[0]*v_n[1][0]*TMath::Sin(1*psi_n[1][0]);
      Float_t Vs_1_y_plus = sumpT[14]*v_n[1][14]*TMath::Sin(1*psi_n[1][14]);
      Float_t Vs_1_x = Vs_1_x_minus + Vs_1_x_plus;
      Float_t Vs_1_y = Vs_1_y_minus + Vs_1_y_plus;
      //if(printDebug)std::cout<<"Vs_1_x = "<<Vs_1_x<<"; Vs_1_y =  "<<Vs_1_y<<std::endl;
      //if(TMath::Abs(Vs_1_x)> v1_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v1 > 100"<<std::endl;
      //if(TMath::Abs(Vs_1_y)> v1_tight && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v1 > 100"<<std::endl;
      if(TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight) hEvents_Diverge_v1_tight->Fill(1);
      if(TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose) hEvents_Diverge_v1_loose->Fill(1);
      Float_t Vs_2_x_minus = sumpT[0]*v_n[2][0]*TMath::Cos(2*psi_n[2][0]);
      Float_t Vs_2_x_plus = sumpT[14]*v_n[2][14]*TMath::Cos(2*psi_n[2][14]);
      Float_t Vs_2_y_minus = sumpT[0]*v_n[2][0]*TMath::Sin(2*psi_n[2][0]);
      Float_t Vs_2_y_plus = sumpT[14]*v_n[2][14]*TMath::Sin(2*psi_n[2][14]);
      Float_t Vs_2_x = Vs_2_x_minus + Vs_2_x_plus;
      Float_t Vs_2_y = Vs_2_y_minus + Vs_2_y_plus;
      //if(printDebug)std::cout<<"Vs_2_x = "<<Vs_2_x<<"; Vs_2_y =  "<<Vs_2_y<<std::endl;
      //if(TMath::Abs(Vs_2_x)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v2 > 140"<<std::endl;
      //if(TMath::Abs(Vs_2_y)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v2 > 140"<<std::endl;
      if(TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight) hEvents_Diverge_v2_tight->Fill(1);
      if(TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose) hEvents_Diverge_v2_loose->Fill(1);
      Float_t Vs_3_x_minus = sumpT[0]*v_n[3][0]*TMath::Cos(3*psi_n[3][0]);
      Float_t Vs_3_x_plus = sumpT[14]*v_n[3][14]*TMath::Cos(3*psi_n[3][14]);
      Float_t Vs_3_y_minus = sumpT[0]*v_n[3][0]*TMath::Sin(3*psi_n[3][0]);
      Float_t Vs_3_y_plus = sumpT[14]*v_n[3][14]*TMath::Sin(3*psi_n[3][14]);
      Float_t Vs_3_x = Vs_3_x_minus + Vs_3_x_plus;
      Float_t Vs_3_y = Vs_3_y_minus + Vs_3_y_plus;
      //if(printDebug)std::cout<<"Vs_3_x = "<<Vs_3_x<<"; Vs_3_y =  "<<Vs_3_y<<std::endl;
      //if(TMath::Abs(Vs_3_x)>120 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v3 > 120"<<std::endl;
      //if(TMath::Abs(Vs_3_y)>120 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v3 > 120"<<std::endl;
      if(TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight) hEvents_Diverge_v3_tight->Fill(1);
      if(TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose) hEvents_Diverge_v3_loose->Fill(1);
      Float_t Vs_4_x_minus = sumpT[0]*v_n[4][0]*TMath::Cos(4*psi_n[4][0]);
      Float_t Vs_4_x_plus = sumpT[14]*v_n[4][14]*TMath::Cos(4*psi_n[4][14]);
      Float_t Vs_4_y_minus = sumpT[0]*v_n[4][0]*TMath::Sin(4*psi_n[4][0]);
      Float_t Vs_4_y_plus = sumpT[14]*v_n[4][14]*TMath::Sin(4*psi_n[4][14]);
      Float_t Vs_4_x = Vs_4_x_minus + Vs_4_x_plus;
      Float_t Vs_4_y = Vs_4_y_minus + Vs_4_y_plus;
      //if(printDebug)std::cout<<"Vs_4_x = "<<Vs_4_x<<"; Vs_4_y =  "<<Vs_4_y<<std::endl;
      //if(TMath::Abs(Vs_4_x)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* cos v4 > 140"<<std::endl;
      //if(TMath::Abs(Vs_4_y)>140 && printDebug) std::cout<<"event "<<jentry<<" has sumpT* sin v4 > 140"<<std::endl;
      if(TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) hEvents_Diverge_v4_tight->Fill(1);
      if(TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) hEvents_Diverge_v4_loose->Fill(1);
      if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
	hEvents_Diverge_tight->Fill(1);
	isDivergeTight = true;
      }
      if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
	hEvents_Diverge_loose->Fill(1);
	isDivergeLoose = true;
      }
      
      
      if(jet55_1 && (eMax_1[0]/pt_1[0]>0.8) && pt_1[0]>80) {
	event_value = evt_1;
	run_value = run_1;
	lumi_value = lumi_1;
	evt_electron_failure[k]->Fill();
      }
      if(jet55_1 && (eMax_1[0]/pt_1[0]<0.4) && pt_1[0]>80) {
	event_value = evt_1;
	run_value = run_1;
	lumi_value = lumi_1;
	evt_electron_good[k]->Fill();
      }
#endif
      
      for(int j = 0;j<nbins_eta;j++){
	for(int g = 0;g<nrefe_1;g++){ // this is the loop for the  Jets we are interested in.  
	  
	  //if((chMax_1[g]/pt_1[g]>0.05) && (chMax_1[g]/pt_1[g]<1.00)){
	  
	  if(eta_1[g]<boundaries_eta[j][0] || eta_1[g]>=boundaries_eta[j][1]) continue;
	  
	  //if(raw_1[g] < 20) continue;
	  
	  // if(printDebug)cout<<"jet variables:  "<<endl;
	  // if(printDebug)cout<<"chSum = "<<chSum_1[g]<<endl;
	  // if(printDebug)cout<<"phSum = "<<phSum_1[g]<<endl;
	  // if(printDebug)cout<<"neSum = "<<neSum_1[g]<<endl;
	  // if(printDebug)cout<<"muSum = "<<muSum_1[g]<<endl;
	  // if(printDebug)cout<<"eSum  = "<<eSum_1[g]<<endl;
	  // if(printDebug)cout<<"jtpt  = "<<pt_1[g]<<endl;
	  // if(printDebug)cout<<"rawpT = "<<raw_1[g]<<endl;
	  // if(printDebug)cout<<"jtpu  = "<<jtpu_1[g]<<endl;
	  // cut3 = (float)(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g] + eSum_1[g] - jtpu_1[g])/(pt_1[g]);
	  // cut1 = (float)chMax_1[g]/(pt_1[g]);
	  // cut2 = (float)neMax_1[g]/TMath::Max(chSum_1[g],neSum_1[g]);
	  // cut2b = (float)phMax_1[g]/TMath::Max(chSum_1[g],neSum_1[g]);
	  // cut4 = (float)(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g] + eSum_1[g])/(0.5*raw_1[g]);
	  // cut5 = (float)neMax_1[g]/(neMax_1[g] + chMax_1[g] + phMax_1[g]);
	  // cut6 = (float)muMax_1[g]/(neMax_1[g] + chMax_1[g] + phMax_1[g]);
	  // cut7 = (float)eMax_1[g]/(neMax_1[g] + chMax_1[g] + phMax_1[g]);
	  // if(printDebug)cout<<"Cut1 value = "<<cut1<<endl;
	  // if(printDebug)cout<<"Cut2 value = "<<cut2<<endl;
	  // if(printDebug)cout<<"Cut3 value = "<<cut3<<endl;
	  // if(printDebug)cout<<"Cut4 value = "<<cut4<<endl;
	  // if(printDebug)cout<<"Cut4 value = "<<cut5<<endl;
	  
	  //hCut3->Fill(cut3);
	  //hCut5->Fill(cut5);
	  //hCut4->Fill(cut4);
	  //hCut2->Fill(cut2);
	  //hCut1->Fill(cut1);
	  
#if 0
	  arrayValues[0] = raw_1[g];
	  arrayValues[1] = pt_1[g];
	  arrayValues[2] = jtpu_1[g];
	  arrayValues[3] = L1_sj36_1;
	  arrayValues[4] = L1_sj36_p_1;
	  arrayValues[5] = L1_sj52_1;
	  arrayValues[6] = L1_sj52_p_1;
	  arrayValues[7] = jet55_1;
	  arrayValues[8] = jet55_p_1;
	  arrayValues[9] = jet65_1;
	  arrayValues[10] = jet65_p_1;
	  arrayValues[11] = jet80_1;
	  arrayValues[12] = jet80_p_1;
	  arrayValues[13] = trgObj_pt_1;
	  arrayValues[14] = centBin;
	  arrayValues[15] = chMax_1[g];
	  arrayValues[16] = chSum_1[g];
	  arrayValues[17] = phMax_1[g];
	  arrayValues[18] = phSum_1[g];
	  arrayValues[19] = neMax_1[g];
	  arrayValues[20] = neSum_1[g];
	  arrayValues[21] = muMax_1[g];
	  arrayValues[22] = muSum_1[g];
	  arrayValues[23] = eMax_1[g];
	  arrayValues[24] = eSum_1[g];
	  arrayValues[25] = trkMax_1[g];
	  arrayValues[26] = trkSum_1[g];
	  jets_ID[k]->Fill(arrayValues);
#endif
	  //if(cut5<0.95 && cut6<0.95 && cut7<0.95 && cut1>0.05 && cut2<0.95 && cut2b<0.95){
	  //if(eSum_1[g]/(chSum_1[g]+neSum_1[g]+phSum_1[g]+muSum_1[g])<0.7){
	  //if((neMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9) && (phMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9) && (chMax_1[g]/pt_1[g]>0.05) && (muMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9) && (chMax_1[g]/(chMax_1[g]+neMax_1[g]+phMax_1[g])<0.9)){
	  //if(1>0){
	  // hpbpb_RecoOverRaw[k][j][centBin]->Fill((Float_t)pt_1[g]/raw_1[g]);
	  // hpbpb_RecoOverRaw[k][j][nbins_cent]->Fill((Float_t)pt_1[g]/raw_1[g]);
	  // hpbpb_RecoOverRaw_rawpt[k][j][centBin]->Fill(raw__1[g],(Float_t)pt_1[g]/raw_1[g]);
	  // hpbpb_RecoOverRaw_rawpt[k][j][nbins_cent]->Fill(raw_1[g],(Float_t)pt_1[g]/raw_1[g]);
	  //get the spectra histograms to make the ratios for the different Jet ID cuts. 0-30%
	  //if(centBin==3 || centBin==4 || centBin==5) continue;
	  if(k==0){
	  
	    hpbpb_chMaxJtpt_jtpt->Fill(pt_1[g],chMax_1[g]/pt_1[g]);
	    hpbpb_eMaxJtpt_jtpt->Fill(pt_1[g],eMax_1[g]/pt_1[g]);
	    hpbpb_eMaxJtpt_chMaxJtpt->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/pt_1[g]);
	    hpbpb_eMaxSumcand_jtpt->Fill(pt_1[g],eMax_1[g]/(chSum_1[g]+phSum_1[g]+neSum_1[g]+muSum_1[g]));
	    hpbpb_eMaxSumcand_chMaxJtpt->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/(chSum_1[g]+phSum_1[g]+neSum_1[g]+muSum_1[g]));
	    for(int a = 0;a<ptSelection;a++) {
	      if(jet80_1){
		if(pt_1[g] >= ptBoundary[a] && pt_1[g] < ptBoundary[a+1]){
		  hpbpb_chMaxJtpt_jtpt_ptselection[2][a]->Fill(pt_1[g],chMax_1[g]/pt_1[g]);
		  hpbpb_phMaxJtpt_jtpt_ptselection[2][a]->Fill(pt_1[g],phMax_1[g]/pt_1[g]);
		  hpbpb_neMaxJtpt_jtpt_ptselection[2][a]->Fill(pt_1[g],neMax_1[g]/pt_1[g]);
		  hpbpb_muMaxJtpt_jtpt_ptselection[2][a]->Fill(pt_1[g],muMax_1[g]/pt_1[g]);
		  hpbpb_eMaxJtpt_jtpt_ptselection[2][a]->Fill(pt_1[g],eMax_1[g]/pt_1[g]);
		  hpbpb_eMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/pt_1[g]);
		  hpbpb_phMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(chMax_1[g]/pt_1[g],phMax_1[g]/pt_1[g]);
		  hpbpb_neMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(chMax_1[g]/pt_1[g],neMax_1[g]/pt_1[g]);
		  hpbpb_muMaxJtpt_chMaxJtpt_ptselection[2][a]->Fill(chMax_1[g]/pt_1[g],muMax_1[g]/pt_1[g]);
		  hpbpb_eMaxSumcand_chMaxJtpt_ptselection[2][a]->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		  hpbpb_eMaxSumcand_eMaxJtpt_ptselection[2][a]->Fill(eMax_1[g]/pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		  hpbpb_eMaxSumcand_jtpt_ptselection[2][a]->Fill(pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		}
	      }
	      if(jet65_1 && !jet80_1){
		if(pt_1[g] >= ptBoundary[a] && pt_1[g] < ptBoundary[a+1]){
		  hpbpb_chMaxJtpt_jtpt_ptselection[1][a]->Fill(pt_1[g],chMax_1[g]/pt_1[g]);
		  hpbpb_phMaxJtpt_jtpt_ptselection[1][a]->Fill(pt_1[g],phMax_1[g]/pt_1[g]);
		  hpbpb_neMaxJtpt_jtpt_ptselection[1][a]->Fill(pt_1[g],neMax_1[g]/pt_1[g]);
		  hpbpb_muMaxJtpt_jtpt_ptselection[1][a]->Fill(pt_1[g],muMax_1[g]/pt_1[g]);
		  hpbpb_eMaxJtpt_jtpt_ptselection[1][a]->Fill(pt_1[g],eMax_1[g]/pt_1[g]);
		  hpbpb_eMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/pt_1[g]);
		  hpbpb_phMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(chMax_1[g]/pt_1[g],phMax_1[g]/pt_1[g]);
		  hpbpb_neMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(chMax_1[g]/pt_1[g],neMax_1[g]/pt_1[g]);
		  hpbpb_muMaxJtpt_chMaxJtpt_ptselection[1][a]->Fill(chMax_1[g]/pt_1[g],muMax_1[g]/pt_1[g]);
		  hpbpb_eMaxSumcand_chMaxJtpt_ptselection[1][a]->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		  hpbpb_eMaxSumcand_eMaxJtpt_ptselection[1][a]->Fill(eMax_1[g]/pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		  hpbpb_eMaxSumcand_jtpt_ptselection[1][a]->Fill(pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		}
	      }
	      if(jet55_1 && !jet65_1 && !jet80_1){
		if(pt_1[g] >= ptBoundary[a] && pt_1[g] < ptBoundary[a+1]){
		  hpbpb_chMaxJtpt_jtpt_ptselection[0][a]->Fill(pt_1[g],chMax_1[g]/pt_1[g]);
		  hpbpb_phMaxJtpt_jtpt_ptselection[0][a]->Fill(pt_1[g],phMax_1[g]/pt_1[g]);
		  hpbpb_neMaxJtpt_jtpt_ptselection[0][a]->Fill(pt_1[g],neMax_1[g]/pt_1[g]);
		  hpbpb_muMaxJtpt_jtpt_ptselection[0][a]->Fill(pt_1[g],muMax_1[g]/pt_1[g]);
		  hpbpb_eMaxJtpt_jtpt_ptselection[0][a]->Fill(pt_1[g],eMax_1[g]/pt_1[g]);
		  hpbpb_eMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/pt_1[g]);
		  hpbpb_phMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(chMax_1[g]/pt_1[g],phMax_1[g]/pt_1[g]);
		  hpbpb_neMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(chMax_1[g]/pt_1[g],neMax_1[g]/pt_1[g]);
		  hpbpb_muMaxJtpt_chMaxJtpt_ptselection[0][a]->Fill(chMax_1[g]/pt_1[g],muMax_1[g]/pt_1[g]);
		  hpbpb_eMaxSumcand_chMaxJtpt_ptselection[0][a]->Fill(chMax_1[g]/pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		  hpbpb_eMaxSumcand_eMaxJtpt_ptselection[0][a]->Fill(eMax_1[g]/pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		  hpbpb_eMaxSumcand_jtpt_ptselection[0][a]->Fill(pt_1[g],eMax_1[g]/(chSum_1[g] + phSum_1[g] + neSum_1[g] + muSum_1[g]));
		}
	      }
	    }
	    if(jet80_1){
	      
		hpbpb_Jet[2][centBin]->Fill(pt_1[g]);
		for(int a = 0;a<5;++a)
		  if(chMax_1[g]/pt_1[g] > (Float_t)(a+1)/100) hpbpb_Jet_chMaxJtpt[2][a][centBin]->Fill(pt_1[g]);
		for(int a = 0;a<9;++a)
		  if(eMax_1[g]/pt_1[g] < (Float_t) (a+1)/10) hpbpb_Jet_eMaxJtpt[2][a][centBin]->Fill(pt_1[g]);
		for(int a = 0;a<3;++a)
		  for(int b = 0;b<3;++b)
		    if(chMax_1[g]/pt_1[g]> (Float_t) (b+1)/100 && eMax_1[g]/pt_1[g] < (Float_t) (a+5)/10 ) hpbpb_Jet_eMaxJtpt_chMaxJtpt[2][a][b][centBin]->Fill(pt_1[g]);
	    }
	    if(jet65_1 && !jet80_1){
		hpbpb_Jet[1][centBin]->Fill(pt_1[g]);
		for(int a = 0;a<5;++a)
		  if(chMax_1[g]/pt_1[g] > (Float_t)(a+1)/100) hpbpb_Jet_chMaxJtpt[1][a][centBin]->Fill(pt_1[g]);
		for(int a = 0;a<9;++a)
		  if(eMax_1[g]/pt_1[g] < (Float_t) (a+1)/10) hpbpb_Jet_eMaxJtpt[1][a][centBin]->Fill(pt_1[g]);
		for(int a = 0;a<3;++a)
		  for(int b = 0;b<3;++b)
		    if(chMax_1[g]/pt_1[g]> (Float_t) (b+1)/100 && eMax_1[g]/pt_1[g] < (Float_t) (a+5)/10 ) hpbpb_Jet_eMaxJtpt_chMaxJtpt[1][a][b][centBin]->Fill(pt_1[g]);
	    }
	    if(jet55_1 && !jet65_1 && !jet80_1){
	
	      hpbpb_Jet[0][centBin]->Fill(pt_1[g],effecPrescl);
		for(int a = 0;a<5;++a)
		  if(chMax_1[g]/pt_1[g] > (Float_t)(a+1)/100) hpbpb_Jet_chMaxJtpt[0][a][centBin]->Fill(pt_1[g],effecPrescl);
		for(int a = 0;a<9;++a)
		  if(eMax_1[g]/pt_1[g] < (Float_t) (a+1)/10) hpbpb_Jet_eMaxJtpt[0][a][centBin]->Fill(pt_1[g],effecPrescl);
		for(int a = 0;a<3;++a)
		  for(int b = 0;b<3;++b)
		    if(chMax_1[g]/pt_1[g]> (Float_t) (b+1)/100 && eMax_1[g]/pt_1[g] < (Float_t) (a+5)/10 ) hpbpb_Jet_eMaxJtpt_chMaxJtpt[0][a][b][centBin]->Fill(pt_1[g],effecPrescl);		
	  
	    }
	    
	  }// check if its r = 0.3
	
#if 0
	  hpbpb_chMax[k][j][centBin]->Fill(chMax_1[g]);
	  hpbpb_phMax[k][j][centBin]->Fill(phMax_1[g]);
	  hpbpb_neMax[k][j][centBin]->Fill(neMax_1[g]);
	  hpbpb_muMax[k][j][centBin]->Fill(muMax_1[g]);
	  hpbpb_eMax[k][j][centBin]->Fill(eMax_1[g]);
	  hpbpb_chSum[k][j][centBin]->Fill(chSum_1[g]);
	  hpbpb_phSum[k][j][centBin]->Fill(phSum_1[g]);
	  hpbpb_neSum[k][j][centBin]->Fill(neSum_1[g]);
	  hpbpb_muSum[k][j][centBin]->Fill(muSum_1[g]);
	  hpbpb_eSum[k][j][centBin]->Fill(eSum_1[g]);
	  // hpbpb_chMax[k][j][nbins_cent]->Fill(chMax_1[g]);
	  // hpbpb_phMax[k][j][nbins_cent]->Fill(phMax_1[g]);
	  // hpbpb_neMax[k][j][nbins_cent]->Fill(neMax_1[g]);
	  // hpbpb_muMax[k][j][nbins_cent]->Fill(muMax_1[g]);
	  // hpbpb_eMax[k][j][nbins_cent]->Fill(eMax_1[g]);
	  // hpbpb_chSum[k][j][nbins_cent]->Fill(chSum_1[g]);
	  // hpbpb_phSum[k][j][nbins_cent]->Fill(phSum_1[g]);
	  // hpbpb_neSum[k][j][nbins_cent]->Fill(neSum_1[g]);
	  // hpbpb_muSum[k][j][nbins_cent]->Fill(muSum_1[g]);
	  // hpbpb_eSum[k][j][nbins_cent]->Fill(eSum_1[g]);
	  if(chMax_1[g]/pt_1[g]>0.02 && eMax_1[g]/pt_1[g]<0.6){
	  hpbpb_chMax_withCut[k][j][centBin]->Fill(chMax_1[g]);
	  hpbpb_phMax_withCut[k][j][centBin]->Fill(phMax_1[g]);
	  hpbpb_neMax_withCut[k][j][centBin]->Fill(neMax_1[g]);
	  hpbpb_muMax_withCut[k][j][centBin]->Fill(muMax_1[g]);
	  hpbpb_eMax_withCut[k][j][centBin]->Fill(eMax_1[g]);
	  hpbpb_chSum_withCut[k][j][centBin]->Fill(chSum_1[g]);
	  hpbpb_phSum_withCut[k][j][centBin]->Fill(phSum_1[g]);
	  hpbpb_neSum_withCut[k][j][centBin]->Fill(neSum_1[g]);
	  hpbpb_muSum_withCut[k][j][centBin]->Fill(muSum_1[g]);
	  hpbpb_eSum_withCut[k][j][centBin]->Fill(eSum_1[g]);
	    //	    hJets->Fill(1);
	    //if(raw_1[g] < 30) continue;
	    if(jet80_1==1 && L1_sj52_1==1) {
	      //if(trgObj_pt_1>=80){
		hpbpb_TrgObj80[k][j][centBin]->Fill(pt_1[g]);
		hpbpb_TrgObj80[k][j][nbins_cent]->Fill(pt_1[g]);
#if 0
		if(jetCounter>=7){
		  hpbpb_TrgObj80_nJet_g7[k][j][centBin]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		  hpbpb_TrgObj80_nJet_g7[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		}
		if(jetCounter<7){
		  hpbpb_TrgObj80_nJet_l7[k][j][centBin]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		  hpbpb_TrgObj80_nJet_l7[k][j][nbins_cent]->Fill(pt_1[g],jet80_p_1*L1_sj52_p_1);
		}
#endif
		// }
#if 0
	      if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
		hpbpb_Jet80[1][1][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[1][1][nbins_cent]->Fill(pt_1[g]);
	      }
	      if(TMath::Abs(Vs_0_x)<v0_tight && TMath::Abs(Vs_0_y)<v0_tight && TMath::Abs(Vs_1_x)<v1_tight && TMath::Abs(Vs_1_y)<v1_tight && TMath::Abs(Vs_2_x)<v2_tight && TMath::Abs(Vs_2_y)<v2_tight && TMath::Abs(Vs_3_x)<v3_tight && TMath::Abs(Vs_3_y)<v3_tight && TMath::Abs(Vs_4_x)<v4_tight && TMath::Abs(Vs_4_y)<v4_tight){
		hpbpb_Jet80[0][1][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[0][1][nbins_cent]->Fill(pt_1[g]);
	      }
	      if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
		hpbpb_Jet80[1][0][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[1][0][nbins_cent]->Fill(pt_1[g]);
	      }
	      if(TMath::Abs(Vs_0_x)<v0_loose && TMath::Abs(Vs_0_y)<v0_loose && TMath::Abs(Vs_1_x)<v1_loose && TMath::Abs(Vs_1_y)<v1_loose && TMath::Abs(Vs_2_x)<v2_loose && TMath::Abs(Vs_2_y)<v2_loose && TMath::Abs(Vs_3_x)<v3_loose && TMath::Abs(Vs_3_y)<v3_loose && TMath::Abs(Vs_4_x)<v4_loose && TMath::Abs(Vs_4_y)<v4_loose){
		hpbpb_Jet80[0][0][centBin]->Fill(pt_1[g]);
		hpbpb_Jet80[0][0][nbins_cent]->Fill(pt_1[g]);
	      }
	      hpbpb_FullJet80[centBin]->Fill(pt_1[g]);
	      hpbpb_FullJet80[nbins_cent]->Fill(pt_1[g]);
#endif
	    }
	    if(jet65_1==1 && L1_sj36_1==1 && jet80_1==0) {
	      // if(trgObj_pt_1>=65 && trgObj_pt_1<80){
	      // 	hpbpb_TrgObj65[k][j][centBin]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	hpbpb_TrgObj65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
		
	      // 	if(jetCounter>=7){
	      // 	  hpbpb_TrgObj65_nJet_g7[k][j][centBin]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	  hpbpb_TrgObj65_nJet_g7[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	}
	      // 	if(jetCounter<7){
	      // 	  hpbpb_TrgObj65_nJet_l7[k][j][centBin]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	  hpbpb_TrgObj65_nJet_l7[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1*L1_sj36_p_1);
	      // 	}
	      // }
	      // hpbpb_Jet65[k][j][centBin]->Fill(pt_1[g],jet65_p_1);
	      // hpbpb_Jet65[k][j][nbins_cent]->Fill(pt_1[g],jet65_p_1);
	    
	      hpbpb_TrgObj65[k][j][centBin]->Fill(pt_1[g]);
	      hpbpb_TrgObj65[k][j][nbins_cent]->Fill(pt_1[g]);
		
#if 0		
		if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
		  hpbpb_Jet65[1][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[1][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_tight && TMath::Abs(Vs_0_y)<v0_tight && TMath::Abs(Vs_1_x)<v1_tight && TMath::Abs(Vs_1_y)<v1_tight && TMath::Abs(Vs_2_x)<v2_tight && TMath::Abs(Vs_2_y)<v2_tight && TMath::Abs(Vs_3_x)<v3_tight && TMath::Abs(Vs_3_y)<v3_tight && TMath::Abs(Vs_4_x)<v4_tight && TMath::Abs(Vs_4_y)<v4_tight){
		  hpbpb_Jet65[0][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[0][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
		  hpbpb_Jet65[1][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[1][0][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_loose && TMath::Abs(Vs_0_y)<v0_loose && TMath::Abs(Vs_1_x)<v1_loose && TMath::Abs(Vs_1_y)<v1_loose && TMath::Abs(Vs_2_x)<v2_loose && TMath::Abs(Vs_2_y)<v2_loose && TMath::Abs(Vs_3_x)<v3_loose && TMath::Abs(Vs_3_y)<v3_loose && TMath::Abs(Vs_4_x)<v4_loose && TMath::Abs(Vs_4_y)<v4_loose){
		  hpbpb_Jet65[0][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet65[0][0][nbins_cent]->Fill(pt_1[g]);
		}
	      
		hpbpb_FullJet65[centBin]->Fill(pt_1[g]);
		hpbpb_FullJet65[nbins_cent]->Fill(pt_1[g]);
#endif
	      
	    }
	    if(jet55_1==1 && L1_sj36_1==1 && jet65_1==0 && jet80_1 == 0) { // passes the jet55 trigger
		hpbpb_TrgObj55[k][j][centBin]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		hpbpb_TrgObj55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
#if 0	      
	      if(trgObj_pt_1>=55 && trgObj_pt_1<65){ // check for the trigger object pt to lie inbetween the two trigger values 
		hpbpb_TrgObj55[k][j][centBin]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		hpbpb_TrgObj55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		
		if(jetCounter>=7){
		  hpbpb_TrgObj55_nJet_g7[k][j][centBin]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		  hpbpb_TrgObj55_nJet_g7[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		}
		if(jetCounter<7){
		  hpbpb_TrgObj55_nJet_l7[k][j][centBin]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		  hpbpb_TrgObj55_nJet_l7[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1*L1_sj36_p_1);
		}
		//if(pt_1[g]>3*trgObj_pt_1){ // get an idea of the event information from these large pt jets 
		//  fHLT_high[k]<<evt_1<<" "<<run_1<<" "<<lumi_1<<" "<<vz_1<<" "<<centBin<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<pt_1[g]<<" "<<raw_1[g]<<" "<<eta_1[g]<<" "<<phi_1[g]<<" "<<endl;
		//}
	      }
	      // hpbpb_Jet55[k][j][centBin]->Fill(pt_1[g],jet55_p_1);
	      // hpbpb_Jet55[k][j][nbins_cent]->Fill(pt_1[g],jet55_p_1);
	      if((jet65_1==0) && (jet80_1==0)){ // this is to just check
		if(TMath::Abs(Vs_0_x)>v0_tight || TMath::Abs(Vs_0_y)>v0_tight || TMath::Abs(Vs_1_x)>v1_tight || TMath::Abs(Vs_1_y)>v1_tight || TMath::Abs(Vs_2_x)>v2_tight || TMath::Abs(Vs_2_y)>v2_tight || TMath::Abs(Vs_3_x)>v3_tight || TMath::Abs(Vs_3_y)>v3_tight || TMath::Abs(Vs_4_x)>v4_tight || TMath::Abs(Vs_4_y)>v4_tight) {
		  hpbpb_Jet55[1][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[1][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_tight && TMath::Abs(Vs_0_y)<v0_tight && TMath::Abs(Vs_1_x)<v1_tight && TMath::Abs(Vs_1_y)<v1_tight && TMath::Abs(Vs_2_x)<v2_tight && TMath::Abs(Vs_2_y)<v2_tight && TMath::Abs(Vs_3_x)<v3_tight && TMath::Abs(Vs_3_y)<v3_tight && TMath::Abs(Vs_4_x)<v4_tight && TMath::Abs(Vs_4_y)<v4_tight){
		  hpbpb_Jet55[0][1][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[0][1][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)>v0_loose || TMath::Abs(Vs_0_y)>v0_loose || TMath::Abs(Vs_1_x)>v1_loose || TMath::Abs(Vs_1_y)>v1_loose || TMath::Abs(Vs_2_x)>v2_loose || TMath::Abs(Vs_2_y)>v2_loose || TMath::Abs(Vs_3_x)>v3_loose || TMath::Abs(Vs_3_y)>v3_loose || TMath::Abs(Vs_4_x)>v4_loose || TMath::Abs(Vs_4_y)>v4_loose) {
		  hpbpb_Jet55[1][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[1][0][nbins_cent]->Fill(pt_1[g]);
		}
		if(TMath::Abs(Vs_0_x)<v0_loose && TMath::Abs(Vs_0_y)<v0_loose && TMath::Abs(Vs_1_x)<v1_loose && TMath::Abs(Vs_1_y)<v1_loose && TMath::Abs(Vs_2_x)<v2_loose && TMath::Abs(Vs_2_y)<v2_loose && TMath::Abs(Vs_3_x)<v3_loose && TMath::Abs(Vs_3_y)<v3_loose && TMath::Abs(Vs_4_x)<v4_loose && TMath::Abs(Vs_4_y)<v4_loose){
		  hpbpb_Jet55[0][0][centBin]->Fill(pt_1[g]);
		  hpbpb_Jet55[0][0][nbins_cent]->Fill(pt_1[g]);
		}
	      
		hpbpb_FullJet55[centBin]->Fill(pt_1[g]);
		hpbpb_FullJet55[nbins_cent]->Fill(pt_1[g]);
	      
	      }
#endif
	    }
	      
	  }//qa cut selection
	
#endif
	}//jet loop
	  
      }//eta bin loop
#endif
    }//nentries_jet55or65 loop
    
  }//radius loop. 
  
  //hEvents_HLT80->Write();
  //hEvents_HLT65->Write();
  //hEvents_HLT55->Write();

  #if 0
  for(int i = 0;i<nbins_cent;++i){
    
    for(int t = 0;t<trigValue-1;t++){
      hpbpb_Jet[t][i]->Write();
      for(int a = 0;a<5;++a){
	hpbpb_Jet_chMaxJtpt[t][a][i]->Write();
      }
      for(int a = 0;a<9;++a){
	hpbpb_Jet_eMaxJtpt[t][a][i]->Write();
      }
      for(int a = 0;a<3;++a)
	for(int b = 0;b<3;++b) { 
	  hpbpb_Jet_eMaxJtpt_chMaxJtpt[t][a][b][i]->Write();
	}
    }
  }
  hpbpb_chMaxJtpt_jtpt->Write();
  hpbpb_eMaxJtpt_jtpt->Write();
  hpbpb_eMaxJtpt_chMaxJtpt->Write();
  hpbpb_eMaxSumcand_jtpt->Write();
  hpbpb_eMaxSumcand_chMaxJtpt->Write();
  for(int a = 0;a<ptSelection;++a){
    for(int t = 0;t<trigValue-1;++t){
      hpbpb_eMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_phMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_neMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_muMaxJtpt_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_eMaxSumcand_chMaxJtpt_ptselection[t][a]->Write();
      hpbpb_eMaxSumcand_jtpt_ptselection[t][a]->Write();
      hpbpb_chMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_phMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_neMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_muMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_eMaxJtpt_jtpt_ptselection[t][a]->Write();
      hpbpb_eMaxSumcand_eMaxJtpt_ptselection[t][a]->Write();
    }
  }
  
  for(int k = 0;k<no_radius;k++){
    for(int j = 0;j<nbins_eta;j++){
      
      for(int i = 0;i<nbins_cent;i++){
	
	//jets_ID[p]->Write();
	hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj80[k][j][i]);
	hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj65[k][j][i]);
	hpbpb_TrgObjComb[k][j][i]->Add(hpbpb_TrgObj55[k][j][i]);
	hpbpb_TrgObjComb[k][j][i]->Write();
	if(printDebug)hpbpb_TrgObjComb[k][j][i]->Print();
	hpbpb_TrgObj80[k][j][i]->Write();
	if(printDebug)hpbpb_TrgObj80[k][j][i]->Print();
	hpbpb_TrgObj65[k][j][i]->Write();
	if(printDebug)hpbpb_TrgObj65[k][j][i]->Print();
	hpbpb_TrgObj55[k][j][i]->Write();
	if(printDebug)hpbpb_TrgObj55[k][j][i]->Print();
	
	// hpbpb_RecoOverRaw_jtpt[k][j][i]->Write();
	// if(printDebug)hpbpb_RecoOverRaw_jtpt[k][j][i]->Print();
	// hpbpb_RecoOverRaw[k][j][i]->Write();
	// if(printDebug)hpbpb_RecoOverRaw[k][j][i]->Print();
	
	hpbpb_chMax[k][j][i]->Write();
	hpbpb_phMax[k][j][i]->Write();
	hpbpb_neMax[k][j][i]->Write();
	hpbpb_muMax[k][j][i]->Write();
	hpbpb_eMax[k][j][i]->Write();
	hpbpb_chSum[k][j][i]->Write();
	hpbpb_phSum[k][j][i]->Write();
	hpbpb_neSum[k][j][i]->Write();
	hpbpb_muSum[k][j][i]->Write();
	hpbpb_eSum[k][j][i]->Write();
	hpbpb_chMax_withCut[k][j][i]->Write();
	hpbpb_phMax_withCut[k][j][i]->Write();
	hpbpb_neMax_withCut[k][j][i]->Write();
	hpbpb_muMax_withCut[k][j][i]->Write();
	hpbpb_eMax_withCut[k][j][i]->Write();
	hpbpb_chSum_withCut[k][j][i]->Write();
	hpbpb_phSum_withCut[k][j][i]->Write();
	hpbpb_neSum_withCut[k][j][i]->Write();
	hpbpb_muSum_withCut[k][j][i]->Write();
	hpbpb_eSum_withCut[k][j][i]->Write();
      }
    }
    //evt_electron_failure[k]->Write();
    //evt_electron_good[k]->Write();
  }
  
  hpbpb_cent[0]->Write();
  hpbpb_vz[0]->Write();
  hpbpb_vx[0]->Write();
  hpbpb_vy[0]->Write();
#endif
  
  hCentEvents->Write();
  /*
  hEvents->Write();
  hEvents_pCES->Write();
  hEvents_pHBHE->Write();
  hEvents_vz15->Write();
  hEvents_supernova->Write();
  hEvents_eta2->Write();
  */
  f.Write();
  f.Close();

 
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;

}

