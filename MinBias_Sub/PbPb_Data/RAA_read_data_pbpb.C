// Raghav Kunnawalkam Elayavalli
// June 5th 2014
// CERN
// for questions or comments: raghav.k.e at CERN dot CH

// 
// read all the MC files for PbPb and pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 
// 

// July 19 - all pp histograms will have 2D arrays with [radius][eta_bin]. the PbPb histograms will be defined by 3D arrays with [radius][eta_bin][centrality]. 
// July 20 - the loop structure(s) are defined as follows for the several histograms 
//            Radius    : iteration variable: k;       number of iterations: no_radius;                                 values for the radii: list_radius
//            Eta bins  : iteration variable: j;       number of iterations: nbins_eta;                                 values of the bins  : boundaries_eta
//            Centrality: iteration variable: i;       number of iterations: nbins_cent +1;                             values of the bins  : boundaries_cent (till nbins_cent) + 0-200 (the whole range) for the final iteration. 
//            p_T Hats  : iteration variable: h;       number of iterations: nbins_pthat (PbPb) and nbinsPP_pthat (pp); values of the bins  : boundaries_pthat (PbPb) and boundariesPP_pthat (pp)  
//            jets      : iteration variable: g;       number of iterations: no of jets in Data[k][h];
//            p_T       : defined just below as nbins_pt with 39 bins. to match our NLO and jet RpA analysis bins. 

// Oct 23 - removed the cuts from the MC -> like the noisefilter etc... 

// Nov 4th - added the supernova event cut rejection based on the no of hits in the pixel. 

// Dec 9th - going to PU for the Jet RAA. 

// Dec 17th - changing the file list to smaller 50k files on which JEC were derived to check for PF electron problems, requested by Marguerite.

// Jan 13th 2015 - adding in the official pp mc (from Dragos) 
//               - this is going to be a bit tricky since each file is split up into 4 smaller files. so each pthat will have a TChain!


// Feb 12th - cleaned up the macro to make it usable (hopefuly) by others.

// Jun 22th - going back to the HiForest with the trees from Pawan's for the event selection cuts and PF electron cuts. 

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

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 3; 
//static const int list_radius[no_radius] = {2,3,4};

static const int nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 
static const Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
static const int trigValue = 4;
static const char trigName [trigValue][256] = {"HLT55","HLT65","HLT80","Combined"};
static const Float_t effecPrescl = 2.047507;
static const char * etaWidth = (char*)"20_eta_20";

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

struct Jet{
  int id;
  float pt;
};
bool compare_pt(Jet jet1, Jet jet2);
bool compare_pt(Jet jet1, Jet jet2){
  return jet1.pt > jet2.pt;
}

using namespace std;

void RAA_read_data_pbpb(int startfile = 0,
		      int endfile = 1,
		      int radius = 3,
		      std::string kFoname="test_output.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;
  std::string infile_Select;

  infile_Forest = "jetRAA_PbPb_data_forest.txt";
  infile_Select = Form("jetRAA_PbPb_data_akPu%d_select.txt",radius);
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  std::ifstream instr_Select(infile_Select.c_str(),std::ifstream::in);
  std::string filename_Select;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
    instr_Select>>filename_Select;
  }

  const int N = 5; //6

  TChain * jetpbpb[N];
  TChain * evt_select, * jet_select; 

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("akPu%dPFJetAnalyzer",radius);
  dir[3] = "akPu3CaloJetAnalyzer";
  dir[4] = "hiEvtAnalyzer";
  // dir[4] = "hltobject";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "t",
    "HiTree"
    // , "jetObjTree"
  };

  for(int t = 0;t<N;t++){
    jetpbpb[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  evt_select = new TChain(Form("akPu%dJetAnalyzer/evtTree",radius));
  jet_select = new TChain(Form("akPu%dJetAnalyzer/jetTree",radius));

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;
    instr_Select>>filename_Select;

    if(printDebug)cout<<"HiForest filename = "<<filename_Forest.c_str()<<endl;

    jetpbpb[0]->Add(filename_Forest.c_str());
    jetpbpb[1]->Add(filename_Forest.c_str());
    jetpbpb[2]->Add(filename_Forest.c_str());
    jetpbpb[3]->Add(filename_Forest.c_str());
    jetpbpb[4]->Add(filename_Forest.c_str());
    jet_select->Add(filename_Select.c_str());
    evt_select->Add(filename_Select.c_str());
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[3]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[4]->GetEntries() << endl;
    if(printDebug)cout << "jet selection file " << jet_select->GetEntries() << endl;
    if(printDebug)cout << "event selection file" << evt_select->GetEntries() << endl;

  }
  
  jetpbpb[2]->AddFriend(jetpbpb[0]);
  jetpbpb[2]->AddFriend(jetpbpb[1]);
  jetpbpb[2]->AddFriend(jetpbpb[4]);
  jetpbpb[3]->AddFriend(jetpbpb[0]);
  jetpbpb[3]->AddFriend(jetpbpb[1]);
  jetpbpb[3]->AddFriend(jetpbpb[4]);
  
  jetpbpb[2]->AddFriend(evt_select);
  jetpbpb[3]->AddFriend(evt_select);

  // Forest files 
  int nref_F;
  float pt_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float chMax_F[1000];
  float trkMax_F[1000];
  float chSum_F[1000];
  float phSum_F[1000];
  float neSum_F[1000];
  float trkSum_F[1000];
  float phMax_F[1000];
  float neMax_F[1000];
  float eMax_F[1000];
  float muMax_F[1000];
  float eSum_F[1000];
  float muSum_F[1000];
  float jtpu_F[1000];
  int jet55_F;
  int jet65_F;
  int jet80_F;
  int L1_sj36_F;
  int L1_sj52_F;
  int L1_sj36_p_F;
  int L1_sj52_p_F;
  int jet55_p_F;
  int jet65_p_F;
  int jet80_p_F;
  float vz_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int hiBin_F;
  int pcollisionEventSelection_F;
  int pHBHENoiseFilter_F;

  float calopt_F[1000];
  jetpbpb[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jetpbpb[4]->SetBranchAddress("evt",&evt_F);
  jetpbpb[4]->SetBranchAddress("run",&run_F);
  jetpbpb[4]->SetBranchAddress("lumi",&lumi_F);
  jetpbpb[4]->SetBranchAddress("hiBin",&hiBin_F);
  jetpbpb[4]->SetBranchAddress("vz",&vz_F);
  jetpbpb[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  jetpbpb[1]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  //jetpbpb[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  //jetpbpb[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpbpb[2]->SetBranchAddress("nref",&nref_F);
  jetpbpb[2]->SetBranchAddress("jtpt",pt_F);
  jetpbpb[2]->SetBranchAddress("jteta",eta_F);
  jetpbpb[2]->SetBranchAddress("jtphi",phi_F);
  jetpbpb[2]->SetBranchAddress("rawpt",rawpt_F);
  jetpbpb[2]->SetBranchAddress("jtpu",jtpu_F);
  jetpbpb[2]->SetBranchAddress("chargedMax",chMax_F);
  jetpbpb[2]->SetBranchAddress("chargedSum",chSum_F);
  jetpbpb[2]->SetBranchAddress("trackMax",trkMax_F);
  jetpbpb[2]->SetBranchAddress("trackSum",trkSum_F);
  jetpbpb[2]->SetBranchAddress("photonMax",phMax_F);
  jetpbpb[2]->SetBranchAddress("photonSum",phSum_F);
  jetpbpb[2]->SetBranchAddress("neutralMax",neMax_F);
  jetpbpb[2]->SetBranchAddress("neutralSum",neSum_F);
  jetpbpb[2]->SetBranchAddress("eSum",eSum_F);
  jetpbpb[2]->SetBranchAddress("eMax",eMax_F);
  jetpbpb[2]->SetBranchAddress("muSum",muSum_F);
  jetpbpb[2]->SetBranchAddress("muMax",muMax_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v1",&jet55_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v1",&jet65_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v1",&jet80_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

  // event tree selection file:
  int hiBin_eS;
  int run_eS;
  int evt_eS;
  int lumi_eS;
  float vz_eS; 
  int isGoodEvent_eS;
  int nref_eS;
  int isMiMatch_eS[1000];
  int isClMatch_eS[1000];
  int isPFElecCut_eS[1000];
  int isTrackCut_eS[1000];
  int isMuCut_eS[1000];
  int index_eS[1000];
  float pt_eS[1000];
  float calopt_eS[1000];
  Double_t weight_eS;

  evt_select->SetBranchAddress("hiBin",&hiBin_eS);
  evt_select->SetBranchAddress("run_value",&run_eS);
  evt_select->SetBranchAddress("evt_value",&evt_eS);
  evt_select->SetBranchAddress("lumi_value",&lumi_eS);
  evt_select->SetBranchAddress("vz",&vz_eS);
  evt_select->SetBranchAddress("weight", &weight_eS);  
  evt_select->SetBranchAddress("isGoodEvt",&isGoodEvent_eS);
  evt_select->SetBranchAddress("nref",&nref_eS);
  evt_select->SetBranchAddress("index", index_eS);
  evt_select->SetBranchAddress("isMlMatch",isMiMatch_eS);
  evt_select->SetBranchAddress("isClMatch",isClMatch_eS);
  evt_select->SetBranchAddress("isPFElecCut",isPFElecCut_eS);
  evt_select->SetBranchAddress("isTrackCut",isTrackCut_eS);
  evt_select->SetBranchAddress("isMuCut",isMuCut_eS);
  evt_select->SetBranchAddress("pfpt",pt_eS);
  evt_select->SetBranchAddress("calopt",calopt_eS);

  // jet tree selection file:
  int hiBin_jS;
  int run_jS;
  int evt_jS;
  int lumi_jS;
  float vz_jS;
  int nref_jS;
  float pt_jS[1000];
  float eta_jS[1000];
  float eMax_jS[1000];
  float chSum_jS[1000];
  float phSum_jS[1000];
  float neSum_jS[1000];
  float muSum_jS[1000];
  float calopt_jS[1000];
  int   isCaloMatch_jS[1000];
  int   isMultiMatch_jS[1000];

  jet_select->SetBranchAddress("hiBin",&hiBin_jS);
  jet_select->SetBranchAddress("run_value",&run_jS);
  jet_select->SetBranchAddress("evt_value",&evt_jS);
  jet_select->SetBranchAddress("lumi_value",&lumi_jS);
  jet_select->SetBranchAddress("vz",&vz_jS);
  jet_select->SetBranchAddress("npf", &nref_jS);  
  jet_select->SetBranchAddress("pfpt", &pt_jS);  
  jet_select->SetBranchAddress("eMax", &eMax_jS);  
  jet_select->SetBranchAddress("pfeta", &eta_jS);  
  jet_select->SetBranchAddress("calopt", &calopt_jS);  
  jet_select->SetBranchAddress("isCaloMatch", &isCaloMatch_jS);  
  jet_select->SetBranchAddress("isMultiMatch", &isMultiMatch_jS);  
  jet_select->SetBranchAddress("chSum", &chSum_jS);  
  jet_select->SetBranchAddress("phSum", &phSum_jS);  
  jet_select->SetBranchAddress("neSum", &neSum_jS);  
  jet_select->SetBranchAddress("muSum", &muSum_jS);  

  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

  TH1F *hpbpb_TrgObj80[nbins_cent];
  TH1F *hpbpb_TrgObj65[nbins_cent];
  TH1F *hpbpb_TrgObj55[nbins_cent];
  TH1F *hpbpb_TrgObjComb[nbins_cent];

  TH1F *hpbpb_raw_TrgObj80[nbins_cent];
  TH1F *hpbpb_raw_TrgObj65[nbins_cent];
  TH1F *hpbpb_raw_TrgObj55[nbins_cent];
  TH1F *hpbpb_raw_TrgObjComb[nbins_cent];

  TH1F *hpbpb_anaBin_TrgObj80[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObj65[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObj55[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObjComb[nbins_cent];

  TH1F *hpbpb_JEC_TrgObj80[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj65[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj55[nbins_cent];
  TH1F *hpbpb_JEC_TrgObjComb[nbins_cent];

  TH1F *hpbpb_Smear_TrgObj80[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj65[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj55[nbins_cent];
  TH1F *hpbpb_Smear_TrgObjComb[nbins_cent];

  TH1F * hcentrality_Jet55[nbins_cent];
  TH1F * hcentrality_Jet55not65not80[nbins_cent];
  TH1F * hcentrality_Jet65[nbins_cent];
  TH1F * hcentrality_Jet65not80[nbins_cent];
  TH1F * hcentrality_Jet80[nbins_cent];
  TH1F * hcentrality_Jet55_Prescl[nbins_cent];
  TH1F * hcentrality_Jet55not65not80_Prescl[nbins_cent];
  TH1F * hcentrality_Jet65_Prescl[nbins_cent];
  TH1F * hcentrality[nbins_cent];


  for(int i = 0;i<nbins_cent;++i){

    hcentrality[i] = new TH1F(Form("hcentrality_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55[i] = new TH1F(Form("hcentrality_Jet55_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55_Prescl[i] = new TH1F(Form("hcentrality_Jet55_Prescl_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55not65not80[i] = new TH1F(Form("hcentrality_Jet55not65not80_cent%d",i),"",200, 0, 200);
    hcentrality_Jet55not65not80_Prescl[i] = new TH1F(Form("hcentrality_Jet55not65not80_Prescl_cent%d",i),"",200, 0, 200);
    hcentrality_Jet65[i] = new TH1F(Form("hcentrality_Jet65_cent%d",i),"",200, 0, 200);
    hcentrality_Jet65_Prescl[i] = new TH1F(Form("hcentrality_Jet65_Prescl_cent%d",i),"",200, 0, 200);
    hcentrality_Jet65not80[i] = new TH1F(Form("hcentrality_Jet65not80_cent%d",i),"",200, 0, 200);
    hcentrality_Jet80[i] = new TH1F(Form("hcentrality_Jet80_cent%d",i),"",200, 0, 200);

    hpbpb_TrgObj80[i] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObj65[i] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObj55[i] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObjComb[i] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_raw_TrgObj80[i] = new TH1F(Form("hpbpb_raw_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_raw_TrgObj65[i] = new TH1F(Form("hpbpb_raw_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_raw_TrgObj55[i] = new TH1F(Form("hpbpb_raw_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_raw_TrgObjComb[i] = new TH1F(Form("hpbpb_raw_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_anaBin_TrgObj80[i] = new TH1F(Form("hpbpb_anaBin_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObj65[i] = new TH1F(Form("hpbpb_anaBin_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObj55[i] = new TH1F(Form("hpbpb_anaBin_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObjComb[i] = new TH1F(Form("hpbpb_anaBin_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    
    hpbpb_JEC_TrgObj80[i] = new TH1F(Form("hpbpb_JEC_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JEC_TrgObj65[i] = new TH1F(Form("hpbpb_JEC_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JEC_TrgObj55[i] = new TH1F(Form("hpbpb_JEC_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JEC_TrgObjComb[i] = new TH1F(Form("hpbpb_JEC_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_Smear_TrgObj80[i] = new TH1F(Form("hpbpb_Smear_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Smear_TrgObj65[i] = new TH1F(Form("hpbpb_Smear_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Smear_TrgObj55[i] = new TH1F(Form("hpbpb_Smear_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Smear_TrgObjComb[i] = new TH1F(Form("hpbpb_Smear_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

  }

  // get the HLT information per forest file.
  TH1F * hJet55_File = new TH1F("hJet55_File","",902, 0, 902);
  TH1F * hJet65_File = new TH1F("hJet65_File","",902, 0, 902);
  TH1F * hJet80_File = new TH1F("hJet80_File","",902, 0, 902);
  TH1F * hJet65not80_File = new TH1F("hJet65not80_File","",902, 0, 902);
  TH1F * hJet55not6580_File = new TH1F("hJet55not6580_File","",902, 0, 902);
  
  // now start the event loop for each file. 
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpbpb[0]->GetEntries();
  Long64_t nGoodEvt = 0;
  if(printDebug) nentries = 10;
  TRandom rnd; 

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    jetpbpb[0]->GetEntry(nEvt);
    jetpbpb[1]->GetEntry(nEvt);
    jetpbpb[2]->GetEntry(nEvt);
    jetpbpb[4]->GetEntry(nEvt);
    jetpbpb[3]->GetEntry(nEvt);
    evt_select->GetEntry(nEvt);

    if(printDebug) cout<<"forest values = "<<hiBin_F<<", "<<evt_F<<", "<<run_F<<", "<<lumi_F<<", "<<vz_F<<endl;
    
    if(!isGoodEvent_eS) continue; 
    
    jet_select->GetEntry(nGoodEvt);
    ++nGoodEvt;

    int cBin = findBin(hiBin_F);//tells us the centrality of the event.
    if(cBin == -1) continue;
    if(cBin == nbins_cent) continue;
    
    hcentrality[cBin]->Fill(hiBin_F);
    if(jet55_F == 1) hcentrality_Jet55[cBin]->Fill(hiBin_F);
    if(jet55_F == 1) hcentrality_Jet55_Prescl[cBin]->Fill(hiBin_F, jet55_p_F);
    if(jet65_F == 1) hcentrality_Jet65[cBin]->Fill(hiBin_F);
    if(jet65_F == 1) hcentrality_Jet65_Prescl[cBin]->Fill(hiBin_F, jet65_p_F);
    if(jet80_F == 1) hcentrality_Jet80[cBin]->Fill(hiBin_F);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80[cBin]->Fill(hiBin_F);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80_Prescl[cBin]->Fill(hiBin_F, jet55_p_F);
    if(jet65_F == 1 && jet80_F == 0) hcentrality_Jet65not80[cBin]->Fill(hiBin_F);

    

    if(jet55_F == 1) hJet55_File->Fill(startfile);
    if(jet65_F == 1) hJet65_File->Fill(startfile);
    if(jet80_F == 1) hJet80_File->Fill(startfile);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hJet55not6580_File->Fill(startfile);
    if(jet65_F == 1 && jet80_F == 0) hJet65not80_File->Fill(startfile);

    if(hiBin_eS != hiBin_F || evt_eS != evt_F || run_eS != run_F || lumi_eS != lumi_F || vz_eS != vz_F) cout<<"ERROR mismatch eS, F"<<endl;
    if(hiBin_eS != hiBin_jS || evt_eS != evt_jS || run_eS != run_jS || lumi_eS != lumi_jS || vz_eS != vz_jS) cout<<"ERROR mismatch eS, jS"<<endl;
    if(hiBin_F != hiBin_jS || evt_F != evt_jS || run_F != run_jS || lumi_F != lumi_jS || vz_F != vz_jS) cout<<"ERROR mismatch F, jS"<<endl;

    if(printDebug) cout<<"hibin hiForest = "<<hiBin_F<<", evtTree = "<<hiBin_eS<<", jetTree = "<<hiBin_jS<<endl;
    if(printDebug) cout<<"evt hiForest   = "<<evt_F<<", evtTree = "<<evt_eS<<", jetTree = "<<evt_jS<<endl;
    if(printDebug) cout<<"lumi hiForest  = "<<lumi_F<<", evtTree = "<<lumi_eS<<", jetTree = "<<lumi_jS<<endl;
    if(printDebug) cout<<"vz hiForest    = "<<vz_F<<", evtTree = "<<vz_eS<<", jetTree = "<<vz_jS<<endl;
    
    if(printDebug) cout<<"nref_F = "<<nref_F<<", nref_eS = "<<nref_eS<<", nref_jS = "<<nref_jS<<endl;

    if(nref_F != nref_eS) cout<<"ERROR mismatch in jet counts"<<endl;

    // //! Sort the jetTree jets according to pT
    // std::vector < Jet > vJet;
    // for(int jet2 = 0; jet2<nref_jS; ++jet2){
    //   //cout <<"\t \t jetTree *** "<< jet2 <<  ", pT " << pt_jS[jet2] <<  ", chSum : "<< chSum_jS[jet2] << endl;
    //   Jet ijet;
    //   ijet.id = jet2;
    //   ijet.pt = pt_jS[jet2];
    //   vJet.push_back(ijet);
    // }
    // std::sort (vJet.begin(), vJet.end(), compare_pt);
    // std::vector < Jet >::const_iterator itr;

    // int jet=0;
    // for(itr=vJet.begin(); itr!=vJet.end(); ++itr, ++jet){

    //   int jetLoc = (*itr).id;
    //   if(isMultiMatch_jS[jetLoc]) {
    // 	++itr;
    // 	jetLoc = (*itr).id;
    // 	if(itr == vJet.end())  break;
    //   }
    //   if(fabs(eta_jS[jetLoc]) > 2) continue;
    //   //if(isPFElecCut_eS[jet] != 1) continue;
    //   // if(isMiMatch_eS[jet]) continue;
    //   if(pt_jS[jetLoc] <15) continue;

    //   bool PFElecCut = false;

    //   Float_t Sumcand = chSum_jS[jetLoc] + phSum_jS[jetLoc] + neSum_jS[jetLoc] + muSum_jS[jetLoc];
    //   if(isCaloMatch_jS[jetLoc] == 1){
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.5 && calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.85 && eMax_jS[jetLoc]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_jS[jetLoc]/pt_jS[jetLoc] - (Float_t)9/7)) PFElecCut = true;
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.85) PFElecCut = true;
    // 	if(calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.5 && eMax_jS[jetLoc]/Sumcand < 0.05) PFElecCut = true;
    //   }
    //   if(isCaloMatch_jS[jetLoc] == 0)
    // 	if(eMax_jS[jetLoc]/Sumcand < 0.05 ) PFElecCut = true;

    //   // if(!PFElecCut) continue;
      
    //   // if(printDebug && (fabs(eta_jS[jet] > 2))) cout<<"jets with |eta| > 2 in jetTree"<<endl;
    //   // if(printDebug && (fabs(eta_F[jet] > 2)))  cout<<"jets with |eta| > 2 in Forest"<<endl;

    Float_t wght = 1; 
      
    //   if(printDebug && index_eS[jet] >= 0 )cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT = "<<calopt_F[index_eS[jet]]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;      // if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;


    for( int jet = 0; jet<nref_F; jet++ ){

      if( fabs(eta_F[jet]) > 2.0 ) continue;
      //if( isPFElecCut_eS[jet] != 1 ) continue;

      //if( chMax_F[jet] < 7 && trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;
      if( trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;

      bool PFElecCut = false;

      Float_t Sumcand = chSum_F[jet] + phSum_F[jet] + neSum_F[jet] + muSum_F[jet];
      if(isClMatch_eS[jet] == 1){
    	if(calopt_eS[jet]/pt_F[jet] > 0.5 && calopt_eS[jet]/pt_F[jet] <= 0.85 && eMax_F[jet]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_eS[jet]/pt_F[jet] - (Float_t)9/7)) PFElecCut = true;
    	if(calopt_eS[jet]/pt_F[jet] > 0.85) PFElecCut = true;
    	if(calopt_eS[jet]/pt_F[jet] <= 0.5 && eMax_F[jet]/Sumcand < 0.05) PFElecCut = true;
      }
      if(isClMatch_eS[jet] == 0)
    	if(eMax_F[jet]/Sumcand < 0.05 ) PFElecCut = true;

      if(isPFElecCut_eS[jet] != PFElecCut && printDebug) 
	cout<<" pf e cut not same, run = "<<run_F<<" "<<run_eS<<", event = "<<evt_F<<" "<<evt_eS<<" , lumi = "<<lumi_F<<" "<<lumi_eS<<endl;
	  
      if(!PFElecCut) continue;

      // also need to cut on the high pT jets from the Jet 55 sample. unmatched PF jets with pfpt > 110 and calopt < 30 including -999, check on their high track content. at the moment dont worry about this.  
    
      float recpt = pt_F[jet];
      
      if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0){
	//if(recpt > 140) continue;
	hpbpb_TrgObj55[cBin]->Fill(recpt, jet55_p_F* wght);
	hpbpb_raw_TrgObj55[cBin]->Fill(rawpt_F[jet], jet55_p_F* wght);
	hpbpb_anaBin_TrgObj55[cBin]->Fill(recpt, jet55_p_F* wght);
	hpbpb_JEC_TrgObj55[cBin]->Fill(recpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), jet55_p_F* wght);
	hpbpb_Smear_TrgObj55[cBin]->Fill(recpt + rnd.Gaus(0,1), jet55_p_F* wght);
      }

      if(jet65_F == 1 && jet80_F == 0){
	//if(recpt > 140) continue;
	hpbpb_TrgObj65[cBin]->Fill(recpt, wght);
	hpbpb_raw_TrgObj65[cBin]->Fill(rawpt_F[jet], wght);
	hpbpb_anaBin_TrgObj65[cBin]->Fill(recpt, wght);
	hpbpb_JEC_TrgObj65[cBin]->Fill(recpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj65[cBin]->Fill(recpt + rnd.Gaus(0,1), wght);
      }
      
      if(jet80_F == 1){
	hpbpb_TrgObj80[cBin]->Fill(recpt, wght);
	hpbpb_raw_TrgObj80[cBin]->Fill(rawpt_F[jet], wght);
	hpbpb_anaBin_TrgObj80[cBin]->Fill(recpt, wght);
	hpbpb_JEC_TrgObj80[cBin]->Fill(recpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj80[cBin]->Fill(recpt + rnd.Gaus(0,1), wght);
      }
      
    }// jet loop
    if(printDebug)cout<<endl;

  }// event loop

  for(int i = 0; i<nbins_cent; ++i){

    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj80[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj65[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj55[i]);
    
    divideBinWidth(hpbpb_TrgObjComb[i]);
    divideBinWidth(hpbpb_TrgObj80[i]);
    divideBinWidth(hpbpb_TrgObj65[i]);
    divideBinWidth(hpbpb_TrgObj55[i]);

    hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj80[i]);
    hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj65[i]);
    hpbpb_raw_TrgObjComb[i]->Add(hpbpb_raw_TrgObj55[i]);

    divideBinWidth(hpbpb_raw_TrgObjComb[i]);
    divideBinWidth(hpbpb_raw_TrgObj80[i]);
    divideBinWidth(hpbpb_raw_TrgObj65[i]);
    divideBinWidth(hpbpb_raw_TrgObj55[i]);

    hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj80[i]);
    hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj65[i]);
    hpbpb_anaBin_TrgObjComb[i]->Add(hpbpb_anaBin_TrgObj55[i]);

    divideBinWidth(hpbpb_anaBin_TrgObjComb[i]);
    divideBinWidth(hpbpb_anaBin_TrgObj80[i]);
    divideBinWidth(hpbpb_anaBin_TrgObj65[i]);
    divideBinWidth(hpbpb_anaBin_TrgObj55[i]);

    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj80[i]);
    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj65[i]);
    hpbpb_Smear_TrgObjComb[i]->Add(hpbpb_Smear_TrgObj55[i]);

    divideBinWidth(hpbpb_Smear_TrgObjComb[i]);

    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj80[i]);
    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj65[i]);
    hpbpb_JEC_TrgObjComb[i]->Add(hpbpb_JEC_TrgObj55[i]);

    divideBinWidth(hpbpb_JEC_TrgObjComb[i]);
    
  }

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
