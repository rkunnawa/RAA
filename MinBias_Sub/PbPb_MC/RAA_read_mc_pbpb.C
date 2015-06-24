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

#define pi 3.14159265

// static const int nbins_pt = 39;
// static const double boundaries_pt[nbins_pt+1] = {
//   3, 4, 5, 7, 9, 12, 
//   15, 18, 21, 24, 28,
//   32, 37, 43, 49, 56,
//   64, 74, 84, 97, 114,
//   133, 153, 174, 196,
//   220, 245, 272, 300, 
//   330, 362, 395, 430,
//   468, 507, 548, 592,
//   638, 686, 1000 
// };

static const int nbins_pt = 14;
static const double boundaries_pt[nbins_pt+1] = {
  49, 56, 64, 74, 84, 97, 114, 133,
  153, 174, 196, 220, 245, 272, 300
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

static const double pthat[10] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 2000};
static const double xsecs[10] = {2.034e-01,
				 1.075e-02,
				 1.025e-03,
				 9.865e-05,
				 1.129e-05,
				 1.465e-06,
				 2.837e-07,
				 5.323e-08,
				 5.934e-09,
				 0.0};

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

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
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

void RAA_read_mc_pbpb(int startfile = 5,
		      int endfile = 6,
		      int radius = 3,
		      std::string kFoname="PbPb_MC_histograms_FromForest_pthatCutonJets_akPu3_20_eta_20_6.root"){
  
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

  infile_Forest = "jetRAA_PbPb_mc_forest.txt";
  infile_Select = Form("jetRAA_PbPb_mc_akPu%d_select.txt",radius);
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

    jetpbpb[0]->Add(filename_Forest.c_str());
    jetpbpb[1]->Add(filename_Forest.c_str());
    jetpbpb[2]->Add(filename_Forest.c_str());
    jetpbpb[3]->Add(filename_Forest.c_str());
    jetpbpb[4]->Add(filename_Forest.c_str());
    jet_select->Add(filename_Select.c_str());
    evt_select->Add(filename_Select.c_str());

    cout<<"filename: "<<filename_Forest<<endl;
    cout<<"pthat = "<<pthat[startfile]<<endl;
    
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

    cout<<"Total number of events loaded in HiForest = "<<jetpbpb[2]->GetEntries()<<endl;

  }
  
  TF1 *fVz=0;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);
  
  
  jetpbpb[2]->AddFriend(jetpbpb[0]);
  jetpbpb[2]->AddFriend(jetpbpb[1]);
  jetpbpb[2]->AddFriend(jetpbpb[4]);
  jetpbpb[3]->AddFriend(jetpbpb[0]);
  jetpbpb[3]->AddFriend(jetpbpb[1]);
  jetpbpb[3]->AddFriend(jetpbpb[4]);

  // get the number of events per pthat bin from the HiForest store that value as an integer.
  double Nevents[9];
  
  for(int iter = 0; iter < 9; ++iter){
    xsecs[iter] = xsecs[iter] - xsecs[iter+1];
    TCut pthatCut(Form("pthat >= %f && pthat < %f",pthat[i],pthat[i+1]));
    Nevents[i] = jetpbpb[2]->GetEntries(pthatCut);
  }
  
  jetpbpb[2]->AddFriend(evt_select);
  jetpbpb[3]->AddFriend(evt_select);
  // Forest files 
  int nref_F;
  float pt_F[1000];
  float refpt_F[1000];
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
  int   subid_F[1000];
  float refdrjt_F[1000];
  float pthat_F;
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

  float calopt_F[1000];
  jetpbpb[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jetpbpb[4]->SetBranchAddress("evt",&evt_F);
  jetpbpb[4]->SetBranchAddress("run",&run_F);
  jetpbpb[4]->SetBranchAddress("lumi",&lumi_F);
  jetpbpb[4]->SetBranchAddress("hiBin",&hiBin_F);
  jetpbpb[4]->SetBranchAddress("vz",&vz_F);
  jetpbpb[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  //jetpbpb[0]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  //jetpbpb[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  //jetpbpb[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpbpb[2]->SetBranchAddress("pthat",&pthat_F);
  jetpbpb[2]->SetBranchAddress("nref",&nref_F);
  jetpbpb[2]->SetBranchAddress("subid",subid_F);
  jetpbpb[2]->SetBranchAddress("refdrjt",refdrjt_F);
  jetpbpb[2]->SetBranchAddress("refpt",refpt_F);
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
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v7",&jet55_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v7",&jet65_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v7",&jet80_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_p_F);
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
  int   isPFElecSel_jS[1000];

  jet_select->SetBranchAddress("hiBin",&hiBin_jS);
  jet_select->SetBranchAddress("run_value",&run_jS);
  jet_select->SetBranchAddress("evt_value",&evt_jS);
  jet_select->SetBranchAddress("lumi_value",&lumi_jS);
  jet_select->SetBranchAddress("vz",&vz_jS);
  jet_select->SetBranchAddress("npf", &nref_jS);  
  jet_select->SetBranchAddress("pfpt", pt_jS);  
  jet_select->SetBranchAddress("eMax", eMax_jS);  
  jet_select->SetBranchAddress("pfeta",eta_jS);  
  jet_select->SetBranchAddress("calopt",calopt_jS);  
  jet_select->SetBranchAddress("isCaloMatch",isCaloMatch_jS);  
  jet_select->SetBranchAddress("isMultiMatch",isMultiMatch_jS);  
  jet_select->SetBranchAddress("isPFElecSel",isPFElecSel_jS);  
  jet_select->SetBranchAddress("chSum",chSum_jS);  
  jet_select->SetBranchAddress("phSum",phSum_jS);  
  jet_select->SetBranchAddress("neSum",neSum_jS);  
  jet_select->SetBranchAddress("muSum",muSum_jS);  

  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

    //get the spectra with the specific trigger object from the different files. 
  TH1F *hpbpb_Jet80_gen[nbins_cent],
    *hpbpb_Jet80_reco[nbins_cent],
    *hpbpb_Jet80_raw[nbins_cent],
    *hpbpb_Jet80_GenSmear[nbins_cent],
    *hpbpb_Jet80_RecoSmear[nbins_cent];
  TH1F *hpbpb_Jet65_gen[nbins_cent],
    *hpbpb_Jet65_reco[nbins_cent],
    *hpbpb_Jet65_raw[nbins_cent],
    *hpbpb_Jet65_GenSmear[nbins_cent],
    *hpbpb_Jet65_RecoSmear[nbins_cent];
  TH1F *hpbpb_Jet55_gen[nbins_cent],
    *hpbpb_Jet55_reco[nbins_cent],
    *hpbpb_Jet55_raw[nbins_cent],
    *hpbpb_Jet55_GenSmear[nbins_cent],
    *hpbpb_Jet55_RecoSmear[nbins_cent],
    * hpbpb_Jet55_reco_5TrigIneff_Smear[nbins_cent],
    * hpbpb_Jet55_reco_10TrigIneff_Smear[nbins_cent],
    * hpbpb_Jet55_gen_5TrigIneff_Smear[nbins_cent],
    * hpbpb_Jet55_gen_10TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_JetComb_gen[nbins_cent],
    *hpbpb_JetComb_reco[nbins_cent],
    *hpbpb_JetComb_raw[nbins_cent],
    *hpbpb_JetComb_GenSmear[nbins_cent],
    *hpbpb_JetComb_RecoSmear[nbins_cent],
    * hpbpb_JetComb_reco_5TrigIneff_Smear[nbins_cent],
    * hpbpb_JetComb_reco_10TrigIneff_Smear[nbins_cent],
    * hpbpb_JetComb_gen_5TrigIneff_Smear[nbins_cent],
    * hpbpb_JetComb_gen_10TrigIneff_Smear[nbins_cent];
  
  TH1F * hpbpb_JetComb_gen2pSmear[nbins_cent],
    * hpbpb_Jet80_gen2pSmear[nbins_cent],
    * hpbpb_Jet65_gen2pSmear[nbins_cent],
    * hpbpb_Jet55_gen2pSmear[nbins_cent];

  TH1F * hpbpb_pthat_Cut[nbins_cent];
  TH1F * hpbpb_pthat_noCut[nbins_cent];

  TH1F *hpbpb_gen[nbins_cent],*hpbpb_reco[nbins_cent];
  TH2F *hpbpb_matrix[nbins_cent];
  TH2F *hpbpb_matrix_HLT[nbins_cent];
  TH2F *hpbpb_Trans_matrix_HLT[nbins_cent];
  TH2F *hpbpb_matrix_HLT_GenSmear[nbins_cent];
  TH2F *hpbpb_matrix_HLT_gen2pSmear[nbins_cent];
  TH2F *hpbpb_matrix_HLT_RecoSmear[nbins_cent];
  TH2F *hpbpb_matrix_HLT_BothSmear[nbins_cent];
  TH2F *hpbpb_matrix_HLT_5TrigIneff_Smear[nbins_cent];
  TH2F *hpbpb_matrix_HLT_10TrigIneff_Smear[nbins_cent];
  
  TH2F *hpbpb_mcclosure_matrix[nbins_cent];
  TH2F *hpbpb_mcclosure_matrix_HLT[nbins_cent];
  TH2F *hpbpb_mcclosure_Trans_matrix_HLT[nbins_cent];
  TH2F *hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[nbins_cent];  
  TH2F *hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[nbins_cent];  
  TH1F *hpbpb_mcclosure_data[nbins_cent];
  TH1F *hpbpb_mcclosure_data_train[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_data[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet80_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet65_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_data_train[nbins_cent];
  // TH1F *hpbpb_mcclosure_JetComb_data_train_5TrigIneff_Smear[nbins_cent];
  // TH1F *hpbpb_mcclosure_JetComb_data_train_10TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet80_data_train[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet65_data_train[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_data_train[nbins_cent];
  // TH1F *hpbpb_mcclosure_Jet55_data_train_5TrigIneff_Smear[nbins_cent];
  // TH1F *hpbpb_mcclosure_Jet55_data_train_10TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet80_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet65_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[nbins_cent];

  TH2F *hpbpb_anaBin_matrix_HLT[nbins_cent];
  TH2F *hpbpb_anaBin_Trans_matrix_HLT[nbins_cent];
  TH1F *hpbpb_anaBin_Jet80_gen[nbins_cent],*hpbpb_anaBin_Jet80_reco[nbins_cent];
  TH1F *hpbpb_anaBin_Jet65_gen[nbins_cent],*hpbpb_anaBin_Jet65_reco[nbins_cent];
  TH1F *hpbpb_anaBin_Jet55_gen[nbins_cent],*hpbpb_anaBin_Jet55_reco[nbins_cent];
  TH1F *hpbpb_anaBin_JetComb_gen[nbins_cent],*hpbpb_anaBin_JetComb_reco[nbins_cent];

  
  for(int i = 0;i<nbins_cent;++i){

    hpbpb_pthat_Cut[i] = new TH1F(Form("pthat_Cut_cent%d",i),"",1000,0,1000);
    hpbpb_pthat_noCut[i] = new TH1F(Form("pthat_noCut_cent%d",i),"",1000,0,1000);
    
    hpbpb_gen[i] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    //cout<<"A"<<endl;
    hpbpb_reco[i] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("Reco jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    //cout<<"B"<<endl;
    hpbpb_matrix[i] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    
    hpbpb_matrix_HLT[i] = new TH2F(Form("hpbpb_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    hpbpb_Trans_matrix_HLT[i] = new TH2F(Form("hpbpb_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    hpbpb_matrix_HLT_GenSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_GenSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    hpbpb_matrix_HLT_gen2pSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_gen2pSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    hpbpb_matrix_HLT_RecoSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_RecoSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    hpbpb_matrix_HLT_BothSmear[i] = new TH2F(Form("hpbpb_matrix_HLT_BothSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    hpbpb_matrix_HLT_5TrigIneff_Smear[i] = new TH2F(Form("hpbpb_matrix_HLT_5TrigIneffSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);
    hpbpb_matrix_HLT_10TrigIneff_Smear[i] = new TH2F(Form("hpbpb_matrix_HLT_10TrigIneffSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000,1000,0,1000);

    hpbpb_anaBin_matrix_HLT[i] = new TH2F(Form("hpbpb_anaBin_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    hpbpb_anaBin_Trans_matrix_HLT[i] = new TH2F(Form("hpbpb_anaBin_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    

    hpbpb_mcclosure_matrix[i] = new TH2F(Form("hpbpb_mcclosure_matrix_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    hpbpb_mcclosure_matrix_HLT[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Trans_matrix_HLT[i] = new TH2F(Form("hpbpb_mcclosure_Trans_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);
    //cout<<"C"<<endl;
    hpbpb_mcclosure_data[i] =new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_data_train[i] =new TH1F(Form("hpbpb_mcclosure_data_train_R%d_%s_cent%d",radius,etaWidth,i),Form("data_train for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_JetComb_data[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet80_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 80  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet65_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 65  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet55_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);

    hpbpb_mcclosure_JetComb_data_train[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_train_R%d_%s_cent%d",radius,etaWidth,i),Form("data_train for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet80_data_train[i] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_train_R%d_%s_cent%d",radius,etaWidth,i),Form("data_train for unfolding mc closure test trigger 80  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet65_data_train[i] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_train_R%d_%s_cent%d",radius,etaWidth,i),Form("data_train for unfolding mc closure test trigger 65  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet55_data_train[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_train_R%d_%s_cent%d",radius,etaWidth,i),Form("data_train for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);

    hpbpb_mcclosure_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_JetComb_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet80_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet80_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet65_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet65_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 65 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet55_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_10TrigIneff_Smear[_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    
    hpbpb_JetComb_gen[i] = new TH1F(Form("hpbpb_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JetComb_reco[i] = new TH1F(Form("hpbpb_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet80_gen[i] = new TH1F(Form("hpbpb_Jet80_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet80_reco[i] = new TH1F(Form("hpbpb_Jet80_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet65_gen[i] = new TH1F(Form("hpbpb_Jet65_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet65_reco[i] = new TH1F(Form("hpbpb_Jet65_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_gen[i] = new TH1F(Form("hpbpb_Jet55_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_reco[i] = new TH1F(Form("hpbpb_Jet55_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
						
    hpbpb_JetComb_GenSmear[i] = new TH1F(Form("hpbpb_JetComb_GenSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JetComb_gen2pSmear[i] = new TH1F(Form("hpbpb_JetComb_gen2pSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("gen2pSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JetComb_RecoSmear[i] = new TH1F(Form("hpbpb_JetComb_RecoSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("RecoSmear jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet80_GenSmear[i] = new TH1F(Form("hpbpb_Jet80_GenSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet80_gen2pSmear[i] = new TH1F(Form("hpbpb_Jet80_gen2pSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("gen2pSmear refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet80_RecoSmear[i] = new TH1F(Form("hpbpb_Jet80_RecoSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("RecoSmear jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet65_GenSmear[i] = new TH1F(Form("hpbpb_Jet65_GenSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet65_gen2pSmear[i] = new TH1F(Form("hpbpb_Jet65_gen2pSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("gen2pSmear refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet65_RecoSmear[i] = new TH1F(Form("hpbpb_Jet65_RecoSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("RecoSmear jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_GenSmear[i] = new TH1F(Form("hpbpb_Jet55_GenSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_gen2pSmear[i] = new TH1F(Form("hpbpb_Jet55_gen2pSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("gen2pSmear refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_RecoSmear[i] = new TH1F(Form("hpbpb_Jet55_RecoSmear_R%d_%s_cent%d",radius,etaWidth,i),Form("RecoSmear jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_JetComb_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_gen_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JetComb_reco_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_reco_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_gen_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_gen_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_reco_5TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_reco_5TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_JetComb_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_gen_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_JetComb_reco_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_JetComb_reco_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_gen_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_gen_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_reco_10TrigIneff_Smear[i] = new TH1F(Form("hpbpb_Jet55_reco_10TrigIneff_Smear_R%d_%s_cent%d",radius,etaWidth,i),Form("GenSmear refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_anaBin_JetComb_gen[i] = new TH1F(Form("hpbpb_anaBin_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_JetComb_reco[i] = new TH1F(Form("hpbpb_anaBin_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet80_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet80_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet80_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet80_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet65_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet65_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet65_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet65_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet55_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet55_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet55_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet55_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);

    hpbpb_JetComb_raw[i] = new TH1F(Form("hpbpb_JetComb_raw_R%d_%s_cent%d",radius,etaWidth,i),Form("raw jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet80_raw[i] = new TH1F(Form("hpbpb_Jet80_raw_R%d_%s_cent%d",radius,etaWidth,i),Form("raw jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet65_raw[i] = new TH1F(Form("hpbpb_Jet65_raw_R%d_%s_cent%d",radius,etaWidth,i),Form("raw jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_Jet55_raw[i] = new TH1F(Form("hpbpb_Jet55_raw_R%d_%s_cent%d",radius,etaWidth,i),Form("raw jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
  }
  
  // now start the event loop for each file. 
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpbpb[0]->GetEntries();
  Long64_t nGoodEvt = 0;
  //if(printDebug) nentries = 10;
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

    //if(printDebug) cout<<"forest values = "<<hiBin_F<<", "<<evt_F<<", "<<run_F<<", "<<lumi_F<<", "<<vz_F<<endl;

    // if(pcollisionEventSelection_F==0) continue; 
    // if(fabs(vz_F)>15) continue;
    if(!isGoodEvent_eS) continue; 
    
    jet_select->GetEntry(nGoodEvt);
    ++nGoodEvt;
 
    int cBin = findBin(hiBin_F);//tells us the centrality of the event. 
    if(cBin==-1 || cBin==nbins_cent) continue;

    if(hiBin_eS != hiBin_F || evt_eS != evt_F || run_eS != run_F || lumi_eS != lumi_F || vz_eS != vz_F) cout<<"ERROR mismatch eS, F"<<endl;
    if(hiBin_eS != hiBin_jS || evt_eS != evt_jS || run_eS != run_jS || lumi_eS != lumi_jS || vz_eS != vz_jS) cout<<"ERROR mismatch eS, jS"<<endl;
    if(hiBin_F != hiBin_jS || evt_F != evt_jS || run_F != run_jS || lumi_F != lumi_jS || vz_F != vz_jS) cout<<"ERROR mismatch F, jS"<<endl;

    
    if(printDebug) cout<<"hibin hiForest = "<<hiBin_F<<", evtTree = "<<hiBin_eS<<", jetTree = "<<hiBin_jS<<endl;
    if(printDebug) cout<<"evt hiForest   = "<<evt_F<<", evtTree = "<<evt_eS<<", jetTree = "<<evt_jS<<endl;
    if(printDebug) cout<<"lumi hiForest  = "<<lumi_F<<", evtTree = "<<lumi_eS<<", jetTree = "<<lumi_jS<<endl;
    if(printDebug) cout<<"vz hiForest    = "<<vz_F<<", evtTree = "<<vz_eS<<", jetTree = "<<vz_jS<<endl;
    
    if(printDebug) cout<<"nref_F = "<<nref_F<<", nref_eS = "<<nref_eS<<", nref_jS = "<<nref_jS<<endl;

    if(nref_F != nref_eS) cout<<"ERROR mismatch in jet counts"<<endl;

    hpbpb_pthat_noCut[cBin]->Fill(pthat_F, weight_eS);

    // // if(pthat_F < pthat[startfile] || pthat_F >= pthat[endfile]) continue;
    // double cweight = exp(- pow(varbin+1.11957e+01,2) / pow(1.34120e+01,2) / 2);

    // double weight = 1;
    // weight = cweight * fVz->Eval(vz_F);    
    // double scale = xsecs[]
    
    // hpbpb_pthat_Cut[cBin]->Fill(pthat_F, weight);
    
    /*    //! Sort the jetTree jets according to pT, from the evtTree now. 
    std::vector < Jet > vJet;
    for(int jet2 = 0; jet2<nref_jS; ++jet2){
      if(printDebug) cout <<"\t \t jetTree *** "<< jet2 <<  ", pT " << pt_jS[jet2] <<", eta : "<< eta_jS[jet2] <<  ", chSum : "<< chSum_jS[jet2]  << endl;
      //if(printDebug) cout <<"\t \t evtTree *** "<< jet2 <<  ", pT " << pt_eS[jet2]<< endl;
      Jet ijet;
      ijet.id = jet2;
      ijet.pt = pt_jS[jet2];
      vJet.push_back(ijet);
    }
    std::sort (vJet.begin(), vJet.end(), compare_pt);
    std::vector < Jet >::const_iterator itr;

    // after sorting the jet array according to pT
    if(printDebug){
      cout<<"**** AFTER  SORTING ****"<<endl;
      for(itr=vJet.begin(); itr!=vJet.end(); ++itr){
	int jetLoc = (*itr).id;
	cout <<"\t \t jetTree *** "<< jetLoc <<  ", pT " << pt_jS[jetLoc]<< endl;
      }
    }
    int jet=0;
    if(printDebug) cout<<"Number of jets = "<< vJet.size()<<endl;
    for(itr=vJet.begin(); itr!=vJet.end(); ++itr, ++jet){

      bool PFElecCut = false;
      int jetLoc = (*itr).id;
      //if(isMiMatch_eS[jetLoc]) ++itr;
      if(isMultiMatch_jS[jetLoc]) {
	++itr;
	jetLoc = (*itr).id;
	if(itr == vJet.end())  break;
      }
      Float_t Sumcand = chSum_jS[jetLoc] + phSum_jS[jetLoc] + neSum_jS[jetLoc] + muSum_jS[jetLoc];
      
      if(fabs(eta_jS[jetLoc]) > 2) continue;
      // //if(isPFElecCut_eS[jet] != 1) continue;
      if(pt_jS[jetLoc] > 2 * pthat_F) continue;
      int refid = -1;
      refid = pfrefidx_jS[jetLoc];
      //if(refid < 0) continue;
      
      //if(printDebug) cout<<"pfrefidx_jS[jetLoc] = "<<pfrefidx_jS[jetLoc]<<", refdrjt_jS[refid] = "<<refdrjt_jS[refid]<<endl;

      Float_t delR = deltaR(eta_F[jet], phi_F[jet], refeta_F[jet], refphi_F[jet]);
      if(subid_F[jet] != 0 || refpt_F[jet] < 0 || delR > 0.2 || refpt_F[jet] > 2 * pthat_F) continue;

      if(pt_jS[jetLoc] != pt_F[jet]) {cout<<"Jets NOT Matched!!! ERROR!!!"<<endl; break;}
      
      //cout<<"refid = "<<refid<<", jetLoc = "<<jetLoc<<endl;

      if(isCaloMatch_jS[jetLoc] == 1){
	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.5 && calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.85 && eMax_jS[jetLoc]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_jS[jetLoc]/pt_jS[jetLoc] - (Float_t)9/7)) PFElecCut = true;
	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.85) PFElecCut = true;
	if(calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.5 && eMax_jS[jetLoc]/Sumcand < 0.05) PFElecCut = true;
      }else if(isCaloMatch_jS[jetLoc] == 0)
	if(eMax_jS[jetLoc]/Sumcand < 0.05 ) PFElecCut = true;       

      if(!PFElecCut) continue;

      if(printDebug) cout<<"jet = "<<jet<<", pt_F[jet] = "<<pt_F[jet]<<", jetLoc = "<<jetLoc<<", pt_jS[jetLoc] = "<<pt_jS[jetLoc]<<", calopt_jS[jetLoc] = "<<calopt_jS[jetLoc]<<", eMax/sumcand = "<<(Float_t)eMax_jS[jetLoc]/Sumcand<<", isCaloMatch_jS[jetLoc] = "<<isCaloMatch_jS[jetLoc]<<", PFEleccut = "<<PFElecCut<<", subid = "<<subid_jS[jetLoc]<<", refpt_F[jet] = "<< refpt_F[jet] <<", refid = "<< refid <<", " <<endl;

      // if(printDebug && (fabs(eta_jS[jet] > 2))) cout<<"jets with |eta| > 2 in jetTree"<<endl;
      // if(printDebug && (fabs(eta_F[jet] > 2)))  cout<<"jets with |eta| > 2 in Forest"<<endl;

      //if(printDebug && index_eS[jet] >= 0 )cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT [indes_eS[jet]] = "<<calopt_F[index_eS[jet]]<<", Calo pT[jetLoc] "<<calopt_F[jetLoc]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl; 
      //if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT [index_eS[jet]] = "<<calopt_F[index_eS[jet]]<<", Calo pT[jetLoc] "<<calopt_F[jetLoc]<< ", calopt_jS[jetLoc] = "<<calopt_jS[jetLoc]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl; 

      // if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;
      */
    for( int jet = 0; jet<nref_F; jet++ ){

      //cout  << "  subid : "  << subid_F[jet] << " refpt : "  << refpt_F[jet] << " pfpt :  " << pt_F[jet] << endl;

      if( fabs(eta_F[jet]) > 2.0 ) continue;
      if( subid_F[jet] != 0 || refpt_F[jet] <= 0) continue;
      if( isPFElecCut_eS[jet] != 1 ) continue;
      // if( (fabs(refpt_F[jet]) > 3.0 * pthat_F) || (pt_F[jet]/refpt_F[jet] > 2.0) ) continue;
      //if( refpt_F[jet] > 2.0 * pthat[startfile] || pt_F[jet] > 2.0 * pthat[startfile]) continue;
      
      Float_t genpt = refpt_F[jet];
      Float_t recpt = pt_F[jet];
      Float_t rawpt = rawpt_F[jet];
      //if(genpt < 0) continue;
      
      hpbpb_gen[cBin]->Fill(genpt, weight_eS);
      hpbpb_reco[cBin]->Fill(recpt, weight_eS);
      hpbpb_matrix[cBin]->Fill(genpt, recpt, weight_eS);
      
      if(nEvt%2==0){
	hpbpb_mcclosure_data[cBin]->Fill(recpt, weight_eS);
      }else {
	hpbpb_mcclosure_matrix[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_mcclosure_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_mcclosure_data_train[cBin]->Fill(recpt, weight_eS);
      }

      if(jet55_F == 1 && jet65_F==0 && jet80_F == 0){
	hpbpb_Jet55_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_Jet55_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet55_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_Jet55_gen_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_Jet55_gen_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.10/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_Jet55_reco_5TrigIneff_Smear[cBin]->Fill(recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_Jet55_reco_10TrigIneff_Smear[cBin]->Fill(recpt * (1. + 0.10/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);

	hpbpb_matrix_HLT_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS);

	hpbpb_Jet55_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_Jet55_RecoSmear[cBin]->Fill(recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet55_raw[cBin]->Fill(rawpt, weight_eS);
	hpbpb_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_Trans_matrix_HLT[cBin]->Fill( recpt, genpt, weight_eS);
	hpbpb_matrix_HLT_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS);
	hpbpb_matrix_HLT_RecoSmear[cBin]->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_matrix_HLT_BothSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS);
	
	hpbpb_anaBin_Jet55_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_anaBin_Jet55_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_anaBin_Trans_matrix_HLT[cBin]->Fill( recpt,genpt, weight_eS);

	if(nEvt%2==0){
	  hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[cBin]->Fill(recpt* (1. + 0.5/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[cBin]->Fill(recpt* (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_Jet55_data[cBin]->Fill(recpt, weight_eS);
	}else {
	  hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	  hpbpb_mcclosure_Trans_matrix_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[cBin]->Fill(genpt* (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_Jet55_gen[cBin]->Fill(genpt, weight_eS);
	  hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[cBin]->Fill(genpt* (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_Jet55_data_train[cBin]->Fill(recpt, weight_eS);
	}
	
	
      }// jet55 selection

      if(jet65_F == 1 && jet80_F == 0){
	hpbpb_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_Jet65_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet65_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	
	hpbpb_matrix_HLT_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS);

	hpbpb_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_Jet65_RecoSmear[cBin]->Fill(recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet65_raw[cBin]->Fill(rawpt, weight_eS);
	hpbpb_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_Trans_matrix_HLT[cBin]->Fill( recpt, genpt, weight_eS);
	hpbpb_matrix_HLT_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS);
	hpbpb_matrix_HLT_RecoSmear[cBin]->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_matrix_HLT_BothSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS);
	
	hpbpb_anaBin_Jet65_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_anaBin_Jet65_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_anaBin_Trans_matrix_HLT[cBin]->Fill( recpt,genpt, weight_eS);

	if(nEvt%2==0){
	  hpbpb_mcclosure_Jet65_data[cBin]->Fill(recpt, weight_eS);
	}else {
	  hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_Trans_matrix_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  hpbpb_mcclosure_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	  hpbpb_mcclosure_Jet65_data_train[cBin]->Fill(recpt, weight_eS);
	  hpbpb_mcclosure_Jet65_gen[cBin]->Fill(genpt, weight_eS);

	}
	
      }// jet65 selection

      if(jet80_F == 1){
	hpbpb_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_Jet80_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet80_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	hpbpb_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);


	hpbpb_matrix_HLT_gen2pSmear[cBin]->Fill(genpt * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), recpt, weight_eS);

	hpbpb_Jet80_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_Jet80_RecoSmear[cBin]->Fill(recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_Jet80_raw[cBin]->Fill(rawpt, weight_eS);
	hpbpb_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_Trans_matrix_HLT[cBin]->Fill( recpt, genpt, weight_eS);
	hpbpb_matrix_HLT_GenSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt, weight_eS);
	hpbpb_matrix_HLT_RecoSmear[cBin]->Fill(genpt, recpt + rnd.Gaus(0,1), weight_eS);
	hpbpb_matrix_HLT_BothSmear[cBin]->Fill(genpt + rnd.Gaus(0,1), recpt + rnd.Gaus(0,1), weight_eS);
	
	hpbpb_anaBin_Jet80_gen[cBin]->Fill(genpt, weight_eS);
	hpbpb_anaBin_Jet80_reco[cBin]->Fill(recpt, weight_eS);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	hpbpb_anaBin_Trans_matrix_HLT[cBin]->Fill( recpt,genpt, weight_eS);

	if(nEvt%2==0){
	  hpbpb_mcclosure_Jet80_data[cBin]->Fill(recpt, weight_eS);
	}else {
	  hpbpb_mcclosure_matrix_HLT_5TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.05/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_matrix_HLT_10TrigIneff_Smear[cBin]->Fill(genpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), recpt * (1. + 0.1/nbins_cent*(nbins_cent-cBin)), weight_eS);
	  hpbpb_mcclosure_Trans_matrix_HLT[cBin]->Fill(recpt, genpt, weight_eS);
	  hpbpb_mcclosure_matrix_HLT[cBin]->Fill(genpt, recpt, weight_eS);
	  hpbpb_mcclosure_Jet80_data_train[cBin]->Fill(recpt, weight_eS);
	  hpbpb_mcclosure_Jet80_gen[cBin]->Fill(genpt, weight_eS);

	}
	
      }// jet80 selection

      
    }// jet loop
    if(printDebug)cout<<endl;

  }// event loop

  for(int i = 0; i<nbins_cent; ++i){
    
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet55_data[i]);
    
    divideBinWidth(hpbpb_mcclosure_JetComb_data[i]);
        

    hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_data_5TrigIneff_Smear[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_data_5TrigIneff_Smear[i]);
    

    hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_data_10TrigIneff_Smear[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_data_10TrigIneff_Smear[i]);
    
    
    hpbpb_mcclosure_JetComb_data_train[i]->Add(hpbpb_mcclosure_Jet80_data_train[i]);
    hpbpb_mcclosure_JetComb_data_train[i]->Add(hpbpb_mcclosure_Jet65_data_train[i]);
    hpbpb_mcclosure_JetComb_data_train[i]->Add(hpbpb_mcclosure_Jet55_data_train[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_data_train[i]);
    
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);
    
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet55_gen[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_gen[i]);

    hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);
    hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    
    hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_gen_10TrigIneff_Smear[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_gen_10TrigIneff_Smear[i]);

    hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);
    hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_mcclosure_Jet55_gen_5TrigIneff_Smear[i]);
    

    divideBinWidth(hpbpb_mcclosure_JetComb_gen_5TrigIneff_Smear[i]);
    
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet80_reco[i]);
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet65_reco[i]);
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet55_reco[i]);
    

    hpbpb_JetComb_reco_5TrigIneff_Smear[i]->Add(hpbpb_Jet80_reco[i]);
    hpbpb_JetComb_reco_5TrigIneff_Smear[i]->Add(hpbpb_Jet65_reco[i]);
    hpbpb_JetComb_reco_5TrigIneff_Smear[i]->Add(hpbpb_Jet55_reco_5TrigIneff_Smear[i]);

    
    hpbpb_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_Jet80_gen[i]);
    hpbpb_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_Jet65_gen[i]);
    hpbpb_JetComb_gen_5TrigIneff_Smear[i]->Add(hpbpb_Jet55_gen_5TrigIneff_Smear[i]);

    
    hpbpb_JetComb_reco_10TrigIneff_Smear[i]->Add(hpbpb_Jet80_reco[i]);
    hpbpb_JetComb_reco_10TrigIneff_Smear[i]->Add(hpbpb_Jet65_reco[i]);
    hpbpb_JetComb_reco_10TrigIneff_Smear[i]->Add(hpbpb_Jet55_reco_10TrigIneff_Smear[i]);

    
    hpbpb_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_Jet80_gen[i]);
    hpbpb_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_Jet65_gen[i]);
    hpbpb_JetComb_gen_10TrigIneff_Smear[i]->Add(hpbpb_Jet55_gen_10TrigIneff_Smear[i]);

    hpbpb_JetComb_RecoSmear[i]->Add(hpbpb_Jet80_RecoSmear[i]);
    
    hpbpb_JetComb_RecoSmear[i]->Add(hpbpb_Jet65_RecoSmear[i]);
    hpbpb_JetComb_RecoSmear[i]->Add(hpbpb_Jet55_RecoSmear[i]);

    hpbpb_JetComb_GenSmear[i]->Add(hpbpb_Jet80_GenSmear[i]);
    hpbpb_JetComb_GenSmear[i]->Add(hpbpb_Jet65_GenSmear[i]);
    
    hpbpb_JetComb_GenSmear[i]->Add(hpbpb_Jet55_GenSmear[i]);

    hpbpb_JetComb_gen2pSmear[i]->Add(hpbpb_Jet80_gen2pSmear[i]);
    hpbpb_JetComb_gen2pSmear[i]->Add(hpbpb_Jet65_gen2pSmear[i]);
    hpbpb_JetComb_gen2pSmear[i]->Add(hpbpb_Jet55_gen2pSmear[i]);
    
    
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet80_gen[i]);
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet65_gen[i]);
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet55_gen[i]);

    
    divideBinWidth(hpbpb_JetComb_gen[i]);
    divideBinWidth(hpbpb_JetComb_reco[i]);

    divideBinWidth(hpbpb_JetComb_GenSmear[i]);
    
    divideBinWidth(hpbpb_JetComb_gen2pSmear[i]);
    divideBinWidth(hpbpb_JetComb_RecoSmear[i]);
    
    divideBinWidth(hpbpb_JetComb_raw[i]);
    
    divideBinWidth(hpbpb_reco[i]);
    divideBinWidth(hpbpb_gen[i]);

    hpbpb_anaBin_JetComb_reco[i]->Add(hpbpb_anaBin_Jet80_reco[i]);
    hpbpb_anaBin_JetComb_reco[i]->Add(hpbpb_anaBin_Jet65_reco[i]);
    
    hpbpb_anaBin_JetComb_reco[i]->Add(hpbpb_anaBin_Jet55_reco[i]);
    
    hpbpb_anaBin_JetComb_gen[i]->Add(hpbpb_anaBin_Jet80_gen[i]);
    hpbpb_anaBin_JetComb_gen[i]->Add(hpbpb_anaBin_Jet65_gen[i]);
    hpbpb_anaBin_JetComb_gen[i]->Add(hpbpb_anaBin_Jet55_gen[i]);
    

    divideBinWidth(hpbpb_anaBin_JetComb_gen[i]);
    divideBinWidth(hpbpb_anaBin_JetComb_reco[i]);

  }

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
