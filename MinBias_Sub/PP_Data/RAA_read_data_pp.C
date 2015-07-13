// Raghav Kunnawalkam Elayavalli
// June 5th 2014
// CERN
// for questions or comments: raghav.k.e at CERN dot CH

// 
// read all the MC files for Pp and pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 
// 

// July 19 - all pp histograms will have 2D arrays with [radius][eta_bin]. the Pp histograms will be defined by 3D arrays with [radius][eta_bin][centrality]. 
// July 20 - the loop structure(s) are defined as follows for the several histograms 
//            Radius    : iteration variable: k;       number of iterations: no_radius;                                 values for the radii: list_radius
//            Eta bins  : iteration variable: j;       number of iterations: nbins_eta;                                 values of the bins  : boundaries_eta
//            Centrality: iteration variable: i;       number of iterations: nbins_cent +1;                             values of the bins  : boundaries_cent (till nbins_cent) + 0-200 (the whole range) for the final iteration. 
//            p_T Hats  : iteration variable: h;       number of iterations: nbins_pthat (Pp) and nbinsPP_pthat (pp); values of the bins  : boundaries_pthat (Pp) and boundariesPP_pthat (pp)  
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
static const char trigName [trigValue][256] = {"HLT40","HLT60","HLT80","Combined"};
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

using namespace std;

void RAA_read_data_pp(int startfile = 8,
		      int endfile = 9,
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

  infile_Forest = "jetRAA_pp_data_forest.txt";
  infile_Select = Form("jetRAA_pp_data_ak%d_select.txt",radius);
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

  TChain * jetpp[N];
  TChain * evt_select, * jet_select; 

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%dPFJetAnalyzer",radius);
  dir[3] = "ak3CaloJetAnalyzer";
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
    jetpp[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  evt_select = new TChain(Form("ak%dJetAnalyzer/evtTree",radius));
  jet_select = new TChain(Form("ak%dJetAnalyzer/jetTree",radius));

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;
    instr_Select>>filename_Select;

    jetpp[0]->Add(filename_Forest.c_str());
    jetpp[1]->Add(filename_Forest.c_str());
    jetpp[2]->Add(filename_Forest.c_str());
    jetpp[3]->Add(filename_Forest.c_str());
    jetpp[4]->Add(filename_Forest.c_str());
    jet_select->Add(filename_Select.c_str());
    evt_select->Add(filename_Select.c_str());
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[3]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpp[4]->GetEntries() << endl;
    if(printDebug)cout << "jet selection file " << jet_select->GetEntries() << endl;
    if(printDebug)cout << "event selection file" << evt_select->GetEntries() << endl;

  }
  
  jetpp[2]->AddFriend(jetpp[0]);
  jetpp[2]->AddFriend(jetpp[1]);
  jetpp[2]->AddFriend(jetpp[4]);
  jetpp[3]->AddFriend(jetpp[0]);
  jetpp[3]->AddFriend(jetpp[1]);
  jetpp[3]->AddFriend(jetpp[4]);
  
  jetpp[2]->AddFriend(evt_select);
  jetpp[3]->AddFriend(evt_select);

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
  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  float vz_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int pcollisionEventSelection_F;

  float calopt_F[1000];
  jetpp[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jetpp[4]->SetBranchAddress("evt",&evt_F);
  jetpp[4]->SetBranchAddress("run",&run_F);
  jetpp[4]->SetBranchAddress("lumi",&lumi_F);
  jetpp[4]->SetBranchAddress("vz",&vz_F);
  jetpp[1]->SetBranchAddress("pPAcollisionEventSelectionPA",&pcollisionEventSelection_F);
  //jetpp[0]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  //jetpp[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  //jetpp[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpp[2]->SetBranchAddress("nref",&nref_F);
  jetpp[2]->SetBranchAddress("jtpt",&pt_F);
  jetpp[2]->SetBranchAddress("jteta",&eta_F);
  jetpp[2]->SetBranchAddress("jtphi",&phi_F);
  jetpp[2]->SetBranchAddress("rawpt",&rawpt_F);
  jetpp[2]->SetBranchAddress("jtpu",&jtpu_F);
  jetpp[2]->SetBranchAddress("chargedMax",&chMax_F);
  jetpp[2]->SetBranchAddress("chargedSum",&chSum_F);
  jetpp[2]->SetBranchAddress("trackMax",&trkMax_F);
  jetpp[2]->SetBranchAddress("trackSum",&trkSum_F);
  jetpp[2]->SetBranchAddress("photonMax",&phMax_F);
  jetpp[2]->SetBranchAddress("photonSum",&phSum_F);
  jetpp[2]->SetBranchAddress("neutralMax",&neMax_F);
  jetpp[2]->SetBranchAddress("neutralSum",&neSum_F);
  jetpp[2]->SetBranchAddress("eSum",&eSum_F);
  jetpp[2]->SetBranchAddress("eMax",&eMax_F);
  jetpp[2]->SetBranchAddress("muSum",&muSum_F);
  jetpp[2]->SetBranchAddress("muMax",&muMax_F);
  jetpp[0]->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40_F);
  jetpp[0]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p_F);
  jetpp[0]->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60_F);
  jetpp[0]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p_F);
  jetpp[0]->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80_F);
  jetpp[0]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p_F);

  
  // event tree selection file:
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
  Double_t weight_eS;

  evt_select->SetBranchAddress("run_value",&run_eS);
  evt_select->SetBranchAddress("evt_value",&evt_eS);
  evt_select->SetBranchAddress("lumi_value",&lumi_eS);
  evt_select->SetBranchAddress("vz",&vz_eS);
  evt_select->SetBranchAddress("nref",&nref_eS);
  evt_select->SetBranchAddress("index", &index_eS);
  evt_select->SetBranchAddress("isGoodEvt",&isGoodEvent_eS);
  evt_select->SetBranchAddress("isMlMatch",&isMiMatch_eS);
  evt_select->SetBranchAddress("isClMatch", &isClMatch_eS);
  evt_select->SetBranchAddress("isPFElecCut",&isPFElecCut_eS);
  evt_select->SetBranchAddress("isTrackCut",&isTrackCut_eS);
  evt_select->SetBranchAddress("isMuCut",&isMuCut_eS);
  evt_select->SetBranchAddress("weight", &weight_eS);  


  // // jet tree selection file:
  // int run_jS;
  // int evt_jS;
  // int lumi_jS;
  // float vz_jS;
  // int nref_jS;
  // float pt_jS[1000];
  // float eta_jS[1000];
  // float eMax_jS[1000];

  // jet_select->SetBranchAddress("run_value",&run_jS);
  // jet_select->SetBranchAddress("evt_value",&evt_jS);
  // jet_select->SetBranchAddress("lumi_value",&lumi_jS);
  // jet_select->SetBranchAddress("vz",&vz_jS);
  // jet_select->SetBranchAddress("npf", &nref_jS);  
  // jet_select->SetBranchAddress("pfpt", &pt_jS);  
  // jet_select->SetBranchAddress("eMax", &eMax_jS);  
  // jet_select->SetBranchAddress("pfeta", &eta_jS);  

  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

  TH1F *hpp_TrgObj80;
  TH1F *hpp_TrgObj60;
  TH1F *hpp_TrgObj40;
  TH1F *hpp_TrgObjComb;

  TH1F *hpp_JEC_TrgObj80;
  TH1F *hpp_JEC_TrgObj60;
  TH1F *hpp_JEC_TrgObj40;
  TH1F *hpp_JEC_TrgObjComb;
  
  TH1F *hpp_Smear_TrgObj80;
  TH1F *hpp_Smear_TrgObj60;
  TH1F *hpp_Smear_TrgObj40;
  TH1F *hpp_Smear_TrgObjComb;
  
  TH1F *hpp_anaBin_TrgObj80;
  TH1F *hpp_anaBin_TrgObj60;
  TH1F *hpp_anaBin_TrgObj40;
  TH1F *hpp_anaBin_TrgObjComb;

    hpp_TrgObj80 = new TH1F(Form("hpp_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_TrgObj60 = new TH1F(Form("hpp_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_TrgObj40 = new TH1F(Form("hpp_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_TrgObjComb = new TH1F(Form("hpp_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000,0,1000);

  hpp_JEC_TrgObj80 = new TH1F(Form("hpp_JEC_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_JEC_TrgObj60 = new TH1F(Form("hpp_JEC_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_JEC_TrgObj40 = new TH1F(Form("hpp_JEC_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_JEC_TrgObjComb = new TH1F(Form("hpp_JEC_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000,0,1000);

  hpp_Smear_TrgObj80 = new TH1F(Form("hpp_Smear_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_Smear_TrgObj60 = new TH1F(Form("hpp_Smear_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_Smear_TrgObj40 = new TH1F(Form("hpp_Smear_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),1000,0,1000);
  hpp_Smear_TrgObjComb = new TH1F(Form("hpp_Smear_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),1000,0,1000);

  hpp_anaBin_TrgObj80 = new TH1F(Form("hpp_anaBin_HLT80_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_TrgObj60 = new TH1F(Form("hpp_anaBin_HLT60_R%d_%s",radius,etaWidth),Form("Spectra from  Jet 60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_TrgObj40 = new TH1F(Form("hpp_anaBin_HLT40_R%d_%s",radius,etaWidth),Form("Spectra from Jet 40 && !jet60 && !jet80 R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);
  hpp_anaBin_TrgObjComb = new TH1F(Form("hpp_anaBin_HLTComb_R%d_%s",radius,etaWidth),Form("Trig Combined Spectra R%d %s ",radius,etaWidth),nbins_pt, boundaries_pt);


  // define the residual correction factors
  TF1 * fResidual = new TF1("fResidual","1 - [0]/pow(x,[1])");
  if(radius == 2){
    fResidual->SetParameter(0, -0.05756);
    fResidual->SetParameter(1,  0.42750);
  }
  if(radius == 3){
    fResidual->SetParameter(0, -0.66229);
    fResidual->SetParameter(1,  1.02119);
  }
  if(radius == 4){
    fResidual->SetParameter(0, -1.28178);
    fResidual->SetParameter(1,  1.17348);
  }
  
  // now start the event loop for each file. 
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpp[0]->GetEntries();
  Long64_t nGoodEvt = 0;
  if(printDebug) nentries = 10;
  TRandom rnd; 

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    jetpp[0]->GetEntry(nEvt);
    jetpp[1]->GetEntry(nEvt);
    jetpp[2]->GetEntry(nEvt);
    jetpp[4]->GetEntry(nEvt);
    jetpp[3]->GetEntry(nEvt);
    evt_select->GetEntry(nEvt);

    //if(printDebug) cout<<"forest values = "<<hiBin_F<<", "<<evt_F<<", "<<run_F<<", "<<lumi_F<<", "<<vz_F<<endl;
    
    // if(pcollisionEventSelection_F==0) continue; 
    // if(fabs(vz_F)>15) continue;
    if(!isGoodEvent_eS) continue; 
    
    jet_select->GetEntry(nGoodEvt);
    ++nGoodEvt;

    //if(hiBin_eS != hiBin_F || evt_eS != evt_F || run_eS != run_F || lumi_eS != lumi_F || vz_eS != vz_F) cout<<"ERROR mismatch eS, F"<<endl;
    //if(hiBin_eS != hiBin_jS || evt_eS != evt_jS || run_eS != run_jS || lumi_eS != lumi_jS || vz_eS != vz_jS) cout<<"ERROR mismatch eS, jS"<<endl;
    //if(hiBin_F != hiBin_jS || evt_F != evt_jS || run_F != run_jS || lumi_F != lumi_jS || vz_F != vz_jS) cout<<"ERROR mismatch F, jS"<<endl;

    //if(printDebug) cout<<"hibin hiForest = "<<hiBin_F<<", evtTree = "<<hiBin_eS<<", jetTree = "<<hiBin_jS<<endl;
    //if(printDebug) cout<<"evt hiForest   = "<<evt_F<<", evtTree = "<<evt_eS<<", jetTree = "<<evt_jS<<endl;
    //if(printDebug) cout<<"lumi hiForest  = "<<lumi_F<<", evtTree = "<<lumi_eS<<", jetTree = "<<lumi_jS<<endl;
    //if(printDebug) cout<<"vz hiForest    = "<<vz_F<<", evtTree = "<<vz_eS<<", jetTree = "<<vz_jS<<endl;
    
    if(printDebug) cout<<"nref_F = "<<nref_F<<", nref_eS = "<<nref_eS<<endl;

    if(nref_F != nref_eS) cout<<"ERROR mismatch in jet counts"<<endl;

    for(int jet = 0; jet<nref_F; ++jet){

      if(fabs(eta_F[jet]) > 2) continue;
      
      if(pt_F[jet] > 40 && pt_F[jet] <= 180) 
	pt_F[jet] = fResidual->Eval(pt_F[jet]) * pt_F[jet];
    
      if(jet40_F == 1 && jet60_F==0 && jet80_F==0) {
	hpp_TrgObj40->Fill(pt_F[jet], jet40_p_F);
	hpp_JEC_TrgObj40->Fill(pt_F[jet]*1.005, jet40_p_F);
	hpp_Smear_TrgObj40->Fill(pt_F[jet]+ rnd.Gaus(0,1), jet40_p_F);
	hpp_anaBin_TrgObj40->Fill(pt_F[jet], jet40_p_F);
      }
    
      if(jet60_F == 1 && jet80_F==0) {
	hpp_TrgObj60->Fill(pt_F[jet]);
	hpp_JEC_TrgObj60->Fill(pt_F[jet]*1.005);
	hpp_Smear_TrgObj60->Fill(pt_F[jet]+rnd.Gaus(0,1));
	hpp_anaBin_TrgObj60->Fill(pt_F[jet]);
      }

      if(jet80_F == 1) {
    	hpp_TrgObj80->Fill(pt_F[jet]);
	hpp_JEC_TrgObj80->Fill(pt_F[jet]*1.005);
	hpp_Smear_TrgObj80->Fill(pt_F[jet]+rnd.Gaus(0,1));
	hpp_anaBin_TrgObj80->Fill(pt_F[jet]);
      }
      
    }// jet loop
    if(printDebug)cout<<endl;

  }// event loop


  hpp_TrgObjComb->Add(hpp_TrgObj80);
  hpp_TrgObjComb->Add(hpp_TrgObj60);
  hpp_TrgObjComb->Add(hpp_TrgObj40);

  divideBinWidth(hpp_TrgObjComb);
  divideBinWidth(hpp_TrgObj80);
  divideBinWidth(hpp_TrgObj60);
  divideBinWidth(hpp_TrgObj40);
  
  hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj80);
  hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj60);
  hpp_JEC_TrgObjComb->Add(hpp_JEC_TrgObj40);

  divideBinWidth(hpp_JEC_TrgObjComb);
  divideBinWidth(hpp_JEC_TrgObj80);
  divideBinWidth(hpp_JEC_TrgObj60);
  divideBinWidth(hpp_JEC_TrgObj40);

  hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj80);
  hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj60);
  hpp_Smear_TrgObjComb->Add(hpp_Smear_TrgObj40);

  divideBinWidth(hpp_Smear_TrgObjComb);
  divideBinWidth(hpp_Smear_TrgObj80);
  divideBinWidth(hpp_Smear_TrgObj60);
  divideBinWidth(hpp_Smear_TrgObj40);

  hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj80);
  hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj60);
  hpp_anaBin_TrgObjComb->Add(hpp_anaBin_TrgObj40);

  divideBinWidth(hpp_anaBin_TrgObjComb);
  divideBinWidth(hpp_anaBin_TrgObj80);
  divideBinWidth(hpp_anaBin_TrgObj60);
  divideBinWidth(hpp_anaBin_TrgObj40);
  
  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
