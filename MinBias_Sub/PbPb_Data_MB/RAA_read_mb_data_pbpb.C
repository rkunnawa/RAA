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

float deltaphi(float phi1, float phi2)
{
  float pi=TMath::Pi();
 
  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}



Double_t getCentWeight(int bin){
  
  TH1F * hCentBin = new TH1F("hCentBin","",200,0,200);
  hCentBin->SetBinContent(1, 0.0424662);
  hCentBin->SetBinContent(2, 0.0425109);
  hCentBin->SetBinContent(3, 0.0432228);
  hCentBin->SetBinContent(4, 0.0398556);
  hCentBin->SetBinContent(5, 0.0407084);
  hCentBin->SetBinContent(6, 0.0330999);
  hCentBin->SetBinContent(7, 0.0347667);
  hCentBin->SetBinContent(8, 0.0330316);
  hCentBin->SetBinContent(9, 0.0332238);
  hCentBin->SetBinContent(10, 0.0289776);
  hCentBin->SetBinContent(11, 0.0285985);
  hCentBin->SetBinContent(12, 0.028012);
  hCentBin->SetBinContent(13, 0.0298808);
  hCentBin->SetBinContent(14, 0.0250504);
  hCentBin->SetBinContent(15, 0.0228468);
  hCentBin->SetBinContent(16, 0.0261552);
  hCentBin->SetBinContent(17, 0.0230045);
  hCentBin->SetBinContent(18, 0.0219975);
  hCentBin->SetBinContent(19, 0.0224328);
  hCentBin->SetBinContent(20, 0.0191567);
  hCentBin->SetBinContent(21, 0.0186773);
  hCentBin->SetBinContent(22, 0.0187965);
  hCentBin->SetBinContent(23, 0.0167002);
  hCentBin->SetBinContent(24, 0.0201528);
  hCentBin->SetBinContent(25, 0.0170583);
  hCentBin->SetBinContent(26, 0.015632);
  hCentBin->SetBinContent(27, 0.0128937);
  hCentBin->SetBinContent(28, 0.0142035);
  hCentBin->SetBinContent(29, 0.013629);
  hCentBin->SetBinContent(30, 0.0111897);
  hCentBin->SetBinContent(31, 0.0115882);
  hCentBin->SetBinContent(32, 0.0113245);
  hCentBin->SetBinContent(33, 0.0121414);
  hCentBin->SetBinContent(34, 0.00921685);
  hCentBin->SetBinContent(35, 0.00907957);
  hCentBin->SetBinContent(36, 0.00835121);
  hCentBin->SetBinContent(37, 0.00827403);
  hCentBin->SetBinContent(38, 0.0076378);
  hCentBin->SetBinContent(39, 0.00657421);
  hCentBin->SetBinContent(40, 0.00681596);
  hCentBin->SetBinContent(41, 0.00611289);
  hCentBin->SetBinContent(42, 0.00575227);
  hCentBin->SetBinContent(43, 0.00537931);
  hCentBin->SetBinContent(44, 0.00525833);
  hCentBin->SetBinContent(45, 0.00544258);
  hCentBin->SetBinContent(46, 0.00440711);
  hCentBin->SetBinContent(47, 0.00386776);
  hCentBin->SetBinContent(48, 0.00382164);
  hCentBin->SetBinContent(49, 0.00364078);
  hCentBin->SetBinContent(50, 0.00359276);
  hCentBin->SetBinContent(51, 0.00308663);
  hCentBin->SetBinContent(52, 0.00335243);
  hCentBin->SetBinContent(53, 0.00307847);
  hCentBin->SetBinContent(54, 0.00357993);
  hCentBin->SetBinContent(55, 0.00290412);
  hCentBin->SetBinContent(56, 0.00247298);
  hCentBin->SetBinContent(57, 0.00238719);
  hCentBin->SetBinContent(58, 0.00239202);
  hCentBin->SetBinContent(59, 0.00236255);
  hCentBin->SetBinContent(60, 0.00177576);
  hCentBin->SetBinContent(61, 0.00183271);
  hCentBin->SetBinContent(62, 0.00152992);
  hCentBin->SetBinContent(63, 0.00156318);
  hCentBin->SetBinContent(64, 0.00169289);
  hCentBin->SetBinContent(65, 0.0016862);
  hCentBin->SetBinContent(66, 0.00156901);
  hCentBin->SetBinContent(67, 0.00154053);
  hCentBin->SetBinContent(68, 0.00127574);
  hCentBin->SetBinContent(69, 0.00121699);
  hCentBin->SetBinContent(70, 0.00132383);
  hCentBin->SetBinContent(71, 0.00121196);
  hCentBin->SetBinContent(72, 0.00105861);
  hCentBin->SetBinContent(73, 0.000991154);
  hCentBin->SetBinContent(74, 0.00108989);
  hCentBin->SetBinContent(75, 0.00123267);
  hCentBin->SetBinContent(76, 0.00101845);
  hCentBin->SetBinContent(77, 0.000888768);
  hCentBin->SetBinContent(78, 0.000943736);
  hCentBin->SetBinContent(79, 0.000920097);
  hCentBin->SetBinContent(80, 0.000663581);
  hCentBin->SetBinContent(81, 0.00070089);
  hCentBin->SetBinContent(82, 0.000817824);
  hCentBin->SetBinContent(83, 0.000745819);
  hCentBin->SetBinContent(84, 0.00061516);
  hCentBin->SetBinContent(85, 0.000517464);
  hCentBin->SetBinContent(86, 0.000559364);
  hCentBin->SetBinContent(87, 0.000774445);
  hCentBin->SetBinContent(88, 0.000600886);
  hCentBin->SetBinContent(89, 0.000497935);
  hCentBin->SetBinContent(90, 0.000555008);
  hCentBin->SetBinContent(91, 0.000655111);
  hCentBin->SetBinContent(92, 0.000562511);
  hCentBin->SetBinContent(93, 0.000385392);
  hCentBin->SetBinContent(94, 0.000377753);
  hCentBin->SetBinContent(95, 0.0006512);
  hCentBin->SetBinContent(96, 0.000378736);
  hCentBin->SetBinContent(97, 0.000374605);
  hCentBin->SetBinContent(98, 0.000388205);
  hCentBin->SetBinContent(99, 0.000356458);
  hCentBin->SetBinContent(100, 0.000327789);
  hCentBin->SetBinContent(101, 0.000284052);
  hCentBin->SetBinContent(102, 0.000423391);
  hCentBin->SetBinContent(103, 0.000347779);
  hCentBin->SetBinContent(104, 0.000385622);
  hCentBin->SetBinContent(105, 0.000349123);
  hCentBin->SetBinContent(106, 0.000230372);
  hCentBin->SetBinContent(107, 0.000337507);
  hCentBin->SetBinContent(108, 0.00025075);
  hCentBin->SetBinContent(109, 0.000243498);
  hCentBin->SetBinContent(110, 0.000422895);
  hCentBin->SetBinContent(111, 0.000286571);
  hCentBin->SetBinContent(112, 0.000291885);
  hCentBin->SetBinContent(113, 0.000177669);
  hCentBin->SetBinContent(114, 0.000217067);
  hCentBin->SetBinContent(115, 0.000378736);
  hCentBin->SetBinContent(116, 0.000270756);
  hCentBin->SetBinContent(117, 0.000260194);
  hCentBin->SetBinContent(118, 0.00021788);
  hCentBin->SetBinContent(119, 0.000197295);
  hCentBin->SetBinContent(120, 0.000144125);
  hCentBin->SetBinContent(121, 0.000220356);
  hCentBin->SetBinContent(122, 0.000267575);
  hCentBin->SetBinContent(123, 0.000252491);
  hCentBin->SetBinContent(124, 0.000201993);
  hCentBin->SetBinContent(125, 0.000176576);
  hCentBin->SetBinContent(126, 0.000242391);
  hCentBin->SetBinContent(127, 0.000156381);
  hCentBin->SetBinContent(128, 0.00019787);
  hCentBin->SetBinContent(129, 0.000192374);
  hCentBin->SetBinContent(130, 0.000151495);
  hCentBin->SetBinContent(131, 0.000183321);
  hCentBin->SetBinContent(132, 0.000108533);
  hCentBin->SetBinContent(133, 0.000257407);
  hCentBin->SetBinContent(134, 0.000257407);
  hCentBin->SetBinContent(135, 0.000234572);
  hCentBin->SetBinContent(136, 0.000204119);
  hCentBin->SetBinContent(137, 0.000215992);
  hCentBin->SetBinContent(138, 0.000129275);
  hCentBin->SetBinContent(139, 0.000106937);
  hCentBin->SetBinContent(140, 0.000167166);
  hCentBin->SetBinContent(141, 0.000307798);
  hCentBin->SetBinContent(142, 0.000302989);
  hCentBin->SetBinContent(143, 7.94725e-05);
  hCentBin->SetBinContent(144, 0.000397363);
  hCentBin->SetBinContent(145, 0.000170698);
  hCentBin->SetBinContent(146, 0.000274405);
  hCentBin->SetBinContent(147, 0.000153899);
  hCentBin->SetBinContent(148, 0.000339348);
  hCentBin->SetBinContent(149, 0);
  hCentBin->SetBinContent(150, 0.000269324);
  hCentBin->SetBinContent(151, 0.000285166);
  hCentBin->SetBinContent(152, 4.75277e-05);
  hCentBin->SetBinContent(153, 0.000177359);
  hCentBin->SetBinContent(154, 0.000165267);
  hCentBin->SetBinContent(155, 0.000346273);
  hCentBin->SetBinContent(156, 0.000142583);
  hCentBin->SetBinContent(157, 0);
  hCentBin->SetBinContent(158, 0.000382723);
  hCentBin->SetBinContent(159, 0.000484782);
  hCentBin->SetBinContent(160, 9.32274e-05);
  hCentBin->SetBinContent(161, 0.000359098);
  hCentBin->SetBinContent(162, 0.000230849);
  hCentBin->SetBinContent(163, 0.000605978);
  hCentBin->SetBinContent(164, 0.000161594);
  hCentBin->SetBinContent(165, 0.000142583);
  hCentBin->SetBinContent(166, 0.000346273);
  hCentBin->SetBinContent(167, 0.000186455);
  hCentBin->SetBinContent(168, 0.000323188);
  hCentBin->SetBinContent(169, 0.000151495);
  hCentBin->SetBinContent(170, 0.000201993);
  hCentBin->SetBinContent(171, 0.000969565);
  hCentBin->SetBinContent(172, 0.000201993);
  hCentBin->SetBinContent(173, 0);
  hCentBin->SetBinContent(174, 0);
  hCentBin->SetBinContent(175, 0);
  hCentBin->SetBinContent(176, 0.00037291);
  hCentBin->SetBinContent(177, 0.000269324);
  hCentBin->SetBinContent(178, 0);
  hCentBin->SetBinContent(179, 0);
  hCentBin->SetBinContent(180, 0);
  hCentBin->SetBinContent(181, 0.000484782);
  hCentBin->SetBinContent(182, 0.000484782);
  hCentBin->SetBinContent(183, 0.000346273);
  hCentBin->SetBinContent(184, 0);
  hCentBin->SetBinContent(185, 0);
  hCentBin->SetBinContent(186, 0);
  hCentBin->SetBinContent(187, 0);
  hCentBin->SetBinContent(188, 0);
  hCentBin->SetBinContent(189, 0.000403985);
  hCentBin->SetBinContent(190, 0);
  hCentBin->SetBinContent(191, 0);
  hCentBin->SetBinContent(192, 0);
  hCentBin->SetBinContent(193, 0);
  hCentBin->SetBinContent(194, 0);
  hCentBin->SetBinContent(195, 0);
  hCentBin->SetBinContent(196, 0);
  hCentBin->SetBinContent(197, 0);
  hCentBin->SetBinContent(198, 0);
  hCentBin->SetBinContent(199, 0);
  hCentBin->SetBinContent(200, 0);

  Double_t value  = hCentBin->GetBinContent(bin);
  delete hCentBin;
  return value;

  
}


using namespace std;

void RAA_read_mb_data_pbpb(int startfile = 8,
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

  infile_Forest = "jetRAA_PbPb_MB_data_forest.txt";
  infile_Select = Form("jetRAA_PbPb_MB_data_akPu%d_select.txt",radius);
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
  int jetMB_F;
  int jetMB_p_F;
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
  jetpbpb[0]->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1",&jetMB_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1_Prescl",&jetMB_p_F);
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
  Double_t weight_eS;

  evt_select->SetBranchAddress("hiBin",&hiBin_eS);
  evt_select->SetBranchAddress("run_value",&run_eS);
  evt_select->SetBranchAddress("evt_value",&evt_eS);
  evt_select->SetBranchAddress("lumi_value",&lumi_eS);
  evt_select->SetBranchAddress("vz",&vz_eS);
  evt_select->SetBranchAddress("nref",&nref_eS);
  evt_select->SetBranchAddress("index", index_eS);
  evt_select->SetBranchAddress("isGoodEvt",&isGoodEvent_eS);
  evt_select->SetBranchAddress("isMlMatch",isMiMatch_eS);
  evt_select->SetBranchAddress("isClMatch", isClMatch_eS);
  evt_select->SetBranchAddress("isPFElecCut",isPFElecCut_eS);
  evt_select->SetBranchAddress("isTrackCut",isTrackCut_eS);
  evt_select->SetBranchAddress("isMuCut",isMuCut_eS);
  evt_select->SetBranchAddress("weight", &weight_eS);  

  // jet tree selection file:
  int hiBin_jS;
  int run_jS;
  int evt_jS;
  int lumi_jS;
  float vz_jS;
  int nref_jS;
  float pt_jS[1000];
  float eta_jS[1000];
  float phi_jS[1000];
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
  jet_select->SetBranchAddress("pfphi", &phi_jS);
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
  
  //get the spectra with the specific trigger object from the different files. 
  TH1F *hpbpb_TrgObjMB[nbins_cent];
  TH1F *hpbpb_TrgObjMBwoLJ[nbins_cent];
  TH1F *hpbpb_TrgObjMBLJ[nbins_cent];
  // TH1F *hpbpb_TrgObjMBwoHT[nbins_cent];
  
  TH2F *hdphiptcent[nbins_cent];
  //TH2F *hdphiptMBwoHLTcent[nbins_cent];
  //  TH1F *hBin_[nbins_cent];
  TH1F *hEvent_Vz_[nbins_cent];
  
  TH1F * hEvent_Vz = new TH1F("hEvent_Vz","Primary Vertex Z",400,-20,20);
  TH1F * hBin = new TH1F("hBin","Centrality Bins",200,0,200);

  //TH1F *heta = new TH1F("heta","eta distribution",100,-2.5,2.5);
  //TH1F *hphi = new TH1F("hphi","phi distribution",100,-3.5,3.5);
  TH1F *hdphi = new TH1F("hdphi","delta phi distribution",70,0,3.5);
  TH2F *hdphipt=new TH2F("hdphipt","pt vs delta phi distribution",100,0,200,70,0,3.5);
  //TH1F *hdptratio = new TH1F("hdptratio","pt ratio distribution",100,0,5);

    
  for(int i = 0;i<nbins_cent;++i){
    hEvent_Vz_[i] = new TH1F(Form("hVz_cent_%d",i),Form("HVz Values %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),400,-20,20);
      
   hpbpb_TrgObjMB[i] = new TH1F(Form("hpbpb_HLTMB_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    // hpbpb_TrgObjMBArray[i] = new TH1F(Form("hpbpb_HLTMBArray_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB ArrayR%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjMBwoLJ[i] = new TH1F(Form("hpbpb_HLTMBwoLJ_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB without LJ R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjMBLJ[i] = new TH1F(Form("hpbpb_HLTMBLJ_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB Leading Jet R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);    
    // hpbpb_TrgObjMBwoHT[i] = new TH1F(Form("hpbpb_HLTMBwoHT_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB without HT R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hdphiptcent[i]=new TH2F(Form("hdphipt_cent%d",i),Form("hdphipt%2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200,160,0,3.2);
    // hdphiptwoHLTcent[i]=new TH2F(Form("hdphiptwoHLT_cent%d",i),Form("hdphipt%2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),100,0,200,70,0,3.5);
  }// centrality bin loop
  
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
    
    // if(pcollisionEventSelection_F==0) continue; 
    // if(fabs(vz_F)>15) continue;
    if(!isGoodEvent_eS) continue; 
    
    jet_select->GetEntry(nGoodEvt);
    ++nGoodEvt;

    int cBin = findBin(hiBin_F);//tells us the centrality of the event. 
    if(cBin==-1 || cBin==nbins_cent) continue;

    if(jetMB_F == 1) {

      hBin->Fill(hiBin_F);
      hEvent_Vz->Fill(vz_F);
      hEvent_Vz_[cBin]->Fill(vz_F);
      
    }

    Double_t weight_cent = getCentWeight(hiBin_F);
    
    if(hiBin_eS != hiBin_F || evt_eS != evt_F || run_eS != run_F || lumi_eS != lumi_F || vz_eS != vz_F) cout<<"ERROR mismatch eS, F"<<endl;
    if(hiBin_eS != hiBin_jS || evt_eS != evt_jS || run_eS != run_jS || lumi_eS != lumi_jS || vz_eS != vz_jS) cout<<"ERROR mismatch eS, jS"<<endl;
    if(hiBin_F != hiBin_jS || evt_F != evt_jS || run_F != run_jS || lumi_F != lumi_jS || vz_F != vz_jS) cout<<"ERROR mismatch F, jS"<<endl;

    if(printDebug) cout<<"hibin hiForest = "<<hiBin_F<<", evtTree = "<<hiBin_eS<<", jetTree = "<<hiBin_jS<<endl;
    if(printDebug) cout<<"evt hiForest   = "<<evt_F<<", evtTree = "<<evt_eS<<", jetTree = "<<evt_jS<<endl;
    if(printDebug) cout<<"lumi hiForest  = "<<lumi_F<<", evtTree = "<<lumi_eS<<", jetTree = "<<lumi_jS<<endl;
    if(printDebug) cout<<"vz hiForest    = "<<vz_F<<", evtTree = "<<vz_eS<<", jetTree = "<<vz_jS<<endl;
    
    if(printDebug) cout<<"nref_F = "<<nref_F<<", nref_eS = "<<nref_eS<<", nref_jS = "<<nref_jS<<endl;

    if(nref_F != nref_eS) cout<<"ERROR mismatch in jet counts"<<endl;

    //! Sort the jetTree jets according to pT
    std::vector < Jet > vJet;
    for(int jet2 = 0; jet2<nref_jS; ++jet2){
      //cout <<"\t \t jetTree *** "<< jet2 <<  ", pT " << pt_jS[jet2] <<  ", chSum : "<< chSum_jS[jet2] << endl;
      Jet ijet;
      ijet.id = jet2;
      ijet.pt = pt_jS[jet2];
      vJet.push_back(ijet);
    }
    std::sort (vJet.begin(), vJet.end(), compare_pt);
    std::vector < Jet >::const_iterator itr;

    std::vector < float > phi;
    std::vector < float > pt;

    int jet=0;
    for(itr=vJet.begin(); itr!=vJet.end(); ++itr, ++jet){

      int jetLoc = (*itr).id;
      if(isMultiMatch_jS[jetLoc]) {
	++itr;
	jetLoc = (*itr).id;
	if(itr == vJet.end())  break;
      }
      if(fabs(eta_jS[jetLoc]) > 2) continue;
      //if(isPFElecCut_eS[jet] != 1) continue;
      // if(isMiMatch_eS[jet]) continue;
      if(pt_jS[jetLoc] <15) continue;

      bool PFElecCut = false;

      Float_t Sumcand = chSum_jS[jetLoc] + phSum_jS[jetLoc] + neSum_jS[jetLoc] + muSum_jS[jetLoc];
      if(isCaloMatch_jS[jetLoc] == 1){
	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.5 && calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.85 && eMax_jS[jetLoc]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_jS[jetLoc]/pt_jS[jetLoc] - (Float_t)9/7)) PFElecCut = true;
	if(calopt_jS[jetLoc]/pt_jS[jetLoc] > 0.85) PFElecCut = true;
	if(calopt_jS[jetLoc]/pt_jS[jetLoc] <= 0.5 && eMax_jS[jetLoc]/Sumcand < 0.05) PFElecCut = true;
      }
      if(isCaloMatch_jS[jetLoc] == 0)
	if(eMax_jS[jetLoc]/Sumcand < 0.05 ) PFElecCut = true;

      // if(!PFElecCut) continue;
      
      // if(printDebug && (fabs(eta_jS[jet] > 2))) cout<<"jets with |eta| > 2 in jetTree"<<endl;
      // if(printDebug && (fabs(eta_F[jet] > 2)))  cout<<"jets with |eta| > 2 in Forest"<<endl;

      Float_t wght = weight_cent; 
      
      if(printDebug && index_eS[jet] >= 0 )cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT = "<<calopt_F[index_eS[jet]]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;      // if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;

      float recpt = pt_jS[jetLoc];

      pt.push_back(recpt);
      phi.push_back(phi_jS[jetLoc]);

      hpbpb_TrgObjMB[cBin]->Fill(recpt, wght);

      if(itr != vJet.begin()){
	hpbpb_TrgObjMBwoLJ[cBin]->Fill(recpt, wght);
      }
            
    }// jet loop

    Float_t wght = weight_cent;
    
    if(pt.size() >=2) {
      float delphi = deltaphi (phi[0], phi[1]);
      hdphi->Fill(delphi, wght);
      hdphipt->Fill(pt[0], delphi, wght);
      hdphiptcent[cBin]->Fill(pt[0], delphi, wght);

      hpbpb_TrgObjMBLJ[cBin]->Fill(pt[0], wght);
      
    }
    
    if(printDebug)cout<<endl;

  }// event loop

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
