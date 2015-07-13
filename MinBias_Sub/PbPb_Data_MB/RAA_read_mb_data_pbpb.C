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
  hCentBin->SetBinContent(1, 6.53615);
  hCentBin->SetBinContent(2, 5.82917);
  hCentBin->SetBinContent(3, 5.51244);
  hCentBin->SetBinContent(4, 5.23009);
  hCentBin->SetBinContent(5, 5.0226);
  hCentBin->SetBinContent(6, 4.82387);
  hCentBin->SetBinContent(7, 4.65832);
  hCentBin->SetBinContent(8, 4.46387);
  hCentBin->SetBinContent(9, 4.31041);
  hCentBin->SetBinContent(10, 4.13481);
  hCentBin->SetBinContent(11, 3.97255);
  hCentBin->SetBinContent(12, 3.82127);
  hCentBin->SetBinContent(13, 3.66531);
  hCentBin->SetBinContent(14, 3.55197);
  hCentBin->SetBinContent(15, 3.37913);
  hCentBin->SetBinContent(16, 3.2234);
  hCentBin->SetBinContent(17, 3.1127);
  hCentBin->SetBinContent(18, 2.96605);
  hCentBin->SetBinContent(19, 2.83942);
  hCentBin->SetBinContent(20, 2.71658);
  hCentBin->SetBinContent(21, 2.58391);
  hCentBin->SetBinContent(22, 2.4721);
  hCentBin->SetBinContent(23, 2.36754);
  hCentBin->SetBinContent(24, 2.24715);
  hCentBin->SetBinContent(25, 2.12635);
  hCentBin->SetBinContent(26, 2.03428);
  hCentBin->SetBinContent(27, 1.92792);
  hCentBin->SetBinContent(28, 1.83357);
  hCentBin->SetBinContent(29, 1.73052);
  hCentBin->SetBinContent(30, 1.64517);
  hCentBin->SetBinContent(31, 1.56986);
  hCentBin->SetBinContent(32, 1.47046);
  hCentBin->SetBinContent(33, 1.38777);
  hCentBin->SetBinContent(34, 1.30793);
  hCentBin->SetBinContent(35, 1.23835);
  hCentBin->SetBinContent(36, 1.17554);
  hCentBin->SetBinContent(37, 1.10455);
  hCentBin->SetBinContent(38, 1.03512);
  hCentBin->SetBinContent(39, 0.971748);
  hCentBin->SetBinContent(40, 0.914222);
  hCentBin->SetBinContent(41, 0.858677);
  hCentBin->SetBinContent(42, 0.812601);
  hCentBin->SetBinContent(43, 0.759309);
  hCentBin->SetBinContent(44, 0.709874);
  hCentBin->SetBinContent(45, 0.669402);
  hCentBin->SetBinContent(46, 0.62571);
  hCentBin->SetBinContent(47, 0.596001);
  hCentBin->SetBinContent(48, 0.555403);
  hCentBin->SetBinContent(49, 0.521659);
  hCentBin->SetBinContent(50, 0.489485);
  hCentBin->SetBinContent(51, 0.458461);
  hCentBin->SetBinContent(52, 0.434022);
  hCentBin->SetBinContent(53, 0.408448);
  hCentBin->SetBinContent(54, 0.39129);
  hCentBin->SetBinContent(55, 0.368376);
  hCentBin->SetBinContent(56, 0.34887);
  hCentBin->SetBinContent(57, 0.328423);
  hCentBin->SetBinContent(58, 0.31229);
  hCentBin->SetBinContent(59, 0.29595);
  hCentBin->SetBinContent(60, 0.283978);
  hCentBin->SetBinContent(61, 0.267748);
  hCentBin->SetBinContent(62, 0.255697);
  hCentBin->SetBinContent(63, 0.244729);
  hCentBin->SetBinContent(64, 0.232192);
  hCentBin->SetBinContent(65, 0.220106);
  hCentBin->SetBinContent(66, 0.209342);
  hCentBin->SetBinContent(67, 0.196523);
  hCentBin->SetBinContent(68, 0.191517);
  hCentBin->SetBinContent(69, 0.183065);
  hCentBin->SetBinContent(70, 0.173901);
  hCentBin->SetBinContent(71, 0.164899);
  hCentBin->SetBinContent(72, 0.159027);
  hCentBin->SetBinContent(73, 0.151787);
  hCentBin->SetBinContent(74, 0.144932);
  hCentBin->SetBinContent(75, 0.139798);
  hCentBin->SetBinContent(76, 0.132586);
  hCentBin->SetBinContent(77, 0.127619);
  hCentBin->SetBinContent(78, 0.123863);
  hCentBin->SetBinContent(79, 0.119335);
  hCentBin->SetBinContent(80, 0.112229);
  hCentBin->SetBinContent(81, 0.108363);
  hCentBin->SetBinContent(82, 0.104996);
  hCentBin->SetBinContent(83, 0.101085);
  hCentBin->SetBinContent(84, 0.0952725);
  hCentBin->SetBinContent(85, 0.0922169);
  hCentBin->SetBinContent(86, 0.0900867);
  hCentBin->SetBinContent(87, 0.0860919);
  hCentBin->SetBinContent(88, 0.0820118);
  hCentBin->SetBinContent(89, 0.0800862);
  hCentBin->SetBinContent(90, 0.0754007);
  hCentBin->SetBinContent(91, 0.0739053);
  hCentBin->SetBinContent(92, 0.0720229);
  hCentBin->SetBinContent(93, 0.0685444);
  hCentBin->SetBinContent(94, 0.0671057);
  hCentBin->SetBinContent(95, 0.0643488);
  hCentBin->SetBinContent(96, 0.0615439);
  hCentBin->SetBinContent(97, 0.0601598);
  hCentBin->SetBinContent(98, 0.0579824);
  hCentBin->SetBinContent(99, 0.0557811);
  hCentBin->SetBinContent(100, 0.0529431);
  hCentBin->SetBinContent(100, 0.0529431);
  hCentBin->SetBinContent(101, 0.0512911);
  hCentBin->SetBinContent(102, 0.0503678);
  hCentBin->SetBinContent(103, 0.0476715);
  hCentBin->SetBinContent(104, 0.0461001);
  hCentBin->SetBinContent(105, 0.0447851);
  hCentBin->SetBinContent(106, 0.0443317);
  hCentBin->SetBinContent(107, 0.041507);
  hCentBin->SetBinContent(108, 0.0411886);
  hCentBin->SetBinContent(109, 0.0387893);
  hCentBin->SetBinContent(110, 0.038585);
  hCentBin->SetBinContent(111, 0.0367217);
  hCentBin->SetBinContent(112, 0.0369704);
  hCentBin->SetBinContent(113, 0.0358748);
  hCentBin->SetBinContent(114, 0.0341547);
  hCentBin->SetBinContent(115, 0.0344083);
  hCentBin->SetBinContent(116, 0.0342802);
  hCentBin->SetBinContent(117, 0.0336188);
  hCentBin->SetBinContent(118, 0.0312598);
  hCentBin->SetBinContent(119, 0.031283);
  hCentBin->SetBinContent(120, 0.0315736);
  hCentBin->SetBinContent(121, 0.0313181);
  hCentBin->SetBinContent(122, 0.0296581);
  hCentBin->SetBinContent(123, 0.0307534);
  hCentBin->SetBinContent(124, 0.0301216);
  hCentBin->SetBinContent(125, 0.0289574);
  hCentBin->SetBinContent(126, 0.0288705);
  hCentBin->SetBinContent(127, 0.029107);
  hCentBin->SetBinContent(128, 0.0291955);
  hCentBin->SetBinContent(129, 0.0283446);
  hCentBin->SetBinContent(130, 0.0277929);
  hCentBin->SetBinContent(131, 0.0280661);
  hCentBin->SetBinContent(132, 0.0273745);
  hCentBin->SetBinContent(133, 0.0268596);
  hCentBin->SetBinContent(134, 0.0280446);
  hCentBin->SetBinContent(135, 0.0264748);
  hCentBin->SetBinContent(136, 0.0277904);
  hCentBin->SetBinContent(137, 0.0272159);
  hCentBin->SetBinContent(138, 0.0262413);
  hCentBin->SetBinContent(139, 0.027097);
  hCentBin->SetBinContent(140, 0.0261707);
  hCentBin->SetBinContent(141, 0.027576);
  hCentBin->SetBinContent(142, 0.0264902);
  hCentBin->SetBinContent(143, 0.0265227);
  hCentBin->SetBinContent(144, 0.0267556);
  hCentBin->SetBinContent(145, 0.0261648);
  hCentBin->SetBinContent(146, 0.027463);
  hCentBin->SetBinContent(147, 0.0272955);
  hCentBin->SetBinContent(148, 0.0263357);
  hCentBin->SetBinContent(149, 0.0260893);
  hCentBin->SetBinContent(150, 0.0277501);
  hCentBin->SetBinContent(151, 0.0262621);
  hCentBin->SetBinContent(152, 0.0260947);
  hCentBin->SetBinContent(153, 0.0282848);
  hCentBin->SetBinContent(154, 0.0270858);
  hCentBin->SetBinContent(155, 0.0257997);
  hCentBin->SetBinContent(156, 0.0262174);
  hCentBin->SetBinContent(157, 0.0228979);
  hCentBin->SetBinContent(158, 0.0262223);
  hCentBin->SetBinContent(159, 0.0277375);
  hCentBin->SetBinContent(160, 0.027118);
  hCentBin->SetBinContent(161, 0.0250594);
  hCentBin->SetBinContent(162, 0.0257246);
  hCentBin->SetBinContent(163, 0.0261892);
  hCentBin->SetBinContent(164, 0.0273813);
  hCentBin->SetBinContent(165, 0.0260768);
  hCentBin->SetBinContent(166, 0.0245042);
  hCentBin->SetBinContent(167, 0.0258077);
  hCentBin->SetBinContent(168, 0.024115);
  hCentBin->SetBinContent(169, 0.0274052);
  hCentBin->SetBinContent(170, 0.0263141);
  hCentBin->SetBinContent(171, 0.0251874);
  hCentBin->SetBinContent(172, 0.0270671);
  hCentBin->SetBinContent(173, 0.0294731);
  hCentBin->SetBinContent(174, 0.0233753);
  hCentBin->SetBinContent(175, 0.0225367);
  hCentBin->SetBinContent(176, 0.0272513);
  hCentBin->SetBinContent(177, 0.0248946);
  hCentBin->SetBinContent(178, 0.0241804);
  hCentBin->SetBinContent(179, 0.0235321);
  hCentBin->SetBinContent(180, 0.0250837);
  hCentBin->SetBinContent(181, 0.0250435);
  hCentBin->SetBinContent(182, 0.0237727);
  hCentBin->SetBinContent(183, 0.0280326);
  hCentBin->SetBinContent(184, 0.0248065);
  hCentBin->SetBinContent(185, 0.0194356);
  hCentBin->SetBinContent(186, 0.0272491);
  hCentBin->SetBinContent(187, 0.0226124);
  hCentBin->SetBinContent(188, 0.019056);
  hCentBin->SetBinContent(189, 0.020624);
  hCentBin->SetBinContent(190, 0.0227978);
  hCentBin->SetBinContent(191, 0.0151329);
  hCentBin->SetBinContent(192, 0.0154592);
  hCentBin->SetBinContent(193, 0.0166911);
  hCentBin->SetBinContent(194, 0.0167844);
  hCentBin->SetBinContent(195, 0.0159046);
  hCentBin->SetBinContent(196, 0.0254608);
  hCentBin->SetBinContent(197, 0.0200352);
  hCentBin->SetBinContent(198, 0.0239546);
  hCentBin->SetBinContent(199, 0);
  hCentBin->SetBinContent(200, 0);

  
  // hCentBin->SetBinContent(1, 0.0424662);
  // hCentBin->SetBinContent(2, 0.0425109);
  // hCentBin->SetBinContent(3, 0.0432228);
  // hCentBin->SetBinContent(4, 0.0398556);
  // hCentBin->SetBinContent(5, 0.0407084);
  // hCentBin->SetBinContent(6, 0.0330999);
  // hCentBin->SetBinContent(7, 0.0347667);
  // hCentBin->SetBinContent(8, 0.0330316);
  // hCentBin->SetBinContent(9, 0.0332238);
  // hCentBin->SetBinContent(10, 0.0289776);
  // hCentBin->SetBinContent(11, 0.0285985);
  // hCentBin->SetBinContent(12, 0.028012);
  // hCentBin->SetBinContent(13, 0.0298808);
  // hCentBin->SetBinContent(14, 0.0250504);
  // hCentBin->SetBinContent(15, 0.0228468);
  // hCentBin->SetBinContent(16, 0.0261552);
  // hCentBin->SetBinContent(17, 0.0230045);
  // hCentBin->SetBinContent(18, 0.0219975);
  // hCentBin->SetBinContent(19, 0.0224328);
  // hCentBin->SetBinContent(20, 0.0191567);
  // hCentBin->SetBinContent(21, 0.0186773);
  // hCentBin->SetBinContent(22, 0.0187965);
  // hCentBin->SetBinContent(23, 0.0167002);
  // hCentBin->SetBinContent(24, 0.0201528);
  // hCentBin->SetBinContent(25, 0.0170583);
  // hCentBin->SetBinContent(26, 0.015632);
  // hCentBin->SetBinContent(27, 0.0128937);
  // hCentBin->SetBinContent(28, 0.0142035);
  // hCentBin->SetBinContent(29, 0.013629);
  // hCentBin->SetBinContent(30, 0.0111897);
  // hCentBin->SetBinContent(31, 0.0115882);
  // hCentBin->SetBinContent(32, 0.0113245);
  // hCentBin->SetBinContent(33, 0.0121414);
  // hCentBin->SetBinContent(34, 0.00921685);
  // hCentBin->SetBinContent(35, 0.00907957);
  // hCentBin->SetBinContent(36, 0.00835121);
  // hCentBin->SetBinContent(37, 0.00827403);
  // hCentBin->SetBinContent(38, 0.0076378);
  // hCentBin->SetBinContent(39, 0.00657421);
  // hCentBin->SetBinContent(40, 0.00681596);
  // hCentBin->SetBinContent(41, 0.00611289);
  // hCentBin->SetBinContent(42, 0.00575227);
  // hCentBin->SetBinContent(43, 0.00537931);
  // hCentBin->SetBinContent(44, 0.00525833);
  // hCentBin->SetBinContent(45, 0.00544258);
  // hCentBin->SetBinContent(46, 0.00440711);
  // hCentBin->SetBinContent(47, 0.00386776);
  // hCentBin->SetBinContent(48, 0.00382164);
  // hCentBin->SetBinContent(49, 0.00364078);
  // hCentBin->SetBinContent(50, 0.00359276);
  // hCentBin->SetBinContent(51, 0.00308663);
  // hCentBin->SetBinContent(52, 0.00335243);
  // hCentBin->SetBinContent(53, 0.00307847);
  // hCentBin->SetBinContent(54, 0.00357993);
  // hCentBin->SetBinContent(55, 0.00290412);
  // hCentBin->SetBinContent(56, 0.00247298);
  // hCentBin->SetBinContent(57, 0.00238719);
  // hCentBin->SetBinContent(58, 0.00239202);
  // hCentBin->SetBinContent(59, 0.00236255);
  // hCentBin->SetBinContent(60, 0.00177576);
  // hCentBin->SetBinContent(61, 0.00183271);
  // hCentBin->SetBinContent(62, 0.00152992);
  // hCentBin->SetBinContent(63, 0.00156318);
  // hCentBin->SetBinContent(64, 0.00169289);
  // hCentBin->SetBinContent(65, 0.0016862);
  // hCentBin->SetBinContent(66, 0.00156901);
  // hCentBin->SetBinContent(67, 0.00154053);
  // hCentBin->SetBinContent(68, 0.00127574);
  // hCentBin->SetBinContent(69, 0.00121699);
  // hCentBin->SetBinContent(70, 0.00132383);
  // hCentBin->SetBinContent(71, 0.00121196);
  // hCentBin->SetBinContent(72, 0.00105861);
  // hCentBin->SetBinContent(73, 0.000991154);
  // hCentBin->SetBinContent(74, 0.00108989);
  // hCentBin->SetBinContent(75, 0.00123267);
  // hCentBin->SetBinContent(76, 0.00101845);
  // hCentBin->SetBinContent(77, 0.000888768);
  // hCentBin->SetBinContent(78, 0.000943736);
  // hCentBin->SetBinContent(79, 0.000920097);
  // hCentBin->SetBinContent(80, 0.000663581);
  // hCentBin->SetBinContent(81, 0.00070089);
  // hCentBin->SetBinContent(82, 0.000817824);
  // hCentBin->SetBinContent(83, 0.000745819);
  // hCentBin->SetBinContent(84, 0.00061516);
  // hCentBin->SetBinContent(85, 0.000517464);
  // hCentBin->SetBinContent(86, 0.000559364);
  // hCentBin->SetBinContent(87, 0.000774445);
  // hCentBin->SetBinContent(88, 0.000600886);
  // hCentBin->SetBinContent(89, 0.000497935);
  // hCentBin->SetBinContent(90, 0.000555008);
  // hCentBin->SetBinContent(91, 0.000655111);
  // hCentBin->SetBinContent(92, 0.000562511);
  // hCentBin->SetBinContent(93, 0.000385392);
  // hCentBin->SetBinContent(94, 0.000377753);
  // hCentBin->SetBinContent(95, 0.0006512);
  // hCentBin->SetBinContent(96, 0.000378736);
  // hCentBin->SetBinContent(97, 0.000374605);
  // hCentBin->SetBinContent(98, 0.000388205);
  // hCentBin->SetBinContent(99, 0.000356458);
  // hCentBin->SetBinContent(100, 0.000327789);
  // hCentBin->SetBinContent(101, 0.000284052);
  // hCentBin->SetBinContent(102, 0.000423391);
  // hCentBin->SetBinContent(103, 0.000347779);
  // hCentBin->SetBinContent(104, 0.000385622);
  // hCentBin->SetBinContent(105, 0.000349123);
  // hCentBin->SetBinContent(106, 0.000230372);
  // hCentBin->SetBinContent(107, 0.000337507);
  // hCentBin->SetBinContent(108, 0.00025075);
  // hCentBin->SetBinContent(109, 0.000243498);
  // hCentBin->SetBinContent(110, 0.000422895);
  // hCentBin->SetBinContent(111, 0.000286571);
  // hCentBin->SetBinContent(112, 0.000291885);
  // hCentBin->SetBinContent(113, 0.000177669);
  // hCentBin->SetBinContent(114, 0.000217067);
  // hCentBin->SetBinContent(115, 0.000378736);
  // hCentBin->SetBinContent(116, 0.000270756);
  // hCentBin->SetBinContent(117, 0.000260194);
  // hCentBin->SetBinContent(118, 0.00021788);
  // hCentBin->SetBinContent(119, 0.000197295);
  // hCentBin->SetBinContent(120, 0.000144125);
  // hCentBin->SetBinContent(121, 0.000220356);
  // hCentBin->SetBinContent(122, 0.000267575);
  // hCentBin->SetBinContent(123, 0.000252491);
  // hCentBin->SetBinContent(124, 0.000201993);
  // hCentBin->SetBinContent(125, 0.000176576);
  // hCentBin->SetBinContent(126, 0.000242391);
  // hCentBin->SetBinContent(127, 0.000156381);
  // hCentBin->SetBinContent(128, 0.00019787);
  // hCentBin->SetBinContent(129, 0.000192374);
  // hCentBin->SetBinContent(130, 0.000151495);
  // hCentBin->SetBinContent(131, 0.000183321);
  // hCentBin->SetBinContent(132, 0.000108533);
  // hCentBin->SetBinContent(133, 0.000257407);
  // hCentBin->SetBinContent(134, 0.000257407);
  // hCentBin->SetBinContent(135, 0.000234572);
  // hCentBin->SetBinContent(136, 0.000204119);
  // hCentBin->SetBinContent(137, 0.000215992);
  // hCentBin->SetBinContent(138, 0.000129275);
  // hCentBin->SetBinContent(139, 0.000106937);
  // hCentBin->SetBinContent(140, 0.000167166);
  // hCentBin->SetBinContent(141, 0.000307798);
  // hCentBin->SetBinContent(142, 0.000302989);
  // hCentBin->SetBinContent(143, 7.94725e-05);
  // hCentBin->SetBinContent(144, 0.000397363);
  // hCentBin->SetBinContent(145, 0.000170698);
  // hCentBin->SetBinContent(146, 0.000274405);
  // hCentBin->SetBinContent(147, 0.000153899);
  // hCentBin->SetBinContent(148, 0.000339348);
  // hCentBin->SetBinContent(149, 0);
  // hCentBin->SetBinContent(150, 0.000269324);
  // hCentBin->SetBinContent(151, 0.000285166);
  // hCentBin->SetBinContent(152, 4.75277e-05);
  // hCentBin->SetBinContent(153, 0.000177359);
  // hCentBin->SetBinContent(154, 0.000165267);
  // hCentBin->SetBinContent(155, 0.000346273);
  // hCentBin->SetBinContent(156, 0.000142583);
  // hCentBin->SetBinContent(157, 0);
  // hCentBin->SetBinContent(158, 0.000382723);
  // hCentBin->SetBinContent(159, 0.000484782);
  // hCentBin->SetBinContent(160, 9.32274e-05);
  // hCentBin->SetBinContent(161, 0.000359098);
  // hCentBin->SetBinContent(162, 0.000230849);
  // hCentBin->SetBinContent(163, 0.000605978);
  // hCentBin->SetBinContent(164, 0.000161594);
  // hCentBin->SetBinContent(165, 0.000142583);
  // hCentBin->SetBinContent(166, 0.000346273);
  // hCentBin->SetBinContent(167, 0.000186455);
  // hCentBin->SetBinContent(168, 0.000323188);
  // hCentBin->SetBinContent(169, 0.000151495);
  // hCentBin->SetBinContent(170, 0.000201993);
  // hCentBin->SetBinContent(171, 0.000969565);
  // hCentBin->SetBinContent(172, 0.000201993);
  // hCentBin->SetBinContent(173, 0);
  // hCentBin->SetBinContent(174, 0);
  // hCentBin->SetBinContent(175, 0);
  // hCentBin->SetBinContent(176, 0.00037291);
  // hCentBin->SetBinContent(177, 0.000269324);
  // hCentBin->SetBinContent(178, 0);
  // hCentBin->SetBinContent(179, 0);
  // hCentBin->SetBinContent(180, 0);
  // hCentBin->SetBinContent(181, 0.000484782);
  // hCentBin->SetBinContent(182, 0.000484782);
  // hCentBin->SetBinContent(183, 0.000346273);
  // hCentBin->SetBinContent(184, 0);
  // hCentBin->SetBinContent(185, 0);
  // hCentBin->SetBinContent(186, 0);
  // hCentBin->SetBinContent(187, 0);
  // hCentBin->SetBinContent(188, 0);
  // hCentBin->SetBinContent(189, 0.000403985);
  // hCentBin->SetBinContent(190, 0);
  // hCentBin->SetBinContent(191, 0);
  // hCentBin->SetBinContent(192, 0);
  // hCentBin->SetBinContent(193, 0);
  // hCentBin->SetBinContent(194, 0);
  // hCentBin->SetBinContent(195, 0);
  // hCentBin->SetBinContent(196, 0);
  // hCentBin->SetBinContent(197, 0);
  // hCentBin->SetBinContent(198, 0);
  // hCentBin->SetBinContent(199, 0);
  // hCentBin->SetBinContent(200, 0);

  Double_t value  = hCentBin->GetBinContent(bin);
  delete hCentBin;
  return value;

  
}


using namespace std;

void RAA_read_mb_data_pbpb(int startfile = 8,
			   int endfile = 9,
			   int radius = 3,
			   std::string kFoname="test_output.root",
			   int ptCut = 15){
  
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
  int hiNpix_F;
  int hiBin_F;
  int pcollisionEventSelection_F;
  int pHBHENoiseFilter_F;

  float calopt_F[1000];
  jetpbpb[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jetpbpb[4]->SetBranchAddress("evt",&evt_F);
  jetpbpb[4]->SetBranchAddress("run",&run_F);
  jetpbpb[4]->SetBranchAddress("lumi",&lumi_F);
  jetpbpb[4]->SetBranchAddress("hiBin",&hiBin_F);
  jetpbpb[4]->SetBranchAddress("hiNpix",&hiNpix_F);
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

  float pt_jS[1000];
  float rawpt_jS[1000];
  jet_select->SetBranchAddress("pfpt", pt_jS);  
  jet_select->SetBranchAddress("pfrawpt", rawpt_jS);  
  
  // Declare the output File and the necessary histograms after that:
  // std::string outdir="";
  // std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

  TH1F * hNumerator_80[2][nbins_cent+1], * hNumerator_65[2][nbins_cent+1], * hNumerator_55[2][nbins_cent+1];
  TH1F * hDenominator[2][nbins_cent+1], * hNumerator_Sum[2][nbins_cent+1];

  TH1F * hcentrality_Jet55[nbins_cent];
  TH1F * hcentrality_Jet55not65not80[nbins_cent];
  TH1F * hcentrality_Jet65[nbins_cent];
  TH1F * hcentrality_Jet65not80[nbins_cent];
  TH1F * hcentrality_JetMB_Jet65not80[nbins_cent];
  TH1F * hcentrality_JetMB_Jet80[nbins_cent];
  TH1F * hcentrality_Jet80[nbins_cent];
  TH1F * hcentrality_Jet55_Prescl[nbins_cent];
  TH1F * hcentrality_Jet55not65not80_Prescl[nbins_cent];
  TH1F * hcentrality_Jet65_Prescl[nbins_cent];
  TH1F * hcentrality[nbins_cent];
  TH1F * hcentrality_JetMB[nbins_cent];
  TH1F * hcentrality_JetMB_Jet55[nbins_cent];
  TH1F * hcentrality_JetMB_Jet55_Prescl[nbins_cent];
  TH1F * hcentrality_JetMB_Jet55not65not80[nbins_cent];
  TH1F * hcentrality_JetMB_Jet55not65not80_Prescl[nbins_cent];
  TH1F * hcentrality_JetMB_Jet65[nbins_cent];
  TH1F * hcentrality_JetMB_Jet65_Prescl[nbins_cent];


  //TH1F * hNumerator_HLT80[2][nbins_cent+1], * hNumerator_HLT65[2][nbins_cent+1], * hNumerator_HLT55[2][nbins_cent+1];
  //TH1F * hDenominator[2][nbins_cent+1];

  for(int i = 0; i<nbins_cent+1; i++){

    if( i != nbins_cent ){

      hcentrality[i] = new TH1F(Form("hcentrality_cent%d",i),"",200, 0, 200);
      hcentrality_Jet55[i] = new TH1F(Form("hcentrality_Jet55_cent%d",i),"",200, 0, 200);
      hcentrality_Jet55_Prescl[i] = new TH1F(Form("hcentrality_Jet55_Prescl_cent%d",i),"",200, 0, 200);
      hcentrality_Jet55not65not80[i] = new TH1F(Form("hcentrality_Jet55not65not80_cent%d",i),"",200, 0, 200);
      hcentrality_Jet55not65not80_Prescl[i] = new TH1F(Form("hcentrality_Jet55not65not80_Prescl_cent%d",i),"",200, 0, 200);
      hcentrality_Jet65[i] = new TH1F(Form("hcentrality_Jet65_cent%d",i),"",200, 0, 200);
      hcentrality_Jet65_Prescl[i] = new TH1F(Form("hcentrality_Jet65_Prescl_cent%d",i),"",200, 0, 200);
      hcentrality_Jet65not80[i] = new TH1F(Form("hcentrality_Jet65not80_cent%d",i),"",200, 0, 200);
      hcentrality_Jet80[i] = new TH1F(Form("hcentrality_Jet80_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB[i] = new TH1F(Form("hcentrality_JetMB_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet55[i] = new TH1F(Form("hcentrality_JetMB_Jet55_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet55_Prescl[i] = new TH1F(Form("hcentrality_JetMB_Jet55_Prescl_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet55not65not80[i] = new TH1F(Form("hcentrality_JetMB_Jet55not65not80_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet55not65not80_Prescl[i] = new TH1F(Form("hcentrality_JetMB_Jet55not65not80_Prescl_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet65[i] = new TH1F(Form("hcentrality_JetMB_Jet65_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet65_Prescl[i] = new TH1F(Form("hcentrality_JetMB_Jet65_Prescl_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet65not80[i] = new TH1F(Form("hcentrality_JetMB_Jet65not80_cent%d",i),"",200, 0, 200);
      hcentrality_JetMB_Jet80[i] = new TH1F(Form("hcentrality_JetMB_Jet80_cent%d",i),"",200, 0, 200);
      
    }
    
    hNumerator_80[0][i] = new TH1F(Form("hNumerator_80_noCut_cent%d",i),"",140, 0, 140);
    hNumerator_65[0][i] = new TH1F(Form("hNumerator_65_noCut_cent%d",i),"",140, 0, 140);
    hNumerator_55[0][i] = new TH1F(Form("hNumerator_55_noCut_cent%d",i),"",140, 0, 140);
    hNumerator_Sum[0][i] = new TH1F(Form("hNumerator_Sum_noCut_cent%d",i),"",140, 0, 140);
    hDenominator[0][i] = new TH1F(Form("hDenominator_noCut_cent%d",i),"",140, 0, 140);
    //hRatio_Sum[0][1] = new TH1F(Form("hRatio_Sum_noCut_cent%d",i),"",140,0,140);

    hNumerator_80[1][i] = new TH1F(Form("hNumerator_80_Cut_cent%d",i),"",140, 0, 140);
    hNumerator_65[1][i] = new TH1F(Form("hNumerator_65_Cut_cent%d",i),"",140, 0, 140);
    hNumerator_55[1][i] = new TH1F(Form("hNumerator_55_Cut_cent%d",i),"",140, 0, 140);
    hNumerator_Sum[1][i] = new TH1F(Form("hNumerator_Sum_Cut_cent%d",i),"",140, 0, 140);
    hDenominator[1][i] = new TH1F(Form("hDenominator_Cut_cent%d",i),"",140, 0, 140);
    //hRatio_Sum[1][1] = new TH1F(Form("hRatio_Sum_noCut_cent%d",i),"",140,0,140);

  }
  
  //get the spectra with the specific trigger object from the different files. 
  TH1F *hpbpb_TrgObjMB[nbins_cent];
  TH1F *hpbpb_TrgObjMBwoLJ[nbins_cent];
  TH1F *hpbpb_TrgObjMBwoLJSbJ[nbins_cent];
  TH1F *hpbpb_TrgObjMBLJ[nbins_cent];
  // TH1F *hpbpb_TrgObjMBwoHT[nbins_cent];
  
  TH2F *hdphiptcent[nbins_cent];
  //TH2F *hdphiptMBwoHLTcent[nbins_cent];
  //  TH1F *hBin_[nbins_cent];
  TH1F *hEvent_Vz_[nbins_cent];

  TH2F * hpbpb_forest_jeccheck[nbins_cent];
  TH2F * hpbpb_trees_jeccheck[nbins_cent];
 
  TH1F * hEvent_Vz = new TH1F("hEvent_Vz","Primary Vertex Z",400,-20,20);
  //TH1F * hBin = new TH1F("hBin","Centrality Bins",200,0,200);

  //TH1F *heta = new TH1F("heta","eta distribution",100,-2.5,2.5);
  //TH1F *hphi = new TH1F("hphi","phi distribution",100,-3.5,3.5);
  TH1F *hdphi = new TH1F("hdphi","delta phi distribution",70,0,3.5);
  TH2F *hdphipt=new TH2F("hdphipt","pt vs delta phi distribution",100,0,200,70,0,3.5);
  //TH1F *hdptratio = new TH1F("hdptratio","pt ratio distribution",100,0,5);

    
  for(int i = 0;i<nbins_cent;++i){

    hpbpb_forest_jeccheck[i] = new TH2F(Form("hpbpb_forest_jeccheck_cent%d",i),"reco/raw vs raw pt",300,0,300,50,0,5);
    hpbpb_trees_jeccheck[i] = new TH2F(Form("hpbpb_trees_jeccheck_cent%d",i),"reco/raw vs raw pt",300,0,300,50,0,5);
    
    hEvent_Vz_[i] = new TH1F(Form("hVz_cent_%d",i),Form("HVz Values %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),400,-20,20);
      
   hpbpb_TrgObjMB[i] = new TH1F(Form("hpbpb_HLTMB_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    // hpbpb_TrgObjMBArray[i] = new TH1F(Form("hpbpb_HLTMBArray_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB ArrayR%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObjMBwoLJ[i] = new TH1F(Form("hpbpb_HLTMBwoLJ_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB without LJ R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObjMBwoLJSbJ[i] = new TH1F(Form("hpbpb_HLTMBwoLJSbJ_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB without LJ and SbJ R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObjMBLJ[i] = new TH1F(Form("hpbpb_HLTMBLJ_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from MB Leading Jet R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);    
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
    
    if(pcollisionEventSelection_F == 0 && pHBHENoiseFilter_F == 0) continue; 
    if(fabs(vz_F)>15) continue;
    // if(!isGoodEvent_eS) continue; 

    int jetCounter = 0;
    
    for(int g = 0;g<nref_F;g++){
      
      if(eta_F[g]>=-2 && eta_F[g]<2){ //to select inside 
	
	if(pt_F[g]>=50) jetCounter++;
	
      }//eta selection cut
      
    }// jet loop
    
    // apply the correct supernova selection cut rejection here: 
    if(hiNpix_F > 38000 - 500*jetCounter){
      if(printDebug) cout<<"removed this supernova event"<<endl;
      continue;
    }    

    // jet_select->GetEntry(nGoodEvt);
    // ++nGoodEvt;

    

    int cBin = findBin(hiBin_F);//tells us the centrality of the event. 
    if(cBin == -1) continue;

    hcentrality[cBin]->Fill(hiBin_F);
    if(jet55_F == 1) hcentrality_Jet55[cBin]->Fill(hiBin_F);
    if(jet55_F == 1) hcentrality_Jet55_Prescl[cBin]->Fill(hiBin_F, jet55_p_F);
    if(jet65_F == 1) hcentrality_Jet65[cBin]->Fill(hiBin_F);
    if(jet65_F == 1) hcentrality_Jet65_Prescl[cBin]->Fill(hiBin_F, jet65_p_F);
    if(jet80_F == 1) hcentrality_Jet80[cBin]->Fill(hiBin_F);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80[cBin]->Fill(hiBin_F);
    if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_Jet55not65not80_Prescl[cBin]->Fill(hiBin_F, jet55_p_F);
    if(jet65_F == 1 && jet80_F == 0) hcentrality_Jet65not80[cBin]->Fill(hiBin_F);
    if(jetMB_F) {
      hcentrality_JetMB[cBin]->Fill(hiBin_F);
      if(jet55_F == 1) hcentrality_JetMB_Jet55[cBin]->Fill(hiBin_F);
      if(jet55_F == 1) hcentrality_JetMB_Jet55_Prescl[cBin]->Fill(hiBin_F, jet55_p_F);
      if(jet65_F == 1) hcentrality_JetMB_Jet65[cBin]->Fill(hiBin_F);
      if(jet65_F == 1) hcentrality_JetMB_Jet65_Prescl[cBin]->Fill(hiBin_F, jet65_p_F);
      if(jet80_F == 1) hcentrality_JetMB_Jet80[cBin]->Fill(hiBin_F);
      if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_JetMB_Jet55not65not80[cBin]->Fill(hiBin_F);
      if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hcentrality_JetMB_Jet55not65not80_Prescl[cBin]->Fill(hiBin_F, jet55_p_F);
      if(jet65_F == 1 && jet80_F == 0) hcentrality_JetMB_Jet65not80[cBin]->Fill(hiBin_F);
    }

    if(cBin==nbins_cent) continue;

    if(jetMB_F == 0) continue; 

    //hBin->Fill(hiBin_F);
    hEvent_Vz->Fill(vz_F);
    hEvent_Vz_[cBin]->Fill(vz_F);
    
    Double_t weight_cent = getCentWeight(hiBin_F);
    
    // if(hiBin_eS != hiBin_F || evt_eS != evt_F || run_eS != run_F || lumi_eS != lumi_F || vz_eS != vz_F) cout<<"ERROR mismatch eS, F"<<endl;
    // if(hiBin_eS != hiBin_jS || evt_eS != evt_jS || run_eS != run_jS || lumi_eS != lumi_jS || vz_eS != vz_jS) cout<<"ERROR mismatch eS, jS"<<endl;
    // if(hiBin_F != hiBin_jS || evt_F != evt_jS || run_F != run_jS || lumi_F != lumi_jS || vz_F != vz_jS) cout<<"ERROR mismatch F, jS"<<endl;

    // if(printDebug) cout<<"hibin hiForest = "<<hiBin_F<<", evtTree = "<<hiBin_eS<<", jetTree = "<<hiBin_jS<<endl;
    // if(printDebug) cout<<"evt hiForest   = "<<evt_F<<", evtTree = "<<evt_eS<<", jetTree = "<<evt_jS<<endl;
    // if(printDebug) cout<<"lumi hiForest  = "<<lumi_F<<", evtTree = "<<lumi_eS<<", jetTree = "<<lumi_jS<<endl;
    // if(printDebug) cout<<"vz hiForest    = "<<vz_F<<", evtTree = "<<vz_eS<<", jetTree = "<<vz_jS<<endl;
    
    // if(printDebug) cout<<"nref_F = "<<nref_F<<", nref_eS = "<<nref_eS<<", nref_jS = "<<nref_jS<<endl;

    if(nref_F != nref_eS) cout<<"ERROR mismatch in jet counts"<<endl;
    std::vector < float > phi;
    std::vector < float > pt;
    
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

    // fill the trigger Turnon histograms
    //if(fabs(eta_F[0] < 2.0)){

      // do the jet ID cut]
      bool PFElecCut1 = false;
      Float_t Sumcand1 = chSum_F[0] + phSum_F[0] + neSum_F[0] + muSum_F[0];
      if(isClMatch_eS[0] == 1){
    	if(calopt_eS[0]/pt_F[0] > 0.5 && calopt_eS[0]/pt_F[0] <= 0.85 && eMax_F[0]/Sumcand1 < ((Float_t)18/7 *(Float_t)calopt_eS[0]/pt_F[0] - (Float_t)9/7)) PFElecCut1 = true;
    	if(calopt_eS[0]/pt_F[0] > 0.85) PFElecCut1 = true;
    	if(calopt_eS[0]/pt_F[0] <= 0.5 && eMax_F[0]/Sumcand1 < 0.05) PFElecCut1 = true;
      }
      if(isClMatch_eS[0] == 0)
    	if(eMax_F[0]/Sumcand1 < 0.05 ) PFElecCut1 = true;

      if(PFElecCut1){
	hDenominator[1][cBin]->Fill(pt_F[0]);
	if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hNumerator_55[1][cBin]->Fill(pt_F[0], jet55_p_F);
	if(jet65_F == 1 && jet80_F == 0) hNumerator_65[1][cBin]->Fill(pt_F[0], jet65_p_F);
	if(jet80_F == 1) hNumerator_80[1][cBin]->Fill(pt_F[0], jet80_p_F);
      }
      
      //}// eta cut 
      
    for( int jet = 0; jet<nref_F; jet++ ){
      
      if( fabs(eta_F[jet]) > 2.0 ) continue;

      //if( chMax_F[jet] < 7 && trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;
      if( trkMax_F[jet] < 7 && neMax_F[jet] < 8 ) continue;

      //if( isPFElecCut_eS[jet] != 1 ) continue;

      bool PFElecCut = false;

      Float_t Sumcand = chSum_F[jet] + phSum_F[jet] + neSum_F[jet] + muSum_F[jet];
      if(isClMatch_eS[jet] == 1){
    	if(calopt_eS[jet]/pt_F[jet] > 0.5 && calopt_eS[jet]/pt_F[jet] <= 0.85 && eMax_F[jet]/Sumcand < ((Float_t)18/7 *(Float_t)calopt_eS[jet]/pt_F[jet] - (Float_t)9/7)) PFElecCut = true;
    	if(calopt_eS[jet]/pt_F[jet] > 0.85) PFElecCut = true;
    	if(calopt_eS[jet]/pt_F[jet] <= 0.5 && eMax_F[jet]/Sumcand < 0.05) PFElecCut = true;
      }
      if(isClMatch_eS[jet] == 0)
    	if(eMax_F[jet]/Sumcand < 0.05 ) PFElecCut = true;

      if(!PFElecCut) continue;

      hDenominator[0][cBin]->Fill(pt_F[jet]);
      if(jet55_F == 1 && jet65_F == 0 && jet80_F == 0) hNumerator_55[0][cBin]->Fill(pt_F[jet], jet55_p_F);
      if(jet65_F == 1 && jet80_F == 0) hNumerator_65[0][cBin]->Fill(pt_F[jet], jet65_p_F);
      if(jet80_F == 1) hNumerator_80[0][cBin]->Fill(pt_F[jet], jet80_p_F);
      
      // if(!PFElecCut) continue;
      
      // if(printDebug && (fabs(eta_jS[jet] > 2))) cout<<"jets with |eta| > 2 in jetTree"<<endl;
      // if(printDebug && (fabs(eta_F[jet] > 2)))  cout<<"jets with |eta| > 2 in Forest"<<endl;
      
      Float_t wght = weight_cent; 
      //if(printDebug && index_eS[jet] >= 0 )cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", Calo pT = "<<calopt_F[index_eS[jet]]<<", onFly flag calculation = "<<PFElecCut<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;      // if(printDebug)cout<<jet<<", hiForest pT = "<<pt_F[jet]<<", jetTree pT = "<<pt_jS[jetLoc]<<", electronCut = "<<isPFElecCut_eS[jetLoc]<<", eMax from hiForest = "<<eMax_F[jet]<<", eMax from jet Tree = "<<eMax_jS[jetLoc]<<endl;

      float recpt = pt_F[jet];

      if(recpt <= ptCut) continue;

      hpbpb_forest_jeccheck[cBin]->Fill((Float_t)  rawpt_F[jet], recpt/rawpt_F[jet]);
      hpbpb_trees_jeccheck[cBin]->Fill((Float_t)  rawpt_jS[jet], pt_jS[jet]/rawpt_jS[jet]);

      pt.push_back(recpt);
      phi.push_back(phi_F[jet]);

      hpbpb_TrgObjMB[cBin]->Fill(recpt, wght);

      // if(jet > 0)
      // 	hpbpb_TrgObjMBwoLJ[cBin]->Fill(recpt, wght);
      // if(jet > 1)
      // 	hpbpb_TrgObjMBwoLJSbJ[cBin]->Fill(recpt, wght);
      
            
    }// jet loop

    Float_t wght = weight_cent;
    
    for(int j = 1; j<pt.size(); ++j){
      if(pt[j] < 35) continue; 
      float delphi = deltaphi (phi[0], phi[j]);
      hdphi->Fill(delphi, wght);
      hdphipt->Fill(pt[0], delphi, wght);
      hdphiptcent[cBin]->Fill(pt[0], delphi, wght);
    }
    
    if(pt.size() != 0){
      hpbpb_TrgObjMBLJ[cBin]->Fill(pt[0], wght);
      for(int j = 0; j<pt.size(); ++j){
	if(j >= 1) hpbpb_TrgObjMBwoLJ[cBin]->Fill(pt[j], wght);
	if(j >= 2) hpbpb_TrgObjMBwoLJSbJ[cBin]->Fill(pt[j], wght);      
      }
    }
    
    if(printDebug)cout<<endl;

    pt.clear();
    phi.clear();

  }// event loop

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
