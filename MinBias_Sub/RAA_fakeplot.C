
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
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include "Headers/plot.h"


static const int nbins_cent = 6;
static const Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
//now we have to multiply by 5, since centrality goes from 0-200. 

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

static const int nbins_pt = 30;
static const double boundaries_pt[nbins_pt+1] = {  3, 4, 5, 7, 9, 12, 15, 18, 21, 24, 28,  32, 37, 43, 49, 56,  64, 74, 84, 97, 114,  133, 153, 174, 196,  220, 245, 300, 330, 362, 395};

void GetCentWeight(TH1F *hCentWeight)
{
  hCentWeight->SetBinContent(1,6.7634);
  hCentWeight->SetBinContent(2,8.9638);
  hCentWeight->SetBinContent(3,8.32666);
  hCentWeight->SetBinContent(4,8.85033);
  hCentWeight->SetBinContent(5,7.48557);
  hCentWeight->SetBinContent(6,7.07842);
  hCentWeight->SetBinContent(7,7.17439);
  hCentWeight->SetBinContent(8,7.39451);
  hCentWeight->SetBinContent(9,7.03825);
  hCentWeight->SetBinContent(10,8.04546);
  hCentWeight->SetBinContent(11,6.32941);
  hCentWeight->SetBinContent(12,4.84289);
  hCentWeight->SetBinContent(13,7.05322);
  hCentWeight->SetBinContent(14,5.6361);
  hCentWeight->SetBinContent(15,5.59001);
  hCentWeight->SetBinContent(16,4.90395);
  hCentWeight->SetBinContent(17,5.13768);
  hCentWeight->SetBinContent(18,4.98226);
  hCentWeight->SetBinContent(19,3.76756);
  hCentWeight->SetBinContent(20,4.44141);
  hCentWeight->SetBinContent(21,4.01054);
  hCentWeight->SetBinContent(22,3.29702);
  hCentWeight->SetBinContent(23,3.21606);
  hCentWeight->SetBinContent(24,3.60559);
  hCentWeight->SetBinContent(25,3.36325);
  hCentWeight->SetBinContent(26,2.6244);
  hCentWeight->SetBinContent(27,3.17479);
  hCentWeight->SetBinContent(28,2.6614);
  hCentWeight->SetBinContent(29,2.1703);
  hCentWeight->SetBinContent(30,2.5898);
  hCentWeight->SetBinContent(31,2.56079);
  hCentWeight->SetBinContent(32,2.4732);
  hCentWeight->SetBinContent(33,2.02533);
  hCentWeight->SetBinContent(34,2.03333);
  hCentWeight->SetBinContent(35,1.75553);
  hCentWeight->SetBinContent(36,1.61111);
  hCentWeight->SetBinContent(37,1.55114);
  hCentWeight->SetBinContent(38,1.63142);
  hCentWeight->SetBinContent(39,1.57826);
  hCentWeight->SetBinContent(40,1.45919);
  hCentWeight->SetBinContent(41,1.5094);
  hCentWeight->SetBinContent(42,1.32344);
  hCentWeight->SetBinContent(43,1.36376);
  hCentWeight->SetBinContent(44,1.00481);
  hCentWeight->SetBinContent(45,0.924943);
  hCentWeight->SetBinContent(46,1.04942);
  hCentWeight->SetBinContent(47,0.976863);
  hCentWeight->SetBinContent(48,0.864318);
  hCentWeight->SetBinContent(49,0.772694);
  hCentWeight->SetBinContent(50,0.864502);
  hCentWeight->SetBinContent(51,0.829401);
  hCentWeight->SetBinContent(52,0.75895);
  hCentWeight->SetBinContent(53,0.702901);
  hCentWeight->SetBinContent(54,0.545314);
  hCentWeight->SetBinContent(55,0.571537);
  hCentWeight->SetBinContent(56,0.447445);
  hCentWeight->SetBinContent(57,0.565426);
  hCentWeight->SetBinContent(58,0.411747);
  hCentWeight->SetBinContent(59,0.381474);
  hCentWeight->SetBinContent(60,0.389027);
  hCentWeight->SetBinContent(61,0.345071);
  hCentWeight->SetBinContent(62,0.39263);
  hCentWeight->SetBinContent(63,0.340061);
  hCentWeight->SetBinContent(64,0.363676);
  hCentWeight->SetBinContent(65,0.342351);
  hCentWeight->SetBinContent(66,0.311693);
  hCentWeight->SetBinContent(67,0.2503);
  hCentWeight->SetBinContent(68,0.258714);
  hCentWeight->SetBinContent(69,0.269137);
  hCentWeight->SetBinContent(70,0.278774);
  hCentWeight->SetBinContent(71,0.254098);
  hCentWeight->SetBinContent(72,0.198022);
  hCentWeight->SetBinContent(73,0.217276);
  hCentWeight->SetBinContent(74,0.233573);
  hCentWeight->SetBinContent(75,0.233931);
  hCentWeight->SetBinContent(76,0.205296);
  hCentWeight->SetBinContent(77,0.187256);
  hCentWeight->SetBinContent(78,0.206262);
  hCentWeight->SetBinContent(79,0.192841);
  hCentWeight->SetBinContent(80,0.174259);
  hCentWeight->SetBinContent(81,0.157487);
  hCentWeight->SetBinContent(82,0.1807);
  hCentWeight->SetBinContent(83,0.135957);
  hCentWeight->SetBinContent(84,0.143054);
  hCentWeight->SetBinContent(85,0.158412);
  hCentWeight->SetBinContent(86,0.158663);
  hCentWeight->SetBinContent(87,0.130637);
  hCentWeight->SetBinContent(88,0.105144);
  hCentWeight->SetBinContent(89,0.109533);
  hCentWeight->SetBinContent(90,0.115536);
  hCentWeight->SetBinContent(91,0.103691);
  hCentWeight->SetBinContent(92,0.0988995);
  hCentWeight->SetBinContent(93,0.0899957);
  hCentWeight->SetBinContent(94,0.091202);
  hCentWeight->SetBinContent(95,0.0947045);
  hCentWeight->SetBinContent(96,0.0990303);
  hCentWeight->SetBinContent(97,0.074485);
  hCentWeight->SetBinContent(98,0.0904833);
  hCentWeight->SetBinContent(99,0.0745771);
  hCentWeight->SetBinContent(100,0.0746246);
  hCentWeight->SetBinContent(101,0.0666776);
  hCentWeight->SetBinContent(102,0.0631808);
  hCentWeight->SetBinContent(103,0.0645528);
  hCentWeight->SetBinContent(104,0.0721828);
  hCentWeight->SetBinContent(105,0.0640522);
  hCentWeight->SetBinContent(106,0.0544978);
  hCentWeight->SetBinContent(107,0.0602298);
  hCentWeight->SetBinContent(108,0.052432);
  hCentWeight->SetBinContent(109,0.0499806);
  hCentWeight->SetBinContent(110,0.05452);
  hCentWeight->SetBinContent(111,0.0456856);
  hCentWeight->SetBinContent(112,0.0464227);
  hCentWeight->SetBinContent(113,0.0389109);
  hCentWeight->SetBinContent(114,0.0429926);
  hCentWeight->SetBinContent(115,0.0423068);
  hCentWeight->SetBinContent(116,0.0436439);
  hCentWeight->SetBinContent(117,0.032317);
  hCentWeight->SetBinContent(118,0.0351724);
  hCentWeight->SetBinContent(119,0.0378572);
  hCentWeight->SetBinContent(120,0.0356574);
  hCentWeight->SetBinContent(121,0.0300515);
  hCentWeight->SetBinContent(122,0.0294732);
  hCentWeight->SetBinContent(123,0.0279459);
  hCentWeight->SetBinContent(124,0.0275134);
  hCentWeight->SetBinContent(125,0.0274872);
  hCentWeight->SetBinContent(126,0.0262874);
  hCentWeight->SetBinContent(127,0.0228082);
  hCentWeight->SetBinContent(128,0.0268362);
  hCentWeight->SetBinContent(129,0.0235638);
  hCentWeight->SetBinContent(130,0.019708);
  hCentWeight->SetBinContent(131,0.0203582);
  hCentWeight->SetBinContent(132,0.0191097);
  hCentWeight->SetBinContent(133,0.0169256);
  hCentWeight->SetBinContent(134,0.018112);
  hCentWeight->SetBinContent(135,0.0175009);
  hCentWeight->SetBinContent(136,0.0144258);
  hCentWeight->SetBinContent(137,0.0155731);
  hCentWeight->SetBinContent(138,0.0135958);
  hCentWeight->SetBinContent(139,0.0129593);
  hCentWeight->SetBinContent(140,0.0134124);
  hCentWeight->SetBinContent(141,0.0102854);
  hCentWeight->SetBinContent(142,0.00902376);
  hCentWeight->SetBinContent(143,0.00938477);
  hCentWeight->SetBinContent(144,0.00979958);
  hCentWeight->SetBinContent(145,0.00981297);
  hCentWeight->SetBinContent(146,0.00830205);
  hCentWeight->SetBinContent(147,0.00828065);
  hCentWeight->SetBinContent(148,0.0075616);
  hCentWeight->SetBinContent(149,0.00721783);
  hCentWeight->SetBinContent(150,0.00742391);
  hCentWeight->SetBinContent(151,0.00668121);
  hCentWeight->SetBinContent(152,0.00490303);
  hCentWeight->SetBinContent(153,0.00689083);
  hCentWeight->SetBinContent(154,0.00620564);
  hCentWeight->SetBinContent(155,0.00501006);
  hCentWeight->SetBinContent(156,0.00467418);
  hCentWeight->SetBinContent(157,0.00358751);
  hCentWeight->SetBinContent(158,0.0043082);
  hCentWeight->SetBinContent(159,0.00353042);
  hCentWeight->SetBinContent(160,0.00356054);
  hCentWeight->SetBinContent(161,0.00277187);
  hCentWeight->SetBinContent(162,0.00259774);
  hCentWeight->SetBinContent(163,0.0026294);
  hCentWeight->SetBinContent(164,0.00266786);
  hCentWeight->SetBinContent(165,0.00251157);
  hCentWeight->SetBinContent(166,0.00218918);
  hCentWeight->SetBinContent(167,0.00229047);
  hCentWeight->SetBinContent(168,0.00178743);
  hCentWeight->SetBinContent(169,0.00182462);
  hCentWeight->SetBinContent(170,0.00204086);
  hCentWeight->SetBinContent(171,0.00189708);
  hCentWeight->SetBinContent(172,0.00203718);
  hCentWeight->SetBinContent(173,0.0020711);
  hCentWeight->SetBinContent(174,0.00180765);
  hCentWeight->SetBinContent(175,0.00159439);
  hCentWeight->SetBinContent(176,0.00216191);
  hCentWeight->SetBinContent(177,0.00136735);
  hCentWeight->SetBinContent(178,0.00182475);
  hCentWeight->SetBinContent(179,0.00160661);
  hCentWeight->SetBinContent(180,0.00138471);
  hCentWeight->SetBinContent(181,0.00156103);
  hCentWeight->SetBinContent(182,0.00200855);
  hCentWeight->SetBinContent(183,0.0023071);
  hCentWeight->SetBinContent(184,0.00211314);
  hCentWeight->SetBinContent(185,0.00155022);
  hCentWeight->SetBinContent(186,0.00204334);
  hCentWeight->SetBinContent(187,0.00180985);
  hCentWeight->SetBinContent(188,0.00165799);
  hCentWeight->SetBinContent(189,0.00253497);
  hCentWeight->SetBinContent(190,0.00271872);
  hCentWeight->SetBinContent(191,0.00223219);
  hCentWeight->SetBinContent(192,0.00272361);
  hCentWeight->SetBinContent(193,0.00296343);
  hCentWeight->SetBinContent(194,0.00455219);
  hCentWeight->SetBinContent(195,0.00947736);
  hCentWeight->SetBinContent(196,0.0159602);
  hCentWeight->SetBinContent(197,0.0463495);
  hCentWeight->SetBinContent(198,0.156464);
  hCentWeight->SetBinContent(199,0);
  hCentWeight->SetBinContent(200,0);
}

using namespace std;


void RAA_fakeplot(char* etaWidth = (char*)"20_eta_20", Int_t radius = 3, Int_t etaLow = 20, Int_t etaHigh = 20){
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // char * etaWidth = (char*)Form("%d_eta_%d",etaLow, etaHigh);
  cout<<"etaWidth = "<<etaWidth<<endl;

  bool isSymm = false;
  if(etaLow == etaHigh) isSymm = true;
  
  TF1 *fVz=0;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);

  // the cut is a 3 step cut based on the different value of the calopt/pfpt - copy the following lines into your loop (with the corresponding branch address set)
  // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) hGood->Fill();
  // if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) ) hGood->Fill();
  // if(calopt/pfpt > 0.85 & eMax/Sumcand > 0.9) hGood->Fill();
  
  TFile * fData, * fMC; 

  fData = TFile::Open("/mnt/hadoop/cms/store/user/pawan/ntuples/JetRaa_akPu234_PbPb_Data.root");
  fMC = TFile::Open("/mnt/hadoop/cms/store/user/pawan/ntuples/JetRaa_akPu234_PbPb_MC.root");

  TTree * Data_matched= (TTree*)fData->Get(Form("akPu%dJetAnalyzer/matchedJets",radius));
  TTree * Data_unmatched = (TTree*)fData->Get(Form("akPu%dJetAnalyzer/unmatchedPFJets",radius));

  TTree * MC_matched = (TTree*)fMC->Get(Form("akPu%dJetAnalyzer/matchedJets",radius));
  TTree * MC_unmatched = (TTree*)fMC->Get(Form("akPu%dJetAnalyzer/unmatchedPFJets",radius));

  TH1F * hMC_Jet55_noCut = new TH1F("hMC_Jet55_noCut","data from matched jets without any jet ID cut",400,0,400);
  TH1F * hMC_Jet55_CutA = new TH1F("hMC_Jet55_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",400,0,400);
  TH1F * hMC_Jet55_CutA_rej = new TH1F("hMC_Jet55_CutA_rej","data from matched jets rejected by Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",400,0,400);

  TH1F * hMC_Jet65_noCut = new TH1F("hMC_Jet65_noCut","data from matched jets without any jet ID cut",400,0,400);
  TH1F * hMC_Jet65_CutA = new TH1F("hMC_Jet65_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",400,0,400);
  TH1F * hMC_Jet65_CutA_rej = new TH1F("hMC_Jet65_CutA_rej","data from matched jets rejected with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",400,0,400);

  TH1F * hMC_Jet80_noCut = new TH1F("hMC_Jet80_noCut","data from matched jets without any jet ID cut",400,0,400);
  TH1F * hMC_Jet80_CutA = new TH1F("hMC_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",400,0,400);
  TH1F * hMC_Jet80_CutA_rej = new TH1F("hMC_Jet80_CutA_rej","data from matched jets rejected with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",400,0,400);

  TH1F * hData_Jet55_noCut = new TH1F("hData_Jet55_noCut","data from matched jets without any jet ID cut",400,0,400);
  TH1F * hData_Jet55_CutA = new TH1F("hData_Jet55_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",400,0,400);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet55_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet55_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet55_CutB->Fill()} 
  TH1F * hData_Jet55_CutA_rej = new TH1F("hData_Jet55_CutA_rej","",400,0,400);

  TH1F * hData_Jet65_noCut = new TH1F("hData_Jet65_noCut","data from matched jets without any jet ID cut",400,0,400);
  TH1F * hData_Jet65_CutA = new TH1F("hData_Jet65_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",400,0,400);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet65_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet65_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet65_CutB->Fill()} 
  TH1F * hData_Jet65_CutA_rej = new TH1F("hData_Jet65_CutA_rej","",400,0,400);
  
  TH1F * hData_Jet80_noCut = new TH1F("hData_Jet80_noCut","data from matched jets without any jet ID cut",400,0,400);
  TH1F * hData_Jet80_CutA = new TH1F("hData_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",400,0,400);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet80_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet80_CutA->Fill()}  
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet80_CutB->Fill()} 
  TH1F * hData_Jet80_CutA_rej = new TH1F("hData_Jet80_CutA_rej","",400,0,400);

  TH1F * hData_unmatched_Jet80_noCut = new TH1F("hData_unmatched_Jet80_noCut","",400,0,400);
  TH1F * hData_unmatched_Jet80_CutA = new TH1F("hData_unmatched_Jet80_CutA","",400,0,400);
  TH1F * hData_unmatched_Jet80_CutA_rej = new TH1F("hData_unmatched_Jet80_CutA_rej","",400,0,400);
  TH1F * hData_unmatched_Jet65_noCut = new TH1F("hData_unmatched_Jet65_noCut","",400,0,400);
  TH1F * hData_unmatched_Jet65_CutA = new TH1F("hData_unmatched_Jet65_CutA","",400,0,400);
  TH1F * hData_unmatched_Jet65_CutA_rej = new TH1F("hData_unmatched_Jet65_CutA_rej","",400,0,400);
  TH1F * hData_unmatched_Jet55_noCut = new TH1F("hData_unmatched_Jet55_noCut","",400,0,400);
  TH1F * hData_unmatched_Jet55_CutA = new TH1F("hData_unmatched_Jet55_CutA","",400,0,400);
  TH1F * hData_unmatched_Jet55_CutA_rej = new TH1F("hData_unmatched_Jet55_CutA_rej","",400,0,400);

  TH1F * hMC_unmatched_Jet80_noCut = new TH1F("hMC_unmatched_Jet80_noCut","",400,0,400);
  TH1F * hMC_unmatched_Jet80_CutA = new TH1F("hMC_unmatched_Jet80_CutA","",400,0,400);
  TH1F * hMC_unmatched_Jet80_CutA_rej = new TH1F("hMC_unmatched_Jet80_CutA_rej","",400,0,400);
  TH1F * hMC_unmatched_Jet65_noCut = new TH1F("hMC_unmatched_Jet65_noCut","",400,0,400);
  TH1F * hMC_unmatched_Jet65_CutA = new TH1F("hMC_unmatched_Jet65_CutA","",400,0,400);
  TH1F * hMC_unmatched_Jet65_CutA_rej = new TH1F("hMC_unmatched_Jet65_CutA_rej","",400,0,400);
  TH1F * hMC_unmatched_Jet55_noCut = new TH1F("hMC_unmatched_Jet55_noCut","",400,0,400);
  TH1F * hMC_unmatched_Jet55_CutA = new TH1F("hMC_unmatched_Jet55_CutA","",400,0,400);
  TH1F * hMC_unmatched_Jet55_CutA_rej = new TH1F("hMC_unmatched_Jet55_CutA_rej","",400,0,400);

  //get the spectra with the specific trigger object from the different files. 
  TH1F *hpbpb_Jet80_gen[nbins_cent],*hpbpb_Jet80_reco[nbins_cent];
  TH1F *hpbpb_Jet65_gen[nbins_cent],*hpbpb_Jet65_reco[nbins_cent];
  TH1F *hpbpb_Jet55_gen[nbins_cent],*hpbpb_Jet55_reco[nbins_cent];
  TH1F *hpbpb_JetComb_gen[nbins_cent],*hpbpb_JetComb_reco[nbins_cent];

  TH1F * hpbpb_MC_noCut[nbins_cent];
  TH1F * hpbpb_MC_Comb_noCut[nbins_cent];
  TH1F * hpbpb_MC_Jet80_noCut[nbins_cent];
  TH1F * hpbpb_MC_Jet65_noCut[nbins_cent];
  TH1F * hpbpb_MC_Jet55_noCut[nbins_cent];
  TH1F * hpbpb_Data_Comb_noCut[nbins_cent];
  TH1F * hpbpb_Data_Jet55_noCut[nbins_cent];
  TH1F * hpbpb_Data_Jet65_noCut[nbins_cent];
  TH1F * hpbpb_Data_Jet80_noCut[nbins_cent];

  TH1F *hpbpb_gen[nbins_cent],*hpbpb_reco[nbins_cent];
  TH2F *hpbpb_matrix[nbins_cent];
  TH2F *hpbpb_matrix_HLT[nbins_cent];
  TH2F *hpbpb_mcclosure_matrix[nbins_cent];
  TH2F *hpbpb_mcclosure_matrix_HLT[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_data[nbins_cent];
  TH1F *hpbpb_mcclosure_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet80_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet65_data[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_data[nbins_cent];
  TH1F *hpbpb_mcclosure_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_JetComb_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet80_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet65_gen[nbins_cent];
  TH1F *hpbpb_mcclosure_Jet55_gen[nbins_cent];

  TH1F *hpbpb_TrgObj80[nbins_cent];
  TH1F *hpbpb_TrgObj65[nbins_cent];
  TH1F *hpbpb_TrgObj55[nbins_cent];
  TH1F *hpbpb_TrgObjComb[nbins_cent];

  TH1F *hpbpb_anaBin_TrgObj80[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObj65[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObj55[nbins_cent];
  TH1F *hpbpb_anaBin_TrgObjComb[nbins_cent];
  TH2F *hpbpb_anaBin_matrix_HLT[nbins_cent];
  TH1F *hpbpb_anaBin_Jet80_gen[nbins_cent],*hpbpb_anaBin_Jet80_reco[nbins_cent];
  TH1F *hpbpb_anaBin_Jet65_gen[nbins_cent],*hpbpb_anaBin_Jet65_reco[nbins_cent];
  TH1F *hpbpb_anaBin_Jet55_gen[nbins_cent],*hpbpb_anaBin_Jet55_reco[nbins_cent];
  TH1F *hpbpb_anaBin_JetComb_gen[nbins_cent],*hpbpb_anaBin_JetComb_reco[nbins_cent];

  TH1F *hpbpb_JEC_TrgObj80[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj65[nbins_cent];
  TH1F *hpbpb_JEC_TrgObj55[nbins_cent];
  TH1F *hpbpb_JEC_TrgObjComb[nbins_cent];

  TH1F *hpbpb_Smear_TrgObj80[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj65[nbins_cent];
  TH1F *hpbpb_Smear_TrgObj55[nbins_cent];
  TH1F *hpbpb_Smear_TrgObjComb[nbins_cent];

  TH1F * hpbpb_fake_JetComb[nbins_cent], * hpbpb_fake_Jet80[nbins_cent], * hpbpb_fake_Jet65[nbins_cent], * hpbpb_fake_Jet55[nbins_cent];  
  
  
  for(int i = 0;i<nbins_cent;++i){
    //cout<<"cent bin = "<<i<<endl;

    hpbpb_fake_JetComb[i] = new TH1F(Form("hpbpb_fake_JetComb_cent%d",i),"",400,0,400);
    hpbpb_fake_Jet80[i] = new TH1F(Form("hpbpb_fake_Jet80_cent%d",i),"",400,0,400);
    hpbpb_fake_Jet65[i] = new TH1F(Form("hpbpb_fake_Jet65_cent%d",i),"",400,0,400);
    hpbpb_fake_Jet55[i] = new TH1F(Form("hpbpb_fake_Jet55_cent%d",i),"",400,0,400);

    hpbpb_MC_noCut[i] = new TH1F(Form("hpbpb_MC_noCut_cent%d",i),"",400,0,400);
    
    hpbpb_MC_Jet55_noCut[i] = new TH1F(Form("hpbpb_MC_Jet55_noCut_cent%d",i),"",400,0,400);
    hpbpb_MC_Jet65_noCut[i] = new TH1F(Form("hpbpb_MC_Jet65_noCut_cent%d",i),"",400,0,400);
    hpbpb_MC_Jet80_noCut[i] = new TH1F(Form("hpbpb_MC_Jet80_noCut_cent%d",i),"",400,0,400);
    hpbpb_MC_Comb_noCut[i] = new TH1F(Form("hpbpb_MC_Comb_noCut_cent%d",i),"",400,0,400);

    hpbpb_Data_Jet55_noCut[i] = new TH1F(Form("hpbpb_Data_Jet55_noCut_cent%d",i),"",400,0,400);
    hpbpb_Data_Jet65_noCut[i] = new TH1F(Form("hpbpb_Data_Jet65_noCut_cent%d",i),"",400,0,400);
    hpbpb_Data_Jet80_noCut[i] = new TH1F(Form("hpbpb_Data_Jet80_noCut_cent%d",i),"",400,0,400);
    hpbpb_Data_Comb_noCut[i] = new TH1F(Form("hpbpb_Data_Comb_noCut_cent%d",i),"",400,0,400);

    hpbpb_TrgObj80[i] = new TH1F(Form("hpbpb_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObj65[i] = new TH1F(Form("hpbpb_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObj55[i] = new TH1F(Form("hpbpb_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_TrgObjComb[i] = new TH1F(Form("hpbpb_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);

    hpbpb_anaBin_TrgObj80[i] = new TH1F(Form("hpbpb_anaBin_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObj65[i] = new TH1F(Form("hpbpb_anaBin_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObj55[i] = new TH1F(Form("hpbpb_anaBin_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_TrgObjComb[i] = new TH1F(Form("hpbpb_anaBin_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    
    hpbpb_JEC_TrgObj80[i] = new TH1F(Form("hpbpb_JEC_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_JEC_TrgObj65[i] = new TH1F(Form("hpbpb_JEC_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_JEC_TrgObj55[i] = new TH1F(Form("hpbpb_JEC_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_JEC_TrgObjComb[i] = new TH1F(Form("hpbpb_JEC_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);

    hpbpb_Smear_TrgObj80[i] = new TH1F(Form("hpbpb_Smear_HLT80_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Smear_TrgObj65[i] = new TH1F(Form("hpbpb_Smear_HLT65_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from  Jet 65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Smear_TrgObj55[i] = new TH1F(Form("hpbpb_Smear_HLT55_R%d_%s_cent%d",radius,etaWidth,i),Form("Spectra from Jet 55 && !jet65 && !jet80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Smear_TrgObjComb[i] = new TH1F(Form("hpbpb_Smear_HLTComb_R%d_%s_cent%d",radius,etaWidth,i),Form("Trig Combined Spectra R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);

    hpbpb_gen[i] = new TH1F(Form("hpbpb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    //cout<<"A"<<endl;
    hpbpb_reco[i] = new TH1F(Form("hpbpb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("Reco jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    //cout<<"B"<<endl;
    hpbpb_matrix[i] = new TH2F(Form("hpbpb_matrix_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400,400,0,400);
    
    
    hpbpb_matrix_HLT[i] = new TH2F(Form("hpbpb_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400,400,0,400);

    hpbpb_anaBin_matrix_HLT[i] = new TH2F(Form("hpbpb_anaBin_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix refpt jtpt from trigger addition R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt,nbins_pt, boundaries_pt);

    hpbpb_mcclosure_matrix[i] = new TH2F(Form("hpbpb_mcclosure_matrix_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400,400,0,400);
    hpbpb_mcclosure_matrix_HLT[i] = new TH2F(Form("hpbpb_mcclosure_matrix_HLT_R%d_%s_cent%d",radius,etaWidth,i),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400,400,0,400);
    //cout<<"C"<<endl;
    hpbpb_mcclosure_data[i] =new TH1F(Form("hpbpb_mcclosure_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_JetComb_data[i] = new TH1F(Form("hpbpb_mcclosure_JetComb_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger combined  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_Jet80_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet80_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 80  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_Jet65_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet65_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 65  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_Jet55_data[i] = new TH1F(Form("hpbpb_mcclosure_Jet55_data_R%d_%s_cent%d",radius,etaWidth,i),Form("data for unfolding mc closure test trigger 55  R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);

    hpbpb_mcclosure_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_JetComb_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_JetComb_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_Jet80_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet80_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_Jet65_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet65_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 65 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_mcclosure_Jet55_gen[i] = new TH1F(Form("hpbpb_mcclosure_gen_Jet55_R%d_%s_cent%d",radius,etaWidth,i),Form("gen spectra for unfolding mc closure test trigger 55 R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);

    
    hpbpb_JetComb_gen[i] = new TH1F(Form("hpbpb_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_JetComb_reco[i] = new TH1F(Form("hpbpb_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Jet80_gen[i] = new TH1F(Form("hpbpb_Jet80_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Jet80_reco[i] = new TH1F(Form("hpbpb_Jet80_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Jet65_gen[i] = new TH1F(Form("hpbpb_Jet65_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Jet65_reco[i] = new TH1F(Form("hpbpb_Jet65_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Jet55_gen[i] = new TH1F(Form("hpbpb_Jet55_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);
    hpbpb_Jet55_reco[i] = new TH1F(Form("hpbpb_Jet55_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),400,0,400);

    hpbpb_anaBin_JetComb_gen[i] = new TH1F(Form("hpbpb_anaBin_JetComb_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_JetComb_reco[i] = new TH1F(Form("hpbpb_anaBin_JetComb_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from HLT trigger combined R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet80_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet80_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet80_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet80_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet65_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet65_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet65_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet65_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet55_gen[i] = new TH1F(Form("hpbpb_anaBin_Jet55_gen_R%d_%s_cent%d",radius,etaWidth,i),Form("Gen refpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);
    hpbpb_anaBin_Jet55_reco[i] = new TH1F(Form("hpbpb_anaBin_Jet55_reco_R%d_%s_cent%d",radius,etaWidth,i),Form("reco jtpt from Jet55 && !Jet65 && !Jet80 trigger R%d %s %2.0f - %2.0f cent",radius,etaWidth,5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt, boundaries_pt);

    
	
  }// centrality bin loop
  


  // Define all the histograms necessary for the analysis: 
  
  // 1 - Data, 2 - MC
  Float_t pfpt_1, pfpt_2;
  Float_t pfrefpt_2;
  Float_t calopt_1, calopt_2;
  Float_t eMax_1, eMax_2;
  Float_t chMax_1, chMax_2;
  Float_t chSum_1, chSum_2;
  Float_t phSum_1, phSum_2;
  Float_t neSum_1, neSum_2;
  Float_t muSum_1, muSum_2;
  Int_t jet55_1, jet65_1, jet80_1;
  Int_t jet55_p_1, jet65_p_1, jet80_p_1;
  Int_t jet55_2, jet65_2, jet80_2;
  Int_t jet55_p_2;
  Float_t weight;
  Int_t subid_2;
  Int_t hiBin_1, hiBin_2;
  Float_t eta_1, eta_2;
  Float_t vz_2; 
  Float_t rawpt_1, rawpt_2; 

  Data_matched->SetBranchAddress("calopt",&calopt_1);
  Data_matched->SetBranchAddress("pfpt",&pfpt_1);
  Data_matched->SetBranchAddress("eMax",&eMax_1);
  Data_matched->SetBranchAddress("chMax",&chMax_1);
  Data_matched->SetBranchAddress("chSum",&chSum_1);
  Data_matched->SetBranchAddress("phSum",&phSum_1);
  Data_matched->SetBranchAddress("neSum",&neSum_1);
  Data_matched->SetBranchAddress("muSum",&muSum_1);
  Data_matched->SetBranchAddress("hiBin",&hiBin_1);
  Data_matched->SetBranchAddress("jet55",&jet55_1);
  Data_matched->SetBranchAddress("jet65",&jet65_1);
  Data_matched->SetBranchAddress("jet80",&jet80_1);
  Data_matched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_matched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_matched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  Data_matched->SetBranchAddress("pfeta",&eta_1);
  Data_matched->SetBranchAddress("pfrawpt",&rawpt_1);
  
  Data_unmatched->SetBranchAddress("pfpt",&pfpt_1);
  Data_unmatched->SetBranchAddress("eMax",&eMax_1);
  Data_unmatched->SetBranchAddress("chMax",&chMax_1);
  Data_unmatched->SetBranchAddress("chSum",&chSum_1);
  Data_unmatched->SetBranchAddress("phSum",&phSum_1);
  Data_unmatched->SetBranchAddress("neSum",&neSum_1);
  Data_unmatched->SetBranchAddress("muSum",&muSum_1);
  Data_unmatched->SetBranchAddress("hiBin",&hiBin_1);
  Data_unmatched->SetBranchAddress("jet55",&jet55_1);
  Data_unmatched->SetBranchAddress("jet65",&jet65_1);
  Data_unmatched->SetBranchAddress("jet80",&jet80_1);
  Data_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_unmatched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_unmatched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  Data_unmatched->SetBranchAddress("pfeta",&eta_1);
  Data_unmatched->SetBranchAddress("pfrawpt",&rawpt_1);

  MC_matched->SetBranchAddress("calopt",&calopt_2);
  MC_matched->SetBranchAddress("pfpt",&pfpt_2);
  MC_matched->SetBranchAddress("eMax",&eMax_2);
  MC_matched->SetBranchAddress("chMax",&chMax_2);
  MC_matched->SetBranchAddress("chSum",&chSum_2);
  MC_matched->SetBranchAddress("phSum",&phSum_2);
  MC_matched->SetBranchAddress("neSum",&neSum_2);
  MC_matched->SetBranchAddress("muSum",&muSum_2);
  MC_matched->SetBranchAddress("hiBin",&hiBin_2);
  MC_matched->SetBranchAddress("refpt",&pfrefpt_2);
  MC_matched->SetBranchAddress("jet55",&jet55_2);
  MC_matched->SetBranchAddress("jet65",&jet65_2);
  MC_matched->SetBranchAddress("jet80",&jet80_2);
  MC_matched->SetBranchAddress("weight", &weight);
  MC_matched->SetBranchAddress("subid", &subid_2);
  MC_matched->SetBranchAddress("jet55_prescl",&jet55_p_2);
  MC_matched->SetBranchAddress("pfeta",&eta_2);
  MC_matched->SetBranchAddress("vz",&vz_2);
  MC_matched->SetBranchAddress("pfrawpt",&rawpt_2);

  MC_unmatched->SetBranchAddress("pfpt",&pfpt_2);
  MC_unmatched->SetBranchAddress("eMax",&eMax_2);
  MC_unmatched->SetBranchAddress("chMax",&chMax_2);
  MC_unmatched->SetBranchAddress("chSum",&chSum_2);
  MC_unmatched->SetBranchAddress("phSum",&phSum_2);
  MC_unmatched->SetBranchAddress("neSum",&neSum_2);
  MC_unmatched->SetBranchAddress("muSum",&muSum_2);
  MC_unmatched->SetBranchAddress("hiBin",&hiBin_2);
  MC_unmatched->SetBranchAddress("refpt",&pfrefpt_2);
  MC_unmatched->SetBranchAddress("jet55",&jet55_2);
  MC_unmatched->SetBranchAddress("jet65",&jet65_2);
  MC_unmatched->SetBranchAddress("jet80",&jet80_2);
  MC_unmatched->SetBranchAddress("weight", & weight);
  MC_unmatched->SetBranchAddress("subid", &subid_2);
  MC_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_2);
  MC_unmatched->SetBranchAddress("pfeta",&eta_2);
  MC_unmatched->SetBranchAddress("vz",&vz_2);
  MC_unmatched->SetBranchAddress("pfrawpt",&rawpt_2);

  
  // Lets implement the correction factors from the fake rate study in the MC to the Data.
  // the method is to take the plots from Marguerite and read off the values for each radii, centrality and pT bin.
  // this correction factor is actually (1 - fake rate) 


  //3,4 ,5 , 7, 9,12,15,18,21,24,28 ,32 , 37  ,  43 , 49  ,  56 ,  64 ,74,84,97,114,133,153,174,196,220,245,300...    
  Float_t R2nCorr[nbins_cent][nbins_pt] = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.75, 0.94, 0.985, 0.99, 0.995, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //0-5%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.80, 0.96, 0.985, 0.99, 0.995, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //5-10%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.91, 0.95, 0.985, 0.99, 0.995, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //10-30%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.92, 0.95, 0.985, 0.99, 0.995, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //30-50%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.92, 0.95, 0.985, 0.99, 0.995, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //50-70%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.91, 0.95, 0.985, 0.99, 0.995, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} //70-90%
  };
  Float_t R3nCorr[nbins_cent][nbins_pt] = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.58, 0.85, 0.99, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //0-5%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.64, 0.91, 0.99, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //5-10%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.75, 0.97, 0.995, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //10-30%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.88, 0.99, 0.995, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //30-50%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.94, 0.99, 0.995, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //50-70%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.95, 0.99, 0.995, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} //70-90%
  };
  Float_t R4nCorr[nbins_cent][nbins_pt] = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.36, 0.52, 0.86, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //0-5%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.39, 0.60, 0.91, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //5-10%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.45, 0.78, 0.99, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //10-30%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.69, 0.96, 0.99, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //30-50%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.91, 0.99, 0.995, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //50-70%
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.97, 0.99, 0.995, 1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} //70-90%
  };
  
  long entries = Data_matched->GetEntries();

#if 0  
  // data loop
  //entries = 1000;
  Float_t Jet55_prescl = 2.0475;
  // get the random value for smear systematics, TRandom rnd, value per jet = rnd.Gaus(0,1);
  TRandom rnd; 
  TH1F * htest = new TH1F("htest","",nbins_pt, boundaries_pt);
  cout<<"matched Data ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry ){
    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    Data_matched->GetEntry(nentry);
    Int_t cBin = findBin(hiBin_1);
    if(cBin == -1 || cBin >= nbins_cent) continue;
    if(isSymm && TMath::Abs(eta_1) > (Float_t)etaHigh/10) continue;       
    if(!isSymm && (TMath::Abs(eta_1) < (Float_t)etaLow/10 || TMath::Abs(eta_1) > (Float_t)etaHigh/10)) continue;
       
    Float_t Sumcand = chSum_1 + phSum_1 + neSum_1 + muSum_1;
    Float_t wght = 1; 
    // if(etaBoundary == 1.0){ 
    //   if(radius == 2) wght = R2nCorr[cBin][htest->FindBin(pfpt_1)];
    //   if(radius == 3) wght = R3nCorr[cBin][htest->FindBin(pfpt_1)];
    //   if(radius == 4) wght = R4nCorr[cBin][htest->FindBin(pfpt_1)];  
    // }
    if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {
      
      hData_Jet55_noCut->Fill(pfpt_1, Jet55_prescl* wght);
      hpbpb_Data_Jet55_noCut[cBin]->Fill(pfpt_1, Jet55_prescl* wght);
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) {
	hData_Jet55_CutA->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_anaBin_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl* wght);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl* wght);
      }
      if(calopt_1/pfpt_1 > 0.85){
	hData_Jet55_CutA->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_anaBin_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl* wght);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl* wght);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand < 0.05) {
	hData_Jet55_CutA->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_anaBin_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl* wght);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl* wght);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl* wght);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet55_CutA_rej->Fill(pfpt_1, Jet55_prescl* wght);
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) hData_Jet55_CutA_rej->Fill(pfpt_1, Jet55_prescl* wght);
      
    }
    
    if(jet65_1 == 1 && jet80_1 == 0 ) {
      
      hData_Jet65_noCut->Fill(pfpt_1, wght);
      hpbpb_Data_Jet65_noCut[cBin]->Fill(pfpt_1, wght);
	    
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)){
	hData_Jet65_CutA->Fill(pfpt_1, wght);
	hpbpb_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
      }
      if(calopt_1/pfpt_1 > 0.85) {
	hData_Jet65_CutA->Fill(pfpt_1, wght);
	hpbpb_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand < 0.05) {
	hData_Jet65_CutA->Fill(pfpt_1, wght);
	hpbpb_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet65_CutA_rej->Fill(pfpt_1, wght);
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) hData_Jet65_CutA_rej->Fill(pfpt_1, wght);
      
    }
    if(jet80_1 == 1) {
    
      hData_Jet80_noCut->Fill(pfpt_1, wght);
      hpbpb_Data_Jet80_noCut[cBin]->Fill(pfpt_1, wght);   
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand < ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) {
	hData_Jet80_CutA->Fill(pfpt_1, wght);
	hpbpb_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
      }
      if(calopt_1/pfpt_1 > 0.85){
	hData_Jet80_CutA->Fill(pfpt_1, wght);
	hpbpb_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand < 0.05){
	hData_Jet80_CutA->Fill(pfpt_1, wght);
	hpbpb_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
      }
      if(calopt_1/pfpt_1 <= 0.5 && eMax_1/Sumcand >= 0.05) hData_Jet80_CutA_rej->Fill(pfpt_1, wght);
      if(calopt_1/pfpt_1 > 0.5 && calopt_1/pfpt_1 <= 0.85 && eMax_1/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_1/pfpt_1 - (Float_t)9/7)) hData_Jet80_CutA_rej->Fill(pfpt_1, wght);
      
    }
    
  }// data ntuple loop
  // data unmatched loop:
  entries = Data_unmatched->GetEntries();
  //entries = 1000;
  cout<<"Unmatched Data ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry ){
    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    Data_unmatched->GetEntry(nentry);
    Int_t cBin = findBin(hiBin_1);
    if(cBin == -1 || cBin >= nbins_cent) continue;
    
    if(isSymm && TMath::Abs(eta_1) > (Float_t)etaHigh/10) continue;       
    if(!isSymm && (TMath::Abs(eta_1) < (Float_t)etaLow/10 || TMath::Abs(eta_1) > (Float_t)etaHigh/10)) continue;
    
    Float_t wght = 1; 
    // if(etaBoundary == 1.0){
    //   if(radius == 2) wght = R2nCorr[cBin][htest->FindBin(pfpt_1)];
    //   if(radius == 3) wght = R3nCorr[cBin][htest->FindBin(pfpt_1)];
    //   if(radius == 4) wght = R4nCorr[cBin][htest->FindBin(pfpt_1)];    
    // }
    Float_t Sumcand = chSum_1 + phSum_1 + neSum_1 + muSum_1;
    if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {
    
      hData_unmatched_Jet55_noCut->Fill(pfpt_1, Jet55_prescl*wght);
      hpbpb_Data_Jet55_noCut[cBin]->Fill(pfpt_1, Jet55_prescl*wght);
      if(eMax_1/Sumcand < 0.05 ){
	hpbpb_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl*wght);
	hpbpb_anaBin_TrgObj55[cBin]->Fill(pfpt_1, Jet55_prescl*wght);
	hData_unmatched_Jet55_CutA->Fill(pfpt_1, Jet55_prescl*wght);
	hpbpb_JEC_TrgObj55[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), Jet55_prescl*wght);
	hpbpb_Smear_TrgObj55[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), Jet55_prescl*wght);
      }
      else hData_unmatched_Jet55_CutA_rej->Fill(pfpt_1, Jet55_prescl*wght);
      
    }
    if(jet65_1 == 1 && jet80_1 == 0 ) {
      hData_unmatched_Jet65_noCut->Fill(pfpt_1, wght);
      hpbpb_Data_Jet65_noCut[cBin]->Fill(pfpt_1, wght);
      if(eMax_1/Sumcand < 0.05  ){
	hpbpb_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj65[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj65[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj65[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
	hData_unmatched_Jet65_CutA->Fill(pfpt_1, wght);}
      else hData_unmatched_Jet65_CutA_rej->Fill(pfpt_1, wght);
      
    }
    if(jet80_1 == 1) {
    
      hData_unmatched_Jet80_noCut->Fill(pfpt_1, wght);
      hpbpb_Data_Jet80_noCut[cBin]->Fill(pfpt_1, wght);
      if(eMax_1/Sumcand < 0.05  ){
	hpbpb_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_anaBin_TrgObj80[cBin]->Fill(pfpt_1, wght);
	hpbpb_JEC_TrgObj80[cBin]->Fill(pfpt_1 * (1. + 0.02/nbins_cent*(nbins_cent-cBin)), wght);
	hpbpb_Smear_TrgObj80[cBin]->Fill(pfpt_1 + rnd.Gaus(0,1), wght);
	hData_unmatched_Jet80_CutA->Fill(pfpt_1, wght);}
      else hData_unmatched_Jet80_CutA_rej->Fill(pfpt_1, wght);
      
    }
    
  }// data ntuple loop
#endif

  //! Centrality re-weighting 
  TH1F *hCentWeight = new TH1F("hCentWeight","Centrality weight",200,0,200);
  GetCentWeight(hCentWeight);

  entries = MC_matched->GetEntries();
  //entries = 1000;
  // MC loop
  cout<<" looking at matched MC ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    MC_matched->GetEntry(nentry);

    if(isSymm && TMath::Abs(eta_2) > (Float_t)etaHigh/10) continue;       
    if(!isSymm && (TMath::Abs(eta_2) < (Float_t)etaLow/10 || TMath::Abs(eta_2) > (Float_t)etaHigh/10)) continue;
    
    Int_t cBin = findBin(hiBin_2);
    if(cBin == -1 || cBin >= nbins_cent) continue;    
    
    Float_t cent_weight = hCentWeight->GetBinContent(hiBin_2);
    Float_t vz_weight = fVz->Eval(vz_2);
    
    Float_t Sumcand = chSum_2 + phSum_2 + neSum_2 + muSum_2;

    if(subid_2!=0){
      
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)){
	if(jet55_2==1 && jet65_2==0 && jet80_2==0) hpbpb_fake_Jet55[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet65_2==1 && jet80_2==0) hpbpb_fake_Jet65[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet80_2==1) hpbpb_fake_Jet80[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
      }
      if(calopt_2/pfpt_2 > 0.85) {
	if(jet55_2==1 && jet65_2==0 && jet80_2==0) hpbpb_fake_Jet55[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet65_2==1 && jet80_2==0) hpbpb_fake_Jet65[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet80_2==1) hpbpb_fake_Jet80[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
	if(jet55_2==1 && jet65_2==0 && jet80_2==0) hpbpb_fake_Jet55[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet65_2==1 && jet80_2==0) hpbpb_fake_Jet65[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet80_2==1) hpbpb_fake_Jet80[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
      }
    }

#if 0
    if(subid_2 != 0) continue;
    if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)){
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
    }
    if(calopt_2/pfpt_2 > 0.85) {
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
    }
    if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
    }
    hpbpb_MC_noCut[cBin]->Fill(pfrefpt_2, weight);
    
    if(jet55_2 == 1 && jet65_2==0 && jet80_2 == 0){
      
      hMC_Jet55_noCut->Fill(pfrefpt_2, weight);
      hpbpb_MC_Jet55_noCut[cBin]->Fill(pfrefpt_2,weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)){
	hMC_Jet55_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
      }
      if(calopt_2/pfpt_2 > 0.85) {
	hMC_Jet55_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
	hMC_Jet55_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet55_CutA_rej->Fill(pfrefpt_2, weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) hMC_Jet55_CutA_rej->Fill(pfrefpt_2, weight);
      
    }
    
    if(jet65_2 == 1 && jet80_2 == 0){
      hMC_Jet65_noCut->Fill(pfrefpt_2, weight);
      hpbpb_MC_Jet65_noCut[cBin]->Fill(pfrefpt_2,weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) {
	hMC_Jet65_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
      }
      if(calopt_2/pfpt_2 > 0.85) {
	hMC_Jet65_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
	hMC_Jet65_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet65_CutA_rej->Fill(pfrefpt_2, weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) hMC_Jet65_CutA_rej->Fill(pfrefpt_2, weight);
    }
    
    if(jet80_2 == 1){
      hMC_Jet80_noCut->Fill(pfrefpt_2, weight);
      hpbpb_MC_Jet80_noCut[cBin]->Fill(pfrefpt_2,weight);
 
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand < ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) {
	hMC_Jet80_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	
      }
      if(calopt_2/pfpt_2 > 0.85) {
	hMC_Jet80_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
      }
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand < 0.05) {
	hMC_Jet80_CutA->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hpbpb_anaBin_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
      }
      
      if(calopt_2/pfpt_2 <= 0.5 && eMax_2/Sumcand >= 0.05) hMC_Jet80_CutA_rej->Fill(pfrefpt_2, weight);
      if(calopt_2/pfpt_2 > 0.5 && calopt_2/pfpt_2 <= 0.85 && eMax_2/Sumcand >= ((Float_t)18/7 *(Float_t)calopt_2/pfpt_2 - (Float_t)9/7)) hMC_Jet80_CutA_rej->Fill(pfrefpt_2, weight);
    }
#endif
    
  }// mc ntuple loop


  entries = MC_unmatched->GetEntries();
  //entries = 1000;
  // MC loop
  cout<<" looking at unmatched MC ntuple"<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%100000 == 0) cout<<nentry<<"/"<<entries<<endl;
    MC_unmatched->GetEntry(nentry);

    if(isSymm && TMath::Abs(eta_2) > (Float_t)etaHigh/10) continue;       
    if(!isSymm && (TMath::Abs(eta_2) < (Float_t)etaLow/10 || TMath::Abs(eta_2) > (Float_t)etaHigh/10)) continue;
    
    Int_t cBin = findBin(hiBin_2);
    if(cBin == -1 || cBin >= nbins_cent) continue;    
    
    Float_t cent_weight = hCentWeight->GetBinContent(hiBin_2);
    Float_t vz_weight = fVz->Eval(vz_2);
    Float_t Sumcand = chSum_2 + phSum_2 + neSum_2 + muSum_2;

    if(subid_2!=0){
      
      if(eMax_2/Sumcand < 0.05 ){
	if(jet55_2==1 && jet65_2==0 && jet80_2==0) hpbpb_fake_Jet55[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet65_2==1 && jet80_2==0) hpbpb_fake_Jet65[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
	if(jet80_2==1) hpbpb_fake_Jet80[cBin]->Fill(pfpt_2, cent_weight*vz_weight);
      }
    }

#if 0
    if(subid_2 != 0) continue;
    if(eMax_2/Sumcand < 0.05  ){
      hpbpb_gen[cBin]->Fill(pfrefpt_2, weight);
      hpbpb_reco[cBin]->Fill(pfpt_2, weight);
      hpbpb_matrix[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
    }
    hpbpb_MC_noCut[cBin]->Fill(pfrefpt_2, weight);
    if(jet55_2 == 1 && jet65_2==0 && jet80_2 == 0){
      hpbpb_MC_Jet55_noCut[cBin]->Fill(pfrefpt_2,weight);
      hMC_unmatched_Jet55_noCut->Fill(pfrefpt_2, weight);
      if(eMax_2/Sumcand < 0.05  ){
	hpbpb_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
     	hMC_unmatched_Jet55_CutA->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet55_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet55_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
      }
      else hMC_unmatched_Jet55_CutA_rej->Fill(pfrefpt_2, weight);
      
    }
    
    if(jet65_2 == 1 && jet80_2 == 0){
      hpbpb_MC_Jet65_noCut[cBin]->Fill(pfrefpt_2,weight);
      hMC_unmatched_Jet65_noCut->Fill(pfrefpt_2, weight);
      if(eMax_2/Sumcand < 0.05  ){
	hpbpb_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hMC_unmatched_Jet65_CutA->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet65_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet65_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
      }
      else hMC_unmatched_Jet65_CutA_rej->Fill(pfrefpt_2);
      
    }
    
    if(jet80_2 == 1){
      hpbpb_MC_Jet80_noCut[cBin]->Fill(pfrefpt_2,weight);
      hMC_unmatched_Jet80_noCut->Fill(pfrefpt_2, weight);
      if(eMax_2/Sumcand < 0.05  ){
	hpbpb_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
	hMC_unmatched_Jet80_CutA->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet80_gen[cBin]->Fill(pfrefpt_2, weight);
	hpbpb_anaBin_Jet80_reco[cBin]->Fill(pfpt_2, weight);
	hpbpb_anaBin_matrix_HLT[cBin]->Fill(pfrefpt_2, pfpt_2, weight);
      }else hMC_unmatched_Jet80_CutA_rej->Fill(pfrefpt_2, weight);
      
    }
#endif
    
  }// mc unmatched  ntuple loop

  TFile fout(Form("../../Output/Pawan_ntuple_PbPb_MC_subidNot0_fake_spectra_JetID_CutA_finebins_%s_R0p%d.root",etaWidth,radius),"RECREATE");
  fout.cd();

  //Int_t TotalEvents = 1973871; 
  // total events per centrality bin; 
  //Float_t TotalEvents[nbins_cent] =      {1.237110e05, 1.081720e05, 4.441970e05, 4.63280e05, 4.85332e05, 5.8576e05};
  Float_t TotalEvents[nbins_cent] =      {6.147599e05, 4.678544e05, 6.17779e05, 7.98961e04, 1.59163e04, 2.366699e03};
  Float_t TotalEvents_Date[nbins_cent] = {9.054905e06, 5.240196e06, 6.678498e06, 8.09324e05, 1.58906e05, 2.4066e04};
  Float_t Lumi = 145.156; //in inverse micro barns 

  for(int i = 0;i<nbins_cent;++i){

    hpbpb_fake_JetComb[i]->Add(hpbpb_fake_Jet80[i]);
    hpbpb_fake_JetComb[i]->Add(hpbpb_fake_Jet65[i]);
    hpbpb_fake_JetComb[i]->Add(hpbpb_fake_Jet55[i]);

    hpbpb_fake_JetComb[i]->Scale(1./TotalEvents[i]);
    //hpbpb_fake_JetComb[i]->Scale(Lumi*1e6);
    hpbpb_fake_Jet80[i]->Scale(1./TotalEvents[i]);
    //hpbpb_fake_Jet80[i]->Scale(Lumi*1e6);
    hpbpb_fake_Jet65[i]->Scale(1./TotalEvents[i]);
    //hpbpb_fake_Jet65[i]->Scale(Lumi*1e6);
    hpbpb_fake_Jet55[i]->Scale(1./TotalEvents[i]);
    //hpbpb_fake_Jet55[i]->Scale(Lumi*1e6);

    hpbpb_MC_Comb_noCut[i]->Add(hpbpb_MC_Jet80_noCut[i]);
    hpbpb_MC_Comb_noCut[i]->Add(hpbpb_MC_Jet65_noCut[i]);
    hpbpb_MC_Comb_noCut[i]->Add(hpbpb_MC_Jet55_noCut[i]);

    divideBinWidth(hpbpb_MC_Comb_noCut[i]);
    divideBinWidth(hpbpb_MC_Jet80_noCut[i]);
    divideBinWidth(hpbpb_MC_Jet65_noCut[i]);
    divideBinWidth(hpbpb_MC_Jet55_noCut[i]);

    hpbpb_Data_Comb_noCut[i]->Add(hpbpb_Data_Jet80_noCut[i]);
    hpbpb_Data_Comb_noCut[i]->Add(hpbpb_Data_Jet65_noCut[i]);
    hpbpb_Data_Comb_noCut[i]->Add(hpbpb_Data_Jet55_noCut[i]);

    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj80[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj65[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj55[i]);

    divideBinWidth(hpbpb_TrgObjComb[i]);
    divideBinWidth(hpbpb_TrgObj80[i]);
    divideBinWidth(hpbpb_TrgObj65[i]);
    divideBinWidth(hpbpb_TrgObj55[i]);

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
    
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet80_data[i]);
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet65_data[i]);
    hpbpb_mcclosure_JetComb_data[i]->Add(hpbpb_mcclosure_Jet55_data[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_data[i]);
    
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet80_gen[i]);
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet65_gen[i]);
    hpbpb_mcclosure_JetComb_gen[i]->Add(hpbpb_mcclosure_Jet55_gen[i]);

    divideBinWidth(hpbpb_mcclosure_JetComb_gen[i]);	
	
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet80_reco[i]);
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet65_reco[i]);
    hpbpb_JetComb_reco[i]->Add(hpbpb_Jet55_reco[i]);
    
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet80_gen[i]);
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet65_gen[i]);
    hpbpb_JetComb_gen[i]->Add(hpbpb_Jet55_gen[i]);

    divideBinWidth(hpbpb_JetComb_gen[i]);
    divideBinWidth(hpbpb_JetComb_reco[i]);
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



  for(int i = 0; i<nbins_cent;++i){

    hpbpb_MC_noCut[i]->Write();
    
    hpbpb_MC_Comb_noCut[i]->Write();
    hpbpb_Data_Comb_noCut[i]->Write();
    
    hpbpb_fake_JetComb[i]->Write();
    hpbpb_fake_Jet80[i]->Write();
    hpbpb_fake_Jet65[i]->Write();
    hpbpb_fake_Jet55[i]->Write();

    hpbpb_TrgObjComb[i]->Write();
    hpbpb_TrgObj80[i]->Write();
    hpbpb_TrgObj65[i]->Write();
    hpbpb_TrgObj55[i]->Write();

    hpbpb_anaBin_TrgObjComb[i]->Write();
    hpbpb_anaBin_TrgObj80[i]->Write();
    hpbpb_anaBin_TrgObj65[i]->Write();
    hpbpb_anaBin_TrgObj55[i]->Write();


    hpbpb_JEC_TrgObjComb[i]->Write();
    hpbpb_JEC_TrgObj80[i]->Write();
    hpbpb_JEC_TrgObj65[i]->Write();
    hpbpb_JEC_TrgObj55[i]->Write();

    hpbpb_Smear_TrgObjComb[i]->Write();
    hpbpb_Smear_TrgObj80[i]->Write();
    hpbpb_Smear_TrgObj65[i]->Write();
    hpbpb_Smear_TrgObj55[i]->Write();
    
    hpbpb_matrix_HLT[i]->Write();
    hpbpb_anaBin_matrix_HLT[i]->Write();

    hpbpb_mcclosure_matrix_HLT[i]->Write();
    hpbpb_mcclosure_JetComb_data[i]->Write();
    hpbpb_mcclosure_Jet80_data[i]->Write();
    hpbpb_mcclosure_Jet65_data[i]->Write();
    hpbpb_mcclosure_Jet55_data[i]->Write();
    hpbpb_mcclosure_JetComb_gen[i]->Write();    
    hpbpb_mcclosure_Jet80_gen[i]->Write();
    hpbpb_mcclosure_Jet65_gen[i]->Write();
    hpbpb_mcclosure_Jet55_gen[i]->Write();

    hpbpb_JetComb_reco[i]->Write();
    hpbpb_Jet80_reco[i]->Write();
    hpbpb_Jet65_reco[i]->Write();
    hpbpb_Jet55_reco[i]->Write();
    hpbpb_JetComb_gen[i]->Write();
    hpbpb_Jet80_gen[i]->Write();
    hpbpb_Jet65_gen[i]->Write();
    hpbpb_Jet55_gen[i]->Write();

    hpbpb_anaBin_JetComb_reco[i]->Write();
    hpbpb_anaBin_Jet80_reco[i]->Write();
    hpbpb_anaBin_Jet65_reco[i]->Write();
    hpbpb_anaBin_Jet55_reco[i]->Write();
    hpbpb_anaBin_JetComb_gen[i]->Write();
    hpbpb_anaBin_Jet80_gen[i]->Write();
    hpbpb_anaBin_Jet65_gen[i]->Write();
    hpbpb_anaBin_Jet55_gen[i]->Write();
    
    hpbpb_reco[i]->Write();
    hpbpb_gen[i]->Write();
    hpbpb_matrix[i]->Write();

  }


#if 0
  // add the unmatched histograms to the matched ones to get the final cut efficiency
  hData_Jet55_noCut->Add(hData_unmatched_Jet55_noCut);
  hData_Jet65_noCut->Add(hData_unmatched_Jet65_noCut);
  hData_Jet80_noCut->Add(hData_unmatched_Jet80_noCut);
  
  hData_Jet55_CutA->Add(hData_unmatched_Jet55_CutA);
  hData_Jet65_CutA->Add(hData_unmatched_Jet65_CutA);
  hData_Jet80_CutA->Add(hData_unmatched_Jet80_CutA);
  hData_Jet55_CutA_rej->Add(hData_unmatched_Jet55_CutA_rej);
  hData_Jet65_CutA_rej->Add(hData_unmatched_Jet65_CutA_rej);
  hData_Jet80_CutA_rej->Add(hData_unmatched_Jet80_CutA_rej);
  hMC_Jet55_noCut->Add(hMC_unmatched_Jet55_noCut);
  hMC_Jet65_noCut->Add(hMC_unmatched_Jet65_noCut);
  hMC_Jet80_noCut->Add(hMC_unmatched_Jet80_noCut);
  
  hMC_Jet55_CutA->Add(hMC_unmatched_Jet55_CutA);
  hMC_Jet65_CutA->Add(hMC_unmatched_Jet65_CutA);
  hMC_Jet80_CutA->Add(hMC_unmatched_Jet80_CutA);
  hMC_Jet55_CutA_rej->Add(hMC_unmatched_Jet55_CutA_rej);
  hMC_Jet65_CutA_rej->Add(hMC_unmatched_Jet65_CutA_rej);
  hMC_Jet80_CutA_rej->Add(hMC_unmatched_Jet80_CutA_rej);
  
  hData_Jet65_noCut->Write();
  hData_Jet65_CutA->Write();
  // hData_Jet65_CutB->Write();
  TH1F * hData_Jet65_CutA_eff = (TH1F*)hData_Jet65_CutA->Clone("hData_Jet65_CutA_eff");
  hData_Jet65_CutA_eff->Divide(hData_Jet65_noCut);
  hData_Jet65_CutA_eff->Write();
  // TH1F * hData_Jet65_CutB_eff = (TH1F*)hData_Jet65_CutB->Clone("hData_Jet65_CutB_eff");
  // hData_Jet65_CutB_eff->Divide(hData_Jet65_noCut);
  // hData_Jet65_CutB_eff->Write();
  hData_Jet65_CutA_rej->Write();
  // hData_Jet65_CutB_rej->Write();
  TH1F * hData_Jet65_CutA_rej_eff = (TH1F*)hData_Jet65_CutA_rej->Clone("hData_Jet65_CutA_rej_eff");
  hData_Jet65_CutA_rej_eff->Divide(hData_Jet65_noCut);
  hData_Jet65_CutA_rej_eff->Write();
  // hData_Jet65_CutB_rej_eff = (TH1F*)hData_Jet65_CutB_rej->Clone("hData_Jet65_CutB_rej_eff");
  // hData_Jet65_CutB_rej_eff->Divide(hData_Jet65_noCut);
  // hData_Jet65_CutB_rej_eff->Write();
  hData_Jet55_noCut->Write();
  hData_Jet55_CutA->Write();
  // hData_Jet55_CutB->Write();
  TH1F * hData_Jet55_CutA_eff = (TH1F*)hData_Jet55_CutA->Clone("hData_Jet55_CutA_eff");
  hData_Jet55_CutA_eff->Divide(hData_Jet55_noCut);
  hData_Jet55_CutA_eff->Write();
  // TH1F * hData_Jet55_CutB_eff = (TH1F*)hData_Jet55_CutB->Clone("hData_Jet55_CutB_eff");
  // hData_Jet55_CutB_eff->Divide(hData_Jet55_noCut);
  // hData_Jet55_CutB_eff->Write();
  hData_Jet55_CutA_rej->Write();
  // hData_Jet55_CutB_rej->Write();
  TH1F * hData_Jet55_CutA_rej_eff = (TH1F*)hData_Jet55_CutA_rej->Clone("hData_Jet55_CutA_rej_eff");
  hData_Jet55_CutA_rej_eff->Divide(hData_Jet55_noCut);
  hData_Jet55_CutA_rej_eff->Write();
  // hData_Jet55_CutB_rej_eff = (TH1F*)hData_Jet55_CutB_rej->Clone("hData_Jet55_CutB_rej_eff");
  // hData_Jet55_CutB_rej_eff->Divide(hData_Jet55_noCut);
  // hData_Jet55_CutB_rej_eff->Write();
  hData_Jet80_noCut->Write();
  hData_Jet80_CutA->Write();
  //hData_Jet80_CutB->Write();
  TH1F * hData_Jet80_CutA_eff = (TH1F*)hData_Jet80_CutA->Clone("hData_Jet80_CutA_eff");
  hData_Jet80_CutA_eff->Divide(hData_Jet80_noCut);
  hData_Jet80_CutA_eff->Write();
  // TH1F * hData_Jet80_CutB_eff = (TH1F*)hData_Jet80_CutB->Clone("hData_Jet80_CutB_eff");
  // hData_Jet80_CutB_eff->Divide(hData_Jet80_noCut);
  // hData_Jet80_CutB_eff->Write();
  hData_Jet80_CutA_rej->Write();
  // hData_Jet80_CutB_rej->Write();
  TH1F * hData_Jet80_CutA_rej_eff = (TH1F*)hData_Jet80_CutA_rej->Clone("hData_Jet80_CutA_rej_eff");
  hData_Jet80_CutA_rej_eff->Divide(hData_Jet80_noCut);
  hData_Jet80_CutA_rej_eff->Write();
  // hData_Jet80_CutB_rej_eff = (TH1F*)hData_Jet80_CutB_rej->Clone("hData_Jet80_CutB_rej_eff");
  // hData_Jet80_CutB_rej_eff->Divide(hData_Jet80_noCut);
  // hData_Jet80_CutB_rej_eff->Write();
  
  
  hMC_Jet80_noCut->Write();
  hMC_Jet80_CutA->Write();
  // hMC_Jet80_CutB->Write();
  TH1F * hMC_Jet80_CutA_eff = (TH1F*)hMC_Jet80_CutA->Clone("hMC_Jet80_CutA_eff");
  hMC_Jet80_CutA_eff->Divide(hMC_Jet80_noCut);
  hMC_Jet80_CutA_eff->Write();
  // TH1F * hMC_Jet80_CutB_eff = (TH1F*)hMC_Jet80_CutB->Clone("hMC_Jet80_CutB_eff");
  // hMC_Jet80_CutB_eff->Divide(hMC_Jet80_noCut);
  // hMC_Jet80_CutB_eff->Write();
  
  hMC_Jet65_noCut->Write();
  hMC_Jet65_CutA->Write();
  // hMC_Jet65_CutB->Write();
  TH1F * hMC_Jet65_CutA_eff = (TH1F*)hMC_Jet65_CutA->Clone("hMC_Jet65_CutA_eff");
  hMC_Jet65_CutA_eff->Divide(hMC_Jet65_noCut);
  hMC_Jet65_CutA_eff->Write();
  // TH1F * hMC_Jet65_CutB_eff = (TH1F*)hMC_Jet65_CutB->Clone("hMC_Jet65_CutB_eff");
  // hMC_Jet65_CutB_eff->Divide(hMC_Jet65_noCut);
  // hMC_Jet65_CutB_eff->Write();
  
  hMC_Jet55_noCut->Write();
  hMC_Jet55_CutA->Write();
  // hMC_Jet55_CutB->Write();
  TH1F * hMC_Jet55_CutA_eff = (TH1F*)hMC_Jet55_CutA->Clone("hMC_Jet55_CutA_eff");
  hMC_Jet55_CutA_eff->Divide(hMC_Jet55_noCut);
  hMC_Jet55_CutA_eff->Write();
  // TH1F * hMC_Jet55_CutB_eff = (TH1F*)hMC_Jet55_CutB->Clone("hMC_Jet55_CutB_eff");
  // hMC_Jet55_CutB_eff->Divide(hMC_Jet55_noCut);
  // hMC_Jet55_CutB_eff->Write();
  // save the unmatched histograms as well:
  hData_unmatched_Jet80_noCut->Write();
  hData_unmatched_Jet80_CutA->Write();
  hData_unmatched_Jet80_CutA_rej->Write();
  hData_unmatched_Jet65_noCut->Write();
  hData_unmatched_Jet65_CutA->Write();
  hData_unmatched_Jet65_CutA_rej->Write();
  hData_unmatched_Jet55_noCut->Write();
  hData_unmatched_Jet55_CutA->Write();
  hData_unmatched_Jet55_CutA_rej->Write();
  hMC_unmatched_Jet80_noCut->Write();
  hMC_unmatched_Jet80_CutA->Write();
  hMC_unmatched_Jet80_CutA_rej->Write();
  hMC_unmatched_Jet65_noCut->Write();
  hMC_unmatched_Jet65_CutA->Write();
  hMC_unmatched_Jet65_CutA_rej->Write();
  hMC_unmatched_Jet55_noCut->Write();
  hMC_unmatched_Jet55_CutA->Write();
  hMC_unmatched_Jet55_CutA_rej->Write();
  
  
  TCanvas * cJet80_CutEfficiency_Jet80 = new TCanvas("cJet80_CutEfficiency_Jet80","",1000,800);
  cJet80_CutEfficiency_Jet80->Divide(2,1);
  cJet80_CutEfficiency_Jet80->cd(1);
  hData_Jet80_CutA_eff->Rebin(5);
  hData_Jet80_CutA_eff->Scale(1./5);
  hData_Jet80_CutA_eff->SetMarkerColor(kRed);
  hData_Jet80_CutA_eff->SetMarkerStyle(24);
  hData_Jet80_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet80_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet80_CutA_eff->SetXTitle(Form("akPu%dPF p_{T}",radius));
  hData_Jet80_CutA_eff->SetTitle("Data");
  hData_Jet80_CutA_eff->SetYTitle("Jet80_Cut efficiency");
  hData_Jet80_CutA_eff->Draw();
  // hData_Jet80_CutB_eff->Rebin(20);
  // hData_Jet80_CutB_eff->Scale(1./20);
  // hData_Jet80_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet80_CutB_eff->SetMarkerStyle(33);
  // hData_Jet80_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet80_CutB_eff->Draw("same");
  cJet80_CutEfficiency_Jet80->cd(2);
  hMC_Jet80_CutA_eff->Rebin(5);
  hMC_Jet80_CutA_eff->Scale(1./5);
  hMC_Jet80_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet80_CutA_eff->SetMarkerStyle(24);
  hMC_Jet80_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet80_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet80_CutA_eff->SetXTitle(Form("akPu%dPF ref p_{T}",radius));
  hMC_Jet80_CutA_eff->SetTitle("MC");
  hMC_Jet80_CutA_eff->SetYTitle("Jet80_Cut efficiency");
  hMC_Jet80_CutA_eff->Draw();
  // hMC_Jet80_CutB_eff->Rebin(20);
  // hMC_Jet80_CutB_eff->Scale(1./20);
  // hMC_Jet80_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet80_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet80_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet80_CutB_eff->Draw("same");
  cJet80_CutEfficiency_Jet80->SaveAs(Form("PbPb_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_R0p%d_zoomed.pdf",radius),"RECREATE");
  // TCanvas * cCutRejection_Jet80 = new TCanvas("cCutRejection_Jet80","",1000,800);
  // hData_Jet80_CutA_rej->Rebin(5);
  // hData_Jet80_CutA_rej->Scale(1./5);
  // hData_Jet80_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet80_CutA_rej->SetMarkerStyle(24);
  // hData_Jet80_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet80_CutA_rej->SetAxisRange(0,1.2,"Y");
  // hData_Jet80_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  // hData_Jet80_CutA_rej->SetTitle("Data");
  // hData_Jet80_CutA_rej->SetYTitle("Jet80_Cut Rejection");
  // hData_Jet80_CutA_rej->Draw();
  // // hData_Jet80_CutB_rej->Rebin(20);
  // // hData_Jet80_CutB_rej->Scale(1./20);
  // // hData_Jet80_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet80_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet80_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet80_CutB_rej->Draw("same");
  // cCutRejection_Jet80->SaveAs("PbPb_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection.pdf","RECREATE");
  TCanvas * cJet55_CutEfficiency_Jet55 = new TCanvas("cJet55_CutEfficiency_Jet55","",1000,800);
  cJet55_CutEfficiency_Jet55->Divide(2,1);
  cJet55_CutEfficiency_Jet55->cd(1);
  hData_Jet55_CutA_eff->Rebin(5);
  hData_Jet55_CutA_eff->Scale(1./5);
  hData_Jet55_CutA_eff->SetMarkerColor(kRed);
  hData_Jet55_CutA_eff->SetMarkerStyle(24);
  hData_Jet55_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet55_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet55_CutA_eff->SetXTitle(Form("akPu%dPF p_{T}",radius));
  hData_Jet55_CutA_eff->SetTitle("Data");
  hData_Jet55_CutA_eff->SetYTitle("Jet55_Cut efficiency");
  hData_Jet55_CutA_eff->Draw();
  // hData_Jet55_CutB_eff->Rebin(20);
  // hData_Jet55_CutB_eff->Scale(1./20);
  // hData_Jet55_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet55_CutB_eff->SetMarkerStyle(33);
  // hData_Jet55_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet55_CutB_eff->SetAxisRange(0,1.2,"Y");
  // hData_Jet55_CutB_eff->Draw("same");
  cJet55_CutEfficiency_Jet55->cd(2);
  hMC_Jet55_CutA_eff->Rebin(5);
  hMC_Jet55_CutA_eff->Scale(1./5);
  hMC_Jet55_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet55_CutA_eff->SetMarkerStyle(24);
  hMC_Jet55_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet55_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet55_CutA_eff->SetXTitle(Form("akPu%dPF ref p_{T}",radius));
  hMC_Jet55_CutA_eff->SetTitle("MC");
  hMC_Jet55_CutA_eff->SetYTitle("Jet55_Cut efficiency");
  hMC_Jet55_CutA_eff->Draw();
  // hMC_Jet55_CutB_eff->Rebin(20);
  // hMC_Jet55_CutB_eff->Scale(1./20);
  // hMC_Jet55_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet55_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet55_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet55_CutB_eff->Draw("same");
  cJet55_CutEfficiency_Jet55->SaveAs(Form("PbPb_YetkinCuts_Jet55_noJet65noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched_R0p%d_zoomed.pdf",radius),"RECREATE");
  // TCanvas * cCutRejection_Jet55 = new TCanvas("cCutRejection_Jet55","",1000,800);
  // hData_Jet55_CutA_rej->Rebin(5);
  // hData_Jet55_CutA_rej->Scale(1./5);
  // hData_Jet55_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet55_CutA_rej->SetMarkerStyle(24);
  // hData_Jet55_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet55_CutA_rej->SetAxisRange(0,1.2,"Y");
  // hData_Jet55_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  // hData_Jet55_CutA_rej->SetTitle("Data");
  // hData_Jet55_CutA_rej->SetYTitle("Jet55_Cut Rejection");
  // hData_Jet55_CutA_rej->Draw();
  // // hData_Jet55_CutB_rej->Rebin(20);
  // // hData_Jet55_CutB_rej->Scale(1./20);
  // // hData_Jet55_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet55_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet55_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet55_CutB_rej->Draw("same");
  // cCutRejection_Jet55->SaveAs("PbPb_YetkinCuts_Jet55_noJet65noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");
  TCanvas * cJet65_CutEfficiency_Jet65 = new TCanvas("cJet65_CutEfficiency_Jet65","",1000,800);
  cJet65_CutEfficiency_Jet65->Divide(2,1);
  cJet65_CutEfficiency_Jet65->cd(1);
  hData_Jet65_CutA_eff->Rebin(5);
  hData_Jet65_CutA_eff->Scale(1./5);
  hData_Jet65_CutA_eff->SetMarkerColor(kRed);
  hData_Jet65_CutA_eff->SetMarkerStyle(24);
  hData_Jet65_CutA_eff->SetAxisRange(30,300,"X");
  hData_Jet65_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hData_Jet65_CutA_eff->SetXTitle(Form("matched akPu%dPF p_{T}",radius));
  hData_Jet65_CutA_eff->SetTitle("Data");
  hData_Jet65_CutA_eff->SetYTitle("Jet65_Cut efficiency");
  hData_Jet65_CutA_eff->Draw();
  // hData_Jet65_CutB_eff->Rebin(20);
  // hData_Jet65_CutB_eff->Scale(1./20);
  // hData_Jet65_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet65_CutB_eff->SetMarkerStyle(33);
  // hData_Jet65_CutB_eff->SetAxisRange(30,300,"X");
  // hData_Jet65_CutB_eff->Draw("same");
  cJet65_CutEfficiency_Jet65->cd(2);
  hMC_Jet65_CutA_eff->Rebin(5);
  hMC_Jet65_CutA_eff->Scale(1./5);
  hMC_Jet65_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet65_CutA_eff->SetMarkerStyle(24);
  hMC_Jet65_CutA_eff->SetAxisRange(30,300,"X");
  hMC_Jet65_CutA_eff->SetAxisRange(0.85,1.1,"Y");
  hMC_Jet65_CutA_eff->SetXTitle(Form("matched akPu%dPF ref p_{T}",radius));
  hMC_Jet65_CutA_eff->SetTitle("MC");
  hMC_Jet65_CutA_eff->SetYTitle("Jet65_Cut efficiency");
  hMC_Jet65_CutA_eff->Draw();
  // hMC_Jet65_CutB_eff->Rebin(20);
  // hMC_Jet65_CutB_eff->Scale(1./20);
  // hMC_Jet65_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet65_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet65_CutB_eff->SetAxisRange(30,300,"X");
  // hMC_Jet65_CutB_eff->Draw("same");
  cJet65_CutEfficiency_Jet65->SaveAs(Form("PbPb_YetkinCuts_Jet65_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched_R0p%d_zoomed.pdf",radius),"RECREATE");
  // TCanvas * cCutRejection_Jet65 = new TCanvas("cCutRejection_Jet65","",1000,800);
  // hData_Jet65_CutA_rej->Rebin(5);
  // hData_Jet65_CutA_rej->Scale(1./5);
  // hData_Jet65_CutA_rej->SetMarkerColor(kRed);
  // hData_Jet65_CutA_rej->SetMarkerStyle(24);
  // hData_Jet65_CutA_rej->SetAxisRange(30,300,"X");
  // hData_Jet65_CutA_rej->SetAxisRange(0, 1.2,"Y");
  // hData_Jet65_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  // hData_Jet65_CutA_rej->SetTitle("Data");
  // hData_Jet65_CutA_rej->SetYTitle("Jet65_Cut Rejection");
  // hData_Jet65_CutA_rej->Draw();
  // // hData_Jet65_CutB_rej->Rebin(20);
  // // hData_Jet65_CutB_rej->Scale(1./20);
  // // hData_Jet65_CutB_rej->SetMarkerColor(kBlack);
  // // hData_Jet65_CutB_rej->SetMarkerStyle(33);
  // // hData_Jet65_CutB_rej->SetAxisRange(30,300,"X");
  // // hData_Jet65_CutB_rej->Draw("same");
  // cCutRejection_Jet65->SaveAs("PbPb_YetkinCuts_Jet65_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");
  // plot the trigger combination from this, and the total cut efficiency:
  
  TH1F * hData_Jet80 = (TH1F*)hData_Jet80_CutA->Clone("hData_Jet80");
  TH1F * hData_Jet65 = (TH1F*)hData_Jet65_CutA->Clone("hData_Jet65");
  TH1F * hData_Jet55 = (TH1F*)hData_Jet55_CutA->Clone("hData_Jet55");
  TH1F * hData_Combined = (TH1F*)hData_Jet80->Clone("hData_Combined");
  hData_Combined->Add(hData_Jet65);
  hData_Combined->Add(hData_Jet55);
  hData_Combined->Print("base");
  
  TH1F * hData_noCut_Jet80 = (TH1F*)hData_Jet80_noCut->Clone("hData_noCut_Jet80");
  TH1F * hData_noCut_Jet65 = (TH1F*)hData_Jet65_noCut->Clone("hData_noCut_Jet65");
  TH1F * hData_noCut_Jet55 = (TH1F*)hData_Jet55_noCut->Clone("hData_noCut_Jet55");
  TH1F * hData_noCut_Combined = (TH1F*)hData_noCut_Jet80->Clone("hData_noCut_Combined");
  hData_noCut_Combined->Add(hData_noCut_Jet65);
  hData_noCut_Combined->Add(hData_noCut_Jet55);
  
  hData_noCut_Combined->Print("base");
  
  TH1F * hMC_Jet80 = (TH1F*)hMC_Jet80_CutA->Clone("hMC_Jet80");
  TH1F * hMC_Jet65 = (TH1F*)hMC_Jet65_CutA->Clone("hMC_Jet65");
  TH1F * hMC_Jet55 = (TH1F*)hMC_Jet55_CutA->Clone("hMC_Jet55");
  TH1F * hMC_Combined = (TH1F*)hMC_Jet80->Clone("hMC_Combined");
  hMC_Combined->Add(hMC_Jet65);
  hMC_Combined->Add(hMC_Jet55);
  hMC_Combined->Print("base");
  
  TH1F * hMC_noCut_Jet80 = (TH1F*)hMC_Jet80_noCut->Clone("hMC_noCut_Jet80");
  TH1F * hMC_noCut_Jet65 = (TH1F*)hMC_Jet65_noCut->Clone("hMC_noCut_Jet65");
  TH1F * hMC_noCut_Jet55 = (TH1F*)hMC_Jet55_noCut->Clone("hMC_noCut_Jet55");
  TH1F * hMC_noCut_Combined = (TH1F*)hMC_noCut_Jet80->Clone("hMC_noCut_Combined");
  hMC_noCut_Combined->Add(hMC_noCut_Jet65);
  hMC_noCut_Combined->Add(hMC_noCut_Jet55);
  hMC_noCut_Combined->Print("base");
  
  TH1F * hData_Combined_Efficiency = (TH1F*)hData_Combined->Clone("hData_Combined_Efficiency");
  hData_Combined_Efficiency->Divide(hData_noCut_Combined);
  hData_Combined_Efficiency->Print("base");
  
  TH1F * hMC_Combined_Efficiency = (TH1F*)hMC_Combined->Clone("hMC_Combined_Efficiency");
  hMC_Combined_Efficiency->Divide(hMC_noCut_Combined);
  hMC_Combined_Efficiency->Print("base");
  
  TCanvas * cCombinedEff = new TCanvas("cCombinedEff","",800,600);
  hData_Combined_Efficiency->SetXTitle("Jet p_{T}");
  hData_Combined_Efficiency->SetYTitle("Combined Jet ID cut efficiency");
  hData_Combined_Efficiency->SetMarkerStyle(20);
  hData_Combined_Efficiency->SetMarkerColor(kBlack);
  hData_Combined_Efficiency->Rebin(10);
  hData_Combined_Efficiency->Scale(1/10);
  hData_Combined_Efficiency->SetAxisRange(30,350,"X");
  hData_Combined_Efficiency->Draw();
  
  hMC_Combined_Efficiency->SetMarkerStyle(24);
  hMC_Combined_Efficiency->SetMarkerColor(kRed);
  hMC_Combined_Efficiency->Rebin(10);
  hMC_Combined_Efficiency->Scale(1/10);
  hMC_Combined_Efficiency->Draw("same");
  cCombinedEff->SaveAs(Form("Combined_trigger_efficiency_YetkinCut_R0p%d.pdf",radius),"RECREATE");
  
  TCanvas * cTriggerCombination = new TCanvas("cTriggerCombination","",800,600);
  cTriggerCombination->SetLogy();
  hData_Combined->SetMarkerColor(kBlack);
  hData_Combined->SetMarkerStyle(25);
  hData_Combined->SetAxisRange(30,350,"X");
  //hData_Combined->SetXTitle("");
  //hData_Combined->Draw();
  hData_Jet80->SetMarkerColor(kRed);
  hData_Jet80->SetMarkerStyle(20);
  hData_Jet80->Draw("same");
  
  hData_Jet65->SetMarkerColor(kBlue);
  hData_Jet65->SetMarkerStyle(20);
  hData_Jet65->Draw("same");
  hData_Jet55->SetMarkerColor(kGreen);
  hData_Jet55->SetMarkerStyle(20);
  hData_Jet55->Draw("same");
  //drawText()
  cTriggerCombination->SaveAs(Form("TriggerCombination_YetkinCuts_R0p%d.pdf",radius),"RECREATE");
#endif

}
