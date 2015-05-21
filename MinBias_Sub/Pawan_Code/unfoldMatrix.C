
#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TMatrixD.h>

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <utility>

#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldResponse.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldBayes.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldSvd.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/RooUnfoldBinByBin.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/prior.h"
#include "/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/Headers/bayesianUnfold.h"


using namespace std;

#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#endif

//! constants
#define iYear 2015

#define pi 3.14159265

#define  ketacut    2.0
#define  kptrawcut  0.0
#define  kptrecocut 0.0
#define  kdelrmatch 0.2
#define  kdelrcut   0.3
#define  kvzcut     15.0

const int nIter=5;

void AddInputFiles(TChain */*ch*/, string /*iname*/, string /*inputTree*/);

int GetPhiBin(float /*phi*/);
int GetEtaBin(float /*eta*/);
int GetPtBin(float /*pt*/);
int GetCentBin(int /*hiBin*/);
float delphi(float /*phi1*/, float /*phi2*/);
float deltaR(float /*eta1*/, float /*phi1*/, 
	      float /*eta2*/, float /*phi2*/);

double GetXsec(double /*maxpthat*/);
void GetCentWeight(TH1F */*hCentWeight*/);
TH2F* CorrelationHist (const TMatrixD& /*cov*/, const char* /*name*/, const char* /*title*/,
		       const unsigned /*nbins_rec*/, const Double_t */*bins*/);

struct Jet{
  int id;
  float pt;
  float eta;
  float phi;
};

bool compare_pt(Jet jet1, Jet jet2);
bool compare_pt(Jet jet1, Jet jet2){
  return jet1.pt > jet2.pt;
}

typedef std::pair< Jet, Jet > CaloPFJetPair;
struct CompareMatchedJets {
  //! Calo-PF match
  bool operator()(const CaloPFJetPair &A1, const CaloPFJetPair &A2){
    
    Jet cj1 = A1.first;  //! CaloJet 1st pair
    Jet pf1 = A1.second; //! PFJet   1st pair
    
    Jet cj2 = A2.first;  //! CaloJet 2nd pair
    Jet pf2 = A2.second; //! PFJet   2nd pair

    float delr1 = deltaR(cj1.eta, cj1.phi, pf1.eta, pf1.phi);
    float delr2 = deltaR(cj2.eta, cj2.phi, pf2.eta, pf2.phi);

    //float delpt1 = fabs(cj1.pt - pf1.pt);
    //float delpt2 = fabs(cj2.pt - pf2.pt);
    
    return ((delr1 < delr2) && (cj1.pt > cj2.pt));
  }
};

typedef std::multiset< CaloPFJetPair, CompareMatchedJets > CaloPFMatchedJets;
typedef std::multiset< CaloPFJetPair >::iterator CPFItr;
//typedef std::multiset< Jet >::value_type MatchedJet;



const int ncen=8;
const char *cdir [ncen] = {"05","510","1030","3050","5070","7090","90100","pp"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","90-100%","pp"};

const int ncand=5;
const char *ccand[ncand] = {"h^{#pm}","e^{#pm}","#mu^{#pm}","#gamma","h0"};


//! pt binning
const double ptbins[] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 300, 
  330, 362, 395,
  430, 468, 507, 
  548
}; 
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

const double etabins[] = {-2.000, -1.4000, -0.4500, 0.000, 0.4500, 1.400, 2.000};
const int neta = sizeof(etabins)/sizeof(double) - 1;
const double phibins[] = {-3.141,-2.100,-1.500,-0.800,-0.300, 
 			  0.300,0.800,1.500,2.100,3.141
};
const int nphi = sizeof(phibins)/sizeof(double) - 1;

double xsec[12][3] ={{2.034e-01,15.  ,30.},   //! 15     0   
		     {1.075e-02,30.  ,50.},   //! 30     1   
		     {1.025e-03,50.  ,80.},   //! 50     2   
		     {9.865e-05,80.  ,120.},  //! 80     3  
		     {1.129e-05,120. ,170.},  //! 120    4  
		     {1.465e-06,170. ,220.},  //! 170    5  
		     {2.837e-07,220. ,280.},  //! 220    6  
		     {5.323e-08,280. ,370.},  //! 280    7  
		     {5.934e-09,370. ,460.},  //! 370    8 
		     {8.125e-10,460. ,540.},  //! 460    9
		     {1.467e-10,540. ,9999.}, //! 540    10
		     {0.0000000,9999.,9999.}  //         11 
};


//! PuPtMin values 
// akPu1PFJets.puPtMin = 10
// akPu2PFJets.puPtMin = 10
// akPu3PFJets.puPtMin = 15
// akPu4PFJets.puPtMin = 20
// akPu5PFJets.puPtMin = 25
// akPu6PFJets.puPtMin = 30
// akPu6PFJets.puPtMin = 35


TStopwatch timer;

int unfoldMatrix(std::string kAlgName="akPu4")
{
  
  timer.Start();

  gSystem->Load("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Headers/RooUnfold-1.1.1/libRooUnfold.so");

  std::cout << std::endl;
  bool printDebug=false;

  std::string kSpecies="pbpb";
  std::string kDataset="mc";
  std::string kFileList="";
  std::string kFoname="Outputhist_unfold_"+kAlgName+"PF.root"; 
  double kMaxpthat=0.0;
		   
  double lumi_scale=1.0;
  // if( kSpecies == "pbpb" ){
  //   lumi_scale = 1./(145.156 * 1e6);
  // } else {
  //   lumi_scale = 1./(5.3 * 1e9);
  // }
  
  int iterations=4;
  int mincen=0;
  int maxcen=ncen-1;
  if( kSpecies == "pp" ){
    mincen=ncen-1;
    maxcen=ncen;
  }


  //! Vertex re-weighting 
  TF1 *fVz=0;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);
  
  //! Centrality re-weighting 
  TH1F *hCentWeight = new TH1F("hCentWeight","Centrality weight",200,0,200);
  GetCentWeight(hCentWeight);

  //! Create output file
  std::string outdir="";
  std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");

  //! 
  //! Define histograms here
  //fout->mkdir(Form("%sJetAnalyzer",kAlgName.c_str()));
  //fout->cd(Form("%sJetAnalyzer"   ,kAlgName.c_str()));

  TH1F *hEvents_Total = new TH1F("hEvents_Total","Total # of events ",10,0.,1.);
  hEvents_Total->Sumw2();
  TH1F *hEvents_pCollEvent = new TH1F("hEvents_pCollEvent","# of events ",10,0.,1.);
  hEvents_pCollEvent->Sumw2();
  TH1F *hEvents_pHBHENoise = new TH1F("hEvents_pHBHENoise","# of events ",10,0.,1.);
  hEvents_pHBHENoise->Sumw2();
  TH1F *hEvents_Vzcut = new TH1F("hEvents_Vzcut","# of events ",10,0.,1.);
  hEvents_Vzcut->Sumw2();
  TH1F *hEvents_Cent = new TH1F("hEvents_Cent","Cent # of events ",10,0.,1.);
  hEvents_Cent->Sumw2();
  TH1F *hEvents_bad = new TH1F("hEvents_bad","bad # of events ",10,0.,1.);
  hEvents_bad->Sumw2();
  TH1F *hEvents_supernova = new TH1F("hEvents_supernova","supernova # of events ",10,0.,1.);
  hEvents_supernova->Sumw2();
  TH1F *hEvents_nopfcalo = new TH1F("hEvents_nopfcalo","nopfcalo # of events ",10,0.,1.);
  hEvents_nopfcalo->Sumw2();
  TH1F *hEvents_maxpthat = new TH1F("hEvents_maxpthat","maxpthat # of events ",10,0.,1.);
  hEvents_maxpthat->Sumw2();

  TH1F *hEvents_jet40 = new TH1F("hEvents_jet40","# of events jet40",10,0.,1.);
  hEvents_jet40->Sumw2();
  TH1F *hEvents_jet60 = new TH1F("hEvents_jet60","# of events jet60",10,0.,1.);
  hEvents_jet60->Sumw2();
  TH1F *hEvents_jet55 = new TH1F("hEvents_jet55","# of events jet55",10,0.,1.);
  hEvents_jet55->Sumw2();
  TH1F *hEvents_jet65 = new TH1F("hEvents_jet65","# of events jet65",10,0.,1.);
  hEvents_jet65->Sumw2();
  TH1F *hEvents_jet80 = new TH1F("hEvents_jet80","# of events jet80",10,0.,1.);
  hEvents_jet80->Sumw2();

  TH1F *hEvents_jet40_nojet60_nojet80 = new TH1F("hEvents_jet40_nojet60_nojet80","# of events jet40 && !jet60 && !jet80",10,0.,1.);
  hEvents_jet40_nojet60_nojet80->Sumw2();
  TH1F *hEvents_jet60_nojet80 = new TH1F("hEvents_jet60_nojet80","# of events jet60 && !jet80",10,0.,1.);
  hEvents_jet60_nojet80->Sumw2();

  TH2F *hEvents_jet40_prescl = new TH2F("hEvents_jet40_prescl","prescaled # of events jet40",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet40_prescl->Sumw2();
  TH2F *hEvents_jet60_prescl = new TH2F("hEvents_jet60_prescl","prescaled # of events jet60",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet60_prescl->Sumw2();				
  TH2F *hEvents_jet55_prescl = new TH2F("hEvents_jet55_prescl","prescaled # of events jet55",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet55_prescl->Sumw2();				
  TH2F *hEvents_jet65_prescl = new TH2F("hEvents_jet65_prescl","prescaled # of events jet65",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet65_prescl->Sumw2();				
  TH2F *hEvents_jet80_prescl = new TH2F("hEvents_jet80_prescl","prescaled # of events jet80",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet80_prescl->Sumw2();
  
  TH1F *hEvents_jet55_nojet65_nojet80 = new TH1F("hEvents_jet55_nojet65_nojet80","# of events jet55 && !jet65 && !jet80",10,0.,1.);
  hEvents_jet55_nojet65_nojet80->Sumw2();
  TH1F *hEvents_jet65_nojet80 = new TH1F("hEvents_jet65_nojet80","# of events jet60 && !jet80",10,0.,1.);
  hEvents_jet65_nojet80->Sumw2();

  RooUnfoldBayes *unf_bayes=0;

  //! 2D histograms for unfolding
  TH2F *hmatrix[2][ncen], *hmatrix_m[2][ncen], *hmatrix_um[2][ncen];
  TH1F *hgen   [2][ncen], *hgen_m   [2][ncen], *hgen_um   [2][ncen];
  TH1F *hrec   [2][ncen], *hrec_m   [2][ncen], *hrec_um   [2][ncen];
  TH1F *hunf   [2][ncen], *hunf_m   [2][ncen], *hunf_um   [2][ncen];
  RooUnfoldResponse *hresp[2][ncen], *hresp_m[2][ncen], *hresp_um[2][ncen];
  TH2F *hresponse[2][ncen], *hresponse_m[2][ncen], *hresponse_um[2][ncen];
  TH2F *hcov[2][ncen], *hcov_m[2][ncen], *hcov_um[2][ncen];

  //! fine bin in pT
  TH2F *hmatrix_f[2][ncen], *hmatrix_m_f[2][ncen], *hmatrix_um_f[2][ncen];
  TH1F *hgen_f   [2][ncen], *hgen_m_f   [2][ncen], *hgen_um_f   [2][ncen];
  TH1F *hrec_f   [2][ncen], *hrec_m_f   [2][ncen], *hrec_um_f   [2][ncen];
  TH1F *hunf_f   [2][ncen], *hunf_m_f   [2][ncen], *hunf_um_f   [2][ncen];
  RooUnfoldResponse *hresp_f[2][ncen], *hresp_m_f[2][ncen], *hresp_um_f[2][ncen];
  TH2F *hresponse_f[2][ncen], *hresponse_m_f[2][ncen], *hresponse_um_f[2][ncen];
  TH2F *hcov_f[2][ncen], *hcov_m_f[2][ncen], *hcov_um_f[2][ncen];

  ///! check
  TH2F *hmatrix_c[2][ncen], *hmatrix_m_c[2][ncen], *hmatrix_um_c[2][ncen];
  TH1F *hgen_c[2][ncen], *hgen_m_c[2][ncen], *hgen_um_c[2][ncen];
  TH1F *hrec_c[2][ncen], *hrec_m_c[2][ncen], *hrec_um_c[2][ncen];
  TH1F *hunf_c[2][ncen], *hunf_m_c[2][ncen], *hunf_um_c[2][ncen];

  TH2F *hmatrix_c_f[2][ncen], *hmatrix_m_c_f[2][ncen], *hmatrix_um_c_f[2][ncen];
  TH1F *hgen_c_f [2][ncen], *hgen_m_c_f [2][ncen], *hgen_um_c_f [2][ncen];
  TH1F *hrec_c_f [2][ncen], *hrec_m_c_f [2][ncen], *hrec_um_c_f [2][ncen];
  TH1F *hunf_c_f [2][ncen], *hunf_m_c_f [2][ncen], *hunf_um_c_f [2][ncen];



  //! 0 : wo jetid; 1 : wjetid
  for(int i=0; i<2; i++){
    for(int ic=0; ic<ncen; ic++){

      hmatrix[i][ic] = new TH2F(Form("hmatrix_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("hmatrix_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			     nbins,ptbins,nbins,ptbins);
      hmatrix[i][ic]->Sumw2();
      hmatrix_m[i][ic] = new TH2F(Form("hmatrix_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("matxhed hmatrix_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			       nbins,ptbins,nbins,ptbins);
      hmatrix_m[i][ic]->Sumw2();
      hmatrix_um[i][ic] = new TH2F(Form("hmatrix_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("un-matched hmatrix_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				nbins,ptbins,nbins,ptbins);
      hmatrix_um[i][ic]->Sumw2();

      hmatrix_f[i][ic] = new TH2F(Form("hmatrix_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("hmatrix_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			       600,0.,600.,600,0.,600.);
      hmatrix_f[i][ic]->Sumw2();
      hmatrix_m_f[i][ic] = new TH2F(Form("hmatrix_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("matxhed hmatrix_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				 600,0.,600.,600,0.,600.);
      hmatrix_m_f[i][ic]->Sumw2();
      hmatrix_um_f[i][ic] = new TH2F(Form("hmatrix_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("un-matched hmatrix_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				  600,0.,600.,600,0.,600.);
      hmatrix_um_f[i][ic]->Sumw2();


      ///
      hmatrix_c[i][ic] = new TH2F(Form("hmatrix_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("hmatrix_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			     nbins,ptbins,nbins,ptbins);
      hmatrix_c[i][ic]->Sumw2();
      hmatrix_m_c[i][ic] = new TH2F(Form("hmatrix_m_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("matxhed hmatrix_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			       nbins,ptbins,nbins,ptbins);
      hmatrix_m_c[i][ic]->Sumw2();
      hmatrix_um_c[i][ic] = new TH2F(Form("hmatrix_um_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("un-matched hmatrix_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				nbins,ptbins,nbins,ptbins);
      hmatrix_um_c[i][ic]->Sumw2();

      hmatrix_c_f[i][ic] = new TH2F(Form("hmatrix_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("hmatrix_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			       600,0.,600.,600,0.,600.);
      hmatrix_c_f[i][ic]->Sumw2();
      hmatrix_m_c_f[i][ic] = new TH2F(Form("hmatrix_m_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("matxhed hmatrix_m_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				 600,0.,600.,600,0.,600.);
      hmatrix_m_c_f[i][ic]->Sumw2();
      hmatrix_um_c_f[i][ic] = new TH2F(Form("hmatrix_um_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic), Form("un-matched hmatrix_um_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				  600,0.,600.,600,0.,600.);
      hmatrix_um_c_f[i][ic]->Sumw2();


      hgen[i][ic] = new TH1F(Form("hgen_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("hgen_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			  nbins,ptbins);
      hgen[i][ic]->Sumw2();
      hgen_m[i][ic] = new TH1F(Form("hgen_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched hgen_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			    nbins,ptbins);
      hgen_m[i][ic]->Sumw2();
      hgen_um[i][ic] = new TH1F(Form("hgen_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("un-matched hgen_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			     nbins,ptbins);
      hgen_um[i][ic]->Sumw2();
      hrec[i][ic] = new TH1F(Form("hrec_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("hrec_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			  nbins,ptbins);
      hrec[i][ic]->Sumw2();
      hrec_m[i][ic] = new TH1F(Form("hrec_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched hrec_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			    nbins,ptbins);
      hrec_m[i][ic]->Sumw2();
      hrec_um[i][ic] = new TH1F(Form("hrec_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("unmatched hrec_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			     nbins,ptbins);
      hrec_um[i][ic]->Sumw2();
      

      hgen_c[i][ic] = new TH1F(Form("hgen_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("hgen_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			  nbins,ptbins);
      hgen_c[i][ic]->Sumw2();
      hgen_m_c[i][ic] = new TH1F(Form("hgen_m_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched hgen_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			    nbins,ptbins);
      hgen_m_c[i][ic]->Sumw2();
      hgen_um_c[i][ic] = new TH1F(Form("hgen_um_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("un-matched hgen_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			     nbins,ptbins);
      hgen_um_c[i][ic]->Sumw2();

      hrec_c[i][ic] = new TH1F(Form("hrec_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("check hrec_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			    nbins,ptbins);
      hrec_c[i][ic]->Sumw2();
      
      hrec_m_c[i][ic] = new TH1F(Form("hrec_m_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched check hrec_m_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			      nbins,ptbins);
      hrec_m_c[i][ic]->Sumw2();
      hrec_um_c[i][ic] = new TH1F(Form("hrec_um_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("un matched check hrec_um_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			       nbins,ptbins);
      hrec_um_c[i][ic]->Sumw2();
      
      
      ////// _f fine bin in pT
      
      hgen_f[i][ic] = new TH1F(Form("hgen_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("hgen_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			    600,0.,600.);
      hgen_f[i][ic]->Sumw2();
      hgen_m_f[i][ic] = new TH1F(Form("hgen_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched hgen_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			      600,0.,600.);
      hgen_m_f[i][ic]->Sumw2();
      hgen_um_f[i][ic] = new TH1F(Form("hgen_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("un-matched hgen_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			       600,0.,600.);
      hgen_um_f[i][ic]->Sumw2();
      hrec_f[i][ic] = new TH1F(Form("hrec_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("hrec_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			    600,0.,600.);
      hrec_f[i][ic]->Sumw2();
      hrec_m_f[i][ic] = new TH1F(Form("hrec_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched hrec_m_f%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			      600,0.,600.);
      hrec_m_f[i][ic]->Sumw2();
      hrec_um_f[i][ic] = new TH1F(Form("hrec_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("unmatched hrec_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			       600,0.,600.);
      hrec_um_f[i][ic]->Sumw2();
      
      hrec_c_f[i][ic] = new TH1F(Form("hrec_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("check hrec_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			      600,0.,600.);
      hrec_c_f[i][ic]->Sumw2();
      
      hrec_m_c_f[i][ic] = new TH1F(Form("hrec_m_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched check hrec_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				600,0.,600.);
      hrec_m_c_f[i][ic]->Sumw2();
      hrec_um_c_f[i][ic] = new TH1F(Form("hrec_um_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("un matched check hrec_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				 600,0.,600.);
      hrec_um_c_f[i][ic]->Sumw2();

      hgen_c_f[i][ic] = new TH1F(Form("hgen_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("hgen_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				 600,0.,600.);
      hgen_c_f[i][ic]->Sumw2();
      hgen_m_c_f[i][ic] = new TH1F(Form("hgen_m_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("matched hgen_c_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
				   600,0.,600.);
      hgen_m_c_f[i][ic]->Sumw2();
      hgen_um_c_f[i][ic] = new TH1F(Form("hgen_um_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),Form("un-matched hgen_c_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(), i,ic),
			     nbins,ptbins);
      hgen_um_c_f[i][ic]->Sumw2();

      //! RooUnfold matrix
      hresp[i][ic] = new RooUnfoldResponse(hrec[i][ic], hgen[i][ic], Form("hresp_%s_%s_%d_%d", kSpecies.c_str(), kAlgName.c_str(), i, ic));
      hresp_m[i][ic] = new RooUnfoldResponse(hrec_m[i][ic], hgen_m[i][ic], Form("hresp_m_%s_%s_%d_%d", kSpecies.c_str(), kAlgName.c_str(), i, ic));
      hresp_um[i][ic] = new RooUnfoldResponse(hrec_um[i][ic], hgen_um[i][ic], Form("hresp_um_%s_%s_%d_%d", kSpecies.c_str(), kAlgName.c_str(), i, ic));

      hresp_f[i][ic] = new RooUnfoldResponse(hrec_f[i][ic], hgen_f[i][ic], Form("hresp_f_%s_%s_%d_%d", kSpecies.c_str(), kAlgName.c_str(), i, ic));
      hresp_m_f[i][ic] = new RooUnfoldResponse(hrec_m_f[i][ic], hgen_m_f[i][ic], Form("hresp_m_f_%s_%s_%d_%d", kSpecies.c_str(), kAlgName.c_str(), i, ic));
      hresp_um_f[i][ic] = new RooUnfoldResponse(hrec_um_f[i][ic], hgen_um_f[i][ic], Form("hresp_um_f_%s_%s_%d_%d", kSpecies.c_str(), kAlgName.c_str(), i, ic));
    }
  }    
  fout->cd("../");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"Initialized the histograms " <<std::endl;

  const int kFiles=9;
  std::string inputFile[kFiles]{
    "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat15_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat30_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat50_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat120_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat170_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat280_Track9_Jet30_matchEqR_merged_forest_0.root",
      "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat370_Track9_Jet30_matchEqR_merged_forest_0.root" };

  const double maxpthat[kFiles]={30.0,50.0,80.0,120.0,170.0,220.0,280.0,370.0,9999};

  if( kSpecies == "pbpb" && kDataset == "mc" ){
    xsec[8][2]=9999; //! 370 in PbPb
    xsec[9][0]=0.00000;  xsec[9][1]=9999;  xsec[9][2]=9999;
  }

  
  
  for(int iFile=0; iFile<kFiles; iFile++){
    kFileList = inputFile[iFile];
    kMaxpthat = maxpthat [iFile];
    
    //tr_jet = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
    TChain *tch_pfjet = new TChain(Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
    AddInputFiles(tch_pfjet,kFileList,Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
    cout <<" # of events in PFJet   Tree : " <<  tch_pfjet->GetEntries() <<endl;
    
    //tr_jet = (TTree*)fin->Get("akPu3CaloJetAnalyzer/t");
    //TChain *tch_calojet = new TChain(Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
    //AddInputFiles(tch_calojet,kFileList,Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
    TChain *tch_calojet=0;
    if( kSpecies == "pbpb" ){
      tch_calojet = new TChain("akPu3CaloJetAnalyzer/t");
      AddInputFiles(tch_calojet,kFileList,"akPu3CaloJetAnalyzer/t");
      cout <<" # of events in CaloJet Tree : " <<  tch_calojet->GetEntries() <<endl;
    }else{
      tch_calojet = new TChain("ak3CaloJetAnalyzer/t");
      AddInputFiles(tch_calojet,kFileList,"ak3CaloJetAnalyzer/t");
      cout <<" # of events in CaloJet Tree : " <<  tch_calojet->GetEntries() <<endl;
    }
    //tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
    TChain *tch_ev = new TChain("hiEvtAnalyzer/HiTree");
    AddInputFiles(tch_ev,kFileList,"hiEvtAnalyzer/HiTree");
    cout <<" # of events in Event   Tree : " <<  tch_ev->GetEntries() <<endl;
    
    //tr_hlt = (TTree*)fin->Get("hltanalysis/HltTree");
    TChain *tch_hlt = new TChain("hltanalysis/HltTree");
    AddInputFiles(tch_hlt,kFileList,"hltanalysis/HltTree");
    cout <<" # of events in HLT     Tree : " <<  tch_hlt->GetEntries() <<endl;
    
    //tr_skim = (TTree*)fin->Get("skimanalysis/HltTree");
    TChain *tch_skim = new TChain("skimanalysis/HltTree");
    AddInputFiles(tch_skim,kFileList,"skimanalysis/HltTree");
    cout <<" # of events in Skim    Tree : " <<  tch_skim->GetEntries() <<endl;
    
    //tr_trobj = (TTree*)fin->Get("hltobject/jetObjTree");  
    // TChain *tch_trgobj = new TChain("hltobject/jetObjTree");
    // AddInputFiles(tch_trgobj,kFileList,"hltobject/jetObjTree");
    // cout <<" # of events in TrigObj Tree : " <<  tch_trgobj->GetEntries() <<endl;
    // cout <<endl;
    
    
    //! Event Tree
    int run_value;
    int evt_value;
    int lumi_value;
    int hiNpix;
    int hiBin;
    float vz;
    
    tch_ev->SetBranchAddress("run",&run_value);  
    tch_ev->SetBranchAddress("evt",&evt_value);  
    tch_ev->SetBranchAddress("lumi",&lumi_value);  
    tch_ev->SetBranchAddress("hiBin",&hiBin);
    tch_ev->SetBranchAddress("hiNpix",&hiNpix);  
    tch_ev->SetBranchAddress("vz",&vz);
    
    
    //! HLT tree
    //! HI
    int jet55;
    int jet55_prescl;
    int jet65;
    int jet65_prescl;
    int jet80;
    int jet80_prescl;
    //! PP
    int jet40;
    int jet40_prescl;
    int jet60;
    int jet60_prescl;
    
    if ( kSpecies == "pbpb" ){
      if( kDataset == "data" ){
	tch_hlt->SetBranchAddress("HLT_HIJet55_v1",&jet55);
	tch_hlt->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_prescl);
	tch_hlt->SetBranchAddress("HLT_HIJet65_v1",&jet65);
	tch_hlt->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_prescl);
	tch_hlt->SetBranchAddress("HLT_HIJet80_v1",&jet80);
	tch_hlt->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_prescl);
      }else{
	tch_hlt->SetBranchAddress("HLT_HIJet55_v7",&jet55);
	tch_hlt->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_prescl);
	tch_hlt->SetBranchAddress("HLT_HIJet65_v7",&jet65);
	tch_hlt->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_prescl);
	tch_hlt->SetBranchAddress("HLT_HIJet80_v7",&jet80);
	tch_hlt->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_prescl);
      }
    }else {
      tch_hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40);
      tch_hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_prescl);
      tch_hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60);
      tch_hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_prescl);
      tch_hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80);
      tch_hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_prescl);
    }
    // int L1_sj36_1;
    // int L1_sj36_p_1;
    // int L1_sj52_1;
    // int L1_sj52_p_1;
    // tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_1);
    // tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_1);
    // tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_1);
    // tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_1);
    
    
    //! Skim Tree
    int pcollisionEventSelection;
    int pHBHENoiseFilter;
    if( kSpecies == "pbpb" )tch_skim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    else tch_skim->SetBranchAddress("pPAcollisionEventSelectionPA",&pcollisionEventSelection);
    tch_skim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  
    //! Trigger object tree
    // float trgObj_id;
    // float trgObj_pt;
    // float trgObj_eta;
    // float trgObj_phi;
    // tch_trgobj->SetBranchAddress("id",&trgObj_id);
    // tch_trgobj->SetBranchAddress("pt",&trgObj_pt);
    // tch_trgobj->SetBranchAddress("eta",&trgObj_eta);
    // tch_trgobj->SetBranchAddress("phi",&trgObj_phi);

    
    //! CaloJet tree
    //Declaration of leaves types
    int   nref_calo;
    float jtpt_calo[1000];
    float rawpt_calo[1000];
    float jtpu_calo[1000];
    float jteta_calo[1000];
    float jtphi_calo[1000];
    float hcalSum_calo[1000];
    float ecalSum_calo[1000];
  
    int sid_calo[1000];
    float refpt_calo[1000];
    float refeta_calo[1000];
    float refphi_calo[1000];
    float refdrjt_calo[1000];

    tch_calojet->SetBranchAddress("nref",&nref_calo);
    tch_calojet->SetBranchAddress("rawpt",rawpt_calo);
    tch_calojet->SetBranchAddress("jtpt" ,jtpt_calo);
    tch_calojet->SetBranchAddress("jtpu" ,jtpu_calo);
    tch_calojet->SetBranchAddress("jteta",jteta_calo);
    tch_calojet->SetBranchAddress("jtphi",jtphi_calo);
    tch_calojet->SetBranchAddress("hcalSum",hcalSum_calo);
    tch_calojet->SetBranchAddress("ecalSum",ecalSum_calo);
    if( kDataset == "mc" ){
      tch_calojet->SetBranchAddress("subid" ,sid_calo);
      tch_calojet->SetBranchAddress("refpt" ,refpt_calo);
      tch_calojet->SetBranchAddress("refeta",refeta_calo);
      tch_calojet->SetBranchAddress("refphi",refphi_calo);
      tch_calojet->SetBranchAddress("refdrjt",refdrjt_calo);
    }


    //! PFJet tree
    //Declaration of leaves types
    int   nref;
    float jtpt[1000];
    float rawpt[1000];
    float jteta[1000];
    float jtpu [1000];
    float jtphi[1000];
    float neutralSum[1000];
    float chargedSum[1000];
    float photonSum[1000];
    float eleSum[1000];
    float muonSum[1000];
    float neutralMax[1000];
    float chargedMax[1000];
    float photonMax[1000];
    float eleMax[1000];
    float muonMax[1000];
    float hcalSum_pf[1000];
    float ecalSum_pf[1000];

    float pthat;
    int   sid[1000];
    float pfrefpt[1000];
    float pfrefeta[1000];
    float pfrefphi[1000];
    float pfrefdrjt[1000];

    tch_pfjet->SetBranchAddress("nref",&nref);
    tch_pfjet->SetBranchAddress("rawpt",rawpt);
    tch_pfjet->SetBranchAddress("jtpt" ,jtpt);
    tch_pfjet->SetBranchAddress("jtpu" ,jtpu);
    tch_pfjet->SetBranchAddress("jteta",jteta);
    tch_pfjet->SetBranchAddress("jtphi",jtphi);
    tch_pfjet->SetBranchAddress("neutralSum",neutralSum);
    tch_pfjet->SetBranchAddress("chargedSum",chargedSum);
    tch_pfjet->SetBranchAddress("photonSum",photonSum);
    tch_pfjet->SetBranchAddress("eSum",eleSum);
    tch_pfjet->SetBranchAddress("muSum",muonSum);
    tch_pfjet->SetBranchAddress("neutralMax",neutralMax);
    tch_pfjet->SetBranchAddress("chargedMax",chargedMax);
    tch_pfjet->SetBranchAddress("photonMax",photonMax);
    tch_pfjet->SetBranchAddress("eMax",eleMax);
    tch_pfjet->SetBranchAddress("muMax",muonMax);
    tch_pfjet->SetBranchAddress("hcalSum",hcalSum_pf);
    tch_pfjet->SetBranchAddress("ecalSum",ecalSum_pf);

    if( kDataset == "mc" ){
      tch_pfjet->SetBranchAddress("pthat",&pthat);    
      tch_pfjet->SetBranchAddress("subid" ,sid);
      tch_pfjet->SetBranchAddress("refpt" ,pfrefpt);
      tch_pfjet->SetBranchAddress("refeta",pfrefeta);
      tch_pfjet->SetBranchAddress("refphi",pfrefphi);
      tch_pfjet->SetBranchAddress("refdrjt",pfrefdrjt);
    }
  
    tch_pfjet->AddFriend(tch_ev);
    tch_pfjet->AddFriend(tch_hlt);
    tch_pfjet->AddFriend(tch_skim);
    //tch_pfjet->AddFriend(tch_trgobj);
    tch_pfjet->AddFriend(tch_calojet);

    //! Disable branches 
    //! Jet Tree
    tch_pfjet->SetBranchStatus("*",0,0);
    tch_pfjet->SetBranchStatus("nref" ,1,0);
    tch_pfjet->SetBranchStatus("rawpt",1,0);
    tch_pfjet->SetBranchStatus("jtpt" ,1,0);
    tch_pfjet->SetBranchStatus("jtpu" ,1,0);
    tch_pfjet->SetBranchStatus("jteta",1,0);
    tch_pfjet->SetBranchStatus("jtphi",1,0);
    tch_pfjet->SetBranchStatus("neutralSum",1,0);
    tch_pfjet->SetBranchStatus("chargedSum",1,0);
    tch_pfjet->SetBranchStatus("photonSum",1,0);
    tch_pfjet->SetBranchStatus("eSum",1,0);
    tch_pfjet->SetBranchStatus("muSum",1,0);
    tch_pfjet->SetBranchStatus("neutralMax",1,0);
    tch_pfjet->SetBranchStatus("chargedMax",1,0);
    tch_pfjet->SetBranchStatus("photonMax",1,0);
    tch_pfjet->SetBranchStatus("eMax",1,0);
    tch_pfjet->SetBranchStatus("muMax",1,0);
    tch_pfjet->SetBranchStatus("hcalSum",1,0);
    tch_pfjet->SetBranchStatus("ecalSum",1,0);

    tch_pfjet->SetBranchStatus("run",1,0);
    tch_pfjet->SetBranchStatus("evt",1,0);
    tch_pfjet->SetBranchStatus("lumi",1,0);
    tch_pfjet->SetBranchStatus("hiNpix",1,0);
    tch_pfjet->SetBranchStatus("hiBin",1,0);
    tch_pfjet->SetBranchStatus("vz",1,0);

    // tch_pfjet->SetBranchStatus("id",1,0);
    // tch_pfjet->SetBranchStatus("pt",1,0);
    // tch_pfjet->SetBranchStatus("eta",1,0);
    // tch_pfjet->SetBranchStatus("phi",1,0);

    if( kSpecies == "pbpb" ){
      if ( kDataset == "data" ){
	tch_pfjet->SetBranchStatus("HLT_HIJet55_v1",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet55_v1_Prescl",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet65_v1",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet65_v1_Prescl",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet80_v1",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet80_v1_Prescl",1,0);
	// tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND",1,0);
	// tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND_Prescl",1,0);
	// tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND",1,0);
	// tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND_Prescl",1,0);
      }else{
	tch_pfjet->SetBranchStatus("HLT_HIJet55_v7",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet55_v7_Prescl",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet65_v7",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet65_v7_Prescl",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet80_v7",1,0);
	tch_pfjet->SetBranchStatus("HLT_HIJet80_v7_Prescl",1,0);
      }
      tch_pfjet->SetBranchStatus("pcollisionEventSelection",1,0);
    }else{
      tch_pfjet->SetBranchStatus("HLT_PAJet40_NoJetID_v1",1,0);
      tch_pfjet->SetBranchStatus("HLT_PAJet40_NoJetID_v1_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_PAJet60_NoJetID_v1",1,0);
      tch_pfjet->SetBranchStatus("HLT_PAJet60_NoJetID_v1_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
      tch_pfjet->SetBranchStatus("HLT_PAJet80_NoJetID_v1_Prescl",1,0);
      tch_pfjet->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
    }
    tch_pfjet->SetBranchStatus("pHBHENoiseFilter",1,0);
  
    if( kDataset == "mc"){
      tch_pfjet->SetBranchStatus("pthat",1,0);    
      tch_pfjet->SetBranchStatus("subid" ,1,0);
      tch_pfjet->SetBranchStatus("refpt" ,1,0);
      tch_pfjet->SetBranchStatus("refeta",1,0);
      tch_pfjet->SetBranchStatus("refphi",1,0);
      tch_pfjet->SetBranchStatus("refdrjt",1,0);
    }
  


    std::cout<<"\t"<<std::endl;
    std::cout<<"\t"<<std::endl;
    std::cout<<"**************************************************** "<<std::endl;
    std::cout<<Form("Dataset : %s, Species : %s, Jet Algorithm :  %s ",kDataset.c_str(), kSpecies.c_str(), kAlgName.c_str())<<std::endl;
    std::cout<<Form("Outfile : %s",outfile.c_str())<<std::endl;
    std::cout<<Form("vertex z (c.m.) cut : %0.3f ",kvzcut)<<std::endl;
    std::cout<<Form("Reco pT cut : %0.3f ; Reco eta cut : %0.3f ",kptrecocut, ketacut)<<std::endl;
    std::cout<<Form("CALO-PF jet matching delta R cut   : %0.3f ",kdelrmatch)<<std::endl;
    std::cout<<"**************************************************** "<<std::endl;
    std::cout<<"\t"<<std::endl;
    std::cout<<"\t"<<std::endl;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    Long64_t nbytes=0;
    Long64_t nentries = tch_pfjet->GetEntries();
    std::cout<<Form("# of entries in TTree for %s %s : ",kAlgName.c_str(),kSpecies.c_str())<<nentries<<std::endl;

    double wxs=1.;
    if( kDataset == "mc" ){
      TEventList* el = new TEventList("el","el");
      stringstream selection; selection<<"pthat<="<<kMaxpthat;
      tch_pfjet->Draw(">>el",selection.str().c_str());
      double fentries = (double)el->GetN();
      std::cout<<"tree entries  :  "<<kAlgName.c_str()<<" algorithm : " << nentries<<" elist: "<< fentries <<std::endl;
      delete el;
      double tmpw = GetXsec(kMaxpthat);
      wxs = tmpw/(fentries/100000.);
    }

    //! Start event loop
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += tch_pfjet->GetEntry(i);
      
      if(printDebug && i==5000)break;
      //if(i==100)break;
      
      int iFill=1;
      if( i%2==1 )iFill=0;

      int trigFill=0;
      if( kSpecies == "pp" ){
	if( jet40==1 && jet60==0 && jet80==0 )trigFill=1;
	if( jet60==1 && jet80==0 )trigFill=1;      
	if( jet80==1 )trigFill=1;      
      }else{
	if( jet55==1 && jet65==0 && jet80==0 )trigFill=1;
	if( jet65==1 && jet80==0 )trigFill=1;      
	if( jet80==1 )trigFill=1;      
      }

      float rndm=gRandom->Rndm();

      hEvents_Total->Fill(rndm);
      if( kSpecies == "pp" ){
	if(jet40)hEvents_jet40->Fill(rndm);
	if(jet40_prescl)hEvents_jet40_prescl->Fill(jet40_prescl,rndm);
	if(jet60)hEvents_jet60->Fill(rndm);
	if(jet60_prescl)hEvents_jet60_prescl->Fill(jet60_prescl,rndm);
	if(jet80)hEvents_jet80->Fill(rndm);
	if(jet80_prescl)hEvents_jet80_prescl->Fill(jet80_prescl,rndm);

	if(jet40==1 && jet60==0 && jet80==0)hEvents_jet40_nojet60_nojet80->Fill(rndm);
	if(jet60==1 && jet80==0)hEvents_jet60_nojet80->Fill(rndm);
      }else{
	if(jet55)hEvents_jet55->Fill(rndm);
	if(jet55_prescl)hEvents_jet55_prescl->Fill(jet55_prescl,rndm);
	if(jet65)hEvents_jet65->Fill(rndm);
	if(jet65_prescl)hEvents_jet65_prescl->Fill(jet65_prescl,rndm);
	if(jet80)hEvents_jet80->Fill(rndm);
	if(jet80_prescl)hEvents_jet80_prescl->Fill(jet80_prescl,rndm);

	if(jet55==1 && jet65==0 && jet80==0)hEvents_jet55_nojet65_nojet80->Fill(rndm);
	if(jet65==1 && jet80==0)hEvents_jet65_nojet80->Fill(rndm);
      }

      if( kDataset == "data" ){
	if(pcollisionEventSelection)hEvents_pCollEvent->Fill(rndm);
	if(pcollisionEventSelection && pHBHENoiseFilter)hEvents_pHBHENoise->Fill(rndm);
	if(pcollisionEventSelection && pHBHENoiseFilter && fabs(vz)<kvzcut )hEvents_Vzcut->Fill(rndm);
	if(pcollisionEventSelection==0 || pHBHENoiseFilter==0 || fabs(vz) > kvzcut){
	  hEvents_bad->Fill(rndm);
	  continue;
	}
      }else if( kDataset == "mc" ){//! HBHENoiseFilter
	if(pcollisionEventSelection)hEvents_pCollEvent->Fill(rndm);
	if(pcollisionEventSelection)hEvents_pHBHENoise->Fill(rndm);
	if(pcollisionEventSelection && fabs(vz)<kvzcut )hEvents_Vzcut->Fill(rndm);
	if(pcollisionEventSelection==0 || fabs(vz) > kvzcut){
	  hEvents_bad->Fill(rndm);
	  continue;
	}
      }
    
      int iCent = -1;
      if( kSpecies == "pp" )iCent=ncen-1;
      else iCent = GetCentBin(hiBin);
      if(iCent<0 || iCent>=ncen)continue;

      if ( kSpecies == "pbpb" ){
	//! SuperNovae events use calo jets
	//int lowjetCounter=0;
	int jetCounter=0;
	for(int g = 0;g<nref_calo;g++){
	  if( fabs(jteta_calo[g]) < 2. && jtpt_calo[g]>=50. ){ //to select inside
	    jetCounter++;
	  }//eta selection cut
	}// jet loop
	// apply the correct supernova selection cut rejection here:
	if( hiNpix > 38000 - 500*jetCounter ){
	  hEvents_supernova->Fill(rndm);
	  continue;
	}
      }
      if( nref==0 && nref_calo==0 ){
	hEvents_nopfcalo->Fill(rndm);
	continue;
      }

      if( kDataset == "mc" && pthat > kMaxpthat ){
	hEvents_maxpthat->Fill(rndm);
	continue;
      }
      hEvents_Cent->Fill(rndm);


      if(printDebug)std::cout << "------------------------------------Start Event # " << i <<"------------------------------------------------------------------ " << std::endl;
      if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<< "  # of Calo jets " <<nref_calo<<std::endl;
      
      int pfjets=0;
      int calojets=0;
      
      std::vector <Jet> vPFJets, vCaloJets;
      std::vector <int> pfid(nref), caloid(nref_calo);

    
      if(printDebug)std::cout << " PF jets : " << std::endl;

      for(int pj=0; pj<nref; pj++){ //! PFjet
      
	if( rawpt[pj] < kptrawcut || jtpt[pj] < kptrecocut ) continue;
	if( fabs(jteta[pj]) > ketacut ) continue;
      
	Jet pfj;
	pfj.id  = pj;
	pfj.eta = jteta[pj];
	pfj.phi = jtphi[pj];
	pfj.pt  = jtpt [pj];

	if(printDebug){
	  if( kDataset == "mc" )std::cout <<"\t" << pj << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " phi : " << jtphi[pj] << " subid : " << sid[pj] << " dr : " << pfrefdrjt[pj] << std::endl;
	  else std::cout <<"\t"<< pj << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " phi : " << jtphi[pj] << std::endl;
	}
	vPFJets.push_back(pfj);
	pfjets++;
      }    

      if(printDebug){
	std::cout << std::endl;
	std::cout << " Calo jets : " << std::endl;
      }
      for(int cj=0; cj<nref_calo; cj++){ //! CaloJet
      
	if( rawpt_calo[cj] < kptrawcut || jtpt_calo[cj] < kptrecocut) continue;
	if( fabs(jteta_calo[cj]) > ketacut ) continue;

	Jet clj;
	clj.id  = cj;
	clj.eta = jteta_calo[cj];
	clj.phi = jtphi_calo[cj];
	clj.pt  = jtpt_calo[cj];

	if(printDebug){
	  if( kDataset == "mc" )std::cout <<"\t" << cj << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " phi : " << jtphi_calo[cj] << " subid : " << sid_calo[cj] << " dr : " << refdrjt_calo[cj] << std::endl;
	  else  std::cout <<"\t" << cj << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " phi : " << jtphi_calo[cj] << std::endl;
	}
	vCaloJets.push_back(clj);
	calojets++;
      }//! calo jet loop
      if(printDebug)std::cout << std::endl;
      if(pfjets==0 && calojets==0){
	if(printDebug){
	  std::cout <<" XXXXXXXXXXX  No Calo and PF jets passed the cuts " << std::endl;
	  std::cout << "------------------------------------End Event # " << i <<"------------------------------------------------------------------ " << "\n"<<std::endl;
	}
	hEvents_nopfcalo->Fill(rndm);
	continue;
      }

      bool onlyCalo   = (pfjets==0 && calojets >0) ? true : false;
      bool onlyPF     = (pfjets>0  && calojets==0) ? true : false;
      bool bothPFCalo = (pfjets>0  && calojets >0) ? true : false;

      int matchedJets=0;
      int unmatchedPFJets=0;
      int unmatchedCaloJets=0;

      double tmpwt = 1.;
      if( kDataset == "mc" ){
	double wvz  = fVz->Eval(vz);
	double wcen = hCentWeight->GetBinContent(hCentWeight->FindBin(hiBin));
	if ( kSpecies == "pbpb" )tmpwt = (wxs*wvz*wcen*lumi_scale);
	else tmpwt = wxs*lumi_scale;
      }else tmpwt = 1.;

      if(printDebug)std::cout <<" Total ==>  # of PF jets : " << pfjets << ", # of calojets : "  << calojets <<"\n"<<std::endl;

      std::vector < Jet >::const_iterator iJet, jJet;

      if( onlyPF ){
      
	for(iJet = vPFJets.begin(); iJet != vPFJets.end(); ++iJet){ //! PFjet

	  int pj = (*iJet).id; 

	  Float_t Sumcand = chargedSum[pj] + neutralSum[pj] + photonSum[pj] + muonSum[pj];
	  bool wJetId  = (eleMax[pj]/Sumcand < 0.05);
	  if( kSpecies == "pp" && kDataset=="mc" )wJetId=true;

	  if( sid[pj] != 0)continue;
	  if( pfrefpt[pj] > 2.5*pthat ) continue;
	
	  if( trigFill ) {
	    if( iFill ){
	      if( fabs(pfrefdrjt[pj]) < kdelrcut ){
		//! w/o jet id
		hgen   [0][iCent]->Fill(pfrefpt[pj],tmpwt);
		hrec   [0][iCent]->Fill(jtpt[pj],tmpwt);
		hmatrix[0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		
		hgen_f   [0][iCent]->Fill(pfrefpt[pj],tmpwt);
		hrec_f   [0][iCent]->Fill(jtpt[pj],tmpwt);
		hmatrix_f[0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		hresp    [0][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);
		hresp_f  [0][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);
		
		hgen_um   [0][iCent]->Fill(pfrefpt[pj],tmpwt);
		hrec_um   [0][iCent]->Fill(jtpt[pj],tmpwt);
		hmatrix_um[0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		
		hgen_um_f   [0][iCent]->Fill(pfrefpt[pj],tmpwt);
		hrec_um_f   [0][iCent]->Fill(jtpt[pj],tmpwt);
		hmatrix_um_f[0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		hresp_um    [0][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);
		hresp_um_f  [0][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);

		if( wJetId ){
		  hgen   [1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  hrec   [1][iCent]->Fill(jtpt[pj],tmpwt);
		  hmatrix[1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  
		  hgen_f   [1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  hrec_f   [1][iCent]->Fill(jtpt[pj],tmpwt);
		  hmatrix_f[1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  hresp    [1][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);
		  hresp_f  [1][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);
		  
		  hgen_um   [1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  hrec_um   [1][iCent]->Fill(jtpt[pj],tmpwt);
		  hmatrix_um[1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  
		  hgen_um_f   [1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  hrec_um_f   [1][iCent]->Fill(jtpt[pj],tmpwt);
		  hmatrix_um_f[1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  hresp_um    [1][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);
		  hresp_um_f  [1][iCent]->Fill(jtpt[pj], pfrefpt[pj], tmpwt);
		}
	      }else{
		hresp     [0][iCent]->Miss(pfrefpt[pj], tmpwt);
		hresp_um  [0][iCent]->Miss(pfrefpt[pj], tmpwt);
		hresp_f   [0][iCent]->Miss(pfrefpt[pj], tmpwt);
		hresp_um_f[0][iCent]->Miss(pfrefpt[pj], tmpwt);
		if( wJetId ){
		  hresp     [1][iCent]->Miss(pfrefpt[pj], tmpwt);
		  hresp_um  [1][iCent]->Miss(pfrefpt[pj], tmpwt);
		  hresp_f   [1][iCent]->Miss(pfrefpt[pj], tmpwt);
		  hresp_um_f[1][iCent]->Miss(pfrefpt[pj], tmpwt);
		}
	      }
	    }else{
	      if( fabs(pfrefdrjt[pj]) < kdelrcut ){
		hmatrix_c  [0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		hmatrix_c_f[0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		hmatrix_um_c[0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		hmatrix_um_c_f[0][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		
		hgen_c     [0][iCent]->Fill(pfrefpt[pj],tmpwt);
		hgen_c_f   [0][iCent]->Fill(pfrefpt[pj],tmpwt);
		hgen_um_c  [0][iCent]->Fill(pfrefpt[pj],tmpwt);
		hgen_um_c_f[0][iCent]->Fill(pfrefpt[pj],tmpwt);
		
		hrec_c     [0][iCent]->Fill(jtpt[pj],tmpwt);
		hrec_c_f   [0][iCent]->Fill(jtpt[pj],tmpwt);
		hrec_um_c  [0][iCent]->Fill(jtpt[pj],tmpwt);
		hrec_um_c_f[0][iCent]->Fill(jtpt[pj],tmpwt);
		if( wJetId ){
		  hmatrix_c  [1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  hmatrix_c_f[1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  hmatrix_um_c[1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  hmatrix_um_c_f[1][iCent]->Fill(jtpt[pj],pfrefpt[pj],tmpwt);
		  
		  hgen_c     [1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  hgen_c_f   [1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  hgen_um_c  [1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  hgen_um_c_f[1][iCent]->Fill(pfrefpt[pj],tmpwt);
		  
		  hrec_c     [1][iCent]->Fill(jtpt[pj],tmpwt);
		  hrec_c_f   [1][iCent]->Fill(jtpt[pj],tmpwt);
		  hrec_um_c  [1][iCent]->Fill(jtpt[pj],tmpwt);
		  hrec_um_c_f[1][iCent]->Fill(jtpt[pj],tmpwt);
		}
	      }
	    }
	  }	
	  if(printDebug)std::cout <<" unmatched PF jets w ncalo=0 : " << unmatchedPFJets << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " phi : " << jtphi[pj] << " dr : " << pfrefdrjt[pj] <<std::endl;
	  unmatchedPFJets++;
	}
      }
      else if( onlyCalo ){
	// for(iJet = vCaloJets.begin(); iJet != vCaloJets.end(); ++iJet){ //! Calojet
	
	//   int cj = (*iJet).id; 	
	
	//   if(printDebug)std::cout <<" unmatched CALO jets w npf=0 : " << unmatchedCaloJets << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " phi : " << jtphi_calo[cj] << " dr : " << refdrjt_calo[cj] << std::endl;
	//   unmatchedCaloJets++;
	// }
      }
      else if( bothPFCalo ){

	CaloPFMatchedJets mCaloPFMatchedJets;
	for(iJet = vCaloJets.begin(); iJet != vCaloJets.end(); ++iJet){ //! Calojet      
	
	  for(jJet = vPFJets.begin(); jJet != vPFJets.end(); ++jJet){ //! PFjet
	  
	    mCaloPFMatchedJets.insert(std::make_pair(*iJet,*jJet));
	  
	  }//! calo jet loop
	}//! PF jet loop
      
	CPFItr itr;
	//! Matched jets (PF jet matched to Calo jet)
	for(itr = mCaloPFMatchedJets.begin(); itr != mCaloPFMatchedJets.end(); ++itr){

	  CaloPFJetPair jetpair = (*itr);
	  Jet clj = jetpair.first;
	  Jet pfj = jetpair.second;

	  float delr  = deltaR(clj.eta, clj.phi, pfj.eta, pfj.phi);
	  //float delpt = fabs(clj.pt - pfj.pt);
	  if( delr < kdelrmatch && caloid[clj.id]==0 ){
	  
	    Float_t Sumcand = chargedSum[pfj.id] + neutralSum[pfj.id] + photonSum[pfj.id] + muonSum[pfj.id];
	    Float_t ePFSel  = (18./7.*jtpt_calo[clj.id]/jtpt[pfj.id]) - 9./7.;
	    bool wJetId     = (jtpt_calo[clj.id]/jtpt[pfj.id]  > 0.50) && (jtpt_calo[clj.id]/jtpt[pfj.id] <= 0.85) && (eleMax[pfj.id]/Sumcand < ePFSel);

	    if( kSpecies == "pp" && kDataset=="mc" )wJetId=true;
	  
	    if( trigFill ) {
	      if( iFill ){
		if( pfrefpt[pfj.id] < 2.5*pthat ) {
		  if( sid[pfj.id] == 0 ){
		    if ( fabs(pfrefdrjt[pfj.id]) < kdelrcut ){
		      
		      //! w/o jet id
		      hgen[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
		      hgen_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
		      hgen_m[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec_m[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix_m[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		    
		      hgen_m_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec_m_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix_m_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);

		      hresp  [0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
		      hresp_m[0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    

		      hresp_f  [0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
		      hresp_m_f[0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    
		    
		      if( wJetId ){
			hgen[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
			hgen_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
			hgen_m[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_m[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix_m[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
			hgen_m_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_m_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix_m_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);

			hresp  [1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
			hresp_m[1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    
			hresp_f  [1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
			hresp_m_f[1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    
		      }
		    }else{
		      hresp[0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
		      hresp_m[0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
		      hresp_f[0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
		      hresp_m_f[0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
		      if( wJetId ){
			hresp[1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
			hresp_m[1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
			hresp_f[1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
			hresp_m_f[1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
		      }
		    }
		  }
		}
	      }else{
		if( pfrefpt[pfj.id] < 2.5*pthat ) {
		  if( sid[pfj.id] == 0 ){
		    if( fabs(pfrefdrjt[pfj.id]) < kdelrcut ){ 
		      hmatrix_c[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      hmatrix_c_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      hmatrix_m_c[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      hmatrix_m_c_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);

		      hgen_c[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hgen_c_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hgen_m_c[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hgen_m_c_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		    
		      hrec_c[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hrec_c_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hrec_m_c[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hrec_m_c_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      if( wJetId ){
			hmatrix_c[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hmatrix_c_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hmatrix_m_c[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hmatrix_m_c_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hgen_c[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hgen_c_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hgen_m_c[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hgen_m_c_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_c[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hrec_c_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hrec_m_c[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hrec_m_c_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      }
		    }
		  }
		}
	      }
	    }
	    if(printDebug)std::cout <<"\t *** Matched" << " id : " << pfj.id << " " << clj.id << " pfpt :  " << pfj.pt << " calopt : " << clj.pt << " spf : " << sid[pfj.id] << " scj : " << sid [clj.id] << " dr : " << pfrefdrjt[pfj.id] << std::endl;
	    //if( sid[pfj.id]!=0 && pfpt>20 )std::cout <<"\t *** Matched" << " pthat : " << pthat << " refpt : " << refpt << " pfpt :  " << pfj.pt << " calopt : " << clj.pt << " spf : " << sid[pfj.id] << " scj : " << sid [clj.id] << std::endl;	  
	  
	    matchedJets++;	  
	    pfid  [pfj.id] = 1;
	    caloid[clj.id] = 1;
	  }
	}//! iJet
      
	// //! Unmatched jets
	for(itr = mCaloPFMatchedJets.begin(); itr != mCaloPFMatchedJets.end(); ++itr){
	  CaloPFJetPair jetpair = (*itr);

	  Jet clj = jetpair.first;
	  Jet pfj = jetpair.second;

	  //float delr  = deltaR(pfj.eta, pfj.phi, clj.eta, clj.phi);
	  //float delpt = fabs(pfj.pt - clj.pt);
	  //if( pfid[pfj.id]==1 || caloid[clj.id]==1 )continue;
      
	  if( pfid[pfj.id] == 0 ){//! Unmatched PF jet

	    Float_t Sumcand = chargedSum[pfj.id] + neutralSum[pfj.id] + photonSum[pfj.id] + muonSum[pfj.id];
	    bool wJetId  = (eleMax[pfj.id]/Sumcand < 0.05);
	  
	    if( trigFill ) {
	      if( iFill ){
		if( pfrefpt[pfj.id] < 2.5*pthat) {
		  if( sid[pfj.id] == 0 ){ 
		    if( fabs(pfrefdrjt[pfj.id]) < kdelrcut ){ 
		      
		      //! w/o jet id
		      hgen[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		    
		      hgen_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		    
		      hgen_um[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec_um[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix_um[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		    
		      hgen_um_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hrec_um_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hmatrix_um_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);

		      hresp  [0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
		      hresp_um[0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    
		      hresp_f  [0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
		      hresp_um_f[0][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    
		    
		      if( wJetId ){
			hgen[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
			hgen_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
			hgen_um[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_um[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix_um[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      
			hgen_um_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_um_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hmatrix_um_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hresp  [1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
			hresp_um[1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    
			hresp_f  [1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);
			hresp_um_f[1][iCent]->Fill(jtpt[pfj.id], pfrefpt[pfj.id], tmpwt);		    
		      }
		    }else{
		      hresp   [0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
		      hresp_um[0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);		    		    
		      hresp_f   [0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
		      hresp_um_f[0][iCent]->Miss(pfrefpt[pfj.id], tmpwt);		    		    
		      if( wJetId ){
			hresp   [1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
			hresp_um[1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);		    
			hresp_f   [1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);
			hresp_um_f[1][iCent]->Miss(pfrefpt[pfj.id], tmpwt);		    
		      }
		    }
		  }
		}
	      }else{
		if( pfrefpt[pfj.id] < 2.5*pthat) {
		  if( sid[pfj.id] == 0 ){ 
		    if( fabs(pfrefdrjt[pfj.id]) < kdelrcut ){ 
		      hmatrix_c[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      hmatrix_c_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      hmatrix_um_c[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
		      hmatrix_um_c_f[0][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);

		      hgen_c[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hgen_c_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hgen_um_c[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		      hgen_um_c_f[0][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
		    
		      hrec_c[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hrec_c_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hrec_um_c[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      hrec_um_c_f[0][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      if( wJetId ){
			hmatrix_c[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hmatrix_c_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hmatrix_um_c[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			hmatrix_um_c_f[1][iCent]->Fill(jtpt[pfj.id],pfrefpt[pfj.id],tmpwt);
			
			hgen_c[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hgen_c_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hgen_um_c[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hgen_um_c_f[1][iCent]->Fill(pfrefpt[pfj.id],tmpwt);
			hrec_c[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hrec_c_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hrec_um_c[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
			hrec_um_c_f[1][iCent]->Fill(jtpt[pfj.id],tmpwt);
		      }
		    }
		  }
		}
	      }
	    }
	    if(printDebug)std::cout <<"\t %%% UnMatched PF" << " id  : " << pfj.id << " pfpt :  " << pfj.pt << " subid : " << sid[pfj.id] << " dr : " << pfrefdrjt[pfj.id] << std::endl;
	  
	    unmatchedPFJets++;	  	  
	    pfid  [pfj.id] = 1;	  
	  }
      
	  if( caloid[clj.id] == 0 ){//! Unmatched Calo jet

	    if(printDebug)std::cout <<"\t XXX UnMatched Calo" << " id : " << clj.id  << " calopt : " << clj.pt << " subid : " << sid_calo [clj.id] << " dr : " << refdrjt_calo[clj.id] << std::endl;
	
	    unmatchedCaloJets++;	  	  
	    caloid[clj.id] = 1;	  
	  }
	}
      }//! bothPFCalo
    
      if(printDebug){
	std::cout << std::endl;
	if( bothPFCalo    )std::cout<<" ****** Both PFCalo Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
	else if( onlyCalo )std::cout<<" ****** Only Calo   Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
	else if( onlyPF   )std::cout<<" ****** Only PF     Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
	std::cout << " " 
		  <<" All ==> PF : " << nref <<", CALO : "<< nref_calo << ";"  
		  <<" After cut ==> PF : " << pfjets << ", Calo : "  << calojets << ";"
		  <<" mCALOPF : "<< matchedJets <<", unmPF : "<<  unmatchedPFJets <<", unmCalo : "<<  unmatchedCaloJets 
		  <<std::endl;
	std::cout << "------------------------------------End Event # " << i <<"------------------------------------------------------------------ " << "\n"<<std::endl;
      }
    }//! event loop ends
  }//! kFiles


  for(int i=0; i<2; i++){
    for(int ic=mincen; ic<maxcen; ic++){
      unf_bayes = new RooUnfoldBayes(hresp[i][ic], hrec[i][ic], iterations);
      hunf[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf[i][ic]->SetName(Form("hunf_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hcov[i][ic] = CorrelationHist(unf_bayes->Ereco((RooUnfold::ErrorTreatment)2), Form("hcov_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic),
				    "Unfolded covariance matrix",hrec[i][ic]->GetNbinsX(), hrec[i][ic]->GetXaxis()->GetXbins()->GetArray());
      unf_bayes = new RooUnfoldBayes(hresp_m[i][ic], hrec_m[i][ic], iterations);
      hunf_m[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_m[i][ic]->SetName(Form("hunf_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hcov_m[i][ic] = CorrelationHist(unf_bayes->Ereco((RooUnfold::ErrorTreatment)2), Form("hcov_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic),
				      "matched Unfolded covariance matrix",hrec_m[i][ic]->GetNbinsX(), hrec_m[i][ic]->GetXaxis()->GetXbins()->GetArray());
      unf_bayes = new RooUnfoldBayes(hresp_um[i][ic], hrec_um[i][ic], iterations);
      hunf_um[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_um[i][ic]->SetName(Form("hunf_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hcov_um[i][ic] = CorrelationHist(unf_bayes->Ereco((RooUnfold::ErrorTreatment)2), Form("hcov_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic),
				      "un matched Unfolded covariance matrix",hrec_um[i][ic]->GetNbinsX(), hrec_um[i][ic]->GetXaxis()->GetXbins()->GetArray());
      

      //! Now check the unfolding on the rest of the data
      unf_bayes = new RooUnfoldBayes(hresp[i][ic], hrec_c[i][ic], iterations);
      hunf_c[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_c[i][ic]->SetName(Form("hunf_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      unf_bayes = new RooUnfoldBayes(hresp_m[i][ic], hrec_m_c[i][ic], iterations);
      hunf_m_c[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_m_c[i][ic]->SetName(Form("hunf_m_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      unf_bayes = new RooUnfoldBayes(hresp_um[i][ic], hrec_um_c[i][ic], iterations);
      hunf_um_c[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_um_c[i][ic]->SetName(Form("hunf_um_c_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));


      hresponse[i][ic] = (TH2F*)hresp[i][ic]->HresponseNoOverflow();
      hresponse[i][ic]->SetName(Form("hresponse_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hresponse_m[i][ic] = (TH2F*)hresp_m[i][ic]->HresponseNoOverflow();
      hresponse_m[i][ic]->SetName(Form("hresponse_m_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hresponse_um[i][ic] = (TH2F*)hresp_um[i][ic]->HresponseNoOverflow();
      hresponse_um[i][ic]->SetName(Form("hresponse_um_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));

      //! finer bins in pT
      unf_bayes = new RooUnfoldBayes(hresp_f[i][ic], hrec_f[i][ic], iterations);
      hunf_f[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_f[i][ic]->SetName(Form("hunf_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hcov_f[i][ic] = CorrelationHist(unf_bayes->Ereco((RooUnfold::ErrorTreatment)2), Form("hcov_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic),
				    "Unfolded covariance matrix",hrec_f[i][ic]->GetNbinsX(), hrec_f[i][ic]->GetXaxis()->GetXbins()->GetArray());
      unf_bayes = new RooUnfoldBayes(hresp_m_f[i][ic], hrec_m_f[i][ic], iterations);
      hunf_m_f[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_m_f[i][ic]->SetName(Form("hunf_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hcov_m_f[i][ic] = CorrelationHist(unf_bayes->Ereco((RooUnfold::ErrorTreatment)2), Form("hcov_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic),
				    "matched Unfolded covariance matrix",hrec_m_f[i][ic]->GetNbinsX(), hrec_m_f[i][ic]->GetXaxis()->GetXbins()->GetArray());
      unf_bayes = new RooUnfoldBayes(hresp_um_f[i][ic], hrec_um_f[i][ic], iterations);
      hunf_um_f[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_um_f[i][ic]->SetName(Form("hunf_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hcov_um_f[i][ic] = CorrelationHist(unf_bayes->Ereco((RooUnfold::ErrorTreatment)2), Form("hcov_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic),
				    "un matched Unfolded covariance matrix",hrec_um_f[i][ic]->GetNbinsX(), hrec_um_f[i][ic]->GetXaxis()->GetXbins()->GetArray());

      //! Now check the unfolding on the rest of the data
      unf_bayes = new RooUnfoldBayes(hresp_f[i][ic], hrec_c_f[i][ic], iterations);
      hunf_c_f[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_c_f[i][ic]->SetName(Form("hunf_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      unf_bayes = new RooUnfoldBayes(hresp_m_f[i][ic], hrec_m_c_f[i][ic], iterations);
      hunf_m_c_f[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_m_c_f[i][ic]->SetName(Form("hunf_m_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      unf_bayes = new RooUnfoldBayes(hresp_um_f[i][ic], hrec_um_c_f[i][ic], iterations);
      hunf_um_c_f[i][ic] = (TH1F*)unf_bayes->Hreco();
      hunf_um_c_f[i][ic]->SetName(Form("hunf_um_c_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));


      hresponse_f[i][ic] = (TH2F*)hresp_f[i][ic]->HresponseNoOverflow();
      hresponse_f[i][ic]->SetName(Form("hresponse_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hresponse_m_f[i][ic] = (TH2F*)hresp_m_f[i][ic]->HresponseNoOverflow();
      hresponse_m_f[i][ic]->SetName(Form("hresponse_m_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));
      hresponse_um_f[i][ic] = (TH2F*)hresp_um_f[i][ic]->HresponseNoOverflow();
      hresponse_um_f[i][ic]->SetName(Form("hresponse_um_f_%s_%s_%d_%d",kSpecies.c_str(), kAlgName.c_str(),i,ic));

    }
  }

  //! Write to output file
  fout->mkdir(Form("%sJetAnalyzer",kAlgName.c_str()));
  fout->cd(Form("%sJetAnalyzer", kAlgName.c_str()));
  for(int i=0; i<2; i++){
    for(int ic=mincen; ic<maxcen; ic++){
      hmatrix[i][ic]->Write(); 
      hmatrix_m[i][ic]->Write();
      hmatrix_um[i][ic]->Write();
      hgen[i][ic]->Write();
      hgen_m[i][ic]->Write();
      hgen_um[i][ic]->Write();
      hrec[i][ic]->Write();
      hrec_m[i][ic]->Write();
      hrec_um[i][ic]->Write();

      hmatrix_c[i][ic]->Write(); 
      hmatrix_m_c[i][ic]->Write();
      hmatrix_um_c[i][ic]->Write();
      hgen_c[i][ic]->Write();
      hgen_m_c[i][ic]->Write();
      hgen_um_c[i][ic]->Write();
      hrec_c[i][ic]->Write();
      hrec_m_c[i][ic]->Write();
      hrec_um_c[i][ic]->Write();

      hcov[i][ic]->Write();
      hcov_m[i][ic]->Write();
      hcov_um[i][ic]->Write();

      hcov_f[i][ic]->Write();
      hcov_m_f[i][ic]->Write();
      hcov_um_f[i][ic]->Write();

      hunf[i][ic]->Write();
      hunf_m[i][ic]->Write();
      hunf_um[i][ic]->Write();

      hunf_c[i][ic]->Write();
      hunf_m_c[i][ic]->Write();
      hunf_um_c[i][ic]->Write();

      hresponse[i][ic]->Write();
      hresponse_m[i][ic]->Write();
      hresponse_um[i][ic]->Write();

      hresp[i][ic]->Write();
      hresp_m[i][ic]->Write();
      hresp_um[i][ic]->Write();

      hmatrix_f[i][ic]->Write(); 
      hmatrix_m_f[i][ic]->Write();
      hmatrix_um_f[i][ic]->Write();
      hgen_f[i][ic]->Write();
      hgen_m_f[i][ic]->Write();
      hgen_um_f[i][ic]->Write();
      hrec_f[i][ic]->Write();
      hrec_m_f[i][ic]->Write();
      hrec_um_f[i][ic]->Write();

      hmatrix_c_f[i][ic]->Write(); 
      hmatrix_m_c_f[i][ic]->Write();
      hmatrix_um_c_f[i][ic]->Write();
      hgen_c_f[i][ic]->Write();
      hgen_m_c_f[i][ic]->Write();
      hgen_um_c_f[i][ic]->Write();
      hrec_c_f[i][ic]->Write();
      hrec_m_c_f[i][ic]->Write();
      hrec_um_c_f[i][ic]->Write();
      
      hunf_f[i][ic]->Write();
      hunf_m_f[i][ic]->Write();
      hunf_um_f[i][ic]->Write();

      hunf_c_f[i][ic]->Write();
      hunf_m_c_f[i][ic]->Write();
      hunf_um_c_f[i][ic]->Write();

      hresponse_f[i][ic]->Write();
      hresponse_m_f[i][ic]->Write();
      hresponse_um_f[i][ic]->Write();

      hresp_f[i][ic]->Write();
      hresp_m_f[i][ic]->Write();
      hresp_um_f[i][ic]->Write();
    }
  }
  //fout->Write();
  fout->cd("../");
  fout->Close();

  //! Check
  timer.Stop();
  double rtime  = timer.RealTime();
  double ctime  = timer.CpuTime();
  
  std::cout<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  return 1;
}

int GetEtaBin(float eta)
{
  for(int ix=0;ix<neta;ix++){
    if(eta>=etabins[ix] && eta<etabins[ix+1]){
      return ix;
    }
  }
  return -1;
}
int GetPhiBin(float phi)
{
  for(int ix=0;ix<nphi;ix++){
    if(phi>=phibins[ix] && phi<phibins[ix+1]){
      return ix;
    }
  }
  return -1;
}

int GetPtBin(float pt)
{
  for(int ix=0;ix<nbins;ix++){
    if(pt>=ptbins[ix] && pt<ptbins[ix+1]){
      return ix;
    }
  }
  return -1;
}
float delphi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  dphi = fabs(atan2(sin(dphi),cos(dphi)));
  return dphi;
}
int GetCentBin(int bin)
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
double GetXsec(double maxpt)
{
  //std::cout << " GetXsec() ::: max pt hat : " << maxpt << std::endl;                                                                          
  double effxsec=0;
  for(int i=0; i<11; i++){
    if(fabs(maxpt - xsec[i][2]) < 1e-08){
      effxsec = xsec[i][0] - xsec[i+1][0];
      //effxsec = xsec[i][0];                                                                                                                   
      //std::cout <<"\t \t  effective xsec : " << effxsec << "\t"<<  xsec[i][0] << "\t pthat : "<< xsec[i][1] << std::endl;                     
      return effxsec;
    }
  }
  return  1;
}
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


void AddInputFiles(TChain *tch, string inputname, string inputTree)
{
  //cout << " Add input files " << tch->GetName() << "  " << inputname.c_str() <<endl;
  //cout<<endl;
  std::stringstream tmpstr(inputname.c_str());
  std::string segment;
  std::vector<string> infiles;

  while(std::getline(tmpstr, segment, ',')){
    infiles.push_back(segment);
  }

  int im=0;
  for(std::vector<string>::iterator it=infiles.begin(); it!=infiles.end(); ++it, im++){
    std::string sFile = (*it);
    //   cout <<im<<"  "<< sFile.c_str() << endl;
    string stree = sFile+"/"+inputTree;
    tch->Add(stree.c_str());
  }
  //  cout<<endl;
}
TH2F* CorrelationHist (const TMatrixD& cov, const char* name, const char* title,
		       const unsigned nbins_rec, const Double_t* bins)
{
  Int_t nb= cov.GetNrows();
  TH2F* h= new TH2F (name, title, nbins_rec, bins, nbins_rec, bins);
  h->SetAxisRange (-1.0, 1.0, "Z");
  for(int i=0; i < nb; i++)
    for(int j=0; j < nb; j++) {
      Float_t Viijj= cov(i,i)*cov(j,j);
      if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
    }
  return h;
}
