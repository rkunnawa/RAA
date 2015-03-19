//
// Raghav Kunnawalkam Elayavalli
// Jan 15th 2014
// Rutgers, raghav.k.e at CERN dot CH
//
// Macro to do calo-pf jet correlation (spacial) to study the PF electron issues - this would tell us if PF jets from a jet55 trigger has a calo jet as a background. following is how to do that: 
//load the forest data files. and then load the calo jets first followed by the pf jet to search for spacial correlation. Once we find the pf jet which is in the same spacial location (which is found by calculating the delta R for all pf jets for one calo jet and the jet which is less than the value of delta R given above becomes it correlated jet), the ratio of the pT of the calo jet and the pf jet (correlated) is calculated. and that distribution is expected to be a gaussian with the outliers clearly saying there is something wrong with the candidates like electron etc...  
//


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

#define pi 3.14159265

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

static const int nbins_cent = 6;
static const double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
static const char centWidth[nbins_cent+1][256] = {"0 - 5 %","5 - 10 %","10 - 30 %","30 - 50 %","50 - 70 %","70 - 90 %","0 - 100 %"};
static const int no_radius = 1;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {3};
static const int TrigValue = 4;
static const char TrigName [TrigValue][256] = {"HLT55","HLT65","HLT80","FullDataset"};

int findBin(int hiBin){
  int binNo = -1;

  for(int i = 0;i<nbins_cent;i++){
    if(hiBin>=5*boundaries_cent[i] && hiBin<5*boundaries_cent[i+1]) {
      binNo = i;
      break;
    }
  }
  return binNo;
}

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

double Calc_deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  double dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}

using namespace std;

void RAA_calo_pf_JetCorrelation(int startfile = 0, int endfile = 1, int radius=3, char *algo = "Pu", int deltaR=2/*which i will divide by 10 later when using*/, Float_t CALOPTCUT = 30.0, Float_t PFPTCUT = 30.0, char *dataset = "Data"){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  timer.Start();

  TDatime date;

  cout<<"Running for Algo = "<<algo<<" "<<endl;
  
  bool printDebug = false;

  // Since this has to run on the HiForest files, it is best to run them as condor jobs similar to the PbPb data read macro, along with the jet trees which get their branch address set. 
  
  std::string infile1;
  if(dataset == "Data")infile1 = "jetRAA_PbPb_data_forest.txt";
  if(dataset == "MC")infile1 = "jetRAA_PbPb_mc_forest.txt";
  if(dataset == "MinBiasUPC")infile1 = "jetRAA_MinBiasUPC_forest.txt";
  
  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;

  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }
  
  const int N = 5;
  
  TChain *jetpbpb1[N][no_radius];

  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%dCaloJetAnalyzer",algo,list_radius[k]);
    dir[3][k] = Form("ak%s%dPFJetAnalyzer",algo,list_radius[k]);
    dir[4][k] = "hiEvtAnalyzer";
    //dir[5][k] = "hltobject";
    //dir[6][k] = "pfcandAnalyzer";
  }
  
  
  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "t",
    "HiTree"
    //"jetObjTree",
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
	
	jetpbpb1[t][k]->Add(filename1.c_str());
	if(printDebug)cout << "Tree loaded  " << string(dir[t][k]+"/"+trees[t]).data() << endl;
	if(printDebug)cout << "Entries : " << jetpbpb1[t][k]->GetEntries() << endl;
		
      }// tree loop ends
      
    }// radius loop ends
    
  }// file loop ends

  for(int k = 0;k<no_radius;k++){
    jetpbpb1[2][k]->AddFriend(jetpbpb1[0][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[1][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[4][k]);
    //jetpbpb1[2][k]->AddFriend(jetpbpb1[5][k]);
    
    jetpbpb1[3][k]->AddFriend(jetpbpb1[0][k]);
    jetpbpb1[3][k]->AddFriend(jetpbpb1[1][k]);
    jetpbpb1[3][k]->AddFriend(jetpbpb1[4][k]);
    //jetpbpb1[3][k]->AddFriend(jetpbpb1[5][k]);
    
  }// radius loop ends
 
  // Ok this should now work. 

  cout<<"total no of entries in the combined forest files = "<<jetpbpb1[2][0]->GetEntries()<<endl;
  //  if(printDebug)cout<<"total no of entries in the Jet80 Tree     = "<<jetpbpb2[2][0]->GetEntries()<<endl;

  //set the branch addresses:
  // jet tree 1 - Calo 
  int nrefe_1;
  float pt_1[1000];
  float raw_1[1000];
  float refpt_1[1000];
  float eta_1[1000];
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
  float hcalSum_1[1000];
  float ecalSum_1[1000];
  int subid_1[1000];
  
  // jet tree 2 - PF
  int nrefe_2;
  float pt_2[1000];
  float raw_2[1000];
  float refpt_2[1000];
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
  float hcalSum_2[1000];
  float ecalSum_2[1000];
  int subid_2[1000];
  
  // event tree
  int evt_1;
  int run_1;
  int lumi_1;
  int hiBin_1;
  float pthat_1;
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

  for(int k = 0;k<no_radius;k++){
    //set the branch addresses:  - one of the most boring parts of the code: 
    jetpbpb1[2][k]->SetBranchAddress("evt",&evt_1);
    jetpbpb1[2][k]->SetBranchAddress("run",&run_1);
    jetpbpb1[2][k]->SetBranchAddress("lumi",&lumi_1);
    jetpbpb1[2][k]->SetBranchAddress("hiBin",&hiBin_1);
    if(dataset=="MC") jetpbpb1[2][k]->SetBranchAddress("pthat",&pthat_1);
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
    jetpbpb1[2][k]->SetBranchAddress("ecalSum",&ecalSum_1);
    jetpbpb1[2][k]->SetBranchAddress("hcalSum",&hcalSum_1);
    
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
    jetpbpb1[3][k]->SetBranchAddress("ecalSum",&ecalSum_2);
    jetpbpb1[3][k]->SetBranchAddress("hcalSum",&hcalSum_2);
    
    if(dataset=="MC"){

      jetpbpb1[2][k]->SetBranchAddress("subid",&subid_1);
      jetpbpb1[3][k]->SetBranchAddress("subid",&subid_2);
      jetpbpb1[2][k]->SetBranchAddress("refpt",&refpt_1);
      jetpbpb1[3][k]->SetBranchAddress("refpt",&refpt_2);

    }
    
    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);

    // dont forget the include this for data and comment out the second set which is for MC. 
    if(dataset == "MC"){
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v7",&jet55_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_p_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v7",&jet65_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_p_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v7",&jet80_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_p_1);
    }
    if(dataset=="Data" || dataset=="MinBiasUPC"){
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v1",&jet55_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v1",&jet65_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v1",&jet80_1);
      jetpbpb1[2][k]->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_1);
    }
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_1);
    jetpbpb1[2][k]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_1);

    // jetpbpb1[2][k]->SetBranchAddress("id",&trgObj_id_1);
    // jetpbpb1[2][k]->SetBranchAddress("pt",&trgObj_pt_1);
    // jetpbpb1[2][k]->SetBranchAddress("eta",&trgObj_eta_1);
    // jetpbpb1[2][k]->SetBranchAddress("phi",&trgObj_phi_1);
    // jetpbpb1[2][k]->SetBranchAddress("mass",&trgObj_mass_1);
    
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

  }// radius loop
  cout<<"after branch declaration"<<endl;

  // // declare the histograms:
  // TH1F *hCaloPFCorr[TrigValue][nbins_cent+1], *hCalo[TrigValue][nbins_cent+1], *hPF[TrigValue][nbins_cent+1], *hRatio[TrigValue][nbins_cent+1];
  // TH2F * hCaloPFpt[TrigValue][nbins_cent+1];
  // TH2F * hCaloPFCorr_pt[TrigValue][nbins_cent+1];
  // TH2F * hDeltaR_deltapT[TrigValue][nbins_cent+1];
  // // TH1F * hDeltaR[TrigValue][nbins_cent+1];
  // // TH1F * hDeltapT[TrigValue][nbins_cent+1];

  // for(int a = 0;a<TrigValue;a++){
  //   for(int i = 0;i<nbins_cent+1;i++){
  //     hCaloPFCorr[a][i] = new TH1F(Form("hCaloPFCorr_%s_cent%d",TrigName[a],i),Form("Ratio Calo jet pT to correlated PF jet for %s %s",TrigName[a],centWidth[i]),200,0,10);
  //     hCaloPFCorr[a][i]->SetXTitle("PF Jet pT / Calo Jet PT");
  //     hCaloPFCorr[a][i]->SetYTitle("Counts");
      
  //     hCaloPFpt[a][i] = new TH2F(Form("hCaloPFpt_%s_cent%d",TrigName[a],i), Form("Matched Calo jet pT vs PF Jet pT %s %s", TrigName[a], centWidth[i]), nbins_pt,boundaries_pt, nbins_pt, boundaries_pt);
  //     hCaloPFpt[a][i]->SetXTitle("Calo jet pT");
  //     hCaloPFpt[a][i]->SetYTitle("Matched PF jet pT");

  //     hCaloPFCorr_pt[a][i] = new TH2F(Form("hCaloPFCorr_pt_%s_cent%d",TrigName[a],i), Form(" Ratio of PF/Calo vs calo pT %s %s", TrigName[a], centWidth[i]), nbins_pt,boundaries_pt, 200,0,10);
  //     hCaloPFCorr_pt[a][i]->SetXTitle("Calo jet pT");
  //     hCaloPFCorr_pt[a][i]->SetYTitle("PF Jet pT / Calo Jet PT");

  //     hCalo[a][i] = new TH1F(Form("hCalo_%s_cent%d",TrigName[a],i),Form("Calo jet spectra %s %s",TrigName[a],centWidth[i]),nbins_pt,boundaries_pt);
  //     hCalo[a][i]->SetXTitle("Jet p_{T} (GeV/c)");
  //     hCalo[a][i]->SetYTitle("Counts");

  //     hPF[a][i] = new TH1F(Form("hPF_%s_cent%d",TrigName[a],i),Form("PF jet spectra %s %s",TrigName[a],centWidth[i]),nbins_pt,boundaries_pt);
  //     hPF[a][i]->SetXTitle("Jet p_{T} (GeV/c)");
  //     hPF[a][i]->SetYTitle("Counts");

  //     hDeltaR_deltapT[a][i] = new TH2F(Form("hDeltaR_deltapT_%s_cent%d",TrigName[a],i), Form(" delta R for matched jets vs delta pT %s %s", TrigName[a], centWidth[i]), 200,0,1, 200,0,200);
  //     // hDeltaR[a][i] = new TH1F(Form("hPF_%s_cent%d",TrigName[a],i),Form("PF jet spectra %s %s",TrigName[a],centWidth[i]),nbins_pt,boundaries_pt);

  //   }
  // }

  cout<<"after histogram declaration"<<endl;

  Float_t calopt, pfpt, deltar, chMax, phMax, neMax, muMax, eMax, chSum, phSum, neSum, muSum, eSum, hcalSum, ecalSum;
  Int_t hiBin, jet80, jet80_prescl, jet65, jet65_prescl, jet55, jet55_prescl; 
  Int_t evt_value;
  Int_t run_value;
  Int_t lumi_value;
  Float_t subid, pfrawpt, calorawpt, pfrefpt, calorefpt, pfjtpu, calojtpu, vz, pthat;
  
  TTree* matchJets = new TTree("matchedJets","Ntuple containing important information about matched jets");
  matchJets->Branch("calopt",&calopt,"calopt/F");   matchJets->Branch("phSum",&phSum,"phSum/F");
  matchJets->Branch("pfpt",&pfpt,"pfpt/F");         matchJets->Branch("neSum",&neSum,"neSum/F");
  matchJets->Branch("deltar",&deltar,"deltar/F");   matchJets->Branch("muSum",&muSum,"muSum/F");
  matchJets->Branch("chMax",&chMax,"chMax/F");      matchJets->Branch("eSum",&eSum,"eSum/F");
  matchJets->Branch("phMax",&phMax,"phMax/F");      matchJets->Branch("hiBin",&hiBin,"hiBin/I");
  matchJets->Branch("neMax",&neMax,"neMax/F");      matchJets->Branch("jet80",&jet80,"jet80/I");
  matchJets->Branch("muMax",&muMax,"muMax/F");      matchJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  matchJets->Branch("eMax",&eMax,"eMax/F");         matchJets->Branch("jet65",&jet65,"jet65/I");
  matchJets->Branch("chSum",&chSum,"chSum/F");      matchJets->Branch("jet65_prescl",&jet65_prescl,"jet65_prescl/I");
  matchJets->Branch("jet55",&jet55,"jet55/I");      matchJets->Branch("jet55_prescl",&jet55_prescl,"jet55_prescl/I");
  matchJets->Branch("hcalSum",&hcalSum,"hcalSum/F");matchJets->Branch("ecalSum",&ecalSum,"ecalSum/F");
  matchJets->Branch("run_value",&run_value,"run_value/I");
  matchJets->Branch("evt_value",&evt_value,"evt_value/I");
  matchJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  if(dataset=="MC") matchJets->Branch("subid",&subid,"subid/I"); matchJets->Branch("pfrawpt",&pfrawpt,"pfrawpt/F");
  if(dataset=="MC") matchJets->Branch("pfrefpt",&pfrefpt,"pfrefpt/F"); matchJets->Branch("pfjtpu",&pfjtpu,"pfjtpu/F");
  if(dataset=="MC") matchJets->Branch("calorefpt",&calorefpt,"calorefpt/F");
  if(dataset=="MC") matchJets->Branch("pthat",&pthat,"pthat/F"); if(dataset=="MC") matchJets->Branch("vz",&vz,"vz/F");
  matchJets->Branch("calorawpt",&calorawpt,"calorawpt/F");
  matchJets->Branch("calojtpu",&calojtpu,"calojtpu/F");
  
  TTree* unmatchPFJets = new TTree("unmatchedPFJets","Ntuple containing important information about unmatched PF jets");
  unmatchPFJets->Branch("phSum",&phSum,"phSum/F");
  unmatchPFJets->Branch("pfpt",&pfpt,"pfpt/F");         unmatchPFJets->Branch("neSum",&neSum,"neSum/F");
  unmatchPFJets->Branch("deltar",&deltar,"deltar/F");   unmatchPFJets->Branch("muSum",&muSum,"muSum/F");
  unmatchPFJets->Branch("chMax",&chMax,"chMax/F");      unmatchPFJets->Branch("eSum",&eSum,"eSum/F");
  unmatchPFJets->Branch("phMax",&phMax,"phMax/F");      unmatchPFJets->Branch("hiBin",&hiBin,"hiBin/I");
  unmatchPFJets->Branch("neMax",&neMax,"neMax/F");      unmatchPFJets->Branch("jet80",&jet80,"jet80/I");
  unmatchPFJets->Branch("muMax",&muMax,"muMax/F");      unmatchPFJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  unmatchPFJets->Branch("eMax",&eMax,"eMax/F");         unmatchPFJets->Branch("jet65",&jet65,"jet65/I");
  unmatchPFJets->Branch("chSum",&chSum,"chSum/F");      unmatchPFJets->Branch("jet65_prescl",&jet65_prescl,"jet65_prescl/I");
  unmatchPFJets->Branch("jet55",&jet55,"jet55/I");      unmatchPFJets->Branch("jet55_prescl",&jet55_prescl,"jet55_prescl/I");
  unmatchPFJets->Branch("hcalSum",&hcalSum,"hcalSum/F");unmatchPFJets->Branch("ecalSum",&ecalSum,"ecalSum/F");
  unmatchPFJets->Branch("run_value",&run_value,"run_value/I");
  unmatchPFJets->Branch("evt_value",&evt_value,"evt_value/I");
  unmatchPFJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  if(dataset=="MC") unmatchPFJets->Branch("subid",&subid,"subid/I"); unmatchPFJets->Branch("pfrawpt",&pfrawpt,"pfrawpt/F");
  if(dataset=="MC") unmatchPFJets->Branch("pfrefpt",&pfrefpt,"pfrefpt/F"); unmatchPFJets->Branch("pfjtpu",&pfjtpu,"pfjtpu/F");
  if(dataset=="MC") unmatchPFJets->Branch("pthat",&pthat,"pthat/F"); if(dataset=="MC") unmatchPFJets->Branch("vz",&vz,"vz/F");
  
  TTree* unmatchCaloJets = new TTree("unmatchedCaloJets","Ntuple containing important information about unmatched Calo jets");
  unmatchCaloJets->Branch("phSum",&phSum,"phSum/F");
  unmatchCaloJets->Branch("calopt",&pfpt,"pfpt/F");       unmatchCaloJets->Branch("neSum",&neSum,"neSum/F");
  unmatchCaloJets->Branch("deltar",&deltar,"deltar/F");   unmatchCaloJets->Branch("muSum",&muSum,"muSum/F");
  unmatchCaloJets->Branch("chMax",&chMax,"chMax/F");      unmatchCaloJets->Branch("eSum",&eSum,"eSum/F");
  unmatchCaloJets->Branch("phMax",&phMax,"phMax/F");      unmatchCaloJets->Branch("hiBin",&hiBin,"hiBin/I");
  unmatchCaloJets->Branch("neMax",&neMax,"neMax/F");      unmatchCaloJets->Branch("jet80",&jet80,"jet80/I");
  unmatchCaloJets->Branch("muMax",&muMax,"muMax/F");      unmatchCaloJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  unmatchCaloJets->Branch("eMax",&eMax,"eMax/F");         unmatchCaloJets->Branch("jet65",&jet65,"jet65/I");
  unmatchCaloJets->Branch("chSum",&chSum,"chSum/F");      unmatchCaloJets->Branch("jet65_prescl",&jet65_prescl,"jet65_prescl/I");
  unmatchCaloJets->Branch("jet55",&jet55,"jet55/I");      unmatchCaloJets->Branch("jet55_prescl",&jet55_prescl,"jet55_prescl/I");
  unmatchCaloJets->Branch("hcalSum",&hcalSum,"hcalSum/F");unmatchCaloJets->Branch("ecalSum",&ecalSum,"ecalSum/F");
  unmatchCaloJets->Branch("run_value",&run_value,"run_value/I");
  unmatchCaloJets->Branch("evt_value",&evt_value,"evt_value/I");
  unmatchCaloJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  if(dataset=="MC") unmatchCaloJets->Branch("subid",&subid,"subid/I"); 
  if(dataset=="MC") unmatchCaloJets->Branch("calorefpt",&calorefpt,"calorefpt/F");
  if(dataset=="MC") unmatchCaloJets->Branch("pthat",&pthat,"pthat/F"); if(dataset=="MC") unmatchCaloJets->Branch("vz",&vz,"vz/F");
  unmatchCaloJets->Branch("calorawpt",&calorawpt,"calorawpt/F");
  unmatchCaloJets->Branch("calojtpu",&calojtpu,"calojtpu/F");

  
  //for the Jet 55 spectra: 
  Float_t effecPrescl = 2.047507;
  
  // declare the 2d calo and pf candidate vectors, deltaR_calovsPF is going to be a 3d vector like so:
  // calo jet on x axis, pf jet on y axis, and delta R, calopT, pfpT on z axis.
  // we also need to add in the trigger information here since we need to know what trigger the matching comes from. maybe i can run that later as an added check. for now just to see if the algorithm is working should try to run things. 
  //vector<vector<double> > caloJet;
  //vector<vector<double> > pfJet;

  // start the loop process.

  for(int k = 0;k<no_radius;k++){

    Long64_t nentries = jetpbpb1[2][k]->GetEntries();
    //nentries = 20;
    
    for(Long64_t nentry = 0; nentry<nentries;nentry++){
      if(printDebug)cout<<"event no = "<<nentry<<endl;
      for(int t = 0;t<N;t++)  jetpbpb1[t][k]->GetEntry(nentry);
      
      int centBin = findBin(hiBin_1);
      if(centBin==-1) continue;
      //if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 
      if(pcollisionEventSelection_1==0) continue; 
      if(dataset=="Data" || dataset=="MinBiasUPC") if(pcollisionEventSelection_1==0) continue;

      if(fabs(vz_1)>15) continue;
      //cout<<"passed the selection"<<endl;

      int jetCounter = 0;//counts jets which are going to be used in the supernova cut rejection. 

      for(int g = 0;g<nrefe_1;g++){
	    
	if(eta_1[g]>=-2 && eta_1[g]<2){ //to select inside 
	      
	  if(pt_1[g]>=50) jetCounter++;
	  
	}//eta selection cut
	
      }// jet loop
      
      // apply the correct supernova selection cut rejection here: 
      if(hiNpix_1 > 38000 - 500*jetCounter){
       	if(printDebug) cout<<"removed this supernova event"<<endl;
      	continue;
      }

      // if(jet80_1!=0)cout<<"jet80 = "<<jet80_1<<" jet80_prescl = "<<jet80_p_1<<endl;
      // if(jet65_1!=0)cout<<"jet65 = "<<jet65_1<<" jet65_prescl = "<<jet65_p_1<<endl;
      // if(jet55_1!=0)cout<<"jet55 = "<<jet55_1<<" jet55_prescl = "<<jet55_p_1<<endl;

      // start doing the search for the match. - best thing to do would be to create a 2D match delta R matrix with each calo jet and pf jet. Once thats done - find the smallest entry in that matrix. the i,j of that smallest entry are matched. now remove the row i and column j and then we have a new distance matrix. where we need to find the smallest element again. keep doing this till we have either no rows or no columns.  
      // declare the necessary variables:
      Float_t deltaRCaloPF = 0;
      Float_t calojet_eta = 0;
      Float_t calojet_phi = 0;
      Float_t calojet_pt = 0;
      Float_t pfjet_eta = 0;
      Float_t pfjet_phi = 0;
      Float_t pfjet_pt = 0;
      
      vector<vector<double> > matchedCaloPFJet;
      vector<vector<vector<double> > > deltaR_calovsPF;
      int calosize = 0;
      Float_t deltapT = 0;
      int calomatchcounter = 0;
      
      for(int g = 0;g<nrefe_1;g++){
	
	calojet_eta = eta_1[g];
	calojet_phi = phi_1[g];
	calojet_pt = pt_1[g];
	
	if(calojet_pt < CALOPTCUT) continue;
	if(TMath::Abs(calojet_eta) > 2.0) continue;
	
	int pfmatchcounter = 0;
	deltaR_calovsPF.push_back(vector<vector<double> > ());
	
	calosize++;
	
	for(int j = 0;j<nrefe_2;j++){

	  pfjet_eta = eta_2[j];
	  pfjet_phi = phi_2[j];
	  pfjet_pt = pt_2[j];

	  if(pfjet_pt < PFPTCUT) continue;
	  if(TMath::Abs(pfjet_eta) > 2.0) continue;

	  deltaRCaloPF = Calc_deltaR(calojet_eta, calojet_phi, pfjet_eta, pfjet_phi);
	  
	  // if(deltaRCaloPF > (Float_t)deltaR/10) continue;
	  deltapT = TMath::Abs(calojet_pt - pfjet_pt);

	  deltaR_calovsPF[calomatchcounter].push_back(vector<double> ());
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(deltaRCaloPF); // 0 - delta R
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(deltapT); // 1 - delta pT	  
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(g); // 2 - calo counter 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(calojet_pt); // 3 - calo jet pT 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(j); // 4 - pf counter 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(pfjet_pt); // 5 - pf jet pT
	  // this will have the candidate variables for the matched jets. 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(chMax_2[j]); // 6 - chMax 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(phMax_2[j]); // 7 - phMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(neMax_2[j]); // 8 - neMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(muMax_2[j]); // 9 - muMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(eMax_2[j]); // 10 - eMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(chSum_2[j]); // 11 - chSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(phSum_2[j]); // 12 - phSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(neSum_2[j]); // 13 - neSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(muSum_2[j]); // 14 - muSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(eSum_2[j]); // 15 - eSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(0.0); // 16 - this is the variable which will tell me if a jet is matched.
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(hcalSum_2[j]); // 17 - hcalSum 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(ecalSum_2[j]); // 18 - hcalSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(raw_2[j]); // 19 - raw pt
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(jtpu_2[j]); // 20 - jtpu
	  //all the following candidate information is for the calo jets 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(chMax_1[g]); // 21 - chMax 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(phMax_1[g]); // 22 - phMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(neMax_1[g]); // 23 - neMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(muMax_1[g]); // 24 - muMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(eMax_1[g]); // 25 - eMax
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(chSum_1[g]); // 26 - chSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(phSum_1[g]); // 27 - phSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(neSum_1[g]); // 28 - neSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(muSum_1[g]); // 29 - muSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(eSum_1[g]); // 30 - eSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(hcalSum_1[g]); // 31 - hcalSum 
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(ecalSum_1[g]); // 32 - hcalSum
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(raw_1[g]); // 33 - raw pt
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(jtpu_1[g]); // 34 - jtpu

	  if(dataset=="MC") {
	    deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(refpt_1[g]); // 35 calo ref pt 
	    deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(refpt_2[j]); // 36 pf ref pt 
	    deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(subid_2[j]); // 37 - subid only for MC
	    
	  }
	  
	  ++pfmatchcounter;

	}// pf jet loop
	++calomatchcounter;
	
      }// calo jet loop

      // set the event variables before proceeding to the matching. 
      hiBin = hiBin_1;
      jet80 = jet80_1;
      jet80_prescl = jet80_p_1;
      jet65 = jet65_1;
      jet65_prescl = jet65_p_1;
      jet55 = jet55_1;
      jet55_prescl = jet55_p_1;
      run_value = run_1;
      evt_value = evt_1;
      lumi_value = lumi_1;
      vz = vz_1;
      if(dataset=="MC")pthat = pthat_1;

      if(printDebug)cout<<deltaR_calovsPF.size()<<endl;
      if(printDebug)for(int a = 0;a<deltaR_calovsPF.size();++a) cout<<deltaR_calovsPF[a].size()<<" ";
      if(printDebug)cout<<endl<<"going to small matching"<<endl;
      
      // now that we have the 2D matrix, lets find the smallest delta R element from that and fill in the value of the 

      Float_t smallDeltaR = 10;
      Int_t small_calo = 0;
      Int_t small_pf = 0;

      for(int c = 0;c<deltaR_calovsPF.size();++c){
	// if(printDebug)cout<<"going through all rows in the  matrix "<<c<<endl;
	for(int a = 0;a<deltaR_calovsPF.size();++a){
	  // if(printDebug)cout<<"calo jet iteration "<<a<<endl;
	  for(int b = 0;b<deltaR_calovsPF[a].size();++b){
	    // if(printDebug)cout<<"pf jet iteration "<<b<<endl;
	    
	    if(deltaR_calovsPF[a][b][16]==1){ 
	      break;
	    }
	    if(smallDeltaR > deltaR_calovsPF[a][b][0]){
	      
	      smallDeltaR = deltaR_calovsPF[a][b][0];
	      small_calo = a;
	      small_pf = b;
	      
	    }
	  }
	}
	
	if(smallDeltaR > (Float_t)deltaR/10 || deltaR_calovsPF[small_calo][small_pf][16] == 1) continue;
	
	calopt = deltaR_calovsPF[small_calo][small_pf][3];
	pfpt = deltaR_calovsPF[small_calo][small_pf][5];
	if(dataset=="MC"){
	  pfrefpt = deltaR_calovsPF[small_calo][small_pf][36];
	  calorefpt = deltaR_calovsPF[small_calo][small_pf][35];
	  subid = deltaR_calovsPF[small_calo][small_pf][37]; 
	}
	pfrawpt = deltaR_calovsPF[small_calo][small_pf][19]; 
	calorawpt = deltaR_calovsPF[small_calo][small_pf][33]; 
	pfjtpu = deltaR_calovsPF[small_calo][small_pf][20]; 
	calojtpu = deltaR_calovsPF[small_calo][small_pf][34];
	deltar = deltaR_calovsPF[small_calo][small_pf][0];
	chMax = deltaR_calovsPF[small_calo][small_pf][6];
	phMax = deltaR_calovsPF[small_calo][small_pf][7];
	neMax = deltaR_calovsPF[small_calo][small_pf][8];
	muMax = deltaR_calovsPF[small_calo][small_pf][9];
	eMax = deltaR_calovsPF[small_calo][small_pf][10];
	chSum = deltaR_calovsPF[small_calo][small_pf][11];
	phSum = deltaR_calovsPF[small_calo][small_pf][12];
	neSum = deltaR_calovsPF[small_calo][small_pf][13];
	muSum = deltaR_calovsPF[small_calo][small_pf][14];
	eSum = deltaR_calovsPF[small_calo][small_pf][15];

	hcalSum = deltaR_calovsPF[small_calo][small_pf][17];
	ecalSum = deltaR_calovsPF[small_calo][small_pf][18];
	
	matchJets->Fill();
	
	//now we have the smallest value, lets get the delta R of that particular matched jets 
	
	// if(jet80_1){
	//   hCaloPFCorr[2][centBin]->Fill((Float_t) pfjet_pt/calojet_pt);
	//   hCaloPFpt[2][centBin]->Fill(calojet_pt,pfjet_pt);
	//   hCaloPFCorr_pt[2][centBin]->Fill(calojet_pt,(Float_t)pfjet_pt/calojet_pt);
	//   hCalo[2][centBin]->Fill(calojet_pt);
	//   hPF[2][centBin]->Fill(pfjet_pt);	    
	//   hDeltaR_deltapT[2][centBin]->Fill(deltaRCaloPF,deltapT);
	// }
	// if(jet65_1 && !jet80_1){
	//   hCaloPFCorr[1][centBin]->Fill((Float_t) pfjet_pt/calojet_pt);
	//   hCaloPFpt[1][centBin]->Fill(calojet_pt,pfjet_pt);
	//   hCaloPFCorr_pt[1][centBin]->Fill(calojet_pt,(Float_t)pfjet_pt/calojet_pt);
	//   hCalo[1][centBin]->Fill(calojet_pt);
	//   hPF[1][centBin]->Fill(pfjet_pt);	    
	//   hDeltaR_deltapT[1][centBin]->Fill(deltaRCaloPF,deltapT);
	// }
	// if(jet55_1 && !jet65_1 && !jet80_1){
	//   hCaloPFCorr[0][centBin]->Fill((Float_t) pfjet_pt/calojet_pt);
	//   hCaloPFpt[0][centBin]->Fill(calojet_pt,pfjet_pt);
	//   hCaloPFCorr_pt[0][centBin]->Fill(calojet_pt,(Float_t)pfjet_pt/calojet_pt);
	//   hCalo[0][centBin]->Fill(calojet_pt);
	//   hPF[0][centBin]->Fill(pfjet_pt);	    
	//   hDeltaR_deltapT[0][centBin]->Fill(deltaRCaloPF,deltapT);
	// }
	
	smallDeltaR = 10;
	for(int b = 0;b<deltaR_calovsPF[small_calo].size();++b){
	  deltaR_calovsPF[small_calo][b][16] = 1; // by setting this value you effectively remove that calo and pf jet for further matching // removes the column for getting matched. 
	}
	for(int a = 0;a<deltaR_calovsPF.size();++a){
	  deltaR_calovsPF[a][small_pf][16] = 1; // this removes the whole row from getting matched later. 
	}

      }// running it for the number of calo jets: 

      if(printDebug)cout<<"now going to find unmatched pf jets"<<endl;
      // ok Now lets find the un-matched jets and fill the necessary unmatched ntuple:

      if(deltaR_calovsPF.size() == 0) continue;
      for(int b = 0;b<deltaR_calovsPF[0].size();++b){
	if(printDebug)cout<<"pf jet iteration "<<b<<endl;

	if(deltaR_calovsPF[0][b][16] == 1) continue;

	pfpt = deltaR_calovsPF[0][b][5];
	deltar = deltaR_calovsPF[0][b][0];
	if(dataset=="MC"){
	  pfrefpt = deltaR_calovsPF[0][b][36];
	  subid = deltaR_calovsPF[0][b][37]; 
	}
	pfrawpt = deltaR_calovsPF[0][b][19]; 
	pfjtpu = deltaR_calovsPF[0][b][20]; 
	chMax = deltaR_calovsPF[0][b][6];
	phMax = deltaR_calovsPF[0][b][7];
	neMax = deltaR_calovsPF[0][b][8];
	muMax = deltaR_calovsPF[0][b][9];
	eMax = deltaR_calovsPF[0][b][10];
	chSum = deltaR_calovsPF[0][b][11];
	phSum = deltaR_calovsPF[0][b][12];
	neSum = deltaR_calovsPF[0][b][13];
	muSum = deltaR_calovsPF[0][b][14];
	eSum = deltaR_calovsPF[0][b][15];

	hcalSum = deltaR_calovsPF[0][b][17];
	ecalSum = deltaR_calovsPF[0][b][18];
	
	unmatchPFJets->Fill();
	
      }// unmatched PF jets
      
      if(printDebug)cout<<"now going to find unmatched calo jets"<<endl;

      if(deltaR_calovsPF[0].size() == 0) continue;
      
      for(int a = 0;a<deltaR_calovsPF.size();++a){
	if(printDebug)cout<<"calo jet iteration "<<a<<endl;

	if(deltaR_calovsPF[a][0][16] == 1) continue;

	calopt = deltaR_calovsPF[a][0][3];
	deltar = deltaR_calovsPF[a][0][0];
	if(dataset=="MC"){
	  calorefpt = deltaR_calovsPF[a][0][35];
	  subid = deltaR_calovsPF[a][0][37]; 
	}
	calorawpt = deltaR_calovsPF[a][0][33]; 
	calojtpu = deltaR_calovsPF[a][0][34]; 
	chMax = deltaR_calovsPF[a][0][21];
	phMax = deltaR_calovsPF[a][0][22];
	neMax = deltaR_calovsPF[a][0][23];
	muMax = deltaR_calovsPF[a][0][24];
	eMax = deltaR_calovsPF[a][0][25];
	chSum = deltaR_calovsPF[a][0][26];
	phSum = deltaR_calovsPF[a][0][27];
	neSum = deltaR_calovsPF[a][0][28];
	muSum = deltaR_calovsPF[a][0][29];
	eSum = deltaR_calovsPF[a][0][30];

	hcalSum = deltaR_calovsPF[a][0][31];
	ecalSum = deltaR_calovsPF[a][0][32];
	
	unmatchCaloJets->Fill();

      }

    }// event loop 

  }// radius loop

  //TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_data_calo_pf_jet_correlation_deltaR_0p%d_ak%s%d_%d_%d.root",deltaR,algo,radius,date.GetDate(),endfile),"RECREATE");

  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_%s_calo_pf_jet_correlation_deltaR_0p%d_ak%s%d_%d_%d.root",dataset,deltaR,algo,radius,date.GetDate(),endfile),"RECREATE");
  f.cd();

  matchJets->Write();
  matchJets->Print();
  unmatchPFJets->Write();
  unmatchPFJets->Print();
  unmatchCaloJets->Write();
  unmatchCaloJets->Print();
  //unmatch_CaloJets->Write();
  //unmatch_CaloJets->Print();
  
  // for(int i = 0;i<nbins_cent;++i){

  //   for(int a = 0;a<TrigValue-1;++a){

  //     hCaloPFCorr[a][i]->Write();
  //     hCaloPFCorr[a][i]->Print();
  //     hCalo[a][i]->Write();
  //     hCalo[a][i]->Print();
  //     hPF[a][i]->Write();
  //     hPF[a][i]->Print();
  //     hCaloPFpt[a][i]->Write();
  //     hCaloPFpt[a][i]->Print();
  //     hCaloPFCorr_pt[a][i]->Write();
  //     hCaloPFCorr_pt[a][i]->Print();
  //     hDeltaR_deltapT[a][i]->Write();
  //     hDeltaR_deltapT[a][i]->Print();
  //   }
    
  // }

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}
