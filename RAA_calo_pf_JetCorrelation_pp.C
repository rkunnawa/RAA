//
// Raghav Kunnawalkam Elayavalli
// March 26th 2014
// Rutgers, raghav.k.e at CERN dot CH
//
// modified version of the current macro which does this matching for the pp, only need to change a few things 
// Macro to do calo-pf jet correlation (spacial) to study the PF electron issues - this would tell us if PF jets from a jet55 trigger has a calo jet as a background. following is how to do that: 
//load the forest data files. and then load the calo jets first followed by the pf jet to search for spacial correlation. Once we find the pf jet which is in the same spacial location (which is found by calculating the delta R for all pf jets for one calo jet and the jet which is less than the value of delta R given above becomes it correlated jet), the ratio of the pT of the calo jet and the pf jet (correlated) is calculated. and that distribution is expected to be a gaussian with the outliers clearly saying there is something wrong with the candidates like electron etc...  
//


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
#include <TEventList.h>
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

static const int nbins_pt = 31;
static const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 272, 300, 
  330, 362, 395
};

static const int no_radius = 1;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {3};
static const int TrigValue = 4;
static const char TrigName [TrigValue][256] = {"HLT40","HLT60","HLT80","FullDataset"};


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

void RAA_calo_pf_JetCorrelation_pp(int startfile = 10, int endfile = 11, int radius=2, int deltaR=2/*which i will divide by 10 later when using*/, Float_t CALOPTCUT = 30.0, Float_t PFPTCUT = 30.0, char *dataset = "MC", char *etaWidth = "n20_eta_p20"){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  timer.Start();

  TDatime date;
  
  bool printDebug = false;

  // Since this has to run on the HiForest files, it is best to run them as condor jobs similar to the Pp data read macro, along with the jet trees which get their branch address set. 
  
  std::string infile1;
  if(dataset == "Data")infile1 = "jetRAA_pp_data_forest.txt";
  if(dataset == "MC")infile1 = "jetRAA_pp_mc_forest.txt";

  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;

  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;

  if(dataset=="MC"){
    for(int ifile = 0;ifile<4*startfile;ifile++){
      instr1>>filename1;
    }
  }
  if(dataset=="Data"){
    for(int ifile = 0;ifile<startfile;ifile++){
      instr1>>filename1;
    }
  }

  
  const int N = 5;
  
  TChain *jetpp1[N][no_radius];

  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%dCaloJetAnalyzer",3);
    dir[3][k] = Form("ak%dPFJetAnalyzer",radius);
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
      jetpp1[t][k] = new TChain(string(dir[t][k]+"/"+trees[t]).data());
    }//tree loop ends
  }// radius loop ends
  
  if(dataset == "Data"){
    for(int ifile = startfile;ifile<endfile;ifile++){
    
      instr1>>filename1;
      if(printDebug)cout<<"File: "<<filename1<<endl;
      for(int k = 0;k<no_radius;k++){

	for(int t = 0;t<N;t++){
	
	  jetpp1[t][k]->Add(filename1.c_str());
	  if(printDebug)cout << "Tree loaded  " << string(dir[t][k]+"/"+trees[t]).data() << endl;
	  if(printDebug)cout << "Entries : " << jetpp1[t][k]->GetEntries() << endl;
		
	}// tree loop ends
      
      }// radius loop ends
    
    }// file loop ends

  }

  // these different file loops for MC is necessary because of eachfile being split into 4 files except the last pthat 540 which has 5 files. 
  if(dataset == "MC" && startfile <10){
    for(int ifile = 4*startfile;ifile<4*endfile;ifile++){
    
      instr1>>filename1;
      if(printDebug)cout<<"File: "<<filename1<<endl;
      for(int k = 0;k<no_radius;k++){

	for(int t = 0;t<N;t++){
	
	  jetpp1[t][k]->Add(filename1.c_str());
	  if(printDebug)cout << "Tree loaded  " << string(dir[t][k]+"/"+trees[t]).data() << endl;
	  if(printDebug)cout << "Entries : " << jetpp1[t][k]->GetEntries() << endl;
		
	}// tree loop ends
      
      }// radius loop ends
    
    }// file loop ends

  }

  if(dataset == "MC" && startfile == 10){
    for(int ifile = 4*startfile;ifile<4*endfile+1;ifile++){
    
      instr1>>filename1;
      if(printDebug)cout<<"File: "<<filename1<<endl;
      for(int k = 0;k<no_radius;k++){

	for(int t = 0;t<N;t++){
	
	  jetpp1[t][k]->Add(filename1.c_str());
	  if(printDebug)cout << "Tree loaded  " << string(dir[t][k]+"/"+trees[t]).data() << endl;
	  if(printDebug)cout << "Entries : " << jetpp1[t][k]->GetEntries() << endl;
		
	}// tree loop ends
      
      }// radius loop ends
    
    }// file loop ends

  }
  
  for(int k = 0;k<no_radius;k++){
    jetpp1[2][k]->AddFriend(jetpp1[0][k]);
    jetpp1[2][k]->AddFriend(jetpp1[1][k]);
    jetpp1[2][k]->AddFriend(jetpp1[4][k]);
    //jetpp1[2][k]->AddFriend(jetpp1[5][k]);
    
    jetpp1[3][k]->AddFriend(jetpp1[0][k]);
    jetpp1[3][k]->AddFriend(jetpp1[1][k]);
    jetpp1[3][k]->AddFriend(jetpp1[4][k]);
    //jetpp1[3][k]->AddFriend(jetpp1[5][k]);
    
  }// radius loop ends
 
  // Ok this should now work. 

  cout<<"total no of entries in the combined forest files = "<<jetpp1[2][0]->GetEntries()<<endl;
  //  if(printDebug)cout<<"total no of entries in the Jet80 Tree     = "<<jetpp2[2][0]->GetEntries()<<endl;

  // get the centrality weight, vz weight and the scale from the cross section, pt weighting. 
  static const Int_t nbins_pthat = 11;
  Double_t xsection[nbins_pthat+1] = {2.034e-01, 1.075e-02, 1.025e-03,  9.865e-05, 1.129e-05, 1.465e-06, 2.837e-07, 5.323e-08, 5.934e-09, 8.125e-10, 1.468e-10, 0};
  Double_t boundaries_pthat[nbins_pthat+1] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 2000};
  TH1F * hpthatBin = new TH1F("hpthatBin","",nbins_pthat, boundaries_pthat);

  // Vertex reweighting for pp
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
  
  
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/pp_%s_closure_histogram_noJetID_deltaR_0p%d_ak%d_%d_%d.root",dataset,deltaR,radius,date.GetDate(),endfile),"RECREATE");
  f.cd();


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
  int pPAcollisionEventSelectionPA_1;
  int pHBHENoiseFilter_1;
  int pprimaryvertexFilter_1;
  int pVertexFilterCutGplus_1;

  // trigger tree
  int jet40_1;
  int jet60_1;
  int jet80_1;
  int jet40_p_1;
  int jet60_p_1;
  int jet80_p_1;

  for(int k = 0;k<no_radius;k++){
    //set the branch addresses:  - one of the most boring parts of the code: 
    jetpp1[2][k]->SetBranchAddress("evt",&evt_1);
    jetpp1[2][k]->SetBranchAddress("run",&run_1);
    jetpp1[2][k]->SetBranchAddress("lumi",&lumi_1);
    if(dataset=="MC") jetpp1[2][k]->SetBranchAddress("pthat",&pthat_1);
    jetpp1[2][k]->SetBranchAddress("vz",&vz_1);
    jetpp1[2][k]->SetBranchAddress("vx",&vx_1);
    jetpp1[2][k]->SetBranchAddress("vy",&vy_1);
    jetpp1[2][k]->SetBranchAddress("hiNpix",&hiNpix_1);
    jetpp1[2][k]->SetBranchAddress("hiNtracks",&hiNtracks_1);
    jetpp1[2][k]->SetBranchAddress("hiHFminus",&hiHFminus_1);
    jetpp1[2][k]->SetBranchAddress("hiHF",&hiHF_1);
    jetpp1[2][k]->SetBranchAddress("hiHFplus",&hiHFplus_1);
    jetpp1[2][k]->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_1);
    jetpp1[2][k]->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_1);
    jetpp1[2][k]->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA_1);
    jetpp1[2][k]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_1);
    //jetpp1[2][k]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_1);
    //jetpp1[2][k]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_1);
  
    jetpp1[2][k]->SetBranchAddress("nref",&nrefe_1);
    jetpp1[2][k]->SetBranchAddress("jtpt",&pt_1);
    jetpp1[2][k]->SetBranchAddress("jteta",&eta_1);
    jetpp1[2][k]->SetBranchAddress("jtphi",&phi_1);
    jetpp1[2][k]->SetBranchAddress("rawpt",&raw_1);
    jetpp1[2][k]->SetBranchAddress("jtpu",&jtpu_1);
    jetpp1[2][k]->SetBranchAddress("chargedMax",&chMax_1);
    jetpp1[2][k]->SetBranchAddress("chargedSum",&chSum_1);
    jetpp1[2][k]->SetBranchAddress("trackMax",&trkMax_1);
    jetpp1[2][k]->SetBranchAddress("trackSum",&trkSum_1);
    jetpp1[2][k]->SetBranchAddress("photonMax",&phMax_1);
    jetpp1[2][k]->SetBranchAddress("photonSum",&phSum_1);
    jetpp1[2][k]->SetBranchAddress("neutralMax",&neMax_1);
    jetpp1[2][k]->SetBranchAddress("neutralSum",&neSum_1);
    jetpp1[2][k]->SetBranchAddress("eSum",&eSum_1);
    jetpp1[2][k]->SetBranchAddress("eMax",&eMax_1);
    jetpp1[2][k]->SetBranchAddress("muSum",&muSum_1);
    jetpp1[2][k]->SetBranchAddress("muMax",&muMax_1);
    jetpp1[2][k]->SetBranchAddress("ecalSum",&ecalSum_1);
    jetpp1[2][k]->SetBranchAddress("hcalSum",&hcalSum_1);
    
    jetpp1[3][k]->SetBranchAddress("nref",&nrefe_2);
    jetpp1[3][k]->SetBranchAddress("jtpt",&pt_2);
    jetpp1[3][k]->SetBranchAddress("jteta",&eta_2);
    jetpp1[3][k]->SetBranchAddress("jtphi",&phi_2);
    jetpp1[3][k]->SetBranchAddress("rawpt",&raw_2);
    jetpp1[3][k]->SetBranchAddress("jtpu",&jtpu_2);
    jetpp1[3][k]->SetBranchAddress("chargedMax",&chMax_2);
    jetpp1[3][k]->SetBranchAddress("chargedSum",&chSum_2);
    jetpp1[3][k]->SetBranchAddress("trackMax",&trkMax_2);
    jetpp1[3][k]->SetBranchAddress("trackSum",&trkSum_2);
    jetpp1[3][k]->SetBranchAddress("photonMax",&phMax_2);
    jetpp1[3][k]->SetBranchAddress("photonSum",&phSum_2);
    jetpp1[3][k]->SetBranchAddress("neutralMax",&neMax_2);
    jetpp1[3][k]->SetBranchAddress("neutralSum",&neSum_2);
    jetpp1[3][k]->SetBranchAddress("eSum",&eSum_2);
    jetpp1[3][k]->SetBranchAddress("eMax",&eMax_2);
    jetpp1[3][k]->SetBranchAddress("muSum",&muSum_2);
    jetpp1[3][k]->SetBranchAddress("muMax",&muMax_2);
    jetpp1[3][k]->SetBranchAddress("ecalSum",&ecalSum_2);
    jetpp1[3][k]->SetBranchAddress("hcalSum",&hcalSum_2);
    
    if(dataset=="MC"){

      jetpp1[2][k]->SetBranchAddress("subid",&subid_1);
      jetpp1[3][k]->SetBranchAddress("subid",&subid_2);
      jetpp1[2][k]->SetBranchAddress("refpt",&refpt_1);
      jetpp1[3][k]->SetBranchAddress("refpt",&refpt_2);

    }
    

    jetpp1[2][k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40_1);
    jetpp1[2][k]->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_p_1);
    jetpp1[2][k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60_1);
    jetpp1[2][k]->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_p_1);
    jetpp1[2][k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80_1);
    jetpp1[2][k]->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_p_1);

    // jetpp1[2][k]->SetBranchAddress("id",&trgObj_id_1);
    // jetpp1[2][k]->SetBranchAddress("pt",&trgObj_pt_1);
    // jetpp1[2][k]->SetBranchAddress("eta",&trgObj_eta_1);
    // jetpp1[2][k]->SetBranchAddress("phi",&trgObj_phi_1);
    // jetpp1[2][k]->SetBranchAddress("mass",&trgObj_mass_1);
    
    /*
    jetpp1[2][k]->SetBranchAddress("nPFpart", &nPFpart);
    jetpp1[2][k]->SetBranchAddress("pfId", pfId);
    jetpp1[2][k]->SetBranchAddress("pfPt", pfPt);
    jetpp1[2][k]->SetBranchAddress("pfVsPtInitial", pfVsPtInitial);
    jetpp1[2][k]->SetBranchAddress("pfVsPt", pfVsPt);
    jetpp1[2][k]->SetBranchAddress("pfEta", pfEta);
    jetpp1[2][k]->SetBranchAddress("pfPhi", pfPhi);
    jetpp1[2][k]->SetBranchAddress("pfArea", pfArea);
    jetpp1[2][k]->SetBranchAddress("vn",&v_n);
    jetpp1[2][k]->SetBranchAddress("psin",&psi_n);
    jetpp1[2][k]->SetBranchAddress("sumpt",&sumpT);
    */

  }// radius loop
  cout<<"after branch declaration"<<endl;

  cout<<"after histogram declaration"<<endl;

  Float_t calopt, pfpt, deltar, chMax, phMax, neMax, muMax, eMax, chSum, phSum, neSum, muSum, eSum, hcalSum, ecalSum;
  Int_t  jet80, jet80_prescl, jet60, jet60_prescl, jet40, jet40_prescl; 
  Int_t evt_value;
  Int_t run_value;
  Int_t lumi_value;
  Float_t weight;
  Float_t subid, pfrawpt, calorawpt, pfrefpt, calorefpt, pfjtpu, calojtpu, vz, pthat, caloeta, calophi, pfeta, pfphi;
  
  TTree* matchJets = new TTree("matchedJets","Ntuple containing important information about matched jets");
  matchJets->Branch("calopt",&calopt,"calopt/F");   matchJets->Branch("phSum",&phSum,"phSum/F");
  matchJets->Branch("pfpt",&pfpt,"pfpt/F");         matchJets->Branch("neSum",&neSum,"neSum/F");
  matchJets->Branch("deltar",&deltar,"deltar/F");   matchJets->Branch("muSum",&muSum,"muSum/F");
  matchJets->Branch("chMax",&chMax,"chMax/F");      matchJets->Branch("eSum",&eSum,"eSum/F");
  matchJets->Branch("phMax",&phMax,"phMax/F"); 
  matchJets->Branch("neMax",&neMax,"neMax/F");      matchJets->Branch("jet80",&jet80,"jet80/I");
  matchJets->Branch("muMax",&muMax,"muMax/F");      matchJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  matchJets->Branch("eMax",&eMax,"eMax/F");         matchJets->Branch("jet60",&jet60,"jet60/I");
  matchJets->Branch("chSum",&chSum,"chSum/F");      matchJets->Branch("jet60_prescl",&jet60_prescl,"jet60_prescl/I");
  matchJets->Branch("jet40",&jet40,"jet40/I");      matchJets->Branch("jet40_prescl",&jet40_prescl,"jet40_prescl/I");
  matchJets->Branch("hcalSum",&hcalSum,"hcalSum/F");matchJets->Branch("ecalSum",&ecalSum,"ecalSum/F");
  matchJets->Branch("run_value",&run_value,"run_value/I");
  matchJets->Branch("evt_value",&evt_value,"evt_value/I");
  matchJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  matchJets->Branch("caloeta", &caloeta, "caloeta/F");
  matchJets->Branch("calophi", &calophi, "calophi/F");
  matchJets->Branch("pfeta", &pfeta, "pfeta/F");
  matchJets->Branch("pfphi", &pfphi, "pfphi/F");
  if(dataset=="MC") matchJets->Branch("subid",&subid,"subid/I"); matchJets->Branch("pfrawpt",&pfrawpt,"pfrawpt/F");
  if(dataset=="MC") matchJets->Branch("pfrefpt",&pfrefpt,"pfrefpt/F"); matchJets->Branch("pfjtpu",&pfjtpu,"pfjtpu/F");
  if(dataset=="MC") matchJets->Branch("calorefpt",&calorefpt,"calorefpt/F");
  if(dataset=="MC") matchJets->Branch("pthat",&pthat,"pthat/F"); if(dataset=="MC") matchJets->Branch("vz",&vz,"vz/F");
  if(dataset=="MC") matchJets->Branch("weight",&weight,"weight/F");
  matchJets->Branch("calorawpt",&calorawpt,"calorawpt/F");
  matchJets->Branch("calojtpu",&calojtpu,"calojtpu/F");
  
  TTree* unmatchPFJets = new TTree("unmatchedPFJets","Ntuple containing important information about unmatched PF jets");
  unmatchPFJets->Branch("phSum",&phSum,"phSum/F");
  unmatchPFJets->Branch("pfpt",&pfpt,"pfpt/F");         unmatchPFJets->Branch("neSum",&neSum,"neSum/F");
  unmatchPFJets->Branch("deltar",&deltar,"deltar/F");   unmatchPFJets->Branch("muSum",&muSum,"muSum/F");
  unmatchPFJets->Branch("chMax",&chMax,"chMax/F");      unmatchPFJets->Branch("eSum",&eSum,"eSum/F");
  unmatchPFJets->Branch("phMax",&phMax,"phMax/F");
  unmatchPFJets->Branch("neMax",&neMax,"neMax/F");      unmatchPFJets->Branch("jet80",&jet80,"jet80/I");
  unmatchPFJets->Branch("muMax",&muMax,"muMax/F");      unmatchPFJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  unmatchPFJets->Branch("eMax",&eMax,"eMax/F");         unmatchPFJets->Branch("jet60",&jet60,"jet60/I");
  unmatchPFJets->Branch("chSum",&chSum,"chSum/F");      unmatchPFJets->Branch("jet60_prescl",&jet60_prescl,"jet60_prescl/I");
  unmatchPFJets->Branch("jet40",&jet40,"jet40/I");      unmatchPFJets->Branch("jet40_prescl",&jet40_prescl,"jet40_prescl/I");
  unmatchPFJets->Branch("hcalSum",&hcalSum,"hcalSum/F");unmatchPFJets->Branch("ecalSum",&ecalSum,"ecalSum/F");
  unmatchPFJets->Branch("run_value",&run_value,"run_value/I");
  unmatchPFJets->Branch("evt_value",&evt_value,"evt_value/I");
  unmatchPFJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  unmatchPFJets->Branch("pfeta", &pfeta, "pfeta/F");
  unmatchPFJets->Branch("pfphi", &pfphi, "pfphi/F");
  if(dataset=="MC") unmatchPFJets->Branch("subid",&subid,"subid/I"); unmatchPFJets->Branch("pfrawpt",&pfrawpt,"pfrawpt/F");
  if(dataset=="MC") unmatchPFJets->Branch("pfrefpt",&pfrefpt,"pfrefpt/F"); unmatchPFJets->Branch("pfjtpu",&pfjtpu,"pfjtpu/F");
  if(dataset=="MC") unmatchPFJets->Branch("weight",&weight,"weight/F");
  if(dataset=="MC") unmatchPFJets->Branch("pthat",&pthat,"pthat/F"); if(dataset=="MC") unmatchPFJets->Branch("vz",&vz,"vz/F");
  
  TTree* unmatchCaloJets = new TTree("unmatchedCaloJets","Ntuple containing important information about unmatched Calo jets");
  unmatchCaloJets->Branch("phSum",&phSum,"phSum/F");
  unmatchCaloJets->Branch("calopt",&calopt,"calopt/F");   unmatchCaloJets->Branch("neSum",&neSum,"neSum/F");
  unmatchCaloJets->Branch("deltar",&deltar,"deltar/F");   unmatchCaloJets->Branch("muSum",&muSum,"muSum/F");
  unmatchCaloJets->Branch("chMax",&chMax,"chMax/F");      unmatchCaloJets->Branch("eSum",&eSum,"eSum/F");
  unmatchCaloJets->Branch("phMax",&phMax,"phMax/F");     
  unmatchCaloJets->Branch("neMax",&neMax,"neMax/F");      unmatchCaloJets->Branch("jet80",&jet80,"jet80/I");
  unmatchCaloJets->Branch("muMax",&muMax,"muMax/F");      unmatchCaloJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  unmatchCaloJets->Branch("eMax",&eMax,"eMax/F");         unmatchCaloJets->Branch("jet60",&jet60,"jet60/I");
  unmatchCaloJets->Branch("chSum",&chSum,"chSum/F");      unmatchCaloJets->Branch("jet60_prescl",&jet60_prescl,"jet60_prescl/I");
  unmatchCaloJets->Branch("jet40",&jet40,"jet40/I");      unmatchCaloJets->Branch("jet40_prescl",&jet40_prescl,"jet40_prescl/I");
  unmatchCaloJets->Branch("hcalSum",&hcalSum,"hcalSum/F");unmatchCaloJets->Branch("ecalSum",&ecalSum,"ecalSum/F");
  unmatchCaloJets->Branch("run_value",&run_value,"run_value/I");
  unmatchCaloJets->Branch("evt_value",&evt_value,"evt_value/I");
  unmatchCaloJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  unmatchCaloJets->Branch("caloeta", &caloeta, "caloeta/F");
  unmatchCaloJets->Branch("calophi", &calophi, "calophi/F");
  if(dataset=="MC") unmatchCaloJets->Branch("subid",&subid,"subid/I"); 
  if(dataset=="MC") unmatchCaloJets->Branch("calorefpt",&calorefpt,"calorefpt/F");
  if(dataset=="MC") unmatchCaloJets->Branch("pthat",&pthat,"pthat/F"); if(dataset=="MC") unmatchCaloJets->Branch("vz",&vz,"vz/F");
  if(dataset=="MC") unmatchCaloJets->Branch("weight",&weight,"weight/F");
  unmatchCaloJets->Branch("calorawpt",&calorawpt,"calorawpt/F");
  unmatchCaloJets->Branch("calojtpu",&calojtpu,"calojtpu/F");


  TH1F * hEvents = new TH1F("hEvents","",20,0,10);

  // define the histograms 

  // define the histograms necessary for that MC closure test.
  TH2F *hpp_mcclosure_matrix_HLT;
  TH2F *hpp_mcclosure_matrix;
  //TH2F *hpp_response;
  TH1F *hpp_mcclosure_JetComb_data;
  TH1F *hpp_mcclosure_data;
  TH1F *hpp_mcclosure_Jet80_data;
  TH1F *hpp_mcclosure_Jet60_data;
  TH1F *hpp_mcclosure_Jet40_data;
  TH1F *hpp_mcclosure_gen;
  TH1F *hpp_mcclosure_JetComb_gen;
  TH1F *hpp_mcclosure_Jet80_gen;
  TH1F *hpp_mcclosure_Jet60_gen;
  TH1F *hpp_mcclosure_Jet40_gen;
  
  hpp_mcclosure_matrix_HLT = new TH2F(Form("hpp_mcclosure_matrix_HLT_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  hpp_mcclosure_matrix = new TH2F(Form("hpp_mcclosure_matrix_R%d_%s",radius,etaWidth),Form("Matrix for mcclosure refpt jtpt from Jet triggers R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  //cout<<"C"<<endl;
  hpp_mcclosure_data = new TH1F(Form("hpp_mcclosure_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_JetComb_data = new TH1F(Form("hpp_mcclosure_JetComb_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger combined  R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_Jet80_data = new TH1F(Form("hpp_mcclosure_Jet80_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 80  R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_Jet60_data = new TH1F(Form("hpp_mcclosure_Jet60_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 60  R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_Jet40_data = new TH1F(Form("hpp_mcclosure_Jet40_data_R%d_%s",radius,etaWidth),Form("data for unfolding mc closure test trigger 40  R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);

  hpp_mcclosure_gen = new TH1F(Form("hpp_mcclosure_gen_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_JetComb_gen = new TH1F(Form("hpp_mcclosure_gen_JetComb_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger combined R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_Jet80_gen = new TH1F(Form("hpp_mcclosure_gen_Jet80_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 80 R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_Jet60_gen = new TH1F(Form("hpp_mcclosure_gen_Jet60_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 60 R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);
  hpp_mcclosure_Jet40_gen = new TH1F(Form("hpp_mcclosure_gen_Jet40_R%d_%s",radius,etaWidth),Form("gen spectra for unfolding mc closure test trigger 40 R%d %s ",radius,etaWidth),nbins_pt,boundaries_pt);

  
  // declare the 2d calo and pf candidate vectors, deltaR_calovsPF is going to be a 3d vector like so:
  // calo jet on x axis, pf jet on y axis, and delta R, calopT, pfpT on z axis.
  // we also need to add in the trigger information here since we need to know what trigger the matching comes from. maybe i can run that later as an added check. for now just to see if the algorithm is working should try to run things. 
  //vector<vector<double> > caloJet;
  //vector<vector<double> > pfJet;

  // start the loop process

  for(int k = 0;k<no_radius;k++){

    Long64_t nentries = jetpp1[2][k]->GetEntries();
    //Long64_t nentries = 100;

    Long64_t fentries= 1;
    if(dataset=="MC"){
      TEventList *el = new TEventList("el","el");
      double pthat_upper = boundaries_pthat[startfile + 1];
      stringstream selection; selection<<"pthat<"<<pthat_upper;
      
      jetpp1[2][k]->Draw(">>el",selection.str().c_str());
      fentries = el->GetN();
      delete el;
    }
    
    for(Long64_t nentry = 0; nentry<nentries;nentry++){
      if(printDebug)cout<<"event no = "<<nentry<<endl;
      for(int t = 0;t<N;t++)  jetpp1[t][k]->GetEntry(nentry);
      
      hEvents->Fill(0);

      //int centBin = findBin(hiBin_1);
      //if(centBin==-1) continue;
      //if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 
      if(pPAcollisionEventSelectionPA_1==0) continue; 
      hEvents->Fill(1);
      if(dataset=="Data") if(pHBHENoiseFilter_1==0 ) continue;
      hEvents->Fill(2);

      if(fabs(vz_1)>15) continue;
      //cout<<"passed the selection"<<endl;

      hEvents->Fill(3);

      weight = 1;

      if(dataset=="MC"){
	
	int pthatBin = hpthatBin->FindBin(pthat_1);
	
	double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/fentries;

	Float_t weight_vz = fVzPP->Eval(vz_1);
	weight = scale*weight_vz;

      }

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

	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(pfjet_eta); // 35 - pf eta
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(pfjet_phi); // 36 - pf phi
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(calojet_eta); // 37 - calo eta
	  deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(calojet_phi); // 38 - calo phi
	  
	  if(dataset=="MC") {
	    deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(refpt_1[g]); // 39 calo ref pt 
	    deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(refpt_2[j]); // 40 pf ref pt 
	    deltaR_calovsPF[calomatchcounter][pfmatchcounter].push_back(subid_2[j]); // 41 - subid only for MC
	    
	  }
	  
	  
	  ++pfmatchcounter;

	}// pf jet loop
	++calomatchcounter;
	
      }// calo jet loop

      // set the event variables before proceeding to the matching. 
      jet80 = jet80_1;
      jet80_prescl = jet80_p_1;
      jet60 = jet60_1;
      jet60_prescl = jet60_p_1;
      jet40 = jet40_1;
      jet40_prescl = jet40_p_1;
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
	  pfrefpt = deltaR_calovsPF[small_calo][small_pf][40];
	  calorefpt = deltaR_calovsPF[small_calo][small_pf][39];
	  subid = deltaR_calovsPF[small_calo][small_pf][41]; 
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
	pfeta = deltaR_calovsPF[small_calo][small_pf][35];
	pfphi = deltaR_calovsPF[small_calo][small_pf][36];
	caloeta = deltaR_calovsPF[small_calo][small_pf][37];
	calophi = deltaR_calovsPF[small_calo][small_pf][38];
	hcalSum = deltaR_calovsPF[small_calo][small_pf][17];
	ecalSum = deltaR_calovsPF[small_calo][small_pf][18];
	
	//matchJets->Fill();
	
	if(c==0) hEvents->Fill(4);
	  
	if(dataset=="MC"){
	  // fill in the matched histograms with the Jet ID cuts:
	  if(subid!=0) continue;
	
	  Float_t Sumcand = chSum + phSum + neSum + muSum;

	  //if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < ((Float_t)18/7 *(Float_t)calopt/pfpt - (Float_t)9/7)){
	    if(nentry%2==1) {
	      hpp_mcclosure_matrix->Fill(pfrefpt, pfpt, weight);
	      hpp_mcclosure_gen->Fill(pfrefpt, weight);
	    }
	    if(nentry%2==1) {
	      hpp_mcclosure_data->Fill(pfpt, weight);
	    }
	    //}
	  
	  if(jet40 == 1 && jet60==0 && jet80 == 0){

	    //if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < ((Float_t)18/7 *(Float_t)calopt/pfpt - (Float_t)9/7)){
	      if(nentry%2==1) {
		hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
		hpp_mcclosure_Jet40_gen->Fill(pfrefpt, weight);
	      }
	      if(nentry%2==1) {
		hpp_mcclosure_Jet40_data->Fill(pfpt, weight);
	      }
	      //}

	    // if(calopt/pfpt > 0.85) {
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
	    // 	hpp_mcclosure_Jet40_gen->Fill(pfrefpt, weight);
	    //   }
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_Jet40_data->Fill(pfpt, weight);
	    //   }
	    // }

	    // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) {
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
	    // 	hpp_mcclosure_Jet40_gen->Fill(pfrefpt, weight);
	    //   }
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_Jet40_data->Fill(pfpt, weight);
	    //   }
	    // }
	  
	  }

	  if(jet60 == 1 && jet80 == 0){

	    //if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < ((Float_t)18/7 *(Float_t)calopt/pfpt - (Float_t)9/7)){
	      if(nentry%2==1) {
		hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
		hpp_mcclosure_Jet60_gen->Fill(pfrefpt, weight);
	      }
	      if(nentry%2==1) {
		hpp_mcclosure_Jet60_data->Fill(pfpt, weight);
	      }
	      //}

	    // if(calopt/pfpt > 0.85) {
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
	    // 	hpp_mcclosure_Jet60_gen->Fill(pfrefpt, weight);
	    //   }
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_Jet60_data->Fill(pfpt, weight);
	    //   }
	    // }

	    // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) {
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
	    // 	hpp_mcclosure_Jet60_gen->Fill(pfrefpt, weight);
	    //   }
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_Jet60_data->Fill(pfpt, weight);
	    //   }
	    // }
	  
	  
	  }

	  if(jet80==1){

	    //if(calopt/pfpt > 0.5 && calopt/pfpt <= 0.85 && eMax/Sumcand < ((Float_t)18/7 *(Float_t)calopt/pfpt - (Float_t)9/7)){
	      if(nentry%2==1) {
		hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
		hpp_mcclosure_Jet80_gen->Fill(pfrefpt, weight);
	      }
	      if(nentry%2==1) {
		hpp_mcclosure_Jet80_data->Fill(pfpt, weight);
	      }
	      //}

	    // if(calopt/pfpt > 0.85) {
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
	    // 	hpp_mcclosure_Jet80_gen->Fill(pfrefpt, weight);
	    //   }
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_Jet80_data->Fill(pfpt, weight);
	    //   }
	    // }

	    // if(calopt/pfpt <= 0.5 && eMax/Sumcand < 0.05) {
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
	    // 	hpp_mcclosure_Jet80_gen->Fill(pfrefpt, weight);
	    //   }
	    //   if(nentry%2==1) {
	    // 	hpp_mcclosure_Jet80_data->Fill(pfpt, weight);
	    //   }
	    // }
	  
	  
	  }

	}


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
	  pfrefpt = deltaR_calovsPF[0][b][40];
	  subid = deltaR_calovsPF[0][b][41]; 
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
	pfeta = deltaR_calovsPF[0][b][35];
	pfphi = deltaR_calovsPF[0][b][36];
	hcalSum = deltaR_calovsPF[0][b][17];
	ecalSum = deltaR_calovsPF[0][b][18];
	
	//unmatchPFJets->Fill();

	if(b==0) hEvents->Fill(5);

	if(dataset=="MC"){
	  // fill in the matched histograms with the Jet ID cuts:
	  if(subid!=0) continue;
	
	  Float_t Sumcand = chSum + phSum + neSum + muSum;

	  //if(eMax/Sumcand < 0.05) {
	    if(nentry%2==1) {
	      hpp_mcclosure_matrix->Fill(pfrefpt, pfpt, weight);
	      hpp_mcclosure_gen->Fill(pfrefpt, weight);
	    }
	    if(nentry%2==1) {
	      hpp_mcclosure_data->Fill(pfpt, weight);
	    }
	    //}
	  
	  
	  if(jet40 == 1 && jet60 == 0 && jet80 == 0){

	    //if(eMax/Sumcand < 0.05) {
	      if(nentry%2==1) {
		hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
		hpp_mcclosure_Jet40_gen->Fill(pfrefpt, weight);
	      }
	      if(nentry%2==1) {
		hpp_mcclosure_Jet40_data->Fill(pfpt, weight);
	      }
	      //}
	  
	  }

	  if(jet60 == 1 && jet80 == 0){

	    //if(eMax/Sumcand < 0.05) {
	      if(nentry%2==1) {
		hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
		hpp_mcclosure_Jet60_gen->Fill(pfrefpt, weight);
	      }
	      if(nentry%2==1) {
		hpp_mcclosure_Jet60_data->Fill(pfpt, weight);
	      }
	      //}
	  
	  
	  }

	  if(jet80==1){

	    //if(eMax/Sumcand < 0.05) {
	      if(nentry%2==1) {
		hpp_mcclosure_matrix_HLT->Fill(pfrefpt, pfpt, weight);
		hpp_mcclosure_Jet80_gen->Fill(pfrefpt, weight);
	      }
	      if(nentry%2==1) {
		hpp_mcclosure_Jet80_data->Fill(pfpt, weight);
	      }
	      //}
	  
	  
	  }
	}

	
      }// unmatched PF jets
      
      if(printDebug)cout<<"now going to find unmatched calo jets"<<endl;

      if(deltaR_calovsPF[0].size() == 0) continue;
      
      for(int a = 0;a<deltaR_calovsPF.size();++a){
	if(printDebug)cout<<"calo jet iteration "<<a<<endl;

	if(deltaR_calovsPF[a][0][16] == 1) continue;

	calopt = deltaR_calovsPF[a][0][3];
	deltar = deltaR_calovsPF[a][0][0];
	if(dataset=="MC"){
	  calorefpt = deltaR_calovsPF[a][0][39];
	  subid = deltaR_calovsPF[a][0][41]; 
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
	caloeta = deltaR_calovsPF[a][0][37];
	calophi = deltaR_calovsPF[a][0][38];
	hcalSum = deltaR_calovsPF[a][0][31];
	ecalSum = deltaR_calovsPF[a][0][32];
	
	//unmatchCaloJets->Fill();

      }

    }// event loop 

  }// radius loop


  // matchJets->Write();
  // matchJets->Print();
  // unmatchPFJets->Write();
  // unmatchPFJets->Print();
  // unmatchCaloJets->Write();
  // unmatchCaloJets->Print();


  hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet80_data);
  hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet60_data);
  hpp_mcclosure_JetComb_data->Add(hpp_mcclosure_Jet40_data);

  hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet80_gen);
  hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet60_gen);
  hpp_mcclosure_JetComb_gen->Add(hpp_mcclosure_Jet40_gen);

  hpp_mcclosure_matrix_HLT->Write();
  hpp_mcclosure_matrix->Write();
  hpp_mcclosure_JetComb_data->Write();
  hpp_mcclosure_data->Write();
  hpp_mcclosure_Jet80_data->Write();
  hpp_mcclosure_Jet60_data->Write();
  hpp_mcclosure_Jet40_data->Write();
  hpp_mcclosure_JetComb_gen->Write();    
  hpp_mcclosure_gen->Write();    
  hpp_mcclosure_Jet80_gen->Write();
  hpp_mcclosure_Jet60_gen->Write();
  hpp_mcclosure_Jet40_gen->Write();
    

  hEvents->Write();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}
