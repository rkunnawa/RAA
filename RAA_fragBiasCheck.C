// Raghav Kunnawalkam Elayavalli
// Feb 17th 2015
// Rutgers

//
// Macro to check the fragmentation lower pT bias similar to fig.4 done by Alice. So its the ratio of the no of jets with a leading track pT cut / no track pt cut. 
//
//

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

static const int nbins_cent = 6;
static Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 5 to get your actual centrality
static Double_t ncoll[nbins_cent+1] = { 1660, 1310, 745, 251, 62.8, 10.8 ,362.24}; //last one is for 0-200 bin. 

static const int no_radius = 3;//necessary for the RAA analysis  
static const int list_radius[no_radius] = {2,3,4};

static const Float_t effecPrescl = 2.047507;

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

int findBin(int hiBin){
  int binNo = 0;

  for(int i = 0;i<nbins_cent;i++){
    if(hiBin>=5*boundaries_cent[i] && hiBin<5*boundaries_cent[i+1]) {
      binNo = i;
      break;
    }
  }

  return binNo;
}


using namespace std;

void RAA_fragBiasCheck(int startfile = 0, int endfile = 1, char *algo="Pu", char *jet_type="PF"){

  
  TH1::SetDefaultSumw2();
  //gStyle->SetOptStat(0);
  
  TStopwatch timer;
  timer.Start();
  
  TDatime date;
  
  cout<<"Running for Algo = "<<algo<<" "<<jet_type<<endl;
  
  bool printDebug = false;

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
  
  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }
  
  const int N = 7;
  
  TChain *jetpbpb1[N][no_radius];
    
  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%d%sJetAnalyzer",algo,list_radius[k],jet_type);
    dir[3][k] = "hiEvtAnalyzer";
    dir[4][k] = "hltobject";
    dir[5][k] = "pfcandAnalyzer";
    dir[6][k] = "pfTowers";
  }

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree",
    "jetObjTree",
    "pfTree",
    "tower"
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
    jetpbpb1[2][k]->AddFriend(jetpbpb1[3][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[4][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[5][k]);
    jetpbpb1[2][k]->AddFriend(jetpbpb1[6][k]);
  }// radius loop ends

  cout<<"total no of entries in the combined forest files = "<<jetpbpb1[2][0]->GetEntries()<<endl;

  //file 1: 
  // jet tree 1
  int nrefe_1;
  float pt_1[1000];
  //float old_pt3[1000];
  float raw_1[1000];
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
  // trigger tree
  int L1_sj36_1;
  int L1_sj52_1;
  int jet55_1;
  int jet65_1;
  int jet80_1;
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

  Int_t nTower;
  Float_t eTower[NOBJECT_MAX];
  Float_t etTower[NOBJECT_MAX];
  Float_t etaTower[NOBJECT_MAX];
  Float_t phiTower[NOBJECT_MAX];
  
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

    jetpbpb1[5][k]->SetBranchAddress("nPFpart", &nPFpart);
    jetpbpb1[5][k]->SetBranchAddress("pfId", pfId);
    jetpbpb1[5][k]->SetBranchAddress("pfPt", pfPt);
    jetpbpb1[5][k]->SetBranchAddress("pfVsPtInitial", pfVsPtInitial);
    jetpbpb1[5][k]->SetBranchAddress("pfVsPt", pfVsPt);
    jetpbpb1[5][k]->SetBranchAddress("pfEta", pfEta);
    jetpbpb1[5][k]->SetBranchAddress("pfPhi", pfPhi);
    jetpbpb1[5][k]->SetBranchAddress("pfArea", pfArea);
    jetpbpb1[5][k]->SetBranchAddress("vn",&v_n);
    jetpbpb1[5][k]->SetBranchAddress("psin",&psi_n);
    jetpbpb1[5][k]->SetBranchAddress("sumpt",&sumpT);

    jetpbpb1[6][k]->SetBranchAddress("n",&nTower);
    jetpbpb1[6][k]->SetBranchAddress("e",&eTower);
    jetpbpb1[6][k]->SetBranchAddress("et",&etTower);
    jetpbpb1[6][k]->SetBranchAddress("eta",&etaTower);
    jetpbpb1[6][k]->SetBranchAddress("phi",&phiTower);

  }//radius loop

  // Declare the output file:
  TFile fout(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_fragBiasCheck_ak%s%s_%d_%d.root",algo,jet_type,date.GetDate(),endfile),"RECREATE");
  fout.cd();
  
  TH1F *hJets_noTrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_3TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_5TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_7TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];
  TH1F *hJets_10TrkpTCut[no_radius][nbins_eta][nbins_cent][trigValue];

  for(int k = 0;k<no_radius;k++){
    for(int l = 0;l<trigValue;l++){
      for(int j = 0;j<nbins_eta;j++){
	for(int i = 0;i<nbins_cent;i++){
	  hJets_noTrkpTCut[k][j][i][l] = new TH1F(Form("hJets_noTrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with no fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_3TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_3TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 3GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_5TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_5TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 5GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_7TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_7TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 7GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  hJets_10TrkpTCut[k][j][i][l] = new TH1F(Form("hJets_10TrkpTCut_%s_R%d_%s_cent%d",trigName[l],list_radius[k],etaWidth[j],i),Form("Jet Spectra with 10GeV fragmentation cut %s R%d %s %2.0f - %2.0f cent",trigName[l],list_radius[k],etaWidth[j],5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
	  
	}// cent loop
	
      }// eta loop

    }// trigger loop

    cout<<"Running Data for R = "<<list_radius[k]<<endl;

    Long64_t nentries_jet55or65 = jetpbpb1[2][k]->GetEntries();
    //if(printDebug)cout<<"nentries_jet55or65or80 = "<<nentries_jet55or65<<endl;
    cout<<"nentries_jet55or65or80 = "<<nentries_jet55or65<<endl;
    if(printDebug) nentries_jet55or65 = 2;

    for(int jentry = 0;jentry<nentries_jet55or65;jentry++){
      jetpbpb1[0][k]->GetEntry(jentry);
      jetpbpb1[1][k]->GetEntry(jentry);
      jetpbpb1[2][k]->GetEntry(jentry);    
      jetpbpb1[3][k]->GetEntry(jentry);
      jetpbpb1[4][k]->GetEntry(jentry);
      jetpbpb1[5][k]->GetEntry(jentry);

      //if(printDebug && jentry%100000==0)cout<<"Jet 55or65 file"<<endl;
      //if(printDebug)cout<<jentry<<": event = "<<evt_1<<"; run = "<<run_1<<"; lumi = "<<lumi_1<<endl;
      int centBin = findBin(hiBin_1);//tells us the centrality of the event.

      if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 
      if(fabs(vz_1)>15) continue;
int jetCounter = 0;//counts jets which are going to be used in the supernova cut rejection. 
      
     
      for(int j = 0;j<nbins_eta;j++){
	  
	for(int g = 0;g<nrefe_1;g++){
	    
	  if(eta_1[g]>=boundaries_eta[j][0] && eta_1[g]<boundaries_eta[j][1]){
	      
	    if(pt_1[g]>=50) jetCounter++;
	    
	  }//eta selection cut
	  
	}// jet loop
	
      }//eta bins loop
      
      // apply the correct supernova selection cut rejection here: 
      if(hiNpix_1 > 38000 - 500*jetCounter){
       	if(printDebug) cout<<"removed this supernova event"<<endl;
      	continue;
      }

      for(int j = 0;j<nbins_eta;j++){

	for(int g = 0;g<nrefe_1;g++){ // this is the loop for the  Jets we are interested in.  

	  vector <Float_t> inJetPFcand_pT;
	  
	  if(eta_1[g]<boundaries_eta[j][0] || eta_1[g]>=boundaries_eta[j][1]) continue;

	  //if(chMax_1[g]/pt_1[g]<0.05) continue; // jet ID cut, but the fragmentation would probably get rid of this anyway. 

	  // lets start the fragmentation check:

	  
	  // have to run through all the particle flow candidates:
	  for(int ipf = 0; ipf<nPFpart; ipf++){

	    if(TMath::Sqrt((eta_1[g] - pfEta[ipf])*(eta_1[g] - pfEta[ipf]) + (phi_1[g] - pfPhi[ipf])*(phi_1[g] - pfPhi[ipf])) <= (float) list_radius[k]/10){

	      inJetPFcand_pT.push_back(pfPt[ipf]);
	      
	    }// checking for candidates to be inside the jet cone 

	  } // pf flow candidate loop

	  // Now find the pT of the largest candidate inside the jet.
	  float large = inJetPFcand_pT[0];
	  for(int a = 0;a<inJetPFcand_pT.size()-1;a++){
	    if(large < inJetPFcand_pT[a]) large = inJetPFcand_pT[a];	      
	  }
	  	  
	  /*
	  // have to run through all the particle flow candidates:
	  for(int ipf = 0; ipf<nTower; ipf++){

	    if(TMath::Sqrt((eta_1[g] - etaTower[ipf])*(eta_1[g] - etaTower[ipf]) + (phi_1[g] - phiTower[ipf])*(phi_1[g] - phiTower[ipf])) <= (float) list_radius[k]/10){

	      inJetPFcand_pT.push_back(etTower[ipf]);
	      
	    }// checking for candidates to be inside the jet cone 

	  } // pf flow candidate loop

	  // Now find the pT of the smallest candidate inside the jet.
	  float large = inJetPFcand_pT[0];
	  for(int a = 1;a<inJetPFcand_pT.size();a++){
	    if(large < inJetPFcand_pT[a])
	      large = inJetPFcand_pT[a];
	  }
	  */
	  // test checking the sum of the pf candidate to equal the jet pt
	  if(printDebug) cout<<"jtpt = "<<pt_1[g]<<"large candidate = "<<large<<endl;
	  
	  // now that we have the smallest candidate we can check for the different fragmentation biases:
	 	  
	  if(jet80_1==1 && L1_sj52_1==1) {

	    hJets_noTrkpTCut[k][j][centBin][2]->Fill(pt_1[g]);
	    if(large > 3) hJets_3TrkpTCut[k][j][centBin][2]->Fill(pt_1[g]);
	    if(large > 5) hJets_5TrkpTCut[k][j][centBin][2]->Fill(pt_1[g]);
	    if(large > 7) hJets_7TrkpTCut[k][j][centBin][2]->Fill(pt_1[g]);
	    if(large > 10) hJets_10TrkpTCut[k][j][centBin][2]->Fill(pt_1[g]);
	      
	  }

	  if(jet65_1==1 && L1_sj36_1==1){

	    hJets_noTrkpTCut[k][j][centBin][1]->Fill(pt_1[g]);
	    if(large > 3) hJets_3TrkpTCut[k][j][centBin][1]->Fill(pt_1[g]);
	    if(large > 5) hJets_5TrkpTCut[k][j][centBin][1]->Fill(pt_1[g]);
	    if(large > 7) hJets_7TrkpTCut[k][j][centBin][1]->Fill(pt_1[g]);
	    if(large > 10) hJets_10TrkpTCut[k][j][centBin][1]->Fill(pt_1[g]);

	  }

	  if(jet55_1==1 && jet65_1==0 && jet80_1==0 && L1_sj36_1==1){
	    
	    hJets_noTrkpTCut[k][j][centBin][0]->Fill(pt_1[g],effecPrescl);
	    if(large > 3) hJets_3TrkpTCut[k][j][centBin][0]->Fill(pt_1[g],effecPrescl);
	    if(large > 5) hJets_5TrkpTCut[k][j][centBin][0]->Fill(pt_1[g],effecPrescl);
	    if(large > 7) hJets_7TrkpTCut[k][j][centBin][0]->Fill(pt_1[g],effecPrescl);
	    if(large > 10) hJets_10TrkpTCut[k][j][centBin][0]->Fill(pt_1[g],effecPrescl);

	  }
	
	  for(int t = 0;t<trigValue-1;t++) {

	    hJets_noTrkpTCut[k][j][centBin][trigValue-1]->Add(hJets_noTrkpTCut[k][j][centBin][t]);
	    hJets_3TrkpTCut[k][j][centBin][trigValue-1]->Add(hJets_3TrkpTCut[k][j][centBin][t]);
	    hJets_5TrkpTCut[k][j][centBin][trigValue-1]->Add(hJets_5TrkpTCut[k][j][centBin][t]);
	    hJets_7TrkpTCut[k][j][centBin][trigValue-1]->Add(hJets_7TrkpTCut[k][j][centBin][t]);
	    hJets_10TrkpTCut[k][j][centBin][trigValue-1]->Add(hJets_10TrkpTCut[k][j][centBin][t]);
	    
	  }

	}// jet loop

      }// eta loop

    }// event loop

  }// radius loop

  
  fout.Write();
  fout.Close();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;

  
}// macro end 
