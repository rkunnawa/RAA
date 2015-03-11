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
  int binNo = 0;

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

using namespace std;

void RAA_calo_pf_JetCorrelation(int startfile = 0, int endfile = 1, int radius=3, char *algo = "Pu", int deltaR=2/*which i will divide by 10 later when using*/){

  TH1::SetDefaultSumw2();

  TStopwatch timer;
  timer.Start();

  TDatime date;

  cout<<"Running for Algo = "<<algo<<" "<<endl;
  
  bool printDebug = false;

  // Since this has to run on the HiForest files, it is best to run them as condor jobs similar to the PbPb data read macro, along with the jet trees which get their branch address set. 
  
  std::string infile1;
  infile1 = "jetRAA_PbPb_data_forest.txt";
  
  std::ifstream instr1(infile1.c_str(),std::ifstream::in);
  std::string filename1;

  cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr1>>filename1;
  }
  
  const int N = 6;
  
  TChain *jetpbpb1[N][no_radius];

  string dir[N][no_radius];
  
  for(int k = 0;k<no_radius;k++){
    dir[0][k] = "hltanalysis";
    dir[1][k] = "skimanalysis";
    dir[2][k] = Form("ak%s%dCaloJetAnalyzer",algo,list_radius[k]);
    dir[3][k] = Form("ak%s%dPFJetAnalyzer",algo,list_radius[k]);
    dir[4][k] = "hiEvtAnalyzer";
    dir[5][k] = "hltobject";
    //dir[6][k] = "pfcandAnalyzer";
  }
  
  
  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "t",
    "HiTree",
    "jetObjTree",
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
    jetpbpb1[2][k]->AddFriend(jetpbpb1[5][k]);
    
    jetpbpb1[3][k]->AddFriend(jetpbpb1[0][k]);
    jetpbpb1[3][k]->AddFriend(jetpbpb1[1][k]);
    jetpbpb1[3][k]->AddFriend(jetpbpb1[4][k]);
    jetpbpb1[3][k]->AddFriend(jetpbpb1[5][k]);
    
  }// radius loop ends
 
  // Ok this should now work. 

  cout<<"total no of entries in the combined forest files = "<<jetpbpb1[2][0]->GetEntries()<<endl;
  //  if(printDebug)cout<<"total no of entries in the Jet80 Tree     = "<<jetpbpb2[2][0]->GetEntries()<<endl;

  //set the branch addresses:
  // jet tree 1 - Calo 
  int nrefe_1;
  float pt_1[1000];
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
  
  // jet tree 2 - PF
  int nrefe_2;
  float pt_2[1000];
  float raw_2[1000];
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
    
    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
    //jetpbpb1[2][k]->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
    //jetpbpb1[2][k]->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
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

  // declare the histograms:
  TH1F *hCaloPFCorr[TrigValue][nbins_cent+1], *hCalo[TrigValue][nbins_cent+1], *hPF[TrigValue][nbins_cent+1], *hRatio[TrigValue][nbins_cent+1];
  TH2F * hCaloPFpt[TrigValue][nbins_cent+1];
  TH2F * hCaloPFCorr_pt[TrigValue][nbins_cent+1];
  TH2F * hDeltaR_deltapT[TrigValue][nbins_cent+1];
  // TH1F * hDeltaR[TrigValue][nbins_cent+1];
  // TH1F * hDeltapT[TrigValue][nbins_cent+1];

  for(int a = 0;a<TrigValue;a++){
    for(int i = 0;i<nbins_cent+1;i++){
      hCaloPFCorr[a][i] = new TH1F(Form("hCaloPFCorr_%s_cent%d",TrigName[a],i),Form("Ratio Calo jet pT to correlated PF jet for %s %s",TrigName[a],centWidth[i]),200,0,10);
      hCaloPFCorr[a][i]->SetXTitle("PF Jet pT / Calo Jet PT");
      hCaloPFCorr[a][i]->SetYTitle("Counts");
      
      hCaloPFpt[a][i] = new TH2F(Form("hCaloPFpt_%s_cent%d",TrigName[a],i), Form("Matched Calo jet pT vs PF Jet pT %s %s", TrigName[a], centWidth[i]), nbins_pt,boundaries_pt, nbins_pt, boundaries_pt);
      hCaloPFpt[a][i]->SetXTitle("Calo jet pT");
      hCaloPFpt[a][i]->SetYTitle("Matched PF jet pT");

      hCaloPFCorr_pt[a][i] = new TH2F(Form("hCaloPFCorr_pt_%s_cent%d",TrigName[a],i), Form(" Ratio of PF/Calo vs calo pT %s %s", TrigName[a], centWidth[i]), nbins_pt,boundaries_pt, 200,0,10);
      hCaloPFCorr_pt[a][i]->SetXTitle("Calo jet pT");
      hCaloPFCorr_pt[a][i]->SetYTitle("PF Jet pT / Calo Jet PT");

      hCalo[a][i] = new TH1F(Form("hCalo_%s_cent%d",TrigName[a],i),Form("Calo jet spectra %s %s",TrigName[a],centWidth[i]),nbins_pt,boundaries_pt);
      hCalo[a][i]->SetXTitle("Jet p_{T} (GeV/c)");
      hCalo[a][i]->SetYTitle("Counts");

      hPF[a][i] = new TH1F(Form("hPF_%s_cent%d",TrigName[a],i),Form("PF jet spectra %s %s",TrigName[a],centWidth[i]),nbins_pt,boundaries_pt);
      hPF[a][i]->SetXTitle("Jet p_{T} (GeV/c)");
      hPF[a][i]->SetYTitle("Counts");

      hDeltaR_deltapT[a][i] = new TH2F(Form("hDeltaR_deltapT_%s_cent%d",TrigName[a],i), Form(" delta R for matched jets vs delta pT %s %s", TrigName[a], centWidth[i]), 200,0,1, 200,0,200);
      // hDeltaR[a][i] = new TH1F(Form("hPF_%s_cent%d",TrigName[a],i),Form("PF jet spectra %s %s",TrigName[a],centWidth[i]),nbins_pt,boundaries_pt);

    }
  }

  cout<<"after histogram declaration"<<endl;
  TNtuple * matchJets = new TNtuple("matchJets","Ntuple containing important information about matched jets","calopt:pfpt:deltaR:chMax:phMax:neMax:muMax:eMax:chSum:phSum:neSum:muSum:eSum");
  
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
    //nentries = 10;
    
    for(Long64_t nentry = 0; nentry<nentries;nentry++){
      //cout<<"event no = "<<nentry<<endl;
      for(int t = 0;t<N;t++)  jetpbpb1[t][k]->GetEntry(nentry);
      
      int centBin = findBin(hiBin_1);

      if(pHBHENoiseFilter_1==0 || pcollisionEventSelection_1==0) continue; 

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

      for(int g = 0;g<nrefe_1;g++){
	
	calojet_eta = eta_1[g];
	calojet_phi = phi_1[g];
	calojet_pt = pt_1[g];
	
	if(calojet_pt < 10) continue;
	
	int pfmatchcounter = 0;
	deltaR_calovsPF.push_back(vector<vector<double> > ());
	
	calosize++;
	
	for(int j = 0;j<nrefe_2;j++){

	  //if(chMax_2[j]/pt_2[j]<0.02 || eMax_2[j]/pt_2[j]>0.6) continue;

	  pfjet_eta = eta_2[j];
	  pfjet_phi = phi_2[j];
	  pfjet_pt = pt_2[j];

	  if(pfjet_pt < 10) continue;

	  deltaRCaloPF = (Float_t)TMath::Sqrt((calojet_eta - pfjet_eta)*(calojet_eta - pfjet_eta) + (calojet_phi - pfjet_phi)*(calojet_phi - pfjet_phi));
	  // if(deltaRCaloPF > (Float_t)deltaR/10) continue;
	  deltapT = TMath::Abs(calojet_pt - pfjet_pt);

	  deltaR_calovsPF[g].push_back(vector<double> ());
	  deltaR_calovsPF[g][pfmatchcounter].push_back(deltaRCaloPF); // 0 - delta R
	  deltaR_calovsPF[g][pfmatchcounter].push_back(deltapT); // 1 - delta pT	  
	  deltaR_calovsPF[g][pfmatchcounter].push_back(g); // 2 - calo counter 
	  deltaR_calovsPF[g][pfmatchcounter].push_back(calojet_pt); // 3 - calo jet pT 
	  deltaR_calovsPF[g][pfmatchcounter].push_back(j); // 4 - pf counter 
	  deltaR_calovsPF[g][pfmatchcounter].push_back(pfjet_pt); // 5 - pf jet pT
	  // this will have the candidate variables for the matched jets. 
	  deltaR_calovsPF[g][pfmatchcounter].push_back(chMax_2[j]); // 6 - chMax 
	  deltaR_calovsPF[g][pfmatchcounter].push_back(phMax_2[j]); // 7 - phMax
	  deltaR_calovsPF[g][pfmatchcounter].push_back(neMax_2[j]); // 8 - neMax
	  deltaR_calovsPF[g][pfmatchcounter].push_back(muMax_2[j]); // 9 - muMax
	  deltaR_calovsPF[g][pfmatchcounter].push_back(eMax_2[j]); // 10 - eMax
	  deltaR_calovsPF[g][pfmatchcounter].push_back(chSum_2[j]); // 11 - chSum
	  deltaR_calovsPF[g][pfmatchcounter].push_back(phSum_2[j]); // 12 - phSum
	  deltaR_calovsPF[g][pfmatchcounter].push_back(neSum_2[j]); // 13 - neSum
	  deltaR_calovsPF[g][pfmatchcounter].push_back(muSum_2[j]); // 14 - muSum
	  deltaR_calovsPF[g][pfmatchcounter].push_back(eSum_2[j]); // 15 - eSum


	  ++pfmatchcounter;

	}// pf jet loop
	
      }// calo jet loop

      // now that we have the 2D matrix, lets find the smallest delta R element from that and fill in the value of the 

      Float_t smallDeltaR = 10;
      Int_t small_calo = 0;
      Int_t small_pf = 0;

      for(int c = 0;c<deltaR_calovsPF.size();++c){
	//cout<<"going through all rows in the  matrix "<<c<<endl;
	for(int a = 0;a<deltaR_calovsPF.size();++a){
	  //cout<<"calo jet iteration "<<a<<endl;
	  for(int b = 0;b<deltaR_calovsPF[a].size();++b){
	    //cout<<"pf jet iteration "<<b<<endl;

	    //smallDeltaR = 10;
	    if(deltaR_calovsPF[a][b][0] == 100) break;
	    if(smallDeltaR > deltaR_calovsPF[a][b][0]){
	     
	      smallDeltaR = deltaR_calovsPF[a][b][0];
	      small_calo = a;
	      small_pf = b;
	      
	    }
	  }
	  //cout<<endl;
	}
	
	if(smallDeltaR > (Float_t)deltaR/10) continue;
	
	calojet_pt = deltaR_calovsPF[small_calo][small_pf][3];
	pfjet_pt = deltaR_calovsPF[small_calo][small_pf][5];
	deltapT = deltaR_calovsPF[small_calo][small_pf][1];
	deltaRCaloPF = deltaR_calovsPF[small_calo][small_pf][0];

	matchJets->Fill(calojet_pt,pfjet_pt,deltaRCaloPF,deltaR_calovsPF[small_calo][small_pf][6],
			deltaR_calovsPF[small_calo][small_pf][7],deltaR_calovsPF[small_calo][small_pf][8],
			deltaR_calovsPF[small_calo][small_pf][9],deltaR_calovsPF[small_calo][small_pf][10],
			deltaR_calovsPF[small_calo][small_pf][11],deltaR_calovsPF[small_calo][small_pf][12],
			deltaR_calovsPF[small_calo][small_pf][13],deltaR_calovsPF[small_calo][small_pf][14],
			deltaR_calovsPF[small_calo][small_pf][15]
			);

	//now we have the smallest value, lets get the delta R of that particular matched jets 

	if(jet80_1){
	  hCaloPFCorr[2][centBin]->Fill((Float_t) pfjet_pt/calojet_pt);
	  hCaloPFpt[2][centBin]->Fill(calojet_pt,pfjet_pt);
	  hCaloPFCorr_pt[2][centBin]->Fill(calojet_pt,(Float_t)pfjet_pt/calojet_pt);
	  hCalo[2][centBin]->Fill(calojet_pt);
	  hPF[2][centBin]->Fill(pfjet_pt);	    
	  hDeltaR_deltapT[2][centBin]->Fill(deltaRCaloPF,deltapT);
	}
	if(jet65_1 && !jet80_1){
	  hCaloPFCorr[1][centBin]->Fill((Float_t) pfjet_pt/calojet_pt);
	  hCaloPFpt[1][centBin]->Fill(calojet_pt,pfjet_pt);
	  hCaloPFCorr_pt[1][centBin]->Fill(calojet_pt,(Float_t)pfjet_pt/calojet_pt);
	  hCalo[1][centBin]->Fill(calojet_pt);
	  hPF[1][centBin]->Fill(pfjet_pt);	    
	  hDeltaR_deltapT[1][centBin]->Fill(deltaRCaloPF,deltapT);
	}
	if(jet55_1 && !jet65_1 && !jet80_1){
	  hCaloPFCorr[0][centBin]->Fill((Float_t) pfjet_pt/calojet_pt);
	  hCaloPFpt[0][centBin]->Fill(calojet_pt,pfjet_pt);
	  hCaloPFCorr_pt[0][centBin]->Fill(calojet_pt,(Float_t)pfjet_pt/calojet_pt);
	  hCalo[0][centBin]->Fill(calojet_pt);
	  hPF[0][centBin]->Fill(pfjet_pt);	    
	  hDeltaR_deltapT[0][centBin]->Fill(deltaRCaloPF,deltapT);
	}
	smallDeltaR = 10;
	deltaR_calovsPF[small_calo][small_pf][0] = 100; // by setting this value you effectively remove that calo and pf jet for further matching 

      }// running it for the number of calo jets: 

    }// event loop 

  }// radius loop

  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_calo_pf_jet_correlation_deltaR_0p%d_ak%s%d_%d_%d.root",deltaR,algo,radius,date.GetDate(),endfile),"RECREATE");
  f.cd();
  matchJets->Write();
  matchJets->Print();
  for(int i = 0;i<nbins_cent;++i){

    for(int a = 0;a<TrigValue-1;++a){

      hCaloPFCorr[a][i]->Write();
      hCaloPFCorr[a][i]->Print();
      hCalo[a][i]->Write();
      hCalo[a][i]->Print();
      hPF[a][i]->Write();
      hPF[a][i]->Print();
      hCaloPFpt[a][i]->Write();
      hCaloPFpt[a][i]->Print();
      hCaloPFCorr_pt[a][i]->Write();
      hCaloPFCorr_pt[a][i]->Print();
      hDeltaR_deltapT[a][i]->Write();
      hDeltaR_deltapT[a][i]->Print();
    }
    
  }

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}
