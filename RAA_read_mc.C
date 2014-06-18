// Raghav Kunnawalkam Elayavalli
// June 5th 2014
// CERN

// 
// read all the MC files for PbPb and pp and make the required histograms for the analysis. 
// need to follow the same cuts used in the data analysis here as well. 
// 

#include <iostream>
#include <stdio.h>
#include <fstream>

  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

//static const int nAlgos = 9;
//static const int BinLabelN = 11;
//remember to change this to run akPu3PF for pPb and akVs3PF for Pbpb datasets. or just create a separate header file which will be way easier. 
//static const char *algoName[nAlgos] = { "", "icPu5", "akPu2PF", "akPu3PF", "akPu4PF", "akPu5PF" , "akPu2Calo", "akPu3Calo", "akPu4Calo" };
//static const char *algoNamePP[nAlgos] = { "", "icPu5", "ak2PF", "ak3PF", "ak4PF", "ak5PF" , "ak2Calo", "ak3Calo", "ak4Calo" };
//static const char *algoNameGen[nAlgos] = { "", "icPu5", "akPu2PF", "akVs3PF", "akPu4PF", "akPu2PF", "akPu3PF", "akPu4PF" };
//static const char *BinLabel[BinLabelN] = {"100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180", "180-200", "200-240","240-300" };


// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}

class JetData
{
public:
  JetData(char *fileName, char *jetTree, char *genJetTree, bool loadGenJet = 0,bool isPbPb = 0) {
		cout <<"Open "<<fileName<<endl;
		tFile = new TFile(fileName,"read");
		tEvt = (TTree*)tFile->Get("hiEvtAnalyzer/HiTree");
		tSkim = (TTree*)tFile->Get("skimanalysis/HltTree");
		tJet = (TTree*)tFile->Get(jetTree);
		tJet->SetBranchAddress("jtpt" , jtpt );
		tJet->SetBranchAddress("trackMax" , trackMax );
		tJet->SetBranchAddress("chargedMax",chargedMax);
		tJet->SetBranchAddress("chargedSum",chargedSum);
		tJet->SetBranchAddress("neutralMax",neutralMax);
		tJet->SetBranchAddress("neutralSum",neutralSum);
		tJet->SetBranchAddress("refpt", refpt);
		tJet->SetBranchAddress("nref" ,&njets);
		tJet->SetBranchAddress("jteta", jteta);
		tJet->SetBranchAddress("jtm",jtmass);
		tJet->SetBranchAddress("pthat",&pthat);
		if (loadGenJet) tGenJet = (TTree*)tFile->Get(genJetTree);
		if (loadGenJet) tGenJet->SetBranchAddress("ngen" ,&ngen);
		if (loadGenJet) tGenJet->SetBranchAddress("genpt", genpt);
		if (loadGenJet) tGenJet->SetBranchAddress("gensubid", gensubid);
		tEvt->SetBranchAddress("hiBin",&bin);
		tEvt->SetBranchAddress("vz",&vz);
		tSkim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
		if(isPbPb) tSkim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
		else tSkim->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
		tJet->AddFriend(tEvt);
		tJet->AddFriend(tSkim);
	};
	TFile *tFile;
	TTree *tJet;
	TTree *tGenJet;
	TTree *tEvt;
	TTree* tSkim;
	float jtpt[1000];
	float refpt[1000];
	float jteta[1000];
	float jtmass[1000];
	float trackMax[1000];
	float chargedMax[1000];
        float neutralMax[1000];
  float chargedSum[1000];
  float neutralSum[1000];
	float genpt[1000];
	int gensubid[1000];
	float vz;
	float pthat;
	int njets;
	int ngen;
	int bin;     
	int pHBHENoiseFilter;
	int pPAcollisionEventSelectionPA;
	int pcollisionEventSelection;
};


using namespace std;

void RAA_read_mc(int radius = 3, char *algo = "Pu"){
   
  TStopwatch timer;
  timer.Start();

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  cout<<"Running for Radius = "<< radius<<" and Algorithm is "<<algo<<endl;
 

  const int nbins_pthat = 5;
  Double_t boundaries_pthat[nbins_pthat+1];
  char *fileName_pthat[nbins_pthat+1];
  Double_t xsection[nbins_pthat+1];
  Double_t entries[nbins_pthat];

  boundaries_pthat[0]=30;
  fileName_pthat[0] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet30_Embedded_FOREST_STARTHI53_LV1_CMSSW_5_3_16_v0_Track8_Jet27_mergev1/0.root";
  xsection[0]= 1.079e-02;
  entries[0] = 274742;
  
  boundaries_pthat[1]=50;
  fileName_pthat[1] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet50_Embedded_FOREST_STARTHI53_LV1_CMSSW_5_3_16_v0_Track8_Jet27_mergedv1/0.root";
  xsection[1]= 1.021e-03;
  entries[1] = 227776;
  
  boundaries_pthat[2]=80;
  fileName_pthat[2] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root";
  xsection[2]= 9.913e-05;
  entries[2] = 209137;
  
  boundaries_pthat[3]=100;
  fileName_pthat[3] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHat_v0/0.root";
  xsection[3]= 3.069e-05 ;
  entries[3] = 188218;
  
  boundaries_pthat[4]=120;
  fileName_pthat[4] = "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0/0.root";
  xsection[4]= 1.128e-05;
  entries[4] = 467004;


  xsection[5] = 0;
  boundaries_pthat[5]=1000;

  // Vertex & centrality reweighting for PbPb
  TF1 *fVz;
  TF1* fCentralityWeight;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);
  fCentralityWeight = new TF1("fCentralityWeight","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
  fCentralityWeight->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04);


  const int nbinsPP_pthat = 8;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  char *fileNamePP_pthat[nbinsPP_pthat+1];
  Double_t xsectionPP[nbinsPP_pthat+1];
  
  boundariesPP_pthat[0]=15;
  fileNamePP_pthat[0]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt15/HiForest_v81_merged01/pt15_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[0]= 1.079e-02;
  xsectionPP[0]= 2.034e-01;
  
  boundariesPP_pthat[1]=30;
  fileNamePP_pthat[1]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt30/HiForest_v81_merged01/pt30_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[1]= 1.021e-03;
  xsectionPP[1]= 1.075e-02;
  
  boundariesPP_pthat[2]=50;
  fileNamePP_pthat[2]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt50/HiForest_v81_merged01/pt50_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[2]= 9.913e-05;
  xsectionPP[2]= 1.025e-03;
  
  boundariesPP_pthat[3]=80;
  fileNamePP_pthat[3]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt80/HiForest_v81_merged01/pt80_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[3]= 1.128e-05;
  xsectionPP[3]= 9.865e-05;
  
  boundariesPP_pthat[4]=120;
  fileNamePP_pthat[4]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt120/HiForest_v81_merged01/pt120_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[4]= 1.470e-06;
  xsectionPP[4]= 1.129e-05;
  
  boundariesPP_pthat[5]=170;
  fileNamePP_pthat[5]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt170/HiForest_v81_merged01/pt170_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[5]= 5.310e-07;
  xsectionPP[5]= 1.465e-06;
  
  boundariesPP_pthat[6]=220;
  fileNamePP_pthat[6]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt220/HiForest_v81_merged01/pt220_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[6]= 1.192e-07;	
  xsectionPP[6]= 2.837e-07;
  
  boundariesPP_pthat[7]=280;
  fileNamePP_pthat[7]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod22/Signal_Pythia_pt280/HiForest_v81_merged01/pt280_pp2013_P01_prod22_v81_merged_forest_0.root";
  //xsectionPP[7]= 3.176e-08;
  xsectionPP[7]= 5.323e-08;
  
  xsectionPP[8] = 0;
  boundariesPP_pthat[8]=1000;
  
  
  static const int nbins_cent = 6;
  Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
  //now we have to multiply by 5, since centrality goes from 0-200. 
  Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };

  TH1F *hpbpb_gen[nbins_cent+1],*hpbpb_reco[nbins_cent+1];
  TH2F *hpbpb_matrix[nbins_cent+1];
  //TH2F *hpbpb_response[nbins_cent+1];
  TH1F *hpbpb_mcclosure_data[nbins_cent+1];

  for(int i = 0;i<nbins_cent;i++){

    hpbpb_gen[i] = new TH1F(Form("hpbpb_gen_cent%d",i),Form("Gen refpt %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_reco[i] = new TH1F(Form("hpbpb_reco_cent%d",i),Form("Reco jtpt %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    hpbpb_matrix[i] = new TH2F(Form("hpbpb_matrix_cent%d",i),Form("Matrix refpt jtpt %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
    hpbpb_mcclosure_data[i] = new TH1F(Form("hpbpb_mcclosure_data_cent%d",i),Form("data for unfolding mc closure test %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt);
    //hpbpb_response[i] = new TH2F(Form("hpbpb_response_cent%d",i),Form("response jtpt refpt %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  }
  
  hpbpb_gen[nbins_cent] = new TH1F(Form("hpbpb_gen_cent%d",nbins_cent),"Gen refpt 0-200 cent",nbins_pt,boundaries_pt);
  hpbpb_reco[nbins_cent] = new TH1F(Form("hpbpb_reco_cent%d",nbins_cent),"Reco jtpt 0-200 cent",nbins_pt,boundaries_pt);
  hpbpb_matrix[nbins_cent] = new TH2F(Form("hpbpb_matrix_cent%d",nbins_cent),"Matrix refpt jtpt 0-200 cent",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  hpbpb_mcclosure_data[nbins_cent] = new TH1F(Form("hpbpb_mcclosure_data_cent%d",nbins_cent),"data for unfolding mc closure test 0-200 cent",nbins_pt,boundaries_pt);
  //hpbpb_response[nbins_cent] = new TH2F(Form("hpbpb_response_cent%d",nbins_cent),"response jtpt refpt 0-200 cent",1000,0,1000,1000,0,1000);
  

  TH1F* hpp_gen = new TH1F("hpp_gen","gen refpt",nbins_pt,boundaries_pt);
  TH1F* hpp_reco = new TH1F("hpp_reco","reco jtpt",nbins_pt,boundaries_pt);
  TH2F* hpp_matrix = new TH2F("hpp_matrix","matrix refpt jtpt",nbins_pt,boundaries_pt,nbins_pt,boundaries_pt);
  //TH2F* hpp_response = new TH2F("hpp_response","response jtpt refpt",1000,0,1000,1000,0,1000);
  TH1F* hpp_mcclosure_data = new TH1F("hpp_mcclosure_data","data for unfolding mc closure test pp",nbins_pt,boundaries_pt);

  TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
  TH1F *hCentData = new TH1F("hCentData","",40,0,40);
  TH1F *hCentMC = new TH1F("hCentMC","",40,0,40);
	
  TH1F *hVzData = new TH1F("hVzData","",60,-15,15);
  TH1F *hVzMC = new TH1F("hVzMC","",60,-15,15);
  TH1F *hVzPPData = new TH1F("hVzPPData","",60,-15,15);
  TH1F *hVzPPMC = new TH1F("hVzPPMC","",60,-15,15);
  hCent->Sumw2();
  hCentData->Sumw2();
  hCentMC->Sumw2();
  hVzData->Sumw2();
  hVzMC->Sumw2();
  hVzPPData->Sumw2();
  hVzPPMC->Sumw2();
  
  // Setup jet data branches
  JetData *data[nbins_pthat]; 
  JetData *dataPP[nbinsPP_pthat];
  
  for (int i=0;i<nbins_pthat;i++) data[i] = new JetData(fileName_pthat[i],Form("ak%s%dPFJetAnalyzer/t",algo,radius),Form("ak%s%dPFJetAnalyzer/t",algo,radius),0,1);	
  for (int i=0;i<nbinsPP_pthat;i++) dataPP[i] = new JetData(fileNamePP_pthat[i],Form("ak%dPFJetAnalyzer/t",radius),Form("ak%dPFJetAnalyzer/t",radius),0,0);

  
  TH1F *hPtHat = new TH1F("hPtHat","",nbins_pthat,boundaries_pthat);
  TH1F *hPtHatRaw = new TH1F("hPtHatRaw","",nbins_pthat,boundaries_pthat);
  TH1F *hPtHatPP = new TH1F("hPtHatPP","",nbinsPP_pthat,boundariesPP_pthat);
  TH1F *hPtHatRawPP = new TH1F("hPtHatRawPP","",nbinsPP_pthat,boundariesPP_pthat);
  
  
  cout<<"reading all the pbpb mc files"<<endl;
  for (int i=0;i<nbins_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbins_pthat,boundaries_pthat);
    data[i]->tJet->Project("hPtHatTmp","pthat");
    hPtHatRaw->Add(hPtHatTmp);
    delete hPtHatTmp;
  }
  
  
  cout<<"reading all the pp mc files"<<endl;
  for (int i=0;i<nbinsPP_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbinsPP_pthat,boundariesPP_pthat);
    dataPP[i]->tJet->Project("hPtHatTmp","pthat");
    hPtHatRawPP->Add(hPtHatTmp);
    delete hPtHatTmp;
  }
  
  hPtHatRawPP->Print("base");
  
  // fill PbPb MC 
  
  for (int i=0;i<nbins_pthat;i++) {
    if (xsection[i]==0) continue;
    cout <<"Loading pthat"<<boundaries_pthat[i]
	 <<" sample, cross section = "<<xsection[i]
	 << Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[i],boundaries_pthat[i+1])<<endl;
    cout<<data[i]->tJet->GetEntries()<<endl;
    for (Long64_t jentry2=0; jentry2<data[i]->tJet->GetEntries();jentry2++) {
      //for (Long64_t jentry2=0; jentry2<100;jentry2++) {

      //cout<<"hi"<<endl;
      data[i]->tEvt->GetEntry(jentry2);
      data[i]->tJet->GetEntry(jentry2);
      //data[i]->tGenJet->GetEntry(jentry2);
      if(data[i]->pthat<boundaries_pthat[i] || data[i]->pthat>boundaries_pthat[i+1]) continue;
      //remember this cut is there because there was some rediculous values of pthats of -1 in the private production forests. 
	
      //if(jentry2%100==0)cout<<"pthat of that event = "<<data[i]->pthat<<endl;
      
      int pthatBin = hPtHat->FindBin(data[i]->pthat);
      //if(jentry2%100==0)cout<<"pthatBin = "<<pthatBin<<endl;
      
      //cout<<xsection[pthatBin-1]-xsection[pthatBin]<<endl;
      //cout<<"nentries = "<<hPtHatRaw->GetBinContent(pthatBin)<<endl;
      double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/hPtHatRaw->GetBinContent(pthatBin);
      //double scale = (double)(xsection[pthatBin-1]-xsection[pthatBin])/entries[i];
      
      if(fabs(data[i]->vz)>15) continue;
      int cBin = hCent->FindBin(data[i]->bin)-1;
      //int cBin = nbins_cent-1;
      double weight_cent=1;
      double weight_pt=1;
      double weight_vz=1;
      
      //weight_cent = fCentralityWeight->Eval(data[i]->bin);
      weight_vz = fVz->Eval(data[i]->vz);
      hCentMC->Fill(data[i]->bin,scale*weight_cent*weight_vz);
      hVzMC->Fill(data[i]->vz,scale*weight_cent*weight_vz);
      if (cBin>=nbins_cent) continue;
      if (cBin==-1) continue;
      hPtHat->Fill(data[i]->pthat,scale*weight_cent*weight_vz);

      if(scale*weight_cent*weight_vz <=0 ) cout<<"RED FLAG RED FLAG RED FLAG"<<endl;
      
      //cout<<"scale = "<<scale<<endl;
      
      /*
	int hasLeadingJet = 0;
	for (int k= 0; k < data[i]->njets; k++) { 
	if ( data[i]->jteta[k]  > 2. || data[i]->jteta[k] < -2. ) continue;
	if ( data[i]->jtpt[k]>100) {
	hasLeadingJet = 1;
	}
	break;
				 
	}
	if (hasLeadingJet == 0) continue;
      */

      for (int k= 0; k < data[i]->njets; k++) { 
	//int subEvt=-1;
	if ( data[i]->refpt[k]  < 30. ) continue;
	if ( data[i]->jteta[k]  > 2. || data[i]->jteta[k] < -2. ) continue;
	if ( data[i]->chargedMax[k]/data[i]->jtpt[k]<0.01) continue;
	if ( data[i]->neutralMax[k]/TMath::Max(data[i]->chargedSum[k],data[i]->neutralSum[k]) < 0.975)continue;

	//for (int l= 0; l< data[i]->ngen;l++) {
	//  if (data[i]->refpt[k]==data[i]->genpt[l]) {
	//    subEvt = data[i]->gensubid[l];
	//    break;
	//  } 
	//}
	//if (subEvt!=0) continue;
	//if (uhist[cBin]->hMeasMatch!=0) {
	//   int ptBinNumber = uhist[cBin]->hMeasMatch->FindBin(data[i]->jtpt[k]);
	//   int ratio = uhist[cBin]->hMeasMatch->GetBinContent(ptBinNumber);
	//if (ratio!=0) weight_pt = 1./ratio;
	//}
	//if (!isMC||jentry2<data[i]->tJet->GetEntries()/2.) {
	//cout<<"going to fill the histograms now"<<endl;
	//cout<<"fvz = "<<weight_vz<<endl;
	
	//hpbpb_response[cBin]->Fill(data[i]->jtpt[k],data[i]->refpt[k],scale*weight_vz);
	hpbpb_matrix[cBin]->Fill(data[i]->refpt[k],data[i]->jtpt[k],scale*weight_vz);
	hpbpb_gen[cBin]->Fill(data[i]->refpt[k],scale*weight_vz);
	hpbpb_reco[cBin]->Fill(data[i]->jtpt[k],scale*weight_vz);
	
	//hpbpb_response[nbins_cent]->Fill(data[i]->jtpt[k],data[i]->refpt[k],scale*weight_vz);
	hpbpb_matrix[nbins_cent]->Fill(data[i]->refpt[k],data[i]->jtpt[k],scale*weight_vz);
	hpbpb_gen[nbins_cent]->Fill(data[i]->refpt[k],scale*weight_vz);
	hpbpb_reco[nbins_cent]->Fill(data[i]->jtpt[k],scale*weight_vz);

	if (jentry2>data[i]->tJet->GetEntries()/2.) {
	  hpbpb_mcclosure_data[cBin]->Fill(data[i]->jtpt[k],scale*weight_vz);
	  hpbpb_mcclosure_data[nbins_cent]->Fill(data[i]->jtpt[k],scale*weight_vz);
	}
	//  uhist[cBin]-> hGen->Fill(data[i]->refpt[k],scale*weight_vz);   
	//  uhist[cBin]-> hMeas->Fill(data[i]->jtpt[k],scale*weight_vz);  	 
	//uhist[cBin]-> hMeasJECSys->Fill(data[i]->jtpt[k]*(1.+0.02/nbins_cent*(nbins_cent-i)),scale*weight_cent*weight_pt*weight_vz); 
	
	
	//}
      }//njets loop
      //uhist[cBin]->hGen->Print("base");
      
    }//nentry loop
    
  }//ptbins loop
  
  
  // Vertex reweighting for pp
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
  
  // fill pp MC
  for (int i=0;i<nbinsPP_pthat;i++) {
    if (xsectionPP[i]==0) continue;
    //float scale=(xsectionPP[i]-xsectionPP[i+1])/dataPP[i]->tJet->GetEntries(Form("pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1])); 
    cout <<"Loading PP pthat"<<boundariesPP_pthat[i]
	 <<" sample, cross section = "<<xsectionPP[i]
	 << Form(" pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1])<<endl;
    //cout<<""<<endl;
    for (Long64_t jentry2=0; jentry2<dataPP[i]->tJet->GetEntries();jentry2++) {
      //for (Long64_t jentry2=0; jentry2<10;jentry2++) {
      dataPP[i]->tEvt->GetEntry(jentry2);
      dataPP[i]->tJet->GetEntry(jentry2);
      //dataPP[i]->tGenJet->GetEntry(jentry2);
      //if(dataPP[i]->pthat<boundariesPP_pthat[i] || dataPP[i]->pthat>boundariesPP_pthat[i+1]) continue;
      //if(dataPP[i]->bin<=28) continue;
      int pthatBin = hPtHatPP->FindBin(dataPP[i]->pthat);
      float scalepp = (xsectionPP[pthatBin-1]-xsectionPP[pthatBin])/hPtHatRawPP->GetBinContent(pthatBin);
      if(fabs(dataPP[i]->vz)>15) continue;
      double weight_cent=1;
      double weight_pt=1;
      double weight_vz=1;
      if(!dataPP[i]->pPAcollisionEventSelectionPA || !dataPP[i]->pHBHENoiseFilter) continue;
      
      weight_vz = fVzPP->Eval(dataPP[i]->vz);
      //if (weight_vz>5||weight_vz<0.5) cout <<dataPP[i]->vz<<" "<<weight_vz<<endl;
      //weight_vz = 1;
      hPtHatPP->Fill(dataPP[i]->pthat,scalepp*weight_vz);
      int hasLeadingJet = 0;
      hVzPPMC->Fill(dataPP[i]->vz,scalepp*weight_vz);
      /*
	for (int k= 0; k < dataPP[i]->njets; k++) { 
	if ( dataPP[i]->jteta[k]  > 2. || dataPP[i]->jteta[k] < -2. ) continue;
	if ( dataPP[i]->jtpt[k]>100) {
	hasLeadingJet = 1;
	}
	break;
	
	}
	if (hasLeadingJet == 0) continue;
      */

      for (int k= 0; k < dataPP[i]->njets; k++) { 
	int subEvt=-1;
	if ( dataPP[i]->refpt[k]  < 30. ) continue;
	if ( dataPP[i]->jteta[k]  > 2. || dataPP[i]->jteta[k] < -2. ) continue;
	if ( dataPP[i]->chargedMax[k]/dataPP[i]->jtpt[k]<0.01) continue;
	if ( dataPP[i]->neutralMax[k]/TMath::Max(dataPP[i]->chargedSum[k],dataPP[i]->neutralSum[k]) < 0.975)continue;
	//if ( dataPP[i]->neu)

	//if (uhist[nbins_cent]->hMeasMatch!=0) {
	//   int ptBinNumber = uhist[nbins_cent]->hMeasMatch->FindBin(dataPP[i]->jtpt[k]);
	//   int ratio = uhist[nbins_cent]->hMeasMatch->GetBinContent(ptBinNumber);
	//if (ratio!=0) weight_pt = 1./ratio;
	//}
					
	//if (!isMC||jentry2<dataPP[i]->tJet->GetEntries()/2.) {
	//if(!isMC){
	//hpp_response->Fill(dataPP[i]->jtpt[k],dataPP[i]->refpt[k],scalepp*weight_vz);
	hpp_matrix->Fill(dataPP[i]->refpt[k],dataPP[i]->jtpt[k],scalepp*weight_vz);
	hpp_gen->Fill(dataPP[i]->refpt[k],scalepp*weight_vz);   
	hpp_reco->Fill(dataPP[i]->jtpt[k],scalepp*weight_vz);
	
	//}	  
	if (jentry2>dataPP[i]->tJet->GetEntries()/2.)
	  hpp_mcclosure_data->Fill(dataPP[i]->jtpt[k],scale*weight_vz);
	
	//uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt[k],scale*weight_vz);   
	//uhist[nbins_cent]-> hMeas->Fill(dataPP[i]->jtpt[k],scale*weight_vz); 
	//}
      }//njet loop     
    }//nentry loop
  }//ptbins loop
  
  TDatime date;

  //declare the output file 
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Output/PbPb_pp_mc_ak%d_%s_%d_chMax_12003cut.root",radius,algo,date.GetDate()),"RECREATE");
  f.cd();

  for(int i = 0;i<=nbins_cent;i++){
    
    //hpbpb_gen[i] = (TH1F*)hpbpb_gen[i]->Rebin(nbins_pt,Form("hpbpb_gen_cent%d",i),boundaries_pt);
    divideBinWidth(hpbpb_gen[i]);
    //hpbpb_reco[i] = (TH1F*)hpbpb_reco[i]->Rebin(nbins_pt,Form("hpbpb_reco_cent%d",i),boundaries_pt);
    divideBinWidth(hpbpb_reco[i]);

    divideBinWidth(hpbpb_mcclosure_data[i]);

    //hpbpb_matrix[i] = (TH2F*)hpbpb_matrix[i]->Rebin(nbins_pt,Form("Matrix  refpt jtpt %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),boundaries_pt);
    //hpbpb_response[i] = (TH2F*)hpbpb_response[i]->Rebin(nbins_pt,Form("Response jtpt refpt %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),boundaries_pt);

  }

  divideBinWidth(hpp_gen);
  divideBinWidth(hpp_reco);
  divideBinWidth(hpp_mcclosure_data);

  

  for(int i = 0;i<=nbins_cent;i++){
    
    hpbpb_gen[i]->Write();
    hpbpb_gen[i]->Print("base");
    hpbpb_reco[i]->Write();
    hpbpb_reco[i]->Print("base");
    hpbpb_matrix[i]->Write();
    hpbpb_matrix[i]->Print("base");
    hpbpb_mcclosure_data[i]->Write();
    hpbpb_mcclosure_data[i]->Print("base");
    //hpbpb_response[i]->Write();
    //hpbpb_response[i]->Print("base");
  }

  hpp_gen->Write();
  hpp_gen->Print("base");
  hpp_reco->Write();
  hpp_reco->Print("base");
  hpp_matrix->Write();
  hpp_matrix->Print("base");
  hpp_mcclosure_data->Write();
  hpp_mcclosure_data->Print("base");
  //hpp_response->Write();
  //hpp_response->Print("base");
  
  f.Write();
  f.Close();
  
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<timer.CpuTime()<<endl;
  cout<<"Real time (min) = "<<timer.RealTime()<<endl;
  
}
