// Raghav Kunnawalkam Elayavalli
// June 4th 2014
// CERN

// 
// read all the data files and create the required histograms needed for analysis. I'll create a new macro to read only the MC files. 
// 
// most of this macro is taken from the merge_pbpb_pp_HLT.C inside WORK/CMSSW_6_0_0/src/RAA/



#include <iostream>
#include <stdio.h>
#include <fstream>

  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

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

void RAA_read_data(int radius = 3, char *algo = "Vs"){

  
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  // number convension:
  // 0 - MB
  // 1 - 55 or 65
  // 2 - 80 or 95 
  // 80 is the unprescaled trigger - yes
  // 
  
  //data files - PbPb 
  //TFile *fpbpb0_old = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root");
  //this one is actually new, but im not going to change the title of the TFile for now. 
  TFile *fpbpb0_old = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_lumi2_FOREST_TRY2merged/0.root");

  TFile *fpbpb1_old = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet55or65_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21/0.root");
  TFile *fpbpb1 = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet55or65_GR_R_53_LV6_25Mar2014_0200CET_Track8_Jet26_TRY2_full/0.root");

  TFile *fpbpb2_old = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21/0.root");
  TFile *fpbpb2 = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_25Mar2014_0200CET_Track8_Jet26_full/0.root");

  // data files - pp 
  TFile *fpp1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet80_v2.root");
  TFile *fpp2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet40_v2.root");

  // grab the trees from the data files. 
  TTree *jetpbpb0_old = (TTree*)fpbpb0_old->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb1 = (TTree*)fpbpb1->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb1_old = (TTree*)fpbpb1_old->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb2 = (TTree*)fpbpb2->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb2_old = (TTree*)fpbpb2_old->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));

  TTree *evtpbpb0_old = (TTree*)fpbpb0_old->Get("hiEvtAnalyzer/HiTree");
  TTree *evtpbpb1 = (TTree*)fpbpb1->Get("hiEvtAnalyzer/HiTree");
  TTree *evtpbpb1_old = (TTree*)fpbpb1_old->Get("hiEvtAnalyzer/HiTree");
  TTree *evtpbpb2 = (TTree*)fpbpb2->Get("hiEvtAnalyzer/HiTree");  
  TTree *evtpbpb2_old = (TTree*)fpbpb2_old->Get("hiEvtAnalyzer/HiTree");

  TTree* hltpbpb0_old = (TTree*)fpbpb0_old->Get("hltanalysis/HltTree");
  TTree* hltpbpb1 = (TTree*)fpbpb1->Get("hltanalysis/HltTree");
  TTree* hltpbpb1_old = (TTree*)fpbpb1_old->Get("hltanalysis/HltTree");
  TTree* hltpbpb2 = (TTree*)fpbpb2->Get("hltanalysis/HltTree");
  TTree* hltpbpb2_old = (TTree*)fpbpb2_old->Get("hltanalysis/HltTree");

  TTree* skmpbpb0_old = (TTree*)fpbpb0_old->Get("skimanalysis/HltTree");
  TTree* skmpbpb1 = (TTree*)fpbpb1->Get("skimanalysis/HltTree");
  TTree* skmpbpb1_old = (TTree*)fpbpb1_old->Get("skimanalysis/HltTree");
  TTree* skmpbpb2 = (TTree*)fpbpb2->Get("skimanalysis/HltTree");
  TTree* skmpbpb2_old = (TTree*)fpbpb2_old->Get("skimanalysis/HltTree");

  TTree* hltobjpbpb0_old = (TTree*)fpbpb0_old->Get("hltobject/jetObjTree");
  TTree* hltobjpbpb1 = (TTree*)fpbpb1->Get("hltobject/jetObjTree");
  TTree* hltobjpbpb1_old = (TTree*)fpbpb1_old->Get("hltobject/jetObjTree");
  TTree* hltobjpbpb2 = (TTree*)fpbpb2->Get("hltobject/jetObjTree");
  TTree* hltobjpbpb2_old = (TTree*)fpbpb2_old->Get("hltobject/jetObjTree");

  jetpbpb0_old->AddFriend(evtpbpb0_old);
  jetpbpb1->AddFriend(evtpbpb1);
  jetpbpb1_old->AddFriend(evtpbpb1_old);
  jetpbpb2->AddFriend(evtpbpb2);
  jetpbpb2_old->AddFriend(evtpbpb2_old);

  jetpbpb0_old->AddFriend(hltpbpb0_old);
  jetpbpb1->AddFriend(hltpbpb1);
  jetpbpb1_old->AddFriend(hltpbpb1_old);
  jetpbpb2->AddFriend(hltpbpb2);
  jetpbpb2_old->AddFriend(hltpbpb2_old);

  jetpbpb0_old->AddFriend(skmpbpb0_old);
  jetpbpb1->AddFriend(skmpbpb1);
  jetpbpb1_old->AddFriend(skmpbpb1_old);
  jetpbpb2->AddFriend(skmpbpb2);
  jetpbpb2_old->AddFriend(skmpbpb2_old);

  jetpbpb0_old->AddFriend(hltobjpbpb0_old);
  jetpbpb1->AddFriend(hltobjpbpb1);
  jetpbpb1_old->AddFriend(hltobjpbpb1_old);
  jetpbpb2->AddFriend(hltobjpbpb2);
  jetpbpb2_old->AddFriend(hltobjpbpb2_old);


  //do it for the pp - need to check up on this. 
  TTree *jetpp1_v2 = (TTree*)fpp1_v2->Get(Form("jetR%d",radius));
  TTree *jetpp2_v2 = (TTree*)fpp2_v2->Get(Form("jetR%d",radius));

  TTree *evtpp1_v2 = (TTree*)fpp1_v2->Get("evt");
  TTree *evtpp2_v2 = (TTree*)fpp2_v2->Get("evt");

  jetpp1_v2->AddFriend(evtpp1_v2);
  jetpp2_v2->AddFriend(evtpp2_v2);


  //get all the pp spectra here: 
  TCut pp3 = "abs(eta)<2&&jet40&&!jet60&&!jet80&&(chMax/pt)>0.01&&(neMax/TMath::Max(chSum,neSum))>=0.975";
  
  TH1F *hpp1 = new TH1F("hpp1","Spectra from Jet 80",1000,0,1000);
  TH1F *hpp2 = new TH1F("hpp2","Spectra from Jet 60 & !Jet80",1000,0,1000);
  TH1F *hpp3 = new TH1F("hpp3","Spectra from Jet 40 & !Jet60 & !Jet80",1000,0,1000);
  TH1F *hppComb = new TH1F("hppComb","Combined spectra according to 12003 method",1000,0,1000);
  
  //get the prescl factor information. 
  //Float_t presclpbpb3 = (Float_t)jetpbpb1_v2->GetEntries("jet80")/jetpbpb1_v2->GetEntries("jet55&&jet80");
  //cout<<"pbpb prescl3 = "<<presclpbpb3<<endl;//1.99871
  //Float_t presclpp3 = (Float_t)jetpp1_v2->GetEntries("jet80")/jetpp1_v2->GetEntries("jet40&&jet80");
  //cout<<"pp prescl3 = "<<presclpp3<<endl; //9.24968

  //root [9] (Float_t)jet->GetEntries("HLT_HIJet80_v1")/jet->GetEntries("HLT_HIJet80_v1&&HLT_HIJet55_v1")
  //(double)2.34995051108819153e+00
  //ive commented this below - to just check for the pbpb histograms to load. 
  
  // include the neutralMax/ max(chargedSum, neutralSum)>0.975 cut here

  // because whatever cut we use in PbPb, we have to use in pp. 

  jetpp1_v2->Project("hpp1","pt","abs(eta)<2&&jet80&&(chMax/pt)>0.01&&(neMax/TMath::Max(chSum,neSum))>=0.975");
  hpp1->Print("base");
 
  jetpp2_v2->Project("hpp2","pt","abs(eta)<2&&jet60&&!jet80&&(chMax/pt)>0.01&&(neMax/TMath::Max(chSum,neSum))>=0.975");
  hpp2->Print("base");

  jetpp2_v2->Project("hpp3","pt","9.25038"*pp3);
  // 9.25038 was the value. 
  //jetpp2_v2->Project("hpp3","pt","jet40_p"*pp3);
  hpp3->Print("base");
 
  
  hpp1->Scale(1./5300e6);//pp lumi
  hpp2->Scale(1./5300e6);
  hpp3->Scale(1./5300e6);

  hpp1->Scale(1./4);//delta eta
  hpp2->Scale(1./4);
  hpp3->Scale(1./4);

  hppComb->Add(hpp1,1);
  hppComb->Add(hpp2,1);
  hppComb->Add(hpp3,1);
  hppComb->Print("base");

  hppComb = (TH1F*)hppComb->Rebin(nbins_pt,"hppComb",boundaries_pt);
  hpp3 = (TH1F*)hpp3->Rebin(nbins_pt,"hpp3",boundaries_pt);
  hpp2 = (TH1F*)hpp2->Rebin(nbins_pt,"hpp2",boundaries_pt);
  hpp1 = (TH1F*)hpp1->Rebin(nbins_pt,"hpp1",boundaries_pt);

  divideBinWidth(hppComb);
  divideBinWidth(hpp1);
  divideBinWidth(hpp2);
  divideBinWidth(hpp3);

  //these were for doing it from the forests directly without the proper JEC's 
  //add the centrality cuts: 

  const int nbins_cent = 6;
  Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
  Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };

  //const int nbins_cent = 1;
  //Double_t boundaries_cent[nbins_cent+1] = {0,40};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
  //Double_t ncoll[nbins_cent] = { 362.24};


  //Double_t jet55or65CentWeight[nbins_cent] = {0.3734,0.2509,0.3222,0.0352,0.0066,0.0010}; //total = 0.9983
  //Double_t jet80or90CentWeight[nbins_cent] = {0.2334,0.1764,0.3601,0.1117,0.0283,0.0054}; //total = 0.9153


  // ok so this is pretty important here: 
  // the structure of the macro realies heavily on the centrality loop. so histograms arrays from 0 to nbins_cent-1 will have the spectra for the pbpb at different centrality classes. the one at nbins_cent contains the spectra for the 0-200 bin. so its the full spectra. this was done as a cross check rather than just adding the other histograms. 
  
  
  TCut pbpb0 = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIZeroBiasPizel_SingleTrack_v1&&chargedMax/jtpt>0.01";//this is just for the MB file. not really used here so far. 

  TCut pbpb1[nbins_cent+1];
  TCut pbpb2[nbins_cent+1];
  TCut pbpb3[nbins_cent+1];

  TCut pbpb80[nbins_cent+1];
  TCut pbpb65[nbins_cent+1];
  TCut pbpb55[nbins_cent+1];

  TH1F *hpbpb1[nbins_cent+1],*hpbpb2[nbins_cent+1],*hpbpb3[nbins_cent+1];
  TH1F *hpbpbComb[nbins_cent+1];
  //TH1F* htest = new TH1F("htest","",1000,0,1000);
  TH1F *hpbpb_80[nbins_cent+1],*hpbpb_65[nbins_cent+1],*hpbpb_55[nbins_cent+1]; //histos to check the separate spectra, weighted by event by event prescl 
  // I should also add the trigger objects merging method. 

  //old way of finding trigger turn on which didnt work since we dont have a good MB sample. 

  TH1F* hTurnon80_old = new TH1F("hTurnon80_old","",150,0,150);
  TH1F* hTurnon65_old = new TH1F("hTurnon65_old","",150,0,150);
  TH1F* hTurnon55_old = new TH1F("hTurnon55_old","",150,0,150);

  TH1F* hTriggerMerged_old = new TH1F("hTriggerMerged_old","",150,0,150);
  TH1F* htest80_old = new TH1F("htest80_old","",150,0,150);
  TH1F* htest65_old = new TH1F("htest65_old","",150,0,150);
  TH1F* htest55_old = new TH1F("htest55_old","",150,0,150);

  // check the trigger turn on curve from the MB file. 
  TH1F* hMB_old = new TH1F("hMB_old","",150,0,150);
  
  //jetpbpb2->Project("htest80","jtpt","HLT_HIJet80_v1");
  //htest80->Print("base");
  //jetpbpb1->Project("htest65","jtpt","HLT_HIJet65_v1&&!HLT_HIJet80_v1");
  //htest65->Print("base");
  //jetpbpb1->Project("htest55","jtpt","HLT_HIJet55_v1_Prescl*(HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1)");
  //htest55->Print("base");
  
  jetpbpb0_old->Project("hTurnon80_old","jtpt","HLT_HIJet80_v1_Prescl*HLT_HIJet80_v1");
  jetpbpb0_old->Project("hTurnon65_old","jtpt","HLT_HIJet65_v1_Prescl*HLT_HIJet65_v1");
  jetpbpb0_old->Project("hTurnon55_old","jtpt","HLT_HIJet55_v1_Prescl*HLT_HIJet55_v1");

  TCut MB_prescl = "HLT_HIMinBiasHfOrBSC_v1_Prescl*HLT_HIMinBiasHfOrBSC_v1";
  jetpbpb0_old->Project("hMB_old","jtpt","30"*MB_prescl);

  //hTurnon80->Print("base");
  //hTurnon65->Print("base");
  //hTurnon55->Print("base");

  //hTriggerMerged->Add(htest80);
  //hTriggerMerged->Add(htest65);
  //hTriggerMerged->Add(htest55);

  hTurnon80_old->Divide(hMB_old);
  hTurnon65_old->Divide(hMB_old);
  hTurnon55_old->Divide(hMB_old);

  
  //centrality loop for the pbpb files/histograms 
  for(int i = 0;i<nbins_cent;i++){

    cout<<"centrality boundary = "<<boundaries_cent[i]*5<<" - "<<boundaries_cent[i+1]*5<<endl;

    /*
      Ok so i screwed up the convention here. the tcuts and the histogram has the number appendage going from 
      1 - 80
      2 - 65
      3 - 55

      but thats the reverse for the jet trees. 
      2 - HLT_80 or HLT_95
      1 - HLT_55 or HLT_65

      so be careful when moving through the code. 
    */

    //list of cuts we are using to remove all the fake jets. 
    // chMax/jtpt>0.01 - studied
    // neMax/jtpt>0.01 - just what im trying out now - needs to be studied
    // neutralMax/TMath::Max(chargedSum,neutralSum)<0.975 used by kurt in 12003. 
    // 

    // lets talk about the actual events/jets which are irritating here. the so called supernova events. 
    // they come from smaller jet triggers but still have larger jet pt. 

    pbpb1[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&neutralMax/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb2[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet65_v1&&!HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&neutralMax/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb3[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&neutralMax/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);

    pbpb80[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&neutralMax/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb65[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet65_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&neutralMax/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);
    pbpb55[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet55_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%2.0f&&hiBin<%2.0f&&neutralMax/TMath::Max(chargedSum,neutralSum)<0.975",5*boundaries_cent[i],5*boundaries_cent[i+1]);

    hpbpb1[i] = new TH1F(Form("hpbpb1_cent%d",i),Form("spectra from Jet 80 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    //hpbpb1[i]->Print("base");
    hpbpb2[i] = new TH1F(Form("hpbpb2_cent%d",i),Form("spectra from jet 65 & !jet 80 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb3[i] = new TH1F(Form("hpbpb3_cent%d",i),Form("spectra from jet 55 & !jet 65 & !jet 80 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpbComb[i] = new TH1F(Form("hpbpbComb_cent%d",i),Form("Spectra Combined using 12003 method %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    hpbpb_80[i] = new TH1F(Form("hpbpb_80_cent%d",i),Form("Spectra from Jet80 trigger alone %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_65[i] = new TH1F(Form("hpbpb_65_cent%d",i),Form("Spectra from Jet65 trigger alone %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_55[i] = new TH1F(Form("hpbpb_55_cent%d",i),Form("Spectra from Jet55 trigger alone %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

    jetpbpb2->Project(Form("hpbpb1_cent%d",i),"jtpt",pbpb1[i]);
    hpbpb1[i]->Print("base");
    //divideBinWidth(hpbpb1);
    
    jetpbpb1->Project(Form("hpbpb2_cent%d",i),"jtpt",pbpb2[i]);
    hpbpb2[i]->Print("base");
    //divideBinWidth(hpbpb2);
    
    jetpbpb1->Project(Form("hpbpb3_cent%d",i),"jtpt","HLT_HIJet55_v1_Prescl"*pbpb3[i]);
    //jetpbpb1->Project("hpbpb3","jtpt","2.34995"*pbpb3);
    hpbpb3[i]->Print("base");
    //divideBinWidth(hpbpb3);

    jetpbpb2->Project(Form("hpbpb_80_cent%d",i),"jtpt","HLT_HIJet80_v1_Prescl"*pbpb80[i]);
    hpbpb_80[i]->Print("base");
    jetpbpb1->Project(Form("hpbpb_65_cent%d",i),"jtpt","HLT_HIJet65_v1_Prescl"*pbpb65[i]);
    hpbpb_65[i]->Print("base");
    jetpbpb1->Project(Form("hpbpb_55_cent%d",i),"jtpt","HLT_HIJet55_v1_Prescl"*pbpb55[i]);
    hpbpb_55[i]->Print("base");

    //scale the PbPb histograms before adding them
    //we have to scale them according to the lumi of the Jet80 file. 
    // HLT file  |   Lumi inverse micro barns 
    // HLT_80    |   149.382 
    // HLT_65    |   3.195
    // HLT_55    |   2.734
    // 

    // the files which im using now is only a fraction of events of that:
    // for the PbPb 55 or 65 file its 0.977
    // for the PbPb 80 file its 0.304
    // now using the full sample file 
    
    hpbpb1[i]->Scale(1./149.382e6);//respective lumi seen by the trigger all in inverse micro barns 
    hpbpb2[i]->Scale(1./3.195e6);
    hpbpb3[i]->Scale(1./2.734e6);
    
    hpbpb1[i]->Scale(1./4);//delta eta
    hpbpb2[i]->Scale(1./4);
    hpbpb3[i]->Scale(1./4);

    //hpbpb1[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i])/jet80or95CentWeight[i]);//centrality bin width and scaling by the centrality events fraction. 
    //hpbpb2[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i])/jet55or65CentWeight[i]);
    //hpbpb3[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i])/jet55or65CentWeight[i]);

    hpbpb1[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));//centrality bin width 
    hpbpb2[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));
    hpbpb3[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));

    //might have to end up adding a centrality weight - the ratio of events per centrality class. 
    
    //add the histograms  
    hpbpbComb[i]->Add(hpbpb1[i]);
    hpbpbComb[i]->Add(hpbpb2[i]);
    hpbpbComb[i]->Add(hpbpb3[i]);
    hpbpbComb[i]->Print("base");

    hpbpbComb[i] = (TH1F*)hpbpbComb[i]->Rebin(nbins_pt,Form("hpbpbComb_cent%d",i),boundaries_pt);
    hpbpb3[i] = (TH1F*)hpbpb3[i]->Rebin(nbins_pt,Form("hpbpb3_cent%d",i),boundaries_pt);
    hpbpb2[i] = (TH1F*)hpbpb2[i]->Rebin(nbins_pt,Form("hpbpb2_cent%d",i),boundaries_pt);
    hpbpb1[i] = (TH1F*)hpbpb1[i]->Rebin(nbins_pt,Form("hpbpb1_cent%d",i),boundaries_pt);

    divideBinWidth(hpbpbComb[i]);
    divideBinWidth(hpbpb3[i]);
    divideBinWidth(hpbpb2[i]);
    divideBinWidth(hpbpb1[i]);
    
  }

  // doing it the 0-200 centrality bin 

  pbpb1[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01";
  pbpb2[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet65_v1&&!HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01";
  pbpb3[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01";
  
  pbpb80[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet80_v1&&(chargedMax/jtpt>0.01)";
  pbpb65[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet65_v1&&(chargedMax/jtpt>0.01";
  pbpb55[nbins_cent] = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet55_v1&&(chargedMax/jtpt)>0.01";

  hpbpb1[nbins_cent] = new TH1F(Form("hpbpb1_cent%d",nbins_cent),"Spectra from Jet80 0-200 cent",1000,0,1000);
  //hpbpb1[i]->Print("base");
  hpbpb2[nbins_cent] = new TH1F(Form("hpbpb2_cent%d",nbins_cent),"Spectra from Jet 65 & !Jet80 0-200 cent",1000,0,1000);
  hpbpb3[nbins_cent] = new TH1F(Form("hpbpb3_cent%d",nbins_cent),"Spectra from Jet 55 & !Jet65 & !Jet80 0-200 cent",1000,0,1000);
  hpbpbComb[nbins_cent] = new TH1F(Form("hpbpbComb_cent%d",nbins_cent),"Combined Jet spectra 12003 method 0-200 cent",1000,0,1000);
  
  hpbpb_80[nbins_cent] = new TH1F(Form("hpbpb_80_cent%d",nbins_cent),"Spectra from Jet 80 alone 0-200 cent",1000,0,1000);
  hpbpb_65[nbins_cent] = new TH1F(Form("hpbpb_65_cent%d",nbins_cent),"Spectra from Jet 65 alone 0-200 cent",1000,0,1000);
  hpbpb_55[nbins_cent] = new TH1F(Form("hpbpb_55_cent%d",nbins_cent),"Spectra from Jet 55 alone 0-200 cent",1000,0,1000);
  
  jetpbpb2->Project(Form("hpbpb1_cent%d",nbins_cent),"jtpt",pbpb1[nbins_cent]);
  hpbpb1[nbins_cent]->Print("base");
  //divideBinWidth(hpbpb1);
    
  jetpbpb1->Project(Form("hpbpb2_cent%d",nbins_cent),"jtpt",pbpb2[nbins_cent]);
  hpbpb2[nbins_cent]->Print("base");
  //divideBinWidth(hpbpb2);
  
  jetpbpb1->Project(Form("hpbpb3_cent%d",nbins_cent),"jtpt","HLT_HIJet55_v1_Prescl"*pbpb3[nbins_cent]);
  //jetpbpb1->Project("hpbpb3","jtpt","2.34995"*pbpb3);
  hpbpb3[nbins_cent]->Print("base");
  //divideBinWidth(hpbpb3);
  
  //following histograms are for the trigger turnon curve. no cuts apart from the trigger selection. 
  jetpbpb2->Project(Form("hpbpb_80_cent%d",nbins_cent),"jtpt","HLT_HIJet80_v1_Prescl*HLT_HIJet80_v1");
  hpbpb_80[nbins_cent]->Print("base");
  jetpbpb1->Project(Form("hpbpb_65_cent%d",nbins_cent),"jtpt","HLT_HIJet65_v1_Prescl*HLT_HIJet65_v1");
  hpbpb_65[nbins_cent]->Print("base");
  jetpbpb1->Project(Form("hpbpb_55_cent%d",nbins_cent),"jtpt","HLT_HIJet55_v1_Prescl*HLT_HIJet55_v1");
  hpbpb_55[nbins_cent]->Print("base");
  
  //scale the PbPb histograms before adding them
  //we have to scale them according to the lumi of the Jet80 file. 
  // HLT file  |   Lumi inverse micro barns 
  // HLT_80    |   149.382 
  // HLT_65    |   3.195
  // HLT_55    |   2.734
  // 
  
  // the files which im using now is only a fraction of events of that:
  // for the PbPb 55 or 65 file its 0.977
  // for the PbPb 80 file its 0.304
  // now using the full sample file 
  
  hpbpb1[nbins_cent]->Scale(1./149.382e6);//respective lumi seen by the trigger all in inverse micro barns 
  hpbpb2[nbins_cent]->Scale(1./3.195e6);
  hpbpb3[nbins_cent]->Scale(1./2.734e6);
  
  hpbpb1[nbins_cent]->Scale(1./4);//delta eta
  hpbpb2[nbins_cent]->Scale(1./4);
  hpbpb3[nbins_cent]->Scale(1./4);
  
  //might have to end up adding a centrality weight - the ratio of events per centrality class. 
  
  //add the histograms  
  hpbpbComb[nbins_cent]->Add(hpbpb1[nbins_cent]);
  hpbpbComb[nbins_cent]->Add(hpbpb2[nbins_cent]);
  hpbpbComb[nbins_cent]->Add(hpbpb3[nbins_cent]);
  hpbpbComb[nbins_cent]->Print("base");
  
  hpbpbComb[nbins_cent] = (TH1F*)hpbpbComb[nbins_cent]->Rebin(nbins_pt,Form("hpbpbComb_cent%d",nbins_cent),boundaries_pt);
  hpbpb3[nbins_cent] = (TH1F*)hpbpb3[nbins_cent]->Rebin(nbins_pt,Form("hpbpb3_cent%d",nbins_cent),boundaries_pt);
  hpbpb2[nbins_cent] = (TH1F*)hpbpb2[nbins_cent]->Rebin(nbins_pt,Form("hpbpb2_cent%d",nbins_cent),boundaries_pt);
  hpbpb1[nbins_cent] = (TH1F*)hpbpb1[nbins_cent]->Rebin(nbins_pt,Form("hpbpb1_cent%d",nbins_cent),boundaries_pt);
  
  divideBinWidth(hpbpbComb[nbins_cent]);
  divideBinWidth(hpbpb3[nbins_cent]);
  divideBinWidth(hpbpb2[nbins_cent]);
  divideBinWidth(hpbpb1[nbins_cent]);
  
  //ok now we have the spectra for the 0-200% centrality. 
  
  
  // do the trigger object merging here: 
  // this has to be done in the event loop which means that we have to get the 
  // create the trees and set the branch address
  // jet tree

  // similarly here 0 - MB file, 1 - 55or65, 2 - 80or95
  
  //file 0:
  // jet tree
  int nrefe_0;
  float pt_0[1000];
  //float old_pt3[1000];
  float raw_0[1000];
  float eta_0[1000];
  float eta_0_CM[1000];
  float phi_0[1000];
  float chMax_0[1000];
  float trkMax_0[1000];
  float chSum_0[1000];
  float phSum_0[1000];
  float neSum_0[1000];
  float trkSum_0[1000];
  float phMax_0[1000];
  float neMax_0[1000];

  // event tree
  int evt_0;
  int run_0;
  int lumi_0;
  int hiBin_0;
  float vx_0;
  float vy_0;
  float vz_0;
  int hiNtracks_0;
  float hiHFminus_0;
  float hiHFplus_0;
  float hiHFplusEta4_0;
  float hiHFminusEta4_0;
  int pcollisionEventSelection_0;
  int pHBHENoiseFilter_0;
  int pprimaryvertexFilter_0;
  int pVertexFilterCutGplus_0;

  // trigger tree
  int L1_MB_0;
  int L1_MB_p_0;
  int jetMB_0;
  int jet55_0;
  int jet65_0;
  int jet80_0;
  int jetMB_p_0;
  int jet55_p_0;
  int jet65_p_0;
  int jet80_p_0;

  // trigger object tree - this contains the maximum value of the particular trigger object. 
  float trgObj_id_0;
  float trgObj_pt_0;
  float trgObj_eta_0;
  float trgObj_phi_0;
  float trgObj_mass_0;

  //file 1: 
  // jet tree
  int nrefe_1;
  float pt_1[1000];
  //float old_pt3[1000];
  float raw_1[1000];
  float eta_1[1000];
  float eta_1_CM[1000];
  float phi_1[1000];
  float chMax_1[1000];
  float trkMax_1[1000];
  float chSum_1[1000];
  float phSum_1[1000];
  float neSum_1[1000];
  float trkSum_1[1000];
  float phMax_1[1000];
  float neMax_1[1000];

  // event tree
  int evt_1;
  int run_1;
  int lumi_1;
  int hiBin_1;
  float vx_1;
  float vy_1;
  float vz_1;
  int hiNtracks_1;
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
  int jetMB_1;
  int jet55_1;
  int jet65_1;
  int jet80_1;
  int jetMB_p_1;
  int jet55_p_1;
  int jet65_p_1;
  int jet80_p_1;


  // trigger object tree - this contains the maximum value of the particular trigger object. 
  float trgObj_id_1;
  float trgObj_pt_1;
  float trgObj_eta_1;
  float trgObj_phi_1;
  float trgObj_mass_1;
  
  //file 2: 
  // jet tree
  int nrefe_2;
  float pt_2[1000];
  //float old_pt3[1000];
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

  // event tree
  int evt_2;
  int run_2;
  int lumi_2;
  int hiBin_2;
  float vx_2;
  float vy_2;
  float vz_2;
  int hiNtracks_2;
  float hiHFminus_2;
  float hiHFplus_2;
  float hiHFplusEta4_2;
  float hiHFminusEta4_2;
  int pcollisionEventSelection_2;
  int pHBHENoiseFilter_2;
  int pprimaryvertexFilter_2;
  int pVertexFilterCutGplus_2;

  // trigger tree
  int L1_MB_2;
  int L1_MB_p_2;
  int jetMB_2;
  int jet55_2;
  int jet65_2;
  int jet80_2;
  int jetMB_p_2;
  int jet55_p_2;
  int jet65_p_2;
  int jet80_p_2;


  // trigger object tree - this contains the maximum value of the particular trigger object. 
  float trgObj_id_2;
  float trgObj_pt_2;
  float trgObj_eta_2;
  float trgObj_phi_2;
  float trgObj_mass_2;
  
  
  //set the branch addresses:  - one of the most boring parts of the code: 
  jetpbpb0_old->SetBranchAddress("evt",&evt_0);
  jetpbpb0_old->SetBranchAddress("run",&run_0);
  jetpbpb0_old->SetBranchAddress("lumi",&lumi_0);
  jetpbpb0_old->SetBranchAddress("hiBin",&hiBin_0);
  jetpbpb0_old->SetBranchAddress("vz",&vz_0);
  jetpbpb0_old->SetBranchAddress("vx",&vx_0);
  jetpbpb0_old->SetBranchAddress("vy",&vy_0);
  jetpbpb0_old->SetBranchAddress("hiNtracks",&hiNtracks_0);
  jetpbpb0_old->SetBranchAddress("hiHFminus",&hiHFminus_0);
  jetpbpb0_old->SetBranchAddress("hiHFplus",&hiHFplus_0);
  jetpbpb0_old->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_0);
  jetpbpb0_old->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_0);
  jetpbpb0_old->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_0);
  jetpbpb0_old->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_0);
  //jetpbpb0_old->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_0);
  //jetpbpb0_old->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_0);
  
  jetpbpb0_old->SetBranchAddress("nref",&nrefe_0);
  jetpbpb0_old->SetBranchAddress("jtpt",&pt_0);
  jetpbpb0_old->SetBranchAddress("jteta",&eta_0);
  jetpbpb0_old->SetBranchAddress("jtphi",&phi_0);
  jetpbpb0_old->SetBranchAddress("rawpt",&raw_0);
  jetpbpb0_old->SetBranchAddress("chargedMax",&chMax_0);
  jetpbpb0_old->SetBranchAddress("chargedSum",&chSum_0);
  jetpbpb0_old->SetBranchAddress("trackMax",&trkMax_0);
  jetpbpb0_old->SetBranchAddress("trackSum",&trkSum_0);
  jetpbpb0_old->SetBranchAddress("photonMax",&phMax_0);
  jetpbpb0_old->SetBranchAddress("photonSum",&phSum_0);
  jetpbpb0_old->SetBranchAddress("neutralMax",&neMax_0);
  jetpbpb0_old->SetBranchAddress("neutralSum",&neSum_0);
  /*
    jetTree->SetBranchAddress("nTrk",&nTrack);
    jetTree->SetBranchAddress("trkPt",&trkPt);
    jetTree->SetBranchAddress("trkEta",&trkEta);
    jetTree->SetBranchAddress("trkPhi",&trkPhi);
    jetTree->SetBranchAddress("highPurity",&highPurity);
    jetTree->SetBranchAddress("trkDz1",&trkDz1);
    jetTree->SetBranchAddress("trkDzError1",&trkDzError1);
    jetTree->SetBranchAddress("trkDxy1",&trkDxy1);
    jetTree->SetBranchAddress("trkDxyError1",&trkDxyError1);
  */
  //jetpbpb0_old->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_0);
  //jetpbpb0_old->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_0);
  //jetpbpb0_old->SetBranchAddress("L1_ZeroBias",&L1_MB_0);
  //jetpbpb0_old->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_0);
  jetpbpb0_old->SetBranchAddress("HLT_HIJet55_v1",&jet55_0);
  jetpbpb0_old->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_0);
  jetpbpb0_old->SetBranchAddress("HLT_HIJet65_v1",&jet65_0);
  jetpbpb0_old->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_0);
  jetpbpb0_old->SetBranchAddress("HLT_HIJet80_v1",&jet80_0);
  jetpbpb0_old->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_0);
  /*
  jetpbpb0_old->SetBranchAddress("id",&trgObj_id_0);
  jetpbpb0_old->SetBranchAddress("pt",&trgObj_pt_0);
  jetpbpb0_old->SetBranchAddress("eta",&trgObj_eta_0);
  jetpbpb0_old->SetBranchAddress("phi",&trgObj_phi_0);
  jetpbpb0_old->SetBranchAddress("mass",&trgObj_mass_0);
  */

  //set the branch addresses:  - one of the most boring parts of the code: 
  jetpbpb1->SetBranchAddress("evt",&evt_1);
  jetpbpb1->SetBranchAddress("run",&run_1);
  jetpbpb1->SetBranchAddress("lumi",&lumi_1);
  jetpbpb1->SetBranchAddress("hiBin",&hiBin_1);
  jetpbpb1->SetBranchAddress("vz",&vz_1);
  jetpbpb1->SetBranchAddress("vx",&vx_1);
  jetpbpb1->SetBranchAddress("vy",&vy_1);
  jetpbpb1->SetBranchAddress("hiNtracks",&hiNtracks_1);
  jetpbpb1->SetBranchAddress("hiHFminus",&hiHFminus_1);
  jetpbpb1->SetBranchAddress("hiHFplus",&hiHFplus_1);
  jetpbpb1->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_1);
  jetpbpb1->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_1);
  jetpbpb1->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_1);
  jetpbpb1->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_1);
  //jetpbpb1->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_1);
  //jetpbpb1->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_1);
  
  jetpbpb1->SetBranchAddress("nref",&nrefe_1);
  jetpbpb1->SetBranchAddress("jtpt",&pt_1);
  jetpbpb1->SetBranchAddress("jteta",&eta_1);
  jetpbpb1->SetBranchAddress("jtphi",&phi_1);
  jetpbpb1->SetBranchAddress("rawpt",&raw_1);
  jetpbpb1->SetBranchAddress("chargedMax",&chMax_1);
  jetpbpb1->SetBranchAddress("chargedSum",&chSum_1);
  jetpbpb1->SetBranchAddress("trackMax",&trkMax_1);
  jetpbpb1->SetBranchAddress("trackSum",&trkSum_1);
  jetpbpb1->SetBranchAddress("photonMax",&phMax_1);
  jetpbpb1->SetBranchAddress("photonSum",&phSum_1);
  jetpbpb1->SetBranchAddress("neutralMax",&neMax_1);
  jetpbpb1->SetBranchAddress("neutralSum",&neSum_1);

  //jetpbpb1->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_1);
  //jetpbpb1->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_1);
  //jetpbpb1->SetBranchAddress("L1_ZeroBias",&L1_MB_1);
  //jetpbpb1->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet55_v1",&jet55_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet65_v1",&jet65_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet80_v1",&jet80_1);
  jetpbpb1->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_1);

  jetpbpb1->SetBranchAddress("id",&trgObj_id_1);
  jetpbpb1->SetBranchAddress("pt",&trgObj_pt_1);
  jetpbpb1->SetBranchAddress("eta",&trgObj_eta_1);
  jetpbpb1->SetBranchAddress("phi",&trgObj_phi_1);
  jetpbpb1->SetBranchAddress("mass",&trgObj_mass_1);

  //set the branch addresses:  - one of the most boring parts of the code: 
  jetpbpb2->SetBranchAddress("evt",&evt_2);
  jetpbpb2->SetBranchAddress("run",&run_2);
  jetpbpb2->SetBranchAddress("lumi",&lumi_2);
  jetpbpb2->SetBranchAddress("hiBin",&hiBin_2);
  jetpbpb2->SetBranchAddress("vz",&vz_2);
  jetpbpb2->SetBranchAddress("vx",&vx_2);
  jetpbpb2->SetBranchAddress("vy",&vy_2);
  jetpbpb2->SetBranchAddress("hiNtracks",&hiNtracks_2);
  jetpbpb2->SetBranchAddress("hiHFminus",&hiHFminus_2);
  jetpbpb2->SetBranchAddress("hiHFplus",&hiHFplus_2);
  jetpbpb2->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4_2);
  jetpbpb2->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4_2);
  jetpbpb2->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_2);
  jetpbpb2->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_2);
  //jetpbpb2->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_2);
  //jetpbpb2->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_2);
  
  jetpbpb2->SetBranchAddress("nref",&nrefe_2);
  jetpbpb2->SetBranchAddress("jtpt",&pt_2);
  jetpbpb2->SetBranchAddress("jteta",&eta_2);
  jetpbpb2->SetBranchAddress("jtphi",&phi_2);
  jetpbpb2->SetBranchAddress("rawpt",&raw_2);
  jetpbpb2->SetBranchAddress("chargedMax",&chMax_2);
  jetpbpb2->SetBranchAddress("chargedSum",&chSum_2);
  jetpbpb2->SetBranchAddress("trackMax",&trkMax_2);
  jetpbpb2->SetBranchAddress("trackSum",&trkSum_2);
  jetpbpb2->SetBranchAddress("photonMax",&phMax_2);
  jetpbpb2->SetBranchAddress("photonSum",&phSum_2);
  jetpbpb2->SetBranchAddress("neutralMax",&neMax_2);
  jetpbpb2->SetBranchAddress("neutralSum",&neSum_2);

  //jetpbpb2->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&jetMB_2);
  //jetpbpb2->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1_Prescl",&jetMB_p_2);
  //jetpbpb2->SetBranchAddress("L1_ZeroBias",&L1_MB_2);
  //jetpbpb2->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p_2);
  jetpbpb2->SetBranchAddress("HLT_HIJet55_v1",&jet55_2);
  jetpbpb2->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_2);
  jetpbpb2->SetBranchAddress("HLT_HIJet65_v1",&jet65_2);
  jetpbpb2->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_2);
  jetpbpb2->SetBranchAddress("HLT_HIJet80_v1",&jet80_2);
  jetpbpb2->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_2);

  jetpbpb2->SetBranchAddress("id",&trgObj_id_2);
  jetpbpb2->SetBranchAddress("pt",&trgObj_pt_2);
  jetpbpb2->SetBranchAddress("eta",&trgObj_eta_2);
  jetpbpb2->SetBranchAddress("phi",&trgObj_phi_2);
  jetpbpb2->SetBranchAddress("mass",&trgObj_mass_2);
  
  //now that we have all the branch addresses set up we can start going through the loop to look at the trigger objects 
  
  //before we go and do the trigger object merging lets do the text file information here: 
  //we need one text file per spectra, which means we need 3 -> one for each trigger object 
  //and we also need one for the high pt jets with low trigger objects 
  // and these files will have the following structure: 
  // run lumi evt HLTobjpt HLTobjeta HLTobjphi hibin jtpt jteta jtphi

  ofstream fHLT_80,fHLT_65,fHLT_55,fHLT_high;
  fHLT_80.open(Form("pbpb_%s_R%d_max_trgObj_80_jets.txt",algo,radius));
  fHLT_65.open(Form("pbpb_%s_R%d_80_max_trgObj_65_jets.txt",algo,radius));
  fHLT_55.open(Form("pbpb_%s_R%d_65_max_trgObj_55_jets.txt",algo,radius));
  fHLT_high.open(Form("pbpb_%s_R%d_abnormal_highpt_jets_maxtrigObj_jtpt_less_0.3.txt",algo,radius));
  // the actual cut we are going to be using is 
  // TMath::Max(chargedMax,neutralMax)/(TMath::Max(chargedSum,neutralSum))<0.975
  // which is pretty much the same as doing the one mentioned above.  
  // ok the cut mentioned above is not the one we are using for the analysis now. We have to check if that cut actually removes a lot of the fake events for us to have a meaningful combined jet spectra. 
  

  //neutralMax/TMath::Max(chargedSum,neutralSum)<0.975 im going to add this from Kurt used in 12003 here to see if it sort of helps. 


  //declare the histograms needed for the hpbpb_TrigObj80, hpbpb_TrigObj65, hpbpb_TrigObj55, hpbpb_TrigObjMB, hpbpb_TrigComb;
  TH1F *hpbpb_TrgObj80[nbins_cent+1], *hpbpb_TrgObj65[nbins_cent+1], *hpbpb_TrgObj55[nbins_cent+1], *hpbpb_TrgObjMB[nbins_cent+1], *hpbpb_TrgObjComb[nbins_cent+1];

  for(int i = 0;i<nbins_cent;i++){

    hpbpb_TrgObj80[i] = new TH1F(Form("hpbpb_TrgObj80_cent%d",i),Form("Spectra from Trig Object > 80 and Jet 80 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObj65[i] = new TH1F(Form("hpbpb_TrgObj65_cent%d",i),Form("Spectra from Trig Object > 65 and < 80 and Jet 65 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObj55[i] = new TH1F(Form("hpbpb_TrgObj55_cent%d",i),Form("Spectra from Trig Object > 55 and < 65 and Jet 55 %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);
    hpbpb_TrgObjMB[i] = new TH1F(Form("hpbpb_TrgObjMB_cent%d",i),"",1000,0,1000);
    hpbpb_TrgObjComb[i] = new TH1F(Form("hpbpb_TrgObjComb_cent%d",i),Form("Trig Combined Spectra using 14007 method %2.0f - %2.0f cent",5*boundaries_cent[i],5*boundaries_cent[i+1]),1000,0,1000);

  }

  hpbpb_TrgObj80[nbins_cent] = new TH1F("hpbpb_TrgObj80_6","Spectra from Trig Object > 80 and Jet 80 0-200 cent",1000,0,1000);
  hpbpb_TrgObj65[nbins_cent] = new TH1F("hpbpb_TrgObj65_6","Spectra from Trig Object > 65 and < 80 and Jet 65 0-200 cent",1000,0,1000);
  hpbpb_TrgObj55[nbins_cent] = new TH1F("hpbpb_TrgObj55_6","Spectra from Trig Object > 55 and < 65 and Jet 55 0-200 cent",1000,0,1000);
  hpbpb_TrgObjMB[nbins_cent] = new TH1F("hpbpb_TrgObjMB_6","",1000,0,1000);
  hpbpb_TrgObjComb[nbins_cent] = new TH1F("hpbpb_TrgObjComb_6","Trig combined spectra using 14007 method 0-200 cent",1000,0,1000);
    
  
  // loop for the jetpbpb1 tree 
  Long64_t nentries_jet55or65 = jetpbpb1->GetEntries();
  cout<<"nentries_jet55or65 = "<<nentries_jet55or65<<endl;
  
  for(int ij = 0;ij<nentries_jet55or65;ij++){
  //for(int ij = 0;ij<10;ij++){
    
    jetpbpb1->GetEntry(ij);
    if(ij%100000==0)cout<<"Jet 55or65 right"<<endl;
    if(ij%100000==0)cout<<ij<<": event = "<<evt_1<<"; run = "<<run_1<<endl;

    // before the cuts get the spectra for the 0-200 centrality bin without any cuts. for trigger turn on curve    
    for(int j = 0;j<nrefe_1;j++){
      
      if(jet55_1&&trgObj_pt_1>=55&&trgObj_pt_1<65) hpbpb_TrgObj55[nbins_cent]->Fill(pt_1[j],jet55_p_1);
	
      if(jet65_1&&trgObj_pt_1>=65&&trgObj_pt_1<80) hpbpb_TrgObj65[nbins_cent]->Fill(pt_1[j],jet65_p_1);
    }
    
    //include the cuts for analysis. 
    if(!pHBHENoiseFilter_1 || !pcollisionEventSelection_1) continue;
    
    if(fabs(vz_1)>15) continue;
    
    if(!jet55_1 && !jet65_1) continue;
    
    if(nbins_cent==1){
      
      //cout<<"inside 0-1"<<endl;
      for(int j = 0;j<nrefe_1;j++){
	  
	if(fabs(eta_1[j])>2) continue;
	  
	if(chMax_1[j]/pt_1[j]<0.01) continue;	
	if(neMax/TMath::Max(chSum,neSum)>=0.975) continue;
	  
	if(jet55_1&&trgObj_pt_1>=55&&trgObj_pt_1<65){ hpbpb_TrgObj55[nbins_cent-1]->Fill(pt_1[j],jet55_p_1);
	  //fHLT_55<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	}
	
	if(jet65_1&&trgObj_pt_1>=65&&trgObj_pt_1<80){ hpbpb_TrgObj65[nbins_cent-1]->Fill(pt_1[j],jet65_p_1);
	  //fHLT_65<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	}
	
	if(TMath::Max(chMax_1[j],neMax_1[j])/(TMath::Max(chSum_1[j],neSum_1[j]))>0.975)
	  fHLT_high<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	
      } // jet loop
	
    }else if(nbins_cent==6){
      
      // checking the centrality class statement
      if(hiBin_1>=5*boundaries_cent[0]&&hiBin_1<5*boundaries_cent[1]){
	//cout<<"inside 0-1"<<endl;
	for(int j = 0;j<nrefe_1;j++){
	  
	  if(fabs(eta_1[j])>2) continue;
	  
	  if(chMax_1[j]/pt_1[j]<0.01) continue;	
	  //here is where we have to change the cuts to remove the fake jets. here and above where we combined triggers using the 12003 method. 
	  if(neMax_1[j]/TMath::Max(chSum_1[j],neSum_1[j])>=0.975) continue;


	  if(jet55_1&&trgObj_pt_1>=55&&trgObj_pt_1<65){ 
	    hpbpb_TrgObj55[0]->Fill(pt_1[j],jet55_p_1);
	    //fHLT_55<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	  
	  if(jet65_1&&trgObj_pt_1>=65&&trgObj_pt_1<80){
	    hpbpb_TrgObj65[0]->Fill(pt_1[j],jet65_p_1);
	    //fHLT_65<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	  
	} // jet loop
	
      }
      if(hiBin_1>=5*boundaries_cent[1]&&hiBin_1<5*boundaries_cent[2]){
	//cout<<"inside 1-2"<<endl;
	
	for(int j = 0;j<nrefe_1;j++){
	  
	  if(fabs(eta_1[j])>2) continue;
	  
	  if(chMax_1[j]/pt_1[j]<0.01) continue;

	  if(neMax_1[j]/TMath::Max(chSum_1[j],neSum_1[j])>=0.975) continue;
	  
	  if(jet55_1&&trgObj_pt_1>=55&&trgObj_pt_1<65){ hpbpb_TrgObj55[1]->Fill(pt_1[j],jet55_p_1);
	    //fHLT_55<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	  
	  if(jet65_1&&trgObj_pt_1>=65&&trgObj_pt_1<80){ hpbpb_TrgObj65[1]->Fill(pt_1[j],jet65_p_1);
	    //fHLT_65<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	  
	} // jet loop
	
      }
      if(hiBin_1>=5*boundaries_cent[2]&&hiBin_1<5*boundaries_cent[3]){
	//cout<<"inside 2-3"<<endl;
	for(int j = 0;j<nrefe_1;j++){
	  
	  if(fabs(eta_1[j])>2) continue;
	  
	  if(chMax_1[j]/pt_1[j]<0.01) continue;

	  if(neMax_1[j]/TMath::Max(chSum_1[j],neSum_1[j])>=0.975) continue;
	  
	  if(jet55_1&&trgObj_pt_1>=55&&trgObj_pt_1<65){ hpbpb_TrgObj55[2]->Fill(pt_1[j],jet55_p_1);
	    //fHLT_55<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	  
	  if(jet65_1&&trgObj_pt_1>=65&&trgObj_pt_1<80){ hpbpb_TrgObj65[2]->Fill(pt_1[j],jet65_p_1);
	    //fHLT_65<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	} // jet loop
	
      }
      if(hiBin_1>=5*boundaries_cent[3]&&hiBin_1<5*boundaries_cent[4]){
	//cout<<"inside 3-4"<<endl;
	
	for(int j = 0;j<nrefe_1;j++){
	
	  if(fabs(eta_1[j])>2) continue;
	
	  if(chMax_1[j]/pt_1[j]<0.01) continue;

	  if(neMax_1[j]/TMath::Max(chSum_1[j],neSum_1[j])>=0.975) continue;
		
	  if(jet55_1 && trgObj_pt_1>=55 && trgObj_pt_1<65){ hpbpb_TrgObj55[3]->Fill(pt_1[j],jet55_p_1);
	    //fHLT_55<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	
	  if(jet65_1 && trgObj_pt_1>=65 && trgObj_pt_1<80){ hpbpb_TrgObj65[3]->Fill(pt_1[j],jet65_p_1);
	    //fHLT_65<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	} // jet loop
      
      }
      if(hiBin_1>=5*boundaries_cent[4] && hiBin_1<5*boundaries_cent[5]){
	//cout<<"inside 4-5"<<endl;

	for(int j = 0;j<nrefe_1;j++){
	
	  if(fabs(eta_1[j])>2) continue;	
	  if(chMax_1[j]/pt_1[j]<0.01) continue;

	  if(neMax_1[j]/TMath::Max(chSum_1[j],neSum_1[j])>=0.975) continue;
		
	  if(jet55_1 && trgObj_pt_1>=55 && trgObj_pt_1<65){ hpbpb_TrgObj55[4]->Fill(pt_1[j],jet55_p_1);
	    //fHLT_55<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	
	  if(jet65_1 && trgObj_pt_1>=65 && trgObj_pt_1<80){ hpbpb_TrgObj65[4]->Fill(pt_1[j],jet65_p_1);
	    //fHLT_65<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	} // jet loop
      
      }
      if(hiBin_1>=5*boundaries_cent[5] && hiBin_1<5*boundaries_cent[6]){
	//cout<<"inside 5-6"<<endl;

	for(int j = 0;j<nrefe_1;j++){
	
	  if(fabs(eta_1[j])>2) continue;
	
	  if(chMax_1[j]/pt_1[j]<0.01) continue;

	  if(neMax_1[j]/TMath::Max(chSum_1[j],neSum_1[j])>=0.975) continue;
		
	  if(jet55_1 && trgObj_pt_1>=55 && trgObj_pt_1<65){ hpbpb_TrgObj55[5]->Fill(pt_1[j],jet55_p_1);
	    //fHLT_55<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	
	  if(jet65_1 && trgObj_pt_1>=65 && trgObj_pt_1<80){ hpbpb_TrgObj65[5]->Fill(pt_1[j],jet65_p_1);
	    //fHLT_65<<run_1<<" "<<evt_1<<" "<<lumi_1<<" "<<trgObj_pt_1<<" "<<trgObj_eta_1<<" "<<trgObj_phi_1<<" "<<hiBin_1<<" "<<pt_1[j]<<" "<<eta_1[j]<<" "<<phi_1[j]<<endl;
	  }
	} // jet loop
      
      }

    }
    
  } // nentries_jet55or65 loop
  
  
    //loop for jetpbpb2 tree
  Long64_t nentries_jet80or95 = jetpbpb2->GetEntries();
  cout<<"nentries_jet80or95 = "<<nentries_jet80or95<<endl;
  
  for(int ij = 0;ij<nentries_jet80or95;ij++){
  //for(int ij = 0;ij<10;ij++){
    
    jetpbpb2->GetEntry(ij);
    //if(ij%100000==0)cout<<"Jet 80 or 95 file"<<endl;
    if(ij%100000==0)cout<<ij<<": event = "<<evt_2<<", run = "<<run_2<<endl;

    for(int j = 0;j<nrefe_2;j++){

      if(jet80_2 && trgObj_pt_2>=80) hpbpb_TrgObj80[nbins_cent]->Fill(pt_2[j],jet80_p_2);

    }
    
    if(!pHBHENoiseFilter_2 || !pcollisionEventSelection_2) continue;
    
    if(fabs(vz_2)>15) continue;
    
    if(!jet80_2) continue;

    if(nbins_cent==1){
      //cout<<"inside 0-1"<<endl;
      for(int j = 0;j<nrefe_2;j++){
	
	if(fabs(eta_2[j])>2) continue;
	
	if(chMax_2[j]/pt_2[j]<0.01) continue;

	if(neMax_2[j]/TMath::Max(chSum_2[j],neSum_2[j])>=0.975) continue;

	
	if(jet80_2 && trgObj_pt_2>=80){ hpbpb_TrgObj80[nbins_cent-1]->Fill(pt_2[j],jet80_p_2);
	  //fHLT_80<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;
	}
	
	if(TMath::Max(chMax_2[j],neMax_2[j])/(TMath::Max(chSum_2[j],neSum_2[j]))>0.975)
	  fHLT_high<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;

      } // jet loop
      
    }else if(nbins_cent==6){
      
      // checking the centrality class statement
      if(hiBin_2>=5*boundaries_cent[0] && hiBin_2<5*boundaries_cent[1]){
	//cout<<"inside 0-1"<<endl;
	for(int j = 0;j<nrefe_2;j++){
	  
	  if(fabs(eta_2[j])>2) continue;
	  
	  if(chMax_2[j]/pt_2[j]<0.01) continue;

	  if(neMax_2[j]/TMath::Max(chSum_2[j],neSum_2[j])>=0.975) continue;
	
	  if(jet80_2 && trgObj_pt_2>=80){ hpbpb_TrgObj80[0]->Fill(pt_2[j],jet80_p_2);
	    //fHLT_80<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;
	  }
	
	} // jet loop
      
      }
      if(hiBin_2>=5*boundaries_cent[1] && hiBin_2<5*boundaries_cent[2]){
	//cout<<"inside 1-2"<<endl;
	for(int j = 0;j<nrefe_2;j++){
	
	  if(fabs(eta_2[j])>2) continue;
	
	  if(chMax_2[j]/pt_2[j]<0.01) continue;

	  if(neMax_2[j]/TMath::Max(chSum_2[j],neSum_2[j])>=0.975) continue;
	       	
	  if(jet80_2 && trgObj_pt_2>=80){ hpbpb_TrgObj80[1]->Fill(pt_2[j],jet80_p_2);
	    //fHLT_80<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;
	  }
	
	} // jet loop
      
      }
      if(hiBin_2>=5*boundaries_cent[2] && hiBin_2<5*boundaries_cent[3]){
	//cout<<"inside 2-3"<<endl;
	for(int j = 0;j<nrefe_2;j++){
	
	  if(fabs(eta_2[j])>2) continue;
	
	  if(chMax_2[j]/pt_2[j]<0.01) continue;

	  if(neMax_2[j]/TMath::Max(chSum_2[j],neSum_2[j])>=0.975) continue;
	       	
	  if(jet80_2 && trgObj_pt_2>=80){ hpbpb_TrgObj80[2]->Fill(pt_2[j],jet80_p_2);
	    //fHLT_80<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;
	  }
	
	} // jet loop
      
      }
      if(hiBin_2>=5*boundaries_cent[3] && hiBin_2<5*boundaries_cent[4]){
	//cout<<"inside 3-4"<<endl;
	for(int j = 0;j<nrefe_2;j++){
	
	  if(fabs(eta_2[j])>2) continue;
	
	  if(chMax_2[j]/pt_2[j]<0.01) continue;

	  if(neMax_2[j]/TMath::Max(chSum_2[j],neSum_2[j])>=0.975) continue;
	       	
	  if(jet80_2 && trgObj_pt_2>=80){ hpbpb_TrgObj80[3]->Fill(pt_2[j],jet80_p_2);
	    //fHLT_80<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;
	  }
	
	} // jet loop
      
      }
      if(hiBin_2>=5*boundaries_cent[4] && hiBin_2<5*boundaries_cent[5]){
	//cout<<"inside 4-4"<<endl;
	for(int j = 0;j<nrefe_2;j++){
	
	  if(fabs(eta_2[j])>2) continue;
	
	  if(chMax_2[j]/pt_2[j]<0.01) continue;

	  if(neMax_2[j]/TMath::Max(chSum_2[j],neSum_2[j])>=0.975) continue;
	       	
	  if(jet80_2 && trgObj_pt_2>=80){ hpbpb_TrgObj80[4]->Fill(pt_2[j],jet80_p_2);
	    //fHLT_80<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;
	  }
	
	} // jet loop
      
      }
      if(hiBin_2>=5*boundaries_cent[5] && hiBin_2<5*boundaries_cent[6]){
	//cout<<"inside 5-5"<<endl;
	for(int j = 0;j<nrefe_2;j++){
	
	  if(fabs(eta_2[j])>2) continue;
	
	  if(chMax_2[j]/pt_2[j]<0.01) continue;

	  if(neMax_2[j]/TMath::Max(chSum_2[j],neSum_2[j])>=0.975) continue;
	       	
	  if(jet80_2 && trgObj_pt_2>=80){ hpbpb_TrgObj80[5]->Fill(pt_2[j],jet80_p_2);
	    //fHLT_80<<run_2<<" "<<evt_2<<" "<<lumi_2<<" "<<trgObj_pt_2<<" "<<trgObj_eta_2<<" "<<trgObj_phi_2<<" "<<hiBin_2<<" "<<pt_2[j]<<" "<<eta_2[j]<<" "<<phi_2[j]<<endl;
	  }
	
	} // jet loop
      
      }

    }
    
  } // nentries_jet80or95 loop

  fHLT_80.close();
  fHLT_65.close();
  fHLT_55.close();
  fHLT_high.close();

  //add the hpbpb_TrgObjComb[nbins_cent] for the trigger turn on curve here. remember no cuts for this histogram 
  hpbpb_TrgObjComb[nbins_cent]->Add(hpbpb_TrgObj80[nbins_cent]);
  hpbpb_TrgObjComb[nbins_cent]->Add(hpbpb_TrgObj65[nbins_cent]);
  hpbpb_TrgObjComb[nbins_cent]->Add(hpbpb_TrgObj55[nbins_cent]);
  

  for(int i = 0;i<nbins_cent;i++){

    hpbpb_TrgObj80[i]->Scale(1./149.382e6);//respective lumi seen by the trigger all in inverse micro barns 
    hpbpb_TrgObj65[i]->Scale(1./3.195e6);
    hpbpb_TrgObj55[i]->Scale(1./2.734e6);
    
    hpbpb_TrgObj80[i]->Scale(1./4);//delta eta
    hpbpb_TrgObj65[i]->Scale(1./4);
    hpbpb_TrgObj55[i]->Scale(1./4);

    if(i!=nbins_cent){
      hpbpb_TrgObj80[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));//centrality bin width 
      hpbpb_TrgObj65[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));
      hpbpb_TrgObj55[i]->Scale(1./0.005/(5*boundaries_cent[i+1]-5*boundaries_cent[i]));
    }

    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj80[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj65[i]);
    hpbpb_TrgObjComb[i]->Add(hpbpb_TrgObj55[i]);

    hpbpb_TrgObjComb[i] = (TH1F*)hpbpb_TrgObjComb[i]->Rebin(nbins_pt,Form("hpbpb_TrgObjComb_cent%d",i),boundaries_pt);
    hpbpb_TrgObj80[i] = (TH1F*)hpbpb_TrgObj80[i]->Rebin(nbins_pt,Form("hpbpb_TrfObj80_cent%d",i),boundaries_pt);
    hpbpb_TrgObj65[i] = (TH1F*)hpbpb_TrgObj65[i]->Rebin(nbins_pt,Form("hpbpb_TrgObj65_cent%d",i),boundaries_pt);
    hpbpb_TrgObj55[i] = (TH1F*)hpbpb_TrgObj55[i]->Rebin(nbins_pt,Form("hpbpb_TrgObj55_cent%d",i),boundaries_pt);

    divideBinWidth(hpbpbComb[i]);
    divideBinWidth(hpbpb3[i]);
    divideBinWidth(hpbpb2[i]);
    divideBinWidth(hpbpb1[i]);
    
    hpbpb_TrgObjComb[i]->Print("base");
    hpbpb_TrgObj80[i]->Print("base");
    hpbpb_TrgObj65[i]->Print("base");
    hpbpb_TrgObj55[i]->Print("base");

  }


 

  TDatime date;

  //declare the output file
  TFile f(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Output/PbPb_pp_data_ak%d_%s_cent%d_chMax_12003cut.root",radius,algo,date.GetDate()),"RECREATE");
  
  f.cd();
  
  for(int i = 0;i<=nbins_cent;i++){
    hpbpb1[i]->Write();
    hpbpb2[i]->Write();
    hpbpb3[i]->Write();
    hpbpbComb[i]->Write();
    hpbpb_TrgObjComb[i]->Write();
    hpbpb_TrgObj80[i]->Write();
    hpbpb_TrgObj65[i]->Write();
    hpbpb_TrgObj55[i]->Write();
    hpbpb_80[i]->Write();
    hpbpb_65[i]->Write();
    hpbpb_55[i]->Write();
  }
  //hPPComb_bins->Write();
  hpp1->Write();
  hpp2->Write();
  hpp3->Write();
  hppComb->Write();
  
  f.Write();

  f.Close();




}
