// Raghav Kunnawalkam Elayavalli
// July 1st 2014
// CERN

//
// Macro to readin the NLO files and compare it with the pp data we have. 
// 

// I need to decide if this macro will also plot the reuqired histograms or would i have to create a new macro for that
// going by what ive done before i think ill create a new macro. 

// July 11 - we also have the NP corrections here /net/hisrv0001/home/rkunnawa/WORK/NP_factors/test/txtfiles

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

  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};


using namespace std;

static const int dir = 50;

void RAA_nlo_compare(int radius = 3, int energy = 2760){

	TStopwatch timer;
	timer.Start();

	TDatime date;

	gStyle->SetOptStat(0);

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();

	cout<<"Starting the comparison of the NLO with data"<<endl;
	cout<<"Running for Energy = "<<energy<<" and Radius = "<<radius<<endl;

	//nlo files: right now they are not for all the eta ranges. they are only present 
	TFile *fNLO_nnpdf = TFile::Open("nlo_files/input_rootfiles/fnl4350a_nnpdf21-nlo_aspdf_new.root");
  	TFile *fNLO_cteq = TFile::Open("nlo_files/input_rootfiles/fnl4350a_cteq66-nlo_aspdf_all_new.root");
  	TFile *fNLO_ct10n = TFile::Open("nlo_files/input_rootfiles/fnl4350a_ct10n-nlo_aspdf_new.root");
  	TFile *fNLO_hera = TFile::Open("nlo_files/input_rootfiles/fnl4350a_hera15all-nlo_aspdf_new.root");

  	//npc files: 
  	// for now only take the central eta file for the different radii. then we can take all the rest. 

	char dirName[dir][256] = {
    "ak2GenJetSpectrum_n10_p10","ak2GenJetSpectrum_n20_p20","ak2GenJetSpectrum_n25_n20","ak2GenJetSpectrum_n20_n15",
    "ak2GenJetSpectrum_n15_n10","ak2GenJetSpectrum_n10_n05","ak2GenJetSpectrum_n05_p05","ak2GenJetSpectrum_p05_p10",
    "ak2GenJetSpectrum_p10_p15","ak2GenJetSpectrum_p15_p20",
    "ak3GenJetSpectrum_n10_p10","ak3GenJetSpectrum_n20_p20","ak3GenJetSpectrum_n25_n20","ak3GenJetSpectrum_n20_n15",
    "ak3GenJetSpectrum_n15_n10","ak3GenJetSpectrum_n10_n05","ak3GenJetSpectrum_n05_p05","ak3GenJetSpectrum_p05_p10",
    "ak3GenJetSpectrum_p10_p15","ak3GenJetSpectrum_p15_p20",
    "ak4GenJetSpectrum_n10_p10","ak4GenJetSpectrum_n20_p20","ak4GenJetSpectrum_n25_n20","ak4GenJetSpectrum_n20_n15",
    "ak4GenJetSpectrum_n15_n10","ak4GenJetSpectrum_n10_n05","ak4GenJetSpectrum_n05_p05","ak4GenJetSpectrum_p05_p10",
    "ak4GenJetSpectrum_p10_p15","ak4GenJetSpectrum_p15_p20",
    "ak5GenJetSpectrum_n10_p10","ak5GenJetSpectrum_n20_p20","ak5GenJetSpectrum_n25_n20","ak5GenJetSpectrum_n20_n15",
    "ak5GenJetSpectrum_n15_n10","ak5GenJetSpectrum_n10_n05","ak5GenJetSpectrum_n05_p05","ak5GenJetSpectrum_p05_p10",
    "ak5GenJetSpectrum_p10_p15","ak5GenJetSpectrum_p15_p20",
    "ak7GenJetSpectrum_n10_p10","ak7GenJetSpectrum_n20_p20","ak7GenJetSpectrum_n25_n20","ak7GenJetSpectrum_n20_n15",
    "ak7GenJetSpectrum_n15_n10","ak7GenJetSpectrum_n10_n05","ak7GenJetSpectrum_n05_p05","ak7GenJetSpectrum_p05_p10",
    "ak7GenJetSpectrum_p10_p15","ak7GenJetSpectrum_p15_p20"};

    char etaWidth[dir][256] = {"n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
	"n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
	"p10_eta_p15","p15_eta_p20",
	"n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
	"n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
	"p10_eta_p15","p15_eta_p20",
	"n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
	"n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
	"p10_eta_p15","p15_eta_p20",
	"n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
	"n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
	"p10_eta_p15","p15_eta_p20",
	"n10_eta_p10","n20_eta_p20","n25_eta_n20","n20_eta_n15",
	"n15_eta_n10","n10_eta_n05","n05_eta_p05","p05_eta_p10",
	"p10_eta_p15","p15_eta_p20"
	};

	char radius_lable[dir][256] = {
	"R2","R2","R2","R2","R2","R2","R2","R2","R2","R2",
	"R3","R3","R3","R3","R3","R3","R3","R3","R3","R3",
	"R4","R4","R4","R4","R4","R4","R4","R4","R4","R4",
	"R5","R5","R5","R5","R5","R5","R5","R5","R5","R5",
	"R7","R7","R7","R7","R7","R7","R7","R7","R7","R7"};

	ifstream fin_txt[dir];

	for(int i = 0;i<dir;i++){

		ostringstream filename;
		filename<<"/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Macros/RAA/nlo_files/input_np_txtfiles/NPC_ak_"<<radius_lable[i]<<etaWidth[i]<<"_energy"<<energy<<".txt";
		fin_txt[i].open(filename.str());

	}



}









































