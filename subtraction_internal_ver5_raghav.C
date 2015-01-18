// Raghav Kunnawalkam Elayavalli
// Dec 2nd 2014
// Rutgers

//
// Macro to study the UE tables and plot the RMS value of the Sum$(PFVsPtInitial) as a function of the eta. 
// based on Yue Shi's macros ver0 to ver4. 
// 

// Dec 4th - checking for the divergences in the polynomials compared to the different flow candidates. 
//         - the values for the divergences are given here:  https://twiki.cern.ch/twiki/pub/CMS/HiVsValidation/4Dec2014_VsValidation_divergence.pdf 
//         - And these are compared with Yue Shi's macro subtraction_internal_v4.cc which plots the polynomials 
//         - these above plots were made from the training code. we are going to check the data and the MC for such a difference. 
//         - this is the x axis for cos:  vn = sumpt[0]*vn[n][0]*cos(n*Psi_[n][0]) + sumpt[14]*vn[n][14]*cos(n*Psi_[n][14])

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TMath.h>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "Headers/plot.h"

#define NOBJECT_MAX 16384

double ue_predictor_pf[3][15][5][2][82];
// value information are given here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/VoronoiUEList#Useful_information_on_the_tables
double ue_interpolation_pf0[15][344];
// still not sure what these 344 are? 
double ue_interpolation_pf1[15][344];
double ue_interpolation_pf2[15][82];
// for 82 look at the previous link. 

size_t pf_id_reduce(const Int_t pf_id)
{
  // Particle::pdgId_ PFCandidate::particleId_
  // PFCandidate::ParticleType Particle
  // 0           0  X          unknown, or dummy 
  // +211, -211  1  h          charged hadron 
  // +11, -11    2  e          electron 
  // +13, -13    3  mu         muon 
  // 22          4  gamma      photon 
  // 130         5  h0         neutral hadron 
  // 130         6  h_HF       hadronic energy in an HF tower 
  // 22          7  egamma_HF  electromagnetic energy in an HF tower

  if (pf_id == 4) {
    return 1;
  }
  else if (pf_id >= 5 && pf_id <= 7) {
    return 2;
  }

  return 0;
}

void subtraction_internal_ver5_raghav(const char *filename = "/Users/keraghav/WORK/RAA/Output/PbPb_data_bad_events.root", const int data = 1, const int calorimetric = 0){

  gStyle->SetPalette(55);
  bool printDebug = false;
  bool printDebug2 = false;
  
  static const size_t nfourier = 5;

  const char *root_tree_name = calorimetric ?
    "rechitanalyzer/tower" : "pfcandAnalyzer/pfTree";

  TFile *root_file = TFile::Open(filename);
  TTree *root_tree = dynamic_cast<TTree *>(gDirectory->Get(root_tree_name));

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

  if (calorimetric) {
    root_tree->SetBranchAddress("n", &nPFpart);
    root_tree->SetBranchAddress("et", pfPt);
    root_tree->SetBranchAddress("eta", pfEta);
    root_tree->SetBranchAddress("phi", pfPhi);
    root_tree->SetBranchAddress("vsArea", pfArea);
  }
  else {
    root_tree->SetBranchAddress("nPFpart", &nPFpart);
    root_tree->SetBranchAddress("pfId", pfId);
    root_tree->SetBranchAddress("pfPt", pfPt);
    root_tree->SetBranchAddress("pfVsPtInitial", pfVsPtInitial);
    root_tree->SetBranchAddress("pfVsPt", pfVsPt);
    root_tree->SetBranchAddress("pfEta", pfEta);
    root_tree->SetBranchAddress("pfPhi", pfPhi);
    root_tree->SetBranchAddress("pfArea", pfArea);
    root_tree->SetBranchAddress("vn",&v_n);
    root_tree->SetBranchAddress("psin",&psi_n);
    root_tree->SetBranchAddress("sumpt",&sumpT);
  }

  TTree *hiTree = dynamic_cast<TTree *>(gDirectory->Get("hiEvtAnalyzer/HiTree"));

  Int_t run;
  Int_t lumi;
  Float_t vx;
  Float_t vy;
  Float_t vz;
  Int_t hiBin;
  Float_t hiHF;
  Float_t hiHFplus;
  Float_t hiHFminus;
  Float_t hiZDC;
  Float_t hiZDCplus;
  Float_t hiZDCminus;
  Float_t hiHFhit;
  Float_t hiHFhitPlus;
  Float_t hiHFhitMinus;
  Float_t hiET;
  Float_t hiEE;
  Float_t hiEB;
  Float_t hiEEplus;
  Float_t hiEEminus;
  Int_t hiNpix;
  Int_t hiNpixelTracks;
  Int_t hiNtracks;
  Int_t hiNtracksPtCut;
  Int_t hiNtracksEtaCut;
  Int_t hiNtracksEtaPtCut;
  Int_t hiNevtPlane;
  Float_t hiEvtPlanes[38];

  // Set branch addresses.
  hiTree->SetBranchAddress("run",&run);
  hiTree->SetBranchAddress("lumi",&lumi);
  hiTree->SetBranchAddress("vx",&vx);
  hiTree->SetBranchAddress("vy",&vy);
  hiTree->SetBranchAddress("vz",&vz);
  hiTree->SetBranchAddress("hiBin",&hiBin);
  hiTree->SetBranchAddress("hiHF",&hiHF);
  hiTree->SetBranchAddress("hiHFplus",&hiHFplus);
  hiTree->SetBranchAddress("hiHFminus",&hiHFminus);
  hiTree->SetBranchAddress("hiZDC",&hiZDC);
  hiTree->SetBranchAddress("hiZDCplus",&hiZDCplus);
  hiTree->SetBranchAddress("hiZDCminus",&hiZDCminus);
  hiTree->SetBranchAddress("hiHFhit",&hiHFhit);
  hiTree->SetBranchAddress("hiHFhitPlus",&hiHFhitPlus);
  hiTree->SetBranchAddress("hiHFhitMinus",&hiHFhitMinus);
  hiTree->SetBranchAddress("hiET",&hiET);
  hiTree->SetBranchAddress("hiEE",&hiEE);
  hiTree->SetBranchAddress("hiEB",&hiEB);
  hiTree->SetBranchAddress("hiEEplus",&hiEEplus);
  hiTree->SetBranchAddress("hiEEminus",&hiEEminus);
  hiTree->SetBranchAddress("hiNpix",&hiNpix);
  hiTree->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks);
  hiTree->SetBranchAddress("hiNtracks",&hiNtracks);
  hiTree->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut);
  hiTree->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut);
  hiTree->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut);
  hiTree->SetBranchAddress("hiNevtPlane",&hiNevtPlane);
  hiTree->SetBranchAddress("hiEvtPlanes",hiEvtPlanes);

  const char *tree_name = "akVs3PFJetAnalyzer/t";

  TTree *t = dynamic_cast<TTree *>(gDirectory->Get(tree_name));

  Int_t evt;
  // Float_t b;
  Int_t nref;
  Float_t rawpt[NOBJECT_MAX];
  Float_t jtpt[NOBJECT_MAX];
  Float_t jteta[NOBJECT_MAX];
  Float_t jty[NOBJECT_MAX];
  Float_t jtphi[NOBJECT_MAX];
  Float_t jtpu[NOBJECT_MAX];

  // Set branch addresses.
  t->SetBranchAddress("evt", &evt);
  // t->SetBranchAddress("b", &b);
  t->SetBranchAddress("nref", &nref);
  t->SetBranchAddress("rawpt", rawpt);
  t->SetBranchAddress("jtpt", jtpt);
  t->SetBranchAddress("jteta", jteta);
  t->SetBranchAddress("jty", jty);
  t->SetBranchAddress("jtphi", jtphi);
  t->SetBranchAddress("jtpu", jtpu);

  std::ifstream in_stream(data ? (calorimetric ? "ue_calibrations_calo_data.txt" : "ue_calibrations_pf_data.txt") : (calorimetric ? "ue_calibrations_calo_mc.txt" : "ue_calibrations_pf_mc.txt"));
  std::string line;
  size_t index = 0;
  const size_t nline_predictor = 3 * 15 * (1 + (5 - 1) * 2) * 82;

  while (std::getline(in_stream, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream line_stream(line);
    double val;
    int bin0, bin1, bin2, bin3, bin4;

    if (index < nline_predictor) {
      line_stream >> bin0 >> bin1 >> bin2 >> bin3 >> bin4 >> val;
      ue_predictor_pf[bin0][bin1][bin2][bin3][bin4] = val;
    }
    else if (index < nline_predictor + sizeof(ue_interpolation_pf0) / sizeof(double)) {
      line_stream >> bin0 >> bin1 >> val;
      ue_interpolation_pf0[bin0][bin1] = val;
    }
    else if (index < nline_predictor + sizeof(ue_interpolation_pf0) / sizeof(double) + sizeof(ue_interpolation_pf1) / sizeof(double)) {
      line_stream >> bin0 >> bin1 >> val;
      ue_interpolation_pf1[bin0][bin1] = val;
    }
    else if (index < nline_predictor + sizeof(ue_interpolation_pf0) / sizeof(double) + sizeof(ue_interpolation_pf1) / sizeof(double) + sizeof(ue_interpolation_pf2) / sizeof(double)) {
      line_stream >> bin0 >> bin1 >> val;
      ue_interpolation_pf2[bin0][bin1] = val;
    }
    index++;
  }

  static const size_t nreduced_id = 3;

  static const size_t nedge_pseudorapidity = 15 + 1;
  static const double edge_pseudorapidity[nedge_pseudorapidity] = {-5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522, 0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191 };

  // there are the bins which we might end up moving to! 
  // static const size_t nedge_pseudorapidity = 19 + 1;
  // static const double edge_pseudorapidity[nedge_pseudorapidity] = {-5.191, -2.964, -2.853 ,-2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522, 0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 2.853, 2.964,5.191 };



  const std::vector<double> edge_pseudorapidity_v(edge_pseudorapidity, edge_pseudorapidity + nedge_pseudorapidity);

  static const size_t ncms_hcal_edge_pseudorapidity = 82 + 1;
  static const double cms_hcal_edge_pseudorapidity[
						   ncms_hcal_edge_pseudorapidity] = {
    -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013,
    -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853,
    -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830,
    -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
    -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609,
    -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
    0.000,
    0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,
    0.696,  0.783,  0.879,  0.957,  1.044,  1.131,  1.218,
    1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,
    1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  2.853,
    2.964,  3.139,  3.314,  3.489,  3.664,  3.839,  4.013,
    4.191,  4.363,  4.538,  4.716,  4.889,  5.191
  };

  const std::vector<double> cms_hcal_edge_pseudorapidity_v(cms_hcal_edge_pseudorapidity, cms_hcal_edge_pseudorapidity + ncms_hcal_edge_pseudorapidity);

  static const size_t ncms_ecal_edge_pseudorapidity = 344 + 1;
  double cms_ecal_edge_pseudorapidity[
				      ncms_ecal_edge_pseudorapidity];

  for (size_t i = 0; i < ncms_ecal_edge_pseudorapidity; i++) {
    cms_ecal_edge_pseudorapidity[i] =
      i * (2 * 2.9928 /
	   (ncms_ecal_edge_pseudorapidity - 1)) -
      2.9928;
  };

  const std::vector<double> cms_ecal_edge_pseudorapidity_v(cms_ecal_edge_pseudorapidity, cms_ecal_edge_pseudorapidity + ncms_ecal_edge_pseudorapidity);

  // lets do a simple check on the ue_predictor_pf to see where v0 (bulk) and v2 blow up. 
 
  for(int i = 0;i<3;i++){
    for(int j = 0;j<15;j++){
      for(int k = 0;k<2;k++){
	for(int l = 0;l<82;l++){
	  //std::cout<<ue_predictor_pf[i][j][0][k][l]<<std::endl;
	  /*
	  if(TMath::Abs(ue_predictor_pf[i][j][0][k][l])>1e6){
	    if(printDebug)std::cout<<i<<" "<<j<<" 0 "<<k<<" "<<l<<" = "<<ue_predictor_pf[i][j][0][k][l]<<std::endl; 
	  }
	  
	  if(TMath::Abs(ue_predictor_pf[i][j][1][k][l])>1e6){
	    if(printDebug)std::cout<<i<<" "<<j<<" 1 "<<k<<" "<<l<<" = "<<ue_predictor_pf[i][j][1][k][l]<<std::endl; 
	  }
	  
	  if(TMath::Abs(ue_predictor_pf[i][j][2][k][l])>1e6){
	    if(printDebug)std::cout<<i<<" "<<j<<" 2 "<<k<<" "<<l<<" = "<<ue_predictor_pf[i][j][2][k][l]<<std::endl; 
	  }
	  
	  if(TMath::Abs(ue_predictor_pf[i][j][3][k][l])>1e6){
	    if(printDebug)std::cout<<i<<" "<<j<<" 3 "<<k<<" "<<l<<" = "<<ue_predictor_pf[i][j][3][k][l]<<std::endl; 
	  }
	  
	  if(TMath::Abs(ue_predictor_pf[i][j][4][k][l])>1e6){
	    if(printDebug)std::cout<<i<<" "<<j<<" 4 "<<k<<" "<<l<<" = "<<ue_predictor_pf[i][j][4][k][l]<<std::endl; 
	  }
	  */
	}
      }
    }
  }

  size_t nentries = root_tree->GetEntries();
  std::cout<<"No of Entries = "<<nentries<<endl;
  // nentries = 2;

  for(size_t i = 0; i<nentries; i++){

    root_tree->GetEntry(i);
    hiTree->GetEntry(i);
    t->GetEntry(i);
    
    
    // Event collective Fourier components, per particle flow ID
    // group. Note that since by heavy ion convention, dN/dphi =
    // v0/(2 pi) + v0 v1/pi cos(phi - Psi_RP) + v0 v2/pi cos(2(phi
    // - Psi_RP)) + ..., and orthonormal relation for the Fourier
    // basis f0 = v0, f1 = v0 v1, ..., if f is the Fourier and v0
    // the phi-averaged 1/(2 pi) dpT/dy

    double perp_fourier[nedge_pseudorapidity - 1][nreduced_id][nfourier][2];

    for (size_t k = 1; k < nedge_pseudorapidity; k++) {
      for (size_t l = 0; l < nreduced_id; l++) {
	for (size_t m = 0; m < nfourier; m++) {
	  for (size_t re_or_im = 0; re_or_im < 2;
	       re_or_im++) {
	    perp_fourier[k - 1][l][m][re_or_im] = 0;
	  }
	}
      }
    }

    memset(perp_fourier, 0,
	   (nedge_pseudorapidity - 1) * nreduced_id * nfourier * 2 *
	   sizeof(double));

    for (Int_t j = 0; j < nPFpart; j++) {
      size_t reduced_id = pf_id_reduce(pfId[j]);

      for (size_t k = 1; k < nedge_pseudorapidity; k++) {
	if (pfEta[j] >= edge_pseudorapidity[k - 1] &&
	    pfEta[j] < edge_pseudorapidity[k]) {
	  for (size_t l = 0; l < nfourier; l++) {
	    perp_fourier[k - 1][reduced_id][l][0] +=
	      pfPt[j] * cos(l * pfPhi[j]);
	    perp_fourier[k - 1][reduced_id][l][1] +=
	      pfPt[j] * sin(l * pfPhi[j]);
	  }
	}
      }
    }

    
    // Event selection
  
    static const size_t nfeature = 2 * nfourier - 1;
    double feature[nfeature];

    // Scale factor to get 95% of the coefficient below 1.0

    std::vector<double> scale(nfourier, 1.0 / 200.0);

    if (nfourier >= 1) {
      scale[0] = 1.0 / 5400.0;
    }
    if (nfourier >= 2) {
      scale[1] = 1.0 / 130.0;
    }
    if (nfourier >= 3) {
      scale[2] = 1.0 / 220.0;
    }

    feature[0] = 0;
    for (size_t l = 0; l < nreduced_id; l++) {
      feature[0] += scale[0] *
	(perp_fourier[0                       ][l][0][0] +
	 perp_fourier[nedge_pseudorapidity - 2][l][0][0]);
    }
    if (i == 6 || i == 11 || i == 12  || i == 14)
      //fprintf(stderr, "%s:%d: HF bulk: %f scaled, %f GeV/c unscaled\n", __FILE__, __LINE__, feature[0], feature[0] / scale[0]);
    for (size_t k = 1; k < nfourier; k++) {
      feature[2 * k - 1] = 0;
      for (size_t l = 0; l < nreduced_id; l++) {
	feature[2 * k - 1] += scale[k] *
	  (perp_fourier[0                       ][l][k][0] +
	   perp_fourier[nedge_pseudorapidity - 2][l][k][0]);
      }
      if (i == 6 || i == 11 || i == 12  || i == 14)
	//fprintf(stderr, "%s:%d: HF order %lu flow: cos %f scaled, %f GeV/c unsacaled\n", __FILE__, __LINE__, k, feature[2 * k - 1], feature[2 * k - 1] / scale[k]);
      feature[2 * k] = 0;
      for (size_t l = 0; l < nreduced_id; l++) {
	feature[2 * k] += scale[k] *
	  (perp_fourier[0                       ][l][k][1] +
	   perp_fourier[nedge_pseudorapidity - 2][l][k][1]);
      }
      if (i == 6 || i == 11 || i == 12  || i == 14)
	//fprintf(stderr, "%s:%d: HF order %lu flow: sin %f scaled, %f GeV/c unsacaled\n", __FILE__, __LINE__, k, feature[2 * k], feature[2 * k] / scale[k]);
    }

#if 0
    const double event_plane = atan2(feature[4], feature[3]);
    const double v2 =
      sqrt(feature[3] * feature[3] +
	   feature[4] * feature[4]) / feature[0];
#endif

    const double max_feature = std::max(-(*std::min_element(feature, feature + nfeature)), *std::max_element(feature, feature + nfeature));

    if (i == 6 || i == 11 || i == 12  || i == 14)
      fprintf(stderr, "%s:%d: %f %f\n", __FILE__, __LINE__, hiBin * 0.5, sqrt(feature[3] * feature[3] + feature[4] * feature[4]));

    for (Int_t k = 0; k < nPFpart; k++) {
      int predictor_index = -1;
      int interpolation_index = -1;
      double density = 0;

      if (pfEta[k] >= edge_pseudorapidity[0] &&
	  pfEta[k] < edge_pseudorapidity[nedge_pseudorapidity - 1]) {
	std::vector<double>::const_iterator p = std::lower_bound(edge_pseudorapidity_v.begin(), edge_pseudorapidity_v.end(), pfEta[k]);

	predictor_index = (p - edge_pseudorapidity_v.begin()) - 1;
      }

      for (size_t j = 0; j < nreduced_id; j++) {
	if (j == 2) {
	  // HCAL
	  if (pfEta[k] >=
	      cms_hcal_edge_pseudorapidity[0] &&
	      pfEta[k] <
	      cms_hcal_edge_pseudorapidity[ncms_hcal_edge_pseudorapidity - 1]) {
	    std::vector<double>::const_iterator p = std::lower_bound(cms_hcal_edge_pseudorapidity_v.begin(), cms_hcal_edge_pseudorapidity_v.end(), pfEta[k]);

	    interpolation_index = (p - cms_hcal_edge_pseudorapidity_v.begin()) - 1;
	  }
	}
	else {
	  // Tracks or ECAL clusters
	  if (pfEta[k] >=
	      cms_ecal_edge_pseudorapidity[0] &&
	      pfEta[k] <
	      cms_ecal_edge_pseudorapidity[ncms_ecal_edge_pseudorapidity - 1]) {
	    std::vector<double>::const_iterator p = std::lower_bound(cms_ecal_edge_pseudorapidity_v.begin(), cms_ecal_edge_pseudorapidity_v.end(), pfEta[k]);

	    interpolation_index = (p - cms_ecal_edge_pseudorapidity_v.begin()) - 1;
	  }
	}
	
	
	if (predictor_index >= 0 && interpolation_index >= 0) {
	  // Calculate the aggregated prediction and
	  // interpolation for the pseudorapidity segment

	  const double azimuth = pfPhi[k];
	  const double (*p)[2][82] =
	    ue_predictor_pf[j][predictor_index];
	  double pred = 0;

	  for (size_t l = 0; l < nfourier; l++) {
	    for (size_t m = 0; m < 2; m++) {
	      float u = p[l][m][0];

	      for (size_t n = 0; n < 2 * nfourier - 1; n++) {
		//if (max_feature < 0.9 || (l == 0 && n == 0) || (l == 2 && (n == 3 || n == 4))) {
		if (true || (l == 0 && n == 0)) {
		  u += (((((((((p[l][m][9 * n + 9]) *
			       feature[n] +
			       p[l][m][9 * n + 8]) *
			      feature[n] +
			      p[l][m][9 * n + 7]) *
			     feature[n] +
			     p[l][m][9 * n + 6]) *
			    feature[n] +
			    p[l][m][9 * n + 5]) *
			   feature[n] +
			   p[l][m][9 * n + 4]) *
			  feature[n] +
			  p[l][m][9 * n + 3]) *
			 feature[n] +
			 p[l][m][9 * n + 2]) *
			feature[n] +
			p[l][m][9 * n + 1]) *
		    feature[n];
		}
		//else if (n == 2 * l - 1 || n == 2 * l) {
		else if ((l == 2 && (n == 3 || n == 4))) {
		  u += p[l][m][9 * n + 1] * feature[n];
		}
	      }
	    
#if 0
	      // This looks at a specific flow component and see how the polynomial is evaluated
	      if (j == 0 && predictor_index == 3 && l == 0 && m == 0) {
		// if(printDebug)fprintf(stderr, "%s:%d: %f %f\n", __FILE__, __LINE__, perp_fourier[0][2][2][0], perp_fourier[nedge_pseudorapidity - 2][2][2][1]);
		if(printDebug)fprintf(stderr, "%s:%d: << %f %f %f %f %f %f %f\n", __FILE__, __LINE__, feature[0], feature[1], feature[2], feature[3], feature[4], u, perp_fourier[predictor_index][j][l][m]);
	      }
#endif
	    
	      pred += u * (l == 0 ? 1.0 : 2.0) *
		(m == 0 ? cos(l * azimuth) :
		 sin(l * azimuth));
	      

	    }
	  }
     
	  //if(printDebug)std::cout<<"pred = "<<pred<<std::endl;

	  double interp;

	  if (j == 0) {
	    interp =
	      ue_interpolation_pf0[predictor_index][
						    interpolation_index];
	  }
	  else if (j == 1) {
	    interp =
	      ue_interpolation_pf1[predictor_index][
						    interpolation_index];
	  }
	  else if (j == 2) {
	    interp =
	      ue_interpolation_pf2[predictor_index][
						    interpolation_index];
	  }

	  // Interpolate down to the finely binned
	  // pseudorapidity

	  density += pred /
	    (2.0 * M_PI *
	     (edge_pseudorapidity[predictor_index + 1] -
	      edge_pseudorapidity[predictor_index])) *
	    interp;

	}
      }
    
      // Prints the subtracted density * area
      // fprintf(stderr, "%s:%d: %.8e %.8e %.8e %.8e %.8e\n", __FILE__, __LINE__, hiBin * 0.5, pfEta[k], pfPhi[k], pfPt[k], density * pfArea[k]);
      // root_histogram0.Fill(pfEta[k], pfPhi[k], pfPt[k] - density * pfArea[k]);

      if(printDebug && TMath::Abs(density * pfArea[k])>10000 && TMath::Abs(pfEta[k])<2) 
	std::cout<<"eta = "<<pfEta[k]<<"; phi = "<<pfPhi[k]<<"; density*pfArea[k] = "<<density * pfArea[k]<<std::endl;
    
    }
  
  }
  

  // declare the histograms we want here: 
  // Sum$(PfVsPtInitial) as a function of the 15 eta bins. -> 15 histograms. since i want the RMS value for each of them. 
  // same for PFVsPt 
  // want to see when bulk = v0+v2 over v0 diverges. or when v2 diverges. or what is the contribution of v3 and v4 for those eevents.   
  TH1F *hSumPfVsPtInitial[nedge_pseudorapidity-1], *hSumPfVsPt[nedge_pseudorapidity-1], *hSumPfPt[nedge_pseudorapidity-1];
  TH2F *hPfVsPtInitial_eta, *hPfVsPt_eta, *hPfPt_eta;
  TH2F *hRMSSumPFVsPtIni_eta;
  TH2F *hRMSSumPFVsPt_eta;
  TH2F *hRMSSumPFPt_eta;
  TH1F *hetabin_test;

  for(int j = 0;j<nedge_pseudorapidity-1;j++){

    hSumPfVsPtInitial[j] = new TH1F(Form("hSumPfVsPtInitial_%d",j),"",10000,-10000,10000);
    hSumPfVsPt[j] = new TH1F(Form("hSumPfVsPt_%d",j),"",10000,-10000,10000);
    hSumPfPt[j] = new TH1F(Form("hSumPfPt_%d",j),"",10000,-10000,10000);
    
  }

  hPfVsPtInitial_eta = new TH2F("hPfVsPtInitial_eta","",10000,-10000,10000,60,-6,+6);
  hPfVsPt_eta = new TH2F("hPfVsPt_eta","",10000,-10000,10000,60,-6,+6);
  hPfPt_eta = new TH2F("hPfPt_eta","",10000,-10000,10000,60,-6,+6);
  hRMSSumPFVsPtIni_eta = new TH2F("hRMSSumPFVsPtIni_eta","",nedge_pseudorapidity-1,edge_pseudorapidity,10000,-1000,3000);
  hRMSSumPFVsPt_eta = new TH2F("hRMSSumPFVsPt_eta","",nedge_pseudorapidity-1,edge_pseudorapidity,10000,-1000,3000);
  hRMSSumPFPt_eta = new TH2F("hRMSSumPFPt_eta","",nedge_pseudorapidity-1,edge_pseudorapidity,10000,-1000,3000);
  hetabin_test = new TH1F("hetabin_test","",nedge_pseudorapidity-1,edge_pseudorapidity);

  if(printDebug)nentries = 2;
  //start the event loop
  for(size_t i = 0; i<nentries; i++){

    root_tree->GetEntry(i);
    hiTree->GetEntry(i);
    t->GetEntry(i);

    Float_t Vs_0_x_minus = sumpT[0]*v_n[0][0]*TMath::Cos(0*psi_n[0][0]);
    Float_t Vs_0_x_plus = sumpT[14]*v_n[0][14]*TMath::Cos(0*psi_n[0][14]);
    Float_t Vs_0_y_minus = sumpT[0]*v_n[0][0]*TMath::Sin(0*psi_n[0][0]);
    Float_t Vs_0_y_plus = sumpT[14]*v_n[0][14]*TMath::Sin(0*psi_n[0][14]);
    Float_t Vs_0_x = Vs_0_x_minus + Vs_0_x_plus;
    Float_t Vs_0_y = Vs_0_y_minus + Vs_0_y_plus;
    if(printDebug)std::cout<<"Vs_0_x = "<<Vs_0_x<<"; Vs_0_y =  "<<Vs_0_y<<std::endl;
    if(TMath::Abs(Vs_0_x)>5200) std::cout<<"event "<<i<<" has sumpT* cos v0 > 5200"<<std::endl;
    if(TMath::Abs(Vs_0_y)>5200) std::cout<<"event "<<i<<" has sumpT* sin v0 > 5200"<<std::endl;

    Float_t Vs_1_x_minus = sumpT[0]*v_n[1][0]*TMath::Cos(2*psi_n[1][0]);
    Float_t Vs_1_x_plus = sumpT[14]*v_n[1][14]*TMath::Cos(2*psi_n[1][14]);
    Float_t Vs_1_y_minus = sumpT[0]*v_n[1][0]*TMath::Sin(2*psi_n[1][0]);
    Float_t Vs_1_y_plus = sumpT[14]*v_n[1][14]*TMath::Sin(2*psi_n[1][14]);
    Float_t Vs_1_x = Vs_1_x_minus + Vs_1_x_plus;
    Float_t Vs_1_y = Vs_1_y_minus + Vs_1_y_plus;
    if(printDebug)std::cout<<"Vs_1_x = "<<Vs_1_x<<"; Vs_1_y =  "<<Vs_1_y<<std::endl;
    if(TMath::Abs(Vs_1_x)>100) std::cout<<"event "<<i<<" has sumpT* cos v1 > 100"<<std::endl;
    if(TMath::Abs(Vs_1_y)>100) std::cout<<"event "<<i<<" has sumpT* sin v1 > 100"<<std::endl;

    Float_t Vs_2_x_minus = sumpT[0]*v_n[2][0]*TMath::Cos(2*psi_n[2][0]);
    Float_t Vs_2_x_plus = sumpT[14]*v_n[2][14]*TMath::Cos(2*psi_n[2][14]);
    Float_t Vs_2_y_minus = sumpT[0]*v_n[2][0]*TMath::Sin(2*psi_n[2][0]);
    Float_t Vs_2_y_plus = sumpT[14]*v_n[2][14]*TMath::Sin(2*psi_n[2][14]);
    Float_t Vs_2_x = Vs_2_x_minus + Vs_2_x_plus;
    Float_t Vs_2_y = Vs_2_y_minus + Vs_2_y_plus;
    if(printDebug)std::cout<<"Vs_2_x = "<<Vs_2_x<<"; Vs_2_y =  "<<Vs_2_y<<std::endl;
    if(TMath::Abs(Vs_2_x)>140) std::cout<<"event "<<i<<" has sumpT* cos v2 > 140"<<std::endl;
    if(TMath::Abs(Vs_2_y)>140) std::cout<<"event "<<i<<" has sumpT* sin v2 > 140"<<std::endl;

    Float_t Vs_3_x_minus = sumpT[0]*v_n[3][0]*TMath::Cos(2*psi_n[3][0]);
    Float_t Vs_3_x_plus = sumpT[14]*v_n[3][14]*TMath::Cos(2*psi_n[3][14]);
    Float_t Vs_3_y_minus = sumpT[0]*v_n[3][0]*TMath::Sin(2*psi_n[3][0]);
    Float_t Vs_3_y_plus = sumpT[14]*v_n[3][14]*TMath::Sin(2*psi_n[3][14]);
    Float_t Vs_3_x = Vs_3_x_minus + Vs_3_x_plus;
    Float_t Vs_3_y = Vs_3_y_minus + Vs_3_y_plus;
    if(printDebug)std::cout<<"Vs_3_x = "<<Vs_3_x<<"; Vs_3_y =  "<<Vs_3_y<<std::endl;
    if(TMath::Abs(Vs_3_x)>120) std::cout<<"event "<<i<<" has sumpT* cos v3 > 120"<<std::endl;
    if(TMath::Abs(Vs_3_y)>120) std::cout<<"event "<<i<<" has sumpT* sin v3 > 120"<<std::endl;

    Float_t Vs_4_x_minus = sumpT[0]*v_n[4][0]*TMath::Cos(2*psi_n[4][0]);
    Float_t Vs_4_x_plus = sumpT[14]*v_n[4][14]*TMath::Cos(2*psi_n[4][14]);
    Float_t Vs_4_y_minus = sumpT[0]*v_n[4][0]*TMath::Sin(2*psi_n[4][0]);
    Float_t Vs_4_y_plus = sumpT[14]*v_n[4][14]*TMath::Sin(2*psi_n[4][14]);
    Float_t Vs_4_x = Vs_4_x_minus + Vs_4_x_plus;
    Float_t Vs_4_y = Vs_4_y_minus + Vs_4_y_plus;
    if(printDebug)std::cout<<"Vs_4_x = "<<Vs_4_x<<"; Vs_4_y =  "<<Vs_4_y<<std::endl;
    if(TMath::Abs(Vs_3_x)>140) std::cout<<"event "<<i<<" has sumpT* cos v4 > 140"<<std::endl;
    if(TMath::Abs(Vs_3_y)>140) std::cout<<"event "<<i<<" has sumpT* sin v4 > 140"<<std::endl;

    // declare the variables to hold the values per event per eta bin. 
    Float_t SumPfVsPtInitial[nedge_pseudorapidity-1], SumPfVsPt[nedge_pseudorapidity-1], SumPfPt[nedge_pseudorapidity-1];
    for(int k = 0;k<nedge_pseudorapidity-1;k++) {SumPfVsPtInitial[k] = 0; SumPfVsPt[k] = 0; SumPfPt[k] = 0;}
    if(printDebug2)std::cout<<"entry = "<< i<<endl;

    for(int ij = 0; ij < nPFpart; ij++){
      if(printDebug2)std::cout<<"pfcand = "<<ij<<endl;

      for(int k = 0;k<nedge_pseudorapidity-1;k++){
	if(printDebug2)std::cout<<"eta bin = "<<k<<endl;
	if(printDebug2)std::cout<<"pfeta = "<<pfEta[ij]<<endl;
	if(printDebug2)std::cout<<"eta boundaries = "<<edge_pseudorapidity[k]<<" "<<edge_pseudorapidity[k+1]<<endl;
	if(pfEta[ij]>= edge_pseudorapidity[k] && pfEta[ij]< edge_pseudorapidity[k+1]){
	  if(printDebug2)std::cout<<"inside eta bin"<<endl;
	  hetabin_test->Fill(edge_pseudorapidity[k]);
	  SumPfVsPtInitial[k] += pfVsPtInitial[ij];
	  SumPfVsPt[k] += pfVsPt[ij];
	  SumPfPt[k] += pfPt[ij];
	  if(printDebug2)std::cout<<SumPfVsPtInitial[k]<<endl;
	}// eta selection statement 

      }// eta bin loop

      hPfVsPtInitial_eta->Fill(pfVsPtInitial[ij],pfEta[ij]);
      hPfVsPt_eta->Fill(pfVsPt[ij],pfEta[ij]);
      hPfPt_eta->Fill(pfPt[ij],pfEta[ij]);

    }// pfcand loop   

    for(int k = 0;k<nedge_pseudorapidity-1;k++){
      
      hSumPfVsPtInitial[k]->Fill(SumPfVsPtInitial[k]);
      hSumPfVsPt[k]->Fill(SumPfVsPt[k]);
      hSumPfPt[k]->Fill(SumPfPt[k]);
      
    }// eta bin loop  
    
  }// entry loop
  
  //get the RMS value from the histograms at this level. 
  for(int k = 0;k<nedge_pseudorapidity-1;k++){
    hRMSSumPFVsPtIni_eta->Fill(edge_pseudorapidity[k],hSumPfVsPtInitial[k]->GetRMS(1));
    hRMSSumPFVsPt_eta->Fill(edge_pseudorapidity[k],hSumPfVsPt[k]->GetRMS(1));
    hRMSSumPFPt_eta->Fill(edge_pseudorapidity[k],hSumPfPt[k]->GetRMS(1));
    if(printDebug2)std::cout<<"RMS value of pfVsPtInitial bin"<<k<<" = "<<hSumPfVsPtInitial[k]->GetRMS(1)<<endl;
    if(printDebug2)std::cout<<"RMS value of pfVsPt bin"<<k<<" = "<<hSumPfVsPt[k]->GetRMS(1)<<endl;
    if(printDebug2)std::cout<<"RMS value of pfPt bin"<<k<<" = "<<hSumPfPt[k]->GetRMS(1)<<endl;
  }

  TCanvas *c1 = new TCanvas("c1","",600,400);
  makeMultiPanelCanvas(c1,3,1,0.0,0.0,0.2,0.15,0.07);
  c1->cd(2);
  hRMSSumPFVsPtIni_eta->SetTitle("RMS value of Sum$(PfVsPtInitial)");
  hRMSSumPFVsPtIni_eta->SetYTitle("GeV/c");
  hRMSSumPFVsPtIni_eta->SetXTitle("Detector #eta bins");
  hRMSSumPFVsPtIni_eta->Draw("col");
  drawText("After Sub, before Equ",0.1,0.2,15);
  c1->cd(3);
  hRMSSumPFVsPt_eta->SetTitle("RMS value of Sum$(PfVsPt)");
  hRMSSumPFVsPt_eta->SetXTitle("Detector #eta bins");
  hRMSSumPFVsPt_eta->Draw("col");
  drawText("Final; after Equ",0.2,0.2,15);
  c1->cd(1);
  hRMSSumPFPt_eta->SetTitle("RMS value of Sum$(PfPt)");
  hRMSSumPFPt_eta->SetXTitle("Detector #eta bins");
  hRMSSumPFPt_eta->Draw("col");
  drawText("Before Subtraction",0.25,0.2,15);
  c1->SaveAs("RMSvalue_SumPt_etaBins_bad_data.pdf","RECREATE");

  TFile f("bad_data_sumPFVsPT_histos.root","RECREATE");  
  f.cd();
  for(int k = 0;k<nedge_pseudorapidity-1;k++){
    hSumPfVsPtInitial[k]->Write();
    hSumPfVsPt[k]->Write();
    hSumPfPt[k]->Write();
  }
  hetabin_test->Write();
  hPfVsPtInitial_eta->Write();
  hPfVsPt_eta->Write();
  hPfPt_eta->Write();
  hRMSSumPFVsPtIni_eta->Write();
  hRMSSumPFVsPt_eta->Write();
  hRMSSumPFPt_eta->Write();
  f.Write();
  f.Close();


}
