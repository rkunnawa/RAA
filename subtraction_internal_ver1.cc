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
#include <TH2D.h>
#include <TMath.h>

#define NOBJECT_MAX 16384

double ue_predictor_pf[3][15][5][2][82];
double ue_interpolation_pf0[15][344];
double ue_interpolation_pf1[15][344];
double ue_interpolation_pf2[15][82];

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

void subtraction_internal_ver1(const char *filename = "/Users/keraghav/WORK/RAA/Output/PbPb_data_bad_events.root", const int data = 1, const int calorimetric = 0)
{
	gStyle->SetPalette(55);

	static const size_t nfourier = 5;

	const char *root_tree_name = calorimetric ?
		"rechitanalyzer/tower" : "pfcandAnalyzer/pfTree";

	TFile *root_file = TFile::Open(filename);
	TTree *root_tree = dynamic_cast<TTree *>(gDirectory->Get(root_tree_name));

	Int_t nPFpart;
	Int_t pfId[NOBJECT_MAX];
	Float_t pfPt[NOBJECT_MAX];
	Float_t pfEta[NOBJECT_MAX];
	Float_t pfPhi[NOBJECT_MAX];
	Float_t pfArea[NOBJECT_MAX];

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
		root_tree->SetBranchAddress("pfEta", pfEta);
		root_tree->SetBranchAddress("pfPhi", pfPhi);
		root_tree->SetBranchAddress("pfArea", pfArea);
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
	static const double edge_pseudorapidity[nedge_pseudorapidity] = {
		-5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522,
		0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191
	};

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

	size_t nentries = root_tree->GetEntries();
	nentries = 20;

	TCanvas canvas0("canvas0", "", 1440, 480);

	const double left_margin = 0.0657133;
	const double right_margin = 0.0023535;
	const double top_margin = 0.05 / 1440 * 480;
	const double bottom_margin = 0.109076 * 1.17;
	const size_t npanel_x = 3;
	const size_t npanel_y = 1;

	std::vector<TPad *> pad;

	for (size_t i = 0; i < 3; i++) {
		canvas0.cd();

		char buf[4096];

		snprintf(buf, 4096, "pad%lu", i);

		pad.push_back(
			new TPad(buf, "",
					 (1 - left_margin - right_margin) / npanel_x * i,
					 (1 - top_margin - bottom_margin) / npanel_y * 0,
					 left_margin + right_margin +
					 (1 - left_margin - right_margin) / npanel_x * (i + 1),
					 top_margin + bottom_margin +
					 (1 - top_margin - bottom_margin) / npanel_y * 1));
		pad.back()->SetLeftMargin(
			left_margin / (left_margin + right_margin +
						   (1 - left_margin - right_margin) / npanel_x));
		pad.back()->SetRightMargin(
			right_margin / (left_margin + right_margin +
							(1 - left_margin - right_margin) / npanel_x));
		pad.back()->SetTopMargin(
			top_margin / (top_margin + bottom_margin +
						  (1 - top_margin - bottom_margin) / npanel_y));
		pad.back()->SetBottomMargin(
			bottom_margin / (top_margin + bottom_margin +
							 (1 - top_margin - bottom_margin) / npanel_y));
		pad.back()->SetFillStyle(0);
		pad.back()->SetFillColor(0);
		pad.back()->Draw();
		pad.back()->Modified();
	}

	TH2D root_histogram0("root_histogram0", "", ncms_hcal_edge_pseudorapidity - 1, cms_hcal_edge_pseudorapidity, 36, -TMath::Pi(), TMath::Pi());
	TH2D root_histogram1("root_histogram1", "", ncms_hcal_edge_pseudorapidity - 1, cms_hcal_edge_pseudorapidity, 36, -TMath::Pi(), TMath::Pi());
	TH2D root_histogram2("root_histogram1", "", ncms_hcal_edge_pseudorapidity - 1, cms_hcal_edge_pseudorapidity, 36, -TMath::Pi(), TMath::Pi());

	for (size_t i = 0; i < nentries; i++) {
		root_tree->GetEntry(i);
		hiTree->GetEntry(i);
		t->GetEntry(i);

		root_histogram0.Reset();
		root_histogram1.Reset();
		root_histogram2.Reset();

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
		fprintf(stderr, "%s:%d: %f %f\n", __FILE__, __LINE__, feature[0], feature[0] / scale[0]);
		for (size_t k = 1; k < nfourier; k++) {
			feature[2 * k - 1] = 0;
			for (size_t l = 0; l < nreduced_id; l++) {
			feature[2 * k - 1] += scale[k] *
				(perp_fourier[0                       ][l][k][0] +
				 perp_fourier[nedge_pseudorapidity - 2][l][k][0]);
			}
			fprintf(stderr, "%s:%d: %f %f %f\n", __FILE__, __LINE__, perp_fourier[0                       ][2][k][0], perp_fourier[nedge_pseudorapidity - 2][2][k][0], feature[2 * k - 1]);
			feature[2 * k] = 0;
			for (size_t l = 0; l < nreduced_id; l++) {
			feature[2 * k] += scale[k] *
				(perp_fourier[0                       ][l][k][1] +
				 perp_fourier[nedge_pseudorapidity - 2][l][k][1]);
			}
			fprintf(stderr, "%s:%d: %f %f %f\n", __FILE__, __LINE__, perp_fourier[0                       ][2][k][1], perp_fourier[nedge_pseudorapidity - 2][2][k][1], feature[2 * k]);
		}

#if 0
		const double event_plane = atan2(feature[4], feature[3]);
		const double v2 =
			sqrt(feature[3] * feature[3] +
				 feature[4] * feature[4]) / feature[0];
#endif

		fprintf(stderr, "%s:%d: %f %f\n", __FILE__, __LINE__, hiBin * 0.5, sqrt(feature[3] * feature[3] +
				 feature[4] * feature[4]));

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

#if 0
							// This looks at a specific flow component and see how the polynomial is evaluated
							if (j == 0 && predictor_index == 3 && l == 0 && m == 0) {
								//fprintf(stderr, "%s:%d: %f %f\n", __FILE__, __LINE__, perp_fourier[0][2][2][0], perp_fourier[nedge_pseudorapidity - 2][2][2][1]);
								fprintf(stderr, "%s:%d: << %f %f %f %f %f %f %f\n", __FILE__, __LINE__, feature[0], feature[1], feature[2], feature[3], feature[4], u, perp_fourier[predictor_index][j][l][m]);
							}
#endif

							pred += u * (l == 0 ? 1.0 : 2.0) *
								(m == 0 ? cos(l * azimuth) :
								 sin(l * azimuth));
						}
					}

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
			//fprintf(stderr, "%s:%d: %.8e %.8e %.8e %.8e %.8e\n", __FILE__, __LINE__, hiBin * 0.5, pfEta[k], pfPhi[k], pfPt[k], density * pfArea[k]);
			root_histogram0.Fill(pfEta[k], pfPhi[k], pfPt[k] - density * pfArea[k]);
			root_histogram1.Fill(pfEta[k], pfPhi[k], pfPt[k]);
		}


		for (int kx = 1; kx < root_histogram2.GetXaxis()->GetNbins() + 1; kx++) {
			for (int ky = 1; ky < root_histogram2.GetYaxis()->GetNbins() + 1; ky++) {
				double pfEta_k = root_histogram2.GetXaxis()->GetBinCenter(kx);
				double pfPhi_k = root_histogram2.GetYaxis()->GetBinCenter(ky);

			int predictor_index = -1;
			int interpolation_index = -1;
			double density = 0;

			if (pfEta_k >= edge_pseudorapidity[0] &&
				pfEta_k < edge_pseudorapidity[nedge_pseudorapidity - 1]) {
				std::vector<double>::const_iterator p = std::lower_bound(edge_pseudorapidity_v.begin(), edge_pseudorapidity_v.end(), pfEta_k);

				predictor_index = (p - edge_pseudorapidity_v.begin()) - 1;
			}

			for (size_t j = 0; j < nreduced_id; j++) {
				if (j == 2) {
					// HCAL
					if (pfEta_k >=
						cms_hcal_edge_pseudorapidity[0] &&
						pfEta_k <
						cms_hcal_edge_pseudorapidity[ncms_hcal_edge_pseudorapidity - 1]) {
						std::vector<double>::const_iterator p = std::lower_bound(cms_hcal_edge_pseudorapidity_v.begin(), cms_hcal_edge_pseudorapidity_v.end(), pfEta_k);

						interpolation_index = (p - cms_hcal_edge_pseudorapidity_v.begin()) - 1;
					}
				}
				else {
					// Tracks or ECAL clusters
					if (pfEta_k >=
						cms_ecal_edge_pseudorapidity[0] &&
						pfEta_k <
						cms_ecal_edge_pseudorapidity[ncms_ecal_edge_pseudorapidity - 1]) {
						std::vector<double>::const_iterator p = std::lower_bound(cms_ecal_edge_pseudorapidity_v.begin(), cms_ecal_edge_pseudorapidity_v.end(), pfEta_k);

						interpolation_index = (p - cms_ecal_edge_pseudorapidity_v.begin()) - 1;
					}
				}

				if (predictor_index >= 0 && interpolation_index >= 0) {
					// Calculate the aggregated prediction and
					// interpolation for the pseudorapidity segment

					const double azimuth = pfPhi_k;
					const double (*p)[2][82] =
						ue_predictor_pf[j][predictor_index];
					double pred = 0;

					for (size_t l = 0; l < nfourier; l++) {
						for (size_t m = 0; m < 2; m++) {
							float u = p[l][m][0];

							for (size_t n = 0; n < 2 * nfourier - 1; n++) {
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

#if 0
							// This looks at a specific flow component and see how the polynomial is evaluated
							if (j == 0 && predictor_index == 3 && l == 0 && m == 0) {
								//fprintf(stderr, "%s:%d: %f %f\n", __FILE__, __LINE__, perp_fourier[0][2][2][0], perp_fourier[nedge_pseudorapidity - 2][2][2][1]);
								fprintf(stderr, "%s:%d: << %f %f %f %f %f %f %f\n", __FILE__, __LINE__, feature[0], feature[1], feature[2], feature[3], feature[4], u, perp_fourier[predictor_index][j][l][m]);
							}
#endif

							pred += u * (l == 0 ? 1.0 : 2.0) *
								(m == 0 ? cos(l * azimuth) :
								 sin(l * azimuth));
						}
					}

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
			root_histogram2.SetBinContent(root_histogram2.FindBin(pfEta_k, pfPhi_k), density * root_histogram2.GetXaxis()->GetBinWidth(root_histogram2.GetXaxis()->FindFixBin(pfEta_k)) * root_histogram2.GetYaxis()->GetBinWidth(root_histogram2.GetYaxis()->FindFixBin(pfPhi_k)));
			}
		}

		pad[0]->cd();

		pad[0]->SetRightMargin(0.125);
		pad[0]->SetLeftMargin(0.0625);
		pad[0]->SetBottomMargin(0.0625);

		root_histogram1.SetMinimum(-100);
		root_histogram1.SetMaximum(100);
		root_histogram1.SetContour(252);
		root_histogram1.Draw("colz");

		pad[0]->Modified();

		pad[1]->cd();

		pad[1]->SetRightMargin(0.125);
		pad[1]->SetLeftMargin(0.0625);
		pad[1]->SetBottomMargin(0.0625);

		root_histogram2.SetMinimum(-100);
		root_histogram2.SetMaximum(100);
		root_histogram2.SetContour(252);
		root_histogram2.Draw("colz");

		pad[2]->cd();

		pad[2]->SetRightMargin(0.125);
		pad[2]->SetLeftMargin(0.0625);
		pad[2]->SetBottomMargin(0.0625);

		root_histogram0.SetMinimum(-100);
		root_histogram0.SetMaximum(100);
		root_histogram0.SetContour(252);
		root_histogram0.Draw("colz");

		pad[2]->Modified();

		char buf[4096];

		snprintf(buf, 4096, "subtraction_data_ver1_%lu.png", i);

		canvas0.SaveAs(buf);
	}

	root_file->Close();

	gSystem->Exit(0);
}
