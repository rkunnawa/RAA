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
#include <TMath.h>
#include <TH2D.h>
#include <TF1.h>

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

void subtraction_internal_ver4(const int data = 1, const int calorimetric = 0, const size_t predictor_index = 8)
{
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

	TCanvas canvas0("canvas0", "", 960, 960);

	const double left_margin = 0.0657133;
	const double right_margin = 0.0023535;
	const double top_margin = 0.05 / 960 * 960;
	const double bottom_margin = 0.109076 * 1.17;
	const size_t npanel_x = 5;
	const size_t npanel_y = 5;

	std::vector<TPad *> pad;

	for (size_t i = 0; i < 5; i++) {
		for (size_t j = 0; j < 5; j++) {
			canvas0.cd();

			char buf[4096];

			snprintf(buf, 4096, "pad%lu", i * 5 + j);

			pad.push_back(new TPad(buf, "",
					(1 - left_margin - right_margin) / npanel_x * i,
					(1 - top_margin - bottom_margin) / npanel_y * (4 - j),
					left_margin + right_margin +
					(1 - left_margin - right_margin) / npanel_x * (i + 1),
					top_margin + bottom_margin +
					(1 - top_margin - bottom_margin) / npanel_y * (5 - j)));
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
	}

	canvas0.SetRightMargin(0.125);
	canvas0.SetLeftMargin(0.0625);
	canvas0.SetBottomMargin(0.0625);

	std::vector<TF1 *> function;
	std::vector<TH1D *> root_histogram;

	const size_t nfourier = 5;
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

	for (size_t i = 0; i < 5; i++) {
		for (size_t j = 0; j < 5; j++) {
			pad[i * 5 + j]->cd();

			const size_t l = j;
			const size_t function_size_begin = function.size();

			for (size_t k = 0; k < 3; k++) {
				const double (*p)[2][82] =
					ue_predictor_pf[k][predictor_index];

				for (size_t n = (i == 0 ? 0 : 2 * i - 1);
					 n < (i == 0 ? 1 : 2 * i + 1); n++) {
					for (size_t m = 0; m < 2; m++) {
						char buf[4096];

						snprintf(buf, 4096, "function%lu", function.size());
						function.push_back(new TF1(buf, "[0]+[1]*(x*[10])+[2]*(x*[10])^2+[3]*(x*[10])^3+[4]*(x*[10])^4+[5]*(x*[10])^5+[6]*(x*[10])^6+[7]*(x*[10])^7+[8]*(x*[10])^8+[9]*(x*[10])^9", 0, 2 / scale[i]));
						function.back()->SetParameter(0, p[l][m][0]);
						for (size_t o = 1; o <= 9; o++) {
							function.back()->SetParameter(o, p[l][m][9 * n + o]);
						}
						function.back()->SetParameter(10, scale[i]);
						if (k == 0) {
							function.back()->SetLineColor(kRed);
						}
						else if (k == 1) {
							function.back()->SetLineColor(kBlue);
						}
						else if (k == 2) {
							function.back()->SetLineColor(kGreen);
						}
						if (m != 0) {
							function.back()->SetLineStyle(2);
						}
						function.back()->SetLineWidth(2);
					}
				}
			}

			char buf[4096];

			snprintf(buf, 4096, "root_histogram%lu", root_histogram.size());

			root_histogram.push_back(new TH1D(buf, "", 1, 0, 2 / scale[i]));

			if (j == 0) {
				root_histogram.back()->GetYaxis()->SetRangeUser(-50, 250);
			}
			else {
				root_histogram.back()->GetYaxis()->SetRangeUser(-50, 50);
			}
			root_histogram.back()->Draw("");

			for (std::vector<TF1 *>::iterator iterator = function.begin() + function_size_begin; iterator != function.end(); iterator++) {
				(*iterator)->Draw("same");
			}

			pad[i * 5 + j]->Modified();
		}
	}

	char buf[4096];

	snprintf(buf, 4096, "polynom_%d.png", 0);
	canvas0.SaveAs(buf);

	gSystem->Exit(0);
}
