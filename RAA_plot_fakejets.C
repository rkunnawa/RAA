// Raghav Kunnawalkam Elayavalli
// June 20 2014
// CERN

// 
// Macro making plots related to the fake jets issue coming from lower p_T triggers. 
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


TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname)
{
	TH1F *hF = (TH1F*)h->Clone(fHistname);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
		hF->SetBinContent(i,var);
		hF->SetBinError(i,0);
	}
	return hF;
}


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



void drawText(const char *text, float xp, float yp, int size){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(63);
	tex->SetTextSize(size);
	tex->SetTextColor(kBlack);
	tex->SetLineWidth(1);
	//tex->SetTextFont(42);
	tex->SetNDC();
	tex->Draw();
}


void putCMSPrel(double x, double y, double size){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Preliminary");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}


void putCMSSim(double x, double y, double size){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Simulation");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}


TLegend *myLegend(double x1,double y1,double x2, double y2)
{
	TLegend *leg = new TLegend(x1,y1,x2,y2);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	return leg; 

}

// Remove bins with error > central value
void cleanup(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val1 = h->GetBinContent(i);
		double valErr1 = h->GetBinError(i);
		if (valErr1>=val1) {
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}   

}

using namespace std;


void RAA_plot_fakejets(int radius = 3, char *algo = "Vs"){

  TH1::SetDefaultSumw2();
  
  TStopwatch timer;
  gStyle->SetOptStat(0);

  // get the file
  TFile *fin = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Output/PbPb_Jet55or65_ak%d_%s_fakejet_histos_combined_v5.root",radius,algo));

  TH1F *Jet55 = (TH1F*)fin->Get("hJet55");
  TH1F *Jet55_QA1 = (TH1F*)fin->Get("hJet55_QA1");
  TH1F *Jet55_QA2 = (TH1F*)fin->Get("hJet55_QA2");
  TH1F *Jet55_QA3 = (TH1F*)fin->Get("hJet55_QA3");
  TH1F *Jet55_QA1_2 = (TH1F*)fin->Get("hJet55_QA1_2");
  TH1F *Jet65 = (TH1F*)fin->Get("hJet65");
  TH1F *Jet65_QA1 = (TH1F*)fin->Get("hJet65_QA1");
  TH1F *Jet65_QA2 = (TH1F*)fin->Get("hJet65_QA2");
  TH1F *Jet65_QA3 = (TH1F*)fin->Get("hJet65_QA3");
  TH1F *Jet65_QA1_2 = (TH1F*)fin->Get("hJet65_QA1_2");
  TH1F *Jet55Fake = (TH1F*)fin->Get("hJet55Fake");
  TH1F *Jet55_only = (TH1F*)fin->Get("hJet55_only");
  TH1F *Jet65_only = (TH1F*)fin->Get("hJet65_only");
  
  TH1F *Jet55_trg = (TH1F*)fin->Get("hJet55_trg");
  TH1F *Jet55_trg_QA1 = (TH1F*)fin->Get("hJet55_trg_QA1");
  TH1F *Jet55_trg_QA2 = (TH1F*)fin->Get("hJet55_trg_QA1");
  TH1F *Jet55_trg_QA3 = (TH1F*)fin->Get("hJet55_trg_QA1");
  TH1F *Jet55_trg_QA1_2 = (TH1F*)fin->Get("hJet55_trg_QA1");
  TH1F *Jet65_trg = (TH1F*)fin->Get("hJet65_trg_QA1");
  TH1F *Jet65_trg_QA1 = (TH1F*)fin->Get("hJet65_trg_QA1");
  TH1F *Jet65_trg_QA2 = (TH1F*)fin->Get("hJet65_trg_QA1");
  TH1F *Jet65_trg_QA3 = (TH1F*)fin->Get("hJet65_trg_QA1");
  TH1F *Jet65_trg_QA1_2 = (TH1F*)fin->Get("hJet65_trg_QA1");

  TH1F *Jet55_trg_QA1_3 = (TH1F*)fin->Get("hJet55_trg_QA1_3"); 
  TH1F *Jet65_trg_QA1_3 = (TH1F*)fin->Get("hJet65_trg_QA1_3"); 
  TH1F *Jet55_QA1_3 = (TH1F*)fin->Get("hJet55_QA1_3"); 
  TH1F *Jet65_QA1_3 = (TH1F*)fin->Get("hJet65_QA1_3");

  

  // list of plots to make - 
  // 1) Jet55 all the QA cuts in one plot. with trigger and with object - (2) 
  // 2) Jet66 same as above.  (2)
  // 3) look at Jet80 plots.  (maybe 2)
  // 4) particular fake plot - would have to decide on how to classify them. 

  TDatime date;
  
  //plot 1
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogy();

  Jet55->SetMarkerColor(kBlue);
  Jet55->SetMarkerStyle(20);
  Jet55->SetXTitle("Jet p_{T} (GeV/c)");
  Jet55->SetYTitle("counts/(pt width)");
  
  Jet55->Draw();

  c1->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_HLT_HIJet65_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");

  //plot 2
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();

  Jet65->SetMarkerColor(kBlue);
  Jet65->SetMarkerStyle(20);
  Jet65->SetXTitle("Jet p_{T} (GeV/c)");
  Jet65->SetYTitle("counts/(pt width)");
  
  Jet65->Draw();

  c2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_HLT_HIJet65_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 3
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetLogy();
  
  Jet55->SetMarkerColor(kBlack);
  Jet55->SetMarkerStyle(20);
  Jet55->SetXTitle("Jet p_{T} (GeV/c)");
  Jet55->SetYTitle("counts/(pt width)");
  
  Jet55->Draw();

  Jet55_QA1->SetMarkerColor(2);
  Jet55_QA1->SetMarkerStyle(20);
  Jet55_QA1->Draw("same");

  Jet55_QA2->SetMarkerColor(3);
  Jet55_QA2->SetMarkerStyle(20);
  Jet55_QA2->Draw("same");

  Jet55_QA3->SetMarkerColor(4);
  Jet55_QA3->SetMarkerStyle(20);
  Jet55_QA3->Draw("same");

  Jet55_QA1_2->SetMarkerColor(5);
  Jet55_QA1_2->SetMarkerStyle(20);
  Jet55_QA1_2->Draw("same");

  Jet55_QA1_3->SetMarkerColor(6);
  Jet55_QA1_3->SetMarkerStyle(20);
  Jet55_QA1_3->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.1,0.8,20);
  drawText("HLT_HIJet55_v1",0.3,0.7,20);

  TLegend *title3 = myLegend(0.34,0.65,0.65,0.75);
  title3->AddEntry(Jet55,"NO QA","pl");
  title3->AddEntry(Jet55_QA1,"QA1","pl");
  title3->AddEntry(Jet55_QA2,"QA2","pl");
  title3->AddEntry(Jet55_QA3,"QA3","pl");
  title3->AddEntry(Jet55_QA1_2,"QA1 & QA2","pl");
  title3->AddEntry(Jet55_QA1_3,"QA1 & QA3","pl");
  title3->Draw();

  c3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2_src/Plots/PbPb_Jet55or65_ak%d_%s_Jet55_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  

  //plot 4
  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetLogy();
  
  Jet65->SetMarkerColor(kBlack);
  Jet65->SetMarkerStyle(20);
  Jet65->SetXTitle("Jet p_{T} (GeV/c)");
  Jet65->SetYTitle("counts/(pt width)");
  
  Jet65->Draw();

  Jet65_QA1->SetMarkerColor(2);
  Jet65_QA1->SetMarkerStyle(20);
  Jet65_QA1->Draw("same");

  Jet65_QA2->SetMarkerColor(3);
  Jet65_QA2->SetMarkerStyle(20);
  Jet65_QA2->Draw("same");

  Jet65_QA3->SetMarkerColor(4);
  Jet65_QA3->SetMarkerStyle(20);
  Jet65_QA3->Draw("same");

  Jet65_QA1_2->SetMarkerColor(5);
  Jet65_QA1_2->SetMarkerStyle(20);
  Jet65_QA1_2->Draw("same");

  Jet65_QA1_3->SetMarkerColor(6);
  Jet65_QA1_3->SetMarkerStyle(20);
  Jet65_QA1_3->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.1,0.8,20);
  drawText("HLT_HIJet65_v1",0.3,0.7,20);

  TLegend *title4 = myLegend(0.34,0.65,0.65,0.75);
  title4->AddEntry(Jet65,"NO QA","pl");
  title4->AddEntry(Jet65_QA1,"QA1","pl");
  title4->AddEntry(Jet65_QA2,"QA2","pl");
  title4->AddEntry(Jet65_QA3,"QA3","pl");
  title4->AddEntry(Jet65_QA1_2,"QA1 & QA2","pl");
  title4->AddEntry(Jet65_QA1_3,"QA1 & QA3","pl");
  title4->Draw();

  c4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2_src/Plots/PbPb_Jet55or65_ak%d_%s_Jet65_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 5
  TCanvas *c5 = new TCanvas("c5","",800,600);
  c5->SetLogy();
  
  Jet55_trg->SetMarkerColor(kBlack);
  Jet55_trg->SetMarkerStyle(20);
  Jet55_trg->SetXTitle("Jet p_{T} (GeV/c)");
  Jet55_trg->SetYTitle("counts/(pt width)");
  
  Jet55_trg->Draw();

  Jet55_trg_QA1->SetMarkerColor(2);
  Jet55_trg_QA1->SetMarkerStyle(20);
  Jet55_trg_QA1->Draw("same");

  Jet55_trg_QA2->SetMarkerColor(3);
  Jet55_trg_QA2->SetMarkerStyle(20);
  Jet55_trg_QA2->Draw("same");

  Jet55_trg_QA3->SetMarkerColor(4);
  Jet55_trg_QA3->SetMarkerStyle(20);
  Jet55_trg_QA3->Draw("same");

  Jet55_trg_QA1_2->SetMarkerColor(5);
  Jet55_trg_QA1_2->SetMarkerStyle(20);
  Jet55_trg_QA1_2->Draw("same");

  Jet55_trg_QA1_3->SetMarkerColor(6);
  Jet55_trg_QA1_3->SetMarkerStyle(20);
  Jet55_trg_QA1_3->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.3,0.8,20);
  drawText("HLT_HIJet55_v1, 55<= trigObject p_T <65",0.3,0.7,20);

  TLegend *title5 = myLegend(0.34,0.65,0.65,0.75);
  title5->AddEntry(Jet55_trg,"NO QA","pl");
  title5->AddEntry(Jet55_trg_QA1,"QA1","pl");
  title5->AddEntry(Jet55_trg_QA2,"QA2","pl");
  title5->AddEntry(Jet55_trg_QA3,"QA3","pl");
  title5->AddEntry(Jet55_trg_QA1_2,"QA1 & QA2","pl");
  title5->AddEntry(Jet55_trg_QA1_3,"QA1 & QA3","pl");
  title5->Draw();

  c5->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2_src/Plots/PbPb_Jet55or65_ak%d_%s_Jet55_trg_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");

  //plot 6
  TCanvas *c6 = new TCanvas("c6","",800,600);
  c6->SetLogy();
  
  Jet65_trg->SetMarkerColor(kBlack);
  Jet65_trg->SetMarkerStyle(20);
  Jet65_trg->SetXTitle("Jet p_{T} (GeV/c)");
  Jet65_trg->SetYTitle("counts/(pt width)");
  
  Jet65_trg->Draw();

  Jet65_trg_QA1->SetMarkerColor(2);
  Jet65_trg_QA1->SetMarkerStyle(20);
  Jet65_trg_QA1->Draw("same");

  Jet65_trg_QA2->SetMarkerColor(3);
  Jet65_trg_QA2->SetMarkerStyle(20);
  Jet65_trg_QA2->Draw("same");

  Jet65_trg_QA3->SetMarkerColor(4);
  Jet65_trg_QA3->SetMarkerStyle(20);
  Jet65_trg_QA3->Draw("same");

  Jet65_trg_QA1_2->SetMarkerColor(5);
  Jet65_trg_QA1_2->SetMarkerStyle(20);
  Jet65_trg_QA1_2->Draw("same");

  Jet65_trg_QA1_3->SetMarkerColor(6);
  Jet65_trg_QA1_3->SetMarkerStyle(20);
  Jet65_trg_QA1_3->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.3,0.8,20);
  drawText("HLT_HIJet65_v1, 65<= trigObj p_T <80",0.3,0.7,20);

  TLegend *title6 = myLegend(0.34,0.65,0.65,0.75);
  title6->AddEntry(Jet65_trg,"NO QA","pl");
  title6->AddEntry(Jet65_trg_QA1,"QA1","pl");
  title6->AddEntry(Jet65_trg_QA2,"QA2","pl");
  title6->AddEntry(Jet65_trg_QA3,"QA3","pl");
  title6->AddEntry(Jet65_trg_QA1_2,"QA1 & QA2","pl");
  title6->AddEntry(Jet65_trg_QA1_3,"QA1 & QA3","pl");
  title6->Draw();

  c6->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2_src/Plots/PbPb_Jet55or65_ak%d_%s_Jet65_trg_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");

}
