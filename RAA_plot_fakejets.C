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


void putCMSPrel(double x=0.1, double y=0.9, double size=0.04){
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
  TFile *fin = TFile::Open(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Output/PbPb_Jet55or65_ak%d_%s_fakejet_histos_combined_v8.root",radius,algo));

  TH1F *Jet55 = (TH1F*)fin->Get("hJet55");
  TH1F *Jet55_QA1 = (TH1F*)fin->Get("hJet55_QA1");
  TH1F *Jet55_QA2_a = (TH1F*)fin->Get("hJet55_QA2_a");
  TH1F *Jet55_QA2_b = (TH1F*)fin->Get("hJet55_QA2_b");
  TH1F *Jet55_QA3 = (TH1F*)fin->Get("hJet55_QA3");
  TH1F *Jet55_QA4 = (TH1F*)fin->Get("hJet55_QA4");
  TH1F *Jet55_QA1_2b = (TH1F*)fin->Get("hJet55_QA1_2b");
  TH1F *Jet65 = (TH1F*)fin->Get("hJet65");
  TH1F *Jet65_QA1 = (TH1F*)fin->Get("hJet65_QA1");
  TH1F *Jet65_QA2_a = (TH1F*)fin->Get("hJet65_QA2_a");
  TH1F *Jet65_QA2_b = (TH1F*)fin->Get("hJet65_QA2_b");
  TH1F *Jet65_QA3 = (TH1F*)fin->Get("hJet65_QA3");
  TH1F *Jet65_QA4 = (TH1F*)fin->Get("hJet65_QA4");
  TH1F *Jet65_QA1_2b = (TH1F*)fin->Get("hJet65_QA1_2b");
  TH1F *Jet55Fake = (TH1F*)fin->Get("hJet55Fake");
  TH1F *Jet55_only = (TH1F*)fin->Get("hJet55_only");
  TH1F *Jet65_only = (TH1F*)fin->Get("hJet65_only");
  
  TH1F *Jet55_trg = (TH1F*)fin->Get("hJet55_trg");
  TH1F *Jet55_trg_QA1 = (TH1F*)fin->Get("hJet55_trg_QA1");
  TH1F *Jet55_trg_QA2_a = (TH1F*)fin->Get("hJet55_trg_QA2_a");
  TH1F *Jet55_trg_QA2_b = (TH1F*)fin->Get("hJet55_trg_QA2_b");
  TH1F *Jet55_trg_QA3 = (TH1F*)fin->Get("hJet55_trg_QA3");
  TH1F *Jet55_trg_QA4 = (TH1F*)fin->Get("hJet55_trg_QA4");
  TH1F *Jet55_trg_QA1_2b = (TH1F*)fin->Get("hJet55_trg_QA1_2b");
  TH1F *Jet65_trg = (TH1F*)fin->Get("hJet65_trg_QA1");
  TH1F *Jet65_trg_QA1 = (TH1F*)fin->Get("hJet65_trg_QA1");
  TH1F *Jet65_trg_QA2_a = (TH1F*)fin->Get("hJet65_trg_QA2_a");
  TH1F *Jet65_trg_QA2_b = (TH1F*)fin->Get("hJet65_trg_QA2_b");
  TH1F *Jet65_trg_QA3 = (TH1F*)fin->Get("hJet65_trg_QA3");
  TH1F *Jet65_trg_QA4 = (TH1F*)fin->Get("hJet65_trg_QA4");
  TH1F *Jet65_trg_QA1_2b = (TH1F*)fin->Get("hJet65_trg_QA1_2b");

  TH1F *Jet55_trg_QA1_3 = (TH1F*)fin->Get("hJet55_trg_QA1_3"); 
  TH1F *Jet65_trg_QA1_3 = (TH1F*)fin->Get("hJet65_trg_QA1_3");

  //TH1F *Jet55_trg_QA2_3 = (TH1F*)fin
 
  TH1F *Jet55_QA1_3 = (TH1F*)fin->Get("hJet55_QA1_3"); 
  TH1F *Jet65_QA1_3 = (TH1F*)fin->Get("hJet65_QA1_3");

  

  // list of plots to make - 
  // 1) Jet55 all the QA cuts in one plot. with trigger and with object - (2) 
  // 2) Jet66 same as above.  (2)
  // 3) look at Jet80 plots.  (maybe 2)
  // 4) particular fake plot - would have to decide on how to classify them. 

  // make them all ratio plots which can show us the difference between cuts. 

  TDatime date;
  
  //plot 1
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogy();

  Jet55->SetMarkerColor(kBlue);
  Jet55->SetMarkerStyle(20);
  Jet55->SetTitle(" ");
  Jet55->SetXTitle("Jet p_{T} (GeV/c)");
  Jet55->SetYTitle("counts/(pt width)");
  //Jet55->SetAxisRange(10,630,"X");
  
  Jet55->Draw();

  Jet55_only->SetMarkerColor(kBlack);
  Jet55_only->SetMarkerStyle(25);
  Jet55_only->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoiseFilter, |vz|<15 && |#eta|<2",0.2,0.8,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);

  
  putCMSPrel();
  
  TLegend *title1 = myLegend(0.2,0.6,0.8,0.8);
  title1->AddEntry(Jet55,"HLT_HIJet55","pl");
  title1->AddEntry(Jet55_only,"HLT_HIJet55 && !HLT_HIJet65 && !HLT_HIJet80","pl");
  title1->SetTextSize(0.03);
  title1->Draw();

  c1->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_HLT_HIJet55_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 2
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();

  Jet65->SetMarkerColor(kBlue);
  Jet65->SetMarkerStyle(20);
  Jet65->SetTitle(" "); 
  Jet65->SetXTitle("Jet p_{T} (GeV/c)");
  Jet65->SetYTitle("counts/(pt width)");
  //Jet65->SetAxisRange(10,630,"X");
  
  Jet65->Draw();
  
  Jet65_only->SetMarkerColor(kBlack);
  Jet65_only->SetMarkerStyle(25);
  Jet65_only->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoiseFilter, |vz|<15 && |#eta|<2",0.2,0.8,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);

  putCMSPrel();

  TLegend *title2 = myLegend(0.3,0.6,0.85,0.8);
  title2->AddEntry(Jet65,"HLT_HIJet65","pl");
  title2->AddEntry(Jet65_only,"HLT_HIJet65 && !HLT_HIJet80","pl");
  title2->SetTextSize(0.03);
  title2->Draw();

  c2->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_HLT_HIJet65_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 3
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetLogy();
  
  Jet55->SetMarkerColor(kBlack);
  Jet55->SetMarkerStyle(31);
  Jet55->SetTitle(" ");
  Jet55->SetXTitle("Jet p_{T} (GeV/c)");
  Jet55->SetYTitle("counts/(pt width)");
  //Jet55->SetAxisRange(10,630,"X");
  Jet55->Draw();

  Jet55_QA1->SetMarkerColor(2);
  Jet55_QA1->SetMarkerStyle(24);
  Jet55_QA1->Draw("same");

  Jet55_QA2->SetMarkerColor(3);
  Jet55_QA2->SetMarkerStyle(25);
  Jet55_QA2->Draw("same");

  Jet55_QA3->SetMarkerColor(4);
  Jet55_QA3->SetMarkerStyle(27);
  Jet55_QA3->Draw("same");

  Jet55_QA1_2->SetMarkerColor(5);
  Jet55_QA1_2->SetMarkerStyle(28);
  Jet55_QA1_2->Draw("same");

  Jet55_QA1_3->SetMarkerColor(6);
  Jet55_QA1_3->SetMarkerStyle(30);
  Jet55_QA1_3->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.2,0.8,16);
  drawText("HLT_HIJet55_v1",0.3,0.85,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);

  putCMSPrel();

  TLegend *title3 = myLegend(0.45,0.3,0.85,0.8);
  title3->AddEntry(Jet55,"NO QA","pl");
  title3->AddEntry(Jet55_QA1,"QA1 = #frac{chMax}{jtpt}>0.01","pl");
  title3->AddEntry(Jet55_QA2,"QA2 = #frac{Max(chMax,neMax)}{Max(chSum,neSum)}<0.975","pl");
  title3->AddEntry(Jet55_QA3,"QA3 = #frac{jtpt}{trgObjpt}<3","pl");
  title3->AddEntry(Jet55_QA1_2,"QA1 & QA2","pl");
  title3->AddEntry(Jet55_QA1_3,"QA1 & QA3","pl");
  title3->SetTextSize(0.03);
  title3->Draw();

  c3->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_Jet55_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  

  //plot 4
  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetLogy();
  
  Jet65->SetMarkerColor(kBlack);
  Jet65->SetMarkerStyle(31);
  Jet65->SetTitle(" ");
  Jet65->SetXTitle("Jet p_{T} (GeV/c)");
  Jet65->SetYTitle("counts/(pt width)");
  //Jet65->SetAxisRange(10,630,"X");
  Jet65->Draw();

  Jet65_QA1->SetMarkerColor(2);
  Jet65_QA1->SetMarkerStyle(24);
  Jet65_QA1->Draw("same");

  Jet65_QA2->SetMarkerColor(3);
  Jet65_QA2->SetMarkerStyle(25);
  Jet65_QA2->Draw("same");

  Jet65_QA3->SetMarkerColor(4);
  Jet65_QA3->SetMarkerStyle(27);
  Jet65_QA3->Draw("same");

  Jet65_QA1_2->SetMarkerColor(5);
  Jet65_QA1_2->SetMarkerStyle(28);
  Jet65_QA1_2->Draw("same");

  Jet65_QA1_3->SetMarkerColor(6);
  Jet65_QA1_3->SetMarkerStyle(30);
  Jet65_QA1_3->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.2,0.8,16);
  drawText("HLT_HIJet65_v1",0.3,0.85,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);

  putCMSPrel();

  TLegend *title4 = myLegend(0.45,0.3,0.85,0.8);
  title4->AddEntry(Jet65,"NO QA","pl");
  title4->AddEntry(Jet65_QA1,"QA1 = #frac{chMax}{jtpt}>0.01","pl");
  title4->AddEntry(Jet65_QA2,"QA2 = #frac{Max(chMax,neMax)}{Max(chSum,neSum)}<0.975","pl");
  title4->AddEntry(Jet65_QA3,"QA3 = #frac{jtpt}{trgObjpt}<3","pl");
  title4->AddEntry(Jet65_QA1_2,"QA1 & QA2","pl");
  title4->AddEntry(Jet65_QA1_3,"QA1 & QA3","pl");
  title4->SetTextSize(0.03);
  title4->Draw();

  c4->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_Jet65_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 5
  TCanvas *c5 = new TCanvas("c5","",800,600);
  c5->SetLogy();
  
  Jet55_trg->SetMarkerColor(kBlack);
  Jet55_trg->SetMarkerStyle(31);
  Jet55_trg->SetTitle(" ");
  Jet55_trg->SetXTitle("Jet p_{T} (GeV/c)");
  Jet55_trg->SetYTitle("counts/(pt width)");
  //Jet55_trg->SetAxisRange(10,630,"X");
  Jet55_trg->Draw();

  Jet55_trg_QA1->SetMarkerColor(2);
  Jet55_trg_QA1->SetMarkerStyle(24);
  Jet55_trg_QA1->Draw("same");

  Jet55_trg_QA2->SetMarkerColor(3);
  Jet55_trg_QA2->SetMarkerStyle(25);
  Jet55_trg_QA2->Draw("same");

  Jet55_trg_QA3->SetMarkerColor(4);
  Jet55_trg_QA3->SetMarkerStyle(27);
  Jet55_trg_QA3->Draw("same");

  Jet55_trg_QA1_2->SetMarkerColor(5);
  Jet55_trg_QA1_2->SetMarkerStyle(28);
  Jet55_trg_QA1_2->Draw("same");

  Jet55_trg_QA1_3->SetMarkerColor(6);
  Jet55_trg_QA1_3->SetMarkerStyle(30);
  Jet55_trg_QA1_3->Draw("same");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.2,0.8,16);
  drawText("HLT_HIJet55_v1, 55<= trigObject p_{T} <65",0.3,0.85,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);

  putCMSPrel();

  TLegend *title5 = myLegend(0.45,0.3,0.85,0.8);
  title5->AddEntry(Jet55_trg,"NO QA","pl");
  title5->AddEntry(Jet55_trg_QA1,"QA1 = #frac{chMax}{jtpt}>0.01","pl");
  title5->AddEntry(Jet55_trg_QA2,"QA2 = #frac{Max(chMax,neMax)}{Max(chSum,neSum)}<0.975","pl");
  title5->AddEntry(Jet55_trg_QA3,"QA3 = #frac{jtpt}{trgObjpt}<3","pl");
  title5->AddEntry(Jet55_trg_QA1_2,"QA1 & QA2","pl");
  title5->AddEntry(Jet55_trg_QA1_3,"QA1 & QA3","pl");
  title5->SetTextSize(0.03);
  title5->Draw();

  c5->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_Jet55_trg_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 6
  TCanvas *c6 = new TCanvas("c6","",800,600);
  c6->SetLogy();
  
  Jet65_trg->SetMarkerColor(kBlack);
  Jet65_trg->SetMarkerStyle(31);
  Jet65_trg->SetTitle(" ");
  Jet65_trg->SetXTitle("Jet p_{T} (GeV/c)");
  Jet65_trg->SetYTitle("counts/(pt width)");
  //Jet65_trg->SetAxisRange(10,630,"X");
  Jet65_trg->Draw();
  //Jet65_trg->Print("base");

  Jet65_trg_QA1->SetMarkerColor(2);
  Jet65_trg_QA1->SetMarkerStyle(24);
  Jet65_trg_QA1->Draw("same");
  //Jet65_trg_QA1->Print("base");

  Jet65_trg_QA2->SetMarkerColor(3);
  Jet65_trg_QA2->SetMarkerStyle(25);
  Jet65_trg_QA2->Draw("same");
  //Jet65_trg_QA2->Print("base");

  Jet65_trg_QA3->SetMarkerColor(4);
  Jet65_trg_QA3->SetMarkerStyle(27);
  Jet65_trg_QA3->Draw("same");
  //Jet65_trg_QA3->Print("base");

  Jet65_trg_QA1_2->SetMarkerColor(5);
  Jet65_trg_QA1_2->SetMarkerStyle(28);
  Jet65_trg_QA1_2->Draw("same");
  //Jet65_trg_QA1_2->Print("base");

  Jet65_trg_QA1_3->SetMarkerColor(6);
  Jet65_trg_QA1_3->SetMarkerStyle(30);
  Jet65_trg_QA1_3->Draw("same");
  //Jet65_trg_QA1_3->Print("base");

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.2,0.8,16);
  drawText("HLT_HIJet65_v1, 65<= trigObj p_{T} <80",0.3,0.85,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);

  putCMSPrel();

  TLegend *title6 = myLegend(0.45,0.3,0.85,0.8);
  title6->AddEntry(Jet65_trg,"NO QA","pl");
  title6->AddEntry(Jet65_trg_QA1,"QA1 = #frac{chMax}{jtpt}>0.01","pl");
  title6->AddEntry(Jet65_trg_QA2,"QA2 = #frac{Max(chMax,neMax)}{Max(chSum,neSum)}<0.975","pl");
  title6->AddEntry(Jet65_trg_QA3,"QA3 = #frac{jtpt}{trgObjpt}<3","pl");
  title6->AddEntry(Jet65_trg_QA1_2,"QA1 & QA2","pl");
  title6->AddEntry(Jet65_trg_QA1_3,"QA1 & QA3","pl");
  title5->SetTextSize(0.03);
  title6->Draw();

  c6->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_Jet65_trg_diff_QAcuts_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 7
  TCanvas *c7 = new TCanvas("c7","",800,600);
  c7->SetLogy();

  Jet55Fake->SetMarkerStyle(20);
  Jet55Fake->SetMarkerColor(kBlack);
  Jet55Fake->SetXTitle("Jet p_{T} (GeV/c)");
  Jet55Fake->SetYTitle("counts/(pt width)");
  Jet55Fake->SetTitle(" ");

  Jet55Fake->Draw();

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.2,0.8,16);
  drawText("HLT_HIJet55_v1 && !HLT_HIJet65_v1 && !HLT_HIJet80, Jet p_{T} > 80",0.2,0.85,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);

  putCMSPrel();

  c7->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_Jet55_only_large_pt_eventSel_%d.pdf",radius,algo,date.GetDate()),"RECREATE");


  //plot 8
  TCanvas *c8 = new TCanvas("c8","",800,600);
  
  //make the ratio plots with QA#/noQA. 
  TH1F *Ratio_Jet55_trg_QA1 = (TH1F*)Jet55_trg_QA1->Clone("Ratio_Jet55_trg_QA1");
  Ratio_Jet55_trg_QA1->Divide(Jet55_trg);
  TH1F *Ratio_Jet55_trg_QA2 = (TH1F*)Jet55_trg_QA2->Clone("Ratio_Jet55_trg_QA2");
  Ratio_Jet55_trg_QA2->Divide(Jet55_trg);
  TH1F *Ratio_Jet55_trg_QA3 = (TH1F*)Jet55_trg_QA3->Clone("Ratio_Jet55_trg_QA3");
  Ratio_Jet55_trg_QA3->Divide(Jet55_trg);
  TH1F *Ratio_Jet55_trg_QA1_2 = (TH1F*)Jet55_trg_QA1_2->Clone("Ratio_Jet55_trg_QA1_2");
  Ratio_Jet55_trg_QA1_2->Divide(Jet55_trg);
  TH1F *Ratio_Jet55_trg_QA1_3 = (TH1F*)Jet55_trg_QA1_3->Clone("Ratio_Jet55_trg_QA1_3");
  Ratio_Jet55_trg_QA1_3->Divide(Jet55_trg);
  //TH1F *Ratio_Jet55_trg_QA3_2 = (TH1F*)Jet55_trg_QA1->Clone("Ratio_Jet55_trg_QA1");
  //Ratio_Jet55_trg_QA1->Divide(Jet55_trg);

  Ratio_Jet55_trg_QA1->SetMarkerColor(2);
  Ratio_Jet55_trg_QA1->SetMarkerStyle(24);//24,25,27,28,30
  Ratio_Jet55_trg_QA1->SetTitle(" ");
  Ratio_Jet55_trg_QA1->SetXTitle("Jet p_{T} (GeV/c)");
  Ratio_Jet55_trg_QA1->SetYTitle("#frac{}{}");
  Ratio_Jet55_trg_QA1->Draw();

  Ratio_Jet55_trg_QA2->SetMarkerColor(3);
  Ratio_Jet55_trg_QA2->SetMarkerStyle(25);
  Ratio_Jet55_trg_QA2->Draw("same");

  Ratio_Jet55_trg_QA3->SetMarkerColor(4);
  Ratio_Jet55_trg_QA3->SetMarkerStyle(27);
  Ratio_Jet55_trg_QA3->Draw("same");

  Ratio_Jet55_trg_QA1_2->SetMarkerColor(5);
  Ratio_Jet55_trg_QA1_2->SetMarkerStyle(28);
  Ratio_Jet55_trg_QA1_2->Draw("same");

  Ratio_Jet55_trg_QA1_3->SetMarkerColor(6);
  Ratio_Jet55_trg_QA1_3->SetMarkerStyle(30);
  Ratio_Jet55_trg_QA1_3->Draw("same");

  TLegend *title8 = myLegend(0.45,0.3,0.85,0.8);
  title8->AddEntry(Ratio_Jet55_trg_QA1,"QA1","pl");
  title8->AddEntry(Ratio_Jet55_trg_QA2,"QA2","pl");
  title8->AddEntry(Ratio_Jet55_trg_QA3,"QA2","pl");
  title8->AddEntry(Ratio_Jet55_trg_QA1_2,"QA2","pl");
  title8->AddEntry(Ratio_Jet55_trg_QA1_3,"QA2","pl");
  title8->SetTextSize(0.03);
  title8->Draw();

  drawText("pcollisionEventSelection, pHBHENoisefilter, |vz|<15 & |#eta|<2",0.2,0.8,16);
  putCMSPrel();
  drawText("HLT_HIJet55_V1, 55<=trgObjpt<65",0.2,0.85,16);
  drawText(Form("ak%s%dPF Jets",algo,radius),0.7,0.9,16);
  
  c8->SaveAs(Form("/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_8_HI_patch2/src/Plots/PbPb_Jet55or65_ak%d_%s_Jet55_trig_QA_ratio_%d",radius,algo,date.GetDate()));
  

}
