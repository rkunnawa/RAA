// Raghav Kunnawalkam Elayavalli
// Plotting the turn on curve from MinBias
//

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TDatime.h>
//#include "Headers/plot.h"

void CorrectBinWidth(TH1 *h1)
{
  for(int ix=1; ix<=h1->GetNbinsX(); ix++){
    float val    = h1->GetBinContent(ix);
    float valErr = h1->GetBinError(ix);
    float width  = h1->GetBinWidth(ix);
    if(val!=0){
      h1->SetBinContent(ix, val/width);
      h1->SetBinError(ix, valErr/width);
    }
  }
}
TLegend *getLegend(double x1, double y1, double x2, double y2)
{
  TLegend *leg = new TLegend(x1,y1,x2,y2,NULL,"BRNDC");
  leg->SetHeader("");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(1001);
  return leg;
}

void RAA_plot_triggerTurnon_v2() {

  TH1::SetDefaultSumw2();
  //gStyle->SetOptStat(0);

  TDatime date; 

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

  TFile * fin = TFile::Open("../../Output/PbPb_MinBiasUPC_supernova30_nojetid_prescl_akPuPF_20150420.root");
  
  TH1F * hDenominator  = (TH1F*)fin->Get("hDenominator_R3_cent0");
  TH1F * hNumerator_80  = (TH1F*)fin->Get("hNumerator_80_R3_cent0");
  TH1F * hNumerator_65  = (TH1F*)fin->Get("hNumerator_65_R3_cent0");
  TH1F * hNumerator_55  = (TH1F*)fin->Get("hNumerator_55_R3_cent0");

  CorrectBinWidth(hDenominator);
  CorrectBinWidth(hNumerator_55);
  CorrectBinWidth(hNumerator_65);
  CorrectBinWidth(hNumerator_80);

  TCanvas *c2 = new TCanvas("c2","Jet pT",808,635);
  c2->cd();
  gPad->SetLogy();
  hDenominator->Draw("p");
  hNumerator_55->Draw("psame");
  hNumerator_65->Draw("psame");
  hNumerator_80->Draw("psame");

  TGraphAsymmErrors *gr_turnon_jet55 = new TGraphAsymmErrors(hNumerator_55, hDenominator,"cl=0.683 b(1,1) mode");
  gr_turnon_jet55->SetTitle("");
  gr_turnon_jet55->GetYaxis()->SetTitle("Trigger turnon");
  gr_turnon_jet55->GetXaxis()->SetTitle("akPu3PF Reco Jet p_{T} (GeV/c)");
  gr_turnon_jet55->SetMarkerStyle(20);
  gr_turnon_jet55->SetMarkerColor(8);
  gr_turnon_jet55->SetLineColor(8);

  TGraphAsymmErrors *gr_turnon_jet65 = new TGraphAsymmErrors(hNumerator_65, hDenominator,"cl=0.683 b(1,1) mode");
  gr_turnon_jet65->SetMarkerStyle(20);
  gr_turnon_jet65->SetMarkerColor(4);
  gr_turnon_jet65->SetLineColor(4);

  TGraphAsymmErrors *gr_turnon_jet80 = new TGraphAsymmErrors(hNumerator_80, hDenominator,"cl=0.683 b(1,1) mode");
  gr_turnon_jet80->SetMarkerStyle(20);
  gr_turnon_jet80->SetMarkerColor(6);
  gr_turnon_jet80->SetLineColor(6);

  TLegend *leg = getLegend(0.60,0.15,0.85,0.35);
  leg->AddEntry(gr_turnon_jet55,"HLT_HIJet55","p");
  leg->AddEntry(gr_turnon_jet65,"HLT_HIJet65","p");
  leg->AddEntry(gr_turnon_jet80,"HLT_HIJet80","p");
  leg->SetTextSize(0.04);

  TCanvas *c3 = new TCanvas("c3","Turn on",808,635);
  c3->cd();
  gr_turnon_jet55->SetMaximum(1.05);
  gr_turnon_jet55->SetMinimum(0.0);
  gr_turnon_jet55->Draw("ap");
  gr_turnon_jet65->Draw("psame");
  gr_turnon_jet80->Draw("psame");
  leg->Draw();

  c3->SaveAs(Form("../../Plots/trigger_turnon_hlt_pbpb_supernova30_nojetid_prescl_akPu3PF_cent0_%d.pdf",date.GetDate()),"RECREATE");

  TH1F * hRatio_Jet80 = (TH1F*)hNumerator_80->Clone("hRatio_Jet80");
  hRatio_Jet80->Divide(hDenominator);
  TH1F * hRatio_Jet65 = (TH1F*)hNumerator_65->Clone("hRatio_Jet65");
  hRatio_Jet65->Divide(hDenominator);
  TH1F * hRatio_Jet55 = (TH1F*)hNumerator_55->Clone("hRatio_Jet55");
  hRatio_Jet55->Divide(hDenominator);

  TLine *line = new TLine(20,1,145,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas * c4 = new TCanvas("c4","",800,600);
  c4->cd();
  hRatio_Jet80->SetMarkerStyle(20);
  hRatio_Jet80->SetMarkerColor(6);
  hRatio_Jet80->SetYTitle("Trigger Contribution");
  hRatio_Jet80->SetXTitle("akPu3PF Reco Jet p_{T} (GeV/c)");
  hRatio_Jet80->SetAxisRange(0,1.2,"Y");
  hRatio_Jet80->SetAxisRange(20,140,"X");
  hRatio_Jet80->Draw();

  hRatio_Jet65->SetMarkerStyle(20);
  hRatio_Jet65->SetMarkerColor(4);
  hRatio_Jet65->SetAxisRange(20,100,"X");
  hRatio_Jet65->Draw("same");

  hRatio_Jet55->SetMarkerStyle(20);
  hRatio_Jet55->SetMarkerColor(8);
  hRatio_Jet55->SetAxisRange(20,90,"X");
  hRatio_Jet55->Draw("same");
  hRatio_Jet80->Draw("same");

  line->Draw();
  leg->Draw();
  
  c4->SaveAs(Form("../../Plots/trigger_spectra_ratio_supernov30_nojetid_prescl_akPu3PF_cent0_cent%d.pdf",date.GetDate()),"RECREATE");

}
