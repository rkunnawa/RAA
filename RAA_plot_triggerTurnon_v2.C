// Raghav Kunnawalkam Elayavalli
// Plotting the turn on curve from MinBias
//

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
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

  TFile * fin = TFile::Open("../../Output/PbPb_MinBiasUPC_akPuCalo_20150318.root");
  
  TH1F * hDenominator  = (TH1F*)fin->Get("hDenominator_R2");
  TH1F * hNumerator_80  = (TH1F*)fin->Get("hNumerator_80_R2");
  TH1F * hNumerator_65  = (TH1F*)fin->Get("hNumerator_65_R2");
  TH1F * hNumerator_55  = (TH1F*)fin->Get("hNumerator_55_R2");

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
  gr_turnon_jet55->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
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

  TLegend *leg = getLegend(0.40,0.30,0.75,0.55);
  leg->AddEntry(gr_turnon_jet55,"HIHLT_Jet55","p");
  leg->AddEntry(gr_turnon_jet65,"HIHLT_Jet65","p");
  leg->AddEntry(gr_turnon_jet80,"HIHLT_Jet80","p");

  TCanvas *c3 = new TCanvas("c3","Turn on",808,635);
  c3->cd();
  gr_turnon_jet55->SetMaximum(1.05);
  gr_turnon_jet55->SetMinimum(0.0);
  gr_turnon_jet55->Draw("ap");
  gr_turnon_jet65->Draw("psame");
  gr_turnon_jet80->Draw("psame");
  leg->Draw();

  c3->SaveAs("../../Plots/trigger_turnon_hlt_pbpb_akPu3Calo.pdf","RECREATE");







}
