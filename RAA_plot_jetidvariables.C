{

  TH1::SetDefaultSumw2();
  gSystem->Load("Headers/plot.h");
  gStyle->SetOptStat(0);
  
  TFile * fData = TFile::Open("../../Output/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150320.root");
  TTree * Data_matched = (TTree*)fData->Get("matchedJets");
  TTree * Data_unmatched = (TTree*)fData->Get("unmatchedPFJets");

  TFile * fMC = TFile::Open("../../Output/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150326.root");
  TTree * MC_matched = (TTree*)fMC->Get("matchedJets");
  TTree * MC_unmatched = (TTree*)fMC->Get("unmatchedPFJets");

  const int nbins_cent = 6;
  const int boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  
  TH1F * hData_chMaxJtpt[nbins_cent], * hData_phMaxJtpt[nbins_cent], * hData_neMaxJtpt[nbins_cent], * hData_muMaxJtpt[nbins_cent], * hData_eMaxJtpt[nbins_cent];
  TH1F * hData_chSumJtpt[nbins_cent], * hData_phSumJtpt[nbins_cent], * hData_neSumJtpt[nbins_cent], * hData_muSumJtpt[nbins_cent], * hData_eSumJtpt[nbins_cent];  
  
  TH1F * hMC_chMaxJtpt[nbins_cent], * hMC_phMaxJtpt[nbins_cent], * hMC_neMaxJtpt[nbins_cent], * hMC_muMaxJtpt[nbins_cent], * hMC_eMaxJtpt[nbins_cent];
  TH1F * hMC_chSumJtpt[nbins_cent], * hMC_phSumJtpt[nbins_cent], * hMC_neSumJtpt[nbins_cent], * hMC_muSumJtpt[nbins_cent], * hMC_eSumJtpt[nbins_cent];  

  TCanvas * cchMaxJtpt[nbins_cent];
  
  for(int i = 0 ; i < nbins_cent ; ++i){

    hData_chMaxJtpt[i] = new TH1F(Form("hData_chMaxJtpt_cent%d",i),Form("Data chMax/jtpt in centrality bin %d",i),100,0,1);
    hMC_chMaxJtpt[i] = new TH1F(Form("hMC_chMaxJtpt_cent%d",i),Form("MC chMax/jtpt in centrality bin ",i),100,0,1);
    
    Data_matched->Draw(Form("chMax/pfpt>>hData_chMaxJtpt_cent%d",i),Form(""),"");
    MC_matched->Draw(Form("chMax/pfpt>>hMC_chMaxJtpt_cent%d",i),"weight*pfpt","");
    
    
    
    cchMaxJtpt[i] = new TCanvas(Form("cchMaxJtpt_cent%d",i),"",800,600);
    cchMaxJtpt[i]->SetLogy();
    hMC_chMaxJtpt[i]->SetMarkerColor(kRed);
    hMC_chMaxJtpt[i]->SetMarkerStyle(25);
    hMC_chMaxJtpt[i]->Print("base");
    hMC_chMaxJtpt[i]->SetTitle(" ");
    hMC_chMaxJtpt[i]->Draw();

    hData_chMaxJtpt[i]->SetMarkerStyle(24);
    hData_chMaxJtpt[i]->SetMarkerColor(kBlack);
    hData_chMaxJtpt[i]->DrawNormalized("same");

    TLegend * leg = myLegend(0.15,0.15,0.3,0.3);
    leg->AddEntry(hMC_chMaxJtpt[1],"MC","pl");
    leg->AddEntry(hData_chMaxJtpt[1],"Data","pl");
    leg->SetTextSize(0.04);
    leg->Draw();
    drawText(Form("%2.0f - %2.0f % ", 2.5 * boundaries_cent[i], 2.5 * boundaries_cent[i+1]), 0.2,0.4,16);
    
    cchMaxJtpt[i]->SaveAs(Form("chMaxJtpt_cent%d.pdf",i),"RECREATE");

  }

}
