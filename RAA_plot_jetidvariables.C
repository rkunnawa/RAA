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


  for(int i = 0 ; i < nbins_cent ; ++i){

    hData_chMaxJtpt[i] = new TH1F(Form("hData_chMaxJtpt_cent%d",i),Form("Data chMax/jtpt in centrality bin %d",i),100,0,1);
    hMC_chMaxJtpt[i] = new TH1F(Form("hMC_chMaxJtpt_cent%d",i),Form("MC chMax/jtpt in centrality bin ",i),100,0,1);

  }
  
  // 1 - Data, 2 - MC
  Float_t pfpt_1, pfpt_2;
  Float_t pfrefpt_2;
  Float_t calopt_1, calopt_2;
  Float_t chMax_1, chMax_2;
  Float_t phMax_1, phMax_2;
  Float_t neMax_1, neMax_2;
  Float_t muMax_1, muMax_2;
  Float_t eMax_1, eMax_2;
  Float_t chSum_1, chSum_2;
  Float_t phSum_1, phSum_2;
  Float_t neSum_1, neSum_2;
  Float_t muSum_1, muSum_2;
  Float_t eSum_1, eSum_2; 
  Int_t jet55_1, jet65_1, jet80_1;
  Int_t jet55_p_1, jet65_p_1, jet80_p_1;
  Int_t jet55_2, jet65_2, jet80_2;
  Int_t jet55_p_2, jet65_p_2, jet80_p_2;
  Int_t hiBin_1, hiBin_2;
  Float_t weight_2;

  Data_matched->SetBranchAddress("calopt",&calopt_1);
  Data_matched->SetBranchAddress("pfpt",&pfpt_1);
  Data_matched->SetBranchAddress("eMax",&eMax_1);
  Data_matched->SetBranchAddress("phMax",&phMax_1);
  Data_matched->SetBranchAddress("muMax",&muMax_1);
  Data_matched->SetBranchAddress("neMax",&neMax_1);
  Data_matched->SetBranchAddress("chMax",&chMax_1);
  Data_matched->SetBranchAddress("chSum",&chSum_1);
  Data_matched->SetBranchAddress("phSum",&phSum_1);
  Data_matched->SetBranchAddress("neSum",&neSum_1);
  Data_matched->SetBranchAddress("muSum",&muSum_1);
  Data_matched->SetBranchAddress("eSum",&eSum_1);
  Data_matched->SetBranchAddress("hiBin",&hiBin_1);
  Data_matched->SetBranchAddress("jet55",&jet55_1);
  Data_matched->SetBranchAddress("jet65",&jet65_1);
  Data_matched->SetBranchAddress("jet80",&jet80_1);
  Data_matched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_matched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_matched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  
  Data_unmatched->SetBranchAddress("pfpt",&pfpt_1);
  Data_unmatched->SetBranchAddress("eMax",&eMax_1);
  Data_unmatched->SetBranchAddress("phMax",&phMax_1);
  Data_unmatched->SetBranchAddress("muMax",&muMax_1);
  Data_unmatched->SetBranchAddress("neMax",&neMax_1);
  Data_unmatched->SetBranchAddress("chMax",&chMax_1);
  Data_unmatched->SetBranchAddress("chSum",&chSum_1);
  Data_unmatched->SetBranchAddress("phSum",&phSum_1);
  Data_unmatched->SetBranchAddress("neSum",&neSum_1);
  Data_unmatched->SetBranchAddress("muSum",&muSum_1);
  Data_unmatched->SetBranchAddress("eSum",&eSum_1);
  Data_unmatched->SetBranchAddress("hiBin",&hiBin_1);
  Data_unmatched->SetBranchAddress("jet55",&jet55_1);
  Data_unmatched->SetBranchAddress("jet65",&jet65_1);
  Data_unmatched->SetBranchAddress("jet80",&jet80_1);
  Data_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_unmatched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_unmatched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  
  MC_matched->SetBranchAddress("calopt",&calopt_2);
  MC_matched->SetBranchAddress("pfpt",&pfpt_2);
  MC_matched->SetBranchAddress("eMax",&eMax_2);
  MC_matched->SetBranchAddress("phMax",&phMax_2);
  MC_matched->SetBranchAddress("muMax",&muMax_2);
  MC_matched->SetBranchAddress("neMax",&neMax_2);
  MC_matched->SetBranchAddress("chMax",&chMax_2);
  MC_matched->SetBranchAddress("chSum",&chSum_2);
  MC_matched->SetBranchAddress("phSum",&phSum_2);
  MC_matched->SetBranchAddress("neSum",&neSum_2);
  MC_matched->SetBranchAddress("muSum",&muSum_2);
  MC_matched->SetBranchAddress("eSum",&eSum_2);
  MC_matched->SetBranchAddress("hiBin",&hiBin_2);
  MC_matched->SetBranchAddress("weight",&weight_2);
  MC_matched->SetBranchAddress("pfrefpt",&pfrefpt_2);
  MC_matched->SetBranchAddress("jet55",&jet55_2);
  MC_matched->SetBranchAddress("jet65",&jet65_2);
  MC_matched->SetBranchAddress("jet80",&jet80_2);
  
  MC_unmatched->SetBranchAddress("pfpt",&pfpt_2);
  MC_unmatched->SetBranchAddress("eMax",&eMax_2);
  MC_unmatched->SetBranchAddress("phMax",&phMax_2);
  MC_unmatched->SetBranchAddress("muMax",&muMax_2);
  MC_unmatched->SetBranchAddress("neMax",&neMax_2);
  MC_unmatched->SetBranchAddress("chMax",&chMax_2);
  MC_unmatched->SetBranchAddress("chSum",&chSum_2);
  MC_unmatched->SetBranchAddress("phSum",&phSum_2);
  MC_unmatched->SetBranchAddress("neSum",&neSum_2);
  MC_unmatched->SetBranchAddress("muSum",&muSum_2);
  MC_unmatched->SetBranchAddress("eSum",&eSum_2);
  MC_unmatched->SetBranchAddress("hiBin",&hiBin_2);
  MC_unmatched->SetBranchAddress("weight",&weight_2);
  MC_unmatched->SetBranchAddress("pfrefpt",&pfrefpt_2);
  MC_unmatched->SetBranchAddress("jet55",&jet55_2);
  MC_unmatched->SetBranchAddress("jet65",&jet65_2);
  MC_unmatched->SetBranchAddress("jet80",&jet80_2);
  
  
  // data loop
  long entries = Data_matched->GetEntries();
  //entries = 1;
  cout<<"matched Data ntuple "<<endl;
  for(int nentry = 0; nentry < entries; ++nentry){
    
    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    Data_matched->GetEntry(nentry);



  }

  
  // data unmatched loop:
  entries = Data_unmatched->GetEntries();
  //entries = 1;
  cout<<"Unmatched Data ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    Data_unmatched->GetEntry(nentry);


  }


  
  entries = MC_matched->GetEntries();
  //entries = 1;
  // MC loop
  cout<<" looking at matched MC ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    MC_matched->GetEntry(nentry);


  }


  
  entries = MC_unmatched->GetEntries();
  //entries = 1;
  // MC loop
  cout<<" looking at unmatched MC ntuple"<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    MC_unmatched->GetEntry(nentry);


  }

  
  
  
  TCanvas * cchMaxJtpt[nbins_cent];
  
  for(int i = 0 ; i < nbins_cent ; ++i){

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
