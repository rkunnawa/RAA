{

  TH1::SetDefaultSumw2();
  
  gSystem->Load("Headers/plot.h");

  TFile * fData = TFile::Open("../../Output/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150320.root");
  TTree * Data_matched = (TTree*)fData->Get("matchedJets");
  TTree * Data_unmatched = (TTree*)fData->Get("unmatchedPFJets");

  TFile * fMC = TFile::Open("../../Output/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150319_9.root");
  TTree * MC_matched = (TTree*)fMC->Get("matchedJets");
  TTree * MC_unmatched = (TTree*)fMC->Get("unmatchedPFJets");

  TH1F * hMC_Jet55_noCut = new TH1F("hMC_Jet55_noCut","data from matched jets without any jet ID cut",1000,0,1000);
  TH1F * hMC_Jet55_CutA = new TH1F("hMC_Jet55_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",1000,0,1000);
  TH1F * hMC_Jet55_CutB = new TH1F("hMC_Jet55_CutB","data from matched jets with Jet ID cut: slant line calopt/pfpt = 0.5, eMax/pfpt = 0 to calopt/pfpt = 1.25, eMax/pfpt = 1.1",1000,0,1000);

  TH1F * hMC_Jet65_noCut = new TH1F("hMC_Jet65_noCut","data from matched jets without any jet ID cut",1000,0,1000);
  TH1F * hMC_Jet65_CutA = new TH1F("hMC_Jet65_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",1000,0,1000);
  TH1F * hMC_Jet65_CutB = new TH1F("hMC_Jet65_CutB","data from matched jets with Jet ID cut: slant line calopt/pfpt = 0.5, eMax/pfpt = 0 to calopt/pfpt = 1.25, eMax/pfpt = 1.1",1000,0,1000);

  TH1F * hMC_Jet80_noCut = new TH1F("hMC_Jet80_noCut","data from matched jets without any jet ID cut",1000,0,1000);
  TH1F * hMC_Jet80_CutA = new TH1F("hMC_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.4 calopt/pfpt from eMax/Sumcand 0 till 0.9  and then calopt/pfpt > 0.85",1000,0,1000);
  TH1F * hMC_Jet80_CutB = new TH1F("hMC_Jet80_CutB","data from matched jets with Jet ID cut: slant line calopt/pfpt = 0.5, eMax/pfpt = 0 to calopt/pfpt = 1.25, eMax/pfpt = 1.1",1000,0,1000);

  TH1F * hData_Jet55_noCut = new TH1F("hData_Jet55_noCut","data from matched jets without any jet ID cut",1000,0,1000);
  TH1F * hData_Jet55_CutA = new TH1F("hData_Jet55_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",1000,0,1000);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet55_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet55_CutA->Fill()}  
  TH1F * hData_Jet55_CutB = new TH1F("hData_Jet55_CutB","data from matched jets with Jet ID cut: slant line calopt/pfpt = 0.5, eMax/pfpt = 0 to calopt/pfpt = 1.25, eMax/pfpt = 1.1",1000,0,1000);
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet55_CutB->Fill()} 
  TH1F * hData_Jet55_CutA_rej = new TH1F("hData_Jet55_CutA_rej","",1000,0,1000);
  TH1F * hData_Jet55_CutB_rej = new TH1F("hData_Jet55_CutB_rej","",1000,0,1000);

  TH1F * hData_Jet65_noCut = new TH1F("hData_Jet65_noCut","data from matched jets without any jet ID cut",1000,0,1000);
  TH1F * hData_Jet65_CutA = new TH1F("hData_Jet65_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",1000,0,1000);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet65_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet65_CutA->Fill()}  
  TH1F * hData_Jet65_CutB = new TH1F("hData_Jet65_CutB","data from matched jets with Jet ID cut: slant line calopt/pfpt = 0.5, eMax/pfpt = 0 to calopt/pfpt = 1.25, eMax/pfpt = 1.1",1000,0,1000);
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet65_CutB->Fill()} 
  TH1F * hData_Jet65_CutA_rej = new TH1F("hData_Jet65_CutA_rej","",1000,0,1000);
  TH1F * hData_Jet65_CutB_rej = new TH1F("hData_Jet65_CutB_rej","",1000,0,1000);
  
  TH1F * hData_Jet80_noCut = new TH1F("hData_Jet80_noCut","data from matched jets without any jet ID cut",1000,0,1000);
  TH1F * hData_Jet80_CutA = new TH1F("hData_Jet80_CutA","data from matched jets with Jet ID cut: slant line from 0.5 calopt/pfpt from eMax/Sumcand 0 till 1 and then calopt/pfpt > 0.85",1000,0,1000);
  // this cut is a combination of two: 
  // if(eMax/SumCand < 0.9 && ( eMax/Sumcand < (18/7 calopt/pfpt - 9/7) ) ) { hData_Jet80_CutA->Fill()}
  // if { eMax/Sumcand >0.9 && calopt/pfpt > 0.85} { hData_Jet80_CutA->Fill()}  
  TH1F * hData_Jet80_CutB = new TH1F("hData_Jet80_CutB","data from matched jets with Jet ID cut: slant line calopt/pfpt = 0.5, eMax/pfpt = 0 to calopt/pfpt = 1.25, eMax/pfpt = 1.1",1000,0,1000);
  // this cut is if( eMax/pfpt < (22/15 * calopt/pfpt - 11/15) ) {hData_Jet80_CutB->Fill()} 
  TH1F * hData_Jet80_CutA_rej = new TH1F("hData_Jet80_CutA_rej","",1000,0,1000);
  TH1F * hData_Jet80_CutB_rej = new TH1F("hData_Jet80_CutB_rej","",1000,0,1000);

  TH1F * hData_unmatched_Jet80_noCut = new TH1F("hData_unmatched_Jet80_noCut","",1000,0,1000);
  TH1F * hData_unmatched_Jet80_CutA = new TH1F("hData_unmatched_Jet80_CutA","",1000,0,1000);
  TH1F * hData_unmatched_Jet80_CutA_rej = new TH1F("hData_unmatched_Jet80_CutA_rej","",1000,0,1000);
  TH1F * hData_unmatched_Jet65_noCut = new TH1F("hData_unmatched_Jet65_noCut","",1000,0,1000);
  TH1F * hData_unmatched_Jet65_CutA = new TH1F("hData_unmatched_Jet65_CutA","",1000,0,1000);
  TH1F * hData_unmatched_Jet65_CutA_rej = new TH1F("hData_unmatched_Jet65_CutA_rej","",1000,0,1000);
  TH1F * hData_unmatched_Jet55_noCut = new TH1F("hData_unmatched_Jet55_noCut","",1000,0,1000);
  TH1F * hData_unmatched_Jet55_CutA = new TH1F("hData_unmatched_Jet55_CutA","",1000,0,1000);
  TH1F * hData_unmatched_Jet55_CutA_rej = new TH1F("hData_unmatched_Jet55_CutA_rej","",1000,0,1000);

  TH1F * hMC_unmatched_Jet80_noCut = new TH1F("hMC_unmatched_Jet80_noCut","",1000,0,1000);
  TH1F * hMC_unmatched_Jet80_CutA = new TH1F("hMC_unmatched_Jet80_CutA","",1000,0,1000);
  TH1F * hMC_unmatched_Jet80_CutA_rej = new TH1F("hMC_unmatched_Jet80_CutA_rej","",1000,0,1000);
  TH1F * hMC_unmatched_Jet65_noCut = new TH1F("hMC_unmatched_Jet65_noCut","",1000,0,1000);
  TH1F * hMC_unmatched_Jet65_CutA = new TH1F("hMC_unmatched_Jet65_CutA","",1000,0,1000);
  TH1F * hMC_unmatched_Jet65_CutA_rej = new TH1F("hMC_unmatched_Jet65_CutA_rej","",1000,0,1000);
  TH1F * hMC_unmatched_Jet55_noCut = new TH1F("hMC_unmatched_Jet55_noCut","",1000,0,1000);
  TH1F * hMC_unmatched_Jet55_CutA = new TH1F("hMC_unmatched_Jet55_CutA","",1000,0,1000);
  TH1F * hMC_unmatched_Jet55_CutA_rej = new TH1F("hMC_unmatched_Jet55_CutA_rej","",1000,0,1000);
  
  
  // 1 - Data, 2 - MC
  Float_t pfpt_1, pfpt_2;
  Float_t pfrefpt_2;
  Float_t calopt_1, calopt_2;
  Float_t eMax_1, eMax_2;
  Float_t chMax_1, chMax_2;
  Float_t chSum_1, chSum_2;
  Float_t phSum_1, phSum_2;
  Float_t neSum_1, neSum_2;
  Float_t muSum_1, muSum_2;
  Int_t jet55_1, jet65_1, jet80_1;
  Int_t jet55_p_1, jet65_p_1, jet80_p_1;
  Int_t jet55_2, jet65_2, jet80_2;
  Int_t jet55_p_2, jet65_p_2, jet80_p_2;

  Data_matched->SetBranchAddress("calopt",&calopt_1);
  Data_matched->SetBranchAddress("pfpt",&pfpt_1);
  Data_matched->SetBranchAddress("eMax",&eMax_1);
  Data_matched->SetBranchAddress("chMax",&chMax_1);
  Data_matched->SetBranchAddress("chSum",&chSum_1);
  Data_matched->SetBranchAddress("phSum",&phSum_1);
  Data_matched->SetBranchAddress("neSum",&neSum_1);
  Data_matched->SetBranchAddress("muSum",&muSum_1);
  Data_matched->SetBranchAddress("jet55",&jet55_1);
  Data_matched->SetBranchAddress("jet65",&jet65_1);
  Data_matched->SetBranchAddress("jet80",&jet80_1);
  Data_matched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_matched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_matched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  
  Data_unmatched->SetBranchAddress("pfpt",&pfpt_1);
  Data_unmatched->SetBranchAddress("eMax",&eMax_1);
  Data_unmatched->SetBranchAddress("chMax",&chMax_1);
  Data_unmatched->SetBranchAddress("chSum",&chSum_1);
  Data_unmatched->SetBranchAddress("phSum",&phSum_1);
  Data_unmatched->SetBranchAddress("neSum",&neSum_1);
  Data_unmatched->SetBranchAddress("muSum",&muSum_1);
  Data_unmatched->SetBranchAddress("jet55",&jet55_1);
  Data_unmatched->SetBranchAddress("jet65",&jet65_1);
  Data_unmatched->SetBranchAddress("jet80",&jet80_1);
  Data_unmatched->SetBranchAddress("jet55_prescl",&jet55_p_1);
  Data_unmatched->SetBranchAddress("jet65_prescl",&jet65_p_1);
  Data_unmatched->SetBranchAddress("jet80_prescl",&jet80_p_1);
  
  MC_matched->SetBranchAddress("calopt",&calopt_2);
  MC_matched->SetBranchAddress("pfpt",&pfpt_2);
  MC_matched->SetBranchAddress("eMax",&eMax_2);
  MC_matched->SetBranchAddress("chMax",&chMax_2);
  MC_matched->SetBranchAddress("chSum",&chSum_2);
  MC_matched->SetBranchAddress("phSum",&phSum_2);
  MC_matched->SetBranchAddress("neSum",&neSum_2);
  MC_matched->SetBranchAddress("muSum",&muSum_2);
  MC_matched->SetBranchAddress("pfrefpt",&pfrefpt_2);
  MC_matched->SetBranchAddress("jet55",&jet55_2);
  MC_matched->SetBranchAddress("jet65",&jet65_2);
  MC_matched->SetBranchAddress("jet80",&jet80_2);
  
  MC_unmatched->SetBranchAddress("pfpt",&pfpt_2);
  MC_unmatched->SetBranchAddress("eMax",&eMax_2);
  MC_unmatched->SetBranchAddress("chMax",&chMax_2);
  MC_unmatched->SetBranchAddress("chSum",&chSum_2);
  MC_unmatched->SetBranchAddress("phSum",&phSum_2);
  MC_unmatched->SetBranchAddress("neSum",&neSum_2);
  MC_unmatched->SetBranchAddress("muSum",&muSum_2);
  MC_unmatched->SetBranchAddress("pfrefpt",&pfrefpt_2);
  MC_unmatched->SetBranchAddress("jet55",&jet55_2);
  MC_unmatched->SetBranchAddress("jet65",&jet65_2);
  MC_unmatched->SetBranchAddress("jet80",&jet80_2);
  
  
  // data loop
  long entries = Data_matched->GetEntries();
  //entries = 1;
  cout<<"matched Data ntuple "<<endl;
  
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    Data_matched->GetEntry(nentry);
    
    Float_t Sumcand = chSum_1 + phSum_1 + neSum_1 + muSum_1;
    
    if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {
      
      hData_Jet55_noCut->Fill(pfpt_1, jet55_p_1);
      
      if(eMax_1/Sumcand < 0.9 && ( eMax_1/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) && chMax_1/pfpt_1 > 0.02  )
	hData_Jet55_CutA->Fill(pfpt_1, jet55_p_1);
      if(eMax_1/Sumcand >=0.9 && calopt_1/pfpt_1 > 0.85 && chMax_1/pfpt_1 > 0.02) 
	hData_Jet55_CutA->Fill(pfpt_1, jet55_p_1);
      
      // if(eMax_1/pfpt_1 < (22/15 * (Float_t)calopt_1/pfpt_1 - 11/15)) 
      // 	hData_Jet55_CutB->Fill(pfpt_1);

      if(eMax_1/Sumcand < 0.9 && ( eMax_1/Sumcand > (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) && chMax_1/pfpt_1 < 0.02 )
	hData_Jet55_CutA_rej->Fill(pfpt_1, jet55_p_1);
      if(eMax_1/Sumcand >0.9 && calopt_1/pfpt_1 < 0.85 && chMax_1/pfpt_1 < 0.02)
	hData_Jet55_CutA_rej->Fill(pfpt_1, jet55_p_1);
      
      // if(eMax_1/pfpt_1 > (22/15 * (Float_t)calopt_1/pfpt_1 - 11/15)) 
      // 	hData_Jet55_CutB_rej->Fill(pfpt_1);
      
    }
    
    if(jet65_1 == 1 && jet80_1 == 0 ) {
      
      hData_Jet65_noCut->Fill(pfpt_1);
      
      if(eMax_1/Sumcand < 0.9 && ( eMax_1/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) && chMax_1/pfpt_1 > 0.02 )
	hData_Jet65_CutA->Fill(pfpt_1);
      if(eMax_1/Sumcand >=0.9 && calopt_1/pfpt_1 > 0.85 && chMax_1/pfpt_1 > 0.02 ) 
	hData_Jet65_CutA->Fill(pfpt_1);
      
      // if(eMax_1/pfpt_1 < (22/15 * (Float_t)calopt_1/pfpt_1 - 11/15)) 
      // 	hData_Jet65_CutB->Fill(pfpt_1);
    
      if(eMax_1/Sumcand < 0.9 && ( eMax_1/Sumcand > (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7)) && chMax_1/pfpt_1 < 0.02 )
	hData_Jet65_CutA_rej->Fill(pfpt_1);
      if(eMax_1/Sumcand >0.9 && calopt_1/pfpt_1 < 0.85 && chMax_1/pfpt_1 < 0.02 )
	hData_Jet65_CutA_rej->Fill(pfpt_1);

      // if(eMax_1/pfpt_1 > (22/15 * (Float_t)calopt_1/pfpt_1 - 11/15)) 
      // 	hData_Jet65_CutB_rej->Fill(pfpt_1);

    }

    if(jet80_1 == 1) {
    
      hData_Jet80_noCut->Fill(pfpt_1);

      if(eMax_1/Sumcand < 0.9 && ( eMax_1/Sumcand < (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7) ) && chMax_1/pfpt_1 > 0.02 )
	hData_Jet80_CutA->Fill(pfpt_1);
      if(eMax_1/Sumcand >=0.9 && calopt_1/pfpt_1 > 0.85 && chMax_1/pfpt_1 > 0.02 ) 
	hData_Jet80_CutA->Fill(pfpt_1);

      // if(eMax_1/pfpt_1 < (22/15 * (Float_t)calopt_1/pfpt_1 - 11/15)) 
      // 	hData_Jet80_CutB->Fill(pfpt_1);
    
      if(eMax_1/Sumcand < 0.9 && ( eMax_1/Sumcand > (18/7 *(Float_t)calopt_1/pfpt_1 - 9/7) )&& chMax_1/pfpt_1 < 0.02 )
	hData_Jet80_CutA_rej->Fill(pfpt_1);
      if(eMax_1/Sumcand >0.9 && calopt_1/pfpt_1 < 0.85 && chMax_1/pfpt_1 < 0.02 )
	hData_Jet80_CutA_rej->Fill(pfpt_1);

      // if(eMax_1/pfpt_1 > (22/15 * (Float_t)calopt_1/pfpt_1 - 11/15)) 
      // 	hData_Jet80_CutB_rej->Fill(pfpt_1);
      
    }
    
  }// data ntuple loop

  // data unmatched loop:
  entries = Data_unmatched->GetEntries();
  //entries = 1;
  cout<<"Unmatched Data ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry ){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    Data_unmatched->GetEntry(nentry);
    
    Float_t Sumcand = chSum_1 + phSum_1 + neSum_1 + muSum_1;

    if(jet55_1 == 1 && jet65_1 == 0 && jet80_1 == 0 ) {
    
      hData_unmatched_Jet55_noCut->Fill(pfpt_1, jet55_p_1);

      if(eMax_1/Sumcand < 0.05 && chMax_1/pfpt_1 > 0.02 ) hData_unmatched_Jet55_CutA->Fill(pfpt_1, jet55_p_1);
      else hData_unmatched_Jet55_CutA_rej->Fill(pfpt_1, jet55_p_1);
      
    }

    if(jet65_1 == 1 && jet80_1 == 0 ) {

      hData_unmatched_Jet65_noCut->Fill(pfpt_1);

      if(eMax_1/Sumcand < 0.05 && chMax_1/pfpt_1 > 0.02 ) hData_unmatched_Jet65_CutA->Fill(pfpt_1);
      else hData_unmatched_Jet65_CutA_rej->Fill(pfpt_1);
      
    }

    if(jet80_1 == 1) {
    
      hData_unmatched_Jet80_noCut->Fill(pfpt_1);

      if(eMax_1/Sumcand < 0.05 && chMax_1/pfpt_1 > 0.02 ) hData_unmatched_Jet80_CutA->Fill(pfpt_1);
      else hData_unmatched_Jet80_CutA_rej->Fill(pfpt_1);
      
    }
    
  }// data ntuple loop



  entries = MC_matched->GetEntries();
  //entries = 1;
  // MC loop
  cout<<" looking at matched MC ntuple "<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    MC_matched->GetEntry(nentry);
    
    Float_t Sumcand = chSum_2 + phSum_2 + neSum_2 + muSum_2;

    if(jet55_2 == 1 && jet65_2==0 && jet80_2 == 0){
      
      hMC_Jet55_noCut->Fill(pfpt_2);

      if(eMax_2/Sumcand < 0.9 && ( eMax_2/Sumcand < (18/7 *(Float_t)calopt_2/pfpt_2 - 9/7) ) && chMax_2/pfpt_2 > 0.02 )
	hMC_Jet55_CutA->Fill(pfpt_2);
      if(eMax_2/Sumcand >0.9 && calopt_2/pfpt_2 > 0.85 && chMax_2/pfpt_2 > 0.02 ) 
	hMC_Jet55_CutA->Fill(pfpt_2);

      // if(eMax_2/pfpt_2 < (22/15 * (Float_t)calopt_2/pfpt_2 - 11/15)) 
      // 	hMC_Jet55_CutB->Fill(pfpt_2);

    }

    
    if(jet65_2 == 1 && jet80_2 == 0){
      
      hMC_Jet65_noCut->Fill(pfpt_2);

      if(eMax_2/Sumcand < 0.9 && ( eMax_2/Sumcand < (18/7 *(Float_t)calopt_2/pfpt_2 - 9/7) ) && chMax_2/pfpt_2 > 0.02 )
	hMC_Jet65_CutA->Fill(pfpt_2);
      if(eMax_2/Sumcand >0.9 && calopt_2/pfpt_2 > 0.85 && chMax_2/pfpt_2 > 0.02 ) 
	hMC_Jet65_CutA->Fill(pfpt_2);

      // if(eMax_2/pfpt_2 < (22/15 * (Float_t)calopt_2/pfpt_2 - 11/15)) 
      // 	hMC_Jet65_CutB->Fill(pfpt_2);

    }

    
    if(jet80_2 == 1){
      
      hMC_Jet80_noCut->Fill(pfpt_2);

      if(eMax_2/Sumcand < 0.9 && ( eMax_2/Sumcand < (18/7 *(Float_t)calopt_2/pfpt_2 - 9/7) ) && chMax_2/pfpt_2 > 0.02 )
	hMC_Jet80_CutA->Fill(pfpt_2);
      if(eMax_2/Sumcand >0.9 && calopt_2/pfpt_2 > 0.85 && chMax_2/pfpt_2 > 0.02 ) 
	hMC_Jet80_CutA->Fill(pfpt_2);

      // if(eMax_2/pfpt_2 < (22/15 * (Float_t)calopt_2/pfpt_2 - 11/15)) 
      // 	hMC_Jet80_CutB->Fill(pfpt_2);

    }
    
    
  }// mc ntuple loop


  entries = MC_unmatched->GetEntries();
  //entries = 1;
  // MC loop
  cout<<" looking at unmatched MC ntuple"<<endl;
  for(long nentry = 0; nentry < entries; ++nentry){

    if(nentry%10000 == 0) cout<<" nentry = "<<nentry<<endl;
    MC_unmatched->GetEntry(nentry);
    
    Float_t Sumcand = chSum_2 + phSum_2 + neSum_2 + muSum_2;

    if(jet55_2 == 1 && jet65_2==0 && jet80_2 == 0){
      
      hMC_unmatched_Jet55_noCut->Fill(pfpt_2);
      if(eMax_2/Sumcand < 0.05 && chMax_2/pfpt_2 > 0.02 ) hMC_unmatched_Jet55_CutA->Fill(pfpt_2);
      //else hMC_unmatched_Jet55_CutA_rej->Fill(pfpt_2);
      
    }

    
    if(jet65_2 == 1 && jet80_2 == 0){
      
      hMC_unmatched_Jet65_noCut->Fill(pfpt_2);
      if(eMax_2/Sumcand < 0.05 && chMax_2/pfpt_2 > 0.02 ) hMC_unmatched_Jet65_CutA->Fill(pfpt_2);
      //else hMC_unmatched_Jet65_CutA_rej->Fill(pfpt_2);
      
    }

    
    if(jet80_2 == 1){
      
      hMC_unmatched_Jet80_noCut->Fill(pfpt_2);
      if(eMax_2/Sumcand < 0.05 && chMax_2/pfpt_2 > 0.02 ) hMC_unmatched_Jet80_CutA->Fill(pfpt_2);
      //else hMC_unmatched_Jet80_CutA_rej->Fill(pfpt_2);
      
    }
    
    
  }// mc unmatched  ntuple loop

  


  TFile fout("PbPb_CutEfficiency_YetkinCuts_matched_slantedlinecalopfpt_addingunmatched_exclusionhighertriggers_eMaxSumcand_A_chMaxJtpt.root","RECREATE");

  // add the unmatched histograms to the matched ones to get the final cut efficiency
  hData_Jet55_noCut->Add(hData_unmatched_Jet55_noCut);
  hData_Jet65_noCut->Add(hData_unmatched_Jet65_noCut);
  hData_Jet80_noCut->Add(hData_unmatched_Jet80_noCut);
  
  hData_Jet55_CutA->Add(hData_unmatched_Jet55_CutA);
  hData_Jet65_CutA->Add(hData_unmatched_Jet65_CutA);
  hData_Jet80_CutA->Add(hData_unmatched_Jet80_CutA);

  hData_Jet55_CutA_rej->Add(hData_unmatched_Jet55_CutA_rej);
  hData_Jet65_CutA_rej->Add(hData_unmatched_Jet65_CutA_rej);
  hData_Jet80_CutA_rej->Add(hData_unmatched_Jet80_CutA_rej);

  hMC_Jet55_noCut->Add(hMC_unmatched_Jet55_noCut);
  hMC_Jet65_noCut->Add(hMC_unmatched_Jet65_noCut);
  hMC_Jet80_noCut->Add(hMC_unmatched_Jet80_noCut);
  
  hMC_Jet55_CutA->Add(hMC_unmatched_Jet55_CutA);
  hMC_Jet65_CutA->Add(hMC_unmatched_Jet65_CutA);
  hMC_Jet80_CutA->Add(hMC_unmatched_Jet80_CutA);

  //hMC_Jet55_CutA_rej->Add(hMC_unmatched_Jet55_CutA_rej);
  //hMC_Jet65_CutA_rej->Add(hMC_unmatched_Jet65_CutA_rej);
  //hMC_Jet80_CutA_rej->Add(hMC_unmatched_Jet80_CutA_rej);


  hData_Jet65_noCut->Write();
  hData_Jet65_CutA->Write();
  // hData_Jet65_CutB->Write();
  TH1F * hData_Jet65_CutA_eff = (TH1F*)hData_Jet65_CutA->Clone("hData_Jet65_CutA_eff");
  hData_Jet65_CutA_eff->Divide(hData_Jet65_noCut);
  hData_Jet65_CutA_eff->Write();
  // TH1F * hData_Jet65_CutB_eff = (TH1F*)hData_Jet65_CutB->Clone("hData_Jet65_CutB_eff");
  // hData_Jet65_CutB_eff->Divide(hData_Jet65_noCut);
  // hData_Jet65_CutB_eff->Write();
  hData_Jet65_CutA_rej->Write();
  // hData_Jet65_CutB_rej->Write();
  hData_Jet65_CutA_rej_eff = (TH1F*)hData_Jet65_CutA_rej->Clone("hData_Jet65_CutA_rej_eff");
  hData_Jet65_CutA_rej_eff->Divide(hData_Jet65_noCut);
  hData_Jet65_CutA_rej_eff->Write();
  // hData_Jet65_CutB_rej_eff = (TH1F*)hData_Jet65_CutB_rej->Clone("hData_Jet65_CutB_rej_eff");
  // hData_Jet65_CutB_rej_eff->Divide(hData_Jet65_noCut);
  // hData_Jet65_CutB_rej_eff->Write();

  hData_Jet55_noCut->Write();
  hData_Jet55_CutA->Write();
  // hData_Jet55_CutB->Write();
  TH1F * hData_Jet55_CutA_eff = (TH1F*)hData_Jet55_CutA->Clone("hData_Jet55_CutA_eff");
  hData_Jet55_CutA_eff->Divide(hData_Jet55_noCut);
  hData_Jet55_CutA_eff->Write();
  // TH1F * hData_Jet55_CutB_eff = (TH1F*)hData_Jet55_CutB->Clone("hData_Jet55_CutB_eff");
  // hData_Jet55_CutB_eff->Divide(hData_Jet55_noCut);
  // hData_Jet55_CutB_eff->Write();
  hData_Jet55_CutA_rej->Write();
  // hData_Jet55_CutB_rej->Write();
  hData_Jet55_CutA_rej_eff = (TH1F*)hData_Jet55_CutA_rej->Clone("hData_Jet55_CutA_rej_eff");
  hData_Jet55_CutA_rej_eff->Divide(hData_Jet55_noCut);
  hData_Jet55_CutA_rej_eff->Write();
  // hData_Jet55_CutB_rej_eff = (TH1F*)hData_Jet55_CutB_rej->Clone("hData_Jet55_CutB_rej_eff");
  // hData_Jet55_CutB_rej_eff->Divide(hData_Jet55_noCut);
  // hData_Jet55_CutB_rej_eff->Write();

  hData_Jet80_noCut->Write();
  hData_Jet80_CutA->Write();
  //hData_Jet80_CutB->Write();
  TH1F * hData_Jet80_CutA_eff = (TH1F*)hData_Jet80_CutA->Clone("hData_Jet80_CutA_eff");
  hData_Jet80_CutA_eff->Divide(hData_Jet80_noCut);
  hData_Jet80_CutA_eff->Write();
  // TH1F * hData_Jet80_CutB_eff = (TH1F*)hData_Jet80_CutB->Clone("hData_Jet80_CutB_eff");
  // hData_Jet80_CutB_eff->Divide(hData_Jet80_noCut);
  // hData_Jet80_CutB_eff->Write();
  hData_Jet80_CutA_rej->Write();
  // hData_Jet80_CutB_rej->Write();
  hData_Jet80_CutA_rej_eff = (TH1F*)hData_Jet80_CutA_rej->Clone("hData_Jet80_CutA_rej_eff");
  hData_Jet80_CutA_rej_eff->Divide(hData_Jet80_noCut);
  hData_Jet80_CutA_rej_eff->Write();
  // hData_Jet80_CutB_rej_eff = (TH1F*)hData_Jet80_CutB_rej->Clone("hData_Jet80_CutB_rej_eff");
  // hData_Jet80_CutB_rej_eff->Divide(hData_Jet80_noCut);
  // hData_Jet80_CutB_rej_eff->Write();
  
  
  hMC_Jet80_noCut->Write();
  hMC_Jet80_CutA->Write();
  // hMC_Jet80_CutB->Write();
  TH1F * hMC_Jet80_CutA_eff = (TH1F*)hMC_Jet80_CutA->Clone("hMC_Jet80_CutA_eff");
  hMC_Jet80_CutA_eff->Divide(hMC_Jet80_noCut);
  hMC_Jet80_CutA_eff->Write();
  // TH1F * hMC_Jet80_CutB_eff = (TH1F*)hMC_Jet80_CutB->Clone("hMC_Jet80_CutB_eff");
  // hMC_Jet80_CutB_eff->Divide(hMC_Jet80_noCut);
  // hMC_Jet80_CutB_eff->Write();
  
  hMC_Jet65_noCut->Write();
  hMC_Jet65_CutA->Write();
  // hMC_Jet65_CutB->Write();
  TH1F * hMC_Jet65_CutA_eff = (TH1F*)hMC_Jet65_CutA->Clone("hMC_Jet65_CutA_eff");
  hMC_Jet65_CutA_eff->Divide(hMC_Jet65_noCut);
  hMC_Jet65_CutA_eff->Write();
  // TH1F * hMC_Jet65_CutB_eff = (TH1F*)hMC_Jet65_CutB->Clone("hMC_Jet65_CutB_eff");
  // hMC_Jet65_CutB_eff->Divide(hMC_Jet65_noCut);
  // hMC_Jet65_CutB_eff->Write();
  
  hMC_Jet55_noCut->Write();
  hMC_Jet55_CutA->Write();
  // hMC_Jet55_CutB->Write();
  TH1F * hMC_Jet55_CutA_eff = (TH1F*)hMC_Jet55_CutA->Clone("hMC_Jet55_CutA_eff");
  hMC_Jet55_CutA_eff->Divide(hMC_Jet55_noCut);
  hMC_Jet55_CutA_eff->Write();
  // TH1F * hMC_Jet55_CutB_eff = (TH1F*)hMC_Jet55_CutB->Clone("hMC_Jet55_CutB_eff");
  // hMC_Jet55_CutB_eff->Divide(hMC_Jet55_noCut);
  // hMC_Jet55_CutB_eff->Write();

  // save the unmatched histograms as well:
  hData_unmatched_Jet80_noCut->Write();
  hData_unmatched_Jet80_CutA->Write();
  hData_unmatched_Jet80_CutA_rej->Write();
  hData_unmatched_Jet65_noCut->Write();
  hData_unmatched_Jet65_CutA->Write();
  hData_unmatched_Jet65_CutA_rej->Write();
  hData_unmatched_Jet55_noCut->Write();
  hData_unmatched_Jet55_CutA->Write();
  hData_unmatched_Jet55_CutA_rej->Write();

  hMC_unmatched_Jet80_noCut->Write();
  hMC_unmatched_Jet80_CutA->Write();
  hMC_unmatched_Jet80_CutA_rej->Write();
  hMC_unmatched_Jet65_noCut->Write();
  hMC_unmatched_Jet65_CutA->Write();
  hMC_unmatched_Jet65_CutA_rej->Write();
  hMC_unmatched_Jet55_noCut->Write();
  hMC_unmatched_Jet55_CutA->Write();
  hMC_unmatched_Jet55_CutA_rej->Write();
  
  
  TCanvas * cJet80_CutEfficiency_Jet80 = new TCanvas("cJet80_CutEfficiency_Jet80","",1000,800);
  cJet80_CutEfficiency_Jet80->Divide(2,1);
  cJet80_CutEfficiency_Jet80->cd(1);

  hData_Jet80_CutA_eff->Rebin(5);
  hData_Jet80_CutA_eff->Scale(1./5);
  hData_Jet80_CutA_eff->SetMarkerColor(kRed);
  hData_Jet80_CutA_eff->SetMarkerStyle(24);
  hData_Jet80_CutA_eff->SetAxisRange(20,600,"X");
  hData_Jet80_CutA_eff->SetAxisRange(0,1.2,"Y");
  hData_Jet80_CutA_eff->SetXTitle("matched akPu3PF p_{T}");
  hData_Jet80_CutA_eff->SetTitle("Data");
  hData_Jet80_CutA_eff->SetYTitle("Jet80_Cut efficiency");
  hData_Jet80_CutA_eff->Draw();
  // hData_Jet80_CutB_eff->Rebin(20);
  // hData_Jet80_CutB_eff->Scale(1./20);
  // hData_Jet80_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet80_CutB_eff->SetMarkerStyle(33);
  // hData_Jet80_CutB_eff->SetAxisRange(20,600,"X");
  // hData_Jet80_CutB_eff->Draw("same");

  cJet80_CutEfficiency_Jet80->cd(2);
  hMC_Jet80_CutA_eff->Rebin(5);
  hMC_Jet80_CutA_eff->Scale(1./5);
  hMC_Jet80_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet80_CutA_eff->SetMarkerStyle(24);
  hMC_Jet80_CutA_eff->SetAxisRange(20,600,"X");
  hMC_Jet80_CutA_eff->SetAxisRange(0,1.2,"Y");
  hMC_Jet80_CutA_eff->SetXTitle("matched akPu3PF p_{T}");
  hMC_Jet80_CutA_eff->SetTitle("MC");
  hMC_Jet80_CutA_eff->SetYTitle("Jet80_Cut efficiency");
  hMC_Jet80_CutA_eff->Draw();
  // hMC_Jet80_CutB_eff->Rebin(20);
  // hMC_Jet80_CutB_eff->Scale(1./20);
  // hMC_Jet80_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet80_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet80_CutB_eff->SetAxisRange(20,600,"X");
  // hMC_Jet80_CutB_eff->Draw("same");

  cJet80_CutEfficiency_Jet80->SaveAs("PbPb_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff.pdf","RECREATE");

  TCanvas * cCutRejection_Jet80 = new TCanvas("cCutRejection_Jet80","",1000,800);

  hData_Jet80_CutA_rej->Rebin(5);
  hData_Jet80_CutA_rej->Scale(1./5);
  hData_Jet80_CutA_rej->SetMarkerColor(kRed);
  hData_Jet80_CutA_rej->SetMarkerStyle(24);
  hData_Jet80_CutA_rej->SetAxisRange(20,600,"X");
  hData_Jet80_CutA_rej->SetAxisRange(0,1.2,"Y");
  hData_Jet80_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  hData_Jet80_CutA_rej->SetTitle("Data");
  hData_Jet80_CutA_rej->SetYTitle("Jet80_Cut Rejection");
  hData_Jet80_CutA_rej->Draw();
  // hData_Jet80_CutB_rej->Rebin(20);
  // hData_Jet80_CutB_rej->Scale(1./20);
  // hData_Jet80_CutB_rej->SetMarkerColor(kBlack);
  // hData_Jet80_CutB_rej->SetMarkerStyle(33);
  // hData_Jet80_CutB_rej->SetAxisRange(20,600,"X");
  // hData_Jet80_CutB_rej->Draw("same");

  cCutRejection_Jet80->SaveAs("PbPb_YetkinCuts_Jet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection.pdf","RECREATE");

  TCanvas * cJet55_CutEfficiency_Jet55 = new TCanvas("cJet55_CutEfficiency_Jet55","",1000,800);
  cJet55_CutEfficiency_Jet55->Divide(2,1);
  cJet55_CutEfficiency_Jet55->cd(1);

  hData_Jet55_CutA_eff->Rebin(5);
  hData_Jet55_CutA_eff->Scale(1./5);
  hData_Jet55_CutA_eff->SetMarkerColor(kRed);
  hData_Jet55_CutA_eff->SetMarkerStyle(24);
  hData_Jet55_CutA_eff->SetAxisRange(20,600,"X");
  hData_Jet55_CutA_eff->SetAxisRange(0,1.2,"Y");
  hData_Jet55_CutA_eff->SetXTitle("matched akPu3PF p_{T}");
  hData_Jet55_CutA_eff->SetTitle("Data");
  hData_Jet55_CutA_eff->SetYTitle("Jet55_Cut efficiency");
  hData_Jet55_CutA_eff->Draw();
  // hData_Jet55_CutB_eff->Rebin(20);
  // hData_Jet55_CutB_eff->Scale(1./20);
  // hData_Jet55_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet55_CutB_eff->SetMarkerStyle(33);
  // hData_Jet55_CutB_eff->SetAxisRange(20,600,"X");
  // hData_Jet55_CutB_eff->SetAxisRange(0,1.2,"Y");
  // hData_Jet55_CutB_eff->Draw("same");

  cJet55_CutEfficiency_Jet55->cd(2);
  hMC_Jet55_CutA_eff->Rebin(5);
  hMC_Jet55_CutA_eff->Scale(1./5);
  hMC_Jet55_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet55_CutA_eff->SetMarkerStyle(24);
  hMC_Jet55_CutA_eff->SetAxisRange(20,600,"X");
  hMC_Jet55_CutA_eff->SetAxisRange(0,1.2,"Y");
  hMC_Jet55_CutA_eff->SetXTitle("matched akPu3PF p_{T}");
  hMC_Jet55_CutA_eff->SetTitle("MC");
  hMC_Jet55_CutA_eff->SetYTitle("Jet55_Cut efficiency");
  hMC_Jet55_CutA_eff->Draw();
  // hMC_Jet55_CutB_eff->Rebin(20);
  // hMC_Jet55_CutB_eff->Scale(1./20);
  // hMC_Jet55_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet55_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet55_CutB_eff->SetAxisRange(20,600,"X");
  // hMC_Jet55_CutB_eff->Draw("same");

  cJet55_CutEfficiency_Jet55->SaveAs("PbPb_YetkinCuts_Jet55_noJet65noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched.pdf","RECREATE");

  TCanvas * cCutRejection_Jet55 = new TCanvas("cCutRejection_Jet55","",1000,800);

  hData_Jet55_CutA_rej->Rebin(5);
  hData_Jet55_CutA_rej->Scale(1./5);
  hData_Jet55_CutA_rej->SetMarkerColor(kRed);
  hData_Jet55_CutA_rej->SetMarkerStyle(24);
  hData_Jet55_CutA_rej->SetAxisRange(20,600,"X");
  hData_Jet55_CutA_rej->SetAxisRange(0,1.2,"Y");
  hData_Jet55_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  hData_Jet55_CutA_rej->SetTitle("Data");
  hData_Jet55_CutA_rej->SetYTitle("Jet55_Cut Rejection");
  hData_Jet55_CutA_rej->Draw();
  // hData_Jet55_CutB_rej->Rebin(20);
  // hData_Jet55_CutB_rej->Scale(1./20);
  // hData_Jet55_CutB_rej->SetMarkerColor(kBlack);
  // hData_Jet55_CutB_rej->SetMarkerStyle(33);
  // hData_Jet55_CutB_rej->SetAxisRange(20,600,"X");
  // hData_Jet55_CutB_rej->Draw("same");

  cCutRejection_Jet55->SaveAs("PbPb_YetkinCuts_Jet55_noJet65noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");

  TCanvas * cJet65_CutEfficiency_Jet65 = new TCanvas("cJet65_CutEfficiency_Jet65","",1000,800);
  cJet65_CutEfficiency_Jet65->Divide(2,1);
  cJet65_CutEfficiency_Jet65->cd(1);

  hData_Jet65_CutA_eff->Rebin(5);
  hData_Jet65_CutA_eff->Scale(1./5);
  hData_Jet65_CutA_eff->SetMarkerColor(kRed);
  hData_Jet65_CutA_eff->SetMarkerStyle(24);
  hData_Jet65_CutA_eff->SetAxisRange(20,600,"X");
  hData_Jet65_CutA_eff->SetAxisRange(0,1.2,"Y");
  hData_Jet65_CutA_eff->SetXTitle("matched akPu3PF p_{T}");
  hData_Jet65_CutA_eff->SetTitle("Data");
  hData_Jet65_CutA_eff->SetYTitle("Jet65_Cut efficiency");
  hData_Jet65_CutA_eff->Draw();
  // hData_Jet65_CutB_eff->Rebin(20);
  // hData_Jet65_CutB_eff->Scale(1./20);
  // hData_Jet65_CutB_eff->SetMarkerColor(kBlack);
  // hData_Jet65_CutB_eff->SetMarkerStyle(33);
  // hData_Jet65_CutB_eff->SetAxisRange(20,600,"X");
  // hData_Jet65_CutB_eff->Draw("same");

  cJet65_CutEfficiency_Jet65->cd(2);
  hMC_Jet65_CutA_eff->Rebin(5);
  hMC_Jet65_CutA_eff->Scale(1./5);
  hMC_Jet65_CutA_eff->SetMarkerColor(kRed);
  hMC_Jet65_CutA_eff->SetMarkerStyle(24);
  hMC_Jet65_CutA_eff->SetAxisRange(20,600,"X");
  hMC_Jet65_CutA_eff->SetAxisRange(0,1.2,"Y");
  hMC_Jet65_CutA_eff->SetXTitle("matched akPu3PF p_{T}");
  hMC_Jet65_CutA_eff->SetTitle("MC");
  hMC_Jet65_CutA_eff->SetYTitle("Jet65_Cut efficiency");
  hMC_Jet65_CutA_eff->Draw();
  // hMC_Jet65_CutB_eff->Rebin(20);
  // hMC_Jet65_CutB_eff->Scale(1./20);
  // hMC_Jet65_CutB_eff->SetMarkerColor(kBlack);
  // hMC_Jet65_CutB_eff->SetMarkerStyle(33);
  // hMC_Jet65_CutB_eff->SetAxisRange(20,600,"X");
  // hMC_Jet65_CutB_eff->Draw("same");

  cJet65_CutEfficiency_Jet65->SaveAs("PbPb_YetkinCuts_Jet65_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_Eff_including_unmatched.pdf","RECREATE");

  TCanvas * cCutRejection_Jet65 = new TCanvas("cCutRejection_Jet65","",1000,800);

  hData_Jet65_CutA_rej->Rebin(5);
  hData_Jet65_CutA_rej->Scale(1./5);
  hData_Jet65_CutA_rej->SetMarkerColor(kRed);
  hData_Jet65_CutA_rej->SetMarkerStyle(24);
  hData_Jet65_CutA_rej->SetAxisRange(20,600,"X");
  hData_Jet65_CutA_rej->SetAxisRange(0, 1.2,"Y");
  hData_Jet65_CutA_rej->SetXTitle("matched akPu3PF p_{T}");
  hData_Jet65_CutA_rej->SetTitle("Data");
  hData_Jet65_CutA_rej->SetYTitle("Jet65_Cut Rejection");
  hData_Jet65_CutA_rej->Draw();
  // hData_Jet65_CutB_rej->Rebin(20);
  // hData_Jet65_CutB_rej->Scale(1./20);
  // hData_Jet65_CutB_rej->SetMarkerColor(kBlack);
  // hData_Jet65_CutB_rej->SetMarkerStyle(33);
  // hData_Jet65_CutB_rej->SetAxisRange(20,600,"X");
  // hData_Jet65_CutB_rej->Draw("same");

  cCutRejection_Jet65->SaveAs("PbPb_YetkinCuts_Jet65_noJet80_eMaxSumcand_A_chMaxJtpt_calopfpt_rejection_including_unmatched.pdf","RECREATE");


  // plot the trigger combination from this, and the total cut efficiency:
  
  TH1F * hData_Jet80 = (TH1F*)hData_Jet80_CutA->Clone("hData_Jet80");
  TH1F * hData_Jet65 = (TH1F*)hData_Jet65_CutA->Clone("hData_Jet65");
  TH1F * hData_Jet55 = (TH1F*)hData_Jet55_CutA->Clone("hData_Jet55");

  TH1F * hData_Combined = (TH1F*)hData_Jet80->Clone("hData_Combined");
  hData_Combined->Add(hData_Jet65);
  hData_Combined->Add(hData_Jet55);

  TH1F * hData_noCut_Jet80 = (TH1F*)hData_Jet80_noCut->Clone("hData_noCut_Jet80");
  TH1F * hData_noCut_Jet65 = (TH1F*)hData_Jet65_noCut->Clone("hData_noCut_Jet65");
  TH1F * hData_noCut_Jet55 = (TH1F*)hData_Jet55_noCut->Clone("hData_noCut_Jet55");

  TH1F * hData_noCut_Combined = (TH1F*)hData_noCut_Jet80->Clone("hData_Combined");
  hData_noCut_Combined->Add(hData_noCut_Jet65);
  hData_noCut_Combined->Add(hData_noCut_Jet55);
  
  TH1F * hMC_Jet80 = (TH1F*)hMC_Jet80_CutA->Clone("hMC_Jet80");
  TH1F * hMC_Jet65 = (TH1F*)hMC_Jet65_CutA->Clone("hMC_Jet65");
  TH1F * hMC_Jet55 = (TH1F*)hMC_Jet55_CutA->Clone("hMC_Jet55");

  TH1F * hMC_Combined = (TH1F*)hMC_Jet80->Clone("hMC_Combined");
  hMC_Combined->Add(hMC_Jet65);
  hMC_Combined->Add(hMC_Jet55);
  
  TH1F * hMC_noCut_Jet80 = (TH1F*)hMC_Jet80_noCut->Clone("hMC_noCut_Jet80");
  TH1F * hMC_noCut_Jet65 = (TH1F*)hMC_Jet65_noCut->Clone("hMC_noCut_Jet65");
  TH1F * hMC_noCut_Jet55 = (TH1F*)hMC_Jet55_noCut->Clone("hMC_noCut_Jet55");

  TH1F * hMC_noCut_Combined = (TH1F*)hMC_noCut_Jet80->Clone("hMC_Combined");
  hMC_noCut_Combined->Add(hMC_noCut_Jet65);
  hMC_noCut_Combined->Add(hMC_noCut_Jet55);
  
  TH1F * hData_Combined_Efficiency = (TH1F*)hData_Combined->Clone("hData_Combined_Efficiency");
  hData_Combined_Efficiency->Divide(hData_noCut_Combined);
  TH1F * hMC_Combined_Efficiency = (TH1F*)hMC_Combined->Clone("hMC_Combined_Efficiency");
  hMC_Combined_Efficiency->Divide(hMC_noCut_Combined);

  TCanvas * cCombinedEff = new TCanvas("cCombinedEff","",800,600);
  hData_Combined_Efficiency->SetXTitle("Jet p_{T}");
  hData_Combined_Efficiency->SetYTitle("Combined Jet ID cut efficiency");
  hData_Combined_Efficiency->SetMarkerStyle(20);
  hData_Combined_Efficiency->SetMarkerColor(kBlack);
  hData_Combined_Efficiency->Rebin(10);
  hData_Combined_Efficiency->Scale(1/10);
  hData_Combined_Efficiency->SetAxisRange(30,350,"X");
  hData_Combined_Efficiency->Draw();

  hMC_Combined_Efficiency->SetMarkerStyle(24);
  hMC_Combined_Efficiency->SetMarkerColor(kRed);
  hMC_Combined_Efficiency->Rebin(10);
  hMC_Combined_Efficiency->Scale(1/10);
  hMC_Combined_Efficiency->Draw("same");

  cCombinedEff->SaveAs("Combined_trigger_efficiency_YetkinCut_chMaxJtpt.pdf","RECREATE");
  
  TCanvas * cTriggerCombination = new TCanvas("cTriggerCombination","",800,600);
  cTriggerCombination->SetLogy();

  hData_Combined->SetMarkerColor(kBlack);
  hData_Combined->SetMarkerStyle(25);
  hData_Combined->SetAxisRange(30,350,"X");
  hData_Combined->Draw();

  hData_Jet80->SetMarkerColor(kRed);
  hData_Jet80->SetMarkerStyle(20);
  hData_Jet80->Draw("same");
  
  hData_Jet65->SetMarkerColor(kBlue);
  hData_Jet65->SetMarkerStyle(20);
  hData_Jet65->Draw("same");

  hData_Jet55->SetMarkerColor(kGreen);
  hData_Jet55->SetMarkerStyle(20);
  hData_Jet55->Draw("same");

  cTriggerCombination->SaveAs("TriggerCombination_YetkinCuts_chMaxJtpt.pdf","RECREATE");


}
