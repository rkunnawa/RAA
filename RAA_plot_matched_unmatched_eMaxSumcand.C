{

  gSystem->Load("Headers/plot.h");

  TFile * fData = TFile::Open("../../Output/PbPb_Data_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150320.root");
  TTree * Data_matched = (TTree*)fData->Get("matchedJets");

  TFile * fMC = TFile::Open("../../Output/PbPb_MC_calo_pf_jet_correlation_deltaR_0p2_akPu3_20150319_9.root");
  TTree * MC_matched = (TTree*)fMC->Get("matchedJets");

  int ptSelection = 17;
  int ptBoundary[ptSelection+1] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
  
  TH2F* hData_eMaxSumcand_calopfpt_ptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_calopfpt_ptselection[ptSelection];
  TH2F* hMC_eMaxSumcand_calopfpt_refptselection[ptSelection];
  TCanvas * ceMaxSumcand [ptSelection];

  TLine * CutA_1 = new TLine(0.5, 0, 0.85, 0.9);
  CutA_1->SetLineStyle(2);
  CutA_1->SetLineWidth(2);
  CutA_1->SetLineColor(kRed);
  TLine * CutA_2 = new TLine(0.85, 0.9, 0.85, 1.5);
  CutA_2->SetLineStyle(2);
  CutA_2->SetLineWidth(2);
  CutA_2->SetLineColor(kRed);
  

  for(int a = 0;a<ptSelection; ++a) {

    hData_eMaxSumcand_calopfpt_ptselection[a] = new TH2F(Form("hData_eMaxSumcand_calopfpt_ptselection_%d",a),"",100, 0, 5, 100, 0, 5);
    hMC_eMaxSumcand_calopfpt_ptselection[a] = new TH2F(Form("hMC_eMaxSumcand_calopfpt_ptselection_%d",a),"",100, 0, 5, 100, 0, 5);
    hMC_eMaxSumcand_calopfpt_refptselection[a] = new TH2F(Form("hMC_eMaxSumcand_calopfpt_refptselection_%d",a),"",100, 0, 5, 100, 0, 5);

    Data_matched->Draw(Form("eMax/(pfpt-eMax):calopt/(pfpt-eMax)>>hData_eMaxSumcand_calopfpt_ptselection_%d",a),Form("pfpt > %d && pfpt < %d && jet55", ptBoundary[a], ptBoundary[a+1]),"goff");
    MC_matched->Draw(Form("eMax/(pfpt-eMax):calopt/(pfpt-eMax)>>hMC_eMaxSumcand_calopfpt_ptselection_%d",a),Form("pfpt > %d && pfpt < %d && jet55", ptBoundary[a], ptBoundary[a+1]),"goff");
    MC_matched->Draw(Form("eMax/(pfpt-eMax):calopt/(pfpt-eMax)>>hMC_eMaxSumcand_calopfpt_refptselection_%d",a),Form("pfpt > %d && pfpt < %d && pfrefpt > %d && pfrefpt < %d && jet55",  ptBoundary[a], ptBoundary[a+1],  ptBoundary[a], ptBoundary[a+1]),"goff");

    ceMaxSumcand[a] = new TCanvas(Form("ceMaxSumcand_%d",a),"",1000,800);
    ceMaxSumcand[a]->Divide(3,1);
    ceMaxSumcand[a]->cd(1);
    ceMaxSumcand[a]->cd(1)->SetLogz();
    hData_eMaxSumcand_calopfpt_ptselection[a]->SetXTitle("calopt/(pfpt-eMax)");
    hData_eMaxSumcand_calopfpt_ptselection[a]->SetYTitle("eMax/(pfpt-eMax)");
    hData_eMaxSumcand_calopfpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < pfpt < %d", ptBoundary[a], ptBoundary[a+1]),0.2,0.7,14);
    drawText("Data",0.25,0.8,14);
    CutA_1->Draw();
    CutA_2->Draw();
    
    ceMaxSumcand[a]->cd(2);
    ceMaxSumcand[a]->cd(2)->SetLogz();
    hMC_eMaxSumcand_calopfpt_ptselection[a]->SetXTitle("calopt/(pfpt-eMax)");
    hMC_eMaxSumcand_calopfpt_ptselection[a]->SetYTitle("eMax/(pfpt-eMax)");
    hMC_eMaxSumcand_calopfpt_ptselection[a]->Draw("colz");
    drawText(Form("%d < pfpt < %d", ptBoundary[a], ptBoundary[a+1]),0.2,0.7,14);
    drawText("MC",0.25,0.8,14);
    CutA_1->Draw();
    CutA_2->Draw();

    ceMaxSumcand[a]->cd(3);
    ceMaxSumcand[a]->cd(3)->SetLogz();
    hMC_eMaxSumcand_calopfpt_refptselection[a]->SetXTitle("calopt/(pfpt-eMax)");
    hMC_eMaxSumcand_calopfpt_refptselection[a]->SetYTitle("eMax/(pfpt-eMax)");
    hMC_eMaxSumcand_calopfpt_refptselection[a]->Draw("colz");
    drawText(Form("%d < pfpt, pfrefpt < %d", ptBoundary[a], ptBoundary[a+1]),0.2,0.7,14);
    drawText("MC",0.25,0.8,14);
    CutA_1->Draw();
    CutA_2->Draw();

    ceMaxSumcand[a]->SaveAs(Form("eMax_OverPfptMinus_eMax_calopfptMinuseMax_jet55_showingCut_%d_ptrange_%d.pdf",ptBoundary[a], ptBoundary[a+1]),"RECREATE");
    
  }

  





}
