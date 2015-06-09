{
  TFile *_file1 = TFile::Open("PbPb_data_vz_cent_akPuPF_20150508_130.root");
  TFile *_file0 = TFile::Open("PbPbMinBias_data_vz_cent_akPuPF_20150508_130.root");

  hMinBiasbin=(TH1F*)_file0->Get("hBinEvents");
  hHLTbin=(TH1F*)_file1->Get("hBinEvents");
  hHLTbin->Draw();
  
  hratio=(TH1F*) hHLTbin->Clone();
  
  // hratio->Draw();
  hratio->Divide(hMinBiasbin);
  hratio->SetLineColor(2);
float integral= hratio->Integral(0,180);
 hratio->Scale(1./integral);
  hratio->Draw("same");


  TFile f("MinBiasHLTRatioWeightsVZ.root","RECREATE");
  f.cd();
  hratio->Write();
  f->Close();

  //hHLTbin->Divide(hMinBiasbin);



}
