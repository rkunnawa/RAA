//
// Example macro for SVD unfolding with RooUnfold
// Run macro in aliroot
// Macro needs as input root file which contains measured spectrum, response matrix, prior and regularization parameter beta
// Unfolding is done with different priors
// Result is compared to chi2 minimization
// Author: Marta Verweij (marta.verweij@cern.ch)
//

//define colors in which different centrality bins must be plotted
static  int      myDarkRed     = TColor::GetColor(128,0,0);
static  int      myDarkGreen   = TColor::GetColor(0,128,0);
static  int      myDarkBlue    = TColor::GetColor(0,0,128);

//-----------------------------------------------------------------------------------------------------------------------
void Load() {
  printf("Load libraries\n");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  // gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics");
}


//_______________________________________________________________________________________________________________________
TMatrixT<double> *CalculatePearsonCoefficients(TMatrixT<double> *covmat) {

  TMatrixT<double> *pearsonCoefs = (TMatrixT<double>*)covmat->Clone("pearsonCoefs");
  //  pearsonCoefs->Clear();

  Int_t nrows = covmat->GetNrows();
  Int_t ncols = covmat->GetNcols();

  Double_t pearson = 0.;

  for(int row = 0; row<nrows; row++) {
    for(int col = 0; col<ncols; col++) {

      pearson = covmat(row,col)/TMath::Sqrt(covmat(row,row)*covmat(col,col));
      //      cout << "(" << row << "," << col << ") = " << pearson << endl;
      pearsonCoefs(row,col) = pearson;
    }
  }

  return pearsonCoefs;

}

//_______________________________________________________________________________________________________________________
void doSVDUnfold(TString strInput = "input/UnfoldingHists.root", Int_t kregDraw = 3, Bool_t useToy = kTRUE) {

  //if useToy==kTRUE, toy experiments are performed to calculate covariance matrix -> adviced for SVD

  Load();
  gStyle->SetOptStat(0);
  // gStyle->SetOptFit(01);
  gStyle->SetOptTitle(0);

  gSystem->Load("$ROOUNFOLD/libRooUnfold.so"); //adjust to your local path

  const Int_t nKregMin = 2;
  const Int_t nKregMax = 7;
  if(kregDraw>nKregMax) nKregMax=kregDraw+1;

  //RooUnfold::ErrorTreatment
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
  if(useToy)  errorTreatment = RooUnfold::kCovToy;
  cout << "errorTreatment: " << errorTreatment << endl;

  Double_t nEvents = 1e4*1e3;

  TFile *fInput = new TFile(strInput.Data());
  TH1D *hUnfoldedChi2 = fInput->Get("hUnfolded");
  TH1D *hFoldedChi2 = fInput->Get("hFolded");
  TH1D *fh1RawJetsCurrent = fInput->Get("fh1RawJetsCurrent");
  TH1D *fh1RawJetsCurrentUnfBinning = fInput->Get("fh1RawJetsCurrentUnfBinning");

  TH1D *fh1RawJetsCurrentUnfBinningOrig =  (TH1D*)fh1RawJetsCurrentUnfBinning->Clone("fh1RawJetsCurrentUnfBinningOrig");

  //Get response matrix: x-axis: generated/true; y-axis: reconstructed/measured
  TH2D *h2ResponseMatrix =  (TH2D*)fInput->Get("h2ResponseMatrixCombinedOrig");
  if(!h2ResponseMatrix) h2ResponseMatrix = (TH2D*)fInput->Get("h2ResponseMatrixOrig");
  TH2D *h2ResponseMatrixNorm = (TH2D*)fInput->Get("h2ResponseMatrixCombined");
  if(!h2ResponseMatrixNorm) h2ResponseMatrixNorm = (TH2D*)fInput->Get("h2ResponseMatrix");
  TH1D *hEfficiency = (TH1D*)fInput->Get("hEfficiencyCombined");
  if(!hEfficiency) hEfficiency = (TH1D*)fInput->Get("hEfficiency");

  //ptmin = minimum pT of measured spectrum
  //ptmax = maximum pt of measured spectrum
  //range of unfolded spectrum given by x axis of response matrix
  Double_t ptmin = h2ResponseMatrixOrig->GetYaxis()->GetXmin();
  Double_t ptmax = h2ResponseMatrixOrig->GetYaxis()->GetXmax();

  //Smooth tail measured spectrum with power law fit -> used as measured prior
  TF1 *fPower = new TF1("fPower","[0]*TMath::Power(x,-([1]))",0.,200.);
  fh1RawJetsCurrentUnfBinning->Fit(fPower,"0","",60.,110.);

  for(int bin=1; bin<=fh1RawJetsCurrentUnfBinning->GetNbinsX(); bin++) {

    if(fh1RawJetsCurrentUnfBinning->GetBinCenter(bin)>110.) {

      fh1RawJetsCurrentUnfBinning->SetBinContent(bin,fPower->Integral(fh1RawJetsCurrentUnfBinning->GetXaxis()->GetBinLowEdge(bin),fh1RawJetsCurrentUnfBinning->GetXaxis()->GetBinUpEdge(bin))/fh1RawJetsCurrentUnfBinning->GetXaxis()->GetBinWidth(bin));
    }
   
  }

  //Plot priors
  TCanvas *c31 = new TCanvas("c31","c31 Priors",800,800);
  gPad->SetLogy();

  hUnfoldedChi2->SetLineColor(myDarkGreen);
  hUnfoldedChi2->SetLineWidth(3);
  hUnfoldedChi2->SetMarkerStyle(23);
  hUnfoldedChi2->SetMarkerColor(hUnfoldedChi2->GetLineColor());
  hUnfoldedChi2->SetLineWidth(2);
  hUnfoldedChi2->DrawCopy();

  fh1RawJetsCurrentUnfBinning->SetLineColor(myDarkRed);
  fh1RawJetsCurrentUnfBinning->SetLineWidth(3);
  fh1RawJetsCurrentUnfBinning->SetMarkerStyle(20);
  fh1RawJetsCurrentUnfBinning->SetMarkerColor(myDarkRed);
  fh1RawJetsCurrentUnfBinning->SetLineWidth(2);
  fh1RawJetsCurrentUnfBinning->DrawCopy("same");

  fh1RawJetsCurrentUnfBinningOrig->SetLineColor(1);
  fh1RawJetsCurrentUnfBinningOrig->SetLineWidth(3);
  fh1RawJetsCurrentUnfBinningOrig->SetMarkerStyle(24);
  fh1RawJetsCurrentUnfBinningOrig->SetMarkerColor(1);
  fh1RawJetsCurrentUnfBinningOrig->SetLineWidth(2);
  fh1RawJetsCurrentUnfBinningOrig->DrawCopy("same");

  fPower->SetLineColor(4);
  fPower->Draw("same");

  TLegend *leg31 = new TLegend(0.52,0.65,0.88,0.88,"Priors");
  leg31->SetFillColor(10);
  leg31->SetFillStyle(0);
  leg31->SetBorderSize(0);
  leg31->SetTextFont(gStyle->GetTextFont());
  leg31->SetTextSize(gStyle->GetTextSize()*0.6);
  leg31->AddEntry(hUnfoldedChi2,"#chi^{2} unfolded","lp");
  leg31->AddEntry(fh1RawJetsCurrentUnfBinning,"measured","lp");
  leg31->Draw();

  //Swap x and y axis in response according to how RooUnfold wants it
  TH2 *hResponseMatrixTranspose = GetTransposeResponsMatrix(h2ResponseMatrix);
  TH2 *hResponseMatrixTransposeNorm = GetTransposeResponsMatrix(h2ResponseMatrixNorm);

  //Use measured as prior
  TH2 *hResponseMatrixTransposePrior =   (TH2*)hResponseMatrixTranspose->Clone("hResponseMatrixTransposePrior");
  hResponseMatrixTransposePrior = NormalizeResponseMatrixYaxisWithPrior(hResponseMatrixTransposePrior, fh1RawJetsCurrentUnfBinning);
  TH2 *hResponseMatrixTransposeNormPrior =   (TH2*)hResponseMatrixTransposeNorm->Clone("hResponseMatrixTransposeNormPrior");
  hResponseMatrixTransposeNormPrior = NormalizeResponseMatrixYaxisWithPrior(hResponseMatrixTransposeNormPrior, fh1RawJetsCurrentUnfBinning);

  //Use unfolded solution from Chi2 method as prior
  TH2 *hResponseMatrixTransposePriorChi2 =   (TH2*)hResponseMatrixTranspose->Clone("hResponseMatrixTransposePriorChi2");
  hResponseMatrixTransposePriorChi2 = NormalizeResponseMatrixYaxisWithPrior(hResponseMatrixTransposePriorChi2, hUnfoldedChi2);//
  RooUnfoldResponse responseChi2(0,0,hResponseMatrixTransposePriorChi2,"respCombinedChi2","respCombinedChi2");

  //Plot resonse matrix
  TCanvas *c1RMComb = new TCanvas("c1RMComb","c1RMComb:PtGen vs PtRec",800,800);
  c1RMComb->Divide(2,2);
  c1RMComb->cd(1);
  gPad->SetLogz();
  hResponseMatrixTranspose->Draw("colz");
  c1RMComb->cd(2);
  hEfficiency->Draw();
  c1RMComb->cd(3);
  gPad->SetLogz();
  hResponseMatrixTransposePrior->Draw("colz");

  //Initialize histos for unfolded, refolded and pearson coefficients
  TH1 *hUnfoldedSVDPriorMeas[nKregMax];
  TH1 *hFoldedSVDPriorMeas[nKregMax];
  TH2 *hPearsonSVDPriorMeas[nKregMax] = {0};

  TH1 *hUnfoldedSVDPriorChi2[nKregMax];
  TH1 *hFoldedSVDPriorChi2[nKregMax];
  TH2 *hPearsonSVDPriorChi2[nKregMax] = {0};

  TH1 *hUnfoldedSVDPriorTruth[nKregMax];
  TH1 *hFoldedSVDPriorTruth[nKregMax];
  TH2 *hPearsonSVDPriorTruth[nKregMax] = {0};

  Int_t nRows = Int_t(TMath::Sqrt(nKregMax))+1;
  TCanvas *cPearsonMatrixIter = new TCanvas("cPearsonMatrixIter","cPearsonMatrixIter",600,600);
  cPearsonMatrixIter->Divide(nRows,nRows);

  TCanvas *cPearsonMatrixIterPriorChi2 = new TCanvas("cPearsonMatrixIterPriorChi2","cPearsonMatrixIterPriorChi2",600,600);
  cPearsonMatrixIterPriorChi2->Divide(nRows,nRows);

  //Response for SVD Unfolding
  RooUnfoldResponse responseSVDMeas(0,0,hResponseMatrixTransposePrior,"respCombinedSVDMeas","respCombinedSVDMeas");
  RooUnfoldResponse responseSVDChi2(0,0,hResponseMatrixTransposePriorChi2,"respCombinedSVDChi2","respCombinedSVDChi2");
  
  // TDecompSVD *svd = unfoldSVD.Impl();
  PrintMatrix(responseSVDMeas.Mresponse(),"","response matrix",10);

  c1RMComb->cd(4);
  gPad->SetLogz();
  TH2 *h2ResponseRooUnfold = responseSVDMeas.Hresponse();
  h2ResponseRooUnfold->DrawCopy("colz"); 
  cout << "h2ResponseRooUnfold->GetNbinsX(): " <<  h2ResponseRooUnfold->GetNbinsX() << endl;
  cout << "h2ResponseRooUnfold->GetNbinsY(): " <<  h2ResponseRooUnfold->GetNbinsY() << endl;


  //_______________________________________________________________________________________________________________________
  // HERE THE UNFOLDING FINALLY STARTS
  // unfold with SVD method and store unfolded spectrum with different priors, pearson coefficients for different regularizations (kreg).
  // As prios used: measured spectrum, chi2 unfolded solution.
  //

  //loop for different regularization -> choose best one based on d-vector which is extracted later
  for(Int_t ikreg=nKregMin; ikreg<nKregMax; ikreg++) {

    RooUnfoldSvd      unfoldSVD(&responseSVDMeas, fh1RawJetsCurrent, ikreg);
    hUnfoldedSVDPriorMeas[ikreg-1] = (TH1D*) unfoldSVD.Hreco(errorTreatment);

    //Get covariance matrix and calculate corresponding Pearson coefficients    
    TMatrixD covmat = unfoldSVD.Ereco(errorTreatment);
    TMatrixD *pearson = (TMatrixD*)CalculatePearsonCoefficients(&covmat);
    if(pearson) {
      cout << "print pearson coefficients" << endl;      
      pearson->Print();

      hPearsonSVDPriorMeas[ikreg-1] = new TH2D(*pearson);
      
      gStyle->SetOptTitle(1);
      cPearsonMatrixIter->cd(ikreg);
      gPad->SetRightMargin(0.16);
      
      hPearsonSVDPriorMeas[ikreg-1]->SetTitle(Form("ikreg=%d",ikreg));
      hPearsonSVDPriorMeas[ikreg-1]->SetName("hPearsonSVDPriorMeas[ikreg-1]");
      hPearsonSVDPriorMeas[ikreg-1]->SetMinimum(-1.);
      hPearsonSVDPriorMeas[ikreg-1]->SetMaximum(1.);
      hPearsonSVDPriorMeas[ikreg-1]->GetZaxis()->SetLabelSize(0.06);
      hPearsonSVDPriorMeas[ikreg-1]->DrawCopy("colz");
    }

    //Correct for efficiency and get refolded distribution
    hUnfoldedSVDPriorMeas[ikreg-1]->Divide(hEfficiency);
    hFoldedSVDPriorMeas[ikreg-1] = MultiplyResponseGenerated(hUnfoldedSVDPriorMeas[ikreg-1],h2ResponseMatrixNorm,hEfficiency);    

    //Redo unfolding with unfolded solution from Chi2 method as prior
    RooUnfoldSvd      unfoldSVDChi2(&responseSVDChi2, fh1RawJetsCurrent, ikreg);
    hUnfoldedSVDPriorChi2[ikreg-1] = (TH1D*) unfoldSVDChi2.Hreco(errorTreatment);
    hFoldedSVDPriorChi2[ikreg-1] = responseSVDChi2.ApplyToTruth(hUnfoldedSVDPriorChi2[ikreg-1]);

    //extract Pearson coefficients
    TMatrixD covmat = unfoldSVD.Ereco(errorTreatment);
    TMatrixD *pearson = (TMatrixD*)CalculatePearsonCoefficients(&covmat);
    if(pearson) {
      hPearsonSVDPriorChi2[ikreg-1] = new TH2D(*pearson);

      gStyle->SetOptTitle(1);
      cPearsonMatrixIterPriorChi2->cd(ikreg);
      gPad->SetRightMargin(0.16);

      hPearsonSVDPriorChi2[ikreg-1]->SetTitle(Form("ikreg=%d",ikreg));
      hPearsonSVDPriorChi2[ikreg-1]->SetName("hPearsonSVDPriorChi2[ikreg-1]");
      hPearsonSVDPriorChi2[ikreg-1]->SetMinimum(-1.);
      hPearsonSVDPriorChi2[ikreg-1]->SetMaximum(1.);
      hPearsonSVDPriorChi2[ikreg-1]->GetZaxis()->SetLabelSize(0.06);
      hPearsonSVDPriorChi2[ikreg-1]->DrawCopy("colz");
    }

    //Correct for efficiency
    hUnfoldedSVDPriorChi2[ikreg-1]->Divide(hEfficiency);
  }
  

  //Get and print singular values and d_i vector, print also on screen
  //Note that these do not depend on the regularization. The opposite: they tell you which regularization to use!
  TSVDUnfold *svdUnfold = unfoldSVD.Impl();
  TH1 *hSVal = svdUnfold->GetSV();
  TH1D *hdi = svdUnfold->GetD();
  TH1D  *hSValClone = (TH1D*)hSVal->Clone("hSValClone");
  cout << "print singular values " << endl;
  for(int bin=1; bin<=hSVal->GetNbinsX(); bin++)
    cout << "bin: " << bin << "  SV: " << hSVal->GetBinContent(bin) << endl;

  cout << "print di vector " <<  endl;
  for(int bin=1; bin<=hdi->GetNbinsX(); bin++)
    cout << "i: " << bin << "  di: " << hdi->GetBinContent(bin) << endl;

  //_______________________________________________________________________________________________________________________
  // Make useful plots
  //

  TCanvas *c11=new TCanvas("c11","c11: Singular Values");
  c11->Divide(2,2);
  c11->cd(1);
  gPad->SetLogy();
  hSVal->SetXTitle("singular values");
  hSVal->DrawCopy();
  c11->cd(2);
  gPad->SetLogy();
  hdi->SetXTitle("|d_{i}^{kreg}|");
  hdi->DrawCopy();

  TCanvas *c2=new TCanvas("c2","c2: unfolding SVD regularizations");
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.15);
  TH1F *frame2 = gPad->DrawFrame(0.,1e-7,130.,10.);
  frame2->SetXTitle("p_{T} (GeV/c)");
  frame2->SetYTitle("counts/event");//dN/dp_{T}");
  frame2->GetXaxis()->SetTitleSize(0.06);
  frame2->GetYaxis()->SetTitleSize(0.06);
  frame2->GetXaxis()->SetTitleOffset(0.9);
  frame2->GetYaxis()->SetTitleOffset(1.0);

  gPad->SetLogy();
  fh1RawJetsCurrent->SetLineColor(2);
  fh1RawJetsCurrent->SetLineWidth(3);
  fh1RawJetsCurrent->SetMarkerStyle(20);
  fh1RawJetsCurrent->SetMarkerColor(fh1RawJetsCurrent->GetLineColor());
  fh1RawJetsCurrent->Draw("same");

  TLegend *leg2 = new TLegend(0.6,0.55,0.88,0.88,"SVD unfolding");
  leg2->SetFillColor(10);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(gStyle->GetTextFont());
  leg2->SetTextSize(gStyle->GetTextSize()*0.6);

  leg2->AddEntry(fh1RawJetsCurrent,"Measured","lp");

  for(int ikreg =nKregMin; ikreg<nKregMax; ikreg++) {
    if(ikreg!=kregDraw )hUnfoldedSVDPriorMeas[ikreg-1]->SetLineColor(ikreg);
    hUnfoldedSVDPriorMeas[ikreg-1]->Draw("same");
    leg2->AddEntry(hUnfoldedSVDPriorMeas[ikreg-1],Form("kreg %d",ikreg),"l");
  }
  leg2->Draw();

  TCanvas *cPearsonMatrix = new TCanvas("cPearsonMatrix","cPearsonMatrix",600,600);
  gPad->SetRightMargin(0.16);
  if(hPearsonSVDPriorMeas[kregDraw-1]) {
    hPearsonSVDPriorMeas[kregDraw-1]->SetTitle("hPearsonSVDPriorMeas");
    hPearsonSVDPriorMeas[kregDraw-1]->SetName("hPearsonSVDPriorMeas");
    hPearsonSVDPriorMeas[kregDraw-1]->SetMinimum(-1.);
    hPearsonSVDPriorMeas[kregDraw-1]->SetMaximum(1.);
    hPearsonSVDPriorMeas[kregDraw-1]->GetZaxis()->SetLabelSize(0.06);
    hPearsonSVDPriorMeas[kregDraw-1]->DrawCopy("colz");
  } 
  

  gStyle->SetOptTitle(0);

  cout << "done SVD unfolding" << endl;
  //_______________________________________________________________________________________________________________________
  
  //_______________________________________________________________________________________________________________________
  //
  //Plot distributions: unfolded, refolded spectra, pearson coefficients
  //

  TCanvas *c1=new TCanvas("c1","c1: unfolding");
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.15);
  TH1F *frame1 = gPad->DrawFrame(0.,1e-7,hUnfoldedSVDPriorMeas[kregDraw-1]->GetXaxis()->GetXmax(),10.);
  frame1->SetXTitle("p_{T} (GeV/c)");
  frame1->SetYTitle("counts/event");//dN/dp_{T}");
  frame1->GetXaxis()->SetTitleSize(0.06);
  frame1->GetYaxis()->SetTitleSize(0.06);
  frame1->GetXaxis()->SetTitleOffset(0.9);
  frame1->GetYaxis()->SetTitleOffset(1.0);

  gPad->SetLogy();
  fh1RawJetsCurrent->SetLineColor(2);
  fh1RawJetsCurrent->SetLineWidth(3);
  fh1RawJetsCurrent->SetMarkerStyle(20);
  fh1RawJetsCurrent->SetMarkerColor(fh1RawJetsCurrent->GetLineColor());
  fh1RawJetsCurrent->Draw("same");

  hUnfoldedSVDPriorMeas[kregDraw-1] ->SetLineColor(kOrange);
  hUnfoldedSVDPriorMeas[kregDraw-1] ->SetLineWidth(3);
  hUnfoldedSVDPriorMeas[kregDraw-1] ->SetMarkerStyle(22);
  hUnfoldedSVDPriorMeas[kregDraw-1] ->SetMarkerColor(hUnfoldedSVDPriorMeas[kregDraw-1]->GetLineColor());
  hUnfoldedSVDPriorMeas[kregDraw-1] ->SetLineWidth(2);
  hUnfoldedSVDPriorMeas[kregDraw-1] ->Draw("same");

  hFoldedSVDPriorMeas[kregDraw-1]->SetLineColor(kBlue-3);
  hFoldedSVDPriorMeas[kregDraw-1]->SetLineWidth(2);
  hFoldedSVDPriorMeas[kregDraw-1]->SetMarkerStyle(33);
  hFoldedSVDPriorMeas[kregDraw-1]->SetLineWidth(2);
  hFoldedSVDPriorMeas[kregDraw-1]->Draw("same");

  hFoldedChi2->SetLineColor(1);
  hFoldedChi2->SetLineWidth(2);
  hFoldedChi2->SetMarkerStyle(33);
  hFoldedChi2->Draw("same hist");

  hUnfoldedChi2->SetLineColor(myDarkGreen);
  hUnfoldedChi2->SetLineWidth(3);
  hUnfoldedChi2->SetMarkerStyle(23);
  hUnfoldedChi2->SetMarkerColor(hUnfoldedChi2->GetLineColor());
  hUnfoldedChi2->SetLineWidth(2);
  hUnfoldedChi2->Draw("same");

  TLegend *leg1 = new TLegend(0.6,0.6,0.88,0.88);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextFont(gStyle->GetTextFont());
  leg1->SetTextSize(gStyle->GetTextSize()*0.6);

  leg1->AddEntry(fh1RawJetsCurrent,"Measured","lp");
  leg1->AddEntry(hUnfoldedChi2,"Unfolded #chi^{2}","lp");
  leg1->AddEntry(hFoldedChi2,"Refolded #chi^{2}","l");

  leg1->AddEntry(hUnfoldedSVDPriorMeas[kregDraw-1],Form("Unfolded SVD kreg=%d",kregDraw),"lp");
  leg1->AddEntry(hFoldedSVDPriorMeas[kregDraw-1],Form("Refolded SVD kreg=%d",kregDraw),"l");

  leg1->Draw();

  TCanvas *c35 = new TCanvas("c35","c35:Ratio SVD/chi2",800,600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.15);
  TH1F *frame35 = gPad->DrawFrame(ptmin,0.,ptmax,2.);
  frame35->SetXTitle("p_{T} (GeV/c)");
  frame35->SetYTitle("SVD/#chi^{2}");
  frame35->GetXaxis()->SetTitleSize(0.06);
  frame35->GetYaxis()->SetTitleSize(0.06);
  frame35->GetXaxis()->SetTitleOffset(0.9);
  frame35->GetYaxis()->SetTitleOffset(1.0);

  TGraphErrors *grSVDPriorMeasChi2 = divide_histos(hUnfoldedSVDPriorMeas[kregDraw-1],hUnfoldedChi2);
  grSVDPriorMeasChi2->SetMarkerStyle(20);
  grSVDPriorMeasChi2->SetMarkerColor(6);
  grSVDPriorMeasChi2->Draw("pz");

  TGraphErrors *grSVDPriorChi2Chi2 = divide_histos(hUnfoldedSVDPriorChi2[kregDraw-1],hUnfoldedChi2);
  grSVDPriorChi2Chi2->SetMarkerStyle(21);
  grSVDPriorChi2Chi2->SetMarkerColor(2);
  grSVDPriorChi2Chi2->Draw("pz");

  TLegend *leg35 = new TLegend(0.6,0.7,0.88,0.88);
  leg35->SetFillColor(10);
  leg35->SetFillStyle(0);
  leg35->SetBorderSize(0);
  leg35->SetTextFont(gStyle->GetTextFont());
  leg35->SetTextSize(gStyle->GetTextSize()*0.6);

  leg35->AddEntry(grSVDPriorMeasChi2,"prior=measured","pl");
  leg35->AddEntry(grSVDPriorChi2Chi2,"prior=#chi^{2} solution","pl");

  leg35->Draw();

  TCanvas *c36 = new TCanvas("c36","c36:Ratio SVDPriorMeasd/chi2 iterations",800,600);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.15);
  TH1F *frame36 = gPad->DrawFrame(ptmin,0.,ptmax,2.);
  frame36->SetXTitle("p_{T} (GeV/c)");
  frame36->SetYTitle("SVD/#chi^{2}");
  frame36->GetXaxis()->SetTitleSize(0.06);
  frame36->GetYaxis()->SetTitleSize(0.06);
  frame36->GetXaxis()->SetTitleOffset(0.9);
  frame36->GetYaxis()->SetTitleOffset(1.0);

  TGraphErrors *grSVDPriorMeasChi2Iter[nKregMax];

  TLegend *leg36 = new TLegend(0.6,0.7,0.88,0.88);
  leg36->SetNColumns(2);
  leg36->SetFillColor(10);
  leg36->SetFillStyle(0);
  leg36->SetBorderSize(0);
  leg36->SetTextFont(gStyle->GetTextFont());
  leg36->SetTextSize(gStyle->GetTextSize()*0.6);

  for(int ikreg =nKregMin; ikreg<nKregMax; ikreg++) {
    grSVDPriorMeasChi2Iter[ikreg-1] = divide_histos(hUnfoldedSVDPriorMeas[ikreg-1],hUnfoldedChi2);
    grSVDPriorMeasChi2Iter[ikreg-1]->SetMarkerStyle(20);
    grSVDPriorMeasChi2Iter[ikreg-1]->SetMarkerColor(hUnfoldedSVDPriorMeas[ikreg-1]->GetLineColor());
    grSVDPriorMeasChi2Iter[ikreg-1]->SetLineColor(hUnfoldedSVDPriorMeas[ikreg-1]->GetLineColor());
    grSVDPriorMeasChi2Iter[ikreg-1]->Draw("pz");
    leg36->AddEntry(grSVDPriorMeasChi2Iter[ikreg-1],Form("k_{reg} %d",ikreg),"pl");
  }
  leg36->Draw();

  TFile *histos = new TFile("SVDSpectra.root","RECREATE");
  hdi->Write("hdi");
  hSVal->Write("hSVal");

  for(int iter =nKregMin; iter<nKregMax; iter++) {
    hUnfoldedSVDPriorMeas[iter-1]->Write(Form("hUnfoldedSVDPriorMeasKReg%d",iter));
    hFoldedSVDPriorMeas[iter-1]->Write(Form("hFoldedSVDPriorMeasKReg%d",iter));
    if(hPearsonSVDPriorMeas[iter-1]) hPearsonSVDPriorMeas[iter-1]->Write(Form("hPearsonSVDPriorMeasKReg%d",iter));

    hUnfoldedSVDPriorChi2[iter-1]->Write(Form("hUnfoldedSVDPriorChi2KReg%d",iter));
    if(hPearsonSVDPriorChi2[iter-1]) hPearsonSVDPriorChi2[iter-1]->Write(Form("hPearsonSVDPriorChi2KReg%d",iter));

    grSVDPriorMeasChi2Iter[iter-1]->Write(Form("grSVDPriorMeasChi2KReg%d",iter));
  }
  hUnfoldedChi2->Write("hUnfoldedChi2");

  histos->Write();
  histos->Close();
}

//-----------------------------------------------------------------------------------------------------------------------
TGraphErrors* divide_histos(TH1 *h1 = 0x0, TH1* h2 = 0x0, Double_t xmax=-1., Bool_t bNoErrorh2 = kFALSE) {

  TGraphErrors *gr = new TGraphErrors();

  float binCent = 0.;
  int j = 0;
  float ratio = 0.;
  float error2 = 0.;
  float binWidth = 0.;
  for(int i=1; i<=h1->GetNbinsX(); i++) {
    binCent = h1->GetXaxis()->GetBinCenter(i);
    if(xmax>0. && binCent>xmax) continue;
    j = h2->FindBin(binCent);
    binWidth = h1->GetXaxis()->GetBinWidth(i);
    if(h2->GetBinContent(j)>0.) {
      ratio = h1->GetBinContent(i)/h2->GetBinContent(j);
      
      Double_t A = 1./h2->GetBinContent(j)*h1->GetBinError(i);
      Double_t B = 0.;
      if(!bNoErrorh2) {
	if(h2->GetBinError(j)>0.) {
	  B = -1.*h1->GetBinContent(i)/(h2->GetBinContent(j)*h2->GetBinContent(j))*h2->GetBinError(j);
	  error2 = A*A + B*B;
	}
	else error2 = A*A;
      }
      else error2 = A*A;

      gr->SetPoint(gr->GetN(),binCent,ratio);
      gr->SetPointError(gr->GetN()-1,0.5*binWidth,TMath::Sqrt(error2));
    }
  }

  return gr;
} 

//-----------------------------------------------------------------------------------------------------------------------
PrintMatrix(const TMatrixD& m, const char* format,
	    const char* name, Int_t cols_per_sheet)
{
  // Print the matrix as a table of elements.
  // Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
  // By default the format "%11.4g" is used to print one element.
  // One can specify an alternative format with eg
  //  format ="%6.2f  "

  if (!m.IsValid()) {
    m.Error("PrintMatrix","%s is invalid",name);
    return;
  }

  const Int_t ncols  = m.GetNcols();
  const Int_t nrows  = m.GetNrows();
  const Int_t collwb = m.GetColLwb();
  const Int_t rowlwb = m.GetRowLwb();

  if (!(format && format[0])) format= "%11.4g ";
  char topbar[1000];
  snprintf(topbar,1000,format,123.456789);
  Int_t nch = strlen(topbar)+1;
  if (nch > 18) nch = 18;
  char ftopbar[20];
  for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
  Int_t nk = 1 + Int_t(log10(ncols));
  snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
  Int_t nch2 = strlen(ftopbar);
  for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
  ftopbar[nch] = '|';
  ftopbar[nch+1] = 0;

  printf("\n%dx%d %s is as follows",nrows,ncols,name);

  if (cols_per_sheet <= 0) {
    cols_per_sheet = 5;
    if (nch <= 8) cols_per_sheet =10;
  }
  nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
  for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
  topbar[nk] = 0;
  for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
    printf("\n\n     |");
    for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
      printf(ftopbar,j+collwb-1);
    printf("\n%s\n",topbar);
    if (m.GetNoElements() <= 0) continue;
    for (Int_t i = 1; i <= nrows; i++) {
      printf("%4d |",i+rowlwb-1);
      for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
	printf(format,m(i+rowlwb-1,j+collwb-1));

      printf("\n");
    }
  }
  printf("\n");
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* GetTransposeResponsMatrix(TH2 *h2RM) {
  //
  // Transpose response matrix
  //

  Double_t *ybinsArray = new Double_t[h2RM->GetNbinsY()+1];
  for(int i=1; i<=h2RM->GetNbinsY(); i++) {
    ybinsArray[i-1] = h2RM->GetYaxis()->GetBinLowEdge(i);
    Printf("i=%d  %f   %f",i,ybinsArray[i-1],h2RM->GetYaxis()->GetBinLowEdge(i));
  }
  ybinsArray[h2RM->GetNbinsY()] = h2RM->GetYaxis()->GetBinUpEdge(h2RM->GetNbinsY());

  Double_t *xbinsArray = new Double_t[h2RM->GetNbinsX()+1];
  for(int i=1; i<=h2RM->GetNbinsX(); i++) 
    xbinsArray[i-1] = h2RM->GetXaxis()->GetBinLowEdge(i);
  xbinsArray[h2RM->GetNbinsX()] = h2RM->GetXaxis()->GetBinUpEdge(h2RM->GetNbinsX());

  //Fill tranposed response matrix
  TH2D *h2RMTranspose = new TH2D("h2RMTranspose","h2RMTranspose",h2RM->GetNbinsY(),ybinsArray,h2RM->GetNbinsX(),xbinsArray);
  for (Int_t ibin = 1; ibin <= h2RMTranspose->GetNbinsX(); ibin++) {
    for (Int_t jbin = 1; jbin <= h2RMTranspose->GetNbinsY(); jbin++) {
      h2RMTranspose->SetBinContent(ibin,jbin, h2RM->GetBinContent(jbin,ibin));
    }
  }
  return h2RMTranspose;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* NormalizeResponseMatrixYaxisWithPrior(TH2 *h2RM, TH1 *hPrior) {
  //
  // Normalize such that the Y projection is the prior
  //

  double intPrior = hPrior->Integral();//"width");
  for (Int_t jbin = 1; jbin <= h2RM->GetNbinsY(); jbin++) {
       for (Int_t ibin = 1; ibin <= h2RM->GetNbinsX(); ibin++) {
	double content = h2RM->GetBinContent(ibin,jbin);
	h2RM->SetBinContent(ibin,jbin,hPrior->GetBinContent(jbin)/hPrior->GetBinWidth(jbin)/intPrior*content);
    }
  }
  return h2RM;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH1D* MultiplyResponseGenerated(TH1 *hGen, TH2 *hResponse,TH1 *hEfficiency,Bool_t bDrawSlices = kFALSE) {

  //
  // Multiply hGen with response matrix to obtain refolded spectrum
  // Efficiency must be given.
  //

  if(!hEfficiency) {
    printf("Setting efficiency to 1 \n");
    hEfficiency = (TH1D*)hGen->Clone("hEfficiency");
    hEfficiency->Reset();
    for(int i=1; i<=hEfficiency->GetNbinsX(); i++) hEfficiency->SetBinContent(i,1.);
  }

  //For response
  //x-axis: generated
  //y-axis: reconstructed

  TH1D *hRec = hResponse->ProjectionY("hRec");
  hRec->Sumw2();
  hRec->Reset();
  hRec->SetTitle("hRec");
  hRec->SetName("hRec");

  for(int irec=1; irec<=hRec->GetNbinsX(); irec++)
    hRec->SetBinContent(irec,0);

  TH1D *hRecGenBin = 0x0;
  TCanvas *cSlices = 0x0;
  if(bDrawSlices) {
    cSlices = new TCanvas("cSlices","cSlices:Slices",400,400);
    gPad->SetLogy();
  }

  Double_t yieldMC = 0.;
  Double_t yieldMCerror = 0.;
  Double_t sumYield = 0.;
  const Int_t nbinsRec = hRec->GetNbinsX()+1;
  Double_t sumError2[999] = {0.};
  Double_t eff = 0.;

  for(int igen=1; igen<=hGen->GetNbinsX(); igen++) {
    //get pTMC
    sumYield = 0.;
    eff = hEfficiency->GetBinContent(igen);
    yieldMC = hGen->GetBinContent(igen)*eff;
    yieldMCerror = hGen->GetBinError(igen)*eff;

    if(bDrawSlices) {
      hRecGenBin = hResponse->ProjectionY(Form("hRecGenBin%d",igen));
      hRecGenBin->Sumw2();
      hRecGenBin->Reset();
      hRecGenBin->SetTitle(Form("hRecGenBin%d",igen));
      hRecGenBin->SetName(Form("hRecGenBin%d",igen));
    }

    for(int irec=1; irec<=hRec->GetNbinsX(); irec++) {
      hRec->AddBinContent(irec,yieldMC*hResponse->GetBinContent(igen,irec));
      sumYield+=hResponse->GetBinContent(igen,irec);
      Double_t B = hResponse->GetBinContent(igen,irec)*yieldMCerror;
      sumError2[irec-1] += B*B;

      if(bDrawSlices)
	hRecGenBin->SetBinContent(irec,yieldMC*hResponse->GetBinContent(igen,irec));
    }
    if(bDrawSlices) {
      cSlices->cd();
      hRecGenBin->SetLineColor(igen);
      if(igen==1) hRecGenBin->DrawCopy();      
      else hRecGenBin->DrawCopy("same");
    }
    if(hRecGenBin) delete hRecGenBin;
    //    cout << "igen: " << igen << "\tpTMC: " << hGen->GetXaxis()->GetBinCenter(igen) << "\teff:" << eff << "\tsumYield: " << sumYield << endl;
  }
  
  for(int i=0; i<nbinsRec; i++) {
    if(sumError2[i]>0.)
      hRec->SetBinError(i+1,TMath::Sqrt(sumError2[i]));
  }
  return hRec;
}

