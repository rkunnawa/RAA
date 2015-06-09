{

  // For centrality with no other FWLite object

  //gSystem->Load("libDataFormatsHeavyIonEvent");
  //gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  //gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
//gROOT->ProcessLine(".x betterColors.C");	

   gStyle->SetErrorX(0);
   gStyle->SetPalette(1,0);
   gStyle->SetPadColor(0);
   gStyle->SetPadBorderSize(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameLineColor(0);
   gStyle->SetTitleColor(0);
   gStyle->SetTitleBorderSize(0); 
   gStyle->SetPalette(1,0); 
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetPadColor(0);
   gStyle->SetPadBorderSize(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameLineColor(0);
   gStyle->SetTextFont(42);//Def 62
   gStyle->SetLabelFont(42,"XYZ");
   gStyle->SetTitleFont(42,"XYZ");
   gStyle->SetTitleColor(kBlue);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetTitleXSize(.055);
   gStyle->SetTitleYSize(.055);
   gStyle->SetTitleXOffset(1.);
   gStyle->SetTitleYOffset(1.5);
   gStyle->SetLabelSize(0.045,"XYZ");
   gStyle->SetLabelOffset(0.01,"X");
   gStyle->SetLabelOffset(0.01,"Y");
   gStyle->SetTitleColor(1,"XYZ");
   gStyle->SetHistFillColor(1);
   gStyle->SetHistFillStyle(0);
   gStyle->SetHistLineColor(1);
   gStyle->SetHistLineStyle(0);
   gStyle->SetHistLineWidth(3);
   gStyle->SetHistLineWidth(1);
   gStyle->SetEndErrorSize(0);
   gStyle->SetErrorX(0);  
   gStyle->SetMarkerStyle(20);
   //gStyle->SetMarkerSize(1.25);
   gStyle->SetMarkerSize(1.2);
   gStyle->SetOptFit(1111);
   gStyle->SetStatColor(0);
   gStyle->SetStatBorderSize(1);
   gStyle->SetOptTitle(0);//if zero get rid of title
   gStyle->SetTitleXOffset(1);
   gStyle->SetTitleFillColor(0);
   gStyle->SetOptStat(0);
   gStyle->SetPadLeftMargin(0.17);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadTopMargin(0.03);
   gStyle->SetPadRightMargin(0.15);
   gROOT->ForceStyle();

   TText *tx= new TLatex(100,
			 100,
			 "CMS Preliminary");
   tx->SetTextAlign(22);
//tx->Draw();

}

