#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <TVector3.h> 
#include <stdio.h>
#include <string.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TLegend.h>
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

#define PI 3.14159;



///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////
void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, 
			  const Int_t rows, const Float_t leftOffset=0.,
                          const Float_t bottomOffset=0., 
			  const Float_t leftMargin=0.2, 
			  const Float_t bottomMargin=0.2,
                          const Float_t edge=0.05);

void drawText(const char *text, float xp, float yp);
void drawDum(float min, float max, double drawXLabel);
TH1D *make1DplotHIforest(TFile *fileIN, const char* file, const char* xVar, int& xBins ,double& xMin, double& xMax, TCut xSel, TCut genSel, int& pad);

//--------------------------------------------------------------
// drawPatch() is a crazy way of removing 0 in the second and third 
// pad which is partially shown due to no margin between the pads
// if anybody has a better way of doing it let me know! - Andre
//--------------------------------------------------------------
void drawPatch(float x1, float y1, float x2, float y2); 
//---------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////
//         TOOL BOX
//////////////////////////////////////////////////////////////////////
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
  if (canv==0) {
    Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];
  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
		      (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			(1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
			   Xlow[i],Ylow[j],Xup[i],Yup[j]);
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
      else pad[i][j]->SetLeftMargin(0);

      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);

      if(j==0) pad[i][j]->SetTopMargin(edge);
      else pad[i][j]->SetTopMargin(0);

      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);

      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}

void rescaleBins(TH1& h){
  for (int i =1; i<=h.GetNbinsX(); i++){//Skip bin=0 since it it is the underflow bin
    double oldBin = h.GetBinContent(i);
    double oldErr = h.GetBinError(i);
    //cout<<"i : "<<i<<" low edge: "<<h.GetBinLowEdge(i)<<" up edge: "<<(h.GetBinLowEdge(i)+h.GetBinWidth(i))<<" width: "<<h.GetBinWidth(i)<<" old content: "<<oldBin<<"  new content : "<<oldBin/(h.GetBinWidth(i))<<endl;
    h.SetBinContent(i,oldBin/(h.GetBinWidth(i)));
    h.SetBinError(i,oldErr/(h.GetBinWidth(i)));
  }
  return;

}
void drawText(const char *text, float xp, float yp){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(17);
  //tex->SetTextSize(0.05);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}

void format1Dhisto(TH1& h1, double Ymax, double Ymin, double& col, double& Mstyle, double& fill, double& style, const char* titx, const char* tity ){
  //void format1Dhisto(TH1& h1, string& xTitle, double Ymax, double Ymin){

  //h1.SetTitle(";XXXX;XXXX");
  if(Ymax!=-1 && Ymin!=-1) h1.GetYaxis()->SetRangeUser(Ymin, Ymax);
  //if(Ymax==-1 && Ymin!=-1) h1.GetYaxis()->SetMinimum(Ymin);
  h1.SetMarkerColor(col);
  h1.SetMarkerStyle(Mstyle);
  h1.SetLineColor(col);
  h1.SetFillColor(fill);
  h1.SetFillStyle(style);
  h1.SetMarkerSize(0.8);
  h1.GetXaxis()->SetTitle(titx);
  h1.GetYaxis()->SetTitle(tity);
  h1.GetXaxis()->CenterTitle();
  h1.GetYaxis()->CenterTitle();
  //cout<<"The title is : "<<tit<<endl;

  return;
}

void format1Dhisto(TH2& h2, double Ymax, double Ymin, double& col, double& Mstyle, double& fill, double& style, const char* titx, const char* tity ){
  //void format1Dhisto(TH1& h1, string& xTitle, double Ymax, double Ymin){

  //h2.SetTitle(";XXXX;XXXX");
  if(Ymax!=-1 && Ymin!=-1) h2.GetYaxis()->SetRangeUser(Ymax, Ymin);
  //if(Ymax==-1 && Ymin!=-1) h2.GetYaxis()->SetMinimum(Ymin);
  //h2.SetMarkerColor(col);
  //h2.SetMarkerStyle(Mstyle);
  //h2.SetLineColor(col);
  //h2.SetFillColor(fill);
  //h2.SetFillStyle(style);
  h2.GetXaxis()->SetTitle(titx);
  h2.GetYaxis()->SetTitle(tity);
  h2.GetXaxis()->CenterTitle();
  h2.GetYaxis()->CenterTitle();
  //cout<<"The title is : "<<tit<<endl;

  return;
}


void format1Dhisto(TH1& h1, double Ymax, double Ymin, double& col, double& fill, double& style, const char* titx, const char* tity ){
  //void format1Dhisto(TH1& h1, string& xTitle, double Ymax, double Ymin){

  //h1.SetTitle(";XXXX;XXXX");
  if(Ymax!=-1 && Ymin!=-1) h1.GetYaxis()->SetRangeUser(Ymax, Ymin);
  //if(Ymax==-1 && Ymin!=-1) h1.GetYaxis()->SetMinimum(Ymin);
  h1.SetMarkerColor(col);
  h1.SetMarkerStyle();
  h1.SetMarkerColor();
  h1.SetLineColor(col);
  h1.SetFillColor(fill);
  h1.SetFillStyle(style);
  h1.GetXaxis()->SetTitle(titx);
  h1.GetYaxis()->SetTitle(tity);
  h1.GetXaxis()->CenterTitle();
  h1.GetYaxis()->CenterTitle();
  //cout<<"The title is : "<<tit<<endl;

  return;
}
void format1Dhisto(TH1& h1, double Ymax, double Ymin, double& col, const char* titx, const char* tity ){
  //void format1Dhisto(TH1& h1, string& xTitle, double Ymax, double Ymin){

  //h1.SetTitle(";XXXX;XXXX");
  if(Ymax!=-1 && Ymin!=-1) h1.GetYaxis()->SetRangeUser(Ymax, Ymin);
  //if(Ymax==-1 && Ymin!=-1) h1.GetYaxis()->SetMinimum(Ymin);
  h1.SetMarkerColor(col);
  h1.SetLineColor(col);
  h1.GetXaxis()->SetTitle(titx);
  h1.GetYaxis()->SetTitle(tity);
  //cout<<"The title is : "<<tit<<endl;

  return;
}
void format1Dhisto(TH1& h1, double Ymax, double Ymin, double& col){
  //void format1Dhisto(TH1& h1, string& xTitle, double Ymax, double Ymin){

  //h1.SetTitle(";XXXX;XXXX");
  if(Ymax!=-1 && Ymin!=-1) h1.GetYaxis()->SetRangeUser(Ymax, Ymin);
  //if(Ymax==-1 && Ymin!=-1) h1.GetYaxis()->SetMinimum(Ymin);
  h1.SetMarkerColor(col);
  h1.SetLineColor(col);

  return;
}
void scaleToBinWidth(TH1& h1){
  double nBins = h1.GetXaxis()->GetNbins();
  for (int iBin = h1.GetXaxis()->GetFirst(); iBin <= h1.GetXaxis()->GetLast(); iBin++){
    double iW = h1.GetBinWidth(iBin);
    double old = h1.GetBinContent(iBin);
    h1.SetBinContent(iBin,old/iW);
  }




}



void LoopCorrRandomCone_MacroV6(){

  gROOT->ProcessLine(".x jorge_rootlogon.C");

  //  gROOT->ProcessLine(".x betterColors.C");

  double trkPtCut = 0;//2,0
  int i_tkPtCut = (int)trkPtCut;
  const char* trackPtCut = Form("min track pT= %4.2f",trkPtCut);
  bool doPrint = true;
  bool doData = true;
  bool doMC = false;
  bool do_c2 = false;
  bool do_c1 = true;
  TTree *akTreeMC[4];
  TTree *akTreeData[4];
  TTree *akTreepAData[1];
  TTree *akTreePbPbData[1];
  Double_t sigmaarray[6][3];
  Double_t sigmaarrayerr[6][3];
  if(doMC){
    // TFile *ak3MCFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak3_pA_HYDJET.root"));
    // TFile *ak4MCFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak4_pA_HYDJET.root"));
    // TFile *ak5MCFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak5_pA_HYDJET.root"));

    TFile *ak3MCFile = TFile::Open(Form("test_randomcone_MC_akVs3PF_20141007.root"));
    TFile *ak4MCFile = TFile::Open(Form("test_randomcone_MC_akVs4PF_20141007.root"));
    TFile *ak5MCFile = TFile::Open(Form("test_randomcone_MC_akVs5PF_20141007.root"));
    
    akTreeMC[0]   = (TTree*)ak3MCFile->Get("nt");
    akTreeMC[1]   = (TTree*)ak4MCFile->Get("nt");
    akTreeMC[2]   = (TTree*)ak5MCFile->Get("nt");
    //akTreeMC[3]   = (TTree*)ak5MCFile->Get("nt");
  }
  if(doData){
    // TFile *ak3dataFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak3_pA_DATA.root"));
    // TFile *ak4dataFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak4_pA_DATA.root"));
    // TFile *ak5dataFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak5_pA_DATA.root"));

    TFile *ak3dataFile = TFile::Open(Form("test_randomcone_data_akVs3PF_20141007.root"));
    TFile *ak4dataFile = TFile::Open(Form("test_randomcone_data_akVs4PF_20141007.root"));
    TFile *ak5dataFile = TFile::Open(Form("test_randomcone_data_akVs5PF_20141007.root"));

    akTreeData[0]   = (TTree*)ak3dataFile->Get("nt");
    akTreeData[1]   = (TTree*)ak4dataFile->Get("nt");
    akTreeData[2]   = (TTree*)ak5dataFile->Get("nt");
    //akTreeData[3]   = (TTree*)ak5dataFile->Get("nt");
    //TFile *ak4pADataFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak4_pA_DATA.root"));
    //TFile *ak4PbPbDataFile = TFile::Open(Form("randomCones_TkpTCut0_ak4_PbPb_DATA.root"));
    //TFile *ak4PbPbDataFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/randomCones_TkpTCut0_ak4_PbPb_DATA.root"));
    //TFile *ak4PbPbDataFile = TFile::Open(Form("/mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3/dataAKSkimNtupleRandomConeRings_v4_TkpTCut0_ak3dataMB.root"));

    //akTreepAData[0]   = (TTree*)ak4pADataFile->Get("nt");
    //akTreePbPbData[0]   = (TTree*)ak4PbPbDataFile->Get("nt");
  }
  const char* centLabel[] = {  "0005"   ,     "0510"       ,     "1030"        ,    "3050"         ,     "5070"        ,     "7090"       };
  const char* cent[] =      {  "0-5%"   ,     "5-10%"      ,    "10-30%"       ,   "30-50%"        ,    "50-70%"       ,     "70-90%"     };
  const char* centCuts[] =  {  "bin<=5" , "bin>5 && bin<=15" , "bin>15 && bin<=55" , "bin>55 && bin<=95", "bin>95 && bin<=135","bin>135 && bin<=175" };
 
  //we might want to look at smaller centrality bins as well 

  //   const char* centLabel[] = {  "002.5"   ,     "2.505"       ,     "057.5"        ,    "7.510"         ,     "1012.5"        ,     "12.515"       };
  //const char* cent[] =      {  "0-2.5%"   ,     "2.5-5%"      ,    "5-7.5%"       ,   "7.5-10%"        ,    "10-12.5%"       ,     "12.5-15%"     };
  //const char* centCuts[] =  {  "bin<=0" , "bin>0 &&bin<=1" , "bin>1 &&bin<=2" , "bin>2 &&bin<=3", "bin>3 &&bin<=4","bin>4&&bin<=5" };


  //const char* centLabel[] = {  "002.5"   ,     "2.505"       ,     "057.5"        ,    "7.510"         ,     "1012.5"        ,     "12.515"       };
  //const char* cent[] =      {  "70-75%"   ,     "75-80%"      ,    "80-85%"       ,   "85-90%"        ,    "90-95%"       ,     "95-100%"     };
  //const char* centCuts[] =  {  "bin>28&&bin<=30" , "bin>30 &&bin<=32" , "bin>32 &&bin<=34" , "bin>34 &&bin<=36", "bin>36 &&bin<=38","bin>38&&bin<=40" };


  const int nCent = 6;
 
  cout<<"going to draw specific histograms"<<endl;

  if (do_c2){

    const int nAlgos = 1;
    TH1D * ranConepAData[nCent][nAlgos];
    TH1D * ranConePbPbData[nCent][nAlgos];

    TLegend* leg2[nCent];
    const char* var     = "ranPFsumEt";
    const char* varLabel[nAlgos] = { "R=0.4" };
    const char* hType[2] = {"pAData    ","PbPbData"};
   
    double marker [2] = {20,25};//,21,25};//22,26};
    double color  [2] = {1,2};//{1,2,4,8};//,2,4};
    double fill   [nAlgos] = {0};//0,0,0,0};//,0,0};
    double meanpAData[nCent][nAlgos];
    double meanErrpAData[nCent][nAlgos]; 
    double meanPbPbData[nCent][nAlgos];
    double meanErrPbPbData[nCent][nAlgos];
   
    const char* xTitle = "sumPF E_{T}";
    const char* yTitle = "normalized counts";
    double Ymax = 10;
    double Ymin = 0.01;
    int c2nBins = 20;
    double c2xMax = 5;
    double c2xMin = 0;
    double xAxisBins[21] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0,1.4,1.8,2.2,2.6, 3.0,3.4,3.8,4.2,4.6,5.0};

    TCanvas *c2 = new TCanvas("c2","c2",1000,800);
    makeMultiPanelCanvas(c2,3,2,0.0,0.0,0.2,0.15,0.07);
    //Loop over the 6 centrality sections
    for (int i=0; i<nCent; i++){
      cout<<"division "<<i<<endl;
      leg2[i] = new TLegend(0.41,0.70,0.73,0.90);//top right
      leg2[i]->SetFillColor(0);
      leg2[i]->SetTextSize(0.028);
      leg2[i]->SetBorderSize(0);
      c2->cd(i+1);
      c2->cd(i+1)->SetLogy();
      for (int ir = 0; ir<=nAlgos-1; ir++){
	cout<<"algo "<<ir<<endl;
	//ranConepAData[i][ir] = new TH1D(Form("sumPFch_%s_%s_pAData",centLabel[i],varLabel[ir]),Form("sumPFch_%s_%s_pAData",centLabel[i],varLabel[ir]),c2nBins,c2xMin,c2xMax);
	ranConepAData[i][ir] = new TH1D(Form("sumPFch_%s_%s_pAData",centLabel[i],varLabel[ir]),Form("sumPFch_%s_%s_pAData",centLabel[i],varLabel[ir]),c2nBins,xAxisBins);
	ranConepAData[i][ir]->Sumw2();
	akTreepAData[ir]->Draw(Form("%s>>sumPFch_%s_%s_pAData",var,centLabel[i],varLabel[ir]),Form("%s && %s!=0",centCuts[i],var),"goff");
	//akTreepAData[ir]->Draw(Form("%s>>sumPFch_%s_%s_pAData",var,centLabel[i],varLabel[ir]),Form("%s",centCuts[i]),"goff");
	double a = ranConepAData[i][ir]->Integral(1,c2nBins);
	if(a!=0) ranConepAData[i][ir]->Scale(1/(a));
	meanpAData[i][ir] = ranConepAData[i][ir]->GetMean();
	meanErrpAData[i][ir] = ranConepAData[i][ir]->GetMeanError();
	format1Dhisto(*ranConepAData[i][ir],-1,-1,color[0],marker[0],color[0],fill[ir],xTitle,yTitle);
	ranConepAData[i][ir]->SetMaximum(1);
	ranConepAData[i][ir]->SetMinimum(0.001);
       
	ranConePbPbData[i][ir] = new TH1D(Form("sumPFch_%s_%s_PbPbData",centLabel[i],varLabel[ir]),Form("sumPFch_%s_%s_PbPbData",centLabel[i],varLabel[ir]),c2nBins,c2xMin,c2xMax);
	ranConePbPbData[i][ir]->Sumw2();
	akTreePbPbData[ir]->Draw(Form("%s>>sumPFch_%s_%s_PbPbData",var,centLabel[i],varLabel[ir]),Form("%s && %s!=0",centCuts[i],var),"goff");
	//akTreePbPbData[ir]->Draw(Form("%s>>sumPFch_%s_%s_PbPbData",var,centLabel[i],varLabel[ir]),Form("%s",centCuts[i]),"goff");
	double b =ranConePbPbData[i][ir]->Integral(1,c2nBins); 
	if(b!=0)ranConePbPbData[i][ir]->Scale(1/(b));
	meanPbPbData[i][ir] = ranConePbPbData[i][ir]->GetMean();
	meanErrPbPbData[i][ir] = ranConePbPbData[i][ir]->GetMeanError();
	format1Dhisto(*ranConePbPbData[i][ir],-1,-1,color[1],marker[1],color[1],fill[ir],xTitle,yTitle);
	ranConePbPbData[i][ir]->SetMaximum(1);
	ranConePbPbData[i][ir]->SetMinimum(0.001);
      

	leg2[i]->AddEntry(ranConepAData[i][ir],Form("%s %s [ mean: %5.2f #pm %5.2f ]",varLabel[ir],hType[0],meanpAData[i][ir],meanErrpAData[i][ir]),"lp");
	leg2[i]->AddEntry(ranConePbPbData[i][ir],Form("%s %s [ mean: %5.2f #pm %5.2f ]",varLabel[ir],hType[1],meanPbPbData[i][ir],meanErrPbPbData[i][ir]),"lp");
 
	ranConePbPbData[i][ir]->Draw(""); 
	leg2[i]->Draw();
	drawText(cent[i], 0.23, 0.83);
	ranConepAData[i][ir]->Draw("same"); 
      }

    }
  }



  if (do_c1){

    cout<<"drawing C1"<<endl;

    const int nAlgos = 3;
    TLegend* leg1[nCent];
    TH1D * ranConeMC[nCent][nAlgos];
    TH1D * ranConeData[nCent][nAlgos];
   
    const char* var     = "ranPFsumEt";
    const char* varLabel[nAlgos] = { "R=0.3","R=0.4","R=0.5" };
    const char* hType[2] = {"MC    ","DataMB"};
   
    double marker [2] = {20,25};//,21,25};//22,26};
    double color  [nAlgos] = {1,2,4};//,2,4};
    double fill   [nAlgos] = {0,0,0};//,0,0};
    double meanMC[nCent][nAlgos];
    double meanErrMC[nCent][nAlgos]; 
    double meanData[nCent][nAlgos];
    double meanErrData[nCent][nAlgos];
   
    const char* xTitle = "sumPF E_{T}";
    const char* yTitle = "normalized counts";
    double Ymax = 10;
    double Ymin = 0.01;
    int c1nBins = 20;
    double c1xMax = 5;
    double c1xMin = 0;
    double c1FixMinX = 0.2;
  
    double xAxisBins[21] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0,1.4,1.8,2.2,2.6, 3.0,3.4,3.8,4.2,4.6,5.0};
    TCanvas *c1MC = new TCanvas("c1MC","c1MC",1000,800);
    makeMultiPanelCanvas(c1MC,3,2,0.0,0.0,0.2,0.15,0.07);
    //Loop over the 6 centrality sections
    for (int i=0; i<nCent; i++){
      cout<<"centrality = "<<i<<endl;
      leg1[i] = new TLegend(0.41,0.70,0.73,0.90);//top right
      leg1[i]->SetFillColor(0);
      leg1[i]->SetTextSize(0.028);
      leg1[i]->SetBorderSize(0);
      c1MC->cd(i+1);
      c1MC->cd(i+1)->SetLogy();
      for (int ir = 0; ir<=nAlgos-1; ir++){
	if(doMC){
	  cout<<"algo "<<ir<<endl;
	  //ranConeMC[i][ir] = new TH1D(Form("sumPFch_%s_%s_MC",centLabel[i],varLabel[ir]),Form("sumPFch_%s_%s_MC",centLabel[i],varLabel[ir]),c1nBins,c1xMin,c1xMax);
	  ranConeMC[i][ir] = new TH1D(Form("sumPFch_%s_%s_MC",centLabel[i],varLabel[ir]),Form("sumPFch_%s_%s_MC",centLabel[i],varLabel[ir]),c1nBins,xAxisBins);
	  ranConeMC[i][ir]->Sumw2();
	  akTreeMC[ir]->Draw(Form("%s>>sumPFch_%s_%s_MC",var,centLabel[i],varLabel[ir]),Form("%s && %s>%f",centCuts[i],var,c1FixMinX),"goff");
	  cout<<" check1: content bin[7]: "<<ranConeMC[i][ir]->GetBinContent(7)<<"   err bin[7]: "<<ranConeMC[i][ir]->GetBinError(i)<<" relative Error: "<<(ranConeMC[i][ir]->GetBinError(i))/(ranConeMC[i][ir]->GetBinContent(7))<<endl;
	  ranConeMC[i][ir]->Scale(1/(ranConeMC[i][ir]->Integral(1,c1nBins)));
	  cout<<" check2(after histo scale): content bin[7]: "<<ranConeMC[i][ir]->GetBinContent(7)<<"   err bin[7]: "<<ranConeMC[i][ir]->GetBinError(i)<<" relative Error: "<<(ranConeMC[i][ir]->GetBinError(i))/(ranConeMC[i][ir]->GetBinContent(7))<<endl;
	  rescaleBins(*ranConeMC[i][ir]);
	  cout<<" check3(after bin scale): content bin[7]: "<<ranConeMC[i][ir]->GetBinContent(7)<<"   err bin[7]: "<<ranConeMC[i][ir]->GetBinError(i)<<" relative Error: "<<(ranConeMC[i][ir]->GetBinError(i))/(ranConeMC[i][ir]->GetBinContent(7))<<endl;
	  meanMC[i][ir] = ranConeMC[i][ir]->GetMean();
	  meanErrMC[i][ir] = ranConeMC[i][ir]->GetMeanError();
	  format1Dhisto(*ranConeMC[i][ir],-1,-1,color[ir],marker[0],color[ir],fill[ir],xTitle,yTitle);
	  ranConeMC[i][ir]->SetMaximum(1);
	  ranConeMC[i][ir]->SetMinimum(0.03);
	}
	if(doData ){
	  //ranConeData[i][ir] = new TH1D(Form("sumPFch_%s_%s_data",centLabel[i],varLabel[ir]),Form("sumPFch_%s_%s_data",centLabel[i],varLabel[ir]),c1nBins,c1xMin,c1xMax);
	  ranConeData[i][ir] = new TH1D(Form("sumPFch_%s_%s_data",centLabel[i],varLabel[ir]),Form("sumPFch_%s_%s_data",centLabel[i],varLabel[ir]),c1nBins,xAxisBins);
	  ranConeData[i][ir]->Sumw2();
	  akTreeData[ir]->Draw(Form("%s>>sumPFch_%s_%s_data",var,centLabel[i],varLabel[ir]),Form("%s && %s>%f",centCuts[i],var,c1FixMinX),"goff");
	  ranConeData[i][ir]->Scale(1/(ranConeData[i][ir]->Integral(1,c1nBins)));
	  rescaleBins(*ranConeData[i][ir]);
	  meanData[i][ir] = ranConeData[i][ir]->GetMean();
	  meanErrData[i][ir] = ranConeData[i][ir]->GetMeanError();
	  format1Dhisto(*ranConeData[i][ir],-1,-1,color[ir],marker[1],color[ir],fill[ir],xTitle,yTitle);
	  ranConeData[i][ir]->SetMaximum(1);
	  ranConeData[i][ir]->SetMinimum(0.03);
	}
	if( doMC ){
	  TF1 *fMC=new TF1("fMC","[0]*exp(-0.5*((x-[1])/[2])^2)",c1FixMinX ,c1xMax);
	  fMC->SetParameters(0.3,meanMC[i][ir],10);
	  fMC->FixParameter(1,meanMC[i][ir]);
	  ranConeMC[i][ir]->Fit("fMC","0");
	  Double_t sigmafitMC=fMC->GetParameter(2);
	  cout<<"sigmafitMC  "<<sigmafitMC<<endl;
	  cout<<"i   "<<i<<" ir  "<<ir<<endl;
	 
	  if (ir==0) ranConeMC[i][ir]->Draw("");
	  if (ir!=0) ranConeMC[i][ir]->Draw("same");	
	  if (ir==0){ leg1[i]->AddEntry("",Form("min track pT cut: %2.1f",trkPtCut),"");}	  
	  leg1[i]->AddEntry(ranConeMC[i][ir],Form("%s %s [ mean: %5.2f #pm %5.2f ]",varLabel[ir],hType[0],meanMC[i][ir],meanErrMC[i][ir]),"lp");
	  leg1[i]->Draw();	      
	  drawText(cent[i], 0.23, 0.83);
	  if (doData ){
	    //	   TF1 *f1=new TF1("f1","gaus",4,150);
	    TF1 *f1=new TF1("f1","[0]*exp(-0.5*((x-[1])/[2])^2)",4,150);
	    f1->SetParameters(0.3,meanData[i][ir],10);
	    f1->FixParameter(1,meanData[i][ir]);
	   
	    ranConeData[i][ir]->Fit("f1","0");
	    Double_t sigmafit=f1->GetParameter(2);
	    cout<<"sigma  "<<sigmafit<<endl;
	    cout<<"i   "<<i<<" ir  "<<ir<<endl;
	    sigmaarray[i][ir]=sigmafit;
	    sigmaarrayerr[i][ir]=f1->GetParError(2);
	    cout<<"test "<<sigmaarray[i][ir]<<endl;
	    ranConeData[i][ir]->Draw("same");
	    TF1 *f2=ranConeData[i][ir]->GetFunction("f1");
	    f2->Draw("same");
	    leg1[i]->AddEntry(ranConeData[i][ir],Form("%s %s [ mean: %5.2f #pm %5.2f ]",varLabel[ir],hType[1],meanData[i][ir],meanErrData[i][ir]),"lp");
	    leg1[i]->Draw();	      
	    drawText(cent[i], 0.23, 0.83);
	 
	  }

	  /*	 if (ir==0) ranConeMC[i][ir]->Draw("same");
		 if (ir!=0) ranConeMC[i][ir]->Draw("same");	
		 if (ir==0){      leg1[i]->AddEntry("",Form("min track pT cut: %2.1f",trkPtCut),"");}
       
		 leg1[i]->AddEntry(ranConeMC[i][ir],Form("%s %s [ mean: %5.2f #pm %5.2f ]",varLabel[ir],hType[0],meanMC[i][ir],meanErrMC[i][ir]),"lp");
		 leg1[i]->Draw();	      
		 drawText(cent[i], 0.23, 0.83);
	  */
	}
      }
    }  
    if (doPrint) c1MC->Print("JetRaa_PFSumEt_PtCut1.pdf");
    if (doPrint) c1MC->Print("JetRaa_PFSumEt_PtCut1.root"); 

    for (int i=0; i<nCent; i++){

      for(int ir=0; ir<3; ir++){
	cout<<"cent  "<<i<<" algos  "<<ir<<"    "<<sigmaarray[i][ir]<<"    "<<sigmaarrayerr[i][ir]<<"\n"<<endl;
	//cout<<"cent  "<<i<<" algos  "<<ir<<"    "<<sigmaarray[i][ir]<<"\n"<<endl;
      }
     
    }
   

  }
  return;
 
}
