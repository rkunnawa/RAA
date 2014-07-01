// Raghav Kunnawalkam Elayavalli
// June 11 2014
// CERN

// RAA - plotting macro. pretty much all the neccessary plots will be made here in this macro. 
// will take input from the read macro and the analysis macro. Maybe not from the read macro - i will have to decide on that. 

// June 25th - just started working on the macro. will plot the iteration systematics, MC closure test, create setup for Unfolded vs measured, RAA and Normalized Response matrix plots as well. 

// July 1st - finished the macro. it would have been very easy and efficient to put them all in 2 segments, 
//            one for PbPb (with the centrality loop) and one for pp. 
//          - But ive decided to keep it this way since we can easily delete any segments which are complete and remake any plots individually. 

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <fstream>
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

  
static const int nbins_pt = 29;
static const double boundaries_pt[nbins_pt+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}


TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname)
{
	TH1F *hF = (TH1F*)h->Clone(fHistname);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
		hF->SetBinContent(i,var);
		hF->SetBinError(i,0);
	}
	return hF;
}


// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}



void drawText(const char *text, float xp, float yp, int size){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(63);
	tex->SetTextSize(size);
	tex->SetTextColor(kBlack);
	tex->SetLineWidth(1);
	//tex->SetTextFont(42);
	tex->SetNDC();
	tex->Draw();
}


void putCMSPrel(double x=0.1, double y=0.9, double size=0.04){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Preliminary");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}


void putCMSSim(double x=0.1, double y=0.9, double size=0.04){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Simulation");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}


TLegend *myLegend(double x1,double y1,double x2, double y2)
{
	TLegend *leg = new TLegend(x1,y1,x2,y2);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	return leg; 

}

// Remove bins with error > central value
void cleanup(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val1 = h->GetBinContent(i);
		double valErr1 = h->GetBinError(i);
		if (valErr1>=val1) {
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}   

}

void makeHistTitle(TH1 *h,char *title, char *xTitle, char *yTitle, int color = -1, bool centerTitle = 1)
{
	h->SetTitle(title);
	h->SetXTitle(xTitle);
	h->SetYTitle(yTitle);
	
	if (centerTitle) {
		h->GetXaxis()->CenterTitle();
		h->GetYaxis()->CenterTitle();
		
	}
	
	if (color!=-1) {
		h->SetLineColor(color);
		h->SetMarkerColor(color);
	}
	
	h->GetYaxis()->SetNdivisions(610); 
	h->GetXaxis()->SetNdivisions(505);

	
	h->GetYaxis()->SetLabelFont(43);
	h->GetYaxis()->SetTitleFont(43);
	h->GetYaxis()->SetLabelSize(20);
	h->GetYaxis()->SetTitleSize(22);
	h->GetYaxis()->SetTitleOffset(2.6);
	
	h->GetXaxis()->SetLabelFont(43);
	h->GetXaxis()->SetTitleFont(43);
	h->GetXaxis()->SetLabelSize(20);
	h->GetXaxis()->SetTitleSize(22);
	h->GetXaxis()->SetTitleOffset(3.1);
	
	h->GetXaxis()->SetNoExponent();
	h->GetXaxis()->SetMoreLogLabels();
	
	h->GetXaxis()->SetTitleOffset(2.4);
	
}

// draw envelope using a systematic uncertainty histogram
TH1F* drawEnvelope(TH1F *h,char *opt,int color = kGray,int fillStyle = 0, int fillColor = 0,double shift = 0)
{
	TH1F *hClone = (TH1F*) h->Clone(Form("%s_mirror",h->GetTitle()));
	TH1F *hMirror = (TH1F*) h->Clone(Form("%s_mirror2",h->GetTitle()));
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val = h->GetBinContent(i);
		hMirror->SetBinContent(i,1-fabs(val-1)+shift);
		hClone->SetBinContent(i,val+shift);
	}
	
	//   hMirror->SetLineStyle(2);
	//   h->SetLineStyle(2);
	hMirror->SetLineColor(color);
	hMirror->SetFillColor(fillColor);
	hMirror->SetFillStyle(fillStyle);
	hClone->SetLineColor(color);
	hClone->SetFillColor(fillColor);
	hClone->SetFillStyle(fillStyle);
	hClone->Draw(opt);
	hMirror->Draw(opt);
	return hMirror;
}

void makeMultiPanelCanvasWithGap(TCanvas*& canv,
								 const Int_t columns,
								 const Int_t rows,
								 const Float_t leftOffset,
								 const Float_t bottomOffset,
								 const Float_t leftMargin,
								 const Float_t bottomMargin,
								 const Float_t edge, const Float_t asyoffset) {
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
	
	//PadHeight = 0.5*PadHeight;
	
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
		if(i==rows-2){
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
		}else{
			Ylow[i] = Yup[rows-1] + ct*PadHeight;
			Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
			//Yup[i] = 0.2*Yup[i];
		}
		ct++;
	}
	
	TString padName;
	for(Int_t i=0;i<columns;i++) {
		for(Int_t j=0;j<rows;j++) {
			canv->cd();
			padName = Form("p_%d_%d",i,j);
			//pad[i][j] = new TPad(padName.Data(),padName.Data(),
			//Xlow[i],Ylow[j],Xup[i],Yup[j]);
			// this is hacked version to create aysmmetric pads around low 
			if(j==0){
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j]-asyoffset,Xup[i],Yup[j]);
			}else{
				pad[i][j] = new TPad(padName.Data(),padName.Data(),
									 Xlow[i],Ylow[j],Xup[i],Yup[j]-asyoffset);
			}
			
			
			if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
			else pad[i][j]->SetLeftMargin(0);
			
			if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
			else pad[i][j]->SetRightMargin(0);
			
			if(j==0) pad[i][j]->SetTopMargin(edge);
			//else pad[i][j]->SetTopMargin(0.01);
			else pad[i][j]->SetTopMargin(0.02);
			
			//if(j==0) pad[i][j]->SetTopMargin(edge*3.5);
			//else pad[i][j]->SetTopMargin(0.0);
			
			if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
			else pad[i][j]->SetBottomMargin(0.15);
			
			pad[i][j]->Draw();
			pad[i][j]->cd();
			pad[i][j]->SetNumber(columns*j+i+1);
			
		}
	}
}

void prepareNcollUnc(int nbins, float maxpt=300.){
	
	int fillsty = 1001;
	
	const int n = nbins;
	
	double xvalue[n];
	double yvalue[n];
	double xerror[n];
	double yerror1[n], yerror2[n], yerror3[n], yerror4[n], yerror5[n], yerror6[n];
	


  for(int i=0;i<nbins;i++){

    xvalue[i] = 300.1 + 1.2*(double)i, yvalue[i]=1.0, xerror[i]=0.0;  
    
      // TAA
    yerror1[i]=0.0409, yerror2[i]=0.0459, yerror3[i]=0.0578, yerror4[i]=0.0944, yerror5[i]=0.143, yerror6[i]=0.176;
 
    
    // add 6% error 
    yerror1[i]=TMath::Sqrt(yerror1[i]*yerror1[i]+0.06*0.06);
    yerror2[i]=TMath::Sqrt(yerror2[i]*yerror2[i]+0.06*0.06);
    yerror3[i]=TMath::Sqrt(yerror3[i]*yerror3[i]+0.06*0.06);
    yerror4[i]=TMath::Sqrt(yerror4[i]*yerror4[i]+0.06*0.06);
    yerror5[i]=TMath::Sqrt(yerror5[i]*yerror5[i]+0.06*0.06);
    yerror6[i]=TMath::Sqrt(yerror6[i]*yerror6[i]+0.06*0.06);
    cout<<"TAA + Lumi uncertainty = "<<yerror1[i]<<endl;   
  
  }
  
  
 // int ci = 29;
  int ci = 15;

  tTAAerr[0] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror1);
  tTAAerr[0]->SetFillColor(ci);
  tTAAerr[0]->SetLineColor(ci);
  tTAAerr[0]->SetFillStyle(fillsty);

  tTAAerr[1] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror2);
  tTAAerr[1]->SetFillColor(ci);
  tTAAerr[1]->SetFillStyle(fillsty);

  tTAAerr[2] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror3);
  tTAAerr[2] ->SetFillColor(ci);
  tTAAerr[2] ->SetFillStyle(fillsty);

  tTAAerr[3]  = new TGraphErrors(n,xvalue,yvalue,xerror,yerror4);
  tTAAerr[3]->SetFillColor(ci);
  tTAAerr[3]->SetFillStyle(fillsty);

  tTAAerr[4] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror5);
  tTAAerr[4]->SetFillColor(ci);
  tTAAerr[4]->SetFillStyle(fillsty);

  tTAAerr[5] = new TGraphErrors(n,xvalue,yvalue,xerror,yerror6);
  tTAAerr[5]->SetFillColor(ci);
  tTAAerr[5]->SetFillStyle(fillsty);
  


}
/*
void DrawNpartTAABand(){
	double xvalueNpart[6];
	double yerrorNpart[6];
	xvalueNpart[0] = 381.29; xvalueNpart[1] = 329.41; xvalueNpart[2] = 224.28;
  	xvalueNpart[3] = 108.12; xvalueNpart[4] = 42.04;  xvalueNpart[5] = 11.43;
  	yerrorNpart[0]=0.0409, yerrorNpart[1]=0.0459, yerrorNpart[2]=0.0578, yerrorNpart[3]=0.0944, yerrorNpart[4]=0.143, yerrorNpart[5]=0.176;
  	
  	int ci = 30;
  	
		for (int i=0;i<6;i++) {
		
			TBox *b = new TBox(xvalueNpart[i]-5,1.-yerrorNpart[i]/2,xvalueNpart[i]+5,1.+yerrorNpart[i]/2);
			b->SetFillColor(ci);
			b->SetFillStyle(3001);
			b->SetLineColor(ci);
			b->Draw();
		}
 
 
 }
 
*/
 

void dumpDatatoTxt(const char *centbin,TH1F *h, TH1F *hsys, TH1F *htotStat, const char *txtfile){
	ofstream outf(txtfile,ios::out);
	for(int ix=1;ix<=h->GetNbinsX();ix++){
		double pt = h->GetBinCenter(ix);
		double val = h->GetBinContent(ix);
		double Uncorerr = h->GetBinError(ix);
		double syserr = hsys->GetBinContent(ix)-1;
		double totStaterr = htotStat->GetBinError(ix);

		outf<<setprecision(0) << fixed <<pt<<"\t" << setprecision(3) << fixed <<val<<"\t" << setprecision(5) << fixed << totStaterr<<"\t"<< setprecision(4) << fixed << syserr*val << endl;
	}

	outf.close();
}


TGraphErrors *ShiftGraph(TGraphErrors* pGraph, Double_t pNumber){
    // shifts a graph by the absolute value in the argument
    TGraphErrors *pGraphtmp;
    for (Int_t i=0;i<pGraph->GetN();i++){
	Double_t x,y;
	Double_t yerr;
	pGraph->GetPoint(i,x,y);
	yerr = pGraph->GetErrorY(i);
	x = x+pNumber;
	pGraphtmp->SetPoint(i,x,y);
	pGraphtmp->SetPointError(i,yerr);
	}
	return pGraphtmp;
}


TGraphErrors* HistToTgraphShift(TH1F* hist,Double_t pNumber){

  TGraphErrors *pGraphtmp;
  int nbins = hist->GetNbinsX();

  const int nlines = nbins;

  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    pterr[i] = 0;
  }

  pGraphtmp = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  
  for (Int_t i=0;i<pGraphtmp->GetN();i++){
	Double_t x,y;
	Double_t yerr;
	pGraphtmp->GetPoint(i,x,y);
	yerr = pGraphtmp->GetErrorY(i);
	x = x+pNumber;
	pGraphtmp->SetPoint(i,x,y);
	pGraphtmp->SetPointError(i,yerr);
	}
  return pGraphtmp;
}


void DrawPanelLabel(int i){
 	TLatex *tex; 
	
 	 if(i==0) tex = new TLatex(0.07,0.9,"(f)");
	 if(i==1) tex = new TLatex(0.07,0.9,"(e)");
	 if(i==2) tex = new TLatex(0.19,0.9,"(d)");
	 if(i==3) tex = new TLatex(0.07,0.9,"(c)");
	 if(i==4) tex = new TLatex(0.07,0.9,"(b)");
	 if(i==5) tex = new TLatex(0.19,0.9,"(a)");
	 tex->SetTextFont(63);
	 tex->SetTextSize(18);
	 tex->SetTextColor(kBlack);
	 tex->SetLineWidth(1);
	 tex->SetNDC();

	 
 	tex->Draw();
	
 }


using namespace std;

void RAA_plot(int radius = 3, char *algo = "Vs"){

  TStopwatch timer;
  timer.Start();

  TDatime date;

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int nbins_cent = 6;
  double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
  double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};

  TFile *fin = TFile::Open(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Output/PbPb_pp_unfo_ak%d_%s_cent%d_chMax_12003cut.root",radius,algo,date.GetDate()));

  TH1F *dPbPb_TrgComb[nbins_cent+1], *dPbPb_Comb[nbins_cent+1], *dPbPb_Trg80[nbins_cent+1], *dPbPb_Trg65[nbins_cent+1], *dPbPb_Trg55[nbins_cent+1], *dPbPb_1[nbins_cent+1], *dPbPb_2[nbins_cent+1], *dPbPb_3[nbins_cent+1], *dPbPb_80[nbins_cent+1], *dPbPb_65[nbins_cent+1], *dPbPb_55[nbins_cent+1];
  
  TH1F *mPbPb_Gen[nbins_cent+1], *mPbPb_Reco[nbins_cent+1];
  TH2F *mPbPb_Matrix[nbins_cent+1], *mPbPb_Response[nbins_cent+1], *mPbPb_ResponseNorm[nbins_cent+1];
  TH1F *mPbPb_mcclosure_data[nbins_cent+1];
  
  TH1F *dPP_1, *dPP_2, *dPP_3, *dPP_Comb;
  
  //TH1F *mPP_Gen, *mPP_Reco;
  TH2F *mPP_Matrix, *mPP_Response;
  TH2F *mPP_ResponseNorm;
  //TH1F *mPP_mcclosure_data;
  
  const int Iterations = 20; //for unfolding systematics. 
  const int BayesIter = 4;
  TH1F *uPbPb_Bayes[nbins_cent+1], *uPbPb_BinByBin[nbins_cent+1]; 
  TH1F *uPbPb_BayesianIter[Iterations][nbins_cent+1];
  TH1F *uPbPb_MC_Bayes[nbins_cent];
  TH1F *uPbPb_MC_BayesianIter[Iterations][nbins_cent+1];
  TH1F *uPbPb_MC_BinByBin[nbins_cent+1];

  TH1F *uPP_Bayes, *uPP_BinByBin;
  TH1F *uPP_BayesianIter[Iterations];
  TH1F *uPP_MC_Bayes, *uPP_MC_BinByBin;
  TH1F *uPP_MC_BayesianIter[Iterations];

  TH1F *RAA_measured[nbins_cent+1];
  TH1F *RAA_binbybin[nbins_cent+1];
  TH1F *RAA_bayesian[nbins_cent+1];

  for(int i = 0;i<=nbins_cent;i++){
    
    cout<<"cent = "<<i<<endl;
    dPbPb_TrgComb[i] = (TH1F*)fin->Get(Form("hPbPb_TrgComb_cent%d",i));
    mPbPb_Reco[i] = (TH1F*)fin->Get(Form("hPbPb_reco_cent%d",i));
    mPbPb_Gen[i] = (TH1F*)fin->Get(Form("hPbPb_gen_cent%d",i));
    mPbPb_ResponseNorm[i] = (TH2F*)fin->Get(Form("mPbPb_ResponseNorm_cent%d",i));
    mPbPb_Response[i] = (TH2F*)fin->Get(Form("mPbPb_Response_cent%d",i));
    uPbPb_Bayes[i] = (TH1F*)fin->Get(Form("uPbPb_Bayes_cent%d",i));
    uPbPb_MC_Bayes[i] = (TH1F*)fin->Get(Form("uPbPb_Bayes_MC_cent%d",i));
    uPbPb_BinByBin[i] = (TH1F*)fin->Get(Form("uPbPb_BinByBin_cent%d",i));
    uPbPb_MC_BinByBin[i] = (TH1F*)fin->Get(Form("uPbPb_MC_BinByBin_cent%d",i));
    RAA[i] = (TH1F*)fin->Get(Form("RAA_cent%d",i));
    mPbPb_mcclosure_data[i] = (TH1F*)fin->Get(Form("mPbPb_mcclosure_data_cent%d",i));

    for(int j = 0;j<Iterations;j++){
    	uPbPb_BayesianIter[j][i] = (TH1F*)fin->Get(Form("uPbPb_BayesianIter%d_cent%d",j,i));
    	uPbPb_MC_BayesianIter[j][i] = (TH1F*)fin->Get(Form("uPbPb_MC_BayesianIter%d_cent%d",j,i));
    }

  }

  dPP_Comb = (TH1F*)fin->Get("hppComb");
  mPP_ResponseNorm = (TH1F*)fin->Get("mPP_ResponseNorm");	
  mPP_mcclosure_data = (TH1F*)fin->Get("mPP_mcclosure_data");

  uPP_Bayes = (TH1F*)fin->Get("uPP_Bayes");
  uPP_BinByBin = (TH1F*)fin->Get("uPP_BinByBin");
  uPP_MC_Bayes = (TH1F*)fin->Get("uPP_MC_Bayes");
  uPP_MC_BinByBin = (TH1F*)fin->Get("uPP_MC_BinByBin");
  mPP_Gen = (TH1F*)fin->Get("hpp_gen");
  mPP_Reco = (TH1F*)fin->Get("hpp_reco");

  for(int i = 0;i<Iterations;i++){

  	uPP_BayesianIter[i] = (TH1F*)fin->Get(Form("uPP_BayesianIter%d",i));
  	uPP_MC_BayesianIter[i] = (TH1F*)fin->Get(Form("uPP_MC_BayesianIter%d",i));

  }

  //Ok now that we have loaded all the histograms we need - lets start making the plots 

  // line at 1
  TLine *line = new TLine(50,1,300,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 1 - PbPb iteration systematics. 
  // this will be a 3 by 2 panel plot showing bayesian for PbPb. per centrality bin. 
  // divide different unfolding iterations with iteration 4 - the nominal one. 
  TCanvas *cIterSysPbPb = new TCanvas("cIterSysPbPb","PbPb Iteration systematics",1200,800);
  makeMultiPanelCanvasWithGap(cIterSysPbPb,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
 
  for(int i = 0;i<nbins_cent;i++){
  	cIterSysPbPb->cd(nbins_cent-i);

  	TLegend *PbPb_itersys = myLegend(0.53,0.65,0.85,0.9);
  	line->Draw();

  	drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

  	for(int j = 2;j<7;j++){
  		
  		uPbPb_BayesianIter[j][i]->Divide(uPbPb_Bayes);
  		uPbPb_BayesianIter[j][i]->SetMarkerStyle(2);
  		uPbPb_BayesianIter[j][i]->SetMarkerColor(j);
  		if(j==2){
  			makeHistTitle(uPbPb_BayesianIter[j][i]," ","Jet p_{T} (GeV/c)","Ratio (Unfolded/Nominal)");
  			uPbPb_BayesianIter[j][i]->Draw();
  		}else {
  			uPbPb_BayesianIter[j][i]->Draw("same");
  		}
  		if(i==0) PbPb_itersys->AddEntry(uPbPb_BayesianIter[j][i],Form("Iteration %d",j),"pl");

  	}


  }
  PbPb_itersys->Draw();

  cIterSysPbPb->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Particle Flow Jets R=0.%d",algo, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cIterSysPbPb->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PbPb_unfoldng_iteration_systematics_ak%d_%s_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 2 - pp iteration systematics 
  // this is just one panel plot showing bayesian for pp. 
  TCanvas *cIterSysPP = new TCanvas("cIterSysPP","PP Iteration systematics",600,400);

  TLegend *PP_itersys = myLegend(0.53,0.65,0.85,0.9);
  line->Draw();

  for(int i = 2;i<7;i++){
  	uPP_BayesianIter[i]->Divide(uPP_Bayes);
  	uPP_BayesianIter[i]->SetMarkerStyle(2);
  	uPP_BayesianIter[i]->SetMarkerColor(i);
  	if(i==2){
  		makeHistTitle(uPP_BayesianIter[i]," ","Jet p_{T} (GeV/c)","Ratio (unfolded/Nominal)");
  		uPP_BayesianIter[i]->Draw();
  	}else uPP_BayesianIter[i]->Draw("same");
  	PP_itersys->AddEntry(uPP_BayesianIter[i],Form("Iteration %d",i),"pl");
  }

  PP_itersys->Draw();

  cIterSysPbPb->cd(1);
  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d", radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cIterSysPP->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PP_unfoldng_iteration_systematics_ak%d_%s_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 3 - RAA 
  // again this will be a 6 panel plot. showing measured, unfolded Bayesian, and unfolded Bin By Bin methods. 
  TCanvas *cRAA = new TCanvas("cRAA","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TLegend *tRAA = myLegend(0.53,0.65,0.85,0.9);

  for(int i = 0;i<nbins_cent;i++){

  	cRAA->cd(nbins_cent-i);
  	drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

  	RAA_measured[i]->SetMarkerColor(kBlack);
  	RAA_measured[i]->SetMarkerStyle(24);
  	RAA_measured[i]->Draw();

  	RAA_bayesian[i]->SetMarkerColor(kRed);
  	RAA_bayesian[i]->SetMarkerStyle(33);
  	RAA_bayesian[i]->Draw("same");

  	RAA_binbybin[i]->SetMarkerStyle(29);
  	RAA_binbybin[i]->SetMarkerColor(kBlue);
  	RAA_binbybin[i]->DraW("same");

  	line->DraW();

  }

  tRAA->AddEntry(RAA_measured[0],"No Unfolding","pl");
  tRAA->AddEntry(RAA_bayesian[0],"Bayesian","pl");
  tRAA->AddEntry(RAA_binbybin[0],"BinbyBin","pl");

  cRAA->cd(0);
  tRAA->Draw();
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Particle Flow Jets R=0.%d",algo, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cRAA->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/RAA_ak%d_%s_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 4 - PbPb MC closure test 
  // this will also be a 6 panel plot showing bayesian iteration 4 and binbybin divided by measured which is actually the MC
  TCanvas *cPbPbMCclosure = new TCanvas("cPbPbMCclosure","PbPb MC closure test",1200,800);
  makeMultiPanelCanvasWithGap(cPbPbMCclosure,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TH1F *hMCClosurePbPb_Meas[nbins_cent+1], *hMCClosurePbPb_Bayesian, *hMCClosurePbPb_BinByBin[nbins_cent+1];
  for(int i = 0;<nbins_cent;i++){

  	cPbPbMCclosure->cd(nbins_cent-i);
  	drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);
  	hMCClosurePbPb_Meas[i] = (TH1F*)mPbPb_Reco[i]->Clone(Form("hMCClosurePbPb_Meas_cent%d",i));
  	hMCClosurePbPb_Bayesian[i] = (TH1F*)uPbPb_MC_Bayes[i]->Clone(Form("hMCClosurePbPb_Bayesian_cent%d",i));
  	hMCClosurePbPb_BinByBin[i] = (TH1F*)uPbPb_MC_BinByBin[i]->Clone(Form("hMCClosurePbPb_BinByBin_cent%d",i));

  	hMCClosurePbPb_Meas[i]->Divide(mPbPb_mcclosure_data[i]);
  	hMCClosurePbPb_Bayesian[i]->Divide(mPbPb_mcclosure_data[i]);
  	hMCClosurePbPb_BinByBin[i]->Divide(mPbPb_mcclosure_data[i]);

  	makeHistTitle(hMCClosurePbPb_Meas[i]," ","Jet p_{T} (GeV/c)","Reco/Truth");
  	hMCClosurePbPb_Meas[i]->SetMarkerStyle(24);
  	hMCClosurePbPb_Meas[i]->SetMarkerColor(kBlack);
  	hMCClosurePbPb_Meas[i]->Draw();

  	hMCClosurePbPb_Bayesian[i]->SetMarkerStyle(33);
  	hMCClosurePbPb_Bayesian[i]->SetMarkerColor(kRed);
  	hMCClosurePbPb_Bayesian[i]->Draw("same");

  	hMCClosurePbPb_BinByBin[i]->SetMarkerColor(kBlue);
  	hMCClosurePbPb_BinByBin[i]->SetMarkerStyle(29);
  	hMCClosurePbPb_BinbyBin[i]->Draw("same");

  	line->Draw();

  }

  cPbPbMCclosure->cd(1);
  TLegend *pbpbmcclosure = myLegend(0.53,0.65,0.85,0.9);
  pbpbmcclosure->AddEntry(hMCClosurePbPb_Meas[0],"PbPb no unfolding","pl");
  pbpbmcclosure->AddEntry(hMCClosurePbPb_Bayesian[0],"PbPb Bayesian 4 Iter","pl");
  pbpbmcclosure->AddEntry(hMCClosurePbPb_BinbyBin[0],"PbPb BinbyBin","pl");
  pbpbmcclosure->Draw();
  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Particle Flow Jets R=0.%d",algo, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPbPbMCclosure>SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PbPb_unfoldng_mc_closure_test_ak%d_%s_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 5 - PP MC closure test
  TCanvas * cPPMCclosure = new TCanvas("cPPMCclosure","PP MC closure test",600,400);
  TH1F *hMCClosurePP_Meas = (TH1F*)mPP_Reco->Clone("hMCClosurePP_Meas");
  TH1F *hMCClosurePP_Bayesian = (TH1F*)uPP_MC_Bayes->Clone("hMCClosurePP_Bayesian");
  TH1F *hMCClosurePP_BinbyBin = (TH1F*)uPP_MC_BinByBin->Clone("hMCClosurePP_BinbyBin");

  hMCClosurePP_Bayesian->Divide(mPP_mcclosure_data);
  hMCClosurePP_Meas->Divide(mPP_mcclosure_data);
  hMCClosurePP_BinbyBin->Divide(mPP_mcclosure_data);

  makeHistTitle(hMCClosurePP_Meas," ","Jet p_{T} (GeV/c)","Reco/Truth");
  hMCClosurePP_Meas->SetAxisRange(50,300,"X");
  hMCClosurePP_Meas->SetAxisRange(0,2,"Y");

  hMCClosurePP_Meas->SetMarkerStyle(24);
  hMCClosurePP_Meas->SetMarkerColor(kBlack);
  hMCClosurePP_Meas->Draw();

  hMCClosurePP_Bayesian->SetMarkerStyle(33);
  hMCClosurePP_Bayesian->SetMarkerColor(kRed);
  hMCClosurePP_Bayesian->Draw("same");

  hMCClosurePP_BinbyBin->SetMarkerStyle(29);
  hMCClosurePP_BinbyBin->SetMarkerColor(kBlue);
  hMCClosurePP_BinbyBin->Draw("same");

  line->Draw();

  TLegend *ppmcclosure = myLegend(0.53,0.65,0.85,0.9);
  ppmcclosure->AddEntry(hMCClosurePP_Meas,"pp no unfolding","pl");
  ppmcclosure->AddEntry(hMCClosurePP_Bayesian,"pp Bayesian","pl");
  ppmcclosure->AddEntry(hMCClosurePP_BinbyBin,"pp BinbyBin","pl");
  ppmcclosure->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,22);
  
  cPPMCclosure>SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PP_unfoldng_mc_closure_test_ak%d_%d.pdf",radius,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 6 - data vs MC for PbPb 
  // i want to make it like ratio plots for each centrality class. and plot the unfolded ratio as well. 
  // so i want a line at 1, and 2 overlayed plots showing measured/Gen ratio and unfo/gen ratio. 
  TCanvas *cPbPb_data_vs_mc = new TCanvas("cPbPb_data_vs_mc","PbPb data vs mc",1200,800);
  makeMultiPanelCanvasWithGap(cPbPb_data_vs_mc,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  TH1F *PbPb_meas_vs_mc[nbins_cent+1];
  TH1F *PbPb_unfo_vs_mc[nbins_cent+1];

  for(int i = 0;i<nbins_cent;i++){

  	cPbPb_data_vs_mc->cd(nbins_cent-i);

  	PbPb_meas_vs_mc[i] = (TH1F*)dPbPb_TrgComb[i]->Clone("PbPb_meas_vs_mc");
  	PbPb_unfo_vs_mc[i] = (TH1F*)uPbPb_Bayes[i]->Clone("PbPb_unfo_vs_mc");

  	PbPb_meas_vs_mc[i]->Divide(mPbPb_Reco[i]);
  	PbPb_unfo_vs_mc[i]->Divide(mPbPb_Reco[i]);

  	makeHistTitle(PbPb_meas_vs_mc[i]," ","Jet p_{T} (GeV/c)", "data/mc");
  	PbPb_meas_vs_mc[i]->SetMarkerStyle(25);
  	PbPb_meas_vs_mc[i]->SetMarkerColor(kBlack);
  	PbPb_meas_vs_mc[i]->Draw();

  	PbPb_unfo_vs_mc[i]->SetMarkerStyle(27);
  	PbPb_unfo_vs_mc[i]->SetMarkerColor(kBlue);
  	PbPb_unfo_vs_mc[i]->Draw("same");

  	line->Draw();

  	drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

  }

  cPbPb_data_vs_mc->cd(1);
  TLegend *PbPb_datavsmc = myLegend(0.53,0.65,0.85,0.9);
  PbPb_datavsmc->AddEntry(PbPb_meas_vs_mc,"no unfolding","pl");
  PbPb_datavsmc->AddEntry(PbPb_unfo_vs_mc,"Bayesian 4 Iter","pl");
  PbPb_datavsmc->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} %s Particle Flow Jets R=0.%d",algo, radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPbPb_data_vs_mc->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PbPb_data_vs_mc_ak%d_%s_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot7 - data vs MC for pp
  // made similar to above. 
  TCanvas *cPP_data_vs_mc = new TCanvas("cPP_data_vs_mc","PP data vs MC",600,400);

  TH1F *PP_meas_vs_mc = (TH1F*)dPP_Comb->Clone("PP_meas_vs_mc");
  TH1F *PP_unfo_vs_mc = (TH1F*)uPP_Bayes->Clone("PP_unfo_vs_mc");

  PP_meas_vs_mc->Divide(mPP_Reco);
  PP_unfo_vs_mc->Divide(mPP_Reco);

  makeHistTitle(PP_meas_vs_mc," ","Jet p_{T} (GeV/c)","data/mc");
  PP_meas_vs_mc->SetMarkerStyle(25);
  PP_meas_vs_mc->SetMarkerColor(kBlack);
  PP_meas_vs_mc->Draw();
  PP_unfo_vs_mc->SetMarkerStyle(27);
  PP_unfo_vs_mc->SetMarkerColor(kBlue);
  PP_unfo_vs_mc->Draw("same");

  line->Draw();

  TLegend *PP_datavsmc = myLegend(0.53,0.65,0.85,0.9);
  PP_datavsmc->AddEntry(PP_meas_vs_mc,"no unfolding","pl");
  PP_datavsmc->AddEntry(PP_unfo_vs_mc,"Bayesian 4 Iter","pl");
  PP_datavsmc->SetTextSize(0.02);
  PP_datavsmc->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);

  cPP_data_vs_mc->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PP_data_vs_mc_ak%d_%d.pdf",radius,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot8 - normalized Response Matrix for PbPb. 
  TCanvas *cPbPb_NormResMat = new TCanvas("cPbPb_NormResMat","Normalized Response Matrix for PbPb",1200,800);
  makeMultiPanelCanvasWithGap(cPbPb_NormResMat,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0;i<nbins_cent;i++){
  	cPbPb_NormResMat->cd(nbins_cent-i);

  	makeHistTitle(mPbPb_ResponseNorm[i]," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");
  	drawText(Form("%2.0f-%2.0f%%",5*boundaries_cent[i],5*boundaries_cent[i+1]),0.8,0.9,20);

  	mPbPb_ResponseNorm[i]->SetAxisRange(1e-10,1,"Z");
  	mPbPb_ResponseNorm[i]->Draw("colz");

  }
  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow %s Jets R=0.%d",algo,radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPbPb_NormResMat->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PbPb_normalized_response_matrix_ak%d_%s_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot9 - normalized response matrix for pp
  TCanvas *cPP_NormResMat = new TCanvas("cPP_NormResMat","Normalized Response Matrix fr PP",600,400);

  makeHistTitle(mPP_ResponseNorm," ","Gen p_{T} (GeV/c)","Reco p_{T} (GeV/c)");

  mPP_ResponseNorm->SetAxisRange(1e-10,1,"Z");
  mPP_ResponseNorm->Draw("colz");

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPP_NormResMat->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PP_normalized_response_matrix_ak%d_%d.pdf",radius,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot10 - Cross section for PbPb - in proper units per centrality bin.  
  TCanvas *cPbPb_sigma = new TCanvas("cPbPb_sigma","PbPb inclusive jet invariant cross section",1200,800);
  makeMultiPanelCanvasWithGap(cPbPb_sigma,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

  for(int i = 0;i<nbins_cent;i++){

  	cPbPb_sigma->cd(nbins_cent-i);
  	cPbPb_sigma->cd(nbins_cent-i)->SetLogy();

  	makeHistTitle(dPbPb_TrgComb[i]," ","Jet p_{T} (GeV/c)","#frac{d^2 #sigma}{d#eta dp_{T}} micro barns");
  	dPbPb_TrgComb[i]->SetMarkerStyle(24);
  	dPbPb_TrgComb[i]->SetMarkerColor(kBlack);
  	dPbPb_TrgComb[i]->Draw();

  	uPbPb_Bayes[i]->SetMarkerStyle(33);
  	uPbPb_Bayes[i]->SetMarkerColor(kRed);
  	uPbPb_Bayes[i]->DraW("same");

  	uPbPb_BinByBin[i]->SetMarkerStyle(29);
  	uPbPb_BinByBin[i]->SetMarkerColor(kBlue);
  	uPbPb_BinByBin[i]->Draw("same");

  }

  cPbPb_sigma->cd(1);
  TLegend *PbPb_sigma = myLegend(0.53,0.65,0.85,0.9);
  PbPb_sigma->AddEntry(dPbPb_TrgComb[0],"No Unfolding","pl");
  PbPb_sigma->AddEntry(uPbPb_Bayes[0],"Bayesian 4 Iter","pl");
  PbPb_sigma->AddEntry(uPbPb_BinByBin[0],"BinByBin","pl");
  PbPb_sigma->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow %s Jets R=0.%d",algo,radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPbPb_sigma->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PbPb_invariant_cross_section_ak%d_%s_%d.pdf",radius,algo,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  //plot 11 - cross section for PP 
  TCanvas *cPP_sigma = new TCanvas("cPP_sigma","PP inclusive jet invariant cross section",600,400);
  cPP_sigma->SetLogy();
  
  makeHistTitle(dPP_Comb,"","Jet p_{T} (GeV/c)","#frac{d^2 #sigma}{d#eta dp_{T}} nano barns");
  dPP_Comb->SetMarkerStyle(24);
  dPP_Comb->SetMarkerColor(kBlack);
  dPP_Comb->Draw();

  uPP_Bayes->SetMarkerColor(kRed);
  uPP_Bayes->SetMarkerStyle(33);
  uPP_Bayes->Draw("same");

  uPP_BinByBin->SetMarkerStyle(29);
  uPP_BinByBin->SetMarkerColor(kBlue);
  uPP_BinByBin->Draw("same");

  TLegend *PP_sigma = myLegend(0.53,0.65,0.85,0.9);
  PP_sigma->AddEntry(dPP_Comb,"No Unfolding","pl");
  PP_sigma->AddEntry(uPP_Bayes,"Bayesian 4 Iter","pl");
  PP_sigma->AddEntry(uPP_BinByBin,"BinByBin","pl");
  PP_sigma->Draw();

  putCMSPrel();
  drawText(Form("Anti-k_{T} Particle Flow Jets R=0.%d",radius),0.2,0.23,20);
  drawText("|#eta|<2, |vz|<15",0.6,0.31,20);
  cPbPb_sigma->SaveAs(Form("/afs/cern.ch/work/r/rkunnawa/WORK/RAA/CMSSW_5_3_18/src/Plots/PP_invariant_cross_section_ak%d_%d.pdf",radius,date.GetDate()),"RECREATE");
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;


}





































