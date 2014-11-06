// Standard library
#include <math.h>
#include <iostream>
#include <fstream>

// ROOT Library
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TCut.h>
#include <TPad.h>

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

void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.425,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.425);
  c->GetPad(2)->SetBottomMargin(0.3);
  c->GetPad(2)->SetGridy(1);
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


void putCMSPrel(double x=0.17, double y=0.96, double size=0.04){
  TLatex *tex=0;
  tex = new TLatex(x,y,"CMS Preliminary");
  tex->SetTextSize(size);
  tex->SetLineWidth(2);
  tex->SetNDC();
  tex->Draw();
}


void putCMSSim(double x=0.17, double y=0.96, double size=0.04){
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

void makeMultiPanelCanvasWithoutGap(TCanvas*& canv,
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


/*
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

*/

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
