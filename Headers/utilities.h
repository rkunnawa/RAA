
#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <TRandom.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TCanvas.h"
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TBox.h>
#include "TF1.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

using namespace std;


// TAA 

TGraphErrors *tTAAerr[6]={0};
TGraphErrors *tTAAerrNpart=0;

const int nbins_cent = 6;
double boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};
double ncoll[nbins_cent+1] = {1660,1310,745,251,62.8,10.8,362.24};
double npart[nbins_cent+1] = {389.84, 307.65, 223.95, 107.5, 41.65, 11.55, 112.9};


const int nbins_pt = 30;
const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12, 
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 300, 
  330, 362, 395 
};

class SysData
{
 public:
  SysData() {
    for (int i=0;i<=nbins_cent;i++) {
      hSys[i]     = new TH1F(Form("hSys_cent%d",i), Form("Totalsys_cent%d",i),nbins_pt, boundaries_pt);
      hSysGeneral[i]= new TH1F(Form("hSysGeneral_cent%d",i), Form("TotalsysGeneral_cent%d",i),nbins_pt, boundaries_pt);
      hSysJEC[i]  = new TH1F(Form("hSysJEC_cent%d",i), Form("JECsys_cent%d",i),nbins_pt, boundaries_pt);
      hSysEff[i]  = new TH1F(Form("hSysEff_cent%d",i), Form("Effsys_cent%d",i),nbins_pt, boundaries_pt);
      hSysSmear[i]  = new TH1F(Form("hSysSmear_cent%d",i), Form("Smearsys_cent%d",i),nbins_pt, boundaries_pt);
      hSysIter[i] = new TH1F(Form("hSysIter_cent%d",i), Form("Itersys_cent%d",i),nbins_pt, boundaries_pt);
      hSysJetID[i] = new TH1F(Form("hSysJetID_cent%d",i), Form("JetID_sys_cent%d",i), nbins_pt, boundaries_pt);
      hSys[i]->SetLineColor(kGray);
      hSysJEC[i]->SetLineColor(4);
      hSysSmear[i]->SetLineColor(kGreen+1);
      hSysIter[i]->SetLineColor(2);
      hSysJetID[i]->SetLineColor(kGreen+2);
    }  
  }
  TH1F *hSys[nbins_cent+1];
  TH1F *hSysGeneral[nbins_cent+1];
  TH1F *hSysJEC[nbins_cent+1];
  TH1F *hSysEff[nbins_cent+1];
  TH1F *hSysSmear[nbins_cent+1];
  TH1F *hSysIter[nbins_cent+1];
  TH1F *hSysNoise[nbins_cent+1];
  TH1F *hSysJetID[nbins_cent+1];
	
  void calcTotalSys(int i) {
    TF1 *fNoise = new TF1("f","1+0.3*0.16*abs(1-([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x))");
    fNoise->SetParameters(0.9521,0.001105,-9.397e-6,3.32e-8,-5.618e-11);
    hSysNoise[i] = functionHist(fNoise,hSys[i],Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetName(Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetLineColor(6);
    for (int j=1;j<=hSys[i]->GetNbinsX();j++) {
      double effSys = 0.01;
      double jetidSys = 0.02;
      hSysEff[i]->SetBinContent(j,1+effSys);
      hSysSmear[i]->SetBinContent(j,1.02);
      hSysJetID[i]->SetBinContent(j, 1+jetidSys);
      double JECSys = hSysJEC[i]->GetBinContent(j)-1;
      double SmearSys = hSysSmear[i]->GetBinContent(j)-1;
      double IterSys = hSysIter[i]->GetBinContent(j)-1; 
      double NoiseSys = hSysNoise[i]->GetBinContent(j)-1;
      cout <<effSys<<" "<<JECSys<<" "<<IterSys<<endl;
      double totalSys = sqrt( effSys * effSys +
			      JECSys * JECSys +
			      SmearSys * SmearSys +
			      NoiseSys * NoiseSys +
			      IterSys* IterSys +
			      jetidSys * jetidSys
			      );
      hSys[i]->SetBinContent(j,totalSys+1);	
      hSys[i]->SetLineWidth(2);				
    }
  }
	
  void calcTotalSysNoUnfolding(int i) {
    TF1 *fNoise = new TF1("f","1+0.3*0.16*abs(1-([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x))");
    fNoise->SetParameters(0.9521,0.001105,-9.397e-6,3.32e-8,-5.618e-11);
    hSysNoise[i] = functionHist(fNoise,hSys[i],Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetName(Form("hSysNoise_cent%d",i));
    hSysNoise[i]->SetLineColor(6);
    for (int j=1;j<=hSysGeneral[i]->GetNbinsX();j++) {
      double effSys = 0.01;
      double jetidSys = 0.02;
      hSysEff[i]->SetBinContent(j,1+effSys);
      hSysSmear[i]->SetBinContent(j,1.02);
      hSysJetID[i]->SetBinContent(j, 1+jetidSys);
      double JECSys = hSysJEC[i]->GetBinContent(j)-1;
      double SmearSys = hSysSmear[i]->GetBinContent(j)-1;
      double NoiseSys = hSysNoise[i]->GetBinContent(j)-1; 
      double totalSys = sqrt( effSys * effSys +
			      JECSys * JECSys +
			      SmearSys * SmearSys +
			      NoiseSys * NoiseSys +
			      jetidSys * jetidSys
			      );
      hSysGeneral[i]->SetBinContent(j,totalSys+1);	
      hSysGeneral[i]->SetLineWidth(2);				
    }
  }
  
  void Draw(TH1F *h,int i, int color) {

    Int_t beginning = h->FindBin(60); Int_t end = h->FindBin(299);

    if(i==0){ beginning = h->FindBin(100); }
    if(i==5 || i==4 || i==3) { end = h->FindBin(240);}

    for (int j=beginning;j<=end;j++) {
      double val = h->GetBinContent(j);
      double err = hSys[i]->GetBinContent(j)-1;
      cout << "Sys Value Check" <<val<<" "<<err<<" "<<h->GetBinLowEdge(j)<<" "<<val*(1-err)<<" "<<h->GetBinLowEdge(j+1)<<" "<<val*(1+err)<<endl;
      TBox *b = new TBox(h->GetBinLowEdge(j),val*(1-err),h->GetBinLowEdge(j+1),val*(1+err));
      //b->SetFillColor(kGray);
      b->SetFillStyle(0);
      //b->SetLineColor(kGray);
			
      
      //***********For Gunther's Color Systematics Band Peference
      b->SetFillColor(color);
      b->SetLineColor(color);
      b->Draw();
    }
  }
	
  void DrawTGraph(TGraphErrors *h,int i) {
    double xv;
    double val;
    for(int j=0;j<h->GetN();j++){	
      h->GetPoint(j,xv,val);
      double err = hSysGeneral[i]->GetBinContent(j+1)-1;
      //cout <<"value" <<val<<" "<<err<<" "<<hSysGeneral[i]->GetBinLowEdge(j+1)<<" "<<val*(1-err)<<" "<<hSysGeneral[i]->GetBinLowEdge(j+2)<<" "<<val*(1+err)<<endl;
      TBox *b_ = new TBox(hSysGeneral[i]->GetBinLowEdge(j+1),val*(1-err),hSysGeneral[i]->GetBinLowEdge(j+2),val*(1+err));
      b_->SetFillColor(kGray);
      b_->SetFillStyle(1001);
      b_->SetLineColor(kGray);
      b_->Draw();
    }
  }
	
  void DrawUnfoErr(TH1F *h,int i) {
    for (int j=1;j<=hSysIter[i]->GetNbinsX();j++) {
      double val = h->GetBinContent(j);
      double err = hSysIter[i]->GetBinContent(j)-1;
      //cout <<"value" << val<<" "<<err<<" "<<h->GetBinLowEdge(j)<<" "<<val*(1-err)<<" "<<h->GetBinLowEdge(j+1)<<" "<<val*(1+err)<<endl;
      TBox *b = new TBox(h->GetBinLowEdge(j),val*(1-err),h->GetBinLowEdge(j+1),val*(1+err));
      b->SetFillColor(29);
      b->SetFillStyle(1001);
      b->SetLineColor(29);
      b->Draw();
    }
  }
	
  
  void DrawNpartSys(double yvalue,int i,double xvalue, int binno) {
    double yerrorNpart[6]= {0.0409, 0.0459,0.0578,0.0944, 0.143, 0.176 };
    double err = hSys[i]->GetBinContent(hSys[i]->FindBin(binno))-1;
    TBox *b = new TBox(xvalue-10.,yvalue*(1-err-yerrorNpart[i]),xvalue+10.,yvalue*(1+err+yerrorNpart[i]));
    //cout << "value " << yvalue<<" err   "<<err<<" xvalue  "<<xvalue<<" "<<yvalue*(1-err)<<" "<<yvalue*(1+err)<<endl;
    b->SetFillColor(kGray);
    b->SetFillStyle(1001);
    b->SetLineColor(kGray);
    
    //***********For Gunther's Color Systematics Band Peference
    //b->SetFillColor(5);
    //b->SetLineColor(5);
			
    b->Draw();
		
  }
  
  
  void DrawComponent(int i) {
    calcTotalSys(i);
    TH1D *h = new TH1D(Form("hSysTmp_cent%d",i),"",nbins_pt, boundaries_pt);
    makeHistTitle(h,"","Jet p_{T} (GeV/c)","Systematic uncertainty");
    h->SetAxisRange(-0.25,0.4,"Y");
    h->SetAxisRange(50,299,"X");
    h->Draw();
    TH1F* sys = drawEnvelope(hSys[i],"same",hSys[i]->GetLineColor(),1001,hSys[i]->GetLineColor(),-1);
    TH1F* sysIter = drawEnvelope(hSysIter[i],"same",hSysIter[i]->GetLineColor(),3004,hSysIter[i]->GetLineColor(),-1);
    TH1F* sysJEC = drawEnvelope(hSysJEC[i],"same",hSysJEC[i]->GetLineColor(),3005,hSysJEC[i]->GetLineColor(),-1);
    TH1F* sysJetID = drawEnvelope(hSysJetID[i],"same",hSysJetID[i]->GetLineColor(),3005,hSysJetID[i]->GetLineColor(),-1);
    TH1F* sysSmear =  drawEnvelope(hSysSmear[i],"same",hSysSmear[i]->GetLineColor(),3001,hSysSmear[i]->GetLineColor(),-1);
    TH1F* sysEff = drawEnvelope(hSysEff[i],"same",hSysEff[i]->GetLineColor(),3002,hSysEff[i]->GetLineColor(),-1);
    TH1F* sysNoise = drawEnvelope(hSysNoise[i],"same",hSysNoise[i]->GetLineColor(),3001,hSysNoise[i]->GetLineColor(),-1);
    TLine *l = new TLine(h->GetBinLowEdge(1),0,h->GetBinLowEdge(h->GetNbinsX()+1),0);
    l->Draw();
    TLine *l2 = new TLine(h->GetBinLowEdge(1),-0.25,h->GetBinLowEdge(1),0.4);
    l2->Draw();
    TLegend *leg = myLegend(0.52,0.6,0.95,0.93);
    leg->SetTextSize(0.043);
    leg->AddEntry(sys,"Total Systematics","f");
    leg->AddEntry(sysIter,"Unfolding","f");
    leg->AddEntry(sysJEC,"Jet Energy Scale","f");
    leg->AddEntry(sysEff,"Jet Trigger Efficiency","f");
    leg->AddEntry(sysJetID, "Jet ID efficiency","f");
    leg->AddEntry(sysSmear,"UE fluctuation","f");
    leg->AddEntry(sysNoise,"HCAL Noise","f");
    if (i==nbins_cent-1)leg->Draw();
  }
	
 
};

// Remove error 
void removeError(TH1F *h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      h->SetBinError(i,0);
    }   
	
}

// Remove Zero
void removeZero(TH1 *h)
{
  double min = 0;
  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i)>min&&h->GetBinContent(i)>0)
      min = h->GetBinContent(i);
  }

  for(int i = 1;i<h->GetNbinsX();i++){
    if(h->GetBinContent(i) == 0){
      h->SetBinContent(i,min/10.);
      h->SetBinError(i,min/10.);
    }
  }
}


// make systematic histogram
void checkMaximumSys(TH1F *hSys, TH1F *h, int opt=0,double minVal = 1)
{
  if (h->GetNbinsX()!=hSys->GetNbinsX()) {
    cout <<"ERROR! Different NBins in subroutine checkMaximumSys!"<<endl;
  } else {
    double val = minVal;
    for (int i=1;i<=h->GetNbinsX();i++) {
      //cout <<i<<" "<<val<<" "<<hSys->GetBinContent(i)<<" "<<h->GetBinContent(i)<<endl;
      if (h->GetBinContent(i)==0) continue;
      if (opt==0) val=minVal;
      if (fabs(hSys->GetBinContent(i))>val) val = fabs(hSys->GetBinContent(i));
      if (fabs(h->GetBinContent(i)-1)+1>val) val=fabs(h->GetBinContent(i)-1)+1;
      hSys->SetBinContent(i,val);
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
 

void dumpDatatoTxt(const char *centbin,TH1F *h, TH1F *hsys, TH1F *htotStat, const char *txtfile)
{
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
