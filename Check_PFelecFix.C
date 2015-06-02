//
// Raghav Kunnawalkam Elayavalli
// Rutgers University, June 2nd 2015 
//

//
// Macro to load in 4 differnet hiForest files and plot the quantities that correspond to pf electron issue with the fix (3 differnet kinds of fixes) and compares them to each other. 
// We will be plotting 3 main quantities
// 1) eMax vs jetpt 
// 2) eMax/jetpt vs jtpt
// 3) eMax/(chSum+neSum+muSum+phSum) vs jtpt 
//


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <cstdlib>
#include <cmath>


static const int nbins_cent = 6;
static Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};

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

void Check_PFelecFix(int radius = 3,
		     char * algo = (char*)"PF",
		     char * bkgsub = (char*)"Pu"){

  
  TH1::SetDefaultSumw2();
  //gStyle->SetOptStat(0);
  
  TStopwatch timer;
  timer.Start();
  
  TDatime date;

  char * fileType[4][256] = {"badFile","Fix_1","Fix_2","Fix_3"};
  
  // get the input hiForest files, these are going to be in an array 0 - bad fine, 1 - fix_1, 2 - fix_2, 3 - fix_3; 
  TFile * fIn[4];

  // input the filen names
  fIn[0] = TFile::Open("");
  fIn[1] = TFile::Open("");
  fIn[2] = TFile::Open("");
  fIn[3] = TFile::Open("");
  
  // get the jet trees from the necessary files (these events already passed the event quality cuts so only hav to get the jet Tree).
  TTree * jet[4];
   
  // get the histograms from these files for the different centrality bins and the fixes. 
  TH1F * heMax[4][nbins_cent], * heSum[4][nbins_cent];
  TH2F * heMax_vs_jtpt[4][nbins_cent], * heMaxJtpt_vs_jtpt[4][nbins_cent], * heMaxSumcand_vs_jtpt[4][nbins_cent];


  for(int k = 0; k<4; ++k){

    jet[k] = (TTree*)fIn[k]->Get(Form("ak%s%d%sJetAnalyzer/t",bkgsub,radius,algo));

    for(int i = 0; i<nbins_cent; ++i){

      heMax[k][i] = new TH1F(Form("heMax_%s_cent%d",fileType,i),"",200,0,200);
      heMax_vs_jtpt[k][i] = new TH1F(Form("heMax_vs_jtpt_%s_cent0",fileType,i),"",400,0,400,100,0,200);
      heMaxJtpt_vs_jtpt[k][i] = new TH1F(Form("heMaxJtpt_vs_jtpt_%s_cent0",fileType,i),"",400,0,400,100,0,10);
      heMaxSumcand_vs_jtpt[k][i] = new TH1F(Form("heMaxSumcand_vs_jtpt_%s_cent0",fileType,i),"",400,0,400,100,0,10);

      jet[k]->Draw(Form("eMax>>heMax_%s_cent%d",fileType,i),Form("hiBin>=%d && hiBin<%d",5 * boundaries_cent[i], 5 * boundaries_cent[i+1]),"goff");
      jet[k]->Draw(Form("eMax:jtpt>>heMax_vs_jtpt_%s_cent%d",fileType,i),Form("hiBin>=%d && hiBin<%d",5 * boundaries_cent[i], 5 * boundaries_cent[i+1]),"goff");
      jet[k]->Draw(Form("eMax/jtpt:jtpt>>heMaxJtpt_vs_jtpt_%s_cent%d",fileType,i),Form("hiBin>=%d && hiBin<%d",5 * boundaries_cent[i], 5 * boundaries_cent[i+1]),"goff");
      jet[k]->Draw(Form("eMax/(chSum+neSum+muSum+phSum):jtpt>>heMaxSumcand_vs_jtpt_%s_cent%d",fileType,i),Form("hiBin>=%d && hiBin<%d",5 * boundaries_cent[i], 5 * boundaries_cent[i+1]),"goff");
    
    }

  }

  // now that we have the histograms, lets do the plotting part here.
  // we need one canvas for each plotted variable.
  // at the moment im only going to plot the most central events. 0 <= hiBin < 5

  TCanvas * ceMax, *ceMax_vs_jtpt, *ceMaxJtpt_vs_jtpt, *ceMaxSumcand_vs_jtpt;

  ceMax = new TCanvas("ceMax","",1200,1000);
  ceMax->Divide(4,1);
  ceMax_vs_jtpt = new TCanvas("ceMax_vs_jtpt","",1200,1000);
  ceMax_vs_jtpt->Divide(4,1);
  ceMaxJtpt_vs_jtpt = new TCanvas("ceMaxJtpt_vs_jtpt","",1200,1000);
  ceMaxJtpt_vs_jtpt->Divide(4,1);
  ceMaxSumcand_vs_jtpt = new TCanvas("ceMaxSumcand_vs_jtpt","",1200,1000);
  ceMaxSumcand_vs_jtpt->Divide(4,1);

  int cent = 0; // change this to draw other centrality classes 
  
  for(int k = 0; k<4; ++k){

    ceMax->cd(k+1);
    heMax[k][cent]->Draw();
    drawText(Form("%s %2.0f-%2.0f%",fileType,2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.2,0.2,14);

    ceMax_vs_jtpt->cd(k+1);
    heMax_vs_jtpt[k][cent]->Draw("colz");
    drawText(Form("%s %2.0f-%2.0f%",fileType,2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.2,0.2,14);

    ceMaxJtpt_vs_jtpt->cd(k+1);
    heMaxJtpt_vs_jtpt[k][cent]->Draw("colz");
    drawText(Form("%s %2.0f-%2.0f%",fileType,2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.2,0.2,14);

    ceMaxSumcand_vs_jtpt->cd(k+1);
    heMax_vs_jtpt[k][cent]->Draw("colz");
    drawText(Form("%s %2.0f-%2.0f%",fileType,2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),0.2,0.2,14);

  }

  ceMax->SaveAs(Form("PbPb_eMaxvariable_cent%d_ak%s%d%s_%d.pdf",cent,bkgsub,radius,algo,date.GetDate()),"RECREATE");
  ceMax_vs_jtpt->SaveAs(Form("PbPb_eMax_vs_jtpt_cent%d_ak%s%d%s_%d.pdf",cent,bkgsub,radius,algo,date.GetDate()),"RECREATE");
  ceMaxJtpt_vs_jtpt->SaveAs(Form("PbPb_eMaxOverJtpt_vs_jtpt_cent%d_ak%s%d%s_%d.pdf",cent,bkgsub,radius,algo,date.GetDate()),"RECREATE");
  ceMaxSumcand_vs_jtpt->SaveAs(Form("PbPb_eMaxOver_SumCandidates_without_eMax_vs_jtpt_cent%d_ak%s%d%s_%d.pdf",cent,bkgsub,radius,algo,date.GetDate()),"RECREATE");
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(float)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(float)timer.RealTime()/60<<endl;



}
