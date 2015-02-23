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
void scaleToBinWidth(TH1 &h1){
  double nBins = h1.GetXaxis()->GetNbins();
  for (int iBin = h1.GetXaxis()->GetFirst(); iBin <= h1.GetXaxis()->GetLast(); iBin++){
    double iW = h1.GetBinWidth(iBin);
    double old = h1.GetBinContent(iBin);
    h1.SetBinContent(iBin,old/iW);
  }


}

