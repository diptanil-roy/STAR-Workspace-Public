#include "levy.h"
void getfitsys_xx(const int errmode = 2)
{

  // gROOT->LoadMacro("~/work/functions/LevyFunction.C");

  const float    D0Mass   = 1.86483;
  const float    MKSize   = 1.7;

  TFile *f = new TFile("./all.root");//from Yifei// using the new pp data point
  TGraphAsymmErrors *gr_Dnew_err = (TGraphAsymmErrors *)f->Get("pp_Dnew_err");
  TGraphAsymmErrors *gr_Dnew_sys = (TGraphAsymmErrors *)f->Get("pp_Dnew_sys");
  f->Close();

  TFile *f2 = new TFile("./pp_D_TBWFit.root");//from yifei// BW fit for pp
  TF1* fcTBW = (TF1*)f2->Get("fcTBW")->Clone("fcTBW");
  fcTBW->SetName("fcTBW");
  f2->Close();

  TFile *f3 = new TFile("./B_and_D0_ptSpectra.root");//from Xiaolong// 
  TH1F* hD0_fonll = (TH1F*)f3->Get("hD0_fonll")->Clone("hD0_fonll");
  TH1F* hD0Min_fonll = (TH1F*)f3->Get("hD0Min_fonll")->Clone("hD0Min_fonll");
  TH1F* hD0Max_fonll = (TH1F*)f3->Get("hD0Max_fonll")->Clone("hD0Max_fonll");
  hD0_fonll->SetDirectory(0);
  hD0Min_fonll->SetDirectory(0);
  hD0Max_fonll->SetDirectory(0);
  //pp12
  TGraphErrors *gD0_pp_Run12_err = (TGraphErrors *)f3->Get("gD0_pp_Run12_err");//run12 w/o D*/D0
  TGraphErrors *gD0_pp_Run12_sys = (TGraphErrors *)f3->Get("gD0_pp_Run12_sys");//run12 w/o D*/D0
  TF1* fD0_pp_Run12 = (TF1*)f3->Get("fD0_pp_Run12")->Clone("fD0_pp_Run12");
  f3->Close();

  TFile *f33 = new TFile("./D0DstarRatio.root");//from Yifei.. Dstar over D0 ratio with pt dependece 
  TF1* funDDR = (TF1*)f33->Get("funDDR")->Clone("funDDR");
  f33->Close();

  for(int i=0; i<gD0_pp_Run12_err->GetN(); i++) {
    float x    = gD0_pp_Run12_err->GetX()[i];
    float y    = gD0_pp_Run12_err->GetY()[i];
    gD0_pp_Run12_err->GetY()[i] *= funDDR->Eval(x);
  }

  TFile *f4 = new TFile("./D_PYTHIA_yifei.root");//from yifei// BW fit for pp
  TH1F* hTrgDPt_PYT = (TH1F*)f4->Get("hTrgDPt")->Clone("hTrgDPt_PYT");
  hTrgDPt_PYT->SetDirectory(0);
  TH1F* hTrgD_PYT= (TH1F*)f4->Get("hTrgD")->Clone("hTrgD_PYT");
  hTrgD_PYT->SetDirectory(0);
  f4->Close();
  hTrgDPt_PYT->Scale(0.1);//abr scale
  hTrgD_PYT->Scale(0.1);//abr scale

  TFile *f5 = new TFile("./D_PYTHIA_mine.root");//from yifei// BW fit for pp
  TH1F* hD0Pt_PYT = (TH1F*)f5->Get("hD0Pt")->Clone("hD0Pt_PYT");
  hD0Pt_PYT->SetDirectory(0);
  TH1F* hD0_PYT= (TH1F*)f5->Get("hD0")->Clone("hD0_PYT");
  hD0_PYT->SetDirectory(0);
  f5->Close();
  hD0Pt_PYT->Scale(0.1);//abr scale
  hD0_PYT->Scale(0.1);//abr scale

  TGraphAsymmErrors *gr_Dmesonnew = (TGraphAsymmErrors *)gr_Dnew_err->Clone();
  gr_Dmesonnew->SetName("gr_Dmesonnew");
  TGraphErrors *grnew = new TGraphErrors(gr_Dmesonnew->GetN());
  grnew->SetName("grnew_Fitsys");
  
  const Int_t NCL = 500;
  TGraphErrors *grnewExtend = new TGraphErrors(NCL);// extend to 10 GeV
  grnewExtend->SetName("grnewExtend_Fitsys");
  for (int i=0; i<NCL; i++) grnewExtend->SetPoint(i, 0 + (float)i/50., 0);

  for(int i=0; i<gr_Dnew_err->GetN(); i++) {
    float x    = gr_Dnew_err->GetX()[i];
    float y    = gr_Dnew_err->GetY()[i];
    float err  = gr_Dnew_err->GetEYhigh()[i];
    float sysl = gr_Dnew_sys->GetEYlow()[i];
    float sysh = gr_Dnew_sys->GetEYhigh()[i];
    float el, eh;
    if(errmode == 0) { el = err; eh = err; }
    else if(errmode == 1) { el = sysl; eh = sysh; }
    else { el = sqrt(err*err + sysl*sysl); eh = sqrt(err*err + sysh*sysh); }
    gr_Dmesonnew->SetPointError(i,0.,0.,el,eh);
    grnew->SetPoint(i,x,y);
    grnew->SetPointError(i,0.,eh);
  }

  // TF1 *Dplnewfun = new TF1("Dplnewfun","[0]*(2*([2]-1.)*([2]-2.))/3.1415926/([2]-3.)**2/[1]**2/(1.0+2*x/[1]/([2]-3.0))**[2]",0.7,6.);
  TF1 *Dplnewfun = new TF1("Dplnewfun","[0]*(2*([2]-1.)*([2]-2.))/3.1415926/([2]-3.)**2/[1]**2/(1.0+2*x/[1]/([2]-3.0))**[2]",0,10.);
  Dplnewfun->SetParameters(0.13,1.1,20.0);

  TF1 *Levynew_pp = new TF1("Levynew_pp",LevyFcn,0.,10.,4);
  Levynew_pp->SetParameters(1e-2,17.751178,0.179449,D0Mass);
  Levynew_pp->FixParameter(3,D0Mass);
  Levynew_pp->SetLineWidth(2);
  Levynew_pp->SetLineStyle(7);
  Levynew_pp->SetLineColor(kCyan+1);

   TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   c1->SetLogy();
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);
   c1->SetTopMargin(0.02);
   c1->SetRightMargin(0.02);

   double x1 = 0.;
   // double x2 = 6.2;
   double x2 = 10.2;
   // double y1 = 5e-9*56;
   double y1 = 5e-10*56;
   double y2 = 5e-2*56;
   TH1 *d0 = new TH1D("d0","",1,x1,x2);
   d0->SetMinimum(y1);
   d0->SetMaximum(y2);
   d0->GetXaxis()->SetNdivisions(208);
   d0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   d0->GetXaxis()->SetTitleOffset(1);
   d0->GetXaxis()->SetTitleSize(0.065);
   d0->GetXaxis()->SetLabelOffset(0.01);
   d0->GetXaxis()->SetLabelSize(0.045);
   d0->GetXaxis()->SetLabelFont(42);
   d0->GetXaxis()->SetTitleFont(42);
   d0->GetYaxis()->SetNdivisions(203);
   // d0->GetYaxis()->SetTitle("(d^{2}#sigma^{c#bar{c}})/(2#pip_{T}dp_{T}dy) [mb/(GeV/c)^{2}]");
   d0->GetYaxis()->SetTitle("(d^{2}#sigma^{D^{0}})/(2#pip_{T}dp_{T}dy) [mb/(GeV/c)^{2}]");
   d0->GetYaxis()->SetTitleOffset(1.1);
   d0->GetYaxis()->SetTitleSize(0.065);
   d0->GetYaxis()->SetLabelOffset(0.005);
   d0->GetYaxis()->SetLabelSize(0.045);
   d0->GetYaxis()->SetLabelFont(42);
   d0->GetYaxis()->SetTitleFont(42);
   d0->Draw();

   TLine *l1 = new TLine(x1,y1,x2,y1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(x1,y2,x2,y2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(x1,y1,x1,y2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(x2,y1,x2,y2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   const float pteL = 0.1;
   const float yeL  = 1.1;

   // gr_Dmesonnew->Fit(Dplnewfun,"NOR");
   gr_Dmesonnew->Fit(Dplnewfun,"NO","",0.3,6);
   TFitResultPtr rplnew = gr_Dmesonnew->Fit(Dplnewfun,"S");
   // cout << "CHECK : chi2  "<< rplnew->Chi2() << " , ndf "<< rplnew->Ndf()  <<endl;
   Dplnewfun->SetLineWidth(2);
   Dplnewfun->SetLineColor(2);
   Dplnewfun->SetLineStyle(1);
   // Dplnewfun->SetRange(0.3,6.);
   // Dplnewfun->SetRange(0.,10.);
   Dplnewfun->Draw("same");

   gr_Dmesonnew->Fit(Levynew_pp,"NOR");
   Levynew_pp->SetLineColor(2);
   Levynew_pp->Draw("same");
   // TFitResultPtr rlvnew = gr_Dmesonnew->Fit(Levynew_pp,"S");

   TFitResultPtr rlvnew = grnew->Fit(Levynew_pp,"S");

   cout << "CHECK : chi2  "<< rlvnew->Chi2() << " , ndf "<< rlvnew->Ndf()  <<endl;
   int tmpn = grnew->GetN();
   double *tmpx = grnew->GetX();
   double tmpci[NCL]={0};
   rlvnew->GetConfidenceIntervals(tmpn, 1, 1, tmpx, tmpci, 0.6827,false);
   // for(int i =0;i<NCL;i++) cout << tmpci[i] << " , "; cout << endl;
   // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grnew,0.68);
   // vector< double > aa ; aa = rlvnew->GetConfidenceIntervals(0.68,false);
   // (ROOT::Fit::FitResult)rlvnew->GetConfidenceIntervals(0.68,false);
   // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grnew,0.68);
   //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(grnew);
   //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(grnew,0.9);
   grnew->SetFillColor(kRed-8);
   // grnew->Draw("e3same");

   // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grnewExtend,0.68);
   // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grnewExtend,0.68);
   int tmpn = grnewExtend->GetN();
   double *tmpx = grnewExtend->GetX();
   double tmpci[NCL]={0};
   double tmpci2[NCL]={0};//pt>6GeV
   rlvnew->GetConfidenceIntervals(tmpn, 1, 1, tmpx, tmpci2, 0.6827,true);
   rlvnew->GetConfidenceIntervals(tmpn, 1, 1, tmpx, tmpci, 0.6827,false);
   cout << " HERE start chercking ........ " << endl; 
   for(int i =0;i<NCL;i++) cout << tmpci[i]/Levynew_pp->Eval(grnewExtend->GetX()[i]) << " , "; cout << endl;
   cout << " FINISH chercking  !!!!! " << endl; 
  for (int i=0; i<NCL; i++) grnewExtend->SetPoint(i, 0 + (float)i/50., Levynew_pp->Eval(0 + (float)i/50));
  for (int i=0; i<NCL; i++) grnewExtend->SetPointError(i, 0, tmpci[i]);
  // // special pt<1 || pt>6 gev, using normalized error
  for (int i=NCL/10.*6.0; i<NCL; i++) grnewExtend->SetPointError(i, 0, tmpci2[i]);
  // for (int i=NCL/10.*5.5; i<NCL/10.*6.0; i++) grnewExtend->SetPointError(i, 0, 0.5* (tmpci[i] + tmpci2[i]));
  for (int i=NCL/10.*5; i<NCL/10.*6.0; i++) grnewExtend->SetPointError(i, 0, 0.5* (tmpci[i] + tmpci2[i]));
  // // for (int i =0 ; i<NCL/10.*0.6; i++) grnewExtend->SetPointError(i, 0, tmpci2[i]);
   grnewExtend->SetFillColor(kRed-8);
   grnewExtend->Draw("e3same");

   grnew->SetFillColor(kRed);
   // grnew->Draw("e3same");

   for(int i=0; i<gr_Dnew_sys->GetN(); i++) {
     float pt  = gr_Dnew_sys->GetX()[i];
     float y   = gr_Dnew_sys->GetY()[i];
     float yel = gr_Dnew_sys->GetEYlow()[i];
     float yeh = gr_Dnew_sys->GetEYhigh()[i];

     TLine *ll = new TLine(pt-pteL,y+yeh,pt+pteL,y+yeh);
     ll->SetLineWidth(2);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt-pteL,y-yel,pt+pteL,y-yel);
     ll->SetLineWidth(2);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt-pteL,(y+yeh)/yeL,pt-pteL,y+yeh);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt+pteL,(y+yeh)/yeL,pt+pteL,y+yeh);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt-pteL,(y-yel)*yeL,pt-pteL,y-yel);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt+pteL,(y-yel)*yeL,pt+pteL,y-yel);
     ll->SetLineColor(2);
     ll->Draw();
   }

   gr_Dmesonnew->SetLineWidth(2);
   gr_Dmesonnew->SetLineColor(2);
   gr_Dmesonnew->SetMarkerStyle(21);
   gr_Dmesonnew->SetMarkerSize(MKSize);
   gr_Dmesonnew->SetMarkerColor(2);
   gr_Dmesonnew->Draw("psame");

   Levynew_pp->Draw("same");
   Dplnewfun->Draw("same");

   fcTBW->Draw("same");

   hD0_fonll->SetLineColor(hD0Min_fonll->GetLineColor());
   hD0Max_fonll->SetLineColor(hD0Min_fonll->GetLineColor());
   hD0_fonll->Draw("same");
   hD0Max_fonll->Draw("same");
   hD0Min_fonll->Draw("same");

   gD0_pp_Run12_err->SetMarkerStyle(24);
   gD0_pp_Run12_err->Draw("Psame");
   // fD0_pp_Run12->Draw("same");

   // hTrgDPt_PYT->Draw("same");
   // hTrgD_PYT->Draw("same");
   
   // hD0Pt_PYT->Draw("same");
   hD0_PYT->Draw("same");

  TLegend* legend = new TLegend(0.2,0.2,0.42,0.43);
  legend->SetFillStyle(0);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->SetTextFont(132);
  legend->AddEntry(gr_Dmesonnew,"run09 w D*/D^{0}","pl");
  legend->AddEntry(gD0_pp_Run12_err,"run12 ","pl");
  legend->Draw();

  TLegend* legend = new TLegend(0.6,0.6,0.92,0.83);
  legend->SetFillStyle(0);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->SetTextFont(132);
  legend->SetNColumns(2);
  legend->AddEntry(Levynew_pp,"Levy","pl");
  legend->AddEntry(fcTBW,"TBW","pl");
  legend->AddEntry(hD0_PYT,"PYT","pl");
  legend->AddEntry(Dplnewfun,"pow-law","pl");
  legend->AddEntry(hD0_fonll,"fonll","pl");
  legend->Draw("same");

  c1->SaveAs("pp_baseLine1.eps");

  // return;

  c1->cd();
  c1->Update();

  TGraphErrors *grnewExtendRatio = (TGraphErrors*)grnewExtend->Clone("grnewExtendRatio");
  grnewExtendRatio->SetName("grnewExtendRatio_sys");
  TGraphAsymmErrors *gr_DmesonnewRatio = (TGraphAsymmErrors *)gr_Dmesonnew->Clone("gr_DmesonnewRatio");
  gr_DmesonnewRatio->SetName("gr_DmesonnewRatio");
  TGraphErrors *gD0_pp_Run12_errRatio = (TGraphErrors *)gD0_pp_Run12_err->Clone("gD0_pp_Run12_errRatio");
  gD0_pp_Run12_errRatio->SetName("gD0_pp_Run12_errRatio");

   for(int i=0; i<grnewExtendRatio->GetN(); i++)
   {
     float x ; x = grnewExtendRatio->GetX()[i];
     float y ; y = grnewExtendRatio->GetY()[i];
     grnewExtendRatio->GetY()[i] /= y;
     grnewExtendRatio->GetEY()[i] /= y;
   }
   cout << " HERE start chercking ........ " << endl; 
   for(int i=0; i<grnewExtendRatio->GetN(); i++)
   {
     cout << grnewExtendRatio->GetEY()[i] << " , ";
   }
   cout << endl << " FINISH chercking  !!!!! " << endl; 

   for(int i=0; i<gr_DmesonnewRatio->GetN(); i++)
   {
     float x  ; x  = gr_DmesonnewRatio->GetX()[i];
     float y0 ; y0 = Levynew_pp->Eval(x);
     float y  ; y  = gr_DmesonnewRatio->GetY()[i];
     float yl ; yl = gr_DmesonnewRatio->GetEYlow()[i];
     float yh ; yh = gr_DmesonnewRatio->GetEYhigh()[i];
     gr_DmesonnewRatio->GetY()[i]/=y0;
     gr_DmesonnewRatio->GetEYlow()[i]/=y0;
     gr_DmesonnewRatio->GetEYhigh()[i]/=y0;
   }

   for(int i=0; i<gD0_pp_Run12_errRatio->GetN(); i++)
   {
     float x  ; x  = gD0_pp_Run12_errRatio->GetX()[i];
     float y0 ; y0 = Levynew_pp->Eval(x);
     float y  ; y  = gD0_pp_Run12_errRatio->GetY()[i];
     float ey ; ey = gD0_pp_Run12_errRatio->GetEY()[i];
     gD0_pp_Run12_errRatio->GetY()[i]/=y0;
     gD0_pp_Run12_errRatio->GetEY()[i]/=y0;
   }

   double x1 = 0.;
   double x2 = 10.2;
   double y1 = -0.4;
   double y2 = 2.1;
   // double y1 = 0.1;
   // double y2 = 1.9;
   TH1 *d00 = new TH1D("d00","",1,x1,x2);
   d00->SetMinimum(y1);
   d00->SetMaximum(y2);
   d00->GetXaxis()->SetNdivisions(208);
   d00->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   d00->GetXaxis()->SetTitleOffset(1);
   d00->GetXaxis()->SetTitleSize(0.065);
   d00->GetXaxis()->SetLabelOffset(0.01);
   d00->GetXaxis()->SetLabelSize(0.045);
   d00->GetXaxis()->SetLabelFont(42);
   d00->GetXaxis()->SetTitleFont(42);
   d00->GetYaxis()->SetNdivisions(203);
   d00->GetYaxis()->SetTitle("Ratio");
   d00->GetYaxis()->SetTitleOffset(1.1);
   d00->GetYaxis()->SetTitleSize(0.065);
   d00->GetYaxis()->SetLabelOffset(0.005);
   d00->GetYaxis()->SetLabelSize(0.045);
   d00->GetYaxis()->SetLabelFont(42);
   d00->GetYaxis()->SetTitleFont(42);
   d00->Draw();

   TLine *l1 = new TLine(x1,y1,x2,y1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(x1,y2,x2,y2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(x1,y1,x1,y2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(x2,y1,x2,y2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   c1->SetLogy(0);

   d00->Draw();
   grnewExtendRatio->Draw("Psame");
   gr_DmesonnewRatio->Draw("Psame");
   gD0_pp_Run12_errRatio->Draw("Psame");
   grnewExtendRatio->GetYaxis()->SetRangeUser(-2,2);
   gr_DmesonnewRatio->GetYaxis()->SetRangeUser(-2,2);
   gD0_pp_Run12_errRatio->GetYaxis()->SetRangeUser(-2,2);

   TLine *l4 = new TLine(x1,1,x2,1);
   l4->SetLineWidth(2);
   l4->SetLineStyle(2);
   l4->SetLineColor(2);
   l4->Draw("same");

   TLine *l4 = new TLine(x1,1.25,x2,1.25);
   l4->SetLineWidth(2);
   l4->SetLineStyle(2);
   l4->SetLineColor(2);
   l4->Draw("same");

   TLine *l4 = new TLine(x1,0.75,x2,0.75);
   l4->SetLineWidth(2);
   l4->SetLineStyle(2);
   l4->SetLineColor(2);
   l4->Draw("same");

   TLatex *latex = new TLatex(0.2,0.85,"Levy Fitting Band");
   latex->SetNDC();
   latex->SetTextFont(42);
   latex->SetTextSize(0.04);
   latex->SetTextColor(1);
   latex->Draw("same");

   c1->SaveAs("pp_baseLine_Ratio.eps");

   // return;

   c1->cd();
   c1->Clear();
   c1->Update();
   c1->SetLogy(1);
   d0->Draw();

   TLine *l1 = new TLine(x1,y1,x2,y1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(x1,y2,x2,y2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(x1,y1,x1,y2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(x2,y1,x2,y2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   Dplnewfun->Draw("same");
   Levynew_pp->Draw("same");
   grnewExtend->Draw("e3same");

   for(int i=0; i<gr_Dnew_sys->GetN(); i++) {
     float pt  = gr_Dnew_sys->GetX()[i];
     float y   = gr_Dnew_sys->GetY()[i];
     float yel = gr_Dnew_sys->GetEYlow()[i];
     float yeh = gr_Dnew_sys->GetEYhigh()[i];

     TLine *ll = new TLine(pt-pteL,y+yeh,pt+pteL,y+yeh);
     ll->SetLineWidth(2);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt-pteL,y-yel,pt+pteL,y-yel);
     ll->SetLineWidth(2);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt-pteL,(y+yeh)/yeL,pt-pteL,y+yeh);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt+pteL,(y+yeh)/yeL,pt+pteL,y+yeh);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt-pteL,(y-yel)*yeL,pt-pteL,y-yel);
     ll->SetLineColor(2);
     ll->Draw();
     ll = new TLine(pt+pteL,(y-yel)*yeL,pt+pteL,y-yel);
     ll->SetLineColor(2);
     ll->Draw();
   }

   gr_Dmesonnew->Draw("psame");

   Levynew_pp->Draw("same");

   fcTBW->Draw("same");

   Dplnewfun->Draw("same");
   // hD0_fonll->Draw("same");
   // hD0Max_fonll->Draw("same");
   // hD0Min_fonll->Draw("same");

   // gD0_pp_Run12_err->SetMarkerStyle(24);
   // gD0_pp_Run12_err->Draw("Psame");
   // fD0_pp_Run12->Draw("same");

   // hTrgDPt_PYT->Draw("same");
   // hTrgD_PYT->Draw("same");
   
   // hD0Pt_PYT->Draw("same");
   // hD0_PYT->Draw("same");

  TLegend* legend = new TLegend(0.2,0.2,0.42,0.43);
  legend->SetFillStyle(0);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->SetTextFont(132);
  legend->AddEntry(gr_Dmesonnew,"run09 w D*/D^{0}","pl");
  // legend->AddEntry(gD0_pp_Run12_err,"run12 w/o","pl");
  legend->Draw();

  TLegend* legend = new TLegend(0.6,0.6,0.92,0.83);
  legend->SetFillStyle(0);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->SetTextFont(132);
  legend->SetNColumns(2);
  legend->AddEntry(Levynew_pp,"Levy","pl");
  legend->AddEntry(fcTBW,"TBW","pl");
  // legend->AddEntry(hD0_PYT,"PYT","pl");
  legend->AddEntry(Dplnewfun,"pow-law","pl");
  // legend->AddEntry(hD0_fonll,"fonll","pl");
  legend->Draw("same");

  c1->SaveAs("pp_baseLine2.eps");

  // return;

  //next is multiply pt for plotting
  //
  TF1 *DplnewfunPt = new TF1("DplnewfunPt","[0]*(2*([2]-1.)*([2]-2.))/3.1415926/([2]-3.)**2/[1]**2/(1.0+2*x/[1]/([2]-3.0))**[2]*x",0,10.);
  DplnewfunPt->SetParameters(Dplnewfun->GetParameters());
  DplnewfunPt->SetLineColor(Dplnewfun->GetLineColor());
  DplnewfunPt->SetLineStyle(Dplnewfun->GetLineStyle());
  TF1 *Levynew_ppPt = new TF1("Levynew_ppPt",LevyFcnPt,0.,10.,4);
  Levynew_ppPt->SetParameters(Levynew_pp->GetParameters());
  Levynew_ppPt->SetLineColor(Levynew_pp->GetLineColor());
  Levynew_ppPt->SetLineStyle(Levynew_pp->GetLineStyle());
  // TF1* fcTBWPt = new TF1("fcTBWPt","fcTBW*x");
  // TF1* fcTBWPt = new TF1("fcTBWPt","",0,10);
  // fcTBWPt->SetParameters(fcTBW->GetParameters());
  // fcTBWPt->SetLineColor(fcTBW->GetLineColor());
  // fcTBWPt->SetLineStyle(fcTBW->GetLineStyle());
  TGraphErrors *grnewExtendPt = new TGraphErrors(NCL);
  for(int i=0; i<grnewExtend->GetN(); i++)
  {
    grnewExtendPt->SetPoint(i, grnewExtend->GetX()[i], grnewExtend->GetY()[i] * grnewExtend->GetX()[i]) ;
    grnewExtendPt->SetPointError(i, 0, grnewExtend->GetEY()[i] * grnewExtend->GetX()[i]) ;
  }
  grnewExtendPt->SetMarkerColor(grnewExtend->GetMarkerColor());
  grnewExtendPt->SetMarkerSize(grnewExtend->GetMarkerSize());
  grnewExtendPt->SetLineColor(grnewExtend->GetLineColor());
  grnewExtendPt->SetLineStyle(grnewExtend->GetLineStyle());
  grnewExtendPt->SetFillColor(grnewExtend->GetFillColor());
  grnewExtendPt->SetFillStyle(grnewExtend->GetFillStyle());

  TGraphAsymmErrors *gr_DmesonnewPt = new TGraphAsymmErrors(NCL);
  for(int i=0; i<gr_Dmesonnew->GetN(); i++)
  {
    gr_DmesonnewPt->SetPoint(i, gr_Dmesonnew->GetX()[i], gr_Dmesonnew->GetY()[i] * gr_Dmesonnew->GetX()[i]) ;
    gr_DmesonnewPt->SetPointError(i, 0,0, gr_Dmesonnew->GetEYlow()[i] * gr_Dmesonnew->GetX()[i], gr_Dmesonnew->GetEYhigh()[i] * gr_Dmesonnew->GetX()[i]) ;
  }
  gr_DmesonnewPt->SetLineColor(gr_Dmesonnew->GetLineColor());
  gr_DmesonnewPt->SetLineStyle(gr_Dmesonnew->GetLineStyle());
  gr_DmesonnewPt->SetMarkerStyle(gr_Dmesonnew->GetMarkerStyle());
  gr_DmesonnewPt->SetMarkerColor(gr_Dmesonnew->GetMarkerColor());
  gr_DmesonnewPt->SetMarkerSize(gr_Dmesonnew->GetMarkerSize());

  TGraphErrors *gD0_pp_Run12_errPt = new TGraphErrors(gD0_pp_Run12_err->GetN());
  for(int i=0; i<gD0_pp_Run12_err->GetN(); i++) {
    float x    = gD0_pp_Run12_err->GetX()[i];
    float y    = gD0_pp_Run12_err->GetY()[i] * x;
    float yerr    = gD0_pp_Run12_err->GetEY()[i] * x;
    gD0_pp_Run12_errPt->SetPoint(i,x,y);
    gD0_pp_Run12_errPt->SetPointError(i,0,yerr);
  }
  gD0_pp_Run12_errPt->SetLineColor(gD0_pp_Run12_err->GetLineColor());
  gD0_pp_Run12_errPt->SetLineStyle(gD0_pp_Run12_err->GetLineStyle());
  gD0_pp_Run12_errPt->SetMarkerStyle(gD0_pp_Run12_err->GetMarkerStyle());
  gD0_pp_Run12_errPt->SetMarkerColor(gD0_pp_Run12_err->GetMarkerColor());
  gD0_pp_Run12_errPt->SetMarkerSize(gD0_pp_Run12_err->GetMarkerSize());

   double x1 = 0.;
   double x2 = 6.2;
   double y1 = 5e-10;
   double y2 = 2.2e-2;

   c1->cd();
   c1->Clear();
   c1->Update();
   c1->SetLogy(0);
   TH1 *d01 = new TH1D("d01","",1,x1,x2);
   d01->SetMinimum(y1);
   d01->SetMaximum(y2);
   d01->GetXaxis()->SetNdivisions(208);
   d01->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   d01->GetXaxis()->SetTitleOffset(1);
   d01->GetXaxis()->SetTitleSize(0.065);
   d01->GetXaxis()->SetLabelOffset(0.01);
   d01->GetXaxis()->SetLabelSize(0.045);
   d01->GetXaxis()->SetLabelFont(42);
   d01->GetXaxis()->SetTitleFont(42);
   d01->GetYaxis()->SetNdivisions(203);
   d01->GetYaxis()->SetTitle("(d^{2}#sigma^{D^{0}})/(2#pidp_{T}dy) [mb/(GeV/c)^{1}]");
   d01->GetYaxis()->SetTitleOffset(1.1);
   d01->GetYaxis()->SetTitleSize(0.065);
   d01->GetYaxis()->SetLabelOffset(0.005);
   d01->GetYaxis()->SetLabelSize(0.045);
   d01->GetYaxis()->SetLabelFont(42);
   d01->GetYaxis()->SetTitleFont(42);
   d01->Draw();

   TLine *l1 = new TLine(x1,y1,x2,y1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(x1,y2,x2,y2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(x1,y1,x1,y2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(x2,y1,x2,y2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   grnewExtendPt->Draw("e3same");
   DplnewfunPt->Draw("same");
   Levynew_ppPt->Draw("same");

   gr_DmesonnewPt->Draw("psame");
   gD0_pp_Run12_errPt->Draw("psame");

   // fcTBWPt->Draw("same");

  TLegend* legend = new TLegend(0.6,0.3,0.82,0.53);
  legend->SetFillStyle(0);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->SetTextFont(132);
  legend->AddEntry(gr_DmesonnewPt,"run09 w D*/D^{0}","pl");
  legend->AddEntry(gD0_pp_Run12_errPt,"run12 w/o","pl");
  legend->Draw();

  TLegend* legend = new TLegend(0.6,0.6,0.92,0.83);
  legend->SetFillStyle(0);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->SetTextFont(132);
  legend->SetNColumns(2);
  legend->AddEntry(Levynew_ppPt,"LevyPt","pl");
  // legend->AddEntry(fcTBWPt,"BWPt","pl");
  legend->AddEntry(DplnewfunPt,"pow-lawPt","pl");
  legend->Draw("same");

  c1->SaveAs("pp_baseLine3Pt.eps");

  // return;

  //This is the one going to use for the final sys. uncertainty calculation
  TGraphErrors *gr_pp_sys= new TGraphErrors(NCL);
  gr_pp_sys->SetName("gr_pp_sys");
  TGraphAsymmErrors *gr_pp_sysAsy= new TGraphAsymmErrors(NCL);
  gr_pp_sysAsy->SetName("gr_pp_sysAsy");

  //This is the one use for the final relative sys. uncertainty calculation, divided to levy baseline
  TGraphErrors *gr_pp_sysRatio= new TGraphErrors(NCL);
  gr_pp_sysRatio->SetName("gr_pp_sysRatio");

  TGraphAsymmErrors *gr_pp_sysAsyRatio= new TGraphAsymmErrors(NCL);
  gr_pp_sysAsyRatio->SetName("gr_pp_sysAsyRatio");

  // for (int i=0; i<NCL; i++)
  // {
  //   float x    = grnewExtend->GetX()[i];
  //   float y    = grnewExtend->GetY()[i];
  //   float err  = grnewExtend->GetEY()[i];
  //   float err2 = fabs(y - fcTBW->Eval(x));
  //   float sys  = sqrt(err*err + err2*err2);
  //   gr_pp_sys->SetPoint(i, x, y);
  //   gr_pp_sys->SetPointError(i,0.,sys);
  // }

  for (int i=0; i<NCL; i++)
  {
    float x    = grnewExtend->GetX()[i];
    float y    = grnewExtend->GetY()[i];
    float err  = grnewExtend->GetEY()[i];
    float err2 = fabs(y - fcTBW->Eval(x));
    float err3 = fabs(y - Dplnewfun->Eval(x));

    err3/= sqrt(3.);

    float sys;
    float sysAsy1;
    float sysAsy2;

    float sysMax;
    //take maximum difference as sys
    float tmp = err > err2 ? err : err2;
    sysMax = tmp > err3 ? tmp : err3;

    //take as additional sys. source
    // //Ratio3 plots
    // sys = sqrt(pow(err,2)+pow(err2,2)+pow(err3,2));
    //
    // remove TBW part
    sys = sqrt(pow(err,2)+pow(err3,2));

    // if(i>= 50 && i<=300)
    // {
      sysAsy1 = sys;
      sysAsy2 = sys;
    // }
    if(i < 1.0 * 50) 
    {
      sysAsy1 = err;// pt <1 
      sysAsy2 = sys;
    }
    else if(i > 6. * 50) 
    {
      sysAsy1 = sys;//err3;//err;//sys;// pt >6 
      sysAsy2 = err;//err3;//err;// pt>6
    }

    //Next is kind of tricky, prefer to using consitent method for whole pt range as above
    //
    // if( x < gr_Dnew_err->GetX()[0])
    // {
    //   sys = err > err2 ? err : err2;
    // }
    // if( x > gr_Dnew_err->GetX()[0] && x < gr_Dnew_err->GetX()[5] )
    // {
    //   sys = err;
    // }
    // if( x > gr_Dnew_err->GetX()[5])
    // {
    //   sys = err > err2 ? err : err2;
    //   sys = err > err3 ? err : err3;
    // }

    gr_pp_sys->SetPoint(i, x, y);
    gr_pp_sys->SetPointError(i,0.,sys);

    gr_pp_sysAsy->SetPoint(i, x, y);
    gr_pp_sysAsy->SetPointError(i,0.,0,sysAsy1,sysAsy2);

    // for( int ii=0; ii < gr_Dnew_err->GetN(); ii++)
    // {
    //   if(fabs(x-gr_Dnew_err->GetX()[ii])<= 5./NCL) cout << x << " , " << y << " , " << err << " , " << sys << endl;
    // }
  }

  gr_Dnew_err->Print("ALL");

  // gr_pp_sys->Draw("e3same");
  gr_pp_sys->SetFillColor(kRed-8);
  gr_pp_sys->SetFillStyle(3020);

  gr_pp_sysAsy->SetFillColor(kRed-8);
  gr_pp_sysAsy->SetFillStyle(3020);

  for (int i=0; i<NCL; i++)
  {
    float x    = gr_pp_sys->GetX()[i];
    float y    = gr_pp_sys->GetY()[i];
    float err  = gr_pp_sys->GetEY()[i];
    float errlw  = gr_pp_sysAsy->GetEYlow()[i];
    float errup  = gr_pp_sysAsy->GetEYhigh()[i];

    gr_pp_sysRatio->SetPoint(i, x, 1.0);
    gr_pp_sysRatio->SetPointError(i,0.,err/y);

    gr_pp_sysAsyRatio->SetPoint(i, x, 1.0);
    gr_pp_sysAsyRatio->SetPointError(i,0.,0,errlw/y,errup/y);

  }


   double x1 = 0.;
   double x2 = 10.2;
   double y1 = -0.4;
   double y2 = 2.1;
   c1->cd();
   c1->SetLogy(0);
   TH1 *d00 = new TH1D("d00","",1,x1,x2);
   d00->SetMinimum(y1);
   d00->SetMaximum(y2);
   d00->GetXaxis()->SetNdivisions(208);
   d00->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   d00->GetXaxis()->SetTitleOffset(1);
   d00->GetXaxis()->SetTitleSize(0.065);
   d00->GetXaxis()->SetLabelOffset(0.01);
   d00->GetXaxis()->SetLabelSize(0.045);
   d00->GetXaxis()->SetLabelFont(42);
   d00->GetXaxis()->SetTitleFont(42);
   d00->GetYaxis()->SetNdivisions(203);
   d00->GetYaxis()->SetTitle("Ratio");
   d00->GetYaxis()->SetTitleOffset(1.1);
   d00->GetYaxis()->SetTitleSize(0.065);
   d00->GetYaxis()->SetLabelOffset(0.005);
   d00->GetYaxis()->SetLabelSize(0.045);
   d00->GetYaxis()->SetLabelFont(42);
   d00->GetYaxis()->SetTitleFont(42);
   d00->Draw();
   d00->SetMinimum(y1);//-2
   d00->SetMaximum(y2);//+3
   d00->Draw();
   gr_pp_sysRatio->Draw("Psame");
   // grnewExtendRatio->Draw("Psame");
   gr_DmesonnewRatio->Draw("Psame");
   gD0_pp_Run12_errRatio->Draw("Psame");

   TLine *l4 = new TLine(x1,1,x2,1);
   l4->SetLineWidth(2);
   l4->SetLineStyle(2);
   l4->SetLineColor(2);
   l4->Draw("same");

   TLine *l4 = new TLine(x1,1.25,x2,1.25);
   l4->SetLineWidth(2);
   l4->SetLineStyle(2);
   l4->SetLineColor(2);
   l4->Draw("same");

   TLine *l4 = new TLine(x1,0.75,x2,0.75);
   l4->SetLineWidth(2);
   l4->SetLineStyle(2);
   l4->SetLineColor(2);
   l4->Draw("same");

   TLatex *latex = new TLatex(0.2,0.85,"Including diff. fitFun to sys. Band");
   latex->SetNDC();
   latex->SetTextFont(42);
   latex->SetTextSize(0.04);
   latex->SetTextColor(1);
   latex->Draw("same");

   // c1->SaveAs("pp_baseLine_Ratio2.eps");
   c1->SaveAs("pp_baseLine_Ratio3.eps");

   TFile *fout = new TFile("out_ppsys.root","recreate");
   fout->cd();
   Dplnewfun->Write();
   Levynew_pp->Write();
   gr_Dnew_err->Write();
   gr_Dnew_sys->Write();
   // grnew->Write();
   grnewExtend->Write();
   fcTBW->Write();
   hD0_fonll->Write();
   hD0Min_fonll->Write();
   hD0Max_fonll->Write();
   hD0Pt_PYT->Write();
   gr_DmesonnewRatio->Write();
   gD0_pp_Run12_errRatio->Write();
   grnewExtendRatio->Write();
   gr_pp_sys->Write();
   gr_pp_sysAsy->Write();
   gr_pp_sysRatio->Write();
   gr_pp_sysAsyRatio->Write();
   fout->Close();

}


