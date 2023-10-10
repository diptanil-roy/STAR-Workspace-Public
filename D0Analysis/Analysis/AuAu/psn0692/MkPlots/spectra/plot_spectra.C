#include "myFunction.h"
#include "myConst.h"
#include "StyleUtilities.h"

TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
   TLatex *latex = new TLatex(x, y, text);
   latex->SetNDC();
   latex->SetTextFont(textFont);
   latex->SetTextSize(textSize);
   latex->SetTextColor(colorIndex);
   latex->Draw("same");
   return latex;
}

TLine* drawLine(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor)
{
   TLine *l1 = new TLine(xlow, ylow, xup, yup);
   l1->SetLineWidth(lineWidth);
   l1->SetLineColor(lineColor);
   l1->Draw("same");
   return l1;
}

void drawLines(Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineColor)
{
   drawLine(xlow, ylow, xup, ylow, lineWidth, lineColor);
   drawLine(xlow, yup, xup, yup, lineWidth, lineColor);
   drawLine(xlow, ylow, xlow, yup, lineWidth, lineColor);
   drawLine(xup, ylow, xup, yup, lineWidth, lineColor);
}

void setpad(TPad *pad, float left, float right, float top, float bottom)
{
   pad->SetFillColor(10);
   pad->SetBorderMode(0);
   pad->SetBorderSize(0);
   pad->SetFrameFillColor(10);
   pad->SetFrameBorderMode(0);
   pad->SetFrameBorderSize(0);
   pad->SetLeftMargin(left);
   pad->SetRightMargin(right);
   pad->SetTopMargin(top);
   pad->SetBottomMargin(bottom);
}

void plot_spectra()
{
   globalSetting();
   char dir[250];
   char name[250];
   char title[250];
   char buf[1024];
   //TString CMD = 0;
   char CMD[250];
   TH1F* h0;

   sprintf(dir, "pic");
   sprintf(CMD, "[ -d %s ] || mkdir -p %s", dir, dir);
   // gSystem->Exec(CMD);

   // const int ncent = 6;
   // const char nameCent[ncent][250] = {"0-80%", "0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
   // const char nameCent1[ncent][250] = {"0-80% [#times10]", "0-10% [#times1]", "10-20% [/3]", "20-40% [/6]", "40-60% [/12]", "60-80% [/24]"};
   // const char nameCentXL[ncent][250] = {"0_80", "0_10", "10_20", "20_40", "40_60", "60_80"};
   // const float scale[ncent] = {10., 1., 1 / 3., 1 / 6., 1. / 12, 1. / 24};
   
   const int ncent = 5;
   const char nameCent1[ncent][250] = {"0-10%(/1)", "10-20%(/2)", "20-40%(/4)", "40-60%(/8)", "60-80%(/16)"};
   const char nameCentXL[ncent][250] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
   const float scale[ncent] = {1., 1 / 2., 1 / 4., 1. / 8, 1. / 16};

   // const int ncent = 5;
   // const char nameCent[ncent][250] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
   // const char nameCent1[ncent][250] = {"0-10%", "10-20% [/3]", "20-40% [/6]", "40-60% [/12]", "60-80% [/24]"};
   // const char nameCentXL[ncent][250] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
   // const float scale[ncent] = {1., 1/3., 1/6., 1./12, 1./24};


   // const int ncent = 4;
   // const char nameCent[ncent][250] = {"0-10%", "10-40%", "40-80%", "0-80%"};
   // const char nameCent1[ncent][250] = {"0-10%", "10-40% [/3]", "40-80% [/12]", "0-80% [/24]"};
   // const char nameCentXL[ncent][250] = {"0_10", "10_40", "40_80", "0_80"};
   // const float scale[ncent] = {1., 1/3., 1./12, 1./24};

   //Read spectra
   //1. from xiaolong
   TGraphErrors* gD0err_xl[ncent];
   TGraphErrors* gD0sys_xl[ncent];
   TF1* fLevy[ncent];
   TFile* fin1 = new TFile("D0_Spectra_Run14HFT.root");
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_err_%s", nameCentXL[icent]));
      gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_sys_%s", nameCentXL[icent]));
      fLevy[icent] = (TF1*)fin1->Get(Form("flevy_%s", nameCentXL[icent]));
   }
   fin1->Close();

   //scale
   for (int icent = 0; icent < ncent; icent++)
   {
      ScaleGraph(gD0err_xl[icent], scale[icent]);
      ScaleGraph(gD0sys_xl[icent], scale[icent]);
      fLevy[icent]->SetParameter(0, fLevy[icent]->GetParameter(0)*scale[icent]);
   }

   //Set for Draw
   float markerSize = 2.0;
   float lineWidth = 2;
   float markerSizeScale[ncent] = {0.75, 1.,1., 1., 1.2};
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent]->SetMarkerStyle(MARKERSTYLE[icent]);
      gD0err_xl[icent]->SetMarkerSize(markerSize*markerSizeScale[icent]);

      gD0err_xl[icent]->SetMarkerColor(COLOR[icent]);
      gD0err_xl[icent]->SetLineWidth(lineWidth);
      gD0err_xl[icent]->SetLineColor(COLOR[icent]);

      fLevy[icent]->SetLineColor(COLOR[icent]);
      fLevy[icent]->SetLineWidth(lineWidth);
      fLevy[icent]->SetLineStyle(7);
   }

   //plot
   // TCanvas* c1 = new TCanvas("c1", "A Canvas", 10, 10, 600, 800);
   // setPad(c1);
   // c1->SetGridx(0);
   // c1->SetGridy(0);

   // TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,550,550);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   c1->SetLogy();
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);

   // c1->SetGridx();
   // c1->SetGridy();

   c1->SetLeftMargin(0.14);
   c1->SetBottomMargin(0.15);
   c1->SetTopMargin(0.02);
   c1->SetRightMargin(0.03);

   // setPad(c1, 0.13, 0.1, 0.06, 0.13, 0, 1);
   // setGraphicsStyle();
   // TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 500, 600);
   // gStyle->SetOptFit(0);
   // gStyle->SetOptStat(0);
   // gStyle->SetEndErrorSize(2);
   // gStyle->SetTitle(0);

   // setpad(c1, 0.14, 0.02, 0.01, 0.12);
   double x1 = 0;
   // double x2 = 9.1;
   double x2 = 12.1;
   // double y1 = 1.0e-8;
   double y1 = 1.1e-9;
   double y2 = .3e0;
   c1->cd(1)->SetLogy(1);
   // gPad->SetPad(0., 0.0, 1., 1.);
   // setpad(gPad,0.13, 0.1, 0.06, 0.13);

   h0 = new TH1F("", "", 1, x1, x2);
   // h0= new TH1F("","",1,0,10);
   h0->Draw();
   // h0->GetXaxis()->CenterTitle();
   // setHisto(h0, "", "p_{T} (GeV/c)", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
   h0->SetMinimum(y1);
   h0->SetMaximum(y2);
   h0->GetXaxis()->SetNdivisions(208);
   h0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h0->GetXaxis()->SetTitleOffset(1.0);
   h0->GetXaxis()->SetTitleSize(0.057);
   h0->GetXaxis()->SetLabelOffset(0.01);
   h0->GetXaxis()->SetLabelSize(0.045);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetYaxis()->SetNdivisions(204);
   h0->GetYaxis()->SetTitle("d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
   h0->GetYaxis()->SetTitleOffset(1.24);
   h0->GetYaxis()->SetTitleSize(0.053);
   h0->GetYaxis()->SetLabelOffset(0.005);
   h0->GetYaxis()->SetLabelSize(0.045);
   h0->GetYaxis()->SetLabelFont(42);
   h0->GetYaxis()->SetTitleFont(42);
   h0->SetLineWidth(2);
   h0->Draw();
   for (int icent = 0; icent < ncent; icent++)
   {
      // if( icent !=2 ) fLevy[icent]->Draw("same");
      fLevy[icent]->Draw("same");
      fLevy[icent]->SetRange(0,9);
      fLevy[icent]->SetLineColor(1);
      fLevy[icent]->SetLineStyle(icent%2+1);
      // gD0err_xl[icent]->Draw("psame");
     int kMarkerStyle = 20;
     if(icent%2==1) kMarkerStyle = 24;
     gD0err_xl[icent]->SetMarkerStyle(kMarkerStyle);
     gD0err_xl[icent]->SetMarkerColor(1);
     gD0err_xl[icent]->SetLineColor(1);
     gD0err_xl[icent]->SetMarkerSize(1.5);
     gD0err_xl[icent]->SetLineWidth(2);
     gD0err_xl[icent]->Draw("psame");
      //draw systematic error
      const float sysw = 0.15;
      for (int i = 0; i < gD0sys_xl[icent]->GetN(); i++)
      {
         const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(gD0err_xl[icent]->GetLineColor());
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(gD0err_xl[icent]->GetLineColor());
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(gD0err_xl[icent]->GetLineColor());
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(gD0err_xl[icent]->GetLineColor());
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(gD0err_xl[icent]->GetLineColor());
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(gD0err_xl[icent]->GetLineColor());
         lv->Draw("same");
      }
   }
   drawLines(x1, y1, x2, y2, 3, 1);

   // sprintf(buf, "STAR Au+Au @ 200 GeV");
   sprintf(buf, "Au+Au #sqrt{s_{NN}} = 200 GeV");
   drawLatex(0.48, 0.88, buf, 42, 0.05, 1);

  TLatex *tex = new TLatex(7.6, y2/100, "D^{0} |y| < 1");
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->Draw("same");

   for(int icent=0;icent<ncent;icent++) {
     if(icent == 0 ) 
     {
       TLatex *tex = new TLatex(9.2, gD0sys_xl[icent]->GetY()[10]*0.9, Form("%s",nameCent1[icent]));
     }
     else if(icent >0 && icent < ncent -1 ) 
     {
       TLatex *tex = new TLatex(9.2, gD0sys_xl[icent]->GetY()[10]*0.6, Form("%s",nameCent1[icent]));
     }
     else if(icent == ncent -1 ) 
     {
       TLatex *tex = new TLatex(9.2, gD0sys_xl[icent]->GetY()[9]*0.1, Form("%s",nameCent1[icent]));
     }
     tex->SetTextFont(42);
     tex->SetTextSize(0.04);
     tex->Draw("same");
   }

   TLine *l1 = new TLine(1.3, 1.5e-8, 1.8, 1.5e-8);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineColor(1);
   l1->Draw("same");
   TLine *l1 = new TLine(1.3, 1.e-8, 1.8, 1.e-8);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");
  TLatex *tex = new TLatex(2.0, 1.0e-8, "Levy Fit");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->Draw("same");


   gPad->SetLogy();
   c1->SaveAs("D0_spectra.eps");
   c1->SaveAs("D0_spectra.pdf");
   c1->SaveAs("D0_spectra.gif");
   c1->SaveAs("D0_spectra.png");
}
