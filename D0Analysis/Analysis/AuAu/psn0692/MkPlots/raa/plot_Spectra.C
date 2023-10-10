#include "myFunction.h"
#include "myConst.h"
#include "levy.h"
// #include "StyleUtilities.h"

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

void plot_Spectra()
{
   globalSetting();
   char dir[250];
   char name[250];
   char title[250];
   char buf[1024];
   //TString CMD = 0;
   char CMD[250];
   TLegend* legend1;
   TLegend* legend2;

   sprintf(dir, "pic");
   sprintf(CMD, "[ -d %s ] || mkdir -p %s", dir, dir);
   // gSystem->Exec(CMD);

   // const int ncent = 5;
   // const char nameCent[ncent][250] = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%"};
   // const char nameCentXL[ncent][250] = {"0_10", "10_20", "20_40", "40_60", "60_80"};
   // const float scale[ncent] = {1., 1., 1., 1., 1.};
   // float Nbin[ncent] = {938.80170, 579.89409, 288.35051, 91.37100, 21.37396};
   // float NbinErr[ncent] = {26.28048, 28.86320, 30.39279, 21.05409, 8.93878};
   const int ncent = 3;
   const char nameCent[ncent][250] = {"0-10%", "10-40%", "40-80%"};
   const char nameCentXL[ncent][250] = {"0_10", "10_40", "40_80"};
   const int scale[ncent] = {20, 5, 2};
   const float scale2[ncent] = {1., 1., 4 };
   float Nbin[ncent] = {938.80170, 386.08527, 56.99229};
   float NbinErr[ncent] = {26.28048, 29.86811, 14.85246 };


   //Read spectra
   //1. from xiaolong
   TGraphErrors* gD0err_xl[ncent];
   TGraphErrors* gD0sys_xl[ncent];
   TF1* fLevy[ncent];
   // TFile* fin1 = new TFile("D0RAA_Run14HFT.root");
   TFile* fin1 = new TFile("D0_Spectra_Run14HFT.root");
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_err_%s", nameCentXL[icent]));
      gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_sys_%s", nameCentXL[icent]));
      ScaleGraph(gD0err_xl[icent], scale[icent]);
      ScaleGraph(gD0sys_xl[icent], scale[icent]);
   }
   fin1->Close();


   //1. from longzhou
   TGraphErrors* gD0err_lz[ncent];
   TGraphErrors* gD0sys_lz[ncent];
   // TFile* fin2 = new TFile("D0_RAA_Publish_Corrected.root");
   // TFile* fin2 = new TFile("D0_Spectra_Publish_Corrected.root");
   TFile* fin2 = new TFile("all_yifei.root");
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_lz[icent] = (TGraphErrors*)fin2->Get(Form("cen_%s_err", nameCentXL[icent]));
      gD0sys_lz[icent] = (TGraphErrors*)fin2->Get(Form("cen_%s_sys", nameCentXL[icent]));

     ScaleGraph(gD0err_lz[icent], scale2[icent]);
     ScaleGraph(gD0sys_lz[icent], scale2[icent]);
   }
   fin2->Close();

   gD0err_lz[1]->RemovePoint(6);
   gD0sys_lz[1]->RemovePoint(6);
   

   TF1* myLevyFcn[ncent];
   for (int icent = 0; icent < ncent; icent++)
   {
   myLevyFcn[icent]= new TF1(Form("myLevyFcn_%d",icent), LevyFcn, 0., 7, 4);
   myLevyFcn[icent]->SetParameters(5.26, 15.89, 0.31, 1.865);
   myLevyFcn[icent]->SetLineWidth(2);
   myLevyFcn[icent]->SetLineStyle(2);
   }

   //Set for Draw
   float markerSize = 2.0;
   float lineWidth = 2;
   float markerSizeScale[ncent + 1] = {0.85, 0.75, 1., 1., };
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent]->SetMarkerStyle(20);
      gD0err_xl[icent]->SetMarkerSize(1.5);
      gD0err_xl[icent]->SetLineWidth(2);
      gD0err_xl[icent]->SetMarkerColor(1);
      gD0err_xl[icent]->SetLineColor(1);

      gD0sys_xl[icent]->SetMarkerStyle(20);
      gD0sys_xl[icent]->SetMarkerSize(1.5);
      gD0sys_xl[icent]->SetLineWidth(2);
      gD0sys_xl[icent]->SetMarkerColor(1);
      gD0sys_xl[icent]->SetLineColor(1);

      gD0err_lz[icent]->SetMarkerStyle(24);
      gD0err_lz[icent]->SetMarkerSize(1.3);
      gD0err_lz[icent]->SetLineWidth(2);
      gD0err_lz[icent]->SetMarkerColor(4);
      gD0err_lz[icent]->SetLineColor(4);

      gD0sys_lz[icent]->SetMarkerStyle(24);
      gD0sys_lz[icent]->SetMarkerSize(1.3);
      gD0sys_lz[icent]->SetLineWidth(2);
      gD0sys_lz[icent]->SetMarkerColor(4);
      gD0sys_lz[icent]->SetLineColor(4);
   }

   TGraphErrors* gD0err_xl_Ratio[ncent];
   TGraphErrors* gD0sys_xl_Ratio[ncent];
   TGraphErrors* gD0err_lz_Ratio[ncent];
   TGraphErrors* gD0sys_lz_Ratio[ncent];
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl_Ratio[icent] = (TGraphErrors*)gD0err_xl[icent]->Clone(Form("gD0_err_Ratio_%s", nameCentXL[icent]));
      gD0sys_xl_Ratio[icent] = (TGraphErrors*)gD0sys_xl[icent]->Clone(Form("gD0_sys_Ratio_%s", nameCentXL[icent]));
     gD0err_lz_Ratio[icent] = (TGraphErrors*)gD0err_lz[icent]->Clone(Form("gD0_err_Ratio_%s", nameCentXL[icent]));
     gD0sys_lz_Ratio[icent] = (TGraphErrors*)gD0sys_lz[icent]->Clone(Form("gD0_sys_Ratio_%s", nameCentXL[icent]));
   }

   gStyle->UseCurrentStyle();

   ifstream inData;
   ofstream outData;

   //plot
   gStyle->Reset("plain");
   // TCanvas* c1 = new TCanvas("c1", "A Canvas", 10, 10, 600, 800);
   TCanvas* c1 = new TCanvas("c1", "A Canvas", 10, 10, 600, 900);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   gStyle->SetTitle(0);
   // setPad(c1);
   c1->SetGridx(0);
   c1->SetGridy(0);
   // setpad(c1, 0.14, 0.02, 0.01, 0.12);
   setpad(c1, 0., 0.02, 0.01, 0.12);

   double x1 = 0;
   double x2 = 8.6;
   double y1 = 7e-7;
   double y2 = 5e1;

   float small = 0;
   c1->Divide(1, 4, small, small);

   c1->cd(1)->SetLogy(1);
   // gPad->SetPad(0., 0.68, 1., 0.98);
   // gPad->SetPad(0., 0.5, 1., 0.98);
   gPad->SetPad(0., 0.55, 1., 1.);
   setpad(gPad, 0.16, 0.05, 0.05, 0.02);
   gPad->SetTickx();
   gPad->SetTicky(0);

   c1->cd(1)->cd();
   TH1* h00 = new TH1F("", "", 1, x1, x2);
   // setHisto(h00, "", "p_{T} (GeV/c)", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) [(GeV/c)^{-2}]" );
   setHisto(h00, "", "", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) [(GeV/c)^{-2}]" );
   h00->GetYaxis()->SetLabelFont(42);
   h00->GetYaxis()->SetTitleFont(42);
   h00->GetYaxis()->SetNdivisions(503);
   h00->Draw("e");
   h00->SetMinimum(y1);
   h00->SetMaximum(y2);
   h00->GetYaxis()->CenterTitle();
   // cout << h00->GetYaxis()->GetTitleOffset() << endl;
   h00->GetYaxis()->SetTitleSize(0.065);
   h00->GetYaxis()->SetTitleOffset(1.0);
   h00->GetYaxis()->SetLabelSize(0.06);
   h00->GetYaxis()->SetLabelOffset(0.01);
   h00->GetXaxis()->SetLabelSize(0);
   // h00->GetListOfFunctions()->FindObject("stats")->Delete();

   const float sysw = 0.08;
   // legend1 = new TLegend(0.81, 0.25, 0.9, 0.80);
   legend1 = new TLegend(0.73, 0.45, 0.93, 0.7);
   legend1->SetFillStyle(0);
   legend1->SetFillColor(10);
   legend1->SetBorderSize(0);
   legend1->SetTextSize(0.06);
   legend1->SetTextFont(42);
   legend1->AddEntry(gD0err_xl[0], "2014", "p");
   legend1->AddEntry(gD0err_lz[0], "2010/11", "p");
   legend1->Draw("same");

   for (int icent = 0; icent < ncent ; icent++)
   {
      // if (icent != 0) continue;
      gD0err_xl[icent]->Draw("psame");
      gD0err_lz[icent]->Draw("psame");
      //draw systematic error
      gD0err_xl[icent]->Fit(myLevyFcn[icent],"R");

      const float sysw = 0.15;
      for (int i = 0; i < gD0sys_xl[icent]->GetN()-1; i++)
      {
         const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(1);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(1);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
      }

      for (int i = 0; i < gD0sys_lz[icent]->GetN(); i++)
      {
         const float sysl = gD0sys_lz[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(4);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
      }


      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.78, 0.75, buf, 62, 0.075, 1);
   TLatex *tex = new TLatex(0.4, 3e-6, Form("(a)"));
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   // tex->SetNDC(kTRUE);
   tex->Draw("same");

   }

   legend1 = new TLegend(0.3, 0.1, 0.4, 0.2);
   legend1->SetFillStyle(0);
   legend1->SetFillColor(10);
   legend1->SetBorderSize(0);
   legend1->SetTextSize(0.055);
   legend1->SetTextFont(42);
   legend1->AddEntry(myLevyFcn[0], "Levy", "l");
   legend1->Draw("same");

   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   drawLine(x1, y2, x2, y2, 3, 1);

   sprintf(buf, "Au+Au #sqrt{s_{NN}} = 200 GeV");
   drawLatex(0.52, 0.80, buf, 42, 0.07, 1);

   // TLine *l1 = new TLine(x1, 1.3, x2, 1.3);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   // l1->Draw("same");

   for(int icent=0;icent<ncent;icent++) {
     TLatex *tex = new TLatex(6.9, gD0sys_xl[icent]->GetY()[9]*0.7, Form("%s (#times%d)",nameCent[icent],scale[icent]));
     tex->SetTextFont(42);
     tex->SetTextSize(0.05);
     tex->Draw("same");
   }
   

   // TGraphErrors* gD0err_xl_Ratio[ncent];
   // TGraphErrors* gD0sys_xl_Ratio[ncent];
   // TGraphErrors* gD0err_lz_Ratio[ncent];
   // TGraphErrors* gD0sys_lz_Ratio[ncent];

   for (int icent = 0; icent < ncent ; icent++)
   {
     for (int i=0;i<gD0err_xl_Ratio[icent]->GetN();i++)
     {
     gD0err_xl_Ratio[icent]->GetY()[i] /= myLevyFcn[icent]->Eval(gD0err_xl_Ratio[icent]->GetX()[i]);//bin shift
     gD0err_xl_Ratio[icent]->GetEY()[i] /= myLevyFcn[icent]->Eval(gD0err_xl_Ratio[icent]->GetX()[i]);//bin shift
     gD0sys_xl_Ratio[icent]->GetY()[i] /= myLevyFcn[icent]->Eval(gD0sys_xl_Ratio[icent]->GetX()[i]);//bin shift
     gD0sys_xl_Ratio[icent]->GetEY()[i] /= myLevyFcn[icent]->Eval(gD0sys_xl_Ratio[icent]->GetX()[i]);//bin shift
     }

     for (int i=0;i<gD0err_lz_Ratio[icent]->GetN();i++)
     {
     gD0err_lz_Ratio[icent]->GetY()[i] /= myLevyFcn[icent]->Eval(gD0err_lz_Ratio[icent]->GetX()[i]);//bin shift
     gD0err_lz_Ratio[icent]->GetEY()[i] /= myLevyFcn[icent]->Eval(gD0err_lz_Ratio[icent]->GetX()[i]);//bin shift
     gD0sys_lz_Ratio[icent]->GetY()[i] /= myLevyFcn[icent]->Eval(gD0sys_lz_Ratio[icent]->GetX()[i]);//bin shift
     gD0sys_lz_Ratio[icent]->GetEY()[i] /= myLevyFcn[icent]->Eval(gD0sys_lz_Ratio[icent]->GetX()[i]);//bin shift
     }
   }


   c1->cd(2)->SetLogy(0);
   c1->cd(2);
   // gPad->SetPad(0., 0.48, 1., 0.68);
   gPad->SetPad(0., 0.39, 1., 0.55);
   setpad(gPad, 0.16, 0.05, 0.01, 0.0);
   gPad->SetTickx();
   gPad->SetTicky(0);

   TH1F* h01 = new TH1F("", "", 1, x1, x2);
   // setHisto(h01, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   // setHisto(h01, "", "p_{T} (GeV/c)", "Ratio");
   setHisto(h01, "", "p_{T} (GeV/c)", "");
   h01->GetYaxis()->SetLabelFont(42);
   h01->GetYaxis()->SetTitleFont(42);
   // y1 = 0.38;
   // y2 = 1.72;
   y1 = 0.1;
   y2 = 1.9;
   h01->Draw();
   h01->SetMinimum(y1);
   h01->SetMaximum(y2);
   h01->GetYaxis()->CenterTitle();
   h01->GetYaxis()->SetTitleSize(0.14);
   h01->GetYaxis()->SetTitleOffset(0.5);
   h01->GetYaxis()->SetLabelSize(0.15);
   h01->GetYaxis()->SetNdivisions(505);
   h01->GetXaxis()->SetTickSize(0.06);
   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 0) continue;
      gD0err_xl_Ratio[icent]->Draw("psame");
      gD0err_lz_Ratio[icent]->Draw("psame");
      //draw systematic error

      const float sysw = 0.15;
      for (int i = 0; i < gD0sys_xl_Ratio[icent]->GetN()-1; i++)
      {
         const float sysl = gD0sys_xl_Ratio[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(1);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(1);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
      }

      for (int i = 0; i < gD0sys_lz_Ratio[icent]->GetN(); i++)
      {
         const float sysl = gD0sys_lz_Ratio[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(4);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
      }

   sprintf(buf, Form("%s%%", nameCent[icent]));
   drawLatex(0.83, 0.2, buf, 42, 0.15, 1);

   TLatex *tex = new TLatex(0.4, 0.2, Form("(b) "));
   tex->SetTextFont(42);
   tex->SetTextSize(0.15);
   // tex->SetNDC(kTRUE);
   tex->Draw("same");

   }
   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");

   c1->cd(3)->SetLogy(0);
   // gPad->SetPad(0., 0.28, 1., 0.48);
   gPad->SetPad(0., 0.23, 1., 0.39);
   setpad(gPad, 0.16, 0.05, 0.0, 0.0);
   gPad->SetTickx();
   gPad->SetTicky(0);

   TH1F* h02 = new TH1F("", "", 1, x1, x2);
   // setHisto(h02, "", "p_{T} (GeV/c)", "Ratio");
   setHisto(h02, "", "p_{T} (GeV/c)", "");
   h02->GetYaxis()->SetLabelFont(42);
   h02->GetYaxis()->SetTitleFont(42);
   // y1 = 0.38;
   // y2 = 1.72;
   y1 = 0.1;
   y2 = 1.9;
   h02->Draw();
   h02->SetMinimum(y1);
   h02->SetMaximum(y2);
   h02->GetYaxis()->CenterTitle();
   h02->GetYaxis()->SetTitleSize(0.13);
   h02->GetYaxis()->SetTitleOffset(0.5);
   h02->GetYaxis()->SetLabelSize(0.15);
   h02->GetXaxis()->SetTickSize(0.06);
   h02->GetYaxis()->SetNdivisions(505);
   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 1) continue;
      gD0err_xl_Ratio[icent]->Draw("psame");
      gD0err_lz_Ratio[icent]->Draw("psame");
      //draw systematic error

      const float sysw = 0.15;
      for (int i = 0; i < gD0sys_xl_Ratio[icent]->GetN()-1; i++)
      {
         const float sysl = gD0sys_xl_Ratio[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(1);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(1);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
      }

      for (int i = 0; i < gD0sys_lz_Ratio[icent]->GetN(); i++)
      {
         const float sysl = gD0sys_lz_Ratio[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(4);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
      }

      sprintf(buf, Form("%s%%", nameCent[icent]));
      drawLatex(0.83, 0.2, buf, 42, 0.145, 1);
   TLatex *tex = new TLatex(0.4, 0.2, Form("(c)", nameCent[icent]));
   tex->SetTextFont(42);
   tex->SetTextSize(0.15);
   // tex->SetNDC(kTRUE);
   tex->Draw("same");

   }
   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");


   c1->cd(4)->SetLogy(0);
   // gPad->SetPad(0., 0.28, 1., 0.48);
   gPad->SetPad(0., 0, 1., 0.2304);
   setpad(gPad, 0.16, 0.05, 0.0, 0.3);
   gPad->SetTickx();
   gPad->SetTicky(0);

   TH1F* h03 = new TH1F("", "", 1, x1, x2);
   // setHisto(h03, "", "p_{T} (GeV/c)", "Ratio");
   setHisto(h03, "", "p_{T} (GeV/c)", "");
   h03->GetYaxis()->SetLabelFont(42);
   h03->GetYaxis()->SetTitleFont(42);
   // y1 = 0.38;
   // y2 = 1.72;
   y1 = 0.1;
   y2 = 1.9;
   h03->Draw();
   h03->SetMinimum(y1);
   h03->SetMaximum(y2);
   h03->GetYaxis()->CenterTitle();
   h03->GetYaxis()->SetTitleSize(0.145);
   h03->GetYaxis()->SetTitleOffset(0.5);
   h03->GetYaxis()->SetLabelSize(0.115);
   h03->GetXaxis()->SetTickSize(0.06);
   h03->GetYaxis()->SetNdivisions(505);
   h03->GetXaxis()->SetLabelSize(0.13);
   h03->GetXaxis()->SetTitleSize(0.12);
   h03->GetXaxis()->CenterTitle();
   h03->GetXaxis()->SetTitleOffset(1.1);
   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 2) continue;
      gD0err_xl_Ratio[icent]->Draw("psame");
      gD0err_lz_Ratio[icent]->Draw("psame");
      //draw systematic error

      const float sysw = 0.15;
      for (int i = 0; i < gD0sys_xl_Ratio[icent]->GetN()-1; i++)
      {
         const float sysl = gD0sys_xl_Ratio[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(1);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(1);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] - gD0sys_xl_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] - sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i], gD0sys_xl_Ratio[icent]->GetX()[i] + sysw, gD0sys_xl_Ratio[icent]->GetY()[i] + gD0sys_xl_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
      }

      for (int i = 0; i < gD0sys_lz_Ratio[icent]->GetN(); i++)
      {
         const float sysl = gD0sys_lz_Ratio[icent]->GetY()[i] * 0.05;
         TLine *llw = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(4);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] - gD0sys_lz_Ratio[icent]->GetEY()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] - sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i], gD0sys_lz_Ratio[icent]->GetX()[i] + sysw, gD0sys_lz_Ratio[icent]->GetY()[i] + gD0sys_lz_Ratio[icent]->GetEY()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
      }

      sprintf(buf, Form("%s%%", nameCent[icent]));
      drawLatex(0.83, 0.35, buf, 42, 0.105, 1);
   TLatex *tex = new TLatex(0.4, 0.2, Form("(d)", nameCent[icent]));
   tex->SetTextFont(42);
   tex->SetTextSize(0.11);
   // tex->SetNDC(kTRUE);
   tex->Draw("same");

   }
   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");

   c1->cd();
   mpad = new TPad(Form("mpad"), "", 0.09, 0.03, 0.0, 0.5);
   // setpad(mpad, -0.5, 0., 0., 0.);
   // mpad->SetFillStyle(4000);
   mpad->SetFillStyle(1001);
   mpad->SetFillColor(10);
   mpad->SetBorderMode(0);
   mpad->SetBorderSize(0);
   mpad->SetFrameFillColor(10);
   mpad->SetFrameBorderMode(0);
   mpad->SetFrameBorderSize(0);
   mpad->SetFrameLineWidth(1);
   mpad->SetTickx();
   mpad->SetTicky(0);
   mpad->Draw();
   mpad->cd();
   sprintf(buf, Form("Ratio"));
   TLatex *tex = new TLatex(0.8, 0.5, buf);
   tex->SetTextFont(42);
   tex->SetTextSize(0.47);
   tex->SetTextAngle(90);
   tex->SetNDC();
   tex->Draw("same");

   c1->cd();
   c1->Update();

   c1->SaveAs("D0_compareSpectra_run10.eps");
   c1->SaveAs("D0_compareSpectra_run10.pdf");
   c1->SaveAs("D0_compareSpectra_run10.gif");
   c1->SaveAs("D0_compareSpectra_run10.png");
}
