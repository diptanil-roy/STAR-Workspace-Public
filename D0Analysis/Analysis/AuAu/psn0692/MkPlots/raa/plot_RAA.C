#include "myFunction.h"
#include "myConst.h"
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

void plot_RAA()
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
   const float scale[ncent] = {1., 1., 1.};
   float Nbin[ncent] = {938.80170, 386.08527, 56.99229};
   float NbinErr[ncent] = {26.28048, 29.86811, 14.85246 };


   //Read spectra
   //1. from xiaolong
   TGraphErrors* gD0err_xl[ncent];
   TGraphErrors* gD0err_xl_cp1[ncent];
   TGraphErrors* gD0err_xl_cp2[ncent];
   TGraphErrors* gD0sys_xl[ncent];
   TGraphAsymmErrors* gD0sys_xl_combine[ncent];
   TGraphAsymmErrors* gD0_RAA_pperr[ncent];
   TF1* fLevy[ncent];
   TFile* fin1 = new TFile("D0RAA_Run14HFT.root");
   // TFile* fin1 = new TFile("D0RAA_Run14HFT.root.0430");
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("D0_RAA_err_%s", nameCentXL[icent]));
      gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("D0_RAA_sys_%s", nameCentXL[icent]));
      gD0_RAA_pperr[icent] = (TGraphAsymmErrors*)fin1->Get(Form("D0_RAA_pperr_%s", nameCentXL[icent]));
      gD0sys_xl_combine[icent] = (TGraphAsymmErrors*)fin1->Get(Form("D0_RAA_sys_combine_%s", nameCentXL[icent]));
   }
   fin1->Close();

   //1. from longzhou
   TGraphErrors* gD0err_lz[ncent];
   TGraphErrors* gD0sys_lz[ncent];
   TGraphAsymmErrors* gD0sys_lz2[ncent];
   // TFile* fin2 = new TFile("D0_RAA_Publish_Corrected.root");
   // TFile* fin2 = new TFile("all_yifei.root");
   // TFile* fin2 = new TFile("D0_RAA_Publish_Corrected_0115.root");//with new sys
   // TFile* fin2 = new TFile("D0_RAA_AsymError.root");//with new sys
   TFile* fin2 = new TFile("D0_RAA_AsymError1.root");//with new sys
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_lz[icent] = (TGraphErrors*)fin2->Get(Form("RAA_%sc_err", nameCentXL[icent]));
      // gD0sys_lz[icent] = (TGraphErrors*)fin2->Get(Form("RAA_%sc_sys", nameCentXL[icent]));
      gD0sys_lz2[icent] = (TGraphAsymmErrors*)fin2->Get(Form("RAA_%sc_sys", nameCentXL[icent]));
   }
   fin2->Close();


   //Set for Draw
   float markerSize = 2.0;
   float lineWidth = 2;
   float markerSizeScale[ncent + 1] = {0.85, 0.75, 1., 1., };
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent]->SetMarkerStyle(20);
      // gD0err_xl[icent]->SetMarkerStyle(27);
      gD0err_xl[icent]->SetMarkerSize(1.5);
      gD0err_xl[icent]->SetLineWidth(2);
      gD0err_xl[icent]->SetMarkerColor(1);
      gD0err_xl[icent]->SetLineColor(1);

      // gD0sys_xl[icent]->SetMarkerStyle(20);
      // gD0sys_xl[icent]->SetMarkerSize(1.5);
      // gD0sys_xl[icent]->SetLineWidth(2);
      // gD0sys_xl[icent]->SetMarkerColor(1);
      // gD0sys_xl[icent]->SetLineColor(1);

      gD0sys_xl_combine[icent]->SetMarkerStyle(20);
      gD0sys_xl_combine[icent]->SetMarkerSize(1.5);
      gD0sys_xl_combine[icent]->SetLineWidth(2);
      gD0sys_xl_combine[icent]->SetMarkerColor(1);
      gD0sys_xl_combine[icent]->SetLineColor(1);

      // gD0err_lz[icent]->SetMarkerStyle(24);
      // gD0err_lz[icent]->SetMarkerSize(1.3);
      gD0err_lz[icent]->SetMarkerStyle(27);
      gD0err_lz[icent]->SetMarkerSize(1.7);
      gD0err_lz[icent]->SetLineWidth(2);
      gD0err_lz[icent]->SetMarkerColor(4);
      gD0err_lz[icent]->SetLineColor(4);

      // gD0sys_lz[icent]->SetMarkerStyle(24);
      // gD0sys_lz[icent]->SetMarkerSize(1.3);
      // gD0sys_lz[icent]->SetLineWidth(2);
      // gD0sys_lz[icent]->SetMarkerColor(4);
      // gD0sys_lz[icent]->SetLineColor(4);

      gD0sys_lz2[icent]->SetMarkerStyle(24);
      gD0sys_lz2[icent]->SetMarkerSize(1.3);
      gD0sys_lz2[icent]->SetLineWidth(2);
      gD0sys_lz2[icent]->SetMarkerColor(4);
      gD0sys_lz2[icent]->SetLineColor(4);
   }

   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl_cp1[icent] = (TGraphErrors*)gD0err_xl[icent]->Clone(Form("D0_RAA_err_%s_cp1", nameCentXL[icent]));
      gD0err_xl_cp2[icent] = (TGraphErrors*)gD0err_xl[icent]->Clone(Form("D0_RAA_err_%s_cp2", nameCentXL[icent]));

      gD0err_xl_cp1[icent]->RemovePoint(10);
      gD0err_xl_cp1[icent]->RemovePoint(9);
      gD0err_xl_cp1[icent]->RemovePoint(1);
      gD0err_xl_cp1[icent]->RemovePoint(0);

      gD0err_xl_cp2[icent]->SetMarkerStyle(24);
      gD0err_xl_cp2[icent]->RemovePoint(8);
      gD0err_xl_cp2[icent]->RemovePoint(7);
      gD0err_xl_cp2[icent]->RemovePoint(6);
      gD0err_xl_cp2[icent]->RemovePoint(5);
      gD0err_xl_cp2[icent]->RemovePoint(4);
      gD0err_xl_cp2[icent]->RemovePoint(3);
      gD0err_xl_cp2[icent]->RemovePoint(2);

   }

   gStyle->UseCurrentStyle();

   ifstream inData;
   ofstream outData;

   //============SHANSHAN LBT D-meson RAA 0_5% ============
   ifstream inLBT("./Models/LBT/RAA_HMnoHI-new.dat");
   float LBT_pt_cen[20],LBT_RAA_cen[20];  
   for(int i=0; i<20; i++) { inLBT >> LBT_pt_cen[i] >> LBT_RAA_cen[i]; }

   //============DUKE LGV D-meson RAA  ============
   ifstream inLGV("./Models/DUKELGV/D_meson_spectra_AuAu200_forSTAR_DukeLGV_reNorm.dat");
   inLGV.getline(buf, 1024);
   inLGV.getline(buf, 1024);
   inLGV.getline(buf, 1024);
   float LGV_pt_cen[20],LGV_RAA_cen[20], LGV_pp[20], LGV_tmp[20];  
   // # pT[GeV/c] pp(charmH) pp(D-meson) 0-10%      10-20%     20-40%     40-60%
   for(int i=0; i<20; i++) { 
     inLGV >> LGV_pt_cen[i] >> LGV_tmp[i] >> LGV_pp[i] >> LGV_RAA_cen[i]>> LGV_tmp[i] >> LGV_tmp[i] >> LGV_tmp[i] ;
     LGV_RAA_cen[i] /= LGV_pp[i];
   }

   //plot
   gStyle->Reset("plain");
   // TCanvas* c1 = new TCanvas("c1", "A Canvas", 10, 10, 600, 800);
   TCanvas* c1 = new TCanvas("c1", "A Canvas", 10, 10, 600, 600);
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
   double x2 = 9.6;
   double y1 = 0.02;
   double y2 = 1.65;
   // double y2 = 2.05;

   float small = 0;
   c1->Divide(1, 3, small, small);

   c1->cd(1)->SetLogy(0);
   // gPad->SetPad(0., 0.68, 1., 0.98);
   gPad->SetPad(0., 0.63, 1., 0.98);
   setpad(gPad, 0.16, 0.05, 0.1, 0.01);
   gPad->SetTickx();
   gPad->SetTicky(0);

   c1->cd(1)->cd();
   TH1* h00 = new TH1F("", "", 1, x1, x2);
   // setHisto(h00, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   setHisto(h00, "", "p_{T} (GeV/c)", "" );
   h00->GetYaxis()->SetLabelFont(42);
   h00->GetYaxis()->SetTitleFont(42);
   h00->Draw("e");
   h00->SetMinimum(y1);
   h00->SetMaximum(y2);
   h00->GetYaxis()->CenterTitle();
   // cout << h00->GetYaxis()->GetTitleOffset() << endl;
   h00->GetYaxis()->SetTitleSize(0.1);
   h00->GetYaxis()->SetTitleOffset(0.65);
   h00->GetYaxis()->SetLabelSize(0.1);
   // h00->GetListOfFunctions()->FindObject("stats")->Delete();

   const float sysw = 0.08;
   // legend1 = new TLegend(0.81, 0.25, 0.9, 0.80);
   legend1 = new TLegend(0.78, 0.32, 0.93, 0.8);
   legend1->SetFillStyle(0);
   legend1->SetFillColor(10);
   legend1->SetBorderSize(0);
   legend1->SetTextSize(0.08);
   legend1->SetTextFont(42);
   legend1->AddEntry(gD0err_xl[0], "Run14", "p");
   legend1->AddEntry(gD0err_lz[0], "Run10/11", "p");
   legend1->Draw("same");

   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 0) continue;
      
      const float sysw = 0.15;
      // for (int i = 0; i < gD0_RAA_pperr[icent]->GetN(); i++)
      // {
      //   TBox *bx = new TBox(gD0_RAA_pperr[icent]->GetX()[i]-sysw,gD0_RAA_pperr[icent]->GetY()[i]- gD0_RAA_pperr[icent]->GetEYlow()[i],gD0_RAA_pperr[icent]->GetX()[i]+sysw,gD0_RAA_pperr[icent]->GetY()[i]+gD0_RAA_pperr[icent]->GetEYhigh()[i]);
      //   bx->SetLineColor(18);
      //   bx->SetFillColor(18);
      //   bx->SetLineWidth(1);
      //   bx->Draw("same");
      // }

      // gD0err_xl[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");
      gD0err_lz[icent]->Draw("psame");
      //draw systematic error

      const float sysw = 0.15;
      // for (int i = 0; i < gD0sys_xl[icent]->GetN(); i++)
      // {
      //    const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(1);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(1);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      // }

      for (int i = 0; i < gD0sys_xl_combine[icent]->GetN(); i++)
      {
         // const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
         const float sysl = 0.06;
         TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(1);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(1);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
      }


      // for (int i = 0; i < gD0sys_lz[icent]->GetN(); i++)
      // {
      //    const float sysl = gD0sys_lz[icent]->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(4);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(4);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      // }

      for (int i = 0; i < gD0sys_lz2[icent]->GetN(); i++)
      {
         // const float sysl = gD0sys_lz2[icent]->GetY()[i] * 0.05;
         const float sysl = 0.06;
         TLine *llw = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(4);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
      }

      // gD0err_xl[icent]->Draw("psame");
      gD0err_lz[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");

      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.78, 0.75, buf, 62, 0.075, 1);
   TLatex *tex = new TLatex(7.8, 0.43, Form("(a)  %s%", nameCent[icent]));
   tex->SetTextFont(42);
   tex->SetTextSize(0.09);
   // tex->SetNDC(kTRUE);
   tex->Draw("same");

   }
   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   drawLine(x1, y2, x2, y2, 3, 1);

   sprintf(buf, "Au+Au #sqrt{s_{NN}} = 200 GeV");
   drawLatex(0.22, 0.76, buf, 42, 0.11, 1);

   //LBT_pt_cen RAA 0-10%
   TGraph *gLBTraa = new TGraph(20,LBT_pt_cen,LBT_RAA_cen); 
   gLBTraa->SetLineColor(4);
   gLBTraa->SetLineStyle(4);
   gLBTraa->SetLineWidth(2);
   // gLBTraa->Draw("c");

   //LGV_pt_cen RAA 0-10%
   TGraph *gLGVraa = new TGraph(20,LGV_pt_cen,LGV_RAA_cen); 
   gLGVraa->SetLineColor(kAzure+10);
   gLGVraa->SetLineStyle(10);
   gLGVraa->SetLineWidth(2);
   // gLGVraa->Draw("c");

   legend2 = new TLegend(0.64, 0.32, 0.73, 0.8);
   legend2->SetFillStyle(0);
   legend2->SetFillColor(10);
   legend2->SetBorderSize(0);
   legend2->SetTextSize(0.08);
   legend2->SetTextFont(42);
   legend2->AddEntry(gLBTraa,"LBT","l");
   legend2->AddEntry(gLGVraa,"Duke","l");
   // legend2->Draw("same");

   // TLine *l1 = new TLine(x1, 1.3, x2, 1.3);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");

   TBox *bx = new TBox(9.3,(1.-NbinErr[0]/Nbin[0]),9.45,(1.+NbinErr[0]/Nbin[0]));
   bx->SetFillColor(kGreen-6);
   bx->Draw();
   const float ppnormsys = sqrt(0.081*0.081+0.052*0.052); //pp normalization 0.081%
   TBox *bx = new TBox(9.45,(1.-ppnormsys),9.58,(1.+ppnormsys));
   bx->SetFillColor(kGreen-2);
   bx->Draw();

   c1->cd(2)->SetLogy(0);
   c1->cd(2);
   // gPad->SetPad(0., 0.48, 1., 0.68);
   gPad->SetPad(0., 0.36, 1., 0.63);
   setpad(gPad, 0.16, 0.05, 0.03, 0.01);
   gPad->SetTickx();
   gPad->SetTicky(0);

   TH1F* h0 = new TH1F("", "", 1, x1, x2);
   // setHisto(h0, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   setHisto(h0, "", "p_{T} (GeV/c)", "");
   h0->GetYaxis()->SetLabelFont(42);
   h0->GetYaxis()->SetTitleFont(42);
   // y1 = 0.3;
   // y2 = 1.44;
   // y2 = 2.05;
   y2 = 1.65;
   h0->Draw();
   h0->SetMinimum(y1);
   h0->SetMaximum(y2);
   h0->GetYaxis()->CenterTitle();
   h0->GetYaxis()->SetTitleSize(0.13);
   h0->GetYaxis()->SetTitleOffset(0.5);
   h0->GetYaxis()->SetLabelSize(0.13);
   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 1) continue;

      const float sysw = 0.15;
      // for (int i = 0; i < gD0_RAA_pperr[icent]->GetN(); i++)
      // {
      //   TBox *bx = new TBox(gD0_RAA_pperr[icent]->GetX()[i]-sysw,gD0_RAA_pperr[icent]->GetY()[i]- gD0_RAA_pperr[icent]->GetEYlow()[i],gD0_RAA_pperr[icent]->GetX()[i]+sysw,gD0_RAA_pperr[icent]->GetY()[i]+gD0_RAA_pperr[icent]->GetEYhigh()[i]);
      //   bx->SetLineColor(18);
      //   bx->SetFillColor(18);
      //   bx->SetLineWidth(1);
      //   bx->Draw("same");
      // }
      // gD0err_xl[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");
      gD0err_lz[icent]->Draw("psame");
      //draw systematic error

      const float sysw = 0.15;
      // for (int i = 0; i < gD0sys_xl[icent]->GetN(); i++)
      // {
      //    const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(1);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(1);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      // }

      for (int i = 0; i < gD0sys_xl_combine[icent]->GetN(); i++)
      {
         // const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
         const float sysl = 0.06;
         TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(1);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(1);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
      }

      // for (int i = 0; i < gD0sys_lz[icent]->GetN(); i++)
      // {
      //    const float sysl = gD0sys_lz[icent]->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(4);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(4);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      // }
      
      for (int i = 0; i < gD0sys_lz2[icent]->GetN() -1; i++)
      {
         // const float sysl = gD0sys_lz2[icent]->GetY()[i] * 0.05;
         const float sysl = 0.06;
         TLine *llw = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(4);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
      }

        TArrow *arw1 = new TArrow(gD0sys_lz2[icent]->GetX()[6], gD0sys_lz2[icent]->GetY()[6] + gD0sys_lz2[icent]->GetEYhigh()[6], gD0sys_lz2[icent]->GetX()[6], gD0sys_lz2[icent]->GetY()[6] + gD0sys_lz2[icent]->GetEYhigh()[6] - 0.7, 0.02, ">");
         arw1->SetLineWidth(2);
         arw1->SetLineColor(4);
         arw1->Draw(">");
         TLine *lhi = new TLine(gD0sys_lz2[icent]->GetX()[6] - sysw, gD0sys_lz2[icent]->GetY()[6] + gD0sys_lz2[icent]->GetEYhigh()[6], gD0sys_lz2[icent]->GetX()[6] + sysw, gD0sys_lz2[icent]->GetY()[6] + gD0sys_lz2[icent]->GetEYhigh()[6]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");

      // gD0err_xl[icent]->Draw("psame");
      gD0err_lz[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");

      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.83, 0.38, buf, 42, 0.12, 1);
   TLatex *tex = new TLatex(7.8, 1.2, Form("(b)  %s%", nameCent[icent]));
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

   TBox *bx = new TBox(9.3,(1.-NbinErr[1]/Nbin[1]),9.45,(1.+NbinErr[1]/Nbin[1]));
   bx->SetFillColor(kGreen-6);
   bx->Draw();
   const float ppnormsys = sqrt(0.081*0.081+0.052*0.052); //pp normalization 0.081%
   TBox *bx = new TBox(9.45,(1.-ppnormsys),9.58,(1.+ppnormsys));
   bx->SetFillColor(kGreen-2);
   bx->Draw();

   c1->cd(3)->SetLogy(0);
   // gPad->SetPad(0., 0.28, 1., 0.48);
   gPad->SetPad(0., 0.0, 1., 0.36);
   setpad(gPad, 0.16, 0.05, 0.03, 0.23);
   gPad->SetTickx();
   gPad->SetTicky(0);

   // y2 = 2.05;
   y2 = 1.65;
   TH1F* h02 = new TH1F("", "", 1, x1, x2);
   // setHisto(h02, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   setHisto(h02, "", "p_{T} (GeV/c)", "");
   h02->GetYaxis()->SetLabelFont(42);
   h02->GetYaxis()->SetTitleFont(42);
   h02->Draw();
   h02->SetMinimum(y1);
   h02->SetMaximum(y2);
   h02->GetYaxis()->CenterTitle();
   h02->GetYaxis()->SetTitleSize(0.095);
   h02->GetYaxis()->SetTitleOffset(0.68);
   h02->GetYaxis()->SetLabelSize(0.11);
   h02->GetXaxis()->SetLabelFont(42);
   h02->GetXaxis()->SetTitleFont(42);
   h02->GetXaxis()->CenterTitle();
   h02->GetXaxis()->SetTitleSize(0.1);
   h02->GetXaxis()->SetTitleOffset(1.1);
   h02->GetXaxis()->SetLabelSize(0.10);
   h02->GetXaxis()->SetLabelOffset(0.025);
   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 2) continue;

      const float sysw = 0.15;
      // for (int i = 0; i < gD0_RAA_pperr[icent]->GetN(); i++)
      // {
      //   TBox *bx = new TBox(gD0_RAA_pperr[icent]->GetX()[i]-sysw,gD0_RAA_pperr[icent]->GetY()[i]- gD0_RAA_pperr[icent]->GetEYlow()[i],gD0_RAA_pperr[icent]->GetX()[i]+sysw,gD0_RAA_pperr[icent]->GetY()[i]+gD0_RAA_pperr[icent]->GetEYhigh()[i]);
      //   bx->SetLineColor(18);
      //   bx->SetFillColor(18);
      //   bx->SetLineWidth(1);
      //   bx->Draw("same");
      // }
      // gD0err_xl[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");
      gD0err_lz[icent]->Draw("psame");

      const float sysw = 0.15;
      // for (int i = 0; i < gD0sys_xl[icent]->GetN(); i++)
      // {
      //    const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(1);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(1);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      // }

      for (int i = 0; i < gD0sys_xl_combine[icent]->GetN(); i++)
      {
         // const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
         const float sysl = 0.06;
         TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(1);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(1);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl_combine[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl_combine[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(1);
         lv->Draw("same");
      }

      // for (int i = 0; i < gD0sys_lz[icent]->GetN(); i++)
      // {
      //    const float sysl = gD0sys_lz[icent]->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(4);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(4);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] - gD0sys_lz[icent]->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] - sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i], gD0sys_lz[icent]->GetX()[i] + sysw, gD0sys_lz[icent]->GetY()[i] + gD0sys_lz[icent]->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(4);
      //    lv->Draw("same");
      // }

      for (int i = 0; i < gD0sys_lz2[icent]->GetN(); i++)
      {
         // const float sysl = gD0sys_lz2[icent]->GetY()[i] * 0.05;
         const float sysl = 0.06;
         TLine *llw = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i]);
         llw->SetLineWidth(2);
         llw->SetLineColor(4);
         llw->Draw("same");
         TLine *lhi = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i]);
         lhi->SetLineWidth(2);
         lhi->SetLineColor(4);
         lhi->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] - gD0sys_lz2[icent]->GetEYlow()[i] + sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] - sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
         TLine *lv = new TLine(gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i], gD0sys_lz2[icent]->GetX()[i] + sysw, gD0sys_lz2[icent]->GetY()[i] + gD0sys_lz2[icent]->GetEYhigh()[i] - sysl);
         lv->SetLineWidth(2);
         lv->SetLineColor(4);
         lv->Draw("same");
      }

      // gD0err_xl[icent]->Draw("psame");
      gD0err_lz[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");
      
      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.83, 0.56, buf, 42, 0.085, 1);
   TLatex *tex = new TLatex(7.8, 1.18, Form("(c)  %s%", nameCent[icent]));
   tex->SetTextFont(42);
   tex->SetTextSize(0.08);
   // tex->SetNDC(kTRUE);
   tex->Draw("same");

   }
   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   drawLine(x1, y1, x2, y1, 3, 1);

   sprintf(buf, "Au+Au #sqrt{s_{NN}} = 200 GeV");
   // drawLatex(0.18, 0.88, buf, 42, 0.045, 1);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");

   TBox *bx = new TBox(9.3,(1.-NbinErr[2]/Nbin[2]),9.45,(1.+NbinErr[2]/Nbin[2]));
   bx->SetFillColor(kGreen-6);
   bx->Draw();
   const float ppnormsys = sqrt(0.081*0.081+0.052*0.052); //pp normalization 0.081%
   TBox *bx = new TBox(9.45,(1.-ppnormsys),9.58,(1.+ppnormsys));
   bx->SetFillColor(kGreen-2);
   bx->Draw();

   c1->cd();
   mpad = new TPad(Form("mpad"), "", 0.06, 0.0, 0.0, 1.0);
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
   sprintf(buf, Form("R_{AA}", nameCent[ncent - 1]));
   TLatex *tex = new TLatex(0.8, 0.5, buf);
   tex->SetTextFont(42);
   tex->SetTextSize(0.8);
   tex->SetTextAngle(90);
   tex->SetNDC();
   tex->Draw("same");

   c1->cd();
   c1->Update();

   // c1->SaveAs("D0_RAA_model.eps");
   // c1->SaveAs("D0_RAA_model.pdf");
   // c1->SaveAs("D0_RAA_model.gif");
   // c1->SaveAs("D0_RAA_model.png");
   c1->SaveAs("D0_RAA.eps");
   c1->SaveAs("D0_RAA.pdf");
   c1->SaveAs("D0_RAA.gif");
   c1->SaveAs("D0_RAA.png");
}
