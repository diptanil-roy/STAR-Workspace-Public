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

void plot_RAA_LHC()
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
   // TFile* fin2 = new TFile("D0_RAA_Publish_Corrected.root");
   // TFile* fin2 = new TFile("all_yifei.root");
   TFile* fin2 = new TFile("D0_RAA_Publish_Corrected_0115.root");//with new sys
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_lz[icent] = (TGraphErrors*)fin2->Get(Form("RAA_%s_err", nameCentXL[icent]));
      gD0sys_lz[icent] = (TGraphErrors*)fin2->Get(Form("RAA_%s_sys", nameCentXL[icent]));
   }
   fin2->Close();


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

      gD0sys_xl_combine[icent]->SetMarkerStyle(20);
      gD0sys_xl_combine[icent]->SetMarkerSize(1.5);
      gD0sys_xl_combine[icent]->SetLineWidth(2);
      gD0sys_xl_combine[icent]->SetMarkerColor(1);
      gD0sys_xl_combine[icent]->SetLineColor(1);

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

   // RAA of D of LHC ALICE 0-10% D-average  - PRC ...
   const Int_t n_D_RAA_LHC = 10;
   double x_D_RAA_LHC[n_D_RAA_LHC] = { 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 10.0, 14.0, 20.0, 30.0};
   double y_D_RAA_LHC[n_D_RAA_LHC] = { 0.695, 0.694, 0.400, 0.250, 0.219, 0.172, 0.149, 0.180, 0.223, 0.352};
   double ye_D_RAA_LHC[n_D_RAA_LHC] = { 0.210, 0.079, 0.034, 0.023, 0.021, 0.016, 0.015, 0.030, 0.042, 0.082};
   double yesl_D_RAA_LHC[n_D_RAA_LHC] = { 0.303, 0.289, 0.118, 0.068, 0.054, 0.041, 0.035, 0.045, 0.080, 0.165};
   double yesr_D_RAA_LHC[n_D_RAA_LHC] = { 0.257, 0.199, 0.088, 0.056, 0.049, 0.037, 0.033, 0.044, 0.066, 0.120};
   TGraphErrors *gr_D_RAA_LHC = new TGraphErrors(n_D_RAA_LHC, x_D_RAA_LHC, y_D_RAA_LHC, 0, ye_D_RAA_LHC);

   // RAA of pi0 at RHIC PHEHIX 0-10%?
   const Int_t n_pi_RAA_RHIC = 23;
   double x_pi_RAA_RHIC[n_pi_RAA_RHIC] = {1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11, 13, 15, 17, 19};
   double y_pi_RAA_RHIC[n_pi_RAA_RHIC] = {0.364, 0.4123, 0.3958, 0.3159, 0.2556, 0.2272, 0.2103, 0.2051, 0.1855, 0.1854, 0.1911, 0.207, 0.1964, 0.1875, 0.1856, 0.1901, 0.2004, 0.22, 0.2099, 0.1663, 0.2892, 0.3739, 0.3449};
   double ye_pi_RAA_RHIC[n_pi_RAA_RHIC] = {0.005983, 0.005426, 0.005159, 0.004425, 0.004123, 0.004208, 0.003988, 0.004284, 0.004507, 0.005381, 0.006178, 0.005652, 0.006401, 0.007132, 0.00837, 0.009373, 0.0117, 0.01542, 0.01134, 0.02028, 0.04825, 0.09209, 0.1869};
   double yes_pi_RAA_RHIC[n_pi_RAA_RHIC] = {0.03572, 0.04197, 0.04287, 0.03558, 0.02952, 0.02641, 0.02547, 0.02483, 0.0224, 0.02232, 0.02299, 0.02482, 0.02355, 0.02251, 0.02237, 0.02302, 0.02452, 0.02719, 0.02812, 0.02704, 0.05437, 0.08636, 0.1534};
   TGraphErrors *gr_pi_RAA_RHIC = new TGraphErrors(n_pi_RAA_RHIC, x_pi_RAA_RHIC, y_pi_RAA_RHIC, 0, ye_pi_RAA_RHIC);

   // RAA of h at LHC ALICE 0-5%
   const Int_t n_h_RAA_LHC = 51;
   double x_h_RAA_LHC[n_h_RAA_LHC] = { 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725,
                                       0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 1.25, 1.35, 1.45,
                                       1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 2.5, 2.7, 2.9,
                                       3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 5.25, 5.75, 6.25,
                                       6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5,
                                       17.0, 19.0
                                     };
   double y_h_RAA_LHC[n_h_RAA_LHC] = { 0.1996, 0.2054, 0.2248, 0.2285, 0.2397, 0.2514, 0.2636, 0.2758, 0.2872,
                                       0.2992, 0.3093, 0.32, 0.3305, 0.3412, 0.3529, 0.3697, 0.3861, 0.3979, 0.4126,
                                       0.4192, 0.4229, 0.4272, 0.4344, 0.4344, 0.4268, 0.4174, 0.3983, 0.3798, 0.3565,
                                       0.3301, 0.3041, 0.2817, 0.2597, 0.2383, 0.2143, 0.1808, 0.1603, 0.1476, 0.145,
                                       0.1389, 0.151, 0.163, 0.1751, 0.1932, 0.2107, 0.2114, 0.2427, 0.2114, 0.2651,
                                       0.3086, 0.3476
                                     };
   double ye_h_RAA_LHC[n_h_RAA_LHC] = { 2.0E-4, 2.0E-4, 3.0E-4, 3.0E-4, 3.0E-4, 3.0E-4, 4.0E-4, 4.0E-4, 5.0E-4,
                                        5.0E-4, 6.0E-4, 6.0E-4, 7.0E-4, 8.0E-4, 6.0E-4, 7.0E-4, 9.0E-4, 0.001, 0.0011,
                                        0.0013, 0.0014, 0.0015, 0.0017, 0.0019, 0.0015, 0.0017, 0.0018, 0.0016, 0.0015,
                                        0.0014, 0.0013, 0.0013, 0.0013, 0.0014, 0.0012, 0.0013, 0.0014, 0.0016, 0.0019,
                                        0.0023, 0.0022, 0.0032, 0.0045, 0.0062, 0.0082, 0.0103, 0.0136, 0.0154, 0.0207,
                                        0.0203, 0.029
                                      };
   double yes_h_RAA_LHC[n_h_RAA_LHC] = { 0.0170, 0.0211, 0.0185, 0.0244, 0.0255, 0.0267, 0.0279, 0.0295, 0.0312,
                                         0.0327, 0.0344, 0.0364, 0.0382, 0.0395, 0.0408, 0.0429, 0.0451, 0.0464, 0.0476,
                                         0.0487, 0.0491, 0.0491, 0.0453, 0.0445, 0.0443, 0.0429, 0.0428, 0.0407, 0.0382,
                                         0.0353, 0.0326, 0.0301, 0.0278, 0.0255, 0.0229, 0.0194, 0.0172, 0.0159, 0.0156,
                                         0.0149, 0.0163, 0.0177, 0.0192, 0.0214, 0.0236, 0.0240, 0.0279, 0.0247, 0.0314,
                                         0.0373, 0.0431
                                       };
   TGraphErrors *gr_h_RAA_LHC = new TGraphErrors(n_h_RAA_LHC, x_h_RAA_LHC, y_h_RAA_LHC, 0, ye_h_RAA_LHC);

   // RAA of non-prompt Jpsi at LHC CMS 0-100%
   const Int_t n_Bjpsi_RAA_LHC = 4;
   double x_Bjpsi_RAA_LHC[] = {7.31, 8.96, 11.30, 16.53};
   double y_Bjpsi_RAA_LHC[] = {0.520, 0.430, 0.430, 0.339};
   double ye_Bjpsi_RAA_LHC[] = {0.121, 0.081, 0.092, 0.072};
   double yes_Bjpsi_RAA_LHC[] = {0.061, 0.050, 0.048, 0.041};
   TGraphErrors *gr_Bjpsi_RAA_LHC = new TGraphErrors(n_Bjpsi_RAA_LHC, x_Bjpsi_RAA_LHC, y_Bjpsi_RAA_LHC, 0, ye_Bjpsi_RAA_LHC);


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


   gStyle->UseCurrentStyle();

   ifstream inData;
   ofstream outData;

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
   double x2 = 10.6;
   double y1 = 0.05;
   double y2 = 1.65;

   float small = 0;
   c1->Divide(1, 2, small, small);

   c1->cd(1)->SetLogy(0);
   // gPad->SetPad(0., 0.68, 1., 0.98);
   gPad->SetPad(0., 0.5, 1., 0.98);
   setpad(gPad, 0.16, 0.05, 0.1, 0.01);
   gPad->SetTickx();
   gPad->SetTicky(0);

   c1->cd(1)->cd();
   TH1* h00 = new TH1F("", "", 1, x1, x2);
   // setHisto(h00, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   setHisto(h00, "", "p_{T} (GeV/c)", "");
   h00->GetYaxis()->SetLabelFont(42);
   h00->GetYaxis()->SetTitleFont(42);
   h00->Draw("e");
   h00->SetMinimum(y1);
   h00->SetMaximum(y2);
   h00->GetYaxis()->CenterTitle();
   // cout << h00->GetYaxis()->GetTitleOffset() << endl;
   h00->GetYaxis()->SetTitleSize(0.09);
   h00->GetYaxis()->SetTitleOffset(0.65);
   h00->GetYaxis()->SetLabelSize(0.08);
   // h00->GetListOfFunctions()->FindObject("stats")->Delete();

   const float sysw = 0.08;
   // legend1 = new TLegend(0.81, 0.25, 0.9, 0.80);
   legend1 = new TLegend(0.6, 0.4, 0.93, 0.7);
   legend1->SetFillStyle(0);
   legend1->SetFillColor(10);
   legend1->SetBorderSize(0);
   legend1->SetTextSize(0.07);
   legend1->SetTextFont(42);
   // legend1->AddEntry(gD0err_xl[0], "Run14", "p");
   // legend1->AddEntry(gD0err_lz[0], "Run10/11", "p");
   legend1->Draw("same");
   legend1->AddEntry(gD0err_xl[0],"D^{0} 0-10% STAR","p");
   legend1->AddEntry(gr_D_RAA_LHC,"D 0-10% ALICE","p");
   // legend1->AddEntry(gr_pi_RAA_RHIC,"#pi^{0} 0-10% PHENIX","p");
   // legend1->AddEntry(gr_h_RAA_LHC,"h^{#pm} 0-5% ALICE","p");

   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 0) continue;

      const float sysw = 0.15;
      // for (int i = 0; i < gD0_RAA_pperr[icent]->GetN(); i++)
      // {
      //    TBox *bx = new TBox(gD0_RAA_pperr[icent]->GetX()[i] - sysw, gD0_RAA_pperr[icent]->GetY()[i] - gD0_RAA_pperr[icent]->GetEYlow()[i], gD0_RAA_pperr[icent]->GetX()[i] + sysw, gD0_RAA_pperr[icent]->GetY()[i] + gD0_RAA_pperr[icent]->GetEYhigh()[i]);
      //    bx->SetLineColor(18);
      //    bx->SetFillColor(18);
      //    bx->SetLineWidth(1);
      //    bx->Draw("same");
      // }

      for (int i = 0; i < n_pi_RAA_RHIC; i++)
      {
         double xl = x_pi_RAA_RHIC[i] * 0.95;
         double xr = x_pi_RAA_RHIC[i] * 1.05;
         double yl = y_pi_RAA_RHIC[i] - yes_pi_RAA_RHIC[i];
         double yr = y_pi_RAA_RHIC[i] + yes_pi_RAA_RHIC[i];
         TBox *box = new TBox(xl, yl, xr, yr);
         box->SetFillColor(18);
         box->SetLineColor(18);
         // box->Draw("same");
      }
      for (int i = 0; i < n_h_RAA_LHC; i++)
      {
         double xl = x_h_RAA_LHC[i] * 0.95;
         double xr = x_h_RAA_LHC[i] * 1.05;
         double yl = y_h_RAA_LHC[i] - yes_h_RAA_LHC[i];
         double yr = y_h_RAA_LHC[i] + yes_h_RAA_LHC[i];
         TBox *box = new TBox(xl, yl, xr, yr);
         box->SetFillColor(18);
         box->SetLineColor(18);
         // box->Draw("same");
      }

      gr_pi_RAA_RHIC->SetMarkerStyle(24);
      gr_pi_RAA_RHIC->SetMarkerSize(1.2);
      gr_pi_RAA_RHIC->SetLineWidth(2);
      gr_pi_RAA_RHIC->SetMarkerColor(1);
      gr_pi_RAA_RHIC->SetLineColor(1);
      // gr_pi_RAA_RHIC->Draw("psame");

      gr_h_RAA_LHC->SetMarkerStyle(25);
      gr_h_RAA_LHC->SetMarkerSize(1.2);
      gr_h_RAA_LHC->SetLineWidth(2);
      gr_h_RAA_LHC->SetMarkerColor(1);
      gr_h_RAA_LHC->SetLineColor(1);
      // gr_h_RAA_LHC->Draw("psame");

      for (int i = 0; i < n_D_RAA_LHC - 3; i++)
      {
         double xl = x_D_RAA_LHC[i] - 0.1;
         double xr = x_D_RAA_LHC[i] + 0.1;
         double yl = y_D_RAA_LHC[i] - yesl_D_RAA_LHC[i];
         double yr = y_D_RAA_LHC[i] + yesr_D_RAA_LHC[i];
         TBox *box = new TBox(xl, yl, xr, yr);
         box->SetFillColor(16);
         box->SetLineColor(16);
         box->Draw("same");
      }

      gr_D_RAA_LHC->SetMarkerStyle(21);
      gr_D_RAA_LHC->SetMarkerSize(1.5);
      gr_D_RAA_LHC->SetMarkerColor(2);
      gr_D_RAA_LHC->SetLineColor(2);
      gr_D_RAA_LHC->SetLineWidth(2);
      gr_D_RAA_LHC->Draw("psame");

      // gD0err_xl[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");
      // gD0err_lz[icent]->Draw("psame");
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
         const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
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

      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");

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


      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.78, 0.75, buf, 62, 0.075, 1);
      TLatex *tex = new TLatex(8., 0.4, Form("(a)  %s%", nameCent[icent]));
      tex->SetTextFont(42);
      tex->SetTextSize(0.09);
      // tex->SetNDC(kTRUE);
      // tex->Draw("same");

   }

   //LBT_pt_cen RAA 0-10%
   TGraph *gLBTraa = new TGraph(20,LBT_pt_cen,LBT_RAA_cen); 
   gLBTraa->SetLineColor(4);
   gLBTraa->SetLineStyle(4);
   gLBTraa->SetLineWidth(2);
   gLBTraa->Draw("c");

   //LGV_pt_cen RAA 0-10%
   TGraph *gLGVraa = new TGraph(20,LGV_pt_cen,LGV_RAA_cen); 
   gLGVraa->SetLineColor(kAzure+10);
   gLGVraa->SetLineStyle(10);
   gLGVraa->SetLineWidth(2);
   gLGVraa->Draw("c");

   legend2 = new TLegend(0.48, 0.4, 0.60, 0.7);
   legend2->SetFillStyle(0);
   legend2->SetFillColor(10);
   legend2->SetBorderSize(0);
   legend2->SetTextSize(0.08);
   legend2->SetTextFont(42);
   legend2->AddEntry(gLBTraa,"LBT","l");
   legend2->AddEntry(gLGVraa,"Duke","l");
   legend2->Draw("same");


   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   drawLine(x1, y2, x2, y2, 3, 1);

   sprintf(buf, "Au+Au #sqrt{s_{NN}} = 200 GeV");
   drawLatex(0.22, 0.76, buf, 42, 0.1, 1);

   // TLine *l1 = new TLine(x1, 1.3, x2, 1.3);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");

   TBox *bx = new TBox(10.3, (1. - NbinErr[0] / Nbin[0]), 10.45, (1. + NbinErr[0] / Nbin[0]));
   bx->SetFillColor(kGreen - 6);
   bx->Draw();
   const float ppnormsys = sqrt(0.081 * 0.081 + 0.052 * 0.052); //pp normalization 0.081%
   TBox *bx = new TBox(10.45, (1. - ppnormsys), 10.58, (1. + ppnormsys));
   bx->SetFillColor(kGreen - 2);
   bx->Draw();

   sprintf(buf, "(a)");
   drawLatex(0.88, 0.3, buf, 42, 0.08, 1);

   c1->cd(2)->SetLogy(0);
   c1->cd(2);
   gPad->SetPad(0., 0.0, 1., 0.5);
   setpad(gPad, 0.16, 0.05, 0.03, 0.2);
   gPad->SetTickx();
   gPad->SetTicky(0);

   TH1F* h0 = new TH1F("", "", 1, x1, x2);
   // setHisto(h0, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   setHisto(h0, "", "p_{T} (GeV/c)", "");
   h0->GetYaxis()->SetLabelFont(42);
   h0->GetYaxis()->SetTitleFont(42);
   // y1 = 0.3;
   y2 = 1.44;
   h0->Draw();
   h0->SetMinimum(y1);
   h0->SetMaximum(y2);
   h0->GetYaxis()->CenterTitle();
   h0->GetYaxis()->SetTitleSize(0.11);
   h0->GetYaxis()->SetTitleOffset(0.55);
   h0->GetYaxis()->SetLabelSize(0.09);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetXaxis()->CenterTitle();
   h0->GetXaxis()->SetTitleSize(0.08);
   h0->GetXaxis()->SetTitleOffset(1.2);
   h0->GetXaxis()->SetLabelSize(0.08);
   h0->GetXaxis()->SetLabelOffset(0.025);

   legend2 = new TLegend(0.6, 0.55, 0.93, 0.85);
   legend2->SetFillStyle(0);
   legend2->SetFillColor(10);
   legend2->SetBorderSize(0);
   legend2->SetTextSize(0.07);
   legend2->SetTextFont(42);
   legend2->Draw("same");
   legend2->AddEntry(gr_h_RAA_LHC,"h^{#pm} 0-5% ALICE","p");
   legend2->AddEntry(gr_pi_RAA_RHIC,"#pi^{0} 0-10% PHENIX","p");

   for (int icent = 0; icent < ncent ; icent++)
   {
      if (icent != 0) continue;

      for (int i = 0; i < n_pi_RAA_RHIC - 5; i++)
      {
         double xl = x_pi_RAA_RHIC[i] * 0.95;
         double xr = x_pi_RAA_RHIC[i] * 1.05;
         double yl = y_pi_RAA_RHIC[i] - yes_pi_RAA_RHIC[i];
         double yr = y_pi_RAA_RHIC[i] + yes_pi_RAA_RHIC[i];
         TBox *box = new TBox(xl, yl, xr, yr);
         box->SetFillColor(18);
         box->SetLineColor(18);
         box->Draw("same");
      }
      for (int i = 0; i < n_h_RAA_LHC - 7; i++)
      {
         double xl = x_h_RAA_LHC[i] * 0.95;
         double xr = x_h_RAA_LHC[i] * 1.05;
         double yl = y_h_RAA_LHC[i] - yes_h_RAA_LHC[i];
         double yr = y_h_RAA_LHC[i] + yes_h_RAA_LHC[i];
         TBox *box = new TBox(xl, yl, xr, yr);
         box->SetFillColor(18);
         box->SetLineColor(18);
         box->Draw("same");
      }

      gr_pi_RAA_RHIC->SetMarkerStyle(24);
      gr_pi_RAA_RHIC->SetMarkerSize(1.2);
      gr_pi_RAA_RHIC->SetLineWidth(2);
      gr_pi_RAA_RHIC->SetMarkerColor(1);
      gr_pi_RAA_RHIC->SetLineColor(1);
      gr_pi_RAA_RHIC->Draw("psame");

      gr_h_RAA_LHC->SetMarkerStyle(25);
      gr_h_RAA_LHC->SetMarkerSize(1.2);
      gr_h_RAA_LHC->SetLineWidth(2);
      gr_h_RAA_LHC->SetMarkerColor(1);
      gr_h_RAA_LHC->SetLineColor(1);
      gr_h_RAA_LHC->Draw("psame");

      for (int i = 0; i < n_D_RAA_LHC; i++)
      {
         double xl = x_D_RAA_LHC[i] - 0.1;
         double xr = x_D_RAA_LHC[i] + 0.1;
         double yl = y_D_RAA_LHC[i] - yesl_D_RAA_LHC[i];
         double yr = y_D_RAA_LHC[i] + yesr_D_RAA_LHC[i];
         TBox *box = new TBox(xl, yl, xr, yr);
         box->SetFillColor(16);
         box->SetLineColor(16);
         // box->Draw("same");
      }

      gr_D_RAA_LHC->SetMarkerStyle(21);
      gr_D_RAA_LHC->SetMarkerSize(1.5);
      gr_D_RAA_LHC->SetMarkerColor(2);
      gr_D_RAA_LHC->SetLineColor(2);
      gr_D_RAA_LHC->SetLineWidth(2);
      // gr_D_RAA_LHC->Draw("psame");

      const float sysw = 0.15;
      // for (int i = 0; i < gD0_RAA_pperr[icent]->GetN(); i++)
      // {
      //    TBox *bx = new TBox(gD0_RAA_pperr[icent]->GetX()[i] - sysw, gD0_RAA_pperr[icent]->GetY()[i] - gD0_RAA_pperr[icent]->GetEYlow()[i], gD0_RAA_pperr[icent]->GetX()[i] + sysw, gD0_RAA_pperr[icent]->GetY()[i] + gD0_RAA_pperr[icent]->GetEYhigh()[i]);
      //    bx->SetLineColor(18);
      //    bx->SetFillColor(18);
      //    bx->SetLineWidth(1);
      //    bx->Draw("same");
      // }
      // gD0err_xl[icent]->Draw("psame");
      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");
      // gD0err_lz[icent]->Draw("psame");
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
         const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
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

      gr_pi_RAA_RHIC->Draw("psame");
      gr_h_RAA_LHC->Draw("psame");

      gD0err_xl_cp1[icent]->Draw("psame");
      gD0err_xl_cp2[icent]->Draw("psame");

      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.83, 0.38, buf, 42, 0.12, 1);
      TLatex *tex = new TLatex(8., 0.5, Form("(b)  %s%", nameCent[icent]));
      tex->SetTextFont(42);
      tex->SetTextSize(0.11);
      // tex->SetNDC(kTRUE);
      // tex->Draw("same");

   }
   legend2->Draw("same");

   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   drawLine(x1, y1, x2, y1, 3, 1);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");

   TBox *bx = new TBox(10.3, (1. - NbinErr[1] / Nbin[1]), 10.45, (1. + NbinErr[1] / Nbin[1]));
   bx->SetFillColor(kGreen - 6);
   bx->Draw();
   const float ppnormsys = sqrt(0.081 * 0.081 + 0.052 * 0.052); //pp normalization 0.081%
   TBox *bx = new TBox(10.45, (1. - ppnormsys), 10.58, (1. + ppnormsys));
   bx->SetFillColor(kGreen - 2);
   bx->Draw();

   sprintf(buf, "(b)");
   drawLatex(0.88, 0.4, buf, 42, 0.08, 1);

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

   c1->SaveAs("D0_RAA_LHC.eps");
   c1->SaveAs("D0_RAA_LHC.pdf");
   c1->SaveAs("D0_RAA_LHC.gif");
   c1->SaveAs("D0_RAA_LHC.png");
}
