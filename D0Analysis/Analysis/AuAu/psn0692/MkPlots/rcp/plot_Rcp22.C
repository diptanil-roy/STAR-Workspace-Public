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

void plot_Rcp22()
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

   const int ncent = 4;
   const char nameCent[ncent][250] = {"0-10%", "10-20%", "20-40%", "40-60%"};
   const char nameCentXL[ncent][250] = {"0_10", "10_20", "20_40", "40_60"};
   const float scale[ncent] = {1., 1., 1., 1.};
   float Nbin[ncent] = {938.80170, 579.89409, 288.35051, 91.37100};
   float NbinErr[ncent] = {26.28048, 28.86320, 30.39279, 21.05409};


   //Read spectra
   //1. from xiaolong
   TGraphErrors* gD0err_xl[ncent];
   TGraphErrors* gD0sys_xl[ncent];
   TF1* fLevy[ncent];
   TFile* fin1 = new TFile("D0_Rcp2_pTshift.root");
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_Rcp_err_%s", nameCentXL[icent]));
      gD0sys_xl[icent] = (TGraphErrors*)fin1->Get(Form("gD0_Rcp_sys_%s", nameCentXL[icent]));
   }
   TGraphErrors* gD0Rcp2_sys_vtx = (TGraphErrors*)fin1->Get(Form("gD0Rcp2_sys_vtx_40_60"));
   fin1->Close();
   gD0Rcp2_sys_vtx->RemovePoint(10);

   //scale
   for (int icent = 0; icent < ncent; icent++)
   {
      // ScaleGraph(gD0err_xl[icent], scale[icent]);
      // ScaleGraph(gD0sys_xl[icent], scale[icent]);
   }

   //Set for Draw
   float markerSize = 2.0;
   float lineWidth = 2;
   float markerSizeScale[ncent + 1] = {0.85, 0.75, 1., 1., 1.};
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0err_xl[icent]->SetMarkerStyle(MARKERSTYLE[icent + 1]);
      gD0err_xl[icent]->SetMarkerSize(markerSize * markerSizeScale[icent + 1]);
      gD0err_xl[icent]->SetLineWidth(lineWidth);
      gD0err_xl[icent]->SetMarkerColor(COLOR[icent + 1]);
      gD0err_xl[icent]->SetLineColor(COLOR[icent + 1]);

      gD0sys_xl[icent]->SetMarkerStyle(MARKERSTYLE[icent + 1]);
      gD0sys_xl[icent]->SetMarkerSize(markerSize * markerSizeScale[icent + 1]);
      gD0sys_xl[icent]->SetLineWidth(lineWidth);
      gD0sys_xl[icent]->SetMarkerColor(COLOR[icent + 1]);
      gD0sys_xl[icent]->SetLineColor(COLOR[icent + 1]);

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
   }

   gStyle->UseCurrentStyle();

   ifstream inData;
   double tmp;
   double pi1_x[29], pi1_y[29], pi1_ysta[29], pi1_ysys[29];
   inData.open(Form("./Utilities/PionsRcp_0_12_40_80.txt"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 29; i++)
   {
      inData >> pi1_x[i] >> pi1_y[i] >> pi1_ysta[i] >> pi1_ysys[i];
      cout << pi1_x[i] << " , " << pi1_y[i] << " , " << pi1_ysta[i] << " , " << pi1_ysys[i] << endl;
   }
   inData.close();

   grpi1 = new TGraphErrors(29, pi1_x, pi1_y, 0, pi1_ysta);
   grpi1sys = new TGraphErrors(29, pi1_x, pi1_y, 0, pi1_ysys);
   grpi1->SetMarkerSize(1.5);
   grpi1->SetMarkerStyle(20 + 7);
   grpi1->SetLineWidth(2);

   double pi2_x[29], pi2_y[29], pi2_ysta[29], pi2_ysys[29];
   inData.open(Form("./Utilities/PionsRcp_20_42_40_80.txt"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 29; i++)
   {
      inData >> pi2_x[i] >> pi2_y[i] >> pi2_ysta[i] >> pi2_ysys[i];
      cout << pi2_x[i] << " , " << pi2_y[i] << " , " << pi2_ysta[i] << " , " << pi2_ysys[i] << endl;
   }
   inData.close();

   grpi2 = new TGraphErrors(29, pi2_x, pi2_y, 0, pi2_ysta);
   grpi2sys = new TGraphErrors(29, pi2_x, pi2_y, 0, pi2_ysys);
   grpi2->SetMarkerSize(1.5);
   grpi2->SetMarkerStyle(20 + 7);
   grpi2->SetLineWidth(2);

   double p1_x[26], p1_y[26], p1_ysta[26], p1_ysys[26];
   inData.open(Form("./Utilities/ProtonsRcp_0_12_40_80.txt"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 26; i++)
   {
      inData >> p1_x[i] >> p1_y[i] >> p1_ysta[i] >> p1_ysys[i];
      cout << p1_x[i] << " , " << p1_y[i] << " , " << p1_ysta[i] << " , " << p1_ysys[i] << endl;
   }
   inData.close();

   grp1 = new TGraphErrors(26, p1_x, p1_y, 0, p1_ysta);
   grp1sys = new TGraphErrors(26, p1_x, p1_y, 0, p1_ysys);
   grp1->SetMarkerSize(1.5);
   grp1->SetMarkerStyle(20 + 8);
   grp1->SetLineWidth(2);
   // grp1->SetMarkerColor(kAzure);
   // grp1->SetLineColor(kAzure);
   grp1->SetMarkerColor(1);
   grp1->SetLineColor(1);

   double p2_x[26], p2_y[26], p2_ysta[26], p2_ysys[26];
   inData.open(Form("./Utilities/ProtonsRcp_20_40_40_80.txt"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 26; i++)
   {
      inData >> p2_x[i] >> p2_y[i] >> p2_ysta[i] >> p2_ysys[i];
      cout << p2_x[i] << " , " << p2_y[i] << " , " << p2_ysta[i] << " , " << p2_ysys[i] << endl;
   }
   inData.close();

   grp2 = new TGraphErrors(26, p2_x, p2_y, 0, p2_ysta);
   grp2sys = new TGraphErrors(26, p2_x, p2_y, 0, p2_ysys);
   grp2->SetMarkerSize(1.5);
   grp2->SetMarkerStyle(20 + 8);
   grp2->SetLineWidth(2);
   // grp2->SetMarkerColor(kAzure);
   // grp2->SetLineColor(kAzure);
   grp2->SetMarkerColor(1);
   grp2->SetLineColor(1);

   //calculated Rcp pi
   TGraphErrors* cal_grpi1[ncent-1];
   double cal_pi1_x[ncent-1][29], cal_pi1_y[ncent-1][29], cal_pi1_ysta[ncent-1][29], cal_pi1_ysys[ncent-1][29];
   for(int ic =0;ic<ncent-1; ic++)
   {
     if(ic==0)
   inData.open(Form("./Utilities/forRcp/AuAu_piplus_Rcp_0_12_40_60.txt",nameCentXL[ic]));
     else
   inData.open(Form("./Utilities/forRcp/AuAu_piplus_Rcp_%s_40_60.txt",nameCentXL[ic]));

   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 29; i++)
   {
      inData >> cal_pi1_x[ic][i] >> cal_pi1_y[ic][i] >> cal_pi1_ysta[ic][i];
      cout << cal_pi1_x[ic][i] << " , " << cal_pi1_y[ic][i] << " , " << cal_pi1_ysta[ic][i] << endl;
   }
   inData.close();
   cal_grpi1[ic] = new TGraphErrors(29, cal_pi1_x[ic], cal_pi1_y[ic], 0, cal_pi1_ysta[ic]);
   cal_grpi1[ic]->SetMarkerSize(1.5);
   cal_grpi1[ic]->SetMarkerStyle(20 + 7);
   cal_grpi1[ic]->SetLineWidth(2);
   cal_grpi1[ic]->SetMarkerColor(1);
   cal_grpi1[ic]->SetLineColor(1);
   }

   //calculated Rcp ks
   TGraphErrors* cal_grKs1[ncent-1];
   double cal_Ks1_x[ncent-1][21], cal_Ks1_y[ncent-1][21], cal_Ks1_ysta[ncent-1][21], cal_Ks1_ysys[ncent-1][21];
   for(int ic =0;ic<ncent-1; ic++)
   {
     if(ic==0)
   inData.open(Form("./Utilities/forRcp/AuAu_K0s_Rcp_0_5_40_60.txt",nameCentXL[ic]));
     else
   inData.open(Form("./Utilities/forRcp/AuAu_K0s_Rcp_%s_40_60.txt",nameCentXL[ic]));

   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 21; i++)
   {
      inData >> cal_Ks1_x[ic][i] >> cal_Ks1_y[ic][i] >> cal_Ks1_ysta[ic][i];
      cout << cal_Ks1_x[ic][i] << " , " << cal_Ks1_y[ic][i] << " , " << cal_Ks1_ysta[ic][i] << endl;
   }
   inData.close();
   cal_grKs1[ic] = new TGraphErrors(21, cal_Ks1_x[ic], cal_Ks1_y[ic], 0, cal_Ks1_ysta[ic]);
   cal_grKs1[ic]->SetMarkerSize(1.0);
   cal_grKs1[ic]->SetMarkerStyle(20 + 5);
   cal_grKs1[ic]->SetLineWidth(2);
   cal_grKs1[ic]->SetMarkerColor(1);
   cal_grKs1[ic]->SetLineColor(1);
   }

   //calculated Rcp phi
   TGraphErrors* cal_grPhi1[ncent-1];
   double cal_Phi1_x[ncent-1][14], cal_Phi1_y[ncent-1][14], cal_Phi1_ysta[ncent-1][14], cal_Phi1_ysys[ncent-1][14];
   for(int ic =0;ic<ncent-1; ic++)
   {
   //   if(ic==0)
   // inData.open(Form("./Utilities/forRcp/AuAu_K0s_Rcp_0_5_60_80.txt",nameCentXL[ic]));
     // else
   inData.open(Form("./Utilities/forRcp/AuAu_phi_Rcp_%s_40_60.txt",nameCentXL[ic]));

   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 14; i++)
   {
      inData >> cal_Phi1_x[ic][i] >> cal_Phi1_y[ic][i] >> cal_Phi1_ysta[ic][i];
      cout << cal_Phi1_x[ic][i] << " , " << cal_Phi1_y[ic][i] << " , " << cal_Phi1_ysta[ic][i] << endl;
   }
   inData.close();
   cal_grPhi1[ic] = new TGraphErrors(14, cal_Phi1_x[ic], cal_Phi1_y[ic], 0, cal_Phi1_ysta[ic]);
   cal_grPhi1[ic]->SetMarkerSize(1.0);
   cal_grPhi1[ic]->SetMarkerStyle(20 + 4);
   cal_grPhi1[ic]->SetLineWidth(2);
   cal_grPhi1[ic]->SetMarkerColor(1);
   cal_grPhi1[ic]->SetLineColor(1);
   }


   //Shanshan D
   double x_Duke[5][20], y_Duke[5][20], Rcp_Duke[5][20];
   inData.open(Form("./Utilities/forSTAR-dNdpT/dNdpT_D_cen-00-10.dat"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 20; i++)
   {
      inData >> x_Duke[0][i] >> y_Duke[0][i] >> tmp >> tmp;
      cout << x_Duke[0][i] << " , " << y_Duke[0][i] << "  " << endl;
   }
   inData.close();

   inData.open(Form("./Utilities/forSTAR-dNdpT/dNdpT_D_cen-10-20.dat"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 20; i++)
   {
      inData >> x_Duke[1][i] >> y_Duke[1][i] >> tmp >> tmp;
      cout << x_Duke[1][i] << " , " << y_Duke[1][i] << "  " << endl;
   }
   inData.close();

   inData.open(Form("./Utilities/forSTAR-dNdpT/dNdpT_D_cen-20-40.dat"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 20; i++)
   {
      inData >> x_Duke[2][i] >> y_Duke[2][i] >> tmp >> tmp;
      cout << x_Duke[2][i] << " , " << y_Duke[2][i] << "  " << endl;
   }
   inData.close();

   inData.open(Form("./Utilities/forSTAR-dNdpT/dNdpT_D_cen-40-60.dat"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 20; i++)
   {
      inData >> x_Duke[3][i] >> y_Duke[3][i] >> tmp >> tmp;
      cout << x_Duke[3][i] << " , " << y_Duke[3][i] << "  " << endl;
   }
   inData.close();

   inData.open(Form("./Utilities/forSTAR-dNdpT/dNdpT_D_cen-40-80.dat"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   inData.getline(buf, 1024);
   cout << buf << endl;
   for (int i = 0; i < 20; i++)
   {
      inData >> x_Duke[4][i] >> y_Duke[4][i] >> tmp >> tmp;
      cout << x_Duke[4][i] << " , " << y_Duke[4][i] << "  " << endl;
   }
   inData.close();

   for(int ic =0; ic<5; ic++)
   {
     for (int i = 0; i < 20; i++)
     {
       Rcp_Duke[ic][i]= y_Duke[ic][i]/y_Duke[3][i];
       cout << "Rcp = " << x_Duke[ic][i] << " , " << y_Duke[ic][i] << " , " << y_Duke[4][i] << " , "<< Rcp_Duke[ic][i] << "  " << endl;
     }
   }

   TGraph* gr_Duke[5];// this is from Shanshan Cao
   for(int ic =0; ic<5; ic++)
   {
     gr_Duke[ic] = new TGraph(20, x_Duke[ic], Rcp_Duke[ic]);
     gr_Duke[ic]->SetLineWidth(3);
     gr_Duke[ic]->SetLineColor(1);
     gr_Duke[ic]->SetLineStyle(3);
   }

   //NEXT is from DUKELGV group
   double x_DukeLGV[20], y_DukeLGV[4][20], Rcp_DukeLGV[4][20];
   inData.open(Form("./Utilities/Duke/D_meson_spectra_AuAu200_forSTAR_DukeLGV_reNorm.dat"));
   if (inData.good())
   {
      cout << " file opened OK" << endl;
   }
   else
   {
      cout << " bad opening " << endl;
   }
   for(int i=0;i<3;i++)
   {
   inData.getline(buf, 1024);
   cout << buf << endl;
   }
   for (int i = 0; i < 20; i++)
   {
      inData >> x_DukeLGV[i] >> tmp >> tmp >> y_DukeLGV[0][i] >> y_DukeLGV[1][i]  >> y_DukeLGV[2][i] >> y_DukeLGV[3][i] ;
      cout << x_DukeLGV[i] << " , " << y_Duke[0][i] << "  " << y_Duke[1][i] << "  " << y_Duke[2][i] << "  " <<y_Duke[3][i] << "  " << endl;
   }
   inData.close();

   for(int ic =0; ic<4; ic++)
   {
     for (int i = 0; i < 20; i++)
     {
       Rcp_DukeLGV[ic][i]= y_DukeLGV[ic][i]/y_DukeLGV[3][i];
       cout << "Rcp = " << x_DukeLGV[i] << " , " << y_DukeLGV[ic][i] << " , " << y_DukeLGV[3][i] << " , "<< Rcp_Duke[ic][i] << "  " << endl;
     }
   }

   TGraph* gr_DukeLGV[4];// this is from Duke LGV
   for(int ic =0; ic<4; ic++)
   {
     gr_DukeLGV[ic] = new TGraph(20, x_DukeLGV, Rcp_DukeLGV[ic]);
     gr_DukeLGV[ic]->SetLineWidth(3);
     gr_DukeLGV[ic]->SetLineColor(4);
     gr_DukeLGV[ic]->SetLineStyle(2);
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
   double x2 = 8.6;
   double y1 = 0.05;
   double y2 = 1.65;

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
   gD0Rcp2_sys_vtx->Draw("e3");
   gD0Rcp2_sys_vtx->SetLineColor(16);
   gD0Rcp2_sys_vtx->SetFillColor(16);

   // grpi1->Draw("Psame");
   grpi1->RemovePoint(29);
   grpi1->RemovePoint(28);
   grpi1->RemovePoint(27);
   grpi1->RemovePoint(26);
   // grp1->Draw("Psame");
   gr_Duke[0]->Draw("csame");
   gr_DukeLGV[0]->Draw("csame");
   // cal_grpi1[0]->Draw("Psame");
   cal_grKs1[0]->RemovePoint(21);
   cal_grKs1[0]->RemovePoint(20);
   cal_grKs1[0]->RemovePoint(19);
   // cal_grKs1[0]->Draw("Psame");
   // cal_grPhi1[0]->Draw("Psame");

   const float sysw = 0.08;
   // for (int i = 0; i < grpi1sys->GetN() - 2; i++)
   // {
   //    const float sysl = grpi1sys->GetY()[i] * 0.05;
   //    TLine *llw = new TLine(grpi1sys->GetX()[i] - sysw, grpi1sys->GetY()[i] - grpi1sys->GetEY()[i], grpi1sys->GetX()[i] + sysw, grpi1sys->GetY()[i] - grpi1sys->GetEY()[i]);
   //    llw->SetLineWidth(2);
   //    llw->SetLineColor(1);
   //    llw->Draw("same");
   //    TLine *lhi = new TLine(grpi1sys->GetX()[i] - sysw, grpi1sys->GetY()[i] + grpi1sys->GetEY()[i], grpi1sys->GetX()[i] + sysw, grpi1sys->GetY()[i] + grpi1sys->GetEY()[i]);
   //    lhi->SetLineWidth(2);
   //    lhi->SetLineColor(1);
   //    lhi->Draw("same");
   //    TLine *lv = new TLine(grpi1sys->GetX()[i] - sysw, grpi1sys->GetY()[i] - grpi1sys->GetEY()[i], grpi1sys->GetX()[i] - sysw, grpi1sys->GetY()[i] - grpi1sys->GetEY()[i] + sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   //    TLine *lv = new TLine(grpi1sys->GetX()[i] + sysw, grpi1sys->GetY()[i] - grpi1sys->GetEY()[i], grpi1sys->GetX()[i] + sysw, grpi1sys->GetY()[i] - grpi1sys->GetEY()[i] + sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   //    TLine *lv = new TLine(grpi1sys->GetX()[i] - sysw, grpi1sys->GetY()[i] + grpi1sys->GetEY()[i], grpi1sys->GetX()[i] - sysw, grpi1sys->GetY()[i] + grpi1sys->GetEY()[i] - sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   //    TLine *lv = new TLine(grpi1sys->GetX()[i] + sysw, grpi1sys->GetY()[i] + grpi1sys->GetEY()[i], grpi1sys->GetX()[i] + sysw, grpi1sys->GetY()[i] + grpi1sys->GetEY()[i] - sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   // }

   const float sysw = 0.08;
   // for (int i = 0; i < grp1sys->GetN() - 2; i++)
   // {
   //    const float sysl = grp1sys->GetY()[i] * 0.05;
   //    TLine *llw = new TLine(grp1sys->GetX()[i] - sysw, grp1sys->GetY()[i] - grp1sys->GetEY()[i], grp1sys->GetX()[i] + sysw, grp1sys->GetY()[i] - grp1sys->GetEY()[i]);
   //    llw->SetLineWidth(2);
   //    llw->SetLineColor(1);
   //    llw->Draw("same");
   //    TLine *lhi = new TLine(grp1sys->GetX()[i] - sysw, grp1sys->GetY()[i] + grp1sys->GetEY()[i], grp1sys->GetX()[i] + sysw, grp1sys->GetY()[i] + grp1sys->GetEY()[i]);
   //    lhi->SetLineWidth(2);
   //    lhi->SetLineColor(1);
   //    lhi->Draw("same");
   //    TLine *lv = new TLine(grp1sys->GetX()[i] - sysw, grp1sys->GetY()[i] - grp1sys->GetEY()[i], grp1sys->GetX()[i] - sysw, grp1sys->GetY()[i] - grp1sys->GetEY()[i] + sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   //    TLine *lv = new TLine(grp1sys->GetX()[i] + sysw, grp1sys->GetY()[i] - grp1sys->GetEY()[i], grp1sys->GetX()[i] + sysw, grp1sys->GetY()[i] - grp1sys->GetEY()[i] + sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   //    TLine *lv = new TLine(grp1sys->GetX()[i] - sysw, grp1sys->GetY()[i] + grp1sys->GetEY()[i], grp1sys->GetX()[i] - sysw, grp1sys->GetY()[i] + grp1sys->GetEY()[i] - sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   //    TLine *lv = new TLine(grp1sys->GetX()[i] + sysw, grp1sys->GetY()[i] + grp1sys->GetEY()[i], grp1sys->GetX()[i] + sysw, grp1sys->GetY()[i] + grp1sys->GetEY()[i] - sysl);
   //    lv->SetLineWidth(2);
   //    lv->SetLineColor(1);
   //    lv->Draw("same");
   // }

   // legend1 = new TLegend(0.81, 0.25, 0.9, 0.75);
   legend1 = new TLegend(0.78, 0.3, 0.93, 0.85);
   legend1->SetFillStyle(0);
   legend1->SetFillColor(10);
   legend1->SetBorderSize(0);
   legend1->SetTextSize(0.08);
   legend1->SetTextFont(42);
   // legend1->AddEntry(grpi1, "pi 0-12%/40-80%", "p");
   // legend1->AddEntry(grp1, "p 0-12%/40-80%", "p");
   // legend1->AddEntry(grpi1, "#pi, 0-12%", "p");
   // legend1->AddEntry(cal_grKs1[0], "K_{s}^{0}, 0-5%", "p");
   // legend1->AddEntry(cal_grPhi1[0], "#phi, 0-10%", "p");
   // legend1->AddEntry(grp1, "p, 0-12%", "p");
   // legend1->AddEntry(gD0err_xl[0], "D^{0}, 0-10%", "p");
   legend1->AddEntry(gr_Duke[0], "LBT", "pl");
   legend1->AddEntry(gr_DukeLGV[0], "Duke", "pl");
   legend1->Draw("same");

   for (int icent = 0; icent < ncent - 1; icent++)
   {
      if (icent != 0) continue;
      gD0err_xl[icent]->Draw("psame");
      //draw systematic error

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

      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.78, 0.75, buf, 62, 0.075, 1);
   TLatex *tex = new TLatex(0.7, 0.2, Form("(a)   %s%", nameCent[icent]));
   tex->SetTextFont(42);
   tex->SetTextSize(0.085);
   // tex->SetNDC(kTRUE);
   tex->Draw("same");

   }
   // drawLines(x1, y1, x2, y2, 2, 1);
   drawLine(x1, y1, x1, y2, 3, 1);
   drawLine(x2, y1, x2, y2, 3, 1);
   drawLine(x1, y2, x2, y2, 3, 1);

   sprintf(buf, "Au+Au #sqrt{s_{NN}} = 200 GeV");
   drawLatex(0.22, 0.76, buf, 42, 0.11, 1);

   // TLine *l1 = new TLine(x1, 1.3, x2, 1.3);
   TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   l1->SetLineWidth(2);
   l1->SetLineColor(1);
   l1->SetLineStyle(2);
   l1->Draw("same");

   TBox *bx = new TBox(8.3,(1.-NbinErr[0]/Nbin[0]),8.45,(1.+NbinErr[0]/Nbin[0]));
   bx->SetFillColor(kGreen-6);
   bx->Draw();
   TBox *bx = new TBox(8.45,(1.-NbinErr[3]/Nbin[3]),8.58,(1.+NbinErr[3]/Nbin[3]));
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
   y2 = 1.44;
   h0->Draw();
   h0->SetMinimum(y1);
   h0->SetMaximum(y2);
   h0->GetYaxis()->CenterTitle();
   h0->GetYaxis()->SetTitleSize(0.13);
   h0->GetYaxis()->SetTitleOffset(0.5);
   h0->GetYaxis()->SetLabelSize(0.13);
   gD0Rcp2_sys_vtx->Draw("e3");
   gr_Duke[1]->Draw("csame");
   gr_DukeLGV[1]->Draw("csame");
   cal_grpi1[1]->RemovePoint(29);
   cal_grpi1[1]->RemovePoint(28);
   cal_grpi1[1]->RemovePoint(27);
   cal_grpi1[1]->RemovePoint(26);
   // cal_grpi1[1]->Draw("Psame");
   cal_grKs1[1]->RemovePoint(21);
   cal_grKs1[1]->RemovePoint(20);
   cal_grKs1[1]->RemovePoint(19);
   // cal_grKs1[1]->Draw("Psame");
   // cal_grPhi1[1]->Draw("Psame");
   for (int icent = 0; icent < ncent - 1; icent++)
   {
      if (icent != 1) continue;
      gD0err_xl[icent]->Draw("psame");
      //draw systematic error

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
      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.83, 0.38, buf, 42, 0.12, 1);
   TLatex *tex = new TLatex(0.7, 0.2, Form("(b)   %s%", nameCent[icent]));
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

   TBox *bx = new TBox(8.3,(1.-NbinErr[1]/Nbin[1]),8.45,(1.+NbinErr[1]/Nbin[1]));
   bx->SetFillColor(kGreen-6);
   bx->Draw();
   TBox *bx = new TBox(8.45,(1.-NbinErr[3]/Nbin[3]),8.58,(1.+NbinErr[3]/Nbin[3]));
   bx->SetFillColor(kGreen-2);
   bx->Draw();
   c1->cd(3)->SetLogy(0);
   // gPad->SetPad(0., 0.28, 1., 0.48);
   gPad->SetPad(0., 0.0, 1., 0.36);
   setpad(gPad, 0.16, 0.05, 0.03, 0.23);
   gPad->SetTickx();
   gPad->SetTicky(0);

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
   gD0Rcp2_sys_vtx->Draw("e3");
   gr_Duke[2]->Draw("csame");
   gr_DukeLGV[2]->Draw("csame");
   grpi2->RemovePoint(29);
   grpi2->RemovePoint(28);
   grpi2->RemovePoint(27);
   grpi2->RemovePoint(26);
   // cal_grpi1[1]->RemovePoint(29);
   // cal_grpi1[1]->RemovePoint(28);
   // cal_grpi1[1]->RemovePoint(27);
   // cal_grpi1[1]->RemovePoint(26);
   // cal_grpi1[1]->Draw("Psame");
   cal_grKs1[1]->RemovePoint(21);
   cal_grKs1[1]->RemovePoint(20);
   cal_grKs1[1]->RemovePoint(19);
   // cal_grKs1[1]->Draw("Psame");
   // cal_grPhi1[1]->Draw("Psame");
   for (int icent = 0; icent < ncent - 1; icent++)
   {
      if (icent != 2) continue;
      gD0err_xl[icent]->Draw("psame");
      //draw systematic error

      // grpi2->Draw("Psame");
      // grp2->Draw("Psame");

      const float sysw = 0.08;
      // for (int i = 0; i < grpi2sys->GetN() - 2; i++)
      // {
      //    const float sysl = grpi2sys->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(grpi2sys->GetX()[i] - sysw, grpi2sys->GetY()[i] - grpi2sys->GetEY()[i], grpi2sys->GetX()[i] + sysw, grpi2sys->GetY()[i] - grpi2sys->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(1);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(grpi2sys->GetX()[i] - sysw, grpi2sys->GetY()[i] + grpi2sys->GetEY()[i], grpi2sys->GetX()[i] + sysw, grpi2sys->GetY()[i] + grpi2sys->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(1);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(grpi2sys->GetX()[i] - sysw, grpi2sys->GetY()[i] - grpi2sys->GetEY()[i], grpi2sys->GetX()[i] - sysw, grpi2sys->GetY()[i] - grpi2sys->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(grpi2sys->GetX()[i] + sysw, grpi2sys->GetY()[i] - grpi2sys->GetEY()[i], grpi2sys->GetX()[i] + sysw, grpi2sys->GetY()[i] - grpi2sys->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(grpi2sys->GetX()[i] - sysw, grpi2sys->GetY()[i] + grpi2sys->GetEY()[i], grpi2sys->GetX()[i] - sysw, grpi2sys->GetY()[i] + grpi2sys->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(grpi2sys->GetX()[i] + sysw, grpi2sys->GetY()[i] + grpi2sys->GetEY()[i], grpi2sys->GetX()[i] + sysw, grpi2sys->GetY()[i] + grpi2sys->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      // }

      const float sysw = 0.08;
      // for (int i = 0; i < grp2sys->GetN() - 2; i++)
      // {
      //    const float sysl = grp2sys->GetY()[i] * 0.05;
      //    TLine *llw = new TLine(grp2sys->GetX()[i] - sysw, grp2sys->GetY()[i] - grp2sys->GetEY()[i], grp2sys->GetX()[i] + sysw, grp2sys->GetY()[i] - grp2sys->GetEY()[i]);
      //    llw->SetLineWidth(2);
      //    llw->SetLineColor(1);
      //    llw->Draw("same");
      //    TLine *lhi = new TLine(grp2sys->GetX()[i] - sysw, grp2sys->GetY()[i] + grp2sys->GetEY()[i], grp2sys->GetX()[i] + sysw, grp2sys->GetY()[i] + grp2sys->GetEY()[i]);
      //    lhi->SetLineWidth(2);
      //    lhi->SetLineColor(1);
      //    lhi->Draw("same");
      //    TLine *lv = new TLine(grp2sys->GetX()[i] - sysw, grp2sys->GetY()[i] - grp2sys->GetEY()[i], grp2sys->GetX()[i] - sysw, grp2sys->GetY()[i] - grp2sys->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(grp2sys->GetX()[i] + sysw, grp2sys->GetY()[i] - grp2sys->GetEY()[i], grp2sys->GetX()[i] + sysw, grp2sys->GetY()[i] - grp2sys->GetEY()[i] + sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(grp2sys->GetX()[i] - sysw, grp2sys->GetY()[i] + grp2sys->GetEY()[i], grp2sys->GetX()[i] - sysw, grp2sys->GetY()[i] + grp2sys->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      //    TLine *lv = new TLine(grp2sys->GetX()[i] + sysw, grp2sys->GetY()[i] + grp2sys->GetEY()[i], grp2sys->GetX()[i] + sysw, grp2sys->GetY()[i] + grp2sys->GetEY()[i] - sysl);
      //    lv->SetLineWidth(2);
      //    lv->SetLineColor(1);
      //    lv->Draw("same");
      // }

      legend1 = new TLegend(0.5, 0.28, 0.9, 0.38);
      legend1->SetFillStyle(0);
      legend1->SetFillColor(10);
      legend1->SetBorderSize(0);
      legend1->SetTextSize(0.12);
      legend1->SetTextFont(132);
      legend1->SetNColumns(2);
      legend1->AddEntry(grpi2, "#pi 20-42%", "p");
      // legend1->AddEntry(grp1, "p 20-42%", "p");
      // legend1->Draw("same");


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
      sprintf(buf, Form("%s%%", nameCent[icent]));
      // drawLatex(0.83, 0.6, buf, 42, 0.085, 1);
   TLatex *tex = new TLatex(0.7, 0.2, Form("(c)   %s%", nameCent[icent]));
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

   TBox *bx = new TBox(8.3,(1.-NbinErr[2]/Nbin[2]),8.45,(1.+NbinErr[2]/Nbin[2]));
   bx->SetFillColor(kGreen-6);
   bx->Draw();
   TBox *bx = new TBox(8.45,(1.-NbinErr[3]/Nbin[3]),8.58,(1.+NbinErr[3]/Nbin[3]));
   bx->SetFillColor(kGreen-2);
   bx->Draw();
   // c1->cd(4)->SetLogy(0);
   // gPad->SetPad(0., 0., 1., 0.28);
   // setpad(gPad, 0.16, 0.05, 0.03, 0.28);
   // gPad->SetTickx(1);
   // gPad->SetTicky(1);
   //
   // TH1F* h1 = new TH1F("", "", 1, x1, x2);
   // setHisto(h1, "", "p_{T} (GeV/c)", Form("R_{cp} (/%s)", nameCent[ncent - 1]));
   // h1->Draw();
   // h1->GetYaxis()->SetLabelFont(42);
   // h1->GetYaxis()->SetTitleFont(42);
   // h1->GetXaxis()->SetLabelFont(42);
   // h1->GetXaxis()->SetTitleFont(42);
   // h1->SetMinimum(y1);
   // h1->SetMaximum(y2);
   // h1->GetXaxis()->CenterTitle();
   // h1->GetYaxis()->CenterTitle();
   // h1->GetXaxis()->SetTitleSize(0.11);
   // h1->GetYaxis()->SetTitleSize(0.095);
   // h1->GetXaxis()->SetTitleOffset(1.1);
   // h1->GetYaxis()->SetTitleOffset(0.7);
   // h1->GetXaxis()->SetLabelSize(0.085);
   // h1->GetXaxis()->SetLabelOffset(0.025);
   // h1->GetYaxis()->SetLabelSize(0.07);
   // // gr_Duke[3]->Draw("csame");
   // for (int icent = 0; icent < ncent - 1; icent++)
   // {
   //    if (icent != 3) continue;
   //    gD0err_xl[icent]->Draw("psame");
   //    //draw systematic error
   //
   //    const float sysw = 0.15;
   //    for (int i = 0; i < gD0sys_xl[icent]->GetN(); i++)
   //    {
   //       const float sysl = gD0sys_xl[icent]->GetY()[i] * 0.05;
   //       TLine *llw = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i]);
   //       llw->SetLineWidth(2);
   //       llw->SetLineColor(1);
   //       llw->Draw("same");
   //       TLine *lhi = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i]);
   //       lhi->SetLineWidth(2);
   //       lhi->SetLineColor(1);
   //       lhi->Draw("same");
   //       TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
   //       lv->SetLineWidth(2);
   //       lv->SetLineColor(1);
   //       lv->Draw("same");
   //       TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] - gD0sys_xl[icent]->GetEY()[i] + sysl);
   //       lv->SetLineWidth(2);
   //       lv->SetLineColor(1);
   //       lv->Draw("same");
   //       TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] - sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
   //       lv->SetLineWidth(2);
   //       lv->SetLineColor(1);
   //       lv->Draw("same");
   //       TLine *lv = new TLine(gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i], gD0sys_xl[icent]->GetX()[i] + sysw, gD0sys_xl[icent]->GetY()[i] + gD0sys_xl[icent]->GetEY()[i] - sysl);
   //       lv->SetLineWidth(2);
   //       lv->SetLineColor(1);
   //       lv->Draw("same");
   //    }
   //    sprintf(buf, Form("%s%%", nameCent[icent]));
   //    drawLatex(0.83, 0.53, buf, 42, 0.075, 1);
   //
   // }
   // // drawLines(x1, y1, x2, y2, 2.5, 1);
   // drawLine(x1, y1, x1, y2, 3, 1);
   // drawLine(x2, y1, x2, y2, 3, 1);
   // drawLine(x1, y1, x2, y1, 3, 1);
   //
   // sprintf(buf, "Au+Au #sqrt{s_{NN}} = 200 GeV");
   // // drawLatex(0.18, 0.88, buf, 42, 0.045, 1);
   // TLine *l1 = new TLine(x1, 1.0, x2, 1.0);
   // l1->SetLineWidth(2);
   // l1->SetLineColor(1);
   // l1->SetLineStyle(2);
   // l1->Draw("same");
   //
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
   sprintf(buf, Form("R_{cp} (/%s%%)", nameCent[ncent - 1]));
   TLatex *tex = new TLatex(0.8, 0.45, buf);
   tex->SetTextFont(42);
   tex->SetTextSize(0.7);
   tex->SetTextAngle(90);
   tex->SetNDC();
   tex->Draw("same");

   c1->SaveAs("D0_Rcp22.eps");
   c1->SaveAs("D0_Rcp22.pdf");
   c1->SaveAs("D0_Rcp22.gif");
   c1->SaveAs("D0_Rcp22.png");
}
