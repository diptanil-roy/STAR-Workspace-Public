#include "../anaCuts.h"
#include "../myFunction.h"
#include "../myConst.h"

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

void plot_sys()
{

   char name[1024];
   TGraphErrors* grEff[ncent];
   TGraphErrors* grEff1[ncent];
   TGraphErrors* grEff2[ncent];
   TGraphErrors* grEff3[ncent];
   TGraphErrors* grEff4[ncent];

   TGraphErrors* grEffRatio1[ncent];
   TGraphErrors* grEffRatio2[ncent];
   TGraphErrors* grEffRatio3[ncent];
   TGraphErrors* grEffRatio4[ncent];

    float y[ncent][npt], yerr[ncent][npt];
    float y1[ncent][npt], y1err[ncent][npt];
    float y2[ncent][npt], y2err[ncent][npt];
    float y3[ncent][npt], y3err[ncent][npt];
    float y4[ncent][npt], y4err[ncent][npt];
    float yRatio1[ncent][npt], yRatio1err[ncent][npt];
    float yRatio2[ncent][npt], yRatio2err[ncent][npt];
    float yRatio3[ncent][npt], yRatio3err[ncent][npt];
    float yRatio4[ncent][npt], yRatio4err[ncent][npt];

    float ptMean[npt];
    float ptErr[npt];
    for(int ipt =0; ipt<npt ; ipt++)
    {
      ptMean[ipt] = 0.5*(nptbin[ipt] + nptbin[ipt+1]);
      ptErr[ipt] = nptbin[ipt+1]- ptMean[ipt];
    }

    ifstream in;
    ofstream out;

   for (int icent = 0; icent < ncent; icent++)
   {
      // grEff[icent] = new TGraph(Form("../default/data/re_yield_%s.csv", nameCent1[icent]), "%lg %lg", ",");
      // grEff1[icent] = new TGraph(Form("data/yield_%s.csv", nameCent1[icent]), "%lg %lg", ",");
      // grEff2[icent] = new TGraph(Form("../ptCut2/data/yield_%s.csv", nameCent1[icent]), "%lg %lg", ",");
      
     gStyle->SetOptStat(0);
        in.open(Form("../default/data/re_yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> y[icent][ipt] >> yerr[icent][ipt];
        in.close();

        in.open(Form("../topoCut1/data/yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> y1[icent][ipt] >> y1err[icent][ipt];
        in.close();

        in.open(Form("../topoCut2/data/yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> y2[icent][ipt] >> y2err[icent][ipt];
        in.close();

        in.open(Form("../ptCut1/data/yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> y3[icent][ipt] >> y3err[icent][ipt];
        in.close();

        in.open(Form("../ptCut2/data/yield_%s.txt",nameCent1[icent]));
        for(int ipt=0; ipt<npt; ipt++) in >> y4[icent][ipt] >> y4err[icent][ipt];
        in.close();

        grEff[icent] = new TGraphErrors(npt,ptMean,y[icent],0,yerr[icent]);
        grEff1[icent] = new TGraphErrors(npt,ptMean,y1[icent],0,y1err[icent]);
        grEff2[icent] = new TGraphErrors(npt,ptMean,y2[icent],0,y2err[icent]);
        grEff3[icent] = new TGraphErrors(npt,ptMean,y3[icent],0,y3err[icent]);
        grEff4[icent] = new TGraphErrors(npt,ptMean,y4[icent],0,y4err[icent]);

        for(int ipt=0; ipt<npt; ipt++) 
        {
          yRatio1[icent][ipt]  = y1[icent][ipt]/ y[icent][ipt];
          yRatio1err[icent][ipt]  = y1err[icent][ipt]/ y[icent][ipt];
          yRatio2[icent][ipt] = y2[icent][ipt]/y[icent][ipt];
          yRatio2err[icent][ipt] = y2err[icent][ipt]/y[icent][ipt];
          yRatio3[icent][ipt] = y3[icent][ipt]/y[icent][ipt];
          yRatio3err[icent][ipt] = y3err[icent][ipt]/y[icent][ipt];
          yRatio4[icent][ipt] = y4[icent][ipt]/y[icent][ipt];
          yRatio4err[icent][ipt] = y4err[icent][ipt]/y[icent][ipt];
        }
        grEffRatio1[icent] = new TGraphErrors(npt,ptMean,yRatio1[icent],0,yRatio1err[icent]);
        grEffRatio2[icent] = new TGraphErrors(npt,ptMean,yRatio2[icent],0,yRatio2err[icent]);
        grEffRatio3[icent] = new TGraphErrors(npt,ptMean,yRatio3[icent],0,yRatio3err[icent]);
        grEffRatio4[icent] = new TGraphErrors(npt,ptMean,yRatio4[icent],0,yRatio4err[icent]);

   }

   // float ymin[ncent] = 1e3;
   // float ymax[ncent] = 0.99e7;

   //plot
   TCanvas* c1 = new TCanvas("c1", "A Canvas", 10, 10, 800, 800);
   setPad(c1);

   float small = 0;
   c1->Divide(1, 2, small, small);

   for (int icent = 0; icent < ncent; icent++)
   {

      c1->cd(1);
      gPad->SetPad(0., 0.3, 1.0, 1.);
      setpad(gPad, 0.12, 0.1, 0.06, 0.0);

      float ymin = 1e3;
      float ymax = 0.99e7;
      h0 = new TH1F("", "", 1, 0, 8);
      // h0= new TH1F("","",1,0,10);
      h0->Draw();
      ymin = grEff[icent]->GetY()[npt-1] / 2;
      ymax = grEff[icent]->GetY()[1] * 2;
      h0->SetMinimum(ymin),
      h0->SetMaximum(ymax);
      // setHisto(h0,"","p_{T} (GeV/c)", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
      // setHisto(h0, "", "p_{T} (GeV/c)", "dN/#varepsilon");
      setHisto(h0, "", "p_{T} (GeV/c)", "dN/#epsilon");
      legend1 = new TLegend(0.68, 0.43, 0.9, 0.93, Form("%s", nameCent[icent]));
      legend1->SetFillStyle(0);
      legend1->SetFillColor(10);
      legend1->SetBorderSize(0);
      legend1->SetTextSize(0.04);
      legend1->SetTextFont(132);
      grEff[icent]->SetMarkerStyle(24);
      grEff[icent]->SetMarkerColor(1);
      grEff[icent]->SetLineColor(1);
      grEff1[icent]->SetMarkerStyle(24);
      grEff1[icent]->SetMarkerColor(2);
      grEff1[icent]->SetLineColor(2);
      grEff2[icent]->SetMarkerStyle(24);
      grEff2[icent]->SetMarkerColor(4);
      grEff2[icent]->SetLineColor(4);
      grEff3[icent]->SetMarkerStyle(24);
      grEff3[icent]->SetMarkerColor(kGreen+1);
      grEff3[icent]->SetLineColor(kGreen+1);
      grEff4[icent]->SetMarkerStyle(24);
      grEff4[icent]->SetMarkerColor(kOrange-1);
      grEff4[icent]->SetLineColor(kOrange-1);
      grEff[icent]->Draw("psame");
      grEff1[icent]->Draw("psame");
      grEff2[icent]->Draw("psame");
      grEff3[icent]->Draw("psame");
      grEff4[icent]->Draw("psame");
      legend1->AddEntry(grEff[icent], "default", "p");
      legend1->AddEntry(grEff1[icent], "eff 50%", "p");
      legend1->AddEntry(grEff2[icent], "eff 150%", "p");
      legend1->AddEntry(grEff3[icent], "pT cut1", "p");
      legend1->AddEntry(grEff4[icent], "pT cut2", "p");
      legend1->Draw("same");
      gPad->SetLogy();

      drawLines(0,ymin,8,ymax,3,1);

      c1->cd(2);
      gPad->SetPad(0., 0., 1.0, 0.3);
      setpad(gPad, 0.12, 0.1, 0.0, 0.3);

      // float ymin = 0.7;
      // float ymax = 1.3;
      float ymin = 0.8;
      float ymax = 1.2;
      h0 = new TH1F("", "", 1, 0, 8);
      h0->Draw();
      h0->SetMinimum(ymin),
      h0->SetMaximum(ymax);
      // setHisto(h0,"","p_{T} (GeV/c)", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
      setHisto(h0, "", "p_{T} (GeV/c)", "Ratio");
      c1->cd(2)->SetGridx();
      c1->cd(2)->SetGridy();
      c1->cd(2)->SetTickx();
      c1->cd(2)->SetTicky();
      h0->GetXaxis()->SetLabelSize(0.08);
      h0->GetYaxis()->SetLabelSize(0.08);
      h0->GetXaxis()->SetTitleSize(0.09);
      h0->GetYaxis()->SetTitleSize(0.09);
      
      grEffRatio1[icent]->SetMarkerStyle(24);
      grEffRatio1[icent]->SetMarkerColor(2);
      grEffRatio1[icent]->SetLineColor(2);
      grEffRatio2[icent]->SetMarkerStyle(24);
      grEffRatio2[icent]->SetMarkerColor(4);
      grEffRatio2[icent]->SetLineColor(4);
      grEffRatio3[icent]->SetMarkerStyle(24);
      grEffRatio3[icent]->SetMarkerColor(kGreen+1);
      grEffRatio3[icent]->SetLineColor(kGreen+1);
      grEffRatio4[icent]->SetMarkerStyle(24);
      grEffRatio4[icent]->SetMarkerColor(kOrange-1);
      grEffRatio4[icent]->SetLineColor(kOrange-1);
      grEffRatio1[icent]->Draw("psame");
      grEffRatio2[icent]->Draw("psame");
      grEffRatio3[icent]->Draw("psame");
      grEffRatio4[icent]->Draw("psame");
      gPad->SetLogy(0);
      drawLines(0,ymin,8,ymax,3,1);

      sprintf(name, Form("./pic/D0_spectra_%s.eps", nameCent1[icent]));
      c1->SaveAs(name);
   }


}
