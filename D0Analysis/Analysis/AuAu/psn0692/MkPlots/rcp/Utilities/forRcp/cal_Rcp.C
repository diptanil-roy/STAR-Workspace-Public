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

void cal_Rcp()
{
   char dir[250];
   char name[250];
   char title[250];
   char buf[1024];
   //TString CMD = 0;
   char CMD[250];
   TLegend* legend1;

   sprintf(dir, "pic");
   sprintf(CMD, "[ -d %s ] || mkdir -p %s", dir, dir);
   // gSystem->Exec(CMD);


   const int  ncent = 5;
   const char nameCent[ncent][250] = {"0-12%", "10-20%", "20-40%", "40-60%", "60-80%"};
   const char nameCentXL[ncent][250] = {"0_12", "10_20", "20_40", "40_60", "60_80"};
   const float scale[ncent] = {1., 1., 1., 1., 1.};
   const float Nbin[ncent] = {900.3, 591.3, 294.17, 93.6, 21.2};
   const float NbinErr[ncent] = {+71.38, +51.9, +40.64, +17.5, +6.6};

   ifstream inData;
   ofstream outData;

   double tmp;
   double pi1_x[ncent][29][2], pi1_y[ncent][29][2], pi1_ysta[ncent][29][2], pi1_ysys[ncent][29][2];//charge
   double Rcp1_x[ncent][29], Rcp1_y[ncent][29], Rcp1_ysta[ncent][29], Rcp1_ysys[ncent][29];
   for (int ic = 0; ic < ncent ; ic++)
   {
      inData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_piplus_spectra_%s.txt", nameCentXL[ic]));
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
         inData >> pi1_x[ic][i][0] >> pi1_y[ic][i][0] >> pi1_ysta[ic][i][0] >> pi1_ysys[ic][i][0];
         cout << pi1_x[ic][i][0] << " , " << pi1_y[ic][i][0] << " , " << pi1_ysta[ic][i][0] << " , " << pi1_ysys[ic][i][0] << endl;
      }
      inData.close();

      inData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_piminus_spectra_%s.txt", nameCentXL[ic]));
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
         inData >> pi1_x[ic][i][1] >> pi1_y[ic][i][1] >> pi1_ysta[ic][i][1] >> pi1_ysys[ic][i][1];
         cout << pi1_x[ic][i][1] << " , " << pi1_y[ic][i][1] << " , " << pi1_ysta[ic][i][1] << " , " << pi1_ysys[ic][i][1] << endl;
      }
      inData.close();
   }

   // cal Rcp with 60-80%
   for (int ic = 0; ic < ncent-1 ; ic++)
   {
      outData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_piplus_Rcp_%s_60_80.txt", nameCentXL[ic]));
      outData << Form("pi Rcp %s/60-80%%, pt y ystat",nameCentXL[ic]) << endl;
      for (int i = 0; i < 29; i++)
      {
        Rcp1_y[ic][i] = ((pi1_y[ic][i][0]+ pi1_y[ic][i][1])/Nbin[ic])/((pi1_y[ncent-1][i][0]+pi1_y[ncent-1][i][1])/Nbin[ncent-1]);
        Rcp1_ysta[ic][i] = Rcp1_y[ic][i]*sqrt(pow(pi1_ysta[ic][i][0]/pi1_y[ic][i][0],2)+pow(pi1_ysta[ncent-1][i][0]/pi1_y[ncent-1][i][0],2));
        outData << pi1_x[ic][i][0] << "  " << Rcp1_y[ic][i] << "  " << Rcp1_ysta[ic][i] << endl;
      }
      outData.close();
   }

   // cal Rcp with 40-60%
   for (int ic = 0; ic < ncent-2 ; ic++)
   {
      outData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_piplus_Rcp_%s_40_60.txt", nameCentXL[ic]));
      outData << Form("pi Rcp %s/40-60%%, pt y ystat",nameCentXL[ic]) << endl;
      for (int i = 0; i < 29; i++)
      {
        Rcp1_y[ic][i] = ((pi1_y[ic][i][0]+ pi1_y[ic][i][1])/Nbin[ic])/((pi1_y[ncent-2][i][0]+pi1_y[ncent-2][i][1])/Nbin[ncent-2]);
        Rcp1_ysta[ic][i] = Rcp1_y[ic][i]*sqrt(pow(pi1_ysta[ic][i][0]/pi1_y[ic][i][0],2)+pow(pi1_ysta[ncent-2][i][0]/pi1_y[ncent-2][i][0],2));
        outData << pi1_x[ic][i][0] << "  " << Rcp1_y[ic][i] << "  " << Rcp1_ysta[ic][i] << endl;
      }
      outData.close();
   }



}
