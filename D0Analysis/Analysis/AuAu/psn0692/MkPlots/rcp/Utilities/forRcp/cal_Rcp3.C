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

void cal_Rcp3()
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


   const int  ncent = 9;
   const char nameCent[ncent][250] =  {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","0-10%","0-5%"}; 
   const char nameCentXL[ncent][250] = {"70_80","60_70","50_60","40_50","30_40","20_30","10_20","0_10","0_5"};
   const float Nbin[ncent] = {12.91047, 29.32408, 62.25545, 120.65795, 215.95050, 360.58912, 579.89409, 938.80170, 1042.75372}; 
   //40-60%/ 60-80% // 91.37100, 21.37396
   // const float NbinErr[ncent] = {};

   ifstream inData;
   ofstream outData;

   double tmp;
   double pi1_x[ncent][14], pi1_xerr[ncent][14], pi1_y[ncent][14], pi1_ysta[ncent][14], pi1_ysys[ncent][14];//charge
   double Rcp1_x[ncent][14], Rcp1_y[ncent][14], Rcp1_ysta[ncent][14], Rcp1_ysys[ncent][14];
   for (int ic = 0; ic < ncent; ic++)
   {
      inData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_phi_MTspectra_%s.txt", nameCentXL[ic]));
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
         inData >> pi1_x[ic][i] >> pi1_xerr[ic][i]>> pi1_y[ic][i] >> pi1_ysta[ic][i] ;
         cout << pi1_x[ic][i] << " , "<< pi1_xerr[ic][i] << " , " << pi1_y[ic][i] << " , " << pi1_ysta[ic][i] << endl;
      }
      inData.close();

   }

   // cal Spectra with 60-80%
   outData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_phi_MTspectra_60_80.txt"));
   outData << Form("phi spectra  60-80%%, pt mt y ystat",nameCentXL[ic]) << endl;
   for (int i = 0; i < 14; i++)
   {
     double p1_y = (pi1_y[0][i] + pi1_y[1][i])/2.;
     double p1_ysta = p1_y*sqrt(pow(pi1_ysta[0][i]/pi1_y[0][i],2)+pow(pi1_ysta[1][i]/pi1_y[1][i],2));
     outData  << pi1_x[0][i] << "  "<< pi1_xerr[0][i] << "  " << p1_y << "  " << p1_ysta << endl;
   }
   outData.close();

   // cal Spectra with 40-60%
   outData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_phi_MTspectra_40_60.txt"));
   outData << Form("phi spectra  40-60%%, pt mt y ystat",nameCentXL[ic]) << endl;
   for (int i = 0; i < 14; i++)
   {
     double p1_y = (pi1_y[2][i] + pi1_y[3][i])/2.;
     double p1_ysta = p1_y*sqrt(pow(pi1_ysta[2][i]/pi1_y[2][i],2)+pow(pi1_ysta[3][i]/pi1_y[3][i],2));
     outData  << pi1_x[0][i] << "  "<< pi1_xerr[0][i] << "  " << p1_y << "  " << p1_ysta << endl;
   }
   outData.close();

   // cal Spectra with 20-40%
   outData.open(Form("/Users/guannanxie/work/D0Analysis/Physics2/SL16d/paperPlots/Utilities/forRcp/AuAu_phi_MTspectra_20_40.txt"));
   outData << Form("phi spectra  20-40%%, pt mt y ystat",nameCentXL[ic]) << endl;
   for (int i = 0; i < 14; i++)
   {
     double p1_y = (pi1_y[4][i] + pi1_y[5][i])/2.;
     double p1_ysta = p1_y*sqrt(pow(pi1_ysta[4][i]/pi1_y[4][i],2)+pow(pi1_ysta[5][i]/pi1_y[5][i],2));
     outData  << pi1_x[0][i] << "  "<< pi1_xerr[0][i] << "  " << p1_y << "  " << p1_ysta << endl;
   }
   outData.close();

}
