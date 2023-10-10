#include "../../myFunction.h"
#include "../../myConst.h"
#include "../../anaCuts.h"

void find_Assym_Band(double &raa1, double &raa2, double &lowEdge, double &highEdge, double auau, double auauErr, double pp, double ppErr)
{
   const Int_t Nevt = 1e6;//1e7
   const Double_t Err = 0.683;
   const Double_t PDF_lowEdge = (1. - Err) / 2.;
   const Double_t PDF_highEdge = 1. - (1. - Err) / 2.;
   const Double_t PDF_median = 0.5;

   // if ((auauErr - auau) > 1e-8) auauErr = auau;
   // if ((ppErr - pp) > 1e-8) ppErr = pp;

   TH1D *hRAA = new TH1D(Form("RAA"), "", 2000, -4., 4.0);
   TH1D *hAA = new TH1D(Form("AA"), "", 2000, auau - 5 * auauErr, auau + 5 * auauErr);
   TH1D *hPP = new TH1D(Form("PP"), "", 2000, pp - 5 * ppErr, pp + 5 * ppErr);

   cout << " Histograms declared. Start processing ... " << endl;
   TRandom3 *gRandom = new TRandom3();
   for (int i = 0; i < Nevt; i++)
   {
      if (i % 100000 == 0)
         cout << " Processing " << i << "-th event ... " << endl;
      double pp_a = gRandom->Gaus(pp, ppErr);
      double AA_a = gRandom->Gaus(auau, auauErr);
      if(pp_a < 0. ) pp_a=1e-16;
      double RAA = AA_a / pp_a;
      hRAA->Fill(RAA);
      hAA->Fill(AA_a);
      hPP->Fill(pp_a);
   }

   double  r = (auauErr / auau) / (ppErr / pp);
   double  mu = 0.;
   double  mu_low = 0.;
   double  mu_high = 999.;
   double  mu_gaus = auau / pp;

   double int_frac = 0.;
   double int_frac_1 = 0.;
   for (int ip = 0; ip < 2000; ip++)
   {
      int_frac += hRAA->GetBinContent(ip + 1) / (double)Nevt;
      if (ip > 0) int_frac_1 += hRAA->GetBinContent(ip) / (double)Nevt;
      if (int_frac > PDF_lowEdge && int_frac_1 < PDF_lowEdge) // linear interpolation
      {
         double x1 = hRAA->GetBinCenter(ip);
         double x2 = hRAA->GetBinCenter(ip + 1);
         mu_low = x1 + (PDF_lowEdge - int_frac_1) / (int_frac - int_frac_1) * (x2 - x1);
         cout << " Found low bound = " << mu_low;
      }
      if (int_frac > PDF_median && int_frac_1 < PDF_median) // linear interpolation
      {
         double x1 = hRAA->GetBinCenter(ip);
         double x2 = hRAA->GetBinCenter(ip + 1);
         mu = x1 + (PDF_median - int_frac_1) / (int_frac - int_frac_1) * (x2 - x1);
         cout << "\t median = " << mu;
      }
      if (int_frac > PDF_highEdge && int_frac_1 < PDF_highEdge) // linear interpolation
      {
         double x1 = hRAA->GetBinCenter(ip);
         double x2 = hRAA->GetBinCenter(ip + 1);
         mu_high = x1 + (PDF_highEdge - int_frac_1) / (int_frac - int_frac_1) * (x2 - x1);
         cout << "\t high bound = " << mu_high << endl;
      }
   }

   raa1 = mu;
   raa2 = mu_gaus;
   lowEdge = mu_low;
   highEdge = mu_high;

   delete hRAA;
   delete hAA;
   delete hPP;
   delete gRandom;
}


void write_RAA()
{
   globalSetting();

   //define Nbin

   //read D0 yield in AuAu
   TGraphErrors* gD0_Run14HFT_err[ncent];
   TGraphErrors* gD0_Run14HFT_sys[ncent];
    TFile* fin = new TFile("../../ptShift/D0_Spectra_Run14HFT.root");
   for (int icent = 0; icent < ncent; icent++)
   {
      gD0_Run14HFT_err[icent] = (TGraphErrors*)fin->Get(Form("gD0_err_%s", nameCent1[icent]));
      gD0_Run14HFT_sys[icent] = (TGraphErrors*)fin->Get(Form("gD0_sys_%s", nameCent1[icent]));
   }
   fin->Close();


   //read D0 yield in pp
   // TFile* fin1 = new TFile("data/B_and_D0_ptSpectra.root");
   // // TF1* fppbase = (TF1*)fin1->Get("fD0_pp_Run09");
   // // TF1* fppbase_lw = (TF1*)fin1->Get("fD0_pp_lw_Run09");
   // // TF1* fppbase_up = (TF1*)fin1->Get("fD0_pp_up_Run09");
   // TF1* fppbase = (TF1*)fin1->Get("fD0_pp_Run12");
   // TF1* fppbase_lw = (TF1*)fin1->Get("fD0_pp_lw_Run12");
   // TF1* fppbase_up = (TF1*)fin1->Get("fD0_pp_up_Run12");
   // fin1->Close();

   //read D0 yield in pp
   // TFile* fin1 = new TFile("data/processed_all_yifei.root");
   // TFile* fin1 = new TFile("ppsys.root");
   TFile* fin1 = new TFile("out_ppsys.root");
   TF1* fppbase = (TF1*)fin1->Get("Levynew_pp");// Yifei's fitting as default base line
   TGraphErrors* gr_pp_sys = (TGraphErrors*)fin1->Get("gr_pp_sys");//taking of fitting sys and BW fit difference as sys
   TGraphAsymmErrors * gr_pp_sysAsy= (TGraphAsymmErrors*)fin1->Get("gr_pp_sysAsy");//taking of fitting sys and BW fit difference as sys
   fin1->Close();

   //calculate pp baseline
   float ppbase[npt], ppbase_lw[npt], ppbase_up[npt];
   // const float factor = 0.565 / 42.;
   const float factor = 0.61 / 42.;
   for (int i = 0; i < npt; i++)
   {
      float pt = gD0_Run14HFT_err[0]->GetX()[i];
      ppbase[i] = fppbase->Eval(pt) * factor;

      int NCL = gr_pp_sys->GetN();
      int tmpN = floor(pt / (10. / NCL));
      // cout << "pt= " << pt << " , iPoint = " << tmpN << endl;
      // float ppbaseErr_lw = fabs(gr_pp_sysAsy->GetEYlow()[tmpN] * factor);
      float ppbaseErr_lw = fabs(gr_pp_sys->GetEY()[tmpN] * factor);
      // ppbaseErr_lw = sqrt(pow(ppbaseErr_lw, 2) - pow(0.06 * ppbase[i], 2));
      ppbaseErr_lw = sqrt(pow(ppbaseErr_lw, 2) );
      ppbase_lw[i] = ppbase[i] - ppbaseErr_lw; //cancel TPC track sys error

      // float ppbaseErr_up =  fabs(gr_pp_sysAsy->GetEYhigh()[tmpN] * factor);
      float ppbaseErr_up =  fabs(gr_pp_sys->GetEY()[tmpN] * factor);
      // ppbaseErr_up = sqrt(pow(ppbaseErr_up, 2) - pow(0.06 * ppbase[i], 2));
      ppbaseErr_up = sqrt(pow(ppbaseErr_up, 2) );
      ppbase_up[i] = ppbase[i] + ppbaseErr_up; //cancel TPC track sys error
   }

   // cout << "IS HERE OK ? "<< endl;
   // return;
   //caculate pp baseline error band
   float ppErr_up[npt], ppErr_lw[npt];
   for (int i = 0; i < npt; i++)
   {

      ppErr_lw[i] = (ppbase[i] - ppbase_lw[i]);// /ppbase[i];
      ppErr_up[i] = (ppbase_up[i] - ppbase[i]);// /ppbase[i];
     //0 0.5 13.2409 +/- 5.17613 (39.09%)  +/- 7.97007 (60.19%) 9.50339 (71.77%)
     if(i==0) { ppErr_up[i] = ppErr_lw[i] = ppbase[i]*0.7177;}
      cout << " check here for pp baseline : "<< ppbase[i] << "\t" <<ppbase_lw[i] << "\t" << ppbase_up[i] << endl;
   }

   //RAA and its errors
   float pt_mean[ncent][npt], y[ncent][npt], yerr[ncent][npt], ysys[ncent][npt];
   float raa[ncent][npt], raaErr[ncent][npt], raaSys[ncent][npt]; //ratio err
   float raamedian[ncent][npt], raaSysLow[ncent][npt], raaSysHigh[ncent][npt]; //sys combine
   float raa_pp_up[ncent][npt], raa_pp_lw[ncent][npt]; //pp err
   ofstream out("D0RAA_Run14HFT.txt");
   for (int icent = 0; icent < ncent; icent++)
   {
      out << nameCent[icent] << endl;
      out << "pT \t RAA \t Statistical error \t Systematic error \t pp error (low) \t pp error (up)" << endl;
      for (int i = 0; i < npt; i++)
      {
         pt_mean[icent][i] = gD0_Run14HFT_err[icent]->GetX()[i];
         y[icent][i] = gD0_Run14HFT_err[icent]->GetY()[i];
         yerr[icent][i] = gD0_Run14HFT_err[icent]->GetEY()[i];
         ysys[icent][i] = sqrt(pow(gD0_Run14HFT_sys[icent]->GetEY()[i], 2)); //cancel TPC track sys error
         //cout << y[icent][i] << "\t" << ppbase[i] << endl;
         raa[icent][i] = y[icent][i] / ppbase[i] / NbinMean[icent];
         raaErr[icent][i] = raa[icent][i] * yerr[icent][i] / y[icent][i];
         raaSys[icent][i] = raa[icent][i] * ysys[icent][i] / y[icent][i];
         //
         // // raa_pp_lw[icent][i] = raa[icent][i] * ppErr_lw[i]/ppbase[i]; //pp is in denominator
         // // raa_pp_up[icent][i] = raa[icent][i] * ppErr_up[i]/ppbase[i];
         //
         // // raa_pp_up[icent][i] = y[icent][i]/ (ppbase[i]-ppErr_lw[i]) / NbinMean[icent];
         // raa_pp_lw[icent][i] = raa[icent][i] - y[icent][i] / (ppErr_up[i] + ppbase[i]) / NbinMean[icent];
         // raa_pp_up[icent][i] = y[icent][i] / (ppbase[i] - ppErr_lw[i]) / NbinMean[icent] - raa[icent][i];
         // raa_pp_up[icent][i] = raa_pp_lw[icent][i];
         //
         // method1 calculate invidual raa up and low, add auau quadratic
         float raa_lw =  y[icent][i] / (ppErr_up[i] + ppbase[i]) / NbinMean[icent];
         float raa_up =  y[icent][i] / (ppbase[i] - ppErr_lw[i]) / NbinMean[icent];
         cout << " ================= start new pt bin ========"<<endl;
         cout << " method1 y : ppErr_up : ppbase : ppErr_lw : "<< ppErr_up[i]<< " , " << ppbase[i] <<" , " << ppErr_lw[i]  << endl;
         cout << " central : "<< raa[icent][i] << " , raa_lw : "<< raa_lw << " , raa_up : " << raa_up << endl;

         cout << " ppbase[i] : "<< ppbase[i] << " , ppErr_up[i] : "<< ppErr_up[i] << " , ppErr_lw[i] : "<< ppErr_lw[i] <<endl;
         cout << " raa[icent][i] : "<< raa[icent][i] << " , raa_lw : "<< raa_lw << " , raa_up : " << raa_up <<  " , raaSys[icent][i] : "<< raaSys[icent][i] << endl;
         raa_pp_lw[icent][i] = sqrt(pow(raa[icent][i] - raa_lw, 2) + pow(raaSys[icent][i], 2));
         raa_pp_up[icent][i] = sqrt(pow(raa_up - raa[icent][i], 2) + pow(raaSys[icent][i], 2));
         // raa_pp_lw[icent][i] = raa[icent][i] - sqrt(pow(raa[icent][i] - raa_lw, 2) + pow(raaSys[icent][i], 2));
         // raa_pp_up[icent][i] = raa[icent][i] + sqrt(pow(raa_up - raa[icent][i], 2) + pow(raaSys[icent][i], 2));


         cout << "input  \t"<< y[icent][i] << "\t" << ysys[icent][i] << "\t"<< ppbase[i]*NbinMean[icent] << "\t"<< ppErr_up[i]*NbinMean[icent]<<endl;
         //method2 sampling 
         double tmpraa1, tmpraa2, tmpraa_lw, tmpraa_up;
         find_Assym_Band(tmpraa1, tmpraa2, tmpraa_lw, tmpraa_up, y[icent][i], ysys[icent][i], ppbase[i]*NbinMean[icent], ppErr_up[i]*NbinMean[icent]);
         // raamedian[icent][i] = tmpraa1;
         raamedian[icent][i] = tmpraa2;
         raaSysLow[icent][i] = tmpraa2 - tmpraa_lw;
         raaSysHigh[icent][i] = tmpraa_up - tmpraa2;
         if(i==0) { raaSysLow[icent][i] = raa_pp_lw[icent][i]; raaSysHigh[icent][i] = raa_pp_up[icent][i]; };

         cout << "tmpraa1 : " << tmpraa1 << " , tmpraa2 : " << tmpraa2 << " , tmpraa_lw : " << tmpraa_lw<< " , tmpraa_up :" << tmpraa_up<< endl;
         cout << "tmpraa1 : " << tmpraa1 << " , tmpraa2 : " << tmpraa2 << " , tmpraa_lw : " << 1- tmpraa_lw/tmpraa2<< " , tmpraa_up :" << tmpraa_up/tmpraa2 -1 << endl;
         // if(i<2 || i>=npt -2)
         // if(i<2)
         // {
         //   raaSysLow[icent][i] = sqrt(pow(raa_pp_lw[icent][i],2)+pow(raaSys[icent][i],2));
         //   raaSysHigh[icent][i] = sqrt(pow(raa_pp_up[icent][i],2)+pow(raaSys[icent][i],2));
         // }
         //
         // ------- specific lowest and higest pt bin
         // if(i<2)
         // {
         //   find_Assym_Band(tmpraa1, tmpraa2, tmpraa_lw, tmpraa_up, y[icent][i], ysys[icent][i], ppbase[i]*NbinMean[icent], ppErr_lw[i]*NbinMean[icent]);
         //   raaSysHigh[icent][i] = tmpraa_up - tmpraa1;
         // }
         //
         // if(i>npt-2)
         // {
         //   find_Assym_Band(tmpraa1, tmpraa2, tmpraa_lw, tmpraa_up, y[icent][i], ysys[icent][i], ppbase[i]*NbinMean[icent], ppErr_lw[i]*NbinMean[icent]);
         //   raaSysHigh[icent][i] = tmpraa_up - tmpraa1;
         // }

         cout << "icent = " << icent << ": " << Form("%s",nameCent[icent]) << " ; ipt =" << i << endl;

         cout << endl <<  " raa : " << raa[icent][i] << " , auauErr : " << raaErr[icent][i] << " , auauSys: " << raaSys[icent][i] << endl;
         cout << "ppbase : " << ppbase[i] << " , ppbase_lw : " << ppbase_lw[i]/ppbase[i] << " ppbase_up :" << ppbase_up[i]/ppbase[i] << endl;
         cout << "raa_pp_lw : " << raa_pp_lw[icent][i] << " , raa_pp_up : " << raa_pp_up[icent][i] << endl;

         // cout << "auau : " << y[icent][i] / NbinMean[icent] << " , auau err: " <<  ysys[icent][i] / NbinMean[icent] << " , pp :" << ppbase[i] << " , pp er : " << ppErr_up[i] << endl;
         // cout << "tmpraa1 : " << tmpraa1 << " , tmpraa2 : " << tmpraa2 << " , tmpraa_lw : " << tmpraa_lw << " , tmpraa_up :" << tmpraa_up << endl;
         cout << "tmpraa1 : " << tmpraa1 << " , tmpraa2 : " << tmpraa2 << " , raaSysLow[icent][i] : " << raaSysLow[icent][i]<< " , raaSysHigh[icent][i] :" << raaSysHigh[icent][i]<< endl;

         out << pt_mean[icent][i] << "\t" << raa[icent][i] << "\t" << raaErr[icent][i] << "\t" << raaSys[icent][i] << "\t" << raa_pp_lw[icent][i] << "\t" << raa_pp_up[icent][i] << endl;
      }
      out << endl;
   }
   out.close();
   TGraphErrors* gRAA[ncent];
   TGraphErrors* gRAA_sys[ncent];
   TGraphAsymmErrors* gRAA_sys_combine[ncent];
   TGraphAsymmErrors* gRAA_pp[ncent];
   for (int icent = 0; icent < ncent; icent++)
   {
      gRAA[icent] = new TGraphErrors(npt, pt_mean[icent], raa[icent], 0, raaErr[icent]);
      gRAA_sys[icent] = new TGraphErrors(npt, pt_mean[icent], raa[icent], 0, raaSys[icent]);
      gRAA_pp[icent] = new TGraphAsymmErrors(npt, pt_mean[icent], raa[icent], 0, 0, raa_pp_lw[icent], raa_pp_up[icent]);
      gRAA_sys_combine[icent] = new TGraphAsymmErrors(npt, pt_mean[icent], raamedian[icent], 0, 0, raaSysLow[icent], raaSysHigh[icent]);
   }

   //Write
   TFile* fout = new TFile("D0RAA_Run14HFT.root", "RECREATE");
   fout->cd();
   for (int icent = 0; icent < ncent; icent++)
   {
      gRAA[icent]->SetMarkerStyle(MARKERSTYLE[0]);
      gRAA[icent]->SetMarkerColor(COLOR[0]);
      gRAA[icent]->SetMarkerSize(2.5);
      gRAA[icent]->Write(Form("D0_RAA_err_%s", nameCent1[icent]));
      gRAA_sys[icent]->SetMarkerStyle(MARKERSTYLE[0]);
      gRAA_sys[icent]->SetMarkerColor(COLOR[0]);
      gRAA_sys[icent]->SetMarkerSize(2.5);
      gRAA_sys[icent]->Write(Form("D0_RAA_sys_%s", nameCent1[icent]));
      gRAA_pp[icent]->SetMarkerStyle(MARKERSTYLE[0]);
      gRAA_pp[icent]->SetMarkerColor(COLOR[0]);
      gRAA_pp[icent]->SetMarkerSize(2.5);
      gRAA_pp[icent]->Write(Form("D0_RAA_pperr_%s", nameCent1[icent]));
      gRAA_sys_combine[icent]->SetMarkerStyle(MARKERSTYLE[0]);
      gRAA_sys_combine[icent]->SetMarkerColor(COLOR[0]);
      gRAA_sys_combine[icent]->SetMarkerSize(2.5);
      gRAA_sys_combine[icent]->Write(Form("D0_RAA_sys_combine_%s", nameCent1[icent]));
   }
   fout->Close();
}
