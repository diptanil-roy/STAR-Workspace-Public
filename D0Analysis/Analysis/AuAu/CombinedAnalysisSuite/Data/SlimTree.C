Double_t GetDeltaPhi(Double_t phia, Double_t phib)
{
    Double_t pi = TMath::Pi();

    if (phia < 0)         phia += 2*pi;
    else if (phia > 2*pi) phia -= 2*pi;
    if (phib < 0)         phib += 2*pi;
    else if (phib > 2*pi) phib -= 2*pi;
    Double_t dphi = phib - phia;

    if (dphi < -1.0*pi)  dphi += 2*pi;
    else if (dphi > pi)  dphi -= 2*pi;

    // phi is between  -pi < phi < pi                                                                                                                                                                                                                                       
    return dphi;
}

Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

const Int_t PTBINS = 5;
Double_t edges[PTBINS+1] = {1.0, 2.0, 3.0, 4.0, 5.0, 10.0};

TH1D *hD0Pt = new TH1D("hD0Pt", "hD0Pt", PTBINS, edges);

TFile *eff = new TFile("EssentialFiles/D0Eff.root", "READ");
TH1F *hEff0 = (TH1F *)eff->Get("hEff0");
TH1F *hEff1 = (TH1F *)eff->Get("hEff1");
TH1F *hEff2 = (TH1F *)eff->Get("hEff2");
TH1F *hEff3 = (TH1F *)eff->Get("hEff3");
TH1F *hEff4 = (TH1F *)eff->Get("hEff4");

TH1F *hEff0_Sys = (TH1F *)hEff0->Clone("hEff0_Sys");
TH1F *hEff1_Sys = (TH1F *)hEff1->Clone("hEff1_Sys");
TH1F *hEff2_Sys = (TH1F *)hEff2->Clone("hEff2_Sys");
TH1F *hEff3_Sys = (TH1F *)hEff3->Clone("hEff3_Sys");
TH1F *hEff4_Sys = (TH1F *)hEff4->Clone("hEff4_Sys");

TFile *Eff_Sys = new TFile("EssentialFiles/Sys4Matt.root", "READ");
TGraph *gSys_cent0_err0 = (TGraph *)Eff_Sys->Get("gSys_cent0_err0");
TGraph *gSys_cent0_err1 = (TGraph *)Eff_Sys->Get("gSys_cent0_err1");
TGraph *gSys_cent0_err2 = (TGraph *)Eff_Sys->Get("gSys_cent0_err2");
TGraph *gSys_cent0_err3 = (TGraph *)Eff_Sys->Get("gSys_cent0_err3");
TGraph *gSys_cent0_err4 = (TGraph *)Eff_Sys->Get("gSys_cent0_err4");
TGraph *gSys_cent0_err5 = (TGraph *)Eff_Sys->Get("gSys_cent0_err5");
TGraph *gSys_cent0_err6 = (TGraph *)Eff_Sys->Get("gSys_cent0_err6");
TGraph *gSys_cent0_err7 = (TGraph *)Eff_Sys->Get("gSys_cent0_err7");
TGraph *gSys_cent0_err8 = (TGraph *)Eff_Sys->Get("gSys_cent0_err8");
TGraph *gSys_cent0_err9 = (TGraph *)Eff_Sys->Get("gSys_cent0_err9");

TGraph *gSys_cent5_err0 = (TGraph *)Eff_Sys->Get("gSys_cent5_err0");
TGraph *gSys_cent5_err1 = (TGraph *)Eff_Sys->Get("gSys_cent5_err1");
TGraph *gSys_cent5_err2 = (TGraph *)Eff_Sys->Get("gSys_cent5_err2");
TGraph *gSys_cent5_err3 = (TGraph *)Eff_Sys->Get("gSys_cent5_err3");
TGraph *gSys_cent5_err4 = (TGraph *)Eff_Sys->Get("gSys_cent5_err4");
TGraph *gSys_cent5_err5 = (TGraph *)Eff_Sys->Get("gSys_cent5_err5");
TGraph *gSys_cent5_err6 = (TGraph *)Eff_Sys->Get("gSys_cent5_err6");
TGraph *gSys_cent5_err7 = (TGraph *)Eff_Sys->Get("gSys_cent5_err7");
TGraph *gSys_cent5_err8 = (TGraph *)Eff_Sys->Get("gSys_cent5_err8");
TGraph *gSys_cent5_err9 = (TGraph *)Eff_Sys->Get("gSys_cent5_err9");

TGraph *gSys_cent6_err0 = (TGraph *)Eff_Sys->Get("gSys_cent6_err0");
TGraph *gSys_cent6_err1 = (TGraph *)Eff_Sys->Get("gSys_cent6_err1");
TGraph *gSys_cent6_err2 = (TGraph *)Eff_Sys->Get("gSys_cent6_err2");
TGraph *gSys_cent6_err3 = (TGraph *)Eff_Sys->Get("gSys_cent6_err3");
TGraph *gSys_cent6_err4 = (TGraph *)Eff_Sys->Get("gSys_cent6_err4");
TGraph *gSys_cent6_err5 = (TGraph *)Eff_Sys->Get("gSys_cent6_err5");
TGraph *gSys_cent6_err6 = (TGraph *)Eff_Sys->Get("gSys_cent6_err6");
TGraph *gSys_cent6_err7 = (TGraph *)Eff_Sys->Get("gSys_cent6_err7");
TGraph *gSys_cent6_err8 = (TGraph *)Eff_Sys->Get("gSys_cent6_err8");
TGraph *gSys_cent6_err9 = (TGraph *)Eff_Sys->Get("gSys_cent6_err9");

TH1F *D0Sys_010 = new TH1F("D0Sys_010", "D0Sys_010", PTBINS, edges);
TH1F *D0Sys_1040 = new TH1F("D0Sys_1040", "D0Sys_1040", PTBINS, edges);
TH1F *D0Sys_4080 = new TH1F("D0Sys_480", "D0Sys_4080", PTBINS, edges);

TFile *doub = new TFile("EssentialFiles/MisPID_SB_Final.root", "READ");
TGraph *gDoub0 = (TGraph *)doub->Get("DoubleCounting_Cen_0_SB");
TGraph *gDoub1 = (TGraph *)doub->Get("DoubleCounting_Cen_1_SB");
TGraph *gDoub2 = (TGraph *)doub->Get("DoubleCounting_Cen_2_SB");
TGraph *gDoub3 = (TGraph *)doub->Get("DoubleCounting_Cen_3_SB");

double getEff(int mCen, float mPt)
{
    double eff = 0;

    // Note I am using the efficiency histograms with "_Sys" here. If running applyWeights() only, this doesn't matter
    // If running doSys(), this will shuffle the efficinecy with a Gaussian each iteration in the _Sys histograms.
    // To not have to change the code here, _Sys is used all the time but is equivelent to original histograms in nominal running
    if (mCen >= 0 && mCen < 10)
        eff = hEff0_Sys->GetBinContent(hEff0_Sys->FindBin(mPt));
    else if (mCen >= 10 && mCen < 20)
        eff = hEff1_Sys->GetBinContent(hEff1_Sys->FindBin(mPt));
    else if (mCen >= 20 && mCen < 40)
        eff = hEff2_Sys->GetBinContent(hEff2_Sys->FindBin(mPt));
    else if (mCen >= 40 && mCen < 60)
        eff = hEff3_Sys->GetBinContent(hEff3_Sys->FindBin(mPt));
    else if (mCen >= 60 && mCen < 80)
        eff = (2. / 3.) * hEff4_Sys->GetBinContent(hEff4_Sys->FindBin(mPt)); // efficiency is from paper which has scale factor
    return eff;
}
double getDoubleCount(int mCen, float mPt)
{
    double dcount = 1;
    if (mCen == 0 && mCen == 80)
        dcount = gDoub0->Eval(mPt);
    else if (mCen >= 0 && mCen < 10)
        dcount = gDoub1->Eval(mPt);
    else if (mCen >= 10 && mCen < 40)
        dcount = gDoub2->Eval(mPt);
    else if (mCen >= 40 && mCen < 80)
        dcount = gDoub3->Eval(mPt);
    return dcount;
}

void Method(TString mode = "US", TH1D *hD0Mass[PTBINS] = NULL, int cenlow = 0, int cenhigh = 80){
    // TFile *f_D = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/D0Jets_WithCS_JetTree_May20.root");
    // TChain *t1 = (TChain*)f_D->Get("JetTree/D0Jets");//JetTreeSigBg/SignalPlusBackground");

  TFile *f_D;
  f_D = new TFile("/Volumes/WorkDrive/STAR-Workspace/D0Analysis/Analysis/AuAu/Files/NewDataSample.root");


  TChain *t1;
  if (mode == "US")t1 = (TChain*)f_D->Get("JetTree_Standard_UnlikeSign/D0Jets");//JetTreeSigBg/SignalPlusBackground");
  else if (mode == "LS")t1 = (TChain*)f_D->Get("JetTree_Standard_LikeSign/D0Jets");//JetTreeSigBg/SignalPlusBackground");

  double binning[6] = {0,0.05,0.1,0.2,0.3,0.4};

  float JetArea;
  float D0Mass;float D0Pt;float D0Phi; float D0Eta;
  float KaonPt;float KaonPhi; float KaonEta;float PionPt;float PionPhi; float PionEta;
  float JetCorrPt; float JetCSEta; float JetEta; float JetCSPhi; float JetPhi;float JetPt;
  float JetHighestTrackPt;
  float Centrality;
  float Weight;
  float TrackPt[10000];
  float TrackEta[10000];
  float TrackPhi[10000];
  float TrackID[10000];
  float TrackPx[10000];
  float TrackPy[10000];
  float TrackPz[10000];
  int JetNConst;
  float _mM; float _mPt; float _mEta; float _mJPt; float _mJPtArea; float _mR; float _mZ;float _mHPt;float _mJEta; float _mZ_Area; float _m_gRefMultCorr;
  int _mCen;float _mWeight;float _mJetArea;int RunID;int EventId; float jetptfromarea; int jetnconstfromarea;
  int _mConstArea; int _mConstCS;

  t1->SetBranchStatus("*",0);
  t1->SetBranchStatus("D0Mass",1);
  t1->SetBranchStatus("KaonPt",1);
  t1->SetBranchStatus("KaonPhi",1);
  t1->SetBranchStatus("KaonEta",1);
  t1->SetBranchStatus("PionPt",1);
  t1->SetBranchStatus("PionPhi",1);
  t1->SetBranchStatus("PionEta",1);
  t1->SetBranchStatus("Centrality",1);
  t1->SetBranchStatus("gRefMult", 1);
  t1->SetBranchStatus("Weight",1);
  t1->SetBranchStatus("RunID",1);
  t1->SetBranchStatus("EventId",1);
  t1->SetBranchStatus("JetNConst",1);
  t1->SetBranchStatus("JetEta",1);
  t1->SetBranchStatus("TrackID", 1);
  t1->SetBranchStatus("TrackPt", 1);
  t1->SetBranchStatus("TrackEta", 1);
  t1->SetBranchStatus("TrackPhi", 1);
  t1->SetBranchStatus("TrackPx", 1);
  t1->SetBranchStatus("TrackPy", 1);
  t1->SetBranchStatus("TrackPz", 1);
  t1->SetBranchStatus("TrackCharge", 1);
  
  t1->SetBranchAddress( "D0Mass" , &D0Mass );
  t1->SetBranchAddress( "KaonPt" , &KaonPt );
  t1->SetBranchAddress( "KaonPhi" , &KaonPhi );
  t1->SetBranchAddress( "KaonEta" , &KaonEta );
  t1->SetBranchAddress( "PionPt" , &PionPt );
  t1->SetBranchAddress( "PionPhi" , &PionPhi );
  t1->SetBranchAddress( "PionEta" , &PionEta );
  t1->SetBranchAddress( "Centrality" , &Centrality );
  t1->SetBranchAddress( "gRefMult" , &_m_gRefMultCorr );
  t1->SetBranchAddress( "Weight" , &Weight );
  t1->SetBranchAddress( "RunID" , &RunID );
  t1->SetBranchAddress( "EventId" , &EventId );
  t1->SetBranchAddress( "JetNConst" , &JetNConst );
  t1->SetBranchAddress( "JetEta" , &JetEta );
  t1->SetBranchAddress("TrackID", TrackID);
  t1->SetBranchAddress("TrackPt", TrackPt);
  t1->SetBranchAddress("TrackEta", TrackEta);
  t1->SetBranchAddress("TrackPhi", TrackPhi);
  t1->SetBranchAddress("TrackPx", TrackPx);
  t1->SetBranchAddress("TrackPy", TrackPy);
  t1->SetBranchAddress("TrackPz", TrackPz);
//   t1->SetBranchAddress("TrackCharge", TrackCharge);

  int loop = t1->GetEntries();//   

  int TotalEventsEncountered = 0;

  for(int i =0;i<loop;i++){
      if(i%10000==0)cout << "On "<< i << " out of " << loop << " " << float(i)/loop*100 << "%" << endl;
      t1->GetEntry(i);
       
      int d0index = -99;
      if (JetNConst == 0) continue;

    for (int itrk = 0; itrk < JetNConst; itrk++){
        if (TrackID[itrk] == 421) {
            d0index = itrk;
            break;
        }
    }

    if (D0Mass < 1.75 || D0Mass > 2.02) continue;
    if (abs(JetEta) > 0.6) continue;
    if (PionPt <= 0.6 || KaonPt <= 0.6) continue;
    //   if(D0Mass<1.75 || D0Mass>2.02)continue;
      if(PionPt<0.6)continue;
      if(KaonPt<0.6)continue;

      TVector3 p1( KaonPt * cos(KaonPhi),KaonPt * sin(KaonPhi),KaonPt * sinh(KaonEta));
      TVector3 p2( PionPt * cos(PionPhi),PionPt * sin(PionPhi),PionPt * sinh(PionEta));
      TVector3 p3;

      p3 = p1 + p2;
      D0Pt = p3.Pt();
      D0Eta = p3.Eta();
      D0Phi = p3.Phi();

      if (D0Pt < 1.0 || D0Pt > 10.0) continue;

    //   cout << D0Pt << "\t" << TrackPt[d0index] << endl;
     
      int ccen = -1;
      if(Centrality == 0)ccen = 70;
      if(Centrality == 1)ccen = 60;
      if(Centrality == 2)ccen = 50;
      if(Centrality == 3)ccen = 40;
      if(Centrality == 4)ccen = 30;
      if(Centrality == 5)ccen = 20;
      if(Centrality == 6)ccen = 10;
      if(Centrality == 7)ccen = 5;
      if(Centrality == 8)ccen = 0;
      if (ccen < cenlow || ccen >= cenhigh) continue;

      TotalEventsEncountered++;
      
      double eff = getEff(ccen, D0Pt);
      double dcount = getDoubleCount(ccen, D0Pt);
      int ptbin = hD0Pt->FindBin(D0Pt);
      for (int ptbin = 0; ptbin < PTBINS; ptbin++){
        if (D0Pt >= edges[ptbin]){
            hD0Mass[ptbin]->Fill(D0Mass, (1.-dcount)*Weight/eff);
        }
      }

    //   hD0Mass[ptbin-1]->Fill(D0Mass, (1.-dcount)*Weight/eff);
  }
  cout << "Total Events Encountered in " << mode.Data() << "\t" << TotalEventsEncountered << endl;
//   return hD0Mass;
}

void SlimTree(){
    TFitter::SetPrecision(0.1);

    TH1D *US[4][PTBINS];
    TH1D *LS[4][PTBINS];

    for (int cent = 0; cent < 4; cent++){
        for (int i = 0; i < PTBINS; i++){
            US[cent][i] = new TH1D(Form("US_%s_%i_%i_%i", "US", cent, (int)edges[i], 10), Form("US_%s_%i_%i_%i", "US", cent, (int)edges[i], 10), 40, 1.75, 2.02);
            LS[cent][i] = new TH1D(Form("LS_%s_%i_%i_%i", "LS", cent, (int)edges[i], 10), Form("LS_%s_%i_%i_%i", "LS", cent, (int)edges[i], 10), 40, 1.75, 2.02);
        }
    }

    Method("US", US[0], 0, 80);
    Method("LS", LS[0], 0, 80);

    Method("US", US[1], 0, 10);
    Method("LS", LS[1], 0, 10);

    Method("US", US[2], 10, 40);
    Method("LS", LS[2], 10, 40);

    Method("US", US[3], 40, 80);
    Method("LS", LS[3], 40, 80);
    
    TF1 *background[4][PTBINS];
    TF1 *fitfunc[4][PTBINS];
    TF1 *backgroundfit[4][PTBINS];

    TCanvas *c[4];

    TH1D *D0PtFromFit[4];
    for (int cent = 0; cent < 4; cent++){
        D0PtFromFit[cent] = new TH1D(Form("D0PtFromFit_%i", cent), Form("D0PtFromFit_%i", cent), PTBINS, edges);
    }

    double dx = US[0][0]->GetBinWidth(1);

    for (int cent = 0; cent < 4; cent++){
        cout << "Centrality Bin = " << cent << endl;
        for (int i = 0; i < PTBINS; i++){
            background[cent][i] = new TF1(Form("background_%i_%i", cent, i), "[0]+[1]*x+[2]*x**2+[3]*x**3", 1.7, 2.10);
            background[cent][i]->SetLineColor(1);
            fitfunc[cent][i] = new TF1(Form("fitfunc_%i_%i", cent, i), "gausn + pol3(3)", 1.7, 2.10);
            fitfunc[cent][i]->SetLineColor(6);
            fitfunc[cent][i]->SetLineWidth(5);

            
            LS[cent][i]->Fit(Form("background_%i_%i", cent, i), "RLQ0");

            backgroundfit[cent][i] = (TF1 *)LS[cent][i]->GetFunction(Form("background_%i_%i", cent, i));

            fitfunc[cent][i]->SetParLimits(0, 0, 200000000);

            fitfunc[cent][i]->SetParameter(1, 1.865);
            fitfunc[cent][i]->SetParLimits(1, 1.86, 1.87);
            fitfunc[cent][i]->SetParLimits(2, 0.01, 0.03);

            fitfunc[cent][i]->SetParameter(3, backgroundfit[cent][i]->GetParameter(0));
            fitfunc[cent][i]->SetParameter(4, backgroundfit[cent][i]->GetParameter(1));
            fitfunc[cent][i]->SetParameter(5, backgroundfit[cent][i]->GetParameter(2));
            fitfunc[cent][i]->SetParameter(6, backgroundfit[cent][i]->GetParameter(3));

            // cout << "STARTING FIT" << endl;

            int status = US[cent][i]->Fit(Form("fitfunc_%i_%i", cent, i), "RLQ0");

            double S_whole = fitfunc[cent][i]->GetParameter(0)/dx;
            double S_whole_err = fitfunc[cent][i]->GetParError(0)/dx;

            // D0PtFromFit[cent]->SetBinContent(i+1, S_whole);
            // D0PtFromFit[cent]->SetBinError(i+1, S_whole_err);

            // cout << "STATUS = " << status << endl;
            cout << edges[i] << " < pT,D0 [GeV/c] < " << 10 << " =========> " << S_whole << "\t" << "pm " << S_whole_err << endl;
        }

        c[cent] = new TCanvas(Form("D0Mass_%i", cent), Form("D0Mass_%i", cent), 2100, 700);
        c[cent]->Divide(PTBINS);

        for (int i = 0; i < PTBINS; i++){
            c[cent]->cd(i+1);
            // US[cent][i]->GetYaxis()->SetRangeUser(0, 1000000);
            US[cent][i]->Draw("SAME");
            fitfunc[cent][i]->Draw("SAME");
            LS[cent][i]->SetLineColor(kRed);
            LS[cent][i]->Draw("SAME");
            backgroundfit[cent][i]->Draw("SAME");
        }
    }

    // TCanvas *d = new TCanvas("D0 pT From Fit", "D0 pT From Fit", 2100, 700);
    // d->cd();
    // for (int cent = 0; cent < 4; cent++){
    //     D0PtFromFit[cent]->SetLineColor(cent+1);
    //     D0PtFromFit[cent]->Draw("SAME");
    // }

    TH1D *D0PtFromSPlot[4];

    TFile *SPlotFits = new TFile("Aug17_2023/Histograms_D01_10GeV_RecoJetPt_1_1000.root");
    SPlotFits->cd();
    D0PtFromSPlot[0] = (TH1D *)gDirectory->Get("D0Pt");
    D0PtFromSPlot[1] = (TH1D *)gDirectory->Get("D0Pt_0_10");
    D0PtFromSPlot[2] = (TH1D *)gDirectory->Get("D0Pt_10_40");
    D0PtFromSPlot[3] = (TH1D *)gDirectory->Get("D0Pt_40_80");

    for (int cent = 0; cent < 4; cent++){
        for (int bin = 0; bin < PTBINS; bin++){
            double splot = D0PtFromSPlot[cent]->Integral(D0PtFromSPlot[cent]->FindBin(edges[bin]), D0PtFromSPlot[cent]->FindBin(10));
            double fit = fitfunc[cent][bin]->GetParameter(0)/dx;
            double deviation = (fit - splot)*100./splot;

            // cout << "Centrality " << cent << "\t" << edges[bin] << " < pT,D0 [GeV/c] < " << 10 << "\t" << splot << "\t" << fit << "\t" << deviation << " % " << endl;
            // cout << Form("Centrality %i \t %i < pT,D0 < 10 < %.2f %", cent, (int)edges[bin], deviation) << endl;
            cout << Form("%.2f, \t", abs(deviation));
        }
        cout << endl;
    }

    // TCanvas *e = new TCanvas("D0 pT Compared", "D0 pT Compared", 2100, 700);
    // e->Divide(4);
    // for (int cent = 0; cent < 4; cent++){
    //     e->cd(cent+1);
    //     gPad->SetLogy();
    //     D0PtFromFit[cent]->SetLineColor(cent+1);
    //     D0PtFromSPlot[cent]->SetLineColor(cent+1);
    //     D0PtFromFit[cent]->SetMarkerColor(cent+1);
    //     D0PtFromSPlot[cent]->SetMarkerColor(cent+1);
    //     D0PtFromFit[cent]->SetMarkerStyle(20);
    //     D0PtFromSPlot[cent]->SetMarkerStyle(24);

    //     D0PtFromFit[cent]->Draw("SAME");
    //     D0PtFromSPlot[cent]->Draw("SAME");
    // }

    // TH1D *Variation[4];

    // for (int cent = 0; cent < 4; cent++){
    //     Variation[cent] = (TH1D *)D0PtFromFit[cent]->Clone(Form("Variation_%i", cent));
    //     Variation[cent]->Reset();
    //     Variation[cent]->SetNameTitle(Form("Variation_%i", cent), Form("Variation_%i", cent));
    //     for (int i = 1; i <= Variation[cent]->GetNbinsX(); i++){
    //         double bincontent = (D0PtFromFit[cent]->GetBinContent(i) - D0PtFromSPlot[cent]->GetBinContent(i))*100/D0PtFromSPlot[cent]->GetBinContent(i);
    //         cout << "Bin " << i << "\t" << "\t" << bincontent << endl;

    //         Variation[cent]->SetBinContent(i, bincontent);
    //     }
    // }

    // TCanvas *d2 = new TCanvas("D0 Yield Variation", "D0 Yield Variation", 2100, 700);
    // d2->Divide(4);    
    // for (int cent = 0; cent < 4; cent++){
    //     d2->cd(cent+1);
    //     Variation[cent]->SetLineColor(cent+1);
    //     Variation[cent]->Draw("HIST");

    //     cout << "Centrality " << cent << "\t" << Variation[cent]->GetMean() << "\t" << Variation[cent]->GetMeanError() << endl;
    // }
}
