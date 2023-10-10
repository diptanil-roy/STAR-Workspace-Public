using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void PlotVsD0pT(TString BaseDirName = "Aug14_FONLL", TString ExtraTag = "WeighedByData", int SUPERITERATION = 0, int iteration = 4){
    TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);


	double taa[3] = {941.23714*1.0318440e+08, 391.35550*3.2123506e+08, 56.62475*4.6679240e+08};

    int color[5] = {kRed, kGreen-2, kBlue, kBlack, kCyan};

    int ptstartpoint[5] = {1, 2, 3, 4, 5};

    TString JetPtXaxisName = "p_{T, Jet}^{D^{0}} [GeV/#it{c}]";
    TString JetZXaxisName = "z_{Jet}^{D^{0}}";
    TString JetdRXaxisName = "#DeltaR";

    TString JetPtYaxisName = "#frac{1}{N_{Evt}} #frac{d^{2}N}{p_{T}dp_{T, Jet}d#eta} [GeV/#it{c}]^{-2}";
    TString JetZYaxisName = "#frac{1}{N_{Evt}} #frac{d^{2}N}{zdz_{Jet}d#eta}";
    TString JetdRYaxisName = "#frac{1}{N_{Jet}} #frac{#DeltaN_{Jet}}{#DeltaR}";

    TString RCPYaxisName = "R_{CP}";
    TString RMPYaxisName = "R_{MP}";

    TString RatiodRYaxisName = Form("Ratio(%s)", JetdRYaxisName.Data());

    TH1D *JetpT[5][3];
    TH1D *JetZ[5][3];
    TH1D *JetdR[5][3];

    TH1D *RCP_Pt[5][3];
    TH1D *RCP_Z[5][3];
    TH1D *Ratio_dR[5][3];

    if (ExtraTag != "") ExtraTag = TString("_") + ExtraTag;

    cout << ExtraTag << endl;

    for (int i = 0; i < 5; i++){ 
        TString DirName = Form("%s_%iGeV_Data%s%i", BaseDirName.Data(), ptstartpoint[i], ExtraTag.Data(), iteration);
        TString FileName = Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins);
        cout << "<=== " << ptstartpoint[i] <<  " < pT [GeV] < 10 from file " << FileName.Data() << " ===>" << endl;

        TFile *f = new TFile(FileName.Data());
        f->cd();

        for (int cent = 0; cent < 3; cent++){
            JetpT[i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded Wide p_{T} Cent = %i SI = %i Iter = %i", cent, 0, iteration));
            // cout << JetpT[i][cent]->GetName() << endl;
            SetName(JetpT[i][cent], Form("%s D0pT = %i GeV", JetpT[i][cent]->GetName(), ptstartpoint[i]));
            SetAxisTitles(JetpT[i][cent], JetPtXaxisName.Data(), JetPtYaxisName.Data());
            JetZ[i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded Wide Z Cent = %i SI = %i Iter = %i", cent, 0, iteration));
            // cout << JetZ[i][cent]->GetName() << endl;
            SetName(JetZ[i][cent], Form("%s D0pT = %i GeV", JetZ[i][cent]->GetName(), ptstartpoint[i]));
            SetAxisTitles(JetZ[i][cent], JetZXaxisName.Data(), JetZYaxisName.Data());
            JetdR[i][cent] = (TH1D *)gDirectory->Get(Form("Unfolded Wide dR Cent = %i SI = %i Iter = %i", cent, 0, iteration));
            // cout << JetZ[i][cent]->GetName() << endl;
            SetName(JetdR[i][cent], Form("%s D0pT = %i GeV", JetdR[i][cent]->GetName(), ptstartpoint[i]));
            SetAxisTitles(JetdR[i][cent], JetdRXaxisName.Data(), JetdRYaxisName.Data());
            JetpT[i][cent]->Scale(1./taa[cent]);
            JetZ[i][cent]->Scale(1./taa[cent]);
            ProcessSpectra(JetpT[i][cent]);
            ProcessSpectra(JetZ[i][cent]);
            JetdR[i][cent]->Scale(1./JetdR[i][cent]->Integral());
            DivideByBinWidth(JetdR[i][cent]);
            // cout << "Done" << endl;
            SetColor(JetpT[i][cent], color[i]);
            SetColor(JetZ[i][cent], color[i]);
            SetColor(JetdR[i][cent], color[i]);

            RCP_Pt[i][cent] = (TH1D *)JetpT[i][cent]->Clone();
            SetName(RCP_Pt[i][cent], Form("%s RCP", JetpT[i][cent]->GetName()));
            RCP_Z[i][cent] = (TH1D *)JetZ[i][cent]->Clone();
            SetName(RCP_Z[i][cent], Form("%s RCP", JetZ[i][cent]->GetName()));
            Ratio_dR[i][cent] = (TH1D *)JetdR[i][cent]->Clone();
            SetName(Ratio_dR[i][cent], Form("%s Ratio", JetdR[i][cent]->GetName()));
        }

        RCP_Pt[i][0]->Divide(RCP_Pt[i][2]);
        RCP_Pt[i][1]->Divide(RCP_Pt[i][2]);
        RCP_Z[i][0]->Divide(RCP_Z[i][2]);
        RCP_Z[i][1]->Divide(RCP_Z[i][2]);
        Ratio_dR[i][0]->Divide(Ratio_dR[i][2]);
        Ratio_dR[i][1]->Divide(Ratio_dR[i][2]);
    }

    TCanvas *pT = new TCanvas("pT", "pT", 1800, 600);
    pT->Divide(3);

    TCanvas *Z = new TCanvas("Z", "Z", 1800, 600);
    Z->Divide(3);

    TCanvas *dR = new TCanvas("dR", "dR", 1800, 600);
    dR->Divide(3);

    TCanvas *RCPpT = new TCanvas("RCPpT", "RCPpT", 1800, 600);
    RCPpT->Divide(2);

    TCanvas *RCPZ = new TCanvas("RCPZ", "RCPZ", 1800, 600);
    RCPZ->Divide(2);

    TCanvas *Ratio_dR_Canvas = new TCanvas("Ratio_dR_Canvas", "Ratio_dR_Canvas", 1800, 600);
    Ratio_dR_Canvas->Divide(2);

    auto legend1 = new TLegend(0.5, 0.6, 0.9, 0.9);
    auto legend2 = new TLegend(0.6, 0.1, 0.9, 0.4);
    auto legend3 = new TLegend(0.5, 0.6, 0.9, 0.9);

    for (int cent = 0; cent < 3; cent++){
        pT->cd(cent+1);
        gPad->SetLeftMargin(0.15);
        for (int i = 0; i < 5; i++){
            JetpT[i][cent]->GetXaxis()->SetRangeUser(5.0, 20.0);
            JetpT[i][cent]->GetYaxis()->SetRangeUser(1e-12, 1e-6);
            
            JetpT[i][cent]->GetXaxis()->SetTitleOffset(1.2);
            JetpT[i][cent]->GetXaxis()->SetTitleSize(0.035);
            JetpT[i][cent]->GetYaxis()->SetTitleOffset(1.8);
            JetpT[i][cent]->GetYaxis()->SetTitleSize(0.035);
            
            JetpT[i][cent]->Draw("SAME 9");
            if (cent == 0) legend1->AddEntry(JetpT[i][cent], Form("%i < p_{T}^{D^{0}} [GeV] < %i", ptstartpoint[i], 10), "lp");
        }
        legend1->SetTextSize(0.035);
        legend1->Draw("same");
        gPad->SetLogy();
        Z->cd(cent+1);
        gPad->SetLeftMargin(0.15);
        for (int i = 0; i < 5; i++){
            JetZ[i][cent]->GetYaxis()->SetRangeUser(1e-10, 1e-4);
            JetZ[i][cent]->GetXaxis()->SetTitleOffset(1.2);
            JetZ[i][cent]->GetYaxis()->SetTitleOffset(2.0);

            JetZ[i][cent]->Draw("same");

            if (cent == 0) legend2->AddEntry(JetZ[i][cent], Form("%i < p_{T}^{D^{0}} [GeV] < %i", ptstartpoint[i], 10), "lp");
        }
        legend2->Draw("same");
        gPad->SetLogy();
        dR->cd(cent+1);
        gPad->SetLeftMargin(0.15);
        for (int i = 0; i < 5; i++){
            JetdR[i][cent]->GetYaxis()->SetRangeUser(1e-6, 1e3);
            JetdR[i][cent]->GetXaxis()->SetTitleOffset(1.2);
            JetdR[i][cent]->GetYaxis()->SetTitleOffset(2.0);

            JetdR[i][cent]->Draw("same");

            if (cent == 0) legend3->AddEntry(JetdR[i][cent], Form("%i < p_{T}^{D^{0}} [GeV] < %i", ptstartpoint[i], 10), "lp");
        }
        legend3->Draw("same");
        gPad->SetLogy();
    }

    auto legend4 = new TLegend(0.1, 0.6, 0.4, 0.9);
    auto legend5 = new TLegend(0.6, 0.6, 0.9, 0.9);
    auto legend6 = new TLegend(0.5, 0.6, 0.9, 0.9);

    for (int centbin = 0; centbin < 2; centbin++){
        RCPpT->cd(centbin+1);
        for (int i = 0; i < 5; i++){
            RCP_Pt[i][centbin]->GetXaxis()->SetRangeUser(5.0, 20.0);
            RCP_Pt[i][centbin]->GetYaxis()->SetRangeUser(0.0, 2.5);
            RCP_Pt[i][centbin]->GetXaxis()->SetTitleOffset(1.2);
            RCP_Pt[i][centbin]->GetYaxis()->SetTitleOffset(1.2);
            if (centbin == 0)SetAxisTitles(RCP_Pt[i][centbin], JetPtXaxisName.Data(), RCPYaxisName.Data());
            else SetAxisTitles(RCP_Pt[i][centbin], JetPtXaxisName.Data(), RMPYaxisName.Data());
            RCP_Pt[i][centbin]->Draw("same");

            if (centbin == 0) legend4->AddEntry(RCP_Pt[i][centbin], Form("%i < p_{T}^{D^{0}} [GeV] < %i", ptstartpoint[i], 10), "lp");
        }
        legend4->SetTextSize(0.035);
        legend4->Draw("same");
        RCPZ->cd(centbin+1);
        for (int i = 0; i < 5; i++){
            RCP_Z[i][centbin]->GetYaxis()->SetRangeUser(-0.5, 2.5);
            RCP_Z[i][centbin]->GetXaxis()->SetTitleOffset(1.2);
            RCP_Z[i][centbin]->GetYaxis()->SetTitleOffset(1.2);
            if (centbin == 0)SetAxisTitles(RCP_Z[i][centbin], JetZXaxisName.Data(), RCPYaxisName.Data());
            else SetAxisTitles(RCP_Z[i][centbin], JetZXaxisName.Data(), RMPYaxisName.Data());
            RCP_Z[i][centbin]->Draw("same");

            if (centbin == 0) legend5->AddEntry(RCP_Z[i][centbin], Form("%i < p_{T}^{D^{0}} [GeV] < %i", ptstartpoint[i], 10), "lp");
        }
        legend5->SetTextSize(0.035);
        legend5->Draw("same");
        Ratio_dR_Canvas->cd(centbin+1);
        for (int i = 0; i < 5; i++){
            Ratio_dR[i][centbin]->GetYaxis()->SetRangeUser(0.0, 2.5);
            Ratio_dR[i][centbin]->GetXaxis()->SetTitleOffset(1.2);
            Ratio_dR[i][centbin]->GetYaxis()->SetTitleOffset(1.2);
            if (centbin == 0)SetAxisTitles(Ratio_dR[i][centbin], JetdRXaxisName.Data(), RatiodRYaxisName.Data());
            else SetAxisTitles(Ratio_dR[i][centbin], JetdRXaxisName.Data(), RatiodRYaxisName.Data());
            Ratio_dR[i][centbin]->Draw("same");

            if (centbin == 0) legend6->AddEntry(Ratio_dR[i][centbin], Form("%i < p_{T}^{D^{0}} [GeV] < %i", ptstartpoint[i], 10), "lp");
        }
        legend6->SetTextSize(0.035);
        legend6->Draw("same");
    }


}