#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TObject.h"
#include "TClonesArray.h"
// #include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include <algorithm>

#include "../plotstyles.h"

using namespace std;

pair <int, int> GetCentBin(int centbin){
    pair <int, int> s;

    if (centbin == 0) { s.first = 1; s.second = 9; } //  0-80%
    if (centbin == 1) { s.first = 8; s.second = 9; } //  0-10%
    if (centbin == 2) { s.first = 5; s.second = 7; } // 10-40%
    if (centbin == 3) { s.first = 1; s.second = 4; } // 40-80%
    if (centbin == 4) { s.first = 7; s.second = 7; } // 10-20%
    if (centbin == 5) { s.first = 5; s.second = 6; } // 20-40%
    if (centbin == 6) { s.first = 3; s.second = 4; } // 40-60%

    return s;
}

pair <int, int> GetD0PtBin(int d0bin){
    pair <int, int> s;

    int nBinsD0Pt = 11;

    if (d0bin == 0)                      {s.first = 1; s.second = nBinsD0Pt;} // 0 - 10 GeV/c
    if (d0bin >=1 && d0bin <= nBinsD0Pt) {s.first = d0bin; s.second = d0bin;} // Usual Bin Definitions
    if (d0bin == 12)                     {s.first = 3; s.second = nBinsD0Pt;} // 1 - 10 GeV/c
    if (d0bin == 13)                     {s.first = 3; s.second = 8;}         // 1 - 5 GeV/c
    if (d0bin == 14)                     {s.first = 9; s.second = nBinsD0Pt;} // 5 - 10 GeV/c
    return s;
}

THn* SetRangeForTHn(THn *h, int centbin = 0, int d0bin = 0){

    THn *n = (THnF *)h->Clone();

    pair<int, int> d0ptbin  = GetD0PtBin(d0bin);
    pair<int, int> cent = GetCentBin(centbin);

    n->GetAxis(2)->SetRange(d0ptbin.first, d0ptbin.second);
    n->GetAxis(0)->SetRange(cent.first, cent.second);

    return n;   
}

THn *Trim(THn *h, int centbin = 0, int d0bin = 0){

    THn *n1;

    n1 = (THnF *)SetRangeForTHn(h, centbin, d0bin);

    return n1;
}

THn *TrimMassWindow(THn *h, double masslow = 1.8, double masshigh = 1.92){

    THn *n1 = (THnF *)h->Clone();

    int binlow = n1->GetAxis(3)->FindBin(masslow);
    int binhigh = n1->GetAxis(3)->FindBin(masshigh);

    cout << binlow << "\t" << binhigh << endl; 

    n1->GetAxis(3)->SetRange(binlow, binhigh);

    return n1;
}

TH1F* D0JetPtFromFit(int centbin = 0, int d0bin = 0, double signalbg[] = NULL, double bg[] = NULL, double signal[] = NULL, double signalerr[] = NULL){

    TH1::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/D0JetYield.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    // cout << k1->GetEntries() << "\t" << k2->GetEntries() << endl;

    THn *t1 = (THn *)Trim(k1, centbin, d0bin);
    THn *t2 = (THn *)Trim(k2, centbin, d0bin);

    // cout << t1->GetEntries() << "\t" << t2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.92);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.92);

    THn *h1 = (THn *)t1->Clone();
    THn *h2 = (THn *)t2->Clone();

    // h1->GetAxis(3)->SetRangeUser(1.8, 1.92);
    // h2->GetAxis(3)->SetRangeUser(1.8, 1.92);

    // cout << h1->GetEntries() << "\t" << h2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.7, 2.1);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.7, 2.1);

    TH1F *jetpt   = (TH1F *)h1->Projection(1, "E");
    TH1F *jetptbg = (TH1F *)h2->Projection(1, "E");

    // cout << jetpt->GetEntries() << "\t" << jetptbg->GetEntries() << endl;

    jetpt->SetNameTitle(Form("jetpt_%i_%i", centbin, d0bin), Form("jetpt_%i_%i", centbin, d0bin));
    jetptbg->SetNameTitle(Form("jetptbg_%i_%i", centbin, d0bin), Form("jetptbg_%i_%i", centbin, d0bin));

    const int nBinsCent = 3;
    TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80"};

    double lsbackground = jetptbg->Integral();

    // cout << bg[d0bin] << "\t" << lsbackground << endl;
    jetptbg->Scale(bg[d0bin]/lsbackground);

    // cout << jetpt->GetBinContent(jetpt->GetNbinsX()) << "\t" << jetptbg->GetBinContent(jetpt->GetNbinsX()) << endl;

    jetpt->Add(jetptbg, -1);

    for (int i = 1; i<= jetpt->GetNbinsX(); i++){
        if (jetpt->GetBinContent(i) < 0) {jetpt->SetBinContent(i, 0); jetpt->SetBinError(i, 0);}
    }

    return jetpt;
}

void CompareTwoDifferentEfficiencyMethods(){

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const int nBinsCent = 3;
    const int nBinsD0Pt = 14;

    TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80"};

    double yield[nBinsCent + 1][nBinsD0Pt + 1];
    double yielderr[nBinsCent + 1][nBinsD0Pt + 1];
    double signalbg[nBinsCent + 1][nBinsD0Pt + 1];
    double bg[nBinsCent + 1][nBinsD0Pt + 1];

    ifstream YieldFiles;

    for (int i = 0; i <= 0; i++){
        cout << centbinsname[i] << endl;
        YieldFiles.open(Form("../D0JetYield/InvMass/Yield/%s.txt", centbinsname[i].Data()), ios::binary);

        for (int d0ptbin = 0; d0ptbin <= nBinsD0Pt; d0ptbin++){
            YieldFiles >> yield[i][d0ptbin] >> yielderr[i][d0ptbin] >> signalbg[i][d0ptbin] >> bg[i][d0ptbin];
        }
        YieldFiles.close();
    }

    TH1F *HistWithoutEfficiencyCorrection[1][nBinsD0Pt + 1];

    for (int centbin = 0; centbin <= 0; centbin++){
        for (int d0bin = 0; d0bin <= nBinsD0Pt; d0bin++){
            HistWithoutEfficiencyCorrection[0][d0bin] = (TH1F *)D0JetPtFromFit(centbin, d0bin, signalbg[centbin], bg[centbin], yield[centbin], yielderr[centbin]);
        }
    }

    TFile* fout = new TFile("../D0ReconstructionEfficiency/effplotsredone.root");
    TH1F *effinusualbins = (TH1F *)fout->Get("heffBinD0_7");
    TH1F *effinbiggerbins = (TH1F *)fout->Get("heffBinForD0Jet_7");

    const int nBinsJetPt = 10;
    double JetPtBins[nBinsJetPt+1] = {3., 4.,5.,7.,9.,11.,13.,15.,20.,30.,50.};

    TH1F *Histogram1 = new TH1F("h1", "h1", nBinsJetPt, JetPtBins);
    TH1F *Histogram2 = new TH1F("h2", "h2", nBinsJetPt, JetPtBins);

    cout << Histogram1->GetEntries() << "\t" << Histogram2->GetEntries() << endl;

    // Histogram 1 is in the usual efficiency bins
    for (int d0bin = 3; d0bin <= 8; d0bin++){
        TH1F *temp = new TH1F(*Histogram1);
        for (int i = 1; i <= temp->GetNbinsX(); i++){
            temp->SetBinContent(i, effinusualbins->GetBinContent(d0bin));
            temp->SetBinError(i, effinusualbins->GetBinError(d0bin));
        }

        TH1F *temp2 = (TH1F *)HistWithoutEfficiencyCorrection[0][d0bin]->Clone();

        temp2->Divide(temp);

        Histogram1->Add(temp2);
    }

    const int myD0BinsForEfficiencyCorrection[7] = {3,4,5,6,7,8,14};

    // Histogram 2 is in the new efficiency bins
    for (int d0bin = 0; d0bin < 6; d0bin++){
        TH1F *temp = new TH1F(*Histogram2);

        int usuald0bintobecalled = myD0BinsForEfficiencyCorrection[d0bin];

        for (int i = 1; i <= temp->GetNbinsX(); i++){
            temp->SetBinContent(i, effinbiggerbins->GetBinContent(d0bin + 3));
            temp->SetBinError(i, effinbiggerbins->GetBinError(d0bin + 3));
        }

        TH1F *temp2 = (TH1F *)HistWithoutEfficiencyCorrection[0][usuald0bintobecalled]->Clone();

        temp2->Divide(temp);

        Histogram2->Add(temp2);
    }

    TCanvas *c[2];

    c[0] = new TCanvas("c0", "c0", 5, 5, 500, 500);
    c[0]->cd();
    Histogram1->SetLineColor(1);
    Histogram1->SetMarkerColor(1);
    Histogram1->SetMarkerStyle(24);
    Histogram1->SetNameTitle("Usual Efficiency Bins", "Usual Efficiency Bins");

    Histogram2->SetLineColor(2);
    Histogram2->SetMarkerColor(2);
    Histogram2->SetMarkerStyle(24);
    Histogram2->SetNameTitle("New Efficiency Bins", "New Efficiency Bins");

    // Histogram1->GetYaxis()->SetRangeUser(0.1, pow(10, 5));
    Histogram1->GetXaxis()->SetTitle("p_{T, Jet} [GeV/#it{c}]");
    Histogram1->GetYaxis()->SetTitle("N_{Jet}");

    Histogram1->Draw();
    Histogram2->Draw("SAME");
    gPad->SetLogy();
    gPad->BuildLegend();

    c[1] = new TCanvas("c1", "c1", 5, 5, 500, 500);
    c[1]->cd();

    TH1F *Ratio = (TH1F *)Histogram2->Clone("Ratio");
    Ratio->Divide(Histogram1);
    Ratio->GetYaxis()->SetRangeUser(0.5, 2.5);
    Ratio->GetXaxis()->SetTitle("p_{T, Jet} [GeV/#it{c}]");
    Ratio->GetYaxis()->SetTitle("Ratio");
    Ratio->Draw();

    TF1 *one = new TF1("One", "1", 0, 100);
    one->SetLineColor(kBlack);
    one->SetLineStyle(3);
    one->SetLineWidth(1);

    one->Draw("SAME");

    MakeGeneralPlots(10, "p_{T, Jet} [GeV/#it{c}]", "N_{Jet}", "pTSpectraWithoutEfficiency",
        HistWithoutEfficiencyCorrection[0][3], HistWithoutEfficiencyCorrection[0][4], 
        HistWithoutEfficiencyCorrection[0][5], HistWithoutEfficiencyCorrection[0][6], 
        HistWithoutEfficiencyCorrection[0][7], HistWithoutEfficiencyCorrection[0][8],
        HistWithoutEfficiencyCorrection[0][9], HistWithoutEfficiencyCorrection[0][10],
        HistWithoutEfficiencyCorrection[0][11], HistWithoutEfficiencyCorrection[0][14]);

    // MakeGeneralPlots(1, "p_{T, Jet} [GeV/#it{c}]", "N_{Jet}", "pTSpectraWithoutEfficiency",
    //     HistWithoutEfficiencyCorrection[0][14]);

}