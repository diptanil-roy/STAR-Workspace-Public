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

    // cout << binlow << "\t" << binhigh << endl; 

    n1->GetAxis(3)->SetRange(binlow, binhigh);

    return n1;
}

TH1F* D0JetPt(int centbin = 0, int d0bin = 0, double bg = 0., double bgerr = 0., double signal = 0., double signalerr = 0.){

    TH1::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/SignalCount/D0JetsDistribution.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    // cout << k1->GetEntries() << "\t" << k2->GetEntries() << endl;

    THn *t1 = (THn *)Trim(k1, centbin, d0bin);
    THn *t2 = (THn *)Trim(k2, centbin, d0bin);

    // cout << t1->GetEntries() << "\t" << t2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.94);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.94);

    THn *h1 = (THn *)TrimMassWindow(t1, 1.7, 2.1);
    THn *h2 = (THn *)TrimMassWindow(t2, 1.7, 2.1);


    // cout << h1->GetEntries() << "\t" << h2->GetEntries() << endl;

    TH1F *jetpt   = (TH1F *)h1->Projection(1, "E");
    TH1F *jetptbg = (TH1F *)h2->Projection(1, "E");

    cout << jetpt->GetNbinsX() << endl;

    jetpt->SetNameTitle(Form("jetpt_%i_%i", centbin, d0bin), Form("jetpt_%i_%i", centbin, d0bin));
    jetptbg->SetNameTitle(Form("jetptbg_%i_%i", centbin, d0bin), Form("jetptbg_%i_%i", centbin, d0bin));

    // if (centbin == 0) cout << "Last Bin: " << d0bin << "\t" << jetpt->GetBinContent(jetpt->GetNbinsX()) << "\t" << jetptbg->GetBinContent(jetptbg->GetNbinsX()) << endl;

    double lsbackground = jetptbg->Integral();

    if (lsbackground != 0) jetptbg->Scale(1./lsbackground);

    for (int i = 1; i <= jetptbg->GetNbinsX(); i++){
        // if (centbin == 3) cout << bg << "\t" << bgerr << endl;
        double bincontent = jetptbg->GetBinContent(i)*bg;
        double binerror;
        if (bincontent != 0.) binerror = bincontent*(TMath::Sqrt(pow(jetptbg->GetBinError(i)/jetptbg->GetBinContent(i), 2) + pow(bgerr/bg, 2)));
        else binerror = 0;
        // if (centbin == 3) cout << bincontent << "\t" << binerror << endl;
        jetptbg->SetBinContent(i, bincontent);
        jetptbg->SetBinError(i, binerror);
    }


    jetpt->Add(jetptbg, -1);

    for (int i = 1; i<= jetpt->GetNbinsX(); i++){
        if (jetpt->GetBinContent(i) < 0) {jetpt->SetBinContent(i, 0.); jetpt->SetBinError(i, 0.);}
    }

    return jetpt;
}

TH2F* D0JetPtvdR(int centbin = 0, int d0bin = 0, double bg = 0., double bgerr = 0., double signal = 0., double signalerr = 0.){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/SignalCount/D0JetsDistribution.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    // cout << k1->GetEntries() << "\t" << k2->GetEntries() << endl;

    THn *t1 = (THn *)Trim(k1, centbin, d0bin);
    THn *t2 = (THn *)Trim(k2, centbin, d0bin);

    // cout << t1->GetEntries() << "\t" << t2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.94);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.94);

    THn *h1 = (THn *)TrimMassWindow(t1, 1.7, 2.1);
    THn *h2 = (THn *)TrimMassWindow(t2, 1.7, 2.1);

    // cout << h1->GetEntries() << "\t" << h2->GetEntries() << endl;

    TH2F *jetptvdR   = (TH2F *)h1->Projection(4, 1, "E");
    TH2F *jetptvdRbg = (TH2F *)h2->Projection(4, 1, "E");

    // cout << jetpt->GetEntries() << "\t" << jetptbg->GetEntries() << endl;

    jetptvdR->SetNameTitle(Form("jetptvdR_%i_%i", centbin, d0bin), Form("jetptvdR_%i_%i", centbin, d0bin));
    jetptvdRbg->SetNameTitle(Form("jetptvdRbg_%i_%i", centbin, d0bin), Form("jetptvdRbg_%i_%i", centbin, d0bin));

    // if (centbin == 0) cout << "Last Bin: " << d0bin << "\t" << jetpt->GetBinContent(jetpt->GetNbinsX()) << "\t" << jetptbg->GetBinContent(jetptbg->GetNbinsX()) << endl;

    double lsbackground = jetptvdRbg->Integral();

    if (lsbackground != 0) jetptvdRbg->Scale(1./lsbackground);

    for (int i = 1; i <= jetptvdRbg->GetNbinsX(); i++){
        for (int j = 1; j <= jetptvdRbg->GetNbinsY(); j++){
            // if (centbin == 3) cout << bg << "\t" << bgerr << endl;
            double bincontent = jetptvdRbg->GetBinContent(i,j)*bg;
            double binerror;
            if (bincontent != 0.) binerror = bincontent*(TMath::Sqrt(pow(jetptvdRbg->GetBinError(i,j)/jetptvdRbg->GetBinContent(i,j), 2) + pow(bgerr/bg, 2)));
            else binerror = 0;
            // if (centbin == 3) cout << bincontent << "\t" << binerror << endl;
            jetptvdRbg->SetBinContent(i, j, bincontent);
            jetptvdRbg->SetBinError(i, j, binerror);
        }
    }

    jetptvdR->Add(jetptvdRbg, -1);

    for (int i = 1; i<= jetptvdR->GetNbinsX(); i++){
        for (int j = 1; j<= jetptvdR->GetNbinsX(); j++){
            if (jetptvdR->GetBinContent(i,j) < 0) {jetptvdR->SetBinContent(i,j, 0.); jetptvdR->SetBinError(i,j, 0.);}
        }
    }

    return jetptvdR;
}


TH2F* D0JetPtvZ(int centbin = 0, int d0bin = 0, double bg = 0., double bgerr = 0., double signal = 0., double signalerr = 0.){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/SignalCount/D0JetsDistribution.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    // cout << k1->GetEntries() << "\t" << k2->GetEntries() << endl;

    THn *t1 = (THn *)Trim(k1, centbin, d0bin);
    THn *t2 = (THn *)Trim(k2, centbin, d0bin);

    // cout << t1->GetEntries() << "\t" << t2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.94);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.94);

    THn *h1 = (THn *)TrimMassWindow(t1, 1.7, 2.1);
    THn *h2 = (THn *)TrimMassWindow(t2, 1.7, 2.1);


    // cout << h1->GetEntries() << "\t" << h2->GetEntries() << endl;

    TH2F *jetptvZ   = (TH2F *)h1->Projection(5, 1, "E");
    TH2F *jetptvZbg = (TH2F *)h2->Projection(5, 1, "E");

    // cout << jetpt->GetEntries() << "\t" << jetptbg->GetEntries() << endl;

    jetptvZ->SetNameTitle(Form("jetptvZ_%i_%i", centbin, d0bin), Form("jetptvZ_%i_%i", centbin, d0bin));
    jetptvZbg->SetNameTitle(Form("jetptvZbg_%i_%i", centbin, d0bin), Form("jetptvZbg_%i_%i", centbin, d0bin));

    // if (centbin == 0) cout << "Last Bin: " << d0bin << "\t" << jetpt->GetBinContent(jetpt->GetNbinsX()) << "\t" << jetptbg->GetBinContent(jetptbg->GetNbinsX()) << endl;

    double lsbackground = jetptvZbg->Integral();

    if (lsbackground != 0) jetptvZbg->Scale(1./lsbackground);

    for (int i = 1; i <= jetptvZbg->GetNbinsX(); i++){
        for (int j = 1; j <= jetptvZbg->GetNbinsY(); j++){
            // if (centbin == 3) cout << bg << "\t" << bgerr << endl;
            double bincontent = jetptvZbg->GetBinContent(i,j)*bg;
            double binerror;
            if (bincontent != 0.) binerror = bincontent*(TMath::Sqrt(pow(jetptvZbg->GetBinError(i,j)/jetptvZbg->GetBinContent(i,j), 2) + pow(bgerr/bg, 2)));
            else binerror = 0;
            // if (centbin == 3) cout << bincontent << "\t" << binerror << endl;
            jetptvZbg->SetBinContent(i, j, bincontent);
            jetptvZbg->SetBinError(i, j, binerror);
        }
    }

    jetptvZ->Add(jetptvZbg, -1);

    for (int i = 1; i<= jetptvZ->GetNbinsX(); i++){
        for (int j = 1; j<= jetptvZ->GetNbinsX(); j++){
            if (jetptvZ->GetBinContent(i,j) < 0) {jetptvZ->SetBinContent(i,j, 0.); jetptvZ->SetBinError(i,j, 0.);}
        }
    }

    return jetptvZ;
}

TH2F* D0PtvZ(int centbin = 0, int d0bin = 0, double bg = 0., double bgerr = 0., double signal = 0., double signalerr = 0.){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/SignalCount/D0JetsDistribution.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    // cout << k1->GetEntries() << "\t" << k2->GetEntries() << endl;

    THn *t1 = (THn *)Trim(k1, centbin, d0bin);
    THn *t2 = (THn *)Trim(k2, centbin, d0bin);

    // cout << t1->GetEntries() << "\t" << t2->GetEntries() << endl;

    THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.94);
    THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.94);

    // cout << h1->GetEntries() << "\t" << h2->GetEntries() << endl;

    TH2F *d0ptvZ   = (TH2F *)h1->Projection(5, 2, "E");
    TH2F *d0ptvZbg = (TH2F *)h2->Projection(5, 2, "E");

    // cout << jetpt->GetEntries() << "\t" << jetptbg->GetEntries() << endl;

    d0ptvZ->SetNameTitle(Form("D0ptvZ_%i_%i", centbin, d0bin), Form("D0ptvZ_%i_%i", centbin, d0bin));
    d0ptvZbg->SetNameTitle(Form("D0ptvZbg_%i_%i", centbin, d0bin), Form("D0ptvZbg_%i_%i", centbin, d0bin));

    // if (centbin == 0) cout << "Last Bin: " << d0bin << "\t" << jetpt->GetBinContent(jetpt->GetNbinsX()) << "\t" << jetptbg->GetBinContent(jetptbg->GetNbinsX()) << endl;

    double lsbackground = d0ptvZbg->Integral();

    if (lsbackground != 0) d0ptvZbg->Scale(1./lsbackground);

    for (int i = 1; i <= d0ptvZbg->GetNbinsX(); i++){
        for (int j = 1; j <= d0ptvZbg->GetNbinsY(); j++){
            // if (centbin == 3) cout << bg << "\t" << bgerr << endl;
            double bincontent = d0ptvZbg->GetBinContent(i,j)*bg;
            double binerror;
            if (bincontent != 0.) binerror = bincontent*(TMath::Sqrt(pow(d0ptvZbg->GetBinError(i,j)/d0ptvZbg->GetBinContent(i,j), 2) + pow(bgerr/bg, 2)));
            else binerror = 0;
            // if (centbin == 3) cout << bincontent << "\t" << binerror << endl;
            d0ptvZbg->SetBinContent(i, j, bincontent);
            d0ptvZbg->SetBinError(i, j, binerror);
        }
    }

    d0ptvZ->Add(d0ptvZbg, -1);

    for (int i = 1; i<= d0ptvZ->GetNbinsX(); i++){
        for (int j = 1; j<= d0ptvZ->GetNbinsX(); j++){
            if (d0ptvZ->GetBinContent(i,j) < 0) {d0ptvZ->SetBinContent(i,j, 0.); d0ptvZ->SetBinError(i,j, 0.);}
        }
    }

    return d0ptvZ;
}

void MakeD0JetSpectra(){
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const int nBinsCent = 4;
    const int nBinsD0Pt = 8;

    TString centbinsname[nBinsCent] = {"0-80", "0-10", "10-40", "40-80"};

    double yield[nBinsCent][nBinsD0Pt];
    double yielderr[nBinsCent][nBinsD0Pt];
    double background[nBinsCent][nBinsD0Pt];
    double backgrounderr[nBinsCent][nBinsD0Pt];

    ifstream YieldFiles;

    for (int i = 0; i < nBinsCent; i++){
        cout << centbinsname[i] << endl;
        YieldFiles.open(Form("../D0JetYield/SignalCount/RooYield/%s.txt", centbinsname[i].Data()), ios::binary);

        for (int d0ptbin = 0; d0ptbin < nBinsD0Pt; d0ptbin++){
            YieldFiles >> yield[i][d0ptbin] >> yielderr[i][d0ptbin] >> background[i][d0ptbin] >> backgrounderr[i][d0ptbin];

            // if (i == 3) cout << background[i][d0ptbin] << "\t" << backgrounderr[i][d0ptbin] << endl;
        }
        YieldFiles.close();
    }

    TFile* fout = new TFile("../D0ReconstructionEfficiency/effplotsredone.root");
    // TH1F *effinusualbins = (TH1F *)fout->Get("heffBinD0_7");
    TH1F *effinbiggerbins[nBinsCent];

    const int binnumberingforefficiencyhistograms[nBinsCent] = {7, 0, 5, 6};

    for (int i = 0; i < nBinsCent; i++){
        effinbiggerbins[i] = (TH1F *)fout->Get(Form("heffBinForD0Jet_%i", binnumberingforefficiencyhistograms[i]));
        effinbiggerbins[i]->SetDirectory(0);
    }

    fout->Close();

    TH1F *D0PtSpectra[nBinsCent];  //Plots to Save

    const int nBinsD0Spectra = 7;
    double D0PtBins[nBinsD0Spectra + 1] = {1., 1.5, 2., 2.5, 3., 4., 5., 10};

    for (int i = 0; i < nBinsCent; i++){
        D0PtSpectra[i] = new TH1F(Form("D0PtSpectra_%s", centbinsname[i].Data()), Form("D0PtSpectra_%s", centbinsname[i].Data()), nBinsD0Spectra, D0PtBins);

        TH1F *temp = new TH1F(*D0PtSpectra[i]);

        for (int bin = 1; bin <= D0PtSpectra[i]->GetNbinsX(); bin++){
            D0PtSpectra[i]->SetBinContent(bin, yield[i][bin]);
            D0PtSpectra[i]->SetBinError(bin, yielderr[i][bin]);

            temp->SetBinContent(bin, effinbiggerbins[i]->GetBinContent(bin + 2));
            temp->SetBinError(bin, effinbiggerbins[i]->GetBinError(bin + 2));
        }

        D0PtSpectra[i]->Divide(temp);
    }

     MakeGeneralPlots(4, "p_{T, D^{0}} [GeV/#it{c}]", "N_{D^{0}}", "D0SpectraWithEfficiency",
                     D0PtSpectra[0], D0PtSpectra[1], D0PtSpectra[2], D0PtSpectra[3]);

    TH1F *JetPtWithoutEfficiencyCorrection[nBinsCent][nBinsD0Pt];
    TH1F *JetPtSpectra[nBinsCent][3]; // Cent: 0-80, 0-10, 10-40, 40-80 || D0Pt: 1-10 GeV, 1-4 GeV, 4-10 GeV //Plots to Save

    const int nBinsJetPt = 16;
    double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0.,1.,2., 3.,4.,5.,7.,9.,11.,13.,15.,20.,30., 50.};

    // const int nBinsJetPt = 11;
    // double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0., 1., 2., 3., 4., 5.,9.,15., 30.};

    // const int nBinsJetPt = 12;
    // double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0., 1., 2., 3., 4., 5., 8., 12., 16., 30.};

    // const int nBinsJetPt = 11;
    // double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0., 1., 2., 3., 4., 5., 8., 10., 30.};

    // const int nBinsJetPt = 11;
    // double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0., 1., 2., 3., 4., 5., 8., 10., 30.};

    // const int nBinsJetPt = 7;
    // double JetPtBins[nBinsJetPt+1] = {-30., -10., 0., 3., 5., 7., 10., 30.};

    // const int nBinsJetPt = 9;
    // double JetPtBins[nBinsJetPt+1] = {-30., -10., 0., 3., 5., 7., 10., 15., 20., 30.};

    // const int nBinsJetPt = 13;
    // double JetPtBins[nBinsJetPt+1] = {0., 1., 2., 3.,4.,5.,7.,9.,11.,13.,15.,20.,30., 50.};

    int allowedD0PtBinsForRooFit[] = {0, 3, 4, 5, 6, 7, 8, 14};

    for (int i = 0; i < nBinsCent; i++){
        for (int d0ptbin = 0; d0ptbin < nBinsD0Pt; d0ptbin++){
            JetPtWithoutEfficiencyCorrection[i][d0ptbin] = (TH1F *)D0JetPt(i, allowedD0PtBinsForRooFit[d0ptbin], background[i][d0ptbin], backgrounderr[i][d0ptbin], yield[i][d0ptbin], yielderr[i][d0ptbin]);
        } 
    }

    for (int i = 0; i < nBinsCent; i++){

        JetPtSpectra[i][0] = new TH1F(Form("JetPt_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-10"), Form("JetPt_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-10"), nBinsJetPt, JetPtBins);
        JetPtSpectra[i][1] = new TH1F(Form("JetPt_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-4"), Form("JetPt_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-4"), nBinsJetPt, JetPtBins);
        JetPtSpectra[i][2] = new TH1F(Form("JetPt_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "4-10"), Form("JetPt_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "4-10"), nBinsJetPt, JetPtBins);

        for (int d0ptbin = 1; d0ptbin < nBinsD0Pt; d0ptbin++){ // Not adding Bin 0 because that's 1-10 GeV
            
            TH1F *temp = new TH1F(*JetPtWithoutEfficiencyCorrection[i][d0ptbin]);

            for (int jetptbin = 1; jetptbin <= temp->GetNbinsX(); jetptbin++){
                // temp->SetBinContent(jetptbin, effinbiggerbins[i]->GetBinContent(d0ptbin + 2)*temp->GetBinWidth(jetptbin)); //(+2 because the first two bins are not counted)
                // temp->SetBinError(jetptbin, effinbiggerbins[i]->GetBinError(d0ptbin + 2)*temp->GetBinWidth(jetptbin));

                temp->SetBinContent(jetptbin, effinbiggerbins[i]->GetBinContent(d0ptbin + 2));
                temp->SetBinError(jetptbin, effinbiggerbins[i]->GetBinError(d0ptbin + 2));
            }

            TH1F *temp2 = (TH1F *)JetPtWithoutEfficiencyCorrection[i][d0ptbin]->Clone();
            temp2->Divide(temp);

            JetPtSpectra[i][0]->Add(temp2);
            if (d0ptbin >= 1 && d0ptbin <= 5) JetPtSpectra[i][1]->Add(temp2);
            else if (d0ptbin >= 6 && d0ptbin <= 7) JetPtSpectra[i][2]->Add(temp2);
        }
    }

    MakeGeneralPlots(3, "p_{T, Jet} [GeV/#it{c}]", "N_{Jet}", "pTSpectraWithEfficiency_1-10",
                     JetPtSpectra[1][0], JetPtSpectra[2][0], JetPtSpectra[3][0]);
    MakeGeneralPlots(3, "p_{T, Jet} [GeV/#it{c}]", "N_{Jet}", "pTSpectraWithEfficiency_1-4",
                     JetPtSpectra[1][1], JetPtSpectra[2][1], JetPtSpectra[3][1]);
    MakeGeneralPlots(3, "p_{T, Jet} [GeV/#it{c}]", "N_{Jet}", "pTSpectraWithEfficiency_4-10",
                     JetPtSpectra[1][2], JetPtSpectra[2][2], JetPtSpectra[3][2]);

    const Int_t nBinsdR = 4;
    Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.2, 0.4};

    TH2F *JetPtvdRWithoutEfficiencyCorrection[nBinsCent][nBinsD0Pt];
    TH2F *JetPtvdRSpectra[nBinsCent][3]; // Cent: 0-80, 0-10, 10-40, 40-80 || D0Pt: 1-10 GeV, 1-4 GeV, 4-10 GeV //Plots to Save

    for (int i = 0; i < nBinsCent; i++){
        for (int d0ptbin = 0; d0ptbin < nBinsD0Pt; d0ptbin++){
            JetPtvdRWithoutEfficiencyCorrection[i][d0ptbin] = (TH2F *)D0JetPtvdR(i, allowedD0PtBinsForRooFit[d0ptbin], background[i][d0ptbin], backgrounderr[i][d0ptbin], yield[i][d0ptbin], yielderr[i][d0ptbin]);
        } 
    }

    for (int i = 0; i < nBinsCent; i++){

        JetPtvdRSpectra[i][0] = new TH2F(Form("JetPtvdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-10"), Form("JetPtvdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-10"), nBinsJetPt, JetPtBins, nBinsdR, dRBins);
        JetPtvdRSpectra[i][1] = new TH2F(Form("JetPtvdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-4"), Form("JetPtvdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-4"), nBinsJetPt, JetPtBins, nBinsdR, dRBins);
        JetPtvdRSpectra[i][2] = new TH2F(Form("JetPtvdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "4-10"), Form("JetPtvdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "4-10"), nBinsJetPt, JetPtBins, nBinsdR, dRBins);

        for (int d0ptbin = 1; d0ptbin < nBinsD0Pt; d0ptbin++){ // Not adding Bin 0 because that's 1-10 GeV
            
            TH2F *temp = new TH2F(*JetPtvdRWithoutEfficiencyCorrection[i][d0ptbin]);

            for (int jetptbin = 1; jetptbin <= temp->GetNbinsX(); jetptbin++){
                for (int drbin = 1; drbin <= temp->GetNbinsY(); drbin++){
                    temp->SetBinContent(jetptbin, drbin, effinbiggerbins[i]->GetBinContent(d0ptbin + 2)); //(+2 because the first two bins are not counted)
                    temp->SetBinError(jetptbin, drbin, effinbiggerbins[i]->GetBinError(d0ptbin + 2));
                }
            }

            TH2F *temp2 = (TH2F *)JetPtvdRWithoutEfficiencyCorrection[i][d0ptbin]->Clone();
            temp2->Divide(temp);

            JetPtvdRSpectra[i][0]->Add(temp2);
            if (d0ptbin >= 1 && d0ptbin <= 5) JetPtvdRSpectra[i][1]->Add(temp2);
            else if (d0ptbin >= 6 && d0ptbin <= 7) JetPtvdRSpectra[i][2]->Add(temp2);
        }
    }

    TString finalD0BinNames[3] = {"1-10", "1-4", "4-10"};

    for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
        TCanvas *c1 = new TCanvas("c1", "c1", 5, 5, 1000, 500);
        c1->Divide(4);
        for (int i = 0; i < nBinsCent; i++){
            c1->cd(i + 1);
            JetPtvdRSpectra[i][finalD0Bin]->Draw("COLZ");
            JetPtvdRSpectra[i][finalD0Bin]->GetYaxis()->SetTitle("#DeltaR");
            JetPtvdRSpectra[i][finalD0Bin]->GetXaxis()->SetTitle("p_{T, Jet} [GeV/#it{c}]");
            gPad->SetLogz();
        }

        c1->SaveAs(Form("./JetPtvdR_%s.pdf", finalD0BinNames[finalD0Bin].Data()));
        delete c1;
    }



    // const int nBinsZ = 17;
    // double ZBins[nBinsZ+1] = {-10.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 10.0};

    // TH2F *JetPtvZWithoutEfficiencyCorrection[nBinsCent][nBinsD0Pt];
    // TH2F *JetPtvZSpectra[nBinsCent][3]; // Cent: 0-80, 0-10, 10-40, 40-80 || D0Pt: 1-10 GeV, 1-4 GeV, 4-10 GeV //Plots to Save

    // for (int i = 0; i < nBinsCent; i++){
    //     for (int d0ptbin = 0; d0ptbin < nBinsD0Pt; d0ptbin++){
    //         JetPtvZWithoutEfficiencyCorrection[i][d0ptbin] = (TH2F *)D0JetPtvZ(i, allowedD0PtBinsForRooFit[d0ptbin], background[i][d0ptbin], backgrounderr[i][d0ptbin], yield[i][d0ptbin], yielderr[i][d0ptbin]);
    //     } 
    // }

    // for (int i = 0; i < nBinsCent; i++){

    //     JetPtvZSpectra[i][0] = new TH2F(Form("JetPtvZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-10"), Form("JetPtvZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-10"), nBinsJetPt, JetPtBins, nBinsZ, ZBins);
    //     JetPtvZSpectra[i][1] = new TH2F(Form("JetPtvZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-4"), Form("JetPtvZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "1-4"), nBinsJetPt, JetPtBins, nBinsZ, ZBins);
    //     JetPtvZSpectra[i][2] = new TH2F(Form("JetPtvZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "4-10"), Form("JetPtvZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), "4-10"), nBinsJetPt, JetPtBins, nBinsZ, ZBins);

    //     for (int d0ptbin = 1; d0ptbin < nBinsD0Pt; d0ptbin++){ // Not adding Bin 0 because that's 1-10 GeV
            
    //         TH2F *temp = new TH2F(*JetPtvZWithoutEfficiencyCorrection[i][d0ptbin]);

    //         for (int jetptbin = 1; jetptbin <= temp->GetNbinsX(); jetptbin++){
    //             for (int zbin = 1; zbin <= temp->GetNbinsY(); zbin++){
    //                 temp->SetBinContent(jetptbin, zbin, effinbiggerbins[i]->GetBinContent(d0ptbin + 2)); //(+2 because the first two bins are not counted)
    //                 temp->SetBinError(jetptbin, zbin, effinbiggerbins[i]->GetBinError(d0ptbin + 2));
    //             }
    //         }

    //         TH2F *temp2 = (TH2F *)JetPtvZWithoutEfficiencyCorrection[i][d0ptbin]->Clone();
    //         temp2->Divide(temp);

    //         JetPtvZSpectra[i][0]->Add(temp2);
    //         if (d0ptbin >= 1 && d0ptbin <= 5) JetPtvZSpectra[i][1]->Add(temp2);
    //         else if (d0ptbin >= 6 && d0ptbin <= 7) JetPtvZSpectra[i][2]->Add(temp2);
    //     }
    // }

    // for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
    //     TCanvas *c1 = new TCanvas("c1", "c1", 5, 5, 1000, 500);
    //     c1->Divide(4);
    //     for (int i = 0; i < nBinsCent; i++){
    //         c1->cd(i + 1);
    //         JetPtvZSpectra[i][finalD0Bin]->Draw("COLZ");
    //         JetPtvZSpectra[i][finalD0Bin]->GetYaxis()->SetTitle("Z");
    //         JetPtvZSpectra[i][finalD0Bin]->GetXaxis()->SetTitle("p_{T, Jet} [GeV/#it{c}]");
    //         gPad->SetLogz();
    //     }

    //     c1->SaveAs(Form("./JetPtvZ_%s.pdf", finalD0BinNames[finalD0Bin].Data()));
    //     delete c1;
    // }


    TH1F *dR[nBinsCent][3]; //Plots to Save
    // TH1F *dZ[nBinsCent][3]; //Plots to Save

    
    for (int i = 0; i < nBinsCent; i++){
        for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
            dR[i][finalD0Bin] = (TH1F *)JetPtvdRSpectra[i][finalD0Bin]->ProjectionY();
            // dZ[i][finalD0Bin] = (TH1F *)JetPtvZSpectra[i][finalD0Bin]->ProjectionY();

            dR[i][finalD0Bin]->SetNameTitle(Form("JetdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), finalD0BinNames[finalD0Bin].Data()), Form("JetdR_Cent_%s_D0Pt_%s", centbinsname[i].Data(), finalD0BinNames[finalD0Bin].Data()));
            // dZ[i][finalD0Bin]->SetNameTitle(Form("JetdZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), finalD0BinNames[finalD0Bin].Data()), Form("JetdZ_Cent_%s_D0Pt_%s", centbinsname[i].Data(), finalD0BinNames[finalD0Bin].Data()));
        }
    }

    // TH2F *D0PtvZWithoutEfficiencyCorrection[nBinsCent][nBinsD0Pt];
    // TH2F *D0PtvZSpectra[nBinsCent]; // Cent: 0-80, 0-10, 10-40, 40-80 || D0Pt: 1-10 GeV, 1-4 GeV, 4-10 GeV //Plots to Save

    // for (int i = 0; i < nBinsCent; i++){
    //     for (int d0ptbin = 0; d0ptbin < nBinsD0Pt; d0ptbin++){
    //         D0PtvZWithoutEfficiencyCorrection[i][d0ptbin] = (TH2F *)D0PtvZ(i, allowedD0PtBinsForRooFit[d0ptbin], background[i][d0ptbin], backgrounderr[i][d0ptbin], yield[i][d0ptbin], yielderr[i][d0ptbin]);
    //     } 
    // }

    // for (int i = 0; i < nBinsCent; i++){

    //     D0PtvZSpectra[i]= new TH2F(Form("D0PtvZ_Cent_%s", centbinsname[i].Data()), Form("D0PtvZ_Cent_%s", centbinsname[i].Data()), nBinsD0Spectra, D0PtBins, nBinsZ, ZBins);

    //     TH2F *temp = new TH2F(*D0PtvZSpectra[i]);

    //     for (int d0ptbin = 1; d0ptbin <= nBinsD0Spectra; d0ptbin++){ // Not adding Bin 0 because that's 1-10 GeV 
    //         for (int zbin = 1; zbin <= nBinsZ; zbin++){
    //             temp->SetBinContent(d0ptbin, zbin, effinbiggerbins[i]->GetBinContent(d0ptbin + 2)); //(+2 because the first two bins are not counted)
    //             temp->SetBinError(d0ptbin, zbin, effinbiggerbins[i]->GetBinError(d0ptbin + 2));
    //         }
    //     }

    //     TH2F *temp2 = (TH2F *)D0PtvZSpectra[i]->Clone();
    //     temp2->Divide(temp);

    //     D0PtvZSpectra[i]->Add(temp2);
    // }

    cout << "I am filling the outputFile" << endl;

    TFile *outputFile = new TFile("./D0JetsFoldedHistograms.root", "RECREATE");
    outputFile->cd();

    for (int i = 0; i < nBinsCent; i++){
        D0PtSpectra[i]->Write();
        // D0PtvZSpectra[i]->Write();
        for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
            JetPtSpectra[i][finalD0Bin]->Write();
        }
        for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
            JetPtvdRSpectra[i][finalD0Bin]->Write();  
        }
        // for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
        //     JetPtvZSpectra[i][finalD0Bin]->Write();  
        // }
        for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
            dR[i][finalD0Bin]->Write();  
        }
        // for (int finalD0Bin = 0; finalD0Bin < 3; finalD0Bin++){
        //     dZ[i][finalD0Bin]->Write();  
        // }
    }
    outputFile->Close();

}

void MakeD0JetPlots(){
    MakeD0JetSpectra();
}