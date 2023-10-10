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

TH1F* D0JetPt(int centbin = 0, int d0bin = 0){

    TH1::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/D0JetYield.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    THn *t1 = (THn *)Trim(k1, centbin, d0bin);
    THn *t2 = (THn *)Trim(k2, centbin, d0bin);

    THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.92);
    THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.92);

    TH1F *jetpt   = (TH1F *)h1->Projection(1, "E");
    TH1F *jetptbg = (TH1F *)h2->Projection(1, "E");

    jetpt->Add(jetptbg, -1);

    cout << jetpt->Integral() << endl;

    // jetpt->Scale(1./jetpt->Integral());

    return jetpt;
}

TH1F* D0JetPtFromFit(int centbin = 0, int d0bin = 0, double signalbg[] = NULL, double bg[] = NULL, double signal[] = NULL, double signalerr[] = NULL){

    TH1::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/D0JetYield.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    cout << k1->GetEntries() << "\t" << k2->GetEntries() << endl;

    THn *t1 = (THn *)Trim(k1, centbin, 0);
    THn *t2 = (THn *)Trim(k2, centbin, 0);

    cout << t1->GetEntries() << "\t" << t2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.92);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.92);

    THn *h1 = (THn *)t1->Clone();
    THn *h2 = (THn *)t2->Clone();

    // h1->GetAxis(3)->SetRangeUser(1.8, 1.92);
    // h2->GetAxis(3)->SetRangeUser(1.8, 1.92);

    cout << h1->GetEntries() << "\t" << h2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.7, 2.1);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.7, 2.1);

    TH1F *jetpt   = (TH1F *)h1->Projection(1, "E");
    TH1F *jetptbg = (TH1F *)h2->Projection(1, "E");

    cout << jetpt->GetEntries() << "\t" << jetptbg->GetEntries() << endl;

    jetptbg->SetNameTitle("jetptbg", "jetptbg");

    const int nBinsCent = 3;
    TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80"};

    double lsbackground = jetptbg->Integral();

    cout << bg[d0bin] << "\t" << lsbackground << endl;
    jetptbg->Scale(bg[d0bin]/lsbackground);

    jetpt->Add(jetptbg, -1);

    return jetpt;
}

TH1F* D0JetPtFromRooFit(int centbin = 0, int d0bin = 0, double bg[] = NULL, double bgerr[] = NULL, double signal[] = NULL, double signalerr[] = NULL){

    TH1::SetDefaultSumw2();

    TFile *f2 =  TFile::Open("../D0JetYield/D0JetYield.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZ_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0JetCentJetPtD0PtMassDeltaRZLikeSign_Standard");

    cout << k1->GetEntries() << "\t" << k2->GetEntries() << endl;

    THn *t1 = (THn *)Trim(k1, centbin, 0);
    THn *t2 = (THn *)Trim(k2, centbin, 0);

    cout << t1->GetEntries() << "\t" << t2->GetEntries() << endl;

    // THn *h1 = (THn *)TrimMassWindow(t1, 1.8, 1.92);
    // THn *h2 = (THn *)TrimMassWindow(t2, 1.8, 1.92);
    THn *h1 = (THn *)t1->Clone();
    THn *h2 = (THn *)t2->Clone();

    cout << h1->GetEntries() << "\t" << h2->GetEntries() << endl;

    TH1F *jetpt   = (TH1F *)h1->Projection(1, "E");
    TH1F *jetptbg = (TH1F *)h2->Projection(1, "E");

    cout << jetpt->GetEntries() << "\t" << jetptbg->GetEntries() << endl;

    jetptbg->SetNameTitle("jetptbg", "jetptbg");

    const int nBinsCent = 3;
    TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80"};

    double lsbackground = jetptbg->Integral();

    cout << bg[d0bin] << "\t" << lsbackground << endl;
    jetptbg->Scale(1./lsbackground);

    //bg[d0bin]

    for (int i = 1; i <= jetptbg->GetNbinsX(); i++){
        double bincontent = jetptbg->GetBinContent(i)*bg[d0bin];
        double binerror = bincontent*(TMath::Sqrt(pow(jetptbg->GetBinError(i)/jetptbg->GetBinContent(i), 2) + pow(bgerr[d0bin]/bg[d0bin], 2)));
        jetptbg->SetBinContent(i, bincontent);
        jetptbg->SetBinError(i, binerror);
    }

    jetpt->Add(jetptbg, -1);

    return jetpt;
}

void CompareErrorsBetweenMethods(){

    TH1::SetDefaultSumw2();

    TFile *f = new TFile("D0JetPt.root");
    f->cd();

    cout << "HERE" << endl;

    TH1F *BinnedRootFit = (TH1F *)gDirectory->Get("JetPtFromFit_0-80_ 0");
    TH1F *UnBinnedRootFit = (TH1F *)gDirectory->Get("JetPtFromRooFit_0-80_ 0");

    BinnedRootFit->SetDirectory(0);
    UnBinnedRootFit->SetDirectory(0);

    f->Close();

    cout << "My Plots are called" << endl;

    TFile *g = new TFile("Matt_Histograms3_NoEff.root");
    g->cd();

    TH1F *sPlot = (TH1F *)gDirectory->Get("JetPt");
    sPlot->SetDirectory(0);

    g->Close();

    cout << "Matt's Plots are called" << endl;

    TH1F *ErrBinned = new TH1F(*BinnedRootFit);
    TH1F *ErrUnBinned = new TH1F(*UnBinnedRootFit);
    TH1F *ErrsPlot = new TH1F(*sPlot);

    for (int i = 1; i <= ErrBinned->GetNbinsX(); i++){
        ErrBinned->SetBinContent(i, BinnedRootFit->GetBinError(i));
        ErrBinned->SetBinError(i, 0.);
    }
    for (int i = 1; i <= ErrUnBinned->GetNbinsX(); i++){
        ErrUnBinned->SetBinContent(i, UnBinnedRootFit->GetBinError(i));
        ErrUnBinned->SetBinError(i, 0.);
    }
    for (int i = 1; i <= ErrsPlot->GetNbinsX(); i++){
        ErrsPlot->SetBinContent(i, sPlot->GetBinError(i));
        ErrsPlot->SetBinError(i, 0.);
    }

    ErrBinned->SetMarkerColor(1);
    ErrBinned->SetMarkerStyle(25);

    ErrUnBinned->SetMarkerColor(1);
    ErrUnBinned->SetMarkerStyle(24);

    ErrsPlot->SetMarkerColor(2);
    ErrsPlot->SetMarkerStyle(20);


    ErrBinned->SetNameTitle("Binned Fit", "Binned Fit");
    ErrUnBinned->SetNameTitle("UnBinned Fit", "UnBinned Fit");
    ErrsPlot->SetNameTitle("sPlot", "sPlot");

    TCanvas *c = new TCanvas("c", "c", 5, 5, 600, 600);
    c->cd();

    ErrBinned->Draw("P");
    ErrUnBinned->Draw("P SAME");
    ErrsPlot->Draw("P SAME");

    ErrBinned->GetYaxis()->SetTitle("Error");
    ErrBinned->GetXaxis()->SetTitle("p_{T, jet} [GeV/#it{c}]");

    auto legend = new TLegend(0.85,0.65,0.52,0.88);
    // legend->SetHeader("InvMass","C"); // option "C" allows to center the header
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(ErrBinned,"Binned Fit","p");
    legend->AddEntry(ErrUnBinned, "Unbinned Fit", "p");
    legend->AddEntry(ErrsPlot,"sPlot","p");
    legend->Draw("SAME");

    TPaveText *data_range = new TPaveText(0.20,0.15,0.40,0.30, "NDC");
    data_range->SetFillColor(0); // text is black on white
    data_range->SetTextSize(0.03); 
    data_range->SetTextAlign(12);

    auto data_range_text0 = data_range->AddText(Form("Centrality #in [%s]", " 0 - 80"));
    auto data_range_text1 = data_range->AddText(Form("p_{T, D^{0}} #in [%s] GeV/#it{c}", "1 - 10"));

    data_range->Draw("SAME");

    // gPad->BuildLegend();
    gPad->SetRightMargin(0.02);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);

    gPad->SetLogy();

    c->SaveAs("../Plots/ErrorsFromDifferentMethods.pdf");


}

void D0JetPtSpectra(){

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const int nBinsCent = 3;
    const int nBinsD0Pt = 14;

    TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80"};

    int colorpallete[] = {1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18};

    TFile *f = new TFile("D0JetPt.root", "RECREATE");
    f->cd();

    double yield[nBinsCent + 1][nBinsD0Pt + 1];
    double yielderr[nBinsCent + 1][nBinsD0Pt + 1];
    double signalbg[nBinsCent + 1][nBinsD0Pt + 1];
    double bg[nBinsCent + 1][nBinsD0Pt + 1];

    double rooyield[nBinsCent + 1][nBinsD0Pt + 1];
    double rooyielderr[nBinsCent + 1][nBinsD0Pt + 1];
    double roobg[nBinsCent + 1][nBinsD0Pt + 1];
    double roobgerr[nBinsCent + 1][nBinsD0Pt + 1];

    ifstream YieldFiles;

    // for (int i = 0; i <= nBinsCent; i++){
    //     cout << centbinsname[i] << endl;
    //     YieldFiles.open(Form("../D0JetYield/InvMass/Yield/%s.txt", centbinsname[i].Data()), ios::binary);

    //     for (int d0ptbin = 0; d0ptbin <= nBinsD0Pt; d0ptbin++){
    //         YieldFiles >> yield[i][d0ptbin] >> yielderr[i][d0ptbin] >> signalbg[i][d0ptbin] >> bg[i][d0ptbin];
    //     }
    //     YieldFiles.close();
    // }

    for (int i = 0; i <= 0; i++){
        cout << centbinsname[i] << endl;
        YieldFiles.open(Form("../D0JetYield/InvMass/Yield/%s.txt", centbinsname[i].Data()), ios::binary);

        for (int d0ptbin = 0; d0ptbin <= nBinsD0Pt; d0ptbin++){
            YieldFiles >> yield[i][d0ptbin] >> yielderr[i][d0ptbin] >> signalbg[i][d0ptbin] >> bg[i][d0ptbin];
        }
        YieldFiles.close();
    }

    for (int i = 0; i <= 0; i++){
        cout << centbinsname[i] << endl;
        YieldFiles.open(Form("../D0JetYield/InvMass/RooYield/%s.txt", centbinsname[i].Data()), ios::binary);

        for (int d0ptbin = 0; d0ptbin <= nBinsD0Pt; d0ptbin++){
            YieldFiles >> rooyield[i][d0ptbin] >> rooyielderr[i][d0ptbin] >> roobg[i][d0ptbin] >> roobgerr[i][d0ptbin];
        }
        YieldFiles.close();
    }

    TCanvas *MethodsComparison = new TCanvas("Jet pT Spectra (Different methods)", "Jet pT Spectra (Different methods)", 5, 5, 600, 600);
    TFile *g = new TFile("Matt_Histograms3_NoEff.root");
    g->cd();
    TH1F *sPlot = (TH1F *)gDirectory->Get("JetPt");
    sPlot->SetNameTitle("sPlot", "sPlot");
    sPlot->SetDirectory(0);
    g->Close();
    

    auto legend = new TLegend(0.88,0.65,0.55,0.88);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(sPlot,"sPlot","lep");

    TPaveText *data_range = new TPaveText(0.20,0.15,0.40,0.30, "NDC");
    data_range->SetFillColor(0); // text is black on white
    data_range->SetTextSize(0.03); 
    data_range->SetTextAlign(12);


    for (int centbin = 0; centbin <= 0; centbin++){
        for (int d0bin = 0; d0bin <= 0; d0bin++){

            

            TH1F *h = (TH1F *)D0JetPtFromFit(centbin, d0bin, signalbg[centbin], bg[centbin], yield[centbin], yielderr[centbin]);
            // TH1F *h = (TH1F *)D0JetPtFromFit(centbin, d0bin, NULL, bg[centbin], yield[centbin], yielderr[centbin]);
            h->SetNameTitle(Form("JetPtFromFit_%s_ %i", centbinsname[centbin].Data(), d0bin), Form("JetPtFromFit_%s_%i", centbinsname[centbin].Data(), d0bin));
            h->SetLineColor(1);
            h->SetMarkerColor(1);
            h->SetMarkerStyle(25);

            

            TH1F *g = (TH1F *)D0JetPt(centbin, d0bin);
            g->SetNameTitle(Form("JetPtFromBinCount_%s_ %i", centbinsname[centbin].Data(), d0bin), Form("JetPtFromBinCount_%s_%i", centbinsname[centbin].Data(), d0bin));
            g->SetLineColor(1);
            g->SetMarkerColor(1);
            g->SetMarkerStyle(24);

            

            TH1F *k = (TH1F *)D0JetPtFromRooFit(centbin, d0bin, roobg[centbin], roobgerr[centbin], rooyield[centbin], rooyielderr[centbin]);
            k->SetNameTitle(Form("JetPtFromRooFit_%s_ %i", centbinsname[centbin].Data(), d0bin), Form("JetPtFromRooFit_%s_%i", centbinsname[centbin].Data(), d0bin));
            k->SetLineColor(1);
            k->SetMarkerColor(1);
            k->SetMarkerStyle(24);

            sPlot->SetMarkerStyle(20);
            sPlot->SetMarkerColor(2);
            sPlot->SetLineColor(2);

            legend->AddEntry(k,"Unbinned Fit","lep");
            legend->AddEntry(h,"Binned Fit","lep");
            legend->AddEntry(g,"Bin Count","lep");

            auto data_range_text0 = data_range->AddText(Form("Centrality #in [%s]", centbinsname[centbin].Data()));
            auto data_range_text1 = data_range->AddText(Form("p_{T, D^{0}} #in [%s] GeV/#it{c}", "1 - 10"));
            
            g->GetYaxis()->SetRangeUser(0.1, pow(10, 5));
            g->GetXaxis()->SetTitle("p_{T, Jet} [GeV/#it{c}]");
            g->GetYaxis()->SetTitle("N_{Jet}");

            g->Draw();
            h->Draw("SAME");
            k->Draw("SAME");
            sPlot->Draw("SAME");
            legend->Draw("SAME");
            data_range->Draw("SAME");

            gPad->SetLogy();

            f->cd();
            g->Write();
            h->Write();
            k->Write();
        }
    }

    MethodsComparison->SaveAs("../Plots/JetPtSpectraBetweenMethods.pdf");

    f->Close();

    CompareErrorsBetweenMethods();

}



// void D0JetR(){

// }

// void D0JetZ(){

// }