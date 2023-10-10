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

    n->GetAxis(1)->SetRange(d0ptbin.first, d0ptbin.second);
    n->GetAxis(0)->SetRange(cent.first, cent.second);

    return n;   
}

THn *Trim(THn *h, int centbin = 0, int d0bin = 0){

    THn *n1;

    n1 = (THnF *)SetRangeForTHn(h, centbin, d0bin);

    return n1;
}

vector<double> invmassplotter(int centbin = 0, int d0bin = 0){ // Signal SignalErr SigBg Bg

    TH1::SetDefaultSumw2();

	cout << "========================================================" << endl;

	cout << "D0Bin = " << d0bin << endl;

	double d0yield = 0;

	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    string CentralityBins_Name[5] = {"0-10", "10-20", "20-40", "40-60", "60-80"};

    TFile *f2 =  TFile::Open("./D0JetYield.root");
    f2->cd();

    THn *k1 = (THnF *)gDirectory->Get("hD0CentPtEtaMDphiDaug_Standard");
    THn *k2 = (THnF *)gDirectory->Get("hD0CentPtEtaMDphiDaugLikeSign_Standard");

    cout << k1->GetEntries() << endl;
    cout << k2->GetEntries() << endl;

    // k1->GetAxis(5)->SetRange(1,4); //|eta| < 1
    // k2->GetAxis(5)->SetRange(1,4);

    int daughter1ptbin = k1->GetAxis(2)->FindBin(0.6+1e-6);
    int daughter2ptbin = k1->GetAxis(4)->FindBin(0.6+1e-6);

    cout << daughter1ptbin << "\t" << daughter2ptbin << endl;

    // k1->GetAxis(2)->SetRange(daughter1ptbin, 4);
    // k1->GetAxis(4)->SetRange(daughter2ptbin, 4);

    // k2->GetAxis(2)->SetRange(daughter1ptbin, 4);
    // k2->GetAxis(4)->SetRange(daughter2ptbin, 4);

    THn *h1 = (THn*) Trim(k1, centbin, d0bin);
    THn *h2 = (THn*) Trim(k2, centbin, d0bin);


    TH1F *fullinvmass = (TH1F *)h1->Projection(3, "E");
    TH1F *fullinvmassbg = (TH1F *)h2->Projection(3, "E");

    TString c0_filename = "./InvMass/";
    c0_filename += "D0_InvMass_Distribution_AuAu";
    c0_filename += "_Cent_";
    c0_filename += centbin;
    c0_filename += "_D0Pt_";
    c0_filename += d0bin;

    int rebincounter[] = {1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1};
  
    TCanvas *c0 = new TCanvas(c0_filename.Data(), c0_filename.Data(), 5,5,800,800);

    c0->cd();
    gPad->SetRightMargin(0.02);
    gPad->SetLeftMargin(0.12);

    TH1F *signalSelectedD0Invmass;
    TH1F *backgroundLSSelectedD0Invmass;

    signalSelectedD0Invmass = (TH1F *)fullinvmass->Clone();
    backgroundLSSelectedD0Invmass = (TH1F *)fullinvmassbg->Clone();
    
    signalSelectedD0Invmass->Rebin(rebincounter[d0bin]);
    signalSelectedD0Invmass->Draw();
    backgroundLSSelectedD0Invmass->Rebin(rebincounter[d0bin]);
    backgroundLSSelectedD0Invmass->Draw("SAME");

    signalSelectedD0Invmass->SetLineColor(4);
    signalSelectedD0Invmass->SetMarkerColor(4);

    signalSelectedD0Invmass->SetMarkerStyle(24);
    backgroundLSSelectedD0Invmass->SetMarkerStyle(20);
    backgroundLSSelectedD0Invmass->SetLineColor(2);
    backgroundLSSelectedD0Invmass->SetMarkerColor(2);


    signalSelectedD0Invmass->GetXaxis()->SetRangeUser(1.7, 2.10);
    signalSelectedD0Invmass->GetYaxis()->SetTitle("Counts/(20 MeV/c^{2})");
    signalSelectedD0Invmass->GetXaxis()->SetTitle("m (GeV/c^{2})");

    Int_t binlow = signalSelectedD0Invmass->FindBin(1.82);
    Int_t binhigh = signalSelectedD0Invmass->FindBin(1.91);

    bool dofit = kFALSE;

    for (int bin = binlow; bin < binhigh; bin++){
        if (signalSelectedD0Invmass->GetBinContent(bin) >  10) dofit = kTRUE;
        // dummy_splusb += signalSelectedD0Invmass->GetBinContent(bin);
        // dummy_b += backgroundLSSelectedD0Invmass->GetBinContent(bin);

        if (dofit) break;
    }

    ////////////// FIT Section //////////////////////

    TF1 *background = new TF1("background", "[0]+[1]*x+[2]*x**2+[3]*x**3", 1.7, 2.10);
    background->SetLineColor(1);
    // TF1 *fitfunc = new TF1("fitfunc", "signal + background(4)", 1.6, 2.1);
    TF1 *fitfunc = new TF1("fitfunc", "gausn + pol3(3)", 1.7, 2.10);
    fitfunc->SetLineColor(6);
    fitfunc->SetLineWidth(5);

    TH1D *invmass = (TH1D *)signalSelectedD0Invmass->Clone();
    // invmass->GetXaxis()->SetRangeUser(1.7, 2.10);

    TH1D *invmassbg = (TH1D *)backgroundLSSelectedD0Invmass->Clone();

    Double_t SplusB_whole = 0;
    Double_t B_whole = 0;
    Double_t S_whole = 0;
    Double_t S_whole_err = 0;

    double dx = invmass->GetXaxis()->GetBinWidth(1);
    cout << dx << endl;
    cout << invmass->GetEntries() << endl;
    cout << invmassbg->GetEntries() << endl;

    if (dofit) {
    	cout << "STARTING FIT" << endl;
        invmassbg->Fit("background", "RLQ0");

        TF1 *backgroundfit = invmassbg->GetFunction("background");

        fitfunc->SetParLimits(0, 0, 200000000);

        fitfunc->SetParameter(1, 1.865);
        fitfunc->SetParLimits(1, 1.86, 1.87);
        fitfunc->SetParLimits(2, 0.01, 0.03);

        fitfunc->SetParameter(3, backgroundfit->GetParameter(0));
        fitfunc->SetParameter(4, backgroundfit->GetParameter(1));
        fitfunc->SetParameter(5, backgroundfit->GetParameter(2));
        fitfunc->SetParameter(6, backgroundfit->GetParameter(3));

        int status = invmass->Fit("fitfunc", "RLQ");

        if (fitfunc->GetChisquare()/fitfunc->GetNDF() > 1.5) status = 1;

        double fitloops = 0;

        while (status != 0){
        if (fitloops >= 1000) break;
        	// cout << fitloops << "\t";
            fitloops += 1;
            fitfunc->SetParameter(0, fitfunc->GetParameter(0));
            fitfunc->SetParLimits(0, 0, 200000000);
            fitfunc->SetParameter(1, fitfunc->GetParameter(1));
            fitfunc->SetParameter(1, 1.865);
            fitfunc->SetParLimits(1, 1.86, 1.87);
            fitfunc->SetParLimits(2, 0.01, 0.03);
            // fitfunc->SetParameter(2, 0.2);
            fitfunc->SetParameter(3, fitfunc->GetParameter(3));
            fitfunc->SetParameter(4, fitfunc->GetParameter(4));
            fitfunc->SetParameter(5, fitfunc->GetParameter(5));
            fitfunc->SetParameter(6, fitfunc->GetParameter(6));

            status = invmass->Fit("fitfunc", "RLQ");

            if (fitfunc->GetChisquare()/fitfunc->GetNDF() > 2.0) status = 1;
        }

        // cout << endl;
        // cout << "Status = " << status << endl;


        TF1 *truebackground = new TF1("truebackground", "[0] + [1]*x + [2]*x**2 + [3]*x**3", 1.7, 2.10);
        truebackground->SetParameter(0, fitfunc->GetParameter(3));
        truebackground->SetParameter(1, fitfunc->GetParameter(4));
        truebackground->SetParameter(2, fitfunc->GetParameter(5));
        truebackground->SetParameter(3, fitfunc->GetParameter(6));

        TF1 *signal = new TF1("signal", "gausn", 1.7, 2.10);
        signal->SetParameter(0, fitfunc->GetParameter(0));
        signal->SetParameter(1, fitfunc->GetParameter(1));
        signal->SetParameter(2, fitfunc->GetParameter(2));

        signal->SetLineStyle(9);
        signal->SetLineColorAlpha(kBlue, 0.35);
        signal->SetFillStyle(1001);
        signal->SetFillColorAlpha(kBlue, 0.20);

        S_whole = fitfunc->GetParameter(0)/dx;
        S_whole_err = fitfunc->GetParError(0)/dx;

        for (int i = invmass->FindBin(1.7); i <= invmass->FindBin(2.10); i++){
            SplusB_whole += fitfunc->Eval(invmass->GetBinCenter(i));
            B_whole += truebackground->Eval(invmass->GetBinCenter(i));
        }

        invmass->GetYaxis()->SetMaxDigits(3);

        invmass->SetLineColor(4);
        invmass->SetMarkerColor(4);

        invmass->SetMarkerStyle(24);
        invmassbg->SetMarkerStyle(20);

        invmassbg->SetLineColor(2);
        invmassbg->SetMarkerColor(2);

        gStyle->SetTitleAlign(23);
        gStyle->SetTitleFontSize(0.03);

        invmass->SetTitle("#sqrt{s_{NN}} = 200 GeV AuAu");
        invmass->GetYaxis()->SetTitle("Counts/(20 MeV/#it{c}^{2})");
        invmass->GetXaxis()->SetTitle("m_{K#pi} [GeV/#it{c}^{2}]");
        invmass->GetYaxis()->SetRangeUser(0, 2*invmass->GetMaximum());
        // if (directorytostart == 5 && !useefficiency_correction) invmass->GetYaxis()->SetRangeUser(0, 50);
        invmass->GetYaxis()->SetNdivisions(505);
        invmass->Draw("EP");
        invmassbg->Draw("EP SAME");

        truebackground->Draw("E4 SAME");

        signal->Draw("E4C SAME");

        auto legend = new TLegend(0.90,0.65,0.57,0.88);
        // legend->SetHeader("InvMass","C"); // option "C" allows to center the header
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        legend->AddEntry(invmass,"Data","lep");
        legend->AddEntry(fitfunc, "Fit", "l");
        legend->AddEntry(invmassbg,"Combinatorial","lep");
        // legend->AddEntry(background, "Background Fit", "l");
        legend->AddEntry(signal, "Signal", "f");
        legend->Draw("SAME");

        TPaveText *pt = new TPaveText(0.15,0.6,0.55,0.85, "NDC"); // NDC sets coords
                                              // relative to pad dimensions
        // pt->SetFillColor(0); // text is black on white
        // pt->SetBorderColor(1);
        pt->SetTextSize(0.03); 
        pt->SetTextAlign(12);

        auto pt_text0 = pt->AddText(Form("M_{D^{0} + #bar{D^{0}}} = %1.3f GeV/#it{c}^{2}", fitfunc->GetParameter(1)));
        auto pt_text1 = pt->AddText(Form("#sigma = %1.4f GeV/#it{c}^{2} ; S = %4.0f", fitfunc->GetParameter(2), S_whole));
        auto pt_text2 = pt->AddText(Form("S/B = %1.2f ; S/#sqrt{S + B} = %2.2f", S_whole/B_whole, S_whole/TMath::Sqrt(SplusB_whole)));
        auto pt_text3 = pt->AddText(Form("#frac{#chi^{2}}{NDF} = %1.2f", fitfunc->GetChisquare()/fitfunc->GetNDF()));

        pt->Draw("SAME");

    }

    else{
        for (int bin = invmass->FindBin(1.7); bin <= invmass->FindBin(2.10); bin++){

			S_whole += (signalSelectedD0Invmass->GetBinContent(bin) - backgroundLSSelectedD0Invmass->GetBinContent(bin));
            S_whole_err += pow(invmass->GetBinError(bin), 2);

            SplusB_whole += signalSelectedD0Invmass->GetBinContent(bin);
            B_whole += backgroundLSSelectedD0Invmass->GetBinContent(bin);
        }
    }

    cout << "Total D0s = " << S_whole << "\t" << fitfunc->GetParError(0) << "\t" << SplusB_whole << "\t" <<  B_whole << endl;
    // cout << "Projected D0s = " << S_whole*875.0*pow(10,6)/totalentries << endl;

    TPaveText *data_range = new TPaveText(0.60,0.45,0.80,0.60, "NDC");
    data_range->SetFillColor(0); // text is black on white
    data_range->SetTextSize(0.03); 
    data_range->SetTextAlign(12);

    TString d0ptbins[15] = {"0.0 - 10.0", "0.0 - 0.5", "0.5 - 1.0", "1.0 - 1.5", "1.5 - 2.0", "2.0 -2.5", "2.5 - 3.0", "3.0 - 4.0", "4.0 - 5.0", "5.0 - 6.0", "6.0 - 8.0", "8.0 - 10.0", "1.0 - 10.0", "1.0 - 5.0", "5.0 - 10.0"};
    TString centbins[7] = {"0 - 80 \%", "0 - 10 \%", "10 - 40 \%", "40 - 80 \%", "10 - 20 \%", "20 - 40 \%", "40 - 60 \%"};
    

    auto data_range_text0 = data_range->AddText(Form("Centrality #in [%s]", centbins[centbin].Data()));
    auto data_range_text2 = data_range->AddText(Form("p_{T,K#pi} [GeV/#it{c}] #in %s", d0ptbins[d0bin].Data()));
    // auto data_range_text3 = data_range->AddText("|#eta| < 1");

    data_range->Draw("SAME");

    c0_filename += ".pdf";

    c0->SaveAs(c0_filename.Data());

    delete c0;

    vector<double> s;
    s.clear();
    s.push_back(S_whole);
    s.push_back(S_whole_err);
    s.push_back(SplusB_whole);
    s.push_back(B_whole);

    // cout << "Counts: " << S_whole << "\t" << S_whole_err;

    // cout.precision(17);
    // cout << "0-80 \% events = " << totalentries << endl;

    return s;
}

void D0JetCounts(){

    const int nBinsCent = 6;
    const int nBinsD0Pt = 14;

    TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80", "10-20", "20-40", "40-60"};

    ofstream YieldFiles[nBinsCent + 1];

    for (int centbin = 0; centbin <= 0; centbin++){

        YieldFiles[centbin].open(Form("./InvMass/Yield/%s.txt", centbinsname[centbin].Data()));

        for (int d0bin = 0; d0bin <= nBinsD0Pt; d0bin++){
            vector<double> s = invmassplotter(centbin, d0bin);

            YieldFiles[centbin] << s[0] << "\t" << s[1] << "\t" << s[2] << "\t" << s[3] << endl; 
        }

        YieldFiles[centbin].close();
    }
}