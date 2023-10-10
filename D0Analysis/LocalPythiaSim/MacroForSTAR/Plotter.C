#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include <TLorentzVector.h>
#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
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
// #include "Pythia8/Pythia.h"
using namespace std;
#endif

void Plotter(){


	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TStopwatch timer;

  	timer.Start();

  	TFile *f1 = new TFile("PDFSet2.root");
    f1->cd();
    TH1F *h = (TH1F *)gDirectory->Get("hD0Pt");

    TFile *f2 = new TFile("PDFSet13.root");
    f2->cd();
    TH1F *k = (TH1F *)gDirectory->Get("hD0Pt");

    TFile *f3 = new TFile("PurePythia.root");
    f3->cd();
    TH1F *l = (TH1F *)gDirectory->Get("hD0Pt");

    TFile *f4 = new TFile("PythiaSTARTune.root");
    f4->cd();
    TH1F *m = (TH1F *)gDirectory->Get("hD0Pt");
    

    TH1F *FONLL = new TH1F("FONLL","FONLL",60,0,30);
    TH1F *MinFONLL = new TH1F("MIN FONLL","MIN FONLL",60,0,30);
    TH1F *MaxFONLL = new TH1F("MAX FONLL","MAX FONLL",60,0,30);
    ifstream data1("FONLL.txt");
    int cnt=1;
    if(data1.is_open()){
      while(!data1.eof()){
        double x;double y;double minerr; double maxerr; double tmp;
        data1 >> x >> y >> minerr >> maxerr >> tmp >> tmp >> tmp >> tmp;
        FONLL->SetBinContent(cnt,y*pow(10, -9)/x);
        MinFONLL->SetBinContent(cnt,minerr*pow(10, -9)/x);
        MaxFONLL->SetBinContent(cnt,maxerr*pow(10, -9)/x);

        // FONLL->SetBinContent(cnt,y/(7.5330e+07));
        // MinFONLL->SetBinContent(cnt,minerr/(4.4410e+07));
        // MaxFONLL->SetBinContent(cnt,maxerr/(1.7240e+08));

        cnt++;
      }
    }

    FONLL->SetLineColor(kRed);
    FONLL->SetMarkerColor(kRed);
    FONLL->SetMarkerStyle(20);
    FONLL->SetMarkerSize(1);

    // MinFONLL->Scale(1./MinFONLL->Integral());
    MinFONLL->SetLineColor(kGreen-2);
    MinFONLL->SetMarkerColor(kGreen-2);
    MinFONLL->SetMarkerStyle(24);

    // MaxFONLL->Scale(1./MaxFONLL->Integral());
    MaxFONLL->SetLineColor(kBlack);
    MaxFONLL->SetMarkerColor(kBlack);
    MaxFONLL->SetMarkerStyle(24);

    MaxFONLL->GetYaxis()->SetTitle("#frac{1}{p_{T}} #frac{d#sigma}{dp_{T}} [mb/GeV]");
    MaxFONLL->GetXaxis()->SetTitle("p_{T}");

    double ptxaxis_published[5] = {1.57, 2.45, 3.44, 4.45, 5.45};
    double dsigmadpt_published[5] = {0.106089296, 0.03107972, 0.00756112, 0.0016320464, 0.00065097852};
    double dsigmadpt_fonll[5] = {1.0515E-01, 3.5115E-02, 9.5505E-03, 2.8380E-03, 9.9716E-04};

    double dsigmaptdpt_published[5]; 
    double dsigmaptdpt_fonll[5]; 

    for (int i = 0; i < 5; i++){
      dsigmaptdpt_published[i] = dsigmadpt_published[i]/ptxaxis_published[i];
      dsigmaptdpt_fonll[i] = dsigmadpt_fonll[i]/ptxaxis_published[i];
    }

    TGraph *publishedspectra = new TGraph(6, ptxaxis_published, dsigmaptdpt_published);
    TGraph *publishedfonll = new TGraph(6, ptxaxis_published, dsigmaptdpt_fonll);
    publishedspectra->SetLineColor(kBlack);
    publishedspectra->SetMarkerColor(kBlack);
    publishedspectra->SetMarkerStyle(30);
    publishedspectra->SetMarkerSize(2);

    publishedfonll->SetLineColor(kRed);
    publishedfonll->SetMarkerColor(kRed);
    publishedfonll->SetMarkerStyle(26);



    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetMarkerStyle(20);

    h->Scale(1./double(1000000));
    // h->Scale(1./105426.);
    h->Scale(1./(0.565*0.5));
    h->Scale(30.);

    k->SetMarkerColor(kRed);
    k->SetLineColor(kRed);
    k->SetMarkerStyle(20);

    k->Scale(1./double(1000000));
    // h->Scale(1./105426.);
    k->Scale(1./(0.565*0.5));
    k->Scale(30.);

    l->SetMarkerColor(kBlue);
    l->SetLineColor(kBlue);
    l->SetMarkerStyle(20);

    l->Scale(1./double(1000000));
    // h->Scale(1./105426.);
    l->Scale(1./(0.565*0.5));
    l->Scale(30.);

    m->SetMarkerColor(kViolet+1);
    m->SetLineColor(kViolet+1);
    m->SetMarkerStyle(29);
    m->SetMarkerSize(2);

    m->Scale(1./double(459909));
    // h->Scale(1./105426.);
    m->Scale(1./(0.565*0.5));
    m->Scale(30.);

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    c->cd();
    gPad->SetLogy();

    MaxFONLL->GetYaxis()->SetRangeUser(5*pow(10, -7), 2);
    MaxFONLL->GetXaxis()->SetRangeUser(0.5, 10);

    MaxFONLL->Draw("C SAME");
    MaxFONLL->SetFillColor(kGreen);
    MaxFONLL->SetLineColor(kGreen);

    MinFONLL->Draw("C SAME");
    MinFONLL->SetFillColor(10);
    MinFONLL->SetLineColor(kGreen);

    h->Draw("EP SAME");
    k->Draw("EP SAME");
    l->Draw("EP SAME");
    m->Draw("EP SAME");

    publishedspectra->Draw("P SAME");
    // publishedfonll->Draw("P SAME");

    auto legend = new TLegend(0.5,0.7,0.9,0.9);
  	legend->AddEntry(MaxFONLL,"FONLL","f");
  	legend->AddEntry(h,"PYTHIA D^{0}/0.565 (CTEQ 5L)","p");
  	legend->AddEntry(k,"PYTHIA D^{0}/0.565 (MONASH (k_{T} changed))","p");
  	legend->AddEntry(l,"PYTHIA D^{0}/0.565 (Out-of-box)","p");
  	legend->AddEntry(m,"PYTHIA D^{0}/0.565 (STAR Tune)","p");
  	legend->AddEntry(publishedspectra,"D^{0} Run 09","p");

  	legend->Draw("SAME");
}