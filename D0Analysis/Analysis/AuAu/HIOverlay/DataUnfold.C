#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"
// #include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
#ifndef __CINT__
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
#include "THn.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TPaveText.h"
// #include "StJetTreeStruct.h"
#include <vector>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"


#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"

#pragma link C++ class vector<int> +;

using namespace std;
// using namespace RooFit;

// #pragma link C++ class StJetTreeStruct+;

// #pragma link C++ class vector<float> +;
// #pragma link C++ class vector<vector<float> >+;
// #pragma link C++ class vector<int> +;
// #pragma link C++ class vector<vector<int> >+;
#endif

#include "BinDef.h"

void Method(int iteration, int numberofiterations = 4){

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile *g2 = new TFile("../SPlotFrameWork/ApplyWeights/Histograms3_D01_10GeV_With2DJetPtZ.root");
    g2->cd();

    TH2D *RecoJetPtvZ[3];
    
    RecoJetPtvZ[0] = (TH2D *)gDirectory->Get("ZPt_0_10");
    RecoJetPtvZ[1] = (TH2D *)gDirectory->Get("ZPt_10_40");
    RecoJetPtvZ[2] = (TH2D *)gDirectory->Get("ZPt_40_80");    


    TFile *g;
    g = new TFile(Form("Closure/ResponsePlots_%i_%i.root", iteration,numberofiterations));
    g->cd();

    RooUnfoldResponse *responseJetPtvZ[3];
    RooUnfoldResponse *responseJetPt[3];

    for (int i = 0; i < 3; i++){
        g->GetObject(Form("ResponseJetPt_%i", i), responseJetPt[i]);
        g->GetObject(Form("ResponseJetPtvZ_%i", i), responseJetPtvZ[i]);
    }

    TFile *fin1 = new TFile("../Files/D0Yield_Feb22.root");
    fin1->cd("D0Tagger");
    TH1F *hCent = (TH1F *)gDirectory->Get("hCentralityWeightedAfterCuts");

    // for(int i =1;i<hCent->GetNbinsX()+1;i++)cout << hCent->GetBinContent(i) << endl;
    float nevts_0_10  = hCent->Integral(1, 2);
    float nevts_10_40 = hCent->Integral(3, 8);
    float nevts_40_80 = hCent->Integral(9, 16);

    double s[3];
    s[0] = 941.23714*nevts_0_10; //(1048*1.38188e+06+838*1.22039e+06)/(1.38188e+06+1.22039e+06);                                                                                                                              
    s[1] = 391.35550*nevts_10_40;//571 + 351 + 206;                                                                                                                                                                           
    s[2] = 56.62475*nevts_40_80;

    cout << s[0] << "\t" << s[1] << "\t" << s[2] << endl;

    TH2D *UnfoldedJetPtvZ[3];

    RooUnfoldBayes unfold0 (responseJetPtvZ[0], RecoJetPtvZ[0], numberofiterations);
    RooUnfoldBayes unfold1 (responseJetPtvZ[1], RecoJetPtvZ[1], numberofiterations);
    RooUnfoldBayes unfold2 (responseJetPtvZ[2], RecoJetPtvZ[2], numberofiterations);

    UnfoldedJetPtvZ[0] = (TH2D *)unfold0.Hreco();
    UnfoldedJetPtvZ[1] = (TH2D *)unfold1.Hreco();
    UnfoldedJetPtvZ[2] = (TH2D *)unfold2.Hreco();

    TH1D *RecoJetPt[3];
    TH1D *RecoJetZ[3];

    TH1D *UnfoldedJetPt[3];
    TH1D *UnfoldedJetZ[3];

    TH1D *RatioJetPt[3];
    TH1D *RatioJetZ[3];

    for (int i = 0; i < 3; i++){
        RecoJetPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)RecoJetPtvZ[i]->ProjectionX());
        RecoJetPt[i]->SetNameTitle(Form("RecoJetPt_%i_%i", i, iteration), Form("RecoJetPt_%i_%i", i, iteration));
        RecoJetZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)RecoJetPtvZ[i]->ProjectionY());
        RecoJetZ[i]->SetNameTitle(Form("RecoJetZ_%i_%i", i, iteration), Form("RecoJetZ_%i_%i", i, iteration));

        UnfoldedJetPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)UnfoldedJetPtvZ[i]->ProjectionX());
        UnfoldedJetPt[i]->SetNameTitle(Form("UnfoldedJetPt_%i_%i", i, iteration), Form("UnfoldedJetPt_%i_%i", i, iteration));
        UnfoldedJetZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)UnfoldedJetPtvZ[i]->ProjectionY());
        UnfoldedJetZ[i]->SetNameTitle(Form("UnfoldedJetZ_%i_%i", i, iteration), Form("UnfoldedJetZ_%i_%i", i, iteration));
    }

    TPaveText *pt = new TPaveText(0.20,0.65,0.30,0.85, "NDC"); // NDC sets coords
                                       // relative to pad dimensions
    pt->SetFillStyle(4000);
    pt->SetTextSize(0.04); 
    pt->SetTextAlign(12);

    auto pt_text0 = pt->AddText(Form("SuperIter = %i, Iter = %i", iteration,numberofiterations));

    TCanvas *b = new TCanvas(Form("JetPt_%i", iteration), Form("JetPt_%i", iteration), 2000, 800);
    b->Divide(3);
    for (int i = 0; i < 3; i++){
        b->cd(i+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);

        RecoJetPt[i]->Draw("P SAME");
        RecoJetPt[i]->SetLineColor(kRed);
        RecoJetPt[i]->SetMarkerColor(kRed);
        RecoJetPt[i]->SetMarkerStyle(20);

        UnfoldedJetPt[i]->Draw("P SAME");
        UnfoldedJetPt[i]->SetLineColor(kBlack);
        UnfoldedJetPt[i]->SetMarkerColor(kBlack);
        UnfoldedJetPt[i]->SetMarkerStyle(20);

        RecoJetPt[i]->GetYaxis()->SetTitle("Jet Spectra");
        RecoJetPt[i]->GetXaxis()->SetTitle("p_{T,jet} [GeV/#it{c}]");

        pt->Draw("SAME");
    }

    b->SaveAs(Form("Data/UnfoldedPtSpectra_%i_%i.pdf", iteration, numberofiterations));

    delete b;

    TCanvas *c = new TCanvas(Form("JetZ_%i", iteration), Form("JetZ_%i", iteration), 2000, 800);
    c->Divide(3);
    for (int i = 0; i < 3; i++){
        c->cd(i+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);

        RecoJetZ[i]->Draw("P SAME");
        RecoJetZ[i]->SetLineColor(kRed);
        RecoJetZ[i]->SetMarkerColor(kRed);
        RecoJetZ[i]->SetMarkerStyle(20);

        UnfoldedJetZ[i]->Draw("P SAME");
        UnfoldedJetZ[i]->SetLineColor(kBlack);
        UnfoldedJetZ[i]->SetMarkerColor(kBlack);
        UnfoldedJetZ[i]->SetMarkerStyle(kBlack);

        RecoJetZ[i]->GetYaxis()->SetTitle("Z Count");
        RecoJetZ[i]->GetXaxis()->SetTitle("Z");

        pt->Draw("SAME");
    }

    c->SaveAs(Form("Data/UnfoldedZSpectra_%i_%i.pdf", iteration, numberofiterations));

    delete c;

    RatioJetPt[0] = (TH1D *)UnfoldedJetPt[0]->Clone(Form("RCPJetPt_%i_%i", 0, iteration)); //Central/Peripheral
    RatioJetPt[1] = (TH1D *)UnfoldedJetPt[1]->Clone(Form("RCPJetPt_%i_%i", 1, iteration)); //MidCentral/Peripheral
    RatioJetPt[2] = (TH1D *)UnfoldedJetPt[0]->Clone(Form("RCPJetPt_%i_%i", 2, iteration)); //Central/MidCentral

    RatioJetPt[0]->SetNameTitle(Form("RCPJetPt_%i_%i", 0, iteration), Form("RCPJetPt_%i_%i", 0, iteration));
    RatioJetPt[1]->SetNameTitle(Form("RCPJetPt_%i_%i", 1, iteration), Form("RCPJetPt_%i_%i", 1, iteration));
    RatioJetPt[2]->SetNameTitle(Form("RCPJetPt_%i_%i", 2, iteration), Form("RCPJetPt_%i_%i", 2, iteration));

    RatioJetPt[0]->Divide(UnfoldedJetPt[2]);
    RatioJetPt[1]->Divide(UnfoldedJetPt[2]);
    RatioJetPt[2]->Divide(UnfoldedJetPt[1]);


    RatioJetZ[0] = (TH1D *)UnfoldedJetZ[0]->Clone(Form("RCPJetZ_%i_%i", 0, iteration)); //Central/Peripheral
    RatioJetZ[1] = (TH1D *)UnfoldedJetZ[1]->Clone(Form("RCPJetZ_%i_%i", 1, iteration)); //MidCentral/Peripheral
    RatioJetZ[2] = (TH1D *)UnfoldedJetZ[0]->Clone(Form("RCPJetZ_%i_%i", 2, iteration)); //Central/MidCentral

    RatioJetZ[0]->SetNameTitle(Form("RCPJetZ_%i_%i", 0, iteration), Form("RCPJetZ_%i_%i", 0, iteration));
    RatioJetZ[1]->SetNameTitle(Form("RCPJetZ_%i_%i", 1, iteration), Form("RCPJetZ_%i_%i", 1, iteration));
    RatioJetZ[2]->SetNameTitle(Form("RCPJetZ_%i_%i", 2, iteration), Form("RCPJetZ_%i_%i", 2, iteration));

    RatioJetZ[0]->Divide(UnfoldedJetZ[2]);
    RatioJetZ[1]->Divide(UnfoldedJetZ[2]);
    RatioJetZ[2]->Divide(UnfoldedJetZ[1]);

    TFile *f = new TFile(Form("Data/DataPlots_%i.root", numberofiterations), "UPDATE");
    f->cd();
    for (int i = 0; i < 3; i++){

        RecoJetPt[i]->Write();
        RecoJetZ[i]->Write();

        UnfoldedJetPt[i]->Write();
        UnfoldedJetZ[i]->Write();

        RatioJetPt[i]->Write();
        RatioJetZ[i]->Write();
    }
    f->Close();

    for (int i = 0; i < 3; i++){

        delete RecoJetPt[i];
        delete RecoJetZ[i];

        delete UnfoldedJetPt[i];
        delete UnfoldedJetZ[i];

        delete RatioJetPt[i];
        delete RatioJetZ[i];
    }

    g->Close();

    g2->Close();
}


void DataUnfold(const int NSUPERITERATIONS = 2, int numberofiterations = 4){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
    int colors[200];

    for (int i = 0; i < 6; i++){
        for (int j = -9; j <= 10; j++){
            colors[i*20 + (j+9)] = col[i] + j; 
        }
    }

    TFile *out = new TFile(Form("Data/DataPlots_%i.root", numberofiterations), "RECREATE");
    out->Close();

    for (int i = 0; i <= NSUPERITERATIONS; i++){
        cout << "Number of Super Iteration = " << i << endl;
        Method(i, numberofiterations);
    }
}