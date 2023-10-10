#ifndef __CINT__
#include <stdexcept>
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

void ResponseMaker(int iteration, int numberofiterations=4){
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // double numberofentries = 10000000;
    double numberofentries = 100;

    TFile *f = new TFile("Plots.root");
    TH1D *NEntries[2];
    TH1D *VzMC[3][2];
    TH1D *MCZ[3][2];
    TH1D *RecoZ[3][2];
    THn  *JetPtZ[3][2]; //Wider range of Z and Jet Pt here.

    TH2D *Resp2MC[3][2];
    TH2D *Resp2Reco[3][2];

    RooUnfoldResponse *response[3][2];
    RooUnfoldResponse *response2[3][2];

    THn  *response2THn[3][2]; //Wider range of Z and Jet Pt here.



    double ratioofeventsthatpassveto[2] = {1070./1000000., 11114./1000000.}; // For scaling, we need to make sure that vetoed events are also included.
    double crosssection[2] = {64.84/66.809, 1.969/66.809};

    for (int bin = 0; bin < 2; bin++){
        f->cd(Form("pthatbin_%i", bin));
        // gROOT->ProcessLine(".ls");
        NEntries[bin] = (TH1D *)gDirectory->Get("NEvents");

        for (int i = 0; i < 3; i++){
          VzMC[i][bin] = (TH1D *)gDirectory->Get(Form("VzMC_%i", i));

          double scale = crosssection[bin]*ratioofeventsthatpassveto[bin]/VzMC[i][bin]->GetEntries();
          cout << scale << endl;

          MCZ[i][bin] = (TH1D *)gDirectory->Get(Form("MCZ_%i", i));
          RecoZ[i][bin] = (TH1D *)gDirectory->Get(Form("RecoZ_%i", i));
          JetPtZ[i][bin] = (THnD *)gDirectory->Get(Form("hMCJetPtMCZRecoJetPtRecoZ_%i", i));

          gDirectory->GetObject(Form("ResponseJetPtvZ_%i", i), response[i][bin]);
          

          gDirectory->GetObject(Form("Response2Reco_%i", i), response2[i][bin]);
          Resp2MC[i][bin] = (TH2D *)gDirectory->Get(Form("resp2MC_%i", i));
          Resp2Reco[i][bin] = (TH2D *)gDirectory->Get(Form("resp2Reco_%i", i));
          response2THn[i][bin] = (THnD *)gDirectory->Get(Form("Response2THn_%i", i));
          
          MCZ[i][bin]->Scale(scale);
          RecoZ[i][bin]->Scale(scale);
          JetPtZ[i][bin]->Scale(scale);

          response[i][bin]->Scale(scale);


          response2[i][bin]->Scale(scale);
          Resp2MC[i][bin]->Scale(scale);
          Resp2Reco[i][bin]->Scale(scale);
          response2THn[i][bin]->Scale(scale);

        }
    }

    for (int i = 0; i < 3; i++){
        MCZ[i][0]->Add(MCZ[i][1]);
        RecoZ[i][0]->Add(RecoZ[i][1]);
        JetPtZ[i][0]->Add(JetPtZ[i][1]);
        response[i][0]->Add(*response[i][1]);

        response2[i][0]->Add(*response2[i][1]);
        Resp2MC[i][0]->Add(Resp2MC[i][1]);
        Resp2Reco[i][0]->Add(Resp2Reco[i][1]);
        response2THn[i][0]->Add(response2THn[i][1]);

    }

    cout << "OG Response Matrices Made" << endl;

    TH2D *UnfoldedJetPtvZ[3];

    RooUnfoldBayes unfold2 (response2[0][0], Resp2Reco[0][0], numberofiterations);

    UnfoldedJetPtvZ[2] = (TH2D *)unfold2.Hreco();
    
    TH2D *tmp[3];

    TCanvas *b = new TCanvas("b", "b", 1600, 800);
    b->Divide(4);

    b->cd(1);
    gPad->SetLogz();
    tmp[0] = (TH2D *)response2[0][0]->Hresponse();
    tmp[0]->Draw("COLZ");

    TH1D *tmp6 = (TH1D *)response2[0][0]->Hmeasured();
    TH1D *tmp7 = (TH1D *)response2[0][0]->Hfakes();
    TH1D *tmp8 = (TH1D *)response2[0][0]->Htruth();

    b->cd(2);
    gPad->SetLogz();
    tmp6->Draw("COLZ");

    b->cd(3);
    gPad->SetLogz();
    tmp7->Draw("COLZ");

    b->cd(4);
    gPad->SetLogz();
    tmp8->Draw("COLZ");

    TCanvas *b2 = new TCanvas("b2", "b2", 1600, 800);
    b2->Divide(3);

    b2->cd(1);
    gPad->SetLogz();
    Resp2MC[0][0]->Draw("COLZ");

    b2->cd(2);
    gPad->SetLogz();
    Resp2Reco[0][0]->Draw("COLZ");

    b2->cd(3);
    gPad->SetLogz();
    UnfoldedJetPtvZ[2]->Draw("COLZ");



    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2);
    c->cd(1);
    gPad->SetLogy();

    TH1D *tmp1 = (TH1D *)UnfoldedJetPtvZ[2]->ProjectionX();

    TH1D *tmpptorg = (TH1D *)Resp2MC[0][0]->ProjectionX();
    TH1D *tmpptreco = (TH1D *)Resp2Reco[0][0]->ProjectionX();

    tmp1->GetYaxis()->SetRangeUser(pow(10, -20), 1);

    tmp1->Draw();
    tmpptorg->Draw("SAME");
    tmpptreco->Draw("SAME");


    tmp1->SetLineColor(kBlack);
    tmp1->SetMarkerColor(kBlack);
    tmp1->SetMarkerSize(1.5);
    tmp1->SetMarkerStyle(24);

    tmpptorg->SetLineColor(kGreen-2);
    tmpptorg->SetMarkerColor(kGreen-2);
    tmpptorg->SetMarkerStyle(20);

    tmpptreco->SetLineColor(kBlue-2);
    tmpptreco->SetMarkerColor(kBlue-2);
    tmpptreco->SetMarkerStyle(20);


    auto legend1 = new TLegend(0.6,0.7,0.85,0.9);
    legend1->AddEntry(tmp1,"Regular Histogram","l");
    legend1->AddEntry(tmpptorg,"Truth","l");
    legend1->AddEntry(tmpptreco,"Detector","l");
    legend1->Draw("SAME");   

    c->cd(2);
    gPad->SetLogy();

    TH1D *tmp3     = (TH1D *)UnfoldedJetPtvZ[2]->ProjectionY();

    TH1D *tmpzorg  = (TH1D *)Resp2MC[0][0]->ProjectionY();
    TH1D *tmpzreco = (TH1D *)Resp2Reco[0][0]->ProjectionY();

    tmp3->GetYaxis()->SetRangeUser(pow(10, -20), 1);
    
    tmp3->Draw();
    tmpzorg->Draw("SAME");
    tmpzreco->Draw("SAME");

    tmp3->SetLineColor(kBlack);
    tmp3->SetMarkerColor(kBlack);
    tmp3->SetMarkerSize(1.5);
    tmp3->SetMarkerStyle(24);

    tmpzorg->SetLineColor(kGreen-2);
    tmpzorg->SetMarkerColor(kGreen-2);
    tmpzorg->SetMarkerStyle(20);

    tmpzreco->SetLineColor(kBlue-2);
    tmpzreco->SetMarkerColor(kBlue-2);
    tmpzreco->SetMarkerStyle(20);

    auto legend2 = new TLegend(0.6,0.7,0.85,0.9);
    legend2->AddEntry(tmp3,"Regular Histogram","lp");
    legend2->AddEntry(tmpzorg,"Truth","lp");
    legend2->AddEntry(tmpzreco,"Reco","lp");
    legend2->Draw("SAME");    

    TH1D *h = (TH1D *)tmpzorg->Clone("Z");
    for (int i = 1; i <= tmpzorg->GetNbinsX(); i++){
        h->SetBinContent(i,tmpzorg->Integral());
    }

    h->Divide(tmpzorg);

    TFile *tmpzfile = new TFile("TMPZFile.root", "RECREATE");
    tmpzfile->cd();
    h->Write();
    tmpzfile->Close();
}

void TestForSampling(){
    ResponseMaker(0, 4);
}