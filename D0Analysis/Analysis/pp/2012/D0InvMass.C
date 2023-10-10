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
#include <TString.h>

#pragma link C++ class vector<int> +;

using namespace std;

double M_PION_PLUS = 0.139570;
double M_KAON_PLUS = 0.493677;

pair<TH1D *, TH1D *> Method(TString TreeName = "D0USTree"){

	TFile *inputfile = new TFile("ppQA.root");
	TDirectory *QADir = (TDirectory *)inputfile->Get("PIDQA_TrackTree");
	TTree *TrackTree = (TTree *)QADir->Get(TreeName.Data());

	// Declaration of leaf types
	Float_t         D0Mass;
	Float_t         D0Pt;
	Float_t         D0Eta;
	Float_t         D0Phi;
	Float_t         PionPt;
	Float_t         PionEta;
	Float_t         PionPhi;
	Float_t         PionDCA;
	Float_t         PionNHitsFit;
	Float_t         PionNSigmaPion;
	Float_t         PionNSigmaKaon;
	Float_t         PionTofBeta;
	Float_t         PionTofYLocal;
	Float_t         KaonPt;
	Float_t         KaonEta;
	Float_t         KaonPhi;
	Float_t         KaonDCA;
	Float_t         KaonNHitsFit;
	Float_t         KaonNSigmaPion;
	Float_t         KaonNSigmaKaon;
	Float_t         KaonTofBeta;
	Float_t         KaonTofYLocal;


	TrackTree->SetBranchAddress("D0Mass", &D0Mass);
	TrackTree->SetBranchAddress("D0Pt", &D0Pt);
	TrackTree->SetBranchAddress("D0Eta", &D0Eta);
	TrackTree->SetBranchAddress("D0Phi", &D0Phi);
	TrackTree->SetBranchAddress("PionPt", &PionPt);
	TrackTree->SetBranchAddress("PionEta", &PionEta);
	TrackTree->SetBranchAddress("PionPhi", &PionPhi);
	TrackTree->SetBranchAddress("PionDCA", &PionDCA);
	TrackTree->SetBranchAddress("PionNHitsFit", &PionNHitsFit);
	TrackTree->SetBranchAddress("PionNSigmaPion", &PionNSigmaPion);
	TrackTree->SetBranchAddress("PionNSigmaKaon", &PionNSigmaKaon);
	TrackTree->SetBranchAddress("PionTofBeta", &PionTofBeta);
	TrackTree->SetBranchAddress("PionTofYLocal", &PionTofYLocal);
	TrackTree->SetBranchAddress("KaonPt", &KaonPt);
	TrackTree->SetBranchAddress("KaonEta", &KaonEta);
	TrackTree->SetBranchAddress("KaonPhi", &KaonPhi);
	TrackTree->SetBranchAddress("KaonDCA", &KaonDCA);
	TrackTree->SetBranchAddress("KaonNHitsFit", &KaonNHitsFit);
	TrackTree->SetBranchAddress("KaonNSigmaPion", &KaonNSigmaPion);
	TrackTree->SetBranchAddress("KaonNSigmaKaon", &KaonNSigmaKaon);
	TrackTree->SetBranchAddress("KaonTofBeta", &KaonTofBeta);
	TrackTree->SetBranchAddress("KaonTofYLocal", &KaonTofYLocal);

	int nentries = TrackTree->GetEntriesFast();

	TH1D *hD0Mass = new TH1D("hD0Mass", "hD0Mass", 2500, 1.7, 2.1);
	TH1D *hD0Pt = new TH1D("hD0Pt", "hD0Pt", 20, 0.0, 10.0);

	for (int event = 0; event < nentries;  event++){

  		if (event%500000==0)cout << event << " Events done." << endl;

   		TrackTree->GetEntry(event);

   		if (D0Pt < 1.0) continue;
   		if (abs(D0Eta) > 1.0) continue;

   		double pionp = PionPt*TMath::CosH(PionEta);
   		double kaonp = KaonPt*TMath::CosH(KaonEta);

   		float oneOverBetaExpectedpion = sqrt(M_PION_PLUS*M_PION_PLUS / pionp / pionp + 1);
   		float oneOverBetaExpectedkaon = sqrt(M_KAON_PLUS*M_KAON_PLUS / kaonp / kaonp + 1);

   		// if (abs(PionNSigmaPion) > 2) continue;
	   	// if (abs(KaonNSigmaKaon) > 2) continue;
	   	// if (abs(PionNSigmaPion) > abs(PionNSigmaKaon)) continue;
	   	// if (abs(KaonNSigmaKaon) > abs(KaonNSigmaPion)) continue;
	   	// if (abs(KaonNSigmaPion) < 2) continue;
	   	// if (abs(PionNSigmaKaon) < 2) continue;
	   	if (PionTofBeta < -1.9 || PionTofBeta > 2.1) continue;
	   	if (abs(PionTofYLocal) > 1.8) continue;
	   	if (abs(KaonTofYLocal) > 1.8) continue;
	   	if (PionNHitsFit < 20) continue;
	   	if (KaonNHitsFit < 20) continue;

	   	hD0Mass->Fill(D0Mass);

	   	if (D0Mass > 1.83 && D0Mass < 1.89) hD0Pt->Fill(D0Pt);
   	}

   	TH1D *h = (TH1D *)hD0Mass->Clone();
   	TH1D *g = (TH1D *)hD0Pt->Clone();

   	delete hD0Mass;
   	delete hD0Pt;

   	pair<TH1D *, TH1D *>s(h, g);

   	return s;
}

void SaveTheD0MassHistograms(){
	TFile *out = new TFile("out.root", "RECREATE");

	pair<TH1D *, TH1D *>s1 = Method("D0USTree");
    pair<TH1D *, TH1D *>s2 = Method("D0LSTree");

    out->cd();

    TH1D *US = (TH1D *)s1.first->Clone("US");
	TH1D *LS = (TH1D *)s2.first->Clone("LS");

    US->Write();
    LS->Write();
    out->Close();
}

pair<TH1D *, TH1D *> ReadTheD0MassHistograms(){
	TFile *out = new TFile("out.root");
	out->cd();

	TH1D *US = (TH1D *)out->Get("US");
	TH1D *LS = (TH1D *)out->Get("LS");

	pair<TH1D *, TH1D *>s(US, LS);

	return s;
}

void FitMethod(){

	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // gStyle->SetOptFit(111111);

    // pair<TH1D *, TH1D *>s1 = Method("D0USTree");
    // pair<TH1D *, TH1D *>s2 = Method("D0LSTree");

    pair<TH1D *, TH1D *>s = ReadTheD0MassHistograms();

	TH1D *h1 = (TH1D *)s.first;
	TH1D *h2 = (TH1D *)s.second;

	int lowbin = h2->FindBin(2.0);
	int highbin = h2->FindBin(2.1);

	double usscalefactor = h1->Integral(lowbin, highbin);
	double lsscalefactor = h2->Integral(lowbin, highbin);

	// lowbin = h2->FindBin(1.70);
	// highbin = h2->FindBin(1.80);

	// usscalefactor += h1->Integral(lowbin, highbin);
	// lsscalefactor += h2->Integral(lowbin, highbin);

	cout << usscalefactor << "\t" << lsscalefactor << endl;

	h2->Scale(usscalefactor/lsscalefactor);

	h1->Rebin(125);
	h2->Rebin(125);

	TH1D *US = (TH1D *)h1->Clone("US");
	TH1D *LS = (TH1D *)h2->Clone("LS");

	US->SetLineColor(1);
	US->SetMarkerStyle(20);
	US->SetMarkerColor(1);
	LS->SetLineColor(1);
	LS->SetMarkerStyle(24);
	LS->SetMarkerColor(1);

	h1->Add(h2, -1);

	h1->GetXaxis()->SetTitle("m_{K#pi} [GeV/#it{c}^{2}]");
	h1->GetYaxis()->SetTitle("Counts");

	TCanvas *c = new TCanvas("c", "c", 700, 500);
	c->cd();
	// US->Draw();
	// LS->Draw("SAME");
	h1->Draw("EP SAME");

	TF1 *bffunc = new TF1("bffunc", "pol2", 1.73, 2.05);
	bffunc->SetParameter(3, 2);
    bffunc->SetParameter(4, -1);
    bffunc->SetParameter(5, -1);

    int bfstatus = h1->Fit("bffunc", "L");

    cout << "BF Status = " << bfstatus << endl;
    cout << bffunc->GetParameter(0) << "\t" << bffunc->GetParameter(1) << "\t" << bffunc->GetParameter(2) << endl;

	// TF1 *fitfunc = new TF1("fitfunc", "signal + background(4)", 1.6, 2.1);
    TF1 *fitfunc = new TF1("fitfunc", "gaus + pol2(3)", 1.73, 2.05);
    fitfunc->SetLineColor(2);
    fitfunc->SetLineWidth(5);

    fitfunc->SetParLimits(0, 1, 5000);
    fitfunc->FixParameter(1, 1.865);
    // fitfunc->SetParLimits(1, 1.863, 1.867);
    fitfunc->SetParLimits(2, 0.00001, 0.1);
    fitfunc->SetParameter(3, bffunc->GetParameter(0));
    fitfunc->SetParameter(4, bffunc->GetParameter(1));
    fitfunc->SetParameter(5, bffunc->GetParameter(2));

    int status = h1->Fit("fitfunc", "");

    cout << "Fit Status = " << status << endl;

    double chi2byndf = fitfunc->GetChisquare()/fitfunc->GetNDF();
    double mean = fitfunc->GetParameter(1);
    double sigma = fitfunc->GetParameter(2);

    TPaveText *pt = new TPaveText(0.65,0.65,0.85,0.85, "NDC"); // NDC sets coords
    pt->SetFillColor(0); // text is black on white
    pt->SetTextSize(0.03); 
    pt->SetTextAlign(12);
    auto pt_text0 = pt->AddText("US - Scaled LS");
    auto pt_text1 = pt->AddText(Form("p_{T, D^{0}} > 1 GeV/#it{c}"));
    auto pt_text2 = pt->AddText(Form("#D^{0} = %.1f #pm %.1f", fitfunc->GetParameter(0), fitfunc->GetParError(0)));

    pt->Draw("SAME");

	 //    TF1 *signal = new TF1("signal", "gausn", 1.7, 2.1);
	 //    signal->SetParameter(0, fitfunc->GetParameter(0));
	 //    signal->SetParameter(1, fitfunc->GetParameter(1));
	 //    signal->SetParameter(2, fitfunc->GetParameter(2));
	 //    // // signal->SetParameter(3, fitfunc->GetParameter(3));

	 //    signal->SetParError(0, fitfunc->GetParError(0));
	 //    signal->SetParError(1, fitfunc->GetParError(1));
	 //    signal->SetParError(2, fitfunc->GetParError(2));

	 //    TCanvas *d = new TCanvas("d", "d", 500, 500);
		// d->cd();

		// TH1D *g1 = (TH1D *)s1.second;
		// TH1D *g2 = (TH1D *)s2.second;

		// // g2->Scale(usscalefactor/lsscalefactor);

		// g1->Add(g2, -1);

		// g1->Draw();

	 //    
}

void D0InvMass(){
	SaveTheD0MassHistograms();
	FitMethod();
}