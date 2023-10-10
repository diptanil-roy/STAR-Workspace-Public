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
#include "THnSparse.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
// #include "StJetTreeStruct.h"
#include <vector>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

// #include "RooUnfold/src/RooUnfoldResponse.h"
// #include "RooUnfold/src/RooUnfoldBayes.h"
// #include "RooUnfold/src/RooUnfoldSvd.h"

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


void Method(int SUPERITERATION = 1, int iteration = 0){

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);


	TFile *f = new TFile(Form("MCMCUnf/MCMCResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", SUPERITERATION, iteration, 200, 30));
	f->cd();
	THnSparseF *hJet[3];
	THnSparseF *hJetZNorm[3];
	THnSparseF *hJetZNormPtNorm[3];
	TH2D *hMCJetPtvRecoJetPt[3];

	// const char *opt = '0';

	for (int cent = 0; cent < 3; cent++){
		hJet[cent] = (THnSparseF *)f->Get(Form("hJet_%i", cent));
		hJetZNorm[cent] = (THnSparseF *)f->Get(Form("hJetZNorm_%i", cent));
		hJetZNormPtNorm[cent] = (THnSparseF *)f->Get(Form("hJetZNormPtNorm_%i", cent));

		cout << hJet[cent]->GetEntries() << "\t" << hJetZNorm[cent]->GetEntries() << "\t" << hJetZNormPtNorm[cent]->GetEntries() << endl;
		hMCJetPtvRecoJetPt[cent] = (TH2D *)hJetZNormPtNorm[cent]->Projection(2, 0);
	}

	TH1D *TruthZ[3];
  	TH1D *TruthPt[3];

	for (int cent = 0; cent < 3; cent++){
		TruthZ[cent] = (TH1D *)hJet[cent]->Projection(1, "E");
		TruthPt[cent] = (TH1D *)hJet[cent]->Projection(0, "E");
	}

	TH1D *TruthZZNorm[3];
	TH1D *TruthPtZNorm[3];

	for (int cent = 0; cent < 3; cent++){
		TruthZZNorm[cent] = (TH1D *)hJetZNorm[cent]->Projection(1, "E");
		TruthPtZNorm[cent] = (TH1D *)hJetZNorm[cent]->Projection(0, "E");
	}

	TH1D *TruthZZNormPtNorm[3];
	TH1D *TruthPtZNormPtNorm[3];

	for (int cent = 0; cent < 3; cent++){
		TruthZZNormPtNorm[cent] = (TH1D *)hJetZNormPtNorm[cent]->Projection(1, "E");
		TruthPtZNormPtNorm[cent] = (TH1D *)hJetZNormPtNorm[cent]->Projection(0, "E");

		cout << TruthPtZNormPtNorm[cent]->Integral() << endl;
	}

	TCanvas *d[3];

	for (int cent = 0; cent < 3; cent++){
		d[cent] = new TCanvas(Form("Testing the weighting mechanism %i", cent), Form("Testing the weighting mechanism %i", cent), 1200, 600);
		d[cent]->Divide(2);
		d[cent]->cd(1);
		gPad->SetLogy();
		SetAxisTitles(TruthPt[cent], "p_{T,Jet} [GeV/#it{c}]", "Normalised Yield");
		SetColor(TruthPt[cent], kBlue, 24);
		SetName(TruthPt[cent], "PYTHIA");
		SetColor(TruthPtZNorm[cent], kRed, 24);
		SetName(TruthPtZNorm[cent], "Z norm (2D)");
		SetColor(TruthPtZNormPtNorm[cent], kBlack, 20);
		SetName(TruthPtZNormPtNorm[cent], "Z norm + PYTHIA pT");
		double scale = TruthPt[cent]->Integral();
		TruthPt[cent]->Scale(1./TruthPt[cent]->Integral());
		TruthPtZNorm[cent]->Scale(1./TruthPtZNorm[cent]->Integral());
		TruthPtZNormPtNorm[cent]->Scale(1./TruthPtZNormPtNorm[cent]->Integral());
		TruthPt[cent]->Draw();
		TruthPtZNorm[cent]->Draw("SAME");
		TruthPtZNormPtNorm[cent]->Draw("SAME");
		gPad->BuildLegend();
		// d[cent]->Update();

		d[cent]->cd(2);
		gPad->SetLogy();
		SetAxisTitles(TruthZ[cent], "Z", "Normalised Yield");
		SetColor(TruthZ[cent], kBlue, 24);
		SetName(TruthZ[cent], "PYTHIA");
		SetColor(TruthZZNorm[cent], kRed, 24);
		SetName(TruthZZNorm[cent], "Z norm (2D)");
		SetColor(TruthZZNormPtNorm[cent], kBlack, 20);
		SetName(TruthZZNormPtNorm[cent], "Z norm + PYTHIA pT");
		scale = TruthZ[cent]->Integral();
		TruthZ[cent]->Scale(1./TruthZ[cent]->Integral());
		TruthZZNorm[cent]->Scale(1./TruthZZNorm[cent]->Integral());
		TruthZZNormPtNorm[cent]->Scale(1./TruthZZNormPtNorm[cent]->Integral());
		TruthZ[cent]->Draw();
		TruthZZNorm[cent]->Draw("SAME");
		TruthZZNormPtNorm[cent]->Draw("SAME");
		gPad->BuildLegend();

		// d[cent]->Update();
	}

	TCanvas *Ratio = new TCanvas("Ratio", "Ratio", 1200, 600);
	Ratio->Divide(2);

	TH1D *rTruthPt[3]; //Ratio
	TH1D *rTruthZ[3]; //Ratio

	for (int i = 0; i < 3; i++){
		// TruthPtZNormPtNorm[i]->Rebin(5);
		// TruthPt[i]->Rebin(5);

		rTruthPt[i] = (TH1D *)TruthPtZNormPtNorm[i]->Clone();
		rTruthPt[i]->Divide(TruthPt[i]);
		SetName(rTruthPt[i], Form("p_{T} Ratio %i", i));
		SetColor(rTruthPt[i], col[i*2]);

		rTruthZ[i] = (TH1D *)TruthZZNormPtNorm[i]->Clone();
		rTruthZ[i]->Divide(TruthZ[i]);
		SetName(rTruthZ[i], Form("Z Ratio %i", i));
		SetColor(rTruthZ[i], col[i*2]);
	}

	Ratio->cd(1);
	for (int i = 0; i < 3; i++){
		rTruthPt[i]->GetYaxis()->SetRangeUser(0, 2);
		SetAxisTitles(rTruthPt[i], "p_{T,Jet} [GeV/#it{c}]", "Ratio");
		rTruthPt[i]->Draw("SAME");
	}
	gPad->BuildLegend();

	Ratio->cd(2);
	for (int i = 0; i < 3; i++){
		rTruthZ[i]->Draw("SAME");
		SetAxisTitles(rTruthZ[i], "Z", "Ratio");
	}
	gPad->BuildLegend();
}

void Plot(){
	Method(0, 3);
}