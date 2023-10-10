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

#include "../BinDef.h"
#include "NewBinDef.h"

void MakeMCComparisonPlots(int iteration = 3, TString DirName = "MCMCUnf", int NSUPERITERATIONS = 50){

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	// gStyle->SetOptTitle(0);

	// const int NSUPERITERATIONS = 50;

	TH2D *UnmodifiedMeasured[3];

	TH1D *UnmodifiedMeasuredPt[3];
	TH1D *UnmodifiedMeasuredZ[3];


	TH2D *Truth[3];
	TH2D *Measured[3][NSUPERITERATIONS];
	TH2D *Unfolded[3][NSUPERITERATIONS];

	TH1D *TruthPt[3];
	TH1D *TruthZ[3];

	TH1D *MeasuredPt[3][NSUPERITERATIONS];
	TH1D *MeasuredZ[3][NSUPERITERATIONS];


	TH1D *UnfoldedPt[3][NSUPERITERATIONS];
	TH1D *UnfoldedZ[3][NSUPERITERATIONS];

	TCanvas *PlottedCanvas[3][NSUPERITERATIONS];

	TString DirectoryName[3] = {"Central", "MidCentral", "Peripheral"};
	
	for (int SUPERITERATION = 1; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
		TFile *f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins));
		cout << f->GetName() << endl;
		f->cd();

		for (int cent = 0; cent < 3; cent++){
			if (SUPERITERATION == 1) {
				Truth[cent] = (TH2D *)f->Get(Form("Truth_%i", cent));
				Truth[cent]->SetDirectory(0);
				TruthPt[cent] = (TH1D *)Truth[cent]->ProjectionX();
				TruthZ[cent] = (TH1D *)Truth[cent]->ProjectionY();
			}

			Measured[cent][SUPERITERATION-1] = (TH2D *)f->Get(Form("Measured_%i", cent));
			SetName(Measured[cent][SUPERITERATION-1], Form("Measured_%i_%i", SUPERITERATION, cent));
			Measured[cent][SUPERITERATION-1]->SetDirectory(0);
			MeasuredPt[cent][SUPERITERATION-1] = (TH1D *)Measured[cent][SUPERITERATION-1]->ProjectionX();
			MeasuredZ[cent][SUPERITERATION-1] = (TH1D *)Measured[cent][SUPERITERATION-1]->ProjectionY();
			MeasuredPt[cent][SUPERITERATION-1]->SetDirectory(0);
			MeasuredZ[cent][SUPERITERATION-1]->SetDirectory(0);
			Unfolded[cent][SUPERITERATION-1] = (TH2D *)f->Get(Form("Unfolded_%i", cent));
			SetName(Unfolded[cent][SUPERITERATION-1], Form("Unfolded_%i_%i", SUPERITERATION, cent));
			Unfolded[cent][SUPERITERATION-1]->SetDirectory(0);
			UnfoldedPt[cent][SUPERITERATION-1] = (TH1D *)Unfolded[cent][SUPERITERATION-1]->ProjectionX();
			UnfoldedZ[cent][SUPERITERATION-1] = (TH1D *)Unfolded[cent][SUPERITERATION-1]->ProjectionY();
			UnfoldedPt[cent][SUPERITERATION-1]->SetDirectory(0);
			UnfoldedZ[cent][SUPERITERATION-1]->SetDirectory(0);

			PlottedCanvas[cent][SUPERITERATION-1] = (TCanvas *)f->Get(Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, cent));
			PlottedCanvas[cent][SUPERITERATION-1]->Draw();
			// // PlottedCanvas[cent][SUPERITERATION-1]->SetDirectory(0);
			PlottedCanvas[cent][SUPERITERATION-1]->SetCanvasSize(2000, 1000);
			PlottedCanvas[cent][SUPERITERATION-1]->SetTitle(Form("SUPERITERATION #%i", SUPERITERATION));

			TPad *pad5 = new TPad("all","all",0,0,1,1);
			pad5->SetFillStyle(4000);  // transparent
			pad5->Draw();
			pad5->cd();
			TLatex *lat = new TLatex();
			lat->DrawLatexNDC(.4,.95,Form("SUPERITERATION #%i", SUPERITERATION));

			PlottedCanvas[cent][SUPERITERATION-1]->Update();
			PlottedCanvas[cent][SUPERITERATION-1]->SaveAs(Form("%s/%s/Plots_Step_%i_Iter_%i_Cent_%i.pdf", DirName.Data(), DirectoryName[cent].Data(), SUPERITERATION, iteration, cent));
			PlottedCanvas[cent][SUPERITERATION-1]->SaveAs(Form("%s/%s/Plots_Step_%i_Iter_%i_Cent_%i.png", DirName.Data(), DirectoryName[cent].Data(), SUPERITERATION, iteration, cent));
			delete PlottedCanvas[cent][SUPERITERATION-1];
		}

		f->Close();

		cout << Unfolded[0][SUPERITERATION-1]->Integral() << endl;
	}

	TFile *g = new TFile(Form("%s/FinalResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, iteration, njpt_gen_bins, nz_gen_bins));
	cout << g->GetName() << endl;
	g->cd();
	// g->ls();

	for (int cent = 0; cent < 3; cent++){
		UnmodifiedMeasured[cent] = (TH2D *)g->Get(Form("Measured_%i", cent));
		SetName(UnmodifiedMeasured[cent], Form("Measured_%i_%i", 0, cent));
		UnmodifiedMeasured[cent]->SetDirectory(0);
		UnmodifiedMeasuredPt[cent] = (TH1D *)UnmodifiedMeasured[cent]->ProjectionX();
		UnmodifiedMeasuredZ[cent] = (TH1D *)UnmodifiedMeasured[cent]->ProjectionY();
		UnmodifiedMeasuredPt[cent]->SetDirectory(0);
		UnmodifiedMeasuredZ[cent]->SetDirectory(0);
	}

	cout << "TAA Called" << endl;

	TH1D *RatioPT[NSUPERITERATIONS+1][3];
	TH1D *RatioZ[NSUPERITERATIONS+1][3];

	for (int SUPERITERATION = 1; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
		for (int cent = 0; cent < 3; cent++){

			// cout << "Cent = " << cent << endl;
			RatioPT[SUPERITERATION-1][cent] = (TH1D *)UnfoldedPt[cent][SUPERITERATION-1]->Clone(Form("RCP_pT_%i_%i", SUPERITERATION, cent));
			RatioPT[SUPERITERATION-1][cent]->Scale(1./RatioPT[SUPERITERATION-1][cent]->Integral());
			RatioZ[SUPERITERATION-1][cent] = (TH1D *)UnfoldedZ[cent][SUPERITERATION-1]->Clone(Form("RCP_Z_%i_%i", SUPERITERATION, cent));
			RatioZ[SUPERITERATION-1][cent]->Scale(1./RatioZ[SUPERITERATION-1][cent]->Integral());
		}
		RatioPT[SUPERITERATION-1][0]->Divide(RatioPT[SUPERITERATION-1][2]);
		RatioZ[SUPERITERATION-1][0]->Divide(RatioZ[SUPERITERATION-1][2]);

		RatioPT[SUPERITERATION-1][1]->Divide(RatioPT[SUPERITERATION-1][2]);
		RatioZ[SUPERITERATION-1][1]->Divide(RatioZ[SUPERITERATION-1][2]);
	}
	

	// g->Close();

	int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
	int colors[12];

	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 2; j++){
			if (j == 0) colors[i*2 + j] = col[i] - 5;
			else if (j == 1) colors[i*2 + j] = col[i] + 5;
		}
	}


	TCanvas *d = new TCanvas ("RCP", "RCP", 2000, 1000);
	d->Divide(2);
	for (int SUPERITERATION = 1; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
		d->cd(1);
		SetColor(RatioPT[SUPERITERATION-1][0], colors[SUPERITERATION-1]);
		RatioPT[SUPERITERATION-1][0]->GetYaxis()->SetRangeUser(0,2);
		RatioPT[SUPERITERATION-1][0]->Draw("EP SAME");
		
		d->cd(2);
		SetColor(RatioZ[SUPERITERATION-1][0], colors[SUPERITERATION-1]);
		RatioZ[SUPERITERATION-1][0]->Draw("EP SAME");
		RatioZ[SUPERITERATION-1][0]->GetYaxis()->SetRangeUser(0,2);
	}
	d->cd(1);
	gPad->BuildLegend(0.55,0.7,0.93,0.9);
	d->cd(2);
	gPad->BuildLegend(0.55,0.7,0.93,0.9);

	d->SaveAs(Form("%s/RCP.pdf", DirName.Data()));

	int color[3] = {kGreen-2, kBlue, kBlack};

	TH1D *Chi2[3];
	for (int cent = 0; cent < 3; cent++){
		Chi2[cent] = new TH1D(Form("Chi2_%i", cent), Form("Chi2_%i", cent), NSUPERITERATIONS, 0.5, NSUPERITERATIONS+0.5);

		for (int SUPERITERATION = 1; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
			Chi2[cent]->SetBinContent(SUPERITERATION, CalculateChi2(Unfolded[cent][SUPERITERATION-1], Truth[cent]));
		}
		SetColor(Chi2[cent], color[cent], 20 + cent);
	}

	TCanvas *c[5];
	for (int i = 0; i < 5; i++){
		c[i] = new TCanvas(Form("Canvas_%i", i), Form("Canvas_%i", i), 2000, 1000);
	}
	
	c[0]->cd();
	gPad->SetLogy();
	for (int cent = 0; cent < 3; cent++){
		Chi2[cent]->Draw("HIST SAME");
	}

	c[0]->SaveAs(Form("%s/Chi2.pdf", DirName.Data()));

	for (int siter = 1; siter <= 3; siter++){
		c[siter]->Divide(2);
		c[siter]->cd(1);
		gPad->SetLogy();
		for (int cent = 0; cent < 3; cent++){

			MeasuredPt[cent][siter-1]->Divide(UnmodifiedMeasuredPt[cent]);
			SetColor(MeasuredPt[cent][siter-1], color[cent], 20);
			MeasuredPt[cent][siter-1]->Draw("SAME");
			

		}
		
		c[siter]->cd(2);
		gPad->SetLogy();
		for (int cent = 0; cent < 3; cent++){
			MeasuredZ[cent][siter-1]->Divide(UnmodifiedMeasuredZ[cent]);
			SetColor(MeasuredZ[cent][siter-1], color[cent], 20);
			MeasuredZ[cent][siter-1]->Draw("SAME");
			
		}
		

		c[siter]->SaveAs(Form("%s/MeasuredDistComparison_Step%i_Iter_%i.pdf", DirName.Data(), siter, iteration));
	}
}
