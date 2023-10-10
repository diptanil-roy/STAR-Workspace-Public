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

void MakeDataComparisonPlots(int iteration = 3, TString DirName = "DataMCUnf", bool oldresp = kFALSE, int NSUPERITERATIONS = 50){

	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	// gStyle->SetOptTitle(0);


	TH2D *Truth[3];
	TH2D *Unfolded[3][NSUPERITERATIONS+1];

	TH1D *TruthPt[3];
	TH1D *TruthZ[3];

	TH1D *UnfoldedPt[3][NSUPERITERATIONS+1];
	TH1D *UnfoldedZ[3][NSUPERITERATIONS+1];

	TCanvas *PlottedCanvas[3][NSUPERITERATIONS+1];

	TString DirectoryName[3] = {"Central", "MidCentral", "Peripheral"};
	
	for (int SUPERITERATION = 0; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
		TFile *f;
		if (oldresp)f = new TFile(Form("%s/OldOutput_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins));
		else f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins));
		cout << f->GetName() << endl;
		f->cd();

		for (int cent = 0; cent < 3; cent++){
			// if (SUPERITERATION == 1) {
			// 	Truth[cent] = (TH2D *)f->Get(Form("Truth_%i", cent));
			// 	Truth[cent]->SetDirectory(0);
			// 	TruthPt[cent] = (TH1D *)Truth[cent]->ProjectionX();
			// 	TruthZ[cent] = (TH1D *)Truth[cent]->ProjectionY();
			// }
			Unfolded[cent][SUPERITERATION] = (TH2D *)f->Get(Form("Unfolded_%i", cent));
			SetName(Unfolded[cent][SUPERITERATION], Form("Unfolded_%i_%i", SUPERITERATION, cent));
			Unfolded[cent][SUPERITERATION]->SetDirectory(0);
			UnfoldedPt[cent][SUPERITERATION] = (TH1D *)Unfolded[cent][SUPERITERATION]->ProjectionX();
			UnfoldedPt[cent][SUPERITERATION]->SetDirectory(0);
			UnfoldedZ[cent][SUPERITERATION] = (TH1D *)Unfolded[cent][SUPERITERATION]->ProjectionY();
			UnfoldedZ[cent][SUPERITERATION]->SetDirectory(0);

			PlottedCanvas[cent][SUPERITERATION] = (TCanvas *)f->Get(Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, cent));
			PlottedCanvas[cent][SUPERITERATION]->Draw();
			// // PlottedCanvas[cent][SUPERITERATION]->SetDirectory(0);
			PlottedCanvas[cent][SUPERITERATION]->SetCanvasSize(2000, 1000);
			PlottedCanvas[cent][SUPERITERATION]->SetTitle(Form("SUPERITERATION #%i", SUPERITERATION));

			TPad *pad5 = new TPad("all","all",0,0,1,1);
			pad5->SetFillStyle(4000);  // transparent
			pad5->Draw();
			pad5->cd();
			TLatex *lat = new TLatex();
			lat->DrawLatexNDC(.4,.95,Form("SUPERITERATION #%i", SUPERITERATION));

			PlottedCanvas[cent][SUPERITERATION]->Update();
			PlottedCanvas[cent][SUPERITERATION]->SaveAs(Form("%s/%s/Plots_Step_%i_Iter_%i_Cent_%i.pdf", DirName.Data(), DirectoryName[cent].Data(), SUPERITERATION, iteration, cent));
			PlottedCanvas[cent][SUPERITERATION]->SaveAs(Form("%s/%s/Plots_Step_%i_Iter_%i_Cent_%i.png", DirName.Data(), DirectoryName[cent].Data(), SUPERITERATION, iteration, cent));
			delete PlottedCanvas[cent][SUPERITERATION];
		}

		f->Close();

		cout << Unfolded[0][SUPERITERATION]->Integral() << endl;
	}

	int color[3] = {kGreen-2, kBlue, kBlack};

	double taa[3] = {941.23714*1.0318440e+08, 391.35550*3.2123506e+08, 56.62475*4.6679240e+08};

	cout << "TAA Called" << endl;

	TH1D *RatioPT[NSUPERITERATIONS+1][3];
	TH1D *RatioZ[NSUPERITERATIONS+1][3];

	TH1D *tmpPT[NSUPERITERATIONS+1][3];
	TH1D *tmpZ[NSUPERITERATIONS+1][3];

	for (int SUPERITERATION = 0; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
		for (int cent = 0; cent < 3; cent++){

			tmpPT[SUPERITERATION][cent] = (TH1D *)UnfoldedPt[cent][SUPERITERATION]->Clone(Form("tmp_pT_%i_%i", SUPERITERATION, cent));
			// tmpPT[SUPERITERATION][cent]->Scale(1./tmpPT[SUPERITERATION][cent]->Integral());
			// tmpPT[SUPERITERATION][cent]->Scale(1./taa[cent]);
			tmpZ[SUPERITERATION][cent] = (TH1D *)UnfoldedZ[cent][SUPERITERATION]->Clone(Form("tmp_Z_%i_%i", SUPERITERATION, cent));
			// tmpZ[SUPERITERATION][cent]->Scale(1./tmpZ[SUPERITERATION][cent]->Integral());
			// tmpZ[SUPERITERATION][cent]->Scale(1./taa[cent]);

			// cout << "Cent = " << cent << endl;
			RatioPT[SUPERITERATION][cent] = (TH1D *)UnfoldedPt[cent][SUPERITERATION]->Clone(Form("RCP_pT_%i_%i", SUPERITERATION, cent));
			RatioPT[SUPERITERATION][cent]->Scale(1./taa[cent]);
			RatioZ[SUPERITERATION][cent] = (TH1D *)UnfoldedZ[cent][SUPERITERATION]->Clone(Form("RCP_Z_%i_%i", SUPERITERATION, cent));
			RatioZ[SUPERITERATION][cent]->Scale(1./taa[cent]);
		}
		RatioPT[SUPERITERATION][0]->Divide(RatioPT[SUPERITERATION][2]);
		RatioZ[SUPERITERATION][0]->Divide(RatioZ[SUPERITERATION][2]);

		RatioPT[SUPERITERATION][1]->Divide(RatioPT[SUPERITERATION][2]);
		RatioZ[SUPERITERATION][1]->Divide(RatioZ[SUPERITERATION][2]);
	}

	cout << "Ratio Plots Made" << endl;

	int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
	int colors[12];

	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 2; j++){
			if (j == 0) colors[i*2 + j] = col[i] - 5;
			else if (j == 1) colors[i*2 + j] = col[i] + 5;
		}
	}

	int centcolor[3] = {kGreen-2, kBlue, kBlack};

	TCanvas *c = new TCanvas ("RCP", "RCP", 2000, 1000);
	c->Divide(2);
	for (int SUPERITERATION = 0; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
		c->cd(1);
		SetColor(RatioPT[SUPERITERATION][0], colors[SUPERITERATION]);
		RatioPT[SUPERITERATION][0]->GetYaxis()->SetRangeUser(0,1);
		RatioPT[SUPERITERATION][0]->Draw("EP SAME");
		
		c->cd(2);
		SetColor(RatioZ[SUPERITERATION][0], colors[SUPERITERATION]);
		RatioZ[SUPERITERATION][0]->Draw("EP SAME");
		RatioZ[SUPERITERATION][0]->GetYaxis()->SetRangeUser(0,5);
	}
	c->cd(1);
	gPad->BuildLegend(0.55,0.7,0.93,0.9);
	c->cd(2);
	gPad->BuildLegend(0.55,0.7,0.93,0.9);

	c->SaveAs(Form("%s/RCP.pdf", DirName.Data()));

	TCanvas *d[NSUPERITERATIONS+1];

	for (int SUPERITERATION = 0; SUPERITERATION <= NSUPERITERATIONS; SUPERITERATION++){
		d[SUPERITERATION] = new TCanvas (Form("Pt and Z %i", SUPERITERATION), Form("Pt and Z %i", SUPERITERATION), 2000, 1000);
		d[SUPERITERATION]->Divide(2);

		for (int cent = 2; cent >= 0; cent--){
			d[SUPERITERATION]->cd(1);
			gPad->SetLogy();

			SetColor(tmpPT[SUPERITERATION][cent], centcolor[cent]);
			tmpPT[SUPERITERATION][cent]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 6));
			tmpPT[SUPERITERATION][cent]->Draw("EP SAME");
			
			d[SUPERITERATION]->cd(2);
			gPad->SetLogy();
			SetColor(tmpZ[SUPERITERATION][cent], centcolor[cent]);
			tmpZ[SUPERITERATION][cent]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 6));
			tmpZ[SUPERITERATION][cent]->Draw("EP SAME");
		}

		d[SUPERITERATION]->cd(1);
		gPad->BuildLegend(0.55,0.7,0.93,0.9);
		d[SUPERITERATION]->cd(2);
		gPad->BuildLegend(0.55,0.7,0.93,0.9);

		d[SUPERITERATION]->SaveAs(Form("%s/PTandZ_Step_%i_Iter_%i.pdf", DirName.Data(), SUPERITERATION, iteration));
	}
	

}
