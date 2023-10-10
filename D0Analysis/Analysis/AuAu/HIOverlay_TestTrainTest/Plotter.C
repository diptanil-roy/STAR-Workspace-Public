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
// #include "StJetTreeStruct.h"
#include <vector>
using namespace std;
#endif

#include "BinDef.h"


void Plotter(){
	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

	int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
	int colors[200];

	for (int i = 0; i < 6; i++){
	    for (int j = -9; j <= 10; j++){
	        colors[i*20 + (j+9)] = col[i] + j; 
	    }
	}

	vector<TH1D *> Pt[3];
	vector<TH1D *> Z[3];

	TH1D *MeasuredPt[3];
	TH1D *MeasuredZ[3];

	TH1D *MCPt[3];
	TH1D *MCZ[3];

	TH2D *Response[3];


	TFile *Out = new TFile("OutFile.root");
	Out->cd();

	for (int iternum = 1; iternum <= 33; iternum++){

	  cout << "Iter # " << iternum << endl; 
	  Out->cd(Form("Iteration_%i", iternum*3));

	  for (int bin = 0; bin < 3; bin++){

	  	TH1D *h = (TH1D *)gDirectory->Get(Form("UnfoldedPt_%i", bin));

	  	Pt[bin].push_back((TH1D *)h->Rebin(nBinsJetPtForPlotting, Form("Iteration_%i", iternum*3), JetPtBinsForPlot));
	  	Z[bin].push_back((TH1D *)gDirectory->Get(Form("UnfoldedZ_%i", bin)));

	  	if (iternum == 1) MeasuredPt[bin] = (TH1D *)gDirectory->Get(Form("MeasuredPt_%i", bin));
	  	if (iternum == 1) MeasuredZ[bin] = (TH1D *)gDirectory->Get(Form("MeasuredZ_%i", bin));

	  	if (iternum == 1) MCPt[bin] = (TH1D *)((TH1D *)gDirectory->Get(Form("MCPt_%i", bin)))->Rebin(nBinsJetPtForPlotting, "True Pt", JetPtBinsForPlot);
	  	if (iternum == 1) MCZ[bin] = (TH1D *)gDirectory->Get(Form("MCZ_%i", bin));

	  	if (iternum == 1) Response[bin] = (TH2D *)gDirectory->Get(Form("Response2NewDef_%i", bin));
	  }
	}

	TCanvas *c[3];
	TCanvas *d[3];
	for (int i = 0; i < 3; i++){
		c[i] = new TCanvas(Form("c_%i", i), Form("c_%i", i), 1000, 1000);
		c[i]->cd();
		// gPad->SetLogy();

		for (int iternum = 1; iternum <= 33; iternum++){
			// Pt[i][iternum-1]->Rebin(nBinsJetPtForPlotting, Form("Iteration_%i", iternum*3), JetPtBinsForPlot);
			Pt[i][iternum-1]->Divide(MCPt[i]);
			// Pt[i][iternum-1]->GetYaxis()->SetRangeUser(7*pow(10, 1), 7*pow(10, 7));
			Pt[i][iternum-1]->GetYaxis()->SetRangeUser(0.1, 3.3);
			Pt[i][iternum-1]->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
			Pt[i][iternum-1]->GetYaxis()->SetTitle("#frac{Unfolded}{True}");
			Pt[i][iternum-1]->GetYaxis()->SetTitleOffset(0.9);
			Pt[i][iternum-1]->Draw("SAME");
			Pt[i][iternum-1]->SetLineColor(colors[iternum*3]);
			Pt[i][iternum-1]->SetMarkerColor(colors[iternum*3]);
			SetName(Pt[i][iternum-1], Form("Iteration_%i", iternum*3));
		}

		// SetName(MeasuredPt[i], "Measured Pt");
		// MeasuredPt[i]->SetLineWidth(3);
		// MeasuredPt[i]->SetLineColor(kBlack);
		// MeasuredPt[i]->SetMarkerColor(kBlack);
		// MeasuredPt[i]->Draw("HIST SAME");

		SetName(MCPt[i], "True Pt");
		MCPt[i]->Divide(MCPt[i]);
		MCPt[i]->SetLineWidth(3);
		MCPt[i]->SetLineColor(kRed);
		MCPt[i]->SetMarkerColor(kRed);
		MCPt[i]->Draw("HIST SAME");

		gPad->BuildLegend();

		SetName(MCZ[i], "True Z");
		MCZ[i]->Rebin(2);

		d[i] = new TCanvas(Form("d_%i", i), Form("d_%i", i), 1000, 1000);
		d[i]->cd();
		// gPad->SetLogy();

		for (int iternum = 1; iternum <= 33; iternum++){
			Z[i][iternum-1]->Rebin(2);
			Z[i][iternum-1]->Divide(MCZ[i]);
			Z[i][iternum-1]->Draw("SAME");
			// Z[i][iternum-1]->GetYaxis()->SetRangeUser(7*pow(10, 2), 7*pow(10, 5));
			Z[i][iternum-1]->GetXaxis()->SetTitle("Z");
			Z[i][iternum-1]->GetYaxis()->SetTitle("#frac{Unfolded}{True}");
			Z[i][iternum-1]->GetYaxis()->SetTitleOffset(0.9);
			Z[i][iternum-1]->GetYaxis()->SetRangeUser(0.1, 3.3);
			Z[i][iternum-1]->SetLineColor(colors[iternum*3]);
			Z[i][iternum-1]->SetMarkerColor(colors[iternum*3]);
			SetName(Z[i][iternum-1], Form("Iteration_%i", iternum*3));
		}

		// SetName(MeasuredZ[i], "Measured Z");
		// MeasuredZ[i]->SetLineWidth(3);
		// MeasuredZ[i]->SetLineColor(kBlack);
		// MeasuredZ[i]->SetMarkerColor(kBlack);
		// MeasuredZ[i]->Draw("HIST SAME");

		MCZ[i]->Divide(MCZ[i]);
		MCZ[i]->SetLineWidth(3);
		MCZ[i]->SetLineColor(kRed);
		MCZ[i]->SetMarkerColor(kRed);
		MCZ[i]->Draw("HIST SAME");

		gPad->BuildLegend();

	}

	TCanvas *e = new TCanvas("e", "e", 1000, 1000);
	e->Divide(3);
	for (int i = 0; i < 3; i++){
		e->cd(i+1);
		gPad->SetLogz();
		Response[i]->Draw("COLZ");
	}
}