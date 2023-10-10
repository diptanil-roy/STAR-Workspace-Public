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

#pragma link C++ class vector<int> +;

using namespace std;
// using namespace RooFit;

// #pragma link C++ class StJetTreeStruct+;

// #pragma link C++ class vector<float> +;
// #pragma link C++ class vector<vector<float> >+;
// #pragma link C++ class vector<int> +;
// #pragma link C++ class vector<vector<int> >+;
#endif


void CrossSectionMerger(){
	TH1::SetDefaultSumw2();
  	gStyle->SetOptStat(0);
  	gStyle->SetOptTitle(0);

	TFile *f[3];
	f[0] = new TFile("st_physics_15094070_raw_0000007_0_3_NoD0Decay.root");
	f[1] = new TFile("st_physics_15094070_raw_0000007_3_-1_NoD0Decay.root");
	f[2] = new TFile("st_physics_15094070_raw_0000007_0_-1_NoD0Decay.root");

	double crosssection[3] = {0.103, 0.006, 0.109};

	// double crosssection[3] = {0.0991233, 0.00589973, 0.10502303};

	TH1F *jetpt[3];
	TH1F *d0pt[3];

	for (int i = 0; i < 3; i++){
		jetpt[i] = (TH1F *)f[i]->Get("hJetPt");
		d0pt[i]  = (TH1F *)f[i]->Get("hD0Pt");

		// jetpt[i]->Scale(crosssection[i]/(jetpt[i]->Integral()*crosssection[2]));
		// d0pt[i]->Scale(crosssection[i]/(d0pt[i]->Integral()*crosssection[2]));

		jetpt[i]->Scale(crosssection[i]/(1000000.*crosssection[2]));
		d0pt[i]->Scale(crosssection[i]/(1000000.*crosssection[2]));

		// jetpt[i]->Scale(crosssection[i]/crosssection[2]);
		// d0pt[i]->Scale(crosssection[i]/crosssection[2]);
	}

	TH1F *mergedjetpt = (TH1F *)jetpt[0]->Clone("MergedJetPt");
	mergedjetpt->Add(jetpt[1]);
	mergedjetpt->SetTitle("Merged #hat{p}_{T} bins");
	jetpt[2]->SetTitle("No #hat{p}_{T} cut");

	TH1F *mergedd0pt = (TH1F *)d0pt[0]->Clone("MergedD0Pt");
	mergedd0pt->Add(d0pt[1]);
	mergedd0pt->SetTitle("Merged #hat{p}_{T} bins");
	d0pt[2]->SetTitle("No #hat{p}_{T} cut");


	TCanvas *c = new TCanvas("c", "c", 1900, 600);
	c->Divide(2);

	c->cd(1);
	gPad->SetLogy();

	mergedjetpt->Draw("EP");
	jetpt[2]->Draw("EP SAME");

	mergedjetpt->GetXaxis()->SetRangeUser(0, 30);
	mergedjetpt->GetXaxis()->SetTitle("p_{T,Jet}");
	mergedjetpt->GetYaxis()->SetTitle("#frac{1}{N_{events}} N_{Jets}");
	mergedjetpt->SetLineColor(kGreen - 2);
	mergedjetpt->SetMarkerColor(kGreen - 2);
	mergedjetpt->SetMarkerStyle(20);

	jetpt[2]->SetLineColor(kBlack);
	jetpt[2]->SetMarkerColor(kBlack);
	jetpt[2]->SetMarkerStyle(24);
	gPad->BuildLegend();

	c->cd(2);
	gPad->SetLogy();

	mergedd0pt->Draw("EP");
	d0pt[2]->Draw("EP SAME");

	mergedd0pt->GetXaxis()->SetRangeUser(0, 30);
	mergedd0pt->GetXaxis()->SetTitle("p_{T,D^{0}}");
	mergedd0pt->GetYaxis()->SetTitle("#frac{1}{N_{events}} N_{D^{0}}");
	mergedd0pt->SetLineColor(kGreen - 2);
	mergedd0pt->SetMarkerColor(kGreen - 2);
	mergedd0pt->SetMarkerStyle(20);

	d0pt[2]->SetLineColor(kBlack);
	d0pt[2]->SetMarkerColor(kBlack);
	d0pt[2]->SetMarkerStyle(24);
	gPad->BuildLegend();

	c->SaveAs("MergedPt.pdf");
}