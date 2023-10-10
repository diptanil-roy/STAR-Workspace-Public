#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TGaxis.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "THn.h"
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
#include "TUnfold.h"
#include "TRandom3.h"
#include "TLegend.h"
#include <algorithm>

#include <cstdarg>

TString yaxisforratioplots = "#frac{Unfolded}{PYTHIA 8}";

TString xaxisforptspectra = "p_{T,jet} [GeV/#it{c}]";
TString yaxisforptspectra = "#frac{dN_{jet}}{dp_{T,jet}} [GeV/#it{c}]^{-1}";


TString xaxisfordRspectra = "#Delta#it{R}";
TString yaxisfordRspectra = "#frac{1}{N_{jet}} #frac{#DeltaN_{jet}}{#Delta#it{R}}";


void Spectra(){

	// std::vector<pair<int, int>> v;

	// pair<int, int> s1(1,2);
	// v.push_back(s1);

	// pair<int, int> s2(2,3);
	// v.push_back(s2);

	// cout << v.size() << endl;

	gROOT->ProcessLine(".x ./myStyle.C");

	TFile *f = new TFile(Form("../Response2022/UnfoldedPlots_Closure_%f.root", 5.0));
	TH1D *CentralJetPt = (TH1D *)f->Get("Central_pT_Closure/UnfoldedJetPt");
	TH1D *MidCentralJetPt = (TH1D *)f->Get("MidCentral_pT_Closure/UnfoldedJetPt");
	TH1D *PeripheralJetPt = (TH1D *)f->Get("Peripheral_pT_Closure/UnfoldedJetPt");

	CentralJetPt->Scale(1./CentralJetPt->Integral());
	MidCentralJetPt->Scale(1./MidCentralJetPt->Integral());
	PeripheralJetPt->Scale(1./PeripheralJetPt->Integral());

	CentralJetPt->SetNameTitle("Central/Peripheral", "Central/Peripheral");
	MidCentralJetPt->SetNameTitle("MidCentral/Peripheral", "MidCentral/Peripheral");

	CentralJetPt->Divide(PeripheralJetPt);
	MidCentralJetPt->Divide(PeripheralJetPt);

	CentralJetPt->SetLineColor(kRed);
	CentralJetPt->SetMarkerColor(kRed);
	CentralJetPt->SetMarkerStyle(20);

	MidCentralJetPt->SetLineColor(kBlue);
	MidCentralJetPt->SetMarkerColor(kBlue);
	MidCentralJetPt->SetMarkerStyle(21);

	CentralJetPt->Draw();
	MidCentralJetPt->Draw("SAME");

	gPad->BuildLegend();


}

void DeltaR(){

	// std::vector<pair<int, int>> v;

	// pair<int, int> s1(1,2);
	// v.push_back(s1);

	// pair<int, int> s2(2,3);
	// v.push_back(s2);

	// cout << v.size() << endl;

	gROOT->ProcessLine(".x ./myStyle.C");

	TFile *f = new TFile(Form("../Response2022/UnfoldedPlots_Closure_%f.root", 5.0));
	TH1D *CentralJetPt = (TH1D *)f->Get("Central_deltaR_Closure_2D/UnfoldedJetdR");
	TH1D *MidCentralJetPt = (TH1D *)f->Get("MidCentral_deltaR_Closure_2D/UnfoldedJetdR");
	TH1D *PeripheralJetPt = (TH1D *)f->Get("Peripheral_deltaR_Closure_2D/UnfoldedJetdR");

	// CentralJetPt->Scale(1./CentralJetPt->Integral());
	// MidCentralJetPt->Scale(1./MidCentralJetPt->Integral());
	// PeripheralJetPt->Scale(1./PeripheralJetPt->Integral());

	CentralJetPt->SetNameTitle("Central/Peripheral", "Central/Peripheral");
	MidCentralJetPt->SetNameTitle("MidCentral/Peripheral", "MidCentral/Peripheral");

	CentralJetPt->Divide(PeripheralJetPt);
	MidCentralJetPt->Divide(PeripheralJetPt);

	CentralJetPt->SetLineColor(kRed);
	CentralJetPt->SetMarkerColor(kRed);
	CentralJetPt->SetMarkerStyle(20);

	MidCentralJetPt->SetLineColor(kBlue);
	MidCentralJetPt->SetMarkerColor(kBlue);
	MidCentralJetPt->SetMarkerStyle(21);

	CentralJetPt->Draw();
	MidCentralJetPt->Draw("SAME");

	gPad->BuildLegend();


}


void RatioOfSpectra(){
	Spectra();
	// DeltaR();
}


