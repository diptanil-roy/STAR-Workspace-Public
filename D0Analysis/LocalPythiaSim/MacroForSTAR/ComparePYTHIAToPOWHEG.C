#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include <TLorentzVector.h>

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

using namespace std;

void ComparePYTHIAToPOWHEG(){

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

	TFile *f = new TFile("st_physics_15094070_raw_0000007_3_0_NoD0Decay.root");
	TFile *g = new TFile("st_physics_15094070_raw_0000007_3_0_POWHEG.root");

	TH1F *JetPt[2];
	TH1F *D0Pt[2];
	TH1F *TrackPt[2];
	TH1F *Z[2];

	JetPt[0]   = (TH1F *)f->Get("hJetPt");
	D0Pt[0]    = (TH1F *)f->Get("hD0Pt");
	TrackPt[0] = (TH1F *)f->Get("hTrackPt");
	Z[0]       = (TH1F *)f->Get("histZ");

	JetPt[1]   = (TH1F *)g->Get("hJetPt");
	D0Pt[1]    = (TH1F *)g->Get("hD0Pt");
	TrackPt[1] = (TH1F *)g->Get("hTrackPt");
	Z[1]       = (TH1F *)g->Get("histZ");

	TCanvas *c[4];
	for (int i = 0; i < 4; i++){
		c[i] = new TCanvas(Form("c_%i", i), Form("c_%i", i), 800, 800);
	}

	c[0]->cd();
	JetPt[0]->SetNameTitle("PYTHIA", "PYTHIA");
	JetPt[1]->SetNameTitle("PYTHIA+POWHEG", "PYTHIA+POWHEG");
	JetPt[0]->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
	JetPt[0]->GetYaxis()->SetTitle("Counts");

	JetPt[0]->Rebin(2);
	JetPt[1]->Rebin(2);

	JetPt[0]->SetLineColor(1);
	JetPt[1]->SetLineColor(2);

	JetPt[0]->GetXaxis()->SetRangeUser(0, 40);
	JetPt[0]->GetXaxis()->SetNdivisions(505);
	gPad->SetLogy();

	JetPt[0]->Draw();
	JetPt[1]->Draw("SAME");
	gPad->BuildLegend();

	c[1]->cd();
	D0Pt[0]->SetNameTitle("PYTHIA", "PYTHIA");
	D0Pt[1]->SetNameTitle("PYTHIA+POWHEG", "PYTHIA+POWHEG");
	D0Pt[0]->GetXaxis()->SetTitle("p_{T,D^{0}} [GeV/#it{c}]");
	D0Pt[0]->GetYaxis()->SetTitle("Counts");

	D0Pt[0]->Rebin(2);
	D0Pt[1]->Rebin(2);

	D0Pt[0]->SetLineColor(1);
	D0Pt[1]->SetLineColor(2);

	D0Pt[0]->GetXaxis()->SetRangeUser(0, 20);
	D0Pt[0]->GetXaxis()->SetNdivisions(505);
	gPad->SetLogy();

	D0Pt[0]->Draw();
	D0Pt[1]->Draw("SAME");
	gPad->BuildLegend();

	c[2]->cd();
	TrackPt[0]->SetNameTitle("PYTHIA", "PYTHIA");
	TrackPt[1]->SetNameTitle("PYTHIA+POWHEG", "PYTHIA+POWHEG");
	TrackPt[0]->GetXaxis()->SetTitle("p_{T,Track} [GeV/#it{c}]");
	TrackPt[0]->GetYaxis()->SetTitle("Counts");

	TrackPt[0]->Rebin(2);
	TrackPt[1]->Rebin(2);

	TrackPt[0]->SetLineColor(1);
	TrackPt[1]->SetLineColor(2);

	TrackPt[0]->GetXaxis()->SetRangeUser(0, 20);
	TrackPt[0]->GetXaxis()->SetNdivisions(505);
	gPad->SetLogy();

	TrackPt[0]->Draw();
	TrackPt[1]->Draw("SAME");
	gPad->BuildLegend();

	c[3]->cd();
	Z[0]->SetNameTitle("PYTHIA", "PYTHIA");
	Z[1]->SetNameTitle("PYTHIA+POWHEG", "PYTHIA+POWHEG");
	Z[0]->GetXaxis()->SetTitle("Z");
	Z[0]->GetYaxis()->SetTitle("Counts");

	Z[0]->Rebin(2);
	Z[1]->Rebin(2);

	Z[0]->SetLineColor(1);
	Z[1]->SetLineColor(2);

	Z[0]->GetXaxis()->SetRangeUser(0, 1.2);
	Z[0]->GetXaxis()->SetNdivisions(505);
	gPad->SetLogy();

	Z[0]->Draw();
	Z[1]->Draw("SAME");
	gPad->BuildLegend();

}
