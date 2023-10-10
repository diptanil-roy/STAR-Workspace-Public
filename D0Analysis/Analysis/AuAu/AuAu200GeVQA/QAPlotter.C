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
using namespace std;
#endif

#include "../BinDef.h"

void Method1(){
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	// gStyle->SetOptTitle(0);
	TFile *f = new TFile("/Volumes/WorkDrive/AuAu200GeVQA/outSL18f/AuAu200GeVQA_SL18f.root");
	TFile *g = new TFile("/Volumes/WorkDrive/AuAu200GeVQA/out/AuAu200GeVQA_SL22c.root");

	f->cd("TrackClusterQA");
	TH1F *Centrality = (TH1F *)gDirectory->Get("fHistCentrality");
	TH1F *Multiplicity   = (TH1F *)gDirectory->Get("fHistMultiplicity");

	TH2F *dEdxVp = (TH2F *)gDirectory->Get("dEdXvp");
	TH2F *m2vp = (TH2F *)gDirectory->Get("m2vp");

	TH1F *TrackPt = (TH1F *)gDirectory->Get("fHistNTrackvsPt");
	TH1F *TrackEta = (TH1F *)gDirectory->Get("fHistNTrackvsEta");
	TH1F *TrackPhi = (TH1F *)gDirectory->Get("fHistNTrackvsPhi");

	TrackPt->Rebin(2);
	TrackPt->Scale(1./Centrality->Integral());
	TrackEta->Scale(1./Centrality->Integral());
	TrackPhi->Scale(1./Centrality->Integral());

	TH1F *TrackHFTPt = (TH1F *)gDirectory->Get("fHistHFTNTrackvsPt");
	TH1F *TrackHFTEta = (TH1F *)gDirectory->Get("fHistHFTNTrackvsEta");
	TH1F *TrackHFTPhi = (TH1F *)gDirectory->Get("fHistHFTNTrackvsPhi");

	TrackHFTPt->Rebin(2);
	TrackHFTPt->Scale(1./Centrality->Integral());
	TrackHFTEta->Scale(1./Centrality->Integral());
	TrackHFTPhi->Scale(1./Centrality->Integral());

	TH2F *PicoTowerVsProjectedTower = (TH2F *)gDirectory->Get("PicoVsProjectedTowerID");
	TH2F *PicoTowerVsProjectedTowerPt[6];
	for (int i = 0; i < 6; i++){
		PicoTowerVsProjectedTowerPt[i] = (TH2F *)gDirectory->Get(Form("PicoVsProjectedTowerID_%i", i));
	}

	TH1D *TowerE = (TH1D *)gDirectory->Get("fHistNTowervsE");
	TH1D *TowerPhi = (TH1D *)gDirectory->Get("fHistNTowervsPhi");

	f->cd("D0Tagger");
	TH2D *VzVsVzVPD = (TH2D *)gDirectory->Get("hEventVzvsVzvpd");
	TH2D *VxVsVy = (TH2D *)gDirectory->Get("hEventVxvsVy");

	g->cd("TrackClusterQA");
	TH1F *CentralityN = (TH1F *)gDirectory->Get("fHistCentrality");
	TH1F *MultiplicityN   = (TH1F *)gDirectory->Get("fHistMultiplicity");

	TH2F *dEdxVpN = (TH2F *)gDirectory->Get("dEdXvp");
	TH2F *m2vpN = (TH2F *)gDirectory->Get("m2vp");

	TH1F *TrackPtN = (TH1F *)gDirectory->Get("fHistNTrackvsPt");
	TH1F *TrackEtaN = (TH1F *)gDirectory->Get("fHistNTrackvsEta");
	TH1F *TrackPhiN = (TH1F *)gDirectory->Get("fHistNTrackvsPhi");

	TrackPtN->Rebin(2);
	TrackPtN->Scale(1./CentralityN->Integral());
	TrackEtaN->Scale(1./CentralityN->Integral());
	TrackPhiN->Scale(1./CentralityN->Integral());

	TH1F *TrackHFTPtN = (TH1F *)gDirectory->Get("fHistHFTNTrackvsPt");
	TH1F *TrackHFTEtaN = (TH1F *)gDirectory->Get("fHistHFTNTrackvsEta");
	TH1F *TrackHFTPhiN = (TH1F *)gDirectory->Get("fHistHFTNTrackvsPhi");

	TrackHFTPtN->Rebin(2);
	TrackHFTPtN->Scale(1./CentralityN->Integral());
	TrackHFTEtaN->Scale(1./CentralityN->Integral());
	TrackHFTPhiN->Scale(1./CentralityN->Integral());

	TH2F *PicoTowerVsProjectedTowerN = (TH2F *)gDirectory->Get("PicoVsProjectedTowerID");
	TH2F *PicoTowerVsProjectedTowerPtN[6];
	for (int i = 0; i < 6; i++){
		PicoTowerVsProjectedTowerPtN[i] = (TH2F *)gDirectory->Get(Form("PicoVsProjectedTowerID_%i", i));
	}

	TH1D *TowerEN = (TH1D *)gDirectory->Get("fHistNTowervsE");
	TH1D *TowerPhiN = (TH1D *)gDirectory->Get("fHistNTowervsPhi");


	g->cd("D0Tagger");
	TH2D *VzVsVzVPDN = (TH2D *)gDirectory->Get("hEventVzvsVzvpd");
	TH2D *VxVsVyN = (TH2D *)gDirectory->Get("hEventVxvsVy");


	TCanvas *c[12];
	for (int i = 0; i < 12; i++){
		c[i] = new TCanvas(Form("old_%i", i), Form("old_%i", i), 1000, 1000);
	}


	TCanvas *d[15];
	for (int i = 0; i < 15; i++){
		d[i] = new TCanvas(Form("new_%i", i), Form("new_%i", i), 1000, 1000);
	}

	TH1D *old = new TH1D("SL18f", "SL18f", 10, 0, 1);
	SetColor(old, kBlack, 20);
	TH1D *newlib = new TH1D("SL22c", "SL22c", 10, 0, 1);
	SetColor(newlib, kRed, 24);

	auto legend = new TLegend(0.5,0.7,0.88,0.9);
	legend->AddEntry(old,"SL18f","lep");
	legend->AddEntry(newlib,"SL22c","lep");


	c[0]->cd();
	// gPad->SetLogy();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	m2vp->Draw("COLZ");
	m2vp->RebinX(10);
	m2vp->RebinY(5);
	m2vp->Scale(1./Centrality->Integral());
	m2vp->GetXaxis()->SetRangeUser(0, 5);
	m2vp->GetYaxis()->SetRangeUser(0, 2);
	m2vp->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
	m2vp->GetYaxis()->SetTitle("m^{2} [GeV^{2}]");
	m2vp->GetXaxis()->SetTitleOffset(1.2);
	m2vp->GetYaxis()->SetTitleOffset(1.3);
	SetName(m2vp, "m^{2} vs p_{T} from TOF (SL18f)");

	c[0]->SaveAs("QAPlots/m2vspT_SL18f.pdf");

	d[0]->cd();
	// gPad->SetLogy();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	m2vpN->Draw("COLZ");
	m2vpN->RebinX(10);
	m2vpN->RebinY(5);
	m2vpN->Scale(1./CentralityN->Integral());
	m2vpN->GetXaxis()->SetRangeUser(0, 5);
	m2vpN->GetYaxis()->SetRangeUser(0, 2);
	m2vpN->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
	m2vpN->GetYaxis()->SetTitle("m^{2} [GeV^{2}]");
	m2vpN->GetXaxis()->SetTitleOffset(1.2);
	m2vpN->GetYaxis()->SetTitleOffset(1.3);
	SetName(m2vpN, "m^{2} vs p_{T} from TOF (SL22c)");

	d[0]->SaveAs("QAPlots/m2vspT_SL22c.pdf");

	c[1]->cd();
	TH1D *m2 = (TH1D *)m2vp->ProjectionY();
	TH1D *m2N = (TH1D *)m2vpN->ProjectionY();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	SetName(m2, "m^{2} from TOF (SL18f)");
	SetColor(m2, kBlack, 20);
	SetName(m2N, "m^{2} from TOF (SL22c)");
	SetColor(m2N, kRed, 24);
	m2->Draw("E SAME");
	m2N->Draw("E SAME");
	
	legend->Draw("SAME");

	c[1]->SaveAs("QAPlots/m2_SL18fvSL22c.pdf");

	c[2]->cd();
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	dEdxVp->Draw("COLZ");
	dEdxVp->RebinX(2);
	// dEdxVp->RebinY(5);
	dEdxVp->Scale(1./Centrality->Integral());
	dEdxVp->GetXaxis()->SetRangeUser(0, 5);
	dEdxVp->GetYaxis()->SetRangeUser(0.5, 10);
	dEdxVp->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
	dEdxVp->GetYaxis()->SetTitle("#frac{dE}{dx} [GeV/cm]");
	dEdxVp->GetXaxis()->SetTitleOffset(1.2);
	dEdxVp->GetYaxis()->SetTitleOffset(1.3);
	dEdxVp->SetTitle("dE/dx vs p_{T} from TPC");

	c[2]->SaveAs("QAPlots/dEdx_SL18f.pdf");

	d[1]->cd();
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	dEdxVpN->Draw("COLZ");
	dEdxVpN->RebinX(2);
	// dEdxVp->RebinY(5);
	dEdxVpN->Scale(1./CentralityN->Integral());
	dEdxVpN->GetXaxis()->SetRangeUser(0, 5);
	dEdxVpN->GetYaxis()->SetRangeUser(0.5, 10);
	dEdxVpN->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
	dEdxVpN->GetYaxis()->SetTitle("#frac{dE}{dx} [GeV/cm]");
	dEdxVpN->GetXaxis()->SetTitleOffset(1.2);
	dEdxVpN->GetYaxis()->SetTitleOffset(1.3);
	dEdxVpN->SetTitle("dE/dx vs p_{T} from TPC");

	d[1]->SaveAs("QAPlots/dEdx_SL22c.pdf");


	c[3]->cd();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TH1D *RatioPt = (TH1D *)TrackHFTPt->Clone("HFT/TPC p_{T}");
	RatioPt->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
	RatioPt->GetYaxis()->SetTitle("#frac{HFT}{TPC}");
	RatioPt->Divide(TrackPt);
	RatioPt->Draw("E");
	RatioPt->GetXaxis()->SetTitleOffset(1.2);
	RatioPt->GetYaxis()->SetTitleOffset(1.3);
	SetName(RatioPt, "HFT/TPC p_{T}");
	SetColor(RatioPt, kBlack, 20);

	TH1D *RatioPtN = (TH1D *)TrackHFTPtN->Clone("HFT/TPC p_{T}");
	RatioPtN->Divide(TrackPtN);
	RatioPtN->Draw("E SAME");
	SetName(RatioPtN, "HFT/TPC p_{T}");
	SetColor(RatioPtN, kRed, 24);
	legend->Draw("SAME");

	c[3]->SaveAs("QAPlots/HFTowerTPCPt_SL18fvSL22c.pdf");

	d[2]->cd();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TH1D *RatioNewOldPt = (TH1D *)RatioPtN->Clone("HFT/TPC p_{T} SL22c/SL18f");
	RatioNewOldPt->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
	RatioNewOldPt->GetYaxis()->SetTitle("#frac{SL22c}{SL18f}");
	RatioNewOldPt->Divide(RatioPt);
	RatioNewOldPt->Draw("E SAME");
	RatioNewOldPt->GetXaxis()->SetTitleOffset(1.2);
	RatioNewOldPt->GetYaxis()->SetTitleOffset(1.3);
	RatioNewOldPt->GetYaxis()->SetRangeUser(0.9, 1.1);
	RatioNewOldPt->GetYaxis()->SetNdivisions(505);
	SetName(RatioNewOldPt, "HFT/TPC p_{T} SL22c/SL18f");
	SetColor(RatioNewOldPt, kRed, 24);

	d[2]->SaveAs("QAPlots/HFTowerTPCPtDoubleRatio_SL18fvSL22c.pdf");

	c[4]->cd();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TH1D *RatioEta = (TH1D *)TrackHFTEta->Clone("HFT/TPC #eta");
	RatioEta->Divide(TrackEta);
	RatioEta->Draw("E");
	RatioEta->GetXaxis()->SetTitle("#eta");
	RatioEta->GetYaxis()->SetTitle("#frac{HFT}{TPC}");
	RatioEta->GetXaxis()->SetTitleOffset(1.2);
	RatioEta->GetYaxis()->SetTitleOffset(1.3);
	RatioEta->GetYaxis()->SetRangeUser(0, 0.95);
	RatioEta->GetYaxis()->SetNdivisions(505);
	SetName(RatioEta, "HFT/TPC #eta");
	SetColor(RatioEta, kBlack, 20);

	TH1D *RatioEtaN = (TH1D *)TrackHFTEtaN->Clone("HFT/TPC #eta");
	RatioEtaN->Divide(TrackEtaN);
	RatioEtaN->Draw("E SAME");
	SetName(RatioEtaN, "HFT/TPC #eta");
	RatioEtaN->GetXaxis()->SetTitleOffset(1.2);
	RatioEtaN->GetYaxis()->SetTitleOffset(1.3);
	SetColor(RatioEtaN, kRed, 24);
	legend->Draw("SAME");

	c[4]->SaveAs("QAPlots/HFTowerTPCeta_SL18fvSL22c.pdf");

	d[3]->cd();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TH1D *RatioNewOldEta = (TH1D *)RatioEtaN->Clone("HFT/TPC #eta SL22c/SL18f");
	RatioNewOldEta->GetXaxis()->SetTitle("#eta");
	RatioNewOldEta->GetYaxis()->SetTitle("#frac{SL22c}{SL18f}");
	RatioNewOldEta->Divide(RatioEta);
	RatioNewOldEta->Draw("E SAME");
	RatioNewOldEta->GetXaxis()->SetTitleOffset(1.2);
	RatioNewOldEta->GetYaxis()->SetTitleOffset(1.3);
	RatioNewOldEta->GetYaxis()->SetRangeUser(0.9, 1.1);
	RatioNewOldEta->GetYaxis()->SetNdivisions(505);
	SetName(RatioNewOldEta, "HFT/TPC #eta SL22c/SL18f");
	SetColor(RatioNewOldEta, kRed, 24);

	d[3]->SaveAs("QAPlots/HFTowerTPCetaDoubleRatio_SL18fvSL22c.pdf");


	c[5]->cd();
	TrackHFTPhi->Rebin(2);
	TrackPhi->Rebin(2);
	TH1D *RatioPhi = (TH1D *)TrackHFTPhi->Clone("HFT/TPC #phi");
	RatioPhi->Divide(TrackPhi);
	RatioPhi->Draw("E");
	SetName(RatioPhi, "HFT/TPC #phi");
	SetColor(RatioPhi, kBlack, 20);

	TrackHFTPhiN->Rebin(3);
	TrackPhiN->Rebin(3);
	TH1D *RatioPhiN = (TH1D *)TrackHFTPhiN->Clone("HFT/TPC #phi");
	RatioPhiN->Divide(TrackPhiN);
	RatioPhiN->Draw("E SAME");
	SetName(RatioPhiN, "HFT/TPC #phi");
	SetColor(RatioPhiN, kRed, 24);
	legend->Draw("SAME");

	d[4]->cd();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	PicoTowerVsProjectedTowerN->GetXaxis()->SetNdivisions(505);
	PicoTowerVsProjectedTowerN->GetYaxis()->SetNdivisions(505);
	SetAxisTitles(PicoTowerVsProjectedTowerN, "PicoDst Tower ID", "Projected Tower ID");
	SetName(PicoTowerVsProjectedTowerN, "Tower IDs Comparison");
	PicoTowerVsProjectedTowerN->GetXaxis()->SetTitleOffset(1.3);
	PicoTowerVsProjectedTowerN->GetYaxis()->SetTitleOffset(1.5);
	PicoTowerVsProjectedTowerN->Draw("COLZ");

	d[4]->SaveAs("QAPlots/TowerIDCalc_SL18fvSL22c.pdf");
	d[4]->SaveAs("QAPlots/TowerIDCalc_SL18fvSL22c.jpg");

	c[6]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Multiplicity->Scale(1./Centrality->Integral());
	MultiplicityN->Scale(1./CentralityN->Integral());
	Multiplicity->Draw("E");
	Multiplicity->GetYaxis()->SetRangeUser(9*pow(10, -8), 8);
	SetName(Multiplicity, "Multiplicity");
	SetColor(Multiplicity, kBlack, 20);
	MultiplicityN->Draw("E SAME");
	SetColor(MultiplicityN, kRed, 24);
	Multiplicity->GetXaxis()->SetTitle("<N_{Part}>");
	Multiplicity->GetXaxis()->SetNdivisions(505);
	legend->Draw("SAME");

	c[6]->SaveAs("QAPlots/Multiplicity_SL18fvSL22c.pdf");

	c[11]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	Centrality->Scale(1./Centrality->Integral());
	CentralityN->Scale(1./CentralityN->Integral());
	Centrality->Draw("E");
	CentralityN->GetYaxis()->SetRangeUser(9*pow(10, -8), 8);
	SetName(Centrality, "Centrality");
	SetColor(Centrality, kBlack, 20);
	CentralityN->Draw("E SAME");
	SetColor(CentralityN, kRed, 24);
	Centrality->GetXaxis()->SetTitle("Centrality (%)");
	Centrality->GetXaxis()->SetNdivisions(505);
	legend->Draw("SAME");

	c[11]->SaveAs("QAPlots/Centrality_SL18fvSL22c.pdf");

	c[7]->cd();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	VzVsVzVPD->Draw("COLZ");
	VzVsVzVPD->RebinX(2);
	VzVsVzVPD->RebinY(2);
	VzVsVzVPD->Scale(1./Centrality->Integral());
	VzVsVzVPD->GetZaxis()->SetRangeUser(pow(10, -8), pow(10, -1));
	VzVsVzVPD->GetXaxis()->SetTitle("V_{z} [cm]");
	VzVsVzVPD->GetYaxis()->SetTitle("V_{z}^{VPD} [cm]");
	VzVsVzVPD->GetXaxis()->SetTitleOffset(1.2);
	VzVsVzVPD->GetYaxis()->SetTitleOffset(1.3);
	VzVsVzVPD->SetTitle("V_{z} vs V_{z}^{VPD}");

	c[7]->SaveAs("QAPlots/VzVsVPD_SL18f.pdf");

	d[5]->cd();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	VzVsVzVPDN->Draw("COLZ");
	VzVsVzVPDN->RebinX(2);
	VzVsVzVPDN->RebinY(2);
	VzVsVzVPDN->Scale(1./CentralityN->Integral());
	VzVsVzVPDN->GetZaxis()->SetRangeUser(pow(10, -8), pow(10, -1));
	VzVsVzVPDN->GetXaxis()->SetTitle("V_{z} [cm]");
	VzVsVzVPDN->GetYaxis()->SetTitle("V_{z}^{VPD} [cm]");
	VzVsVzVPDN->GetXaxis()->SetTitleOffset(1.2);
	VzVsVzVPDN->GetYaxis()->SetTitleOffset(1.3);
	VzVsVzVPDN->SetTitle("V_{z} vs V_{z}^{VPD}");

	d[5]->SaveAs("QAPlots/VzVsVPD_SL22c.pdf");

	d[6]->cd();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TH2D *hVzVzVPDRatio = (TH2D *)VzVsVzVPDN->Clone("Ratio");
	hVzVzVPDRatio->Divide(VzVsVzVPD);
	hVzVzVPDRatio->Draw("COLZ");
	hVzVzVPDRatio->GetZaxis()->SetRangeUser(0.1, 2);

	c[8] ->cd();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	VxVsVy->Draw("COLZ");
	VxVsVy->Scale(1./Centrality->Integral());
	VxVsVy->GetZaxis()->SetRangeUser(pow(10, -8), pow(10, -1));
	VxVsVy->GetXaxis()->SetTitle("V_{x} [cm]");
	VxVsVy->GetYaxis()->SetTitle("V_{y} [cm]");
	VxVsVy->GetXaxis()->SetTitleOffset(1.2);
	VxVsVy->GetYaxis()->SetTitleOffset(1.3);
	VxVsVy->SetTitle("V_{x} vs V_{y}");

	c[8]->SaveAs("QAPlots/VxVy_SL18f.pdf");

	d[7]->cd();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	gPad->SetLogz();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	VxVsVyN->Draw("COLZ");
	VxVsVyN->Scale(1./CentralityN->Integral());
	VxVsVyN->GetZaxis()->SetRangeUser(pow(10, -8), pow(10, -1));
	VxVsVyN->GetXaxis()->SetTitle("V_{x} [cm]");
	VxVsVyN->GetYaxis()->SetTitle("V_{y} [cm]");
	VxVsVyN->GetXaxis()->SetTitleOffset(1.2);
	VxVsVyN->GetYaxis()->SetTitleOffset(1.3);
	VxVsVyN->SetTitle("V_{x} vs V_{y}");

	d[7]->SaveAs("QAPlots/VxVy_SL22c.pdf");


	c[9]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TowerE->Scale(1./Centrality->Integral());
	TowerEN->Scale(1./CentralityN->Integral());
	TowerE->Draw("E");
	TowerEN->Draw("E SAME");
	SetAxisTitles(TowerE, "Energy [GeV]", "Arb. Units");
	SetName(TowerE, "Tower Energy");
	SetColor(TowerE, kBlack, 20);
	SetColor(TowerEN, kRed, 24);
	TowerE->GetXaxis()->SetTitleOffset(1.2);
	TowerE->GetYaxis()->SetTitleOffset(1.4);
	legend->Draw("SAME");

	c[9]->SaveAs("QAPlots/TowerE_SL18fvSL22c.pdf");

	c[10]->cd();
	gPad->SetLogy();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TowerPhi->Scale(1./Centrality->Integral());
	TowerPhiN->Scale(1./CentralityN->Integral());
	TowerPhi->Draw("E");
	TowerPhiN->Draw("E SAME");
	TowerPhi->GetYaxis()->SetRangeUser(2*pow(10, -4), 999);
	SetAxisTitles(TowerPhi, "#phi", "Arb. Units");
	SetName(TowerPhi, "Tower ID vs Phi");
	SetColor(TowerPhi, kBlack, 20);
	SetColor(TowerPhiN, kRed, 24);
	TowerPhi->GetXaxis()->SetTitleOffset(1.2);
	TowerPhi->GetYaxis()->SetTitleOffset(1.4);
	legend->Draw("SAME");

	c[10]->SaveAs("QAPlots/TowerPhi_SL18fvSL22c.pdf");

	TString PicoPlotsName[6] = {"0.2 < p_{T} [GeV/#it{c}] < 0.5", 
								"0.5 < p_{T} [GeV/#it{c}] < 1.0",
								"1.0 < p_{T} [GeV/#it{c}] < 2.0",
								"2.0 < p_{T} [GeV/#it{c}] < 5.0",
								"5.0 < p_{T} [GeV/#it{c}] < 10.0",
								"10.0 < p_{T} [GeV/#it{c}] < 40.0"};
 
	for (int i = 8; i < 14; i++){
		d[i]->cd();
		gStyle->SetPalette(kInvertedDarkBodyRadiator);
		gPad->SetLogz();
		gPad->SetRightMargin(0.125);
		gPad->SetLeftMargin(0.125);
		gPad->SetTopMargin(0.1);
		gPad->SetBottomMargin(0.125);
		PicoTowerVsProjectedTowerPtN[i-8]->GetXaxis()->SetNdivisions(505);
		PicoTowerVsProjectedTowerPtN[i-8]->GetYaxis()->SetNdivisions(505);
		SetAxisTitles(PicoTowerVsProjectedTowerPtN[i-8], "PicoDst Tower ID", "Projected Tower ID");
		SetName(PicoTowerVsProjectedTowerPtN[i-8], PicoPlotsName[i-8].Data());
		PicoTowerVsProjectedTowerPtN[i-8]->GetXaxis()->SetTitleOffset(1.3);
		PicoTowerVsProjectedTowerPtN[i-8]->GetYaxis()->SetTitleOffset(1.5);
		PicoTowerVsProjectedTowerPtN[i-8]->Draw("COLZ");

		d[i]->SaveAs(Form("QAPlots/TowerIDCalc_%i_SL18fvSL22c.pdf", i-8));
		d[i]->SaveAs(Form("QAPlots/TowerIDCalc_%i_SL18fvSL22c.jpg", i-8));
	}
}

void Method2(int filenumber = 0){
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TFile *f = new TFile(Form("/Volumes/WorkDrive/AuAu200GeVQA/outsmallQA/SL18f_%i.root", filenumber));
	f->cd("TrackClusterQA");
	TTree *oldtree = (TTree *)gDirectory->Get("EventTree");

	TFile *g = new TFile(Form("/Volumes/WorkDrive/AuAu200GeVQA/outsmallQA/SL22c_%i.root", filenumber));
	g->cd("TrackClusterQA");
	TTree *newtree = (TTree *)gDirectory->Get("EventTree");

	int runidold;
	int eventidold;
	int hfttrackold;
	int tpctrackold;

	oldtree->SetBranchAddress("RunId", &runidold);
	oldtree->SetBranchAddress("EventId", &eventidold);
	oldtree->SetBranchAddress("HFTTracks", &hfttrackold);
	oldtree->SetBranchAddress("TPCTracks", &tpctrackold);

	int runidnew;
	int eventidnew;
	int hfttracknew;
	int tpctracknew;

	newtree->SetBranchAddress("RunId", &runidnew);
	newtree->SetBranchAddress("EventId", &eventidnew);
	newtree->SetBranchAddress("HFTTracks", &hfttracknew);
	newtree->SetBranchAddress("TPCTracks", &tpctracknew);

	if (oldtree->GetEntries() != newtree->GetEntries()){
		cout << "Check this file" << endl;
		return;
	}

	oldtree->GetEntry(0);

	int RUNNUM = runidold;

	double mineventid = eventidold;

	cout << mineventid << endl;

	oldtree->GetEntry(oldtree->GetEntries()-1);

	double maxeventid = eventidold;

	cout << maxeventid << endl;

	const int nbins = maxeventid - mineventid + 1;

	TH1D *HFTTPCTrackRatioPerEvent_Old = new TH1D("HFTTPCTrackRatioPerEvent_SL18f", "HFTTPCTrackRatioPerEvent_SL18f", nbins, mineventid - 0.5, maxeventid + 0.5);
	TH1D *HFTTPCTrackRatioPerEvent_New = new TH1D("HFTTPCTrackRatioPerEvent_SL22c", "HFTTPCTrackRatioPerEvent_SL22c", nbins, mineventid - 0.5, maxeventid + 0.5);

	for (int i = 0; i < oldtree->GetEntries(); i++){
		if (i%10!=0) continue;
		oldtree->GetEntry(i);
		double ratio = (double)hfttrackold/(double)tpctrackold;
		HFTTPCTrackRatioPerEvent_Old->Fill(eventidold, ratio);
	}

	for (int i = 0; i < newtree->GetEntries(); i++){
		if (i%10!=0) continue;
		newtree->GetEntry(i);
		double ratio = (double)hfttracknew/(double)tpctracknew;
		HFTTPCTrackRatioPerEvent_New->Fill(eventidnew, ratio);
	}

	TCanvas *c = new TCanvas("Ratio_Event_Event", "Ratio_Event_Event", 1000, 600);
	c->cd();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	SetName(HFTTPCTrackRatioPerEvent_Old, "SL18f");
	SetColor(HFTTPCTrackRatioPerEvent_Old, kBlack, 20);
	SetAxisTitles(HFTTPCTrackRatioPerEvent_Old, "Event ID", "#frac{HFT}{TPC}");

	SetName(HFTTPCTrackRatioPerEvent_New, "SL22c");
	SetColor(HFTTPCTrackRatioPerEvent_New, kRed, 24);

	HFTTPCTrackRatioPerEvent_Old->GetYaxis()->SetRangeUser(0.1, 1.5);
	HFTTPCTrackRatioPerEvent_Old->Draw("HIST P");
	HFTTPCTrackRatioPerEvent_New->Draw("HIST P SAME");

	auto legend = new TLegend(0.5,0.7,0.88,0.9);
	legend->SetHeader(Form("Run ID = %i", RUNNUM));
	legend->AddEntry(HFTTPCTrackRatioPerEvent_Old,"SL18f","p");
	legend->AddEntry(HFTTPCTrackRatioPerEvent_New,"SL22c","p");
	legend->Draw("SAME");

	c->SaveAs(Form("QAPlots/HFTvsTPC_PerEvent_SL18f_SL22c_File_%i.pdf", filenumber));

	TCanvas *d = new TCanvas("DoubleRatio_Event_Event", "DoubleRatio_Event_Event", 1000, 600);
	d->cd();
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.125);
	TH1D *Ratio = (TH1D *)HFTTPCTrackRatioPerEvent_New->Clone("SL22c/SL18f");
	Ratio->Divide(HFTTPCTrackRatioPerEvent_Old);
	Ratio->GetYaxis()->SetRangeUser(0.94, 1.06);
	SetAxisTitles(Ratio, "Event ID", "#frac{SL22c}{SL18f}");
	Ratio->Draw("HIST P");

	auto legend2 = new TLegend(0.5,0.7,0.88,0.9);
	legend2->SetHeader(Form("Run ID = %i", RUNNUM));
	legend2->Draw("SAME");

	d->SaveAs(Form("QAPlots/HFTvsTPCRatio_PerEvent_SL18f_SL22c_File_%i.pdf", filenumber));

	delete c;
	delete d;

	f->Close();
	g->Close();

}

void QAPlotter(){
	// Method1();

	for (int i = 0; i < 15; i++){
		Method2(i);
	}
}