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

void ppQA(){
	TFile *inputfile = new TFile("ppQA.root");
	TDirectory *QADir = (TDirectory *)inputfile->Get("PIDQA_TrackTree");
	TTree *TrackTree = (TTree *)QADir->Get("D0USTree");

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
  	cout << nentries << endl;

  	TH2F *hntofsigmapionvp = new TH2F("hntofsigmapionvp", "NTOFSigma Pion v p", 200, 0, 2, 200, -4, 4);
  	TH2F *hntofsigmakaonvp = new TH2F("hntofsigmakaonvp", "NTOFSigma Kaon v p", 200, 0, 2, 200, -4, 4);

  	TH2F *hntpcsigmapionvp = new TH2F("hntpcsigmapionvp", "NTPCSigma Pion v p", 200, 0, 2, 200, -4, 4);
  	TH2F *hntpcsigmakaonvp = new TH2F("hntpcsigmakaonvp", "NTPCSigma Kaon v p", 200, 0, 2, 200, -4, 4);

  	TH2F *hpionsigmapionvsigmakaon = new TH2F("hpionsigmapionvsigmakaon", "Pion Sigma Pion v Kaon", 300, -5, 5, 300, -5, 5);
  	TH2F *hkaonsigmakaonvsigmapion = new TH2F("hkaonsigmakaonvsigmapion", "Kaon Sigma Kaon v Pion", 300, -5, 5, 300, -5, 5);

  	TH1F *hD0Mass = new TH1F("hD0Mass", "hD0Mass", 100, 1.8, 1.9);

  	for (int event = 0; event < nentries;  event++){

  		if (event%10000==0)cout << event << " Events done." << endl;

   		TrackTree->GetEntry(event);

   		if (D0Pt < 1.0) continue;

   		double pionp = PionPt*TMath::CosH(PionEta);
   		double kaonp = KaonPt*TMath::CosH(KaonEta);

   		float oneOverBetaExpectedpion = sqrt(M_PION_PLUS*M_PION_PLUS / pionp / pionp + 1);
   		float oneOverBetaExpectedkaon = sqrt(M_KAON_PLUS*M_KAON_PLUS / kaonp / kaonp + 1);

   		if (PionTofBeta > 0) hntofsigmapionvp->Fill(pionp, (1./PionTofBeta - oneOverBetaExpectedpion)/0.012);
   		if (KaonTofBeta > 0) hntofsigmakaonvp->Fill(kaonp, (1./KaonTofBeta - oneOverBetaExpectedkaon)/0.012);

	   	if (abs(PionNSigmaPion) > 2) continue;
	   	if (abs(KaonNSigmaKaon) > 2) continue;
	   	// if (abs(PionNSigmaPion) > abs(PionNSigmaKaon)) continue;
	   	// if (abs(KaonNSigmaKaon) > abs(KaonNSigmaPion)) continue;
	   	// if (abs(KaonNSigmaPion) < 2) continue;
	   	// if (abs(PionNSigmaKaon) < 2) continue;
	   	if (abs(PionTofYLocal) > 1.8) continue;
	   	if (abs(KaonTofYLocal) > 1.8) continue;
	   	if (PionNHitsFit < 20) continue;
	   	if (KaonNHitsFit < 20) continue;

	   	hntpcsigmapionvp->Fill(pionp, PionNSigmaPion);
	   	hntpcsigmakaonvp->Fill(kaonp, KaonNSigmaKaon);

	   	hpionsigmapionvsigmakaon->Fill(PionNSigmaPion, PionNSigmaKaon);
	   	hkaonsigmakaonvsigmapion->Fill(KaonNSigmaPion, KaonNSigmaKaon);

	   	hD0Mass->Fill(D0Mass);
   	}

   	TCanvas *c = new TCanvas("c", "c", 1200, 500);
   	c->Divide(2);
   	c->cd(1);
   	gPad->SetLogz();
   	hntofsigmapionvp->Draw("COLZ");
   	c->cd(2);
   	gPad->SetLogz();
   	hntofsigmakaonvp->Draw("COLZ");

   	TCanvas *d = new TCanvas("d", "d", 1200, 500);
   	d->Divide(2);
   	d->cd(1);
   	gPad->SetLogz();
   	hntpcsigmapionvp->Draw("COLZ");
   	d->cd(2);
   	gPad->SetLogz();
   	hntpcsigmakaonvp->Draw("COLZ");

   	TCanvas *e = new TCanvas("e", "e", 1200, 500);
   	e->Divide(2);
   	e->cd(1);
   	gPad->SetLogz();
   	hpionsigmapionvsigmakaon->Draw("COLZ");
   	e->cd(2);
   	gPad->SetLogz();
   	hkaonsigmakaonvsigmapion->Draw("COLZ");

   	TCanvas *f = new TCanvas("f", "f", 1200, 500);
   	f->cd();
   	hD0Mass->Draw();
}

