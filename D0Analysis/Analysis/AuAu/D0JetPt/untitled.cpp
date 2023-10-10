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
#include "TRandom3.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
// #include "StJetTreeStruct.h"
#include <vector>

#pragma link C++ class vector<int> +;

using namespace std;

// #pragma link C++ class StJetTreeStruct+;

// #pragma link C++ class vector<float> +;
// #pragma link C++ class vector<vector<float> >+;
// #pragma link C++ class vector<int> +;
// #pragma link C++ class vector<vector<int> >+;
#endif

Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

Double_t standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}


void Writer(TString CentralityBin = "Central"){

	TStopwatch timer;

	timer.Start();

	gSystem->Load("./RooUnfold/build/libRooUnfold.dylib");

	// Declaring Global variables

	const int maximum = 5000;

	Int_t           mcRunID;
	Int_t           mcEventId;
	Float_t         mcRefMult;
	Float_t         mcCentrality;
	vector<int> *mcTriggers = new vector<int>;
	vector<double>  *mcPrimaryVertex = new vector<double>;
	vector<double>  *mcPrimaryVertexErr = new vector<double>;
	Float_t         mcJetPt;
	Float_t         mcJetCorrPt;
	Float_t         mcJetEta;
	Float_t         mcJetPhi;
	Float_t         mcJetArea;
	Float_t         mcJetRadius;
	Float_t         mcJetE;
	Float_t         mcJetNEF;
	Float_t         mcJetRhoVal;
	Float_t         mcJetHighestTrackPt;
	Int_t           mcJetNConst;
	Float_t         mcD0Mass;
	Float_t         mcD0Pt;
	Float_t         mcPionPt;
	Float_t         mcPionEta;
	Float_t         mcPionPhi;
	Float_t         mcPionCharge;
	Float_t         mcKaonPt;
	Float_t         mcKaonEta;
	Float_t         mcKaonPhi;
	Float_t         mcKaonCharge;
	Float_t         mcTrackID[maximum];   //[numberofconstituents]
	Float_t         mcTrackPt[maximum];   //[numberofconstituents]
	Float_t         mcTrackEta[maximum];   //[numberofconstituents]
	Float_t         mcTrackPhi[maximum];   //[numberofconstituents]
	Float_t         mcTrackPx[maximum];   //[numberofconstituents]
	Float_t         mcTrackPy[maximum];   //[numberofconstituents]
	Float_t         mcTrackPz[maximum];   //[numberofconstituents]
	Float_t         mcTrackCharge[maximum];   //[numberofconstituents]

	Int_t           recoRunID;
	Int_t           recoEventId;
	Float_t         recoRefMult;
	Float_t         recoCentrality;
	vector<int> *recoTriggers = new vector<int>;
	vector<double>  *recoPrimaryVertex = new vector<double>;
	vector<double>  *recoPrimaryVertexErr = new vector<double>;
	Float_t         recoJetPt;
	Float_t         recoJetCorrPt;
	Float_t         recoJetEta;
	Float_t         recoJetPhi;
	Float_t         recoJetArea;
	Float_t         recoJetRadius;
	Float_t         recoJetE;
	Float_t         recoJetNEF;
	Float_t         recoJetRhoVal;
	Float_t         recoJetHighestTrackPt;
	Int_t           recoJetNConst;
	Float_t         recoD0Mass;
	Float_t         recoD0Pt;
	Float_t         recoPionPt;
	Float_t         recoPionEta;
	Float_t         recoPionPhi;
	Float_t         recoPionCharge;
	Float_t         recoKaonPt;
	Float_t         recoKaonEta;
	Float_t         recoKaonPhi;
	Float_t         recoKaonCharge;
	Float_t         recoTrackID[maximum];   //[numberofconstituents]
	Float_t         recoTrackPt[maximum];   //[numberofconstituents]
	Float_t         recoTrackEta[maximum];   //[numberofconstituents]
	Float_t         recoTrackPhi[maximum];   //[numberofconstituents]
	Float_t         recoTrackPx[maximum];   //[numberofconstituents]
	Float_t         recoTrackPy[maximum];   //[numberofconstituents]
	Float_t         recoTrackPz[maximum];   //[numberofconstituents]
	Float_t         recoTrackCharge[maximum];   //[numberofconstituents]


	TFile *filein;

	filein = TFile::Open(Form("../ResponseMatrix/%s/Response_%s.root", CentralityBin.Data(), CentralityBin.Data()));

	TString mcinputstring;
	TString recoinputstring;

	mcinputstring = "SimJetSaverSmear/MCJets";
	recoinputstring = "SimJetSaverSmear/RecoJets";

	// TFile *f = TFile::Open("test.root");

	TTree *mcjettree = (TTree*)filein->Get(mcinputstring.Data());
	TTree *recojettree = (TTree*)filein->Get(recoinputstring.Data());

	mcjettree->SetBranchAddress("RunID", &mcRunID);
	mcjettree->SetBranchAddress("EventId", &mcEventId);
	mcjettree->SetBranchAddress("RefMult", &mcRefMult);
	mcjettree->SetBranchAddress("Centrality", &mcCentrality);
	mcjettree->SetBranchAddress("Triggers", &mcTriggers);
	mcjettree->SetBranchAddress("PrimaryVertex", &mcPrimaryVertex);
	mcjettree->SetBranchAddress("PrimaryVertexErr", &mcPrimaryVertexErr);
	mcjettree->SetBranchAddress("JetPt", &mcJetPt);
	mcjettree->SetBranchAddress("JetCorrPt", &mcJetCorrPt);
	mcjettree->SetBranchAddress("JetEta", &mcJetEta);
	mcjettree->SetBranchAddress("JetPhi", &mcJetPhi);
	mcjettree->SetBranchAddress("JetArea", &mcJetArea);
	mcjettree->SetBranchAddress("JetRadius", &mcJetRadius);
	mcjettree->SetBranchAddress("JetE", &mcJetE);
	mcjettree->SetBranchAddress("JetNEF", &mcJetNEF);
	mcjettree->SetBranchAddress("JetRhoVal", &mcJetRhoVal);
	mcjettree->SetBranchAddress("JetHighestTrackPt", &mcJetHighestTrackPt);
	mcjettree->SetBranchAddress("JetNConst", &mcJetNConst);
	mcjettree->SetBranchAddress("D0Mass", &mcD0Mass);
	mcjettree->SetBranchAddress("D0Pt", &mcD0Pt);
	mcjettree->SetBranchAddress("PionPt", &mcPionPt);
	mcjettree->SetBranchAddress("PionEta", &mcPionEta);
	mcjettree->SetBranchAddress("PionPhi", &mcPionPhi);
	mcjettree->SetBranchAddress("PionCharge", &mcPionCharge);
	mcjettree->SetBranchAddress("KaonPt", &mcKaonPt);
	mcjettree->SetBranchAddress("KaonEta", &mcKaonEta);
	mcjettree->SetBranchAddress("KaonPhi", &mcKaonPhi);
	mcjettree->SetBranchAddress("KaonCharge", &mcKaonCharge);
	mcjettree->SetBranchAddress("TrackID", mcTrackID);
	mcjettree->SetBranchAddress("TrackPt", mcTrackPt);
	mcjettree->SetBranchAddress("TrackEta", mcTrackEta);
	mcjettree->SetBranchAddress("TrackPhi", mcTrackPhi);
	mcjettree->SetBranchAddress("TrackPx", mcTrackPx);
	mcjettree->SetBranchAddress("TrackPy", mcTrackPy);
	mcjettree->SetBranchAddress("TrackPz", mcTrackPz);
	mcjettree->SetBranchAddress("TrackCharge", mcTrackCharge);

	recojettree->SetBranchAddress("RunID", &recoRunID);
	recojettree->SetBranchAddress("EventId", &recoEventId);
	recojettree->SetBranchAddress("RefMult", &recoRefMult);
	recojettree->SetBranchAddress("Centrality", &recoCentrality);
	recojettree->SetBranchAddress("Triggers", &recoTriggers);
	recojettree->SetBranchAddress("PrimaryVertex", &recoPrimaryVertex);
	recojettree->SetBranchAddress("PrimaryVertexErr", &recoPrimaryVertexErr);
	recojettree->SetBranchAddress("JetPt", &recoJetPt);
	recojettree->SetBranchAddress("JetCorrPt", &recoJetCorrPt);
	recojettree->SetBranchAddress("JetEta", &recoJetEta);
	recojettree->SetBranchAddress("JetPhi", &recoJetPhi);
	recojettree->SetBranchAddress("JetArea", &recoJetArea);
	recojettree->SetBranchAddress("JetRadius", &recoJetRadius);
	recojettree->SetBranchAddress("JetE", &recoJetE);
	recojettree->SetBranchAddress("JetNEF", &recoJetNEF);
	recojettree->SetBranchAddress("JetRhoVal", &recoJetRhoVal);
	recojettree->SetBranchAddress("JetHighestTrackPt", &recoJetHighestTrackPt);
	recojettree->SetBranchAddress("JetNConst", &recoJetNConst);
	recojettree->SetBranchAddress("D0Mass", &recoD0Mass);
	recojettree->SetBranchAddress("D0Pt", &recoD0Pt);
	recojettree->SetBranchAddress("PionPt", &recoPionPt);
	recojettree->SetBranchAddress("PionEta", &recoPionEta);
	recojettree->SetBranchAddress("PionPhi", &recoPionPhi);
	recojettree->SetBranchAddress("PionCharge", &recoPionCharge);
	recojettree->SetBranchAddress("KaonPt", &recoKaonPt);
	recojettree->SetBranchAddress("KaonEta", &recoKaonEta);
	recojettree->SetBranchAddress("KaonPhi", &recoKaonPhi);
	recojettree->SetBranchAddress("KaonCharge", &recoKaonCharge);
	recojettree->SetBranchAddress("TrackID", recoTrackID);
	recojettree->SetBranchAddress("TrackPt", recoTrackPt);
	recojettree->SetBranchAddress("TrackEta", recoTrackEta);
	recojettree->SetBranchAddress("TrackPhi", recoTrackPhi);
	recojettree->SetBranchAddress("TrackPx", recoTrackPx);
	recojettree->SetBranchAddress("TrackPy", recoTrackPy);
	recojettree->SetBranchAddress("TrackPz", recoTrackPz);
	recojettree->SetBranchAddress("TrackCharge", recoTrackCharge);

	TFile *weightfile = new TFile("./WeightFile.root");
	weightfile->cd();

	TH2D *D0PtvZ_Peripheral = (TH2D *)weightfile->Get("RatioDataSIMPeripheral");

	TFile *fonllweight = new TFile("./FONLLWeights.root");
	fonllweight->cd();

	TH1D *FONLL = (TH1D *)fonllweight->Get("FONLLWeights");

	const Int_t nBinsdR = 5;
	Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4};

	const int nBinsJetPt = 12;
  double JetPtBins[nBinsJetPt+1] = {0.,1, 2, 3.,4.,5.,7.,9.,11.,13.,15.,20.,30};


  const int nBinsJetPtMC = 12;
  double JetPtBinsMC[nBinsJetPtMC+1] = {0,1,2,3.,4.,5.,7.,9.,11.,13.,15.,20.,30};


	const int nbinsZ = 16;
	double ZBins[nbinsZ+1] = {-2.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

	const int nBinsD0Spectra = 7;
  double D0PtBins[nBinsD0Spectra + 1] = {1., 1.5, 2., 2.5, 3., 4., 5., 10};

	// Dimension defining histogram only

	TH1D *hRecoD0Pt = new TH1D("hRecoD0Pt", "hRecoD0Pt", nBinsD0Spectra, D0PtBins);
	TH1D *hRecoJetPt = new TH1D("hRecoJetPt", "hRecoJetPt", 100, -50, 50);
	TH1D *hRecoJetCorrPt = new TH1D("hRecoJetCorrPt", "hRecoJetCorrPt", 100, -50, 50);

	TH1D *hMCdRBeforeWeight = new TH1D("hMCdRBeforeWeight", "hMCdRBeforeWeight", nBinsdR, dRBins);
	TH1D *hRecodRBeforeWeight = new TH1D("hRecodRBeforeWeight", "hRecodRBeforeWeight", nBinsdR, dRBins);

	TH1D *hMCdRAfterWeight = new TH1D("hMCdRAfterWeight", "hMCdRAfterWeight", nBinsdR, dRBins);
	TH1D *hRecodRAfterWeight = new TH1D("hRecodRAfterWeight", "hRecodRAfterWeight", nBinsdR, dRBins);

	TH1D *hist_measured = new TH1D("hist_measured", "hist_measured", nBinsJetPt, JetPtBins);
	TH1D *hist_truth = new TH1D("hist_truth", "hist_truth", nBinsJetPtMC, JetPtBinsMC);
	TH2F *signalJetPtvdRD0Jet = new TH2F("signalJetPtvdRD0Jet", "signalJetPtvdRD0Jet", nBinsJetPt, JetPtBins, nBinsdR, dRBins);

	// The histograms above this line are not written and not filled.

	TH1D *mchist_z = new TH1D("mchist_z", "mchist_z", nbinsZ, ZBins);
	TH1D *recohist_z = new TH1D("recohist_z", "recohist_z", nbinsZ, ZBins);

	TH2D *mcJetPtvZ = new TH2D("mcJetPtvZ", "mcJetPtvZ", nBinsJetPt, JetPtBins, nbinsZ, ZBins);
	TH2D *recoJetPtvZ = new TH2D("recoJetPtvZ", "recoJetPtvZ", nBinsJetPt, JetPtBins, nbinsZ, ZBins);

	TH1D *hTrueMCPt = new TH1D("hTrueMCPt", "hTrueMCPt", nBinsJetPtMC, JetPtBinsMC);
	TH1D *hTrueRecoPt = new TH1D("hTrueRecoPt", "hTrueRecoPt", nBinsJetPt, JetPtBins);
	TH2D *hTrueRecoMCPt = new TH2D("hTrueRecoMCPt", "hTrueRecoMCPt", nBinsJetPt, JetPtBins, nBinsJetPtMC, JetPtBinsMC);

	TH1D *hTrueMCdR = new TH1D("hTrueMCdR", "hTrueMCdR", nBinsdR, dRBins);
	TH1D *hTrueRecodR = new TH1D("hTrueRecodR", "hTrueRecodR", nBinsdR, dRBins);
	TH2D *hTrueRecoMCdR = new TH2D("hTrueRecoMCdR", "hTrueRecoMCdR", nBinsdR, dRBins, nBinsdR, dRBins);

	TH2D *hMCPtvDeltaR = new TH2D("hMCPtvDeltaR", "hMCPtvDeltaR", nBinsJetPt, JetPtBins, nBinsdR, dRBins);
	TH2D *hRecoPtvDeltaR = new TH2D("hRecoPtvDeltaR", "hRecoPtvDeltaR", nBinsJetPt, JetPtBins, nBinsdR, dRBins);


	int bins[4] = {nBinsJetPt, nBinsdR, nBinsJetPt, nBinsdR};
	THnD *hRecoMCPtvDeltaR = new THnD("hRecoMCPtvDeltaR", "hRecoMCPtvDeltaR", 4, bins, NULL, NULL);
	hRecoMCPtvDeltaR->SetBinEdges(0, JetPtBins);
	hRecoMCPtvDeltaR->SetBinEdges(1, dRBins);
	hRecoMCPtvDeltaR->SetBinEdges(2, JetPtBins);
	hRecoMCPtvDeltaR->SetBinEdges(3, dRBins);

	RooUnfoldResponse* response, *response2D;

	response = new RooUnfoldResponse(hist_measured, hist_truth, "Response", "Response");
	response2D = new RooUnfoldResponse(signalJetPtvdRD0Jet, signalJetPtvdRD0Jet, "Response2D", "Response2D");

	int mcentries = mcjettree->GetEntries();
	int recoentries = recojettree->GetEntries();

	int numberoffakes = 0;
	int numberofmisses = 0;
	int numberofoks = 0;
	
	// for (int event = 0; event < 1000;  event++){
	for (int event = 0; event < mcentries;  event++){
		if (event%10000 == 0) cout << event << " Events analysed." << endl;
		mcjettree->GetEntry(event);
		recojettree->GetEntry(event);

		// Global EVENT CUTS (These events get tossed.)
		if (mcD0Pt <= 1 || mcD0Pt >= 10) continue;
		if (mcJetNConst == 0) continue;
		if (recoJetNConst == 0) continue; 


		// Smear the Detector Level Jet Pt with the background fluctuations. These have been obtained on March 14, 2022, so are the most recent values we have.

		TRandom *r = new TRandom(0);
		if (CentralityBin.CompareTo("Central")==0) recoJetCorrPt = recoJetPt + r->Gaus(0., 5.9);
		else if (CentralityBin.CompareTo("MidCentral")==0) recoJetCorrPt = recoJetPt + r->Gaus(0., 4.6);
		else if (CentralityBin.CompareTo("Peripheral")==0) recoJetCorrPt = recoJetPt + r->Gaus(0., 2.3);

		hRecoJetPt->Fill(recoJetPt);
		hRecoJetCorrPt->Fill(recoJetCorrPt); // We also need to fill a 2D histogram here.

		bool fake = kFALSE;
		bool miss = kFALSE;
		bool ok = kFALSE;

		// Jet Specific cuts defined from this point

		double deltaR = TMath::Sqrt(pow(mcJetEta - recoJetEta, 2) + pow(dPhi(mcJetPhi, recoJetPhi), 2)); // This is the jet matching condition

		double mcdeltaeta, mcdeltaphi, mcdeltar; // This is the D0 Delta R in the Particle Level Jet
		for(int itrk = 0; itrk < mcJetNConst; itrk++) {
			if (mcTrackID[itrk] == 421){
				mcdeltaeta = abs(mcJetEta - mcTrackEta[itrk]);
				mcdeltaphi = dPhi(mcJetPhi, standardPhi(mcTrackPhi[itrk]));
				mcdeltar = TMath::Sqrt(pow(mcdeltaeta, 2) + pow(mcdeltaphi, 2));
			}
		}

		double recodeltaeta, recodeltaphi, recodeltar; // This is the D0 Delta R in the Detector Level Jet
		for(int itrk = 0; itrk < recoJetNConst; itrk++) {
			if (recoTrackID[itrk] == 421){
				recodeltaeta = abs(recoJetEta - recoTrackEta[itrk]);
				recodeltaphi = dPhi(recoJetPhi, standardPhi(recoTrackPhi[itrk]));
				recodeltar = TMath::Sqrt(pow(recodeltaeta, 2) + pow(recodeltaphi, 2));
			}
		}

		// Define OK

		if (mcJetPt > 3. && mcJetPt < 30. && mcdeltar > 0. && mcdeltar < 0.4 && recoJetCorrPt > 3. && recoJetCorrPt < 30. && recodeltar > 0. && recodeltar < 0.4 && recoD0Pt > 1.0 && recoD0Pt < 10.0 && deltaR > 0. && deltaR < 0.4) ok = kTRUE;

		// Define FAKE

		if ((mcJetPt <= 3. || mcJetPt >= 30. || mcdeltar <= 0. || mcdeltar >= 0.4 ) && (recoJetCorrPt > 3. && recoJetCorrPt < 30. && recodeltar > 0. && recodeltar < 0.4 && recoD0Pt > 1.0 && recoD0Pt < 10.0)) fake = kTRUE;

		// Define MISS

		if ((mcJetPt > 3. && mcJetPt < 30. && mcdeltar > 0. && mcdeltar < 0.4) && (recoJetCorrPt <= 3. || recoJetCorrPt <= 30. || recodeltar <= 0. || recodeltar >= 0.4 || recoD0Pt <= 1.0 || recoD0Pt >= 10.0)) miss = kTRUE;

		// Check to make sure there's no event which is both fake/miss, fake/ok, miss/ok

		if ((fake && ok) || (fake && miss) || (miss && ok)) { cout << "Something's wrong! Terminating!" << endl; return; } // Should never be true!

		if (!fake && !miss && !ok) continue;

		if (ok){
			numberofoks++;

			double mcdeltar, mcz, mcpt;
			for(int itrk = 0; itrk < mcJetNConst; itrk++) {
				if (mcTrackID[itrk] == 421){
					mcdeltar = TMath::Sqrt(pow(mcJetEta - mcTrackEta[itrk], 2) + pow(dPhi(mcJetPhi, standardPhi(mcTrackPhi[itrk])), 2));
					double jetpx = mcJetCorrPt*TMath::Cos(mcJetPhi);
	        double jetpy = mcJetCorrPt*TMath::Sin(mcJetPhi);
	        mcz = (jetpx*mcTrackPx[itrk] + jetpy*mcTrackPy[itrk])/pow(mcJetCorrPt, 2);
	        // mcpt = mcTrackPt[itrk];
				}
			}

			double recopt, w;

			double recodeltar, recoz;
			for(int itrk = 0; itrk < recoJetNConst; itrk++) {
				if (recoTrackID[itrk] == 421){
					recodeltar = TMath::Sqrt(pow(recoJetEta - recoTrackEta[itrk], 2) + pow(dPhi(recoJetPhi, standardPhi(recoTrackPhi[itrk])), 2));
					double jetpx = recoJetCorrPt*TMath::Cos(recoJetPhi);
	        double jetpy = recoJetCorrPt*TMath::Sin(recoJetPhi);
	        recoz = (jetpx*recoTrackPx[itrk] + jetpy*recoTrackPy[itrk])/pow(recoJetCorrPt, 2);
	        recopt = recoTrackPt[itrk];
				}
			}

			int binx = D0PtvZ_Peripheral->GetXaxis()->FindBin(recopt);
			int biny = D0PtvZ_Peripheral->GetYaxis()->FindBin(recoz);

			// w = D0PtvZ_Peripheral->GetBinContent(binx, biny);

			int mcbinx = FONLL->GetXaxis()->FindBin(mcJetPt);
			w = FONLL->GetBinContent(mcbinx);

			hRecoD0Pt->Fill(recoD0Pt, w);

			response->Fill(recoJetCorrPt, mcJetPt, w);
			hTrueMCPt->Fill(mcJetPt, w);
	  	hTrueRecoPt->Fill(recoJetCorrPt, w);
	  	hTrueRecoMCPt->Fill(recoJetCorrPt, mcJetPt, w);

	  	response2D->Fill(recoJetCorrPt, recodeltar, mcJetPt, mcdeltar, w);
	  	hMCPtvDeltaR->Fill(mcJetPt, mcdeltar, w);
	  	hRecoPtvDeltaR->Fill(recoJetCorrPt, recodeltar, w);

	  	double delvec[4] = {recoJetCorrPt, recodeltar, mcJetPt, mcdeltar};
	  	hRecoMCPtvDeltaR->Fill(delvec, w);

	  	mchist_z->Fill(mcz, w);
	  	recohist_z->Fill(recoz, w);
	  	mcJetPtvZ->Fill(mcJetPt, mcz, w);
	  	recoJetPtvZ->Fill(recoJetCorrPt, recoz, w);

	  	hMCdRBeforeWeight->Fill(mcdeltar);
	  	hMCdRAfterWeight->Fill(mcdeltar, w);

	  	hRecodRBeforeWeight->Fill(recodeltar);
	  	hRecodRAfterWeight->Fill(recodeltar, w);

	  	hTrueMCdR->Fill(mcdeltar, w);
	  	hTrueRecodR->Fill(recodeltar, w);
	  	hTrueRecoMCdR->Fill(recodeltar, mcdeltar, w);
	  	
		}

		if (fake){

			numberoffakes++;

			double recopt, w;

			double recodeltar, recoz;
			for(int itrk = 0; itrk < recoJetNConst; itrk++) {
				if (recoTrackID[itrk] == 421){
					recodeltar = TMath::Sqrt(pow(recoJetEta - recoTrackEta[itrk], 2) + pow(standardPhi(recoJetPhi) - standardPhi(recoTrackPhi[itrk]), 2));
					double jetpx = recoJetCorrPt*TMath::Cos(recoJetPhi);
	        double jetpy = recoJetCorrPt*TMath::Sin(recoJetPhi);
	        recoz = (jetpx*recoTrackPx[itrk] + jetpy*recoTrackPy[itrk])/pow(recoJetCorrPt, 2);
	        recopt = recoTrackPt[itrk];
				}
			}

			int binx = D0PtvZ_Peripheral->GetXaxis()->FindBin(recopt);
			int biny = D0PtvZ_Peripheral->GetYaxis()->FindBin(recoz);

			// w = D0PtvZ_Peripheral->GetBinContent(binx, biny);
			int mcbinx = FONLL->GetXaxis()->FindBin(mcJetPt);
			w = FONLL->GetBinContent(mcbinx);

			hRecoD0Pt->Fill(recoD0Pt, w);
			response->Fake(recoJetCorrPt, w);
			hTrueRecoPt->Fill(recoJetCorrPt, w);

			response2D->Fake(recoJetCorrPt, recodeltar, w);
			hRecoPtvDeltaR->Fill(recoJetCorrPt, recodeltar, w);

			recohist_z->Fill(recoz, w);
			recoJetPtvZ->Fill(recoJetCorrPt, recoz, w);

			hRecodRBeforeWeight->Fill(recodeltar);
	  	hRecodRAfterWeight->Fill(recodeltar, w);

	  	hTrueRecodR->Fill(recodeltar, w);
		}

		if (miss){

			numberofmisses++;

			double mcdeltar, mcz, mcpt;

			double w;

			for(int itrk = 0; itrk < mcJetNConst; itrk++) {
				if (mcTrackID[itrk] == 421){
					mcdeltar = TMath::Sqrt(pow(mcJetEta - mcTrackEta[itrk], 2) + pow(standardPhi(mcJetPhi) - standardPhi(mcTrackPhi[itrk]), 2));
					double jetpx = mcJetCorrPt*TMath::Cos(mcJetPhi);
	        double jetpy = mcJetCorrPt*TMath::Sin(mcJetPhi);
	        mcz = (jetpx*mcTrackPx[itrk] + jetpy*mcTrackPy[itrk])/pow(mcJetCorrPt, 2);
	        // mcpt = mcTrackPt[itrk];
				}
			}

			int mcbinx = FONLL->GetXaxis()->FindBin(mcJetPt);
			w = FONLL->GetBinContent(mcbinx);

			hRecoD0Pt->Fill(recoD0Pt, w);
			response->Miss(mcJetPt, w);
			hTrueMCPt->Fill(mcJetPt, w);

			response2D->Miss(mcJetPt, mcdeltar, w);
	  	hMCPtvDeltaR->Fill(mcJetPt, mcdeltar, w);

	  	mchist_z->Fill(mcz, w);
	  	mcJetPtvZ->Fill(mcJetPt, mcz, w);

	  	hMCdRBeforeWeight->Fill(mcdeltar);
	  	hMCdRAfterWeight->Fill(mcdeltar, w);

	  	hTrueMCdR->Fill(mcdeltar, w);
		}

	}

	cout << "OK:MISS:FAKE = " << numberofoks << "\t" << numberofmisses << "\t" << numberoffakes << endl;

	TFile *f = new TFile(Form("ResponseMatrix_Weighted_%s.root", CentralityBin.Data()), "RECREATE");
	f->cd();
	f->WriteObject(response, "Response");
	f->WriteObject(response2D, "Response2D");

	hTrueMCPt->Write();
	hTrueRecoPt->Write();
	hTrueRecoMCPt->Write();

	hMCPtvDeltaR->Write();
	hRecoPtvDeltaR->Write();
	hRecoMCPtvDeltaR->Write();

	mchist_z->Write();
	recohist_z->Write();
	mcJetPtvZ->Write();
	recoJetPtvZ->Write();

	hRecoD0Pt->Write();

	hRecoJetPt->Write();
	hRecoJetCorrPt->Write();

	hMCdRBeforeWeight->Write();
	hMCdRAfterWeight->Write();

	hRecodRBeforeWeight->Write();
	hRecodRAfterWeight->Write();

	hTrueMCdR->Write();
	hTrueRecodR->Write();
	hTrueRecoMCdR->Write();

	f->Close();

	mcjettree->Reset();
	recojettree->Reset();
	filein->Close();

}

void ResponseMatricesForTesting(){
	Writer("Central");
	Writer("MidCentral");
	Writer("Peripheral");
}


