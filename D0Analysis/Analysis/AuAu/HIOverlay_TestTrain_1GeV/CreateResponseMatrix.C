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

#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"

#pragma link C++ class vector<int> +;

using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

void Method(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf"){

	TString inputFileName = Form("%s/MCMCResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins);
	TFile *f = new TFile(inputFileName.Data(), "UPDATE");
	f->cd();

	THnSparseF *hJet[3];
	THnSparseF *hJetZNorm[3];
	THnSparseF *hJetZNormPtNorm[3];

	for (int cent = 0; cent < 3; cent++){
		hJet[cent] = (THnSparseF *)f->Get(Form("hJet_%i", cent));
		hJetZNorm[cent] = (THnSparseF *)f->Get(Form("hJetZNorm_%i", cent));
		hJetZNormPtNorm[cent] = (THnSparseF *)f->Get(Form("hJetZNormPtNorm_%i", cent));

		cout << hJet[cent]->GetEntries() << "\t" << hJetZNorm[cent]->GetEntries() << "\t" << hJetZNormPtNorm[cent]->GetEntries() << endl;
	}

	const Int_t ndim = 4;
  	Int_t dim[ndim];
	for(Int_t i = 0; i<ndim; i++) dim[i] = i;

	Int_t iPtTrue = 0;
	Int_t iZTrue = 1;
	Int_t iPtDet = 2;
	Int_t iZDet = 3;

	TH2D *fReco[3]; //Measured
	TH2D *fMC[3]; //Truth

	RooUnfoldResponse *resp[3]; //Response

	TH2D *fMeas[3];
	TH2D *fTrue[3];
	TH2D *fMiss[3];

	for (int cent = 0; cent < 3; cent++){
		fReco[cent] = (TH2D *)hJetZNormPtNorm[cent]->Projection(iZDet, iPtDet, "E");
		fMC[cent] = (TH2D *)hJetZNormPtNorm[cent]->Projection(iZTrue, iPtTrue, "E");

		fMeas[cent] = new TH2D(Form("fMeas_Cent_%i", cent), Form("fMeas_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);
		fTrue[cent] = new TH2D(Form("fTrue_Cent_%i", cent), Form("fTrue_Cent_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
		fMiss[cent] = new TH2D(Form("fMiss_Cent_%i", cent), Form("fMiss_Cent_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

		resp[cent] = new RooUnfoldResponse(Form("Resp_%i", cent), Form("Resp_%i", cent));
		resp[cent]->Setup(fMeas[cent], fTrue[cent]); //Setup Response Matrix Definition
	}

	//Fill Detector Info
	for (int cent = 0; cent < 3; cent++){
		for ( int ibinx = 1; ibinx <= fMeas[cent]->GetNbinsX(); ibinx++ ){
			Double_t xlow  = fMeas[cent]->GetXaxis()->GetBinLowEdge(ibinx);
			Double_t xhigh = fMeas[cent]->GetXaxis()->GetBinUpEdge(ibinx);
			Int_t ptbinlow = fReco[cent]->GetXaxis()->FindBin(xlow+0.000001);
			Int_t ptbinhigh = fReco[cent]->GetXaxis()->FindBin(xhigh-0.000001);

			for ( int ibiny = 1; ibiny <= fMeas[cent]->GetNbinsY(); ibiny++ ){
				Double_t ylow  = fMeas[cent]->GetYaxis()->GetBinLowEdge(ibiny);
				Double_t yhigh = fMeas[cent]->GetYaxis()->GetBinUpEdge(ibiny);
				Int_t zbinlow = fReco[cent]->GetYaxis()->FindBin(ylow+0.000001);
				Int_t zbinhigh = fReco[cent]->GetYaxis()->FindBin(yhigh-0.000001);

				Double_t err = 0;
				Double_t content = fReco[cent]->IntegralAndError(ptbinlow, ptbinhigh, zbinlow, zbinhigh, err);
				// cout << "Pt vs Z = " << ibinx << "\t" << ibiny << "\t" << ptbinlow << "\t" << ptbinhigh << "\t" << zbinlow << "\t" << zbinhigh << "\t" << content << endl;

				fMeas[cent]->SetBinContent(ibinx, ibiny, content);
				fMeas[cent]->SetBinError(ibinx, ibiny, err);
			}
		}
	}

	//Fill Particle Info
	for (int cent = 0; cent < 3; cent++){
		for ( int ibinx = 1; ibinx <= fTrue[cent]->GetNbinsX(); ibinx++ ){
			Double_t xlow  = fTrue[cent]->GetXaxis()->GetBinLowEdge(ibinx);
			Double_t xhigh = fTrue[cent]->GetXaxis()->GetBinUpEdge(ibinx);
			Int_t ptbinlow = fMC[cent]->GetXaxis()->FindBin(xlow+0.000001);
			Int_t ptbinhigh = fMC[cent]->GetXaxis()->FindBin(xhigh-0.000001);

			for ( int ibiny = 1; ibiny <= fTrue[cent]->GetNbinsY(); ibiny++ ){
				Double_t ylow  = fTrue[cent]->GetYaxis()->GetBinLowEdge(ibiny);
				Double_t yhigh = fTrue[cent]->GetYaxis()->GetBinUpEdge(ibiny);
				Int_t zbinlow = fMC[cent]->GetYaxis()->FindBin(ylow+0.000001);
				Int_t zbinhigh = fMC[cent]->GetYaxis()->FindBin(yhigh-0.000001);

				Double_t err = 0;
				Double_t content = fMC[cent]->IntegralAndError(ptbinlow, ptbinhigh, zbinlow, zbinhigh, err);
				// cout << "Pt vs Z = " << ibinx << "\t" << ibiny << "\t" << ptbinlow << "\t" << ptbinhigh << "\t" << zbinlow << "\t" << zbinhigh << "\t" << content << endl;

				fTrue[cent]->SetBinContent(ibinx, ibiny, content);
				fTrue[cent]->SetBinError(ibinx, ibiny, err);
			}
		}
	}

	//Fill RooUnfoldResponse Object

	int* coord = new int[ndim];

	for (int cent = 0; cent < 3; cent++){
		int nbin = hJetZNormPtNorm[cent]->GetNbins();

		cout << "Bins = " << nbin << endl;

		for (int bin = 0; bin < nbin; bin++){
			Double_t w = hJetZNormPtNorm[cent]->GetBinContent(bin, coord);
			Double_t pttrue = hJetZNormPtNorm[cent]->GetAxis(0)->GetBinCenter(coord[0]);
			Double_t ztrue = hJetZNormPtNorm[cent]->GetAxis(1)->GetBinCenter(coord[1]);
			Double_t ptdet =  hJetZNormPtNorm[cent]->GetAxis(2)->GetBinCenter(coord[2]); 
			Double_t zdet =  hJetZNormPtNorm[cent]->GetAxis(3)->GetBinCenter(coord[3]);

			if (zdet >= z_low && zdet <= z_high
			&&	ztrue >= z_gen_low && ztrue <= z_gen_high
			&& ptdet >= jetpt_low && ptdet <= jetpt_high
			&& pttrue >= jetpt_gen_low && pttrue <= jetpt_gen_high
			){
				resp[cent]->Fill(ptdet, zdet, pttrue, ztrue, w);
			}
			else{
				resp[cent]->Miss(pttrue, ztrue, w);
				fMiss[cent]->Fill(pttrue, ztrue, w);
			} 
		}
	}

	delete [] coord;
	
	TString outputFileName = inputFileName;
	outputFileName.ReplaceAll("MCMCResponse_", "FinalResponse_");

	cout << outputFileName.Data() << endl;

	TFile *out = new TFile(outputFileName.Data(), "RECREATE");

	for (int cent = 0; cent < 3; cent++){
		fMeas[cent]->Write();
		fTrue[cent]->Write();
		fMiss[cent]->Write();
		SetName(fReco[cent], Form("Measured_%i", cent));
		fReco[cent]->Write();
		SetName(fMC[cent], Form("Truth_%i", cent));
		fMC[cent]->Write();
		gDirectory->WriteObject(resp[cent], Form("Resp_%i", cent));
	}

	out->Close();
	f->Close();
}

void CreateResponseMatrix(int SUPERITERATION = 1, int iteration = 3, TString DirName = "MCMCUnf"){
	Method(SUPERITERATION, iteration, DirName);
}
