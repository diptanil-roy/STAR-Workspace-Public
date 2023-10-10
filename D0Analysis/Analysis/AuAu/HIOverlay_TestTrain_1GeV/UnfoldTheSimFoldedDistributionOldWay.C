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
#include "TRatioPlot.h"
#include <vector>

#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldSvd.h"

#pragma link C++ class vector<int> +;

using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

void Method(int SUPERITERATION = 0, int iteration = 3, TString DirName = "DataMCUnf", bool oldResponse = kFALSE){
	
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TString RespFileName = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins, 1);
	TFile *RespFile = new TFile(RespFileName.Data());

	cout << RespFileName.Data() << endl;
	RespFile->cd();

	RooUnfoldResponse *resp[3];
	RooUnfoldResponse *resp1D[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
		gDirectory->GetObject(Form("Resp1D_%i", i), resp1D[i]);
	}

	TString FoldedFileName = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins, 2); // Step 0 is always where I store the folded distro to be unfolded
	TFile *FoldedFile = new TFile(FoldedFileName.Data());
	// FoldedFile->cd("Response");
	// gDirectory->ls();

	cout << FoldedFileName.Data() << endl;

	cout << resp[0]->GetDimensionMeasured() << endl;

	TH1D *Measured1D[3];
	TH1D *Truth1D[3];

	Measured1D[0] = (TH1D *)gDirectory->Get("fMeas1D_Cent_0");
	Measured1D[1] = (TH1D *)gDirectory->Get("fMeas1D_Cent_1");
	Measured1D[2] = (TH1D *)gDirectory->Get("fMeas1D_Cent_2");

	Truth1D[0] = (TH1D *)gDirectory->Get("fTrue1D_Cent_0");
	Truth1D[1] = (TH1D *)gDirectory->Get("fTrue1D_Cent_1");
	Truth1D[2] = (TH1D *)gDirectory->Get("fTrue1D_Cent_2");


	cout << "Integral = " << Measured1D[0]->Integral() << "\t" << Measured1D[1]->Integral() << "\t" << Measured1D[2]->Integral() << endl;

	TH2D *Measured[3];
	TH2D *Truth[3];

	Measured[0] = (TH2D *)gDirectory->Get("fMeas_Cent_0");
	Measured[1] = (TH2D *)gDirectory->Get("fMeas_Cent_1");
	Measured[2] = (TH2D *)gDirectory->Get("fMeas_Cent_2");

	Truth[0] = (TH2D *)gDirectory->Get("fTrue_Cent_0");
	Truth[1] = (TH2D *)gDirectory->Get("fTrue_Cent_1");
	Truth[2] = (TH2D *)gDirectory->Get("fTrue_Cent_2");

	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

	RooUnfoldBayes unfoldcent (resp[0], Measured[0], iteration);
	RooUnfoldBayes unfoldmid (resp[1], Measured[1], iteration);
	RooUnfoldBayes unfoldperi (resp[2], Measured[2], iteration);

	TH2D *Unfolded[3];

	Unfolded[0] = (TH2D *)unfoldcent.Hreco(errorTreatment);
	Unfolded[1] = (TH2D *)unfoldmid.Hreco(errorTreatment);
	Unfolded[2] = (TH2D *)unfoldperi.Hreco(errorTreatment);

	RooUnfoldBayes unfold1Dcent (resp1D[0], Measured1D[0], iteration);
	RooUnfoldBayes unfold1Dmid (resp1D[1], Measured1D[1], iteration);
	RooUnfoldBayes unfold1Dperi (resp1D[2], Measured1D[2], iteration);

	TH1D *Unfolded1D[3];

	Unfolded1D[0] = (TH1D *)unfold1Dcent.Hreco(errorTreatment);
	Unfolded1D[1] = (TH1D *)unfold1Dmid.Hreco(errorTreatment);
	Unfolded1D[2] = (TH1D *)unfold1Dperi.Hreco(errorTreatment);

	TH1D *MeasuredPt[3];
	TH1D *MeasuredZ[3];

	TH1D *TruthPt[3];
	TH1D *TruthZ[3];

	TH1D *UnfoldedPt[3];
	TH1D *UnfoldedZ[3];

	TH1D *UnfoldedZ_JPtBins[3][njpt_gen_bins];

	TCanvas *c[3];
	TCanvas *d[3][3];

	TRatioPlot *rp1[3];
	TRatioPlot *rp2[3];
	TRatioPlot *rp3[3];

	for (int i = 0; i < 3; i++){
		MeasuredPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)Measured[i]->ProjectionX());
		MeasuredZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)Measured[i]->ProjectionY());

		UnfoldedPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)Unfolded[i]->ProjectionX());
		UnfoldedZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)Unfolded[i]->ProjectionY());

		TruthPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)Truth[i]->ProjectionX());
		TruthZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)Truth[i]->ProjectionY());

		c[i] = new TCanvas(Form("Closure_Cent_%i", i), Form("Closure_Cent_%i", i), 1800, 800);
		c[i]->Divide(3);

		SetAxisTitles(TruthPt[i], "p_{T,Jet} [Gev/#it{c}]", "Counts");
		SetAxisTitles(UnfoldedPt[i], "p_{T,Jet} [Gev/#it{c}]", "Counts");

		SetName(TruthPt[i], Form("TruthPt_%i", i));
		SetName(UnfoldedPt[i], Form("UnfoldedPt_%i", i));
		SetColor(TruthPt[i], kRed, 24);
		SetColor(UnfoldedPt[i], kBlack, 20);

		c[i]->cd(1);
		rp1[i] = new TRatioPlot(UnfoldedPt[i], TruthPt[i]);
		rp1[i]->Draw();
		rp1[i]->GetUpperRefYaxis()->SetRangeUser(0.5, 5*pow(10, 5));
		rp1[i]->GetUpperPad()->SetLogy(1);
		rp1[i]->GetLowerRefYaxis()->SetRangeUser(0.5, 1.2);
		rp1[i]->GetLowerPad()->SetLogy(0);
		
		SetAxisTitles(TruthZ[i], "z = #frac{p_{T,D^{0}}}{p_{T,Jet}}", "Counts");
		SetAxisTitles(UnfoldedZ[i], "z = #frac{p_{T,D^{0}}}{p_{T,Jet}}", "Counts");

		SetName(TruthZ[i], Form("TruthZ_%i", i));
		SetName(UnfoldedZ[i], Form("UnfoldedZ_%i", i));
		SetColor(TruthZ[i], kRed, 24);
		SetColor(UnfoldedZ[i], kBlack, 20);

		c[i]->cd(2);
		rp2[i] = new TRatioPlot(UnfoldedZ[i], TruthZ[i]);
		rp2[i]->Draw();
		rp2[i]->GetUpperRefYaxis()->SetRangeUser(5*pow(10, 1), 5*pow(10, 7));
		rp2[i]->GetUpperPad()->SetLogy();
		rp2[i]->GetLowerRefYaxis()->SetRangeUser(0.5, 1.2);
		rp2[i]->GetLowerPad()->SetLogy(0);

		Unfolded1D[i] = (TH1D *)ProcessSpectraHistogram(Unfolded1D[i]);
		Truth1D[i] = (TH1D *)ProcessSpectraHistogram(Truth1D[i]);

		SetAxisTitles(Truth1D[i], "p_{T,Jet} [Gev/#it{c}]", "Counts");
		SetAxisTitles(Unfolded1D[i], "p_{T,Jet} [Gev/#it{c}]", "Counts");

		SetName(Truth1D[i], Form("Truth1D_%i", i));
		SetName(Unfolded1D[i], Form("Unfolded1D_%i", i));
		SetColor(Truth1D[i], kBlue, 24);
		SetColor(Unfolded1D[i], kBlack, 20);

		c[i]->cd(3);
		rp3[i] = new TRatioPlot(Unfolded1D[i], Truth1D[i]);
		rp3[i]->Draw();
		rp3[i]->GetUpperRefYaxis()->SetRangeUser(0.5, 5*pow(10, 5));
		rp3[i]->GetUpperPad()->SetLogy();
		rp3[i]->GetLowerRefYaxis()->SetRangeUser(0.5, 1.2);
		rp3[i]->GetLowerPad()->SetLogy(0);

		c[i]->cd(1);
		gPad = (TPad *)rp1[i]->GetUpperPad();
		TH1D *tmp = (TH1D *)Truth1D[i]->Clone();
		SetColor(tmp, kViolet, 25);
		TH1D *tmp2 = (TH1D *)Unfolded1D[i]->Clone();
		SetColor(tmp2, kBlack, 25);
		tmp->Draw("HIST SAME");

		c[i]->SaveAs(Form("%s/Closure_%i.pdf", DirName.Data(), i));
		// delete c[i];
	}

	for (int i = 0; i < 3; i++){
		d[i][0] = new TCanvas(Form("Closure_Cent_%i_2D", i), Form("Closure_Cent_%i_2D", i), 1200, 600);
		d[i][0]->Divide(2);

		d[i][0]->cd(1);
		rp1[i]->Draw();
		rp1[i]->GetUpperRefYaxis()->SetRangeUser(0.5, 5*pow(10, 5));
		rp1[i]->GetUpperPad()->SetLogy(1);
		rp1[i]->GetLowerRefYaxis()->SetRangeUser(0.5, 1.2);
		rp1[i]->GetLowerRefYaxis()->SetNdivisions(505);
		rp1[i]->GetLowerRefXaxis()->SetTitle("p_{T,Jet} [Gev/#it{c}]");
		rp1[i]->GetLowerPad()->SetLogy(0);
		d[i][0]->cd(2);
		rp2[i]->Draw();
		rp2[i]->GetUpperRefYaxis()->SetRangeUser(5*pow(10, 1), 5*pow(10, 7));
		rp2[i]->GetUpperPad()->SetLogy();
		rp2[i]->GetLowerRefYaxis()->SetRangeUser(0.5, 1.2);
		rp2[i]->GetLowerRefXaxis()->SetTitle("z = #frac{p_{T,D^{0}}}{p_{T,Jet}}");
		rp2[i]->GetLowerRefYaxis()->SetNdivisions(504);
		rp2[i]->GetLowerPad()->SetLogy(0);

		d[i][0]->SaveAs(Form("%s/Closure_%i_2D.pdf", DirName.Data(), i));
		d[i][0]->SaveAs(Form("%s/Closure_%i_2D.png", DirName.Data(), i));
		// delete d[i][0];

		d[i][1] = new TCanvas(Form("Closure_Cent_%i_1D", i), Form("Closure_Cent_%i_1D", i), 600, 600);

		d[i][1]->cd();
		rp3[i]->Draw();
		rp3[i]->GetUpperRefYaxis()->SetRangeUser(0.5, 5*pow(10, 5));
		rp3[i]->GetUpperPad()->SetLogy();
		rp3[i]->GetLowerRefYaxis()->SetRangeUser(0.5, 1.2);
		rp3[i]->GetLowerRefYaxis()->SetNdivisions(504);
		rp3[i]->GetLowerRefXaxis()->SetTitle("p_{T,Jet} [Gev/#it{c}]");
		rp3[i]->GetLowerPad()->SetLogy(0);
		
		d[i][1]->SaveAs(Form("%s/Closure_%i_1D.pdf", DirName.Data(), i));
		d[i][1]->SaveAs(Form("%s/Closure_%i_1D.png", DirName.Data(), i));
		// delete d[i][1];
	}

	/*

	TString FinalOutName = RespFileName;
	if (!oldResponse) FinalOutName.ReplaceAll("FinalResponse_", "Output_");
	else FinalOutName.ReplaceAll("OldResponse_", "OldOutput_");

	cout << FinalOutName.Data() << endl;

	TFile *out = new TFile(FinalOutName.Data(), "RECREATE");
	out->cd();

	for (int i = 0; i < 3; i++){
		Measured[i]->Write();
		// Truth[i]->Write();
    	SetName(Unfolded[i], Form("Unfolded_%i", i));
		Unfolded[i]->Write();

		c[i]->Write();
	}

	for (int i = 0; i < 3; i++){
		Measured1D[i]->Write();
		// Truth[i]->Write();
    	SetName(Unfolded1D[i], Form("Unfolded1D_%i", i));
		Unfolded1D[i]->Write();

		d[i]->Write();
	}

	// for (int i = 0; i < 3; i++){
	// 	for (int jptbin = 1; jptbin <= njpt_gen_bins; jptbin++){
	// 		SetName(UnfoldedZ_JPtBins[i][jptbin-1], Form("ZUnf_%i_%i", i, jptbin));
	// 		UnfoldedZ_JPtBins[i][jptbin-1]->Write();
	// 	}
	// }
	out->Close();

	*/
}

void UnfoldTheSimFoldedDistributionOldWay(int SUPERITERATION = 1, int iteration = 3, TString DirName = "SimMCUnf", bool oldResponse = kFALSE){
	Method(SUPERITERATION, iteration, DirName, oldResponse);
}