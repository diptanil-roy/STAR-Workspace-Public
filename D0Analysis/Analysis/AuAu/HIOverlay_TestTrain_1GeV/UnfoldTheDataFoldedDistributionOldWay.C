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

#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldSvd.h"

#pragma link C++ class vector<int> +;

using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

void Method(int SUPERITERATION = 0, int iteration = 3, TString DirName = "DataMCUnf", bool oldResponse = kFALSE, TString ResponseAppend = "Resp1D"){
	TString RespFileName = Form("%s/FinalResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins);
	if (oldResponse) RespFileName = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins);
	TFile *RespFile = new TFile(RespFileName.Data());

	cout << RespFileName.Data() << endl;
	RespFile->cd();

	RooUnfoldResponse *resp[3];
	RooUnfoldResponse *resp1D[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
		gDirectory->GetObject(Form("%s_%i", ResponseAppend.Data(), i), resp1D[i]);
	}

	TString FoldedFileName = Form("../SPlot_Mar2_2023/Histograms3_D01_10GeV_NewBinningOldCuts_Apr3.root"); // Step 0 is always where I store the folded distro to be unfolded
	TFile *FoldedFile = new TFile(FoldedFileName.Data());
	// FoldedFile->cd("Response");
	// gDirectory->ls();

	cout << FoldedFileName.Data() << endl;

	TH1D *Measured1D[3];
	// TH2D *Truth[3];

	Measured1D[0] = (TH1D *)gDirectory->Get("JetPt_0_10");
	Measured1D[1] = (TH1D *)gDirectory->Get("JetPt_10_40");
	Measured1D[2] = (TH1D *)gDirectory->Get("JetPt_40_80");

	cout << "Integral = " << Measured1D[0]->Integral() << "\t" << Measured1D[1]->Integral() << "\t" << Measured1D[2]->Integral() << endl;

	TH2D *Measured[3];
	// TH2D *Truth[3];

	Measured[0] = (TH2D *)gDirectory->Get("ZPt_0_10");
	Measured[1] = (TH2D *)gDirectory->Get("ZPt_10_40");
	Measured[2] = (TH2D *)gDirectory->Get("ZPt_40_80");

	cout << "Integral = " << Measured[0]->Integral() << "\t" << Measured[1]->Integral() << "\t" << Measured[2]->Integral() << endl;

	// for (int i = 0; i < 3; i++){
	// 	Measured[i] = (TH2D *)gDirectory->Get(Form("Measured_%i", i));
	// 	Truth[i] = (TH2D *)gDirectory->Get(Form("Truth_%i", i));
	// 	cout << "Step " << i << endl;
	// 	cout << Measured[i]->GetNbinsX() << "\t" << Truth[i]->GetNbinsX() << endl;

	// }


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


	cout << "Unfolded the distributions" << endl;

	TH1D *MeasuredPt[3];
	TH1D *MeasuredZ[3];

	// TH1D *TruthPt[3];
	// TH1D *TruthZ[3];

	TH1D *UnfoldedPt[3];
	TH1D *UnfoldedZ[3];

	TH1D *UnfoldedZ_JPtBins[3][njpt_gen_bins];

	TCanvas *c[3];
	TCanvas *c2[3];
	TCanvas *d[3];

	for (int i = 0; i < 3; i++){
		MeasuredPt[i] = (TH1D *)Measured[i]->ProjectionX();
		MeasuredZ[i] = (TH1D *)Measured[i]->ProjectionY();

		// TruthPt[i] = (TH1D *)Truth[i]->ProjectionX();
		// TruthZ[i] = (TH1D *)Truth[i]->ProjectionY();

		UnfoldedPt[i] = (TH1D *)Unfolded[i]->ProjectionX();
		UnfoldedZ[i] = (TH1D *)Unfolded[i]->ProjectionY();

		for (int jptbin = 1; jptbin <= njpt_gen_bins; jptbin++){
			TH2D *h = (TH2D *)Unfolded[i]->Clone("tmp");
			h->GetXaxis()->SetRange(jptbin, jptbin);
			UnfoldedZ_JPtBins[i][jptbin-1] = (TH1D *)h->ProjectionY(Form("ZUnf_%i_%i", i, jptbin));	
			// UnfoldedZ_JPtBins[i][jptbin-1] = NULL;
		}

		c[i] = new TCanvas(Form("Plots_2D_Pt_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("Plots_2D_Pt_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 600, 600);
		// c[i]->Divide(2);

		SetName(MeasuredPt[i], Form("MeasuredPt_%i", i));
		// SetName(TruthPt[i], Form("TruthPt_%i", i));
		SetName(UnfoldedPt[i], Form("UnfoldedPt_%i", i));
		SetColor(MeasuredPt[i], kBlue, 24);
		// SetColor(TruthPt[i], kRed, 24);
		SetColor(UnfoldedPt[i], kBlack, 20);

		c[i]->cd(1);
		MeasuredPt[i]->Draw("EP");
		MeasuredPt[i]->GetYaxis()->SetRangeUser(pow(10, -2), pow(10, 5));
		// TruthPt[i]->Draw("EP SAME");
		UnfoldedPt[i]->Draw("EP SAME");


		gPad->SetLogy();

		SetName(MeasuredZ[i], Form("MeasuredZ_%i", i));
		// SetName(TruthZ[i], Form("TruthZ_%i", i));
		SetName(UnfoldedZ[i], Form("UnfoldedZ_%i", i));
		SetColor(MeasuredZ[i], kBlue, 24);
		// SetColor(TruthZ[i], kRed, 24);
		SetColor(UnfoldedZ[i], kBlack, 20);

		c[i]->SaveAs(Form("%s/UnfoldedPlots_Pt_2D_Cent_%i.pdf", DirName.Data(), i));

		c2[i] = new TCanvas(Form("Plots_2D_Z_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("Plots_2D_Z_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 600, 600);
		c2[i]->cd(1);
		MeasuredZ[i]->Draw("EP");
		MeasuredZ[i]->GetYaxis()->SetRangeUser(pow(10, -2), pow(10, 5));
		// TruthZ[i]->Draw("EP SAME");
		UnfoldedZ[i]->Draw("EP SAME");

		gPad->SetLogy();

		c2[i]->SaveAs(Form("%s/UnfoldedPlots_Z_2D_Cent_%i.pdf", DirName.Data(), i));

		d[i] = new TCanvas(Form("1DPlots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("1DPlots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 600, 600);
		d[i]->cd();

		SetName(Measured1D[i], Form("Measured1D_%i", i));
		SetName(Unfolded1D[i], Form("Unfolded1D_%i", i));
		SetColor(Measured1D[i], kBlue, 24);
		SetColor(Unfolded1D[i], kBlack, 20);

		d[i]->cd();
		Measured1D[i]->Draw("EP");
		Measured1D[i]->GetYaxis()->SetRangeUser(pow(10, -2), pow(10, 5));
		Unfolded1D[i]->Draw("EP SAME");

		gPad->SetLogy();

		d[i]->SaveAs(Form("%s/UnfoldedPlots_Pt_1D_Cent_%i.pdf", DirName.Data(), i));
	}

	

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


}

void UnfoldTheDataFoldedDistributionOldWay(int SUPERITERATION = 1, int iteration = 3, TString DirName = "DataMCUnf", bool oldResponse = kFALSE, TString ResponseAppend = "Resp1D"){
	Method(SUPERITERATION, iteration, DirName, oldResponse, ResponseAppend);
}