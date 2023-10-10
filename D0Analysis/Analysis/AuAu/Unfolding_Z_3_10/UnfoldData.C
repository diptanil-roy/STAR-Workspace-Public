R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, bool simpleunfolding = kFALSE){
	TString RespFileName;
	if (!simpleunfolding) RespFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins_var, nz_gen_bins);
	else RespFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins);
	TFile *RespFile = new TFile(RespFileName.Data());
	RespFile->cd();

	RooUnfoldResponse *resp[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
	}

	TFile *DataFile = new TFile(Form("../SPlot_Mar2_2023/Histograms3_D0%i_10GeV_RecoJetPt_%i_%i_NewBinningOldCuts.root", D0pTLow, (int)nbinsjetpt[0], (int)nbinsjetpt[njpt_bins] ), "READ");
	DataFile->cd();

	TH2D *Measured[3];

	Measured[0] = (TH2D *)gDirectory->Get("ZPt_0_10");
	Measured[1] = (TH2D *)gDirectory->Get("ZPt_10_40");
	Measured[2] = (TH2D *)gDirectory->Get("ZPt_40_80");

	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

	RooUnfoldBayes unfoldcent (resp[0], Measured[0], iteration);
	RooUnfoldBayes unfoldmid (resp[1], Measured[1], iteration);
	RooUnfoldBayes unfoldperi (resp[2], Measured[2], iteration);

	TH2D *Unfolded[3];

	Unfolded[0] = (TH2D *)unfoldcent.Hreco(errorTreatment);
	Unfolded[1] = (TH2D *)unfoldmid.Hreco(errorTreatment);
	Unfolded[2] = (TH2D *)unfoldperi.Hreco(errorTreatment);

	TH1D *MeasuredPt[3];
	TH1D *MeasuredZ[3];

	TH1D *UnfoldedPt[3];
	TH1D *UnfoldedZ[3];

	TH1D *UnfoldedZ_JPtBins[3][njpt_gen_bins_var];

	TCanvas *c[3];

	for (int i = 0; i < 3; i++){
		MeasuredPt[i] = (TH1D *)Measured[i]->ProjectionX();
		MeasuredZ[i] = (TH1D *)Measured[i]->ProjectionY();

		UnfoldedPt[i] = (TH1D *)Unfolded[i]->ProjectionX();
		UnfoldedZ[i] = (TH1D *)Unfolded[i]->ProjectionY();

		for (int jptbin = 1; jptbin <= njpt_gen_bins_var; jptbin++){
			TH2D *h = (TH2D *)Unfolded[i]->Clone("tmp");
			h->GetXaxis()->SetRange(jptbin, jptbin);
			UnfoldedZ_JPtBins[i][jptbin-1] = (TH1D *)h->ProjectionY(Form("ZUnf_%i_%i", i, jptbin));	
			// UnfoldedZ_JPtBins[i][jptbin-1] = NULL;
		}

		c[i] = new TCanvas(Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 800, 800);
		c[i]->Divide(2);

		SetName(MeasuredPt[i], Form("MeasuredPt_%i", i));
		SetName(UnfoldedPt[i], Form("UnfoldedPt_%i", i));
		SetColor(MeasuredPt[i], kBlue, 24);
		SetColor(UnfoldedPt[i], kBlack, 20);

		c[i]->cd(1);
		MeasuredPt[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		MeasuredPt[i]->Draw("EP");
		UnfoldedPt[i]->Draw("EP SAME");

		gPad->SetLogy();

		SetName(MeasuredZ[i], Form("MeasuredZ_%i", i));
		SetName(UnfoldedZ[i], Form("UnfoldedZ_%i", i));
		SetColor(MeasuredZ[i], kBlue, 24);
		SetColor(UnfoldedZ[i], kBlack, 20);

		c[i]->cd(2);
		MeasuredZ[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		MeasuredZ[i]->Draw("EP");
		UnfoldedZ[i]->Draw("EP SAME");

		gPad->SetLogy();
	}

	TString FinalOutName = RespFileName;
	FinalOutName.ReplaceAll("Response_", "Output_");
	FinalOutName.ReplaceAll("_Split_1", "");

	TFile *out = new TFile(FinalOutName.Data(), "RECREATE");
	out->cd();

	for (int i = 0; i < 3; i++){
		SetName(Measured[i], Form("Measured_%i", i));
		Measured[i]->Write();
    	SetName(Unfolded[i], Form("Unfolded_%i", i));
		Unfolded[i]->Write();

		MeasuredPt[i]->Write();
		MeasuredZ[i]->Write();

		UnfoldedPt[i]->Write();
		UnfoldedZ[i]->Write();

		c[i]->Write();
	}
	for (int i = 0; i < 3; i++){
		for (int jptbin = 1; jptbin <= njpt_gen_bins_var; jptbin++){
			SetName(UnfoldedZ_JPtBins[i][jptbin-1], Form("ZUnf_%i_%i", i, jptbin));
			UnfoldedZ_JPtBins[i][jptbin-1]->Write();
		}
	}
	out->Close();

}

void UnfoldData(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, bool simpleunfolding = kFALSE){
	Method(DirName, D0pTLow, SUPERITERATION, iteration, simpleunfolding);
}