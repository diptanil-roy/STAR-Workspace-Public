R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, bool simpleunfolding = kFALSE){
	TString RespFileName;
	RespFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins);
	
	cout << RespFileName.Data() << endl;
	TFile *RespFile = new TFile(RespFileName.Data());
	RespFile->cd();

	RooUnfoldResponse *resp[3];
	RooUnfoldResponse *resp1D[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
		gDirectory->GetObject(Form("Resp1D_%i", i), resp1D[i]);
	}

	TFile *DataFile = new TFile(Form("../SPlotUpdated/AreaBased/Histograms3_D0%i_10GeV_RecoJetPt_%i_%i_NewBinningOldCuts.root", D0pTLow, 0, (int)nbinsjetpt[njpt_bins] ), "READ");
	DataFile->cd();

	TH1D *Measured1D[3];

	Measured1D[0] = (TH1D *)gDirectory->Get("JetPt_0_10");
	Measured1D[1] = (TH1D *)gDirectory->Get("JetPt_10_40");
	Measured1D[2] = (TH1D *)gDirectory->Get("JetPt_40_80");

	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

	RooUnfoldBayes unfoldcent1D (resp1D[0], Measured1D[0], iteration);
	RooUnfoldBayes unfoldmid1D (resp1D[1], Measured1D[1], iteration);
	RooUnfoldBayes unfoldperi1D (resp1D[2], Measured1D[2], iteration);

	TH1D *Unfolded1D[3];

	Unfolded1D[0] = (TH1D *)unfoldcent1D.Hreco(errorTreatment);
	Unfolded1D[1] = (TH1D *)unfoldmid1D.Hreco(errorTreatment);
	Unfolded1D[2] = (TH1D *)unfoldperi1D.Hreco(errorTreatment);

	for (int i = 0; i < 3; i++){
		SetName(Measured1D[i], Form("Measured1D_%i", i));
		SetName(Unfolded1D[i], Form("Unfolded1D_%i", i));
		SetColor(Measured1D[i], kBlue, 24);
		SetColor(Unfolded1D[i], kBlack, 20);
	}

	TString FinalOutName = RespFileName;
	FinalOutName.ReplaceAll("Response_", "Output_");
	FinalOutName.ReplaceAll("_Split_1", "");

	TFile *out = new TFile(FinalOutName.Data(), "RECREATE");
	out->cd();

	for (int i = 0; i < 3; i++){
		
		SetName(Measured1D[i], Form("Measured1D_%i", i));
		SetName(Unfolded1D[i], Form("Unfolded1D_%i", i));
		SetColor(Measured1D[i], kBlue, 24);
		SetColor(Unfolded1D[i], kBlack, 20);

		Measured1D[i]->Write();
		Unfolded1D[i]->Write();
	}

	out->Close();

}

void UnfoldData(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, bool simpleunfolding = kFALSE){
	Method(DirName, D0pTLow, SUPERITERATION, iteration, simpleunfolding);
}