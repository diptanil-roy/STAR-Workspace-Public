R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int SUPERITERATION = 1, int iteration = 4, bool mimicDataUncertainty = kFALSE){
	// for (int i = 0; i < njpt_gen_bins_var+1; i++){
	// 	jetpt_var_bin[i] = 5 + (20. - 5.)/njpt_gen_bins_var * i;
	// }

	TString RespFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), SUPERITERATION, 0, njpt_gen_bins_var, nz_gen_bins, 1);
	TFile *RespFile = new TFile(RespFileName.Data());
	cout << RespFileName.Data() << endl;
	RespFile->cd();

	RooUnfoldResponse *resp[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("Resp1D_%i", i), resp[i]);
	}

	TString FoldedFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins, 2); // Step 0 is always where I store the folded distro to be unfolded
	TFile *FoldedFile = new TFile(FoldedFileName.Data());
	cout << FoldedFileName.Data() << endl;

	THnSparseF *hJet[3];

	for (int i = 0; i < 3; i++){
		hJet[i] = (THnSparseF *)gDirectory->Get(Form("hJet_%i", i));
	}

	TH2D *Measured[3];
	TH2D *Truth[3];

	TH1D *Measured1D[3];
	TH1D *Truth1D[3];

	for (int i = 0; i < 3; i++){

		TH2D *tmpMeas = (TH2D *)hJet[i]->Projection(3, 2, "E");
		TH2D *tmpTrue = (TH2D *)hJet[i]->Projection(1, 0, "E");

		TH1D *tmpMeas1D = (TH1D *)hJet[i]->Projection(2, "E");
		TH1D *tmpTrue1D = (TH1D *)hJet[i]->Projection(0, "E");

		cout << i << "\t" << tmpMeas->Integral() << "\t" << tmpTrue->Integral() << endl;

		Truth[i] = new TH2D(Form("Truth_%i", i), Form("Truth_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
		Measured[i] = new TH2D(Form("Measured_%i", i), Form("Measured_%i", i), njpt_bins, nbinsjetpt, nz_bins, nbinsz);

		Truth1D[i] = new TH1D(Form("Truth1D_%i", i), Form("Truth1D_%i", i), njpt_gen_bins_var, jetpt_var_bin);
		Measured1D[i] = new TH1D(Form("Measured1D_%i", i), Form("Measured1D_%i", i), njpt_bins, nbinsjetpt);

		// Truth Distribution
		for (int binx = 1; binx <= Truth1D[i]->GetNbinsX(); binx++){
			double xlow = Truth1D[i]->GetXaxis()->GetBinLowEdge(binx);
			double xup = Truth1D[i]->GetXaxis()->GetBinUpEdge(binx);

			int binxlow = tmpTrue1D->GetXaxis()->FindBin(xlow + 0.00001);
			int binxup = tmpTrue1D->GetXaxis()->FindBin(xup - 0.00001);

			double error = 0.;
			double content = tmpTrue1D->IntegralAndError(binxlow, binxup, error, "");

			Truth1D[i]->SetBinContent(binx, content);
			Truth1D[i]->SetBinError(binx, error);
		}

		// Measured Distribution
		for (int binx = 1; binx <= Measured1D[i]->GetNbinsX(); binx++){
			double xlow = Measured1D[i]->GetXaxis()->GetBinLowEdge(binx);
			double xup = Measured1D[i]->GetXaxis()->GetBinUpEdge(binx);

			int binxlow = tmpMeas1D->GetXaxis()->FindBin(xlow + 0.00001);
			int binxup = tmpMeas1D->GetXaxis()->FindBin(xup - 0.00001);

			double error = 0.;
			double content = tmpMeas1D->IntegralAndError(binxlow, binxup, error);

			Measured1D[i]->SetBinContent(binx, content);
			Measured1D[i]->SetBinError(binx, error);
		}
	}

	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

	RooUnfoldBayes unfoldcent (resp[0], Measured1D[0], iteration);
	RooUnfoldBayes unfoldmid (resp[1], Measured1D[1], iteration);
	RooUnfoldBayes unfoldperi (resp[2], Measured1D[2], iteration);

	TH1D *Unfolded1D[3];

	Unfolded1D[0] = (TH1D *)unfoldcent.Hreco(errorTreatment);
	Unfolded1D[1] = (TH1D *)unfoldmid.Hreco(errorTreatment);
	Unfolded1D[2] = (TH1D *)unfoldperi.Hreco(errorTreatment);

	TCanvas *c[3];

	for (int i = 0; i < 3; i++){
		
		c[i] = new TCanvas(Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 800, 800);
		// c[i]->Divide(2);

		SetName(Measured1D[i], Form("MeasuredPt_%i", i));
		SetName(Truth1D[i], Form("TruthPt_%i", i));
		SetName(Unfolded1D[i], Form("UnfoldedPt_%i", i));
		SetColor(Measured1D[i], kBlue, 24);
		SetColor(Truth1D[i], kRed, 24);
		SetColor(Unfolded1D[i], kBlack, 20);

		c[i]->cd(1);
		Measured1D[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		Measured1D[i]->Draw("EP");
		Truth1D[i]->Draw("EP SAME");
		Unfolded1D[i]->Draw("EP SAME");

		gPad->SetLogy();
	}

	TString FinalOutName = RespFileName;
	FinalOutName.ReplaceAll("Response_", "Output_");
	FinalOutName.ReplaceAll("_Split_1", "");

	TH1D *RatioPt[3];

	TFile *out = new TFile(FinalOutName.Data(), "RECREATE");
	out->cd();

	for (int i = 0; i < 3; i++){
		Measured1D[i]->Write();
		Truth1D[i]->Write();
		Unfolded1D[i]->Write();

		RatioPt[i] = (TH1D *)Unfolded1D[i]->Clone(Form("RatioPt_%i", i));
		RatioPt[i]->Divide(Truth1D[i]);
		RatioPt[i]->Write();

		c[i]->Write();
	}
	
	out->Close();

}

void UnfoldMC(TString DirName = "MCMCUnf", int SUPERITERATION = 1, int iteration = 4, bool mimicDataUncertainty = kFALSE){
	Method(DirName, SUPERITERATION, iteration, mimicDataUncertainty);
}