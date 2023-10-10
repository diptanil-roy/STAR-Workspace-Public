R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void ChangeCountOfHistogramByPercent(TH1 *h, double percent){
	for (int i = 1; i <= h->GetNbinsX(); i++){
		h->SetBinContent(i, h->GetBinContent(i) * (1 + percent/100.));
		h->SetBinError(i, h->GetBinError(i) * (1 + percent/100.));
	}
}

void ChangeCountOfHistogramByPercent(TH2 *h, double percent){
	for (int i = 1; i <= h->GetNbinsX(); i++){
		for (int j = 1; j <= h->GetNbinsY(); j++){
			double oldval = h->GetBinContent(i, j);
			double olderr = h->GetBinError(i, j);
			h->SetBinContent(i, j, oldval * (1 + percent/100.));
			h->SetBinError(i, j, olderr * (1 + percent/100.));
			cout << i << "\t" << j << "\t" << percent << "\t" << oldval << "\t" << olderr << "\t" << h->GetBinContent(i, j) << "\t" << h->GetBinError(i, j) << endl;
		}
	}
}

void RemoveNegativeBins(TH1 *h){
	for (int i = 1; i <= h->GetNbinsX(); i++){
		if (h->GetBinContent(i) < 0.) h->SetBinContent(i, 0.);
	}
}

void RemoveNegativeBins(TH2 *h){
	for (int i = 1; i <= h->GetNbinsX(); i++){
		for (int j = 1; j <= h->GetNbinsY(); j++){
			if (h->GetBinContent(i, j) < 0.) h->SetBinContent(i, j, 0.);
		}
	}
}

void Method(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, bool simpleunfolding = kTRUE, int d0efficiencycorrmode = 0){

	// d0efficiencycorrmode = 0 -> No change
	// d0efficiencycorrmode = 1 -> Req % changes in +ve direction
	// d0efficiencycorrmode = -1 -> Req % changes in -ve direction
	// d0efficiencycorrmode = 2 -> D0 reconstruction % changes in +ve direction
	// d0efficiencycorrmode = -2 -> D0 reconstruction % changes in -ve direction

	double D0CountFluctuationPercent[3][5] =   {{2.58, 	2.61, 	5.36, 	3.28, 	10.46},
												{0.95, 	1.92, 	1.44, 	3.94, 	4.33},
												{0.16, 	0.22, 	2.03, 	0.24, 	2.57}};
	
	TString RespFileName;
	if (!simpleunfolding) RespFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins_var, nz_gen_bins);
	else RespFileName = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins);
	
	TString RespFileName2;
	if (!simpleunfolding) RespFileName2 = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins_var, nz_gen_bins);
	else RespFileName2 = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins);

	cout << RespFileName2.Data() << endl;
	TFile *RespFile2 = new TFile(RespFileName2.Data());
	RespFile2->cd();

	RooUnfoldResponse *resp[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
	}

	cout << RespFileName.Data() << endl;
	TFile *RespFile = new TFile(RespFileName.Data());
	RespFile->cd();

	RooUnfoldResponse *resp1D[3];
	for (int i = 0; i < 3; i++){
		// gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
		gDirectory->GetObject(Form("Resp1D_%i", i), resp1D[i]);
	}

	RooUnfoldResponse *respwide[3];
	for (int i = 0; i < 3; i++){
		// gDirectory->GetObject(Form("Resp_%i", i), resp[i]);
		gDirectory->GetObject(Form("RespWide_%i", i), respwide[i]);
	}

	RooUnfoldResponse *respdR[3];
	for (int i = 0; i < 3; i++){
		gDirectory->GetObject(Form("RespWidedR_%i", i), respdR[i]);
	}

	TFile *DataFile;
	if (d0efficiencycorrmode == 2) DataFile = new TFile(Form("../Data/Aug17_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_%i_NoVCSysUp.root", D0pTLow, D0pTLow, 1000), "READ");
	else if (d0efficiencycorrmode == -2) DataFile = new TFile(Form("../Data/Aug17_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_%i_NoVCSysDown.root", D0pTLow, D0pTLow, 1000), "READ");
	else if (d0efficiencycorrmode == 3) DataFile = new TFile(Form("../Data/Aug17_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_%i_VCSysUp.root", D0pTLow, D0pTLow, 1000), "READ");
	else if (d0efficiencycorrmode == -3) DataFile = new TFile(Form("../Data/Aug17_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_%i_VCSysDown.root", D0pTLow, D0pTLow, 1000), "READ");
	else DataFile = new TFile(Form("../Data/Aug17_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_%i.root", D0pTLow, D0pTLow, 1000), "READ");
	DataFile->cd();

	TH2D *Measured[3];

	Measured[0] = (TH2D *)gDirectory->Get("ZPt_0_10");
	Measured[1] = (TH2D *)gDirectory->Get("ZPt_10_40");
	Measured[2] = (TH2D *)gDirectory->Get("ZPt_40_80");

	TH2D *MeasuredWide[3];

	MeasuredWide[0] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_0_10");
	MeasuredWide[1] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_10_40");
	MeasuredWide[2] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_40_80");

	TH2D *MeasuredWideOld[3];
	for (int cent = 0; cent < 3; cent++){
		MeasuredWideOld[cent] = (TH2D *)MeasuredWide[cent]->Clone();
		SetName(MeasuredWideOld[cent], Form("%s_%s", MeasuredWide[cent]->GetName(), "Old"));
	}

	if (d0efficiencycorrmode == 1){
		for (int i = 0; i < 3; i++){
			cout << i << "\t" << D0CountFluctuationPercent[i][4] << endl;
			ChangeCountOfHistogramByPercent(Measured[i], D0CountFluctuationPercent[i][4]);
			ChangeCountOfHistogramByPercent(MeasuredWide[i], D0CountFluctuationPercent[i][4]);
			RemoveNegativeBins(Measured[i]);
			RemoveNegativeBins(MeasuredWide[i]);
		}
	}
	else if (d0efficiencycorrmode == -1){
		for (int i = 0; i < 3; i++){
			cout << i << "\t" << D0CountFluctuationPercent[i][4] << endl;
			ChangeCountOfHistogramByPercent(Measured[i], -1 * D0CountFluctuationPercent[i][4]);
			ChangeCountOfHistogramByPercent(MeasuredWide[i], -1 * D0CountFluctuationPercent[i][4]);
			RemoveNegativeBins(Measured[i]);
			RemoveNegativeBins(MeasuredWide[i]);
		}
	}

	TH1D *Measured1D[3];

	Measured1D[0] = (TH1D *)MeasuredWide[0]->ProjectionX();
	Measured1D[1] = (TH1D *)MeasuredWide[1]->ProjectionX();
	Measured1D[2] = (TH1D *)MeasuredWide[2]->ProjectionX();

	// Measured1D[0] = (TH1D *)gDirectory->Get("JetPt_Area_Wide_0_10");
	// Measured1D[1] = (TH1D *)gDirectory->Get("JetPt_Area_Wide_10_40");
	// Measured1D[2] = (TH1D *)gDirectory->Get("JetPt_Area_Wide_40_80");

	SetName(Measured1D[0], Form("Measured1D_%i", 0));
	SetName(Measured1D[1], Form("Measured1D_%i", 1));
	SetName(Measured1D[2], Form("Measured1D_%i", 2));

	TH1D *Measured1DOld[3];
	for (int cent = 0; cent < 3; cent++){
		Measured1DOld[cent] = (TH1D *)Measured1D[cent]->Clone();
		SetName(Measured1DOld[cent], Form("%s_%s", Measured1D[cent]->GetName(), "Old"));
	}

	TH2D *Measured2DdRWide[3];

	Measured2DdRWide[0] = (TH2D *)gDirectory->Get("JetPt_R_Area_Wide_0_10");
	Measured2DdRWide[1] = (TH2D *)gDirectory->Get("JetPt_R_Area_Wide_10_40");
	Measured2DdRWide[2] = (TH2D *)gDirectory->Get("JetPt_R_Area_Wide_40_80");

	TH2D *Measured2DdRWideOld[3];

	for (int cent = 0; cent < 3; cent++){
		Measured2DdRWideOld[cent] = (TH2D *)Measured2DdRWide[cent]->Clone();
		SetName(Measured2DdRWideOld[cent], Form("%s_%s", Measured2DdRWide[cent]->GetName(), "Old"));
	}

	if (d0efficiencycorrmode == 1){
		for (int i = 0; i < 3; i++){
			ChangeCountOfHistogramByPercent(Measured2DdRWide[i], D0CountFluctuationPercent[i][4]);
			RemoveNegativeBins(Measured2DdRWide[i]);
		}
	}
	else if (d0efficiencycorrmode == -1){
		for (int i = 0; i < 3; i++){
			ChangeCountOfHistogramByPercent(Measured2DdRWide[i], -1 * D0CountFluctuationPercent[i][4]);
			RemoveNegativeBins(Measured2DdRWide[i]);
		}
	}

	for (int i = 0; i < 3; i++){
		cout << i << "\t" << Measured[i]->Integral() << endl;
		cout << i << "\t" << Measured1D[i]->Integral() << endl;
		cout << i << "\t" << MeasuredWide[i]->Integral() << endl;
		cout << i << "\t" << Measured2DdRWide[i]->Integral() << endl;
	}

	// for (int i = 0; i < 3; i++){
	// 	for (int binx = 0; binx <= Measured[i]->GetNbinsX()+1; binx++){
	// 		for (int biny = 0; biny <= Measured[i]->GetNbinsY()+1; biny++){
	// 			if (Measured[i]->GetBinContent(binx, biny) < 0.) Measured[i]->SetBinContent(binx, biny, 0);
	// 		}
	// 	}

	// 	cout << i << "\t" << Measured[i]->Integral() << endl;
	// }

	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

	RooUnfoldBayes unfoldcent (resp[0], Measured[0], iteration);
	RooUnfoldBayes unfoldmid (resp[1], Measured[1], iteration);
	RooUnfoldBayes unfoldperi (resp[2], Measured[2], iteration);

	RooUnfoldBayes unfold1Dcent (resp1D[0], Measured1D[0], iteration);
	RooUnfoldBayes unfold1Dmid (resp1D[1], Measured1D[1], iteration);
	RooUnfoldBayes unfold1Dperi (resp1D[2], Measured1D[2], iteration);

	RooUnfoldBayes unfoldwidecent (respwide[0], MeasuredWide[0], iteration);
	RooUnfoldBayes unfoldwidemid (respwide[1], MeasuredWide[1], iteration);
	RooUnfoldBayes unfoldwideperi (respwide[2], MeasuredWide[2], iteration);

	RooUnfoldBayes unfoldwidecentdR (respdR[0], Measured2DdRWide[0], iteration);
	RooUnfoldBayes unfoldwidemiddR (respdR[1], Measured2DdRWide[1], iteration);
	RooUnfoldBayes unfoldwideperidR (respdR[2], Measured2DdRWide[2], iteration);

	cout << "Called RooUnfoldBayes Object" << endl;

	TH2D *Unfolded[3];

	Unfolded[0] = (TH2D *)unfoldcent.Hreco(errorTreatment);
	Unfolded[1] = (TH2D *)unfoldmid.Hreco(errorTreatment);
	Unfolded[2] = (TH2D *)unfoldperi.Hreco(errorTreatment);

	TH1D *Unfolded1D[3];

	Unfolded1D[0] = (TH1D *)unfold1Dcent.Hreco(errorTreatment);
	Unfolded1D[1] = (TH1D *)unfold1Dmid.Hreco(errorTreatment);
	Unfolded1D[2] = (TH1D *)unfold1Dperi.Hreco(errorTreatment);

	TH2D *UnfoldedWide[3];

	UnfoldedWide[0] = (TH2D *)unfoldwidecent.Hreco(errorTreatment);
	UnfoldedWide[1] = (TH2D *)unfoldwidemid.Hreco(errorTreatment);
	UnfoldedWide[2] = (TH2D *)unfoldwideperi.Hreco(errorTreatment);

	TH2D *Unfolded2DdRWide[3];

	Unfolded2DdRWide[0] = (TH2D *)unfoldwidecentdR.Hreco(errorTreatment);
	Unfolded2DdRWide[1] = (TH2D *)unfoldwidemiddR.Hreco(errorTreatment);
	Unfolded2DdRWide[2] = (TH2D *)unfoldwideperidR.Hreco(errorTreatment);

	TH1D *MeasuredPt[3];
	TH1D *MeasuredZ[3];
	TH1D *UnfoldedPt[3];
	TH1D *UnfoldedZ[3];

	TH1D *MeasuredPtWide[3];
	TH1D *MeasuredZWide[3];
	TH1D *MeasureddRWide[3];

	TH1D *MeasuredPtWideOld[3];
	TH1D *MeasuredZWideOld[3];
	TH1D *MeasureddRWideOld[3];

	TH1D *UnfoldedPtWide[3];
	TH1D *UnfoldedZWide[3];
	TH1D *UnfoldeddRWide[3];

	TH1D *UnfoldedZ_JPtBins[3][njpt_gen_bins_var];

	TH1D *DiffPt[3];
	TH1D *DiffZ[3];

	cout << "Creating Histograms Here" << endl;

	TCanvas *c[3];
	TCanvas *d[3];
	TCanvas *d2[3];

	for (int i = 0; i < 3; i++){
		MeasuredPt[i] = (TH1D *)Measured[i]->ProjectionX();
		MeasuredZ[i] = (TH1D *)Measured[i]->ProjectionY();

		UnfoldedPt[i] = (TH1D *)Unfolded[i]->ProjectionX();
		UnfoldedZ[i] = (TH1D *)Unfolded[i]->ProjectionY();

		MeasuredPtWide[i] = (TH1D *)MeasuredWide[i]->ProjectionX();
		MeasuredZWide[i] = (TH1D *)MeasuredWide[i]->ProjectionY();
		MeasureddRWide[i] = (TH1D *)Measured2DdRWide[i]->ProjectionY();

		MeasuredPtWideOld[i] = (TH1D *)MeasuredWideOld[i]->ProjectionX();
		MeasuredZWideOld[i] = (TH1D *)MeasuredWideOld[i]->ProjectionY();
		MeasureddRWideOld[i] = (TH1D *)Measured2DdRWideOld[i]->ProjectionY();

		DiffPt[i] = (TH1D *)MeasuredPtWide[i]->Clone();
		DiffPt[i]->Add(MeasuredPtWideOld[i], -1);
		DiffPt[i]->Divide(MeasuredPtWideOld[i]);
		DiffPt[i]->Scale(100);

		DiffZ[i] = (TH1D *)MeasuredZWide[i]->Clone();
		DiffZ[i]->Add(MeasuredZWideOld[i], -1);
		DiffZ[i]->Divide(MeasuredZWideOld[i]);
		DiffZ[i]->Scale(100);

		UnfoldedPtWide[i] = (TH1D *)UnfoldedWide[i]->ProjectionX();
		UnfoldedZWide[i] = (TH1D *)UnfoldedWide[i]->ProjectionY();
		UnfoldeddRWide[i] = (TH1D *)Unfolded2DdRWide[i]->ProjectionY();

		cout << "Made Histos here" << endl;

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

		d[i] = new TCanvas(Form("PlotsWide_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("Plots_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 800, 800);
		d[i]->Divide(3);
		d[i]->cd(1);

		SetName(MeasuredPtWide[i], Form("MeasuredPtWide_%i", i));
		SetName(UnfoldedPtWide[i], Form("UnfoldedPtWide_%i", i));
		SetColor(MeasuredPtWide[i], kBlue, 24);
		SetColor(UnfoldedPtWide[i], kBlack, 20);
		SetColor(MeasuredPtWideOld[i], kGreen-2, 20);
		SetName(Measured1D[i], Form("Measured1D_%i", i));
		SetName(Unfolded1D[i], Form("Unfolded1D_%i", i));
		SetColor(Measured1D[i], kBlue, 25);
		SetColor(Unfolded1D[i], kBlack, 21);
		
		MeasuredPtWide[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		MeasuredPtWide[i]->Draw("EP");
		UnfoldedPtWide[i]->Draw("EP SAME");
		MeasuredPtWideOld[i]->Draw("EP SAME");
		Measured1D[i]->Draw("EP SAME");
		Unfolded1D[i]->Draw("EP SAME");

		gPad->SetLogy();

		SetName(MeasuredZWide[i], Form("MeasuredZWide_%i", i));
		SetName(UnfoldedZWide[i], Form("UnfoldedZWide_%i", i));
		SetColor(MeasuredZWide[i], kBlue, 24);
		SetColor(UnfoldedZWide[i], kBlack, 20);
		SetColor(MeasuredZWideOld[i], kGreen-2, 20);

		d[i]->cd(2);
		MeasuredZWide[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		MeasuredZWide[i]->Draw("EP");
		MeasuredZWideOld[i]->Draw("EP SAME");
		UnfoldedZWide[i]->Draw("EP SAME");

		gPad->SetLogy();

		SetName(MeasureddRWide[i], Form("MeasureddRWide_%i", i));
		SetName(UnfoldeddRWide[i], Form("UnfoldeddRWide_%i", i));
		SetColor(MeasureddRWide[i], kBlue, 24);
		SetColor(UnfoldeddRWide[i], kBlack, 20);

		d[i]->cd(3);
		MeasureddRWide[i]->GetYaxis()->SetRangeUser(pow(10, 0), pow(10, 7));
		MeasureddRWide[i]->Draw("EP");
		UnfoldeddRWide[i]->Draw("EP SAME");

		gPad->SetLogy();

		d2[i] = new TCanvas(Form("Diff_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), Form("Diff_Step_%i_Iter_%i_Cent_%i", SUPERITERATION, iteration, i), 800, 800);
		d2[i]->Divide(2);
		d2[i]->cd(1);
		DiffPt[i]->Draw("HIST");
		d2[i]->cd(2);
		DiffZ[i]->Draw("HIST");

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

		SetName(Measured1D[i], Form("Measured1D_%i", i));
		Measured1D[i]->Write();
		SetName(Unfolded1D[i], Form("Unfolded1D_%i", i));
		Unfolded1D[i]->Write();

		SetName(MeasuredWide[i], Form("MeasuredWide_%i", i));
		MeasuredWide[i]->Write();
    	SetName(UnfoldedWide[i], Form("UnfoldedWide_%i", i));
		UnfoldedWide[i]->Write();

		SetName(Measured2DdRWide[i], Form("Measured2DdRWide_%i", i));
		Measured2DdRWide[i]->Write();
		SetName(Unfolded2DdRWide[i], Form("Unfolded2DdRWide_%i", i));
		Unfolded2DdRWide[i]->Write();

		MeasuredPt[i]->Write();
		MeasuredZ[i]->Write();

		UnfoldedPt[i]->Write();
		UnfoldedZ[i]->Write();

		MeasuredPtWide[i]->Write();
		MeasuredZWide[i]->Write();
		MeasureddRWide[i]->Write();

		UnfoldedPtWide[i]->Write();
		UnfoldedZWide[i]->Write();
		UnfoldeddRWide[i]->Write();

		c[i]->Write();
		d[i]->Write();
	}
	// for (int i = 0; i < 3; i++){
	// 	for (int jptbin = 1; jptbin <= njpt_gen_bins_var; jptbin++){
	// 		SetName(UnfoldedZ_JPtBins[i][jptbin-1], Form("ZUnf_%i_%i", i, jptbin));
	// 		UnfoldedZ_JPtBins[i][jptbin-1]->Write();
	// 	}
	// }
	out->Close();

}

void UnfoldData(TString DirName = "MCMCUnf", int D0pTLow = 1, int SUPERITERATION = 1, int iteration = 4, bool simpleunfolding = kFALSE, int d0efficiencycorrmode = 0){
	Method(DirName, D0pTLow, SUPERITERATION, iteration, simpleunfolding, d0efficiencycorrmode);
}