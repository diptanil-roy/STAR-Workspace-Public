using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

double Chi2Func(TH1D *h1, TH1D *h2){ //h1 is experimental, h2 is truth
	assert (h1->GetNbinsX() == h2->GetNbinsX());

	double chi2 = 0;

	TH1D *hist1 = (TH1D *)h1->Clone("hist1");
	TH1D *hist2 = (TH1D *)h2->Clone("hist2");

	// cout << "Of course they are the same" << endl;

	hist1->Add(hist2, -1);
	hist1->Multiply(hist1);
	hist1->Divide(hist2);



	for (int bin = 1; bin <= hist1->GetNbinsX(); bin++){
		chi2 += hist1->GetBinContent(bin);
	}
	// cout << " Chi2 = " << chi2 << endl;
	return chi2; 
}

double Chi2Func2D(TH2D *h1, TH2D *h2){ //h1 is experimental, h2 is truth
	assert (h1->GetNbinsX() == h2->GetNbinsX());
	assert (h1->GetNbinsY() == h2->GetNbinsY());

	double chi2 = 0;

	TH2D *hist1 = (TH2D *)h1->Clone("hist1");
	TH2D *hist2 = (TH2D *)h2->Clone("hist2");

	// cout << "Of course they are the same" << endl;

	hist1->Add(hist2, -1);
	hist1->Multiply(hist1);
	hist1->Divide(hist2);



	for (int binx = 1; binx <= hist1->GetNbinsX(); binx++){
		for (int biny = 1; biny <= hist1->GetNbinsY(); biny++){
			chi2 += hist1->GetBinContent(binx, biny);
		}
	}
	// cout << " Chi2 = " << chi2 << endl;
	return chi2; 
}



void TrimHistogram(TH1 *h1, TH1 *h2, double xlow, double xhigh){

	for (int bin = 1; bin <= h2->GetNbinsX(); bin++){
		double center = h2->GetXaxis()->GetBinCenter(bin);
		if (center > xlow && center < xhigh) {
			int newbin = h1->GetXaxis()->FindBin(center);
			h1->SetBinContent(newbin, h2->GetBinContent(bin));
			h1->SetBinError(newbin, h2->GetBinError(bin));
		}
	}
	h1->SetMarkerStyle(h2->GetMarkerStyle());
	h1->SetMarkerColor(h2->GetMarkerColor());
	h1->SetLineColor(h2->GetLineColor());
	h1->SetMarkerSize(h2->GetMarkerSize());
	h1->SetLineWidth(h2->GetLineWidth());
	h1->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
	h1->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
	h1->GetXaxis()->SetRangeUser(xlow+0.001, xhigh-0.001);

	cout << h1->GetBinContent(0) << "\t" << h1->GetBinContent(h1->GetNbinsX()+1) << endl;
}

void PlotTheUnfoldedDistributionForMC(TString DirName = "tmp4", int D0pTLow = 5, int SUPERITERATION = 20, int iteration = 4, bool mimicDataUncertainty = kFALSE, bool isAreaBased = kFALSE){

	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	const int nSITER = SUPERITERATION;

	int colors[nSITER+1];
	Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,nSITER+1);
	colors[0] = kBlack;
	for (int i=1; i<=nSITER; i++) colors[i] = FI+i;

	TH2D *Truth[3];
	TH2D *Truth2DdR[3];
	TH2D *Measured[3];
	TH2D *Measured2DdR[3];

	TH2D *Unfolded[3][nSITER+1];
	TH2D *Unfolded2DdR[3][nSITER+1];

	TH1D *Ratio1D[3][nSITER+1];

	for (int i = 0; i <= nSITER; i++){
		TFile *f;
		if (isAreaBased){
			if (mimicDataUncertainty) f = new TFile(Form("%s/Output_AreaBased_DataUnc_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
			else f = new TFile(Form("%s/Output_AreaBased_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
		}
		else{
			if (mimicDataUncertainty) f = new TFile(Form("%s/Output_CS1_DataUnc_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
			else f = new TFile(Form("%s/Output_CS1_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
		}

		cout << f->GetName() << endl;
		f->cd();
		if (i == 0){
			for (int cent = 0; cent < 3; cent++){
				Truth[cent] = (TH2D *)f->Get(Form("Truth_%i", cent));
				Truth2DdR[cent] = (TH2D *)f->Get(Form("Truth2DdR_%i", cent));
				Measured[cent] = (TH2D *)f->Get(Form("Measured_%i", cent));
				Measured2DdR[cent] = (TH2D *)f->Get(Form("Measured2DdR_%i", cent));
				Truth[cent]->SetDirectory(0);
				Truth2DdR[cent]->SetDirectory(0);
				Measured[cent]->SetDirectory(0);
				Measured2DdR[cent]->SetDirectory(0);

				Ratio1D[cent][i] = (TH1D *)f->Get(Form("Non Closure 1D Cent %i SI %i", cent, i));
				Ratio1D[cent][i]->SetDirectory(0);

				// Truth[cent]->GetXaxis()->SetRangeUser(D0pTLow, 30);

				cout << Truth[cent]->GetName() << endl;
				cout << Truth2DdR[cent]->GetName() << "\t" << Truth2DdR[cent]->Integral() << endl;
				cout << Measured[cent]->GetName() << endl;
				cout << Measured2DdR[cent]->GetName() << endl;
			}
		}
		for (int cent = 0; cent < 3; cent++){
			Unfolded[cent][i] = (TH2D *)f->Get(Form("Unfolded_%i", cent));
			SetName(Unfolded[cent][i], Form("Unfolded Cent = %i SI = %i Iter = %i", cent, i, iteration));
			Unfolded[cent][i]->SetDirectory(0);

			Unfolded2DdR[cent][i] = (TH2D *)f->Get(Form("Unfolded2DdR_%i", cent));
			SetName(Unfolded2DdR[cent][i], Form("Unfolded2DdR Cent = %i SI = %i Iter = %i", cent, i, iteration));
			Unfolded2DdR[cent][i]->SetDirectory(0);

			// Unfolded[cent][i]->GetXaxis()->SetRangeUser(D0pTLow, 30);

			// cout << "SI = " << i << "\t" << Unfolded[cent][i]->GetName() << endl;
		}

		f->Close();
	}

	cout << "Imported all histograms" << endl;

	TH1D *TruthPt[3];
	TH1D *MeasuredPt[3];
	TH1D *UnfoldedPt[3][nSITER+1];
	TH1D *Chi2Pt[3];

	TH1D *TruthZ[3];
	TH1D *MeasuredZ[3];
	TH1D *UnfoldedZ[3][nSITER+1];
	TH1D *Chi2Z[3];

	TH1D *TruthdR[3];
	TH1D *MeasureddR[3];
	TH1D *UnfoldeddR[3][nSITER+1];
	TH1D *Chi2dR[3];

	int lowptbin = Unfolded[0][0]->GetXaxis()->FindBin(5);
	int highptbin = Unfolded[0][0]->GetXaxis()->FindBin(20);

	TH1D *CombinedChi2[3];

	for (int cent = 0; cent < 3; cent++){
		Chi2Pt[cent] = new TH1D(Form("Chi2Pt_%i", cent), Form("Chi2Pt_%i", cent), nSITER+1, -0.5, nSITER+0.5);
		Chi2Z[cent] = new TH1D(Form("Chi2Z_%i", cent), Form("Chi2Z_%i", cent), nSITER+1, -0.5, nSITER+0.5);
		Chi2dR[cent] = new TH1D(Form("Chi2dR_%i", cent), Form("Chi2dR_%i", cent), nSITER+1, -0.5, nSITER+0.5);
		CombinedChi2[cent] = new TH1D(Form("CombinedChi2_%i", cent), Form("CombinedChi2_%i", cent), nSITER+1, -0.5, nSITER+0.5);
	}

	cout << "Declared all 1D histograms" << endl;

	for (int cent = 0; cent < 3; cent++){
		TruthPt[cent] = (TH1D *)Truth[cent]->ProjectionX();
		TruthPt[cent]->GetXaxis()->SetRangeUser(5, 20);
		SetName(TruthPt[cent], Form("Truth p_{T} %i", cent));
		MeasuredPt[cent] = (TH1D *)Measured[cent]->ProjectionX();
		SetName(MeasuredPt[cent], Form("Measured p_{T} %i", cent));	
		
		TruthZ[cent] = (TH1D *)Truth[cent]->ProjectionY("", lowptbin, highptbin);
		SetName(TruthZ[cent], Form("Truth Z %i", cent));
		MeasuredZ[cent] = (TH1D *)Measured[cent]->ProjectionY();
		SetName(MeasuredZ[cent], Form("Measured Z %i", cent));

		cout << Truth2DdR[cent]->GetName() << "\t" << Truth2DdR[cent]->Integral() << endl;
		cout << Measured2DdR[cent]->GetName() << "\t" << Measured2DdR[cent]->Integral() << endl;

		TruthdR[cent] = (TH1D *)Truth2DdR[cent]->ProjectionY("", lowptbin, highptbin);
		SetName(TruthdR[cent], Form("Truth dR %i", cent));
		MeasureddR[cent] = (TH1D *)Measured2DdR[cent]->ProjectionY();
		SetName(MeasureddR[cent], Form("Measured dR %i", cent));

		cout << TruthdR[cent]->GetName() << "\t" << TruthdR[cent]->Integral() << endl;
		cout << MeasureddR[cent]->GetName() << "\t" << MeasureddR[cent]->Integral() << endl;

		for (int i = 0; i <= nSITER; i++){
			CombinedChi2[cent]->SetBinContent(i+1, Chi2Func2D(Unfolded[cent][i], Truth[cent]));

			UnfoldedPt[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionX();
			UnfoldedPt[cent][i]->GetXaxis()->SetRangeUser(5, 20);
			SetName(UnfoldedPt[cent][i], Form("Unfolded p_{T} Cent = %i SI = %i Iter = %i", cent, i+1, iteration));
			Chi2Pt[cent]->SetBinContent(i+1, Chi2Func(UnfoldedPt[cent][i], TruthPt[cent]));

			UnfoldedZ[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionY("", lowptbin, highptbin);
			// cout << "Cent = " << cent << "\t" << "SI = " << i << "\t" << Unfolded[cent][i]->GetNbinsX() << "\t" << Unfolded[cent][i]->GetNbinsY() << endl;
			SetName(UnfoldedZ[cent][i], Form("Unfolded Z Cent = %i SI = %i Iter = %i", cent, i+1, iteration));
			Chi2Z[cent]->SetBinContent(i+1, Chi2Func(UnfoldedZ[cent][i], TruthZ[cent]));

			UnfoldeddR[cent][i] = (TH1D *)Unfolded2DdR[cent][i]->ProjectionY("", lowptbin, highptbin);
			SetName(UnfoldeddR[cent][i], Form("Unfolded dR Cent = %i SI = %i Iter = %i", cent, i+1, iteration));
			Chi2dR[cent]->SetBinContent(i+1, Chi2Func(UnfoldeddR[cent][i], TruthdR[cent]));
		}
	}

	cout << "Made 1D Histograms" << endl;

	TCanvas *Pt = new TCanvas("pT", "pT", 2000, 1000);
	Pt->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Pt->cd(cent+1);
		gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			// UnfoldedPt[cent][i]->Rebin(5);
			SetColor(UnfoldedPt[cent][i], colors[i], 20, 1.5);
			UnfoldedPt[cent][i]->Draw("EP SAME");
			UnfoldedPt[cent][i]->GetYaxis()->SetRangeUser(pow(10,-2), pow(10, 8));
			ProcessSpectra(UnfoldedPt[cent][i]);
			UnfoldedPt[cent][i]->GetYaxis()->SetTitle("dN_{Jet}/dp_{T} (GeV/c)^{-1}");
			UnfoldedPt[cent][i]->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
		}

		// TruthPt[cent]->Rebin(5);
		ProcessSpectra(TruthPt[cent]);
		SetColor(TruthPt[cent], kGreen-2, 29, 1.5);
		TruthPt[cent]->Draw("EP SAME");

		gPad->BuildLegend();
	}

	if (isAreaBased){
		if (mimicDataUncertainty) Pt->SaveAs(Form("%s/Plots_AreaBased_DataUnc/UnfoldedpT.pdf", DirName.Data()));
		else Pt->SaveAs(Form("%s/Plots_AreaBased/UnfoldedpT.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) Pt->SaveAs(Form("%s/Plots_CS1_DataUnc/UnfoldedpT.pdf", DirName.Data()));
		else Pt->SaveAs(Form("%s/Plots_CS1/UnfoldedpT.pdf", DirName.Data()));
	}

	TCanvas *Z = new TCanvas("Z", "Z", 2000, 1000);
	Z->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Z->cd(cent+1);
		gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			// UnfoldedZ[cent][i]->Rebin(5);
			SetColor(UnfoldedZ[cent][i], colors[i], 20, 1.5);
			UnfoldedZ[cent][i]->Draw("EP SAME");
			UnfoldedZ[cent][i]->GetYaxis()->SetRangeUser(pow(10,0), pow(10, 8));
			ProcessSpectra(UnfoldedZ[cent][i]);
			UnfoldedZ[cent][i]->GetYaxis()->SetTitle("dN_{Jet}/dz");
			UnfoldedZ[cent][i]->GetXaxis()->SetTitle("z");
			
		}
		ProcessSpectra(TruthZ[cent]);
		SetColor(TruthZ[cent], kGreen-2, 29, 1.5);
		TruthZ[cent]->Draw("EP SAME");

		gPad->BuildLegend();
	}

	// Z->SaveAs(Form("%s/Plots/UnfoldedZ.pdf", DirName.Data()));

	if (isAreaBased){
		if (mimicDataUncertainty) Z->SaveAs(Form("%s/Plots_AreaBased_DataUnc/UnfoldedZ.pdf", DirName.Data()));
		else Z->SaveAs(Form("%s/Plots_AreaBased/UnfoldedZ.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) Z->SaveAs(Form("%s/Plots_CS1_DataUnc/UnfoldedZ.pdf", DirName.Data()));
		else Z->SaveAs(Form("%s/Plots_CS1/UnfoldedZ.pdf", DirName.Data()));
	}

	TCanvas *dR = new TCanvas("dR", "dR", 2000, 1000);
	dR->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		dR->cd(cent+1);
		gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			// UnfoldedZ[cent][i]->Rebin(5);
			SetColor(UnfoldeddR[cent][i], colors[i], 20, 1.5);
			UnfoldeddR[cent][i]->Draw("EP SAME");
			UnfoldeddR[cent][i]->GetYaxis()->SetRangeUser(pow(10,0), pow(10, 8));
			UnfoldeddR[cent][i]->GetXaxis()->SetRangeUser(0, 0.2);
			DivideByBinWidth(UnfoldeddR[cent][i]);
			UnfoldeddR[cent][i]->GetYaxis()->SetTitle("dN_{Jet}/d#DeltaR");
			UnfoldeddR[cent][i]->GetXaxis()->SetTitle("#DeltaR");
			
		}
		DivideByBinWidth(TruthdR[cent]);
		SetColor(TruthdR[cent], kGreen-2, 29, 1.5);
		TruthdR[cent]->Draw("EP SAME");

		gPad->BuildLegend();
	}

	// Z->SaveAs(Form("%s/Plots/UnfoldedZ.pdf", DirName.Data()));

	if (isAreaBased){
		if (mimicDataUncertainty) dR->SaveAs(Form("%s/Plots_AreaBased_DataUnc/UnfoldeddR.pdf", DirName.Data()));
		else dR->SaveAs(Form("%s/Plots_AreaBased/UnfoldeddR.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) dR->SaveAs(Form("%s/Plots_CS1_DataUnc/UnfoldeddR.pdf", DirName.Data()));
		else dR->SaveAs(Form("%s/Plots_CS1/UnfoldeddR.pdf", DirName.Data()));
	}

	TH2D *Ratio[3][nSITER+1];
	for (int cent = 0; cent < 3; cent++){
		for (int i = 0; i <= nSITER; i++){
			Ratio[cent][i] = (TH2D *)Unfolded[cent][i]->Clone();
			SetName(Ratio[cent][i], Form("Non Closure Cent %i", cent));
			cout << "Ratio = " << Ratio[cent][i]->Integral() << endl;
			Ratio[cent][i]->Divide(Truth[cent]);
			
			cout << Ratio[cent][i]->Integral() << endl;
		}
	}

	TH1D *RatioPt[3][nSITER+1];
	TH1D *RatioZ[3][nSITER+1];
	TH1D *RatiodR[3][nSITER+1];

	TCanvas *RatioCanvasPt = new TCanvas("Ratio Canvas pT", "Ratio Canvas pT", 2000, 1000);
	RatioCanvasPt->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		RatioCanvasPt->cd(cent+1);
		// gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			RatioPt[cent][i] = (TH1D *)UnfoldedPt[cent][i]->Clone();
			RatioPt[cent][i]->Divide(TruthPt[cent]);
			RatioPt[cent][i]->GetYaxis()->SetTitle("R_{CP}");

			SetColor(RatioPt[cent][i], colors[i], 20);
			RatioPt[cent][i]->Draw("EP SAME");
			RatioPt[cent][i]->GetYaxis()->SetRangeUser(0, 2);
		}

		// gPad->BuildLegend();
	
		// TH1D *tmp = (TH1D *)GetLineAtOne(TruthPt[cent]);
		// tmp->Draw("P HIST SAME");

		TLine *tmp = (TLine *)GetLineAtOne(TruthPt[cent]);
		tmp->Draw("SAME");
	}

	// RatioCanvasPt->SaveAs(Form("%s/Plots/RatiopT.pdf", DirName.Data()));

	if (isAreaBased){
		if (mimicDataUncertainty) RatioCanvasPt->SaveAs(Form("%s/Plots_AreaBased_DataUnc/RatiopT.pdf", DirName.Data()));
		else RatioCanvasPt->SaveAs(Form("%s/Plots_AreaBased/RatiopT.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) RatioCanvasPt->SaveAs(Form("%s/Plots_CS1_DataUnc/RatiopT.pdf", DirName.Data()));
		else RatioCanvasPt->SaveAs(Form("%s/Plots_CS1/RatiopT.pdf", DirName.Data()));
	}


	TCanvas *RatioCanvasZ = new TCanvas("Ratio Canvas Z", "Ratio Canvas Z", 2000, 1000);
	RatioCanvasZ->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		RatioCanvasZ->cd(cent+1);
		// gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			// cout << "Cent = " << cent << " SI = " << i << " Bins = " << UnfoldedZ[cent][i]->GetNbinsX() << "\t" << TruthZ[cent]->GetNbinsX() << endl;
			RatioZ[cent][i] = (TH1D *)UnfoldedZ[cent][i]->Clone();
			RatioZ[cent][i]->Divide(TruthZ[cent]);
			RatioZ[cent][i]->GetYaxis()->SetTitle("R_{CP}");

			SetColor(RatioZ[cent][i], colors[i], 20);
			RatioZ[cent][i]->Draw("EP SAME");
			RatioZ[cent][i]->GetYaxis()->SetRangeUser(0, 2);
		}

		// gPad->BuildLegend();

		TLine *tmp = (TLine *)GetLineAtOne(TruthZ[cent]);
		tmp->Draw("SAME");
	}

	// RatioCanvasZ->SaveAs(Form("%s/Plots/RatioZ.pdf", DirName.Data()));
	if (isAreaBased){
		if (mimicDataUncertainty) RatioCanvasZ->SaveAs(Form("%s/Plots_AreaBased_DataUnc/RatioZ.pdf", DirName.Data()));
		else RatioCanvasZ->SaveAs(Form("%s/Plots_AreaBased/RatioZ.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) RatioCanvasZ->SaveAs(Form("%s/Plots_CS1_DataUnc/RatioZ.pdf", DirName.Data()));
		else RatioCanvasZ->SaveAs(Form("%s/Plots_CS1/RatioZ.pdf", DirName.Data()));
	}

	///////////////////////////////////////////
	TCanvas *RatioCanvasdR = new TCanvas("Ratio Canvas dR", "Ratio Canvas dR", 2000, 1000);
	RatioCanvasdR->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		RatioCanvasdR->cd(cent+1);
		// gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			// cout << "Cent = " << cent << " SI = " << i << " Bins = " << UnfoldedZ[cent][i]->GetNbinsX() << "\t" << TruthZ[cent]->GetNbinsX() << endl;
			RatiodR[cent][i] = (TH1D *)UnfoldeddR[cent][i]->Clone();
			RatiodR[cent][i]->Divide(TruthdR[cent]);
			RatiodR[cent][i]->GetYaxis()->SetTitle("R_{CP}");

			SetColor(RatiodR[cent][i], colors[i], 20);
			RatiodR[cent][i]->Draw("EP SAME");
			RatiodR[cent][i]->GetYaxis()->SetRangeUser(0, 2);
		}

		// gPad->BuildLegend();

		TLine *tmp = (TLine *)GetLineAtOne(TruthdR[cent]);
		tmp->Draw("SAME");
	}

	// RatioCanvasZ->SaveAs(Form("%s/Plots/RatioZ.pdf", DirName.Data()));
	if (isAreaBased){
		if (mimicDataUncertainty) RatioCanvasdR->SaveAs(Form("%s/Plots_AreaBased_DataUnc/RatiodR.pdf", DirName.Data()));
		else RatioCanvasdR->SaveAs(Form("%s/Plots_AreaBased/RatiodR.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) RatioCanvasdR->SaveAs(Form("%s/Plots_CS1_DataUnc/RatiodR.pdf", DirName.Data()));
		else RatioCanvasdR->SaveAs(Form("%s/Plots_CS1/RatiodR.pdf", DirName.Data()));
	}

	TH1D *UnfoldedPtTrimmed[3][nSITER+1];
	TH1D *UnfoldedZTrimmed[3][nSITER+1];
	TH1D *UnfoldeddRTrimmed[3][nSITER+1];

	TH1D *TruthPtTrimmed[3];
	TH1D *TruthZTrimmed[3];
	TH1D *TruthdRTrimmed[3];

	double plottingpT[7] = {5, 7, 9, 11, 13, 15, 20};
	double plottingdR[4] = {0, 0.05, 0.1, 0.2};

	for (int i = 0; i < 3; i++){
		UnfoldedPtTrimmed[i][0] = new TH1D (Form("%s_Trimmed", UnfoldedPt[i][0]->GetName()), Form("%s_Trimmed", UnfoldedPt[i][0]->GetName()), 6, plottingpT);
		UnfoldeddRTrimmed[i][0] = new TH1D (Form("%s_Trimmed", UnfoldeddR[i][0]->GetName()), Form("%s_Trimmed", UnfoldeddR[i][0]->GetName()), 3, plottingdR);

		TruthPtTrimmed[i] = new TH1D (Form("%s_Trimmed", TruthPt[i]->GetName()), Form("%s_Trimmed", TruthPt[i]->GetName()), 6, plottingpT);
		TruthdRTrimmed[i] = new TH1D (Form("%s_Trimmed", TruthdR[i]->GetName()), Form("%s_Trimmed", TruthdR[i]->GetName()), 3, plottingdR);

		TrimHistogram(UnfoldedPtTrimmed[i][0], UnfoldedPt[i][0], 5, 20);
		TrimHistogram(UnfoldeddRTrimmed[i][0], UnfoldeddR[i][0], 0, 0.2);

		TrimHistogram(TruthPtTrimmed[i], TruthPt[i], 5, 20);
		TrimHistogram(TruthdRTrimmed[i], TruthdR[i], 0, 0.2);
	}


	TRatioPlot *rpPt[3][nSITER+1];
	TRatioPlot *rpZ[3][nSITER+1];
	TRatioPlot *rpdR[3][nSITER+1];

	TCanvas *RatioPlotCanvasPt = new TCanvas("Ratio Plot Canvas pT", "Ratio Plot Canvas pT", 2000, 1000);
	RatioPlotCanvasPt->Divide(3);
	for (int i = 0; i < 3; i++){
		RatioPlotCanvasPt->cd(i+1);
		UnfoldedPt[i][0]->GetXaxis()->SetRange(1,6);
		TruthPtTrimmed[i]->GetXaxis()->SetRangeUser(1,3);
		rpPt[i][0] = new TRatioPlot(UnfoldedPtTrimmed[i][0], TruthPtTrimmed[i]);
		// gPad->SetTicks(0, 1);
		rpPt[i][0]->Draw();
		rpPt[i][0]->GetLowerRefGraph()->SetMinimum(0);
		rpPt[i][0]->GetLowerRefGraph()->SetMaximum(2);
		rpPt[i][0]->GetLowYaxis()->SetNdivisions(505);
		rpPt[i][0]->GetUpperRefYaxis()->SetTitle("Unfolded/Truth");
		rpPt[i][0]->GetUpperRefYaxis()->SetTitleOffset(1.5);
		rpPt[i][0]->GetUpperRefXaxis()->SetRangeUser(0, 20);
		rpPt[i][0]->RangeAxisChanged();
		rpPt[i][0]->GetUpperPad()->SetLogy();
		rpPt[i][0]->GetUpperPad()->cd();
		gPad->BuildLegend();
		RatioPlotCanvasPt->cd(i+1);
		gPad->Modified();
		gPad->Update();
		

	}
	if (isAreaBased){
		if (mimicDataUncertainty) RatioPlotCanvasPt->SaveAs(Form("%s/Plots_AreaBased_DataUnc/RatioPlotPt.pdf", DirName.Data()));
		else RatioPlotCanvasPt->SaveAs(Form("%s/Plots_AreaBased/RatioPlotPt.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) RatioPlotCanvasPt->SaveAs(Form("%s/Plots_CS1_DataUnc/RatioPlotPt.pdf", DirName.Data()));
		else RatioPlotCanvasPt->SaveAs(Form("%s/Plots_CS1/RatioPlotPt.pdf", DirName.Data()));
	}
	
	TCanvas *RatioPlotCanvasZ = new TCanvas("Ratio Plot Canvas Z", "Ratio Plot Canvas Z", 2000, 1000);
	RatioPlotCanvasZ->Divide(3);
	for (int i = 0; i < 3; i++){
		RatioPlotCanvasZ->cd(i+1);
		rpZ[i][0] = new TRatioPlot(UnfoldedZ[i][0], TruthZ[i]);
		gPad->SetTicks(0, 1);
		rpZ[i][0]->Draw();
		rpZ[i][0]->GetLowerRefGraph()->SetMinimum(0);
		rpZ[i][0]->GetLowerRefGraph()->SetMaximum(2);
		rpZ[i][0]->GetUpperRefXaxis()->SetRangeUser(0, 1);
		rpZ[i][0]->GetLowYaxis()->SetNdivisions(505);
		rpZ[i][0]->GetUpperRefYaxis()->SetTitle("Unfolded/Truth");
		rpZ[i][0]->GetUpperRefYaxis()->SetTitleOffset(1.5);
		rpZ[i][0]->GetUpperPad()->SetLogy();
		rpZ[i][0]->GetUpperPad()->cd();
		gPad->BuildLegend();
	}

	if (isAreaBased){
		if (mimicDataUncertainty) RatioPlotCanvasZ->SaveAs(Form("%s/Plots_AreaBased_DataUnc/RatioPlotZ.pdf", DirName.Data()));
		else RatioPlotCanvasZ->SaveAs(Form("%s/Plots_AreaBased/RatioPlotZ.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) RatioPlotCanvasZ->SaveAs(Form("%s/Plots_CS1_DataUnc/RatioPlotZ.pdf", DirName.Data()));
		else RatioPlotCanvasZ->SaveAs(Form("%s/Plots_CS1/RatioPlotZ.pdf", DirName.Data()));
	}



	TCanvas *RatioPlotCanvasdR = new TCanvas("Ratio Plot Canvas dR", "Ratio Plot Canvas dR", 2000, 1000);
	RatioPlotCanvasdR->Divide(3);
	for (int i = 0; i < 3; i++){
		RatioPlotCanvasdR->cd(i+1);
		rpdR[i][0] = new TRatioPlot(UnfoldeddRTrimmed[i][0], TruthdRTrimmed[i]);
		gPad->SetTicks(0, 1);
		rpdR[i][0]->Draw();
		rpdR[i][0]->GetLowerRefGraph()->SetMinimum(0);
		rpdR[i][0]->GetLowerRefGraph()->SetMaximum(2);
		rpdR[i][0]->GetUpperRefXaxis()->SetRangeUser(0, 0.2);
		rpdR[i][0]->GetLowYaxis()->SetNdivisions(505);
		rpdR[i][0]->GetUpperRefYaxis()->SetTitle("Unfolded/Truth");
		rpdR[i][0]->GetUpperRefYaxis()->SetTitleOffset(1.5);
		rpdR[i][0]->GetUpperPad()->SetLogy();
		rpdR[i][0]->GetUpperPad()->cd();
		gPad->BuildLegend();
	}

	if (isAreaBased){
		if (mimicDataUncertainty) RatioPlotCanvasdR->SaveAs(Form("%s/Plots_AreaBased_DataUnc/RatioPlotdR.pdf", DirName.Data()));
		else RatioPlotCanvasdR->SaveAs(Form("%s/Plots_AreaBased/RatioPlotdR.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) RatioPlotCanvasdR->SaveAs(Form("%s/Plots_CS1_DataUnc/RatioPlotdR.pdf", DirName.Data()));
		else RatioPlotCanvasdR->SaveAs(Form("%s/Plots_CS1/RatioPlotdR.pdf", DirName.Data()));
	}

	///////////////////////////////////////////

	TCanvas *Chi2CanvasPt = new TCanvas("Chi2 Canvas pT", "Chi2 Canvas pT", 2000, 1000);
	Chi2CanvasPt->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Chi2CanvasPt->cd(cent+1);
		Chi2Pt[cent]->Draw("HIST SAME");
		// Chi2Pt[cent]->GetYaxis()->SetRangeUser(0, 100);
		gPad->SetLogy();
	}

	// Chi2CanvasPt->SaveAs(Form("%s/Plots/Chi2pT.pdf", DirName.Data()));
	if (isAreaBased){
		if (mimicDataUncertainty) Chi2CanvasPt->SaveAs(Form("%s/Plots_AreaBased_DataUnc/Chi2pT.pdf", DirName.Data()));
		else Chi2CanvasPt->SaveAs(Form("%s/Plots_AreaBased/Chi2pT.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) Chi2CanvasPt->SaveAs(Form("%s/Plots_CS1_DataUnc/Chi2pT.pdf", DirName.Data()));
		else Chi2CanvasPt->SaveAs(Form("%s/Plots_CS1/Chi2pT.pdf", DirName.Data()));
	}

	TCanvas *Chi2CanvasZ = new TCanvas("Chi2 Canvas Z", "Chi2 Canvas Z", 2000, 1000);
	Chi2CanvasZ->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Chi2CanvasZ->cd(cent+1);
		Chi2Z[cent]->Draw("HIST SAME");
		// Chi2Z[cent]->GetYaxis()->SetRangeUser(0, 100);
		gPad->SetLogy();
	}

	// Chi2CanvasZ->SaveAs(Form("%s/Plots/Chi2Z.pdf", DirName.Data()));
	if (isAreaBased){
		if (mimicDataUncertainty) Chi2CanvasZ->SaveAs(Form("%s/Plots_AreaBased_DataUnc/Chi2Z.pdf", DirName.Data()));
		else Chi2CanvasZ->SaveAs(Form("%s/Plots_AreaBased/Chi2Z.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) Chi2CanvasZ->SaveAs(Form("%s/Plots_CS1_DataUnc/Chi2Z.pdf", DirName.Data()));
		else Chi2CanvasZ->SaveAs(Form("%s/Plots_CS1/Chi2Z.pdf", DirName.Data()));
	}

	TCanvas *Chi2Canvas = new TCanvas("Chi2 Canvas", "Chi2 Canvas", 2000, 1000);
	Chi2Canvas->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Chi2Canvas->cd(cent+1);
		CombinedChi2[cent]->Draw("HIST SAME");
		// Chi2Z[cent]->GetYaxis()->SetRangeUser(0, 100);
		gPad->SetLogy();

		cout << CombinedChi2[cent]->GetBinCenter(CombinedChi2[cent]->GetMinimumBin()) << endl;
	}

	// Chi2Canvas->SaveAs(Form("%s/Plots/CombinedChi2.pdf", DirName.Data()));
	if (isAreaBased){
		if (mimicDataUncertainty) Chi2Canvas->SaveAs(Form("%s/Plots_AreaBased_DataUnc/CombinedChi2.pdf", DirName.Data()));
		else Chi2Canvas->SaveAs(Form("%s/Plots_AreaBased/CombinedChi2.pdf", DirName.Data()));
	}
	else{
		if (mimicDataUncertainty) Chi2Canvas->SaveAs(Form("%s/Plots_CS1_DataUnc/CombinedChi2.pdf", DirName.Data()));
		else Chi2Canvas->SaveAs(Form("%s/Plots_CS1/CombinedChi2.pdf", DirName.Data()));
	}


	TFile *SummaryStat;
	if (isAreaBased){
		if (mimicDataUncertainty) SummaryStat = new TFile(Form("%s/Plots_AreaBased_DataUnc/Summary.root", DirName.Data()), "RECREATE");
		else SummaryStat = new TFile(Form("%s/Plots_AreaBased/Summary.root", DirName.Data()), "RECREATE");
	}
	else{
		if (mimicDataUncertainty) SummaryStat = new TFile(Form("%s/Plots_CS1_DataUnc/Summary.root", DirName.Data()), "RECREATE");
		else SummaryStat = new TFile(Form("%s/Plots_CS1/Summary.root", DirName.Data()), "RECREATE");
	}
	SummaryStat->cd();

	for (int cent = 0; cent < 3; cent++){
		TruthPt[cent]->Write(Form("Truth Pt Cent %i", cent));
		TruthZ[cent]->Write(Form("Truth Z Cent %i", cent));
		TruthdR[cent]->Write(Form("Truth dR Cent %i", cent));
		Truth[cent]->Write(Form("Truth Cent %i", cent));

		Chi2Pt[cent]->Write(Form("Chi2 Pt Cent %i", cent));
		Chi2Z[cent]->Write(Form("Chi2 Z Cent %i", cent));
		CombinedChi2[cent]->Write(Form("Combined Chi2 Cent %i", cent));

		Ratio[cent][0]->Write(Form("Non Closure Cent %i", cent));
		cout << Ratio[cent][0]->GetName() << endl;
	}
	for (int cent = 0; cent < 3; cent++){
		for (int i = 0; i <= nSITER; i++){
			UnfoldedPt[cent][i]->Write(Form("Unfolded Pt Cent %i SI %i", cent, i));
			RatioPt[cent][i]->Write(Form("Non Closure Pt Cent %i SI %i", cent, i));
			UnfoldedZ[cent][i]->Write(Form("Unfolded Z Cent %i SI %i", cent, i));
			RatioZ[cent][i]->Write(Form("Non Closure Z Cent %i SI %i", cent, i));
			UnfoldeddR[cent][i]->Write(Form("Unfolded dR Cent %i SI %i", cent, i));
			RatiodR[cent][i]->Write(Form("Non Closure dR Cent %i SI %i", cent, i));
			Unfolded[cent][i]->Write(Form("Unfolded Cent %i SI %i", cent, i));
			Unfolded2DdR[cent][i]->Write(Form("Unfolded2DdR Cent %i SI %i", cent, i));
			
		}
	}
	for (int cent = 0; cent < 3; cent++){
		Ratio1D[cent][0]->Write(Form("Non Closure 1D Cent %i SI %i", cent, 0));
	} 
	cout << "Wrote to file " << SummaryStat->GetName() << endl;
	SummaryStat->Close();
}