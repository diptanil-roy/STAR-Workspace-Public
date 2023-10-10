using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void PlotTheUnfoldedDistributionForAreaData(TString DirName = "tmp4", int D0pTLow = 5, int SUPERITERATION = 20, int iteration = 4, TString NonClosureDir = "", bool isSimpleUnfolding = false){

	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);

  	// Scaling by TAA

	double taa[3] = {941.23714*1.0318440e+08, 391.35550*3.2123506e+08, 56.62475*4.6679240e+08};

	const int nSITER = SUPERITERATION;

	TH2D *Truth[3];
	TH2D *Measured[3];
	TH2D *MeasureddR[3];
	TH1D *Unfolded1D[3][nSITER+1];
	TH2D *Unfolded[3][nSITER+1];
	TH2D *UnfoldeddR[3][nSITER+1];

	TH1D *UnfoldedpTCS[3][nSITER+1];
	TH1D *UnfoldedZCS[3][nSITER+1];
	TH1D *RCPPtCS[3][nSITER+1];
	TH1D *RCPZCS[3][nSITER+1];

	for (int i = 0; i <= nSITER; i++){
		TFile *f;
		
		if (isSimpleUnfolding) f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
		else f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, iteration, njpt_gen_bins_var, nz_gen_bins), "READ");
		
		cout << f->GetName() << endl;
		f->cd();
		if (i == 0){
			for (int cent = 0; cent < 3; cent++){
				Measured[cent] = (TH2D *)f->Get(Form("MeasuredWide_%i", cent));
				Measured[cent]->SetDirectory(0);
				cout << Measured[cent]->GetName() << endl;

				MeasureddR[cent] = (TH2D *)f->Get(Form("MeasureddRWide_%i", cent));
				MeasureddR[cent]->SetDirectory(0);
				cout << MeasureddR[cent]->GetName() << endl;
			}
		}
		for (int cent = 0; cent < 3; cent++){
			Unfolded1D[cent][i] = (TH1D *)f->Get(Form("Unfolded1D_%i", cent));
			SetName(Unfolded1D[cent][i], Form("Unfolded1D Cent = %i SI = %i Iter = %i", cent, i, iteration));
			Unfolded1D[cent][i]->SetDirectory(0);
			cout << Unfolded1D[cent][i]->GetName() << endl;

			Unfolded[cent][i] = (TH2D *)f->Get(Form("UnfoldedWide_%i", cent));
			SetName(Unfolded[cent][i], Form("UnfoldedWide Cent = %i SI = %i Iter = %i", cent, i, iteration));
			Unfolded[cent][i]->SetDirectory(0);
			cout << Unfolded[cent][i]->GetName() << endl;

			UnfoldeddR[cent][i] = (TH2D *)f->Get(Form("Unfolded2DdRWide_%i", cent));
			SetName(UnfoldeddR[cent][i], Form("UnfoldeddRWide Cent = %i SI = %i Iter = %i", cent, i, iteration));
			UnfoldeddR[cent][i]->SetDirectory(0);
			cout << UnfoldeddR[cent][i]->GetName() << endl;

			// Unfolded[cent][i]->GetXaxis()->SetRangeUser(D0pTLow, 30);

			UnfoldedpTCS[cent][i] = (TH1D *)f->Get(Form("Unfolded p_{T} Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
			UnfoldedpTCS[cent][i]->SetDirectory(0);
			UnfoldedZCS[cent][i] = (TH1D *)f->Get(Form("Unfolded Z Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
			UnfoldedZCS[cent][i]->SetDirectory(0);

			cout << "SI = " << i << "\t" << Unfolded[cent][i]->GetName() << endl;
		}

		RCPPtCS[0][i] = (TH1D *)f->Get(Form("R_{CP} p_{T} #SI %i", i));
		RCPPtCS[1][i] = (TH1D *)f->Get(Form("R_{MP} p_{T} #SI %i", i));

		RCPPtCS[0][i]->SetDirectory(0);
		RCPPtCS[1][i]->SetDirectory(0);


		RCPZCS[0][i] = (TH1D *)f->Get(Form("R_{CP} Z #SI %i", i));
		RCPZCS[1][i] = (TH1D *)f->Get(Form("R_{MP} Z #SI %i", i));

		RCPZCS[0][i]->SetDirectory(0);
		RCPZCS[1][i]->SetDirectory(0);

		f->Close();
	}

	cout << "Imported all histograms" << endl;

	int lowptbin = Unfolded[0][0]->GetXaxis()->FindBin(5);
	int highptbin = Unfolded[0][0]->GetXaxis()->FindBin(20);

	TH1D *MeasuredPt[3];
	TH1D *UnfoldedPt[3][nSITER+1];

	TH1D *MeasuredZ[3];
	TH1D *UnfoldedZ[3][nSITER+1];
	TH1D *Unfolded1DdR[3][nSITER+1];

	TH1D *Unfolded1DRebinned[3][nSITER+1];
	TH1D *UnfoldedPtRebinned[3][nSITER+1];
	TH1D *UnfoldedZRebinned[3][nSITER+1];
	TH1D *UnfoldeddRRebinned[3][nSITER+1];

	for (int cent = 0; cent < 3; cent++){
		MeasuredPt[cent] = (TH1D *)Measured[cent]->ProjectionX();
		SetName(MeasuredPt[cent], Form("Measured Wide p_{T} %i", cent));	
		
		MeasuredZ[cent] = (TH1D *)Measured[cent]->ProjectionY();
		SetName(MeasuredZ[cent], Form("Measured Wide Z %i", cent));

		for (int i = 0; i <= nSITER; i++){
			Unfolded1D[cent][i]->GetXaxis()->SetRangeUser(5, 20);

			UnfoldedPt[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionX();
			UnfoldedPt[cent][i]->GetXaxis()->SetRangeUser(5, 20);
			SetName(UnfoldedPt[cent][i], Form("Unfolded Wide p_{T} Cent = %i SI = %i Iter = %i", cent, i, iteration));
			UnfoldedZ[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionY("", lowptbin, highptbin);
			SetName(UnfoldedZ[cent][i], Form("Unfolded Wide Z Cent = %i SI = %i Iter = %i", cent, i, iteration));
			Unfolded1DdR[cent][i] = (TH1D *)UnfoldeddR[cent][i]->ProjectionY("", lowptbin, highptbin);
			SetName(Unfolded1DdR[cent][i], Form("Unfolded Wide dR Cent = %i SI = %i Iter = %i", cent, i, iteration));

			Unfolded1DRebinned[cent][i] = (TH1D *)Unfolded1D[cent][i]->Clone();
			SetName(Unfolded1DRebinned[cent][i], Form("Unfolded1D Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
			UnfoldedPtRebinned[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionX();
			UnfoldedPtRebinned[cent][i]->GetXaxis()->SetRangeUser(5, 20);
			SetName(UnfoldedPtRebinned[cent][i], Form("Unfolded Wide p_{T} Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
			UnfoldedZRebinned[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionY("", lowptbin, highptbin);
			SetName(UnfoldedZRebinned[cent][i], Form("Unfolded Wide Z Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
			UnfoldeddRRebinned[cent][i] = (TH1D *)UnfoldeddR[cent][i]->ProjectionY("", lowptbin, highptbin);
			SetName(UnfoldeddRRebinned[cent][i], Form("Unfolded Wide dR Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
		}
	}

	cout << "Made 1D Histograms" << endl;

	if (NonClosureDir != ""){

		TFile *SummaryFile = new TFile(Form("%s/Summary.root", NonClosureDir.Data()), "READ");
		SummaryFile->cd();

		TH1D *tmp1DNonClosure;
		TH1D *tmpPtNonClosure;
		TH1D *tmpZNonClosure;
		TH1D *tmpdRNonClosure;

		for (int cent = 0; cent < 3; cent++){
			for (int i = 0; i <= nSITER; i++){
				tmp1DNonClosure = (TH1D *)SummaryFile->Get(Form("Non Closure 1D Cent %i SI %i", cent, i));
				tmpPtNonClosure = (TH1D *)SummaryFile->Get(Form("Non Closure Pt Cent %i SI %i", cent, i));
				tmpZNonClosure = (TH1D *)SummaryFile->Get(Form("Non Closure Z Cent %i SI %i", cent, i));
				tmpdRNonClosure = (TH1D *)SummaryFile->Get(Form("Non Closure dR Cent %i SI %i", cent, i));

				Unfolded1DRebinned[cent][i]->Divide(tmp1DNonClosure);
				UnfoldedPtRebinned[cent][i]->Divide(tmpPtNonClosure);
				UnfoldedZRebinned[cent][i]->Divide(tmpZNonClosure);
				UnfoldeddRRebinned[cent][i]->Divide(tmpdRNonClosure);
				cout << "Division Done" << endl;
			}
		}
	}

	// TH1D *UnfoldedPtRebinned[3][nSITER+1];

	TCanvas *Pt = new TCanvas("pT", "pT", 2000, 1000);
	Pt->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Pt->cd(cent+1);
		gPad->SetLogy();

		// ProcessSpectra(MeasuredPt[cent]);
		// MeasuredPt[cent]->Scale(1./taa[cent]);
		// SetColor(MeasuredPt[cent], kRed, 29, 1.5);
		// MeasuredPt[cent]->Draw("EP SAME");

		// MeasuredPt[cent]->GetYaxis()->SetRangeUser(pow(10,-11), pow(10, -3));

		for (int i = 0; i <= nSITER; i++){
			// UnfoldedPt[cent][i]->Rebin(2);
			// UnfoldedPtRebinned[cent][i] = (TH1D *)UnfoldedPt[cent][i]->Rebin(nbins_jpt, Form("%s_Rebinned", UnfoldedPt[cent][i]->GetName()), binning_jpt);
			// UnfoldedPtRebinned[cent][i] = (TH1D *)UnfoldedPt[cent][i]->Clone();
			// ProcessSpectra(UnfoldedPtRebinned[cent][i]);
			// UnfoldedPtRebinned[cent][i]->Scale(1./taa[cent]);
			
			SetColor(UnfoldedPtRebinned[cent][i], i + 1, 20);
			UnfoldedPtRebinned[cent][i]->Draw("EP SAME");
			UnfoldedPtRebinned[cent][i]->GetYaxis()->SetRangeUser(pow(10,-1), pow(10, 5));

			SetColor(Unfolded1DRebinned[cent][i], i + 1, 24);
			Unfolded1DRebinned[cent][i]->Draw("EP SAME");
			Unfolded1DRebinned[cent][i]->GetYaxis()->SetRangeUser(pow(10,-1), pow(10, 5));
		}

		// MeasuredPt[cent]->Rebin(2);

		gPad->BuildLegend();
	}

	Pt->SaveAs(Form("%s/Plots/UnfoldedpT_Wide_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TCanvas *Z = new TCanvas("Z", "Z", 2000, 1000);
	Z->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Z->cd(cent+1);
		gPad->SetLogy();

		// SetColor(MeasuredZ[cent], kRed, 29, 1.5);
		// ProcessSpectra(MeasuredZ[cent]);
		// MeasuredZ[cent]->Scale(1./taa[cent]);
		

		// MeasuredZ[cent]->Draw("EP SAME");

		// MeasuredZ[cent]->GetYaxis()->SetRangeUser(pow(10,-9), pow(10, -1));

		for (int i = 0; i <= nSITER; i++){
			// UnfoldedZ[cent][i]->Rebin(4);
			// UnfoldedZ[cent][i]->Scale(1./taa[cent]);
			SetColor(UnfoldedZRebinned[cent][i], i + 1, 20);
			// ProcessSpectra(UnfoldedZ[cent][i]);
			UnfoldedZRebinned[cent][i]->Draw("EP SAME");
			UnfoldedZRebinned[cent][i]->GetYaxis()->SetRangeUser(pow(10,-1), pow(10, 5));
			
		}

		gPad->BuildLegend();
	}

	Z->SaveAs(Form("%s/Plots/UnfoldedZ_Wide_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TCanvas *dR = new TCanvas("dR", "dR", 2000, 1000);
	dR->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		dR->cd(cent+1);
		gPad->SetLogy();

		// SetColor(MeasureddR[cent], kRed, 29, 1.5);
		// ProcessSpectra(MeasureddR[cent]);
		// MeasureddR[cent]->Scale(1./taa[cent]);
		// MeasureddR[cent]->Draw("EP SAME");

		// MeasureddR[cent]->GetYaxis()->SetRangeUser(pow(10,-11), pow(10, -3));

		for (int i = 0; i <= nSITER; i++){
			// UnfoldeddR[cent][i]->Rebin(2);
			// UnfoldeddR[cent][i]->Scale(1./taa[cent]);
			SetColor(UnfoldeddRRebinned[cent][i], i + 1, 20);
			UnfoldeddRRebinned[cent][i]->Scale(1./UnfoldeddRRebinned[cent][i]->Integral());
			// ProcessSpectra(UnfoldeddR[cent][i]);
			UnfoldeddRRebinned[cent][i]->Draw("EP SAME");
			UnfoldeddRRebinned[cent][i]->GetYaxis()->SetRangeUser(pow(10,-10), pow(10, 0));
		}

		gPad->BuildLegend();
	}

	dR->SaveAs(Form("%s/Plots/UnfoldeddR_Wide_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TH1D *RCPPt[3][nSITER+1];
	TH1D *RCPZ[3][nSITER+1];
	TH1D *RCPdR[3][nSITER+1];

	TCanvas *RatioCanvasPt = new TCanvas("RCP pT", "RCP pT", 2000, 1000);
	RatioCanvasPt->Divide(2);

	RatioCanvasPt->cd(1);

	for (int i = 0; i <= nSITER; i++){
		RCPPt[0][i] = (TH1D *)UnfoldedPtRebinned[0][i]->Clone();
		RCPPt[0][i]->Divide(UnfoldedPtRebinned[2][i]);
		RCPPt[0][i]->Scale(taa[2]/taa[0]);

		SetName(RCPPt[0][i], Form("R_{CP} Wide p_{T} #SI %i", i));

		SetColor(RCPPt[0][i], i + 1, 20);
		RCPPt[0][i]->Draw("EP SAME");
		RCPPt[0][i]->GetYaxis()->SetRangeUser(0, 2);


	}

	TLine *tmppT = (TLine *)GetLineAtOne(RCPPt[0][0]);
	tmppT->Draw("SAME");

	RatioCanvasPt->cd(2);

	for (int i = 0; i <= nSITER; i++){
		RCPPt[1][i] = (TH1D *)UnfoldedPtRebinned[1][i]->Clone();
		RCPPt[1][i]->Divide(UnfoldedPtRebinned[2][i]);
		RCPPt[1][i]->Scale(taa[2]/taa[1]);

		SetName(RCPPt[1][i], Form("R_{MP} Wide p_{T} #SI %i", i));

		SetColor(RCPPt[1][i], i + 1, 20);
		RCPPt[1][i]->Draw("EP SAME");
		RCPPt[1][i]->GetYaxis()->SetRangeUser(0, 2);
	}

	gPad->BuildLegend();

	tmppT->Draw("SAME");

	RatioCanvasPt->SaveAs(Form("%s/Plots/RatiopT_Wide_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TCanvas *RatioCanvasZ = new TCanvas("RCP Z", "RCP Z", 2000, 1000);
	RatioCanvasZ->Divide(2);

	RatioCanvasZ->cd(1);

	for (int i = 0; i <= nSITER; i++){
		RCPZ[0][i] = (TH1D *)UnfoldedZRebinned[0][i]->Clone();
		RCPZ[0][i]->Divide(UnfoldedZRebinned[2][i]);
		RCPZ[0][i]->Scale(taa[2]/taa[0]);

		SetName(RCPZ[0][i], Form("R_{CP} Wide Z #SI %i", i));

		SetColor(RCPZ[0][i], i + 1, 20);
		RCPZ[0][i]->Draw("EP SAME");
		RCPZ[0][i]->GetYaxis()->SetRangeUser(0, 2);
	}

	gPad->BuildLegend();

	TLine *tmpZ = (TLine *)GetLineAtOne(RCPZ[0][0]);
	tmpZ->Draw("SAME");

	RatioCanvasZ->cd(2);

	for (int i = 0; i <= nSITER; i++){
		RCPZ[1][i] = (TH1D *)UnfoldedZRebinned[1][i]->Clone();
		RCPZ[1][i]->Divide(UnfoldedZRebinned[2][i]);
		RCPZ[1][i]->Scale(taa[2]/taa[1]);

		SetName(RCPZ[1][i], Form("R_{MP} Wide Z #SI %i", i));

		SetColor(RCPZ[1][i], i + 1, 20);
		RCPZ[1][i]->Draw("EP SAME");
		RCPZ[1][i]->GetYaxis()->SetRangeUser(0, 2);
	}

	gPad->BuildLegend();

	tmpZ->Draw("SAME");

	RatioCanvasZ->SaveAs(Form("%s/Plots/RatioZ_Wide_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	// Delta R

	TCanvas *RatioCanvasdR = new TCanvas("RCP dR", "RCP dR", 2000, 1000);
	RatioCanvasdR->Divide(2);

	RatioCanvasdR->cd(1);

	for (int i = 0; i <= nSITER; i++){
		RCPdR[0][i] = (TH1D *)UnfoldeddRRebinned[0][i]->Clone();
		RCPdR[0][i]->Divide(UnfoldeddRRebinned[2][i]);

		SetName(RCPdR[0][i], Form("R_{CP} Wide dR #SI %i", i));

		SetColor(RCPdR[0][i], i + 1, 20);
		RCPdR[0][i]->Draw("EP SAME");
		RCPdR[0][i]->GetYaxis()->SetRangeUser(0, 2);
	}

	gPad->BuildLegend();

	TLine *tmpdR = (TLine *)GetLineAtOne(RCPdR[0][0]);
	tmpdR->Draw("SAME");

	RatioCanvasdR->cd(2);

	for (int i = 0; i <= nSITER; i++){
		RCPdR[1][i] = (TH1D *)UnfoldeddRRebinned[1][i]->Clone();
		RCPdR[1][i]->Divide(UnfoldeddRRebinned[2][i]);

		SetName(RCPdR[1][i], Form("R_{MP} Wide dR #SI %i", i));

		SetColor(RCPdR[1][i], i + 1, 20);
		RCPdR[1][i]->Draw("EP SAME");
		RCPdR[1][i]->GetYaxis()->SetRangeUser(0, 2);
	}

	gPad->BuildLegend();

	tmpdR->Draw("SAME");

	RatioCanvasdR->SaveAs(Form("%s/Plots/RatiodR_Wide_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TCanvas *RatioCanvasPtCompared = new TCanvas("pT Area Vs CS", "pT Area Vs CS", 2100, 700);
	RatioCanvasPtCompared->Divide(3);

	auto legend = new TLegend(0.1, 0.1, 0.9, 0.3);

	for (int cent = 0; cent < 3; cent++){
		RatioCanvasPtCompared->cd(cent+1);
		gPad->SetLogy();
		gPad->SetLeftMargin(0.12);
		for (int i = 0; i <= nSITER; i++){
			SetColor(UnfoldedpTCS[cent][i], i + 1, 24);
			UnfoldedpTCS[cent][i]->Draw("EP SAME");
			UnfoldedPtRebinned[cent][i]->Draw("EP SAME");
			UnfoldedpTCS[cent][i]->GetXaxis()->SetRangeUser(5, 20);
			SetAxisTitles(UnfoldedpTCS[cent][i], "p_{T, Jet} [GeV/#it{c}]", "Unfolded Counts");
			if (cent == 2){
				legend->AddEntry(UnfoldedpTCS[cent][i], "Constituent BG Subtracted", "lep");
				legend->AddEntry(UnfoldedPtRebinned[cent][i], "AreaBased BG Subtracted", "lep");
			}
		}
		if (cent == 2) legend->Draw("same");
	}

	// gPad->BuildLegend();

	RatioCanvasPtCompared->SaveAs(Form("%s/Plots/AreaVsCS_pT.pdf", DirName.Data()));

	TCanvas *RatioCanvasZCompared = new TCanvas("Z Area Vs CS", "Z Area Vs CS", 2100, 700);
	RatioCanvasZCompared->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		RatioCanvasZCompared->cd(cent+1);
		gPad->SetLogy();

		for (int i = 0; i <= nSITER; i++){
			SetColor(UnfoldedZCS[cent][i], i + 1, 24);
			UnfoldedZCS[cent][i]->Draw("EP SAME");
			UnfoldedZRebinned[cent][i]->Draw("EP SAME");
		}
	}

	gPad->BuildLegend();

	RatioCanvasZCompared->SaveAs(Form("%s/Plots/AreaVsCS_Z.pdf", DirName.Data()));

	TCanvas *RatioCanvasRCPPt = new TCanvas("RCP pT Area Vs CS", "RCP pT Area Vs CS", 2100, 700);
	RatioCanvasRCPPt->Divide(2);

	for (int pad = 0; pad < 2; pad++){
		RatioCanvasRCPPt->cd(pad+1);
		SetColor(RCPPtCS[pad][0], 1, 24);
		RCPPtCS[pad][0]->Draw("EP SAME");
		RCPPt[pad][0]->Draw("EP SAME");
	}

	gPad->BuildLegend();

	RatioCanvasRCPPt->SaveAs(Form("%s/Plots/AreaVsCS_RCP_pT.pdf", DirName.Data()));

	TCanvas *RatioCanvasRCPZ = new TCanvas("RCP Z Area Vs CS", "RCP Z Area Vs CS", 2100, 700);
	RatioCanvasRCPZ->Divide(2);

	for (int pad = 0; pad < 2; pad++){
		RatioCanvasRCPZ->cd(pad+1);
		SetColor(RCPZCS[pad][0], 1, 24);
		RCPZCS[pad][0]->Draw("EP SAME");
		RCPZ[pad][0]->Draw("EP SAME");
	}

	gPad->BuildLegend();

	RatioCanvasRCPZ->SaveAs(Form("%s/Plots/AreaVsCS_RCP_Z.pdf", DirName.Data()));

	TCanvas *RatioCanvasRCPdR = new TCanvas("RCP dR Area", "RCP dR Area", 2100, 700);
	RatioCanvasRCPdR->Divide(2);

	for (int pad = 0; pad < 2; pad++){
		RatioCanvasRCPdR->cd(pad+1);
		SetColor(RCPdR[pad][0], 1, 24);
		RCPdR[pad][0]->Draw("EP");
	}

	gPad->BuildLegend();

	RatioCanvasRCPdR->SaveAs(Form("%s/Plots/RCP_dR.pdf", DirName.Data()));

	for (int i = 0; i <= nSITER; i++){
		TFile *f;
		if (isSimpleUnfolding) f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "UPDATE");
		else f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, iteration, njpt_gen_bins_var, nz_gen_bins), "UPDATE");
		Unfolded1D[0][i]->Write();
		Unfolded1D[1][i]->Write();
		Unfolded1D[2][i]->Write();

		Unfolded1DRebinned[0][i]->Write();
		Unfolded1DRebinned[1][i]->Write();
		Unfolded1DRebinned[2][i]->Write();

		UnfoldedPt[0][i]->Write();
		UnfoldedPt[1][i]->Write();
		UnfoldedPt[2][i]->Write();
		UnfoldedZ[0][i]->Write();
		UnfoldedZ[1][i]->Write();
		UnfoldedZ[2][i]->Write();
		Unfolded1DdR[0][i]->Write();
		Unfolded1DdR[1][i]->Write();
		Unfolded1DdR[2][i]->Write();

		UnfoldedPtRebinned[0][i]->Write();
		UnfoldedPtRebinned[1][i]->Write();
		UnfoldedPtRebinned[2][i]->Write();
		UnfoldedZRebinned[0][i]->Write();
		UnfoldedZRebinned[1][i]->Write();
		UnfoldedZRebinned[2][i]->Write();
		UnfoldeddRRebinned[0][i]->Write();
		UnfoldeddRRebinned[1][i]->Write();
		UnfoldeddRRebinned[2][i]->Write();

		RCPPt[0][i]->Write();
		RCPPt[1][i]->Write();
		RCPZ[0][i]->Write();
		RCPZ[1][i]->Write();
		RCPdR[0][i]->Write();
		RCPdR[1][i]->Write();
		
		f->Close();
	}


}