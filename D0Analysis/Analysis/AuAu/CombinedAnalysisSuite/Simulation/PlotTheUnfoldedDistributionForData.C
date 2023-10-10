using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void PlotTheUnfoldedDistributionForData(TString DirName = "tmp4", int D0pTLow = 5, int SUPERITERATION = 20, int iteration = 4, TString NonClosureDir = "", bool isSimpleUnfolding = true){

	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);

  	// Scaling by TAA

	double taa[3] = {941.23714*1.0318440e+08, 391.35550*3.2123506e+08, 56.62475*4.6679240e+08};

	const int nSITER = SUPERITERATION;

	TH2D *Truth[3];
	TH2D *Measured[3];

	TH2D *Unfolded[3][nSITER+1];

	for (int i = 0; i <= nSITER; i++){
		TFile *f;
		
		if (isSimpleUnfolding) f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
		else f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, iteration, njpt_gen_bins_var, nz_gen_bins), "READ");
		
		cout << f->GetName() << endl;
		f->cd();
		if (i == 0){
			for (int cent = 0; cent < 3; cent++){
				Measured[cent] = (TH2D *)f->Get(Form("Measured_%i", cent));
				Measured[cent]->SetDirectory(0);
				cout << Measured[cent]->GetName() << endl;
			}
		}
		for (int cent = 0; cent < 3; cent++){
			Unfolded[cent][i] = (TH2D *)f->Get(Form("Unfolded_%i", cent));
			SetName(Unfolded[cent][i], Form("Unfolded Cent = %i SI = %i Iter = %i", cent, i, iteration));
			Unfolded[cent][i]->SetDirectory(0);

			Unfolded[cent][i]->GetXaxis()->SetRangeUser(D0pTLow, 30);

			cout << "SI = " << i << "\t" << Unfolded[cent][i]->GetName() << endl;
		}

		f->Close();
	}

	// if (NonClosureDir != ""){

	// 	TFile *SummaryFile = new TFile(Form("%s/Summary.root", NonClosureDir.Data()), "READ");
	// 	SummaryFile->cd();

	// 	TH2D *tmpNonClosure;
	// 		for (int cent = 0; cent < 3; cent++){
	// 			for (int i = 0; i <= nSITER; i++){
	// 				tmpNonClosure = (TH2D *)SummaryFile->Get(Form("Non Closure Cent %i", cent));
	// 				Unfolded[cent][i]->Divide(tmpNonClosure);
	// 			}
	// 	}
	// }

	cout << "Imported all histograms" << endl;

	int lowptbin = Unfolded[0][0]->GetXaxis()->FindBin(5);
	int highptbin = Unfolded[0][0]->GetXaxis()->FindBin(20);

	TH1D *MeasuredPt[3];
	TH1D *UnfoldedPt[3][nSITER+1];

	TH1D *MeasuredZ[3];
	TH1D *UnfoldedZ[3][nSITER+1];

	TH1D *UnfoldedPtRebinned[3][nSITER+1];
	TH1D *UnfoldedZRebinned[3][nSITER+1];

	for (int cent = 0; cent < 3; cent++){
		MeasuredPt[cent] = (TH1D *)Measured[cent]->ProjectionX();
		SetName(MeasuredPt[cent], Form("Measured p_{T} %i", cent));	
		
		MeasuredZ[cent] = (TH1D *)Measured[cent]->ProjectionY();
		SetName(MeasuredZ[cent], Form("Measured Z %i", cent));

		for (int i = 0; i <= nSITER; i++){
			UnfoldedPt[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionX();
			UnfoldedPt[cent][i]->GetXaxis()->SetRangeUser(5, 30);
			SetName(UnfoldedPt[cent][i], Form("Unfolded p_{T} Cent = %i SI = %i Iter = %i", cent, i, iteration));
			UnfoldedZ[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionY("", lowptbin, highptbin);
			SetName(UnfoldedZ[cent][i], Form("Unfolded Z Cent = %i SI = %i Iter = %i", cent, i, iteration));
			cout << "Division Done" << endl;

			UnfoldedPtRebinned[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionX();
			UnfoldedPtRebinned[cent][i]->GetXaxis()->SetRangeUser(5, 30);
			SetName(UnfoldedPtRebinned[cent][i], Form("Unfolded p_{T} Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
			UnfoldedZRebinned[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionY("", lowptbin, highptbin);
			SetName(UnfoldedZRebinned[cent][i], Form("Unfolded Z Cent = %i SI = %i Iter = %i Closed", cent, i, iteration));
		}
	}

	cout << "Made 1D Histograms" << endl;

	if (NonClosureDir != ""){

		TFile *SummaryFile = new TFile(Form("%s/Summary.root", NonClosureDir.Data()), "READ");
		SummaryFile->cd();

		TH1D *tmpPtNonClosure;
		TH1D *tmpZNonClosure;

		for (int cent = 0; cent < 3; cent++){
			for (int i = 0; i <= nSITER; i++){
				tmpPtNonClosure = (TH1D *)SummaryFile->Get(Form("Non Closure Pt Cent %i SI %i", cent, i));
				tmpZNonClosure = (TH1D *)SummaryFile->Get(Form("Non Closure Z Cent %i SI %i", cent, i));

				UnfoldedPtRebinned[cent][i]->Divide(tmpPtNonClosure);
				UnfoldedZRebinned[cent][i]->Divide(tmpZNonClosure);
			}
		}
	}

	

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
			// UnfoldedPtRebinned[cent][i] = (TH1D *)UnfoldedPt[cent][i]->Clone();
			// ProcessSpectra(UnfoldedPtRebinned[cent][i]);
			// UnfoldedPtRebinned[cent][i]->Scale(1./taa[cent]);
			
			SetColor(UnfoldedPtRebinned[cent][i], i + 1, 20);
			UnfoldedPtRebinned[cent][i]->Draw("EP SAME");
			UnfoldedPtRebinned[cent][i]->GetYaxis()->SetRangeUser(pow(10,-1), pow(10, 5));
		}

		// MeasuredPt[cent]->Rebin(2);

		gPad->BuildLegend();
	}

	Pt->SaveAs(Form("%s/Plots/UnfoldedpT_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

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

	Z->SaveAs(Form("%s/Plots/UnfoldedZ_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TH1D *RCPPt[3][nSITER+1];
	TH1D *RCPZ[3][nSITER+1];

	TCanvas *RatioCanvasPt = new TCanvas("RCP pT", "RCP pT", 2000, 1000);
	RatioCanvasPt->Divide(2);

	RatioCanvasPt->cd(1);

	for (int i = 0; i <= nSITER; i++){
		RCPPt[0][i] = (TH1D *)UnfoldedPtRebinned[0][i]->Clone();
		RCPPt[0][i]->Divide(UnfoldedPtRebinned[2][i]);
		RCPPt[0][i]->Scale(taa[2]/taa[0]);

		SetName(RCPPt[0][i], Form("R_{CP} p_{T} #SI %i", i));

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

		SetName(RCPPt[1][i], Form("R_{MP} p_{T} #SI %i", i));

		SetColor(RCPPt[1][i], i + 1, 20);
		RCPPt[1][i]->Draw("EP SAME");
		RCPPt[1][i]->GetYaxis()->SetRangeUser(0, 2);
	}

	gPad->BuildLegend();

	tmppT->Draw("SAME");

	RatioCanvasPt->SaveAs(Form("%s/Plots/RatiopT_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TCanvas *RatioCanvasZ = new TCanvas("RCP Z", "RCP Z", 2000, 1000);
	RatioCanvasZ->Divide(2);

	RatioCanvasZ->cd(1);

	for (int i = 0; i <= nSITER; i++){
		RCPZ[0][i] = (TH1D *)UnfoldedZRebinned[0][i]->Clone();
		RCPZ[0][i]->Divide(UnfoldedZRebinned[2][i]);
		RCPZ[0][i]->Scale(taa[2]/taa[0]);

		SetName(RCPZ[0][i], Form("R_{CP} Z #SI %i", i));

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

		SetName(RCPZ[1][i], Form("R_{MP} Z #SI %i", i));

		SetColor(RCPZ[1][i], i + 1, 20);
		RCPZ[1][i]->Draw("EP SAME");
		RCPZ[1][i]->GetYaxis()->SetRangeUser(0, 2);
	}

	gPad->BuildLegend();

	tmpZ->Draw("SAME");

	RatioCanvasZ->SaveAs(Form("%s/Plots/RatioZ_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	for (int i = 0; i <= nSITER; i++){
		TFile *f;
		if (isSimpleUnfolding) f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, 0, njpt_gen_bins_var, nz_gen_bins), "UPDATE");
		else f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, iteration, njpt_gen_bins_var, nz_gen_bins), "UPDATE");
		UnfoldedPt[0][i]->Write();
		UnfoldedPt[1][i]->Write();
		UnfoldedPt[2][i]->Write();
		UnfoldedZ[0][i]->Write();
		UnfoldedZ[1][i]->Write();
		UnfoldedZ[2][i]->Write();
		UnfoldedPtRebinned[0][i]->Write();
		UnfoldedPtRebinned[1][i]->Write();
		UnfoldedPtRebinned[2][i]->Write();
		UnfoldedZRebinned[0][i]->Write();
		UnfoldedZRebinned[1][i]->Write();
		UnfoldedZRebinned[2][i]->Write();
		RCPPt[0][i]->Write();
		RCPPt[1][i]->Write();
		RCPZ[0][i]->Write();
		RCPZ[1][i]->Write();
		f->Close();
	}


}