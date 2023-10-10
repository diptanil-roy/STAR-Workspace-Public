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

void PlotTheUnfoldedDistributionForMC(TString DirName = "tmp4", int SUPERITERATION = 20, int iteration = 4){

	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);

	const int nSITER = SUPERITERATION;

	int colors[nSITER+1];
	Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,nSITER+1);
	colors[0] = kBlack;
	for (int i=1; i<=nSITER; i++) colors[i] = FI+i;

	TH2D *Truth[3];
	TH2D *Measured[3];

	TH2D *Unfolded[3][nSITER+1];

	for (int i = 0; i <= nSITER; i++){
		TFile *f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), i, iteration, njpt_gen_bins_var, nz_gen_bins), "READ");
		
		cout << f->GetName() << endl;
		f->cd();
		if (i == 0){
			for (int cent = 0; cent < 3; cent++){
				Truth[cent] = (TH2D *)f->Get(Form("Truth_%i", cent));
				Measured[cent] = (TH2D *)f->Get(Form("Measured_%i", cent));
				Truth[cent]->SetDirectory(0);
				Measured[cent]->SetDirectory(0);

				cout << Truth[cent]->GetName() << endl;
				cout << Measured[cent]->GetName() << endl;
			}
		}
		for (int cent = 0; cent < 3; cent++){
			Unfolded[cent][i] = (TH2D *)f->Get(Form("Unfolded_%i", cent));
			SetName(Unfolded[cent][i], Form("Unfolded Cent = %i SI = %i Iter = %i", cent, i, iteration));
			Unfolded[cent][i]->SetDirectory(0);

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

	TH1D *CombinedChi2[3];

	for (int cent = 0; cent < 3; cent++){
		Chi2Pt[cent] = new TH1D(Form("Chi2Pt_%i", cent), Form("Chi2Pt_%i", cent), nSITER+1, -0.5, nSITER+0.5);
		Chi2Z[cent] = new TH1D(Form("Chi2Z_%i", cent), Form("Chi2Z_%i", cent), nSITER+1, -0.5, nSITER+0.5);
		CombinedChi2[cent] = new TH1D(Form("CombinedChi2_%i", cent), Form("CombinedChi2_%i", cent), nSITER+1, -0.5, nSITER+0.5);
	}

	for (int cent = 0; cent < 3; cent++){
		TruthPt[cent] = (TH1D *)Truth[cent]->ProjectionX();
		SetName(TruthPt[cent], Form("Truth p_{T} %i", cent));
		MeasuredPt[cent] = (TH1D *)Measured[cent]->ProjectionX();
		SetName(MeasuredPt[cent], Form("Measured p_{T} %i", cent));	
		
		TruthZ[cent] = (TH1D *)Truth[cent]->ProjectionY();
		SetName(TruthZ[cent], Form("Truth Z %i", cent));
		MeasuredZ[cent] = (TH1D *)Measured[cent]->ProjectionY();
		SetName(MeasuredZ[cent], Form("Measured Z %i", cent));

		for (int i = 0; i <= nSITER; i++){
			CombinedChi2[cent]->SetBinContent(i+1, Chi2Func2D(Unfolded[cent][i], Truth[cent]));

			UnfoldedPt[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionX();
			SetName(UnfoldedPt[cent][i], Form("Unfolded p_{T} Cent = %i SI = %i Iter = %i", cent, i+1, iteration));
			Chi2Pt[cent]->SetBinContent(i+1, Chi2Func(UnfoldedPt[cent][i], TruthPt[cent]));

			UnfoldedZ[cent][i] = (TH1D *)Unfolded[cent][i]->ProjectionY();
			// cout << "Cent = " << cent << "\t" << "SI = " << i << "\t" << Unfolded[cent][i]->GetNbinsX() << "\t" << Unfolded[cent][i]->GetNbinsY() << endl;
			SetName(UnfoldedZ[cent][i], Form("Unfolded Z Cent = %i SI = %i Iter = %i", cent, i+1, iteration));
			Chi2Z[cent]->SetBinContent(i+1, Chi2Func(UnfoldedZ[cent][i], TruthZ[cent]));
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
			SetColor(UnfoldedPt[cent][i], colors[i], 20);
			UnfoldedPt[cent][i]->Draw("EP SAME");
			UnfoldedPt[cent][i]->GetYaxis()->SetRangeUser(pow(10,0), pow(10, 7));
		}

		// TruthPt[cent]->Rebin(5);

		SetColor(TruthPt[cent], kBlack, 29);
		TruthPt[cent]->Draw("EP SAME");

		gPad->BuildLegend();
	}

	Pt->SaveAs(Form("%s/Plots/UnfoldedpT.pdf", DirName.Data()));

	TCanvas *Z = new TCanvas("Z", "Z", 2000, 1000);
	Z->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Z->cd(cent+1);
		gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			// UnfoldedZ[cent][i]->Rebin(5);
			SetColor(UnfoldedZ[cent][i], colors[i], 20);
			UnfoldedZ[cent][i]->Draw("EP SAME");
			UnfoldedZ[cent][i]->GetYaxis()->SetRangeUser(pow(10,0), pow(10, 7));
			
		}

		SetColor(TruthZ[cent], kBlack, 29);
		TruthZ[cent]->Draw("EP SAME");

		gPad->BuildLegend();
	}

	Z->SaveAs(Form("%s/Plots/UnfoldedZ.pdf", DirName.Data()));

	TH1D *RatioPt[3][nSITER+1];
	TH1D *RatioZ[3][nSITER+1];

	TCanvas *RatioCanvasPt = new TCanvas("Ratio Canvas pT", "Ratio Canvas pT", 2000, 1000);
	RatioCanvasPt->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		RatioCanvasPt->cd(cent+1);
		// gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			RatioPt[cent][i] = (TH1D *)UnfoldedPt[cent][i]->Clone();
			RatioPt[cent][i]->Divide(TruthPt[cent]);

			SetColor(RatioPt[cent][i], colors[i], 20);
			RatioPt[cent][i]->Draw("EP SAME");
			RatioPt[cent][i]->GetYaxis()->SetRangeUser(0, 2);
		}

		gPad->BuildLegend();
	
		// TH1D *tmp = (TH1D *)GetLineAtOne(TruthPt[cent]);
		// tmp->Draw("P HIST SAME");

		TLine *tmp = (TLine *)GetLineAtOne(TruthPt[cent]);
		tmp->Draw("SAME");
	}

	RatioCanvasPt->SaveAs(Form("%s/Plots/RatiopT.pdf", DirName.Data()));

	TCanvas *RatioCanvasZ = new TCanvas("Ratio Canvas Z", "Ratio Canvas Z", 2000, 1000);
	RatioCanvasZ->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		RatioCanvasZ->cd(cent+1);
		// gPad->SetLogy();
		for (int i = 0; i <= nSITER; i++){
			// cout << "Cent = " << cent << " SI = " << i << " Bins = " << UnfoldedZ[cent][i]->GetNbinsX() << "\t" << TruthZ[cent]->GetNbinsX() << endl;
			RatioZ[cent][i] = (TH1D *)UnfoldedZ[cent][i]->Clone();
			RatioZ[cent][i]->Divide(TruthZ[cent]);

			SetColor(RatioZ[cent][i], colors[i], 20);
			RatioZ[cent][i]->Draw("EP SAME");
			RatioZ[cent][i]->GetYaxis()->SetRangeUser(0, 5);
		}

		gPad->BuildLegend();

		TLine *tmp = (TLine *)GetLineAtOne(TruthZ[cent]);
		tmp->Draw("SAME");
	}

	RatioCanvasZ->SaveAs(Form("%s/Plots/RatioZ.pdf", DirName.Data()));

	TCanvas *Chi2CanvasPt = new TCanvas("Chi2 Canvas pT", "Chi2 Canvas pT", 2000, 1000);
	Chi2CanvasPt->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Chi2CanvasPt->cd(cent+1);
		Chi2Pt[cent]->Draw("HIST SAME");
		// Chi2Pt[cent]->GetYaxis()->SetRangeUser(0, 100);
		gPad->SetLogy();
	}

	Chi2CanvasPt->SaveAs(Form("%s/Plots/Chi2pT.pdf", DirName.Data()));

	TCanvas *Chi2CanvasZ = new TCanvas("Chi2 Canvas Z", "Chi2 Canvas Z", 2000, 1000);
	Chi2CanvasZ->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Chi2CanvasZ->cd(cent+1);
		Chi2Z[cent]->Draw("HIST SAME");
		// Chi2Z[cent]->GetYaxis()->SetRangeUser(0, 100);
		gPad->SetLogy();
	}

	Chi2CanvasZ->SaveAs(Form("%s/Plots/Chi2Z.pdf", DirName.Data()));

	TCanvas *Chi2Canvas = new TCanvas("Chi2 Canvas", "Chi2 Canvas", 2000, 1000);
	Chi2Canvas->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Chi2Canvas->cd(cent+1);
		CombinedChi2[cent]->Draw("HIST SAME");
		// Chi2Z[cent]->GetYaxis()->SetRangeUser(0, 100);
		gPad->SetLogy();

		cout << CombinedChi2[cent]->GetBinCenter(CombinedChi2[cent]->GetMinimumBin()) << endl;
	}

	Chi2Canvas->SaveAs(Form("%s/Plots/CombinedChi2.pdf", DirName.Data()));

	TFile *SummaryStat = new TFile(Form("%s/Plots/Summary.root", DirName.Data()), "RECREATE");
	SummaryStat->cd();

	for (int cent = 0; cent < 3; cent++){
		TruthPt[cent]->Write(Form("Truth Pt Cent %i", cent));
		TruthZ[cent]->Write(Form("Truth Z Cent %i", cent));
		Truth[cent]->Write(Form("Truth Cent %i", cent));

		Chi2Pt[cent]->Write(Form("Chi2 Pt Cent %i", cent));
		Chi2Z[cent]->Write(Form("Chi2 Z Cent %i", cent));
		CombinedChi2[cent]->Write(Form("Combined Chi2 Cent %i", cent));
	}
	for (int cent = 0; cent < 3; cent++){
		for (int i = 0; i <= nSITER; i++){
			UnfoldedPt[cent][i]->Write(Form("Unfolded Pt Cent %i SI %i", cent, i));
			RatioPt[cent][i]->Write(Form("Non Closure Pt Cent %i SI %i", cent, i));
			UnfoldedZ[cent][i]->Write(Form("Unfolded Z Cent %i SI %i", cent, i));
			RatioZ[cent][i]->Write(Form("Non Closure Z Cent %i SI %i", cent, i));
			Unfolded[cent][i]->Write(Form("Unfolded Cent %i SI %i", cent, i));
		}
	} 
	cout << "Wrote to file " << SummaryStat->GetName() << endl;
	SummaryStat->Close();
}