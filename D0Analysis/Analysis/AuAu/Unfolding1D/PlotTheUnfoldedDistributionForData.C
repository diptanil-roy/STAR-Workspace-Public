using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void PlotTheUnfoldedDistributionForData(TString DirName = "tmp4", int SUPERITERATION = 20, int iteration = 4, TString NonClosureDir = "", bool isSimpleUnfolding = false){

	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);

  	// Scaling by TAA

	double taa[3] = {941.23714*1.0318440e+08, 391.35550*3.2123506e+08, 56.62475*4.6679240e+08};

	const int nSITER = SUPERITERATION;

	TH2D *Truth[3];
	TH2D *Measured[3];

	TH1D *Measured1D[3];
	TH1D *Unfolded1D[3];

	TFile *f;
		
	f = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins), "READ");
	
	cout << f->GetName() << endl;
	f->cd();

	for (int cent = 0; cent < 3; cent++){
		Measured1D[cent] = (TH1D *)f->Get(Form("Measured1D_%i", cent));
		Measured1D[cent]->SetDirectory(0);
		cout << Measured1D[cent]->GetName() << endl;

		Unfolded1D[cent] = (TH1D *)f->Get(Form("Unfolded1D_%i", cent));
		SetName(Unfolded1D[cent], Form("Unfolded 1D Cent = %i SI = %i Iter = %i", cent, 0, iteration));
		Unfolded1D[cent]->SetDirectory(0);

		cout << Unfolded1D[cent]->GetName() << endl;
	}

	cout << "Imported all histograms" << endl;

	cout << "Made 1D Histograms" << endl;

	TH1D *UnfoldedPtRebinned[3];

	TCanvas *Pt = new TCanvas("pT", "pT", 2000, 1000);
	Pt->Divide(3);

	for (int cent = 0; cent < 3; cent++){
		Pt->cd(cent+1);
		gPad->SetLogy();

		UnfoldedPtRebinned[cent] = (TH1D *)Unfolded1D[cent]->Rebin(nbins_jpt, Form("%s_Rebinned", Unfolded1D[cent]->GetName()), binning_jpt);
		// ProcessSpectra(UnfoldedPtRebinned[cent]);
		UnfoldedPtRebinned[cent]->Draw("EP SAME");

		gPad->BuildLegend();
	}

	Pt->SaveAs(Form("%s/Plots/UnfoldedpT_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TH1D *RCPPt[3];

	TCanvas *RatioCanvasPt = new TCanvas("RCP pT", "RCP pT", 2000, 1000);
	RatioCanvasPt->Divide(2);

	RatioCanvasPt->cd(1);

	RCPPt[0] = (TH1D *)UnfoldedPtRebinned[0]->Clone();
	RCPPt[0]->Divide(UnfoldedPtRebinned[2]);
	RCPPt[0]->Scale(taa[2]/taa[0]);

	SetName(RCPPt[0], Form("R_{CP} p_{T} #SI %i", 0));

	SetColor(RCPPt[0], 1, 20);
	RCPPt[0]->Draw("EP SAME");
	RCPPt[0]->GetYaxis()->SetRangeUser(0, 1);

	TLine *tmppT = (TLine *)GetLineAtOne(RCPPt[0]);
	tmppT->Draw("SAME");

	RatioCanvasPt->cd(2);

	RCPPt[1] = (TH1D *)UnfoldedPtRebinned[1]->Clone();
	RCPPt[1]->Divide(UnfoldedPtRebinned[2]);
	RCPPt[1]->Scale(taa[2]/taa[1]);

	SetName(RCPPt[1], Form("R_{MP} p_{T} #SI %i", 0));

	SetColor(RCPPt[1], 1, 20);
	RCPPt[1]->Draw("EP SAME");
	RCPPt[1]->GetYaxis()->SetRangeUser(0, 1);

	tmppT->Draw("SAME");

	gPad->BuildLegend();

	RatioCanvasPt->SaveAs(Form("%s/Plots/RatiopT_%s.pdf", DirName.Data(), (NonClosureDir != "") ? "NC" : "Reg"));

	TFile *g;

	g = new TFile(Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins), "UPDATE");

	UnfoldedPtRebinned[0]->Write();
	UnfoldedPtRebinned[1]->Write();
	UnfoldedPtRebinned[2]->Write();

	RCPPt[0]->Write();
	RCPPt[1]->Write();

	g->Close();

}