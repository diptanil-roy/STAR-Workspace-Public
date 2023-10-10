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

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

void ComparePtSmearinGenPtBin(){

	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);

  	const int filenum = 3;

  	int color[filenum] = {kViolet, kAzure+10, kRed};

  	int marker[filenum] = {20, 24, 20};

  	TH1D *dPt[filenum][3][nbins_jpt];

  	TString fname[filenum];

  	fname[0] = "Jan26_FONLL_MC_CombinedFile4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
  	fname[1] = "Feb14_FONLL_HI4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_7.root";
  	fname[2] = "Jan26_FONLL_HI_CombinedFile_pTSmearedFromMC4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";

  	TString LegAppend[filenum];
	LegAppend[0] = "MC (FONLL)";
	LegAppend[1] = "HI (FONLL)";
	LegAppend[2] = "HI (FONLL MC)";

  	for (int fi = 0; fi < filenum; fi++){
		// cout << RFile[fi].Data() << endl;
		TFile *f = new TFile(fname[fi].Data());
		f->cd();
		for (int i = 0; i < 3; i++){
			for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
				dPt[fi][i][ptbin] = (TH1D *)gDirectory->Get(Form("hDiffJetPt_%i_%i", i, ptbin));
				SetName(dPt[fi][i][ptbin], Form("hDiffJetPt_%i_%i_%i", fi, i, ptbin));
				dPt[fi][i][ptbin]->SetDirectory(0);
			}
		}
		f->Close();
	}
	

	TCanvas *c[3];
	TLegend *l[3][nbins_jpt];

	TF1 * f = new TF1("f","2.*gaus(x,[0],[1],[2])*ROOT::Math::normal_cdf([3]*x,1,0)", -20, 30);
	

	TF1 *g[3][nbins_jpt];

	for (int cent = 0; cent < 3; cent++){
		c[cent] = new TCanvas(Form("c_%i", cent), Form("c_%i", cent), 3000, 400);
		c[cent]->Divide(nbins_jpt);

		// cout << "Canvas " << cent << endl;

		for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
			c[cent]->cd(ptbin+1);
			gPad->SetLogy();

			l[cent][ptbin] =  new TLegend(0.4,0.8,0.9,0.9);
			// cout << "Pad " << ptbin << endl;

			for (int fi = 0; fi < 2; fi++){
				SetColor(dPt[fi][cent][ptbin], color[fi], marker[fi]);
				dPt[fi][cent][ptbin]->GetXaxis()->SetRangeUser(-30, 40);
				dPt[fi][cent][ptbin]->GetYaxis()->SetRangeUser(0.01, 5*pow(10, 6));
				dPt[fi][cent][ptbin]->Draw("EP SAME");
				l[cent][ptbin]->AddEntry(dPt[fi][cent][ptbin],Form("%s = %.2f, %.2f", LegAppend[fi].Data(), dPt[fi][cent][ptbin]->GetMean(), dPt[fi][cent][ptbin]->GetRMS()),"p");
				if (fi == 0) SetName(dPt[fi][cent][ptbin], Form("Gen p_{T} #in  [%.1f, %.1f], Cent = %i", binning_jpt[ptbin], binning_jpt[ptbin+1], cent));
			}

			l[cent][ptbin]->SetTextSize(0.03);
			l[cent][ptbin]->Draw("SAME");
		}
		c[cent]->SaveAs(Form("JetQA/DeltapT_%i.pdf", cent));
	}

	/*
	TCanvas *d[3];

	double skewfactor[3][nbins_jpt] = {{0.1, 0.1, 0.1, 0.05, 0.05, 0.05},
					 				   {0.2, 0.2, 0.2, 0.1, 0.1, 0.1},
					                   {0.6, 0.5, 0.4, 0.4, 0.2, 0.1}
					                  };

	for (int cent = 0; cent < 3; cent++){
		d[cent] = new TCanvas(Form("d_%i", cent), Form("d_%i", cent), 3000, 400);
		d[cent]->Divide(nbins_jpt);

		for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
			d[cent]->cd(ptbin+1);
			gPad->SetLogy();

			g[cent][ptbin] = (TF1 *)f->Clone(Form("FitForm_%i_%i", cent, ptbin));
			g[cent][ptbin]->SetParameters(10000, dPt[0][cent][ptbin]->GetMean(), dPt[0][cent][ptbin]->GetRMS(), skewfactor[cent][ptbin]);
			// g[cent][ptbin]->SetParLimits(3, -0.8, 0.8);

			dPt[0][cent][ptbin]->Fit(Form("FitForm_%i_%i", cent, ptbin), "ELM");
			dPt[0][cent][ptbin]->Draw();
		}

		d[cent]->SaveAs(Form("JetQA/DeltapT_Fit_%i.pdf", cent));
	}

	TFile *PtSmear = new TFile("PtSmear_Fit.root", "RECREATE");
	PtSmear->cd();

	for (int cent = 0; cent < 3; cent++){
		for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
			dPt[0][cent][ptbin]->Write();
			g[cent][ptbin]->Write();
		}
	}

	PtSmear->Close();

	*/
}
