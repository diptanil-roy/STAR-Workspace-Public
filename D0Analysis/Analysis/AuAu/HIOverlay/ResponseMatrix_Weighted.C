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


#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"

#pragma link C++ class vector<int> +;

using namespace std;
// using namespace RooFit;

// #pragma link C++ class StJetTreeStruct+;

// #pragma link C++ class vector<float> +;
// #pragma link C++ class vector<vector<float> >+;
// #pragma link C++ class vector<int> +;
// #pragma link C++ class vector<vector<int> >+;
#endif


#include "BinDef.h"

void ResponseMatrix_Weighted(double power = 0.0){
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	double numberofentries = 10000000;
	// double numberofentries = 100;

	TFile *f = new TFile("Plots.root");
	TH1D *NEntries[2];
	TH1D *VzMC[3][2];
	TH1D *MCZ[3][2];
	TH1D *RecoZ[3][2];
	THn  *Z[3][2];
	TH2D *ResponsePt[3][2];

	TH1D *MCD0JetPtBeforeFiducialCuts[3][11][2];
	TH1D *MCD0JetPtAfterFiducialCuts[3][11][2];

	TH1D *RecoD0JetPtBeforeFiducialCuts[3][11][2];
	TH1D *RecoD0JetPtAfterFiducialCuts[3][11][2];

	double ratioofeventsthatpassveto[2] = {1070./1000000., 11114./1000000.}; // For scaling, we need to make sure that vetoed events are also included.
  	double crosssection[2] = {64.84/66.809, 1.969/66.809};

  	for (int bin = 0; bin < 2; bin++){
		f->cd(Form("pthatbin_%i", bin));
		// gROOT->ProcessLine(".ls");
		NEntries[bin] = (TH1D *)gDirectory->Get("NEvents");

		for (int i = 0; i < 3; i++){
		  VzMC[i][bin] = (TH1D *)gDirectory->Get(Form("VzMC_%i", i));

		  double scale = crosssection[bin]*ratioofeventsthatpassveto[bin]/VzMC[i][bin]->GetEntries();
		  cout << scale << endl;

		  MCZ[i][bin] = (TH1D *)gDirectory->Get(Form("MCZ_%i", i));
		  RecoZ[i][bin] = (TH1D *)gDirectory->Get(Form("RecoZ_%i", i));
		  Z[i][bin] = (THnD *)gDirectory->Get(Form("hD0MCZD0PtJetPtRecoD0PtJetPt_%i", i));
		  ResponsePt[i][bin] = (TH2D *)gDirectory->Get(Form("hResponsePt_%i", i));

		  for (int cut = 0; cut < 11; cut++){
		  	MCD0JetPtBeforeFiducialCuts[i][cut][bin] = (TH1D *)gDirectory->Get(Form("hMCD0JetPtBeforeFiducialCuts_%i_%i", i, cut));
		  	MCD0JetPtAfterFiducialCuts[i][cut][bin] = (TH1D *)gDirectory->Get(Form("hMCD0JetPtAfterFiducialCuts_%i_%i", i, cut));
		  }

		  for (int cut = 0; cut < 1; cut++){
		  	RecoD0JetPtBeforeFiducialCuts[i][cut][bin] = (TH1D *)gDirectory->Get(Form("hRecoD0JetPtBeforeFiducialCuts_%i_%i", i, cut));
		  	RecoD0JetPtAfterFiducialCuts[i][cut][bin] = (TH1D *)gDirectory->Get(Form("hRecoD0JetPtAfterFiducialCuts_%i_%i", i, cut));
		  }
		  
		  MCZ[i][bin]->Scale(scale);
		  RecoZ[i][bin]->Scale(scale);
		  Z[i][bin]->Scale(scale);
		  ResponsePt[i][bin]->Scale(scale);
		  for (int cut = 0; cut < 11; cut++){
		  	MCD0JetPtBeforeFiducialCuts[i][cut][bin]->Scale(scale);
		  	MCD0JetPtAfterFiducialCuts[i][cut][bin]->Scale(scale);
		  }

		  for (int cut = 0; cut < 1; cut++){
		  	RecoD0JetPtBeforeFiducialCuts[i][cut][bin]->Scale(scale);
		  	RecoD0JetPtAfterFiducialCuts[i][cut][bin]->Scale(scale);
		  }
		}
	}

	for (int i = 0; i < 3; i++){
		MCZ[i][0]->Add(MCZ[i][1]);
		RecoZ[i][0]->Add(RecoZ[i][1]);
		Z[i][0]->Add(Z[i][1]);
		ResponsePt[i][0]->Add(ResponsePt[i][1]);
		for (int cut = 0; cut < 11; cut++){
			MCD0JetPtBeforeFiducialCuts[i][cut][0]->Add(MCD0JetPtBeforeFiducialCuts[i][cut][1]);
			MCD0JetPtAfterFiducialCuts[i][cut][0]->Add(MCD0JetPtAfterFiducialCuts[i][cut][1]);
		}
		for (int cut = 0; cut < 1; cut++){
			RecoD0JetPtBeforeFiducialCuts[i][cut][0]->Add(RecoD0JetPtBeforeFiducialCuts[i][cut][1]);
			RecoD0JetPtAfterFiducialCuts[i][cut][0]->Add(RecoD0JetPtAfterFiducialCuts[i][cut][1]);
		}
	}

	TH1D *EfficiencyRatioAtGenLevel[3];

	for (int i = 0; i < 3; i++){
		EfficiencyRatioAtGenLevel[i] = (TH1D *)MCD0JetPtAfterFiducialCuts[i][6][0]->Clone(Form("EfficiencyRatioAtGenLevel_%i", i));
		EfficiencyRatioAtGenLevel[i]->SetNameTitle(Form("EfficiencyRatioAtGenLevel_%i", i), Form("EfficiencyRatioAtGenLevel_%i", i));
		EfficiencyRatioAtGenLevel[i]->Divide(MCD0JetPtBeforeFiducialCuts[i][6][0]);
	}

	TH1D *EfficiencyRatioAtRecoLevel[3]; //This estimates the contribution of particle level jets with |eta| < 0.6 in the detector level within fiducial acceptance

	for (int i = 0; i < 3; i++){
		EfficiencyRatioAtRecoLevel[i] = (TH1D *)RecoD0JetPtAfterFiducialCuts[i][0][0]->Clone(Form("EfficiencyRatioAtRecoLevel_%i", i));
		EfficiencyRatioAtRecoLevel[i]->SetNameTitle(Form("EfficiencyRatioAtRecoLevel_%i", i), Form("EfficiencyRatioAtRecoLevel_%i", i));
		EfficiencyRatioAtRecoLevel[i]->Divide(RecoD0JetPtBeforeFiducialCuts[i][0][0]);
	}

	THn  *ZResampled[3];

	int nBinsDaug[nDimDaug] = {nBinsMCZ, nBinsMCZ, nBinsMCD0Pt, nBinsMCJetPt, nBinsRecoD0Pt, nBinsRecoJetPt, nBinsEta, nBinsEta}; //MC Z, Reco Z, MC D0 Pt, MC Jet Pt, Reco D0 Pt, Reco Jet Pt

	TF1 *f1 = new TF1("f1",Form("pow(1-x, %.1f)", power),0, 1);
	TH1D *RequiredZ = (TH1D *)f1->GetHistogram();
	// cout << RequiredZ->GetNbinsX() << endl;
	RequiredZ->Rebin(10);
	// RequiredZ->Scale(RecoZ[0][0]->Integral()/RequiredZ->Integral());
	// RequiredZ->SetBinContent(1, 0);

	// TCanvas *a = new TCanvas("a", "a", 800, 800);
	// gPad->SetLogy();
	// RequiredZ->Draw();
	// MCZ[0][0]->SetLineColor(kBlack);
	// MCZ[0][0]->SetMarkerColor(kBlack);
	// MCZ[0][0]->SetNameTitle("MC", "MC");
	// RecoZ[0][0]->Draw("SAME");
	// RecoZ[0][0]->SetNameTitle("Reco", "Reco");
	// RequiredZ->GetYaxis()->SetRangeUser(pow(10, -9), pow(10, 0));

	TFile *g2 = new TFile("../SPlotFrameWork/ApplyWeights/Histograms3_D01_10GeV.root");
    g2->cd();

    TH1D *DataZ[3];
    DataZ[0] = (TH1D *)gDirectory->Get("Z_0_10");
    DataZ[1] = (TH1D *)gDirectory->Get("Z_10_40");
    DataZ[2] = (TH1D *)gDirectory->Get("Z_40_80");

    TH1D *ratio[3];

    TCanvas *a = new TCanvas("Z", "Z", 1600,400);
	a->Divide(3);
	for (int i = 0; i < 3; i++){
		a->cd(i+1);
		gPad->SetLogy();
		
		TH1D *tmp = (TH1D *)RequiredZ->Clone("Z_Reweighted");
		tmp->Scale(RecoZ[i][0]->Integral()/tmp->Integral());
		
		MCZ[i][0]->Draw();
		MCZ[i][0]->GetYaxis()->SetRangeUser(pow(10, -9), pow(10, 0));

		if (power < 100) tmp->Draw("SAME");

		RecoZ[i][0]->Draw("SAME");
		MCZ[i][0]->SetLineColor(kBlue);
		MCZ[i][0]->SetMarkerColor(kBlue);
		MCZ[i][0]->SetMarkerSize(20);

		RecoZ[i][0]->SetLineColor(kBlack);
		RecoZ[i][0]->SetMarkerColor(kBlack);
		RecoZ[i][0]->SetMarkerSize(20);

		TH1D *tmp2 = (TH1D *)DataZ[i]->Clone("Data");
		tmp2->GetYaxis()->SetRangeUser(0,1);
		tmp2->Scale(RecoZ[i][0]->Integral()/tmp2->Integral());
		tmp2->Draw("SAME");
		tmp2->SetLineColor(kGreen-2);
		tmp2->SetMarkerColor(kGreen-2);

		// ratio[i] = (TH1D *)tmp2->Clone(Form("Ratio_%i", i));
		// ratio[i]->Divide(RecoZ[i][0]);

		ratio[i] = (TH1D *)RecoZ[i][0]->Clone(Form("Ratio_%i", i));
		ratio[i]->Divide(tmp2);
	}

	a->SaveAs(Form("RecoZReweight/Z_pow_%.1f.pdf", power));

	delete a;

	TCanvas *a0 = new TCanvas("ZRatio", "ZRatio", 400,400);
	a0->cd();
	for (int i = 0; i < 3; i++){
		gPad->SetLogy();
		ratio[i]->Draw("SAME");
		ratio[i]->SetLineColor(i+1);
		ratio[i]->SetMarkerColor(i+1);
	}

	a0->SaveAs(Form("RecoZReweight/RatioZDatavReco_pow_%.1f.pdf", power));

	delete a0;

	THnD *MCRecoD0PtJetPt[3][10];
	const int dimcount[7] = {0,2,3,4,5,6,7};
	for (int cent = 0; cent < 3; cent++){
		for (int i = 1; i <= 10; i++){
			THnD *tmp = (THnD *)Z[cent][0]->Clone();
			tmp->GetAxis(1)->SetRange(i, i);
			MCRecoD0PtJetPt[cent][i-1] = (THnD *)tmp->Projection(7, dimcount);
		}
	}

	for (int i = 0; i < 3; i++){
		ZResampled[i] = new THnF(Form("ZResampled_%i", i), Form("ZResampled_%i", i), nDimDaug, nBinsDaug, NULL, NULL);

	    ZResampled[i]->SetBinEdges(0, ZBins);
	    ZResampled[i]->SetBinEdges(1, ZBins);
	    ZResampled[i]->SetBinEdges(2, MCD0PtBins);
	    ZResampled[i]->SetBinEdges(3, MCJetPtBins);
	    ZResampled[i]->SetBinEdges(4, RecoD0PtBins);
	    ZResampled[i]->SetBinEdges(5, RecoJetPtBins);
	    ZResampled[i]->SetBinEdges(6, EtaBins);
    	ZResampled[i]->SetBinEdges(7, EtaBins);
	}

	TH1D *DataZRescaled[3];
	for (int i = 0; i < 3; i++){
		DataZRescaled[i] = (TH1D *)DataZ[i]->Clone();
		DataZRescaled[i]->Scale(RecoZ[i][0]->Integral()/DataZ[i]->Integral());
	}


	TRandom3 *r = new TRandom3(0);

	for (int cent = 0; cent < 3; cent++){
		for (int entries = 0; entries < numberofentries; entries++){		

			double randZ;
			if (power == 200.) randZ = DataZRescaled[2]->GetRandom(r);
			else if (power == 300.) randZ = DataZRescaled[cent]->GetRandom(r);
			else{
				randZ = RequiredZ->GetRandom(r);
			}
			int bin = RecoZ[cent][0]->FindBin(randZ);
			double val[7];
			MCRecoD0PtJetPt[cent][bin - 1]->GetRandom(val);

			double responsefiller[8];
			responsefiller[0] = val[0];
			responsefiller[1] = randZ;
			responsefiller[2] = val[1];
			responsefiller[3] = val[2];
			responsefiller[4] = val[3];
			responsefiller[5] = val[4];
			responsefiller[6] = val[5];
			responsefiller[7] = val[6];

			ZResampled[cent]->Fill(responsefiller);
		}
	}

	TFile *g = new TFile(Form("RecoZReweight/ResponsePlots_Weighted_pow_%.1f.root", power), "RECREATE");
	g->cd();
	for (int cent = 0; cent < 3; cent++){
		MCZ[cent][0]->Write();
		Z[cent][0]->Write();
		EfficiencyRatioAtGenLevel[cent]->Write();
		EfficiencyRatioAtRecoLevel[cent]->Write();
		ZResampled[cent]->Write();
	}
	g->Close();


	cout << "Resampled Histograms are filled " << endl;


	TCanvas *b = new TCanvas("Response", "Response", 2000, 600);
	b->Divide(3);
	for (int i = 0; i < 3; i++){
		b->cd(i+1);
		gPad->SetLogz(); 
		THnD *ntmp = (THnD *)ZResampled[i]->Clone();
		cout << "NTMP" << endl;
		// ntmp->Scale(numberofentries/ntmp->GetEntries());
		ntmp->GetAxis(6)->SetRange(2,2);
		ntmp->GetAxis(7)->SetRange(2,2);

		TH2D *tmp = (TH2D *)ntmp->Projection(3,5);
		cout << "TMP" << endl;
		tmp->Draw("COLZ");
		tmp->GetYaxis()->SetRangeUser(1,30);
		tmp->GetXaxis()->SetRangeUser(1,30);
		// tmp->GetZaxis()->SetRangeUser(pow(10, -10), pow(10, -4));

		cout << "Resampled Histograms are drawn " << endl;
	}

	b->SaveAs(Form("RecoZReweight/Response_pow_%.1f.pdf", power));

	// delete b;
}