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
#include "TRandom3.h"
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

TString CutNames[11] = {  "1 < Reco D0 pT < 10", 
                          "1 < Reco Jet pT < 30", 
                          "Reco Jet |Eta| < 0.6", 
                          "0 < Reco Z < 1",
                          "1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6",
                          "1 < Reco D0 pT < 10 && 0 < Reco Z < 1 && Reco Jet |Eta| < 0.6",
                          "1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6 && 0 < Reco Z < 1",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6 && 1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6 && 1 < Reco D0 pT < 10 && 0 < Reco Z < 1 && Reco Jet |Eta| < 0.6",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6 && 1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6 && 0 < Reco Z < 1"
                        };


void IterativeResponseMatrix(int iteration = 1){
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	if (iteration == 1) {
		TFile *ZFile = new TFile("ZFile.root", "RECREATE");

		TH1D *ZForPrior = new TH1D(Form("ZForPrior_%i", iteration), Form("ZForPrior_%i", iteration), nBinsForZHist, ZBinsForZHist);

		for (int bin = ZForPrior->FindBin(0); bin <= ZForPrior->FindBin(1); bin++){
			ZForPrior->SetBinContent(bin, 1);
		}
		ZForPrior->Write();
		ZFile->Close();
	}

	double numberofentries = 10000000;
	// double numberofentries = 100;

	TFile *f = new TFile("Plots.root");
	TH1D *NEntries[2];
	TH1D *VzMC[3][2];
	TH1D *MCZ[3][2];
	TH1D *RecoZ[3][2];
	THn  *Z[3][2];
	THn  *JetPtZ[3][2]; //Wider range of Z and Jet Pt here.
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

		  JetPtZ[i][bin] = (THnD *)gDirectory->Get(Form("hMCJetPtMCZRecoJetPtRecoZ_%i", i));
		  
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

		  JetPtZ[i][bin]->Scale(scale);
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
		JetPtZ[i][0]->Add(JetPtZ[i][1]);
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

	int nBinsZHist[nDimZHist] = {nBinsJetPtForZHist, nBinsForZHist, nBinsJetPtForZHist, nBinsForZHist};

	THn *JetPtZResampled[3];

	RooUnfoldResponse *responseJetPtvZ[3];
	RooUnfoldResponse *responseJetPt[3];

	TH2D *tmpMC = new TH2D("tmpMC", "tmpMC", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);
	TH2D *tmpReco = new TH2D("tmpReco", "tmpReco", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

	TH1D *tmpMC1D = new TH1D("tmpMC1D", "tmpMC1D", nBinsJetPtForZHist, JetPtBinsForZHist);
	TH1D *tmpReco1D = new TH1D("tmpReco1D", "tmpReco1D", nBinsJetPtForZHist, JetPtBinsForZHist);

	for (int i = 0; i < 3; i++){
		JetPtZResampled[i] = new THnF(Form("JetPtZResampled_%i", i), Form("JetPtZResampled_%i", i), nDimZHist, nBinsZHist, NULL, NULL);

		JetPtZResampled[i]->SetBinEdges(0, JetPtBinsForZHist); //MCJetPtBins
	    JetPtZResampled[i]->SetBinEdges(1, ZBinsForZHist); //MCZBins
	    JetPtZResampled[i]->SetBinEdges(2, JetPtBinsForZHist); //RecoJetPtBins
	    JetPtZResampled[i]->SetBinEdges(3, ZBinsForZHist); //RecoZBins

		responseJetPtvZ[i] = new RooUnfoldResponse(tmpMC, tmpReco, Form("ResponseJetPtvZ_%i", i), Form("ResponseJetPtvZ_%i", i));
		responseJetPt[i] = new RooUnfoldResponse(tmpMC1D, tmpReco1D, Form("ResponseJetPt_%i", i), Form("ResponseJetPt_%i", i));

	}

	TRandom3 *r = new TRandom3(0);

	TFile *FileWithZPrior = new TFile("ZFile.root", "UPDATE");
	TH1D *PriorZ = (TH1D *)FileWithZPrior->Get(Form("ZForPrior_%i", iteration));

	THnD *MCJetPtRecoJetPtRecoZ[3][nBinsForZHist];
	const int dimcount[3] = {0,2,3};
	for (int cent = 0; cent < 3; cent++){
		for (int i = 1; i <= nBinsForZHist; i++){
			THnD *tmp = (THnD *)JetPtZ[cent][0]->Clone();
			tmp->GetAxis(1)->SetRange(i, i);
			MCJetPtRecoJetPtRecoZ[cent][i-1] = (THnD *)tmp->Projection(3, dimcount);
		}
	}

	for (int cent = 0; cent < 3; cent++){
		for (int entries = 0; entries < numberofentries; entries++){		
			double v[3];
			double randZ;
			randZ = PriorZ->GetRandom(r);

			int bin = PriorZ->FindBin(randZ);

			MCJetPtRecoJetPtRecoZ[cent][bin-1]->GetRandom(v);
			double vtofill[4] = {v[0], randZ, v[1], v[2]};
			JetPtZResampled[cent]->Fill(vtofill);
			responseJetPtvZ[cent]->Fill(vtofill[2], vtofill[3], vtofill[0], vtofill[1]);
			responseJetPt[cent]->Fill(vtofill[2], vtofill[0]);
		}
	}

	TFile *g2 = new TFile("../SPlotFrameWork/ApplyWeights/Histograms3_D01_10GeV_With2DJetPtZ.root");
    g2->cd();

    TH2D *JetPtvZ[3];
    JetPtvZ[0] = (TH2D *)gDirectory->Get("ZPt_0_10");
    JetPtvZ[1] = (TH2D *)gDirectory->Get("ZPt_10_40");
    JetPtvZ[2] = (TH2D *)gDirectory->Get("ZPt_40_80");

	RooUnfoldBayes unfold0 (responseJetPtvZ[0], JetPtvZ[0], 4);
    RooUnfoldBayes unfold1 (responseJetPtvZ[1], JetPtvZ[1], 4);
    RooUnfoldBayes unfold2 (responseJetPtvZ[2], JetPtvZ[2], 4);

    TH2D *UnfoldedJetPtvZ[3];

    UnfoldedJetPtvZ[0] = (TH2D *)unfold0.Hreco();
    UnfoldedJetPtvZ[1] = (TH2D *)unfold1.Hreco();
    UnfoldedJetPtvZ[2] = (TH2D *)unfold2.Hreco();

    FileWithZPrior->cd();

    TH1D *UnfoldedJetPt[3];
    TH1D *UnfoldedZ[3];
    for (int i = 0; i < 3; i++){
        UnfoldedJetPt[i] = (TH1D *)UnfoldedJetPtvZ[i]->ProjectionX();
        UnfoldedZ[i] = (TH1D *)UnfoldedJetPtvZ[i]->ProjectionY();
        UnfoldedZ[i]->Write();
    }

    TFile *fin1 = new TFile("../Files/D0Yield_Feb22.root");
    fin1->cd("D0Tagger");
    TH1F *hCent = (TH1F *)gDirectory->Get("hCentralityWeightedAfterCuts");

    for(int i =1;i<hCent->GetNbinsX()+1;i++)cout << hCent->GetBinContent(i) << endl;
    float nevts_0_10  = hCent->Integral(1, 2);
    float nevts_10_40 = hCent->Integral(3, 8);
    float nevts_40_80 = hCent->Integral(9, 16);

    double s[3];
    s[0] = 941.23714*nevts_0_10; //(1048*1.38188e+06+838*1.22039e+06)/(1.38188e+06+1.22039e+06);                                                                                                                              
    s[1] = 391.35550*nevts_10_40;//571 + 351 + 206;                                                                                                                                                                           
    s[2] = 56.62475*nevts_40_80;

    TCanvas *b = new TCanvas(Form("Unfolded_%i", iteration), Form("Unfolded_%.if", iteration), 800, 800);
    b->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    for (int i = 0; i < 3; i++){
        UnfoldedJetPt[i]->Scale(1./s[i]);
        UnfoldedJetPt[i]->Draw("P SAME");
        UnfoldedJetPt[i]->SetLineColor(1+i);
        UnfoldedJetPt[i]->SetMarkerColor(1+i);
        UnfoldedJetPt[i]->SetMarkerStyle(20);
        UnfoldedJetPt[i]->GetYaxis()->SetTitle("Jet Spectra");
        UnfoldedJetPt[i]->GetXaxis()->SetTitle("p_{T,jet} [GeV/#it{c}]");
   	}

    b->SaveAs(Form("IterativeRecoPtvZ/UnfoldedPtSpectra_%i.pdf", iteration));

    TCanvas *c = new TCanvas(Form("UnfoldedZ_%i", iteration), Form("UnfoldedZ_%i", iteration), 800, 800);
    c->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);

    TH1D *tmp[3];
    for (int i = 0; i < 3; i++){
    	tmp[i] = (TH1D *)UnfoldedZ[i]->Clone(Form("TmpZ_%i", i));
        tmp[i]->Scale(1./s[i]);
        tmp[i]->Draw("P SAME");
        tmp[i]->SetLineColor(1+i);
        tmp[i]->SetMarkerColor(1+i);
        tmp[i]->SetMarkerStyle(20);
        tmp[i]->GetYaxis()->SetTitle("Z Spectra");
        tmp[i]->GetXaxis()->SetTitle("Z");
   	}

    c->SaveAs(Form("IterativeRecoPtvZ/UnfoldedZSpectra_%i.pdf", iteration));

    FileWithZPrior->Close();
}