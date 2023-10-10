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
#include "TRandom3.h"

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


void Method(int iteration){
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


	TCanvas *a = new TCanvas("a", "a", 800, 800);
	gPad->SetLogy();
	MCZ[0][0]->Draw();
	MCZ[0][0]->GetYaxis()->SetRangeUser(pow(10, -9), pow(10, 0));

	TCanvas *b0 = new TCanvas("b", "b", 800, 800);
	gPad->SetLogy();
	RecoZ[0][0]->Draw();
	RecoZ[0][0]->GetYaxis()->SetRangeUser(pow(10, -9), pow(10, 0));

	b0->SaveAs(Form("Closure/RZ_pow_%i_%i.pdf", iteration, 100));

	delete a;
	delete b0;

	THn  *ZResampled[3];

	int nBinsDaug[nDimDaug] = {nBinsMCZ, nBinsMCZ, nBinsMCD0Pt, nBinsMCJetPt, nBinsRecoD0Pt, nBinsRecoJetPt, nBinsEta, nBinsEta}; //MC Z, Reco Z, MC D0 Pt, MC Jet Pt, Reco D0 Pt, Reco Jet Pt

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

	delete gRandom;
	gRandom = new TRandom3(0);

	for (int cent = 0; cent < 3; cent++){
		for (int entries = 0; entries < numberofentries; entries++){		
			double v[7];
			Z[cent][0]->GetRandom(v);
			ZResampled[cent]->Fill(v);
		}
	}

	int nBinsZHist[nDimZHist] = {nBinsJetPtForZHist, nBinsForZHist, nBinsJetPtForZHist, nBinsForZHist};

	THn *JetPtZResampled[3];

	RooUnfoldResponse *responseJetPtvZ[3];
	RooUnfoldResponse *responseJetPt[3];

	TH2D *tmpMC = new TH2D("tmpMC", "tmpMC", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);
	TH2D *tmpReco = new TH2D("tmpReco", "tmpReco", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

	TH1D *tmpMC1D = new TH1D("tmpMC1D", "tmpMC1D", nBinsJetPtForZHist, JetPtBinsForZHist);
	TH1D *tmpReco1D = new TH1D("tmpReco1D", "tmpReco1D", nBinsJetPtForZHist, JetPtBinsForZHist);

	TH2D *MCJetPtvZ[3];
	TH2D *RecoJetPtvZ[3];

	for (int i = 0; i < 3; i++){
		JetPtZResampled[i] = new THnF(Form("JetPtZResampled_%i", i), Form("JetPtZResampled_%i", i), nDimZHist, nBinsZHist, NULL, NULL);
		MCJetPtvZ[i] = new TH2D(Form("MCJetPtvZ_%i", i), Form("MCJetPtvZ_%i", i), nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);
		RecoJetPtvZ[i] = new TH2D(Form("RecoJetPtvZ_%i", i), Form("RecoJetPtvZ_%i", i), nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

		JetPtZResampled[i]->SetBinEdges(0, JetPtBinsForZHist); //MCJetPtBins
	    JetPtZResampled[i]->SetBinEdges(1, ZBinsForZHist); //MCZBins
	    JetPtZResampled[i]->SetBinEdges(2, JetPtBinsForZHist); //RecoJetPtBins
	    JetPtZResampled[i]->SetBinEdges(3, ZBinsForZHist); //RecoZBins

		responseJetPtvZ[i] = new RooUnfoldResponse(tmpMC, tmpReco, Form("ResponseJetPtvZ_%i", i), Form("ResponseJetPtvZ_%i", i));
		responseJetPt[i] = new RooUnfoldResponse(tmpMC1D, tmpReco1D, Form("ResponseJetPt_%i", i), Form("ResponseJetPt_%i", i));

	}


	for (int cent = 0; cent < 3; cent++){
		for (int entries = 0; entries < numberofentries; entries++){		
			double v[4];
			JetPtZ[cent][0]->GetRandom(v);
			JetPtZResampled[cent]->Fill(v);
			MCJetPtvZ[cent]->Fill(v[0], v[1]);
			RecoJetPtvZ[cent]->Fill(v[2], v[3]);
			responseJetPtvZ[cent]->Fill(v[2], v[3], v[0], v[1]);
			responseJetPt[cent]->Fill(v[2], v[0]);
		}
	}

	TFile *g = new TFile(Form("Closure/ResponsePlots_%i.root", iteration), "RECREATE");
	g->cd();
	for (int cent = 0; cent < 3; cent++){
		MCZ[cent][0]->Write();
		Z[cent][0]->Write();
		EfficiencyRatioAtGenLevel[cent]->Write();
		EfficiencyRatioAtRecoLevel[cent]->Write();
		ZResampled[cent]->Write();
		JetPtZResampled[cent]->Write();
		MCJetPtvZ[cent]->Write();
		RecoJetPtvZ[cent]->Write();
		g->WriteObject(responseJetPt[cent], Form("ResponseJetPt_%i", cent));
		g->WriteObject(responseJetPtvZ[cent], Form("ResponseJetPtvZ_%i", cent));
	}
	g->Close();


	cout << "Resampled Histograms are filled " << endl;


	TCanvas *b = new TCanvas("Response", "Response", 2000, 600);
	b->Divide(3);
	for (int i = 0; i < 3; i++){
		b->cd(i+1);
		gPad->SetLogz(); 
		THnD *ntmp = (THnD *)ZResampled[i]->Clone(Form("ntmp_%i_%i", iteration, i));
		cout << "NTMP" << endl;
		// ntmp->Scale(numberofentries/ntmp->GetEntries());
		ntmp->GetAxis(6)->SetRange(2,2);
		ntmp->GetAxis(7)->SetRange(2,2);

		TH2D *tmp = (TH2D *)ntmp->Projection(3,5);
		cout << "TMP" << endl;
		tmp->Draw("COLZ");
		tmp->GetYaxis()->SetRangeUser(1, 30);
		tmp->GetXaxis()->SetRangeUser(1, 30);
		// tmp->GetZaxis()->SetRangeUser(pow(10, -10), pow(10, -4));

		cout << "Resampled Histograms are drawn " << endl;
	}

	b->SaveAs(Form("Closure/Response_pow_%i_%i.pdf", iteration, 100));

	delete b;
	gStyle->SetOptTitle(1);

	TCanvas *c = new TCanvas("Efficiency Factors", "Efficiency Factors", 2000, 800);
	c->cd();
	c->Divide(6,2);
	for (int cut = 0; cut < 12; cut++){
		c->cd(cut + 1);
		// TPaveText *pt = new TPaveText(0.60,0.25,0.80,0.35, "NDC"); // NDC sets coords
		if (cut < 11){
			for (int i = 0; i < 3; i++){
				TH1D *tmp = (TH1D *)MCD0JetPtAfterFiducialCuts[i][cut][0]->Clone(Form("tmp_%i_%i", iteration, i));
				tmp->SetNameTitle(Form("Cut %i", cut+1), Form("Cut %i", cut+1));
				tmp->Divide(MCD0JetPtBeforeFiducialCuts[i][cut][0]);
				tmp->Draw("SAME");
				tmp->SetLineColor(i+1);
				tmp->GetYaxis()->SetRangeUser(0, 1.2);
			}
		}

		else{
			for (int i = 0; i < 3; i++){
				TH1D *tmp = (TH1D *)RecoD0JetPtAfterFiducialCuts[i][0][0]->Clone(Form("tmp_%i_%i", iteration, i));
				tmp->SetNameTitle(Form("Cut %i", cut+1), Form("Cut %i", cut+1));
				tmp->Divide(RecoD0JetPtBeforeFiducialCuts[i][0][0]);
				tmp->Draw("SAME");
				tmp->SetLineColor(i+1);
				tmp->GetYaxis()->SetRangeUser(0, 1.2);
			}
		}
	}

	c->SaveAs(Form("Closure/AllEfficiencyPlotsQA_pow_%i_%i.pdf", iteration, 100));

	delete c;

	TCanvas *d = new TCanvas("d", "d", 2000, 1800);
	d->Divide(3,3);
	TH2D *tmp1, *tmp2, *tmp3;
	
	for (int i = 0; i < 3; i++){
		d->cd(i*3+1);
		gPad->SetLogz();

		THnD *n2tmp = (THnD *)JetPtZResampled[i]->Clone(Form("n2tmp_%i_%i", iteration, i));

		tmp1 = (TH2D *)n2tmp->Projection(1,0);
		tmp1->Draw("COLZ");
		// tmp1->GetZaxis()->SetRangeUser(pow(10, -9), pow(10,-2));

		d->cd(i*3+2);
		gPad->SetLogz();
		tmp2 = (TH2D *)n2tmp->Projection(3,2);
		tmp2->Draw("COLZ");
		// tmp2->GetZaxis()->SetRangeUser(pow(10, -9), pow(10,-2));

		d->cd(i*3+3);
		gPad->SetLogz();
		tmp3 = (TH2D *)n2tmp->Projection(3,1);
		tmp3->Draw("COLZ");
		// tmp3->GetZaxis()->SetRangeUser(pow(10, -9), pow(10,-2));

	}

	d->SaveAs(Form("Closure/JetPtvZResponse_pow_%i_%i.pdf", iteration, 100));

	delete d;

	for (int i = 0; i < 3; i++){
		delete ZResampled[i];
		delete JetPtZResampled[i];
		delete MCJetPtvZ[i];
		delete RecoJetPtvZ[i];
	}

	delete tmp1;
	delete tmp2;
	delete tmp3;
	
}

void ResponseMatrix(){
	Method(1);
	Method(2);
}