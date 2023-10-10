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

TH1D *ProcessSpectraHistogram(TH1D *h){
    TH1D *R = (TH1D *)h->Clone();
    for(int i = 1;i<h->GetNbinsX()+1;i++){
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        double center = fabs(R->GetBinCenter(i));

        R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.035);
        R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.035);
    }

    return R;
}

const int NEntries = 1000000;

void Z1(TH1D *h, double smearfactor, TH2D *hResponse, TH1D *hMC, TH1D *hReco, TH1D *hZ, TH1D *hMCD0Pt){
	TRandom3 *r = new TRandom3(0);
	for (int ziter = 0; ziter < 10; ziter++){
		int i = 0;
		while(i < NEntries){
			double d0pt = h->GetRandom();
			double zval = r->Uniform(ziter*0.1, (ziter+1)*0.1);
			double mcjetpt = d0pt/zval;
			if (mcjetpt < 1 || mcjetpt > 20) continue;
			i++;
			double recojetpt = mcjetpt + r->Gaus(0, smearfactor);
			if (recojetpt > 3 && recojetpt < 20 && mcjetpt > 3 && mcjetpt < 20) hResponse->Fill(recojetpt, mcjetpt);
			if (mcjetpt > 3 && mcjetpt < 20) hMC->Fill(mcjetpt);
			if (recojetpt > 3 && recojetpt < 20) hReco->Fill(recojetpt);
			hZ->Fill(zval);
			hMCD0Pt->Fill(d0pt);
		}
	}
}

void Z2(TH1D *h, double smearfactor, TH2D *hResponse, TH1D *hMC, TH1D *hReco, TH1D *hZ, TH1D *hMCD0Pt){
	TRandom3 *r = new TRandom3(0);
	for (int ziter = 0; ziter < 10; ziter++){
		int i = 0;
		while(i < (ziter + 1)*NEntries){
			double d0pt = h->GetRandom();
			double zval = r->Uniform(ziter*0.1, (ziter+1)*0.1);
			double mcjetpt = d0pt/zval;
			if (mcjetpt < 1 || mcjetpt > 20) continue;
			i++;
			double recojetpt = mcjetpt + r->Gaus(0, smearfactor);
			if (recojetpt > 3 && recojetpt < 20 && mcjetpt > 3 && mcjetpt < 20) hResponse->Fill(recojetpt, mcjetpt);
			if (mcjetpt > 3 && mcjetpt < 20) hMC->Fill(mcjetpt);
			if (recojetpt > 3 && recojetpt < 20) hReco->Fill(recojetpt);
			hZ->Fill(zval);
			hMCD0Pt->Fill(d0pt);
		}
	}
}

void Z3(TH1D *h, double smearfactor, TH2D *hResponse, TH1D *hMC, TH1D *hReco, TH1D *hZ, TH1D *hMCD0Pt){
	TRandom3 *r = new TRandom3(0);
	for (int ziter = 0; ziter < 10; ziter++){
		int i = 0;
		while(i < NEntries/(ziter + 1)){
			double d0pt = h->GetRandom();
			double zval = r->Uniform(ziter*0.1, (ziter+1)*0.1);
			double mcjetpt = d0pt/zval;
			if (mcjetpt < 1 || mcjetpt > 20) continue;
			i++;
			double recojetpt = mcjetpt + r->Gaus(0, smearfactor);
			if (recojetpt > 3 && recojetpt < 20 && mcjetpt > 3 && mcjetpt < 20) hResponse->Fill(recojetpt, mcjetpt);
			if (mcjetpt > 3 && mcjetpt < 20) hMC->Fill(mcjetpt);
			if (recojetpt > 3 && recojetpt < 20) hReco->Fill(recojetpt);
			hZ->Fill(zval);
			hMCD0Pt->Fill(d0pt);
		}
	}
}

TH1D *Method(int method){

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	gSystem->Load("RooUnfold/libRooUnfold");    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

	TFile *f = new TFile("Plots.root");
	TH1D *NEntries[2];
	TH1D *VzMC[3][2];
	TH1D *D0Pt[3][2];

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

		  D0Pt[i][bin] = (TH1D *)gDirectory->Get(Form("MCD0Pt_%i", i));
		}
	}
	D0Pt[0][0]->Add(D0Pt[0][1]);
	TH1D *h = (TH1D *)D0Pt[0][0]->Clone();

	TH2D *Response[3];
	TH1D *MC[3];
	TH1D *Reco[3];
	TH1D *Z[3];
	TH1D *MCD0Pt[3];

	const int nBinsMCJetPt = 8;
	double MCJetPtBins[nBinsMCJetPt + 1]         = {1, 3, 5, 7, 9, 11, 13, 15, 20};
	const int nBinsRecoJetPt = 8;
	double RecoJetPtBins[nBinsRecoJetPt + 1]     = {1, 3, 5, 7, 9, 11, 13, 15, 20};

	double smearfactors[3] = {6.69, 4.60, 2.16};

	for (int i = 0; i < 3; i++){
		Response[i] = new TH2D(Form("Response_%i", i), Form("Response_%i", i), nBinsRecoJetPt, RecoJetPtBins, nBinsMCJetPt, MCJetPtBins);
		MC[i] = new TH1D(Form("MC_%i", i), Form("MC_%i", i), nBinsMCJetPt, MCJetPtBins);
		Reco[i] = new TH1D(Form("Reco_%i", i), Form("Reco_%i", i), nBinsRecoJetPt, RecoJetPtBins);
		Z[i] = new TH1D(Form("SampleRecoZ_%i", i), Form("SampleRecoZ_%i", i), 10, 0, 1);
		MCD0Pt[i] = new TH1D(Form("SampleMCD0Pt_%i", i), Form("SampleMCD0Pt_%i", i), 20, 0., 10.);

		if (method == 1)Z1(h, smearfactors[i], Response[i], MC[i], Reco[i], Z[i], MCD0Pt[i]);
		else if (method == 2)Z2(h, smearfactors[i], Response[i], MC[i], Reco[i], Z[i], MCD0Pt[i]);
		else if (method == 3)Z3(h, smearfactors[i], Response[i], MC[i], Reco[i], Z[i], MCD0Pt[i]);
	}

	RooUnfoldResponse *response[3];

    for (int i = 0; i < 3; i++){
		response[i] = new RooUnfoldResponse(Reco[i],MC[i],Response[i],Form("Response_%i", i),Form("Response_%i", i));
    }

	TFile *g2 = new TFile("../SPlotFrameWork/ApplyWeights/Histograms3_D01_10GeV_ForTest.root");
    g2->cd();

    TH1D *JetPt[3];
    JetPt[0] = (TH1D *)gDirectory->Get("JetPt_0_10");
    JetPt[1] = (TH1D *)gDirectory->Get("JetPt_10_40");
    JetPt[2] = (TH1D *)gDirectory->Get("JetPt_40_80");

	TH1D *UnfoldedJetPt[3];

	RooUnfoldBayes unfold0 (response[0], JetPt[0], 4);
    RooUnfoldBayes unfold1 (response[1], JetPt[1], 4);
    RooUnfoldBayes unfold2 (response[2], JetPt[2], 4);

    UnfoldedJetPt[0] = (TH1D *)ProcessSpectraHistogram((TH1D*) unfold0.Hreco());
    UnfoldedJetPt[1] = (TH1D *)ProcessSpectraHistogram((TH1D*) unfold1.Hreco());
    UnfoldedJetPt[2] = (TH1D *)ProcessSpectraHistogram((TH1D*) unfold2.Hreco());

    TFile *fin1 = new TFile("../Files/D0Yield_Feb22.root");
    fin1->cd("D0Tagger");
    TH1F *hCent = (TH1F *)gDirectory->Get("hCentralityWeightedAfterCuts");

    for(int i =1;i<hCent->GetNbinsX()+1;i++)cout << hCent->GetBinContent(i) << endl;
    float nevts_0_10 = hCent->Integral(1, 2);
    float nevts_10_40 = hCent->Integral(3, 8);
    float nevts_40_80 = hCent->Integral(9, 16);

    double s[3];
    s[0] = 941.23714*nevts_0_10; //(1048*1.38188e+06+838*1.22039e+06)/(1.38188e+06+1.22039e+06);                                                                                                                              
    s[1] = 391.35550*nevts_10_40;//571 + 351 + 206;                                                                                                                                                                           
    s[2] = 56.62475*nevts_40_80;

    TCanvas *a = new TCanvas(Form("Response_%i", method), Form("Response_%i", method), 1600, 800);
	a->Divide(3);
	for (int i = 0; i < 3; i++){
		a->cd(i+1);
		gPad->SetLogz(); 
		Response[i]->Draw("COLZ");
	}

	TCanvas *a0 = new TCanvas(Form("D0Pt_%i", method), Form("D0Pt_%i", method), 1600, 800);
	a0->Divide(3);
	for (int i = 0; i < 3; i++){
		a0->cd(i+1);
		gPad->SetLogy(); 
		MCD0Pt[i]->Draw();
	}

	TCanvas *a1 = new TCanvas(Form("Z_%i", method), Form("Z_%i", method), 1600, 800);
	a1->Divide(3);
	for (int i = 0; i < 3; i++){
		a1->cd(i+1);
		gPad->SetLogy(); 
		Z[i]->Draw();
	}

	TCanvas *a2 = new TCanvas(Form("JetPt_%i", method), Form("JetPt_%i", method), 1600, 800);
	a2->Divide(3);
	for (int i = 0; i < 3; i++){
		a2->cd(i+1);
		gPad->SetLogy(); 
		TH1D *tmp = (TH1D *)ProcessSpectraHistogram(MC[i]);
		tmp->Draw();
		UnfoldedJetPt[i]->Draw("P SAME");
	}


    TCanvas *b = new TCanvas(Form("Unfolded_%i", method), Form("Unfolded_%i", method), 1600, 800);
    b->cd();
    for (int i = 0; i < 3; i++){
        gPad->SetLogy();
        UnfoldedJetPt[i]->Scale(1./s[i]);
        UnfoldedJetPt[i]->Draw("P SAME");
        UnfoldedJetPt[i]->SetLineColor(1+i);
        UnfoldedJetPt[i]->GetYaxis()->SetTitle("Jet Spectra");
        UnfoldedJetPt[i]->GetXaxis()->SetTitle("p_{T,jet} [GeV/#it{c}]");
        // cout << "Resampled Histograms are drawn " << endl;
    }

    TH1D *Ratio[3];
    Ratio[0] = (TH1D *)UnfoldedJetPt[0]->Clone();
    Ratio[1] = (TH1D *)UnfoldedJetPt[1]->Clone();
    Ratio[2] = (TH1D *)UnfoldedJetPt[0]->Clone();

    Ratio[0]->Divide(UnfoldedJetPt[2]);
    Ratio[1]->Divide(UnfoldedJetPt[2]);
    Ratio[2]->Divide(UnfoldedJetPt[1]);

    TCanvas *c = new TCanvas(Form("Ratio_%i", method), Form("Ratio_%i", method), 1600, 800);
    c->cd();
    for (int i = 0; i < 2; i++){
        gPad->SetLogy();
        Ratio[i]->Draw("P SAME");
        Ratio[i]->GetYaxis()->SetTitle("R_{CP}");
        Ratio[i]->GetXaxis()->SetTitle("p_{T,jet} [GeV/#it{c}]");
        // Ratio[i]->SetLineColor(1+i);
        // cout << "Resampled Histograms are drawn " << endl;
    }

    TH1D *k = (TH1D *)Ratio[0]->Clone();
    return k;

}


void Test(){
	TH1D *R[3];
	for (int i = 0; i < 3; i++){
		R[i] = (TH1D *)Method(i+1);
	}

	int colors[3] = {kBlack, kRed, kGreen - 2};
	TString Name[3] = {"Flat", "~z", "~1/z"};
	TCanvas *a = new TCanvas("RCP For Different Z", "RCP For Different Z", 800, 800);
	a->cd();
	gPad->SetLogy();
	for (int i = 0; i < 3; i++){
		R[i]->SetLineColor(colors[i]);
		R[i]->SetMarkerColor(colors[i]);
		R[i]->SetMarkerStyle(20);
		R[i]->SetNameTitle(Name[i].Data(), Name[i].Data());
		R[i]->Draw("SAME");
	}

	gPad->BuildLegend();


}