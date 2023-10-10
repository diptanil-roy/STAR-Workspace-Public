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
#include "TPaveText.h"
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

TH2D *ResponseFunc(THnD* Z){
	THnD *tmp = (THnD *)Z->Clone();
	tmp->GetAxis(6)->SetRange(2,2); //For response, everything is within |eta| < 0.6 (MC and Reco)
	tmp->GetAxis(7)->SetRange(2,2);
	// tmp->GetAxis(4)->SetRange(5,12);

	TH2D *h = (TH2D *)tmp->Projection(3,5);

	return h;
}

TH1D *MCFunc(THnD* Z){
	THnD *tmp = (THnD *)Z->Clone();
	tmp->GetAxis(6)->SetRange(2,2); //For truth, we only want all mc jets within |eta| < 0.6

	TH1D *h = (TH1D *)tmp->Projection(3);

	return h;
}

TH1D *RecoFunc(THnD* Z){
	THnD *tmp = (THnD *)Z->Clone();
	tmp->GetAxis(7)->SetRange(2,2);  //For detector, we only want all reco jets within |eta| < 0.6
	// tmp->GetAxis(4)->SetRange(5,12); //For detector, we only want all reco jets within 0 < pTCorr < 20 GeV/c

	TH1D *h = (TH1D *)tmp->Projection(5);

	return h;
}

TH1D *PreProcessFoldedHistogram(TH1D *h, TH1D *k = NULL){
    TH1D *R = (TH1D *)h->Clone();
    if (k) {R->Multiply(k); cout << "Corrected by reco efficiency" << endl;}

    return R;
}

TH1D *ProcessSpectraHistogram(TH1D *h, TH1D *k = NULL){
    TH1D *R = (TH1D *)h->Clone();
    for(int i = 1;i<h->GetNbinsX()+1;i++){
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        double center = fabs(R->GetBinCenter(i));

        R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.035);
        R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.035);
    }

    if (k) {R->Divide(k); cout << "Corrected by gen efficiency" << endl;}

    return R;
}

TH1D *Method(double method){
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    THnD *ZResampled[3];
    TH2D *Response[3];
    TH1D *MCJetPt[3];
    TH1D *RecoJetPt[3];

    TH1D *EfficiencyRatioForGen[3];
    TH1D *EfficiencyRatioForReco[3];

    TFile *g;

    if (method == 100)g = new TFile("RecoZReweight/ResponsePlots.root");
    else g = new TFile(Form("RecoPtvZ/ResponsePlots_Weighted_pow_%.1f.root", method));
	g->cd();

    RooUnfoldResponse *responseJetPtvZ[3];
    RooUnfoldResponse *responseJetPt[3];

	for (int i = 0; i < 3; i++){
        EfficiencyRatioForGen[i] = (TH1D *)gDirectory->Get(Form("EfficiencyRatioAtGenLevel_%i", i));
        EfficiencyRatioForReco[i] = (TH1D *)gDirectory->Get(Form("EfficiencyRatioAtRecoLevel_%i", i));
		ZResampled[i] = (THnD *)gDirectory->Get(Form("ZResampled_%i", i));
        g->GetObject(Form("ResponseJetPt_%i", i), responseJetPt[i]);
        g->GetObject(Form("ResponseJetPtvZ_%i", i), responseJetPtvZ[i]);
	}

    TFile *g2 = new TFile("../SPlotFrameWork/ApplyWeights/Histograms3_D01_10GeV_With2DJetPtZ.root");
    g2->cd();

    TH1D *JetPt1D[3];
    JetPt1D[0] = (TH1D *)gDirectory->Get("JetPt_0_10");
    JetPt1D[1] = (TH1D *)gDirectory->Get("JetPt_10_40");
    JetPt1D[2] = (TH1D *)gDirectory->Get("JetPt_40_80");

    TH2D *JetPtvZ[3];
    JetPtvZ[0] = (TH2D *)gDirectory->Get("ZPt_0_10");
    JetPtvZ[1] = (TH2D *)gDirectory->Get("ZPt_10_40");
    JetPtvZ[2] = (TH2D *)gDirectory->Get("ZPt_40_80");

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

    cout << s[0] << "\t" << s[1] << "\t" << s[2] << endl;


    TH1D *JetPt[3];
    TH1D *JetPtRecoCorrected[3];
    TCanvas *a = new TCanvas(Form("JetPtRecoEfficiencyCorrected_%.1f", method), Form("JetPtRecoEfficiencyCorrected_%.1f", method), 800, 800);
    a->cd();
    gPad->SetLogy();
    for (int i = 0; i < 3; i++){
        JetPt[i] = (TH1D *)JetPtvZ[i]->ProjectionX();
        JetPt[i]->Draw("SAME");
        // JetPtRecoCorrected[i] = (TH1D *)PreProcessFoldedHistogram(JetPt[i], EfficiencyRatioForReco[i]);
        JetPt[i]->SetLineColor(i+1);
        JetPt[i]->SetMarkerColor(i+1);
        JetPt[i]->SetMarkerStyle(20);

        // JetPtRecoCorrected[i]->SetLineColor(i+1);
        // JetPtRecoCorrected[i]->SetMarkerColor(i+1);
        // JetPtRecoCorrected[i]->SetMarkerStyle(24);
        // JetPtRecoCorrected[i]->Draw("SAME");
    }


    TH2D *UnfoldedJetPtvZ[3];
    TH1D *UnfoldedJetPt1D[3];

    RooUnfoldBayes unfold0 (responseJetPtvZ[0], JetPtvZ[0], 4);
    RooUnfoldBayes unfold1 (responseJetPtvZ[1], JetPtvZ[1], 4);
    RooUnfoldBayes unfold2 (responseJetPtvZ[2], JetPtvZ[2], 4);

    RooUnfoldBayes unfold1D0 (responseJetPt[0], JetPt1D[0], 4);
    RooUnfoldBayes unfold1D1 (responseJetPt[1], JetPt1D[1], 4);
    RooUnfoldBayes unfold1D2 (responseJetPt[2], JetPt1D[2], 4);

    UnfoldedJetPtvZ[0] = (TH2D *)unfold0.Hreco();
    UnfoldedJetPtvZ[1] = (TH2D *)unfold1.Hreco();
    UnfoldedJetPtvZ[2] = (TH2D *)unfold2.Hreco();

    UnfoldedJetPt1D[0] = (TH1D *)unfold1D0.Hreco();
    UnfoldedJetPt1D[1] = (TH1D *)unfold1D1.Hreco();
    UnfoldedJetPt1D[2] = (TH1D *)unfold1D2.Hreco();

    TH1D *UnfoldedJetPt[3];
    TH1D *UnfoldedZ[3];
    for (int i = 0; i < 3; i++){
        UnfoldedJetPt[i] = (TH1D *)UnfoldedJetPtvZ[i]->ProjectionX();
        UnfoldedZ[i] = (TH1D *)UnfoldedJetPtvZ[i]->ProjectionY();
    }


    TPaveText *pt = new TPaveText(0.20,0.65,0.30,0.85, "NDC"); // NDC sets coords
                                       // relative to pad dimensions
    pt->SetFillStyle(4000);
    pt->SetTextSize(0.04); 
    pt->SetTextAlign(12);

    if (method == 100)auto pt_text0 = pt->AddText("PYTHIA Z");
    else if (method == 200)auto pt_text0 = pt->AddText("Peripheral Z (unc.)");
    else if (method == 300)auto pt_text0 = pt->AddText("Centrality Z (unc.)");
    else auto pt_text0 = pt->AddText(Form("#frac{dN}{dz} ~ (1-z)^{%.1f}", method));

    TCanvas *b = new TCanvas(Form("Unfolded_%.1f", method), Form("Unfolded_%.1f", method), 800, 800);
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

        UnfoldedJetPt1D[i]->Scale(1./s[i]);
        UnfoldedJetPt1D[i]->Draw("P SAME");
        UnfoldedJetPt1D[i]->SetLineColor(1+i);
        UnfoldedJetPt1D[i]->SetMarkerColor(1+i);
        UnfoldedJetPt1D[i]->SetMarkerStyle(24);
        // cout << "Resampled Histograms are drawn " << endl;
    }
    pt->Draw("SAME");

    b->SaveAs(Form("RecoPtvZ/UnfoldedSpectra_%.1f.pdf", method));

    // delete b;

    TH1D *Ratio[3];
    TH1D *Ratio1D[3];
    Ratio[0] = (TH1D *)UnfoldedJetPt[0]->Clone();
    Ratio[1] = (TH1D *)UnfoldedJetPt[1]->Clone();
    Ratio[2] = (TH1D *)UnfoldedJetPt[0]->Clone();

    Ratio1D[0] = (TH1D *)UnfoldedJetPt1D[0]->Clone();
    Ratio1D[1] = (TH1D *)UnfoldedJetPt1D[1]->Clone();
    Ratio1D[2] = (TH1D *)UnfoldedJetPt1D[0]->Clone();

    Ratio[0]->Divide(UnfoldedJetPt[2]);
    Ratio[1]->Divide(UnfoldedJetPt[2]);
    Ratio[2]->Divide(UnfoldedJetPt[1]);

    Ratio1D[0]->Divide(UnfoldedJetPt1D[2]);
    Ratio1D[1]->Divide(UnfoldedJetPt1D[2]);
    Ratio1D[2]->Divide(UnfoldedJetPt1D[1]);

    TCanvas *c = new TCanvas(Form("Ratio_%.1f", method), Form("Ratio_%.1f", method), 800, 800);
    c->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    for (int i = 0; i < 2; i++){        
        Ratio[i]->GetYaxis()->SetRangeUser(pow(10, -1), pow(10,2));
        Ratio[i]->Draw("P SAME");
        Ratio1D[i]->Draw("P SAME");
        Ratio[i]->GetYaxis()->SetTitle("R_{CP}");
        Ratio[i]->GetXaxis()->SetTitle("p_{T,jet} [GeV/#it{c}]");
        // cout << "Resampled Histograms are drawn " << endl;
    }
    pt->Draw("SAME");

    c->SaveAs(Form("RecoPtvZ/RCP_%.1f.pdf", method));

    TCanvas *d = new TCanvas(Form("UnfoldedZ_%.1f", method), Form("UnfoldedZ_%.1f", method), 800, 800);
    d->cd();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    for (int i = 0; i < 3; i++){
        UnfoldedZ[i]->Scale(1./s[i]);
        UnfoldedZ[i]->Draw("P SAME");
        UnfoldedZ[i]->SetLineColor(1+i);
        UnfoldedZ[i]->GetYaxis()->SetTitle("Z Spectra");
        UnfoldedZ[i]->GetXaxis()->SetTitle("Z");
        // cout << "Resampled Histograms are drawn " << endl;
    }
    pt->Draw("SAME");

    d->SaveAs(Form("RecoPtvZ/UnfoldedZ_%.1f.pdf", method));

    TH1D *h = (TH1D *)Ratio[0]->Clone();
    return h;
}


void JetPtvsZ(){

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

	TH1D *R1[13];
	TH1D *R2[13];
	// int col[13] = {1, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84};
	// double methodname[13] = {100.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, -0.5, -0.8, -1.0, -1.2, -1.5, -2.0};
	// double methodname[13] = {100.0, -2.0, -1.5, -1.2, -1.0, -0.8, -0.5, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
	int col[13] = {1, 51, 84};
	double methodname[3] = {100.0};

	for (int i = 0; i < 1; i++){
		R1[i] = (TH1D *)Method(methodname[i]);
	}

	TCanvas *c 	= new TCanvas("RatioOfRatios", "RatioOfRatios", 800, 800);
	TH1D *DoubleRatioOfRCP[13];
	c->cd();
	gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
	for (int i = 0; i < 1; i++){
		// R1[i]->Draw("SAME");
		DoubleRatioOfRCP[i] = (TH1D *)R1[i]->Clone();
		DoubleRatioOfRCP[i]->Divide(R1[0]);
		DoubleRatioOfRCP[i]->SetLineColor(col[i]);
		DoubleRatioOfRCP[i]->SetMarkerColor(col[i]);
		DoubleRatioOfRCP[i]->SetMarkerStyle(20);
		DoubleRatioOfRCP[i]->SetMarkerSize(1.2);

        if (i == 0){
            for (int bin = 1; bin < DoubleRatioOfRCP[i]->GetNbinsX(); bin++){
                DoubleRatioOfRCP[i]->SetBinError(bin, 0);
                // if (bin == DoubleRatioOfRCP[i]->GetNbinsX()) DoubleRatioOfRCP[i]->SetBinContent(bin, 0);
            }
        }

		DoubleRatioOfRCP[i]->GetYaxis()->SetTitle("Ratio of RCP");
		DoubleRatioOfRCP[i]->GetYaxis()->SetRangeUser(-1, 5);

		if (i==0)DoubleRatioOfRCP[i]->SetNameTitle("PYTHIA", "PYTHIA");
		// else if (i==1)DoubleRatioOfRCP[i]->SetNameTitle("Peripheral Reweight", "Peripheral Reweight");
		else if (i==1)DoubleRatioOfRCP[i]->SetNameTitle("Cent Reweight", "Cent Reweight");
		// else DoubleRatioOfRCP[i]->SetNameTitle(Form("~(1-z)^%.1f", methodname[i]), Form("~(1-z)^%.1f", methodname[i]));
		DoubleRatioOfRCP[i]->Draw("SAME"); 
	}

	gPad->BuildLegend(0.55,0.65,0.85,0.85);

    c->SaveAs(Form("RecoPtvZ/RatioofRCPs.pdf"));
}
