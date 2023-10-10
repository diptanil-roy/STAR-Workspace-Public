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

void PlotterForComparingSingleParticleWithHIOverlay(){

	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
  	gStyle->SetOptStat(0);

	const int filenum = 7;

	TH1D *dPt[3][filenum];
	TH1D *dEta[3][filenum];
	TH1D *dPhi[3][filenum];

	TH2D *Unfolded2D[filenum][3];
	TH1D *UnfoldedPt2D[filenum][3];
	TH1D *UnfoldedZ2D[filenum][3];
	TH1D *UnfoldedPt1D[filenum][3];

	TH1D *TruthPt1D[filenum][3];

	TH1D *Truth;

	TH1D *RCPPt[filenum];
	TH1D *RCPPt1D[filenum];
	TH1D *RCPZ[filenum];

	TString RFile[filenum];


	RFile[0] = "Dec26_SinglePart/SingleParticleEmbedding.root";
	RFile[1] = "Jan26_FONLL_MC_dPtvPt4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	RFile[2] = "Jan26_FONLL_HI4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	// RFile[3] = "Jan26_FONLL_MC_CombinedFile4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	RFile[3] = "Mar2_FONLL_MC_Final_1D4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	// RFile[4] = "Jan26_FONLL_HI_CombinedFile4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	RFile[4] = "Feb15_FONLL_HI_CombinedFile4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	RFile[5] = "Feb15_FONLL_MC_Checking2D4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_8.root";
	// RFile[6] = "Feb15_FONLL_HI_Checking2D4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_8.root";
	// RFile[6] = "Feb20_FONLL_HI_WiderpTBins4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_7.root";
	RFile[6] = "Mar2_FONLL_HI_Final_2D4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	// RFile[7] = "Dec9_DataHI_RecoMaxTrackPtCut1MaxTowerPtCut13/OldResponse_Step_1_IterParam_3_njptbin_10_nzbin_10.root";
	// // RFile[8] = "Dec9_DataMC_PtEtaPhiChangeFromSP3/OldResponse_Step_1_IterParam_3_njptbin_10_nzbin_10.root";
	// RFile[8] = "Jan8_FONLL_MC4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	// RFile[9] = "Dec9_DataMC_EtaPhiChangeFromHI3/OldResponse_Step_1_IterParam_3_njptbin_10_nzbin_10.root";
	// RFile[10] = "Dec13_DataMC_PtChangeFromSP3/OldResponse_Step_1_IterParam_3_njptbin_10_nzbin_10.root";
	// RFile[11] = "Dec9_DataMC_PtEtaPhiChangeFromHI3/OldResponse_Step_1_IterParam_3_njptbin_10_nzbin_10.root";
	// RFile[12] = "Dec29_SmearFromMC_DataHIUnf_Updated4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";
	// RFile[13] = "Jan5_DataHIUnf_NewHIResponse_PtCutOff34/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root";

	TString OFile[filenum];

	for (int fi = 1; fi < filenum; fi++){
		OFile[fi] = RFile[fi];
		OFile[fi].ReplaceAll("OldResponse", "OldOutput");
		// cout << OFile[fi].Data() << endl;
	}

	TString LegAppend[filenum];
	LegAppend[0] = "SP Width";
	LegAppend[1] = "MC (FONLL), Iter = 4";
	LegAppend[2] = "HI (FONLL), Iter = 4";
	LegAppend[3] = "MC (FONLL, Iter = 4)";
	LegAppend[4] = "HI (FONLL, Iter = 4)";
	LegAppend[5] = "MC (FONLL), Iter = 4 (Checking 2D)";
	LegAppend[6] = "HI (FONLL, Iter = 4)";
	// LegAppend[7] = "HI Width (pT, ET Max = 1)";
	// LegAppend[8] = "MC Width (dp_{T} from SP, d#eta, d#phi from HI)";
	// LegAppend[9] = "MC Width (d#eta, d#phi from HI)";
	// LegAppend[10] = "MC Width (dp_{T} from SP)";
	// LegAppend[11] = "MC Width (dp_{T}, d#eta, d#phi from HI)";
	// LegAppend[12] = "HI Width Smeared from MC";
	// LegAppend[13] = "New HI (pT, ET Max = 3)";

	for (int fi = 0; fi < filenum; fi++){
		// cout << RFile[fi].Data() << endl;
		TFile *f = new TFile(RFile[fi].Data());
		f->cd();
		for (int i = 0; i < 3; i++){
			dPt[i][fi] = (TH1D *)gDirectory->Get(Form("hDiffJetPt_%i", i));
			dEta[i][fi] = (TH1D *)gDirectory->Get(Form("hDiffJetEta_%i", i));
			dPhi[i][fi] = (TH1D *)gDirectory->Get(Form("hDiffJetPhi_%i", i));

			dPt[i][fi]->SetDirectory(0);
			dEta[i][fi]->SetDirectory(0);
			dPhi[i][fi]->SetDirectory(0);
		}
		f->Close();
	}

	// cout << "Imported All the Plots" << endl;

	double naa[3] = {1.0318440e+08, 3.2123506e+08, 4.6679240e+08};
	double taa[3] = {941.23714, 391.35550, 56.62475};
	
	TH1D *tmp[3];

	for (int fi = 1; fi < filenum; fi++){
		TFile *f = new TFile(OFile[fi].Data());
		f->cd();

		for (int i = 0; i < 3; i++){

			UnfoldedPt1D[fi][i] = (TH1D *)gDirectory->Get(Form("Unfolded1D_%i", i));
			UnfoldedPt1D[fi][i]->SetDirectory(0);
			UnfoldedPt1D[fi][i]->Scale(1./naa[i]);
			// cout << UnfoldedPt1D[fi][i]->Integral() << endl;

			Unfolded2D[fi][i] = (TH2D *)gDirectory->Get(Form("Unfolded_%i", i));
			UnfoldedPt2D[fi][i] = (TH1D *)Unfolded2D[fi][i]->ProjectionX();
			UnfoldedZ2D[fi][i] = (TH1D *)Unfolded2D[fi][i]->ProjectionY(Form("UnfoldedZ2D_%i", i), UnfoldedPt2D[fi][i]->FindBin(7.001), UnfoldedPt2D[fi][i]->FindBin(19.999));
			Unfolded2D[fi][i]->SetDirectory(0);
			UnfoldedPt2D[fi][i]->SetDirectory(0);
			UnfoldedZ2D[fi][i]->SetDirectory(0);
			UnfoldedPt2D[fi][i]->Scale(1./naa[i]);
			UnfoldedZ2D[fi][i]->Scale(1./naa[i]);

			UnfoldedPt2D[fi][i] = (TH1D *)UnfoldedPt2D[fi][i]->Rebin(nbins_jpt, Form("Unfolded2DRebinned_%i", i), binning_jpt);
			UnfoldedPt1D[fi][i] = (TH1D *)UnfoldedPt1D[fi][i]->Rebin(nbins_jpt, Form("Unfolded1DRebinned_%i", i), binning_jpt);

			UnfoldedPt1D[fi][i] = (TH1D *)ProcessSpectraHistogram(UnfoldedPt1D[fi][i]);
			UnfoldedPt1D[fi][i]->SetDirectory(0);

			UnfoldedPt2D[fi][i] = (TH1D *)ProcessSpectraHistogram(UnfoldedPt2D[fi][i]);
			UnfoldedPt2D[fi][i]->SetDirectory(0);

			UnfoldedZ2D[fi][i] = (TH1D *)ProcessSpectraHistogram(UnfoldedZ2D[fi][i]);
			UnfoldedZ2D[fi][i]->SetDirectory(0);

			// if (fi == 6){
			// 	TruthPt1D[fi][i] = (TH1D *)gDirectory->Get(Form("fTrue1D_Cent_%i", i));
			// 	cout << TruthPt1D[fi][i]->Integral() << endl;
			// 	TruthPt1D[fi][i]->SetDirectory(0);
			// 	TruthPt1D[fi][i]->Scale(1./TruthPt1D[fi][i]->Integral());
			// 	TruthPt1D[fi][i]->SetDirectory(0);
			// 	TruthPt1D[fi][i] = (TH1D *)ProcessSpectraHistogram(TruthPt1D[fi][i]);
			// 	TruthPt1D[fi][i]->SetDirectory(0);
			// }
		}

		// cout << UnfoldedPt2D[fi][0]->GetName() << "\t" << UnfoldedPt2D[fi][0]->Integral() << endl;

		RCPPt[fi] = (TH1D *)UnfoldedPt2D[fi][0]->Clone(Form("RCP_pT_1_0"));
		RCPPt[fi]->Divide(UnfoldedPt2D[fi][2]);
		RCPPt[fi]->Scale(taa[2]/taa[0]);
		
		// RCPPt1D[fi] = (TH1D *)gDirectory->Get(Form("RCP1D_pT_1_0"));
		RCPZ[fi] = (TH1D *)gDirectory->Get(Form("RCP_Z_1_0"));

		RCPZ[fi] = (TH1D *)UnfoldedZ2D[fi][0]->Clone(Form("RCP_Z_1_0"));
		RCPZ[fi]->Divide(UnfoldedZ2D[fi][2]);
		RCPZ[fi]->Scale(taa[2]/taa[0]);


		RCPPt1D[fi] = (TH1D *)UnfoldedPt1D[fi][0]->Clone(Form("RCP1D_pT_1_0"));
		RCPPt1D[fi]->Divide(UnfoldedPt1D[fi][2]);
		RCPPt1D[fi]->Scale(taa[2]/taa[0]);

		RCPPt[fi]->SetDirectory(0);
		RCPPt1D[fi]->SetDirectory(0);
		RCPZ[fi]->SetDirectory(0);
		f->Close();
	}

	for (int fi = 1; fi < filenum; fi++){
		TFile *f = new TFile(RFile[fi].Data());
		f->cd();

		for (int i = 0; i < 3; i++){
			if (fi == 6){
				TruthPt1D[fi][i] = (TH1D *)gDirectory->Get(Form("fTrue1D_Cent_%i", i));
				cout << TruthPt1D[fi][i]->Integral() << endl;
				TruthPt1D[fi][i]->SetDirectory(0);
				TruthPt1D[fi][i]->Scale(1./TruthPt1D[fi][i]->Integral(TruthPt1D[fi][i]->FindBin(5.01), TruthPt1D[fi][i]->FindBin(19.99)));
				TruthPt1D[fi][i]->SetDirectory(0);
				TruthPt1D[fi][i] = (TH1D *)ProcessFONLLHistogram(TruthPt1D[fi][i]);
				TruthPt1D[fi][i]->SetDirectory(0);
			}
		}
		f->Close();
	}

	Truth = (TH1D *)TruthPt1D[6][0]->Clone();
	SetName(Truth, "FONLL");

	Truth->SetBinContent(1, 5.2420e+05);
	Truth->SetBinContent(2, 9.2870e+04);
	Truth->SetBinContent(3, 2.1330e+04);
	Truth->SetBinContent(4, 5.8870e+03);
	Truth->SetBinContent(5, 1.8600e+03);
	Truth->SetBinContent(6, 9.6000e+02);

	Truth->Scale(1./6.4710e+05);
	Truth = (TH1D *)ProcessFONLLHistogram(Truth);
	cout << "Unfolded Pt 1D = " << UnfoldedPt1D[6][2]->Integral() << endl;

	for (int i = 0; i < 3; i++){
		TruthPt1D[6][i] = (TH1D *)Truth->Clone(Form("FONLL_%i", i));
		TruthPt1D[6][i]->Scale(UnfoldedPt1D[6][i]->Integral()/Truth->Integral());
	}

	// Truth->Scale(UnfoldedPt1D[6][2]->Integral()/Truth->Integral());

	// TH1D *RCPPt1DForNumber8Rebinned = (TH1D *)UnfoldedPt1D[8][0]->Clone("RCPPt1DForNumber8Rebinned");
	// RCPPt1DForNumber8Rebinned->Divide(UnfoldedPt1D[8][2]);
	// // RCPPt1DForNumber8Rebinned->SetLineColor(kGreen-2);
	// // RCPPt1DForNumber8Rebinned->SetMarkerColor(kGreen-2);
	// RCPPt1D[8] = (TH1D *)RCPPt1DForNumber8Rebinned->Clone();


	TFile *SpectraFile = new TFile("/Volumes/WorkDrive/work/2022/PreliminaryPlots/April2/JetRCP.root");
	SpectraFile->cd();
	TH1D *Central    = (TH1D *)gDirectory->Get("JetRCP_0_10");
	Central->SetFillColorAlpha(kBlack, 0);
	Central->SetLineColor(kBlack);
	Central->SetMarkerColor(kBlack);
	Central->SetMarkerStyle(29);
	Central->SetMarkerSize(1.5);

	TFile *SpectraFileNewBinning = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/RCP/JetRCP_NewBinning.root");
	SpectraFileNewBinning->cd();
	TH1D *CentralNewBinning;
	// CentralNewBinning    = (TH1D *)gDirectory->Get("JetRCP_0_10");
	// CentralNewBinning->SetFillColorAlpha(kGray+2, 0);
	// CentralNewBinning->SetLineColor(kGray+2);
	// CentralNewBinning->SetMarkerColor(kGray+2);
	// CentralNewBinning->SetMarkerStyle(24);
	// CentralNewBinning->SetMarkerSize(1.5);

	int idx = 4;
    int method1 = 30;

    TH1F *hData_Old[3];

	TFile *fin1 = new TFile(Form("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/hists/Unfolded_Hist_Cen_0_10_Method_%i_Indx_%i.root",method1,idx),"READ");
    hData_Old[0] = (TH1F*)ProcessSpectraHistogram((TH1F*)fin1->Get(Form("hData_0_10_%i_%i",method1,idx)));
    TFile *fin2 = new TFile(Form("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/hists/Unfolded_Hist_Cen_10_40_Method_%i_Indx_%i.root",method1,idx),"READ");
    hData_Old[1] = (TH1F*)ProcessSpectraHistogram((TH1F*)fin2->Get(Form("hData_10_40_%i_%i",method1,idx)));
    TFile *fin3 = new TFile(Form("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/hists/Unfolded_Hist_Cen_40_80_Method_%i_Indx_%i.root",method1,idx),"READ");
    hData_Old[2] = (TH1F*)ProcessSpectraHistogram((TH1F*)fin3->Get(Form("hData_40_80_%i_%i",method1,idx)));

    int method2 = 33;

    TH1F *hData_Old_NewBinning[3];
    TH1F *hData_Old_Rebinned[3];

	TFile *fin4 = new TFile(Form("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/hists/Unfolded_Hist_Cen_0_10_Method_%i_Indx_%i.root",method2,idx),"READ");
    hData_Old_NewBinning[0] = (TH1F*)fin4->Get(Form("hData_0_10_%i_%i",method2,idx));
    TFile *fin5 = new TFile(Form("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/hists/Unfolded_Hist_Cen_10_40_Method_%i_Indx_%i.root",method2,idx),"READ");
    hData_Old_NewBinning[1] = (TH1F*)fin5->Get(Form("hData_10_40_%i_%i",method2,idx));
    TFile *fin6 = new TFile(Form("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/hists/Unfolded_Hist_Cen_40_80_Method_%i_Indx_%i.root",method2,idx),"READ");
    hData_Old_NewBinning[2] = (TH1F*)fin6->Get(Form("hData_40_80_%i_%i",method2,idx));


    for (int i = 0; i < 3; i++){
    	hData_Old[i]->Scale(1./naa[i]);
    	hData_Old_Rebinned[i] = (TH1F*)ProcessSpectraHistogram((TH1F *)hData_Old_NewBinning[i]->Rebin(nbins_jpt, Form("hData_%i", i), binning_jpt));
    	hData_Old_Rebinned[i]->Scale(1./naa[i]);
    }

    CentralNewBinning = (TH1D *)hData_Old_Rebinned[0]->Clone("JetRCP_0_10_NewBin");
    CentralNewBinning->Divide(hData_Old_Rebinned[2]);
    CentralNewBinning->Scale(taa[2]/taa[0]);

    CentralNewBinning->SetFillColorAlpha(kGray+2, 0);
	CentralNewBinning->SetLineColor(kGray+2);
	CentralNewBinning->SetMarkerColor(kGray+2);
	CentralNewBinning->SetMarkerStyle(45);
	CentralNewBinning->SetMarkerSize(1.5);

	int color[14] = {kBlack, kOrange+10, kPink-10, kViolet, kAzure+10, kTeal+10, kGreen-2, kRed+2, kMagenta+2, kBlue+2, kCyan+2, kGreen+2, kGreen, kOrange};

	for (int cs = 0; cs < filenum; cs++){
		for (int i = 0; i < 3; i++){
			dPt[i][cs]->Scale(1./dPt[i][cs]->Integral());
			dEta[i][cs]->Scale(1./dEta[i][cs]->Integral());
			dPhi[i][cs]->Scale(1./dPhi[i][cs]->Integral());

			SetColor(dPt[i][cs], color[cs]);
			SetColor(dEta[i][cs], color[cs]);
			SetColor(dPhi[i][cs], color[cs]);
		}
	}

	gStyle->SetOptTitle(0);

	TCanvas *c[3];
	TLegend *l[3][3];

	TString centrality[3] = {"Central", "MidCentral", "Peripheral"};

	for (int i = 0; i < 3; i++){
		c[i] = new TCanvas(Form("c_%s", centrality[i].Data()), Form("c_%s", centrality[i].Data()), 1200, 400);
		c[i]->Divide(3);

		c[i]->cd(1);
		l[0][i] =  new TLegend(0.1,0.8,0.9,0.9);
		gPad->SetLogy();

		for (int fi = 0; fi < filenum; fi++){
			if (fi != 3 && fi != 6) continue;
			dPt[i][fi]->Draw("EP SAME");
			dPt[i][fi]->GetXaxis()->SetRangeUser(-30, 30);
			dPt[i][fi]->GetXaxis()->SetTitle("#Delta p_{T} = p_{T,Jet}^{Reco} - p_{T,Jet}^{Gen} (GeV/#it{c})");
			dPt[i][fi]->GetYaxis()->SetTitle("arb.units");
			dPt[i][fi]->GetYaxis()->SetRangeUser(pow(10, -8), pow(10, 3));
			l[0][i]->AddEntry(dPt[i][fi],Form("%s = %.2f, %.2f", LegAppend[fi].Data(), dPt[i][fi]->GetMean(), dPt[i][fi]->GetRMS()),"p");
		}
		l[0][i]->SetTextSize(0.03);
		l[0][i]->Draw("SAME");

		c[i]->cd(2);
		l[1][i] =  new TLegend(0.1,0.8,0.9,0.9);
		gPad->SetLogy();

		for (int fi = 0; fi < filenum; fi++){
			if (fi != 3 && fi != 6) continue;
			dEta[i][fi]->Draw("EP SAME");
			dEta[i][fi]->GetXaxis()->SetRangeUser(0, 0.6);
			dEta[i][fi]->GetXaxis()->SetTitle("#Delta #eta = |#eta_{Jet}^{Reco} - #eta_{Jet}^{Gen}|");
			dEta[i][fi]->GetYaxis()->SetTitle("arb.units");
			dEta[i][fi]->GetYaxis()->SetRangeUser(pow(10, -8), pow(10, 3));
			l[1][i]->AddEntry(dEta[i][fi],Form("%s = %.2f, %.2f", LegAppend[fi].Data(), dEta[i][fi]->GetMean(), dEta[i][fi]->GetRMS()),"p");
		}
		l[1][i]->SetTextSize(0.03);
		l[1][i]->Draw("SAME");

		c[i]->cd(3);
		l[2][i] =  new TLegend(0.1,0.8,0.9,0.9);
		gPad->SetLogy();

		for (int fi = 0; fi < filenum; fi++){
			if (fi != 3 && fi != 6) continue;
			dPhi[i][fi]->Draw("EP SAME");
			dPhi[i][fi]->GetXaxis()->SetRangeUser(0, 0.6);
			dPhi[i][fi]->GetXaxis()->SetTitle("#Delta #phi = |#phi_{Jet}^{Reco} - #phi_{Jet}^{Gen}|");
			dPhi[i][fi]->GetYaxis()->SetTitle("arb.units");
			dPhi[i][fi]->GetYaxis()->SetRangeUser(pow(10, -8), pow(10, 3));
			l[2][i]->AddEntry(dPhi[i][fi],Form("%s =  %.2f, %.2f", LegAppend[fi].Data(), dPhi[i][fi]->GetMean(), dPhi[i][fi]->GetRMS()),"p");
		}
		l[2][i]->SetTextSize(0.03);
		l[2][i]->Draw("SAME");

		c[i]->SaveAs(Form("JetQA/%s.pdf", centrality[i].Data()));

		delete c[i];
	}

	gStyle->SetOptTitle(0);

	// TCanvas *e = new TCanvas("Pt Spectra", "Pt Spectra", 1200, 400);
	

	TCanvas *e[3]; 
	TLegend *l3[3];

	for (int i = 0; i < 3; i++){
		e[i] = new TCanvas(Form("Pt Spectra %i", i), Form("Pt Spectra %i", i), 1000, 800);
		e[i]->cd();

		l3[i] =  new TLegend(0.5,0.7,0.95,0.9, "", "NDC");

		gPad->SetLogy();
		gPad->SetLeftMargin(0.25);
		gPad->SetRightMargin(0.05);
		gPad->SetBottomMargin(0.15);

		hData_Old[i]->Draw("HIST E SAME");
		hData_Old[i]->SetFillColorAlpha(kBlack, 0);
		hData_Old[i]->SetLineColor(kBlack);
		hData_Old[i]->SetMarkerColor(kBlack);
		hData_Old[i]->SetMarkerStyle(29);
		hData_Old[i]->SetMarkerSize(1.5);
		hData_Old[i]->GetXaxis()->SetNdivisions(505);
		hData_Old[i]->GetYaxis()->SetNdivisions(504);
		hData_Old[i]->GetXaxis()->SetTitle("p_{T,Jet}^{D^{0}} [GeV/#it{c}]");
		hData_Old[i]->GetYaxis()->SetTitle("#frac{1}{2#pi#it{N}_{evt}} #frac{d^{2}#it{N}_{Jet}}{#it{p}_{T,jet} d#it{p}_{T,jet} d#it{#eta}} [GeV/#it{c}]^{-2}");
		hData_Old[i]->GetYaxis()->SetTitleSize(0.055);
		hData_Old[i]->GetXaxis()->SetTitleSize(0.055);
		hData_Old[i]->GetYaxis()->SetLabelSize(0.06);
		hData_Old[i]->GetXaxis()->SetLabelSize(0.06);
		hData_Old[i]->GetYaxis()->SetTitleOffset(2.0);
		hData_Old[i]->GetXaxis()->SetTitleOffset(1.1);
		hData_Old[i]->GetXaxis()->SetRangeUser(5, 15);
		hData_Old[i]->GetYaxis()->SetRangeUser(5*pow(10, -10), 2*pow(10, -3));

		l3[i]->AddEntry(hData_Old[i], "QM Result", "p");

		// hData_Old_Rebinned[i]->Draw("HIST E SAME");
		// hData_Old_Rebinned[i]->SetFillColorAlpha(kBlack, 0);
		// hData_Old_Rebinned[i]->SetLineColor(kGray+2);
		// hData_Old_Rebinned[i]->SetMarkerColor(kGray+2);
		// hData_Old_Rebinned[i]->SetMarkerStyle(45);
		// hData_Old_Rebinned[i]->SetMarkerSize(1.5);

		for (int fi = 1; fi < filenum; fi++){
			if (fi != 3 && fi != 6) continue;
			// if (fi != 4) continue;
			UnfoldedPt1D[fi][i]->Draw("EP SAME");
			UnfoldedPt1D[fi][i]->GetXaxis()->SetRangeUser(5, 15);
			UnfoldedPt1D[fi][i]->GetYaxis()->SetRangeUser(5*pow(10, -10), 2*pow(10, -3));
			SetColor(UnfoldedPt1D[fi][i], color[fi]);
			l3[i]->AddEntry(UnfoldedPt1D[fi][i],Form("%s", LegAppend[fi].Data()),"p");
		}

		

		for (int gi = 1; gi < filenum; gi++){
			if (gi != 6) continue;
			cout << UnfoldedPt2D[gi][i]->GetName() << "\t" << UnfoldedPt2D[gi][i]->Integral() << endl;
			UnfoldedPt2D[gi][i]->Draw("EP SAME");
			UnfoldedPt2D[gi][i]->GetXaxis()->SetRangeUser(5, 15);
			UnfoldedPt2D[gi][i]->GetYaxis()->SetRangeUser(5*pow(10, -10), 2*pow(10, -3));
			SetColor(UnfoldedPt2D[gi][i], color[gi], 29);
			l3[i]->AddEntry(UnfoldedPt2D[gi][i],Form("%s (2D)", LegAppend[gi].Data()),"p");
		}

		

		// Truth->Draw("EP SAME");
		// SetColor(Truth, kRed, 23);
		// if (i == 2)l3->AddEntry(Truth,Form("mFONLL"),"p");

		// for (int ti = 1; ti < filenum; ti++){
		// 	if (ti != 6) continue;
		// 	cout << TruthPt1D[ti][i]->GetName() << "\t" << TruthPt1D[ti][i]->Integral() << endl;
		// 	TruthPt1D[ti][i]->Draw("EP SAME");
		// 	TruthPt1D[ti][i]->GetXaxis()->SetRangeUser(5, 15);
		// 	TruthPt1D[ti][i]->GetYaxis()->SetRangeUser(5*pow(10, -10), 2*pow(10, -3));
		// 	SetColor(TruthPt1D[ti][i], kRed, 23);
		// 	if (i == 2)l3->AddEntry(TruthPt1D[ti][i],Form("FONLL"),"p");
		// }

		// if (i == 2)l3->AddEntry(hData_Old[i], "QM Result","p");
		// if (i == 2)l3->AddEntry(hData_Old_Rebinned[i], "QM Macro New Sample","p");

		l3[i]->SetTextSize(0.03); 
		l3[i]->Draw("SAME");

		e[i]->SaveAs(Form("JetQA/PtSpectraFrom1D_%i.pdf", i));
		e[i]->SaveAs(Form("JetQA/PtSpectraFrom1D_%i.png", i));

		delete e[i];
	}

	

	TLegend *l2;
	l2 = new TLegend(0.5,0.7,0.9,0.9);

	TCanvas *d2 = new TCanvas("RCP1D", "RCP1D", 800, 800);
	d2->cd();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);

	for (int fi = 1; fi < filenum; fi++){
		if (fi != 3 && fi != 6) continue;
		// if (fi != 4) continue;
		SetColor(RCPPt1D[fi], color[fi]);
		RCPPt1D[fi]->GetYaxis()->SetTitle("R_{CP}");
		RCPPt1D[fi]->GetXaxis()->SetTitle("p_{T,Jet}^{D^{0}} [GeV/#it{c}]");
		RCPPt1D[fi]->Draw("EP SAME");

		RCPPt1D[fi]->GetYaxis()->SetRangeUser(0, 1.4);
		RCPPt1D[fi]->GetXaxis()->SetRangeUser(5, 15);
		RCPPt1D[fi]->GetYaxis()->SetTitleSize(0.055);
		RCPPt1D[fi]->GetXaxis()->SetTitleSize(0.055);
		RCPPt1D[fi]->GetYaxis()->SetLabelSize(0.06);
		RCPPt1D[fi]->GetXaxis()->SetLabelSize(0.06);
		RCPPt1D[fi]->GetYaxis()->SetTitleOffset(1.1);
		RCPPt1D[fi]->GetXaxis()->SetTitleOffset(1.1);
		RCPPt1D[fi]->GetXaxis()->SetNdivisions(505);
		RCPPt1D[fi]->GetYaxis()->SetNdivisions(505);
		
		
		l2->AddEntry(RCPPt1D[fi], Form("%s", LegAppend[fi].Data()));
	}

	Central->Draw("EP SAME");
	// CentralNewBinning->Draw("EP SAME");
	l2->AddEntry(Central, "QM Result", "lp");
	// l2->AddEntry(CentralNewBinning, "QM Macro New Sample", "lp");

	Central->GetYaxis()->SetTitle("R_{CP}");
	Central->GetYaxis()->SetRangeUser(0, 1.2);

	
	for (int gi = 1; gi < filenum; gi++){
		if (gi != 6) continue;
		RCPPt[gi]->Draw("EP SAME");
		RCPPt[gi]->GetYaxis()->SetRangeUser(0, 1.4);
		SetColor(RCPPt[gi], color[gi], 29);
		l2->AddEntry(RCPPt[gi],Form("%s (2D)", LegAppend[gi].Data()));
	}
	

	l2->Draw("SAME");
	d2->SaveAs("JetQA/RCP1D.pdf");
	d2->SaveAs("JetQA/RCP1D.png");
	delete d2;

	// TCanvas *d = new TCanvas("RCP", "RCP", 2000, 800);
	// d->Divide(2);
	
	// for (int fi = 1; fi < filenum; fi++){
	// 	if (fi != 3 && fi != 6) continue;
	// 	d->cd(1);
	// 	SetColor(RCPPt[fi], color[fi]);
	// 	RCPPt[fi]->GetYaxis()->SetRangeUser(0,1.2);
	// 	RCPPt[fi]->Draw("EP SAME");

	// 	d->cd(2);
	// 	SetColor(RCPZ[fi], color[fi]);
	// 	RCPZ[fi]->GetYaxis()->SetRangeUser(0,2);
	// 	RCPZ[fi]->Draw("EP SAME");
	// 	l2->AddEntry(RCPZ[fi], Form("%s", LegAppend[fi].Data()));
	// }

	// d->cd(3);
	// l2->Draw("SAME");

	// d->SaveAs("JetQA/RCP.pdf");

	TLegend *l5;
	l5 = new TLegend(0.5,0.7,0.9,0.9);

	TCanvas *d3 = new TCanvas("RCPZ", "RCPZ", 800, 800);
	d3->cd();
	gPad->SetLeftMargin(0.2);
	gPad->SetBottomMargin(0.15);

	for (int gi = 1; gi < filenum; gi++){
		if (gi != 6) continue;
		// if (fi != 4) continue;
		SetColor(RCPZ[gi], color[gi]);
		RCPZ[gi]->GetYaxis()->SetTitle("R_{CP}");
		RCPZ[gi]->GetXaxis()->SetTitle("z");
		RCPZ[gi]->Draw("EP SAME");

		RCPZ[gi]->GetYaxis()->SetRangeUser(0, 1.4);
		RCPZ[gi]->GetXaxis()->SetRangeUser(0.3, 1);
		RCPZ[gi]->GetYaxis()->SetTitleSize(0.055);
		RCPZ[gi]->GetXaxis()->SetTitleSize(0.055);
		RCPZ[gi]->GetYaxis()->SetLabelSize(0.06);
		RCPZ[gi]->GetXaxis()->SetLabelSize(0.06);
		RCPZ[gi]->GetYaxis()->SetTitleOffset(1.1);
		RCPZ[gi]->GetXaxis()->SetTitleOffset(1.1);
		RCPZ[gi]->GetXaxis()->SetNdivisions(505);
		RCPZ[gi]->GetYaxis()->SetNdivisions(505);
		
		
		l5->AddEntry(RCPZ[gi], Form("%s", LegAppend[gi].Data()));
	}

	l5->Draw("SAME");
	d3->SaveAs("JetQA/RCPZ.pdf");
	d3->SaveAs("JetQA/RCPZ.png");
	delete d3;

	TCanvas *ez[3]; 
	TLegend *l6[3];

	for (int i = 0; i < 3; i++){
		ez[i] = new TCanvas(Form("Z Spectra %i", i), Form("Z Spectra %i", i), 1000, 800);
		ez[i]->cd();

		l6[i] =  new TLegend(0.5,0.7,0.95,0.9, "", "NDC");

		gPad->SetLogy();
		gPad->SetLeftMargin(0.25);
		gPad->SetRightMargin(0.05);
		gPad->SetBottomMargin(0.15);

		for (int gi = 1; gi < filenum; gi++){
			if (gi != 6) continue;
			// cout << UnfoldedPt2D[gi][i]->GetName() << "\t" << UnfoldedPt2D[gi][i]->Integral() << endl;
			UnfoldedZ2D[gi][i]->GetXaxis()->SetNdivisions(505);
			UnfoldedZ2D[gi][i]->GetYaxis()->SetNdivisions(504);
			UnfoldedZ2D[gi][i]->GetXaxis()->SetTitle("z");
			UnfoldedZ2D[gi][i]->GetYaxis()->SetTitle("#frac{1}{2#pi#it{N}_{evt}} #frac{d^{2}#it{N}_{Jet}}{z dz d#it{#eta}}");
			UnfoldedZ2D[gi][i]->GetYaxis()->SetTitleSize(0.055);
			UnfoldedZ2D[gi][i]->GetXaxis()->SetTitleSize(0.055);
			UnfoldedZ2D[gi][i]->GetYaxis()->SetLabelSize(0.06);
			UnfoldedZ2D[gi][i]->GetXaxis()->SetLabelSize(0.06);
			UnfoldedZ2D[gi][i]->GetYaxis()->SetTitleOffset(2.0);
			UnfoldedZ2D[gi][i]->GetXaxis()->SetTitleOffset(1.1);
			UnfoldedZ2D[gi][i]->Draw("EP SAME");
			UnfoldedZ2D[gi][i]->GetXaxis()->SetRangeUser(0.3, 1);
			UnfoldedZ2D[gi][i]->GetYaxis()->SetRangeUser(5*pow(10, -8), 2*pow(10, -1));
			SetColor(UnfoldedZ2D[gi][i], color[gi], 20);
			l6[i]->AddEntry(UnfoldedZ2D[gi][i],Form("%s (2D)", LegAppend[gi].Data()),"p");
		}

		

		l6[i]->SetTextSize(0.03); 
		l6[i]->Draw("SAME");

		ez[i]->SaveAs(Form("JetQA/ZSpectraFrom2D_%i.pdf", i));
		ez[i]->SaveAs(Form("JetQA/ZSpectraFrom2D_%i.png", i));

		delete ez[i];
	}
}