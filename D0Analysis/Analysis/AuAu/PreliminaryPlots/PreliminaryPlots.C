#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TObject.h"
#include "TClonesArray.h"
// #include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
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
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include <algorithm>
#include <TString.h>

using namespace std;

TString xaxisforptspectra = "p_{T,jet} [GeV/#it{c}]";
TString yaxisforptspectra = "#frac{dN_{jet}}{dp_{T,jet}} [GeV/#it{c}]^{-1}";

TString xaxisfordRspectra = "#Delta#it{R}";
TString yaxisfordRspectra = "#frac{1}{N_{jet}} #frac{#DeltaN_{jet}}{#Delta#it{R}}";

TString DirectoryName = "./";



void JetPtSpectra(){

	gROOT->ProcessLine(".x ./myStyle.C");

	// myStyle();

	TFile *SpectraFile = new TFile("April2/JetSpectra.root");

	TH1D *Central    = (TH1D *)gDirectory->Get("JetSpectra_0_10");
	TH1D *MidCentral = (TH1D *)gDirectory->Get("JetSpectra_10_40");
	TH1D *Peripheral = (TH1D *)gDirectory->Get("JetSpectra_40_80");

	TH1D *Central_Sys    = (TH1D *)gDirectory->Get("JetSpectra_0_10_Sys");
	TH1D *MidCentral_Sys = (TH1D *)gDirectory->Get("JetSpectra_10_40_Sys");
	TH1D *Peripheral_Sys = (TH1D *)gDirectory->Get("JetSpectra_40_80_Sys");

	double centralscalefactor = 10000;

	Central->Scale(centralscalefactor);
	Central_Sys->Scale(centralscalefactor);

	TH1D *tmp1 = (TH1D *)Central_Sys->Clone("tmp1");
	
	Central_Sys->Draw("E2");
	Central_Sys->SetFillColorAlpha(kGreen-2, 0.35);
	Central_Sys->SetFillStyle(3356);

	Central_Sys->GetXaxis()->SetTitle("#it{p}_{T,jet} [GeV/#it{c}]");
	Central_Sys->GetYaxis()->SetTitle("#frac{1}{2#pi#it{N}_{evt}} #frac{d^{2}#it{N}_{Jet}}{#it{p}_{T,jet} d#it{p}_{T,jet} d#it{#eta}} [GeV/#it{c}]^{-2}");


	// tmp1->Draw("E5 SAME");
	// tmp1->SetFillColorAlpha(kGreen-2, 0.30);
	Central->Draw("SAME");
	Central->SetLineColor(kGreen-2);
	Central->SetMarkerColor(kGreen-2);
	Central_Sys->GetYaxis()->SetNdivisions(503);
	Central_Sys->GetYaxis()->SetRangeUser(9*pow(10, -10), 3*pow(10, 3));
	Central_Sys->GetXaxis()->SetRangeUser(4, 15);
	Central_Sys->GetYaxis()->SetMaxDigits(2);

	Central->SetMarkerStyle(20);
	Central->SetMarkerSize(1.5);

	double midcentralscalefactor = 100.;

	MidCentral->Scale(midcentralscalefactor);
	MidCentral_Sys->Scale(midcentralscalefactor);

	TH1D *tmp2 = (TH1D *)MidCentral_Sys->Clone("tmp2");
	
	MidCentral_Sys->Draw("E2 SAME");
	MidCentral_Sys->SetFillColorAlpha(kBlue, 0.35);
	// tmp2->Draw("E5 SAME");
	// tmp2->SetFillColorAlpha(kBlue, 0.30);
	MidCentral_Sys->SetFillStyle(3365);
	MidCentral->Draw("SAME");
	MidCentral->SetLineColor(kBlue);
	MidCentral->SetMarkerColor(kBlue);
	MidCentral_Sys->GetYaxis()->SetRangeUser(pow(10, -17), pow(10, -5));

	MidCentral->SetMarkerStyle(21);
	MidCentral->SetMarkerSize(1.5);

	double peripheralscalefactor = 1.;

	Peripheral->Scale(peripheralscalefactor);
	Peripheral_Sys->Scale(peripheralscalefactor);

	Peripheral_Sys->SetBinContent(10, 0);
	Peripheral_Sys->SetBinError(10, 0);
	Peripheral->SetBinContent(10, 0);
	Peripheral->SetBinError(10, 0);

	TH1D *tmp3 = (TH1D *)Peripheral_Sys->Clone("tmp3");
	
	Peripheral_Sys->Draw("E2 SAME");
	Peripheral_Sys->SetFillColorAlpha(kBlack, 0.35);
	// tmp3->Draw("E5 SAME");
	// tmp3->SetFillColorAlpha(kBlack, 0.30);
	Peripheral_Sys->SetFillStyle(3395);
	Peripheral->Draw("SAME");
	Peripheral->SetLineColor(kBlack);
	Peripheral->SetMarkerColor(kBlack);
	Peripheral_Sys->GetYaxis()->SetRangeUser(9*pow(10, -11), 5*pow(10, 2));

	Peripheral_Sys->SetMarkerStyle(107);
	Peripheral_Sys->SetMarkerSize(1.5);
	Peripheral->SetMarkerStyle(107);
	Peripheral->SetMarkerSize(1.5);

	gPad->SetLogy();

	auto legend = new TLegend(0.90,0.67,0.65,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetTextFont(132);
    legend->AddEntry(Central,"0-10\% (x10^{4})","p");
    legend->AddEntry(MidCentral, "10-40\% (x10^{2})", "p");
    legend->AddEntry(Peripheral,"40-80\%","p");
    // legend->AddEntry(Central,"Stat. Unc.","ep");
    legend->AddEntry(Peripheral_Sys,"Sys. Unc.","f");

    legend->Draw("SAME");

	TPaveText *data = new TPaveText(0.25,0.71,0.55,0.86, "NDC");
	data->SetFillColor(0);
	data->SetTextAlign(12);
	data->SetTextSize(0.035);
	data->SetBorderSize(0);
	auto system1 = data->AddText("Au+Au #sqrt{#it{s}_{NN}} = 200 GeV");
	// system1->SetTextSize(0.035);
	data->AddText("anti-#it{k}_{T}, #it{R} = 0.4");
	data->AddText("#it{p}_{T,D^{0}} > 5 GeV/#it{c}");

	data->Draw("SAME");

    TLatex *myLabel = new TLatex();
	myLabel->SetTextFont(132);
	myLabel->SetTextColor(kBlack);
	myLabel->SetTextSize(0.04);
	myLabel->SetTextAlign(12);	
    myLabel->SetNDC(kTRUE);
    myLabel->SetText(0.27, 0.88, "#bf{STAR} #it{Preliminary}");
    myLabel->Draw("SAME");

    TLatex *ZDist = new TLatex();
	ZDist->SetTextFont(132);
	ZDist->SetTextColor(kBlack);
	ZDist->SetTextSize(0.04);
	ZDist->SetTextAlign(12);	
    ZDist->SetNDC(kTRUE);
    ZDist->SetText(0.25, 0.17, "*Fragmentation from PYTHIA 8");
    ZDist->Draw("SAME");

    gPad->SaveAs("JetPtSpectra_D0Pt5GeV_AuAu200GeV_Run14_STAR.pdf");

}


void RCP1040(){
	gROOT->ProcessLine(".x ./myStyle.C");
	gStyle->SetPadLeftMargin(0.15);

	TFile *Analysis2018 = new TFile("../psn0692/HEPData-ins1711377-v2-root.root");
    Analysis2018->cd("RCP (40-60%)");

    TGraph *ana2018[3];

    for (int i = 0; i < 3; i++){
        ana2018[i] = (TGraph *)gDirectory->Get(Form("Graph1D_y%i", i+1));
        // ana2018[i]->SetDirectory(0);
        ana2018[i]->SetTitle("Published D^{0}");
        ana2018[i]->SetLineColor(1);
        ana2018[i]->SetMarkerColor(1);
        // ana2018[i]->SetMarkerStyle(20);
    }
    
    // Analysis2018->Close();


    // TGraph *RCP0010_1020 = (TGraph *)ana2018[0]->Clone();

    const int bins = ana2018[0]->GetN();

    double x[bins];
    double rcp0010[bins];
    double rcp1020[bins];
    double rcp2040[bins];

    for (int i = 1; i <= bins; i++){
    	x[i-1] = ana2018[0]->GetPointX(i);
    	rcp0010[i-1] = ana2018[0]->GetPointY(i);
    	rcp1020[i-1] = ana2018[1]->GetPointY(i);
    	rcp2040[i-1] = ana2018[2]->GetPointY(i);

    	cout << rcp0010[i-1] << "\t" << rcp1020[i-1] << "\t" << rcp2040[i-1] << endl;

    	rcp1020[i-1] = rcp0010[i-1]/rcp1020[i-1];
    	rcp2040[i-1] = rcp0010[i-1]/rcp2040[i-1];
    }

    TGraph *RCP0010_1020 = new TGraph(bins, x, rcp1020);
    TGraph *RCP0010_2040 = new TGraph(bins, x, rcp2040);

    RCP0010_1020->SetLineColor(kBlue);
    RCP0010_2040->SetLineColor(kGreen-3);
    RCP0010_1020->SetMarkerColor(kBlue);
    RCP0010_2040->SetMarkerColor(kGreen-3);

	// myStyle();

	TFile *SpectraFile = new TFile("JetSpectra.root");

	TH1D *Central    = (TH1D *)gDirectory->Get("JetSpectra_0_10");
	TH1D *MidCentral = (TH1D *)gDirectory->Get("JetSpectra_10_40");
	// TH1D *Peripheral = (TH1D *)gDirectory->Get("JetRCP_40_80");

	TH1D *Central_Sys    = (TH1D *)gDirectory->Get("JetSpectra_0_10_Sys");
	TH1D *MidCentral_Sys = (TH1D *)gDirectory->Get("JetSpectra_10_40_Sys");
	// TH1D *Peripheral_Sys = (TH1D *)gDirectory->Get("JetRCP_40_80_Sys");

	double centralscalefactor = 1./941.23714;

	Central->Scale(centralscalefactor);
	Central_Sys->Scale(centralscalefactor);

	TH1D *tmp1 = (TH1D *)Central_Sys->Clone("tmp1");
	
	Central_Sys->Draw("E2");
	Central_Sys->SetFillColorAlpha(kGreen-2, 0.35);
	Central_Sys->SetFillStyle(3356);
	Central_Sys->SetTitleOffset(1.1, "Y");
	Central_Sys->SetMarkerStyle(20);


	Central_Sys->GetXaxis()->SetTitle("#it{p}_{T, Jet} [GeV/#it{c}]");
	Central_Sys->GetYaxis()->SetTitle("R_{CP} (/10-40\%)");
	Central_Sys->SetNdivisions(505,"y");
	// Central_Sys->SetLabelOffset(1.1, "Y");

	Central_Sys->SetBinContent(10, 0);
	Central_Sys->SetBinError(10,0);
	Central->SetBinContent(10, 0);
	Central->SetBinError(10,0);

	// tmp1->Draw("E5 SAME");
	// tmp1->SetFillColorAlpha(kGreen-2, 0.30);
	Central->Draw("SAME");
	Central->SetLineColor(kGreen-2);
	Central->SetMarkerColor(kGreen-2);
	Central->SetMarkerStyle(20);
	Central->SetMarkerSize(1.5);
	Central_Sys->GetYaxis()->SetRangeUser(0, 3);
	Central_Sys->GetXaxis()->SetRangeUser(4, 15);
	Central_Sys->GetYaxis()->SetMaxDigits(2);

	RCP0010_1020->Draw("P SAME");
	RCP0010_2040->Draw("P SAME");



	double midcentralscalefactor = 1./391.35550;

	MidCentral->Scale(midcentralscalefactor);
	MidCentral_Sys->Scale(midcentralscalefactor);

	Central->Divide(MidCentral);
	Central_Sys->Divide(MidCentral_Sys);

	TH1D *tmp2 = (TH1D *)MidCentral_Sys->Clone("tmp2");
	
	// MidCentral_Sys->Draw("E2 SAME");
	MidCentral_Sys->SetFillColorAlpha(kBlue, 0.35);
	MidCentral_Sys->SetMarkerStyle(21);
	// tmp2->Draw("E5 SAME");
	// tmp2->SetFillColorAlpha(kBlue, 0.30);
	MidCentral_Sys->SetFillStyle(3356);
	// MidCentral->Draw("SAME");
	MidCentral->SetLineColor(kBlue);
	MidCentral->SetMarkerColor(kBlue);
	MidCentral->SetMarkerStyle(21);
	MidCentral->SetMarkerSize(1.5);
	// MidCentral_Sys->GetYaxis()->SetRangeUser(pow(10, -17), pow(10, -5));

	MidCentral_Sys->SetBinContent(10, 0);
	MidCentral_Sys->SetBinError(10,0);
	MidCentral->SetBinContent(10, 0);
	MidCentral->SetBinError(10,0);

	TF1 *one = new TF1("One", "1", 0, 100);
	one->SetLineColorAlpha(kBlack, 1.);
	one->SetLineStyle(9);

	one->Draw("SAME");

	TH1D *CentralNCollError = new TH1D("CentralNCollError", "CentralNCollError", 1, 4.4, 4.6);
	CentralNCollError->SetBinContent(1, 1.);
	CentralNCollError->SetBinError(1, TMath::Sqrt(pow(30.21318/391.35550, 2) + pow(26.27357/941.23714, 2)));

	CentralNCollError->SetFillColor(kGreen-2);
	CentralNCollError->SetFillStyle(1001);
	CentralNCollError->SetMarkerStyle(1);
	CentralNCollError->SetMarkerColor(kGreen-2);
	CentralNCollError->Draw("E2SAME");

	TH1D *MidCentralNCollError = new TH1D("MidCentralNCollError", "MidCentralNCollError", 1, 4.6, 4.8);
	MidCentralNCollError->SetBinContent(1, 1.);
	MidCentralNCollError->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(30.21318/391.35550, 2)));

	MidCentralNCollError->SetFillColor(kBlue);
	MidCentralNCollError->SetFillStyle(1001);
	MidCentralNCollError->SetMarkerStyle(1);
	MidCentralNCollError->SetMarkerColor(kBlue);
	// MidCentralNCollError->Draw("E2SAME");
	// gPad->SetLogy();

	auto legend = new TLegend(0.90,0.67,0.53,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(132);
    legend->AddEntry(Central,"0-10\%","p");
    // legend->AddEntry(MidCentral, "10-40\%", "p");
    // legend->AddEntry(Peripheral,"40-80\%","p");
    // legend->AddEntry(Central,"Stat. Unc.","ep");
    legend->AddEntry(Central_Sys,"Sys. Unc. (Uncorrelated)","f");
    legend->AddEntry(CentralNCollError, "0-10\% <N_{Coll}> Unc.", "f");
    // legend->AddEntry(MidCentralNCollError, "10-40\% <N_{Coll}> Unc.", "f");

    legend->Draw("SAME");

	TPaveText *data = new TPaveText(0.18,0.71,0.50,0.86, "NDC");
	data->SetFillColor(0);
	data->SetTextAlign(12);
	data->SetTextSize(0.04);
	data->SetBorderSize(0);
	auto system1 = data->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
	// system1->SetTextSize(0.035);
	data->AddText("anti-#it{k}_{T}, R = 0.4");
	data->AddText("#it{p}_{T, D^{0}} > 5 GeV/#it{c}");

	data->Draw("SAME");

    TLatex *myLabel = new TLatex();
	myLabel->SetTextFont(132);
	myLabel->SetTextColor(kBlack);
	myLabel->SetTextSize(0.04);
	myLabel->SetTextAlign(12);	
    myLabel->SetNDC(kTRUE);
    myLabel->SetText(0.18, 0.60, "#bf{STAR} #it{Preliminary}");
    myLabel->Draw("SAME");

    gPad->SaveAs("RCP1040_D0Pt5GeV_AuAu200GeV_Run14_STAR.pdf");
}


void RCP(){
	gROOT->ProcessLine(".x ./myStyle.C");
	gStyle->SetPadLeftMargin(0.15);

	// myStyle();

	TFile *SpectraFile = new TFile("April2/JetRCP.root");

	TH1D *Central    = (TH1D *)gDirectory->Get("JetRCP_0_10");
	TH1D *MidCentral = (TH1D *)gDirectory->Get("JetRCP_10_40");
	// TH1D *Peripheral = (TH1D *)gDirectory->Get("JetRCP_40_80");

	TH1D *Central_Sys    = (TH1D *)gDirectory->Get("JetRCP_0_10_Sys");
	TH1D *MidCentral_Sys = (TH1D *)gDirectory->Get("JetRCP_10_40_Sys");
	// TH1D *Peripheral_Sys = (TH1D *)gDirectory->Get("JetRCP_40_80_Sys");

	double centralscalefactor = 1;

	Central->Scale(centralscalefactor);
	Central_Sys->Scale(centralscalefactor);

	TH1D *tmp1 = (TH1D *)Central_Sys->Clone("tmp1");
	
	Central_Sys->Draw("E2");
	Central_Sys->SetFillColorAlpha(kGreen-2, 0.35);
	Central_Sys->SetFillStyle(3356);
	Central_Sys->SetTitleOffset(1.1, "Y");
	Central_Sys->SetMarkerStyle(20);


	Central_Sys->GetXaxis()->SetTitle("#it{p}_{T,jet} [GeV/#it{c}]");
	Central_Sys->GetYaxis()->SetTitle("R^{*}_{CP} (/40-80\%)");
	Central_Sys->SetNdivisions(505,"y");
	// Central_Sys->SetLabelOffset(1.1, "Y");

	Central_Sys->SetBinContent(10, 0);
	Central_Sys->SetBinError(10,0);
	Central->SetBinContent(10, 0);
	Central->SetBinError(10,0);



	// tmp1->Draw("E5 SAME");
	// tmp1->SetFillColorAlpha(kGreen-2, 0.30);
	Central->Draw("SAME");
	Central->SetLineColor(kGreen-2);
	Central->SetMarkerColor(kGreen-2);
	Central->SetMarkerStyle(20);
	Central->SetMarkerSize(1.5);
	Central_Sys->GetYaxis()->SetRangeUser(0, 2.9);
	Central_Sys->GetXaxis()->SetRangeUser(4, 15);
	Central_Sys->GetYaxis()->SetMaxDigits(2);


	double midcentralscalefactor = 1.;

	MidCentral->Scale(midcentralscalefactor);
	MidCentral_Sys->Scale(midcentralscalefactor);

	TH1D *tmp2 = (TH1D *)MidCentral_Sys->Clone("tmp2");
	
	MidCentral_Sys->Draw("E2 SAME");
	MidCentral_Sys->SetFillColorAlpha(kBlue, 0.35);
	MidCentral_Sys->SetMarkerStyle(21);
	// tmp2->Draw("E5 SAME");
	// tmp2->SetFillColorAlpha(kBlue, 0.30);
	MidCentral_Sys->SetFillStyle(3365);
	MidCentral->Draw("SAME");
	MidCentral->SetLineColor(kBlue);
	MidCentral->SetMarkerColor(kBlue);
	MidCentral->SetMarkerStyle(21);
	MidCentral->SetMarkerSize(1.5);
	// MidCentral_Sys->GetYaxis()->SetRangeUser(pow(10, -17), pow(10, -5));

	MidCentral_Sys->SetBinContent(10, 0);
	MidCentral_Sys->SetBinError(10,0);
	MidCentral->SetBinContent(10, 0);
	MidCentral->SetBinError(10,0);

	TF1 *one = new TF1("One", "1", 0, 100);
	one->SetLineColorAlpha(kBlack, 1.);
	one->SetLineStyle(9);

	one->Draw("SAME");

	TH1D *CentralNCollError = new TH1D("CentralNCollError", "CentralNCollError", 1, 4.4, 4.6);
	CentralNCollError->SetBinContent(1, 1.);
	CentralNCollError->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(26.27357/941.23714, 2)));

	CentralNCollError->SetFillColor(kGreen-2);
	CentralNCollError->SetFillStyle(1001);
	CentralNCollError->SetMarkerStyle(1);
	CentralNCollError->SetMarkerColor(kGreen-2);
	CentralNCollError->Draw("E2SAME");

	TH1D *MidCentralNCollError = new TH1D("MidCentralNCollError", "MidCentralNCollError", 1, 4.6, 4.8);
	MidCentralNCollError->SetBinContent(1, 1.);
	MidCentralNCollError->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(30.21318/391.35550, 2)));

	MidCentralNCollError->SetFillColor(kBlue);
	MidCentralNCollError->SetFillStyle(1001);
	MidCentralNCollError->SetMarkerStyle(1);
	MidCentralNCollError->SetMarkerColor(kBlue);
	MidCentralNCollError->Draw("E2SAME");
	// gPad->SetLogy();

	auto legend = new TLegend(0.90,0.67,0.53,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(132);
    legend->AddEntry(Central,"0-10\%","p");
    legend->AddEntry(MidCentral, "10-40\%", "p");
    // legend->AddEntry(Peripheral,"40-80\%","p");
    // legend->AddEntry(Central,"Stat. Unc.","ep");
    legend->AddEntry(MidCentral_Sys,"Sys. Unc. (Uncorrelated)","f");
    legend->AddEntry(CentralNCollError, "0-10\% <N_{Coll}> Unc.", "f");
    legend->AddEntry(MidCentralNCollError, "10-40\% <N_{Coll}> Unc.", "f");

    legend->Draw("SAME");

	TPaveText *data = new TPaveText(0.18,0.71,0.50,0.86, "NDC");
	data->SetFillColor(0);
	data->SetTextAlign(12);
	data->SetTextSize(0.035);
	data->SetBorderSize(0);
	auto system1 = data->AddText("Au+Au #sqrt{#it{s}_{NN}} = 200 GeV");
	// system1->SetTextSize(0.035);
	data->AddText("anti-#it{k}_{T}, #it{R} = 0.4");
	data->AddText("#it{p}_{T,D^{0}} > 5 GeV/#it{c}");

	data->Draw("SAME");

    TLatex *myLabel = new TLatex();
	myLabel->SetTextFont(132);
	myLabel->SetTextColor(kBlack);
	myLabel->SetTextSize(0.04);
	myLabel->SetTextAlign(12);	
    myLabel->SetNDC(kTRUE);
    myLabel->SetText(0.20, 0.88, "#bf{STAR} #it{Preliminary}");
    myLabel->Draw("SAME");

    TLatex *ZDist = new TLatex();
	ZDist->SetTextFont(132);
	ZDist->SetTextColor(kBlack);
	ZDist->SetTextSize(0.04);
	ZDist->SetTextAlign(12);	
    ZDist->SetNDC(kTRUE);
    ZDist->SetText(0.45, 0.20, "*Fragmentation from PYTHIA 8");
    ZDist->Draw("SAME");

    gPad->SaveAs("RCP_D0Pt5GeV_AuAu200GeV_Run14_STAR.pdf");
}

void DeltaR(){
	gROOT->ProcessLine(".x ./myStyle.C");

	// myStyle();

	TFile *SpectraFile = new TFile("April2/DeltaR.root");

	TH1D *Central    = (TH1D *)gDirectory->Get("DR_0_10");
	TH1D *MidCentral = (TH1D *)gDirectory->Get("DR_10_40");
	TH1D *Peripheral = (TH1D *)gDirectory->Get("DR_40_80");

	TH1D *Central_Sys    = (TH1D *)gDirectory->Get("DR_0_10_Sys");
	TH1D *MidCentral_Sys = (TH1D *)gDirectory->Get("DR_10_40_Sys");
	TH1D *Peripheral_Sys = (TH1D *)gDirectory->Get("DR_40_80_Sys");

	double centralscalefactor = 1;

	Central->Scale(centralscalefactor);
	Central_Sys->Scale(centralscalefactor);

	TH1D *tmp1 = (TH1D *)Central_Sys->Clone("tmp1");
	
	Central_Sys->Draw("E2");
	Central_Sys->SetFillColorAlpha(kGreen-2, 0.35);
	Central_Sys->SetFillStyle(3356);

	Central_Sys->GetXaxis()->SetTitle("r");
	Central_Sys->GetYaxis()->SetTitle("#frac{1}{#it{N}_{jet}} #frac{d#it{N}_{jet}}{dr}");
	Central_Sys->SetNdivisions(504,"y");


	// tmp1->Draw("E5 SAME");
	// tmp1->SetFillColorAlpha(kGreen-2, 0.30);
	Central->Draw("SAME");
	Central->SetLineColor(kGreen-2);
	Central->SetMarkerColor(kGreen-2);
	Central_Sys->GetYaxis()->SetRangeUser(5*pow(10, -3), 9*pow(10, 3));
	Central_Sys->GetXaxis()->SetRangeUser(0, 0.2);
	Central_Sys->GetYaxis()->SetMaxDigits(2);


	double midcentralscalefactor = 1.;

	MidCentral->Scale(midcentralscalefactor);
	MidCentral_Sys->Scale(midcentralscalefactor);

	TH1D *tmp2 = (TH1D *)MidCentral_Sys->Clone("tmp2");
	
	MidCentral_Sys->Draw("E2 SAME");
	MidCentral_Sys->SetFillColorAlpha(kBlue, 0.35);
	// tmp2->Draw("E5 SAME");
	// tmp2->SetFillColorAlpha(kBlue, 0.30);
	MidCentral_Sys->SetFillStyle(3365);
	MidCentral->Draw("SAME");
	MidCentral->SetLineColor(kBlue);
	MidCentral->SetMarkerColor(kBlue);
	// MidCentral_Sys->GetYaxis()->SetRangeUser(pow(10, -17), pow(10, -5));

	double peripheralscalefactor = 1.;

	Peripheral->Scale(peripheralscalefactor);
	Peripheral_Sys->Scale(peripheralscalefactor);

	TH1D *tmp3 = (TH1D *)Peripheral_Sys->Clone("tmp3");
	
	Peripheral_Sys->Draw("E2 SAME");
	Peripheral_Sys->SetFillColorAlpha(kBlack, 0.35);
	// tmp3->Draw("E5 SAME");
	// tmp3->SetFillColorAlpha(kBlack, 0.30);
	Peripheral_Sys->SetFillStyle(3395);
	Peripheral->Draw("SAME");
	Peripheral->SetLineColor(kBlack);
	Peripheral->SetMarkerColor(kBlack);
	// Peripheral_Sys->GetYaxis()->SetRangeUser(pow(10, -17), pow(10, -5));

	for (int i = 1; i <= MidCentral->GetNbinsX(); i++){
		cout << "Bin " << i << " = " << MidCentral->GetBinError(i) << "\t" << MidCentral_Sys->GetBinError(i) << endl;
		cout << "Bin " << i << " = " << Peripheral->GetBinError(i) << "\t" << Peripheral_Sys->GetBinError(i) << endl;
	}

	Central->SetMarkerStyle(20);
	Central->SetMarkerSize(1.5);

	MidCentral->SetMarkerStyle(21);
	MidCentral->SetMarkerSize(1.5);
	Peripheral_Sys->SetMarkerStyle(107);
	Peripheral_Sys->SetMarkerSize(1.5);
	Peripheral->SetMarkerStyle(107);
	Peripheral->SetMarkerSize(1.5);

	gPad->SetLogy();

	auto legend = new TLegend(0.65,0.67,0.90,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetTextFont(132);
    legend->AddEntry(Central,"0-10\%","p");
    legend->AddEntry(MidCentral, "10-40\%", "p");
    legend->AddEntry(Peripheral,"40-80\%","p");
    // legend->AddEntry(Central,"Stat. Unc.","ep");
    legend->AddEntry(Peripheral_Sys,"Sys. Unc.","f");

    legend->Draw("SAME");

	TPaveText *data = new TPaveText(0.25,0.65,0.55,0.86, "NDC");
	data->SetFillColor(0);
	data->SetTextAlign(12);
	data->SetTextSize(0.035);
	data->SetBorderSize(0);
	auto system1 = data->AddText("Au+Au #sqrt{#it{s}_{NN}} = 200 GeV");
	// system1->SetTextSize(0.035);
	data->AddText("anti-#it{k}_{T}, #it{R} = 0.4");
	data->AddText("#it{p}_{T,D^{0}} > 5 GeV/#it{c}");
	data->AddText("5 < #it{p}_{T,jet} [GeV/#it{c}] < 15");

	data->Draw("SAME");

    TLatex *myLabel = new TLatex();
	myLabel->SetTextFont(132);
	myLabel->SetTextColor(kBlack);
	myLabel->SetTextSize(0.04);
	myLabel->SetTextAlign(12);	
    myLabel->SetNDC(kTRUE);
    myLabel->SetText(0.27, 0.88, "#bf{STAR} #it{Preliminary}");
    myLabel->Draw("SAME");

    TLatex *ZDist = new TLatex();
	ZDist->SetTextFont(132);
	ZDist->SetTextColor(kBlack);
	ZDist->SetTextSize(0.035);
	ZDist->SetTextAlign(12);	
    ZDist->SetNDC(kTRUE);
    ZDist->SetText(0.25, 0.25, "#splitline{*Fragmentation}{from PYTHIA 8}");
    ZDist->Draw("SAME");

    gPad->SaveAs("D0DeltaR_D0Pt5GeV_AuAu200GeV_Run14_STAR.pdf");

}

void DeltaRRatio(){
	gROOT->ProcessLine(".x ./myStyle.C");
	gStyle->SetPadLeftMargin(0.15);

	// myStyle();

	TFile *SpectraFile = new TFile("April2/DRRatios.root");

	TH1D *Central    = (TH1D *)gDirectory->Get("DR_Ratio_0_10");
	TH1D *MidCentral = (TH1D *)gDirectory->Get("DR_Ratio_10_40");
	// TH1D *Peripheral = (TH1D *)gDirectory->Get("JetRCP_40_80");

	TH1D *Central_Sys    = (TH1D *)gDirectory->Get("DR_Ratio_0_10_Sys");
	TH1D *MidCentral_Sys = (TH1D *)gDirectory->Get("DR_Ratio_10_40_Sys");
	// TH1D *Peripheral_Sys = (TH1D *)gDirectory->Get("JetRCP_40_80_Sys");

	double centralscalefactor = 1;

	Central->Scale(centralscalefactor);
	Central_Sys->Scale(centralscalefactor);

	TH1D *tmp1 = (TH1D *)Central_Sys->Clone("tmp1");
	
	Central_Sys->Draw("E2");
	Central_Sys->SetFillColorAlpha(kGreen-2, 0.35);
	Central_Sys->SetFillStyle(3356);
	Central_Sys->SetTitleOffset(1.1, "Y");

	Central_Sys->GetXaxis()->SetTitle("r");
	Central_Sys->GetYaxis()->SetTitle("Ratio* #scale[2.0]{(}#frac{1}{#it{N}_{jet}} #frac{d#it{N}_{jet}}{dr}#scale[2.0]{)}");
	Central_Sys->SetNdivisions(504,"y");
	// Central_Sys->SetLabelOffset(1.1, "Y");

	// tmp1->Draw("E5 SAME");
	// tmp1->SetFillColorAlpha(kGreen-2, 0.30);
	Central->Draw("SAME");
	Central->SetLineColor(kGreen-2);
	Central->SetMarkerColor(kGreen-2);
	Central_Sys->GetYaxis()->SetRangeUser(0.7, 2.2);
	Central_Sys->GetXaxis()->SetRangeUser(0, 0.2);
	Central_Sys->GetYaxis()->SetMaxDigits(2);


	double midcentralscalefactor = 1.;

	MidCentral->Scale(midcentralscalefactor);
	MidCentral_Sys->Scale(midcentralscalefactor);

	TH1D *tmp2 = (TH1D *)MidCentral_Sys->Clone("tmp2");
	
	MidCentral_Sys->Draw("E2 SAME");
	MidCentral_Sys->SetFillColorAlpha(kBlue, 0.35);
	// tmp2->Draw("E5 SAME");
	// tmp2->SetFillColorAlpha(kBlue, 0.30);
	MidCentral_Sys->SetFillStyle(3365);
	MidCentral->Draw("SAME");
	MidCentral->SetLineColor(kBlue);
	MidCentral->SetMarkerColor(kBlue);
	// MidCentral_Sys->GetYaxis()->SetRangeUser(pow(10, -17), pow(10, -5));

	for (int i = 1; i <= MidCentral->GetNbinsX(); i++){
		cout << "Bin " << i << " = " << MidCentral->GetBinError(i) << "\t" << MidCentral_Sys->GetBinError(i) << endl;
	}

	TF1 *one = new TF1("One", "1", 0, 100);
	one->SetLineColorAlpha(kBlack, 1.);
	one->SetLineStyle(9);

	one->Draw("SAME");

	Central->SetMarkerStyle(20);
	Central->SetMarkerSize(1.5);

	MidCentral->SetMarkerStyle(21);
	MidCentral->SetMarkerSize(1.5);

	// TH1D *CentralNCollError = new TH1D("CentralNCollError", "CentralNCollError", 1, 0, 0.005);
	// CentralNCollError->SetBinContent(1, 1.);
	// CentralNCollError->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(26.27357/941.23714, 2)));

	// CentralNCollError->SetFillColor(kGreen-3);
	// CentralNCollError->SetFillStyle(1001);
	// CentralNCollError->SetMarkerStyle(1);
	// CentralNCollError->SetMarkerColor(kGreen-3);
	// CentralNCollError->Draw("E2SAME");

	// gPad->SetLogy();

	auto legend = new TLegend(0.90,0.67,0.53,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(132);
    legend->AddEntry(Central,"(0-10\%)/(40-80\%)","p");
    legend->AddEntry(MidCentral, "(10-40\%)/(40-80\%)", "p");
    // legend->AddEntry(Peripheral,"40-80\%","p");
    // legend->AddEntry(Central,"Stat. Unc.","ep");
    legend->AddEntry(MidCentral_Sys,"Sys. Unc. (Uncorrelated)","f");
    // legend->AddEntry(CentralNCollError, "<N_{Coll}> Unc.", "f");

    legend->Draw("SAME");

	TPaveText *data = new TPaveText(0.18,0.71,0.55,0.86, "NDC");
	data->SetFillColor(0);
	data->SetTextAlign(12);
	data->SetTextSize(0.035);
	data->SetBorderSize(0);
	auto system1 = data->AddText("Au+Au #sqrt{#it{s}_{NN}} = 200 GeV");
	// system1->SetTextSize(0.035);
	data->AddText("anti-#it{k}_{T}, #it{R} = 0.4");
	data->AddText("#it{p}_{T,D^{0}} > 5 GeV/#it{c}");
	data->AddText("5 < #it{p}_{T,jet} [GeV/#it{c}] < 15");

	data->Draw("SAME");

    TLatex *myLabel = new TLatex();
	myLabel->SetTextFont(132);
	myLabel->SetTextColor(kBlack);
	myLabel->SetTextSize(0.04);
	myLabel->SetTextAlign(12);	
    myLabel->SetNDC(kTRUE);
    myLabel->SetText(0.20, 0.88, "#bf{STAR} #it{Preliminary}");
    myLabel->Draw("SAME");

    TLatex *ZDist = new TLatex();
	ZDist->SetTextFont(132);
	ZDist->SetTextColor(kBlack);
	ZDist->SetTextSize(0.035);
	ZDist->SetTextAlign(12);	
    ZDist->SetNDC(kTRUE);
    ZDist->SetText(0.20, 0.20, "#splitline{*Fragmentation}{from PYTHIA 8}");
    ZDist->Draw("SAME");

    gPad->SaveAs("DeltaRRatio_D0Pt5GeV_AuAu200GeV_Run14_STAR.pdf");
}


void dRComparisonWithPythia(){
	gROOT->ProcessLine(".x ./myStyle.C");
	TFile *f = new TFile("../Response2022/ResponseMatrix_Prelim_Central.root");
	TH1F *Pythia = (TH1F *)f->Get("hTrueMCdR");
	Pythia->Scale(1./Pythia->Integral());

	for (int i = 1; i <= Pythia->GetNbinsX(); i++){
		Pythia->SetBinContent(i, Pythia->GetBinContent(i)/Pythia->GetBinWidth(i));
		Pythia->SetBinError(i, Pythia->GetBinError(i)/Pythia->GetBinWidth(i));
	}

	TFile *SpectraFile = new TFile("DeltaR.root");
	TH1D *Peripheral = (TH1D *)gDirectory->Get("DR_40_80");
	TH1D *Peripheral_Sys = (TH1D *)gDirectory->Get("DR_40_80_Sys");

	Peripheral_Sys->SetFillColorAlpha(kBlack, 0.35);
	Peripheral_Sys->SetFillStyle(3356);
	Peripheral->SetLineColor(kBlack);
	Peripheral->SetMarkerColor(kBlack);

	Pythia->SetLineColor(kGreen-3);
	Pythia->SetMarkerColor(kGreen-3);

	Peripheral_Sys->GetYaxis()->SetTitle("#frac{1}{N_{Jet}} #frac{dN_{Jet}}{dr}");
	Peripheral_Sys->GetXaxis()->SetTitle("r");
	// Peripheral_Sys->Draw("E2");
	// Peripheral->Draw("SAME");
	// Pythia->Draw("SAME");

	TH1D *RatioPythia = (TH1D *)Peripheral->Clone();
	TH1D *RatioPythia_Sys = (TH1D *)Peripheral_Sys->Clone();

	for (int i = 1; i <= Peripheral->GetNbinsX(); i++){
		RatioPythia->SetBinContent(i, Peripheral->GetBinContent(i)/Pythia->GetBinContent(i));
		RatioPythia->SetBinError(i, TMath::Sqrt(pow(Peripheral->GetBinError(i), 2) + pow(Pythia->GetBinError(i), 2)));

		RatioPythia_Sys->SetBinContent(i, RatioPythia->GetBinContent(i));
		// RatioPythia_Sys->SetBinError(i, Peripheral_Sys->GetBinError(i));
	}

	RatioPythia_Sys->Draw("E2");
	RatioPythia->Draw("SAME");

	gPad->SetLogy();
}


void JetPtComparisonWithPythia(){

	gROOT->ProcessLine(".x ./myStyle.C");


	TFile *f = new TFile("../Response2022/ResponseMatrix_Prelim_Central.root");
	TH1F *Pythia = (TH1F *)f->Get("hTrueMCPt");

	TFile *SpectraFile = new TFile("JetSpectra.root");

	TH1D *Peripheral = (TH1D *)gDirectory->Get("JetSpectra_40_80");
	TH1D *Peripheral_Sys = (TH1D *)gDirectory->Get("JetSpectra_40_80_Sys");

	Peripheral_Sys->GetXaxis()->SetTitle("#it{p}_{T, Jet} [GeV/#it{c}]");
	Peripheral_Sys->GetYaxis()->SetTitle("#frac{1}{2#piN_{Evt}} #frac{d^{2}N_{Jet}}{#it{p}_{T} d#it{p}_{T} d#eta} [GeV/#it{c}]^{-2}");

	Pythia->Scale(1./Pythia->Integral());

	TH1F *Weightfactor = (TH1F *)Pythia->Clone("tmp");
	Weightfactor->Divide(Peripheral);

	Pythia->Divide(Weightfactor);

	// Pythia->Scale(56.62475);

	Peripheral_Sys->SetMarkerStyle(107);
	Peripheral_Sys->SetMarkerSize(1.5);
	Peripheral->SetMarkerStyle(107);
	Peripheral->SetMarkerSize(1.5);


	Pythia->SetMarkerStyle(29);
	Pythia->SetMarkerSize(2.0);
	Pythia->SetMarkerColor(kBlack);

	Pythia->Draw("SAME");
	Peripheral_Sys->Draw("E2 SAME");
	Peripheral->Draw("SAME");
	
}

void InvariantMassComparison(){

	TH1::SetDefaultSumw2();

	gROOT->ProcessLine(".x ./myStyle.C");
	TFile *SpectraFile = new TFile("../InvMassJets/D0Spectra.root");
	TH1D *InvMass = (TH1D *)SpectraFile->Get("Spectra_0-80");
	TFile *sPlotFile = new TFile("Histograms3_Test.root");
	TH1D *sPlot = (TH1D *)sPlotFile->Get("D0Pt");

	TFile *fin = new TFile("../D0Yield_Feb22.root");
    fin->cd("D0Tagger");
    TH1F *hCent = (TH1F *)gDirectory->Get("hCentralityWeightedAfterCuts");

    double NEvents = hCent->Integral(1, 16);


	const int nSpectraBinsD0Pt = 11;
    double D0SpectraPtBins[nSpectraBinsD0Pt + 1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    Double_t mids[nSpectraBinsD0Pt] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 7.0, 9.0};

    // double yieldunbinned[nSpectraBinsD0Pt]    = {0., 0., 1842.33, 4568.98, 6539.62, 5971.52, 8201.72, 3204.77, 1234.3, 548.535, 74.7504};
    // double yieldunbinnederr[nSpectraBinsD0Pt] = {0., 0., 77.3385, 119.491, 131.101, 114.621, 138.507, 104.923, 50.0427, 357.791, 16.6071};

    TFile *eff = new TFile("../D0ReconstructionEfficiency/effplotsredone.root");
    int efficiencycurvetopullfromthefile = 7; // The efficiency curves have a different arrangement of centrality bins than I do.
    TH1F *eff_hist = (TH1F *)eff->Get(Form("heffBinD0_%i", efficiencycurvetopullfromthefile));	


	ifstream YieldFiles;
	YieldFiles.open(Form("./RooYield/%s.txt", "0-80"), ios::binary);

	double yieldunbinned[nSpectraBinsD0Pt];
	double yieldunbinnederr[nSpectraBinsD0Pt];
	double tmp[nSpectraBinsD0Pt];

	for (int d0ptbin = 0; d0ptbin < 10; d0ptbin++){
        YieldFiles >> yieldunbinned[d0ptbin] >> yieldunbinnederr[d0ptbin] >> tmp[d0ptbin] >> tmp[d0ptbin];
    }
    YieldFiles.close();

    for (int d0ptbin = 0; d0ptbin < 10; d0ptbin++){
    	cout << yieldunbinned[d0ptbin] << "\t" << yieldunbinnederr[d0ptbin] << endl;
    }

    for (int d0ptbin = 3; d0ptbin <= nSpectraBinsD0Pt; d0ptbin++){
        sPlot->SetBinContent(d0ptbin, sPlot->GetBinContent(d0ptbin));
        sPlot->SetBinError(d0ptbin, sPlot->GetBinError(d0ptbin));

        InvMass->SetBinContent(d0ptbin, yieldunbinned[d0ptbin - 2]);
        InvMass->SetBinError(d0ptbin, yieldunbinnederr[d0ptbin - 2]);
    }

    // InvMass->Divide(eff_hist);

    sPlot->Scale(1./NEvents);
    InvMass->Scale(1./NEvents);

    InvMass->SetLineColor(kBlue);
	InvMass->SetMarkerColor(kBlue);

	// InvMass->Draw();

	sPlot->SetLineColor(kGreen-2);
	sPlot->SetMarkerColor(kGreen-2);

	// sPlot->Draw("SAME");

	sPlot->Divide(InvMass);
	sPlot->Draw();
	sPlot->GetYaxis()->SetNdivisions(504);
	sPlot->GetYaxis()->SetTitle("#frac{sPlot}{Like Sign Background Subtraction}");
	sPlot->GetYaxis()->SetRangeUser(0.7, 1.5);

	auto legend = new TLegend(0.90,0.67,0.53,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(132);
    legend->AddEntry(sPlot,"Ratio","p");
 

    TPaveText *data = new TPaveText(0.25,0.71,0.50,0.86, "NDC");
	data->SetFillColor(0);
	data->SetTextAlign(12);
	data->SetTextSize(0.04);
	data->SetBorderSize(0);
	auto system1 = data->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
	auto system2 = data->AddText("1 < #it{p}_{T, K#pi} [GeV/#it{c}] < 10");

	data->Draw();

	TF1 *one = new TF1("One", "1", 0, 100);
	one->SetLineColorAlpha(kBlack, 1.);
	one->SetLineStyle(9);

	one->Draw("SAME");

	gPad->SetLogy(0);

}

void SingleParticleEmbeddingPlot(){
	gROOT->ProcessLine(".x ./myStyle.C");

	TFile *SPFile = new TFile("BKGHistogramsFINAL.root");

	TH1D *Central    = (TH1D *)gDirectory->Get("hBKG_0_10");
	TH1D *MidCentral    = (TH1D *)gDirectory->Get("hBKG_10_40");
	TH1D *Peripheral    = (TH1D *)gDirectory->Get("hBKG_40_80");
	// TH1D *MidCentral = (TH1D *)gDirectory->Get("DR_Ratio_10_40");

	Central->SetLineColor(kGreen-2);
	Central->SetMarkerColor(kGreen-2);
	Central->Scale(1./Central->Integral());

	Central->GetYaxis()->SetRangeUser(2*pow(10, -6), 3);
	Central->GetXaxis()->SetRangeUser(-20, 20);
	Central->GetXaxis()->SetTitle("#Delta#it{p}_{T, SP Jet} [GeV/#it{c}]");
	Central->GetYaxis()->SetTitle("#frac{1}{N_{Jet}} #frac{dN_{jet}}{d#it{p}_{T}} [GeV/#it{c}]^{-1}");

	MidCentral->SetLineColor(kBlue);
	MidCentral->SetMarkerColor(kBlue);
	MidCentral->Scale(1./MidCentral->Integral());

	Peripheral->SetLineColor(kBlack);
	Peripheral->SetMarkerColor(kBlack);
	Peripheral->Scale(1./Peripheral->Integral());

	Central->Draw();
	MidCentral->Draw("SAME");
	Peripheral->Draw("SAME");

	gPad->SetLogy();

	auto legend = new TLegend(0.90,0.67,0.70,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetTextFont(132);
    legend->AddEntry(Central,"0-10\%","p");
    legend->AddEntry(MidCentral, "10-40\%", "p");
    legend->AddEntry(Peripheral,"40-80\%","p");

    legend->Draw("SAME");

    gPad->SaveAs("SingleParticleEmbeddingPlot_pp200GeV.pdf");
}

void CCNU(){
	gROOT->ProcessLine(".x ./myStyle.C");

	const int nrbinsfortheory = 8;
	double drbinsfortheory[nrbinsfortheory + 1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
	double drbinsmidfortheory[nrbinsfortheory ] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375};
	double ratiotheory[nrbinsfortheory] = {0.95405, 1.03273, 1.09413, 1.12471, 1.28143, 1.3403, 1.26103, 1.27431};

	// TGraph *theory = new TGraph(nrbinsfortheory, drbinsmidfortheory, ratiotheory);
	// theory->SetMarkerStyle(29);
	// theory->SetMarkerSize(2.0);

	TH1D *theory = new TH1D("Theory", "Theory", nrbinsfortheory, drbinsfortheory);
	for (int i = 1; i <= nrbinsfortheory; i++){
		theory->SetBinContent(i,ratiotheory[i - 1]);
	}


	TF1 *one = new TF1("One", "1", 0, 100);
	one->SetLineColorAlpha(kBlack, 0.3);
	one->SetLineStyle(9);


	theory->SetLineColor(kRed);
	theory->SetMarkerColor(kRed);
	theory->SetLineWidth(6);
	theory->SetLineStyle(kDashed);
	theory->GetYaxis()->SetNdivisions(504);
	theory->GetYaxis()->SetTitleOffset(1.3);

	theory->GetYaxis()->SetRangeUser(0.7, 2.2);
	theory->GetXaxis()->SetTitle("#DeltaR");
	theory->GetYaxis()->SetTitle("Au+Au (0-10\%)/p+p");

	theory->Draw();
	one->Draw("SAME");

	auto legend = new TLegend(0.90,0.67,0.70,0.90);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->SetTextFont(132);
    legend->AddEntry(theory,"CCNU","l");

    legend->Draw("SAME");

    TPaveText *data = new TPaveText(0.25,0.60,0.55,0.86, "NDC");
	data->SetFillColor(0);
	data->SetTextAlign(12);
	data->SetTextSize(0.04);
	data->SetBorderSize(0);
	auto system1 = data->AddText("Au+Au #sqrt{#it{s}_{NN}} = 200 GeV");
	auto system2 = data->AddText("p+p #sqrt{#it{s}} = 200 GeV");
	data->AddText("anti-#it{k}_{T}, #it{R} = 0.4");
	data->AddText("#it{p}_{T,D^{0}} > 5 GeV/#it{c}");


	// system1->SetTextSize(0.035);
	

	data->Draw("SAME");

    gPad->SaveAs("CCNU_pp200GeV.pdf");
}

void PreliminaryPlots(){
	// JetPtSpectra();
	// DeltaR();
	// RCP();
	DeltaRRatio();
	// InvariantMassComparison();
	// SingleParticleEmbeddingPlot();
	// CCNU();
	// RCP1040();
	// ComparisonWithPythia();
	// dRComparisonWithPythia();
}