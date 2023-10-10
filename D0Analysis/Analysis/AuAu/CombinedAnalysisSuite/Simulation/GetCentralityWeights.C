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
#include "TRandom3.h"
// #include "StJetTreeStruct.h"
#include <vector>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

#pragma link C++ class vector<int> +;

using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

void GetCentralityWeights(int D0pTLow = 1){
	TH1::SetDefaultSumw2();

	// TFile *f1 = new TFile("CentWeight/Response_Step_0_IterParam_0_njptbin_6_nzbin_6.root");
	TFile *f1 = new TFile(Form("CentWeight_D0pT_%iGeV/Response_Step_0_IterParam_0_njptbin_%i_nzbin_%i.root", D0pTLow, njpt_gen_bins_var, nz_gen_bins));
	// TFile *f1 = new TFile("Dec30_DataHIUnf_Updated_WithCentWeight4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root ");
	TFile *f2 = new TFile(Form("../Data/Jul18_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_%i.root", D0pTLow, D0pTLow, 1000));

	TH1D *HICent = (TH1D *)f1->Get("hCentrality");
	TH1D *DataCent = (TH1D *)f2->Get("centrality");

	TH1D *HIgRef = (TH1D *)f1->Get("hgRefMultCorr");
	TH1D *DatagRef = (TH1D *)f2->Get("JetRefMult");

	HIgRef->Scale(1./HIgRef->Integral());
	DatagRef->Scale(1./DatagRef->Integral());

	const int ncentbin = 9;
	double centbin[ncentbin + 1] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

	TH1D *h1 = (TH1D *)HICent->Rebin(ncentbin, "HI Centrality", centbin);
	TH1D *h2 = (TH1D *)DataCent->Rebin(ncentbin, "Data Centrality", centbin);
	TH1D *h3 = (TH1D *)DataCent->Clone("Cent");
	h1->SetNameTitle("HI_NoWeights", "HI_NoWeights");
	// h1->SetNameTitle("HI_Weights", "HI_Weights");
	h2->SetNameTitle("Data", "Data");
	h3->Scale(1./h3->Integral());

	TCanvas *c = new TCanvas("c", "c", 800, 800);
	c->cd();
	// h1->Scale(1./h1->Integral());
	// h2->Scale(1./h2->Integral());

	TH1D *k1 = (TH1D *)h1->Clone("HI");
	TH1D *k2 = (TH1D *)h2->Clone("Data");

	k1->Scale(1./k1->Integral());
	k2->Scale(1./k2->Integral());

	SetColor(k1, kRed, 20);
	SetColor(k2, kBlack, 24);

	k1->GetYaxis()->SetRangeUser(0, 0.4);

	k1->Draw("HIST");
	k2->Draw("EP SAME");
	// h3->Draw("EP SAME");


	double val9 = DataCent->GetBinContent(DataCent->FindBin(70));
	double val8 = DataCent->GetBinContent(DataCent->FindBin(60));
	double val7 = DataCent->GetBinContent(DataCent->FindBin(50));
	double val6 = DataCent->GetBinContent(DataCent->FindBin(40));
	double val5 = DataCent->GetBinContent(DataCent->FindBin(30));
	double val4 = DataCent->GetBinContent(DataCent->FindBin(20));
	double val3 = DataCent->GetBinContent(DataCent->FindBin(10));
	double val2 = DataCent->GetBinContent(DataCent->FindBin(5));
	double val1 = DataCent->GetBinContent(DataCent->FindBin(0));    

	double sval9 = HICent->GetBinContent(HICent->FindBin(70)); // 70-80
	double sval8 = HICent->GetBinContent(HICent->FindBin(60));
	double sval7 = HICent->GetBinContent(HICent->FindBin(50));
	double sval6 = HICent->GetBinContent(HICent->FindBin(40));
	double sval5 = HICent->GetBinContent(HICent->FindBin(30));
	double sval4 = HICent->GetBinContent(HICent->FindBin(20));
	double sval3 = HICent->GetBinContent(HICent->FindBin(10));
	double sval2 = HICent->GetBinContent(HICent->FindBin(5));
	double sval1 = HICent->GetBinContent(HICent->FindBin(0)); // 0-5



	cout << "w9 " << val9<< endl;
	cout << "w8 " << val8<< endl;
	cout << "w7 " << val7<< endl;
	cout << "w6 " << val6<< endl;
	cout << "w5 " << val5<< endl;
	cout << "w4 " << val4<< endl;
	cout << "w3 " << val3<< endl;
	cout << "w2 " << val2<< endl;
	cout << "w1 " << val1<< endl;




	double norm = val9+val8+val7+val6;
	double snorm = sval9+sval8+sval7+sval6;
	double w9 = (val9/norm)/(sval9/snorm);
	double w8 = (val8/norm)/(sval8/snorm);
	double w7 = (val7/norm)/(sval7/snorm);
	double w6 = (val6/norm)/(sval6/snorm);


	norm= val3+val4+val5;
	snorm = sval3+sval4+sval5;
	double w5 = (val5/norm)/(sval5/snorm);
	double w4 = (val4/norm)/(sval4/snorm);
	double w3 = (val3/norm)/(sval3/snorm);

	norm= val1+val2;
	snorm = sval1+sval2;
	double w2 = (val2/norm)/(sval2/snorm);
	double w1 = (val1/norm)/(sval1/snorm);

	cout << "w9 " << w9<< endl;
	cout << "w8 " << w8<< endl;
	cout << "w7 " << w7<< endl;
	cout << "w6 " << w6<< endl;
	cout << "w5 " << w5<< endl;
	cout << "w4 " << w4<< endl;
	cout << "w3 " << w3<< endl;
	cout << "w2 " << w2<< endl;
	cout << "w1 " << w1<< endl;

	TH1D *Weight = (TH1D *)h2->Clone("Centrality Weights");
	Weight->SetBinContent(1, w1);
	Weight->SetBinContent(2, w2);
	Weight->SetBinContent(3, w3);
	Weight->SetBinContent(4, w4);
	Weight->SetBinContent(5, w5);
	Weight->SetBinContent(6, w6);
	Weight->SetBinContent(7, w7);
	Weight->SetBinContent(8, w8);
	Weight->SetBinContent(9, w9);

	// TH1D *Weight = (TH1D *)k2->Clone("Centrality Weights");
	// Weight->Divide(k1);

	for (int i = 1; i <= Weight->GetNbinsX(); i++){
		Weight->SetBinError(i, 0);
	}

	TH1D *WeightFromRefCorr = (TH1D *)DatagRef->Clone("Centrality Weights From RefCorr");
	SetName(WeightFromRefCorr, "Centrality Weights From RefCorr");
	WeightFromRefCorr->Divide(HIgRef);

	for (int i = 1; i <= WeightFromRefCorr->GetNbinsX(); i++){
		WeightFromRefCorr->SetBinError(i, 0);
	}

	TCanvas *d = new TCanvas("d", "d", 800, 800);
	d->cd();
	Weight->Draw("HIST");

	// TF1 *Glauber = new TF1("Glauber", "[0]*([1]*(x-[2])*(x-[2]))/(1 + exp((x-[3])/[4]))", 0, 600);

	// // TF1 *Glauber = 
	// Glauber->SetParameters(1,1,1,1,1);

	TCanvas *e= new TCanvas("e", "e", 800, 800);
	e->cd();
	WeightFromRefCorr->Draw();
	// WeightFromRefCorr->Fit("Glauber", "LLR");
	

	TFile *f = new TFile(Form("CentWeight_D0pT_%iGeV/CentWeight.root", D0pTLow), "RECREATE");
	f->cd();
	Weight->Write();
	WeightFromRefCorr->Write();
	h1->Write();
	h2->Write();
	f->Close();

	cout << f->GetName() << endl;
}
