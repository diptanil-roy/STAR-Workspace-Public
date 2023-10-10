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

void Test(){
	TFile *f1 = new TFile("../SPlotFrameWork/ApplyWeights/Histograms3_D05_10GeV_NewBinningOldCuts.root");
	TFile *f2 = new TFile("Jan26_FONLL_HI4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root");

	TH1D *h1 = (TH1D *)f1->Get("centrality");
	TH1D *h2 = (TH1D *)f2->Get("hCentrality");

	TH1D *k1 = (TH1D *)h2->Clone();
	TH1D *k2 = (TH1D *)h2->Clone();

	double acent0 = h1->Integral(1, 2);
	double acent1 = h1->Integral(3, 8);
	double acent2 = h1->Integral(9, 16);

	double bcent0 = h2->Integral(1, 2);
	double bcent1 = h2->Integral(3, 5);
	double bcent2 = h2->Integral(6, 9);

	k1->SetBinContent(1, h1->GetBinContent(1)/acent0);
	k1->SetBinContent(2, h1->GetBinContent(2)/acent0);
	k1->SetBinContent(3, h1->GetBinContent(3)/acent1);
	k1->SetBinContent(4, h1->GetBinContent(5)/acent1);
	k1->SetBinContent(5, h1->GetBinContent(7)/acent1);
	k1->SetBinContent(6, h1->GetBinContent(9)/acent2);
	k1->SetBinContent(7, h1->GetBinContent(11)/acent2);
	k1->SetBinContent(8, h1->GetBinContent(13)/acent2);
	k1->SetBinContent(9, h1->GetBinContent(15)/acent2);

	k2->SetBinContent(1, h2->GetBinContent(1)/bcent0);
	k2->SetBinContent(2, h2->GetBinContent(2)/bcent0);
	k2->SetBinContent(3, h2->GetBinContent(3)/bcent1);
	k2->SetBinContent(4, h2->GetBinContent(4)/bcent1);
	k2->SetBinContent(5, h2->GetBinContent(5)/bcent1);
	k2->SetBinContent(6, h2->GetBinContent(6)/bcent2);
	k2->SetBinContent(7, h2->GetBinContent(7)/bcent2);
	k2->SetBinContent(8, h2->GetBinContent(8)/bcent2);
	k2->SetBinContent(9, h2->GetBinContent(9)/bcent2);

	SetColor(k1, kRed, 20);
	SetColor(k2, kBlack, 24);

	k1->Draw();
	k2->Draw("SAME");
}