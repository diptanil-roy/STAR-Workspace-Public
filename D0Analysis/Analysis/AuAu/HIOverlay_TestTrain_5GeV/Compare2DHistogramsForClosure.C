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
#include "TRatioPlot.h"
#include <vector>

#pragma link C++ class vector<int> +;

#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldSvd.h"

using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

void Compare2DHistogramsForClosure(){
	TFile *f1 = new TFile("Feb13_FONLL_HI_Fixing2D4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_7_Split_1.root");
	TFile *f2 = new TFile("Feb13_FONLL_HI_Fixing2D4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_7_Split_2.root");

	RooUnfoldResponse *resp1 = (RooUnfoldResponse*)f1->Get("Resp_0");
	resp1->SetNameTitle("Resp_0_Split1", "Resp_0_Split1");
	RooUnfoldResponse *resp2 = (RooUnfoldResponse*)f2->Get("Resp_0");
	resp2->SetNameTitle("Resp_0_Split2", "Resp_0_Split2");

	TH2D *Meas1 = (TH2D *)resp1->Hmeasured();
	SetName(Meas1, "f_Meas_Cent_0_Split1");
	TH2D *Meas2 = (TH2D *)resp2->Hmeasured();
	SetName(Meas1, "f_Meas_Cent_0_Split2");

	TH1D *MeasPt1 = (TH1D *)Meas1->ProjectionX();
	TH1D *MeasPt2 = (TH1D *)Meas2->ProjectionX();

	TH1D *MeasZ1 = (TH1D *)Meas1->ProjectionY();
	TH1D *MeasZ2 = (TH1D *)Meas2->ProjectionY();

	SetColor(MeasPt1, kRed, 20);
	SetColor(MeasPt2, kBlue, 24);

	SetColor(MeasZ1, kRed, 20);
	SetColor(MeasZ2, kBlue, 24);

	TCanvas *c = new TCanvas("c", "c", 1000, 1000);
	c->cd();
	MeasPt1->Draw();
	MeasPt2->Draw("SAME");
	gPad->SetLogy();

	TCanvas *d = new TCanvas("d", "d", 1000, 1000);
	d->cd();
	MeasZ1->Draw();
	MeasZ2->Draw("SAME");
	gPad->SetLogy();
}