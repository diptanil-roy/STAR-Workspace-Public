#ifndef __CINT__
#include <stdexcept>
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
using namespace std;
#endif

#include "../BinDef.h"

void Tower(){
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	TFile *f = new TFile("/Users/diptanilroy/Downloads/test.root");
	f->cd("TowerMapper");

	TH2F *hTower = (TH2F *)gDirectory->Get("hBEMCEtaPhi");

	TCanvas *c = new TCanvas("c", "c", 6000, 3000);
	c->cd();
	hTower->SetMarkerSize(0.4);
	hTower->Draw("TEXT0");
}