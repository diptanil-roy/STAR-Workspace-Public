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

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "../RooUnfold/src/RooUnfoldSvd.h"

#pragma link C++ class vector<int> +;

using namespace std;

#endif

#include "../BinDef.h"

TRandom3 *r = new TRandom3(0);

RooUnfoldResponse *MakeResponseMatrixFromPrior(int numberofentries = 5000000, TH1D *JetPt = NULL, TH1D *ZDist = NULL, TH1D *DeltaPhi = NULL, TH1D *JetPtDiff[] = NULL, TH1D *D0PtDiff[] = NULL, TH1D *RecoDeltaPhi = NULL, TString Name = "", double MCD0PtCut = 5.0){
	TH2D *tmpMC = new TH2D("tmpMC", "tmpMC", nBinsJetPtForMCZHist, JetPtBinsForMCZHist, nBinsForMCZHist, ZBinsForMCZHist);
  	TH2D *tmpReco = new TH2D("tmpReco", "tmpReco", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

  	TH1D *tmpMCPt = new TH1D("tmpMCPt", "tmpMCPt", extendedJetPtBins, extendedJetPt);

  	TH1D *tmpMCD0Pt = new TH1D("tmpMCD0Pt", "tmpMCPt", nBinsMCD0Pt, MCD0PtBins);

	RooUnfoldResponse *resp = new RooUnfoldResponse(tmpReco, tmpMC, Name.Data(), Name.Data());

	TH2D *respPt = new TH2D("RespPt", "RespPt", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsJetPtForMCZHist, JetPtBinsForMCZHist);
	TH2D *respZ = new TH2D("RespZ", "RespZ", nBinsForZHist, ZBinsForZHist, nBinsForMCZHist, ZBinsForMCZHist);
	TH2D *respD0Pt = new TH2D("RespD0Pt", "RespD0Pt", nBinsRecoD0Pt, RecoD0PtBins, nBinsMCD0Pt, MCD0PtBins);

	TH1D *hMCJetPt = new TH1D("hMCJetPt", "hMCJetPt", nBinsJetPtForMCZHist, JetPtBinsForMCZHist);
	TH1D *hRecoJetPt = new TH1D("hRecoJetPt", "hRecoJetPt", nBinsJetPtForZHist, JetPtBinsForZHist);

	TH1D *hMCD0Pt = new TH1D("hMCD0Pt", "hMCD0Pt", nBinsMCD0Pt, MCD0PtBins);
	TH1D *hRecoD0Pt = new TH1D("hRecoD0Pt", "hRecoD0Pt", nBinsRecoD0Pt, RecoD0PtBins);

	TH1D *hMCZ = new TH1D("hMCZ", "hMCZ", nBinsForMCZHist, ZBinsForMCZHist);
	TH1D *hRecoZ = new TH1D("hRecoZ", "hRecoZ", nBinsForZHist, ZBinsForZHist);

	TH1D *hDeltaPhi = new TH1D("hDeltaPhi", "hDeltaPhi", 1000, -0.1, 0.5);

	int count = 0;

	while (count < numberofentries){
		double mcjetpt = JetPt->GetRandom(r);
		double mcz = ZDist->GetRandom(r);
		// cout << mcz << endl;
		double maxD0Pt = mcz*mcjetpt/TMath::Cos(0.4); //Because jets are of radius 0.4
		if (maxD0Pt < MCD0PtCut) continue;

		double mcD0Pt = 0.0;
		double deltaphi = -99;

		int howmanyevents = 0;
		while(mcD0Pt < MCD0PtCut){
			deltaphi = r->Uniform(0, 0.4);
			mcD0Pt = mcz*mcjetpt/TMath::Cos(deltaphi);
			howmanyevents++;

			// cout << mcz << "\t" << mcjetpt << "\t" << deltaphi << "\t" << mcD0Pt << endl;
			// if (howmanyevents == 100000) {cout << howmanyevents << " events looked at" << endl; break;}
		}		

		// double mcD0Pt = mcz*mcjetpt;
		if (mcjetpt < 5) continue;
		if (mcD0Pt > mcjetpt) continue;
		if (mcD0Pt < MCD0PtCut || mcD0Pt > 10) continue;

	    int D0HistBin = tmpMCD0Pt->FindBin(mcD0Pt) - 1;
	    int histbin = tmpMCPt->FindBin(mcjetpt) - 1;

	    double recoD0Pt = mcD0Pt + D0PtDiff[D0HistBin]->GetRandom(r);
	    double recoJetPt = mcjetpt + JetPtDiff[D0HistBin]->GetRandom(r);

	    bool isRecoD0Pt = recoD0Pt > MCD0PtCut && recoD0Pt < 10;
	    bool isRecoJetPt = recoJetPt > 0 && recoJetPt < 50; // 78% of jets with pT > 5 GeV is captured within this range

	    if (!isRecoD0Pt || !isRecoJetPt ) continue;

	    double recodeltaphi = RecoDeltaPhi->GetRandom(r);

	    double recoz = recoD0Pt*TMath::Cos(recodeltaphi)/recoJetPt;

	    resp->Fill(recoJetPt, recoz, mcjetpt, mcz);

	    respPt->Fill(recoJetPt, mcjetpt);
	    respZ->Fill(recoz, mcz);
	    respD0Pt->Fill(recoD0Pt, mcD0Pt);

	    hDeltaPhi->Fill(deltaphi);

	    hMCJetPt->Fill(mcjetpt);
	    hMCD0Pt->Fill(mcD0Pt);
	    hMCZ->Fill(mcz);

	    hRecoJetPt->Fill(recoJetPt);
	    hRecoD0Pt->Fill(recoD0Pt);
	    hRecoZ->Fill(recoz);

	    count++;
	}

	TCanvas *c[3];
    for (int i = 0; i < 3; i++){
    	c[i] = new TCanvas(Form("c_%i", i), Form("c_%i", i), 1000, 1000);
    	gPad->SetLogz();
    }
    c[0]->cd();
    respPt->Draw("COLZ");
    c[1]->cd();
    respZ->Draw("COLZ");
    c[2]->cd();
    respD0Pt->Draw("COLZ");

    TCanvas *d[3];
    for (int i = 0; i < 3; i++){
    	d[i] = new TCanvas(Form("d_%i", i), Form("d_%i", i), 1000, 1000);
    	gPad->SetLogy();
    }

    d[0]->cd();
    hRecoJetPt->Draw("SAME");
    hMCJetPt->Draw("SAME");
    JetPt->Draw("SAME");
    SetColor(hRecoJetPt, kRed);
    SetColor(hMCJetPt, kBlack);
    JetPt->Scale(hMCJetPt->Integral()/JetPt->Integral());
    SetColor(JetPt, kGreen-2);

    d[1]->cd();
    hRecoZ->Draw("SAME");
    hMCZ->Draw("SAME");
    ZDist->Draw("SAME");
    SetColor(hRecoZ, kRed);
    SetColor(hMCZ, kBlack);
    SetColor(ZDist, kGreen-2);
    ZDist->Scale(hMCZ->Integral()/ZDist->Integral());


    d[2]->cd();
    hRecoD0Pt->Draw("SAME");
    hMCD0Pt->Draw("SAME");
    SetColor(hRecoD0Pt, kRed);
    SetColor(hMCD0Pt, kBlack);

    TCanvas *e = new TCanvas(Form("e"), Form("e"), 1000, 1000);
    e->cd();
    gPad->SetLogy();
    hDeltaPhi->Draw();

	return resp;
}

void ResponseMaker(){
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();


	TFile *f = new TFile("SmearFactors.root");

	TH1D *MCJetPt;
	TH1D *MCJetZ;
	TH1D *MCDeltaPhi;

	TH1D *MCJetPtDiff[3][nBinsJetPtForMCZHist];
	TH1D *MCD0PtDiff[nBinsMCD0Pt];

	TH1D *RecoDeltaPhi[3];

	MCJetPt = (TH1D *)f->Get("MCPt");
	MCJetZ  = (TH1D *)f->Get("MCZ");
	MCDeltaPhi = (TH1D *)f->Get("DeltaPhi_3"); // We only take the sum total Delta Phi because there is no centrality dependence

	for (int cent = 0; cent < 3; cent++){
		for (int i = 0; i < nBinsJetPtForMCZHist; i++ ){
			MCJetPtDiff[cent][i] = (TH1D *)f->Get(Form("Cent %i, Pt in [%.1f, %.1f]", cent, extendedJetPt[i], extendedJetPt[i+1]));
		}
	}

	for (int i = 0; i < nBinsMCD0Pt; i++){
		MCD0PtDiff[i] = (TH1D *)f->Get(Form("D0 Pt in [%.1f, %.1f]", MCD0PtBins[i], MCD0PtBins[i+1]));
	}

	for (int cent = 0; cent < 3; cent++){
		RecoDeltaPhi[cent] = (TH1D *)f->Get(Form("RecoDeltaPhi_%i", cent));
	}

	TH1D *FlatZ = new TH1D("FlatZ", "FlatZ", nBinsMCZ, ZBins);
	for (int zbins = 1; zbins <= nBinsMCZ; zbins++){
		FlatZ->SetBinContent(zbins, 1.);
	}

	RooUnfoldResponse *resp = (RooUnfoldResponse *)MakeResponseMatrixFromPrior(50000000, MCJetPt, FlatZ, MCDeltaPhi, MCJetPtDiff[0], MCD0PtDiff, RecoDeltaPhi[0], "Central", 1.0);
}

