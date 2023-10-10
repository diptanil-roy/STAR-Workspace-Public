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
// #ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
// #include <cstdlib>
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

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void MakeFONLLHistogram(double lowptcutoff = 5.0, double highptcutoff = 20.0){
	TString FileName = "FONLL.txt";
	TH1D *FONLL = new TH1D(Form("FONLL"), Form("FONLL"), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

	double pt, ptweight;
	ifstream myfile (FileName.Data());
	if (myfile.is_open())
	{
		while ( !myfile.eof() )
		{
			myfile >> pt >> ptweight;
			// cout << x << "\t" << y << endl;
			if (pt < lowptcutoff || pt > highptcutoff) continue;
			int bin = FONLL->FindBin(pt+0.0001);
			FONLL->SetBinContent(bin, ptweight);
		}
	}

	FONLL->Scale(1./FONLL->Integral());
	FONLL->Draw();
	gPad->SetLogy();

	TFile *f = new TFile(Form("FONLL_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "RECREATE");
	f->cd();
	FONLL->Write();
	f->Close();
}
