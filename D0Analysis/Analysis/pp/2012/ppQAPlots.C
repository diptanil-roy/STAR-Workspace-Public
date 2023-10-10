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

void ppQAPlots(){

	TFile *inputfile = new TFile("ppQA.root");
	TDirectory *PlotsDir = (TDirectory *)inputfile->Get("PIDQA");

	TH2F *normalised_invbetavpT_tof_ka = (TH2F *)PlotsDir->Get("normalised_invbetavpT_tof_ka");

	const int nbins = normalised_invbetavpT_tof_ka->GetNbinsX();
	double p[nbins];
	double lowres[nbins];
	double midres[nbins];
	double highres[nbins];

	for (int i = 1; i <= nbins; i++){
		double fres =  0.884 + 0.0174/pow(normalised_invbetavpT_tof_ka->GetXaxis()->GetBinCenter(i) + 0.0839, 4.23);
		double fpos = 0.0316 + 0.00137/pow(normalised_invbetavpT_tof_ka->GetXaxis()->GetBinCenter(i) + 0.101, 6.89);

		p[i - 1] = normalised_invbetavpT_tof_ka->GetXaxis()->GetBinCenter(i);
		lowres[i - 1] = -2.0*fres + fpos;
		midres[i - 1] = fpos;
		highres[i - 1] = 2.0*fres + fpos; 
	}

	TGraph *low = new TGraph(nbins, p, lowres);
	TGraph *mid = new TGraph(nbins, p, midres);
	TGraph *high = new TGraph(nbins, p, highres);

	normalised_invbetavpT_tof_ka->Draw("COLZ");
	normalised_invbetavpT_tof_ka->GetXaxis()->SetRangeUser(0.1, 2);
	gPad->SetLogz();

	low->SetLineColor(2);
	high->SetLineColor(2);
	mid->SetLineColor(1);
	low->Draw("C SAME");
	mid->Draw("C SAME");
	high->Draw("C SAME");
}