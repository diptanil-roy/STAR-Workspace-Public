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

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

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

double Chi2Func(TH1D *h1, TH1D *h2){ //h1 is experimental, h2 is truth
	assert (h1->GetNbinsX() == h2->GetNbinsX());

	double chi2 = 0;

	TH1D *hist1 = (TH1D *)h1->Clone("hist1");
	TH1D *hist2 = (TH1D *)h2->Clone("hist2");

	// cout << "Of course they are the same" << endl;

	hist1->Add(hist2, -1);
	hist1->Multiply(hist1);
	hist1->Divide(hist2);



	for (int bin = 1; bin <= hist1->GetNbinsX(); bin++){
		chi2 += hist1->GetBinContent(bin);
	}
	// cout << " Chi2 = " << chi2 << endl;
	return chi2; 
}

void Method(){

	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

	TFile *fin1 = new TFile("../Files/D0Yield_Feb22.root");
	fin1->cd("D0Tagger");
	TH1F *hCent = (TH1F *)gDirectory->Get("hCentralityWeightedAfterCuts");

	for(int i =1;i<hCent->GetNbinsX()+1;i++)cout << hCent->GetBinContent(i) << endl;
	float nevts_0_10  = hCent->Integral(1, 2);
	float nevts_10_40 = hCent->Integral(3, 8);
	float nevts_40_80 = hCent->Integral(9, 16);

	cout << "NEvents = " << nevts_0_10 << "\t" << nevts_10_40 << "\t" << nevts_40_80 << endl;

	double s[3];
	s[0] = 941.23714*nevts_0_10; //(1048*1.38188e+06+838*1.22039e+06)/(1.38188e+06+1.22039e+06);                                                                                                                              
	s[1] = 391.35550*nevts_10_40;//571 + 351 + 206;                                                                                                                                                                           
	s[2] = 56.62475*nevts_40_80;

	cout << "s = " << s[0] << "\t" << s[1] << "\t" << s[2] << endl;

   double taa[3] = {941.23714, 391.35550, 56.62475};

	int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
	int colors[200];

	for (int i = 0; i < 6; i++){
	    for (int j = -9; j <= 10; j++){
	        colors[i*20 + (j+9)] = col[i] + j; 
	    }
	}

	TH1D *Ratio[3];

	TFile *fin2 = new TFile("Inefficiency.root");
	for (int i = 0; i < 3; i++){
		Ratio[i] = (TH1D *)fin2->Get(Form("Ratio_%i", i));
	}

	// Superiteration need to be one loop for each iteration

	TH1D *MCJetPt[3];
	TH1D *MCJetZ[3];

	const int NLOWITERATIONS = 3;
	const int NITERATIONS = 7;
	const int NSUPERITERATIONS = 40;

	TH1D *PriorMCPt[NITERATIONS][NSUPERITERATIONS + 2][3];
	TH1D *PriorMCZ[NITERATIONS][NSUPERITERATIONS + 2][3];

	TH1D *UnfoldedJetPt[NITERATIONS][NSUPERITERATIONS + 2][3];
	TH1D *UnfoldedJetZ[NITERATIONS][NSUPERITERATIONS + 2][3];

	TH1D *RatioJetPt[NITERATIONS][NSUPERITERATIONS + 2][3];
	TH1D *RatioJetZ[NITERATIONS][NSUPERITERATIONS + 2][3];

	TH1D *Weight[NITERATIONS][NSUPERITERATIONS + 2][3];


	TH1D *VzMC[NITERATIONS][3];

	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
		// cout << "====================================== ITER = " << iter << endl;
		for (int siter = 0; siter <= NSUPERITERATIONS; siter++){
			TFile *f = new TFile(Form("Response_%i_%i.root", siter, iter));
			for (int cent = 0; cent < 3; cent++){

				if (siter == 0) VzMC[iter-3][cent] = (TH1D *)f->Get(Form("pthatbin_0/VzMC_%i", cent));

				if (siter == 0) {
					MCJetPt[cent] = (TH1D *)ProcessSpectraHistogram((TH1D *)f->Get(Form("PriorMCPt_%i", cent)));
					MCJetZ[cent] = (TH1D *)ProcessSpectraHistogram((TH1D *)f->Get(Form("PriorMCZ_%i", cent)));

					MCJetPt[cent]->SetDirectory(0);
					MCJetZ[cent]->SetDirectory(0);

					MCJetPt[cent]->SetNameTitle(Form("PythiaPt_%i", cent), Form("PythiaPt_%i", cent));
					MCJetZ[cent]->SetNameTitle(Form("PythiaZ_%i", cent), Form("PythiaZ_%i", cent));

					MCJetPt[cent]->Scale(1./VzMC[iter-3][cent]->Integral());
					MCJetZ[cent]->Scale(1./VzMC[iter-3][cent]->Integral());
				}

				PriorMCPt[iter-3][siter][cent] = (TH1D *)f->Get(Form("PriorMCPt_%i", cent));
				PriorMCZ[iter-3][siter][cent] = (TH1D *)f->Get(Form("PriorMCZ_%i", cent));

				// UnfoldedJetPt[iter-3][siter][cent] = (TH1D *)f->Get(Form("UnfoldedPt_%i", cent));
				UnfoldedJetPt[iter-3][siter][cent] = (TH1D *)ProcessSpectraHistogram((TH1D *)f->Get(Form("UnfoldedPt_%i", cent)));
				cout << "Dividing here " << UnfoldedJetPt[iter-3][siter][cent]->GetNbinsX() << "\t" << Ratio[cent]->GetNbinsX() << endl;

				// for (int i = 1; i <= UnfoldedJetPt[iter-3][siter][cent]->GetNbinsX(); i++){
				// 	if (Ratio[cent]->GetBinContent(i) != 0) {
				// 		UnfoldedJetPt[iter-3][siter][cent]->SetBinContent(i, UnfoldedJetPt[iter-3][siter][cent]->GetBinContent(i)/Ratio[cent]->GetBinContent(i));
				// 		UnfoldedJetPt[iter-3][siter][cent]->SetBinError(i, UnfoldedJetPt[iter-3][siter][cent]->GetBinError(i)/Ratio[cent]->GetBinContent(i));
				// 	}
				// 	else UnfoldedJetPt[iter-3][siter][cent]->SetBinContent(i,0);
				// }
				// UnfoldedJetPt[iter-3][siter][cent]->Divide(Ratio[cent]);
				// UnfoldedJetPt[iter-3][siter][cent]->Scale(1./VzMC[iter-3][cent]->Integral());
				UnfoldedJetPt[iter-3][siter][cent]->Scale(1./s[cent]);
				UnfoldedJetPt[iter-3][siter][cent]->SetNameTitle(Form("UnfoldedPt_%i_%i_%i", iter, siter, cent), Form("UnfoldedPt_%i_%i_%i", iter, siter, cent));
				
				// UnfoldedJetZ[iter-3][siter][cent] = (TH1D *)f->Get(Form("UnfoldedZ_%i", cent));
				UnfoldedJetZ[iter-3][siter][cent] = (TH1D *)ProcessSpectraHistogram((TH1D *)f->Get(Form("UnfoldedZ_%i", cent)));
				// UnfoldedJetZ[iter-3][siter][cent]->Divide(Ratio[cent]);
				// UnfoldedJetZ[iter-3][siter][cent]->Scale(1./VzMC[iter-3][cent]->Integral());
				UnfoldedJetZ[iter-3][siter][cent]->Scale(1./s[cent]);
				UnfoldedJetZ[iter-3][siter][cent]->SetNameTitle(Form("UnfoldedZ_%i_%i_%i", iter, siter, cent), Form("UnfoldedZ_%i_%i_%i", iter, siter, cent));


				Weight[iter-3][siter][cent] = (TH1D *)f->Get(Form("Weight_%i", cent));
				if (!Weight[iter-3][siter][cent]) cout << "WTF" << endl;
				Weight[iter-3][siter][cent]->SetNameTitle(Form("Weight_%i_%i_%i", iter, siter, cent), Form("Weight_%i_%i_%i", iter, siter, cent));

				cout << iter << "\t" << siter << "\t" << cent << "\t" << Weight[iter-3][siter][cent]->GetName() << endl;
				// UnfoldedJetPt[iter-3][siter][cent]->Scale(1./taa[cent]);
				// UnfoldedJetZ[iter-3][siter][cent]->Scale(1./taa[cent]);

				// cout << iter << "\t" << siter << "\t" << UnfoldedJetZ[iter-3][siter][cent]->GetName() << endl;

				PriorMCPt[iter-3][siter][cent]->SetDirectory(0);
				PriorMCZ[iter-3][siter][cent]->SetDirectory(0);
				UnfoldedJetPt[iter-3][siter][cent]->SetDirectory(0);
				UnfoldedJetZ[iter-3][siter][cent]->SetDirectory(0);
				Weight[iter-3][siter][cent]->SetDirectory(0);
				VzMC[iter-3][cent]->SetDirectory(0);

				// cout << iter << "\t" << siter << "\t" << UnfoldedJetZ[iter-3][siter][cent]->GetName() << endl;

				// RatioJetPt[iter-3][siter][cent]->SetDirectory(0);
				// RatioJetZ[iter-3][siter][cent]->SetDirectory(0);

				// UnfoldedJetPt[iter-3][siter][cent]->Add(MCJetPt[cent], -1);
				// UnfoldedJetZ[iter-3][siter][cent]->Add(MCJetZ[cent], -1);
			}
			f->Close();
			// cout << iter << "\t" << siter << "\t" << UnfoldedJetZ[iter-3][siter][0]->GetName() << endl;
			// cout << iter << "\t" << siter << "\t" << UnfoldedJetZ[iter-3][siter][1]->GetName() << endl;
			// cout << iter << "\t" << siter << "\t" << UnfoldedJetZ[iter-3][siter][2]->GetName() << endl;
		}
	}
	
	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
		for (int siter = 0; siter <= NSUPERITERATIONS; siter++){
			RatioJetPt[iter-3][siter][0] = (TH1D *)UnfoldedJetPt[iter-3][siter][0]->Clone();
			RatioJetPt[iter-3][siter][1] = (TH1D *)UnfoldedJetPt[iter-3][siter][1]->Clone();
			RatioJetPt[iter-3][siter][2] = (TH1D *)UnfoldedJetPt[iter-3][siter][0]->Clone();
			RatioJetPt[iter-3][siter][0]->Divide(UnfoldedJetPt[iter-3][siter][2]);
			RatioJetPt[iter-3][siter][1]->Divide(UnfoldedJetPt[iter-3][siter][2]);
			RatioJetPt[iter-3][siter][2]->Divide(MCJetPt[0]);

			RatioJetPt[iter-3][siter][0]->SetNameTitle(Form("RCP_%i_%i_%i", iter, siter, 0), Form("RCP_%i_%i_%i", iter, siter, 0));
			RatioJetPt[iter-3][siter][1]->SetNameTitle(Form("RCP_%i_%i_%i", iter, siter, 1), Form("RCP_%i_%i_%i", iter, siter, 1));
			RatioJetPt[iter-3][siter][2]->SetNameTitle(Form("RCP_%i_%i_%i", iter, siter, 2), Form("RCP_%i_%i_%i", iter, siter, 2));
		}
	}

	/*
	
	cout << "This is fine" << endl;

	TCanvas *c[3];
	for (int cent = 0; cent < 3; cent++){
		c[cent] = new TCanvas (Form("AllZPlots_%i", cent), Form("AllZPlots_%i", cent), 2000, 1000);
		c[cent]->Divide(6,3);
		for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
			// cout << "====================================== ITER = " << iter << endl;
			c[cent]->cd(iter - 3 + 1);
			gPad->SetLogy();
			for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
				// cout << iter << "\t" << siter << "\t" << UnfoldedJetZ[iter-3][siter][cent]->GetName() << endl;
				UnfoldedJetZ[iter-3][siter][cent]->GetXaxis()->SetRangeUser(0,1);
				// UnfoldedJetZ[iter-3][siter][cent]->GetYaxis()->SetNameTitle(Form("#frac{1}{N} #frac{dN}{dz}"), Form("#frac{1}{N} #frac{dN}{dz}"));
				UnfoldedJetZ[iter-3][siter][cent]->GetXaxis()->SetNameTitle(Form("z"), Form("z"));
				UnfoldedJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter*3],1);
            	// else UnfoldedJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter/2],0.5);
            	UnfoldedJetZ[iter-3][siter][cent]->Draw("L SAME");
            	if (siter == NSUPERITERATIONS)UnfoldedJetZ[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));

            }
			// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
            // MCJetZ[cent]->GetXaxis()->SetRangeUser(0,1);
            // MCJetZ[cent]->SetLineColor(kBlack);
            // MCJetZ[cent]->Draw("HIST SAME");
            gPad->BuildLegend();
        }
    }

    cout << "This is fine" << endl;
 //    TCanvas *c2[3];
	// for (int cent = 0; cent < 3; cent++){
	// 	c2[cent] = new TCanvas (Form("RatioZPlots_%i", cent), Form("RatioZPlots_%i", cent), 2000, 1000);
	// 	c2[cent]->Divide(6,3);
	// 	for (int iter = 3; iter <= NITERATIONS; iter++){
	// 		c2[cent]->cd(iter - 3 + 1);
	// 		gPad->SetLogy();
	// 		for (int siter = NSUPERITERATIONS; siter >= NSUPERITERATIONS-2; siter--){
	// 			RatioJetZ[iter-3][siter][cent]->GetXaxis()->SetRangeUser(0,1);
	// 			RatioJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter],1);
 //            	// else RatioJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter/2],0.5);
 //            	if(siter == 0)RatioJetZ[iter-3][siter][cent]->SetLineColorAlpha(kBlack,1);
 //            	RatioJetZ[iter-3][siter][cent]->Draw("L SAME");
 //            }
 //            // MCJetZ[cent]->SetLineColor(kBlack);
 //            // MCJetZ[cent]->Draw("HIST SAME");
 //            gPad->BuildLegend();
 //        }
 //    }
    

   TCanvas *d[3];
	for (int cent = 0; cent < 3; cent++){
		d[cent] = new TCanvas (Form("AllPtPlots_%i", cent), Form("AllPtPlots_%i", cent), 2000, 1000);
		d[cent]->Divide(6,3);
		for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
			d[cent]->cd(iter - 3 + 1);
			gPad->SetLogy();
			for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
				UnfoldedJetPt[iter-3][siter][cent]->GetXaxis()->SetRangeUser(1, 30);
				// UnfoldedJetPt[iter-3][siter][cent]->GetYaxis()->SetNameTitle(Form("#frac{1}{N} #frac{dN}{dp_{T,jet}}"), Form("#frac{1}{N} #frac{dN}{dp_{T,jet}}"));
				UnfoldedJetPt[iter-3][siter][cent]->GetXaxis()->SetNameTitle(Form("p_{T,jet}"), Form("p_{T,jet}"));
				UnfoldedJetPt[iter-3][siter][cent]->GetYaxis()->SetRangeUser(1.0*pow(10, -3), pow(10, 7));
				UnfoldedJetPt[iter-3][siter][cent]->SetLineColorAlpha(colors[siter*3],1);
            	UnfoldedJetPt[iter-3][siter][cent]->Draw("L SAME");
            	if (siter == NSUPERITERATIONS)UnfoldedJetPt[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
            }
			// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
            // MCJetPt[cent]->GetXaxis()->SetRangeUser(1,30);
            // MCJetPt[cent]->SetLineColor(kBlack);
            // MCJetPt[cent]->Draw("HIST SAME");
        }
    }

    cout << "This is fine 3" << endl;
    
 //    TCanvas *d2[3];
	// for (int cent = 0; cent < 3; cent++){
	// 	d2[cent] = new TCanvas (Form("RatioPtPlots_%i", cent), Form("RatioPtPlots_%i", cent), 2000, 1000);
	// 	d2[cent]->Divide(6,3);
	// 	for (int iter = 3; iter <= NITERATIONS; iter++){
	// 		d2[cent]->cd(iter - 3 + 1);
	// 		gPad->SetLogy();
	// 		for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
	// 			RatioJetPt[iter-3][siter][cent]->GetXaxis()->SetRangeUser(5, 20);
 //            	RatioJetPt[iter-3][siter][cent]->Draw("L SAME");
 //            }
 //            // MCJetPt[cent]->SetLineColor(kBlack);
 //            // MCJetPt[cent]->Draw("HIST SAME");
 //        }
 //    }
    
    double low = -0.5;
    double high = low + 100.5;

    TH1D *Chi2ForJetPt[NITERATIONS+1][3];
    TH1D *Chi2ForJetZ[NITERATIONS+1][3];

    for (int cent = 0; cent < 3; cent++){
    	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
    		
        	Chi2ForJetPt[iter-3][cent] = new TH1D(Form("Chi2ForJetPt_%i_%i", iter, cent), Form("Chi2ForJetPt_%i_%i", iter, cent), 101, low, high);
        	Chi2ForJetZ[iter-3][cent] = new TH1D(Form("Chi2ForJetZ_%i_%i", iter, cent), Form("Chi2ForJetZ_%i_%i", iter, cent), 101, low, high);
        	
        	cout << iter << "\t" << endl;

        	for (int siter = 0; siter <= NSUPERITERATIONS; siter++){
        		double chi2ptvalue = MCJetPt[cent]->Chi2Test(UnfoldedJetPt[iter-3][siter][cent], "WW");
        		
        		double chi2ptvalue_mine = Chi2Func(UnfoldedJetPt[iter-3][siter][cent], MCJetPt[cent]);

        		cout << chi2ptvalue_mine << endl;

        		Chi2ForJetPt[iter-3][cent]->SetBinContent(siter+1, chi2ptvalue_mine);

        		double chi2zvalue = MCJetZ[cent]->Chi2Test(UnfoldedJetZ[iter-3][siter][cent], "WW");
        		double chi2zvalue_mine = Chi2Func(UnfoldedJetZ[iter-3][siter][cent], MCJetZ[cent]);

        		Chi2ForJetZ[iter-3][cent]->SetBinContent(siter+1, chi2zvalue_mine);

        		if (siter == 0 || siter == NSUPERITERATIONS-1 || siter == NSUPERITERATIONS) cout << cent << "\t" << siter << "\t" << iter << "\t" << UnfoldedJetPt[iter-3][siter][cent]->GetName() << "\t" << chi2ptvalue_mine << "\t" << UnfoldedJetZ[iter-3][siter][cent]->GetName() << "\t" << chi2zvalue_mine << endl;
        	}
        }
    }

    cout << "This is fine 5" << endl;

    TCanvas *e = new TCanvas("ChiSq", "ChiSq", 2000, 1000);
    e->Divide(3,2);

    for (int cent = 0; cent < 3; cent++){
    	e->cd(cent + 1);
    	gPad->SetLogy();
    	auto legend0 = new TLegend(0.6,0.7,0.85,0.9);
    	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
    		Chi2ForJetPt[iter-3][cent]->GetXaxis()->SetTitle("#SI");
    		Chi2ForJetPt[iter-3][cent]->GetYaxis()->SetTitle("#chi^{2}");
    		Chi2ForJetPt[iter-3][cent]->GetXaxis()->SetRangeUser(0, 40);
    		Chi2ForJetPt[iter-3][cent]->GetYaxis()->SetRangeUser(5*pow(10, -6), 5*pow(10, -4));
    		Chi2ForJetPt[iter-3][cent]->SetLineColor(colors[iter*3]);
    		Chi2ForJetPt[iter-3][cent]->Draw("HIST SAME");
    		legend0->AddEntry(Chi2ForJetPt[iter-3][cent],Form("Iteration %i", iter),"l");
    	}

        legend0->Draw("SAME");

    	e->cd(cent + 4);
    	gPad->SetLogy();
    	auto legend1 = new TLegend(0.6,0.7,0.85,0.9);
    	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
    		Chi2ForJetZ[iter-3][cent]->GetXaxis()->SetTitle("#SI");
    		Chi2ForJetZ[iter-3][cent]->GetYaxis()->SetTitle("#chi^{2}");
    		Chi2ForJetZ[iter-3][cent]->GetXaxis()->SetRangeUser(0, 40);
    		Chi2ForJetZ[iter-3][cent]->GetYaxis()->SetRangeUser(pow(10, -2), pow(10,2));
    		Chi2ForJetZ[iter-3][cent]->SetLineColor(colors[iter*3]);
    		Chi2ForJetZ[iter-3][cent]->Draw("HIST SAME");
    		legend1->AddEntry(Chi2ForJetZ[iter-3][cent],Form("Iteration %i", iter),"l");
    	}

        legend1->Draw("SAME");
    }

    cout << "This is fine" << endl;

   TCanvas *d3;
	d3 = new TCanvas (Form("RCP_%i", 0), Form("RCP_%i", 0), 2000, 1000);
	d3->Divide(5);
	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
		// d3->cd();
		d3->cd(iter - 3 + 1);
		// gPad->SetLogy();
		for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
			RatioJetPt[iter-3][siter][0]->GetXaxis()->SetRangeUser(5, 20);
			RatioJetPt[iter-3][siter][0]->GetYaxis()->SetRangeUser(-0.5, 5);
			RatioJetPt[iter-3][siter][0]->SetLineColorAlpha(colors[siter*3],1);
        	RatioJetPt[iter-3][siter][0]->Draw("L SAME");
        	if (siter == NSUPERITERATIONS)RatioJetPt[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
        }
        gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
        gPad->BuildLegend();
   }

   cout << "This is fine" << endl;

   TCanvas *d4;
	d4 = new TCanvas (Form("RCP_%i", 2), Form("RCP_%i", 2), 2000, 1000);
	d4->Divide(5);
	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
		// d3->cd();
		d4->cd(iter - 3 + 1);
		// gPad->SetLogy();
		for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
			RatioJetPt[iter-3][siter][2]->GetXaxis()->SetRangeUser(5, 20);
			RatioJetPt[iter-3][siter][2]->GetYaxis()->SetRangeUser(0.5,5);
			RatioJetPt[iter-3][siter][2]->SetLineColorAlpha(colors[siter*3],1);
        	RatioJetPt[iter-3][siter][2]->Draw("L SAME");
        	cout << "Integral = " << RatioJetPt[iter-3][siter][2]->Integral() << endl;
        	if (siter == NSUPERITERATIONS)RatioJetPt[iter-3][siter][2]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
        }
        gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
        gPad->BuildLegend();
   }

   cout << "This is fine" << endl;

 //    TCanvas *d2;
	// d2 = new TCanvas (Form("RCP_%i", 0), Form("RCP_%i", 0), 2000, 1000);
	// for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
	// 	RatioJetPt[4][siter][0]->GetXaxis()->SetRangeUser(5, 20);
	// 	RatioJetPt[4][siter][0]->GetYaxis()->SetRangeUser(-0.5, 5);
	// 	RatioJetPt[4][siter][0]->SetNameTitle(Form("SI #%i", siter), Form("SI #%i", siter));
	// 	RatioJetPt[4][siter][0]->SetLineColorAlpha(colors[siter],1);
 //    	RatioJetPt[4][siter][0]->Draw("L SAME");
 //    }
 //    gPad->BuildLegend();

    cout << "============================================= WEIGHTS ===================================================" << endl;
    TCanvas *k[3];
	for (int cent = 0; cent < 3; cent++){
		k[cent] = new TCanvas (Form("AllWeights_%i", cent), Form("AllWeights_%i", cent), 2000, 1000);
		k[cent]->Divide(6,3);
		for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
			// cout << "====================================== ITER = " << iter << endl;
			k[cent]->cd(iter - 3 + 1);
			// k[cent]->cd();
			gPad->SetLogy();
			for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
				cout << iter << "\t" << siter << "\t" << cent << "\t" << Weight[iter-3][siter][cent]->GetName() << endl;
				Weight[iter-3][siter][cent]->GetXaxis()->SetRangeUser(0,1);
				// UnfoldedJetZ[iter-3][siter][cent]->GetYaxis()->SetNameTitle(Form("#frac{1}{N} #frac{dN}{dz}"), Form("#frac{1}{N} #frac{dN}{dz}"));
				Weight[iter-3][siter][cent]->GetXaxis()->SetNameTitle(Form("z"), Form("z"));
				Weight[iter-3][siter][cent]->SetLineColorAlpha(colors[siter*3],1);
            	// else UnfoldedJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter/2],0.5);
            	Weight[iter-3][siter][cent]->Draw("L SAME");
            	if (siter == NSUPERITERATIONS)Weight[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));

            }
			// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
            // MCJetZ[cent]->GetXaxis()->SetRangeUser(0,1);
            // MCJetZ[cent]->SetLineColor(kBlack);
            // MCJetZ[cent]->Draw("HIST SAME");
            gPad->BuildLegend();
        }
    }

    */

    /*

    cout << "============================================= WEIGHTS ===================================================" << endl;
    TCanvas *k2[3];
	for (int cent = 0; cent < 1; cent++){
		k2[cent] = new TCanvas (Form("WeightsvsZ_%i", cent), Form("WeightsvsZ_%i", cent), 2000, 1000);
		k2[cent]->Divide(10, 4);
		for (int iter = 3; iter <= 3; iter++){
			// cout << "====================================== ITER = " << iter << endl;
			for (int siter = 1; siter <= NSUPERITERATIONS; siter++){
				k2[cent]->cd(NSUPERITERATIONS - siter + 1);
				gPad->SetLogy();
				cout << iter << "\t" << siter << "\t" << cent << "\t" << Weight[iter-3][siter][cent]->GetName() << endl;

				Weight[iter-3][siter][cent]->SetLineColor(kBlack);
				UnfoldedJetZ[iter-3][siter-1][cent]->SetLineColor(kBlue);
				Weight[iter-3][siter][cent]->GetYaxis()->SetRangeUser(pow(10, -20), pow(10,5));
            	Weight[iter-3][siter][cent]->Draw("L SAME");
            	UnfoldedJetZ[iter-3][siter-1][cent]->Divide(Weight[iter-3][siter][cent]);
            	UnfoldedJetZ[iter-3][siter-1][cent]->Draw("L SAME");

            	if (siter == NSUPERITERATIONS)Weight[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));

            }
			// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
            // MCJetZ[cent]->GetXaxis()->SetRangeUser(0,1);
            // MCJetZ[cent]->SetLineColor(kBlack);
            // MCJetZ[cent]->Draw("HIST SAME");
            gPad->BuildLegend();
        }
    }

    TCanvas *k3[3];
	for (int cent = 0; cent < 1; cent++){
		k3[cent] = new TCanvas (Form("MCZ_%i", cent), Form("MCZ_%i", cent), 2000, 1000);
		// k3[cent]->Divide(1);
		for (int iter = 3; iter <= 3; iter++){
			// cout << "====================================== ITER = " << iter << endl;
			for (int siter = 1; siter <= NSUPERITERATIONS; siter++){
				// k3[cent]->cd(NSUPERITERATIONS - siter + 1);
				gPad->SetLogy();

            	UnfoldedJetZ[iter-3][siter-1][cent]->Draw("L SAME");
            	UnfoldedJetZ[iter-3][siter-1][cent]->SetLineColorAlpha(colors[siter*10],1);
            }
			// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
            // MCJetZ[cent]->GetXaxis()->SetRangeUser(0,1);
            MCJetZ[cent]->SetLineColor(kBlack);
            MCJetZ[cent]->Draw("HIST SAME");
            gPad->BuildLegend();
        }
    }

    TCanvas *k4[3];
	for (int cent = 0; cent < 1; cent++){
		k4[cent] = new TCanvas (Form("TestingEachStep_%i", cent), Form("TestingEachStep_%i", cent), 2000, 1000);
		k4[cent]->Divide(3);
		for (int iter = 3; iter <= 3; iter++){
			
			for (int siter = 1; siter <= 3; siter++){
				k4[cent]->cd(siter);
				gPad->SetLogy();

				PriorMCZ[iter-3][siter][cent]->Draw("L SAME");
				PriorMCZ[iter-3][siter][cent]->SetLineColor(kRed);

				PriorMCZ[iter-3][siter][cent]->GetXaxis()->SetRangeUser(0,1);
				PriorMCZ[iter-3][siter][cent]->GetYaxis()->SetRangeUser(pow(10, -7), pow(10, -3));

				MCJetZ[cent]->SetLineColor(kBlack);
            	MCJetZ[cent]->Draw("HIST SAME");
            
            	UnfoldedJetZ[iter-3][siter][cent]->Draw("L SAME");
            	UnfoldedJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter*10],1);

            	cout << "SITER = " << siter << "\t" << MCJetZ[cent]->Integral() << "\t" << PriorMCZ[iter-3][siter][cent]->Integral() << "\t" << UnfoldedJetZ[iter-3][siter][cent]->Integral() << endl;

            	gPad->BuildLegend();
            }
			// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
            // MCJetZ[cent]->GetXaxis()->SetRangeUser(0,1);
        }
    }

    */

    TCanvas *p1 = new TCanvas("p1", "p1", 600, 600);
    p1->cd();
    TH1D *h1 = (TH1D *)PriorMCZ[0][0][0]->Clone();
    TH1D *h2 = (TH1D *)PriorMCZ[0][1][0]->Clone();
    TH1D *h3 = (TH1D *)PriorMCZ[0][2][0]->Clone();

    SetName(h1, "PYTHIA z");
    SetName(h2, "Flat z");
    SetName(h3, "z close to folded data");

    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
    h3->Scale(1./h3->Integral());

    h1->GetYaxis()->SetRangeUser(0,0.5);
    h1->GetXaxis()->SetRangeUser(0,1);
    h1->Draw("LEP");
    h1->SetLineColor(kBlack);
    h1->SetMarkerColor(kBlack);
    h2->Draw("LEP SAME");
    h2->SetLineColor(kBlue);
    h3->Draw("LEP SAME");
    h3->SetLineColor(kRed);

    gPad->BuildLegend();

    // delete h1;
    // delete h2;
    // delete h3;

    delete p1;


    TCanvas *p2[3];

    int threecolpallete[3] = {kBlack, kBlue, kRed};

    for (int siter = 0; siter < 3; siter++){
    	p2[siter] = new TCanvas(Form("p2_%i", siter), Form("p2_%i", siter), 600, 600);
    	gPad->SetLogy();
    	gPad->SetLeftMargin(0.15);
    	UnfoldedJetPt[0][siter][0]->Draw("LEP");
    	UnfoldedJetPt[0][siter][2]->Draw("LEP SAME");

    	SetName(UnfoldedJetPt[0][siter][0], "Central");
    	SetName(UnfoldedJetPt[0][siter][2], "Peripheral");

    	UnfoldedJetPt[0][siter][0]->GetXaxis()->SetRangeUser(5,15);
    	UnfoldedJetPt[0][siter][0]->GetYaxis()->SetRangeUser(2*pow(10, -13), 5*pow(10, -5));

    	UnfoldedJetPt[0][siter][0]->GetXaxis()->SetTitle("p_{T,Jet}");
    	UnfoldedJetPt[0][siter][0]->GetYaxis()->SetTitle("#frac{1}{N_{Jets}}#frac{dN_{Jets}}{dp_{T}}");

    	UnfoldedJetPt[0][siter][0]->SetLineColor(threecolpallete[siter]);
    	UnfoldedJetPt[0][siter][0]->SetMarkerColor(threecolpallete[siter]);

    	UnfoldedJetPt[0][siter][2]->SetLineColor(threecolpallete[siter]);
    	UnfoldedJetPt[0][siter][2]->SetMarkerColor(threecolpallete[siter]);

    	UnfoldedJetPt[0][siter][0]->SetMarkerStyle(20);
    	UnfoldedJetPt[0][siter][2]->SetMarkerStyle(24);

    	gPad->BuildLegend();

    	p2[siter]->SaveAs(Form("JetPt_%i.pdf", siter));

    	delete p2[siter];
    }

    TCanvas *p3 = new TCanvas("p3", "p3", 600, 600);
    p3->cd();
    RatioJetPt[0][0][0]->Draw("LEP");
    RatioJetPt[0][1][0]->Draw("LEP SAME");
    RatioJetPt[0][2][0]->Draw("LEP SAME");

    SetName(RatioJetPt[0][0][0], "PYTHIA z");
    SetName(RatioJetPt[0][1][0], "Flat z");
    SetName(RatioJetPt[0][2][0], "z close to folded data");

    RatioJetPt[0][0][0]->GetXaxis()->SetTitle("p_{T,Jet}");
    RatioJetPt[0][0][0]->GetYaxis()->SetTitle("R_{CP}");

    RatioJetPt[0][0][0]->GetXaxis()->SetRangeUser(5,15);
    RatioJetPt[0][0][0]->GetYaxis()->SetRangeUser(0,5);

    for (int siter = 0; siter < 3; siter++){
    	RatioJetPt[0][siter][0]->SetLineColor(threecolpallete[siter]);
    	RatioJetPt[0][siter][0]->SetMarkerColor(threecolpallete[siter]);
    	RatioJetPt[0][siter][0]->SetMarkerStyle(20);
    }

    gPad->BuildLegend();

    p3->SaveAs("RCP.pdf");

    delete p3;

    TCanvas *p4 = new TCanvas("p4", "p4", 1200, 600);
    p4->Divide(2);

    h1 = (TH1D *)PriorMCZ[0][0][0]->Clone();
    TH1D *k1 = (TH1D *)PriorMCPt[0][0][0]->Clone();

    h1->Scale(1./h1->Integral());
    k1->Scale(1./k1->Integral());

    h2 = (TH1D *)PriorMCZ[0][1][0]->Clone();
    TH1D *k2 = (TH1D *)PriorMCPt[0][1][0]->Clone();

    h2->Scale(1./h2->Integral());
    k2->Scale(1./k2->Integral());

    SetName(h1, "PYTHIA Z");
    SetName(h2, "Flat Z");

    SetName(k1, "PYTHIA Z");
    SetName(k2, "Flat Z");

    SetColor(h1, kBlack);
    SetColor(k1, kBlack);

    SetColor(h2, kRed);
    SetColor(k2, kRed);

    h1->GetXaxis()->SetRangeUser(0,1);
    h1->GetYaxis()->SetRangeUser(2*pow(10, -3), 0.9);

    h1->GetXaxis()->SetTitle("Z");
    h1->GetYaxis()->SetTitle("arb. units");

    k1->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
    k1->GetYaxis()->SetTitle("arb. units");

    k1->GetXaxis()->SetRangeUser(5,15);
    k1->GetYaxis()->SetRangeUser(2*pow(10, -4), 1.1);

    p4->cd(1);
    gPad->SetLogy();

    h1->Draw();
    h2->Draw("SAME");

    gPad->BuildLegend();

    p4->cd(2);
    gPad->SetLogy();

    k1->Draw();
    k2->Draw("SAME");

    gPad->BuildLegend();

    delete p4;

   TFile *OriginalResponseFile = new TFile("../SPlotFrameWork/ApplyWeights/Histograms3_D01_10GeV_JetPt0_20_With2DJetPtZ.root");
  	OriginalResponseFile->cd();

  	TH2D *FoldedData;

  	FoldedData = (TH2D *)gDirectory->Get("ZPt_0_10");

  	FoldedData->Scale(1./s[0]);

  	TCanvas *p5 = new TCanvas("p5", "p5", 1200, 600);
   p5->Divide(2);

   p5->cd(2);
   gPad->SetLogy();
   gPad->SetLeftMargin(0.12);

  	h1 = (TH1D *)UnfoldedJetPt[0][2][0]->Clone();
  	SetColor(h1, kRed);
  	SetName(h1, "Unfolded");
  	h1->GetYaxis()->SetTitle("arb. units");
  	h1->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");

  	h2 = (TH1D *)ProcessSpectraHistogram((TH1D *)FoldedData->ProjectionX());
  	SetColor(h2, kBlack);
  	SetName(h2, "Folded");

  	h1->GetXaxis()->SetRangeUser(0, 20);

  	h1->Draw();
  	h2->Draw("SAME");

  	gPad->BuildLegend();

  	p5->cd(1);
   gPad->SetLogy();
   gPad->SetLeftMargin(0.12);

  	k1 = PriorMCZ[0][3][0];
  	k1->Scale(1/s[0]);
  	SetColor(k1, kRed);
  	SetName(k1, "Unfolded");
  	k1->GetYaxis()->SetTitle("arb. units");
  	k1->GetXaxis()->SetTitle("Z");

  	k2 = (TH1D *)FoldedData->ProjectionY();
  	SetColor(k2, kBlack);
  	SetName(k2, "Folded");

  	k1->GetXaxis()->SetRangeUser(0, 1);

  	k1->Draw();
  	k2->Draw("SAME");

  	gPad->BuildLegend();

  	// delete p5;
  	
  	// TCanvas *p6 = new TCanvas("p6", "p6", 1200, 600);
   //  p6->Divide(2);

   //  h1 = (TH1D *)PriorMCZ[0][0][0]->Clone();
   //  k1 = (TH1D *)PriorMCPt[0][0][0]->Clone();

   //  h1->Scale(1./h1->Integral());
   //  k1->Scale(1./k1->Integral());

   //  h2 = (TH1D *)PriorMCZ[0][1][0]->Clone();
   //  k2 = (TH1D *)PriorMCPt[0][1][0]->Clone();

   //  h2->Scale(1./h2->Integral());
   //  k2->Scale(1./k2->Integral());

   //  h3 = (TH1D *)PriorMCZ[0][2][0]->Clone();
   //  TH1D *k3 = (TH1D *)PriorMCPt[0][2][0]->Clone();

   //  h3->Scale(1./h3->Integral());
   //  k3->Scale(1./k3->Integral());

   //  SetName(h1, "PYTHIA Z");
   //  SetName(h2, "Flat Z");
   //  SetName(h3, "Z from prev. step");

   //  SetName(k1, "PYTHIA Z");
   //  SetName(k2, "Flat Z");
   //  SetName(k3, "Z from prev. step");

   //  SetColor(h1, kBlack);
   //  SetColor(k1, kBlack);

   //  SetColor(h2, kRed);
   //  SetColor(k2, kRed);

   //  SetColor(h3, kGreen-2);
   //  SetColor(k3, kGreen-2);

   //  h1->GetXaxis()->SetRangeUser(0,1);
   //  h1->GetYaxis()->SetRangeUser(2*pow(10, -3), 0.9);

   //  h1->GetXaxis()->SetTitle("Z");
   //  h1->GetYaxis()->SetTitle("arb. units");

   //  k1->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
   //  k1->GetYaxis()->SetTitle("arb. units");

   //  k1->GetXaxis()->SetRangeUser(5,15);
   //  k1->GetYaxis()->SetRangeUser(2*pow(10, -4), 1.1);

   //  p6->cd(1);
   //  gPad->SetLogy();

   //  h1->Draw();
   //  h2->Draw("SAME");
   //  h3->Draw("SAME");

   //  gPad->BuildLegend();

   //  p6->cd(2);
   //  gPad->SetLogy();

   //  k1->Draw();
   //  k2->Draw("SAME");
   //  k3->Draw("SAME");

   //  gPad->BuildLegend();

    // delete p6;
}

void PlotsComparingSuperIter(){
	Method();
}