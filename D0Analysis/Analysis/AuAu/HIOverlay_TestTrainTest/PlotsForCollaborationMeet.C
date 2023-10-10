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

    double s[3];
    s[0] = 941.23714*nevts_0_10; //(1048*1.38188e+06+838*1.22039e+06)/(1.38188e+06+1.22039e+06);                                                                                                                              
    s[1] = 391.35550*nevts_10_40;//571 + 351 + 206;                                                                                                                                                                           
    s[2] = 56.62475*nevts_40_80;

    double taa[3] = {941.23714, 391.35550, 56.62475};

	int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
	int colors[200];

	for (int i = 0; i < 6; i++){
	    for (int j = -9; j <= 10; j++){
	        colors[i*20 + (j+9)] = col[i] + j; 
	    }
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

				UnfoldedJetPt[iter-3][siter][cent] = (TH1D *)ProcessSpectraHistogram((TH1D *)f->Get(Form("UnfoldedPt_%i", cent)));
				UnfoldedJetPt[iter-3][siter][cent]->Scale(1./VzMC[iter-3][cent]->Integral());
				UnfoldedJetPt[iter-3][siter][cent]->SetNameTitle(Form("UnfoldedPt_%i_%i_%i", iter, siter, cent), Form("UnfoldedPt_%i_%i_%i", iter, siter, cent));
				UnfoldedJetZ[iter-3][siter][cent] = (TH1D *)ProcessSpectraHistogram((TH1D *)f->Get(Form("UnfoldedZ_%i", cent)));
				UnfoldedJetZ[iter-3][siter][cent]->Scale(1./VzMC[iter-3][cent]->Integral());
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

	
	cout << "This is fine" << endl;

	const int numofplots = 11;
	int iterationstoplot[numofplots] = {1,2,4,5,10,15,20,25,30,35,39};

	/*
	TCanvas *c[3];
	for (int cent = 0; cent < 3; cent++){
		c[cent] = new TCanvas (Form("AllZPlots_%i", cent), Form("AllZPlots_%i", cent), 1000, 1000);
		// c[cent]->Divide(6,3);
		for (int iter = NLOWITERATIONS; iter <= NLOWITERATIONS; iter++){
			// cout << "====================================== ITER = " << iter << endl;
			// c[cent]->cd(iter - 3 + 1);
			c[cent]->cd();
			gPad->SetLogy();
			gPad->SetLeftMargin(0.12);
			for (int siter = 0; siter <= NSUPERITERATIONS; siter++){
				if (std::find(iterationstoplot, iterationstoplot + numofplots, siter) == iterationstoplot + numofplots) continue;
				// cout << iter << "\t" << siter << "\t" << UnfoldedJetZ[iter-3][siter][cent]->GetName() << endl;
				UnfoldedJetZ[iter-3][siter][cent]->GetXaxis()->SetRangeUser(0,1);
				UnfoldedJetZ[iter-3][siter][cent]->GetYaxis()->SetRangeUser(2*pow(10,-2), 5);
				UnfoldedJetZ[iter-3][siter][cent]->GetXaxis()->SetTitle(Form("Z"));

				UnfoldedJetZ[iter-3][siter][cent]->GetYaxis()->SetTitle(Form("#frac{1}{N} #frac{dN}{dz}"));

				SetName(UnfoldedJetZ[iter-3][siter][cent], Form("SI # %i", siter));
				SetColor(UnfoldedJetZ[iter-3][siter][cent], colors[siter*3]);
				// UnfoldedJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter*3],1);
            	// else UnfoldedJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter/2],0.5);
            UnfoldedJetZ[iter-3][siter][cent]->Draw("LEP SAME");
            // if (siter == NSUPERITERATIONS)UnfoldedJetZ[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));

            }

            MCJetZ[0]->Draw("HIST SAME");
            SetName(MCJetZ[0], "PYTHIA Z");
            SetColor(MCJetZ[0], kBlack, 1);

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
		d[cent] = new TCanvas (Form("AllPtPlots_%i", cent), Form("AllPtPlots_%i", cent), 1000, 1000);
		// c[cent]->Divide(6,3);
		for (int iter = NLOWITERATIONS; iter <= NLOWITERATIONS; iter++){
			// cout << "====================================== ITER = " << iter << endl;
			// c[cent]->cd(iter - 3 + 1);
			d[cent]->cd();
			gPad->SetLogy();
			gPad->SetLeftMargin(0.15);
			for (int siter = 0; siter <= NSUPERITERATIONS; siter++){
				if (std::find(iterationstoplot, iterationstoplot + numofplots, siter) == iterationstoplot + numofplots) continue;
					UnfoldedJetPt[iter-3][siter][cent]->GetXaxis()->SetRangeUser(5,15);
					UnfoldedJetPt[iter-3][siter][cent]->GetYaxis()->SetRangeUser(2*pow(10,-6), 0.9);
					UnfoldedJetPt[iter-3][siter][cent]->GetXaxis()->SetTitle(Form("p_{T,Jet} [GeV/#it{c}]"));

					UnfoldedJetPt[iter-3][siter][cent]->GetYaxis()->SetTitle(Form("#frac{1}{N} #frac{dN}{dp_{T}}"));

					SetName(UnfoldedJetPt[iter-3][siter][cent], Form("SI # %i", siter));
					SetColor(UnfoldedJetPt[iter-3][siter][cent], colors[siter*3]);
	            UnfoldedJetPt[iter-3][siter][cent]->Draw("LEP SAME");
            }

            MCJetPt[0]->Draw("HIST SAME");
            SetName(MCJetPt[0], "PYTHIA Z");
            SetColor(MCJetPt[0], kBlack, 1);

            gPad->BuildLegend();
        }
   }
    


 //    TCanvas *d[3];
	// for (int cent = 0; cent < 3; cent++){
	// 	d[cent] = new TCanvas (Form("AllPtPlots_%i", cent), Form("AllPtPlots_%i", cent), 2000, 1000);
	// 	d[cent]->Divide(6,3);
	// 	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
	// 		d[cent]->cd(iter - 3 + 1);
	// 		gPad->SetLogy();
	// 		for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
	// 			UnfoldedJetPt[iter-3][siter][cent]->GetXaxis()->SetRangeUser(1, 30);
	// 			// UnfoldedJetPt[iter-3][siter][cent]->GetYaxis()->SetNameTitle(Form("#frac{1}{N} #frac{dN}{dp_{T,jet}}"), Form("#frac{1}{N} #frac{dN}{dp_{T,jet}}"));
	// 			UnfoldedJetPt[iter-3][siter][cent]->GetXaxis()->SetNameTitle(Form("p_{T,jet}"), Form("p_{T,jet}"));
	// 			UnfoldedJetPt[iter-3][siter][cent]->GetYaxis()->SetRangeUser(1.0*pow(10, -8), pow(10, 0));
	// 			UnfoldedJetPt[iter-3][siter][cent]->SetLineColorAlpha(colors[siter*3],1);
 //            	UnfoldedJetPt[iter-3][siter][cent]->Draw("L SAME");
 //            	if (siter == NSUPERITERATIONS)UnfoldedJetPt[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
 //            }
	// 		// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
 //            // MCJetPt[cent]->GetXaxis()->SetRangeUser(1,30);
 //            // MCJetPt[cent]->SetLineColor(kBlack);
 //            // MCJetPt[cent]->Draw("HIST SAME");
 //        }
 //    }

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
    */

	/*
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

    int colorsforchisquareplot[5] = {kBlack, kBlue, kRed, kGreen-2, kPink};

    TCanvas *e = new TCanvas("ChiSq", "ChiSq", 2000, 1000);
    e->Divide(3,2);

    for (int cent = 0; cent < 3; cent++){
    	e->cd(cent + 1);
    	gPad->SetLogy();
    	auto legend0 = new TLegend(0.6,0.7,0.85,0.9);
    	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
    		Chi2ForJetPt[iter-3][cent]->GetXaxis()->SetTitle("#SI");
    		Chi2ForJetPt[iter-3][cent]->GetYaxis()->SetTitle("#chi^{2} for p_{T,Jet}");
    		Chi2ForJetPt[iter-3][cent]->GetXaxis()->SetRangeUser(0, 40);
    		Chi2ForJetPt[iter-3][cent]->GetYaxis()->SetRangeUser(5*pow(10, -8), 5*pow(10, 0));
    		Chi2ForJetPt[iter-3][cent]->SetLineColor(colorsforchisquareplot[iter-3]);
    		Chi2ForJetPt[iter-3][cent]->Draw("HIST SAME");
    		legend0->AddEntry(Chi2ForJetPt[iter-3][cent],Form("Regularisation Parameter %i", iter),"l");
    	}

        legend0->Draw("SAME");

    	e->cd(cent + 4);
    	gPad->SetLogy();
    	auto legend1 = new TLegend(0.6,0.7,0.85,0.9);
    	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
    		Chi2ForJetZ[iter-3][cent]->GetXaxis()->SetTitle("#SI");
    		Chi2ForJetZ[iter-3][cent]->GetYaxis()->SetTitle("#chi^{2} for Z");
    		Chi2ForJetZ[iter-3][cent]->GetXaxis()->SetRangeUser(0, 40);
    		Chi2ForJetZ[iter-3][cent]->GetYaxis()->SetRangeUser(pow(10, -2), 5*pow(10,2));
    		Chi2ForJetZ[iter-3][cent]->SetLineColor(colorsforchisquareplot[iter-3]);
    		Chi2ForJetZ[iter-3][cent]->Draw("HIST SAME");
    		legend1->AddEntry(Chi2ForJetZ[iter-3][cent],Form("Regularisation Parameter %i", iter),"l");
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
			// if (std::find(iterationstoplot, iterationstoplot + numofplots, siter) == iterationstoplot + numofplots) continue;
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
	*/
   TLegend *legend2[5];

   TCanvas *d4;
	d4 = new TCanvas (Form("RCP_%i", 2), Form("RCP_%i", 2), 3000, 600);
	d4->Divide(5);
	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
		// d3->cd();
		d4->cd(iter - 3 + 1);
		// gPad->SetLogy();
		gPad->SetLeftMargin(0.12);
		legend2[iter - 3] = new TLegend(0.15,0.45,0.4,0.9);
		for (int siter = 0; siter <= NSUPERITERATIONS; siter++){
			if (std::find(iterationstoplot, iterationstoplot + numofplots, siter) == iterationstoplot + numofplots) continue;
			cout << iter << "\t" << siter << endl;
			RatioJetPt[iter-3][siter][2]->GetXaxis()->SetRangeUser(5, 20);
			RatioJetPt[iter-3][siter][2]->GetYaxis()->SetRangeUser(0.5,5);

			RatioJetPt[iter-3][siter][2]->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
			RatioJetPt[iter-3][siter][2]->GetYaxis()->SetTitle("Ratio with PYTHIA p_{T,Jet} Spectra");
			// RatioJetPt[iter-3][siter][2]->SetLineColorAlpha(colors[siter*3],1);

			SetName(RatioJetPt[iter-3][siter][2], Form("SI # %i", siter));
			SetColor(RatioJetPt[iter-3][siter][2], colors[siter*3]);

        	RatioJetPt[iter-3][siter][2]->Draw("L SAME");
        	cout << "Integral = " << RatioJetPt[iter-3][siter][2]->Integral() << endl;
        	// if (siter == NSUPERITERATIONS)RatioJetPt[iter-3][siter][2]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
        	legend2[iter - 3]->AddEntry(RatioJetPt[iter-3][siter][2],Form("SI #%i", siter),"lp");
        }
        legend2[iter - 3]->Draw("SAME");
        legend2[iter - 3]->SetTextSize(0.045);
        gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
        // gPad->BuildLegend();
   }

   cout << "This is fine" << endl;


 // //    TCanvas *d2;
	// // d2 = new TCanvas (Form("RCP_%i", 0), Form("RCP_%i", 0), 2000, 1000);
	// // for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
	// // 	RatioJetPt[4][siter][0]->GetXaxis()->SetRangeUser(5, 20);
	// // 	RatioJetPt[4][siter][0]->GetYaxis()->SetRangeUser(-0.5, 5);
	// // 	RatioJetPt[4][siter][0]->SetNameTitle(Form("SI #%i", siter), Form("SI #%i", siter));
	// // 	RatioJetPt[4][siter][0]->SetLineColorAlpha(colors[siter],1);
 // //    	RatioJetPt[4][siter][0]->Draw("L SAME");
 // //    }
 // //    gPad->BuildLegend();

 //    cout << "============================================= WEIGHTS ===================================================" << endl;
 //    TCanvas *k[3];
	// for (int cent = 0; cent < 3; cent++){
	// 	k[cent] = new TCanvas (Form("AllWeights_%i", cent), Form("AllWeights_%i", cent), 2000, 1000);
	// 	k[cent]->Divide(6,3);
	// 	for (int iter = NLOWITERATIONS; iter <= NITERATIONS; iter++){
	// 		// cout << "====================================== ITER = " << iter << endl;
	// 		k[cent]->cd(iter - 3 + 1);
	// 		// k[cent]->cd();
	// 		gPad->SetLogy();
	// 		for (int siter = NSUPERITERATIONS; siter >= 0; siter--){
	// 			cout << iter << "\t" << siter << "\t" << cent << "\t" << Weight[iter-3][siter][cent]->GetName() << endl;
	// 			Weight[iter-3][siter][cent]->GetXaxis()->SetRangeUser(0,1);
	// 			// UnfoldedJetZ[iter-3][siter][cent]->GetYaxis()->SetNameTitle(Form("#frac{1}{N} #frac{dN}{dz}"), Form("#frac{1}{N} #frac{dN}{dz}"));
	// 			Weight[iter-3][siter][cent]->GetXaxis()->SetNameTitle(Form("z"), Form("z"));
	// 			Weight[iter-3][siter][cent]->SetLineColorAlpha(colors[siter*3],1);
 //            	// else UnfoldedJetZ[iter-3][siter][cent]->SetLineColorAlpha(colors[siter/2],0.5);
 //            	Weight[iter-3][siter][cent]->Draw("L SAME");
 //            	if (siter == NSUPERITERATIONS)Weight[iter-3][siter][0]->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));

 //            }
	// 		// gPad->SetTitle(Form("Resolution Parameter #%i", iter - 3 + 1));
 //            // MCJetZ[cent]->GetXaxis()->SetRangeUser(0,1);
 //            // MCJetZ[cent]->SetLineColor(kBlack);
 //            // MCJetZ[cent]->Draw("HIST SAME");
 //            gPad->BuildLegend();
 //        }
 //    }

}

void PlotsForCollaborationMeet(){
	Method();
}