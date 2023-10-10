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


#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"

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

void Method1(TH1D *ChiJetPt[3], TH1D *ChiJetZ[3], double chipt[25][3], double chiz[25][3], int SuperIteration=1, bool doPlots = kTRUE){
	TFile *f[25];


	TH1D *MCJetPt[3][25];
    TH1D *MCJetZ[3][25];

    TH1D *UnfoldedJetPt[3][25];
    TH1D *UnfoldedJetZ[3][25];

    TH1D *RatioJetPt[3][25];
    TH1D *RatioJetZ[3][25];

	for (int iter = 3; iter <= 25; iter++){
		f[iter-3] = new TFile(Form("Closure/ClosurePlots_%i.root", iter));

		for (int i = 0; i < 3; i++){
			MCJetPt[i][iter-3] = (TH1D *)f[iter-3]->Get(Form("MCJetPt_%i_%i", i, SuperIteration));
	        MCJetZ[i][iter-3] = (TH1D *)f[iter-3]->Get(Form("MCJetZ_%i_%i", i, SuperIteration));

	        UnfoldedJetPt[i][iter-3] = (TH1D *)f[iter-3]->Get(Form("UnfoldedJetPt_%i_%i", i, SuperIteration));
	        UnfoldedJetZ[i][iter-3] = (TH1D *)f[iter-3]->Get(Form("UnfoldedJetZ_%i_%i", i, SuperIteration));

	        RatioJetPt[i][iter-3] = (TH1D *)f[iter-3]->Get(Form("RatioJetPt_%i_%i", i, SuperIteration));
	        RatioJetZ[i][iter-3] = (TH1D *)f[iter-3]->Get(Form("RatioJetZ_%i_%i", i, SuperIteration));
	    }
	}

	double low = 0.5;
    double high = low + 25;

    TH1D *Chi2ForJetPt[3];
    TH1D *Chi2ForJetZ[3];

    for (int i = 0; i < 3; i++){
        Chi2ForJetPt[i] = new TH1D(Form("Chi2ForJetPt_%i_%i", i, SuperIteration), Form("Chi2ForJetPt_%i_%i", i, SuperIteration), 25, low, high);
        Chi2ForJetZ[i] = new TH1D(Form("Chi2ForJetZ_%i_%i", i, SuperIteration), Form("Chi2ForJetZ_%i_%i", i, SuperIteration), 25, low, high);
    }

    for (int i = 0; i < 3; i++){
        for (int iter = 3; iter <= 25; iter++){
            // Chi2ForJetPt[i]->SetBinContent(iter, Chi2Test(UnfoldedJetPt[i][iter - 1], MCJetPt[i][iter - 1], "CHI2"));
            // Chi2ForJetZ[i]->SetBinContent(iter, Chi2Test(UnfoldedJetZ[i][iter - 1], MCJetZ[i][iter - 1], "CHI2"));

            // Chi2ForJetPt[i]->SetBinContent(iter, UnfoldedJetPt[i][iter - 3]->Chi2Test(MCJetPt[i][iter - 3], "NORM UU CHI2/NDF"));
            // Chi2ForJetZ[i]->SetBinContent(iter, UnfoldedJetZ[i][iter - 3]->Chi2Test(MCJetZ[i][iter - 3], "NORM UU CHI2/NDF"));

            double chiptvalue = MCJetPt[i][iter - 3]->Chi2Test(UnfoldedJetPt[i][iter - 3], "WW CHI2/NDF");
            double chizvalue = MCJetZ[i][iter - 3]->Chi2Test(UnfoldedJetZ[i][iter - 3], "WW CHI2/NDF");

            Chi2ForJetPt[i]->SetBinContent(iter, chiptvalue);
            Chi2ForJetZ[i]->SetBinContent(iter, chizvalue);

            chipt[iter-3][i] = chiptvalue;
            chiz[iter-3][i] = chizvalue;


        }

        ChiJetPt[i] = (TH1D *) Chi2ForJetPt[i]->Clone();
        ChiJetZ[i] = (TH1D *)Chi2ForJetZ[i]->Clone();

        ChiJetPt[i]->SetDirectory(0);
        ChiJetZ[i]->SetDirectory(0);
    }

    if (doPlots){

        TCanvas *e = new TCanvas(Form("Chi2"), Form("Chi2"), 2000, 800);
        e->Divide(2);

        e->cd(1);
        gPad->SetLogy();
        for (int i = 0; i < 3; i++){
            Chi2ForJetPt[i]->GetYaxis()->SetRangeUser(pow(10, -1), pow(10, 7));
            Chi2ForJetPt[i]->Draw("P SAME");
            Chi2ForJetPt[i]->SetLineColor(i+1);
            Chi2ForJetPt[i]->SetMarkerColor(i+1);
            Chi2ForJetPt[i]->SetMarkerStyle(20);
        }

        e->cd(2);
        gPad->SetLogy();
        for (int i = 0; i < 3; i++){
            Chi2ForJetZ[i]->GetYaxis()->SetRangeUser(pow(10, -1), pow(10, 7));
            Chi2ForJetZ[i]->Draw("P SAME");
            Chi2ForJetZ[i]->SetLineColor(i+1);
            Chi2ForJetZ[i]->SetMarkerColor(i+1);
            Chi2ForJetZ[i]->SetMarkerStyle(20);
        }

        e->SaveAs(Form("Closure/PlotsToUpload/SuperIterChi2_%i.pdf", SuperIteration));

        delete e;

    }

    for (int iter = 3; iter <= 25; iter++){
		f[iter-3]->Close();
	}
}

void Chi2ForSuperIteration(){

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
    int colors[200];

    for (int i = 0; i < 6; i++){
        for (int j = -9; j <= 10; j++){
            colors[i*20 + (j+9)] = col[i] + j; 
        }
    }



	const int NSUPERITERATIONS = 200;

    TH1D *hJetPtChi[NSUPERITERATIONS][3];
    TH1D *hJetZChi[NSUPERITERATIONS][3];

    double chipt[NSUPERITERATIONS][25][3];
    double chiz[NSUPERITERATIONS][25][3];


	for (int siter = 1; siter <= NSUPERITERATIONS; siter++){
		Method1(hJetPtChi[siter-1], hJetZChi[siter-1], chipt[siter-1], chiz[siter-1], siter);
	}

    TCanvas *a = new TCanvas(Form("Chi2Pt"), Form("Chi2Pt"), 2000, 1200);
    a->Divide(3);

    for (int i = 0; i < 3; i++){
        a->cd(i+1);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONS; iter++){
            
            if (iter%2 == 0)hJetPtChi[iter-1][i]->SetLineColorAlpha(colors[iter/2],1);
            else hJetPtChi[iter-1][i]->SetLineColorAlpha(colors[iter/2],0.5);

            hJetPtChi[iter-1][i]->GetYaxis()->SetRangeUser(pow(10, -1), pow(10, 7));
            
            hJetPtChi[iter-1][i]->Draw("HIST SAME");
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(hJetPtChi[0][i],"Iteration 1","l");
            legend->AddEntry(hJetPtChi[NSUPERITERATIONS-1][i],Form("Iteration %i", NSUPERITERATIONS),"l");
            legend->Draw("SAME");
        }
    }

    a->SaveAs(Form("Closure/PlotsToUpload/Chi2Pt_Cent.pdf"));

    delete a;

    TCanvas *b = new TCanvas(Form("Chi2Z"), Form("Chi2Z"), 2000, 1200);
    b->Divide(3);

    for (int i = 0; i < 3; i++){
        b->cd(i+1);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONS; iter++){
            
            if (iter%2 == 0)hJetZChi[iter-1][i]->SetLineColorAlpha(colors[iter/2],1);
            else hJetZChi[iter-1][i]->SetLineColorAlpha(colors[iter/2],0.5);

            hJetZChi[iter-1][i]->GetYaxis()->SetRangeUser(pow(10, 2), pow(10, 7));
            
            hJetZChi[iter-1][i]->Draw("HIST SAME");
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(hJetZChi[0][i],"Iteration 1","l");
            legend->AddEntry(hJetZChi[NSUPERITERATIONS-1][i],Form("Iteration %i", NSUPERITERATIONS),"l");
            legend->Draw("SAME");
        }
    }

    b->SaveAs(Form("Closure/PlotsToUpload/Chi2Z_Cent.pdf"));

    delete b;

    TH2D *ChiSquarePtMapOfHyperParameters[3];
    TH2D *ChiSquareZMapOfHyperParameters[3];

    double lowsuperiter  = 0.5;
    double highsuperiter = lowsuperiter + NSUPERITERATIONS;

    double lowiter = 2.5;
    double highiter = lowiter + 23;


    double niterbins[NSUPERITERATIONS+1];
    for (int binx = 0; binx <= NSUPERITERATIONS; binx++){
        niterbins[binx] = 0.5 + (binx);

        // cout << niterbins[binx-1] << endl;
    }

    const int NITERATIONS = 23;
    double iterbins[NITERATIONS+1];
    for (int biny = 0; biny <= NITERATIONS; biny++){
        iterbins[biny] = 2.5 + (biny);

        cout << iterbins[biny] << endl;
    }


    for (int i = 0; i < 3; i++){
        ChiSquarePtMapOfHyperParameters[i] = new TH2D(Form("ChiSquarePtMapOfHyperParameters_%i", i), Form("ChiSquarePtMapOfHyperParameters_%i", i), NSUPERITERATIONS, niterbins, NITERATIONS, iterbins);
        ChiSquareZMapOfHyperParameters[i] = new TH2D(Form("ChiSquareZMapOfHyperParameters_%i", i), Form("ChiSquareZMapOfHyperParameters_%i", i), NSUPERITERATIONS, niterbins, NITERATIONS, iterbins);
        
        for (int binx = 1; binx <= NSUPERITERATIONS; binx++){
            for (int biny = 1; biny <= 23; biny++){
                ChiSquarePtMapOfHyperParameters[i]->SetBinContent(binx, biny, chipt[binx-1][biny-1][i]);
                ChiSquareZMapOfHyperParameters[i]->SetBinContent(binx, biny, chiz[binx-1][biny-1][i]);
            }
        }

    }

    TCanvas *c = new TCanvas(Form("ChiSquareMapPt"), Form("ChiSquareMapPt"), 2000, 400);
    c->Divide(3,1);
    
    for (int i = 0; i < 3; i++){
        c->cd(i+1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.1);
        ChiSquarePtMapOfHyperParameters[i]->Draw("COLZ");
    }

    c->SaveAs(Form("Closure/PlotsToUpload/ChiSquareMapPt.pdf"));


    TCanvas *d = new TCanvas(Form("ChiSquareMapZ"), Form("ChiSquareMapZ"), 2000, 400);
    d->Divide(3,1);
    
    for (int i = 0; i < 3; i++){
        d->cd(i+1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.1);
        ChiSquareZMapOfHyperParameters[i]->Draw("COLZ");
    }

    d->SaveAs(Form("Closure/PlotsToUpload/ChiSquareMapZ.pdf"));

    double percentagechangept[NSUPERITERATIONS][25][3];
    double percentagechangez[NSUPERITERATIONS][25][3];

    TH1D *hPercentChangePt[NSUPERITERATIONS][3];
    TH1D *hPercentChangeZ[NSUPERITERATIONS][3];


    for (int i = 0; i < 3; i++){
        for (int binx = 1; binx <= NSUPERITERATIONS; binx++){

            hPercentChangePt[binx-1][i] = new TH1D(Form("hPercentChangePt_%i_%i", binx, i), Form("hPercentChangePt_%i_%i", binx, i), NITERATIONS, iterbins);
            hPercentChangeZ[binx-1][i] = new TH1D(Form("hPercentChangeZ_%i_%i", binx, i), Form("hPercentChangeZ_%i_%i", binx, i), NITERATIONS, iterbins);
        }
    }
    

    for (int i = 0; i < 3; i++){
        for (int binx = 1; binx <= NSUPERITERATIONS; binx++){
            for (int biny = 2; biny <= NITERATIONS; biny++){
                // cout << i << "\t" << binx << "\t" << biny << "\t" << hPercentChangePt[binx-1][i]->GetBinCenter(biny) << "\t" << hPercentChangeZ[binx-1][i]->GetBinCenter(biny) << endl;

                double val1 = abs(chipt[binx-1][biny-1][i]-chipt[binx-1][biny-2][i])/chipt[binx-1][biny-2][i];
                double val2 = abs(chiz[binx-1][biny-1][i]-chiz[binx-1][biny-2][i])/chiz[binx-1][biny-2][i];

                // cout << hPercentChangePt[binx][i]->GetBinCenter(biny) << "\t" << percentagechangept[binx][biny][i] << endl;


                hPercentChangePt[binx-1][i]->SetBinContent(biny, val1);
                hPercentChangeZ[binx-1][i]->SetBinContent(biny, val2);
            }
        } 
    } 

    cout << "The error is before this line" << endl;

    TCanvas *e = new TCanvas(Form("ChangeOfChi2Pt"), Form("ChangeOfChi2Pt"), 2000, 1200);
    e->Divide(3);

    for (int i = 0; i < 3; i++){
        e->cd(i+1);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONS; iter++){
            
            if (iter%2 == 0)hPercentChangePt[iter-1][i]->SetLineColorAlpha(colors[iter/2],1);
            else hPercentChangePt[iter-1][i]->SetLineColorAlpha(colors[iter/2],0.5);

            hPercentChangePt[iter-1][i]->GetYaxis()->SetRangeUser(pow(10, -5), pow(10, 0));
            hPercentChangePt[iter-1][i]->Draw("HIST SAME");
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(hPercentChangePt[0][i],"Iteration 1","l");
            legend->AddEntry(hPercentChangePt[NSUPERITERATIONS-1][i],Form("Iteration %i", NSUPERITERATIONS),"l");
            legend->Draw("SAME");
        }
    }

    e->SaveAs(Form("Closure/PlotsToUpload/ChangeOfChi2Pt.pdf"));

    TCanvas *e1 = new TCanvas(Form("ChangeOfChi2Z"), Form("ChangeOfChi2Z"), 2000, 1200);
    e1->Divide(3);

    for (int i = 0; i < 3; i++){
        e1->cd(i+1);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONS; iter++){
            
            if (iter%2 == 0)hPercentChangeZ[iter-1][i]->SetLineColorAlpha(colors[iter/2],1);
            else hPercentChangeZ[iter-1][i]->SetLineColorAlpha(colors[iter/2],0.5);

            hPercentChangeZ[iter-1][i]->GetYaxis()->SetRangeUser(pow(10, -5), pow(10, 0));
            hPercentChangeZ[iter-1][i]->Draw("HIST SAME");
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(hPercentChangeZ[0][i],"Iteration 1","l");
            legend->AddEntry(hPercentChangeZ[NSUPERITERATIONS-1][i],Form("Iteration %i", NSUPERITERATIONS),"l");
            legend->Draw("SAME");
        }
    }

    e1->SaveAs(Form("Closure/PlotsToUpload/ChangeOfChi2Z.pdf"));


    double lpercentagechangept[NITERATIONS][NSUPERITERATIONS][3];
    double lpercentagechangez[NITERATIONS][NSUPERITERATIONS][3];

    TH1D *lPercentChangePt[NITERATIONS][3];
    TH1D *lPercentChangeZ[NITERATIONS][3];


    for (int i = 0; i < 3; i++){
        for (int binx = 1; binx <= NITERATIONS; binx++){

            lPercentChangePt[binx-1][i] = new TH1D(Form("lPercentChangePt_%i_%i", binx, i), Form("lPercentChangePt_%i_%i", binx, i), NSUPERITERATIONS, niterbins);
            lPercentChangeZ[binx-1][i] = new TH1D(Form("lPercentChangeZ_%i_%i", binx, i), Form("lPercentChangeZ_%i_%i", binx, i), NSUPERITERATIONS, niterbins);

            // for (int biny = 1; biny <= NSUPERITERATIONS; biny++){
            //     cout << i << "\t" << binx << "\t" << biny << "\t" << lPercentChangePt[binx-1][i]->GetBinCenter(biny) << "\t" << lPercentChangeZ[binx-1][i]->GetBinCenter(biny) << endl;
            // }
        }
    }
    

    for (int i = 0; i < 3; i++){
        for (int binx = 1; binx <= NITERATIONS; binx++){
            for (int biny = 2; biny <= NSUPERITERATIONS; biny++){
                double val1 = abs(chipt[biny-1][binx-1][i]-chipt[biny-2][binx-1][i])/chipt[biny-2][binx-1][i];
                double val2 = abs(chiz[biny-1][binx-1][i]-chiz[biny-2][binx-1][i])/chiz[biny-2][binx-1][i];

                // cout << i << "\t" << lPercentChangePt[binx-1][i]->GetBinCenter(biny) << "\t" << chipt[biny-2][binx-1][i] << "\t" << chipt[biny-1][binx-1][i] << endl;


                lPercentChangePt[binx-1][i]->SetBinContent(biny, val1);
                lPercentChangeZ[binx-1][i]->SetBinContent(biny, val2);
            }
        } 
    } 

    TCanvas *m = new TCanvas(Form("ChangeOfChi2Pt_SuperIter"), Form("ChangeOfChi2Pt_SuperIter"), 2000, 1200);
    m->Divide(3);

    for (int i = 0; i < 3; i++){
        m->cd(i+1);
        gPad->SetLogy();

        for (int iter = 1; iter <= NITERATIONS; iter++){
            
            lPercentChangePt[iter-1][i]->SetLineColorAlpha((colors[iter]-1)*3,1);

            lPercentChangePt[iter-1][i]->GetYaxis()->SetRangeUser(pow(10, -5), pow(10, 0));
            lPercentChangePt[iter-1][i]->Draw("HIST SAME");
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(lPercentChangePt[0][i],"Iteration 1","l");
            legend->AddEntry(lPercentChangePt[NITERATIONS-1][i],Form("Iteration %i", NITERATIONS),"l");
            legend->Draw("SAME");
        }
    }

    m->SaveAs(Form("Closure/PlotsToUpload/ChangeOfChi2Pt_SuperIter.pdf"));

    TCanvas *m1 = new TCanvas(Form("ChangeOfChi2Z_SuperIter"), Form("ChangeOfChi2Z_SuperIter"), 2000, 1200);
    m1->Divide(3);

    for (int i = 0; i < 3; i++){
        m1->cd(i+1);
        gPad->SetLogy();

        for (int iter = 1; iter <= NITERATIONS; iter++){
            
            lPercentChangeZ[iter-1][i]->SetLineColorAlpha((colors[iter]-1)*3,1);

            lPercentChangeZ[iter-1][i]->GetYaxis()->SetRangeUser(pow(10, -5), pow(10, 0));
            lPercentChangeZ[iter-1][i]->Draw("HIST SAME");
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(lPercentChangeZ[0][i],"Iteration 1","l");
            legend->AddEntry(lPercentChangeZ[NITERATIONS-1][i],Form("Iteration %i", NITERATIONS),"l");
            legend->Draw("SAME");
        }
    }

    m1->SaveAs(Form("Closure/PlotsToUpload/ChangeOfChi2Z_SuperIter.pdf"));

    TFile *f = new TFile("Closure/PlotsToUpload/SummaryPlots.root", "RECREATE");
    f->cd();

    e->Write();
    e1->Write();
    m->Write();
    m1->Write();

    f->Close();

    

}


