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
#include "TPaveText.h"
// #include "StJetTreeStruct.h"
#include <vector>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

using namespace std;
#endif

#include "BinDef.h"


void Method(int NSUPERITERATIONS=200, int numberofiterations=3){

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
    int colors[200];

    for (int i = 0; i < 6; i++){
        for (int j = -9; j <= 10; j++){
            colors[i*20 + (j+9)] = col[i] + j; 
        }
    }

    TH1D *RecoJetPt[3][NSUPERITERATIONS+1];
    TH1D *RecoJetZ[3][NSUPERITERATIONS+1];

	TH1D *UnfoldedJetPt[3][NSUPERITERATIONS+1];
    TH1D *UnfoldedJetZ[3][NSUPERITERATIONS+1];

    TH1D *RatioJetPt[3][NSUPERITERATIONS+1];
    TH1D *RatioJetZ[3][NSUPERITERATIONS+1];

    TFile *f = new TFile(Form("Data/DataPlots_%i.root", numberofiterations));
    f->cd();

    int NSUPERITERATIONSAVAILABLE = NSUPERITERATIONS;
    
    for (int i = 0; i < 3; i++){
        for (int iter = 0; iter <= NSUPERITERATIONS; iter++){

            // if (!gDirectory->GetListOfKeys()->Contains(Form("RatioJetPt_%i_%i", i, iter))) {NSUPERITERATIONSAVAILABLE = iter-1; break;}

            RecoJetPt[i][iter] = (TH1D *)f->Get(Form("RecoJetPt_%i_%i", i, iter));
            RecoJetZ[i][iter] = (TH1D *)f->Get(Form("RecoJetZ_%i_%i", i, iter));

            UnfoldedJetPt[i][iter] = (TH1D *)f->Get(Form("UnfoldedJetPt_%i_%i", i, iter));
            UnfoldedJetZ[i][iter] = (TH1D *)f->Get(Form("UnfoldedJetZ_%i_%i", i, iter));

            RatioJetPt[i][iter] = (TH1D *)f->Get(Form("RCPJetPt_%i_%i", i, iter));
            RatioJetZ[i][iter] = (TH1D *)f->Get(Form("RCPJetZ_%i_%i", i, iter));
        }
    }


    TCanvas *c = new TCanvas(Form("PtandZ"), Form("PtandZ"), 2000, 800);
    c->Divide(3,2);
    
    for (int i = 0; i < 3; i++){
        c->cd(1 + i);
        gPad->SetLogy();

        for (int iter = 0; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            
            if (iter%2 == 0)UnfoldedJetPt[i][iter]->SetLineColorAlpha(colors[iter/2],1);
            else UnfoldedJetPt[i][iter]->SetLineColorAlpha(colors[iter/2],0.5);
            UnfoldedJetPt[i][iter]->SetMarkerColor(colors[iter]);

            UnfoldedJetPt[i][iter]->GetYaxis()->SetRangeUser(pow(10,-1), pow(10, 7));

            // if ((iter-1)%3 != 0) continue;

            UnfoldedJetPt[i][iter]->Draw("HIST SAME");

            // RatioJetPt[i][iter-1]->SetLineColor(col[iter - 1]);
            // RatioJetPt[i][iter-1]->SetMarkerColor(col[iter - 1]);
        }

        RecoJetPt[i][0]->Draw("HIST SAME");
        RecoJetPt[i][0]->SetLineColor(kBlack);

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(UnfoldedJetPt[i][0],"Iteration 0","l");
            legend->AddEntry(UnfoldedJetPt[i][NSUPERITERATIONSAVAILABLE],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->AddEntry(RecoJetPt[i][0], Form("Measured"));
            legend->Draw("SAME");
        }
    }

    for (int i = 0; i < 3; i++){
        c->cd(4 + i);
        gPad->SetLogy();

        for (int iter = 0; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            

            if (iter%2 == 0)UnfoldedJetZ[i][iter]->SetLineColorAlpha(colors[iter/2],1);
            else UnfoldedJetZ[i][iter]->SetLineColorAlpha(colors[iter/2],0.5);

            UnfoldedJetZ[i][iter]->GetYaxis()->SetRangeUser(pow(10,4), pow(10, 9));
            
            UnfoldedJetZ[i][iter]->Draw("HIST SAME");

        }

        RecoJetZ[i][0]->Draw("HIST SAME");
        RecoJetZ[i][0]->SetLineColor(kBlack);

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(UnfoldedJetZ[i][0],"Iteration 0","l");
            legend->AddEntry(UnfoldedJetZ[i][NSUPERITERATIONSAVAILABLE],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->AddEntry(RecoJetZ[i][0], Form("Measured"));
            legend->Draw("SAME");
        }
    }

    c->SaveAs(Form("Data/PlotsToUpload/UnfoldedSpectraCompared_%i.pdf", numberofiterations));


    TCanvas *d = new TCanvas(Form("Ratio"), Form("Ratio"), 2000, 800);
    d->Divide(2,2);
    
    for (int i = 0; i < 2; i++){
        d->cd(1 + i);
        // gPad->SetLogy();

        for (int iter = 0; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            
            if (iter%2 == 0)RatioJetPt[i][iter]->SetLineColorAlpha(colors[iter/2],1);
            else RatioJetPt[i][iter]->SetLineColorAlpha(colors[iter/2],0.5);

            RatioJetPt[i][iter]->Draw("HIST SAME");

            RatioJetPt[i][iter]->GetYaxis()->SetRangeUser(-2, 20);
            RatioJetPt[i][iter]->GetXaxis()->SetRangeUser(0, 20);

            // cout << i << "\t" << iter << endl;

        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(RatioJetPt[i][0],"Iteration 0","l");
            legend->AddEntry(RatioJetPt[i][NSUPERITERATIONSAVAILABLE],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->Draw("SAME");
        }
    }

    
    for (int i = 0; i < 2; i++){
        d->cd(3 + i);
        // gPad->SetLogy();

        for (int iter = 0; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            
            if (iter%2 == 0)RatioJetZ[i][iter]->SetLineColorAlpha(colors[iter/2],1);
            else RatioJetZ[i][iter]->SetLineColorAlpha(colors[iter/2],0.5);
            
            RatioJetZ[i][iter]->Draw("HIST SAME");

            RatioJetZ[i][iter]->GetYaxis()->SetRangeUser(-2, 20);
            RatioJetZ[i][iter]->GetXaxis()->SetRangeUser(0,1);

            // cout << i << "\t" << iter << endl;
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(RatioJetZ[i][0],"Iteration 0","l");
            legend->AddEntry(RatioJetZ[i][NSUPERITERATIONSAVAILABLE],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->Draw("SAME");
        }
    }


    d->SaveAs(Form("Data/PlotsToUpload/RatioSpectraCompared_%i.pdf", numberofiterations));

    delete c;

    delete d;

    for (int i = 0; i < 3; i++){
        for (int iter = 0; iter <= NSUPERITERATIONS; iter++){
            delete UnfoldedJetPt[i][iter];
            delete UnfoldedJetZ[i][iter];

            delete RatioJetPt[i][iter];
            delete RatioJetZ[i][iter];
        }
    }

    f->Close();

}

void DataPlots(){
    for (int i = 18; i <= 25; i++){
        Method(200, i);
    }
}