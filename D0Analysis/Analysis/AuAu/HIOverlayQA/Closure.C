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

void ResponseMaker(int iteration, int numberofiterations=4){
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    double numberofentries = 10000000;
    // double numberofentries = 100;

    TFile *f = new TFile("Plots.root");
    TH1D *NEntries[2];
    TH1D *VzMC[3][2];
    TH1D *MCZ[3][2];
    TH1D *RecoZ[3][2];
    THn  *JetPtZ[3][2]; //Wider range of Z and Jet Pt here.


    double ratioofeventsthatpassveto[2] = {1070./1000000., 11114./1000000.}; // For scaling, we need to make sure that vetoed events are also included.
    double crosssection[2] = {64.84/66.809, 1.969/66.809};

    for (int bin = 0; bin < 2; bin++){
        f->cd(Form("pthatbin_%i", bin));
        // gROOT->ProcessLine(".ls");
        NEntries[bin] = (TH1D *)gDirectory->Get("NEvents");

        for (int i = 0; i < 3; i++){
          VzMC[i][bin] = (TH1D *)gDirectory->Get(Form("VzMC_%i", i));

          double scale = crosssection[bin]*ratioofeventsthatpassveto[bin]/VzMC[i][bin]->GetEntries();
          cout << scale << endl;

          MCZ[i][bin] = (TH1D *)gDirectory->Get(Form("MCZ_%i", i));
          RecoZ[i][bin] = (TH1D *)gDirectory->Get(Form("RecoZ_%i", i));
          JetPtZ[i][bin] = (THnD *)gDirectory->Get(Form("hMCJetPtMCZRecoJetPtRecoZ_%i", i));
          
          MCZ[i][bin]->Scale(scale);
          RecoZ[i][bin]->Scale(scale);
          JetPtZ[i][bin]->Scale(scale);
        }
    }

    for (int i = 0; i < 3; i++){
        MCZ[i][0]->Add(MCZ[i][1]);
        RecoZ[i][0]->Add(RecoZ[i][1]);
        JetPtZ[i][0]->Add(JetPtZ[i][1]);
    }

    if (iteration == 1) {
        TFile *ZFile = new TFile(Form("Closure/ZFile_%i.root", numberofiterations), "RECREATE");

        TH1D *ZForPrior[3];

        for (int i = 0; i < 3; i++){
            ZForPrior[i]= new TH1D(Form("ZForPrior_%i_%i", i, iteration), Form("ZForPrior_%i_%i", i, iteration), nBinsForZHist, ZBinsForZHist);

            for (int bin = ZForPrior[i]->FindBin(0.01); bin <= ZForPrior[i]->FindBin(0.99); bin++){
                ZForPrior[i]->SetBinContent(bin, 1);
            }

            ZForPrior[i]->Write();
        }

        
        ZFile->Close();
    }

    delete gRandom;
    gRandom = new TRandom3(0);

    THnD *MCJetPtRecoJetPtRecoZ[3][nBinsForZHist];
    const int dimcount[3] = {0,2,3};
    for (int cent = 0; cent < 3; cent++){
        for (int i = 1; i <= nBinsForZHist; i++){
            THnD *tmp = (THnD *)JetPtZ[cent][0]->Clone();
            tmp->GetAxis(1)->SetRange(i, i);
            MCJetPtRecoJetPtRecoZ[cent][i-1] = (THnD *)tmp->Projection(3, dimcount);
        }
    }

    int nBinsZHist[nDimZHist] = {nBinsJetPtForZHist, nBinsForZHist, nBinsJetPtForZHist, nBinsForZHist};

    THn *JetPtZResampled[3];

    RooUnfoldResponse *responseJetPtvZ[3];
    RooUnfoldResponse *responseJetPt[3];

    TH2D *tmpMC = new TH2D("tmpMC", "tmpMC", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);
    TH2D *tmpReco = new TH2D("tmpReco", "tmpReco", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

    TH1D *tmpMC1D = new TH1D("tmpMC1D", "tmpMC1D", nBinsJetPtForZHist, JetPtBinsForZHist);
    TH1D *tmpReco1D = new TH1D("tmpReco1D", "tmpReco1D", nBinsJetPtForZHist, JetPtBinsForZHist);

    TH2D *MCJetPtvZ[3];
    TH2D *RecoJetPtvZ[3];

    for (int i = 0; i < 3; i++){
        JetPtZResampled[i] = new THnF(Form("JetPtZResampled_%i", i), Form("JetPtZResampled_%i", i), nDimZHist, nBinsZHist, NULL, NULL);
        MCJetPtvZ[i] = new TH2D(Form("MCJetPtvZ_%i", i), Form("MCJetPtvZ_%i", i), nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);
        RecoJetPtvZ[i] = new TH2D(Form("RecoJetPtvZ_%i", i), Form("RecoJetPtvZ_%i", i), nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

        JetPtZResampled[i]->SetBinEdges(0, JetPtBinsForZHist); //MCJetPtBins
        JetPtZResampled[i]->SetBinEdges(1, ZBinsForZHist); //MCZBins
        JetPtZResampled[i]->SetBinEdges(2, JetPtBinsForZHist); //RecoJetPtBins
        JetPtZResampled[i]->SetBinEdges(3, ZBinsForZHist); //RecoZBins

        responseJetPtvZ[i] = new RooUnfoldResponse(tmpMC, tmpReco, Form("ResponseJetPtvZ_%i", i), Form("ResponseJetPtvZ_%i", i));
        responseJetPt[i] = new RooUnfoldResponse(tmpMC1D, tmpReco1D, Form("ResponseJetPt_%i", i), Form("ResponseJetPt_%i", i));

    }

    TFile *FileWithZPrior;
    if (iteration != 0)FileWithZPrior = new TFile(Form("Closure/ZFile_%i.root", numberofiterations));
    TH1D *PriorZ[3];

    for (int i = 0; i < 3; i++){
        if (iteration != 0) PriorZ[i] = (TH1D *)FileWithZPrior->Get(Form("ZForPrior_%i_%i", i, iteration));
    }

    for (int cent = 0; cent < 3; cent++){
        for (int entries = 0; entries < numberofentries; entries++){

            if (iteration == 0){
                double v[4];
                JetPtZ[cent][0]->GetRandom(v);
                JetPtZResampled[cent]->Fill(v);
                MCJetPtvZ[cent]->Fill(v[0], v[1]);
                RecoJetPtvZ[cent]->Fill(v[2], v[3]);
                responseJetPtvZ[cent]->Fill(v[2], v[3], v[0], v[1]);
                responseJetPt[cent]->Fill(v[2], v[0]);
            }
            else{        
                double v[3];
                double randZ;
                randZ = PriorZ[cent]->GetRandom();
                int bin = PriorZ[cent]->FindBin(randZ);
                MCJetPtRecoJetPtRecoZ[cent][bin-1]->GetRandom(v);
                double vtofill[4] = {v[0], randZ, v[1], v[2]};
                JetPtZResampled[cent]->Fill(vtofill);
                MCJetPtvZ[cent]->Fill(vtofill[0], vtofill[1]);
                RecoJetPtvZ[cent]->Fill(vtofill[2], vtofill[3]);
                responseJetPtvZ[cent]->Fill(vtofill[2], vtofill[3], vtofill[0], vtofill[1]);
                responseJetPt[cent]->Fill(vtofill[2], vtofill[0]);
            }
        }
    }

    TFile *g = new TFile(Form("Closure/ResponsePlots_%i_%i.root", iteration,numberofiterations), "RECREATE");
    g->cd();
    for (int cent = 0; cent < 3; cent++){
        JetPtZResampled[cent]->Write();
        MCJetPtvZ[cent]->Write();
        RecoJetPtvZ[cent]->Write();
        g->WriteObject(responseJetPt[cent], Form("ResponseJetPt_%i", cent));
        g->WriteObject(responseJetPtvZ[cent], Form("ResponseJetPtvZ_%i", cent));
    }
    g->Close();

    if (iteration != 0)FileWithZPrior->Close();


    cout << "Resampled Histograms are filled " << endl;

    gStyle->SetOptTitle(1);


    TCanvas *d = new TCanvas("d", "d", 2000, 1800);
    d->Divide(3,3);
    TH2D *tmp1, *tmp2, *tmp3;
    
    for (int i = 0; i < 3; i++){
        d->cd(i*3+1);
        gPad->SetLogz();

        THnD *n2tmp = (THnD *)JetPtZResampled[i]->Clone(Form("n2tmp_%i_%i", iteration, i));

        tmp1 = (TH2D *)n2tmp->Projection(1,0);
        tmp1->Draw("COLZ");
        // tmp1->GetZaxis()->SetRangeUser(pow(10, -9), pow(10,-2));

        d->cd(i*3+2);
        gPad->SetLogz();
        tmp2 = (TH2D *)n2tmp->Projection(3,2);
        tmp2->Draw("COLZ");
        // tmp2->GetZaxis()->SetRangeUser(pow(10, -9), pow(10,-2));

        d->cd(i*3+3);
        gPad->SetLogz();
        tmp3 = (TH2D *)n2tmp->Projection(3,1);
        tmp3->Draw("COLZ");
        // tmp3->GetZaxis()->SetRangeUser(pow(10, -9), pow(10,-2));

    }

    d->SaveAs(Form("Closure/JetPtvZResponse_pow_%i_%i.pdf", iteration, numberofiterations));

    delete d;

    for (int i = 0; i < 3; i++){
        delete JetPtZResampled[i];
        delete MCJetPtvZ[i];
        delete RecoJetPtvZ[i];
    }

    delete tmp1;
    delete tmp2;
    delete tmp3;
    
    f->Close();
}

void ClosurePlotMaker(int iteration, int numberofiterations = 4){
	TH1::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFile *g2 = new TFile(Form("Closure/ResponsePlots_%i_%i.root", 0,numberofiterations));
    g2->cd();

    TH2D *MCJetPtvZ[3];
    TH2D *RecoJetPtvZ[3];
    for (int i = 0; i < 3; i++){
        MCJetPtvZ[i] = (TH2D *)gDirectory->Get(Form("MCJetPtvZ_%i", i));
        RecoJetPtvZ[i] = (TH2D *)gDirectory->Get(Form("RecoJetPtvZ_%i", i));
    }

    TFile *g;
    g = new TFile(Form("Closure/ResponsePlots_%i_%i.root", iteration,numberofiterations));
    g->cd();

    RooUnfoldResponse *responseJetPtvZ[3];
    RooUnfoldResponse *responseJetPt[3];

    for (int i = 0; i < 3; i++){
        g->GetObject(Form("ResponseJetPt_%i", i), responseJetPt[i]);
        g->GetObject(Form("ResponseJetPtvZ_%i", i), responseJetPtvZ[i]);
    }

    double s[3];
    s[0] = 10000000; //(1048*1.38188e+06+838*1.22039e+06)/(1.38188e+06+1.22039e+06);                                                                                                                              
    s[1] = 10000000;//571 + 351 + 206;                                                                                                                                                                           
    s[2] = 10000000;

    TH2D *UnfoldedJetPtvZ[3];

    RooUnfoldBayes unfold0 (responseJetPtvZ[0], RecoJetPtvZ[0], numberofiterations);
    RooUnfoldBayes unfold1 (responseJetPtvZ[1], RecoJetPtvZ[1], numberofiterations);
    RooUnfoldBayes unfold2 (responseJetPtvZ[2], RecoJetPtvZ[2], numberofiterations);

    UnfoldedJetPtvZ[0] = (TH2D *)unfold0.Hreco();
    UnfoldedJetPtvZ[1] = (TH2D *)unfold1.Hreco();
    UnfoldedJetPtvZ[2] = (TH2D *)unfold2.Hreco();

    TH1D *PriorZ[3];

    TFile *FileWithZPrior = new TFile(Form("Closure/ZFile_%i.root", numberofiterations), "UPDATE");
    FileWithZPrior->cd();
    for (int i = 0; i < 3; i++){
        PriorZ[i] = (TH1D *)UnfoldedJetPtvZ[i]->ProjectionY();
        PriorZ[i]->SetNameTitle(Form("ZForPrior_%i_%i", i, iteration + 1), Form("ZForPrior_%i_%i", i, iteration + 1));
        PriorZ[i]->Write();
    }
    FileWithZPrior->Close();

    TH1D *MCJetPt[3];
    TH1D *MCJetZ[3];

    TH1D *RecoJetPt[3];
    TH1D *RecoJetZ[3];

    TH1D *UnfoldedJetPt[3];
    TH1D *UnfoldedJetZ[3];

    TH1D *RatioJetPt[3];
    TH1D *RatioJetZ[3];

    for (int i = 0; i < 3; i++){
        MCJetPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)MCJetPtvZ[i]->ProjectionX());
        MCJetPt[i]->SetNameTitle(Form("MCJetPt_%i_%i", i, iteration), Form("MCJetPt_%i_%i", i, iteration));
        MCJetZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)MCJetPtvZ[i]->ProjectionY());
        MCJetZ[i]->SetNameTitle(Form("MCJetZ_%i_%i", i, iteration), Form("MCJetZ_%i_%i", i, iteration));

        RecoJetPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)RecoJetPtvZ[i]->ProjectionX());
        RecoJetPt[i]->SetNameTitle(Form("RecoJetPt_%i_%i", i, iteration), Form("RecoJetPt_%i_%i", i, iteration));
        RecoJetZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)RecoJetPtvZ[i]->ProjectionY());
        RecoJetZ[i]->SetNameTitle(Form("RecoJetZ_%i_%i", i, iteration), Form("RecoJetZ_%i_%i", i, iteration));

        UnfoldedJetPt[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)UnfoldedJetPtvZ[i]->ProjectionX());
        UnfoldedJetPt[i]->SetNameTitle(Form("UnfoldedJetPt_%i_%i", i, iteration), Form("UnfoldedJetPt_%i_%i", i, iteration));
        UnfoldedJetZ[i] = (TH1D *)ProcessSpectraHistogram((TH1D *)UnfoldedJetPtvZ[i]->ProjectionY());
        UnfoldedJetZ[i]->SetNameTitle(Form("UnfoldedJetZ_%i_%i", i, iteration), Form("UnfoldedJetZ_%i_%i", i, iteration));
    }

    TCanvas *b = new TCanvas(Form("JetPt_%i", iteration), Form("JetPt_%i", iteration), 2000, 800);
    b->Divide(3);
    for (int i = 0; i < 3; i++){
        b->cd(i+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        MCJetPt[i]->Draw("P SAME");
        MCJetPt[i]->SetLineColor(kBlue);
        MCJetPt[i]->SetMarkerColor(kBlue);
        MCJetPt[i]->SetMarkerStyle(20);

        RecoJetPt[i]->Draw("P SAME");
        RecoJetPt[i]->SetLineColor(kRed);
        RecoJetPt[i]->SetMarkerColor(kRed);
        RecoJetPt[i]->SetMarkerStyle(20);

        UnfoldedJetPt[i]->Draw("P SAME");
        UnfoldedJetPt[i]->SetLineColor(kBlack);
        UnfoldedJetPt[i]->SetMarkerColor(kBlack);
        UnfoldedJetPt[i]->SetMarkerStyle(20);

        MCJetPt[i]->GetYaxis()->SetTitle("Jet Spectra");
        MCJetPt[i]->GetXaxis()->SetTitle("p_{T,jet} [GeV/#it{c}]");
    }
    // pt->Draw("SAME");

    b->SaveAs(Form("Closure/UnfoldedPtSpectra_%i_%i.pdf", iteration, numberofiterations));

    delete b;

    TCanvas *c = new TCanvas(Form("JetZ_%i", iteration), Form("JetZ_%i", iteration), 2000, 800);
    c->Divide(3);
    for (int i = 0; i < 3; i++){
        c->cd(i+1);
        gPad->SetLogy();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.1);
        MCJetZ[i]->Draw("P SAME");
        MCJetZ[i]->SetLineColor(kBlue);
        MCJetZ[i]->SetMarkerColor(kBlue);
        MCJetZ[i]->SetMarkerStyle(20);

        RecoJetZ[i]->Draw("P SAME");
        RecoJetZ[i]->SetLineColor(kRed);
        RecoJetZ[i]->SetMarkerColor(kRed);
        RecoJetZ[i]->SetMarkerStyle(20);

        UnfoldedJetZ[i]->Draw("P SAME");
        UnfoldedJetZ[i]->SetLineColor(kBlack);
        UnfoldedJetZ[i]->SetMarkerColor(kBlack);
        UnfoldedJetZ[i]->SetMarkerStyle(kBlack);

        MCJetZ[i]->GetYaxis()->SetTitle("Z Count");
        MCJetZ[i]->GetXaxis()->SetTitle("Z");
    }
    // pt->Draw("SAME");

    c->SaveAs(Form("Closure/UnfoldedZSpectra_%i_%i.pdf", iteration, numberofiterations));

    delete c;

    for (int i = 0; i < 3; i++){
        RatioJetPt[i] = (TH1D *)UnfoldedJetPt[i]->Clone(Form("RatioJetPt_%i_%i", i, iteration));
        RatioJetPt[i]->Divide(MCJetPt[i]);

        RatioJetZ[i] = (TH1D *)UnfoldedJetZ[i]->Clone(Form("RatioJetZ_%i_%i", i, iteration));
        RatioJetZ[i]->Divide(MCJetZ[i]);
    }

    TFile *f = new TFile(Form("Closure/ClosurePlots_%i.root", numberofiterations), "UPDATE");
    f->cd();
    for (int i = 0; i < 3; i++){
        MCJetPt[i]->Write();
        MCJetZ[i]->Write();

        RecoJetPt[i]->Write();
        RecoJetZ[i]->Write();

        UnfoldedJetPt[i]->Write();
        UnfoldedJetZ[i]->Write();

        RatioJetPt[i]->Write();
        RatioJetZ[i]->Write();
    }
    f->Close();

    g->Close();

    g2->Close();

}


void Closure(const int NSUPERITERATIONS = 2, int numberofiterations = 4){

    gStyle->SetPalette(kVisibleSpectrum);

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


    // double lowwavelength = 400.;
    // double highwavelength = 700.;

    // double numcolors = 200;

    // for (int i = 0; i < numcolors; i++){
    //     double wavelength = lowwavelength + (highwavelength - lowwavelength)*double(i)/numcolors;
    //     double r, g, b;

    //     spectral_color(r,g,b, wavelength);

    //     Int_t ci = TColor::GetFreeColorIndex();
    //     TColor *color = new TColor(ci, r,g,b);
    //     colors[i] = ci;
    // }

 //    TFile *out = new TFile(Form("Closure/ClosurePlots_%i.root", numberofiterations), "RECREATE");
 //    out->Close();

	// ResponseMaker(0, numberofiterations);

 //    for (int i = 1; i <= NSUPERITERATIONS; i++){
 //        cout << "Number of Super Iteration = " << i << endl;
 //        ResponseMaker(i, numberofiterations);
 //        ClosurePlotMaker(i, numberofiterations);
 //    }

    TH1D *MCJetPt[3][NSUPERITERATIONS];
    TH1D *MCJetZ[3][NSUPERITERATIONS];

    TH1D *UnfoldedJetPt[3][NSUPERITERATIONS];
    TH1D *UnfoldedJetZ[3][NSUPERITERATIONS];

    TH1D *RatioJetPt[3][NSUPERITERATIONS];
    TH1D *RatioJetZ[3][NSUPERITERATIONS];

    TFile *f = new TFile(Form("Closure/ClosurePlots_%i.root", numberofiterations));
    f->cd();

    int NSUPERITERATIONSAVAILABLE = NSUPERITERATIONS;
    
    for (int i = 0; i < 3; i++){
        for (int iter = 1; iter <= NSUPERITERATIONS; iter++){

            if (!gDirectory->GetListOfKeys()->Contains(Form("RatioJetPt_%i_%i", i, iter))) {NSUPERITERATIONSAVAILABLE = iter-1; break;}

            MCJetPt[i][iter-1] = (TH1D *)f->Get(Form("MCJetPt_%i_%i", i, iter));
            MCJetZ[i][iter-1] = (TH1D *)f->Get(Form("MCJetZ_%i_%i", i, iter));

            UnfoldedJetPt[i][iter-1] = (TH1D *)f->Get(Form("UnfoldedJetPt_%i_%i", i, iter));
            UnfoldedJetZ[i][iter-1] = (TH1D *)f->Get(Form("UnfoldedJetZ_%i_%i", i, iter));

            RatioJetPt[i][iter-1] = (TH1D *)f->Get(Form("RatioJetPt_%i_%i", i, iter));
            RatioJetZ[i][iter-1] = (TH1D *)f->Get(Form("RatioJetZ_%i_%i", i, iter));
        }
    }


    TCanvas *c = new TCanvas(Form("PtandZ"), Form("PtandZ"), 2000, 800);
    c->Divide(3,2);
    
    for (int i = 0; i < 3; i++){
        c->cd(1 + i);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            
            if (iter%2 == 0)UnfoldedJetPt[i][iter-1]->SetLineColorAlpha(colors[iter/2],1);
            else UnfoldedJetPt[i][iter-1]->SetLineColorAlpha(colors[iter/2],0.5);
            UnfoldedJetPt[i][iter-1]->SetMarkerColor(colors[iter]);

            UnfoldedJetPt[i][iter-1]->GetYaxis()->SetRangeUser(pow(10,-1), pow(10, 7));

            // if ((iter-1)%3 != 0) continue;

            UnfoldedJetPt[i][iter-1]->Draw("HIST SAME");

            // RatioJetPt[i][iter-1]->SetLineColor(col[iter - 1]);
            // RatioJetPt[i][iter-1]->SetMarkerColor(col[iter - 1]);
        }

        MCJetPt[0][i]->SetLineColor(kBlack);
        MCJetPt[0][i]->Draw("HIST SAME");

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(UnfoldedJetPt[i][0],"Iteration 1","l");
            legend->AddEntry(UnfoldedJetPt[i][NSUPERITERATIONSAVAILABLE-1],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->AddEntry(MCJetPt[0][i],Form("Truth"),"l");
            legend->Draw("SAME");
        }
    }

    for (int i = 0; i < 3; i++){
        c->cd(4 + i);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            

            if (iter%2 == 0)UnfoldedJetZ[i][iter-1]->SetLineColorAlpha(colors[iter/2],1);
            else UnfoldedJetZ[i][iter-1]->SetLineColorAlpha(colors[iter/2],0.5);

            UnfoldedJetZ[i][iter-1]->GetYaxis()->SetRangeUser(pow(10,4), pow(10, 9));


            // UnfoldedJetZ[i][iter-1]->SetLineColor(colors[iter]);
            // UnfoldedJetZ[i][iter-1]->SetMarkerColor(colors[iter]);

            // if ((iter-1)%3 != 0) continue;
            
            UnfoldedJetZ[i][iter-1]->Draw("HIST SAME");

            // RatioJetZ[i][iter-1]->SetLineColor(col[iter - 1]);
            // RatioJetZ[i][iter-1]->SetMarkerColor(col[iter - 1]);
        }

        MCJetZ[0][i]->SetLineColor(kBlack);
        MCJetZ[0][i]->Draw("HIST SAME");

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(UnfoldedJetZ[i][0],"Iteration 1","l");
            legend->AddEntry(UnfoldedJetZ[i][NSUPERITERATIONSAVAILABLE-1],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->AddEntry(MCJetZ[0][i],Form("Truth"),"l");
            legend->Draw("SAME");
        }
    }

    c->SaveAs(Form("Closure/PlotsToUpload/UnfoldedSpectraCompared_%i.pdf", numberofiterations));


    TCanvas *d = new TCanvas(Form("Ratio"), Form("Ratio"), 2000, 800);
    d->Divide(3,2);
    
    for (int i = 0; i < 3; i++){
        d->cd(1 + i);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            
            if (iter%2 == 0)RatioJetPt[i][iter-1]->SetLineColorAlpha(colors[iter/2],1);
            else RatioJetPt[i][iter-1]->SetLineColorAlpha(colors[iter/2],0.5);
            // if ((iter-1)%3 != 0) continue;

            RatioJetPt[i][iter-1]->Draw("HIST SAME");

            // RatioJetPt[i][iter-1]->SetLineColor(col[iter - 1]);
            // RatioJetPt[i][iter-1]->SetMarkerColor(col[iter - 1]);
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(RatioJetPt[i][0],"Iteration 1","l");
            legend->AddEntry(RatioJetPt[i][NSUPERITERATIONSAVAILABLE-1],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->Draw("SAME");
        }
    }

    for (int i = 0; i < 3; i++){
        d->cd(4 + i);
        gPad->SetLogy();

        for (int iter = 1; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            
            if (iter%2 == 0)RatioJetZ[i][iter-1]->SetLineColorAlpha(colors[iter/2],1);
            else RatioJetZ[i][iter-1]->SetLineColorAlpha(colors[iter/2],0.5);

            // if ((iter-1)%3 != 0) continue;
            
            RatioJetZ[i][iter-1]->Draw("HIST SAME");

            // RatioJetZ[i][iter-1]->SetLineColor(col[iter - 1]);
            // RatioJetZ[i][iter-1]->SetMarkerColor(col[iter - 1]);
        }

        if (i == 0){
            auto legend = new TLegend(0.6,0.7,0.85,0.9);
            legend->AddEntry(RatioJetPt[i][0],"Iteration 1","l");
            legend->AddEntry(RatioJetPt[i][NSUPERITERATIONSAVAILABLE-1],Form("Iteration %i", NSUPERITERATIONSAVAILABLE),"l");
            legend->Draw("SAME");
        }
    }

    d->SaveAs(Form("Closure/PlotsToUpload/RatioSpectraCompared_%i.pdf", numberofiterations));

    // gSystem->Exec(Form("mkdir Closure/RatioJetPt_%i", numberofiterations));

    // TCanvas *RatioJetPtGif[3][NSUPERITERATIONS];
    // TLegend *l[3][NSUPERITERATIONS];
    // for (int i = 0; i < 3; i++){
    //     for (int iter = 1; iter <= NSUPERITERATIONS; iter++){
    //         RatioJetPtGif[i][iter - 1] = new TCanvas(Form("RatioJetPtGif_%i_%i", i, iter), Form("RatioJetPtGif_%i_%i", i, iter), 600, 600);
    //         gPad->SetLogy(0);

    //         TH1D *tmp = (TH1D *)RatioJetPt[i][iter - 1]->Clone();
    //         tmp->Draw("HIST SAME");
    //         tmp->SetLineColor(kBlack);
    //         tmp->GetYaxis()->SetRangeUser(0,12);
    //         l[i][iter-1] = new TLegend(0.6,0.7,0.85,0.9);
    //         l[i][iter-1]->AddEntry(tmp,Form("Iteration %i", iter),"l");
    //         l[i][iter-1]->Draw("SAME");

    //         RatioJetPtGif[i][iter - 1]->SaveAs(Form("Closure/RatioJetPt_%i/RatioJetPtGif_%i_%i.pdf", numberofiterations, i, iter));

    //         delete RatioJetPtGif[i][iter - 1];
    //     }
    // }

    // gSystem->Exec(Form("mkdir Closure/RatioJetZ_%i", numberofiterations));

    // TCanvas *RatioJetZGif[3][NSUPERITERATIONS];
    // TLegend *m[3][NSUPERITERATIONS];
    // for (int i = 0; i < 3; i++){
    //     for (int iter = 1; iter <= NSUPERITERATIONS; iter++){
    //         RatioJetZGif[i][iter - 1] = new TCanvas(Form("RatioJetZGif_%i_%i", i, iter), Form("RatioJetZGif_%i_%i", i, iter), 600, 600);
    //         gPad->SetLogy(0);
    //         TH1D *tmp = (TH1D *)RatioJetZ[i][iter - 1]->Clone();
    //         tmp->Draw("HIST SAME");
    //         tmp->SetLineColor(kBlack);
    //         tmp->GetYaxis()->SetRangeUser(0,12);
    //         m[i][iter-1] = new TLegend(0.6,0.7,0.85,0.9);
    //         m[i][iter-1]->AddEntry(tmp,Form("Iteration %i", iter),"l");
    //         m[i][iter-1]->Draw("SAME");

    //         RatioJetZGif[i][iter - 1]->SaveAs(Form("Closure/RatioJetZ_%i/RatioJetZGif_%i_%i.pdf", numberofiterations, i, iter));

    //         delete RatioJetZGif[i][iter - 1];
    //     }
    // }

    double low = 0.5;
    double high = low + NSUPERITERATIONS;

    TH1D *Chi2ForJetPt[3];
    TH1D *Chi2ForJetZ[3];

    for (int i = 0; i < 3; i++){
        Chi2ForJetPt[i] = new TH1D(Form("Chi2ForJetPt_%i", i), Form("Chi2ForJetPt_%i", i), NSUPERITERATIONS, low, high);
        Chi2ForJetZ[i] = new TH1D(Form("Chi2ForJetZ_%i", i), Form("Chi2ForJetZ_%i", i), NSUPERITERATIONS, low, high);
    }

    for (int i = 0; i < 3; i++){
        for (int iter = 1; iter <= NSUPERITERATIONSAVAILABLE; iter++){
            // Chi2ForJetPt[i]->SetBinContent(iter, Chi2Test(UnfoldedJetPt[i][iter - 1], MCJetPt[i][iter - 1], "CHI2"));
            // Chi2ForJetZ[i]->SetBinContent(iter, Chi2Test(UnfoldedJetZ[i][iter - 1], MCJetZ[i][iter - 1], "CHI2"));

            Chi2ForJetPt[i]->SetBinContent(iter, UnfoldedJetPt[i][iter - 1]->Chi2Test(MCJetPt[i][iter - 1], "WW CHI2/NDF"));
            Chi2ForJetZ[i]->SetBinContent(iter, UnfoldedJetZ[i][iter - 1]->Chi2Test(MCJetZ[i][iter - 1], "WW CHI2/NDF"));
        }
    }

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

    e->SaveAs(Form("Closure/PlotsToUpload/Chi2_%i.pdf", numberofiterations));

    f->Close();

    
}
