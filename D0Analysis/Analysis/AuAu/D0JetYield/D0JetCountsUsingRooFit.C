#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TAxis.h"
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

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooConstVar.h"

using namespace std;

using namespace RooFit;

Double_t fline(Double_t *x, Double_t *par)
{
    if (x[0] > 1.8 && x[0] < 1.92) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}


Double_t gline(Double_t *x, Double_t *par)
{
    if (x[0] > 1.8 && x[0] < 1.92) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

vector<double> invmassplotter(const int centbin = 0, const int d0bin = 0, TFile *file = 0x0){


    double d0yield = 0;

	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);


    RooRealVar invmassLS("invmassLS", "invmassLS", 1.8, 1.92);
    // RooRealVar invmassUS("invmassUS", "invmassUS", 1.8, 1.92);
    RooRealVar invmassSig("invmassSig", "invmassSig", 1.8, 1.92);

    TString SignalUS = "SigBg_D0Pt_";
    TString LSBg = "LSBg_D0Pt_";

    TString usfilename = "./InvMass/TxtFiles/Sigbg_";
    TString lsfilename = "./InvMass/TxtFiles/LSbg_";

    RooDataSet* Sigdata = RooDataSet::read(Form("%s%i_%i.txt", usfilename.Data(), centbin, d0bin), invmassSig) ;
    RooDataSet* LSdata = RooDataSet::read(Form("%s%i_%i.txt", lsfilename.Data(), centbin, d0bin), invmassLS) ;

    TH1F *US;
    TH1F *LS;

    // US = (TH1F *)Sigdata->createHistogram(Form("UnlikeSign_%s_%s", cent_name[directorytostart].Data(), d0pt_name[d0ptbintostart].Data()), invmassSig, Binning(200, 1.7, 2.02));
    // LS = (TH1F *)LSdata->createHistogram(Form("LikeSign_%s_%s", cent_name[directorytostart].Data(), d0pt_name[d0ptbintostart].Data()), invmassLS, Binning(200, 1.7, 2.02));

    // TH1F *Signal = (TH1F *)US->Clone(Form("Signal_%s_%s", cent_name[directorytostart].Data(), d0pt_name[d0ptbintostart].Data()));
    // Signal->Add(LS, -1);

    // for (int bin = 0; bin < Signal->GetNbinsX(); bin++){
    //     if (Signal->GetBinContent(bin) < 0) {
    //         Signal->SetBinContent(bin, 0);
    //         Signal->SetBinError(bin, 0);
    //     }
    // }

    //Constructing the Polynomial Function For LS Data

    RooRealVar p1("p1","p1", 2, -10, 10);
    RooRealVar p2("p2","p2", 1, -10, 10);

    RooPolynomial lspol2("LSpol2", "LSpol2", invmassLS, RooArgList(p1, p2));
    RooRealVar lsbkg_yield("lsbkg_yield", "LS Background Yield", LSdata->sumEntries(), 1, 1000000000);


    // Fit the LS Background here

    RooArgList shape;
    RooArgList yield;

    shape.add(lspol2); yield.add(lsbkg_yield);
    RooAddPdf LSBackgroundPDF("LSBackgroundPDF", "Background Fit", shape, yield);

    RooFitResult *bgfit = LSBackgroundPDF.fitTo(*LSdata, Extended(kTRUE), PrintEvalErrors(-1), Save(), "Q");

    RooPlot *BgPlot = invmassLS.frame(Title("LS Inv Mass"));
    LSdata->plotOn(BgPlot, MarkerSize(0.9));
    LSBackgroundPDF.plotOn(BgPlot, LineColor(kRed));

    RooArgList fitparam = bgfit->floatParsFinal(); //1st param is yield, second param is p1, third is p2

    cout << p1.getVal() << endl;

    // cout << val << endl;

    // Constructing the Polynomial Function For US Data

    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    double n1startvalue, n2startvalue;

    if (bgfit->status() > 0) {n1startvalue = 2; n2startvalue = 1;}
    else{n1startvalue = p1.getVal(); n2startvalue = p2.getVal();}

    RooRealVar n1("n1","n1", n1startvalue, -10, 10);
    RooRealVar n2("n2","n2", n2startvalue, -10, 10);

    RooPolynomial pol2("pol2", "pol2", invmassSig, RooArgList(n1, n2));
    RooRealVar bkg_yield("bkg_yield", "Background Yield", 1000, 1., 10000000);

    //Constructing the Signal Function

    RooRealVar mean("mean", "mean of gaussian", 1.865, 1.86, 1.87);
    RooRealVar sigma("sigma", "width of gaussian", 0.005, 0.03);

    RooGaussian sig("sig","Signal component", invmassSig, mean, sigma) ;
    RooRealVar sig_yield("sig_yield", "Signal Yield", 10 ,10.,10000000.);

    RooArgList finalshape;
    RooArgList finalyield;

    RooAddPdf TotalPDF("TotalPDF", "Total Fit", RooArgList(sig, pol2), RooArgList(sig_yield, bkg_yield));

    RooFitResult *fit = TotalPDF.fitTo(*Sigdata, Extended(kTRUE), PrintEvalErrors(-1), Save(), "Q");

    RooPlot *Plot = invmassSig.frame(Bins(50), Title("m_{K#pi} [GeV/#it{c}^2]"));
    Sigdata->plotOn(Plot, MarkerSize(0.9));
    // LSdata->plotOn(Plot, MarkerSize(0.9));

    // if (fit->status() == 0){
        
    // }

    TotalPDF.plotOn(Plot, Components(pol2), LineColor(kBlack));
    TotalPDF.plotOn(Plot, Components(sig), LineColor(kBlue));
    TotalPDF.plotOn(Plot, LineColor(kRed));


    cout << n1startvalue << "\t" << n2startvalue << endl;
    cout << "Sig-Bkg = " << Sigdata->numEntries()-LSdata->numEntries() << endl;
    cout << "Signal = " << sig_yield.getVal() << endl;
    cout << "Chi-Square = " <<  Plot->chiSquare() << endl;
    cout << "Status = " << fit->status() << endl;

    double D0Yield = 0;
    double D0YieldErr = 0;
    double BackgroundYield = 0;
    double BackgroundYieldErr = 0;

    if (fit->status() == 0 && Plot->chiSquare() >= 0.7){
        D0Yield = sig_yield.getVal();
        D0YieldErr = sig_yield.getError();
        BackgroundYield = bkg_yield.getVal();
        BackgroundYieldErr = bkg_yield.getError();
    }

    else{
        D0Yield = Sigdata->sumEntries("invmassSig > 1.8 && invmassSig < 1.92");

        D0YieldErr = Sigdata->sumEntries("invmassSig > 1.8 && invmassSig < 1.92");

        D0Yield -= LSdata->sumEntries("invmassLS > 1.8 && invmassLS < 1.92");

        D0YieldErr += LSdata->sumEntries("invmassLS > 1.8 && invmassLS < 1.92");
        // for (int bin = Signal->FindBin(1.83); bin <= Signal->FindBin(1.87); bin++){
        //     D0Yield += (Signal->GetBinContent(bin));
        //     D0YieldErr += pow(Signal->GetBinError(bin), 2);
        // }

        D0YieldErr = TMath::Sqrt(D0YieldErr);

        BackgroundYield += LSdata->sumEntries("invmassLS > 1.8 && invmassLS < 1.92");

        BackgroundYieldErr = TMath::Sqrt(BackgroundYield);
    }

    // Signal Histograms

    string CentralityBins_Name[5] = {"Cent_0-10", "Cent_10-20", "Cent_20-40", "Cent_40-60", "Cent_60-80"};

    TString d0ptbins[15] = {"0.0 - 10.0", "0.0 - 0.5", "0.5 - 1.0", "1.0 - 1.5", "1.5 - 2.0", "2.0 -2.5", "2.5 - 3.0", "3.0 - 4.0", "4.0 - 5.0", "5.0 - 6.0", "6.0 - 8.0", "8.0 - 10.0", "1.0 - 10.0", "1.0 - 5.0", "5.0 - 10.0"};
    TString centbins[7] = {"0 - 80 \%", "0 - 10 \%", "10 - 40 \%", "40 - 80 \%", "10 - 20 \%", "20 - 40 \%", "40 - 60 \%"};

    TString c0_filename = "./InvMass/";
    c0_filename += "D0_RooFit_Distribution_AuAu";
    c0_filename += "_Cent_";
    c0_filename += centbin;
    c0_filename += "_D0Pt_";
    c0_filename += d0bin;

    // ofstream txfile (TxtFileName.Data(), ios::out | ios::binary );

    // txfile << D0Yield << "\t" << D0YieldErr << endl;
    // txfile << Plot->chiSquare() << endl;

    cout << "================ Pull Started =============" << endl;
    RooHist* hresid = Plot->residHist() ;
    RooHist* hpull = Plot->pullHist();
    RooPlot* pullframe = invmassSig.frame(Bins(50),Title("Pull Distribution")) ;
    pullframe->addPlotable(hpull,"P") ;
  
    TCanvas *c0 = new TCanvas(c0_filename.Data(), c0_filename.Data(), 5,5,800,800);

    c0->SetLeftMargin(0.15);
    c0->SetRightMargin(0.05);
    c0->cd();

    // pullframe->Draw();
    if (fit->status() != 0) pullframe->Draw();
    else Plot->Draw();

    TPaveText *pt = new TPaveText(0.63,0.61,0.83,0.71, "NDC"); // NDC sets coords
                                       // relative to pad dimensions
    pt->SetFillStyle(4000);
    pt->SetTextSize(0.025); 
    pt->SetTextAlign(12);

    if (fit->status() == 0) auto pt_text0 = pt->AddText(Form("%s", "FIT Converged"));
    else auto pt_text00 = pt->AddText(Form("%s", "Bad FIT"));
    auto pt_text1 = pt->AddText(Form("S = %4.0f #pm %.0f", D0Yield, D0YieldErr));
    if (fit->status() == 0) auto pt_text2 = pt->AddText(Form("#frac{#chi^{2}}{NDF} = %1.1f", Plot->chiSquare()));


    pt->Draw("SAME");

    // TPaveText *data_range = new TPaveText(0.63,0.45,0.83,0.60, "NDC");
    // data_range->SetFillStyle(4000); // text is black on white
    // data_range->SetTextSize(0.025); 
    // data_range->SetTextAlign(12);

    // auto data_range_text0 = data_range->AddText(Form("%s - %s %%", cent_name[directorytostart].Data(), cent_name[directorytoend].Data()));
    // auto data_range_text2 = data_range->AddText(Form("%s < p_{T,K#pi} < %s (GeV/#it{c})", d0pt_label[d0ptbintostart].c_str(), d0pt_label[d0ptbintostop].c_str()));
    // auto data_range_text3 = data_range->AddText("|#eta| < 1");

    // data_range->Draw("SAME");


    c0_filename += ".pdf";

    c0->SaveAs(c0_filename.Data());

    delete c0;

    // delete US;
    // delete LS;

    // cout << "Total D0s = " << S_whole << endl;

    // txfile.close();

    cout << "D0 Final Yield = " << D0Yield << "\t" << D0YieldErr << "\t" << BackgroundYield << "\t" << BackgroundYieldErr << endl;

    std::vector<double> s;

    s.push_back(D0Yield);
    s.push_back(D0YieldErr);
    s.push_back(BackgroundYield);
    s.push_back(BackgroundYieldErr);

    return s;
}

void D0JetCountsUsingRooFit(){
	const int nBinsCent = 6;
  	const int nBinsD0Pt = 11;

  	TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80", "10-20", "20-40", "40-60"};

  	ofstream YieldFiles[nBinsCent + 1];

	for (int centbin = 0; centbin < 1; centbin++){

        YieldFiles[centbin].open(Form("./InvMass/RooYield/%s.txt", centbinsname[centbin].Data()));

        for (int d0bin = 0; d0bin <= 0; d0bin++){
            std::vector<double> s = invmassplotter(centbin, d0bin);

            YieldFiles[centbin] << s[0] << "\t" << s[1] << "\t" << s[2] << endl; 
        }

        YieldFiles[centbin].close();
    }
}