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

vector<double> invmassplotter(const int centbin = 0, const int d0bin = 0){


    double d0yield = 0;

	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    int allowedD0PtBinsForRooFit[] = {0, 3, 4, 5, 6, 7, 8, 14}; // Refer to the D0 pT definitions above

    RooRealVar invmassSig("invmassSig", "invmassSig", 1.7, 2.1);
    RooRealVar invmassSigWeight("invmassSigWeight", "invmassSigWeight", 0, 10);

    RooDataSet* data1 = RooDataSet::read(Form("./SignalCount/TxtFiles/Sigbg_CentBin_%i_D0PtBin_%i.txt", centbin, allowedD0PtBinsForRooFit[d0bin]), RooArgList(invmassSig, invmassSigWeight), "Q") ;
    RooDataSet *Sigdata = new RooDataSet(data1->GetName(),data1->GetTitle(),data1,*data1->get(),0,invmassSigWeight.GetName()) ;

    RooRealVar invmassBg("invmassBg", "invmassBg", 1.7, 2.1);
    RooRealVar invmassBgWeight("invmassBgWeight", "invmassBgWeight", 0, 10);

    RooDataSet* data2 = RooDataSet::read(Form("./SignalCount/TxtFiles/LSbg_CentBin_%i_D0PtBin_%i.txt", centbin, allowedD0PtBinsForRooFit[d0bin]), RooArgList(invmassBg, invmassBgWeight), "Q") ;
    RooDataSet *LSdata = new RooDataSet(data2->GetName(),data2->GetTitle(),data2,*data2->get(),0,invmassBgWeight.GetName()) ;

    //Constructing the Polynomial Function For LS Data

    RooRealVar p1("p1","p1", 2, -10, 10);
    RooRealVar p2("p2","p2", 1, -10, 10);

    RooPolynomial lspol2("LSpol2", "LSpol2", invmassBg, RooArgList(p1, p2));
    RooRealVar lsbkg_yield("lsbkg_yield", "LS Background Yield", LSdata->sumEntries(), 1, 1000000000);

    // Fit the LS Background here

    RooArgList shape;
    RooArgList yield;

    shape.add(lspol2); yield.add(lsbkg_yield);
    RooAddPdf LSBackgroundPDF("LSBackgroundPDF", "Background Fit", shape, yield);

    RooFitResult *bgfit = LSBackgroundPDF.fitTo(*LSdata, Extended(kTRUE), PrintEvalErrors(-1), Save(), SumW2Error(kTRUE), "Q");

    RooPlot *BgPlot = invmassBg.frame(Title("LS Inv Mass"));
    LSdata->plotOn(BgPlot, MarkerSize(0.9));
    LSBackgroundPDF.plotOn(BgPlot, LineColor(kRed));

    RooArgList fitparam = bgfit->floatParsFinal(); //1st param is yield, second param is p1, third is p2

    cout << "Background Fit Parameters = " << p1.getVal() << "\t" << p2.getVal() << endl;

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
    RooRealVar sigma("sigma", "width of gaussian", 0.005, 0.02);

    RooGaussian sig("sig","Signal component", invmassSig, mean, sigma) ;
    RooRealVar sig_yield("sig_yield", "Signal Yield", 10 ,10.,10000000.);

    RooArgList finalshape;
    RooArgList finalyield;

    RooAddPdf TotalPDF("TotalPDF", "Total Fit", RooArgList(sig, pol2), RooArgList(sig_yield, bkg_yield));

    RooFitResult *fit = TotalPDF.fitTo(*Sigdata, Extended(kTRUE), PrintEvalErrors(-1), Save(), SumW2Error(kTRUE), "Q");

    RooPlot *Plot = invmassSig.frame(Bins(70), Title("m_{K#pi} [GeV/#it{c}^2]"));
    Plot->GetYaxis()->SetTitle("Counts (per 50 MeV/#it{c}^{2})");
    Plot->GetXaxis()->SetTitle("m_{K#pi} [GeV/#it{c}^{2}]");
    if (d0bin == 0) Plot->GetYaxis()->SetRangeUser(1, 1000);
    else if (centbin == 0 && d0bin >= 1 && d0bin <= 5) Plot->GetYaxis()->SetRangeUser(1, 1200);
    else if (centbin != 3 && d0bin >= 1 && d0bin <= 5) Plot->GetYaxis()->SetRangeUser(1, 400);
    else  Plot->GetYaxis()->SetRangeUser(1, 120);
    Plot->GetYaxis()->SetNdivisions(505);

    Sigdata->plotOn(Plot, Name("invmass"), MarkerSize(0.9));
    TotalPDF.plotOn(Plot, Name("invmassbg"), Components(pol2), LineColor(kBlack));
    TotalPDF.plotOn(Plot, Name("signal"), Components(sig), DrawOption("E4C"), LineColor(kBlue));
    TotalPDF.plotOn(Plot, Name("fitfunc"), LineColor(kRed));

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
        D0Yield = Sigdata->sumEntries("invmassSig > 1.8 && invmassSig < 1.94");

        D0YieldErr = Sigdata->sumEntries("invmassSig > 1.8 && invmassSig < 1.94");

        D0Yield -= LSdata->sumEntries("invmassBg > 1.8 && invmassBg < 1.94");

        D0YieldErr += LSdata->sumEntries("invmassBg > 1.8 && invmassBg < 1.94");
        // for (int bin = Signal->FindBin(1.83); bin <= Signal->FindBin(1.87); bin++){
        //     D0Yield += (Signal->GetBinContent(bin));
        //     D0YieldErr += pow(Signal->GetBinError(bin), 2);
        // }

        D0YieldErr = TMath::Sqrt(D0YieldErr);

        BackgroundYield += LSdata->sumEntries("invmassBg > 1.8 && invmassBg < 1.94");

        BackgroundYieldErr = TMath::Sqrt(BackgroundYield);
    }

    // Signal Histograms

    TString CentralityBins_Name[4] = {"Cent_0-80", "Cent_0-10", "Cent_10-40", "Cent_40-80"};

    // TString d0ptbins[8] = {"0.0 - 10.0", "1.0 - 1.5", "1.5 - 2.0", "2.0 - 2.5", "2.5 - 3.0", "3.0 - 4.0", "2.5 - 3.0", "3.0 - 4.0", "4.0 - 5.0", "5.0 - 6.0", "6.0 - 8.0", "8.0 - 10.0", "1.0 - 10.0", "1.0 - 5.0", "5.0 - 10.0"};
    TString centbins_name[4] = {"0 - 80 \%", "0 - 10 \%", "10 - 40 \%", "40 - 80 \%"};
    TString d0ptbins_name[8] = {"1.0 - 10.0", "1.0 - 1.5", "1.5 - 2.0", "2.0 - 2.5", "2.5 - 3.0", "3.0 - 4.0", "4.0 - 5.0", "5.0 - 10.0"};

    TString c0_filename = "./SignalCount/Plots/";
    c0_filename += "D0_RooFit_Distribution_AuAu";
    c0_filename += "_Cent_";
    c0_filename += centbin;
    c0_filename += "_D0Pt_";
    c0_filename += d0bin;

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

    auto legend = new TLegend(0.90,0.65,0.57,0.88);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    legend->AddEntry("invmass","Data","lep");
    legend->AddEntry("fitfunc", "Fit (Gaus + pol2)", "l");
    legend->AddEntry("invmassbg","Combinatorial","l");
    legend->AddEntry("signal", "Signal", "l");
    legend->Draw("SAME");

    TPaveText *pt = new TPaveText(0.20,0.65,0.55,0.85, "NDC"); // NDC sets coords
                                       // relative to pad dimensions
    pt->SetFillStyle(4000);
    pt->SetTextSize(0.025); 
    pt->SetTextAlign(12);

    auto pt_text0 = pt->AddText(Form("M_{D^{0} + #bar{D^{0}}} = %1.3f #pm %1.3f GeV/#it{c}^{2}", mean.getVal(), sigma.getVal()));
    auto pt_text1 = pt->AddText(Form("S = %4.0f #pm %.0f", D0Yield, D0YieldErr));
    auto pt_text2 = pt->AddText(Form("S/B = %1.2f, #frac{#chi^{2}}{NDF} = %1.2f", D0Yield/BackgroundYield, Plot->chiSquare()));
    // auto pt_text3 = pt->AddText(Form("#frac{#chi^{2}}{NDF} = %1.2f", Plot->chiSquare()));


    // if (fit->status() == 0) auto pt_text0 = pt->AddText(Form("%s", "FIT Converged"));
    // else auto pt_text00 = pt->AddText(Form("%s", "Bad FIT"));
    // auto pt_text1 = pt->AddText(Form("S = %4.0f #pm %.0f", D0Yield, D0YieldErr));
    // if (fit->status() == 0) auto pt_text2 = pt->AddText(Form("#frac{#chi^{2}}{NDF} = %1.1f", Plot->chiSquare()));


    pt->Draw("SAME");

    TPaveText *data_range = new TPaveText(0.65,0.45,0.85,0.60, "NDC");
    data_range->SetFillStyle(4000); // text is black on white
    data_range->SetTextSize(0.025); 
    data_range->SetTextAlign(12);

    auto data_range_text0 = data_range->AddText(Form("%s", centbins_name[centbin].Data()));
    auto data_range_text2 = data_range->AddText(Form("p_{T,K#pi} #in %s [GeV/#it{c}]", d0ptbins_name[d0bin].Data()));
    auto data_range_text3 = data_range->AddText("|#eta| < 1");

    data_range->Draw("SAME");


    c0_filename += ".pdf";

    c0->SaveAs(c0_filename.Data());

    delete c0;

    cout << "D0 Final Yield = " << D0Yield << "\t" << D0YieldErr << "\t" << BackgroundYield << "\t" << BackgroundYieldErr << endl;

    std::vector<double> s;

    s.push_back(D0Yield);
    s.push_back(D0YieldErr);
    s.push_back(BackgroundYield);
    s.push_back(BackgroundYieldErr);

    return s;
}

void test(){
    RooRealVar invmassSig("invmassSig", "invmassSig", 1.8, 1.94);
    RooRealVar invmassSigWeight("invmassSigWeight", "invmassSigWeight", 0, 10);

    RooDataSet* Sigdata = RooDataSet::read("./SignalCount/TxtFiles/Sigbg_CentBin_0_D0PtBin_0.txt", RooArgList(invmassSig, invmassSigWeight), "Q") ;
    RooDataSet *wdata = new RooDataSet(Sigdata->GetName(),Sigdata->GetTitle(),Sigdata,*Sigdata->get(),0,invmassSigWeight.GetName()) ;

    Sigdata->Print();
    wdata->Print();

    if (Sigdata->isWeighted()) cout << "Something's wrong" << endl;
    if (wdata->isWeighted()) cout << "Things are OK" << endl;
}

void D0JetSignalCount(){
	const int nBinsCent = 4;
  	const int nBinsD0Pt = 8;

  	TString centbinsname[nBinsCent + 1] = {"0-80", "0-10", "10-40", "40-80"};

  	ofstream YieldFiles[nBinsCent + 1];

	for (int centbin = 0; centbin < 4; centbin++){

        YieldFiles[centbin].open(Form("./SignalCount/RooYield/%s.txt", centbinsname[centbin].Data()));

        for (int d0bin = 0; d0bin < nBinsD0Pt; d0bin++){
            std::vector<double> s = invmassplotter(centbin, d0bin);

            YieldFiles[centbin] << s[0] << "\t" << s[1] << "\t" << s[2] << "\t" << s[3] << endl; 
        }

        YieldFiles[centbin].close();
    }

    // test();
}