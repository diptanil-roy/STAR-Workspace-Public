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

double invmassplotter(const int directorytostart = 0, const int directorytoend = 1, const int d0ptbintostart = 0, const int d0ptbintostop = 11, double efficiencyvaluetouse = 1., TFile *file = 0x0, TH1F *hJetPt = 0x0, TH1F *hJetPt_Eff = 0x0, const int jetptbintostart = 1, const int jetptbintostop = 8){


    double d0yield = 0;

	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);


    RooRealVar invmassLS("invmassLS", "invmassLS", 1.8, 1.92);
    RooRealVar invmassUS("invmassUS", "invmassUS", 1.8, 1.92);
    RooRealVar invmassSig("invmassSig", "invmassSig", 1.8, 1.92);

    TString TxtFileDirectory = "Feb6/InvMass/TxtFiles/";

    TString SignalUS = "SigBg_D0Pt_";
    TString USBg = "USBg_D0Pt_";
    TString LSBg = "LSBg_D0Pt_";

    const Int_t PTBINS = 11;
    TString d0pt_name[PTBINS+1] = {"00", "05", "10", "15", "20", "25", "30", "40", "50", "60", "80", "100"};

    const Int_t CENTBINS = 5;
    TString cent_name[CENTBINS+1] = {"0", "10", "20", "40", "60", "80"};

    TString signalusfilename = TxtFileDirectory + SignalUS + cent_name[directorytostart] + "_" + d0pt_name[d0ptbintostart] + ".txt";
    TString usbgfilename = TxtFileDirectory + USBg + cent_name[directorytostart] + "_" + d0pt_name[d0ptbintostart] + ".txt";
    TString lsbgfilename = TxtFileDirectory + LSBg + cent_name[directorytostart] + "_" + d0pt_name[d0ptbintostart] + ".txt";

    RooDataSet* Sigdata = RooDataSet::read(signalusfilename.Data(), invmassSig) ;
    RooDataSet* USdata = RooDataSet::read(usbgfilename.Data(), invmassSig) ;
    Sigdata->append(*USdata);
    RooDataSet* LSdata = RooDataSet::read(lsbgfilename.Data(), invmassLS) ;


    TH1F *US;
    TH1F *LS;

    US = (TH1F *)Sigdata->createHistogram(Form("UnlikeSign_%s_%s", cent_name[directorytostart].Data(), d0pt_name[d0ptbintostart].Data()), invmassSig, Binning(200, 1.7, 2.02));
    LS = (TH1F *)LSdata->createHistogram(Form("LikeSign_%s_%s", cent_name[directorytostart].Data(), d0pt_name[d0ptbintostart].Data()), invmassLS, Binning(200, 1.7, 2.02));

    TH1F *Signal = (TH1F *)US->Clone(Form("Signal_%s_%s", cent_name[directorytostart].Data(), d0pt_name[d0ptbintostart].Data()));
    Signal->Add(LS, -1);

    for (int bin = 0; bin < Signal->GetNbinsX(); bin++){
        if (Signal->GetBinContent(bin) < 0) {
            Signal->SetBinContent(bin, 0);
            Signal->SetBinError(bin, 0);
        }
    }

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

    RooRealVar mean("mean", "mean of gaussian", 1.865, 1.86, 1.88);
    RooRealVar sigma("sigma", "width of gaussian", 0.005, 0.5);

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

    if (fit->status() == -999999){
        D0Yield = sig_yield.getVal();
        D0YieldErr = sig_yield.getError();
    }

    else{
        D0Yield = Sigdata->sumEntries("invmassSig > 1.83 && invmassSig < 1.89");

        D0YieldErr = Sigdata->sumEntries("invmassSig > 1.83 && invmassSig < 1.89");

        D0Yield -= LSdata->sumEntries("invmassLS > 1.83 && invmassLS < 1.89");

        D0YieldErr += LSdata->sumEntries("invmassLS > 1.83 && invmassLS < 1.89");
        // for (int bin = Signal->FindBin(1.83); bin <= Signal->FindBin(1.87); bin++){
        //     D0Yield += (Signal->GetBinContent(bin));
        //     D0YieldErr += pow(Signal->GetBinError(bin), 2);
        // }

        D0YieldErr = TMath::Sqrt(D0YieldErr);
    }

    // Signal Histograms

    string CentralityBins_Name[5] = {"Cent_0-10", "Cent_10-20", "Cent_20-40", "Cent_40-60", "Cent_60-80"};

    const Int_t NBINS = 5;
    Double_t edges[NBINS+1] = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4};

    const int nfbins = 11;
    double fbins[nfbins+1] = {3.,4.,5.,7.,9.,11.,13.,15.,20.,30.,40., 50};

    // string cent_name[6] = {"0", "10", "20", "40", "60", "80"};
    string jetpt_name[9] = {"-0", "3", "5", "10", "15", "20", "30", "50", "100"};
    // string d0pt_name[12] = {"00", "05", "10", "15", "20", "25", "30", "40", "50", "60", "80", "100"};
    string d0pt_label[12] = {"0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "4.0", "5.0", "6.0", "8.0", "10.0"};


    TString c0_filename = "Feb6/InvMass/";
    c0_filename += CentralityBins_Name[directorytostart] + "/";
    c0_filename += "D0InvMass";
    c0_filename += "_Cent_";
    c0_filename += cent_name[directorytostart] + "-" + cent_name[directorytoend];
    c0_filename += "_JetPt_";
    c0_filename += jetpt_name[jetptbintostart] + "-" + jetpt_name[jetptbintostop];
    c0_filename += "_D0Pt_";
    c0_filename += d0pt_name[d0ptbintostart] + "-" + d0pt_name[d0ptbintostop];

    TString TxtFileName = c0_filename;
    TxtFileName += ".txt";
    // Save out the final parameters, chi-square, signal and associated errors

    //1st Row Signal: Yield YieldErr
    //2nd Row Signal: ReducedChiSquare

    ofstream txfile (TxtFileName.Data(), ios::out | ios::binary );

    txfile << D0Yield << "\t" << D0YieldErr << endl;
    txfile << Plot->chiSquare() << endl;

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

    pt->Draw("SAME");

    TPaveText *data_range = new TPaveText(0.63,0.45,0.83,0.60, "NDC");
    data_range->SetFillStyle(4000); // text is black on white
    data_range->SetTextSize(0.025); 
    data_range->SetTextAlign(12);

    auto data_range_text0 = data_range->AddText(Form("%s - %s %%", cent_name[directorytostart].Data(), cent_name[directorytoend].Data()));
    auto data_range_text2 = data_range->AddText(Form("%s < p_{T,K#pi} < %s (GeV/#it{c})", d0pt_label[d0ptbintostart].c_str(), d0pt_label[d0ptbintostop].c_str()));
    auto data_range_text3 = data_range->AddText("|#eta| < 1");

    data_range->Draw("SAME");


    c0_filename += ".pdf";

    c0->SaveAs(c0_filename.Data());

    delete c0;

    // delete US;
    // delete LS;

    // cout << "Total D0s = " << S_whole << endl;

    txfile.close();

    
    return 0;

}

void D0YieldUsingRooFit(){

    // sumplotter();
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    const Int_t PTBINS = 11;
    Double_t edges[PTBINS+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 15.0};
    Double_t mids[PTBINS] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 7.0, 9.0};

    int colorpallete[5] = {1, 2, 3, 4, 6};
    int openmarkerpallete[5] = {26, 25, 24, 27, 32};
    int closedmarkerpallete[5] = {22, 21, 20, 33, 23};

    const char* cent_name[5] = {"0 - 10", "10 - 20", "20 - 40", "40 - 60", "60 - 80"};
    string centbinedges_name[6] = {"0", "10", "20", "40", "60", "80"};
    string d0pt_name[12] = {"00", "05", "10", "15", "20", "25", "30", "40", "50", "60", "80", "100"};

    TH1F *D0Yield_Without_Efficiency_Correction[5];
    TH1F *D0Yield_With_Efficiency_Correction[5];

    TH1F *D0Yield_Without_Efficiency_Correction_CentBin[5];
    TH1F *D0Yield_With_Efficiency_Correction_CentBin[5];

    TH1F *D0JetPt_Without_Efficiency[5][12];
    TH1F *D0JetPt_With_Efficiency[5][12];

    for (int centbin = 0; centbin < 5; centbin++){
        for (int ptbin = 0; ptbin <= PTBINS; ptbin++){
            D0JetPt_Without_Efficiency[centbin][ptbin - 1] = new TH1F(Form("hJetPt_%i_%i", centbin, ptbin),Form("hJetPt_%i_%i", centbin, ptbin),150, -50, 100);
            D0JetPt_With_Efficiency[centbin][ptbin - 1] = new TH1F(Form("hJetPt_Eff_%i_%i", centbin, ptbin),Form("hJetPt_Eff_%i_%i", centbin, ptbin),150, -50, 100);
        }
    }

    // efficiencycorrection->Close();

    TFile *f = new TFile("Feb6/InvMass/D0Yield.root", "RECREATE");
    cout << "Here" << endl;

    double d0yields_in_centbins[5][PTBINS] = {0};

    for (int centbin = 0; centbin < 1; centbin++){
        for (int ptbin = 1; ptbin <= 1; ptbin++){

            d0yields_in_centbins[centbin][ptbin - 1] = invmassplotter(centbin, centbin + 1, ptbin - 1, ptbin, 1, f, D0JetPt_Without_Efficiency[centbin][ptbin - 1], D0JetPt_With_Efficiency[centbin][ptbin - 1], 1, 6);
        }
    }

    // double d0yield = invmassplotter(0, 1, 2, 3, 1, f, 0x0, 0x0, 1, 6);

}