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
#include "TRandom3.h"
#include "TLegend.h"
// #include "StJetTreeStruct.h"
#include <vector>

#endif


void spectral_color(double &r,double &g,double &b,double l) // RGB <0,1> <- lambda l <400,700> [nm]
{
    double t;  r=0.0; g=0.0; b=0.0;
         if ((l>=400.0)&&(l<410.0)) { t=(l-400.0)/(410.0-400.0); r=    +(0.33*t)-(0.20*t*t); }
    else if ((l>=410.0)&&(l<475.0)) { t=(l-410.0)/(475.0-410.0); r=0.14         -(0.13*t*t); }
    else if ((l>=545.0)&&(l<595.0)) { t=(l-545.0)/(595.0-545.0); r=    +(1.98*t)-(     t*t); }
    else if ((l>=595.0)&&(l<650.0)) { t=(l-595.0)/(650.0-595.0); r=0.98+(0.06*t)-(0.40*t*t); }
    else if ((l>=650.0)&&(l<700.0)) { t=(l-650.0)/(700.0-650.0); r=0.65-(0.84*t)+(0.20*t*t); }
         if ((l>=415.0)&&(l<475.0)) { t=(l-415.0)/(475.0-415.0); g=             +(0.80*t*t); }
    else if ((l>=475.0)&&(l<590.0)) { t=(l-475.0)/(590.0-475.0); g=0.8 +(0.76*t)-(0.80*t*t); }
    else if ((l>=585.0)&&(l<639.0)) { t=(l-585.0)/(639.0-585.0); g=0.84-(0.84*t)           ; }
         if ((l>=400.0)&&(l<475.0)) { t=(l-400.0)/(475.0-400.0); b=    +(2.20*t)-(1.50*t*t); }
    else if ((l>=475.0)&&(l<560.0)) { t=(l-475.0)/(560.0-475.0); b=0.7 -(     t)+(0.30*t*t); }
}


const int nDimDaug = 8;

const int nBinsMCZ = 10;
double ZBins[nBinsMCZ + 1]                   = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
const int nBinsMCD0Pt = 9;
double MCD0PtBins[nBinsMCD0Pt + 1]           = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
const int nBinsMCJetPt = 9;
double MCJetPtBins[nBinsMCJetPt + 1]         = {1, 3, 5, 7, 9, 11, 13, 15, 20, 30};
const int nBinsRecoD0Pt = 9;
double RecoD0PtBins[nBinsRecoD0Pt + 1]       = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
const int nBinsRecoJetPt = 9;
double RecoJetPtBins[nBinsRecoJetPt + 1]     = {1, 3, 5, 7, 9, 11, 13, 15, 20, 30};
const int nBinsEta = 3;
double EtaBins[nBinsEta + 1]                 = {-1.0, -0.6, 0.6, 1.0};



const int nDimZHist = 4;

const int nBinsJetPtForZHist = 13;
// double JetPtBinsForZHist[nBinsJetPtForZHist + 1] = {-20, -10, -5, 0, 1, 3, 5, 7, 9, 11, 13, 15, 20, 30};
double JetPtBinsForZHist[nBinsJetPtForZHist + 1] = {-7, -5, -2, 0, 1, 3, 5, 7, 9, 11, 13, 15, 20, 30};

const int nBinsForZHist = 13;
// double ZBinsForZHist[nBinsForZHist + 1] = {-2, -0.001, 0.5, 1.001, 3};

double ZBinsForZHist[nBinsForZHist + 1] = {-2., -1.5, -1., -0.5, -0.001, 0.2, 0.4, 0.6, 0.8, 1.001, 1.5, 2., 2.5, 3.};


TH1D *ProcessSpectraHistogram(TH1D *h, TH1D *k = NULL){
    TH1D *R = (TH1D *)h->Clone();
    for(int i = 1;i<h->GetNbinsX()+1;i++){
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        double center = fabs(R->GetBinCenter(i));

        R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.035);
        R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.035);
    }

    if (k) {R->Divide(k); cout << "Corrected by gen efficiency" << endl;}

    return R;
}


double CalculateChi2(TH1D *Unfolded, TH1D *MC){
    if (Unfolded->GetNbinsX() != MC->GetNbinsX()) {throw invalid_argument( "Histograms do not have the same binning." );}

    double chi2 = 0;

    for (int i = 1; i <= MC->GetNbinsX(); i++){
        if (MC->GetBinContent(i) != 0) chi2 += pow((Unfolded->GetBinContent(i) - MC->GetBinContent(i)), 2)/MC->GetBinContent(i);
    }

    return chi2; 
}

