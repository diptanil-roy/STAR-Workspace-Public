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


Double_t Red[3]    = { 1.00, 0.00, 0.00};
Double_t Green[3]  = { 0.00, 1.00, 0.00};
Double_t Blue[3]   = { 1.00, 0.00, 1.00};
Double_t Length[3] = { 0.00, 0.50, 1.00 };

// void spectral_color(double &r,double &g,double &b,double l) // RGB <0,1> <- lambda l <400,700> [nm]
// {
//     double t;  r=0.0; g=0.0; b=0.0;
//          if ((l>=400.0)&&(l<410.0)) { t=(l-400.0)/(410.0-400.0); r=    +(0.33*t)-(0.20*t*t); }
//     else if ((l>=410.0)&&(l<475.0)) { t=(l-410.0)/(475.0-410.0); r=0.14         -(0.13*t*t); }
//     else if ((l>=545.0)&&(l<595.0)) { t=(l-545.0)/(595.0-545.0); r=    +(1.98*t)-(     t*t); }
//     else if ((l>=595.0)&&(l<650.0)) { t=(l-595.0)/(650.0-595.0); r=0.98+(0.06*t)-(0.40*t*t); }
//     else if ((l>=650.0)&&(l<700.0)) { t=(l-650.0)/(700.0-650.0); r=0.65-(0.84*t)+(0.20*t*t); }
//          if ((l>=415.0)&&(l<475.0)) { t=(l-415.0)/(475.0-415.0); g=             +(0.80*t*t); }
//     else if ((l>=475.0)&&(l<590.0)) { t=(l-475.0)/(590.0-475.0); g=0.8 +(0.76*t)-(0.80*t*t); }
//     else if ((l>=585.0)&&(l<639.0)) { t=(l-585.0)/(639.0-585.0); g=0.84-(0.84*t)           ; }
//          if ((l>=400.0)&&(l<475.0)) { t=(l-400.0)/(475.0-400.0); b=    +(2.20*t)-(1.50*t*t); }
//     else if ((l>=475.0)&&(l<560.0)) { t=(l-475.0)/(560.0-475.0); b=0.7 -(     t)+(0.30*t*t); }
// }

int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};


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

const int nBinsJetPtForMCZHist = 9;
double JetPtBinsForMCZHist[nBinsJetPtForMCZHist + 1] = {5, 6, 7, 8, 9, 10, 11, 13, 15, 20};

const int nBinsJetPtForZHist = 10;
double JetPtBinsForZHist[nBinsJetPtForZHist + 1] = {0, 1, 3, 5, 7, 9, 11, 13, 15, 20, 30};

const int nBinsJetPtForPlotting = 6;
double JetPtBinsForPlot[nBinsJetPtForPlotting + 1] = {5, 7, 9, 11, 13, 15, 20};

const int nBinsForMCZHist = 10;
double ZBinsForMCZHist[nBinsForMCZHist + 1] = {0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

const int nBinsForZHist = 9;
double ZBinsForZHist[nBinsForZHist + 1] = {0., 0.2, 0.4, 0.6, 0.8, 1.001, 1.5, 2., 2.5, 3.};

void ProcessSpectra(TH1 *R){
  for(int i = 1;i<R->GetNbinsX()+1;i++){
    double val = R->GetBinContent(i);
    double er = R->GetBinError(i);
    double width = R->GetBinWidth(i);
    double center = fabs(R->GetBinCenter(i));

    R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.039/2);
    R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.039/2);
  }
}


TH1 *ProcessSpectraHistogram(TH1 *h, TH1D *k = NULL){
    TH1D *R = (TH1D *)h->Clone();
    for(int i = 1;i<h->GetNbinsX()+1;i++){
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        double center = fabs(R->GetBinCenter(i));

        R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.039/2);
        R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.039/2);
    }

    if (k) {R->Divide(k); cout << "Corrected by gen efficiency" << endl;}

    return R;
}

TH1 *ProcessFONLLHistogram(TH1 *h, TH1D *k = NULL){
    TH1D *R = (TH1D *)h->Clone();
    for(int i = 1;i<h->GetNbinsX()+1;i++){
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        double center = fabs(R->GetBinCenter(i));

        R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center);
        R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center);
    }

    if (k) {R->Divide(k); cout << "Corrected by gen efficiency" << endl;}

    return R;
}

TH1 *SetName(TH1 *h, TString Name){
  h->SetNameTitle(Name.Data(), Name.Data());
  return h;
}

TH2 *SetName(TH2 *h, TString Name){
  h->SetNameTitle(Name.Data(), Name.Data());
  return h;
}

TH1 *SetAxisTitles(TH1 *h, TString xaxis = "", TString yaxis = ""){
  h->GetXaxis()->SetTitle(xaxis.Data());
  h->GetYaxis()->SetTitle(yaxis.Data());
  return h;
}

TH1 *SetColor(TH1 *h, int color, int marker = 20, double markersize = 1.2){
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(markersize);
  return h;
}

TH1D *Rebin(TH1D *h, TString Name, int nbins, double *bins){
  TH1D *k = (TH1D *)h->Rebin(nbins, Name.Data(), bins);  
  return k;
}

// TH1D *GetLineAtOne(TH1D *h){
//   TH1D *one = (TH1D *)h->Clone();
//   SetName(one, Form("%s_one", h->GetName()));

//   for (int i = 1; i <= h->GetNbinsX(); i++){
//     one->SetBinContent(i, 1.);
//     one->SetBinError(i, 0);
//   }

//   one = (TH1D *)SetColor(one, kBlack, 29, 1.2);
//   return one;
// }

TLine *GetLineAtOne(TH1D *h){
  double low = h->GetBinLowEdge(1);
  double high = h->GetBinLowEdge(h->GetNbinsX()) + h->GetBinWidth(h->GetNbinsX());

  TLine *l = new TLine(low, 1, high, 1);
  l->SetLineColor(kBlack);
  return l;
}

double CalculateChi2(TH1D *Unfolded, TH1D *MC){
    if (Unfolded->GetNbinsX() != MC->GetNbinsX()) {throw invalid_argument( "Histograms do not have the same binning." );}

    double chi2 = 0;

    for (int i = 1; i <= MC->GetNbinsX(); i++){
        if (MC->GetBinContent(i) != 0) chi2 += pow((Unfolded->GetBinContent(i) - MC->GetBinContent(i)), 2)/MC->GetBinContent(i);
    }

    return chi2; 
}


const double R = 0.4;
const double deltar = 0.05;
const int numberofbins = R/deltar;


// Important functions that are used repeatedly in the calculations

// function to calculate relative phi between 2 objects and shift between 0 and 2pi 
//___________________________________________________________________________________________
Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

int GetSignOfDeltaPhi(Double_t phi1, Double_t phi2){ //it's always returning the sign of phi1 - phi2
  if (phi2 > 1.5*(TMath::Pi()) && phi2 < 2*(TMath::Pi()) && phi1 > 0 && phi1 < 0.5*(TMath::Pi())) return 1;
  else return ( (phi1 - phi2 > 0) - (phi1 - phi2 < 0) );
}

Double_t standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard > 2*(TMath::Pi())) phi_standard-=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dEta(Double_t eta1, Double_t eta2) {
  Double_t deltaEta;
  deltaEta = abs(eta1 - eta2);

  return deltaEta;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dR(Double_t delphi, Double_t deleta) {
  Double_t dRad;
  dRad = TMath::Sqrt(pow(delphi,2) + pow(deleta,2));

  return dRad;
}

// function to calculate invariant mass of a pair of objects
//___________________________________________________________________________________________
Double_t Mass(Double_t m1, Double_t m2, Double_t E1, Double_t E2, Double_t px1, Double_t px2, Double_t py1, Double_t py2, Double_t pz1, Double_t pz2) {
  Double_t m;
  m = TMath::Sqrt(pow(m1,2) + pow(m2,2) + 2*(E1*E2 - px1*px2 - py1*py2 - pz1*pz2));

  return m;
}

// function to calculate invariant mass of a pair of objects
//___________________________________________________________________________________________
Double_t pTforD0(Double_t px1, Double_t px2, Double_t py1, Double_t py2) {
  Double_t m;
  m = TMath::Sqrt(pow(px1 + px2, 2) + pow(py1 + py2, 2));

  return m;
}

Double_t pX(Double_t pT, Double_t phi){
  Double_t px;
  px = pT*TMath::Cos(phi);

  return px;
}

Double_t pY(Double_t pT, Double_t phi){
  Double_t py;
  py = pT*TMath::Sin(phi);

  return py;
}

Double_t p(Double_t E, Double_t m){
   
  Double_t pmag;
  pmag = TMath::Sqrt(pow(E,2) - pow(m,2));

  return pmag;
}

Double_t p(Double_t px, Double_t py, Double_t pz){
   
  Double_t pmag;
  pmag = TMath::Sqrt(pow(px,2) + pow(py,2) + pow(pz, 2));

  return pmag;
}

Double_t pZ(Double_t E, Double_t m, Double_t eta){
   
  Double_t pz;
  pz = p(E,m)*(TMath::Exp(2*eta)-1)/(TMath::Exp(2*eta)+1);

  return pz;
}



