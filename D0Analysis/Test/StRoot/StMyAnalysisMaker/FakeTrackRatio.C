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
#include "TRandom3.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
// #include "StJetTreeStruct.h"
#include <vector>
#include <algorithm>
#include <cstdarg>

#pragma link C++ class vector<int> +;

using namespace std;

// #pragma link C++ class StJetTreeStruct+;

// #pragma link C++ class vector<float> +;
// #pragma link C++ class vector<vector<float> >+;
// #pragma link C++ class vector<int> +;
// #pragma link C++ class vector<vector<int> >+;
#endif



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

Double_t standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dEta(Double_t eta1, Double_t eta2) {
  Double_t deltaEta;
  deltaEta = eta1 - eta2;

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

int GetD0PtBin(Double_t pt){
  if (pt > 0 && pt < 1.0) return 0;
  if (pt >= 1.0 && pt < 3.0) return 1;
  if (pt >= 3.0 && pt < 5.0) return 2;
  if (pt >= 5.0 && pt < 7.0) return 3;
  if (pt >= 7.0 && pt < 10.0) return 4;
  if (pt >= 10.0 && pt < 15.0) return 5;
  if (pt >= 15.0 && pt < 20.0) return 6;
  if (pt >= 20.0 && pt < 30.0) return 7;
  if (pt >= 30.0) return 8;
}

int GetJetPtBin(Double_t pt){
  int jetbin = -99;

  if (pt < 0.0) jetbin = 0;
  else if (pt >= 0.0 && pt < 5.0) jetbin = 1;
  else if (pt >= 5.0 && pt < 10.0) jetbin = 2;
  else if (pt >= 10.0 && pt < 15.0) jetbin = 3;
  else if (pt >= 15.0 && pt < 20.0) jetbin = 4;
  else if (pt >= 20.0 && pt < 30.0) jetbin = 5;
  else if (pt >= 30.0 && pt < 50.0) jetbin = 6;
  else jetbin = 7;

  return jetbin;
}

void MakeGeneralPlots(int num, TString xtitle = "", TString ytitle = "", TString Name = "Test", ...){

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); 

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  c->SetLogy();

  va_list valist;
  va_start(valist, num);

  int binlow = 1;
  int binhigh = 1000000;

  int binlow_old = 1000000;
  int binhigh_old = 1;

  double y_high = pow(10, 10);
  double y_low = pow(10, -5);

  double y_old_high = y_low;
  double y_old_low = y_high;

  for (int i = 0; i < num; i++){

    TH1D *g =  va_arg(valist, TH1D *);

    TH1D *h = (TH1D *)g->Clone();

    if (i == 0) {

        h->GetXaxis()->SetTitle(xtitle.Data());
        h->GetYaxis()->SetTitle(ytitle.Data());
        h->GetXaxis()->SetTitleOffset(1.2);
        h->GetYaxis()->SetTitleOffset(1.4);
    }

    binlow = h->FindFirstBinAbove(0);
    binhigh = h->FindLastBinAbove(0);

    binlow = ( binlow - 5 > 1) ? binlow - 5 : 1;
    binhigh = ( binhigh + 5 < h->GetNbinsX() ) ? binhigh + 5 : h->GetNbinsX();

    binlow_old = (binlow_old < binlow) ? binlow_old : binlow;
    binhigh_old = (binhigh_old > binhigh) ? binhigh_old: binhigh;

    h->GetXaxis()->SetRange(binlow_old, binhigh_old);

    y_high = h->GetMaximum()*9;
    // y_low = ( h->GetMinimum() > 0) ? h->GetMinimum(pow(10, -5)) : pow(10, -5);

    y_low = h->GetMinimum(pow(10, -5))*0.1;

    y_old_high = (y_old_high > y_high) ? y_old_high : y_high;
    y_old_low = (y_old_low < y_low) ? y_old_low : y_low;

    h->GetYaxis()->SetRangeUser(y_low, y_old_high);
    h->GetYaxis()->SetNdivisions(505);
    h->GetXaxis()->SetNdivisions(505);
    h->SetTitle(h->GetName());
    h->SetLineColor(i+1);
    h->SetMarkerColor(i+1);
    h->SetMarkerStyle(20+i);

    h->Draw("SAME lep");

      // leg->AddEntry(h, h->GetName(), "lep");
     
  }

  va_end(valist);

  TLegend *leg = gPad->BuildLegend(0.65, 0.7, 0.85, 0.85);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);

  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);

  TString Directory = "Plots_BMeson/";
  c->SaveAs((Directory+Name+".pdf").Data());

  delete c;

}

// void JetAnalysis(){}

void Writer(TString CentBin = "4080", bool smear = kFALSE, bool useZweights = kTRUE, TString Closure = ""){

  TStopwatch timer;

  timer.Start();

  gSystem->Load("/Users/diptanilroy/Desktop/Y2020/STAR/STARAnalysis/2021/EMBEDDING/RAGHAVSFRAMEWORK/ResponsePlots/RooUnfold/libRooUnfold.so");

  // Declaring Global variables

  const int maximum = 5000;

  int runid;
  int eventid;
  float refmult;
  float centrality;
  vector <unsigned int> triggers;
  vector <double> primaryvertex;
  vector <double> primaryvertexerror;

  float hardpt;
  float hardeta;
  float hardphi;

  float jetpt;
  float jetcorrectedpt;
  float jeteta;
  float jetphi;
  float jetarea;
  float jetradius;
  float jetenergy;
  float jetnef;
  float fRhoValforjet;
  float jethighesttrackpt;
  int numberofconstituents;

  float mTrackID[MaxTrack];
  float mTrackPt[MaxTrack];
  float mTrackEta[MaxTrack];
  float mTrackPhi[MaxTrack];
  float mTrackPx[MaxTrack];
  float mTrackPy[MaxTrack];
  float mTrackPz[MaxTrack];
  float mTrackCharge[MaxTrack];

  TFile *filein;

  filein = TFile::Open("FakeJets.root");

  TString inputstring = "FakeTree/FakeJetTree";

  // TFile *f = TFile::Open("test.root");

  TTree *mcjettree = (TTree*)filein->Get(inputstring.Data());

  // // jettree->SetBranchStatus("Triggers", 0);
  // jettree->SetBranchStatus("PrimaryVertex", 0);
  // jettree->SetBranchStatus("PrimaryVertexErr", 0);

  mcjettree->SetBranchAddress("RunID", &runid);
  mcjettree->SetBranchAddress("EventId", &eventid);
  mcjettree->SetBranchAddress("RefMult", &refmult);
  mcjettree->SetBranchAddress("Centrality", &centrality);
  mcjettree->SetBranchAddress("Triggers", &triggers);
  mcjettree->SetBranchAddress("PrimaryVertex", &primaryvertex);
  mcjettree->SetBranchAddress("PrimaryVertexErr", &primaryvertexerror);
  mcjettree->SetBranchAddress("HardPt", &hardpt);
  mcjettree->SetBranchAddress("HardEta", &hardeta);
  mcjettree->SetBranchAddress("HardPhi", &hardphi);
  mcjettree->SetBranchAddress("JetPt", &jetpt);
  mcjettree->SetBranchAddress("JetCorrPt", &jetcorrectedpt);
  mcjettree->SetBranchAddress("JetEta", &jeteta);
  mcjettree->SetBranchAddress("JetPhi", &jetphi);
  mcjettree->SetBranchAddress("JetArea", &jetarea);
  mcjettree->SetBranchAddress("JetRadius", &jetradius);
  mcjettree->SetBranchAddress("JetE", &jetenergy);
  mcjettree->SetBranchAddress("JetNEF", &jetnef);
  mcjettree->SetBranchAddress("JetRhoVal", &fRhoValforjet);
  mcjettree->SetBranchAddress("JetHighestTrackPt", &jethighesttrackpt);
  mcjettree->SetBranchAddress("JetNConst", &numberofconstituents);
  
  mcjettree->SetBranchAddress("TrackID", mTrackID);
  mcjettree->SetBranchAddress("TrackPt", mTrackPt);
  mcjettree->SetBranchAddress("TrackEta", mTrackEta);
  mcjettree->SetBranchAddress("TrackPhi", mTrackPhi);
  mcjettree->SetBranchAddress("TrackPx", mTrackPx);
  mcjettree->SetBranchAddress("TrackPy", mTrackPy);
  mcjettree->SetBranchAddress("TrackPz", mTrackPz);
  mcjettree->SetBranchAddress("TrackCharge", mTrackCharge);

  const int nfbins = 11;
  double fbins[nfbins+1] = {3.,4.,5.,7.,9.,11.,13.,15.,20.,30.,40., 50};

  TH1D *CentralD0Pt = new TH1D("CentralD0Pt", "CentralD0Pt", nfbins, fbins);
  TH1D *MidCentralD0Pt = new TH1D("MidCentralD0Pt", "MidCentralD0Pt", nfbins, fbins);
  TH1D *PeripheralD0Pt = new TH1D("PeripheralD0Pt", "PeripheralD0Pt", nfbins, fbins);

  for (int event = 0; event < mcjettree->GetEntries();  event++){

    if (JetCorrPt > 3.0) continue;

    if (JetNConst < 2.0) continue;

    if (Centrality < 10) CentralD0Pt->Fill(hardpt);
    else if (Centrality >= 10 && Centrality < 40) MidCentralD0Pt->Fill(hardpt);
    else if (Centrality >= 40 && Centrality < 80) PeripheralD0Pt->Fill(hardpt);
  }

  CentralD0Pt->Scale(1./CentralD0Pt->Integral());
  MidCentralD0Pt->Scale(1./MidCentralD0Pt->Integral());
  PeripheralD0Pt->Scale(1./PeripheralD0Pt->Integral());

  MakeGeneralPlots(3, "p_{T}^{D^0} [GeV/#it{c}]", "dN^{D^0}", "FakeTracksInJet", CentralD0Pt, MidCentralD0Pt, PeripheralD0Pt);

}