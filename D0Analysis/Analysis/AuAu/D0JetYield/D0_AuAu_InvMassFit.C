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
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
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

// int GetD0PtBin(Double_t pt){
//   if (pt > 0 && pt < 1.0) return 0;
//   if (pt >= 1.0 && pt < 3.0) return 1;
//   if (pt >= 3.0 && pt < 5.0) return 2;
//   if (pt >= 5.0 && pt < 7.0) return 3;
//   if (pt >= 7.0 && pt < 10.0) return 4;
//   if (pt >= 10.0 && pt < 15.0) return 5;
//   if (pt >= 15.0 && pt < 20.0) return 6;
//   if (pt >= 20.0 && pt < 30.0) return 7;
//   if (pt >= 30.0) return 8;
// }

int GetD0PtBin(Double_t pt){
  if (pt > 0 && pt < 0.5) return 0;
  if (pt >= 0.5 && pt < 1.0) return 1;
  if (pt >= 1.0 && pt < 1.5) return 2;
  if (pt >= 1.5 && pt < 2.0) return 3;
  if (pt >= 2.0 && pt < 2.5) return 4;
  if (pt >= 2.5 && pt < 3.0) return 5;
  if (pt >= 3.0 && pt < 4.0) return 6;
  if (pt >= 4.0 && pt < 5.0) return 7;
  if (pt >= 5.0 && pt < 6.0) return 8;
  if (pt >= 6.0 && pt < 8.0) return 9;
  if (pt >= 8.0) return 10;
}

int GetJetPtBin(Double_t pt){
  int jetbin = -99;

  if (pt < 3.0) jetbin = 0;
  else if (pt >= 3.0 && pt < 5.0) jetbin = 1;
  else if (pt >= 5.0 && pt < 10.0) jetbin = 2;
  else if (pt >= 10.0 && pt < 15.0) jetbin = 3;
  else if (pt >= 15.0 && pt < 20.0) jetbin = 4;
  else if (pt >= 20.0 && pt < 30.0) jetbin = 5;
  else if (pt >= 30.0 && pt < 50.0) jetbin = 6;
  else jetbin = 7;

  return jetbin;
}

int GetCentBin(Double_t Cent){
  int centbin = -99;

  if (Cent >= 0 && Cent < 10.0) centbin = 0;
  else if (Cent >= 10.0 && Cent < 20.0) centbin = 1;
  else if (Cent >= 20.0 && Cent < 40.0) centbin = 2;
  else if (Cent >= 40.0 && Cent < 60.0) centbin = 3;
  else if (Cent >= 60.0 && Cent < 80.0) centbin = 4;
  else centbin = -99;

  return centbin;
}

// void JetAnalysis(){}

void JetAnalysis(TString cutlevel = "Standard", TString Sign = "UnlikeSign", TFile *outfile = 0x0){

  TStopwatch timer;

  timer.Start();

  // Declaring Global variables

  int nentries;

  const int maximum = 5000;

  Int_t           RunID;
  Int_t           EventId;
  Float_t         RefMult;
  Float_t         Centrality;
  Float_t         Weight;
  vector<int> *Triggers = new vector<int>;
  vector<double>  *PrimaryVertex = new vector<double>;
  vector<double>  *PrimaryVertexErr = new vector<double>;
  Float_t         JetPt;
  Float_t         JetCorrPt;
  Float_t         JetEta;
  Float_t         JetPhi;
  Float_t         JetArea;
  Float_t         JetRadius;
  Float_t         JetE;
  Float_t         JetNEF;
  Float_t         JetRhoVal;
  Float_t         JetHighestTrackPt;
  Int_t           JetNConst;
  Float_t         D0Mass;
  Float_t         D0Pt;
  Float_t         D0Eta;
  Float_t         D0Phi;
  Float_t         PionPt;
  Float_t         PionEta;
  Float_t         PionPhi;
  Float_t         PionCharge;
  Float_t         KaonPt;
  Float_t         KaonEta;
  Float_t         KaonPhi;
  Float_t         KaonCharge;
  Float_t         TrackID[maximum];   //[numberofconstituents]
  Float_t         TrackPt[maximum];   //[numberofconstituents]
  Float_t         TrackEta[maximum];   //[numberofconstituents]
  Float_t         TrackPhi[maximum];   //[numberofconstituents]
  Float_t         TrackPx[maximum];   //[numberofconstituents]
  Float_t         TrackPy[maximum];   //[numberofconstituents]
  Float_t         TrackPz[maximum];   //[numberofconstituents]
  Float_t         TrackCharge[maximum];   //[numberofconstituents]


  TFile *filein;
  filein = TFile::Open("../JetTrees/Feb28_JetTree.root");
  
  TString inputstring;

  inputstring += Form("JetTree_%s_%s/D0Jets", cutlevel.Data(), Sign.Data());

  TTree *jettree = (TTree*)filein->Get(inputstring.Data());

  jettree->SetBranchAddress("RunID", &RunID);
  jettree->SetBranchAddress("EventId", &EventId);
  jettree->SetBranchAddress("RefMult", &RefMult);
  jettree->SetBranchAddress("Centrality", &Centrality);
  jettree->SetBranchAddress("Weight", &Weight);
  jettree->SetBranchAddress("Triggers", &Triggers);
  jettree->SetBranchAddress("PrimaryVertex", &PrimaryVertex);
  jettree->SetBranchAddress("PrimaryVertexErr", &PrimaryVertexErr);
  jettree->SetBranchAddress("JetPt", &JetPt);
  jettree->SetBranchAddress("JetCorrPt", &JetCorrPt);
  jettree->SetBranchAddress("JetEta", &JetEta);
  jettree->SetBranchAddress("JetPhi", &JetPhi);
  jettree->SetBranchAddress("JetArea", &JetArea);
  jettree->SetBranchAddress("JetRadius", &JetRadius);
  jettree->SetBranchAddress("JetE", &JetE);
  jettree->SetBranchAddress("JetNEF", &JetNEF);
  jettree->SetBranchAddress("JetRhoVal", &JetRhoVal);
  jettree->SetBranchAddress("JetHighestTrackPt", &JetHighestTrackPt);
  jettree->SetBranchAddress("JetNConst", &JetNConst);
  jettree->SetBranchAddress("D0Mass", &D0Mass);
  jettree->SetBranchAddress("PionPt", &PionPt);
  jettree->SetBranchAddress("PionEta", &PionEta);
  jettree->SetBranchAddress("PionPhi", &PionPhi);
  jettree->SetBranchAddress("PionCharge", &PionCharge);
  jettree->SetBranchAddress("KaonPt", &KaonPt);
  jettree->SetBranchAddress("KaonEta", &KaonEta);
  jettree->SetBranchAddress("KaonPhi", &KaonPhi);
  jettree->SetBranchAddress("KaonCharge", &KaonCharge);
  jettree->SetBranchAddress("TrackID", TrackID);
  jettree->SetBranchAddress("TrackPt", TrackPt);
  jettree->SetBranchAddress("TrackEta", TrackEta);
  jettree->SetBranchAddress("TrackPhi", TrackPhi);
  jettree->SetBranchAddress("TrackPx", TrackPx);
  jettree->SetBranchAddress("TrackPy", TrackPy);
  jettree->SetBranchAddress("TrackPz", TrackPz);
  jettree->SetBranchAddress("TrackCharge", TrackCharge);


  ////////////////////////////////// Make Jet Shape Histograms /////////////////////////////////


  nentries = jettree->GetEntries();
  cout << nentries << endl;

  // Defining the D0 histogram

  const int nDimDaug = 6;
  const int nBinsCent = 9;
  double CentBins[nBinsCent + 1]           = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  const int nBinsD0Pt = 11;
  double D0PtBins[nBinsD0Pt + 1]           = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
  const int nBinsDaugPt = 4;
  double DaugPtBins[nBinsDaugPt + 1]       = {0.3, 0.6, 1.2, 5.0, 10.0};
  const int nBinsRapidity = 4;
  double D0RapidityBins[nBinsRapidity + 1] = {-1.0, -0.5, 0.0, 0.5, 1.0};

  const int nBinsMass = 100;
  double D0MassBins[nBinsMass + 1];
  double lowmasslimit = 1.7;
  double highmasslimit = 2.1;

  for (int i = 0; i <= nBinsMass; i++){
    D0MassBins[i] = lowmasslimit + (highmasslimit - lowmasslimit)*float(i)/float(nBinsMass);
  }
  int nBinsDaug[nDimDaug] = {nBinsCent, nBinsD0Pt, nBinsDaugPt, nBinsMass, nBinsDaugPt, nBinsRapidity}; //cent, pt*charge, daughterpt1, m, daughterpt2*charge, rapidity

  THnF *hD0CentPtEtaMDphiDaug;
  TString histogramname;
  if (Sign.CompareTo("UnlikeSign") == 0) histogramname = "hD0CentPtEtaMDphiDaug";
  else if (Sign.CompareTo("LikeSign") == 0) histogramname = "hD0CentPtEtaMDphiDaugLikeSign";

  hD0CentPtEtaMDphiDaug = new THnF(Form("%s_%s", histogramname.Data(), cutlevel.Data()), Form("%s_%s", histogramname.Data(), cutlevel.Data()), nDimDaug, nBinsDaug, NULL, NULL);
  hD0CentPtEtaMDphiDaug->SetBinEdges(0, CentBins);
  hD0CentPtEtaMDphiDaug->SetBinEdges(1, D0PtBins);
  hD0CentPtEtaMDphiDaug->SetBinEdges(2, DaugPtBins);
  hD0CentPtEtaMDphiDaug->SetBinEdges(3, D0MassBins);
  hD0CentPtEtaMDphiDaug->SetBinEdges(4, DaugPtBins);
  hD0CentPtEtaMDphiDaug->SetBinEdges(5, D0RapidityBins);

  // Defining the Jet related Histogram

  // Includes Centrality, Jet Pt, D0 Pt, D0 Mass, Delta R, Z

  const Int_t nDimJet = 6;

  const Int_t nBinsdR = 8;
  Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8};

  const int nBinsJetPt = 10;
  double JetPtBins[nBinsJetPt+1] = {3., 4.,5.,7.,9.,11.,13.,15.,20.,30.,50.};

  const int nbinsZ = 17;
  double ZBins[nbinsZ+1] = {-10.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 10.0};

  int nBinsJet[nDimJet] = {nBinsCent, nBinsJetPt, nBinsD0Pt, nBinsMass, nBinsdR, nbinsZ};

  THnF *hD0JetCentJetPtD0PtDeltaRZ;
  TString jethistogramname;

  if (Sign.CompareTo("UnlikeSign") == 0) jethistogramname = "hD0JetCentJetPtD0PtMassDeltaRZ";
  else if (Sign.CompareTo("LikeSign") == 0) jethistogramname = "hD0JetCentJetPtD0PtMassDeltaRZLikeSign";

  hD0JetCentJetPtD0PtDeltaRZ = new THnF(Form("%s_%s", jethistogramname.Data(), cutlevel.Data()), Form("%s_%s", jethistogramname.Data(), cutlevel.Data()), nDimJet, nBinsJet, NULL, NULL);
  hD0JetCentJetPtD0PtDeltaRZ->SetBinEdges(0, CentBins);
  hD0JetCentJetPtD0PtDeltaRZ->SetBinEdges(1, JetPtBins);
  hD0JetCentJetPtD0PtDeltaRZ->SetBinEdges(2, D0PtBins);
  hD0JetCentJetPtD0PtDeltaRZ->SetBinEdges(3, D0MassBins);
  hD0JetCentJetPtD0PtDeltaRZ->SetBinEdges(4, dRBins);
  hD0JetCentJetPtD0PtDeltaRZ->SetBinEdges(5, ZBins);

  double ptbinsmid[] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 7.0, 9.0};

  /// Double counting correction files

  TFile *dc = TFile::Open("../psn0692/Analysis/DoubleCount/MisPID_SB_Final.root");
  TGraph *dc_corr[5];

  for (int i = 0; i < 5; i++){
    dc_corr[i] = (TGraph *)dc->Get(Form("DoubleCounting_Cen_%i_SB",i));
  }

  double goodevents = 0;

  for (int event = 0; event < nentries;  event++){
    jettree->GetEntry(event);

    int d0index = -99;

    for (int itrk = 0; itrk < JetNConst; itrk++){
      if (TrackID[itrk] == 421) {
        d0index = itrk;
        break;
      }
    }

    if (D0Mass < 1.70 || D0Mass > 2.10) continue;
    // if (JetCorrPt < 3.0) continue;

    // if (D0Mass < 1.8 || D0Mass > 1.92) continue;
    if (JetCorrPt < 3.0) continue;
    if (PionPt <= 0.6 || KaonPt <= 0.6) continue;
    if (TrackPt[d0index] < 1.0) continue;

    double deltaD0phi = dPhi(JetPhi, standardPhi(TrackPhi[d0index]));
    double deltaD0eta = dEta(JetEta, TrackEta[d0index]);
    double deltaD0R = dR(deltaD0phi, deltaD0eta);

    double JetPx = JetCorrPt*TMath::Cos(JetPhi);
    double JetPy = JetCorrPt*TMath::Sin(JetPhi);
    double z = (JetPx*TrackPx[d0index] + JetPy*TrackPy[d0index])/pow(JetCorrPt, 2);

    /// Double counting correction bit

    int dc_correction_bin;
    double dc_correction_val;

    if (Centrality == 7 || Centrality == 8) dc_correction_bin = 0; //  0 - 10%
    if (Centrality == 6)                    dc_correction_bin = 1; // 10 - 20%
    if (Centrality == 4 || Centrality == 5) dc_correction_bin = 2; // 20 - 40%
    if (Centrality == 2 || Centrality == 3) dc_correction_bin = 3; // 40 - 60%
    if (Centrality == 0 || Centrality == 1) dc_correction_bin = 4; // 60 - 80%

    int d0bin  = GetD0PtBin(TrackPt[d0index]);

    dc_correction_val = (1.0 - dc_corr[dc_correction_bin]->Eval(ptbinsmid[d0bin]));

    double toFillDaug[6] = { Centrality + 0.5, TrackPt[d0index], PionPt, D0Mass, KaonPt, TrackEta[d0index] };
    double toFillJet[6] = { Centrality + 0.5, JetCorrPt, TrackPt[d0index], D0Mass, deltaD0R, z };

    hD0CentPtEtaMDphiDaug->Fill(toFillDaug, Weight*dc_correction_val);
    hD0JetCentJetPtD0PtDeltaRZ->Fill(toFillJet, Weight*dc_correction_val); // There might be a call to increase stat later, for now, this ad-hoc way is good enough to exclude low pt kaons and pions

    goodevents++;

  }

  cout << "Events = " << goodevents << endl;

  outfile->cd();
  hD0CentPtEtaMDphiDaug->Write();
  hD0JetCentJetPtD0PtDeltaRZ->Write();

  jettree->Reset();
  filein->Close();

  delete hD0CentPtEtaMDphiDaug;

  timer.Stop();
  cout << "Real Time Used: " << timer.RealTime()/60 << "m" << endl;
}

void Method(){
  const char *outfile;

  outfile = "./D0JetYield.root";

  TFile *histos = new TFile(outfile, "RECREATE");

  JetAnalysis("Standard", "UnlikeSign", histos);
  JetAnalysis("Standard", "LikeSign", histos);

  histos->Close();
  
}

void D0_AuAu_InvMassFit(){
  Method();
}
