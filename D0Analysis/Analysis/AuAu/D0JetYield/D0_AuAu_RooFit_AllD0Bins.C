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

int GetD0PtBinFor2018Analysis(Double_t pt){
  int bin = -99;
  if (pt > 0 && pt < 0.5) bin = 0;
  if (pt >= 0.5 && pt < 1.0) bin = 1;
  if (pt >= 1.0 && pt < 1.5) bin = 2;
  if (pt >= 1.5 && pt < 2.0) bin = 3;
  if (pt >= 2.0 && pt < 2.5) bin = 4;
  if (pt >= 2.5 && pt < 3.0) bin = 5;
  if (pt >= 3.0 && pt < 4.0) bin = 6;
  if (pt >= 4.0 && pt < 5.0) bin = 7;
  if (pt >= 5.0 && pt < 6.0) bin = 8;
  if (pt >= 6.0 && pt < 8.0) bin = 9;
  if (pt >= 8.0 && pt <= 10.0) bin = 10;
  return bin;
}

int GetCentBinFor2018Analysis(Double_t Cent){ //This is pretty much only going to useful to call on the efficiency corrections from the old analysis
  int centbin = -99;

  if (Cent >= 7 && Cent <= 8) centbin = 0; //0 - 10%
  else if (Cent >= 6 && Cent <= 6) centbin = 1; // 10 - 20 %
  else if (Cent >= 4 && Cent <= 5) centbin = 2; // 20 - 40 %
  else if (Cent >= 2 && Cent <= 3) centbin = 3; // 40 - 60 %
  else if (Cent >= 0 && Cent <= 1) centbin = 4; // 60 - 80 %
  return centbin;
}

int GetD0PtBin(Double_t pt){
  int bin = -99;

  if (pt > 0 && pt < 0.5) bin = 1;
  if (pt >= 0.5 && pt < 1.0) bin = 2;
  if (pt >= 1.0 && pt < 1.5) bin = 3;
  if (pt >= 1.5 && pt < 2.0) bin = 4;
  if (pt >= 2.0 && pt < 2.5) bin = 5;
  if (pt >= 2.5 && pt < 3.0) bin = 6;
  if (pt >= 3.0 && pt < 4.0) bin = 7;
  if (pt >= 4.0 && pt < 5.0) bin = 8;
  if (pt >= 5.0 && pt <= 10.0) bin = 14; // The numbering is used to keep everything consistent across the analyses.

  return bin;
}

int GetCentBin(Double_t Cent){ 
  int centbin = -99;

  if (Cent >= 7 && Cent <= 8) centbin = 1; //0 - 10%
  else if (Cent >= 4 && Cent <= 6) centbin = 2; // 10 - 40 %
  else if (Cent >= 0 && Cent <= 3) centbin = 3; // 40 - 80 % 
  return centbin; // The numbering is used to keep everything consistent across the analyses.
}


void JetAnalysis(TString cutlevel = "Standard", TString Sign = "UnlikeSign", TFile *outfile = 0x0, double fiducialptcut = 2.0){

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
  filein = TFile::Open("../JetTrees/DataWithNewEtaCut.root");
  
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

  double ptbinsmid[] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 7.0, 9.0};

  /// Double counting correction files

  TFile *dc = TFile::Open("../psn0692/Analysis/DoubleCount/MisPID_SB_Final.root");
  TGraph *dc_corr[5];

  for (int i = 0; i < 5; i++){
    dc_corr[i] = (TGraph *)dc->Get(Form("DoubleCounting_Cen_%i_SB",i));
  }

  // Bin Definition For The THn

  const int nBinsCent = 9;
  double CentBins[nBinsCent + 1]           = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  const int nBinsD0Pt = 11;
  double D0PtBins[nBinsD0Pt + 1]           = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

  const int nBinsMass = 100;
  double D0MassBins[nBinsMass + 1];
  double lowmasslimit = 1.7;
  double highmasslimit = 2.1;

  for (int i = 0; i <= nBinsMass; i++){
    D0MassBins[i] = lowmasslimit + (highmasslimit - lowmasslimit)*float(i)/float(nBinsMass);
  }

  const Int_t nBinsdR = 4;
  Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.2, 0.4};

  const int nBinsJetPt = 16;
  double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0., 1., 2., 3.,4.,5.,7.,9.,11.,13.,15.,20.,30.,50.};

  // const int nBinsJetPt = 9;
  // double JetPtBins[nBinsJetPt+1] = {-30., -10., 0., 3., 5., 7., 10., 15., 20., 30.};

  // const int nBinsJetPt = 11;
  // double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0., 1., 2., 3., 4., 5., 8., 10., 30.};

  // const int nBinsJetPt = 13;
  // double JetPtBins[nBinsJetPt+1] = {0., 1., 2., 3.,4.,5.,7.,9.,11.,13.,15.,20.,30., 50.};

  // const int nbins_jpt =12;
  // double binning_jpt[nbins_jpt+1] = {0,1,2,3,4,5,7,9,11,13,15,20,30};

  const int nbinsZ = 17;
  double ZBins[nbinsZ+1] = {-10.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 10.0};

  const Int_t nDimJet = 6;

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

  // for (int event = 0; event < 100;  event++){
  for (int event = 0; event < nentries;  event++){
    jettree->GetEntry(event);

    // if (event%1000 == 0) cout << "Event #" << event << endl;

    int d0index = -99;

    for (int itrk = 0; itrk < JetNConst; itrk++){
      if (TrackID[itrk] == 421) {
        d0index = itrk;
        break;
      }
    }
    // cout << d0bin << endl;

    if (D0Mass < 1.7 || D0Mass > 2.10) continue;
    // if (JetCorrPt < 3.0) continue;
    if (abs(JetEta) > 0.6) continue;
    // if (JetCorrPt >= 30.0) continue;
    if (PionPt <= 0.6 || KaonPt <= 0.6) continue;
    if (TrackPt[d0index] < fiducialptcut) continue;
    // if (abs(TrackEta[d0index]) > 0.6) continue;

    int d0bin2018 = GetD0PtBinFor2018Analysis(TrackPt[d0index]);
    int centbin2018 = GetCentBinFor2018Analysis(Centrality);

    /// Double counting correction bit
    double dc_correction_val;
    dc_correction_val = (1.0 - dc_corr[centbin2018]->Eval(ptbinsmid[d0bin2018]));

    double weightofevents = Weight*dc_correction_val;

    // if (TrackPt)

    double deltaD0phi = dPhi(JetPhi, standardPhi(TrackPhi[d0index]));
    double deltaD0eta = dEta(JetEta, TrackEta[d0index]);
    double deltaD0R = dR(deltaD0phi, deltaD0eta);

    double JetPx = JetCorrPt*TMath::Cos(JetPhi);
    double JetPy = JetCorrPt*TMath::Sin(JetPhi);
    double z = (JetPx*TrackPx[d0index] + JetPy*TrackPy[d0index])/pow(JetCorrPt, 2);

    // if (z > 1.2) continue;

    double toFillJet[6] = { Centrality + 0.5, JetCorrPt, TrackPt[d0index], D0Mass, deltaD0R, z };

    hD0JetCentJetPtD0PtDeltaRZ->Fill(toFillJet, weightofevents);

    int centbin = GetCentBin(Centrality);
    int d0bin = GetD0PtBin(TrackPt[d0index]);

    TString filename;
    if (Sign.CompareTo("UnlikeSign") == 0) filename = "./SignalCount/TxtFiles/Sigbg";
    else if (Sign.CompareTo("LikeSign") == 0) filename = "./SignalCount/TxtFiles/LSbg";

    ofstream file;
    file.open(Form("%s_CentBin_%i_D0PtBin_%i.txt", filename.Data(), 0, 0), ios::out | ios::app);
    file << D0Mass << "\t" << weightofevents << endl;
    file.close();

    ofstream file2;
    file2.open(Form("%s_CentBin_%i_D0PtBin_%i.txt", filename.Data(), 0, d0bin), ios::out | ios::app);
    file2 << D0Mass << "\t" << weightofevents << endl;
    file2.close();

    ofstream file3;
    file3.open(Form("%s_CentBin_%i_D0PtBin_%i.txt", filename.Data(), centbin, 0), ios::out | ios::app);
    file3 << D0Mass << "\t" << weightofevents << endl;
    file3.close();

    ofstream file4;
    file4.open(Form("%s_CentBin_%i_D0PtBin_%i.txt", filename.Data(), centbin, d0bin), ios::out | ios::app);
    file4 << D0Mass << "\t" << weightofevents << endl;
    file4.close();
  }

  outfile->cd();
  hD0JetCentJetPtD0PtDeltaRZ->Write();

  jettree->Reset();
  filein->Close();

  delete hD0JetCentJetPtD0PtDeltaRZ;

  timer.Stop();
  cout << "Real Time Used: " << timer.RealTime()/60 << "m" << endl;
}

void Method(double fiducialptcut = 2.0){
  
  const int nBinsCent = 4;
  const int nBinsD0PtForRooFit = 8;

  ofstream sigbg[nBinsCent][nBinsD0PtForRooFit];
  ofstream lsbg[nBinsCent][nBinsD0PtForRooFit];

  int allowedD0PtBinsForRooFit[] = {0, 3, 4, 5, 6, 7, 8, 14}; // Refer to the D0 pT definitions above

  for (int i = 0; i < nBinsCent; i++){
    for (int j = 0; j < nBinsD0PtForRooFit; j++){

      // cout << i << "\t" << j << endl;
      sigbg[i][j].open(Form("./SignalCount/TxtFiles/Sigbg_CentBin_%i_D0PtBin_%i.txt", i, allowedD0PtBinsForRooFit[j]), ios::out | ios::binary);
      sigbg[i][j].close();
      lsbg[i][j].open(Form("./SignalCount/TxtFiles/LSbg_CentBin_%i_D0PtBin_%i.txt", i, allowedD0PtBinsForRooFit[j]), ios::out | ios::binary);
      lsbg[i][j].close();
    }
  }

  const char *outfile;

  outfile = "./SignalCount/D0JetsDistribution.root";

  TFile *histos = new TFile(outfile, "RECREATE");

  cout << "Starting Analysis" << endl;

  JetAnalysis("Standard", "UnlikeSign", histos, fiducialptcut);

  JetAnalysis("Standard", "LikeSign", histos, fiducialptcut);

  histos->Close();
  
}

void D0_AuAu_RooFit_AllD0Bins(double fiducialptcut = 2.0){
  Method(fiducialptcut);
}
