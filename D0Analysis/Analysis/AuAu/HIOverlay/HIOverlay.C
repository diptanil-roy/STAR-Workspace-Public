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


void Method(int pthatbin = 0){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TFile *f;
  if (pthatbin == 0)f = new TFile("pt_0_3.root");
  else if (pthatbin == 1)f = new TFile("pt_3_inf.root");
  f->cd("HIJetSaver");

  TTree *JetTree = (TTree *)gDirectory->Get("Jets");
  // TTree *RecoJetTree = (TTree *)gDirectory->Get("RecoJets");

  Float_t         Centrality;
  vector<double>  *MCPrimaryVertex = new vector<double>;
  vector<double>  *RecoPrimaryVertex = new vector<double>;
  Float_t         MCJetPt;
  Float_t         MCJetEta;
  Float_t         MCJetPhi;
  Float_t         MCJetArea;
  Float_t         MCJetE;
  Int_t           MCJetNConst;
  Float_t         MCD0Pt;
  Float_t         MCD0Eta;
  Float_t         MCD0Phi;
  Float_t         MCPionPt;
  Float_t         MCPionEta;
  Float_t         MCPionPhi;
  Float_t         MCKaonPt;
  Float_t         MCKaonEta;
  Float_t         MCKaonPhi;
  Float_t         RecoJetPt;
  Float_t         RecoJetCorrPt;
  Float_t         RecoJetEta;
  Float_t         RecoJetPhi;
  Float_t         RecoJetArea;
  Float_t         RecoJetE;
  Float_t         RecoJetRhoVal;
  Int_t           RecoJetNConst;
  Float_t         RecoD0Pt;
  Float_t         RecoD0Eta;
  Float_t         RecoD0Phi;
  Float_t         RecoPionPt;
  Float_t         RecoPionEta;
  Float_t         RecoPionPhi;
  Float_t         RecoKaonPt;
  Float_t         RecoKaonEta;
  Float_t         RecoKaonPhi;

  JetTree->SetBranchAddress("Centrality", &Centrality);
  JetTree->SetBranchAddress("MCPrimaryVertex", &MCPrimaryVertex);
  JetTree->SetBranchAddress("RecoPrimaryVertex", &RecoPrimaryVertex);
  JetTree->SetBranchAddress("MCJetPt", &MCJetPt);
  JetTree->SetBranchAddress("MCJetEta", &MCJetEta);
  JetTree->SetBranchAddress("MCJetPhi", &MCJetPhi);
  JetTree->SetBranchAddress("MCJetArea", &MCJetArea);
  JetTree->SetBranchAddress("MCJetE", &MCJetE);
  JetTree->SetBranchAddress("MCJetNConst", &MCJetNConst);
  JetTree->SetBranchAddress("MCD0Pt", &MCD0Pt);
  JetTree->SetBranchAddress("MCD0Eta", &MCD0Eta);
  JetTree->SetBranchAddress("MCD0Phi", &MCD0Phi);
  JetTree->SetBranchAddress("MCPionPt", &MCPionPt);
  JetTree->SetBranchAddress("MCPionEta", &MCPionEta);
  JetTree->SetBranchAddress("MCPionPhi", &MCPionPhi);
  JetTree->SetBranchAddress("MCKaonPt", &MCKaonPt);
  JetTree->SetBranchAddress("MCKaonEta", &MCKaonEta);
  JetTree->SetBranchAddress("MCKaonPhi", &MCKaonPhi);
 
  JetTree->SetBranchAddress("RecoJetPt", &RecoJetPt);
  JetTree->SetBranchAddress("RecoJetCorrPt", &RecoJetCorrPt);
  JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
  JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
  JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
  JetTree->SetBranchAddress("RecoJetE", &RecoJetE);
  JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
  JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConst);
  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);


  TH1D *NEvents;
  TH1D *hVzMC[3];
  TH1D *hVzReco[3];
  TH1D *hVz[3];
  TH1D *hDeltaPt[3];
  TH1D *hDeltaEta[3];
  TH1D *hDeltaPhi[3];

  TH2D *hD0Pt[3];
  TH2D *hD0Eta[3];
  TH2D *hD0Phi[3];

  TH1D *hDeltaD0Pt;

  TH1D *hMCD0Pt[3];
  TH1D *hMCD0CrossSection[3];
  TH1D *hMCD0PtAverage[3];
  TH1D *hMCD0JetPt[3];
  TH1D *hMCD0PtWeighted[3];
  TH1D *hMCD0JetPtWeighted[3];

  TString CutNames[11] = {"1 < Reco D0 pT < 10", 
                          "1 < Reco Jet pT < 30", 
                          "Reco Jet |Eta| < 0.6", 
                          "0 < Reco Z < 1",
                          "1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6",
                          "1 < Reco D0 pT < 10 && 0 < Reco Z < 1 && Reco Jet |Eta| < 0.6",
                          "1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6 && 0 < Reco Z < 1",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6 && 1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6 && 1 < Reco D0 pT < 10 && 0 < Reco Z < 1 && Reco Jet |Eta| < 0.6",
                          "MC Jet |Eta| > 0.6 && Reco Jet |Eta| < 0.6 && 1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6 && 0 < Reco Z < 1"
                        };

  TString CutNamesForReco[11] = {"Reco Jet |Eta| < 0.6 && 1 < Reco D0 pT < 10 && 1 < Reco Jet pT < 30 && Reco Jet |Eta| < 0.6 && 0 < Reco Z < 1 && MC Jet |Eta| > 0.6"};

  TH1D *hMCD0JetPtBeforeFiducialCuts[3][11];
  TH1D *hMCD0JetPtAfterFiducialCuts[3][11];

  TH1D *hRecoD0JetPtBeforeFiducialCuts[3][11];
  TH1D *hRecoD0JetPtAfterFiducialCuts[3][11];

  TH1D *hMCZ[3];
  TH1D *hRecoZ[3];
  //Saving out a THn to sample jet candidates with specific values of MC z. 
  //The idea is to create response matrices with specific distributions of z. 
  //My assumption is that the underlying z distribution should not alter the response matrix drastically, 
  //should the entire phase space be sampled sensibly
  THn  *hD0MCZD0PtJetPtRecoD0PtJetPt[3]; //MC Z, MC D0 Pt, MC Jet Pt, Reco D0 Pt, Reco Jet Pt

  int nBinsDaug[nDimDaug] = {nBinsMCZ, nBinsMCZ, nBinsMCD0Pt, nBinsMCJetPt, nBinsRecoD0Pt, nBinsRecoJetPt, nBinsEta, nBinsEta}; //MC Z, Reco Z, MC D0 Pt, MC Jet Pt, Reco D0 Pt, Reco Jet Pt
  
  THn  *hMCJetPtMCZRecoJetPtRecoZ[3];

  int nBinsZHist[nDimZHist] = {nBinsJetPtForZHist, nBinsForZHist, nBinsJetPtForZHist, nBinsForZHist};


  TH2D *hResponsePt[3];
  TH2D *hResponsePtWeighed[3];


  RooUnfoldResponse *responseJetPtvZ[3];
  RooUnfoldResponse *responseJetPt[3];

  RooUnfoldResponse *resp2[3];

  TH2D *tmpMC = new TH2D("tmpMC", "tmpMC", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);
  TH2D *tmpReco = new TH2D("tmpReco", "tmpReco", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

  TH1D *tmpMC1D = new TH1D("tmpMC1D", "tmpMC1D", nBinsJetPtForZHist, JetPtBinsForZHist);
  TH1D *tmpReco1D = new TH1D("tmpReco1D", "tmpReco1D", nBinsJetPtForZHist, JetPtBinsForZHist);

  TH2D *resp2MC[3];
  TH2D *resp2Reco[3];
  THn *resp2THn[3];

  for (int i = 0; i < 3; i++){
    hVzMC[i] = new TH1D(Form("VzMC_%i", i), Form("VzMC_%i", i), 200, -20, 20); 
    hVzReco[i] = new TH1D(Form("VzReco_%i", i), Form("VzReco_%i", i), 200, -20, 20); 

    hVz[i] = new TH1D(Form("Vz_%i", i), Form("Vz_%i", i), 200, -20, 20); 
    hDeltaPt[i] = new TH1D(Form("DeltaPt_%i", i), Form("DeltaPt_%i", i), 100, -50, 50); 
    hDeltaEta[i] = new TH1D(Form("DeltaEta_%i", i), Form("DeltaEta_%i", i), 20, -1, 1);
    hDeltaPhi[i] = new TH1D(Form("DeltaPhi_%i", i), Form("DeltaPhi_%i", i), 100, -4, 4);

    hD0Pt[i] = new TH2D(Form("hD0Pt_%i", i), Form("hD0Pt_%i", i), 100, 0, 20, 100, 0, 20); //Reco v MC
    hD0Eta[i] = new TH2D(Form("hD0Eta_%i", i), Form("hD0Eta_%i", i), 100, -2, 2, 100, -2, 2);
    hD0Phi[i] = new TH2D(Form("hD0Phi_%i", i), Form("hD0Phi_%i", i), 100, -10, 10, 100, -10, 10);

    hMCD0Pt[i] = new TH1D(Form("MCD0Pt_%i", i), Form("MCD0Pt_%i", i), 20, 0., 10.); 
    hMCD0CrossSection[i] = new TH1D(Form("MCD0CrossSection_%i", i), Form("MCD0CrossSection_%i", i), 20, 0., 10.); 
    hMCD0PtAverage[i] = new TH1D(Form("MCD0PtAvg_%i", i), Form("MCD0PtAvg_%i", i), 20, 0., 10.); 
    hMCD0JetPt[i] = new TH1D(Form("MCD0JetPt_%i", i), Form("MCD0JetPt_%i", i), 100, 0., 50.); 
    hMCD0PtWeighted[i] = new TH1D(Form("MCD0PtWeighted_%i", i), Form("MCD0PtWeighted_%i", i), 20, 0., 10.); 
    hMCD0JetPtWeighted[i] = new TH1D(Form("MCD0JetPtWeighted_%i", i), Form("MCD0JetPtWeighted_%i", i), 100, 0., 50.);

    for (int cut = 0; cut < 11; cut++){
      hMCD0JetPtBeforeFiducialCuts[i][cut] = new TH1D(Form("hMCD0JetPtBeforeFiducialCuts_%i_%i", i, cut), Form("hMCD0JetPtBeforeFiducialCuts_%i_%i", i, cut), nBinsMCJetPt, MCJetPtBins);
      hMCD0JetPtAfterFiducialCuts[i][cut] = new TH1D(Form("hMCD0JetPtAfterFiducialCuts_%i_%i", i, cut), Form("hMCD0JetPtAfterFiducialCuts_%i_%i", i, cut), nBinsMCJetPt, MCJetPtBins);

      hRecoD0JetPtBeforeFiducialCuts[i][cut] = new TH1D(Form("hRecoD0JetPtBeforeFiducialCuts_%i_%i", i, cut), Form("hRecoD0JetPtBeforeFiducialCuts_%i_%i", i, cut), nBinsRecoJetPt, RecoJetPtBins);
      hRecoD0JetPtAfterFiducialCuts[i][cut] = new TH1D(Form("hRecoD0JetPtAfterFiducialCuts_%i_%i", i, cut), Form("hRecoD0JetPtAfterFiducialCuts_%i_%i", i, cut), nBinsRecoJetPt, RecoJetPtBins);
    }
    
    hMCZ[i] = new TH1D(Form("MCZ_%i", i), Form("MCZ_%i", i), nBinsMCZ, ZBins); 
    hRecoZ[i] = new TH1D(Form("RecoZ_%i", i), Form("RecoZ_%i", i), nBinsMCZ, ZBins); 

    hResponsePt[i] = new TH2D(Form("hResponsePt_%i", i), Form("hResponsePt_%i", i), 100, -50, 50, 100, -50, 50); //Reco v MC
    hResponsePtWeighed[i] = new TH2D(Form("hResponsePtWeighed_%i", i), Form("hResponsePtWeighed_%i", i), 100, -50, 50, 100, -50, 50); //Reco v MC
    
    responseJetPtvZ[i] = new RooUnfoldResponse(tmpMC, tmpReco, Form("ResponseJetPtvZ_%i", i), Form("ResponseJetPtvZ_%i", i));
    responseJetPt[i] = new RooUnfoldResponse(tmpMC1D, tmpReco1D, Form("ResponseJetPt_%i", i), Form("ResponseJetPt_%i", i));

    resp2[i] = new RooUnfoldResponse(tmpMC, tmpReco, Form("Response2NewDef_%i", i), Form("Response2NewDef_%i", i));
    resp2MC[i] = new TH2D(Form("resp2MC_%i", i), Form("resp2MC_%i", i), nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);
    resp2Reco[i] = new TH2D(Form("resp2Reco_%i", i), Form("resp2Reco_%i", i), nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

    resp2THn[i] = new THnF(Form("Response2THn_%i", i), Form("Response2THn_%i", i), nDimZHist, nBinsZHist, NULL, NULL);

    resp2THn[i]->SetBinEdges(0, JetPtBinsForZHist); //MCJetPtBins
    resp2THn[i]->SetBinEdges(1, ZBinsForZHist); //MCZBins
    resp2THn[i]->SetBinEdges(2, JetPtBinsForZHist); //RecoJetPtBins
    resp2THn[i]->SetBinEdges(3, ZBinsForZHist); //RecoZBins

    hD0MCZD0PtJetPtRecoD0PtJetPt[i] = new THnF(Form("hD0MCZD0PtJetPtRecoD0PtJetPt_%i", i), Form("hD0MCZD0PtJetPtRecoD0PtJetPt_%i", i), nDimDaug, nBinsDaug, NULL, NULL);

    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(0, ZBins); //MCBins
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(1, ZBins); //RecoBins
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(2, MCD0PtBins);
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(3, MCJetPtBins);
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(4, RecoD0PtBins);
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(5, RecoJetPtBins);
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(6, EtaBins);
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->SetBinEdges(7, EtaBins);

    hMCJetPtMCZRecoJetPtRecoZ[i] = new THnF(Form("hMCJetPtMCZRecoJetPtRecoZ_%i", i), Form("hMCJetPtMCZRecoJetPtRecoZ_%i", i), nDimZHist, nBinsZHist, NULL, NULL);

    hMCJetPtMCZRecoJetPtRecoZ[i]->SetBinEdges(0, JetPtBinsForZHist); //MCJetPtBins
    hMCJetPtMCZRecoJetPtRecoZ[i]->SetBinEdges(1, ZBinsForZHist); //MCZBins
    hMCJetPtMCZRecoJetPtRecoZ[i]->SetBinEdges(2, JetPtBinsForZHist); //RecoJetPtBins
    hMCJetPtMCZRecoJetPtRecoZ[i]->SetBinEdges(3, ZBinsForZHist); //RecoZBins

  }

  NEvents = new TH1D("NEvents", "NEvents", 3, -0.5, 2.5);
  hDeltaD0Pt = new TH1D("hDeltaD0Pt", "hDeltaD0Pt", 100, -20, 20);

  TFile *fonllweight = new TFile("./FONLLWeights_Mar26.root");
  fonllweight->cd();

  TH1D *FONLL = (TH1D *)fonllweight->Get("FONLLWeights");

  cout << JetTree->GetEntries() << endl;

  int nentries = JetTree->GetEntries();

  int goodentries = 0;

  TFile *tmpzfile = new TFile("TMPZFile.root");
  TH1D *tmpZ = (TH1D *)tmpzfile->Get("Z");

  for (int i = 0; i < nentries; i++){

    if (i%100000 ==0) cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    // cout << Centrality << "\t" << MCJetPt << "\t" << RecoJetCorrPt << endl;

    // if (abs(MCJetEta) > 1.) continue;
    if (abs(RecoPrimaryVertex->at(2)) > 6) continue;
    if (abs(MCPrimaryVertex->at(2)) > 6) continue;

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

    if (centhistogramtofill < 0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);

    // if (mcz > 1.0) cout << MCJetPt << "\t" << MCD0Pt << mcz << endl;

    double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);


    double w = tmpZ->GetBinContent(tmpZ->FindBin(mcz));

    //Defining response matrix using fakes and misses to check closure

    bool isMCD0Pt = MCD0Pt > 1 && MCD0Pt < 10;
    bool isRecoD0Pt = RecoD0Pt > 1 && RecoD0Pt < 10;
    bool isMCZ   = mcz > 0 && mcz < 1;
    bool isRecoZ = recoz > 0 && recoz < 1;
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 30;
    bool isRecoJetPt = RecoJetCorrPt > -7 && RecoJetCorrPt < 30;
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;

    // Response

    if (isMCD0Pt && isRecoD0Pt && isMCJetPt && isRecoJetPt && isMCJetEta && isRecoJetEta) {
      resp2[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
      resp2MC[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      resp2Reco[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);

      double vtofill[4] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
      resp2THn[centhistogramtofill]->Fill(vtofill, w);
    }

    // Fake

    if ((!isMCD0Pt || !isMCJetPt || !isMCJetEta) && isRecoD0Pt && isRecoJetPt && isRecoJetEta) {
      resp2[centhistogramtofill]->Fake(RecoJetCorrPt, recoz, w);
      resp2Reco[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
    }

    // Miss

    if (isMCD0Pt && isMCJetPt && isMCJetEta && (!isRecoD0Pt || !isRecoJetPt || !isRecoJetEta)) {
      resp2[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      resp2MC[centhistogramtofill]->Fill(MCJetPt, mcz, w);
    }

    if (MCD0Pt > 10) continue;
    if (MCJetNConst <= 1) continue;

    // Here I am going to fill up the cut histograms for QA

    if (abs(MCJetEta)<0.6 && MCJetPt >= 1. && MCJetPt <= 30){
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][0]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][1]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][2]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][3]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][4]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][5]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][6]->Fill(MCJetPt);

      if (RecoD0Pt > 1.0 && RecoD0Pt < 10.0)                                                                                                     hMCD0JetPtAfterFiducialCuts[centhistogramtofill][0]->Fill(MCJetPt);
      if (RecoJetCorrPt > 1.0 && RecoJetCorrPt < 30.0)                                                                                           hMCD0JetPtAfterFiducialCuts[centhistogramtofill][1]->Fill(MCJetPt);
      if (abs(RecoJetEta) < 0.6)                                                                                                                 hMCD0JetPtAfterFiducialCuts[centhistogramtofill][2]->Fill(MCJetPt);
      if (recoz > 0 && recoz < 1)                                                                                                                hMCD0JetPtAfterFiducialCuts[centhistogramtofill][3]->Fill(MCJetPt);
      if ((RecoD0Pt > 1.0 && RecoD0Pt < 10.0)&&(RecoJetCorrPt > 1.0 && RecoJetCorrPt < 30.0)&&(abs(RecoJetEta) < 0.6))                           hMCD0JetPtAfterFiducialCuts[centhistogramtofill][4]->Fill(MCJetPt);
      if ((RecoD0Pt > 1.0 && RecoD0Pt < 10.0)&&(recoz > 0 && recoz < 1)&&(abs(RecoJetEta) < 0.6))                                                hMCD0JetPtAfterFiducialCuts[centhistogramtofill][5]->Fill(MCJetPt);
      if ((RecoD0Pt > 1.0 && RecoD0Pt < 10.0)&&(RecoJetCorrPt > 1.0 && RecoJetCorrPt < 30.0)&&(abs(RecoJetEta) < 0.6)&&(recoz > 0 && recoz < 1)) hMCD0JetPtAfterFiducialCuts[centhistogramtofill][6]->Fill(MCJetPt);
    }

    if (abs(MCJetEta)>=0.6 && MCJetPt >= 1. && MCJetPt <= 30){
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][7]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][8]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][9]->Fill(MCJetPt);
      hMCD0JetPtBeforeFiducialCuts[centhistogramtofill][10]->Fill(MCJetPt);
      if (abs(RecoJetEta) < 0.6)                                                                                                                 hMCD0JetPtAfterFiducialCuts[centhistogramtofill][7]->Fill(MCJetPt);
      if ((RecoD0Pt > 1.0 && RecoD0Pt < 10.0)&&(RecoJetCorrPt > 1.0 && RecoJetCorrPt < 30.0)&&(abs(RecoJetEta) < 0.6))                           hMCD0JetPtAfterFiducialCuts[centhistogramtofill][8]->Fill(MCJetPt);
      if ((RecoD0Pt > 1.0 && RecoD0Pt < 10.0)&&(recoz > 0 && recoz < 1)&&(abs(RecoJetEta) < 0.6))                                                hMCD0JetPtAfterFiducialCuts[centhistogramtofill][9]->Fill(MCJetPt);
      if ((RecoD0Pt > 1.0 && RecoD0Pt < 10.0)&&(RecoJetCorrPt > 1.0 && RecoJetCorrPt < 30.0)&&(abs(RecoJetEta) < 0.6)&&(recoz > 0 && recoz < 1)) hMCD0JetPtAfterFiducialCuts[centhistogramtofill][10]->Fill(MCJetPt);
    }

    if ((RecoD0Pt > 1.0 && RecoD0Pt < 10.0)&&(RecoJetCorrPt > 1.0 && RecoJetCorrPt < 30.0)&&(abs(RecoJetEta) < 0.6)&&(recoz > 0 && recoz < 1)){
      hRecoD0JetPtBeforeFiducialCuts[centhistogramtofill][0]->Fill(RecoJetCorrPt);
      if (abs(MCJetEta)<0.6 && MCJetPt >= 1. && MCJetPt <= 30)                                                                                  hRecoD0JetPtAfterFiducialCuts[centhistogramtofill][0]->Fill(RecoJetCorrPt);
    }

    // Done filling up cut histograms for QA

    if (abs(MCJetEta)> 0.6) continue; //This is a little sus. We will come back to this.
    if (MCJetPt < 1.) continue;
    if (MCJetPt > 30.) continue;
    

    // hMCD0JetPtBeforeFiducialCuts[centhistogramtofill]->Fill(MCJetPt);
    if (abs(RecoJetEta)> 0.6) continue;

    double zhisttofill[4] = {MCJetPt, mcz, RecoJetCorrPt, recoz};

    if (recoz > -2 && recoz < 3. && RecoJetCorrPt > -20. && RecoJetCorrPt < 30.) {

      hMCJetPtMCZRecoJetPtRecoZ[centhistogramtofill]->Fill(zhisttofill);
      responseJetPtvZ[centhistogramtofill]->Fill(zhisttofill[2], zhisttofill[3], zhisttofill[0], zhisttofill[1]);
      responseJetPt[centhistogramtofill]->Fill(zhisttofill[2], zhisttofill[0]);

    }

    // Different definitions of Response Matrices

    if (RecoJetCorrPt < 1.) continue;
    if (RecoJetCorrPt > 30.) continue;
    if (recoz < 0. || recoz > 1.) continue;

    double deltaR = dR(dPhi(MCJetPhi, RecoJetPhi), dEta(MCJetEta, RecoJetEta));
    // if (deltaR > 0.4) continue;

    double dVz = RecoPrimaryVertex->at(2) - MCPrimaryVertex->at(2);
    double deltaptofjets = RecoJetCorrPt - MCJetPt;
    double deltaetaofjets = RecoJetEta - MCJetEta;
    double deltaphiofjets = dPhi(MCJetPhi, RecoJetPhi);

    // double w;
    // int mcbinx = FONLL->GetXaxis()->FindBin(MCJetPt);
    // w = FONLL->GetBinContent(mcbinx);

    // hMCD0JetPtAfterFiducialCuts[centhistogramtofill]->Fill(MCJetPt);

    // if (abs(dVz) < 3 || abs(dVz) > 6) continue;
    // if (abs(dVz) > 6) continue;

    goodentries++;

    // cout << deltaptofjets << endl;


    TVector3 RecoD0, MCD0;
    RecoD0.SetPtEtaPhi(RecoD0Pt, RecoD0Eta, RecoD0Phi);
    MCD0.SetPtEtaPhi(MCD0Pt, MCD0Eta, MCD0Phi);
    // hDeltaD0Pt->Fill(sqrt(pow(RecoD0.X() - MCD0.X(), 2) + pow(RecoD0.Y() - MCD0.Y(), 2)));
    hDeltaD0Pt->Fill(RecoD0.Pt() - MCD0.Pt());

    NEvents->Fill(centhistogramtofill);

    hVzMC[centhistogramtofill]->Fill(MCPrimaryVertex->at(2));
    hVzReco[centhistogramtofill]->Fill(RecoPrimaryVertex->at(2));
    hVz[centhistogramtofill]->Fill(dVz);

    hDeltaPt[centhistogramtofill]->Fill(deltaptofjets);
    hDeltaEta[centhistogramtofill]->Fill(deltaetaofjets);
    hDeltaPhi[centhistogramtofill]->Fill(deltaphiofjets);

    hD0Pt[centhistogramtofill]->Fill(RecoD0Pt, MCD0Pt);
    hD0Eta[centhistogramtofill]->Fill(RecoD0Eta, MCD0Eta);
    hD0Phi[centhistogramtofill]->Fill(RecoD0Phi, MCD0Phi);

    hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt, 1./MCD0Pt);
    hMCD0CrossSection[centhistogramtofill]->Fill(MCD0Pt, 1./MCD0Pt);
    hMCD0PtAverage[centhistogramtofill]->Fill(MCD0Pt, MCD0Pt);
    hMCD0JetPt[centhistogramtofill]->Fill(MCJetPt);
    // hMCD0PtWeighted[centhistogramtofill]->Fill(MCD0Pt, w);
    // hMCD0JetPtWeighted[centhistogramtofill]->Fill(MCJetPt, w);

    hMCZ[centhistogramtofill]->Fill(mcz);
    hRecoZ[centhistogramtofill]->Fill(recoz);
    // cout << mcz << endl;

    double toFill[8] = {mcz, recoz, MCD0Pt, MCJetPt, RecoD0Pt, RecoJetCorrPt, MCJetEta, RecoJetEta };

    hD0MCZD0PtJetPtRecoD0PtJetPt[centhistogramtofill]->Fill(toFill);

    hResponsePt[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt);
    // hResponsePtWeighed[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);
  }

  cout << "Number of entries processed = " << goodentries << "\t" << double(goodentries)/double(nentries)*100 << endl;

  TFile *outfile = new TFile("Plots.root", "UPDATE");
  outfile->mkdir(Form("pthatbin_%i", pthatbin));
  outfile->cd(Form("pthatbin_%i", pthatbin));

  NEvents->Write();
  hDeltaD0Pt->Write();

  for (int i = 0; i < 3; i++){
    hVzMC[i]->Write();
    hVzReco[i]->Write();
    hVz[i]->Write();

    hDeltaPt[i]->Write();
    hDeltaEta[i]->Write();
    hDeltaPhi[i]->Write();

    hD0Pt[i]->Write();
    hD0Eta[i]->Write();
    hD0Phi[i]->Write();

    hMCD0Pt[i]->Write();
    hMCD0CrossSection[i]->Write();
    hMCD0PtAverage[i]->Write();
    hMCD0JetPt[i]->Write();
    hMCD0PtWeighted[i]->Write();
    hMCD0JetPtWeighted[i]->Write();

    for (int cut = 0; cut < 11; cut++){
      hMCD0JetPtBeforeFiducialCuts[i][cut]->Write();
      hMCD0JetPtAfterFiducialCuts[i][cut]->Write();

      hRecoD0JetPtBeforeFiducialCuts[i][cut]->Write();
      hRecoD0JetPtAfterFiducialCuts[i][cut]->Write();
    }
    

    hMCZ[i]->Write();
    hRecoZ[i]->Write();
    hD0MCZD0PtJetPtRecoD0PtJetPt[i]->Write();

    hMCJetPtMCZRecoJetPtRecoZ[i]->Write();

    hResponsePt[i]->Write();
    hResponsePtWeighed[i]->Write();

    gDirectory->WriteObject(responseJetPt[i], Form("ResponseJetPt_%i", i));
    gDirectory->WriteObject(responseJetPtvZ[i], Form("ResponseJetPtvZ_%i", i));

    gDirectory->WriteObject(resp2[i], Form("Response2Reco_%i", i));
    resp2MC[i]->Write();
    resp2Reco[i]->Write();

    resp2THn[i]->Write();
  }

  outfile->Close();

  JetTree->Reset();
  f->Close();
}

void PlotMethod(){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *f = new TFile("Plots.root");
  TH1D *NEntries[2];
  TH1D *VzMC[3][2];
  TH1D *MCD0Pt[3][2];
  TH1D *MCD0CrossSection[3][2];
  TH1D *MCD0PtAverage[3][2];
  TH1D *MCD0JetPt[3][2];
  TH1D *MCZ[3][2];
  TH1D *RecoZ[3][2];
  TH2D *ResponsePt[3][2];

  // double entries[2] = {8701140., 27064505.};

  // [0-3]:   177037/1000000
  // [3-inf]: 364262/1000000

  // double ratioofeventsthatpassveto[2] = {177037./1000000., 364262./1000000.}; // For scaling, we need to make sure that vetoed events are also included.
  // double crosssection[2] = {0.103/0.109, 0.006/0.109};

  // double ratioofeventsthatpassveto[2] = {100./405., 100./284.}; // For scaling, we need to make sure that vetoed events are also included.
  // double crosssection[2] = {0.099, 0.006};

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

      MCD0Pt[i][bin] = (TH1D *)gDirectory->Get(Form("MCD0Pt_%i", i));
      MCD0CrossSection[i][bin] = (TH1D *)gDirectory->Get(Form("MCD0CrossSection_%i", i));
      MCD0PtAverage[i][bin] = (TH1D *)gDirectory->Get(Form("MCD0PtAvg_%i", i));
      MCD0JetPt[i][bin] = (TH1D *)gDirectory->Get(Form("MCD0JetPt_%i", i));
      MCZ[i][bin] = (TH1D *)gDirectory->Get(Form("MCZ_%i", i));
      RecoZ[i][bin] = (TH1D *)gDirectory->Get(Form("RecoZ_%i", i));
      ResponsePt[i][bin] = (TH2D *)gDirectory->Get(Form("hResponsePt_%i", i));

      // cout << "Imported the histograms" << endl;

      MCD0Pt[i][bin]->Scale(scale);
      MCD0CrossSection[i][bin]->Scale(scale);
      MCD0PtAverage[i][bin]->Scale(scale);
      MCD0JetPt[i][bin]->Scale(scale);
      MCZ[i][bin]->Scale(scale);
      RecoZ[i][bin]->Scale(scale);
      ResponsePt[i][bin]->Scale(scale);
    }
  }

  double branchratio = 0.0389;

  for (int i = 0; i < 3; i++){
    MCD0Pt[i][0]->Add(MCD0Pt[i][1]);
    MCD0CrossSection[i][0]->Add(MCD0CrossSection[i][1]);
    MCD0PtAverage[i][0]->Add(MCD0PtAverage[i][1]);
    MCD0JetPt[i][0]->Add(MCD0JetPt[i][1]);
    MCZ[i][0]->Add(MCZ[i][1]);
    RecoZ[i][0]->Add(RecoZ[i][1]);
    ResponsePt[i][0]->Add(ResponsePt[i][1]);
    MCD0Pt[i][0]->Scale(30./(0.565*0.5)); // This ensures the output is in millibarn
    MCD0CrossSection[i][0]->Scale(30./(2*2.*TMath::Pi()*branchratio*0.5)); // This ensures the output is in millibarn
  }

  // for (int i = 0; i < 3; i++){
  //   MCD0PtAverage[i][0]->Divide(MCD0Pt[i][0]);
  //   for (int bin = 1; bin <= MCD0PtAverage[i][0]->GetNbinsX(); bin++){
  //     cout << bin << "\t" <<  MCD0PtAverage[i][0]->GetBinCenter(bin) << "\t" << MCD0PtAverage[i][0]->GetBinContent(bin) << endl;
  //   }
  // }

  
  TH1F *FONLL = new TH1F("FONLL","FONLL",60,0,30);
  TH1F *MinFONLL = new TH1F("MIN FONLL","MIN FONLL",60,0,30);
  TH1F *MaxFONLL = new TH1F("MAX FONLL","MAX FONLL",60,0,30);
  ifstream data1("FONLL_New.txt");
  int cnt=1;
  if(data1.is_open()){
      while(!data1.eof()){
          double x;double y;double minerr; double maxerr; double tmp;
          data1 >> x >> y >> minerr >> maxerr >> tmp >> tmp >> tmp >> tmp;
          FONLL->SetBinContent(cnt,y*pow(10, -9)/x);
          MinFONLL->SetBinContent(cnt,minerr*pow(10, -9)/x);
          MaxFONLL->SetBinContent(cnt,maxerr*pow(10, -9)/x);

          // FONLL->SetBinContent(cnt,y/(7.5330e+07));
          // MinFONLL->SetBinContent(cnt,minerr/(4.4410e+07));
          // MaxFONLL->SetBinContent(cnt,maxerr/(1.7240e+08));

          cnt++;
      }
  }

  // for(int i = 1; i <= 61; i++){
  //   cout << i << "\t" << FONLL->GetBinCenter(i) << endl;
  // }

  // FONLL->Scale(1./FONLL->Integral());
  FONLL->SetLineColor(kRed);
  FONLL->SetMarkerColor(kRed);
  FONLL->SetMarkerStyle(20);
  FONLL->SetMarkerSize(1);

  // MinFONLL->Scale(1./MinFONLL->Integral());
  MinFONLL->SetLineColor(kGreen-2);
  MinFONLL->SetMarkerColor(kGreen-2);
  MinFONLL->SetMarkerStyle(24);

  // MaxFONLL->Scale(1./MaxFONLL->Integral());
  MaxFONLL->SetLineColor(kBlack);
  MaxFONLL->SetMarkerColor(kBlack);
  MaxFONLL->SetMarkerStyle(24);

  MaxFONLL->GetYaxis()->SetTitle("#frac{1}{p_{T}} #frac{d#sigma}{dp_{T}} [mb/GeV]");
  MaxFONLL->GetXaxis()->SetTitle("p_{T}");

  double ptxaxis_published[5] = {1.57, 2.45, 3.44, 4.45, 5.45};
  double dsigmadpt_published[5] = {0.106089296, 0.03107972, 0.00756112, 0.0016320464, 0.00065097852};
  double dsigmadpt_fonll[5] = {1.0515E-01, 3.5115E-02, 9.5505E-03, 2.8380E-03, 9.9716E-04};

  double dsigmaptdpt_published[5];
  double dsigmaptdpt_fonll[5];

  for (int i = 0; i < 5; i++){
      dsigmaptdpt_published[i] = dsigmadpt_published[i]/ptxaxis_published[i];
      dsigmaptdpt_fonll[i] = dsigmadpt_fonll[i]/ptxaxis_published[i];
    }

  TGraph *publishedspectra = new TGraph(6, ptxaxis_published, dsigmaptdpt_published);
  TGraph *publishedfonll = new TGraph(6, ptxaxis_published, dsigmaptdpt_fonll);
  publishedspectra->SetLineColor(kRed);
  publishedspectra->SetMarkerColor(kRed);
  publishedspectra->SetMarkerStyle(22);

  publishedfonll->SetLineColor(kRed);
  publishedfonll->SetMarkerColor(kRed);
  publishedfonll->SetMarkerStyle(26);

  auto legend = new TLegend(0.5,0.7,0.9,0.9);
  legend->AddEntry(MaxFONLL,"FONLL","f");
  legend->AddEntry(MCD0Pt[0][0],"PYTHIA D^{0}/0.565","p");
  legend->AddEntry(publishedspectra,"D^{0} Run 09","p");
  // legend->AddEntry(MaxFONLL,"FONLL","f");


  TCanvas *a = new TCanvas("D0 Pt", "D0 Pt", 1600, 400);
  a->Divide(3);
  for (int i = 0; i < 3; i++){
    a->cd(i+1);
    gPad->SetLogy();
    MaxFONLL->GetYaxis()->SetRangeUser(5*pow(10, -6), 2);
    MaxFONLL->GetXaxis()->SetRangeUser(0.5, 10);
    // MCD0CrossSection[i][0]->Draw();
    

    MaxFONLL->Draw("C SAME");
    MaxFONLL->SetFillColor(kGreen);
    MaxFONLL->SetLineColor(kGreen);

    MinFONLL->Draw("C SAME");
    MinFONLL->SetFillColor(10);
    MinFONLL->SetLineColor(kGreen);

    // FONLL->Draw("P SAME");
    
    publishedspectra->Draw("P SAME");
    // publishedfonll->Draw("P SAME");

    MCD0Pt[i][0]->Draw("EP SAME");
    MCD0Pt[i][0]->SetMarkerColor(kBlack);
    MCD0Pt[i][0]->SetLineColor(kBlack);
    MCD0Pt[i][0]->SetMarkerStyle(20);
    // MCD0Pt[i][0]->SetMarkerSize(2);

    legend->Draw("SAME");     
  }

  a->SaveAs("D0Pt.pdf");

  
  TCanvas *b = new TCanvas("D0 Jet Pt", "D0 Jet Pt", 1600, 400);
  b->Divide(3);
  for (int i = 0; i < 3; i++){
    b->cd(i+1);
    gPad->SetLogy();
    MCD0JetPt[i][0]->Draw();
  }

  b->SaveAs("JetPt.pdf");

  TCanvas *c = new TCanvas("Response", "Response", 1600, 400);
  c->Divide(3);
  for (int i = 0; i < 3; i++){
    c->cd(i+1);
    gPad->SetLogz();
    ResponsePt[i][0]->Draw("COLZ");
    ResponsePt[i][0]->GetYaxis()->SetRangeUser(0, 50);
    ResponsePt[i][0]->GetZaxis()->SetRangeUser(pow(10, -10), pow(10, -4));
  }

  c->SaveAs("Response.pdf");

  double integral = MCZ[0][0]->Integral();
  // cout << integral << endl;
  TF1 *f1 = new TF1("f1","pow(1-x, 4)*0.00109084",0, 1);
  

  TCanvas *d = new TCanvas("Z", "Z", 1600,400);
  d->Divide(3);
  for (int i = 0; i < 3; i++){
    d->cd(i+1);
    gPad->SetLogy();
    f1->SetMarkerColor(kBlack);
    f1->SetMarkerStyle(25);
    
    MCZ[i][0]->Draw("SAME");
    f1->Draw("SAME");
    RecoZ[i][0]->Draw("SAME");
    MCZ[i][0]->SetLineColor(kRed);
    MCZ[i][0]->SetMarkerColor(kRed);
    MCZ[i][0]->SetMarkerSize(20);

    RecoZ[i][0]->SetLineColor(kBlack);
    RecoZ[i][0]->SetMarkerColor(kBlack);
    RecoZ[i][0]->SetMarkerSize(20);
  }

  
}

void HIOverlay(){
    
  TFile *outfile = new TFile("Plots.root", "RECREATE");
  outfile->Close();
  Method(0);
  Method(1);

  PlotMethod();
}