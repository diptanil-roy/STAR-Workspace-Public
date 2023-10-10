#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"
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

const double Mpion = 0.139570;
const double Mkaon = 0.493677;
const double Mproton = 0.938272;

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

void EfficiencyPlotter(){

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0); 


  // TFile *infile = TFile::Open("/Users/diptanilroy/Desktop/Y2020/STAR/STARAnalysis/2021/EMBEDDING/RAGHAVSFRAMEWORK/ResponsePlots/SampleEfficiency.root");
  TFile *infile = TFile::Open("EmbeddingQA.root");
  // TFile *infile = TFile::Open("/Users/diptanilroy/Desktop/Y2020/STAR/STARAnalysis/2021/2022/Simulation/PythiaSIM_3_0_NoD0Decay_Feb7.root");
  // infile->cd("DummyMaker");

  TH1F *hD0Mass = (TH1F *)gDirectory->Get("hD0Mass");

  TH1F *hMCD0KaonPt = (TH1F *)gDirectory->Get("hMCD0KaonPt");
  TH1F *hMCD0PionPt = (TH1F *)gDirectory->Get("hMCD0PionPt");

  TH1F *hRecoD0KaonPt = (TH1F *)gDirectory->Get("hRecoD0KaonPt");
  TH1F *hRecoD0PionPt = (TH1F *)gDirectory->Get("hRecoD0PionPt");

  TH2F *hRecoVMCD0KaonPt = (TH2F *)gDirectory->Get("hRecoVMCD0KaonPt");
  TH2F *hRecoVMCD0PionPt = (TH2F *)gDirectory->Get("hRecoVMCD0PionPt");

  TH1F *hMCD0KaonPhi = (TH1F *)gDirectory->Get("hMCD0KaonPhi");
  TH1F *hMCD0PionPhi = (TH1F *)gDirectory->Get("hMCD0PionPhi");

  TH1F *hRecoD0KaonPhi = (TH1F *)gDirectory->Get("hRecoD0KaonPhi");
  TH1F *hRecoD0PionPhi = (TH1F *)gDirectory->Get("hRecoD0PionPhi");

  TH2F *hRecoVMCD0KaonPhi = (TH2F *)gDirectory->Get("hRecoVMCD0KaonPhi");
  TH2F *hRecoVMCD0PionPhi = (TH2F *)gDirectory->Get("hRecoVMCD0PionPhi");

  TH1F *hMCD0KaonEta = (TH1F *)gDirectory->Get("hMCD0KaonEta");
  TH1F *hMCD0PionEta = (TH1F *)gDirectory->Get("hMCD0PionEta");

  TH1F *hRecoD0KaonEta = (TH1F *)gDirectory->Get("hRecoD0KaonEta");
  TH1F *hRecoD0PionEta = (TH1F *)gDirectory->Get("hRecoD0PionEta");

  TH2F *hRecoVMCD0KaonEta = (TH2F *)gDirectory->Get("hRecoVMCD0KaonEta");
  TH2F *hRecoVMCD0PionEta = (TH2F *)gDirectory->Get("hRecoVMCD0PionEta");

  hMCD0KaonPt->Rebin(2);
  hMCD0PionPt->Rebin(2);
  hRecoD0KaonPt->Rebin(2);
  hRecoD0PionPt->Rebin(2);
  hMCD0KaonPhi->Rebin(2);
  hMCD0PionPhi->Rebin(2);
  hRecoD0KaonPhi->Rebin(2);
  hRecoD0PionPhi->Rebin(2);
  hMCD0KaonEta->Rebin(2);
  hMCD0PionEta->Rebin(2);
  hRecoD0KaonEta->Rebin(2);
  hRecoD0PionEta->Rebin(2);

  ////////////////////////// All Kaon Pion Histograms //////////////////////////////////////////

  TH1F *hMCKaonPt = (TH1F *)gDirectory->Get("hMCKaonPt");
  TH1F *hMCPionPt = (TH1F *)gDirectory->Get("hMCPionPt");

  TH1F *hRecoKaonPt = (TH1F *)gDirectory->Get("hRecoKaonPt");
  TH1F *hRecoPionPt = (TH1F *)gDirectory->Get("hRecoPionPt");

  TH2F *hRecoVMCKaonPt = (TH2F *)gDirectory->Get("hRecoVMCKaonPt");
  TH2F *hRecoVMCPionPt = (TH2F *)gDirectory->Get("hRecoVMCPionPt");

  TH1F *hMCKaonPhi = (TH1F *)gDirectory->Get("hMCKaonPhi");
  TH1F *hMCPionPhi = (TH1F *)gDirectory->Get("hMCPionPhi");

  TH1F *hRecoKaonPhi = (TH1F *)gDirectory->Get("hRecoKaonPhi");
  TH1F *hRecoPionPhi = (TH1F *)gDirectory->Get("hRecoPionPhi");

  TH2F *hRecoVMCKaonPhi = (TH2F *)gDirectory->Get("hRecoVMCKaonPhi");
  TH2F *hRecoVMCPionPhi = (TH2F *)gDirectory->Get("hRecoVMCPionPhi");

  TH1F *hMCKaonEta = (TH1F *)gDirectory->Get("hMCKaonEta");
  TH1F *hMCPionEta = (TH1F *)gDirectory->Get("hMCPionEta");

  TH1F *hRecoKaonEta = (TH1F *)gDirectory->Get("hRecoKaonEta");
  TH1F *hRecoPionEta = (TH1F *)gDirectory->Get("hRecoPionEta");

  TH2F *hRecoVMCKaonEta = (TH2F *)gDirectory->Get("hRecoVMCKaonEta");
  TH2F *hRecoVMCPionEta = (TH2F *)gDirectory->Get("hRecoVMCPionEta");

  hMCKaonPt->Rebin(2);
  hMCPionPt->Rebin(2);
  hRecoKaonPt->Rebin(2);
  hRecoPionPt->Rebin(2);
  hMCKaonPhi->Rebin(2);
  hMCPionPhi->Rebin(2);
  hRecoKaonPhi->Rebin(2);
  hRecoPionPhi->Rebin(2);
  hMCKaonEta->Rebin(2);
  hMCPionEta->Rebin(2);
  hRecoKaonEta->Rebin(2);
  hRecoPionEta->Rebin(2);

  ////////////////////////// All Kaon Pion (TPC) PID Histograms //////////////////////////////////////////

  TH1F *hRecoTPCKaonPt = (TH1F *)gDirectory->Get("hRecoTPCKaonPt");
  TH1F *hRecoTPCPionPt = (TH1F *)gDirectory->Get("hRecoTPCPionPt");

  TH2F *hRecoVMCTPCKaonPt = (TH2F *)gDirectory->Get("hRecoVMCTPCKaonPt");
  TH2F *hRecoVMCTPCPionPt = (TH2F *)gDirectory->Get("hRecoVMCTPCPionPt");

  TH1F *hRecoTPCKaonPhi = (TH1F *)gDirectory->Get("hRecoTPCKaonPhi");
  TH1F *hRecoTPCPionPhi = (TH1F *)gDirectory->Get("hRecoTPCPionPhi");

  TH2F *hRecoVMCTPCKaonPhi = (TH2F *)gDirectory->Get("hRecoVMCTPCKaonPhi");
  TH2F *hRecoVMCTPCPionPhi = (TH2F *)gDirectory->Get("hRecoVMCTPCPionPhi");

  TH1F *hRecoTPCKaonEta = (TH1F *)gDirectory->Get("hRecoTPCKaonEta");
  TH1F *hRecoTPCPionEta = (TH1F *)gDirectory->Get("hRecoTPCPionEta");

  TH2F *hRecoVMCTPCKaonEta = (TH2F *)gDirectory->Get("hRecoVMCTPCKaonEta");
  TH2F *hRecoVMCTPCPionEta = (TH2F *)gDirectory->Get("hRecoVMCTPCPionEta");

  hRecoTPCKaonPt->Rebin(2);
  hRecoTPCPionPt->Rebin(2);
  hRecoTPCKaonPhi->Rebin(2);
  hRecoTPCPionPhi->Rebin(2);
  hRecoTPCKaonEta->Rebin(2);
  hRecoTPCPionEta->Rebin(2);

  ////////////////////////// All Kaon Pion (TPC + TOF) PID Histograms //////////////////////////////////////////

  TH1F *hRecoTOFKaonPt = (TH1F *)gDirectory->Get("hRecoTOFKaonPt");
  TH1F *hRecoTOFPionPt = (TH1F *)gDirectory->Get("hRecoTOFPionPt");

  TH2F *hRecoVMCTOFKaonPt = (TH2F *)gDirectory->Get("hRecoVMCTOFKaonPt");
  TH2F *hRecoVMCTOFPionPt = (TH2F *)gDirectory->Get("hRecoVMCTOFPionPt");

  TH1F *hRecoTOFKaonPhi = (TH1F *)gDirectory->Get("hRecoTOFKaonPhi");
  TH1F *hRecoTOFPionPhi = (TH1F *)gDirectory->Get("hRecoTOFPionPhi");

  TH2F *hRecoVMCTOFKaonPhi = (TH2F *)gDirectory->Get("hRecoVMCTOFKaonPhi");
  TH2F *hRecoVMCTOFPionPhi = (TH2F *)gDirectory->Get("hRecoVMCTOFPionPhi");

  TH1F *hRecoTOFKaonEta = (TH1F *)gDirectory->Get("hRecoTOFKaonEta");
  TH1F *hRecoTOFPionEta = (TH1F *)gDirectory->Get("hRecoTOFPionEta");

  TH2F *hRecoVMCTOFKaonEta = (TH2F *)gDirectory->Get("hRecoVMCTOFKaonEta");
  TH2F *hRecoVMCTOFPionEta = (TH2F *)gDirectory->Get("hRecoVMCTOFPionEta");

  hRecoTOFKaonPt->Rebin(2);
  hRecoTOFPionPt->Rebin(2);
  hRecoTOFKaonPhi->Rebin(2);
  hRecoTOFPionPhi->Rebin(2);
  hRecoTOFKaonEta->Rebin(2);
  hRecoTOFPionEta->Rebin(2);

  ////////////////////////// All Kaon Pion (TPC + HFT) PID Histograms //////////////////////////////////////////

  TH1F *hRecoHFTKaonPt = (TH1F *)gDirectory->Get("hRecoHFTKaonPt");
  TH1F *hRecoHFTPionPt = (TH1F *)gDirectory->Get("hRecoHFTPionPt");

  TH2F *hRecoVMCHFTKaonPt = (TH2F *)gDirectory->Get("hRecoVMCHFTKaonPt");
  TH2F *hRecoVMCHFTPionPt = (TH2F *)gDirectory->Get("hRecoVMCHFTPionPt");

  TH1F *hRecoHFTKaonPhi = (TH1F *)gDirectory->Get("hRecoHFTKaonPhi");
  TH1F *hRecoHFTPionPhi = (TH1F *)gDirectory->Get("hRecoHFTPionPhi");

  TH2F *hRecoVMCHFTKaonPhi = (TH2F *)gDirectory->Get("hRecoVMCHFTKaonPhi");
  TH2F *hRecoVMCHFTPionPhi = (TH2F *)gDirectory->Get("hRecoVMCHFTPionPhi");

  TH1F *hRecoHFTKaonEta = (TH1F *)gDirectory->Get("hRecoHFTKaonEta");
  TH1F *hRecoHFTPionEta = (TH1F *)gDirectory->Get("hRecoHFTPionEta");

  TH2F *hRecoVMCHFTKaonEta = (TH2F *)gDirectory->Get("hRecoVMCHFTKaonEta");
  TH2F *hRecoVMCHFTPionEta = (TH2F *)gDirectory->Get("hRecoVMCHFTPionEta");

  hRecoHFTKaonPt->Rebin(2);
  hRecoHFTPionPt->Rebin(2);
  hRecoHFTKaonPhi->Rebin(2);
  hRecoHFTPionPhi->Rebin(2);
  hRecoHFTKaonEta->Rebin(2);
  hRecoHFTPionEta->Rebin(2);

  TH1F *RecoByMC[5][2][3];
  TString RatioType[5] = {"RecoD0ByMCD0", "RecoByMC", "RecoTPCByMC", "RecoTOFByMC", "RecoHFTvRecoTPC"};
  TString ParticleType[2] = {"Kaon", "Pion"};
  TString PlotType[3] = {"pT", "eta", "phi"};
  int numberofbins[3] = {16, 20, 40};
  double lowerbound[3] = {-0.5, -1.0, -10.};
  double upperbound[3] = {15.5, 1.0, 10.};

  for (int i = 0; i < 5; i++){
    for (int part = 0; part < 2; part++){
      for (int j = 0; j < 3; j++){
        TString Name = RatioType[i] + "_" + ParticleType[part] + "_" + PlotType[j];
        RecoByMC[i][part][j] = new TH1F(Name.Data(), Name.Data(), numberofbins[j], lowerbound[j], upperbound[j]);
      }
    }
  }

  //// D0 Reco Bois /////
  RecoByMC[0][0][0]->Add(hRecoD0KaonPt);
  RecoByMC[0][0][1]->Add(hRecoD0KaonEta);
  RecoByMC[0][0][2]->Add(hRecoD0KaonPhi);

  RecoByMC[0][1][0]->Add(hRecoD0PionPt);
  RecoByMC[0][1][1]->Add(hRecoD0PionEta);
  RecoByMC[0][1][2]->Add(hRecoD0PionPhi);

  //// All Reco Bois Geant Pid ////

  RecoByMC[1][0][0]->Add(hRecoKaonPt);
  RecoByMC[1][0][1]->Add(hRecoKaonEta);
  RecoByMC[1][0][2]->Add(hRecoKaonPhi);

  RecoByMC[1][1][0]->Add(hRecoPionPt);
  RecoByMC[1][1][1]->Add(hRecoPionEta);
  RecoByMC[1][1][2]->Add(hRecoPionPhi);

  //// All Reco Bois TPC PID /////

  RecoByMC[2][0][0]->Add(hRecoTPCKaonPt);
  RecoByMC[2][0][1]->Add(hRecoTPCKaonEta);
  RecoByMC[2][0][2]->Add(hRecoTPCKaonPhi);

  RecoByMC[2][1][0]->Add(hRecoTPCPionPt);
  RecoByMC[2][1][1]->Add(hRecoTPCPionEta);
  RecoByMC[2][1][2]->Add(hRecoTPCPionPhi);

  //// All Reco Bois TOF PID /////

  RecoByMC[3][0][0]->Add(hRecoTOFKaonPt);
  RecoByMC[3][0][1]->Add(hRecoTOFKaonEta);
  RecoByMC[3][0][2]->Add(hRecoTOFKaonPhi);

  RecoByMC[3][1][0]->Add(hRecoTOFPionPt);
  RecoByMC[3][1][1]->Add(hRecoTOFPionEta);
  RecoByMC[3][1][2]->Add(hRecoTOFPionPhi);

  //// All Reco Bois HFT PID /////

  RecoByMC[4][0][0]->Add(hRecoHFTKaonPt);
  RecoByMC[4][0][1]->Add(hRecoHFTKaonEta);
  RecoByMC[4][0][2]->Add(hRecoHFTKaonPhi);

  RecoByMC[4][1][0]->Add(hRecoHFTPionPt);
  RecoByMC[4][1][1]->Add(hRecoHFTPionEta);
  RecoByMC[4][1][2]->Add(hRecoHFTPionPhi);


  //// Divide by appropriate MC histo

  RecoByMC[0][0][0]->Divide(hMCD0KaonPt);
  RecoByMC[0][0][1]->Divide(hMCD0KaonEta);
  RecoByMC[0][0][2]->Divide(hMCD0KaonPhi);

  RecoByMC[0][1][0]->Divide(hMCD0PionPt);
  RecoByMC[0][1][1]->Divide(hMCD0PionEta);
  RecoByMC[0][1][2]->Divide(hMCD0PionPhi);

  //// All Reco Bois Geant Pid ////

  RecoByMC[1][0][0]->Divide(hMCKaonPt);
  RecoByMC[1][0][1]->Divide(hMCKaonEta);
  RecoByMC[1][0][2]->Divide(hMCKaonPhi);

  RecoByMC[1][1][0]->Divide(hMCPionPt);
  RecoByMC[1][1][1]->Divide(hMCPionEta);
  RecoByMC[1][1][2]->Divide(hMCPionPhi);

  //// All Reco Bois TPC PID /////

  RecoByMC[2][0][0]->Divide(hMCKaonPt);
  RecoByMC[2][0][1]->Divide(hMCKaonEta);
  RecoByMC[2][0][2]->Divide(hMCKaonPhi);

  RecoByMC[2][1][0]->Divide(hMCPionPt);
  RecoByMC[2][1][1]->Divide(hMCPionEta);
  RecoByMC[2][1][2]->Divide(hMCPionPhi);

  //// All Reco Bois TOF PID /////

  RecoByMC[3][0][0]->Divide(hMCKaonPt);
  RecoByMC[3][0][1]->Divide(hMCKaonEta);
  RecoByMC[3][0][2]->Divide(hMCKaonPhi);

  RecoByMC[3][1][0]->Divide(hMCPionPt);
  RecoByMC[3][1][1]->Divide(hMCPionEta);
  RecoByMC[3][1][2]->Divide(hMCPionPhi);

  //// All Reco Bois HFT PID /////

  RecoByMC[4][0][0]->Divide(hRecoTPCKaonPt);
  RecoByMC[4][0][1]->Divide(hRecoTPCKaonEta);
  RecoByMC[4][0][2]->Divide(hRecoTPCKaonPhi);

  RecoByMC[4][1][0]->Divide(hRecoTPCPionPt);
  RecoByMC[4][1][1]->Divide(hRecoTPCPionEta);
  RecoByMC[4][1][2]->Divide(hRecoTPCPionPhi);


  TCanvas *c[5];
  for (int i = 0; i < 5; i++){
    c[i] = new TCanvas(RatioType[i].Data(), RatioType[i].Data(), 800, 800);
    c[i]->Divide(3,2);
    c[i]->cd(1);
    RecoByMC[i][0][0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    RecoByMC[i][0][0]->GetYaxis()->SetRangeUser(0, 2);
    RecoByMC[i][0][0]->Draw();
    c[i]->cd(2);
    RecoByMC[i][0][1]->GetXaxis()->SetTitle("#eta");
    RecoByMC[i][0][1]->GetYaxis()->SetRangeUser(0, 2);
    RecoByMC[i][0][1]->Draw();
    c[i]->cd(3);
    RecoByMC[i][0][2]->GetXaxis()->SetTitle("#phi");
    RecoByMC[i][0][2]->GetYaxis()->SetRangeUser(0, 2);
    RecoByMC[i][0][2]->Draw();

    c[i]->cd(4);
    RecoByMC[i][1][0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    RecoByMC[i][1][0]->GetYaxis()->SetRangeUser(0, 2);
    RecoByMC[i][1][0]->Draw();
    c[i]->cd(5);
    RecoByMC[i][1][1]->GetXaxis()->SetTitle("#eta");
    RecoByMC[i][1][1]->GetYaxis()->SetRangeUser(0, 2);
    RecoByMC[i][1][1]->Draw();
    c[i]->cd(6);
    RecoByMC[i][1][2]->GetXaxis()->SetTitle("#phi");
    RecoByMC[i][1][2]->GetYaxis()->SetRangeUser(0, 2);
    RecoByMC[i][1][2]->Draw();

  }

  TCanvas *mass = new TCanvas("InvMass", "InvMass", 800, 800);

  TF1 *fitfunc = new TF1("fitfunc", "gaus", 1.7, 2.02);
  fitfunc->SetLineColor(6);
  fitfunc->SetLineWidth(5);

  fitfunc->SetParLimits(0, 0, 5000);
  fitfunc->SetParameter(1, 1.865);
  fitfunc->SetParLimits(1, 1.85, 1.88);
  fitfunc->SetParLimits(2, 0., 0.2);

  int status = hD0Mass->Fit("fitfunc", "RLQ");

  // if (fitfunc->GetChisquare()/fitfunc->GetNDF() > 1.5) status = 1;

  // cout << status << endl;

  // double fitloops = 0;

  // while (status != 0){
  //   if (fitloops >= 50) break;
  //       fitloops += 1;
  //       fitfunc->SetParameter(0, fitfunc->GetParameter(0));
  //       fitfunc->SetParLimits(0, 0, 5000);
  //       fitfunc->SetParameter(1, fitfunc->GetParameter(1));
  //       fitfunc->SetParameter(2, fitfunc->GetParameter(2));
  //       fitfunc->SetParLimits(1, 1.86, 1.87);
  //       fitfunc->SetParLimits(2, 0., 0.2);

  //       status = hD0Mass->Fit("fitfunc", "RLQ");

  //       if (fitfunc->GetChisquare()/fitfunc->GetNDF() > 1.5) status = 1;

  //       cout << status << endl;

  //   }

  hD0Mass->Draw();

  TPaveText *pt = new TPaveText(0.12,0.65,0.52,0.85, "NDC"); // NDC sets coords
                                              // relative to pad dimensions
  pt->SetTextSize(0.03); 
  pt->SetTextAlign(12);

  auto pt_text0 = pt->AddText(Form("D0 + #bar{D0} M = %1.3f GeV/c^{2}", fitfunc->GetParameter(1)));
  auto pt_text1 = pt->AddText(Form("#sigma = %1.4f GeV/c^{2}", fitfunc->GetParameter(2)));
  auto pt_text3 = pt->AddText(Form("#frac{#chi^{2}}{NDF} = %1.2f", fitfunc->GetChisquare()/fitfunc->GetNDF()));

  pt->Draw("SAME");




  // TCanvas *c = new TCanvas("c", "c", 800, 800);
  // c->Divide(3,1);
  // c->cd(1);
  // hHFTPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  // hHFTPt->GetYaxis()->SetRangeUser(0, 2);
  // hHFTPt->Draw();

  // c->cd(2);
  // hHFTEta->GetXaxis()->SetTitle("#eta");
  // hHFTEta->GetYaxis()->SetRangeUser(0, 2);
  // hHFTEta->Draw();

  // c->cd(3);
  // hHFTPhi->GetXaxis()->SetTitle("#phi");
  // hHFTPhi->GetYaxis()->SetRangeUser(0, 2);
  // hHFTPhi->Draw();

  
}