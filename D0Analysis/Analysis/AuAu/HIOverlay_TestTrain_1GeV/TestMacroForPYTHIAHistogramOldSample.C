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
#include "THnSparse.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TRandom3.h"
// #include "StJetTreeStruct.h"
#include <vector>

using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

void TestMacroForPYTHIAHistogramOldSample(double lowptcutoff = 5.0, double highptcutoff = 20.0){
	TFile *f = new TFile("/Volumes/WorkDrive/work/2022/Response2022/Response_Mar24.root");
	f->cd("SimJetSaverSmear_Central");

	cout << gDirectory->GetName() << endl;

  TTree *JetTree = (TTree *)gDirectory->Get("MCJets");
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
  Float_t         MCD0Z;
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
  Float_t         RecoD0Z;
  Float_t         RecoD0Pt;
  Float_t         RecoD0Eta;
  Float_t         RecoD0Phi;
  Float_t         RecoPionPt;
  Float_t         RecoPionEta;
  Float_t         RecoPionPhi;
  Float_t         RecoKaonPt;
  Float_t         RecoKaonEta;
  Float_t         RecoKaonPhi;

  // JetTree->SetBranchAddress("Centrality", &Centrality);
  // JetTree->SetBranchAddress("MCPrimaryVertex", &MCPrimaryVertex);
  JetTree->SetBranchAddress("JetPt", &MCJetPt);
  JetTree->SetBranchAddress("JetEta", &MCJetEta);
  JetTree->SetBranchAddress("JetPhi", &MCJetPhi);
  JetTree->SetBranchAddress("JetArea", &MCJetArea);
  JetTree->SetBranchAddress("JetE", &MCJetE);
  JetTree->SetBranchAddress("JetNConst", &MCJetNConst);
  // JetTree->SetBranchAddress("MCD0Z", &MCD0Z);
  JetTree->SetBranchAddress("PionPt", &MCPionPt);
  JetTree->SetBranchAddress("PionEta", &MCPionEta);
  JetTree->SetBranchAddress("PionPhi", &MCPionPhi);
  JetTree->SetBranchAddress("KaonPt", &MCKaonPt);
  JetTree->SetBranchAddress("KaonEta", &MCKaonEta);
  JetTree->SetBranchAddress("KaonPhi", &MCKaonPhi);
 
  // JetTree->SetBranchAddress("RecoJetPt", &RecoJetPt);
  // JetTree->SetBranchAddress("RecoJetCorrPt", &RecoJetCorrPt);
  // JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
  // JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
  // JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
  // JetTree->SetBranchAddress("RecoJetE", &RecoJetE);
  // // JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
  // JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConst);
  // // JetTree->SetBranchAddress("RecoD0Z", &RecoD0Z);
  // JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  // JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  // JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  // JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  // JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  // JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  // JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  // JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  // JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  int nEntries = JetTree->GetEntries();

  const int nptbins = 6;
  double ptbins[nptbins + 1] = {5, 7, 9, 11, 13, 15, 20};

  TH1D *PYTHIAPt = new TH1D("PYTHIA pT", "PYTHIA pT", njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
  TH1D *PYTHIAZ = new TH1D("PYTHIA Z", "PYTHIA Z", nz_gen_bins, z_gen_low, z_gen_high);

  TH2D *PYTHIA = new TH2D("PYTHIA", "PYTHIA", njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

  TH1D *FONLLPt = new TH1D("FONLL pT", "FONLL pT", njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
  TH1D *FONLLZ = new TH1D("FONLL Z", "FONLL Z", nz_gen_bins, z_gen_low, z_gen_high);

  // TH1D *FONLLCurve = new TH1D("FONLLCurve", "FONLLCurve", nptbins, ptbins);

  TH2D *FONLLWeighted = new TH2D("FONLL Weighted", "FONLL Weighted", njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

  for (int i = 0; i < nEntries; i++){
  	if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    // cout << "Read Entry " << i << endl;
    JetTree->GetEntry(i);

    double MCPionPx = MCPionPt*TMath::Cos(MCPionPhi);
    double MCPionPy = MCPionPt*TMath::Sin(MCPionPhi);
    double MCPionPz = MCPionPt*TMath::SinH(MCPionEta);

    double MCKaonPx = MCKaonPt*TMath::Cos(MCKaonPhi);
    double MCKaonPy = MCKaonPt*TMath::Sin(MCKaonPhi);
    double MCKaonPz = MCKaonPt*TMath::SinH(MCKaonEta);

    TVector3 Pion(MCPionPx, MCPionPy, MCPionPz);
    TVector3 Kaon(MCKaonPx, MCKaonPy, MCKaonPz);

    TVector3 D0;
    D0 = Pion + Kaon;

    MCD0Pt = D0.Pt();
    MCD0Phi = D0.Phi();
    MCD0Eta = D0.PseudoRapidity();

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isMCD0Eta = abs(MCD0Eta) < 1.;
    bool isMCPionEta = abs(MCPionEta) < 1.;
    bool isMCKaonEta = abs(MCKaonEta) < 1.;
    bool isPionPt = MCPionPt > 0.6;
    bool isKaonPt = MCKaonPt > 0.6;

    if (!isMCD0Pt) continue;
    if (!isMCD0Eta) continue;
    // if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (!isMCJetEta) continue;
    // if (!isMCPionEta) continue;
    // if (!isMCKaonEta) continue;
    // if (!isPionPt) continue;
    // if (!isKaonPt) continue;
    // if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

    PYTHIAPt->Fill(MCJetPt);
    PYTHIAZ->Fill(mcz);

    PYTHIA->Fill(MCJetPt, mcz);
  }

  PYTHIAPt->Scale(1./PYTHIAPt->Integral());
  PYTHIAZ->Scale(1./PYTHIAZ->Integral());
  PYTHIA->Scale(1./PYTHIA->Integral());

  TString FileName = "FONLL.txt";

  // double pt, ptweight;
  // ifstream myfile (FileName.Data());
  // if (myfile.is_open())
  // {
  //   while ( !myfile.eof() )
  //   {
  //     myfile >> pt >> ptweight;
  //     // cout << x << "\t" << y << endl;
  //     int bin = FONLLPt->FindBin(pt+0.0001);
  //     FONLLCurve->SetBinContent(bin, ptweight);
  //   }
  // }

  // TH1D *Weight = (TH1D *)FONLLCurve->Clone("FONLLvPYTHIAWeights");
  // Weight->Divide(PYTHIAPt);

  TFile *FONLL = new TFile(Form("FONLL_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "READ");
  TH1D  *FONLLCurve = (TH1D *)FONLL->Get("FONLL");

  TH1D *Weight = (TH1D *)FONLLCurve->Clone("FONLLvPYTHIAWeights");
  Weight->Divide(PYTHIAPt);

  for (int i = 0; i < nEntries; i++){
  	if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    // cout << "Read Entry " << i << endl;
    JetTree->GetEntry(i);

    double MCPionPx = MCPionPt*TMath::Cos(MCPionPhi);
    double MCPionPy = MCPionPt*TMath::Sin(MCPionPhi);
    double MCPionPz = MCPionPt*TMath::SinH(MCPionEta);

    double MCKaonPx = MCKaonPt*TMath::Cos(MCKaonPhi);
    double MCKaonPy = MCKaonPt*TMath::Sin(MCKaonPhi);
    double MCKaonPz = MCKaonPt*TMath::SinH(MCKaonEta);

    TVector3 Pion(MCPionPx, MCPionPy, MCPionPz);
    TVector3 Kaon(MCKaonPx, MCKaonPy, MCKaonPz);

    TVector3 D0;
    D0 = Pion + Kaon;

    MCD0Pt = D0.Pt();
    MCD0Phi = D0.Phi();
    MCD0Eta = D0.PseudoRapidity();

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isMCD0Eta = abs(MCD0Eta) < 1.;
    bool isMCPionEta = abs(MCPionEta) < 1.;
    bool isMCKaonEta = abs(MCKaonEta) < 1.;
    bool isPionPt = MCPionPt > 0.6;
    bool isKaonPt = MCKaonPt > 0.6;

    if (!isMCD0Pt) continue;
    if (!isMCD0Eta) continue;
    // if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (!isMCJetEta) continue;
    // if (!isMCPionEta) continue;
    // if (!isMCKaonEta) continue;
    // if (!isPionPt) continue;
    // if (!isKaonPt) continue;
    // if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

    double weight = Weight->GetBinContent(Weight->FindBin(MCJetPt));

    FONLLPt->Fill(MCJetPt, weight);
    FONLLZ->Fill(mcz, weight);

    FONLLWeighted->Fill(MCJetPt, mcz, weight);
  }

  FONLLPt->Scale(1./FONLLPt->Integral());
  FONLLZ->Scale(1./FONLLZ->Integral());
  FONLLWeighted->Scale(1./FONLLWeighted->Integral());

  TH2F *Weight2D = (TH2F *)PYTHIA->Clone("Weight2D");
  Weight2D->Divide(FONLLWeighted);

  TCanvas *d = new TCanvas("d", "d", 1200, 600);
  Weight2D->Draw("COLZ");

  TH1D *PYTHIAWPt = new TH1D("PYTHIA W pT", "PYTHIA W pT", njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
  TH1D *PYTHIAWZ = new TH1D("PYTHIA W Z", "PYTHIA W Z", nz_gen_bins, z_gen_low, z_gen_high);
  TH2D *PYTHIAZWeighted = new TH2D("PYTHIA Z Weighted", "PYTHIA Z Weighted", njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

  for (int i = 0; i < nEntries; i++){
    if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    // cout << "Read Entry " << i << endl;
    JetTree->GetEntry(i);

    double MCPionPx = MCPionPt*TMath::Cos(MCPionPhi);
    double MCPionPy = MCPionPt*TMath::Sin(MCPionPhi);
    double MCPionPz = MCPionPt*TMath::SinH(MCPionEta);

    double MCKaonPx = MCKaonPt*TMath::Cos(MCKaonPhi);
    double MCKaonPy = MCKaonPt*TMath::Sin(MCKaonPhi);
    double MCKaonPz = MCKaonPt*TMath::SinH(MCKaonEta);

    TVector3 Pion(MCPionPx, MCPionPy, MCPionPz);
    TVector3 Kaon(MCKaonPx, MCKaonPy, MCKaonPz);

    TVector3 D0;
    D0 = Pion + Kaon;

    MCD0Pt = D0.Pt();
    MCD0Phi = D0.Phi();
    MCD0Eta = D0.PseudoRapidity();

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isMCD0Eta = abs(MCD0Eta) < 1.;
    bool isMCPionEta = abs(MCPionEta) < 1.;
    bool isMCKaonEta = abs(MCKaonEta) < 1.;
    bool isPionPt = MCPionPt > 0.6;
    bool isKaonPt = MCKaonPt > 0.6;

    if (!isMCD0Pt) continue;
    if (!isMCD0Eta) continue;
    // if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (!isMCJetEta) continue;
    // if (!isMCPionEta) continue;
    // if (!isMCKaonEta) continue;
    // if (!isPionPt) continue;
    // if (!isKaonPt) continue;
    // if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

    int bin = Weight2D->GetBin(MCJetPt, mcz);
    double w = Weight2D->GetBinContent(bin);

    PYTHIAWPt->Fill(MCJetPt, w);
    PYTHIAWZ->Fill(mcz, w);
    PYTHIAZWeighted->Fill(MCJetPt, mcz, w);
  }

  TCanvas *c = new TCanvas("c", "c", 1200, 600);
  c->Divide(5);
  c->cd(1);
  gPad->SetLogy();
  // PYTHIAPt->Draw();
  // FONLLPt->Draw("SAME");
  // SetColor(FONLLPt, kRed);
  PYTHIAWPt->Draw("EP SAME");
  SetColor(PYTHIAWPt, kGreen-2);
  c->cd(2);
  gPad->SetLogy();
  // PYTHIAZ->Draw();
  // FONLLZ->Draw("SAME");
  // SetColor(FONLLZ, kRed);
  PYTHIAWZ->Draw("EP SAME");
  SetColor(PYTHIAWZ, kGreen-2);
  c->cd(3);
  gPad->SetLogz();
  PYTHIA->Draw("COLZ");
  PYTHIA->GetZaxis()->SetRangeUser(pow(10, -6),1);
  c->cd(4);
  gPad->SetLogz();
  FONLLWeighted->Draw("COLZ");
  FONLLWeighted->GetZaxis()->SetRangeUser(pow(10, -6),1);
  c->cd(5);
  gPad->SetLogy();
  Weight->Draw();

  TFile *g = new TFile(Form("PYTHIA_Test_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "RECREATE");
  g->cd();
  PYTHIAPt->Write();
  PYTHIAZ->Write();
  PYTHIA->Write();
  Weight->Write();
  g->Close();
}