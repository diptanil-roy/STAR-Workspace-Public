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

#pragma link C++ class vector<int> +;

using namespace std;
// using namespace RooFit;

// #pragma link C++ class StJetTreeStruct+;

// #pragma link C++ class vector<float> +;
// #pragma link C++ class vector<vector<float> >+;
// #pragma link C++ class vector<int> +;
// #pragma link C++ class vector<vector<int> >+;
#endif

#include "../BinDef.h"

void SmearPlotter(int pthatbin = 0, int SUPERITERATION = 0, int iteration = 3){

  cout << "File Reader Method" << endl;

  int col[6] = {kViolet, kAzure, kTeal, kSpring, kOrange, kPink};
  int colors[200];

  for (int i = 0; i < 6; i++){
      for (int j = -9; j <= 10; j++){
          colors[i*20 + (j+9)] = col[i] + j; 
      }
  }

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TFile *f;
  if (pthatbin == 1)f = new TFile("pt_0_3.root");
  else if (pthatbin == 0)f = new TFile("pt_3_inf.root");

  cout << f->GetName() <<  endl;
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

  JetTree->SetBranchAddress("Centrality", &Centrality);
  JetTree->SetBranchAddress("MCPrimaryVertex", &MCPrimaryVertex);
  JetTree->SetBranchAddress("RecoPrimaryVertex", &RecoPrimaryVertex);
  JetTree->SetBranchAddress("MCJetPt", &MCJetPt);
  JetTree->SetBranchAddress("MCJetEta", &MCJetEta);
  JetTree->SetBranchAddress("MCJetPhi", &MCJetPhi);
  JetTree->SetBranchAddress("MCJetArea", &MCJetArea);
  JetTree->SetBranchAddress("MCJetE", &MCJetE);
  JetTree->SetBranchAddress("MCJetNConst", &MCJetNConst);
  JetTree->SetBranchAddress("MCD0Z", &MCD0Z);
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
  JetTree->SetBranchAddress("RecoD0Z", &RecoD0Z);
  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  cout << "Here" << endl;

  TH1D *PtDiff[3][extendedJetPtBins+1];

  TH1D *D0PtDiff[3][nBinsMCD0Pt+1];

  TH1D *D0PtDiffAllCent[nBinsMCD0Pt+1];

  TH2D *RecoPtvsMCPt[3];

  TH2D *RecoD0PtvsMCD0Pt[3];

  TH1D *DeltaPhi[4];

  TH1D *RecoDeltaPhi[4];

  TH1D *tmpMCPt = new TH1D("tmpMCPt", "tmpMCPt", extendedJetPtBins, extendedJetPt);

  TH1D *tmpMCD0Pt = new TH1D("tmpMCD0Pt", "tmpMCPt", nBinsMCD0Pt, MCD0PtBins);

  TH1D *MCPt = new TH1D("MCPt", "MCJetPt", extendedJetPtBins, extendedJetPt);
  TH1D *MCZ = new TH1D("MCZ", "MCJetZ", nBinsMCZ, ZBins);

  for (int cent = 0; cent < 3; cent++){
    RecoPtvsMCPt[cent] = new TH2D(Form("RecoPtvsMCPt_%i", cent), Form("RecoPtvsMCPt_%i", cent), extendedRecoJetPtBins, extendedRecoJetPt, extendedJetPtBins, extendedJetPt);

    RecoD0PtvsMCD0Pt[cent] = new TH2D(Form("RecoD0PtvsMCD0Pt_%i", cent), Form("RecoD0PtvsMCD0Pt_%i", cent), nBinsRecoD0Pt, RecoD0PtBins, nBinsMCD0Pt, MCD0PtBins);

    DeltaPhi[cent] = new TH1D(Form("DeltaPhi_%i", cent), Form("DeltaPhi_%i", cent), 1000, -0.1, 0.5);

    RecoDeltaPhi[cent] = new TH1D(Form("RecoDeltaPhi_%i", cent), Form("RecoDeltaPhi_%i", cent), 1000, -0.1, 0.5);

    for (int i = 0; i < extendedJetPtBins; i++){
      PtDiff[cent][i] = new TH1D(Form("PtDiff_%i_%i", cent, i), Form("PtDiff_%i_%i", cent, i), 160, -40, 40);
    }

    for (int i = 0; i < nBinsMCD0Pt; i++){
      D0PtDiff[cent][i] = new TH1D(Form("D0PtDiff_%i_%i", cent, i), Form("D0PtDiff_%i_%i", cent, i), 500, -5, 5);
    }

    PtDiff[cent][extendedJetPtBins] = new TH1D(Form("PtDiff_%i_%i", cent, extendedJetPtBins), Form("PtDiff_%i_%i", cent, extendedJetPtBins), 160, -40, 40);

    D0PtDiff[cent][nBinsMCD0Pt] = new TH1D(Form("D0PtDiff_%i_%i", cent, nBinsMCD0Pt), Form("D0PtDiff_%i_%i", cent, nBinsMCD0Pt), 500, -5, 5);
  }

  DeltaPhi[3] = new TH1D(Form("DeltaPhi_%i", 3), Form("DeltaPhi_%i", 3), 1000, -0.1, 0.5);

  RecoDeltaPhi[3] = new TH1D(Form("RecoDeltaPhi_%i", 3), Form("RecoDeltaPhi_%i", 3), 1000, -0.1, 0.5);

  for (int i = 0; i < nBinsMCD0Pt; i++){
    D0PtDiffAllCent[i] = new TH1D(Form("D0PtDiffAllCent_%i", i), Form("D0PtDiffAllCent_%i", i), 500, -5, 5);
  }

  D0PtDiffAllCent[nBinsMCD0Pt] = new TH1D(Form("D0PtDiffAllCent_%i", nBinsMCD0Pt), Form("D0PtDiffAllCent_%i", nBinsMCD0Pt), 500, -5, 5);

  int nentries = JetTree->GetEntries();

  for (int i = 0; i < nentries; i++){

    if (i%1000000 == 0) cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

    if (centhistogramtofill < 0) continue;

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 5 && RecoD0Pt < 10;
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 30;
    bool isRecoJetPt = RecoJetCorrPt > 0 && RecoJetCorrPt < 30; // 78% of jets with pT > 5 GeV is captured within this range
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;

    int D0HistBin = tmpMCD0Pt->FindBin(MCD0Pt) - 1;
    D0PtDiff[centhistogramtofill][D0HistBin]->Fill(RecoD0Pt - MCD0Pt);
    D0PtDiff[centhistogramtofill][nBinsMCD0Pt]->Fill(RecoD0Pt - MCD0Pt);

    D0PtDiffAllCent[D0HistBin]->Fill(RecoD0Pt - MCD0Pt);
    D0PtDiffAllCent[nBinsMCD0Pt]->Fill(RecoD0Pt - MCD0Pt);

    RecoD0PtvsMCD0Pt[centhistogramtofill]->Fill(RecoD0Pt, MCD0Pt);

    if (!isMCD0Pt || !isMCJetPt) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);

    MCPt->Fill(MCJetPt);
    MCZ->Fill(mcz);

    DeltaPhi[centhistogramtofill]->Fill(dPhi(MCJetPhi, MCD0Phi));
    DeltaPhi[3]->Fill(dPhi(MCJetPhi, MCD0Phi));

    RecoDeltaPhi[centhistogramtofill]->Fill(dPhi(RecoJetPhi, RecoD0Phi));
    RecoDeltaPhi[3]->Fill(dPhi(RecoJetPhi, RecoD0Phi));

    int histbin = tmpMCPt->FindBin(MCJetPt) - 1;
    // cout << MCJetPt << "\t" << histbin << endl;
    PtDiff[centhistogramtofill][histbin]->Fill(RecoJetCorrPt - MCJetPt);
    PtDiff[centhistogramtofill][extendedJetPtBins]->Fill(RecoJetCorrPt - MCJetPt);
    RecoPtvsMCPt[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt);
  }

  cout << "Number of entries = " << nentries << endl;

  TCanvas *c[3];

  for (int cent = 0; cent < 3; cent++){
    c[cent] = new TCanvas(Form("c_%i", cent), Form("c_%i", cent), 1000, 1000);
    c[cent]->cd();
    gPad->SetLogy();
    for (int i = 0; i < extendedJetPtBins; i++){
      SetName(PtDiff[cent][i], Form("Cent %i, Pt in [%.1f, %.1f]", cent, extendedJetPt[i], extendedJetPt[i+1]));
      SetColor(PtDiff[cent][i], colors[i*8]);
      cout << "Cent PtBin Integral Mean RMS = " << cent << "\t" << i << "\t" << PtDiff[cent][i]->Integral() << "\t" << PtDiff[cent][i]->GetMean() << "\t" << PtDiff[cent][i]->GetRMS() << endl;
      PtDiff[cent][i]->Scale(1./PtDiff[cent][i]->Integral());
      PtDiff[cent][i]->Draw("SAME");
    }
    gPad->BuildLegend();
  }

  TCanvas *d = new TCanvas("d", "d", 1000, 1000);
  for (int cent = 0; cent < 3; cent++){
    d->cd();
    cout << "Cent Integral Mean RMS = " << cent << "\t" << PtDiff[cent][extendedJetPtBins]->Integral() << "\t" << PtDiff[cent][extendedJetPtBins]->GetMean() << "\t" << PtDiff[cent][extendedJetPtBins]->GetRMS() << endl;
    SetName(PtDiff[cent][extendedJetPtBins], Form("Cent %i, Pt in [%.1f, %.1f]", cent, 5., 30.));
    PtDiff[cent][extendedJetPtBins]->Scale(1./PtDiff[cent][extendedJetPtBins]->Integral());
    SetColor(PtDiff[cent][extendedJetPtBins], col[cent*2]);
    PtDiff[cent][extendedJetPtBins]->Draw("SAME");
    gPad->SetLogy();
  }

  TCanvas *e[3];

  for (int cent = 0; cent < 3; cent++){
    e[cent] = new TCanvas(Form("e_%i", cent), Form("e_%i", cent), 1000, 1000);
    e[cent]->cd();
    gPad->SetLogz();
    RecoPtvsMCPt[cent]->Draw("COLZ");
  }

  TCanvas *fg[3];

  for (int cent = 0; cent < 3; cent++){
    fg[cent] = new TCanvas(Form("f_%i", cent), Form("f_%i", cent), 1000, 1000);
    fg[cent]->cd();
    gPad->SetLogy();
    for (int i = 0; i < nBinsMCD0Pt; i++){
      SetName(D0PtDiff[cent][i], Form("D0 Cent %i, Pt in [%.1f, %.1f]", cent, MCD0PtBins[i], MCD0PtBins[i+1]));
      SetColor(D0PtDiff[cent][i], colors[i*8]);
      cout << "Cent PtBin Integral Mean RMS = " << cent << "\t" << i << "\t" << D0PtDiff[cent][i]->Integral() << "\t" << D0PtDiff[cent][i]->GetMean() << "\t" << D0PtDiff[cent][i]->GetRMS() << endl;
      D0PtDiff[cent][i]->Scale(1./D0PtDiff[cent][i]->Integral());
      D0PtDiff[cent][i]->Draw("SAME");
    }
    gPad->BuildLegend();
  }

  TCanvas *g = new TCanvas("g", "g", 1000, 1000);
  for (int cent = 0; cent < 3; cent++){
    g->cd();
    SetName(D0PtDiff[cent][nBinsMCD0Pt], Form("D0 Cent %i, Pt in [%.1f, %.1f]", cent, 1., 10.));
    cout << "Cent Integral Mean RMS = " << cent << "\t" << D0PtDiff[cent][nBinsMCD0Pt]->Integral() << "\t" << D0PtDiff[cent][nBinsMCD0Pt]->GetMean() << "\t" << D0PtDiff[cent][nBinsMCD0Pt]->GetRMS() << endl;
    D0PtDiff[cent][nBinsMCD0Pt]->Scale(1./D0PtDiff[cent][nBinsMCD0Pt]->Integral());
    SetColor(D0PtDiff[cent][nBinsMCD0Pt], col[cent*2]);
    D0PtDiff[cent][nBinsMCD0Pt]->Draw("SAME");
    gPad->SetLogy();
  }

  TCanvas *h[3];

  for (int cent = 0; cent < 3; cent++){
    h[cent] = new TCanvas(Form("h_%i", cent), Form("h_%i", cent), 1000, 1000);
    h[cent]->cd();
    gPad->SetLogz();
    RecoD0PtvsMCD0Pt[cent]->Draw("COLZ");
  }

  TCanvas *ha = new TCanvas("ha", "ha", 1000, 1000);
  ha->cd();
  gPad->SetLogy();
  for (int i = 0; i <= nBinsMCD0Pt; i++){
    cout << "Cent Integral Mean RMS = " << i << "\t" << D0PtDiffAllCent[i]->Integral() << "\t" << D0PtDiffAllCent[i]->GetMean() << "\t" << D0PtDiffAllCent[i]->GetRMS() << endl;
    if (i!=nBinsMCD0Pt) SetName(D0PtDiffAllCent[i], Form("D0 Pt in [%.1f, %.1f]", MCD0PtBins[i], MCD0PtBins[i+1]));
    else SetName(D0PtDiffAllCent[i], Form("D0 Pt in [%.1f, %.1f]", 1., 10.));
    SetColor(D0PtDiffAllCent[i], colors[i*8]);
    D0PtDiffAllCent[i]->Scale(1./D0PtDiffAllCent[i]->Integral());
    D0PtDiffAllCent[i]->Draw("SAME");
  }
  gPad->BuildLegend();

  TCanvas *hb = new TCanvas("hb", "hb", 1000, 1000);
  hb->cd();
  gPad->SetLogy();
  for (int i = 0; i < 4; i++){
    DeltaPhi[i]->Scale(1./DeltaPhi[i]->Integral());
    SetColor(DeltaPhi[i], col[i]);
    DeltaPhi[i]->Draw("SAME");
  }
  gPad->BuildLegend();

  TCanvas *hc = new TCanvas("hc", "hc", 1000, 1000);
  hc->cd();
  gPad->SetLogy();
  for (int i = 0; i < 4; i++){
    RecoDeltaPhi[i]->Scale(1./RecoDeltaPhi[i]->Integral());
    SetColor(RecoDeltaPhi[i], col[i]);
    RecoDeltaPhi[i]->Draw("SAME");
  }
  gPad->BuildLegend();

  TFile *out = new TFile("SmearFactors.root", "RECREATE");
  out->cd();
  for (int cent = 0; cent < 3; cent++){
    for (int i = 0; i <= extendedJetPtBins; i++){
      PtDiff[cent][i]->Write();
    }
  }

  for (int cent = 0; cent < 3; cent++){
    RecoPtvsMCPt[cent]->Write();
  }

  for (int cent = 0; cent < 3; cent++){
    for (int i = 0; i <= nBinsMCD0Pt; i++){
      D0PtDiff[cent][i]->Write();
    }
  }

  for (int i = 0; i <= nBinsMCD0Pt; i++){
    D0PtDiffAllCent[i]->Write();
  }

  for (int cent = 0; cent < 3; cent++){
    RecoD0PtvsMCD0Pt[cent]->Write();
  }

  for (int cent = 0; cent < 4; cent++){
    DeltaPhi[cent]->Write();
  }

  for (int cent = 0; cent < 4; cent++){
    RecoDeltaPhi[cent]->Write();
  }

  MCPt->Write();
  MCZ->Write();

  out->Close();
}
