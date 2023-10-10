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
// #ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
// #include <cstdlib>
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
// #include "StJetTreeStruct.h"
#include <vector>

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void NonPtInefficiency(TString mode = "QM"){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);

  //// Importing all the files below this line

  // TFile *f = new TFile("HI5GeV_WithMaxTrackTower.root");
  TFile *f = new TFile("HIResponse_MCIncluded_Feb6.root");

  cout << f->GetName() <<  endl;
  f->cd("HIJetSaver");

  TTree *JetTree = (TTree *)gDirectory->Get("Jets");
  // TTree *RecoJetTree = (TTree *)gDirectory->Get("RecoJets");

  Float_t         Centrality;
  Float_t         CWeight;
  vector<double>  *MCPrimaryVertex = new vector<double>;
  vector<double>  *RecoPrimaryVertex = new vector<double>;
  Float_t         RecoMaxTrackPt;
  Float_t         RecoMaxTowerEtBeforeHC;
  Float_t         RecoMaxTowerEtAfterHC;
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
  JetTree->SetBranchAddress("Weight", &CWeight);
  JetTree->SetBranchAddress("MCPrimaryVertex", &MCPrimaryVertex);
  JetTree->SetBranchAddress("RecoPrimaryVertex", &RecoPrimaryVertex);
  JetTree->SetBranchAddress("RecoMaxTrackPt", &RecoMaxTrackPt);
  JetTree->SetBranchAddress("RecoMaxTowerEtBeforeHC", &RecoMaxTowerEtBeforeHC);
  JetTree->SetBranchAddress("RecoMaxTowerEtAfterHC", &RecoMaxTowerEtAfterHC);
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
  
  if (mode == "HI"){  
    JetTree->SetBranchAddress("RecoJetPt", &RecoJetPt);
    JetTree->SetBranchAddress("RecoJetCorrPt", &RecoJetCorrPt);
    JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
    JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
    JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
    JetTree->SetBranchAddress("RecoJetE", &RecoJetE);
    JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
    JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConst);
  }
  else if (mode == "QM"){
    JetTree->SetBranchAddress("RecoJetPtFromPYTHIA", &RecoJetPt);
    JetTree->SetBranchAddress("RecoJetEtaFromPYTHIA", &RecoJetEta);
    JetTree->SetBranchAddress("RecoJetPhiFromPYTHIA", &RecoJetPhi);
    JetTree->SetBranchAddress("RecoJetAreaFromPYTHIA", &RecoJetArea);
    JetTree->SetBranchAddress("RecoJetEFromPYTHIA", &RecoJetE);
    JetTree->SetBranchAddress("RecoJetNConstFromPYTHIA", &RecoJetNConst);
  }
  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  TFile* FONLL = new TFile("FONLL_Pt_5_20.root","READ");
  TH1F * FONLLCurve = (TH1F*) FONLL->Get("FONLL");

  TFile *PYTHIA = new TFile("PYTHIA_Pt_5_20.root", "READ");
  TH1F  *PYTHIAPtCurve = (TH1F *) PYTHIA->Get("PYTHIA pT");
  TH1F  *PYTHIAZCurve  = (TH1F *) PYTHIA->Get("PYTHIA Z");
  TH2F  *PYTHIA2D      = (TH2F *) PYTHIA->Get("PYTHIA");
  TH1F  *FONLLvPYTHIAWeights = (TH1F *) PYTHIA->Get("FONLLWeights");

  TFile* CentWeights = new TFile("Cent_Weight.root");
  TH1F * CENTWeight = (TH1F *)CentWeights->Get("Centrality Weights");

  TFile *pTdepdpTSmearing = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/Pt_Res_FINAL_Jan26_2023.root");
  TFile *EtaResFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/Eta_Res_FINAL_Jan26_2023.root");
  TFile *PhiResFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/Phi_Res_FINAL_Jan26_2023.root");

  TH1F *dPt[3][40]; // This is the dPt vs Pt generated from Embed_April2.root
  TH1F *Eta_Mean[3];
  TH1F *Eta_Sigma[3];
  TH1F *Phi_Mean[3];
  TH1F *Phi_Sigma[3];

  for (int i = 0; i < 40; i++){
    dPt[0][i] = (TH1F *)pTdepdpTSmearing->Get(Form("dPt_Cent_%i_%i_Pt_%i_%i", 0, 10, i, i+1));
  }
  Eta_Mean[0] = (TH1F *)EtaResFile->Get(Form("eta_mean_pt_%i_%i", 0, 10));
  Eta_Sigma[0] = (TH1F *)EtaResFile->Get(Form("eta_sigma_pt_%i_%i", 0, 10));

  Phi_Mean[0] = (TH1F *)PhiResFile->Get(Form("phi_mean_pt_%i_%i", 0, 10));
  Phi_Sigma[0] = (TH1F *)PhiResFile->Get(Form("phi_sigma_pt_%i_%i", 0, 10));
  
  for (int i = 0; i < 40; i++){
    dPt[1][i] = (TH1F *)pTdepdpTSmearing->Get(Form("dPt_Cent_%i_%i_Pt_%i_%i", 10, 40, i, i+1));
  }
  Eta_Mean[1] = (TH1F *)EtaResFile->Get(Form("eta_mean_pt_%i_%i", 10, 40));
  Eta_Sigma[1] = (TH1F *)EtaResFile->Get(Form("eta_sigma_pt_%i_%i", 10, 40));

  Phi_Mean[1] = (TH1F *)PhiResFile->Get(Form("phi_mean_pt_%i_%i", 10, 40));
  Phi_Sigma[1] = (TH1F *)PhiResFile->Get(Form("phi_sigma_pt_%i_%i", 10, 40));
  
  for (int i = 0; i < 40; i++){
    dPt[2][i] = (TH1F *)pTdepdpTSmearing->Get(Form("dPt_Cent_%i_%i_Pt_%i_%i", 40, 80, i, i+1));
  }

  Eta_Mean[2] = (TH1F *)EtaResFile->Get(Form("eta_mean_pt_%i_%i", 40, 80));
  Eta_Sigma[2] = (TH1F *)EtaResFile->Get(Form("eta_sigma_pt_%i_%i", 40, 80));

  Phi_Mean[2] = (TH1F *)PhiResFile->Get(Form("phi_mean_pt_%i_%i", 40, 80));
  Phi_Sigma[2] = (TH1F *)PhiResFile->Get(Form("phi_sigma_pt_%i_%i", 40, 80));

  //// Imported all the files above this line

  TH1D *hPtMissedDueToEta[3];
  TH1D *hPtAll[3];
  TH1D *hPtRatio[3];

  TH1D *hEtaMissedDueToEta[3];
  TH1D *hEtaAll[3];
  TH1D *hEtaRatio[3];

  TH2D *hPtEtaMissedDueToEta[3];
  TH2D *hPtEtaAll[3];
  TH2D *hPtEtaRatio[3];

  TH1D *hPtMissedDueToDeltaR[3];
  TH1D *hPtDueToDeltaRRatio[3];

  const int nbins_jpt = 6;
  double binning_jpt[nbins_jpt+1] = {5,7,9,11,13,15,20};

  for (int i = 0; i < 3; i++){
    hPtMissedDueToEta[i] = new TH1D(Form("hPtMissedDueToEta_%i", i), Form("hPtMissedDueToEta_%i", i), nbins_jpt, binning_jpt);
    hPtAll[i] = new TH1D(Form("hPtAll_%i", i), Form("hPtAll_%i", i), nbins_jpt, binning_jpt);

    hPtMissedDueToDeltaR[i] = new TH1D(Form("hPtMissedDueToDeltaR_%i", i), Form("hPtMissedDueToDeltaR_%i", i), nbins_jpt, binning_jpt);

    hEtaMissedDueToEta[i] = new TH1D(Form("hEtaMissedDueToEta_%i", i), Form("hEtaMissedDueToEta_%i", i), 20, -1, 1);
    hEtaAll[i] = new TH1D(Form("hEtaAll_%i", i), Form("hEtaAll_%i", i), 20, -1, 1);

    hPtEtaMissedDueToEta[i] = new TH2D(Form("hPtEtaMissedDueToEta_%i", i), Form("hPtEtaMissedDueToEta_%i", i), nbins_jpt, binning_jpt, 20, -1, 1);
    hPtEtaAll[i] = new TH2D(Form("hPtEtaAll_%i", i), Form("hPtEtaAll_%i", i), nbins_jpt, binning_jpt, 20, -1, 1);
  }

  TH1D *tmppT = new TH1D("tmppT", "tmppT", 40, 0, 40);

  int nentries = JetTree->GetEntries();

  int lowlimit = 0;
  int highlimit = nentries;

  TRandom* r = new TRandom(0);

  for (int i = lowlimit; i < highlimit; i++){

    JetTree->GetEntry(i);

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

    if (centhistogramtofill < 0) continue;

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 5 && RecoD0Pt < 10;
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    // if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);
    double MCDeltaR = dR(dEta(MCJetEta, MCD0Eta), dPhi(MCJetPhi, MCD0Phi));

    if (mcz == 1.0) mcz = 0.999; // Padding the boundaries

    if (mode == "QM"){
      int ptbin = tmppT->FindBin(RecoJetPt);
      RecoJetCorrPt = RecoJetPt + dPt[centhistogramtofill][ptbin-1]->GetRandom();
      double eta_mean_val = Eta_Mean[centhistogramtofill]->GetBinContent(ptbin);
      double eta_sigma_val = Eta_Sigma[centhistogramtofill]->GetBinContent(ptbin);
      double phi_mean_val = Phi_Mean[centhistogramtofill]->GetBinContent(ptbin);
      double phi_sigma_val = Phi_Sigma[centhistogramtofill]->GetBinContent(ptbin);
      RecoJetEta    = RecoJetEta + r->Gaus(eta_mean_val, eta_sigma_val);
      RecoJetPhi    = RecoJetPhi + r->Gaus(phi_mean_val, phi_sigma_val);
    }

    if(RecoJetPhi>=2.*TMath::Pi()) RecoJetPhi = RecoJetPhi - 2.*TMath::Pi();
    if(RecoJetPhi<0) RecoJetPhi = 2.*TMath::Pi() + RecoJetPhi;

    double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);
    double RecoDeltaR = dR(dEta(RecoJetEta, RecoD0Eta), dPhi(RecoJetPhi, RecoD0Phi));

    bool isMCJetEta = abs(MCJetEta) < 0.601;
    bool isMCDeltaR = MCDeltaR < 0.4;

    bool isRecoJetPt = RecoJetCorrPt > 3 && RecoJetCorrPt < 30;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoDeltaR = RecoDeltaR < 0.4;

    bool isOK   = isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta;
    bool isMISS = !isOK && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOK && isRecoJetPt && isRecoJetEta;

    if (!isOK && !isMISS && !isFAKE) {cout << "Doesn't fit any criteria" << endl; continue;}

    int centbin = CENTWeight->FindBin(Centrality);

    double wcent = CENTWeight->GetBinContent(centbin);
    double fonllweight = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));

    // double w = wcent*fonllweight;
    double w = wcent;
    if (mode == "QM") w = 1;
    // double w = 1.;
    // if (!isRecoJetPt) continue;

    if (!isRecoJetEta){
      hPtMissedDueToEta[centhistogramtofill]->Fill(MCJetPt, w);
      hEtaMissedDueToEta[centhistogramtofill]->Fill(MCJetEta, w);
      hPtEtaMissedDueToEta[centhistogramtofill]->Fill(MCJetPt, MCJetEta, w);
    }

    if (RecoJetNConst==0){
      hPtMissedDueToDeltaR[centhistogramtofill]->Fill(MCJetPt, w);
      // cout << Form("%.2f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\n", MCD0Pt, MCJetPt, MCJetEta, RecoD0Pt, RecoJetCorrPt, RecoJetEta);
    }

    hPtAll[centhistogramtofill]->Fill(MCJetPt, w);
    hEtaAll[centhistogramtofill]->Fill(MCJetEta, w);
    hPtEtaAll[centhistogramtofill]->Fill(MCJetPt, MCJetEta, w);
  }

  TCanvas *c = new TCanvas("c", "c", 2000, 1000);
  c->Divide(3, 2);
  for (int cent = 0; cent < 3; cent++){
    c->cd(1 + cent);
    hPtEtaMissedDueToEta[cent]->Draw("COLZ");
    gPad->SetLogz();

    c->cd(4 + cent);
    hPtEtaAll[cent]->Draw("COLZ");
    gPad->SetLogz(); 
  }

  int color[3] = {kGreen-2, kBlue, kBlack};
  int marker[3] = {20, 21, 20};

  for (int cent = 0; cent < 3; cent++){
    hPtRatio[cent] = (TH1D *)hPtMissedDueToEta[cent]->Clone();
    hPtRatio[cent]->Divide(hPtAll[cent]);
    SetColor(hPtRatio[cent], color[cent], marker[cent]);
    hEtaRatio[cent] = (TH1D *)hEtaMissedDueToEta[cent]->Clone();
    hEtaRatio[cent]->Divide(hEtaAll[cent]);
    SetColor(hEtaRatio[cent], color[cent], marker[cent]);
    hPtEtaRatio[cent] = (TH2D *)hPtEtaMissedDueToEta[cent]->Clone();
    hPtEtaRatio[cent]->Divide(hPtEtaAll[cent]);

    hPtDueToDeltaRRatio[cent] = (TH1D *)hPtMissedDueToDeltaR[cent]->Clone();
    hPtDueToDeltaRRatio[cent]->Divide(hPtAll[cent]);
  }

  TCanvas *d = new TCanvas("d", "d", 2000, 1000);
  d->Divide(2);
  d->cd(1);
  for (int cent = 0; cent < 3; cent++){
    hPtRatio[cent]->Draw("SAME");
    hPtRatio[cent]->GetYaxis()->SetRangeUser(0.01, 0.1);
  }

  d->cd(2);
  for (int cent = 0; cent < 3; cent++){
    hEtaRatio[cent]->Draw("SAME");
  }

  TCanvas *e = new TCanvas("e", "e", 2000, 1000);
  e->Divide(3);
  for (int cent = 0; cent < 3; cent++){
    e->cd(cent+1);
    hPtEtaRatio[cent]->Draw("COLZ");
    gPad->SetLogz();
  }

  TCanvas *e2 = new TCanvas("e2", "e2", 2000, 1000);
  e2->Divide(3);
  for (int cent = 0; cent < 3; cent++){
    e2->cd(cent+1);
    hPtDueToDeltaRRatio[cent]->Draw();
  }

  TFile *g = new TFile(Form("JetQA/NonPtInefficiency_%s.root", mode.Data()), "RECREATE");
  g->cd();

  for (int cent = 0; cent < 3; cent++){
    // hPtMissedDueToEta[cent]->Write();
    // hPtAll[cent]->Write();
    hPtRatio[cent]->Write();

    // hEtaMissedDueToEta[cent]->Write();
    // hEtaAll[cent]->Write();
    hEtaRatio[cent]->Write();

    // hPtEtaMissedDueToEta[cent]->Write();
    // hPtEtaAll[cent]->Write();
    hPtEtaRatio[cent]->Write();
  }

  g->Close();
}