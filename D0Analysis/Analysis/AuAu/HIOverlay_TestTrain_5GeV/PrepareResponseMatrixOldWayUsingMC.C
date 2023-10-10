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

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldSvd.h"

#pragma link C++ class vector<int> +;

using namespace std;

#endif

#include "BinDef.h"
#include "NewBinDef.h"

TH1D *hZ[3];
TH2D *hJetPtvZ[3];
THnSparseF *hJet[3];
THnSparseF *hJetZNorm[3];
THnSparseF *hJetZNormPtNorm[3];
TH2D *Weight[3];
TH2D *fWeight[3];
RooUnfoldResponse *resp[3]; //Response
TH2D *fMeas[3];
TH2D *fTrue[3];
TH2D *fMiss[3];

RooUnfoldResponse *resp1D[3]; //Response
TH2D *fResp1D[3];
TH1D *fMeas1D[3];
TH1D *fTrue1D[3];
TH1D *fMiss1D[3];
TH1D *fFake1D[3];

TH1D *fMeasZ1D[3];
TH1D *fTrueZ1D[3];
TH1D *fMissZ1D[3];

TH1D *hDiffJetPt[3]; //Reco - MC
TH1D *hDiffJetEta[3]; //Reco - MC
TH1D *hDiffJetPhi[3]; //Reco - MC

TH1D *hMCD0Pt[3]; //Reco - MC
TH1D *hRecoD0Pt[3]; //Reco - MC

void Method(int cent = 0, int iteration = 3, TString DirName = "MCMCUnf"){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TFile *f = new TFile("MCResponse.root");

  cout << f->GetName() <<  endl;
  if (cent == 0) f->cd("SimJetSaverSmear_Central");
  else if (cent == 1) f->cd("SimJetSaverSmear_MidCentral");
  else if (cent == 2) f->cd("SimJetSaverSmear_Peripheral");

  cout << gDirectory->GetName() << endl;

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

  // JetTree->SetBranchAddress("Centrality", &Centrality);
  JetTree->SetBranchAddress("MCPrimaryVertex", &MCPrimaryVertex);
  JetTree->SetBranchAddress("MCJetPt", &MCJetPt);
  JetTree->SetBranchAddress("MCJetEta", &MCJetEta);
  JetTree->SetBranchAddress("MCJetPhi", &MCJetPhi);
  JetTree->SetBranchAddress("MCJetArea", &MCJetArea);
  JetTree->SetBranchAddress("MCJetE", &MCJetE);
  JetTree->SetBranchAddress("MCJetNConst", &MCJetNConst);
  // JetTree->SetBranchAddress("MCD0Z", &MCD0Z);
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
  // JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
  JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConst);
  // JetTree->SetBranchAddress("RecoD0Z", &RecoD0Z);
  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  TFile *PYTHIA = new TFile("PYTHIA_Pt_5_20.root", "READ");
  TH2F  *PYTHIA2D      = (TH2F *) PYTHIA->Get("PYTHIA");
  TH1F  *FONLLvPYTHIAWeights = (TH1F *) PYTHIA->Get("FONLLvPYTHIAWeights");

  cout << "PYTHIA Integrals = " << PYTHIA2D->Integral() << "\t" << FONLLvPYTHIAWeights->Integral() << endl;

  TFile *FONLL = new TFile("FONLL_Pt_5_20.root", "READ");
  TH1D  *FONLLCurve = (TH1D *)FONLL->Get("FONLL");

  cout << "FONLL Integrals = " << FONLLCurve->Integral() << endl;

  TFile *SmearFile = new TFile("TestSinglePartEmbedding/SingleParticleEmbedding.root", "READ");
  SmearFile->cd();
  TH1D *Smear[3];

  for (int i = 0; i < 3; i++){
    Smear[i] = (TH1D *)gDirectory->Get(Form("hDiffJetPtScaled_%i", i));
    SetName(Smear[i], Form("Smear_%i", i));
  }


  cout << JetTree->GetEntries() << endl;

  for (int i = cent; i < cent + 1; i++){
    hZ[i] = new TH1D(Form("hZ_%i", i), Form("hZ_%i", i), nz_gen_bins, z_gen_low, z_gen_high);
    hJetPtvZ[i] = new TH2D(Form("hJetPtvZ_%i", i), Form("hJetPtvZ_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
    hJet[i] = new THnSparseF(Form("hJet_%i", i), Form("hJet_%i", i), ndim, nbins, xmin, xmax);
    hJetZNorm[i] = new THnSparseF(Form("hJetZNorm_%i", i), Form("hJetZNorm_%i", i), ndim, nbins, xmin, xmax);
    hJetZNormPtNorm[i] = new THnSparseF(Form("hJetZNormPtNorm_%i", i), Form("hJetZNormPtNorm_%i", i), ndim, nbins, xmin, xmax);

    fMeas[i] = new TH2D(Form("fMeas_Cent_%i", i), Form("fMeas_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);
    fTrue[i] = new TH2D(Form("fTrue_Cent_%i", i), Form("fTrue_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
    fMiss[i] = new TH2D(Form("fMiss_Cent_%i", i), Form("fMiss_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

    resp[i] = new RooUnfoldResponse(Form("Resp_%i", i), Form("Resp_%i", i));
    resp[i]->Setup(fMeas[i], fTrue[i]); //Setup Response Matrix Definition

    // fMeas1D[i] = new TH1D(Form("fMeas1D_Cent_%i", i), Form("fMeas1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    // fTrue1D[i] = new TH1D(Form("fTrue1D_Cent_%i", i), Form("fTrue1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
    // fMiss1D[i] = new TH1D(Form("fMiss1D_Cent_%i", i), Form("fMiss1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
    // fFake1D[i] = new TH1D(Form("fFake1D_Cent_%i", i), Form("fFake1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    // fResp1D[i] = new TH2D(Form("fResp1D_Cent_%i", i), Form("fResp1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

    fMeas1D[i] = new TH1D(Form("fMeas1D_Cent_%i", i), Form("fMeas1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    fTrue1D[i] = new TH1D(Form("fTrue1D_Cent_%i", i), Form("fTrue1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
    fMiss1D[i] = new TH1D(Form("fMiss1D_Cent_%i", i), Form("fMiss1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

    fMeasZ1D[i] = new TH1D(Form("fMeasZ1D_Cent_%i", i), Form("fMeasZ1D_Cent_%i", i), nz_bins, z_low, z_high);
    fTrueZ1D[i] = new TH1D(Form("fTrueZ1D_Cent_%i", i), Form("fTrueZ1D_Cent_%i", i), nz_gen_bins, z_gen_low, z_gen_high);
    fMissZ1D[i] = new TH1D(Form("fMissZ1D_Cent_%i", i), Form("fMissZ1D_Cent_%i", i), nz_gen_bins, z_gen_low, z_gen_high);

    resp1D[i] = new RooUnfoldResponse(Form("Resp1D_%i", i), Form("Resp1D_%i", i));
    resp1D[i]->Setup(fMeas1D[i], fTrue1D[i]); //Setup Response Matrix Definition

    hDiffJetPt[i] = new TH1D(Form("hDiffJetPt_%i", i), Form("hDiffJetPt_%i", i), 100, -50, 50);
    hDiffJetEta[i] = new TH1D(Form("hDiffJetEta_%i", i), Form("hDiffJetEta_%i", i), 40, -1., 1.);
    hDiffJetPhi[i] = new TH1D(Form("hDiffJetPhi_%i", i), Form("hDiffJetPhi_%i", i), 40, -1., 1.);

    hMCD0Pt[i] = new TH1D(Form("hMCD0Pt_%i", i), Form("hMCD0Pt_%i", i), 20, 0, 10);
    hRecoD0Pt[i] = new TH1D(Form("hRecoD0Pt_%i", i), Form("hRecoD0Pt_%i", i), 20, 0, 10);

    hJet[i]->Sumw2();
    hJetZNorm[i]->Sumw2();
    hJetZNormPtNorm[i]->Sumw2();

    hJet[i]->CalculateErrors();
    hJetZNorm[i]->CalculateErrors();
    hJetZNormPtNorm[i]->CalculateErrors();

    Weight[i] = new TH2D(Form("Weight_%i", i), Form("Weight_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
    fWeight[i] = new TH2D(Form("fWeight_%i", i), Form("fWeight_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

    for (int binx = 1; binx <= Weight[i]->GetNbinsX(); binx++ ){
      TH2D *tmp1 = (TH2D *)PYTHIA2D->Clone();
      tmp1->GetXaxis()->SetRange(binx, binx);
      TH1D *tmp2 = (TH1D *)tmp1->ProjectionY();
      tmp2->Scale(1./tmp2->Integral());
      for (int biny = 1; biny <= Weight[i]->GetNbinsY(); biny++ ){          
        Weight[i]->SetBinContent(binx, biny, tmp2->GetBinContent(biny));
      }
    }

    cout << "Weight Integrals = " << Weight[i]->Integral() << endl;
  }
  
  cout << JetTree->GetEntries() << endl;

  int nentries = JetTree->GetEntries();

  int lowlimit = 0;
  int highlimit = nentries;

  int counter = 0;

  double recooutside[3] = {0};

  double etasmearfactors[3] = {0.11, 0.09, 0.05}; // Eta Smear Factors From Heavy Ion Overlay
  double phismearfactors[3] = {0.15, 0.13 ,0.09}; // Phi Smear Factors From Heavy Ion Overlay

  TRandom3 *r = new TRandom3(0);

  for (int i = lowlimit; i < highlimit; i++){
    
    if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    // cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    int centhistogramtofill = cent;

    if (centhistogramtofill < 0) continue;

    // MCD0Phi = MCTrackPhi[MCD0Index];
    // RecoD0Phi = RecoTrackPhi[RecoD0Index];

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    // double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

    // if (mcz > 1.0) cout << MCJetPt << "\t" << MCD0Pt << mcz << endl;
    if (mcz == 1.0) mcz = 0.999; // Padding the boundaries

    // RecoJetCorrPt = RecoJetPt + r->Gaus(0, ptsmearfactors[centhistogramtofill]);
    RecoJetCorrPt = RecoJetPt + Smear[centhistogramtofill]->GetRandom();
    // RecoJetCorrPt = MCJetPt + r->Gaus(ptoffsetfactors[centhistogramtofill], ptsmearfactors[centhistogramtofill]);
    // RecoJetEta    = RecoJetEta + (RecoJetEta - MCJetEta)/(abs(RecoJetEta - MCJetEta)) * abs(r->Gaus(0, etasmearfactors[centhistogramtofill]));
    // RecoJetPhi    = RecoJetPhi + (RecoJetPhi - MCJetPhi)/(abs(RecoJetPhi - MCJetPhi)) * abs(r->Gaus(0, phismearfactors[centhistogramtofill]));

    RecoJetEta    = RecoJetEta + r->Gaus(0, etasmearfactors[centhistogramtofill]);
    RecoJetPhi    = RecoJetPhi + r->Gaus(0, phismearfactors[centhistogramtofill]);

    if(RecoJetPhi>=2.*TMath::Pi()) RecoJetPhi = RecoJetPhi - 2.*TMath::Pi();
    if(RecoJetPhi<0) RecoJetPhi = 2.*TMath::Pi() + RecoJetPhi;

    double RecoJetPx = RecoJetCorrPt*TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    // double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);
    double recoz = RecoD0Pt/RecoJetCorrPt * TMath::Cos(RecoJetPhi - RecoD0Phi);

    double RecoDeltaR = dR(dEta(RecoJetEta, RecoD0Eta), dPhi(RecoJetPhi, RecoD0Phi));
    // cout << MCJetPx << "\t" << MCD0Px << "\t" << RecoJetPx << "\t" << RecoD0Px << endl;
    // cout << MCD0Pt << "\t" << MCJetPt << "\t" << MCJetPhi << "\t" << mcz << "\t" << RecoD0Pt << "\t" << RecoJetCorrPt << "\t" << RecoJetPhi << "\t" << recoz << endl;

    // if (recoz > z_low-0.001 && recoz <= z_low) recoz = z_low + 0.001; // Padding the boundaries
    // if (recoz == 1.0) recoz = 0.99;

    //Defining response matrix using fakes and misses to check closure

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 5 && RecoD0Pt < 10;
    bool isMCZ   = mcz > z_gen_low && mcz < z_gen_high; 
    bool isRecoZ = recoz > z_low && recoz < z_high; 
    // bool isMCJetPt = MCJetPt > 5 && MCJetPt < 30;
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    // bool isRecoJetPt = RecoJetCorrPt > -7 && RecoJetCorrPt < 30; // This is used because 99% of jets with pT > 5 GeV is captured within this range
    // bool isRecoJetPt = RecoJetCorrPt > 3 && RecoJetCorrPt < 30; // 78% of jets with pT > 5 GeV is captured within this range
    bool isRecoJetPt = RecoJetCorrPt > 3;
    // bool isRecoJetPt = kTRUE;
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoDeltaR = RecoDeltaR < 0.4;

    // bool isConstGreaterThanJet = (RecoJetCorrPt - RecoD0Pt > 0);

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (RecoJetNConst==0) continue;
    
    // hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt);
    // hDiffJetEta[centhistogramtofill]->Fill(RecoJetEta - MCJetEta);
    // hDiffJetPhi[centhistogramtofill]->Fill(RecoJetPhi - MCJetPhi);

    // hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt);
    // hRecoD0Pt[centhistogramtofill]->Fill(RecoD0Pt);

    // isRecoZ = kTRUE;
    // isMCZ = kTRUE;

    bool isOK = isRecoZ && isMCZ && isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta && isRecoDeltaR;
    bool isMISS = !isOK && isMCZ && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOK && isRecoZ && isRecoJetPt && isRecoJetEta && isRecoDeltaR;

    double w = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));
    // double w = 1.;

    // cout << w << endl;

    if (!isOK && !isMISS) continue;
    // hCentrality->Fill(Centrality, w);

    hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt, w);
    hDiffJetEta[centhistogramtofill]->Fill(RecoJetEta - MCJetEta, w);
    hDiffJetPhi[centhistogramtofill]->Fill(RecoJetPhi - MCJetPhi, w);

    hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt, w);
    hRecoD0Pt[centhistogramtofill]->Fill(RecoD0Pt, w);

    double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
    hJet[centhistogramtofill]->Fill(fvalue, w);
  }

  // This fixes the Z distribution to PYTHIA (this is necessary when there are centrality weights in the picture.)

  // Changing the z distribution here
  for (int gbinx = 1; gbinx <= njpt_gen_bins; gbinx++ ){
    for (int gbiny = 1; gbiny <= nz_gen_bins; gbiny++ ){
      double scalefactor = 0;
      double valintegral = 0;
      // scalefactor = Weight[cent]->GetBinContent(gbinx, gbiny);

      for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
        for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){
          int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
          int bin = hJet[cent]->GetBin(xbin);
          valintegral += hJet[cent]->GetBinContent(bin);
        }
      }

      if (valintegral > 0) {
        scalefactor = Weight[cent]->GetBinContent(gbinx, gbiny)/valintegral;
      }

      // cout << scalefactor << endl;

      fWeight[cent]->SetBinContent(gbinx, gbiny, scalefactor);

      for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
        for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){

          int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
          int binZnorm = hJetZNorm[cent]->GetBin(xbin);
          int bin = hJet[cent]->GetBin(xbin);
          
          double content = hJetZNorm[cent]->GetBinContent(binZnorm);
          double err = hJetZNorm[cent]->GetBinError(binZnorm);

          double newbincontent = content + hJet[cent]->GetBinContent(bin)*scalefactor;
          double newbinerr = err + hJet[cent]->GetBinError(bin)*scalefactor;
          // hJetZNorm[cent]->SetBinContent(binZnorm, newbincontent);
          // hJetZNorm[cent]->SetBinError(binZnorm, newbinerr);

          double mid[ndim];
          mid[0] = hJetZNorm[cent]->GetAxis(0)->GetBinCenter(gbinx);
          mid[1] = hJetZNorm[cent]->GetAxis(1)->GetBinCenter(gbiny);
          mid[2] = hJetZNorm[cent]->GetAxis(2)->GetBinCenter(rbinx);
          mid[3] = hJetZNorm[cent]->GetAxis(3)->GetBinCenter(rbiny);

          hJetZNorm[cent]->Fill(mid, hJet[cent]->GetBinContent(bin)*scalefactor);
          hJetZNorm[cent]->SetBinError(binZnorm, newbinerr);
        }
      }
    }
  }

  // This fixes the pT distribution to either FONLL or PYTHIA based on what we want.
  for (int gbinx = 1; gbinx <= njpt_gen_bins; gbinx++){
    double scalefactor = 0;
    double origIntegral = 0;
    double scalIntegral = 0;

    // origIntegral = PYTHIAPtCurve->GetBinContent(gbinx);
    

    for (int gbiny = 1; gbiny <= nz_gen_bins; gbiny++ ){
      for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
        for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){

          int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
          int bin = hJet[cent]->GetBin(xbin);
          int binZnorm = hJetZNorm[cent]->GetBin(xbin);
          cout << gbinx << "\t" << gbiny << "\t" << rbinx << "\t" << rbiny << "\t" << hJet[cent]->GetBinContent(bin) << endl;
          origIntegral+=hJet[cent]->GetBinContent(bin);
          scalIntegral+=hJetZNorm[cent]->GetBinContent(binZnorm);
        }
      }
    }

    cout << "Jet pT Bin = " << gbinx << "\t" << origIntegral << "\t" << scalIntegral << endl;

    scalefactor = (scalIntegral > 0) ? origIntegral/scalIntegral : 0;
    // cout << scalefactor << endl;

    for (int gbiny = 1; gbiny <= nz_gen_bins; gbiny++ ){
      for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
        for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){

          int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
          int binZnorm = hJetZNorm[cent]->GetBin(xbin);
          int binZnormPtnorm = hJetZNormPtNorm[cent]->GetBin(xbin);
          
          double content = hJetZNormPtNorm[cent]->GetBinContent(binZnormPtnorm);
          double err = hJetZNormPtNorm[cent]->GetBinError(binZnormPtnorm);

          double newbincontent = content + hJetZNorm[cent]->GetBinContent(binZnorm)*scalefactor;
          double newbinerr = err + hJetZNorm[cent]->GetBinError(binZnorm)*scalefactor;
          hJetZNormPtNorm[cent]->SetBinContent(binZnormPtnorm, newbincontent);
          hJetZNormPtNorm[cent]->SetBinError(binZnormPtnorm, newbinerr);
        }
      }
    }
  }

  //Fill RooUnfoldResponse Object

  int* coord = new int[ndim];

  int nbin = hJetZNormPtNorm[cent]->GetNbins();

  cout << "Bins = " << nbin << endl;

  for (int bin = 0; bin < nbin; bin++){
    Double_t w = hJetZNormPtNorm[cent]->GetBinContent(bin, coord);
    Double_t pttrue = hJetZNormPtNorm[cent]->GetAxis(0)->GetBinCenter(coord[0]);
    Double_t ztrue = hJetZNormPtNorm[cent]->GetAxis(1)->GetBinCenter(coord[1]);
    Double_t ptdet =  hJetZNormPtNorm[cent]->GetAxis(2)->GetBinCenter(coord[2]); 
    Double_t zdet =  hJetZNormPtNorm[cent]->GetAxis(3)->GetBinCenter(coord[3]);

    if (zdet >= z_low && zdet <= z_high
    &&  ztrue >= z_gen_low && ztrue <= z_gen_high
    && ptdet >= 3 && ptdet <= 30
    && pttrue >= 5 && pttrue <= 20
    ){
      resp[cent]->Fill(ptdet, zdet, pttrue, ztrue, w);
      resp1D[cent]->Fill(ptdet, pttrue, w);
      fTrue[cent]->Fill(pttrue, ztrue, w);
      fTrue1D[cent]->Fill(pttrue, w);
      fTrueZ1D[cent]->Fill(ztrue, w);

      fMeas[cent]->Fill(ptdet, zdet, w);
      fMeas1D[cent]->Fill(ptdet, w);
      fMeasZ1D[cent]->Fill(zdet, w);
    }
    else{
      resp[cent]->Miss(pttrue, ztrue, w);
      resp1D[cent]->Miss(pttrue, w);
      fMiss[cent]->Fill(pttrue, ztrue, w);
      fMiss1D[cent]->Fill(pttrue, w);

      fTrue[cent]->Fill(pttrue, ztrue, w);
      fTrue1D[cent]->Fill(pttrue, w);
      fTrueZ1D[cent]->Fill(ztrue, w);
    } 
  }

  delete [] coord;

  // cout << "Entries = " << fMiss1D[cent]->Integral() << "\t" << fFake1D[cent]->Integral() << "\t" << fResp1D[cent]->Integral() << endl;

  TString filename;
  filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins);
  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "UPDATE");
  outfile->cd();

  hZ[cent]->Write();
  hJetPtvZ[cent]->Write();
  hJet[cent]->Write();
  hJetZNorm[cent]->Write();
  hJetZNormPtNorm[cent]->Write();

  fMeas[cent]->Write();
  fTrue[cent]->Write();
  fMiss[cent]->Write();

  fMeas1D[cent]->Write();
  fTrue1D[cent]->Write();
  fMiss1D[cent]->Write();

  fMeasZ1D[cent]->Write();
  fTrueZ1D[cent]->Write();
  fMissZ1D[cent]->Write();

  gDirectory->WriteObject(resp[cent], Form("Resp_%i", cent));
  gDirectory->WriteObject(resp1D[cent], Form("Resp1D_%i", cent));

  hDiffJetPt[cent]->Write();
  hDiffJetEta[cent]->Write();
  hDiffJetPhi[cent]->Write();

  hMCD0Pt[cent]->Write();
  hRecoD0Pt[cent]->Write();

  Weight[cent]->Write();
  fWeight[cent]->Write();
  
  outfile->Close();
}

void PrepareResponseMatrixOldWayUsingMC(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf"){
  TString filename;
  filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins);
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->Close();

  for (int cent = 0; cent < 3; cent++){
    Method(cent, iteration, DirName);
  }
}