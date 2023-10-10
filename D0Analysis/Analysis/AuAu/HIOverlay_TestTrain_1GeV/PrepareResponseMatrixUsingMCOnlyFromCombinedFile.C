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
RooUnfoldResponse *resp[3]; //Response
TH2D *fMeas[3];
TH2D *fTrue[3];
TH2D *fMiss[3];
TH2D *fFake[3];

TH2D *hMeasured[3];
TH2D *hTrue[3];
TH1D *hMeasured1D[3];
TH1D *hTrue1D[3];

RooUnfoldResponse *resp1D[3]; //Response
TH2D *fResp1D[3];
TH1D *fMeas1D[3];
TH1D *fTrue1D[3];
TH1D *fMiss1D[3];
TH1D *fFake1D[3];

TH2D *hMiss1D[3]; //These misses are because of eta-phi mismatch and not a pT mismatch

TH1D *hDiffJetPt[3]; //Reco - MC
TH1D *hDiffJetPtInGenPtBins[3][nbins_jpt]; //Reco - MC
TH1D *hDiffJetEta[3]; //Reco - MC
TH1D *hDiffJetEtaBeforeSmear[3]; //Reco - MC
TH1D *hDiffJetEtaAfterSmear[3]; //Reco - MC
TH1D *hDiffJetPhi[3]; //Reco - MC

TH1D *hMCD0Pt[3]; //Reco - MC
TH1D *hRecoD0Pt[3]; //Reco - MC


void Method(int iteration = 3, TString DirName = "MCMCUnf", int mode = 0){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TFile *f = new TFile("HIResponse_MCIncluded_Feb6.root");

  cout << f->GetName() <<  endl;
  f->cd("HIJetSaver");
  
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

  JetTree->SetBranchAddress("Centrality", &Centrality);
  JetTree->SetBranchAddress("MCPrimaryVertex", &MCPrimaryVertex);
  JetTree->SetBranchAddress("RecoPrimaryVertex", &RecoPrimaryVertex);
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
 
  JetTree->SetBranchAddress("RecoJetPtFromPYTHIA", &RecoJetPt);
  JetTree->SetBranchAddress("RecoJetEtaFromPYTHIA", &RecoJetEta);
  JetTree->SetBranchAddress("RecoJetPhiFromPYTHIA", &RecoJetPhi);
  JetTree->SetBranchAddress("RecoJetAreaFromPYTHIA", &RecoJetArea);
  JetTree->SetBranchAddress("RecoJetEFromPYTHIA", &RecoJetE);
  JetTree->SetBranchAddress("RecoJetNConstFromPYTHIA", &RecoJetNConst);
  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);


  // TFile *SmearFile = new TFile("TestSinglePartEmbedding/SingleParticleEmbedding.root");
  TFile *SmearFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/BKGHistogramsFINALDec14_22.root");
  SmearFile->cd();
  TH1D *Smear[3];

  Smear[0] = (TH1D *)gDirectory->Get("hBKG_0_10_Weighted");
  Smear[1] = (TH1D *)gDirectory->Get("hBKG_10_40_Weighted");
  Smear[2] = (TH1D *)gDirectory->Get("hBKG_40_80_Weighted");
  // for (int i = 0; i < 3; i++){
  //   Smear[i] = (TH1D *)gDirectory->Get(Form("hDiffJetPtScaled_%i", i));
  //   SetName(Smear[i], Form("Smear_%i", i));
  // }


  // TFile *ZFile = new TFile("ZFile.root", "UPDATE");
  // TH1D *ZWanted = (TH1D *)ZFile->Get(Form("Z_%i_%i", SUPERITERATION, iteration));

  cout << JetTree->GetEntries() << endl;

  for (int i = 0; i < 3; i++){
    hZ[i] = new TH1D(Form("hZ_%i", i), Form("hZ_%i", i), nz_gen_bins, z_gen_low, z_gen_high);
    hJetPtvZ[i] = new TH2D(Form("hJetPtvZ_%i", i), Form("hJetPtvZ_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
    hJet[i] = new THnSparseF(Form("hJet_%i", i), Form("hJet_%i", i), ndim, nbins, xmin, xmax);
    hJetZNorm[i] = new THnSparseF(Form("hJetZNorm_%i", i), Form("hJetZNorm_%i", i), ndim, nbins, xmin, xmax);
    hJetZNormPtNorm[i] = new THnSparseF(Form("hJetZNormPtNorm_%i", i), Form("hJetZNormPtNorm_%i", i), ndim, nbins, xmin, xmax);

    fMeas[i] = new TH2D(Form("fMeas_Cent_%i", i), Form("fMeas_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);
    fTrue[i] = new TH2D(Form("fTrue_Cent_%i", i), Form("fTrue_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
    fMiss[i] = new TH2D(Form("fMiss_Cent_%i", i), Form("fMiss_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
    fFake[i] = new TH2D(Form("fFake_Cent_%i", i), Form("fFake_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);

    resp[i] = new RooUnfoldResponse(Form("Resp_%i", i), Form("Resp_%i", i));
    resp[i]->Setup(fMeas[i], fTrue[i]); //Setup Response Matrix Definition

    fMeas1D[i] = new TH1D(Form("fMeas1D_Cent_%i", i), Form("fMeas1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    fTrue1D[i] = new TH1D(Form("fTrue1D_Cent_%i", i), Form("fTrue1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
    fMiss1D[i] = new TH1D(Form("fMiss1D_Cent_%i", i), Form("fMiss1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
    fFake1D[i] = new TH1D(Form("fFake1D_Cent_%i", i), Form("fFake1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    fResp1D[i] = new TH2D(Form("fResp1D_Cent_%i", i), Form("fResp1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

    resp1D[i] = new RooUnfoldResponse(Form("Resp1D_%i", i), Form("Resp1D_%i", i));
    resp1D[i]->Setup(fMeas1D[i], fTrue1D[i]); //Setup Response Matrix Definition

    hMiss1D[i] = new TH2D(Form("hMiss1D_Cent_%i", i), Form("hMiss1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

    hDiffJetPt[i] = new TH1D(Form("hDiffJetPt_%i", i), Form("hDiffJetPt_%i", i), 100, -50, 50);
    hDiffJetEta[i] = new TH1D(Form("hDiffJetEta_%i", i), Form("hDiffJetEta_%i", i), 40, -1., 1.);
    hDiffJetEtaBeforeSmear[i] = new TH1D(Form("hDiffJetEtaBeforeSmear_%i", i), Form("hDiffJetEtaBeforeSmear_%i", i), 40, -1., 1.);
    hDiffJetEtaAfterSmear[i] = new TH1D(Form("hDiffJetEtaAfterSmear_%i", i), Form("hDiffJetEtaAfterSmear_%i", i), 40, -1., 1.);
    hDiffJetPhi[i] = new TH1D(Form("hDiffJetPhi_%i", i), Form("hDiffJetPhi_%i", i), 40, -1., 1.);

    for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
      hDiffJetPtInGenPtBins[i][ptbin] = new TH1D(Form("hDiffJetPt_%i_%i", i, ptbin), Form("hDiffJetPt_%i_%i", i, ptbin), 100, -50, 50);
    }

    hMCD0Pt[i] = new TH1D(Form("hMCD0Pt_%i", i), Form("hMCD0Pt_%i", i), 20, 0, 10);
    hRecoD0Pt[i] = new TH1D(Form("hRecoD0Pt_%i", i), Form("hRecoD0Pt_%i", i), 20, 0, 10);

    hJet[i]->Sumw2();
    hJetZNorm[i]->Sumw2();
    hJetZNormPtNorm[i]->Sumw2();

    hJet[i]->CalculateErrors();
    hJetZNorm[i]->CalculateErrors();
    hJetZNormPtNorm[i]->CalculateErrors();

    Weight[i] = new TH2D(Form("Weight_%i", i), Form("Weight_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

    for (int binx = 1; binx <= Weight[i]->GetNbinsX(); binx++ ){
      for (int biny = 1; biny <= Weight[i]->GetNbinsY(); biny++ ){
        Weight[i]->SetBinContent(binx, biny, 1.);
      }
    }
  }
  
  cout << JetTree->GetEntries() << endl;

  int nentries = JetTree->GetEntries();

  int lowlimit = 0;
  int highlimit = nentries;
  // int highlimit = 100;

  // nentries = 100;

  // if (SUPERITERATION == 0){lowlimit = 0; highlimit = nentries;}
  // else {lowlimit = nentries/2; highlimit = nentries;}

  int counter = 0;

  double recooutside[3] = {0};

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

  
  double etasmearfactors[3] = {0.11, 0.09, 0.05}; // Eta Smear Factors From Heavy Ion Overlay
  double phismearfactors[3] = {0.15, 0.13 ,0.09}; // Phi Smear Factors From Heavy Ion Overlay

  TRandom* r = new TRandom(0);
  gRandom->SetSeed(0);

  TFile *PYTHIA = new TFile("PYTHIA_Pt_5_20.root", "READ");
  // TFile *PYTHIA = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/FONLLWeights_New.root", "READ");
  TH1F  *FONLLvPYTHIAWeights = (TH1F *) PYTHIA->Get("FONLLWeights");

  cout << "Mode == " << mode << endl;

  TH1D *tmppT = new TH1D("tmppT", "tmppT", 40, 0, 40);

  TH1D *plottingpT = new TH1D("plottingpT", "plottingpT", nbins_jpt, binning_jpt);

  for (int i = lowlimit; i < highlimit; i++){
    
    if (mode == 1) {if (i%2 != 0) continue;}
    else if (mode == 2) {if (i%2 == 0) continue;}

    
    // cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

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

    int ptbin = tmppT->FindBin(RecoJetPt);
    // RecoJetCorrPt = RecoJetPt + r->Gaus(0, ptsmearfactors[centhistogramtofill]);
    RecoJetCorrPt = RecoJetPt + dPt[centhistogramtofill][ptbin-1]->GetRandom();
    // RecoJetCorrPt = MCJetPt + r->Gaus(ptoffsetfactors[centhistogramtofill], ptsmearfactors[centhistogramtofill]);
    // RecoJetEta    = RecoJetEta + (RecoJetEta - MCJetEta)/(abs(RecoJetEta - MCJetEta)) * abs(r->Gaus(0, etasmearfactors[centhistogramtofill]));
    // RecoJetPhi    = RecoJetPhi + (RecoJetPhi - MCJetPhi)/(abs(RecoJetPhi - MCJetPhi)) * abs(r->Gaus(0, phismearfactors[centhistogramtofill]));
    double eta_mean_val = Eta_Mean[centhistogramtofill]->GetBinContent(ptbin);
    double eta_sigma_val = Eta_Sigma[centhistogramtofill]->GetBinContent(ptbin);

    double phi_mean_val = Phi_Mean[centhistogramtofill]->GetBinContent(ptbin);
    double phi_sigma_val = Phi_Sigma[centhistogramtofill]->GetBinContent(ptbin);


    // RecoJetEta    = RecoJetEta + r->Gaus(0, etasmearfactors[centhistogramtofill]);
    // RecoJetPhi    = RecoJetPhi + r->Gaus(0, phismearfactors[centhistogramtofill]);

    double oldRecoJetEta = RecoJetEta;

    RecoJetEta    = RecoJetEta + r->Gaus(eta_mean_val, eta_sigma_val);
    RecoJetPhi    = RecoJetPhi + r->Gaus(phi_mean_val, phi_sigma_val);

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

    bool isMCD0Pt = MCD0Pt > 4 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 4 && RecoD0Pt < 10;
    bool isMCZ   = mcz > z_gen_low && mcz < z_gen_high; 
    bool isRecoZ = recoz > z_low && recoz < z_high; 
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    bool isRecoJetPt = RecoJetCorrPt > 3 && RecoJetCorrPt < 30; // 78% of jets with pT > 5 GeV is captured within this range
    bool isMCJetEta = abs(MCJetEta) < 0.601;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoDeltaR = RecoDeltaR < 0.4;

    // bool isConstGreaterThanJet = (RecoJetCorrPt - RecoD0Pt > 0);

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (RecoJetNConst==0) continue;
    // if (!isMCJetPt) continue;

    // if (RecoJetCorrPt >= jetpt_low && RecoJetCorrPt <= jetpt_high
    //   && MCJetPt >= jetpt_gen_low && MCJetPt <= jetpt_gen_high
    // ){
    //   resp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt);
    // }
    // else{
    //   resp1D[centhistogramtofill]->Miss(MCJetPt);
    //   fMiss1D[centhistogramtofill]->Fill(MCJetPt);
    // }

    isRecoZ = kTRUE;
    isMCZ = kTRUE;

    // bool isOK = isRecoZ && isMCZ && isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta && isRecoDeltaR;
    // bool isMISS = !isOK && isMCZ && isMCJetPt && isMCJetEta;
    // bool isFAKE = !isOK && isRecoZ && isRecoJetPt && isRecoJetEta && isRecoDeltaR;

    bool isOK   = isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta;
    bool isMISS = !isOK && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOK && isRecoJetPt && isRecoJetEta;

    double w = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));
    // double w = 1.;

    if (!isOK && !isMISS && !isFAKE) continue;

    // cout << MCJetPt << "\t" << MCJetEta << "\t" << MCJetPhi << "\t" << RecoJetCorrPt << "\t" << RecoJetEta << "\t" << RecoJetPhi << endl;

    // cout << Form("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", MCJetPt, MCJetEta, MCJetPhi, RecoJetCorrPt, RecoJetEta, RecoJetPhi);

    hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt);
    hDiffJetEta[centhistogramtofill]->Fill(dEta(RecoJetEta, MCJetEta));
    hDiffJetEtaBeforeSmear[centhistogramtofill]->Fill(oldRecoJetEta - MCJetEta);
    hDiffJetEtaAfterSmear[centhistogramtofill]->Fill(RecoJetEta - MCJetEta);
    hDiffJetPhi[centhistogramtofill]->Fill(dPhi(RecoJetPhi, MCJetPhi));

    int genptbin = plottingpT->FindBin(MCJetPt);

    hDiffJetPtInGenPtBins[centhistogramtofill][genptbin-1]->Fill(RecoJetCorrPt - MCJetPt);

    hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt, w);
    hRecoD0Pt[centhistogramtofill]->Fill(RecoD0Pt, w);

    if (isOK){
      resp[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);
      fResp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);

      fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);
    }
    else if (isMISS){
      resp[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      fMiss[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Miss(MCJetPt, w);
      fMiss1D[centhistogramtofill]->Fill(MCJetPt, w);

      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);
    }
    else if (isFAKE){
      resp[centhistogramtofill]->Fake(RecoJetCorrPt, recoz, w);
      resp1D[centhistogramtofill]->Fake(RecoJetCorrPt, w);

      fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);

      fFake[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fFake1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);
    }

    if (isMISS) hMiss1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);

    // if (!isConstGreaterThanJet) continue;
    // if (MCJetNConst <= 2 ) continue;
    // if (mcz > 0.88) cout << MCJetPt << "\t" << MCD0Pt << "\t" << MCJetNConst << "\t" << mcz << endl;

    if (!isMCJetPt) continue;
    if (!isRecoJetPt) continue;
    if (!isMCJetEta) continue;
    if (!isRecoJetEta) continue;
    if (!isMCZ) continue;

    if (!isRecoZ) recooutside[centhistogramtofill]++;

    double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
    hJet[centhistogramtofill]->Fill(fvalue);
    hZ[centhistogramtofill]->Fill(mcz);
    hJetPtvZ[centhistogramtofill]->Fill(MCJetPt, mcz);
  }

  cout << "Entries = " << fMeas[0]->Integral() << "\t" << fTrue[0]->Integral() << "\t" << fMeas1D[0]->Integral() << "\t" << fTrue1D[0]->Integral() << endl;
  cout << "Entries = " << fMeas[1]->Integral() << "\t" << fTrue[1]->Integral() << "\t" << fMeas1D[1]->Integral() << "\t" << fTrue1D[1]->Integral() << endl;
  cout << "Entries = " << fMeas[2]->Integral() << "\t" << fTrue[2]->Integral() << "\t" << fMeas1D[2]->Integral() << "\t" << fTrue1D[2]->Integral() << endl;

  TString filename;

  switch (mode){
    case 1:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    case 2:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    default:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins);
  }

  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "UPDATE");
  outfile->cd();

  for (int cent = 0; cent < 3; cent++){
    hZ[cent]->Write();
    hJetPtvZ[cent]->Write();
    hJet[cent]->Write();

    fMeas[cent]->Write();
    fTrue[cent]->Write();
    fFake[cent]->Write();

    fMeas1D[cent]->Write();
    fTrue1D[cent]->Write();

    fMiss[cent]->Write();

    fResp1D[cent]->Write();
    fMiss1D[cent]->Write();
    fFake1D[cent]->Write();

    gDirectory->WriteObject(resp[cent], Form("Resp_%i", cent));
    gDirectory->WriteObject(resp1D[cent], Form("Resp1D_%i", cent));

    hMiss1D[cent]->Write();

    hDiffJetPt[cent]->Write();
    hDiffJetEta[cent]->Write();
    hDiffJetPhi[cent]->Write();

    hDiffJetEtaBeforeSmear[cent]->Write();
    hDiffJetEtaAfterSmear[cent]->Write();

    hMCD0Pt[cent]->Write();
    hRecoD0Pt[cent]->Write();
  }

  for (int cent = 0; cent < 3; cent++){
    for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
      hDiffJetPtInGenPtBins[cent][ptbin]->Write();
    }
  }

  outfile->Close();
}

void PrepareResponseMatrixUsingMCOnlyFromCombinedFile(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", int mode = 0){

  TString filename;

  switch (mode){
    case 1:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    case 2:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    default:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins);
  }
  
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->Close();

  Method(iteration, DirName, mode);

  // for (int cent = 0; cent < 3; cent++){
    
  // }

  // for (int cent = 0; cent < 1; cent++){
  //   Method(cent, iteration, DirName, mode);
  // }
}