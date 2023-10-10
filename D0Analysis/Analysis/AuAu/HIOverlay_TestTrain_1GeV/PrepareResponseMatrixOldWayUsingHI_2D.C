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

#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/diptanilroy/ROOT_INSTALL/RooUnfold/src/RooUnfoldSvd.h"

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", bool makeInitialDistributionDifferent = kFALSE, int mode = 0){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);


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

  // TFile *ZFile = new TFile("ZFile.root", "UPDATE");
  // TH1D *ZWanted = (TH1D *)ZFile->Get(Form("Z_%i_%i", SUPERITERATION, iteration));

  // // TFile *SmearFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/BKGHistogramsFINALDec14_22.root");
  // TFile *SmearFile = new TFile("Jan26_FONLL_MC_CombinedFile4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root");
  // SmearFile->cd();
  // TH1D *PtSmear[3][nbins_jpt];
  // // TH1D *Smear[3];
  // TH1D *EtaSmear[3];
  // TH1D *PhiSmear[3];

  // // Smear[0] = (TH1D *)gDirectory->Get("hBKG_0_10_Weighted");
  // // Smear[1] = (TH1D *)gDirectory->Get("hBKG_10_40_Weighted");
  // // Smear[2] = (TH1D *)gDirectory->Get("hBKG_40_80_Weighted");
  // for (int i = 0; i < 3; i++){
  //   for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
  //     PtSmear[i][ptbin] = (TH1D *)gDirectory->Get(Form("hDiffJetPt_%i_%i", i, ptbin));
  //     cout << PtSmear[i][ptbin]->Integral() << endl;
  //     SetName(PtSmear[i][ptbin], Form("Smear_%i_%i", i, ptbin));
  //   }

  //   // Smear[i] = (TH1D *)gDirectory->Get(Form("hDiffJetPt_%i", i));
  //   EtaSmear[i] = (TH1D *)gDirectory->Get(Form("hDiffJetEta_%i", i));
  //   PhiSmear[i] = (TH1D *)gDirectory->Get(Form("hDiffJetPhi_%i", i));
    
  //   SetName(EtaSmear[i], Form("EtaSmear_%i", i));
  //   SetName(PhiSmear[i], Form("PhiSmear_%i", i));
  // }

  // TF1 *PtSmearFit[3][nbins_jpt];

  // TFile *PtSmearFitFile = new TFile("PtSmear_Fit.root");
  // PtSmearFitFile->cd();

  // for (int i = 0; i < 3; i++){
  //   for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
  //     PtSmearFit[i][ptbin] = (TF1 *)gDirectory->Get(Form("FitForm_%i_%i", i, ptbin));
  //   }
  // }

  TFile* FONLL = new TFile("FONLL_Pt_5_20.root","READ");
  TH1F * FONLLCurve = (TH1F*) FONLL->Get("FONLL");

  TFile *PYTHIA = new TFile("PYTHIA_Pt_5_20.root", "READ");
  TH1F  *PYTHIAPtCurve = (TH1F *) PYTHIA->Get("PYTHIA pT");
  TH1F  *PYTHIAZCurve  = (TH1F *) PYTHIA->Get("PYTHIA Z");
  TH2F  *PYTHIA2D      = (TH2F *) PYTHIA->Get("PYTHIA");
  TH1F  *FONLLvPYTHIAWeights = (TH1F *) PYTHIA->Get("FONLLWeights");

  TFile* CentWeights = new TFile("Cent_Weight.root");
  TH1F * CENTWeight = (TH1F *)CentWeights->Get("Centrality Weights");

  ///// Everything we need to set up response matrices is imported before this line. /////

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

  RooUnfoldResponse *resp1D[3]; //Response
  RooUnfoldResponse *resp1DFONLL[3]; //Response Trying out the separation of centrality and fonll weights

  TH2D *fMeas2D[3];
  TH1D *fMeas1D[3];
  TH1D *fTrue1D[3];
  TH1D *fMiss1D[3];
  TH1D *fFake1D[3];

  TH1D *fMeasZ1D[3];
  TH1D *fTrueZ1D[3];
  TH1D *fMissZ1D[3];

  TH2D *hResp1D[3];
  TH2D *hMiss1D[3]; //These misses are because of eta-phi mismatch and not a pT mismatch
  TH2D *hFake1D[3]; //These fakes are because of eta-phi mismatch and not a pT mismatch

  TH2D *hResp1DFONLL[3];
  TH2D *hMiss1DFONLL[3]; //These misses are because of eta-phi mismatch and not a pT mismatch
  TH2D *hFake1DFONLL[3]; //These fakes are because of eta-phi mismatch and not a pT mismatch

  // TH2D *hRecoJetPtvsMCJetPt[3];
  TH1D *hDiffJetPt[3]; //Reco - MC
  TH1D *hDiffJetPtInGenPtBins[3][nbins_jpt]; //Reco - MC
  TH1D *hDiffJetEta[3]; //Reco - MC
  TH1D *hDiffJetPhi[3]; //Reco - MC

  TH1D *hMCD0Pt[3]; //Reco - MC
  TH1D *hRecoD0Pt[3]; //Reco - MC
  // TH1D *hDiffJetPhi[3]; //Reco - MC

  TH1D *hCentrality;

  TH1D *hRecoJetCorrPt[3];
  TH1D *hRecoZ[3];
  TH2D *hRecoJetCorrPtvZ[3];

  TH1D *hRecoJetCorrPt_Smeared[3];
  TH1D *hRecoZ_Smeared[3];
  TH2D *hRecoJetCorrPtvZ_Smeared[3]; // These histograms are generated by smearing MCJetPt with distributions from MC

  TH2D *OldRecoVsNewRecoPt[3];

  cout << JetTree->GetEntries() << endl;

  for (int i = 0; i < 3; i++){
    hZ[i] = new TH1D(Form("hZ_%i", i), Form("hZ_%i", i), nz_gen_bins, z_gen_low, z_gen_high);
    hJetPtvZ[i] = new TH2D(Form("hJetPtvZ_%i", i), Form("hJetPtvZ_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_low, z_gen_high);
    hJet[i] = new THnSparseF(Form("hJet_%i", i), Form("hJet_%i", i), ndim, nbins, xmin, xmax);
    hJetZNorm[i] = new THnSparseF(Form("hJetZNorm_%i", i), Form("hJetZNorm_%i", i), ndim, nbins, xmin, xmax);
    hJetZNormPtNorm[i] = new THnSparseF(Form("hJetZNormPtNorm_%i", i), Form("hJetZNormPtNorm_%i", i), ndim, nbins, xmin, xmax);

    fMeas[i] = new TH2D(Form("fMeas_Cent_%i", i), Form("fMeas_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);
    fTrue[i] = new TH2D(Form("fTrue_Cent_%i", i), Form("fTrue_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_low, z_gen_high);
    fMiss[i] = new TH2D(Form("fMiss_Cent_%i", i), Form("fMiss_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_low, z_gen_high);
    fFake[i] = new TH2D(Form("fFake_Cent_%i", i), Form("fFake_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);

    resp[i] = new RooUnfoldResponse(Form("Resp_%i", i), Form("Resp_%i", i));
    resp[i]->Setup(fMeas[i], fTrue[i]); //Setup Response Matrix Definition

    fMeas1D[i] = new TH1D(Form("fMeas1D_Cent_%i", i), Form("fMeas1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    fTrue1D[i] = new TH1D(Form("fTrue1D_Cent_%i", i), Form("fTrue1D_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin);
    fMiss1D[i] = new TH1D(Form("fMiss1D_Cent_%i", i), Form("fMiss1D_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin);
    fFake1D[i] = new TH1D(Form("fFake1D_Cent_%i", i), Form("fFake1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    
    fMeas2D[i] = new TH2D(Form("fMeas2D_Cent_%i", i), Form("fMeas2D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, 100, 0, 10);

    fMeasZ1D[i] = new TH1D(Form("fMeasZ1D_Cent_%i", i), Form("fMeasZ1D_Cent_%i", i), 100, 0, 10);
    fTrueZ1D[i] = new TH1D(Form("fTrueZ1D_Cent_%i", i), Form("fTrueZ1D_Cent_%i", i), 100, 0, 10);
    fMissZ1D[i] = new TH1D(Form("fMissZ1D_Cent_%i", i), Form("fMissZ1D_Cent_%i", i), 100, 0, 10);

    hDiffJetPt[i] = new TH1D(Form("hDiffJetPt_%i", i), Form("hDiffJetPt_%i", i), 100, -50, 50);
    hDiffJetEta[i] = new TH1D(Form("hDiffJetEta_%i", i), Form("hDiffJetEta_%i", i), 40, -1., 1.);
    hDiffJetPhi[i] = new TH1D(Form("hDiffJetPhi_%i", i), Form("hDiffJetPhi_%i", i), 40, -1., 1.);

    for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
      hDiffJetPtInGenPtBins[i][ptbin] = new TH1D(Form("hDiffJetPt_%i_%i", i, ptbin), Form("hDiffJetPt_%i_%i", i, ptbin), 100, -50, 50);
    }

    hMCD0Pt[i] = new TH1D(Form("hMCD0Pt_%i", i), Form("hMCD0Pt_%i", i), 20, 0, 10);
    hRecoD0Pt[i] = new TH1D(Form("hRecoD0Pt_%i", i), Form("hRecoD0Pt_%i", i), 20, 0, 10);

    hRecoJetCorrPt[i] = new TH1D(Form("hRecoJetCorrPt_%i", i), Form("hRecoJetCorrPt_%i", i), 100, -50, 50);
    hRecoZ[i] = new TH1D(Form("hRecoZ_%i", i), Form("hRecoZ_%i", i), 100, -5, 5);
    hRecoJetCorrPtvZ[i] = new TH2D(Form("hRecoJetCorrPtvZ_%i", i), Form("hRecoJetCorrPtvZ_%i", i), 100, -50, 50, 100, -5, 5);

    hRecoJetCorrPt_Smeared[i] = new TH1D(Form("hRecoJetCorrPt_Smeared_%i", i), Form("hRecoJetCorrPt_Smeared_%i", i), 100, -50, 50);
    hRecoZ_Smeared[i] = new TH1D(Form("hRecoZ_Smeared_%i", i), Form("hRecoZ_Smeared_%i", i), 100, -5, 5);
    hRecoJetCorrPtvZ_Smeared[i] = new TH2D(Form("hRecoJetCorrPtvZ_Smeared_%i", i), Form("hRecoJetCorrPtvZ_Smeared_%i", i), 100, -50, 50, 100, -5, 5);

    resp1D[i] = new RooUnfoldResponse(Form("Resp1D_%i", i), Form("Resp1D_%i", i));
    resp1D[i]->Setup(fMeas1D[i], fTrue1D[i]); //Setup Response Matrix Definition

    resp1DFONLL[i] = new RooUnfoldResponse(Form("Resp1D_New_%i", i), Form("Resp1D_New_%i", i));
    resp1DFONLL[i]->Setup(fMeas1D[i], fTrue1D[i]); //Setup Response Matrix Definition

    hResp1D[i] = new TH2D(Form("hResp1D_Cent_%i", i), Form("hResp1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins_var, jetpt_var_bin);
    hMiss1D[i] = new TH2D(Form("hMiss1D_Cent_%i", i), Form("hMiss1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins_var, jetpt_var_bin);
    hFake1D[i] = new TH2D(Form("hFake1D_Cent_%i", i), Form("hFake1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins_var, jetpt_var_bin);

    hResp1DFONLL[i] = new TH2D(Form("hResp1DFONLL_Cent_%i", i), Form("hResp1DFONLL_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins_var, jetpt_var_bin);
    hMiss1DFONLL[i] = new TH2D(Form("hMiss1DFONLL_Cent_%i", i), Form("hMiss1DFONLL_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins_var, jetpt_var_bin);
    hFake1DFONLL[i] = new TH2D(Form("hFake1DFONLL_Cent_%i", i), Form("hFake1DFONLL_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins_var, jetpt_var_bin);

    OldRecoVsNewRecoPt[i] = new TH2D(Form("OldRecoVsNewRecoPt_%i",i), Form("OldRecoVsNewRecoPt_%i",i), 100, -50, 50, 100, -50, 50);

    hJet[i]->Sumw2();
    hJetZNorm[i]->Sumw2();
    hJetZNormPtNorm[i]->Sumw2();

    hJet[i]->CalculateErrors();
    hJetZNorm[i]->CalculateErrors();
    hJetZNormPtNorm[i]->CalculateErrors();

    Weight[i] = new TH2D(Form("Weight_%i", i), Form("Weight_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_low, z_gen_high);

    if (!makeInitialDistributionDifferent){ // This is with PYTHIA pT and Z distribution in mind.
      for (int binx = 1; binx <= Weight[i]->GetNbinsX(); binx++ ){
        TH2D *tmp1 = (TH2D *)PYTHIA2D->Clone();
        tmp1->GetXaxis()->SetRange(binx, binx);
        TH1D *tmp2 = (TH1D *)tmp1->ProjectionY();
        tmp2->Scale(1./tmp2->Integral());
        for (int biny = 1; biny <= Weight[i]->GetNbinsY(); biny++ ){          
          Weight[i]->SetBinContent(binx, biny, tmp2->GetBinContent(biny));
        }
      }
    }
    else{ // Different Z Distribution to unfold
      for (int binx = 1; binx <= Weight[i]->GetNbinsX(); binx++ ){
        for (int biny = 1; biny <= Weight[i]->GetNbinsY(); biny++ ){
          // double w = 1. - pow(Weight[i]->GetYaxis()->GetBinCenter(biny), 2);
          // double w = 1. - pow(Weight[i]->GetYaxis()->GetBinCenter(biny) - 0.4, 2); // Moved the peak z value to 0.4 to see if it can be mimicked
          double w = 1. - pow(Weight[i]->GetYaxis()->GetBinCenter(biny) - 0.5, 3); // Moved the peak z value to 0.5 to see if it can be mimicked
          Weight[i]->SetBinContent(binx, biny, w);
        }
      }
    }
  }

  hCentrality = new TH1D("hCentrality","hCentrality", ncentbin, centbin);
  
  cout << JetTree->GetEntries() << endl;

  int nentries = JetTree->GetEntries();

  int lowlimit = 0;
  int highlimit = nentries;

  // nentries = 100;

  // if (SUPERITERATION == 0){lowlimit = 0; highlimit = nentries;}
  // else {lowlimit = nentries/2; highlimit = nentries;}

  double ptoffsetfactors[3]  = {0.89, 0.71, 0.20};
  double ptsmearfactors[3]  = {6.1184648, 4.4642124, 2.1123038};

  int counter = 0;

  double recooutside[3] = {0};

  TH1D *plottingpT = new TH1D("plottingpT", "plottingpT", nbins_jpt, binning_jpt);

  cout << "Mode == " << mode << endl;

  for (int i = lowlimit; i < highlimit; i++){
    // if (!dodata){
    //   if (SUPERITERATION == 0){if (i%2 == 0) continue;}
    //   else {if (i%2 != 0) continue;}
    // }

    // int midpoint = (highlimit - lowlimit)/2.;

    if (mode == 1) {if (i%2 != 0) continue;}
    else if (mode == 2) {if (i%2 == 0) continue;}

    if (mode == 1) {if (i%1000000 == 0) cout << "Read Entry " << i << endl;}
    else {if (i%1000000 == 1) cout << "Read Entry " << i << endl;}
    // cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

    if (centhistogramtofill < 0) continue;

    // if (MCPrimaryVertex->size() < 3 || RecoPrimaryVertex->size() < 3) continue;

    bool isMCD0Pt = MCD0Pt > 1 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 1 && RecoD0Pt < 10;
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);

    // if (mcz > 1.0) cout << MCJetPt << "\t" << MCD0Pt << mcz << endl;
    if (mcz == 1.0) mcz = 0.999; // Padding the boundaries

    bool isMCZ   = mcz > z_gen_low && mcz < z_gen_high; 

    if (!isMCZ) continue;

    // double smearfactors[3] = {5.8, 4.4, 1.875};
    // RecoJetCorrPt = RecoJetPt + gRandom->Gaus(0, smearfactors[centhistogramtofill]);

    double oldrecoz = (RecoD0Pt*TMath::Cos(RecoJetPhi - RecoD0Phi))/RecoJetCorrPt;

    hRecoJetCorrPt[centhistogramtofill]->Fill(RecoJetCorrPt);
    hRecoZ[centhistogramtofill]->Fill(oldrecoz);
    hRecoJetCorrPtvZ[centhistogramtofill]->Fill(RecoJetCorrPt, oldrecoz);

    int genptbin = plottingpT->FindBin(MCJetPt);

    double oldRecoJetCorrPt = RecoJetCorrPt;

    // RecoJetCorrPt = MCJetPt + PtSmearFit[centhistogramtofill][genptbin-1]->GetRandom(-20, 30);

    // if (centhistogramtofill == 0){if (RecoJetCorrPt - MCJetPt > 30) continue;}
    // if (centhistogramtofill == 1){if (RecoJetCorrPt - MCJetPt > 24) continue;}
    // if (centhistogramtofill == 2){if (RecoJetCorrPt - MCJetPt > 20) continue;}

    OldRecoVsNewRecoPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt, oldRecoJetCorrPt - MCJetPt);

    // RecoJetEta = MCJetEta + ((RecoJetEta - MCJetEta > 0) - (RecoJetEta - MCJetEta < 0))*abs(EtaSmear[centhistogramtofill]->GetRandom());
    // RecoJetPhi = MCJetPhi + GetSignOfDeltaPhi(RecoJetPhi, MCJetPhi)*PhiSmear[centhistogramtofill]->GetRandom();

    if(RecoJetPhi>=2.*TMath::Pi()) RecoJetPhi = RecoJetPhi - 2.*TMath::Pi();
    if(RecoJetPhi<0) RecoJetPhi = 2.*TMath::Pi() + RecoJetPhi;



    // if (centhistogramtofill == 0) RecoJetCorrPt = RecoJetCorrPt + (1.95-0.87);
    // if (centhistogramtofill == 1) RecoJetCorrPt = RecoJetCorrPt + (1.64-0.72);
    // if (centhistogramtofill == 2) RecoJetCorrPt = RecoJetCorrPt + (1.07-0.24);

    double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);

    hRecoJetCorrPt_Smeared[centhistogramtofill]->Fill(RecoJetCorrPt);
    hRecoZ_Smeared[centhistogramtofill]->Fill(recoz);
    hRecoJetCorrPtvZ_Smeared[centhistogramtofill]->Fill(RecoJetCorrPt, recoz);
    double RecoDeltaR = dR(dEta(RecoJetEta, RecoD0Eta), dPhi(RecoJetPhi, RecoD0Phi));

    // if (recoz > z_low-0.001 && recoz <= z_low) recoz = z_low + 0.001; // Padding the boundaries
    // if (recoz == 1.0) recoz = 0.99;

    //Defining response matrix using fakes and misses to check closure

    // bool isMCJetPt = MCJetPt > 5 && MCJetPt < 30;
    // bool isRecoJetPt = RecoJetCorrPt > -7 && RecoJetCorrPt < 30; // This is used because 99% of jets with pT > 5 GeV is captured within this range
    // bool isRecoJetPt = RecoJetCorrPt > -50 && RecoJetCorrPt < 50; // 78% of jets with pT > 5 GeV is captured within this range
    bool isRecoJetPt = RecoJetCorrPt > jetpt_low && RecoJetCorrPt < jetpt_high;
    bool isMCJetEta = abs(MCJetEta) < 0.601;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoDeltaR = RecoDeltaR < 0.4;

    // double rdca = TMath::Sqrt(pow(RecoPrimaryVertex->at(0) - MCPrimaryVertex->at(0), 2) + pow(RecoPrimaryVertex->at(1) - MCPrimaryVertex->at(1), 2) +  pow(RecoPrimaryVertex->at(2) - MCPrimaryVertex->at(2), 2));
    // if (rdca > 1.5) continue;
    // if (abs(RecoPrimaryVertex->at(2) - MCPrimaryVertex->at(2)) > 0.5) continue;

    // for (int i = 0; i < RecoPrimaryVertex->size(); i++){
    //   cout << RecoPrimaryVertex->at(i) << "\t" << MCPrimaryVertex->at(i) << endl;
    // }

    // cout << "Entry = " << i << "\t" << MCPrimaryVertex->at(0) << "\t" << MCPrimaryVertex->at(1) << "\t" << MCPrimaryVertex->at(2) << endl; 

    // bool isConstGreaterThanJet = (RecoJetCorrPt - RecoD0Pt > 0);

    // if (!isMCJetPt) continue;
    // if (!isRecoJetPt) continue;
    // if (!isMCJetEta) continue;
    // if (!isRecoJetEta) continue;
    // if (!isMCZ) continue;
    // double cutoff = 3.0;

    // if (RecoMaxTrackPt > cutoff) continue;
    // if (RecoMaxTowerEtAfterHC > cutoff) continue;

    int centbin = CENTWeight->FindBin(Centrality);

    double wcent = CENTWeight->GetBinContent(centbin);
    double fonllweight = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));
    // double wcent = 1.;

    // isRecoZ = kTRUE;
    // isMCZ = kTRUE;

    bool isOK   = isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta;
    bool isMISS = !isOK && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOK && isRecoJetPt && isRecoJetEta;

    // Ideally, we should decouple the centrality weights from the FONLL weights. I am not sure the decoupling changes anything, I don't know why it would, to be honest.
    // However, ideally this is the only thing left to check.

    double w = wcent*fonllweight;
    // double w = fonllweight;
    // double w = 1.;
    // if (abs(RecoJetCorrPt - MCJetPt) > 10) continue;

    if (!isOK && !isMISS && !isFAKE) continue;

    bool isRecoZ = recoz > z_low && recoz < z_high;

    if (isOK){
      fMeasZ1D[centhistogramtofill]->Fill(recoz, w);
      fMeas2D[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
    }

    // if (!isRecoZ) continue;

    hCentrality->Fill(Centrality, w);

    // if (abs(RecoJetPhi - MCJetPhi) > 0.6) cout << Form("MC Reco Delta = %.2f\t%.2f\t%.2f\t%.2f\n", MCJetPhi, RecoJetPhi, RecoJetPhi - MCJetPhi, dPhi(MCJetPhi, RecoJetPhi));

    hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt, wcent);
    hDiffJetEta[centhistogramtofill]->Fill(dEta(RecoJetEta, MCJetEta), wcent);
    hDiffJetPhi[centhistogramtofill]->Fill(dPhi(RecoJetPhi, MCJetPhi), wcent);

    hDiffJetPtInGenPtBins[centhistogramtofill][genptbin-1]->Fill(RecoJetCorrPt - MCJetPt, wcent);

    hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt, w);
    hRecoD0Pt[centhistogramtofill]->Fill(RecoD0Pt, w);

    double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
    hJet[centhistogramtofill]->Fill(fvalue, w);

    if (isOK)   hResp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, wcent);
    if (isMISS) hMiss1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, wcent);
    if (isFAKE) hFake1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, wcent); // This never happens in the current definition.

    if (isOK) {
      resp[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);

      fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);
    }
    else if (isMISS) {
      resp[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Miss(MCJetPt, w);

      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);
    }
    else if (isFAKE) {
      resp[centhistogramtofill]->Fake(RecoJetCorrPt, recoz, w);
      resp1D[centhistogramtofill]->Fake(RecoJetCorrPt, w);

      fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);

      fFake[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fFake1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);
    }
  }

  for (int cent = 0; cent < 3; cent++){
    cout << "Fakes = " << fFake[cent]->Integral() << endl; 
  }

  // FONLLising the weights
  
  for (int cent = 0; cent < 3; cent++){
    for (int binx = 1; binx <= hResp1D[cent]->GetNbinsX(); binx++){
      for (int biny = 1; biny <= hResp1D[cent]->GetNbinsY(); biny++){
        double binycenter = hResp1D[cent]->GetYaxis()->GetBinCenter(biny);
        double bincontent = hResp1D[cent]->GetBinContent(binx, biny);
        double binerror   = hResp1D[cent]->GetBinError(binx, biny);

        double fonllweight = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(binycenter));

        hResp1DFONLL[cent]->SetBinContent(binx, biny, bincontent*fonllweight);
        hResp1DFONLL[cent]->SetBinError(binx, biny, binerror*fonllweight);
      }
    }

    for (int binx = 0; binx <= hMiss1D[cent]->GetNbinsX()+1; binx++){
      for (int biny = 1; biny <= hMiss1D[cent]->GetNbinsY(); biny++){
        double binycenter = hMiss1D[cent]->GetYaxis()->GetBinCenter(biny);
        double bincontent = hMiss1D[cent]->GetBinContent(binx, biny);
        double binerror   = hMiss1D[cent]->GetBinError(binx, biny);

        double fonllweight = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(binycenter));

        hMiss1DFONLL[cent]->SetBinContent(binx, biny, bincontent*fonllweight);
        hMiss1DFONLL[cent]->SetBinError(binx, biny, binerror*fonllweight);
      }
    }
  }

  // Filling Response

  for (int cent = 0; cent < 3; cent++){
    for (int binx = 1; binx <= hResp1DFONLL[cent]->GetNbinsX(); binx++){
      for (int biny = 1; biny <= hResp1DFONLL[cent]->GetNbinsY(); biny++){
        double binxcenter = hResp1DFONLL[cent]->GetXaxis()->GetBinCenter(binx);
        double binycenter = hResp1DFONLL[cent]->GetYaxis()->GetBinCenter(biny);
        double bincontent = hResp1DFONLL[cent]->GetBinContent(binx, biny);
        double binerror   = hResp1DFONLL[cent]->GetBinError(binx, biny);

        resp1DFONLL[cent]->Fill(binxcenter, binycenter, bincontent);
      }
    }
  }

  // Filling Miss

  // // Testing to see if the missed distribution from the QM method remedies the mismatch between QM and HI

  // TFile *MissFileFromMC = new TFile("Jan26_FONLL_MC_CombinedFile4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root");
  // for (int cent = 0; cent < 3; cent++){
  //   hMiss1DFONLL[cent] = (TH2D *)MissFileFromMC->Get(Form("hMiss1D_Cent_%i", cent));
  // }

  for (int cent = 0; cent < 3; cent++){
    for (int biny = 1; biny <= hMiss1DFONLL[cent]->GetNbinsY(); biny++){
      double bincontent = 0.;
      double binycenter = hMiss1DFONLL[cent]->GetYaxis()->GetBinCenter(biny);

      for (int binx = 0; binx <= hMiss1D[cent]->GetNbinsX()+1; binx++){
        bincontent += hMiss1DFONLL[cent]->GetBinContent(binx, biny);
      }
        
      resp1DFONLL[cent]->Miss(binycenter, bincontent);
    }
  }

  cout << "Unaccepted = " << recooutside[0] << "\t" << recooutside[1] << "\t" << recooutside[2] << endl;

  TString filename;
  switch (mode){
    case 1:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins, mode);
    break;
    case 2:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    default:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins);
  }

  cout << "Centrality Integral = " << hCentrality->Integral() << endl;

  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->cd();

  for (int cent = 0;  cent < 3; cent++){
    hZ[cent]->Write();
    hJetPtvZ[cent]->Write();
    hJet[cent]->Write();
    hJetZNorm[cent]->Write();
    hJetZNormPtNorm[cent]->Write();

    fMeas[cent]->Write();
    fTrue[cent]->Write();
    fMiss[cent]->Write();
    fFake[cent]->Write();

    fMeas1D[cent]->Write();
    fTrue1D[cent]->Write();
    fMiss1D[cent]->Write();
    fFake1D[cent]->Write();

    fMeas2D[cent]->Write();
    fMeasZ1D[cent]->Write();
    fTrueZ1D[cent]->Write();
    fMissZ1D[cent]->Write();

    gDirectory->WriteObject(resp[cent], Form("Resp_%i", cent));
    gDirectory->WriteObject(resp1D[cent], Form("Resp1D_%i", cent));
    gDirectory->WriteObject(resp1DFONLL[cent], Form("Resp1D_FONLL_%i", cent));

    hDiffJetPt[cent]->Write();
    hDiffJetEta[cent]->Write();
    hDiffJetPhi[cent]->Write();

    hMCD0Pt[cent]->Write();
    hRecoD0Pt[cent]->Write();

    Weight[cent]->Write();

    hRecoJetCorrPt[cent]->Write();
    hRecoZ[cent]->Write();
    hRecoJetCorrPtvZ[cent]->Write();

    hRecoJetCorrPt_Smeared[cent]->Write();
    hRecoZ_Smeared[cent]->Write();
    hRecoJetCorrPtvZ_Smeared[cent]->Write();

    OldRecoVsNewRecoPt[cent]->Write();
  }

  for (int cent = 0; cent < 3; cent++){
    for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
      hDiffJetPtInGenPtBins[cent][ptbin]->Write();
    }
  }

  hCentrality->Write();
  // c->Write();
  // Ratio->Write();

  outfile->Close(); 

}

void PrepareResponseMatrixOldWayUsingHI_2D(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", bool makeInitialDistributionDifferent = kFALSE, int mode = 0){
  Method(SUPERITERATION, iteration, DirName, makeInitialDistributionDifferent, mode);
}