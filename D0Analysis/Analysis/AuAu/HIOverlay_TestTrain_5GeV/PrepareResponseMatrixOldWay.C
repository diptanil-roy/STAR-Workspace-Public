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

void Method(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", bool makeInitialDistributionDifferent = kFALSE, bool dodata = kFALSE){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);


  // TFile *f = new TFile("HI5GeV_WithMaxTrackTower.root");
  TFile *f = new TFile("Jan4_HIResponse.root");

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

  // TFile *SmearFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/BKGHistogramsFINALDec14_22.root");
  TFile *SmearFile = new TFile("Dec28_DataMC_PtChangeFromNewSP_OldCuts4/OldResponse_Step_1_IterParam_4_njptbin_15_nzbin_10.root");
  SmearFile->cd();
  TH1D *Smear[3];

  // Smear[0] = (TH1D *)gDirectory->Get("hBKG_0_10_Weighted");
  // Smear[1] = (TH1D *)gDirectory->Get("hBKG_10_40_Weighted");
  // Smear[2] = (TH1D *)gDirectory->Get("hBKG_40_80_Weighted");
  for (int i = 0; i < 3; i++){
    Smear[i] = (TH1D *)gDirectory->Get(Form("hDiffJetPt_%i", i));
    SetName(Smear[i], Form("Smear_%i", i));
  }

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

  RooUnfoldResponse *resp1D[3]; //Response
  TH1D *fMeas1D[3];
  TH1D *fTrue1D[3];
  TH1D *fMiss1D[3];

  // TH2D *hRecoJetPtvsMCJetPt[3];
  TH1D *hDiffJetPt[3]; //Reco - MC
  TH1D *hDiffJetEta[3]; //Reco - MC
  TH1D *hDiffJetPhi[3]; //Reco - MC

  TH1D *hMCD0Pt[3]; //Reco - MC
  TH1D *hRecoD0Pt[3]; //Reco - MC
  // TH1D *hDiffJetPhi[3]; //Reco - MC

  TH1D *hCentrality;

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

    resp[i] = new RooUnfoldResponse(Form("Resp_%i", i), Form("Resp_%i", i));
    resp[i]->Setup(fMeas[i], fTrue[i]); //Setup Response Matrix Definition

    fMeas1D[i] = new TH1D(Form("fMeas1D_Cent_%i", i), Form("fMeas1D_Cent_%i", i), njpt_bins, jetpt_low, jetpt_high);
    fTrue1D[i] = new TH1D(Form("fTrue1D_Cent_%i", i), Form("fTrue1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
    fMiss1D[i] = new TH1D(Form("fMiss1D_Cent_%i", i), Form("fMiss1D_Cent_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

    hDiffJetPt[i] = new TH1D(Form("hDiffJetPt_%i", i), Form("hDiffJetPt_%i", i), 100, -50, 50);
    hDiffJetEta[i] = new TH1D(Form("hDiffJetEta_%i", i), Form("hDiffJetEta_%i", i), 40, -1., 1.);
    hDiffJetPhi[i] = new TH1D(Form("hDiffJetPhi_%i", i), Form("hDiffJetPhi_%i", i), 40, -1., 1.);

    hMCD0Pt[i] = new TH1D(Form("hMCD0Pt_%i", i), Form("hMCD0Pt_%i", i), 20, 0, 10);
    hRecoD0Pt[i] = new TH1D(Form("hRecoD0Pt_%i", i), Form("hRecoD0Pt_%i", i), 20, 0, 10);

    resp1D[i] = new RooUnfoldResponse(Form("Resp1D_%i", i), Form("Resp1D_%i", i));
    resp1D[i]->Setup(fMeas1D[i], fTrue1D[i]); //Setup Response Matrix Definition

    hJet[i]->Sumw2();
    hJetZNorm[i]->Sumw2();
    hJetZNormPtNorm[i]->Sumw2();

    hJet[i]->CalculateErrors();
    hJetZNorm[i]->CalculateErrors();
    hJetZNormPtNorm[i]->CalculateErrors();

    Weight[i] = new TH2D(Form("Weight_%i", i), Form("Weight_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

    if (SUPERITERATION == 0){
      if (!makeInitialDistributionDifferent){ //Regular
        for (int binx = 1; binx <= Weight[i]->GetNbinsX(); binx++ ){
          for (int biny = 1; biny <= Weight[i]->GetNbinsY(); biny++ ){
            Weight[i]->SetBinContent(binx, biny, 0.);
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
    else if (SUPERITERATION == 1){
      for (int binx = 1; binx <= Weight[i]->GetNbinsX(); binx++ ){
        for (int biny = 1; biny <= Weight[i]->GetNbinsY(); biny++ ){
          Weight[i]->SetBinContent(binx, biny, 1.);
        }
      }
    }
  }

  hCentrality = new TH1D("hCentrality","hCentrality",20,0,100);
  
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

  TFile* FONLL = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/FONLLWeights_New.root","READ");
  TH1F * FONLLWeights = (TH1F*) FONLL->Get("FONLLWeights");

  TFile* CentWeights = new TFile("Cent_Weight.root");
  TH1F * CENTWeight = (TH1F *)CentWeights->Get("Centrality Weights");

  for (int i = lowlimit; i < highlimit; i++){
    if (!dodata){
      if (SUPERITERATION == 0){if (i%2 == 0) continue;}
      else {if (i%2 != 0) continue;}
    }

    if (SUPERITERATION > 0) {if (i%1000000 == 0) cout << "Read Entry " << i << endl;}
    else if (i%1000000 == 1) {cout << "Read Entry " << i << endl;}
    // cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

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
    if (mcz == 1.0) mcz = 0.999; // Padding the boundaries

    // double smearfactors[3] = {5.8, 4.4, 1.875};
    // RecoJetCorrPt = RecoJetPt + gRandom->Gaus(0, smearfactors[centhistogramtofill]);

    // RecoJetCorrPt = MCJetPt + Smear[centhistogramtofill]->GetRandom();

    // if (centhistogramtofill == 0) RecoJetCorrPt = RecoJetCorrPt + (1.95-0.87);
    // if (centhistogramtofill == 1) RecoJetCorrPt = RecoJetCorrPt + (1.64-0.72);
    // if (centhistogramtofill == 2) RecoJetCorrPt = RecoJetCorrPt + (1.07-0.24);

    double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);

    double RecoDeltaR = dR(dEta(RecoJetEta, RecoD0Eta), dPhi(RecoJetPhi, RecoD0Phi));

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
    // bool isRecoJetPt = RecoJetCorrPt > -50 && RecoJetCorrPt < 50; // 78% of jets with pT > 5 GeV is captured within this range
    bool isRecoJetPt = RecoJetCorrPt > 3 && RecoJetCorrPt < 30;
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoDeltaR = RecoDeltaR < 0.4;

    // bool isConstGreaterThanJet = (RecoJetCorrPt - RecoD0Pt > 0);

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (RecoJetNConst==0) continue;
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
    // double wcent = 1.;

    hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt, wcent);
    hDiffJetEta[centhistogramtofill]->Fill(RecoJetEta - MCJetEta, wcent);
    hDiffJetPhi[centhistogramtofill]->Fill(RecoJetPhi - MCJetPhi, wcent);

    hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt);
    hRecoD0Pt[centhistogramtofill]->Fill(RecoD0Pt);

    // isRecoZ = kTRUE;
    // isMCZ = kTRUE;

    bool isOK   = isRecoZ && isMCZ && isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta && isRecoDeltaR;
    bool isMISS = !isOK && isMCZ && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOK && isRecoZ && isRecoJetPt && isRecoJetEta && isRecoDeltaR;

    hCentrality->Fill(Centrality, wcent);

    // double w = FONLLWeights->GetBinContent(FONLLWeights->FindBin(MCJetPt))*wcent;
    double w = wcent;

    // if (abs(RecoJetCorrPt - MCJetPt) > 10) continue;

    if (isOK){
      resp[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);
    }
    else if (isMISS){
      resp[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      fMiss[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Miss(MCJetPt, w);
      fMiss1D[centhistogramtofill]->Fill(MCJetPt, w);
    }

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

  cout << "Unaccepted = " << recooutside[0] << "\t" << recooutside[1] << "\t" << recooutside[2] << endl;

  TString filename;
  filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins);
  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->cd();

  hCentrality->Write();

  for (int cent = 0;  cent < 3; cent++){
    hZ[cent]->Write();
    hJetPtvZ[cent]->Write();
    hJet[cent]->Write();

    fMeas[cent]->Write();
    fTrue[cent]->Write();
    fMiss[cent]->Write();

    fMeas1D[cent]->Write();

    gDirectory->WriteObject(resp[cent], Form("Resp_%i", cent));
    gDirectory->WriteObject(resp1D[cent], Form("Resp1D_%i", cent));

    hDiffJetPt[cent]->Write();
    hDiffJetEta[cent]->Write();
    hDiffJetPhi[cent]->Write();

    hMCD0Pt[cent]->Write();
    hRecoD0Pt[cent]->Write();
  }
  // Ratio->Write();

  outfile->Close(); 

}

void PrepareResponseMatrixOldWay(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", bool makeInitialDistributionDifferent = kFALSE, bool dodata = kFALSE){
  Method(SUPERITERATION, iteration, DirName, makeInitialDistributionDifferent, dodata);
}