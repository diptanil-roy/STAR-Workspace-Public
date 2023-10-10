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
// #include "StJetTreeStruct.h"
#include <vector>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

// #include "RooUnfold/src/RooUnfoldResponse.h"
// #include "RooUnfold/src/RooUnfoldBayes.h"
// #include "RooUnfold/src/RooUnfoldSvd.h"

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
#include "NewBinDef.h"

void Method(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", bool makeInitialDistributionDifferent = kFALSE, bool dodata = kFALSE){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  // TFile *f = new TFile("HIResponse.root");
  TFile *f = new TFile("HIResponse_MCIncluded_Feb6.root");

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
  JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
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

  // TFile *ZFile = new TFile("ZFile.root", "UPDATE");
  // TH1D *ZWanted = (TH1D *)ZFile->Get(Form("Z_%i_%i", SUPERITERATION, iteration));

  TH1D *hZ[3];
  TH2D *hJetPtvZ[3];
  THnSparseF *hJet[3];
  THnSparseF *hJetZNorm[3];
  THnSparseF *hJetZNormPtNorm[3];
  TH2D *Weight[3];
  cout << JetTree->GetEntries() << endl;

  for (int i = 0; i < 3; i++){
    hZ[i] = new TH1D(Form("hZ_%i", i), Form("hZ_%i", i), nz_gen_bins, z_gen_low, z_gen_high);
    hJetPtvZ[i] = new TH2D(Form("hJetPtvZ_%i", i), Form("hJetPtvZ_%i", i), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
    hJet[i] = new THnSparseF(Form("hJet_%i", i), Form("hJet_%i", i), ndim, nbins, xmin, xmax);
    hJetZNorm[i] = new THnSparseF(Form("hJetZNorm_%i", i), Form("hJetZNorm_%i", i), ndim, nbins, xmin, xmax);
    hJetZNormPtNorm[i] = new THnSparseF(Form("hJetZNormPtNorm_%i", i), Form("hJetZNormPtNorm_%i", i), ndim, nbins, xmin, xmax);

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
    
  if (SUPERITERATION > 1){
    TString laststepfilename;
    laststepfilename = Form("%s/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION-1, iteration, njpt_gen_bins, nz_gen_bins);
    cout << laststepfilename.Data() << endl;
    TFile *laststep = new TFile(laststepfilename.Data());
    laststep->cd();
    // laststep->ls();
    TH1D *ZUnf[3];
    TH1D *ZUnf_JetPt[3][njpt_gen_bins];

    for (int i = 0; i < 3; i++){
      cout << Weight[i]->GetName() << endl;
      ZUnf[i] = (TH1D *)gDirectory->Get(Form("UnfoldedZ_%i", i));
      // ZUnf[i]->SetDirectory(0);
      for (int jptbin = 1; jptbin <= njpt_gen_bins; jptbin++){
        ZUnf_JetPt[i][jptbin-1] = (TH1D *)gDirectory->Get(Form("ZUnf_%i_%i", i, jptbin));
        // cout << ZUnf_JetPt[i][jptbin-1]->Integral() << endl;
        ZUnf_JetPt[i][jptbin-1]->Scale(1./ZUnf_JetPt[i][jptbin-1]->Integral());
        for (int zbin = 1; zbin <= nz_gen_bins; zbin++){
          // cout << i << "\t" << jptbin << "\t" << zbin << "\t" << ZUnf_JetPt[i][jptbin-1]->GetBinContent(zbin) << endl;
          Weight[i]->SetBinContent(jptbin, zbin, ZUnf_JetPt[i][jptbin-1]->GetBinContent(zbin));
        }
        // ZUnf_JetPt[i][jptbin-1]->SetDirectory(0);
      }
    }
    // laststep->Close();
  } // I will include the definitions for superiteration > 1 from here. Need to do the full step once to know what histograms I need to import.


  cout << JetTree->GetEntries() << endl;

  int nentries = JetTree->GetEntries();

  int lowlimit = 0;
  int highlimit = nentries;

  // nentries = 100;

  // if (SUPERITERATION == 0){lowlimit = 0; highlimit = nentries;}
  // else {lowlimit = nentries/2; highlimit = nentries;}

  int counter = 0;

  double recooutside[3] = {0};

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

    double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);

    if (recoz > z_low-0.001 && recoz <= z_low) recoz = z_low + 0.001; // Padding the boundaries
    // if (recoz == 1.0) recoz = 0.99;

    //Defining response matrix using fakes and misses to check closure

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 5 && RecoD0Pt < 10;
    bool isMCZ   = mcz > z_gen_low && mcz < z_gen_high; 
    bool isRecoZ = recoz > z_low && recoz < z_high; 
    // bool isMCJetPt = MCJetPt > 5 && MCJetPt < 30;
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    // bool isRecoJetPt = RecoJetCorrPt > -7 && RecoJetCorrPt < 30; // This is used because 99% of jets with pT > 5 GeV is captured within this range
    bool isRecoJetPt = RecoJetCorrPt > 3 && RecoJetCorrPt < 30; // 78% of jets with pT > 5 GeV is captured within this range
    bool isMCJetEta = abs(MCJetEta) < 0.601;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    if (RecoJetNConst==0) continue;

    // bool isConstGreaterThanJet = (RecoJetCorrPt - RecoD0Pt > 0);

    // if (!isMCJetPt) continue;
    // // if (!isRecoJetPt) continue;
    // if (!isMCJetEta) continue;
    // if (!isRecoJetEta) continue;
    // if (!isMCZ) continue;

    bool isOK   = isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta;
    bool isMISS = !isOK && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOK && isRecoJetPt && isRecoJetEta;

    if (!isOK && !isMISS && !isFAKE) continue;

    // if (!isConstGreaterThanJet) continue;
    // if (MCJetNConst <= 2 ) continue;
    // if (mcz > 0.88) cout << MCJetPt << "\t" << MCD0Pt << "\t" << MCJetNConst << "\t" << mcz << endl;

    if (!isRecoJetPt || !isRecoZ) recooutside[centhistogramtofill]++;


    double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
    hJet[centhistogramtofill]->Fill(fvalue);
    hZ[centhistogramtofill]->Fill(mcz);
    hJetPtvZ[centhistogramtofill]->Fill(MCJetPt, mcz);
  }

  cout << "Unaccepted = " << recooutside[0] << "\t" << recooutside[1] << "\t" << recooutside[2] << endl;

  TH1D *TruthZ[3];
  TH1D *TruthPt[3];

  for (int cent = 0; cent < 3; cent++){
    TruthZ[cent] = (TH1D *)hJet[cent]->Projection(1, "E");
    TruthPt[cent] = (TH1D *)hJet[cent]->Projection(0, "E");
  }

  // TH2D *PythiaDist[3];
  // TCanvas *c[3];

  // for (int cent = 0; cent < 3; cent++){
  //   PythiaDist[cent] = (TH2D *)hJet[cent]->Projection(1,0, "E");
  //   Weight[cent]->Divide(PythiaDist[cent]);
  //   c[cent] = new TCanvas(Form("c_%i", cent), Form("c_%i", cent), 800, 800);
  //   c[cent]->cd();
  //   Weight[cent]->Draw("COLZ");
  // }

  // for (int cent = 0; cent < 3; cent++){
  //   hJetZNorm[cent] = (THnD *)hJet[cent]->Clone(Form("hJetZNorm_%i", cent));
  //   for (int gbinx = 0; gbinx <= njpt_gen_bins; gbinx++ ){
  //     for (int gbiny = 0; gbiny <= nz_gen_bins; gbiny++ ){
  //       for (int rbinx = 0; rbinx <= njpt_bins; rbinx++ ){
  //         for (int rbiny = 0; rbiny <= nz_bins; rbiny++ ){
  //           int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
  //           int binZnorm = hJetZNorm[cent]->GetBin(xbin);
  //           hJetZNorm[cent]->SetBinContent(binZnorm, 0);
  //           hJetZNorm[cent]->SetBinError(binZnorm, 0);
  //         }
  //       }
  //     }
  //   }
  // }



  // Changing the z distribution here
  for (int cent = 0; cent < 3; cent++){
    for (int gbinx = 0; gbinx <= njpt_gen_bins; gbinx++ ){
      for (int gbiny = 0; gbiny <= nz_gen_bins; gbiny++ ){
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
          if (SUPERITERATION == 0 && !makeInitialDistributionDifferent) scalefactor = 1;
          else scalefactor = Weight[cent]->GetBinContent(gbinx, gbiny)/valintegral;
        }
        Weight[cent]->SetBinContent(gbinx, gbiny, scalefactor);

        for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
          for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){

            int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
            int binZnorm = hJetZNorm[cent]->GetBin(xbin);
            int bin = hJet[cent]->GetBin(xbin);
            
            double content = hJetZNorm[cent]->GetBinContent(binZnorm);
            double err = hJetZNorm[cent]->GetBinError(binZnorm);

            double newbincontent = content + hJet[cent]->GetBinContent(bin)*scalefactor;
            double newbinerr = err + hJet[cent]->GetBinError(bin)*scalefactor;
            hJetZNorm[cent]->SetBinContent(binZnorm, newbincontent);
            hJetZNorm[cent]->SetBinError(binZnorm, newbinerr);

            // double mid[ndim];
            // mid[0] = hJetZNorm[cent]->GetAxis(0)->GetBinCenter(gbinx);
            // mid[1] = hJetZNorm[cent]->GetAxis(1)->GetBinCenter(gbiny);
            // mid[2] = hJetZNorm[cent]->GetAxis(2)->GetBinCenter(rbinx);
            // mid[3] = hJetZNorm[cent]->GetAxis(3)->GetBinCenter(rbiny);

            // hJetZNorm[cent]->Fill(mid, hJet[cent]->GetBinContent(bin)*scalefactor);
            // hJetZNorm[cent]->SetBinError(binZnorm, newbinerr);
          }
        }
      }
    }
  }

  // cout << "Changed Z Distribution" << endl;

  TCanvas *c[3];
  for (int cent = 0; cent < 3; cent++){
    c[cent] = new TCanvas(Form("c_%i", cent), Form("c_%i", cent), 800, 800);
    c[cent]->cd();
    // Weight[cent]->Draw("COLZ");
    hJetPtvZ[cent]->Draw("COLZ");
  }

  // cout << "Changed Z Distribution" << endl;

  TH1D *TruthZZNorm[3];
  TH1D *TruthPtZNorm[3];

  for (int cent = 0; cent < 3; cent++){
    TruthZZNorm[cent] = (TH1D *)hJetZNorm[cent]->Projection(1, "E");
    TruthPtZNorm[cent] = (TH1D *)hJetZNorm[cent]->Projection(0, "E");
  }

  double ent[3][njpt_gen_bins+2] = {0};

  // for (int cent = 0; cent < 3; cent++){
  //   hJetZNormPtNorm[cent] = (THnD *)hJet[cent]->Clone(Form("hJetZNormPtNorm_%i", cent));
  //   for (int gbinx = 0; gbinx <= njpt_gen_bins; gbinx++ ){
  //     for (int gbiny = 0; gbiny <= nz_gen_bins; gbiny++ ){
  //       for (int rbinx = 0; rbinx <= njpt_bins; rbinx++ ){
  //         for (int rbiny = 0; rbiny <= nz_bins; rbiny++ ){
  //           int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
  //           int binZnorm = hJetZNormPtNorm[cent]->GetBin(xbin);
  //           hJetZNormPtNorm[cent]->SetBinContent(binZnorm, 0);
  //           hJetZNormPtNorm[cent]->SetBinError(binZnorm, 0);
  //         }
  //       }
  //     }
  //   }
  // } 

  // Changing the pT distribution to go back to the original distribution
  for (int cent = 0; cent < 3; cent++){
    for (int gbinx = 0; gbinx <= njpt_gen_bins; gbinx++ ){
      double scalefactor = 0;
      double origIntegral = 0;
      double scalIntegral = 0;
      for (int gbiny = 0; gbiny <= nz_gen_bins; gbiny++ ){
        for (int rbinx = 0; rbinx <= njpt_bins+1; rbinx++ ){
          for (int rbiny = 0; rbiny <= nz_bins+1; rbiny++ ){

            int xbin[ndim] = {gbinx, gbiny, rbinx, rbiny};
            int binZnorm = hJetZNorm[cent]->GetBin(xbin);
            int bin = hJet[cent]->GetBin(xbin);

            origIntegral+=hJet[cent]->GetBinContent(bin);
            scalIntegral+=hJetZNorm[cent]->GetBinContent(bin);
            
          }
        }
      }
      if (scalIntegral > 0) {
        if (SUPERITERATION == 0 && !makeInitialDistributionDifferent) scalefactor = 1.;
        else scalefactor = origIntegral/scalIntegral;
      }

      for (int gbiny = 0; gbiny <= nz_gen_bins; gbiny++ ){
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

            // double mid[ndim];
            // mid[0] = hJetZNormPtNorm[cent]->GetAxis(0)->GetBinCenter(gbinx);
            // mid[1] = hJetZNormPtNorm[cent]->GetAxis(1)->GetBinCenter(gbiny);
            // mid[2] = hJetZNormPtNorm[cent]->GetAxis(2)->GetBinCenter(rbinx);
            // mid[3] = hJetZNormPtNorm[cent]->GetAxis(3)->GetBinCenter(rbiny);

            // hJetZNormPtNorm[cent]->Fill(mid, hJetZNorm[cent]->GetBinContent(binZnorm)*scalefactor);

            ent[cent][gbinx]+=hJetZNorm[cent]->GetBinContent(binZnorm)*scalefactor;
            // hJetZNormPtNorm[cent]->SetBinError(binZnormPtnorm, newbinerr);
          }
        }
      }
    }
  } //End of centrality loop

  // for (int cent = 0; cent < 3; cent++){
  //   cout << hJet[cent]->GetEntries() << "\t" << hJetZNorm[cent]->GetEntries() << "\t" << hJetZNormPtNorm[cent]->GetEntries() << endl;
  // }

  // for (int cent = 0; cent < 3; cent++){
  //   for (int gbinx = 0; gbinx <= njpt_gen_bins; gbinx++ ){
  //     cout << "Cent " << cent << "Bin X " << gbinx << ": " << ent[cent][gbinx] << endl;
  //   }
  // } 

  TH1D *TruthZZNormPtNorm[3];
  TH1D *TruthPtZNormPtNorm[3];

  for (int cent = 0; cent < 3; cent++){
    TruthZZNormPtNorm[cent] = (TH1D *)hJetZNormPtNorm[cent]->Projection(1, "E");
    TruthPtZNormPtNorm[cent] = (TH1D *)hJetZNormPtNorm[cent]->Projection(0, "E");

    cout << TruthPtZNormPtNorm[cent]->Integral() << endl;
  }

  for (int cent = 0; cent < 3; cent++){
    cout << hJet[cent]->GetEntries() << "\t" << hJetZNorm[cent]->GetEntries() << "\t" << hJetZNormPtNorm[cent]->GetEntries() << endl;
  }

  TCanvas *d[3];

  for (int cent = 0; cent < 3; cent++){
    d[cent] = new TCanvas(Form("Testing the weighting mechanism %i", cent), Form("Testing the weighting mechanism %i", cent), 1200, 600);
    d[cent]->Divide(2);
    d[cent]->cd(1);
    gPad->SetLogy();
    SetAxisTitles(TruthPt[cent], "p_{T,Jet} [GeV/#it{c}]", "Normalised Yield");
    SetColor(TruthPt[cent], kBlue, 24);
    SetName(TruthPt[cent], "PYTHIA");
    SetColor(TruthPtZNorm[cent], kRed, 24);
    SetName(TruthPtZNorm[cent], "Z norm (2D)");
    SetColor(TruthPtZNormPtNorm[cent], kBlack, 20);
    SetName(TruthPtZNormPtNorm[cent], "Z norm + PYTHIA pT");
    double scale = TruthPt[cent]->Integral();
    TruthPt[cent]->Scale(1./TruthPt[cent]->Integral());
    TruthPtZNorm[cent]->Scale(1./TruthPtZNorm[cent]->Integral());
    TruthPtZNormPtNorm[cent]->Scale(1./TruthPtZNormPtNorm[cent]->Integral());
    TruthPt[cent]->Draw();
    TruthPtZNorm[cent]->Draw("SAME");
    TruthPtZNormPtNorm[cent]->Draw("SAME");
    gPad->BuildLegend();
    // d[cent]->Update();

    d[cent]->cd(2);
    gPad->SetLogy();
    SetAxisTitles(TruthZ[cent], "Z", "Normalised Yield");
    SetColor(TruthZ[cent], kBlue, 24);
    SetName(TruthZ[cent], "PYTHIA");
    SetColor(TruthZZNorm[cent], kRed, 24);
    SetName(TruthZZNorm[cent], "Z norm (2D)");
    SetColor(TruthZZNormPtNorm[cent], kBlack, 20);
    SetName(TruthZZNormPtNorm[cent], "Z norm + PYTHIA pT");
    scale = TruthZ[cent]->Integral();
    TruthZ[cent]->Scale(1./TruthZ[cent]->Integral());
    TruthZZNorm[cent]->Scale(1./TruthZZNorm[cent]->Integral());
    TruthZZNormPtNorm[cent]->Scale(1./TruthZZNormPtNorm[cent]->Integral());
    TruthZ[cent]->Draw();
    TruthZZNorm[cent]->Draw("SAME");
    TruthZZNormPtNorm[cent]->Draw("SAME");
    gPad->BuildLegend();

    // d[cent]->Update();
  }

  TCanvas *Ratio = new TCanvas("Ratio", "Ratio", 1200, 600);
  Ratio->Divide(2);
  
  TH1D *rTruthPt[3]; //Ratio
  TH1D *rTruthZ[3]; //Ratio

  for (int i = 0; i < 3; i++){
    rTruthPt[i] = (TH1D *)TruthPtZNormPtNorm[i]->Clone();
    rTruthPt[i]->Divide(TruthPt[i]);
    SetName(rTruthPt[i], Form("p_{T} Ratio %i", i));
    SetColor(rTruthPt[i], col[i*2]);

    rTruthZ[i] = (TH1D *)TruthZZNormPtNorm[i]->Clone();
    rTruthZ[i]->Divide(TruthZ[i]);
    SetName(rTruthZ[i], Form("Z Ratio %i", i));
    SetColor(rTruthZ[i], col[i*2]);
  }

  Ratio->cd(1);
  for (int i = 0; i < 3; i++){
    rTruthPt[i]->GetYaxis()->SetRangeUser(0, 2);
    rTruthPt[i]->Draw("SAME");
  }
  gPad->BuildLegend();

  Ratio->cd(2);
  for (int i = 0; i < 3; i++){
    rTruthZ[i]->Draw("SAME");
  }
  gPad->BuildLegend();


  TH2D *InitialDist[3];
  TH2D *FinalDist[3];
  TCanvas *e[3];

  for (int cent = 0; cent < 3; cent++){
    InitialDist[cent] = (TH2D *)hJet[cent]->Projection(1,0, "E");  
    FinalDist[cent] = (TH2D *)hJetZNormPtNorm[cent]->Projection(1,0, "E");
    e[cent] = new TCanvas(Form("e_%i", cent), Form("e_%i", cent), 800, 800);
    e[cent]->cd();

    FinalDist[cent]->Scale(1./FinalDist[cent]->Integral());
    InitialDist[cent]->Scale(1./InitialDist[cent]->Integral());
    FinalDist[cent]->Divide(InitialDist[cent]);
    SetName(FinalDist[cent], Form("Final Weights %i", cent));
    FinalDist[cent]->Draw("COLZ");
  }



  TString filename;
  filename = Form("%s/MCMCResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), SUPERITERATION, iteration, njpt_gen_bins, nz_gen_bins);
  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->cd();

  for (int cent = 0;  cent < 3; cent++){
    hZ[cent]->Write();
    hJetPtvZ[cent]->Write();
    hJet[cent]->Write();
    hJetZNorm[cent]->Write();
    hJetZNormPtNorm[cent]->Write();

    c[cent]->Write();
    d[cent]->Write();
    e[cent]->Write();
  }
  Ratio->Write();

  outfile->Close(); 

}

void PrepareResponseMatrix(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", bool makeInitialDistributionDifferent = kFALSE, bool dodata = kFALSE){
  Method(SUPERITERATION, iteration, DirName, makeInitialDistributionDifferent, dodata);
}