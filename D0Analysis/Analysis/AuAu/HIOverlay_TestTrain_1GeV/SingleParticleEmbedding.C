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

void Method(TString DirName = "SP"){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TFile *f = new TFile("/Volumes/WorkDrive/work/2022/Response2022/Embed_April2.root");

  cout << f->GetName() <<  endl;
  f->cd("JetTree_Embed");

  TTree *JetTree = (TTree *)gDirectory->Get("D0Jets");
  // TTree *RecoJetTree = (TTree *)gDirectory->Get("RecoJets");

  Float_t         Centrality;
  vector<double>  *MCPrimaryVertex = new vector<double>;
  vector<double>  *RecoPrimaryVertex = new vector<double>;
  Float_t         Weight;
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

  JetTree->SetBranchStatus("*", 0);
  JetTree->SetBranchStatus("Centrality", 1);
  JetTree->SetBranchStatus("Weight", 1);
  JetTree->SetBranchStatus("JetCorrPt", 1);
  JetTree->SetBranchStatus("JetEta", 1);
  JetTree->SetBranchStatus("JetPhi", 1);
  JetTree->SetBranchStatus("PionPt", 1);
  JetTree->SetBranchStatus("PionEta", 1);
  JetTree->SetBranchStatus("PionPhi", 1);

  JetTree->SetBranchAddress("Centrality", &Centrality);
  JetTree->SetBranchAddress("Weight", &Weight);
  JetTree->SetBranchAddress("JetCorrPt", &RecoJetCorrPt);
  JetTree->SetBranchAddress("JetEta", &RecoJetEta);
  JetTree->SetBranchAddress("JetPhi", &RecoJetPhi);
  JetTree->SetBranchAddress("PionPt", &MCJetPt);
  JetTree->SetBranchAddress("PionEta", &MCJetEta);
  JetTree->SetBranchAddress("PionPhi", &MCJetPhi);

  int nentries = JetTree->GetEntries();

  TH1D *hDiffJetPt[3]; //Reco - MC
  TH1D *hDiffJetEta[3]; //Reco - MC
  TH1D *hDiffJetPhi[3]; //Reco - MC

  TH1D *hDiffJetPtWeighed[3]; //Reco - MC
  TH1D *hDiffJetEtaWeighed[3]; //Reco - MC
  TH1D *hDiffJetPhiWeighed[3]; //Reco - MC

  TH1D *hMCD0Pt[3]; //Reco - MC
  TH1D *hRecoD0Pt[3]; //Reco - MC

  TH1D *hMCJetPt[3];
  TH1D *hMCJetPtWeighed[3];

  TH1F *hDiffJetPtNew[9];
  TH1F *hDiffJetPtNewWeighed[9];

  cout << JetTree->GetEntries() << endl;

  for (int i = 0; i < 3; i++){
    hDiffJetPt[i] = new TH1D(Form("hDiffJetPt_%i", i), Form("hDiffJetPt_%i", i), 100, -50, 50);
    hDiffJetEta[i] = new TH1D(Form("hDiffJetEta_%i", i), Form("hDiffJetEta_%i", i), 40, -1., 1.);
    hDiffJetPhi[i] = new TH1D(Form("hDiffJetPhi_%i", i), Form("hDiffJetPhi_%i", i), 40, -1., 1.);

    hDiffJetPtWeighed[i] = new TH1D(Form("hDiffJetPtWeighed_%i", i), Form("hDiffJetPtWeighed_%i", i), 100, -50, 50);
    hDiffJetEtaWeighed[i] = new TH1D(Form("hDiffJetEtaWeighed_%i", i), Form("hDiffJetEtaWeighed_%i", i), 40, -1., 1.);
    hDiffJetPhiWeighed[i] = new TH1D(Form("hDiffJetPhiWeighed_%i", i), Form("hDiffJetPhiWeighed_%i", i), 40, -1., 1.);

    hMCD0Pt[i] = new TH1D(Form("hMCD0Pt_%i", i), Form("hMCD0Pt_%i", i), 20, 0, 10);
    hRecoD0Pt[i] = new TH1D(Form("hRecoD0Pt_%i", i), Form("hRecoD0Pt_%i", i), 20, 0, 10);

    hMCJetPt[i] = new TH1D(Form("hMCJetPt_%i", i), Form("hMCJetPt_%i", i), 80, 0, 40);
    hMCJetPtWeighed[i] = new TH1D(Form("hMCJetPtWeighed_%i", i), Form("hMCJetPtWeighed_%i", i), 80, 0, 40);
  }

  for (int i = 0; i < 9; i++){
    hDiffJetPtNew[i] = new TH1F(Form("hDiffJetPtNew_%i", i), Form("hDiffJetPtNew_%i", i), 100, -50, 50);
    hDiffJetPtNewWeighed[i] = new TH1F(Form("hDiffJetPtNewWeighed_%i", i), Form("hDiffJetPtNewWeighed_%i", i), 100, -50, 50);
  }

  for (int i = 0; i < nentries; i++){
    JetTree->GetEntry(i);
    if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    if (Centrality < 0) continue;

    int centhistogramtofill = -99;
    if (Centrality > 6 && Centrality <= 8) centhistogramtofill = 0;
    else if (Centrality > 3 && Centrality <= 6) centhistogramtofill = 1;
    else if (Centrality >= 0 && Centrality <= 3) centhistogramtofill = 2;


    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;

    if (!isRecoJetEta) continue;

    Weight = 1.;

    hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt, Weight);
    hDiffJetEta[centhistogramtofill]->Fill(RecoJetEta - MCJetEta, Weight);
    hDiffJetPhi[centhistogramtofill]->Fill(RecoJetPhi - MCJetPhi, Weight);

    hDiffJetPtWeighed[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt, Weight*1.0/(MCJetPt));
    hDiffJetEtaWeighed[centhistogramtofill]->Fill(RecoJetEta - MCJetEta, Weight*1.0/(MCJetPt));
    hDiffJetPhiWeighed[centhistogramtofill]->Fill(RecoJetPhi - MCJetPhi, Weight*1.0/(MCJetPt));

    hDiffJetPtNew[(int)Centrality]->Fill(RecoJetCorrPt - MCJetPt, Weight);
    hDiffJetPtNewWeighed[(int)Centrality]->Fill(RecoJetCorrPt - MCJetPt, Weight*1.0/(MCJetPt));

    hMCJetPt[centhistogramtofill]->Fill(MCJetPt, Weight);
    hMCJetPtWeighed[centhistogramtofill]->Fill(MCJetPt, Weight*1.0/(MCJetPt));
  }

  TFile* fin = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/ApplyWeights/Histograms3_D05GeV.root");
  TH1F*centrality = (TH1F*)fin->Get("centrality");
  double val9 = centrality->GetBinContent(centrality->FindBin(70));
  double val8 = centrality->GetBinContent(centrality->FindBin(60));
  double val7 = centrality->GetBinContent(centrality->FindBin(50));
  double val6 = centrality->GetBinContent(centrality->FindBin(40));
  double val5 = centrality->GetBinContent(centrality->FindBin(30));
  double val4 = centrality->GetBinContent(centrality->FindBin(20));
  double val3 = centrality->GetBinContent(centrality->FindBin(10));
  double val2 = centrality->GetBinContent(centrality->FindBin(5));
  double val1 = centrality->GetBinContent(centrality->FindBin(0));    

  double sval9 = hDiffJetPtNew[0]->Integral(); // 70-80
  double sval8 = hDiffJetPtNew[1]->Integral();
  double sval7 = hDiffJetPtNew[2]->Integral();
  double sval6 = hDiffJetPtNew[3]->Integral();
  double sval5 = hDiffJetPtNew[4]->Integral();
  double sval4 = hDiffJetPtNew[5]->Integral();
  double sval3 = hDiffJetPtNew[6]->Integral();
  double sval2 = hDiffJetPtNew[7]->Integral();
  double sval1 = hDiffJetPtNew[8]->Integral(); // 0-5



  cout << "w9 " << val9<< endl;
  cout << "w8 " << val8<< endl;
  cout << "w7 " << val7<< endl;
  cout << "w6 " << val6<< endl;
  cout << "w5 " << val5<< endl;
  cout << "w4 " << val4<< endl;
  cout << "w3 " << val3<< endl;
  cout << "w2 " << val2<< endl;
  cout << "w1 " << val1<< endl;




  double norm = val9+val8+val7+val6;
  double snorm = sval9+sval8+sval7+sval6;
  double w9 = (val9/norm)/(sval9/snorm);
  double w8 = (val8/norm)/(sval8/snorm);
  double w7 = (val7/norm)/(sval7/snorm);
  double w6 = (val6/norm)/(sval6/snorm);
 

  norm= val3+val4+val5;
  snorm = sval3+sval4+sval5;
  double w5 = (val5/norm)/(sval5/snorm);
  double w4 = (val4/norm)/(sval4/snorm);
  double w3 = (val3/norm)/(sval3/snorm);
  
  norm= val1+val2;
  snorm = sval1+sval2;
  double w2 = (val2/norm)/(sval2/snorm);
  double w1 = (val1/norm)/(sval1/snorm);

  cout << "w9 " << w9<< endl;
  cout << "w8 " << w8<< endl;
  cout << "w7 " << w7<< endl;
  cout << "w6 " << w6<< endl;
  cout << "w5 " << w5<< endl;
  cout << "w4 " << w4<< endl;
  cout << "w3 " << w3<< endl;
  cout << "w2 " << w2<< endl;
  cout << "w1 " << w1<< endl;

  TH1F *hDiffJetPtScaled[3];

  hDiffJetPtScaled[0] = (TH1F *)hDiffJetPtNew[8]->Clone("hDiffJetPtScaled_0");
  hDiffJetPtScaled[0]->Scale(w1);
  hDiffJetPtScaled[0]->Add(hDiffJetPtNew[7], w2);

  hDiffJetPtScaled[1] = (TH1F *)hDiffJetPtNew[6]->Clone("hDiffJetPtScaled_1");
  hDiffJetPtScaled[1]->Scale(w3);
  hDiffJetPtScaled[1]->Add(hDiffJetPtNew[5], w4);
  hDiffJetPtScaled[1]->Add(hDiffJetPtNew[4], w5);

  hDiffJetPtScaled[2] = (TH1F *)hDiffJetPtNew[3]->Clone("hDiffJetPtScaled_2");
  hDiffJetPtScaled[2]->Scale(w6);
  hDiffJetPtScaled[2]->Add(hDiffJetPtNew[2], w7);
  hDiffJetPtScaled[2]->Add(hDiffJetPtNew[1], w8);
  hDiffJetPtScaled[2]->Add(hDiffJetPtNew[0], w9);


  TString filename;
  filename = Form("%s/SingleParticleEmbedding.root", DirName.Data());
  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->cd();

  for (int cent = 0;  cent < 3; cent++){
    hDiffJetPt[cent]->Write();
    hDiffJetEta[cent]->Write();
    hDiffJetPhi[cent]->Write();

    hDiffJetPtWeighed[cent]->Write();
    hDiffJetEtaWeighed[cent]->Write();
    hDiffJetPhiWeighed[cent]->Write();

    hMCJetPt[cent]->Write();
    hMCJetPtWeighed[cent]->Write();

    hDiffJetPtScaled[cent]->Write();
    // hMCD0Pt[cent]->Write();
    // hRecoD0Pt[cent]->Write();
  }
  // Ratio->Write();

  outfile->Close();

  TCanvas *c = new TCanvas("c", "c", 2000, 600);
  c->Divide(3);

  for (int cent = 0;  cent < 3; cent++){
    c->cd(cent + 1);
    gPad->SetLogy();
    hMCJetPtWeighed[cent]->Draw("HIST SAME");
    hMCJetPtWeighed[cent]->Scale(1./hMCJetPtWeighed[cent]->Integral());
    SetColor(hMCJetPtWeighed[cent], kRed);
    hMCJetPt[cent]->Draw("SAME");
    hMCJetPt[cent]->Scale(1./hMCJetPt[cent]->Integral());
  }

  TCanvas *d = new TCanvas("d", "d", 2000, 600);
  d->Divide(3);

  for (int cent = 0;  cent < 3; cent++){
    d->cd(cent + 1);
    gPad->SetLogy();
    hDiffJetPt[cent]->Draw();
    // hDiffJetPt[cent]->Scale(1./hDiffJetPt[cent]->Integral());
    hDiffJetPtScaled[cent]->Draw("HIST SAME");
    // hDiffJetPtScaled[cent]->Scale(1./hDiffJetPtScaled[cent]->Integral());
    SetColor(hDiffJetPtScaled[cent], kRed);
  }

}

void SingleParticleEmbedding(TString DirName){
  Method(DirName);
}