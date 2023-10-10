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

#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"

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

void Method(int pthatbin = 0, int SUPERITERATION = 0, int iteration = 3){

  cout << "File Reader Method" << endl;

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

  // TFile *ZFile = new TFile("ZFile.root", "UPDATE");
  // TH1D *ZWanted = (TH1D *)ZFile->Get(Form("Z_%i_%i", SUPERITERATION, iteration));

  // Calling all Histograms that will be saved out

  TH1D *NEvents;
  TH1D *NEventsWeighted;
  TH1D *hVzMC[3];

  int nBinsZHist[nDimZHist] = {nBinsJetPtForZHist, nBinsForZHist, nBinsJetPtForZHist, nBinsForZHist};


  RooUnfoldResponse *resp2[3];
  RooUnfoldResponse *newresp[3];
  RooUnfoldResponse *newresp2[3];

  TH2D *tmpMC = new TH2D("tmpMC", "tmpMC", nBinsJetPtForMCZHist, JetPtBinsForMCZHist, nBinsForMCZHist, ZBinsForMCZHist);
  TH2D *tmpReco = new TH2D("tmpReco", "tmpReco", nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

  TH1D *tmpMC1D = new TH1D("tmpMC1D", "tmpMC1D", nBinsJetPtForZHist, JetPtBinsForZHist);
  TH1D *tmpReco1D = new TH1D("tmpReco1D", "tmpReco1D", nBinsJetPtForZHist, JetPtBinsForZHist);

  const int ynbins = nBinsJetPtForMCZHist*nBinsForMCZHist; // MC
  const int xnbins = nBinsJetPtForZHist*nBinsForZHist; // Reco

  // cout  << "Bins = " << xnbins << "\t" << ynbins << endl;

  // Temporary Matrices to find bins

  TH1D *hMCPt = new TH1D("hMCPt", "hMCPt", nBinsJetPtForMCZHist, JetPtBinsForMCZHist);
  TH1D *hMCZ = new TH1D("hMCZ", "hMCZ", nBinsForMCZHist, ZBinsForMCZHist);
  TH1D *hRecoPt = new TH1D("hRecoPt", "hRecoPt", nBinsJetPtForZHist, JetPtBinsForZHist);
  TH1D *hRecoZ = new TH1D("hRecoZ", "hRecoZ", nBinsForZHist, ZBinsForZHist);

  TH2D *NewResponseMatrix[3];
  TH1D *NewMC[3];
  TH1D *NewReco[3];

  // TH2D *NewResponseMatrix = new TH2D("NewResponseMatrix", xnbins, 0, xnbins, ynbins, 0, ynbins);
  // TH1D *NewMC = new TH1D("NewMC", "NewMC", ynbins, 0, ynbins);
  // TH1D *NewReco = new TH1D("NewReco", "NewReco", xnbins, 0, xnbins);

  TH2D *resp2MC[3];
  TH2D *resp2Reco[3];
  THn  *resp2THn[3];

  double xbins[xnbins+1] = {};
  double ybins[ynbins+1] = {};

  for (int i = 0; i <= xnbins; i++){
    xbins[i] = i+1;
    // cout << xbins[i] << endl; 
  }
  for (int i = 0; i <= ynbins; i++){
    ybins[i] = i+1; 
  }

  for (int i = 0; i < 3; i++){
    hVzMC[i] = new TH1D(Form("VzMC_%i", i), Form("VzMC_%i", i), 200, -20, 20); 

    resp2[i] = new RooUnfoldResponse(tmpReco, tmpMC, Form("Response2NewDef_%i", i), Form("Response2NewDef_%i", i));
    resp2MC[i] = new TH2D(Form("resp2MC_%i", i), Form("resp2MC_%i", i), nBinsJetPtForMCZHist, JetPtBinsForMCZHist, nBinsForMCZHist, ZBinsForMCZHist);
    resp2Reco[i] = new TH2D(Form("resp2Reco_%i", i), Form("resp2Reco_%i", i), nBinsJetPtForZHist, JetPtBinsForZHist, nBinsForZHist, ZBinsForZHist);

    NewResponseMatrix[i] = new TH2D(Form("NewResponseMatrix_%i", i), Form("NewResponseMatrix_%i", i), xnbins, xbins, ynbins, ybins);
    NewMC[i] = new TH1D(Form("NewMC_%i",i), Form("NewMC_%i",i), ynbins, ybins);
    NewReco[i] = new TH1D(Form("NewReco_%i",i), Form("NewReco_%i",i), xnbins, xbins);

    newresp[i] = new RooUnfoldResponse(NewReco[i], NewMC[i], Form("NewResp_%i", i), Form("NewResp_%i", i));
  }

  NEvents = new TH1D("NEvents", "NEvents", 3, -0.5, 2.5);
  NEventsWeighted = new TH1D("NEventsWeighted", "NEventsWeighted", 3, -0.5, 2.5);

  cout << JetTree->GetEntries() << endl;

  int nentries = JetTree->GetEntries();

  TH1D *Weight[3];
  for (int i = 0; i < 3; i++){
    Weight[i] = new TH1D(Form("Weight_%i", i), Form("Weight_%i", i), nBinsForMCZHist, ZBinsForMCZHist);

    for (int bin = Weight[i]->FindBin(0.001); bin <= Weight[i]->FindBin(0.999); bin++){
      Weight[i]->SetBinContent(bin, 1.);
    }
  }

  TFile *f1;
  TFile *f2;

  TH1D *TrueZDist[3];
  TH1D *PriorZ[3];

  if (SUPERITERATION > 0){
    f1 = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", 0, iteration));
    f1->cd();

    TH1D *TrueZDist[3];
    for (int i = 0; i < 3; i++){
      TrueZDist[i] = (TH1D *)f1->Get(Form("MCZ_%i", i));
    }

    f2 = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION-1, iteration));
    f2->cd();
    TH1D *PriorZ[3];

    for (int i = 0; i < 3; i++){
      PriorZ[i] = (TH1D *)f2->Get(Form("UnfoldedZ_%i", i));
      PriorZ[i]->Scale(1./PriorZ[i]->Integral());
      // TrueZDist[i]->Scale(PriorZ[i]->Integral()/TrueZDist[i]->Integral());
      TrueZDist[i]->Scale(1./TrueZDist[i]->Integral());
    }

    // cout << "Problem Doesn't Start Here" << endl;
    if (SUPERITERATION == 1){
      for (int i = 0; i < 3; i++){
        cout << "True Z Integral = " <<  TrueZDist[i]->Integral() << endl;
        for (int bin = Weight[i]->FindBin(0.001); bin <= Weight[i]->FindBin(0.999); bin++){
          if (bin >= Weight[i]->FindBin(0.001) && bin <= Weight[i]->FindBin(0.999))
            // Weight[i]->SetBinContent(bin, TrueZDist[i]->Integral()/(Weight[i]->FindBin(0.999) - Weight[i]->FindBin(0.001)));
            Weight[i]->SetBinContent(bin, 1.0);
            // Weight[i]->SetBinContent(bin, PriorZ[i]->Integral());
          else
            Weight[i]->SetBinContent(bin, 0.);
        }
        // Weight[i]->Scale(1./Weight[i]->Integral());
        Weight[i]->Divide(TrueZDist[i]);
        cout << "Weight Z  = " <<  Weight[i]->Integral() << endl;
        // Weight[i]->Scale(TrueZDist[i]->Integral()/Weight[i]->Integral());
        // cout << "Weight Z  = " <<  Weight[i]->Integral() << endl;
      }  
    }
    else{
      TFile *f2 = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION-1, iteration));
      cout << f2->GetName() << endl;
      f2->cd();
      for (int i = 0; i < 3; i++){
        Weight[i] = (TH1D *)PriorZ[i]->Clone();
        // Weight[i]->Scale(TrueZDist[i]->Integral()/Weight[i]->Integral());
        Weight[i]->Divide(TrueZDist[i]);
        // Weight[i]->Scale(PriorZ[i]->Integral()/Weight[i]->Integral());
        Weight[i]->Draw("SAME");
        Weight[i]->SetNameTitle(Form("Weight_%i", i), Form("Weight_%i", i));

        cout << "Weight Z  = " <<  Weight[i]->Integral() << endl;
      }
      // cout << "Problem Doesn't Start Here" << endl;      
    }
  }

  // for (int i = 0; i < 3; i++){
  //   for (int bin = 1; bin <= Weight[i]->GetNbinsX(); bin++){
  //     cout << Weight[i]->GetBinContent(bin) << endl;
  //   }
  // }

  int lowlimit = 0;
  int highlimit = 100;

  // nentries = 100;

  if (SUPERITERATION == 0){lowlimit = 0; highlimit = nentries;}
  else {lowlimit = nentries/2; highlimit = nentries;}

  int counter = 0;

  for (int i = lowlimit; i < highlimit; i++){

    if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    // cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    // cout << Centrality << "\t" << MCJetPt << "\t" << RecoJetCorrPt << endl;

    // if (abs(MCJetEta) > 1.) continue;
    // if (abs(RecoPrimaryVertex->at(2)) > 6) continue;
    // if (abs(MCPrimaryVertex->at(2)) > 6) continue;

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

    double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);


    //Defining response matrix using fakes and misses to check closure

    bool isMCD0Pt = MCD0Pt > 5 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 5 && RecoD0Pt < 10;
    bool isMCZ   = mcz > 0 && mcz < 1; // Not used 
    bool isRecoZ = recoz > 0 && recoz < 1; // Not used
    // bool isMCJetPt = MCJetPt > 5 && MCJetPt < 30;
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    // bool isRecoJetPt = RecoJetCorrPt > -7 && RecoJetCorrPt < 30; // This is used because 99% of jets with pT > 5 GeV is captured within this range
    bool isRecoJetPt = RecoJetCorrPt > 0 && RecoJetCorrPt < 30; // 78% of jets with pT > 5 GeV is captured within this range
    bool isMCJetEta = abs(MCJetEta) < 0.6;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;

    double w = Weight[centhistogramtofill]->GetBinContent(Weight[centhistogramtofill]->FindBin(mcz));

    if ( !isMCJetPt || !isRecoJetPt || !isMCJetEta || !isRecoJetEta || !isMCD0Pt || !isRecoD0Pt ) continue;

    // Response

    int mcptbin = hMCPt->FindBin(MCJetPt);
    int mczbin  = hMCZ->FindBin(mcz);
    int recoptbin = hRecoPt->FindBin(RecoJetCorrPt);
    int recozbin = hRecoZ->FindBin(recoz);

    int xbin = (recoptbin - 1)*nBinsJetPtForZHist + recozbin;
    int ybin = (mcptbin - 1)*nBinsJetPtForMCZHist + mczbin;

    // cout << "Info: " << mcptbin << "\t" << mczbin << "\t" << recoptbin << "\t" << recozbin << "\t" << xbin << "\t" << ybin << endl;

    newresp[centhistogramtofill]->Fill(xbin, ybin, w);
    NewResponseMatrix[centhistogramtofill]->Fill(xbin, ybin, w);
    NewMC[centhistogramtofill]->Fill(ybin, w);
    NewReco[centhistogramtofill]->Fill(xbin, w);

    if (isMCD0Pt && isRecoD0Pt && isMCJetPt && isRecoJetPt && isMCJetEta && isRecoJetEta) {
      if (SUPERITERATION == 0){
        if (i < nentries/2){
          resp2MC[centhistogramtofill]->Fill(MCJetPt, mcz, w);
          resp2Reco[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
        }
        else{
          resp2[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
        }
      }

      else{
        resp2[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
        resp2MC[centhistogramtofill]->Fill(MCJetPt, mcz, w);
        resp2Reco[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      }
    }
    counter++;
    // Fake

    if ((!isMCD0Pt || !isMCJetPt || !isMCJetEta) && isRecoD0Pt && isRecoJetPt && isRecoJetEta) {
      cout << "HALT 1" << endl;
      return;
      resp2[centhistogramtofill]->Fake(RecoJetCorrPt, recoz, w);
      resp2Reco[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
    }

    // Miss

    if (isMCD0Pt && isMCJetPt && isMCJetEta && (!isRecoD0Pt || !isRecoJetPt || !isRecoJetEta)) {
      cout << "HALT 2" << endl;
      return;
      resp2[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      resp2MC[centhistogramtofill]->Fill(MCJetPt, mcz, w);
    }

    NEvents->Fill(centhistogramtofill);
    NEventsWeighted->Fill(centhistogramtofill, w);

    // hVzMC[centhistogramtofill]->Fill(MCPrimaryVertex->at(2), w); // Do we need the w here? I guess so, because the events need to be weighted.
    hVzMC[centhistogramtofill]->Fill(MCPrimaryVertex->at(2), w);
  }

  for (int cent = 0; cent < 3; cent++){
    for (int binx = 1; binx <= NewResponseMatrix[cent]->GetNbinsX(); binx++){
      for (int biny = 1; biny <= NewResponseMatrix[cent]->GetNbinsY(); biny++){
        double mcbinval = NewMC[cent]->GetBinContent(binx);
        double mcbinerr = NewMC[cent]->GetBinError(binx);
        double respval = NewResponseMatrix[cent]->GetBinContent(binx, biny);
        double resperr = NewResponseMatrix[cent]->GetBinError(binx, biny);

        double finalval = (mcbinval != 0) ? respval/mcbinval : 0;
        double finalerr = (mcbinval != 0 && respval != 0) ? TMath::Sqrt(pow(mcbinerr/mcbinval,2) + pow(resperr/respval, 2)) : 0;

        NewResponseMatrix[cent]->SetBinContent(binx, biny, finalval);
        NewResponseMatrix[cent]->SetBinError(binx, biny, finalerr);
      }
    }

    newresp2[cent] = new RooUnfoldResponse(NewReco[cent], NewMC[cent], NewResponseMatrix[cent], Form("NewResp2_%i", cent), Form("NewResp2_%i", cent));
  }

  cout << "Counted " << counter << " events" << endl;

  TFile *outfile = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION, iteration), "UPDATE");
  outfile->mkdir(Form("pthatbin_%i", pthatbin));
  outfile->cd(Form("pthatbin_%i", pthatbin));

  NEvents->Write();
  NEventsWeighted->Write();

  for (int i = 0; i < 3; i++){
    hVzMC[i]->Write();
    gDirectory->WriteObject(resp2[i], Form("Response2Reco_%i", i));
    resp2MC[i]->Write();
    resp2Reco[i]->Write();

    gDirectory->WriteObject(newresp[i], Form("NewResponse_%i", i));
    gDirectory->WriteObject(newresp2[i], Form("NewResponseNormalised_%i", i));
    NewResponseMatrix[i]->Write();
    NewMC[i]->Write();
    NewReco[i]->Write();
  }


  if (pthatbin == 0){
    outfile->cd();
    for (int i = 0; i < 3; i++){
      Weight[i]->Write();
    }
  }
  cout << "SI # " << SUPERITERATION << " Iter = " << iteration << "\t" << Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION, iteration) << endl;

  outfile->Close();
  JetTree->Reset();
  f->Close();

  if (SUPERITERATION > 0){
    f2->Close();
    f1->Close();  
  }

  cout << "pTHatBin = " << Form("pthatbin_%i", pthatbin) << endl;
}


void ResponseMaker(int SUPERITERATION = 0, int iteration = 0){
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  cout << "SI # " << SUPERITERATION << " Iter = " << iteration << "\t" << Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION, iteration) << endl;
    
  TFile *outfile = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION, iteration), "RECREATE");
  outfile->Close();
  Method(0, SUPERITERATION, iteration);
  // Method(1, SUPERITERATION, iteration);

  cout << "Response Maker" << endl;

  TFile *responseoutfile = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION, iteration), "UPDATE");
  responseoutfile->cd();

  TH1D *NEntries[2];
  TH1D *VzMC[3][2];

  TH2D *Resp2MC[3][2];
  TH2D *Resp2Reco[3][2];

  RooUnfoldResponse *response2[3][2];
  
  RooUnfoldResponse *newresponse[3][2];
  RooUnfoldResponse *newresponsenorm[3][2];
  TH1D *NewResp2MC[3][2];
  TH1D *NewResp2Reco[3][2];

  // double ratioofeventsthatpassveto[2] = {1070./1000000., 11114./1000000.}; // For scaling, we need to make sure that vetoed events are also included.
  double ratioofeventsthatpassveto[2] = {11114./1000000., 1070./1000000.}; // For scaling, we need to make sure that vetoed events are also included.
  double crosssection[2] = {1.969/66.809, 64.84/66.809};

  for (int bin = 0; bin < 1; bin++){
    responseoutfile->cd(Form("pthatbin_%i", bin));
    // gROOT->ProcessLine(".ls");
    NEntries[bin] = (TH1D *)gDirectory->Get("NEvents");

    for (int i = 0; i < 3; i++){
      VzMC[i][bin] = (TH1D *)gDirectory->Get(Form("VzMC_%i", i));

      // double scale = crosssection[bin]*ratioofeventsthatpassveto[bin]/VzMC[i][bin]->Integral();
      double scale = 1.;

      gDirectory->GetObject(Form("Response2Reco_%i", i), response2[i][bin]);
      Resp2MC[i][bin] = (TH2D *)gDirectory->Get(Form("resp2MC_%i", i));
      Resp2Reco[i][bin] = (TH2D *)gDirectory->Get(Form("resp2Reco_%i", i));

      response2[i][bin]->Scale(scale);
      Resp2MC[i][bin]->Scale(scale);
      Resp2Reco[i][bin]->Scale(scale);


      gDirectory->GetObject(Form("NewResponse_%i", i), newresponse[i][bin]);
      gDirectory->GetObject(Form("NewResponseNormalised_%i", i), newresponsenorm[i][bin]);
      NewResp2MC[i][bin] = (TH1D *)gDirectory->Get(Form("NewMC_%i", i));
      NewResp2Reco[i][bin] = (TH1D *)gDirectory->Get(Form("NewReco_%i", i));

    }
  }

  cout << "Got all the new responses" << endl;

  for (int i = 0; i < 3; i++){
    // response2[i][0]->Add(*response2[i][1]);
    response2[i][0]->SetNameTitle(Form("Response_%i", i), Form("Response_%i", i));
    // Resp2MC[i][0]->Add(Resp2MC[i][1]);
    Resp2MC[i][0]->SetNameTitle(Form("MC_%i", i), Form("MC_%i", i));
    // Resp2Reco[i][0]->Add(Resp2Reco[i][1]);
    Resp2Reco[i][0]->SetNameTitle(Form("Reco_%i", i), Form("Reco_%i", i));
  }

  for (int i = 0; i < 3; i++){
    // response2[i][0]->Add(*response2[i][1]);
    newresponse[i][0]->SetNameTitle(Form("Response2_%i", i), Form("Response2_%i", i));
    newresponsenorm[i][0]->SetNameTitle(Form("ResponseNormalised_%i", i), Form("ResponseNormalised_%i", i));
    // Resp2MC[i][0]->Add(Resp2MC[i][1]);
    NewResp2MC[i][0]->SetNameTitle(Form("MC2_%i", i), Form("MC2_%i", i));
    // Resp2Reco[i][0]->Add(Resp2Reco[i][1]);
    NewResp2Reco[i][0]->SetNameTitle(Form("Reco2_%i", i), Form("Reco2_%i", i));
  }

  
  responseoutfile->cd();

  for (int i = 0; i < 3; i++){
    response2[i][0]->Write();
    Resp2MC[i][0]->Write();
    Resp2Reco[i][0]->Write();
  }

  for (int i = 0; i < 3; i++){
    newresponse[i][0]->Write();
    newresponsenorm[i][0]->Write();
    NewResp2MC[i][0]->Write();
    NewResp2Reco[i][0]->Write();
  }

  responseoutfile->Close();
}

void UnfoldMaker(int SUPERITERATION = 0, int iteration = 3){
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  cout << "Unfold Maker" << endl;

  TFile *responseoutfile = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", SUPERITERATION, iteration), "UPDATE");
  responseoutfile->cd();

  RooUnfoldResponse *Response[3];
  RooUnfoldResponse *Response2[3];
  RooUnfoldResponse *Response2Norm[3];

  TH2D *PriorReco[3];
  TH2D *PriorMC[3];

  TH1D *PriorReco2[3];
  TH1D *PriorMC2[3];

  cout << "What happened " << endl;

  for (int i = 0; i < 3; i++){
    gDirectory->GetObject(Form("Response_%i", i), Response[i]);
    PriorMC[i] = (TH2D *)gDirectory->Get(Form("MC_%i", i));
    PriorMC[i]->SetNameTitle(Form("PriorMC_%i", i), Form("PriorMC_%i", i));

    gDirectory->GetObject(Form("Response2_%i", i), Response2[i]);
    gDirectory->GetObject(Form("ResponseNormalised_%i", i), Response2Norm[i]);
    PriorMC2[i] = (TH1D *)gDirectory->Get(Form("MC2_%i", i));
    PriorMC2[i]->SetNameTitle(Form("PriorMC2_%i", i), Form("PriorMC2_%i", i));
  }

  cout << "What happened " << endl;

  TH1D *PriorMCPt[3];
  TH1D *PriorMCZ[3];

  // TH1D *PriorMCPt2[3];
  // TH1D *PriorMCZ2[3];

  for (int i = 0; i < 3; i++){
    PriorMCPt[i] = (TH1D *)PriorMC[i]->ProjectionX();
    PriorMCPt[i]->SetNameTitle(Form("PriorMCPt_%i", i), Form("PriorMCPt_%i", i));
    PriorMCZ[i]  = (TH1D *)PriorMC[i]->ProjectionY();
    PriorMCZ[i]->SetNameTitle(Form("PriorMCZ_%i", i), Form("PriorMCZ_%i", i));

    // PriorMCPt2[i] = (TH1D *)PriorMC2[i]->ProjectionX();
    // PriorMCPt2[i]->SetNameTitle(Form("PriorMCPt2_%i", i), Form("PriorMCPt2_%i", i));
    // PriorMCZ2[i]  = (TH1D *)PriorMC2[i]->ProjectionY();
    // PriorMCZ2[i]->SetNameTitle(Form("PriorMCZ2_%i", i), Form("PriorMCZ2_%i", i));
  }

  cout << "What happened " << endl;

  TFile *OriginalResponseFile = new TFile(Form("TestTrain_JetPt0GeV/Response_%i_%i.root", 0, iteration));
  OriginalResponseFile->cd();

  TH2D *Reco[3];
  TH2D *MC[3];

  for (int i = 0; i < 3; i++){
    Reco[i] = (TH2D *)gDirectory->Get(Form("Reco_%i", i));
    MC[i] = (TH2D *)gDirectory->Get(Form("MC_%i", i));
  }

  TH1D *Reco2[3];
  TH1D *MC2[3];

  for (int i = 0; i < 3; i++){
    Reco2[i] = (TH1D *)gDirectory->Get(Form("Reco2_%i", i));
    MC2[i]   = (TH1D *)gDirectory->Get(Form("MC2_%i", i));
  }

  TH1D *Fakes[3];

  TFile *fin2 = new TFile("Inefficiency.root");
  for (int i = 0; i < 3; i++){
    Fakes[i] = (TH1D *)fin2->Get(Form("Fakes_%i", i));
  }

  // for (int cent = 0; cent < 3; cent++){
  //   for (int binx = 1; binx <= Reco[cent]->GetNbinsX(); binx++){
  //     double weight =  Fakes[cent]->GetBinContent(binx);
  //     for (int biny = 1; biny <= Reco[cent]->GetNbinsY(); biny++){
  //       double val = (1.-weight)*Reco[cent]->GetBinContent(binx, biny);
  //       double err = (1.-weight)*Reco[cent]->GetBinError(binx, biny);

  //       Reco[cent]->SetBinContent(binx, biny, val);
  //       Reco[cent]->SetBinError(binx, biny, err);
  //     }
  //   }
  // }

  TH2D *UnfoldedJetPtvZ[3];

  cout << "Unfolding OG Response" << endl;

  RooUnfoldBayes unfold0 (Response[0], Reco[0], iteration);
  RooUnfoldBayes unfold1 (Response[1], Reco[1], iteration);
  RooUnfoldBayes unfold2 (Response[2], Reco[2], iteration);

  UnfoldedJetPtvZ[0] = (TH2D *)unfold0.Hreco();
  UnfoldedJetPtvZ[1] = (TH2D *)unfold1.Hreco();
  UnfoldedJetPtvZ[2] = (TH2D *)unfold2.Hreco();

  TH1D *MCPt[3];
  TH1D *UnfoldedPt[3];

  TH1D *MCZ[3];
  TH1D *UnfoldedZ[3];

  TH1D *RatioPt[3];
  TH1D *RatioZ[3];

  for (int i = 0; i < 3; i++){
    MCPt[i] = (TH1D *)MC[i]->ProjectionX();
    UnfoldedPt[i] = (TH1D *)UnfoldedJetPtvZ[i]->ProjectionX();

    MCPt[i]->SetNameTitle(Form("MCPt_%i", i), Form("MCPt_%i", i));
    UnfoldedPt[i]->SetNameTitle(Form("UnfoldedPt_%i", i), Form("UnfoldedPt_%i", i));

    RatioPt[i] = (TH1D *)UnfoldedPt[i]->Clone();
    RatioPt[i]->Divide(MCPt[i]);
    RatioPt[i]->SetNameTitle(Form("RatioPt_%i", i), Form("RatioPt_%i", i));

    MCZ[i] = (TH1D *)MC[i]->ProjectionY();
    UnfoldedZ[i] = (TH1D *)UnfoldedJetPtvZ[i]->ProjectionY();

    MCZ[i] -> SetNameTitle(Form("MCZ_%i", i), Form("MCZ_%i", i));
    UnfoldedZ[i] -> SetNameTitle(Form("UnfoldedZ_%i", i), Form("UnfoldedZ_%i", i));

    RatioZ[i] = (TH1D *)UnfoldedZ[i]->Clone();
    RatioZ[i]->Divide(MCZ[i]);
    RatioZ[i]->SetNameTitle(Form("RatioZ_%i", i), Form("RatioZ_%i", i));
  }

  /*

  TH1D *UnfoldedJetPtvZ2[3];

  cout << "Unfolding New Response" << endl;

  RooUnfoldBayes unfold02 (Response2[0], Reco2[0], iteration);
  RooUnfoldBayes unfold12 (Response2[1], Reco2[1], iteration);
  RooUnfoldBayes unfold22 (Response2[2], Reco2[2], iteration);

  UnfoldedJetPtvZ2[0] = (TH1D *)unfold02.Hreco();
  UnfoldedJetPtvZ2[1] = (TH1D *)unfold12.Hreco();
  UnfoldedJetPtvZ2[2] = (TH1D *)unfold22.Hreco();

  TH1D *UnfoldedJetPtvZ2Norm[3];

  cout << "Unfolding New Response Normalised" << endl;

  RooUnfoldBayes unfold02Norm (Response2Norm[0], Reco2[0], iteration);
  RooUnfoldBayes unfold12Norm (Response2Norm[1], Reco2[1], iteration);
  RooUnfoldBayes unfold22Norm (Response2Norm[2], Reco2[2], iteration);

  UnfoldedJetPtvZ2Norm[0] = (TH1D *)unfold02Norm.Hreco();
  UnfoldedJetPtvZ2Norm[1] = (TH1D *)unfold12Norm.Hreco();
  UnfoldedJetPtvZ2Norm[2] = (TH1D *)unfold22Norm.Hreco();

  for (int i = 0; i < 3; i++){
    UnfoldedJetPtvZ2[i]->SetNameTitle(Form("UnfoldedJetPtvZ2_%i", i), Form("UnfoldedJetPtvZ2_%i", i));
    UnfoldedJetPtvZ2Norm[i]->SetNameTitle(Form("UnfoldedJetPtvZ2Norm_%i", i), Form("UnfoldedJetPtvZ2Norm_%i", i));
  }

  */

  responseoutfile->cd();

  for (int i = 0; i < 3; i++){
    PriorMCPt[i]->Write();
    MCPt[i]->Write();
    UnfoldedPt[i]->Write();
    RatioPt[i]->Write();

    PriorMCZ[i]->Write();
    MCZ[i]->Write();
    UnfoldedZ[i]->Write();
    RatioZ[i]->Write();
  }

  /*
  for (int i = 0; i < 3; i++){
    UnfoldedJetPtvZ2[i]->Write();
    UnfoldedJetPtvZ2Norm[i]->Write();
  }
  */
  responseoutfile->Close();
  OriginalResponseFile->Close();  
}

void HIOverlay(int iteration = 3){

  for (int superiter = 0; superiter <= 200; superiter++){
    ResponseMaker(superiter, iteration);
    UnfoldMaker(superiter, iteration);
  }
}