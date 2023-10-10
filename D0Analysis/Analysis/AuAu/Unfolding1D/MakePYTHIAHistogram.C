#include "BinDef.h"
#include "NewBinDef.h"

using namespace std;

void MakePYTHIAHistogram(double lowptcutoff = 1.0, double highptcutoff = 20.0){
	TFile *f = new TFile("/Volumes/WorkDrive/April10_MCResponse_NoJetPtCutOff/MCResponse_April10_NoJetpTCutOff.root");
	f->cd("SimJetSaverSmear_Central");

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

  int nEntries = JetTree->GetEntries();

  TH1D *PYTHIAPt = new TH1D("PYTHIA pT", "PYTHIA pT", 76, lowptcutoff, highptcutoff);
  TH1D *PYTHIAZ = new TH1D("PYTHIA Z", "PYTHIA Z", nz_gen_bins, z_gen_low, z_gen_high);

  TH2D *PYTHIA = new TH2D("PYTHIA", "PYTHIA", 76, lowptcutoff, highptcutoff, nz_gen_bins, z_gen_low, z_gen_high);

  TH1D *FONLLPt = new TH1D("FONLL pT", "FONLL pT", 76, lowptcutoff, highptcutoff);
  TH1D *FONLLZ = new TH1D("FONLL Z", "FONLL Z", nz_gen_bins, z_gen_low, z_gen_high);

  TH2D *FONLLWeighted = new TH2D("FONLL Weighted", "FONLL Weighted", 76, lowptcutoff, highptcutoff, nz_gen_bins, z_gen_low, z_gen_high);

  for (int i = 0; i < nEntries; i++){
  	if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    // cout << "Read Entry " << i << endl;
    JetTree->GetEntry(i);

    bool isMCD0Pt = MCD0Pt > 1 && MCD0Pt < 10;
    bool isMCJetPt = MCJetPt > lowptcutoff && MCJetPt < jetpt_var_bin[njpt_gen_bins_var];

    if (!isMCD0Pt) continue;
    // if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    // if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

    bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue;

    PYTHIAPt->Fill(MCJetPt);
    PYTHIAZ->Fill(mcz);

    PYTHIA->Fill(MCJetPt, mcz);
  }

  PYTHIAPt->Scale(1./PYTHIAPt->Integral());
  PYTHIAZ->Scale(1./PYTHIAZ->Integral());
  PYTHIA->Scale(1./PYTHIA->Integral());

  TFile *FONLL = new TFile(Form("FONLL_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "READ");
  TH1D  *FONLLCurve = (TH1D *)FONLL->Get("FONLL");

  TH1D *Weight = (TH1D *)FONLLCurve->Clone("FONLLWeights");
  Weight->Divide(PYTHIAPt);

  for (int i = 0; i < nEntries; i++){
  	if (i%1000000 == 0) cout << "Read Entry " << i << endl;
    // cout << "Read Entry " << i << endl;
    JetTree->GetEntry(i);

    bool isMCD0Pt = MCD0Pt > 1 && MCD0Pt < 10;
    bool isMCJetPt = MCJetPt > lowptcutoff && MCJetPt < jetpt_var_bin[njpt_gen_bins_var];

    if (!isMCD0Pt) continue;
    // if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    // if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);
    
    bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue;

    double weight = Weight->GetBinContent(Weight->FindBin(MCJetPt));

    FONLLPt->Fill(MCJetPt, weight);
    FONLLZ->Fill(mcz, weight);

    FONLLWeighted->Fill(MCJetPt, mcz, weight);
  }

  FONLLPt->Scale(1./FONLLPt->Integral());
  FONLLZ->Scale(1./FONLLZ->Integral());
  FONLLWeighted->Scale(1./FONLLWeighted->Integral());

  // FONLLWeighted->Fit()

  TCanvas *c = new TCanvas("c", "c", 1200, 600);
  c->Divide(5);
  c->cd(1);
  gPad->SetLogy();
  PYTHIAPt->Draw();
  FONLLPt->Draw("SAME");
  SetColor(FONLLPt, kRed);
  c->cd(2);
  gPad->SetLogy();
  PYTHIAZ->Draw();
  FONLLZ->Draw("SAME");
  SetColor(FONLLZ, kRed);
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

  TFile *g = new TFile(Form("PYTHIA_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "RECREATE");
  g->cd();
  PYTHIAPt->Write();
  PYTHIAZ->Write();
  PYTHIA->Write();
  Weight->Write();
  g->Close();
}