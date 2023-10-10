#include "BinDef.h"
#include "NewBinDef.h"

using namespace std;

void MakeLIDOHistogram(double lowptcutoff = 1.0, double highptcutoff = 20.0){

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

	TFile *f = new TFile("HIOverlay_HFJets_WithCS_Jun20_2023_pthat_3_inf.root");
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
  Float_t         RecoJetPtFromArea;
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
 
  JetTree->SetBranchAddress("RecoJetPt", &RecoJetPtFromArea);
  // JetTree->SetBranchAddress("RecoJetCorrPt", &RecoJetCorrPt);
  JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
  JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
  // JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
  // JetTree->SetBranchAddress("RecoJetE", &RecoJetE);
  // // JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
  // JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConst);
  // // JetTree->SetBranchAddress("RecoD0Z", &RecoD0Z);
  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  // JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  // JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  // JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  // JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  // JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  // JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  int nEntries = JetTree->GetEntries();

  int nptbins = (highptcutoff - lowptcutoff)*4;

  int nz_bins_theory = 9;
    double nbinsz_theory[10] = {0., 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  TH1D *PYTHIAPt = new TH1D("PYTHIA pT", "PYTHIA pT", njpt_gen_bins_var, jetpt_var_bin);
  // TH1D *PYTHIAPt = new TH1D("PYTHIA pT", "PYTHIA pT", 40, 0, 40);
  TH1D *PYTHIAZ = new TH1D("PYTHIA Z", "PYTHIA Z", nz_bins_theory, nbinsz_theory);
  TH2D *PYTHIA = new TH2D("PYTHIA", "PYTHIA", njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);

  TH1D *LIDOPt = new TH1D("LIDO pT", "LIDO pT", njpt_gen_bins_var, jetpt_var_bin);
  TH1D *LIDOZ = new TH1D("LIDO Z", "LIDO Z", nz_bins_theory, nbinsz_theory);
  TH2D *LIDOWeighted = new TH2D("LIDO Weighted", "LIDO Weighted", njpt_gen_bins_var, jetpt_var_bin, nz_bins_theory, nbinsz_theory);
  
  TH1D *PYTHIAPtNormal = new TH1D("PYTHIANormal pT", "PYTHIANormal pT", nptbins, lowptcutoff, highptcutoff);
  TH1D *PYTHIAZNormal = new TH1D("PYTHIAZNormal Z", "PYTHIAZNormal Z", nz_bins_theory, nbinsz_theory);

  for (int i = 0; i < 1000000; i++){
  	if (i%1000000 == 0) cout << "Read Entry " << i << "\r" << flush;
    // cout << "Read Entry " << i << endl;
    JetTree->GetEntry(i);

    bool isMCD0Pt = MCD0Pt > lowptcutoff && MCD0Pt < 10;
    bool isMCJetPt = MCJetPt > lowptcutoff && MCJetPt < highptcutoff;

    if (!isMCD0Pt) continue;
    // if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
    // if (RecoJetNConst==0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

    if (mcz >=1) mcz = 0.9999;

    bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue;

    PYTHIAPt->Fill(MCJetPt);
    PYTHIAZ->Fill(mcz);

    PYTHIAPtNormal->Fill(MCJetPt);
    if (MCJetPt > 5.) PYTHIAZNormal->Fill(mcz);

    PYTHIA->Fill(MCJetPt, mcz);
  }
  cout << endl;
  
  TH1D *LIDOCurvePt[3];
  TH1D *LIDOCurveZ[3];

  TFile *FinalPlots = new TFile("Plots/FinalPlots.root", "READ");
  FinalPlots->cd(Form("D0pT_%i_10", (int)lowptcutoff));
  for (int cent = 0; cent < 3; cent++){
    LIDOCurvePt[cent] = (TH1D *)gDirectory->Get(Form("TheoryD0JetPt_%i_%i", cent, (int)lowptcutoff));
    for (int i = 0; i < LIDOCurvePt[cent]->GetNbinsX(); i++){
      LIDOCurvePt[cent]->SetBinContent(i+1, LIDOCurvePt[cent]->GetBinContent(i+1)*LIDOCurvePt[cent]->GetBinCenter(i+1));
      LIDOCurvePt[cent]->SetBinError(i+1, LIDOCurvePt[cent]->GetBinError(i+1)*LIDOCurvePt[cent]->GetBinCenter(i+1));
    }
    LIDOCurveZ[cent] = (TH1D *)gDirectory->Get(Form("TheoryD0JetZ_%i_%i", cent, (int)lowptcutoff));
    LIDOCurvePt[cent]->Scale(1./LIDOCurvePt[cent]->Integral());
    LIDOCurvePt[cent]->SetDirectory(0);
    LIDOCurveZ[cent]->SetDirectory(0);
  }
  DivideByBinWidth(PYTHIAPt);
  PYTHIAPt->Scale(1./PYTHIAPt->Integral(5,20));

  TF1 *f1 = new TF1("PowerLaw", "[0]*x*(1 + [1]*x)**(-1.0*[2])", 1, 20);
  f1->SetParameters(2, 1, 3);
  LIDOCurvePt[0]->Fit(f1, "RLMES0", "SAME", 5.0, 20.0);

  TF1 *f2 = new TF1("PowerLaw", "[0]*x*(1 + [1]*x)**(-1.0*[2])", 1, 20);
  f2->SetParameters(2, 1, 3);
  LIDOCurvePt[2]->Fit(f2, "RLMES0", "SAME", 5.0, 20.0);

  TCanvas *d = new TCanvas("d", "d", 1800, 900);
  d->Divide(2,1);
  d->cd(1);
  gPad->SetLogy();
  SetColor(LIDOCurvePt[0], kRed, 20);
  PYTHIAPt->GetXaxis()->SetRangeUser(1, 20);
  PYTHIAPt->GetYaxis()->SetRangeUser(1e-6, 1e2);
  PYTHIAPt->Draw("HIST SAME");
  LIDOCurvePt[0]->Draw("EP SAME");
  f1->Draw("SAME");
  d->cd(2);
  gPad->SetLogy();
  SetColor(LIDOCurvePt[2], kRed, 20);
  PYTHIAPt->Draw("HIST SAME");
  LIDOCurvePt[2]->Draw("EP SAME");
  f2->Draw("SAME");

  

  TH1D *LIDONormal[3];
  TH1D *LIDOZNormal[3]; 
  for (int cent = 0; cent < 3; cent++){
    LIDONormal[cent] = new TH1D(Form("LIDONormal pT %i", cent), Form("LIDONormal pT %i", cent), nptbins, lowptcutoff, highptcutoff);
    LIDOZNormal[cent] = new TH1D(Form("LIDOZNormal Z %i", cent), Form("LIDOZNormal Z %i", cent), nz_bins_theory, nbinsz_theory);
  
    for (int i = 0; i < nptbins; i++){
        if (cent == 1) continue;
        if (cent == 0)LIDONormal[cent]->SetBinContent(i+1, f1->Eval(LIDONormal[cent]->GetBinCenter(i+1)));
        if (cent == 2)LIDONormal[cent]->SetBinContent(i+1, f2->Eval(LIDONormal[cent]->GetBinCenter(i+1)));
    }
  }

    PYTHIAPtNormal->Scale(1./PYTHIAPtNormal->Integral());
    PYTHIAZNormal->Scale(1./PYTHIAZNormal->Integral());

  LIDONormal[1] = (TH1D *)LIDONormal[0]->Clone("LIDONormal 1");

  TH1D *Weight[3];

    for (int i = 0; i < 3; i++){
        LIDONormal[i]->Scale(1./LIDONormal[i]->Integral());
        Weight[i] = (TH1D *)LIDONormal[i]->Clone(Form("Weight %i", i));
        Weight[i]->Divide(PYTHIAPtNormal);
    }
    TCanvas *e = new TCanvas("e", "e", 1800, 900);
    e->cd();
    gPad->SetLogy();
    SetColor(Weight[0], kRed, 20);
    Weight[0]->Draw("HIST SAME");
    Weight[2]->Draw("HIST SAME");

    for (int i = 0; i < 1000000; i++){
        if (i%1000000 == 0) cout << "Read Entry " << i << "\r" << flush;
        // cout << "Read Entry " << i << endl;
        JetTree->GetEntry(i);

        bool isMCD0Pt = MCD0Pt > lowptcutoff && MCD0Pt < 10;
        bool isMCJetPt = MCJetPt > lowptcutoff && MCJetPt < highptcutoff;

        if (!isMCD0Pt) continue;
        // if (!isRecoD0Pt) continue;
        if (!isMCJetPt) continue;
        // if (RecoJetNConst==0) continue;

        double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
        double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
        double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
        double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
        double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

        if (mcz >=1) mcz = 0.9999;

        bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
        if (!isMCZ) continue;

        int centhistogramtofill = -99;
        if (Centrality < 10) centhistogramtofill = 0;
        else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
        else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

        if (centhistogramtofill < 0) continue;

        if(MCJetPt > 5) LIDOZNormal[centhistogramtofill]->Fill(mcz, Weight[centhistogramtofill]->GetBinContent(Weight[centhistogramtofill]->FindBin(MCJetPt)));
    }
  cout << endl;

  TH1D *WeightZ[3];

  for (int cent = 0; cent < 3; cent++){
        for (int i = 0; i < LIDOZNormal[cent]->GetNbinsX(); i++){
            LIDOZNormal[cent]->SetBinContent(i+1, LIDOZNormal[cent]->GetBinContent(i+1)*LIDOZNormal[cent]->GetBinCenter(i+1));
            LIDOZNormal[cent]->SetBinError(i+1, LIDOZNormal[cent]->GetBinError(i+1)*LIDOZNormal[cent]->GetBinCenter(i+1));
        }
        LIDOZNormal[cent]->Scale(1./LIDOZNormal[cent]->Integral());
        WeightZ[cent] = (TH1D *)LIDOZNormal[cent]->Clone(Form("Weight Z %i", cent));
        WeightZ[cent]->Divide(PYTHIAZNormal);
    }

    TH1D *LIDOPtNormalWeighted[3];
    TH1D *LIDOZNormalWeighted[3];

    for (int cent = 0; cent < 3; cent++){
        LIDOPtNormalWeighted[cent] = new TH1D(Form("LIDOPtNormalWeighted pT %i", cent), Form("LIDOPtNormalWeighted pT %i", cent), nptbins, lowptcutoff, highptcutoff);
        LIDOZNormalWeighted[cent] = new TH1D(Form("LIDOZNormalWeighted Z %i", cent), Form("LIDOZNormalWeighted Z %i", cent), nz_bins_theory, nbinsz_theory);
    }

    for (int i = 0; i < 1000000; i++){
        if (i%1000000 == 0) cout << "Read Entry " << i << "\r" << flush;
        // cout << "Read Entry " << i << endl;
        JetTree->GetEntry(i);

        bool isMCD0Pt = MCD0Pt > lowptcutoff && MCD0Pt < 10;
        bool isMCJetPt = MCJetPt > lowptcutoff && MCJetPt < highptcutoff;

        if (!isMCD0Pt) continue;
        // if (!isRecoD0Pt) continue;
        if (!isMCJetPt) continue;
        // if (RecoJetNConst==0) continue;

        double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
        double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
        double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
        double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
        double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

        if (mcz >=1) mcz = 0.9999;

        bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
        if (!isMCZ) continue;

        int centhistogramtofill = -99;
        if (Centrality < 10) centhistogramtofill = 0;
        else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
        else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

        if (centhistogramtofill < 0) continue;

        double w = Weight[centhistogramtofill]->GetBinContent(Weight[centhistogramtofill]->FindBin(MCJetPt));
        w *= WeightZ[centhistogramtofill]->GetBinContent(WeightZ[centhistogramtofill]->FindBin(mcz));

        LIDOPtNormalWeighted[centhistogramtofill]->Fill(MCJetPt, w);
        if(MCJetPt > 5) LIDOZNormalWeighted[centhistogramtofill]->Fill(mcz, w);
    }
  cout << endl;

    for (int cent = 0; cent < 3; cent++){
        for (int i = 0; i < LIDOPtNormalWeighted[cent]->GetNbinsX(); i++){
            LIDOPtNormalWeighted[cent]->SetBinContent(i+1, LIDOPtNormalWeighted[cent]->GetBinContent(i+1)*LIDOPtNormalWeighted[cent]->GetBinCenter(i+1));
            LIDOPtNormalWeighted[cent]->SetBinError(i+1, LIDOPtNormalWeighted[cent]->GetBinError(i+1)*LIDOPtNormalWeighted[cent]->GetBinCenter(i+1));
        }
        LIDOPtNormalWeighted[cent]->Scale(1./LIDOPtNormalWeighted[cent]->Integral());
    }

    for (int cent = 0; cent < 3; cent++){
        for (int i = 0; i < LIDOZNormalWeighted[cent]->GetNbinsX(); i++){
            LIDOZNormalWeighted[cent]->SetBinContent(i+1, LIDOZNormalWeighted[cent]->GetBinContent(i+1)*LIDOZNormalWeighted[cent]->GetBinCenter(i+1));
            LIDOZNormalWeighted[cent]->SetBinError(i+1, LIDOZNormalWeighted[cent]->GetBinError(i+1)*LIDOZNormalWeighted[cent]->GetBinCenter(i+1));
        }
        LIDOZNormalWeighted[cent]->Scale(1./LIDOZNormalWeighted[cent]->Integral());
    }

    for (int cent = 0; cent < 3; cent++){
        for (int i = 0; i < LIDOCurveZ[cent]->GetNbinsX(); i++){
            LIDOCurveZ[cent]->SetBinContent(i+1, LIDOCurveZ[cent]->GetBinContent(i+1)*LIDOCurveZ[cent]->GetBinCenter(i+1));
            LIDOCurveZ[cent]->SetBinError(i+1, LIDOCurveZ[cent]->GetBinError(i+1)*LIDOCurveZ[cent]->GetBinCenter(i+1));
        }
        LIDOCurveZ[cent]->Scale(1./LIDOCurveZ[cent]->Integral());
        LIDOZNormal[cent]->Scale(1./LIDOZNormal[cent]->Integral());
    }

    TCanvas *e1 = new TCanvas("e1", "e1", 1800, 900);
    e1->Divide(2,1);
    e1->cd(1);
    gPad->SetLogy();
    SetColor(LIDONormal[0], kRed, 20);
    SetColor(LIDOPtNormalWeighted[0], kBlue, 20);
    SetColor(PYTHIAPtNormal, kBlack, 20);
    LIDONormal[0]->Draw("EP SAME");
    LIDOPtNormalWeighted[0]->Draw("EP SAME");
    PYTHIAPtNormal->Draw("EP SAME");
    e1->cd(2);
    gPad->SetLogy();
    SetColor(LIDONormal[2], kRed, 20);
    SetColor(LIDOPtNormalWeighted[2], kBlue, 20);
    LIDONormal[2]->Draw("EP SAME");
    LIDOPtNormalWeighted[2]->Draw("EP SAME");
    PYTHIAPtNormal->Draw("EP SAME");

    TCanvas *e2 = new TCanvas("e2", "e2", 1800, 900);
    e2->Divide(2,1);
    e2->cd(1);
    gPad->SetLogy();
    SetColor(LIDOCurveZ[0], kRed, 20);
    SetColor(LIDOZNormalWeighted[0], kBlue, 20);
    SetColor(PYTHIAZNormal, kBlack, 20);
    LIDOCurveZ[0]->Draw("EP SAME");
    LIDOZNormalWeighted[0]->Draw("EP SAME");
    PYTHIAZNormal->Draw("EP SAME");
    e2->cd(2);
    gPad->SetLogy();
    SetColor(LIDOCurveZ[2], kRed, 20);
    SetColor(LIDOZNormalWeighted[2], kBlue, 20);
    LIDOCurveZ[2]->Draw("EP SAME");
    LIDOZNormalWeighted[2]->Draw("EP SAME");
    PYTHIAZNormal->Draw("EP SAME");

    TFile *WeightFile = new TFile("LIDOHistogram_1_20.root", "RECREATE");
    WeightFile->cd();
    for (int cent = 0; cent < 3; cent++){
        Weight[cent]->Write();
    }


}