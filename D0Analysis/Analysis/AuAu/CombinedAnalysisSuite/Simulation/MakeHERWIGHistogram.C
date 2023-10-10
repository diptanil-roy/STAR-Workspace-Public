#include "BinDef.h"
#include "NewBinDef.h"

using namespace std;

void MakeHERWIGHistogram(double lowptcutoff = 1.0, double highptcutoff = 30.0){

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *herwig = new TFile("HERWIG_D0Jets.root");
  herwig->cd("TagD0Event");
  TH2D *HERWIGJetPtVsZ = (TH2D *)gDirectory->Get(Form("D0ZJetPt_D0PtAbove%i", (int)lowptcutoff));

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
 
  // JetTree->SetBranchAddress("RecoJetPt", &RecoJetPt);
  // JetTree->SetBranchAddress("RecoJetCorrPt", &RecoJetCorrPt);
  // JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
  // JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
  // JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
  // JetTree->SetBranchAddress("RecoJetE", &RecoJetE);
  // // JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
  // JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConst);
  // // JetTree->SetBranchAddress("RecoD0Z", &RecoD0Z);
  // JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  // JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  // JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  // JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  // JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  // JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  // JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  // JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  // JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  int nEntries = JetTree->GetEntries();

  int nptbins = (highptcutoff - lowptcutoff)*4;

  TH1D *PYTHIAPt = new TH1D("PYTHIA pT", "PYTHIA pT", 40, 0.0, 40.0);
  TH1D *PYTHIAZ = new TH1D("PYTHIA Z", "PYTHIA Z", 10, 0.0, 1.0);

  TH2D *PYTHIA = new TH2D("PYTHIA", "PYTHIA", 40, 0.0, 40.0, 10, 0.0, 1.0);

  TH1D *HERWIGPt = new TH1D("HERWIG pT", "HERWIG pT", 40, 0.0, 40.0);
  TH1D *HERWIGZ = new TH1D("HERWIG Z", "HERWIG Z", 10, 0.0, 1.0);

  TH2D *HERWIG = new TH2D("HERWIG", "HERWIG", 40, 0.0, 40.0, 10, 0.0, 1.0);

  TH1D *FitPt = new TH1D("HERWIG Fit pT", "HERWIG FONLL pT", 40, 0.0, 40.0);
  TH1D *FitZ = new TH1D("HERWIG Fit Z", "HERWIG FONLL Z", 10, 0.0, 1.0);

  TH2D *FitWeighted = new TH2D("HERWIG Fit Weighted", "HERWIG Fit Weighted", 40, 0.0, 40.0, 10, 0.0, 1.0);

  TF1 *f1 = new TF1("PowerLaw", "[0]*x*(1 + [1]*x)**(-1.0*[2])", lowptcutoff, highptcutoff);
  
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

    PYTHIA->Fill(MCJetPt, mcz);
  }
  cout << endl;

  TH2D *Weight = (TH2D *)HERWIGJetPtVsZ->Clone("HERWIGWeights");
  cout << Weight->GetBinCenter(0) << " " << Weight->GetNbinsY() << endl;
  cout << PYTHIA->GetNbinsX() << " " << PYTHIA->GetNbinsY() << endl;
  SetName(Weight, "HERWIGWeights");
  Weight->Divide(PYTHIA);

  cout << "Weights Integral = " << Weight->Integral() << endl;

  Weight->Draw("COLZ");

  /*
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

    bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue;

    double w = Weight->GetBinContent(Weight->FindBin(MCJetPt, mcz));

    HERWIGPt->Fill(MCJetPt, w);
    HERWIGZ->Fill(mcz, w);

    HERWIG->Fill(MCJetPt, mcz, w);
  }
  cout << endl;
  cout << "Start Fit" << endl;
  f1->SetParameters(5000, 1, 1);
  HERWIGPt->Fit(f1, "RLME0S", "SAME", 5.0, 20.0);

  TH1D *FitHist = new TH1D("HERWIG Fit Histogram", "HERWIG Fit Histogram", 40, 0.0, 40.0);

  for (int i = 1; i <= FitHist->GetNbinsX(); i++ ){
    double pt = FitHist->GetBinCenter(i);
    double weight = f1->Eval(pt);
    FitHist->SetBinContent(i, weight);
  }

  TCanvas *d2 = new TCanvas("d2", "d2", 1800, 900);
  d2->cd();
  d2->SetLogy();
  SetColor(FitHist, kGreen-2, 20);
  TH1D *FitHistClone = (TH1D *)FitHist->Clone("FitHistClone");
  FitHistClone->Draw("EP");
  FitHistClone->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
  FitHistClone->GetYaxis()->SetTitle("Counts");
  TH1D *HerwigPtClone = (TH1D *)HERWIGPt->Clone("HerwigPtClone");
  HerwigPtClone->Draw("HIST SAME");
  auto legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->AddEntry("HerwigPtClone", "HERWIG", "l");
  legend->AddEntry("FitHistClone", "Power Law Fit", "p");
  legend->Draw("SAME");
  d2->SaveAs("HERWIG_Pt_Fit.pdf");

  FitHist->Scale(1./FitHist->Integral());
  HERWIGPt->Scale(1./HERWIGPt->Integral());
  HERWIGZ->Scale(1./HERWIGZ->Integral());
  HERWIG->Scale(1./HERWIG->Integral());

  TH1D *WeightFromFit = (TH1D *)FitHist->Clone("FitWeights");
  WeightFromFit->Divide(PYTHIAPt);

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
    
    bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue;

    double weightfromfit = WeightFromFit->GetBinContent(WeightFromFit->FindBin(MCJetPt));

    FitPt->Fill(MCJetPt, weightfromfit);
    if (MCJetPt > 5.0) FitZ->Fill(mcz, weightfromfit);
    FitWeighted->Fill(MCJetPt, mcz, weightfromfit);
  }
  cout << endl;

  FitPt->Scale(1./FitPt->Integral());
  FitZ->Scale(1./FitZ->Integral());
  FitWeighted->Scale(1./FitWeighted->Integral());

  // FONLLWeighted->Fit()

  TCanvas *c = new TCanvas("c", "c", 1800, 900);
  c->Divide(3);
  c->cd(1);
  gPad->SetLogy();
  PYTHIAPt->GetYaxis()->SetRangeUser(pow(10, -8), 1);
  PYTHIAPt->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
  PYTHIAPt->GetYaxis()->SetTitle("a.u.");
  SetColor(HERWIGPt, kRed, 20);
  SetColor(FitPt, kGreen-2, 20);
  PYTHIAPt->Draw("HIST");
  HERWIGPt->Draw("EP SAME");
  FitPt->Draw("EP SAME");
  auto legend2 = new TLegend(0.5,0.6,0.9,0.9);
  legend2->AddEntry("PYTHIAPt", "PYTHIA", "l");
  legend2->AddEntry("HERWIGPt", "HERWIGPt", "p");
  legend2->AddEntry("FitPt", "Power Law Fit for HERWIG", "p");
  legend2->Draw("SAME");
  c->cd(2);
  gPad->SetLogy();
  PYTHIAZ->Draw("HIST");
  PYTHIAZ->GetXaxis()->SetTitle("z");
  PYTHIAZ->GetYaxis()->SetTitle("a.u.");
  HERWIGZ->Draw("SAME");
  FitZ->Draw("SAME");
  SetColor(HERWIGZ, kRed);
  SetColor(FitZ, kGreen-2);
  c->cd(3);
  gPad->SetLogy();
  Weight->Draw();
  Weight->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
  WeightFromFit->Draw("SAME");
  SetColor(Weight, kRed);
  SetColor(WeightFromFit, kGreen-2);

  c->SaveAs("HERWIG_Weighted.pdf");

  TCanvas *d = new TCanvas("d", "d", 1800, 900);
  d->Divide(3);
  d->cd(1);
  gPad->SetLogz();
  PYTHIA->Draw("COLZ");
  PYTHIA->GetZaxis()->SetRangeUser(pow(10, -6),1);
  d->cd(2);
  gPad->SetLogz();
  HERWIG->Draw("COLZ");
  HERWIG->GetZaxis()->SetRangeUser(pow(10, -6),1);
  d->cd(3);
  gPad->SetLogz();
  FitWeighted->Draw("COLZ");
  FitWeighted->GetZaxis()->SetRangeUser(pow(10, -6),1);
  
  */
//   TFile *g = new TFile(Form("HERWIG_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "RECREATE");
//   g->cd();
//   PYTHIAPt->Write();
//   PYTHIAZ->Write();
//   PYTHIA->Write();
//   Weight->Write();
//   WeightFromFit->Write();
//   g->Close();

}