#include "BinDef.h"
#include "NewBinDef.h"

using namespace std;

void MakePYTHIAHistogram(double lowptcutoff = 1.0, double highptcutoff = 30.0){

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

  TH1D *PYTHIAPt = new TH1D("PYTHIA pT", "PYTHIA pT", nptbins, lowptcutoff, highptcutoff);
  // TH1D *PYTHIAPt = new TH1D("PYTHIA pT", "PYTHIA pT", 40, 0, 40);
  TH1D *PYTHIAZ = new TH1D("PYTHIA Z", "PYTHIA Z", nz_gen_bins, z_gen_bin);

  TH2D *PYTHIA = new TH2D("PYTHIA", "PYTHIA", nptbins, lowptcutoff, highptcutoff, nz_gen_bins, z_gen_bin);

  TH1D *FONLLPt = new TH1D("FONLL pT", "FONLL pT", nptbins, lowptcutoff, highptcutoff);
  TH1D *FONLLZ = new TH1D("FONLL Z", "FONLL Z", nz_gen_bins, z_gen_bin);

  TH1D *FitPt = new TH1D("PYTHIA Fit pT", "PYTHIA FONLL pT", nptbins, lowptcutoff, highptcutoff);
  TH1D *FitZ = new TH1D("PYTHIA Fit Z", "PYTHIA FONLL Z", nz_gen_bins, z_gen_bin);

  TH1D *DataWeighedPt[3];
  TH1D *DataWeighedZ[3];
  TH1D *DataWeighedPtFit[3];
  TH1D *DataWeighedZFit[3];
  TH2D *DataWeighed[3];

  for (int i = 0; i < 3; i++){
    DataWeighedPt[i] = new TH1D(Form("DataWeighed pT %i", i), Form("DataWeighed pT %i", i), nptbins, lowptcutoff, highptcutoff);
    DataWeighedZ[i] = new TH1D(Form("DataWeighed Z %i", i), Form("DataWeighed Z %i", i), nz_gen_bins, z_gen_bin);
    DataWeighedPtFit[i] = new TH1D(Form("DataWeighed Fit pT %i", i), Form("DataWeighed Fit pT %i", i), nptbins, lowptcutoff, highptcutoff);
    DataWeighedZFit[i] = new TH1D(Form("DataWeighed Fit Z %i", i), Form("DataWeighed Fit Z %i", i), nz_gen_bins, z_gen_bin);
    DataWeighed[i] = new TH2D(Form("Data Weighted %i", i), Form("Data Weighted %i", i), nptbins, lowptcutoff, highptcutoff, nz_gen_bins, z_gen_bin);
  }
  
  TH2D *FONLLWeighted = new TH2D("FONLL Weighted", "FONLL Weighted", nptbins, lowptcutoff, highptcutoff, nz_gen_bins, z_gen_bin);
  TH2D *FitWeighted = new TH2D("PYTHIA Fit Weighted", "PYTHIA Fit Weighted", nptbins, lowptcutoff, highptcutoff, nz_gen_bins, z_gen_bin);


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

  // PYTHIAPt->Scale(1./PYTHIAPt->Integral());

  int low = PYTHIAPt->GetXaxis()->FindBin(5.0);
  int high = PYTHIAPt->GetXaxis()->FindBin(30.0);

  f1->SetParameters(10000, 1, 3);
  // f1->SetParLimits(0, 10000, 100000);
  // f1->SetParLimits(2, 0, 10);
  PYTHIAPt->Fit(f1, "RLME0S", "SAME", 7.0, 20.0);

  TH1D *FitHist = new TH1D("PYTHIA Fit Histogram", "PYTHIA Fit Histogram", nptbins, lowptcutoff, highptcutoff);

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
  TH1D *PythiaPtClone = (TH1D *)PYTHIAPt->Clone("PythiaPtClone");
  PythiaPtClone->Draw("HIST SAME");
  auto legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->AddEntry("PythiaPtClone", "PYTHIA", "l");
  legend->AddEntry("FitHistClone", "Power Law Fit", "p");
  legend->Draw("SAME");
  d2->SaveAs(Form("PYTHIA_Pt_Fit_%i.pdf", (int)lowptcutoff));

  FitHist->Scale(1./FitHist->Integral());
  PYTHIAPt->Scale(1./PYTHIAPt->Integral());
  PYTHIAZ->Scale(1./PYTHIAZ->Integral());
  PYTHIA->Scale(1./PYTHIA->Integral());
  
  TFile *FONLL = new TFile(Form("FONLL_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "READ");
  TH1D  *FONLLCurve = (TH1D *)FONLL->Get("FONLL");

  TH1D *Weight = (TH1D *)FONLLCurve->Clone("FONLLWeights");
  Weight->Divide(PYTHIAPt);

  cout << Weight->Integral() << endl;

  TH1D *WeightFromFit = (TH1D *)FitHist->Clone("FitWeights");
  WeightFromFit->Divide(PYTHIAPt);

  // Aug14_FONLL_1GeV_Data_WeighedByData4/GEANTvDATA_1.root
  TFile *DataWeightFile = new TFile(Form("Aug14_FONLL_%iGeV_Data_WeighedByData4/GEANTvDATA_%i.root", (int)lowptcutoff, (int)lowptcutoff), "READ");
  DataWeightFile->cd();
  TH2D *DataWeight[3];

  for (int i = 0; i < 3; i++){
    DataWeight[i] = (TH2D *)gDirectory->Get(Form("DataWeight_%i", i));
    DataWeight[i]->SetDirectory(0);
  }

  for (int i = 0; i < 5000000; i++){
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

    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);

    double RecoJetPxWide = RecoJetPtFromArea *TMath::Cos(RecoJetPhi);
    double RecoJetPyWide = RecoJetPtFromArea*TMath::Sin(RecoJetPhi);
    double recozwide = (RecoJetPxWide*RecoD0Px + RecoJetPyWide*RecoD0Py)/pow(RecoJetPtFromArea, 2);

    if (mcz >=1) mcz = 0.9999;
    
    bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue;

    double weight = Weight->GetBinContent(Weight->FindBin(MCJetPt));

    double weightfromfit = WeightFromFit->GetBinContent(WeightFromFit->FindBin(MCJetPt));

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

    if (centhistogramtofill < 0) continue;

    FONLLPt->Fill(MCJetPt, weight);
    if (MCJetPt > 5) FONLLZ->Fill(mcz, weight);
    FONLLWeighted->Fill(MCJetPt, mcz, weight);

    FitPt->Fill(MCJetPt, weightfromfit);
    if (MCJetPt > 5) FitZ->Fill(mcz, weightfromfit);
    FitWeighted->Fill(MCJetPt, mcz, weightfromfit);

    int measuredptbin = DataWeight[centhistogramtofill]->GetXaxis()->FindBin(RecoJetPtFromArea);
    int measuredzbin = DataWeight[centhistogramtofill]->GetYaxis()->FindBin(recozwide);
    double w3 = DataWeight[centhistogramtofill]->GetBinContent(measuredptbin, measuredzbin);

    DataWeighedPt[centhistogramtofill]->Fill(MCJetPt, weight*w3);
    if (MCJetPt > 5) DataWeighedZ[centhistogramtofill]->Fill(mcz, weight*w3);
    DataWeighed[centhistogramtofill]->Fill(MCJetPt, mcz, weight*w3);

  }
  cout << endl;

  TH1D *WeightFromDataFit[3];
  for (int i = 0; i < 3; i++){
    WeightFromDataFit[i] = (TH1D *)DataWeighedPt[i]->Clone(Form("DataWeighedPt_%i", i));
    SetName(WeightFromDataFit[i], Form("DataWeighed_PtWeight_%i", i));
    WeightFromDataFit[i]->Divide(PYTHIAPt);
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

    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);

    double RecoJetPxWide = RecoJetPtFromArea *TMath::Cos(RecoJetPhi);
    double RecoJetPyWide = RecoJetPtFromArea*TMath::Sin(RecoJetPhi);
    double recozwide = (RecoJetPxWide*RecoD0Px + RecoJetPyWide*RecoD0Py)/pow(RecoJetPtFromArea, 2);

    if (mcz >=1) mcz = 0.9999;
    
    bool isMCZ = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue;

    double weight = Weight->GetBinContent(Weight->FindBin(MCJetPt));

    double weightfromfit = WeightFromFit->GetBinContent(WeightFromFit->FindBin(MCJetPt));

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

    if (centhistogramtofill < 0) continue;

    DataWeighedPtFit[centhistogramtofill]->Fill(MCJetPt, WeightFromDataFit[centhistogramtofill]->GetBinContent(WeightFromDataFit[centhistogramtofill]->FindBin(MCJetPt)));
    if (MCJetPt > 5)DataWeighedZFit[centhistogramtofill]->Fill(mcz, WeightFromDataFit[centhistogramtofill]->GetBinContent(WeightFromDataFit[centhistogramtofill]->FindBin(MCJetPt)));

  }
  cout << endl;



  int lowbin = PYTHIAPt->GetXaxis()->FindBin(5.0);
  int highbin = PYTHIAPt->GetXaxis()->FindBin(30.0);
  double integral = PYTHIAPt->Integral(lowbin, highbin);

  // cout << FONLLPt->Integral(lowbin, highbin) << "\t" << FitPt->Integral(lowbin, highbin) << "\t" << DataWeighedPt->Integral(lowbin, highbin) << endl;

  FONLLPt->Scale(integral/FONLLPt->Integral(lowbin, highbin));
  FONLLZ->Scale(integral/FONLLZ->Integral());
  FONLLWeighted->Scale(1./FONLLWeighted->Integral());

  FitPt->Scale(integral/FitPt->Integral(lowbin, highbin));
  FitZ->Scale(integral/FitZ->Integral());
  FitWeighted->Scale(1./FitWeighted->Integral());

  for (int cent = 0; cent < 3; cent++){
    DataWeighedPt[cent]->Scale(integral/DataWeighedPt[cent]->Integral(lowbin, highbin));
    DataWeighedPtFit[cent]->Scale(integral/DataWeighedPtFit[cent]->Integral(lowbin, highbin));
    DataWeighedZ[cent]->Scale(integral/DataWeighedZ[cent]->Integral());
    DataWeighedZFit[cent]->Scale(integral/DataWeighedZFit[cent]->Integral());
    DataWeighed[cent]->Scale(1./DataWeighed[cent]->Integral());
  }
  
  // FONLLWeighted->Fit()

  TLegend *legend2 = new TLegend(0.5,0.6,0.9,0.9);
  

  TCanvas *c = new TCanvas("c", "c", 1800, 900);
  c->Divide(3);
  for (int cent = 0; cent < 3; cent++){
    c->cd(cent+1);
    gPad->SetLogy();
    PYTHIAPt->GetYaxis()->SetRangeUser(pow(10, -8), 1);
    PYTHIAPt->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
    PYTHIAPt->GetYaxis()->SetTitle("a.u.");
    SetColor(FONLLPt, kRed, 20);
    SetColor(FitPt, kGreen-2, 20);
    SetColor(DataWeighedPt[cent], kBlack, 24);
    SetColor(DataWeighedPtFit[cent], kBlue, 22);
    FONLLPt->Draw("EP SAME");
    FitPt->Draw("EP SAME");
    DataWeighedPt[cent]->Draw("EP SAME");
    // DataWeighedPtFit[cent]->Draw("EP SAME");
    PYTHIAPt->Draw("HIST SAME");
    if (cent == 0){
      legend2->AddEntry(PYTHIAPt, "PYTHIA", "lp");
      legend2->AddEntry(FONLLPt, "FONLL", "lp");
      legend2->AddEntry(FitPt, "Power Law Fit", "lp");
      legend2->AddEntry(DataWeighedPt[cent], "Data Weighted", "lp");
      // legend2->AddEntry(DataWeighedPtFit[cent], "Data Weighted Fit", "lp");
    }
    legend2->Draw("SAME");
  }

  c->SaveAs(Form("PYTHIA_Pt_D0pT_%i.pdf", (int)lowptcutoff));

  TLegend *legend3 = new TLegend(0.6,0.1,0.9,0.4);

  TCanvas *c2 = new TCanvas("c2", "c2", 1800, 900);
  c2->Divide(3);
  for (int cent = 0; cent < 3; cent++){
    c2->cd(cent+1);
    gPad->SetLogy();
    PYTHIAZ->GetYaxis()->SetRangeUser(pow(10, -4), pow(10, -2));
    PYTHIAZ->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
    PYTHIAZ->GetYaxis()->SetTitle("a.u.");
    SetColor(FONLLZ, kRed, 20);
    SetColor(FitZ, kGreen-2, 20);
    SetColor(DataWeighedZ[cent], kBlack, 24);
    SetColor(DataWeighedZFit[cent], kBlue, 22);
    // PYTHIAZ->Draw("HIST SAME");
    FONLLZ->Draw("EP SAME");
    FitZ->Draw("EP SAME");
    DataWeighedZ[cent]->Draw("EP SAME");
    // DataWeighedZFit[cent]->Draw("EP SAME");
    
    if (cent == 0){
      legend3->AddEntry(PYTHIAZ, "PYTHIA", "lp");
      legend3->AddEntry(FONLLZ, "FONLL", "lp");
      legend3->AddEntry(FitZ, "Power Law Fit", "lp");
      legend3->AddEntry(DataWeighedZ[cent], "Data Weighted", "lp");
      // legend3->AddEntry(DataWeighedZFit[cent], "Data Weighted Fit", "lp");
    }
    legend3->Draw("SAME");
  }

  c2->SaveAs(Form("PYTHIA_Z_D0pT_%i.pdf", (int)lowptcutoff));
  // c->cd(2);
  // gPad->SetLogy();
  // PYTHIAZ->Draw("HIST");
  // PYTHIAZ->GetXaxis()->SetTitle("z");
  // PYTHIAZ->GetYaxis()->SetTitle("a.u.");
  // FONLLZ->Draw("SAME");
  // FitZ->Draw("SAME");
  // DataWeighedZ->Draw("SAME");
  // SetColor(FONLLZ, kRed);
  // SetColor(FitZ, kGreen-2);
  // SetColor(DataWeighedZ, kBlack);
  // c->cd(3);
  // gPad->SetLogy();
  // Weight->Draw();
  // Weight->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
  // WeightFromFit->Draw("SAME");
  // SetColor(Weight, kRed);
  // SetColor(WeightFromFit, kGreen-2);

  c->SaveAs("PYTHIA_Weighted.pdf");

  TCanvas *d = new TCanvas("d", "d", 1800, 900);
  d->Divide(4);
  d->cd(1);
  gPad->SetLogz();
  PYTHIA->Draw("COLZ");
  PYTHIA->GetZaxis()->SetRangeUser(pow(10, -6),1);
  d->cd(2);
  gPad->SetLogz();
  FONLLWeighted->Draw("COLZ");
  FONLLWeighted->GetZaxis()->SetRangeUser(pow(10, -6),1);
  d->cd(3);
  gPad->SetLogz();
  FitWeighted->Draw("COLZ");
  FitWeighted->GetZaxis()->SetRangeUser(pow(10, -6),1);
  d->cd(4);
  gPad->SetLogz();
  DataWeighed[0]->Draw("COLZ");
  DataWeighed[0]->GetZaxis()->SetRangeUser(pow(10, -6),1);
  

  TFile *g = new TFile(Form("PYTHIA_Pt_%i_%i.root", (int)lowptcutoff, (int)highptcutoff), "RECREATE");
  g->cd();
  PYTHIAPt->Write();
  PYTHIAZ->Write();
  PYTHIA->Write();
  Weight->Write();
  WeightFromFit->Write();
  WeightFromDataFit[0]->Write();
  WeightFromDataFit[1]->Write();
  WeightFromDataFit[2]->Write();
  g->Close();

}