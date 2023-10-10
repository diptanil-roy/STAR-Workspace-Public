R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int D0pTLow = 1, int mode = 0, bool isForCentralityWeight = false, int priormode = 1){
  // PriorMode
  // 0 -> FONLL (default)
  // 1 -> FONLL + Reweighed with data distros
  // 2 -> PYTHIA
  // 3 -> HERWIG

  cout << D0pTLow << "\t" << mode  << "\t" << isForCentralityWeight << "\t" << priormode << endl; 

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  // TFile *f = new TFile("HIOverlay_CS_May20.root");
  // TFile *f = new TFile("HIOverlay_HFJets_WithCS_Jun6_2023_pthat_3_inf.root");
  TFile *f = new TFile("HIOverlay_HFJets_WithCS_Jun20_2023_pthat_3_inf.root");
  // TFile *f = new TFile("/Volumes/WorkDrive/Jan26_2023_PYTHIA/HIResponse_MCIncluded_Feb6.root");

  cout << f->GetName() <<  endl;
  f->cd("HIJetSaver");

  TTree *JetTree = (TTree *)gDirectory->Get("Jets");
  // TTree *RecoJetTree = (TTree *)gDirectory->Get("RecoJets");

  Float_t         Centrality;
  Float_t         gRefMult;
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
  Float_t         RecoJetPtFromArea;
  Int_t           RecoJetNConstFromArea;

  JetTree->SetBranchStatus("RecoJet*FromPYTHIA*", '0');

  JetTree->SetBranchAddress("Centrality", &Centrality);
  JetTree->SetBranchAddress("gRefMult", &gRefMult);
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

  JetTree->SetBranchAddress("RecoJetPtCS", &RecoJetCorrPt);
  JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
  JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
  // JetTree->SetBranchAddress("RecoJetAreaCS", &RecoJetArea);
  JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
  JetTree->SetBranchAddress("RecoJetECS", &RecoJetE);
  // JetTree->SetBranchAddress("RecoJetRhoValCS", &RecoJetRhoVal);
  JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
  JetTree->SetBranchAddress("RecoJetNConstCS", &RecoJetNConst);

  // JetTree->SetBranchAddress("RecoJetCorrPt", &RecoJetPtFromArea);
  JetTree->SetBranchAddress("RecoJetPt", &RecoJetPtFromArea);
  // JetTree->SetBranchAddress("RecoJetCorrPt", &RecoJetCorrPt);
  // JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
  // JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
  // JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
  // JetTree->SetBranchAddress("RecoJetE", &RecoJetE);
  // JetTree->SetBranchAddress("RecoJetRhoVal", &RecoJetRhoVal);
  JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConstFromArea);

  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  TString FONLLFileName = Form("FONLL_Pt_%i_30.root", 1);
  TString PYTHIAFileName = Form("PYTHIA_Pt_%i_30.root", 1);

  TFile* FONLL = new TFile(FONLLFileName.Data(),"READ");
  // TFile* FONLL = new TFile("FONLL_Pt_5_20.root","READ");

  TH1F * FONLLCurve = (TH1F*) FONLL->Get("FONLL");

  TFile *PYTHIA = new TFile(PYTHIAFileName.Data(), "READ");
  // TFile *PYTHIA = new TFile("PYTHIA_Pt_5_20.root", "READ");

  TH1F  *PYTHIAPtCurve = (TH1F *) PYTHIA->Get("PYTHIA pT");
  TH1F  *PYTHIAZCurve  = (TH1F *) PYTHIA->Get("PYTHIA Z");
  TH2F  *PYTHIA2D      = (TH2F *) PYTHIA->Get("PYTHIA");
  TH1F  *FONLLvPYTHIAWeights = (TH1F *) PYTHIA->Get("FONLLWeights");
  cout << "FONLLvPYTHIAWeights->GetNbinsX() = " << FONLLvPYTHIAWeights->GetNbinsX() << endl;

  TH1F * CENTWeight = NULL;
  TH1F * gREFMULTWeight = NULL;

  if (!isForCentralityWeight){
    TString CentWeightFileName = Form("CentWeight_D0pT_%iGeV/CentWeight.root", D0pTLow);
    TFile* CentWeights = new TFile(CentWeightFileName.Data(), "READ");
    CENTWeight = (TH1F *)CentWeights->Get("Centrality Weights");
    CENTWeight->SetDirectory(0);
    gREFMULTWeight = (TH1F *)CentWeights->Get("Centrality Weights From RefCorr");
    gREFMULTWeight->SetDirectory(0);
  }

  TH1D *DataWeightPt[3];
  TH1D *DataWeightZ[3];
  TH2D *DataWeight[3];

  if (priormode == 1){
    TFile *DataWeightFile = new TFile(Form("%s/GEANTvDATA_%i.root", DirName.Data(), D0pTLow), "READ");
    DataWeightFile->cd();
    for (int i = 0; i < 3; i++){
      DataWeightPt[i] = (TH1D *)gDirectory->Get(Form("DataWeightPt_%i", i));
      DataWeightZ[i] = (TH1D *)gDirectory->Get(Form("DataWeightZ_%i", i));
      DataWeight[i] = (TH2D *)gDirectory->Get(Form("DataWeight_%i", i));
      DataWeightPt[i]->SetDirectory(0);
      DataWeightZ[i]->SetDirectory(0);
      DataWeight[i]->SetDirectory(0);
    }
  }

  ///// Everything we need to set up response matrices is imported before this line. /////

  TH1D *hZ[3];
  TH2D *hJetPtvZ[3];
  THnSparseF *hJetWide[3];
  THnSparseF *hJet[3];
  THnSparseF *hJ[3];

  TH2D *Weight[3];
  RooUnfoldResponse *resp[3]; //Response
  RooUnfoldResponse *respwide[3]; //Response
  RooUnfoldResponse *resp1D[3]; //Response
  // RooUnfoldResponse *respdR[3]; //Response

  TH2D *fMeas[3];
  TH2D *fTrue[3];
  TH2D *fMiss[3];
  TH2D *fFake[3];

  TH2D *fMeasWide[3];
  TH2D *fFakeWide[3];

  TH2D *fMissZNorm[3];
  TH2D *fMissZNormPtNorm[3];

  TH2D *fMeas2D[3];
  TH1D *fMeas1D[3];
  TH1D *fTrue1D[3];
  TH1D *fMiss1D[3];
  TH1D *fFake1D[3];

  TH1D *fMeasZ1D[3];
  TH1D *fTrueZ1D[3];
  TH1D *fMissZ1D[3];

  // TH2D *hRecoJetPtvsMCJetPt[3];
  TH1D *hDiffJetPt[3]; //Reco - MC
  TH1D *hDiffJetPtInGenPtBins[3][njpt_gen_bins_var]; //Reco - MC
  TH1D *hDiffJetEta[3]; //Reco - MC
  TH1D *hDiffJetPhi[3]; //Reco - MC

  TH1D *hMCD0Pt[3]; //Reco - MC
  TH1D *hRecoD0Pt[3]; //Reco - MC
  // TH1D *hDiffJetPhi[3]; //Reco - MC

  TH2D *hJetPtAreaVsJetPtCS[3];
  TH2D *hJetAreaVsJetNConstArea[3];
  TH2D *hJetAreaVsJetNConstCS[3];
  TH2D *hJetPtAreaVsJetNConstArea[3];
  TH2D *hJetPtCSVsJetNConstCS[3];

  TH2D *hGendRvJetPt[3];
  TH1D *hGendR[3];

  TH2D *hdRvJetPt[3];
  TH1D *hdR[3];

  TH1D *hCentrality;

  TH1D *hgRefMultCorr;

  TH1D *hRecoD0Mass;

  TH2D *hCentVsJetNConst;

  // THnSparseF *masterTHn = new THnSparseF("masterTHn", "masterTHn", ndim, )

  cout << JetTree->GetEntries() << endl;

  for (int i = 0; i < 3; i++){
    hZ[i] = new TH1D(Form("hZ_%i", i), Form("hZ_%i", i), nz_gen_bins, z_gen_bin);
    hJetPtvZ[i] = new TH2D(Form("hJetPtvZ_%i", i), Form("hJetPtvZ_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);

    hJ[i] = new THnSparseF(Form("hJ_%i", i), Form("hJ_%i", i), ndim, nbins, NULL, NULL);
    hJet[i] = new THnSparseF(Form("hJet_%i", i), Form("hJet_%i", i), ndim, nbins, NULL, NULL);
    hJetWide[i] = new THnSparseF(Form("hJetWide_%i", i), Form("hJetWide_%i", i), ndim, nbins_wide, NULL, NULL);

    hJ[i]->SetBinEdges(0, jetpt_var_bin);
    hJ[i]->SetBinEdges(1, z_gen_bin);
    hJ[i]->SetBinEdges(2, nbinsjetpt);
    hJ[i]->SetBinEdges(3, nbinsz);

    hJet[i]->SetBinEdges(0, jetpt_var_bin);
    hJet[i]->SetBinEdges(1, z_gen_bin);
    hJet[i]->SetBinEdges(2, nbinsjetpt);
    hJet[i]->SetBinEdges(3, nbinsz);

    hJetWide[i]->SetBinEdges(0, jetpt_var_bin);
    hJetWide[i]->SetBinEdges(1, z_gen_bin);
    hJetWide[i]->SetBinEdges(2, nbinsjetpt_wide);
    hJetWide[i]->SetBinEdges(3, nbinsz_wide);

    fMeas[i] = new TH2D(Form("fMeas_Cent_%i", i), Form("fMeas_Cent_%i", i), njpt_bins, nbinsjetpt, nz_bins, nbinsz);
    fTrue[i] = new TH2D(Form("fTrue_Cent_%i", i), Form("fTrue_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
    fMiss[i] = new TH2D(Form("fMiss_Cent_%i", i), Form("fMiss_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
    fFake[i] = new TH2D(Form("fFake_Cent_%i", i), Form("fFake_Cent_%i", i), njpt_bins, nbinsjetpt, nz_bins, nbinsz);

    fMeasWide[i] = new TH2D(Form("fMeasWide_Cent_%i", i), Form("fMeasWide_Cent_%i", i), njpt_bins_wide, nbinsjetpt_wide, nz_bins_wide, nbinsz_wide);
    fFakeWide[i] = new TH2D(Form("fFakeWide_Cent_%i", i), Form("fFakeWide_Cent_%i", i), njpt_bins_wide, nbinsjetpt_wide, nz_bins_wide, nbinsz_wide);

    fMeas1D[i] = new TH1D(Form("fMeas1D_Cent_%i", i), Form("fMeas1D_Cent_%i", i), njpt_bins, nbinsjetpt);
    fTrue1D[i] = new TH1D(Form("fTrue1D_Cent_%i", i), Form("fTrue1D_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin);
    fMiss1D[i] = new TH1D(Form("fMiss1D_Cent_%i", i), Form("fMiss1D_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin);
    fFake1D[i] = new TH1D(Form("fFake1D_Cent_%i", i), Form("fFake1D_Cent_%i", i), njpt_bins, nbinsjetpt);

    // cout << "Here" << endl;
    resp[i] = new RooUnfoldResponse(Form("Resp_%i", i), Form("Resp_%i", i));
    resp[i]->Setup(fMeas[i], fTrue[i]); //Setup Response Matrix Definition

    respwide[i] = new RooUnfoldResponse(Form("RespWide_%i", i), Form("RespWide_%i", i));
    respwide[i]->Setup(fMeasWide[i], fTrue[i]); //Setup Response Matrix Definition

    resp1D[i] = new RooUnfoldResponse(Form("Resp1D_%i", i), Form("Resp1D_%i", i));
    resp1D[i]->Setup(fMeas1D[i], fTrue1D[i]); //Setup Response Matrix Definition

    fMissZNorm[i] = new TH2D(Form("fMissZNorm_Cent_%i", i), Form("fMissZNorm_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
    fMissZNormPtNorm[i] = new TH2D(Form("fMissZNormPtNorm_Cent_%i", i), Form("fMissZNormPtNorm_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
    
    fMeas2D[i] = new TH2D(Form("fMeas2D_Cent_%i", i), Form("fMeas2D_Cent_%i", i), njpt_bins, nbinsjetpt, 100, 0, 10);

    fMeasZ1D[i] = new TH1D(Form("fMeasZ1D_Cent_%i", i), Form("fMeasZ1D_Cent_%i", i), 100, 0, 10);
    fTrueZ1D[i] = new TH1D(Form("fTrueZ1D_Cent_%i", i), Form("fTrueZ1D_Cent_%i", i), 100, 0, 10);
    fMissZ1D[i] = new TH1D(Form("fMissZ1D_Cent_%i", i), Form("fMissZ1D_Cent_%i", i), 100, 0, 10);

    hDiffJetPt[i] = new TH1D(Form("hDiffJetPt_%i", i), Form("hDiffJetPt_%i", i), 100, -50, 50);
    hDiffJetEta[i] = new TH1D(Form("hDiffJetEta_%i", i), Form("hDiffJetEta_%i", i), 40, -1., 1.);
    hDiffJetPhi[i] = new TH1D(Form("hDiffJetPhi_%i", i), Form("hDiffJetPhi_%i", i), 40, -1., 1.);

    for (int ptbin = 0; ptbin < njpt_gen_bins_var; ptbin++){
      hDiffJetPtInGenPtBins[i][ptbin] = new TH1D(Form("hDiffJetPt_%i_%i", i, ptbin), Form("hDiffJetPt_%i_%i", i, ptbin), 100, -50, 50);
    }

    hJetPtAreaVsJetPtCS[i] = new TH2D(Form("hJetPtAreaVsJetPtCS_%i", i), Form("hJetPtAreaVsJetPtCS_%i", i), 38, -25.5, 50.5, 38, -25.5, 50.5);

    hJetAreaVsJetNConstArea[i] = new TH2D(Form("hJetAreaVsJetNConstArea_%i", i), Form("hJetAreaVsJetNConstArea_%i", i), 100, 0., 1., 151, -0.5, 150.5);
    hJetAreaVsJetNConstCS[i] = new TH2D(Form("hJetAreaVsJetNConstCS_%i", i), Form("hJetAreaVsJetNConstCS_%i", i), 100, 0., 1., 151, -0.5, 150.5);

    hJetPtAreaVsJetNConstArea[i] = new TH2D(Form("hJetPtAreaVsJetNConstArea_%i", i), Form("hJetPtAreaVsJetNConstArea_%i", i), njpt_bins_wide, nbinsjetpt_wide, 151, -0.5, 150.5);
    hJetPtCSVsJetNConstCS[i] = new TH2D(Form("hJetPtCSVsJetNConstCS_%i", i), Form("hJetPtCSVsJetNConstCS_%i", i), njpt_bins, nbinsjetpt, 151, -0.5, 150.5);

    hMCD0Pt[i] = new TH1D(Form("hMCD0Pt_%i", i), Form("hMCD0Pt_%i", i), 20, 0, 10);
    hRecoD0Pt[i] = new TH1D(Form("hRecoD0Pt_%i", i), Form("hRecoD0Pt_%i", i), 20, 0, 10);

    hGendRvJetPt[i] = new TH2D(Form("hGendRvJetPt_%i", i), Form("hGendRvJetPt_%i", i), ndrbins, drbins, njpt_gen_bins_var, jetpt_var_bin);
    hGendR[i] = new TH1D(Form("hGendR_%i", i), Form("hGendR_%i", i), ndrbins, drbins);

    hdRvJetPt[i] = new TH2D(Form("hdRvJetPt_%i", i), Form("hdRvJetPt_%i", i), ndrbins, drbins, njpt_bins_wide, nbinsjetpt_wide);
    hdR[i] = new TH1D(Form("hdR_%i", i), Form("hdR_%i", i), ndrbins, drbins);

    hJet[i]->Sumw2();
    hJet[i]->CalculateErrors();

    hJ[i]->Sumw2();
    hJ[i]->CalculateErrors();
  }

hCentrality = new TH1D("hCentrality","hCentrality", ncentbin, centbin);

hgRefMultCorr = new TH1D("hgRefMultCorr","hgRefMultCorr", 100, 0, 800);

hRecoD0Mass = new TH1D("Reco D0 Mass", "Reco D0 Mass", 1000, 1, 3);

hCentVsJetNConst = new TH2D("Cent vs Jet Constituents", "Cent vs Jet Constituents", ncentbin, centbin, njetconstbin, jetconstbin);

cout << JetTree->GetEntries() << endl;

int nentries = JetTree->GetEntries();

int lowlimit = 0;
int highlimit = nentries;

int counter = 0;

TH1D *plottingpT = new TH1D("plottingpT", "plottingpT", njpt_gen_bins_var, jetpt_var_bin);

cout << "Mode == " << mode << endl;

cout << "Reading events from " << lowlimit << " to " << highlimit << endl;

cout << "================ Limits =================" << endl;
cout << Form("%i < pT,D0 [GeV/c] < %i", D0pTLow, 10) << endl;
cout << Form("%.1f < MC pT,Jet [GeV/c] < %.1f", 5., jetpt_var_bin[njpt_gen_bins_var]) << endl;
cout << Form("%.1f < Reco pT,Jet [GeV/c] < %.1f", (float)D0pTLow, nbinsjetpt[njpt_bins]) << endl;
cout << Form("%.1f < MC Z,Jet < %.1f", z_gen_bin[0], z_gen_bin[nz_gen_bins]) << endl;
cout << Form("%.1f < Reco Z,Jet < %.1f", nbinsz[0], nbinsz[nz_bins]) << endl;
cout << "=========================================" << endl;

if (countcutoff > 0){
  cout << "Count cutoff is " << countcutoff << endl;
  for (int i = lowlimit; i < highlimit; i++){

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

      bool isMCD0Pt = MCD0Pt > D0pTLow && MCD0Pt < 10;
      bool isRecoD0Pt = RecoD0Pt > D0pTLow && RecoD0Pt < 10;
      bool isMCJetPt = MCJetPt > 5 && MCJetPt < jetpt_var_bin[njpt_gen_bins_var];

      if (!isMCD0Pt) continue;
      if (!isRecoD0Pt) continue;
      if (!isMCJetPt) continue;
      if (RecoJetNConst==0) continue;

      double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
      double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
      double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
      double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
      double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);

      // if (mcz < 0.3) cout << MCJetPt << "\t" << MCD0Pt << mcz << endl;
      if (mcz >= 1.0) mcz = 0.999; // Padding the boundaries

      bool isMCZ   = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
      if (!isMCZ) continue; 

      int genptbin = plottingpT->FindBin(MCJetPt);

      if(RecoJetPhi>=2.*TMath::Pi()) RecoJetPhi = RecoJetPhi - 2.*TMath::Pi();
      if(RecoJetPhi<0) RecoJetPhi = 2.*TMath::Pi() + RecoJetPhi;

      // if (RecoJetCorrPt < RecoD0Pt) RecoJetCorrPt = RecoD0Pt;

      double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
      double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
      double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
      double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
      double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);

      TLorentzVector RecoPion;
      TLorentzVector RecoKaon;
      TLorentzVector RecoD0;

      RecoPion.SetPtEtaPhiM(RecoPionPt, RecoPionEta, RecoPionPhi, 0.13957);
      RecoKaon.SetPtEtaPhiM(RecoKaonPt, RecoKaonEta, RecoKaonPhi, 0.49368);
      RecoD0 = RecoPion + RecoKaon;

      hRecoD0Mass->Fill(RecoD0.M());

      // if (RecoD0.M() < 1.75 || RecoD0.M() > 2.02) continue;

      if (recoz >= 1.0) recoz = 0.999; // Padding the boundaries

      // cout << recoz << endl;
      double RecoDeltaR = dR(dEta(RecoJetEta, RecoD0Eta), dPhi(RecoJetPhi, RecoD0Phi));

      // if (RecoJetCorrPt < nbinsjetpt[0] || RecoJetCorrPt > nbinsjetpt[njpt_bins] || recoz < nbinsz[0] || recoz > nbinsz[nz_bins]) cout << "RecoJetCorrPt = " << RecoJetCorrPt << "\t" << "RecoD0Pt = " << RecoD0Pt << "\t" << "recoz = " << recoz << endl;

      //Defining response matrix using fakes and misses to check closure
      bool isRecoJetPt = RecoJetCorrPt >= nbinsjetpt[0] && RecoJetCorrPt <= nbinsjetpt[njpt_bins];
      // bool isRecoJetPt = RecoJetPtFromArea > 3 && RecoJetPtFromArea < 30;
      bool isMCJetEta = abs(MCJetEta) < 0.601;
      bool isRecoJetEta = abs(RecoJetEta) < 0.6;
      bool isRecoDeltaR = RecoDeltaR < 0.4;

      int centbin = CENTWeight->FindBin(Centrality);

      double wcent = CENTWeight->GetBinContent(centbin);
      double fonllweight = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));

      if (mcz > z_gen_bin[nz_gen_bins]) mcz = z_gen_bin[nz_gen_bins] + z_gen_binwidth_high*0.5;
      if (mcz < z_gen_bin[0]) mcz = z_gen_bin[0] - z_gen_binwidth_low*0.5;
      if (mcz < mczcutoff) continue;

      // if (!isRecoJetPt) continue;

      bool isOK   = isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta;
      bool isMISS = !isOK && isMCJetPt && isMCJetEta;
      bool isFAKE = !isOK && isRecoJetPt && isRecoJetEta;

      double w = wcent*fonllweight;
      // double w = 1.;
      
      if (!isOK && !isMISS && !isFAKE) continue;

      if (recoz > nbinsz[nz_bins]) {
        cout << recoz << endl;
        recoz = nbinsz[nz_bins] + z_binwidth_high*0.5;
      }
      if (recoz < nbinsz[0]) {
        cout << recoz << endl;
        recoz = nbinsz[0] - z_binwidth_low*0.5;
      }

      if (RecoJetCorrPt < nbinsjetpt[0]) {
        cout << RecoJetCorrPt << endl;
        RecoJetCorrPt = nbinsjetpt[0] - jetpt_binwidth_low*0.5;
      }    
      if (RecoJetCorrPt > nbinsjetpt[njpt_bins]) {
        cout << RecoJetCorrPt << endl;
        RecoJetCorrPt = nbinsjetpt[njpt_bins] + jetpt_binwidth_high*0.5;
      }
      bool isRecoZ = recoz > nbinsz[0] && recoz < nbinsz[nz_bins];

      if (isOK) {
          double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
          hJ[centhistogramtofill]->Fill(fvalue, 1); //Main body of the response matrix is the only thing I fill in the main THnSparses. The misses are accounted for separately.
      }

      else if (isMISS) {

          double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};

          if (!isRecoJetEta) {
            if (abs(RecoJetCorrPt - nbinsjetpt[0]) < abs(RecoJetCorrPt - nbinsjetpt[njpt_bins])){
              fvalue[2] = nbinsjetpt[0] - jetpt_binwidth_low*0.5;
              if (abs(recoz - nbinsz[0]) < abs(recoz - nbinsz[nz_bins])) fvalue[3] = nbinsz[0] - z_binwidth_low*0.5;
              else fvalue[3] = nbinsz[nz_bins] + z_binwidth_high*0.5;
            }
            else{
              fvalue[2] = nbinsjetpt[njpt_bins] + jetpt_binwidth_high*0.5;
              if (abs(recoz - nbinsz[0]) < abs(recoz - nbinsz[nz_bins])) fvalue[3] = nbinsz[0] - z_binwidth_low*0.5;
              else fvalue[3] = nbinsz[nz_bins] + z_binwidth_high*0.5;
            }
          }

          hJ[centhistogramtofill]->Fill(fvalue, 1); //Misses of the response matrix. Let's see if this works.
      }
  }

}
cout << "Low D0 pT = " << D0pTLow << " GeV" << endl;

double okcount[3] = {0};
double misscount[3] = {0};
double fakecount[3] = {0};
double singlejetcount[3] = {0};
double singlemcjetcount[3] = {0};

for (int i = lowlimit; i < highlimit; i++){

    if (mode == 1) {if (i%2 != 0) continue;}
    else if (mode == 2) {if (i%2 == 0) continue;}

    if (mode == 1) {if (i%1000000 == 0) cout << "Read Entry " << i << "\r" << flush;}
    else {if (i%1000000 == 1) cout << "Read Entry " << i << "\r" << flush;}

    // if (mode == 1) {if (i%1000000 == 0) cout << "Read Entry " << i << "\n";}
    // else {if (i%1000000 == 1) cout << "Read Entry " << i << "\n";}

    // cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    int centhistogramtofill = -99;
    if (Centrality < 10) centhistogramtofill = 0;
    else if (Centrality >= 10 && Centrality < 40) centhistogramtofill = 1;
    else if (Centrality >= 40 && Centrality <= 80) centhistogramtofill = 2;

    if (centhistogramtofill < 0) continue;

    // if (MCPrimaryVertex->size() < 3 || RecoPrimaryVertex->size() < 3) continue;

    if (MCJetPt < MCD0Pt && MCJetNConst > 0) MCJetPt = MCD0Pt;

    bool isMCD0Pt = MCD0Pt > D0pTLow && MCD0Pt < 10; // We wanna check 1 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > D0pTLow && RecoD0Pt < 10;
    bool isMCJetPt = MCJetPt > 5. && MCJetPt < jetpt_var_bin[njpt_gen_bins_var];

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    // if (!isMCJetPt) continue;
    if (RecoJetNConst == 0) continue;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);
    double MCDeltaR = dR(dEta(MCJetEta, MCD0Eta), dPhi(MCJetPhi, MCD0Phi));
    // double mcz = MCD0Pt/MCJetPt;

    // if (mcz > 1.0) cout << MCJetPt << "\t" << MCD0Pt << mcz << endl;
    if (mcz >= 1.0) mcz = 0.999; // Padding the boundaries

    bool isMCZ   = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    if (!isMCZ) continue; 

    int genptbin = plottingpT->FindBin(MCJetPt);

    if(RecoJetPhi>=2.*TMath::Pi()) RecoJetPhi = RecoJetPhi - 2.*TMath::Pi();
    if(RecoJetPhi<0) RecoJetPhi = 2.*TMath::Pi() + RecoJetPhi;

    // if (RecoJetCorrPt < RecoD0Pt) RecoJetCorrPt = RecoD0Pt; // This is to avoid negative values in the response matrix
    // if (RecoJetCorrPt > 5 && RecoJetPtFromArea > RecoD0Pt) RecoJetCorrPt = RecoJetPtFromArea; // Just to make sure the differences are artifacts of differences of pT distributions
    
    if (RecoJetCorrPt < RecoD0Pt) RecoJetCorrPt = RecoD0Pt;
    // if (RecoJetCorrPt > 9.) RecoJetCorrPt = RecoJetPtFromArea;

    // if (RecoJetPtFromArea < nbinsjetpt_wide[0] && RecoJetNConst != 1)

    double RecoJetPx = RecoJetCorrPt *TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);
    // double recoz = RecoD0Pt/RecoJetCorrPt;
    double RecoDeltaR = dR(dEta(RecoJetEta, RecoD0Eta), dPhi(RecoJetPhi, RecoD0Phi));

    double RecoJetPxWide = RecoJetPtFromArea *TMath::Cos(RecoJetPhi);
    double RecoJetPyWide = RecoJetPtFromArea*TMath::Sin(RecoJetPhi);
    double recozwide = (RecoJetPxWide*RecoD0Px + RecoJetPyWide*RecoD0Py)/pow(RecoJetPtFromArea, 2);

    TLorentzVector RecoPion;
    TLorentzVector RecoKaon;
    TLorentzVector RecoD0;

    RecoPion.SetPtEtaPhiM(RecoPionPt, RecoPionEta, RecoPionPhi, 0.13957);
    RecoKaon.SetPtEtaPhiM(RecoKaonPt, RecoKaonEta, RecoKaonPhi, 0.49368);
    RecoD0 = RecoPion + RecoKaon;

    // if (RecoD0.M() < 1.75 || RecoD0.M() > 2.02) continue;

    if (recoz >= 1.0) recoz = 0.999; // Padding the boundaries
  
    //Defining response matrix using fakes and misses to check closure
    bool isRecoJetPt = RecoJetCorrPt >= D0pTLow && RecoJetCorrPt <= nbinsjetpt[njpt_bins];
    bool isRecoJetPtWide = RecoJetPtFromArea >= nbinsjetpt_wide[0] && RecoJetPtFromArea <= nbinsjetpt_wide[njpt_bins_wide];
    // bool isRecoJetPt = RecoJetPtFromArea > -25 && RecoJetPtFromArea < 50;
    bool isMCJetEta = abs(MCJetEta) < 0.601;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoDeltaR = RecoDeltaR < 0.4;
    bool isRecoArea = RecoJetArea > 0.4 && RecoJetArea < 0.7;
    bool isRecoJetNConst = RecoJetNConst >= 1;

    /*
    if (isRecoJetPt && !isRecoJetPtWide) {
      cout << "Case 1 " << RecoJetCorrPt << endl; 
      cout << "DANGER!" << endl; 
      return;
    }
    if (!isRecoJetPt && isRecoJetPtWide) {
      cout << "Case 2" << "\t" << RecoD0Pt << "\t" << RecoJetNConst << "\t"  << RecoJetCorrPt << "\t" <<  RecoJetPtFromArea << endl; 
      cout << "DANGER!" << endl; 
      return;
    }
    */

    // if (!isRecoJetPt) continue;

    // int centbin = CENTWeight->FindBin(Centrality);
    // double wcent = CENTWeight->GetBinContent(centbin);


    double fonllweight = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));
    if (priormode == 2) fonllweight = 1.0;

    if (mcz < mczcutoff) continue;

    bool isOK   = isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta && isRecoJetNConst;
    bool isOKWIDE   = isRecoJetPtWide && isMCJetPt && isRecoJetEta && isMCJetEta && isRecoJetNConst;
    // if (isOK && !isOKWIDE) {cout << "DANGER!" << endl; return;}
    // if (!isOK && isOKWIDE) {cout << "DANGER!" << endl; return;}
    bool isMISS = !isOKWIDE && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOKWIDE && isRecoJetPtWide && isRecoJetEta;
    // bool isFAKE = !isOKWIDE && isRecoJetPtWide;

    double w;
    if (!isForCentralityWeight){
      int grefmultbin = gREFMULTWeight->FindBin(gRefMult);
      double wcent = gREFMULTWeight->GetBinContent(grefmultbin);
      w = wcent*fonllweight*CWeight;
    }
    else w = fonllweight*CWeight;

    // double w = 1.;
    
    if (!isOKWIDE && !isMISS && !isFAKE) continue;

    /// CS
    if (recoz > nbinsz[nz_bins]) {
      recoz = nbinsz[nz_bins] - z_binwidth_high*0.5;
    }
    if (recoz < nbinsz[0]) {
      recoz = nbinsz[0] + z_binwidth_low*0.5;
    }

    if (RecoJetCorrPt < nbinsjetpt[0]) {
      RecoJetCorrPt = nbinsjetpt[0] + jetpt_binwidth_low*0.5;
    }
        
    if (RecoJetCorrPt > nbinsjetpt[njpt_bins]) {
      RecoJetCorrPt = nbinsjetpt[njpt_bins] - jetpt_binwidth_high*0.5;
    }

    /// Area
    if (recozwide > nbinsz_wide[nz_bins_wide]) {
      recozwide = nbinsz_wide[nz_bins_wide] - 0.001;
    }
    if (recozwide < nbinsz_wide[0]) {
      recozwide = nbinsz_wide[0] + 0.001;
    }

    if (RecoJetPtFromArea < nbinsjetpt_wide[0]) {
      RecoJetPtFromArea = nbinsjetpt_wide[0] + 0.001;
    }    
    if (RecoJetPtFromArea > nbinsjetpt_wide[njpt_bins_wide]) {
      RecoJetPtFromArea = nbinsjetpt_wide[njpt_bins_wide] - 0.001;
    }

    bool isRecoZ = recoz > nbinsz[0] && recoz < nbinsz[nz_bins];

    if (priormode == 1){
      int measuredptbin = fMeasWide[centhistogramtofill]->GetXaxis()->FindBin(RecoJetPtFromArea);
      int measuredzbin = fMeasWide[centhistogramtofill]->GetYaxis()->FindBin(recozwide);
      double w1 = DataWeightPt[centhistogramtofill]->GetBinContent(measuredptbin);
      double w2 = DataWeightZ[centhistogramtofill]->GetBinContent(measuredzbin);
      double w3 = DataWeight[centhistogramtofill]->GetBinContent(measuredptbin, measuredzbin);
      w = w*w3;
    }

    if (isOKWIDE){
      fMeasZ1D[centhistogramtofill]->Fill(recoz, w);
      fMeas2D[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
    }

    hCentrality->Fill(Centrality, w);
    hgRefMultCorr->Fill(gRefMult, w);

    hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt, w);
    hDiffJetEta[centhistogramtofill]->Fill(dEta(RecoJetEta, MCJetEta), w);
    hDiffJetPhi[centhistogramtofill]->Fill(dPhi(RecoJetPhi, MCJetPhi), w);

    hDiffJetPtInGenPtBins[centhistogramtofill][genptbin-1]->Fill(RecoJetCorrPt - MCJetPt, w);

    hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt, w);
    hRecoD0Pt[centhistogramtofill]->Fill(RecoD0Pt, w);

    if (isOKWIDE) {
      double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
      double fvaluewide[ndim] = {MCJetPt, mcz, RecoJetPtFromArea, recozwide};

      if (RecoJetNConst == 1){
        singlejetcount[centhistogramtofill]+=w;
      }

      if (MCJetNConst == 1){
        singlemcjetcount[centhistogramtofill]+=w;
      }

      int binJ = hJ[centhistogramtofill]->GetBin(fvalue);

      if (countcutoff > 0) {if (hJ[centhistogramtofill]->GetBinContent(binJ) <= countcutoff) continue;}

      int tmpbin = hJet[centhistogramtofill]->Fill(fvalue, w); //Main body of the response matrix is the only thing I fill in the main THnSparses. The misses are accounted for separately.
      int tmpbinwide = hJetWide[centhistogramtofill]->Fill(fvaluewide, w); //Main body of the response matrix is the only thing I fill in the main THnSparses. The misses are accounted for separately.

      int* coord = new int[ndim];

      Double_t bincount = hJet[centhistogramtofill]->GetBinContent(tmpbin, coord);
			Double_t pttrue = hJet[centhistogramtofill]->GetAxis(0)->GetBinCenter(coord[0]);
			Double_t ztrue = hJet[centhistogramtofill]->GetAxis(1)->GetBinCenter(coord[1]);
			Double_t ptdet =  hJet[centhistogramtofill]->GetAxis(2)->GetBinCenter(coord[2]); 
			Double_t zdet =  hJet[centhistogramtofill]->GetAxis(3)->GetBinCenter(coord[3]);

      if (pttrue < jetpt_var_bin[0] || pttrue > jetpt_var_bin[njpt_gen_bins_var]) cout << fvalue[0] << "\t" << MCD0Pt << "\t" << fvalue[2] << "\t" << RecoD0Pt << endl;
      if (ztrue < z_gen_bin[0] || ztrue > z_gen_bin[nz_gen_bins]) cout << mcz << "\t" << MCD0Pt << "\t" << MCJetPt << "\t" << fvalue[1] << endl;

      // resp[centhistogramtofill]->Fill(ptdet, zdet, pttrue, ztrue, w);
      resp[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
      respwide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Fill(ptdet, pttrue, w);

      fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);

      fMeasWide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, w);

      fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);

      hJetPtAreaVsJetPtCS[centhistogramtofill]->Fill(RecoJetPtFromArea, RecoJetCorrPt, w);

      hJetAreaVsJetNConstArea[centhistogramtofill]->Fill(RecoJetArea, RecoJetNConstFromArea, w);
      hJetAreaVsJetNConstCS[centhistogramtofill]->Fill(RecoJetArea, RecoJetNConst, w);

      hJetPtAreaVsJetNConstArea[centhistogramtofill]->Fill(RecoJetPtFromArea, RecoJetNConstFromArea, w);
      hJetPtCSVsJetNConstCS[centhistogramtofill]->Fill(RecoJetCorrPt, RecoJetNConst, w);

      hCentVsJetNConst->Fill(Centrality, RecoJetNConst, w);

      hGendRvJetPt[centhistogramtofill]->Fill(MCDeltaR, MCJetPt, w);
      hGendR[centhistogramtofill]->Fill(MCDeltaR, w);  

      hdRvJetPt[centhistogramtofill]->Fill(RecoDeltaR, RecoJetPtFromArea, w);
      hdR[centhistogramtofill]->Fill(RecoDeltaR, w);

      okcount[centhistogramtofill] += w;
    }

    else if (isMISS) {

      double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
      double fvaluewide[ndim] = {MCJetPt, mcz, RecoJetPtFromArea, recozwide};

      int binJ = hJ[centhistogramtofill]->GetBin(fvalue);

      if (countcutoff > 0) {if (hJ[centhistogramtofill]->GetBinContent(binJ) <= countcutoff) continue;}

      if (RecoJetNConst == 1){
        singlejetcount[centhistogramtofill]+=w;
      }

      if (MCJetNConst == 1){
        singlemcjetcount[centhistogramtofill]+=w;
      }

      if (!isRecoJetEta) {
        if (abs(RecoJetCorrPt - nbinsjetpt[0]) < abs(RecoJetCorrPt - nbinsjetpt[njpt_bins])){
            fvalue[2] = nbinsjetpt[0] - jetpt_binwidth_low*0.5;
            if (abs(recoz - nbinsz[0]) < abs(recoz - nbinsz[nz_bins])) fvalue[3] = nbinsz[0] - z_binwidth_low*0.5;
            else fvalue[3] = nbinsz[nz_bins] + z_binwidth_high*0.5;
        }
        else{
            fvalue[2] = nbinsjetpt[njpt_bins] + jetpt_binwidth_high*0.5;
            if (abs(recoz - nbinsz[0]) < abs(recoz - nbinsz[nz_bins])) fvalue[3] = nbinsz[0] - z_binwidth_low*0.5;
            else fvalue[3] = nbinsz[nz_bins] + z_binwidth_high*0.5;
        }

        if (abs(RecoJetPtFromArea - nbinsjetpt_wide[0]) < abs(RecoJetPtFromArea - nbinsjetpt_wide[njpt_bins_wide])){
            fvalue[2] = nbinsjetpt_wide[0] - jetpt_binwidth_low_wide*0.5;
            if (abs(recozwide - nbinsz_wide[0]) < abs(recozwide - nbinsz_wide[nz_bins_wide])) fvaluewide[3] = nbinsz_wide[0] - z_binwidth_low_wide*0.5;
            else fvaluewide[3] = nbinsz_wide[nz_bins_wide] + z_binwidth_high_wide*0.5;
        }
        else{
            fvaluewide[2] = nbinsjetpt_wide[njpt_bins_wide] + jetpt_binwidth_high_wide*0.5;
            if (abs(recozwide - nbinsz_wide[0]) < abs(recozwide - nbinsz_wide[nz_bins_wide])) fvaluewide[3] = nbinsz_wide[0] - z_binwidth_low_wide*0.5;
            else fvaluewide[3] = nbinsz_wide[nz_bins_wide] + z_binwidth_high_wide*0.5;
        }
      }

      // hJet[centhistogramtofill]->Fill(fvalue, w); //Misses of the response matrix. Let's see if this works.

      int tmpbin = hJet[centhistogramtofill]->Fill(fvalue, w); //Main body of the response matrix is the only thing I fill in the main THnSparses. The misses are accounted for separately.
      int tmpbinwide = hJetWide[centhistogramtofill]->Fill(fvaluewide, w); //Main body of the response matrix is the only thing I fill in the main THnSparses. The misses are accounted for separately.

      int* coord = new int[ndim];

      Double_t bincount = hJet[centhistogramtofill]->GetBinContent(tmpbin, coord);
      Double_t pttrue = hJet[centhistogramtofill]->GetAxis(0)->GetBinCenter(coord[0]);
      Double_t ztrue = hJet[centhistogramtofill]->GetAxis(1)->GetBinCenter(coord[1]);
      Double_t ptdet =  hJet[centhistogramtofill]->GetAxis(2)->GetBinCenter(coord[2]); 
      Double_t zdet =  hJet[centhistogramtofill]->GetAxis(3)->GetBinCenter(coord[3]);

      if (pttrue < jetpt_var_bin[0] || pttrue > jetpt_var_bin[njpt_gen_bins_var]) cout << fvalue[0] << "\t" << MCD0Pt << "\t" << fvalue[2] << "\t" << RecoD0Pt << endl;
      if (ztrue < z_gen_bin[0] || ztrue > z_gen_bin[nz_gen_bins]) cout << mcz << "\t" << MCD0Pt << "\t" << MCJetPt << "\t" << fvalue[1] << endl;


      fMiss[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fMiss1D[centhistogramtofill]->Fill(MCJetPt, w);

      resp[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      respwide[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Miss(pttrue, w);

      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);

      hGendRvJetPt[centhistogramtofill]->Fill(MCDeltaR, MCJetPt, w);
      hGendR[centhistogramtofill]->Fill(MCDeltaR, w); 

      misscount[centhistogramtofill] += w;
    }

  else if (isFAKE) { // This never happens in the current definition.
    // cout << "FAKE == " << MCD0Pt << "\t" << MCJetPt << "\t" << RecoD0Pt << "\t" << RecoJetCorrPt << endl;

    double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
    double fvaluewide[ndim] = {MCJetPt, mcz, RecoJetPtFromArea, recozwide};

    int binJ = hJ[centhistogramtofill]->GetBin(fvalue);

    if (countcutoff > 0) {if (hJ[centhistogramtofill]->GetBinContent(binJ) <= countcutoff) continue;}

    if (RecoJetNConst == 1){
      singlejetcount[centhistogramtofill]+=w;
    }

    if (MCJetNConst == 1){
      singlemcjetcount[centhistogramtofill]+=w;
    }

      // hJet[centhistogramtofill]->Fill(fvalue, w); //Misses of the response matrix. Let's see if this works.

    int tmpbin = hJet[centhistogramtofill]->Fill(fvalue, w); //Main body of the response matrix is the only thing I fill in the main THnSparses. The misses are accounted for separately.
    int tmpbinwide = hJetWide[centhistogramtofill]->Fill(fvaluewide, w); //Main body of the response matrix is the only thing I fill in the main THnSparses. The misses are accounted for separately.


    resp[centhistogramtofill]->Fake(RecoJetCorrPt, recoz, w);
    respwide[centhistogramtofill]->Fake(RecoJetPtFromArea, recozwide, w);
    resp1D[centhistogramtofill]->Fake(RecoJetCorrPt, w);

    fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
    fMeasWide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, w);
    fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);

    fFake[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
    fFakeWide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, w);
    fFake1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);

    fakecount[centhistogramtofill] += w;
  }
}

cout << endl;

for (int i = 0; i < 3; i++){
  cout << "Centrality = " << i << "\t" << okcount[i] + misscount[i] << "\t" << okcount[i] << "\t" << misscount[i] << "\t" << fakecount[i] << endl;
}

for (int i = 0; i < 3; i++){
  cout << "Centrality = " << i << "\t" << fTrue[i]->Integral() << "\t" << fMeas[i]->Integral() << "\t" << fMiss[i]->Integral() << endl;
}

for (int i = 0; i < 3; i++){
  cout << "Wide pT Central = " << i << "\t" << fTrue[i]->Integral() << "\t" << fMeasWide[i]->Integral() << "\t" << fMiss[i]->Integral() << endl;
}

for (int i = 0; i < 3; i++){
  cout << "MC = " << i << "\t" << singlemcjetcount[i] << "\t" << fTrue[i]->Integral() << "\t" << singlemcjetcount[i]/fTrue[i]->Integral() << endl;
  cout << "Reco = " << i << "\t" << singlejetcount[i] << "\t" << fMeas[i]->Integral() << "\t" << singlejetcount[i]/fMeas[i]->Integral() << endl;
}


// cout << "Unaccepted = " << recooutside[0] << "\t" << recooutside[1] << "\t" << recooutside[2] << endl;

TString filename;
switch (mode){
    case 1:
    filename = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins, mode);
    break;
    case 2:
    filename = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins, mode);
    break;
    default:
    filename = Form("%s/Response_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 0, 0, njpt_gen_bins_var, nz_gen_bins);
}

  cout << "Centrality Integral = " << hCentrality->Integral() << endl;

  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->cd();

  for (int cent = 0; cent < 3; cent++){
      for (int ptbin = 0; ptbin < njpt_gen_bins_var; ptbin++){
        hDiffJetPtInGenPtBins[cent][ptbin]->Write();
    }
  }

  for (int cent = 0;  cent < 3; cent++){
    hZ[cent]->Write();
    hJetPtvZ[cent]->Write();
    
    fMeas[cent]->Write();
    fTrue[cent]->Write();
    fMiss[cent]->Write();
    fFake[cent]->Write();

    fMeasWide[cent]->Write();
    fFakeWide[cent]->Write();

    fMeas1D[cent]->Write();
    fTrue1D[cent]->Write();
    fMiss1D[cent]->Write();
    fMissZ1D[cent]->Write();
    fFake1D[cent]->Write();

    fMeas2D[cent]->Write();
    fMeasZ1D[cent]->Write();
    fTrueZ1D[cent]->Write();
    fMissZ1D[cent]->Write();

    gDirectory->WriteObject(resp[cent], Form("Resp_%i", cent));
    gDirectory->WriteObject(respwide[cent], Form("RespWide_%i", cent));
    gDirectory->WriteObject(resp1D[cent], Form("Resp1D_%i", cent));

    hDiffJetPt[cent]->Write();
    hDiffJetEta[cent]->Write();
    hDiffJetPhi[cent]->Write();

    hMCD0Pt[cent]->Write();
    hRecoD0Pt[cent]->Write();

    hJetPtAreaVsJetPtCS[cent]->Write();
    hJetAreaVsJetNConstArea[cent]->Write();
    hJetAreaVsJetNConstCS[cent]->Write();
    hJetPtAreaVsJetNConstArea[cent]->Write();
    hJetPtCSVsJetNConstCS[cent]->Write();

    hGendRvJetPt[cent]->Write();
    hGendR[cent]->Write();

    hdRvJetPt[cent]->Write();
    hdR[cent]->Write();
    // Stuff for SuperIteration Starts Here

    hJ[cent] = (THnSparseF*)hJet[cent]->Clone(Form("hJ_%i", cent));

    hJ[cent]->Write();
    hJet[cent]->Write();
    hJetWide[cent]->Write();
  }

  hCentVsJetNConst->Write();

  hCentrality->Write();
  hgRefMultCorr->Write();
  hRecoD0Mass->Write();

  outfile->Close();

}

void PrepareResponseHistograms(TString DirName = "MCMCUnf", int D0pTLow = 1, int mode = 0, bool isForCentralityWeight = false, int priormode = 0){
  // gSystem->ListLibraries();

  Method(DirName, D0pTLow, mode, isForCentralityWeight, priormode);
}