R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
// #include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int D0pTLow = 1, int mode = 0, bool isForCentralityWeight = false, int priormode = 1){
  // PriorMode
  // 0 -> FONLL (default)
  // 1 -> FONLL + Reweighed with data distros
  // 2 -> PYTHIA
  // 3 -> PYTHIA With Fit
  // 4 -> FONLL (default) With Tracking Efficiency
  // 5 -> PYTHIA With Weights Derived From DATA

  cout << D0pTLow << "\t" << mode  << "\t" << isForCentralityWeight << "\t" << priormode << endl; 

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH1::StatOverflows();

const int njpt_gen_bins_var = 10; //x2 MC
double jetpt_var_bin[njpt_gen_bins_var+1] = {1, 3, 5,7,9,11,13,15,20,25,30};
double jetpt_gen_binwidth_low = 2.;
double jetpt_gen_binwidth_high = 5.;

const int nz_gen_bins = 7; //x2 MC
double z_gen_bin[nz_gen_bins+1] = {0., 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.};
double z_gen_binwidth_low = 0.1;
double z_gen_binwidth_high = 0.1;

const int njpt_bins_wide = 16;
double nbinsjetpt_wide[njpt_bins_wide + 1] = {-10, -5, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 40, 50};

double jetpt_binwidth_low_wide = 5.;
double jetpt_binwidth_high_wide = 10.;

const int nz_bins_wide = 24;
double nbinsz_wide[nz_bins_wide + 1] = {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000};

double z_binwidth_low_wide = 500.;
double z_binwidth_high_wide = 500.;

  
  // TFile *f = new TFile("HIOverlay_CS_May20.root");
  // TFile *f = new TFile("HIOverlay_HFJets_WithCS_Jun6_2023_pthat_3_inf.root");
  TFile *f;
  if (priormode != 4)f = new TFile("HIOverlay_HFJets_WithCS_Jun20_2023_pthat_3_inf.root");
  else f = new TFile("HIOverlay_HFJets_WithCS_TrackEff_Aug17_2023_pthat_3_inf.root");
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
  TH1F  *FITvPYTHIAWeights = (TH1F *) PYTHIA->Get("FitWeights");

  TH1F  *DataWeighedFit[3];
  for (int i = 0; i < 3; i++){
    DataWeighedFit[i] = (TH1F *) PYTHIA->Get(Form("DataWeighed_PtWeight_%i", i));
  }

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

  TH2D *Weight[3];
  RooUnfoldResponse *respwide[3]; //Response
  RooUnfoldResponse *resp1D[3]; //Response
  
  TH2D *fTrue[3];
  TH2D *fMiss[3];
  TH2D *fMeasWide[3];
  TH2D *fFakeWide[3];

  TH1D *fMeas1D[3];
  TH1D *fTrue1D[3];
  TH1D *fMiss1D[3];
  TH1D *fFake1D[3];

  // THnSparseF *masterTHn = new THnSparseF("masterTHn", "masterTHn", ndim, )

  cout << JetTree->GetEntries() << endl;

  for (int i = 0; i < 3; i++){


    fTrue[i] = new TH2D(Form("fTrue_Cent_%i", i), Form("fTrue_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
    fMiss[i] = new TH2D(Form("fMiss_Cent_%i", i), Form("fMiss_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin, nz_gen_bins, z_gen_bin);
    fMeasWide[i] = new TH2D(Form("fMeasWide_Cent_%i", i), Form("fMeasWide_Cent_%i", i), njpt_bins_wide, nbinsjetpt_wide, nz_bins_wide, nbinsz_wide);
    fFakeWide[i] = new TH2D(Form("fFakeWide_Cent_%i", i), Form("fFakeWide_Cent_%i", i), njpt_bins_wide, nbinsjetpt_wide, nz_bins_wide, nbinsz_wide);

    fMeas1D[i] = new TH1D(Form("fMeas1D_Cent_%i", i), Form("fMeas1D_Cent_%i", i), njpt_bins_wide, nbinsjetpt_wide);
    fTrue1D[i] = new TH1D(Form("fTrue1D_Cent_%i", i), Form("fTrue1D_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin);
    fMiss1D[i] = new TH1D(Form("fMiss1D_Cent_%i", i), Form("fMiss1D_Cent_%i", i), njpt_gen_bins_var, jetpt_var_bin);
    fFake1D[i] = new TH1D(Form("fFake1D_Cent_%i", i), Form("fFake1D_Cent_%i", i), njpt_bins_wide, nbinsjetpt_wide);

    // cout << "Here" << endl;
    respwide[i] = new RooUnfoldResponse(Form("RespWide_%i", i), Form("RespWide_%i", i));
    respwide[i]->Setup(fMeasWide[i], fTrue[i]); //Setup Response Matrix Definition

    resp1D[i] = new RooUnfoldResponse(Form("Resp1D_%i", i), Form("Resp1D_%i", i));
    resp1D[i]->Setup(fMeas1D[i], fTrue1D[i]); //Setup Response Matrix Definition

  }

cout << JetTree->GetEntries() << endl;

int nentries = JetTree->GetEntries();

int lowlimit = 0;
int highlimit = 1000000;

int counter = 0;

TH1D *plottingpT = new TH1D("plottingpT", "plottingpT", njpt_gen_bins_var, jetpt_var_bin);

cout << "Mode == " << mode << endl;

cout << "Reading events from " << lowlimit << " to " << highlimit << endl;

cout << "================ Limits =================" << endl;
cout << Form("%i < pT,D0 [GeV/c] < %i", D0pTLow, 10) << endl;
cout << Form("%.1f < MC pT,Jet [GeV/c] < %.1f", (float)D0pTLow, jetpt_var_bin[njpt_gen_bins_var]) << endl;
cout << Form("%.1f < Reco pT,Jet [GeV/c] < %.1f", (float)nbinsjetpt_wide[0], nbinsjetpt_wide[njpt_bins_wide]) << endl;
cout << Form("%.1f < MC Z,Jet < %.1f", z_gen_bin[0], z_gen_bin[nz_gen_bins]) << endl;
cout << Form("%.1f < Reco Z,Jet < %.1f", nbinsz_wide[0], nbinsz_wide[nz_bins_wide]) << endl;
cout << "=========================================" << endl;


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

    bool isMCD0Pt = MCD0Pt > D0pTLow && MCD0Pt < 10; // We wanna check 1 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > D0pTLow && RecoD0Pt < 10;
    bool isMCJetPt = MCJetPt > D0pTLow && MCJetPt < jetpt_var_bin[njpt_gen_bins_var];

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (!isMCJetPt) continue;
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
    if (!isMCZ) {
        cout << "MCZ Not Satistied" << endl;
        return;
        // continue;
    } 

    int genptbin = plottingpT->FindBin(MCJetPt);

    if(RecoJetPhi>=2.*TMath::Pi()) RecoJetPhi = RecoJetPhi - 2.*TMath::Pi();
    if(RecoJetPhi<0) RecoJetPhi = 2.*TMath::Pi() + RecoJetPhi;

    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    double RecoJetPxWide = RecoJetPtFromArea *TMath::Cos(RecoJetPhi);
    double RecoJetPyWide = RecoJetPtFromArea*TMath::Sin(RecoJetPhi);
    double recozwide = (RecoJetPxWide*RecoD0Px + RecoJetPyWide*RecoD0Py)/pow(RecoJetPtFromArea, 2);

    TLorentzVector RecoPion;
    TLorentzVector RecoKaon;
    TLorentzVector RecoD0;

    RecoPion.SetPtEtaPhiM(RecoPionPt, RecoPionEta, RecoPionPhi, 0.13957);
    RecoKaon.SetPtEtaPhiM(RecoKaonPt, RecoKaonEta, RecoKaonPhi, 0.49368);
    RecoD0 = RecoPion + RecoKaon;


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
  
    //Defining response matrix using fakes and misses to check closure
    bool isRecoJetPtWide = RecoJetPtFromArea >= nbinsjetpt_wide[0] && RecoJetPtFromArea <= nbinsjetpt_wide[njpt_bins_wide];
    bool isMCJetEta = abs(MCJetEta) < 0.601;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoArea = RecoJetArea > 0.4 && RecoJetArea < 0.7;
    bool isRecoJetNConst = RecoJetNConst >= 1;

    double priorweight = 1.;
    double fonllweight = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));
    double fitweight = FONLLvPYTHIAWeights->GetBinContent(FITvPYTHIAWeights->FindBin(MCJetPt));
    double dataweighedweight = DataWeighedFit[centhistogramtofill]->GetBinContent(DataWeighedFit[centhistogramtofill]->FindBin(MCJetPt));
    if (priormode == 0 || priormode == 1 || priormode == 4) priorweight = fonllweight;
    if (priormode == 2) priorweight = 1.0;
    if (priormode == 3) priorweight = fitweight;
    if (priormode == 5) priorweight = dataweighedweight;

    bool isOKWIDE   = isRecoJetPtWide && isMCJetPt && isRecoJetEta && isMCJetEta;
    bool isMISS = !isOKWIDE && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOKWIDE && isRecoJetPtWide && isRecoJetEta;

    double w;
    w = priorweight*CWeight;

    if (!isOKWIDE && !isMISS && !isFAKE) continue;

    if (priormode == 1){
      int measuredptbin = fMeasWide[centhistogramtofill]->GetXaxis()->FindBin(RecoJetPtFromArea);
      int measuredzbin = fMeasWide[centhistogramtofill]->GetYaxis()->FindBin(recozwide);
      double w1 = DataWeightPt[centhistogramtofill]->GetBinContent(measuredptbin);
      double w2 = DataWeightZ[centhistogramtofill]->GetBinContent(measuredzbin);
      double w3 = DataWeight[centhistogramtofill]->GetBinContent(measuredptbin, measuredzbin);
      w = w*w3;
    }

    if (isOKWIDE) {
      respwide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Fill(RecoJetPtFromArea, MCJetPt, w);

      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);

      fMeasWide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, w);
      fMeas1D[centhistogramtofill]->Fill(RecoJetPtFromArea, w);
    }

    else if (isMISS) {
      respwide[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Miss(MCJetPt, w);
      
      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);

      fMiss[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fMiss1D[centhistogramtofill]->Fill(MCJetPt, w);

      misscount[centhistogramtofill] += w;
    }

    else if (isFAKE) { // This never happens in the current definition.
        // cout << "FAKE == " << MCD0Pt << "\t" << MCJetPt << "\t" << RecoD0Pt << "\t" << RecoJetCorrPt << endl;
        respwide[centhistogramtofill]->Fake(RecoJetPtFromArea, recozwide, w);
        resp1D[centhistogramtofill]->Fake(RecoJetPtFromArea, w);

        fMeasWide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, w);
        fMeas1D[centhistogramtofill]->Fill(RecoJetPtFromArea, w);

        fFakeWide[centhistogramtofill]->Fill(RecoJetPtFromArea, recozwide, w);
        fFake1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);

        fakecount[centhistogramtofill] += w;
    }
}

cout << endl;

TFile *DataFile = new TFile(Form("../Data/Aug17_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_%i.root", D0pTLow, D0pTLow, 1000), "READ");
DataFile->cd();

TH2D *MeasuredWide[3];

MeasuredWide[0] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_0_10");
MeasuredWide[1] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_10_40");
MeasuredWide[2] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_40_80");

TH1D *Measured1D[3];

Measured1D[0] = (TH1D *)MeasuredWide[0]->ProjectionX();
Measured1D[1] = (TH1D *)MeasuredWide[1]->ProjectionX();
Measured1D[2] = (TH1D *)MeasuredWide[2]->ProjectionX();

RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
RooUnfold::ErrorTreatment errorTreatment1D = RooUnfold::kErrors;

int iteration = 4;

RooUnfoldBayes unfold1Dcent (resp1D[0], Measured1D[0], iteration);
RooUnfoldBayes unfold1Dmid (resp1D[1], Measured1D[1], iteration);
RooUnfoldBayes unfold1Dperi (resp1D[2], Measured1D[2], iteration);

RooUnfoldBayes unfoldwidecent (respwide[0], MeasuredWide[0], iteration+2);
RooUnfoldBayes unfoldwidemid (respwide[1], MeasuredWide[1], iteration+2);
RooUnfoldBayes unfoldwideperi (respwide[2], MeasuredWide[2], iteration+2);

TH1D *Unfolded1D[3];

Unfolded1D[0] = (TH1D *)unfold1Dcent.Hreco(errorTreatment1D);
Unfolded1D[1] = (TH1D *)unfold1Dmid.Hreco(errorTreatment1D);
Unfolded1D[2] = (TH1D *)unfold1Dperi.Hreco(errorTreatment1D);

TH2D *UnfoldedWide[3];

UnfoldedWide[0] = (TH2D *)unfoldwidecent.Hreco(errorTreatment);
UnfoldedWide[1] = (TH2D *)unfoldwidemid.Hreco(errorTreatment);
UnfoldedWide[2] = (TH2D *)unfoldwideperi.Hreco(errorTreatment);

cout << "Measured Integral" << endl;
for (int i = 0; i < 3; i++){
    cout << i << "\t" << Measured1D[i]->Integral(0,  Measured1D[i]->GetNbinsX()+1) << endl;
    cout << i << "\t" << MeasuredWide[i]->Integral() << endl;
}

cout << "Fake Integral" << endl;
for (int i = 0; i < 3; i++){
    cout << i << "\t" << fFake1D[i]->Integral(0,  fFake1D[i]->GetNbinsX()+1) << endl;
    cout << i << "\t" << fFakeWide[i]->Integral() << endl;
}

cout << "Miss Integral" << endl;
for (int i = 0; i < 3; i++){
    cout << i << "\t" << fMiss1D[i]->Integral(0,  fMiss1D[i]->GetNbinsX()+1) << endl;
    cout << i << "\t" << fMiss[i]->Integral() << endl;
}

cout << "True Integral" << endl;
for (int i = 0; i < 3; i++){
    cout << i << "\t" << fTrue1D[i]->Integral(0,  fTrue1D[i]->GetNbinsX()+1) << endl;
    cout << i << "\t" << fTrue[i]->Integral() << endl;

    cout << "% of misses = " << fMiss1D[i]->Integral(0,  fMiss1D[i]->GetNbinsX()+1)/fTrue1D[i]->Integral(0,  fTrue1D[i]->GetNbinsX()+1) << endl;
}

cout << "Unfolded Integral" << endl;
for (int i = 0; i < 3; i++){
    cout << i << "\t" << Unfolded1D[i]->Integral(0, Unfolded1D[i]->GetNbinsX()+1) << endl;
    cout << i << "\t" << UnfoldedWide[i]->Integral() << endl;
}

TH1D *UnfoldedPt[3];
TH1D *UnfoldedZ[3];

UnfoldedPt[0] = (TH1D *)UnfoldedWide[0]->ProjectionX();
UnfoldedPt[1] = (TH1D *)UnfoldedWide[1]->ProjectionX();
UnfoldedPt[2] = (TH1D *)UnfoldedWide[2]->ProjectionX();

UnfoldedZ[0] = (TH1D *)UnfoldedWide[0]->ProjectionY();
UnfoldedZ[1] = (TH1D *)UnfoldedWide[1]->ProjectionY();
UnfoldedZ[2] = (TH1D *)UnfoldedWide[2]->ProjectionY();

TCanvas *c = new TCanvas("c", "c", 1200, 600);
c->Divide(3);

for (int i = 0; i < 3; i++){
    c->cd(i+1);
    gPad->SetLogy();
    SetColor(Unfolded1D[i], kRed, 20);
    SetColor(UnfoldedPt[i], kCyan, 21);
    Unfolded1D[i]->GetXaxis()->SetRangeUser(5,20);
    Unfolded1D[i]->Draw("EP");
    UnfoldedPt[i]->Draw("EP SAME");
}

int lowptbin = Unfolded1D[0]->FindBin(5.0);
int highptbin = Unfolded1D[0]->FindBin(20.0);

TH1D *VariationPt[3];

for (int i = 0; i < 3; i++){
    VariationPt[i] = (TH1D *)UnfoldedPt[i]->Clone();
    VariationPt[i]->Reset();
    for (int j = lowptbin; j <= highptbin; j++){
        double relerr = (Unfolded1D[i]->GetBinContent(j) - UnfoldedPt[i]->GetBinContent(j))*100./UnfoldedPt[i]->GetBinContent(j);
        VariationPt[i]->SetBinContent(j, relerr);
        VariationPt[i]->SetBinError(j, 0);
    }
}

TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
c2->Divide(3);
for (int i = 0; i < 3; i++){
    c2->cd(i+1);
    // gPad->SetLogy();
    VariationPt[i]->GetXaxis()->SetRangeUser(5, 20);
    VariationPt[i]->GetYaxis()->SetRangeUser(-100, 100);
    SetColor(VariationPt[i], kRed, 20);
    VariationPt[i]->Draw("EP");
}


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

  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->cd();

  for (int cent = 0;  cent < 3; cent++){

    fTrue[cent]->Write();
    fTrue[cent]->ProjectionX()->Write(Form("TruePt_%i", cent));
    fTrue[cent]->ProjectionY()->Write(Form("TrueZ_%i", cent));
    fMiss[cent]->Write();
    fMiss[cent]->ProjectionX()->Write(Form("MissPt_%i", cent));
    fMiss[cent]->ProjectionY()->Write(Form("MissZ_%i", cent));
    fMeasWide[cent]->Write();
    fMeasWide[cent]->ProjectionX()->Write(Form("MeasPt_%i", cent));
    fMeasWide[cent]->ProjectionY()->Write(Form("MeasZ_%i", cent));
    fFakeWide[cent]->Write();
    fFakeWide[cent]->ProjectionX()->Write(Form("FakePt_%i", cent));
    fFakeWide[cent]->ProjectionY()->Write(Form("FakeZ_%i", cent));

    fMeas1D[cent]->Write();
    fTrue1D[cent]->Write();
    fMiss1D[cent]->Write();
    fFake1D[cent]->Write();

    gDirectory->WriteObject(respwide[cent], Form("RespWide_%i", cent));
    gDirectory->WriteObject(resp1D[cent], Form("Resp1D_%i", cent));
  }

  outfile->Close();

}

void Test(TString DirName = "tmp", int D0pTLow = 1, int mode = 0, bool isForCentralityWeight = false, int priormode = 0){
  // gSystem->ListLibraries();
  
  gSystem->mkdir(DirName.Data(), true);
  Method(DirName, D0pTLow, mode, isForCentralityWeight, priormode);
}