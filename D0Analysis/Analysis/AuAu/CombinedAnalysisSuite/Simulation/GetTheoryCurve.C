R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void Method(TString DirName = "MCMCUnf", int D0pTLow = 1, int mode = 0, bool isForCentralityWeight = false, int priormode = 1){
  // PriorMode
  // 0 -> FONLL (default)
  // 1 -> FONLL + Reweighed with data distros
  // 2 -> PYTHIA
  // 3 -> PYTHIA With Fit
  // 4 -> FONLL (default) With Tracking Efficiency
  // 5 -> PYTHIA With Weights Derived From DATA
  // 6 -> FONLL (default) D0 Yield Varied In Data In Positive Direction
  // 7 -> FONLL (default) D0 Yield Varied In Data In Negative Direction
  // 8 -> FONLL (default) D0 Reconstruction Efficiency Without Vertex Correction Varied In Data In Positive Direction
  // 9 -> FONLL (default) D0 Reconstruction Efficiency Without Vertex Correction Varied In Data In Negative Direction
  // 10 -> FONLL (default) D0 Reconstruction Efficiency With Vertex Correction Varied In Data In Positive Direction
  // 11 -> FONLL (default) D0 Reconstruction Efficiency With Vertex Correction Varied In Data In Negative Direction
  // 12 -> LIDO (This is underdeveloped. Only available for 1-10 GeV/c for central and peripheral)

  cout << D0pTLow << "\t" << mode  << "\t" << isForCentralityWeight << "\t" << priormode << endl; 

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
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

  TH1D *LIDOPtCurveWeight[3];

  TFile *LIDOWeightFile = new TFile("LIDOHistogram_1_20.root");
  for (int i = 0; i < 3; i++){
    LIDOPtCurveWeight[i] = (TH1D *) LIDOWeightFile->Get(Form("Weight %i", i));
    cout << LIDOPtCurveWeight[i]->GetName() << "\t" << LIDOPtCurveWeight[i]->GetNbinsX() << endl;
    LIDOPtCurveWeight[i]->SetDirectory(0);
  }

  
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
  RooUnfoldResponse *respdR[3]; //Response
  // RooUnfoldResponse *resp1DdR[3]; //Response

  TH2D *fMeas[3];
  TH2D *fTrue[3];
  TH2D *fMiss[3];
  TH2D *fFake[3];

  TH2D *fMeasWide[3];
  TH2D *fFakeWide[3];

  TH2D *fMeasdRWide[3];
  TH2D *fTruedR[3];
  TH2D *fMissdR[3];
  TH2D *fFakedRWide[3];

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

  TH1D *hMCJetPt[3];
  TH1D *hMCJetZ[3];
  TH1D *hMCJetdR[3];

  // THnSparseF *masterTHn = new THnSparseF("masterTHn", "masterTHn", ndim, )

  cout << JetTree->GetEntries() << endl;

  for (int i = 0; i < 3; i++){
    hMCJetPt[i] = new TH1D(Form("hMCJetPt_%i", i), Form("hMCJetPt_%i", i), njpt_gen_bins_var, jetpt_var_bin);
    hMCJetZ[i] = new TH1D(Form("hMCJetZ_%i", i), Form("hMCJetZ_%i", i), nz_gen_bins, z_gen_bin);
    hMCJetdR[i] = new TH1D(Form("hMCJetdR_%i", i), Form("hMCJetdR_%i", i), ndrbins, drbins);
  }
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
cout << Form("%.1f < MC pT,Jet [GeV/c] < %.1f", (float)D0pTLow, jetpt_var_bin[njpt_gen_bins_var]) << endl;
cout << Form("%.1f < Reco pT,Jet [GeV/c] < %.1f", (float)D0pTLow, nbinsjetpt[njpt_bins]) << endl;
cout << Form("%.1f < MC Z,Jet < %.1f", z_gen_bin[0], z_gen_bin[nz_gen_bins]) << endl;
cout << Form("%.1f < Reco Z,Jet < %.1f", nbinsz[0], nbinsz[nz_bins]) << endl;
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
    bool isMCJetPt = MCJetPt > 5. && MCJetPt < 20.;

    if (!isMCD0Pt) continue;
    if (!isMCJetPt) continue;

    double fitweight = FONLLvPYTHIAWeights->GetBinContent(FITvPYTHIAWeights->FindBin(MCJetPt));
    double priorweight = fitweight;

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);
    double MCDeltaR = dR(dEta(MCJetEta, MCD0Eta), dPhi(MCJetPhi, MCD0Phi));
    // double mcz = MCD0Pt/MCJetPt;

    // if (mcz > 1.0) cout << MCJetPt << "\t" << MCD0Pt << mcz << endl;
    if (mcz >= 1.0) mcz = 0.999; // Padding the boundaries

    // bool isMCZ   = mcz > z_gen_bin[0] && mcz < z_gen_bin[nz_gen_bins];
    // if (!isMCZ) continue; 

    hMCJetPt[centhistogramtofill]->Fill(MCJetPt, priorweight);
    hMCJetZ[centhistogramtofill]->Fill(mcz, priorweight);
    hMCJetdR[centhistogramtofill]->Fill(MCDeltaR, priorweight);
    
}

cout << endl;



// cout << "Unaccepted = " << recooutside[0] << "\t" << recooutside[1] << "\t" << recooutside[2] << endl;

TString filename;
filename = "PYTHIA_TheoryCurve_Jets.root";

  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->cd();


  for (int cent = 0;  cent < 3; cent++){
    hMCJetPt[cent]->Write();
    hMCJetZ[cent]->Write();
    hMCJetdR[cent]->Write();
  }

  outfile->Close();

}

void GetTheoryCurve(TString DirName = "MCMCUnf", int D0pTLow = 1, int mode = 0, bool isForCentralityWeight = false, int priormode = 0){
  // gSystem->ListLibraries();

  Method(DirName, D0pTLow, mode, isForCentralityWeight, priormode);
}