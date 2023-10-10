R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

TH1D *hZ[3];
TH2D *hJetPtvZ[3];
THnSparseF *hJet[3];
THnSparseF *hJetZNorm[3];
THnSparseF *hJetZNormPtNorm[3];
TH2D *Weight[3];
RooUnfoldResponse *resp[3]; //Response
TH2D *fMeas[3];
TH2D *fTrue[3];
TH2D *fMiss[3];
TH2D *fFake[3];

TH2D *hMeasured[3];
TH2D *hTrue[3];
TH1D *hMeasured1D[3];
TH1D *hTrue1D[3];

RooUnfoldResponse *resp1D[3]; //Response
TH2D *fResp1D[3];
TH1D *fMeas1D[3];
TH1D *fTrue1D[3];
TH1D *fMiss1D[3];
TH1D *fFake1D[3];

TH1D *fMissRatio[3];
TH1D *fFakeRatio[3];

TH2D *hMiss1D[3]; //These misses are because of eta-phi mismatch and not a pT mismatch

TH1D *hDiffJetPt[3]; //Reco - MC
TH1D *hDiffJetPtInGenPtBins[3][nbins_jpt]; //Reco - MC
TH1D *hDiffJetEta[3]; //Reco - MC
TH1D *hDiffJetEtaBeforeSmear[3]; //Reco - MC
TH1D *hDiffJetEtaAfterSmear[3]; //Reco - MC
TH1D *hDiffJetPhi[3]; //Reco - MC

TH1D *hMCD0Pt[3]; //Reco - MC
TH1D *hRecoD0Pt[3]; //Reco - MC

TH1D *Smear[3];

TH1F *dPt[3][40]; // This is the dPt vs Pt generated from Embed_April2.root
TH1F *Eta_Mean[3];
TH1F *Eta_Sigma[3];
TH1F *Phi_Mean[3];
TH1F *Phi_Sigma[3];


void Method(int cent = 0, int iteration = 3, TString DirName = "MCMCUnf", int mode = 0){

  cout << "File Reader Method" << endl;

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TFile *f = new TFile("/Volumes/WorkDrive/April10_MCResponse_NoJetPtCutOff/MCResponse_April10_NoJetpTCutOff.root");
  cout << f->GetName() <<  endl;
  if (cent == 0) f->cd("SimJetSaverSmear_Central");
  else if (cent == 1) f->cd("SimJetSaverSmear_MidCentral");
  else if (cent == 2) f->cd("SimJetSaverSmear_Peripheral");
  
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
  // JetTree->SetBranchAddress("RecoPrimaryVertex", &RecoPrimaryVertex);
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
  JetTree->SetBranchAddress("RecoJetEta", &RecoJetEta);
  JetTree->SetBranchAddress("RecoJetPhi", &RecoJetPhi);
  JetTree->SetBranchAddress("RecoJetArea", &RecoJetArea);
  JetTree->SetBranchAddress("RecoJetE", &RecoJetE);
  JetTree->SetBranchAddress("RecoJetNConst", &RecoJetNConst);
  JetTree->SetBranchAddress("RecoD0Pt", &RecoD0Pt);
  JetTree->SetBranchAddress("RecoD0Eta", &RecoD0Eta);
  JetTree->SetBranchAddress("RecoD0Phi", &RecoD0Phi);
  JetTree->SetBranchAddress("RecoPionPt", &RecoPionPt);
  JetTree->SetBranchAddress("RecoPionEta", &RecoPionEta);
  JetTree->SetBranchAddress("RecoPionPhi", &RecoPionPhi);
  JetTree->SetBranchAddress("RecoKaonPt", &RecoKaonPt);
  JetTree->SetBranchAddress("RecoKaonEta", &RecoKaonEta);
  JetTree->SetBranchAddress("RecoKaonPhi", &RecoKaonPhi);

  TFile *SmearFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/BKGHistogramsFINALDec14_22.root");
  SmearFile->cd();

  Smear[cent] = (TH1D *)gDirectory->Get("hBKG_0_10_Weighted");
  Smear[cent] = (TH1D *)gDirectory->Get("hBKG_10_40_Weighted");
  Smear[cent] = (TH1D *)gDirectory->Get("hBKG_40_80_Weighted");

  TFile *pTdepdpTSmearing = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/Pt_Res_FINAL_Jan26_2023.root");
  TFile *EtaResFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/Eta_Res_FINAL_Jan26_2023.root");
  TFile *PhiResFile = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Background/Phi_Res_FINAL_Jan26_2023.root");

  for (int i = 0; i < 40; i++){
    dPt[0][i] = (TH1F *)pTdepdpTSmearing->Get(Form("dPt_Cent_%i_%i_Pt_%i_%i", 0, 10, i, i+1));

    // cout << dPt[0][i]->Integral() << endl;
  }
  Eta_Mean[0] = (TH1F *)EtaResFile->Get(Form("eta_mean_pt_%i_%i", 0, 10));
  Eta_Sigma[0] = (TH1F *)EtaResFile->Get(Form("eta_sigma_pt_%i_%i", 0, 10));

  Phi_Mean[0] = (TH1F *)PhiResFile->Get(Form("phi_mean_pt_%i_%i", 0, 10));
  Phi_Sigma[0] = (TH1F *)PhiResFile->Get(Form("phi_sigma_pt_%i_%i", 0, 10));

  // cout << Eta_Sigma[0]->Integral() << endl;

  for (int i = 0; i < 40; i++){
    dPt[1][i] = (TH1F *)pTdepdpTSmearing->Get(Form("dPt_Cent_%i_%i_Pt_%i_%i", 10, 40, i, i+1));
  }
  Eta_Mean[1] = (TH1F *)EtaResFile->Get(Form("eta_mean_pt_%i_%i", 10, 40));
  Eta_Sigma[1] = (TH1F *)EtaResFile->Get(Form("eta_sigma_pt_%i_%i", 10, 40));

  Phi_Mean[1] = (TH1F *)PhiResFile->Get(Form("phi_mean_pt_%i_%i", 10, 40));
  Phi_Sigma[1] = (TH1F *)PhiResFile->Get(Form("phi_sigma_pt_%i_%i", 10, 40));

  // cout << Eta_Sigma[1]->Integral() << endl;

  for (int i = 0; i < 40; i++){
    dPt[2][i] = (TH1F *)pTdepdpTSmearing->Get(Form("dPt_Cent_%i_%i_Pt_%i_%i", 40, 80, i, i+1));
  }

  Eta_Mean[2] = (TH1F *)EtaResFile->Get(Form("eta_mean_pt_%i_%i", 40, 80));
  Eta_Sigma[2] = (TH1F *)EtaResFile->Get(Form("eta_sigma_pt_%i_%i", 40, 80));

  Phi_Mean[2] = (TH1F *)PhiResFile->Get(Form("phi_mean_pt_%i_%i", 40, 80));
  Phi_Sigma[2] = (TH1F *)PhiResFile->Get(Form("phi_sigma_pt_%i_%i", 40, 80));

  // cout << Eta_Sigma[2]->Integral() << endl;


  cout << JetTree->GetEntries() << endl;

  hZ[cent] = new TH1D(Form("hZ_%i", cent), Form("hZ_%i", cent), nz_gen_bins, z_gen_low, z_gen_high);
  hJetPtvZ[cent] = new TH2D(Form("hJetPtvZ_%i", cent), Form("hJetPtvZ_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
  hJet[cent] = new THnSparseF(Form("hJet_%i", cent), Form("hJet_%i", cent), ndim, nbins, xmin, xmax);
  hJetZNorm[cent] = new THnSparseF(Form("hJetZNorm_%i", cent), Form("hJetZNorm_%i", cent), ndim, nbins, xmin, xmax);
  hJetZNormPtNorm[cent] = new THnSparseF(Form("hJetZNormPtNorm_%i", cent), Form("hJetZNormPtNorm_%i", cent), ndim, nbins, xmin, xmax);

  fMeas[cent] = new TH2D(Form("fMeas_Cent_%i", cent), Form("fMeas_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);
  fTrue[cent] = new TH2D(Form("fTrue_Cent_%i", cent), Form("fTrue_Cent_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
  fMiss[cent] = new TH2D(Form("fMiss_Cent_%i", cent), Form("fMiss_Cent_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);
  fFake[cent] = new TH2D(Form("fFake_Cent_%i", cent), Form("fFake_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high, nz_bins, z_low, z_high);

  resp[cent] = new RooUnfoldResponse(Form("Resp_%i", cent), Form("Resp_%i", cent));
  resp[cent]->Setup(fMeas[cent], fTrue[cent]); //Setup Response Matrix Definition

  fMeas1D[cent] = new TH1D(Form("fMeas1D_Cent_%i", cent), Form("fMeas1D_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high);
  fTrue1D[cent] = new TH1D(Form("fTrue1D_Cent_%i", cent), Form("fTrue1D_Cent_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
  fMiss1D[cent] = new TH1D(Form("fMiss1D_Cent_%i", cent), Form("fMiss1D_Cent_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
  fMissRatio[cent] = new TH1D(Form("fMissRatio_Cent_%i", cent), Form("fMissRatio_Cent_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);
  fFake1D[cent] = new TH1D(Form("fFake1D_Cent_%i", cent), Form("fFake1D_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high);
  fFakeRatio[cent] = new TH1D(Form("fFakeRatio_Cent_%i", cent), Form("fFakeRatio_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high);
  fResp1D[cent] = new TH2D(Form("fResp1D_Cent_%i", cent), Form("fResp1D_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

  resp1D[cent] = new RooUnfoldResponse(Form("Resp1D_%i", cent), Form("Resp1D_%i", cent));
  resp1D[cent]->Setup(fMeas1D[cent], fTrue1D[cent]); //Setup Response Matrix Definition

  hMiss1D[cent] = new TH2D(Form("hMiss1D_Cent_%i", cent), Form("hMiss1D_Cent_%i", cent), njpt_bins, jetpt_low, jetpt_high, njpt_gen_bins, jetpt_gen_low, jetpt_gen_high);

  hDiffJetPt[cent] = new TH1D(Form("hDiffJetPt_%i", cent), Form("hDiffJetPt_%i", cent), 100, -50, 50);
  hDiffJetEta[cent] = new TH1D(Form("hDiffJetEta_%i", cent), Form("hDiffJetEta_%i", cent), 40, -1., 1.);
  hDiffJetEtaBeforeSmear[cent] = new TH1D(Form("hDiffJetEtaBeforeSmear_%i", cent), Form("hDiffJetEtaBeforeSmear_%i", cent), 40, -1., 1.);
  hDiffJetEtaAfterSmear[cent] = new TH1D(Form("hDiffJetEtaAfterSmear_%i", cent), Form("hDiffJetEtaAfterSmear_%i", cent), 40, -1., 1.);
  hDiffJetPhi[cent] = new TH1D(Form("hDiffJetPhi_%i", cent), Form("hDiffJetPhi_%i", cent), 40, -1., 1.);

  for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
    hDiffJetPtInGenPtBins[cent][ptbin] = new TH1D(Form("hDiffJetPt_%i_%i", cent, ptbin), Form("hDiffJetPt_%i_%i", cent, ptbin), 100, -50, 50);
  }

  hMCD0Pt[cent] = new TH1D(Form("hMCD0Pt_%i", cent), Form("hMCD0Pt_%i", cent), 20, 0, 10);
  hRecoD0Pt[cent] = new TH1D(Form("hRecoD0Pt_%i", cent), Form("hRecoD0Pt_%i", cent), 20, 0, 10);

  hJet[cent]->Sumw2();
  hJetZNorm[cent]->Sumw2();
  hJetZNormPtNorm[cent]->Sumw2();

  hJet[cent]->CalculateErrors();
  hJetZNorm[cent]->CalculateErrors();
  hJetZNormPtNorm[cent]->CalculateErrors();

  Weight[cent] = new TH2D(Form("Weight_%i", cent), Form("Weight_%i", cent), njpt_gen_bins, jetpt_gen_low, jetpt_gen_high, nz_gen_bins, z_gen_low, z_gen_high);

  for (int binx = 1; binx <= Weight[cent]->GetNbinsX(); binx++ ){
    for (int biny = 1; biny <= Weight[cent]->GetNbinsY(); biny++ ){
      Weight[cent]->SetBinContent(binx, biny, 1.);
    }
  }
  
  cout << JetTree->GetEntries() << endl;

  int nentries = JetTree->GetEntries();

  int lowlimit = 0;
  int highlimit = nentries;

  int counter = 0;

  double recooutside[3] = {0};

  double etasmearfactors[3] = {0.11, 0.09, 0.05}; // Eta Smear Factors From Heavy Ion Overlay
  double phismearfactors[3] = {0.15, 0.13 ,0.09}; // Phi Smear Factors From Heavy Ion Overlay

  TRandom* r = new TRandom(0);
  gRandom->SetSeed(0);

  TFile *PYTHIA = new TFile("PYTHIA_Pt_5_20.root", "READ");
  // TFile *PYTHIA = new TFile("/Volumes/WorkDrive/MattsOldFramework_Changed/Unfold/FONLLWeights_New.root", "READ");
  TH1F  *FONLLvPYTHIAWeights = (TH1F *) PYTHIA->Get("FONLLWeights");

  cout << "Mode == " << mode << endl;

  TH1D *tmppT = new TH1D("tmppT", "tmppT", 40, 0, 40);

  TH1D *plottingpT = new TH1D("plottingpT", "plottingpT", nbins_jpt, binning_jpt);

  for (int i = lowlimit; i < highlimit; i++){
    
    if (mode == 1) {if (i%2 != 0) continue;}
    else if (mode == 2) {if (i%2 == 0) continue;}    
    // cout << "Read Entry " << i << endl;

    JetTree->GetEntry(i);

    if (MCJetPt < 1) continue;

    int centhistogramtofill = cent;

    if (centhistogramtofill < 0) continue;

    // MCD0Phi = MCTrackPhi[MCD0Index];
    // RecoD0Phi = RecoTrackPhi[RecoD0Index];

    double MCJetPx = MCJetPt*TMath::Cos(MCJetPhi);
    double MCJetPy = MCJetPt*TMath::Sin(MCJetPhi);
    double MCD0Px = MCD0Pt*TMath::Cos(MCD0Phi);
    double MCD0Py = MCD0Pt*TMath::Sin(MCD0Phi);
    // double mcz = (MCJetPx*MCD0Px + MCJetPy*MCD0Py)/pow(MCJetPt, 2);
    double mcz = MCD0Pt/MCJetPt * TMath::Cos(MCJetPhi - MCD0Phi);

    // if (mcz > 1.0) cout << MCJetPt << "\t" << MCD0Pt << mcz << endl;

    if (mcz == 1.0) mcz = 0.999; // Padding the boundaries
    int ptbin = tmppT->FindBin(RecoJetPt);
    
    // cout << ptbin << "\t" << RecoJetPt << endl;
    RecoJetCorrPt = RecoJetPt + dPt[centhistogramtofill][ptbin-1]->GetRandom();

    double eta_mean_val = Eta_Mean[centhistogramtofill]->GetBinContent(ptbin);
    double eta_sigma_val = Eta_Sigma[centhistogramtofill]->GetBinContent(ptbin);

    double phi_mean_val = Phi_Mean[centhistogramtofill]->GetBinContent(ptbin);
    double phi_sigma_val = Phi_Sigma[centhistogramtofill]->GetBinContent(ptbin);

    // cout << "Eta Phi = " << eta_mean_val << eta_sigma_val << phi_mean_val << phi_sigma_val << endl;

    double oldRecoJetEta = RecoJetEta;

    RecoJetEta    = RecoJetEta + r->Gaus(eta_mean_val, eta_sigma_val);
    RecoJetPhi    = RecoJetPhi + r->Gaus(phi_mean_val, phi_sigma_val);

    if(RecoJetPhi>=2.*TMath::Pi()) RecoJetPhi = RecoJetPhi - 2.*TMath::Pi();
    if(RecoJetPhi<0) RecoJetPhi = 2.*TMath::Pi() + RecoJetPhi;

    double RecoJetPx = RecoJetCorrPt*TMath::Cos(RecoJetPhi);
    double RecoJetPy = RecoJetCorrPt*TMath::Sin(RecoJetPhi);
    double RecoD0Px = RecoD0Pt*TMath::Cos(RecoD0Phi);
    double RecoD0Py = RecoD0Pt*TMath::Sin(RecoD0Phi);
    // double recoz = (RecoJetPx*RecoD0Px + RecoJetPy*RecoD0Py)/pow(RecoJetCorrPt, 2);
    double recoz = RecoD0Pt/RecoJetCorrPt * TMath::Cos(RecoJetPhi - RecoD0Phi);

    double RecoDeltaR = dR(dEta(RecoJetEta, RecoD0Eta), dPhi(RecoJetPhi, RecoD0Phi));
    // cout << MCJetPx << "\t" << MCD0Px << "\t" << RecoJetPx << "\t" << RecoD0Px << endl;
    // cout << MCD0Pt << "\t" << MCJetPt << "\t" << MCJetPhi << "\t" << mcz << "\t" << RecoD0Pt << "\t" << RecoJetCorrPt << "\t" << RecoJetPhi << "\t" << recoz << endl;

    // if (recoz > z_low-0.001 && recoz <= z_low) recoz = z_low + 0.001; // Padding the boundaries
    // if (recoz == 1.0) recoz = 0.99;

    //Defining response matrix using fakes and misses to check closure

    bool isMCD0Pt = MCD0Pt > 1 && MCD0Pt < 10; // We wanna check 5 GeV in this folder.
    bool isRecoD0Pt = RecoD0Pt > 1 && RecoD0Pt < 10;
    bool isMCZ   = mcz > z_gen_low && mcz < z_gen_high; 
    bool isRecoZ = recoz > z_low && recoz < z_high; 
    bool isMCJetPt = MCJetPt > 5 && MCJetPt < 20;
    bool isRecoJetPt = RecoJetCorrPt > 3 && RecoJetCorrPt < 30; // 78% of jets with pT > 5 GeV is captured within this range
    bool isMCJetEta = abs(MCJetEta) < 0.601;
    bool isRecoJetEta = abs(RecoJetEta) < 0.6;
    bool isRecoDeltaR = RecoDeltaR < 0.4;

    // bool isConstGreaterThanJet = (RecoJetCorrPt - RecoD0Pt > 0);

    if (!isMCD0Pt) continue;
    if (!isRecoD0Pt) continue;
    if (RecoJetNConst==0) continue;
  
    isRecoZ = kTRUE;
    isMCZ = kTRUE;

    // bool isOK = isRecoZ && isMCZ && isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta && isRecoDeltaR;
    // bool isMISS = !isOK && isMCZ && isMCJetPt && isMCJetEta;
    // bool isFAKE = !isOK && isRecoZ && isRecoJetPt && isRecoJetEta && isRecoDeltaR;

    bool isOK   = isRecoJetPt && isMCJetPt && isRecoJetEta && isMCJetEta;
    bool isMISS = !isOK && isMCJetPt && isMCJetEta;
    bool isFAKE = !isOK && isRecoJetPt && isRecoJetEta;

    // double w = FONLLvPYTHIAWeights->GetBinContent(FONLLvPYTHIAWeights->FindBin(MCJetPt));
    double w = 1.;

    // cout << "W = " << w << endl;

    if (!isOK && !isMISS && !isFAKE) continue;

    // cout << MCJetPt << "\t" << MCJetEta << "\t" << MCJetPhi << "\t" << RecoJetCorrPt << "\t" << RecoJetEta << "\t" << RecoJetPhi << endl;

    // cout << Form("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", MCJetPt, MCJetEta, MCJetPhi, RecoJetCorrPt, RecoJetEta, RecoJetPhi);

    hDiffJetPt[centhistogramtofill]->Fill(RecoJetCorrPt - MCJetPt);
    hDiffJetEta[centhistogramtofill]->Fill(dEta(RecoJetEta, MCJetEta));
    hDiffJetEtaBeforeSmear[centhistogramtofill]->Fill(oldRecoJetEta - MCJetEta);
    hDiffJetEtaAfterSmear[centhistogramtofill]->Fill(RecoJetEta - MCJetEta);
    hDiffJetPhi[centhistogramtofill]->Fill(dPhi(RecoJetPhi, MCJetPhi));

    // int genptbin = plottingpT->FindBin(MCJetPt);

    // hDiffJetPtInGenPtBins[centhistogramtofill][genptbin-1]->Fill(RecoJetCorrPt - MCJetPt);

    hMCD0Pt[centhistogramtofill]->Fill(MCD0Pt, w);
    hRecoD0Pt[centhistogramtofill]->Fill(RecoD0Pt, w);

    if (isOK){
      resp[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);
      fResp1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);

      fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);
    }
    else if (isMISS){
      resp[centhistogramtofill]->Miss(MCJetPt, mcz, w);
      fMiss[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      resp1D[centhistogramtofill]->Miss(MCJetPt, w);
      fMiss1D[centhistogramtofill]->Fill(MCJetPt, w);

      fTrue[centhistogramtofill]->Fill(MCJetPt, mcz, w);
      fTrue1D[centhistogramtofill]->Fill(MCJetPt, w);
    }
    else if (isFAKE){
      resp[centhistogramtofill]->Fake(RecoJetCorrPt, recoz, w);
      resp1D[centhistogramtofill]->Fake(RecoJetCorrPt, w);

      fMeas[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fMeas1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);

      fFake[centhistogramtofill]->Fill(RecoJetCorrPt, recoz, w);
      fFake1D[centhistogramtofill]->Fill(RecoJetCorrPt, w);
    }



    if (isMISS) hMiss1D[centhistogramtofill]->Fill(RecoJetCorrPt, MCJetPt, w);

    if (!isMCJetPt) continue;
    if (!isRecoJetPt) continue;
    if (!isMCJetEta) continue;
    if (!isRecoJetEta) continue;
    if (!isMCZ) continue;

    if (!isRecoZ) recooutside[centhistogramtofill]++;

    double fvalue[ndim] = {MCJetPt, mcz, RecoJetCorrPt, recoz};
    hJet[centhistogramtofill]->Fill(fvalue);
    hZ[centhistogramtofill]->Fill(mcz);
    hJetPtvZ[centhistogramtofill]->Fill(MCJetPt, mcz);
  }

  TString filename;

  switch (mode){
    case 1:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    case 2:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    default:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins);
  }

  cout << filename.Data() << endl;
  TFile *outfile = new TFile(filename.Data(), "UPDATE");
  outfile->cd();

  hZ[cent]->Write();
  hJetPtvZ[cent]->Write();
  hJet[cent]->Write();

  fMeas[cent]->Write();
  fTrue[cent]->Write();
  fFake[cent]->Write();

  fMeas1D[cent]->Write();
  fTrue1D[cent]->Write();

  fMiss[cent]->Write();

  fResp1D[cent]->Write();
  fMiss1D[cent]->Write();
  fFake1D[cent]->Write();

  fMissRatio[cent]->Add(fMiss1D[cent]);
  fMissRatio[cent]->Divide(fTrue1D[cent]);

  fMissRatio[cent]->Write();

  fFakeRatio[cent]->Add(fFake1D[cent]);
  fFakeRatio[cent]->Divide(fMeas1D[cent]);

  fFakeRatio[cent]->Write();

  gDirectory->WriteObject(resp[cent], Form("Resp_%i", cent));
  gDirectory->WriteObject(resp1D[cent], Form("Resp1D_%i", cent));

  hMiss1D[cent]->Write();

  hDiffJetPt[cent]->Write();
  hDiffJetEta[cent]->Write();
  hDiffJetPhi[cent]->Write();

  hDiffJetEtaBeforeSmear[cent]->Write();
  hDiffJetEtaAfterSmear[cent]->Write();

  hMCD0Pt[cent]->Write();
  hRecoD0Pt[cent]->Write();

  // for (int cent = 0; cent < 3; cent++){
  //   for (int ptbin = 0; ptbin < nbins_jpt; ptbin++){
  //     hDiffJetPtInGenPtBins[cent][ptbin]->Write();
  //   }
  // }

  outfile->Close();
}

void PrepareEfficiencyHistogramsUsingMC(int SUPERITERATION = 0, int iteration = 3, TString DirName = "MCMCUnf", int mode = 0){

  TString filename;

  switch (mode){
    case 1:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    case 2:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i_Split_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins, mode);
      break;
    default:
      filename = Form("%s/OldResponse_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", DirName.Data(), 1, iteration, njpt_gen_bins, nz_gen_bins);
  }
  
  TFile *outfile = new TFile(filename.Data(), "RECREATE");
  outfile->Close();

  // Method(iteration, DirName, mode);

  for (int cent = 0; cent < 3; cent++){
    Method(cent, iteration, DirName, mode);
  }
}