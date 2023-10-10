// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// This code is set as an AnalysisMaker task, where it can perform:
// 1) jet analysis
// 	- tagging
// 	- jet-hadron correlations
//      - jet shape analysis
// 	- mixed events: use of an event pool to mix triggers with
//      - Rho (underlying event) subtraction to jets
//      - leading jet and subleading tag
//      - event plane calculation with BBC, ZDC, TPC
//      - event plane corrections with BBC, ZDC, TPC
//      - access to jet constituents
//      - general QA
//
// can get a pointer to:
// 1) collection of jets  	
// 2) event wise rho parameter
// 3) jet constituents (4 vectors)
// 4) leading + subleading jets
// 5) event plane calculation for various pt bins
// 6) externally read in Event Pool for event mixing
// ################################################################

#include "StJetHadronCorrelations.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <sstream>
#include <fstream>

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h" // NEW name
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"  // NEW (OLD: StPicoEmcPidTraits.h)

// jet-framework includes
#include "StEventPlaneMaker.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPoolManager.h"
#include "StFemtoTrack.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StJetHadronCorrelations)

//______________________________________________________________________________
StJetHadronCorrelations::StJetHadronCorrelations(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", bool mDoComments = kFALSE, double minJetPt = 1.0, double trkbias = 0.15, const char* jetMakerName = "", const char* rhoMakerName = "")
  : StJetFrameworkPicoBase(name)
{
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;
  doppAnalysis = kFALSE;
  doJetShapeAnalysis = kFALSE;
  fJetAnalysisJetType = kLeadingJets;  // see header, LJ, SubLJ, inclusive
  doRequireAjSelection = kFALSE;
  fCorrJetPt = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  doWriteTrackQAHist = kTRUE;
  doWriteJetQAHist = kTRUE;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fLeadingJet = 0x0; fSubLeadingJet = 0x0; fExcludeLeadingJetsFromFit = 1.0;
  fTrackWeight = 1; //StJetFrameworkPicoBase::kPtLinear2Const5Weight; // see StJetFrameworkPicoBase::EPtrackWeightType 
  fEventPlaneMaxTrackPtCut = 5.0;
  fTPCEPmethod = 1;
  fHistCentBinMin = 0;
  fHistCentBinMax = 9;               // 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
  fHistZvertBinMin = 0;
  fHistZvertBinMax = 20;             // (-40, 40) 4cm bins
  TPC_PSI2 = -999;
  TPCA_PSI2 = -999; // subevent A
  TPCB_PSI2 = -999; // subevent B
  BBC_PSI2 = -999; ZDC_PSI2 = -999; 
  BBC_PSI1 = -999; ZDC_PSI1 = -999;
  PSI2 = -999;
  RES = -999;
  TPC_raw_comb = 0.; TPC_raw_neg = 0.; TPC_raw_pos = 0.;
  BBC_raw_comb = 0.; BBC_raw_east = 0.; BBC_raw_west = 0.;
  ZDC_raw_comb = 0.; ZDC_raw_east = 0.; ZDC_raw_west = 0.;
  fPoolMgr = 0x0;
  fJets = 0x0;
  fRunNumber = 0;
  fEPTPCn = 0.; fEPTPCp = 0.;
  fEPTPC = 0.; fEPBBC = 0.; fEPZDC = 0.;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  mOutNameEP = "";
  mOutNameQA = "";
  mOutNameME = "";
  doPrintEventCounter = kFALSE;
  fDoEffCorr = kFALSE;
  doEventPlaneRes = kFALSE;
  doTPCptassocBin = kFALSE;
  fTPCptAssocBin = -99;
  doUseMainEPAngle = kFALSE;  // kTRUE: use 0.2-2.0 GeV charged tracks for event plane
  doRejectBadRuns = kFALSE;
  fMinPtJet = minJetPt;
  fJetConstituentCut = 2.0;
  fTrackBias = trkbias;
  fTowerBias = 0.2;
  fJetRad = 0.4;
  fLeadJetPtMin = 20.0; fSubLeadJetPtMin = 10.0;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0; fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0; fTowerPhiMaxCut = 2.0*TMath::Pi();
  fDoEventMixing = 0; fMixingTracks = 50000; fNMIXtracks = 5000; fNMIXevents = 5;
  fCentBinSize = 5; fCentBinSizeJS = 10; fReduceStatsCent = -1;
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
  //mVertex = 0x0;
  zVtx = 0.0;
  fDoFilterPtMixEvents = kFALSE;
  fDoUseMultBins = kFALSE;
  doUseEPBins = kFALSE;
  fnEPBins = 6;
  doIgnoreExternalME = kTRUE;
  fEmcTriggerEventType = 0; fMBEventType = 2; fMixingEventType = 0;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4800; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }
  doComments = mDoComments;
  fhnJH = 0x0;
  fhnMixedEvents = 0x0;
  fhnCorr = 0x0;
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  fEventPlaneMakerName = "";
  mBaseMaker = 0x0;
  mCentMaker = 0x0;
  //fEventPoolOutputList = {{}}; // FIXME - does this make sense as is?
  fUsePtBinnedEventPool = kFALSE;
  fCheckEventNumberInMixedEvent = kFALSE;
  fListOfPools = 0x0;
}
//
//
//__________________________________________________________________________________
StJetHadronCorrelations::~StJetHadronCorrelations()
{ /*  */
  // destructor
  if(hdEPReactionPlaneFnc) delete hdEPReactionPlaneFnc;
  if(hdEPEventPlaneFncN2)  delete hdEPEventPlaneFncN2;
  if(hdEPEventPlaneFncP2)  delete hdEPEventPlaneFncP2;
  if(hdEPEventPlaneFnc2)   delete hdEPEventPlaneFnc2;
  if(hdEPEventPlaneClass)  delete hdEPEventPlaneClass;
  if(hReactionPlaneFnc)    delete hReactionPlaneFnc;
  if(hEventPlaneFncN2)     delete hEventPlaneFncN2;
  if(hEventPlaneFncP2)     delete hEventPlaneFncP2;
  if(hEventPlaneFnc2)      delete hEventPlaneFnc2;
  if(hEventPlaneClass)     delete hEventPlaneClass;

  if(hEventPlane)   delete hEventPlane;
  if(fHistEPTPCn)   delete fHistEPTPCn;
  if(fHistEPTPCp)   delete fHistEPTPCp;
  if(fHistEPBBC)    delete fHistEPBBC;
  if(fHistEPZDC)    delete fHistEPZDC;
  if(hEventZVertex) delete hEventZVertex;
  if(hCentrality)   delete hCentrality;
  if(hMultiplicity) delete hMultiplicity;
  if(hStats)        delete hStats;
  if(hRhovsCent)    delete hRhovsCent;
  for(int i=0; i<5; i++) { 
    if(hdEPtrk[i])  delete hdEPtrk[i]; 
  }
  for(int i=0; i<9; i++){ // centrality
    if(hTrackPhi[i]) delete hTrackPhi[i];
    if(hTrackEta[i]) delete hTrackEta[i];
    if(hTrackPt[i])  delete hTrackPt[i];
  }
  if(hTrackEtavsPhi) delete hTrackEtavsPhi;

  if(hJetPt)             delete hJetPt;
  if(hJetCorrPt)         delete hJetCorrPt;
  if(hJetLeadingPt)      delete hJetLeadingPt;
  if(hJetSubLeadingPt)   delete hJetSubLeadingPt;
  if(hJetLeadingPtAj)    delete hJetLeadingPtAj;
  if(hJetSubLeadingPtAj) delete hJetSubLeadingPtAj;
  if(hJetDiJetAj)        delete hJetDiJetAj;
  if(hJetE)              delete hJetE;
  if(hJetEta)            delete hJetEta;
  if(hJetPhi)            delete hJetPhi;
  if(hJetNEF)            delete hJetNEF;
  if(hJetArea)           delete hJetArea;
  if(hJetTracksPt)       delete hJetTracksPt;
  if(hJetTracksPhi)      delete hJetTracksPhi;
  if(hJetTracksEta)      delete hJetTracksEta;
  if(hJetTracksZ)        delete hJetTracksZ;
  if(hJetPtvsArea)       delete hJetPtvsArea;
  if(hJetEventEP)        delete hJetEventEP;
  if(hJetPhivsEP)        delete hJetPhivsEP;

  if(hJetPtIn)       delete hJetPtIn;
  if(hJetPhiIn)      delete hJetPhiIn;
  if(hJetEtaIn)      delete hJetEtaIn;
  if(hJetEventEPIn)  delete hJetEventEPIn;
  if(hJetPhivsEPIn)  delete hJetPhivsEPIn;
  if(hJetPtMid)      delete hJetPtMid;
  if(hJetPhiMid)     delete hJetPhiMid;
  if(hJetEtaMid)     delete hJetEtaMid;
  if(hJetEventEPMid) delete hJetEventEPMid;
  if(hJetPhivsEPMid) delete hJetPhivsEPMid;
  if(hJetPtOut)      delete hJetPtOut;
  if(hJetPhiOut)     delete hJetPhiOut;
  if(hJetEtaOut)     delete hJetEtaOut;
  if(hJetEventEPOut) delete hJetEventEPOut;
  if(hJetPhivsEPOut) delete hJetPhivsEPOut;

  if(fHistJetHEtaPhi)        delete fHistJetHEtaPhi;
  if(fHistEventSelectionQA)  delete fHistEventSelectionQA;
  if(fHistEventSelectionQAafterCuts) delete fHistEventSelectionQAafterCuts;
  if(hTriggerIds)            delete hTriggerIds;
  if(hEmcTriggers)           delete hEmcTriggers;
  if(hBadTowerFiredTrigger)  delete hBadTowerFiredTrigger;
  if(hMixEvtStatZVtx)        delete hMixEvtStatZVtx;
  if(hMixEvtStatCent)        delete hMixEvtStatCent;
  if(hMixEvtStatZvsCent)     delete hMixEvtStatZvsCent;
  if(hTriggerEvtStatZVtx)    delete hTriggerEvtStatZVtx;
  if(hTriggerEvtStatCent)    delete hTriggerEvtStatCent;
  if(hTriggerEvtStatZvsCent) delete hTriggerEvtStatZvsCent;
  if(hMBvsMult)              delete hMBvsMult;
  if(hMB5vsMult)             delete hMB5vsMult;
  if(hMB30vsMult)            delete hMB30vsMult;
  if(hHTvsMult)              delete hHTvsMult;
  if(hNMixEvents)            delete hNMixEvents;

  if(hTPCvsBBCep) delete hTPCvsBBCep;
  if(hTPCvsZDCep) delete hTPCvsZDCep;
  if(hBBCvsZDCep) delete hBBCvsZDCep;

  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        if(fProfJetV2[k][j][i]) delete fProfJetV2[k][j][i];
      }
    }
  }

  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      if(fProfV2Resolution[i]) delete fProfV2Resolution[i];
      if(fProfV3Resolution[i]) delete fProfV3Resolution[i];
      if(fProfV4Resolution[i]) delete fProfV4Resolution[i];
      if(fProfV5Resolution[i]) delete fProfV5Resolution[i];
    }
  }

  // delete sparses
  if(fhnJH)          delete fhnJH;
  if(fhnMixedEvents) delete fhnMixedEvents;
  if(fhnCorr)        delete fhnCorr;

  // clear and delete objects
//  fJets->Clear(); delete fJets;
//  fRho->Clear(); delete fRho; 
//  fPoolMgr->Clear(); delete fPoolMgr;

  // Clear unnecessary pools before saving - FIXME
  fPoolMgr->ClearPools();
  if(fListOfPools) delete fListOfPools;
}
//
// initialize objects & set up
//_________________________________________________________________________________________
Int_t StJetHadronCorrelations::Init() {
  //StJetFrameworkPicoBase::Init();

  // Initialize output list of event pools - TEST
  if (fListOfPools != NULL){
        delete fListOfPools;
        fListOfPools = NULL;
  }
  if (!fListOfPools){
        fListOfPools = new TList();
        fListOfPools->SetOwner(kTRUE);
  }

  // initialize the histograms
  DeclareHistograms();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);

  return kStOK;
}
//
//
//______________________________________________________________________________________
Int_t StJetHadronCorrelations::Finish() { 
  //  Summarize the run.
  cout << "StJetHadronCorrelations::Finish()\n";

  // Write event pool manager object to file and close it
  if(mOutNameME!="") {
    TFile *fOutME = new TFile(mOutNameME.Data(), "RECREATE");
    fOutME->cd();
    fListOfPools->Write();  // write pools to file if we have them and want to
    fOutME->Write();
    fOutME->Close();
  }

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    //fout->ls();
    fout->cd();
    fout->mkdir(fAnalysisMakerName);
    fout->cd(fAnalysisMakerName);
    WriteHistograms();
   
    fout->cd();
    fout->Write();
    fout->Close();
  }

/*
  //  Write QA histos to file and close it.
  if(mOutNameQA!=""  && fJetShapePtAssocBin < 5) {
    TFile *fQAout = new TFile(mOutNameQA.Data(), "UPDATE");
    fQAout->cd();

    // track QA
    if(doWriteTrackQAHist && (fJetShapePtAssocBin < 5)) {
      fQAout->mkdir(Form("TrackQA_bin%i", fTPCptAssocBin));
      fQAout->cd(Form("TrackQA_bin%i", fTPCptAssocBin));
      WriteTrackQAHistograms();
      fQAout->cd();
    }

    // jet QA
    if(doWriteJetQAHist && (fJetShapePtAssocBin < 5)) {
      fQAout->mkdir(Form("JetEPQA_bin%i", fTPCptAssocBin));
      fQAout->cd(Form("JetEPQA_bin%i", fTPCptAssocBin));
      WriteJetEPQAHistograms();
      fQAout->cd();
    }

    fQAout->Write();
    fQAout->Close();
  }
*/

  cout<<"End of StJetHadronCorrelations::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}
//
// function to declare histograms and set up global objects
//__________________________________________________________________________________________
void StJetHadronCorrelations::DeclareHistograms() {
  // constants
  double pi = 1.0*TMath::Pi();

  // binning for cent histograms
  int nHistCentBins = 0;
  if(fCentBinSize == 10) nHistCentBins = 10;
  if(fCentBinSize ==  5) nHistCentBins = 20;

  // binning for mult histograms:    pp : AuAu
  double kHistMultMax = (doppAnalysis) ? 100. : 800.;
  int kHistMultBins = (doppAnalysis) ? 100 : 400;

  // QA histos
  hdEPReactionPlaneFnc = new TH1F("hdEPReactionPlaneFnc", "jets relative to EP from reaction plane function", 3, 0.0, 0.5*pi);
  hdEPEventPlaneFncN2 = new TH1F("hdEPEventPlaneFncN2", "jets relative to EP from event plane function SUB1 meth2", 3, 0.0, 0.5*pi);
  hdEPEventPlaneFncP2 = new TH1F("hdEPEventPlaneFncP2", "jets relative to EP from event plane function SUB2 meth2", 3, 0.0, 0.5*pi);
  hdEPEventPlaneFnc2 = new TH1F("hdEPEventPlaneFnc2", "jets relative to EP from event plane function meth2", 3, 0.0, 0.5*pi);
  hdEPEventPlaneClass = new TH1F("hdEPEventPlaneClass", "jets relative to EP from event plane class", 3, 0.0, 0.5*pi);
  hReactionPlaneFnc = new TH1F("hReactionPlaneFnc", "Event plane distribution from reaction plane function", 72, 0.0, 1.0*pi);
  hEventPlaneFncN2 = new TH1F("hEventPlaneFncN2", "Event plane distribution from event plane function SUB1 meth2", 72, 0.0, 1.0*pi);
  hEventPlaneFncP2 = new TH1F("hEventPlaneFncP2", "Event plane distribution from event plane function SUB2 meth2", 72, 0.0, 1.0*pi);
  hEventPlaneFnc2 = new TH1F("hEventPlaneFnc2", "Event plane distribution from event plane function meth2", 72, 0.0, 1.0*pi);
  hEventPlaneClass = new TH1F("hEventPlaneClass", "Event plane distribution from event plane class", 72, 0.0, 1.0*pi);

  hEventPlane = new TH1F("hEventPlane", "Event plane distribution", 72, 0.0, 1.0*pi);
  fHistEPTPCn = new TH2F("fHistEPTPCn", "", 20, 0., 100., 72, -pi, pi);
  fHistEPTPCp = new TH2F("fHistEPTPCp", "", 20, 0., 100., 72, -pi, pi);
  fHistEPBBC = new TH2F("fHistEPBBC", "", 20, 0., 100., 72, -pi, pi);
  fHistEPZDC = new TH2F("fHistEPZDC", "", 20, 0., 100., 72, -pi, pi);
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 100, -50, 50);
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100); 
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  hStats = new TH1F("hStats", "QA stats for cuts and objects", 25, 0.5, 25.5);
  hStats->GetXaxis()->SetBinLabel(1, "total events");
  hStats->GetXaxis()->SetBinLabel(2, "good run event");
  hStats->GetXaxis()->SetBinLabel(3, "passed max trk pt cut");
  hStats->GetXaxis()->SetBinLabel(4, "passed z-vtx cut");
  hStats->GetXaxis()->SetBinLabel(5, "passed cent cut");
  hStats->GetXaxis()->SetBinLabel(6, "post cent selection");
  hStats->GetXaxis()->SetBinLabel(7, "have jet object");
  hStats->GetXaxis()->SetBinLabel(8, "have jets");
  hStats->GetXaxis()->SetBinLabel(9, "have rho object");
  hStats->GetXaxis()->SetBinLabel(10, "");
  hStats->GetXaxis()->SetBinLabel(11, "do jet shape analysis");
  hStats->GetXaxis()->SetBinLabel(12, "");
  hStats->GetXaxis()->SetBinLabel(13, "");
  hStats->GetXaxis()->SetBinLabel(14, "");
  hStats->GetXaxis()->SetBinLabel(15, "");
  hStats->GetXaxis()->SetBinLabel(16, "");
  hStats->GetXaxis()->SetBinLabel(17, "");
  hStats->GetXaxis()->SetBinLabel(18, "");
  hStats->GetXaxis()->SetBinLabel(19, "");
  hStats->GetXaxis()->SetBinLabel(20, "");

  hRhovsCent = new TH2F("hRhovsCent", "#rho vs centrality", 20, 0, 100, 200, 0, 200);

  for(int i=0; i<5; i++) { // pt bins
    hdEPtrk[i] = new TH1F(Form("hdEPtrk%i", i), Form("tracks relative to event plane, p_{T} bin=%i", i), 72, 0, 0.5*pi);
  }

  // track phi distribution for centrality
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i] = new TH1F(Form("hTrackPhi%d", i), Form("track distribution vs #phi, centr%d", i), 144, 0, 2.0*pi);
    hTrackEta[i] = new TH1F(Form("hTrackEta%d", i), Form("track distribution vs #eta, centr%d", i), 40, -1.0, 1.0);
    hTrackPt[i] = new TH1F(Form("hTrackPt%d", i), Form("track distribution vs p_{T}, centr%d", i), 120, 0., 30.0);
  }
  hTrackEtavsPhi = new TH2F(Form("hTrackEtavsPhi"), Form("track distribution: #eta vs #phi"), 144, 0, 2.0*pi, 40, -1.0, 1.0);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);
  hJetLeadingPt = new TH1F("hJetLeadingPt", "Leading Jet p_{T}", 100, 0, 100);
  hJetSubLeadingPt = new TH1F("hJetSubLeadingPt", "SubLeading Jet p_{T}", 100, 0, 100);
  hJetLeadingPtAj = new TH1F("hJetLeadingPtAj", "Leading Jet p_{T} with Aj cut", 100, 0, 100);
  hJetSubLeadingPtAj = new TH1F("hJetSubLeadingPtAj", "SubLeading Jet p_{T} with Aj cut", 100, 0, 100);
  hJetDiJetAj = new TH1F("hJetDiJetAj", "DiJet imbalance: Aj", 100, 0., 1.);
  hJetE = new TH1F("hJetE", "Jet energy distribution", 100, 0, 100);
  hJetEta = new TH1F("hJetEta", "Jet #eta distribution", 24, -1.2, 1.2);
  hJetPhi = new TH1F("hJetPhi", "Jet #phi distribution", 72, 0.0, 2.0*pi);
  hJetNEF = new TH1F("hJetNEF", "Jet NEF", 100, 0, 1);
  hJetArea = new TH1F("hJetArea", "Jet Area", 100, 0, 1);
  hJetTracksPt = new TH1F("hJetTracksPt", "Jet track constituent p_{T}", 120, 0, 30.0);
  hJetTracksPhi = new TH1F("hJetTracksPhi", "Jet track constituent #phi", 72, 0, 2.0*pi);
  hJetTracksEta = new TH1F("hJetTracksEta", "Jet track constituent #eta", 56, -1.4, 1.4);
  hJetTracksZ = new TH1F("hJetTracksZ", "Jet track fragmentation function", 144, 0, 1.44);
  hJetPtvsArea = new TH2F("hJetPtvsArea", "Jet p_{T} vs Jet area", 100, 0, 100, 100, 0, 1);
  hJetEventEP = new TH1F("hJetEventEP", "no of jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEP = new TH2F("hJetPhivsEP", "Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);

  // jet EP QA histos
  hJetPtIn = new TH1F("hJetPtIn", "no of jets in-plane vs p_{T}", 100, 0.0, 100);
  hJetPhiIn = new TH1F("hJetPhiIn", "no of jets in-plane vs #phi", 72, 0.0, 2.0*pi);
  hJetEtaIn = new TH1F("hJetEtaIn", "no of jets in-plane vs #eta", 40, -1.0, 1.0);
  hJetEventEPIn = new TH1F("hJetEventEPIn", "no of in-plane jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEPIn = new TH2F("hJetPhivsEPIn", "in-plane Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);
  hJetPtMid = new TH1F("hJetPtMid", "no of jets mid-plane vs p_{T}", 100, 0.0, 100);
  hJetPhiMid = new TH1F("hJetPhiMid", "no of jets mid-plane vs #phi", 72, 0.0, 2.0*pi);
  hJetEtaMid = new TH1F("hJetEtaMid", "no of jets mid-plane vs #eta", 40, -1.0, 1.0);
  hJetEventEPMid = new TH1F("hJetEventEPMid", "no of mid-plane jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEPMid = new TH2F("hJetPhivsEPMid", "mid-plane Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);
  hJetPtOut = new TH1F("hJetPtOut", "no of jets out-of-plane vs p_{T}", 100, 0.0, 100);
  hJetPhiOut = new TH1F("hJetPhiOut", "no of jets out-of-plane vs #phi", 72, 0.0, 2.0*pi);
  hJetEtaOut = new TH1F("hJetEtaOut", "no of jets out-of-plane vs #eta", 40, -1.0, 1.0);
  hJetEventEPOut = new TH1F("hJetEventEPOut", "no of out-of-plane jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEPOut = new TH2F("hJetPhivsEPOut", "out-of-plane Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);

  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi", "Jet-hadron #Delta#eta-#Delta#phi", 72, -1.8, 1.8, 72, -0.5*pi, 1.5*pi);

  // Event Selection QA histo + event QA
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
  hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
  hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);
  hBadTowerFiredTrigger = new TH1F("hBadTowerFiredTrigger", "# bad tower fired trigger vs tower Id", 4800, 0.5, 4800.5);
  hMixEvtStatZVtx = new TH1F("hMixEvtStatZVtx", "no of events in pool vs zvtx", 20, -40.0, 40.0);
  hMixEvtStatCent = new TH1F("hMixEvtStatCent", "no of events in pool vs Centrality", nHistCentBins, 0, 100);
  hMixEvtStatZvsCent = new TH2F("hMixEvtStatZvsCent", "no of events: zvtx vs Centality", nHistCentBins, 0, 100, 20, -40.0, 40.0);
  hTriggerEvtStatZVtx = new TH1F("hTriggerEvtStatZVtx", "no of trigger events used vs zvtx", 20, -40.0, 40.0);
  hTriggerEvtStatCent = new TH1F("hTriggerEvtStatCent", "no of trigger events used vs Centrality", nHistCentBins, 0, 100);
  hTriggerEvtStatZvsCent = new TH2F("hTriggerEvtStatZvsCent", "no of trigger events used: zvtx vs Centality", nHistCentBins, 0, 100, 20, -40.0, 40.0);
  hMBvsMult = new TH1F("hMBvsMult", "# MB events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hMB5vsMult = new TH1F("hMB5vsMult", "# MB5 events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hMB30vsMult = new TH1F("hMB30vsMult", "# MB30 events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hHTvsMult = new TH1F("hHTvsMult", "# HT events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hNMixEvents = new TH1F("hNMixEvents", "number of mixing events", 200, 0, 200);

  //// res_cen=new TProfile("res_cen","res vs. cen",10,0,10,-2,2);
  // 2-D event plane differences - (updated ranges on February 7th)
  hTPCvsBBCep = new TH2F("hTPCvsBBCep", "TPC vs BBC 2nd order event plane", 144, 0.*pi, 1.0*pi, 144, 0.*pi, 1.*pi);
  hTPCvsZDCep = new TH2F("hTPCvsZDCep", "TPC vs ZDC 2nd order event plane", 144, 0.*pi, 1.0*pi, 144, 0.*pi, 1.*pi);
  hBBCvsZDCep = new TH2F("hBBCvsZDCep", "BBC vs ZDC 2nd order event plane", 144, 0.*pi, 1.0*pi, 144, 0.*pi, 1.*pi);

  // jet v2 histogram - primarily used for EP resolution unfolding
  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        fProfJetV2[k][j][i] = new TProfile(Form("fProfJetV2_%i_%i_%i", k, j, i), Form("fProfJetV2_%i_%i_%i", k, j, i), 11, -0.5, 10.5);
      }
    }
  }

  if(doEventPlaneRes){
    // Reaction Plane resolution as function of centrality - corrected for 2nd order event plane
    for (Int_t i=0; i<9; i++){
      // 2nd order correction
      fProfV2Resolution[i] = new TProfile(Form("fProfV2Resolution_%i", i), Form("fProfV2Resolution_%i", i), 25, 0.5, 25.5);
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(2(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(2(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(2(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(2(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(2(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(2(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(2(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(2(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(2(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(2(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(2(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(2(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(2(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(2(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(2(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(2(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 3rd order correction
      fProfV3Resolution[i] = new TProfile(Form("fProfV3Resolution_%i", i), Form("fProfV3Resolution_%i", i), 25, 0.5, 25.5);
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(3(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(3(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(3(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(3(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(3(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(3(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(3(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(3(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(3(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(3(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(3(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(3(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(3(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(3(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(3(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(3(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(3(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(3(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 4th order correction
      fProfV4Resolution[i] = new TProfile(Form("fProfV4Resolution_%i", i), Form("fProfV4Resolution_%i", i), 25, 0.5, 25.5);
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(4(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(4(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(4(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(4(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(4(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(4(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(4(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(4(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(4(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(4(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(4(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(4(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(4(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(4(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(4(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(4(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(4(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(4(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 5th order correction
      fProfV5Resolution[i] = new TProfile(Form("fProfV5Resolution_%i", i), Form("fProfV5Resolution_%i", i), 25, 0.5, 25.5);
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(5(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(5(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(5(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(5(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(5(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(5(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(5(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(5(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(5(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(5(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(5(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(5(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(5(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(5(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(5(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(5(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(5(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(5(#Psi_{ZDC1} - #Psi_{TPCp}))>");
    }
  }

  // set up jet-hadron sparse
  UInt_t bitcodeMESE = 0; // bit coded, see GetDimParams() below
  bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; // | 1<<8 | 1<<9 | 1<<10;
  fhnJH = NewTHnSparseF("fhnJH", bitcodeMESE);

  // =================================================================================================
  // set up centrality bins for mixed events
  // for pp we need mult bins for event mixing. Create binning here, to also make a histogram from it

  // Setup for Au-Au collisions: cent bin size can only be 5 or 10% bins
  // centrality bins for mixed events
  int nCentralityBinsJS = 100 + 1;
  double multJS = 1.0;
  if(fCentBinSizeJS==5){ // will be most commonly used
    nCentralityBinsJS = 20;
    multJS = 5.0;
  } else if(fCentBinSizeJS==10){
    nCentralityBinsJS = 10;
    multJS = 10.0;
  } else if(fCentBinSizeJS==20){
    nCentralityBinsJS = 5;
    multJS = 20.0;
  }

  // this is temp as the above and various other implementation attempts would not work for both cases
  // need to look into this, but a few hours already is too much.  Don't want to have to have this hard-coded
  Int_t nCentBins = 8;
  Double_t* centralityBins = GenerateFixedBinArray(nCentBins, 0., 8.);

  // for AuAu data
  // +1 to accomodate the fact that we define bins rather than array entries
  const int nMultBinsJS = 26;  // Alt-2 - 27 values, 26 ranges
  Double_t multBinsJS[nMultBinsJS + 1] = {0, 10,15,21,31,42,53,66,   80, 95, 112, 130, 149, 169, 190, 212, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469, 800};
  Double_t *multiplicityBinsJS = multBinsJS;

  // cent bins for AuAu data
  Int_t nCentBinsJS = 20;
  Double_t* centralityBinsJSnew = GenerateFixedBinArray(nCentBinsJS, 0., 20.); 

  // for pp data
  const int nMultBinsJSpp = 7;
  //Double_t multBinsJSpp[nMultBinsJSpp + 1] = {0.0, 4., 9, 15, 25, 35, 55, 100.0, 500.0};  // 8 (9)
  Double_t multBinsJSpp[nMultBinsJSpp + 1] = {0.0, 4.0, 6.0, 8.0, 10.0, 13.0, 30., 100.};   // 7 (8)
  Double_t *multiplicityBinsJSpp = multBinsJSpp;

  // z-vertex bins for mixed events
  Int_t nZvBins  = 20; // 4 cm wide, 40 for 2 cm wide
  Double_t* zvBins = GenerateFixedBinArray(nZvBins, -40., 40.); // min/max doesn't matter as data is cut zmin/zmax

  // event plane bins for mixed events
  Int_t nEPBins = 6;   // default 6 from 0-180 degrees, (0 - pi) range
  if(fnEPBins != 6) nEPBins = fnEPBins;
  Double_t* epBins = GenerateFixedBinArray(nEPBins, 0., 1.0*pi);

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implementation of AliEventPoolManager
  if(doIgnoreExternalME) {
    if(fDoUseMultBins) {
      if(!doppAnalysis && !doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins); // not pp
      if(!doppAnalysis &&  doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins); // not pp
      if( doppAnalysis) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJSpp, (Double_t*)multiplicityBinsJSpp, nZvBins, (Double_t*)zvBins); // is pp

    } else { // centrality binning  
      if(!doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsJS, (Double_t*)centralityBinsJSnew, nZvBins, (Double_t*)zvBins);
      if( doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsJS, (Double_t*)centralityBinsJSnew, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins);

    } // correlation analysis setup
  }

  // ================================================================================================
  // TEMP testing 
  if(!doIgnoreExternalME) {
    Int_t nPtBins = 1;
    Double_t defaultPtBins[2] = {-9999., 9999.};
    Double_t* ptBins = defaultPtBins;

    // check if fPoolMgr is set, else construct object
    if(!fPoolMgr) {
      cout<<"Don't have pre-existing fPoolMgr object, We will now create it!"<<endl;
      fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins, nEPBins, epBins, nPtBins, ptBins);
    } else {
      cout<<"User has provided fPoolMgr object.. goodluck!"<<endl;
    }

    // If some bins of the pool should be saved, fEventPoolOutputList must be given using AddEventPoolsToOutput()
    // Note that this is in principle also possible, if an external poolmanager was given
    for(UInt_t i = 0; i < fEventPoolOutputList.size(); i++) {
      Double_t minCent = fEventPoolOutputList[i][0];
      Double_t maxCent = fEventPoolOutputList[i][1];
      Double_t minZvtx = fEventPoolOutputList[i][2];
      Double_t maxZvtx = fEventPoolOutputList[i][3];
      Double_t minPsi2 = fEventPoolOutputList[i][4];
      Double_t maxPsi2 = fEventPoolOutputList[i][5];
      Double_t minPt   = fEventPoolOutputList[i][6];
      Double_t maxPt   = fEventPoolOutputList[i][7];

      fPoolMgr->SetSaveFlag(minCent, maxCent, minZvtx, maxZvtx, minPsi2, maxPsi2, minPt, maxPt);
    }

    // Basic checks and printing of pool properties
    fPoolMgr->Validate();

    // save to output if requested
    if(fEventPoolOutputList.size()) fListOfPools->Add(fPoolMgr);

  } // external use switch
  // ================================================================================================

  // set up event mixing sparse
  bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; // | 1<<8 | 1<<9;
  fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", bitcodeMESE);

  // jet counter for normalizations in correlation analysis
  UInt_t bitcodeCorr = 0; // bit coded, see GetDimparamsCorr() below
  bitcodeCorr = 1<<0 | 1<<1 | 1<<2 | 1<<3; // | 1<<4;
  fhnCorr = NewTHnSparseFCorr("fhnCorr", bitcodeCorr);

  // Switch on Sumw2 for all histos - (except profiles)
  SetSumw2();
}
//
// write track QA histograms
//_____________________________________________________________________________
void StJetHadronCorrelations::WriteTrackQAHistograms() {
  // track phi distribution for centrality
  for(int i=0; i<5; i++) { hdEPtrk[i]->Write(); }
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i]->Write();
    hTrackEta[i]->Write();
    hTrackPt[i]->Write();
  }
  hTrackEtavsPhi->Write();
}
//
// write Jet event plane QA histograms
//______________________________________________________________________________
void StJetHadronCorrelations::WriteJetEPQAHistograms() {
  hdEPReactionPlaneFnc->Write();
  hdEPEventPlaneFncN2->Write();
  hdEPEventPlaneFncP2->Write();
  hdEPEventPlaneFnc2->Write();
  hdEPEventPlaneClass->Write();
  hReactionPlaneFnc->Write();
  hEventPlaneFncN2->Write();
  hEventPlaneFncP2->Write();
  hEventPlaneFnc2->Write();
  hEventPlaneClass->Write();

  hJetPtIn->Write();
  hJetPhiIn->Write();
  hJetEtaIn->Write();
  hJetEventEPIn->Write();
  hJetPhivsEPIn->Write();
  hJetPtMid->Write();
  hJetPhiMid->Write();
  hJetEtaMid->Write();
  hJetEventEPMid->Write();
  hJetPhivsEPMid->Write();
  hJetPtOut->Write();
  hJetPhiOut->Write();
  hJetEtaOut->Write();
  hJetEventEPOut->Write();
  hJetPhivsEPOut->Write();
}
//
// write histograms
//_____________________________________________________________________________
void StJetHadronCorrelations::WriteHistograms() {
  // default histos
  hEventPlane->Write();
  fHistEPTPCn->Write();
  fHistEPTPCp->Write();
  fHistEPBBC->Write();
  fHistEPZDC->Write();
  hEventZVertex->Write();
  hCentrality->Write();
  hMultiplicity->Write();
  hStats->Write();
  hRhovsCent->Write();

  // jet histos
  hJetPt->Write();
  hJetCorrPt->Write();
  hJetLeadingPt->Write();
  hJetSubLeadingPt->Write();
  hJetLeadingPtAj->Write();
  hJetSubLeadingPtAj->Write();
  hJetDiJetAj->Write();
  hJetE->Write();
  hJetEta->Write();
  hJetPhi->Write();
  hJetNEF->Write();
  hJetArea->Write();
  hJetTracksPt->Write();
  hJetTracksPhi->Write();
  hJetTracksEta->Write();
  hJetTracksZ->Write();
  hJetPtvsArea->Write();
  hJetEventEP->Write();
  hJetPhivsEP->Write();

  fHistJetHEtaPhi->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();
  hBadTowerFiredTrigger->Write();
  hMixEvtStatZVtx->Write();
  hMixEvtStatCent->Write();
  hMixEvtStatZvsCent->Write();
  hTriggerEvtStatZVtx->Write();
  hTriggerEvtStatCent->Write();
  hTriggerEvtStatZvsCent->Write();
  hMBvsMult->Write();
  hMB5vsMult->Write();
  hMB30vsMult->Write();
  hHTvsMult->Write();
  hNMixEvents->Write();

  fhnJH->Write();
  fhnMixedEvents->Write();
  fhnCorr->Write();

  // (perhaps temp - save resolution hists to main output file instead of event plane calibration file)
  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      fProfV2Resolution[i]->Write();
      fProfV3Resolution[i]->Write();
      fProfV4Resolution[i]->Write();
      fProfV5Resolution[i]->Write();
    }
  }

  // event plane method: jet v2 histograms
  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        fProfJetV2[k][j][i]->Write();
      }
    }
  }

}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StJetHadronCorrelations::Clear(Option_t *opt) {
  fJets->Clear();
}
// 
//  This method is called every event.
//_____________________________________________________________________________
Int_t StJetHadronCorrelations::Make() {
  hStats->Fill(1);

  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // constants
  const double pi = 1.0*TMath::Pi();
  //StMemStat::PrintMem("MyAnalysisMaker at beginning of make");

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent 
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
    hStats->Fill(2);
  }

  // cut event on max track pt > 35.0 GeV (lowered to 30 on Oct 25, 2018)
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;
  hStats->Fill(3);

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut - the Aj analysis cut on (-40, 40) for reference
  // cut on (-30, 30) when using NEW centrality definitions - perhaps cut on (-28, 28)
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  hEventZVertex->Fill(zVtx);
  hStats->Fill(4);

  // get the fill, and event ID
  int fillId = mPicoEvent->fillId();
  int eventId = mPicoEvent->eventId();
  if(fDebugLevel == kDebugGeneralEvt) cout<<"RunID = "<<fRunNumber<<"  fillID = "<<fillId<<"  eventID = "<<eventId<<endl; // what is eventID?

  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 
  hStats->Fill(5);
 
  // bin-age to use for mixed event and sparses
  Int_t centbin10 = GetCentBin10(centbin);
  double centBinToUse;
  if(fCentBinSize==10)       { centBinToUse = (double)centbin10 * 10.0;
  } else if(fCentBinSize==5) { centBinToUse = (double)centbin * 5.0; }

  // centrality / multiplicity histograms
  // event activity - compensate for pp or AuAu
  double kEventActivity = (doppAnalysis) ? (double)grefMult : refCorr2;
  hMultiplicity->Fill(kEventActivity);
  hCentrality->Fill(fCentralityScaled);

  // to limit filling unused entries in sparse, only fill for certain centrality ranges
  // ranges can be different than functional cent bin setter
  Int_t cbin = -1;
  // this is actually not used since the below line does the cut:  if(fRequireCentSelection)
  if (centbin>-1 && centbin < 2)    cbin = 1; // 0-10%
  else if (centbin>1 && centbin<4)  cbin = 2; // 10-20%
  else if (centbin>3 && centbin<6)  cbin = 3; // 20-30%
  else if (centbin>5 && centbin<10) cbin = 4; // 30-50%
  else if (centbin>9 && centbin<16) cbin = 5; // 50-80%
  else cbin = -99;

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }
  hStats->Fill(6);
  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggersHist(hEmcTriggers);

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds(); 
  if(fDebugLevel == kDebugEmcTrigger) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", "; 
  }
  if(fDebugLevel == kDebugEmcTrigger) cout<<endl;

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30); 
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);
  bool fRunForMB = kFALSE;  // used to differentiate pp and AuAu
  if(doppAnalysis)  fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

  // fill arrays for towers that fired trigger
  FillTowerTriggersArr();

  // switches for Jet and Event Plane analysis
  Bool_t doJetAnalysis = kFALSE; // set false by default
  Bool_t doEPAnalysis = kFALSE;  // set false by default

  // if we have trigger: perform jet analysis
  if(fHaveEmcTrigger) { doJetAnalysis = kTRUE; }

  // if we have trigger && AuAu dataset: run event plane analysis
  if(fHaveEmcTrigger && (fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200 || fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB || fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200)) {
    doEPAnalysis = kTRUE;
  }
  // ======================== end of Triggers ============================= //

  // ================= JetMaker ================ //
  // get JetMaker
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get JetMaker collection of jets
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn;
  }

  // fill stats histo to do QA
  hStats->Fill(7);
  if(fJets->GetEntries() > 0) hStats->Fill(8);

  // check if bad/dead towers fired trigger and kill event if true
  //if(DidBadTowerFireTrigger()) { return kStOK; }

  // ================================================================================

  // ============== RhoMaker =============== //
  // get RhoMaker from event: old names "StRho_JetsBG", "OutRho", "StMaker#0"
  RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  if(!fRho) {
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;
  }

  // get rho/area value from rho object     fRho->ls("");
  //double value = GetRhoValue(fRhoMakerName);
  fRhoVal = fRho->GetVal();
  hRhovsCent->Fill(centbin*5.0, fRhoVal);
  hStats->Fill(9);
  if(fDebugLevel == kDebugRhoEstimate) cout<<"   fRhoVal = "<<fRhoVal<<"   Correction = "<<1.0*TMath::Pi()*fJetRad*fJetRad*fRhoVal<<endl;

  // =========== Leading and Subleading Jets ============= //
  // cache the leading + subleading jets within acceptance
  if(fCorrJetPt) {
    fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
  } else {
    fLeadingJet = GetLeadingJet(fJetMakerName);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
  }

  // fill leading and subleading pt spectra
  if(fLeadingJet)    hJetLeadingPt->Fill(fLeadingJet->Pt());
  if(fSubLeadingJet) hJetSubLeadingPt->Fill(fSubLeadingJet->Pt());   

  // get dijet imbalance Aj -  z = if(condition) then(?) <do this> else(:) <do this>  
  double fDiJetAj = (fLeadingJet && fSubLeadingJet) ? GetDiJetAj(fLeadingJet, fSubLeadingJet, fRho, fCorrJetPt) : -999.;
  if(fDiJetAj > 0)  hJetDiJetAj->Fill(fDiJetAj);

  // =========== Event Plane Angle ============= //
  // basic method - not used for anything..
  double rpAngle = GetReactionPlane();
  hEventPlane->Fill(rpAngle);

  // ============== EventPlaneMaker =============== //
  // get StEventPlaneMaker from event
  double angle2, angle2n, angle2p, angle1, angle4;
  StEventPlaneMaker *EventPlaneMaker0, *EventPlaneMaker1, *EventPlaneMaker2, *EventPlaneMaker3, *EventPlaneMaker4;
  double tpc2EP_bin0, tpc2EP_bin1, tpc2EP_bin2, tpc2EP_bin3, tpc2EP_bin4;
  if(doEPAnalysis) {
    if(!doTPCptassocBin) {
      EventPlaneMaker = static_cast<StEventPlaneMaker*>(GetMaker(fEventPlaneMakerName));
      const char *fEventPlaneMakerNameCh = fEventPlaneMakerName;
      if(!EventPlaneMaker) {
        LOG_WARN << Form(" No %s! Skip! ", fEventPlaneMakerNameCh) << endm;
        return kStWarn;
      }

      // set event plane
      double tpc2EP = (double)EventPlaneMaker->GetTPCEP();
      TPC_PSI2 = tpc2EP;

    } else { // pt-dependent bin mode
      const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
      EventPlaneMaker0 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 0)));
      EventPlaneMaker1 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 1)));
      EventPlaneMaker2 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 2)));
      EventPlaneMaker3 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 3)));
      EventPlaneMaker4 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 4)));

      // check for requested EventPlaneMaker pointer
      if((fTPCptAssocBin == 0) && (!EventPlaneMaker0)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 1) && (!EventPlaneMaker1)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 2) && (!EventPlaneMaker2)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 3) && (!EventPlaneMaker3)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 4) && (!EventPlaneMaker4)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}

      // get event plane angle for different pt bins
      // could also write this as:  tpc2EP_bin = (EventPlaneMaker) ? (double)EventPlaneMaker->GetTPCEP() : -999;
      if(EventPlaneMaker0) { tpc2EP_bin0 = (double)EventPlaneMaker0->GetTPCEP(); } else { tpc2EP_bin0 = -999; }
      if(EventPlaneMaker1) { tpc2EP_bin1 = (double)EventPlaneMaker1->GetTPCEP(); } else { tpc2EP_bin1 = -999; }
      if(EventPlaneMaker2) { tpc2EP_bin2 = (double)EventPlaneMaker2->GetTPCEP(); } else { tpc2EP_bin2 = -999; }
      if(EventPlaneMaker3) { tpc2EP_bin3 = (double)EventPlaneMaker3->GetTPCEP(); } else { tpc2EP_bin3 = -999; }
      if(EventPlaneMaker4) { tpc2EP_bin4 = (double)EventPlaneMaker4->GetTPCEP(); } else { tpc2EP_bin4 = -999; }

      // assign global event plane to selected pt-dependent bin
      if(fTPCptAssocBin == 0) TPC_PSI2 = tpc2EP_bin0;
      if(fTPCptAssocBin == 1) TPC_PSI2 = tpc2EP_bin1;
      if(fTPCptAssocBin == 2) TPC_PSI2 = tpc2EP_bin2;
      if(fTPCptAssocBin == 3) TPC_PSI2 = tpc2EP_bin3;
      if(fTPCptAssocBin == 4) TPC_PSI2 = tpc2EP_bin4;
      //cout<<"EP angles:  "<<tpc2EP_bin0<<"   "<<tpc2EP_bin1<<"   "<<tpc2EP_bin2<<"   "<<tpc2EP_bin3<<"   "<<tpc2EP_bin4<<endl;
    }

    // fill histogram with angle from GetReactionPlane() function
    angle1 = GetReactionPlane();
    hReactionPlaneFnc->Fill(angle1);

    // fill histogram with angle from GetEventPlane() function
    GetEventPlane(kFALSE, 2, fTPCEPmethod, 2.0, fTPCptAssocBin); // 2nd to last param not used (ptcut)
    angle2     = fEPTPC;
    angle2n    = fEPTPCn;
    angle2p    = fEPTPCp;
    hEventPlaneFncN2->Fill(angle2n);
    hEventPlaneFncP2->Fill(angle2p);
    hEventPlaneFnc2->Fill(angle2);

    // fill histogram with angle from StEventPlaneMaker class
    angle4 = TPC_PSI2;
    hEventPlaneClass->Fill(angle4);
  }

  // FIXME FIXME - March 3, 2019 - this is outdated at this point, need to recode below chunk if doing ptbin approach
  // switch to require specific trigger (for Event Plane corrections + Resolution)
  if(doEPAnalysis) {
    // compare BBC, ZDC, TPC event planes
    // only truely relevant for STEP3 - when both recentering and shifting corrections are read in
    hTPCvsBBCep->Fill(BBC_PSI2, TPC_PSI2);
    hTPCvsZDCep->Fill(ZDC_PSI2, TPC_PSI2);
    hBBCvsZDCep->Fill(ZDC_PSI2, BBC_PSI2);
 
    // calculate / fill event plane resolution histograms
    if(doEventPlaneRes){
      // corrected event planes used for Resolution
      // MAKE sure for meaningful results to be in mode: STEP3
      CalculateEventPlaneResolution(BBC_PSI2, ZDC_PSI2, TPC_PSI2, TPCA_PSI2, TPCB_PSI2, BBC_PSI1, ZDC_PSI1);
    }
  } // have Emc HT trigger - process event plane calculation/ corrections / resolutions
  // ===================================================================================

  // ========================== Jet Shape Analysis ===================================== //
  int jsret = -99;
  if(doJetShapeAnalysis) {
    hStats->Fill(11);

    // declare pool pointer
    StEventPool *pool = 0x0;

    // require event mixing
    if(fDoEventMixing > 0) {
      // convert back to integer bins for mixed event pool - 10% bins (0, 7), 5% bins (0, 15)
      Int_t mixcentbin = TMath::Floor(fCentralityScaled / fCentBinSizeJS);
      //cout<<"fCentralityScaled: "<<fCentralityScaled<<"  fCentBinSizeJS: "<<fCentBinSizeJS<<"  mixcentbin: "<<mixcentbin<<"  zVtx: "<<zVtx<<endl;
  
      // for an angle (0-pi)    
      const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
      StEventPlaneMaker *EPMaker = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 4)));
      double psi2 = (EPMaker) ? (double)EPMaker->GetTPCEP() : -999;

      double psibin = -1.;
      if(psi2 >=       0.0 && psi2 <  1.*pi/6.0) psibin = 0.5; // psibin = 1;
      if(psi2 >= 1.*pi/6.0 && psi2 <  2.*pi/6.0) psibin = 1.5; // psibin = 2;
      if(psi2 >= 2.*pi/6.0 && psi2 <  4.*pi/6.0) psibin = 2.5; // psibin = 3;
      if(psi2 >= 4.*pi/6.0 && psi2 <  6.*pi/6.0) psibin = 1.5; // psibin = 2;
      if(psi2 >= 5.*pi/6.0 && psi2 <= 6.*pi/6.0) psibin = 0.5; // psibin = 1;
      //cout<<"psi2 = "<<psi2<<"  psibin = "<<psibin<<endl;
    
      // initialize event pools - different cases for each dataset
      if(fDoUseMultBins) {
        if(!doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx);
        //if( doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx, psibin);
        if( doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx, psi2);
      } else {
        if(!doUseEPBins) pool = fPoolMgr->GetEventPool(mixcentbin, zVtx);
        if( doUseEPBins) pool = fPoolMgr->GetEventPool(mixcentbin, zVtx, psi2); // psibin
      }

      // check if pool exists
      if(!pool) {
        Form("No pool found for centrality = %f, zVtx = %f", (double)mixcentbin, zVtx);
        return kTRUE;
      }

      //cout<<"nPsiIndex: "<<pool->PsiBinIndex()<<"   zvertexbinIndex: "<<pool->ZvtxBinIndex()<<"   multbinindex: "<<pool->MultBinIndex()<<endl;
    }

    // initialize some variables
    double ljpttemp = 0., sljpttemp = 0.;
    bool isBackToBack = kFALSE;

    // check for back-to-back jets
    if(fLeadingJet && fSubLeadingJet) { 
      double BackToBackPhi = fLeadingJet->Phi() - fSubLeadingJet->Phi() - pi;
      if(BackToBackPhi < 0.0) BackToBackPhi += 2.0*pi;
      if(BackToBackPhi < 0.4) isBackToBack = kTRUE;

      // get lj and sublj pt
      if(fCorrJetPt) {
        ljpttemp = fLeadingJet->Pt() - fLeadingJet->Area()*fRhoVal;
        sljpttemp = fSubLeadingJet->Pt() - fLeadingJet->Area()*fRhoVal;
      } else {
        ljpttemp = fLeadingJet->Pt();
        sljpttemp = fSubLeadingJet->Pt(); 
      }

    } // lj + sublj

    // Aj selection for Jet Shape Analysis: LJ pt > 20, SubLJ pt > 10 defaults
    bool doAjSelection = (isBackToBack && ljpttemp > fLeadJetPtMin && sljpttemp > fSubLeadJetPtMin) ? kTRUE : kFALSE;
    if(doAjSelection) hJetLeadingPtAj->Fill(ljpttemp);      // leading jet w/ Aj cut
    if(doAjSelection) hJetSubLeadingPtAj->Fill(sljpttemp);  // sub-leading jet w/ Aj cut

    // Triggered events and leading/subleading jets - do Jet Shape Analysis
    // check for back to back jets: must have leading + subleading jet, subleading jet must be > 10 GeV, subleading jet must be within 0.4 of pi opposite of leading jet
    if(doRequireAjSelection) {
      for(int ptbin=0; ptbin<9; ptbin++) {
        if(doAjSelection && fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet)       cout<<endl; //
        if(doAjSelection && fHaveEmcTrigger && fJetAnalysisJetType == kSubLeadingJets && fSubLeadingJet) cout<<endl; //
      }
    } else { // don't require back-to-back jets meeting Aj criteria
      for(int ptbin=0; ptbin<9; ptbin++) {
        if(fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet)       cout<<endl; //
        if(fHaveEmcTrigger && fJetAnalysisJetType == kSubLeadingJets && fSubLeadingJet) cout<<endl; //
      }
    }

    // jet shape - case for: inclusive jets
    if(fHaveEmcTrigger && fJetAnalysisJetType == kInclusiveJets) {

      // loop over Jets in the event
      for(int ijet = 0; ijet < fJets->GetEntries(); ijet++) {  // JET LOOP
        // get pointer to jet
        StJet *jet = static_cast<StJet*>(fJets->At(ijet));
        if(!jet) continue;

        // loop over pt associated bins
        for(int ptbin = 0; ptbin < 9; ptbin++) {
          if(doRequireAjSelection && doAjSelection) { cout<<endl;
          } else { cout<<endl;  }

        }

      } // loop over jets
    }   // inclusive jet case

    // use only tracks from MB events
    //if(fDoEventMixing > 0 && fRunForMB && (!fHaveEmcTrigger)) { // kMB5 or kMB30 - AuAu, kMB - pp (excluding HT)
    if(fDoEventMixing > 0 && fRunForMB) { // kMB5 or kMB30 - AuAu, kMB - pp (don't exclude HT)
      // update pool: create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      pool->UpdatePool(CloneAndReduceTrackList());

      // fill QA histo's
      hMBvsMult->Fill(kEventActivity);    // MB5 || MB30 AuAu and MB pp
      if(fHaveMB5event)  hMB5vsMult->Fill(refCorr2);   // MB5
      if(fHaveMB30event) hMB30vsMult->Fill(refCorr2);  // MB30
      hMixEvtStatZVtx->Fill(zVtx);
      hMixEvtStatCent->Fill(centBinToUse);
      hMixEvtStatZvsCent->Fill(centBinToUse, zVtx);
    } // MB 

    // only return if doing jet shape analysis
    FillEventTriggerQA(fHistEventSelectionQAafterCuts);

    // only return if doing jet shape analysis
    return kStOK;  
  }

  // run Track QA and fill histograms
  if((doWriteTrackQAHist) && (doJetAnalysis)) TrackQA();

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();
  if(fDebugLevel == kDebugGeneralEvt) {
    cout<<"njets = "<<njets<<"  ntracks = "<<ntracks<<"  nglobaltracks = "<<nglobaltracks<<"  refCorr2 = "<<refCorr2<<"  grefMult = "<<grefMult<<"  centbin = "<<centbin<<endl;
  }

  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables
  Int_t ijethi = -1;
  Double_t highestjetpt = 0.0;
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    // Run - Trigger Selection to process jets from
    if(!doJetAnalysis) continue;

    // get pointer to jets
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    // get some jet parameters
    double jetarea = jet->Area();
    double jetpt = jet->Pt();
    double corrjetpt = jet->Pt() - jetarea*fRhoVal;
    double jetE = jet->E();
    double jetEta = jet->Eta();
    double jetPhi = jet->Phi();    
    double jetNEF = jet->NEF();
    //double dEP = RelativeEPJET(jet->Phi(), rpAngle);         // difference between jet and EP
    if(!doppAnalysis) {
      double dEPmethod1      = RelativeEPJET(jetPhi, angle1);       // from Reaction Plane function
      double dEPmethod2p     = RelativeEPJET(jetPhi, angle2p);      // from Event Plane function: sub event 1 - method 2
      double dEPmethod2n     = RelativeEPJET(jetPhi, angle2n);      // from Event Plane function: sub event 2 - method 2
      double dEPmethod2      = RelativeEPJET(jetPhi, angle2);       // from Event Plane function - method 2
      double dEPmethod4      = RelativeEPJET(jetPhi, angle4);       // from Event Plane class
      hdEPReactionPlaneFnc->Fill(dEPmethod1);
      hdEPEventPlaneFncN2->Fill(dEPmethod2n);
      hdEPEventPlaneFncP2->Fill(dEPmethod2p);
      hdEPEventPlaneFnc2->Fill(dEPmethod2);
      hdEPEventPlaneClass->Fill(dEPmethod4);
    }

    // get pt dependent event plane calculation
    // FIXME, double check, might want to code this nicer
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, TPC_PSI2) : 0; // CORRECTED event plane angle - STEP3
    if(doTPCptassocBin && !doppAnalysis) {
      // z = if(condition) then(?) <do this> else(:) <do this>  
      double dEP0 = (EventPlaneMaker0) ? RelativeEPJET(jetPhi, tpc2EP_bin0) : -999;
      double dEP1 = (EventPlaneMaker1) ? RelativeEPJET(jetPhi, tpc2EP_bin1) : -999;
      double dEP2 = (EventPlaneMaker2) ? RelativeEPJET(jetPhi, tpc2EP_bin2) : -999;
      double dEP3 = (EventPlaneMaker3) ? RelativeEPJET(jetPhi, tpc2EP_bin3) : -999;
      double dEP4 = (EventPlaneMaker4) ? RelativeEPJET(jetPhi, tpc2EP_bin4) : -999;
      if(fTPCptAssocBin == 0) dEP = dEP0;
      if(fTPCptAssocBin == 1) dEP = dEP1;
      if(fTPCptAssocBin == 2) dEP = dEP2;
      if(fTPCptAssocBin == 3) dEP = dEP3;
      if(fTPCptAssocBin == 4) dEP = dEP4;

      if(fDebugLevel == kDebugJetvsEPtype) cout<<"jetPhi = "<<jetPhi<<"  RP = "<<angle1<<"  bin0 = "<<tpc2EP_bin0<<"  bin1 = "<<tpc2EP_bin1<<"  bin2 = "<<tpc2EP_bin2<<"  bin3 = "<<tpc2EP_bin3<<"  bin4 = "<<tpc2EP_bin4<<endl;
      if(fDebugLevel == kDebugJetvsEPtype) cout<<"dRP = "<<dEP1<<"  dEP0: "<<dEP0<<"  dEP1: "<<dEP1<<"  dEP2: "<<dEP2<<"  dEP3: "<<dEP3<<"  dEP4: "<<dEP4<<endl;
    }

    // this is a double check - should not happen, but if it does -> kill the job so it can be fixed - most likely in the runPico macro
    if(dEP < -900) return kStFatal;

    // some threshold cuts
    if(fCorrJetPt) {  // background subtracted jet pt
      if(corrjetpt < fMinPtJet) continue;
    } else { if(jetpt < fMinPtJet) continue; }
    if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) continue;

    // TODO check that jet contains a tower that fired the trigger
    if(!DidTowerConstituentFireTrigger(jet)) { continue; }
    // =======================================================================================================================================

    // loop over constituent tracks
    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);      

      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
      if(!trk){ continue; }

      // track momentum vector
      TVector3 mTrkMom;
      if(doUsePrimTracks) {
        if(!(trk->isPrimary())) return kFALSE; // check if primary
        // get primary track vector
        mTrkMom = trk->pMom();
      } else {
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield);
      }

      // track variables
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      double pt = mTrkMom.Perp();
      double px = mTrkMom.x();
      double py = mTrkMom.y();
      double pz = mTrkMom.z();
      double jetZ = jet->GetZ(px, py, pz);

      // shift angle (0, 2*pi) 
      if(phi < 0.0)    phi += 2.0*pi;
      if(phi > 2.0*pi) phi -= 2.0*pi;

      // fill jet track constituent histograms
      hJetTracksPt->Fill(pt);
      hJetTracksPhi->Fill(phi);
      hJetTracksEta->Fill(eta);
      hJetTracksZ->Fill(jetZ);
    } // constituent track loop

    // fill some jet histograms
    hJetPt->Fill(jetpt);
    hJetCorrPt->Fill(corrjetpt);
    hJetE->Fill(jetE);
    hJetEta->Fill(jetEta);
    hJetPhi->Fill(jetPhi);
    hJetNEF->Fill(jetNEF);
    hJetArea->Fill(jetarea);
    hJetPtvsArea->Fill(jetpt, jetarea);
    hJetEventEP->Fill(TPC_PSI2);
    hJetPhivsEP->Fill(jetPhi, TPC_PSI2);

    // fill some jet QA plots for each orientation
    if(dEP >= 0.0*pi/6.0 && dEP < 1.0*pi/6.0) {
      hJetPtIn->Fill(jetpt);
      hJetPhiIn->Fill(jetPhi);
      hJetEtaIn->Fill(jetEta);
      hJetEventEPIn->Fill(TPC_PSI2);
      hJetPhivsEPIn->Fill(jetPhi, TPC_PSI2);
    } else if(dEP >= 1.0*pi/6.0 && dEP < 2.0*pi/6.0) {
      hJetPtMid->Fill(jetpt);
      hJetPhiMid->Fill(jetPhi);
      hJetEtaMid->Fill(jetEta);
      hJetEventEPMid->Fill(TPC_PSI2);
      hJetPhivsEPMid->Fill(jetPhi, TPC_PSI2);
    } else if(dEP >= 2.0*pi/6.0 && dEP <= 3.0*pi/6.0) {
      hJetPtOut->Fill(jetpt);
      hJetPhiOut->Fill(jetPhi);
      hJetEtaOut->Fill(jetEta);
      hJetEventEPOut->Fill(TPC_PSI2);
      hJetPhivsEPOut->Fill(jetPhi, TPC_PSI2);
    }

    // get nTracks and maxTrackPt
    double maxtrackpt = jet->GetMaxTrackPt();
    int NtrackConstit = jet->GetNumberOfTracks();
    if(doComments) cout<<"Jet# = "<<ijet<<"  JetPt = "<<jetpt<<"  JetE = "<<jetE<<"   MaxtrackPt = "<<maxtrackpt<<"  NtrackConstit = "<<NtrackConstit<<endl;

    // get highest Pt jet in event (leading jet)
    if(highestjetpt < jetpt){
      ijethi = ijet;
      highestjetpt = jetpt;
    }

    // the below track and cluster cut is already done
    //if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fTowerBias)){
    // set up and fill Jet-Hadron trigger jets THnSparse
    // ====================================================================================
    double jetptselected;
    if(fCorrJetPt) { jetptselected = corrjetpt; 
    } else { jetptselected = jetpt; }

    // fill jet sparse for trigger normalization
    Double_t CorrEntries[4] = {centBinToUse, jetptselected, dEP, zVtx};
    if(fReduceStatsCent > 0) {
      if(cbin == fReduceStatsCent) fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    } else fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    // ====================================================================================
    //} // check on max track and tower pt/E

    // track loop inside jet loop - loop over ALL tracks in PicoDst
    for(int itrack = 0; itrack < ntracks; itrack++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

      // get momentum vector of track - global or primary track
      TVector3 mTrkMom;
      if(doUsePrimTracks) {
        // get primary track vector
        mTrkMom = trk->pMom();
      } else {
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield);
      }

      // track variables
      double pt = mTrkMom.Perp();
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      short charge = trk->charge();

      // get jet - track relations
      //deta = eta - jetEta;               // eta betweeen hadron and jet
      double deta = jetEta - eta;               // eta betweeen jet and hadron
      double dphijh = RelativePhi(jetPhi, phi); // angle between jet and hadron

      // === June 1, 2018 - these cuts should not be here!
/*
      // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
      // when doing event plane calculation via pt assoc bin
      if(doTPCptassocBin) {
        if(fTPCptAssocBin == 0) { if((pt < 0.20) || (pt >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 1) { if((pt < 0.50) || (pt >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 2) { if((pt < 1.00) || (pt >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 3) { if((pt < 1.50) || (pt >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 4) { if((pt < 2.00) || (pt >= 30.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
        if(fTPCptAssocBin == 5) { if((pt < 2.00) || (pt >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 6) { if((pt < 3.00) || (pt >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 7) { if((pt < 4.00) || (pt >= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
      }
*/

      // fill jet sparse 
      double triggerEntries[8] = {centBinToUse, jetptselected, pt, deta, dphijh, dEP, zVtx, (double)charge};
      double trefficiency = 1.0;
      //if(fDoEventMixing) {
        if(fReduceStatsCent > 0) {
          if(cbin == fReduceStatsCent) fhnJH->Fill(triggerEntries, 1.0/trefficiency);    // fill Sparse Histo with trigger entries
        } else fhnJH->Fill(triggerEntries, 1.0/trefficiency);
      //}

      fHistJetHEtaPhi->Fill(deta,dphijh);                          // fill jet-hadron  eta--phi distributio
    // =====================================================================================
    } // track loop

  } // jet loop

// ***************************************************************************************************************
// ******************************** Event MIXING *****************************************************************
// ***************************************************************************************************************
  if(fDebugLevel == kDebugMixedEvents) cout<<"StMyAnMaker event# = "<<EventCounter()<<"  Centbin = "<<centbin<<"  zVtx = "<<zVtx<<endl;

  // Prepare to do event mixing
  if(fDoEventMixing > 0){
    // event mixing

    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    // mix jets from triggered events with tracks from MB events
    // get the trigger bit, need to change trigger bits between different runs

    //Double_t mycentbin = (Double_t)centbin + 0.001;
    // mixed event centbin - FIXME
    Int_t mixcentbin;
    if(fCentBinSize==10)       { mixcentbin = centbin10;    // 10% bins (0,  7)
    } else if(fCentBinSize==5) { mixcentbin = centbin; }    //  5% bins (0, 15)
    //cout<<"mixcentbin = "<<mixcentbin<<"  centbin = "<<centbin<<" centbin10 = "<<centbin10<<"  zvtx = "<<zVtx<<endl;

    // initialize event pools
    StEventPool *pool = fPoolMgr->GetEventPool(mixcentbin, zVtx); //FIXME? 
    if(!pool) {
      Form("No pool found for centrality = %i, zVtx = %f", mixcentbin, zVtx); // FIXME if cent changes to double
      return kTRUE; // FIXME
    }

    if(fDebugLevel == kDebugMixedEvents) cout<<"NtracksInPool = "<<pool->NTracksInPool()<<"  CurrentNEvents = "<<pool->GetCurrentNEvents()<<endl;

    // initialize background tracks array
    TObjArray *bgTracks;

  // do event mixing when Signal Jet is part of event with a HT1 or HT2 or HT3 trigger firing
  if(doJetAnalysis) { // trigger type requested was fired for this event - do mixing
    if(pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {

      // loop over jets (passing cuts - set by jet maker)
      for(int ijet = 0; ijet < njets; ijet++) {
        // leading jet - why was this a double?
        Int_t leadjet = 0;
        if(ijet == ijethi) leadjet = 1; //FIXME for leading jet

        // get jet pointer
        StJet *jet = static_cast<StJet*>(fJets->At(ijet));
        if(!jet) continue;

        // get some get parameters of jets for mixing
        double Mixjetarea = jet->Area();
        double Mixjetpt = jet->Pt();
        double Mixcorrjetpt = Mixjetpt - Mixjetarea*fRhoVal;
        double MixjetEta = jet->Eta();
        double MixjetPhi = jet->Phi();
        //double dMixEP = RelativeEPJET(jet->Phi(), rpAngle);         // difference between jet and EP
        double dMixEP = (!doppAnalysis) ? RelativeEPJET(jet->Phi(), TPC_PSI2) : 0; // CORRECTED event plane angle - STEP3
        if(doTPCptassocBin && !doppAnalysis) {
          // z = if(condition) then(?) <do this> else(:) <do this>  
          double dMixEP0 = (EventPlaneMaker0) ? RelativeEPJET(MixjetPhi, tpc2EP_bin0) : -999;
          double dMixEP1 = (EventPlaneMaker1) ? RelativeEPJET(MixjetPhi, tpc2EP_bin1) : -999;
          double dMixEP2 = (EventPlaneMaker2) ? RelativeEPJET(MixjetPhi, tpc2EP_bin2) : -999;
          double dMixEP3 = (EventPlaneMaker3) ? RelativeEPJET(MixjetPhi, tpc2EP_bin3) : -999;
          double dMixEP4 = (EventPlaneMaker4) ? RelativeEPJET(MixjetPhi, tpc2EP_bin4) : -999;
          if(fTPCptAssocBin == 0) dMixEP = dMixEP0;
          if(fTPCptAssocBin == 1) dMixEP = dMixEP1;
          if(fTPCptAssocBin == 2) dMixEP = dMixEP2;
          if(fTPCptAssocBin == 3) dMixEP = dMixEP3;
          if(fTPCptAssocBin == 4) dMixEP = dMixEP4;
        }

        // this is a double check - should not happen, but if it does -> kill the job so it can be fixed - most likely in the runPico macro
        if(dMixEP < -900) return kStFatal;

        // some threshold cuts - do mixing only if we have a jet meeting our pt threshold and bias
        if(fCorrJetPt) {
          if(Mixcorrjetpt < fMinPtJet) continue;
        } else { if(Mixjetpt < fMinPtJet) continue; }
        //if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) continue;
   	//if (!AcceptJet(jet)) continue;  // acceptance cuts done to jet in JetMaker

        // get number of current events in pool
        int nMix = pool->GetCurrentNEvents();

        // check that jet contains a tower that fired the trigger
        if( !DidTowerConstituentFireTrigger(jet) ) continue;

        // Fill for biased jet triggers only
        if((jet->GetMaxTrackPt() > fTrackBias) || (jet->GetMaxTowerEt() > fTowerBias)) {
          // Fill mixed-event histos here: loop over nMix events
          for(int jMix = 0; jMix < nMix; jMix++) {
            // get jMix'th event
            bgTracks = pool->GetEvent(jMix);

            //TObjArray *bgTracks = pool->GetEvent(jMix);
            const Int_t Nbgtrks = bgTracks->GetEntries();

            // loop over background (mixed event) tracks
            for(int ibg = 0; ibg < Nbgtrks; ibg++) {
              // get Femto track pointer
              StFemtoTrack *trk = static_cast<StFemtoTrack*>(bgTracks->At(ibg));
              if(!trk) continue;

              // mixed track variables
              double Mixphi = trk->Phi();
              double Mixeta = trk->Eta();
              double Mixpt = trk->Pt();
              short Mixcharge = trk->Charge();

              // shift angle (0, 2*pi) 
              if(Mixphi < 0.0)    Mixphi += 2.0*pi;
              if(Mixphi > 2.0*pi) Mixphi -= 2.0*pi;
              if(fDebugLevel == kDebugMixedEvents) cout<<"itrack = "<<ibg<<"  phi = "<<Mixphi<<"  eta = "<<Mixeta<<"  pt = "<<Mixpt<<"  q = "<<Mixcharge<<endl;

              // get jet - track relations
              //double deta = eta - jetEta;               // eta betweeen hadron and jet
              double dMixeta = MixjetEta - Mixeta;               // eta betweeen jet and hadron
              double dMixphijh = RelativePhi(MixjetPhi, Mixphi); // angle between jet and hadron

              // print tracks outside of acceptance somehow
              if(fDebugLevel == kDebugMixedEvents) if(TMath::Abs(dMixeta > 1.6)) cout<<"DELTA ETA out of bounds... deta = "<<dMixeta<<"   iTrack = "<<ibg<<"  jetEta = "<<MixjetEta<<"  trk eta = "<<Mixeta<<endl;

              // calculate single particle tracking efficiency of mixed events for correlations (-999)
              double mixefficiency = 1.0;
              // FIXME mixefficiency = EffCorrection(part->Eta(), part->Pt(), fDoEffCorr);                           

              // select which jet pt to use for filling
              double Mixjetptselected;
              if(fCorrJetPt) { Mixjetptselected = Mixcorrjetpt;
              } else { Mixjetptselected = Mixjetpt; }

              // === June 1, 2018 - these cuts should not be here!
/*
              // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
              // when doing event plane calculation via pt assoc bin
              if(doTPCptassocBin) {
                if(fTPCptAssocBin == 0) { if((Mixpt < 0.20) || (Mixpt >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 1) { if((Mixpt < 0.50) || (Mixpt >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 2) { if((Mixpt < 1.00) || (Mixpt >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 3) { if((Mixpt < 1.50) || (Mixpt >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 4) { if((Mixpt < 2.00) || (Mixpt >= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
                if(fTPCptAssocBin == 5) { if((Mixpt < 2.00) || (Mixpt >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 6) { if((Mixpt < 3.00) || (Mixpt >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 7) { if((Mixpt < 4.00) || (Mixpt >= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
              }
*/

              // create / fill mixed event sparse  (centbin*5.0)
              double triggerEntries[8] = {centBinToUse, Mixjetptselected, Mixpt, dMixeta, dMixphijh, dMixEP, zVtx, (double)Mixcharge};

              if(fReduceStatsCent > 0) {
                if(cbin == fReduceStatsCent) fhnMixedEvents->Fill(triggerEntries, 1./(nMix*mixefficiency));
              } else fhnMixedEvents->Fill(triggerEntries, 1./(nMix*mixefficiency));   // fill Sparse histo of mixed events

            } // end of background track loop
          } // end of filling mixed-event histo's:  jth mix event loop
        } // end of check for biased jet triggers
      } // end of jet loop
    } // end of check for pool being ready
  } // end EMC triggered loop

    // use only tracks from MB (and Semi-Central) events
    ///if(fMixingEventType) { //kMB) {
    //if(fRunForMB && (!fHaveEmcTrigger)) { // kMB or kMB30 - AuAu, kMB - pp, TODO probably want to use this line in future, may not matter
    //if(fRunForMB) { // kMB or kMB30 (don't exclude HT)
    if(fHaveMBevent) { // kMB - FIXME, this is how it was set up for correlations
      // update pool - create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      pool->UpdatePool(CloneAndReduceTrackList());

      // fill QA histo's
      hMixEvtStatZVtx->Fill(zVtx);
      hMixEvtStatCent->Fill(centBinToUse);
      hMixEvtStatZvsCent->Fill(centBinToUse, zVtx);
    } // MB event - updating pool

  } // end of event mixing

  // event counter at end of maker
  mInputEventCounter++;

  // fill Event Trigger QA
  //FillEventTriggerQA(fHistEventSelectionQAafterCuts);
  //StMemStat::PrintMem("MyAnalysisMaker at end of make");

  return kStOK;
}
//
//______________________________________________________________________
THnSparse* StJetHadronCorrelations::NewTHnSparseF(const char* name, UInt_t entries)
{
   // generate new THnSparseF, axes are defined in GetDimParams()
   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){
         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }

      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF
//
//________________________________________________________________________
void StJetHadronCorrelations::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

   case 0:
      label = "centrality 5% bin";
      if(fCentBinSize==10) {
        nbins = 10;
      } else {
        nbins = 20; //16;
      } 
      xmin = 0.;
      xmax = 100.; //16.;     
      break;

   case 1:
      if(fCorrJetPt) { // correct jet pt
        label = "Jet Corrected p_{T}";
        nbins = 30;
        xmin = -50.;
        xmax = 100.;
      } else { // don't correct jet pt
        label = "Jet p_{T}";
        nbins = 20;
        xmin = 0.;
        xmax = 100.;
      }
      break;

   case 2:
      label = "Track p_{T}";
      nbins = 80; 
      xmin = 0.;
      xmax = 20.;
      break;

   case 3:
      label = "Relative Eta";
      nbins = 72; // 48
      xmin = -1.8;
      xmax = 1.8;
      break;

   case 4: 
      label = "Relative Phi";
      nbins = 72;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;

   case 5:
      label = "Relative angle of Jet and Reaction Plane";
      nbins = 3; // (12) 72
      xmin = 0;
      xmax = 0.5*pi;
      break;

   case 6:
      label = "z-vertex";
      nbins = 20; // 10
      xmin = -40; //-10
      xmax =  40; //+10
      break;

   case 7:
      label = "track charge";
      nbins = 3;
      xmin = -1.5;
      xmax = 1.5;
      break;

   case 8:
      label = "leading jet";
      nbins = 3;
      xmin = -0.5;
      xmax = 2.5;
      break;

   case 9: // need to update
      label = "leading track";
      nbins = 10;
      xmin = 0;
      xmax = 50;
      break; 

   } // end of switch
} // end of getting dim-params
//
//______________________________________________________________________
THnSparse* StJetHadronCorrelations::NewTHnSparseFCorr(const char* name, UInt_t entries) {
  // generate new THnSparseD, axes are defined in GetDimParamsD()
  Int_t count = 0;
  UInt_t tmp = entries;
  while(tmp!=0){
    count++;
    tmp = tmp &~ -tmp;  // clear lowest bit
  }

  TString hnTitle(name);
  const Int_t dim = count;
  Int_t nbins[dim];
  Double_t xmin[dim];
  Double_t xmax[dim];

  Int_t i=0;
  Int_t c=0;
  while(c<dim && i<32){
    if(entries&(1<<i)){
      TString label("");
      GetDimParamsCorr(i, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label.Data());
      c++;
    }

    i++;
  }
  hnTitle += ";";

  return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF
//
//______________________________________________________________________________________________
void StJetHadronCorrelations::GetDimParamsCorr(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   //stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

    case 0:
      label = "centrality 5% bin";
      if(fCentBinSize==10) {
        nbins = 10;
      } else {
        nbins = 20; //16;
      }
      xmin = 0.;
      xmax = 100.; //16.;
      break;

    case 1:
      if(fCorrJetPt) { // correct jet pt
        label = "Jet Corrected p_{T}";
        nbins = 30;
        xmin = -50.;
        xmax = 100.;
      } else { // don't correct jet pt
        label = "Jet p_{T}";
        nbins = 20;
        xmin = 0.;
        xmax = 100.;
      }
      break;

    case 2:
      label = "Relative angle: Jet and Reaction Plane";
      nbins = 3; // (12) 48
      xmin = 0.;
      xmax = 0.5*pi;
      break;

    case 3:
      label = "Z-vertex";
      nbins = 20;
      xmin = -40.;
      xmax = 40.;
      break;

    case 4: // may delete this case
      label = "Jet p_{T} corrected with Rho";
      nbins = 50; // 250
      xmin = -50.;
      xmax = 200.;  
      break;

   }// end of switch
} // end of Correction (ME) sparse
//
// From CF event mixing code PhiCorrelations
// clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
//____________________________________________________________________________
TClonesArray* StJetHadronCorrelations::CloneAndReduceTrackList()
{
  // create array for Femto tracks
  TClonesArray *tracksClone = new TClonesArray("StFemtoTrack");
//  tracksClone->SetName("tracksClone");
//  tracksClone->SetOwner(kTRUE);

  // construct variables, get # of tracks
  int nMixTracks = mPicoDst->numberOfTracks();
  int iterTrk = 0;

  // loop over tracks
  for(int i = 0; i < nMixTracks; i++) { 
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      if(!(trk->isPrimary())) continue; // check if primary
      // get primary track vector
      mTrkMom = trk->pMom();
    } else {
      // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield);
    }

    // track variables - used with alt method below
    double pt = mTrkMom.Perp();
    //double phi = mTrkMom.Phi();
    //double eta = mTrkMom.PseudoRapidity();
    //short charge = trk->charge();

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    // this is TEMP, it will filter track by the pt bin used for analysis
    if(doTPCptassocBin && fDoFilterPtMixEvents) {
      if(fTPCptAssocBin == 0) { if((pt < 0.20) || (pt >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 1) { if((pt < 0.50) || (pt >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 2) { if((pt < 1.00) || (pt >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 3) { if((pt < 1.50) || (pt >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 4) { if((pt < 2.00) || (pt >= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(fTPCptAssocBin == 5) { if((pt < 2.00) || (pt >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 6) { if((pt < 3.00) || (pt >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 7) { if((pt < 4.00) || (pt >= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }

    // create StFemtoTracks out of accepted tracks - light-weight object for mixing
    //  StFemtoTrack *t = new StFemtoTrack(pt, eta, phi, charge);
    StFemtoTrack *t = new StFemtoTrack(trk, Bfield, mVertex, doUsePrimTracks);
    if(!t) continue;

    // add light-weight tracks passing cuts to TClonesArray
    ((*tracksClone)[iterTrk]) =  t;

    //delete t;
    ++iterTrk;
  } // end of looping through tracks

  return tracksClone;
}
//
//
//_________________________________________________________________________
TH1* StJetHadronCorrelations::FillEmcTriggersHist(TH1* h) {
  // set bin labels
  h->GetXaxis()->SetBinLabel(1, "HT0");
  h->GetXaxis()->SetBinLabel(2, "HT1");
  h->GetXaxis()->SetBinLabel(3, "HT2");
  h->GetXaxis()->SetBinLabel(4, "HT3");
  h->GetXaxis()->SetBinLabel(5, "JP0");
  h->GetXaxis()->SetBinLabel(6, "JP1");
  h->GetXaxis()->SetBinLabel(7, "JP2");
  h->GetXaxis()->SetBinLabel(10, "Any");
  h->LabelsOption("v");  // set x-axis labels vertically
  //h->LabelsDeflate("X");

  // number of Emcal Triggers
  for(int i = 0; i < 8; i++) { fEmcTriggerArr[i] = 0; }
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  //if(fDebugLevel == kDebugEmcTrigger) { cout<<"nEmcTrigger = "<<nEmcTrigger<<endl; }

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT0 = emcTrig->isHT0();
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();
    bool isJP0 = emcTrig->isJP0();
    bool isJP1 = emcTrig->isJP1();
    bool isJP2 = emcTrig->isJP2();

    // print some EMCal Trigger info: exclude JP2, JP1
    if(emcTrig->flag() != 112 && emcTrig->flag() != 96) { // FIXME - here for test May 30, 2019
      if(fDebugLevel == kDebugEmcTrigger) {
        cout<<"i = "<<i<<"  id = "<<emcTrig->id()<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc();
        cout<<"  isHT0: "<<isHT0<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3;
        cout<<"  isJP0: "<<isJP0<<"  isJP1: "<<isJP1<<"  isJP2: "<<isJP2<<endl;
      }
    }

    // fill for valid triggers
    if(isHT0) { h->Fill(1); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(isHT1) { h->Fill(2); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(isHT2) { h->Fill(3); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(isHT3) { h->Fill(4); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(isJP0) { h->Fill(5); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(isJP1) { h->Fill(6); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(isJP2) { h->Fill(7); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }

  // kAny trigger - filled once per event
  h->Fill(10); 

  return h;
}
//
// general function to get local reaction plane calculation
//________________________________________________________________________
Double_t StJetHadronCorrelations::GetReactionPlane() { 
  // initialize some variables and constants
  TVector2 mQ;
  double mQx = 0., mQy = 0.;
  int order = 2;
  double pi = 1.0*TMath::Pi();

  // leading jet check and removal
  Double_t excludeInEta = -999;
  if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from EP estimate
    if(fLeadingJet) excludeInEta = fLeadingJet->Eta();
  }

  // loop over tracks
  int nTrack = mPicoDst->numberOfTracks();
  for(int i = 0; i < nTrack; i++) {
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 5.0 GeV
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // check for leading jet removal
    if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit) )) continue;

    // configure track weight when performing Q-vector summation
    double trackweight;
    if(fTrackWeight == kNoWeight) {
      trackweight = 1.0;
    } else if(fTrackWeight == kPtLinearWeight) {
      trackweight = pt;
    } else if(fTrackWeight == kPtLinear2Const5Weight) {
      if(pt <= 2.0) trackweight = pt;
      if(pt >  2.0) trackweight = 2.0;
    } else {
      // nothing choosen, so don't use weight
      trackweight = 1.0;
    }

    // sum up q-vectors
    mQx += trackweight * cos(phi * order);
    mQy += trackweight * sin(phi * order);
  }

  // set q-vector components 
  mQ.Set(mQx, mQy);
  double psi = mQ.Phi() / order;

  return psi;
}
//
// Set the bin errors on histograms
// __________________________________________________________________________________
void StJetHadronCorrelations::SetSumw2() {
  // set sum weights
  hdEPReactionPlaneFnc->Sumw2();
  hdEPEventPlaneFncN2->Sumw2();
  hdEPEventPlaneFncP2->Sumw2();
  hdEPEventPlaneFnc2->Sumw2();
  hdEPEventPlaneClass->Sumw2();
  hReactionPlaneFnc->Sumw2();
  hEventPlaneFncN2->Sumw2();
  hEventPlaneFncP2->Sumw2();
  hEventPlaneFnc2->Sumw2();
  hEventPlaneClass->Sumw2();

  hEventPlane->Sumw2();
  fHistEPTPCn->Sumw2();
  fHistEPTPCp->Sumw2();
  fHistEPBBC->Sumw2();
  fHistEPZDC->Sumw2();
  //hEventZVertex->Sumw2();
  //hCentrality->Sumw2();
  //hMultiplicity->Sumw2();
  //hStats->Sumw2();
  //hRhovsCent->Sumw2();
  for(int i=0; i<5; i++) { hdEPtrk[i]->Sumw2(); }
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i]->Sumw2();
    hTrackEta[i]->Sumw2();
    hTrackPt[i]->Sumw2();
  }
  hTrackEtavsPhi->Sumw2();

  hJetPt->Sumw2();
  hJetCorrPt->Sumw2();
  hJetLeadingPt->Sumw2();
  hJetSubLeadingPt->Sumw2();
  hJetLeadingPtAj->Sumw2();
  hJetSubLeadingPtAj->Sumw2();
  hJetDiJetAj->Sumw2();
  hJetE->Sumw2();
  hJetEta->Sumw2();
  hJetPhi->Sumw2();
  hJetNEF->Sumw2();
  hJetArea->Sumw2();
  hJetTracksPt->Sumw2();
  hJetTracksPhi->Sumw2();
  hJetTracksEta->Sumw2();
  hJetTracksZ->Sumw2();
  hJetPtvsArea->Sumw2();
  hJetEventEP->Sumw2();
  hJetPhivsEP->Sumw2();

  hJetPtIn->Sumw2();
  hJetPhiIn->Sumw2();
  hJetEtaIn->Sumw2();
  hJetEventEPIn->Sumw2();
  hJetPhivsEPIn->Sumw2();
  hJetPtMid->Sumw2();
  hJetPhiMid->Sumw2();
  hJetEtaMid->Sumw2();
  hJetEventEPMid->Sumw2();
  hJetPhivsEPMid->Sumw2();
  hJetPtOut->Sumw2();
  hJetPhiOut->Sumw2();
  hJetEtaOut->Sumw2();
  hJetEventEPOut->Sumw2();
  hJetPhivsEPOut->Sumw2();

  fHistJetHEtaPhi->Sumw2();
  //fHistEventSelectionQA->Sumw2();
  //fHistEventSelectionQAafterCuts->Sumw2();
  //hTriggerIds->Sumw2();
  //hEmcTriggers->Sumw2();
  //hBadTowerFiredTrigger->Sumw2()
  hMixEvtStatZVtx->Sumw2();
  hMixEvtStatCent->Sumw2();
  hMixEvtStatZvsCent->Sumw2();
  hTriggerEvtStatZVtx->Sumw2();
  hTriggerEvtStatCent->Sumw2();
  hTriggerEvtStatZvsCent->Sumw2();
  hMBvsMult->Sumw2();
  hMB5vsMult->Sumw2();
  hMB30vsMult->Sumw2();
  hHTvsMult->Sumw2();
  hNMixEvents->Sumw2();

  hTPCvsBBCep->Sumw2();
  hTPCvsZDCep->Sumw2();
  hBBCvsZDCep->Sumw2();

/*
  // don't set Sumw2 for TProfile plots
  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        fProfJetV2[k][j][i]->Sumw2();
      }
    }
  }
*/

/*
  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      fProfV2Resolution[i]->Sumw2();
      fProfV3Resolution[i]->Sumw2();
      fProfV4Resolution[i]->Sumw2();
      fProfV5Resolution[i]->Sumw2();
    }
  }
*/

  fhnJH->Sumw2();
  fhnMixedEvents->Sumw2();
  fhnCorr->Sumw2();
}
//
//
// ________________________________________________________________________________________
void StJetHadronCorrelations::GetEventPlane(Bool_t flattenEP, Int_t n, Int_t method, Double_t ptcut, Int_t ptbin)
{ 
  // local variables
  TVector2 mQtpcn, mQtpcp, mQtpc;
  double mQtpcnx = 0., mQtpcny = 0., mQtpcpx = 0., mQtpcpy = 0., mQtpcX = 0., mQtpcY = 0.;
  int order = n;
  double pi = 1.0*TMath::Pi();
  int ntracksNEG = 0, ntracksPOS = 0;

  // leading jet check and removal
  double excludeInEta = -999., excludeInPhi = -999.;
  double excludeInEtaSub = -999., excludeInPhiSub = -999.;
  if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from EP estimate
    // check for leading jet
    if(fLeadingJet) {
      excludeInEta = fLeadingJet->Eta();
      excludeInPhi = fLeadingJet->Phi();
      //cout<<"leading: pt = "<<fLeadingJet->Pt()<<"  eta = "<<fLeadingJet->Eta()<<"  phi = "<<fLeadingJet->Phi()<<endl;
    }

    // check for subleading jet
    if(fSubLeadingJet) {
      excludeInEtaSub = fSubLeadingJet->Eta();
      excludeInPhiSub = fSubLeadingJet->Phi();
      //cout<<"subleading: pt = "<<fSubLeadingJet->Pt()<<"  eta = "<<fSubLeadingJet->Eta()<<"  phi = "<<fSubLeadingJet->Phi()<<endl;
    }
  } // leading jets

  // loop over tracks
  TRandom3 *rand = new TRandom3();
  int nTOT = 0, nA = 0, nB = 0;
  int nTrack = mPicoDst->numberOfTracks();
  for(int i = 0; i < nTrack; i++) {
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 5.0 GeV
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin - flipped syntax back in April 2020
    if(doTPCptassocBin) {
      if(ptbin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
      if(ptbin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(ptbin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(ptbin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(ptbin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(ptbin == 5) { if((pt > 2.00) && (pt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(ptbin == 6) { if((pt > 3.00) && (pt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(ptbin == 7) { if((pt > 4.00) && (pt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }

    // remove strip only when we have a leading jet
    // Method1: kRemoveEtaStrip
    if(fTPCEPmethod == kRemoveEtaStrip){
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && 
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    } else if(fTPCEPmethod == kRemoveEtaPhiCone){
      // remove cone (in eta and phi) around leading jet
      // Method2: kRemoveEtaPhiCone - FIXME found bug May25
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    } else if(fTPCEPmethod == kRemoveLeadingJetConstituents){
      // remove tracks above 2 GeV in cone around leading jet
      // Method3: kRemoveLeadingJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
    } else if(fTPCEPmethod == kRemoveEtaStripLeadSub){
      // remove strip only when we have a leading + subleading jet
      // Method4: kRemoveEtaStripLeadSub
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEtaSub) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    } else if(fTPCEPmethod == kRemoveEtaPhiConeLeadSub){
      // remove cone (in eta and phi) around leading + subleading jet
      // Method5: kRemoveEtaPhiConeLeadSub
      double deltaR    = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi-excludeInPhiSub)*(phi-excludeInPhiSub));
      if((fLeadingJet)    && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR    < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaRSub < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    } else if(fTPCEPmethod == kRemoveLeadingSubJetConstituents){
      // remove tracks above 2 GeV in cone around leading + subleading jet
      // Method6: kRemoveLeadingSubJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi - excludeInPhiSub)*(phi - excludeInPhiSub));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaRSub < fJetRad)) continue;
    } else {
      // DO NOTHING! nothing is removed...
    }

    // configure track weight when performing Q-vector summation
    double trackweight;
    if(fTrackWeight == kNoWeight) {
      trackweight = 1.0;
    } else if(fTrackWeight == kPtLinearWeight) {
      trackweight = pt;
    } else if(fTrackWeight == kPtLinear2Const5Weight) {
      if(pt <= 2.0) trackweight = pt;
      if(pt >  2.0) trackweight = 2.0;
    } else {
      // nothing choosen, so don't use weight
      trackweight = 1.0;
    }

    // generate random distribution from 0 -> 1: and split subevents for [0,0.5] and [0.5, 1]
    double randomNum = rand->Rndm();
    ////double randomNum = gRandom->Rndm();  // > 0.5?

    // split up Q-vectors into 2 random TPC sub-events (A and B)
    if(doTPCptassocBin) {
      if(randomNum >= 0.5) {
        // sum up q-vectors for random event A
        mQtpcpx += trackweight * cos(phi * order);
        mQtpcpy += trackweight * sin(phi * order);
        nA++;
      }
      if(randomNum < 0.5) {
        // sum up q-vectors for random event B
        mQtpcnx += trackweight * cos(phi * order);
        mQtpcny += trackweight * sin(phi * order);
        nB++;
      }
    }  // pt-dependent mode

    // non-pt dependent mode
    // split up Q-vectors to +/- eta regions
    if(!doTPCptassocBin) {
      if(eta < 0) {
        // sum up q-vectors on negative eta side
        mQtpcnx += trackweight * cos(phi * order);
        mQtpcny += trackweight * sin(phi * order);
        ntracksNEG++;
      }
      if(eta >= 0) {
        // sum up q-vectors on positive eta side
        mQtpcpx += trackweight * cos(phi * order);
        mQtpcpy += trackweight * sin(phi * order);
        ntracksPOS++;
      }
    }  

    // combined TPC event plane q-vectors
    mQtpcX += trackweight * cos(phi * order);
    mQtpcY += trackweight * sin(phi * order);
    nTOT++;

  } // end of track loop

  // test statements
  //cout<<"ntracksNEG = "<<ntracksNEG<<"  ntracksPOS = "<<ntracksPOS<<endl;
  //cout<<"mQtpcpx = "<<mQtpcpx<<"  mQtpcpy = "<<mQtpcpy<<"  mQtpcnx = "<<mQtpcnx<<"  mQtpcny = "<<mQtpcny<<endl;
  //cout<<"nTOT = "<<nTOT<<"  nA = "<<nA<<"  nB = "<<nB<<endl;

  // set q-vector components 
  mQtpcn.Set(mQtpcnx, mQtpcny);
  mQtpcp.Set(mQtpcpx, mQtpcpy);
  mQtpc.Set(mQtpcX, mQtpcY);

  // test..
  TVector2 mQtpcComb;
  mQtpcComb.Set(mQtpcnx + mQtpcpx, mQtpcny + mQtpcpy);
  fEPTPCcomb = mQtpcComb.Phi() / order;

  // Calculate the Event Plane
  fEPTPCn = mQtpcn.Phi() / order;
  fEPTPCp = mQtpcp.Phi() / order;
  fEPTPC = mQtpc.Phi() / order;

  // make sure event plane is contrained from [0, pi]
  if((fEPTPCn > pi) || (fEPTPCn < 0)) cout<<"fEPTPCn out of range..... "<<fEPTPCn<<endl;
  if((fEPTPCp > pi) || (fEPTPCp < 0)) cout<<"fEPTPCp out of range..... "<<fEPTPCp<<endl;
  if((fEPTPC  > pi) || (fEPTPC  < 0)) cout<<"fEPTPC out of range..... "<<fEPTPC<<endl;
  if(fEPTPCn > pi) fEPTPCn -= pi;
  if(fEPTPCn <  0) fEPTPCn += pi;
  if(fEPTPCp > pi) fEPTPCp -= pi;
  if(fEPTPCp <  0) fEPTPCp += pi;
  if(fEPTPC > pi)  fEPTPC  -= pi;
  if(fEPTPC <  0)  fEPTPC  += pi;

  // combine x-y vectors for neg and pos Eta ranges (cross-check)
  //double tpc2 = (0.5*TMath::ATan2(mQtpcY, mQtpcX));

  // standard event plane distributions as function of centrality
  fHistEPTPCn->Fill(fCentralityScaled, fEPTPCn);
  fHistEPTPCp->Fill(fCentralityScaled, fEPTPCp);
}
//
// Fill event plane resolution histograms
//_____________________________________________________________________________
void StJetHadronCorrelations::CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP, Double_t bbc1, Double_t zdc1)
{ 
    // fill the profiles for the resolution parameters
    // R2 resolution for 2nd order event plane
    fProfV2Resolution[ref9]->Fill(2., TMath::Cos(2.*(bbc - tpc)));
    fProfV2Resolution[ref9]->Fill(3., TMath::Cos(2.*(bbc - zdc)));
    fProfV2Resolution[ref9]->Fill(4., TMath::Cos(2.*(tpc - bbc)));  // bin2
    fProfV2Resolution[ref9]->Fill(5., TMath::Cos(2.*(tpc - zdc)));  
    fProfV2Resolution[ref9]->Fill(6., TMath::Cos(2.*(zdc - tpc)));  // bin6
    fProfV2Resolution[ref9]->Fill(7., TMath::Cos(2.*(zdc - bbc)));  // bin3
    fProfV2Resolution[ref9]->Fill(8., TMath::Cos(2.*(bbc - tpcN)));
    fProfV2Resolution[ref9]->Fill(9., TMath::Cos(2.*(bbc - tpcP)));
    fProfV2Resolution[ref9]->Fill(10., TMath::Cos(2.*(zdc - tpcN)));
    fProfV2Resolution[ref9]->Fill(11., TMath::Cos(2.*(zdc - tpcP)));
    fProfV2Resolution[ref9]->Fill(12., TMath::Cos(2.*(tpcP - tpcN)));
    fProfV2Resolution[ref9]->Fill(17., TMath::Cos(2.*(bbc1 - tpc)));
    fProfV2Resolution[ref9]->Fill(18., TMath::Cos(2.*(bbc1 - tpcN)));
    fProfV2Resolution[ref9]->Fill(19., TMath::Cos(2.*(bbc1 - tpcP)));
    fProfV2Resolution[ref9]->Fill(20., TMath::Cos(2.*(bbc1 - zdc1)));
    fProfV2Resolution[ref9]->Fill(21., TMath::Cos(2.*(zdc1 - tpc)));
    fProfV2Resolution[ref9]->Fill(22., TMath::Cos(2.*(zdc1 - tpcN)));
    fProfV2Resolution[ref9]->Fill(23., TMath::Cos(2.*(zdc1 - tpcP)));

    // R3 resolution for 2nd order event plane
    fProfV3Resolution[ref9]->Fill(2., TMath::Cos(3.*(bbc - tpc)));
    fProfV3Resolution[ref9]->Fill(3., TMath::Cos(3.*(bbc - zdc)));
    fProfV3Resolution[ref9]->Fill(4., TMath::Cos(3.*(tpc - bbc)));
    fProfV3Resolution[ref9]->Fill(5., TMath::Cos(3.*(tpc - zdc)));
    fProfV3Resolution[ref9]->Fill(6., TMath::Cos(3.*(zdc - tpc)));
    fProfV3Resolution[ref9]->Fill(7., TMath::Cos(3.*(zdc - bbc)));
    fProfV3Resolution[ref9]->Fill(8., TMath::Cos(3.*(bbc - tpcN)));
    fProfV3Resolution[ref9]->Fill(9., TMath::Cos(3.*(bbc - tpcP)));
    fProfV3Resolution[ref9]->Fill(10., TMath::Cos(3.*(zdc - tpcN)));
    fProfV3Resolution[ref9]->Fill(11., TMath::Cos(3.*(zdc - tpcP)));
    fProfV3Resolution[ref9]->Fill(12., TMath::Cos(3.*(tpcP - tpcN)));
    fProfV3Resolution[ref9]->Fill(17., TMath::Cos(3.*(bbc1 - tpc)));
    fProfV3Resolution[ref9]->Fill(18., TMath::Cos(3.*(bbc1 - tpcN)));
    fProfV3Resolution[ref9]->Fill(19., TMath::Cos(3.*(bbc1 - tpcP)));
    fProfV3Resolution[ref9]->Fill(20., TMath::Cos(3.*(bbc1 - zdc1)));
    fProfV3Resolution[ref9]->Fill(21., TMath::Cos(3.*(zdc1 - tpc)));
    fProfV3Resolution[ref9]->Fill(22., TMath::Cos(3.*(zdc1 - tpcN)));
    fProfV3Resolution[ref9]->Fill(23., TMath::Cos(3.*(zdc1 - tpcP)));

    // R4 resolution for 2nd order event plane
    fProfV4Resolution[ref9]->Fill(2., TMath::Cos(4.*(bbc - tpc)));
    fProfV4Resolution[ref9]->Fill(3., TMath::Cos(4.*(bbc - zdc)));
    fProfV4Resolution[ref9]->Fill(4., TMath::Cos(4.*(tpc - bbc)));
    fProfV4Resolution[ref9]->Fill(5., TMath::Cos(4.*(tpc - zdc)));
    fProfV4Resolution[ref9]->Fill(6., TMath::Cos(4.*(zdc - tpc)));
    fProfV4Resolution[ref9]->Fill(7., TMath::Cos(4.*(zdc - bbc)));
    fProfV4Resolution[ref9]->Fill(8., TMath::Cos(4.*(bbc - tpcN)));
    fProfV4Resolution[ref9]->Fill(9., TMath::Cos(4.*(bbc - tpcP)));
    fProfV4Resolution[ref9]->Fill(10., TMath::Cos(4.*(zdc - tpcN)));
    fProfV4Resolution[ref9]->Fill(11., TMath::Cos(4.*(zdc - tpcP)));
    fProfV4Resolution[ref9]->Fill(12., TMath::Cos(4.*(tpcP - tpcN)));
    fProfV4Resolution[ref9]->Fill(17., TMath::Cos(4.*(bbc1 - tpc)));
    fProfV4Resolution[ref9]->Fill(18., TMath::Cos(4.*(bbc1 - tpcN)));
    fProfV4Resolution[ref9]->Fill(19., TMath::Cos(4.*(bbc1 - tpcP)));
    fProfV4Resolution[ref9]->Fill(20., TMath::Cos(4.*(bbc1 - zdc1)));
    fProfV4Resolution[ref9]->Fill(21., TMath::Cos(4.*(zdc1 - tpc)));
    fProfV4Resolution[ref9]->Fill(22., TMath::Cos(4.*(zdc1 - tpcN)));
    fProfV4Resolution[ref9]->Fill(23., TMath::Cos(4.*(zdc1 - tpcP)));

    // R5 resolution for 2nd order event plane
    fProfV5Resolution[ref9]->Fill(2., TMath::Cos(5.*(bbc - tpc)));
    fProfV5Resolution[ref9]->Fill(3., TMath::Cos(5.*(bbc - zdc)));
    fProfV5Resolution[ref9]->Fill(4., TMath::Cos(5.*(tpc - bbc)));
    fProfV5Resolution[ref9]->Fill(5., TMath::Cos(5.*(tpc - zdc)));
    fProfV5Resolution[ref9]->Fill(6., TMath::Cos(5.*(zdc - tpc)));
    fProfV5Resolution[ref9]->Fill(7., TMath::Cos(5.*(zdc - bbc)));
    fProfV5Resolution[ref9]->Fill(8., TMath::Cos(5.*(bbc - tpcN)));
    fProfV5Resolution[ref9]->Fill(9., TMath::Cos(5.*(bbc - tpcP)));
    fProfV5Resolution[ref9]->Fill(10., TMath::Cos(5.*(zdc - tpcN)));
    fProfV5Resolution[ref9]->Fill(11., TMath::Cos(5.*(zdc - tpcP)));
    fProfV5Resolution[ref9]->Fill(12., TMath::Cos(5.*(tpcP - tpcN)));
    fProfV5Resolution[ref9]->Fill(17., TMath::Cos(5.*(bbc1 - tpc)));
    fProfV5Resolution[ref9]->Fill(18., TMath::Cos(5.*(bbc1 - tpcN)));
    fProfV5Resolution[ref9]->Fill(19., TMath::Cos(5.*(bbc1 - tpcP)));
    fProfV5Resolution[ref9]->Fill(20., TMath::Cos(5.*(bbc1 - zdc1)));
    fProfV5Resolution[ref9]->Fill(21., TMath::Cos(5.*(zdc1 - tpc)));
    fProfV5Resolution[ref9]->Fill(22., TMath::Cos(5.*(zdc1 - tpcN)));
    fProfV5Resolution[ref9]->Fill(23., TMath::Cos(5.*(zdc1 - tpcP)));

    // for the resolution of the combined vzero event plane, use two tpc halves as uncorrelated subdetectors
} 
//
// function to calculate: event plane chi
//_____________________________________________________________________________
Double_t StJetHadronCorrelations::CalculateEventPlaneChi(Double_t res) {
  // return chi for given resolution to combine event plane estimates from two subevents
  // see Phys. Rev. C no. CS6346 (http://arxiv.org/abs/nucl-ex/9805001)
  Double_t chi(2.), delta(1.), con((TMath::Sqrt(TMath::Pi()))/(2.*TMath::Sqrt(2)));
  for (Int_t i(0); i < 15; i++) {
    chi = ((con*chi*TMath::Exp(-chi*chi/4.)*(TMath::BesselI0(chi*chi/4.)+TMath::BesselI1(chi*chi/4.))) < res) ? chi + delta : chi - delta;
    delta = delta / 2.;
  }

  return chi;
}
//
// track QA function to fill some histograms with track information
//_____________________________________________________________________________
void StJetHadronCorrelations::TrackQA()
{
  // get # of tracks
  int nTrack = mPicoDst->numberOfTracks();
  double pi = 1.0*TMath::Pi();

  // loop over all tracks
  for(int i = 0; i < nTrack; i++) {
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // for heavy ion collisions, get EP and fill QA plots
    // get angle between tracks and event plane
    double dEPtrk = (!doppAnalysis) ? RelativeEPJET(phi, TPC_PSI2) : 0;
    //cout<<"dEPtrk = "<<dEPtrk<<"  phi = "<<phi<<"  EP2 = "<<TPC_PSI2<<endl;
    if((pt > 0.20) && (pt <= 0.5)) hdEPtrk[0]->Fill(dEPtrk);
    if((pt > 0.50) && (pt <= 1.0)) hdEPtrk[1]->Fill(dEPtrk);
    if((pt > 1.00) && (pt <= 1.5)) hdEPtrk[2]->Fill(dEPtrk);
    if((pt > 1.50) && (pt <= 2.0)) hdEPtrk[3]->Fill(dEPtrk);
    if((pt > 2.00) && (pt <= 20.)) hdEPtrk[4]->Fill(dEPtrk);

    // fill other track QA plots
    hTrackPhi[ref9]->Fill(phi);
    hTrackEta[ref9]->Fill(eta);
    hTrackPt[ref9]->Fill(pt);
    hTrackEtavsPhi->Fill(phi, eta);

  }
}
//
// function filling trigger arrays with booleans set to kTRUE for fired triggers
//_________________________________________________________________________
void StJetHadronCorrelations::FillTowerTriggersArr() {
  // tower - HT trigger types array: zero these out - so they are refreshed for each event
  for(int i = 0; i < 4800; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }

  // get number of Emc triggers
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get pointer to trigger
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // emc trigger parameters
    int emcTrigID = emcTrig->id();

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();

    // set trigger level type
    if(isHT1) fTowerToTriggerTypeHT1[emcTrigID] = kTRUE;
    if(isHT2) fTowerToTriggerTypeHT2[emcTrigID] = kTRUE;
    if(isHT3) fTowerToTriggerTypeHT3[emcTrigID] = kTRUE;

    //cout<<"i = "<<i<<"  EmcTrigID = "<<emcTrigID<<"  adc = "<<emcTrig->adc()<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3<<endl;
  }

}
//
// function to require that a jet constituent tower fired a HT trigger
//___________________________________________________________________________________________
Bool_t StJetHadronCorrelations::DidTowerConstituentFireTrigger(StJet *jet) {  
  // tower constituent fired trigger
  Bool_t mFiredTrigger = kFALSE;

  // loop over constituent towers
  for(int itow = 0; itow < jet->GetNumberOfClusters(); itow++) {
    int towerIndex = jet->ClusterAt(itow);

    // get tower pointer
    StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towerIndex));
    if(!tow) continue;
    
    // tower ID: get from index of array shifted by +1
    int towID = towerIndex + 1;
    if(towID < 0) continue;

    // change flag to true if jet tower fired trigger
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && fTowerToTriggerTypeHT1[towerIndex]) mFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && fTowerToTriggerTypeHT2[towerIndex]) mFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && fTowerToTriggerTypeHT3[towerIndex]) mFiredTrigger = kTRUE;
  
  } // tower constituent loop

  return mFiredTrigger;
}  
//
// FIXME TODO - still working on how this function will perform.. May 31, 2019
// function to check if a bad tower fired the events HT trigger
//___________________________________________________________________________________________
Bool_t StJetHadronCorrelations::DidBadTowerFireTrigger() {
  // bad/dead tower fired trigger
  bool mBadTowerFiredTrigger = kFALSE;
  bool mBadTowerFiring = kFALSE;

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // update switch
    mBadTowerFiredTrigger = kFALSE;

    // get tower pointer
    StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tow) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: get from index of array shifted by +1
    int towID = itow + 1;
    int towIndex = towID - 1;
    if(towID < 0) continue; // double check these aren't still in the event list

    // check if tower is bad or dead - functions return kTRUE if ok and kFALSE if NOT dead 
    // isTowerOK(towID) = kTRUE if tower is OK and not dead
    // isTowerDead(towID) = kTRUE if tower is dead, = kFALSE if tower is NOT dead
//    bool isTowOk = (mBaseMaker->IsTowerOK(towID) && !mBaseMaker->IsTowerDead(towID)); // FIXME
    bool isTowOk = (mBaseMaker->IsTowerOK(towID)); // FIXME

    // change flag to true if jet tower fired trigger
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && fTowerToTriggerTypeHT1[towIndex] && !isTowOk) mBadTowerFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && fTowerToTriggerTypeHT2[towIndex] && !isTowOk) mBadTowerFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && fTowerToTriggerTypeHT3[towIndex] && !isTowOk) mBadTowerFiredTrigger = kTRUE;

    // for bad tower firings: fill histo, inform user, return kTRUE
    if(mBadTowerFiredTrigger == kTRUE) {
      hBadTowerFiredTrigger->Fill(towID);
      mBadTowerFiring = kTRUE;
      cout<<"Bad tower fired trigger, towID: "<<towID<<endl;
      //break;
    }

  } // tower loop

  // kTRUE: bad tower fired trigger
  //return mBadTowerFiredTrigger;
  return mBadTowerFiring;
}
//
// function to calculate Jet v2 using the event plane method
//____________________________________________________________________________________________
void StJetHadronCorrelations::GetJetV2(StJet *jet, Double_t EPangle, Int_t ptAssocBin)
{
  // get centrality bin
  int centBin = (!doppAnalysis) ? Get4CentBin(fCentralityScaled) : 0; // this secures it will still run for pp
  if(centBin < 0) return; // kStOK;

  // get jet pt bin and value
  double jetPt = -99.;
  if(fCorrJetPt) { jetPt = jet->Pt() - jet->Area()*fRhoVal;
  } else { jetPt = jet->Pt(); }
  int jetPtBin = GetJetPtBin(jetPt);
  if(jetPtBin < 0) return; // kStOK;

  // require tower and or track bias for jet
  //if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) return kStOK;

  // check that jet contains a tower that fired the trigger
  //if(!DidTowerConstituentFireTrigger(jet)) return kStOK;

  // jet parameters
  double jetPhi = jet->Phi();

  // get event plane bin
  double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, EPangle) : -99.; // CORRECTED event plane angle - STEP3
  int EPBin = (!doppAnalysis) ? GetJetEPBin(dEP) : 0; // this secures it will still run for pp, and use EP bin = 0
  if(EPBin < 0) return; // kStOK;

  // fill histogram
  fProfJetV2[jetPtBin][centBin][EPBin]->Fill(ptAssocBin, TMath::Cos(2.*(jetPhi - EPangle)));
  fProfJetV2[jetPtBin][centBin][3]->Fill(ptAssocBin, TMath::Cos(2.*(jetPhi - EPangle)));

}
//
// Function to add event pools to output file
//____________________________________________________________________
void StJetHadronCorrelations::AddEventPoolsToOutput(Double_t minCent, Double_t maxCent, Double_t minZvtx, Double_t maxZvtx, Double_t minPsi2, Double_t maxPsi2, Double_t minPt, Double_t maxPt)
{
  // create vector of vectors
  std::vector<Double_t> binVec;
  binVec.push_back(minCent);
  binVec.push_back(maxCent);
  binVec.push_back(minZvtx);
  binVec.push_back(maxZvtx);
  binVec.push_back(minPsi2);
  binVec.push_back(maxPsi2);
  binVec.push_back(minPt);
  binVec.push_back(maxPt);

  // add vectors to output list
  fEventPoolOutputList.push_back(binVec);
}
