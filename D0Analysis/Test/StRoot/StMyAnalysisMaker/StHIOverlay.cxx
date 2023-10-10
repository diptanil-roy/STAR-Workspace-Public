// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StHIOverlay.h"
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


// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StHIOverlay)

//________________________________________________________________________
StHIOverlay::StHIOverlay(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* filename = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StHIOverlay::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fCorrJetPt = kFALSE;
  fMinPtJet = 0.0;
  fTrackBias = 0.0;
  fTowerBias = 0.0;
  fJetRad = 0.4;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;
  fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0;  fTowerPhiMaxCut = 2.0*TMath::Pi();
  fCentralityScaled = 0.; ref16 = -99; ref9 = -99; // FIXME - maybe not make global
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  fRho = 0x0;
  fRhoVal = 0;
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;

  fMCFileListName = filename;

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  mBemcGeom = 0x0;

  fMinJetTrackPt = 0.2; fMaxJetTrackPt = 30.0;
  fJetTrackEtaMin = -1.0; fJetTrackEtaMax = 1.0;
  fJetTrackPhiMin = 0.0; fJetTrackPhiMax = 2.0*TMath::Pi();
  fJetTrackDCAcut = 3.0;
  fJetTracknHitsFit = 15;
  fJetTracknHitsRatio = 0.52;

  fJetTowerEMin = 0.2; fJetTowerEMax = 100.0;
  fJetTowerEtaMin = -1.0; fJetTowerEtaMax = 1.0;
  fJetTowerPhiMin = 0.0; fJetTowerPhiMax = 2.0*TMath::Pi();
  mTowerEnergyMin = 0.2;

  fSetNumberOfEvents = 0;
  fNumberofeventsoverlayed = 0;

  for (int i = 0; i < 100; i++){
    fMCEventTracks[i].clear();
    fMCEventTowers[i].clear();
    fRecoEventTracks[i].clear();
    fRecoEventTowers[i].clear();
    fMCD0Information[i] = {};
    fRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0,0,0);
  }

}

//
//________________________________________________________________________
StHIOverlay::~StHIOverlay()
{ /*  */

}

//
//________________________________________________________________________
Int_t StHIOverlay::Init() {
  StJetFrameworkPicoBase::Init();

  // position object for Emc
  mBemcGeom = StEmcGeom::instance("bemc");
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  ifstream filelistforMCEvents(fMCFileListName.Data());

  if (!filelistforMCEvents.is_open()) {
    LOG_ERROR << "No MC File List! Exiting!" << endm; 
    return kStOk;
  }

  string line;

  // cout << "Files read in: " << endl;

  while(getline(filelistforMCEvents,line))
  {
    TString s(line);

    // cout << s << endl;
    filenamesforHIOverlay.push_back(s);  
  }

  TFile f("/star/u/droy1/Y2019/STAR/Momentum_resolution_SL16d.root");
  fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPion");
  fKaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaon");
  fProtonMomResolution = (TF1*)f.Get("fProton")->Clone("fProton");

  TFile effweight("/star/u/droy1/Y2019/STAR/EffWeightsInCentralityBins.root");

  fPionWeight[0] = (TGraph *)effweight.Get("Pion_0_10");
  fKaonWeight[0] = (TGraph *)effweight.Get("Kaon_0_10");
  fProtonWeight[0] = (TGraph *)effweight.Get("Proton_0_10");
  fAProtonWeight[0] = (TGraph *)effweight.Get("AProton_0_10");

  fPionWeight[1] = (TGraph *)effweight.Get("Pion_10_40");
  fKaonWeight[1] = (TGraph *)effweight.Get("Kaon_10_40");
  fProtonWeight[1] = (TGraph *)effweight.Get("Proton_10_40");
  fAProtonWeight[1] = (TGraph *)effweight.Get("AProton_10_40");

  fPionWeight[2] = (TGraph *)effweight.Get("Pion_40_80");
  fKaonWeight[2] = (TGraph *)effweight.Get("Kaon_40_80");
  fProtonWeight[2] = (TGraph *)effweight.Get("Proton_40_80");
  fAProtonWeight[2] = (TGraph *)effweight.Get("AProton_40_80");

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StHIOverlay::Finish() { 
  cout << "StHIOverlay::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StHIOverlay::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StHIOverlay::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StHIOverlay::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  for (int i = 0; i < 100; i++){
    fMCEventTracks[i].clear();
    fMCEventTowers[i].clear();
    fRecoEventTracks[i].clear();
    fRecoEventTowers[i].clear();
    fMCD0Information[i] = {};
    fRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0,0,0);
  }

  fNumberofeventsoverlayed = 0;

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
  }

  //for(int i = 1; i<4801; i++) {  if(!mBaseMaker->IsTowerOK(i))  cout<<"tower: "<<i<<" is not good!!"<<endl;  }
  // can check for bad towers like this (if needed):
  // bool isTowOK = mBaseMaker->IsTowerOK(towerID);
  // ===========================================================================================

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  double zVtx_VPD = mPicoEvent->vzVpd();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


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

  // cout << "Centrality = " << grefMult << "\t" << ref16 << "\t" << cent16 << "\t" << fCentralityScaled << endl;


  if (fCentralityScaled >= 0 && fCentralityScaled < 10) centralitybinforefficiency = 0;
  else if (fCentralityScaled >= 10 && fCentralityScaled < 40) centralitybinforefficiency = 1;
  else if (fCentralityScaled >= 40 && fCentralityScaled < 80) centralitybinforefficiency = 2;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // cout << "1" << endl;

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // cout << "2" << endl;

  // ============================ end of CENTRALITY ============================== //

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5) return kStOK;

  // cout << "3" << endl;

  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // cout << "4" << endl;

  int arrMB5_Run14[]  = {450005, 450015, 450025, 450050, 450060};

  bool matchMB = kFALSE;

  for(int i = 0; i < sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14); i++) {
    if(mPicoEvent->isTrigger(arrMB5_Run14[i])) matchMB = kTRUE;
    if(matchMB) break;
  }

  if (!matchMB) return kStOk;

  // cout << "5" << endl;

  if (abs(zVtx) > 6.) return kStOk;

  // cout << "6" << endl;

  if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.) return kStOk;

  // cout << "7" << endl;

  if (abs(zVtx - zVtx_VPD) > 3) return kStOk;

  // cout << "8" << endl;

  // cout << "Centrality = " << fCentralityScaled << endl;
  // pick a random file
  // PickARandomFile();

  // run Towers:
  SampleMCEvents();

  return kStOK;
}

void StHIOverlay::PickARandomFile(){


}


void StHIOverlay::SampleMCEvents(){

  TRandom3 *r1 = new TRandom3(0);
  int filenumber = r1->Integer(filenamesforHIOverlay.size());

  // cout << filenumber << "\t" << filenamesforHIOverlay[filenumber] << endl;

  TFile *f = new TFile(filenamesforHIOverlay[filenumber].Data());

  fMCPico = (TTree *)f->Get("PicoDst");
  // fMCPico->SetDirectory(0);

  double pi0mass = Pico::mMass[0]; // GeV

  int nentries = fMCPico->GetEntriesFast();

  

  // cout << "Entries in PicoDst = " << nentries << endl;

  // fMCPico->Print();

  TRandom3 *r2 = new TRandom3(0);

  vector<int> randomlisttoevents;

  randomlisttoevents.clear();

  // for (int i = 0; i < fSetNumberOfEvents; i++){ //This has been FIXED to be able to automate the number of MC events for overlay studies
  //   randomlisttoevents.push_back(r2->Integer(nentries));
  // }

  int numberofmceventstoconvulate = fSetNumberOfEvents;

  if (centralitybinforefficiency == 0) numberofmceventstoconvulate = 2*fSetNumberOfEvents;
  else if (centralitybinforefficiency == 1) numberofmceventstoconvulate = fSetNumberOfEvents;
  else if (centralitybinforefficiency == 2) numberofmceventstoconvulate = 3*fSetNumberOfEvents;

  while (randomlisttoevents.size()!=numberofmceventstoconvulate){
    int num = r2->Integer(nentries);
    if (std::find(randomlisttoevents.begin(), randomlisttoevents.end(), num) == randomlisttoevents.end()) { //New entry every time...
      randomlisttoevents.push_back(num);
    }
  }

  // for (int i = 0; i < fSetNumberOfEvents; i++){cout << randomlisttoevents[i] << endl;}

  // Structure from PicoDst that we need to mimic its functionalities.
  // I reckon this will be the most memory intensive process of the whole operation.
  // There is significant I/O happening here.

  // Investigate ways to cut I/O here.

  fMCPico->SetMakeClass(1);

  // Fixed size dimensions of array or collections stored in the TTree if any.
  static constexpr Int_t kMaxEvent = 1;
  static constexpr Int_t kMaxTrack = 107;
  static constexpr Int_t kMaxEmcTrigger = 89;
  static constexpr Int_t kMaxMtdTrigger = 1;
  static constexpr Int_t kMaxBTowHit = 4800;
  static constexpr Int_t kMaxBTofHit = 83;
  static constexpr Int_t kMaxMtdHit = 4;
  static constexpr Int_t kMaxBbcHit = 1;
  static constexpr Int_t kMaxEpdHit = 1;
  static constexpr Int_t kMaxFmsHit = 1;
  static constexpr Int_t kMaxEmcPidTraits = 7;
  static constexpr Int_t kMaxBTofPidTraits = 31;
  static constexpr Int_t kMaxMtdPidTraits = 1;
  static constexpr Int_t kMaxTrackCovMatrix = 107;
  static constexpr Int_t kMaxBEmcSmdEHit = 1;
  static constexpr Int_t kMaxBEmcSmdPHit = 1;
  static constexpr Int_t kMaxETofHit = 1;
  static constexpr Int_t kMaxETofPidTraits = 1;
  static constexpr Int_t kMaxMcVertex = 172;
  static constexpr Int_t kMaxMcTrack = 446;

  // Declaration of leaf types
  Int_t           Event_;
  Int_t           Event_mRunId[kMaxEvent];   //[Event_]
  Int_t           Event_mEventId[kMaxEvent];   //[Event_]
  UShort_t        Event_mFillId[kMaxEvent];   //[Event_]
  Float_t         Event_mBField[kMaxEvent];   //[Event_]
  Int_t           Event_mTime[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexX[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexY[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexZ[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexErrorX[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexErrorY[kMaxEvent];   //[Event_]
  Float_t         Event_mPrimaryVertexErrorZ[kMaxEvent];   //[Event_]
  Float_t         Event_mRanking[kMaxEvent];   //[Event_]
  UShort_t        Event_mNBEMCMatch[kMaxEvent];   //[Event_]
  UShort_t        Event_mNBTOFMatch[kMaxEvent];   //[Event_]
  vector<unsigned int> Event_mTriggerIds[kMaxEvent];
  UShort_t        Event_mRefMultFtpcEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultFtpcWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultNeg[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultPos[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2NegEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2PosEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2NegWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult2PosWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult3NegEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult3PosEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult3NegWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult3PosWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult4NegEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult4PosEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult4NegWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMult4PosWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfNegEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfPosEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfNegWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mRefMultHalfPosWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mGRefMult[kMaxEvent];   //[Event_]
  UShort_t        Event_mNumberOfGlobalTracks[kMaxEvent];   //[Event_]
  UShort_t        Event_mbTofTrayMultiplicity[kMaxEvent];   //[Event_]
  UShort_t        Event_mNHitsHFT[kMaxEvent][4];   //[Event_]
  UChar_t         Event_mNVpdHitsEast[kMaxEvent];   //[Event_]
  UChar_t         Event_mNVpdHitsWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mNTofT0[kMaxEvent];   //[Event_]
  Float_t         Event_mVzVpd[kMaxEvent];   //[Event_]
  UInt_t          Event_mZDCx[kMaxEvent];   //[Event_]
  UInt_t          Event_mBBCx[kMaxEvent];   //[Event_]
  Float_t         Event_mBackgroundRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcBlueBackgroundRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcYellowBackgroundRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcEastRate[kMaxEvent];   //[Event_]
  Float_t         Event_mBbcWestRate[kMaxEvent];   //[Event_]
  Float_t         Event_mZdcEastRate[kMaxEvent];   //[Event_]
  Float_t         Event_mZdcWestRate[kMaxEvent];   //[Event_]
  UShort_t        Event_mZdcSumAdcEast[kMaxEvent];   //[Event_]
  UShort_t        Event_mZdcSumAdcWest[kMaxEvent];   //[Event_]
  UShort_t        Event_mZdcSmdEastHorizontal[kMaxEvent][8];   //[Event_]
  UShort_t        Event_mZdcSmdEastVertical[kMaxEvent][8];   //[Event_]
  UShort_t        Event_mZdcSmdWestHorizontal[kMaxEvent][8];   //[Event_]
  UShort_t        Event_mZdcSmdWestVertical[kMaxEvent][8];   //[Event_]
  UShort_t        Event_mBbcAdcEast[kMaxEvent][24];   //[Event_]
  UShort_t        Event_mBbcAdcWest[kMaxEvent][24];   //[Event_]
  UChar_t         Event_mHighTowerThreshold[kMaxEvent][4];   //[Event_]
  UChar_t         Event_mJetPatchThreshold[kMaxEvent][4];   //[Event_]
  UChar_t         Event_mBunchCrossId[kMaxEvent];   //[Event_]
  UShort_t        Event_mETofHitMultiplicity[kMaxEvent];   //[Event_]
  UShort_t        Event_mETofDigiMultiplicity[kMaxEvent];   //[Event_]
  UShort_t        Event_mNumberOfPrimaryTracks[kMaxEvent];   //[Event_]
  UShort_t        Event_mZdcUnAttenuated[kMaxEvent][2];   //[Event_]
  Int_t           Track_;
  UShort_t        Track_mId[kMaxTrack];   //[Track_]
  UShort_t        Track_mChi2[kMaxTrack];   //[Track_]
  Float_t         Track_mPMomentumX[kMaxTrack];   //[Track_]
  Float_t         Track_mPMomentumY[kMaxTrack];   //[Track_]
  Float_t         Track_mPMomentumZ[kMaxTrack];   //[Track_]
  Float_t         Track_mGMomentumX[kMaxTrack];   //[Track_]
  Float_t         Track_mGMomentumY[kMaxTrack];   //[Track_]
  Float_t         Track_mGMomentumZ[kMaxTrack];   //[Track_]
  Float_t         Track_mOriginX[kMaxTrack];   //[Track_]
  Float_t         Track_mOriginY[kMaxTrack];   //[Track_]
  Float_t         Track_mOriginZ[kMaxTrack];   //[Track_]
  Float16_t       Track_mDedx[kMaxTrack];   //[Track_]
  Float16_t       Track_mDedxError[kMaxTrack];   //[Track_]
  Char_t          Track_mNHitsFit[kMaxTrack];   //[Track_]
  UChar_t         Track_mNHitsMax[kMaxTrack];   //[Track_]
  UChar_t         Track_mNHitsDedx[kMaxTrack];   //[Track_]
  Short_t         Track_mNSigmaPion[kMaxTrack];   //[Track_]
  Short_t         Track_mNSigmaKaon[kMaxTrack];   //[Track_]
  Short_t         Track_mNSigmaProton[kMaxTrack];   //[Track_]
  Short_t         Track_mNSigmaElectron[kMaxTrack];   //[Track_]
  UInt_t          Track_mTopologyMap[kMaxTrack][2];   //[Track_]
  Short_t         Track_mBEmcPidTraitsIndex[kMaxTrack];   //[Track_]
  Short_t         Track_mBTofPidTraitsIndex[kMaxTrack];   //[Track_]
  Short_t         Track_mMtdPidTraitsIndex[kMaxTrack];   //[Track_]
  Short_t         Track_mETofPidTraitsIndex[kMaxTrack];   //[Track_]
  Short_t         Track_mBEmcMatchedTowerIndex[kMaxTrack];   //[Track_]
  ULong64_t       Track_mTopoMap_iTpc[kMaxTrack];   //[Track_]
  UShort_t        Track_mIdTruth[kMaxTrack];   //[Track_]
  UShort_t        Track_mQATruth[kMaxTrack];   //[Track_]
  Char_t          Track_mVertexIndex[kMaxTrack];   //[Track_]
  Int_t           EmcTrigger_;
  UChar_t         EmcTrigger_mFlag[kMaxEmcTrigger];   //[EmcTrigger_]
  UShort_t        EmcTrigger_mId[kMaxEmcTrigger];   //[EmcTrigger_]
  UShort_t        EmcTrigger_mAdc[kMaxEmcTrigger];   //[EmcTrigger_]
  vector<unsigned short> EmcTrigger_mSmdE[kMaxEmcTrigger];
  vector<unsigned short> EmcTrigger_mSmdP[kMaxEmcTrigger];
  Int_t           MtdTrigger_;
  UShort_t        MtdTrigger_mVpdTacSum[kMaxMtdTrigger];   //[MtdTrigger_]
  UInt_t          MtdTrigger_mTHUBtime[kMaxMtdTrigger][2];   //[MtdTrigger_]
  UShort_t        MtdTrigger_mQTtacSum[kMaxMtdTrigger][8][8];   //[MtdTrigger_]
  UShort_t        MtdTrigger_mMT101Tac[kMaxMtdTrigger][8][2];   //[MtdTrigger_]
  UChar_t         MtdTrigger_mMT101Id[kMaxMtdTrigger][8][2];   //[MtdTrigger_]
  UInt_t          MtdTrigger_mTF201TriggerBit[kMaxMtdTrigger];   //[MtdTrigger_]
  Char_t          MtdTrigger_mShouldHaveRejectEvent[kMaxMtdTrigger];   //[MtdTrigger_]
  Int_t           BTowHit_;
  UShort_t        BTowHit_mAdc[kMaxBTowHit];   //[BTowHit_]
  Short_t         BTowHit_mE[kMaxBTowHit];   //[BTowHit_]
  Int_t           BTofHit_;
  Short_t         BTofHit_mId[kMaxBTofHit];   //[BTofHit_]
  Int_t           MtdHit_;
  Short_t         MtdHit_mgChannel[kMaxMtdHit];   //[MtdHit_]
  UChar_t         MtdHit_mTriggerFlag[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mLeadingEdgeTime_first[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mLeadingEdgeTime_second[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mTrailingEdgeTime_first[kMaxMtdHit];   //[MtdHit_]
  Float_t         MtdHit_mTrailingEdgeTime_second[kMaxMtdHit];   //[MtdHit_]
  Int_t           BbcHit_;
  Short_t         BbcHit_mId[kMaxBbcHit];   //[BbcHit_]
  Int_t           BbcHit_mQTdata[kMaxBbcHit];   //[BbcHit_]
  Int_t           EpdHit_;
  Short_t         EpdHit_mId[kMaxEpdHit];   //[EpdHit_]
  Int_t           EpdHit_mQTdata[kMaxEpdHit];   //[EpdHit_]
  Float_t         EpdHit_mnMIP[kMaxEpdHit];   //[EpdHit_]
  Int_t           FmsHit_;
  UShort_t        FmsHit_mChannelDetectorId[kMaxFmsHit];   //[FmsHit_]
  UShort_t        FmsHit_mAdc[kMaxFmsHit];   //[FmsHit_]
  Int_t           EmcPidTraits_;
  Short_t         EmcPidTraits_mTrackIndex[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBemcId[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBemcAdc0[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBemcE0[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBemcE[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBemcZDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBemcPhiDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
  UChar_t         EmcPidTraits_mBemcSmdNEta[kMaxEmcPidTraits];   //[EmcPidTraits_]
  UChar_t         EmcPidTraits_mBemcSmdNPhi[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBtowId[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Char_t          EmcPidTraits_mBtowId23[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBtowE[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBtowE2[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBtowE3[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBtowEtaDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Short_t         EmcPidTraits_mBtowPhiDist[kMaxEmcPidTraits];   //[EmcPidTraits_]
  Int_t           BTofPidTraits_;
  Short_t         BTofPidTraits_mTrackIndex[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofCellId[kMaxBTofPidTraits];   //[BTofPidTraits_]
  UChar_t         BTofPidTraits_mBTofMatchFlag[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Float_t         BTofPidTraits_mBTof[kMaxBTofPidTraits];   //[BTofPidTraits_]
  UShort_t        BTofPidTraits_mBTofBeta[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofYLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofZLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofHitPosX[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofHitPosY[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mBTofHitPosZ[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mNSigmaElectron[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mNSigmaPion[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mNSigmaKaon[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Short_t         BTofPidTraits_mNSigmaProton[kMaxBTofPidTraits];   //[BTofPidTraits_]
  Int_t           MtdPidTraits_;
  Short_t         MtdPidTraits_mTrackIndex[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Short_t         MtdPidTraits_mMtdHitIndex[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Char_t          MtdPidTraits_mMatchFlag[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Short_t         MtdPidTraits_mDeltaY[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Short_t         MtdPidTraits_mDeltaZ[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Float_t         MtdPidTraits_mDeltaTimeOfFlight[kMaxMtdPidTraits];   //[MtdPidTraits_]
  UShort_t        MtdPidTraits_mBeta[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Short_t         MtdPidTraits_mMtdHitChan[kMaxMtdPidTraits];   //[MtdPidTraits_]
  Int_t           TrackCovMatrix_;
  Float16_t       TrackCovMatrix_mImp[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
  Float16_t       TrackCovMatrix_mZ[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
  Float16_t       TrackCovMatrix_mPsi[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
  Float16_t       TrackCovMatrix_mPti[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
  Float16_t       TrackCovMatrix_mTan[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
  Float16_t       TrackCovMatrix_mCurv[kMaxTrackCovMatrix];   //[TrackCovMatrix_]
  Float16_t       TrackCovMatrix_mSigma[kMaxTrackCovMatrix][5];   //[TrackCovMatrix_]
  Float16_t       TrackCovMatrix_mCorr[kMaxTrackCovMatrix][10];   //[TrackCovMatrix_]
  Int_t           BEmcSmdEHit_;
  Short_t         BEmcSmdEHit_mId[kMaxBEmcSmdEHit];   //[BEmcSmdEHit_]
  Short_t         BEmcSmdEHit_mAdc[kMaxBEmcSmdEHit];   //[BEmcSmdEHit_]
  Float_t         BEmcSmdEHit_mEnergy[kMaxBEmcSmdEHit];   //[BEmcSmdEHit_]
  Int_t           BEmcSmdPHit_;
  Short_t         BEmcSmdPHit_mId[kMaxBEmcSmdPHit];   //[BEmcSmdPHit_]
  Short_t         BEmcSmdPHit_mAdc[kMaxBEmcSmdPHit];   //[BEmcSmdPHit_]
  Float_t         BEmcSmdPHit_mEnergy[kMaxBEmcSmdPHit];   //[BEmcSmdPHit_]
  Int_t           ETofHit_;
  UChar_t         ETofHit_mGeomId[kMaxETofHit];   //[ETofHit_]
  Short_t         ETofHit_mLocalX[kMaxETofHit];   //[ETofHit_]
  Short_t         ETofHit_mLocalY[kMaxETofHit];   //[ETofHit_]
  UChar_t         ETofHit_mClusterSize[kMaxETofHit];   //[ETofHit_]
  Float_t         ETofHit_mLeadingEdgeTime[kMaxETofHit];   //[ETofHit_]
  UShort_t        ETofHit_mTimeOverThreshold[kMaxETofHit];   //[ETofHit_]
  Int_t           ETofPidTraits_;
  Short_t         ETofPidTraits_mTrackIndex[kMaxETofPidTraits];   //[ETofPidTraits_]
  Short_t         ETofPidTraits_mHitIndex[kMaxETofPidTraits];   //[ETofPidTraits_]
  Char_t          ETofPidTraits_mMatchFlag[kMaxETofPidTraits];   //[ETofPidTraits_]
  Float_t         ETofPidTraits_mTimeOfFlight[kMaxETofPidTraits];   //[ETofPidTraits_]
  UShort_t        ETofPidTraits_mBeta[kMaxETofPidTraits];   //[ETofPidTraits_]
  Short_t         ETofPidTraits_mDeltaX[kMaxETofPidTraits];   //[ETofPidTraits_]
  Short_t         ETofPidTraits_mDeltaY[kMaxETofPidTraits];   //[ETofPidTraits_]
  Short_t         ETofPidTraits_mCrossingX[kMaxETofPidTraits];   //[ETofPidTraits_]
  Short_t         ETofPidTraits_mCrossingY[kMaxETofPidTraits];   //[ETofPidTraits_]
  Short_t         ETofPidTraits_mCrossingZ[kMaxETofPidTraits];   //[ETofPidTraits_]
  Int_t           McVertex_;
  Int_t           McVertex_mId[kMaxMcVertex];   //[McVertex_]
  UShort_t        McVertex_mNoDaughters[kMaxMcVertex];   //[McVertex_]
  Int_t           McVertex_mIdParTrk[kMaxMcVertex];   //[McVertex_]
  Int_t           McVertex_mIsInterm[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mTime[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mVx[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mVy[kMaxMcVertex];   //[McVertex_]
  Float_t         McVertex_mVz[kMaxMcVertex];   //[McVertex_]
  Int_t           McTrack_;
  UShort_t        McTrack_mId[kMaxMcTrack];   //[McTrack_]
  Int_t           McTrack_mGePid[kMaxMcTrack];   //[McTrack_]
  Char_t          McTrack_mCharge[kMaxMcTrack];   //[McTrack_]
  UChar_t         McTrack_mHits[kMaxMcTrack][22];   //[McTrack_]
  Float_t         McTrack_mPx[kMaxMcTrack];   //[McTrack_]
  Float_t         McTrack_mPy[kMaxMcTrack];   //[McTrack_]
  Float_t         McTrack_mPz[kMaxMcTrack];   //[McTrack_]
  Float_t         McTrack_mE[kMaxMcTrack];   //[McTrack_]
  Bool_t          McTrack_mIsFromShower[kMaxMcTrack];   //[McTrack_]
  Short_t         McTrack_mIdVtxStart[kMaxMcTrack];   //[McTrack_]
  Short_t         McTrack_mIdVtxStop[kMaxMcTrack];   //[McTrack_]
  Short_t         McTrack_mIdVtxItrmd[kMaxMcTrack];   //[McTrack_]

  // List of branches
  TBranch        *b_Event_;   //!
  TBranch        *b_Event_mRunId;   //!
  TBranch        *b_Event_mEventId;   //!
  TBranch        *b_Event_mFillId;   //!
  TBranch        *b_Event_mBField;   //!
  TBranch        *b_Event_mTime;   //!
  TBranch        *b_Event_mPrimaryVertexX;   //!
  TBranch        *b_Event_mPrimaryVertexY;   //!
  TBranch        *b_Event_mPrimaryVertexZ;   //!
  TBranch        *b_Event_mPrimaryVertexErrorX;   //!
  TBranch        *b_Event_mPrimaryVertexErrorY;   //!
  TBranch        *b_Event_mPrimaryVertexErrorZ;   //!
  TBranch        *b_Event_mRanking;   //!
  TBranch        *b_Event_mNBEMCMatch;   //!
  TBranch        *b_Event_mNBTOFMatch;   //!
  TBranch        *b_Event_mTriggerIds;   //!
  TBranch        *b_Event_mRefMultFtpcEast;   //!
  TBranch        *b_Event_mRefMultFtpcWest;   //!
  TBranch        *b_Event_mRefMultNeg;   //!
  TBranch        *b_Event_mRefMultPos;   //!
  TBranch        *b_Event_mRefMult2NegEast;   //!
  TBranch        *b_Event_mRefMult2PosEast;   //!
  TBranch        *b_Event_mRefMult2NegWest;   //!
  TBranch        *b_Event_mRefMult2PosWest;   //!
  TBranch        *b_Event_mRefMult3NegEast;   //!
  TBranch        *b_Event_mRefMult3PosEast;   //!
  TBranch        *b_Event_mRefMult3NegWest;   //!
  TBranch        *b_Event_mRefMult3PosWest;   //!
  TBranch        *b_Event_mRefMult4NegEast;   //!
  TBranch        *b_Event_mRefMult4PosEast;   //!
  TBranch        *b_Event_mRefMult4NegWest;   //!
  TBranch        *b_Event_mRefMult4PosWest;   //!
  TBranch        *b_Event_mRefMultHalfNegEast;   //!
  TBranch        *b_Event_mRefMultHalfPosEast;   //!
  TBranch        *b_Event_mRefMultHalfNegWest;   //!
  TBranch        *b_Event_mRefMultHalfPosWest;   //!
  TBranch        *b_Event_mGRefMult;   //!
  TBranch        *b_Event_mNumberOfGlobalTracks;   //!
  TBranch        *b_Event_mbTofTrayMultiplicity;   //!
  TBranch        *b_Event_mNHitsHFT;   //!
  TBranch        *b_Event_mNVpdHitsEast;   //!
  TBranch        *b_Event_mNVpdHitsWest;   //!
  TBranch        *b_Event_mNTofT0;   //!
  TBranch        *b_Event_mVzVpd;   //!
  TBranch        *b_Event_mZDCx;   //!
  TBranch        *b_Event_mBBCx;   //!
  TBranch        *b_Event_mBackgroundRate;   //!
  TBranch        *b_Event_mBbcBlueBackgroundRate;   //!
  TBranch        *b_Event_mBbcYellowBackgroundRate;   //!
  TBranch        *b_Event_mBbcEastRate;   //!
  TBranch        *b_Event_mBbcWestRate;   //!
  TBranch        *b_Event_mZdcEastRate;   //!
  TBranch        *b_Event_mZdcWestRate;   //!
  TBranch        *b_Event_mZdcSumAdcEast;   //!
  TBranch        *b_Event_mZdcSumAdcWest;   //!
  TBranch        *b_Event_mZdcSmdEastHorizontal;   //!
  TBranch        *b_Event_mZdcSmdEastVertical;   //!
  TBranch        *b_Event_mZdcSmdWestHorizontal;   //!
  TBranch        *b_Event_mZdcSmdWestVertical;   //!
  TBranch        *b_Event_mBbcAdcEast;   //!
  TBranch        *b_Event_mBbcAdcWest;   //!
  TBranch        *b_Event_mHighTowerThreshold;   //!
  TBranch        *b_Event_mJetPatchThreshold;   //!
  TBranch        *b_Event_mBunchCrossId;   //!
  TBranch        *b_Event_mETofHitMultiplicity;   //!
  TBranch        *b_Event_mETofDigiMultiplicity;   //!
  TBranch        *b_Event_mNumberOfPrimaryTracks;   //!
  TBranch        *b_Event_mZdcUnAttenuated;   //!
  TBranch        *b_Track_;   //!
  TBranch        *b_Track_mId;   //!
  TBranch        *b_Track_mChi2;   //!
  TBranch        *b_Track_mPMomentumX;   //!
  TBranch        *b_Track_mPMomentumY;   //!
  TBranch        *b_Track_mPMomentumZ;   //!
  TBranch        *b_Track_mGMomentumX;   //!
  TBranch        *b_Track_mGMomentumY;   //!
  TBranch        *b_Track_mGMomentumZ;   //!
  TBranch        *b_Track_mOriginX;   //!
  TBranch        *b_Track_mOriginY;   //!
  TBranch        *b_Track_mOriginZ;   //!
  TBranch        *b_Track_mDedx;   //!
  TBranch        *b_Track_mDedxError;   //!
  TBranch        *b_Track_mNHitsFit;   //!
  TBranch        *b_Track_mNHitsMax;   //!
  TBranch        *b_Track_mNHitsDedx;   //!
  TBranch        *b_Track_mNSigmaPion;   //!
  TBranch        *b_Track_mNSigmaKaon;   //!
  TBranch        *b_Track_mNSigmaProton;   //!
  TBranch        *b_Track_mNSigmaElectron;   //!
  TBranch        *b_Track_mTopologyMap;   //!
  TBranch        *b_Track_mBEmcPidTraitsIndex;   //!
  TBranch        *b_Track_mBTofPidTraitsIndex;   //!
  TBranch        *b_Track_mMtdPidTraitsIndex;   //!
  TBranch        *b_Track_mETofPidTraitsIndex;   //!
  TBranch        *b_Track_mBEmcMatchedTowerIndex;   //!
  TBranch        *b_Track_mTopoMap_iTpc;   //!
  TBranch        *b_Track_mIdTruth;   //!
  TBranch        *b_Track_mQATruth;   //!
  TBranch        *b_Track_mVertexIndex;   //!
  TBranch        *b_EmcTrigger_;   //!
  TBranch        *b_EmcTrigger_mFlag;   //!
  TBranch        *b_EmcTrigger_mId;   //!
  TBranch        *b_EmcTrigger_mAdc;   //!
  TBranch        *b_EmcTrigger_mSmdE;   //!
  TBranch        *b_EmcTrigger_mSmdP;   //!
  TBranch        *b_MtdTrigger_;   //!
  TBranch        *b_MtdTrigger_mVpdTacSum;   //!
  TBranch        *b_MtdTrigger_mTHUBtime;   //!
  TBranch        *b_MtdTrigger_mQTtacSum;   //!
  TBranch        *b_MtdTrigger_mMT101Tac;   //!
  TBranch        *b_MtdTrigger_mMT101Id;   //!
  TBranch        *b_MtdTrigger_mTF201TriggerBit;   //!
  TBranch        *b_MtdTrigger_mShouldHaveRejectEvent;   //!
  TBranch        *b_BTowHit_;   //!
  TBranch        *b_BTowHit_mAdc;   //!
  TBranch        *b_BTowHit_mE;   //!
  TBranch        *b_BTofHit_;   //!
  TBranch        *b_BTofHit_mId;   //!
  TBranch        *b_MtdHit_;   //!
  TBranch        *b_MtdHit_mgChannel;   //!
  TBranch        *b_MtdHit_mTriggerFlag;   //!
  TBranch        *b_MtdHit_mLeadingEdgeTime_first;   //!
  TBranch        *b_MtdHit_mLeadingEdgeTime_second;   //!
  TBranch        *b_MtdHit_mTrailingEdgeTime_first;   //!
  TBranch        *b_MtdHit_mTrailingEdgeTime_second;   //!
  TBranch        *b_BbcHit_;   //!
  TBranch        *b_BbcHit_mId;   //!
  TBranch        *b_BbcHit_mQTdata;   //!
  TBranch        *b_EpdHit_;   //!
  TBranch        *b_EpdHit_mId;   //!
  TBranch        *b_EpdHit_mQTdata;   //!
  TBranch        *b_EpdHit_mnMIP;   //!
  TBranch        *b_FmsHit_;   //!
  TBranch        *b_FmsHit_mChannelDetectorId;   //!
  TBranch        *b_FmsHit_mAdc;   //!
  TBranch        *b_EmcPidTraits_;   //!
  TBranch        *b_EmcPidTraits_mTrackIndex;   //!
  TBranch        *b_EmcPidTraits_mBemcId;   //!
  TBranch        *b_EmcPidTraits_mBemcAdc0;   //!
  TBranch        *b_EmcPidTraits_mBemcE0;   //!
  TBranch        *b_EmcPidTraits_mBemcE;   //!
  TBranch        *b_EmcPidTraits_mBemcZDist;   //!
  TBranch        *b_EmcPidTraits_mBemcPhiDist;   //!
  TBranch        *b_EmcPidTraits_mBemcSmdNEta;   //!
  TBranch        *b_EmcPidTraits_mBemcSmdNPhi;   //!
  TBranch        *b_EmcPidTraits_mBtowId;   //!
  TBranch        *b_EmcPidTraits_mBtowId23;   //!
  TBranch        *b_EmcPidTraits_mBtowE;   //!
  TBranch        *b_EmcPidTraits_mBtowE2;   //!
  TBranch        *b_EmcPidTraits_mBtowE3;   //!
  TBranch        *b_EmcPidTraits_mBtowEtaDist;   //!
  TBranch        *b_EmcPidTraits_mBtowPhiDist;   //!
  TBranch        *b_BTofPidTraits_;   //!
  TBranch        *b_BTofPidTraits_mTrackIndex;   //!
  TBranch        *b_BTofPidTraits_mBTofCellId;   //!
  TBranch        *b_BTofPidTraits_mBTofMatchFlag;   //!
  TBranch        *b_BTofPidTraits_mBTof;   //!
  TBranch        *b_BTofPidTraits_mBTofBeta;   //!
  TBranch        *b_BTofPidTraits_mBTofYLocal;   //!
  TBranch        *b_BTofPidTraits_mBTofZLocal;   //!
  TBranch        *b_BTofPidTraits_mBTofHitPosX;   //!
  TBranch        *b_BTofPidTraits_mBTofHitPosY;   //!
  TBranch        *b_BTofPidTraits_mBTofHitPosZ;   //!
  TBranch        *b_BTofPidTraits_mNSigmaElectron;   //!
  TBranch        *b_BTofPidTraits_mNSigmaPion;   //!
  TBranch        *b_BTofPidTraits_mNSigmaKaon;   //!
  TBranch        *b_BTofPidTraits_mNSigmaProton;   //!
  TBranch        *b_MtdPidTraits_;   //!
  TBranch        *b_MtdPidTraits_mTrackIndex;   //!
  TBranch        *b_MtdPidTraits_mMtdHitIndex;   //!
  TBranch        *b_MtdPidTraits_mMatchFlag;   //!
  TBranch        *b_MtdPidTraits_mDeltaY;   //!
  TBranch        *b_MtdPidTraits_mDeltaZ;   //!
  TBranch        *b_MtdPidTraits_mDeltaTimeOfFlight;   //!
  TBranch        *b_MtdPidTraits_mBeta;   //!
  TBranch        *b_MtdPidTraits_mMtdHitChan;   //!
  TBranch        *b_TrackCovMatrix_;   //!
  TBranch        *b_TrackCovMatrix_mImp;   //!
  TBranch        *b_TrackCovMatrix_mZ;   //!
  TBranch        *b_TrackCovMatrix_mPsi;   //!
  TBranch        *b_TrackCovMatrix_mPti;   //!
  TBranch        *b_TrackCovMatrix_mTan;   //!
  TBranch        *b_TrackCovMatrix_mCurv;   //!
  TBranch        *b_TrackCovMatrix_mSigma;   //!
  TBranch        *b_TrackCovMatrix_mCorr;   //!
  TBranch        *b_BEmcSmdEHit_;   //!
  TBranch        *b_BEmcSmdEHit_mId;   //!
  TBranch        *b_BEmcSmdEHit_mAdc;   //!
  TBranch        *b_BEmcSmdEHit_mEnergy;   //!
  TBranch        *b_BEmcSmdPHit_;   //!
  TBranch        *b_BEmcSmdPHit_mId;   //!
  TBranch        *b_BEmcSmdPHit_mAdc;   //!
  TBranch        *b_BEmcSmdPHit_mEnergy;   //!
  TBranch        *b_ETofHit_;   //!
  TBranch        *b_ETofHit_mGeomId;   //!
  TBranch        *b_ETofHit_mLocalX;   //!
  TBranch        *b_ETofHit_mLocalY;   //!
  TBranch        *b_ETofHit_mClusterSize;   //!
  TBranch        *b_ETofHit_mLeadingEdgeTime;   //!
  TBranch        *b_ETofHit_mTimeOverThreshold;   //!
  TBranch        *b_ETofPidTraits_;   //!
  TBranch        *b_ETofPidTraits_mTrackIndex;   //!
  TBranch        *b_ETofPidTraits_mHitIndex;   //!
  TBranch        *b_ETofPidTraits_mMatchFlag;   //!
  TBranch        *b_ETofPidTraits_mTimeOfFlight;   //!
  TBranch        *b_ETofPidTraits_mBeta;   //!
  TBranch        *b_ETofPidTraits_mDeltaX;   //!
  TBranch        *b_ETofPidTraits_mDeltaY;   //!
  TBranch        *b_ETofPidTraits_mCrossingX;   //!
  TBranch        *b_ETofPidTraits_mCrossingY;   //!
  TBranch        *b_ETofPidTraits_mCrossingZ;   //!
  TBranch        *b_McVertex_;   //!
  TBranch        *b_McVertex_mId;   //!
  TBranch        *b_McVertex_mNoDaughters;   //!
  TBranch        *b_McVertex_mIdParTrk;   //!
  TBranch        *b_McVertex_mIsInterm;   //!
  TBranch        *b_McVertex_mTime;   //!
  TBranch        *b_McVertex_mVx;   //!
  TBranch        *b_McVertex_mVy;   //!
  TBranch        *b_McVertex_mVz;   //!
  TBranch        *b_McTrack_;   //!
  TBranch        *b_McTrack_mId;   //!
  TBranch        *b_McTrack_mGePid;   //!
  TBranch        *b_McTrack_mCharge;   //!
  TBranch        *b_McTrack_mHits;   //!
  TBranch        *b_McTrack_mPx;   //!
  TBranch        *b_McTrack_mPy;   //!
  TBranch        *b_McTrack_mPz;   //!
  TBranch        *b_McTrack_mE;   //!
  TBranch        *b_McTrack_mIsFromShower;   //!
  TBranch        *b_McTrack_mIdVtxStart;   //!
  TBranch        *b_McTrack_mIdVtxStop;   //!
  TBranch        *b_McTrack_mIdVtxItrmd;   //!

  // fMCPico->SetBranchStatus("*", 0);

  fMCPico->SetBranchAddress("Event", &Event_, &b_Event_);
  fMCPico->SetBranchAddress("Event.mRunId", Event_mRunId, &b_Event_mRunId);
  fMCPico->SetBranchAddress("Event.mEventId", Event_mEventId, &b_Event_mEventId);
  fMCPico->SetBranchAddress("Event.mFillId", Event_mFillId, &b_Event_mFillId);
  fMCPico->SetBranchAddress("Event.mBField", Event_mBField, &b_Event_mBField);
  fMCPico->SetBranchAddress("Event.mTime", Event_mTime, &b_Event_mTime);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexX", Event_mPrimaryVertexX, &b_Event_mPrimaryVertexX);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexY", Event_mPrimaryVertexY, &b_Event_mPrimaryVertexY);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexZ", Event_mPrimaryVertexZ, &b_Event_mPrimaryVertexZ);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexErrorX", Event_mPrimaryVertexErrorX, &b_Event_mPrimaryVertexErrorX);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexErrorY", Event_mPrimaryVertexErrorY, &b_Event_mPrimaryVertexErrorY);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexErrorZ", Event_mPrimaryVertexErrorZ, &b_Event_mPrimaryVertexErrorZ);
  fMCPico->SetBranchAddress("Event.mRanking", Event_mRanking, &b_Event_mRanking);
  fMCPico->SetBranchAddress("Event.mNBEMCMatch", Event_mNBEMCMatch, &b_Event_mNBEMCMatch);
  fMCPico->SetBranchAddress("Event.mNBTOFMatch", Event_mNBTOFMatch, &b_Event_mNBTOFMatch);
  fMCPico->SetBranchAddress("Event.mTriggerIds", Event_mTriggerIds, &b_Event_mTriggerIds);
  fMCPico->SetBranchAddress("Event.mRefMultFtpcEast", Event_mRefMultFtpcEast, &b_Event_mRefMultFtpcEast);
  fMCPico->SetBranchAddress("Event.mRefMultFtpcWest", Event_mRefMultFtpcWest, &b_Event_mRefMultFtpcWest);
  fMCPico->SetBranchAddress("Event.mRefMultNeg", Event_mRefMultNeg, &b_Event_mRefMultNeg);
  fMCPico->SetBranchAddress("Event.mRefMultPos", Event_mRefMultPos, &b_Event_mRefMultPos);
  fMCPico->SetBranchAddress("Event.mRefMult2NegEast", Event_mRefMult2NegEast, &b_Event_mRefMult2NegEast);
  fMCPico->SetBranchAddress("Event.mRefMult2PosEast", Event_mRefMult2PosEast, &b_Event_mRefMult2PosEast);
  fMCPico->SetBranchAddress("Event.mRefMult2NegWest", Event_mRefMult2NegWest, &b_Event_mRefMult2NegWest);
  fMCPico->SetBranchAddress("Event.mRefMult2PosWest", Event_mRefMult2PosWest, &b_Event_mRefMult2PosWest);
  fMCPico->SetBranchAddress("Event.mRefMult3NegEast", Event_mRefMult3NegEast, &b_Event_mRefMult3NegEast);
  fMCPico->SetBranchAddress("Event.mRefMult3PosEast", Event_mRefMult3PosEast, &b_Event_mRefMult3PosEast);
  fMCPico->SetBranchAddress("Event.mRefMult3NegWest", Event_mRefMult3NegWest, &b_Event_mRefMult3NegWest);
  fMCPico->SetBranchAddress("Event.mRefMult3PosWest", Event_mRefMult3PosWest, &b_Event_mRefMult3PosWest);
  fMCPico->SetBranchAddress("Event.mRefMult4NegEast", Event_mRefMult4NegEast, &b_Event_mRefMult4NegEast);
  fMCPico->SetBranchAddress("Event.mRefMult4PosEast", Event_mRefMult4PosEast, &b_Event_mRefMult4PosEast);
  fMCPico->SetBranchAddress("Event.mRefMult4NegWest", Event_mRefMult4NegWest, &b_Event_mRefMult4NegWest);
  fMCPico->SetBranchAddress("Event.mRefMult4PosWest", Event_mRefMult4PosWest, &b_Event_mRefMult4PosWest);
  fMCPico->SetBranchAddress("Event.mRefMultHalfNegEast", Event_mRefMultHalfNegEast, &b_Event_mRefMultHalfNegEast);
  fMCPico->SetBranchAddress("Event.mRefMultHalfPosEast", Event_mRefMultHalfPosEast, &b_Event_mRefMultHalfPosEast);
  fMCPico->SetBranchAddress("Event.mRefMultHalfNegWest", Event_mRefMultHalfNegWest, &b_Event_mRefMultHalfNegWest);
  fMCPico->SetBranchAddress("Event.mRefMultHalfPosWest", Event_mRefMultHalfPosWest, &b_Event_mRefMultHalfPosWest);
  fMCPico->SetBranchAddress("Event.mGRefMult", Event_mGRefMult, &b_Event_mGRefMult);
  fMCPico->SetBranchAddress("Event.mNumberOfGlobalTracks", Event_mNumberOfGlobalTracks, &b_Event_mNumberOfGlobalTracks);
  fMCPico->SetBranchAddress("Event.mbTofTrayMultiplicity", Event_mbTofTrayMultiplicity, &b_Event_mbTofTrayMultiplicity);
  fMCPico->SetBranchAddress("Event.mNHitsHFT[4]", Event_mNHitsHFT, &b_Event_mNHitsHFT);
  fMCPico->SetBranchAddress("Event.mNVpdHitsEast", Event_mNVpdHitsEast, &b_Event_mNVpdHitsEast);
  fMCPico->SetBranchAddress("Event.mNVpdHitsWest", Event_mNVpdHitsWest, &b_Event_mNVpdHitsWest);
  fMCPico->SetBranchAddress("Event.mNTofT0", Event_mNTofT0, &b_Event_mNTofT0);
  fMCPico->SetBranchAddress("Event.mVzVpd", Event_mVzVpd, &b_Event_mVzVpd);
  fMCPico->SetBranchAddress("Event.mZDCx", Event_mZDCx, &b_Event_mZDCx);
  fMCPico->SetBranchAddress("Event.mBBCx", Event_mBBCx, &b_Event_mBBCx);
  fMCPico->SetBranchAddress("Event.mBackgroundRate", Event_mBackgroundRate, &b_Event_mBackgroundRate);
  fMCPico->SetBranchAddress("Event.mBbcBlueBackgroundRate", Event_mBbcBlueBackgroundRate, &b_Event_mBbcBlueBackgroundRate);
  fMCPico->SetBranchAddress("Event.mBbcYellowBackgroundRate", Event_mBbcYellowBackgroundRate, &b_Event_mBbcYellowBackgroundRate);
  fMCPico->SetBranchAddress("Event.mBbcEastRate", Event_mBbcEastRate, &b_Event_mBbcEastRate);
  fMCPico->SetBranchAddress("Event.mBbcWestRate", Event_mBbcWestRate, &b_Event_mBbcWestRate);
  fMCPico->SetBranchAddress("Event.mZdcEastRate", Event_mZdcEastRate, &b_Event_mZdcEastRate);
  fMCPico->SetBranchAddress("Event.mZdcWestRate", Event_mZdcWestRate, &b_Event_mZdcWestRate);
  fMCPico->SetBranchAddress("Event.mZdcSumAdcEast", Event_mZdcSumAdcEast, &b_Event_mZdcSumAdcEast);
  fMCPico->SetBranchAddress("Event.mZdcSumAdcWest", Event_mZdcSumAdcWest, &b_Event_mZdcSumAdcWest);
  fMCPico->SetBranchAddress("Event.mZdcSmdEastHorizontal[8]", Event_mZdcSmdEastHorizontal, &b_Event_mZdcSmdEastHorizontal);
  fMCPico->SetBranchAddress("Event.mZdcSmdEastVertical[8]", Event_mZdcSmdEastVertical, &b_Event_mZdcSmdEastVertical);
  fMCPico->SetBranchAddress("Event.mZdcSmdWestHorizontal[8]", Event_mZdcSmdWestHorizontal, &b_Event_mZdcSmdWestHorizontal);
  fMCPico->SetBranchAddress("Event.mZdcSmdWestVertical[8]", Event_mZdcSmdWestVertical, &b_Event_mZdcSmdWestVertical);
  fMCPico->SetBranchAddress("Event.mBbcAdcEast[24]", Event_mBbcAdcEast, &b_Event_mBbcAdcEast);
  fMCPico->SetBranchAddress("Event.mBbcAdcWest[24]", Event_mBbcAdcWest, &b_Event_mBbcAdcWest);
  fMCPico->SetBranchAddress("Event.mHighTowerThreshold[4]", Event_mHighTowerThreshold, &b_Event_mHighTowerThreshold);
  fMCPico->SetBranchAddress("Event.mJetPatchThreshold[4]", Event_mJetPatchThreshold, &b_Event_mJetPatchThreshold);
  fMCPico->SetBranchAddress("Event.mBunchCrossId", Event_mBunchCrossId, &b_Event_mBunchCrossId);
  fMCPico->SetBranchAddress("Event.mETofHitMultiplicity", Event_mETofHitMultiplicity, &b_Event_mETofHitMultiplicity);
  fMCPico->SetBranchAddress("Event.mETofDigiMultiplicity", Event_mETofDigiMultiplicity, &b_Event_mETofDigiMultiplicity);
  fMCPico->SetBranchAddress("Event.mNumberOfPrimaryTracks", Event_mNumberOfPrimaryTracks, &b_Event_mNumberOfPrimaryTracks);
  fMCPico->SetBranchAddress("Event.mZdcUnAttenuated[2]", Event_mZdcUnAttenuated, &b_Event_mZdcUnAttenuated);
  fMCPico->SetBranchAddress("Track", &Track_, &b_Track_);
  fMCPico->SetBranchAddress("Track.mId", Track_mId, &b_Track_mId);
  fMCPico->SetBranchAddress("Track.mChi2", Track_mChi2, &b_Track_mChi2);
  fMCPico->SetBranchAddress("Track.mPMomentumX", Track_mPMomentumX, &b_Track_mPMomentumX);
  fMCPico->SetBranchAddress("Track.mPMomentumY", Track_mPMomentumY, &b_Track_mPMomentumY);
  fMCPico->SetBranchAddress("Track.mPMomentumZ", Track_mPMomentumZ, &b_Track_mPMomentumZ);
  fMCPico->SetBranchAddress("Track.mGMomentumX", Track_mGMomentumX, &b_Track_mGMomentumX);
  fMCPico->SetBranchAddress("Track.mGMomentumY", Track_mGMomentumY, &b_Track_mGMomentumY);
  fMCPico->SetBranchAddress("Track.mGMomentumZ", Track_mGMomentumZ, &b_Track_mGMomentumZ);
  fMCPico->SetBranchAddress("Track.mOriginX", Track_mOriginX, &b_Track_mOriginX);
  fMCPico->SetBranchAddress("Track.mOriginY", Track_mOriginY, &b_Track_mOriginY);
  fMCPico->SetBranchAddress("Track.mOriginZ", Track_mOriginZ, &b_Track_mOriginZ);
  fMCPico->SetBranchAddress("Track.mDedx", Track_mDedx, &b_Track_mDedx);
  fMCPico->SetBranchAddress("Track.mDedxError", Track_mDedxError, &b_Track_mDedxError);
  fMCPico->SetBranchAddress("Track.mNHitsFit", Track_mNHitsFit, &b_Track_mNHitsFit);
  fMCPico->SetBranchAddress("Track.mNHitsMax", Track_mNHitsMax, &b_Track_mNHitsMax);
  fMCPico->SetBranchAddress("Track.mNHitsDedx", Track_mNHitsDedx, &b_Track_mNHitsDedx);
  fMCPico->SetBranchAddress("Track.mNSigmaPion", Track_mNSigmaPion, &b_Track_mNSigmaPion);
  fMCPico->SetBranchAddress("Track.mNSigmaKaon", Track_mNSigmaKaon, &b_Track_mNSigmaKaon);
  fMCPico->SetBranchAddress("Track.mNSigmaProton", Track_mNSigmaProton, &b_Track_mNSigmaProton);
  fMCPico->SetBranchAddress("Track.mNSigmaElectron", Track_mNSigmaElectron, &b_Track_mNSigmaElectron);
  fMCPico->SetBranchAddress("Track.mTopologyMap[2]", Track_mTopologyMap, &b_Track_mTopologyMap);
  fMCPico->SetBranchAddress("Track.mBEmcPidTraitsIndex", Track_mBEmcPidTraitsIndex, &b_Track_mBEmcPidTraitsIndex);
  fMCPico->SetBranchAddress("Track.mBTofPidTraitsIndex", Track_mBTofPidTraitsIndex, &b_Track_mBTofPidTraitsIndex);
  fMCPico->SetBranchAddress("Track.mMtdPidTraitsIndex", Track_mMtdPidTraitsIndex, &b_Track_mMtdPidTraitsIndex);
  fMCPico->SetBranchAddress("Track.mETofPidTraitsIndex", Track_mETofPidTraitsIndex, &b_Track_mETofPidTraitsIndex);
  fMCPico->SetBranchAddress("Track.mBEmcMatchedTowerIndex", Track_mBEmcMatchedTowerIndex, &b_Track_mBEmcMatchedTowerIndex);
  fMCPico->SetBranchAddress("Track.mTopoMap_iTpc", Track_mTopoMap_iTpc, &b_Track_mTopoMap_iTpc);
  fMCPico->SetBranchAddress("Track.mIdTruth", Track_mIdTruth, &b_Track_mIdTruth);
  fMCPico->SetBranchAddress("Track.mQATruth", Track_mQATruth, &b_Track_mQATruth);
  fMCPico->SetBranchAddress("Track.mVertexIndex", Track_mVertexIndex, &b_Track_mVertexIndex);
  fMCPico->SetBranchAddress("EmcTrigger", &EmcTrigger_, &b_EmcTrigger_);
  fMCPico->SetBranchAddress("EmcTrigger.mFlag", EmcTrigger_mFlag, &b_EmcTrigger_mFlag);
  fMCPico->SetBranchAddress("EmcTrigger.mId", EmcTrigger_mId, &b_EmcTrigger_mId);
  fMCPico->SetBranchAddress("EmcTrigger.mAdc", EmcTrigger_mAdc, &b_EmcTrigger_mAdc);
  fMCPico->SetBranchAddress("EmcTrigger.mSmdE", EmcTrigger_mSmdE, &b_EmcTrigger_mSmdE);
  fMCPico->SetBranchAddress("EmcTrigger.mSmdP", EmcTrigger_mSmdP, &b_EmcTrigger_mSmdP);
  fMCPico->SetBranchAddress("MtdTrigger", &MtdTrigger_, &b_MtdTrigger_);
  fMCPico->SetBranchAddress("MtdTrigger.mVpdTacSum", MtdTrigger_mVpdTacSum, &b_MtdTrigger_mVpdTacSum);
  fMCPico->SetBranchAddress("MtdTrigger.mTHUBtime[2]", MtdTrigger_mTHUBtime, &b_MtdTrigger_mTHUBtime);
  fMCPico->SetBranchAddress("MtdTrigger.mQTtacSum[8][8]", MtdTrigger_mQTtacSum, &b_MtdTrigger_mQTtacSum);
  fMCPico->SetBranchAddress("MtdTrigger.mMT101Tac[8][2]", MtdTrigger_mMT101Tac, &b_MtdTrigger_mMT101Tac);
  fMCPico->SetBranchAddress("MtdTrigger.mMT101Id[8][2]", MtdTrigger_mMT101Id, &b_MtdTrigger_mMT101Id);
  fMCPico->SetBranchAddress("MtdTrigger.mTF201TriggerBit", MtdTrigger_mTF201TriggerBit, &b_MtdTrigger_mTF201TriggerBit);
  fMCPico->SetBranchAddress("MtdTrigger.mShouldHaveRejectEvent", MtdTrigger_mShouldHaveRejectEvent, &b_MtdTrigger_mShouldHaveRejectEvent);
  fMCPico->SetBranchAddress("BTowHit", &BTowHit_, &b_BTowHit_);
  fMCPico->SetBranchAddress("BTowHit.mAdc", BTowHit_mAdc, &b_BTowHit_mAdc);
  fMCPico->SetBranchAddress("BTowHit.mE", BTowHit_mE, &b_BTowHit_mE);
  fMCPico->SetBranchAddress("BTofHit", &BTofHit_, &b_BTofHit_);
  fMCPico->SetBranchAddress("BTofHit.mId", BTofHit_mId, &b_BTofHit_mId);
  fMCPico->SetBranchAddress("MtdHit", &MtdHit_, &b_MtdHit_);
  fMCPico->SetBranchAddress("MtdHit.mgChannel", MtdHit_mgChannel, &b_MtdHit_mgChannel);
  fMCPico->SetBranchAddress("MtdHit.mTriggerFlag", MtdHit_mTriggerFlag, &b_MtdHit_mTriggerFlag);
  fMCPico->SetBranchAddress("MtdHit.mLeadingEdgeTime.first", MtdHit_mLeadingEdgeTime_first, &b_MtdHit_mLeadingEdgeTime_first);
  fMCPico->SetBranchAddress("MtdHit.mLeadingEdgeTime.second", MtdHit_mLeadingEdgeTime_second, &b_MtdHit_mLeadingEdgeTime_second);
  fMCPico->SetBranchAddress("MtdHit.mTrailingEdgeTime.first", MtdHit_mTrailingEdgeTime_first, &b_MtdHit_mTrailingEdgeTime_first);
  fMCPico->SetBranchAddress("MtdHit.mTrailingEdgeTime.second", MtdHit_mTrailingEdgeTime_second, &b_MtdHit_mTrailingEdgeTime_second);
  fMCPico->SetBranchAddress("BbcHit", &BbcHit_, &b_BbcHit_);
  fMCPico->SetBranchAddress("BbcHit.mId", &BbcHit_mId, &b_BbcHit_mId);
  fMCPico->SetBranchAddress("BbcHit.mQTdata", &BbcHit_mQTdata, &b_BbcHit_mQTdata);
  fMCPico->SetBranchAddress("EpdHit", &EpdHit_, &b_EpdHit_);
  fMCPico->SetBranchAddress("EpdHit.mId", &EpdHit_mId, &b_EpdHit_mId);
  fMCPico->SetBranchAddress("EpdHit.mQTdata", &EpdHit_mQTdata, &b_EpdHit_mQTdata);
  fMCPico->SetBranchAddress("EpdHit.mnMIP", &EpdHit_mnMIP, &b_EpdHit_mnMIP);
  fMCPico->SetBranchAddress("FmsHit", &FmsHit_, &b_FmsHit_);
  fMCPico->SetBranchAddress("FmsHit.mChannelDetectorId", &FmsHit_mChannelDetectorId, &b_FmsHit_mChannelDetectorId);
  fMCPico->SetBranchAddress("FmsHit.mAdc", &FmsHit_mAdc, &b_FmsHit_mAdc);
  fMCPico->SetBranchAddress("EmcPidTraits", &EmcPidTraits_, &b_EmcPidTraits_);
  fMCPico->SetBranchAddress("EmcPidTraits.mTrackIndex", EmcPidTraits_mTrackIndex, &b_EmcPidTraits_mTrackIndex);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcId", EmcPidTraits_mBemcId, &b_EmcPidTraits_mBemcId);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcAdc0", EmcPidTraits_mBemcAdc0, &b_EmcPidTraits_mBemcAdc0);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcE0", EmcPidTraits_mBemcE0, &b_EmcPidTraits_mBemcE0);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcE", EmcPidTraits_mBemcE, &b_EmcPidTraits_mBemcE);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcZDist", EmcPidTraits_mBemcZDist, &b_EmcPidTraits_mBemcZDist);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcPhiDist", EmcPidTraits_mBemcPhiDist, &b_EmcPidTraits_mBemcPhiDist);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcSmdNEta", EmcPidTraits_mBemcSmdNEta, &b_EmcPidTraits_mBemcSmdNEta);
  fMCPico->SetBranchAddress("EmcPidTraits.mBemcSmdNPhi", EmcPidTraits_mBemcSmdNPhi, &b_EmcPidTraits_mBemcSmdNPhi);
  fMCPico->SetBranchAddress("EmcPidTraits.mBtowId", EmcPidTraits_mBtowId, &b_EmcPidTraits_mBtowId);
  fMCPico->SetBranchAddress("EmcPidTraits.mBtowId23", EmcPidTraits_mBtowId23, &b_EmcPidTraits_mBtowId23);
  fMCPico->SetBranchAddress("EmcPidTraits.mBtowE", EmcPidTraits_mBtowE, &b_EmcPidTraits_mBtowE);
  fMCPico->SetBranchAddress("EmcPidTraits.mBtowE2", EmcPidTraits_mBtowE2, &b_EmcPidTraits_mBtowE2);
  fMCPico->SetBranchAddress("EmcPidTraits.mBtowE3", EmcPidTraits_mBtowE3, &b_EmcPidTraits_mBtowE3);
  fMCPico->SetBranchAddress("EmcPidTraits.mBtowEtaDist", EmcPidTraits_mBtowEtaDist, &b_EmcPidTraits_mBtowEtaDist);
  fMCPico->SetBranchAddress("EmcPidTraits.mBtowPhiDist", EmcPidTraits_mBtowPhiDist, &b_EmcPidTraits_mBtowPhiDist);
  fMCPico->SetBranchAddress("BTofPidTraits", &BTofPidTraits_, &b_BTofPidTraits_);
  fMCPico->SetBranchAddress("BTofPidTraits.mTrackIndex", BTofPidTraits_mTrackIndex, &b_BTofPidTraits_mTrackIndex);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofCellId", BTofPidTraits_mBTofCellId, &b_BTofPidTraits_mBTofCellId);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofMatchFlag", BTofPidTraits_mBTofMatchFlag, &b_BTofPidTraits_mBTofMatchFlag);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTof", BTofPidTraits_mBTof, &b_BTofPidTraits_mBTof);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofBeta", BTofPidTraits_mBTofBeta, &b_BTofPidTraits_mBTofBeta);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofYLocal", BTofPidTraits_mBTofYLocal, &b_BTofPidTraits_mBTofYLocal);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofZLocal", BTofPidTraits_mBTofZLocal, &b_BTofPidTraits_mBTofZLocal);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofHitPosX", BTofPidTraits_mBTofHitPosX, &b_BTofPidTraits_mBTofHitPosX);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofHitPosY", BTofPidTraits_mBTofHitPosY, &b_BTofPidTraits_mBTofHitPosY);
  fMCPico->SetBranchAddress("BTofPidTraits.mBTofHitPosZ", BTofPidTraits_mBTofHitPosZ, &b_BTofPidTraits_mBTofHitPosZ);
  fMCPico->SetBranchAddress("BTofPidTraits.mNSigmaElectron", BTofPidTraits_mNSigmaElectron, &b_BTofPidTraits_mNSigmaElectron);
  fMCPico->SetBranchAddress("BTofPidTraits.mNSigmaPion", BTofPidTraits_mNSigmaPion, &b_BTofPidTraits_mNSigmaPion);
  fMCPico->SetBranchAddress("BTofPidTraits.mNSigmaKaon", BTofPidTraits_mNSigmaKaon, &b_BTofPidTraits_mNSigmaKaon);
  fMCPico->SetBranchAddress("BTofPidTraits.mNSigmaProton", BTofPidTraits_mNSigmaProton, &b_BTofPidTraits_mNSigmaProton);
  fMCPico->SetBranchAddress("MtdPidTraits", &MtdPidTraits_, &b_MtdPidTraits_);
  fMCPico->SetBranchAddress("MtdPidTraits.mTrackIndex", MtdPidTraits_mTrackIndex, &b_MtdPidTraits_mTrackIndex);
  fMCPico->SetBranchAddress("MtdPidTraits.mMtdHitIndex", MtdPidTraits_mMtdHitIndex, &b_MtdPidTraits_mMtdHitIndex);
  fMCPico->SetBranchAddress("MtdPidTraits.mMatchFlag", MtdPidTraits_mMatchFlag, &b_MtdPidTraits_mMatchFlag);
  fMCPico->SetBranchAddress("MtdPidTraits.mDeltaY", MtdPidTraits_mDeltaY, &b_MtdPidTraits_mDeltaY);
  fMCPico->SetBranchAddress("MtdPidTraits.mDeltaZ", MtdPidTraits_mDeltaZ, &b_MtdPidTraits_mDeltaZ);
  fMCPico->SetBranchAddress("MtdPidTraits.mDeltaTimeOfFlight", MtdPidTraits_mDeltaTimeOfFlight, &b_MtdPidTraits_mDeltaTimeOfFlight);
  fMCPico->SetBranchAddress("MtdPidTraits.mBeta", MtdPidTraits_mBeta, &b_MtdPidTraits_mBeta);
  fMCPico->SetBranchAddress("MtdPidTraits.mMtdHitChan", MtdPidTraits_mMtdHitChan, &b_MtdPidTraits_mMtdHitChan);
  fMCPico->SetBranchAddress("TrackCovMatrix", &TrackCovMatrix_, &b_TrackCovMatrix_);
  fMCPico->SetBranchAddress("TrackCovMatrix.mImp", TrackCovMatrix_mImp, &b_TrackCovMatrix_mImp);
  fMCPico->SetBranchAddress("TrackCovMatrix.mZ", TrackCovMatrix_mZ, &b_TrackCovMatrix_mZ);
  fMCPico->SetBranchAddress("TrackCovMatrix.mPsi", TrackCovMatrix_mPsi, &b_TrackCovMatrix_mPsi);
  fMCPico->SetBranchAddress("TrackCovMatrix.mPti", TrackCovMatrix_mPti, &b_TrackCovMatrix_mPti);
  fMCPico->SetBranchAddress("TrackCovMatrix.mTan", TrackCovMatrix_mTan, &b_TrackCovMatrix_mTan);
  fMCPico->SetBranchAddress("TrackCovMatrix.mCurv", TrackCovMatrix_mCurv, &b_TrackCovMatrix_mCurv);
  fMCPico->SetBranchAddress("TrackCovMatrix.mSigma[5]", TrackCovMatrix_mSigma, &b_TrackCovMatrix_mSigma);
  fMCPico->SetBranchAddress("TrackCovMatrix.mCorr[10]", TrackCovMatrix_mCorr, &b_TrackCovMatrix_mCorr);
  fMCPico->SetBranchAddress("BEmcSmdEHit", &BEmcSmdEHit_, &b_BEmcSmdEHit_);
  fMCPico->SetBranchAddress("BEmcSmdEHit.mId", &BEmcSmdEHit_mId, &b_BEmcSmdEHit_mId);
  fMCPico->SetBranchAddress("BEmcSmdEHit.mAdc", &BEmcSmdEHit_mAdc, &b_BEmcSmdEHit_mAdc);
  fMCPico->SetBranchAddress("BEmcSmdEHit.mEnergy", &BEmcSmdEHit_mEnergy, &b_BEmcSmdEHit_mEnergy);
  fMCPico->SetBranchAddress("BEmcSmdPHit", &BEmcSmdPHit_, &b_BEmcSmdPHit_);
  fMCPico->SetBranchAddress("BEmcSmdPHit.mId", &BEmcSmdPHit_mId, &b_BEmcSmdPHit_mId);
  fMCPico->SetBranchAddress("BEmcSmdPHit.mAdc", &BEmcSmdPHit_mAdc, &b_BEmcSmdPHit_mAdc);
  fMCPico->SetBranchAddress("BEmcSmdPHit.mEnergy", &BEmcSmdPHit_mEnergy, &b_BEmcSmdPHit_mEnergy);
  fMCPico->SetBranchAddress("ETofHit", &ETofHit_, &b_ETofHit_);
  fMCPico->SetBranchAddress("ETofHit.mGeomId", &ETofHit_mGeomId, &b_ETofHit_mGeomId);
  fMCPico->SetBranchAddress("ETofHit.mLocalX", &ETofHit_mLocalX, &b_ETofHit_mLocalX);
  fMCPico->SetBranchAddress("ETofHit.mLocalY", &ETofHit_mLocalY, &b_ETofHit_mLocalY);
  fMCPico->SetBranchAddress("ETofHit.mClusterSize", &ETofHit_mClusterSize, &b_ETofHit_mClusterSize);
  fMCPico->SetBranchAddress("ETofHit.mLeadingEdgeTime", &ETofHit_mLeadingEdgeTime, &b_ETofHit_mLeadingEdgeTime);
  fMCPico->SetBranchAddress("ETofHit.mTimeOverThreshold", &ETofHit_mTimeOverThreshold, &b_ETofHit_mTimeOverThreshold);
  fMCPico->SetBranchAddress("ETofPidTraits", &ETofPidTraits_, &b_ETofPidTraits_);
  fMCPico->SetBranchAddress("ETofPidTraits.mTrackIndex", &ETofPidTraits_mTrackIndex, &b_ETofPidTraits_mTrackIndex);
  fMCPico->SetBranchAddress("ETofPidTraits.mHitIndex", &ETofPidTraits_mHitIndex, &b_ETofPidTraits_mHitIndex);
  fMCPico->SetBranchAddress("ETofPidTraits.mMatchFlag", &ETofPidTraits_mMatchFlag, &b_ETofPidTraits_mMatchFlag);
  fMCPico->SetBranchAddress("ETofPidTraits.mTimeOfFlight", &ETofPidTraits_mTimeOfFlight, &b_ETofPidTraits_mTimeOfFlight);
  fMCPico->SetBranchAddress("ETofPidTraits.mBeta", &ETofPidTraits_mBeta, &b_ETofPidTraits_mBeta);
  fMCPico->SetBranchAddress("ETofPidTraits.mDeltaX", &ETofPidTraits_mDeltaX, &b_ETofPidTraits_mDeltaX);
  fMCPico->SetBranchAddress("ETofPidTraits.mDeltaY", &ETofPidTraits_mDeltaY, &b_ETofPidTraits_mDeltaY);
  fMCPico->SetBranchAddress("ETofPidTraits.mCrossingX", &ETofPidTraits_mCrossingX, &b_ETofPidTraits_mCrossingX);
  fMCPico->SetBranchAddress("ETofPidTraits.mCrossingY", &ETofPidTraits_mCrossingY, &b_ETofPidTraits_mCrossingY);
  fMCPico->SetBranchAddress("ETofPidTraits.mCrossingZ", &ETofPidTraits_mCrossingZ, &b_ETofPidTraits_mCrossingZ);
  fMCPico->SetBranchAddress("McVertex", &McVertex_, &b_McVertex_);
  fMCPico->SetBranchAddress("McVertex.mId", McVertex_mId, &b_McVertex_mId);
  fMCPico->SetBranchAddress("McVertex.mNoDaughters", McVertex_mNoDaughters, &b_McVertex_mNoDaughters);
  fMCPico->SetBranchAddress("McVertex.mIdParTrk", McVertex_mIdParTrk, &b_McVertex_mIdParTrk);
  fMCPico->SetBranchAddress("McVertex.mIsInterm", McVertex_mIsInterm, &b_McVertex_mIsInterm);
  fMCPico->SetBranchAddress("McVertex.mTime", McVertex_mTime, &b_McVertex_mTime);
  fMCPico->SetBranchAddress("McVertex.mVx", McVertex_mVx, &b_McVertex_mVx);
  fMCPico->SetBranchAddress("McVertex.mVy", McVertex_mVy, &b_McVertex_mVy);
  fMCPico->SetBranchAddress("McVertex.mVz", McVertex_mVz, &b_McVertex_mVz);
  fMCPico->SetBranchAddress("McTrack", &McTrack_, &b_McTrack_);
  fMCPico->SetBranchAddress("McTrack.mId", McTrack_mId, &b_McTrack_mId);
  fMCPico->SetBranchAddress("McTrack.mGePid", McTrack_mGePid, &b_McTrack_mGePid);
  fMCPico->SetBranchAddress("McTrack.mCharge", McTrack_mCharge, &b_McTrack_mCharge);
  fMCPico->SetBranchAddress("McTrack.mHits[22]", McTrack_mHits, &b_McTrack_mHits);
  fMCPico->SetBranchAddress("McTrack.mPx", McTrack_mPx, &b_McTrack_mPx);
  fMCPico->SetBranchAddress("McTrack.mPy", McTrack_mPy, &b_McTrack_mPy);
  fMCPico->SetBranchAddress("McTrack.mPz", McTrack_mPz, &b_McTrack_mPz);
  fMCPico->SetBranchAddress("McTrack.mE", McTrack_mE, &b_McTrack_mE);
  fMCPico->SetBranchAddress("McTrack.mIsFromShower", McTrack_mIsFromShower, &b_McTrack_mIsFromShower);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStart", McTrack_mIdVtxStart, &b_McTrack_mIdVtxStart);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStop", McTrack_mIdVtxStop, &b_McTrack_mIdVtxStop);
  fMCPico->SetBranchAddress("McTrack.mIdVtxItrmd", McTrack_mIdVtxItrmd, &b_McTrack_mIdVtxItrmd);


  int counter = 0;
  int eventsconsidered = 0;

  for (Long64_t jentry=0; jentry<randomlisttoevents.size();jentry++) {
  // while(counter < fSetNumberOfEvents){

    // cout << "Events considered = " << eventsconsidered << endl;

    // cout << "Random Number = " << randomlisttoevents[jentry] << "\t" << counter << endl;

    fMCPico->GetEntry(randomlisttoevents[jentry]);
    // fMCPico->GetEntry(randomlisttoevents[eventsconsidered]);

    // MC Entries here
    vector<int> vertexids;
    vector<int> pionids;
    vector<int> kaonids;
    vertexids.clear();
    pionids.clear();
    kaonids.clear();

    // cout << Event_mPrimaryVertexX[0] << "\t" << Event_mPrimaryVertexY[0] << "\t" << Event_mPrimaryVertexZ[0] << endl;

    // Saving the D0 vertices
    for (int mc = 0; mc < McTrack_; mc++){
      if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38){
        TLorentzVector v;
        v.SetPxPyPzE(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], McTrack_mE[mc]);
        // track variables
        double pt = v.Perp();
        double phi = v.Phi();
        if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
        if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
        double eta = v.PseudoRapidity();

        if (pt < 1.0 || pt > 30.0) continue;
        if ((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax) || (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) continue;
        
        cout << "D0 Info = " << pt << "\t" << eta << "\t" << phi << "\t" << McTrack_mIdVtxStop[mc] << endl;

        bool goodvertex = kTRUE;

        int pionid = -99;
        int kaonid = -99;

        for(int daug = 0; daug < McTrack_; daug++){
          if (McTrack_mIdVtxStart[daug] != McTrack_mIdVtxStop[mc]) continue; // We only want the kaon and pion that originated from the D0 we are interested in.

          cout << "Daughters are = " << McTrack_mGePid[daug] << endl;

          if (McTrack_mIdVtxStop[daug] != 0) {goodvertex = kFALSE; break;} // If these are not final state kaon pions, then they are of no use to us.

          if (McTrack_mGePid[daug] != 8 && McTrack_mGePid[daug] != 9 && McTrack_mGePid[daug] != 11 && McTrack_mGePid[daug] != 12)
          { 
            goodvertex = kFALSE; break;
          }

          // Push the pion and kaon ids into the vector
          if (McTrack_mGePid[daug] == 8 || McTrack_mGePid[daug] == 9){
            pionid = daug;
          }
          
          else if (McTrack_mGePid[daug] == 11 || McTrack_mGePid[daug] == 12){
            kaonid = daug;
          } 
          
        }

        // cout << "Vertex Stop = " << McTrack_mIdVtxStop[mc] << endl;
        if (goodvertex){
          if (pionid == -99 || kaonid == -99) { cout << "Something wrong with D0. Exiting." << endl; continue; }
          vertexids.push_back(McTrack_mIdVtxStop[mc]);
          pionids.push_back(pionid);
          kaonids.push_back(kaonid);
        } 
        // cout << "D0 = " << pt << "\t" << eta << "\t" << phi << endl;
      }
    }

    // cout << "Number of D0s = " << vertexids.size() << endl;

    int D0Counter = 0;


    // Loop over each D0
    for (int D0 = 0; D0 < vertexids.size(); D0++){

      // I only want to use TLorentzVector to propagate the information ahead. All the processing is done within this class.
      // This means I need to find a way to make sure I can identify the D0 track when I save out the information.
      // Since only D0 needs to be identified, I am saving the mass information for the D0 as 1.865.
      // All the other tracks are saved out with pi0mass, because ultimately, jet constituents are assumed to have that mass.

      // How do I propagate charge information though?
      // Do I need it? Mostly, nope! In fact, once I separate track and tower, it should be enough.

      // This loop fills the input vector for the MC Side

      // cout << "Number of MC tracks = " << McTrack_ << endl;

      for (int mc = 0; mc < McTrack_; mc++){
        TLorentzVector v;
        if ((McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38)) v.SetXYZM(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], 1.865); // This is a D0
        else v.SetXYZM(McTrack_mPx[mc], McTrack_mPy[mc], McTrack_mPz[mc], pi0mass);

        // track variables
        double pt = v.Perp();
        double phi = v.Phi();
        if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
        if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
        double eta = v.PseudoRapidity();

        // if (McTrack_mGePid[mc] == 8 || McTrack_mGePid[mc] == 9 || McTrack_mGePid[mc] == 11 || McTrack_mGePid[mc] == 12){
        //   cout << "Vertex Start Stop = " << McTrack_mIdVtxStart[mc] << "\t" << McTrack_mIdVtxStop[mc] << endl;
        // }

        if (McTrack_mIdVtxStop[mc] != 0 && (McTrack_mGePid[mc] != 37 && McTrack_mGePid[mc] != 38)) continue; // Only Final State Particles are included in JetMaker. D0s are the only exception

        if ((McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) && pt < 1.0) continue; // We don't need D0s < 1 GeV
        if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) {
          if (McTrack_mIdVtxStop[mc] != vertexids[D0]) continue; // Only consider the current D0
        }

        // if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38){D0Counter++;}
        // Unstable particles which shouldn't make it to the end are discarded by hand. The list provisionally includes:
        /*
          Lambda, Eta, Sigma0, Xi0, Muon, Neutrino, KS0, KL0
        */

        if (McTrack_mGePid[mc] == 4 || McTrack_mGePid[mc] == 5 || McTrack_mGePid[mc] == 6 || McTrack_mGePid[mc] == 10 || McTrack_mGePid[mc] == 16 || McTrack_mGePid[mc]== 17 || McTrack_mGePid[mc] == 18 || McTrack_mGePid[mc] == 20 || McTrack_mGePid[mc] == 22 ) continue;

        // Have to discard Kaons and Pions which come from the mCurvrent D0

        if (McTrack_mIdVtxStart[mc] == vertexids[D0]){
          continue;
        }

        if ((pt < fMinJetTrackPt) || (pt > fMaxJetTrackPt) || (eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax) || (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) continue;
        if (McTrack_mCharge[mc] != 0 || McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) {
          if (McTrack_mGePid[mc] == 37 || McTrack_mGePid[mc] == 38) {
            D0Counter++; 
            // cout << "D0 = " << pt << "\t" << eta << "\t" << phi << endl;
          }
          fMCEventTracks[counter].push_back(v); //Neutral particles which are not D0 are saved in towers.
        }
        else fMCEventTowers[counter].push_back(v); 
      }

      const int numberoftowers = BTowHit_;
      double towerenergy[numberoftowers];
      for (int i = 0; i < numberoftowers; i++){towerenergy[i] = 0.;}

      // This loop fills the input vector for the RECO side
      for (int reco = 0; reco < Track_; reco++){

        TVector3 o;
        o.SetXYZ(Track_mOriginX[reco], Track_mOriginY[reco], Track_mOriginZ[reco]);
        TVector3 g;
        g.SetXYZ(Track_mGMomentumX[reco], Track_mGMomentumY[reco], Track_mGMomentumZ[reco]);

        // track variables
        double pt = g.Perp();
        double phi = g.Phi();
        if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
        if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
        double eta = g.PseudoRapidity();
        double px = g.x();
        double py = g.y();
        double pz = g.z();
        double p = g.Mag();
        double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
        short charge = (Track_mNHitsFit[reco] > 0) ? 1 : -1;

        double dca = (mVertex - o).Mag();
        
        //// Variables For FastSim

        double pt_new = pt;
        double phi_new = phi;
        double eta_new = eta;
        double px_new = px;
        double py_new = py;
        double pz_new = pz;
        double p_new = p;
        double energy_new = energy;
        short charge_new = charge;

        bool mctrackavailable = kTRUE;

        int mcid = Track_mIdTruth[reco];

        if (mcid <= 0) mctrackavailable = kFALSE;

        bool isatrackfromD0 = kFALSE;

        // Here, we have two paths to take. If the track needs replacement, we replace it with the fastsim method that is standardised.
        // Else, the pt, eta, phi are sent as is to the final vector.

        if (mctrackavailable){
          if (McTrack_mIdVtxStart[mcid] == vertexids[D0]) isatrackfromD0 = kTRUE; // Kaons and Pions that come from the current D0 need to be tossed, and replaced by the fast sim version 

          TVector3 mg(McTrack_mPx[mcid], McTrack_mPy[mcid], McTrack_mPz[mcid]);

          double relativesmearing = TMath::Sqrt(pow(mg.Px() - px, 2) + pow(mg.Py() - py, 2))/(mg.Pt());

          double fastsimsmearing;

          int pid = McTrack_mGePid[mcid];

          if     (pid ==8  || pid == 9) fastsimsmearing =  fPionMomResolution->Eval(mg.Pt());// Pion
          else if(pid ==11 || pid == 12)fastsimsmearing =  fKaonMomResolution->Eval(mg.Pt());// Kaon
          else if(pid ==15 || pid == 14)fastsimsmearing =  fProtonMomResolution->Eval(mg.Pt());// Proton
          else fastsimsmearing = fPionMomResolution->Eval(mg.Pt());// Catch all: pions

          if (relativesmearing > 3*fastsimsmearing){
            TVector3 fastsimsmearedmom = FastSimMom(mg, pid);
            pt_new = fastsimsmearedmom.Perp();
            phi_new = fastsimsmearedmom.Phi();
            if(phi_new < 0.0)    phi_new += 2.0*pi;  // force from 0-2pi
            if(phi_new > 2.0*pi) phi_new -= 2.0*pi;  // force from 0-2pi
            eta_new = fastsimsmearedmom.PseudoRapidity();
            px_new = fastsimsmearedmom.x();
            py_new = fastsimsmearedmom.y();
            pz_new = fastsimsmearedmom.z();
            p_new = fastsimsmearedmom.Mag();
            energy_new = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
          }
        }

        // jet track acceptance cuts now
        if(pt_new < fMinJetTrackPt) continue;
        if(pt_new > fMaxJetTrackPt) continue; // 20.0 STAR, 100.0 ALICE
        if((eta_new < fJetTrackEtaMin) || (eta_new > fJetTrackEtaMax)) continue;
        if(phi_new < 0.0)    phi_new += 2.0*pi;  // force from 0-2pi
        if(phi_new > 2.0*pi) phi_new -= 2.0*pi;  // force from 0-2pi
        if((phi_new < fJetTrackPhiMin) || (phi_new > fJetTrackPhiMax)) continue;
            
        // additional quality cuts for tracks
        if(dca > fJetTrackDCAcut)            continue;
        if(abs(Track_mNHitsFit[reco]) < fJetTracknHitsFit)     continue;
        if(abs(Track_mNHitsFit[reco])/Track_mNHitsMax[reco] < fJetTracknHitsRatio) continue;

        // cout << "TRACK = " << px_new << "\t" << py_new << "\t" << pz_new << endl;

        int particleid = -99;

        if      (abs(Track_mNSigmaPion[reco]/1000.) < 2) particleid = 1;
        else if (abs(Track_mNSigmaKaon[reco]/1000.) < 2.) particleid = 2;
        else if (abs(Track_mNSigmaProton[reco]/1000.) < 2.) particleid = 3*charge;

        if (!KeepTrack(particleid, centralitybinforefficiency, pt_new) && !isatrackfromD0) continue; // To match the efficiency, we start tossing random tracks.

        TLorentzVector v;
        v.SetXYZM(px_new, py_new, pz_new, pi0mass);

        if (!isatrackfromD0) fRecoEventTracks[counter].push_back(v);
        //// Fill the Tower Array with Energy Depositions here

        int matchedTowerIndex = abs(GetMatchedBtowID(Track_mBEmcMatchedTowerIndex[reco], g, o, charge_new)) - 1; // towerIndex = towerID - 1
        
        if (matchedTowerIndex >= 0){towerenergy[matchedTowerIndex]+=energy_new;} //This place takes care of all energy depositions due to tracks accepted for jet reco. It also includes the kaon and pion from D0.        
      }

      TVector3 mcKaon(McTrack_mPx[kaonids[D0]], McTrack_mPy[kaonids[D0]], McTrack_mPz[kaonids[D0]]);
      TVector3 mcPion(McTrack_mPx[pionids[D0]], McTrack_mPy[pionids[D0]], McTrack_mPz[pionids[D0]]);

      fMCD0Information[counter] = {mcPion, mcKaon}; 

      TVector3 recoKaon = FastSimMom(mcKaon, 11);
      TVector3 recoPion = FastSimMom(mcPion, 8);

      fRecoD0Information[counter] = {recoPion, recoKaon};

      fOrigin[counter].SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

      // fOrigin[counter].push_back(eventorigin); 

      TVector3 recoD0;
      recoD0 = recoKaon + recoPion;

      TLorentzVector v;
      v.SetXYZM(recoD0.x(), recoD0.y(), recoD0.z(), 1.865); // This is a D0

      fRecoEventTracks[counter].push_back(v);

      // This loop fills the input vector for the towers

      for (int tower = 0; tower < numberoftowers; tower++){
        int towerID = tower + 1;
        if(towerID < 0) continue; // double check these aren't still in the event list

        TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
        double towerPhi = towerPosition.Phi();
        if(towerPhi < 0.0)    towerPhi += 2.0*pi;  // force from 0-2pi
        if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;  // force from 0-2pi
        double towerEta = towerPosition.PseudoRapidity();

        // check for bad (and dead) towers
        bool TowerOK = mBaseMaker->IsTowerOK(towerID);      // kTRUE means GOOD
        bool TowerDead = mBaseMaker->IsTowerDead(towerID);  // kTRUE means BAD
        if(!TowerOK)  { continue; }
        if(TowerDead) { continue; }

        // jet track acceptance cuts njow
        if((towerEta < fJetTowerEtaMin) || (towerEta > fJetTowerEtaMax)) continue;
        if((towerPhi < fJetTowerPhiMin) || (towerPhi > fJetTowerPhiMax)) continue;

        double towerEunCorr = BTowHit_mE[tower]/1000.;  // uncorrected energy
        double towerE = BTowHit_mE[tower]/1000.;        // corrected energy (hadronically - done below)
        double towEtunCorr = towerE / (1.0*TMath::CosH(towerEta));

        // cut on min tower energy after filling histos
        if(towerEunCorr < mTowerEnergyMin) continue; // if we don't have enough E to start with, why mess around

        // =======================================================================
        // HADRONIC CORRECTION
        
        double sumEt = (towerEunCorr - towerenergy[tower])/(1.0*TMath::CosH(towerEta));
        double towerEt = sumEt;
        if(towerEt < mTowerEnergyMin) continue;
        towerE = towerEt  * 1.0*TMath::CosH(towerEta);

        Double_t p = 1.0*TMath::Sqrt(towerE*towerE - pi0mass*pi0mass);

        double posX = towerPosition.x();
        double posY = towerPosition.y();
        double posZ = towerPosition.z();

        Double_t r = TMath::Sqrt(posX*posX + posY*posY + posZ*posZ) ;

        TLorentzVector v;
        v.SetXYZM(p*posX/r, p*posY/r, p*posZ/r, pi0mass);

        // cout << "TOWER = " << v.X() << "\t" << v.Y() << "\t" << v.Z() << endl;

        fRecoEventTowers[counter].push_back(v);
      }

      counter++; // Since one MC event can be invoked multiple times (due to having multiple D0s, this step is necessary.)
    }
    eventsconsidered++;
  }

  fNumberofeventsoverlayed = counter;

  fMCPico->Reset();
  f->Close();

}

int StHIOverlay::GetMatchedBtowID(int trkbemcid, TVector3 gMom, TVector3 org, int charge){
  Double_t bemc_radius = mBemcGeom->Radius();
  // Magnetic field in Tesla 
  Double_t mBField_tesla = Bfield / 10.0; //Check this definition. Magnetic fields are minefields of error in STAR

  // Needed for projection of the track onto the barrel radius
  TVector3 bemc_pos, bemc_mom;

  // BEMC hardware indices 
  Int_t h_m, h_e, h_s = 0;

  // tower index: if no tower can be matched, assign 0
  // picoTrk->setBEmcMatchedTowerIndex(0);
  Int_t tow_id = 0;
  Bool_t close_match = false;

  // int trkbemcid = trk->bemcTowerIndex();

  // Check if the track can be projected onto the current radius
  // if not, track can't be matched.
  // By JetCorr request the global track projection to BEMC is used.
  if ( mEmcPosition->projTrack(&bemc_pos, &bemc_mom, gMom, org, Bfield, charge, bemc_radius) ) {
    // First, examine track eta. If it falls in two regions:
    // 0 < |eta| < etaMin()
    // etaMax() < |eta| < 1.0
    // then shift the eta for the projection slightly into the neighboring tower
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(org, trkbemcid + 1);
    // cout << bemc_pos.Phi() << "\t" << towerPosition.Phi() << "\t" << bemc_pos.PseudoRapidity() << "\t" << towerPosition.PseudoRapidity() << endl;

    if ( fabs(bemc_pos.PseudoRapidity()) < mBemcGeom->EtaMin() ) {
      Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.001;
      // Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.000001;
      Double_t unsigned_theta = 2.0 * atan(exp(-1.0 * unsigned_eta));
      Double_t signed_theta = (bemc_pos.PseudoRapidity() >= 0 ? 1.0 : -1.0) * unsigned_theta;
      bemc_pos.SetTheta(signed_theta);
      close_match = true;
    } 
    else if ( fabs(bemc_pos.PseudoRapidity()) > mBemcGeom->EtaMax() &&
      fabs(bemc_pos.PseudoRapidity()) < 1.0 ) {
      Double_t unsigned_eta = mBemcGeom->EtaMax() - 0.001;
      // Double_t unsigned_eta = mBemcGeom->EtaMax() - 0.000001;
      Double_t unsigned_theta = 2.0 * atan(exp(-1.0 * unsigned_eta));
      Double_t signed_theta = (bemc_pos.PseudoRapidity() >= 0 ? 1.0 : -1.0) * unsigned_theta;
      bemc_pos.SetTheta(signed_theta);
      close_match = true;
    }


    // Get the BEMC hardware location in (m, e, s) and translate to id
    // If StEmcGeom::getBin() != 0: track was not matched to a tower.
    // Its outside of the BEMC eta range (> 1.0).

    if ( mBemcGeom->getBin(bemc_pos.Phi(),bemc_pos.PseudoRapidity(),h_m,h_e,h_s) == 0 ) {
      // If StEmcGeom::getId() == 0: the track was matched successfully. Otherwise, 
      // the track was not matched to a tower at this radius, the track was projected
      // into the gap between modules in phi. 
      if ( h_s != -1 ) {
        mBemcGeom->getId(h_m,h_e,h_s,tow_id);
        if (close_match) {
          return -1*tow_id;
        }

        else{
          return tow_id;
        }
      }

      // Track fell in between modules in phi. We will find which module it is closer
      // to by shifting phi slightly.
      else {
        // Value of the "dead space" per module in phi:
        // 2*pi/60 (amount of azimuth covered per module)
        // 2*0.0495324 (active size of module)

        Double_t dphi = (TMath::Pi() / 60.0) - 0.0495324;
        // Shift the projected phi by dphi in positive and negative directions
        // if we look for the projection for both of these, only one should give
        // a tower id, and the other should still be in the inter-tower space

        TVector3 bemc_pos_shift_pos(bemc_pos); 
        bemc_pos_shift_pos.SetPhi(bemc_pos_shift_pos.Phi() + dphi);
        TVector3 bemc_pos_shift_neg(bemc_pos); 
        bemc_pos_shift_neg.SetPhi(bemc_pos_shift_neg.Phi() - dphi);

        if ( mBemcGeom->getBin(bemc_pos_shift_pos.Phi(),bemc_pos_shift_pos.PseudoRapidity(),h_m,h_e,h_s) == 0 && h_s != -1 ){
          mBemcGeom->getId(h_m,h_e,h_s,tow_id);        
          return -1*tow_id;
        }

        else if( mBemcGeom->getBin(bemc_pos_shift_neg.Phi(),bemc_pos_shift_neg.PseudoRapidity(),h_m,h_e,h_s) == 0 && h_s != -1 ){
          mBemcGeom->getId(h_m,h_e,h_s,tow_id);
          return -1*tow_id;
        } 
      }
    }
  }

  return tow_id;

}


TVector3 StHIOverlay::FastSimMom(TVector3 p, int pid)
{
    float pt = p.Perp();
    float pt1 = pt;
    // if(pt1>2) pt1 = 2;//Used for high pt-hat bin smearing test
    if(pt1>10) pt1 = 10;//Used for high pt-hat bin smearing test
    float sPt = -1;
    
    if(pid ==8 || pid == 9)sPt = gRandom->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Pion
    else if(pid ==11 || pid == 12)sPt = gRandom->Gaus(pt, pt * fKaonMomResolution->Eval(pt1));// Kaon
    else if(pid ==15 || pid == 14)sPt = gRandom->Gaus(pt, pt * fProtonMomResolution->Eval(pt1));// Proton
    else sPt = gRandom->Gaus(pt, pt * fPionMomResolution->Eval(pt1));// Catch all: pions

    TVector3 smearedmom(sPt * cos(p.Phi()), sPt * sin(p.Phi()), sPt * sinh(p.PseudoRapidity()));
    
    return smearedmom;
}


Bool_t StHIOverlay::KeepTrack(int particleid, int centralitybin, double pt)
{
  bool keeptrack = kTRUE;

  TRandom3 *r = new TRandom3(0);
  double rando = r->Rndm();

  // cout << "Random = " << rando << "\t" << fPionWeight->Eval(pt) << "\t" << fKaonWeight->Eval(pt) << "\t" << fProtonWeight->Eval(pt) << "\t" << fAProtonWeight->Eval(pt) << endl;

  if (particleid == 1 || particleid == -99) keeptrack = (rando > fPionWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE; //Either a pion
  else if (particleid == 2) keeptrack = (rando > fKaonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == 3) keeptrack = (rando > fProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == -3) keeptrack = (rando > fAProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;

  return keeptrack;
}

