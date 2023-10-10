// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StTestClass.h"
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
#include "StFJWrapper.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"


// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StTestClass)

//________________________________________________________________________
StTestClass::StTestClass(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StTestClass::fRunFlagEnum
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
  // fjw("HIMC", "HIMC");

  // fMCFileListName = filename;

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  mBemcGeom = 0x0;

  if (!name) return;
  SetName(name);
}

//
//________________________________________________________________________
StTestClass::~StTestClass()
{ /*  */

}

//
//________________________________________________________________________
Int_t StTestClass::Init() {

  StJetFrameworkPicoBase::Init();

  // position object for Emc
  mBemcGeom = StEmcGeom::instance("bemc");
  mEmcPosition = new StEmcPosition2();

  int numberofcentbins = 101;

  hCentralitySoft = new TH1F("hCentralitySoft", "hCentralitySoft", numberofcentbins, -0.5, 100.5);
  hCentrality     = new TH1F("hCentrality", "hCentrality", numberofcentbins, -0.5, 100.5);

  hSoftEventMaxPt = new TH1F("hSoftEventMaxPt", "hSoftEventMaxPt", 60, 0, 30);
  hMaxPt          = new TH1F("hMaxPt", "hMaxPt", 60, 0, 30);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTestClass::Finish() { 
  cout << "StTestClass::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    hCentralitySoft->Write();
    hCentrality->Write();
    hSoftEventMaxPt->Write();
    hMaxPt->Write();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StTestClass::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTestClass::Clear(Option_t *opt) {
  // fjw->Clear();
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTestClass::Make() {
  // fjw->Clear();
  // fJets->Delete();

  // fJetsArr.clear();
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // for (int i = 0; i < 100; i++){
  //   fMCEventTracks[i].clear();
  //   fMCEventTowers[i].clear();
  //   fRecoEventTracks[i].clear();
  //   fRecoEventTowers[i].clear();
  //   fMCD0Information[i] = {};
  //   fRecoD0Information[i] = {};
  //   fOrigin[i].SetXYZ(0,0,0);
  // }

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

  // double maxpt = GetMaxTrackPt();

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


  // if (centralitybinforefficiency != 0) return kStOk; ///// Remove this after the test

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

  // We only want events with no other hard scattering. This might skew things weirdly, but it's still worth a try

  hCentrality->Fill(fCentralityScaled);
  hMaxPt->Fill(GetMaxTrackPt());

  cout << "Max Track & Tower = " << GetMaxTrackPt() << "\t" << GetMaxTowerEt() << endl;

  if (GetMaxTrackPt() > 5.) return kStOk;

  hCentralitySoft->Fill(fCentralityScaled);
  hSoftEventMaxPt->Fill(GetMaxTrackPt());

  return kStOK;
}

Bool_t StTestClass::IsItASoftEvent(){
  bool softevent = kTRUE;

  const Int_t ntracks = mPicoDst->numberOfTracks();

  for(unsigned short itrk = 0; itrk < ntracks; itrk++){
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
    if(!trk){ continue; }

    TVector3 mTrkMom;
    mTrkMom = trk->gMom(mVertex, Bfield);

    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    if (abs(eta) > 1.) continue;
    if (pt > 5.) {
      softevent = kFALSE;
      break;
    }
  }

  return softevent;
}
