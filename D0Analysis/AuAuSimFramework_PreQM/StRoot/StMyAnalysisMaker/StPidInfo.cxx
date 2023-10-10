// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StPidInfo.h"
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
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"

// KFParticle includes
#include "StRoot/StEvent/StDcaGeometry.h"
#include "StRoot/StarRoot/KFParticle.h"

// Bichsel includes
#include "StBichsel/Bichsel.h"

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

ClassImp(StPidInfo)

//________________________________________________________________________
StPidInfo::StPidInfo(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  doUsePrimTracks = kTRUE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StPidInfo::fRunFlagEnum
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
  fEventZVtxMinCut = -60.0; fEventZVtxMaxCut = 60.0;
  fMaxEventTrackPt = 50.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 50.0;
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
  fJetMakerName = jetMakerName;
  //fRhoMakerName = rhoMakerName;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  for(int i=0; i<10; i++) {numberofevents[i] = 0; numberoftracks[i] = 0;}
  //Mass cut
  fInvMassSignal1 = 1.84;
  fInvMassSignal2 = 1.89;

  fInvMassULBg1 = 1.70;
  fInvMassULBg2 = 2.08;

  fInvMassLSBg1 = 1.70;
  fInvMassLSBg2 = 2.08;

}

//
//________________________________________________________________________
StPidInfo::~StPidInfo()
{ 
  //if (gPtvpPt) delete gPtvpPt;
  if (dEdXvpT) delete dEdXvpT;
  if (invbetavpT) delete invbetavpT;
  // if (event_cuts) delete event_cuts;
  // if (track_cuts) delete track_cuts;

  if (cuthistogram_event) delete cuthistogram_event;
  if (cuthistogram_track) delete cuthistogram_track;
  //gPtvpPt->Write();
  // if (dEdXvpT) delete dEdXvpT;
  if (dEdXvp) delete dEdXvp;

  if (dEdXvpT_pion) delete dEdXvpT_pion;
  if (dEdXvpT_kaon) delete dEdXvpT_kaon;
  if (dEdXvpT_proton) delete dEdXvpT_proton;
  if (dEdXvpT_electron) delete dEdXvpT_electron;

  if (dEdXvp_pion) delete dEdXvp_pion;
  if (dEdXvp_kaon) delete dEdXvp_kaon;
  if (dEdXvp_proton) delete dEdXvp_proton;
  if (dEdXvp_electron) delete dEdXvp_electron;

  if (dEdXthvp_pi) delete dEdXthvp_pi;
  if (dEdXthvp_ka) delete dEdXthvp_ka;
  if (dEdXthvp_pr) delete dEdXthvp_pr;

  if (z_pi) delete z_pi;
  if (z_ka) delete z_ka;
  if (z_pr) delete z_pr;

  // if (invbetavpT) delete invbetavpT;
  if (invbetavpT_tof) delete invbetavpT_tof;
  if (normalised_invbetavpT_tof_pi) delete normalised_invbetavpT_tof_pi;
  if (normalised_invbetavpT_tof_ka) delete normalised_invbetavpT_tof_ka;
  if (normalised_invbetavpT_tof_pr) delete normalised_invbetavpT_tof_pr;

  if (mvpT) delete mvpT;
  if (EvP) delete EvP;

  if (hJetShaped0) delete hJetShaped0;
  if (hJetShaped0BgUL) delete hJetShaped0BgUL;
  if (hJetShaped0BgLS) delete hJetShaped0BgLS;

  if (hJetDistd0) delete hJetDistd0;
  if (hJetDistd0BgUL) delete hJetDistd0BgUL;
  if (hJetDistd0BgLS) delete hJetDistd0BgLS;

}

//________________________________________________________________________
Int_t StPidInfo::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StPidInfo::Finish() { 

  cout << "StPidInfo::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "RECREATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StPidInfo::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StPidInfo::DeclareHistograms() {

  //gPtvpPt = new TH2F("gPtvpPt", "gPtvpPt", 1000, 0, 20, 1000, 0, 20);

  cuthistogram_event = new TH1F("cuthistogram_event", "cuthistogram_event", 10, 0, 10);

  cuthistogram_track = new TH1F("cuthistogram_track", "cuthistogram_track", 10, 0, 10);

  dEdXvpT = new TH2F("dEdXvpT", "dEdXvpT", 1000, 0, 5, 1000, 0, 20);

  dEdXvp = new TH2F("dEdXvp", "dEdXvp", 1000, 0, 5, 1000, 0, 20);

  dEdXvpT_pion = new TH2F("dEdXvpT_pion", "dEdXvpT_pion", 1000, 0, 5, 1000, 0, 20);

  dEdXvpT_kaon = new TH2F("dEdXvpT_kaon", "dEdXvpT_kaon", 1000, 0, 5, 1000, 0, 20);

  dEdXvpT_proton = new TH2F("dEdXvpT_proton", "dEdXvpT_proton", 1000, 0, 5, 1000, 0, 20);

  dEdXvpT_electron = new TH2F("dEdXvpT_electron", "dEdXvpT_electron", 1000, 0, 5, 1000, 0, 20);

  dEdXvp_pion = new TH2F("dEdXvp_pion", "dEdXvp_pion", 1000, 0, 5, 1000, 0, 20);

  dEdXvp_kaon = new TH2F("dEdXvp_kaon", "dEdXvp_kaon", 1000, 0, 5, 1000, 0, 20);

  dEdXvp_proton = new TH2F("dEdXvp_proton", "dEdXvp_proton", 1000, 0, 5, 1000, 0, 20);

  dEdXvp_electron = new TH2F("dEdXvp_electron", "dEdXvp_electron", 1000, 0, 5, 1000, 0, 20);

  dEdXthvp_pi = new TH2F("dEdXthvp_pi", "dEdXthvp_pi", 1000, 0, 5, 1000, 0, 20);

  dEdXthvp_ka = new TH2F("dEdXthvp_ka", "dEdXthvp_ka", 1000, 0, 5, 1000, 0, 20);

  dEdXthvp_pr = new TH2F("dEdXthvp_pr", "dEdXthvp_pr", 1000, 0, 5, 1000, 0, 20);

  z_pi = new TH2F("z_pi", "z_pi", 1000, 0, 5, 2000, -20, 20);

  z_ka = new TH2F("z_ka", "z_ka", 1000, 0, 5, 2000, -20, 20);

  z_pr = new TH2F("z_pr", "z_pr", 1000, 0, 5, 2000, -20, 20);

  invbetavpT = new TH2F("invbetavpT", "invbetavpT", 1000, 0, 5, 1000, 0, 5);

  invbetavpT_tof = new TH2F("invbetavpT_tof", "invbetavpT_tof", 1000, 0, 5, 1000, 0, 5);

  normalised_invbetavpT_tof_pi = new TH2F("normalised_invbetavpT_tof_pi", "normalised_invbetavpT_tof_pi", 1000, 0, 5, 2000, -20, 20);

  normalised_invbetavpT_tof_ka = new TH2F("normalised_invbetavpT_tof_ka", "normalised_invbetavpT_tof_ka", 1000, 0, 5, 2000, -20, 20);

  normalised_invbetavpT_tof_pr = new TH2F("normalised_invbetavpT_tof_pr", "normalised_invbetavpT_tof_pr", 1000, 0, 5, 2000, -20, 20);

  mvpT = new TH2F("mvp", "mvp", 2000, -5, 5, 1000, -1.5, 1.5);

  EvP = new TH2F("EvP", "EvP", 1000, 0, 10, 1000, 0, 10);

  invmass = new TH1F("invmass", "invmass", 1000, 0.5, 2.5);

  kfmass = new TH1F("kfmass", "kfmass", 1000, 0.5, 2.5);

  // Jetshape histograms 

  hDeltaEta = new TH1F("hDeltaEta","Relative Eta of D0 and Jets from IV", 100, -1, 1); // Histogram for Delta Eta

  hDeltaPhi = new TH1F("hDeltaPhi","Relative Phi of D0 and Jets from IV", 100, 0, 2); // Histogram for Delta Phi\

  hRad = new TH1F("hRad","Rad of D0 and Jets from IV", 100, 0, 10); // Histogram for dR

  hJetShaped0 = new TH1F("hJetShaped0", "hJetShaped0", numberofbins, 0, R);

  hJetShaped0BgUL = new TH1F("hJetShaped0BgUL", "hJetShaped0BgUL", numberofbins, 0, R);

  hJetShaped0BgLS = new TH1F("hJetShaped0BgLS", "hJetShaped0BgLS", numberofbins, 0, R);

  hJetDistd0 = new TH1F("hJetDistd0", "hJetDistd0", numberofbins, 0, R);

  hJetDistd0BgUL =  new TH1F("hJetDistd0BgUL", "hJetDistd0BgUL", numberofbins, 0, R);

  hJetDistd0BgLS = new TH1F("hJetDistd0BgLS", "hJetDistd0BgLS", numberofbins, 0, R);
}
//
// write histograms
//_____________________________________________________________________________
void StPidInfo::WriteHistograms() {

  const char *event_cuts[10] = {"Total Events", "|V_{z} - V_{z(VPD)}| < 6 cm", "nBEMCMatch > 0", "nBTOFMatch > 0", "|V_{z}| < 60 cm", "", "", "", "", ""};
  const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 20", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

  for (int i=1; i <= 10; i++){
    cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
    cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
    cuthistogram_track->SetBinContent(i, numberoftracks[i-1]);
    cuthistogram_track->GetXaxis()->SetBinLabel(i, track_cuts[i-1]);
  }
  // cuthistogram_event->GetXaxis()->LabelsOption(“v”);
  // cuthistogram_track->GetXaxis()->LabelsOption(“v”);

  cuthistogram_event->Write();
  cuthistogram_track->Write();
  //gPtvpPt->Write();
  dEdXvpT->Write();
  dEdXvp->Write();

  dEdXvpT_pion->Write();
  dEdXvpT_kaon->Write();
  dEdXvpT_proton->Write();
  dEdXvpT_electron->Write();

  dEdXvp_pion->Write();
  dEdXvp_kaon->Write();
  dEdXvp_proton->Write();
  dEdXvp_electron->Write();

  dEdXthvp_pi->Write();
  dEdXthvp_ka->Write();
  dEdXthvp_pr->Write();

  z_pi->Write();
  z_ka->Write();
  z_pr->Write();

  invbetavpT->Write();
  invbetavpT_tof->Write();
  normalised_invbetavpT_tof_pi->Write();
  normalised_invbetavpT_tof_ka->Write();
  normalised_invbetavpT_tof_pr->Write();

  mvpT->Write();
  EvP->Write();

  invmass->Write();
  kfmass->Write();

  hDeltaPhi->Write();
  hDeltaEta->Write();
  hRad->Write();

  hJetDistd0->Write();
  hJetShaped0->Write();

  hJetDistd0BgUL->Write();
  hJetShaped0BgUL->Write();

  hJetDistd0BgLS->Write();
  hJetShaped0BgLS->Write();

}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StPidInfo::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StPidInfo::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

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

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // fill histograms
  // hCentrality->Fill(fCentralityScaled);
  // hMultiplicity->Fill(refCorr2);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  //FillEmcTriggers();

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
  }
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<endl;

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);
  // ======================== end of Triggers ============================= //

  numberofevents[0]++;

  if (abs(mPicoEvent->primaryVertex().z()- mPicoEvent->vzVpd() > 6)) return kStOk;
  numberofevents[1]++;
   
  if (mPicoEvent->nBEMCMatch()==0) return kStOk;
  numberofevents[2]++;

  if (mPicoEvent->nBTOFMatch() == 0) return kStOk;
  numberofevents[3]++;

  if (abs(mPicoEvent->primaryVertex().z() > 60)) return kStOk;
  numberofevents[4]++;

  // =========================== JetMaker =============================== //
  // get JetMaker pointer
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get jet collection associated with JetMaker
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn;
  }

  // // ============================= RhoMaker ============================== //
  // // get RhoMaker pointer from event: old names "StRho_JetsBG"
  // RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  // const char *fRhoMakerNameCh = fRhoMakerName;
  // if(!RhoMaker) {
  //   LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
  //   return kStWarn;
  // }

  // // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  // fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  // if(!fRho) {
  //   LOG_WARN << Form("Couldn't get fRho object! ") << endm;
  //   return kStWarn;    
  // } 
  
  // // get rho/area value from rho object     fRho->ls("");
  // fRhoVal = fRho->GetVal();
  // =======================================================================

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  RunJets();
  //TestTracks();
  ProcessJetForJetShape();

  return kStOK;
}

void StPidInfo::TestTracks()
{

  const Int_t ntracks = mPicoDst->numberOfTracks();

  // for (int i = 0; i < ntracks; i++){

    // cout <<  "Track Number " << i << endl;

    // StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(i));
    // if(!trk1){ continue; }

    // TVector3 mTrk1Org, mTrk1Mom;

    // mTrk1Org = trk1->origin();

    // if(doUsePrimTracks) {
    //   mTrk1Mom = trk1->pMom();
    // } else {
    //   mTrk1Mom = trk1->gMom(mVertex, 0.0);
    // }
    // // cout << "Track exists" << endl;

    // cout << "========================================================================================================" << endl;

    // cout << "Track number = " << i << "\t" << "Track ID = " << trk1->id() << endl;

    // trk1->Print();
    // cout << "========================================================================================================" << endl;

    // cout << "x= "<< mTrk1Org.X() << " y= "<< mTrk1Org.Y() << " z= " << mTrk1Org.Z() << " px= " << mTrk1Mom.X() << " py= " << mTrk1Mom.Y() << " pz= " << mTrk1Mom.Z() << endl;

    // StPicoTrackCovMatrix *cov1 = static_cast<StPicoTrackCovMatrix*>(mPicoDst->trackCovMatrix(i));

    // if (cov1){
    //   cout << "Cov exists." << endl;
    //   cov1->Print();

    //   cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    //   Double_t params[6];
    //   Double_t covariance[21];

    //   ProcessTrackForKF(trk1, cov1, params, covariance);
    // }
    // else continue;

    

    // Bool_t badcov = cov1->isBadCovMatrix();

    // if (!badcov) {cout << "Cov exists." << endl;}

    //IsD0KFParticle(trk1);
    // int tof_loc = trk1->bTofPidTraitsIndex();

    // cout << "TOF LOC is " << tof_loc << endl;

    // if (tof_loc < 0) continue;

    // StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));

    // if (tofpointer){
    //   cout << "TOF track found" << endl;
    //   double invbeta_from_tof = tofpointer->btofBeta();
    //   cout << "Beta is " << invbeta_from_tof << endl;
    // }

  // }

  for (int i = 0; i < ntracks; i++){
    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!trk1){ continue; }
    for (int j = i + 1; j < ntracks; j++){
      if (j >= ntracks) break;
      // cout <<  "Track Number " << i << "\t" << j << endl;
      StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(j));
      if(!trk2){ continue; }

      StPicoTrackCovMatrix *cov1 = static_cast<StPicoTrackCovMatrix*>(mPicoDst->trackCovMatrix(i));
      StPicoTrackCovMatrix *cov2 = static_cast<StPicoTrackCovMatrix*>(mPicoDst->trackCovMatrix(j));

      bool D0, D0BgUnlike, D0BgLike;
      double mass;
      TVector3 mom3;

      


      IsD0(trk1, trk2, D0, D0BgUnlike, D0BgLike, mass, mom3);
      // if (D0) cout << "D0" << "\t" << mass << endl; 
      // else if (D0BgUnlike) cout << "D0BgUnLike" << "\t" << mass << endl; 
      // else if (D0BgLike) cout << "D0BgLike" << "\t" << mass << endl;
      // else cout << "None" << "\t" << mass << endl;
      IsD0KFParticle(trk1, trk2, cov1, cov2);



      // cout << "===================================================================================" << endl;
    }
  }

  // StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(1));

  // cout << trk1->charge() << endl;

  // KFParticle d1, d2, d3, d4, d5; 
  
  // Double_t params1[6] = {0.0,0.0,0.0,5.0,0.0,1.0};
  // Double_t params2[6] = {0.0,0.0,0.0,0.0,0.0,1.0};
  // Double_t covar1[21] = {0};


  // d1.Create(params1, covar1, 1, 321);
  // d2.Create(params1, covar1, 1, -321);
  // d3.Create(params1, covar1, -1, 321);
  // d4.Create(params2, covar1, -1, -211);

  // d5 = KFParticle(d1,d4);

  // cout << d1.GetMass() << "\t" << d2.GetMass() << "\t" << d3.GetMass() << "\t" << d4.GetMass() << "\t" << d5.GetMass() << endl;
}

//________________________________________________________________________
void StPidInfo::RunJets()
{

  if(fCorrJetPt) {
    fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
  } else {
    fLeadingJet = GetLeadingJet(fJetMakerName);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
  }

  
  Int_t njets = fJets->GetEntries();
  //if (njets!=0) cout << njets << endl;

  Int_t d0CandidateJet, d0BgCandidateULJet, d0BgCandidateLSJet;

  for (int ijet = 0; ijet < njets; ijet++){

    Bool_t D0Candidate = kFALSE;
    Bool_t D0BgCandidateUL = kFALSE;
    Bool_t D0BgCandidateLS = kFALSE;

    // get jet pointer
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    //if (abs(jet->Eta()) > 1.) continue;
    if ((jet->Pt() < 20.)) continue;

    numberoftracks[0]+=jet->GetNumberOfTracks();

    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
      if(!trk){ cout << "No track" << endl; continue; }
      FillPidHistograms(trk);
    }
    //cout << ijet << "\t" << jet->GetNumberOfTracks() << endl;

    // for(int itrk1 = 0; itrk1 < jet->GetNumberOfTracks(); itrk1++) {
    //   int trackid1 = jet->TrackAt(itrk1);

    //   StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(trackid1));
    //   if(!trk1){ cout << "No track" << endl; continue; }

    //   for(int itrk2 = itrk1 + 1; itrk2 < jet->GetNumberOfTracks(); itrk2++) {
    //     int trackid2 = jet->TrackAt(itrk2);

    //     StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(trackid2));
    //     if(!trk2){ cout << "No track" << endl; continue; }

    //     //cout << "Track 1 = " << itrk1  << "\t" << "Track 2 = " << itrk2 << endl;

        
    //     bool D0, D0BgUnlike, D0BgLike;
    //     double mass;
    //     TVector3 mom3;

    //     //cout << "Step before IsDo" << endl;

    //     IsD0(trk1, trk2, D0, D0BgUnlike, D0BgLike, mass, mom3);

    //     if (D0) D0Candidate = kTRUE;
    //     if (D0BgUnlike) D0BgCandidateUL = kTRUE;
    //     if (D0BgLike) D0BgCandidateLS = kTRUE;

    //     if ((D0Candidate) && (D0BgCandidateUL) && (D0BgCandidateLS)) break;
    //   }
    //   if ((D0Candidate) && (D0BgCandidateUL) && (D0BgCandidateLS)) break;
    // }

    // if (D0Candidate){ d0CandidateID.push_back(ijet); cout << "Found a D0 Jet" << endl;}
    // if (D0BgCandidateUL){ d0BgCandidateULID.push_back(ijet); cout << "Found a UL Jet" << endl;}
    // if (D0BgCandidateLS){ d0BgCandidateLSID.push_back(ijet); cout << "Found a LS Jet" << endl;}
  }
}

Bool_t StPidInfo::IsAnAcceptableTrack(StPicoTrack *trk){
  if (trk->nHitsDedx() < 20) return kFALSE;
  if (float(trk->nHitsDedx())/float(trk->nHitsMax()) < 0.52) return kFALSE;
  if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 3.) return kFALSE; 

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, 0.0);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  double px = mTrkMom.x();
  double py = mTrkMom.y();
  double pz = mTrkMom.z();
  double p = mTrkMom.Mag();
  short charge = trk->charge();

  // bool bemctrack = trk->isBemcMatchedTrack();
  // if ((!bemctrack)) return kFALSE;

  bool toftrack = trk->isTofTrack();
  if ((!toftrack) && (pt < 1.6)) return kFALSE;
  if (charge == 0) return kFALSE;

  return kTRUE;
}

void StPidInfo::FillPidHistograms(StPicoTrack *trk){
  if (!IsAnAcceptableTrack(trk)) return;

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, 0.0);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  double px = mTrkMom.x();
  double py = mTrkMom.y();
  double pz = mTrkMom.z();
  double p = mTrkMom.Mag();
  short charge = trk->charge();

  // dedx info
  double dedx = trk->dEdx();
  double dedxresolution = trk->dEdxError();

  // bichsel function approximation
  double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
  double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
  double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

  //z - variables
  double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
  double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
  double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

  dEdXvpT->Fill(pt, dedx);
  z_pi->Fill(pt, zpi);
  z_ka->Fill(pt, zka);
  z_pr->Fill(pt, zpr);

  // int tower_loc = trk->bemcTowerIndex();
  int tof_loc = trk->bTofPidTraitsIndex();

  if (tof_loc >= 0)
  {
    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));

    if (tofpointer){
      // cout << "TOF track found" << endl;
      double invbeta_from_tof = tofpointer->btofBeta();
      if (invbeta_from_tof > 0.0001) {
        // cout << "Beta " << invbeta_from_tof << "\n"; 
        invbeta_from_tof = 1/invbeta_from_tof;

        invbetavpT_tof->Fill(pt*charge, invbeta_from_tof);

        double energy = invbeta_from_tof*p;
        double masssq = (pow(invbeta_from_tof, 2) - 1)*pow(p,2);

        mvpT->Fill(p*charge, masssq);

        double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
        double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
        double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

        double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.012;
        double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.012;
        double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr)/0.012;

        normalised_invbetavpT_tof_pi->Fill(pt, normalisedinvbeta_for_pi);
        normalised_invbetavpT_tof_ka->Fill(pt, normalisedinvbeta_for_ka);
        normalised_invbetavpT_tof_pr->Fill(pt, normalisedinvbeta_for_pr);
      }
    }
  }
}

void StPidInfo::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){

  if(!IsAnAcceptableTrack(trk)){pid = 0; m = 0.; e = 0.; return;}
  else{

    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      mTrkMom = trk->pMom();
    } else {
      mTrkMom = trk->gMom(mVertex, 0.0);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();
    double px = mTrkMom.x();
    double py = mTrkMom.y();
    double pz = mTrkMom.z();
    double p = mTrkMom.Mag();
    short charge = trk->charge();

    double energy = 0.;
    double mass = 0.;

    // int tower_loc = trk->bemcTowerIndex();
    int tof_loc = trk->bTofPidTraitsIndex();

    if (tof_loc >= 0)
    {
      StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
      // if (tof) {cout << "TOF track found" << endl;}
      // tof matched tracks
      if (tofpointer){
        // cout << "TOF track found" << endl;
        double invbeta_from_tof = tofpointer->btofBeta();
        if (invbeta_from_tof > 0.0001) {
          // cout << "Beta " << invbeta_from_tof << "\n"; 
          invbeta_from_tof = 1/invbeta_from_tof;
          energy = invbeta_from_tof*p;
          if (pow(energy, 2) - pow(p, 2) < 0) {pid = 0; m = 0.; e = 0; return;}
          mass = TMath::Sqrt(pow(invbeta_from_tof, 2) - 1)*p;

          double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
          double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
          double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

          double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.012;
          double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.012;
          double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr)/0.012;

          double f_res = 0.884 + 0.0174/pow((p + 0.0839), 4.23);
          double f_pos = 0.0316 + 0.00137/pow((p + 0.101), 6.89);

          if ((normalisedinvbeta_for_ka > -2*f_res + f_pos) && (normalisedinvbeta_for_ka < 2*f_res + f_pos) && (charge == 1)) {pid = 2; m = mass; e = energy; return;}
          else if ((normalisedinvbeta_for_ka > -2*f_res + f_pos) && (normalisedinvbeta_for_ka < 2*f_res + f_pos) && (charge == -1)) {pid = -2; m = mass; e = energy; return;}
          else if ((normalisedinvbeta_for_pi > -1.9) && (normalisedinvbeta_for_pi < 2.1) && (charge == 1)) {pid = 1; m = mass; e = energy; return;}
          else if ((normalisedinvbeta_for_pi > -1.9) && (normalisedinvbeta_for_pi < 2.1) && (charge == -1)) {pid = -1; m = mass; e = energy; return;}
          else if ((normalisedinvbeta_for_pr > -1.9) && (normalisedinvbeta_for_pr < 2.1) && (charge == 1)) {pid = 3; m = mass; e = energy; return;}
          else if ((normalisedinvbeta_for_pr > -1.9) && (normalisedinvbeta_for_pr < 2.1) && (charge == -1)) {pid = -3; m = mass; e = energy; return;}
          else {pid = 0; m = 0.; e = 0; return;}
        }

        else {pid = 0; m = 0.; e = 0; return;}
      }

      else {pid = 0; m = 0.; e = 0; return;}
    }
    else {
      if (pt > 1.6){
        // dedx info
        double dedx = trk->dEdx();
        double dedxresolution = trk->dEdxError();

        // bichsel function approximation
        double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
        double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
        double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

        //z - variables
        double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
        double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
        double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

        if ((zpi > -2) && (zpi < 2)){
          energy = TMath::Sqrt(pow(p,2) + pow(Mpion, 2));
          if (charge == 1) {pid = 1; m = Mpion; e = energy; return;}
          if (charge == -1) {pid = -1; m = Mpion; e = energy; return;}
        }

        if ((zka > -2) && (zka < 2)){
          energy = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2));
          if (charge == 1) {pid = 2; m = Mkaon; e = energy; return;}
          if (charge == -1) {pid = -2; m = Mkaon; e = energy; return;}
        }

        if ((zpr > -2) && (zpr < 2)){
          energy = TMath::Sqrt(pow(p,2) + pow(Mproton, 2));
          if (charge == 1) {pid = 3; m = Mproton; e = energy; return;}
          if (charge == -1) {pid = -3; m = Mproton; e = energy; return;}
        }

        else {pid = 0; m = 0.; e = 0; return;}
      }

      else {pid = 0; m = 0.; e = 0; return;}
    }
    // if (tofpointer) delete tofpointer;
  }
}

void StPidInfo::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum){

  int particle1, particle2;
  double mass1, mass2;
  double energy1, energy2;

  // cout << "IsWhatParticle called" << endl;

  IsWhatParticle(trk1, particle1, mass1, energy1);
  IsWhatParticle(trk2, particle2, mass2, energy2);

  TVector3 mTrk1Mom, mTrk2Mom;

  if ((mass1<0.000001) || (mass2<0.000001) || (energy1<0.000001) || (energy2<0.000001)){
    particle=0;
    invmass=0.;
    momentum = mTrk1Mom + mTrk2Mom;
    return;
  }
  else{
    if(doUsePrimTracks) {
      mTrk1Mom = trk1->pMom();
    } else {
      mTrk1Mom = trk1->gMom(mVertex, 0.0);
    }

    if(doUsePrimTracks) {
      mTrk2Mom = trk2->pMom();
    } else {
      mTrk2Mom = trk2->gMom(mVertex, 0.0);
    }

    particle = particle1*particle2;
    invmass = TMath::Sqrt(pow(mass1,2) + pow(mass2,2) + 2*(energy1*energy2 - mTrk1Mom*mTrk2Mom));
    momentum = mTrk1Mom + mTrk2Mom;
  }
}

void StPidInfo::IsD0(StPicoTrack *trk1, StPicoTrack *trk2, bool &D0, bool &D0BgUnlike, bool &D0BgLike, double &mass, TVector3 &mom3){

  int particle;

  InvariantMass(trk1, trk2, particle, mass, mom3);

  D0 = kFALSE; D0BgUnlike = kFALSE; D0BgLike = kFALSE;

  if(particle == -2){
    // cout << "===================================================================================" << endl;

    // cout << "InvariantMass is " << mass << endl;

    invmass->Fill(mass);

    if ((mass > fInvMassSignal1) && (mass < fInvMassSignal2)){
      D0 = kTRUE;
      // cout << "InvariantMass is " << mass << endl;
    }
    if (((mass > fInvMassULBg1) && (mass < fInvMassSignal1)) || ((mass > fInvMassSignal2) && (mass < fInvMassULBg2))){
      D0BgUnlike = kTRUE;
      // cout << "InvariantMass is " << mass << endl;
    }
  }

  if(particle == 2){
    // cout << "===================================================================================" << endl;

    // cout << "InvariantMass is " << mass << endl;
    if ((mass > fInvMassLSBg1) && (mass < fInvMassLSBg2)){
      D0BgLike = kTRUE;
      // cout << "InvariantMass is " << mass << endl;
    }
  }
}

void StPidInfo::ProcessTrackForKF(StPicoTrack *trk, StPicoTrackCovMatrix *picocov, Double_t params[6], Double_t cov[21]){
  const Float_t *paramfromcov = picocov->params();
  const Float_t *sigmafromcov = picocov->sigmas();
  const Float_t *corrfromcov = picocov->correlations();

  static StDcaGeometry dcaGeometryPi;

  Float_t err[15] = {sigmafromcov[0], 
                     corrfromcov[0], sigmafromcov[1],
                     corrfromcov[1], corrfromcov[2], sigmafromcov[2],
                     corrfromcov[3], corrfromcov[4], corrfromcov[5], sigmafromcov[3],
                     corrfromcov[6], corrfromcov[7], corrfromcov[8], corrfromcov[9], sigmafromcov[4]};

  dcaGeometryPi.set(paramfromcov, err);
  dcaGeometryPi.GetXYZ(params, cov);
}


void StPidInfo::IsD0KFParticle(StPicoTrack *trk1, StPicoTrack *trk2, StPicoTrackCovMatrix *cov1, StPicoTrackCovMatrix *cov2){
  Double_t params1[6], params2[6];
  Double_t covar1[21], covar2[21];

  int charge1 = trk1->charge();
  int charge2 = trk2->charge();

  Int_t pid1, pid2;

  Double_t m1, m2, e1, e2;

  IsWhatParticle(trk1, pid1, m1, e1);
  IsWhatParticle(trk2, pid2, m2, e2);

  Int_t particle = pid1*pid2;

  if (particle!=-2){return;} 

  ProcessTrackForKF(trk1, cov1, params1, covar1);
  ProcessTrackForKF(trk2, cov2, params2, covar2);

  Int_t pdg1, pdg2;

  switch(pid1){
    case 1: pdg1 = 211;     break;
    case 2: pdg1 = 321;     break;
    // case 3: pdg1 = 2212;    break;
    case -1:  pdg1 = -211;   break;
    case -2:  pdg1 = -321;   break;
    // case -3:  pdg1 = -2212;  break;
  }

  switch(pid2){
    case 1: pdg2 = 211;     break;
    case 2: pdg2 = 321;     break;
    // case 3: pdg2 = 2212;    break;
    case -1:  pdg2 = -211;   break;
    case -2:  pdg2 = -321;   break;
    // case -3:  pdg2 = -2212;  break;
  }

  KFParticle d1, d2; 
  

  d1.Create(params1, covar1, charge1, pdg1);
  d2.Create(params2, covar2, charge2, pdg2);

  KFParticle mother = KFParticle(d1, d2);

  double mass1, mass1err;
  d1.GetMass(mass1, mass1err);
  double mass2, mass2err;
  d2.GetMass(mass2, mass2err);

  double mass_mom, mass_mom_err;
  mother.GetMass(mass_mom, mass_mom_err);
  //mother.SetMassConstraint(1.865);

  kfmass->Fill(mass_mom);
}


void StPidInfo::ProcessJetForJetShape(){

  Int_t njets = fJets->GetEntries();

  float diffd0dist[numberofbins];
  float diffd0distBgUL[numberofbins];
  float diffd0distBgLS[numberofbins];

  float diffd0shape[numberofbins];
  float diffd0shapeBgUL[numberofbins];
  float diffd0shapeBgLS[numberofbins];

  // Jets which have a D0 candidate

  //for (Int_t iJetD0 = 0; iJetD0 = d0CandidateID.size(); iJetD0++){
  for (int ijet = 0; ijet < njets; ijet++){
    // get jet pointer
    //StJet *jet = static_cast<StJet*>(fJets->At(d0CandidateID[iJetD0]));
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;
    numberoftracks[0]+=jet->GetNumberOfTracks();

    double jetpt = jet->Pt();
    double jeteta = jet->Eta();
    double jetphi = jet->Phi();

    if ((jetpt < 20.)) continue;

    for (int bin = 0; bin < numberofbins; bin++){
      diffd0dist[bin] = 0;
      diffd0distBgUL[bin] = 0;
      diffd0distBgLS[bin] = 0;

      diffd0shape[bin] = 0;
      diffd0shapeBgUL[bin] = 0;
      diffd0shapeBgLS[bin] = 0;
    }

    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {

      int trackid = jet->TrackAt(itrk);

      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
      if(!trk){ cout << "No track" << endl; continue; }

      TVector3 mTrkMom;
      if(doUsePrimTracks) {
        mTrkMom = trk->pMom();
      } else {
        mTrkMom = trk->gMom(mVertex, 0.0);
      }

      // track variables
      double pt = mTrkMom.Perp();
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      double px = mTrkMom.x();
      double py = mTrkMom.y();
      double pz = mTrkMom.z();
      double p = mTrkMom.Mag();
      short charge = trk->charge();

      double delphiD = dPhi(jetphi, phi);
      double deletaD = dEta(jeteta, eta);
      double delRD = dR(delphiD, deletaD);

      hDeltaPhi->Fill(delphiD);
      hDeltaEta->Fill(deletaD);
      hRad->Fill(delRD);

      for (int bin = 0; bin < numberofbins; bin++){
        if (delRD >= bin*deltar && delRD < (bin+1)*deltar){
          diffd0dist[bin]++;
          diffd0shape[bin] += pt/jetpt;
        }
      } // end of bin loop
    } // end of track loop

    for (int bin = 0; bin < numberofbins; bin++){
      hJetDistd0->Fill((2*bin+1)*deltar/2, diffd0dist[bin]);
      hJetShaped0->Fill((2*bin+1)*deltar/2, diffd0shape[bin]);
    }

    // cout << "Jet # " << ijet << "#oftracks " << jet->GetNumberOfTracks() << endl;

    // for (int bin = 0; bin < numberofbins; bin++){
    //   cout << "Bin# " << bin << "\t" << diffd0dist[bin] << "\t" << diffd0shape[bin] << endl;
    // }
  } //end of jet loop

  // // Jets which have a D0 Bg from UL 

  // for (Int_t iJetD0 = 0; iJetD0 = d0BgCandidateULID.size(); iJetD0++){
  //   // get jet pointer
  //   StJet *jet = static_cast<StJet*>(fJets->At(d0BgCandidateULID[iJetD0]));
  //   if(!jet) continue;
  //   numberoftracks[0]+=jet->GetNumberOfTracks();

  //   double jetpt = jet->Pt();
  //   double jeteta = jet->Eta();
  //   double jetphi = jet->Phi();

  //   for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {

  //     int trackid = jet->TrackAt(itrk);

  //     StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
  //     if(!trk){ cout << "No track" << endl; continue; }

  //     TVector3 mTrkMom;
  //     if(doUsePrimTracks) {
  //       mTrkMom = trk->pMom();
  //     } else {
  //       mTrkMom = trk->gMom(mVertex, 0.0);
  //     }

  //     // track variables
  //     double pt = mTrkMom.Perp();
  //     double phi = mTrkMom.Phi();
  //     double eta = mTrkMom.PseudoRapidity();
  //     double px = mTrkMom.x();
  //     double py = mTrkMom.y();
  //     double pz = mTrkMom.z();
  //     double p = mTrkMom.Mag();
  //     short charge = trk->charge();

  //     double delphiD = dPhi(jetphi, phi);
  //     double deletaD = dEta(jeteta, eta);
  //     double delRD = dR(delphiD, deletaD);

  //     for (int bin = 0; bin < numberofbins; bin++){
  //       if (delRD >= bin*deltar && delRD < (bin+1)*deltar){
  //         diffd0distBgUL[bin]++;
  //         diffd0shapeBgUL[bin]+= pt/jetpt;
  //       }
  //     } // end of bin loop
  //   } // end of track loop
  // } //end of jet loop

  // // Jets which have a D0 Bg from LS

  // for (Int_t iJetD0 = 0; iJetD0 = d0BgCandidateLSID.size(); iJetD0++){
  //   // get jet pointer
  //   StJet *jet = static_cast<StJet*>(fJets->At(d0BgCandidateLSID[iJetD0]));
  //   if(!jet) continue;
  //   numberoftracks[0]+=jet->GetNumberOfTracks();

  //   double jetpt = jet->Pt();
  //   double jeteta = jet->Eta();
  //   double jetphi = jet->Phi();

  //   for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {

  //     int trackid = jet->TrackAt(itrk);

  //     StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
  //     if(!trk){ cout << "No track" << endl; continue; }

  //     TVector3 mTrkMom;
  //     if(doUsePrimTracks) {
  //       mTrkMom = trk->pMom();
  //     } else {
  //       mTrkMom = trk->gMom(mVertex, 0.0);
  //     }

  //     // track variables
  //     double pt = mTrkMom.Perp();
  //     double phi = mTrkMom.Phi();
  //     double eta = mTrkMom.PseudoRapidity();
  //     double px = mTrkMom.x();
  //     double py = mTrkMom.y();
  //     double pz = mTrkMom.z();
  //     double p = mTrkMom.Mag();
  //     short charge = trk->charge();

  //     double delphiD = dPhi(jetphi, phi);
  //     double deletaD = dEta(jeteta, eta);
  //     double delRD = dR(delphiD, deletaD);

  //     for (int bin = 0; bin < numberofbins; bin++){
  //       if (delRD >= bin*deltar && delRD < (bin+1)*deltar){
  //         diffd0distBgLS[bin]++;
  //         diffd0shapeBgLS[bin] += pt/jetpt;
  //       }
  //     } // end of bin loop
  //   } // end of track loop
  // } //end of jet loop

  // for (int bin = 0; bin < numberofbins; bin++){
  //   hJetDistd0->Fill((2*bin+1)*deltar/2, diffd0dist[bin]);
  //   hJetShaped0->Fill((2*bin+1)*deltar/2, diffd0shape[bin]);

  //   hJetDistd0BgUL->Fill((2*bin+1)*deltar/2, diffd0distBgUL[bin]);
  //   hJetShaped0BgUL->Fill((2*bin+1)*deltar/2, diffd0shapeBgUL[bin]);

  //   hJetDistd0BgLS->Fill((2*bin+1)*deltar/2, diffd0distBgLS[bin]);
  //   hJetShaped0BgLS->Fill((2*bin+1)*deltar/2, diffd0shapeBgLS[bin]);
  //  }
} //end of function

// jet function

// void StPidInfo::WriteInfo(){

//   std::ofstream normalisedinvbetaVp_pi;
//   normalisedinvbetaVp_pi.open("normalisedinvbetaVp_pi.csv", ios::out | ios::app | ios::binary);
//   normalisedinvbetaVp_pi << track_p << "\t" << normalisedinvbeta_for_pi << endl;
//   normalisedinvbetaVp_pi.close();

//   std::ofstream normalisedinvbetaVp_ka;
//   normalisedinvbetaVp_ka.open("normalisedinvbetaVp_ka.csv", ios::out | ios::app | ios::binary);
//   normalisedinvbetaVp_ka << track_p << "\t" << normalisedinvbeta_for_ka << endl;
//   normalisedinvbetaVp_ka.close();

// }
