// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StTagD0TPCEvents.h"
#include "StMemStat.h"
#include "phys_constants.h"
#include <limits>
#include "math.h"
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
#include "TLorentzVector.h"


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
// #include "StDcaGeometry.h"
#include "StPhysicalHelix.hh"
#include "StThreeVectorF.hh"


// #include "StarClassLibrary/StRedefinedHelix.hh"


// // KFParticle includes
// #include "StRoot/StEvent/StDcaGeometry.h"
// #include "StRoot/StarRoot/KFParticle.h"

#include "StBTofUtil/tofPathLength.hh"

// Bichsel includes
#include "StBichsel/Bichsel.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StTagD0TPCEvents)

//________________________________________________________________________
StTagD0TPCEvents::StTagD0TPCEvents(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTagD0TPCEvents::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fMCEventsWithoutCent = kFALSE;
  fCentralitySelectionCut = -99;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fTopoLevel = 1;
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
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  for(int i=0; i<14; i++) {numberofevents[i] = 0.; numberoftracks[i] = 0.;}
  //Mass cut
  fInvMassSignal1 = 1.80;
  fInvMassSignal2 = 1.92;

  fInvMassULBg1 = 1.70;
  fInvMassULBg2 = 2.08;

  fInvMassLSBg1 = 1.70;
  fInvMassLSBg2 = 2.08;

  fd0 = kFALSE;
  fd0BgUS = kFALSE;
  fd0BgLS = kFALSE;

  // fd0TrackIndices.clear();
  // fd0BgUSTrackIndices.clear();
  // fd0BgLSTrackIndices.clear();
}

//
//________________________________________________________________________
StTagD0TPCEvents::~StTagD0TPCEvents()
{ 

  if (hEventZvertex_whole) delete hEventZvertex_whole;
  if (hEventZvertex_VPD) delete hEventZvertex_VPD;
  if (hEventZvertex_diff) delete hEventZvertex_diff;

  if (hCentrality) delete hCentrality;
  if (hMultiplicity) delete hMultiplicity;
  if (cuthistogram_event) delete cuthistogram_event;

  for (int i = 0; i < 11; i++){
    if (hAcceptedPt[i]) delete hAcceptedPt[i];
    if (hAcceptedEta[i]) delete hAcceptedEta[i];
  }

  if (z_pi) delete z_pi;
  if (z_ka) delete z_ka;
  if (normalised_invbetavpT_tof_pi) delete normalised_invbetavpT_tof_pi;
  if (normalised_invbetavpT_tof_ka) delete normalised_invbetavpT_tof_ka;
  if (dEdXvp) delete dEdXvp;
  if (invbetavp) delete invbetavp;

  if (hD0CentPtEtaMDphiDaug) delete hD0CentPtEtaMDphiDaug;
  if (hD0CentPtEtaMDphiDaugLikeSign) delete hD0CentPtEtaMDphiDaugLikeSign;

  if (hNumberOfD0s) delete hNumberOfD0s;
  if (hNumberOfD0BgUS) delete hNumberOfD0BgUS;
  if (hNumberOfD0BgLS) delete hNumberOfD0BgLS;


  if(mEmcPosition) delete mEmcPosition;

}

//________________________________________________________________________
Int_t StTagD0TPCEvents::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTagD0TPCEvents::Finish() { 

  cout << "StTagD0TPCEvents::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),  "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StTagD0TPCEvents::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StTagD0TPCEvents::DeclareHistograms() {

  //gPtvpPt = new TH2F("gPtvpPt", "gPtvpPt", 1000, 0, 20, 1000, 0, 20);

  // binning for cent histograms
  int nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100.;
  }

  hEventZvertex_whole = new TH1F("hEventZvertex_whole", "z-vertex distribution for all events", 200, -100., 100.);
  hEventZvertex_VPD = new TH1F("hEventZvertex_VPD", "z-vertex distribution from VPD", 200, -100., 100.);
  hEventZvertex_diff = new TH1F("hEventZvertex_diff", "z-vertex distribution Difference", 200, -100., 100.);

  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);

  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  cuthistogram_event = new TH1F("cuthistogram_event", "cuthistogram_event", 14, 0, 14);

  // Pt After Track Cuts

  for (int i = 0; i < 11; i++){
    hAcceptedPt[i] = new TH1F(Form("hAcceptedPt_%i", i), Form("hAcceptedPt_%i", i), 1000, 0, 30);
    hAcceptedEta[i] = new TH1F(Form("hAcceptedEta%i", i), Form("hAcceptedEta%i", i), 1000, -2, 2);
  }
  

  z_pi = new TH2F("z_pi", "z_pi", 3000, 0, 30, 2000, -10, 10);
  z_ka = new TH2F("z_ka", "z_ka", 3000, 0, 30, 2000, -10, 10);
  normalised_invbetavpT_tof_pi = new TH2F("normalised_invbetavpT_tof_pi", "normalised_invbetavpT_tof_pi", 3000, 0, 30, 2000, -10, 10);
  normalised_invbetavpT_tof_ka = new TH2F("normalised_invbetavpT_tof_ka", "normalised_invbetavpT_tof_ka", 3000, 0, 30, 2000, -10, 10);

  dEdXvp = new TH2F("dEdXvp", "dEdXvp", 3000, 0, 30, 4000, -20, 20);
  invbetavp = new TH2F("invbetavp", "invbetavp", 3000, 0, 30, 4000, -20, 20);

  const int nDimDaug = 6;
  int nBinsDaug[nDimDaug] = {9, 20, 5, 100, 5, 4};//cent, pt*charge, daughterpt1, m, daughterpt2*charge, rapidity
  double xMinDaug[nDimDaug] = {0, -10, 0, 1.7, 0, -1.0};
  double xMaxDaug[nDimDaug] = {9, 10, 10, 2.1, 10, 1.0};

  hD0CentPtEtaMDphiDaug = new THnF("hD0CentPtEtaMDphiDaug", "hD0CentPtEtaMDphiDaug", nDimDaug, nBinsDaug, xMinDaug, xMaxDaug);
  hD0CentPtEtaMDphiDaugLikeSign = new THnF("hD0CentPtEtaMDphiDaugLikeSign", "hD0CentPtEtaMDphiDaugLikeSign", nDimDaug, nBinsDaug, xMinDaug, xMaxDaug);

  hNumberOfD0s = new TH1F("hNumberOfD0s", "hNumberOfD0s", 21, -0.5, 20.5);
  hNumberOfD0BgUS = new TH1F("hNumberOfD0BgUS", "hNumberOfD0BgUS", 21, -0.5, 20.5);
  hNumberOfD0BgLS = new TH1F("hNumberOfD0BgLS", "hNumberOfD0BgLS", 21, -0.5, 20.5);

}
//
// write histograms
//_____________________________________________________________________________
void StTagD0TPCEvents::WriteHistograms() {

  // Centrality Histograms

  hEventZvertex_whole->Write();
  hEventZvertex_VPD->Write();
  hEventZvertex_diff->Write();

  hCentrality->Write();
  hMultiplicity->Write();

  // Event Cuts Histograms
  const char *event_cuts[14] = {"Total Events", "Good Runs", "Max Track Pt < 30", "|V_{i}| != 0", "|V_{z}| < 30 cm", "MinBias","|V_{z}| < 6 cm", "|V_{r}| < 2 cm ", "|V_{z} - V_{z(VPD)}| < 3 cm", "nBEMCMatch > 0", "nBTOFMatch > 0", "HT1", "HT2", "HT3"};
  const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 20", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

  for (int i=1; i <= 14; i++){
    cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
    cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
    // cuthistogram_track->SetBinContent(i, numberoftracks[i-1]);
    // cuthistogram_track->GetXaxis()->SetBinLabel(i, track_cuts[i-1]);
  }
  cuthistogram_event->Write();
  // cuthistogram_track->Write();

  // Track Cut Histograms
  for (int i = 0; i < 11; i++){
    hAcceptedPt[i]->Write();
    hAcceptedEta[i]->Write();
  }

  // Particle Identification Histograms
  z_pi->Write();
  z_ka->Write();
  normalised_invbetavpT_tof_pi->Write();
  normalised_invbetavpT_tof_ka->Write();

  dEdXvp->Write();
  invbetavp->Write();

  hD0CentPtEtaMDphiDaug->Write();
  hD0CentPtEtaMDphiDaugLikeSign->Write();

  // Pt and Invmass Histograms for KPi pairs after topological cuts

  hNumberOfD0s->Write();
  hNumberOfD0BgUS->Write();
  hNumberOfD0BgLS->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTagD0TPCEvents::Clear(Option_t *opt) {
  // fJets->Clear();
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTagD0TPCEvents::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  fd0 = kFALSE;
  fd0BgUS = kFALSE;
  fd0BgLS = kFALSE;

  //zero out global vectors
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();

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

  numberofevents[0]++;

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  numberofevents[1]++;

  numberofevents[2]++;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // cout << "B = " << Bfield << endl;
  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  // cout << "Vertex" << mVertex.X() << "\t" << mVertex.Y() << "\t" << mVertex.Z() << endl;

  // if (mVertex.x() == 0 || mVertex.y() == 0 || mVertex.z() == 0) return kStOK;

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5) return kStOK;


  numberofevents[3]++;

  double zVtx_VPD = mPicoEvent->vzVpd();

  hEventZvertex_whole->Fill(zVtx);
  hEventZvertex_VPD->Fill(zVtx_VPD);
  hEventZvertex_diff->Fill(zVtx - zVtx_VPD);

  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  numberofevents[4]++;


  // HFT Analysis Specific Cuts For Event

  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CentMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int cent9 = mCentMaker->GetCent9();
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();

  const int nCent = 5;
  const float CentEdge[nCent+1] = { -0.5, 1.5, 3.5, 5.5, 6.5, 8.5 };

  int cent9bin = -99;

  // cout << "The Centrality bin is : " << grefMult << "\t" << cent16 << endl;

  for (int i = 0; i < nCent; i++)
  {
    if ((cent9 >= CentEdge[i]) && (cent9 < CentEdge[i + 1]))
    cent9bin = i;
  }

  // cout << "The Centrality bin is : "  << fCentralityScaled << "\t" << cent9 << "\t" << GetFiveCentBin(fCentralityScaled) << "\t" << cent9bin  << endl;

  // cout << "Event Weight : " << mCentMaker->Weight() << endl;
  //cout << "The centrality bin is " << fCentralityScaled << endl;
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%

  if (!fMCEventsWithoutCent)
  {
    if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

    // cout << "Here" << endl;

    // fill histograms
    hCentrality->Fill(fCentralityScaled);
    hMultiplicity->Fill(refCorr2);

  

    // cut on centrality for analysis before doing anything
    if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

    // ============================ end of CENTRALITY ============================== //

    // ========================= Trigger Info =============================== //
    // looking at the EMCal triggers - used for QA and deciding on HT triggers
    //  FillEmcTriggers();

    // get trigger IDs from PicoEvent class and loop over them
    vector<unsigned int> mytriggers = mPicoEvent->triggerIds();

    int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
    int arrHT1_Run14[] = {450201, 450211, 460201};
    int arrHT2_Run14[] = {450202, 450212, 460202, 460212};
    int arrHT3_Run14[] = {450203, 450213, 460203};

    if (!doppAnalysis){

      bool matchMB = kFALSE;

      for(int i = 0; i < sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14); i++) {
        if(mPicoEvent->isTrigger(arrMB5_Run14[i])) matchMB = kTRUE;
        if(matchMB) break;
      }

      if (!matchMB) return kStOk;

    }

    numberofevents[5]++;
    // ======================== end of Triggers ============================= //

    

    if (abs(zVtx) > 6.) return kStOk;
    numberofevents[6]++;

    // if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.) return kStOK;
    if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) <= 2.) numberofevents[7]++;

    if (abs(zVtx - zVtx_VPD) > 3) return kStOk;
    numberofevents[8]++;
     
    if (mPicoEvent->nBEMCMatch() != 0) numberofevents[9]++;

    if (mPicoEvent->nBTOFMatch() != 0) numberofevents[10]++;

    if (!doppAnalysis){

      bool matchHT1 = kFALSE;
      bool matchHT2 = kFALSE;
      bool matchHT3 = kFALSE;

      for(int i = 0; i < sizeof(arrHT1_Run14)/sizeof(*arrHT1_Run14); i++) {
        if(mPicoEvent->isTrigger(arrHT1_Run14[i])) matchHT1 = kTRUE;
        if(matchHT1) break;
      }

      for(int i = 0; i < sizeof(arrHT2_Run14)/sizeof(*arrHT2_Run14); i++) {
        if(mPicoEvent->isTrigger(arrHT2_Run14[i])) matchHT2 = kTRUE;
        if(matchHT2) break;
      }

      for(int i = 0; i < sizeof(arrHT3_Run14)/sizeof(*arrHT3_Run14); i++) {
        if(mPicoEvent->isTrigger(arrHT3_Run14[i])) matchHT3 = kTRUE;
        if(matchHT3) break;
      }

      if (matchHT1 || matchHT2 || matchHT3) numberofevents[11]++;


      if (matchHT2 || matchHT3) numberofevents[12]++;

     
      if (matchHT3) numberofevents[13]++;
    // cout << "MinBias Event Accepted." << endl;
    }
  }

  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  if (fDebugLevel == 1) cout << "Number of tracks found = " << ntracks << "\t" << nglobaltracks << endl;

  // cout << Mpion << "\t" << M_PION_PLUS << "\t" << Mkaon << "\t" << M_KAON_PLUS << endl;
  RunTracks();
  // ProcessJetForJetShape();
  if (fDebugLevel == 1) TestTracks();
  return kStOK;
}


//__________________________________________________________________________________________
  
void StTagD0TPCEvents::RunTracks(){

  // cout << "Started Run Tracks" << endl;

  // cout << "Using primary tracks : " << std::boolalpha << doUsePrimTracks << endl; 

  const Int_t ntracks = mPicoDst->numberOfTracks();

  for(unsigned short itrk1 = 0; itrk1 < ntracks; itrk1++){
    // get track pointer
    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(itrk1));
    if(!trk1){ continue; }

    if (!IsAnAcceptableTrack(trk1, kTRUE)) continue;

    // FillPidHistograms(trk1);

    TVector3 mTrk1Mom;

    if(doUsePrimTracks) {
      mTrk1Mom = trk1->pMom();
    } else {
      mTrk1Mom = trk1->gMom(mVertex, Bfield);
    }

    for(unsigned short itrk2 = 0; itrk2 < ntracks; itrk2++){
      if (itrk2 >= ntracks) continue;
      if (itrk2 == itrk1) continue;
      // get track pointer
      StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(itrk2));
      if(!trk2){ continue; }

      if (!IsAnAcceptableTrack(trk2, kFALSE)) continue;

      int pid1 = IsWhatParticleNew(trk1);
      int pid2 = IsWhatParticleNew(trk2);

      if (abs(pid1) != 1 || abs(pid2) != 2) continue;

      double mass = InvariantMass(trk1, trk2);

      if ( mass < fInvMassSignal1 || mass > fInvMassSignal2 ) continue;

      if (fDebugLevel == 1) cout << "Tracks: " << itrk1 << "\t" << itrk2 << "\t" << mass << endl;

      TVector3 mTrk2Mom, mResMom;

      if(doUsePrimTracks) {
        mTrk2Mom = trk2->pMom();
      } else {
        mTrk2Mom = trk2->gMom(mVertex, Bfield);
      }

      mResMom = mTrk1Mom + mTrk2Mom;
      
      if (pid1*pid2==-2) {
        if ((mass >= fInvMassSignal1) && (mass <= fInvMassSignal2)) {
          fd0 = kTRUE;

          double toFillDaug[6] = {mCentMaker->GetCent9() + 0.5, mResMom.Perp()*trk1->charge(), mTrk1Mom.Perp(), mass, mTrk2Mom.Perp(), mResMom.PseudoRapidity() };

          fd0TrackIndices.push_back({itrk1, itrk2, mass});
          
          hD0CentPtEtaMDphiDaug->Fill(toFillDaug);

          if (fDebugLevel == 1){
            cout << fAnalysisMakerName << endl;
            cout << "Pion : " << mTrk1Mom.Perp() << "\t" << standardPhi(mTrk1Mom.Phi()) << "\t" << mTrk1Mom.PseudoRapidity() << endl;
            cout << "Kaon : " << mTrk2Mom.Perp() << "\t" << standardPhi(mTrk2Mom.Phi()) << "\t" << mTrk2Mom.PseudoRapidity() << endl;
            cout << "D0 : " << mResMom.Perp() << "\t" << standardPhi(mResMom.Phi()) << "\t" << mResMom.PseudoRapidity() << endl;
          } 
        }
      }
      if (pid1*pid2==2) {
        if ((mass >= fInvMassLSBg1) && (mass <= fInvMassLSBg2)) {
          fd0BgLS = kTRUE;

          double toFillDaug[6] = {mCentMaker->GetCent9() + 0.5, mResMom.Perp()*trk1->charge(), mTrk1Mom.Perp(), mass, mTrk2Mom.Perp(), mResMom.PseudoRapidity() };

          fd0BgLSTrackIndices.push_back({itrk1, itrk2, mass});
        
          hD0CentPtEtaMDphiDaugLikeSign->Fill(toFillDaug);

          if (fDebugLevel == 1){
            cout << fAnalysisMakerName << endl;
            cout << "Pion : " << mTrk1Mom.Perp() << "\t" << standardPhi(mTrk1Mom.Phi()) << "\t" << mTrk1Mom.PseudoRapidity() << endl;
            cout << "Kaon : " << mTrk2Mom.Perp() << "\t" << standardPhi(mTrk2Mom.Phi()) << "\t" << mTrk2Mom.PseudoRapidity() << endl;
            cout << "D0 Bg LS : " << mResMom.Perp() << "\t" << standardPhi(mResMom.Phi()) << "\t" << mResMom.PseudoRapidity() << endl;
          } 
        }
      }
    }
  }

  hNumberOfD0s->Fill(fd0TrackIndices.size());
  hNumberOfD0BgUS->Fill(fd0BgUSTrackIndices.size());
  hNumberOfD0BgLS->Fill(fd0BgLSTrackIndices.size());

  // cout << "Fraction of Mismatched Tracks: " << BadTracks << "\t" << AcceptedTracks << "\t" << BadTracks/AcceptedTracks << endl;
}

Bool_t StTagD0TPCEvents::IsAnAcceptableTrack(StPicoTrack *trk, bool dohistograms = kFALSE){ // This is updated as of Dec 25, 2021 to match the D0 Analysis from 2018. The cuts haven't been changed.

  TVector3 mTrkMom;
  mTrkMom = trk->gMom(mVertex, Bfield);

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  double px = mTrkMom.x();
  double py = mTrkMom.y();
  double pz = mTrkMom.z();
  double p = mTrkMom.Mag();
  short charge = trk->charge();

  if (charge == 0) return kFALSE;

  if (dohistograms){ hAcceptedPt[0]->Fill(pt); hAcceptedEta[0]->Fill(eta); }

  if (pt <= 0.2) return kFALSE; //This is different
  
  if (dohistograms){ hAcceptedPt[1]->Fill(pt); hAcceptedEta[1]->Fill(eta); }

  if (abs(eta) > 1) return kFALSE;

  if (dohistograms){ hAcceptedPt[2]->Fill(pt); hAcceptedEta[2]->Fill(eta); }

  if (trk->nHitsFit() < 20 || trk->nHitsFit() >= 50) return kFALSE; //This is different

  if (dohistograms){ hAcceptedPt[3]->Fill(pt); hAcceptedEta[3]->Fill(eta); }

  if (float(trk->nHitsFit())/float(trk->nHitsMax()) <= 0.52) return kFALSE;

  if (dohistograms){ hAcceptedPt[4]->Fill(pt); hAcceptedEta[4]->Fill(eta); }

  if (trk->nHitsDedx() <= 15) return kFALSE;

  if (dohistograms){ hAcceptedPt[5]->Fill(pt); hAcceptedEta[5]->Fill(eta); }

  if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 1.) return kFALSE;

  if (dohistograms){ hAcceptedPt[6]->Fill(pt); hAcceptedEta[6]->Fill(eta); }


  return kTRUE;
}

void StTagD0TPCEvents::FillPidHistograms(StPicoTrack *trk){

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, Bfield);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double p = mTrkMom.Mag();
  double eta = mTrkMom.PseudoRapidity();

  dEdXvp->Fill(p, trk->dEdx());

  bool toftrack = kFALSE;

  int tof_loc = trk->bTofPidTraitsIndex();

  if (tof_loc >= 0)
  {
    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    if (tofpointer){toftrack = kTRUE;}
  }

  if (toftrack){

    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    // double invbeta_from_tof = tofpointer->btofBeta();
    double invbeta_from_tof = GetTofBeta(trk);
    invbeta_from_tof = 1/invbeta_from_tof;

    invbetavp->Fill(p, invbeta_from_tof);

    double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
    double normalisedinvbeta_for_pi = abs(invbeta_from_tof-norm_invbeta_pi);
    normalised_invbetavpT_tof_pi->Fill(pt, normalisedinvbeta_for_pi);
    double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
    double normalisedinvbeta_for_ka = abs(invbeta_from_tof-norm_invbeta_ka);
    normalised_invbetavpT_tof_ka->Fill(pt, normalisedinvbeta_for_ka);
  }

  else{
      double zpi = trk->nSigmaPion();
      z_pi->Fill(pt, zpi);
      double zka = trk->nSigmaKaon();
      z_ka->Fill(pt, zka);
  }
}

bool StTagD0TPCEvents::IsTpcPion(StPicoTrack *trk){
  double zpi = trk->nSigmaPion();
  if (abs(zpi) < 2.) return true;
  return false;
}

bool StTagD0TPCEvents::IsTpcKaon(StPicoTrack *trk){
  double zka = trk->nSigmaKaon();
  if (abs(zka) < 2.) return true;
  return false;
}

bool StTagD0TPCEvents::IsPion(StPicoTrack *trk){
  if (!IsTpcPion(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_PION_PLUS*M_PION_PLUS / p / p + 1);
  double nsigma = (1./beta - oneOverBetaExpected)/0.013;

  double fpimin = -2.0;
  double fpimax;

  if (p >= 1.6) fpimax = 2.0 ;
  else fpimax = 5.43 - 2.14*p;

  if (nsigma > fpimin && nsigma < fpimax) return true;
  return false;
}

bool StTagD0TPCEvents::IsKaon(StPicoTrack *trk){
  if (!IsTpcKaon(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_KAON_PLUS*M_KAON_PLUS / p / p + 1);
  double nsigma = (1./beta - oneOverBetaExpected)/0.013;

  double fkaonmin;
  if (p >= 1.37) fkaonmin = -2.0 ;
  else fkaonmin = -7.54 + 5.83*p - 1.31*p*p;
  double fkaonmax;
  if (p >= 1.9) fkaonmax = 2.0 ;
  else fkaonmax = 8.69 - 6.02*p + 1.32*p*p;

  if (nsigma > fkaonmin && nsigma < fkaonmax) return true;
  return false;
}

void StTagD0TPCEvents::IsWhatParticleNew(StPicoTrack *trk, int &pid, double &m, double &e){ // This function is being rewritten to match the old analysis code closely. I do not expect differences, but it's nevertheless a good thing to be thorough.
  if(!IsAnAcceptableTrack(trk, kFALSE)){pid = 0; m = 0.; e = 0.; return;}

  pid = 0;
  m = 0.0;
  e = 0.0;

  double p = trk->gMom(mVertex, Bfield).Mag();
  short charge = trk->charge();

  if (IsPion(trk) && IsKaon(trk)) return; 

  if (IsPion(trk)) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
  if (IsKaon(trk)) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}

  return;
}

Int_t StTagD0TPCEvents::IsWhatParticleNew(StPicoTrack *trk){ // Just to get the PID out
  int pid;
  double m; 
  double e;
  IsWhatParticleNew(trk, pid, m, e);
  return pid;
}

void StTagD0TPCEvents::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum){

  int particle1, particle2;
  double mass1, mass2;
  double energy1, energy2;

  // cout << "IsWhatParticle called" << endl;

  IsWhatParticleNew(trk1, particle1, mass1, energy1);
  IsWhatParticleNew(trk2, particle2, mass2, energy2);

  TVector3 mTrk1Mom, mTrk2Mom;

  mTrk1Mom = trk1->gMom(mVertex, Bfield);

  mTrk2Mom = trk2->gMom(mVertex, Bfield);

  TLorentzVector mTrk1, mTrk2;

  mTrk1.SetXYZM(mTrk1Mom.X(), mTrk1Mom.Y(), mTrk1Mom.Z(), mass1);
  mTrk2.SetXYZM(mTrk2Mom.X(), mTrk2Mom.Y(), mTrk2Mom.Z(), mass2);

  TLorentzVector mD0;
  mD0 = mTrk1 + mTrk2;

  particle = particle1*particle2;

  invmass = mD0.M();

  momentum = mD0.Vect();

}

Double_t StTagD0TPCEvents::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2){
  int particle; 
  double invmass; 
  TVector3 momentum;

  InvariantMass(trk1, trk2, particle, invmass, momentum);
  return invmass;
}

float StTagD0TPCEvents::GetTofBeta(StPicoTrack *trk){
  int index2tof = trk->bTofPidTraitsIndex();
    
  float beta = std::numeric_limits<double>::quiet_NaN();

  if (index2tof < 0) return beta;

  StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(index2tof);
      
  if (tofPid)
  {
    float ylocal = tofPid->btofYLocal();
    if (abs(ylocal) >= 1.8) return beta;

    beta = tofPid->btofBeta();
    
    if (beta < 1e-4)
    {
        TVector3 const btofHitPosvec = tofPid->btofHitPos();

        StThreeVectorF const btofHitPos(btofHitPosvec.x(), btofHitPosvec.y(), btofHitPosvec.z()); 
        StThreeVectorF const pVertex(mVertex.x(), mVertex.y(), mVertex.z()); 
        
        StPicoPhysicalHelix helix = trk->helix(Bfield);
        float L = tofPathLength(&pVertex, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
    }

    // cout << tofPid->btofBeta() - beta << endl;
    // if (abs(tofPid->btofBeta() - beta) > 1e-4) cout << "Compare Betas : " << tofPid->btofBeta() << "\t" << beta << endl;
  }

  return beta;

}


// Centrality bin getter
Int_t StTagD0TPCEvents::GetFiveCentBin(Double_t scaledCent) const{
  int centbin = -99;
  // get centrality bin number
  if(scaledCent >= 0 && scaledCent <  10.0)  { centbin = 4; }
  else if(scaledCent >= 10.0 && scaledCent <  20.0)                { centbin = 3; }
  else if(scaledCent >= 20.0 && scaledCent <  40.0)                { centbin = 2; }
  else if(scaledCent >= 40.0 && scaledCent < 60.0)                { centbin = 1; }
  else if(scaledCent >= 60.0 && scaledCent <= 80.0)                { centbin = 0; }

  return centbin;
}

// D0 bin getter
Int_t StTagD0TPCEvents::GetD0PtBin(Double_t D0Pt) const{
  int d0bin = -99;
  // get centrality bin number
  if(D0Pt >= 0 && D0Pt <  0.5)                       { d0bin = 0; }
  else if(D0Pt >= 0.5 && D0Pt <  1.0)                { d0bin = 1; }
  else if(D0Pt >= 1.0 && D0Pt <  2.0)                { d0bin = 2; }
  else if(D0Pt >= 2.0 && D0Pt <  3.0)                { d0bin = 3; }
  else if(D0Pt >= 3.0 && D0Pt <  5.0)                { d0bin = 4; }
  else if(D0Pt >= 5.0 && D0Pt < 10.0)                { d0bin = 5; }

  return d0bin;
}