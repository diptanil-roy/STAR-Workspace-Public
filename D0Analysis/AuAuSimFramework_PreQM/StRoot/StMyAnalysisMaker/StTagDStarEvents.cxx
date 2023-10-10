// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StTagDStarEvents.h"
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
#include "TLorentzVector.h"
#include "TF1.h"


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

// #include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"

// // KFParticle includes
// #include "StRoot/StEvent/StDcaGeometry.h"
// #include "StRoot/StarRoot/KFParticle.h"

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

ClassImp(StTagDStarEvents)

//________________________________________________________________________
StTagDStarEvents::StTagDStarEvents(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTagDStarEvents::fRunFlagEnum
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
  fTrackDCAcut = 2.0;
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

  for(int i=0; i<10; i++) {numberofevents[i] = 0; numberoftracks[i] = 0;}
  //Mass cut
  fInvMassSignal1 = 1.80;
  fInvMassSignal2 = 1.92;

  fInvMassULBg1 = 1.70;
  fInvMassULBg2 = 2.08;

  fInvMassLSBg1 = 1.70;
  fInvMassLSBg2 = 2.08;

  fTagDStarEvents = kFALSE;

  fd0 = kFALSE;
  fd0BgUS = kFALSE;
  fd0BgLS = kFALSE;

  fD0SigTree = {0};
  fD0BgTree = {0};

  fDStarSigTree = {0};
  fDStarBgTree = {0};

  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();

  fdstarTrackIndices.clear();
  fdstarLSTrackIndices.clear();
}

//
//________________________________________________________________________
StTagDStarEvents::~StTagDStarEvents()
{ 

/*
  if (hEventZvertex_whole) delete hEventZvertex_whole;
  if (hEventZvertex_VPD) delete hEventZvertex_VPD;
  if (hEventZvertex_diff) delete hEventZvertex_diff;

  if (hCentrality) delete hCentrality;
  if (hMultiplicity) delete hMultiplicity;
  if (cuthistogram_event) delete cuthistogram_event;

  for (int i = 0; i < 10; i++){
    if (hAcceptedPt[i]) delete hAcceptedPt[i];
    if (hAcceptedEta[i]) delete hAcceptedEta[i];
  }

  if (z_pi) delete z_pi;
  if (z_ka) delete z_ka;
  if (normalised_invbetavpT_tof_pi) delete normalised_invbetavpT_tof_pi;
  if (normalised_invbetavpT_tof_ka) delete normalised_invbetavpT_tof_ka;
  if (dEdXvp) delete dEdXvp;
  if (invbetavp) delete invbetavp;

  if (decaylengthd0US) delete decaylengthd0US;
  if (distancepikUS) delete distancepikUS;
  if (distanced0PVUS) delete distanced0PVUS;
  if (dcakPVUS) delete dcakPVUS;
  if (dcapiPVUS) delete dcapiPVUS;

  if (decaylengthd0LS) delete decaylengthd0LS;
  if (distancepikLS) delete distancepikLS;
  if (distanced0PVLS) delete distanced0PVLS;
  if (dcakPVLS) delete dcakPVLS;
  if (dcapiPVLS) delete dcapiPVLS;

  if (kaonpt) delete kaonpt;
  if (pionpt) delete pionpt;
  if (d0pt) delete d0pt;
  if (kaonpionpt) delete kaonpionpt;
  if (kaoneta) delete kaoneta;
  if (pioneta) delete pioneta;
  if (d0eta) delete d0eta;
  if (kaonphi) delete kaonphi;
  if (pionphi) delete pionphi;
  if (d0phi) delete d0phi;

  if (kaonbgUSpt) delete kaonbgUSpt;
  if (pionbgUSpt) delete pionbgUSpt;
  if (d0bgUSpt) delete d0bgUSpt;
  if (kaonpionbgUSpt) delete kaonpionbgUSpt;
  if (kaonbgUSeta) delete kaonbgUSeta;
  if (pionbgUSeta) delete pionbgUSeta;
  if (d0bgUSeta) delete d0bgUSeta;
  if (kaonbgUSphi) delete kaonbgUSphi;
  if (pionbgUSphi) delete pionbgUSphi;
  if (d0bgUSphi) delete d0bgUSphi;

  if (kaonbgLSpt) delete kaonbgLSpt;
  if (pionbgLSpt) delete pionbgLSpt;
  if (d0bgLSpt) delete d0bgLSpt;
  if (kaonpionbgLSpt) delete kaonpionbgLSpt;
  if (kaonbgLSeta) delete kaonbgLSeta;
  if (pionbgLSeta) delete pionbgLSeta;
  if (d0bgLSeta) delete d0bgLSeta;
  if (kaonbgLSphi) delete kaonbgLSphi;
  if (pionbgLSphi) delete pionbgLSphi;
  if (d0bgLSphi) delete d0bgLSphi;

  if (invmass) delete invmass;
  if (invmassbg) delete invmassbg;

  if (invmassvpt) delete invmassvpt;
  if (invmassbgvpt) delete invmassbgvpt;

  if (invmassDStar) delete invmassDStar;
  if (invmassbgDStar) delete invmassbgDStar;

  if (invmassDStarvpt) delete invmassDStarvpt;
  if (invmassbgDStarvpt) delete invmassbgDStarvpt;

  if (invmassdiffDStarvpt) delete invmassdiffsideDStarvpt;
  if (invmassdiffsideDStarvpt) delete invmassdiffsideDStarvpt;
  if (invmassdiffwrongDStarvpt) delete invmassdiffwrongDStarvpt;

  if (invmassdiff) delete invmassdiff;
  if (invmassdiffside) delete invmassdiffside;
  if (invmassdiffwrong) delete invmassdiffwrong;

  if (hNumberOfD0s) delete hNumberOfD0s;
  if (hNumberOfD0BgUS) delete hNumberOfD0BgUS;
  if (hNumberOfD0BgLS) delete hNumberOfD0BgLS;

  */
  if(mEmcPosition) delete mEmcPosition;

}

//________________________________________________________________________
Int_t StTagDStarEvents::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  if (fTagDStarEvents){
    DeclareTree();
    BookTree();
  }
  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTagDStarEvents::Finish() { 

  cout << "StTagDStarEvents::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),  "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());

    const char *histdir = Form("%s/%s", GetName(), "Histograms");
    fout->mkdir(histdir);
    fout->cd(histdir);
    WriteHistograms();

    if (fTagDStarEvents){
      fout->cd(GetName());
      WriteTree();
    }

    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StTagDStarEvents::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

void StTagDStarEvents::DeclareTree() {

  if (fTagDStarEvents){

    TString d0sigtreename;

    d0sigtreename = "D0UnlikeSign";

    TString d0bgtreename;

    d0bgtreename = "D0LikeSign";

    d0sigtree = new TTree(d0sigtreename.Data(), d0sigtreename.Data());
    d0bgtree = new TTree(d0bgtreename.Data(), d0bgtreename.Data());

  
    TString dstarsigtreename;

    dstarsigtreename = "DStarUnlikeSign";

    TString dstarbgtreename;

    dstarbgtreename = "DStarLikeSign";

    dstarsigtree = new TTree(dstarsigtreename.Data(), dstarsigtreename.Data());
    dstarbgtree = new TTree(dstarbgtreename.Data(), dstarbgtreename.Data());
  }

}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StTagDStarEvents::DeclareHistograms() {

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

  cuthistogram_event = new TH1F("cuthistogram_event", "cuthistogram_event", 10, 0, 10);

  // Pt After Track Cuts

  for (int i = 0; i < 10; i++){
    hAcceptedPt[i] = new TH1F(Form("hAcceptedPt_%i", i), Form("hAcceptedPt_%i", i), 1000, 0, 30);
    hAcceptedEta[i] = new TH1F(Form("hAcceptedEta%i", i), Form("hAcceptedEta%i", i), 1000, -2, 2);
  }
  

  z_pi = new TH2F("z_pi", "z_pi", 3000, 0, 30, 2000, -10, 10);
  z_ka = new TH2F("z_ka", "z_ka", 3000, 0, 30, 2000, -10, 10);
  normalised_invbetavpT_tof_pi = new TH2F("normalised_invbetavpT_tof_pi", "normalised_invbetavpT_tof_pi", 3000, 0, 30, 2000, -10, 10);
  normalised_invbetavpT_tof_ka = new TH2F("normalised_invbetavpT_tof_ka", "normalised_invbetavpT_tof_ka", 3000, 0, 30, 2000, -10, 10);

  dEdXvp = new TH2F("dEdXvp", "dEdXvp", 3000, 0, 30, 4000, -20, 20);
  invbetavp = new TH2F("invbetavp", "invbetavp", 3000, 0, 30, 4000, -20, 20);


  // Topological histograms

  decaylengthd0US = new TH1F("decaylengthd0US", "decaylengthd0US", 10000, 0, 1);
  distancepikUS = new TH1F("distancepikUS", "distancepikUS", 10000, 0, 1);
  distanced0PVUS = new TH1F("distanced0PVUS", "distanced0PVUS", 10000, 0, 1);
  dcakPVUS = new TH1F("dcakPVUS", "dcakPVUS", 10000, 0, 1);
  dcapiPVUS = new TH1F("dcapiPVUS", "dcapiPVUS", 10000, 0, 1);

  decaylengthd0LS = new TH1F("decaylengthd0LS", "decaylengthd0LS", 10000, 0, 1);
  distancepikLS = new TH1F("distancepikLS", "distancepikLS", 10000, 0, 1);
  distanced0PVLS = new TH1F("distanced0PVLS", "distanced0PVLS", 10000, 0, 1);
  dcakPVLS = new TH1F("dcakPVLS", "dcakPVLS", 10000, 0, 1);
  dcapiPVLS = new TH1F("dcapiPVLS", "dcapiPVLS", 10000, 0, 1);

  // hD0US_pT = new TH1F("hD0US_pT", "hD0US_pT", 1000, 0, 50);
  // hD0LS_pT = new TH1F("hD0LS_pT", "hD0LS_pT", 1000, 0, 50);


  //Invariant Mass Plots

  kaonpt = new TH1F("kaonpt", "kaonpt", 1000, 0, 30);
  pionpt = new TH1F("pionpt", "pionpt", 1000, 0, 30);
  d0pt = new TH1F("d0pt", "d0pt", 1000, 0, 30);
  kaonpionpt = new TH2F("kaonpionpt", "kaonpionpt", 1000, 0, 30, 1000, 0, 30);
  kaoneta = new TH1F("kaoneta", "kaoneta", 1000, -1, 1);
  pioneta = new TH1F("pioneta", "pioneta", 1000, -1, 1);
  d0eta = new TH1F("d0eta", "d0eta", 1000, -5, 5);
  kaonphi = new TH1F("kaonphi", "kaonphi", 1000, -10, 10);
  pionphi = new TH1F("pionphi", "pionphi", 1000, -10, 10);
  d0phi = new TH1F("d0phi", "d0phi", 1000, -10, 10);


  kaonbgUSpt = new TH1F("kaonbgUSpt", "kaonbgUSpt", 1000, 0, 30);
  pionbgUSpt = new TH1F("pionbgUSpt", "pionbgUSpt", 1000, 0, 30);
  d0bgUSpt = new TH1F("d0bgUSpt", "d0bgUSpt", 1000, 0, 30);
  kaonpionbgUSpt = new TH2F("kaonpionbgUSpt", "kaonpionbgUSpt", 1000, 0, 30, 1000, 0, 30);
  kaonbgUSeta = new TH1F("kaonbgUSeta", "kaonbgUSeta", 1000, -1, 1);
  pionbgUSeta = new TH1F("pionbgUSeta", "pionbgUSeta", 1000, -1, 1);
  d0bgUSeta = new TH1F("d0bgUSeta", "d0bgUSeta", 1000, -5, 5);
  kaonbgUSphi = new TH1F("kaonbgUSphi", "kaonbgUSphi", 1000, -10, 10);
  pionbgUSphi = new TH1F("pionbgUSphi", "pionbgUSphi", 1000, -10, 10);
  d0bgUSphi = new TH1F("d0bgUSphi", "d0bgUSphi", 1000, -10, 10);


  kaonbgLSpt = new TH1F("kaonbgLSpt", "kaonbgLSpt", 1000, 0, 30);
  pionbgLSpt = new TH1F("pionbgLSpt", "pionbgLSpt", 1000, 0, 30);
  d0bgLSpt = new TH1F("d0bgLSpt", "d0bgLSpt", 1000, 0, 30);
  kaonpionbgLSpt = new TH2F("kaonpionbgLSpt", "kaonpionbgLSpt", 1000, 0, 30, 1000, 0, 30);
  kaonbgLSeta = new TH1F("kaonbgLSeta", "kaonbgLSeta", 1000, -1, 1);
  pionbgLSeta = new TH1F("pionbgLSeta", "pionbgLSeta", 1000, -1, 1);
  d0bgLSeta = new TH1F("d0bgLSeta", "d0bgLSeta", 1000, -5, 5);
  kaonbgLSphi = new TH1F("kaonbgLSphi", "kaonbgLSphi", 1000, -10, 10);
  pionbgLSphi = new TH1F("pionbgLSphi", "pionbgLSphi", 1000, -10, 10);
  d0bgLSphi = new TH1F("d0bgLSphi", "d0bgLSphi", 1000, -10, 10);

  invmass = new TH1F("invmass", "invmass", 1000, 0.5, 2.5);
  invmassbg = new TH1F("invmassbg", "invmassbg", 1000, 0.5, 2.5);

  invmassvpt = new TH2F("invmassvpt", "invmassvpt", 1000, 0.5, 2.5, 50, 0., 50);
  invmassbgvpt = new TH2F("invmassbgvpt", "invmassbgvpt", 1000, 0.5, 2.5, 50, 0., 50.);

  invmassDStar = new TH1F ("invmassDStar", "invmassDStar", 1000, 0.5, 2.5);
  invmassbgDStar = new TH1F("invmassbgDStar", "invmassbgDStar", 1000, 0.5, 2.5);

  invmassDStarvpt = new TH2F("invmassDStarvpt", "invmassDStarvpt", 1000, 0.5, 2.5, 50, 0., 50);
  invmassbgDStarvpt = new TH2F("invmassbgDStarvpt", "invmassbgDStarvpt", 1000, 0.5, 2.5, 50, 0., 50.);

  invmassdiff = new TH1F("invmassdiff", "invmassdiff", 1000, 0.1, 0.2);
  invmassdiffside = new TH1F("invmassdiffside", "invmassdiffside", 1000, 0.1, 0.2);
  invmassdiffwrong = new TH1F("invmassdiffwrong", "invmassdiffwrong", 1000, 0.1, 0.2);

  invmassdiffDStarvpt = new TH2F("invmassdiffDStarvpt", "invmassdiffDStarvpt", 1000, 0.1, 0.2, 50, 0., 50);
  invmassdiffsideDStarvpt = new TH2F("invmassdiffsideDStarvpt", "invmassdiffsideDStarvpt", 1000, 0.1, 0.2, 50, 0., 50.);
  invmassdiffwrongDStarvpt = new TH2F("invmassdiffwrongDStarvpt", "invmassdiffwrongDStarvpt", 1000, 0.1, 0.2, 50, 0., 50.);

  hNumberOfD0s = new TH1F("hNumberOfD0s", "hNumberOfD0s", 21, -0.5, 20.5);
  hNumberOfD0BgUS = new TH1F("hNumberOfD0BgUS", "hNumberOfD0BgUS", 21, -0.5, 20.5);
  hNumberOfD0BgLS = new TH1F("hNumberOfD0BgLS", "hNumberOfD0BgLS", 21, -0.5, 20.5);

}
//
// write histograms
//_____________________________________________________________________________
void StTagDStarEvents::WriteHistograms() {

  // Centrality Histograms

  hEventZvertex_whole->Write();
  hEventZvertex_VPD->Write();
  hEventZvertex_diff->Write();

  hCentrality->Write();
  hMultiplicity->Write();

  // Event Cuts Histograms
  const char *event_cuts[10] = {"Total Events", "MinBias","|V_{z}| < 6 cm", "|V_{z} - V_{z(VPD)}| < 3 cm", "nBEMCMatch > 0", "nBTOFMatch > 0", "HT1", "HT2", "HT3", ""};
  const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 20", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

  for (int i=1; i <= 10; i++){
    cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
    cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
    // cuthistogram_track->SetBinContent(i, numberoftracks[i-1]);
    // cuthistogram_track->GetXaxis()->SetBinLabel(i, track_cuts[i-1]);
  }
  cuthistogram_event->Write();
  // cuthistogram_track->Write();

  // Track Cut Histograms
  for (int i = 0; i < 8; i++){
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

  // Topological Histograms

  decaylengthd0US->Write();
  distancepikUS->Write();
  distanced0PVUS->Write();
  dcakPVUS->Write();
  dcapiPVUS->Write();

  decaylengthd0LS->Write();
  distancepikLS->Write();
  distanced0PVLS->Write();
  dcakPVLS->Write();
  dcapiPVLS->Write();

  // Pt and Invmass Histograms for KPi pairs after topological cuts

  kaonpt->Write();
  pionpt->Write();
  d0pt->Write();
  kaonpionpt->Write();
  kaoneta->Write();
  pioneta->Write();
  d0eta->Write();
  kaonphi->Write();
  pionphi->Write();
  d0phi->Write();

  kaonbgUSpt->Write();
  pionbgUSpt->Write();
  d0bgUSpt->Write();
  kaonpionbgUSpt->Write();
  kaonbgUSeta->Write();
  pionbgUSeta->Write();
  d0bgUSeta->Write();
  kaonbgUSphi->Write();
  pionbgUSphi->Write();
  d0bgUSphi->Write();

  kaonbgLSpt->Write();
  pionbgLSpt->Write();
  d0bgLSpt->Write();
  kaonpionbgLSpt->Write();
  kaonbgLSeta->Write();
  pionbgLSeta->Write();
  d0bgLSeta->Write();
  kaonbgLSphi->Write();
  pionbgLSphi->Write();
  d0bgLSphi->Write();

  invmass->Write();
  invmassbg->Write();

  invmassvpt->Write();
  invmassbgvpt->Write();

  invmassDStar->Write();
  invmassbgDStar->Write();

  invmassDStarvpt->Write();
  invmassbgDStarvpt->Write();

  invmassdiff->Write();
  invmassdiffside->Write();
  invmassdiffwrong->Write();

  invmassdiffDStarvpt->Write();
  invmassdiffsideDStarvpt->Write();
  invmassdiffwrongDStarvpt->Write();

  hNumberOfD0s->Write();
  hNumberOfD0BgUS->Write();
  hNumberOfD0BgLS->Write();
}

void StTagDStarEvents::WriteTree(){
  
  if (fTagDStarEvents){
    d0sigtree->Write();
    d0bgtree->Write();

    dstarsigtree->Write();
    dstarbgtree->Write();
  }
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTagDStarEvents::Clear(Option_t *opt) {
  // fJets->Clear();
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();

  fdstarTrackIndices.clear();
  fdstarLSTrackIndices.clear();
}

void StTagDStarEvents::BookTree()
{
  if (fTagDStarEvents){
    // Branches to save D0 Sig event info
    d0sigtree->Branch("RunID", &fD0SigTree.runid, "runid/I");
    d0sigtree->Branch("EventId", &fD0SigTree.eventid, "eventid/I");
    d0sigtree->Branch("RefMult", &fD0SigTree.refmult, "refmult/F");
    d0sigtree->Branch("Centrality", &fD0SigTree.centrality, "centrality/F");
    d0sigtree->Branch("Triggers", &fD0SigTree.triggers);
    d0sigtree->Branch("PrimaryVertex", &fD0SigTree.primaryvertex);
    d0sigtree->Branch("PrimaryVertexErr", &fD0SigTree.primaryvertexerror);

    d0sigtree->Branch("D0Mass", &fD0SigTree.d0mass, "d0mass/F");
    d0sigtree->Branch("D0Pt", &fD0SigTree.d0pt, "d0pt/F");
    d0sigtree->Branch("D0Eta", &fD0SigTree.d0eta, "d0eta/F");
    d0sigtree->Branch("D0Phi", &fD0SigTree.d0phi, "d0phi/F");

    d0sigtree->Branch("PionPt", &fD0SigTree.pionpt, "pionpt/F");
    d0sigtree->Branch("PionEta", &fD0SigTree.pioneta, "pioneta/F");
    d0sigtree->Branch("PionPhi", &fD0SigTree.pionphi, "pionphi/F");
    d0sigtree->Branch("PionCharge", &fD0SigTree.pioncharge, "pioncharge/F");
    d0sigtree->Branch("PionDCA", &fD0SigTree.piondca, "piondca/F");
    d0sigtree->Branch("PionNSigma", &fD0SigTree.nsigmapion, "nsigmapion/F");

    d0sigtree->Branch("KaonPt", &fD0SigTree.kaonpt, "kaonpt/F");
    d0sigtree->Branch("KaonEta", &fD0SigTree.kaoneta, "kaoneta/F");
    d0sigtree->Branch("KaonPhi", &fD0SigTree.kaonphi, "kaonphi/F");
    d0sigtree->Branch("KaonCharge", &fD0SigTree.kaoncharge, "kaoncharge/F");
    d0sigtree->Branch("KaonDCA", &fD0SigTree.kaondca, "kaondca/F");
    d0sigtree->Branch("KaonNSigma", &fD0SigTree.nsigmakaon, "nsigmakaon/F");

    // Branches to save BG event info
    d0bgtree->Branch("RunID", &fD0BgTree.runid, "runid/I");
    d0bgtree->Branch("EventId", &fD0BgTree.eventid, "eventid/I");
    d0bgtree->Branch("RefMult", &fD0BgTree.refmult, "refmult/F");
    d0bgtree->Branch("Centrality", &fD0BgTree.centrality, "centrality/F");
    d0bgtree->Branch("Triggers", &fD0BgTree.triggers);
    d0bgtree->Branch("PrimaryVertex", &fD0BgTree.primaryvertex);
    d0bgtree->Branch("PrimaryVertexErr", &fD0BgTree.primaryvertexerror);

    d0bgtree->Branch("D0Mass", &fD0BgTree.d0mass, "d0mass/F");
    d0bgtree->Branch("D0Pt", &fD0BgTree.d0pt, "d0pt/F");
    d0bgtree->Branch("D0Eta", &fD0BgTree.d0eta, "d0eta/F");
    d0bgtree->Branch("D0Phi", &fD0BgTree.d0phi, "d0phi/F");

    d0bgtree->Branch("PionPt", &fD0BgTree.pionpt, "pionpt/F");
    d0bgtree->Branch("PionEta", &fD0BgTree.pioneta, "pioneta/F");
    d0bgtree->Branch("PionPhi", &fD0BgTree.pionphi, "pionphi/F");
    d0bgtree->Branch("PionCharge", &fD0BgTree.pioncharge, "pioncharge/F");
    d0bgtree->Branch("PionDCA", &fD0BgTree.piondca, "piondca/F");
    d0bgtree->Branch("PionNSigma", &fD0BgTree.nsigmapion, "nsigmapion/F");

    d0bgtree->Branch("KaonPt", &fD0BgTree.kaonpt, "kaonpt/F");
    d0bgtree->Branch("KaonEta", &fD0BgTree.kaoneta, "kaoneta/F");
    d0bgtree->Branch("KaonPhi", &fD0BgTree.kaonphi, "kaonphi/F");
    d0bgtree->Branch("KaonCharge", &fD0BgTree.kaoncharge, "kaoncharge/F");
    d0bgtree->Branch("KaonDCA", &fD0BgTree.kaondca, "kaondca/F");
    d0bgtree->Branch("KaonNSigma", &fD0BgTree.nsigmakaon, "nsigmakaon/F");

    // Branches to save DStar Sig event info
    dstarsigtree->Branch("RunID", &fDStarSigTree.runid, "runid/I");
    dstarsigtree->Branch("EventId", &fDStarSigTree.eventid, "eventid/I");
    dstarsigtree->Branch("RefMult", &fDStarSigTree.refmult, "refmult/F");
    dstarsigtree->Branch("Centrality", &fDStarSigTree.centrality, "centrality/F");
    dstarsigtree->Branch("Triggers", &fDStarSigTree.triggers);
    dstarsigtree->Branch("PrimaryVertex", &fDStarSigTree.primaryvertex);
    dstarsigtree->Branch("PrimaryVertexErr", &fDStarSigTree.primaryvertexerror);

    dstarsigtree->Branch("D0Mass", &fDStarSigTree.d0mass, "d0mass/F");
    dstarsigtree->Branch("D0Pt", &fDStarSigTree.d0pt, "d0pt/F");
    dstarsigtree->Branch("D0Eta", &fDStarSigTree.d0eta, "d0eta/F");
    dstarsigtree->Branch("D0Phi", &fDStarSigTree.d0phi, "d0phi/F");

    dstarsigtree->Branch("DStarMass", &fDStarSigTree.dstarmass, "dstarmass/F");
    dstarsigtree->Branch("DStarPt", &fDStarSigTree.dstarpt, "dstarpt/F");
    dstarsigtree->Branch("DStarEta", &fDStarSigTree.dstareta, "dstareta/F");
    dstarsigtree->Branch("DStarPhi", &fDStarSigTree.dstarphi, "dstarphi/F");

    dstarsigtree->Branch("PionPt", &fDStarSigTree.pionpt, "pionpt/F");
    dstarsigtree->Branch("PionEta", &fDStarSigTree.pioneta, "pioneta/F");
    dstarsigtree->Branch("PionPhi", &fDStarSigTree.pionphi, "pionphi/F");
    dstarsigtree->Branch("PionCharge", &fDStarSigTree.pioncharge, "pioncharge/F");
    dstarsigtree->Branch("PionDCA", &fDStarSigTree.piondca, "piondca/F");
    dstarsigtree->Branch("PionNSigma", &fDStarSigTree.nsigmapion, "nsigmapion/F");

    dstarsigtree->Branch("KaonPt", &fDStarSigTree.kaonpt, "kaonpt/F");
    dstarsigtree->Branch("KaonEta", &fDStarSigTree.kaoneta, "kaoneta/F");
    dstarsigtree->Branch("KaonPhi", &fDStarSigTree.kaonphi, "kaonphi/F");
    dstarsigtree->Branch("KaonCharge", &fDStarSigTree.kaoncharge, "kaoncharge/F");
    dstarsigtree->Branch("KaonDCA", &fDStarSigTree.kaondca, "kaondca/F");
    dstarsigtree->Branch("KaonNSigma", &fDStarSigTree.nsigmakaon, "nsigmakaon/F");

    dstarsigtree->Branch("SPionPt", &fDStarSigTree.spionpt, "spionpt/F");
    dstarsigtree->Branch("SPionEta", &fDStarSigTree.spioneta, "spioneta/F");
    dstarsigtree->Branch("SPionPhi", &fDStarSigTree.spionphi, "spionphi/F");
    dstarsigtree->Branch("SPionCharge", &fDStarSigTree.spioncharge, "spioncharge/F");
    dstarsigtree->Branch("SPionDCA", &fDStarSigTree.spiondca, "spiondca/F");
    dstarsigtree->Branch("SPionNSigma", &fDStarSigTree.snsigmapion, "snsigmapion/F");

    // Branches to save DStar Bg event info
    dstarbgtree->Branch("RunID", &fDStarBgTree.runid, "runid/I");
    dstarbgtree->Branch("EventId", &fDStarBgTree.eventid, "eventid/I");
    dstarbgtree->Branch("RefMult", &fDStarBgTree.refmult, "refmult/F");
    dstarbgtree->Branch("Centrality", &fDStarBgTree.centrality, "centrality/F");
    dstarbgtree->Branch("Triggers", &fDStarBgTree.triggers);
    dstarbgtree->Branch("PrimaryVertex", &fDStarBgTree.primaryvertex);
    dstarbgtree->Branch("PrimaryVertexErr", &fDStarBgTree.primaryvertexerror);

    dstarbgtree->Branch("D0Mass", &fDStarBgTree.d0mass, "d0mass/F");
    dstarbgtree->Branch("D0Pt", &fDStarBgTree.d0pt, "d0pt/F");
    dstarbgtree->Branch("D0Eta", &fDStarBgTree.d0eta, "d0eta/F");
    dstarbgtree->Branch("D0Phi", &fDStarBgTree.d0phi, "d0phi/F");

    dstarbgtree->Branch("DStarMass", &fDStarBgTree.dstarmass, "dstarmass/F");
    dstarbgtree->Branch("DStarPt", &fDStarBgTree.dstarpt, "dstarpt/F");
    dstarbgtree->Branch("DStarEta", &fDStarBgTree.dstareta, "dstareta/F");
    dstarbgtree->Branch("DStarPhi", &fDStarBgTree.dstarphi, "dstarphi/F");

    dstarbgtree->Branch("PionPt", &fDStarBgTree.pionpt, "pionpt/F");
    dstarbgtree->Branch("PionEta", &fDStarBgTree.pioneta, "pioneta/F");
    dstarbgtree->Branch("PionPhi", &fDStarBgTree.pionphi, "pionphi/F");
    dstarbgtree->Branch("PionCharge", &fDStarBgTree.pioncharge, "pioncharge/F");
    dstarbgtree->Branch("PionDCA", &fDStarBgTree.piondca, "piondca/F");
    dstarbgtree->Branch("PionNSigma", &fDStarBgTree.nsigmapion, "nsigmapion/F");

    dstarbgtree->Branch("KaonPt", &fDStarBgTree.kaonpt, "kaonpt/F");
    dstarbgtree->Branch("KaonEta", &fDStarBgTree.kaoneta, "kaoneta/F");
    dstarbgtree->Branch("KaonPhi", &fDStarBgTree.kaonphi, "kaonphi/F");
    dstarbgtree->Branch("KaonCharge", &fDStarBgTree.kaoncharge, "kaoncharge/F");
    dstarbgtree->Branch("KaonDCA", &fDStarBgTree.kaondca, "kaondca/F");
    dstarbgtree->Branch("KaonNSigma", &fDStarBgTree.nsigmakaon, "nsigmakaon/F");

    dstarbgtree->Branch("SPionPt", &fDStarBgTree.spionpt, "spionpt/F");
    dstarbgtree->Branch("SPionEta", &fDStarBgTree.spioneta, "spioneta/F");
    dstarbgtree->Branch("SPionPhi", &fDStarBgTree.spionphi, "spionphi/F");
    dstarbgtree->Branch("SPionCharge", &fDStarBgTree.spioncharge, "spioncharge/F");
    dstarbgtree->Branch("SPionDCA", &fDStarBgTree.spiondca, "spiondca/F");
    dstarbgtree->Branch("SPionNSigma", &fDStarBgTree.snsigmapion, "snsigmapion/F");

  }
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTagDStarEvents::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  fd0 = kFALSE;
  fd0BgUS = kFALSE;
  fd0BgLS = kFALSE;

  fdstar = kFALSE;
  fdstarLS = kFALSE;

  //zero out global vectors
  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();

  fdstarTrackIndices.clear();
  fdstarLSTrackIndices.clear();

  fD0SigTree = {0};
  fD0BgTree = {0};
  
  fDStarSigTree = {0};
  fDStarBgTree = {0};

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

  // cout << "B = " << Bfield << endl;
  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  double zVtx_VPD = mPicoEvent->vzVpd();

  hEventZvertex_whole->Fill(zVtx);
  hEventZvertex_VPD->Fill(zVtx_VPD);
  hEventZvertex_diff->Fill(zVtx - zVtx_VPD);

  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


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
    int centbin = mCentMaker->GetRef16();
    double refCorr2 = mCentMaker->GetRefCorr2();
    frefCorr2 = refCorr2;
    fCentralityScaled = mCentMaker->GetCentScaled();

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
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"EventTriggers: ";
    for(unsigned int i=0; i<mytriggers.size(); i++) {
      if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
    }
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<endl;

    numberofevents[0]++;

    // check for MB/HT event
    bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
    bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
    bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
    bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);

    bool fRunForMB = kFALSE;
    if(doppAnalysis) fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
    if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

    int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
    int arrHT1_Run14[] = {450201, 450211, 460201};
    int arrHT2_Run14[] = {450202, 450212, 460202, 460212};
    int arrHT3_Run14[] = {450203, 450213, 460203};

    int arrMB_Run12[] = {370001, 370011};
    int arrHT1_Run12[] = {370511, 370546};
    int arrHT2_Run12[] = {370521, 370522, 370531, 370980};

    if (!doppAnalysis){

      bool matchMB = kFALSE;

      for(int i = 0; i < sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14); i++) {
        if(mPicoEvent->isTrigger(arrMB5_Run14[i])) matchMB = kTRUE;
        if(matchMB) break;
      }

      if (!matchMB) return kStOk;

    }

    else{
      bool matchMB = kFALSE;

      for(int i = 0; i < sizeof(arrMB_Run12)/sizeof(*arrMB_Run12); i++) {
        if(mPicoEvent->isTrigger(arrMB_Run12[i])) matchMB = kTRUE;
        if(matchMB) break;
      }

      if (!matchMB) return kStOk;
    }

    numberofevents[1]++;
    // ======================== end of Triggers ============================= //

    

    if (abs(zVtx) > 60) return kStOk; //Set for PP
    numberofevents[2]++;

    if (abs(zVtx - zVtx_VPD) > 6) return kStOk; // Set For PP
    numberofevents[3]++;
     
    if (mPicoEvent->nBEMCMatch() == 0) return kStOk;
    numberofevents[4]++;

    if (mPicoEvent->nBTOFMatch() == 0) return kStOk;
    numberofevents[5]++;

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

      if (matchHT1 || matchHT2 || matchHT3) numberofevents[6]++;


      if (matchHT2 || matchHT3) numberofevents[7]++;

     
      if (matchHT3) numberofevents[8]++;
    // cout << "MinBias Event Accepted." << endl;
    }
  }

  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // cout << "Number of tracks found = " << ntracks << endl;

  RunTracks();
  // ProcessJetForJetShape();
  if (fDebugLevel == 1) TestTracks();
  return kStOK;
}


//__________________________________________________________________________________________
  
void StTagDStarEvents::RunTracks(){

  // cout << "Started Run Tracks" << endl; 

  fD0SigTree = {0};
  fD0BgTree = {0};
  fDStarSigTree = {0};
  fDStarBgTree = {0};

  fd0TrackIndices.clear();
  fd0BgUSTrackIndices.clear();
  fd0BgLSTrackIndices.clear();

  fdstarTrackIndices.clear();
  fdstarLSTrackIndices.clear();

  const Int_t ntracks = mPicoDst->numberOfTracks();

  for(unsigned short itrk1 = 0; itrk1 < ntracks; itrk1++){
    // get track pointer
    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(itrk1));
    if(!trk1){ continue; }

    if (!IsAnAcceptableTrack(trk1, kTRUE)) continue;

    FillPidHistograms(trk1);

    TVector3 mTrk1Mom;

    if(doUsePrimTracks) {
      mTrk1Mom = trk1->pMom();
    } else {
      mTrk1Mom = trk1->gMom(mVertex, Bfield);
    }

    for(unsigned short itrk2 = itrk1 + 1; itrk2 < ntracks; itrk2++){
      if (itrk2 >= ntracks) continue;
      // get track pointer
      StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(itrk2));
      if(!trk2){ continue; }

      if (!IsAnAcceptableTrack(trk2, kFALSE)) continue;

      if (fTopoLevel > 0 && TopologicalCuts(trk1, trk2, kTRUE) != fTopoLevel) continue; // If I set fTopoLevel to 0, it doesn't invoke Topological cuts at all! Sometimes my genius, ..., it generates gravity!

      TVector3 mTrk2Mom, mResMom;

      if(doUsePrimTracks) {
        mTrk2Mom = trk2->pMom();
      } else {
        mTrk2Mom = trk2->gMom(mVertex, Bfield);
      }

      mResMom = mTrk1Mom + mTrk2Mom;


      int pid1 = IsWhatParticle(trk1);
      int pid2 = IsWhatParticle(trk2);

      double mass = InvariantMass(trk1, trk2);

      int trackid1 = trk1->id();
      int trackid2 = trk2->id();

      if (fDebugLevel > 1){cout << "Track IDs : " << itrk1 << "\t" << itrk2 << "\t" << trackid1 << "\t" << trackid2 << endl;}

      // if (pid1*pid2!=-2) continue;

      if ((mass < fInvMassULBg1) || (mass > fInvMassULBg2)) continue; // Only tracks between 1.7, 2.02

      unsigned short piontrack = (abs(pid1) == 1) ? itrk1 : itrk2; 
      unsigned short kaontrack = (abs(pid1) == 2) ? itrk1 : itrk2;

      TVector3 pionMom, kaonMom;
      int charge1, charge2;

      if (piontrack == itrk1) {pionMom = mTrk1Mom; kaonMom = mTrk2Mom; charge1 = trk1->charge(); charge2 = trk2->charge();}
      else if (piontrack == itrk2) {pionMom = mTrk2Mom; kaonMom = mTrk1Mom; charge1 = trk2->charge(); charge2 = trk1->charge();}

      if (abs(pid1*pid2)!=2) continue;

      if (pid1*pid2==-2){

        if (mass >= fInvMassSignal1 && mass <= fInvMassSignal2) {fd0 = kTRUE; fd0TrackIndices.push_back({piontrack, kaontrack, mass});}
        else {fd0BgUS = kTRUE; fd0BgUSTrackIndices.push_back({piontrack, kaontrack, mass});}
  
        invmass->Fill(mass);
        invmassvpt->Fill(mass, mResMom.Perp());

        // cout << "Index Array Sizes = " << fd0TrackIndices.size() << "\t" << fd0BgUSTrackIndices.size() << endl;

        if (fTagDStarEvents){

          // cout << "This should not be printed." << endl;
          fD0SigTree = {0};

          fD0SigTree.runid = mPicoEvent->runId();
          fD0SigTree.eventid = mPicoEvent->eventId();
          fD0SigTree.refmult = frefCorr2;
          fD0SigTree.centrality = fCentralityScaled;
          fD0SigTree.triggers = mPicoEvent->triggerIds();

          vector<double> pv;
          pv.clear();
          pv.push_back(mPicoEvent->primaryVertex().X());
          pv.push_back(mPicoEvent->primaryVertex().Y());
          pv.push_back(mPicoEvent->primaryVertex().Z());

          fD0SigTree.primaryvertex = pv;

          vector<double> pverr;
          pverr.clear();
          pverr.push_back(mPicoEvent->primaryVertexError().X());
          pverr.push_back(mPicoEvent->primaryVertexError().Y());
          pverr.push_back(mPicoEvent->primaryVertexError().Z());

          fD0SigTree.primaryvertexerror = pverr;

          fD0SigTree.d0mass = mass;
          fD0SigTree.d0pt = mResMom.Perp();
          fD0SigTree.d0eta = mResMom.PseudoRapidity();
          fD0SigTree.d0phi = standardPhi(mResMom.Phi());

          fD0SigTree.pionpt = pionMom.Perp();
          fD0SigTree.pioneta = pionMom.PseudoRapidity();
          fD0SigTree.pionphi = standardPhi(pionMom.Phi());
          fD0SigTree.pioncharge = charge1;
          fD0SigTree.piondca = (piontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
          fD0SigTree.nsigmapion = (piontrack == itrk1) ? trk1->nSigmaPion() : trk2->nSigmaPion();

          fD0SigTree.kaonpt = kaonMom.Perp();
          fD0SigTree.kaoneta = kaonMom.PseudoRapidity();
          fD0SigTree.kaonphi = standardPhi(kaonMom.Phi());
          fD0SigTree.kaoncharge = charge2;
          fD0SigTree.kaondca = (kaontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
          fD0SigTree.nsigmakaon = (kaontrack == itrk1) ? trk1->nSigmaKaon() : trk2->nSigmaKaon();

          d0sigtree->Fill();
        }
      }

      if (pid1*pid2==2){
        fd0BgLS = kTRUE;
        fd0BgLSTrackIndices.push_back({piontrack, kaontrack, mass});

        // cout << "Index Array Sizes = " << fd0BgLSTrackIndices.size() << "\n";

        invmassbg->Fill(mass);
        invmassbgvpt->Fill(mass, mResMom.Perp());

        if (fTagDStarEvents){

          fD0BgTree = {0};

          fD0BgTree.runid = mPicoEvent->runId();
          fD0BgTree.eventid = mPicoEvent->eventId();
          fD0BgTree.refmult = frefCorr2;
          fD0BgTree.centrality = fCentralityScaled;
          fD0BgTree.triggers = mPicoEvent->triggerIds();

          vector<double> pv;
          pv.clear();
          pv.push_back(mPicoEvent->primaryVertex().X());
          pv.push_back(mPicoEvent->primaryVertex().Y());
          pv.push_back(mPicoEvent->primaryVertex().Z());

          fD0BgTree.primaryvertex = pv;

          vector<double> pverr;
          pverr.clear();
          pverr.push_back(mPicoEvent->primaryVertexError().X());
          pverr.push_back(mPicoEvent->primaryVertexError().Y());
          pverr.push_back(mPicoEvent->primaryVertexError().Z());

          fD0BgTree.primaryvertexerror = pverr;

          fD0BgTree.d0mass = mass;
          fD0BgTree.d0pt = mResMom.Perp();
          fD0BgTree.d0eta = mResMom.PseudoRapidity();
          fD0BgTree.d0phi = standardPhi(mResMom.Phi());

          fD0BgTree.pionpt = pionMom.Perp();
          fD0BgTree.pioneta = pionMom.PseudoRapidity();
          fD0BgTree.pionphi = standardPhi(pionMom.Phi());
          fD0BgTree.pioncharge = charge1;
          fD0BgTree.piondca = (piontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
          fD0BgTree.nsigmapion = (piontrack == itrk1) ? trk1->nSigmaPion() : trk2->nSigmaPion();

          fD0BgTree.kaonpt = kaonMom.Perp();
          fD0BgTree.kaoneta = kaonMom.PseudoRapidity();
          fD0BgTree.kaonphi = standardPhi(kaonMom.Phi());
          fD0BgTree.kaoncharge = charge2;
          fD0BgTree.kaondca = (kaontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
          fD0BgTree.nsigmakaon = (kaontrack == itrk1) ? trk1->nSigmaKaon() : trk2->nSigmaKaon();

          d0bgtree->Fill();

        }
      }

      if (fTagDStarEvents){

        if ((mass < fInvMassSignal1) || (mass > fInvMassSignal2)) continue; // Only tracks between 1.8, 1.92

        TLorentzVector D0;
        D0.SetPxPyPzE(mResMom.x(), mResMom.y(), mResMom.z(), mass);

        TLorentzVector K;

        if (kaontrack == itrk1) K.SetPxPyPzE(mTrk1Mom.x(), mTrk1Mom.y(), mTrk1Mom.z(), Mkaon);
        else K.SetPxPyPzE(mTrk2Mom.x(), mTrk2Mom.y(), mTrk2Mom.z(), Mkaon);

        K.Boost(-mResMom);

        double costheta = (mResMom * K.Vect())/(mResMom.Mag()*K.Vect().Mag());

        // if (costheta > 0.8) continue;

        for(unsigned short itrk3 = 0; itrk3 < ntracks; itrk3++){
          if (itrk3 == itrk1 || itrk3 == itrk2) continue;

          StPicoTrack *trk3 = static_cast<StPicoTrack*>(mPicoDst->track(itrk3));
          if(!trk3){ continue; }

          TVector3 mTrk3Mom;

          if(doUsePrimTracks) {
            mTrk3Mom = trk3->pMom();
          } else {
            mTrk3Mom = trk3->gMom(mVertex, Bfield);
          }

          short charge3 = trk3->charge();

          if (!IsSoftPion(trk3)) continue;

          short chargeD0Pi = (piontrack == itrk1) ? trk1->charge() : trk2->charge();

          double softness = mResMom.Perp()/mTrk3Mom.Perp();
          if ( softness < 7.0 || softness > 20.0 ) continue;

          double energyD0 = TMath::Sqrt(pow(mResMom.Mag(),2) + pow(mass, 2));
          double energysPion = TMath::Sqrt(pow(mTrk3Mom.Mag(),2) + pow(Mpion, 2));

          double massDStar = TMath::Sqrt(pow(mass,2) + pow(Mpion,2) + 2*(energyD0*energysPion - mResMom*mTrk3Mom));

          double massDstarminusmassD0 = massDStar - mass;

          if (massDstarminusmassD0 > 0.2) continue;

          TVector3 mDStarMom;

          mDStarMom = mResMom + mTrk3Mom;

          if (chargeD0Pi*charge3 == 1){ // Choosing Pions with the same charge as the Pion From D0

            if ((mass >= fInvMassSignal1) && (mass <= fInvMassSignal2)){
              invmassdiff->Fill(massDstarminusmassD0);
              invmassDStar->Fill(massDStar);

              invmassdiffDStarvpt->Fill(massDstarminusmassD0, mDStarMom.Perp());
              invmassDStarvpt->Fill(massDstarminusmassD0, mDStarMom.Perp());

              fdstarTrackIndices.push_back({piontrack, kaontrack, itrk3, mass, massDstarminusmassD0});

              fdstar = kTRUE;

              fDStarSigTree = {0};

              fDStarSigTree.runid = mPicoEvent->runId();
              fDStarSigTree.eventid = mPicoEvent->eventId();
              fDStarSigTree.refmult = frefCorr2;
              fDStarSigTree.centrality = fCentralityScaled;
              fDStarSigTree.triggers = mPicoEvent->triggerIds();

              vector<double> pv;
              pv.clear();
              pv.push_back(mPicoEvent->primaryVertex().X());
              pv.push_back(mPicoEvent->primaryVertex().Y());
              pv.push_back(mPicoEvent->primaryVertex().Z());

              fDStarSigTree.primaryvertex = pv;

              vector<double> pverr;
              pverr.clear();
              pverr.push_back(mPicoEvent->primaryVertexError().X());
              pverr.push_back(mPicoEvent->primaryVertexError().Y());
              pverr.push_back(mPicoEvent->primaryVertexError().Z());

              fDStarSigTree.primaryvertexerror = pverr;

              fDStarSigTree.d0mass = mass;
              fDStarSigTree.d0pt = mResMom.Perp();
              fDStarSigTree.d0eta = mResMom.PseudoRapidity();
              fDStarSigTree.d0phi = standardPhi(mResMom.Phi());

              fDStarSigTree.dstarmass = massDStar;
              fDStarSigTree.dstarpt = mDStarMom.Perp();
              fDStarSigTree.dstareta = mDStarMom.PseudoRapidity();
              fDStarSigTree.dstarphi = standardPhi(mDStarMom.Phi());

              fDStarSigTree.pionpt = pionMom.Perp();
              fDStarSigTree.pioneta = pionMom.PseudoRapidity();
              fDStarSigTree.pionphi = standardPhi(pionMom.Phi());
              fDStarSigTree.pioncharge = charge1;
              fDStarSigTree.piondca = (piontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
              fDStarSigTree.nsigmapion = (piontrack == itrk1) ? trk1->nSigmaPion() : trk2->nSigmaPion();

              fDStarSigTree.kaonpt = kaonMom.Perp();
              fDStarSigTree.kaoneta = kaonMom.PseudoRapidity();
              fDStarSigTree.kaonphi = standardPhi(kaonMom.Phi());
              fDStarSigTree.kaoncharge = charge2;
              fDStarSigTree.kaondca = (kaontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
              fDStarSigTree.nsigmakaon = (kaontrack == itrk1) ? trk1->nSigmaKaon() : trk2->nSigmaKaon();

              fDStarSigTree.spionpt = mTrk3Mom.Perp();
              fDStarSigTree.spioneta = mTrk3Mom.PseudoRapidity();
              fDStarSigTree.spionphi = standardPhi(mTrk3Mom.Phi());
              fDStarSigTree.spioncharge = charge3;
              fDStarSigTree.spiondca = trk3->gDCA(mVertex).Mag();
              fDStarSigTree.snsigmapion = trk3->nSigmaPion();

              dstarsigtree->Fill();
            }

            else{
              invmassdiffside->Fill(massDstarminusmassD0);
              invmassdiffsideDStarvpt->Fill(massDstarminusmassD0, mDStarMom.Perp());
            }

          }

          else if (chargeD0Pi*charge3 == -1){

            fdstarLS = kTRUE;

            fdstarLSTrackIndices.push_back({piontrack, kaontrack, itrk3, mass, massDstarminusmassD0});

            invmassbgDStar->Fill(massDStar);
            invmassbgDStarvpt->Fill(massDStar, mDStarMom.Perp());

            invmassdiffwrong->Fill(massDstarminusmassD0);
            invmassdiffwrongDStarvpt->Fill(massDstarminusmassD0, mDStarMom.Perp());

            fDStarBgTree = {0};

            fDStarBgTree.runid = mPicoEvent->runId();
            fDStarBgTree.eventid = mPicoEvent->eventId();
            fDStarBgTree.refmult = frefCorr2;
            fDStarBgTree.centrality = fCentralityScaled;
            fDStarBgTree.triggers = mPicoEvent->triggerIds();

            vector<double> pv;
            pv.clear();
            pv.push_back(mPicoEvent->primaryVertex().X());
            pv.push_back(mPicoEvent->primaryVertex().Y());
            pv.push_back(mPicoEvent->primaryVertex().Z());

            fDStarBgTree.primaryvertex = pv;

            vector<double> pverr;
            pverr.clear();
            pverr.push_back(mPicoEvent->primaryVertexError().X());
            pverr.push_back(mPicoEvent->primaryVertexError().Y());
            pverr.push_back(mPicoEvent->primaryVertexError().Z());

            fDStarBgTree.primaryvertexerror = pverr;

            fDStarBgTree.d0mass = mass;
            fDStarBgTree.d0pt = mResMom.Perp();
            fDStarBgTree.d0eta = mResMom.PseudoRapidity();
            fDStarBgTree.d0phi = standardPhi(mResMom.Phi());

            fDStarBgTree.dstarmass = massDStar;
            fDStarBgTree.dstarpt = mDStarMom.Perp();
            fDStarBgTree.dstareta = mDStarMom.PseudoRapidity();
            fDStarBgTree.dstarphi = standardPhi(mDStarMom.Phi());

            fDStarBgTree.pionpt = pionMom.Perp();
            fDStarBgTree.pioneta = pionMom.PseudoRapidity();
            fDStarBgTree.pionphi = standardPhi(pionMom.Phi());
            fDStarBgTree.pioncharge = charge1;
            fDStarBgTree.piondca = (piontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
            fDStarBgTree.nsigmapion = (piontrack == itrk1) ? trk1->nSigmaPion() : trk2->nSigmaPion();

            fDStarBgTree.kaonpt = kaonMom.Perp();
            fDStarBgTree.kaoneta = kaonMom.PseudoRapidity();
            fDStarBgTree.kaonphi = standardPhi(kaonMom.Phi());
            fDStarBgTree.kaoncharge = charge2;
            fDStarBgTree.kaondca = (kaontrack == itrk1) ? trk1->gDCA(mVertex).Mag() : trk2->gDCA(mVertex).Mag();
            fDStarBgTree.nsigmakaon = (kaontrack == itrk1) ? trk1->nSigmaKaon() : trk2->nSigmaKaon();

            fDStarBgTree.spionpt = mTrk3Mom.Perp();
            fDStarBgTree.spioneta = mTrk3Mom.PseudoRapidity();
            fDStarBgTree.spionphi = standardPhi(mTrk3Mom.Phi());
            fDStarBgTree.spioncharge = charge3;
            fDStarBgTree.spiondca = trk3->gDCA(mVertex).Mag();
            fDStarBgTree.snsigmapion = trk3->nSigmaPion();

            dstarbgtree->Fill();
          }
        } //Trk3
      }
    } // Trk2
  } // Trk1
}

// void StTagDStarEvents::TestTracks(){

//   if (fd0TrackIndices.size() != 0){
//     cout << "Testing D0 Tracks" << endl;
//     for (int ican = 0; ican < fd0TrackIndices.size(); ican++){
//       StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(fd0TrackIndices[ican][0]));
//       StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(fd0TrackIndices[ican][1]));

//       TVector3 mTrk1Mom, mTrk2Mom, mResMom;

//       if(doUsePrimTracks) {
//         mTrk1Mom = trk1->pMom();
//       } else {
//         mTrk1Mom = trk1->gMom(mVertex, Bfield);
//       }

//       if(doUsePrimTracks) {
//         mTrk2Mom = trk2->pMom();
//       } else {
//         mTrk2Mom = trk2->gMom(mVertex, Bfield);
//       }

//       mResMom = mTrk1Mom + mTrk2Mom;

//       Double_t pt_pion = mTrk1Mom.Perp();
//       Double_t phi_pion = standardPhi(mTrk1Mom.Phi());
//       Double_t eta_pion = mTrk1Mom.PseudoRapidity();

//       Double_t pt_kaon = mTrk2Mom.Perp();
//       Double_t phi_kaon = standardPhi(mTrk2Mom.Phi());
//       Double_t eta_kaon = mTrk2Mom.PseudoRapidity();

//       Double_t phi_mom = standardPhi(mResMom.Phi());
//       Double_t eta_mom = mResMom.PseudoRapidity();
//       Double_t pt_mom = mResMom.Perp();

//       Double_t mass = fd0TrackIndices[ican][2];

//       cout << "Pion : " << pt_pion << "\t" << phi_pion << "\t" << eta_pion << endl;
//       cout << "Kaon : " << pt_kaon << "\t" << phi_kaon << "\t" << eta_kaon << endl;
//       cout << "D0 : " << pt_mom << "\t" << phi_mom << "\t" << eta_mom << endl;
//     }
//   }

//   if (fd0BgUSTrackIndices.size() != 0){
//     cout << "Testing D0 Bg US Tracks" << endl;
//     for (int ican = 0; ican < fd0BgUSTrackIndices.size(); ican++){
//       StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(fd0BgUSTrackIndices[ican][0]));
//       StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(fd0BgUSTrackIndices[ican][1]));

//       TVector3 mTrk1Mom, mTrk2Mom, mResMom;

//       if(doUsePrimTracks) {
//         mTrk1Mom = trk1->pMom();
//       } else {
//         mTrk1Mom = trk1->gMom(mVertex, Bfield);
//       }

//       if(doUsePrimTracks) {
//         mTrk2Mom = trk2->pMom();
//       } else {
//         mTrk2Mom = trk2->gMom(mVertex, Bfield);
//       }

//       mResMom = mTrk1Mom + mTrk2Mom;

//       Double_t pt_pion = mTrk1Mom.Perp();
//       Double_t phi_pion = standardPhi(mTrk1Mom.Phi());
//       Double_t eta_pion = mTrk1Mom.PseudoRapidity();

//       Double_t pt_kaon = mTrk2Mom.Perp();
//       Double_t phi_kaon = standardPhi(mTrk2Mom.Phi());
//       Double_t eta_kaon = mTrk2Mom.PseudoRapidity();

//       Double_t phi_mom = standardPhi(mResMom.Phi());
//       Double_t eta_mom = mResMom.PseudoRapidity();
//       Double_t pt_mom = mResMom.Perp();

//       Double_t mass = fd0BgUSTrackIndices[ican][2];

//       cout << "Pion : " << pt_pion << "\t" << phi_pion << "\t" << eta_pion << endl;
//       cout << "Kaon : " << pt_kaon << "\t" << phi_kaon << "\t" << eta_kaon << endl;
//       cout << "D0 Bg US: " << pt_mom << "\t" << phi_mom << "\t" << eta_mom << endl;
//     }
//   }

//   if (fd0BgLSTrackIndices.size() != 0){
//     cout << "Testing D0 Bg LS Tracks" << endl;
//     for (int ican = 0; ican < fd0BgLSTrackIndices.size(); ican++){
//       StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(fd0BgLSTrackIndices[ican][0]));
//       StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(fd0BgLSTrackIndices[ican][1]));

//       TVector3 mTrk1Mom, mTrk2Mom, mResMom;

//       if(doUsePrimTracks) {
//         mTrk1Mom = trk1->pMom();
//       } else {
//         mTrk1Mom = trk1->gMom(mVertex, Bfield);
//       }

//       if(doUsePrimTracks) {
//         mTrk2Mom = trk2->pMom();
//       } else {
//         mTrk2Mom = trk2->gMom(mVertex, Bfield);
//       }

//       mResMom = mTrk1Mom + mTrk2Mom;

//       Double_t pt_pion = mTrk1Mom.Perp();
//       Double_t phi_pion = standardPhi(mTrk1Mom.Phi());
//       Double_t eta_pion = mTrk1Mom.PseudoRapidity();

//       Double_t pt_kaon = mTrk2Mom.Perp();
//       Double_t phi_kaon = standardPhi(mTrk2Mom.Phi());
//       Double_t eta_kaon = mTrk2Mom.PseudoRapidity();

//       Double_t phi_mom = standardPhi(mResMom.Phi());
//       Double_t eta_mom = mResMom.PseudoRapidity();
//       Double_t pt_mom = mResMom.Perp();

//       Double_t mass = fd0BgLSTrackIndices[ican][2];

//       cout << "Pion : " << pt_pion << "\t" << phi_pion << "\t" << eta_pion << endl;
//       cout << "Kaon : " << pt_kaon << "\t" << phi_kaon << "\t" << eta_kaon << endl;
//       cout << "D0 Bg LS: " << pt_mom << "\t" << phi_mom << "\t" << eta_mom << endl;
//     }
//   }
// }

Bool_t StTagDStarEvents::IsAnAcceptableTrack(StPicoTrack *trk, bool dohistograms = kFALSE){

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, Bfield);
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

  if (charge == 0) return kFALSE;

  if (dohistograms){ hAcceptedPt[0]->Fill(pt); hAcceptedEta[0]->Fill(eta); }

  if (pt < 0.2) return kFALSE;
  
  if (dohistograms){ hAcceptedPt[1]->Fill(pt); hAcceptedEta[1]->Fill(eta); }

  if (abs(eta) > 1) return kFALSE;

  if (dohistograms){ hAcceptedPt[2]->Fill(pt); hAcceptedEta[2]->Fill(eta); }

  if (trk->nHitsDedx() < 15) return kFALSE;

  if (dohistograms){ hAcceptedPt[3]->Fill(pt); hAcceptedEta[3]->Fill(eta); }

  if (float(trk->nHitsDedx())/float(trk->nHitsMax()) < 0.52) return kFALSE;

  if (dohistograms){ hAcceptedPt[4]->Fill(pt); hAcceptedEta[4]->Fill(eta); }

  if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > fTrackDCAcut) return kFALSE;

  if (dohistograms){ hAcceptedPt[5]->Fill(pt); hAcceptedEta[5]->Fill(eta); }

  if (trk->chi2() > 3) return kFALSE;

  if (dohistograms){ hAcceptedPt[6]->Fill(pt); hAcceptedEta[6]->Fill(eta); }

  // HFT Hits only

  if(!doppAnalysis && fTopoLevel > 0){
    if (!trk->isHFTTrack()) return kFALSE;
    if (!trk->hasPxl1Hit() || !trk->hasPxl2Hit() || !trk->hasIstHit()) return kFALSE;
  }

  if (dohistograms){ hAcceptedPt[7]->Fill(pt); hAcceptedEta[7]->Fill(eta); }

  // bool bemctrack = trk->isBemcMatchedTrack();
  // if ((!bemctrack)) return kFALSE;

  return kTRUE;
}


Int_t StTagDStarEvents::TopologicalCuts(StPicoTrack *trk1, StPicoTrack *trk2, bool dohistograms = kFALSE){

  int cutnum = -99;

  int pid1, pid2;
  double m1, m2;
  double e1, e2;

  IsWhatParticle(trk1, pid1, m1, e1);
  IsWhatParticle(trk2, pid2, m2, e2);

  if (abs(pid1*pid2) != 2) return cutnum;

  TVector3 mTrk1Mom;
  if(doUsePrimTracks) {
    mTrk1Mom = trk1->pMom();
  } else {
    mTrk1Mom = trk1->gMom(mVertex, Bfield);
  }

  TVector3 mTrk2Mom;
  if(doUsePrimTracks) {
    mTrk2Mom = trk2->pMom();
  } else {
    mTrk2Mom = trk2->gMom(mVertex, Bfield);
  }

  // to be used for testing with preview II pico production
  StPicoPhysicalHelix kHelix = trk1->helix(Bfield);
  StPicoPhysicalHelix pHelix = trk2->helix(Bfield);

  // move origins of helices to the primary vertex origin
  kHelix.moveOrigin(kHelix.pathLength(mVertex));
  pHelix.moveOrigin(pHelix.pathLength(mVertex));

  // use straight lines approximation to get point of DCA of kaon-pion pair
  TVector3 const kMom = kHelix.momentum(Bfield*(1.e-14));//*(1.e-14);
  TVector3 const pMom = pHelix.momentum(Bfield*(1.e-14));//*(1.e-14);

  // cout << mTrk1Mom.Mag() << "\t" << kMom.Mag() << "\t" << mTrk2Mom.Mag() << "\t" << pMom.Mag() << endl;

  StPicoPhysicalHelix const kStraightLine(kMom, kHelix.origin(), 0, trk1->charge());
  StPicoPhysicalHelix const pStraightLine(pMom, pHelix.origin(), 0, trk2->charge());


  pair<double, double> const ss = kStraightLine.pathLengths(pStraightLine);
  TVector3 const kAtDcaToPion = kStraightLine.at(ss.first);
  TVector3 const pAtDcaToKaon = pStraightLine.at(ss.second);

  // calculate DCA of pion to kaon at their DCA
  Double_t mDcaDaughters = (kAtDcaToPion - pAtDcaToKaon).Mag();

  // calculate Lorentz vector of kaon-pion pair
  TVector3 const kMomAtDca = kHelix.momentumAt(ss.first, Bfield*(1.e-14)); //*(1.e-14);
  TVector3 const pMomAtDca = pHelix.momentumAt(ss.second, Bfield*(1.e-14)); //*(1.e-14);

  TLorentzVector const kFourMom(kMomAtDca, e1);
  TLorentzVector const pFourMom(pMomAtDca, e2);

  TLorentzVector mLorentzVector = kFourMom + pFourMom;

  // calculate cosThetaStar
  //TLorentzVector const kpFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
  //TLorentzVector const kFourMomStar = kFourMom.Boost(kpFourMomReverse);
  //mCosThetaStar = std::cos(kFourMomStar.Vect().Angle(mLorentzVector.Vect()));

  // calculate pointing angle and decay length
  TVector3 const vtxToV0 = (kAtDcaToPion + pAtDcaToKaon) * 0.5 - mVertex;
  Double_t mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());

  //mPointingAngle = vtxToV0.Angle(finalFourVector.vect());

  Double_t mDecayLength = vtxToV0.Mag();

  Double_t perpDcaToVtx = mDecayLength*TMath::Sin(mPointingAngle);

  // calculate DCA of tracks to primary vertex
  Double_t mKaonDca = (kHelix.origin() - mVertex).Mag();
  Double_t mPionDca = (pHelix.origin() - mVertex).Mag();

  if (dohistograms && pid1*pid2==-2) decaylengthd0US->Fill(mDecayLength);
  if (dohistograms && pid1*pid2==-2) distancepikUS->Fill(mDcaDaughters);
  if (dohistograms && pid1*pid2==-2) distanced0PVUS->Fill(perpDcaToVtx);
  if (dohistograms && pid1*pid2==-2) dcakPVUS->Fill(mKaonDca);
  if (dohistograms && pid1*pid2==-2) dcapiPVUS->Fill(mPionDca);

  if (dohistograms && pid1*pid2==2) decaylengthd0LS->Fill(mDecayLength);
  if (dohistograms && pid1*pid2==2) distancepikLS->Fill(mDcaDaughters);
  if (dohistograms && pid1*pid2==2) distanced0PVLS->Fill(perpDcaToVtx);
  if (dohistograms && pid1*pid2==2) dcakPVLS->Fill(mKaonDca);
  if (dohistograms && pid1*pid2==2) dcapiPVLS->Fill(mPionDca);

  // if (mDecayLength < 0.0212) return kFALSE;
  // if (mDcaDaughters > 0.0057) return kFALSE;
  // if (perpDcaToVtx > 0.0038) return kFALSE;
  // if (mKaonDca < 0.0095) return kFALSE;
  // if (mPionDca < 0.0086) return kFALSE;

  if ((mDecayLength > 0.0212) && (mDcaDaughters < 0.0057) && (perpDcaToVtx < 0.0038) && (mKaonDca > 0.0095) && (mPionDca > 0.0086)){
    cutnum = 1;
    return cutnum;
  }

  if ((mDecayLength > 0.5*0.0212) && (mDcaDaughters < 2*0.0057) && (perpDcaToVtx < 2*0.0038) && (mKaonDca > 0.5*0.0095) && (mPionDca > 0.5*0.0086)){
    cutnum = 2;
    return cutnum;
  }

  if ((mDecayLength > 0.25*0.0212) && (mDcaDaughters < 4*0.0057) && (perpDcaToVtx < 4*0.0038) && (mKaonDca > 0.25*0.0095) && (mPionDca > 0.25*0.0086)){
    cutnum = 3;
    return cutnum;
  }

  return cutnum;
  // cout << "Topo Cuts: " << mDecayLength << "\t" << mDcaDaughters << "\t" << perpDcaToVtx << "\t" << mKaonDca << "\t" << mPionDca << endl;
} 


void StTagDStarEvents::FillPidHistograms(StPicoTrack *trk){

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
    double invbeta_from_tof = tofpointer->btofBeta();
    invbeta_from_tof = 1/invbeta_from_tof;

    invbetavp->Fill(p, invbeta_from_tof);

    double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
    double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
    normalised_invbetavpT_tof_pi->Fill(pt, normalisedinvbeta_for_pi);
    double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
    double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.011;
    normalised_invbetavpT_tof_ka->Fill(pt, normalisedinvbeta_for_ka);
  }

  else{
      double zpi = trk->nSigmaPion();
      z_pi->Fill(pt, zpi);
      double zka = trk->nSigmaKaon();
      z_ka->Fill(pt, zka);
  }
}

bool StTagDStarEvents::IsSoftPion(StPicoTrack *trk){
  // if(!IsAnAcceptableTrack(trk, kFALSE)){return kFALSE;}

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, Bfield);
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

  double dca = trk->gDCA(mVertex).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // track acceptance cuts now - after getting 3vector - hardcoded
  if(pt < 0.15) return kFALSE;
  if(pt > fTrackPtMaxCut) return kFALSE; // 20.0 STAR, (increased to 30.0) 100.0 ALICE
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  if((phi < fTrackPhiMinCut) || (phi > fTrackPhiMaxCut)) return kFALSE;

  if(dca > fTrackDCAcut)            return kFALSE;
  if(nHitsFit < fTracknHitsFit)     return kFALSE;
  if(nHitsRatio < fTracknHitsRatio) return kFALSE;

  double dedx = trk->dEdx();
  double dedxresolution = trk->dEdxError();

  // bichsel function approximation
  double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
  double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
  double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

  //z - variables
  // double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
  // double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
  // double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

  double zpi = trk->nSigmaPion();
  double zka = trk->nSigmaKaon();
  double zpr = trk->nSigmaProton();


  if (abs(zpi) < 3) return kTRUE;
  else return kFALSE;

}

void StTagDStarEvents::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){ // NEW PID APPROACH
  if(!IsAnAcceptableTrack(trk, kFALSE)){pid = 0; m = 0.; e = 0.; return;}

  pid = 0;
  m = 0.0;
  e = 0.0;

  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    mTrkMom = trk->pMom();
  } else {
    mTrkMom = trk->gMom(mVertex, Bfield);
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

  double dedx = trk->dEdx();
  double dedxresolution = trk->dEdxError();

  // bichsel function approximation
  double dedxth_pi = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mpion)));
  double dedxth_ka = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mkaon)));
  double dedxth_pr = charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(p/Mproton)));

  //z - variables
  // double zpi = TMath::Log(dedx/dedxth_pi)/dedxresolution;
  // double zka = TMath::Log(dedx/dedxth_ka)/dedxresolution;
  // double zpr = TMath::Log(dedx/dedxth_pr)/dedxresolution;

  double zpi = trk->nSigmaPion();
  double zka = trk->nSigmaKaon();
  double zpr = trk->nSigmaProton();

  if ((abs(zpi) > 3) && (abs(zka) > 3)) {pid = 0; m = 0.; e = 0.; return;} //Wider because this is in pp DStarAnalysis Only

  bool tpc_pion = kFALSE;
  bool tpc_kaon = kFALSE;

  tpc_pion = (abs(zpi) <= abs(zka)) ? kTRUE : kFALSE;
  tpc_kaon = (abs(zpi) > abs(zka)) ? kTRUE : kFALSE;

  bool toftrack = kFALSE;

  int tof_loc = trk->bTofPidTraitsIndex();

  if (tof_loc >= 0)
  {
    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    if (tofpointer){toftrack = kTRUE;}
  }

  if(toftrack){ //Only For Kaons

    bool tof_pion = kFALSE;
    bool tof_kaon = kFALSE;

    StPicoBTofPidTraits *tofpointer = static_cast<StPicoBTofPidTraits*>(mPicoDst->btofPidTraits(tof_loc));
    double invbeta_from_tof = tofpointer->btofBeta();
    invbeta_from_tof = 1/invbeta_from_tof;

    double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
    double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
    double norm_invbeta_pr = TMath::Sqrt(pow(Mproton,2)/pow(p,2) + 1);

    double normalisedinvbeta_for_pi = (invbeta_from_tof-norm_invbeta_pi)/0.011;
    double normalisedinvbeta_for_ka = (invbeta_from_tof-norm_invbeta_ka)/0.011;
    double normalisedinvbeta_for_pr = (invbeta_from_tof-norm_invbeta_pr)/0.011;

    TF1 *PSigmaBand = new TF1("PSigmaBand", "2.0*([0] + [1]/pow((x + [2]), [3]))", 0.2, 10.);
    PSigmaBand->SetParameters(1.015, 0.01884, 0.03166, 3.469);

    if (normalisedinvbeta_for_ka > PSigmaBand->Eval(pt)*(-1.) && normalisedinvbeta_for_ka < PSigmaBand->Eval(pt)*(1.)) tof_kaon = kTRUE;

    if (tof_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
  }

  else{
    if (tpc_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tpc_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}
  }
}


Int_t StTagDStarEvents::IsWhatParticle(StPicoTrack *trk){ // Just to get the PID out
  int pid;
  double m; 
  double e;
  IsWhatParticle(trk, pid, m, e);
  return pid;
}

void StTagDStarEvents::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum){

  int particle1, particle2;
  double mass1, mass2;
  double energy1, energy2;

  // cout << "IsWhatParticle called" << endl;

  IsWhatParticle(trk1, particle1, mass1, energy1);
  IsWhatParticle(trk2, particle2, mass2, energy2);

  TVector3 mTrk1Mom, mTrk2Mom;

  if(doUsePrimTracks) {
    mTrk1Mom = trk1->pMom();
  } else {
    mTrk1Mom = trk1->gMom(mVertex, Bfield);
  }

  if(doUsePrimTracks) {
    mTrk2Mom = trk2->pMom();
  } else {
    mTrk2Mom = trk2->gMom(mVertex, Bfield);
  }

  if ((mass1<0.000001) || (mass2<0.000001) || (energy1<0.000001) || (energy2<0.000001)){
    particle=0;
    invmass=0.;
    momentum = mTrk1Mom + mTrk2Mom;
    return;
  }
  else{
    particle = particle1*particle2;
    invmass = TMath::Sqrt(pow(mass1,2) + pow(mass2,2) + 2*(energy1*energy2 - mTrk1Mom*mTrk2Mom));
    momentum = mTrk1Mom + mTrk2Mom;
  }
}

Double_t StTagDStarEvents::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2){
  int particle; 
  double invmass; 
  TVector3 momentum;

  InvariantMass(trk1, trk2, particle, invmass, momentum);
  return invmass;
}


// Centrality bin getter
Int_t StTagDStarEvents::GetFourCentBin(Double_t scaledCent) const{
  int centbin = -99;
  // get centrality bin number
  if(scaledCent >= 0 && scaledCent <  10.0)  { centbin = 0; }
  else if(scaledCent >= 10.0 && scaledCent <  20.0)                { centbin = 1; }
  else if(scaledCent >= 20.0 && scaledCent <  50.0)                { centbin = 2; }
  else if(scaledCent >= 50.0 && scaledCent <= 80.0)                { centbin = 3; }

  return centbin;
}

