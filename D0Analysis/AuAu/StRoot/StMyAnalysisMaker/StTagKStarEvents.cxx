// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StTagKStarEvents.h"
#include "StMemStat.h"
#include "phys_constants.h"
#include <limits>
#include "math.h"
// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "THn.h"
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

ClassImp(StTagKStarEvents)

//________________________________________________________________________
StTagKStarEvents::StTagKStarEvents(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTagKStarEvents::fRunFlagEnum
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
  fInvMassSignal1 = 1.70;
  fInvMassSignal2 = 2.10;

  fInvMassULBg1 = 1.70;
  fInvMassULBg2 = 2.10;

  fInvMassLSBg1 = 1.70;
  fInvMassLSBg2 = 2.10;

  fd0 = kFALSE;
  fd0Bg = kFALSE;

  fTightd0 = kFALSE;
  fTightd0Bg = kFALSE;

  fLoosed0 = kFALSE;
  fLoosed0Bg = kFALSE;

  fdoQA_Histograms = kTRUE;
  fdoDaug_Histograms = kTRUE;

  fdoSingleParticleEmbedding = kFALSE;
  fSingleParticle.SetPxPyPzE(0, 0, 0, 0);

  // fd0TrackIndices.clear();
  // fd0BgUSTrackIndices.clear();
  // fd0BgLSTrackIndices.clear();
}

//
//________________________________________________________________________
StTagKStarEvents::~StTagKStarEvents()
{ 

  if (fdoQA_Histograms){

    if (hEventZvertex_diff) delete hEventZvertex_diff;
    if (hEventVzvsVzvpd) delete hEventVzvsVzvpd;
    if (hEventVxvsVy) delete hEventVxvsVy;

    if (hCentralityBeforeCuts) delete hCentralityBeforeCuts;
    if (hCentralityAfterCuts) delete hCentralityAfterCuts;
    if (hCentralityWeightedBeforeCuts) delete hCentralityWeightedBeforeCuts;
    if (hCentralityWeightedAfterCuts) delete hCentralityWeightedAfterCuts;
    if (hRefMultiplicity) delete hRefMultiplicity;
    if (hgRefMultiplicity) delete hgRefMultiplicity;

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
  }

  if (fdoDaug_Histograms){

    for (int i = 0 ; i < 3; i++){
      if (decaylengthd0US[i]) delete decaylengthd0US[i];
      if (distancepikUS[i]) delete distancepikUS[i];
      if (distanced0PVUS[i]) delete distanced0PVUS[i];
      if (dcakPVUS[i]) delete dcakPVUS[i];
      if (dcapiPVUS[i]) delete dcapiPVUS[i];
      if (costhetaDVPVUS[i]) delete costhetaDVPVUS[i];

      if (decaylengthd0LS[i]) delete decaylengthd0LS[i];
      if (distancepikLS[i]) delete distancepikLS[i];
      if (distanced0PVLS[i]) delete distanced0PVLS[i];
      if (dcakPVLS[i]) delete dcakPVLS[i];
      if (dcapiPVLS[i]) delete dcapiPVLS[i];
      if (costhetaDVPVLS[i]) delete costhetaDVPVLS[i];
    
      if (hD0CentPtEtaMDphiDaug[i]) delete hD0CentPtEtaMDphiDaug[i];
      if (hD0CentPtEtaMDphiDaugLikeSign[i]) delete hD0CentPtEtaMDphiDaugLikeSign[i];
    }

  }

  if(mEmcPosition) delete mEmcPosition;

}

//________________________________________________________________________
Int_t StTagKStarEvents::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTagKStarEvents::Finish() { 

  cout << "StTagKStarEvents::Finish()\n";

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

  cout<<"End of StTagKStarEvents::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StTagKStarEvents::DeclareHistograms() {

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

  if (fdoQA_Histograms){
    hEventZvertex_diff = new TH1F("hEventZvertex_diff", "z-vertex distribution Difference", 200, -100., 100.);
    hEventVzvsVzvpd = new TH2F("hEventVzvsVzvpd", "Vz vs VPD Vz", 200, -10, 10, 200, -10, 10);
    hEventVxvsVy = new TH2F("hEventVxvsVy", "Vx vs Vy", 200, -10, 10, 200, -10, 10);

    hCentralityBeforeCuts = new TH1F("hCentralityBeforeCuts", "No. events vs centrality (Before Cuts)", nHistCentBins, 0, 100);
    hCentralityAfterCuts = new TH1F("hCentralityAfterCuts", "No. events vs centrality (After Cuts)", nHistCentBins, 0, 100);
    hCentralityWeightedBeforeCuts = new TH1F("hCentralityWeightedBeforeCuts", "No. events vs centrality weighted (Before Cuts)", nHistCentBins, 0, 100);
    hCentralityWeightedAfterCuts = new TH1F("hCentralityWeightedAfterCuts", "No. events vs centrality weighted (After Cuts)", nHistCentBins, 0, 100);
    hRefMultiplicity = new TH1F("hRefMultiplicity", "No. events vs refmultiplicity", kHistMultBins, 0, kHistMultMax);
    hgRefMultiplicity = new TH1F("hgRefMultiplicity", "No. events vs grefmultiplicity", kHistMultBins, 0, kHistMultMax);

    cuthistogram_event = new TH1F("cuthistogram_event", "cuthistogram_event", 14, 0, 14);

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
  }

  if (fdoDaug_Histograms){
    TString cutlevel[3] = {"Standard", "50percent", "150percent"};

    const int nDimDaug = 6;

    const int nBinsCent = 9;
    double CentBins[nBinsCent + 1]           = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const int nBinsD0Pt = 22;
    double D0PtBins[nBinsD0Pt + 1]           = {-10.0, -8.0, -6.0, -5.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    const int nBinsDaugPt = 4;
    double DaugPtBins[nBinsDaugPt + 1]       = {0.3, 0.6, 1.2, 5.0, 10.0};
    const int nBinsRapidity = 4;
    double D0RapidityBins[nBinsRapidity + 1] = {-1.0, -0.5, 0.0, 0.5, 1.0};

    const int nBinsMass = 100;
    double D0MassBins[nBinsMass + 1];
    double lowmasslimit = fInvMassSignal1;
    double highmasslimit = fInvMassSignal2;

    for (int i = 0; i <= nBinsMass; i++){
      D0MassBins[i] = lowmasslimit + (highmasslimit - lowmasslimit)*float(i)/float(nBinsMass);
    }
    int nBinsDaug[nDimDaug] = {nBinsCent, nBinsD0Pt, nBinsDaugPt, nBinsMass, nBinsDaugPt, nBinsRapidity}; //cent, pt*charge, daughterpt1, m, daughterpt2*charge, rapidity

    // Topological histograms

    for (int i = 0; i < 3; i++){
      decaylengthd0US[i] = new TH3F(Form("decaylengthd0US_%s", cutlevel[i].Data()), Form("decaylengthd0US_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      distancepikUS[i] = new TH3F(Form("distancepikUS_%s", cutlevel[i].Data()), Form("distancepikUS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      distanced0PVUS[i] = new TH3F(Form("distanced0PVUS_%s", cutlevel[i].Data()), Form("distanced0PVUS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      dcakPVUS[i] = new TH3F(Form("dcakPVUS_%s", cutlevel[i].Data()), Form("dcakPVUS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      dcapiPVUS[i] = new TH3F(Form("dcapiPVUS_%s", cutlevel[i].Data()), Form("dcapiPVUS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      costhetaDVPVUS[i] = new TH3F(Form("costhetaDVPVUS_%s", cutlevel[i].Data()), Form("costhetaDVPVUS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 500, -1, 1);

      decaylengthd0LS[i] = new TH3F(Form("decaylengthd0LS_%s", cutlevel[i].Data()), Form("decaylengthd0LS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      distancepikLS[i] = new TH3F(Form("distancepikLS_%s", cutlevel[i].Data()), Form("distancepikLS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      distanced0PVLS[i] = new TH3F(Form("distanced0PVLS_%s", cutlevel[i].Data()), Form("distanced0PVLS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      dcakPVLS[i] = new TH3F(Form("dcakPVLS_%s", cutlevel[i].Data()), Form("dcakPVLS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      dcapiPVLS[i] = new TH3F(Form("dcapiPVLS_%s", cutlevel[i].Data()), Form("dcapiPVLS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 300, 0, 0.3);
      costhetaDVPVLS[i] = new TH3F(Form("costhetaDVPVLS_%s", cutlevel[i].Data()), Form("costhetaDVPVLS_%s", cutlevel[i].Data()), 11, -0.5, 10.5, 20, 0, 10, 500, -1, 1);
    
      hD0CentPtEtaMDphiDaug[i] = new THnF(Form("hD0CentPtEtaMDphiDaug_%s", cutlevel[i].Data()), Form("hD0CentPtEtaMDphiDaug_%s", cutlevel[i].Data()), nDimDaug, nBinsDaug, NULL, NULL);
      hD0CentPtEtaMDphiDaugLikeSign[i] = new THnF(Form("hD0CentPtEtaMDphiDaugLikeSign_%s", cutlevel[i].Data()), Form("hD0CentPtEtaMDphiDaugLikeSign_%s", cutlevel[i].Data()), nDimDaug, nBinsDaug, NULL, NULL);

      hD0CentPtEtaMDphiDaug[i]->SetBinEdges(0, CentBins);
      hD0CentPtEtaMDphiDaug[i]->SetBinEdges(1, D0PtBins);
      hD0CentPtEtaMDphiDaug[i]->SetBinEdges(2, DaugPtBins);
      hD0CentPtEtaMDphiDaug[i]->SetBinEdges(3, D0MassBins);
      hD0CentPtEtaMDphiDaug[i]->SetBinEdges(4, DaugPtBins);
      hD0CentPtEtaMDphiDaug[i]->SetBinEdges(5, D0RapidityBins);

      hD0CentPtEtaMDphiDaugLikeSign[i]->SetBinEdges(0, CentBins);
      hD0CentPtEtaMDphiDaugLikeSign[i]->SetBinEdges(1, D0PtBins);
      hD0CentPtEtaMDphiDaugLikeSign[i]->SetBinEdges(2, DaugPtBins);
      hD0CentPtEtaMDphiDaugLikeSign[i]->SetBinEdges(3, D0MassBins);
      hD0CentPtEtaMDphiDaugLikeSign[i]->SetBinEdges(4, DaugPtBins);
      hD0CentPtEtaMDphiDaugLikeSign[i]->SetBinEdges(5, D0RapidityBins);
    }
  }
}

//
// write histograms
//_____________________________________________________________________________
void StTagKStarEvents::WriteHistograms() {

  if (fdoQA_Histograms){
    hEventZvertex_diff->Write();
    hEventVzvsVzvpd->Write();
    hEventVxvsVy->Write();

    hCentralityBeforeCuts->Write();
    hCentralityAfterCuts->Write();
    hCentralityWeightedBeforeCuts->Write();
    hCentralityWeightedAfterCuts->Write();
    hRefMultiplicity->Write();
    hgRefMultiplicity->Write();

    // Event Cuts Histograms
    const char *event_cuts[14] = {"Total Events", "Good Runs", "Max Track Pt < 30", "|V_{i}| != 0", "|V_{z}| < 30 cm", "VPDMB-5","|V_{z}| < 6 cm", "|V_{r}| < 2 cm ", "|V_{z} - V_{z(VPD)}| < 3 cm", "VPDMB-30", "nBEMCMatch && nBTOFMatch > 0", "HT1", "HT2", "HT3"};
    const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 20", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

    for (int i=1; i <= 14; i++){
      cuthistogram_event->SetBinContent(i, numberofevents[i-1]);
      cuthistogram_event->GetXaxis()->SetBinLabel(i, event_cuts[i-1]);
    }
    cuthistogram_event->Write();

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
  }

  // Topological Histograms

  if (fdoDaug_Histograms){
    for (int i = 0; i < 3; i++){
      decaylengthd0US[i]->Write();
      distancepikUS[i]->Write();
      distanced0PVUS[i]->Write();
      dcakPVUS[i]->Write();
      dcapiPVUS[i]->Write();
      costhetaDVPVUS[i]->Write();

      decaylengthd0LS[i]->Write();
      distancepikLS[i]->Write();
      distanced0PVLS[i]->Write();
      dcakPVLS[i]->Write();
      dcapiPVLS[i]->Write();
      costhetaDVPVLS[i]->Write();

      hD0CentPtEtaMDphiDaug[i]->Write();
      hD0CentPtEtaMDphiDaugLikeSign[i]->Write();
    }
  }  
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTagKStarEvents::Clear(Option_t *opt) {
  // fJets->Clear();
  fPionIndices.clear();
  fKaonIndices.clear();

  fd0TrackIndices.clear();
  fd0BgTrackIndices.clear();

  fTightd0TrackIndices.clear();
  fTightd0BgTrackIndices.clear();

  fLoosed0TrackIndices.clear();
  fLoosed0BgTrackIndices.clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTagKStarEvents::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  fd0 = kFALSE;
  fd0Bg = kFALSE;

  fTightd0 = kFALSE;
  fTightd0Bg = kFALSE;

  fLoosed0 = kFALSE;
  fLoosed0Bg = kFALSE;

  //zero out global vectors
  fPionIndices.clear();
  fKaonIndices.clear();

  //zero out global vectors
  fd0TrackIndices.clear();
  fd0BgTrackIndices.clear();

  fTightd0TrackIndices.clear();
  fTightd0BgTrackIndices.clear();

  fLoosed0TrackIndices.clear();
  fLoosed0BgTrackIndices.clear();

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

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // cout << "B = " << Bfield << endl;
  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  double zVtx_VPD = mPicoEvent->vzVpd();

  if (fdoQA_Histograms){
    hEventZvertex_diff->Fill(zVtx - zVtx_VPD);
    hEventVzvsVzvpd->Fill(zVtx, zVtx_VPD);
    hEventVxvsVy->Fill(mVertex.x(), mVertex.y());
  }

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() <= fMaxEventTrackPt){
    numberofevents[2]++;
  }

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
  double weight = mCentMaker->GetWeight();

  if (fdoQA_Histograms){
    hCentralityBeforeCuts->Fill(fCentralityScaled);
    hCentralityWeightedBeforeCuts->Fill(fCentralityScaled, weight);
  }

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5) return kStOK;

  numberofevents[3]++;

  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  numberofevents[4]++;




  if (!fMCEventsWithoutCent)
  {
    if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

    // cout << "Here" << endl;  

    // cut on centrality for analysis before doing anything
    if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

    // ============================ end of CENTRALITY ============================== //

    // ========================= Trigger Info =============================== //
    // looking at the EMCal triggers - used for QA and deciding on HT triggers
    //  FillEmcTriggers();

    // get trigger IDs from PicoEvent class and loop over them
    vector<unsigned int> mytriggers = mPicoEvent->triggerIds();

    int arrMB5_Run14[]  = {450005, 450015, 450025, 450050, 450060};
    int arrMB30_Run14[] = {450008, 450009, 450014, 450018, 450024};
    int arrHT1_Run14[]  = {450201, 450211, 460201};
    int arrHT2_Run14[]  = {450202, 450212, 460202, 460212};
    int arrHT3_Run14[]  = {450203, 450213, 460203};

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

    if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.) return kStOK;
    numberofevents[7]++;

    if (abs(zVtx - zVtx_VPD) > 3) return kStOk;
    numberofevents[8]++;

    if (!doppAnalysis){
      bool matchMB30 = kFALSE;

      for(int i = 0; i < sizeof(arrMB30_Run14)/sizeof(*arrMB30_Run14); i++) {
        if(mPicoEvent->isTrigger(arrMB30_Run14[i])) matchMB30 = kTRUE;
        if(matchMB30) break;
      }

      if (matchMB30) numberofevents[9]++;

    }
     
    if (mPicoEvent->nBEMCMatch() != 0 || mPicoEvent->nBTOFMatch() != 0) numberofevents[10]++;

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

    // fill histograms
    if (fdoQA_Histograms){
      hCentralityAfterCuts->Fill(fCentralityScaled);
      hCentralityWeightedAfterCuts->Fill(fCentralityScaled, weight);
      hRefMultiplicity->Fill(refMult);
      hgRefMultiplicity->Fill(grefMult);
    }
  }

  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  if (fDebugLevel == 1)cout << "Number of tracks found = " << ntracks << "\t" << nglobaltracks << endl;

  RunTracks();

  return kStOK;
}

//__________________________________________________________________________________________
  
void StTagKStarEvents::RunTracks(){

  fd0 = kFALSE;
  fd0Bg = kFALSE;

  fTightd0 = kFALSE;
  fTightd0Bg = kFALSE;

  fLoosed0 = kFALSE;
  fLoosed0Bg = kFALSE;

  fPionIndices.clear();
  fKaonIndices.clear();

  fd0TrackIndices.clear();
  fd0BgTrackIndices.clear();

  fTightd0TrackIndices.clear();
  fTightd0BgTrackIndices.clear();

  fLoosed0TrackIndices.clear();
  fLoosed0BgTrackIndices.clear();

  fSingleParticleVector.clear();

  const Int_t ntracks = mPicoDst->numberOfTracks();

  for(unsigned short itrk = 0; itrk < ntracks; itrk++){
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
    if(!trk){ continue; }

    if (!IsAnAcceptableTrack(trk, fdoQA_Histograms)) continue;

    if(fdoQA_Histograms) FillPidHistograms(trk);

    // If a track is both pion and kaon, it ends up in both the arrays! There needs to be a protection in the main loop. Run over pions first, kaons second.
    if (IsPion(trk)) fPionIndices.push_back(itrk);
    if (IsKaon(trk)) fKaonIndices.push_back(itrk);
  }

  // Main loop (Pion goes first, kaon goes later)
  for (int piontrack = 0; piontrack < fPionIndices.size(); piontrack++){
    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(fPionIndices[piontrack]));
    if(!trk1){ continue; }

    for(int kaontrack = 0; kaontrack < fKaonIndices.size(); kaontrack++){
      StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(fKaonIndices[kaontrack]));
      if(!trk2){ continue; }

      if (fPionIndices[piontrack] == fKaonIndices[kaontrack]) continue; //Avoid pairs where pion and kaon are the same tracks

      double mass = InvariantMass(trk1, trk2);
      if ( mass < fInvMassSignal1 || mass > fInvMassSignal2 ) continue;

      bool topocuts[3];

      if (fTopoLevel == 0){ // Suppose we don't need topo cuts. Shouldn't happen, ever!
        bool tmp;
        TopologicalCutsNew(trk1, trk2, tmp, tmp, tmp);

        topocuts[0] = kTRUE; topocuts[1] = kFALSE; topocuts[2] = kFALSE;
      }

      else{
        TopologicalCutsNew(trk1, trk2, topocuts[0], topocuts[1], topocuts[2]);
      }

      if (!topocuts[0] && !topocuts[1] && !topocuts[2]) continue;
      // cout << mass << endl;

      TVector3 mTrk1Mom, mTrk2Mom, mResMom;
      mTrk1Mom = trk1->gMom(mVertex, Bfield);
      mTrk2Mom = trk2->gMom(mVertex, Bfield);

      mResMom = mTrk1Mom + mTrk2Mom;

      int charge1 = trk1->charge();
      int charge2 = trk2->charge();

      double toFillDaug[6] = {mCentMaker->GetCent9() + 0.5, mResMom.Perp()*trk1->charge(), mTrk1Mom.Perp(), mass, mTrk2Mom.Perp(), mResMom.PseudoRapidity() };

      for (int cut = 0; cut < 3; cut++){
        if (!topocuts[cut]) continue;
        if (charge1*charge2<0){ //Unlike sign pair
          if (cut == 0)      { fd0 = kTRUE;     fd0TrackIndices.push_back({fPionIndices[piontrack], fKaonIndices[kaontrack], mass});}
          else if (cut == 1) { fTightd0 = kTRUE; fTightd0TrackIndices.push_back({fPionIndices[piontrack], fKaonIndices[kaontrack], mass});}
          else if (cut == 2) { fLoosed0 = kTRUE; fLoosed0TrackIndices.push_back({fPionIndices[piontrack], fKaonIndices[kaontrack], mass});}
          if (fdoDaug_Histograms) hD0CentPtEtaMDphiDaug[cut]->Fill(toFillDaug, mCentMaker->GetWeight());

          if (cut == 0 && fdoSingleParticleEmbedding){

            double pt = mResMom.Perp();
            double phi = mResMom.Phi();
            if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
            if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
            double eta = mResMom.PseudoRapidity();
            double px = mResMom.x();
            double py = mResMom.y();
            double pz = mResMom.z();
            double p = mResMom.Mag();

            // double low_ = phi - pi/2.;
            // double high_ = phi + pi/2.;
            // if(low_ < 0.0)    low_ += 2.0*pi;
            //       if(low_ > 2.0*pi) low_ -= 2.0*pi;
            // if(high_ < 0.0)    high_ += 2.0*pi;
            //       if(high_ > 2.0*pi) high_ -= 2.0*pi;
           
            // double phi_ = 0;
            // if(high_>low_)
            //   phi_ = gRandom->Uniform(low_,high_);
            // else
            //   phi_ = gRandom->Uniform(high_,low_);

            // if(phi_ < 0.0)    phi_ += 2.0*pi;  // force from 0-2pi                                                                                                                                                                                           
            // if(phi_ > 2.0*pi) phi_ -= 2.0*pi;  // force from 0-2pi 

            double phi_ = phi + TMath::Pi();
            if(phi_ < 0.0)    phi_ += 2.0*pi;  // force from 0-2pi                                                                                                                                                                                           
            if(phi_ > 2.0*pi) phi_ -= 2.0*pi;  // force from 0-2pi 

            double eta_ = 0;
            if(eta>0) 
              eta_ = gRandom->Uniform(-1,0);
            else
              eta_ = gRandom->Uniform(0,1);

            double pi0mass = Pico::mMass[0]; // GeV

            double pt_ = gRandom->Uniform(0,40);
            double p_ = sqrt(pt_*pt_ + pt_*pt_* TMath::SinH(eta_)* TMath::SinH(eta_));
            double px_ = pt_*TMath::Cos(phi_);
            double py_ = pt_*TMath::Sin(phi_);
            double pz_ = pt_ * TMath::SinH(eta_);
            double pe_ = sqrt(p_*p_+pi0mass*pi0mass);

            // cout << "EMBEDDED PARTICLE = " << px_ << "\t" << py_ << "\t" << pz_ << "\t" << pe_ << endl;

            fSingleParticle.SetPxPyPzE(px_, py_, pz_, pe_);
            fSingleParticleVector.push_back(fSingleParticle);
          }
        }
        else{ //Like sign pair
          if (cut == 0)      { fd0Bg = kTRUE; fd0BgTrackIndices.push_back({fPionIndices[piontrack], fKaonIndices[kaontrack], mass});}
          else if (cut == 1) { fTightd0Bg = kTRUE; fTightd0BgTrackIndices.push_back({fPionIndices[piontrack], fKaonIndices[kaontrack], mass});}
          else if (cut == 2) { fLoosed0Bg = kTRUE; fLoosed0BgTrackIndices.push_back({fPionIndices[piontrack], fKaonIndices[kaontrack], mass});}
          if (fdoDaug_Histograms) hD0CentPtEtaMDphiDaugLikeSign[cut]->Fill(toFillDaug, mCentMaker->GetWeight());
        }
      }
    }
  }
}



Bool_t StTagKStarEvents::IsAnAcceptableTrack(StPicoTrack *trk, bool dohistograms = kFALSE){ // This is updated as of Dec 25, 2021 to match the D0 Analysis from 2018. The cuts haven't been changed.

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

  if (pt < 0.3) return kFALSE; //This is different
  
  if (dohistograms){ hAcceptedPt[1]->Fill(pt); hAcceptedEta[1]->Fill(eta); }

  if (abs(eta) >= 1) return kFALSE;

  if (dohistograms){ hAcceptedPt[2]->Fill(pt); hAcceptedEta[2]->Fill(eta); }

  if (trk->nHitsFit() < 20) return kFALSE; //This is different

  if (dohistograms){ hAcceptedPt[3]->Fill(pt); hAcceptedEta[3]->Fill(eta); }

  if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 1.5) return kFALSE;

  if (trk->isPrimary()) {
    // return kFALSE;
    if (dohistograms){ hAcceptedPt[4]->Fill(pt); hAcceptedEta[4]->Fill(eta); }

    if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() < 1.) { // This is not there in the analysis note. I am a little confused.

      if (dohistograms){ hAcceptedPt[5]->Fill(pt); hAcceptedEta[5]->Fill(eta); }

      if (float(trk->nHitsFit())/float(trk->nHitsMax()) > 0.52) { // This is not there in the analysis note. I am a little confused.

        if (dohistograms){ hAcceptedPt[6]->Fill(pt); hAcceptedEta[6]->Fill(eta); }

      }
    } 
  }

  // HFT Hits only

  if(!doppAnalysis && fTopoLevel > 0){ // This is the same as the old cuts. Hmm.
    if (!trk->isHFTTrack()) return kFALSE;
  }

  if (dohistograms){ hAcceptedPt[7]->Fill(pt); hAcceptedEta[7]->Fill(eta); }

  if(!doppAnalysis && fTopoLevel > 0){ // This is the same as the old cuts. Hmm.
    if (trk->hasPxl1Hit() && trk->hasPxl2Hit()) {
      if (dohistograms){ hAcceptedPt[8]->Fill(pt); hAcceptedEta[8]->Fill(eta); }
      if (trk->hasIstHit()) {
        if (dohistograms){ hAcceptedPt[9]->Fill(pt); hAcceptedEta[9]->Fill(eta); }
      }
    }
  }

  return kTRUE;
}

void StTagKStarEvents::TopologicalCutsNew(StPicoTrack *trk1, StPicoTrack *trk2, bool &standard, bool &tight, bool &loose){ // Always needs pion as first track, kaon as second track

  standard = kFALSE;
  tight = kFALSE;
  loose = kFALSE;

  int pid1 = trk1->charge();
  int pid2 = trk2->charge();

  TVector3 mTrk1Mom = trk1->gMom(mVertex, Bfield);
  TVector3 mTrk2Mom = trk2->gMom(mVertex, Bfield);

  double mass = InvariantMass(trk1, trk2); //Always needs pion as first track, kaon as second track
  // cout << "Mass : " << mass << endl;

  TVector3 mD0Mom;
  mD0Mom = mTrk1Mom + mTrk2Mom;
  double d0pt = mD0Mom.Perp();

  int centbin = GetFiveCentBin(fCentralityScaled);
  int d0bin = GetD0PtBin(d0pt);

  if (fDebugLevel == 1) cout << "Cent and D0 Pt Bin: " << centbin << "\t" << d0bin << endl;
  if (centbin == -99 || d0bin == -99) return;

  // to be used for testing with preview II pico production
  StPicoPhysicalHelix pHelix = trk1->helix(Bfield);
  StPicoPhysicalHelix kHelix = trk2->helix(Bfield);

  pHelix.moveOrigin(pHelix.pathLength(mVertex));
  kHelix.moveOrigin(kHelix.pathLength(mVertex));
    
  if ( (pHelix.at(pHelix.pathLength(mVertex)) - mVertex).Mag() <= 0.0050 || (kHelix.at(kHelix.pathLength(mVertex)) - mVertex).Mag() <= 0.0050) return;

  StPicoPhysicalHelix const pStraightLine(pHelix.momentum(Bfield*(kilogauss)), pHelix.origin(), 0, trk1->charge());
  StPicoPhysicalHelix const kStraightLine(kHelix.momentum(Bfield*(kilogauss)), kHelix.origin(), 0, trk2->charge());

  pair<double, double> const ss = pStraightLine.pathLengths(kStraightLine);
  TVector3 const pAtDcaToKaon = pStraightLine.at(ss.first);
  TVector3 const kAtDcaToPion = kStraightLine.at(ss.second);

  // calculate DCA of pion to kaon at their DCA
  Double_t mDcaDaughters = (pAtDcaToKaon - kAtDcaToPion).Mag();

  TVector3 const pMomAtDca = pHelix.momentumAt(ss.first,  Bfield*(kilogauss)); //*(1.e-14);
  TVector3 const kMomAtDca = kHelix.momentumAt(ss.second, Bfield*(kilogauss)); //*(1.e-14);

  // cout << "======" << endl;
  // cout << pMomAtDca.X() << "\t" << pMomAtDca.Y() << "\t" << pMomAtDca.Z() << endl;
  // cout << mTrk1Mom.X() << "\t" << mTrk1Mom.Y() << "\t" << mTrk1Mom.Z() << endl;
  // cout << kMomAtDca.X() << "\t" << kMomAtDca.Y() << "\t" << kMomAtDca.Z() << endl;
  // cout << mTrk2Mom.X() << "\t" << mTrk2Mom.Y() << "\t" << mTrk2Mom.Z() << endl;

  TLorentzVector pFourMom;
  TLorentzVector kFourMom;

  pFourMom.SetXYZM(pMomAtDca.X(), pMomAtDca.Y(), pMomAtDca.Z(), M_PION_PLUS);
  kFourMom.SetXYZM(kMomAtDca.X(), kMomAtDca.Y(), kMomAtDca.Z(), M_KAON_PLUS);

  TLorentzVector mLorentzVector = pFourMom + kFourMom;

  if (fDebugLevel == 1) cout << "Mass : " << mass << "\t" << mLorentzVector.M() << endl;
  if (fDebugLevel == 1) cout << "Pt : " << mD0Mom.Perp() << "\t" << mLorentzVector.Vect().Perp() << endl;
  if (fDebugLevel == 1) cout << "Rapidity = " << mLorentzVector.Rapidity() << endl;

  if (abs(mLorentzVector.Rapidity()) >= 1.) return;

  // calculate pointing angle and decay length
  TVector3 const vtxToV0 = (pAtDcaToKaon + kAtDcaToPion) * 0.5 - mVertex;
  Double_t mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());

  Double_t point1 = (pAtDcaToKaon - mVertex).Angle(pFourMom.Vect());
  Double_t point2 = (kAtDcaToPion - mVertex).Angle(kFourMom.Vect());

  // cout << "Pointing Angles : " << mPointingAngle << "\t" << point1 << "\t" << point2 << endl;

  Double_t costheta = TMath::Cos(mPointingAngle);

  Double_t mDecayLength = vtxToV0.Mag();

  Double_t perpDcaToVtx = mDecayLength*TMath::Sin(mPointingAngle);

  // calculate DCA of tracks to primary vertex
  Double_t mPionDca = (pHelix.origin() - mVertex).Mag();
  Double_t mKaonDca = (kHelix.origin() - mVertex).Mag();

  if (fTopoLevel == 0) FillTopoHistograms(centbin, d0bin, pid1*pid2, d0pt, mDecayLength, mDcaDaughters, perpDcaToVtx, mKaonDca, mPionDca, costheta);
  else{
    if (FillStandardTopoHistograms(centbin, d0bin, pid1*pid2, d0pt, mDecayLength, mDcaDaughters, perpDcaToVtx, mKaonDca, mPionDca, costheta)) standard = kTRUE;
    if (FillFiftyPercentSignificanceTopoHistograms(centbin, d0bin, pid1*pid2, d0pt, mDecayLength, mDcaDaughters, perpDcaToVtx, mKaonDca, mPionDca, costheta)) tight = kTRUE;
    if (FillOneHundredFiftyPercentSignificanceTopoHistograms(centbin, d0bin, pid1*pid2, d0pt, mDecayLength, mDcaDaughters, perpDcaToVtx, mKaonDca, mPionDca, costheta)) loose = kTRUE;
  }

  return;
}

void StTagKStarEvents::FillTopoHistograms(int centbin, int d0bin, int pid1pid2, double d0pt, double mDecayLength, double mDcaDaughters, double perpDcaToVtx, double mKaonDca, double mPionDca, double costheta){
  
  if (fdoDaug_Histograms){
    if (pid1pid2<0) decaylengthd0US[0]->Fill(centbin, d0pt, mDecayLength);
    if (pid1pid2>0) decaylengthd0LS[0]->Fill(centbin, d0pt, mDecayLength);

    if (pid1pid2<0) distancepikUS[0]->Fill(centbin, d0pt, mDcaDaughters);
    if (pid1pid2>0) distancepikLS[0]->Fill(centbin, d0pt, mDcaDaughters);

    if (pid1pid2<0) distanced0PVUS[0]->Fill(centbin, d0pt, perpDcaToVtx);
    if (pid1pid2>0) distanced0PVLS[0]->Fill(centbin, d0pt, perpDcaToVtx);

    if (pid1pid2<0) dcakPVUS[0]->Fill(centbin, d0pt, mKaonDca);
    if (pid1pid2>0) dcakPVLS[0]->Fill(centbin, d0pt, mKaonDca);

    if (pid1pid2<0) dcapiPVUS[0]->Fill(centbin, d0pt, mPionDca);
    if (pid1pid2>0) dcapiPVLS[0]->Fill(centbin, d0pt, mPionDca);

    if (pid1pid2<0) costhetaDVPVUS[0]->Fill(centbin, d0pt, costheta);
    if (pid1pid2>0) costhetaDVPVLS[0]->Fill(centbin, d0pt, costheta);
  }
}


bool StTagKStarEvents::FillStandardTopoHistograms(int centbin, int d0bin, int pid1pid2, double d0pt, double mDecayLength, double mDcaDaughters, double perpDcaToVtx, double mKaonDca, double mPionDca, double costheta){
  
  bool DoesItPassTheCut = kFALSE;

  const int centbins = 5;
  const int d0bins = 6;
  // default
  // float const cosTheta = 0.95;
  float const kDca[centbins][d0bins] = {
      {0.0106, 0.0106, 0.0069, 0.0068, 0.0050, 0.0050}, // 60-80%
      {0.0140, 0.0100, 0.0075, 0.0072, 0.0060, 0.0050}, // 40-60%
      {0.0151, 0.0102, 0.0104, 0.0099, 0.0063, 0.0050}, // 20-40%
      {0.0145, 0.0113, 0.0094, 0.0089, 0.0069, 0.0050}, // 10-20%
      {0.0138, 0.0109, 0.0082, 0.0094, 0.0076, 0.0054}  // 0-10%
  };
  float const pDca[centbins][d0bins] = {
      {0.0098, 0.0098, 0.0083, 0.0073, 0.0056, 0.0050}, // 60-80%
      {0.0145, 0.0128, 0.0072, 0.0079, 0.0060, 0.0051}, // 40-60%
      {0.0131, 0.0113, 0.0099, 0.0106, 0.0065, 0.0052}, // 20-40%
      {0.0141, 0.0100, 0.0074, 0.0077, 0.0066, 0.0052}, // 10-20%
      {0.0133, 0.0105, 0.0093, 0.0097, 0.0067, 0.0055}  // 0-10%
  };
  float const dcaV0ToPv[centbins][d0bins] = {
      {0.0076, 0.0076, 0.0053, 0.0054, 0.0054, 0.0042}, // 60-80%
      {0.0072, 0.0057, 0.0058, 0.0049, 0.0049, 0.0047}, // 40-60%
      {0.0066, 0.0055, 0.0053, 0.0046, 0.0041, 0.0050}, // 20-40%
      {0.0063, 0.0047, 0.0045, 0.0046, 0.0042, 0.0044}, // 10-20%
      {0.0062, 0.0055, 0.0040, 0.0040, 0.0040, 0.0044}  // 0-10%
  };
  float const dcaDaughters[centbins][d0bins] = {
      {0.0077, 0.0077, 0.0094, 0.0078, 0.0081, 0.0120}, // 60-80%
      {0.0080, 0.0083, 0.0092, 0.0081, 0.0094, 0.0106}, // 40-60%
      {0.0078, 0.0073, 0.0080, 0.0093, 0.0096, 0.0103}, // 20-40%
      {0.0076, 0.0078, 0.0092, 0.0072, 0.0086, 0.0085}, // 10-20%
      {0.0071, 0.0064, 0.0070, 0.0063, 0.0082, 0.0080}  // 0-10%
  };
  float const decayLength[centbins][d0bins] = {
      {0.0175, 0.0175, 0.0187, 0.0178, 0.0184, 0.0187}, // 60-80%
      {0.0171, 0.0196, 0.0210, 0.0187, 0.0190, 0.0214}, // 40-60%
      {0.0178, 0.0206, 0.0221, 0.0209, 0.0219, 0.0240}, // 20-40%
      {0.0172, 0.0215, 0.0252, 0.0232, 0.0236, 0.0237}, // 10-20%
      {0.0100, 0.0199, 0.0227, 0.0232, 0.0236, 0.0255}  // 0-10%
  };

  if (fdoDaug_Histograms){
    if ((mDcaDaughters < dcaDaughters[centbin][d0bin])
     // && (mDecayLength > decayLength[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) decaylengthd0US[0]->Fill(centbin, d0pt, mDecayLength);
      if (pid1pid2>0) decaylengthd0LS[0]->Fill(centbin, d0pt, mDecayLength);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     // && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) distancepikUS[0]->Fill(centbin, d0pt, mDcaDaughters);
      if (pid1pid2>0) distancepikLS[0]->Fill(centbin, d0pt, mDcaDaughters);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     // && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) distanced0PVUS[0]->Fill(centbin, d0pt, perpDcaToVtx);
      if (pid1pid2>0) distanced0PVLS[0]->Fill(centbin, d0pt, perpDcaToVtx);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     // && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) dcakPVUS[0]->Fill(centbin, d0pt, mKaonDca);
      if (pid1pid2>0) dcakPVLS[0]->Fill(centbin, d0pt, mKaonDca);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     // && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) dcapiPVUS[0]->Fill(centbin, d0pt, mPionDca);
      if (pid1pid2>0) dcapiPVLS[0]->Fill(centbin, d0pt, mPionDca);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     // && costheta > 0.95
      )
    {
      if (pid1pid2<0) costhetaDVPVUS[0]->Fill(centbin, d0pt, costheta);
      if (pid1pid2>0) costhetaDVPVLS[0]->Fill(centbin, d0pt, costheta);
    }                                                               
  }

  if ((mDecayLength > decayLength[centbin][d0bin])
   && (mDcaDaughters < dcaDaughters[centbin][d0bin])
   && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
   && (mKaonDca > kDca[centbin][d0bin]) 
   && (mPionDca > pDca[centbin][d0bin]) 
   && (costheta > 0.95)){
    DoesItPassTheCut = kTRUE;
  }


  return DoesItPassTheCut;

}

bool StTagKStarEvents::FillFiftyPercentSignificanceTopoHistograms(int centbin, int d0bin, int pid1pid2, double d0pt, double mDecayLength, double mDcaDaughters, double perpDcaToVtx, double mKaonDca, double mPionDca, double costheta){
  
  bool DoesItPassTheCut = kFALSE;

  const int centbins = 5;
  const int d0bins = 6;
  // tight
  // float const cosTheta = 0.95;
  float const kDca[centbins][d0bins] = {
      {0.0126, 0.0126, 0.0116, 0.0097, 0.0076, 0.0050}, // 60-80%
      {0.0176, 0.0100, 0.0123, 0.0092, 0.0068, 0.0056}, // 40-60%
      {0.0158, 0.0123, 0.0133, 0.0125, 0.0103, 0.0053}, // 20-40%
      {0.0172, 0.0165, 0.0119, 0.0125, 0.0091, 0.0057}, // 10-20%
      {0.0170, 0.0109, 0.0117, 0.0109, 0.0111, 0.0059}  // 0-10%
  };
  float const pDca[centbins][d0bins] = {
      {0.0130, 0.0130, 0.0130, 0.0095, 0.0097, 0.0086}, // 60-80%
      {0.0143, 0.0100, 0.0072, 0.0145, 0.0113, 0.0095}, // 40-60%
      {0.0150, 0.0113, 0.0063, 0.0149, 0.0107, 0.0105}, // 20-40%
      {0.0172, 0.0100, 0.0144, 0.0133, 0.0130, 0.0092}, // 10-20%
      {0.0139, 0.0148, 0.0093, 0.0133, 0.0080, 0.0088}  // 0-10%
  };
  float const dcaV0ToPv[centbins][d0bins] = {
      {0.0058, 0.0058, 0.0051, 0.0036, 0.0036, 0.0031}, // 60-80%
      {0.0056, 0.0045, 0.0050, 0.0037, 0.0029, 0.0026}, // 40-60%
      {0.0062, 0.0060, 0.0037, 0.0038, 0.0030, 0.0026}, // 20-40%
      {0.0056, 0.0049, 0.0039, 0.0036, 0.0033, 0.0028}, // 10-20%
      {0.0055, 0.0053, 0.0037, 0.0027, 0.0026, 0.0025}  // 0-10%
  };
  float const dcaDaughters[centbins][d0bins] = {
      {0.0076, 0.0076, 0.0078, 0.0095, 0.0076, 0.0093}, // 60-80%
      {0.0088, 0.0059, 0.0081, 0.0083, 0.0073, 0.0108}, // 40-60%
      {0.0061, 0.0043, 0.0070, 0.0089, 0.0066, 0.0070}, // 20-40%
      {0.0077, 0.0049, 0.0042, 0.0056, 0.0053, 0.0119}, // 10-20%
      {0.0077, 0.0044, 0.0047, 0.0073, 0.0060, 0.0061}  // 0-10%
  };
  float const decayLength[centbins][d0bins] = {
      {0.0203, 0.0203, 0.0206, 0.0228, 0.0161, 0.0216}, // 60-80%
      {0.0222, 0.0229, 0.0269, 0.0236, 0.0232, 0.0182}, // 40-60%
      {0.0240, 0.0242, 0.0268, 0.0319, 0.0176, 0.0338}, // 20-40%
      {0.0219, 0.0240, 0.0213, 0.0231, 0.0261, 0.0399}, // 10-20%
      {0.0100, 0.0230, 0.0268, 0.0292, 0.0249, 0.0225}  // 0-10%
  };

  if (fdoDaug_Histograms){
    if ((mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) decaylengthd0US[1]->Fill(centbin, d0pt, mDecayLength);
      if (pid1pid2>0) decaylengthd0LS[1]->Fill(centbin, d0pt, mDecayLength);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     // && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) distancepikUS[1]->Fill(centbin, d0pt, mDcaDaughters);
      if (pid1pid2>0) distancepikLS[1]->Fill(centbin, d0pt, mDcaDaughters);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     // && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) distanced0PVUS[1]->Fill(centbin, d0pt, perpDcaToVtx);
      if (pid1pid2>0) distanced0PVLS[1]->Fill(centbin, d0pt, perpDcaToVtx);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     // && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) dcakPVUS[1]->Fill(centbin, d0pt, mKaonDca);
      if (pid1pid2>0) dcakPVLS[1]->Fill(centbin, d0pt, mKaonDca);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     // && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) dcapiPVUS[1]->Fill(centbin, d0pt, mPionDca);
      if (pid1pid2>0) dcapiPVLS[1]->Fill(centbin, d0pt, mPionDca);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     // && costheta > 0.95
      )
    {
      if (pid1pid2<0) costhetaDVPVUS[1]->Fill(centbin, d0pt, costheta);
      if (pid1pid2>0) costhetaDVPVLS[1]->Fill(centbin, d0pt, costheta);
    }  
  }                                                             

  if ((mDecayLength > decayLength[centbin][d0bin])
   && (mDcaDaughters < dcaDaughters[centbin][d0bin])
   && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
   && (mKaonDca > kDca[centbin][d0bin]) 
   && (mPionDca > pDca[centbin][d0bin]) 
   && (costheta > 0.95)){
    DoesItPassTheCut = kTRUE;
  }

  return DoesItPassTheCut;

}

bool StTagKStarEvents::FillOneHundredFiftyPercentSignificanceTopoHistograms(int centbin, int d0bin, int pid1pid2, double d0pt, double mDecayLength, double mDcaDaughters, double perpDcaToVtx, double mKaonDca, double mPionDca, double costheta){
  
  bool DoesItPassTheCut = kFALSE;

  const int centbins = 5;
  const int d0bins = 6;
  // loose
  // float const cosTheta = 0.95;
  float const kDca[centbins][d0bins] = {
      {0.0090, 0.0090, 0.0073, 0.0057, 0.0050, 0.0050}, // 60-80%
      {0.0115, 0.0114, 0.0067, 0.0050, 0.0050, 0.0050}, // 40-60%
      {0.0134, 0.0112, 0.0074, 0.0063, 0.0050, 0.0050}, // 20-40%
      {0.0135, 0.0120, 0.0070, 0.0060, 0.0050, 0.0050}, // 10-20%
      {0.0111, 0.0109, 0.0080, 0.0067, 0.0050, 0.0050}  // 0-10%
  };
  float const pDca[centbins][d0bins] = {
      {0.0098, 0.0098, 0.0080, 0.0058, 0.0050, 0.0050}, // 60-80%
      {0.0124, 0.0114, 0.0072, 0.0050, 0.0050, 0.0050}, // 40-60%
      {0.0111, 0.0113, 0.0089, 0.0062, 0.0050, 0.0050}, // 20-40%
      {0.0108, 0.0100, 0.0074, 0.0069, 0.0050, 0.0050}, // 10-20%
      {0.0116, 0.0105, 0.0093, 0.0072, 0.0050, 0.0050}  // 0-10%
  };
  float const dcaV0ToPv[centbins][d0bins] = {
      {0.0083, 0.0083, 0.0090, 0.0077, 0.0059, 0.0055}, // 60-80%
      {0.0069, 0.0083, 0.0087, 0.0077, 0.0085, 0.0072}, // 40-60%
      {0.0075, 0.0073, 0.0056, 0.0053, 0.0100, 0.0090}, // 20-40%
      {0.0071, 0.0063, 0.0055, 0.0052, 0.0065, 0.0055}, // 10-20%
      {0.0078, 0.0063, 0.0054, 0.0045, 0.0062, 0.0047}  // 0-10%
  };
  float const dcaDaughters[centbins][d0bins] = {
      {0.0097, 0.0097, 0.0123, 0.0092, 0.0093, 0.0092}, // 60-80%
      {0.0095, 0.0084, 0.0100, 0.0108, 0.0130, 0.0113}, // 40-60%
      {0.0081, 0.0090, 0.0078, 0.0094, 0.0137, 0.0121}, // 20-40%
      {0.0082, 0.0086, 0.0099, 0.0102, 0.0121, 0.0105}, // 10-20%
      {0.0098, 0.0087, 0.0088, 0.0078, 0.0100, 0.0101}  // 0-10%
  };
  float const decayLength[centbins][d0bins] = {
      {0.0154, 0.0154, 0.0163, 0.0147, 0.0126, 0.0140}, // 60-80%
      {0.0158, 0.0153, 0.0172, 0.0150, 0.0126, 0.0164}, // 40-60%
      {0.0100, 0.0177, 0.0177, 0.0194, 0.0131, 0.0150}, // 20-40%
      {0.0157, 0.0197, 0.0222, 0.0180, 0.0155, 0.0189}, // 10-20%
      {0.0148, 0.0179, 0.0215, 0.0198, 0.0167, 0.0206}  // 0-10%
  };

  if (fdoDaug_Histograms){

    if ((mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) decaylengthd0US[2]->Fill(centbin, d0pt, mDecayLength);
      if (pid1pid2>0) decaylengthd0LS[2]->Fill(centbin, d0pt, mDecayLength);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     // && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) distancepikUS[2]->Fill(centbin, d0pt, mDcaDaughters);
      if (pid1pid2>0) distancepikLS[2]->Fill(centbin, d0pt, mDcaDaughters);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     // && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) distanced0PVUS[2]->Fill(centbin, d0pt, perpDcaToVtx);
      if (pid1pid2>0) distanced0PVLS[2]->Fill(centbin, d0pt, perpDcaToVtx);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     // && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) dcakPVUS[2]->Fill(centbin, d0pt, mKaonDca);
      if (pid1pid2>0) dcakPVLS[2]->Fill(centbin, d0pt, mKaonDca);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     // && (mPionDca > pDca[centbin][d0bin])
     && costheta > 0.95
      )
    {
      if (pid1pid2<0) dcapiPVUS[2]->Fill(centbin, d0pt, mPionDca);
      if (pid1pid2>0) dcapiPVLS[2]->Fill(centbin, d0pt, mPionDca);
    }

    if ((mDecayLength > decayLength[centbin][d0bin]) 
     && (mDcaDaughters < dcaDaughters[centbin][d0bin]) 
     && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
     && (mKaonDca > kDca[centbin][d0bin]) 
     && (mPionDca > pDca[centbin][d0bin])
     // && costheta > 0.95
      )
    {
      if (pid1pid2<0) costhetaDVPVUS[2]->Fill(centbin, d0pt, costheta);
      if (pid1pid2>0) costhetaDVPVLS[2]->Fill(centbin, d0pt, costheta);
    }
  }                                                               

  if ((mDecayLength > decayLength[centbin][d0bin])
   && (mDcaDaughters < dcaDaughters[centbin][d0bin])
   && (perpDcaToVtx < dcaV0ToPv[centbin][d0bin]) 
   && (mKaonDca > kDca[centbin][d0bin]) 
   && (mPionDca > pDca[centbin][d0bin]) 
   && (costheta > 0.95)){
    DoesItPassTheCut = kTRUE;
  }

  return DoesItPassTheCut;

}

void StTagKStarEvents::FillPidHistograms(StPicoTrack *trk){

  TVector3 mTrkMom;
  mTrkMom = trk->gMom(mVertex, Bfield);

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
      z_pi->Fill(p, zpi);
      double zka = trk->nSigmaKaon();
      z_ka->Fill(p, zka);
  }
}

bool StTagKStarEvents::IsTpcPion(StPicoTrack *trk){
  double zpi = trk->nSigmaPion();
  if (abs(zpi) < 3.) return true;
  return false;
}

bool StTagKStarEvents::IsTpcKaon(StPicoTrack *trk){
  double zka = trk->nSigmaKaon();
  if (abs(zka) < 2.) return true;
  return false;
}

bool StTagKStarEvents::IsPion(StPicoTrack *trk){
  if (!IsTpcPion(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_PION_PLUS*M_PION_PLUS / p / p + 1);
  double nsigma = abs(1./beta - oneOverBetaExpected);
  if (nsigma < 0.03) return true;
  return false;
}

bool StTagKStarEvents::IsKaon(StPicoTrack *trk){
  if (!IsTpcKaon(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_KAON_PLUS*M_KAON_PLUS / p / p + 1);
  double nsigma = abs(1./beta - oneOverBetaExpected);
  if (nsigma < 0.03) return true;
  return false;
}

Double_t StTagKStarEvents::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2){

  TVector3 mTrk1Mom, mTrk2Mom;

  mTrk1Mom = trk1->gMom(mVertex, Bfield);

  mTrk2Mom = trk2->gMom(mVertex, Bfield);

  TLorentzVector mTrk1, mTrk2;

  mTrk1.SetXYZM(mTrk1Mom.X(), mTrk1Mom.Y(), mTrk1Mom.Z(), M_PION_PLUS);
  mTrk2.SetXYZM(mTrk2Mom.X(), mTrk2Mom.Y(), mTrk2Mom.Z(), M_KAON_PLUS);

  TLorentzVector mD0;
  mD0 = mTrk1 + mTrk2;

  return mD0.M();
}

float StTagKStarEvents::GetTofBeta(StPicoTrack *trk){
  int index2tof = trk->bTofPidTraitsIndex();
    
  float beta = std::numeric_limits<double>::quiet_NaN();

  if (index2tof < 0) return beta;

  StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(index2tof);
      
  if (tofPid)
  {
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
Int_t StTagKStarEvents::GetFiveCentBin(Double_t scaledCent) const{
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
Int_t StTagKStarEvents::GetD0PtBin(Double_t D0Pt) const{
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

