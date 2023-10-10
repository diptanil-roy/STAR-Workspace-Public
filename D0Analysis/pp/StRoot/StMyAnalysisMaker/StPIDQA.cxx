// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StPIDQA.h"
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

ClassImp(StPIDQA)

//________________________________________________________________________
StPIDQA::StPIDQA(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StPIDQA::fRunFlagEnum
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
  
  
  fMBEvents = kTRUE;
  fHTEvents = kFALSE;
  //Mass cut
  fInvMassSignal1 = 1.70;
  fInvMassSignal2 = 2.10;

  fInvMassULBg1 = 1.70;
  fInvMassULBg2 = 2.10;

  fInvMassLSBg1 = 1.70;
  fInvMassLSBg2 = 2.10;

  fTrackTree = {0};
  fD0Tree    = {0};

  fD0MassWindow = kTRUE;

  fSaveD0Trees = kFALSE;

  
}

//
//________________________________________________________________________
StPIDQA::~StPIDQA()
{ 

  if (fdoQA_Histograms){

    if (hEventZvertex_diff) delete hEventZvertex_diff;
    if (hEventVzvsVzvpd) delete hEventVzvsVzvpd;
    if (hEventVxvsVy) delete hEventVxvsVy;

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

    if (hD0Mass) delete hD0Mass;
  }

  if(mEmcPosition) delete mEmcPosition;

}

void StPIDQA::DeclareTree() {
  TString ustreename = "D0USTree";
  d0unlikesigntree = new TTree(ustreename.Data(), ustreename.Data());

  TString lstreename = "D0LSTree";
  d0likesigntree = new TTree(lstreename.Data(), lstreename.Data());
}

void StPIDQA::WriteTree() {
  d0unlikesigntree->Write();
  d0likesigntree->Write();
}

void StPIDQA::BookTree()
{
  d0unlikesigntree->Branch("Vz", &fD0Tree.Vz, "Vz/F");
  d0unlikesigntree->Branch("VzVPD", &fD0Tree.VzVPD, "VzVPD/F");
  d0unlikesigntree->Branch("D0Mass", &fD0Tree.D0mass, "D0mass/F");
  d0unlikesigntree->Branch("D0Pt", &fD0Tree.D0pt, "D0pt/F");
  d0unlikesigntree->Branch("D0Eta", &fD0Tree.D0eta, "D0eta/F");
  d0unlikesigntree->Branch("D0Phi", &fD0Tree.D0phi, "D0phi/F");
  d0unlikesigntree->Branch("PionPt", &fD0Tree.pionpt, "pionpt/F");
  d0unlikesigntree->Branch("PionEta", &fD0Tree.pioneta, "pioneta/F");
  d0unlikesigntree->Branch("PionPhi", &fD0Tree.pionphi, "pionphi/F");
  d0unlikesigntree->Branch("PionDCA", &fD0Tree.piondca, "piondca/F");
  d0unlikesigntree->Branch("PionNHitsFit", &fD0Tree.pionnhitsfit, "pionnhitsfit/F");
  d0unlikesigntree->Branch("PionNSigmaPion", &fD0Tree.pionnsigmapion, "pionnsigmapion/F");
  d0unlikesigntree->Branch("PionNSigmaKaon", &fD0Tree.pionnsigmakaon, "pionnsigmakaon/F");
  d0unlikesigntree->Branch("PionTofBeta", &fD0Tree.piontofbeta, "piontofbeta/F");
  d0unlikesigntree->Branch("PionTofYLocal", &fD0Tree.ntofpionylocal, "ntofpionylocal/F");
  d0unlikesigntree->Branch("KaonPt", &fD0Tree.kaonpt, "kaonpt/F");
  d0unlikesigntree->Branch("KaonEta", &fD0Tree.kaoneta, "kaoneta/F");
  d0unlikesigntree->Branch("KaonPhi", &fD0Tree.kaonphi, "kaonphi/F");
  d0unlikesigntree->Branch("KaonDCA", &fD0Tree.kaondca, "kaondca/F");
  d0unlikesigntree->Branch("KaonNHitsFit", &fD0Tree.kaonnhitsfit, "kaonnhitsfit/F");
  d0unlikesigntree->Branch("KaonNSigmaPion", &fD0Tree.kaonnsigmapion, "kaonnsigmapion/F");
  d0unlikesigntree->Branch("KaonNSigmaKaon", &fD0Tree.kaonnsigmakaon, "kaonnsigmakaon/F");
  d0unlikesigntree->Branch("KaonTofBeta", &fD0Tree.kaontofbeta, "kaontofbeta/F");
  d0unlikesigntree->Branch("KaonTofYLocal", &fD0Tree.ntofkaonylocal, "ntofkaonylocal/F");

  d0likesigntree->Branch("Vz", &fD0Tree.Vz, "Vz/F");
  d0likesigntree->Branch("VzVPD", &fD0Tree.VzVPD, "VzVPD/F");
  d0likesigntree->Branch("D0Mass", &fD0Tree.D0mass, "D0mass/F");
  d0likesigntree->Branch("D0Pt", &fD0Tree.D0pt, "D0pt/F");
  d0likesigntree->Branch("D0Eta", &fD0Tree.D0eta, "D0eta/F");
  d0likesigntree->Branch("D0Phi", &fD0Tree.D0phi, "D0phi/F");
  d0likesigntree->Branch("PionPt", &fD0Tree.pionpt, "pionpt/F");
  d0likesigntree->Branch("PionEta", &fD0Tree.pioneta, "pioneta/F");
  d0likesigntree->Branch("PionPhi", &fD0Tree.pionphi, "pionphi/F");
  d0likesigntree->Branch("PionDCA", &fD0Tree.piondca, "piondca/F");
  d0likesigntree->Branch("PionNHitsFit", &fD0Tree.pionnhitsfit, "pionnhitsfit/F");
  d0likesigntree->Branch("PionNSigmaPion", &fD0Tree.pionnsigmapion, "pionnsigmapion/F");
  d0likesigntree->Branch("PionNSigmaKaon", &fD0Tree.pionnsigmakaon, "pionnsigmakaon/F");
  d0likesigntree->Branch("PionTofBeta", &fD0Tree.piontofbeta, "piontofbeta/F");
  d0likesigntree->Branch("PionTofYLocal", &fD0Tree.ntofpionylocal, "ntofpionylocal/F");
  d0likesigntree->Branch("KaonPt", &fD0Tree.kaonpt, "kaonpt/F");
  d0likesigntree->Branch("KaonEta", &fD0Tree.kaoneta, "kaoneta/F");
  d0likesigntree->Branch("KaonPhi", &fD0Tree.kaonphi, "kaonphi/F");
  d0likesigntree->Branch("KaonDCA", &fD0Tree.kaondca, "kaondca/F");
  d0likesigntree->Branch("KaonNHitsFit", &fD0Tree.kaonnhitsfit, "kaonnhitsfit/F");
  d0likesigntree->Branch("KaonNSigmaPion", &fD0Tree.kaonnsigmapion, "kaonnsigmapion/F");
  d0likesigntree->Branch("KaonNSigmaKaon", &fD0Tree.kaonnsigmakaon, "kaonnsigmakaon/F");
  d0likesigntree->Branch("KaonTofBeta", &fD0Tree.kaontofbeta, "kaontofbeta/F");
  d0likesigntree->Branch("KaonTofYLocal", &fD0Tree.ntofkaonylocal, "ntofkaonylocal/F");
}

//________________________________________________________________________
Int_t StPIDQA::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  if (fSaveD0Trees){
    DeclareTree();
    BookTree();
  }

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StPIDQA::Finish() { 

  cout << "StPIDQA::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),  "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();

    if (fSaveD0Trees){
      fout->cd();
      fout->mkdir(Form("%s_%s", GetName(), "TrackTree"));
      fout->cd(Form("%s_%s", GetName(), "TrackTree"));
      WriteTree();
    }

    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StPIDQA::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StPIDQA::DeclareHistograms() {

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
    hEventVzvsVzvpd = new TH2F("hEventVzvsVzvpd", "Vz vs VPD Vz", 2000, -200, 200, 2000, -200, 200);
    hEventVxvsVy = new TH2F("hEventVxvsVy", "Vx vs Vy", 2000, -200, 200, 2000, -200, 200);

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

    hD0Mass   = new TH1F("hD0Mass", "hD0Mass", 1000, 0, 2.5);
  }
}

//
// write histograms
//_____________________________________________________________________________
void StPIDQA::WriteHistograms() {

  if (fdoQA_Histograms){
    hEventZvertex_diff->Write();
    hEventVzvsVzvpd->Write();
    hEventVxvsVy->Write();

    hRefMultiplicity->Write();
    hgRefMultiplicity->Write();

    // Event Cuts Histograms
    const char *event_cuts[14] = {"Total Events", "Good Runs", "Max Track Pt < 30", "|V_{i}| != 0", "|V_{z}| -> No Cut", "VPDMB/HT","|V_{z}| < 100 cm", "|V_{z} - V_{z(VPD)}| < 6 cm", "Ranking > 0", "nBEMCMatch && nBTOFMatch > 0", "", "", "", ""};
    const char *track_cuts[10] = {"Total Tracks", "nHitsDedx > 15", "#frac{nHitsDedx}{nHitsMax} > 0.52", "DCA < 3 cm", "BEMC Hit", "BTOF Hit", "Charged Tracks", "", "", ""};

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

    hD0Mass->Write();
  }
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StPIDQA::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StPIDQA::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  fd0TrackIndices.clear();

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

  if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5) return kStOK;

  numberofevents[3]++;
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  // if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

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

    int arrMB_Run12[]   = {370001, 370011};
    int arrHT_Run12[]   = {370501, 370541, 370542, 370511, 370546, 370521, 370522, 370531, 370980};
    int arrJP_Run12[]   = {370601, 370611, 370621, 370641};
    
    bool matchMB = kFALSE;
    bool matchHT = kFALSE;

    for(int i = 0; i < sizeof(arrMB_Run12)/sizeof(*arrMB_Run12); i++) {
      if(mPicoEvent->isTrigger(arrMB_Run12[i])) matchMB = kTRUE;
      if(matchMB) break;
    }

    for(int i = 0; i < sizeof(arrHT_Run12)/sizeof(*arrHT_Run12); i++) {
      if(mPicoEvent->isTrigger(arrHT_Run12[i])) matchHT = kTRUE;
      if(matchHT) break;
    }

    for(int i = 0; i < sizeof(arrJP_Run12)/sizeof(*arrJP_Run12); i++) {
      if(mPicoEvent->isTrigger(arrJP_Run12[i])) matchHT = kTRUE;
      if(matchHT) break;
    }
   

    if (fMBEvents && !matchMB) return kStOk;
    if (fHTEvents && !matchHT) return kStOk;

    numberofevents[5]++;
    // ======================== end of Triggers ============================= //

    if (abs(zVtx) > 100.) return kStOk;
    numberofevents[6]++;

    if (mPicoEvent->ranking() <= 0) return kStOk;
    numberofevents[7]++;

    if (fdoQA_Histograms){
      hEventZvertex_diff->Fill(zVtx - zVtx_VPD);
      hEventVzvsVzvpd->Fill(zVtx, zVtx_VPD);
      hEventVxvsVy->Fill(mVertex.x(), mVertex.y());
    }

    if (abs(zVtx - zVtx_VPD) > 6) return kStOk;
    numberofevents[8]++;

     
    if (mPicoEvent->nBEMCMatch() != 0 || mPicoEvent->nBTOFMatch() != 0) numberofevents[9]++;

    // fill histograms
    if (fdoQA_Histograms){
      hRefMultiplicity->Fill(refMult);
      hgRefMultiplicity->Fill(grefMult);
    }
  }

  MakeD0s();

  return kStOK;
}

void StPIDQA::MakeD0s(){
  const Int_t ntracks = mPicoDst->numberOfTracks();

  fPionIndices.clear();
  fKaonIndices.clear();

  for(unsigned short itrk = 0; itrk < ntracks; itrk++){
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
    if(!trk){ continue; }

    if(!IsAnAcceptableTrack(trk, kTRUE)) continue;

    FillPidHistograms(trk);

    if (IsPion(trk)) fPionIndices.push_back(itrk);
    if (IsKaon(trk)) fKaonIndices.push_back(itrk);
  } 

  for (int pion = 0; pion < fPionIndices.size(); pion++){
    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(fPionIndices[pion]));
    if(!trk1){ continue; }

    for (int kaon = 0; kaon < fKaonIndices.size(); kaon++){
      StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(fKaonIndices[kaon]));
      if(!trk2){ continue; }

      if (fPionIndices[pion] == fKaonIndices[kaon]) continue; //Don't want a D0 candidate made out of the same track as a pion or kaon
      double mass = InvariantMass(trk1, trk2);

      TVector3 mTrk1Mom, mTrk2Mom, mResMom;
      mTrk1Mom = trk1->gMom(mVertex, Bfield);
      mTrk2Mom = trk2->gMom(mVertex, Bfield);

      mResMom = mTrk1Mom + mTrk2Mom;

      int charge1 = trk1->charge();
      int charge2 = trk2->charge();

      fD0Tree = {0};

      fD0Tree.Vz = mVertex.z();
      fD0Tree.VzVPD = mPicoEvent->vzVpd();
      fD0Tree.D0mass = mass;
      fD0Tree.D0pt = mResMom.Perp();
      fD0Tree.D0eta = mResMom.PseudoRapidity();
      double phi = standardPhi(mResMom.Phi());
      if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
      if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
      fD0Tree.D0phi = phi;

      fD0Tree.pionpt = mTrk1Mom.Perp();
      fD0Tree.pioneta = mTrk1Mom.PseudoRapidity();
      fD0Tree.pionphi = standardPhi(mTrk1Mom.Phi());
      fD0Tree.piondca = trk1->gDCA(mPicoEvent->primaryVertex()).Mag();
      fD0Tree.pionnhitsfit = trk1->nHitsFit();
      fD0Tree.pionnsigmapion = trk1->nSigmaPion();
      fD0Tree.pionnsigmakaon = trk1->nSigmaKaon();

      float pionbeta = GetTofBeta(trk1);
      fD0Tree.piontofbeta = (isnan(pionbeta) || pionbeta < 0) ? -999 : pionbeta;
      fD0Tree.ntofpionylocal = (isnan(pionbeta) || pionbeta < 0) ? -999 : mPicoDst->btofPidTraits(trk1->bTofPidTraitsIndex())->btofYLocal();

      fD0Tree.kaonpt = mTrk2Mom.Perp();
      fD0Tree.kaoneta = mTrk2Mom.PseudoRapidity();
      fD0Tree.kaonphi = standardPhi(mTrk2Mom.Phi());
      fD0Tree.kaondca = trk2->gDCA(mPicoEvent->primaryVertex()).Mag();
      fD0Tree.kaonnhitsfit = trk2->nHitsFit();
      fD0Tree.kaonnsigmapion = trk2->nSigmaPion();
      fD0Tree.kaonnsigmakaon = trk2->nSigmaKaon();

      float kaonbeta = GetTofBeta(trk2);
      fD0Tree.kaontofbeta = (isnan(kaonbeta) || kaonbeta < 0) ? -999 : kaonbeta;
      fD0Tree.ntofkaonylocal = (isnan(kaonbeta) || kaonbeta < 0) ? -999 : mPicoDst->btofPidTraits(trk2->bTofPidTraitsIndex())->btofYLocal();

      if (charge1*charge2 < 0) hD0Mass->Fill(mass);

      if (fD0MassWindow && (mass < 1.7 || mass > 2.1)) continue;

      if (fSaveD0Trees){
        if (charge1*charge2 < 0) d0unlikesigntree->Fill();
        else d0likesigntree->Fill();
      }

      fd0TrackIndices.push_back({fPionIndices[pion], fKaonIndices[kaon]});
    }
  }
}

Bool_t StPIDQA::IsAnAcceptableTrack(StPicoTrack *trk, bool dohistograms = kFALSE){ // This is updated as of Dec 25, 2021 to match the D0 Analysis from 2018. The cuts haven't been changed.

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

  if (pt < 0.2) return kFALSE; //This is different
  
  if (dohistograms){ hAcceptedPt[1]->Fill(pt); hAcceptedEta[1]->Fill(eta); }

  if (abs(eta) > 1) return kFALSE;

  if (dohistograms){ hAcceptedPt[2]->Fill(pt); hAcceptedEta[2]->Fill(eta); }

  if (trk->nHitsFit() < 15) return kFALSE; //This is different

  if (dohistograms){ hAcceptedPt[3]->Fill(pt); hAcceptedEta[3]->Fill(eta); }

  if (float(trk->nHitsFit())/float(trk->nHitsMax()) < 0.52) return kFALSE;
 
  if (dohistograms){ hAcceptedPt[4]->Fill(pt); hAcceptedEta[4]->Fill(eta); }

  if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 2.) return kFALSE;

  if (dohistograms){ hAcceptedPt[5]->Fill(pt); hAcceptedEta[5]->Fill(eta); }

  if (trk->isPrimary()) {
    if (dohistograms){ hAcceptedPt[6]->Fill(pt); hAcceptedEta[6]->Fill(eta); }
  }

  return kTRUE;
}

void StPIDQA::FillPidHistograms(StPicoTrack *trk){

  TVector3 mTrkMom;
  mTrkMom = trk->gMom(mVertex, Bfield);

  // track variables
  double pt = mTrkMom.Perp();
  double p = mTrkMom.Mag();
  double eta = mTrkMom.PseudoRapidity();

  dEdXvp->Fill(p, trk->dEdx());

  double zpi = trk->nSigmaPion();
  z_pi->Fill(p, zpi);
  double zka = trk->nSigmaKaon();
  z_ka->Fill(p, zka);

  float tofbeta = GetTofBeta(trk);
  if (isnan(tofbeta) || tofbeta < 0) return;

  double invbeta_from_tof = 1/tofbeta;
  invbetavp->Fill(p, invbeta_from_tof);
  double norm_invbeta_pi = TMath::Sqrt(pow(Mpion,2)/pow(p,2) + 1);
  double normalisedinvbeta_for_pi = invbeta_from_tof-norm_invbeta_pi;
  normalised_invbetavpT_tof_pi->Fill(pt, normalisedinvbeta_for_pi/0.012);
  double norm_invbeta_ka = TMath::Sqrt(pow(Mkaon,2)/pow(p,2) + 1);
  double normalisedinvbeta_for_ka = invbeta_from_tof-norm_invbeta_ka;
  normalised_invbetavpT_tof_ka->Fill(pt, normalisedinvbeta_for_ka/0.012);
}

float StPIDQA::GetTofBeta(StPicoTrack *trk){
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

bool StPIDQA::IsKaon(StPicoTrack *trk){
  TVector3 mTrkMom;
  mTrkMom = trk->gMom(mVertex, Bfield);
  double p = mTrkMom.Mag();
  double pt = mTrkMom.Perp();

  float tofbeta = GetTofBeta(trk);
  if ((isnan(tofbeta) || tofbeta < 0) && pt < 1.6) return kFALSE;
  if ((isnan(tofbeta) || tofbeta < 0) && pt >= 1.6) {
    if (trk->isBemcMatchedTrack()) return (abs(trk->nSigmaKaon()) < 2);
    else return kFALSE;
  }

  double fres =  0.884 + 0.0174/pow(p + 0.0839, 4.23);
  double fpos = 0.0316 + 0.00137/pow(p + 0.101, 6.89);

  float oneOverBetaExpected = sqrt(M_KAON_PLUS*M_KAON_PLUS / p / p + 1);
  double nsigma = (1./tofbeta - oneOverBetaExpected)/0.012;

  if (nsigma < -2.0*fres + fpos || nsigma > 2.0*fres + fpos) return kFALSE;
  return kTRUE;
}

bool StPIDQA::IsPion(StPicoTrack *trk){
  return (abs(trk->nSigmaPion()) < 2.0);
}

// bool StPIDQA::IsPion(StPicoTrack *trk){ // This is the same as the pp500 version. Instead, I am using all TPC pions
//   TVector3 mTrkMom;
//   mTrkMom = trk->gMom(mVertex, Bfield);
//   double p = mTrkMom.Mag();
//   double pt = mTrkMom.Perp();

//   float tofbeta = GetTofBeta(trk);
//   if ((isnan(tofbeta) || tofbeta < 0) && pt < 1.6) return kFALSE;
//   if ((isnan(tofbeta) || tofbeta < 0) && pt >= 1.6) {
//     if (trk->isBemcMatchedTrack()) return (abs(trk->nSigmaPion()) < 2);
//     else return kFALSE;
//   }

//   float oneOverBetaExpected = sqrt(M_PION_PLUS*M_PION_PLUS / p / p + 1);
//   double nsigma = (1./tofbeta - oneOverBetaExpected)/0.012;

//   if (nsigma < -1.9 || nsigma > 2.1) return kFALSE;
//   return kTRUE;
// }

Double_t StPIDQA::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2){

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
