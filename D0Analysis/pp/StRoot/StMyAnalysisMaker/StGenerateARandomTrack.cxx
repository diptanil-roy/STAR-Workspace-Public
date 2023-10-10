// ################################################################
// Author: Diptanil Roy
// Based on Joel Mazer's framework for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StGenerateARandomTrack.h"
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

ClassImp(StGenerateARandomTrack)

//________________________________________________________________________
StGenerateARandomTrack::StGenerateARandomTrack(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{ 
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StGenerateARandomTrack::fRunFlagEnum
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
StGenerateARandomTrack::~StGenerateARandomTrack()
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
  }

  if(mEmcPosition) delete mEmcPosition;

}

//________________________________________________________________________
Int_t StGenerateARandomTrack::Init() {
  StJetFrameworkPicoBase::Init();

  DeclareHistograms();

  //position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StGenerateARandomTrack::Finish() { 

  cout << "StGenerateARandomTrack::Finish()\n";

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

  cout<<"End of StGenerateARandomTrack::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StGenerateARandomTrack::DeclareHistograms() {

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

  }
}

//
// write histograms
//_____________________________________________________________________________
void StGenerateARandomTrack::WriteHistograms() {

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
  }
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StGenerateARandomTrack::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StGenerateARandomTrack::Make() {
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

  fSingleParticleVector.clear();

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
    if (weight > 0) hCentralityWeightedBeforeCuts->Fill(fCentralityScaled, weight);
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

    if (abs(zVtx) > 6.) return kStOk;
    numberofevents[6]++;

    if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > 2.) return kStOK;
    numberofevents[7]++;

    if (abs(zVtx - zVtx_VPD) > 3) return kStOk;
    numberofevents[8]++;

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

  // cout << "Number of tracks found = " << ntracks << "\t" << nglobaltracks << endl;

  RunTracks();

  return kStOK;
}

//__________________________________________________________________________________________
  
void StGenerateARandomTrack::RunTracks(){

  
  fSingleParticleVector.clear();

  double pi0mass = Pico::mMass[0]; // GeV

  double pt_  = gRandom->Uniform(1,40);
  double phi_ = gRandom->Uniform(0, 2*TMath::Pi());
  double eta_ = gRandom->Uniform(-1, 1);

  double px_ = pt_*TMath::Cos(phi_);
  double py_ = pt_*TMath::Sin(phi_);
  double pz_ = pt_ * TMath::SinH(eta_);
  double p_ = sqrt(px_*px_ + py_*py_ + pz_*pz_);
  double pe_ = sqrt(p_*p_+pi0mass*pi0mass);  

  // cout << pt_ << "\t" << eta_ << "\t" << phi_ << "\t" <<  pe_ << endl;    
  TLorentzVector fSingleParticle;
  fSingleParticle.SetPxPyPzE(px_, py_, pz_, pe_);
  fSingleParticleVector.push_back(fSingleParticle);
}

