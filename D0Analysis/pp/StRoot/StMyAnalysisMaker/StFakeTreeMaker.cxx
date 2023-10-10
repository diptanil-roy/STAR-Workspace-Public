// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StFakeTreeMaker.h"
#include "StMemStat.h"


// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <vector>

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
#include "StTagD0Events.h"
#include "StThrowARandomTrack.h"
#include "StD0EventsJetMaker.h"
#include "StD0Rho.h"


// Bichsel includes
#include "StBichsel/Bichsel.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StFakeTreeMaker)

//________________________________________________________________________
StFakeTreeMaker::StFakeTreeMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  fD0Kind = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StFakeTreeMaker::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;

  fMinCentrality = 0;
  fMaxCentrality = 80;

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
  frefCorr2 = -99.;
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
  fRhoMakerName = rhoMakerName;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  d0TrackIndices.clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};
}

//
//________________________________________________________________________
StFakeTreeMaker::~StFakeTreeMaker()
{ /*  */
  
  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StFakeTreeMaker::Init() {
  StJetFrameworkPicoBase::Init();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);

  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};

  DeclareTree();
  BookTree();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StFakeTreeMaker::Finish() { 
  cout << "StFakeTreeMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    // WriteHistograms();

    WriteTree();

    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StFakeTreeMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}


void StFakeTreeMaker::DeclareTree() {
  TString treename;

  treename = "FakeJetTree";

  jettree = new TTree(treename.Data(), treename.Data());
}

void StFakeTreeMaker::WriteTree() {
  jettree->Write();
}

//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StFakeTreeMaker::Clear(Option_t *opt) {
  fJets->Clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};

}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StFakeTreeMaker::Make() {
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
  frefCorr2 = refCorr2;
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  // FillEmcTriggers();

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

  // =========================== JetMaker =============================== //
  // get JetMaker pointer
  // JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  // const char *fJetMakerNameCh = fJetMakerName;
  // if(!JetMaker) {
  //   LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
  //   return kStWarn;
  // }

  // // get jet collection associated with JetMaker
  // fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  // if(!fJets) {
  //   LOG_WARN << Form(" No fJets object! Skip! ") << endm;
  //   return kStWarn;
  // }

  // get JetMaker


  mRandomTrack = static_cast<StThrowARandomTrack*>(GetMaker("RandomTrack"));
  if(!mRandomTrack) {
    LOG_WARN << " No D0Tagger! Skip! " << endm;
    return kStWarn;
  }

  // ======================== end of Triggers ============================= //

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

  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
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
  fRhoVal = fRho->GetVal();
  // =======================================================================

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();
  //   }
  // }

  // run Jets:
  // RunJets();

  SaveJets(9998);
  SaveJets(9999);

  // // run Tracks:
  // RunTracks();

  // // run Towers:
  // RunTowers();

  return kStOK;
}

//
//
//_____________________________________________________________________________________________

void StFakeTreeMaker::BookTree()
{
  // Branches to save event info
  jettree->Branch("RunID", &fJetTree.runid, "runid/I");
  jettree->Branch("EventId", &fJetTree.eventid, "eventid/I");
  jettree->Branch("RefMult", &fJetTree.refmult, "refmult/F");
  jettree->Branch("Centrality", &fJetTree.centrality, "centrality/F");
  jettree->Branch("Triggers", &fJetTree.triggers);
  jettree->Branch("PrimaryVertex", &fJetTree.primaryvertex);
  jettree->Branch("PrimaryVertexErr", &fJetTree.primaryvertexerror);

  jettree->Branch("HardPt", &fJetTree.hardpt, "hardpt/F");
  jettree->Branch("HardEta", &fJetTree.hardeta, "hardeta/F");
  jettree->Branch("HardPhi", &fJetTree.hardphi, "hardphi/F");

  jettree->Branch("JetPt", &fJetTree.jetpt, "jetpt/F");
  jettree->Branch("JetCorrPt", &fJetTree.jetcorrectedpt, "jetcorrectedpt/F");
  jettree->Branch("JetEta", &fJetTree.jeteta, "jeteta/F");
  jettree->Branch("JetPhi", &fJetTree.jetphi, "jetphi/F");
  jettree->Branch("JetArea", &fJetTree.jetarea, "jetarea/F");
  jettree->Branch("JetRadius", &fJetTree.jetradius, "jetradius/F");
  jettree->Branch("JetE", &fJetTree.jetenergy, "jetenergy/F");
  jettree->Branch("JetNEF", &fJetTree.jetnef, "jetnef/F");
  jettree->Branch("JetRhoVal", &fJetTree.fRhoValforjet, "fRhoValforjet/F");
  jettree->Branch("JetHighestTrackPt", &fJetTree.jethighesttrackpt, "jethighesttrackpt/F");
  jettree->Branch("JetNConst", &fJetTree.numberofconstituents, "numberofconstituents/I");

  // Branches to save jet data info
  jettree->Branch("TrackID", &fJetTree.mTrackID, "mTrackID[numberofconstituents]/F");
  jettree->Branch("TrackPt", &fJetTree.mTrackPt, "mTrackPt[numberofconstituents]/F");
  jettree->Branch("TrackEta", &fJetTree.mTrackEta, "mTrackEta[numberofconstituents]/F");
  jettree->Branch("TrackPhi", &fJetTree.mTrackPhi, "mTrackPhi[numberofconstituents]/F");
  jettree->Branch("TrackPx", &fJetTree.mTrackPx, "mTrackPx[numberofconstituents]/F");
  jettree->Branch("TrackPy", &fJetTree.mTrackPy, "mTrackPy[numberofconstituents]/F");
  jettree->Branch("TrackPz", &fJetTree.mTrackPz, "mTrackPz[numberofconstituents]/F");
  jettree->Branch("TrackCharge", &fJetTree.mTrackCharge, "mTrackCharge[numberofconstituents]/F");
}


void StFakeTreeMaker::SaveJets(Int_t FakeTrackID = 9998)
{
  
  fJetTree = {0};

  TVector3 hardtrack = mRandomTrack->GetRandomHardTrack();

  if (FakeTrackID == 9999) hardtrack = -hardtrack;

  // track variables
  double pt = hardtrack.Perp();
  double phi = hardtrack.Phi();
  if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
  if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
  double eta = hardtrack.PseudoRapidity();

  fJetTree.runid = mPicoEvent->runId();
  fJetTree.eventid = mPicoEvent->eventId();
  fJetTree.refmult = frefCorr2;
  fJetTree.centrality = fCentralityScaled;
  fJetTree.triggers = mPicoEvent->triggerIds();

  vector<double> pv;
  pv.clear();
  pv.push_back(mPicoEvent->primaryVertex().X());
  pv.push_back(mPicoEvent->primaryVertex().Y());
  pv.push_back(mPicoEvent->primaryVertex().Z());

  fJetTree.primaryvertex = pv;

  vector<double> pverr;
  pverr.clear();
  pverr.push_back(mPicoEvent->primaryVertexError().X());
  pverr.push_back(mPicoEvent->primaryVertexError().Y());
  pverr.push_back(mPicoEvent->primaryVertexError().Z());

  fJetTree.primaryvertexerror = pverr;

  fJetTree.hardpt = pt;
  fJetTree.hardeta = eta;
  fJetTree.hardphi = phi;

  Int_t njets = fJets->GetEntries();

  vector<double> FakeJetID;
  FakeJetID.clear();

  // loop over jets
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    // get jet pointer
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

    // get nTracks and maxTrackPt
    // double maxtrackpt = jet->GetMaxTrackPt();
    // double NtrackConstit = jet->GetNumberOfTracks();

    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);

      // if (trackid > 5000) cout << trackid << endl;
      if (trackid == FakeTrackID) {FakeJetID.push_back(ijet); break;}
    }

    if (FakeJetID.size() == 1) break;
  }

  if (FakeJetID.size() == 0){
    jettree->Fill();
    return;
  }

  else{
    StJet *jet = static_cast<StJet*>(fJets->At(FakeJetID[0]));
    

    // get some jet parameters
    double jetarea = jet->Area();
    double jetpt = jet->Pt();
    double corrjetpt = jet->Pt() - jetarea*fRhoVal;
    double jetE = jet->E();
    double jetEta = jet->Eta();
    double jetPhi = jet->Phi();
    double jetNEF = jet->NEF();

    // get nTracks and maxTrackPt
    // double maxtrackpt = jet->GetMaxTrackPt();
    double NtrackConstit = jet->GetNumberOfTracks();

    fJetTree.jetpt = jetpt;
    fJetTree.jetcorrectedpt = corrjetpt;
    fJetTree.jeteta = jetEta;
    fJetTree.jetphi = jetPhi;
    fJetTree.jetarea = jetarea;
    fJetTree.jetradius = 0.4;
    fJetTree.jetenergy = jetE;
    fJetTree.jetnef = jetNEF;
    fJetTree.fRhoValforjet = fRhoVal;
    fJetTree.jethighesttrackpt = -999;
    fJetTree.numberofconstituents = NtrackConstit;

    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);

      // fJetTree.mTrackID[itrk] = trackid;

      // cout << trackid << endl;

      if (trackid != FakeTrackID){
        StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
        if(!trk){ continue; }

        // get momentum vector
        TVector3 mTrkMom;
        if(doUsePrimTracks) {
        if(!(trk->isPrimary())) continue; // check if primary
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
        double px = mTrkMom.x();
        double py = mTrkMom.y();
        double pz = mTrkMom.z();
        short charge = trk->charge();

        // Track level information fill

        fJetTree.mTrackID[itrk] = IsWhatParticle(trk);
        fJetTree.mTrackPt[itrk] = pt;
        fJetTree.mTrackEta[itrk] = eta;
        fJetTree.mTrackPhi[itrk] = standardPhi(phi);
        fJetTree.mTrackPx[itrk] = px;
        fJetTree.mTrackPy[itrk] = py;
        fJetTree.mTrackPz[itrk] = pz;
        fJetTree.mTrackCharge[itrk] = charge;
      }

      else{
        fJetTree.mTrackID[itrk] = FakeTrackID;
        fJetTree.mTrackPt[itrk] = hardtrack.Perp();
        fJetTree.mTrackEta[itrk] = hardtrack.PseudoRapidity();
        fJetTree.mTrackPhi[itrk] = standardPhi(hardtrack.Phi());
        fJetTree.mTrackPx[itrk] = hardtrack.x();
        fJetTree.mTrackPy[itrk] = hardtrack.y();
        fJetTree.mTrackPz[itrk] = hardtrack.z();
        fJetTree.mTrackCharge[itrk] = 0;
      }


    }      

    jettree->Fill();

  }

}

void StFakeTreeMaker::IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e){ // NEW PID APPROACH

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

  if ((abs(zpi) > 2) && (abs(zka) > 2)) {pid = 0; m = 0.; e = 0.; return;}

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

  if(toftrack){

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

    if ((abs(normalisedinvbeta_for_pi) > 2) && (abs(normalisedinvbeta_for_ka) > 2)) {pid = 0; m = 0.; e = 0.; return;}

    tof_pion = (abs(normalisedinvbeta_for_pi) <= abs(normalisedinvbeta_for_ka)) ? kTRUE : kFALSE;
    tof_kaon = (abs(normalisedinvbeta_for_pi) > abs(normalisedinvbeta_for_ka)) ? kTRUE : kFALSE;



    if (tof_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tof_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}

  }

  else{
    if (tpc_pion) {pid = 1*charge; m = Mpion; e = TMath::Sqrt(pow(p,2) + pow(Mpion, 2)); return;}
    else if (tpc_kaon) {pid = 2*charge; m = Mkaon; e = TMath::Sqrt(pow(p,2) + pow(Mkaon, 2)); return;}
    else {pid = 0; m = 0.; e = 0.; return;}
  }
}


Int_t StFakeTreeMaker::IsWhatParticle(StPicoTrack *trk){ // Just to get the PID out
  int pid;
  double m; 
  double e;
  IsWhatParticle(trk, pid, m, e);
  return pid;
}

Double_t StFakeTreeMaker::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}


