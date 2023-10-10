// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StLeadJetChargeCorrelator.h"
#include "StMemStat.h"
#include "phys_constants.h"

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
#include "StEmcPosition2.h"
#include "StCentMaker.h"

#include "StPhysicalHelix.hh"
#include "StThreeVectorF.hh"
#include "StBTofUtil/tofPathLength.hh"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StLeadJetChargeCorrelator)

//________________________________________________________________________
StLeadJetChargeCorrelator::StLeadJetChargeCorrelator(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StLeadJetChargeCorrelator::fRunFlagEnum
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
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  mChargeCorrStruct = {0};

}

//
//________________________________________________________________________
StLeadJetChargeCorrelator::~StLeadJetChargeCorrelator()
{ /*  */
  // destructor
  if(mEmcPosition) delete mEmcPosition;
}

void StLeadJetChargeCorrelator::DeclareTree() {
  TString treename = "ChargeCorrelation";

  TDirectory *writedir;

  if(mOutName!="") {
    fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    writedir = (TDirectory *)fout->Get(GetName());
    // fout->Close();
  }

  jettree = new TTree(treename.Data(), treename.Data());
  jettree->SetDirectory(writedir);
}

//
//________________________________________________________________________
Int_t StLeadJetChargeCorrelator::Init() {
  StJetFrameworkPicoBase::Init();

  // book Histograms
  mChargeCorrStruct = {0};

  DeclareTree();
  BookTree();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StLeadJetChargeCorrelator::Finish() { 
  cout << "StLeadJetChargeCorrelator::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StLeadJetChargeCorrelator::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

void StLeadJetChargeCorrelator::BookTree()
{
  // Branches to save event info
  jettree->Branch("RunID", &mChargeCorrStruct.runid, "runid/I");
  jettree->Branch("EventId", &mChargeCorrStruct.eventid, "eventid/I");
  jettree->Branch("RefMult", &mChargeCorrStruct.refmult, "refmult/F");
  jettree->Branch("gRefMult", &mChargeCorrStruct.grefmult, "grefmult/F");
  jettree->Branch("Centrality", &mChargeCorrStruct.centrality, "centrality/F");
  jettree->Branch("Weight", &mChargeCorrStruct.weight, "weight/F");
  jettree->Branch("IsMB", &mChargeCorrStruct.isMB, "isMB/O");
  jettree->Branch("IsHT1", &mChargeCorrStruct.isHT1, "isHT1/O");
  jettree->Branch("IsHT2", &mChargeCorrStruct.isHT2, "isHT2/O");
  jettree->Branch("IsHT3", &mChargeCorrStruct.isHT3, "isHT3/O");


  jettree->Branch("JetPt", &mChargeCorrStruct.jetpt, "jetpt/F");
  jettree->Branch("JetEta", &mChargeCorrStruct.jeteta, "jeteta/F");
  jettree->Branch("JetPhi", &mChargeCorrStruct.jetphi, "jetphi/F");
  jettree->Branch("JetArea", &mChargeCorrStruct.jetarea, "jetarea/F");
  jettree->Branch("JetE", &mChargeCorrStruct.jetenergy, "jetenergy/F");
  jettree->Branch("JetNEF", &mChargeCorrStruct.jetnef, "jetnef/F");
  jettree->Branch("JetNConst", &mChargeCorrStruct.numberofconstituents, "numberofconstituents/I");

  jettree->Branch("LeadTrackPt", &mChargeCorrStruct.leadtrackpt, "leadtrackpt/F");
  jettree->Branch("LeadTrackEta", &mChargeCorrStruct.leadtracketa, "leadtracketa/F");
  jettree->Branch("LeadTrackPhi", &mChargeCorrStruct.leadtrackphi, "leadtrackphi/F");
  jettree->Branch("LeadTrackCharge", &mChargeCorrStruct.leadtrackcharge, "leadtrackcharge/I");
  jettree->Branch("LeadTrackIsPion", &mChargeCorrStruct.leadtrackispion, "leadtrackispion/O");
  jettree->Branch("LeadTrackIsKaon", &mChargeCorrStruct.leadtrackiskaon, "leadtrackiskaon/O");
  jettree->Branch("LeadTrackIsProton", &mChargeCorrStruct.leadtrackisproton, "leadtrackisproton/O");

  jettree->Branch("SubLeadTrackPt", &mChargeCorrStruct.subleadtrackpt, "subleadtrackpt/F");
  jettree->Branch("SubLeadTrackEta", &mChargeCorrStruct.subleadtracketa, "subleadtracketa/F");
  jettree->Branch("SubLeadTrackPhi", &mChargeCorrStruct.subleadtrackphi, "subleadtrackphi/F");
  jettree->Branch("SubLeadTrackCharge", &mChargeCorrStruct.subleadtrackcharge, "subleadtrackcharge/I");
  jettree->Branch("SubLeadTrackIsPion", &mChargeCorrStruct.subleadtrackispion, "subleadtrackispion/O");
  jettree->Branch("SubLeadTrackIsKaon", &mChargeCorrStruct.subleadtrackiskaon, "subleadtrackiskaon/O");
  jettree->Branch("SubLeadTrackIsProton", &mChargeCorrStruct.subleadtrackisproton, "subleadtrackisproton/O");
}

// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StLeadJetChargeCorrelator::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StLeadJetChargeCorrelator::Make() {

  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  mChargeCorrStruct = {0};
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

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }


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
  bool fRunForMB = kFALSE;  // used to differentiate pp and AuAu
  if(doppAnalysis)  fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

  // TODO - YOU WILL NEED TO USE THIS BOOLEAN's appropriately!!! - TODO

  int arrMB_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
  int arrHT1_Run14[]  = {450201, 450211, 460201};
  int arrHT2_Run14[]  = {450202, 450212, 460202, 460212};
  int arrHT3_Run14[]  = {450203, 450213, 460203};

  matchMB = false;
  matchHT1 = false;
  matchHT2 = false;
  matchHT3 = false;

  if (!doppAnalysis){
    for(int i = 0; i < sizeof(arrMB_Run14)/sizeof(*arrMB_Run14); i++) {
      if(mPicoEvent->isTrigger(arrMB_Run14[i])) matchMB = kTRUE;
      if(matchMB) break;
    }

    if (!matchMB) return kStOk;

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

  }

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

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // run Jets:
  RunJets();

  return kStOK;
}

void StLeadJetChargeCorrelator::RunJets(){

  mChargeCorrStruct = {0}; // clear struct

  Int_t njets = fJets->GetEntries();

  for(int ijet = 0; ijet < njets; ijet++) {
    StJet *jet = static_cast<StJet*>(fJets->At(0));
    if(!jet) continue;

    vector<fastjet::PseudoJet> fConstituents = jet->GetJetConstituents();

    if (fConstituents.size() < 2) continue; // For this analysis, we only need the leading jet with at least 2 constituents

    int trackid1 = fConstituents[0].user_index(); // leading track
    int trackid2 = fConstituents[1].user_index(); // subleading track

    if (trackid1 < 0 || trackid2 < 0) continue; // if either track is neutral, skip

    mChargeCorrStruct = {0}; // clear struct

    // These info will be stored in the tree

    mChargeCorrStruct.runid = mPicoEvent->runId();
    mChargeCorrStruct.eventid = mPicoEvent->eventId();
    mChargeCorrStruct.refmult = mCentMaker->GetrefMult();
    mChargeCorrStruct.grefmult = mCentMaker->GetgrefMult();
    mChargeCorrStruct.centrality = fCentralityScaled;
    mChargeCorrStruct.weight = mCentMaker->GetWeight();

    mChargeCorrStruct.isMB = matchMB;
    mChargeCorrStruct.isHT1 = matchHT1;
    mChargeCorrStruct.isHT2 = matchHT2;
    mChargeCorrStruct.isHT3 = matchHT3;

    mChargeCorrStruct.jetarea = jet->Area();
    mChargeCorrStruct.jetpt = jet->Pt();
    mChargeCorrStruct.jetenergy = jet->E();
    mChargeCorrStruct.jeteta = jet->Eta();
    mChargeCorrStruct.jetphi = jet->Phi();
    mChargeCorrStruct.jetnef = jet->NEF();
    mChargeCorrStruct.numberofconstituents = jet->GetNumberOfTracks();

    // get jet track pointers
    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(trackid1));
    StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(trackid2));
    if(!trk1 || !trk2){ continue; }

    // get momentum vectors

    TVector3 mTrkMom1 = trk1->gMom(mVertex, Bfield);
    TVector3 mTrkMom2 = trk2->gMom(mVertex, Bfield);

    mChargeCorrStruct.leadtrackcharge = trk1->charge(); 
    mChargeCorrStruct.subleadtrackcharge = trk2->charge();

    mChargeCorrStruct.leadtracketa = mTrkMom1.Eta();
    mChargeCorrStruct.subleadtracketa = mTrkMom2.Eta();

    mChargeCorrStruct.leadtrackphi = mTrkMom1.Phi();
    mChargeCorrStruct.subleadtrackphi = mTrkMom2.Phi();

    mChargeCorrStruct.leadtrackpt = mTrkMom1.Pt();
    mChargeCorrStruct.subleadtrackpt = mTrkMom2.Pt();

    mChargeCorrStruct.leadtrackispion = IsPion(trk1);
    mChargeCorrStruct.subleadtrackispion = IsPion(trk2);

    mChargeCorrStruct.leadtrackiskaon = IsKaon(trk1);
    mChargeCorrStruct.subleadtrackiskaon = IsKaon(trk2);

    mChargeCorrStruct.leadtrackisproton = false;
    mChargeCorrStruct.subleadtrackisproton = false;

    jettree->Fill();
  }
}

bool StLeadJetChargeCorrelator::IsTpcPion(StPicoTrack *trk){
  double zpi = trk->nSigmaPion();
  if (abs(zpi) < 3.) return true;
  return false;
}

bool StLeadJetChargeCorrelator::IsTpcKaon(StPicoTrack *trk){
  double zka = trk->nSigmaKaon();
  if (abs(zka) < 2.) return true;
  return false;
}

bool StLeadJetChargeCorrelator::IsPion(StPicoTrack *trk){
  if (!IsTpcPion(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_PION_PLUS*M_PION_PLUS / p / p + 1);
  double nsigma = abs(1./beta - oneOverBetaExpected);
  if (nsigma < 0.03) return true;
  return false;
}

bool StLeadJetChargeCorrelator::IsKaon(StPicoTrack *trk){
  if (!IsTpcKaon(trk)) return false;
  float beta = GetTofBeta(trk);
  if (isnan(beta) || beta < 0) return true;
  double p = trk->gMom(mVertex, Bfield).Mag();
  float oneOverBetaExpected = sqrt(M_KAON_PLUS*M_KAON_PLUS / p / p + 1);
  double nsigma = abs(1./beta - oneOverBetaExpected);
  if (nsigma < 0.03) return true;
  return false;
}

// Get TOF Information

float StLeadJetChargeCorrelator::GetTofBeta(StPicoTrack *trk){
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
