// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StJetTreeUpdatedMaker.h"
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
#include "StGenerateARandomTrack.h"




// Bichsel includes
#include "StBichsel/Bichsel.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StJetTreeUpdatedMaker)

//________________________________________________________________________
StJetTreeUpdatedMaker::StJetTreeUpdatedMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StJetTreeUpdatedMaker::fRunFlagEnum
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

  mD0JetMaker = 0x0;
  mCSJetMaker1 = 0x0;
  mCSJetMaker2 = 0x0;

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  d0TrackIndices.clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};
  fJetTreeCS1 = {0};
  fJetTreeCS2 = {0};

  fEmbedSingleParticle = kFALSE;
}

//
//________________________________________________________________________
StJetTreeUpdatedMaker::~StJetTreeUpdatedMaker()
{ /*  */
  
  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StJetTreeUpdatedMaker::Init() {
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
Int_t StJetTreeUpdatedMaker::Finish() { 
  cout << "StJetTreeUpdatedMaker::Finish()\n";

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

  cout<<"End of StJetTreeUpdatedMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}


void StJetTreeUpdatedMaker::DeclareTree() {
  TString treename = "D0Jets";
  jettree = new TTree(treename.Data(), treename.Data());
}

void StJetTreeUpdatedMaker::WriteTree() {
  jettree->Write();
}

//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StJetTreeUpdatedMaker::Clear(Option_t *opt) {
  fJets->Clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};

}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StJetTreeUpdatedMaker::Make() {
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

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  // if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


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

  mD0Tagger = static_cast<StTagD0Events*>(GetMaker("D0Tagger"));
  if(!mD0Tagger) {
      LOG_WARN << " No D0Tagger! Skip! " << endm;
  return kStWarn;
  }

  if (!mD0Tagger->DoesEventHaveD0()) return kStOK;

  d0TrackIndices = mD0Tagger->GetD0Indices();

  if (d0TrackIndices.size() == 0) return kStOK;

  mD0JetMaker = static_cast<StD0EventsJetMaker*>(GetMaker("D0JetMaker"));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!mD0JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", "D0JetMaker") << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  if(mD0JetMaker) {
    fJetsArr =  mD0JetMaker->GetJets();
  }
  if(fJetsArr.size()==0) return kStOk;

  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
  D0RhoMaker = static_cast<StD0Rho*>(GetMaker("D0Rho"));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!D0RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", "D0Rho") << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  rhovalue = D0RhoMaker->GetArrayofRhoValue();
  if(rhovalue.size()==0) {
    LOG_WARN << Form("Couldn't get rhovalue object! ") << endm;

    for (int i = 0; i < fJetsArr.size(); i++){
      rhovalue.push_back(0);
    }
    // return kStWarn;    
  }

  mCSJetMaker1 = static_cast<StCSJetMaker*>(GetMaker("CSJetMaker1"));
  if (!mCSJetMaker1) {
    LOG_WARN << " No CSJetMaker1! Skip! " << endm;
    return kStWarn;
  }

  fJetsArrCS1 =  mCSJetMaker1->GetJets();
  rhovalueCS1 = mCSJetMaker1->GetArrayofRhoValue();
  sigmavalueCS1 = mCSJetMaker1->GetArrayofSigmaValue();

  mCSJetMaker2 = static_cast<StCSJetMaker*>(GetMaker("CSJetMaker2"));
  if (mCSJetMaker2) {
    // 
    // return kStWarn;

    fJetsArrCS2 =  mCSJetMaker2->GetJets(); 
    rhovalueCS2 = mCSJetMaker2->GetArrayofRhoValue();
    sigmavalueCS2 = mCSJetMaker2->GetArrayofSigmaValue();
  }
  else{
    LOG_WARN << " No CSJetMaker2! Skip! " << endm;
  }

  RunJets();


  return kStOK;
}

//
//
//_____________________________________________________________________________________________

void StJetTreeUpdatedMaker::BookTree()
{
  // Branches to save event info
  jettree->Branch("RunID", &fJetTree.runid, "runid/I");
  jettree->Branch("EventId", &fJetTree.eventid, "eventid/I");
  jettree->Branch("RefMult", &fJetTree.refmult, "refmult/F");
  jettree->Branch("gRefMult", &fJetTree.grefmult, "grefmult/F");
  jettree->Branch("Centrality", &fJetTree.centrality, "centrality/F");
  jettree->Branch("Weight", &fJetTree.weight, "weight/F");
  jettree->Branch("Triggers", &fJetTree.triggers);
  jettree->Branch("PrimaryVertex", &fJetTree.primaryvertex);
  jettree->Branch("PrimaryVertexErr", &fJetTree.primaryvertexerror);

  // Branches to D0 Info

  jettree->Branch("D0Mass", &fJetTree.d0mass, "d0mass/F");
  jettree->Branch("D0Pt", &fJetTree.d0pt, "d0pt/F");
  jettree->Branch("D0Eta", &fJetTree.d0eta, "d0eta/F");
  jettree->Branch("D0Phi", &fJetTree.d0phi, "d0phi/F");

  jettree->Branch("PionPt", &fJetTree.pionpt, "pionpt/F");
  jettree->Branch("PionEta", &fJetTree.pioneta, "pioneta/F");
  jettree->Branch("PionPhi", &fJetTree.pionphi, "pionphi/F");
  jettree->Branch("PionCharge", &fJetTree.pioncharge, "pioncharge/F");

  jettree->Branch("KaonPt", &fJetTree.kaonpt, "kaonpt/F");
  jettree->Branch("KaonEta", &fJetTree.kaoneta, "kaoneta/F");
  jettree->Branch("KaonPhi", &fJetTree.kaonphi, "kaonphi/F");
  jettree->Branch("KaonCharge", &fJetTree.kaoncharge, "kaoncharge/F");

  // Branch to save jet info after area-based background subtraction

  jettree->Branch("JetPt", &fJetTree.jetpt, "jetpt/F");
  jettree->Branch("JetEta", &fJetTree.jeteta, "jeteta/F");
  jettree->Branch("JetPhi", &fJetTree.jetphi, "jetphi/F");

  // If the information above is different from different methods, we got a problem.

  // Branch to save jet info after area-based background subtraction

  jettree->Branch("JetCorrPt", &fJetTree.jetcorrectedpt, "jetcorrectedpt/F");
  jettree->Branch("JetArea", &fJetTree.jetarea, "jetarea/F");
  jettree->Branch("JetE", &fJetTree.jetenergy, "jetenergy/F");
  jettree->Branch("JetNEF", &fJetTree.jetnef, "jetNEF/F");
  jettree->Branch("JetRhoVal", &fJetTree.fRhoValforjet, "fRhoVal/F");
  jettree->Branch("JetNConst", &fJetTree.numberofconstituents, "numberofconstituents/I");

  // Branch to save jet info after CS-based background subtraction (with D0 uncorrected)
  jettree->Branch("JetCSCorrPt", &fJetTreeCS1.jetcorrectedpt, "jetcscorrectedpt/F");
  jettree->Branch("JetCSEta", &fJetTreeCS1.jeteta, "jetcs1eta/F");
  jettree->Branch("JetCSPhi", &fJetTreeCS1.jetphi, "jetcs1phi/F");
  jettree->Branch("JetCSArea", &fJetTreeCS1.jetarea, "jetcsarea/F");
  jettree->Branch("JetCSE", &fJetTreeCS1.jetenergy, "jetcsenergy/F");
  jettree->Branch("JetCSNEF", &fJetTreeCS1.jetnef, "jetcsNEF/F");
  jettree->Branch("JetCSRhoVal", &fJetTreeCS1.fRhoValforjet, "fcsRhoVal/F");
  jettree->Branch("JetCSSigVal", &fJetTreeCS1.fSigmaValforjet, "fcsSigma/F");
  jettree->Branch("JetCSNConst", &fJetTreeCS1.numberofconstituents, "numberofcsconstituents/I");

  // Branch to save jet info after CS-based background subtraction (with KPi part of the full event)
  // In this case, we need to tag D0 jets differently from the rest of the jets
  // In fact, in this case, the jet original pT will also be different.
  
  jettree->Branch("JetCS2CorrPt", &fJetTreeCS2.jetcorrectedpt, "jetcs2correctedpt/F");
  jettree->Branch("JetCS2Eta", &fJetTreeCS2.jeteta, "jetcs2eta/F");
  jettree->Branch("JetCS2Phi", &fJetTreeCS2.jetphi, "jetcs2phi/F");
  jettree->Branch("JetCS2D0DeltaR", &fJetTreeCS2.d0deltar, "jetcs2d0deltar/F");
  jettree->Branch("JetCS2Area", &fJetTreeCS2.jetarea, "jetcs2area/F");
  jettree->Branch("JetCS2E", &fJetTreeCS2.jetenergy, "jetcs2energy/F");
  jettree->Branch("JetCS2NEF", &fJetTreeCS2.jetnef, "jetcs2NEF/F");
  jettree->Branch("JetCS2RhoVal", &fJetTreeCS2.fRhoValforjet, "fcs2RhoVal/F");
  jettree->Branch("JetCS2SigVal", &fJetTreeCS2.fSigmaValforjet, "fcs2Sigma/F");
  jettree->Branch("JetCS2NConst", &fJetTreeCS2.numberofconstituents, "numberofcs2constituents/I");
  
}


void StJetTreeUpdatedMaker::RunJets()
{
  
  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables

  // cout << Form("=============== %s ===============", GetName()) << endl;


  for (int jetcollection = 0; jetcollection < d0TrackIndices.size(); jetcollection++){

    fJetTree = {0};

    // Event level information filler

    fJetTree.runid = mPicoEvent->runId();
    fJetTree.eventid = mPicoEvent->eventId();
    fJetTree.refmult = mCentMaker->GetrefMult();
    fJetTree.grefmult = mCentMaker->GetgrefMult();
    fJetTree.centrality = mCentMaker->GetCent9();
    fJetTree.weight = mCentMaker->GetWeight();
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

    fJetTree.d0mass = d0TrackIndices[jetcollection][2];

    iTrack1 = d0TrackIndices[jetcollection][0];
    iTrack2 = d0TrackIndices[jetcollection][1];

    StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack1));
    if(!trk1){ continue; }

    StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack2));
    if(!trk2){ continue; }

    TVector3 mTrk1Mom, mTrk2Mom;
    mTrk1Mom = trk1->gMom(mVertex, Bfield);
    mTrk2Mom = trk2->gMom(mVertex, Bfield);

    fJetTree.pionpt = mTrk1Mom.Perp();
    fJetTree.pioneta = mTrk1Mom.PseudoRapidity();
    fJetTree.pionphi = standardPhi(mTrk1Mom.Phi());
    fJetTree.pioncharge = trk1->charge();

    fJetTree.kaonpt = mTrk2Mom.Perp();
    fJetTree.kaoneta = mTrk2Mom.PseudoRapidity();
    fJetTree.kaonphi = standardPhi(mTrk2Mom.Phi());
    fJetTree.kaoncharge = trk2->charge();
    
    TVector3 mTrkD0Mom;
    mTrkD0Mom = mTrk1Mom + mTrk2Mom;

    fJetTree.d0pt = mTrkD0Mom.Perp();
    fJetTree.d0eta = mTrkD0Mom.PseudoRapidity();
    fJetTree.d0phi = standardPhi(mTrkD0Mom.Phi());

    // Jet TClonesArray
    fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it

    fJets = fJetsArr[jetcollection];
    fRhoVal = rhovalue[jetcollection];

    Int_t njets = fJets->GetEntries();

    bool d0Jet = kFALSE;

    // Area based background subtracted loop
    // loop over jets

    int D0JetIndex = -1;

    for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
      // get jet pointer
      StJet *jet = static_cast<StJet*>(fJets->At(ijet));
      // if(!jet) continue;

      d0Jet = kFALSE;
      D0JetIndex = -1;

      for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
        int trackid = jet->TrackAt(itrk);      

        if (trackid == 99999) {
          // cout << "D0 Jet Found." << endl;
          d0Jet = kTRUE;
          break;
        }
      }

      if (!d0Jet) continue;

      D0JetIndex = ijet;
      break;
    }

    // cout << GetName() <<  ":: D0JetIndex from Area Based Subtraction = " << D0JetIndex << "\t" << fJetTree.d0pt << endl;

    if (D0JetIndex != -1){

      StJet *jet = static_cast<StJet*>(fJets->At(D0JetIndex));
      
      // Jet level information filler

      fJetTree.jetpt = jet->Pt();
      fJetTree.jetcorrectedpt = jet->Pt() - jet->Area()*fRhoVal;
      fJetTree.jeteta = jet->Eta();
      fJetTree.jetphi = jet->Phi();
      fJetTree.jetarea = jet->Area();
      fJetTree.jetenergy = jet->E();
      fJetTree.jetnef = jet->NEF();
      fJetTree.fRhoValforjet = fRhoVal;
      fJetTree.jethighesttrackpt = jet->GetMaxTrackPt();
      fJetTree.numberofconstituents = jet->GetNumberOfTracks();

    }

    // CS1 Information Filler

    D0JetIndex = -1;
    fJetTreeCS1 = {0};

    fJetsCS1 = fJetsArrCS1[jetcollection];

    int njetsCS1 = fJetsCS1->GetEntries();

    D0JetIndex = (njetsCS1 == 1) ? 0 : -1;
    
    // for(int ijet = 0; ijet < njetsCS1; ijet++) {  // JET LOOP
    //   // get jet pointer
    //   StJet *jet = static_cast<StJet*>(fJetsCS1->At(ijet));
    //   // if(!jet) continue;


    //   D0JetIndex = -1;
    //   d0Jet = kFALSE;

    //   vector<fastjet::PseudoJet> fConstituents = jet->GetJetConstituents();
    //   // cout << Form("%s:: Constituent Size = %i", GetName(), fConstituents.size()) << endl;

    //   for(int itrk = 0; itrk < fConstituents.size(); itrk++) {
    //     int trackid = fConstituents[itrk].user_index();      

    //     if (trackid == 99999) {
    //       // cout << "D0 Jet Found From CS." << endl;

    //       // cout << "NConstituents CS1 = " << fConstituents.size() << "\t" << jet->GetNumberOfTracks() << "\t" << jet->GetNumberOfTowers() << endl;
    //       // cout << "Actual D0 pT = " << fJetTree.d0pt << "\t" << "BG Subtracted D0 pT = " << fConstituents[itrk].perp() << endl;

    //       d0Jet = kTRUE;
    //       break;
    //     }
    //   }

    //   if (!d0Jet) continue;

    //   D0JetIndex = ijet;
    //   break;
    // }

    if (D0JetIndex != -1){
      StJet *jetCS1 = static_cast<StJet*>(fJetsCS1->At(D0JetIndex));

      fJetTreeCS1.jetcorrectedpt = jetCS1->Pt();
      fJetTreeCS1.jeteta = jetCS1->Eta();
      fJetTreeCS1.jetphi = jetCS1->Phi();
      fJetTreeCS1.jetarea = jetCS1->Area();
      fJetTreeCS1.jetenergy = jetCS1->E();
      fJetTreeCS1.jetnef = jetCS1->NEF();
      fJetTreeCS1.fRhoValforjet = rhovalueCS1[jetcollection];
      fJetTreeCS1.fSigmaValforjet = sigmavalueCS1[jetcollection];
      fJetTreeCS1.jethighesttrackpt = jetCS1->GetMaxTrackPt();
      fJetTreeCS1.numberofconstituents = jetCS1->GetNumberOfTracks();

      double deltaR = TMath::Sqrt(pow(jetCS1->Eta() - fJetTree.d0eta, 2) + pow(dPhi(jetCS1->Phi(), fJetTree.d0phi), 2));

      // cout << "Event ID = " << mPicoEvent->eventId() << "\t" << " Index = " << D0JetIndex << " Delta R = " << deltaR << "\t" << " Jet pT = " << fJetTree.jetpt << "\t" << fJetTreeCS1.jetcorrectedpt << "\t" << fJetTree.pionpt << "\t" << fJetTree.kaonpt << endl;

      // if (fJetTreeCS1.jeteta != fJetTree.jeteta) cout << "Jet Eta Different: Halt STOP " << fJetTree.jeteta << "\t" << fJetTreeCS1.jeteta << endl;
    }

    // CS2 Information Filler

    if (mCSJetMaker2){

      D0JetIndex = -1;
      fJetTreeCS2 = {0};

      fJetsCS2 = fJetsArrCS2[jetcollection];

      int njetsCS2 = fJetsCS2->GetEntries();

      // cout << "CSJetMaker2:: Entries = " << njetsCS2 << endl;

      double mindeltaR = 999.0;
    
      for(int ijet = 0; ijet < njetsCS2; ijet++) {  // JET LOOP
        // get jet pointer
        StJet *jet = static_cast<StJet*>(fJetsCS2->At(ijet));
        if(!jet) continue;

        double deltaR = TMath::Sqrt(pow(jet->Eta() - fJetTree.d0eta, 2) + pow(dPhi(jet->Phi(), fJetTree.d0phi), 2));

        if (deltaR < mindeltaR){
          mindeltaR = deltaR;
          D0JetIndex = ijet;
        }
      }

      if (mindeltaR <= 1.0 && D0JetIndex != -1){
        StJet *jetCS2 = static_cast<StJet*>(fJetsCS2->At(D0JetIndex));

        fJetTreeCS2.jetcorrectedpt = jetCS2->Pt();
        fJetTreeCS2.jeteta = jetCS2->Eta();
        fJetTreeCS2.jetphi = jetCS2->Phi();
        fJetTreeCS2.jetarea = jetCS2->Area();
        fJetTreeCS2.d0deltar = mindeltaR;
        fJetTreeCS2.jetenergy = jetCS2->E();
        fJetTreeCS2.jetnef = jetCS2->NEF();
        fJetTreeCS2.fRhoValforjet = rhovalueCS2[jetcollection];
        fJetTreeCS2.fSigmaValforjet = sigmavalueCS2[jetcollection];
        fJetTreeCS2.jethighesttrackpt = jetCS2->GetMaxTrackPt();
        fJetTreeCS2.numberofconstituents = jetCS2->GetNumberOfTracks();
      }
    }

    jettree->Fill();

  }

}

Double_t StJetTreeUpdatedMaker::dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

Double_t StJetTreeUpdatedMaker::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}


