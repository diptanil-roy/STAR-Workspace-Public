// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StJetTreeMaker.h"
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
#include "StPIDQA.h"
#include "StGenerateARandomTrack.h"
#include "StD0EventsJetMaker.h"
#include "StD0Rho.h"

#include "StPhysicalHelix.hh"
#include "StThreeVectorF.hh"

#include "StBTofUtil/tofPathLength.hh"


// Bichsel includes
#include "StBichsel/Bichsel.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StJetTreeMaker)

//________________________________________________________________________
StJetTreeMaker::StJetTreeMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StJetTreeMaker::fRunFlagEnum
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

  fEmbedSingleParticle = kFALSE;
}

//
//________________________________________________________________________
StJetTreeMaker::~StJetTreeMaker()
{ /*  */
  
  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StJetTreeMaker::Init() {
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
Int_t StJetTreeMaker::Finish() { 
  cout << "StJetTreeMaker::Finish()\n";

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

  cout<<"End of StJetTreeMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}


void StJetTreeMaker::DeclareTree() {
  TString ustreename = "D0JetsUSTree";
  d0unlikesigntree = new TTree(ustreename.Data(), ustreename.Data());

  TString lstreename = "D0JetsLSTree";
  d0likesigntree = new TTree(lstreename.Data(), lstreename.Data());
}

void StJetTreeMaker::WriteTree() {
  d0unlikesigntree->Write();
  d0likesigntree->Write();
}

//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StJetTreeMaker::Clear(Option_t *opt) {
  fJets->Clear();
  fJetsArr.clear();
  rhovalue.clear();

  fJetTree = {0};

}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StJetTreeMaker::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;
  fEmbeddedSingleParticleVector.clear();

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
  // if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

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
  
  if (!fEmbedSingleParticle){
    mD0Tagger = static_cast<StPIDQA*>(GetMaker("PIDQA"));
    if(!mD0Tagger) {
      LOG_WARN << " No D0Tagger! Skip! " << endm;
      return kStWarn;
    }

    if (!mD0Tagger->DoesEventHaveD0()) return kStOK;
    d0TrackIndices = mD0Tagger->GetD0Indices();
    if (d0TrackIndices.size() == 0) return kStOK;
  }

  else {
    mGenerateRandomParticle = static_cast<StGenerateARandomTrack*>(GetMaker("RandomParticle"));
    if(!mGenerateRandomParticle) {
      LOG_WARN << " No Random Particle Generator! Skip! " << endm;
      return kStWarn;
    }

    fEmbeddedSingleParticleVector = mGenerateRandomParticle->GetEmbeddedParticle();

    if (fEmbeddedSingleParticleVector.size() != 1) return kStOK;
  }

  D0JetMaker = static_cast<StD0EventsJetMaker*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!D0JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  if(D0JetMaker) {
    fJetsArr =  D0JetMaker->GetJets();
    //fJets->SetName("BGJetsRho");  // name is set by Maker who created it
  }
  if(fJetsArr.size()==0) return kStOk;
  
  if ( !fEmbedSingleParticle ) RunJets();
  else RunJetsForEmbeddedParticle();

  return kStOK;
}

//
//
//_____________________________________________________________________________________________

void StJetTreeMaker::BookTree()
{
  // Branches to save event info
  
  d0unlikesigntree->Branch("Vz", &fJetTree.Vz, "Vz/F");
  d0unlikesigntree->Branch("VzVPD", &fJetTree.VzVPD, "VzVPD/F");
  d0unlikesigntree->Branch("Triggers", &fJetTree.Triggers);
  d0unlikesigntree->Branch("D0Mass", &fJetTree.D0mass, "D0mass/F");
  d0unlikesigntree->Branch("D0Pt", &fJetTree.D0pt, "D0pt/F");
  d0unlikesigntree->Branch("D0Eta", &fJetTree.D0eta, "D0eta/F");
  d0unlikesigntree->Branch("D0Phi", &fJetTree.D0phi, "D0phi/F");
  d0unlikesigntree->Branch("PionPt", &fJetTree.pionpt, "pionpt/F");
  d0unlikesigntree->Branch("PionEta", &fJetTree.pioneta, "pioneta/F");
  d0unlikesigntree->Branch("PionPhi", &fJetTree.pionphi, "pionphi/F");
  d0unlikesigntree->Branch("PionDCA", &fJetTree.piondca, "piondca/F");
  d0unlikesigntree->Branch("PionNHitsFit", &fJetTree.pionnhitsfit, "pionnhitsfit/F");
  d0unlikesigntree->Branch("PionNSigmaPion", &fJetTree.pionnsigmapion, "pionnsigmapion/F");
  d0unlikesigntree->Branch("PionNSigmaKaon", &fJetTree.pionnsigmakaon, "pionnsigmakaon/F");
  d0unlikesigntree->Branch("PionTofBeta", &fJetTree.piontofbeta, "piontofbeta/F");
  d0unlikesigntree->Branch("PionTofYLocal", &fJetTree.ntofpionylocal, "ntofpionylocal/F");
  d0unlikesigntree->Branch("KaonPt", &fJetTree.kaonpt, "kaonpt/F");
  d0unlikesigntree->Branch("KaonEta", &fJetTree.kaoneta, "kaoneta/F");
  d0unlikesigntree->Branch("KaonPhi", &fJetTree.kaonphi, "kaonphi/F");
  d0unlikesigntree->Branch("KaonDCA", &fJetTree.kaondca, "kaondca/F");
  d0unlikesigntree->Branch("KaonNHitsFit", &fJetTree.kaonnhitsfit, "kaonnhitsfit/F");
  d0unlikesigntree->Branch("KaonNSigmaPion", &fJetTree.kaonnsigmapion, "kaonnsigmapion/F");
  d0unlikesigntree->Branch("KaonNSigmaKaon", &fJetTree.kaonnsigmakaon, "kaonnsigmakaon/F");
  d0unlikesigntree->Branch("KaonTofBeta", &fJetTree.kaontofbeta, "kaontofbeta/F");
  d0unlikesigntree->Branch("KaonTofYLocal", &fJetTree.ntofkaonylocal, "ntofkaonylocal/F");

  d0unlikesigntree->Branch("JetPt", &fJetTree.jetpt, "jetpt/F");
  d0unlikesigntree->Branch("JetEta", &fJetTree.jeteta, "jeteta/F");
  d0unlikesigntree->Branch("JetPhi", &fJetTree.jetphi, "jetphi/F");
  d0unlikesigntree->Branch("JetArea", &fJetTree.jetarea, "jetarea/F");
  d0unlikesigntree->Branch("JetE", &fJetTree.jetenergy, "jetenergy/F");
  d0unlikesigntree->Branch("JetNEF", &fJetTree.jetnef, "jetnef/F");
  d0unlikesigntree->Branch("JetNConst", &fJetTree.numberofconstituents, "numberofconstituents/I");

  // Branches to save jet data info
  d0unlikesigntree->Branch("TrackPt", &fJetTree.mTrackPt, "mTrackPt[numberofconstituents]/F");
  d0unlikesigntree->Branch("TrackEta", &fJetTree.mTrackEta, "mTrackEta[numberofconstituents]/F");
  d0unlikesigntree->Branch("TrackPhi", &fJetTree.mTrackPhi, "mTrackPhi[numberofconstituents]/F");
  d0unlikesigntree->Branch("TrackCharge", &fJetTree.mTrackCharge, "mTrackCharge[numberofconstituents]/F");

  d0likesigntree->Branch("Vz", &fJetTree.Vz, "Vz/F");
  d0likesigntree->Branch("VzVPD", &fJetTree.VzVPD, "VzVPD/F");
  d0likesigntree->Branch("Triggers", &fJetTree.Triggers);
  d0likesigntree->Branch("D0Mass", &fJetTree.D0mass, "D0mass/F");
  d0likesigntree->Branch("D0Pt", &fJetTree.D0pt, "D0pt/F");
  d0likesigntree->Branch("D0Eta", &fJetTree.D0eta, "D0eta/F");
  d0likesigntree->Branch("D0Phi", &fJetTree.D0phi, "D0phi/F");
  d0likesigntree->Branch("PionPt", &fJetTree.pionpt, "pionpt/F");
  d0likesigntree->Branch("PionEta", &fJetTree.pioneta, "pioneta/F");
  d0likesigntree->Branch("PionPhi", &fJetTree.pionphi, "pionphi/F");
  d0likesigntree->Branch("PionDCA", &fJetTree.piondca, "piondca/F");
  d0likesigntree->Branch("PionNHitsFit", &fJetTree.pionnhitsfit, "pionnhitsfit/F");
  d0likesigntree->Branch("PionNSigmaPion", &fJetTree.pionnsigmapion, "pionnsigmapion/F");
  d0likesigntree->Branch("PionNSigmaKaon", &fJetTree.pionnsigmakaon, "pionnsigmakaon/F");
  d0likesigntree->Branch("PionTofBeta", &fJetTree.piontofbeta, "piontofbeta/F");
  d0likesigntree->Branch("PionTofYLoccal", &fJetTree.ntofpionylocal, "ntofpionylocal/F");
  d0likesigntree->Branch("KaonPt", &fJetTree.kaonpt, "kaonpt/F");
  d0likesigntree->Branch("KaonEta", &fJetTree.kaoneta, "kaoneta/F");
  d0likesigntree->Branch("KaonPhi", &fJetTree.kaonphi, "kaonphi/F");
  d0likesigntree->Branch("KaonDCA", &fJetTree.kaondca, "kaondca/F");
  d0likesigntree->Branch("KaonNHitsFit", &fJetTree.kaonnhitsfit, "kaonnhitsfit/F");
  d0likesigntree->Branch("KaonNSigmaPion", &fJetTree.kaonnsigmapion, "kaonnsigmapion/F");
  d0likesigntree->Branch("KaonNSigmaKaon", &fJetTree.kaonnsigmakaon, "kaonnsigmakaon/F");
  d0likesigntree->Branch("KaonTofBeta", &fJetTree.kaontofbeta, "kaontofbeta/F");
  d0likesigntree->Branch("KaonTofYLocal", &fJetTree.ntofkaonylocal, "ntofkaonylocal/F");

  d0likesigntree->Branch("JetPt", &fJetTree.jetpt, "jetpt/F");
  d0likesigntree->Branch("JetEta", &fJetTree.jeteta, "jeteta/F");
  d0likesigntree->Branch("JetPhi", &fJetTree.jetphi, "jetphi/F");
  d0likesigntree->Branch("JetArea", &fJetTree.jetarea, "jetarea/F");
  d0likesigntree->Branch("JetE", &fJetTree.jetenergy, "jetenergy/F");
  d0likesigntree->Branch("JetNEF", &fJetTree.jetnef, "jetnef/F");
  d0likesigntree->Branch("JetNConst", &fJetTree.numberofconstituents, "numberofconstituents/I");

  // Branches to save jet data info
  d0likesigntree->Branch("TrackPt", &fJetTree.mTrackPt, "mTrackPt[numberofconstituents]/F");
  d0likesigntree->Branch("TrackEta", &fJetTree.mTrackEta, "mTrackEta[numberofconstituents]/F");
  d0likesigntree->Branch("TrackPhi", &fJetTree.mTrackPhi, "mTrackPhi[numberofconstituents]/F");
  d0likesigntree->Branch("TrackCharge", &fJetTree.mTrackCharge, "mTrackCharge[numberofconstituents]/F");
}


void StJetTreeMaker::RunJets()
{
  
  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables


  for (int jetcollection = 0; jetcollection < fJetsArr.size(); jetcollection++){

    // Jet TClonesArray
    fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it

    fJets = fJetsArr[jetcollection];

    Int_t njets = fJets->GetEntries();

    bool d0Jet = kFALSE;

    // loop over jets
    for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
      // get jet pointer
      StJet *jet = static_cast<StJet*>(fJets->At(ijet));
      if(!jet) continue;

      d0Jet = kFALSE;

      // get some jet parameters
      double jetarea = jet->Area();
      double jetpt = jet->Pt();
      double jetE = jet->E();
      double jetEta = jet->Eta();
      double jetPhi = jet->Phi();
      double jetNEF = jet->NEF();

      // get nTracks and maxTrackPt
      double NtrackConstit = jet->GetNumberOfTracks();

      for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
        int trackid = jet->TrackAt(itrk);      

        if (trackid == 99999) {
          d0Jet = kTRUE;
          break;
        }
      }

      if (d0Jet){

        fJetTree = {0};

        iTrack1 = d0TrackIndices[jetcollection][0];
        iTrack2 = d0TrackIndices[jetcollection][1];

        StPicoTrack *trk1 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack1));
        if(!trk1){ continue; }

        StPicoTrack *trk2 = static_cast<StPicoTrack*>(mPicoDst->track(iTrack2));
        if(!trk2){ continue; }

        // Event level information filler

        double mass = InvariantMass(trk1, trk2);

        TVector3 mTrk1Mom, mTrk2Mom, mResMom;
        mTrk1Mom = trk1->gMom(mVertex, Bfield);
        mTrk2Mom = trk2->gMom(mVertex, Bfield);

        mResMom = mTrk1Mom + mTrk2Mom;

        int charge1 = trk1->charge();
        int charge2 = trk2->charge();

        fJetTree.Vz = mVertex.z();
        fJetTree.VzVPD = mPicoEvent->vzVpd();
        fJetTree.Triggers = mPicoEvent->triggerIds();
        fJetTree.D0mass = mass;
        fJetTree.D0pt = mResMom.Perp();
        fJetTree.D0eta = mResMom.PseudoRapidity();
        double phi = standardPhi(mResMom.Phi());
        if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
        if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
        fJetTree.D0phi = phi;

        fJetTree.pionpt = mTrk1Mom.Perp();
        fJetTree.pioneta = mTrk1Mom.PseudoRapidity();
        fJetTree.pionphi = standardPhi(mTrk1Mom.Phi());
        fJetTree.piondca = trk1->gDCA(mPicoEvent->primaryVertex()).Mag();
        fJetTree.pionnhitsfit = trk1->nHitsFit();
        fJetTree.pionnsigmapion = trk1->nSigmaPion();
        fJetTree.pionnsigmakaon = trk1->nSigmaKaon();

        float pionbeta = GetTofBeta(trk1);
        fJetTree.piontofbeta = (isnan(pionbeta) || pionbeta < 0) ? -999 : pionbeta;
        fJetTree.ntofpionylocal = (isnan(pionbeta) || pionbeta < 0) ? -999 : mPicoDst->btofPidTraits(trk1->bTofPidTraitsIndex())->btofYLocal();

        fJetTree.kaonpt = mTrk2Mom.Perp();
        fJetTree.kaoneta = mTrk2Mom.PseudoRapidity();
        fJetTree.kaonphi = standardPhi(mTrk2Mom.Phi());
        fJetTree.kaondca = trk2->gDCA(mPicoEvent->primaryVertex()).Mag();
        fJetTree.kaonnhitsfit = trk2->nHitsFit();
        fJetTree.kaonnsigmapion = trk2->nSigmaPion();
        fJetTree.kaonnsigmakaon = trk2->nSigmaKaon();

        float kaonbeta = GetTofBeta(trk2);
        fJetTree.kaontofbeta = (isnan(kaonbeta) || kaonbeta < 0) ? -999 : kaonbeta;
        fJetTree.ntofkaonylocal = (isnan(kaonbeta) || kaonbeta < 0) ? -999 : mPicoDst->btofPidTraits(trk2->bTofPidTraitsIndex())->btofYLocal();

        // Jet level information filler

        fJetTree.jetpt = jetpt;
        fJetTree.jeteta = jetEta;
        fJetTree.jetphi = jetPhi;
        fJetTree.jetarea = jetarea;
        fJetTree.jetenergy = jetE;
        fJetTree.jetnef = jetNEF;
        fJetTree.numberofconstituents = NtrackConstit;

        for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
          int trackid = jet->TrackAt(itrk);

          if (trackid != 99999){
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

            fJetTree.mTrackPt[itrk] = pt;
            fJetTree.mTrackEta[itrk] = eta;
            fJetTree.mTrackPhi[itrk] = standardPhi(phi);
            fJetTree.mTrackCharge[itrk] = charge;
          }
        }

        if (charge1*charge2 < 0) d0unlikesigntree->Fill();
        else d0likesigntree->Fill();

      }

    }

  }

}


float StJetTreeMaker::GetTofBeta(StPicoTrack *trk){
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

bool StJetTreeMaker::IsKaon(StPicoTrack *trk){
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

bool StJetTreeMaker::IsPion(StPicoTrack *trk){
  return (abs(trk->nSigmaPion()) < 2.0);
}

Double_t StJetTreeMaker::InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2){

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

Double_t StJetTreeMaker::standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard > 2*(TMath::Pi())) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}


