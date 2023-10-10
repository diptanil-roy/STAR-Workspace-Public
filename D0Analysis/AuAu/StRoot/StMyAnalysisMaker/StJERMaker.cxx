// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StJERMaker.h"
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

ClassImp(StJERMaker)

//________________________________________________________________________
StJERMaker::StJERMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName1 = "", const char* rhoMakerName1 = "", const char* jetMakerName2 = "", const char* rhoMakerName2 = "") : StJetFrameworkPicoBase(name) //StMaker(name),
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
  fRunFlag = 0;       // see StJERMaker::fRunFlagEnum
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
  fJetMakerName = jetMakerName1;
  fRhoMakerName = rhoMakerName1;
  gJetMakerName = jetMakerName2;
  gRhoMakerName = rhoMakerName2;

  gJetMaker = 0x0;
  gJets = 0x0;
  gRhoMaker = 0x0;
  gRho = 0x0;
  gRhoVal = 0;

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

}

//
//________________________________________________________________________
StJERMaker::~StJERMaker()
{ /*  */
  // destructor
  if(hCentrality)  delete hCentrality;
  if(hMultiplicity)delete hMultiplicity;
  if(hJetPt)       delete hJetPt;
  if(hJetCorrPt)   delete hJetCorrPt;

  if(mEmcPosition) delete mEmcPosition;

  for (int i = 0; i < 7; i++){
    if (hJER[i]) delete hJER[i];
    if (hNEF[i]) delete hNEF[i];
  }

  if (hJetCorrPtvKetCorrPt) delete hJetCorrPtvKetCorrPt;
  if (hJetNEFvKetNEF) delete hJetNEFvKetNEF;
}

//
//________________________________________________________________________
Int_t StJERMaker::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms
  DeclareHistograms();

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
Int_t StJERMaker::Finish() { 
  cout << "StJERMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StJERMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StJERMaker::DeclareHistograms() {
  // binning for cent histograms
  int nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100;
  }

  // histograms
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);

  // hNEF = new TH1F("hNEF", "NEF", 100, 0, 1);

  for (int i = 0; i < 7; i++){
    hJER[i] = new TH1F(Form("hJER_%i", i), Form("hJER_%i", i), 1000, -1, 1);
    hNEF[i] = new TH1F(Form("hNEFR_%i", i), Form("hNEFR_%i", i), 1000, -1, 1);
  }


  hJetCorrPtvKetCorrPt = new TH2F("JetPtCompared", "JetPtCompared", 101, -30.5, 70.5, 101, -30.5, 70.5);
  hJetNEFvKetNEF = new TH2F("JetNEFCompared", "JetNEFCompared", 1000, -1, 1, 1000, -1, 1);

  SetSumw2();
}
//
// write histograms
//_____________________________________________________________________________
void StJERMaker::WriteHistograms() {
  // writing of histograms done here
  hCentrality->Write();
  hMultiplicity->Write();
  // hJetPt->Write();
  // hJetCorrPt->Write();
  // hNEF->Write();

  for (int i = 0; i < 7; i++){
    hJER[i]->Write();
    hNEF[i]->Write();
  }

  hJetCorrPtvKetCorrPt->Write();
  hJetNEFvKetNEF->Write();
}

//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StJERMaker::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StJERMaker::Make() {
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
  hCentrality->Fill(fCentralityScaled);
  hMultiplicity->Fill(refCorr2);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggers();

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

  fRhoVal = fRho->GetVal();

  // =========================== JetMaker =============================== //
  // get JetMaker pointer
  gJetMaker = static_cast<StJetMakerTask*>(GetMaker(gJetMakerName));
  const char *gJetMakerNameCh = gJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", gJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get jet collection associated with JetMaker
  gJets = static_cast<TClonesArray*>(gJetMaker->GetJets());
  if(!gJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn;
  }

  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
  gRhoMaker = static_cast<StRho*>(GetMaker(gRhoMakerName));
  const char *gRhoMakerNameCh = gRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", gRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  gRho = static_cast<StRhoParameter*>(gRhoMaker->GetRho());
  if(!gRho) {
    LOG_WARN << Form("Couldn't get gRho object! ") << endm;
    return kStWarn;    
  } 
  
  // get rho/area value from rho object     fRho->ls("");
  gRhoVal = gRho->GetVal();
  // =======================================================================

  // run Jets:
  RunJets();

  return kStOK;
}

//
//
//_____________________________________________________________________________________________
void StJERMaker::RunJets()
{

  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables
  Int_t ijethi = -1;
  Double_t highestjetpt = 0.0;
  Int_t njets = fJets->GetEntries();
  Int_t ojets = gJets->GetEntries();


  // Overall JER
  // JER in bins of 3-5, 5-10, 10-15, 15-20, 20-30, >30

  // Overall NEF
  // NEF in bins of 3-5, 5-10, 10-15, 15-20, 20-30, >30

  // What about jets with no matching?

  // loop over jets
  cout << "Number of Jets = " << njets << "\t" << ojets << endl;
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

    if (corrjetpt < 3) continue;

    bool hardtrack = kFALSE;

    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);      

      // get jet track pointer
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
      if (pt > 3.0) hardtrack = kTRUE;

      if (hardtrack) break;
    } // track constit loop

    double currentdeltar = 10000.;
    // double olddeltar = 10000.;

    cout << "IJET # = " << ijet << endl;

    int jetid = -99;

    for(int jjet = 0; jjet < ojets; jjet++) {  // JET LOOP
      // get jet pointer
      StJet *ket = static_cast<StJet*>(gJets->At(jjet));
      if(!ket) continue;

      // get some jet parameters
      double ketarea = ket->Area();
      double ketpt = ket->Pt();
      double korrjetpt = ket->Pt() - ketarea*gRhoVal;
      double ketE = ket->E();
      double ketEta = ket->Eta();
      double ketPhi = ket->Phi();
      double ketNEF = ket->NEF();

      double deltaR = TMath::Sqrt(pow(dEta(jetEta, ketEta), 2) + pow(dPhi(jetPhi, ketPhi), 2));

      if (deltaR > 0.4) continue;

      cout << "JJet # = " << jjet << "\t" << deltaR << endl;
      if (deltaR < currentdeltar) { currentdeltar = deltaR; jetid = jjet; }
      else {currentdeltar = currentdeltar; jetid = jetid;} 
    }

    if (jetid == -99) continue;

    StJet *ket = static_cast<StJet*>(gJets->At(jetid));
    if(!ket) continue;

    // get some jet parameters
    double ketarea = ket->Area();
    double ketpt = ket->Pt();
    double korrjetpt = ket->Pt() - ketarea*gRhoVal;
    double ketE = ket->E();
    double ketEta = ket->Eta();
    double ketPhi = ket->Phi();
    double ketNEF = ket->NEF();

    double JER = (korrjetpt - corrjetpt)/corrjetpt;
    double NEFR = (ketNEF - jetNEF)/jetNEF;

    double deltaR = TMath::Sqrt(pow(dEta(jetEta, ketEta), 2) + pow(dPhi(jetPhi, ketPhi), 2));

    cout << "JJet ID = " << jetid << "\t" << deltaR << "\t" << corrjetpt << "\t" << korrjetpt << "\t" << jetEta << "\t" << ketEta << "\t" << jetPhi << "\t" << ketPhi << endl;

    if (abs(JER) > 0.8) cout << "************************* CHECK THIS JET ***********************" << endl;


    // cout << "NEFR = " << jetNEF << "\t" << ketNEF << "\t" << NEFR << endl;

    hJetCorrPtvKetCorrPt->Fill(corrjetpt, korrjetpt);
    hJetNEFvKetNEF->Fill(jetNEF, ketNEF);

    if (corrjetpt > 3 && corrjetpt < 5){hJER[0]->Fill(JER); hNEF[0]->Fill(NEFR);}
    if (corrjetpt > 5 && corrjetpt < 10){hJER[1]->Fill(JER); hNEF[1]->Fill(NEFR);}
    if (corrjetpt > 10 && corrjetpt < 15){hJER[2]->Fill(JER); hNEF[2]->Fill(NEFR);}
    if (corrjetpt > 15 && corrjetpt < 20){hJER[3]->Fill(JER); hNEF[3]->Fill(NEFR);}
    if (corrjetpt > 20 && corrjetpt < 30){hJER[4]->Fill(JER); hNEF[4]->Fill(NEFR);}
    if (corrjetpt > 30 && corrjetpt < 50){hJER[5]->Fill(JER); hNEF[5]->Fill(NEFR);}

    if (corrjetpt > 3) hJER[6]->Fill(JER); hNEF[6]->Fill(NEFR);

  }

}

// //
// //
// //________________________________________________________________________
// void StJERMaker::RunTracks()
// {
//   // constants: assume neutral pion mass
//   double pi = 1.0*TMath::Pi();
//   double pi0mass = Pico::mMass[0]; // GeV
//   unsigned int ntracks = mPicoDst->numberOfTracks();

//   // loop over ALL tracks in PicoDst 
//   for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
//     // get track pointer
//     StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
//     if(!trk){ continue; }

//     // acceptance and kinematic quality cuts
//     if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

//     // primary track switch: get momentum vector of track - global or primary track
//     TVector3 mTrkMom;
//     if(doUsePrimTracks) {
//       // get primary track vector
//       mTrkMom = trk->pMom();
//     } else {
//       // get global track vector
//       mTrkMom = trk->gMom(mVertex, Bfield);
//     }

//     // track variables
//     double pt = mTrkMom.Perp();
//     double phi = mTrkMom.Phi();
//     if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
//     if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi
//     double eta = mTrkMom.PseudoRapidity();
//     double px = mTrkMom.x();
//     double py = mTrkMom.y();
//     double pz = mTrkMom.z();
//     short charge = trk->charge();

//     // do some things with tracks here
//     // ......

//     // fill track histograms here
//     // .........

//   } // track loop

// }  // track function

// //
// //
// //________________________________________________________________________
// void StJERMaker::RunTowers()
// {
//   // constants: assume neutral pion mass
//   double pi = 1.0*TMath::Pi();
//   double pi0mass = Pico::mMass[0]; // GeV

//   // looping over clusters - STAR: matching already done, get # of clusters and set variables
//   unsigned int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();

//   // loop over ALL clusters in PicoDst
//   for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
//     StPicoBEmcPidTraits *cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
//     if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

//     // cluster and tower ID
//     int clusID = cluster->bemcId();  // index in bemc point array
//     int towID = cluster->btowId();   // projected tower Id: 1 - 4800
//     int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
//     int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
//     if(towID < 0) continue;

//     // cluster and tower position - from vertex and ID
//     TVector3 towPosition = mEmcPosition->getPosFromVertex(mVertex, towID);
//     double towPhi = towPosition.Phi();
//     double towEta = towPosition.PseudoRapidity();

//     // matched track index
//     int trackIndex = cluster->trackIndex();
//     // get track pointer
//     StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
//     if(!trk) { cout<<"No trk pointer...."<<endl; continue; }
//     if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

//   } // BEmc loop

//   // loop over towers
//   int nTowers = mPicoDst->numberOfBTowHits();
//   for(int itow = 0; itow < nTowers; itow++) {
//     // get tower pointer
//     StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
//     if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

//     // tower ID: tower ID shifted by +1 from array index
//     int towerID = itow + 1;
//     if(towerID < 0) continue; // double check these aren't still in the event list

//     // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
//     TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
//     double towerPhi = towerPosition.Phi();
//     double towerEta = towerPosition.PseudoRapidity();
//     int towerADC = tower->adc();
//     double towerE = tower->energy();
//     double towerEt = towerE / (1.0*TMath::CosH(towerEta)); // THIS should be USED

//     // do stuff with towers and fill histograms here
//     // ........


//   } // tower loop

// }  // run towers function

//
//
// __________________________________________________________________________________
void StJERMaker::SetSumw2() {

  TH1::SetDefaultSumw2();
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  hJetPt->Sumw2();
  hJetCorrPt->Sumw2();
}

//
// Function: get relative phi of jet and track (-0.5pi, 1.5pi)
//________________________________________________________________________
Double_t StJERMaker::RelativePhi(Double_t mphi,Double_t vphi) const
{ // function to calculate relative PHI
  double dphi = mphi-vphi;

  // set dphi to operate on adjusted scale
  if(dphi < -0.5*pi) dphi += 2.0*TMath::Pi();
  if(dphi >  1.5*pi) dphi -= 2.0*TMath::Pi();

  // test check
  if( dphi < -0.5*pi || dphi > 1.5*pi )
    Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName());

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}

//
// Function: calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
//_________________________________________________________________________
Double_t StJERMaker::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{
  Double_t pi = 1.0*TMath::Pi();
  Double_t dphi = 1.0*TMath::Abs(EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi < -pi ){
    dphi = dphi + pi;
  } // this assumes we are doing full jets currently 
 
  if(dphi > 1.5*pi) dphi -= 2.0*pi;
  if((dphi > 1.0*pi) && (dphi < 1.5*pi)) dphi -= 1.0*pi;
  if((dphi > 0.5*pi) && (dphi < 1.0*pi)) dphi -= 1.0*pi;
  dphi = 1.0*TMath::Abs(dphi);

  // test check
  if( dphi < 0.0 || dphi > 0.5*pi ) {
    //Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());
    cout<<"dPhi not in range [0, 0.5*Pi]!"<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}

//
//
//_________________________________________________________________________
void StJERMaker::FillEmcTriggers() {
  // number of Emcal Triggers
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // fill for valid triggers
    if(emcTrig->isHT0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
}
//
// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StJERMaker::DoComparison(int myarr[], int elems) {
  //std::cout << "Length of array = " << (sizeof(myarr)/sizeof(*myarr)) << std::endl;
  bool match = kFALSE;

  // loop over specific physics selection array and compare to specific event trigger
  for(int i=0; i<elems; i++) {
    if(mPicoEvent->isTrigger(myarr[i])) match = kTRUE;
    if(match) break;
  }

  return match;
}
