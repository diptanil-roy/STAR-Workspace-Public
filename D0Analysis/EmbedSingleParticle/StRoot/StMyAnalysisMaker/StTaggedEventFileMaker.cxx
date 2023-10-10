// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StTaggedEventFileMaker.h"
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
#include "TTree.h"
#include "TChain.h"

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

// D0 Includes
#include "StTagD0Events.h"


ClassImp(StTaggedEventFileMaker)

//________________________________________________________________________
StTaggedEventFileMaker::StTaggedEventFileMaker(const char* name, StPicoDstMaker *picoMaker, const char* D0TaggerName = "", const char* outName="") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StTaggedEventFileMaker::fRunFlagEnum
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
  fD0MakerName = D0TaggerName;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

}

//
//________________________________________________________________________
StTaggedEventFileMaker::~StTaggedEventFileMaker()
{ /*  */
  // destructor
  // if(hCentrality)  delete hCentrality;
  // if(hMultiplicity)delete hMultiplicity;
  // if(hJetPt)       delete hJetPt;
  // if(hJetCorrPt)   delete hJetCorrPt;

  if (outtree) delete outtree;
  if(mEmcPosition) delete mEmcPosition;
}

//
//________________________________________________________________________
Int_t StTaggedEventFileMaker::Init() {
  StJetFrameworkPicoBase::Init();

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  outtree = mPicoDstMaker->chain()->CloneTree(0);
  // intree = new TTree();
  // outtree = new TTree();
  // declare histograms
  // DeclareHistograms();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  // fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StTaggedEventFileMaker::Finish() { 
  cout << "StTaggedEventFileMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    outtree->Write();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StTaggedEventFileMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StTaggedEventFileMaker::DeclareHistograms() {
  // binning for cent histograms
  // int nHistCentBins = 20;

  // // binning for mult histograms
  // double kHistMultMax = 800.;
  // int kHistMultBins = 400;

  // // pp specific settings
  // if(doppAnalysis) {
  //   kHistMultMax = 100.;
  //   kHistMultBins = 100.;
  // }

  // // histograms
  // hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
  // hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  // // jet QA histos
  // hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  // hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);
}
//
// write histograms
//_____________________________________________________________________________
void StTaggedEventFileMaker::WriteHistograms() {
  // writing of histograms done here
  // hCentrality->Write();
  // hMultiplicity->Write();
  // hJetPt->Write();
  // hJetCorrPt->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StTaggedEventFileMaker::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StTaggedEventFileMaker::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;



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


  mD0Tagger = static_cast<StTagD0Events*>(GetMaker(fD0MakerName));
  if(!mD0Tagger) {
    LOG_WARN << " No D0Tagger! Skip! " << endm;
    return kStWarn;
  }

  // cout << mPicoDstMaker->chain()->GetEntries() << endl;

  // intree = mPicoDstMaker->tree();

  // if (!intree) cout << "YELL" << endl;
  // cout << mPicoDstMaker->chain()->GetEntries() << endl;

  if (!mD0Tagger->DoesEventHaveD0() && !mD0Tagger->DoesEventHaveD0Bg() 
   && !mD0Tagger->DoesEventHaveTightD0() && !mD0Tagger->DoesEventHaveTightD0Bg() 
   && !mD0Tagger->DoesEventHaveLooseD0() && !mD0Tagger->DoesEventHaveLooseD0Bg()) return kStOK;

  // cout << "Made it this far." << endl;

  outtree->Fill();


  return kStOK;
}
