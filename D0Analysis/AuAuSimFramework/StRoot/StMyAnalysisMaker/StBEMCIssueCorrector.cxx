// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StBEMCIssueCorrector.h"
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


ClassImp(StBEMCIssueCorrector)

//________________________________________________________________________
StBEMCIssueCorrector::StBEMCIssueCorrector(const char* name, StPicoDstMaker *picoMaker, const char* filelistname = "", const char* outName="") : StJetFrameworkPicoBase(name) //StMaker(name),
{

  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  
  mfilelist = filelistname;

}

//
//________________________________________________________________________
StBEMCIssueCorrector::~StBEMCIssueCorrector()
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
Int_t StBEMCIssueCorrector::Init() {
  StJetFrameworkPicoBase::Init();

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  outtree = mPicoDstMaker->chain()->CloneTree(0);
  outtree->SetBranchAddress("Track.mBEmcMatchedTowerIndex", copiedTrack_mBEmcMatchedTowerIndex);

  fVariablesAssignedToTree = new TChain("PicoDst");
  fVariablesAssignedToTree->Add(mfilelist.Data());

  fVariablesAssignedToTree->SetBranchStatus("*",0);
  fVariablesAssignedToTree->SetBranchStatus("Event", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mRunId", 1);
  fVariablesAssignedToTree->SetBranchStatus("Event.mEventId", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track", 1);
  fVariablesAssignedToTree->SetBranchStatus("Track.mId", 1);

  fVariablesAssignedToTree->SetBranchAddress("Event", &Event_);
  fVariablesAssignedToTree->SetBranchAddress("Event.mRunId", Event_mRunId);
  fVariablesAssignedToTree->SetBranchAddress("Event.mEventId", Event_mEventId);
  fVariablesAssignedToTree->SetBranchAddress("Track", &Track_);
  fVariablesAssignedToTree->SetBranchAddress("Track.mId", Track_mId);
  fVariablesAssignedToTree->SetBranchAddress("Track.mBEmcMatchedTowerIndex", Track_mBEmcMatchedTowerIndex);

  fVariablesAssignedToTree->BuildIndex("mRunId", "mEventId");

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StBEMCIssueCorrector::Finish() { 
  cout << "StBEMCIssueCorrector::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    outtree->Write();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StBEMCIssueCorrector::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StBEMCIssueCorrector::Clear(Option_t *opt) {
  // fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StBEMCIssueCorrector::Make() {
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

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();

  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  fVariablesAssignedToTree->GetEntryWithIndex( mPicoEvent->runId(),  mPicoEvent->eventId());

  const Int_t ntracks = mPicoDst->numberOfTracks();

  for(unsigned short itrk = 0; itrk < ntracks; itrk++){
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrk));
    if(!trk){ continue; }

    cout << trk->bemcTowerIndex() << "\t" << Track_mBEmcMatchedTowerIndex[itrk] << endl;

    copiedTrack_mBEmcMatchedTowerIndex[itrk] = Track_mBEmcMatchedTowerIndex[itrk];
  }

  outtree->Fill();


  return kStOK;
}
