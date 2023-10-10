#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoPrescales/StPicoPrescales.h"
#include "StPicoMixedEventMaker.h"
#include "StPicoEventMixer.h"
#include "StD0Hists.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StEventPlane/StEventPlane.h"

#include <vector>

ClassImp(StPicoMixedEventMaker)
// _________________________________________________________
StPicoMixedEventMaker::StPicoMixedEventMaker(char const* name, StPicoDstMaker* picoMaker, StRefMultCorr* grefmultCorrUtil, StEventPlane* eventPlaneMaker,
					     char const* outputBaseFileName,  char const* inputPicoList, int harmonic) :
   StMaker(name), mPicoDstMaker(picoMaker),  mPicoEvent(NULL),
   mGRefMultCorrUtil(grefmultCorrUtil), mEventPlaneMaker(eventPlaneMaker),
   mFailedRunnumber(0), mOuputFileBaseName(outputBaseFileName), 
  mInputFileName(inputPicoList), mEventCounter(0), mHarmonic(harmonic)
{
   mGRefMultCorrUtil->print();
   for (int iVz = 0 ; iVz < 10 ; ++iVz)
   {
      for (int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality)
      {
         for (int iEventPlane = 0 ; iEventPlane < 10 ; ++iEventPlane)
         {
            mPicoEventMixer[iVz][iCentrality][iEventPlane] = NULL;
         }
      }
   }

   //   TH1::AddDirectory(false);

   mOutputFile = new TFile(Form("%s.d0Hists.root", mOuputFileBaseName.Data()), "RECREATE");

   mEventPlaneMaker->setFileOut(mOutputFile);
   mD0Hists = new StD0Hists("picoME",mHarmonic);

   // -- constructor
}

// _________________________________________________________
StPicoMixedEventMaker::~StPicoMixedEventMaker()
{
   delete mGRefMultCorrUtil;
   for (int iVz = 0 ; iVz < 10 ; ++iVz)
   {
      for (int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality)
      {
         for (int iEventPlane = 0 ; iEventPlane < 10 ; ++iEventPlane)
         {
            delete mPicoEventMixer[iVz][iCentrality][iEventPlane];
         }
      }
   }
}

// Method should load Q vector stuff from Hao, needs fixing
// _________________________________________________________
bool StPicoMixedEventMaker::loadEventPlaneCorr(StEventPlane const * mEventPlane)
{
   //needs to implement, will currently break maker
   if (!mEventPlane)
   {
      LOG_WARN << "No EventPlane ! Skipping! " << endm;
      return kFALSE;
   }
   if (!mEventPlane->getAcceptEvent())
   {
      // LOG_WARN << "StPicoMixedEvent::THistograms and TProiles NOT found! shoudl check the input Qvector files From HaoQiu ! Skipping this run! " << endm;
      return kFALSE;
   }
   return kTRUE;
}

// _________________________________________________________
Int_t StPicoMixedEventMaker::Init()
{
   for (int iVz = 0 ; iVz < 10 ; ++iVz)
   {
      for (int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality)
      {
         for (int iEventPlane = 0 ; iEventPlane < 10 ; ++iEventPlane)
         {
            mPicoEventMixer[iVz][iCentrality][iEventPlane] = new StPicoEventMixer(iCentrality, iVz, iEventPlane, mEventPlaneMaker, mD0Hists);
            //         mPicoEventMixer[iVz][iCentrality][iEventPlane]->setEventBuffer(10);
         }
      }
   }
   // mGRefMultCorrUtil = new StRefMultCorr("grefmult");
   // if(!LoadEventPlaneCorr(mRunId)){
   // LOG_WARN << "Event plane calculations unavalable! Skipping"<<endm;
   // return kStOk;
   // }

   // -- reset event to be in a defined state
   //resetEvent();
   return kStOK;
}

// _________________________________________________________
Int_t StPicoMixedEventMaker::Finish()
{
   // -- Inhertited from StMaker
   //    NOT TO BE OVERWRITTEN by daughter class
   //    daughter class should implement FinishHF()
   for (int iVz = 0 ; iVz < 10 ; ++iVz)
   {
      for (int iCentrality = 0 ; iCentrality < 9 ; ++iCentrality)
      {
         for (int iEventPlane = 0 ; iEventPlane < 10 ; ++iEventPlane)
         {
            mPicoEventMixer[iVz][iCentrality][iEventPlane]->finish();
            //delete mPicoEventMixer[iVz][iCentrality];
         }
      }
   }
   cout << "StPicoMixedEventMaker::Finish()" << endl;
   mOutputFile->cd();
   for(int ii = 0 ; ii<mxeCuts::nCutsSets; ++ii){

#ifdef __run_w_DaugHisto__     
     //Daught
     mD0Hists->hD0CentPtEtaMDphiDaug[ii]->Write();
     mD0Hists->hD0CentPtEtaMDphiDaugLikeSign[ii]->Write();
     mD0Hists->hD0CentPtEtaMDphiDaugMixed[ii]->Write();
     mD0Hists->hD0CentPtEtaMDphiDaugLikeSignMixed[ii]->Write();
#endif

     //Eta sub
     // mD0Hists->hD0EtaSubCentPtMDphi[ii]->Write();
     // mD0Hists->hD0EtaSubCentPtMDphiLikeSign[ii]->Write();
     // mD0Hists->hD0EtaSubCentPtMDphiMixed[ii]->Write();
     // mD0Hists->hD0EtaSubCentPtMDphiLikeSignMixed[ii]->Write();
     
     // Eta Gap
     // mD0Hists->hD0CentPtMDphiEtaGap[ii]->Write();
     // mD0Hists->hD0CentPtMDphiEtaGapLikeSign[ii]->Write();
     // mD0Hists->hD0CentPtMDphiEtaGapMixed[ii]->Write();
     // mD0Hists->hD0CentPtMDphiEtaGapLikeSignMixed[ii]->Write();
   }
#ifdef __run_w_QA__
   // QA plots
   mD0Hists->mSE_US_DecayTopology->Write();
   mD0Hists->mSE_US_PointingAngle->Write();
   mD0Hists->mSE_US_DecayL->Write();
   mD0Hists->mSE_US_Dca12->Write();
   mD0Hists->mSE_US_PionDca2Vtx->Write();
   mD0Hists->mSE_US_KaonDca2Vtx->Write();
   mD0Hists->mSE_US_D0Dca2Vtx->Write();

   //
   mD0Hists->mSE_LS_PointingAngle->Write();
   mD0Hists->mSE_LS_DecayL->Write();
   mD0Hists->mSE_LS_Dca12->Write();
   mD0Hists->mSE_LS_PionDca2Vtx->Write();
   mD0Hists->mSE_LS_KaonDca2Vtx->Write();
   mD0Hists->mSE_LS_D0Dca2Vtx->Write();
   //
   mD0Hists->mME_US_DecayTopology->Write();
   mD0Hists->mME_US_PointingAngle->Write();
   mD0Hists->mME_US_DecayL->Write();
   mD0Hists->mME_US_Dca12->Write();
   mD0Hists->mME_US_PionDca2Vtx->Write();
   mD0Hists->mME_US_KaonDca2Vtx->Write();
   mD0Hists->mME_US_D0Dca2Vtx->Write();
#endif
   mOutputFile->Write();
   //   mOutputFile->Close();
   cout << "StPicoMixedEventMaker::Finish() done" << endl;

   return kStOK;
}
// _________________________________________________________
void StPicoMixedEventMaker::Clear(Option_t* opt)
{
}
// _________________________________________________________
Int_t StPicoMixedEventMaker::Make()
{

   if (!mPicoDstMaker)
   {
      LOG_WARN << "No PicoDstMaker! Skipping! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();
   if (!picoDst)
   {
      LOG_WARN << "No picoDst ! Skipping! " << endm;
      return kStWarn;
   }

   // - GRef from Guannan
   if (!mGRefMultCorrUtil)
   {
      LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
      return kStWarn;
   }

   //Load event
   mPicoEvent = picoDst->event();
   if (!mPicoEvent)
   {
      cerr << "Error opening picoDst Event, skip!" << endl;
      return kStWarn;
   }


   mD0Hists->hTotalNumberOfEvents->Fill(0);

   for (int i = 0; i < 32; i++)
     if ( mPicoEvent->isTrigger(i) )
         mD0Hists->hTrigger->Fill(i);

   if(!isMinBiasTrigger()) return kStOk;
   StThreeVectorF pVtx = mPicoEvent->primaryVertex();
   //Remove bad vertices
   mD0Hists->hVzVpdVz->Fill(pVtx.z(), mPicoEvent->vzVpd());
   mD0Hists->hVzDiff->Fill(mPicoEvent->vzVpd() - pVtx.z());
   mD0Hists->hVxy->Fill(pVtx.x(), pVtx.y());
   
   if(isGoodEvent(pVtx))
   {
     mD0Hists->hRefMult->Fill(mPicoEvent->refMult());
     mD0Hists->hGRefMult->Fill(mPicoEvent->grefMult());

     // - GRef from Guannan
     if (!mGRefMultCorrUtil)
     {
       LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
       return kStWarn;
     }

     mGRefMultCorrUtil->init(mPicoEvent->runId());
     mGRefMultCorrUtil->initEvent(mPicoEvent->grefMult(), pVtx.z(), mPicoEvent->ZDCx()) ;
     int const centrality  = mGRefMultCorrUtil->getCentralityBin9();
     float weight = mGRefMultCorrUtil->getWeight();
     mD0Hists->hCentrality->Fill(centrality);
     mD0Hists->hCentralityWeighted->Fill(centrality, weight);
     if (centrality < 0 || centrality > 8) return kStOk;
     // if (centrality < -0.5 || centrality > 1.5) return kStOk;

     int const vz_bin = (int)((6 + pVtx.z()) / 1.2) ;
     if (vz_bin < 0  ||  vz_bin > 9) return kStOk;


     if (mFailedRunnumber != mPicoEvent->runId())
     {
       if (!loadEventPlaneCorr(mEventPlaneMaker))
       {
         LOG_WARN << "Event plane calculations unavalable! Skipping" << endm;
         mFailedRunnumber = picoDst->event()->runId();
         return kStOK;
       }
     }
     else  return kStOK;


     float const eventPlane = mEventPlaneMaker->getEventPlane();
     int const eventPlane_bin = (int)(eventPlane / TMath::Pi() * 10.) ;
     if (eventPlane_bin < 0  ||  eventPlane_bin > 9 || mEventPlaneMaker->eventPlaneStatus()) return kStOk;

     mD0Hists->hCentVzPsi->Fill(centrality, pVtx.z(), eventPlane, weight);

     if (mPicoEventMixer[vz_bin][centrality][eventPlane_bin]->addPicoEvent(picoDst, pVtx, weight))
       mPicoEventMixer[vz_bin][centrality][eventPlane_bin]->mixEvents();
   }

   return kStOk;
}

bool StPicoMixedEventMaker::isGoodEvent(StThreeVectorF const& pVtx) const
{
  return !(fabs(pVtx.x()) < mxeCuts::Verror && fabs(pVtx.y()) < mxeCuts::Verror && fabs(pVtx.z()) < mxeCuts::Verror) &&
    fabs(pVtx.z()) < mxeCuts::maxVz && fabs(pVtx.z() - mPicoEvent->vzVpd()) < mxeCuts::vzVpdVz &&
    sqrt(pow(pVtx.x(), 2) + pow(pVtx.y(), 2)) < mxeCuts::Vrcut;
}
bool StPicoMixedEventMaker::isMinBiasTrigger() const 
{
  for(int ii = 0; ii<mxeCuts::nTrig; ++ii){
    if( mPicoEvent->isTrigger(mxeCuts::mTriggerId[ii]) )
      return true;
  }
  return false;
}
