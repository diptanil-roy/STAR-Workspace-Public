#ifndef StBEMCIssueCorrector_h
#define StBEMCIssueCorrector_h

#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

// ROOT classes
// ROOT classes
class TTree;
class TBranch;
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class THnSparse;
class TProfile;
class TString;
class TVector3;
class TTree;
class TChain;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

// jet-framework classes
class StCentMaker;
class StEmcPosition2;
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;

//D0 Includes
class StTagD0Events;


// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t gMaxEvent = 1;
const Int_t gMaxTrack = 5000;
const Int_t gMaxEmcTrigger = 5000;
const Int_t gMaxMtdTrigger = 5000;
const Int_t gMaxBTowHit = 5000;
const Int_t gMaxBTofHit = 5000;
const Int_t gMaxMtdHit = 5000;
const Int_t gMaxBbcHit = 5000;
const Int_t gMaxEpdHit = 5000;
const Int_t gMaxFmsHit = 5000;
const Int_t gMaxEmcPidTraits = 5000;
const Int_t gMaxBTofPidTraits = 5000;
const Int_t gMaxMtdPidTraits = 5000;
const Int_t gMaxTrackCovMatrix = 5000;

class StBEMCIssueCorrector : public StJetFrameworkPicoBase {
  public:

    StBEMCIssueCorrector(const char* name, StPicoDstMaker *picoMaker, const char* filelistname, const char* outName);
    virtual ~StBEMCIssueCorrector();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // const char *mFilename;

    TString mfilelist;

    TChain *fVariablesAssignedToTree;

    // Declaration of leaf types
    Int_t           Event_;
    Int_t           Event_mRunId[gMaxEvent];   //[Event_]
    Int_t           Event_mEventId[gMaxEvent];   //[Event_]
    UShort_t        Event_mFillId[gMaxEvent];   //[Event_]
    Float_t         Event_mBField[gMaxEvent];   //[Event_]
    Int_t           Event_mTime[gMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexX[gMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexY[gMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexZ[gMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexErrorX[gMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexErrorY[gMaxEvent];   //[Event_]
    Float_t         Event_mPrimaryVertexErrorZ[gMaxEvent];   //[Event_]
    Float_t         Event_mRanking[gMaxEvent];   //[Event_]
    UShort_t        Event_mNBEMCMatch[gMaxEvent];   //[Event_]
    UShort_t        Event_mNBTOFMatch[gMaxEvent];   //[Event_]
    vector<unsigned int> Event_mTriggerIds[gMaxEvent];
    UShort_t        Event_mRefMultFtpcEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultFtpcWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultNeg[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultPos[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2NegEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2PosEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2NegWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult2PosWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3NegEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3PosEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3NegWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult3PosWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4NegEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4PosEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4NegWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMult4PosWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfNegEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfPosEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfNegWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mRefMultHalfPosWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mGRefMult[gMaxEvent];   //[Event_]
    UShort_t        Event_mNumberOfGlobalTracks[gMaxEvent];   //[Event_]
    UShort_t        Event_mbTofTrayMultiplicity[gMaxEvent];   //[Event_]
    UShort_t        Event_mNHitsHFT[gMaxEvent][4];   //[Event_]
    UChar_t         Event_mNVpdHitsEast[gMaxEvent];   //[Event_]
    UChar_t         Event_mNVpdHitsWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mNTofT0[gMaxEvent];   //[Event_]
    Float_t         Event_mVzVpd[gMaxEvent];   //[Event_]
    UInt_t          Event_mZDCx[gMaxEvent];   //[Event_]
    UInt_t          Event_mBBCx[gMaxEvent];   //[Event_]
    Float_t         Event_mBackgroundRate[gMaxEvent];   //[Event_]
    Float_t         Event_mBbcBlueBackgroundRate[gMaxEvent];   //[Event_]
    Float_t         Event_mBbcYellowBackgroundRate[gMaxEvent];   //[Event_]
    Float_t         Event_mBbcEastRate[gMaxEvent];   //[Event_]
    Float_t         Event_mBbcWestRate[gMaxEvent];   //[Event_]
    Float_t         Event_mZdcEastRate[gMaxEvent];   //[Event_]
    Float_t         Event_mZdcWestRate[gMaxEvent];   //[Event_]
    UShort_t        Event_mZdcSumAdcEast[gMaxEvent];   //[Event_]
    UShort_t        Event_mZdcSumAdcWest[gMaxEvent];   //[Event_]
    UShort_t        Event_mZdcSmdEastHorizontal[gMaxEvent][8];   //[Event_]
    UShort_t        Event_mZdcSmdEastVertical[gMaxEvent][8];   //[Event_]
    UShort_t        Event_mZdcSmdWestHorizontal[gMaxEvent][8];   //[Event_]
    UShort_t        Event_mZdcSmdWestVertical[gMaxEvent][8];   //[Event_]
    UShort_t        Event_mBbcAdcEast[gMaxEvent][24];   //[Event_]
    UShort_t        Event_mBbcAdcWest[gMaxEvent][24];   //[Event_]
    UChar_t         Event_mHighTowerThreshold[gMaxEvent][4];   //[Event_]
    UChar_t         Event_mJetPatchThreshold[gMaxEvent][4];   //[Event_]
    Int_t           Track_;
    UShort_t        Track_mId[gMaxTrack];   //[Track_]
    UShort_t        Track_mChi2[gMaxTrack];   //[Track_]
    Float_t         Track_mPMomentumX[gMaxTrack];   //[Track_]
    Float_t         Track_mPMomentumY[gMaxTrack];   //[Track_]
    Float_t         Track_mPMomentumZ[gMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumX[gMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumY[gMaxTrack];   //[Track_]
    Float_t         Track_mGMomentumZ[gMaxTrack];   //[Track_]
    Float_t         Track_mOriginX[gMaxTrack];   //[Track_]
    Float_t         Track_mOriginY[gMaxTrack];   //[Track_]
    Float_t         Track_mOriginZ[gMaxTrack];   //[Track_]
    Float16_t       Track_mDedx[gMaxTrack];   //[Track_]
    Float16_t       Track_mDedxError[gMaxTrack];   //[Track_]
    Char_t          Track_mNHitsFit[gMaxTrack];   //[Track_]
    UChar_t         Track_mNHitsMax[gMaxTrack];   //[Track_]
    UChar_t         Track_mNHitsDedx[gMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaPion[gMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaKaon[gMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaProton[gMaxTrack];   //[Track_]
    Short_t         Track_mNSigmaElectron[gMaxTrack];   //[Track_]
    UInt_t          Track_mTopologyMap[gMaxTrack][2];   //[Track_]
    Short_t         Track_mBEmcPidTraitsIndex[gMaxTrack];   //[Track_]
    Short_t         Track_mBTofPidTraitsIndex[gMaxTrack];   //[Track_]
    Short_t         Track_mBEmcMatchedTowerIndex[gMaxTrack];   //[Track_]

    Short_t         copiedTrack_mBEmcMatchedTowerIndex[gMaxTrack];   //[Track_]

    Short_t         Track_mMtdPidTraitsIndex[gMaxTrack];   //[Track_]
    Int_t           EmcTrigger_;
    UChar_t         EmcTrigger_mFlag[gMaxEmcTrigger];   //[EmcTrigger_]
    UShort_t        EmcTrigger_mId[gMaxEmcTrigger];   //[EmcTrigger_]
    UShort_t        EmcTrigger_mAdc[gMaxEmcTrigger];   //[EmcTrigger_]
    Int_t           MtdTrigger_;
    UShort_t        MtdTrigger_mVpdTacSum[gMaxMtdTrigger];   //[MtdTrigger_]
    UInt_t          MtdTrigger_mTHUBtime[gMaxMtdTrigger][2];   //[MtdTrigger_]
    UShort_t        MtdTrigger_mQTtacSum[gMaxMtdTrigger][8][8];   //[MtdTrigger_]
    UShort_t        MtdTrigger_mMT101Tac[gMaxMtdTrigger][8][2];   //[MtdTrigger_]
    UChar_t         MtdTrigger_mMT101Id[gMaxMtdTrigger][8][2];   //[MtdTrigger_]
    UInt_t          MtdTrigger_mTF201TriggerBit[gMaxMtdTrigger];   //[MtdTrigger_]
    Char_t          MtdTrigger_mShouldHaveRejectEvent[gMaxMtdTrigger];   //[MtdTrigger_]
    Int_t           BTowHit_;
    UShort_t        BTowHit_mAdc[gMaxBTowHit];   //[BTowHit_]
    Short_t         BTowHit_mE[gMaxBTowHit];   //[BTowHit_]
    Int_t           BTofHit_;
    Short_t         BTofHit_mId[gMaxBTofHit];   //[BTofHit_]
    Int_t           MtdHit_;
    Short_t         MtdHit_mgChannel[gMaxMtdHit];   //[MtdHit_]
    UChar_t         MtdHit_mTriggerFlag[gMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mLeadingEdgeTime_first[gMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mLeadingEdgeTime_second[gMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mTrailingEdgeTime_first[gMaxMtdHit];   //[MtdHit_]
    Float_t         MtdHit_mTrailingEdgeTime_second[gMaxMtdHit];   //[MtdHit_]
    Int_t           BbcHit_;
    Short_t         BbcHit_mId[gMaxBbcHit];   //[BbcHit_]
    Int_t           BbcHit_mQTdata[gMaxBbcHit];   //[BbcHit_]
    Int_t           EpdHit_;
    Short_t         EpdHit_mId[gMaxEpdHit];   //[EpdHit_]
    Int_t           EpdHit_mQTdata[gMaxEpdHit];   //[EpdHit_]
    Float_t         EpdHit_mnMIP[gMaxEpdHit];   //[EpdHit_]
    Int_t           FmsHit_;
    UShort_t        FmsHit_mChannelDetectorId[gMaxFmsHit];   //[FmsHit_]
    UShort_t        FmsHit_mAdc[gMaxFmsHit];   //[FmsHit_]
    Int_t           EmcPidTraits_;
    Short_t         EmcPidTraits_mTrackIndex[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcId[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcAdc0[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcE0[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcE[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcZDist[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBemcPhiDist[gMaxEmcPidTraits];   //[EmcPidTraits_]
    UChar_t         EmcPidTraits_mBemcSmdNEta[gMaxEmcPidTraits];   //[EmcPidTraits_]
    UChar_t         EmcPidTraits_mBemcSmdNPhi[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowId[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Char_t          EmcPidTraits_mBtowId23[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowE[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowE2[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowE3[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowEtaDist[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Short_t         EmcPidTraits_mBtowPhiDist[gMaxEmcPidTraits];   //[EmcPidTraits_]
    Int_t           BTofPidTraits_;
    Short_t         BTofPidTraits_mTrackIndex[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofCellId[gMaxBTofPidTraits];   //[BTofPidTraits_]
    UChar_t         BTofPidTraits_mBTofMatchFlag[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Float_t         BTofPidTraits_mBTof[gMaxBTofPidTraits];   //[BTofPidTraits_]
    UShort_t        BTofPidTraits_mBTofBeta[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofYLocal[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofZLocal[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofHitPosX[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofHitPosY[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Short_t         BTofPidTraits_mBTofHitPosZ[gMaxBTofPidTraits];   //[BTofPidTraits_]
    Int_t           MtdPidTraits_;
    Short_t         MtdPidTraits_mTrackIndex[gMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mMtdHitIndex[gMaxMtdPidTraits];   //[MtdPidTraits_]
    Char_t          MtdPidTraits_mMatchFlag[gMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mDeltaY[gMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mDeltaZ[gMaxMtdPidTraits];   //[MtdPidTraits_]
    Float_t         MtdPidTraits_mDeltaTimeOfFlight[gMaxMtdPidTraits];   //[MtdPidTraits_]
    UShort_t        MtdPidTraits_mBeta[gMaxMtdPidTraits];   //[MtdPidTraits_]
    Short_t         MtdPidTraits_mMtdHitChan[gMaxMtdPidTraits];   //[MtdPidTraits_]
    Int_t           TrackCovMatrix_;
    Float16_t       TrackCovMatrix_mImp[gMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mZ[gMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mPsi[gMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mPti[gMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mTan[gMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mCurv[gMaxTrackCovMatrix];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mSigma[gMaxTrackCovMatrix][5];   //[TrackCovMatrix_]
    Float16_t       TrackCovMatrix_mCorr[gMaxTrackCovMatrix][10];   //[TrackCovMatrix_]

  protected:
    
  private:
    // variables
    Int_t                   fRunNumber;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);

    // position object
    StEmcPosition2         *mEmcPosition;

    // D0 tagger objects
    StTagD0Events         *mD0Tagger;               // D0 Tagger object

    // histos
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
 
    // jet histos
    TH1F *hJetPt;//!
    TH1F *hJetCorrPt;//!

    // Trees
    TTree *intree;
    TTree *outtree;

    // bad and dead tower list
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list
    std::set<Int_t>        badRuns;

    // base class pointer object
    StJetFrameworkPicoBase *mBaseMaker;

    // maker names
    TString                fD0MakerName;

    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;


    ClassDef(StBEMCIssueCorrector, 1)
};
#endif
