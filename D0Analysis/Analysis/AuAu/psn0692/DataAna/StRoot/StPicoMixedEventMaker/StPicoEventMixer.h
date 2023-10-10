#ifndef StPicoEventMixer_hh
#define StPicoEventMixer_hh

/* **************************************************
 * Class stores event buffer used in event mixing. Mixing
 * is done automatically once buffer reaches defined maximum.
 * User should rpesonalize mixEvent() method to cosntruct
 * desired background.
 *
 * **************************************************
 *
 * Initial Authors:
 *          **Michael Lomnitz (mrlomnitz@lbl.gov)
 *          Musta Mustafa   (mmustafa@lbl.gov)
 *
 *  ** Code maintainer
 *
 * **************************************************
 */

#include <vector>
#include "TVector2.h"
#include "StThreeVectorF.hh"
#include "StMixerCuts.h"
#include "TopologyCuts.h"

class TTree;
class TH2F;
class StPicoEvent;
class StPicoTrack;
class StPicoDst;
class StMixerTrack;
class StMixerEvent;
class StMixerPair;
class StD0Hists;
class StEventPlane;

class StPicoEventMixer
{
public:
    StPicoEventMixer(int centBin, int vzBin, int psiBin, StEventPlane* eventPlaneMaker, StD0Hists* d0Hists);
    ~StPicoEventMixer();
    bool addPicoEvent(StPicoDst const* picoDst, StThreeVectorF pVertex, float weight = 1);
    void setEventsBufferSize(int bufferSize);
    void mixEvents();
    void finish();
private:
    int mCentBin, mVzBin, mPsiBin;
    StEventPlane* mEventPlaneMaker;
    
    //TTree * ntp_ME;
    std::vector <StMixerEvent*> mEvents;
    std::vector <StMixerEvent*> mFirstEvents; // first events in a job stored to mix with last events
    StD0Hists* mD0Hists;
    unsigned short int mEventsBufferSize;
    float dca1, dca2, dcaDaughters, theta_hs, decayL_hs;
    float pt_hs, mass_hs, eta_hs, phi_hs;
    
    bool isGoodTrack(StPicoTrack const*, StPicoDst const*, StThreeVectorF const& kfVtx) const;
    bool isCloseTrack(StPicoTrack const*, StThreeVectorF const& pVtx) const;
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*) const;
    bool isPion(StPicoTrack const*, StPicoDst const*, StThreeVectorF const& pVtx) const;
    bool isKaon(StPicoTrack const*, StPicoDst const*, StThreeVectorF const& pVtx) const;
    bool isGoodPair(StMixerPair const& pair, topoCuts::TopologicalCuts const& cuts) const;
    int getD0PtIndex(StMixerPair const& pair, std::vector<float> const& edges) const;
    int getCentIndex(int const icent, std::vector<float> const& edges) const;
    float getTofBeta(StPicoTrack const*, StPicoDst const*, StThreeVectorF const& pVtx) const;
};

inline void StPicoEventMixer::setEventsBufferSize(int bufferSize)
{
    mEventsBufferSize = bufferSize;
}

#endif
