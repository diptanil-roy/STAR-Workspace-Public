#include "StMixerTrack.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StThreeVectorF.hh"
#include <limits>

StMixerTrack::StMixerTrack() : mOrigin(StThreeVectorF()), mMom(StThreeVectorF()), mTrackInfo(std::numeric_limits<short>::min())
{
}
StMixerTrack::StMixerTrack(StThreeVectorF const & pVtx, float B, StPicoTrack const& picoTrack, bool isPion, bool isKaon, TVector2 const & q) :
   mTrackInfo(0)
{
   StPhysicalHelixD helix = picoTrack.helix();
   helix.moveOrigin(helix.pathLength(pVtx));
   mOrigin = helix.origin();
   mMom = helix.momentum(B) * kilogauss;
   mq = q;
   if (picoTrack.charge() == 1) mTrackInfo |= 1;
   if (isPion) mTrackInfo |= 2;
   if (isKaon) mTrackInfo |= 4;
}
StMixerTrack::StMixerTrack(StMixerTrack const * t) : mOrigin(t->mOrigin), mMom(t->mMom), mTrackInfo(t->mTrackInfo)
{
}

