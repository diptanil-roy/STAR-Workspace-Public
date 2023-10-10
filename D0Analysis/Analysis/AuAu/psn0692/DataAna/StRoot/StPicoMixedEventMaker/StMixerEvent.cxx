#include "StMixerEvent.h"
#include "StEventPlane/StEventPlane.h"
#include "TVector2.h"
#include <limits>
StMixerEvent::StMixerEvent() :  mVtx(StThreeVectorF()),
   mBField(std::numeric_limits<float>::quiet_NaN())
{
}
StMixerEvent::StMixerEvent(StMixerEvent *t) : mVtx(t->mVtx), mBField(t->mBField),
					      mWeight(t->mWeight), mQ(t->mQ), mQEtaPlusGap005(t->mQEtaPlusGap005), mQEtaMinusGap005(t->mQEtaMinusGap005),
					      mTracks(t->mTracks),
					      mKaonsIds(t->mKaonsIds), mPionsIds(t->mPionsIds)
{
  for(int i=0; i<20; ++i) mQEta[i] = t->mQEta[i];
}
StMixerEvent::StMixerEvent(StThreeVectorF const& vtx, float b, StEventPlane* eventPlaneMaker, float weight) :
  mVtx(vtx), mBField(b), mWeight(weight)
{
   mQ = eventPlaneMaker->Q();
   mQEtaPlusGap005 = eventPlaneMaker->QEtaPlusGap005();
   mQEtaMinusGap005 = eventPlaneMaker->QEtaMinusGap005();
   for (int i = 0; i < 20; ++i)
      mQEta[i] = eventPlaneMaker->QEta(i);
}
void StMixerEvent::addTrack(StMixerTrack const& t)
{
   mTracks.push_back(t);
   return;
}
void StMixerEvent::addPion(int arrayId)
{
   mPionsIds.push_back(arrayId);
   return;
}
void StMixerEvent::addKaon(int arrayId)
{
   mKaonsIds.push_back(arrayId);
   return;
}
TVector2 StMixerEvent::QEtaGap(int iEta, int nEtaGaps) const
{
   TVector2 QEtaGap_(0, 0);
   int iEta_ = iEta;
   if (iEta_ < nEtaGaps) iEta_ = nEtaGaps - 1;
   if (iEta_ > 20 - nEtaGaps) iEta_ = 20 - nEtaGaps;
   for (int i = 0; i < 20; i++)
   {
      if (fabs(i - iEta_) >= nEtaGaps) QEtaGap_ += mQEta[i];
   }
   return QEtaGap_;
}
