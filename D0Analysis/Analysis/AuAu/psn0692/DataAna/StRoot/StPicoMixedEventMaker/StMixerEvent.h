#ifndef StMixerEvent_hh
#define StMixerEvent_hh
#include <math.h>
#include <vector>
#include "TVector2.h"
#include "StThreeVectorF.hh"
#include "StMixerTrack.h"
/* **************************************************
 *
 * Event class used for mixed event buffer, stripped down
 * to minimum information neede to reconstruct the helix
 * and basic track information. Currently include:
 * 1) primVtx
 * 2) B-Field
 * 3) MixerTrack array
 *
 * **************************************************
 *
 *  Initial Authors:
 *         ** Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

class StMixerTrack;
class StEventPlane;

class StMixerEvent
{
public:
   StMixerEvent();
   StMixerEvent(StMixerEvent*);
   StMixerEvent(StThreeVectorF const& vertexPos, float B, StEventPlane* eventPlaneMaker, float weight = 1);
   ~StMixerEvent()
   {
      ;
   };
   void addPion(int);
   void addKaon(int);
   void addTrack(StMixerTrack const&);
   void setPos(float,float,float);
   void setField(float);
   int getNoTracks() const;
   int getNoKaons() const;
   int getNoPions() const;
   int pionId(int counter) const;
   int kaonId(int counter) const;
   StMixerTrack const&  pionAt(int) const;
   StMixerTrack const&  kaonAt(int) const;
   StThreeVectorF const & vertex() const;
   double const field() const;
   float const weight() const;
   TVector2 const Q() const;
   TVector2 const QEtaPlusGap005() const;
   TVector2 const QEtaMinusGap005() const;
   TVector2 QEtaGap(int iEta, int nEtaGaps) const;
private:
   StThreeVectorF mVtx;
   float mBField;
   float mWeight;
   TVector2 mQ;
   TVector2 mQEtaPlusGap005;
   TVector2 mQEtaMinusGap005;
   TVector2 mQEta[20];
   std::vector <StMixerTrack  > mTracks;
   std::vector <int  > mKaonsIds;
   std::vector <int  > mPionsIds;
};
inline void StMixerEvent::setPos(float const vx, float const vy, float const vz)
{
   mVtx = StThreeVectorF(vx, vy, vz);
}
inline void StMixerEvent::setField(float const field)
{
   mBField = field;
}
inline int StMixerEvent::getNoPions() const
{
   return mPionsIds.size();
}
inline int StMixerEvent::getNoKaons() const
{
   return mKaonsIds.size();
}
inline int StMixerEvent::getNoTracks() const
{
   return mTracks.size();
}
inline int StMixerEvent::pionId(int counter) const
{
   return mPionsIds[counter];
}
inline int StMixerEvent::kaonId(int counter) const
{
   return mKaonsIds[counter];
}
inline StMixerTrack const& StMixerEvent::pionAt(int const counter) const
{
   return mTracks[mPionsIds[counter]];
}
inline StMixerTrack const& StMixerEvent::kaonAt(int const counter) const
{
   return mTracks[mKaonsIds[counter]];
}
inline StThreeVectorF const & StMixerEvent::vertex() const
{
   return mVtx;
}
inline double const StMixerEvent::field() const
{
   return mBField;
}
inline float const StMixerEvent::weight() const
{
   return mWeight;
}
inline TVector2 const StMixerEvent::Q() const
{
   return mQ;
}
inline TVector2 const StMixerEvent::QEtaPlusGap005() const
{
   return mQEtaPlusGap005;
}
inline TVector2 const StMixerEvent::QEtaMinusGap005() const
{
   return mQEtaMinusGap005;
}
#endif
