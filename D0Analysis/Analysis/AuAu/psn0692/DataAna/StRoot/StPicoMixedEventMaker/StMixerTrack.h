#ifndef StMixerTrack_hh
#define StMixerTrack_hh
/* **************************************************
 *
 * Track class used for mixed event buffer, stripped down
 * to minimum information neede to reconstruct the helix
 * and basic track information. Currently include:
 * 1) charge
 * 2) isTpcPi & isTofPi
 * 3) isTpcKaon & is TofKaon
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
#include <math.h>
#include "TVector2.h"
#include "StThreeVectorF.hh"

class StPicoTrack;

class StMixerTrack
{
public:
   StMixerTrack();
   StMixerTrack(StMixerTrack const *);
   StMixerTrack(StThreeVectorF const & pVtx, float B, StPicoTrack const& picoTrack, bool isPion, bool isKaon, TVector2 const& q);
   short const getTrackInfo() const;
   int const charge() const ;
   StThreeVectorF const& gMom() const;
   StThreeVectorF const& origin() const;
   TVector2 const q() const;
   ~StMixerTrack()
   {
      ;
   };
private:
   StThreeVectorF mOrigin;
   StThreeVectorF mMom;
   short mTrackInfo;
   TVector2 mq;
   //Removed origin, allt racks shoud me set to 0,0,0
};
inline short const StMixerTrack::getTrackInfo() const
{
   return (mTrackInfo);
}
inline StThreeVectorF const & StMixerTrack::gMom() const
{
   return (mMom) ;
}
inline StThreeVectorF const & StMixerTrack::origin() const
{
   return (mOrigin) ;
}
inline int const StMixerTrack::charge() const
{
   if (mTrackInfo & 1) return 1;
   else return -1;
}
inline TVector2 const StMixerTrack::q() const
{
   return mq;
}
#endif
