#include <limits>

#include "TH3F.h"
#include "THn.h"
#include "phys_constants.h"
#include "StPicoEventMixer.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StEventPlane/StEventPlane.h"
#include "StPicoMixedEventMaker.h"
#include "StMixerEvent.h"
#include "StMixerPair.h"
#include "StD0Hists.h"
#include "StBTofUtil/tofPathLength.hh"

StPicoEventMixer::StPicoEventMixer(int centBin, int vzBin, int psiBin, StEventPlane* eventPlaneMaker, StD0Hists* d0Hists):
mEvents(), mD0Hists(d0Hists)
{
    mCentBin = centBin;
    mVzBin = vzBin;
    mPsiBin = psiBin;
    mEventPlaneMaker = eventPlaneMaker;
    setEventsBufferSize(11);
    
}
StPicoEventMixer::~StPicoEventMixer()
{
    for (size_t i = 0 ; i < mEvents.size() ; i++)
    {
        delete mEvents[i];
    }
}
void StPicoEventMixer::finish()
{
    if (!mFirstEvents.size())
    {
        // cout << "warning: not enough events to mix!   centBin: " << mCentBin << "  vzBin: " << mVzBin << "  psiBin: " << mPsiBin << "  nEvents: " << mEvents.size() << endl;
        mEventsBufferSize = mEvents.size();
        mixEvents();
    }
    for (int i = 0; i < mEventsBufferSize - 1; i++)
    {
        mEvents.push_back(mFirstEvents.at(0));
        mFirstEvents.erase(mFirstEvents.begin());
        mixEvents();
    }
    
}
bool StPicoEventMixer::addPicoEvent(StPicoDst const* const picoDst, StThreeVectorF pVertex, float weight)
{
    int nTracks = picoDst->numberOfTracks();
    StMixerEvent* event = new StMixerEvent(pVertex, picoDst->event()->bField(), mEventPlaneMaker, weight);
    //Event.setNoTracks( nTracks );
    for (int iTrk = 0; iTrk < nTracks; ++iTrk)
    {
        StPicoTrack const* trk = picoDst->track(iTrk);
        bool saveTrack = false;
        bool isPion_ = false;
        bool isKaon_ = false;
        
        if (!isGoodTrack(trk, picoDst, pVertex )  || isCloseTrack(trk, pVertex)) continue;//good track and Not close trak
        if (isPion(trk, picoDst, pVertex))
        {
            event->addPion(event->getNoTracks());
            isPion_ = true;
            saveTrack = true;
        }
        if (isKaon(trk, picoDst, pVertex))
        {
            event->addKaon(event->getNoTracks());
            isKaon_ = true;
            saveTrack = true;
        }
        if (saveTrack == true)
        {
            StMixerTrack mTrack(pVertex, picoDst->event()->bField(), *trk, isPion_, isKaon_, mEventPlaneMaker->q(iTrk));
            event->addTrack(mTrack);
        }
    }
    //   if (event->getNoPions() > 0 ||  event->getNoKaons() > 0)
    //   {
    mEvents.push_back(event);
    //   }
    /*
     else
     {
     delete event;
     return false;
     }
     */
    //Returns true if need to do mixing, false if buffer has space still
    if (mEvents.size() == mEventsBufferSize)
        return true;
    return false;
}
void StPicoEventMixer::mixEvents()
{
    if(!mEvents.size()) return;
    int mHarmonic = mEventPlaneMaker->harmonic();
    //Template for D0 studies
    for (size_t iEvt2 = 0; iEvt2 < mEvents.size(); ++iEvt2)
    {
        if (iEvt2 == 0)
        {
            mD0Hists->hCentVzPsiSameEventNoWeight->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10);
            mD0Hists->hCentVzPsiSameEvent->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10, mEvents[0]->weight());
        }
        else
        {
            mD0Hists->hCentVzPsiMixedNoWeight->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10);
            mD0Hists->hCentVzPsiMixed->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10, mEvents[0]->weight());
        }
        
        for (int iTrk2 = 0; iTrk2 < mEvents[iEvt2]->getNoKaons(); ++iTrk2) // loop over kaons
        {
            for (int iTrk1 = 0; iTrk1 < mEvents[0]->getNoPions(); ++iTrk1) // loop over pions
            {
                if (iEvt2 == 0 && mEvents[0]->pionId(iTrk1) == mEvents[iEvt2]->kaonId(iTrk2)) continue;
                
                StMixerPair pair(mEvents[0]->pionAt(iTrk1), mEvents[iEvt2]->kaonAt(iTrk2),
                                 mxeCuts::pidMass[mxeCuts::kPion], mxeCuts::pidMass[mxeCuts::kKaon],
                                 mEvents[0]->vertex(), mEvents[iEvt2]->vertex(),
                                 mEvents[0]->field());
                
                int charge2 = mEvents[0]->pionAt(iTrk1).charge() * mEvents[iEvt2]->kaonAt(iTrk2).charge();
                
                //Topology histos, fill before checking cuts
#ifdef __run_w_QA__
                if (iEvt2 == 0)
                {
                    if (charge2 < 0) mD0Hists->fillSameEvt_US_QADist(pair, mCentBin, topoCuts::D0Cuts);
                    else mD0Hists->fillSameEvt_LS_QADist(pair, mCentBin, topoCuts::D0Cuts);
                }
                else if (charge2 < 0) mD0Hists->fillMixedEvtQADist(pair, mCentBin, topoCuts::D0Cuts);
#endif
                bool pass_cut_set[3] = {isGoodPair(&pair,topoCuts::D0Cuts),
                    isGoodPair(&pair,topoCuts::D0Cuts_50eff),
                    isGoodPair(&pair,topoCuts::D0Cuts_150eff) };
                if( !pass_cut_set[0] && !pass_cut_set[1] && !pass_cut_set[2]) continue;
                
                TVector2 QSub = mEvents[0]->Q() - mEvents[0]->pionAt(iTrk1).q();
                if (iEvt2 == 0) QSub -= mEvents[iEvt2]->kaonAt(iTrk2).q();
                float dPhi = pair.phi() - QSub.Phi() / mHarmonic;
                while (dPhi < 0) dPhi += (2.0/mHarmonic)*TMath::Pi();
                while (dPhi >= (2.0/mHarmonic)*TMath::Pi()) dPhi -= (2.0/mHarmonic)*TMath::Pi();
                
                double toFill[5] = {mCentBin + 0.5, pair.pt(), pair.eta(), pair.m(), dPhi};
                // double toFillDaug[5] = {mCentBin + 0.5, pair.pt(), mEvents[0]->pionAt(iTrk1).gMom().perp(), pair.m(), mEvents[iEvt2]->kaonAt(iTrk2).gMom().perp()};
                double toFillDaug[7] = {mCentBin + 0.5, pair.pt(), mEvents[0]->pionAt(iTrk1).gMom().perp(), pair.m(), mEvents[iEvt2]->kaonAt(iTrk2).gMom().perp(), pair.lorentzVector().rapidity(), mEvents[0]->pionAt(iTrk1).charge() };
                //Daughter histograms
#ifdef __run_w_DaugHisto__
                for(int ii = 0; ii < mxeCuts::nCutsSets; ++ii){
                    if( !pass_cut_set[ii] ) continue;
                    if (iEvt2 == 0)
                    {
                        if (charge2 < 0)
                        {
                            mD0Hists->hD0CentPtEtaMDphiDaug[ii]->Fill(toFillDaug, mEvents[0]->weight());
                        }
                        else
                        {
                            mD0Hists->hD0CentPtEtaMDphiDaugLikeSign[ii]->Fill(toFillDaug, mEvents[0]->weight());
                        }
                    }
                    else
                    {
                        if (charge2 < 0)
                        {
                            mD0Hists->hD0CentPtEtaMDphiDaugMixed[ii]->Fill(toFillDaug, mEvents[0]->weight());
                        }
                        else
                        {
                            mD0Hists->hD0CentPtEtaMDphiDaugLikeSignMixed[ii]->Fill(toFillDaug, mEvents[0]->weight());
                        }
                    }
                }
#endif
                // Eta sub
                TVector2 QSubEtaSub;
                if( pair.eta() > 0 )
                    QSubEtaSub = mEvents[0]->QEtaMinusGap005();
                else
                    QSubEtaSub = mEvents[0]->QEtaPlusGap005();
                
                float pionEta = mEvents[0]->pionAt(iTrk1).gMom().pseudoRapidity();
                float kaonEta = mEvents[iEvt2]->kaonAt(iTrk2).gMom().pseudoRapidity();
                if( pionEta*pair.eta() < 0 && std::fabs(pionEta) > 0.05  )
                    QSubEtaSub -= mEvents[0]->pionAt(iTrk1).q();
                if( iEvt2==0 && kaonEta*pair.eta() < 0 && std::fabs(kaonEta) > 0.05  )
                    QSubEtaSub -= mEvents[0]->kaonAt(iTrk2).q();
                
                float dPhiEtaSub = pair.phi() - QSubEtaSub.Phi() / mHarmonic;
                while (dPhiEtaSub < 0) dPhiEtaSub += (2.0/mHarmonic)*TMath::Pi();
                while (dPhiEtaSub >= (2.0/mHarmonic)*TMath::Pi()) dPhiEtaSub -= (2.0/mHarmonic)*TMath::Pi();
                for(int ii = 0; ii < mxeCuts::nCutsSets; ++ii){
                    if( !pass_cut_set[ii] ) continue;
                    double toFillEtaSub[5] = {mCentBin + 0.5, pair.pt(), pair.m(), dPhiEtaSub};
                    if (iEvt2 == 0)
                    {
                        if (charge2 < 0)
                        {
                            // mD0Hists->hD0EtaSubCentPtMDphi[ii]->Fill(toFillEtaSub, mEvents[0]->weight());
                        }
                        else
                        {
                            // mD0Hists->hD0EtaSubCentPtMDphiLikeSign[ii]->Fill(toFillEtaSub, mEvents[0]->weight());
                        }
                    }
                    else
                    {
                        if (charge2 < 0)
                        {
                            // mD0Hists->hD0EtaSubCentPtMDphiMixed[ii]->Fill(toFillEtaSub, mEvents[0]->weight());
                        }
                        else
                        {
                            // mD0Hists->hD0EtaSubCentPtMDphiLikeSignMixed[ii]->Fill(toFillEtaSub, mEvents[0]->weight());
                        }
                    }
                }
                
                // Eta gap
                int iEta = (int)(pair.eta() * 10 + 10);
                for (int nEtaGaps = 0; nEtaGaps < 8; ++nEtaGaps)
                {
                    TVector2 QSubEtaGap = mEvents[0]->QEtaGap(iEta, nEtaGaps);
                    int iEta_ = iEta;
                    if (iEta_ < nEtaGaps) iEta_ = nEtaGaps - 1;
                    if (iEta_ > 20 - nEtaGaps) iEta_ = 20 - nEtaGaps;
                    int iEtaPion = (int)(mEvents[0]->pionAt(iTrk1).gMom().pseudoRapidity() * 10 + 10);
                    if (fabs(iEtaPion - iEta_) >= nEtaGaps)
                        QSubEtaGap -= mEvents[0]->pionAt(iTrk1).q();
                    if (iEvt2 == 0)
                    {
                        int iEtaKaon = (int)(mEvents[iEvt2]->kaonAt(iTrk2).gMom().pseudoRapidity() * 10 + 10);
                        if (fabs(iEtaKaon - iEta_) >= nEtaGaps)
                            QSubEtaGap -= mEvents[iEvt2]->kaonAt(iTrk2).q();
                    }
                    if (QSubEtaGap.Mod() == 0)
                    {
                        //cout << "QSubEtaGap.Mod()==0  nEtaGaps: " << nEtaGaps << endl;
                        continue;
                    }
                    float dPhiEtaGap = pair.phi() - QSubEtaGap.Phi() / mHarmonic;
                    while (dPhiEtaGap < 0) dPhiEtaGap += (2.0/mHarmonic)*TMath::Pi();
                    while (dPhiEtaGap >= (2.0/mHarmonic)*TMath::Pi()) dPhiEtaGap -= (2.0/mHarmonic)*TMath::Pi();
                    double toFill[5] = {mCentBin + 0.5, pair.pt(), pair.m(), dPhiEtaGap, 0.1 * nEtaGaps + 0.05};
                    for(int ii = 0; ii < mxeCuts::nCutsSets; ++ii){
                        if( !pass_cut_set[ii] ) continue;
                        if (iEvt2 == 0)
                        {
                            // if (charge2 < 0) mD0Hists->hD0CentPtMDphiEtaGap[ii]->Fill(toFill, mEvents[0]->weight());
                            // else mD0Hists->hD0CentPtMDphiEtaGapLikeSign[ii]->Fill(toFill, mEvents[0]->weight());
                        }
                        else
                        {
                            // if (charge2 < 0) mD0Hists->hD0CentPtMDphiEtaGapMixed[ii]->Fill(toFill, mEvents[0]->weight());
                            // else mD0Hists->hD0CentPtMDphiEtaGapLikeSignMixed[ii]->Fill(toFill, mEvents[0]->weight());
                        }
                    }
                }
                
            } //second event track loop
        } //first event track loop
    } //loop over second events
    
    if (mFirstEvents.size() == static_cast<unsigned short>(mEventsBufferSize - 1))
        delete mEvents[0];
    else
        mFirstEvents.push_back(mEvents[0]);
    
    mEvents.erase(mEvents.begin());
}
bool StPicoEventMixer::isKaon(StPicoTrack const* const trk, StPicoDst const* const picoDst, StThreeVectorF const& pVertex) const
{
    if (!isTpcKaon(trk)) return false;
    float beta = getTofBeta(trk, picoDst, pVertex);
    if (isnan(beta) || beta < 0) return true;
    float p = trk->gMom(pVertex, picoDst->event()->bField()).mag();
    float oneOverBetaExpected = sqrt(M_KAON_PLUS*M_KAON_PLUS / p / p + 1);
    
    return fabs(1. / beta - oneOverBetaExpected) < mxeCuts::tofOneOverBetaDiffPion;
}
bool StPicoEventMixer::isPion(StPicoTrack const* const trk, StPicoDst const* const picoDst, StThreeVectorF const& pVertex) const
{
    if (!isTpcPion(trk)) return false;
    float beta = getTofBeta(trk, picoDst, pVertex);
    if (isnan(beta) || beta < 0) return true;
    float p = trk->gMom(pVertex, picoDst->event()->bField()).mag();
    float oneOverBetaExpected = sqrt(M_PION_PLUS*M_PION_PLUS / p / p + 1);
    
    return fabs(1. / beta - oneOverBetaExpected) < mxeCuts::tofOneOverBetaDiffKaon;
}
bool StPicoEventMixer::isTpcPion(StPicoTrack const * const trk) const
{
    return fabs(trk->nSigmaPion()) < mxeCuts::nSigmaPion;
}
bool StPicoEventMixer::isTpcKaon(StPicoTrack const * const trk) const
{
    return fabs(trk->nSigmaKaon()) < mxeCuts::nSigmaKaon;
}
bool StPicoEventMixer::isGoodTrack(StPicoTrack const * const trk, StPicoDst const* const picoDst,  StThreeVectorF const& kfVtx) const
{
    StThreeVectorF mom = trk->gMom(kfVtx, picoDst->event()->bField());
    return ((!mxeCuts::mRequireHft || trk->isHFTTrack()) &&
            trk->nHitsFit() >= mxeCuts::nHitsFit && trk->gPt() > mxeCuts::minPt && fabs(mom.pseudoRapidity()) <= mxeCuts::Eta);
}
bool StPicoEventMixer::isCloseTrack(StPicoTrack const* const trk, StThreeVectorF const& pVtx) const
{
    StPhysicalHelixD helix = trk->dcaGeometry().helix();
    return (helix.at(helix.pathLength(pVtx)) - pVtx).mag() <= mxeCuts::dca2pVtx;
}
bool StPicoEventMixer::isGoodPair(StMixerPair const& pair, topoCuts::TopologicalCuts const& cuts) const
{
    int ptIndex = getD0PtIndex(pair,cuts.PtEdge);
    int centIndex = getCentIndex(mCentBin,cuts.CentEdge);
    if(ptIndex==-1 || centIndex==-1) return false;
    return (pair.m() > mxeCuts::massMin && pair.m() < mxeCuts::massMax &&
            std::abs(pair.lorentzVector().rapidity()) < cuts.RapidityCut &&
            pair.particle1Dca() > cuts.pDca[centIndex][ptIndex] && pair.particle2Dca() > cuts.kDca[centIndex][ptIndex] &&
            pair.dcaDaughters() < cuts.dcaDaughters[centIndex][ptIndex] &&
            pair.decayLength() > cuts.decayLength[centIndex][ptIndex] &&
            std::cos(pair.pointingAngle()) > cuts.cosTheta[centIndex][ptIndex] &&
            ((pair.decayLength()) * sin(pair.pointingAngle())) < cuts.dcaV0ToPv[centIndex][ptIndex]);
}//----------------------------------------------------------------------------- 
//-----------------------------------------------------------------------------
int StPicoEventMixer::getD0PtIndex(StMixerPair const& pair, std::vector<float> const& edges) const
{
    for (int i = 0; i < edges.size(); i++)
    {
        if ((pair.pt() >= edges[i]) && (pair.pt() < edges[i + 1]))
            return i;
    }
    return - 1;
}
//-----------------------------------------------------------------------------
int StPicoEventMixer::getCentIndex(int const icent, std::vector<float> const& edges) const
{
    for (int i = 0; i < edges.size(); i++)
    {
        if ((icent >= edges[i]) && (icent < edges[i + 1]))
            return i;
    }
    return - 1;
}
//-----------------------------------------------------------------------------
float StPicoEventMixer::getTofBeta(StPicoTrack const* const trk, StPicoDst const* const picoDst, StThreeVectorF const& pVertex) const
{
    int index2tof = trk->bTofPidTraitsIndex();
    
    float beta = std::numeric_limits<float>::quiet_NaN();
    
    if (index2tof >= 0)
    {
        StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2tof);
        
        if (tofPid)
        {
            beta = tofPid->btofBeta();
            
            if (beta < 1e-4)
            {
                StThreeVectorF const btofHitPos = tofPid->btofHitPos();
                
                StPhysicalHelixD helix = trk->helix();
                float L = tofPathLength(&pVertex, &btofHitPos, helix.curvature());
                float tof = tofPid->btof();
                if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
                else beta = std::numeric_limits<float>::quiet_NaN();
            }
        }
    }
    
    return beta;
}
