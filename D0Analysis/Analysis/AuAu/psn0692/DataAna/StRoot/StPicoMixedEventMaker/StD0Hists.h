#ifndef StD0Hists__h
#define StD0Hists__h

/* **************************************************
 *  A class to create and save my D0 analysis histograms.
 *
 *  Authors: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include <string>
#include "TObject.h"

class TH1D;
class TH1F;
class TH2F;
class TH3F;
class THn;
class TProfile;
class TFile;
class TString;

#include "StMixerPair.h"
#include "StMixerCuts.h"
#include "TopologyCuts.h"

class StD0Hists
{
public:
  StD0Hists(std::string fileBaseName, int harmonic = 2);
   virtual ~StD0Hists();
   void closeFile();
#ifdef __run_w_QA__
   void fillSameEvt_US_QADist(StMixerPair const&, int const, topoCuts::TopologicalCuts const&);
   void fillSameEvt_LS_QADist(StMixerPair const&, int const, topoCuts::TopologicalCuts const&);
   void fillMixedEvtQADist(StMixerPair const&, int const, topoCuts::TopologicalCuts const&);
#endif
   //event level qa
   TH1D*      hTotalNumberOfEvents;
   TH2F*      hVzVpdVz;
   TH1F*      hVzDiff;
   TH2F*      hVxy;
   TH1F*      hRefMult;
   TH1F*      hGRefMult;
   TH1F*      hTrigger;
   TH1F*      hCentrality;
   TH1F*      hCentralityWeighted;

   TH2F*      hOneOverBetaDiffKaonP;
   TH2F*      hOneOverBetaDiffPionP;

   //mix event counts
   TH3F* hCentVzPsi;
   TH3F* hCentVzPsiSameEventNoWeight;
   TH3F* hCentVzPsiMixedNoWeight;
   TH3F* hCentVzPsiSameEvent;
   TH3F* hCentVzPsiMixed;

#ifdef __run_w_DaugHisto__
   //add daughter pT
   THn* hD0CentPtEtaMDphiDaug[mxeCuts::nCutsSets];
   THn* hD0CentPtEtaMDphiDaugLikeSign[mxeCuts::nCutsSets];
   THn* hD0CentPtEtaMDphiDaugMixed[mxeCuts::nCutsSets];
   THn* hD0CentPtEtaMDphiDaugLikeSignMixed[mxeCuts::nCutsSets];
#endif

   // Half eta's event plane 
   THn* hD0EtaSubCentPtMDphi[mxeCuts::nCutsSets];
   THn* hD0EtaSubCentPtMDphiLikeSign[mxeCuts::nCutsSets];
   THn* hD0EtaSubCentPtMDphiMixed[mxeCuts::nCutsSets];
   THn* hD0EtaSubCentPtMDphiLikeSignMixed[mxeCuts::nCutsSets];

   // Eta Gap Event Plane
   THn* hD0CentPtMDphiEtaGap[mxeCuts::nCutsSets];
   THn* hD0CentPtMDphiEtaGapLikeSign[mxeCuts::nCutsSets];
   THn* hD0CentPtMDphiEtaGapMixed[mxeCuts::nCutsSets];
   THn* hD0CentPtMDphiEtaGapLikeSignMixed[mxeCuts::nCutsSets];


#ifdef __run_w_QA__
   //QA hists
   THn*  mSE_US_DecayTopology;
   TH3F* mSE_US_PointingAngle;
   TH3F* mSE_US_DecayL;
   TH3F* mSE_US_Dca12;
   TH3F* mSE_US_PionDca2Vtx;
   TH3F* mSE_US_KaonDca2Vtx;
   TH3F* mSE_US_D0Dca2Vtx;
   //
   TH3F* mSE_LS_PointingAngle;
   TH3F* mSE_LS_DecayL;
   TH3F* mSE_LS_Dca12;
   TH3F* mSE_LS_PionDca2Vtx;
   TH3F* mSE_LS_KaonDca2Vtx;
   TH3F* mSE_LS_D0Dca2Vtx;
   //
   THn*  mME_US_DecayTopology;
   TH3F* mME_US_PointingAngle;
   TH3F* mME_US_DecayL;
   TH3F* mME_US_Dca12;
   TH3F* mME_US_PionDca2Vtx;
   TH3F* mME_US_KaonDca2Vtx;
   TH3F* mME_US_D0Dca2Vtx;
#endif
};

#endif
