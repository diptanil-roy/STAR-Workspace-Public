#ifndef StMixerCuts_H
#define StMixerCuts_H
/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Michael Lomnitz (mrlomnitz@lbl.gov)
 *
 * **************************************************
 */
// ---------- QA Switch
#define __run_w_QA__
#define __run_w_DaugHisto__
//
#include "Rtypes.h"
#include <string>
#include <vector>

#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "TGraph.h"
#include "StMixerCuts.h"

namespace mxeCuts
{
   enum mePid {kKaon, kPion, kProton};
   float const pidMass[3] = { M_KAON_PLUS, M_PION_PLUS, M_PROTON};
   int const nTrig = 5;
   int const mTriggerId[nTrig] = {450050, 450060,
				  450005, 450015,
				  450025 };
   //Event
   //UShort_t const triggerWord = 0x1F; //first five bits see http://rnc.lbl.gov/~xdong/SoftHadron/picoDst.html
   UShort_t const minBiasTrigger = 0x1F;
   float const maxVz = 6.0;// cm.
   float const vzVpdVz = 3.0; // 3 cm.
   float const Verror = 1.0e-5; //
   float const Vrcut = 2.0; //

   //Tracking
   float const minPt = 0.3;//1.2//0.6
   float const nHitsFit = 20;
   float const dca2pVtx = 0.0050;
   bool const mRequireHft = true;
   float const Eta = 1.0;
   //PID
   // pions
   float const nSigmaPion = 3.0;
   float const tofOneOverBetaDiffPion = 0.03;
   // kaons
   float const nSigmaKaon = 2.0;
   float const tofOneOverBetaDiffKaon = 0.03;

   float const massMin = 0.;
   float const massMax = 2.5;
   int const massBins = 250;
   //Topology QA
   float const QAmassMin = 1.828;
   float const QAmassMax = 1.892;

   //BAsic Decay QA Topology Cuts
   int    const decayTopologyCentrality = 7; // minimumCentrality
   double const decayTopologyMinDca = 0.0050;
   double const decayTopologyMaxD0Dca2Vtx = 0.0120;
   int const nCutsSets = 3;
   string const cutsSetName[nCutsSets] = {"standard","50efficiency","150efficiency"};
}
#endif
