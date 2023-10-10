// $Id$
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoName".Apppend("_Scaled").
//

#include "StSimD0RhoHIOverlay.h"

// ROOT includes
#include <TClonesArray.h>
#include <TMath.h>
#include "TH2.h"
#include "TH2F.h"
#include "TVector3.h"

// jet-framework includes
#include "StJet.h"
#include "StRhoParameter.h"
// #include "StJetMakerTask.h"
#include "StSimD0EventsHIOverlayJetMaker.h"
#include "StCentMaker.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

// STAR centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"


// ROOT classes
class TH2;
class TH2F;

ClassImp(StSimD0RhoHIOverlay)

//________________________________________________________________________
StSimD0RhoHIOverlay::StSimD0RhoHIOverlay() : StRhoBase("")
{
  fNExclLeadJets = 0;
  fJets = 0x0;
  // fJetsArr = 0x0;
  fHistMultvsRho = 0x0;
  // rhovalue = 0x0;
  mOutName = ""; 
  fJetMakerName = "";
  fRhoMakerName = "";
  fD0Kind = 0;
}

//________________________________________________________________________
StSimD0RhoHIOverlay::StSimD0RhoHIOverlay(const char *name, Bool_t histo, const char *outName, const char *jetMakerName) :
  StRhoBase(name, histo, jetMakerName)
{
  fNExclLeadJets = 0;
  fJets = 0x0;
  // fJetsArr = 0x0;
  fHistMultvsRho = 0x0;
  // rhovalue = 0x0;
  mBaseMaker = 0x0;
  mOutName = outName;
  fJetMakerName = jetMakerName;
  fRhoMakerName = name;
  fD0Kind = 0;
  // Constructor.
  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StSimD0RhoHIOverlay::~StSimD0RhoHIOverlay()
{ /*  */
  // destructor
  if(fHistMultvsRho) delete fHistMultvsRho;

}

//________________________________________________________________________
Int_t StSimD0RhoHIOverlay::Init()
{
  // nothing done - base class should take care of that
  // this in effect inherits from StJetFrameworkPicoBase - check it out!
  StRhoBase::Init();

  // declare histogram
  DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);   

  fJetsArr.clear();
  rhovalue.clear();   

  return kStOk;
}

//________________________________________________________________________
Int_t StSimD0RhoHIOverlay::Finish() {
  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(fRhoMakerName);
    fout->cd(fRhoMakerName);
    ///StSimD0RhoHIOverlayBase::WriteHistograms();
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  return kStOK;
}
//
// Function: declare histograms 
//________________________________________________________________________
void StSimD0RhoHIOverlay::DeclareHistograms() {
    // weird that this next line is needed to remove cppcheck warning
    delete fHistMultvsRho;

    // mult was 150, 0, 1500
    fHistMultvsRho = new TH2F("fHistMultvsRho", "fHistMultvsRho", 160, 0., 800., 100, 0., 100.);
    fHistMultvsRho->GetXaxis()->SetTitle("Charged track multiplicity");
    fHistMultvsRho->GetYaxis()->SetTitle("#rho (GeV/c)/A");
}
//
// Function: write histograms / objects to file
//________________________________________________________________________
void StSimD0RhoHIOverlay::WriteHistograms() {
  fHistMultvsRho->Write();
}

//________________________________________________________________________
void StSimD0RhoHIOverlay::Clear(Option_t *opt) {
//  StSimD0RhoHIOverlayBase::Clear();
//  fJets->Clear();
}

//
// Function that runs the analysis for each event
//________________________________________________________________________
Int_t StSimD0RhoHIOverlay::Make() 
{
  fJetsArr.clear();
  rhovalue.clear();   

  // get PicoDstMaker
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();

  // get run number, check bad runs list if desired (kFALSE if bad)
  int fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get vertex 3 vector and declare variables
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // z-vertex cut - per the Aj analysis (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // get JetMaker
  D0JetMaker = static_cast<StSimD0EventsHIOverlayJetMaker*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!D0JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  if(D0JetMaker) {
    fJetsArr =  D0JetMaker->GetJets();
    //fJets->SetName("BGJetsRho");  // name is set by Maker who created it
  }
  if(fJetsArr.size()==0) return kStOk;

  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  double fCent = mCentMaker->GetCentScaled(); // integer bins scaled up by 5% per
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  cout << "RhoMaker not engaged" << endl;
  cout << cent16 << endl;
  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // get event multiplicity - TODO is this correct? same as that used for centrality
  //const int multiplicity = mPicoDst->numberOfTracks(); // this is total tracks not multiplicity
  const double multiplicity = refCorr2;

  // ============================ end of CENTRALITY ============================== //

  cout << "RhoMaker engaged" << endl;

  // initialize Rho and scaled Rho
  fOutRho->SetVal(0);
  if(fOutRhoScaled) fOutRhoScaled->SetVal(0);

  // cout << "Number of Jet Collections = " << fJetsArr.size() << endl;
  for (int jetcollection = 0; jetcollection < fJetsArr.size(); jetcollection++){

    // cout << "Jet Collection # " << jetcollection << endl;
    fJets->Clear();
    // get number of jets, initialize arrays
    fJets = fJetsArr[jetcollection];
    
    const Int_t Njets = fJetsArr[jetcollection]->GetEntries();

    // cout << " ======================================================================================= " << endl;

    cout << "Number of Jets = " <<  Njets << "\t" << fJets->GetEntries() << endl;

    // for (int jet = 0; jet < Njets; jet++){
    //   StJet *jetc = static_cast<StJet*>(fJets->At(jet));
    //   if(!jetc) { continue; } 
    //   cout << Form("Found jet n. %d, pt = %f, area = %f, constituents = %d", jet, jetc->Pt(), jetc->Area(), jetc->GetNumberOfConstituents()) << endl;
    // }

    Int_t maxJetIds[]   = {-1, -1};
    Float_t maxJetPts[] = { 0,  0};

    // exclude leading jets
    if(fNExclLeadJets > 0) {
      // loop over jets
      for(Int_t ij = 0; ij < Njets; ++ij) {
        // get jet pointer
        
        StJet *jet = static_cast<StJet*>(fJets->At(ij));
        if(!jet) { continue; } 

        //if (!AcceptJet(jet)) continue; // FIXME if wanting additional cuts
        // get some jet parameters
        double jetPt = jet->Pt();
        // some threshold cuts for tests
        if(jetPt < 0) continue;

        // get ID and pt of the leading and sub-leading jet   
        if(jetPt > maxJetPts[0]) {
  	maxJetPts[1] = maxJetPts[0];
  	maxJetIds[1] = maxJetIds[0];
  	maxJetPts[0] = jetPt;
  	maxJetIds[0] = ij;
        } else if (jetPt > maxJetPts[1]) {
  	maxJetPts[1] = jetPt;
  	maxJetIds[1] = ij;
        }
      }

      // only set to remove leading jet 
      if(fNExclLeadJets < 2) {
        maxJetIds[1] = -1;
        maxJetPts[1] = 0;
      }
    }

    static Double_t rhovec[999];
    Int_t NjetAcc = 0;

    // push all jets within selected acceptance into stack
    for(Int_t iJets = 0; iJets < Njets; ++iJets) {
      // excluding lead jets
      if(iJets == maxJetIds[0] || iJets == maxJetIds[1])
        continue;

      // get jet pointer
      StJet *jet = static_cast<StJet*>(fJets->At(iJets));
      if(!jet) { continue; } 

      // NEED TO CHECK FOR DEFAULTS - cuts are done at the jet finder level
      //if(!AcceptJet(jet)) continue; //FIXME if wanting additional cuts
      // get some get parameters
      double jetPt = jet->Pt();
      double jetArea = jet->Area();
      // some threshold cuts for tests
      if(jetPt < 0) continue;
   
      rhovec[NjetAcc] = jetPt / jetArea;
      ++NjetAcc;
    }

    // cout << "Number of accepted jets = " << NjetAcc << endl;

    // when we have accepted Jets - calculate and set rho
    if(NjetAcc > 0) {
      // find median value
      Double_t rho = TMath::Median(NjetAcc, rhovec);
      fOutRho->SetVal(rho);

      rhovalue.push_back(fOutRho->GetVal());

      cout << "Rho = " << fOutRho->GetVal() << endl;

      // cout << "Rho value = " << fOutRho->GetVal() << endl;

      // fill histo
      fHistMultvsRho->Fill(multiplicity, rho);

      // if we want scaled Rho from charged -> ch+ne
      if(fOutRhoScaled) {
        //Double_t rhoScaled = rho * GetScaleFactor(fCent); //Don't need yet, fCent is in 5% bins
        Double_t rhoScaled = rho * 1.0;
        fOutRhoScaled->SetVal(rhoScaled);
      }
    }

    // StRhoBase::FillHistograms();

  }
  return kStOk;
} 
