#ifndef StSimD0RhoHIOverlay_H
#define StSimD0RhoHIOverlay_H

// $Id$
// adapted from the AliROOT class AliAnalysisTaskRho.h

#include "StRhoBase.h"

// additional includes
#include "StMaker.h"

// ROOT classes
class TH2;
class TH2F;
class TClonesArray;

// STAR classes
class StMaker;
class StSimD0EventsHIOverlayJetMaker;

class StSimD0RhoHIOverlay : public StRhoBase {

 public:
  StSimD0RhoHIOverlay();
  StSimD0RhoHIOverlay(const char *name, Bool_t histo=kFALSE, const char* outName="", const char* jetMakerName="");
  //virtual ~StSimD0RhoHIOverlay() {}
  virtual ~StSimD0RhoHIOverlay();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  std::vector<double> GetArrayofRhoValue() {return rhovalue;}

  void    SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

  // set Jets for D0s, D0s Bg US, D0s Bg LS
  void         SetD0Kind(Int_t kind)                      { fD0Kind           = kind  ; }

 protected:
  UInt_t            fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  TClonesArray     *fJets;//!jet collection
  std::vector<TClonesArray *> fJetsArr;

  // set Jets for D0s, D0s Bg US, D0s Bg LS
  Double_t               fD0Kind;                 // D0s, Bg D0s US, LS

 private:
  TH2F             *fHistMultvsRho;//!

  StSimD0EventsHIOverlayJetMaker *D0JetMaker;//!

  std::vector<double> rhovalue; //!

  StSimD0RhoHIOverlay(const StSimD0RhoHIOverlay&);             // not implemented
  StSimD0RhoHIOverlay& operator=(const StSimD0RhoHIOverlay&);  // not implemented
  
  ClassDef(StSimD0RhoHIOverlay, 2); // Rho task
};
#endif
