#ifndef StFJWrapper_H
#define StFJWrapper_H

#if !defined(__CINT__)

// ROOT includes
#include <TMath.h>
#include <TList.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <vector>
#include <TString.h>

// jet includes
#include "FJ_includes.h"
#include "StJetPicoDefinitions.h"

class StFJWrapper
{
 public:
  StFJWrapper(const char *name, const char *title);
  virtual ~StFJWrapper();

  virtual void  AddInputVector (Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index = -99999);
  virtual void  AddInputVector (const fastjet::PseudoJet& vec,                Int_t index = -99999);
  virtual void  AddInputVectors(const std::vector<fastjet::PseudoJet>& vecs,  Int_t offsetIndex = -99999);
  virtual void  AddInputGhost  (Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index = -99999);
  virtual const char *ClassName()                            const { return "StFJWrapper";              }
  virtual void  Clear(const Option_t* /*opt*/ = "");
  virtual void  ClearMemory();
  virtual void  CopySettingsFrom (const StFJWrapper& wrapper);
  virtual void  GetMedianAndSigma(Double_t& median, Double_t& sigma, Int_t remove = 0) const;
  // void          Description();
  fastjet::ClusterSequenceArea*           GetClusterSequence() const   { return fClustSeq;                 }
  fastjet::ClusterSequence*               GetClusterSequenceSA() const { return fClustSeqSA;               }
  fastjet::ClusterSequenceActiveAreaExplicitGhosts* GetClusterSequenceGhosts() const { return fClustSeqActGhosts; }
  const std::vector<fastjet::PseudoJet>&  GetInputVectors()    const { return fInputVectors;               }
  const std::vector<fastjet::PseudoJet>&  GetInputGhosts()     const { return fInputGhosts;                }
  const std::vector<fastjet::PseudoJet>&  GetInclusiveJets()   const { return fInclusiveJets;              }
  const std::vector<fastjet::PseudoJet>&  GetFilteredJets()    const { return fFilteredJets;               }
  std::vector<fastjet::PseudoJet>         GetJetConstituents(UInt_t idx) const;
  std::vector<fastjet::PseudoJet>         GetFilteredJetConstituents(UInt_t idx) const;
  Double_t                                GetMedianUsedForBgSubtraction() const { return fMedUsedForBgSub; }
  const char*                             GetName()            const { return fName;                       }
  const char*                             GetTitle()           const { return fTitle;                      }
  Double_t                                GetJetArea         (UInt_t idx) const;
  fastjet::PseudoJet                      GetJetAreaVector   (UInt_t idx) const;
  Double_t                                GetFilteredJetArea (UInt_t idx) const;
  fastjet::PseudoJet                      GetFilteredJetAreaVector(UInt_t idx) const;
  Double_t                                GetJetSubtractedPt (UInt_t idx) const;
  virtual std::vector<double>             GetSubtractedJetsPts(Double_t median_pt = -1, Bool_t sorted = kFALSE);
  Bool_t                                  GetLegacyMode()            { return fLegacyMode; }
  Bool_t                                  GetDoFilterArea()          { return fDoFilterArea; }
  Double_t                                NSubjettiness(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Int_t Option=0);
  Double32_t                              NSubjettinessDerivativeSub(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Double_t JetR, fastjet::PseudoJet jet, Int_t Option=0);
#ifdef FASTJET_VERSION
/*
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetMass()        const {return fGenSubtractorInfoJetMass        ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetAngularity()  const {return fGenSubtractorInfoJetAngularity  ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetpTD()         const {return fGenSubtractorInfoJetpTD         ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetCircularity() const {return fGenSubtractorInfoJetCircularity ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetSigma2()      const {return fGenSubtractorInfoJetSigma2      ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetConstituent() const {return fGenSubtractorInfoJetConstituent ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetLeSub()       const {return fGenSubtractorInfoJetLeSub       ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJet1subjettiness_kt()       const {return fGenSubtractorInfoJet1subjettiness_kt ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJet2subjettiness_kt()       const {return fGenSubtractorInfoJet2subjettiness_kt ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJet3subjettiness_kt()       const {return fGenSubtractorInfoJet3subjettiness_kt ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetOpeningAngle_kt()       const {return fGenSubtractorInfoJetOpeningAngle_kt ; }
*/
  const std::vector<fastjet::PseudoJet>                      GetConstituentSubtrJets()            const { return fConstituentSubtrJets           ; }
  const std::vector<fastjet::PseudoJet>                      GetGroomedJets()                     const { return fGroomedJets                    ; }
////  Int_t CreateGenSub();          // fastjet::contrib::GenericSubtractor
  Int_t CreateConstituentSub();  // fastjet::contrib::ConstituentSubtractor
  Int_t CreateSoftDrop();
#endif
  virtual std::vector<double>                                GetGRNumerator()                     const { return fGRNumerator                    ; }
  virtual std::vector<double>                                GetGRDenominator()                   const { return fGRDenominator                  ; }
  virtual std::vector<double>                                GetGRNumeratorSub()                  const { return fGRNumeratorSub                 ; }
  virtual std::vector<double>                                GetGRDenominatorSub()                const { return fGRDenominatorSub               ; }

  virtual void RemoveLastInputVector();

  virtual Int_t Run();
  virtual Int_t Filter();
/*
  virtual Int_t DoGenericSubtractionJetMass();
  virtual Int_t DoGenericSubtractionGR(Int_t ijet);
  virtual Int_t DoGenericSubtractionJetAngularity();
  virtual Int_t DoGenericSubtractionJetpTD();
  virtual Int_t DoGenericSubtractionJetCircularity();
  virtual Int_t DoGenericSubtractionJetSigma2();
  virtual Int_t DoGenericSubtractionJetConstituent();
  virtual Int_t DoGenericSubtractionJetLeSub();
  virtual Int_t DoGenericSubtractionJet1subjettiness_kt();
  virtual Int_t DoGenericSubtractionJet2subjettiness_kt();
  virtual Int_t DoGenericSubtractionJet3subjettiness_kt();
  virtual Int_t DoGenericSubtractionJetOpeningAngle_kt();
*/
  virtual Int_t DoConstituentSubtraction();
  virtual Int_t DoSoftDrop();
  
  void SetName(const char* name)        { fName           = name;    }
  void SetTitle(const char* title)      { fTitle          = title;   }
  void SetStrategy(const fastjet::Strategy &strat)                 { fStrategy = strat;  }
  void SetAlgorithm(const fastjet::JetAlgorithm &algor)            { fAlgor    = algor;  }
  void SetRecombScheme(const fastjet::RecombinationScheme &scheme) { fScheme   = scheme; }
  void SetAreaType(const fastjet::AreaType &atype)                 { fAreaType = atype;  }
  void SetNRepeats(Int_t nrepeat)       { fNGhostRepeats  = nrepeat; }
  void SetGhostArea(Double_t gharea)    { fGhostArea      = gharea;  }
  void SetMaxRap(Double_t maxrap)       { fMaxRap         = maxrap;  }
  void SetR(Double_t r)                 { fR              = r;       }
  void SetGridScatter(Double_t gridSc)  { fGridScatter    = gridSc;  }
  void SetKtScatter(Double_t ktSc)      { fKtScatter      = ktSc;    }
  void SetMeanGhostKt(Double_t meankt)  { fMeanGhostKt    = meankt;  }
  void SetPluginAlgor(Int_t plugin)     { fPluginAlgor    = plugin;  }
  void SetUseArea4Vector(Bool_t useA4v) { fUseArea4Vector = useA4v;  }
  void SetupAlgorithmfromOpt(const char *option);
  void SetupAreaTypefromOpt(const char *option);
  void SetupSchemefromOpt(const char *option);
  void SetupStrategyfromOpt(const char *option);
  void SetLegacyMode (Bool_t mode)      { fLegacyMode ^= mode; }
  void SetLegacyFJ();
  void SetUseExternalBkg(Bool_t b, Double_t rho, Double_t rhom) { fUseExternalBkg = b; fRho = rho; fRhom = rhom;}
  void SetRMaxAndStep(Double_t rmax, Double_t dr) {fRMax = rmax; fDRStep = dr; }
  void SetRhoRhom (Double_t rho, Double_t rhom) { fUseExternalBkg = kTRUE; fRho = rho; fRhom = rhom;} // if using rho,rhom then fUseExternalBkg is true
  void SetMinJetPt(Double_t MinPt) {fMinJetPt=MinPt;}

 protected:
  TString                                fName;               //!
  TString                                fTitle;              //!
  std::vector<fastjet::PseudoJet>        fInputVectors;       //!
  std::vector<fastjet::PseudoJet>        fInputGhosts;        //!
  std::vector<fastjet::PseudoJet>        fInclusiveJets;      //!
  std::vector<fastjet::PseudoJet>        fFilteredJets;       //!
  std::vector<double>                    fSubtractedJetsPt;   //!
  std::vector<fastjet::PseudoJet>        fConstituentSubtrJets; //!
  std::vector<fastjet::PseudoJet>        fGroomedJets;        //!
  fastjet::AreaDefinition               *fAreaDef;            //!
  fastjet::VoronoiAreaSpec              *fVorAreaSpec;        //!
  fastjet::GhostedAreaSpec              *fGhostedAreaSpec;    //!
  fastjet::JetDefinition                *fJetDef;             //!
  fastjet::JetDefinition::Plugin        *fPlugin;             //!
#ifndef FASTJET_VERSION
  fastjet::RangeDefinition              *fRange;              //!
#else
  fastjet::Selector                     *fRange;              //!
#endif
  fastjet::ClusterSequenceArea          *fClustSeq;           //!
  fastjet::ClusterSequence              *fClustSeqSA;                //!
  fastjet::ClusterSequenceActiveAreaExplicitGhosts *fClustSeqActGhosts; //!
  fastjet::Strategy                      fStrategy;           //!
  fastjet::JetAlgorithm                  fAlgor;              //!
  fastjet::RecombinationScheme           fScheme;             //!
  fastjet::AreaType                      fAreaType;           //!
  Int_t                                  fNGhostRepeats;      //!
  Double_t                               fGhostArea;	      //!
  Double_t                               fMaxRap;	      //!
  Double_t                               fR;                  //!
  Double_t                               fMinJetPt;
  // no setters for the moment - used default values in the constructor
  Double_t                               fGridScatter;        //!
  Double_t                               fKtScatter;	      //!
  Double_t                               fMeanGhostKt;        //!
  Int_t                                  fPluginAlgor;        //!
  // extra parameters
  Double_t                               fMedUsedForBgSub;    //!
  Bool_t                                 fUseArea4Vector;     //!
  // condition to stop the grooming (rejection of soft splitting) z > fZcut theta^fBeta
  Double_t                               fZcut;               // fZcut = 0.1                
  Double_t                               fBeta;               // fBeta = 0
#ifdef FASTJET_VERSION
  fastjet::JetMedianBackgroundEstimator   *fBkrdEstimator;    //!
  //from contrib package
////  fastjet::contrib::GenericSubtractor     *fGenSubtractor;    //!
  fastjet::contrib::ConstituentSubtractor *fConstituentSubtractor;    //!
  fastjet::contrib::SoftDrop              *fSoftDrop;        //!
/*
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetMass;        //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoGRNum;          //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoGRDen;          //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetAngularity;  //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetpTD;         //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetCircularity; //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetSigma2;      //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetConstituent; //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetLeSub;       //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJet1subjettiness_kt;       //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJet2subjettiness_kt;       //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJet3subjettiness_kt;       //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetOpeningAngle_kt;       //!
*/
#endif
  Bool_t                                   fDoFilterArea;         //!
  Bool_t                                   fLegacyMode;           //!
  Bool_t                                   fUseExternalBkg;       //!
  Double_t                                 fRho;                  //  pT background density
  Double_t                                 fRhom;                 //  mT background density
  Double_t                                 fRMax;             //!
  Double_t                                 fDRStep;           //!
  std::vector<double>                      fGRNumerator;      //!
  std::vector<double>                      fGRDenominator;    //!
  std::vector<double>                      fGRNumeratorSub;   //!
  std::vector<double>                      fGRDenominatorSub; //!

  virtual void   SubtractBackground(const Double_t median_pt = -1);

 private:
  StFJWrapper();
  StFJWrapper(const StFJWrapper& wrapper);
  StFJWrapper& operator = (const StFJWrapper& wrapper);
};
#endif
#endif

#ifdef StFJWrapper_CXX
#undef StFJWrapper_CXX

#if defined __GNUC__
#pragma GCC system_header
#endif

namespace fj = fastjet;

//_________________________________________________________________________________________________
StFJWrapper::StFJWrapper(const char *name, const char *title)
  :
    fName              (name)
  , fTitle             (title)
  , fInputVectors      ( )
  , fInputGhosts       ( )
  , fInclusiveJets     ( )
  , fFilteredJets      ( )
  , fSubtractedJetsPt  ( )
  , fConstituentSubtrJets ( )
  , fSoftDrop          ( )
  , fAreaDef           (0)
  , fVorAreaSpec       (0)
  , fGhostedAreaSpec   (0)
  , fJetDef            (0)
  , fPlugin            (0)
  , fRange             (0)
  , fClustSeq          (0)
  , fClustSeqSA        (0)
  , fClustSeqActGhosts (0)
  , fStrategy          (fj::Best)
  , fAlgor             (fj::kt_algorithm)
  , fScheme            (fj::BIpt_scheme)
  , fAreaType          (fj::active_area)
  , fNGhostRepeats     (1)
  , fGhostArea         (0.005)
  , fMaxRap            (1.)
  , fR                 (0.4)
  , fGridScatter       (1.0)
  , fKtScatter         (0.1)
  , fMeanGhostKt       (1e-100)
  , fPluginAlgor       (0)
  , fMedUsedForBgSub   (0)
  , fUseArea4Vector    (kFALSE)
  , fZcut(0.1)
  , fBeta(0)
#ifdef FASTJET_VERSION
  , fBkrdEstimator     (0)
////  , fGenSubtractor     (0)
  , fConstituentSubtractor (0)
/*
  , fGenSubtractorInfoJetMass ( )
  , fGenSubtractorInfoGRNum ( )
  , fGenSubtractorInfoGRDen ( )
  , fGenSubtractorInfoJetAngularity ( )
  , fGenSubtractorInfoJetpTD ( )
  , fGenSubtractorInfoJetCircularity( )
  , fGenSubtractorInfoJetSigma2()
  , fGenSubtractorInfoJetConstituent ( )
  , fGenSubtractorInfoJetLeSub ( )
  , fGenSubtractorInfoJet1subjettiness_kt ( )
  , fGenSubtractorInfoJet2subjettiness_kt ( )
  , fGenSubtractorInfoJet3subjettiness_kt ( )
  , fGenSubtractorInfoJetOpeningAngle_kt ( )
*/
#endif
  , fDoFilterArea      (false)
  , fLegacyMode        (false)
  , fUseExternalBkg    (false)
  , fRho               (0)
  , fRhom              (0)
  , fRMax(2.)
  , fDRStep(0.04)
  , fGRNumerator()
  , fGRDenominator()
  , fGRNumeratorSub()
  , fGRDenominatorSub()
{
  // Constructor.
}

//_________________________________________________________________________________________________
StFJWrapper::~StFJWrapper()
{
  // Destructor.
  ClearMemory();
}

//_________________________________________________________________________________________________
void StFJWrapper::ClearMemory()
{
  // Destructor.
  if (fAreaDef)           { delete fAreaDef;           fAreaDef         = NULL; }
  if (fVorAreaSpec)       { delete fVorAreaSpec;       fVorAreaSpec     = NULL; }
  if (fGhostedAreaSpec)   { delete fGhostedAreaSpec;   fGhostedAreaSpec = NULL; }
  if (fJetDef)            { delete fJetDef;            fJetDef          = NULL; }
  if (fPlugin)            { delete fPlugin;            fPlugin          = NULL; }
  if (fRange)             { delete fRange;             fRange           = NULL; }
  if (fClustSeq)          { delete fClustSeq;          fClustSeq        = NULL; }
  if (fClustSeqSA)        { delete fClustSeqSA;        fClustSeqSA        = NULL; }
  if (fClustSeqActGhosts) { delete fClustSeqActGhosts; fClustSeqActGhosts = NULL; }
  #ifdef FASTJET_VERSION
  if (fBkrdEstimator)          { delete fBkrdEstimator; fBkrdEstimator = NULL; }
////  if (fGenSubtractor)          { delete fGenSubtractor; fGenSubtractor = NULL; }
  if (fConstituentSubtractor)  { delete fConstituentSubtractor; fConstituentSubtractor = NULL; }
  if (fSoftDrop)          { delete fSoftDrop; fSoftDrop = NULL;}
  #endif
}

//_________________________________________________________________________________________________
void StFJWrapper::CopySettingsFrom(const StFJWrapper& wrapper)
{
  // Copy some settings.
  // You very often want to keep most of the settings
  // but change only the algorithm or R - do it after call to this function

  fStrategy         = wrapper.fStrategy;
  fAlgor            = wrapper.fAlgor;
  fScheme           = wrapper.fScheme;
  fAreaType         = wrapper.fAreaType;
  fNGhostRepeats    = wrapper.fNGhostRepeats;
  fGhostArea        = wrapper.fGhostArea;
  fMaxRap           = wrapper.fMaxRap;
  fR                = wrapper.fR;
  fGridScatter      = wrapper.fGridScatter;
  fKtScatter        = wrapper.fKtScatter;
  fMeanGhostKt      = wrapper.fMeanGhostKt;
  fPluginAlgor      = wrapper.fPluginAlgor;
  fUseArea4Vector   = wrapper.fUseArea4Vector;
  fZcut             = wrapper.fZcut;
  fBeta             = wrapper.fBeta;
  fLegacyMode       = wrapper.fLegacyMode;
  fUseExternalBkg   = wrapper.fUseExternalBkg;
  fRho              = wrapper.fRho;
  fRhom             = wrapper.fRhom;
}

//_________________________________________________________________________________________________
void StFJWrapper::Clear(const Option_t */*opt*/)
{
  // Simply clear the input vectors.
  // Make sure done on every event if the instance is reused
  // Reset the median to zero.

  fInputVectors.clear();
  fInputGhosts.clear();
  fMedUsedForBgSub = 0;

  // for the moment brute force delete everything
  ClearMemory();
}

//_________________________________________________________________________________________________
void StFJWrapper::RemoveLastInputVector()
{
  // Remove last input vector
  fInputVectors.pop_back();
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputVector(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index)
{
  // Make the input pseudojet.
  fastjet::PseudoJet inVec(px, py, pz, E);

  // Salvatore Aiola: not sure why this was done...
  //if (index > -99999) {
  inVec.set_user_index(index);
  //} else {
  //inVec.set_user_index(fInputVectors.size());
  //}

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputVector(const fj::PseudoJet& vec, Int_t index)
{
  // Add an input pseudojet.
  fj::PseudoJet inVec = vec;

  // Salvatore Aiola: not sure why this was done...
  ///if (index > -99999) {
  inVec.set_user_index(index);
  //} else {
  //inVec.set_user_index(fInputVectors.size());
  //}

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputVectors(const std::vector<fj::PseudoJet>& vecs, Int_t offsetIndex)
{
  // Add the input from vector of pseudojets.
  for (UInt_t i = 0; i < vecs.size(); ++i) {
    fj::PseudoJet inVec = vecs[i];
    if (offsetIndex > -99999)
      inVec.set_user_index(fInputVectors.size() + offsetIndex);
    // add to the fj container of input vectors
    fInputVectors.push_back(inVec);
  }
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputGhost(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index)
{
  // Make the input pseudojet.
  fastjet::PseudoJet inVec(px, py, pz, E);

  if (index > -99999) {
    inVec.set_user_index(index);
  } else {
    inVec.set_user_index(fInputGhosts.size());
  }

  // add to the fj container of input vectors
  fInputGhosts.push_back(inVec);
  if (!fDoFilterArea) fDoFilterArea = kTRUE;
}

//_________________________________________________________________________________________________
Double_t StFJWrapper::GetJetArea(UInt_t idx) const
{
  // Get the jet area.
  Double_t retval = -1; // really wrong area..
  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->area(fInclusiveJets[idx]);
  } else {
    __ERROR(Form("Wrong index: %d",idx));
  }
  return retval;
}

//_________________________________________________________________________________________________
Double_t StFJWrapper::GetFilteredJetArea(UInt_t idx) const
{
  // Get the filtered jet area.
  Double_t retval = -1; // really wrong area..
  if (fDoFilterArea && fClustSeqActGhosts && (idx<fFilteredJets.size())) {
    retval = fClustSeqActGhosts->area(fFilteredJets[idx]);
  } else {
    __ERROR(Form("Wrong index: %d",idx));
  }
  return retval;
}

//_________________________________________________________________________________________________
fastjet::PseudoJet StFJWrapper::GetJetAreaVector(UInt_t idx) const
{
  // Get the jet area as vector.
  fastjet::PseudoJet retval;
  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->area_4vector(fInclusiveJets[idx]);
  } else {
    __ERROR(Form("Wrong index: %d",idx));
  }
  return retval;
}

//_________________________________________________________________________________________________
fastjet::PseudoJet StFJWrapper::GetFilteredJetAreaVector(UInt_t idx) const
{
  // Get the jet area as vector.
  fastjet::PseudoJet retval;
  if (fDoFilterArea && fClustSeqActGhosts && (idx<fFilteredJets.size())) {
    retval = fClustSeqActGhosts->area_4vector(fFilteredJets[idx]);
  } else {
    __ERROR(Form("Wrong index: %d",idx));
  }
  return retval;
}

//_________________________________________________________________________________________________
std::vector<double> StFJWrapper::GetSubtractedJetsPts(Double_t median_pt, Bool_t sorted)
{
  // Get subtracted jets pTs, returns vector.
  SubtractBackground(median_pt);

  if (kTRUE == sorted) {
    std::sort(fSubtractedJetsPt.begin(), fSubtractedJetsPt.begin());
  }
  return fSubtractedJetsPt;
}

//_________________________________________________________________________________________________
Double_t StFJWrapper::GetJetSubtractedPt(UInt_t idx) const
{
  // Get subtracted jets pTs, returns Double_t.
  Double_t retval = -99999.; // really wrong pt..
  if ( idx < fSubtractedJetsPt.size() ) {
    retval = fSubtractedJetsPt[idx];
  }
  return retval;
}

//_________________________________________________________________________________________________
std::vector<fastjet::PseudoJet>
StFJWrapper::GetJetConstituents(UInt_t idx) const
{
  // Get jets constituents.
  std::vector<fastjet::PseudoJet> retval;

  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->constituents(fInclusiveJets[idx]);
  } else {
    __ERROR(Form("Wrong index: %d",idx));
  }

  return retval;
}

//_________________________________________________________________________________________________
std::vector<fastjet::PseudoJet>
StFJWrapper::GetFilteredJetConstituents(UInt_t idx) const
{
  // Get jets constituents.
  std::vector<fastjet::PseudoJet> retval;

  if ( idx < fFilteredJets.size() ) {
    if (fClustSeqSA)        retval = fClustSeqSA->constituents(fFilteredJets[idx]);
    if (fClustSeqActGhosts) retval = fClustSeqActGhosts->constituents(fFilteredJets[idx]);
  } else {
    __ERROR(Form("Wrong index: %d",idx));
  }

  return retval;
}

//_________________________________________________________________________________________________
void StFJWrapper::GetMedianAndSigma(Double_t &median, Double_t &sigma, Int_t remove) const
{
  // Get the median and sigma from fastjet.
  // User can also do it on his own because the cluster sequence is exposed (via a getter)
  if (!fClustSeq) {
    __ERROR(Form("Run the jfinder first."));
    return;
  }

  Double_t mean_area = 0;
  try {
    if(0 == remove) {
      fClustSeq->get_median_rho_and_sigma(*fRange, fUseArea4Vector, median, sigma, mean_area);
    }  else {
      std::vector<fastjet::PseudoJet> input_jets = sorted_by_pt(fClustSeq->inclusive_jets());
      input_jets.erase(input_jets.begin(), input_jets.begin() + remove);
      fClustSeq->get_median_rho_and_sigma(input_jets, *fRange, fUseArea4Vector, median, sigma, mean_area);
      input_jets.clear();
    }
  } catch (fj::Error) {
    __WARNING(Form("FJ Exception caught."));
    median = -1.;
    sigma = -1;
  }
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::Run()
{
  // Run the actual jet finder.
  if (fAreaType == fj::voronoi_area) {
    // Rfact - check dependence - default is 1.
    // NOTE: hardcoded variable!
    fVorAreaSpec = new fj::VoronoiAreaSpec(1.);
    fAreaDef     = new fj::AreaDefinition(*fVorAreaSpec);
  } else {
    fGhostedAreaSpec = new fj::GhostedAreaSpec(fMaxRap,
                                               fNGhostRepeats,
                                               fGhostArea,
                                               fGridScatter,
                                               fKtScatter,
                                               fMeanGhostKt);

    fAreaDef = new fj::AreaDefinition(*fGhostedAreaSpec, fAreaType);
  }

  // this is acceptable by fastjet:
#ifndef FASTJET_VERSION
  fRange = new fj::RangeDefinition(fMaxRap - 0.95 * fR);
#else
  fRange = new fj::Selector(fj::SelectorAbsRapMax(fMaxRap - 0.95 * fR));
#endif

  if (fAlgor == fj::plugin_algorithm) {
    if (fPluginAlgor == 0) {
      // SIS CONE ALGOR
      // NOTE: hardcoded split parameter
      Double_t overlap_threshold = 0.75; // NOTE: this actually splits a lot: thr/min(pt1,pt2)
      fPlugin = new fj::SISConePlugin(fR,
                                      overlap_threshold,
                                      0,    //search of stable cones - zero = until no more
                                      1.0); // this should be seed effectively for proto jets
      fJetDef = new fastjet::JetDefinition(fPlugin);
    } else if (fPluginAlgor == 1) {
      // CDF cone
      // NOTE: hardcoded split parameter
      Double_t overlap_threshold = 0.75; // NOTE: this actually splits a lot: thr/min(pt1,pt2)
      fPlugin = new fj::CDFMidPointPlugin(fR,
                                      overlap_threshold,
                                      1.0,    //search of stable cones - zero = until no more
                                      1.0); // this should be seed effectively for proto jets
      fJetDef = new fastjet::JetDefinition(fPlugin);
    } else {
     __ERROR(Form("Unrecognized plugin number!"));
    }
  } else {
    fJetDef = new fj::JetDefinition(fAlgor, fR, fScheme, fStrategy);
  }

  try {
    fClustSeq = new fj::ClusterSequenceArea(fInputVectors, *fJetDef, *fAreaDef);
  } catch (fj::Error) {
    __ERROR(Form("FJ Exception caught."));
    return -1;
  }

  // FJ3 :: Define an JetMedianBackgroundEstimator just in case it will be used
#ifdef FASTJET_VERSION
  fBkrdEstimator     = new fj::JetMedianBackgroundEstimator(fj::SelectorAbsRapMax(fMaxRap));
#endif

  if (fLegacyMode) { SetLegacyFJ(); } // for FJ 2.x even if fLegacyMode is set, SetLegacyFJ is dummy

  // cout << "============================== Input Vectors =================================" << endl; 
  // for (int i = 0; i < fInputVectors.size(); i++){
  //   cout << Form("%i\t%.3f\t%.3f\t%.3f\t%.3f", fInputVectors[i].user_index(), fInputVectors[i].px(), fInputVectors[i].py(), fInputVectors[i].pz(), fInputVectors[i].e()) << endl;
  //   // cout << fInputVectors[i].user_index() << "\t" << fInputVectors[i].px() << "\t" << fInputVectors[i].py() << "\t" << fInputVectors[i].pz() << "\t" << fInputVectors[i].e() << endl;
  // }
  // cout << "============================== End of Input =================================" << endl; 

  // cout << "Jet Def" << fJetDef->description() << endl;
  // cout << "Area Def" << fAreaDef->description() << endl;

  // inclusive jets:
  fInclusiveJets.clear();
  fInclusiveJets = fClustSeq->inclusive_jets(0.0);

  return 0;
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::Filter()
{
//
//  StFJWrapper::Filter
//
  fJetDef = new fj::JetDefinition(fAlgor, fR, fScheme, fStrategy);

  if (fDoFilterArea) {
    if (fInputGhosts.size()>0) {
      try {
        fClustSeqActGhosts = new fj::ClusterSequenceActiveAreaExplicitGhosts(fInputVectors,
                                                                           *fJetDef,
                                                                            fInputGhosts,
                                                                            fGhostArea);
      } catch (fj::Error) {
        __WARNING(Form("FJ Exception caught."));
        return -1;
      }

      fFilteredJets.clear();
      fFilteredJets =  fClustSeqActGhosts->inclusive_jets(0.0);
    } else {
      return -1;
    }
  } else {
    try {
      fClustSeqSA = new fastjet::ClusterSequence(fInputVectors, *fJetDef);
    } catch (fj::Error) {
      __WARNING(Form("FJ Exception caught."));
      return -1;
    }

    fFilteredJets.clear();
    fFilteredJets = fClustSeqSA->inclusive_jets(0.0);
  }

  return 0;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetLegacyFJ()
{
  // This methods enable legacy behaviour (FastJet 2.x) when StROOT is compiled with FastJet 3.x
#ifdef FASTJET_VERSION
    std::cout << "WARNING! Setting FastJet in legacy mode" << std::endl;
    if (fGhostedAreaSpec) { fGhostedAreaSpec->set_fj2_placement(kTRUE); }
     if (fBkrdEstimator) {
      fBkrdEstimator->set_provide_fj2_sigma(kTRUE);
      fBkrdEstimator->set_use_area_4vector(kFALSE);
    }
#endif
}

//_________________________________________________________________________________________________
void StFJWrapper::SubtractBackground(Double_t median_pt)
{
  // Subtract the background (specify the value - see below the meaning).
  // Negative argument means the bg will be determined with the current algorithm
  // this is the default behavior. Zero is allowed
  // Note: user may set the switch for area4vector based subtraction.

  Double_t median    = 0;
  Double_t sigma     = 0;
  Double_t mean_area = 0;

  // clear the subtracted jet pt's vector<double>
  fSubtractedJetsPt.clear();

  // check what was specified (default is -1)
  if (median_pt < 0) {
    try {
      fClustSeq->get_median_rho_and_sigma(*fRange, fUseArea4Vector, median, sigma, mean_area);
    }

    catch (fj::Error) {
      __WARNING(Form("FJ Exception caught."));
      median = -9999.;
      sigma = -1;
      fMedUsedForBgSub = median;
      return;
    }
    fMedUsedForBgSub = median;
  } else {
    // we do not know the sigma in this case
    sigma = -1;
    if (0.0 == median_pt) {
      __WARNING(Form("Median specified for bg subtraction is ZERO: %f \n", median_pt ));
      fMedUsedForBgSub = 0.;
    } else {
      fMedUsedForBgSub = median_pt;
    }
  }

  // subtract:
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    if ( fUseArea4Vector ) {
      // subtract the background using the area4vector
      fj::PseudoJet area4v = fClustSeq->area_4vector(fInclusiveJets[i]);
      fj::PseudoJet jet_sub = fInclusiveJets[i] - area4v * fMedUsedForBgSub;
      fSubtractedJetsPt.push_back(jet_sub.perp()); // here we put only the pt of the jet - note: this can be negative
    } else {
      // subtract the background using scalars
      // fj::PseudoJet jet_sub = fInclusiveJets[i] - area * fMedUsedForBgSub_;
      Double_t area = fClustSeq->area(fInclusiveJets[i]);
      // standard subtraction
      Double_t pt_sub = fInclusiveJets[i].perp() - fMedUsedForBgSub * area;
      fSubtractedJetsPt.push_back(pt_sub); // here we put only the pt of the jet - note: this can be negative
    }
  }
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::DoConstituentSubtraction() {
  // Do constituent subtraction
#ifdef FASTJET_VERSION
  CreateConstituentSub();
  // fConstituentSubtractor->set_alpha(/* double alpha */);
  // fConstituentSubtractor->set_max_deltaR(/* double max_deltaR */);

  // clear constituent subtracted jets
  fConstituentSubtrJets.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::PseudoJet subtracted_jet(0.,0.,0.,0.);
    if(fInclusiveJets[i].perp()>0.)
      subtracted_jet = (*fConstituentSubtractor)(fInclusiveJets[i]);
    fConstituentSubtrJets.push_back(subtracted_jet);
  }
  if(fConstituentSubtractor) { delete fConstituentSubtractor; fConstituentSubtractor = NULL; }

#endif
  return 0;
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::DoSoftDrop() {
  // Do grooming
#ifdef FASTJET_VERSION
  CreateSoftDrop();

  // clear groomed jets
  fGroomedJets.clear();
  //fastjet::Subtractor fjsub (fBkrdEstimator);
  //fSoftDrop->set_subtractor(&fjsub);
  //fSoftDrop->set_input_jet_is_subtracted(false); //??
  
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::PseudoJet groomed_jet(0.,0.,0.,0.);
    if(fInclusiveJets[i].perp()>0.){
      groomed_jet = (*fSoftDrop)(fInclusiveJets[i]);
      if(groomed_jet!=0) fGroomedJets.push_back(groomed_jet);
    }
    
  }
  if(fSoftDrop) { delete fSoftDrop; fSoftDrop = NULL; }

#endif
  return 0;
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::CreateSoftDrop() {
  // Do grooming
  #ifdef FASTJET_VERSION
  if (fSoftDrop) { delete fSoftDrop; } // protect against memory leaks
  
  fSoftDrop   = new fj::contrib::SoftDrop(fBeta,fZcut);
  
  #endif
  return 0;
}

/*
//_________________________________________________________________________________________________
Int_t StFJWrapper::CreateGenSub() {
  // Do generic subtraction for jet mass
  #ifdef FASTJET_VERSION
  if (fGenSubtractor) { delete fGenSubtractor; } // protect against memory leaks

  if (fUseExternalBkg)
    { fGenSubtractor   = new fj::contrib::GenericSubtractor(fRho,fRhom); }
  else
    {
    fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);
    #if FASTJET_VERSION_NUMBER >= 30100
    fGenSubtractor->set_common_bge_for_rho_and_rhom(); // see contrib 1.020 GenericSubtractor.hh line 62
    #endif
    }

  #endif
  return 0;
}
*/

//_________________________________________________________________________________________________
Int_t StFJWrapper::CreateConstituentSub() {
  // Do generic subtraction for jet mass
  #ifdef FASTJET_VERSION
  if (fConstituentSubtractor) { delete fConstituentSubtractor; } // protect against memory leaks

  // see ConstituentSubtractor.hh signatures
  // ConstituentSubtractor(double rho, double rhom=0, double alpha=0, double maxDeltaR=-1)
  if (fUseExternalBkg) { fConstituentSubtractor = new fj::contrib::ConstituentSubtractor(fRho,fRhom); }
  else                 { fConstituentSubtractor = new fj::contrib::ConstituentSubtractor(fBkrdEstimator); }  // FIXME Nov15, 2018 commented back in

  #endif
  return 0;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupAlgorithmfromOpt(const char *option)
{
  // Setup algorithm from char.
  std::string opt(option);

  if (!opt.compare("kt"))                fAlgor    = fj::kt_algorithm;
  if (!opt.compare("antikt"))            fAlgor    = fj::antikt_algorithm;
  if (!opt.compare("cambridge"))         fAlgor    = fj::cambridge_algorithm;
  if (!opt.compare("genkt"))             fAlgor    = fj::genkt_algorithm;
  if (!opt.compare("cambridge_passive")) fAlgor    = fj::cambridge_for_passive_algorithm;
  if (!opt.compare("genkt_passive"))     fAlgor    = fj::genkt_for_passive_algorithm;
  if (!opt.compare("ee_kt"))             fAlgor    = fj::ee_kt_algorithm;
  if (!opt.compare("ee_genkt"))          fAlgor    = fj::ee_genkt_algorithm;
  if (!opt.compare("plugin"))            fAlgor    = fj::plugin_algorithm;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupAreaTypefromOpt(const char *option)
{
  // Setup area type from char.
  std::string opt(option);

  if (!opt.compare("active"))                      fAreaType = fj::active_area;
  if (!opt.compare("invalid"))                     fAreaType = fj::invalid_area;
  if (!opt.compare("active_area_explicit_ghosts")) fAreaType = fj::active_area_explicit_ghosts;
  if (!opt.compare("one_ghost_passive"))           fAreaType = fj::one_ghost_passive_area;
  if (!opt.compare("passive"))                     fAreaType = fj::passive_area;
  if (!opt.compare("voronoi"))                     fAreaType = fj::voronoi_area;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupSchemefromOpt(const char *option)
{
  //
  // setup scheme from char
  //
  std::string opt(option);

  if (!opt.compare("BIpt"))   fScheme   = fj::BIpt_scheme;
  if (!opt.compare("BIpt2"))  fScheme   = fj::BIpt2_scheme;
  if (!opt.compare("E"))      fScheme   = fj::E_scheme;
  if (!opt.compare("pt"))     fScheme   = fj::pt_scheme;
  if (!opt.compare("pt2"))    fScheme   = fj::pt2_scheme;
  if (!opt.compare("Et"))     fScheme   = fj::Et_scheme;
  if (!opt.compare("Et2"))    fScheme   = fj::Et2_scheme;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupStrategyfromOpt(const char *option)
{
  // Setup strategy from char.
  std::string opt(option);

  if (!opt.compare("Best"))            fStrategy = fj::Best;
  if (!opt.compare("N2MinHeapTiled"))  fStrategy = fj::N2MinHeapTiled;
  if (!opt.compare("N2Tiled"))         fStrategy = fj::N2Tiled;
  if (!opt.compare("N2PoorTiled"))     fStrategy = fj::N2PoorTiled;
  if (!opt.compare("N2Plain"))         fStrategy = fj::N2Plain;
  if (!opt.compare("N3Dumb"))          fStrategy = fj::N3Dumb;
  if (!opt.compare("NlnN"))            fStrategy = fj::NlnN;
  if (!opt.compare("NlnN3pi"))         fStrategy = fj::NlnN3pi;
  if (!opt.compare("NlnN4pi"))         fStrategy = fj::NlnN4pi;
  if (!opt.compare("NlnNCam4pi"))      fStrategy = fj::NlnNCam4pi;
  if (!opt.compare("NlnNCam2pi2R"))    fStrategy = fj::NlnNCam2pi2R;
  if (!opt.compare("NlnNCam"))         fStrategy = fj::NlnNCam;
  if (!opt.compare("plugin"))          fStrategy = fj::plugin_strategy;
}

//_______________________________________________________________________________________________ 
Double_t StFJWrapper::NSubjettiness(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Int_t Option){
  //Option 0=Nsubjettiness result, 1=opening angle between axes in Eta-Phi plane, 2=Distance between axes in Eta-Phi plane
  
  fJetDef = new fj::JetDefinition(fAlgor, fR*100, fScheme, fStrategy ); //the *2 is becasue of a handful of jets that end up missing a track for some reason.

  try {
    fClustSeqSA = new fastjet::ClusterSequence(fInputVectors, *fJetDef);
    // ClustSeqSA = new fastjet::ClusterSequenceArea(fInputVectors, *fJetDef, *fAreaDef);
  } catch (fj::Error) {
    __WARNING(Form("FJ Exception caught."));
    return -1;
  }
  fFilteredJets.clear();
  fFilteredJets = fClustSeqSA->inclusive_jets(fMinJetPt-0.1); //becasue this is < not <=
  Double_t Result=-1;
  std::vector<fastjet::PseudoJet> SubJet_Axes;
  fj::PseudoJet SubJet1_Axis;
  fj::PseudoJet SubJet2_Axis;
  if (Algorithm==0){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==1) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==2){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==3) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==4) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==5){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==6){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==7){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==8){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==9){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==10){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::MultiPass_Axes(100), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }

  SubJet1_Axis=SubJet_Axes[0];	
  Double_t SubJet1_Eta=SubJet1_Axis.pseudorapidity();
  Double_t SubJet2_Eta;
  Double_t SubJet1_Phi=SubJet1_Axis.phi();
  if(SubJet1_Phi < -1*TMath::Pi()) SubJet1_Phi += (2*TMath::Pi());
  else if (SubJet1_Phi > TMath::Pi()) SubJet1_Phi -= (2*TMath::Pi());
  Double_t SubJet2_Phi;
  Double_t DeltaPhi=-5;
  if (SubJet_Axes.size()>1){
    SubJet2_Axis=SubJet_Axes[1];
    SubJet2_Eta=SubJet2_Axis.pseudorapidity();
    SubJet2_Phi=SubJet2_Axis.phi();
    if(SubJet2_Phi < -1*TMath::Pi()) SubJet2_Phi += (2*TMath::Pi());
    else if (SubJet2_Phi > TMath::Pi()) SubJet2_Phi -= (2*TMath::Pi());
    DeltaPhi=SubJet1_Phi-SubJet2_Phi;
    if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
    else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  }
    
  if (Option==0) return Result;
  else if (Option==1 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else if (Option==2 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else return -2;
}

//_______________________________________________________________________________________________
Double32_t StFJWrapper::NSubjettinessDerivativeSub(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Double_t JetR, fastjet::PseudoJet jet, Int_t Option){ //For derivative subtraction

  Double_t Result=-1;
  std::vector<fastjet::PseudoJet> SubJet_Axes;
  fj::PseudoJet SubJet1_Axis;
  fj::PseudoJet SubJet2_Axis;
  if (Algorithm==0){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==1) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==2){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==3) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==4) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==5){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==6){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==7){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==8){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==9){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==10){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::MultiPass_Axes(100), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }

  SubJet1_Axis=SubJet_Axes[0];	
  Double_t SubJet1_Eta=SubJet1_Axis.pseudorapidity();
  Double_t SubJet2_Eta;
  Double_t SubJet1_Phi=SubJet1_Axis.phi();
  if(SubJet1_Phi < -1*TMath::Pi()) SubJet1_Phi += (2*TMath::Pi());
  else if (SubJet1_Phi > TMath::Pi()) SubJet1_Phi -= (2*TMath::Pi());
  Double_t SubJet2_Phi;
  Double_t DeltaPhi=-5;
  if (SubJet_Axes.size()>1){
    SubJet2_Axis=SubJet_Axes[1];
    SubJet2_Eta=SubJet2_Axis.pseudorapidity();
    SubJet2_Phi=SubJet2_Axis.phi();
    if(SubJet2_Phi < -1*TMath::Pi()) SubJet2_Phi += (2*TMath::Pi());
    else if (SubJet2_Phi > TMath::Pi()) SubJet2_Phi -= (2*TMath::Pi());
    DeltaPhi=SubJet1_Phi-SubJet2_Phi;
    if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
    else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  }
    
  if (Option==0) return Result;
  else if (Option==1 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else if (Option==2 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else return -2;

}
#endif
