#ifndef StPIDQA_h
#define StPIDQA_h

#include <iostream>
#include <fstream>
#include "StJetFrameworkPicoBase.h"
#include "StTrackTreeStruct.h"
#include "StD0TreeStruct.h"
#include "StDcaGeometry.h"
#include <TLorentzVector.h>
#include <vector>
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class TH3F;
class THn;
class THnSparse;
class TProfile;
class TString;
class TVector3;

// STAR classes
//class StarRoot;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoTrackCovMatrix;
class StDcaGeometry;
// class StPhysicalHelix;
class StRefMultCorr;
// class StarClassLibrary;

// jet-framework classes
class StCentMaker;
class StEmcPosition2;

class StPIDQA : public StJetFrameworkPicoBase {
  public:

    bool testrun = kTRUE;

    StPIDQA(const char* name, StPicoDstMaker *picoMaker, const char* outName);
    virtual ~StPIDQA();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    void    MakeD0s();

    // booking of histograms (optional)
    void    DeclareTree();
    void    BookTree();
    void    WriteTree();

    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    // particle identification
    Bool_t    IsAnAcceptableTrack(StPicoTrack *trk, bool dohistograms);
    Bool_t    IsTrackWithinJet(StJet *jet, StPicoTrack *trk);

    double    InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2);
    bool      IsPion(StPicoTrack *trk);
    bool      IsKaon(StPicoTrack *trk);

    void      FillPidHistograms(StPicoTrack *trk);


    float     GetTofBeta(StPicoTrack *trk);
    StDcaGeometry DcaGeometry(int trackid);

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp;}
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    virtual void            SetTestRun(Bool_t trun)            { testrun = trun; }

    virtual void            SetMCAnalysisWithoutCent()         { fMCEventsWithoutCent = kTRUE   ; }

    virtual void            SaveWideMassBin()                  { fD0MassWindow = kFALSE; } // Saves all KPi pairs instead of the ones between 1.5 and 2.2
    virtual void            SaveD0Trees()                      { fSaveD0Trees = kTRUE; } 
    virtual void            ProcessHTEvents()                  { fMBEvents = kFALSE; fHTEvents = kTRUE;} // Processs HT and JP events     
    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

    // topological cut setter

    virtual void            SetTopoCutsLevel(int l)             {fTopoLevel = l         ; }

    // tower setters
    virtual void            SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
    virtual void            SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
    virtual void            SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)    { fEmcTriggerEventType = te;  }
    virtual void            SetMBEventType(UInt_t mbe)           { fMBEventType = mbe; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Bool_t effcorr)          { fDoEffCorr = effcorr; }

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }

    virtual void            DoNotFillQAHistograms()                    { fdoQA_Histograms = kFALSE; }
    virtual void            DoNotFillDaugHistograms()                  { fdoDaug_Histograms = kFALSE; }

    virtual void            DoSingleParticleEmbedding()                { fdoSingleParticleEmbedding = kTRUE; }

    // mass band
    virtual void            SetSignalRange(Double_t lmass, Double_t umass)  {fInvMassSignal1 = lmass; fInvMassSignal2 = umass;}
    virtual void            SetBackgroundRangeUL(Double_t lmassUL, Double_t umassUL)    {fInvMassULBg1 = lmassUL; fInvMassULBg2 = umassUL;}
    virtual void            SetBackgroundRangeLS(Double_t lmassLS, Double_t umassLS)    {fInvMassLSBg1 = lmassLS; fInvMassLSBg2 = umassLS;}

    //Centrality bin setter
    Int_t                   GetFiveCentBin(Double_t centbin) const;
    Int_t                   GetD0PtBin(Double_t d0bin) const;

    // Standard D0 event getters
    Bool_t                  DoesEventHaveD0()                {return (fd0TrackIndices.size()!=0);}

    vector<TLorentzVector>           GetEmbeddedParticle()               {return fSingleParticleVector;}

    std::vector<std::vector<int>>    GetD0Indices()                      {return fd0TrackIndices;}


  protected:
    void                    RunTracks();
    void                    RunTowers();
    void                    TestTracks();
    void                    RunJets();
    Double_t                RelativePhi(Double_t mphi, Double_t vphi) const;      // relative jet track angle
    Double_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    void                    FillEmcTriggers();                          // EmcTrigger counter histo
    Bool_t                  DoComparison(int myarr[], int elems);
    void                    SetSumw2(); // set errors weights 

    // Topo Cut level
    int                     fTopoLevel;             // Topo cut level (0 is null, 1 is tight, 2 is less tight, 3 is lax)

    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Bool_t                  fDoEffCorr;              // efficiency correction to tracks

    Bool_t                  fD0MassWindow;
    Bool_t                  fSaveD0Trees;

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    Bool_t                  fMCEventsWithoutCent;
    Bool_t                  fdoQA_Histograms;
    Bool_t                  fdoDaug_Histograms;

    Bool_t                  fdoSingleParticleEmbedding;


    Bool_t                  fd0;
    Bool_t                  fd0Bg;

    Bool_t                  fTightd0;
    Bool_t                  fTightd0Bg;

    Bool_t                  fLoosed0;
    Bool_t                  fLoosed0Bg;

    Bool_t                  fMBEvents;
    Bool_t                  fHTEvents;

    TLorentzVector          fSingleParticle;
    vector<TLorentzVector>  fSingleParticleVector;

    // Saves out indices of PiK pairs and the mass (The first is always a pion, the second always a kaon (even if the same particle can be both pion and kaon))
    std::vector<std::vector<int> >  fd0TrackIndices; 

  private:

    // funcs
    // function to calculate relative phi between 2 objects and shift between 0 and 2pi 
    //___________________________________________________________________________________________
    Double_t dPhi(Double_t phi1, Double_t phi2) {
      Double_t deltaPhi;
      deltaPhi = abs(phi1 - phi2); //TODO absolute values
      if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
      if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

      if (deltaPhi > TMath::Pi()) deltaPhi=2*(TMath::Pi()) - deltaPhi;
      return deltaPhi;   // dphi in [0, 2Pi]
    }

    // Standardise phi to 0 to 2 Pi
    Double_t standardPhi(Double_t phi){
        Double_t phi_standard = phi;
        if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
        if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
        return phi_standard;
    }

    // function to calculate relative eta between 2 objects
    //___________________________________________________________________________________________
    Double_t dEta(Double_t eta1, Double_t eta2) {
      Double_t deltaEta;
      deltaEta = eta1 - eta2;

      return deltaEta;
    }

    // function to calculate relative eta between 2 objects
    //___________________________________________________________________________________________
    Double_t dR(Double_t delphi, Double_t deleta) {
      Double_t dRad;
      dRad = TMath::Sqrt(pow(delphi,2) + pow(deleta,2));
    
      return dRad;
    }

    // variables
    Int_t                   fRunNumber;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);

    // position object
    StEmcPosition2         *mEmcPosition;

    StTrackTreeStruct       fTrackTree;
    StD0TreeStruct           fD0Tree;
    TTree                  *d0unlikesigntree;
    TTree                  *d0likesigntree;

    vector<int>             d0CandidateID;
    vector<int>             d0BgCandidateULID;
    vector<int>             d0BgCandidateLSID;

    std::vector<double>     fPionIndices; //Vector of Track Indices corresponding to pion
    std::vector<double>     fKaonIndices; //Vector of Track Indices corresponding to kaon

    const double Mpion = 0.139570;
    const double Mkaon = 0.493677;
    const double Mproton = 0.938272;

    const double R = 0.4;
    const double deltar = 0.05;
    const int numberofbins = 8;
    const int numberofptbins = 9;

    // Vertex Histograms

    TH1F *hEventZvertex_diff;
    TH2F *hEventVzvsVzvpd;
    TH2F *hEventVxvsVy;

    // Event Histograms with and without cuts 

    TH1F *hCentralityBeforeCuts;
    TH1F *hCentralityAfterCuts;
    TH1F *hCentralityWeightedBeforeCuts;
    TH1F *hCentralityWeightedAfterCuts;
    TH1F *hRefMultiplicity;
    TH1F *hgRefMultiplicity;

    TH1F *cuthistogram_event;

    // Track acceptances

    TH1F *hAcceptedPt[11];
    TH1F *hAcceptedEta[11];

    // PID histograms

    TH2F *z_pi;
    TH2F *z_ka;
    TH2F *normalised_invbetavpT_tof_pi;
    TH2F *normalised_invbetavpT_tof_ka;
    TH2F *dEdXvp;
    TH2F *invbetavp;

    TH1F *hD0Mass;

    double fInvMassSignal1;
    double fInvMassSignal2;

    double fInvMassULBg1;
    double fInvMassULBg2;

    double fInvMassLSBg1;
    double fInvMassLSBg2;

    Int_t                   numberofevents[14];
    Int_t                   numberoftracks[10];

    // bad and dead tower list
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list
    std::set<Int_t>        badRuns;

    // base class pointer object
    StJetFrameworkPicoBase *mBaseMaker;

    // maker names
    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;

    ClassDef(StPIDQA, 1)
};
#endif
