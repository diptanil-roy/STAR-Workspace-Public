#ifndef StTagDMesonEvents_h
#define StTagDMesonEvents_h

#include <iostream>
#include <fstream>
#include "StJetFrameworkPicoBase.h"
#include <vector>
#include "StD0TreeStruct.h"
#include "StDStarTreeStruct.h"

class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class THnSparse;
class TProfile;
class TString;
class TVector3;
class TTree;

// STAR classes
//class StarRoot;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
// class StPicoTrackCovMatrix;
class StRefMultCorr;

// jet-framework classes
class StCentMaker;
class StEmcPosition2;

class StTagDMesonEvents : public StJetFrameworkPicoBase {
  public:

    bool testrun = kTRUE;

    StTagDMesonEvents(const char* name, StPicoDstMaker *picoMaker, const char* outName);
    virtual ~StTagDMesonEvents();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    void    DeclareD0Tree();
    void    BookD0Tree();
    void    WriteD0Tree();

    void    DeclareDStarTree();
    void    BookDStarTree();
    void    WriteDStarTree();

    // particle identification
    Bool_t    IsAnAcceptableTrack(StPicoTrack *trk, bool dohistograms);
    Bool_t    IsTrackWithinJet(StJet *jet, StPicoTrack *trk);
    void      IsWhatParticle(StPicoTrack *trk, int &pid, double &m, double &e);
    Int_t     IsWhatParticle(StPicoTrack *trk);
    void      InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2, int &particle, double &invmass, TVector3 &momentum);
    Double_t  InvariantMass(StPicoTrack *trk1, StPicoTrack *trk2);
    void      FillPidHistograms(StPicoTrack *trk);
    Int_t     TopologicalCuts(StPicoTrack *trk1, StPicoTrack *trk2, bool dohistograms);

    bool      IsSoftPion(StPicoTrack *trk);
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

    // mass band
    virtual void            SetSignalRange(Double_t lmass, Double_t umass)  {fInvMassSignal1 = lmass; fInvMassSignal2 = umass;}
    virtual void            SetBackgroundRangeUL(Double_t lmassUL, Double_t umassUL)    {fInvMassULBg1 = lmassUL; fInvMassULBg2 = umassUL;}
    virtual void            SetBackgroundRangeLS(Double_t lmassLS, Double_t umassLS)    {fInvMassLSBg1 = lmassLS; fInvMassLSBg2 = umassLS;}

    //Centrality bin setter
    Int_t                   GetFourCentBin(Double_t centbin) const;

    virtual void            SetDStarAnalysis()                 {fTagDStarEvents = kTRUE;}

    // event getters
    Bool_t                  DoesEventHaveD0()                  {return fd0;}
    Bool_t                  DoesEventHaveD0BgUS()              {return fd0BgUS;}
    Bool_t                  DoesEventHaveD0BgLS()              {return fd0BgLS;}

    Bool_t                  DoesEventHaveDStar()               {return fdstar;}
    Bool_t                  DoesEventHaveDStarLS()             {return fdstarLS;}

    std::vector<std::vector<double>> GetD0Indices()                   {return fd0TrackIndices;}
    std::vector<std::vector<double>> GetD0BgUSIndices()               {return fd0BgUSTrackIndices;}
    std::vector<std::vector<double>> GetD0BgLSIndices()               {return fd0BgLSTrackIndices;}

    std::vector<std::vector<double>> GetDStarIndices()                {return fdstarTrackIndices;}
    std::vector<std::vector<double>> GetDStarLSIndices()              {return fdstarLSTrackIndices;}

  protected:
    void                    GetD0Tracks();
    void                    GetDStarTracks();
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

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    Bool_t                  fMCEventsWithoutCent;

    Bool_t                  fd0;
    Bool_t                  fd0BgUS;
    Bool_t                  fd0BgLS;

    Bool_t                  fdstar;
    Bool_t                  fdstarLS;

    // Saves out indices of KPi pairs and the mass
    std::vector<std::vector<double> >  fd0TrackIndices; 
    std::vector<std::vector<double> >  fd0BgUSTrackIndices;
    std::vector<std::vector<double> >  fd0BgLSTrackIndices;

    // Saves out indices of KPiPi pairs and the mass
    std::vector<std::vector<double> >  fdstarTrackIndices; //Index of Pion Kaon sPion D0Mass DiffMass
    std::vector<std::vector<double> >  fdstarLSTrackIndices;

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

    Double_t Sigma(Double_t *x, Double_t *par) //Parameterisation For TOF Kaon Choosing
    {
       return 2.0*(par[0] + par[1]/pow((x[0] + par[2]), par[3]));
    }
    // variables
    Int_t                   fRunNumber;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);

    // position object
    StEmcPosition2         *mEmcPosition;

    vector<int>             d0CandidateID;
    vector<int>             d0BgCandidateULID;
    vector<int>             d0BgCandidateLSID;

    const double Mpion = 0.139570;
    const double Mkaon = 0.493677;
    const double Mproton = 0.938272;

    const double R = 0.4;
    const double deltar = 0.05;
    const int numberofbins = 8;
    const int numberofptbins = 9;

    Bool_t fTagDStarEvents;

    TH1F *hEventZvertex_whole;
    TH1F *hEventZvertex_VPD;
    TH1F *hEventZvertex_diff;

    TH1F *hCentrality;
    TH1F *hMultiplicity;
    // TH1F *hptdistro;
    // TH1F *hetadistro;
    TH1F *cuthistogram_event;

    TH1F *hAcceptedPt[10];
    TH1F *hAcceptedEta[10];

    TH2F *z_pi;
    TH2F *z_ka;
    TH2F *normalised_invbetavpT_tof_pi;
    TH2F *normalised_invbetavpT_tof_ka;

    TH2F *dEdXvp;
    TH2F *invbetavp;

    TH1F *decaylengthd0US;
    TH1F *distancepikUS;
    TH1F *distanced0PVUS;
    TH1F *dcakPVUS;
    TH1F *dcapiPVUS;

    TH1F *decaylengthd0LS;
    TH1F *distancepikLS;
    TH1F *distanced0PVLS;
    TH1F *dcakPVLS;
    TH1F *dcapiPVLS;


    TH1F *kaonpt;
    TH1F *pionpt;
    TH1F *d0pt;
    TH2F *kaonpionpt;
    TH1F *kaoneta;
    TH1F *pioneta;
    TH1F *d0eta;
    TH1F *kaonphi;
    TH1F *pionphi;
    TH1F *d0phi;

    TH1F *kaonbgUSpt;
    TH1F *pionbgUSpt;
    TH1F *d0bgUSpt;
    TH2F *kaonpionbgUSpt;
    TH1F *kaonbgUSeta;
    TH1F *pionbgUSeta;
    TH1F *d0bgUSeta;
    TH1F *kaonbgUSphi;
    TH1F *pionbgUSphi;
    TH1F *d0bgUSphi;

    TH1F *kaonbgLSpt;
    TH1F *pionbgLSpt;
    TH1F *d0bgLSpt;
    TH2F *kaonpionbgLSpt;
    TH1F *kaonbgLSeta;
    TH1F *pionbgLSeta;
    TH1F *d0bgLSeta;
    TH1F *kaonbgLSphi;
    TH1F *pionbgLSphi;
    TH1F *d0bgLSphi;


    TH1F *invmass;
    TH1F *invmassbg;

    TH2F *invmassvpt;
    TH2F *invmassbgvpt;

    TH1F *invmassDStar;
    TH1F *invmassbgDStar;

    TH2F *invmassDStarvpt;
    TH2F *invmassbgDStarvpt;

    TH1F *invmassdiff;
    TH1F *invmassdiffwrong;

    TH2F *invmassdiffDStarvpt;
    TH2F *invmassdiffwrongDStarvpt;

    TH1F *hNumberOfD0s;
    TH1F *hNumberOfD0BgUS;
    TH1F *hNumberOfD0BgLS;

    TH1F *hNumberOfDStar;
    TH1F *hNumberOfDStarBgLS;
    
    StD0TreeStruct fD0SigTree;
    StD0TreeStruct fD0BgTree;

    StDStarTreeStruct fDStarSigTree;
    StDStarTreeStruct fDStarBgTree;

    TTree *d0sigtree;
    TTree *d0bgtree;

    TTree *dstarsigtree;
    TTree *dstarbgtree;

    double frefCorr2 = -99;

    double fInvMassSignal1;
    double fInvMassSignal2;

    double fInvMassULBg1;
    double fInvMassULBg2;

    double fInvMassLSBg1;
    double fInvMassLSBg2;

    Int_t                   numberofevents[10];
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

    ClassDef(StTagDMesonEvents, 1)
};
#endif
