#ifndef StJetShapeMTPMaker_h
#define StJetShapeMTPMaker_h

#include "StJetFrameworkPicoBase.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

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

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

// jet-framework classes
class StCentMaker;
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;
class StEventPoolManager;
class StEventPool;
class StCentMaker;

class StJetShapeMTPMaker : public StJetFrameworkPicoBase {
  public:

    StJetShapeMTPMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char *jetMakerName, const char *rhoMakerName);
    virtual ~StJetShapeMTPMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();
    
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp;}
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    //virtual void	    SetdoConstituentSubtr(Bool_t c)    { doBGSubs = c;}
    //virtual void            SetSystematicUncType(Int_t a)      { fSysUncType = a; }
    
    // jet setters
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetConstituentCut(Double_t mc)  { fJetConstituentCut= mc;}    // min constituent pt cut
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetMaxTowerEt(Double_t t)        { fTowerBias        = t; }    // tower bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }
    virtual void            SetBadRunListVers(Int_t i)         { fBadRunListVers = i; }
    
    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }
    // tower setters
    virtual void            SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
    virtual void            SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
    virtual void            SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }
    // event selection - setters    virtual void            SetEmcTriggerEventType(UInt_t te)    { fEmcTriggerEventType = te;  }
    virtual void            SetMBEventType(UInt_t mbe)           { fMBEventType = mbe; }
    // efficiency correction setter
    virtual void            SetDoEffCorr(Bool_t effcorr)         { fDoEffCorr = effcorr; }
    //virtual void            SetTrackEfficiencyType(Int_t t)      { fTrackEfficiencyType = t; }
    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }

    // event mixing - setters
    virtual void            SetEventMixing(Int_t yesno)	       { fDoEventMixing=yesno; }
    virtual void            SetMixingTracks(Int_t tracks)      { fMixingTracks = tracks; }
    virtual void            SetNMixedTr(Int_t nmt)             { fNMIXtracks = nmt; }
    virtual void            SetNMixedEvt(Int_t nme)            { fNMIXevents = nme; }
    virtual void            SetCentBinSize(Int_t centbins)     { fCentBinSize = centbins; }
    virtual void            SetCentBinSizeJS(Int_t centbins)   { fCentBinSizeJS = centbins; }
    virtual void            SetDoUseMultBins(Bool_t mult)      { fDoUseMultBins = mult; }
    virtual void            SetdoUseEPBins(Bool_t ep)          { doUseEPBins = ep; }
    virtual void            SetnEPBins(Int_t nep)              { fnEPBins = nep; }    
   
  protected:
    void                    RunTracks();
    void                    RunTowers();
    void                    RunJets(StEventPool *pool);
    void                    FillEmcTriggers();                          // EmcTrigger counter histo
    void                    SetSumw2(); // set errors weights 

    void    		    fillTrackQA(int bgMethod, vector<fastjet::PseudoJet> arr, int ptBin, int centBin);
    
    
    vector<fastjet::PseudoJet> GetEtaRefTracks(Double_t jetEta, Double_t jetPhi);
    
    double                  calculate_energy(double px, double py, double pz);
    double                  calculate_lesub(vector<fastjet::PseudoJet> arr);
    double                  calculate_ptd(vector<fastjet::PseudoJet> arr);
    double                  calculate_girth(vector<fastjet::PseudoJet> arr, double jet_pt, double jet_eta, double jet_phi);

    int                     getPtBin(float pt);
    int                     getCentBin(float cent);
    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Bool_t                  fDoEffCorr;              // efficiency correction to tracks

    // event mixing
    Int_t                   fDoEventMixing;          // switch ON/off event mixing
    Int_t                   fMixingTracks;           // MAX # of mixing tracks to keep in pool, before removing old to add new
    Int_t                   fNMIXtracks;             // MIN # of mixing track in pool before performing mixing
    Int_t                   fNMIXevents;             // MIN # of mixing events in pool before performing mixing
    Int_t                   fCentBinSize;            // centrality bin size of mixed event pools
    Int_t                   fCentBinSizeJS;          // centrality bin size of mixed event pools for jet shape analysis
    Bool_t                  fDoUseMultBins;          // use multiplicity bins instead of centrality bins - used for Jet Shape Analysis
    Bool_t                  doUseEPBins;             // use event plane bins: 0.2-2.0 GeV charged tracks
    Int_t                   fnEPBins;                // number of event plane bins to use for event mixing (0, pi) range

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA
        
    // track efficiency file
    TFile                  *fEfficiencyInputFile;

    // event pool
    TClonesArray           *CloneAndReduceTrackList();
    StEventPoolManager     *fPoolMgr;//!  // event pool Manager object

  private:
    // variables
    Int_t                   fRunNumber;
    Double_t		    refCorr2;
    
    //constants
    static const int nJetQATypes = 4;
    static const int nJetBGSMethods = 3;

    static const int nTrkQATypes = 3;
    static const int nTrkBGSMethods = 3;

    static const int nJetShapeTypes = 3;

    static const int nPtBins = 3;
    static const int nCentBins = 6;
    
    const double jetRadius = 0.4;
    const float fJetEtaCut = 0.6;
    
    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);
    
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
    
    
    
    string QANames[nJetQATypes] = {"Pt", "Eta", "Phi","Area"};
    //string BGSNames[nJetBGSMethods] = {"None", "Constit", "Eta Ref BG", "Eta Ref BG Subtracted"};
    string BGSNames[nJetBGSMethods] = {"None", "Constit", "Eta Ref BG"};
    string PtBinNames[nPtBins] = {"All", "15-25", ">25"};
    string CentBinNames[nCentBins] = {"All", "0-10", "10-20","20-50","50-80","80-100"};
    string JetShapeNames[nJetShapeTypes] = {"Lesub","Girth","PtD"};
    //0-10, 10-20, 20-50, 
    //15-20,20-40,10-15
   
    

    // Jet QA histos
    //types: pt, eta, phi, area
    //BGSubMethods: None, constit, eta ref, eta ref subbed
    //pt bins: all, 15-25, 25-inf
    //cent bins: all, 0-10, 10-30, 30-50, 50-70, 70-100
    //           [type][BGSubMethod][Pt Bin][Cent bin]
    TH1F           *fHistJetQA[nJetQATypes][nJetBGSMethods][nPtBins][nCentBins];

    // Track QA histos
    //types: pt, eta, phi
    //BGSubMethods: None, constit, eta ref
    //pt, cent bins: same as jet QA
    //           [type][BGSubMethod][Pt Bin][Cent bin]
    TH1F           *fHistTrackQA[nTrkQATypes][nTrkBGSMethods][nPtBins][nCentBins];

    // Jet Shape histos
    //types: Lesub, Girth, PtD
    //everything else same as jet QA
    TH1F           *fHistJetShapes[nJetShapeTypes][nJetBGSMethods][nPtBins][nCentBins];
    
    
    
   
    
    
    

    ClassDef(StJetShapeMTPMaker, 1)
};
#endif
