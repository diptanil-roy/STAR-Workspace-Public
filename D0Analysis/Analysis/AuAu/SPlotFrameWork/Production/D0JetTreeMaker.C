// #include "StHFAnalysis.h"
// #include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TStopwatch.h"

class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StRefMultCorr;
// my jet-framework STAR classes
class StJetFrameworkPicoBase;
class StTagD0Events;
class StJetMakerTask;
class StRho;
class StRhoBase;
class StMyAnalysisMaker;


StChain *chain;

Int_t RunYear = 14; 
bool dopp = kFALSE;   // FIXME double check before submitting!
bool doTEST = kTRUE;  // FIXME
bool doMBset = kTRUE;  // HFT dataset (MB triggers for Run14 AuAu)
Bool_t RejectBadRuns = kFALSE;
Bool_t doCentSelection = kFALSE;
Bool_t doBackgroundJets = kFALSE;
Bool_t doConstituentSubtraction = kFALSE; // Remember to turn off this switch if you don't want constituent subtraction
Bool_t d0D0Analysis = kTRUE;
Bool_t makesmallfile = kFALSE;

Double_t ZVtxMin = -30.0;
Double_t ZVtxMax = 30.0;

// constants
const double pi = 1.0*TMath::Pi();


void LoadLibs()
{
  // load fastjet libraries 3.x
  //gSystem->Load("libCGAL"); - not installed 
    // Set $LD_LIBRARY_PATH via spack
/*  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libfastjettools");
  gSystem->Load("libfastjetcontribfragile");
*/
    gSystem->Load("$FASTJET/lib/libfastjet");
    gSystem->Load("$FASTJET/lib/libsiscone");
    gSystem->Load("$FASTJET/lib/libsiscone_spherical");
    gSystem->Load("$FASTJET/lib/libfastjetplugins");
    gSystem->Load("$FASTJET/lib/libfastjettools");
    gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");
  // add include path to use its functionality
  gSystem->AddIncludePath("-I/star/u/droy1/Y2019/STAR/FastJet/fastjet-install/include");

  // load the system libraries - these were defaults
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // these are needed for new / additional classes
  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  // my libraries
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StMyAnalysisMaker");

  gSystem->ListLibraries();
}


void D0JetTreeMaker(const Char_t *inputFile="Run14_P18ih_HPSS_15164046.list", const Char_t *outputFile="test.root", Int_t nEv = 10, const Char_t *fOutJobappend=""){

  TStopwatch timer;

  timer.Start();

  Int_t nEvents = 10;
  if(nEv > 100) nEvents = 100000000;

  // open and close output .root file (so it exist and can be updated by Analysis Tasks)
  TFile *fout = new TFile(outputFile, "RECREATE");
  fout->Close();

  LoadLibs();

    // run enumerators
  enum RunType_t {
      mJobSubmission = 0, // 0th spot - default unless set
      mRun07 =  7, mRun08 = 8,  mRun09 = 9, mRun10 = 10,
      mRun11 = 11, mRun12 = 12, mRun13 = 13,
      mRun14 = 14, mRun15 = 15, mRun16 = 16, mRun17 = 17, mRun18 = 18
  };

    // set up Jet enumerators:
    // jet type
  enum EJetType_t {
      kFullJet,  // tracks + clusters
      kChargedJet,
      kNeutralJet
  };

    // jet algorithm
    enum EJetAlgo_t {
      kt_algorithm                    = 0, // background jets
      antikt_algorithm                = 1, // signal jets
      cambridge_algorithm             = 2,
      genkt_algorithm                 = 3,
      cambridge_for_passive_algorithm = 11,
      genkt_for_passive_algorithm     = 13,
      plugin_algorithm                = 99,
      undefined_jet_algorithm         = 999
    };

    // jet recombination scheme
    enum ERecoScheme_t {
      E_scheme        = 0,
      pt_scheme       = 1, pt2_scheme      = 2,
      Et_scheme       = 3, Et2_scheme      = 4,
      BIpt_scheme     = 5, BIpt2_scheme    = 6,
      WTA_pt_scheme   = 7, WTA_modp_scheme = 8,
      external_scheme = 99
    };

    // jet shape jet type
    enum EJetShapeJetType_t {
      kInclusiveJets,
      kLeadingJets,
      kSubLeadingJets
    };

    // =============================================================================== //
    // over-ride functions
    if(dopp) {
      doCentSelection = kFALSE;  // can't ask for a particular centrality if requesting pp collisions
      doBackgroundJets = kFALSE; // don't do a rho background calculation when using pp collisions   
    }

    if((RunYear == mRun14) && doMBset && doTEST) inputFile = "test.list";
    if((RunYear == mRun14) && !doMBset && doTEST) inputFile = "Run14.list"; //Run_15151042_files.list"; //"testLIST_Run14.list";
    cout<<"inputFileName = "<<inputFile<<endl;

    // centrality global flags - no centrality for pp collisions
    Int_t CentralitySelection;
    Int_t CentralityDefinition;
    if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi; // (NEW - from Nick Aug22, 2019: all lumi) this will be default
    cout<<"Centrality definition: "<<CentralityDefinition<<endl;

    // Run/Event Flag
    Int_t RunFlag;
    if(RunYear == mRun14 && !doMBset) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
    if(RunYear == mRun14 &&  doMBset) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200_MB;
    Bool_t RejectBadRuns = kFALSE; // switch to load and than omit bad runs
    Int_t fBadRunListVers = StJetFrameworkPicoBase::fBadRuns_w_missing_HT; // fBadRuns_w_missing_HT, fBadRuns_wo_missing_HT,

    // trigger flags - update default
    Int_t EmcTriggerEventType; // kIsHT1 or kIsHT2 or kIsHT3
    if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;

    Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default

    Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;  // kTriggerANY, kTriggerMB, kTriggerHT  - only used by JetMaker and EPMaker (set to HT when doing EP corrections)

    // track flags
    bool usePrimaryTracks;
    if(RunYear == mRun14) usePrimaryTracks = kFALSE;  // = kTRUE for Run14, kFALSE for Run16

    // jet type flags
    Int_t fJetType;
    if(RunYear == mRun14) fJetType = kFullJet; // kChargedJet

    double fJetRadius = 0.4;  // 0.4, 0.3, 0.2
    double fJetConstituentCut = 0.2; // correlation analysis: 2.0, jet shape analysis: 2.0 (been using 2.0 for corrections)
    Int_t fJetAnalysisJetType = kInclusiveJets;  // Jet analysis jet types - options: kInclusiveJets, kLeadingJets, kSubLeadingJets (need to set up in analysis when you want to use)
    cout<<"fJetType: "<<fJetType<<endl;

    // FIXME - be aware of which list is used! 
    // tower flags - lists to load for bad towers, see StJetFrameworkPicoBase and below
    Int_t TowerListToUse = 136; // doesn't matter for charged jets
    if(dopp) TowerListToUse = 169;
    // see StJetFrameworkPicoBase:   9992000 - 2 GeV, 9991000 - 1 GeV, 9990200 - 0.2 GeV  (applicable currently for Run12 pp and Run14 AuAu)
    if(fJetConstituentCut == 2.0) TowerListToUse = 9992000;
    if(fJetConstituentCut == 1.0) TowerListToUse = 9991000;
    if(fJetConstituentCut == 0.2) TowerListToUse = 9990200;
    // Run12: 1 - Raghav's list, 102 - my initial list, 169 - new list
    // Run14: 136 - main list (updated version for AuAu 200 GeV Run14), 122 - past used list
    // Run14 P18ih: 999 (initial) 
    cout<<"TowerListUsed: "<<TowerListToUse<<endl;

    // update settings for new centrality definitions - certain productions had settings for z-vertex < 30 when calculating centrality definitions, etc..
    if(CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 ||
       CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30 ||
       CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi
    ) { ZVtxMin = -30.0; ZVtxMax = 30.0; }

    // =============================================================================== //

    StChain* chain = new StChain();

    // create the analysis maker!
    bool doComments = kFALSE;

    const int centralitybins = 7;

    Int_t CentralityBins[centralitybins] = {-1, 5, 10, 20, 30, 40, 80};
    string CentralityBinsName[centralitybins] = {"00", "05", "10", "20", "30", "40", "80"};

  // create the picoMaker maker
  //StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
    StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst");
    picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));

    // create base class maker pointer
    StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
    baseMaker->SetRunFlag(RunFlag);                  // run flag (year)
    baseMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
    baseMaker->SetBadRunListVers(fBadRunListVers);          // switch to select specific bad run version file
    baseMaker->SetBadTowerListVers(TowerListToUse);
    cout<<baseMaker->GetName()<<endl;  // print name of class instance

    // create centrality class maker pointer
    StCentMaker *CentMaker = new StCentMaker("CentMaker", picoMaker, outputFile, doComments);
    CentMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
    CentMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
    CentMaker->SetRunFlag(RunFlag);                         // Run Flag
    CentMaker->SetdoppAnalysis(dopp);                       // pp-analysis switch
    CentMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
    CentMaker->SetUseBBCCoincidenceRate(kFALSE);            // BBC or ZDC (default) rate used?
    CentMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
    CentMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
    cout<<CentMaker->GetName()<<endl;  // print name of class instance


    StTagD0Events *D0Tagger = new StTagD0Events("D0Tagger", picoMaker, outputFile);
    D0Tagger->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
    D0Tagger->SetTurnOnCentSelection(doCentSelection);    // run analysis for specific centrality: BOOLEAN
    D0Tagger->SetCentralityBinCut(CentralitySelection);   // specific centrality range to run: if above is FALSE, this doesn't matter
    D0Tagger->SetdoppAnalysis(dopp);
    D0Tagger->SetSignalRange(1.8, 1.92);
    D0Tagger->SetBackgroundRangeUL(1.7, 2.02);
    D0Tagger->SetBackgroundRangeLS(1.7, 2.02);
    D0Tagger->SetTopoCutsLevel(0);
    D0Tagger->SetDebugLevel(0);
    
    cout<<D0Tagger->GetName()<<endl;  // print name of class instance

    cout << "Starting Signal+Bg Jet Task" << endl;

    StD0EventsJetMaker *jetTask = new StD0EventsJetMaker("JetMakerSigBg", fJetConstituentCut, kTRUE, outputFile);
    jetTask->SetJetType(fJetType);          // jetType
    jetTask->SetJetAlgo(antikt_algorithm);  // jetAlgo
    jetTask->SetRecombScheme(E_scheme); // recombination scheme
    jetTask->SetRadius(fJetRadius);         // jet radius
    jetTask->SetJetsName("Jets");
    jetTask->SetMinJetPt(5.0);             // 15.0 signal jets
    jetTask->SetMaxJetTrackPt(30.0);        // max track constituent
    jetTask->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
    jetTask->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
    jetTask->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
    jetTask->SetGhostArea(0.005);           // ghost area
    jetTask->SetMinJetArea(0.0);            // minimum jet area
    jetTask->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // fiducial eta acceptance
    jetTask->SetJetPhiRange(0,2.0*pi);      // phi acceptance
    jetTask->SetUsePrimaryTracks(usePrimaryTracks);
    jetTask->SetRunFlag(RunFlag);           // run flag      
    jetTask->SetdoppAnalysis(dopp);         // pp switch
    jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax);     // can be tighter for Run16 (-20,20)
    jetTask->SetTurnOnCentSelection(doCentSelection);
    jetTask->SetdoConstituentSubtr(doConstituentSubtraction); 
    jetTask->SetCentralityBinCut(CentralitySelection);
    jetTask->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
    jetTask->SetTriggerToUse(TriggerToUse);
    jetTask->SetRejectBadRuns(RejectBadRuns);    // switch to load and than omit bad runs
    jetTask->SetD0Analysis(d0D0Analysis);
    jetTask->SetD0Kind(0);
    cout << jetTask->GetName() << endl;

    StD0EventsJetMaker *jetTaskBG;
    if(doBackgroundJets) {
      //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
      // task Name, track constituent cut, doHistos, output file name
      jetTaskBG = new StD0EventsJetMaker("JetMakerBGSigBg", fJetConstituentCut, kTRUE, outputFile);
      jetTaskBG->SetJetType(fJetType);          // jetType
      jetTaskBG->SetJetAlgo(kt_algorithm);      // jetAlgo
      jetTaskBG->SetRecombScheme(E_scheme); // recombination scheme
      jetTaskBG->SetRadius(fJetRadius);         // jet radius
      jetTaskBG->SetJetsName("JetsBG");
      jetTaskBG->SetMinJetPt(0.0);
      jetTaskBG->SetMaxJetTrackPt(30.0);
      jetTaskBG->SetMinJetTowerE(fJetConstituentCut);     // inclusive: 0.2
      jetTaskBG->SetHadronicCorrFrac(1.0);      // hadronic correlation fraction 0-1
      jetTaskBG->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
      jetTaskBG->SetGhostArea(0.005);           // ghost area
      jetTaskBG->SetMinJetArea(0.0);            // minimum jet area
      jetTaskBG->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
      jetTaskBG->SetJetPhiRange(0, 2.0*pi);   // 0,pi
      jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
      jetTaskBG->SetRunFlag(RunFlag);
      jetTaskBG->SetdoppAnalysis(dopp);
      jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);        // can be tighter for Run16 (-20,20)
      jetTaskBG->SetTurnOnCentSelection(doCentSelection);
      jetTaskBG->SetCentralityBinCut(CentralitySelection);
      jetTaskBG->SetEmcTriggerEventType(EmcTriggerEventType);// kIsHT1 or kIsHT2 or kIsHT3
      jetTaskBG->SetTriggerToUse(TriggerToUse);
      jetTaskBG->SetdoConstituentSubtr(kFALSE);
      jetTaskBG->SetRejectBadRuns(RejectBadRuns);      // switch to load and than omit bad runs
      jetTaskBG->SetD0Analysis(d0D0Analysis);
      jetTaskBG->SetD0Kind(0);
      cout << jetTaskBG->GetName() << endl;
    }


    bool dohisto = kFALSE;
    // Rho task, and scale it up to include neutral constituents
    StD0Rho *rhoTask;
    // RhoMaker name, doHisto switch, output filename, Background JetMaker name
    if(doBackgroundJets) { rhoTask = new StD0Rho("StRhoSigBg_JetsBG", dohisto, outputFile, "JetMakerBGSigBg");
    } else { rhoTask = new StD0Rho("StRhoSigBg_JetsBG", dohisto, outputFile, "JetMakerSigBg"); }
    rhoTask->SetExcludeLeadJets(2);
    rhoTask->SetOutRhoName("OutRho");
    rhoTask->SetRunFlag(RunFlag);
    rhoTask->SetdoppAnalysis(dopp);
    rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
    rhoTask->SetTurnOnCentSelection(doCentSelection);
    rhoTask->SetCentralityBinCut(CentralitySelection);
    rhoTask->SetRejectBadRuns(RejectBadRuns);     // switch to load and than omit bad runs


    StJetTreeMaker *JetTreeSigBg = new StJetTreeMaker("JetTreeSigBg", picoMaker, outputFile, "JetMakerSigBg", "StRhoSigBg_JetsBG");
    JetTreeSigBg->SetD0Kind(0);
    cout << JetTreeSigBg->GetName() << endl;


    cout << "Starting Unlike Bg Jet Task" << endl;

    StD0EventsJetMaker *jetTaskUL = new StD0EventsJetMaker("JetMakerUnlike", fJetConstituentCut, kTRUE, outputFile);
    jetTaskUL->SetJetType(fJetType);          // jetType
    jetTaskUL->SetJetAlgo(antikt_algorithm);  // jetAlgo
    jetTaskUL->SetRecombScheme(E_scheme); // recombination scheme
    jetTaskUL->SetRadius(fJetRadius);         // jet radius
    jetTaskUL->SetJetsName("Jets");
    jetTaskUL->SetMinJetPt(5.0);             // 15.0 signal jets
    jetTaskUL->SetMaxJetTrackPt(30.0);        // max track constituent
    jetTaskUL->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
    jetTaskUL->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
    jetTaskUL->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
    jetTaskUL->SetGhostArea(0.005);           // ghost area
    jetTaskUL->SetMinJetArea(0.0);            // minimum jet area
    jetTaskUL->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // fiducial eta acceptance
    jetTaskUL->SetJetPhiRange(0,2.0*pi);      // phi acceptance
    jetTaskUL->SetUsePrimaryTracks(usePrimaryTracks);
    jetTaskUL->SetRunFlag(RunFlag);           // run flag      
    jetTaskUL->SetdoppAnalysis(dopp);         // pp switch
    jetTaskUL->SetEventZVtxRange(ZVtxMin, ZVtxMax);     // can be tighter for Run16 (-20,20)
    jetTaskUL->SetTurnOnCentSelection(doCentSelection);
    jetTaskUL->SetdoConstituentSubtr(doConstituentSubtraction); 
    jetTaskUL->SetCentralityBinCut(CentralitySelection);
    jetTaskUL->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
    jetTaskUL->SetTriggerToUse(TriggerToUse);
    jetTaskUL->SetRejectBadRuns(RejectBadRuns);    // switch to load and than omit bad runs
    jetTaskUL->SetD0Analysis(d0D0Analysis);
    jetTaskUL->SetD0Kind(1);
    cout << jetTaskUL->GetName() << endl;

    StD0EventsJetMaker *jetTaskULBG;
    if(doBackgroundJets) {
      //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
      // task Name, track constituent cut, doHistos, output file name
      jetTaskULBG = new StD0EventsJetMaker("JetMakerUnlikeBG", fJetConstituentCut, kTRUE, outputFile);
      jetTaskULBG->SetJetType(fJetType);          // jetType
      jetTaskULBG->SetJetAlgo(kt_algorithm);      // jetAlgo
      jetTaskULBG->SetRecombScheme(E_scheme); // recombination scheme
      jetTaskULBG->SetRadius(fJetRadius);         // jet radius
      jetTaskULBG->SetJetsName("JetsBG");
      jetTaskULBG->SetMinJetPt(0.0);
      jetTaskULBG->SetMaxJetTrackPt(30.0);
      jetTaskULBG->SetMinJetTowerE(fJetConstituentCut);     // inclusive: 0.2
      jetTaskULBG->SetHadronicCorrFrac(1.0);      // hadronic correlation fraction 0-1
      jetTaskULBG->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
      jetTaskULBG->SetGhostArea(0.005);           // ghost area
      jetTaskULBG->SetMinJetArea(0.0);            // minimum jet area
      jetTaskULBG->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
      jetTaskULBG->SetJetPhiRange(0, 2.0*pi);   // 0,pi
      jetTaskULBG->SetUsePrimaryTracks(usePrimaryTracks);
      jetTaskULBG->SetRunFlag(RunFlag);
      jetTaskULBG->SetdoppAnalysis(dopp);
      jetTaskULBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);        // can be tighter for Run16 (-20,20)
      jetTaskULBG->SetTurnOnCentSelection(doCentSelection);
      jetTaskULBG->SetCentralityBinCut(CentralitySelection);
      jetTaskULBG->SetEmcTriggerEventType(EmcTriggerEventType);// kIsHT1 or kIsHT2 or kIsHT3
      jetTaskULBG->SetTriggerToUse(TriggerToUse);
      jetTaskULBG->SetdoConstituentSubtr(kFALSE);
      jetTaskULBG->SetRejectBadRuns(RejectBadRuns);      // switch to load and than omit bad runs
      jetTaskULBG->SetD0Analysis(d0D0Analysis);
      jetTaskULBG->SetD0Kind(1);
      cout << jetTaskULBG->GetName() << endl;
    }


    bool dohisto = kFALSE;
    // Rho task, and scale it up to include neutral constituents
    StD0Rho *rhoTaskUL;
    // RhoMaker name, doHisto switch, output filename, Background JetMaker name
    if(doBackgroundJets) { rhoTaskUL = new StD0Rho("StRhoUnlike_JetsBG", dohisto, outputFile, "JetMakerUnlikeBG");
    } else { rhoTaskUL = new StD0Rho("StRhoUnlike_JetsBG", dohisto, outputFile, "JetMakerUnlike"); }
    rhoTaskUL->SetExcludeLeadJets(2);
    rhoTaskUL->SetOutRhoName("OutRho");
    rhoTaskUL->SetRunFlag(RunFlag);
    rhoTaskUL->SetdoppAnalysis(dopp);
    rhoTaskUL->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
    rhoTaskUL->SetTurnOnCentSelection(doCentSelection);
    rhoTaskUL->SetCentralityBinCut(CentralitySelection);
    rhoTaskUL->SetRejectBadRuns(RejectBadRuns);     // switch to load and than omit bad runs


    // StD0AnalysisJetShape *D0JetShapeUL = new StD0AnalysisJetShape("D0JetShapeUnlike", picoMaker, outputFile, "JetMakerUnlike", "StRhoUnlike_JetsBG");
    // D0JetShapeUL->SetD0Kind(1);
    // cout << D0JetShapeUL->GetName() << endl;

    StJetTreeMaker *JetTreeULBg = new StJetTreeMaker("JetTreeULBg", picoMaker, outputFile, "JetMakerUnlike", "StRhoUnlike_JetsBG");
    JetTreeULBg->SetD0Kind(1);
    cout << JetTreeULBg->GetName() << endl;


    cout << "Starting Like Bg Jet Task" << endl;

    StD0EventsJetMaker *jetTaskLS = new StD0EventsJetMaker("JetMakerLike", fJetConstituentCut, kTRUE, outputFile);
    jetTaskLS->SetJetType(fJetType);          // jetType
    jetTaskLS->SetJetAlgo(antikt_algorithm);  // jetAlgo
    jetTaskLS->SetRecombScheme(E_scheme); // recombination scheme
    jetTaskLS->SetRadius(fJetRadius);         // jet radius
    jetTaskLS->SetJetsName("Jets");
    jetTaskLS->SetMinJetPt(5.0);             // 15.0 signal jets
    jetTaskLS->SetMaxJetTrackPt(30.0);        // max track constituent
    jetTaskLS->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
    jetTaskLS->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
    jetTaskLS->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
    jetTaskLS->SetGhostArea(0.005);           // ghost area
    jetTaskLS->SetMinJetArea(0.0);            // minimum jet area
    jetTaskLS->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // fiducial eta acceptance
    jetTaskLS->SetJetPhiRange(0,2.0*pi);      // phi acceptance
    jetTaskLS->SetUsePrimaryTracks(usePrimaryTracks);
    jetTaskLS->SetRunFlag(RunFlag);           // run flag      
    jetTaskLS->SetdoppAnalysis(dopp);         // pp switch
    jetTaskLS->SetEventZVtxRange(ZVtxMin, ZVtxMax);     // can be tighter for Run16 (-20,20)
    jetTaskLS->SetTurnOnCentSelection(doCentSelection);
    jetTaskLS->SetdoConstituentSubtr(doConstituentSubtraction); 
    jetTaskLS->SetCentralityBinCut(CentralitySelection);
    jetTaskLS->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
    jetTaskLS->SetTriggerToUse(TriggerToUse);
    jetTaskLS->SetRejectBadRuns(RejectBadRuns);    // switch to load and than omit bad runs
    jetTaskLS->SetD0Analysis(d0D0Analysis);
    jetTaskLS->SetD0Kind(2);
    cout << jetTaskLS->GetName() << endl;

    StD0EventsJetMaker *jetTaskLSBG;
    if(doBackgroundJets) {
      //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
      // task Name, track constituent cut, doHistos, output file name
      jetTaskLSBG = new StD0EventsJetMaker("JetMakerLikeBG", fJetConstituentCut, kTRUE, outputFile);
      jetTaskLSBG->SetJetType(fJetType);          // jetType
      jetTaskLSBG->SetJetAlgo(kt_algorithm);      // jetAlgo
      jetTaskLSBG->SetRecombScheme(E_scheme); // recombination scheme
      jetTaskLSBG->SetRadius(fJetRadius);         // jet radius
      jetTaskLSBG->SetJetsName("JetsBG");
      jetTaskLSBG->SetMinJetPt(0.0);
      jetTaskLSBG->SetMaxJetTrackPt(30.0);
      jetTaskLSBG->SetMinJetTowerE(fJetConstituentCut);     // inclusive: 0.2
      jetTaskLSBG->SetHadronicCorrFrac(1.0);      // hadronic correlation fraction 0-1
      jetTaskLSBG->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
      jetTaskLSBG->SetGhostArea(0.005);           // ghost area
      jetTaskLSBG->SetMinJetArea(0.0);            // minimum jet area
      jetTaskLSBG->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
      jetTaskLSBG->SetJetPhiRange(0, 2.0*pi);   // 0,pi
      jetTaskLSBG->SetUsePrimaryTracks(usePrimaryTracks);
      jetTaskLSBG->SetRunFlag(RunFlag);
      jetTaskLSBG->SetdoppAnalysis(dopp);
      jetTaskLSBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);        // can be tighter for Run16 (-20,20)
      jetTaskLSBG->SetTurnOnCentSelection(doCentSelection);
      jetTaskLSBG->SetCentralityBinCut(CentralitySelection);
      jetTaskLSBG->SetEmcTriggerEventType(EmcTriggerEventType);// kIsHT1 or kIsHT2 or kIsHT3
      jetTaskLSBG->SetTriggerToUse(TriggerToUse);
      jetTaskLSBG->SetdoConstituentSubtr(kFALSE);
      jetTaskLSBG->SetRejectBadRuns(RejectBadRuns);      // switch to load and than omit bad runs
      jetTaskLSBG->SetD0Analysis(d0D0Analysis);
      jetTaskLSBG->SetD0Kind(2);
      cout << jetTaskLSBG->GetName() << endl;
    }


    bool dohisto = kFALSE;
    // Rho task, and scale it up to include neutral constituents
    StD0Rho *rhoTaskLS;
    // RhoMaker name, doHisto switch, output filename, Background JetMaker name
    if(doBackgroundJets) { rhoTaskLS = new StD0Rho("StRhoLike_JetsBG", dohisto, outputFile, "JetMakerLikeBG");
    } else { rhoTaskLS = new StD0Rho("StRhoLike_JetsBG", dohisto, outputFile, "JetMakerLike"); }
    rhoTaskLS->SetExcludeLeadJets(2);
    rhoTaskLS->SetOutRhoName("OutRho");
    rhoTaskLS->SetRunFlag(RunFlag);
    rhoTaskLS->SetdoppAnalysis(dopp);
    rhoTaskLS->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
    rhoTaskLS->SetTurnOnCentSelection(doCentSelection);
    rhoTaskLS->SetCentralityBinCut(CentralitySelection);
    rhoTaskLS->SetRejectBadRuns(RejectBadRuns);     // switch to load and than omit bad runs


    StJetTreeMaker *JetTreeLSBg = new StJetTreeMaker("JetTreeLSBg", picoMaker, outputFile, "JetMakerLike", "StRhoLike_JetsBG");
    JetTreeLSBg->SetD0Kind(2);
    cout << JetTreeLSBg->GetName() << endl;


    // initialize chain
    chain->Init();
    cout<<"chain->Init();"<<endl;
    int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if(nEvents > total) nEvents = total;
    //nEvents = 20000;

    for (Int_t i = 0; i < nEvents; i++){
      if (doTEST) {if (i%1000 == 0) cout << "Working on eventNumber " << i << endl;}
      else {if (i%1000 == 0) cout << "Working on eventNumber " << i << endl;}
      // cout << "Working on eventNumber " << i << endl;

      chain->Clear();
      int iret = chain->Make(i);  
      if (iret) { cout << "Bad return code!" << iret << endl; break;}

      total++;  
    }

    chain->Finish();

    delete chain;

  // close output file if open
    if(fout->IsOpen())   fout->Close();
    if (makesmallfile && !doTEST){
      if(fsmalldstout->IsOpen()) fsmalldstout->Close();
    } 
    timer.Stop();
    cout << "Real Time Used: " << timer.RealTime()/60 << "m" << endl;
}
