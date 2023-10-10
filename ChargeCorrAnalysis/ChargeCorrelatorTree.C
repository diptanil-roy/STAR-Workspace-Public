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
bool doTEST = kFALSE;  // FIXME
bool doMBset = kTRUE;  // HFT dataset (MB triggers for Run14 AuAu)
Bool_t RejectBadRuns = kFALSE;
Bool_t doCentSelection = kFALSE;
Bool_t doBackgroundJets = kFALSE;
Bool_t doConstituentSubtraction = kFALSE; // Remember to turn off this switch if you don't want constituent subtraction
Bool_t d0D0Analysis = kTRUE;
Bool_t makesmallfile = kFALSE;

Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

// constants
const double pi = 1.0*TMath::Pi();


void LoadLibs()
{
  // load fastjet libraries 3.x
  //gSystem->Load("libCGAL"); - not installed 
  gSystem->Load("$FASTJETNEW/lib/libfastjet");
  gSystem->Load("$FASTJETNEW/lib/libsiscone");
  gSystem->Load("$FASTJETNEW/lib/libsiscone_spherical");
  gSystem->Load("$FASTJETNEW/lib/libfastjetplugins");
  gSystem->Load("$FASTJETNEW/lib/libfastjettools");
  gSystem->Load("$FASTJETNEW/lib/libfastjetcontribfragile");

  // add include path to use its functionality
  gSystem->AddIncludePath("-I/gpfs01/star/pwg/droy1/SOFTWARE/fastjet-install/include");

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


void ChargeCorrelatorTree(const Char_t *inputFile="Run14_P18ih_HPSS_15164046.list", const Char_t *outputFile="test2.root", Int_t nEv = 10, const Char_t *fOutJobappend=""){

  TStopwatch timer;

  timer.Start();

  Int_t nEvents = 4000;
  if(nEv > 100) nEvents = 100000000;

  // open and close output .root file (so it exist and can be updated by Analysis Tasks)
  TFile *fout = new TFile(outputFile, "RECREATE");
  fout->Close();

  // Making the JetTree

  TString JetTreeFileName = "JetTree_";
  JetTreeFileName += outputFile;

  TFile *JetTreeFile = new TFile(JetTreeFileName.Data(), "RECREATE");
  JetTreeFile->Close();

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

    // input file for tests (based on Run) - updated for new Runs as needed
    if((RunYear == mRun12) && doTEST) inputFile = "Run12_D0.list";
    // if((RunYear == mRun14) && doMBset && doTEST) inputFile = "/gpfs01/star/pwg/droy1/TestCodeDirectory/TestFileList.list";
    if((RunYear == mRun14) && doMBset && doTEST) inputFile = "Run14MB_TestFileList.list";
    if((RunYear == mRun14) && !doMBset && doTEST) inputFile = "Run14.list"; //Run_15151042_files.list"; //"testLIST_Run14.list";
    if((RunYear == mRun16) && doTEST) inputFile = "test_run17124003_files.list";
    if((RunYear == mRun17) && doTEST && dopp) inputFile = "filelist_pp2017.list";
    cout<<"inputFileName = "<<inputFile<<endl;

    // centrality global flags - no centrality for pp collisions
    Int_t CentralitySelection;
    Int_t CentralityDefinition;
    if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // no centrality defintion for pp, just set one 
    // if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi; // (NEW - from Nick Aug22, 2019: all lumi) this will be default
    //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // Run14 P18ih (NEW - from Nick June10, 2019)
    if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id;
    if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id;         // Run16 - option: StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx;
    cout<<"Centrality definition: "<<CentralityDefinition<<endl;

    // Run/Event Flag
    Int_t RunFlag;
    if(RunYear == mRun12 && dopp) RunFlag = StJetFrameworkPicoBase::Run12_pp200;
    if(RunYear == mRun14 && !doMBset) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
    if(RunYear == mRun14 &&  doMBset) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200_MB;
    if(RunYear == mRun16) RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
    if(RunYear == mRun17 && dopp) RunFlag = StJetFrameworkPicoBase::Run17_pp510;
    Bool_t RejectBadRuns = kFALSE; // switch to load and than omit bad runs
    Int_t fBadRunListVers = StJetFrameworkPicoBase::fBadRuns_w_missing_HT; // fBadRuns_w_missing_HT, fBadRuns_wo_missing_HT,

    // trigger flags - update default
    Int_t EmcTriggerEventType; // kIsHT1 or kIsHT2 or kIsHT3
    if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;
    if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;
    // if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kAny;
    if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
    if(RunYear == mRun17) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT3; 
    Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default
    if(RunYear == mRun12) MBEventType = StJetFrameworkPicoBase::kRun12main; // default for Run12 pp
    if(RunYear == mRun17) MBEventType = StJetFrameworkPicoBase::kVPDMB; // default for Run17 pp
    Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;  // kTriggerANY, kTriggerMB, kTriggerHT  - only used by JetMaker and EPMaker (set to HT when doing EP corrections)

    // track flags
    bool usePrimaryTracks;
    if(RunYear == mRun12) usePrimaryTracks = kTRUE;
    if(RunYear == mRun14) usePrimaryTracks = kFALSE;  // = kTRUE for Run14, kFALSE for Run16
    if(RunYear == mRun16) usePrimaryTracks = kFALSE; 
    if(RunYear == mRun17) usePrimaryTracks = kTRUE;  

    // jet type flags
    Int_t fJetType;
    if(RunYear == mRun12) fJetType = kFullJet;
    if(RunYear == mRun14) fJetType = kFullJet; // kChargedJet
    if(RunYear == mRun16) fJetType = kChargedJet;
    if(RunYear == mRun17) fJetType = kFullJet;
    double fJetRadius = 0.4;  // 0.4, 0.3, 0.2
    double fJetConstituentCut = 2.0; // correlation analysis: 2.0, jet shape analysis: 2.0 (been using 2.0 for corrections)
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

    Int_t CentralityBins[centralitybins] = {-1, 5, 10, 20, 30, 50, 80};
    string CentralityBinsName[centralitybins] = {"00", "05", "10", "20", "30", "50", "80"};

  // create the picoMaker maker
  //StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
    StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst");
    // picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));
    picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::VpdOrDefault));

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

    // create JetFinder first (JetMaker)
    // task Name, track constituent cut, doHistos, output file name
    StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", fJetConstituentCut, kTRUE, outputFile);
    jetTask->SetJetType(fJetType);          // jetType
    jetTask->SetJetAlgo(antikt_algorithm);  // jetAlgo
    jetTask->SetRecombScheme(E_scheme); // recombination scheme
    //jetTask->SetRecombScheme(E_scheme); // recomb - this scheme actually doesn't pre-process the 4-vectors during the recombination scheme to set mass to 0 - USED for jet mass
    jetTask->SetRadius(fJetRadius);         // jet radius
    jetTask->SetJetsName("Jets");
    jetTask->SetMinJetPt(5.0);             // min signal jet pt to save to fJets array
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
    jetTask->SetdoConstituentSubtr(kFALSE);      // implement constituent subtractor, if TRUE, don't want to subtract underlying event Rho
    jetTask->SetRejectBadRuns(RejectBadRuns);    // switch to load and than omit bad runs


    StLeadJetChargeCorrelator *jetchargecorr = new StLeadJetChargeCorrelator(Form("JetTree"), picoMaker, JetTreeFileName.Data(), "JetMaker", "");
    cout << jetchargecorr->GetName() << endl;
    

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
    if (JetTreeFile->IsOpen()) JetTreeFile->Close();

    timer.Stop();
    cout << "Real Time Used: " << timer.RealTime()/60 << "m" << endl;
}
