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
Bool_t doBackgroundJets = kTRUE;
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


void D0JetTreeMaker(const Char_t *inputFile="Run14_P18ih_HPSS_15164046.list", const Char_t *outputFile="test2.root", Int_t nEv = 10, const Char_t *fOutJobappend=""){

  TStopwatch timer;

  timer.Start();

  Int_t nEvents = 100;
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

    // input file for tests (based on Run) - updated for new Runs as needed
    if((RunYear == mRun12) && doTEST) inputFile = "Run12_D0.list";
    // if((RunYear == mRun14) && doMBset && doTEST) inputFile = "/gpfs01/star/pwg/droy1/TestCodeDirectory/TestFileList.list";
    if((RunYear == mRun14) && doMBset && doTEST) inputFile = "TestFileList.list";
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


    StTagD0Events *D0Tagger = new StTagD0Events("D0Tagger", picoMaker, outputFile);
    D0Tagger->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
    D0Tagger->SetTurnOnCentSelection(doCentSelection);    // run analysis for specific centrality: BOOLEAN
    D0Tagger->SetCentralityBinCut(CentralitySelection);   // specific centrality range to run: if above is FALSE, this doesn't matter
    D0Tagger->SetdoppAnalysis(dopp);
    D0Tagger->SetSignalRange(1.7, 2.10);
    D0Tagger->SetBackgroundRangeUL(1.7, 2.10);
    D0Tagger->SetBackgroundRangeLS(1.7, 2.10);
    D0Tagger->SetTopoCutsLevel(1);
    D0Tagger->DoNotFillQAHistograms();
    // D0Tagger->DoNotFillDaugHistograms();
    D0Tagger->SetDebugLevel(0);
    
    cout<<D0Tagger->GetName()<<endl;  // print name of class instance

    // Making D0 Jet Trees

    StD0EventsJetMaker *jetTask[3][2]; //3 sets of cuts, US and LS per cut
    StD0EventsJetMaker *jetTaskBG[3][2];
    StD0Rho            *rhoTask[3][2];
    StJetTreeMaker     *jettree[3][2];

    TString CutNames[3] = {"Standard", "Tight", "Loose"};
    TString SignPair[2] = {"UnlikeSign", "LikeSign"};

    // Making the JetTree

    TString JetTreeFileName = "JetTree_";
    JetTreeFileName += outputFile;

    TFile *JetTreeFile = new TFile(JetTreeFileName.Data(), "RECREATE");
    JetTreeFile->Close();

    for (int cut = 0; cut < 1; cut++){
      for (int bin = 0; bin < 2; bin++){
        double signpair;

        if (bin == 0) signpair = -1.;
        else if (bin == 1) signpair = 1.;

        double d0kind = signpair*(cut + 1);

        // cout << d0kind << endl;

        jetTask[cut][bin] = new StD0EventsJetMaker(Form("JetMaker_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()), fJetConstituentCut, kTRUE, outputFile);
        jetTask[cut][bin]->SetJetType(fJetType);          // jetType
        jetTask[cut][bin]->SetJetAlgo(antikt_algorithm);  // jetAlgo
        jetTask[cut][bin]->SetRecombScheme(E_scheme); // recombination scheme
        jetTask[cut][bin]->SetRadius(fJetRadius);         // jet radius
        jetTask[cut][bin]->SetJetsName("Jets");
        jetTask[cut][bin]->SetMinJetPt(0.0);             // 15.0 signal jets
        jetTask[cut][bin]->SetMaxJetTrackPt(30.0);        // max track constituent
        jetTask[cut][bin]->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
        jetTask[cut][bin]->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
        jetTask[cut][bin]->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        jetTask[cut][bin]->SetGhostArea(0.005);           // ghost area
        jetTask[cut][bin]->SetMinJetArea(0.0);            // minimum jet area
        jetTask[cut][bin]->SetJetEtaRange(-1.0 , 1.0); // fiducial eta acceptance
        jetTask[cut][bin]->SetJetPhiRange(0,2.0*pi);      // phi acceptance
        jetTask[cut][bin]->SetUsePrimaryTracks(usePrimaryTracks);
        jetTask[cut][bin]->SetRunFlag(RunFlag);           // run flag      
        jetTask[cut][bin]->SetdoppAnalysis(dopp);         // pp switch
        jetTask[cut][bin]->SetdoConstituentSubtr(doConstituentSubtraction); 
        jetTask[cut][bin]->SetD0Analysis(d0D0Analysis);
        jetTask[cut][bin]->SetD0Kind(d0kind);
        cout << jetTask[cut][bin]->GetName() << endl;

        if(doBackgroundJets) {
          jetTaskBG[cut][bin] = new StD0EventsJetMaker(Form("JetMakerBG_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()), fJetConstituentCut, kTRUE, outputFile);
          jetTaskBG[cut][bin]->SetJetType(fJetType);          // jetType
          jetTaskBG[cut][bin]->SetJetAlgo(kt_algorithm);      // jetAlgo
          jetTaskBG[cut][bin]->SetRecombScheme(E_scheme); // recombination scheme
          jetTaskBG[cut][bin]->SetRadius(fJetRadius);         // jet radius
          jetTaskBG[cut][bin]->SetJetsName("JetsBG");
          jetTaskBG[cut][bin]->SetMinJetPt(0.0);
          jetTaskBG[cut][bin]->SetMaxJetTrackPt(30.0);
          jetTaskBG[cut][bin]->SetMinJetTowerE(fJetConstituentCut);     // inclusive: 0.2
          jetTaskBG[cut][bin]->SetHadronicCorrFrac(1.0);      // hadronic correlation fraction 0-1
          jetTaskBG[cut][bin]->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
          jetTaskBG[cut][bin]->SetGhostArea(0.005);           // ghost area
          jetTaskBG[cut][bin]->SetMinJetArea(0.0);            // minimum jet area
          jetTaskBG[cut][bin]->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
          jetTaskBG[cut][bin]->SetJetPhiRange(0, 2.0*pi);   // 0,pi
          jetTaskBG[cut][bin]->SetUsePrimaryTracks(usePrimaryTracks);
          jetTaskBG[cut][bin]->SetRunFlag(RunFlag);
          jetTaskBG[cut][bin]->SetdoppAnalysis(dopp);
          jetTaskBG[cut][bin]->SetdoConstituentSubtr(kFALSE);
          jetTaskBG[cut][bin]->SetD0Analysis(d0D0Analysis);
          jetTaskBG[cut][bin]->SetD0Kind(d0kind);
          cout << jetTaskBG[cut][bin]->GetName() << endl;
        }

        bool dohisto = kFALSE;

        if(doBackgroundJets) { rhoTask[cut][bin] = new StD0Rho(Form("Rho_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()), dohisto, outputFile, Form("JetMakerBG_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()));
        } else { rhoTask[cut][bin] = new StD0Rho(Form("Rho_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()), dohisto, outputFile, Form("JetMaker_%s_%s", CutNames[cut].Data(), SignPair[bin].Data())); }
        rhoTask[cut][bin]->SetExcludeLeadJets(2);
        rhoTask[cut][bin]->SetOutRhoName("OutRho");
        rhoTask[cut][bin]->SetRunFlag(RunFlag);
        rhoTask[cut][bin]->SetdoppAnalysis(dopp);
        rhoTask[cut][bin]->SetRejectBadRuns(RejectBadRuns);     // switch to load and than omit bad runs

        jettree[cut][bin] = new StJetTreeMaker(Form("JetTree_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()), picoMaker, JetTreeFileName.Data(), Form("JetMaker_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()), Form("Rho_%s_%s", CutNames[cut].Data(), SignPair[bin].Data()));
        jettree[cut][bin]->SetD0Kind(d0kind);
        cout << jettree[cut][bin]->GetName() << endl;
      }
    }

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
