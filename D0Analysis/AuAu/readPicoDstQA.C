// ************************************** //
// Sept28, 2017
// some current notes:
// Run16 AuAu 200 GeV:  P16ij on RCAS has no EmcalTriggers and no clusters
//     ::P16ij has triggers that need to be accessed a different way
//
//

#include <TSystem>

// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
class StRefMultCorr;
// jet-framework STAR classes
class StJetMakerTask;
class StJetFrameworkPicoBase;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();
Double_t ZVtxMin = -40.0;
Double_t ZVtxMax =  40.0;

bool doCentSelection = kFALSE; //kTRUE;
bool dopp = kTRUE; // FIXME kTRUE for pp data
int RunYear = 12;   // FIXME
// kTRUE for local tests, kFALSE for job submission
bool doMBset = kFALSE;  // HFT dataset (MB triggers for Run14 AuAu)
bool doLow   = kFALSE;  // low-luminosity subset
bool doTEST  = kFALSE;  //FIXME FIXME!!!! be aware before submission

StChain *chain;

//__________________________________________________________________________________________________________
void readPicoDstQA(const Char_t *inputFile="Run14_P18ih_SL20d_mid_testFiles.list", const Char_t *outputFile="towerQA.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="_this_is_a_test")
{
//        Int_t nEvents = 1000;
        Int_t nEvents = 10000;
//        Int_t nEvents = 50000;
//        Int_t nEvents = 100000;        
        if(nEv > 100) nEvents = 100000000;

        // run enumerators
        enum RunType_t {
          mJobSubmission = 0, // 0th spot - default unless set
          mRun07 =  7, mRun08 = 8,  mRun09 = 9,  mRun10 = 10, mRun11 = 11, mRun12 = 12, mRun13 = 13,
          mRun14 = 14, mRun15 = 15, mRun16 = 16, mRun17 = 17, mRun18 = 18, mRun19 = 19, mRun20 = 20
        };

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        // =============================================================================== //
        // =============================================================================== //
        // over-ride functions
        if(dopp) {
          doCentSelection = kFALSE;  // can't ask for a particular centrality if requesting pp collisions
        }

        // input file for tests (based on Run) - updated for new Runs as needed
        //if((RunYear == mRun12) && doTEST && dopp)    inputFile = "pp200_production_2012_P12id.SL21d_QA-small.list";  // NFS files missing 
        if((RunYear == mRun12) && doTEST && dopp)    inputFile = "2012pp_testFILES-Aug2022.list";        // test files for pp, updated Aug22
        if((RunYear == mRun14) && doTEST)            inputFile = "Run14_P18ih_SL20d_mid_testFiles.list"  // test files for mid-lumi subset of Run14 AuAu
        if((RunYear == mRun14) && doTEST && doMBset) inputFile = "Run14_P16id_SL18f_MB_test.list";       // test files for HF version with MB triggers
        if((RunYear == mRun14) && doTEST && doLow)   inputFile = "Run14AuAu_P18ih_SL20d_low_test.list";  // test files for low-lumi subset of Run14 AuAu
        //if((RunYear == mRun16) && doTEST)            inputFile = "";
        cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags - no centrality for pp collisions
        Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050;
        Int_t CentralityDefinition;
        if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // no centrality defintion for Run 12 pp 
        if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi; // (NEW - from Nick: Aug16, 2019 set for all lumi)
        //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // Run14 P18ih (NEW - from Nick June10, 2019)
        if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id;         // Run16: StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx;
        cout<<"Centrality definition: "<<CentralityDefinition<<endl;

        // Run/Event Flag
        Int_t RunFlag;
        if(RunYear == mRun14 && !doMBset)          RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
        if(RunYear == mRun14 &&  doMBset)          RunFlag = StJetFrameworkPicoBase::Run14_AuAu200_MB;
        if(RunYear == mRun14 && !doMBset && doLow) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200_low;
        if(RunYear == mRun16)                      RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
        if(RunYear == mRun12 && dopp)              RunFlag = StJetFrameworkPicoBase::Run12_pp200;
        Bool_t RejectBadRuns = kFALSE; // switch to load and than omit bad runs
        Int_t fBadRunListVers = StJetFrameworkPicoBase::fBadRuns_w_missing_HT; // fBadRuns_w_missing_HT, fBadRuns_wo_missing_HT,

        // trigger flags - update default
        Int_t EmcTriggerEventType; // kIsHT1 or kIsHT2 or kIsHT3
        if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;
        if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT2 Run14
        if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
        Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default
        if(RunYear == mRun12) MBEventType = StJetFrameworkPicoBase::kRun12main; // default for Run12 pp
        //Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerHT;   // kTriggerANY, kTriggerMB, kTriggerHT - only used by JetMaker and EPMaker (set to HT when doing EP corrections)
        Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;  // kTriggerANY, kTriggerMB, kTriggerHT

        // FIXME - be aware of which list is used! 
        // tower flags - lists to load for bad towers, see StJetFrameworkPicoBase and below
        Int_t TowerListToUse = 136; // doesn't matter for charged jets - Run14 136-122: jet-hadron, early jet shape - Pt dep lists set below
        if(dopp) TowerListToUse = 169;
        // see StJetFrameworkPicoBase:   9992000 - 2 GeV, 9991000 - 1 GeV, 9990200 - 0.2 GeV  (applicable currently for Run12 pp and Run14 AuAu)
        // if using generic list:
        Double_t fJetConstituentCut = 2.0;
        if(fJetConstituentCut == 2.0) TowerListToUse = 9992000;
        if(fJetConstituentCut == 1.0) TowerListToUse = 9991000;
        if(fJetConstituentCut == 0.2) TowerListToUse = 9990200;
        // Run12: 1 - Raghav's list, 102 - my initial list, 169 - new list
        // Run14: 136 - main list (updated version for AuAu 200 GeV Run14), 122 - past used list
        // Run14 P18ih: 999 (initial) 
        cout<<"TowerListUsed: "<<TowerListToUse<<endl;

        // track flags
        bool usePrimaryTracks;
        if(RunYear == mRun12) usePrimaryTracks = kTRUE;
        if(RunYear == mRun14) usePrimaryTracks = kTRUE;  
        if(RunYear == mRun16) usePrimaryTracks = kFALSE; // don't have primary tracks (at least mid 2018 with that current production)

        // update settings for new centrality definitions
        if(CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 ||
           CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30 ||
           CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi
        ) { ZVtxMin = -30.0; ZVtxMax = 30.0; }

        // =============================================================================== //
        // =============================================================================== //

        // open and close output .root file (so it exist and can be updated by Analysis Tasks)
        TFile *fout = new TFile(outputFile, "RECREATE");
        fout->Close();
	
        // create chain
        StChain* chain = new StChain();

	// create the picoMaker maker:  (PicoIoMode, inputFile, name="picoDst")
	// - Write PicoDst's: PicoIoMode::IoWrite -> StPicoDstMaker::IoWrite
	// - Read  PicoDst's: PicoIoMode::IoRead  -> StPicoDstMaker::IoRead
        StPicoDstMaker *picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile, "picoDst");
        picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));

        // create base class maker pointer
        StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
        baseMaker->SetRunFlag(RunFlag);                  // run flag (year)
        baseMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        baseMaker->SetBadRunListVers(fBadRunListVers);          // switch to select specific bad run version file
        baseMaker->SetBadTowerListVers(TowerListToUse);
        cout<<baseMaker->GetName()<<endl;  // print name of class instance

        // create centrality class maker pointer
        StCentMaker *CentMaker = new StCentMaker("CentMaker", picoMaker, outputFile, kFALSE);
        CentMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
        CentMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
        CentMaker->SetRunFlag(RunFlag);                         // Run Flag
        CentMaker->SetdoppAnalysis(dopp);                       // pp-analysis switch
        CentMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
        CentMaker->SetUseBBCCoincidenceRate(kFALSE);            // BBC or ZDC (default) rate used? - USE ZDC (keep kFALSE)
        CentMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        CentMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        cout<<CentMaker->GetName()<<endl;  // print name of class instance

        // =======================================================================================================
        // QA task - MB events
        StPicoTrackClusterQA *Task = new StPicoTrackClusterQA("TrackClusterQA", kTRUE, outputFile);
        Task->SetTrackPtRange(0.2, 30.0);
        Task->SetTrackPhiRange(0., 2.0*pi);
        Task->SetTrackEtaRange(-1.0, 1.0);
        Task->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task->SetTowerERange(0.2, 100.0);
        Task->SetUsePrimaryTracks(usePrimaryTracks);
        Task->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Task->SetRunFlag(RunFlag);                      // RunFlag
        Task->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task->SetDebugLevel(8);
        Task->SetHadronicCorrFrac(1.0);
        Task->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        Task->SetDoTowerQAforHT(kFALSE);
        Task->SetdoppAnalysis(dopp);
        Task->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs

        // =======================================================================================================
        // QA task - MB events, little to no cuts
        StPicoTrackClusterQA *TasknoCuts = new StPicoTrackClusterQA("TrackClusterQAnocuts", kTRUE, outputFile);
        TasknoCuts->SetTrackPtRange(0.2, 30.0);
        TasknoCuts->SetTrackPhiRange(0., 2.0*pi);
        TasknoCuts->SetTrackEtaRange(-1.0, 1.0);
        TasknoCuts->SetEventZVtxRange(-40., 40.);  // TEST
        TasknoCuts->SetTowerERange(0.2, 100.0);
        TasknoCuts->SetUsePrimaryTracks(usePrimaryTracks);
        TasknoCuts->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        TasknoCuts->SetRunFlag(RunFlag);                      // RunFlag
        TasknoCuts->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        TasknoCuts->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        TasknoCuts->SetHadronicCorrFrac(1.0);
        TasknoCuts->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        TasknoCuts->SetDoTowerQAforHT(kFALSE);
        TasknoCuts->SetdoppAnalysis(dopp);
        TasknoCuts->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        TasknoCuts->SetDebugLevel(99); // TEST
        TasknoCuts->SetMaxEventTrackPt(100.0); // TEST

        // ======================================================================================================= 
        // QA task - HT events - defined by EmcTriggerEventType above
        StPicoTrackClusterQA *TaskQAHT = new StPicoTrackClusterQA("TrackClusterQAHT", kTRUE, outputFile);
        TaskQAHT->SetTrackPtRange(0.2, 30.0);
        TaskQAHT->SetTrackPhiRange(0., 2.0*pi);
        TaskQAHT->SetTrackEtaRange(-1.0, 1.0);
        TaskQAHT->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        TaskQAHT->SetTowerERange(0.2, 100.0);
        TaskQAHT->SetUsePrimaryTracks(usePrimaryTracks);
        TaskQAHT->SetEmcTriggerEventType(EmcTriggerEventType);    // kIsHT1 or kIsHT2 or kIsHT3
        TaskQAHT->SetRunFlag(RunFlag);                      // RunFlag
        TaskQAHT->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        TaskQAHT->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        TaskQAHT->SetHadronicCorrFrac(1.0);
        TaskQAHT->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        TaskQAHT->SetDoTowerQAforHT(kTRUE);   // kFALSE FIXME this was kFALSE on Apr 2021
        TaskQAHT->SetdoppAnalysis(dopp);
        TaskQAHT->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        ////////////TaskQAHT->SetTriggerToUse(trg2use); 

        // =======================================================================================================
        // QA task - compare with Hanseul, MB | HT2 | HT3
        const char *makerName = (dopp) ? "TrackClusterQAMBHT2HT3" : "TrackClusterQAMB30HT2HT3";
        Int_t trg2use = StJetFrameworkPicoBase::kTriggerMB30HT2HT3;
        if(dopp) trg2use = StJetFrameworkPicoBase::kTriggerMBHT2HT3;
        cout<<"maker name: "<<Form("%s", makerName)<<endl;

        StPicoTrackClusterQA *Task3 = new StPicoTrackClusterQA(Form("%s", makerName), kTRUE, outputFile);
        Task3->SetTrackPtRange(0.2, 30.0);
        Task3->SetTrackPhiRange(0., 2.0*pi);
        Task3->SetTrackEtaRange(-1.0, 1.0);
        Task3->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task3->SetTowerERange(0.2, 100.0);
        Task3->SetUsePrimaryTracks(usePrimaryTracks);
        Task3->SetEmcTriggerEventType(EmcTriggerEventType);    // kIsHT1 or kIsHT2 or kIsHT3
        Task3->SetRunFlag(RunFlag);                      // RunFlag
        Task3->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task3->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        Task3->SetHadronicCorrFrac(1.0);
        Task3->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        Task3->SetDoTowerQAforHT(kFALSE);
        Task3->SetdoppAnalysis(dopp);
        Task3->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task3->SetTriggerToUse(trg2use); 

/*
        // QA task - HT1
        StPicoTrackClusterQA *TaskHT1 = new StPicoTrackClusterQA("TrackClusterQAHT1", kTRUE, outputFile);
        TaskHT1->SetTrackPtRange(0.2, 30.0);
        TaskHT1->SetTrackPhiRange(0., 2.0*pi);
        TaskHT1->SetTrackEtaRange(-1.0, 1.0);
        TaskHT1->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        TaskHT1->SetTowerERange(0.2, 100.0);
        TaskHT1->SetUsePrimaryTracks(usePrimaryTracks);
        TaskHT1->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);    // kIsHT1 or kIsHT2 or kIsHT3
        TaskHT1->SetRunFlag(RunFlag);                      // RunFlag
        TaskHT1->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        TaskHT1->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        TaskHT1->SetHadronicCorrFrac(1.0);
        TaskHT1->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        TaskHT1->SetDoTowerQAforHT(kTRUE);
        TaskHT1->SetdoppAnalysis(dopp);
        TaskHT1->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs

        // QA task - HT2
        StPicoTrackClusterQA *TaskHT2 = new StPicoTrackClusterQA("TrackClusterQAHT2", kTRUE, outputFile);
        TaskHT2->SetTrackPtRange(0.2, 30.0);
        TaskHT2->SetTrackPhiRange(0., 2.0*pi);
        TaskHT2->SetTrackEtaRange(-1.0, 1.0);
        TaskHT2->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        TaskHT2->SetTowerERange(0.2, 100.0);
        TaskHT2->SetUsePrimaryTracks(usePrimaryTracks);
        TaskHT2->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT2);    // kIsHT1 or kIsHT2 or kIsHT3
        TaskHT2->SetRunFlag(RunFlag);                      // RunFlag
        TaskHT2->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        TaskHT2->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        TaskHT2->SetHadronicCorrFrac(1.0);
        TaskHT2->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        TaskHT2->SetDoTowerQAforHT(kTRUE);
        TaskHT2->SetdoppAnalysis(dopp);
        TaskHT2->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
*/

        // initialize chain
        chain->Init();
        cout<<"chain->Init();"<<endl;
        int total = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << total << endl;
        if(nEvents>total) nEvents = total;
  
        for (Int_t i=0; i<nEvents; i++){
          if(i%100==0) cout << "Working on eventNumber " << i << endl;

          chain->Clear();
          int iret = chain->Make(i);	
          if (iret) { cout << "Bad return code!" << iret << endl; break;}

          total++;		
	}
	
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	
	delete chain;	

        // close output file if open
        if(fout->IsOpen()) fout->Close();
}

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
  gSystem->AddIncludePath("-I$FASTJET/include");

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

void LoadMacros()
{
}
