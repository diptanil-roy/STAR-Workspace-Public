// ************************************** //
// July 15, 2019
// some current notes:
// - set up as means of using StCentMaker to provide other classes with centrality information
//

#include <TSystem>

// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StRefMultCorr;
class StMyAnalysisMaker;

// my jet-framework STAR classes
class StJetFrameworkPicoBase;
// class StTagD0Events;
class StTagD0MCEvents;
// class StMCD0EventsJetMaker;
class StMCD0JetMaker;
class StSimD0EventsJetMaker;
class StMCRecoJetMatcher;
class StJetMakerTask;
class StRho;
class StRhoBase;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// find kt jets and perform rho subtraction - not needed for constituents 2.0+ GeV or certain analyses
bool doBackgroundJets = kTRUE;

// run analysis for specific centrality bin
bool doCentSelection = kFALSE; // keep false to run over all centralities
Bool_t doConstituentSubtraction = kFALSE;

Bool_t d0D0Analysis = kTRUE;

bool MCSim = kFALSE;


// data parameters
Int_t RunYear = 14; 
// kFALSE when submitting jobs, kTRUE for tests
bool doTEST = kTRUE;   // FIXME double check before submitting!
bool dopp = kFALSE;      // FIXME kTRUE for pp data

Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

// constants
const double pi = 1.0*TMath::Pi();

StChain *chain;

void TestMacro(const Char_t *inputFile="Run12.list", const Char_t *outputFile="Test", Int_t nEv = 1, const Char_t *fOutJobappend="")
{
       Int_t nEvents = 200;
        // Int_t nEvents = 20;
//        Int_t nEvents = 1000000;
        if(nEv > 100) nEvents = 1000000;

        TString out = outputFile;
        out+=".root";

        TFile *outfile = new TFile(out.Data(), "RECREATE");
        outfile->Close();

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

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        // =============================================================================== //
        // over-ride functions
        if(dopp) {
          doCentSelection = kFALSE;  // can't ask for a particular centrality if requesting pp collisions
          doBackgroundJets = kFALSE; // don't do a rho background calculation when using pp collisions   
        }

        // input file for tests (based on Run) - updated for new Runs as needed
        if((RunYear == mRun12) && doTEST) inputFile = "Run12.list";
        //if((RunYear == mRun14) && doTEST) inputFile = "Run14_P18ih_HPSS_15164046.list"; //Run_15151042_files.list"; //"testLIST_Run14.list";
        // if((RunYear == mRun14) && doTEST) inputFile = "/star/u/droy1/Y2019/STAR/Test_PythiaSmallSampleFileList.list";
        // if((RunYear == mRun14) && doTEST) inputFile = "/star/u/droy1/Y2019/STAR/PythiaSmallFileList_15_20.list";

        if((RunYear == mRun14) && doTEST) {
            inputFile = "Run14PicoDst.list";
        }
        // if((RunYear == mRun14) && doTEST) inputFile = "SampleFileListVertex.list";
        if((RunYear == mRun16) && doTEST) inputFile = "test_run17124003_files.list";
        if((RunYear == mRun17) && doTEST && dopp) inputFile = "filelist_pp2017.list";
        cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags - no centrality for pp collisions
        Int_t CentralitySelection;
        Int_t CentralityDefinition;
        if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // no centrality defintion for pp, just set one 
        if(RunYear == mRun14) CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P16id; // (NEW - from Nick Aug22, 2019: all lumi) this will be default
        if(RunYear == mRun14 && MCSim) CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // no centrality defintion for pp, just set one 
        //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // Run14 P18ih (NEW - from Nick June10, 2019)
        if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id;         // Run16 - option: StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx;
        cout<<"Centrality definition: "<<CentralityDefinition<<endl;

        // Run/Event Flag
        Int_t RunFlag;
        if(RunYear == mRun12 && dopp) RunFlag = StJetFrameworkPicoBase::Run12_pp200;
        if(RunYear == mRun14) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
        if(RunYear == mRun16) RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
        if(RunYear == mRun17 && dopp) RunFlag = StJetFrameworkPicoBase::Run17_pp510;
        Bool_t RejectBadRuns = kFALSE; // switch to load and than omit bad runs
        Int_t fBadRunListVers = StJetFrameworkPicoBase::fBadRuns_w_missing_HT; // fBadRuns_w_missing_HT, fBadRuns_wo_missing_HT,

        // trigger flags - update default
        Int_t EmcTriggerEventType; // kIsHT1 or kIsHT2 or kIsHT3
        if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;
        if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;
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

        // open and close output .root file (so it exist and can be updated by Analysis Tasks)
        // TFile *fout = new TFile(outputFile, "RECREATE");
        // fout->Close();

        // create the analysis maker!
        bool doComments = kFALSE;

        // create chain
        StChain* chain = new StChain();

  // create the picoMaker maker
  //StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
        StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst");
        picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::VpdOrDefault));
        // picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));


        // create base class maker pointer
        StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
        baseMaker->SetRunFlag(RunFlag);                  // run flag (year)
        baseMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        baseMaker->SetBadRunListVers(fBadRunListVers);          // switch to select specific bad run version file
        baseMaker->SetBadTowerListVers(TowerListToUse);
        baseMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
        cout<<baseMaker->GetName()<<endl;  // print name of class instance

        // create centrality class maker pointer
        StCentMaker *CentMaker = new StCentMaker("CentMaker", picoMaker, out.Data(), doComments);
        CentMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
        CentMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
        CentMaker->SetRunFlag(RunFlag);                         // Run Flag
        CentMaker->SetdoppAnalysis(dopp);                       // pp-analysis switch
        CentMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
        CentMaker->SetUseBBCCoincidenceRate(kFALSE);            // BBC or ZDC (default) rate used?
        CentMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        CentMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        cout<<CentMaker->GetName()<<endl;  // print name of class instance

        StTestClass *Test = new StTestClass("Test", picoMaker, out.Data());
        cout << Test->GetName() << endl;

        // initialize chain
        chain->Init();
        cout<<"chain->Init();"<<endl;
        int total = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << total << endl;
       // if(nEvents > total) nEvents = total;

        nEvents = TMath::Min(total, nEvents);
  
        for (Int_t i = 0; i < nEvents; i++){
          if (doTEST){
            if (i%1 == 0){
              cout << "****************************************** " << endl;
              cout << "Working on eventNumber " << i << endl;
            }
          }
          else{
            if(i%10 == 0) {
                // cout << "****************************************** " << endl;
                cout << "Working on eventNumber " << i << endl;
            }
          }
          
            // cout << "Working on eventNumber " << i << endl;

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
        if(outfile->IsOpen())   outfile->Close();

        //StMemStat::PrintMem("load StChain");
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

void LoadMacros()
{
}
