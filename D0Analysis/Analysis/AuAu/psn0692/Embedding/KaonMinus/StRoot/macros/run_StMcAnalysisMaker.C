#include <iostream>
#include <string>

void Load();

using namespace std;

void run_StMcAnalysisMaker(const char* file, std::string outFile = "test")
{
   TStopwatch*   stopWatch = new TStopwatch();
   stopWatch->Start();

   //Check STAR Library. Please set SL_version to the original star library used
   // in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl

   string SL_version = "SL16d";
   string env_SL = getenv("STAR");

   if (env_SL.find(SL_version) == string::npos)
   {
      cout << "Environment Star Library does not match the requested library in run_st_etree.C. Exiting..." << endl;
      return;
   }

   // Load shared libraries
   Load();

   // Create chain
   StChain* chain = new StChain;

   // I/O maker
   StIOMaker* ioMaker = new StIOMaker;
   ioMaker->SetFile(file);
   ioMaker->SetIOMode("r");
   ioMaker->SetBranch("*", 0, "0");
   //ioMaker->SetBranch("McEventBranch",0,"r");
   ioMaker->SetBranch("geantBranch", 0, "r");
   ioMaker->SetBranch("eventBranch", 0, "r");

   TString mudstfile = file;

   if(mudstfile.First("$") != -1)
   {
     mudstfile.ReplaceAll("$","");
     mudstfile = getenv(mudstfile.Data());
   }

   mudstfile.ReplaceAll(".event.root", ".MuDst.root");
   mudstfile.ReplaceAll(".geant.root", ".MuDst.root");
   cout << "Reading MuDst file " << mudstfile << endl;
   StMuDstMaker* muDstMaker = new StMuDstMaker(0, 0, "", mudstfile.Data(), "", 100000, "MuDst");

   StMcEventMaker *mcEventMaker = new StMcEventMaker();
   mcEventMaker->doPrintEventInfo = false;
   mcEventMaker->doPrintMemoryInfo = false;

   StAssociationMaker* assoc = new StAssociationMaker;
   assoc->useInTracker();
   assoc->SetDebug();

   //.. see example in CVS: StRoot/macros/mudst/exampleEmc.C
   // Need St_db_Maker for Emc calibration
   // St_db_Maker* dbMk = new St_db_Maker("StarDb", "MySQL:StarDb");
   //dbMk->SetMaxEntryTime(20100301,0);
   //dbMk->SetDateTime(20080101,000001);
   int refMultCorrLoad = gSystem->Load("StRefMultCorr");
   StRefMultCorr* grefmultCorrUtil = NULL;

   if (refMultCorrLoad == -1)
   {
      cout << "StRefMultCorr library is not available" << endl;
   }
   else
   {
      grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id();
      grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
      grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");

      // for (int i = 0; i < 6; ++i)
      // {
      //    cout << i << " " << grefmultCorrUtil->get(i, 0) << endl;
      // }
   }

   // Monte Carlo event maker
   StMcAnalysisMaker* analysis = new StMcAnalysisMaker;
   analysis->setOutFileName(outFile);
   // analysis->fillTpcHitsNtuple();
   analysis->setRefMultCorr(grefmultCorrUtil);

   // Initialize chain
   chain->Init();
   chain->EventLoop(1e6);
   chain->Finish();

   //delete chain;
   stopWatch->Stop();   
   stopWatch->Print();
}

void Load()
{
   gROOT->Macro("loadMuDst.C");
   gROOT->Macro("LoadLogger.C");
   gSystem->Load("StMcEvent");
   gSystem->Load("StMcEventMaker");
   gSystem->Load("StAssociationMaker");
   gSystem->Load("StDbLib.so");
   gSystem->Load("StDbBroker.so");
   gSystem->Load("libglobal_Tables.so");
   gSystem->Load("St_db_Maker.so");
   gSystem->Load("StDetectorDbMaker");
   gSystem->Load("StTpcDb");
   gSystem->Load("StDbUtilities");
   gSystem->Load("StMcEvent");
   gSystem->Load("StMcEventMaker");
   gSystem->Load("StAssociationMaker");
   gSystem->Load("StEmcRawMaker");
   gSystem->Load("StEmcADCtoEMaker");
   gSystem->Load("StPreEclMaker");
   gSystem->Load("StEpcMaker");
   gSystem->Load("StEmcSimulatorMaker");
   gSystem->Load("StEmcUtil");
   gSystem->Load("StEEmcUtil");
   gSystem->Load("StEEmcDbMaker");
   gSystem->Load("StEmcTriggerMaker");
   gSystem->Load("StDaqLib");
   gSystem->Load("StMcAnalysisMaker");
}
