R__ADD_INCLUDE_PATH("/opt/homebrew/Cellar/hepmc2/2.06.11/include/");
R__LOAD_LIBRARY(/opt/homebrew/Cellar/hepmc2/2.06.11/lib/libHepMC.dylib);

// #include "HepMC/IO_GenEvent.h"
#include "/opt/homebrew/Cellar/hepmc2/2.06.11/include/HepMC/GenEvent.h"

#include <iostream>

using namespace HepMC;
using namespace std;

void JewelReaderHepMC2(){
    
    TString InputFile = "/Volumes/WorkDrive/JEWEL/0010_AuAu200GeV_dijet_wRec_1.hepmc";

    TString OutputFile = "/Volumes/WorkDrive/JEWEL/0010_AuAu200GeV_dijet_wRec_1.root";

    TFile *f = new TFile(OutputFile, "RECREATE");

    MuHijing mu;

    TTree *mutree = new TTree("MuTree", "MuTree");
    mutree->SetDirectory(f);

    mutree->Branch("NParticles", &mu.nparticles, "nparticles/I");
    mutree->Branch("TrackPt", &mu.mPt, "mPt[nparticles]/F");
    mutree->Branch("TrackEta", &mu.mEta, "mEta[nparticles]/F");
    mutree->Branch("TrackPhi", &mu.mPhi, "mPhi[nparticles]/F");
    mutree->Branch("TrackEnergy", &mu.mEnergy, "mEnergy[nparticles]/F");
    mutree->Branch("TrackPID", &mu.mPID, "mPID[nparticles]/I");

    // ReaderAsciiHepMC2 text_input (InputFile.Data());

    IO_GenEvent inputFile(InputFile.Data(), ios::in);

    GenEvent* event = new HepMC3::GenEvent();

    int events_parsed = 0;

    while (!inputFile.failed()) {
        inputFile >> event;
        if (inputFile.failed()) break;
    }

    // while( !text_input.failed() ) {

    //     cout << events_parsed << endl;
    //     // if (events_parsed > 2) break;

    //     GenEvent evt(Units::GEV,Units::MM);

    //     text_input.read_event(evt);

    //     if( text_input.failed() ) break;

    //     if( events_parsed == 0 ) {
    //         cout << "First event: " << endl;
    //         Print::listing(evt);
    //     }
        
    //     ++events_parsed;
    // }
}
