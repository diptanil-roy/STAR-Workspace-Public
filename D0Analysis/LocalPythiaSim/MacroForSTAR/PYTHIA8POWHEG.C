#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "Pythia8Plugins/PowhegHooks.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include <TLorentzVector.h>
#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
// #include "Pythia8/Pythia.h"
using namespace std;
#endif
using namespace Pythia8;

// Include Fastjet Classes
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/Selector.hh>

using namespace Pythia8;

// Defining a User Class to store PID info for particles used as input to FastJet
class MyParticleId : public fastjet::PseudoJet::UserInfoBase{
    public:
        MyParticleId(const int & pdg_id_in, const int & part_id_in) : _pdg_id(pdg_id_in), _part_id(part_id_in){}

        int pid() const{
            return _pdg_id;
        }
        int partn() const{
            return _part_id;
        }
    private:
        int _pdg_id;
        int _part_id;
};


// function to calculate invariant mass of a pair of objects
//___________________________________________________________________________________________
Double_t Mass(Double_t m1, Double_t m2, Double_t E1, Double_t E2, Double_t px1, Double_t px2, Double_t py1, Double_t py2, Double_t pz1, Double_t pz2) {
  Double_t m;
  m = TMath::Sqrt(pow(m1,2) + pow(m2,2) + 2*(E1*E2 - px1*px2 - py1*py2 - pz1*pz2));

  return m;
}

Double_t standardPhi(Double_t phi){
  Double_t phi_standard = phi;
  if (phi_standard < 0) phi_standard+=2*(TMath::Pi()); //FIXME
  if (phi_standard < 0) cout << "Something wrong with angle!" << endl;
  return phi_standard;
}

// This function is to count the number of events in the vertex file.
int countlines(const char* VertexFileName = "./VERTEXFULL/st_physics_15094070_raw_0000007.txt"){
    int numberoflines = 0;
    std::string line;
    std::ifstream file(VertexFileName, ios::in);

    while (std::getline(file, line))
        ++numberoflines;

    return numberoflines;
}

void method(Int_t nev = 100, const int random_seed = 18134, const char* VertexFileName = "./VERTEXFULL/st_physics_15094070_raw_0000007.txt", float R = 0.4, Int_t power = -1, Bool_t printjetinfo = kFALSE, Bool_t smallOutput = kTRUE){

  //nev = Number of Events (Equal to the number of events asked for, or the number of events in the vertex file, whichever is less.)
  //Random Seed = set by hand to be able to regenerate the dataset if needed.
  //Vertexfilename = Text files containing the vertex information from picodsts
  //R = Radius of Jets
  //Power = Defines the jet reconstruction algorithm, -1 stands for anti-kT and doesn't need to be changed
  //printjetinfo = For debugging



  TStopwatch timer;

  timer.Start();

  gRandom = new TRandom3();

  // pT HARD Bins defined here
  double pTlimit[2] = {3., 0.}; 

  // FileName for the .root Pythia Events Files, salvaged from the vertex files names to ensure correspondence
  TString filename;

  filename = VertexFileName;
  filename.ReplaceAll(".txt", "");
  filename.ReplaceAll("./VERTEXFULL/Vertex_", "");
  filename += "_" + std::to_string(int(pTlimit[0])) + "_" + std::to_string(int(pTlimit[1]));
  filename += "_POWHEG";

  // FileName for the .txt Pythia Events Files, salvaged from the vertex files names to ensure correspondence
  TString textfilename;
  textfilename = filename + ".txt";
  filename += ".root";

  // Defining a file to save Jets
  TFile *jetfile = new TFile(filename.Data(), "RECREATE");

  cout << textfilename << endl;

  ofstream txfile (textfilename.Data(), ios::out | ios::binary );

  // Create pythia8 object
  TPythia8        *pythia = new TPythia8(); 
  Pythia8::Pythia *pythia8  = pythia->Pythia8();
  TDatabasePDG *pdg = new TDatabasePDG();

  // pythia8->readString("HardQCD:gg2ccbar = on");
  // pythia8->readString("HardQCD:qqbar2ccbar = on");
  // pythia8->readString("Charmonium:all = on");

  pythia8->readFile("main.cmnd");

  // Read in main settings
  int nEvent      = pythia8->settings.mode("Main:numberOfEvents");
  int nError      = pythia8->settings.mode("Main:timesAllowErrors");
  // Read in key POWHEG merging settings
  int vetoMode    = pythia8->settings.mode("POWHEG:veto");
  int MPIvetoMode = pythia8->settings.mode("POWHEG:MPIveto");
  bool loadHooks  = (vetoMode > 0 || MPIvetoMode > 0);

  cout << nEvent << endl;

  // Add in user hooks for shower vetoing
  shared_ptr<PowhegHooks> powhegHooks;
  if (loadHooks) {

    // Set ISR and FSR to start at the kinematical limit
    if (vetoMode > 0) {
      pythia8->readString("TimeShower:pTmaxMatch = 2");
    }

    // Set MPI to start at the kinematical limit
    if (MPIvetoMode > 0) {
      pythia8->readString("MultipartonInteractions:pTmaxMatch = 2");
    }

    powhegHooks = make_shared<PowhegHooks>();
    pythia8->setUserHooksPtr((UserHooksPtr)powhegHooks);
  }

  // Configure Pythia to generate events
  // I will use all HardQCD processes, and make sure the D0s do not decay. Perform a check by invoking isFinal for PID 421

  pythia->ReadString("Random:setSeed = on");
  pythia->ReadString(Form("Random:seed = %i", random_seed));

  if (smallOutput) {
      pythia8->readString("Init:showProcesses = off");
      pythia8->readString("Init:showMultipartonInteractions = off");
      pythia8->readString("Init:showChangedSettings = off");
      pythia8->readString("Init:showChangedParticleData = off");
      pythia8->readString("Next:numberCount = 1000000000");
      pythia8->readString("Next:numberShowInfo = 0");
      pythia8->readString("Next:numberShowProcess = 0");
      pythia8->readString("Next:numberShowEvent = 0");
  }

  


  /*
  D0 Decay Stopped, Will happen in GEANT
  */

  // pythia8->particleData.readString("421:mayDecay = 0");

  // Turning off weak and long-lived particle decays! They will only be decayed by GEANT

  /*
  111:mayDecay = 0
  211:mayDecay = 0
  221:mayDecay = 0
  321:mayDecay = 0
  130:mayDecay = 0
  310:mayDecay = 0
  3122:mayDecay = 0
  3212:mayDecay = 0
  3112:mayDecay = 0
  3222:mayDecay = 0
  3312:mayDecay = 0
  3322:mayDecay = 0
  3334:mayDecay = 0
  */

  pythia8->particleData.readString("111:mayDecay = 0");
  pythia8->particleData.readString("211:mayDecay = 0");
  pythia8->particleData.readString("221:mayDecay = 0");
  pythia8->particleData.readString("321:mayDecay = 0");
  pythia8->particleData.readString("130:mayDecay = 0");
  pythia8->particleData.readString("310:mayDecay = 0");
  pythia8->particleData.readString("3122:mayDecay = 0");
  pythia8->particleData.readString("3212:mayDecay = 0");
  pythia8->particleData.readString("3112:mayDecay = 0");
  pythia8->particleData.readString("3222:mayDecay = 0");
  pythia8->particleData.readString("3312:mayDecay = 0");
  pythia8->particleData.readString("3322:mayDecay = 0");
  pythia8->particleData.readString("3334:mayDecay = 0");

  pythia8->particleData.readString("421:mayDecay = 0");

  // pythia8->particleData.listChanged();

  // Open the vertex file
  fstream vertexfile;
  vertexfile.open(VertexFileName, ios::in);

  // The same macro can be looped over multiple pt hard bins but for now, it is simpler to change the pt hard bins by hand, and have different macros for each bin
  for (int ipT = 0; ipT < 1; ipT++){

    cout << "pT Bin : " << pTlimit[ipT] << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << pTlimit[ipT + 1] << endl; 
    // pythia8->settings.parm("PhaseSpace:pTHatMin", pTlimit[ipT]);
    // pythia8->settings.parm("PhaseSpace:pTHatMax", pTlimit[ipT + 1]);

    TH1F *hJetPt = new TH1F("hJetPt", "hJetPt", 2000, 0.0, 200.0);
    TH1F *hJetEta = new TH1F("hJetEta", "hJetEta", 2000, -10.0, 10.0);
    TH1F *hJetPhi = new TH1F("hJetPhi", "hJetPhi", 2000, -10.0, 10.0);

    TH1F *hTrackPt = new TH1F("hTrackPt", "hTrackPt", 2000, 0.0, 200.0);
    TH1F *hTrackEta = new TH1F("hTrackEta", "hTrackEta", 2000, -10.0, 10.0);
    TH1F *hTrackPhi = new TH1F("hTrackPhi", "hTrackPhi", 2000, -10.0, 10.0);

    TH1F *hD0Pt = new TH1F("hD0Pt", "hD0Pt", 2000, 0.0, 200.0);
    TH1F *hD0Eta = new TH1F("hD0Eta", "hD0Eta", 2000, -10.0, 10.0);
    TH1F *hD0Phi = new TH1F("hD0Phi", "hD0Phi", 2000, -10.0, 10.0);

    TH1F *histZ = new TH1F("histZ", "histZ", 150, 0, 1.5);

    // =========================== Initialising FastJet =====================================

    // Will use E_scheme recombination for this analysis. To choose the others, replace any of these in the fastjet::RecombinationScheme: 
    // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme
    fastjet::Strategy strategy = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; //Change as you need
    fastjet::JetDefinition *jetdefinition = NULL;
    fastjet::JetAlgorithm algorithm;

    if (power == -1)      algorithm = fastjet::antikt_algorithm; // Will use this one for now
    if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
    if (power ==  1)      algorithm = fastjet::kt_algorithm;

    // Defining the Jet
    jetdefinition = new fastjet::JetDefinition(algorithm, R, recombScheme, strategy);

    // Defining Vectors for FastJet inputs
    std::vector <fastjet::PseudoJet> fjInputs; // Will store px, py, pz, E info



    // Initialise Pythia
    // Beam 1 = p
    // Beam 2 = p
    // Ecm =  200 GeV/c^2
    // pythia->Initialize(2212, 2212, 200.);
    // pythia->Initialize();
    pythia8->init();

    cout << "Particles are : " << pythia8->info.idA() << "\t" << pythia8->info.idB() << endl;
    cout << "Energies are : " << pythia8->info.eA() << "\t" << pythia8->info.eB() << endl;
    
    cout << "Pythia initialised." << endl;

    // TClonesArray to access the list of particles generated by Pythia
    TClonesArray *particles = (TClonesArray*)pythia->GetListOfParticles();

    int d0Events = 0; // To count the number of D0 events
    int totalD0s = 0;
    int numberofd0Jetsabove3GeV = 0;

    int requiredgoodd0events = nev;

    float pthat;
    float eventCrossSection;

    cout << "Asked for " << requiredgoodd0events << " events" << endl;

    // EVENT LOOP
    for (int event = 0; event < 10000000000; event++){
      // if (d0Events % 10 == 0) {cout << "D0 Event " << d0Events << endl;}

      if (numberofd0Jetsabove3GeV >= requiredgoodd0events || pythia8->info.atEndOfFile()) break;

      if (event % 1000 == 0) { cout << "Event # " << event << "\t" << d0Events << "\t" << numberofd0Jetsabove3GeV << endl;}
      pythia->GenerateEvent();
      //if (event < 1) pythia->EventListing();

      //Import all particles generated to the TClonesArray defined earlier)
      pythia->ImportParticles(particles, "All");
      Int_t npart = particles->GetEntriesFast();

      //Reset FastJet input for every event
      fjInputs.clear();
      fjInputs.resize(0);


      double weight = pythia8->info.weight();
      // pythia8->info.list();
      pthat = pythia8->info.pTHat();
      eventCrossSection = pythia8->info.sigmaGen();

      if (printjetinfo) cout << "Weight and pTHat is = " << weight << "\t" << eventCrossSection << "\t" << pthat << endl;

      std::vector<std::vector<int> >  fd0TrackIndices;
      fd0TrackIndices.clear(); // Save the D0 track indices 

      bool d0Event = kFALSE; // Boolean to select a D0 event

      for (int part1 = 0; part1 < npart; part1++){
        TParticle *MPart1 = (TParticle *) particles->At(part1);
        int pid1 = MPart1->GetPdgCode();
        Int_t ist = MPart1->GetStatusCode();

        if (MPart1->Pt() < 0.2 || abs(MPart1->Eta()) > 1) continue; // No point in having D0s beyond our acceptance!

        if (ist > 0) {hTrackPt->Fill(MPart1->Pt()); hTrackPhi->Fill(MPart1->Phi()); hTrackEta->Fill(MPart1->Eta());}

        if (abs(pid1) == 421 && MPart1->Pt() > 1.0) {d0Event = kTRUE; break;}
      }

      if (!d0Event) continue;
      d0Events++;

      for (int part = 0; part < npart; part++){
        TParticle *MPart = (TParticle *) particles->At(part);
        int pid = MPart->GetPdgCode();
        Int_t ist = MPart->GetStatusCode();

        if (ist <= 0) continue;

        if (abs(pid) == 421) {totalD0s++; hD0Pt->Fill(MPart->Pt()); hD0Phi->Fill(MPart->Phi()); hD0Eta->Fill(MPart->Eta());}

        if (MPart->Pt() < 0.2 || abs(MPart->Eta()) > 1) continue; // No point in having D0s beyond our acceptance!

        fastjet::PseudoJet p(MPart->Px(), MPart->Py(), MPart->Pz(), MPart->Energy());
        p.set_user_info(new MyParticleId(pid, part));
        // Store as input to Fastjet
        fjInputs.push_back(p);
      }

      vector <fastjet::PseudoJet> inclusiveJets, sortedJets, selectedJets;
      fastjet::ClusterSequence clustSeq(fjInputs, *jetdefinition);
      
      // Extract Inclusive Jets sorted by pT
      inclusiveJets = clustSeq.inclusive_jets(3.); //Only jets with pT > 3 GeV/c are selected
      sortedJets = sorted_by_pt(inclusiveJets);

      fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-0.6, 0.6); // Selects jets with |eta| < 0.6
      fastjet::Selector selector = eta_selector;
      selectedJets = eta_selector(sortedJets);

      // cout << event << "\t" << totalD0s << "\t" << numberofd0Jetsabove3GeV << "\t" << selectedJets.size() << endl;


      for(unsigned jetid = 0; jetid < selectedJets.size(); jetid++){
        // cout<<"jet " <<jetid<<": "<<selectedJets[jetid].pt()<<" "<<selectedJets[jetid].eta()<<" "<<selectedJets[jetid].phi()<<endl;
        // Loop over jet constituents to get some info
        vector <fastjet::PseudoJet> constituents = selectedJets[jetid].constituents();
        vector <fastjet::PseudoJet> sortedconstituents = sorted_by_pt(constituents);

        bool d0Jet = kFALSE;

        for (unsigned j = 0; j < sortedconstituents.size(); j++){
          if (abs(sortedconstituents[j].user_info<MyParticleId>().pid()) == 421 && sortedconstituents[j].pt() > 1.){
              d0Jet = kTRUE; 
          }
          if (d0Jet) { break;}
        }

        if (d0Jet) {
          numberofd0Jetsabove3GeV++;

          hJetPt->Fill(selectedJets[jetid].pt());
          hJetEta->Fill(selectedJets[jetid].eta());
          hJetPhi->Fill(standardPhi(selectedJets[jetid].phi()));

          int d0ID = -99;

          for (unsigned j = 0; j < sortedconstituents.size(); j++){
            if (abs(sortedconstituents[j].user_info<MyParticleId>().pid()) == 421){ d0ID = j; break;}
          }

          Double_t z = (sortedconstituents[d0ID].px()*selectedJets[jetid].px() + sortedconstituents[d0ID].py()*selectedJets[jetid].py())/(pow(selectedJets[jetid].pt(), 2));
          histZ->Fill(z);
        }
      }
    } // END OF EVENT LOOP

    cout << "Number of D0 Events = " << d0Events << endl;

    Double_t weightsum = pythia8->info.weightSum();
    Double_t sigmagen = pythia8->info.sigmaGen();

    jetfile->cd();

    hJetPt->Write();
    hJetEta->Write();
    hJetPhi->Write();


    hD0Pt->Write();
    hD0Eta->Write();
    hD0Phi->Write();

    hTrackPt->Write();
    hTrackEta->Write();
    hTrackPhi->Write();

    histZ->Write();
  }

  jetfile->Close();
  txfile.close();

  timer.Stop();
  cout << "Real Time Used: " << timer.RealTime()/60 << "m" << endl;

}

void PYTHIA8POWHEG(Int_t nev = 20000, const int random_seed = 18134, const char* VertexFileName = "./VERTEXFULL/Vertex_st_physics_15094070_raw_0000007.txt"){
  method(nev, random_seed, VertexFileName);
}

