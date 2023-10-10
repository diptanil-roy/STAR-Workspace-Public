R__LOAD_LIBRARY(libEG);
R__LOAD_LIBRARY(libEGPythia8);

R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/FASTJET/fastjet-install/lib/libfastjet.dylib);


#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
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

struct JetSaver{
  float pt;
  float eta;
  float phi;
  float D0pt;
  float D0eta;
  float D0phi;
  float z;
};

void method(Int_t nev = 100, const int random_seed = 18134, const char* FileName = "./VERTEXFULL/st_physics_15094070_raw_0000007.txt", float R = 0.4, Int_t power = -1, Bool_t printjetinfo = kFALSE, Bool_t smallOutput = kTRUE){

  //nev = Number of Events (Equal to the number of events asked for, or the number of events in the vertex file, whichever is less.)
  //Random Seed = set by hand to be able to regenerate the dataset if needed.
  //Vertexfilename = Text files containing the vertex information from picodsts
  //R = Radius of Jets
  //Power = Defines the jet reconstruction algorithm, -1 stands for anti-kT and doesn't need to be changed
  //printjetinfo = For debugging



  TStopwatch timer;

  timer.Start();

  // Load Pythia 8
  // gSystem->Load("libEG");
  // gSystem->Load("libEGPythia8");

  gRandom = new TRandom3();

  // pT HARD Bins defined here
  double pTlimit[2] = {3, -1}; 

  // FileName for the .root Pythia Events Files, salvaged from the vertex files names to ensure correspondence

  // Defining a file to save Jets
  TFile *jetfile = new TFile(FileName, "RECREATE");

  cout << FileName << endl;

  // Create pythia8 object
  TPythia8        *pythia = new TPythia8(); 
  Pythia8::Pythia *pythia8  = pythia->Pythia8();
  TDatabasePDG *pdg = new TDatabasePDG();

  // Configure Pythia to generate events
  // I will use all HardQCD processes, and make sure the D0s do not decay. Perform a check by invoking isFinal for PID 421

  pythia->ReadString("Random:setSeed = on");
  pythia->ReadString(Form("Random:seed = %i", random_seed));

  pythia->ReadConfigFile("detroit.cmnd");

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

  // pythia8->readString("HardQCD:all = on");

  // pythia8->readString("HardQCD:gg2ccbar = on");
  // pythia8->readString("HardQCD:qqbar2ccbar = on");
  // pythia8->readString("Charmonium:all = on");

  pythia8->readString("HardQCD:gg2bbbar = on");
  pythia8->readString("HardQCD:qqbar2bbbar = on");
  pythia8->readString("Bottomonium:all = on");

  /*
  D0 Decay Stopped, Will happen in GEANT
  */

  pythia8->particleData.readString("421:mayDecay = 0");

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

  // pythia8->particleData.listChanged();

  // The same macro can be looped over multiple pt hard bins but for now, it is simpler to change the pt hard bins by hand, and have different macros for each bin
  for (int ipT = 0; ipT < 1; ipT++){

    cout << "pT Bin : " << pTlimit[ipT] << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << pTlimit[ipT + 1] << endl; 
    pythia8->settings.parm("PhaseSpace:pTHatMin", pTlimit[ipT]);
    pythia8->settings.parm("PhaseSpace:pTHatMax", pTlimit[ipT + 1]);

    JetSaver jet = {};

    TTree *jettree = new TTree("Jets", "Jets");
    jettree->SetDirectory(jetfile);

    jettree->Branch("JetPt", &jet.pt, "JetPt/F");
    jettree->Branch("JetEta", &jet.eta, "JetEta/F");
    jettree->Branch("JetPhi", &jet.phi, "JetPhi/F");
    jettree->Branch("D0Pt", &jet.D0pt, "D0Pt/F");
    jettree->Branch("D0Eta", &jet.D0eta, "D0Eta/F");
    jettree->Branch("D0Phi", &jet.D0phi, "D0Phi/F");
    jettree->Branch("D0Z", &jet.z, "D0Z/F");


    TH1F *hJetPt = new TH1F("hJetPt", "hJetPt", 100, 0.0, 50.0);
    TH1F *hJetEta = new TH1F("hJetEta", "hJetEta", 600, -1.2, 1.2);
    TH1F *hJetPhi = new TH1F("hJetPhi", "hJetPhi", 500, 0, 10.0);

    TH1F *hTrackPt = new TH1F("hTrackPt", "hTrackPt", 100, 0.0, 50.0);
    TH1F *hTrackEta = new TH1F("hTrackEta", "hTrackEta", 600, -1.2, 1.2);
    TH1F *hTrackPhi = new TH1F("hTrackPhi", "hTrackPhi", 500, 0, 10.0);

    TH1F *hD0Pt = new TH1F("hD0Pt", "hD0Pt", 400, 0.0, 20.0);
    TH1F *hD0Eta = new TH1F("hD0Eta", "hD0Eta", 600, -1.2, 1.2);
    TH1F *hD0Phi = new TH1F("hD0Phi", "hD0Phi", 500, 0, 10.0);
    TH1F *histZ = new TH1F("histZ", "histZ", 80, 0., 1.);
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
    pythia->Initialize(2212, 2212, 200.);
    
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
    for (int event = 0; event < nev; event++){

      // if (d0Events >= requiredgoodd0events) break;

      if (event % 1000 == 0) { cout << "Event # " << event << "\t" << d0Events << "\t" << numberofd0Jetsabove3GeV << "...\r" << flush;}
      if (event == nev-1) cout << endl;
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

        if (abs(pid1) == 421 && MPart1->Pt() > 1.0 && MPart1->Pt() < 10.0) {d0Event = kTRUE; break;}
      }

      if (!d0Event) continue;
      d0Events++;

      for (int part = 0; part < npart; part++){
        TParticle *MPart = (TParticle *) particles->At(part);
        int pid = MPart->GetPdgCode();
        Int_t ist = MPart->GetStatusCode();

        if (ist <= 0) continue;

        if (MPart->Pt() < 0.2) continue; // No point in having D0s beyond our acceptance!

        if (MPart->Pt() > 30.0) continue;

        if (abs(MPart->Eta()) > 1) continue;

        if (abs(pid) == 421 && MPart->Pt() > 1.0 && MPart->Pt() < 10.0) {totalD0s++; hD0Pt->Fill(MPart->Pt()); hD0Phi->Fill(MPart->Phi()); hD0Eta->Fill(MPart->Eta());}

        fastjet::PseudoJet p(MPart->Px(), MPart->Py(), MPart->Pz(), MPart->Energy());
        p.set_user_info(new MyParticleId(pid, part));
        // Store as input to Fastjet
        fjInputs.push_back(p);
      }

      vector <fastjet::PseudoJet> inclusiveJets, sortedJets, selectedJets;
      fastjet::ClusterSequence clustSeq(fjInputs, *jetdefinition);
      
      // Extract Inclusive Jets sorted by pT
      inclusiveJets = clustSeq.inclusive_jets(1.0); //Only jets with pT > 3 GeV/c are selected
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

        if (selectedJets[jetid].pt() > 20.) continue;

        for (unsigned j = 0; j < sortedconstituents.size(); j++){
          if (abs(sortedconstituents[j].user_info<MyParticleId>().pid()) == 421 && sortedconstituents[j].pt() > 1. && sortedconstituents[j].pt() < 20.){
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

          if (z == 1) z = 0.999;

          jet.pt = selectedJets[jetid].pt();
          jet.eta = selectedJets[jetid].eta();
          jet.phi = selectedJets[jetid].phi();
          jet.D0pt = sortedconstituents[d0ID].pt();
          jet.D0eta = sortedconstituents[d0ID].eta();
          jet.D0phi = sortedconstituents[d0ID].phi();
          jet.z = z;

          histZ->Fill(z);

          jettree->Fill();
        }
      }
    } // END OF EVENT LOOP

    cout << "Number of D0 Events = " << "\t" << d0Events << "\t"  << numberofd0Jetsabove3GeV << endl;

    Double_t weightsum = pythia8->info.weightSum();
    Double_t sigmagen = pythia8->info.sigmaGen();

    jetfile->cd();

    jettree->Write();

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

  timer.Stop();
  cout << "Real Time Used: " << timer.RealTime()/60 << "m" << endl;

}

void PYTHIA8ForD0Decay(Int_t nev = 2000000, const int random_seed = 2500001, const char* FileName = "PYTHIA_BBBar_FullRange.root"){

  method(nev, random_seed, FileName);
}

