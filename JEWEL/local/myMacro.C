#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "Particle.h"
#include "TVector3.h"

void ConvertHepMC2ToROOT(const char* inputFile, const char* outputFile) {
  // Open the HEPMC2 input file
  std::ifstream file(inputFile);
  if (!file.is_open()) {
    std::cerr << "Failed to open HEPMC2 input file." << std::endl;
    return;
  }

  // Create the ROOT output file
  TFile* rootFile = new TFile(outputFile, "RECREATE");
  TTree* tree = new TTree("Events", "Event Data");

  // Create variables to hold the event and particle information
  Int_t eventNumber;
  std::vector<Particle> particles;

  // Create branches in the TTree
  tree->Branch("EventNumber", &eventNumber, "EventNumber/I");
  tree->Branch("Particles", &particles);

  // Variables for reading lines from the file
  std::string line;
  std::istringstream iss;

  int eventnum = 0;
  double tmp;

  // Read the file line by line
  while (std::getline(file, line)) {
    iss.clear();
    iss.str(line);

    // Check the line type
    char lineType;
    iss >> lineType;

    if (lineType == 'E') {
      // Event line
      iss.ignore(1);  // Ignore the event entry number
      iss >> eventNumber;
      eventnum++;
      if (eventnum%10000==0) cout << "Processing event " << eventnum << "\r" << flush;
      // if (eventnum%1==0) cout << "Processing event " << eventnum << "\r" << endl;
      // if (eventnum>10000) break;
    } 
    else if (lineType == 'P') {
      // Particle line
      Particle particle;
      iss.ignore(1);  // Ignore the particle entry number
      iss >> particle.id;
      iss >> particle.id;
      iss >> particle.px;
      iss >> particle.py;
      iss >> particle.pz;
      iss >> particle.energy;
      iss >> tmp; // mass (Not saving it because pdg mass is available)
      iss >> particle.status;

      // TVector3 p(particle.px, particle.py, particle.pz);
      // if (p.Pt() > 2. && particle.status==1) cout << Form("Particle ID = %i\t%.2f\t%.2f\t%i", particle.id, p.Pt(), p.Eta(), particle.status) << endl;

      // std::cout << "Particle ID = " << particle.id << "\t" << particle.px << "\t" << particle.energy << "\t" << particle.status << endl;
      // cout << Form("Particle ID = %i\t%.2f\t%.2f\t%i", particle.id, particle.px, particle.energy, particle.status) << endl;
      particles.push_back(particle);
    }
  }

  cout << "\n";

  // Fill the TTree
  tree->Fill();

  // Write the TTree to the ROOT file
  rootFile->Write();

  // Clean up
  rootFile->Close();
  delete rootFile;

  std::cout << "Conversion completed successfully!" << std::endl;
}

void myMacro() {
//   const char* inputFile = "smallhepmc.hepmc";
  const char* inputFile = "/Volumes/WorkDrive/JEWEL/0010_AuAu200GeV_dijet_wRec_1.hepmc";
  const char* outputFile = "out.root";

  ConvertHepMC2ToROOT(inputFile, outputFile);
}
