#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "MuHijing.h"
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
  TTree* mutree = new TTree("MuTree", "MuTree");
  mutree->SetDirectory(rootFile);

  // Create variables to hold the event and particle information
  Int_t eventNumber;
  MuHijing mu;

  // Create branches in the TTree
  mutree->Branch("NParticles", &mu.nparticles, "nparticles/I");
  mutree->Branch("TrackPt", mu.mPt, "mPt[nparticles]/F");
  mutree->Branch("TrackEta", mu.mEta, "mEta[nparticles]/F");
  mutree->Branch("TrackPhi", mu.mPhi, "mPhi[nparticles]/F");
  mutree->Branch("TrackEnergy", mu.mEnergy, "mEnergy[nparticles]/F");
  mutree->Branch("TrackPID", mu.mPID, "mPID[nparticles]/I");

  // Variables for reading lines from the file
  std::string line;
  std::istringstream iss;

  int eventnum = 0;
  double tmp;
  int status;

  int id;
  double px,py,pz,energy;

  int goodparticles = 0;

  // Read the file line by line
  while (std::getline(file, line)) {
    iss.clear();
    iss.str(line);

    // Check the line type
    char lineType;
    iss >> lineType;

    if (lineType == 'E') {
      // Fill the TTree
      if (eventnum!=0) {
        // mu.nparticles = goodparticles;

        // cout << eventnum << "\t" << mu.nparticles << endl;
        // for (int i = 0; i < mu.nparticles; i++){
        //   cout << Form("Particle ID = %i\t%.2f\t%.2f\t%.2f", mu.mPID[i], mu.mPt[i], mu.mEta[i], mu.mPhi[i]) << endl;
        // }
        
        mutree->Fill();
      }

      // Event line
      iss.ignore(1);  // Ignore the event entry number
      iss >> eventNumber;

      eventnum++;

      mu.nparticles = 0;
      // goodparticles = 0;
      
      for (int i = 0; i < ConstMax; i++){
        mu.mPt[i] = 0;
        mu.mEta[i] = 0;
        mu.mPhi[i] = 0;
        mu.mEnergy[i] = 0;
        mu.mPID[i] = 0;
      }
      
      // if (eventnum>10) break;
      // if (eventnum%1==0) cout << "Processing event " << eventnum << "\r" << flush;
      if (eventnum%1000==0) cout << "Processing event " << eventnum << endl;
      
    } 
    else if (lineType == 'P') {
      // Particle line
      // Particle particle;
      iss.ignore(1);  // Ignore the particle entry number
      iss >> id;
      iss >> id;
      iss >> px;
      iss >> py;
      iss >> pz;
      iss >> energy;
      iss >> tmp; // mass (Not saving it because pdg mass is available)
      iss >> status;

      TVector3 p(px, py, pz);

      if (status != 1) continue; //Only final state particles saved out

      mu.nparticles++;

      mu.mPt[mu.nparticles-1] = p.Pt();
      mu.mEta[mu.nparticles-1] = p.Eta();
      mu.mPhi[mu.nparticles-1] = p.Phi();
      mu.mEnergy[mu.nparticles-1] = energy;
      mu.mPID[mu.nparticles-1] = id;

      // cout << Form("Particle ID = %i\t%.2f\t%.2f\t%.2f\t%i", mu.mPID[goodparticles-1], mu.mPt[goodparticles-1], mu.mEta[goodparticles-1], mu.mPhi[goodparticles-1], status) << endl;
      // cout << Form("Particle ID = %i\t%.2f\t%.2f\t%i", particle.id, particle.px, particle.energy, particle.status) << endl;
    }
  }

  cout << "\n";

  // Write the TTree to the ROOT file
  mutree->Write();

  // Clean up
  rootFile->Close();
  delete rootFile;

  std::cout << "Conversion completed successfully!" << std::endl;
}

void JEWELToRoot(TString input = "./smallhepmc.hepmc") {
  const char* inputFile = input.Data();
  const char* outputFile = "out.root";

  ConvertHepMC2ToROOT(inputFile, outputFile);
}
