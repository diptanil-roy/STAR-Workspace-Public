#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
// #include "MuHijing.h"
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
  // MuHijing mu;

  int nparticles = 0;
  float mPt[5000];
  float mEta[5000];
  float mPhi[5000];
  float mEnergy[5000];
  int mPID[5000];

  // Create branches in the TTree
  mutree->Branch("NParticles", &nparticles, "NParticles/I");
  mutree->Branch("TrackPt", mPt, "TrackPt[NParticles]/F");
  mutree->Branch("TrackEta", mEta, "TrackEta[NParticles]/F");
  mutree->Branch("TrackPhi", mPhi, "TrackPhi[NParticles]/F");
  mutree->Branch("TrackEnergy", mEnergy, "TrackEnergy[NParticles]/F");
  mutree->Branch("TrackPID", mPID, "TrackPID[NParticles]/I");

  // Variables for reading lines from the file
  std::string line;
  std::istringstream iss;

  int eventnum = 0;
  double tmp;
  int status;

  int id;
  double px,py,pz,energy;

  // Read the file line by line
  while (std::getline(file, line)) {
    iss.clear();
    iss.str(line);

    // Check the line type
    char lineType;
    iss >> lineType;

    if (lineType == 'E') {
      // Fill the TTree
      if (eventnum != 0) {
        // cout << eventnum << "\t" << nparticles << endl;
        // for (int i = 0; i < nparticles; i++){
        //   cout << Form("Particle ID = %i\t%.2f\t%.2f\t%.2f", mPID[i], mPt[i], mEta[i], mPhi[i]) << endl;
        // }
        mutree->Fill();
      }

      // Event line
      iss.ignore(1);  // Ignore the event entry number
      iss >> eventNumber;

      eventnum++;

      nparticles = 0;
      // goodparticles = 0;
      
      for (int i = 0; i < 5000; i++){
        mPt[i] = 0;
        mEta[i] = 0;
        mPhi[i] = 0;
        mEnergy[i] = 0;
        mPID[i] = 0;
      }
      
      // if (eventnum>20000) break;
      if (eventnum%1==0) cout << "Processing event " << eventnum << "\r" << flush;
      // if (eventnum%10000==0) cout << "Processing event " << eventnum << endl;
      
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

      mPt[nparticles] = p.Pt();
      mEta[nparticles] = p.Eta();
      mPhi[nparticles] = p.Phi();
      mEnergy[nparticles] = energy;
      mPID[nparticles] = id;

      nparticles++;

      // cout << Form("Particle ID = %i\t%.2f\t%.2f\t%.2f\t%i", mPID[nparticles-1], mPt[nparticles-1], mEta[nparticles-1], mPhi[nparticles-1], status) << endl;
      // cout << Form("Particle ID = %i\t%.2f\t%.2f\t%i", particle.id, particle.px, particle.energy, particle.status) << endl;
    }
  }

  cout << "\n";

  // Write the TTree to the ROOT file
  // mutree->Print();
  mutree->Write();

  // Clean up
  rootFile->Close();
  delete rootFile;

  std::cout << "Conversion completed successfully!" << std::endl;
}

void JEWELToRoot(TString input = "smallhepmc.hepmc") {
  const char* inputFile = input.Data();
  TString output = input;
  output.ReplaceAll(".hepmc", ".root");
  const char* outputFile = output.Data();

  cout << "Input file: " << inputFile << endl;
  cout << "Output file: " << outputFile << endl;

  ConvertHepMC2ToROOT(inputFile, outputFile);
}
