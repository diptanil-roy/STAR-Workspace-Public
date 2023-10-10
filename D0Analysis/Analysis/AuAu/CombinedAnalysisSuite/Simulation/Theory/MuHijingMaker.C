#define Test_cxx
#include "Test.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "MuHijing.h"

using namespace std;

void Test::Loop()
{

   if (fChain == 0) return;

   MuHijing mu;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   TDatabasePDG *pdg = new TDatabasePDG();

   TFile *out = new TFile("PYTHIA.root", "RECREATE");
   TTree *mutree = new TTree("MuTree", "MuTree");
   mutree->SetDirectory(out);

   mutree->Branch("NParticles", &mu.nparticles, "nparticles/I");
   mutree->Branch("TrackPt", &mu.mPt, "mPt[nparticles]/F");
   mutree->Branch("TrackEta", &mu.mEta, "mEta[nparticles]/F");
   mutree->Branch("TrackPhi", &mu.mPhi, "mPhi[nparticles]/F");
   mutree->Branch("TrackEnergy", &mu.mEnergy, "mEnergy[nparticles]/F");
   mutree->Branch("TrackPID", &mu.mPID, "mPID[nparticles]/I");

   cout << "Total Entries = " << nentries << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (jentry%5000==0) cout << jentry << " events finished." << "\r" << flush;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      mu = {0};

      for (int i = 0; i < ConstMax; i++){
         mu.mPt[i] = 0;
         mu.mEta[i] = 0;
         mu.mPhi[i] = 0;
         mu.mEnergy[i] = 0;
         mu.mPID[i] = 0;
      }

      // cout << mu.impactparameter << "\t" << impactParameter << endl;
      int goodparticles = 0;

      for (int i = 0; i < mParticles_; i++){
         TVector3 p(mParticles_mPx[i], mParticles_mPy[i], mParticles_mPz[i]);

         // cout << mParticles_mStatus[i] << endl;
        //  if (abs(mParticles_mId[i]) == 421) cout << mParticles_mId[i] << "\t" << mParticles_mStatus[i] << endl;

         if (mParticles_mStatus[i] != 1) continue; //Only final state particles saved out

         if (abs(p.Eta()) > 1.) continue;
         if (p.Pt() < 0.2) continue;

         goodparticles++;

         mu.mPt[goodparticles-1] = p.Pt();
         mu.mEta[goodparticles-1] = p.Eta();
         mu.mPhi[goodparticles-1] = p.Phi();
         mu.mEnergy[goodparticles-1] = mParticles_mEnergy[i];
         mu.mPID[goodparticles-1] = mParticles_mId[i];

         // if (mParticles_mEnergy[i] > 30.) cout << jentry << "\t" << impactParameter << "\t" << i << "\t" << mParticles_mId[i] << "\t" << mParticles_mEnergy[i] << "\t" << p.Mag() << "\t" << p.Pz() << "\t" << p.Pt() << "\t" << p.Eta() << "\t" << p.Phi() << endl;
      }

      mu.nparticles = goodparticles;

      // cout << "Entry = " << jentry << "\t" << goodparticles << endl;

      mutree->Fill();
   }

   cout << endl;
   
   mutree->Write();
}
