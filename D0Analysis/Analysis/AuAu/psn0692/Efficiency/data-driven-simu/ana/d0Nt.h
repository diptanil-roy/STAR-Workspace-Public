//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 15 08:55:02 2015 by ROOT version 5.34/09
// from TChain nt/
//////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
#ifndef d0Nt_h
#define d0Nt_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class d0Nt {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         cent;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         pid;
   Float_t         w;
   Float_t         m;
   Float_t         pt;
   Float_t         eta;
   Float_t         y;
   Float_t         phi;
   Float_t         v0x;
   Float_t         v0y;
   Float_t         v0z;
   Float_t         rM;
   Float_t         rPt;
   Float_t         rEta;
   Float_t         rY;
   Float_t         rPhi;
   Float_t         reco;
   Float_t         dca12;
   Float_t         decayLength;
   Float_t         dcaD0ToPv;
   Float_t         cosTheta;
   Float_t         angle12;
   Float_t         cosThetaStar;
   Float_t         kM;
   Float_t         kPt;
   Float_t         kEta;
   Float_t         kY;
   Float_t         kPhi;
   Float_t         kDca;
   Float_t         kRM;
   Float_t         kRPt;
   Float_t         kREta;
   Float_t         kRY;
   Float_t         kRPhi;
   Float_t         kRVx;
   Float_t         kRVy;
   Float_t         kRVz;
   Float_t         kRDca;
   Float_t         kRSDca;
   Float_t         kRDcaXY;
   Float_t         kRDcaZ;
   Float_t         kTpc;
   Float_t         pM;
   Float_t         pPt;
   Float_t         pEta;
   Float_t         pY;
   Float_t         pPhi;
   Float_t         pDca;
   Float_t         pRM;
   Float_t         pRPt;
   Float_t         pREta;
   Float_t         pRY;
   Float_t         pRPhi;
   Float_t         pRVx;
   Float_t         pRVy;
   Float_t         pRVz;
   Float_t         pRDca;
   Float_t         pRSDca;
   Float_t         pRDcaXY;
   Float_t         pRDcaZ;
   Float_t         pTpc;
   Float_t         kHft;
   Float_t         pHft;
   Float_t         primaryPt;
   Float_t         primaryY;
   Float_t         primaryPhi;
   // List of branches
   TBranch        *b_cent;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_w;   //!
   TBranch        *b_m;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_y;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_v0x;   //!
   TBranch        *b_v0y;   //!
   TBranch        *b_v0z;   //!
   TBranch        *b_rM;   //!
   TBranch        *b_rPt;   //!
   TBranch        *b_rEta;   //!
   TBranch        *b_rY;   //!
   TBranch        *b_rPhi;   //!
   TBranch        *b_reco;   //!
   TBranch        *b_dca12;   //!
   TBranch        *b_decayLength;   //!
   TBranch        *b_dcaD0ToPv;   //!
   TBranch        *b_cosTheta;   //!
   TBranch        *b_angle12;   //!
   TBranch        *b_cosThetaStar;   //!
   TBranch        *b_kM;   //!
   TBranch        *b_kPt;   //!
   TBranch        *b_kEta;   //!
   TBranch        *b_kY;   //!
   TBranch        *b_kPhi;   //!
   TBranch        *b_kDca;   //!
   TBranch        *b_kRM;   //!
   TBranch        *b_kRPt;   //!
   TBranch        *b_kREta;   //!
   TBranch        *b_kRY;   //!
   TBranch        *b_kRPhi;   //!
   TBranch        *b_kRVx;   //!
   TBranch        *b_kRVy;   //!
   TBranch        *b_kRVz;   //!
   TBranch        *b_kRDca;   //!
   TBranch        *b_kRSDca;   //!
   TBranch        *b_kRDcaXY;   //!
   TBranch        *b_kRDcaZ;   //!
   TBranch        *b_kTpc;   //!
   TBranch        *b_pM;   //!
   TBranch        *b_pPt;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pY;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pDca;   //!
   TBranch        *b_pRM;   //!
   TBranch        *b_pRPt;   //!
   TBranch        *b_pREta;   //!
   TBranch        *b_pRY;   //!
   TBranch        *b_pRPhi;   //!
   TBranch        *b_pRVx;   //!
   TBranch        *b_pRVy;   //!
   TBranch        *b_pRVz;   //!
   TBranch        *b_pRDca;   //!
   TBranch        *b_pRSDca;   //!
   TBranch        *b_pRDcaXY;   //!
   TBranch        *b_pRDcaZ;   //!
   TBranch        *b_pTpc;   //!
   TBranch        *b_kHft;   //!
   TBranch        *b_pHft;   //!
   TBranch        *b_primaryPt;
   TBranch        *b_primaryY;
   TBranch        *b_primaryPhi;

   d0Nt(std::string file);
   virtual ~d0Nt();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t GetEntries() const { return fChain->GetEntries();}
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef d0Nt_cxx
d0Nt::d0Nt(std::string file) : fChain(0) 
{

  TChain * chain = new TChain("nt","");
  std::string dirFile = file;
  if (dirFile.find(".list")!=std::string::npos)  {
    ifstream inputStream(dirFile.c_str());
    if(!(inputStream.good())) {
      std::cout << "ERROR: Cannot open list file " << dirFile << std::endl;
    }
    char line[512];
    int nFile=0;
    std::string ltest;
    while (inputStream.good()) {
      inputStream.getline(line,512);
      std::string aFile = line;      
      if (inputStream.good() && aFile.find(".root")!=std::string::npos) {
        TFile *ftmp = TFile::Open(line);
        if(ftmp && ftmp->IsOpen() && ftmp->GetNkeys()) {
          std::cout << " Read in root file " << line << std::endl;
          chain->Add(line);
          nFile++;
        }
      }
    }
    std::cout << " Total " << nFile << " files have been read in. " << std::endl;
  } else if (dirFile.find(".root")!=std::string::npos)  {
    chain->Add(dirFile.c_str());
  } else {
    std::cout << " No good input file to read ... " << std::endl;
  }
  //chain->Add(file.c_str());
  TTree* tree = chain;

   Init(tree);
}

d0Nt::~d0Nt()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t d0Nt::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t d0Nt::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void d0Nt::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("cent", &cent, &b_cent);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("pid", &pid, &b_pid);
   fChain->SetBranchAddress("w", &w, &b_w);
   fChain->SetBranchAddress("m", &m, &b_m);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("v0x", &v0x, &b_v0x);
   fChain->SetBranchAddress("v0y", &v0y, &b_v0y);
   fChain->SetBranchAddress("v0z", &v0z, &b_v0z);
   fChain->SetBranchAddress("rM", &rM, &b_rM);
   fChain->SetBranchAddress("rPt", &rPt, &b_rPt);
   fChain->SetBranchAddress("rEta", &rEta, &b_rEta);
   fChain->SetBranchAddress("rY", &rY, &b_rY);
   fChain->SetBranchAddress("rPhi", &rPhi, &b_rPhi);
   fChain->SetBranchAddress("reco", &reco, &b_reco);
   fChain->SetBranchAddress("dca12", &dca12, &b_dca12);
   fChain->SetBranchAddress("decayLength", &decayLength, &b_decayLength);
   fChain->SetBranchAddress("dcaD0ToPv", &dcaD0ToPv, &b_dcaD0ToPv);
   fChain->SetBranchAddress("cosTheta", &cosTheta, &b_cosTheta);
   fChain->SetBranchAddress("angle12", &angle12, &b_angle12);
   fChain->SetBranchAddress("cosThetaStar", &cosThetaStar, &b_cosThetaStar);
   fChain->SetBranchAddress("kM", &kM, &b_kM);
   fChain->SetBranchAddress("kPt", &kPt, &b_kPt);
   fChain->SetBranchAddress("kEta", &kEta, &b_kEta);
   fChain->SetBranchAddress("kY", &kY, &b_kY);
   fChain->SetBranchAddress("kPhi", &kPhi, &b_kPhi);
   fChain->SetBranchAddress("kDca", &kDca, &b_kDca);
   fChain->SetBranchAddress("kRM", &kRM, &b_kRM);
   fChain->SetBranchAddress("kRPt", &kRPt, &b_kRPt);
   fChain->SetBranchAddress("kREta", &kREta, &b_kREta);
   fChain->SetBranchAddress("kRY", &kRY, &b_kRY);
   fChain->SetBranchAddress("kRPhi", &kRPhi, &b_kRPhi);
   fChain->SetBranchAddress("kRVx", &kRVx, &b_kRVx);
   fChain->SetBranchAddress("kRVy", &kRVy, &b_kRVy);
   fChain->SetBranchAddress("kRVz", &kRVz, &b_kRVz);
   fChain->SetBranchAddress("kRDca", &kRDca, &b_kRDca);
   fChain->SetBranchAddress("kRSDca", &kRSDca, &b_kRSDca);
   fChain->SetBranchAddress("kRDcaXY", &kRDcaXY, &b_kRDcaXY);
   fChain->SetBranchAddress("kRDcaZ", &kRDcaZ, &b_kRDcaZ);
   fChain->SetBranchAddress("kTpc", &kTpc, &b_kTpc);
   fChain->SetBranchAddress("pM", &pM, &b_pM);
   fChain->SetBranchAddress("pPt", &pPt, &b_pPt);
   fChain->SetBranchAddress("pEta", &pEta, &b_pEta);
   fChain->SetBranchAddress("pY", &pY, &b_pY);
   fChain->SetBranchAddress("pPhi", &pPhi, &b_pPhi);
   fChain->SetBranchAddress("pDca", &pDca, &b_pDca);
   fChain->SetBranchAddress("pRM", &pRM, &b_pRM);
   fChain->SetBranchAddress("pRPt", &pRPt, &b_pRPt);
   fChain->SetBranchAddress("pREta", &pREta, &b_pREta);
   fChain->SetBranchAddress("pRY", &pRY, &b_pRY);
   fChain->SetBranchAddress("pRPhi", &pRPhi, &b_pRPhi);
   fChain->SetBranchAddress("pRVx", &pRVx, &b_pRVx);
   fChain->SetBranchAddress("pRVy", &pRVy, &b_pRVy);
   fChain->SetBranchAddress("pRVz", &pRVz, &b_pRVz);
   fChain->SetBranchAddress("pRDca", &pRDca, &b_pRDca);
   fChain->SetBranchAddress("pRSDca", &pRSDca, &b_pRSDca);
   fChain->SetBranchAddress("pRDcaXY", &pRDcaXY, &b_pRDcaXY);
   fChain->SetBranchAddress("pRDcaZ", &pRDcaZ, &b_pRDcaZ);
   fChain->SetBranchAddress("pTpc", &pTpc, &b_pTpc);
   fChain->SetBranchAddress("kHft", &kHft, &b_kHft);
   fChain->SetBranchAddress("pHft", &pHft, &b_pHft);
   fChain->SetBranchAddress("primaryPt", &primaryPt, &b_primaryPt);
   fChain->SetBranchAddress("primaryY", &primaryY, &b_primaryY);
   fChain->SetBranchAddress("primaryPhi", &primaryPhi, &b_primaryPhi);
   Notify();
}

Bool_t d0Nt::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void d0Nt::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t d0Nt::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef d0Nt_cxx
