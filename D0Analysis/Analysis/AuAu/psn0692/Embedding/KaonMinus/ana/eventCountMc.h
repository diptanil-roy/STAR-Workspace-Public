//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 18 13:30:55 2016 by ROOT version 5.34/30
// from TTree eventCount/eventCount
// found on file: ../Run14_SL16d_Embedding_KaonMinus_Oct17.McAna.root
//////////////////////////////////////////////////////////

#ifndef eventCountMc_h
#define eventCountMc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class eventCountMc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         runId;
   Float_t         eventId;
   Float_t         mcVx;
   Float_t         mcVy;
   Float_t         mcVz;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         vzVpd;
   Float_t         centrality;
   Float_t         gRefMult;
   Float_t         posRefMult;
   Float_t         negRefMult;
   Float_t         zdc;
   Float_t         bbc;
   Float_t         nMcTracks;
   Float_t         nRTracks;
   Float_t         magField;
   Float_t         t0;
   Float_t         t1;
   Float_t         t2;
   Float_t         t3;
   Float_t         t4;
   Float_t         t5;

   // List of branches
   TBranch        *b_runId;   //!
   TBranch        *b_eventId;   //!
   TBranch        *b_mcVx;   //!
   TBranch        *b_mcVy;   //!
   TBranch        *b_mcVz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vzVpd;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_gRefMult;   //!
   TBranch        *b_posRefMult;   //!
   TBranch        *b_negRefMult;   //!
   TBranch        *b_zdc;   //!
   TBranch        *b_bbc;   //!
   TBranch        *b_nMcTracks;   //!
   TBranch        *b_nRTracks;   //!
   TBranch        *b_magField;   //!
   TBranch        *b_t0;   //!
   TBranch        *b_t1;   //!
   TBranch        *b_t2;   //!
   TBranch        *b_t3;   //!
   TBranch        *b_t4;   //!
   TBranch        *b_t5;   //!

   eventCountMc(TTree *tree=0);
   virtual ~eventCountMc();
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

#ifdef eventCountMc_cxx
eventCountMc::eventCountMc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Run14_SL16d_Embedding_KaonMinus_Oct17.McAna.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../Run14_SL16d_Embedding_KaonMinus_Oct17.McAna.root");
      }
      f->GetObject("eventCount",tree);

   }
   Init(tree);
}

eventCountMc::~eventCountMc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t eventCountMc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t eventCountMc::LoadTree(Long64_t entry)
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

void eventCountMc::Init(TTree *tree)
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

   fChain->SetBranchAddress("runId", &runId, &b_runId);
   fChain->SetBranchAddress("eventId", &eventId, &b_eventId);
   fChain->SetBranchAddress("mcVx", &mcVx, &b_mcVx);
   fChain->SetBranchAddress("mcVy", &mcVy, &b_mcVy);
   fChain->SetBranchAddress("mcVz", &mcVz, &b_mcVz);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("vzVpd", &vzVpd, &b_vzVpd);
   fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
   fChain->SetBranchAddress("gRefMult", &gRefMult, &b_gRefMult);
   fChain->SetBranchAddress("posRefMult", &posRefMult, &b_posRefMult);
   fChain->SetBranchAddress("negRefMult", &negRefMult, &b_negRefMult);
   fChain->SetBranchAddress("zdc", &zdc, &b_zdc);
   fChain->SetBranchAddress("bbc", &bbc, &b_bbc);
   fChain->SetBranchAddress("nMcTracks", &nMcTracks, &b_nMcTracks);
   fChain->SetBranchAddress("nRTracks", &nRTracks, &b_nRTracks);
   fChain->SetBranchAddress("magField", &magField, &b_magField);
   fChain->SetBranchAddress("t0", &t0, &b_t0);
   fChain->SetBranchAddress("t1", &t1, &b_t1);
   fChain->SetBranchAddress("t2", &t2, &b_t2);
   fChain->SetBranchAddress("t3", &t3, &b_t3);
   fChain->SetBranchAddress("t4", &t4, &b_t4);
   fChain->SetBranchAddress("t5", &t5, &b_t5);
   Notify();
}

Bool_t eventCountMc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void eventCountMc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t eventCountMc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef eventCountMc_cxx
