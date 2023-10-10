//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May  3 13:15:31 2022 by ROOT version 6.20/08
// from TTree D0USTree/D0USTree
// found on file: Test.root
//////////////////////////////////////////////////////////

#ifndef Test_h
#define Test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         D0Mass;
   Float_t         D0Pt;
   Float_t         D0Eta;
   Float_t         D0Phi;
   Float_t         PionPt;
   Float_t         PionEta;
   Float_t         PionPhi;
   Float_t         PionDCA;
   Float_t         PionNHitsFit;
   Float_t         PionNSigmaPion;
   Float_t         PionNSigmaKaon;
   Float_t         PionTofBeta;
   Float_t         PionTofYLocal;
   Float_t         KaonPt;
   Float_t         KaonEta;
   Float_t         KaonPhi;
   Float_t         KaonDCA;
   Float_t         KaonNHitsFit;
   Float_t         KaonNSigmaPion;
   Float_t         KaonNSigmaKaon;
   Float_t         KaonTofBeta;
   Float_t         KaonTofYLocal;

   // List of branches
   TBranch        *b_D0mass;   //!
   TBranch        *b_D0pt;   //!
   TBranch        *b_D0eta;   //!
   TBranch        *b_D0phi;   //!
   TBranch        *b_pionpt;   //!
   TBranch        *b_pioneta;   //!
   TBranch        *b_pionphi;   //!
   TBranch        *b_piondca;   //!
   TBranch        *b_pionnhitsfit;   //!
   TBranch        *b_pionnsigmapion;   //!
   TBranch        *b_pionnsigmakaon;   //!
   TBranch        *b_piontofbeta;   //!
   TBranch        *b_ntofpionylocal;   //!
   TBranch        *b_kaonpt;   //!
   TBranch        *b_kaoneta;   //!
   TBranch        *b_kaonphi;   //!
   TBranch        *b_kaondca;   //!
   TBranch        *b_kaonnhitsfit;   //!
   TBranch        *b_kaonnsigmapion;   //!
   TBranch        *b_kaonnsigmakaon;   //!
   TBranch        *b_kaontofbeta;   //!
   TBranch        *b_ntofkaonylocal;   //!

   Test(TTree *tree=0);
   virtual ~Test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Test_cxx
Test::Test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Test.root:/PIDQA_TrackTree");
      dir->GetObject("D0USTree",tree);

   }
   Init(tree);
}

Test::~Test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Test::LoadTree(Long64_t entry)
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

void Test::Init(TTree *tree)
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

   fChain->SetBranchAddress("D0Mass", &D0Mass, &b_D0mass);
   fChain->SetBranchAddress("D0Pt", &D0Pt, &b_D0pt);
   fChain->SetBranchAddress("D0Eta", &D0Eta, &b_D0eta);
   fChain->SetBranchAddress("D0Phi", &D0Phi, &b_D0phi);
   fChain->SetBranchAddress("PionPt", &PionPt, &b_pionpt);
   fChain->SetBranchAddress("PionEta", &PionEta, &b_pioneta);
   fChain->SetBranchAddress("PionPhi", &PionPhi, &b_pionphi);
   fChain->SetBranchAddress("PionDCA", &PionDCA, &b_piondca);
   fChain->SetBranchAddress("PionNHitsFit", &PionNHitsFit, &b_pionnhitsfit);
   fChain->SetBranchAddress("PionNSigmaPion", &PionNSigmaPion, &b_pionnsigmapion);
   fChain->SetBranchAddress("PionNSigmaKaon", &PionNSigmaKaon, &b_pionnsigmakaon);
   fChain->SetBranchAddress("PionTofBeta", &PionTofBeta, &b_piontofbeta);
   fChain->SetBranchAddress("PionTofYLocal", &PionTofYLocal, &b_ntofpionylocal);
   fChain->SetBranchAddress("KaonPt", &KaonPt, &b_kaonpt);
   fChain->SetBranchAddress("KaonEta", &KaonEta, &b_kaoneta);
   fChain->SetBranchAddress("KaonPhi", &KaonPhi, &b_kaonphi);
   fChain->SetBranchAddress("KaonDCA", &KaonDCA, &b_kaondca);
   fChain->SetBranchAddress("KaonNHitsFit", &KaonNHitsFit, &b_kaonnhitsfit);
   fChain->SetBranchAddress("KaonNSigmaPion", &KaonNSigmaPion, &b_kaonnsigmapion);
   fChain->SetBranchAddress("KaonNSigmaKaon", &KaonNSigmaKaon, &b_kaonnsigmakaon);
   fChain->SetBranchAddress("KaonTofBeta", &KaonTofBeta, &b_kaontofbeta);
   fChain->SetBranchAddress("KaonTofYLocal", &KaonTofYLocal, &b_ntofkaonylocal);
   Notify();
}

Bool_t Test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Test_cxx
