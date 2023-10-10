//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 18 13:31:09 2016 by ROOT version 5.34/30
// from TTree tracks/
// found on file: ../Run14_SL16d_Embedding_KaonMinus_Oct17.McAna.root
//////////////////////////////////////////////////////////

#ifndef tracksMc_h
#define tracksMc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class tracksMc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         centrality;
   Float_t         pt;
   Float_t         p;
   Float_t         eta;
   Float_t         y;
   Float_t         phi;
   Float_t         geantId;
   Float_t         eventGenLabel;
   Float_t         startVtxX;
   Float_t         startVtxY;
   Float_t         startVtxZ;
   Float_t         stopVtxX;
   Float_t         stopVtxY;
   Float_t         stopVtxZ;
   Float_t         gPt;
   Float_t         gEta;
   Float_t         gPhi;
   Float_t         nFit;
   Float_t         nMax;
   Float_t         nCom;
   Float_t         nDedx;
   Float_t         dedx;
   Float_t         nSigPi;
   Float_t         nSigK;
   Float_t         dca;
   Float_t         dcaXY;
   Float_t         dcaZ;
   Float_t         hftTopo;
   Float_t         hftTruth;
   Float_t         trkMap0;
   Float_t         trkMap1;

   // List of branches
   TBranch        *b_centrality;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_p;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_y;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_geantId;   //!
   TBranch        *b_eventGenLabel;   //!
   TBranch        *b_startVtxX;   //!
   TBranch        *b_startVtxY;   //!
   TBranch        *b_startVtxZ;   //!
   TBranch        *b_stopVtxX;   //!
   TBranch        *b_stopVtxY;   //!
   TBranch        *b_stopVtxZ;   //!
   TBranch        *b_gPt;   //!
   TBranch        *b_gEta;   //!
   TBranch        *b_gPhi;   //!
   TBranch        *b_nFit;   //!
   TBranch        *b_nMax;   //!
   TBranch        *b_nCom;   //!
   TBranch        *b_nDedx;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_nSigPi;   //!
   TBranch        *b_nSigK;   //!
   TBranch        *b_dca;   //!
   TBranch        *b_dcaXY;   //!
   TBranch        *b_dcaZ;   //!
   TBranch        *b_hftTopo;   //!
   TBranch        *b_hftTruth;   //!
   TBranch        *b_trkMap0;   //!
   TBranch        *b_trkMap1;   //!

   tracksMc(TTree *tree=0);
   virtual ~tracksMc();
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

#ifdef tracksMc_cxx
tracksMc::tracksMc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Run14_SL16d_Embedding_KaonMinus_Oct17.McAna.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../Run14_SL16d_Embedding_KaonMinus_Oct17.McAna.root");
      }
      f->GetObject("tracks",tree);

   }
   Init(tree);
}

tracksMc::~tracksMc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tracksMc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tracksMc::LoadTree(Long64_t entry)
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

void tracksMc::Init(TTree *tree)
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

   fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("geantId", &geantId, &b_geantId);
   fChain->SetBranchAddress("eventGenLabel", &eventGenLabel, &b_eventGenLabel);
   fChain->SetBranchAddress("startVtxX", &startVtxX, &b_startVtxX);
   fChain->SetBranchAddress("startVtxY", &startVtxY, &b_startVtxY);
   fChain->SetBranchAddress("startVtxZ", &startVtxZ, &b_startVtxZ);
   fChain->SetBranchAddress("stopVtxX", &stopVtxX, &b_stopVtxX);
   fChain->SetBranchAddress("stopVtxY", &stopVtxY, &b_stopVtxY);
   fChain->SetBranchAddress("stopVtxZ", &stopVtxZ, &b_stopVtxZ);
   fChain->SetBranchAddress("gPt", &gPt, &b_gPt);
   fChain->SetBranchAddress("gEta", &gEta, &b_gEta);
   fChain->SetBranchAddress("gPhi", &gPhi, &b_gPhi);
   fChain->SetBranchAddress("nFit", &nFit, &b_nFit);
   fChain->SetBranchAddress("nMax", &nMax, &b_nMax);
   fChain->SetBranchAddress("nCom", &nCom, &b_nCom);
   fChain->SetBranchAddress("nDedx", &nDedx, &b_nDedx);
   fChain->SetBranchAddress("dedx", &dedx, &b_dedx);
   fChain->SetBranchAddress("nSigPi", &nSigPi, &b_nSigPi);
   fChain->SetBranchAddress("nSigK", &nSigK, &b_nSigK);
   fChain->SetBranchAddress("dca", &dca, &b_dca);
   fChain->SetBranchAddress("dcaXY", &dcaXY, &b_dcaXY);
   fChain->SetBranchAddress("dcaZ", &dcaZ, &b_dcaZ);
   fChain->SetBranchAddress("hftTopo", &hftTopo, &b_hftTopo);
   fChain->SetBranchAddress("hftTruth", &hftTruth, &b_hftTruth);
   fChain->SetBranchAddress("trkMap0", &trkMap0, &b_trkMap0);
   fChain->SetBranchAddress("trkMap1", &trkMap1, &b_trkMap1);
   Notify();
}

Bool_t tracksMc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tracksMc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tracksMc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tracksMc_cxx
