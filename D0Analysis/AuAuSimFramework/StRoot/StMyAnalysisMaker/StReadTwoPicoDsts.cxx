// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StReadTwoPicoDsts.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMcTrack.h"
#include "StRoot/StPicoEvent/StPicoMcVertex.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StReadTwoPicoDsts)

//________________________________________________________________________
StReadTwoPicoDsts::StReadTwoPicoDsts(const char* name, StPicoDstReader *Reader1, StPicoDstReader *Reader2, const char* vertexfilename = "", const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;

  mPicoDstMakerMC = 0x0;
  mPicoDstMC = 0x0;
  mPicoEventMC = 0x0;

  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StReadTwoPicoDsts::fRunFlagEnum
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fCorrJetPt = kFALSE;
  fMinPtJet = 0.0;
  fTrackBias = 0.0;
  fTowerBias = 0.0;
  fJetRad = 0.4;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;
  fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0;  fTowerPhiMaxCut = 2.0*TMath::Pi();
  fCentralityScaled = 0.; ref16 = -99; ref9 = -99; // FIXME - maybe not make global
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  fRho = 0x0;
  fRhoVal = 0;
  mEmcPosition = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  mVertexFileName = vertexfilename;
  // DataFileName = inFile1;
  // MCFileName = inFile2;

  mPicoDstReader = Reader1;
  mPicoDstReaderMC = Reader2;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

}

//
//________________________________________________________________________
StReadTwoPicoDsts::~StReadTwoPicoDsts()
{ /*  */
  // destructor
  if(hCentrality)  delete hCentrality;
  if(hMultiplicity)delete hMultiplicity;
  if(hJetPt)       delete hJetPt;
  if(hJetCorrPt)   delete hJetCorrPt;

  if(mEmcPosition) delete mEmcPosition;

//______________________Test Code #FIXME ________________________________
  if(hTrackPt) delete hTrackPt;
  if(hTrackEt) delete hTrackEt;

}

//
//________________________________________________________________________
Int_t StReadTwoPicoDsts::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms
  DeclareHistograms();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);


  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StReadTwoPicoDsts::Finish() { 
  cout << "StReadTwoPicoDsts::Finish()\n";
  vertexfile.close();

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }
  
  cout<<"End of StReadTwoPicoDsts::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StReadTwoPicoDsts::DeclareHistograms() {
  // binning for cent histograms
  int nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100.;
  }

  // histograms
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);
  
  //___________________TEST Code #FIXME____________________________________

  
  hTrackPt = new TH2F("hTrackPt", "No. of tracks", 200, -2.0, 2.0, 200, -5.0, 5.0);
  hTrackEt = new TH2F("hTrackEt", "No. of tracks", 200, -2.0, 2.0, 200, -5.0, 5.0);

  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 200, -100., 100.);

  hTPCPt = new TH1F("hTPCPt", "TPC hits vs Pt", 32, -0.5, 15.5);
  hHFTPt = new TH1F("hHFTPt", "HFT hits vs Pt", 32, -0.5, 15.5);
  hMCPt = new TH1F("hMCPt", "MC hits vs Pt", 32, -0.5, 15.5);
  hGenPt = new TH1F("hGenPt", "Gen vs Pt", 32, -0.5, 15.5);

  hTPCEta = new TH1F("hTPCEta", "TPC hits vs Eta", 40, -1.0, 1.0);
  hHFTEta = new TH1F("hHFTEta", "HFT hits vs Eta", 40, -1.0, 1.0);
  hMCEta = new TH1F("hMCEta", "MC hits vs Eta", 40, -1.0, 1.0);
  hGenEta = new TH1F("hGenEta", "Gen vs Eta", 40, -1.0, 1.0);

  hTPCPhi = new TH1F("hTPCPhi", "TPC hits vs Phi", 80, -10, 10);
  hHFTPhi = new TH1F("hHFTPhi", "HFT hits vs Phi", 80, -10, 10);
  hMCPhi = new TH1F("hMCPhi", "MC hits vs Phi", 80, -10, 10);
  hGenPhi = new TH1F("hGenPhi", "Gen vs Phi", 80, -10, 10);


}
//
// write histograms
//_____________________________________________________________________________
void StReadTwoPicoDsts::WriteHistograms() {
  // writing of histograms done here
  // hCentrality->Write();
  // hMultiplicity->Write();
  // hJetPt->Write();
  // hJetCorrPt->Write();
  // hTrackPt->SetOption("LEGO2");
  // hTrackEt->SetOption("LEGO2");
  // hTrackPt->Write();
  // hTrackEt->Write();

  hEventZVertex->Write();

  hTPCPt->Write();
  hTPCEta->Write();
  hTPCPhi->Write();

  hHFTPt->Write();
  hHFTEta->Write();
  hHFTPhi->Write();

  hMCPt->Write();
  hMCEta->Write();
  hMCPhi->Write();

  hGenPt->Write();
  hGenEta->Write();
  hGenPhi->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StReadTwoPicoDsts::Clear(Option_t *opt) {
  fJets->Clear();
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StReadTwoPicoDsts::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // DataFile = TFile::Open(DataFileName);
  // if (!DataFile->IsOpen()){
  //   return kStWarn;
  // }

  // MCFile = TFile::Open(MCFileName);
  // if (!DataFile->IsOpen()){
  //   return kStWarn;
  // }

  // cout << DataFileName << endl;
  // cout << MCFileName << endl;

  // mPicoDstReader = new StPicoDstReader(DataFileName);
  // mPicoDstReaderMC = new StPicoDstReader(MCFileName);

  mPicoDstReader->readPicoEvent(20);
  mPicoDstReaderMC->readPicoEvent(1);

  // // construct PicoDst object from maker
  // mPicoDst = (StPicoDst *)DataFile->Get("PicoDst");
  // if(!mPicoDst) {
  //   LOG_WARN << " No PicoDst! Skip! " << endm;
  //   return kStWarn;
  // }

  // // create pointer to PicoEvent 
  // mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  // if(!mPicoEvent) {
  //   LOG_WARN << " No PicoEvent! Skip! " << endm;
  //   return kStWarn;
  // }

  // cout << mPicoEvent->runId() << "\t" << mPicoEvent->eventId() << endl;

  // // construct PicoDst object from maker
  // mPicoDstMC = (StPicoDst *)MCFile->Get("PicoDst");
  // if(!mPicoDstMC) {
  //   LOG_WARN << " No PicoDst! Skip! " << endm;
  //   return kStWarn;
  // }

  // // create pointer to PicoEvent 
  // mPicoEventMC = static_cast<StPicoEvent*>(mPicoDstMC->event());
  // if(!mPicoEventMC) {
  //   LOG_WARN << " No PicoEvent! Skip! " << endm;
  //   return kStWarn;
  // }

  // cout << mPicoEventMC->runId() << "\t" << mPicoEventMC->eventId() << endl;


  return kStOK;
}

//
//
//_____________________________________________________________________________________________
// void StReadTwoPicoDsts::RunJets()
// {
//   // cache the leading + subleading jets within acceptance
//   // first parameter is Jet Maker name, 2nd is Rho Parameter: fRho
//   if(fCorrJetPt) {
//     fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
//     fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
//   } else {
//     fLeadingJet = GetLeadingJet(fJetMakerName);
//     fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
//   }

//   // ====================== Jet loop below ============================
//   // loop over Jets in the event: initialize some parameter variables
//   Int_t ijethi = -1;
//   Double_t highestjetpt = 0.0;
//   Int_t njets = fJets->GetEntries();
//   if (njets!=0) cout << "Number of jets = " << njets << endl;
//   // loop over jets
//   for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
//     // get jet pointer
//     StJet *jet = static_cast<StJet*>(fJets->At(ijet));
//     if(!jet) continue;

//     // get some jet parameters
//     double jetarea = jet->Area();
//     double jetpt = jet->Pt();
//     double corrjetpt = jet->Pt() - jetarea*fRhoVal;
//     double jetE = jet->E();
//     double jetEta = jet->Eta();
//     double jetPhi = jet->Phi();
//     double jetNEF = jet->NEF();
//     //double jetM = jet->M(); //#FIXME
    

//     // get nTracks and maxTrackPt
//     double maxtrackpt = jet->GetMaxTrackPt();
//     double NtrackConstit = jet->GetNumberOfTracks();

//     // get highest Pt jet in event (leading jet)
//     if(highestjetpt < jetpt){
//       ijethi = ijet;
//       highestjetpt = jetpt;
//     }

//     // fill some basic histos
//     hJetPt->Fill(jetpt);
//     hJetCorrPt->Fill(corrjetpt);

//    //______________ TEST Code #FIXME________________________

//    // hTrackPt->Fill(jetpt, jetEta, jetPhi);
//    // hTrackEt->Fill(jetpt, jetEta, jetPhi);

//     // TEST - when using constituent subtractor
//     vector<fastjet::PseudoJet> fConstituents = jet->GetJetConstituents();
//     for(UInt_t ic = 0; ic < fConstituents.size(); ++ic) {
//       // get user defined index
//       Int_t uid = fConstituents[ic].user_index();
//       double cpt = fConstituents[ic].perp();
//       double ceta = fConstituents[ic].eta();
//       double cphi = fConstituents[ic].phi();
//       cout<<"ic = "<<ic<<", uid = "<<uid<<", cpt = "<<cpt<<", ceta = "<<ceta<<", cphi = "<<cphi<<endl;
//     }

//     cout << ijet << "\t" << jet->GetNumberOfTracks() << endl;

//     // get jet constituents: loop over constituent tracks
//     for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
//       int trackid = jet->TrackAt(itrk); 



//       // get jet track pointer
//       StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
//       if(!trk){ continue; }

//       // get momentum vector
//       TVector3 mTrkMom;
//       if(doUsePrimTracks) {
//         if(!(trk->isPrimary())) continue; // check if primary
//         // get primary track vector
//         mTrkMom = trk->pMom();
//       } else {
//         // get global track vector
//         mTrkMom = trk->gMom(mVertex, Bfield);
//       }

//       // track variables
//       double pt = mTrkMom.Perp();
//       double phi = mTrkMom.Phi();
//       double eta = mTrkMom.PseudoRapidity();
//       double px = mTrkMom.x();
//       double py = mTrkMom.y();
//       double pz = mTrkMom.z();
//       short charge = trk->charge();

//     } // track constit loop

//     // loop over constituents towers
//     for(int itow = 0; itow < jet->GetNumberOfClusters(); itow++) {
//       int towerid = jet->ClusterAt(itow);

//       // get jet tower pointer
//       StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towerid));
//       if(!tow){ continue; }

//       // tower ID shifted by +1 from array index
//       int towID = itow + 1;
    
//     } // tower constit loop

//   } // jet loop

// }

//
//
//________________________________________________________________________
void StReadTwoPicoDsts::RunTracks()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV
  unsigned int ntracks = mPicoDst->numberOfTracks();
  unsigned int nmctracks = mPicoDst->numberOfMcTracks();

  cout << ntracks << "\t" << nmctracks << endl;

  for(unsigned short iTracks = 0; iTracks < nmctracks; iTracks++){
    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(iTracks));
    if(mctrk->idVtxStart() > 1) continue;
    int mcgepid = mctrk->id();
    if (mcgepid == 8 || mcgepid == 9){
      hGenPt->Fill(mctrk->pt());
      hGenEta->Fill(mctrk->eta());
      hGenPhi->Fill(mctrk->fourMomentum().Phi());
    }
  }

  // loop over ALL tracks in PicoDst 
  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // primary track switch: get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = trk->pMom();
    } else {
      // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();
    double px = mTrkMom.x();
    double py = mTrkMom.y();
    double pz = mTrkMom.z();
    short charge = trk->charge();

    int mctrkid = trk->idTruth();

    StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(mctrkid));
    if(!mctrk) continue;
    if(mctrk->idVtxStart() > 1) continue;

    int mcgepid = mctrk->id();
    if (mcgepid == 8 || mcgepid == 9){
      hMCPt->Fill(mctrk->pt());
      hMCEta->Fill(mctrk->eta());
      hMCPhi->Fill(mctrk->fourMomentum().Phi());
    }

    double zpi = trk->nSigmaPion();

    if (abs(zpi) > 2) continue;

    if (pt < 0.2) continue;

    if (abs(eta) > 1) continue;

    if (trk->nHitsDedx() < 15) continue;

    if (float(trk->nHitsDedx())/float(trk->nHitsMax()) < 0.52) continue;

    if (trk->gDCA(mPicoEvent->primaryVertex()).Mag() > 3.) continue;

    if (trk->chi2() > 3) continue;


    hTPCPt->Fill(pt);
    hTPCEta->Fill(eta);
    hTPCPhi->Fill(phi);

    if (trk->isHFTTrack()){
      hHFTPt->Fill(pt);
      hHFTEta->Fill(eta);
      hHFTPhi->Fill(phi);
    }

//     int mctrkid = trk->idTruth();

//     StPicoMcTrack *mctrk = static_cast<StPicoMcTrack*>(mPicoDst->mcTrack(mctrkid));

//     int mcgepid = mctrk->id();
//     if (mcgepid == 8 || mcgepid == 9) cout << "Pions here" << endl;
//     else continue;

//     // if (trk->isHFTTrack()) cout << iTracks << "\t" << pt << "\t" << phi << "\t" << eta << endl;

//     if (trk->isHFTTrack()) cout << "HFT" << endl; 
//     // if (trk->hasIstHit()) cout << "IST" << endl; 
//     // if (trk->hasPxl1Hit() || trk->hasPxl2Hit()) cout << "PXL" << endl; 


// //______________ TEST Code #FIXME________________________

//     hTrackPt->Fill(eta, phi, pt);
//

    // do some things with tracks here
    // ......

    // fill track histograms here
    // .........

  } // track loop

}  // track function

//
//
//________________________________________________________________________
void StReadTwoPicoDsts::RunTowers()
{
  // constants: assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // looping over clusters - STAR: matching already done, get # of clusters and set variables
  unsigned int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();

  // loop over ALL clusters in PicoDst and add to jet //TODO
  for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
    StPicoBEmcPidTraits *cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
    if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

    // cluster and tower ID
    int clusID = cluster->bemcId();  // index in bemc point array
    int towID = cluster->btowId();   // projected tower Id: 1 - 4800
    int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    if(towID < 0) continue;

    // cluster and tower position - from vertex and ID
    TVector3 towPosition = mEmcPosition->getPosFromVertex(mVertex, towID);
    double towPhi = towPosition.Phi();
    double towEta = towPosition.PseudoRapidity();

    // matched track index
    int trackIndex = cluster->trackIndex();
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
    if(!trk) { cout<<"No trk pointer...."<<endl; continue; }
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

  } // BEmc loop

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // get tower pointer
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: tower ID shifted by +1 from array index
    int towerID = itow + 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    double towerPhi = towerPosition.Phi();
    double towerEta = towerPosition.PseudoRapidity();
    int towerADC = tower->adc();
    double towerE = tower->energy();
    double towerEt = towerE / (1.0*TMath::CosH(towerEta)); // THIS should be USED

    // do stuff with towers and fill histograms here
    // ........
    hTrackEt->Fill(towerEta, towerPhi, towerEt);


  } // tower loop

}  // run towers function

//
//
// __________________________________________________________________________________
void StReadTwoPicoDsts::SetSumw2() {
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  hJetPt->Sumw2();
  hJetCorrPt->Sumw2();
}

//
// Function: get relative phi of jet and track (-0.5pi, 1.5pi)
//________________________________________________________________________
Double_t StReadTwoPicoDsts::RelativePhi(Double_t mphi,Double_t vphi) const
{ // function to calculate relative PHI
  double dphi = mphi-vphi;

  // set dphi to operate on adjusted scale
  if(dphi < -0.5*pi) dphi += 2.0*TMath::Pi();
  if(dphi >  1.5*pi) dphi -= 2.0*TMath::Pi();

  // test check
  if( dphi < -0.5*pi || dphi > 1.5*pi )
    Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName());

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}

//
// Function: calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
//_________________________________________________________________________
Double_t StReadTwoPicoDsts::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{
  Double_t pi = 1.0*TMath::Pi();
  Double_t dphi = 1.0*TMath::Abs(EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi < -pi ){
    dphi = dphi + pi;
  } // this assumes we are doing full jets currently 
 
  if(dphi > 1.5*pi) dphi -= 2.0*pi;
  if((dphi > 1.0*pi) && (dphi < 1.5*pi)) dphi -= 1.0*pi;
  if((dphi > 0.5*pi) && (dphi < 1.0*pi)) dphi -= 1.0*pi;
  dphi = 1.0*TMath::Abs(dphi);

  // test check
  if( dphi < 0.0 || dphi > 0.5*pi ) {
    //Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());
    cout<<"dPhi not in range [0, 0.5*Pi]!"<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}

//
//
//_________________________________________________________________________
void StReadTwoPicoDsts::FillEmcTriggers() {
  // number of Emcal Triggers
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // fill for valid triggers
    if(emcTrig->isHT0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
}
//
// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StReadTwoPicoDsts::DoComparison(int myarr[], int elems) {
  //std::cout << "Length of array = " << (sizeof(myarr)/sizeof(*myarr)) << std::endl;
  bool match = kFALSE;

  // loop over specific physics selection array and compare to specific event trigger
  for(int i=0; i<elems; i++) {
    if(mPicoEvent->isTrigger(myarr[i])) match = kTRUE;
    if(match) break;
  }

  return match;
}
