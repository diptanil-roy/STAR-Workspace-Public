/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
 */

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "TStopwatch.h"
// #include "TSystem.h"
// #include "TMemStat.h"

using namespace std;

void setDecayChannels(int const mMode);
bool GetD0FromBDecay(int& kf, TLorentzVector& b, TVector3& v0_D0, TClonesArray& daughters);
void decayAndFill(int const kf, TLorentzVector* b, TVector3 const v0_D0, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
TVector3 getVertex(int centrality);
bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom);
float matchCorrection(int iParticleIndex, int iEtaIndex, int iVzIndex, int iPhiIndex,  int iCent,  float pt);
bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);
bool reconstructD0(int const centrality, TLorentzVector const& mom);
void bookObjects();
void write();
int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);
int getPhiIndexDca(double);

int getPtIndexHftRatio(double);
int getEtaIndexHftRatio(double);
int getVzIndexHftRatio(double);
int getPhiIndexHftRatio(double);

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;
float primaryPt;
float primaryY;
float primaryPhi;

//TF1* fKaonMomResolution = NULL;
//TF1* fPionMomResolution = NULL;
TF1* fKaonPlusMomResolution = NULL;
TF1* fPionPlusMomResolution = NULL;
TF1* fKaonMinusMomResolution = NULL;
TF1* fPionMinusMomResolution = NULL;
TF1* fWeightFunction = NULL;
TGraph* grEff[3];
const Int_t nParticles = 2;
const Int_t nCentHftRatio = 9;

// HFT ratio binning
const Int_t nEtasHftRatio = 10;
const Int_t nVzsHftRatio = 6;
const Int_t nPtBinsHftRatio = 50;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] =
{
    -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
};
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] =
{
    -6.0e4, -4.0e4, -2.0e4, 0.0, 2.0e4, 4.0e4, 6.0e4
};
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] =
{
    0.2, 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 ,
    1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
    2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 , 3.2 , 3.4 , 3.6 , 3.8 ,
    4. , 4.2 , 4.4 , 4.6 , 4.8 , 5.0 , 5.2 , 5.4 , 5.6 , 5.8 ,
    6. , 6.5 , 7.0 , 7.5 , 8.0 , 8.5 , 9.0 , 9.5 , 10. , 10.5,
    11.0 , 11.5 , 12.0
};
const Int_t nPhisHftRatio = 11;
const Double_t PhiEdgeHftRatio[nPhisHftRatio + 1] =
{
    -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector  // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
};

const Int_t nPhisDca = 11;
const Double_t PhiEdgeDca[nPhisDca + 1] =
{
    -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector  // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
};

// DCA binning
int const nVzsDca = 4;
float const VzEdgeDca[nVzsDca + 1] = {   -6.e4, -3.e4, 0, 3.e4, 6.e4};

int const nEtasDca = 5;
float const EtaEdgeDca[nEtasDca + 1] = { -1.0, -0.6, -0.2, 0.2, 0.6, 1.0};

const Int_t nPtBinsDca = 19;
const Double_t ptEdgeDca[nPtBinsDca + 1] =
{
    0.3, 0.4, 0.5,
    0.6,  0.7 , 0.8 , 0.9 ,
    1. , 1.25 , 1.5 , 1.75 , 2.  , 2.25 , 2.5 , 2.75 , 3.0 , 3.5,
    // 3. , 3.5 , 4.  , 4.5 , 5. , 6. , 8.0 , 10. , 12.0
    4.  , 6. , 12.0
};

TH1D* h1Vz[nCentHftRatio];

TH1D* hHftRatioCorrect[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][1];
TH1D* hHftRatio1[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][nCentHftRatio];
// TH1D* h1DcaZ1[nParticles][nEtasDca][nVzsDca][nCentDca][nPtBinsDca];
// TH1D* h1DcaXY1[nParticles][nEtasDca][nVzsDca][nCentDca][nPtBinsDca];
int const nCentDca = 9;
TH2D* h2Dca[nParticles][nEtasDca][nVzsDca][nCentDca][nPtBinsDca];
// TH2D* h2Dca[nParticles][nEtasDca][nVzsDca][nPhisDca][nPtBinsDca];

TH1D* hTpcPiPlus[nCentHftRatio];
TH1D* hTpcPiMinus[nCentHftRatio];
TH1D* hTpcKPlus[nCentHftRatio];
TH1D* hTpcKMinus[nCentHftRatio];

TString outFileName;// = "D0.toyMc.root";
std::pair<float, float> const momentumRange(0, 12);

float const gVzCut = 6.0e4;
float const acceptanceRapidity = 1.2;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCentHftRatio] = {31., 18.1, 12.8, 9.3, 7.2, 5.9, 5., 4.6, 4.};

int mMode = 0;
//const float M_D_0 = 1.86484;
const float M_B_0 = 5.27958;
const float M_B_PLUS = 5.27926;
//============== main  program ==================
void toyMcBtoD(int npart = 10, TString output = "D0.toyMc.root", TString particleName = "D0", bool isCombinB = false)
{
    // TMemStat mem;
    // mem.Enable();
    TStopwatch*   stopWatch = new TStopwatch();
    stopWatch->Start();
    gRandom->SetSeed();
    
    //In this part, if you want, add the random number to choose specific D0, B0, and B+ ratio
    //I plan to combine B0 and B+
    //

    if(particleName.CompareTo("D0")==0)      mMode = 1;
    else if(particleName.CompareTo("B0")==0) mMode = 2;
    else if(particleName.CompareTo("Bpm")==0) mMode = 3;
    else {
        cout << ">>>>>>>Please input the right particle name: D0, B0, Bpm ......<<<<<<<<" << endl;
        exit(1);
    }
    
    outFileName = output;
    if(!output.Contains(".root")) outFileName += ".root";
    bookObjects();
    
    pydecay = TPythia6Decayer::Instance();
    pydecay->Init();
    TPythia6::Instance()->SetMRPY(1, 88158204); //random seed number
    
    TLorentzVector* b_d = new TLorentzVector;
    TClonesArray ptl("TParticle", 100);
    for (int ipart = 0; ipart < npart; ipart++)
    {
        if (npart>1000 && !(ipart % (npart/100)))
            cout << "____________ ipart = " << ipart / static_cast<float>(npart) << " ________________" << endl;
        
        if(isCombinB) {
            float FR_B0 = 0.4;
            float FR_Bplus = 0.4;
            float FR_Sum = FR_B0 + FR_Bplus;
            float FR_Rdm = FR_Sum*gRandom->Rndm();
            if(FR_Rdm<FR_B0) mMode = 2;
            else if(FR_Rdm<FR_B0+FR_Bplus) mMode = 3;
            else continue;  //this doesn't work, because I close some decay branch
        }
        
        setDecayChannels(mMode);
        //see particle mass: http://www.johnmarcampbell.com/lambdaMixingDoc/da/dae/phys__constants_8h_source.html
        if(mMode==1) getKinematics(*b_d, M_D_0);
        else if(mMode==2) getKinematics(*b_d, M_B_0);
        else if(mMode==3) getKinematics(*b_d, M_B_PLUS);
        else {
            cout << ">>>>>>>>>>>No exact input particle<<<<<<<<<<" << endl;
            exit(1);
        }
        
        int KF;
        TVector3 v0_D0(0.,0.,0.);
        if(mMode==1) {
            KF = 421;
            if(gRandom->Rndm()<0.5) KF = -KF;
            decayAndFill(KF, b_d, v0_D0, fWeightFunction->Eval(b_d->Perp()), ptl);
        }
        else {
            if(mMode==2) KF = 511;
            if(mMode==3) KF = 521;
            bool ifGetD0 = GetD0FromBDecay(KF, *b_d, v0_D0, ptl);  //B decay to D0
            if(!ifGetD0) continue;
            decayAndFill(KF, b_d, v0_D0, fWeightFunction->Eval(b_d->Perp()), ptl);
        }
        if (ipart % 6000 == 1) nt->AutoSave("SaveSelf");
    }
    
    write();
    // mem.Show();
    stopWatch->Stop();
    stopWatch->Print();
}

void setDecayChannels(int const mMode)
{
    // Set some stable particle
    TPythia6::Instance()->SetMDCY(122,1,0); //Set D+/D- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(128,1,0); //Set Ds+/Ds- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(102,1,0); //Set pi0 to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(15,1,0); //Set tau+/tau- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(107,1,0); //Set rho+/rho- to be stable particle, not decay
    TPythia6::Instance()->SetMDCY(281,1,0); //Set a_1+/a_1- to be stable particle, not decay
    
    // Set intermediate particle channels -->D0
    // ---not set this, in order to keep the B->D0 branch ratio the same as Pythia
    
    // set D0 decay channel, constrained to D0->k-pi+, D0bar->k+pi-
    for (int idc = 747; idc < 807 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
    TPythia6::Instance()->SetMDME(763, 1, 1);
    
    // set B0 decay channels
    if(mMode==2) {
        for (int idc = 863; idc < 898 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
        for (int idc = 864; idc < 868 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
        for (int idc = 870; idc < 874 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
        TPythia6::Instance()->SetMDME(876, 1, 1);
        for (int idc = 880; idc < 882 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
        for (int idc = 885; idc < 886 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
    }
    
    // set B+ decay channels
    if(mMode==3) {
        for (int idc = 908; idc < 943 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
        for (int idc = 908; idc < 931 + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 1);
    }
}

bool GetD0FromBDecay(int& kf, TLorentzVector& b, TVector3& v0_D0, TClonesArray& daughters)
{
    if(gRandom->Rndm()<0.5) kf = -kf;
    pydecay->Decay(kf, &b);
    pydecay->ImportParticles(&daughters);
    
    TLorentzVector D0Mom;
    int nTrk = daughters.GetEntriesFast();
    //cout << "nTrk = " << nTrk << ", mMode = " << mMode << endl;
    bool ifGetD0 = false;
    for (int iTrk = 0; iTrk < nTrk; ++iTrk)
    {
        TParticle* ptl0_pre = (TParticle*)daughters.At(iTrk);
        //cout << "itrk = " << ptl0_pre->GetPdgCode() << endl;
        switch (ptl0_pre->GetPdgCode())
        {
            case 421:
                ptl0_pre->Momentum(D0Mom);
                v0_D0.SetXYZ(ptl0_pre->Vx() * 1000., ptl0_pre->Vy() * 1000., ptl0_pre->Vz() * 1000.); // converted to μm
                kf = 421;
                ifGetD0 = true;
                break;
            case -421:
                ptl0_pre->Momentum(D0Mom);
                v0_D0.SetXYZ(ptl0_pre->Vx() * 1000., ptl0_pre->Vy() * 1000., ptl0_pre->Vz() * 1000.); // converted to μm
                kf = -421;
                ifGetD0 = true;
                break;
            default:
                break;
        }
    }
    daughters.Clear();
    b = D0Mom;
    return ifGetD0;
}

void decayAndFill(int const kf, TLorentzVector* b, TVector3 const v0_D0, double const weight, TClonesArray& daughters)
{
    pydecay->Decay(kf, b);
    pydecay->ImportParticles(&daughters);
    
    TLorentzVector kMom;
    TLorentzVector pMom;
    TVector3 v00;
    
    int nTrk = daughters.GetEntriesFast();
    for (int iTrk = 0; iTrk < nTrk; ++iTrk)
    {
        TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
        
        switch (abs(ptl0->GetPdgCode()))
        {
            case 321:
                ptl0->Momentum(kMom);
                // v00.SetXYZ(0,0,0);
                v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
                break;
            case 211:
                ptl0->Momentum(pMom);
                break;
            default:
                break;
        }
    }
    daughters.Clear();
    
    TVector3 v0 = v00 + v0_D0;
    fill(kf, b, weight, kMom, pMom, v0);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 v00)
{
    int const centrality = floor(nCentHftRatio * gRandom->Rndm());
    
    TVector3 const vertex = getVertex(centrality);
    // smear primary vertex
    // float const sigmaVertex = sigmaVertexCent[cent];
    // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));
    
    v00 += vertex;
    
    // smear momentum
    //TLorentzVector const kRMom = smearMom(kMom, fKaonMomResolution);
    //TLorentzVector const pRMom = smearMom(pMom, fPionMomResolution);
    TLorentzVector kRMom;
    TLorentzVector pRMom;
    if(kf<0) {
        kRMom = smearMom(kMom, fKaonPlusMomResolution);
        pRMom = smearMom(pMom, fPionMinusMomResolution);
    }
    else {
        kRMom = smearMom(kMom, fKaonMinusMomResolution);
        pRMom = smearMom(pMom, fPionPlusMomResolution);
    }
    
    
    // smear position
    TVector3 const kRPos = smearPosData(1, vertex.z(), centrality, kRMom, v00);
    TVector3 const pRPos = smearPosData(0, vertex.z(), centrality, pRMom, v00);
    // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
    // TVector3 const pRPos = smearPos(pMom, pRMom, v00);
    
    // reconstruct
    TLorentzVector const rMom = kRMom + pRMom;
    float const kDca = dca(kMom.Vect(), v00, vertex);
    float const pDca = dca(pMom.Vect(), v00, vertex);
    float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
    float const kRSDca = dcaSigned(kRMom.Vect(), kRPos, vertex);
    float const kRDcaXY = dcaXY(kRMom.Vect(), kRPos, vertex);
    float const kRDcaZ = dcaZ(kRMom.Vect(), kRPos, vertex);
    float const pRDca = dca(pRMom.Vect(), pRPos, vertex);
    float const pRSDca = dcaSigned(pRMom.Vect(), pRPos, vertex);
    float const pRDcaXY = dcaXY(pRMom.Vect(), pRPos, vertex);
    float const pRDcaZ = dcaZ(pRMom.Vect(), pRPos, vertex);

//TLorentzVector const trMom = kMom + pMom;   
//cout << trMom.M() << endl;
 
    TVector3 v0;
    float const dca12 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, v0);
    float const decayLength = (v0 - vertex).Mag();
    float const dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
    float const D0DcaXY = dcaXY(rMom.Vect(), v0, vertex);
    float const D0DcaZ = dcaZ(rMom.Vect(), v0, vertex);
    float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
    float const angle12 = kRMom.Vect().Angle(pRMom.Vect());
    
    TLorentzVector kRMomRest = kRMom;
    TVector3 beta;
    beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
    kRMomRest.Boost(-beta);
    float const cosThetaStar = rMom.Vect().Unit().Dot(kRMomRest.Vect().Unit());
    
    int const charge = kf > 0 ? 1 : -1;
    // save
    float arr[110];
    int iArr = 0;
    arr[iArr++] = centrality;
    arr[iArr++] = vertex.X();
    arr[iArr++] = vertex.Y();
    arr[iArr++] = vertex.Z();
    arr[iArr++] = getVzIndexDca(vertex.Z());
    
    arr[iArr++] = kf;
    arr[iArr++] = weight;
    arr[iArr++] = b->M();
    arr[iArr++] = b->Perp();
    arr[iArr++] = b->PseudoRapidity();
    arr[iArr++] = b->Rapidity();
    arr[iArr++] = b->Phi();
    arr[iArr++] = v00.X();
    arr[iArr++] = v00.Y();
    arr[iArr++] = v00.Z();
    
    arr[iArr++] = rMom.M();
    arr[iArr++] = rMom.Perp();
    arr[iArr++] = rMom.PseudoRapidity();
    arr[iArr++] = rMom.Rapidity();
    arr[iArr++] = rMom.Phi();
    arr[iArr++] = v0.X();
    arr[iArr++] = v0.Y();
    arr[iArr++] = v0.Z();
    arr[iArr++] = reconstructD0(centrality, rMom);
    
    arr[iArr++] = dca12;
    arr[iArr++] = decayLength;
    arr[iArr++] = dcaD0ToPv;
    arr[iArr++] = D0DcaXY;
    arr[iArr++] = D0DcaZ;
    arr[iArr++] = cosTheta;
    arr[iArr++] = angle12;
    arr[iArr++] = cosThetaStar;
    
    arr[iArr++] = kMom.M();
    arr[iArr++] = kMom.Perp();
    arr[iArr++] = kMom.PseudoRapidity();
    arr[iArr++] = kMom.Rapidity();
    arr[iArr++] = kMom.Phi();
    arr[iArr++] = kDca;
    
    arr[iArr++] = kRMom.M();
    arr[iArr++] = kRMom.Perp();
    arr[iArr++] = kRMom.PseudoRapidity();
    arr[iArr++] = kRMom.Rapidity();
    arr[iArr++] = kRMom.Phi();
    arr[iArr++] = kRPos.X();
    arr[iArr++] = kRPos.Y();
    arr[iArr++] = kRPos.Z();
    arr[iArr++] = kRDca;
    arr[iArr++] = kRSDca;
    arr[iArr++] = kRDcaXY;
    arr[iArr++] = kRDcaZ;
    arr[iArr++] = getEtaIndexDca(kRMom.PseudoRapidity());
    arr[iArr++] = getPtIndexDca(kRMom.Perp());
    arr[iArr++] = tpcReconstructed(1, -1 * charge, centrality, kRMom);
    
    arr[iArr++] = pMom.M();
    arr[iArr++] = pMom.Perp();
    arr[iArr++] = pMom.PseudoRapidity();
    arr[iArr++] = pMom.Rapidity();
    arr[iArr++] = pMom.Phi();
    arr[iArr++] = pDca;
    
    arr[iArr++] = pRMom.M();
    arr[iArr++] = pRMom.Perp();
    arr[iArr++] = pRMom.PseudoRapidity();
    arr[iArr++] = pRMom.Rapidity();
    arr[iArr++] = pRMom.Phi();
    arr[iArr++] = pRPos.X();
    arr[iArr++] = pRPos.Y();
    arr[iArr++] = pRPos.Z();
    arr[iArr++] = pRDca;
    arr[iArr++] = pRSDca;
    arr[iArr++] = pRDcaXY;
    arr[iArr++] = pRDcaZ;
    arr[iArr++] = getEtaIndexDca(pRMom.PseudoRapidity());
    arr[iArr++] = getPtIndexDca(pRMom.Perp());
    arr[iArr++] = tpcReconstructed(0, charge, centrality, pRMom);
    
    arr[iArr++] = matchHft(1, vertex.z(), centrality, kRMom);
    arr[iArr++] = matchHft(0, vertex.z(), centrality, pRMom);

    arr[iArr++] = primaryPt;
    arr[iArr++] = primaryY;
    arr[iArr++] = primaryPhi;
    
    nt->Fill(arr);
}

void getKinematics(TLorentzVector& b, double const mass)
{
    float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
    float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
    float const phi = TMath::TwoPi() * gRandom->Rndm();

    primaryPt = pt;
    primaryY = y;
    primaryPhi = phi;      

    float const mT = sqrt(mass * mass + pt * pt);
    float const pz = mT * sinh(y);
    float const E = mT * cosh(y);
    
    b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
    
    return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 newPos(pos);
    newPos.SetZ(0);
    
    TVector3 newP(p);
    newP.SetZ(0);
    
    TVector3 newVertex(vertex);
    newVertex.SetZ(0);
    
    TVector3 posDiff = newPos - newVertex;
    float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
    return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
    TVector3 posDiff = pos - vertex;
    if (sin(p.Theta()) == 0) return 0;
    else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0)
{
    TVector3 posDiff = pos2 - pos1;
    TVector3 pu1 = p1.Unit();
    TVector3 pu2 = p2.Unit();
    double pu1Pu2 = pu1.Dot(pu2);
    double g = posDiff.Dot(pu1);
    double k = posDiff.Dot(pu2);
    double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
    double s1 = g + s2 * pu1Pu2;
    TVector3 posDca1 = pos1 + pu1 * s1;
    TVector3 posDca2 = pos2 + pu2 * s2;
    v0 = 0.5 * (posDca1 + posDca2);
    return (posDca1 - posDca2).Mag();
}

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
    float const pt = b.Perp();
    float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));
    
    TLorentzVector sMom;
    sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
    return sMom;
}

TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos)
{
    float thetaMCS = 13.6 / mom.Beta() / rMom.P() / 1000 * sqrt(pxlLayer1Thickness / fabs(sin(mom.Theta())));
    float sigmaMCS = thetaMCS * 28000 / fabs(sin(mom.Theta()));
    float sigmaPos = sqrt(pow(sigmaMCS, 2) + pow(sigmaPos0, 2));
    
    return TVector3(gRandom->Gaus(pos.X(), sigmaPos), gRandom->Gaus(pos.Y(), sigmaPos), gRandom->Gaus(pos.Z(), sigmaPos));
}

int getPtIndexDca(double pT)
{
    for (int i = 0; i < nPtBinsDca; i++)
    {
        if ((pT >= ptEdgeDca[i]) && (pT < ptEdgeDca[i + 1]))
            return i;
    }
    return nPtBinsDca - 1 ;
}

int getEtaIndexDca(double Eta)
{
    for (int i = 0; i < nEtasDca; i++)
    {
        if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
            return i;
    }
    return nEtasDca - 1 ;
}

int getVzIndexDca(double Vz)
{
    for (int i = 0; i < nVzsDca; i++)
    {
        if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
            return i;
    }
    return nVzsDca - 1 ;
}

int getPhiIndexDca(double Phi)
{
    for (int i = 0; i < nPhisDca; i++)
    {
        if ((Phi >= PhiEdgeDca[i]) && (Phi < PhiEdgeDca[i + 1]))
            return i;
    }
    return nPhisDca - 1 ;
}

int getPtIndexHftRatio(double pT)
{
    for (int i = 0; i < nPtBinsHftRatio; i++)
    {
        if ((pT >= ptEdgeHftRatio[i]) && (pT < ptEdgeHftRatio[i + 1]))
            return i;
    }
    return nPtBinsHftRatio - 1 ;
}

int getEtaIndexHftRatio(double Eta)
{
    for (int i = 0; i < nEtasHftRatio; i++)
    {
        if ((Eta >= EtaEdgeHftRatio[i]) && (Eta < EtaEdgeHftRatio[i + 1]))
            return i;
    }
    return nEtasHftRatio - 1 ;
}

int getVzIndexHftRatio(double Vz)
{
    for (int i = 0; i < nVzsHftRatio; i++)
    {
        if ((Vz >= VzEdgeHftRatio[i]) && (Vz < VzEdgeHftRatio[i + 1]))
            return i;
    }
    return nVzsHftRatio - 1 ;
}

int getPhiIndexHftRatio(double Phi)
{
    for (int i = 0; i < nPhisHftRatio; i++)
    {
        if ((Phi >= PhiEdgeHftRatio[i]) && (Phi < PhiEdgeHftRatio[i + 1]))
            return i;
    }
    return nPhisHftRatio - 1 ;
}

TVector3 smearPosData(int const iParticleIndex, double const vz, int cent, TLorentzVector const& rMom, TVector3 const& pos)
{
    int const iEtaIndex = getEtaIndexDca(rMom.PseudoRapidity());
    int const iVzIndex = getVzIndexDca(vz);
    // int const iPhiIndex = getPhiIndexDca(rMom.Phi());
    int const iPtIndex = getPtIndexDca(rMom.Perp());
    
    double sigmaPosZ = 0;
    double sigmaPosXY = 0;
    
    // if (cent == 8) cent = 7;
    // All the centrality position smear was based on 0-10% centrality input
    // changed to 0-80%
    
    h2Dca[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom2(sigmaPosXY,sigmaPosZ);
    // h2Dca[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][iPtIndex]->GetRandom2(sigmaPosXY, sigmaPosZ);
    sigmaPosZ *= 1.e4;
    sigmaPosXY *= 1.e4;
    /*if (h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
     {
     do sigmaPosZ = h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
     while (fabs(sigmaPosZ) > 1.e3);
     }
     
     if (h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->ComputeIntegral())
     {
     do sigmaPosXY = h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;
     while (fabs(sigmaPosXY) > 1.e3);
     }
     */
    
    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
    newPos -= momPerp.Unit() * sigmaPosXY;
    
    return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TVector3 getVertex(int const centrality)
{
    double rdmVz;
    
    if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
    else
    {
        do rdmVz = h1Vz[centrality]->GetRandom() * 1e4;
        while (fabs(rdmVz) > gVzCut);
    }
    
    return TVector3(0., 0., rdmVz);
}

bool reconstructD0(int const centrality, TLorentzVector const& mom)
{
    TGraph* gr = NULL;
    
    if (centrality < 4) gr = grEff[0];
    else if (centrality < 7) gr = grEff[1];
    else gr = grEff[2];
    
    return gRandom->Rndm() < gr->Eval(mom.Perp());
}

bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom)
{
    TH1D* h = NULL;
    
    if (iParticleIndex == 0)
    {
        if (charge > 0) h = hTpcPiPlus[cent];
        else h = hTpcPiMinus[cent];
    }
    else
    {
        if (charge > 0) h = hTpcKPlus[cent];
        else h = hTpcKMinus[cent];
    }
    
    int const bin = h->FindBin(mom.Perp());
    
    return gRandom->Rndm() < h->GetBinContent(bin);
}

float matchCorrection(int iParticleIndex, int iEtaIndex, int iVzIndex, int iPhiIndex,  int iCent,  float pt) // correct HFT ratios difference for inclusive and primary tracks
{
   int const bin = hHftRatioCorrect[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][iCent]->FindBin(pt);
   return hHftRatioCorrect[iParticleIndex][iEtaIndex][iVzIndex][iParticleIndex][iCent]->GetBinContent(bin);
}

bool matchHft(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& mom)
{
    int const iEtaIndex = getEtaIndexHftRatio(mom.PseudoRapidity());
    int const iVzIndex = getVzIndexHftRatio(vz);
    int const iPhiIndex = getPhiIndexHftRatio(mom.Phi());

    float AdditionalCorrectFactor = matchCorrection(iParticleIndex, iEtaIndex, iVzIndex, iPhiIndex, 0, mom.Perp());//inclusive / primary
    if(fabs(AdditionalCorrectFactor - 0.) < 1.e-5) AdditionalCorrectFactor = 1.0;// if this correction factor is 0, then use 1 instead
    
    // int const bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->FindBin(mom.Perp());
    int bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->FindBin(mom.Perp());
    if(hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->GetBinContent(bin) < 0.005 && mom.Perp()>4.0) bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->FindBin(4.0);// use 4 GeV as HFT ratio for high pT when there is no statistics
    // return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->GetBinContent(bin);
    if(cent == 0 && mom.Perp()>4.0) //for 70-80%, pt>4GeV. HFT matching ratio is limited by statistics, updated with 60-70%
    {
      return gRandom->Rndm() < (hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent+1]->GetBinContent(bin)/AdditionalCorrectFactor);
    }
    else
    {
      return gRandom->Rndm() < (hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->GetBinContent(bin)/AdditionalCorrectFactor);
    }
}
//___________
void bookObjects()
{
    cout << "Loading input momentum resolution ..." << endl;
    //TFile f("momentum_resolution.root");
    // TFile f("Momentum_resolution_May11.root");
    TFile f("Momentum_resolution_FromHijing.root");// new Momentun resolution 
    //fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPion");
    //fKaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaon");
    // fPionPlusMomResolution = (TF1*)f.Get("fPionPlus")->Clone("fPionPlus");
    // fPionMinusMomResolution = (TF1*)f.Get("fPionMinus")->Clone("fPionMinus");
    // fKaonPlusMomResolution = (TF1*)f.Get("fKaonPlus")->Clone("fKaonPlus");
    // fKaonMinusMomResolution = (TF1*)f.Get("fKaonMinus")->Clone("fKaonMinus");
    fPionPlusMomResolution = (TF1*)f.Get("fPionPlus_3")->Clone("fPionPlus");
    fPionMinusMomResolution = (TF1*)f.Get("fPionPlus_3")->Clone("fPionMinus");
    fKaonPlusMomResolution = (TF1*)f.Get("fKaonMinus_3")->Clone("fKaonPlus");
    fKaonMinusMomResolution = (TF1*)f.Get("fKaonMinus_3")->Clone("fKaonMinus");
    f.Close();
    
    cout << "Loading input spectra ..." << endl;
    TFile fPP("pp200_spectra.root");
    fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("f1Levy");
    fPP.Close();
    
    // TFile fVertex("Vz_Cent_May25.root");
    TFile fVertex("Vz_Cent.root");// new Vz
    
    for (int ii = 0; ii < nCentHftRatio; ++ii)
    {
        h1Vz[ii]      = (TH1D*)(fVertex.Get(Form("mh1Vz_%i", ii)));
        h1Vz[ii]->SetDirectory(0);
    }
    
    fVertex.Close();
    
    std::cout << "Loading HFT ratios Correction factor for inclusive and primary..." << std::endl;
    TFile fHftRatioCorrect("HFT_Ratio_Correction_Hijing_QM.root");
    for (int iParticle = 0; iParticle < nParticles; ++iParticle)
    {
       for (int iCent = 0; iCent < 1 ; ++iCent)
       {
          // HFT ratio
          for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
          {
             for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
             {
                for (int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
                {
                   hHftRatioCorrect[iParticle][iEta][iVz][iPhi][iCent]  = (TH1D*)(fHftRatioCorrect.Get(Form("mhHFTRatio_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent)));
                   hHftRatioCorrect[iParticle][iEta][iVz][iPhi][iCent] ->SetDirectory(0);
                }
             }
          }
       }
    }
    fHftRatioCorrect.Close();

    cout << "Loading input HFT ratios and DCA ..." << endl;
    // TFile fHftRatio1("HFT_Ratio_VsPt_Centrality_Eta_Phi_Vz_Zdcx_170813_combineCharge.root");
    TFile fHftRatio1("HFT_Ratio_VsPt_Centrality_Eta_Phi_Vz_Zdcx_170813_combineCharge_reBin.root");
    TFile fDca1("2DProjection_simCent_NoBinWidth_3D_Dca_VsPt_Centrality_Eta_Phi_Vz_Zdcx_20170813_combineCharge.root");
    
    for (int iParticle = 0; iParticle < nParticles; ++iParticle)
    {
        for (int iCent = 0; iCent < nCentHftRatio; ++iCent)
        {
            // HFT ratio
            for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
            {
                for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
                {
                    for (int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
                    {
                        hHftRatio1[iParticle][iEta][iVz][iPhi][iCent] = (TH1D*)(fHftRatio1.Get(Form("mh1HFT1PtCentPartEtaVzPhiRatio_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent)));
                        hHftRatio1[iParticle][iEta][iVz][iPhi][iCent]->SetDirectory(0);
                    }
                }
            }
        }
        cout << "Finished loading HFT Ratio: " <<  endl;
        
        for(int iCent = 0; iCent < nCentDca; ++iCent)
        {
            // DCA
            for (int iEta = 0; iEta < nEtasDca; ++iEta)
            {
                for (int iVz = 0; iVz < nVzsDca; ++iVz)
                {
                    // for (int iPhi = 0; iPhi < nPhisDca; ++iPhi)
                    for (int iPt = 0; iPt < nPtBinsDca; ++iPt)
                    {
                        h2Dca[iParticle][iEta][iVz][iCent][iPt] = (TH2D*)((fDca1.Get(Form("mh2DcaPtCentPartEtaVzPhi_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iCent, iPt))));
                        h2Dca[iParticle][iEta][iVz][iCent][iPt]->SetDirectory(0);
                    }
                }
            }
        }
        // cout << "Finished loading centrality: " << iCent << endl;
    }
    cout << "Finished loading Dca: " <<  endl;
    
    fHftRatio1.Close();
    fDca1.Close();
    
    cout << " Loading TPC tracking efficiencies " << endl;
    
    TFile fTpcPiPlus("Eff_PionPlus_embedding.root");
    TFile fTpcPiMinus("Eff_PionMinus_embedding.root");
    TFile fTpcKPlus("Eff_KaonPlus_embedding.root");
    TFile fTpcKMinus("Eff_KaonMinus_embedding.root");
    
    for (int iCent = 0; iCent < nCentHftRatio; ++iCent)
    {
        hTpcPiPlus[iCent] = (TH1D*)fTpcPiPlus.Get(Form("h1Ratiocent_%i", iCent));
        hTpcPiPlus[iCent]->SetDirectory(0);
        hTpcPiMinus[iCent] = (TH1D*)fTpcPiMinus.Get(Form("h1Ratiocent_%i", iCent));
        hTpcPiMinus[iCent] ->SetDirectory(0);
        hTpcKPlus[iCent] = (TH1D*)fTpcKPlus.Get(Form("h1Ratiocent_%i", iCent));
        hTpcKPlus[iCent]->SetDirectory(0);
        hTpcKMinus[iCent] = (TH1D*)fTpcKMinus.Get(Form("h1Ratiocent_%i", iCent));
        hTpcKMinus[iCent]->SetDirectory(0);
    }
    
    fTpcPiPlus.Close();
    fTpcPiMinus.Close();
    fTpcKPlus.Close();
    fTpcKMinus.Close();
    
    cout << "Done with loading all files ..." << endl;
    
    grEff[0] = new TGraph("eff_4080.csv", "%lg %lg", ",");
    grEff[1] = new TGraph("eff_1040.csv", "%lg %lg", ",");
    grEff[2] = new TGraph("eff_010.csv", "%lg %lg", ",");
    
    result = new TFile(outFileName.Data(), "recreate");
    result->SetCompressionLevel(1);
    result->cd();
    
    int BufSize = (int)pow(2., 16.);
    // int Split = 1;
    nt = new TNtuple("nt", "", "cent:vx:vy:vz:vzIdx:"
                     "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC D0
                     "rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:reco:" // Rc D0
                     "dca12:decayLength:dcaD0ToPv:dcaXY:dcaZ:cosTheta:angle12:cosThetaStar:" // Rc pair
                     "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
                     "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:kRSDca:kRDcaXY:kRDcaZ:kEtaIdx:kPtIdx:kTpc:" // Rc Kaon
                     "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion1
                     "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:pRSDca:pRDcaXY:pRDcaZ:pEtaIdx:pPtIdx:pTpc:" // Rc Pion1
                     "kHft:pHft:primaryPt:primaryY:primaryPhi", BufSize); //add primary pt, y, phi
    // nt->SetAutoSave(-500000); // autosave every 1 Mbytes
}
//___________
void write()
{
    result->cd();
    nt->Write();
    result->Close();
}
