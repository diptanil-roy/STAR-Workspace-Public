// R__LOAD_LIBRARY(/opt/homebrew/opt/fastjet/lib/libfastjet.dylib);
using namespace std;

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TProfile.h"
// #include "TPythia8.h"
#include "TParticle.h"
// #include "TDatabasePDG.h"
#include <TLorentzVector.h>
#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "THn.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TLegend.h"
#include <vector>


#endif

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/Selector.hh>

// Defining a User Class to store PID info for particles used as input to FastJet
class MyParticleId : public fastjet::PseudoJet::UserInfoBase{
    public:
        MyParticleId(const int & pdg_id_in, const int & part_id_in) : _pdg_id(pdg_id_in), _part_id(part_id_in){}

        int pid() const{
            return _pdg_id;
        }
        int partn() const{
            return _part_id;
        }
    private:
        int _pdg_id;
        int _part_id;
};

void BinLogX(TH1*h)
{
	TAxis *axis = h->GetXaxis();
	int bins = axis->GetNbins();
	Axis_t from = axis->GetXmin();
	Axis_t to = axis->GetXmax();
	Axis_t width = (to - from) / bins;
	Axis_t *new_bins = new Axis_t[bins + 1];

	for (int i = 0; i <= bins; i++) {
		new_bins[i] = TMath::Power(10, from + i * width);
		cout << new_bins[i] << ", ";
	}
	cout << endl;
	axis->Set(bins, new_bins);
	delete[] new_bins;
}

void BinLogY(TProfile *h)
{
	TAxis *axis = h->GetYaxis();
	int bins = axis->GetNbins();
	Axis_t from = axis->GetXmin();
	Axis_t to = axis->GetXmax();
	Axis_t width = (to - from) / bins;
	Axis_t *new_bins = new Axis_t[bins + 1];

	for (int i = 0; i <= bins; i++) {
		new_bins[i] = TMath::Power(10, from + i * width);
		// cout << new_bins[i] << ", ";
	}
	// cout << endl;
	axis->Set(bins, new_bins);
	delete[] new_bins;
}

Double_t dPhi(Double_t phi1, Double_t phi2) {
  Double_t deltaPhi;
  deltaPhi = abs(phi1 - phi2); //TODO absolute values
  if (deltaPhi>(2*TMath::Pi()))  deltaPhi-=2*(TMath::Pi());
  if (deltaPhi<(0*TMath::Pi())) deltaPhi+=2*(TMath::Pi()); 

  if (deltaPhi > TMath::Pi()) deltaPhi= 2*(TMath::Pi()) - deltaPhi;
  return deltaPhi;   // dphi in [0, 2Pi]
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dEta(Double_t eta1, Double_t eta2) {
  Double_t deltaEta;
  deltaEta = abs(eta1 - eta2);

  return deltaEta;
}

// function to calculate relative eta between 2 objects
//___________________________________________________________________________________________
Double_t dR(Double_t delphi, Double_t deleta) {
  Double_t dRad;
  dRad = TMath::Sqrt(pow(delphi,2) + pow(deleta,2));

  return dRad;
}


int main(){

	TString BaseDir = "./";
	// if (!strcmp(argv[3], "")) BaseDir = argv[2];

	cout << BaseDir.Data() << endl;
	
	TH1::SetDefaultSumw2();
  	TH2::SetDefaultSumw2();
	TProfile::SetDefaultSumw2();
  	// gStyle->SetOptStat(0);

  	TString InputFileName;
  	InputFileName = Form("%s/PYTHIA.root", BaseDir.Data());

	if (InputFileName == "") return 0;

	// Import Tree

	TFile *f = new TFile(InputFileName.Data());

	cout << f->GetName() <<  endl;
	f->cd();

	TTree *MuTree = (TTree *)gDirectory->Get("MuTree");
	// TTree *RecoJetTree = (TTree *)gDirectory->Get("RecoJets");

	int numberofbinary;
    float impactparameter;

    int nparticles;
    float mPt[5000];
    float mEta[5000];
    float mPhi[5000];
    float mEnergy[5000];
    int mPID[5000];

  	MuTree->SetBranchAddress("NParticles", &nparticles);
  	MuTree->SetBranchAddress("TrackPt", mPt);
  	MuTree->SetBranchAddress("TrackEta", mEta);
  	MuTree->SetBranchAddress("TrackPhi", mPhi);
  	MuTree->SetBranchAddress("TrackEnergy", mEnergy);
  	MuTree->SetBranchAddress("TrackPID", mPID);

  	int nentries = MuTree->GetEntriesFast();

  	// Define Fastjet

  	fastjet::Strategy strategy = fastjet::Best;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; //Change as you need
    fastjet::JetDefinition *jetdefinition = NULL;
    fastjet::JetAlgorithm algorithm = fastjet::antikt_algorithm; // Will use this one for now
    double R = 0.4;
    // Defining the Jet
    jetdefinition = new fastjet::JetDefinition(algorithm, R, recombScheme, strategy);
    // Defining Vectors for FastJet inputs
    std::vector <fastjet::PseudoJet> fjInputs; // Will store px, py, pz, E info

    // Defining histograms here
    // 4 histograms for 4 kinds of particles (charged, pi, K, p)
    // 2 histograms for unlike and like sign

    const int njpt_gen_bins_var = 6; //x2 MC
    double jetpt_var_bin[njpt_gen_bins_var+1] = {5,7,9,11,13,15,20};

    const int nz_gen_bins = 7; //x2 MC
    double z_gen_bin[nz_gen_bins+1] = {0., 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.};

    const int ndrbins = 5;
    double drbins[ndrbins+1] = {0.0, 0.05, 0.1, 0.2, 0.4, 0.6};

    TH1D *D0Pt;
    TH1D *LeadJetPt;
    TH1D *LeadJetZ;
    TH1D *LeadJetdR;

    D0Pt = new TH1D("D0Pt", "D0Pt", 20, 0, 10);
    LeadJetPt = new TH1D("LeadJetPt", "LeadJetPt", njpt_gen_bins_var, jetpt_var_bin);
    LeadJetZ = new TH1D("LeadJetZ", "LeadJetZ", nz_gen_bins, z_gen_bin);
    LeadJetdR = new TH1D("LeadJetdR", "LeadJetdR", ndrbins, drbins);


  	for (int i = 0; i < nentries; i++){
  		MuTree->GetEntry(i);
  		// cout << numberofbinary << "\t" << impactparameter << endl;
  		//Reset FastJet input for every event
        fjInputs.clear();
        fjInputs.resize(0);

        for (int track = 0; track < nparticles; track++){
        	TVector3 p;
        	p.SetPtEtaPhi(mPt[track], mEta[track], mPhi[track]);

            if (mPt[track] < 0.2 || mPt[track] > 30.0) continue;
            if (abs(mEta[track]) > 1.0) continue;
			
        	fastjet::PseudoJet k(p.X(), p.Y(), p.Z(), mEnergy[track]);
        	k.set_user_info(new MyParticleId(mPID[track], track));
            if (abs(mPID[track]) == 421 && mPt[track] > 1. && mPt[track] < 10. ) D0Pt->Fill(mPt[track]);
	        // Store as input to Fastjet
	        fjInputs.push_back(k);
        }

        fjInputs = sorted_by_pt(fjInputs);

        if (fjInputs.size() == 0) continue;
        // if (fjInputs[0].pt() < 2.0) continue;

        // Cooking Jets Now (Run FastJet Algorithm)
	    vector <fastjet::PseudoJet> inclusiveJets, sortedJets, jets;
	    fastjet::ClusterSequence clustSeq(fjInputs, *jetdefinition);

	    // Extract Inclusive Jets sorted by pT
	    inclusiveJets = clustSeq.inclusive_jets(5.0);
	    sortedJets = sorted_by_pt(inclusiveJets);

	    fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-0.6, 0.6); // Selects |eta| < 0.6
	    jets = eta_selector(sortedJets);

        if (jets.size() == 0) continue;

		for (int jet = 0; jet < jets.size(); jet++){
			if (jets[0].pt() > 20.0) continue; // If Lead Jet is < 10 GeV, we do not need to consider the event.

			vector <fastjet::PseudoJet> constituents = sorted_by_pt(jets[0].constituents());
			
			bool D0Jet = kFALSE;
            int D0Index = -1;

			for (int con = 0; con < constituents.size(); con++){
                if (abs(constituents[con].user_info<MyParticleId>().pid()) == 421) {
                    if (constituents[con].pt() > 1 && constituents[con].pt() < 10.0 && abs(constituents[con].eta()) < 1.0) {
                        D0Jet = kTRUE;
                        D0Index = con;
                        break;
                    }
                }
            }

			if (D0Jet) {
                LeadJetPt->Fill(jets[jet].pt());
                double z = (constituents[D0Index].px()*jets[jet].px() + constituents[D0Index].py()*jets[jet].py())/(jets[jet].px()*jets[jet].px() + jets[jet].py()*jets[jet].py());
                if (z >= 1.0) z = 0.9999;
                LeadJetZ->Fill(z);
                double dR = sqrt(pow(dPhi(constituents[D0Index].phi(), jets[jet].phi()), 2) + pow(dEta(constituents[D0Index].eta(), jets[jet].eta()), 2));
                LeadJetdR->Fill(dR);
            }

		}
    }

    TFile *out = new TFile(Form("%s/D0Jet_Out.root", BaseDir.Data()), "RECREATE");
    out->cd();
    D0Pt->Write();
    LeadJetPt->Write();
    LeadJetZ->Write();
    LeadJetdR->Write();
    out->Close();

	return 1;
}

