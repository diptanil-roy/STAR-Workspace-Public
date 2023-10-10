#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include <TLorentzVector.h>
#ifndef __CINT__
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
// #include "Pythia8/Pythia.h"
using namespace std;
#endif


using namespace Pythia8;

class MyUserHooks : public UserHooks {

  public:

    MyUserHooks(){}

    virtual bool canModifySigma() {return true;}

    virtual bool canVetoAfterHadronization() {return true;}

    virtual bool doVetoAfterHadronization(const Event& event) {
      bool go = false;

      // cout << event.size() << endl;
      for (int i = 0; i < event.size(); i++){
        const Particle& p = event[i];
        // if (!p.isFinal()) continue;
        // if (p.m()!=1.865) continue;
        cout << p.m() << endl;
        if (p.pT()<1.) continue;

        cout << "Found a D0" << endl;
        go = true;
        if (go) break;
      }

      return go;
    }
  private:

};

void PythiaCrossSections(Int_t nev = 1000000, const int random_seed = 181349, const char* VertexFileName = "", float R = 0.4, Int_t power = -1, Bool_t printjetinfo = kFALSE, Bool_t smallOutput = kTRUE){

  //nev = Number of Events (Equal to the number of events asked for, or the number of events in the vertex file, whichever is less.)
  //Random Seed = set by hand to be able to regenerate the dataset if needed.
  //Vertexfilename = Text files containing the vertex information from picodsts
  //R = Radius of Jets
  //Power = Defines the jet reconstruction algorithm, -1 stands for anti-kT and doesn't need to be changed
  //printjetinfo = For debugging



  TStopwatch timer;

  timer.Start();

  // Load Pythia 8
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  gRandom = new TRandom3();

  // pT HARD Bins defined here
  double pTlimit[2] = {0, -1}; 
  // double pTlimit[2] = {3, 0}; 
  // double pTlimit[2] = {0, 3}; 

  // FileName for the .root Pythia Events Files, salvaged from the vertex files names to ensure correspondence
  TString filename;

  filename = "PYTHIA";
  filename += "_" + std::to_string(int(pTlimit[0])) + "_" + std::to_string(int(pTlimit[1]));

  filename += ".root";

  // Defining a file to save Jets
  TFile *jetfile = new TFile(filename.Data(), "RECREATE");

  // Create pythia8 object
  TPythia8        *pythia = new TPythia8(); 
  Pythia8::Pythia *pythia8  = pythia->Pythia8();
  TDatabasePDG *pdg = new TDatabasePDG();

  // Configure Pythia to generate events
  // I will use all HardQCD processes, and make sure the D0s do not decay. Perform a check by invoking isFinal for PID 421

  pythia->ReadString("Random:setSeed = on");
  pythia->ReadString(Form("Random:seed = %i", random_seed));

  if (smallOutput) {
      pythia8->readString("Init:showProcesses = off");
      pythia8->readString("Init:showMultipartonInteractions = off");
      pythia8->readString("Init:showChangedSettings = off");
      pythia8->readString("Init:showChangedParticleData = off");
      pythia8->readString("Next:numberCount = 1000000000");
      pythia8->readString("Next:numberShowInfo = 0");
      pythia8->readString("Next:numberShowProcess = 0");
      pythia8->readString("Next:numberShowEvent = 0");
  }

  // pythia8->readString("HardQCD:gg2ccbar = on");
  // pythia8->readString("HardQCD:qqbar2ccbar = on");
  // pythia8->readString("Charmonium:all = on");
  // pythia8->readString("HardQCD:gg2bbbar = on");
  // pythia8->readString("HardQCD:qqbar2bbbar = on");
  // pythia8->readString("Charmonium:all = on");
  pythia8->readString("HardQCD:all = on");
  // pythia8->readString("SoftQCD:nonDiffractive = on");

  // pythia8->readString("PDF:pSet = 2");
  /// STAR TUNE
  pythia8->readString("MultipartonInteractions:bprofile=2");

  pythia8->readString("MultipartonInteractions:pT0Ref=1.4");
  pythia8->readString("MultipartonInteractions:ecmPow=0.135");
  pythia8->readString("MultipartonInteractions:coreRadius=0.56");
  pythia8->readString("MultipartonInteractions:coreFraction=0.78");
  pythia8->readString("ColourReconnection:range=5.4");

  // pythia8->readString("BeamRemnants:primordialKThard = 0.9");
  // pythia8->readString("BeamRemnants:primordialKTsoft = 0.45");

  // pythia8->readString("SigmaProcess:renormScale2=3");
  // pythia8->readString("SigmaProcess:factorScale2=3");
  // pythia8->readString("SigmaProcess:renormMultFac = 0.25");
  // pythia8->readString("SigmaProcess:factorMultFac = 0.25");

  // pythia8->readString("StringFlav:mesonCvector = 1.5");
  // pythia8->readString("StringFlav:mesonBvector = 3");

  // pythia8->readString("4:m0 = 1.43");
  // pythia8->readString("5:m0 = 4.30");

  /*
  D0 Decay Stopped, Will happen in GEANT
  */

  pythia8->particleData.readString("421:mayDecay = 0");

  // Turning off weak and long-lived particle decays! They will only be decayed by GEANT

  /*
  111:mayDecay = 0
  211:mayDecay = 0
  221:mayDecay = 0
  321:mayDecay = 0
  130:mayDecay = 0
  310:mayDecay = 0
  3122:mayDecay = 0
  3212:mayDecay = 0
  3112:mayDecay = 0
  3222:mayDecay = 0
  3312:mayDecay = 0
  3322:mayDecay = 0
  3334:mayDecay = 0
  */

  pythia8->particleData.readString("111:mayDecay = 0");
  pythia8->particleData.readString("211:mayDecay = 0");
  pythia8->particleData.readString("221:mayDecay = 0");
  pythia8->particleData.readString("321:mayDecay = 0");
  pythia8->particleData.readString("130:mayDecay = 0");
  pythia8->particleData.readString("310:mayDecay = 0");
  pythia8->particleData.readString("3122:mayDecay = 0");
  pythia8->particleData.readString("3212:mayDecay = 0");
  pythia8->particleData.readString("3112:mayDecay = 0");
  pythia8->particleData.readString("3222:mayDecay = 0");
  pythia8->particleData.readString("3312:mayDecay = 0");
  pythia8->particleData.readString("3322:mayDecay = 0");
  pythia8->particleData.readString("3334:mayDecay = 0");

  // auto myUserHooks = make_shared<MyUserHooks>();
  // pythia8->setUserHooksPtr( myUserHooks);

  // pythia8->particleData.listChanged();

  // Open the vertex file
  fstream vertexfile;
  vertexfile.open(VertexFileName, ios::in);

  // The same macro can be looped over multiple pt hard bins but for now, it is simpler to change the pt hard bins by hand, and have different macros for each bin
  // for (int ipT = 0; ipT < 1; ipT++){

    cout << "pT Bin : " << pTlimit[0] << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< " << pTlimit[1] << endl; 
    pythia8->settings.parm("PhaseSpace:pTHatMin", pTlimit[0]);
    pythia8->settings.parm("PhaseSpace:pTHatMax", pTlimit[1]);

    // Initialise Pythia
    // Beam 1 = p
    // Beam 2 = p
    // Ecm =  200 GeV/c^2
    pythia->Initialize(2212, 2212, 200.);
    
    cout << "Pythia initialised." << endl;

    // TClonesArray to access the list of particles generated by Pythia
    TClonesArray *particles = (TClonesArray*)pythia->GetListOfParticles();

    float pthat;
    float eventCrossSection;

    cout << "Asked for " << nev << " events" << endl;

    int d0Events = 0;

    bool makeevents = kTRUE;

    if (makeevents){

      TH1F *hD0Pt = new TH1F("hD0Pt", "hD0Pt", 20, 0.0, 10.0);

      // EVENT LOOP
      for (int event = 0; event < nev; event++){

        if (event % 10000 == 0) { cout << "Event # " << event << "\t" << endl;}
        pythia->GenerateEvent();
        //if (event < 1) pythia->EventListing();

        // pythia->PrintStatistics();

        //Import all particles generated to the TClonesArray defined earlier)
        pythia->ImportParticles(particles, "All");
        Int_t npart = particles->GetEntriesFast();

        double weight = pythia8->info.weight();
        // pythia8->info.list();
        pthat = pythia8->info.pTHat();
        eventCrossSection = pythia8->info.sigmaGen();

        if (printjetinfo) cout << "Weight and pTHat is = " << weight << "\t" << eventCrossSection << "\t" << pthat << endl;

        bool d0Event = kFALSE; // Boolean to select a D0 event

        for (int part1 = 0; part1 < npart; part1++){
          TParticle *MPart1 = (TParticle *) particles->At(part1);
          int pid1 = MPart1->GetPdgCode();
          Int_t ist = MPart1->GetStatusCode();

          // if (MPart1->Pt() < 0.2 || abs(MPart1->Eta()) > 1) continue; // No point in having D0s beyond our acceptance!

          TLorentzVector p;
          p.SetPxPyPzE(MPart1->Px(), MPart1->Py(), MPart1->Pz(), MPart1->Energy());
          if (abs(p.Rapidity()) > 1) continue;

          if (abs(pid1) == 421 && MPart1->Pt() > 1.0) {d0Event = kTRUE; break;}
        }

        if (!d0Event) continue;
        d0Events++;

        for (int part1 = 0; part1 < npart; part1++){
          TParticle *MPart1 = (TParticle *) particles->At(part1);
          int pid1 = MPart1->GetPdgCode();

          TLorentzVector p;
          p.SetPxPyPzE(MPart1->Px(), MPart1->Py(), MPart1->Pz(), MPart1->Energy());

          if (abs(pid1) == 421 && MPart1->Pt() >= 1.0 && abs(p.Rapidity()) < 1) hD0Pt->Fill(MPart1->Pt(), 1./MPart1->Pt());
        }

  	  }

      cout << nev << " events processed, " << d0Events << " D0 Events found." << endl;

      pythia8->stat();
    	cout << "Cross-Section = " << pythia8->info.sigmaGen() << "\t" << pythia8->info.sigmaErr() << endl;


      TFile *f = new TFile("tmp.root", "RECREATE");
      f->cd();
      hD0Pt->Write();
      f->Close();
    }

    TFile *g = new TFile("tmp.root");
    // TFile *g = new TFile("D0Pt.root");
    g->cd();

    TH1F *h = (TH1F *)gDirectory->Get("hD0Pt");
    

    TH1F *FONLL = new TH1F("FONLL","FONLL",60,0,30);
    TH1F *MinFONLL = new TH1F("MIN FONLL","MIN FONLL",60,0,30);
    TH1F *MaxFONLL = new TH1F("MAX FONLL","MAX FONLL",60,0,30);
    ifstream data1("FONLL_New.txt");
    int cnt=1;
    if(data1.is_open()){
      while(!data1.eof()){
        double x;double y;double minerr; double maxerr; double tmp;
        data1 >> x >> y >> minerr >> maxerr >> tmp >> tmp >> tmp >> tmp;
        FONLL->SetBinContent(cnt,y*pow(10, -9)/x);
        MinFONLL->SetBinContent(cnt,minerr*pow(10, -9)/x);
        MaxFONLL->SetBinContent(cnt,maxerr*pow(10, -9)/x);

        // FONLL->SetBinContent(cnt,y/(7.5330e+07));
        // MinFONLL->SetBinContent(cnt,minerr/(4.4410e+07));
        // MaxFONLL->SetBinContent(cnt,maxerr/(1.7240e+08));

        cnt++;
      }
    }

    FONLL->SetLineColor(kRed);
    FONLL->SetMarkerColor(kRed);
    FONLL->SetMarkerStyle(20);
    FONLL->SetMarkerSize(1);

    // MinFONLL->Scale(1./MinFONLL->Integral());
    MinFONLL->SetLineColor(kGreen-2);
    MinFONLL->SetMarkerColor(kGreen-2);
    MinFONLL->SetMarkerStyle(24);

    // MaxFONLL->Scale(1./MaxFONLL->Integral());
    MaxFONLL->SetLineColor(kBlack);
    MaxFONLL->SetMarkerColor(kBlack);
    MaxFONLL->SetMarkerStyle(24);

    MaxFONLL->GetYaxis()->SetTitle("#frac{1}{p_{T}} #frac{d#sigma}{dp_{T}} [mb/GeV]");
    MaxFONLL->GetXaxis()->SetTitle("p_{T}");

    double ptxaxis_published[5] = {1.57, 2.45, 3.44, 4.45, 5.45};
    double dsigmadpt_published[5] = {0.106089296, 0.03107972, 0.00756112, 0.0016320464, 0.00065097852};
    double dsigmadpt_fonll[5] = {1.0515E-01, 3.5115E-02, 9.5505E-03, 2.8380E-03, 9.9716E-04};

    double dsigmaptdpt_published[5]; 
    double dsigmaptdpt_fonll[5]; 

    for (int i = 0; i < 5; i++){
      dsigmaptdpt_published[i] = dsigmadpt_published[i]/ptxaxis_published[i];
      dsigmaptdpt_fonll[i] = dsigmadpt_fonll[i]/ptxaxis_published[i];
    }

    TGraph *publishedspectra = new TGraph(6, ptxaxis_published, dsigmaptdpt_published);
    TGraph *publishedfonll = new TGraph(6, ptxaxis_published, dsigmaptdpt_fonll);
    publishedspectra->SetLineColor(kBlack);
    publishedspectra->SetMarkerColor(kBlack);
    publishedspectra->SetMarkerStyle(22);

    publishedfonll->SetLineColor(kRed);
    publishedfonll->SetMarkerColor(kRed);
    publishedfonll->SetMarkerStyle(26);



    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetMarkerStyle(20);

    h->Scale(1./double(nev));
    // h->Scale(1./105426.);
    h->Scale(1./(0.565*0.5));
    h->Scale(30.);

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->cd();
    gPad->SetLogy();

    MaxFONLL->GetYaxis()->SetRangeUser(5*pow(10, -6), 2);
    MaxFONLL->GetXaxis()->SetRangeUser(0.5, 10);

    MaxFONLL->Draw("C SAME");
    MaxFONLL->SetFillColor(kGreen);
    MaxFONLL->SetLineColor(kGreen);

    MinFONLL->Draw("C SAME");
    MinFONLL->SetFillColor(10);
    MinFONLL->SetLineColor(kGreen);

    FONLL->Draw("P SAME");

    h->Draw("EP SAME");

    publishedspectra->Draw("P SAME");
    publishedfonll->Draw("P SAME");

    c->SaveAs("PythiaSim.pdf");

    
  // }

}



