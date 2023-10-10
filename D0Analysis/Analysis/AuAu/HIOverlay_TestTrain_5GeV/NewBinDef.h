#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"
// #include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
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
// #include "StJetTreeStruct.h"
#include <vector>

#endif


// Calling all Histograms that will be saved out

/*

const int njpt_gen_bins = 15; //x2 MC
double jetpt_gen_low = 5.;
double jetpt_gen_high = 20.;

const int nz_gen_bins = 10; //x2 MC
double z_gen_low = 0.2;
double z_gen_high = 1.;

const int njpt_bins = 10; 
double jetpt_low = 0.;
double jetpt_high = 30.;

const int nz_bins = 20;
double z_low = 0.0;
double z_high = 4.0;

const int ndim = 4;
int nbins[ndim] = {njpt_gen_bins, nz_gen_bins, njpt_bins, nz_bins};
double xmin[ndim] = {jetpt_gen_low, z_gen_low, jetpt_low, z_low};
double xmax[ndim] = {jetpt_gen_high, z_gen_high, jetpt_high, z_high};

const int ncentbin = 9;
double centbin[ncentbin + 1] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

// For plotting and as such
const int nbins_jpt = 6;
double binning_jpt[nbins_jpt+1] = {5,7,9,11,13,15,20};

*/

// Everything above this works for 1D. I am testing things out for 2D now.

const int njpt_gen_bins = 15; //x2 MC
double jetpt_gen_low = 5.;
double jetpt_gen_high = 20.;

const int njpt_gen_bins_var = 6; //x2 MC
double jetpt_var_bin[njpt_gen_bins_var+1] = {5,7,9,11,13,15,20};

const int nz_gen_bins = 10; //x2 MC
double z_gen_low = 0.;
double z_gen_high = 1.;

const int njpt_bins = 9; 
double jetpt_low = 3.;
double jetpt_high = 30.;

const int nz_bins = 5;
double z_low = 0.0;
double z_high = 2.5;

const int ndim = 4;
int nbins[ndim] = {njpt_gen_bins, nz_gen_bins, njpt_bins, nz_bins};
double xmin[ndim] = {jetpt_gen_low, z_gen_low, jetpt_low, z_low};
double xmax[ndim] = {jetpt_gen_high, z_gen_high, jetpt_high, z_high};

const int ncentbin = 9;
double centbin[ncentbin + 1] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

// For plotting and as such
const int nbins_jpt = 6;
double binning_jpt[nbins_jpt+1] = {5,7,9,11,13,15,20};

// const int njpt_gen_bins = 60; //x2 MC
// double jetpt_gen_low = 5.;
// double jetpt_gen_high = 20.;

// const int njpt_gen_bins_var = 6; //x2 MC
// double jetpt_var_bin[njpt_gen_bins_var+1] = {5,7,9,11,13,15,20};

// const int nz_gen_bins = 28; //x2 MC
// double z_gen_low = 0.3;
// double z_gen_high = 1.;

// const int njpt_bins = 9; 
// double jetpt_low = 3.;
// double jetpt_high = 30.;

// const int nz_bins = 5;
// double z_low = 0.0;
// double z_high = 2.5;

// const int ndim = 4;
// int nbins[ndim] = {njpt_gen_bins, nz_gen_bins, njpt_bins, nz_bins};
// double xmin[ndim] = {jetpt_gen_low, z_gen_low, jetpt_low, z_low};
// double xmax[ndim] = {jetpt_gen_high, z_gen_high, jetpt_high, z_high};

// const int ncentbin = 9;
// double centbin[ncentbin + 1] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

// // For plotting and as such
// const int nbins_jpt = 6;
// double binning_jpt[nbins_jpt+1] = {5,7,9,11,13,15,20};

