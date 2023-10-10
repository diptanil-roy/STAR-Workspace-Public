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

double mczcutoff = 0.0;

int countcutoff = 0;

const int njpt_gen_bins_var = 8; //x2 MC
double jetpt_var_bin[njpt_gen_bins_var+1] = {5,7,9,11,13,15,20,25,30};
double jetpt_gen_binwidth_low = 2.;
double jetpt_gen_binwidth_high = 5.;

// const int njpt_gen_bins_var = 8; //x2 MC
// double jetpt_var_bin[njpt_gen_bins_var+1] = {5,7,9,11,13,15,20,25,30};
// double jetpt_gen_binwidth_low = 2.;
// double jetpt_gen_binwidth_high = 5.;

const int nz_gen_bins = 7; //x2 MC
double z_gen_bin[nz_gen_bins+1] = {0., 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.};
double z_gen_binwidth_low = 0.1;
double z_gen_binwidth_high = 0.1;

const int njpt_bins = 16;
double nbinsjetpt[njpt_bins + 1] = {1, 2, 3, 4, 5,6, 7, 8, 9, 10, 12, 15, 18, 21, 24, 30, 50};

double jetpt_binwidth_low = 1.;
double jetpt_binwidth_high = 20.;

// const int njpt_bins = 12;
// double nbinsjetpt[njpt_bins + 1] = {5, 6, 7, 8, 9, 10, 12, 15, 18, 21, 24, 30, 50};

// double jetpt_binwidth_low = 1.;
// double jetpt_binwidth_high = 20.;

const int nz_bins = 9;
double nbinsz[nz_bins + 1] = {0., 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

double z_binwidth_low = 0.2;
double z_binwidth_high = 0.1;

const int njpt_bins_wide = 20;
double nbinsjetpt_wide[njpt_bins_wide + 1] = {-25, -20, -15, -10, -5, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 35, 40, 50};

double jetpt_binwidth_low_wide = 5.;
double jetpt_binwidth_high_wide = 10.;

const int nz_bins_wide = 18;
double nbinsz_wide[nz_bins_wide + 1] = {-100, -50, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 50, 100};

double z_binwidth_low_wide = 50;
double z_binwidth_high_wide = 50;

// These ae for Area Based Wide Range (I want to club all the negative bins together, things that go to negative values should be ashamed of themselves)

// const int njpt_bins_wide = 20;
// double nbinsjetpt_wide[njpt_bins_wide + 1] = {-25, -20, -15, -10, -5, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 35, 40, 50};

// double jetpt_binwidth_low_wide = 1.;
// double jetpt_binwidth_high_wide = 10.;

// const int nz_bins_wide = 11;
// double nbinsz_wide[nz_bins_wide + 1] = {-10., 0, 0.2, 0.4, 0.6, 1.0, 1.5, 2.0, 2.5, 5, 7.5, 10};

// double z_binwidth_low_wide = 1;
// double z_binwidth_high_wide = 2.5;

const int ndim = 4;
int nbins[ndim] = {njpt_gen_bins_var, nz_gen_bins, njpt_bins, nz_bins};
int nbins_wide[ndim] = {njpt_gen_bins_var, nz_gen_bins, njpt_bins_wide, nz_bins_wide};

const int ncentbin = 9;
double centbin[ncentbin + 1] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

// For plotting and as such
const int nbins_jpt = 6;
double binning_jpt[nbins_jpt+1] = {5,7,9,11,13,15,20};

const int njetconstbin = 17;
double jetconstbin[njetconstbin+1] = {1, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80};

const int ndrbins = 8;
double drbins[ndrbins+1] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6};
