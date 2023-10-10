// Jet PT distribtuions
//const int nbins_jpt = 32;
//double binning_jpt[nbins_jpt+1] = {-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,7,9,11,13,15,20,30};
// const int nbins_jpt =13;
// double binning_jpt[nbins_jpt+1] = {0,1,2,3,4,5,7,9,11,13,15,20,30,50};

// const int nbins_jpt =8;
// double binning_jpt[nbins_jpt+1] = {0, 3, 5, 7, 9, 11, 13, 15, 20};
// const int nbins_jpt =9;
// double binning_jpt[nbins_jpt+1] = {1, 3, 5, 7, 9, 11, 13, 15, 20, 30};

const int nbins_jpt = 9;
// double binning_jpt[nbins_jpt+1] = {-7, -5, -2, 0, 1, 3, 5, 7, 9, 11, 13, 15, 20, 30};
// double binning_jpt[nbins_jpt+1] = {0, 1, 3, 5, 7, 9, 11, 13, 15, 20};
double binning_jpt[nbins_jpt+1] = {1, 3, 5, 7, 9, 11, 13, 15, 20, 30};

// const int nbins_z = 34;
// double binning_z[nbins_z+1] = {-100,-10,-4,-2,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.5,4,5,10,100}; 

// const int nbins_z = 10;
// double binning_z[nbins_z+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.};

const int nbins_z = 9;
double binning_z[nbins_z+1] = {0., 0.2, 0.4, 0.6, 0.8, 1.001, 1.5, 2., 2.5, 3.}; 

//const Int_t nBinsdR = 5;
//Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4};

const Int_t PTBINS = 11;
Double_t edges[PTBINS+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

//const Int_t nBinsdR = 7;
//Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 1.0};

const Int_t nBinsdR = 7;
Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 1.0};

const int nBinsJetPt_R = 16;
double JetPtBins_R[nBinsJetPt_R+1] = {-30., -20., -10., 0., 1., 2., 3.,4.,5.,7.,9.,11.,13.,15.,20.,30.,50.};


const int nBinsJetPt = 16;
double JetPtBins[nBinsJetPt+1] = {-30., -20., -10., 0., 1., 2., 3.,4.,5.,7.,9.,11.,13.,15.,20.,30.,50.};


// This is added on Nov 8, 2022.
// const int njpt_bins = 10; 
// double jetpt_low = 0.;
// double jetpt_high = 30.;

// const int nz_bins = 20; // Changed from 40. Let's see what it does.
// double z_low = 0.0;
// double z_high = 4.0;


// These settings work for 2D
// const int njpt_bins = 30; 
// double jetpt_low = 0.;
// double jetpt_high = 30.;

// const int nz_bins = 10; // Changed from 40. Let's see what it does.
// double z_low = 0.0;
// double z_high = 2.5;

// const int njpt_bins = 21;
// double nbinsjetpt[njpt_bins + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 21, 23, 25, 30, 35};

// const int njpt_bins = 18;
// double nbinsjetpt[njpt_bins + 1] = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 21, 23, 25, 30, 35, 50};

// const int nz_bins = 10;
// double nbinsz[nz_bins + 1] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

/// May 17, 2023

// const int njpt_bins = 10;
// double nbinsjetpt[njpt_bins + 1] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30};

const int njpt_bins = 20;
double nbinsjetpt[njpt_bins + 1] = {-25, -20, -15, -10, -5, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 35, 40, 50};

const int nz_bins = 11;
double nbinsz[nz_bins + 1] = {-5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5};

// const int njpt_bins = 12;
// double nbinsjetpt[njpt_bins + 1] = {5, 7, 9, 12, 15, 18, 21, 24, 27, 30, 35, 40, 50};

// const int nz_bins = 5;
// double nbinsz[nz_bins + 1] = {0, 0.2, 0.4, 0.6, 0.8, 1.};
