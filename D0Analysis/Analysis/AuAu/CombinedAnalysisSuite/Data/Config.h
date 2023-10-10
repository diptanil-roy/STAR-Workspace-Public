const int nbins_jpt = 9;
double binning_jpt[nbins_jpt+1] = {1, 3, 5, 7, 9, 11, 13, 15, 20, 30};

const int nbins_z = 9;
double binning_z[nbins_z+1] = {0., 0.2, 0.4, 0.6, 0.8, 1.001, 1.5, 2., 2.5, 3.}; 

//const Int_t nBinsdR = 5;
//Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.2, 0.3, 0.4};

const Int_t PTBINS = 11;
Double_t edges[PTBINS+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

// const Int_t PTBINS = 5;
// Double_t edges[PTBINS+1] = {1.0, 2.0, 3.0, 4.0, 5.0, 10.0};

const Int_t nBinsdR = 5;
Double_t dRBins[nBinsdR+1] = {0.0, 0.05, 0.1, 0.2, 0.4, 0.6};


// These ae for Area Based Wide Range

const int njpt_bins_wide = 16;
double nbinsjetpt_wide[njpt_bins_wide + 1] = {-10, -5, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 40, 50};

const int nz_bins_wide = 24;
double nbinsz_wide[nz_bins_wide + 1] =  {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000};

// // These ae for Area Based Wide Range (I want to club all the negative bins together, things that go to negative values should be ashamed of themselves)

// const int njpt_bins_wide = 20;
// double nbinsjetpt_wide[njpt_bins_wide + 1] = {-25, -20, -15, -10, -5, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 35, 40, 50};

// const int nz_bins_wide = 11;
// double nbinsz_wide[nz_bins_wide + 1] = {-10., 0, 0.2, 0.4, 0.6, 1.0, 1.5, 2.0, 2.5, 5, 7.5, 10};

// These are for CS

const int njpt_bins_cs = 16;
double nbinsjetptcs[njpt_bins_cs + 1] = {1,2,3,4,5, 6, 7, 8, 9, 10, 12, 15, 18, 21, 24, 30, 50};

const int nz_bins_cs = 9;
double nbinszcs[nz_bins_cs + 1] = {0., 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

// These are for Area Based

const int njpt_bins = 9;
double nbinsjetpt[njpt_bins + 1] = {3, 6, 9, 12, 15, 18, 21, 24, 27, 30};

const int nz_bins = 5;
double nbinsz[nz_bins + 1] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};

// Centrality

const int nCentBins = 9;
double CentBins[nCentBins+1] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};

// NJetConstituents

const int njetconstbin = 17;
double jetconstbin[njetconstbin+1] = {1, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80};

