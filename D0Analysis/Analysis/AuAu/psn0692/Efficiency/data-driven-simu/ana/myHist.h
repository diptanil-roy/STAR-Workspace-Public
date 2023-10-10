#ifndef MYHIST_H_INCLUDED
#define MYHIST_H_INCLUDED
#include "TH1F.h"
#include "TH1D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
//hist
//test
TH1D* hD0;
TH1D* hD0Wg;
TH1F* hpt;
TH1F* hptWg;
TH1F* hppt;
TH1F* hpptWg;

//efficiency
TH2F* h2Pt; //default--clean + default topo cuts + daughter pT>1
TH2F* h2PtCut;
TH2F* h2PtCut_acc;
TH2F* h2PtCut_tpc;
TH2F* h2PtCut_pid;
TH2F* h2PtCut_hft;
TH2F* h2PtCut_hftTopo;
TH2F* h2PtCut_clean;
TH2F* h2PtCut_pt1;
TH2F* h2PtCut_pt2;
TH2F* h2PtCut_topoCut1;
TH2F* h2PtCut_topoCut2;
TH2F* h2PtCut_pt1topoCut1;
TH2F* h2PtCut_pt1topoCut2;

TH3F* h3PtCentY;
TH3F* h3PtCentYCut;
TH3F* h3PtCentYCut_clean;
TH3F* h3PtCentYCut_pt1;
TH3F* h3PtCentYCut_pt2;
TH3F* h3PtCentYCut_topoCut1;
TH3F* h3PtCentYCut_topoCut2;
TH3F* h3PtCentYCut_pt1topoCut1;
TH3F* h3PtCentYCut_pt1topoCut2;

//cuts
namespace anaCuts
{
    const int nCent = 5;
    const float CentEdge[nCent+1] = { -0.5, 1.5, 3.5, 5.5, 6.5, 8.5 };
    const char nameCent[nCent][100] = {"60-80%", "40-60%", "20-40%", "10-20%", "0-10%"};

   // int   const nPtBins = 5;
   // float const PtEdge[nPtBins + 1] = {0., 1., 2., 3., 5., 10.}; // for topo cuts
    const int nPtBins = 6;
    // const float PtEdge[nPtBins+1] = {0, 0.5, 1., 2., 3., 5., 10.};
    const float PtEdge[nPtBins+1] = {0, 0.5, 1., 2., 3., 5., 12.};

   //default cuts
   // float const dcaV0ToPv[nPtBins] = {0.0061, 0.0049, 0.0038, 0.0038, 0.0040};
   // float const decayLength[nPtBins] = {0.0145, 0.0181, 0.0212, 0.0247, 0.0259};
   // float const cosTheta[nPtBins] = {0.95, 0.95, 0.95, 0.95, 0.95};
   // float const dcaDaughters[nPtBins] = {0.0084, 0.0066, 0.0057, 0.0050, 0.0060}; //0.0050;
   // float const kDca[nPtBins] = {0.0103, 0.0091, 0.0095, 0.0079, 0.0058};//0.008, // minimum
   // float const pDca[nPtBins] = {0.0110, 0.0111, 0.0086, 0.0081, 0.0062};//0.008
    float const kDca[nCent][nPtBins] = {
        {0.0106, 0.0106, 0.0069, 0.0068, 0.0050, 0.0050}, // 60-80%
        {0.0140, 0.0100, 0.0075, 0.0072, 0.0060, 0.0050}, // 40-60%
        {0.0151, 0.0102, 0.0104, 0.0099, 0.0063, 0.0050}, // 20-40%
        {0.0145, 0.0113, 0.0094, 0.0089, 0.0069, 0.0050}, // 10-20%
        {0.0138, 0.0109, 0.0082, 0.0094, 0.0076, 0.0054}  // 0-10%
    };
    float const pDca[nCent][nPtBins] = {
        {0.0098, 0.0098, 0.0083, 0.0073, 0.0056, 0.0050}, // 60-80%
        {0.0145, 0.0128, 0.0072, 0.0079, 0.0060, 0.0051}, // 40-60%
        {0.0131, 0.0113, 0.0099, 0.0106, 0.0065, 0.0052}, // 20-40%
        {0.0141, 0.0100, 0.0074, 0.0077, 0.0066, 0.0052}, // 10-20%
        {0.0133, 0.0105, 0.0093, 0.0097, 0.0067, 0.0055}  // 0-10%
    };
    float const dcaV0ToPv[nCent][nPtBins] = {
        {0.0076, 0.0076, 0.0053, 0.0054, 0.0054, 0.0042}, // 60-80%
        {0.0072, 0.0057, 0.0058, 0.0049, 0.0049, 0.0047}, // 40-60%
        {0.0066, 0.0055, 0.0053, 0.0046, 0.0041, 0.0050}, // 20-40%
        {0.0063, 0.0047, 0.0045, 0.0046, 0.0042, 0.0044}, // 10-20%
        {0.0062, 0.0055, 0.0040, 0.0040, 0.0040, 0.0044}  // 0-10%
    };
    float const dcaDaughters[nCent][nPtBins] = {
        {0.0077, 0.0077, 0.0094, 0.0078, 0.0081, 0.0120}, // 60-80%
        {0.0080, 0.0083, 0.0092, 0.0081, 0.0094, 0.0106}, // 40-60%
        {0.0078, 0.0073, 0.0080, 0.0093, 0.0096, 0.0103}, // 20-40%
        {0.0076, 0.0078, 0.0092, 0.0072, 0.0086, 0.0085}, // 10-20%
        {0.0071, 0.0064, 0.0070, 0.0063, 0.0082, 0.0080}  // 0-10%
    };
    float const decayLength[nCent][nPtBins] = {
        {0.0175, 0.0175, 0.0187, 0.0178, 0.0184, 0.0187}, // 60-80%
        {0.0171, 0.0196, 0.0210, 0.0187, 0.0190, 0.0214}, // 40-60%
        {0.0178, 0.0206, 0.0221, 0.0209, 0.0219, 0.0240}, // 20-40%
        {0.0172, 0.0215, 0.0252, 0.0232, 0.0236, 0.0237}, // 10-20%
        {0.0100, 0.0199, 0.0227, 0.0232, 0.0236, 0.0255}  // 0-10%
    };
    float const cosTheta[nCent][nPtBins] = {
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 60-80%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 40-60%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 20-40%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 10-20%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}  // 0-10%
    };
    

   //tight cuts
   // float const dcaV0ToPv1[nPtBins] = {0.0044, 0.0036, 0.0031, 0.0026, 0.0032};
   // float const decayLength1[nPtBins] = {0.0144, 0.0204, 0.0242, 0.0245, 0.0300};
   // float const cosTheta1[nPtBins] = {0.95, 0.95, 0.95, 0.95, 0.95};
   // float const dcaDaughters1[nPtBins] = {0.0069, 0.0048, 0.0044, 0.0049, 0.0047}; //0.0050;
   // float const kDca1[nPtBins] = {0.0119, 0.0110, 0.0109, 0.0106, 0.0080};//0.008, // minimum
   // float const pDca1[nPtBins] = {0.0120, 0.0102, 0.0118, 0.0109, 0.0096};//0.008
    float const kDca1[nCent][nPtBins] = {
        {0.0126, 0.0126, 0.0116, 0.0097, 0.0076, 0.0050}, // 60-80%
        {0.0176, 0.0100, 0.0123, 0.0092, 0.0068, 0.0056}, // 40-60%
        {0.0158, 0.0123, 0.0133, 0.0125, 0.0103, 0.0053}, // 20-40%
        {0.0172, 0.0165, 0.0119, 0.0125, 0.0091, 0.0057}, // 10-20%
        {0.0170, 0.0109, 0.0117, 0.0109, 0.0111, 0.0059}  // 0-10%
    };
    float const pDca1[nCent][nPtBins] = {
        {0.0130, 0.0130, 0.0130, 0.0095, 0.0097, 0.0086}, // 60-80%
        {0.0143, 0.0100, 0.0072, 0.0145, 0.0113, 0.0095}, // 40-60%
        {0.0150, 0.0113, 0.0063, 0.0149, 0.0107, 0.0105}, // 20-40%
        {0.0172, 0.0100, 0.0144, 0.0133, 0.0130, 0.0092}, // 10-20%
        {0.0139, 0.0148, 0.0093, 0.0133, 0.0080, 0.0088}  // 0-10%
    };
    float const dcaV0ToPv1[nCent][nPtBins] = {
        {0.0058, 0.0058, 0.0051, 0.0036, 0.0036, 0.0031}, // 60-80%
        {0.0056, 0.0045, 0.0050, 0.0037, 0.0029, 0.0026}, // 40-60%
        {0.0062, 0.0060, 0.0037, 0.0038, 0.0030, 0.0026}, // 20-40%
        {0.0056, 0.0049, 0.0039, 0.0036, 0.0033, 0.0028}, // 10-20%
        {0.0055, 0.0053, 0.0037, 0.0027, 0.0026, 0.0025}  // 0-10%
    };
    float const dcaDaughters1[nCent][nPtBins] = {
        {0.0076, 0.0076, 0.0078, 0.0095, 0.0076, 0.0093}, // 60-80%
        {0.0088, 0.0059, 0.0081, 0.0083, 0.0073, 0.0108}, // 40-60%
        {0.0061, 0.0043, 0.0070, 0.0089, 0.0066, 0.0070}, // 20-40%
        {0.0077, 0.0049, 0.0042, 0.0056, 0.0053, 0.0119}, // 10-20%
        {0.0077, 0.0044, 0.0047, 0.0073, 0.0060, 0.0061}  // 0-10%
    };
    float const decayLength1[nCent][nPtBins] = {
        {0.0203, 0.0203, 0.0206, 0.0228, 0.0161, 0.0216}, // 60-80%
        {0.0222, 0.0229, 0.0269, 0.0236, 0.0232, 0.0182}, // 40-60%
        {0.0240, 0.0242, 0.0268, 0.0319, 0.0176, 0.0338}, // 20-40%
        {0.0219, 0.0240, 0.0213, 0.0231, 0.0261, 0.0399}, // 10-20%
        {0.0100, 0.0230, 0.0268, 0.0292, 0.0249, 0.0225}  // 0-10%
    };
    float const cosTheta1[nCent][nPtBins] = {
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 60-80%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 40-60%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 20-40%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 10-20%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}  // 0-10%
    };

   //loose cuts
   // float const dcaV0ToPv2[nPtBins] = {0.0072, 0.0053, 0.0047, 0.0042, 0.0062};
   // float const decayLength2[nPtBins] = {0.0110, 0.0168, 0.0187, 0.0199, 0.0180};
   // float const cosTheta2[nPtBins] = {0.95, 0.95, 0.95, 0.95, 0.95};
   // float const dcaDaughters2[nPtBins] = {0.0077, 0.0078, 0.0074, 0.0068, 0.0066}; //0.0050;
   // float const kDca2[nPtBins] = {0.0105, 0.0068, 0.0080, 0.0066, 0.0041};//0.008, // minimum
   // float const pDca2[nPtBins] = {0.0092, 0.0078, 0.0086, 0.0065, 0.0047};//0.008
    float const kDca2[nCent][nPtBins] = {
        {0.0090, 0.0090, 0.0073, 0.0057, 0.0050, 0.0050}, // 60-80%
        {0.0115, 0.0114, 0.0067, 0.0050, 0.0050, 0.0050}, // 40-60%
        {0.0134, 0.0112, 0.0074, 0.0063, 0.0050, 0.0050}, // 20-40%
        {0.0135, 0.0120, 0.0070, 0.0060, 0.0050, 0.0050}, // 10-20%
        {0.0111, 0.0109, 0.0080, 0.0067, 0.0050, 0.0050}  // 0-10%
    };
    float const pDca2[nCent][nPtBins] = {
        {0.0098, 0.0098, 0.0080, 0.0058, 0.0050, 0.0050}, // 60-80%
        {0.0124, 0.0114, 0.0072, 0.0050, 0.0050, 0.0050}, // 40-60%
        {0.0111, 0.0113, 0.0089, 0.0062, 0.0050, 0.0050}, // 20-40%
        {0.0108, 0.0100, 0.0074, 0.0069, 0.0050, 0.0050}, // 10-20%
        {0.0116, 0.0105, 0.0093, 0.0072, 0.0050, 0.0050}  // 0-10%
    };
    float const dcaV0ToPv2[nCent][nPtBins] = {
        {0.0083, 0.0083, 0.0090, 0.0077, 0.0059, 0.0055}, // 60-80%
        {0.0069, 0.0083, 0.0087, 0.0077, 0.0085, 0.0072}, // 40-60%
        {0.0075, 0.0073, 0.0056, 0.0053, 0.0100, 0.0090}, // 20-40%
        {0.0071, 0.0063, 0.0055, 0.0052, 0.0065, 0.0055}, // 10-20%
        {0.0078, 0.0063, 0.0054, 0.0045, 0.0062, 0.0047}  // 0-10%
    };
    float const dcaDaughters2[nCent][nPtBins] = {
        {0.0097, 0.0097, 0.0123, 0.0092, 0.0093, 0.0092}, // 60-80%
        {0.0095, 0.0084, 0.0100, 0.0108, 0.0130, 0.0113}, // 40-60%
        {0.0081, 0.0090, 0.0078, 0.0094, 0.0137, 0.0121}, // 20-40%
        {0.0082, 0.0086, 0.0099, 0.0102, 0.0121, 0.0105}, // 10-20%
        {0.0098, 0.0087, 0.0088, 0.0078, 0.0100, 0.0101}  // 0-10%
    };
    float const decayLength2[nCent][nPtBins] = {
        {0.0154, 0.0154, 0.0163, 0.0147, 0.0126, 0.0140}, // 60-80%
        {0.0158, 0.0153, 0.0172, 0.0150, 0.0126, 0.0164}, // 40-60%
        {0.0100, 0.0177, 0.0177, 0.0194, 0.0131, 0.0150}, // 20-40%
        {0.0157, 0.0197, 0.0222, 0.0180, 0.0155, 0.0189}, // 10-20%
        {0.0148, 0.0179, 0.0215, 0.0198, 0.0167, 0.0206}  // 0-10%
    };
    float const cosTheta2[nCent][nPtBins] = {
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 60-80%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 40-60%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 20-40%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}, // 10-20%
      {0.95, 0.95, 0.95, 0.95, 0.95, 0.95}  // 0-10%
    };
}

//define hist
void DefineHist()
{
   TH1::SetDefaultSumw2();
   hD0 = new TH1D("hD0", "hD0", 2, 0, 2);
   hD0Wg = new TH1D("hD0Wg", "hD0Wg", 2, 0, 2);
   hpt = new TH1F("hpt", "hpt", 500, 0, 10);
   hptWg = new TH1F("hptWg", "hptWg", 500, 0, 10);
   hppt = new TH1F("hppt", "hppt", 500, 0, 10);
   hpptWg = new TH1F("hpptWg", "hpptWg", 500, 0, 10);

   h2Pt = new TH2F("h2Pt", "h2Pt", 200, 0, 10, 10, 0, 10);
   h2PtCut = new TH2F("h2PtCut", "h2PtCut", 200, 0, 10, 10, 0, 10);
   h2PtCut_acc = new TH2F("h2PtCut_acc", "h2PtCut_acc", 200, 0, 10, 10, 0, 10);
   h2PtCut_tpc = new TH2F("h2PtCut_tpc", "h2PtCut_tpc", 200, 0, 10, 10, 0, 10);
   h2PtCut_pid = new TH2F("h2PtCut_pid", "h2PtCut_pid", 200, 0, 10, 10, 0, 10);
   h2PtCut_hft = new TH2F("h2PtCut_hft", "h2PtCut_hft", 200, 0, 10, 10, 0, 10);
   h2PtCut_hftTopo = new TH2F("h2PtCut_hftTopo", "h2PtCut_hftTopo", 200, 0, 10, 10, 0, 10);
   h2PtCut_clean = new TH2F("h2PtCut_clean", "h2PtCut_clean", 200, 0, 10, 10, 0, 10);
   h2PtCut_pt1 = new TH2F("h2PtCut_pt1", "h2PtCut_pt1", 200, 0, 10, 10, 0, 10);
   h2PtCut_pt2 = new TH2F("h2PtCut_pt2", "h2PtCut_pt2", 200, 0, 10, 10, 0, 10);
   h2PtCut_topoCut1 = new TH2F("h2PtCut_topoCut1", "h2PtCut_topoCut1", 200, 0, 10, 10, 0, 10);
   h2PtCut_topoCut2 = new TH2F("h2PtCut_topoCut2", "h2PtCut_topoCut2", 200, 0, 10, 10, 0, 10);
   h2PtCut_pt1topoCut1 = new TH2F("h2PtCut_pt1topoCut1", "h2PtCut_pt1topoCut1", 200, 0, 10, 10, 0, 10);
   h2PtCut_pt1topoCut2 = new TH2F("h2PtCut_pt1topoCut2", "h2PtCut_pt1topoCut2", 200, 0, 10, 10, 0, 10);

   h3PtCentY = new TH3F("h3PtCentY", "no Cuts;p_{T};cent;y", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut = new TH3F("h3PtCentYCut", "with all cuts;p_{T};cent;y", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut_topoCut1 = new TH3F("h3PtCentYCut_topoCut1", "with all cuts (50% eff topo.);p_{T};cent;y", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut_topoCut2 = new TH3F("h3PtCentYCut_topoCut2", "with all cuts (150% eff topo.);p_{T};cent;y", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut_clean = new TH3F("h3PtCentYCut_clean", "h3PtCentYCut_clean", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut_pt1 = new TH3F("h3PtCentYCut_pt1", "h3PtCentYCut_pt1", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut_pt2 = new TH3F("h3PtCentYCut_pt2", "h3PtCentYCut_pt2", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut_pt1topoCut1 = new TH3F("h3PtCentYCut_pt1topoCut1", "h3PtCentYCut_pt1topoCut1", 200, 0, 10, 9, 0, 9, 20, -1, 1);
   h3PtCentYCut_pt1topoCut2 = new TH3F("h3PtCentYCut_pt1topoCut2", "h3PtCentYCut_pt1topoCut2", 200, 0, 10, 9, 0, 9, 20, -1, 1);
}
//write hist
void WriteHist(TFile* fout)
{
   fout->cd();
   h3PtCentY->Write();
   h3PtCentYCut->Write();
   h3PtCentYCut_clean->Write();
   h3PtCentYCut_pt1->Write();
   h3PtCentYCut_pt2->Write();
   h3PtCentYCut_topoCut1->Write();
   h3PtCentYCut_topoCut2->Write();
   h3PtCentYCut_pt1topoCut1->Write();
   h3PtCentYCut_pt1topoCut2->Write();

   hD0->Write();
   hD0Wg->Write();
   hpt->Write();
   hptWg->Write();
   hppt->Write();
   hpptWg->Write();

   h2Pt->Write();
   h2PtCut->Write();
   h2PtCut_acc->Write();
   h2PtCut_tpc->Write();
   h2PtCut_pid->Write();
   h2PtCut_hft->Write();
   h2PtCut_hftTopo->Write();
   h2PtCut_clean->Write();
   h2PtCut_pt1->Write();
   h2PtCut_pt2->Write();
   h2PtCut_topoCut1->Write();
   h2PtCut_topoCut2->Write();
   h2PtCut_pt1topoCut1->Write();
   h2PtCut_pt1topoCut2->Write();

}
void DeleteHist()
{
   if (h3PtCentY)
   {
      delete h3PtCentY;
      h3PtCentY = NULL;
   }
   if (h3PtCentYCut)
   {
      delete h3PtCentYCut;
      h3PtCentYCut = NULL;
   }
   if (h3PtCentYCut_clean)
   {
      delete h3PtCentYCut_clean;
      h3PtCentYCut_clean = NULL;
   }
   if (h3PtCentYCut_pt1)
   {
      delete h3PtCentYCut_pt1;
      h3PtCentYCut_pt1 = NULL;
   }
   if (h3PtCentYCut_pt2)
   {
      delete h3PtCentYCut_pt2;
      h3PtCentYCut_pt2 = NULL;
   }
   if (h3PtCentYCut_topoCut1)
   {
      delete h3PtCentYCut_topoCut1;
      h3PtCentYCut_topoCut1 = NULL;
   }
   if (h3PtCentYCut_topoCut2)
   {
      delete h3PtCentYCut_topoCut2;
      h3PtCentYCut_topoCut2 = NULL;
   }
   if (h3PtCentYCut_pt1topoCut1)
   {
      delete h3PtCentYCut_pt1topoCut1;
      h3PtCentYCut_pt1topoCut1 = NULL;
   }
   if (h3PtCentYCut_pt1topoCut2)
   {
      delete h3PtCentYCut_pt1topoCut2;
      h3PtCentYCut_pt1topoCut2 = NULL;
   }
   if (hD0)
   {
      delete hD0;
      hD0 = NULL;
   }
   if (hD0Wg)
   {
      delete hD0Wg;
      hD0Wg = NULL;
   }
   if (hpt)
   {
      delete hpt;
      hpt = NULL;
   }
   if (hptWg)
   {
      delete hptWg;
      hptWg = NULL;
   }
   if (hppt)
   {
      delete hppt;
      hppt = NULL;
   }
   if (hpptWg)
   {
      delete hpptWg;
      hpptWg = NULL;
   }

   if (h2Pt)
   {
      delete h2Pt;
      h2Pt = NULL;
   }
   if (h2PtCut)
   {
      delete h2PtCut;
      h2PtCut = NULL;
   }
   if (h2PtCut_acc)
   {
      delete h2PtCut_acc;
      h2PtCut_acc = NULL;
   }
   if (h2PtCut_tpc)
   {
      delete h2PtCut_tpc;
      h2PtCut_tpc = NULL;
   }
   if (h2PtCut_pid)
   {
      delete h2PtCut_pid;
      h2PtCut_pid = NULL;
   }
   if (h2PtCut_hft)
   {
      delete h2PtCut_hft;
      h2PtCut_hft = NULL;
   }
   if (h2PtCut_hftTopo)
   {
      delete h2PtCut_hftTopo;
      h2PtCut_hftTopo = NULL;
   }
   if (h2PtCut_clean)
   {
      delete h2PtCut_clean;
      h2PtCut_clean = NULL;
   }
   if (h2PtCut_pt1)
   {
      delete h2PtCut_pt1;
      h2PtCut_pt1 = NULL;
   }
   if (h2PtCut_pt2)
   {
      delete h2PtCut_pt2;
      h2PtCut_pt2 = NULL;
   }
   if (h2PtCut_topoCut1)
   {
      delete h2PtCut_topoCut1;
      h2PtCut_topoCut1 = NULL;
   }
   if (h2PtCut_topoCut2)
   {
      delete h2PtCut_topoCut2;
      h2PtCut_topoCut2 = NULL;
   }
   if (h2PtCut_pt1topoCut1)
   {
      delete h2PtCut_pt1topoCut1;
      h2PtCut_pt1topoCut1 = NULL;
   }
   if (h2PtCut_pt1topoCut2)
   {
      delete h2PtCut_pt1topoCut2;
      h2PtCut_pt1topoCut2 = NULL;
   }

}
#endif // VARIABLE_H_INCLUDED
