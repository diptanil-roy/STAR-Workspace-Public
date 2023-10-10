#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include "iostream"
#include <string>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib> //exit function

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#endif

//#include "dataDrivenFastSimulator.h"
#include "d0Nt.h"
#include "myHist.h"
using namespace std;

int getD0PtIndex(float const pt);
int getD0CentIndex(float const cent);
void loadInputFile();

//input 1--D0 pT weight
TF1* fweight;

//input 2--k/pi tof match efficiency
TF1* htof_pip[9];
TF1* htof_kp[9];
TF1* htof_pim[9];
TF1* htof_km[9];

//input 3--k/pi pid eff
TF1* gTpcPID_pi;
TF1* gTpcPID_k;
TF1* gTofPID_pi;
TF1* gTofPID_k;

TH1F* hptWeight = NULL;
int main(int argc, char **argv)
{
   if (argc != 3)
   {
      cout << "Error: the number of arguments is wrong!" << endl;
      cout << "First argument is:\t input root file or file list" << endl;
      cout << "Second argument is:\t output root file name" << endl;
      exit(1);
   }
   std::string inputFile = argv[1];
   std::string outputFileName = argv[2];
   loadInputFile();
   d0Nt* t = new d0Nt(inputFile);

   //double FR1 = 0.4; //b->B+
   //double BR1 = 0.086+0.79;//0.086: B+ -> D0 X; 0.79: B+ -> D0bar X
   //double FR2 = 0.4; //b->B0
   //double BR2 = 0.081+0.474; // 0.081: B0->D0 X, 0.474: B0->D0bar X
   //double R = FR1*BR1 + FR2*BR2;

   if (outputFileName.find(".root", 0) == std::string::npos) outputFileName += ".root";
   TFile* fOut = new TFile(outputFileName.c_str(), "recreate");
   DefineHist();

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;
   int num = 0;
   int num1 = 0;
   int num2 = 0;
   int num3 = 0;
   int num4 = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (!((i + 1) % (nEntries / 10)))
         cout << "____________ ipart = " << i / static_cast<float>(nEntries) << " ________________" << endl;
      num++;

      int centrality = t->cent;
      float primarypt = t->primaryPt;
      float pt = t->pt;
      float rcpt = t->rPt;
      float dca = t->dcaD0ToPv / 1.e4;
      float decayL = t->decayLength / 1.e4;
      float dca12 = t->dca12 / 1.e4;
      float cosTheta = t->cosTheta;
      float ppt = t->pRPt;
      float pp = ppt * cosh(t->pREta);
      float kpt = t->kRPt;
      float kp = kpt * cosh(t->kREta);
      float pdca = t->pRDca / 1.e4;
      float kdca = t->kRDca / 1.e4;

      float weight = fweight->Eval(pt);

      if (fabs(t->rY) < 1.0)
      {
         hpt->Fill(pt);
         hptWg->Fill(pt, weight);
         hppt->Fill(primarypt);
         hpptWg->Fill(primarypt, weight);
      }

      if (t->pid == 421)
      {
         hD0->Fill(0);
         hD0Wg->Fill("D^{0}", weight);
      }
      else if (t->pid == -421)
      {
         hD0->Fill(1);
         hD0Wg->Fill("#bar{D^{0}", weight);
      }
      else
      {
         hD0->Fill(2);
         hD0Wg->Fill(2., weight);
      }

      if (fabs(t->rY) >= 1.0) continue;
      h2Pt->Fill(pt, centrality, weight);
      h3PtCentY->Fill(pt, centrality, t->rY, weight);

      // single eff factor
      bool isGoodTrack_QA = t->kRPt > 0.3 && fabs(t->kREta) < 1 && t->pRPt > 0.3 && fabs(t->pREta) < 1;
      if (isGoodTrack_QA) h2PtCut_acc->Fill(pt, centrality, weight);
      // isGoodTrack_QA = true;
      bool passTpc_QA = t->kTpc > 0 && t->pTpc > 0;
      bool passHft_QA = t->kHft > 0 && t->pHft > 0;
      int ptIndex = getD0PtIndex(t->rPt);
      int centIndex = getD0CentIndex(centrality);
      if (ptIndex < 0) continue;
      bool mtopoCut_QA = (kdca > anaCuts::kDca[centIndex][ptIndex] &&
                          pdca > anaCuts::pDca[centIndex][ptIndex] &&
                          dca12 < anaCuts::dcaDaughters[centIndex][ptIndex] &&
                          decayL > anaCuts::decayLength[centIndex][ptIndex] &&
                          cosTheta > anaCuts::cosTheta[centIndex][ptIndex] &&
                          dca < anaCuts::dcaV0ToPv[centIndex][ptIndex]
                         );

      if (isGoodTrack_QA && passTpc_QA) 
      {
      float dedxEff_tmp = gTpcPID_k->Eval(kpt) * gTpcPID_pi->Eval(ppt); //for kaon <2--0.952, for pion <3--0.997
      float tofMatchEff_k_tmp;
      float tofMatchEff_pi_tmp;
      if (t->pid == 421) 
      {
        tofMatchEff_pi_tmp = htof_pip[centrality]->Eval(ppt);
        tofMatchEff_k_tmp = htof_km[centrality]->Eval(kpt);
      }
      else if (t->pid == -421) 
      {
        tofMatchEff_pi_tmp = htof_pim[centrality]->Eval(ppt);
        tofMatchEff_k_tmp = htof_kp[centrality]->Eval(kpt);
      }
      const float tofMatchEff_tmp = tofMatchEff_pi_tmp * tofMatchEff_k_tmp;

      //tof pid eff.
      const float tofPidEff_k_tmp = gTofPID_k->Eval(kpt);
      const float tofPidEff_pi_tmp = gTofPID_pi->Eval(ppt);
      const float tofPidEff_tmp = tofPidEff_k_tmp * tofPidEff_pi_tmp;

      //combine tof
      const float tof_hybrid_k_tmp = 1.0 - tofMatchEff_k_tmp + tofMatchEff_k_tmp * tofPidEff_k_tmp;
      const float tof_hybrid_pi_tmp = 1.0 - tofMatchEff_pi_tmp + tofMatchEff_pi_tmp * tofPidEff_pi_tmp;
      const float tof_hybrid_tmp = tof_hybrid_pi_tmp * tof_hybrid_k_tmp;
      float pidweight = dedxEff_tmp * tof_hybrid_tmp;

      // if (isGoodTrack_QA && passTpc_QA) h2PtCut_pid->Fill(pt, centrality, weight *pidweight);
      if (isGoodTrack_QA && passTpc_QA) h2PtCut_pid->Fill(pt, centrality, weight *pidweight);
   }


      if (isGoodTrack_QA && passTpc_QA) h2PtCut_tpc->Fill(pt, centrality, weight);
      if (isGoodTrack_QA && passTpc_QA && passHft_QA) h2PtCut_hft->Fill(pt, centrality, weight);
      if (isGoodTrack_QA && passTpc_QA && passHft_QA && mtopoCut_QA) h2PtCut_hftTopo->Fill(pt, centrality, weight);

      //cut
      bool isGoodTrack = t->kRPt > 0.3 && fabs(t->kREta) < 1 && t->pRPt > 0.3 && fabs(t->pREta) < 1;
      if (!isGoodTrack) continue;
      num1++;
      bool passTpc = t->kTpc > 0 && t->pTpc > 0;
      if (!passTpc) continue;
      num2++;
      bool passHft = t->kHft > 0 && t->pHft > 0;
      if (!passHft) continue;
      num3++;
      //int ptIndex = getD0PtIndex(t->rPt);
      //if(ptIndex < 0) continue;
      num4++;

      bool mtopoCut = (kdca > anaCuts::kDca[centIndex][ptIndex] &&
                       pdca > anaCuts::pDca[centIndex][ptIndex] &&
                       dca12 < anaCuts::dcaDaughters[centIndex][ptIndex] &&
                       decayL > anaCuts::decayLength[centIndex][ptIndex] &&
                       cosTheta > anaCuts::cosTheta[centIndex][ptIndex] &&
                       dca < anaCuts::dcaV0ToPv[centIndex][ptIndex]
                      );

      bool mtopoCut1 = (kdca > anaCuts::kDca1[centIndex][ptIndex] &&
                        pdca > anaCuts::pDca1[centIndex][ptIndex] &&
                        dca12 < anaCuts::dcaDaughters1[centIndex][ptIndex] &&
                        decayL > anaCuts::decayLength1[centIndex][ptIndex] &&
                        cosTheta > anaCuts::cosTheta1[centIndex][ptIndex] &&
                        dca < anaCuts::dcaV0ToPv1[centIndex][ptIndex]
                       );

      bool mtopoCut2 = (kdca > anaCuts::kDca2[centIndex][ptIndex] &&
                        pdca > anaCuts::pDca2[centIndex][ptIndex] &&
                        dca12 < anaCuts::dcaDaughters2[centIndex][ptIndex] &&
                        decayL > anaCuts::decayLength2[centIndex][ptIndex] &&
                        cosTheta > anaCuts::cosTheta2[centIndex][ptIndex] &&
                        dca < anaCuts::dcaV0ToPv2[centIndex][ptIndex]
                       );

      bool mPtCut = kpt > 0.6 && ppt > 0.6;
      bool mPtCut1 = kpt > 0.3 && ppt > 0.3;
      bool mPtCut2 = kpt > 0.5 && ppt > 0.5;

      //dedx cut eff.
      // const float dedxEff = 0.952 * 0.997; //for kaon <2--0.952, for pion <3--0.997
      const float dedxEff = gTpcPID_k->Eval(kpt) * gTpcPID_pi->Eval(ppt); //for kaon <2--0.952, for pion <3--0.997

      //tof match eff.
      // const float tofMatchEff_k = kpt < 1.6 ? htof_k[centrality]->GetBinContent(htof_k[centrality]->FindBin(kpt)) : 1.0;
      // const float tofMatchEff_pi = ppt < 1.6 ? htof_pi[centrality]->GetBinContent(htof_pi[centrality]->FindBin(ppt)) : 1.0;
      float tofMatchEff_k;
      float tofMatchEff_pi;
      if (t->pid == 421) 
      {
        tofMatchEff_pi = htof_pip[centrality]->Eval(ppt);
        tofMatchEff_k = htof_km[centrality]->Eval(kpt);
      }
      else if (t->pid == -421) 
      {
        tofMatchEff_pi = htof_pim[centrality]->Eval(ppt);
        tofMatchEff_k = htof_kp[centrality]->Eval(kpt);
      }
      const float tofMatchEff = tofMatchEff_pi * tofMatchEff_k;

      //tof pid eff.
      // const float tofPidEff_k = kp < 1.6 ? gTofPID_k->Eval(kp) : 1.0;
      // const float tofPidEff_pi = pp < 1.6 ? gTofPID_pi->Eval(pp) : 1.0;
      const float tofPidEff_k = gTofPID_k->Eval(kpt);
      const float tofPidEff_pi = gTofPID_pi->Eval(ppt);
      const float tofPidEff = tofPidEff_k * tofPidEff_pi;

      //combine tof
      const float tof_hybrid_k = 1.0 - tofMatchEff_k + tofMatchEff_k * tofPidEff_k;
      const float tof_hybrid_pi = 1.0 - tofMatchEff_pi + tofMatchEff_pi * tofPidEff_pi;
      const float tof_hybrid = tof_hybrid_pi * tof_hybrid_k;

      const float tof_clean_k = tofMatchEff_k * tofPidEff_k;
      const float tof_clean_pi = tofMatchEff_pi * tofPidEff_pi;
      const float tof_clean = tof_clean_k * tof_clean_pi;

      //default --hybrid + default topo cuts + daughter pT>0.3
      if (mtopoCut && mPtCut)
      {
         h2PtCut->Fill(rcpt, centrality, weight * dedxEff * tof_hybrid);
         h3PtCentYCut->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_hybrid);
      }

      //clean pid
      if (mtopoCut && mPtCut)
      {
         h2PtCut_clean->Fill(rcpt, centrality, weight * dedxEff * tof_clean);
         h3PtCentYCut_clean->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_clean);
      }

      //daughter pT>0.3
      if (mtopoCut && mPtCut1)
      {
         h2PtCut_pt1->Fill(rcpt, centrality, weight * dedxEff * tof_hybrid);
         h3PtCentYCut_pt1->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_hybrid);
      }

      //daughter pT>0.5
      if (mtopoCut && mPtCut2)
      {
         h2PtCut_pt2->Fill(rcpt, centrality, weight * dedxEff * tof_hybrid);
         h3PtCentYCut_pt2->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_hybrid);
      }

      //topo cut1 -- tight
      if (mtopoCut1 && mPtCut)
      {
         h2PtCut_topoCut1->Fill(rcpt, centrality, weight * dedxEff * tof_hybrid);
         h3PtCentYCut_topoCut1->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_hybrid);
      }

      //topo cut2 -- loose
      if (mtopoCut2 && mPtCut)
      {
         h2PtCut_topoCut2->Fill(rcpt, centrality, weight * dedxEff * tof_hybrid);
         h3PtCentYCut_topoCut2->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_hybrid);
      }

      //topo cut1 -- tight pt1 0.3
      if (mtopoCut1 && mPtCut1)
      {
         h2PtCut_pt1topoCut1->Fill(rcpt, centrality, weight * dedxEff * tof_hybrid);
         h3PtCentYCut_pt1topoCut1->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_hybrid);
      }

      //topo cut2 -- loose pt1 0.3
      if (mtopoCut2 && mPtCut1)
      {
         h2PtCut_pt1topoCut2->Fill(rcpt, centrality, weight * dedxEff * tof_hybrid);
         h3PtCentYCut_pt1topoCut2->Fill(rcpt, centrality, t->rY, weight * dedxEff * tof_hybrid);
      }

   }
   cout << "before cuts, \tnum = " << num << endl;
   cout << "after track cuts, \tnum = " << num1 << endl;
   cout << "after Embedding eff, \tnum = " << num2 << endl;
   cout << "after Hft ratio, \tnum = " << num3 << endl;
   cout << "not find cut pT bin\tnum = " << num4 << endl;

   //write histograph
   WriteHist(fOut);
   cout << "Write histo completed !!! " << endl;

   //deleta histograph
   if (fweight)
   {
      delete fweight;
      fweight = NULL;
   }
   for (int i = 0; i < 9; i++)
   {
      if (htof_kp[i])
      {
         delete htof_kp[i];
         htof_kp[i] = NULL;
      }
      if (htof_pip[i])
      {
         delete htof_pip[i];
         htof_pip[i] = NULL;
      }
      if (htof_km[i])
      {
         delete htof_km[i];
         htof_km[i] = NULL;
      }
      if (htof_pim[i])
      {
         delete htof_pim[i];
         htof_pim[i] = NULL;
      }
   }
   DeleteHist();

   fOut->Close();
}

int getD0PtIndex(float const pt)
{
   int bin = -1;
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((pt >= anaCuts::PtEdge[i]) && (pt < anaCuts::PtEdge[i + 1]))
         bin = i;
   }
   return bin;
}
int getD0CentIndex(float const cent)
{
   int bin = -1;
   for (int i = 0; i < anaCuts::nCent; i++)
   {
      if ((cent >= anaCuts::CentEdge[i]) && (cent < anaCuts::CentEdge[i + 1]))
         bin = i;
   }
   return bin;
}
void loadInputFile()
{
   cout << "Load D0 pT weight..." << endl;
   // TFile* finWeight = new TFile("AuAu010_weight.root");
   // fweight = (TF1*)finWeight->Get("f1Levy010");
   TFile* finWeight = new TFile("D0_Spectra_Run14HFT.root");
   // fweight = (TF1*)finWeight->Get("flevy_0_10");
   fweight = (TF1*)finWeight->Get("flevy_time_pt_0_10");
   finWeight->Close();
   cout << "End of load D0 weight..." << endl;

   cout << "Load k/pi tof match efficiency..." << endl;
   TFile* fTofMatch = new TFile("tofMatch_fit_Run14_17Jan19.root");
   for (int icent = 0; icent < 9; icent++)
   {
      htof_pip[icent] = (TF1*)fTofMatch->Get(Form("funpip_%i", icent));
      htof_kp[icent] = (TF1*)fTofMatch->Get(Form("funkp_%i", icent));
      htof_pim[icent] = (TF1*)fTofMatch->Get(Form("funpim_%i", icent));
      htof_km[icent] = (TF1*)fTofMatch->Get(Form("funkm_%i", icent));
   }
   fTofMatch->Close();

   cout << "Load k/pi tof/tpc pid efficiency..." << endl;
   TFile* fPIDeff1 = new TFile("pion_PidEff_Ks_170822_use.root");
   gTpcPID_pi = (TF1*)fPIDeff1->Get("fpionNsig_eff");
   gTofPID_pi = (TF1*)fPIDeff1->Get("fpionNsigTof_eff");
   fPIDeff1->Close();
   TFile* fPIDeff2 = new TFile("kaon_PidEff_phi_170822_use.root");
   gTpcPID_k = (TF1*)fPIDeff2->Get("fkaonNsig_eff");
   gTofPID_k = (TF1*)fPIDeff2->Get("fkaonNsigTof_eff");
   fPIDeff2->Close();

   cout << "End of load input..." << endl;
}
