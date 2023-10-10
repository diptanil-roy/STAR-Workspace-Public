/* *********************************************************************
 *
 *  Analysis code to read Embedding files.
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/
#include "iostream"
#include <string>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"

#include "TPDFManager.h"
#include "StyleUtilities.h"
#include "tracksMc.h"
#include "eventCountMc.h"

using namespace std;

int main(int argc, char **argv)
{
   eventCountMc* evtMc = new eventCountMc();
   tracksMc* trkMc = new tracksMc();

   TFile* fOut = new TFile("Eff_KaonMinus_embedding.root", "recreate");

   const Int_t nPtBins = 8;
   const Double_t ptEdge[nPtBins + 1] = {0, 1, 2, 3, 4, 5, 7, 9, 12};
   Double_t pt[nPtBins];
   Double_t pterr[nPtBins];
   for (int i = 0; i < nPtBins; i++)
   {
      pt[i] = 0.5 * (ptEdge[i] + ptEdge[i + 1]);
      pterr[i] = pt[i] - ptEdge[i];
   }

   const Int_t nEtaBins = 5;
   const Double_t EtaEdge[nEtaBins + 1] = { -1., -0.6, -0.2, 0.2, 0.6, 1.};
   Double_t Eta[nEtaBins] = { -0.8, -0.4, 0, 0.4, 0.8};
   Double_t Etaerr[nEtaBins] = {0.2, 0.2 , 0.2, 0.2, 0.2};

   const Int_t nPhiBins = 5;
   const Double_t PhiEdge[nPhiBins + 1] = { -3.14159, -1.88496, -0.628319, 0.628319, 1.88496, 3.14159};
   Double_t Phi[nPhiBins] = { -2.51328, -1.25664, 0, 1.25664, 2.51328};
   Double_t Phierr[nPhiBins] = {1.25664, 1.25664, 1.25664, 1.25664, 1.25664};

   const int nCent = 9;
   const Double_t CentEdge[nCent + 1] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5 };
   Double_t Cent[nCent] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
   Double_t Centerr[nCent] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
   const char cent[nCent][1024] = {"70-80%", "60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%", "5-10%", "0-5%"};

   //From Mc Data

   TH3D* pTEtaPhiMc[nCent];
   TH3D* pTEtaPhiRc[nCent];
   TH2D* McpTgpTRc[nCent];
   TH2D* pTResRc[nCent];
   for (int ic = 0; ic < nCent; ic++)
   {
      pTEtaPhiMc[ic] = new TH3D(Form("pTEtaPhiMc_%d", ic), Form("cent_%s;McpT;#eta;#phi", cent[ic]), 120, 0, 12, 100, -1.5, 1.5, 100, -3.14, 3.14);
      pTEtaPhiRc[ic] = new TH3D(Form("pTEtaPhiRc_%d", ic), Form("cent_%s;McpT;#eta;#phi", cent[ic]), 120, 0, 12, 100, -1.5, 1.5, 100, -3.14, 3.14);
      McpTgpTRc[ic] = new TH2D(Form("McpTgpTRc_%d", ic), Form("cent_%d;McpT;RcpT", cent[ic]), 120, 0, 12, 120, 0, 12);
      pTResRc[ic] = new TH2D(Form("pTResRc_%d", ic), Form("cent_%s;McpT;(RcpT - McpT)/ McpT", cent[ic]), 120, 0, 12, 200, -5, 5);
   }

   Long64_t nEntriesEvtMc = evtMc->GetEntries();
   cout << "nEntriesEvtMc = " << nEntriesEvtMc << endl;
   Long64_t nEntriesTrkMc = trkMc->GetEntries();
   cout << "nEntriesTrkMc = " << nEntriesTrkMc << endl;

   for (Long64_t i = 0; i < evtMc->GetEntries(); ++i)
   {
      evtMc->GetEntry(i);

      if (i && i % 10000 == 0) cout << static_cast<float>(i) / nEntriesEvtMc << endl;
      // std::cout<<evtMc->runId<<std::endl;
   } // end event looping

   for (Long64_t i = 0; i < trkMc->GetEntries(); ++i)
   {
      trkMc->GetEntry(i);

      if (i && i % 100000 == 0) cout << static_cast<float>(i) / nEntriesTrkMc << endl;
      // std::cout<<evtMc->runId<<std::endl;
      int ic = trkMc->centrality;
      if (trkMc->centrality < 0) continue;
      // if (trkMc->gPt != -999 && trkMc->pt != -999 && abs(trkMc->eta) < 1.0)pTEtaPhiMc[ic]->Fill(trkMc->pt, trkMc->eta, trkMc->phi);
      if (trkMc->pt != -999 && abs(trkMc->eta) < 1.0)pTEtaPhiMc[ic]->Fill(trkMc->pt, trkMc->eta, trkMc->phi);
      if (trkMc->gPt != -999 && trkMc->pt != -999 && abs(trkMc->eta) < 1.0 && abs(trkMc->gEta) < 1.0 && trkMc->nFit >= 20 && trkMc->nCom > 10  && trkMc->dca<1.5) pTEtaPhiRc[ic]->Fill(trkMc->pt, trkMc->eta, trkMc->phi);
      if (trkMc->gPt != -999 && trkMc->pt != -999 && abs(trkMc->eta) < 1.0 && abs(trkMc->gEta) < 1.0 && trkMc->nFit >= 20 && trkMc->nCom > 10) McpTgpTRc[ic]->Fill(trkMc->pt, trkMc->gPt);
      if (trkMc->gPt != -999 && trkMc->pt != -999 && abs(trkMc->eta) < 1.0 && abs(trkMc->gEta) < 1.0 && trkMc->nFit >= 20 && trkMc->nCom > 10) pTResRc[ic]->Fill(trkMc->pt, (trkMc->gPt - trkMc->pt) / trkMc->pt * 1.e0);
   } // end event looping

   setGraphicsStyle();

   // setStyle(h2);
   // setStyle(h1, 20, 4);

   // gStyle->SetOptTitle(0);
   gStyle->SetOptStat(1);
   TPDFManager* pdf = new TPDFManager("AuAu200GeV.KaonMinus.Embedding.Efficiency");

   TH1D *hsliceMean[nCent], *hsliceSigma[nCent];
   // TH1D* h1mean = new TH1D("h1mean", "Mc;pT;mean", nPtBins, ptEdge);
   char buf[1024];

   for (int ic = 0; ic < nCent; ++ic)
   {
      pdf->newPage(2, 2, Form("Rc pT vs Mc pT, cent_%s", cent[ic])); // divide the page into 2x1 canvases
      pdf->draw(McpTgpTRc[ic], "colz", false , false, false , true);
      pdf->draw(pTResRc[ic], "colz", false , false, false , true);
      pTResRc[ic]->FitSlicesY();
      hsliceMean[ic] = (TH1D*)gDirectory->Get(Form("pTResRc_%d_1", ic));
      hsliceSigma[ic] = (TH1D*)gDirectory->Get(Form("pTResRc_%d_2", ic));
      setStyle(hsliceMean[ic], 20, 1);
      setStyle(hsliceSigma[ic], 20, 1);
      hsliceMean[ic]->GetYaxis()->SetRangeUser(-0.1, 0.1);
      hsliceMean[ic]->GetYaxis()->SetTitle("Mean");
      pdf->draw(hsliceMean[ic], "", false , false, false , true);
      hsliceSigma[ic]->GetYaxis()->SetRangeUser(0, 0.2);
      hsliceSigma[ic]->GetYaxis()->SetTitle("Sigma");
      hsliceSigma[ic]->SetName(Form("hsliceSigma_%d", ic));
      hsliceSigma[ic]->SetTitle(Form("sigma_cent%s", cent[ic]));
      pdf->draw(hsliceSigma[ic], "", false , false, false , true);
   }

   gStyle->SetOptStat(0);

   TH1D *h1effcent[nCent];
   TH1D *h2effcent[nCent];
   TH1D *h1Ratiocent[nCent];
   TGraphAsymmErrors *g1Ratiocent[nCent];
   const Int_t colorIndex[nCent] = {kBlack, kRed, kBlue, kGreen, kPink, kOrange, kTeal, kViolet,  kAzure};

   pdf->newPage(3, 3, Form("Eff Vs pT :( abs(trkMc->gEta) < 1.0 && trkMc->nFit >= 20 && trkMc->nCom > 10 && dca<1.5), diff cent")); // divide the page into 2x1 canvases
   for (int ic = 0; ic < nCent; ic++)
   {
      h1effcent[ic] = (TH1D*)pTEtaPhiRc[ic]->ProjectionX(Form("h1effcent_%d", ic));
      h1effcent[ic]->SetName(Form("h1effcent_%d", ic));
      h1effcent[ic]->SetTitle(Form("h1effcent_%s", cent[ic]));
      // h1effcent[ic]->Sumw2();
      h2effcent[ic] = (TH1D*)pTEtaPhiMc[ic]->ProjectionX(Form("h2effcent_%d", ic));
      h2effcent[ic]->SetName(Form("h2effcent_%d", ic));
      h2effcent[ic]->SetTitle(Form("h2effcent_%s", cent[ic]));
      // h2effcent[ic]->Sumw2();
      // h1effcent[ic]->Scale(1./h1effcent[ic]->GetEntries());
      g1Ratiocent[ic] = new TGraphAsymmErrors(h1effcent[ic], h2effcent[ic]);
      g1Ratiocent[ic]->SetName(Form("g1Ratiocent_%d", ic));
      g1Ratiocent[ic]->SetTitle(Form("g1 efficeency_%s", cent[ic]));
      g1Ratiocent[ic]->SetLineColor(colorIndex[ic]);
      g1Ratiocent[ic]->GetYaxis()->SetRangeUser(0, 1.2);
      h1Ratiocent[ic] = (TH1D*)h1effcent[ic]->Clone(Form("h1Ratiocent_%d", ic));
      h1Ratiocent[ic]->SetName(Form("h1Ratiocent_%d", ic));
      h1Ratiocent[ic]->SetTitle(Form("h1 efficeency_%s", cent[ic]));
      h1Ratiocent[ic]->Divide(h2effcent[ic]);
      h1Ratiocent[ic]->SetLineColor(colorIndex[ic]);
      h1Ratiocent[ic]->GetYaxis()->SetRangeUser(0, 1.2);
      pdf->newLegend(Form("cent_%s", cent[ic]), 0.5, 0.3, 0.7, 0.5); // create a new legend
      // pdf->draw(h1Ratiocent[ic], "", true, false, false, true);
      pdf->draw(g1Ratiocent[ic], "", true, false, false, true);
   }

   TH1D *h1effeta[nCent][nEtaBins];
   TH1D *h2effeta[nCent][nEtaBins];
   TH1D *h1Ratioeta[nCent][nEtaBins];
   TGraphAsymmErrors *g1Ratioeta[nCent][nEtaBins];

   pdf->newPage(3, 3, Form("Eff Vs pT :( abs(trkMc->gEta) < 1.0 && trkMc->nFit >= 20 && trkMc->nCom > 10 && dca<1.5), diff cent,diff Eta")); // divide the page into 2x1 canvases
   for (int ic = 0; ic < nCent; ic++)
   {
      vector<TH1*> h1Ratioetazvector;
      h1Ratioetazvector.clear();
      vector<TGraphAsymmErrors*> g1Ratioetazvector;
      g1Ratioetazvector.clear();
      for (int ieta = 0; ieta < nEtaBins; ieta++)
      {
         int ibin1 = pTEtaPhiRc[ic]->GetYaxis()->FindBin(EtaEdge[ieta] + 1e-6);
         int ibin2 = pTEtaPhiRc[ic]->GetYaxis()->FindBin(EtaEdge[ieta + 1] - 1e-6);
         int ibin3 = pTEtaPhiMc[ic]->GetYaxis()->FindBin(EtaEdge[ieta] + 1e-6);
         int ibin4 = pTEtaPhiMc[ic]->GetYaxis()->FindBin(EtaEdge[ieta + 1] - 1e-6);
         if (ibin1 != ibin3 || ibin2 != ibin4) std::cout << "something  bins width Error" << std::endl;

         h1effeta[ic][ieta] = (TH1D*)pTEtaPhiRc[ic]->ProjectionX(Form("h1effeta_%d_%d", ic, ieta), ibin1, ibin2, 0, -1);
         h1effeta[ic][ieta]->SetName(Form("h1effeta_%d_%d", ic, ieta));
         h1effeta[ic][ieta]->SetTitle(Form("h1effeta_%s_%d", cent[ic], Eta[ieta]));
         // h1effeta[ic][ieta]->Sumw2();
         h2effeta[ic][ieta] = (TH1D*)pTEtaPhiMc[ic]->ProjectionX(Form("h2effeta_%d_%d", ic, ieta), ibin1, ibin2, 0, -1);
         h2effeta[ic][ieta]->SetName(Form("h2effeta_%d_%d", ic, ieta));
         h2effeta[ic][ieta]->SetTitle(Form("h2effeta_%s_%d", cent[ic], Eta[ieta]));
         // h2effeta[ic][ieta]->Sumw2();
         g1Ratioeta[ic][ieta] = new TGraphAsymmErrors(h1effeta[ic][ieta], h2effeta[ic][ieta]);
         g1Ratioeta[ic][ieta]->SetName(Form("g1effeta_%d_%d", ic, ieta));
         g1Ratioeta[ic][ieta]->SetLineColor(colorIndex[ieta]);
         g1Ratioeta[ic][ieta]->GetYaxis()->SetRangeUser(0, 1.2);
         h1Ratioeta[ic][ieta] = (TH1D*)h1effeta[ic][ieta]->Clone(Form("h1effeta_%d_%d", ic, ieta));
         h1Ratioeta[ic][ieta]->SetName(Form("h1effeta_%d_%d", ic, ieta));
         h1Ratioeta[ic][ieta]->SetTitle(Form("h1effeta_cent%s_eta%2.1f", cent[ic], Eta[ieta]));
         h1Ratioeta[ic][ieta]->Divide(h2effeta[ic][ieta]);
         h1Ratioeta[ic][ieta]->SetLineColor(colorIndex[ieta]);
         h1Ratioeta[ic][ieta]->GetYaxis()->SetRangeUser(0, 1.2);
         pdf->newLegend(Form("cent_%s", cent[ic]), 0.5, 0.3, 0.7, 0.5); // create a new legend
         // h1Ratioetazvector.push_back(h1Ratioeta[ic][ieta]);
         g1Ratioetazvector.push_back(g1Ratioeta[ic][ieta]);
      }
      // pdf->draw(h1Ratioetazvector, "", true, false, false, true);
      pdf->draw(g1Ratioetazvector, "", true, false, false, true);
   }

   TH1D *h1effphi[nCent][nPhiBins];
   TH1D *h2effphi[nCent][nPhiBins];
   TH1D *h1Ratiophi[nCent][nPhiBins];
   TGraphAsymmErrors *g1Ratiophi[nCent][nPhiBins];

   pdf->newPage(3, 3, Form("Eff Vs pT :( abs(trkMc->gEta) < 1.0 && trkMc->nFit >= 20 && trkMc->nCom > 10 && dca<1.5), diff cent,diff Phi")); // divide the page into 2x1 canvases
   for (int ic = 0; ic < nCent; ic++)
   {
      vector<TH1*> h1Ratiophizvector;
      h1Ratiophizvector.clear();
      vector<TGraphAsymmErrors*> g1Ratiophizvector;
      g1Ratiophizvector.clear();
      for (int iphi = 0; iphi < nPhiBins; iphi++)
      {
         int ibin1 = pTEtaPhiRc[ic]->GetZaxis()->FindBin(PhiEdge[iphi] + 1e-6);
         int ibin2 = pTEtaPhiRc[ic]->GetZaxis()->FindBin(PhiEdge[iphi + 1] - 1e-6);
         int ibin3 = pTEtaPhiMc[ic]->GetZaxis()->FindBin(PhiEdge[iphi] + 1e-6);
         int ibin4 = pTEtaPhiMc[ic]->GetZaxis()->FindBin(PhiEdge[iphi + 1] - 1e-6);
         if (ibin1 != ibin3 || ibin2 != ibin4) std::cout << "something  bins width Error" << std::endl;

         h1effphi[ic][iphi] = (TH1D*)pTEtaPhiRc[ic]->ProjectionX(Form("h1effphi_%d_%d", ic, iphi), 0, -1, ibin1, ibin2);
         h1effphi[ic][iphi]->SetName(Form("h1effphi_%d_%d", ic, iphi));
         h1effphi[ic][iphi]->SetTitle(Form("h1effphi_%s_%d", cent[ic], Phi[iphi]));
         // h1effphi[ic][iphi]->Sumw2();
         h2effphi[ic][iphi] = (TH1D*)pTEtaPhiMc[ic]->ProjectionX(Form("h2effphi_%d_%d", ic, iphi), 0, -1, ibin1, ibin2);
         h2effphi[ic][iphi]->SetName(Form("h2effphi_%d_%d", ic, iphi));
         h2effphi[ic][iphi]->SetTitle(Form("h2effphi_%s_%d", cent[ic], Phi[iphi]));
         // h2effphi[ic][iphi]->Sumw2();
         g1Ratiophi[ic][iphi] = new TGraphAsymmErrors(h1effphi[ic][iphi], h2effphi[ic][iphi]);
         g1Ratiophi[ic][iphi]->SetName(Form("g1effphi_%d_%d", ic, iphi));
         g1Ratiophi[ic][iphi]->SetLineColor(colorIndex[iphi]);
         g1Ratiophi[ic][iphi]->GetYaxis()->SetRangeUser(0, 1.2);
         h1Ratiophi[ic][iphi] = (TH1D*)h1effphi[ic][iphi]->Clone(Form("h1effphi_%d_%d", ic, iphi));
         h1Ratiophi[ic][iphi]->SetName(Form("h1effphi_%d_%d", ic, iphi));
         h1Ratiophi[ic][iphi]->SetTitle(Form("h1effphi_cent%s_phi%2.1f", cent[ic], Phi[iphi]));
         h1Ratiophi[ic][iphi]->Divide(h2effphi[ic][iphi]);
         h1Ratiophi[ic][iphi]->SetLineColor(colorIndex[iphi]);
         h1Ratiophi[ic][iphi]->GetYaxis()->SetRangeUser(0, 1.2);
         pdf->newLegend(Form("cent_%s", cent[ic]), 0.5, 0.3, 0.7, 0.5); // create a new legend
         h1Ratiophizvector.push_back(h1Ratiophi[ic][iphi]);
         g1Ratiophizvector.push_back(g1Ratiophi[ic][iphi]);
      }
      // pdf->draw(h1Ratiophizvector, "", true, false, false, true);
      pdf->draw(g1Ratiophizvector, "", true, false, false, true);
   }


   pdf->close();

   fOut->cd();
   for (int ic = 0; ic < nCent; ++ic)
   {
      hsliceSigma[ic]->Write();
   }
   for (int ic = 0; ic < nCent; ++ic)
   {
      h1Ratiocent[ic]->Write();
   }
   for (int ic = 0; ic < nCent; ++ic)
   {
      for (int ieta = 0; ieta < nEtaBins; ieta++)
      {
         h1Ratioeta[ic][ieta]->Write();
      }
   }
   for (int ic = 0; ic < nCent; ++ic)
   {
      for (int iphi = 0; iphi < nPhiBins; iphi++)
      {
         h1Ratiophi[ic][iphi]->Write();
      }
   }

   fOut->Close();

   return 0;
}
