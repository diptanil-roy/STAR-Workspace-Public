#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TObject.h"
#include "TClonesArray.h"
// #include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>
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
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include <algorithm>
#include <TString.h>

#pragma link C++ class vector<int> +;

using namespace std;

double M_PION_PLUS = 0.139570;
double M_KAON_PLUS = 0.493677;

Double_t fline(Double_t *x, Double_t *par)
{
    if (x[0] > 1.8 && x[0] < 1.92) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t gline(Double_t *x, Double_t *par)
{
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}


vector<TH1D *> Method(TString TreeName = "D0JetsUSTree"){

	TFile *inputfile = new TFile("ppJetTreeJP.root");
	TDirectory *QADir = (TDirectory *)inputfile->Get("JetTree");
	TTree *TrackTree = (TTree *)QADir->Get(TreeName.Data());

	// Declaration of leaf types
	Float_t         D0Mass;
	Float_t         D0Pt;
	Float_t         D0Eta;
	Float_t         D0Phi;
	Float_t         PionPt;
	Float_t         PionEta;
	Float_t         PionPhi;
	Float_t         PionDCA;
	Float_t         PionNHitsFit;
	Float_t         PionNSigmaPion;
	Float_t         PionNSigmaKaon;
	Float_t         PionTofBeta;
	Float_t         PionTofYLocal;
	Float_t         KaonPt;
	Float_t         KaonEta;
	Float_t         KaonPhi;
	Float_t         KaonDCA;
	Float_t         KaonNHitsFit;
	Float_t         KaonNSigmaPion;
	Float_t         KaonNSigmaKaon;
	Float_t         KaonTofBeta;
	Float_t         KaonTofYLocal;

	Float_t         JetPt;
	Float_t			JetEta;
	Float_t			JetPhi;
	vector <unsigned int> *Triggers = new vector <unsigned int>;


	TrackTree->SetBranchAddress("D0Mass", &D0Mass);
	TrackTree->SetBranchAddress("D0Pt", &D0Pt);
	TrackTree->SetBranchAddress("D0Eta", &D0Eta);
	TrackTree->SetBranchAddress("D0Phi", &D0Phi);
	TrackTree->SetBranchAddress("PionPt", &PionPt);
	TrackTree->SetBranchAddress("PionEta", &PionEta);
	TrackTree->SetBranchAddress("PionPhi", &PionPhi);
	TrackTree->SetBranchAddress("PionDCA", &PionDCA);
	TrackTree->SetBranchAddress("PionNHitsFit", &PionNHitsFit);
	TrackTree->SetBranchAddress("PionNSigmaPion", &PionNSigmaPion);
	TrackTree->SetBranchAddress("PionNSigmaKaon", &PionNSigmaKaon);
	TrackTree->SetBranchAddress("PionTofBeta", &PionTofBeta);
	// TrackTree->SetBranchAddress("PionTofYLocal", &PionTofYLocal);
	if (TreeName.CompareTo("D0JetsUSTree") == 0) TrackTree->SetBranchAddress("PionTofYLocal", &PionTofYLocal);
	else TrackTree->SetBranchAddress("PionTofYLoccal", &PionTofYLocal);
	TrackTree->SetBranchAddress("KaonPt", &KaonPt);
	TrackTree->SetBranchAddress("KaonEta", &KaonEta);
	TrackTree->SetBranchAddress("KaonPhi", &KaonPhi);
	TrackTree->SetBranchAddress("KaonDCA", &KaonDCA);
	TrackTree->SetBranchAddress("KaonNHitsFit", &KaonNHitsFit);
	TrackTree->SetBranchAddress("KaonNSigmaPion", &KaonNSigmaPion);
	TrackTree->SetBranchAddress("KaonNSigmaKaon", &KaonNSigmaKaon);
	TrackTree->SetBranchAddress("KaonTofBeta", &KaonTofBeta);
	TrackTree->SetBranchAddress("KaonTofYLocal", &KaonTofYLocal);

	TrackTree->SetBranchAddress("JetPt", &JetPt);
	TrackTree->SetBranchAddress("JetEta", &JetEta);
	TrackTree->SetBranchAddress("JetPhi", &JetPhi);
	TrackTree->SetBranchAddress("Triggers", &Triggers);

	int nentries = TrackTree->GetEntriesFast();

	TH1D *hD0Mass = new TH1D("hD0Mass", "hD0Mass", 2500, 1.7, 2.05);

	const int nbins = 4;
	double ptbins[nbins + 1] = {0, 1, 3, 5, 10};
	TH1D *hD0Pt = new TH1D("hD0Pt", "hD0Pt", nbins, ptbins);
	TH1D *hD0PtSB = new TH1D("hD0PtSB", "hD0PtSB", nbins, ptbins);

	const int njetbins = 5;
	double ptjetbins[njetbins + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.};
	TH1D *hJetPt = new TH1D("hJetPt", "hJetPt", njetbins, ptjetbins);
	TH1D *hJetPtSB = new TH1D("hJetPtSB", "hJetPtSB", njetbins, ptjetbins);

	for (int event = 0; event < nentries;  event++){

  		if (event%500000==0)cout << event << " Events done." << endl;

   		TrackTree->GetEntry(event);

   		if (D0Mass < 1.7 || D0Mass > 2.1) continue;

   		if (D0Pt < 1.0 || D0Pt > 5.0) continue;
   		if (abs(D0Eta) > 1.0) continue;

   		if (abs(JetEta) > 0.6) continue;
   		if (JetPt < 3.) continue;

   		if (find(Triggers->begin(), Triggers->end(), 370501) != Triggers->end()) continue;
   		if (find(Triggers->begin(), Triggers->end(), 370541) != Triggers->end()) continue;
   		if (find(Triggers->begin(), Triggers->end(), 370542) != Triggers->end()) continue;

   		if (find(Triggers->begin(), Triggers->end(), 370601) != Triggers->end()) continue;
   		if (find(Triggers->begin(), Triggers->end(), 370611) != Triggers->end()) continue;
   		if (find(Triggers->begin(), Triggers->end(), 370621) != Triggers->end()) continue;
   		if (find(Triggers->begin(), Triggers->end(), 370641) != Triggers->end()) continue;

   		TVector3 D0, Jet;
   		D0.SetPtEtaPhi(D0Pt, D0Eta, D0Phi);
   		Jet.SetPtEtaPhi(JetPt, JetEta, JetPhi);

   		double z = (D0.X()*Jet.X() + D0.Y()*Jet.Y())/(Jet.Perp()*Jet.Perp());

   		double pionp = PionPt*TMath::CosH(PionEta);
   		double kaonp = KaonPt*TMath::CosH(KaonEta);

   		float oneOverBetaExpectedpion = sqrt(M_PION_PLUS*M_PION_PLUS / pionp / pionp + 1);
   		float oneOverBetaExpectedkaon = sqrt(M_KAON_PLUS*M_KAON_PLUS / kaonp / kaonp + 1);

   		// if (abs(PionNSigmaPion) > 2) continue;
	   	// if (abs(KaonNSigmaKaon) > 2) continue;
	   	// if (abs(PionNSigmaPion) > abs(PionNSigmaKaon)) continue;
	   	// if (abs(KaonNSigmaKaon) > abs(KaonNSigmaPion)) continue;
	   	// if (abs(KaonNSigmaPion) < 2) continue;
	   	// if (abs(PionNSigmaKaon) < 2) continue;
	   	if (PionTofBeta < -1.9 || PionTofBeta > 2.1) continue;
	   	if (abs(PionTofYLocal) > 1.8) continue;
	   	if (abs(KaonTofYLocal) > 1.8) continue;
	   	if (PionNHitsFit < 15) continue;
	   	if (KaonNHitsFit < 15) continue;

	   	hD0Mass->Fill(D0Mass);	   	

	   	if ((D0Mass >= 1.7 && D0Mass <= 1.78) || (D0Mass >= 1.94 && D0Mass <= 2.02)) {
	   		hD0PtSB->Fill(D0Pt);
	   		hJetPtSB->Fill(z);
	   	}
	   	else if (D0Mass >= 1.82 && D0Mass <= 1.9) {
	   		hD0Pt->Fill(D0Pt);
	   		hJetPt->Fill(z);
	   	}
   	}

   	TH1D *h = (TH1D *)hD0Mass->Clone();
   	TH1D *g = (TH1D *)hD0Pt->Clone();
   	TH1D *k = (TH1D *)hJetPt->Clone();
   	TH1D *l = (TH1D *)hD0PtSB->Clone();
   	TH1D *m = (TH1D *)hJetPtSB->Clone();

   	delete hD0Mass;
   	delete hD0Pt;
   	delete hJetPt;
   	delete hD0PtSB;
   	delete hJetPtSB;

   	vector<TH1D *>s;
   	s.push_back(h);
   	s.push_back(g);
   	s.push_back(k);
   	s.push_back(l);
   	s.push_back(m);

   	return s;
}

void SaveTheD0MassHistograms(){
	TFile *out = new TFile("out.root", "RECREATE");

	vector<TH1D *>s1 = Method("D0JetsUSTree");
    vector<TH1D *>s2 = Method("D0JetsLSTree");

    out->cd();

    TH1D *US = (TH1D *)s1[0]->Clone("US");
	TH1D *LS = (TH1D *)s2[0]->Clone("LS");

	TH1D *USPt = (TH1D *)s1[1]->Clone("USPt");
	TH1D *LSPt = (TH1D *)s2[1]->Clone("LSPt");

	TH1D *USJetPt = (TH1D *)s1[2]->Clone("USJetPt");
	TH1D *LSJetPt = (TH1D *)s2[2]->Clone("LSJetPt");

	TH1D *USPtSB = (TH1D *)s1[3]->Clone("USPtSB");
	TH1D *LSPtSB = (TH1D *)s2[3]->Clone("LSPtSB");

	TH1D *USJetPtSB = (TH1D *)s1[4]->Clone("USJetPtSB");
	TH1D *LSJetPtSB = (TH1D *)s2[4]->Clone("LSJetPtSB");

    US->Write();
    LS->Write();
    USPt->Write();
    LSPt->Write();
    USJetPt->Write();
    LSJetPt->Write();
    USPtSB->Write();
    LSPtSB->Write();
    USJetPtSB->Write();
    LSJetPtSB->Write();
    out->Close();
}

pair<TH1D *, TH1D *> ReadTheD0MassHistograms(){
	TFile *out = new TFile("out.root");
	out->cd();

	TH1D *US = (TH1D *)out->Get("US");
	TH1D *LS = (TH1D *)out->Get("LS");

	pair<TH1D *, TH1D *>s(US, LS);

	return s;
}

pair<TH1D *, TH1D *> ReadTheD0PtHistograms(){
	TFile *out = new TFile("out.root");
	out->cd();

	TH1D *US = (TH1D *)out->Get("USPt");
	TH1D *USSB = (TH1D *)out->Get("USPtSB");

	pair<TH1D *, TH1D *>s(US, USSB);

	return s;
}

pair<TH1D *, TH1D *> ReadTheJetPtHistograms(){
	TFile *out = new TFile("out.root");
	out->cd();

	TH1D *US = (TH1D *)out->Get("USJetPt");
	TH1D *USSB = (TH1D *)out->Get("USJetPtSB");

	pair<TH1D *, TH1D *>s(US, USSB);

	return s;
}


void FitMethod(){

	TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // gStyle->SetOptFit(111111);

    // pair<TH1D *, TH1D *>s1 = Method("D0USTree");
    // pair<TH1D *, TH1D *>s2 = Method("D0LSTree");

    pair<TH1D *, TH1D *>s = ReadTheD0MassHistograms();
    pair<TH1D *, TH1D *>p = ReadTheD0PtHistograms();
    pair<TH1D *, TH1D *>j = ReadTheJetPtHistograms();

	TH1D *h1 = (TH1D *)s.first;
	TH1D *h2 = (TH1D *)s.second;
	
	TH1D *h1pt = (TH1D *)p.first;
	TH1D *h2pt = (TH1D *)p.second;

	TH1D *h1Jetpt = (TH1D *)j.first;
	TH1D *h2Jetpt = (TH1D *)j.second;

	int lowbin = h2->FindBin(2.0);
	int highbin = h2->FindBin(2.1);

	double usscalefactor = h1->Integral(lowbin, highbin);
	double lsscalefactor = h2->Integral(lowbin, highbin);

	int bin1 = h1->FindBin(1.7);
	int bin2 = h1->FindBin(1.78);
	int bin3 = h1->FindBin(1.82);
	int bin4 = h1->FindBin(1.9);
	int bin5 = h1->FindBin(1.94);
	int bin6 = h1->FindBin(2.02);
	// lowbin = h2->FindBin(1.70);
	// highbin = h2->FindBin(1.80);

	cout << h1->Integral(bin1, bin2) << "\t" << h1->Integral(bin3, bin4) << "\t" << h1->Integral(bin5, bin6) << endl;

	cout << "Fraction = " << h1->Integral(bin3, bin4)/(h1->Integral(bin1, bin2) + h1->Integral(bin5, bin6)) << endl;

	double fraction = h1->Integral(bin3, bin4)/(h1->Integral(bin1, bin2) + h1->Integral(bin5, bin6));

	// usscalefactor += h1->Integral(lowbin, highbin);
	// lsscalefactor += h2->Integral(lowbin, highbin);

	cout << usscalefactor << "\t" << lsscalefactor << endl;

	h2->Scale(usscalefactor/lsscalefactor);
	

	h1->Rebin(100);
	h2->Rebin(100);

	TCanvas *b = new TCanvas("b", "b", 700, 500);
	b->cd();

	TH1D *US = (TH1D *)h1->Clone("US");
	TH1D *LS = (TH1D *)h2->Clone("LS");

	US->GetXaxis()->SetTitle("m_{K#pi} [GeV/#it{c}^{2}]");
	US->GetYaxis()->SetTitle("Counts");

	TF1 *fl = new TF1("fl",fline,1.7, 2.1, 2);
   	fl->SetParameters(2,-1, -1);  

	US->Fit("fl", "0L");

	TF1 *gl = new TF1("gl",gline,1.7, 2.1, 2);
   	gl->SetParameters(fl->GetParameters());
   	US->GetListOfFunctions()->Add(gl);

   	TH1D *USSubbed = (TH1D *)US->Clone("USSubbed");

   	for (int i = 1; i <= USSubbed->GetNbinsX(); i++){
   		USSubbed->SetBinContent(i, USSubbed->GetBinContent(i) - gl->Eval(USSubbed->GetBinCenter(i)));
   		USSubbed->SetBinError(i, USSubbed->GetBinError(i));
   	}  


   	TF1 *fitgaus = new TF1("fitgaus", "gaus + pol2(3)", 1.7, 2.1);
   	fitgaus->SetParLimits(0, 1, 2000);
    fitgaus->SetParameter(1, 1.865);
    fitgaus->SetParLimits(1, 1.863, 1.867);
    fitgaus->SetParLimits(2, 0.00001, 0.1);
    fitgaus->SetParameter(3, fl->GetParameter(0));
    fitgaus->SetParameter(4, fl->GetParameter(1));
    fitgaus->SetParameter(5, fl->GetParameter(2));
    US->Fit("fitgaus", "L");



	// US->Draw("EP");
   	US->Draw("EP");

   	TPaveText *ppt = new TPaveText(0.65,0.65,0.85,0.85, "NDC"); // NDC sets coords
    ppt->SetFillColor(0); // text is black on white
    ppt->SetTextSize(0.03); 
    ppt->SetTextAlign(12);
    auto ppt_text0 = ppt->AddText("US");
    auto ppt_text1 = ppt->AddText(Form("p_{T, D^{0}} > 1 GeV/#it{c}"));
    auto ppt_text2 = ppt->AddText(Form("#D^{0} = %.1f #pm %.1f", fitgaus->GetParameter(0), fitgaus->GetParError(0)));
   	ppt->Draw("SAME");

   	cout << "Integrals = " << h1pt->Integral() << "\t" << h2pt->Integral() << "\t" << h1pt->Integral()/h2pt->Integral() << endl;

   	double background = 0.; 
   	double signalplusbackground = 0.;

   	bin1 = US->FindBin(1.7);
	bin2 = US->FindBin(1.78);
	bin3 = US->FindBin(1.82);
	bin4 = US->FindBin(1.9);
	bin5 = US->FindBin(1.94);
	bin6 = US->FindBin(2.02);

	TF1 *backgroundfunction = new TF1("background", "pol2", 1.7, 2.1);
	backgroundfunction->SetParameter(0, fitgaus->GetParameter(3));
	backgroundfunction->SetParameter(1, fitgaus->GetParameter(4));
	backgroundfunction->SetParameter(2, fitgaus->GetParameter(5));

   	for (int i = bin1; i <= bin2; i++){
   		background += backgroundfunction->Eval(US->GetBinCenter(i));
   	}

   	for (int i = bin5; i <= bin6; i++){
   		background += backgroundfunction->Eval(US->GetBinCenter(i));
   	}

	for (int i = bin3; i <= bin4; i++){
   		signalplusbackground += backgroundfunction->Eval(US->GetBinCenter(i));
   	}   

   	cout << "Fractions = " << signalplusbackground << "\t" << background << "\t" << signalplusbackground/background << endl;

   	// h2pt->Scale(h1pt->Integral()/h2pt->Integral());
	// h2Jetpt->Scale(h1Jetpt->Integral()/h2Jetpt->Integral());	

	h2pt->Scale(h1pt->Integral()/h2pt->Integral());
	h2Jetpt->Scale(h1Jetpt->Integral()/h2Jetpt->Integral());

	US->SetLineColor(1);
	US->SetMarkerStyle(20);
	US->SetMarkerColor(1);
	LS->SetLineColor(1);
	LS->SetMarkerStyle(24);
	LS->SetMarkerColor(1);

	h1->Add(h2, -1);
	// h1pt->Add(h2pt, -1);

	h1->GetXaxis()->SetTitle("m_{K#pi} [GeV/#it{c}^{2}]");
	h1->GetYaxis()->SetTitle("Counts");

	TCanvas *c = new TCanvas("c", "c", 700, 500);
	c->cd();
	// US->Draw("EP");
	// LS->Draw("EP SAME");
	h1->Draw("SAME");

	TF1 *bffunc = new TF1("bffunc", "pol2", 1.73, 2.05);
	bffunc->SetParameter(3, 2);
    bffunc->SetParameter(4, -1);
    bffunc->SetParameter(5, -1);

    int bfstatus = h1->Fit("bffunc", "L");

    cout << "BF Status = " << bfstatus << endl;
    cout << bffunc->GetParameter(0) << "\t" << bffunc->GetParameter(1) << "\t" << bffunc->GetParameter(2) << endl;

	// TF1 *fitfunc = new TF1("fitfunc", "signal + background(4)", 1.6, 2.1);
    TF1 *fitfunc = new TF1("fitfunc", "gaus + pol2(3)", 1.73, 2.05);
    fitfunc->SetLineColor(2);
    fitfunc->SetLineWidth(5);

    fitfunc->SetParLimits(0, 1, 1500);
    fitfunc->FixParameter(1, 1.865);
    // fitfunc->SetParLimits(1, 1.863, 1.867);
    fitfunc->SetParLimits(2, 0.00001, 0.1);
    fitfunc->SetParameter(3, bffunc->GetParameter(0));
    fitfunc->SetParameter(4, bffunc->GetParameter(1));
    fitfunc->SetParameter(5, bffunc->GetParameter(2));

    int status = h1->Fit("fitfunc", "");

    cout << "Fit Status = " << status << endl;

    double chi2byndf = fitfunc->GetChisquare()/fitfunc->GetNDF();
    double mean = fitfunc->GetParameter(1);
    double sigma = fitfunc->GetParameter(2);

    TPaveText *pt = new TPaveText(0.65,0.65,0.85,0.85, "NDC"); // NDC sets coords
    pt->SetFillColor(0); // text is black on white
    pt->SetTextSize(0.03); 
    pt->SetTextAlign(12);
    auto pt_text0 = pt->AddText("US - Scaled LS");
    auto pt_text1 = pt->AddText(Form("p_{T, D^{0}} > 1 GeV/#it{c}"));
    auto pt_text2 = pt->AddText(Form("#D^{0} = %.1f #pm %.1f", fitfunc->GetParameter(0), fitfunc->GetParError(0)));

    pt->Draw("SAME");

    TCanvas *d = new TCanvas("d", "d", 700, 500);
	d->cd();

	h1pt->GetXaxis()->SetTitle("p_{T, D^{0}} [GeV/#it{c}]");
	h1pt->GetYaxis()->SetTitle("Counts");
	h1pt->SetNameTitle("Unlike Sign", "Unlike Sign");
	h2pt->SetNameTitle("Sideband (Scaled)", "Sideband (Scaled)");

	h1pt->SetLineColor(1);
	h1pt->SetMarkerColor(1);
	h2pt->SetLineColor(2);
	h2pt->SetMarkerColor(2);
	h1pt->Draw();
	h2pt->Draw("SAME");
	gPad->BuildLegend();

	TCanvas *e = new TCanvas("e", "e", 700, 500);
	e->cd();

	h1Jetpt->GetXaxis()->SetTitle("z");
	h1Jetpt->GetYaxis()->SetTitle("Counts");
	h1Jetpt->SetNameTitle("Unlike Sign", "Unlike Sign");
	h2Jetpt->SetNameTitle("Sideband (Scaled)", "Sideband (Scaled)");


	h1Jetpt->SetLineColor(1);
	h1Jetpt->SetMarkerColor(1);
	h2Jetpt->SetLineColor(2);
	h2Jetpt->SetMarkerColor(2);
	h1Jetpt->Draw();
	h2Jetpt->Draw("SAME");
	gPad->BuildLegend();

	TCanvas *f = new TCanvas("f", "f", 700, 500);
	f->cd();

	TH1D *h1ptClone = (TH1D *)h1pt->Clone();
	h1ptClone->Add(h2pt, -1);
	h1ptClone->Draw();

	TCanvas *g = new TCanvas("g", "g", 700, 500);
	g->cd();

	TH1D *h1JetptClone = (TH1D *)h1Jetpt->Clone();
	h1JetptClone->Add(h2Jetpt, -1);
	h1JetptClone->Draw();


}

void D0InvMassInJets(){
	SaveTheD0MassHistograms();
	FitMethod();
}