#include "../Config.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TVector3.h"
#include <iostream>
#include "TRandom3.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
using namespace std;
Double_t GetDeltaPhi(Double_t phia, Double_t phib);
void MakeResponseNew(int save = 1, int use_weights = 0, int cen_low = 40, int cen_high = 80){


    TFile* bkg_file = new TFile("BKG_Fluc_Correlations.root");
    THnD* ALL_SMEARING;
    TFile* bkg_file1 = new TFile("BKGHistogramsFINAL.root");
    TH1F* smearing;
    
        
    char outfile[200];

    if(use_weights == 0 ) sprintf(outfile,"ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New_Sample_FINAL.root", cen_low , cen_high);
    else sprintf(outfile,"ResponseMatrices_ppWeighted_%i_%i.root", cen_low , cen_high);
    double sigma = 2.7;
    double cutoff = 30.;
    double lowcutoff = 0.;
    double D0PTCUT = 5.;
    double RECOPTCUT = 3.;
    char inmat[200];
    char inmat1[200];
    TFile* FakeRates = new TFile("FakeRates_FINAL.root");
    TH1F *hFakes;
    TH1F *hFakes_2_5;
    TH1F *hFakes_5_10;
    TH1F *hFakes_10_30;
    if(cen_low==40 && cen_high == 80){
	//sprintf(inmat,"/gpfs01/star/pwg/droy1/SimulationFrameworkOfficial/JetsFromSimulation/Response_Peripheral.root");
	//sprintf(inmat,"/gpfs01/star/pwg/droy1/EMBEDDING_2021_D0Analysis/FullPYTHIA/Response4080_BMeson_Nov4/Response4080.root");  
	sprintf(inmat,"SimJetSaverSmear_Peripheral/MCJets");
	sprintf(inmat1,"SimJetSaverSmear_Peripheral/RecoJets");
	cutoff=30;
	sigma = 1.875;//1.62;//2.3;//2.7;
	//RECOPTCUT=10;
	ALL_SMEARING = (THnD*)bkg_file->Get("Smearing_pT_dpT_eta_phi_40_80");
	smearing = (TH1F*)bkg_file1->Get("hBKG_40_80_Weighted");
	hFakes = (TH1F*)FakeRates->Get("ETA_3");
	hFakes_2_5 = (TH1F*)FakeRates->Get("ETA_3_2_5");
	hFakes_5_10 = (TH1F*)FakeRates->Get("ETA_3_5_10");
	hFakes_10_30 = (TH1F*)FakeRates->Get("ETA_3_10_30");
    }
    else if (cen_low ==10 && cen_high ==40){
	//sprintf(inmat,"/gpfs01/star/pwg/droy1/SimulationFrameworkOfficial/JetsFromSimulation/Response_MidCentral.root");
	//sprintf(inmat,"/gpfs01/star/pwg/droy1/EMBEDDING_2021_D0Analysis/FullPYTHIA/Response1040_BMeson_Nov16/Response1040_BMeson.root");
	sprintf(inmat,"SimJetSaverSmear_MidCentral/MCJets");
	sprintf(inmat1,"SimJetSaverSmear_MidCentral/RecoJets");
	cutoff=30;
	sigma =4.4;//5.5;
	//RECOPTCUT = 16;
	smearing = (TH1F*)bkg_file1->Get("hBKG_10_40_Weighted");
	ALL_SMEARING = (THnD*)bkg_file->Get("Smearing_pT_dpT_eta_phi_10_40");
	hFakes = (TH1F*)FakeRates->Get("ETA_2");
	hFakes_2_5 = (TH1F*)FakeRates->Get("ETA_2_2_5");
        hFakes_5_10 = (TH1F*)FakeRates->Get("ETA_2_5_10");
        hFakes_10_30 = (TH1F*)FakeRates->Get("ETA_2_10_30");
    }
    else if (cen_low == 0 && cen_high ==10){
	//sprintf(inmat,"/gpfs01/star/pwg/droy1/SimulationFrameworkOfficial/JetsFromSimulation/Response_Central.root");
	//sprintf(inmat,"/gpfs01/star/pwg/droy1/EMBEDDING_2021_D0Analysis/FullPYTHIA/Response0010_BMeson_Nov16/Response0010_BMeson.root");
	sprintf(inmat,"SimJetSaverSmear_Central/MCJets");
	sprintf(inmat1,"SimJetSaverSmear_Central/RecoJets");
	cutoff=30;
	sigma = 5.8;//6.5;//5.9;//6.6;
	//RECOPTCUT=20;
	smearing = (TH1F*)bkg_file1->Get("hBKG_0_10_Weighted");
	ALL_SMEARING = (THnD*)bkg_file->Get("Smearing_pT_dpT_eta_phi_0_10");
	hFakes = (TH1F*)FakeRates->Get("ETA_1");
	hFakes_2_5 = (TH1F*)FakeRates->Get("ETA_1_2_5");
        hFakes_5_10 = (TH1F*)FakeRates->Get("ETA_1_5_10");
        hFakes_10_30 = (TH1F*)FakeRates->Get("ETA_1_10_30");
    }


    
    cout <<"\n\n > Loading input matrix " << inmat << " \n > Sigma smearing = " << sigma << endl;
    TFile *f_D = new TFile("/gpfs01/star/pwg_tasks/jetcorr03/Mar24_AllUpdatedFiles/Response_Mar24.root","READ");
    gROOT->ProcessLine(".x ~/myStyle.C");
    //gSystem->Load("RooUnfold/libRooUnfold");    
    TChain *reco = (TChain*)f_D->Get(inmat1);//"SimJetSaverSmear/RecoJets");
    TChain *gen = (TChain*)f_D->Get(inmat);//"SimJetSaverSmear/MCJets");

    TH1F *eta;
    TH1F *phi;
    TF1 *eta_ext;
    TF1 *phi_ext;
    
    TFile* eta_file = new TFile("Eta_Res_FINAL.root","READ");
    if(cen_low==40 && cen_high == 80)eta = (TH1F*)eta_file->Get("eta_sigma_pt_40_80");
    else  if (cen_low ==10 && cen_high ==40) eta = (TH1F*)eta_file->Get("eta_sigma_pt_10_40");
    else  if (cen_low ==0 && cen_high ==10) eta = (TH1F*)eta_file->Get("eta_sigma_pt_0_10");
    if(cen_low==40 && cen_high == 80)eta_ext = (TF1*)eta_file->Get("eta_fit_1");
    else  if (cen_low ==10 && cen_high ==40)eta_ext = (TF1*)eta_file->Get("eta_fit_2");
    else  if (cen_low ==40 && cen_high ==80)eta_ext = (TF1*)eta_file->Get("eta_fit_3");
    TFile* phi_file = new TFile("Phi_Res_FINAL.root","READ");
    if(cen_low==40 && cen_high == 80)phi = (TH1F*)phi_file->Get("phi_sigma_pt_40_80");
    else  if (cen_low ==10 && cen_high ==40) phi = (TH1F*)phi_file->Get("phi_sigma_pt_10_40");
    else  if (cen_low ==0 && cen_high ==10) phi = (TH1F*)phi_file->Get("phi_sigma_pt_0_10");
      if(cen_low==40 && cen_high == 80)phi_ext = (TF1*)phi_file->Get("phi_fit_1");
    else  if (cen_low ==10 && cen_high ==40)phi_ext = (TF1*)phi_file->Get("phi_fit_2");
    else  if (cen_low ==40 && cen_high ==80)phi_ext = (TF1*)phi_file->Get("phi_fit_3");



    TH1D *hReco_det = new TH1D("hReco_det ","hReco_det",nbins_jpt,binning_jpt);
    TH1D *hTruth_det = new TH1D("hTruth_det ","hTruth_det",nbins_jpt,binning_jpt);

    
    TH1D *hReco = new TH1D("hReco","hReco",nbins_jpt,binning_jpt);
    TH1D *hTruth= new TH1D("hTruth","hTruth",nbins_jpt,binning_jpt);
    TH1D *hTruthFine= new TH1D("hTruthFine","hTruthFine",61,-0.25,30.25);
    TH1D *hFake= new TH1D("hFake","hFake",nbins_jpt,binning_jpt);
    TH2D *hRes = new TH2D("hRes","hRes",nbins_jpt,binning_jpt,nbins_jpt,binning_jpt);

    TH1D *hReco_dR = new TH1D("hReco_dR","hReco_dR",nBinsdR,dRBins);
    TH1D *hTruth_dR = new TH1D("hTruth_dR","hTruth_dR",nBinsdR,dRBins);
    TH2D *hRes_dR = new TH2D("hRes_dR","hRes_dR",nBinsdR,dRBins,nBinsdR,dRBins);

    TH1D *hTruth_Eta = new TH1D("hTruth_Eta","hTruth_Eta",150,-1.5,1.5);
    TH1D *hD0Pt = new TH1D("hD0Pt","hD0Pt",100,0,10);
    TH1D *hD0Eta = new TH1D("hD0Eta","hD0Eta",100,-1,1);
    TH1D *hJetEta = new TH1D("hJetEta","hJetEta",100,-1,1);
    TH1D *hTrueJetEta = new TH1D("hTrueJetEta","hTrueJetEta",100,-1,1);    
    TH2D *hRes_det = new TH2D("hRes_det","hRes_det",nbins_jpt,binning_jpt,nbins_jpt,binning_jpt);
    TH2D *hRes_nobkg = new TH2D("hRes_nobkg","hRes_nobkg",nbins_jpt,binning_jpt,nbins_jpt,binning_jpt);
    TH1F *hRecoZ = new TH1F("hRecoZ","hRecoZ",nbins_z,binning_z);
    TH2F *hRecoZPt3 = new TH2F("hRecoZPt3","hRecoZPt3",nbins_jpt,binning_jpt,nbins_z,binning_z);

    TH1F *hGenZ = new TH1F("hGenZ","hGenZ",1000,-100,100);
    TH2F *hGenZPt3 = new TH2F("hGenZPt3","hGenZPt3",100,0,50,20,0,1);
    TH1F *hRecoZ3 = new TH1F("hRecoZ3","hRecoZ3",nbins_z,binning_z);
    TH1F *hRecoZ3_nobkg = new TH1F("hRecoZ3_nobkg","hRecoZ3_nobkg",nbins_z,binning_z);
    TH1F *hRecoZ3_1 = new TH1F("hRecoZ3_1","hRecoZ3_1",nbins_z,binning_z);
    const Int_t PTBINS = 11;
    Double_t edges[PTBINS+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    TH1F *hD0 = new TH1F("hD0","hD0",PTBINS,edges);
    TH1F *hD0_1 = new TH1F("hD0_1","hD0_1",PTBINS,edges);
    TH2F *hD0JPt = new TH2F("hD0JPt","hD0JPt",nbins_jpt,binning_jpt,PTBINS,edges);

    /// Weighting histograms


    TH2F *ResPt = new TH2F("ResPt","ResPt",50,0,50,100,-100,100);
    TH2F *g_Hard2D = new TH2F("g_Hard2D","g_Hard2D",100,0,50,100,0,50);
    TH2F *r_Hard2D = new TH2F("r_Hard2D","r_Hard2D",100,0,50,100,0,50);
    TH2F *g_Hard2D2 = new TH2F("g_Hard2D2","g_Hard2D2",100,0,50,100,0,5);
    TH2F *r_Hard2D2 = new TH2F("r_Hard2D2","r_Hard2D2",100,0,50,100,0,5);

    //RooUnfoldResponse* response;
    //RooUnfoldResponse* response_3;
    //response = new RooUnfoldResponse(hReco,hTruth, "response", "response");
    //response_3 = new RooUnfoldResponse(hReco,hTruth, "response_3", "response_3");

    TH2F *signalJetPtvdRD0Jet = new TH2F("signalJetPtvdRD0Jet", "signalJetPtvdRD0Jet", nBinsJetPt, JetPtBins, nBinsdR, dRBins);
    RooUnfoldResponse *response2D;

    response2D = new RooUnfoldResponse(signalJetPtvdRD0Jet, signalJetPtvdRD0Jet, "Response2D", "Response2D");
    
    float r_TrackCharge[1000];
    float r_TrackPt[1000];
    float r_TrackPx[1000];
    float r_TrackPy[1000];
    float r_TrackPz[1000];
    float r_TrackID[1000];
    float r_D0Mass;
    float r_JetPt;
    float r_JetPhi;
    float r_JetEta;
    float r_R;
    float r_KaonPt;
    float r_KaonPhi;
    float r_KaonEta;
    float r_PionPt;
    float r_PionPhi;
    float r_PionEta;
    int r_JetNConst;
    float r_JetHighestTrackPt;
    reco->SetBranchAddress( "JetNConst" , &r_JetNConst );
    reco->SetBranchAddress( "JetHighestTrackPt" , &r_JetHighestTrackPt );
    //reco->SetBranchAddress( "TrackPt" , &r_TrackPt[r_JetNConst] );
    reco->SetBranchAddress( "TrackPx" , r_TrackPx );
    reco->SetBranchAddress( "TrackPy" , r_TrackPy );
    reco->SetBranchAddress( "TrackPz" , r_TrackPz );
    reco->SetBranchAddress( "TrackID" , r_TrackID );
    //reco->SetBranchAddress( "TrackCharge" , &r_TrackCharge[r_JetNConst] );
    reco->SetBranchAddress( "D0Mass" , &r_D0Mass );
    reco->SetBranchAddress( "JetPt" , &r_JetPt );
    reco->SetBranchAddress( "JetPhi" , &r_JetPhi );
    reco->SetBranchAddress( "JetEta" , &r_JetEta );
    reco->SetBranchAddress( "KaonPt" , &r_KaonPt );
    reco->SetBranchAddress( "KaonPhi" , &r_KaonPhi );
    reco->SetBranchAddress( "PionPt" , &r_PionPt );
    reco->SetBranchAddress( "PionPhi" , &r_PionPhi );
    reco->SetBranchAddress( "PionEta" , &r_PionEta );
    reco->SetBranchAddress( "KaonEta" , &r_KaonEta );
    reco->SetBranchAddress( "JetRadius" , &r_R );
    
    float g_D0Mass;
    float g_JetPt;
    float g_JetPhi;
    float g_KaonPt;
    float g_KaonPhi;
    float g_PionPt;
    float g_PionPhi;
    float g_PionEta;
    float g_KaonEta;
    float g_R;
    float g_JetEta;
    int g_JetNConst;
    float g_JetHighestTrackPt;
    float g_TrackCharge[1000];
    float g_TrackPt[1000];
    float g_TrackPx[1000];
    float g_TrackPy[1000];
    float g_TrackPz[1000];
    float g_TrackID[1000];
    gen->SetBranchAddress( "JetNConst" , &g_JetNConst );
    gen->SetBranchAddress( "JetHighestTrackPt" , &g_JetHighestTrackPt );
    gen->SetBranchAddress( "TrackPt" , &g_TrackPt );
    gen->SetBranchAddress( "TrackPx" , &g_TrackPx );
    gen->SetBranchAddress( "TrackPy" , &g_TrackPy );
    gen->SetBranchAddress( "TrackPz" , &g_TrackPz );
    gen->SetBranchAddress( "TrackID" , &g_TrackID );
//gen->SetBranchAddress( "TrackCharge" , &g_TrackCharge[g_JetNConst] );
    gen->SetBranchAddress( "D0Mass" , &g_D0Mass );
    gen->SetBranchAddress( "JetPt" , &g_JetPt );
    gen->SetBranchAddress( "JetPhi" , &g_JetPhi );
    gen->SetBranchAddress( "JetEta" , &g_JetEta );
    gen->SetBranchAddress( "KaonPt" , &g_KaonPt );
    gen->SetBranchAddress( "KaonPhi" , &g_KaonPhi );
    gen->SetBranchAddress( "PionPt" , &g_PionPt );
    gen->SetBranchAddress( "PionPhi" , &g_PionPhi );
    gen->SetBranchAddress( "PionEta" , &g_PionEta );
    gen->SetBranchAddress( "KaonEta" , &g_KaonEta );
    gen->SetBranchAddress( "JetRadius" , &g_R );
    int loop = reco->GetEntries();

    TFile* FONLL = new TFile("FONLLWeights_New.root","READ");
    TH1F * FONLLWeights = (TH1F*) FONLL->Get("FONLLWeights");
    TF1 * FONLLWeights_HighPt = (TF1*) FONLL->Get("FONLLWeights_HighPt");




    //    TFile* fD0Weights = new TFile("D0Weights.root","READ");
    //D0Weights = (TH1F*)fD0Weights->Get("D0Weights");
    
    TH3F* smearing_3d[nbins_jpt];
    for(int i = 0; i< nbins_jpt; i++){
	ALL_SMEARING->GetAxis(0)->SetRangeUser(binning_jpt[i],binning_jpt[i+1]);
	smearing_3d[i] = (TH3F*)ALL_SMEARING->Projection(1,2,3);
	smearing_3d[i]->SetName(Form("smearing_3d_%i",i));
    }
    




    for(int i =0;i<loop;i++){
	if(i%100000==0)cout << "On "<< i << " out of " << loop << " " << float(i)/loop*100 << "%" << endl;
	reco->GetEntry(i);
	gen->GetEntry(i);
	
	//if(g_JetPt<3)continue;
	//if(g_JetPt>cutoff)continue;
	//if(g_JetPt<lowcutoff)continue;
	//Global Cut
	//if(g_JetNConst<=1)continue;
	if(r_JetNConst==0)continue;
	if(r_PionPt<=0.6)continue;
        if(r_KaonPt<=0.6)continue;
        if(r_PionEta>=1.)continue;
        if(r_KaonEta>=1.)continue;
	if(g_PionPt<=0.6)continue;
	if(g_KaonPt<=0.6)continue;
	if(g_PionEta>=1.)continue;
	if(g_KaonEta>=1.)continue;

	TRandom* r = new TRandom(0);
	
	double _dpt;
	double _deta;
	double _dphi;
	int project_bin = hReco->FindBin(r_JetPt);
	double DPT = 1;//r->Gaus(0,sigma);
        //smearing_3d[project_bin-1]->GetXaxis()->SetRange(DPT-5. , DPT+5.);
	//TH2F* eta_phi = (TH2F*)smearing_3d[project_bin-1]->Project3D("zy");
	//eta_phi->GetRandom2(_deta,_dphi);
	//delete eta_phi;
	//smearing_3d[project_bin-1]->GetRandom3(_dpt,_deta,_dphi);
	//double r_jet_pt = r_JetPt + r->Gaus(0,sigma);  //
	double r_jet_pt = r_JetPt +smearing->GetRandom();//gRandom->Gaus(0,sigma);
	int bin_eta = eta->FindBin(r_JetPt);
	int bin_phi = phi->FindBin(r_JetPt);
	double eta_sigma = eta->GetBinContent(bin_eta);
	double phi_sigma = phi->GetBinContent(bin_phi);
	//e ta_sigma = eta->GetBinContent(bin_eta);
	//phi_sigma = phi->GetBinContent(bin_phi);
		  
	double r_jet_eta = r_JetEta + r->Gaus(0,eta_sigma);
	double r_jet_phi = r_JetPhi + r->Gaus(0,phi_sigma);
	if(r_jet_phi>=2.*TMath::Pi()) r_jet_phi = r_jet_phi - 2.*TMath::Pi();
	if(r_jet_phi<0) r_jet_phi = 2.*TMath::Pi() + r_jet_phi;
	double g_pi_px = g_PionPt*cos(g_PionPhi);
	double g_pi_py = g_PionPt*sin(g_PionPhi);
	double g_k_px = g_KaonPt*cos(g_KaonPhi);
	double g_k_py = g_KaonPt*sin(g_KaonPhi);
	double g_k_pz = g_KaonPt * sinh(g_KaonEta);
	double g_pi_pz = g_PionPt * sinh(g_PionEta);
	double g_d0_px = g_pi_px + g_k_px;
	double g_d0_py = g_pi_py + g_k_py;
	double g_d0_pz = g_pi_pz + g_k_pz;
	TVector3 d0(g_d0_px,g_d0_py,0);
	//Global Cut
	if(sqrt(g_d0_px*g_d0_px+g_d0_py*g_d0_py)<D0PTCUT) continue;
	if(sqrt(g_d0_px*g_d0_px+g_d0_py*g_d0_py)>10.)continue;
	if(sqrt(g_d0_px*g_d0_px+g_d0_py*g_d0_py)>7.)hTrueJetEta->Fill(g_JetEta);
	if(g_JetPt<D0PTCUT)continue;
	double g_jet_px = g_JetPt*cos(g_JetPhi);
	double g_jet_py = g_JetPt*sin(g_JetPhi);

	TVector3 p4;
	p4.SetXYZ(g_d0_px,g_d0_py,g_d0_pz);
	/*for(int idx = 0; idx<g_JetNConst;idx++){
	    if(g_TrackID[idx]==421){
		p4.SetXYZ(g_TrackPx[idx],g_TrackPy[idx],g_TrackPz[idx]);
		break;
	    }
	    }*/
	double g_R = sqrt((p4.PseudoRapidity()-g_JetEta)*(p4.PseudoRapidity()-g_JetEta)+GetDeltaPhi(p4.Phi(),g_JetPhi)*GetDeltaPhi(p4.Phi(),g_JetPhi));
	
	double g_z = (g_jet_px*g_d0_px + g_jet_py*g_d0_py)/TMath::Power(g_JetPt,2);


	double WEIGHT = 1;

	if(g_JetPt<15)WEIGHT = FONLLWeights->GetBinContent(FONLLWeights->FindBin(g_JetPt));//d0pt_weight->GetBinContent(d0pt_weight->FindBin(sqrt(g_d0_px*g_d0_px+g_d0_py*g_d0_py)));
	//else WEIGHT = FONLLWeights->GetBinContent(FONLLWeights->FindBin(12));
	else WEIGHT = FONLLWeights_HighPt->Eval(g_JetPt);
 //TMath::Power(1-g_z,4);
	//if(g_JetPt<3.)WEIGHT*=10;
        if(WEIGHT <=0)WEIGHT=1;//Should not happen anyway                                                                                                                                                                                                    
	//if(fabs(g_JetEta)>=0.6)WEIGHT*=10;
	//if(WEIGHT >300.)WEIGHT=300.;                             
	
	//WEIGHT *= D0Weights->GetBinContent(D0Weights->FindBin(sqrt(g_d0_px*g_d0_px+g_d0_py*g_d0_py)));
	
	
	double WW = WEIGHT;//*TMath::Power((1.+ g_z),2);//(1.+ g_z)*(1.+ g_z)*(1.+ g_z)*(1.+ g_z);
	hTruth_Eta->Fill(g_JetEta,WW);
	
	//If D0 is in jet
	if(r_D0Mass>0){
	    double r_pi_px = r_PionPt*cos(r_PionPhi);
	    double r_pi_py = r_PionPt*sin(r_PionPhi);
	    double r_k_px = r_KaonPt*cos(r_KaonPhi);
	    double r_k_py = r_KaonPt*sin(r_KaonPhi);
	    double r_k_pz = r_KaonPt * sinh(r_KaonEta);
	    double r_pi_pz = r_PionPt * sinh(r_PionEta);
	    double r_d0_px = r_pi_px + r_k_px;
	    double r_d0_py = r_pi_py + r_k_py;
	    double r_d0_pz = r_pi_pz + r_k_pz;
	    double r_d0_pt = sqrt(r_d0_px*r_d0_px+r_d0_py*r_d0_py);
	    double r_jet_px = r_jet_pt*cos(r_JetPhi);
	    double r_jet_py = r_jet_pt*sin(r_JetPhi);
	    double R_jet_px = r_JetPt*cos(r_JetPhi);
            double R_jet_py = r_JetPt*sin(r_JetPhi);
	    double r_z = (r_jet_px*r_d0_px + r_jet_py*r_d0_py)/TMath::Power(r_jet_pt,2);
	    
	    hGenZ->Fill(g_z,WW);
	    hGenZPt3->Fill(g_JetPt,g_z,WW);
	    g_Hard2D->Fill(d0.Pt(),g_JetHighestTrackPt);
	    g_Hard2D2->Fill(g_JetPt,g_JetHighestTrackPt/d0.Pt());
	    double R_z = (R_jet_px*r_d0_px + R_jet_py*r_d0_py)/TMath::Power(r_JetPt,2);
	    TVector3 d0r(r_d0_px,r_d0_py,0);
	    TVector3 p3;
	    p3.SetXYZ(r_d0_px,r_d0_py,r_d0_pz);
	    /*for(int idx = 0; idx<r_JetNConst;idx++){
		if(r_TrackID[idx]==421){
		    p3.SetXYZ(r_TrackPx[idx],r_TrackPy[idx],r_TrackPz[idx]);
		    break;
		}
		}*/
	    double r_R = sqrt((p3.PseudoRapidity()-r_JetEta)*(p3.PseudoRapidity()-r_JetEta)+GetDeltaPhi(p3.Phi(),r_JetPhi)*GetDeltaPhi(p3.Phi(),r_JetPhi));
	    double jet_R = sqrt((g_JetEta-r_JetEta)*(g_JetEta-r_JetEta)+GetDeltaPhi(g_JetPhi,r_JetPhi)*GetDeltaPhi(g_JetPhi,r_JetPhi));
	    
	    double r_R_shift = sqrt((p3.PseudoRapidity()-r_jet_eta)*(p3.PseudoRapidity()-r_jet_eta)+GetDeltaPhi(p3.Phi(),r_jet_phi)*GetDeltaPhi(p3.Phi(),r_jet_phi));

	    hRecoZ->Fill(r_z);
	    //response->Fill(r_jet_pt,g_JetPt,ww);
	    double ww = 1;
	    double ww1 = WEIGHT;//*TMath::Power((1.+ g_z),2);//*(1.+ g_z)*(1.+ g_z)*(1.+ g_z);

	    ResPt->Fill(g_JetPt,(g_JetPt-r_JetPt),ww1);
	    

	    int Jet_flag=-1;
	    //Note:jet_R = radius between particle and detector level jets. NO Eta or Phi smearing applied here!
	    //if(fabs(g_JetEta)>=0.6 && fabs(r_jet_eta)>=0.6)continue;
	    int ISFAKE=0;
	    if(fabs(g_JetEta)>=0.6){
		double val_ = r->Uniform(0,1);
		double comp=-1;
		if(g_JetPt<5)comp = hFakes_2_5->GetBinContent(hFakes_2_5->FindBin(fabs(g_JetEta)));
		else if(g_JetPt>=5 && g_JetPt<10)comp = hFakes_5_10->GetBinContent(hFakes_5_10->FindBin(fabs(g_JetEta)));
		else if(g_JetPt>=10 && g_JetPt<300)comp = hFakes_10_30->GetBinContent(hFakes_10_30->FindBin(fabs(g_JetEta)));
		if(val_ > comp)continue;
		else if(val_ <= comp) ISFAKE=1;
	    }

	    if(r_jet_pt > RECOPTCUT && jet_R <0.4 && fabs(g_JetEta)<0.6 && fabs(r_jet_eta)<0.6) Jet_flag = 1; // OK Jets
	    else if(r_jet_pt > RECOPTCUT && fabs(g_JetEta)>=0.6 && ISFAKE==1) Jet_flag = 2; // Fakes -> REMOVED fabs(r_jet_eta)<0.6 && here for New Fake def.
	    else if(fabs(g_JetEta)<0.6 && (r_jet_pt <= RECOPTCUT || jet_R>=0.4 || fabs(r_jet_eta)>=0.6)) Jet_flag = 3; //Misses
	    else continue;
	    //if(r_jet_pt > 3. && g_JetPt > 1.)Jet_flag = 1; // OK Jets  
	    //else if(r_jet_pt <= 3.) Jet_flag = 3; //Misses   
	    
	    if(Jet_flag ==1){
		hTruth_det->Fill(r_JetPt,ww1);
		hReco_det->Fill(r_jet_pt,ww1);
		hRes_det->Fill(r_jet_pt,r_JetPt,ww1);
	    }
	    else if(Jet_flag == 3){ //Misses     
		hTruth_det->Fill(r_JetPt,ww1);
	    }
	    else if(Jet_flag == 2){
		hReco_det->Fill(r_jet_pt,ww1);
	    }
	    
	    

	    if(Jet_flag ==1){// OK Jets
		if(r_JetPt/g_JetPt<1.)hRecoZ3_1->Fill(r_z,ww1);
		hD0Pt->Fill(r_d0_pt,ww1);
		hD0Eta->Fill(p3.PseudoRapidity(),ww1);
		hJetEta->Fill(r_jet_eta,ww1);
		hRecoZ3->Fill(r_z,ww1);
		hRecoZ3_nobkg->Fill(R_z,ww1);
		hRes->Fill(r_jet_pt,g_JetPt,ww1);
		hRes_dR->Fill(r_R_shift,g_R,ww1);
		hRes_nobkg->Fill(r_JetPt,g_JetPt,ww1);
		hReco->Fill(r_jet_pt,ww1);
		hReco_dR->Fill(r_R_shift,ww1);
		hTruth->Fill(g_JetPt,ww1);
		hTruth_dR->Fill(g_R,ww1);
		hTruthFine->Fill(g_JetPt,ww1);
		hD0->Fill(r_d0_pt,ww1);
		hD0JPt->Fill(r_jet_pt,r_d0_pt,ww1);
		hD0_1->Fill(r_d0_pt,ww1);
		hRecoZPt3->Fill(r_jet_pt,r_z,ww1);
		r_Hard2D->Fill(d0r.Pt(),r_JetHighestTrackPt,ww1);
		r_Hard2D2->Fill(r_JetPt,r_JetHighestTrackPt/d0r.Pt(),ww1);
		response2D->Fill(r_jet_pt, r_R_shift, g_JetPt, g_R, ww1);
		//response_3->Fill(r_jet_pt,g_JetPt,ww1);
		if(g_z<0)cout <<"Negative reco z " << r_JetPhi << " " << d0r.Phi() << endl;	
	    }else if(Jet_flag == 3){ //Misses
		hTruth->Fill(g_JetPt,ww1);
		hTruth_dR->Fill(g_R,ww1);
		hTruthFine->Fill(g_JetPt,ww1);
		hGenZ->Fill(g_z,WW);
		hGenZPt3->Fill(g_JetPt,g_z,WW);
		response2D->Miss( g_JetPt, g_R, WW);
		hD0_1->Fill(r_d0_pt,ww1);
		g_Hard2D->Fill(d0.Pt(),g_JetHighestTrackPt);
		g_Hard2D2->Fill(g_JetPt,g_JetHighestTrackPt/d0.Pt());
		//response_3->Miss(g_JetPt,ww1);  
	    }else if(Jet_flag == 2){ //fakes
		hD0Pt->Fill(r_d0_pt,ww1);
		hD0Eta->Fill(p3.PseudoRapidity(),ww1);
		hJetEta->Fill(r_jet_eta,ww1);
		response2D->Fake(r_jet_pt, r_R_shift,ww1);
		hReco_dR->Fill(r_R_shift,ww1);
		hReco->Fill(r_jet_pt,ww1); 
		hD0->Fill(r_d0_pt,ww1);
		hTruthFine->Fill(g_JetPt,ww1);
		hRes_nobkg->Fill(r_JetPt,g_JetPt,ww1);
		hGenZ->Fill(g_z,WW);
                hGenZPt3->Fill(g_JetPt,g_z,WW);
	    }
	}else{
	    hGenZ->Fill(g_z,WW);
	    hGenZPt3->Fill(g_JetPt,g_z,WW);
	    g_Hard2D->Fill(d0.Pt(),g_JetHighestTrackPt);
	    g_Hard2D2->Fill(g_JetPt,g_JetHighestTrackPt/d0.Pt());
	    //hTruth->Fill(g_JetPt);
	    //response->Miss(g_JetPt);
	    //response_3->Miss(g_JetPt);
	}
    }


    TCanvas *c1 = new TCanvas("c1","c1");
    hRecoZ3->Draw();
    if(save){
	TFile *out = new TFile(outfile,"RECREATE");
	out->cd();
	hD0->Write();
	hD0_1->Write();
	hD0JPt->Write();
	hTruth->Write();
	hTruth_dR->Write();
	hReco_dR->Write();
	hRes_dR->Write();
	hTruthFine->Write();
	hReco->Write();
	hTrueJetEta->Write();
	hTruth_det->Write();
        hReco_det->Write();
	hFake->Write();
	hRecoZ->Write();
	hRecoZ3->Write();
	hRecoZ3_nobkg->Write();
	hRecoZ3_1->Write();
	hGenZ->Write();
	//out->WriteTObject(response, "response");
	//out->WriteTObject(response_3, "response_3");
	hD0Pt->Write();
	hD0Eta->Write();
	hJetEta->Write();
	signalJetPtvdRD0Jet->Write();
	out->WriteObject(response2D, "Response2D");
	hRes->Write();
	hRes_det->Write();
	hRes_nobkg->Write();
	ResPt->Write();
	r_Hard2D->Write();
	g_Hard2D->Write();
	r_Hard2D2->Write();
        g_Hard2D2->Write();
	hGenZPt3->Write();
	hRecoZPt3->Write();
	hTruth_Eta->Write();
    }
}
Double_t GetDeltaPhi(Double_t phia, Double_t phib)
{
    Double_t pi = TMath::Pi();
    
    if (phia < 0)         phia += 2*pi;
    else if (phia > 2*pi) phia -= 2*pi;
    if (phib < 0)         phib += 2*pi;
    else if (phib > 2*pi) phib -= 2*pi;
    Double_t dphi = phib - phia;

    if (dphi < -1.0*pi)  dphi += 2*pi;
    else if (dphi > pi)  dphi -= 2*pi;

    // phi is between  -pi < phi < pi                                                                                                                                                                                                \
                                                                                                                                                                                                                                      
    return dphi;
}
void RUN(){
    MakeResponseNew(1,0,40,80);
    MakeResponseNew(1,0,10,40);
    MakeResponseNew(1,0,0,10);
}
