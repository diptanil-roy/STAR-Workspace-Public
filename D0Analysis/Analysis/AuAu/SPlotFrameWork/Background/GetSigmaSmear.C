#include "TFile.h"
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TChain.h"
#include "TF1.h"
#include "THn.h"
using namespace std;
void GetSigmaSmear(){

    //TFile *f_D = new TFile("/gpfs01/star/pwg_tasks/jetcorr03/tmp/Embed_April1.root");
    TFile *f_D = new TFile("/gpfs01/star/pwg_tasks/jetcorr03/tmp/Embed_April2.root");
    TChain* ch = (TChain*) f_D->Get("JetTree_Embed/D0Jets");
    float gpt; float rpt;float  mCen;float jpt;float geta; int JetNConst;
    float Weight;float Centrality;float jeta; float jphi;float gphi;
    ch->SetBranchAddress( "Weight" , &Weight );
    ch->SetBranchAddress( "PionPt" , &gpt );// Both Pion and Kaon variables have the generated PT
    ch->SetBranchAddress( "PionEta" , &geta );
    ch->SetBranchAddress( "PionPhi" , &gphi );
    ch->SetBranchAddress( "JetCorrPt" , &rpt );// jet PT - rho*A
    ch->SetBranchAddress( "Centrality" , &Centrality );//centrality (0,10,20,etc)                                                                                                                                     
    ch->SetBranchAddress( "JetPt" , &jpt );// jet PT - rho*A                                                                                                                                                                   
    ch->SetBranchAddress( "JetEta" , &jeta );
    ch->SetBranchAddress( "JetPhi" , &jphi );
    ch->SetBranchAddress( "JetNConst" , &JetNConst);
    int loop = ch->GetEntries();
    cout <<"> Filling response matrix " << endl;
    TH3D *PT2D  = new TH3D("PT3D","PT3D",400,-40,40,100,-20,20,10,0,100);
    TH3D *Corr3D_phi  = new TH3D("Corr3D_phi","Corr3D_phi",40,0,40,500,-4,4,10,0,100);
    TH3D *Corr3D_eta  = new TH3D("Corr3D_eta","Corr3D_eta",40,0,40,100,-1,1,10,0,100);

    TH2D *Corr2D  = new TH2D("Corr2D","Corr2D",400,-40,40,10,0,100);
    TH1D *Corr1D  = new TH1D("Corr1D","Corr1D",9,0.5,9.5);


    TH2D *Corr_Eta_Phi = new TH2D("Corr_Eta_Phi","Corr_Eta_Phi",100,-1,1,500,-4,4);
    TH3D *Corr_Pt_R = new TH3D("Corr_Pt_R","Corr_Pt_R",100,-20,20,200,-2,2,10,0,100);
    const int nbins_jpt =13;
    double binning_jpt[nbins_jpt+1] = {0,1,2,3,4,5,7,9,11,13,15,20,30,50};

    int nbins[4] ={nbins_jpt,100,100,100};
    double bins1[101];
    double bins2[101];
    double bins3[101];
    double bins4[101];

    //for(int i = 0;i<nbins[0]+1;i++)bins1[i] = i*50./nbins[0];
    for(int i = 0;i<nbins[1]+1;i++)bins2[i] = -40+i*80./nbins[1];
    for(int i = 0;i<nbins[2]+1;i++)bins3[i] = -1.+i*2./nbins[2];
    for(int i = 0;i<nbins[3]+1;i++)bins4[i] = -4.+i*8./nbins[3];

    for(int i = 0;i<nbins[2]+1;i++) cout << bins3[i] << endl;
    

    THnD *Smearing_pT_dpT_eta_phi_0_10 = new THnD("Smearing_pT_dpT_eta_phi_0_10 ","Smearing_pT_dpT_eta_phi_0_10", 4, nbins, NULL, NULL);
    THnD *Smearing_pT_dpT_eta_phi_10_40 = new THnD("Smearing_pT_dpT_eta_phi_40_80 ","Smearing_pT_dpT_eta_phi_40_80",4, nbins, NULL, NULL);
    THnD *Smearing_pT_dpT_eta_phi_40_80 = new THnD("Smearing_pT_dpT_eta_phi_10_40 ","Smearing_pT_dpT_eta_phi_10_40",4, nbins, NULL, NULL);

    Smearing_pT_dpT_eta_phi_0_10->SetBinEdges(0,binning_jpt);
    Smearing_pT_dpT_eta_phi_0_10->SetBinEdges(1,bins2);
    Smearing_pT_dpT_eta_phi_0_10->SetBinEdges(2,bins3);
    Smearing_pT_dpT_eta_phi_0_10->SetBinEdges(3,bins4);

    Smearing_pT_dpT_eta_phi_10_40->SetBinEdges(0,binning_jpt);
    Smearing_pT_dpT_eta_phi_10_40->SetBinEdges(1,bins2);
    Smearing_pT_dpT_eta_phi_10_40->SetBinEdges(2,bins3);
    Smearing_pT_dpT_eta_phi_10_40->SetBinEdges(3,bins4);

    Smearing_pT_dpT_eta_phi_40_80->SetBinEdges(0,binning_jpt);
    Smearing_pT_dpT_eta_phi_40_80->SetBinEdges(1,bins2);
    Smearing_pT_dpT_eta_phi_40_80->SetBinEdges(2,bins3);
    Smearing_pT_dpT_eta_phi_40_80->SetBinEdges(3,bins4);



    for(int i =0;i<loop;i++){
        ch->GetEntry(i);
	if(i%10000 == 0) cout << "on " << i << " of " << loop << endl;
	
	int ccen = -1;
	if(Centrality == 0)ccen = 70;
	if(Centrality == 1)ccen = 60;
	if(Centrality == 2)ccen = 50;
	if(Centrality == 3)ccen = 40;
	if(Centrality == 4)ccen = 30;
	if(Centrality == 5)ccen = 20;
	if(Centrality == 6)ccen = 10;
	if(Centrality == 7)ccen = 5;
	if(Centrality == 8)ccen = 0;
	
//	if(fabs(geta)>0.6)continue;
        if(fabs(jeta)>0.6)continue;
	//if(fabs(geta)>0.2)continue;
	Weight = 1;	
	PT2D->Fill(rpt-gpt,gpt,ccen+0.1,Weight);
        Corr2D->Fill(jpt-rpt,ccen+0.1,Weight);
	Corr3D_phi->Fill(gpt,jphi-gphi,ccen+0.1,Weight);
	Corr3D_eta->Fill(gpt,jeta-geta,ccen+0.1,Weight);
	Corr_Eta_Phi->Fill(jeta-geta,jphi-gphi);
	Corr_Pt_R->Fill(rpt-gpt,sqrt((geta-jeta)*(geta-jeta)+(gphi-jphi)*(gphi-jphi)),ccen);
	double tofill[4] = {gpt,rpt-gpt,jeta-geta,jphi-gphi};
	if(ccen<10)Smearing_pT_dpT_eta_phi_0_10->Fill(tofill);
	if(ccen>=10 && ccen< 40)Smearing_pT_dpT_eta_phi_10_40->Fill(tofill);
	if(ccen>=40)Smearing_pT_dpT_eta_phi_40_80->Fill(tofill);
    }
	
    TFile* outfile = new TFile("BKG_Fluc_CorrelationsNew.root","RECREATE");
    PT2D->Write();
    Corr3D_phi->Write();
    Corr3D_eta->Write();
    Corr2D->Write();
    Corr_Eta_Phi->Write();
    Corr_Pt_R->Write();
    Smearing_pT_dpT_eta_phi_0_10->Write();
    Smearing_pT_dpT_eta_phi_10_40->Write();
    Smearing_pT_dpT_eta_phi_40_80->Write();

}
