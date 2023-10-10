

#include "TFile.h"
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TChain.h"
#include "TF1.h"

using namespace std;
void EtaProb(int save=1){
    gROOT->ProcessLine(".x ~/myStyle.C");
    //TFile *f_D = new TFile("../Production/output/full_new.root","READ");
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
    TH1D *ETA_1  = new TH1D("ETA_1","ETA_1",50,0.,1.);
    TH1D *ETA_2  = new TH1D("ETA_2","ETA_2",50,0.,1.);
    TH1D *ETA_3  = new TH1D("ETA_3","ETA_3",50,0.,1.);
    TH1D *ETA_1_2_5  = new TH1D("ETA_1_2_5","ETA_1_2_5",50,0.,1.);
    TH1D *ETA_2_2_5  = new TH1D("ETA_2_2_5","ETA_2_2_5",50,0.,1.);
    TH1D *ETA_3_2_5   = new TH1D("ETA_3_2_5","ETA_3_2_5",50,0.,1.);
    TH1D *ETA_1_5_10   = new TH1D("ETA_1_5_10 ","ETA_1_5_10 ",50,0.,1.);
    TH1D *ETA_2_5_10   = new TH1D("ETA_2_5_10 ","ETA_2_5_10 ",50,0.,1.);
    TH1D *ETA_3_5_10  = new TH1D("ETA_3_5_10 ","ETA_3_5_10 ",50,0.,1.);
    TH1D *ETA_1_10_30  = new TH1D("ETA_1_10_30","ETA_1_10_30",50,0.,1.);
    TH1D *ETA_2_10_30  = new TH1D("ETA_2_10_30","ETA_2_10_30",50,0.,1.);
    TH1D *ETA_3_10_30  = new TH1D("ETA_3_10_30","ETA_3_10_30",50,0.,1.);
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
	if(gpt<2.)continue;
	if(fabs(jeta)>=0.6.)continue;
	Weight = 1;	
	if(ccen>=0 && ccen<10)ETA_1->Fill(geta);
	if(ccen>=10 && ccen<40)ETA_2->Fill(geta);
	if(ccen>=40 && ccen<80)ETA_3->Fill(geta);
	if(gpt>=2 && gpt<5){
	    if(ccen>=0 && ccen<10)ETA_1_2_5->Fill(geta);
	    if(ccen>=10 && ccen<40)ETA_2_2_5->Fill(geta);
	    if(ccen>=40 && ccen<80)ETA_3_2_5->Fill(geta);
	}
	else if(gpt>=5 && gpt<10){
            if(ccen>=0 && ccen<10)ETA_1_5_10->Fill(geta);
            if(ccen>=10 && ccen<40)ETA_2_5_10->Fill(geta);
            if(ccen>=40 && ccen<80)ETA_3_5_10->Fill(geta);
	}
	else if(gpt>10 && gpt<1000){
            if(ccen>=0 && ccen<10)ETA_1_10_30->Fill(geta);
            if(ccen>=10 && ccen<40)ETA_2_10_30->Fill(geta);
            if(ccen>=40 && ccen<80)ETA_3_10_30->Fill(geta);
        }
    }

    ETA_1->SetLineColor(kRed);
    ETA_2->SetLineColor(kBlue);
    ETA_3->SetLineColor(kGreen-2);

    TF1 *f1 = new TF1("f1","[0]",0,0.4);
    TF1 *f2 = new TF1("f2","[0]",0,0.4);
    TF1 *f3 = new TF1("f3","[0]",0,0.4);

    TF1 *f1_1 = new TF1("f1_1","[0]",0,0.4);
    TF1 *f2_1 = new TF1("f2_1","[0]",0,0.4);
    TF1 *f3_1 = new TF1("f3_1","[0]",0,0.4);
    TF1 *f1_2 = new TF1("f1_2","[0]",0,0.4);
    TF1 *f2_2 = new TF1("f2_2","[0]",0,0.4);
    TF1 *f3_2 = new TF1("f3_2","[0]",0,0.4);
    TF1 *f1_3 = new TF1("f1_3","[0]",0,0.4);
    TF1 *f2_3 = new TF1("f2_3","[0]",0,0.4);
    TF1 *f3_3 = new TF1("f3_3","[0]",0,0.4);


    ETA_1->Fit(f1,"R");
    ETA_2->Fit(f2,"R");
    ETA_3->Fit(f3,"R");
    
    ETA_1_2_5->Fit(f1_1,"R");
    ETA_2_2_5->Fit(f2_1,"R");
    ETA_3_2_5->Fit(f3_1,"R");

    ETA_1_5_10->Fit(f1_2,"R");
    ETA_2_5_10->Fit(f2_2,"R");
    ETA_3_5_10->Fit(f3_2,"R");

    ETA_1_10_30->Fit(f1_3,"R");
    ETA_2_10_30->Fit(f2_3,"R");
    ETA_3_10_30->Fit(f3_3,"R");

    TCanvas *c11 = new TCanvas("c11","C11");
    ETA_1->DrawClone("hist same");
    ETA_2->DrawClone("hist same");
    ETA_3->DrawClone("hist same");
    f1->Draw("same");
    f2->Draw("same");
    f3->Draw("same");
    


    TCanvas *c1 = new TCanvas("c1","C1");
    ETA_1->DrawNormalized("hist same");
    ETA_2->DrawNormalized("hist same");
    ETA_3->DrawNormalized("hist same");

    double val1 = f1->GetParameter(0);
    double val2= f2->GetParameter(0);
    double val3= f3->GetParameter(0);

    double val1_1 = f1_1->GetParameter(0);
    double val2_1= f2_1->GetParameter(0);
    double val3_1= f3_1->GetParameter(0);

    double val1_2 = f1_2->GetParameter(0);
    double val2_2= f2_2->GetParameter(0);
    double val3_2= f3_2->GetParameter(0);

    double val1_3 = f1_3->GetParameter(0);
    double val2_3= f2_3->GetParameter(0);
    double val3_3= f3_3->GetParameter(0);


    for(int i = 1;i<ETA_1->GetNbinsX()+1;i++){
	if(fabs(ETA_1->GetBinCenter(i))<0.6){
	    ETA_1->SetBinContent(i,1);
	    ETA_2->SetBinContent(i,1);
	    ETA_3->SetBinContent(i,1);
	    ETA_1_2_5->SetBinContent(i,1);
            ETA_2_2_5->SetBinContent(i,1);
            ETA_3_2_5->SetBinContent(i,1);
	    ETA_1_5_10->SetBinContent(i,1);
            ETA_2_5_10->SetBinContent(i,1);
            ETA_3_5_10->SetBinContent(i,1);
	    ETA_1_10_30->SetBinContent(i,1);
            ETA_2_10_30->SetBinContent(i,1);
            ETA_3_10_30->SetBinContent(i,1);
	}
	else{
	    double num1 = ETA_1->GetBinContent(i);
	    double num2= ETA_2->GetBinContent(i);
	    double num3= ETA_3->GetBinContent(i);
	    ETA_1->SetBinContent(i,num1/val1);
	    ETA_2->SetBinContent(i,num2/val2);
	    ETA_3->SetBinContent(i,num3/val3);
	

	    num1 = ETA_1_2_5->GetBinContent(i);
            num2= ETA_2_2_5->GetBinContent(i);
            num3= ETA_3_2_5->GetBinContent(i);
            ETA_1_2_5->SetBinContent(i,num1/val1_1);
            ETA_2_2_5->SetBinContent(i,num2/val2_1);
            ETA_3_2_5->SetBinContent(i,num3/val3_1);

            num1 = ETA_1_5_10->GetBinContent(i);
            num2= ETA_2_5_10->GetBinContent(i);
            num3= ETA_3_5_10->GetBinContent(i);
            ETA_1_5_10->SetBinContent(i,num1/val1_2);
            ETA_2_5_10->SetBinContent(i,num2/val2_2);
            ETA_3_5_10->SetBinContent(i,num3/val3_2);

            num1 = ETA_1_10_30->GetBinContent(i);
            num2= ETA_2_10_30->GetBinContent(i);
            num3= ETA_3_10_30->GetBinContent(i);
            ETA_1_10_30->SetBinContent(i,num1/val1_3);
            ETA_2_10_30->SetBinContent(i,num2/val2_3);
            ETA_3_10_30->SetBinContent(i,num3/val3_3);

	}
	
    }
    TCanvas *c111 = new TCanvas("c111","C111");
    ETA_1->DrawClone("hist same");
    ETA_1_2_5->DrawClone("hist same");
    ETA_1_5_10->DrawClone("hist same");
    ETA_1_10_30->DrawClone("hist same");

    TCanvas *c1111 = new TCanvas("c1111","C1111");
    ETA_1->DrawClone("hist same");
    ETA_2->DrawClone("hist same");
    ETA_3->DrawClone("hist same");
    
    
    
    if(save){
    TFile* outfile = new TFile("FakeRates_FINAL.root","RECREATE");
    ETA_1->Write();
    ETA_2->Write();
    ETA_3->Write();

    ETA_1_2_5->Write();
    ETA_2_2_5->Write();
    ETA_3_2_5->Write();

    ETA_1_5_10->Write();
    ETA_2_5_10->Write();
    ETA_3_5_10->Write();

    ETA_1_10_30->Write();
    ETA_2_10_30->Write();
    ETA_3_10_30->Write();
    }
}
