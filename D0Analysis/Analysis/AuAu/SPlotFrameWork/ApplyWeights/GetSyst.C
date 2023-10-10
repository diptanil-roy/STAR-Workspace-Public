#include "../Config.h"
void DN(TH1D *R, int col, int mar, double sum,int doit, double cutoff);
void GetSyst(){
    //=============
    int method=21;
    int idx=4;
    TFile *fin = new TFile("../../Files/D0Yield_Feb22.root");
    fin->cd("D0Tagger");
    TH1F *hCent = (TH1F *)gDirectory->Get("hCentralityWeightedAfterCuts");


//    for(int i =1;i<hCent->GetNbinsX()+1;i++)cout << hCent->GetBinContent(i) << endl;
    float nevts_0_10 = hCent->Integral(1, 2);
    float nevts_10_40 = hCent->Integral(3, 8);
    float nevts_40_80 = hCent->Integral(9, 16);


    double s1 = nevts_0_10;//941.23714*nevts_0_10;
    double s2 = nevts_10_40;//391.35550*nevts_10_40;
    double s3 = nevts_40_80;//56.62475*nevts_40_80;
//============= 


    gROOT->ProcessLine(".x ../myStyle.C");
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();


    TFile *fin1 = new TFile("Histograms3_D01_3GeV.root","READ");
    TH1D* hData_0_10_1_4 = (TH1D*)fin1->Get("JetPt_0_10");
    TH1D* hData_10_40_1_4 = (TH1D*)fin1->Get("JetPt_10_40");
    TH1D* hData_40_80_1_4 = (TH1D*)fin1->Get("JetPt_40_80");

    int n_sys=100;
    double bin_sys[201];
    for(int i =0;i<n_sys+1;i++)bin_sys[i]= -2 + i*4./n_sys;
    
    
    //TH2D * SYS1 = new TH2D("SYS1","SYS1",PTBINS,edges,n_sys,bin_sys);
    //TH2D * SYS2 = new TH2D("SYS2","SYS2",PTBINS,edges,n_sys,bin_sys);
    //TH2D * SYS3 = new TH2D("SYS3","SYS3",PTBINS,edges,n_sys,bin_sys);
    
    TH2D * SYS1 = new TH2D("SYS1","SYS1",nbins_jpt,binning_jpt,n_sys,bin_sys);
    TH2D * SYS2 = new TH2D("SYS2","SYS2",nbins_jpt,binning_jpt,n_sys,bin_sys);
    TH2D * SYS3 = new TH2D("SYS3","SYS3",nbins_jpt,binning_jpt,n_sys,bin_sys);

    for(int i = 0;i<166;i++){

	TFile *fin11 = new TFile(Form("./systematics_new/Histograms3_D05GeV_%i.root",i),"READ");
	TH1D* temp1 = (TH1D*)fin11->Get(Form("JetPt_0_10"));
	TH1D* temp2 = (TH1D*)fin11->Get(Form("JetPt_10_40"));
	TH1D* temp3 = (TH1D*)fin11->Get(Form("JetPt_40_80"));

	for(int j=1;j<hData_0_10_1_4->GetNbinsX()+1;j++){
	    double val = hData_0_10_1_4->GetBinContent(j);
	    if(val<=0)continue;
	    double val_new = temp1->GetBinContent(j);
	    SYS1->Fill(hData_0_10_1_4->GetBinCenter(j),val_new/val);
	}
	for(int j=1;j<hData_10_40_1_4->GetNbinsX()+1;j++){
            double val = hData_10_40_1_4->GetBinContent(j);
            if(val<=0)continue;
            double val_new = temp2->GetBinContent(j);
            SYS2->Fill(hData_10_40_1_4->GetBinCenter(j),val_new/val);
	} 
	for(int j=1;j<hData_40_80_1_4->GetNbinsX()+1;j++){
            double val = hData_40_80_1_4->GetBinContent(j);
            if(val<=0)continue;
            double val_new = temp3->GetBinContent(j);
            SYS3->Fill(hData_40_80_1_4->GetBinCenter(j),val_new/val);
	}
	fin11->Close();
    }
    
    TH1D* JetSpectra_0_10_Sys_D0  = (TH1D*)hData_0_10_1_4->Clone("JetSpectra_0_10_Sys_D0");
    TH1D* JetSpectra_10_40_Sys_D0 = (TH1D*)hData_10_40_1_4->Clone("JetSpectra_10_40_Sys_D0");
    TH1D* JetSpectra_40_80_Sys_D0 = (TH1D*)hData_40_80_1_4->Clone("JetSpectra_40_80_Sys_D0");


    for(int j=1;j<hData_0_10_1_4->GetNbinsX()+1;j++){

	double val1 = hData_0_10_1_4->GetBinContent(j);
	TH1D* _temp1 = (TH1D*)SYS1->ProjectionY("_temp1",j,j);
	double er1 = _temp1->GetRMS();
	JetSpectra_0_10_Sys_D0->SetBinContent(j,val1);
	JetSpectra_0_10_Sys_D0->SetBinError(j,val1*er1);

	double val2 = hData_10_40_1_4->GetBinContent(j);
	TH1D* _temp2 =(TH1D*)SYS2->ProjectionY("_temp2",j,j);
	double er2 = _temp2->GetRMS();
	JetSpectra_10_40_Sys_D0->SetBinContent(j,val2);
	JetSpectra_10_40_Sys_D0->SetBinError(j,val2*er2);

	double val3 = hData_40_80_1_4->GetBinContent(j);
	TH1D* _temp3 =(TH1D*)SYS3->ProjectionY("_temp3",j,j);
	double er3 = _temp3->GetRMS();
	JetSpectra_40_80_Sys_D0->SetBinContent(j,val3);
	JetSpectra_40_80_Sys_D0->SetBinError(j,val3*er3);
	
	cout << "Rel errors: D0 " << hData_0_10_1_4->GetBinCenter(j) << " " << er1 << " " << er2 << " " << er3 << endl;

    }
    

    DN(hData_0_10_1_4,2,8,1,1,30);
    DN(hData_10_40_1_4,9,8,1,1,30);
    DN(hData_40_80_1_4,32,8,1,1,30);
    DN(JetSpectra_0_10_Sys_D0,2,8,1,1,30);
    DN(JetSpectra_10_40_Sys_D0,9,8,1,1,30);
    DN(JetSpectra_40_80_Sys_D0,32,8,1,1,30);    
    hData_0_10_1_4->Scale(1./s1);
    hData_10_40_1_4->Scale(1./s2);
    hData_40_80_1_4->Scale(1./s3);

    JetSpectra_0_10_Sys_D0->Scale(1./s1);
    JetSpectra_10_40_Sys_D0->Scale(1./s2);
    JetSpectra_40_80_Sys_D0->Scale(1./s3);
    TCanvas *c2 = new TCanvas("c2","c2");
    hData_0_10_1_4->DrawClone("PE X0 same");
    JetSpectra_0_10_Sys_D0->DrawClone("same E2");
    gPad->SetLogy();
    TCanvas *c21 = new TCanvas("c21","c21");
    hData_10_40_1_4->DrawClone("PE X0 same");
    JetSpectra_10_40_Sys_D0->DrawClone("same E2");
    gPad->SetLogy();
    TCanvas *c22 = new TCanvas("c22","c22");
    hData_40_80_1_4->DrawClone("PE X0 same");
    JetSpectra_40_80_Sys_D0->DrawClone("same E2");
    gPad->SetLogy();

}

void getRatio(TH1D *N,TH1D *D){

    for(int i = 1;i<N->GetNbinsX()+1;i++){
	double val1 = N->GetBinContent(i);
	double val2 = D->GetBinContent(i);
	double er1 = N->GetBinError(i)/ N->GetBinContent(i);
	double er2 = D->GetBinError(i)/D->GetBinContent(i);
	if(val1>0 && val2>0){
	    N->SetBinContent(i,val1/val2);
	    N->SetBinError(i,val1/val2*sqrt(er1*er1+er2*er2));
	    
	}
	else{
	    N->SetBinContent(i,-1);
            N->SetBinError(i,0);
	}

    }
} 

void DN(TH1D *R, int col, int mar, double sum,int doit, double cutoff){
    sum = 1;//R->Integral();
    if(doit){
	for(int i = 1;i<R->GetNbinsX()+1;i++){
	    double val = R->GetBinContent(i);
	    double er = R->GetBinError(i);
	    double width = R->GetBinWidth(i);
	    double center = fabs(R->GetBinCenter(i));
	    if(center<cutoff){
		R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.039/2);//Final divide by two is for (D0 + D0-bar)/2.
		R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.039/2.);//Final divide by two is for (D0 + D0-bar)/2.  
	    }
	    else{
		R->SetBinContent(i,0);
		R->SetBinError(i,0);
	    }
	}
    }
    R->SetLineColor(col);
    R->SetMarkerColor(col);
    R->SetFillColor(col);
    R->SetFillStyle(3004);
    R->SetMarkerStyle(mar);
    R->GetYaxis()->SetTitle("(1/2#pi) d^{2}N/d#it{p}_{T}d#it{y} (#it{c}^{2}/GeV^{2})");
    R->GetXaxis()->SetTitle("Jet #it{p}_{T} (GeV/#it{c})");
    //R->GetYaxis()->SetRangeUser(0.000001,R->GetMaximum()*10);
    R->GetYaxis()->SetRangeUser(0.01,5000000);//R->GetMaximum()*10);
    R->GetXaxis()->SetRangeUser(3,30);
}
void getData(TGraph *g,char file[]){
    int cnt=0;
    ifstream data(file);
    if(data.is_open()){
        while(!data.eof()){
            double x;double y;
            data >> x >> y;
            g->SetPoint(cnt,x,y/x/3000/2);cnt++;
        }
    }
}
