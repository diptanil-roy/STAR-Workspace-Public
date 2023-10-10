void DN(TH1F *R, int col, int mar, double sum,int doit);
void Prelim(){
    //=============
    int idx = 4;
    int idx2 = 4;
    int idx3 = 4;
    int method = 23;
    int method2 = 27;
    int method3 = 18;
    int sys = 0;

    char legg1[100]="Gaussian 0-10% D^{0} p_{T}>5 GeV";
    char legg2[100]="Gassian 10-40% D^{0} p_{T}>5 GeV";
    char legg3[100]="Sampling 0-10% D^{0} p_{T}>5 GeV";
    char legg4[100]="Sampling 10-40% D^{0} p_{T}>5 GeV";
    char legg5[100]="0-10% D^{0} p_{T}>4 GeV";
    char legg6[100]="10-40% D^{0} p_{T}>4 GeV";
    //=============
    TFile *fin1 = new TFile("/gpfs01/star/pwg/droy1/D0Jets/D0Yield_Feb22.root");
    fin1->cd("D0Tagger");
    TH1F *hCent = (TH1F *)gDirectory->Get("hCentralityWeightedAfterCuts");

    for(int i =1;i<hCent->GetNbinsX()+1;i++)cout << hCent->GetBinContent(i) << endl;
    float nevts_0_10 = hCent->Integral(1, 2);
    float nevts_10_40 = hCent->Integral(3, 8);
    float nevts_40_80 = hCent->Integral(9, 16);

    gROOT->ProcessLine(".x ~/myStyle.C");
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();


    TFile *fin1 = new TFile(Form("../hists/Unfolded_Hist_Cen_0_10_Method_%i_Indx_%i.root",method,idx),"READ");
    hData_0_10_1_4 = (TH1D*)fin1->Get(Form("hData_0_10_%i_%i",method,idx));
    TFile *fin2 = new TFile(Form("../hists/Unfolded_Hist_Cen_10_40_Method_%i_Indx_%i.root",method,idx),"READ");
    hData_10_40_1_4 = (TH1D*)fin2->Get(Form("hData_10_40_%i_%i",method,idx));
    TFile *fin3 = new TFile(Form("../hists/Unfolded_Hist_Cen_40_80_Method_%i_Indx_%i.root",method,idx),"READ");
    hData_40_80_1_4 = (TH1D*)fin3->Get(Form("hData_40_80_%i_%i",method,idx));
    
    TFile *fin11 = new TFile(Form("../hists/Unfolded_Hist_Cen_0_10_Method_%i_Indx_%i.root",method2,idx2),"READ");
    hData_0_10_2_4 = (TH1D*)fin11->Get(Form("hData_0_10_%i_%i",method2,idx2));
    TFile *fin21 = new TFile(Form("../hists/Unfolded_Hist_Cen_10_40_Method_%i_Indx_%i.root",method2,idx2),"READ");
    hData_10_40_2_4 = (TH1D*)fin21->Get(Form("hData_10_40_%i_%i",method2,idx2));
    TFile *fin31 = new TFile(Form("../hists/Unfolded_Hist_Cen_40_80_Method_%i_Indx_%i.root",method2,idx2),"READ");
    hData_40_80_2_4 = (TH1D*)fin31->Get(Form("hData_40_80_%i_%i",method2,idx2));
        
    TFile *fin111 = new TFile(Form("../hists/Unfolded_Hist_Cen_0_10_Method_%i_Indx_%i.root",method3,idx3),"READ");
    hData_0_10_3_4 = (TH1D*)fin111->Get(Form("hData_0_10_%i_%i",method3,idx3));
    TFile *fin211 = new TFile(Form("../hists/Unfolded_Hist_Cen_10_40_Method_%i_Indx_%i.root",method3,idx3),"READ");
    hData_10_40_3_4 = (TH1D*)fin211->Get(Form("hData_10_40_%i_%i",method3,idx3));
    TFile *fin311 = new TFile(Form("../hists/Unfolded_Hist_Cen_40_80_Method_%i_Indx_%i.root",method3,idx3),"READ");
    hData_40_80_3_4 = (TH1D*)fin311->Get(Form("hData_40_80_%i_%i",method3,idx3));


    DN(hData_0_10_2_4,2,25,1,1,30);
    DN(hData_10_40_2_4,9,25,1,1,30);
    DN(hData_40_80_2_4,32,25,1,1,30);
    
    DN(hData_0_10_1_4,2,8,1,1,30);
    DN(hData_10_40_1_4,9,8,1,1,30);
    DN(hData_40_80_1_4,32,8,1,1,30);
    
    DN(hData_0_10_3_4,2,24,1,1,30);
    DN(hData_10_40_3_4,9,24,1,1,30);
    DN(hData_40_80_3_4,32,24,1,1,30);

    double y1 = 6.97025e+07+6.9921e+07;
    double y2 = 1.35716e+08 +1.25515e+08 + 1.17859e+08;
    double y3 = 4.9484e+07  + 8.81849e+07 +1.08515e+08+1.13374e+08 ;
    double s1 = 941.23714*nevts_0_10; //(1048*1.38188e+06+838*1.22039e+06)/(1.38188e+06+1.22039e+06);                                                                                                                              
    double s2 = 391.35550*nevts_10_40;//571 + 351 + 206;                                                                                                                                                                           
    double s3 = 56.62475*nevts_40_80;
    

    hData_0_10_1_4->Scale(1./s1);
    hData_10_40_1_4->Scale(1./s2);
    hData_40_80_1_4->Scale(1./s3);


    hData_0_10_2_4->Scale(1./s1);
    hData_10_40_2_4->Scale(1./s2);
    hData_40_80_2_4->Scale(1./s3);

    hData_0_10_3_4->Scale(1./s1);
    hData_10_40_3_4->Scale(1./s2);
    hData_40_80_3_4->Scale(1./s3);


    TCanvas *c2 = new TCanvas("c2","c2");
    hData_0_10_1_4->DrawClone("PE X0 same");
    hData_10_40_1_4->DrawClone("PE X0 same");
    hData_40_80_1_4->DrawClone("PE X0 same");
    if(sys >0){
	hData_0_10_2_4->DrawClone("PE X0 same");
	hData_10_40_2_4->DrawClone("PE X0 same");
	hData_40_80_2_4->DrawClone("PE X0 same");
    }
    if(sys == 2){
        hData_0_10_3_4->DrawClone("PE X0 same");
        hData_10_40_3_4->DrawClone("PE X0 same");
        hData_40_80_3_4->DrawClone("PE X0 same");
    }
    gPad->SetLogy();



    hData_0_10_1_4->Divide(hData_40_80_1_4);
    hData_10_40_1_4->Divide(hData_40_80_1_4);
    hData_0_10_2_4->Divide(hData_40_80_2_4);
    hData_10_40_2_4->Divide(hData_40_80_2_4);
    hData_0_10_3_4->Divide(hData_40_80_3_4);
    hData_10_40_3_4->Divide(hData_40_80_3_4);
    /*
    getRatio(hData_0_10_1_4,hData_40_80_1_4);
    getRatio(hData_10_40_1_4,hData_40_80_1_4);
    getRatio(hData_0_10_2_4,hData_40_80_2_4);
    getRatio(hData_10_40_2_4,hData_40_80_2_4);
    getRatio(hData_0_10_1_4,hData_40_80_3_4);
    getRatio(hData_10_40_1_4,hData_40_80_3_4);
    */
    TLegend *leg  = new TLegend(0.7,0.5,0.9,0.9);
    leg->AddEntry(hData_0_10_1_4,legg1,"PE");
    leg->AddEntry(hData_10_40_1_4,legg2,"PE");
    if(sys >0){
        leg->AddEntry(hData_0_10_2_4,legg3,"PE");
        leg->AddEntry(hData_10_40_2_4,legg4,"PE");
    }
    if(sys ==2){
        leg->AddEntry(hData_0_10_3_4,legg5,"PE");
        leg->AddEntry(hData_10_40_3_4,legg6,"PE");
    }

    TCanvas *c3 = new TCanvas("c3","c3");
    hData_0_10_1_4->GetYaxis()->SetTitle("R_{CP}(/40-80\%)");
    hData_0_10_1_4->GetYaxis()->SetRangeUser(0,2);
    hData_0_10_1_4->GetXaxis()->SetRangeUser(2,30);

    hData_0_10_1_4->DrawClone("PE X0 same");
    hData_10_40_1_4->DrawClone("PE X0 same");
    if(sys > 0){
	hData_0_10_2_4->DrawClone("PE  X0 same");
	hData_10_40_2_4->DrawClone("PE X0 same");
    }
    if(sys == 2){
        hData_0_10_3_4->DrawClone("PE  X0 same");
        hData_10_40_3_4->DrawClone("PE X0 same");
    }
    leg->Draw("same");
    TF1 *line = new TF1("line","1",-100,100);line->SetLineStyle(7);
    line->Draw("same");

    TCanvas *c33 = new TCanvas("c33","c33");
    hData_0_10_1_4->GetYaxis()->SetTitle("R_{CP}(0-10%/40-80\%) Ratio");
    hData_0_10_1_4->GetYaxis()->SetRangeUser(0,2);
    hData_0_10_1_4->GetXaxis()->SetRangeUser(2,30);
    hData_0_10_1_4->Divide(hData_0_10_2_4);
    hData_0_10_1_4->DrawClone("hist same");
    line->Draw("same");
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
    double sum = 1;//R->Integral();
    if(doit){
	for(int i = 1;i<R->GetNbinsX()+1;i++){
	    double val = R->GetBinContent(i);
	    double er = R->GetBinError(i);
	    double width = R->GetBinWidth(i);
	    double center = fabs(R->GetBinCenter(i));
	    if(center<cutoff){
		R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.035);
	    R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.035);
	    }
	    else{
		R->SetBinContent(i,0);
		R->SetBinError(i,0);
	    }
	}
    }
    R->SetLineColor(col);
    R->SetMarkerColor(col);
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
