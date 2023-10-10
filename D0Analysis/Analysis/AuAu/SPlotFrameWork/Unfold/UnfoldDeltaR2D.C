#include "../Config.h"
void doSys(){
    for(int i = 0;i<166;i++)
        RUN(1,i);
}
void RUN(int sys=0, int i=0){
    int method = 1;
    if(sys==0){
	UnfoldDeltaR2D(0,10,method,4,1,0,0);
	UnfoldDeltaR2D(10,40,method,4,1,0,0);
	UnfoldDeltaR2D(40,80,method,4,1,0,0);
    }
    else{
	UnfoldDeltaR2D(0,10,method,4,1,1,i);
	UnfoldDeltaR2D(10,40,method,4,1,1,i);
	UnfoldDeltaR2D(40,80,method,4,1,1,i);
    }
}
void UnfoldDeltaR2D(int cen_high=0,int cen_low=10,  int method=1,int index=4,int save=1,int sys=0, int iter=0){
    
    gSystem->Load("RooUnfold/libRooUnfold");    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    // Get data histograms

    if(sys!=1 && method!=26){
	if(method==1 || method ==25)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	if(method==2)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D04GeV.root","READ");
	else if (method==5)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D0PT_3_10.root","READ");
    }else if(sys==1){
	if(method==1)TFile *fin1 = new TFile(Form("../ApplyWeights/systematics_new_ratios/Histograms3_D05GeV_%i.root",iter),"READ");
    }
    else if(method==26){
	TFile *fin1 = new TFile("D0JetsFoldedHistogramsNew.root","READ");
    }
    if(method!=26){
    if     (cen_high==0  && cen_low==10) TH1D *data = (TH1D*) fin1->Get("JetPt_R_0_10");
    else if(cen_high==10 && cen_low==20) TH1D *data = (TH1D*) fin1->Get("JetPt_R_10_20");
    else if(cen_high==20 && cen_low==30) TH1D *data = (TH1D*) fin1->Get("JetPt_R_20_30");
    else if(cen_high==30 && cen_low==40) TH1D *data = (TH1D*) fin1->Get("JetPt_R_30_40");
    else if(cen_high==40 && cen_low==80) TH1D *data = (TH1D*) fin1->Get("JetPt_R_40_80");
    else if(cen_high==10 && cen_low==40) TH1D *data = (TH1D*) fin1->Get("JetPt_R_10_40");
    else if(cen_high==0  && cen_low==80) TH1D *data = (TH1D*) fin1->Get("JetPt_R");
    }
    else{
	if     (cen_high==0  && cen_low==10) TH1D *data = (TH1D*) fin1->Get("JetPtvdR_Cent_0-10_D0Pt_1-10");
        else if(cen_high==40 && cen_low==80) TH1D *data = (TH1D*) fin1->Get("JetPtvdR_Cent_40-80_D0Pt_1-10");
        else if(cen_high==10 && cen_low==40) TH1D *data = (TH1D*) fin1->Get("JetPtvdR_Cent_10-40_D0Pt_1-10");

	
    }
//    TFile *fin1 = new TFile("D0JetsFoldedHistograms_Neil.root");
//    if     (cen_high==0  && cen_low==10) TH1D *data = (TH1D*) fin1->Get("JetPt_Cent_0-10_D0Pt_1-10");
//    else if(cen_high==40 && cen_low==80) TH1D *data = (TH1D*) fin1->Get("JetPt_Cent_40-80_D0Pt_1-10");
//    else if(cen_high==10 && cen_low==40) TH1D *data = (TH1D*) fin1->Get("JetPt_Cent_10-40_D0Pt_1-10");


    // ============== Done data histograms      
    // ============== First get response matrix
    if(method==1 || method==26)TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New_Sample_FINAL_test.root",cen_high,cen_low),"READ");
    else  if(method==25)TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_PYTHIA_D05GeV_Reco3_New_Sample_FINAL_test.root",cen_high,cen_low),"READ");
    else if (method==2) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D04GeV_Reco3_New_Sample_FINAL.root",cen_high,cen_low),"READ");
    else if (method==3) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D03GeV.root",cen_high,cen_low),"READ");
    else if (method==4) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i.root",cen_high,cen_low),"READ");
    else if (method==5) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D0PT_3_10.root",cen_high,cen_low),"READ");
    

    RooUnfoldResponse *response;
    fin->GetObject("Response2D",response);
    TH2F *signalJetPtvdRD0Jet = fin->Get("signalJetPtvdRD0Jet");
    TH1F* hReco = (TH1F*)fin->Get("hReco_dR"); 
    TH1F* hTruth = (TH1F*)fin->Get("hTruth_dR");
    TH2F* hRes = (TH2F*)fin->Get("hRes_dR");
    
//    RooUnfoldResponse *response = new RooUnfoldResponse(hReco,hTruth,hRes,"response","response");
    
    TMatrixD matx = response->Mresponse();
    TH2D *res_matx = new TH2D(matx);
    // ============== Done with response matrix
    // ============== Unfold
    RooUnfoldBayes unfold (response, data, index);
    TH2D* hData2D = (TH2D*) unfold.Hreco();
    hData2D->GetXaxis()->SetRangeUser(5,15);
    TH1D* hData = (TH1D*)hData2D->ProjectionY();
    // =============== Now draw
    if(sys!=1){
	gROOT->ProcessLine(".x ~/myStyle.C"); 
	TLatex lat;
	char label[100];
	sprintf(label,"%i-%i percentile",cen_high,cen_low);
	
	TCanvas *c1 = new TCanvas("c1","c1");
	res_matx->GetYaxis()->SetTitle("True p_{T} Bin Index");
	res_matx->GetXaxis()->SetTitle("Reco. p_{T} Bin Index");
	res_matx->Draw("COLZ");
	lat.DrawLatex(6,36,label);
	gPad->SetLogz();
	
	TCanvas *c2 = new TCanvas("c2","c2");
	hData->GetXaxis()->SetTitle("D^{0} #Delta R");
	hData->GetYaxis()->SetTitle("Counts");
	data->SetMarkerStyle(24);
	hData->SetMarkerColor(kBlue);
	hData->SetLineColor(kBlue);
	TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
	leg->SetHeader(label);
	leg->AddEntry(data,"Data","PE");
	leg->AddEntry(hData,"Unfolded","PE");
	hData->GetYaxis()->SetRangeUser(1,hData->GetMaximum()*5);
	hData->Draw("PE");
	data->Draw("pe same");
	leg->Draw("same");
	gPad->SetLogy();
    }
    if(save){
	char out[100];
	if(sys==0)sprintf(out,"hists/Unfolded_Hist_Cen_%i_%i_Method_%i_Indx_%i_dR2D.root",cen_high,cen_low,method,index);
	else sprintf(out,"hists/sys_ratios/Unfolded_Hist_Cen_%i_%i_Method_%i_Indx_%i_Sys_%i_dR2D.root",cen_high,cen_low,method,index,iter);
	TFile *fout = new TFile(out,"RECREATE");
	char hout[100];
	sprintf(hout,"hData_%i_%i_%i_%i",cen_high,cen_low,method,index);
	hData->Write(hout);
	fout->Close();
    }
}
