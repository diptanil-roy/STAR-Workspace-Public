#include "../Config.h"

void doSys(){
    for(int i = 0;i<166;i++)
	RUN(1,i);
}

void RUN(int sys=0, int i=0){
    int method = 21;
    if(sys==0){
	UnfoldAll(0,10,method,4,1,0,0);
	UnfoldAll(10,40,method,4,1,0,0);
	UnfoldAll(40,80,method,4,1,0,0);
    }
    else{
	UnfoldAll(0,10,method,4,1,1,i);
        UnfoldAll(10,40,method,4,1,1,i);
        UnfoldAll(40,80,method,4,1,1,i);
    }
}
void UnfoldAll(int cen_high=0,int cen_low=10,  int method=21,int index=4,int save=0,int sys=0, int iter=0){
    
    gSystem->Load("RooUnfold/libRooUnfold");    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    // Get data histograms

    if(sys!=1 && method!=26){
	if(method==1 || method==2 || method==3 || method ==4 || method==5)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D03GeV.root","READ");
	else if (method==6 || method==7)TFile *fin1 = new TFile("../ApplyWeights/Histograms7_D05GeV.root","READ");
	else if (method==8)TFile *fin1 = new TFile("../ApplyWeights/HistogramsVary_D02GeV.root","READ");
	else if (method==10)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D02GeV.root","READ");
	else if (method==11)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D03GeV.root","READ");
	else if (method==12)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	else if (method==13)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D02GeV.root","READ");
	else if (method==14)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D03GeV.root","READ");
	else if (method==15)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	else if (method==16)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D02GeV.root","READ");
        else if (method==17)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D03GeV.root","READ");
        else if (method==18)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	else if (method==19)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D04GeV.root","READ");
	else if (method==20)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D03GeV.root","READ");
	else if (method==21)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	else if (method==22)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D04GeV.root","READ");
	else if (method==23)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D02GeV.root","READ");
	else if (method==24)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D04GeV.root","READ");
	else if (method==25)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	else if (method==27)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	else if (method==28)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
	else if (method==29)TFile *fin1 = new TFile("../ApplyWeights/Histograms3_D05GeV.root","READ");
    }else if(sys==1){
	TFile *fin1 = new TFile(Form("../ApplyWeights/systematics_new/Histograms3_D05GeV_%i.root",iter),"READ");
    }else if(method==26){
        TFile *fin1 = new TFile("D0JetsFoldedHistogramsNew.root","READ");
    }
    if(method!=26){
    if     (cen_high==0  && cen_low==10) TH1D *data = (TH1D*) fin1->Get("JetPt_0_10");
    else if(cen_high==10 && cen_low==20) TH1D *data = (TH1D*) fin1->Get("JetPt_10_20");
    else if(cen_high==20 && cen_low==30) TH1D *data = (TH1D*) fin1->Get("JetPt_20_30");
    else if(cen_high==30 && cen_low==40) TH1D *data = (TH1D*) fin1->Get("JetPt_30_40");
    else if(cen_high==40 && cen_low==80) TH1D *data = (TH1D*) fin1->Get("JetPt_40_80");
    else if(cen_high==10 && cen_low==40) TH1D *data = (TH1D*) fin1->Get("JetPt_10_40");
    else if(cen_high==0  && cen_low==80) TH1D *data = (TH1D*) fin1->Get("JetPt");
    }
    else{
	if     (cen_high==0  && cen_low==10)TH1D *data1 = (TH1D*) fin1->Get("JetPt_Cent_0-10_D0Pt_1-10");
	else if(cen_high==40 && cen_low==80) TH1D *data1 = (TH1D*) fin1->Get("JetPt_Cent_40-80_D0Pt_1-10");
	else if(cen_high==10 && cen_low==40) TH1D *data1 = (TH1D*) fin1->Get("JetPt_Cent_10-40_D0Pt_1-10");      
	TH1D *data =new TH1D("data","data",nbins_jpt,binning_jpt);
	for(int i = 1;i<data->GetNbinsX()+1;i++){
	    data->SetBinContent(i,data1->GetBinContent(i+3));
	    data->SetBinError(i,data1->GetBinError(i+3));
	}

    }
//    TFile *fin1 = new TFile("D0JetsFoldedHistograms_Neil.root");
//    if     (cen_high==0  && cen_low==10) TH1D *data = (TH1D*) fin1->Get("JetPt_Cent_0-10_D0Pt_1-10");
//    else if(cen_high==40 && cen_low==80) TH1D *data = (TH1D*) fin1->Get("JetPt_Cent_40-80_D0Pt_1-10");
//    else if(cen_high==10 && cen_low==40) TH1D *data = (TH1D*) fin1->Get("JetPt_Cent_10-40_D0Pt_1-10");


    // ============== Done data histograms      
    // ============== First get response matrix
    if(1){
	if(method==1)TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_D02GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==2) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D02GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==3) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_PeripheralWeights_D03GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==4) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D03GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==5) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D05GeV_Reco3.root",cen_high,cen_low),"READ");
        else if (method==6) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D05GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==7) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D05GeV_Reco7_Scale.root",cen_high,cen_low),"READ");
        else if (method==8) TFile *fin = new TFile(Form("./ResponseMatrices_%i_%i_FONLL_D02GeV_RecoVary.root",cen_high,cen_low),"READ");
	else if (method==10) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D02GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==11) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D03GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==12) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3.root",cen_high,cen_low),"READ");
	else if (method==13) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D02GeV_Reco3_New.root",cen_high,cen_low),"READ");
	else if (method==14) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D03GeV_Reco3_New.root",cen_high,cen_low),"READ");
	else if (method==15) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New.root",cen_high,cen_low),"READ");
	else if (method==16) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_PYTHIA_D02GeV_Reco3_New.root",cen_high,cen_low),"READ");
        else if (method==17) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_PYTHIA_D03GeV_Reco3_New.root",cen_high,cen_low),"READ");
        else if (method==18) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_PYTHIA_D05GeV_Reco3_New.root",cen_high,cen_low),"READ");
	else if (method==19) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D04GeV_Reco3_New_Sample_FINAL.root",cen_high,cen_low),"READ");
        else if (method==20) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D03GeV_Reco3_New_Sample.root",cen_high,cen_low),"READ");
        else if (method==21) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New_Sample_FINAL.root",cen_high,cen_low),"READ");
	else if (method==22) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D04GeV_Reco3_New.root",cen_high,cen_low),"READ");
	else if (method==23) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D02GeV_Reco3_New_Sample_FINAL.root",cen_high,cen_low),"READ");
	else if (method==24) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_PYTHIA_D04GeV_Reco3_New.root",cen_high,cen_low),"READ");
	else if (method==25) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_PYTHIA_D05GeV_Reco3_New_Sample_FINAL.root",cen_high,cen_low),"READ");
	else if (method==26) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New_Sample_FINAL.root",cen_high,cen_low),"READ");
	else if (method==27) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New_Sample_WEIGHT3.root",cen_high,cen_low),"READ");
	else if (method==28) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New_Sample.root",cen_high,cen_low),"READ");
	else if (method==29) TFile *fin = new TFile(Form("./ResponseMatricesNew_%i_%i_FONLL_D05GeV_Reco3_New.root",cen_high,cen_low),"READ");
	TH1F* hReco = (TH1F*)fin->Get("hReco"); 
	TH1F* hTruth = (TH1F*)fin->Get("hTruth");
	TH2F* hRes = (TH2F*)fin->Get("hRes");
	RooUnfoldResponse *response = new RooUnfoldResponse(hReco,hTruth,hRes,"response","response");
    }
    else{
	//Neils response matrix
	if(cen_high==0 && cen_low==10)TFile *ResponseFile = new TFile("ResponseMatrix_Prelim_Central.root","READ");
	if(cen_high==10 && cen_low==40)TFile *ResponseFile = new TFile("ResponseMatrix_Prelim_MidCentral.root","READ");
	if(cen_high==40 && cen_low==80)TFile *ResponseFile = new TFile("ResponseMatrix_Prelim_Peripheral.root","READ");
	TH1D *hMCJetPt = (TH1D*)ResponseFile->Get("hTrueMCPt");
	TH1D *hMCRecoSmearedJetPt = (TH1D*)ResponseFile->Get("hTrueRecoPt")->Clone("hMCRecoSmearedJetPt");
	TH2D *ResponseJetPt = (TH2D*)ResponseFile->Get("hTrueRecoMCPt")->Clone("Response");
	RooUnfoldResponse * response = new RooUnfoldResponse(hMCRecoSmearedJetPt, hMCJetPt, ResponseJetPt,"response","response");
    }

    
    TMatrixD matx = response->Mresponse();
    TH2D *res_matx = new TH2D(matx);
    // ============== Done with response matrix
    // ============== Unfold
    RooUnfoldBayes unfold (response, data, index);
    TH1D* hData = (TH1D*) unfold.Hreco();
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
	hData->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
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
	if(sys==0)sprintf(out,"hists/Unfolded_Hist_Cen_%i_%i_Method_%i_Indx_%i.root",cen_high,cen_low,method,index);
	else sprintf(out,"hists/sys/Unfolded_Hist_Cen_%i_%i_Method_%i_Indx_%i_Sys_%i.root",cen_high,cen_low,method,index,iter);
	TFile *fout = new TFile(out,"RECREATE");
	char hout[100];
	sprintf(hout,"hData_%i_%i_%i_%i",cen_high,cen_low,method,index);
	hData->Write(hout);
	fout->Close();
    }
}
