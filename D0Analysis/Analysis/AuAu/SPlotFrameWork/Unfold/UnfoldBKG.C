#include "../Config.h"

void UnfoldBKG(float cen_high=40,float cen_low=80, int method=100,int index=4,int save=1, int sim =0){
    
    double reco_low=-1000;
    double gen_low = 3;//Always for PYTIHA
    if(method==1 || method==21 )reco_low= -1000;
    if(method==2 || method==22 )reco_low= 3;
    if(method==3)reco_low=-1000;
    if(method==4)reco_low= 3;
    if(method==5 || method == 25)reco_low= 5;
    if(method==9)reco_low= 5;

    gROOT->ProcessLine(".x ~/myStyle.C");  
    gSystem->Load("RooUnfold/libRooUnfold");    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TGraph * FONLL = new TGraph(30);
    getData(FONLL,"c.txt");
    // ============== Load data histograms                                                                                                                                                                                                                              

    cout << "> Loading in data histograms " << endl;
    if(!sim){
	if(method==11)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
	if(method==10)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
	if(method==9)TFile *fin = new TFile("../ApplyWeights/Histograms_3GeV.root","READ");
	if(method==8)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
	if(method==7)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
	if(method==6)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
	if(method==5 || method==25)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
	if(method==3 || method==23)TFile *fin = new TFile("../ApplyWeights/Histograms_3GeV.root","READ");
	if(method==100 || method==21)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
	if(method==2 || method==22)TFile *fin = new TFile("../ApplyWeights/Histograms.root","READ");
        if(method==4 || method==24)TFile *fin = new TFile("../ApplyWeights/Histograms_3GeV.root","READ");
    }
    else TFile *fin = new TFile("ToyHistograms.root","READ");

    if     (cen_high==0  && cen_low==10) TH1D *data = (TH1D*) fin->Get("JetPt_0_10");
    else if(cen_high==10 && cen_low==20) TH1D *data = (TH1D*) fin->Get("JetPt_10_20");
    else if(cen_high==20 && cen_low==30) TH1D *data = (TH1D*) fin->Get("JetPt_20_30");
    else if(cen_high==30 && cen_low==40) TH1D *data = (TH1D*) fin->Get("JetPt_30_40");
    else if(cen_high==40 && cen_low==80) TH1D *data = (TH1D*) fin->Get("JetPt_40_80");
    else if(cen_high==0  && cen_low==80) TH1D *data = (TH1D*) fin->Get("JetPt");
    cout << "> Done loading in data histograms " << endl;

    // ============== Done data histograms      
    // ============== First get response matrix
    for(int i = 1;i<data->GetNbinsX()+1;i++){
	if(data->GetBinCenter(i)<reco_low){
	    data->SetBinContent(i,0);
	    data->SetBinError(i,0);
	}
    }
    TH1D *hReco = new TH1D("hReco","hReco",nbins_jpt,binning_jpt);
    TH1D *hTruth= new TH1D("hTruth","hTruth",nbins_jpt,binning_jpt);
    TH2D *PT2D  = new TH2D("PT2D","PT2D",nbins_jpt,binning_jpt,nbins_jpt,binning_jpt);
    TH2D *Corr2D  = new TH2D("Corr2D","Corr2D",100,-20,20,10,0,100);
    TH1D *hTruth_R= new TH1D("hTruth_R","hTruth_R",nbins_jpt,binning_jpt);

    TFile *fin = new TFile("./Response_NoRecoPtCut.root","READ");
    TH1F * hMCJetPt = (TH1F*)fin->Get("MCRecoJetMatcher/hMCJetPt"); 
    TH1F * hRecoJetPt = (TH1F*)fin->Get("MCRecoJetMatcher/hRecoJetPt");
    TH2F * hMCvRecoJetPt = (TH2F*)fin->Get("MCRecoJetMatcher/hRecovMCJetPt");
    hMCvRecoJetPt->SetName("hMCvRecoJetPt");
    //This part is constructing a pT distribution in each GEN bin of which I can randomly sample later to construct a high-stat response matrix.
    // This really is just a rebining of the provided 2D response to match analysis binning
    TH1D *hSample[nbins_jpt];
    cout <<"> Now doing projections for sampling histograms "<<endl;
    for(int i = 0;i<nbins_jpt;i++){
	if(binning_jpt[i] < hMCvRecoJetPt->GetYaxis()->GetBinLowEdge(1)){
	    hSample[i]= new TH1D(Form("hSample_%i",i),Form("hSample_%i",i),nbins_jpt,binning_jpt);
	}
	else{
	    hMCvRecoJetPt->GetYaxis()->SetRangeUser(binning_jpt[i],binning_jpt[i+1]);
	    hSample[i]= (TH1D*)hMCvRecoJetPt->ProjectionX(Form("hSample_%i",i));
	}
    }    
    hMCvRecoJetPt->GetYaxis()->SetRangeUser(0,100);
    
    TFile *f_D = new TFile("/gpfs01/star/pwg/mkelsey/D0JetAuAu/Production/output/output_1/all.root","READ");
    TChain* ch = (TChain*) f_D->Get("JetTreeEmbed/Embed");
    float gpt; float rpt;float  mCen;float jpt;
    ch->SetBranchAddress( "PionPt" , &gpt );// Both Pion and Kaon variables have the generated PT
    ch->SetBranchAddress( "JetCorrPt" , &rpt );// jet PT - rho*A
    ch->SetBranchAddress( "JetPt" , &jpt );// jet PT - rho*A   
    ch->SetBranchAddress( "Centrality" , &mCen );//centrality (0,10,20,etc)
    int loop = ch->GetEntries()*1;
        
    TH3D *PT3D  = new TH3D("PT3D","PT3D",100,-20,20,100,-20,20,10,0,100);
    cout <<"> Filling response matrix " << endl;
    for(int i =0;i<loop;i++){
        ch->GetEntry(i);
        PT3D->Fill(rpt-gpt,gpt,mCen+0.1);
	Corr2D->Fill(jpt-rpt,mCen+0.1);
    }

    PT3D->GetZaxis()->SetRangeUser(cen_high,cen_low);
    TH2D * temp = (TH2D*)PT3D->Project3D("yx");
    TH1D *h1 = temp->ProjectionX();
    TF1 *fit = new TF1("fit","gaus(0)",-10,10);
    fit->SetParameter(2,5);
    h1->Fit(fit,"R");
    double rms = fit->GetParameter(2);

    cout << "RMS OF PROJECTION " << rms << endl;
    loop = 1000000;
    RooUnfoldResponse response(hReco,hTruth);
    for(int i =0;i<loop;i++){
        
	// Use gRandom for FONLL only
	//double gpt1 = gRandom->Uniform(0,40);
	double gpt1 = hMCJetPt->GetRandom();
	int bin = hTruth->FindBin(gpt1);
	double response_fac = hSample[bin-1]->GetRandom() - gpt1; //(Det. RECO-GEN)
	if(method>20)response_fac=0;
	double rpt1 =  gpt1 + response_fac + rms*gRandom->Gaus(0,1);
	//double rpt1 =  gpt1 + rms*gRandom->Gaus(0,1);

	//Don't need weights if I am sampleing from PYTHIA distribtuion
	//double weight = FONLL->Eval(gpt1)/FONLL->Integral(0,40); 
	double weight = 1;//hMCJetPt->GetBinContent(hMCJetPt->FindBin(gpt1))/hMCJetPt->Integral();
	if(gpt1>gen_low){
	    // reco pt cutt-off so add inefficiency part here
	    if(reco_low>=0){
		hTruth->Fill(gpt1,weight);                                                                                                                                                                                                    
		if(rpt1>reco_low){
		    PT2D->Fill(rpt1,gpt1,weight);
		    hReco->Fill(rpt1,weight);
		    response.Fill(rpt1,gpt1,weight);
		    hTruth_R->Fill(gpt1,weight); 
		}
		else{
		    response.Miss(gpt1,weight);
		}
	    }
	    // no low pt cut-off
	    else{ 
		PT2D->Fill(rpt1,gpt1,weight);
	        response.Fill(rpt1,gpt1,weight);
	    }
	}
	else{
	    cout <<"THIS SHOULD NOT HAPPEN" << endl;
	}
    }
    cout <<"> Done filling response matrix " << endl;
    f_D->Close();
    //RooUnfoldResponse response (hReco,hTruth,PT2D);//,"response","response"); 

    TMatrixD matx = response.Mresponse();
    TH2D *res_matx = new TH2D(matx);
    // ============== Done with response matrix
    // ============== Unfold

    cout << "> Starting the unfolding part " << endl;
    RooUnfoldBayes unfold (&response, data, index);
    //if(method==2)RooUnfoldSvd unfold (&response, data, index);
    TH1D* hData = (TH1D*) unfold.Hreco();

    cout << "> Finished with unfolding part " << endl;

    // =============== Now draw

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
    if(save && !sim){
	char out[100];
	sprintf(out,"hists/Unfolded_Hist_Cen_%i_%i_Method_%i_Indx_%i.root",cen_high,cen_low,method,index);
	TFile *fout = new TFile(out,"RECREATE");
	char hout[100];
	sprintf(hout,"hData_%i_%i_%i_%i",cen_high,cen_low,method,index);
	hData->Write(hout);
	for(int i = 0;i<nbins_jpt;i++)hSample[i]->Write(Form("hSample_%i",i));
	hMCvRecoJetPt->Write();
	hRecoJetPt->Write();
	hMCJetPt->Write();
	PT2D->Write();
	hReco->Write();
	hTruth->Write();
	hTruth_R->Write();
	fout->Close();
    }
    else if(save && sim){
	char out[100];
        sprintf(out,"hists/ToyUnfolded_Hist_Cen_%i_%i_Method_%i_Indx_%i.root",cen_high,cen_low,method,index);
        TFile *fout = new TFile(out,"RECREATE");
        char hout[100];
        sprintf(hout,"hData_%i_%i_%i_%i",cen_high,cen_low,method,index);
        hData->Write(hout);
	hReco->Write();
        hTruth->Write();
        hTruth_R->Write();
	PT2D->Write();
	hMCvRecoJetPt->Write();
        hRecoJetPt->Write();
        hMCJetPt->Write();
        fout->Close();
	
    }
}

void getData(TGraph *g,char file[]){
    int cnt=0;
    ifstream data(file);
    if(data.is_open()){
        while(!data.eof()){
            double x;double y;
            data >> x >> y;
            g->SetPoint(cnt,x,y);cnt++;
        }
    }
}
