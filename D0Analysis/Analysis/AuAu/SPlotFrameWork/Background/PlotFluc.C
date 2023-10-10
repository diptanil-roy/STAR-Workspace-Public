void PlotFluc(){
    gROOT->ProcessLine(".x ~/myStyle.C");
    TFile* outfile = new TFile("BKG_Fluc_Correlations.root");
    PT2D = (TH3F*)outfile->Get("PT3D");
    Corr2D = (TH2F*)outfile->Get("Corr2D");
    temp2 = (TH2F*) PT2D->Project3D("yz");

    temp1 = (TH2F*) PT2D->Project3D("zx");
    temp1->SetName("temp1");

    PT2D->GetYaxis()->SetRangeUser(0,100);
    
    PT2D->GetZaxis()->SetRangeUser(0,9);
    hBKG_0_10 = (TH2F*) PT2D->Project3D("x");hBKG_0_10->SetName("hBKG_0_10");
    PT2D->GetZaxis()->SetRangeUser(10,39);
    hBKG_10_40 = (TH2F*) PT2D->Project3D("x");hBKG_10_40->SetName("hBKG_10_40");
    PT2D->GetZaxis()->SetRangeUser(40,80);
    hBKG_40_80 = (TH2F*) PT2D->Project3D("x");hBKG_40_80->SetName("hBKG_40_80");
    

    PT2D->GetZaxis()->SetRangeUser(0,5);
    hBKG_0_5 = (TH1F*) PT2D->Project3D("x");hBKG_0_5->SetName("hBKG_0_5");
    PT2D->GetZaxis()->SetRangeUser(5,10);
    hBKG_5_10 = (TH1F*) PT2D->Project3D("x");hBKG_5_10->SetName("hBKG_5_10");
    PT2D->GetZaxis()->SetRangeUser(10,20);
    hBKG_10_20 = (TH1F*) PT2D->Project3D("x");hBKG_10_20->SetName("hBKG_10_20");
    PT2D->GetZaxis()->SetRangeUser(20,30);
    hBKG_20_30 = (TH1F*) PT2D->Project3D("x");hBKG_20_30->SetName("hBKG_20_30");
    PT2D->GetZaxis()->SetRangeUser(30,40);
    hBKG_30_40 = (TH1F*) PT2D->Project3D("x");hBKG_30_40->SetName("hBKG_30_40");
    PT2D->GetZaxis()->SetRangeUser(40,50);
    hBKG_40_50 = (TH1F*) PT2D->Project3D("x");hBKG_40_50->SetName("hBKG_40_50");
    PT2D->GetZaxis()->SetRangeUser(50,60);
    hBKG_50_60 = (TH1F*) PT2D->Project3D("x");hBKG_50_60->SetName("hBKG_40_60");
    PT2D->GetZaxis()->SetRangeUser(60,70);
    hBKG_60_70 = (TH1F*) PT2D->Project3D("x");hBKG_60_70->SetName("hBKG_50_70");
    PT2D->GetZaxis()->SetRangeUser(70,80);
    hBKG_70_80 = (TH1F*) PT2D->Project3D("x");hBKG_70_80->SetName("hBKG_70_80");

    TFile* fin = new TFile("../ApplyWeights/Histograms3_D05GeV.root");
    centrality = (TH1F*)fin->Get("centrality");
    double val9 = centrality->GetBinContent(centrality->FindBin(70));
    double val8 = centrality->GetBinContent(centrality->FindBin(60));
    double val7 = centrality->GetBinContent(centrality->FindBin(50));
    double val6 = centrality->GetBinContent(centrality->FindBin(40));
    double val5 = centrality->GetBinContent(centrality->FindBin(30));
    double val4 = centrality->GetBinContent(centrality->FindBin(20));
    double val3 = centrality->GetBinContent(centrality->FindBin(10));
    double val2 = centrality->GetBinContent(centrality->FindBin(5));
    double val1 = centrality->GetBinContent(centrality->FindBin(0));    

    double sval9 = hBKG_70_80->Integral();
    double sval8 = hBKG_60_70->Integral();
    double sval7 = hBKG_50_60->Integral();
    double sval6 = hBKG_40_50->Integral();
    double sval5 = hBKG_30_40->Integral();
    double sval4 = hBKG_20_30->Integral();
    double sval3 = hBKG_10_20->Integral();
    double sval2 = hBKG_5_10->Integral();
    double sval1 = hBKG_0_5->Integral();



    cout << "w9 " << val9<< endl;
    cout << "w8 " << val8<< endl;
    cout << "w7 " << val7<< endl;
    cout << "w6 " << val6<< endl;
    cout << "w5 " << val5<< endl;
    cout << "w4 " << val4<< endl;
    cout << "w3 " << val3<< endl;
    cout << "w2 " << val2<< endl;
    cout << "w1 " << val1<< endl;




    double norm = val9+val8+val7+val6;
    double snorm = sval9+sval8+sval7+sval6;
    double w9 = (val9/norm)/(sval9/snorm);
    double w8 = (val8/norm)/(sval8/snorm);
    double w7 = (val7/norm)/(sval7/snorm);
    double w6 = (val6/norm)/(sval6/snorm);
   

    double norm= val3+val4+val5;
    double snorm = sval3+sval4+sval5;
    double w5 = (val5/norm)/(sval5/snorm);
    double w4 = (val4/norm)/(sval4/snorm);
    double w3 = (val3/norm)/(sval3/snorm);
    
    double norm= val1+val2;
    double snorm = sval1+sval2;
    double w2 = (val2/norm)/(sval2/snorm);
    double w1 = (val1/norm)/(sval1/snorm);

    cout << "w9 " << w9<< endl;
    cout << "w8 " << w8<< endl;
    cout << "w7 " << w7<< endl;
    cout << "w6 " << w6<< endl;
    cout << "w5 " << w5<< endl;
    cout << "w4 " << w4<< endl;
    cout << "w3 " << w3<< endl;
    cout << "w2 " << w2<< endl;
    cout << "w1 " << w1<< endl;


    hBKG_40_50->Scale(w6);
    hBKG_40_50->Add(hBKG_50_60,w7);
    hBKG_40_50->Add(hBKG_50_60,w8);
    hBKG_40_50->Add(hBKG_50_60,w9);

    hBKG_10_20->Scale(w3);
    hBKG_10_20->Add(hBKG_20_30,w4);
    hBKG_10_20->Add(hBKG_30_40,w5);

    hBKG_0_5->Scale(w1);
    hBKG_0_5->Add(hBKG_5_10,w2);

    TFile *save = new TFile("BKGHistogramsFINAL.root","RECREATE");
    hBKG_0_10->Write();
    hBKG_10_40->Write();
    hBKG_40_80->Write();
    hBKG_40_50->Write("hBKG_40_80_Weighted");
    hBKG_10_20->Write("hBKG_10_40_Weighted");
    hBKG_0_5->Write("hBKG_0_10_Weighted");


    PT2D->GetZaxis()->SetRangeUser(40,80);
    temp = (TH2F*) PT2D->Project3D("yx");
    PT2D->GetZaxis()->SetRangeUser(40,50);
    //PT2D->GetYaxis()->SetRangeUser(0,5);
    _temp = (TH2F*) PT2D->Project3D("x");
    PT2D->GetYaxis()->SetRangeUser(5,15);
    __temp = (TH2F*) PT2D->Project3D("x");

    TCanvas *c1 = new TCanvas("c1","c1");
    temp1->Draw("COLZ");
    
    TCanvas *c2 = new TCanvas("c2","c2");
    TF1 *fit = new TF1("fit","gaus(0)",-10,10);
    fit->SetParameter(2,1.7);

    TH1F* h1 = (TH1F*)temp->ProjectionX();
    h1->Fit(fit,"R");
    cout <<"RMS HERE " << h1->GetRMS() <<endl;

    TCanvas *c22 = new TCanvas("c2™","c22");
    h1->DrawNormalized("hist");
    hBKG_40_50->DrawNormalized("same PE");
//    _temp->DrawNormalized("same PE");
    //__temp->DrawNormalized("same PE");
    double rms = fit->GetParameter(2);

    double low = 0;
    double high = 0;
    double max = h1->GetMaximum();
    double _cnt=0;
    for(int i = 1;i<h1->GetNbinsX()+1;i++){
	double val = h1->GetBinContent(i);
	//cout << ">Calc " << max/2. << " " << val << " " << h1->GetBinLowEdge(i)  << " " << _cnt << endl; 
	if(_cnt==2)continue;
        if(_cnt== 1 &&val < max/2.){
            high = h1->GetBinLowEdge(i-1)+h1->GetBinWidth(i-1);
            _cnt++;
        }
	if(_cnt == 0 && val > max/2.){
	    low = h1->GetBinLowEdge(i);
	    _cnt++;
	} 
    }
    cout <<" FWHM " << high - low << " sigma " <<  (high - low)/2.35482  <<" low " << low << " high " << high << endl;


    cout << "RMS " << rms << endl;
    temp2->RebinY();
    temp2->GetXaxis()->SetRangeUser(40,80);
    hh1 =  (TH1F*)temp2->ProjectionY("hh1");
    temp2->GetXaxis()->SetRangeUser(10,40);
    hh2 =  (TH1F*)temp2->ProjectionY("hh2");
    temp2->GetXaxis()->SetRangeUser(0,10);
    hh3 =  (TH1F*)temp2->ProjectionY("hh3");

    hh1->GetXaxis()->SetRangeUser(0,20);
    TCanvas *c3 = new TCanvas("c3","c3");
    hh1->DrawNormalized("PE");
    hh2->SetLineColor(kBlue);
    hh2->DrawNormalized("same  hist");
    hh3->SetLineColor(kRed);
    hh3->DrawNormalized("same  hist ");
}
