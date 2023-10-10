void PhiRes(){
    gROOT->ProcessLine(".x ~/myStyle.C");
    TFile* outfile = new TFile("BKG_Fluc_CorrelationsNew.root");
    PT2D = (TH3F*)outfile->Get("Corr3D_phi");
    temp1 = (TH2F*) PT2D->Project3D("yx");
    temp1->SetName("temp1");
    PT2D->GetZaxis()->SetRangeUser(0,5);
    temp_1 = (TH2F*) PT2D->Project3D("xy");temp_1->SetName("temp_1");
    PT2D->GetZaxis()->SetRangeUser(5,10);
    temp_2 = (TH2F*) PT2D->Project3D("xy");temp_2->SetName("temp_2");
    PT2D->GetZaxis()->SetRangeUser(10,20);
    temp_3 = (TH2F*) PT2D->Project3D("xy");temp_3->SetName("temp_3");
    PT2D->GetZaxis()->SetRangeUser(20,30);
    temp_4 = (TH2F*) PT2D->Project3D("xy");temp_4->SetName("temp_4");
    PT2D->GetZaxis()->SetRangeUser(30,40);
    temp_5 = (TH2F*) PT2D->Project3D("xy");temp_5->SetName("temp_5");
    PT2D->GetZaxis()->SetRangeUser(40,50);
    temp_6 = (TH2F*) PT2D->Project3D("xy");temp_6->SetName("temp_6");
    PT2D->GetZaxis()->SetRangeUser(50,60);
    temp_7 = (TH2F*) PT2D->Project3D("xy");temp_7->SetName("temp_7");
    PT2D->GetZaxis()->SetRangeUser(60,70);
    temp_8 = (TH2F*) PT2D->Project3D("xy");temp_8->SetName("temp_8");
    PT2D->GetZaxis()->SetRangeUser(70,80);
    temp_9 = (TH2F*) PT2D->Project3D("xy");temp_9->SetName("temp_9");



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


    temp_6->Scale(w6);
    temp_6->Add(temp_7,w7);
    temp_6->Add(temp_8,w8);
    temp_6->Add(temp_9,w9);

    temp_3->Scale(w3);
    temp_3->Add(temp_4,w4);
    temp_3->Add(temp_5,w5);

    temp_1->Scale(w1);
    temp_1->Add(temp_2,w2);


    TCanvas *c1 = new TCanvas("c1","c1");
    temp1->Draw("COLZ");
    
    TCanvas *c2 = new TCanvas("c2","c2");
    TF1 *fit = new TF1("fit","gaus(0)",-70,70);
    fit->SetParameter(2,5);

    TH1F *eta_sigma_pt_0_10 = new TH1F("eta_sigma_pt_0_10","eta_sigma_pt_0_10",40,0,40);
    TH1F *eta_sigma_pt_10_40 = new TH1F("eta_sigma_pt_10_40","eta_sigma_pt_10_40",40,0,40);
    TH1F *eta_sigma_pt_40_80 = new TH1F("eta_sigma_pt_40_80","eta_sigma_pt_40_80",40,0,40);






    for(int i = 0;i<eta_sigma_pt_0_10->GetNbinsX()+1;i++){
	
	int bin1= temp_1->GetYaxis()->FindBin(eta_sigma_pt_0_10->GetBinLowEdge(i)+0.0001);
	int bin2= temp_1->GetYaxis()->FindBin(eta_sigma_pt_0_10->GetBinLowEdge(i)+eta_sigma_pt_0_10->GetBinWidth(i));

	TH1F* h1 = (TH1F*)temp_1->ProjectionX("h1",i,i);
	eta_sigma_pt_0_10->SetBinContent(i,h1->GetRMS());
	cout << eta_sigma_pt_0_10->GetBinLowEdge(i) << " " << bin1 << " " << bin2 << " " << h1->GetRMS() << " " << temp_1->GetYaxis()->GetBinLowEdge(bin1) << " " << temp_1->GetYaxis()->GetBinLowEdge(bin2)+ temp_1->GetYaxis()->GetBinWidth(bin2) <<endl;

	TH1F* h1 = (TH1F*)temp_3->ProjectionX("h1",i,i);
	eta_sigma_pt_10_40->SetBinContent(i,h1->GetRMS());
	
	TH1F* h1 = (TH1F*)temp_6->ProjectionX("h1",i,i);
	eta_sigma_pt_40_80->SetBinContent(i,h1->GetRMS());


    }
    eta_sigma_pt_0_10->GetXaxis()->SetTitle("Embedded Jet p_{T} (GeV)");
    eta_sigma_pt_0_10->GetYaxis()->SetTitle("#eta Resolution");
    eta_sigma_pt_0_10->SetLineColor(kRed);
    eta_sigma_pt_10_40->SetLineColor(kBlue);
    eta_sigma_pt_40_80->SetLineColor(kGreen-2);

    TF1 *eta_fit_1 = new TF1("eta_fit_1","[0]+TMath::Power(x,[1])",0,50);
    TF1 *eta_fit_2 = new TF1("eta_fit_2","[0]+TMath::Power(x,[1])",0,50);
    TF1 *eta_fit_3 = new TF1("eta_fit_3","[0]+TMath::Power(x,[1])",0,50);
    eta_sigma_pt_0_10->Fit(eta_fit_1);
    eta_sigma_pt_10_40->Fit(eta_fit_2);
    eta_sigma_pt_40_80->Fit(eta_fit_3);

    TCanvas *c33 = new TCanvas("c33","c33");
    eta_sigma_pt_0_10->Draw("hist");
    eta_sigma_pt_10_40->Draw("hist same");
    eta_sigma_pt_40_80->Draw("hist same");
    eta_fit_1->Draw("same");
    eta_fit_2->Draw("same");
    eta_fit_3->Draw("same");
    PT2D->GetZaxis()->SetRangeUser(40,80);
    hh1 =  (TH1F*)PT2D->Project3D("y");hh1->SetName("hh1");
    PT2D->GetZaxis()->SetRangeUser(10,40);
    hh2 =  (TH1F*)PT2D->Project3D("y");hh2->SetName("hh2");
    PT2D->GetZaxis()->SetRangeUser(0,10);
    hh3 =  (TH1F*)PT2D->Project3D("y");hh3->SetName("hh3");

    TCanvas *c3 = new TCanvas("c3","c3");
    hh1->DrawNormalized("PE");
    hh2->SetLineColor(kBlue);
    hh2->DrawNormalized("same  hist");
    hh3->SetLineColor(kRed);
    hh3->DrawNormalized("same  hist ");
    TFile *fout = new TFile("Phi_Res_FINAL.root","RECREATE");
    eta_sigma_pt_0_10->Write("phi_sigma_pt_0_10");
    eta_sigma_pt_10_40->Write("phi_sigma_pt_10_40");
    eta_sigma_pt_40_80->Write("phi_sigma_pt_40_80");
    eta_fit_1->Write("phi_fit_1");;
    eta_fit_2->Write("phi_fit_2");
    eta_fit_3->Write("phi_fit_3");

}
