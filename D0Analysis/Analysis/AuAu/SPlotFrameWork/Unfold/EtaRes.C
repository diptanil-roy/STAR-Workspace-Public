void EtaRes(){
    gROOT->ProcessLine(".x ~/myStyle.C");
    TFile* outfile = new TFile("BKG_Fluc_Correlations.root");
    PT2D = (TH3F*)outfile->Get("Corr3D_eta");
    temp1 = (TH2F*) PT2D->Project3D("yx");
    temp1->SetName("temp1");
    PT2D->GetZaxis()->SetRangeUser(0,10);
    temp_1 = (TH2F*) PT2D->Project3D("xy");temp_1->SetName("temp_1");
    PT2D->GetZaxis()->SetRangeUser(10,40);
    temp_2 = (TH2F*) PT2D->Project3D("xy");temp_2->SetName("temp_2");
    PT2D->GetZaxis()->SetRangeUser(40,80);
    temp_3 = (TH2F*) PT2D->Project3D("xy");temp_3->SetName("temp_3");

    TCanvas *c1 = new TCanvas("c1","c1");
    temp1->Draw("COLZ");
    
    TCanvas *c2 = new TCanvas("c2","c2");
    TF1 *fit = new TF1("fit","gaus(0)",-70,70);
    fit->SetParameter(2,5);

    TH1F *eta_sigma_pt_0_10 = new TH1F("eta_sigma_pt_0_10","eta_sigma_pt_0_10",40,0,40);
    TH1F *eta_sigma_pt_10_40 = new TH1F("eta_sigma_pt_10_40","eta_sigma_pt_10_40",40,0,40);
    TH1F *eta_sigma_pt_40_80 = new TH1F("eta_sigma_pt_40_80","eta_sigma_pt_40_80",40,0,40);
    for(int i = 0;i<eta_sigma_pt_0_10->GetNbinsX()+1;i++){
	
	int bin1= temp_1->GetYaxis()->FindBin(eta_sigma_pt_0_10->GetBinLowEdge(i));
	int bin2= temp_1->GetYaxis()->FindBin(eta_sigma_pt_0_10->GetBinLowEdge(i)+eta_sigma_pt_0_10->GetBinWidth(i));
	TH1F* h1 = (TH1F*)temp_1->ProjectionX("h1",bin1,bin2);
	eta_sigma_pt_0_10->SetBinContent(i,h1->GetRMS());
    
	TH1F* h1 = (TH1F*)temp_2->ProjectionX("h1",bin1,bin2);
	eta_sigma_pt_10_40->SetBinContent(i,h1->GetRMS());
	
	TH1F* h1 = (TH1F*)temp_3->ProjectionX("h1",bin1,bin2);
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

    TFile *fout = new TFile("Eta_Res_FINAL.root","RECREATE");
    eta_sigma_pt_0_10->Write();
    eta_sigma_pt_10_40->Write();
    eta_sigma_pt_40_80->Write();
    eta_fit_1->Write();;
    eta_fit_2->Write();
    eta_fit_3->Write();
}
