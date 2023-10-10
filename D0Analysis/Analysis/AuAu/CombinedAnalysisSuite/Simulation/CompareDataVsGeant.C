using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void CompareDataVsGeant(TString BaseFolderName = "", int D0pT = 1, TString outputdir = "", TString extratagforoutput = ""){
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    if (extratagforoutput != "") extratagforoutput = "_" + extratagforoutput;

    TFile *geant = new TFile(Form("%s/Response_Step_0_IterParam_0_njptbin_%i_nzbin_%i.root", BaseFolderName.Data(), njpt_gen_bins_var, nz_gen_bins));
    geant->cd();
    // geant->ls();
    THnSparseF *hJetWide[3];
    for (int i = 0; i < 3; i++){
        hJetWide[i] = (THnSparseF *)gDirectory->Get(Form("hJetWide_%i", i));
        hJetWide[i]->SetName(Form("hJetWide_%i", i));
        // hJetWide[i]->SetDirectory(0);
    }

    TFile *data = new TFile(Form("../Data/Aug17_2023/Histograms_D0%i_10GeV_RecoJetPt_%i_1000.root", D0pT, D0pT));
	data->cd();

	TH2D *MeasuredWide[3];

	MeasuredWide[0] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_0_10");
	MeasuredWide[1] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_10_40");
	MeasuredWide[2] = (TH2D *)gDirectory->Get("ZPt_Area_Wide_40_80");

    for (int i = 0; i < 3; i++){
        for (int binx = 0; binx < MeasuredWide[i]->GetNbinsX(); binx++){
            for (int biny = 0; biny < MeasuredWide[i]->GetNbinsY(); biny++){
                if (MeasuredWide[i]->GetBinContent(binx+1, biny+1) < 0){
                    MeasuredWide[i]->SetBinContent(binx+1, biny+1, 0);
                }
            }
        }
    }

    TH1D *hMCPt[3];
    TH1D *hMCZ[3];
    TH1D *hRecoPt[3];
    TH1D *hRecoZ[3];
    TH2D *hReco[3];
    TH1D *hDataPt[3];
    TH1D *hDataZ[3];
    TH1D *hRatioGeantVsDataPt[3];
    TH1D *hRatioGeantVsDataZ[3];
    TH1D *hRatioGeantVsData[3];

    for (int i = 0; i < 3; i++){
        hMCPt[i] = (TH1D *)hJetWide[i]->Projection(0);
        SetName(hMCPt[i], Form("hMCPt_%i", i));
        hMCZ[i] = (TH1D *)hJetWide[i]->Projection(1);
        SetName(hMCZ[i], Form("hMCZ_%i", i));
        hRecoPt[i] = (TH1D *)hJetWide[i]->Projection(2);
        SetName(hRecoPt[i], Form("hRecoPt_%i", i));
        hRecoZ[i] = (TH1D *)hJetWide[i]->Projection(3);
        SetName(hRecoZ[i], Form("hRecoZ_%i", i));
        hReco[i] = (TH2D *)hJetWide[i]->Projection(3, 2);
        SetName(hReco[i], Form("hReco_%i", i));
        hDataPt[i] = (TH1D *)MeasuredWide[i]->ProjectionX();
        SetName(hDataPt[i], Form("hDataPt_%i", i));
        hDataZ[i] = (TH1D *)MeasuredWide[i]->ProjectionY();
        SetName(hDataZ[i], Form("hDataZ_%i", i));
    }

    for (int i = 0; i < 3; i++){
        hMCPt[i]->Scale(1./hMCPt[i]->Integral());
        SetColor(hMCPt[i], kRed);
        hMCZ[i]->Scale(1./hMCZ[i]->Integral());
        SetColor(hMCZ[i], kRed);
        hRecoPt[i]->Scale(1./hRecoPt[i]->Integral());
        SetColor(hRecoPt[i], kBlue);
        hRecoZ[i]->Scale(1./hRecoZ[i]->Integral());
        SetColor(hRecoZ[i], kBlue);
        hReco[i]->Scale(1./hReco[i]->Integral());
        hDataPt[i]->Scale(1./hDataPt[i]->Integral());
        SetColor(hDataPt[i], kBlack);
        hDataZ[i]->Scale(1./hDataZ[i]->Integral());
        SetColor(hDataZ[i], kBlack);

        hReco[i]->Scale(1./hReco[i]->Integral());
        MeasuredWide[i]->Scale(1./MeasuredWide[i]->Integral());
    }

    for (int i = 0; i < 3; i++){
        hRatioGeantVsDataPt[i] = (TH1D *)hDataPt[i]->Clone(Form("hRatioGeantVsDataPt_%i", i));
        hRatioGeantVsDataPt[i]->Divide(hRecoPt[i]);
        hRatioGeantVsDataZ[i] = (TH1D *)hDataZ[i]->Clone(Form("hRatioGeantVsDataZ_%i", i));
        hRatioGeantVsDataZ[i]->Divide(hRecoZ[i]);
        hRatioGeantVsData[i] = (TH1D *)MeasuredWide[i]->Clone(Form("hRatioGeantVsData_%i", i));
        hRatioGeantVsData[i]->Divide(hReco[i]);
    }

    TCanvas *cPt = new TCanvas("cPt", "cPt", 1200, 600);
    cPt->Divide(3, 1);

    for (int i = 0; i < 3; i++){
        cPt->cd(i+1);
        gPad->SetLogy();
        // hMCPt[i]->Draw();
        DivideByBinWidth(hMCPt[i]);
        DivideByBinWidth(hRecoPt[i]);
        DivideByBinWidth(hDataPt[i]);
        hRecoPt[i]->Draw("same");
        hDataPt[i]->Draw("same");
        // hRatioGeantVsDataPt[i]->GetYaxis()->SetRangeUser(0, 2);
        // hRatioGeantVsDataPt[i]->Draw();
    }

    cPt->SaveAs(Form("%s/GEANTvData_pT_%i%s.pdf", outputdir.Data(), D0pT, extratagforoutput.Data()));

    TCanvas *cZ = new TCanvas("cZ", "cZ", 1200, 600);
    cZ->Divide(3, 1);

    for (int i = 0; i < 3; i++){
        cZ->cd(i+1);
        gPad->SetLogy();
        // hMCZ[i]->Draw();
        DivideByBinWidth(hMCZ[i]);
        DivideByBinWidth(hRecoZ[i]);
        DivideByBinWidth(hDataZ[i]);
        // hRecoZ[i]->GetXaxis()->SetRangeUser(-2, 10);
        hRecoZ[i]->Draw("same");
        hDataZ[i]->Draw("same");
        // hRatioGeantVsDataZ[i]->GetYaxis()->SetRangeUser(0, 5);
        // hRatioGeantVsDataZ[i]->GetXaxis()->SetRangeUser(-10, 10);
        // hRatioGeantVsDataZ[i]->Draw();
    }

    cZ->SaveAs(Form("%s/GEANTvData_Z_%i%s.pdf", outputdir.Data(), D0pT, extratagforoutput.Data()));

    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    c->Divide(3, 1);

    for (int i = 0; i < 3; i++){
        c->cd(i+1);
        hRatioGeantVsData[i]->Draw("colz");
    }

    TCanvas *dPt = new TCanvas("dPt", "dPt", 1200, 600);
    dPt->Divide(3, 1);
    for (int i = 0; i < 3; i++){
        dPt->cd(i+1);
        // gPad->SetLogy();

        // hRatioGeantVsDataPt[i]->GetYaxis()->SetRangeUser(0, 2);
        hRatioGeantVsDataPt[i]->Draw();
    }

    dPt->SaveAs(Form("%s/GEANTvData_pTRatio_%i%s.pdf", outputdir.Data(), D0pT, extratagforoutput.Data()));

    TCanvas *dZ = new TCanvas("dZ", "dZ", 1200, 600);
    dZ->Divide(3, 1);

    for (int i = 0; i < 3; i++){
        dZ->cd(i+1);
        // gPad->SetLogy();

        // hRatioGeantVsDataZ[i]->GetYaxis()->SetRangeUser(0, 5);
        // hRatioGeantVsDataZ[i]->GetXaxis()->SetRangeUser(-10, 10);
        hRatioGeantVsDataZ[i]->Draw();
    }

    dZ->SaveAs(Form("%s/GEANTvData_ZRatio_%i%s.pdf", outputdir.Data(), D0pT, extratagforoutput.Data()));


    TFile *WeightFile = new TFile(Form("%s/GEANTvDATA_%i%s.root", outputdir.Data(), D0pT, extratagforoutput.Data()), "RECREATE");
    cout << WeightFile->GetName() << endl;
    WeightFile->cd();
    for (int i = 0; i < 3; i++){
        SetName(hRatioGeantVsDataPt[i], Form("DataWeightPt_%i", i));
        SetName(hRatioGeantVsDataZ[i], Form("DataWeightZ_%i", i));
        SetName(hRatioGeantVsData[i], Form("DataWeight_%i", i));
        hRatioGeantVsDataPt[i]->Write();
        hRatioGeantVsDataZ[i]->Write();
        hRatioGeantVsData[i]->Write();
    }
    WeightFile->Close();
}