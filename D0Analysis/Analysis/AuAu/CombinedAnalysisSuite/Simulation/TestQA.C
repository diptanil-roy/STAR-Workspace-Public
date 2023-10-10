R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

#include "BinDef.h"
#include "NewBinDef.h"

void TestQA(){
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TH2D *Measured[3][3];
    TH1D *MeasuredPt[3][3];
    TH1D *MeasuredZ[3][3];

    TH2D *Truth[3][3];
    TH1D *TruthPt[3][3];
    TH1D *TruthZ[3][3];

    TH2D *Unfolded[3][3];
    TH1D *UnfoldedPt[3][3];
    TH1D *UnfoldedZ[3][3];
    TH1D *Unfolded1D[3][3];
    TH1D *Unfolded1DdR[3][3];

    TFile *FONLL = new TFile("Aug14_FONLL_1GeV_Data4/Output_Step_0_IterParam_0_njptbin_12_nzbin_7.root");
    FONLL->cd();
    for (int cent = 0; cent < 3; cent++){
        Measured[0][cent] = (TH2D *)gDirectory->Get(Form("Measured_%i", cent));
        SetName(Measured[0][cent], Form("%s %s", Measured[0][cent]->GetName(), "FONLL"));
        Measured[0][cent]->Scale(1./Measured[0][cent]->Integral());
        MeasuredPt[0][cent] = (TH1D *)Measured[0][cent]->ProjectionX();
        MeasuredZ[0][cent] = (TH1D *)Measured[0][cent]->ProjectionY();

        Truth[0][cent] = (TH2D *)gDirectory->Get(Form("fTruthWide_Cent_%i", cent));
    }

    TFile *FONLLData = new TFile("Aug14_FONLL_1GeV_Data_WeighedByData4/Response_Step_0_IterParam_0_njptbin_12_nzbin_7.root");
    FONLLData->cd();
    for (int cent = 0; cent < 3; cent++){
        Measured[1][cent] = (TH2D *)gDirectory->Get(Form("Measured_%i", cent));
        SetName(Measured[1][cent], Form("%s %s", Measured[1][cent]->GetName(), "FONLL+Data"));
        Measured[1][cent]->Scale(1./Measured[1][cent]->Integral());
        MeasuredPt[1][cent] = (TH1D *)Measured[1][cent]->ProjectionX();
        MeasuredZ[1][cent] = (TH1D *)Measured[1][cent]->ProjectionY();
    }

    // TFile *PYTHIAPowerLawFit = new TFile("Aug14_PYTHIA_Corrected_1GeV_Data4/Response_Step_0_IterParam_0_njptbin_12_nzbin_7.root");
    // PYTHIAPowerLawFit->cd();
    // for (int cent = 0; cent < 3; cent++){
    //     Measured[2][cent] = (TH2D *)gDirectory->Get(Form("Measured_%i", cent));
    //     SetName(Measured[2][cent], Form("%s %s", Measured[2][cent]->GetName(), "PYTHIA+PowerLawFit"));
    //     Measured[2][cent]->Scale(1./Measured[2][cent]->Integral());
    //     MeasuredPt[2][cent] = (TH1D *)Measured[2][cent]->ProjectionX();
    //     MeasuredZ[2][cent] = (TH1D *)Measured[2][cent]->ProjectionY();
    // }

    TCanvas *c = new TCanvas("c", "c", 1200, 1200);
    c->Divide(3);

    for (int cent = 0; cent < 3; cent++){
        c->cd(cent+1);
        gPad->SetLogy();
        SetColor(MeasuredPt[0][cent], kBlack);
        SetColor(MeasuredPt[1][cent], kRed);
        SetColor(MeasuredPt[2][cent], kBlue);
        MeasuredPt[0][cent]->Draw();
        MeasuredPt[1][cent]->Draw("same");
        MeasuredPt[2][cent]->Draw("same");
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 1200);
    c2->Divide(3);

    for (int cent = 0; cent < 3; cent++){
        c2->cd(cent+1);
        gPad->SetLogy();
        MeasuredZ[0][cent]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        SetColor(MeasuredZ[0][cent], kBlack);
        SetColor(MeasuredZ[1][cent], kRed);
        SetColor(MeasuredZ[2][cent], kBlue);
        MeasuredZ[0][cent]->Draw();
        MeasuredZ[1][cent]->Draw("same");
        MeasuredZ[2][cent]->Draw("same");
    }
    
}