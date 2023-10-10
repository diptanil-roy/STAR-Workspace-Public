using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

bool CompareHistos(TH1 *h1, TH1 *h2){
    if (h1->GetNbinsX() != h2->GetNbinsX()){
        cout << "ERROR: Histograms have different number of bins" << endl;
        return false;
    }
    for (int i = 0; i < h1->GetNbinsX()+1; i++){
        if (h1->GetBinContent(i) < h2->GetBinContent(i)) return false;
    }
    return true;
}

bool CompareHistos(TH2 *h1, TH2 *h2){
    if (h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY()){
        cout << "ERROR: Histograms have different number of bins" << endl;
        return false;
    }
    for (int i = 0; i < h1->GetNbinsX()+1; i++){
        for (int j = 0; j < h1->GetNbinsY()+1; j++){
            if (h1->GetBinContent(i, j) < h2->GetBinContent(i, j)) return false;
        }
    }
    return true;
}

void TmpChecker(){
    TH2D *Measured[5][3];
    TH1D *MeasuredPt[5][3];
    TH1D *MeasuredZ[5][3];

    TH2D *Unfolded[5][3];
    TH1D *UnfoldedPt[5][3];
    TH1D *UnfoldedZ[5][3];

    int colors[5] = {kBlack, kRed, kBlue, kGreen, kMagenta};

    for (int i = 0; i < 5; i++){
        for (int cent = 0; cent < 3; cent++){
            TFile *f = new TFile(Form("Jun25_%iGeV_Data%i/Output_Step_%i_IterParam_%i_njptbin_%i_nzbin_%i.root", i+1, 4, 0, 0, njpt_gen_bins_var, nz_gen_bins));
            f->cd();
            Measured[i][cent] = (TH2D *)gDirectory->Get(Form("MeasuredWide_%i", cent));
            MeasuredPt[i][cent] = (TH1D *)gDirectory->Get(Form("MeasuredPtWide_%i", cent));
            MeasuredZ[i][cent] = (TH1D *)gDirectory->Get(Form("MeasuredZWide_%i", cent));
            SetName(Measured[i][cent], Form("%s_%iGeV", Measured[i][cent]->GetName(), i+1));
            SetName(MeasuredPt[i][cent], Form("%s_%iGeV", MeasuredPt[i][cent]->GetName(), i+1));
            SetName(MeasuredZ[i][cent], Form("%s_%iGeV", MeasuredZ[i][cent]->GetName(), i+1));

            Unfolded[i][cent] = (TH2D *)gDirectory->Get(Form("UnfoldedWide_%i", cent));
            UnfoldedPt[i][cent] = (TH1D *)gDirectory->Get(Form("UnfoldedPtWide_%i", cent));
            UnfoldedZ[i][cent] = (TH1D *)gDirectory->Get(Form("UnfoldedZWide_%i", cent));
            SetName(Unfolded[i][cent], Form("%s_%iGeV", Unfolded[i][cent]->GetName(), i+1));
            SetName(UnfoldedPt[i][cent], Form("%s_%iGeV", UnfoldedPt[i][cent]->GetName(), i+1));
            SetName(UnfoldedZ[i][cent], Form("%s_%iGeV", UnfoldedZ[i][cent]->GetName(), i+1));


            MeasuredPt[i][cent]->SetDirectory(0);
            MeasuredZ[i][cent]->SetDirectory(0);
            Measured[i][cent]->SetDirectory(0);

            UnfoldedPt[i][cent]->SetDirectory(0);
            UnfoldedZ[i][cent]->SetDirectory(0);
            Unfolded[i][cent]->SetDirectory(0);
        }
    }

    for (int i = 0; i < 5; i++){
        for (int cent = 0; cent < 3; cent++){
            // MeasuredPt[i][cent]->Divide(MeasuredPt[4][cent]);
            // MeasuredZ[i][cent]->Divide(MeasuredZ[4][cent]);
            SetColor(MeasuredPt[i][cent], colors[i]);
            SetColor(MeasuredZ[i][cent], colors[i]);
        }
    }

    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    c->Divide(3);
    for (int cent = 0; cent < 3; cent++){
        c->cd(cent+1);
        // MeasuredPt[0][cent]->GetYaxis()->SetRangeUser(-0.5, 100.5);
        MeasuredPt[0][cent]->Draw();
        MeasuredPt[1][cent]->Draw("same");
        MeasuredPt[2][cent]->Draw("same");
        MeasuredPt[3][cent]->Draw("same");
        MeasuredPt[4][cent]->Draw("same");
        gPad->SetLogy();
        gPad->BuildLegend();
    }

    // cout << CompareHistos(MeasuredPt[0][0], MeasuredPt[1][0]) << endl;
    // cout << CompareHistos(MeasuredPt[1][0], MeasuredPt[2][0]) << endl;
    // cout << CompareHistos(MeasuredPt[2][0], MeasuredPt[3][0]) << endl;
    // cout << CompareHistos(MeasuredPt[3][0], MeasuredPt[4][0]) << endl;

    TCanvas *d = new TCanvas("d", "d", 1200, 600);
    d->Divide(3);
    for (int cent = 0; cent < 3; cent++){
        d->cd(cent+1);
        // MeasuredZ[0][cent]->GetYaxis()->SetRangeUser(-0.5, 100.5);
        MeasuredZ[0][cent]->Draw();
        MeasuredZ[1][cent]->Draw("same");
        MeasuredZ[2][cent]->Draw("same");
        MeasuredZ[3][cent]->Draw("same");
        MeasuredZ[4][cent]->Draw("same");
        gPad->SetLogy();
        gPad->BuildLegend();
    }

    // cout << CompareHistos(MeasuredZ[0][0], MeasuredZ[1][0]) << endl;
    // cout << CompareHistos(MeasuredZ[1][0], MeasuredZ[2][0]) << endl;
    // cout << CompareHistos(MeasuredZ[2][0], MeasuredZ[3][0]) << endl;
    // cout << CompareHistos(MeasuredZ[3][0], MeasuredZ[4][0]) << endl;

    //////////////////////////////////////////////////////////////////////////////////////////

    int lowptbin = UnfoldedPt[0][0]->FindBin(13);
    int highptbin = UnfoldedPt[0][0]->FindBin(13);

    cout << UnfoldedPt[0][0]->Integral(lowptbin, highptbin) << "\t" << UnfoldedPt[0][2]->Integral(lowptbin, highptbin) << endl;
    cout << UnfoldedPt[1][0]->Integral(lowptbin, highptbin) << "\t" << UnfoldedPt[1][2]->Integral(lowptbin, highptbin) << endl;
    cout << UnfoldedPt[2][0]->Integral(lowptbin, highptbin) << "\t" << UnfoldedPt[2][2]->Integral(lowptbin, highptbin) << endl;
    cout << UnfoldedPt[3][0]->Integral(lowptbin, highptbin) << "\t" << UnfoldedPt[3][2]->Integral(lowptbin, highptbin) << endl;
    cout << UnfoldedPt[4][0]->Integral(lowptbin, highptbin) << "\t" << UnfoldedPt[4][2]->Integral(lowptbin, highptbin) << endl;


    cout << CompareHistos(Measured[0][0], Measured[1][0]) << "\t" << CompareHistos(Unfolded[0][0], Unfolded[1][0]) << endl;
    cout << CompareHistos(Measured[1][0], Measured[2][0]) << "\t" << CompareHistos(Unfolded[1][0], Unfolded[2][0]) << endl;
    cout << CompareHistos(Measured[2][0], Measured[3][0]) << "\t" << CompareHistos(Unfolded[2][0], Unfolded[3][0]) << endl;
    cout << CompareHistos(Measured[3][0], Measured[4][0]) << "\t" << CompareHistos(Unfolded[3][0], Unfolded[4][0]) << endl;

    cout << CompareHistos(Measured[0][2], Measured[1][2]) << "\t" << CompareHistos(Unfolded[0][2], Unfolded[1][2]) << endl;
    cout << CompareHistos(Measured[1][2], Measured[2][2]) << "\t" << CompareHistos(Unfolded[1][2], Unfolded[2][2]) << endl;
    cout << CompareHistos(Measured[2][2], Measured[3][2]) << "\t" << CompareHistos(Unfolded[2][2], Unfolded[3][2]) << endl;
    cout << CompareHistos(Measured[3][2], Measured[4][2]) << "\t" << CompareHistos(Unfolded[3][2], Unfolded[4][2]) << endl;
    
}