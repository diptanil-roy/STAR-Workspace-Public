R__LOAD_LIBRARY(/Users/diptanilroy/ROOT_INSTALL/RooUnfold/libRooUnfold.so);

using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void QA(){

    TRandom3 *r = new TRandom3(0);

    TH1D *True = new TH1D("True", "True", 10, 0, 1);
    TH1D *Meas = new TH1D("Meas", "Meas", 10, -2, 2);

    TH2D *True2D = new TH2D("True2D", "True2D", 10, 0, 1, 10, 0, 1);
    TH2D *Meas2D = new TH2D("Meas2D", "Meas2D", 10, -2, 2, 10, -5, 5);

    RooUnfoldResponse *resp = new RooUnfoldResponse("Resp", "Resp");
    resp->Setup(Meas, True); //Setup Response Matrix Definition

    RooUnfoldResponse *resp2D = new RooUnfoldResponse("Resp2D", "Resp2D");
    resp2D->Setup(Meas2D, True2D); //Setup Response Matrix Definition

    for (int i = 0; i < 100000; i++){
        double x = r->Rndm();
        double smear = r->Gaus(0, 1);
        double y = x + smear;
        double x2 = 1-x;
        double z = r->Gaus(0, 3);
        double y2 = z/(y*10);
        // double y2 = z;
        True->Fill(x);
        Meas->Fill(y);
        resp->Fill(y,x);

        True2D->Fill(x, x2);
        Meas2D->Fill(y, y2);
        resp2D->Fill(y,y2, x, x2);
    }

    TH1D *TrueTest = new TH1D("TrueTest", "TrueTest", 10, 0, 1);
    TH1D *MeasTest = new TH1D("MeasTest", "MeasTest", 10, -2, 2);

    TH2D *True2DTest = new TH2D("True2DTest", "True2DTest", 10, 0, 1, 10, 0, 1);
    TH2D *Meas2DTest = new TH2D("Meas2DTest", "Meas2DTest", 10, -2, 2, 10, -5, 5);

    for (int i = 0; i < 20000; i++){
        double x = r->Rndm();
        double smear = r->Gaus(0.1, 1.2);
        double y = x + smear;
        double x2 = 1-x;
        double z = r->Gaus(0, 3);
        double y2 = z/(y*10);
        // double y2 = z;
        TrueTest->Fill(x);
        MeasTest->Fill(y);

        True2DTest->Fill(x, x2);
        Meas2DTest->Fill(y, y2);
    }

    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

	RooUnfoldBayes unfold (resp, MeasTest, 4);

    TH1D *Unf = (TH1D *)unfold.Hreco(errorTreatment);

    RooUnfoldBayes unfold2D (resp2D, Meas2DTest, 4);

    TH1D *True2DX = (TH1D *)True2DTest->ProjectionX();
    TH1D *True2DY = (TH1D *)True2DTest->ProjectionY();

    TH2D *Unf2D = (TH2D *)unfold2D.Hreco(errorTreatment);
    TH1D *Unf2DX = (TH1D *)Unf2D->ProjectionX();
    TH1D *Unf2DY = (TH1D *)Unf2D->ProjectionY();

    cout << Unf->Integral() << "\t" << Unf2D->Integral() << endl;

    TCanvas *c = new TCanvas("c", "c", 600, 600);
    c->Divide(2);
    c->cd(1);
    gPad->SetLogy();
    SetColor(Unf, kBlack, 20);
    SetColor(Unf2DX, kBlue, 24);
    SetColor(TrueTest, kGreen-2, 24);
    Unf->Draw();
    Unf2DX->Draw("SAME");
    TrueTest->Draw("SAME");
    True2DX->Draw("SAME");
    c->cd(2);
    gPad->SetLogy();
    SetColor(Unf2DY, kBlue, 24);
    SetColor(True2DY, kGreen-2, 24);
    Unf2DY->Draw();
    True2DY->Draw("SAME");
}