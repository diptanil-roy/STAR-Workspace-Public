void NewPlot(){
    gROOT->ProcessLine(".x ~/myStyle.C");
    TFile *save1 = new TFile("../Unfold/BKGHistogramsNew.root","READ");
    hBKG_0_101 = (TH1F*)save1->Get("hBKG_0_10");
    hBKG_10_401= (TH1F*)save1->Get("hBKG_10_40");
    hBKG_40_801= (TH1F*)save1->Get("hBKG_40_80");
    hBKG_0_101->SetName("hBKG_0_101");
    hBKG_10_401->SetName("hBKG_10_401");
    hBKG_40_801->SetName("hBKG_40_801");
    TFile *save = new TFile("BKGHistogramsFINAL.root","READ");
    hBKG_0_10 = (TH1F*)save->Get("hBKG_0_10_Weighted");
    hBKG_10_40= (TH1F*)save->Get("hBKG_10_40_Weighted");
    hBKG_40_80= (TH1F*)save->Get("hBKG_40_80_Weighted");

    hBKG_0_10->Rebin();
    hBKG_10_40->Rebin();
    hBKG_40_80->Rebin();
    hBKG_0_101->Rebin();
    hBKG_10_401->Rebin();
    hBKG_40_801->Rebin();

    TCanvas *c1 = new TCanvas("c1","c1");
    hBKG_40_80->DrawNormalized("PE");
    hBKG_40_801->DrawNormalized("hist same");


    TCanvas *c2 = new TCanvas("c12","c12");
    hBKG_10_40->DrawNormalized("PE");
    hBKG_10_401->DrawNormalized("hist same");

    TCanvas *c3 = new TCanvas("c13","c3");
    hBKG_0_10->DrawNormalized("PE");
    hBKG_0_101->DrawNormalized("hist same");

}
