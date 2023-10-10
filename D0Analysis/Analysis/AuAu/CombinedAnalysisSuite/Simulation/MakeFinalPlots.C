using namespace std;

#include "BinDef.h"
#include "NewBinDef.h"

void SetPlotProperties(TH1D *h, int color, int marker, double size = 1.5, int style = kSolid, int width = 1){
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineStyle(style);
    h->SetLineWidth(width);
}

void SetPlotProperties(TGraphAsymmErrors *h, int color, int marker, double size = 1.5, int style = kSolid, int width = 1){
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(size);
    h->SetLineColor(color);
    h->SetLineStyle(style);
    h->SetLineWidth(width);
}

void TrimAndScale(TGraphAsymmErrors *h, double scale, double low, double high){
    for (int i = h->GetN()-1; i >= 0; i--){
        double x, y;
        h->GetPoint(i, x, y);
        if (x < low || x > high) h->RemovePoint(i);
        else {
            h->SetPoint(i, x, y*scale);
            h->SetPointEYhigh(i, h->GetErrorYhigh(i)*scale);
            h->SetPointEYlow(i, h->GetErrorYlow(i)*scale);
        }
    }
}

void TrimAndScale(TH1 *h, TGraphAsymmErrors *g, double scale, double low, double high){
    g = new TGraphAsymmErrors(h);
    TrimAndScale(g, scale, low, high);
}

void CanvasPartition(TCanvas *C,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   int Nx = 1;
 
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
 
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin);
//    Float_t hStep  = 0;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   hposl = 0.0;
   hposr = lMargin + hStep;
   hfactor = hposr-hposl;
   hmarl = lMargin;
   hmarr = rMargin;
 
    for (Int_t j=0;j<Ny;j++) {

        if (j==0) {
        vposd = 0.0;
        vposu = bMargin + vStep;
        vfactor = vposu-vposd;
        vmard = bMargin / vfactor;
        vmaru = 0.0;
        } else if (j == Ny-1) {
        vposd = vposu + vSpacing;
        vposu = vposd + vStep + tMargin;
        vfactor = vposu-vposd;
        vmard = 0.0;
        vmaru = tMargin / (vposu-vposd);
        } else {
        vposd = vposu + vSpacing;
        vposu = vposd + vStep;
        vfactor = vposu-vposd;
        vmard = 0.0;
        vmaru = 0.0;
        }

        C->cd(0);

        cout << "i = " << Nx-1 << ", j = " << j << "\t" << hposl << "\t" << vposd << "\t" << hposr << "\t" << vposu << endl;

        auto name = TString::Format("pad_%d_%d",Nx-1,j);
        auto pad = (TPad*) C->FindObject(name.Data());
        if (pad) delete pad;
        pad = new TPad(name.Data(),"",hposl,vposd,hposr,vposu);
        pad->SetLeftMargin(hmarl);
        pad->SetRightMargin(hmarr);
        pad->SetBottomMargin(vmard);
        pad->SetTopMargin(vmaru);

        pad->SetFrameBorderMode(0);
        pad->SetBorderMode(0);
        pad->SetBorderSize(0);

        pad->Draw();
    }
}
 
double XtoPad(double x)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double pw = xu-xl;
   double lm = gPad->GetLeftMargin();
   double rm = gPad->GetRightMargin();
   double fw = pw-pw*lm-pw*rm;
   return (x*fw+pw*lm)/pw;
}
 
double YtoPad(double y)
{
   double xl,yl,xu,yu;
   gPad->GetPadPar(xl,yl,xu,yu);
   double ph = yu-yl;
   double tm = gPad->GetTopMargin();
   double bm = gPad->GetBottomMargin();
   double fh = ph-ph*bm-ph*tm;
   return (y*fh+bm*ph)/ph;
}


void MakeFinalPlots(){

    float lmargin = 0.22;
    float rmargin = 0.05;
    float bmargin = 0.14;
    float tmargin = 0.07;

    gStyle->SetPaperSize(20,26);
    gStyle->SetPadTopMargin(tmargin);
    gStyle->SetPadRightMargin(rmargin); // increase for colz plots
    gStyle->SetPadBottomMargin(bmargin);
    gStyle->SetPadLeftMargin(lmargin);

    gStyle->SetLabelOffset(0.010,"X");
    gStyle->SetLabelOffset(0.010,"Y");

    gStyle->SetTitleOffset(1.1,"X");
    gStyle->SetLabelOffset(1.1,"Y");

    // put tick marks on top and RHS of plots
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // histogram divisions: only 5 in x to avoid label overlaps
    gStyle->SetNdivisions(505,"x");
    gStyle->SetNdivisions(504,"y");
    gStyle->SetOptStat(0);  
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);

    // gROOT->ForceStyle();

    double taa[3] = {941.23714, 391.35550, 56.62475};

    TFile *f = new TFile("Plots/FinalPlots.root");

    TString CentName[4] = {"0-10%", "10-40%", "40-80%", "pp"};

    TGraphAsymmErrors *FinalJetPt[5][3];
    TGraphAsymmErrors *FinalJetZ[5][3];
    TGraphAsymmErrors *FinalJetdR[5][3];

    TGraphAsymmErrors *FinalRCP_Pt[5][2];
    TGraphAsymmErrors *FinalRCP_Z[5][2];
    TGraphAsymmErrors *FinalRCP_dR[5][2];

    TGraphAsymmErrors *FinalJetPt_Sys[5][3];
    TGraphAsymmErrors *FinalJetZ_Sys[5][3];
    TGraphAsymmErrors *FinalJetdR_Sys[5][3];

    TGraphAsymmErrors *FinalRCP_Pt_Sys[5][2];
    TGraphAsymmErrors *FinalRCP_Z_Sys[5][2];
    TGraphAsymmErrors *FinalRCP_dR_Sys[5][2];

    TGraphAsymmErrors *TheoryD0JetPt[5][4];
    TGraphAsymmErrors *TheoryD0JetZ[5][4];
    TGraphAsymmErrors *TheoryD0JetZPerJet[5][4];
    TGraphAsymmErrors *TheoryD0JetdR[5][4];

    TGraphAsymmErrors *TheoryRCP_JetPt[5][2];
    TGraphAsymmErrors *TheoryRCP_JetZ[5][2];
    TGraphAsymmErrors *TheoryRCP_JetdR[5][2];

    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 3; cent++){
            f->cd(Form("D0pT_%i_%i", pT, 10));
            FinalJetPt[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("Jet pT Spectra %s Stat", CentName[cent].Data()));
            FinalJetPt_Sys[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("Jet pT Spectra %s Sys", CentName[cent].Data()));
            FinalJetZ[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("Jet Z Spectra %s Stat", CentName[cent].Data()));
            FinalJetZ_Sys[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("Jet Z Spectra %s Sys", CentName[cent].Data()));
            FinalJetdR[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("Jet dR %s Stat", CentName[cent].Data()));
            FinalJetdR_Sys[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("Jet dR %s Sys", CentName[cent].Data()));
        }

        for (int cent = 0; cent < 4; cent++){
            if (cent == 0 && (pT != 1 && pT != 4)) continue;
            if (cent == 1) continue;
            if (cent == 2 && (pT != 1 && pT != 4)) continue;
            if (cent == 3 && pT != 1) continue;
            TheoryD0JetPt[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("LIDO Jet pT Spectra %s Stat", CentName[cent].Data()));
            TheoryD0JetZ[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("LIDO Jet Z Spectra %s Stat", CentName[cent].Data()));
            TheoryD0JetdR[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("LIDO Jet dR %s Stat", CentName[cent].Data()));
        }
    }
    cout << "Here" << endl;
    for (int pT = 1; pT <= 5; pT++){
        for (int cent = 0; cent < 2; cent++){
            f->cd(Form("D0pT_%i_%i", pT, 10));
            FinalRCP_Pt[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP Jet pT Spectra %s Stat", CentName[cent].Data()));
            FinalRCP_Pt_Sys[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP Jet pT Spectra %s Sys", CentName[cent].Data()));
            FinalRCP_Z[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP Jet Z Spectra %s Stat", CentName[cent].Data()));
            FinalRCP_Z_Sys[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP Jet Z Spectra %s Sys", CentName[cent].Data()));
            FinalRCP_dR[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP Jet dR %s Stat", CentName[cent].Data()));
            FinalRCP_dR_Sys[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP Jet dR %s Sys", CentName[cent].Data()));

            if (cent == 0 && (pT != 1 && pT != 4)) continue;
            if (cent == 1) continue;
            if (cent == 2 && (pT != 1 && pT != 4)) continue;
            TheoryRCP_JetPt[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP LIDO Jet pT Spectra %s Stat", CentName[cent].Data()));
            TheoryRCP_JetZ[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP LIDO Jet Z Spectra %s Stat", CentName[cent].Data()));
            TheoryRCP_JetdR[pT-1][cent] = (TGraphAsymmErrors*)gDirectory->Get(Form("RCP LIDO Jet dR %s Stat", CentName[cent].Data()));
        }
    }
    cout << "Here" << endl;

    // So there will be 15 Plots - one for each pT bin and observable (i.e jet pT, jet Z, jet dR)
    // Each plot will have 3 centralities in 1 pad. Scale central by 10, midcentral by 1, and peripheral by 0.1

    // Central Marker --> 20, Red For Main Plot, Blue For Systematics
    // Midcentral Marker --> 21, Green-2 For Main Plot, Blue For Systematics
    // Peripheral Marker --> 22, Black For Main Plot, Blue For Systematics

    int color[3] = {kRed, kGreen-2, kBlack};
    int marker[3] = {20, 21, 22};
    double scalefactors[3] = {100, 10, 1};

    TFile *PYTHIA_Theory = new TFile("Theory/D0Jet_Out.root");
    TH1D *PYTHIA_JetPt = (TH1D *)gDirectory->Get("LeadJetPt");
    TH1D *PYTHIA_JetZ = (TH1D *)gDirectory->Get("LeadJetZ");
    TH1D *PYTHIA_dR = (TH1D *)gDirectory->Get("LeadJetdR");

    cout << PYTHIA_JetPt->Integral() << "\t" << PYTHIA_JetZ->Integral() << endl;
    PYTHIA_JetPt->Scale(1./36243862.);
    PYTHIA_JetPt->Scale(1.967/64.6921);
    PYTHIA_JetPt->Scale(43.82/64.6921);
    PYTHIA_JetPt->Scale(0.039);
    PYTHIA_JetPt->Scale(taa[2]*scalefactors[2]);
    ProcessSpectra(PYTHIA_JetPt);
    
    PYTHIA_JetZ->Scale(1./36243862.);
    PYTHIA_JetZ->Scale(1.967/64.6921);
    PYTHIA_JetZ->Scale(43.82/64.6921);
    PYTHIA_JetZ->Scale(0.039);
    PYTHIA_JetZ->Scale(taa[2]*scalefactors[2]);
    ProcessSpectra(PYTHIA_JetZ);
    
    PYTHIA_dR->Scale(1./PYTHIA_dR->Integral());
    DivideByBinWidth(PYTHIA_dR);

    TFile *D0PaperRCP = new TFile("D0Paper_HEPData.root");
    D0PaperRCP->cd("RCP (40-60%)");
    TGraphAsymmErrors *D0RCP010 = (TGraphAsymmErrors *)gDirectory->Get("Graph1D_y1");

    cout << "Imported All Histograms" << endl;

    gSystem->mkdir("Plots/Preliminaries");

    TCanvas *PtSpectra[5];
    TLegend *PtSpectraLegend;
    TPaveText *PtSpectraLabel[5];

    for (int pT = 1; pT <= 5; pT++){
        cout << "pT = " << pT << endl;
        PtSpectra[pT-1] = new TCanvas(Form("PtSpectra_%i", pT), Form("PtSpectra_%i", pT), 800, 800);
        if (pT == 1) {
            PtSpectraLegend = new TLegend(0.5, 0.7, 0.9, 0.9);
            PtSpectraLegend->SetBorderSize(0);
            PtSpectraLegend->SetFillStyle(0);
            PtSpectraLegend->SetTextFont(43);
            PtSpectraLegend->SetTextSize(24);
            PtSpectraLegend->AddEntry(FinalJetPt[pT-1][0], "Central (0-10 %)(#times 100)", "lp");
            PtSpectraLegend->AddEntry(FinalJetPt[pT-1][1], "Midcentral (10-40 %)(#times 10)", "lp");
            PtSpectraLegend->AddEntry(FinalJetPt[pT-1][2], "Peripheral (40-80 %)(#times 1)", "lp");
            PtSpectraLegend->AddEntry(FinalJetPt_Sys[pT-1][0], "Sys. Unc.", "f");
        }

        PtSpectraLabel[pT-1] = new TPaveText(0.32, 0.15, 0.62, 0.35, "NDC");
        PtSpectraLabel[pT-1]->SetBorderSize(0);
        PtSpectraLabel[pT-1]->SetFillStyle(0);
        PtSpectraLabel[pT-1]->SetTextFont(43);
        PtSpectraLabel[pT-1]->SetTextSize(24);
        TText *t1 = PtSpectraLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        PtSpectraLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        PtSpectraLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        PtSpectraLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));

        cout << "Added Legend" << endl;

        gPad->SetLogy();
        TrimAndScale(FinalJetPt[pT-1][0], scalefactors[0], 5, 20);
        TrimAndScale(FinalJetPt_Sys[pT-1][0], scalefactors[0], 5, 20);
        TrimAndScale(FinalJetPt[pT-1][1], scalefactors[1], 5, 20);
        TrimAndScale(FinalJetPt_Sys[pT-1][1], scalefactors[1], 5, 20);
        TrimAndScale(FinalJetPt[pT-1][2], scalefactors[2], 5, 20);
        TrimAndScale(FinalJetPt_Sys[pT-1][2], scalefactors[2], 5, 20);
        SetPlotProperties(FinalJetPt[pT-1][0], color[0], marker[0]);
        SetPlotProperties(FinalJetPt[pT-1][1], color[1], marker[1]);
        SetPlotProperties(FinalJetPt[pT-1][2], color[2], marker[2]);
        cout << "Set Plot Properties" << endl;
        FinalJetPt_Sys[pT-1][0]->Draw("A2");
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetLimits(4.7, 20);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetRangeUser(5*1e-11,3*1e-1);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalJetPt_Sys[pT-1][1]->Draw("2 SAME");
        FinalJetPt_Sys[pT-1][2]->Draw("2 SAME");
        FinalJetPt[pT-1][0]->Draw("P SAME");
        FinalJetPt[pT-1][1]->Draw("P SAME");
        FinalJetPt[pT-1][2]->Draw("P SAME");
        PtSpectraLabel[pT-1]->Draw("SAME");
        cout << "Drew Data Plots" << endl;
        SetColor(PYTHIA_JetPt, kMagenta, 23, 1.5, kMagenta);
        PYTHIA_JetPt->Draw("E3 SAME");
        if (pT == 1){
            cout << "LIDO pp Loop" << endl;
            cout << TheoryD0JetPt[pT-1][3]->GetName() << "\t" << TheoryD0JetPt[pT-1][3]->GetN() << "\t" << TheoryD0JetPt[pT-1][3]->Integral() << endl;
            TrimAndScale(TheoryD0JetPt[pT-1][3], taa[2], 5, 20);
            SetPlotProperties(TheoryD0JetPt[pT-1][3], kCyan, 1, 1, kSolid, 3);
            TheoryD0JetPt[pT-1][3]->Draw("3L SAME");
        }
        PtSpectraLegend->Draw("SAME");
        PtSpectra[pT-1]->SaveAs(Form("Plots/Preliminaries/PtSpectra_%i.pdf", pT));
    }

    TCanvas *PtSpectraTAAScaled[5];
    TGraphAsymmErrors *FinalJetPtTAAScaled[5][3];
    TGraphAsymmErrors *FinalJetPtTAAScaledSys[5][3];
    TString JetPtYaxisNameForTAAScaled = "#frac{1}{T_{AA}} #frac{1}{N_{Evt}} #frac{d^{2}N}{p_{T}dp_{T,Jet}d#eta_{Jet}} [GeV/#it{c}]^{-2}";

    for (int pT = 1; pT <= 5; pT++){
        PtSpectraTAAScaled[pT-1] = new TCanvas(Form("PtSpectraTAAScaled_%i", pT), Form("PtSpectraTAAScaled_%i", pT), 800, 800);
        gPad->SetLogy();
        FinalJetPtTAAScaled[pT-1][0] = (TGraphAsymmErrors*)FinalJetPt[pT-1][0]->Clone(Form("%s TAA Scaled", FinalJetPt[pT-1][0]->GetName()));
        FinalJetPtTAAScaled[pT-1][1] = (TGraphAsymmErrors*)FinalJetPt[pT-1][1]->Clone(Form("%s TAA Scaled", FinalJetPt[pT-1][1]->GetName()));
        FinalJetPtTAAScaled[pT-1][2] = (TGraphAsymmErrors*)FinalJetPt[pT-1][2]->Clone(Form("%s TAA Scaled", FinalJetPt[pT-1][2]->GetName()));
        FinalJetPtTAAScaledSys[pT-1][0] = (TGraphAsymmErrors*)FinalJetPt_Sys[pT-1][0]->Clone(Form("%s TAA Scaled", FinalJetPt_Sys[pT-1][0]->GetName()));
        FinalJetPtTAAScaledSys[pT-1][1] = (TGraphAsymmErrors*)FinalJetPt_Sys[pT-1][1]->Clone(Form("%s TAA Scaled", FinalJetPt_Sys[pT-1][1]->GetName()));
        FinalJetPtTAAScaledSys[pT-1][2] = (TGraphAsymmErrors*)FinalJetPt_Sys[pT-1][2]->Clone(Form("%s TAA Scaled", FinalJetPt_Sys[pT-1][2]->GetName()));
        FinalJetPtTAAScaled[pT-1][0]->Scale(1./taa[0]);
        FinalJetPtTAAScaled[pT-1][0]->Scale(1./scalefactors[0]);
        FinalJetPtTAAScaled[pT-1][1]->Scale(1./taa[1]);
        FinalJetPtTAAScaled[pT-1][1]->Scale(1./scalefactors[1]);
        FinalJetPtTAAScaled[pT-1][2]->Scale(1./taa[2]);
        FinalJetPtTAAScaled[pT-1][2]->Scale(1./scalefactors[2]);
        FinalJetPtTAAScaledSys[pT-1][0]->Scale(1./taa[0]);
        FinalJetPtTAAScaledSys[pT-1][0]->Scale(1./scalefactors[0]);
        FinalJetPtTAAScaledSys[pT-1][1]->Scale(1./taa[1]);
        FinalJetPtTAAScaledSys[pT-1][1]->Scale(1./scalefactors[1]);
        FinalJetPtTAAScaledSys[pT-1][2]->Scale(1./taa[2]);
        // FinalJetPtTAAScaledSys[pT-1][2]->Scale(1./scalefactors[2]);
        SetPlotProperties(FinalJetPtTAAScaled[pT-1][0], color[0], marker[0]);
        SetPlotProperties(FinalJetPtTAAScaled[pT-1][1], color[1], marker[1]);
        SetPlotProperties(FinalJetPtTAAScaled[pT-1][2], color[2], marker[2]);
        FinalJetPtTAAScaledSys[pT-1][0]->GetYaxis()->SetTitle(JetPtYaxisNameForTAAScaled.Data());
        FinalJetPtTAAScaledSys[pT-1][0]->Draw("A2");
        FinalJetPtTAAScaledSys[pT-1][0]->GetXaxis()->SetLimits(4.7, 20);
        FinalJetPtTAAScaledSys[pT-1][0]->GetYaxis()->SetRangeUser(5*1e-15,3*1e-5);
        FinalJetPtTAAScaledSys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalJetPtTAAScaledSys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalJetPtTAAScaledSys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalJetPtTAAScaledSys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalJetPtTAAScaledSys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalJetPtTAAScaledSys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalJetPtTAAScaledSys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalJetPtTAAScaledSys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalJetPtTAAScaledSys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalJetPtTAAScaledSys[pT-1][1]->Draw("2 SAME");
        FinalJetPtTAAScaledSys[pT-1][2]->Draw("2 SAME");
        FinalJetPtTAAScaled[pT-1][0]->Draw("P SAME");
        FinalJetPtTAAScaled[pT-1][1]->Draw("P SAME");
        FinalJetPtTAAScaled[pT-1][2]->Draw("P SAME");
        PtSpectraLabel[pT-1]->Draw("SAME");
        PtSpectraLegend->Draw("SAME");
        PtSpectraTAAScaled[pT-1]->SaveAs(Form("Plots/Preliminaries/PtSpectraTAAScaled_%i.pdf", pT));
    }

    TCanvas *PtSpectraWithTheory[5];
    TLegend *PtSpectraWithTheoryLegend;

    for (int pT = 1; pT <= 5; pT++){
        if (pT != 1 && pT != 4) continue;
        PtSpectra[pT-1] = new TCanvas(Form("PtSpectraWithTheory_%i", pT), Form("PtSpectraWithTheory_%i", pT), 800, 800);
        PtSpectraLabel[pT-1] = new TPaveText(0.32, 0.15, 0.62, 0.35, "NDC");
        PtSpectraLabel[pT-1]->SetBorderSize(0);
        PtSpectraLabel[pT-1]->SetFillStyle(0);
        PtSpectraLabel[pT-1]->SetTextFont(43);
        PtSpectraLabel[pT-1]->SetTextSize(24);
        TText *t1 = PtSpectraLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        PtSpectraLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        PtSpectraLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        PtSpectraLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));

        gPad->SetLogy();
        TrimAndScale(TheoryD0JetPt[pT-1][0], scalefactors[0], 5, 20);
        SetPlotProperties(TheoryD0JetPt[pT-1][0], kMagenta, 1, 1, kSolid, 3);
        TheoryD0JetPt[pT-1][0]->SetFillColorAlpha(kMagenta, 0.5);
        TrimAndScale(TheoryD0JetPt[pT-1][2], scalefactors[2], 5, 20);
        SetPlotProperties(TheoryD0JetPt[pT-1][2], kMagenta, 1, 1, kSolid, 3);
        TheoryD0JetPt[pT-1][2]->SetFillColorAlpha(kMagenta, 0.5);
        if (pT == 1) {
            PtSpectraWithTheoryLegend = new TLegend(0.43, 0.7, 0.9, 0.9);
            PtSpectraWithTheoryLegend->SetBorderSize(0);
            PtSpectraWithTheoryLegend->SetFillStyle(0);
            PtSpectraWithTheoryLegend->SetTextFont(43);
            PtSpectraWithTheoryLegend->SetTextSize(24);
            PtSpectraWithTheoryLegend->AddEntry(FinalJetPt[pT-1][0], "Central (0-10 %)(#times 100)", "lp");
            PtSpectraWithTheoryLegend->AddEntry(FinalJetPt[pT-1][1], "Midcentral (10-40 %)(#times 10)", "lp");
            PtSpectraWithTheoryLegend->AddEntry(FinalJetPt[pT-1][2], "Peripheral (40-80 %)(#times 1)", "lp");
            PtSpectraWithTheoryLegend->AddEntry(FinalJetPt_Sys[pT-1][0], "Sys. Unc.", "f");
            PtSpectraWithTheoryLegend->AddEntry(TheoryD0JetPt[pT-1][0], "LIDO (MPI = off) Stat. Unc. Only", "l");
        }
        FinalJetPt_Sys[pT-1][0]->Draw("A2");
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetLimits(4.7, 20);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetRangeUser(5*1e-11,3*1e-1);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalJetPt_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalJetPt_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalJetPt_Sys[pT-1][1]->Draw("2 SAME");
        TheoryD0JetPt[pT-1][0]->Draw("3L SAME");
        FinalJetPt_Sys[pT-1][2]->Draw("2 SAME");
        TheoryD0JetPt[pT-1][2]->Draw("3L SAME");
        FinalJetPt[pT-1][0]->Draw("P SAME");
        FinalJetPt[pT-1][1]->Draw("P SAME");
        FinalJetPt[pT-1][2]->Draw("P SAME");
        
        // TheoryD0JetPt[pT-1][0]->SetErrorX(0);
        
        PtSpectraLabel[pT-1]->Draw("SAME");
        PtSpectraWithTheoryLegend->Draw("SAME");
        PtSpectra[pT-1]->SaveAs(Form("Plots/Preliminaries/PtSpectraWithTheory_%i.pdf", pT));
    }

    TCanvas *ZSpectra[5];
    TLegend *ZSpectraLegend;
    TPaveText *ZSpectraLabel[5];

    for (int pT = 1; pT <= 5; pT++){
        ZSpectra[pT-1] = new TCanvas(Form("ZSpectra_%i", pT), Form("ZSpectra_%i", pT), 800, 800);
        if (pT == 1) {
            ZSpectraLegend = new TLegend(0.22, 0.72, 0.62, 0.92);
            ZSpectraLegend->SetBorderSize(0);
            ZSpectraLegend->SetFillStyle(0);
            ZSpectraLegend->SetTextFont(43);
            ZSpectraLegend->SetTextSize(24);
            ZSpectraLegend->AddEntry(FinalJetZ[pT-1][0], "Central (0-10 %)(#times 100)", "lp");
            ZSpectraLegend->AddEntry(FinalJetZ[pT-1][1], "Midcentral (10-40 %)(#times 10)", "lp");
            ZSpectraLegend->AddEntry(FinalJetZ[pT-1][2], "Peripheral (40-80 %)(#times 1)", "lp");
            ZSpectraLegend->AddEntry(FinalJetZ_Sys[pT-1][0], "Sys. Unc.", "f");
        }

        ZSpectraLabel[pT-1] = new TPaveText(0.55, 0.15, 0.86, 0.4, "NDC");
        ZSpectraLabel[pT-1]->SetBorderSize(0);
        ZSpectraLabel[pT-1]->SetFillStyle(0);
        ZSpectraLabel[pT-1]->SetTextFont(43);
        ZSpectraLabel[pT-1]->SetTextSize(24);
        TText *t1 = ZSpectraLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        ZSpectraLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        ZSpectraLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4,  |#eta_{Jet}| < 0.6"));
        ZSpectraLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));
        ZSpectraLabel[pT-1]->AddText(Form("%i < p_{T,Jet} [GeV/#it{c}] < %i", 5, 20));

        gPad->SetLogy();
        if (pT >= 3){
            TrimAndScale(FinalJetZ[pT-1][0], scalefactors[0], 0.2, 1);
            TrimAndScale(FinalJetZ_Sys[pT-1][0], scalefactors[0], 0.2, 1);
            TrimAndScale(FinalJetZ[pT-1][1], scalefactors[1], 0.2, 1);
            TrimAndScale(FinalJetZ_Sys[pT-1][1], scalefactors[1], 0.2, 1);
            TrimAndScale(FinalJetZ[pT-1][2], scalefactors[2], 0.2, 1);
            TrimAndScale(FinalJetZ_Sys[pT-1][2], scalefactors[2], 0.2, 1);
        }
        else{
            TrimAndScale(FinalJetZ[pT-1][0], scalefactors[0], 0, 1);
            TrimAndScale(FinalJetZ_Sys[pT-1][0], scalefactors[0], 0, 1);
            TrimAndScale(FinalJetZ[pT-1][1], scalefactors[1], 0, 1);
            TrimAndScale(FinalJetZ_Sys[pT-1][1], scalefactors[1], 0, 1);
            TrimAndScale(FinalJetZ[pT-1][2], scalefactors[2], 0, 1);
            TrimAndScale(FinalJetZ_Sys[pT-1][2], scalefactors[2], 0, 1);
        }
        SetPlotProperties(FinalJetZ[pT-1][0], color[0], marker[0]);
        SetPlotProperties(FinalJetZ[pT-1][1], color[1], marker[1]);
        SetPlotProperties(FinalJetZ[pT-1][2], color[2], marker[2]);
        FinalJetZ_Sys[pT-1][0]->GetXaxis()->SetLimits(0, 1.1);
        FinalJetZ_Sys[pT-1][0]->GetYaxis()->SetRangeUser(2*1e-8,9*1e2);
        FinalJetZ_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalJetZ_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalJetZ_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalJetZ_Sys[pT-1][0]->GetYaxis()->SetTitleOffset(1.8);
        FinalJetZ_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalJetZ_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalJetZ_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalJetZ_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalJetZ_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalJetZ_Sys[pT-1][0]->Draw("A2");
        FinalJetZ_Sys[pT-1][1]->Draw("2 SAME");
        FinalJetZ_Sys[pT-1][2]->Draw("2 SAME");
        FinalJetZ[pT-1][0]->Draw("P SAME");
        FinalJetZ[pT-1][1]->Draw("P SAME");
        FinalJetZ[pT-1][2]->Draw("P SAME");
        SetColor(PYTHIA_JetZ, kMagenta, 23, 1.5, kMagenta);
        PYTHIA_JetZ->Draw("E3 SAME");
        if (pT == 1){
            cout << "LIDO pp Loop" << endl;
            cout << TheoryD0JetZ[pT-1][3]->GetName() << "\t" << TheoryD0JetZ[pT-1][3]->GetN() << "\t" << TheoryD0JetZ[pT-1][3]->Integral() << endl;
            TrimAndScale(TheoryD0JetZ[pT-1][3], taa[2], 0, 1);
            SetPlotProperties(TheoryD0JetZ[pT-1][3], kCyan, 1, 1, kSolid, 3);
            TheoryD0JetZ[pT-1][3]->SetFillColorAlpha(kCyan, 0.5);
            TheoryD0JetZ[pT-1][3]->Draw("3L SAME");
        }
        ZSpectraLegend->Draw("SAME");
        ZSpectraLabel[pT-1]->Draw("SAME");
        ZSpectra[pT-1]->SaveAs(Form("Plots/Preliminaries/ZSpectra_%i.pdf", pT));
    }

    TCanvas *ZSpectraTAAScaled[5];
    TGraphAsymmErrors *FinalJetZTAAScaled[5][3];
    TGraphAsymmErrors *FinalJetZTAAScaledSys[5][3];

    TString JetZYaxisNameForTAAScaled = "#frac{1}{T_{AA}} #frac{1}{N_{Evt}} #frac{d^{2}N}{zdz_{Jet}d#eta_{Jet}}";

    for (int pT = 1; pT <= 5; pT++){
        ZSpectraTAAScaled[pT-1] = new TCanvas(Form("ZSpectraTAAScaled_%i", pT), Form("ZSpectraTAAScaled_%i", pT), 800, 800);
        gPad->SetLogy();
        FinalJetZTAAScaled[pT-1][0] = (TGraphAsymmErrors*)FinalJetZ[pT-1][0]->Clone(Form("%s TAA Scaled", FinalJetZ[pT-1][0]->GetName()));
        FinalJetZTAAScaled[pT-1][1] = (TGraphAsymmErrors*)FinalJetZ[pT-1][1]->Clone(Form("%s TAA Scaled", FinalJetZ[pT-1][1]->GetName()));
        FinalJetZTAAScaled[pT-1][2] = (TGraphAsymmErrors*)FinalJetZ[pT-1][2]->Clone(Form("%s TAA Scaled", FinalJetZ[pT-1][2]->GetName()));
        FinalJetZTAAScaledSys[pT-1][0] = (TGraphAsymmErrors*)FinalJetZ_Sys[pT-1][0]->Clone(Form("%s TAA Scaled", FinalJetZ_Sys[pT-1][0]->GetName()));
        FinalJetZTAAScaledSys[pT-1][1] = (TGraphAsymmErrors*)FinalJetZ_Sys[pT-1][1]->Clone(Form("%s TAA Scaled", FinalJetZ_Sys[pT-1][1]->GetName()));
        FinalJetZTAAScaledSys[pT-1][2] = (TGraphAsymmErrors*)FinalJetZ_Sys[pT-1][2]->Clone(Form("%s TAA Scaled", FinalJetZ_Sys[pT-1][2]->GetName()));
        FinalJetZTAAScaled[pT-1][0]->Scale(1./taa[0]);
        FinalJetZTAAScaled[pT-1][0]->Scale(1./scalefactors[0]);
        FinalJetZTAAScaled[pT-1][1]->Scale(1./taa[1]);
        FinalJetZTAAScaled[pT-1][1]->Scale(1./scalefactors[1]);
        FinalJetZTAAScaled[pT-1][2]->Scale(1./taa[2]);
        FinalJetZTAAScaled[pT-1][2]->Scale(1./scalefactors[2]);
        FinalJetZTAAScaledSys[pT-1][0]->Scale(1./taa[0]);
        FinalJetZTAAScaledSys[pT-1][0]->Scale(1./scalefactors[0]);
        FinalJetZTAAScaledSys[pT-1][1]->Scale(1./taa[1]);
        FinalJetZTAAScaledSys[pT-1][1]->Scale(1./scalefactors[1]);
        FinalJetZTAAScaledSys[pT-1][2]->Scale(1./taa[2]);
        FinalJetZTAAScaledSys[pT-1][2]->Scale(1./scalefactors[2]);
        SetPlotProperties(FinalJetZTAAScaled[pT-1][0], color[0], marker[0]);
        SetPlotProperties(FinalJetZTAAScaled[pT-1][1], color[1], marker[1]);
        SetPlotProperties(FinalJetZTAAScaled[pT-1][2], color[2], marker[2]);
        FinalJetZTAAScaledSys[pT-1][0]->GetYaxis()->SetTitle(JetZYaxisNameForTAAScaled.Data());
        FinalJetZTAAScaledSys[pT-1][0]->Draw("A2");
        FinalJetZTAAScaledSys[pT-1][0]->GetXaxis()->SetLimits(0, 1.1);
        FinalJetZTAAScaledSys[pT-1][0]->GetYaxis()->SetRangeUser(2*1e-12,9*1e-2);
        FinalJetZTAAScaledSys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalJetZTAAScaledSys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalJetZTAAScaledSys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalJetZTAAScaledSys[pT-1][0]->GetXaxis()->SetTitleOffset(2.2);
        FinalJetZTAAScaledSys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalJetZTAAScaledSys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalJetZTAAScaledSys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalJetZTAAScaledSys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalJetZTAAScaledSys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalJetZTAAScaledSys[pT-1][1]->Draw("2 SAME");
        FinalJetZTAAScaledSys[pT-1][2]->Draw("2 SAME");
        FinalJetZTAAScaled[pT-1][0]->Draw("P SAME");
        FinalJetZTAAScaled[pT-1][1]->Draw("P SAME");
        FinalJetZTAAScaled[pT-1][2]->Draw("P SAME");
        ZSpectraLabel[pT-1]->Draw("SAME");
        ZSpectraLegend->Draw("SAME");
        ZSpectraTAAScaled[pT-1]->SaveAs(Form("Plots/Preliminaries/ZSpectraTAAScaled_%i.pdf", pT));
    }

    TCanvas *ZSpectraWithTheory[5];
    TLegend *ZSpectraWithTheoryLegend;

    for (int pT = 1; pT <= 5; pT++){
        if (pT != 1 && pT != 4) continue;
        ZSpectraWithTheory[pT-1] = new TCanvas(Form("ZSpectraWithTheory_%i", pT), Form("ZSpectraWithTheory_%i", pT), 800, 800);
        gPad->SetLogy();
        if (pT >= 3){
            TrimAndScale(TheoryD0JetZ[pT-1][0], scalefactors[0], 0.2, 1);
            TrimAndScale(TheoryD0JetZ[pT-1][2], scalefactors[2], 0.2, 1);
        }
        else{
            TrimAndScale(TheoryD0JetZ[pT-1][0], scalefactors[0], 0, 1);
            TrimAndScale(TheoryD0JetZ[pT-1][2], scalefactors[2], 0, 1);
        }
        SetPlotProperties(TheoryD0JetZ[pT-1][0], kMagenta, 1, 1, kSolid, 3);
        TheoryD0JetZ[pT-1][0]->SetFillColorAlpha(kMagenta, 0.5);
        SetPlotProperties(TheoryD0JetZ[pT-1][2], kMagenta, 1, 1, kSolid, 3);
        TheoryD0JetZ[pT-1][2]->SetFillColorAlpha(kMagenta, 0.5);

        if (pT == 1) {
            ZSpectraWithTheoryLegend = new TLegend(0.22, 0.7, 0.62, 0.92);
            ZSpectraWithTheoryLegend->SetBorderSize(0);
            ZSpectraWithTheoryLegend->SetFillStyle(0);
            ZSpectraWithTheoryLegend->SetTextFont(43);
            ZSpectraWithTheoryLegend->SetTextSize(24);
            ZSpectraWithTheoryLegend->AddEntry(FinalJetZ[pT-1][0], "Central (0-10 %)(#times 100)", "lp");
            ZSpectraWithTheoryLegend->AddEntry(FinalJetZ[pT-1][1], "Midcentral (10-40 %)(#times 10)", "lp");
            ZSpectraWithTheoryLegend->AddEntry(FinalJetZ[pT-1][2], "Peripheral (40-80 %)(#times 1)", "lp");
            ZSpectraWithTheoryLegend->AddEntry(FinalJetZ_Sys[pT-1][0], "Sys. Unc.", "f");
            ZSpectraWithTheoryLegend->AddEntry(TheoryD0JetZ[pT-1][0], "LIDO (MPI = off) Stat. Unc. Only", "l");
        }

        FinalJetZ_Sys[pT-1][0]->Draw("A2");
        TheoryD0JetZ[pT-1][0]->Draw("3L SAME");
        FinalJetZ_Sys[pT-1][1]->Draw("2 SAME");
        TheoryD0JetZ[pT-1][2]->Draw("3L SAME");
        FinalJetZ_Sys[pT-1][2]->Draw("2 SAME");
        FinalJetZ[pT-1][0]->Draw("P SAME");
        FinalJetZ[pT-1][1]->Draw("P SAME");
        FinalJetZ[pT-1][2]->Draw("P SAME");
        
        // TheoryD0JetZ[pT-1][0]->SetErrorX(0);
        ZSpectraWithTheoryLegend->Draw("SAME");
        ZSpectraLabel[pT-1]->Draw("SAME");
        ZSpectraWithTheory[pT-1]->SaveAs(Form("Plots/Preliminaries/ZSpectraWithTheory_%i.pdf", pT));
    }

    TCanvas *dRSpectra[5];
    TLegend *dRSpectraLegend;
    TPaveText *dRSpectraLabel[5];

    for (int pT = 1; pT <= 5; pT++){
        dRSpectra[pT-1] = new TCanvas(Form("dRSpectra_%i", pT), Form("dRSpectra_%i", pT), 800, 800);
        if (pT == 1) {
            dRSpectraLegend = new TLegend(0.5, 0.7, 0.9, 0.9);
            dRSpectraLegend->SetBorderSize(0);
            dRSpectraLegend->SetFillStyle(0);
            dRSpectraLegend->SetTextFont(43);
            dRSpectraLegend->SetTextSize(24);
            dRSpectraLegend->AddEntry(FinalJetdR[pT-1][0], "Central (0-10 %)(#times 100)", "lp");
            dRSpectraLegend->AddEntry(FinalJetdR[pT-1][1], "Midcentral (10-40 %)(#times 10)", "lp");
            dRSpectraLegend->AddEntry(FinalJetdR[pT-1][2], "Peripheral (40-80 %)(#times 1)", "lp");
            dRSpectraLegend->AddEntry(FinalJetdR_Sys[pT-1][0], "Sys. Unc.", "f");
        }

        dRSpectraLabel[pT-1] = new TPaveText(0.34, 0.16, 0.56, 0.4, "NDC");
        dRSpectraLabel[pT-1]->SetBorderSize(0);
        dRSpectraLabel[pT-1]->SetFillStyle(0);
        dRSpectraLabel[pT-1]->SetTextFont(43);
        dRSpectraLabel[pT-1]->SetTextSize(24);
        TText *t1 = dRSpectraLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        dRSpectraLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        dRSpectraLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        dRSpectraLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));
        dRSpectraLabel[pT-1]->AddText(Form("%i < p_{T,Jet} [GeV/#it{c}] < %i", 5, 20));


        gPad->SetLogy();
        TrimAndScale(FinalJetdR[pT-1][0], scalefactors[0], 0, 0.2);
        TrimAndScale(FinalJetdR_Sys[pT-1][0], scalefactors[0], 0, 0.2);
        TrimAndScale(FinalJetdR[pT-1][1], scalefactors[1], 0, 0.2);
        TrimAndScale(FinalJetdR_Sys[pT-1][1], scalefactors[1], 0, 0.2);
        TrimAndScale(FinalJetdR[pT-1][2], scalefactors[2], 0, 0.2);
        TrimAndScale(FinalJetdR_Sys[pT-1][2], scalefactors[2], 0, 0.2);
        SetPlotProperties(FinalJetdR[pT-1][0], color[0], marker[0]);
        SetPlotProperties(FinalJetdR[pT-1][1], color[1], marker[1]);
        SetPlotProperties(FinalJetdR[pT-1][2], color[2], marker[2]);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetRangeUser(0, 0.2);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetRangeUser(2*1e-2,9*1e4);
        FinalJetdR_Sys[pT-1][0]->Draw("A2");
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalJetdR_Sys[pT-1][1]->Draw("2 SAME");
        FinalJetdR_Sys[pT-1][2]->Draw("2 SAME");
        FinalJetdR[pT-1][0]->Draw("P SAME");
        FinalJetdR[pT-1][1]->Draw("P SAME");
        FinalJetdR[pT-1][2]->Draw("P SAME");
        SetColor(PYTHIA_dR, kMagenta, 23, 1.5, kMagenta);
        PYTHIA_dR->Draw("E3 SAME");
        if (pT == 1){
            cout << "LIDO pp Loop" << endl;
            cout << TheoryD0JetdR[pT-1][3]->GetName() << "\t" << TheoryD0JetdR[pT-1][3]->GetN() << "\t" << TheoryD0JetdR[pT-1][3]->Integral() << endl;
            TrimAndScale(TheoryD0JetdR[pT-1][3], 1, 0, 0.2);
            SetPlotProperties(TheoryD0JetdR[pT-1][3], kCyan, 1, 1, kSolid, 3);
            TheoryD0JetdR[pT-1][3]->SetFillColorAlpha(kCyan, 0.5);
            TheoryD0JetdR[pT-1][3]->Draw("3L SAME");
        }
        dRSpectraLegend->Draw("SAME");
        dRSpectraLabel[pT-1]->Draw("SAME");
        dRSpectra[pT-1]->SaveAs(Form("Plots/Preliminaries/dRSpectra_%i.pdf", pT));
    }

    TCanvas *dRSpectraWithTheory[5];
    TLegend *dRSpectraWithTheoryLegend;

    for (int pT = 1; pT <= 5; pT++){
        if (pT != 1 && pT != 4) continue;
        dRSpectraWithTheory[pT-1] = new TCanvas(Form("dRSpectraWithTheory_%i", pT), Form("dRSpectraWithTheory_%i", pT), 800, 800);

        gPad->SetLogy();
        TrimAndScale(TheoryD0JetdR[pT-1][0], scalefactors[0], 0, 0.2);
        SetPlotProperties(TheoryD0JetdR[pT-1][0], kMagenta, 1, 1, kSolid, 3);
        TheoryD0JetdR[pT-1][0]->SetFillColorAlpha(kMagenta, 0.5);
        TrimAndScale(TheoryD0JetdR[pT-1][2], scalefactors[2], 0, 0.2);
        SetPlotProperties(TheoryD0JetdR[pT-1][2], kMagenta, 1, 1, kSolid, 3);
        TheoryD0JetdR[pT-1][2]->SetFillColorAlpha(kMagenta, 0.5);

        if (pT == 1) {
            dRSpectraWithTheoryLegend = new TLegend(0.45, 0.7, 0.8, 0.9);
            dRSpectraWithTheoryLegend->SetBorderSize(0);
            dRSpectraWithTheoryLegend->SetFillStyle(0);
            dRSpectraWithTheoryLegend->SetTextFont(43);
            dRSpectraWithTheoryLegend->SetTextSize(24);
            dRSpectraWithTheoryLegend->AddEntry(FinalJetdR[pT-1][0], "Central (0-10 %)(#times 100)", "lp");
            dRSpectraWithTheoryLegend->AddEntry(FinalJetdR[pT-1][1], "Midcentral (10-40 %)(#times 10)", "lp");
            dRSpectraWithTheoryLegend->AddEntry(FinalJetdR[pT-1][2], "Peripheral (40-80 %)(#times 1)", "lp");
            dRSpectraWithTheoryLegend->AddEntry(FinalJetdR_Sys[pT-1][0], "Sys. Unc.", "f");
            dRSpectraWithTheoryLegend->AddEntry(TheoryD0JetdR[pT-1][0], "LIDO (MPI = off) Stat. Unc. Only", "l");
        }

        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetRangeUser(0, 0.2);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetRangeUser(2*1e-2,9*1e4);
        FinalJetdR_Sys[pT-1][0]->Draw("A2");
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalJetdR_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalJetdR_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalJetdR_Sys[pT-1][1]->Draw("2 SAME");
        FinalJetdR_Sys[pT-1][2]->Draw("2 SAME");
        TheoryD0JetdR[pT-1][0]->Draw("3L SAME");
        TheoryD0JetdR[pT-1][2]->Draw("3L SAME");
        FinalJetdR[pT-1][0]->Draw("P SAME");
        FinalJetdR[pT-1][1]->Draw("P SAME");
        FinalJetdR[pT-1][2]->Draw("P SAME");
        // TheoryD0JetdR[pT-1][0]->SetErrorX(0);
        dRSpectraWithTheoryLegend->Draw("SAME");
        dRSpectraLabel[pT-1]->Draw("SAME");
        dRSpectraWithTheory[pT-1]->SaveAs(Form("Plots/Preliminaries/dRSpectraWithTheory_%i.pdf", pT));
    }

    TH1D *CentralNCollError = new TH1D("CentralNCollError", "CentralNCollError", 1, 4.7, 4.95);
	CentralNCollError->SetBinContent(1, 1.);
	CentralNCollError->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(26.27357/941.23714, 2)));
	CentralNCollError->SetFillColor(kCyan-2);
	CentralNCollError->SetFillStyle(1001);
	CentralNCollError->SetMarkerStyle(1);
	CentralNCollError->SetMarkerColor(kCyan-2);

	TH1D *MidCentralNCollError = new TH1D("MidCentralNCollError", "MidCentralNCollError", 1, 4.7, 4.95);
	MidCentralNCollError->SetBinContent(1, 1.);
	MidCentralNCollError->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(30.21318/391.35550, 2)));
	MidCentralNCollError->SetFillColor(kCyan-2);
	MidCentralNCollError->SetFillStyle(1001);
	MidCentralNCollError->SetMarkerStyle(1);
	MidCentralNCollError->SetMarkerColor(kCyan-2);

    TCanvas *RCP_Pt[5];
    TCanvas *RMP_Pt[5];
    TLegend *RCP_PtLegend;
    TLegend *RMP_PtLegend;
    TPaveText *RCP_PtLabel[5];
    TPaveText *RMP_PtLabel[5];

    TLine *PtLineAtOne = new TLine(4.7, 1, 20, 1);
    PtLineAtOne->SetLineColorAlpha(kBlack, 0.8);
    PtLineAtOne->SetLineStyle(kDashed);

    for (int pT = 1; pT <= 5; pT++){
        RCP_Pt[pT-1] = new TCanvas(Form("RCP_Pt_%i", pT), Form("RCP_Pt_%i", pT), 800, 800);
        if (pT == 1) {
            RCP_PtLegend = new TLegend(0.75, 0.7, 0.88, 0.9);
            RCP_PtLegend->SetBorderSize(0);
            RCP_PtLegend->SetFillStyle(0);
            RCP_PtLegend->SetTextFont(43);
            RCP_PtLegend->SetTextSize(24);
            RCP_PtLegend->AddEntry(FinalRCP_Pt[pT-1][0], "Data", "lp");
            RCP_PtLegend->AddEntry(FinalRCP_Pt_Sys[pT-1][0], "Sys. Unc.", "f");
            RCP_PtLegend->AddEntry(CentralNCollError, "T_{AA} Unc.", "f");
        }

        RCP_PtLabel[pT-1] = new TPaveText(0.32, 0.7, 0.6, 0.9, "NDC");
        RCP_PtLabel[pT-1]->SetBorderSize(0);
        RCP_PtLabel[pT-1]->SetFillStyle(0);
        RCP_PtLabel[pT-1]->SetTextFont(43);
        RCP_PtLabel[pT-1]->SetTextSize(24);

        TText *t1 = RCP_PtLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        RCP_PtLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        RCP_PtLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        RCP_PtLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));
        FinalRCP_Pt[pT-1][0]->SetMarkerSize(2.0);
        FinalRCP_Pt_Sys[pT-1][0]->GetYaxis()->SetTitle("R_{CP} #left(#frac{0-10 %}{40-80 %}#right)");
        FinalRCP_Pt_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalRCP_Pt_Sys[pT-1][0]->GetXaxis()->SetRangeUser(4.7, 20.1);
        FinalRCP_Pt_Sys[pT-1][0]->GetYaxis()->SetRangeUser(0,2.3);
        // FinalRCP_Pt_Sys[pT-1][0]->SetMaximum(2.3);
        FinalRCP_Pt_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalRCP_Pt_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalRCP_Pt_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalRCP_Pt_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalRCP_Pt_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalRCP_Pt_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalRCP_Pt_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalRCP_Pt_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalRCP_Pt_Sys[pT-1][0]->Draw("A2 SAME");
        FinalRCP_Pt[pT-1][0]->Draw("P SAME");
        PtLineAtOne->Draw("SAME");
        CentralNCollError->Draw("SAME E2");
        RCP_PtLegend->Draw("SAME");
        RCP_PtLabel[pT-1]->Draw("SAME");
        RCP_Pt[pT-1]->SaveAs(Form("Plots/Preliminaries/RCP_Pt_%i.pdf", pT));

        RMP_Pt[pT-1] = new TCanvas(Form("RMP_Pt_%i", pT), Form("RMP_Pt_%i", pT), 800, 800);
        if (pT == 1) {
            RMP_PtLegend = new TLegend(0.24, 0.15, 0.4, 0.35);
            RMP_PtLegend->SetBorderSize(0);
            RMP_PtLegend->SetFillStyle(0);
            RMP_PtLegend->SetTextFont(43);
            RMP_PtLegend->SetTextSize(24);
            RMP_PtLegend->AddEntry(FinalRCP_Pt[pT-1][1], "Data", "lp");
            RMP_PtLegend->AddEntry(FinalRCP_Pt_Sys[pT-1][1], "Sys. Unc.", "f");
            RMP_PtLegend->AddEntry(MidCentralNCollError, "T_{AA} Unc.", "f");
        }

        RMP_PtLabel[pT-1] = new TPaveText(0.52, 0.15, 0.88, 0.35, "NDC");
        RMP_PtLabel[pT-1]->SetBorderSize(0);
        RMP_PtLabel[pT-1]->SetFillStyle(0);
        RMP_PtLabel[pT-1]->SetTextFont(43);
        RMP_PtLabel[pT-1]->SetTextSize(24);
        TText *t2 = RMP_PtLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t2->SetTextColor(kRed);
        RMP_PtLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        RMP_PtLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        RMP_PtLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));

        FinalRCP_Pt[pT-1][1]->SetMarkerSize(2.0);
        
        if (pT != 5) {
            FinalRCP_Pt[pT-1][1]->SetMaximum(20);
            TrimAndScale(FinalRCP_Pt[pT-1][1], 1., 5, 20);

            FinalRCP_Pt_Sys[pT-1][1]->SetMaximum(20);
            TrimAndScale(FinalRCP_Pt_Sys[pT-1][1], 1., 5, 20);
            
        }
        FinalRCP_Pt_Sys[pT-1][1]->GetYaxis()->SetTitle("R_{CP} #left(#frac{10-40 %}{40-80 %}#right)");
        FinalRCP_Pt_Sys[pT-1][1]->GetYaxis()->CenterTitle(true);
        FinalRCP_Pt_Sys[pT-1][1]->GetXaxis()->SetRangeUser(4.7, 20);
        FinalRCP_Pt_Sys[pT-1][1]->GetYaxis()->SetRangeUser(0, 2.3);
        if (pT != 5) {
            FinalRCP_Pt[pT-1][1]->RemovePoint(FinalRCP_Pt_Sys[pT-1][1]->GetN()-1);
            FinalRCP_Pt_Sys[pT-1][1]->RemovePoint(FinalRCP_Pt_Sys[pT-1][1]->GetN()-1);
        }
        FinalRCP_Pt_Sys[pT-1][1]->GetXaxis()->SetNdivisions(605);
        FinalRCP_Pt_Sys[pT-1][1]->GetYaxis()->SetNdivisions(505);
        FinalRCP_Pt_Sys[pT-1][1]->GetXaxis()->SetTitleOffset(2.1);
        FinalRCP_Pt_Sys[pT-1][1]->GetXaxis()->SetTitleOffset(1.1);
        FinalRCP_Pt_Sys[pT-1][1]->GetXaxis()->SetTitleSize(0.05);
        FinalRCP_Pt_Sys[pT-1][1]->GetYaxis()->SetTitleSize(0.05);
        FinalRCP_Pt_Sys[pT-1][1]->GetXaxis()->SetLabelSize(0.05);
        FinalRCP_Pt_Sys[pT-1][1]->GetYaxis()->SetLabelSize(0.05);
        FinalRCP_Pt_Sys[pT-1][1]->Draw("A20 SAME");
        FinalRCP_Pt[pT-1][1]->Draw("P SAME");
        PtLineAtOne->Draw("SAME");
        MidCentralNCollError->Draw("SAME E2");
        RMP_PtLegend->Draw("SAME");
        RMP_PtLabel[pT-1]->Draw("SAME");
        RMP_Pt[pT-1]->SaveAs(Form("Plots/Preliminaries/RMP_Pt_%i.pdf", pT));
    }

    TCanvas *RCP_Pt_WithTheory[5];
    TLegend *RCP_Pt_WithTheoryLegend;
    TPaveText *RCP_Pt_WithTheoryLabel[5];

    for (int pT = 1; pT <= 5; pT++){
        if (pT != 1 && pT != 4) continue;
        RCP_Pt_WithTheory[pT-1] = new TCanvas(Form("RCP_Pt_WithTheory_%i", pT), Form("RCP_Pt_WithTheory_%i", pT), 800, 800);
        if (pT == 1) {
            RCP_Pt_WithTheoryLegend = new TLegend(0.24, 0.6, 0.52, 0.8);
            RCP_Pt_WithTheoryLegend->SetBorderSize(0);
            RCP_Pt_WithTheoryLegend->SetFillStyle(0);
            RCP_Pt_WithTheoryLegend->SetTextFont(43);
            RCP_Pt_WithTheoryLegend->SetTextSize(24);
            RCP_Pt_WithTheoryLegend->AddEntry(FinalRCP_Pt[pT-1][0], "Data", "lp");
            RCP_Pt_WithTheoryLegend->AddEntry(FinalRCP_Pt_Sys[pT-1][0], "Sys. Unc.", "f");
            RCP_Pt_WithTheoryLegend->AddEntry(CentralNCollError, "T_{AA} Unc.", "f");
            RCP_Pt_WithTheoryLegend->AddEntry(TheoryRCP_JetPt[pT-1][0], "LIDO (MPI = off) Stat. Unc. Only", "l");
        }

        RCP_Pt_WithTheoryLabel[pT-1] = new TPaveText(0.5, 0.7, 0.88, 0.9, "NDC");
        RCP_Pt_WithTheoryLabel[pT-1]->SetBorderSize(0);
        RCP_Pt_WithTheoryLabel[pT-1]->SetFillStyle(0);
        RCP_Pt_WithTheoryLabel[pT-1]->SetTextFont(43);
        RCP_Pt_WithTheoryLabel[pT-1]->SetTextSize(24);

        TText *t1 = RCP_Pt_WithTheoryLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        RCP_Pt_WithTheoryLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        RCP_Pt_WithTheoryLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        RCP_Pt_WithTheoryLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));

        SetPlotProperties(TheoryRCP_JetPt[pT-1][0], kMagenta, 1, 1, kSolid, 3);
        TheoryRCP_JetPt[pT-1][0]->SetFillColorAlpha(kMagenta, 0.5);
        FinalRCP_Pt_Sys[pT-1][0]->Draw("A2");
        FinalRCP_Pt[pT-1][0]->Draw("P SAME");
        PtLineAtOne->Draw("SAME");
        CentralNCollError->Draw("SAME E2");
        RCP_Pt_WithTheoryLabel[pT-1]->Draw("SAME");
        TheoryRCP_JetPt[pT-1][0]->Draw("3 SAME");
        RCP_Pt_WithTheoryLegend->Draw("SAME");
        RCP_Pt_WithTheory[pT-1]->SaveAs(Form("Plots/Preliminaries/RCP_Pt_WithTheory_%i.pdf", pT));
    }

    TH1D *CentralNCollErrorForZ = new TH1D("CentralNCollErrorForZ", "CentralNCollErrorForZ", 1, 1.07, 1.1);
	CentralNCollErrorForZ->SetBinContent(1, 1.);
	CentralNCollErrorForZ->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(26.27357/941.23714, 2)));
	CentralNCollErrorForZ->SetFillColor(kCyan-2);
	CentralNCollErrorForZ->SetFillStyle(1001);
	CentralNCollErrorForZ->SetMarkerStyle(1);
	CentralNCollErrorForZ->SetMarkerColor(kCyan-2);

	TH1D *MidCentralNCollErrorForZ = new TH1D("MidCentralNCollErrorForZ", "MidCentralNCollErrorForZ", 1, 1.07, 1.1);
	MidCentralNCollErrorForZ->SetBinContent(1, 1.);
	MidCentralNCollErrorForZ->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(30.21318/391.35550, 2)));
	MidCentralNCollErrorForZ->SetFillColor(kCyan-2);
	MidCentralNCollErrorForZ->SetFillStyle(1001);
	MidCentralNCollErrorForZ->SetMarkerStyle(1);
	MidCentralNCollErrorForZ->SetMarkerColor(kCyan-2);

    
    TCanvas *RCP_Z[5];
    TCanvas *RMP_Z[5];
    TLegend *RCP_ZLegend;
    TLegend *RMP_ZLegend;
    TPaveText *RCP_ZLabel[5];

    TLine *ZLineAtOne = new TLine(0, 1, 1.1, 1);
    ZLineAtOne->SetLineColorAlpha(kBlack, 0.8);
    ZLineAtOne->SetLineStyle(kDashed);

    for (int pT = 1; pT <= 5; pT++){
        RCP_Z[pT-1] = new TCanvas(Form("RCP_Z_%i", pT), Form("RCP_Z_%i", pT), 800, 800);
        if (pT == 1) {
            RCP_ZLegend = new TLegend(0.24, 0.7, 0.64, 0.9);
            RCP_ZLegend->SetBorderSize(0);
            RCP_ZLegend->SetFillStyle(0);
            RCP_ZLegend->SetTextFont(43);
            RCP_ZLegend->SetTextSize(24);
            RCP_ZLegend->AddEntry(FinalRCP_Z[pT-1][0], "Data", "lp");
            RCP_ZLegend->AddEntry(FinalRCP_Z_Sys[pT-1][0], "Sys. Unc.", "f");
            RCP_ZLegend->AddEntry(CentralNCollErrorForZ, "T_{AA} Unc.", "f");
        }

        RCP_ZLabel[pT-1] = new TPaveText(0.55, 0.68, 0.88, 0.92, "NDC");
        RCP_ZLabel[pT-1]->SetBorderSize(0);
        RCP_ZLabel[pT-1]->SetFillStyle(0);
        RCP_ZLabel[pT-1]->SetTextFont(43);
        RCP_ZLabel[pT-1]->SetTextSize(24);

        TText *t1 = RCP_ZLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        RCP_ZLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        RCP_ZLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        RCP_ZLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));
        RCP_ZLabel[pT-1]->AddText(Form("%i < p_{T,Jet} [GeV/#it{c}] < %i", 5, 20));

        FinalRCP_Z[pT-1][0]->SetMarkerSize(2.0);
        if (pT >= 3){
            TrimAndScale(FinalRCP_Z[pT-1][0], 1, 0.2, 1);
            TrimAndScale(FinalRCP_Z_Sys[pT-1][0], 1, 0.2, 1);
        }
        else{
            TrimAndScale(FinalRCP_Z[pT-1][0], 1, 0, 1);
            TrimAndScale(FinalRCP_Z_Sys[pT-1][0], 1, 0, 1);
        }
        FinalRCP_Z_Sys[pT-1][0]->Draw("A2");
        FinalRCP_Z_Sys[pT-1][0]->GetYaxis()->SetTitle("R_{CP} #left(#frac{0-10 %}{40-80 %}#right)");
        FinalRCP_Z_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalRCP_Z_Sys[pT-1][0]->GetXaxis()->SetLimits(0, 1.1);
        FinalRCP_Z_Sys[pT-1][0]->GetYaxis()->SetRangeUser(0, 2.3);
        FinalRCP_Z_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalRCP_Z_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalRCP_Z_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalRCP_Z_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalRCP_Z_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalRCP_Z_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalRCP_Z_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalRCP_Z_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalRCP_Z[pT-1][0]->Draw("P SAME");
        ZLineAtOne->Draw("SAME");
        CentralNCollErrorForZ->Draw("SAME E2");
        RCP_ZLegend->Draw("SAME");
        RCP_ZLabel[pT-1]->Draw("SAME");
        RCP_Z[pT-1]->SaveAs(Form("Plots/Preliminaries/RCP_Z_%i.pdf", pT));

        RMP_Z[pT-1] = new TCanvas(Form("RMP_Z_%i", pT), Form("RMP_Z_%i", pT), 800, 800);
        if (pT == 1) {
            RMP_ZLegend = new TLegend(0.24, 0.15, 0.64, 0.35);
            RMP_ZLegend->SetBorderSize(0);
            RMP_ZLegend->SetFillStyle(0);
            RMP_ZLegend->SetTextFont(43);
            RMP_ZLegend->SetTextSize(24);
            RMP_ZLegend->AddEntry(FinalRCP_Z[pT-1][1], "Data", "lp");
            RMP_ZLegend->AddEntry(FinalRCP_Z_Sys[pT-1][1], "Sys. Unc.", "f");
            RMP_ZLegend->AddEntry(MidCentralNCollErrorForZ, "T_{AA} Unc.", "f");
        }

        FinalRCP_Z[pT-1][1]->SetMarkerSize(2.0);
        if (pT >= 3){
            TrimAndScale(FinalRCP_Z[pT-1][1], 1, 0.2, 1);
            TrimAndScale(FinalRCP_Z_Sys[pT-1][1], 1, 0.2, 1);
        }
        else{
            TrimAndScale(FinalRCP_Z[pT-1][1], 1, 0, 1);
            TrimAndScale(FinalRCP_Z_Sys[pT-1][1], 1, 0, 1);
        }
        FinalRCP_Z_Sys[pT-1][1]->Draw("A2");
        FinalRCP_Z_Sys[pT-1][1]->GetYaxis()->SetTitle("R_{CP} #left(#frac{10-40 %}{40-80 %}#right)");
        FinalRCP_Z_Sys[pT-1][1]->GetYaxis()->CenterTitle(true);
        FinalRCP_Z_Sys[pT-1][1]->GetXaxis()->SetLimits(0, 1.1);
        FinalRCP_Z_Sys[pT-1][1]->GetYaxis()->SetRangeUser(0, 2.3);
        FinalRCP_Z_Sys[pT-1][1]->GetXaxis()->SetNdivisions(605);
        FinalRCP_Z_Sys[pT-1][1]->GetYaxis()->SetNdivisions(505);
        FinalRCP_Z_Sys[pT-1][1]->GetXaxis()->SetTitleOffset(2.1);
        FinalRCP_Z_Sys[pT-1][1]->GetXaxis()->SetTitleOffset(1.1);
        FinalRCP_Z_Sys[pT-1][1]->GetXaxis()->SetTitleSize(0.05);
        FinalRCP_Z_Sys[pT-1][1]->GetYaxis()->SetTitleSize(0.05);
        FinalRCP_Z_Sys[pT-1][1]->GetXaxis()->SetLabelSize(0.05);
        FinalRCP_Z_Sys[pT-1][1]->GetYaxis()->SetLabelSize(0.05);
        FinalRCP_Z[pT-1][1]->Draw("P SAME");
        ZLineAtOne->Draw("SAME");
        MidCentralNCollErrorForZ->Draw("SAME E2");
        RMP_ZLegend->Draw("SAME");
        RCP_ZLabel[pT-1]->Draw("SAME");
        RMP_Z[pT-1]->SaveAs(Form("Plots/Preliminaries/RMP_Z_%i.pdf", pT));
    }

    TCanvas *RCP_Z_WithTheory[5];
    TLegend *RCP_Z_WithTheoryLegend;

    for (int pT = 1; pT <= 5; pT++){
        if (pT != 1 && pT != 4) continue;
        RCP_Z_WithTheory[pT-1] = new TCanvas(Form("RCP_Z_WithTheory_%i", pT), Form("RCP_Z_WithTheory_%i", pT), 800, 800);
        if (pT == 1) {
            RCP_Z_WithTheoryLegend = new TLegend(0.24, 0.69, 0.64, 0.94);
            RCP_Z_WithTheoryLegend->SetBorderSize(0);
            RCP_Z_WithTheoryLegend->SetFillStyle(0);
            RCP_Z_WithTheoryLegend->SetTextFont(43);
            RCP_Z_WithTheoryLegend->SetTextSize(24);
            RCP_Z_WithTheoryLegend->AddEntry(FinalRCP_Z[pT-1][0], "Data", "lp");
            RCP_Z_WithTheoryLegend->AddEntry(FinalRCP_Z_Sys[pT-1][0], "Sys. Unc.", "f");
            RCP_Z_WithTheoryLegend->AddEntry(CentralNCollErrorForZ, "T_{AA} Unc.", "f");
            RCP_Z_WithTheoryLegend->AddEntry(TheoryRCP_JetZ[pT-1][0], "#splitline{LIDO (MPI = off)}{Stat. Unc. Only}", "l");
        }

        SetPlotProperties(TheoryRCP_JetZ[pT-1][0], kMagenta, 1, 1, kSolid, 3);
        TheoryRCP_JetZ[pT-1][0]->SetFillColorAlpha(kMagenta, 0.5);
        if (pT >= 3){
            TrimAndScale(FinalRCP_Z[pT-1][0], 1, 0.2, 1);
            TrimAndScale(FinalRCP_Z_Sys[pT-1][0], 1, 0.2, 1);
        }
        else{
            TrimAndScale(FinalRCP_Z[pT-1][0], 1, 0, 1);
            TrimAndScale(FinalRCP_Z_Sys[pT-1][0], 1, 0, 1);
        }
        FinalRCP_Z_Sys[pT-1][0]->GetYaxis()->SetRangeUser(0, 2.3);
        FinalRCP_Z_Sys[pT-1][0]->Draw("A2");
        FinalRCP_Z[pT-1][0]->Draw("P SAME");
        ZLineAtOne->Draw("SAME");
        CentralNCollErrorForZ->Draw("SAME E2");
        RCP_Z_WithTheoryLegend->Draw("SAME");
        RCP_ZLabel[pT-1]->Draw("SAME");
        TheoryRCP_JetZ[pT-1][0]->Draw("3 SAME");
        RCP_Z_WithTheory[pT-1]->SaveAs(Form("Plots/Preliminaries/RCP_Z_WithTheory_%i.pdf", pT));
    }

    TCanvas *RCP_dR[5];
    TCanvas *RMP_dR[5];
    TLegend *RCP_dRLegend;
    TLegend *RMP_dRLegend;
    TPaveText *RCP_dRLabel[5];

    TLine *dRLineAtOne = new TLine(0, 1, 0.2, 1);
    dRLineAtOne->SetLineColorAlpha(kBlack, 0.8);
    dRLineAtOne->SetLineStyle(kDashed);

    for (int pT = 1; pT <= 5; pT++){
        RCP_dR[pT-1] = new TCanvas(Form("RCP_dR_%i", pT), Form("RCP_dR_%i", pT), 800, 800);
        if (pT == 1) {
            RCP_dRLegend = new TLegend(0.25, 0.15, 0.55, 0.35);
            RCP_dRLegend->SetBorderSize(0);
            RCP_dRLegend->SetFillStyle(0);
            RCP_dRLegend->SetTextFont(43);
            RCP_dRLegend->SetTextSize(24);
            RCP_dRLegend->AddEntry(FinalRCP_dR[pT-1][0], "Data", "lp");
            RCP_dRLegend->AddEntry(FinalRCP_dR_Sys[pT-1][0], "Sys. Unc.", "f");
        }

        RCP_dRLabel[pT-1] = new TPaveText(0.35, 0.65, 0.55, 0.9, "NDC");
        RCP_dRLabel[pT-1]->SetBorderSize(0);
        RCP_dRLabel[pT-1]->SetFillStyle(0);
        RCP_dRLabel[pT-1]->SetTextFont(43);
        RCP_dRLabel[pT-1]->SetTextSize(24);

        TText *t1 = RCP_dRLabel[pT-1]->AddText(Form("#it{STAR Preliminary}"));
        t1->SetTextColor(kRed);
        RCP_dRLabel[pT-1]->AddText(Form("Au+Au #sqrt{s_{NN}} = 200 GeV"));
        RCP_dRLabel[pT-1]->AddText(Form("Full Jets, anti-k_{T}, R = 0.4, |#eta_{Jet}| < 0.6"));
        RCP_dRLabel[pT-1]->AddText(Form("%i < p_{T,D^{0}} [GeV/#it{c}] < %i", pT, 10));
        RCP_dRLabel[pT-1]->AddText(Form("%i < p_{T,Jet} [GeV/#it{c}] < %i", 5, 20));
    
        FinalRCP_dR[pT-1][0]->SetMarkerSize(2.0);
        FinalRCP_dR_Sys[pT-1][0]->Draw("A2");
        FinalRCP_dR_Sys[pT-1][0]->GetYaxis()->SetTitle("#frac{0-10 %}{40-80 %}");
        FinalRCP_dR_Sys[pT-1][0]->GetYaxis()->CenterTitle(true);
        FinalRCP_dR_Sys[pT-1][0]->SetMaximum(0.2);
        TrimAndScale(FinalRCP_dR[pT-1][0], 1, 0, 0.2);
        TrimAndScale(FinalRCP_dR_Sys[pT-1][0], 1, 0, 0.2);
        FinalRCP_dR_Sys[pT-1][0]->GetYaxis()->SetRangeUser(0.7, 1.5);
        FinalRCP_dR_Sys[pT-1][0]->GetXaxis()->SetRangeUser(0, 0.2);
        FinalRCP_dR_Sys[pT-1][0]->GetXaxis()->SetNdivisions(605);
        FinalRCP_dR_Sys[pT-1][0]->GetYaxis()->SetNdivisions(505);
        FinalRCP_dR_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(2.1);
        FinalRCP_dR_Sys[pT-1][0]->GetXaxis()->SetTitleOffset(1.1);
        FinalRCP_dR_Sys[pT-1][0]->GetXaxis()->SetTitleSize(0.05);
        FinalRCP_dR_Sys[pT-1][0]->GetYaxis()->SetTitleSize(0.05);
        FinalRCP_dR_Sys[pT-1][0]->GetXaxis()->SetLabelSize(0.05);
        FinalRCP_dR_Sys[pT-1][0]->GetYaxis()->SetLabelSize(0.05);
        FinalRCP_dR[pT-1][0]->Draw("P SAME");
        dRLineAtOne->Draw("SAME");
        RCP_dRLegend->Draw("SAME");
        RCP_dRLabel[pT-1]->Draw("SAME");
        RCP_dR[pT-1]->SaveAs(Form("Plots/Preliminaries/RCP_dR_%i.pdf", pT));

        RMP_dR[pT-1] = new TCanvas(Form("RMP_dR_%i", pT), Form("RMP_dR_%i", pT), 800, 800);
        if (pT == 1) {
            RMP_dRLegend = new TLegend(0.25, 0.15, 0.45, 0.35);
            RMP_dRLegend->SetBorderSize(0);
            RMP_dRLegend->SetFillStyle(0);
            RMP_dRLegend->SetTextFont(43);
            RMP_dRLegend->SetTextSize(24);
            RMP_dRLegend->AddEntry(FinalRCP_dR[pT-1][1], "Data", "lp");
            RMP_dRLegend->AddEntry(FinalRCP_dR_Sys[pT-1][1], "Sys. Unc.", "f");
        }

        FinalRCP_dR[pT-1][1]->SetMarkerSize(2.0);
        FinalRCP_dR_Sys[pT-1][1]->Draw("A2");
        FinalRCP_dR_Sys[pT-1][1]->GetYaxis()->SetTitle("#frac{10-40 %}{40-80 %}");
        FinalRCP_dR_Sys[pT-1][1]->GetYaxis()->CenterTitle(true);
        FinalRCP_dR_Sys[pT-1][1]->SetMaximum(0.2);
        TrimAndScale(FinalRCP_dR[pT-1][1], 1, 0, 0.2);
        TrimAndScale(FinalRCP_dR_Sys[pT-1][1], 1, 0, 0.2);
        FinalRCP_dR_Sys[pT-1][1]->GetYaxis()->SetRangeUser(0.7, 1.5);
        FinalRCP_dR_Sys[pT-1][1]->GetXaxis()->SetRangeUser(0, 0.2);
        FinalRCP_dR_Sys[pT-1][1]->GetXaxis()->SetNdivisions(605);
        FinalRCP_dR_Sys[pT-1][1]->GetYaxis()->SetNdivisions(505);
        FinalRCP_dR_Sys[pT-1][1]->GetXaxis()->SetTitleOffset(2.1);
        FinalRCP_dR_Sys[pT-1][1]->GetXaxis()->SetTitleOffset(1.1);
        FinalRCP_dR_Sys[pT-1][1]->GetXaxis()->SetTitleSize(0.05);
        FinalRCP_dR_Sys[pT-1][1]->GetYaxis()->SetTitleSize(0.05);
        FinalRCP_dR_Sys[pT-1][1]->GetXaxis()->SetLabelSize(0.05);
        FinalRCP_dR_Sys[pT-1][1]->GetYaxis()->SetLabelSize(0.05);
        FinalRCP_dR[pT-1][1]->Draw("P SAME");
        dRLineAtOne->Draw("SAME");
        RMP_dRLegend->Draw("SAME");
        RCP_dRLabel[pT-1]->Draw("SAME");
        RMP_dR[pT-1]->SaveAs(Form("Plots/Preliminaries/RMP_dR_%i.pdf", pT));
    }

    TCanvas *RCP_dR_WithTheory[5];
    TLegend *RCP_dR_WithTheoryLegend;

    for (int pT = 1; pT <= 5; pT++){
        if (pT != 1 && pT != 4) continue;
        RCP_dR_WithTheory[pT-1] = new TCanvas(Form("RCP_dR_WithTheory_%i", pT), Form("RCP_dR_WithTheory_%i", pT), 800, 800);
        if (pT == 1) {
            RCP_dR_WithTheoryLegend = new TLegend(0.25, 0.15, 0.45, 0.35);
            RCP_dR_WithTheoryLegend->SetBorderSize(0);
            RCP_dR_WithTheoryLegend->SetFillStyle(0);
            RCP_dR_WithTheoryLegend->SetTextFont(43);
            RCP_dR_WithTheoryLegend->SetTextSize(24);
            RCP_dR_WithTheoryLegend->AddEntry(FinalRCP_dR[pT-1][0], "Data", "lp");
            RCP_dR_WithTheoryLegend->AddEntry(FinalRCP_dR_Sys[pT-1][0], "Sys. Unc.", "f");
            RCP_dR_WithTheoryLegend->AddEntry(TheoryRCP_JetdR[pT-1][0], "LIDO (MPI = off) Stat. Unc. Only", "l");
        }

        SetPlotProperties(TheoryRCP_JetdR[pT-1][0], kMagenta, 1, 1, kSolid, 3);
        TheoryRCP_JetdR[pT-1][0]->SetFillColorAlpha(kMagenta, 0.5);
        TrimAndScale(TheoryRCP_JetdR[pT-1][0], 1, 0, 0.2);
        FinalRCP_dR_Sys[pT-1][0]->GetYaxis()->SetRangeUser(0.7, 1.5);
        FinalRCP_dR_Sys[pT-1][0]->Draw("A2");
        FinalRCP_dR[pT-1][0]->Draw("P SAME");
        dRLineAtOne->Draw("SAME");
        RCP_dR_WithTheoryLegend->Draw("SAME");
        RCP_dRLabel[pT-1]->Draw("SAME");
        TheoryRCP_JetdR[pT-1][0]->Draw("3 SAME");
        RCP_dR_WithTheory[pT-1]->SaveAs(Form("Plots/Preliminaries/RCP_dR_WithTheory_%i.pdf", pT));
    }

    // This part is only for 1 < D0 pT < 10 GeV/c
    TCanvas *RCP_Pt_MultiPlot = new TCanvas("RCP_Pt_MultiPlot", "RCP_Pt_MultiPlot", 800, 800);
    RCP_Pt_MultiPlot->SetFillStyle(0);
    RCP_Pt_MultiPlot->cd(0);

    TPad *pad[2];

    cout << "Starting to draw" << endl;

    TH2D *RCP_Pt_Hist = new TH2D("RCP_Pt_Hist", "RCP_Pt_Hist", 10, 4.7, 20.1, 10, 0.0001, 1.9999);
    RCP_Pt_Hist->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
    RCP_Pt_Hist->GetYaxis()->SetLabelOffset(0.015);
    RCP_Pt_Hist->GetXaxis()->SetNdivisions(605);
    RCP_Pt_Hist->GetYaxis()->SetNdivisions(505);
    RCP_Pt_Hist->GetXaxis()->SetTitleOffset(1.1);
    RCP_Pt_Hist->GetXaxis()->SetTitleSize(0.08);
    RCP_Pt_Hist->GetXaxis()->SetLabelSize(0.08);
    RCP_Pt_Hist->GetYaxis()->SetLabelSize(0.08);

    // double ymidpoint = (gPad->GetUymax() - tMargin - gPad->GetUymin() - bMargin)/2.;
    // double ymidpoint = 0.5;

    pad[0] = new TPad(Form("%s_Pad_%i", RCP_Pt_MultiPlot->GetName(), 1),"", 0., 0., 1., 0.5);
    pad[0]->SetLeftMargin(lmargin);
    pad[0]->SetRightMargin(rmargin);
    pad[0]->SetBottomMargin(0.2);
    pad[0]->SetTopMargin(0);
    pad[0]->SetFrameBorderMode(0);
    pad[0]->SetBorderMode(0);
    pad[0]->SetBorderSize(0);
    pad[0]->Draw();
    pad[0]->SetFillStyle(4000);
    pad[0]->SetFrameFillStyle(4000);
    pad[0]->cd();

    RCP_Pt_Hist->Draw("AXIS");

    FinalRCP_Pt_Sys[0][1]->Draw("2 SAME");
    FinalRCP_Pt[0][1]->Draw("P SAME");
    MidCentralNCollError->Draw("SAME E2");
    PtLineAtOne->Draw("SAME");
    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_Pt_MultiPlot->cd(0);

    pad[1] = new TPad(Form("%s_Pad_%i", RCP_Pt_MultiPlot->GetName(), 2),"",0., 0.5, 1., 1.);
    pad[1]->SetLeftMargin(lmargin);
    pad[1]->SetRightMargin(rmargin);
    pad[1]->SetBottomMargin(0);
    pad[1]->SetTopMargin(0.2);
    pad[1]->SetFrameBorderMode(0);
    pad[1]->SetBorderMode(0);
    pad[1]->SetBorderSize(0);
    pad[1]->Draw();
    pad[1]->SetFillStyle(4000);
    pad[1]->SetFrameFillStyle(4000);
    pad[1]->cd();

    RCP_Pt_Hist->Draw("AXIS");

    FinalRCP_Pt_Sys[0][0]->Draw("2 SAME");
    TheoryRCP_JetPt[0][0]->Draw("3 SAME");
    FinalRCP_Pt[0][0]->Draw("P SAME");
    CentralNCollError->Draw("SAME E2");
    PtLineAtOne->Draw("SAME");

    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_Pt_MultiPlot->cd(0);
    TPad *yaxistitlepad = new TPad("yaxistitlepad","yaxistitlepad",0,0,1,1);
    yaxistitlepad->SetFillStyle(4000);  // transparent
    yaxistitlepad->Draw();
    yaxistitlepad->cd();
    TLatex *lat = new TLatex();
    // lat->SetTextFont(43);
    lat->SetTextAngle(90);
    lat->DrawLatexNDC(.1,.4,"R_{CP}(/40-80%)");
    auto RCP_PtLegend2 = (TLegend*)RCP_Pt_WithTheoryLegend->Clone();
    RCP_PtLegend2->SetTextSize(21);
    RCP_PtLegend2->SetY1(0.72);
    RCP_PtLegend2->SetY2(0.88);
    RCP_PtLegend2->Draw("SAME");
    auto RCP_PtLabel2 = (TPaveText*)RCP_PtLabel[0]->Clone();
    RCP_PtLabel2->SetTextSize(20);
    RCP_PtLabel2->SetX1(0.62);
    RCP_PtLabel2->SetX2(0.9);
    RCP_PtLabel2->SetY1(0.11);
    RCP_PtLabel2->SetY2(0.26);
    RCP_PtLabel2->Draw("SAME");

    TLatex *CentralTag = new TLatex();
    CentralTag->SetTextFont(43);
    CentralTag->SetTextSize(24);
    CentralTag->DrawLatexNDC(0.3, 0.52, "Central (0-10%)");
    TLatex *MidCentralTag = new TLatex();
    MidCentralTag->SetTextFont(43);
    MidCentralTag->SetTextSize(24);
    MidCentralTag->DrawLatexNDC(0.3, 0.12, "MidCentral (10-40%)");

    RCP_Pt_MultiPlot->SaveAs("Plots/Preliminaries/RCP_Pt_MultiPlot_1.pdf");

    TLine *PtLineAtOneExtended = new TLine(0.8, 1, 20, 1);
    PtLineAtOneExtended->SetLineColorAlpha(kBlack, 0.8);
    PtLineAtOneExtended->SetLineStyle(kDashed);

    TH1D *CentralNCollErrorExtended = new TH1D("CentralNCollErrorExtended", "CentralNCollErrorExtended", 1, 0.7, 0.9);
	CentralNCollErrorExtended->SetBinContent(1, 1.);
	CentralNCollErrorExtended->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(26.27357/941.23714, 2)));
	CentralNCollErrorExtended->SetFillColor(kCyan-2);
	CentralNCollErrorExtended->SetFillStyle(1001);
	CentralNCollErrorExtended->SetMarkerStyle(1);
	CentralNCollErrorExtended->SetMarkerColor(kCyan-2);

	TH1D *MidCentralNCollErrorExtended = new TH1D("MidCentralNCollErrorExtended", "MidCentralNCollErrorExtended", 1, 0.7, 0.9);
	MidCentralNCollErrorExtended->SetBinContent(1, 1.);
	MidCentralNCollErrorExtended->SetBinError(1, TMath::Sqrt(pow(13.62653/56.62475, 2) + pow(30.21318/391.35550, 2)));
	MidCentralNCollErrorExtended->SetFillColor(kCyan-2);
	MidCentralNCollErrorExtended->SetFillStyle(1001);
	MidCentralNCollErrorExtended->SetMarkerStyle(1);
	MidCentralNCollErrorExtended->SetMarkerColor(kCyan-2);


    // This part is only for 1 < D0 pT < 10 GeV/c
    TCanvas *RCP_Pt_MultiPlot_WithD0RCP = new TCanvas("RCP_Pt_MultiPlot_WithD0RCP", "RCP_Pt_MultiPlot_WithD0RCP", 800, 800);
    RCP_Pt_MultiPlot_WithD0RCP->SetFillStyle(0);
    RCP_Pt_MultiPlot_WithD0RCP->cd(0);

    TPad *pad4[2];

    cout << "Starting to draw" << endl;

    TH2D *RCP_Pt_Hist_ForD0RCP = new TH2D("RCP_Pt_Hist_ForD0RCP", "RCP_Pt_Hist_ForD0RCP", 10, 0.7, 20.1, 10, 0.0001, 2.1999);
    RCP_Pt_Hist_ForD0RCP->GetXaxis()->SetTitle("p_{T,Jet} [GeV/#it{c}]");
    RCP_Pt_Hist_ForD0RCP->GetYaxis()->SetLabelOffset(0.015);
    RCP_Pt_Hist_ForD0RCP->GetXaxis()->SetNdivisions(605);
    RCP_Pt_Hist_ForD0RCP->GetYaxis()->SetNdivisions(505);
    RCP_Pt_Hist_ForD0RCP->GetXaxis()->SetTitleOffset(1.1);
    RCP_Pt_Hist_ForD0RCP->GetXaxis()->SetTitleSize(0.08);
    RCP_Pt_Hist_ForD0RCP->GetXaxis()->SetLabelSize(0.08);
    RCP_Pt_Hist_ForD0RCP->GetYaxis()->SetLabelSize(0.08);

    // double ymidpoint = (gPad->GetUymax() - tMargin - gPad->GetUymin() - bMargin)/2.;
    double ymidpoint = 0.5;

    pad4[0] = new TPad(Form("%s_Pad_%i", RCP_Pt_MultiPlot_WithD0RCP->GetName(), 1),"", 0., 0., 1., 0.5);
    pad4[0]->SetLeftMargin(lmargin);
    pad4[0]->SetRightMargin(rmargin);
    pad4[0]->SetBottomMargin(0.2);
    pad4[0]->SetTopMargin(0);
    pad4[0]->SetFrameBorderMode(0);
    pad4[0]->SetBorderMode(0);
    pad4[0]->SetBorderSize(0);
    pad4[0]->Draw();
    pad4[0]->SetFillStyle(4000);
    pad4[0]->SetFrameFillStyle(4000);
    pad4[0]->cd();

    RCP_Pt_Hist_ForD0RCP->Draw("AXIS");

    FinalRCP_Pt_Sys[0][1]->Draw("2 SAME");
    FinalRCP_Pt[0][1]->Draw("P SAME");
    MidCentralNCollErrorExtended->Draw("SAME E2");
    PtLineAtOneExtended->Draw("SAME");
    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_Pt_MultiPlot_WithD0RCP->cd(0);

    pad4[1] = new TPad(Form("%s_Pad_%i", RCP_Pt_MultiPlot_WithD0RCP->GetName(), 2),"",0., 0.5, 1., 1.);
    pad4[1]->SetLeftMargin(lmargin);
    pad4[1]->SetRightMargin(rmargin);
    pad4[1]->SetBottomMargin(0);
    pad4[1]->SetTopMargin(0.2);
    pad4[1]->SetFrameBorderMode(0);
    pad4[1]->SetBorderMode(0);
    pad4[1]->SetBorderSize(0);
    pad4[1]->Draw();
    pad4[1]->SetFillStyle(4000);
    pad4[1]->SetFrameFillStyle(4000);
    pad4[1]->cd();

    RCP_Pt_Hist_ForD0RCP->Draw("AXIS");

    TrimAndScale(D0RCP010, 1, 1, 10);
    TrimAndScale(FinalRCP_Pt_Sys[0][0], 1, 5, 20);
    TrimAndScale(FinalRCP_Pt[0][0], 1, 5, 20);
    // SetPlotProperties(D0RCP010, kBlack, 1, 1, kSolid, 3);
    SetPlotProperties(D0RCP010, kBlack, 20);

    FinalRCP_Pt_Sys[0][0]->Draw("2 SAME");
    TheoryRCP_JetPt[0][0]->Draw("3 SAME");
    FinalRCP_Pt[0][0]->Draw("P SAME");
    D0RCP010->Draw("P SAME");
    CentralNCollErrorExtended->Draw("SAME E2");
    PtLineAtOneExtended->Draw("SAME");

    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_Pt_MultiPlot_WithD0RCP->cd(0);
    TPad *yaxistitlepad4 = new TPad("yaxistitlepad4","yaxistitlepad4",0,0,1,1);
    yaxistitlepad4->SetFillStyle(4000);  // transparent
    yaxistitlepad4->Draw();
    yaxistitlepad4->cd();
    TLatex *lat4 = new TLatex();
    // lat->SetTextFont(43);
    lat4->SetTextAngle(90);
    lat4->DrawLatexNDC(.1,.4,"R_{CP}(/40-80%)");
    auto RCP_PtLegend3 = (TLegend*)RCP_Pt_WithTheoryLegend->Clone();
    RCP_PtLegend3->SetTextSize(18);
    RCP_PtLegend3->SetY1(0.72);
    RCP_PtLegend3->SetY2(0.9);
    RCP_PtLegend3->AddEntry(D0RCP010, "D^{0} R_{CP} (0-10%)/(40-60%) [PRC 2019]", "lp");
    RCP_PtLegend3->Draw("SAME");
    auto RCP_PtLabel3 = (TPaveText*)RCP_PtLabel[0]->Clone();
    RCP_PtLabel3->SetTextSize(18);
    RCP_PtLabel3->SetX1(0.3);
    RCP_PtLabel3->SetX2(0.5);
    RCP_PtLabel3->SetY1(0.3);
    RCP_PtLabel3->SetY2(0.48);
    RCP_PtLabel3->Draw("SAME");

    TLatex *CentralTag4 = new TLatex();
    CentralTag4->SetTextFont(43);
    CentralTag4->SetTextSize(24);
    CentralTag4->DrawLatexNDC(0.3, 0.52, "Central (0-10%)");
    TLatex *MidCentralTag4 = new TLatex();
    MidCentralTag4->SetTextFont(43);
    MidCentralTag4->SetTextSize(24);
    MidCentralTag4->DrawLatexNDC(0.3, 0.12, "MidCentral (10-40%)");

    RCP_Pt_MultiPlot_WithD0RCP->SaveAs("Plots/Preliminaries/RCP_Pt_MultiPlot_WithD0RCP_1.pdf");

    // This part is only for 1 < D0 pT < 10 GeV/c
    TCanvas *RCP_Z_MultiPlot = new TCanvas("RCP_Z_MultiPlot", "RCP_Z_MultiPlot", 800, 800);
    RCP_Z_MultiPlot->SetFillStyle(0);
    RCP_Z_MultiPlot->cd(0);

    TPad *pad2[2];

    TH2D *RCP_Z_Hist = new TH2D("RCP_Z_Hist", "RCP_Z_Hist", 10, 0.0, 1.1, 10, 0.0001, 2.49999);
    RCP_Z_Hist->GetXaxis()->SetTitle(FinalRCP_Z_Sys[1][0]->GetXaxis()->GetTitle());
    RCP_Z_Hist->GetYaxis()->SetLabelOffset(0.015);
    RCP_Z_Hist->GetXaxis()->SetNdivisions(605);
    RCP_Z_Hist->GetYaxis()->SetNdivisions(505);
    RCP_Z_Hist->GetXaxis()->SetTitleOffset(1.1);
    RCP_Z_Hist->GetXaxis()->SetTitleSize(0.08);
    RCP_Z_Hist->GetXaxis()->SetLabelSize(0.08);
    RCP_Z_Hist->GetYaxis()->SetLabelSize(0.08);


    pad2[0] = new TPad(Form("%s_Pad_%i", RCP_Z_MultiPlot->GetName(), 1),"", 0., 0., 1., 0.5);
    pad2[0]->SetLeftMargin(lmargin);
    pad2[0]->SetRightMargin(rmargin);
    pad2[0]->SetBottomMargin(0.205);
    pad2[0]->SetTopMargin(0);
    pad2[0]->SetFrameBorderMode(0);
    pad2[0]->SetBorderMode(0);
    pad2[0]->SetBorderSize(0);
    pad2[0]->Draw();
    pad2[0]->SetFillStyle(4000);
    pad2[0]->SetFrameFillStyle(4000);
    pad2[0]->cd();
    
    RCP_Z_Hist->Draw("AXIS");

    FinalRCP_Z_Sys[0][1]->Draw("2 SAME");
    FinalRCP_Z[0][1]->Draw("P SAME");
    MidCentralNCollErrorForZ->Draw("SAME E2");
    ZLineAtOne->Draw("SAME");
    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_Z_MultiPlot->cd(0);

    pad2[1] = new TPad(Form("%s_Pad_%i", RCP_Z_MultiPlot->GetName(), 2),"",0., 0.5, 1., 1.);
    pad2[1]->SetLeftMargin(lmargin);
    pad2[1]->SetRightMargin(rmargin);
    pad2[1]->SetBottomMargin(0);
    pad2[1]->SetTopMargin(0.205);
    pad2[1]->SetFrameBorderMode(0);
    pad2[1]->SetBorderMode(0);
    pad2[1]->SetBorderSize(0);
    pad2[1]->Draw();
    pad2[1]->SetFillStyle(4000);
    pad2[1]->SetFrameFillStyle(4000);
    pad2[1]->cd();

    RCP_Z_Hist->Draw("AXIS");

    FinalRCP_Z_Sys[0][0]->Draw("2 SAME");
    TheoryRCP_JetZ[0][0]->Draw("3 SAME");
    FinalRCP_Z[0][0]->Draw("P SAME");
    CentralNCollErrorForZ->Draw("SAME E2");
    ZLineAtOne->Draw("SAME");
    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_Z_MultiPlot->cd(0);

    TPad *yaxistitlepad2 = new TPad("yaxistitlepad2","yaxistitlepad2",0,0,1,1);
    yaxistitlepad2->SetFillStyle(4000);  // transparent
    yaxistitlepad2->Draw();
    yaxistitlepad2->cd();
    
    TLatex *lat2 = new TLatex();
    // lat2->SetTextFont(43);
    lat2->SetTextAngle(90);
    lat2->DrawLatexNDC(.1,.4,"R_{CP}(/40-80%)");
    auto RCP_ZLegend2 = (TLegend*)RCP_PtLegend2->Clone();
    RCP_ZLegend2->SetTextSize(21);
    RCP_ZLegend2->SetY1(0.75);
    RCP_ZLegend2->SetY2(0.9);
    RCP_ZLegend2->Draw("SAME");

    auto RCP_ZLabel2 = (TPaveText*)RCP_ZLabel[0]->Clone();
    RCP_ZLabel2->SetTextSize(20);
    RCP_ZLabel2->SetX1(0.61);
    RCP_ZLabel2->SetX2(0.9);
    RCP_ZLabel2->SetY1(0.3);
    RCP_ZLabel2->SetY2(0.5);
    RCP_ZLabel2->Draw("SAME");

    TLatex *CentralTag2 = new TLatex();
    CentralTag2->SetTextFont(43);
    CentralTag2->SetTextSize(24);
    CentralTag2->DrawLatexNDC(0.3, 0.52, "Central (0-10%)");
    TLatex *MidCentralTag2 = new TLatex();
    MidCentralTag2->SetTextFont(43);
    MidCentralTag2->SetTextSize(24);
    MidCentralTag2->DrawLatexNDC(0.3, 0.12, "MidCentral (10-40%)");

    RCP_Z_MultiPlot->SaveAs("Plots/Preliminaries/RCP_Z_MultiPlot_1.pdf");

    // This part is only for 1 < D0 pT < 10 GeV/c

    TCanvas *RCP_dR_MultiPlot = new TCanvas("RCP_dR_MultiPlot", "RCP_dR_MultiPlot", 800, 800);
    RCP_dR_MultiPlot->SetFillStyle(0);
    RCP_dR_MultiPlot->cd(0);

    TPad *pad3[2];

    TH2D *RCP_dR_Hist = new TH2D("RCP_dR_Hist", "RCP_dR_Hist", 10, 0.0, 0.2, 10, 0.700001, 1.4999);
    RCP_dR_Hist->GetXaxis()->SetTitle(FinalRCP_dR_Sys[1][0]->GetXaxis()->GetTitle());
    RCP_dR_Hist->GetYaxis()->SetLabelOffset(0.015);
    RCP_dR_Hist->GetXaxis()->SetNdivisions(605);
    RCP_dR_Hist->GetYaxis()->SetNdivisions(505);
    RCP_dR_Hist->GetXaxis()->SetTitleOffset(1.1);
    RCP_dR_Hist->GetXaxis()->SetTitleSize(0.08);
    RCP_dR_Hist->GetXaxis()->SetLabelSize(0.08);
    RCP_dR_Hist->GetYaxis()->SetLabelSize(0.08);

    pad3[0] = new TPad(Form("%s_Pad_%i", RCP_dR_MultiPlot->GetName(), 1),"", 0., 0., 1., 0.5);
    pad3[0]->SetLeftMargin(lmargin);
    pad3[0]->SetRightMargin(rmargin);
    pad3[0]->SetBottomMargin(0.2);
    pad3[0]->SetTopMargin(0);
    pad3[0]->SetFrameBorderMode(0);
    pad3[0]->SetBorderMode(0);
    pad3[0]->SetBorderSize(0);
    pad3[0]->Draw();
    pad3[0]->SetFillStyle(4000);
    pad3[0]->SetFrameFillStyle(4000);
    pad3[0]->cd();

    RCP_dR_Hist->Draw("AXIS");

    FinalRCP_dR_Sys[0][1]->Draw("2 SAME");
    FinalRCP_dR[0][1]->Draw("P SAME");
    MidCentralNCollError->Draw("SAME E2");
    dRLineAtOne->Draw("SAME");
    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_dR_MultiPlot->cd(0);

    pad3[1] = new TPad(Form("%s_Pad_%i", RCP_dR_MultiPlot->GetName(), 2),"",0., 0.5, 1., 1.);
    pad3[1]->SetLeftMargin(lmargin);
    pad3[1]->SetRightMargin(rmargin);
    pad3[1]->SetBottomMargin(0);
    pad3[1]->SetTopMargin(0.2);
    pad3[1]->SetFrameBorderMode(0);
    pad3[1]->SetBorderMode(0);
    pad3[1]->SetBorderSize(0);
    pad3[1]->Draw();
    pad3[1]->SetFillStyle(4000);
    pad3[1]->SetFrameFillStyle(4000);
    pad3[1]->cd();

    RCP_dR_Hist->Draw("AXIS");

    FinalRCP_dR_Sys[0][0]->Draw("2 SAME");
    TheoryRCP_JetdR[0][0]->Draw("3 SAME");
    FinalRCP_dR[0][0]->Draw("P SAME");
    CentralNCollError->Draw("SAME E2");
    dRLineAtOne->Draw("SAME");
    // gPad->SetTickx(0);
    gPad->Modified();
    gPad->Update();

    RCP_dR_MultiPlot->cd(0);

    TPad *yaxistitlepad3 = new TPad("yaxistitlepad3","yaxistitlepad3",0,0,1,1);
    yaxistitlepad3->SetFillStyle(4000);  // transparent
    yaxistitlepad3->Draw();
    yaxistitlepad3->cd();

    TLatex *lat3 = new TLatex();
    // lat3->SetTextFont(43);
    lat3->SetTextAngle(90);
    lat3->DrawLatexNDC(.1,.4,"Ratio(/40-80%)");
    auto RCP_dRLegend2 = (TLegend*)RCP_PtLegend2->Clone();
    RCP_dRLegend2->SetTextSize(21);
    RCP_dRLegend2->SetY1(0.75);
    RCP_dRLegend2->SetY2(0.9);
    RCP_dRLegend2->Draw("SAME");
    auto RCP_dRLabel2 = (TPaveText*)RCP_dRLabel[0]->Clone();
    RCP_dRLabel2->SetTextSize(20);
    RCP_dRLabel2->SetY1(0.28);
    RCP_dRLabel2->SetY2(0.48);
    RCP_dRLabel2->SetX1(0.6);
    RCP_dRLabel2->SetX2(0.9);
    RCP_dRLabel2->Draw("SAME");

    TLatex *CentralTag3 = new TLatex();
    CentralTag3->SetTextFont(43);
    CentralTag3->SetTextSize(24);
    CentralTag3->DrawLatexNDC(0.3, 0.52, "Central (0-10%)");
    TLatex *MidCentralTag3 = new TLatex();
    MidCentralTag3->SetTextFont(43);
    MidCentralTag3->SetTextSize(24);
    MidCentralTag3->DrawLatexNDC(0.3, 0.12, "MidCentral (10-40%)");

    RCP_dR_MultiPlot->SaveAs("Plots/Preliminaries/RCP_dR_MultiPlot_1.pdf");

}